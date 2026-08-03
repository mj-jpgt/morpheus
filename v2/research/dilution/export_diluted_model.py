"""Stage 2: re-export a trained V2 checkpoint with tissue-diluted patch bags.

Stage 1 dilutes the unweighted ``concat(mean, std)`` global skip path, which is
structurally unable to down-weight a non-tumour patch.  This script asks the
complementary question: does the learned attention (``SoftSlotPool`` at the
local/slide/patient levels) recover any of what that path loses?

Everything except the bag contents is held fixed -- same checkpoint, same split,
same targets, same export code path -- so the only difference between levels is
the tissue fed in.

Normal patches are attached to the patient's EXISTING slide identifiers, in
proportion to each slide's tumour patch count, rather than to a synthetic extra
slide.  Whole-slide sampling of an unannotated cohort mixes tumour and normal
inside one slide; giving the normal patches their own slide token would hand the
hierarchy a free segregation cue that a real external cohort would not provide.

One artifact per dilution level, named ``<prefix>_dilution_<level>.npz``, so
``run_calibra`` reports each level as its own ``method`` with the canonical state
names intact.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import torch

from morpheus.src.training.train_bio_query_former import load_bio_query_data
from morpheus.v2.contracts import trainable_state_names
from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.preflight import restrict_cohort_to_split, split_file_digest
from morpheus.v2.runner import (UncappedHoptimusBatches, _standardize_clinical,
                                attach_external_programme_targets, attach_v2_targets)

from build_dilution_artifact import load_normal_pool, normal_draw_order, patient_rng


class DilutedHoptimusBatches(UncappedHoptimusBatches):
    """``UncappedHoptimusBatches`` with a fixed non-tumour fraction spliced in."""

    def __init__(self, *args, pool: np.ndarray, donor_lookup, level: float, seed: int, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.pool, self.donor_lookup, self.level, self.draw_seed = pool, donor_lookup, float(level), int(seed)
        self.drawn: dict[str, int] = {}

    def __iter__(self):
        for batch in super().__iter__():
            if self.level <= 0.0:
                yield batch
                continue
            indices = batch["indices"].numpy()
            patches, mask, slide_ids = (batch["patches"].numpy(), batch["patch_mask"].numpy(),
                                        batch["slide_ids"].numpy())
            bags = []
            for row, index in enumerate(indices):
                patient = str(self.data.patient_ids[index])
                valid = mask[row]
                tumour, codes = patches[row][valid], slide_ids[row][valid]
                needed = int(round(self.level / (1.0 - self.level) * len(tumour)))
                rng = patient_rng(self.draw_seed, patient)
                donors = self.donor_lookup(patient, str(self.data.cancers[index]))
                rows = normal_draw_order(donors, rng, needed)
                self.drawn[patient] = int(len(rows))
                if len(rows) == 0:
                    bags.append((tumour, codes))
                    continue
                unique, counts = np.unique(codes, return_counts=True)
                assigned = rng.choice(unique, size=len(rows), p=counts / counts.sum())
                bags.append((np.concatenate([tumour, self.pool[rows]], axis=0),
                             np.concatenate([codes, assigned])))
            width = max(len(bag) for bag, _ in bags)
            new_patches = np.zeros((len(bags), width, patches.shape[-1]), dtype=np.float32)
            new_mask = np.zeros((len(bags), width), dtype=bool)
            new_slides = np.zeros((len(bags), width), dtype=np.int64)
            for row, (bag, codes) in enumerate(bags):
                new_patches[row, :len(bag)] = bag
                new_mask[row, :len(bag)] = True
                new_slides[row, :len(bag)] = codes
            batch["patches"] = torch.from_numpy(new_patches)
            batch["patch_mask"] = torch.from_numpy(new_mask)
            batch["slide_ids"] = torch.from_numpy(new_slides)
            batch.pop("coordinates", None)
            batch.pop("coordinate_present", None)
            yield batch


def build_model(checkpoint, data, device: str):
    manifest = checkpoint.get("manifest", {})
    dims = manifest.get("modal_dims", {})
    cfg = V2ModelConfig(**manifest["model_config"])
    head = (manifest.get("run_configuration") or {}).get("programme_head_dim") or data._v2_programmes.shape[1]
    model = TumorStateV2(cfg, clinical_dim=int(dims.get("clinical", 0)) or None,
                         snv_dim=int(dims.get("snv", 0)) or None, cnv_dim=int(dims.get("cnv", 0)) or None,
                         programme_dim=int(head)).to(device)
    model.load_state_dict(checkpoint["model"])
    return model.eval(), manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--data-config", required=True)
    parser.add_argument("--split-file", required=True)
    parser.add_argument("--programme-targets", default="")
    parser.add_argument("--normal-staging", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--levels", default="0.0,0.10,0.20,0.40,0.60")
    parser.add_argument("--token-budget", type=int, default=16384)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    levels = [float(v) for v in args.levels.split(",")]
    checkpoint = torch.load(args.checkpoint, map_location=args.device, weights_only=False)
    trained_split = (checkpoint.get("manifest", {}).get("split_manifest") or {})
    if trained_split.get("split_digest") and trained_split["split_digest"] != split_file_digest(args.split_file):
        raise ValueError("--split-file does not match the split this checkpoint trained on")
    data = load_bio_query_data(args.data_config, args.split_file, wsi_mode="hoptimus_patch")
    if (trained_split.get("cohort_restriction") or {}).get("enabled"):
        data, _ = restrict_cohort_to_split(data, args.split_file)
    manifest = checkpoint.get("manifest", {})
    fit_mask = (np.asarray(data.split).astype(str) != "test"
                if manifest.get("fit_population") == "development_train_val"
                else np.asarray(data.split).astype(str) == "train")
    if (manifest.get("programme_targets", {}) or {}).get("source") == "hallmark":
        attach_v2_targets(data, fit_mask)
    else:
        attach_external_programme_targets(data, args.programme_targets, fit_mask)
    if int((manifest.get("modal_dims") or {}).get("clinical", 0)) > 0:
        _standardize_clinical(data)
    model, manifest = build_model(checkpoint, data, args.device)

    pool, slides = load_normal_pool(Path(args.normal_staging))
    print(f"[pool] {len(slides)} normal slides, {len(pool)} patches", flush=True)
    trained = tuple(manifest.get("trained_states", ("wsi_identity", "rna_identity")))
    widths = {name: 256 for name in ("wsi_identity", "rna_identity", "wsi_biology", "rna_biology",
                                     "full_identity", "full_biology", "full_patient")}
    n = len(data.patient_ids)

    for level in levels:
        arrays = {key: np.zeros((n, width), np.float32) for key, width in widths.items() if key in trained}
        loader = DilutedHoptimusBatches(data, np.arange(n), args.token_budget, 917, shuffle=False,
                                        include_clinical=int((manifest.get("modal_dims") or {}).get("clinical", 0)) > 0,
                                        pool=pool, donor_lookup=lambda patient, cancer: slides,
                                        level=level, seed=args.seed)
        with torch.no_grad():
            for batch in loader:
                indices = batch.pop("indices").numpy()
                batch = {k: v.to(args.device, non_blocking=True) for k, v in batch.items()}
                views = {"wsi": model(batch, "wsi"), "rna": model(batch, "rna"), "full": model(batch, "full")}
                values = {"wsi_identity": views["wsi"]["z_identity"], "rna_identity": views["rna"]["z_identity"],
                          "wsi_biology": views["wsi"]["z_biology"], "rna_biology": views["rna"]["z_biology"],
                          "full_identity": views["full"]["z_identity"], "full_biology": views["full"]["z_biology"],
                          "full_patient": views["full"]["z_patient"]}
                for key, value in values.items():
                    if key in arrays:
                        arrays[key][indices] = value.cpu().float().numpy()
        output = Path(f"{args.output_prefix}_dilution_{int(round(level * 100)):03d}.npz")
        output.parent.mkdir(parents=True, exist_ok=True)
        level_manifest = dict(manifest)
        level_manifest["dilution"] = {
            "level": level, "definition": "fraction_of_final_bag_that_is_non_tumour",
            "normal_patches_drawn": int(sum(loader.drawn.values())),
            "normal_pool_slides": int(len(slides)), "normal_pool_patches": int(len(pool)),
            "slide_id_assignment": "existing_patient_slides_weighted_by_tumour_count",
            "seed": args.seed, "checkpoint": str(Path(args.checkpoint).resolve()),
            "checkpoint_sha256": sha256(Path(args.checkpoint).read_bytes()).hexdigest()}
        np.savez_compressed(output, patient_ids=np.asarray(data.patient_ids), split=np.asarray(data.split),
                            cancers=np.asarray(data.cancers),
                            trained_states=np.asarray(trainable_state_names(*trained)),
                            artifact_version=np.asarray(int(manifest.get("artifact_version", 4))),
                            manifest_json=np.asarray(json.dumps(level_manifest, sort_keys=True, default=str)),
                            **arrays)
        print(f"[ok] level {level:g} -> {output} (normal patches {sum(loader.drawn.values())})", flush=True)


if __name__ == "__main__":
    main()
