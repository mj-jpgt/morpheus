"""Export frozen V2 representations in the shared patient-artifact schema."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import torch
from morpheus.src.training.train_bio_query_former import load_bio_query_data
from .preflight import restrict_cohort_to_split, split_file_digest
from .model import TumorStateV2, V2ModelConfig
from .runner import (UncappedHoptimusBatches, _attach_numeric_table,
                     _standardize_clinical, attach_mlp_clip_anchor, attach_v2_targets,
                     attach_external_programme_targets)
from .contracts import trainable_state_names


def _programme_head_dim(manifest: dict, data) -> int:
    """Head width the CHECKPOINT was built with, not the raw target width.

    Training widens the programme head to --programme-head-dim (default 256) and
    pads the target; loading with the raw width (50 Hallmark / 128 PBS) raises
    "size mismatch for programme_mean.weight" and every export fails after the
    GPU time is already spent.
    """
    recorded = (manifest.get("run_configuration") or {}).get("programme_head_dim")
    if recorded:
        return int(recorded)
    return int(getattr(data, "_v2_programme_head_dim", 0) or data._v2_programmes.shape[1])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data-config", required=True); p.add_argument("--split-file", required=True)
    p.add_argument("--checkpoint", required=True); p.add_argument("--output", required=True)
    p.add_argument("--token-budget", type=int, default=32768); p.add_argument("--hidden-dim", type=int, default=512)
    p.add_argument("--layers", type=int, default=4); p.add_argument("--heads", type=int, default=8); p.add_argument("--device", default="cuda")
    p.add_argument("--semantic-dim", type=int, default=0)
    p.add_argument("--mlp-clip-anchor", default="")
    p.add_argument("--programme-targets", default="",
                   help="immutable PBS target artifact required to export a PBS-supervised checkpoint")
    p.add_argument("--snv-features", default=""); p.add_argument("--cnv-features", default="")
    p.add_argument("--diagnostics-output", default="", help="optional JSON sidecar for anchor/slot diagnostics")
    a = p.parse_args()
    checkpoint = torch.load(a.checkpoint, map_location=a.device, weights_only=False)
    manifest = checkpoint.get("manifest", {})
    # The cohort comes from the CHECKPOINT, never from a flag. Exporting against a
    # different split than the model trained on scores a different set of patients,
    # and a paired comparison would not notice because BOTH arms would be wrong.
    # The PBS branch below catches it incidentally via its target fit digest; the
    # Hallmark branch -- which is the whole D1 path -- has no other check.
    trained_split = manifest.get("split_manifest") or {}
    if trained_split.get("split_digest") and trained_split["split_digest"] != split_file_digest(a.split_file):
        raise ValueError(f"--split-file does not match the split this checkpoint was trained on "
                         f"(expected digest {trained_split['split_digest'][:16]}...)")
    data = load_bio_query_data(a.data_config, a.split_file, wsi_mode="hoptimus_patch")
    if (trained_split.get("cohort_restriction") or {}).get("enabled"):
        data, _ = restrict_cohort_to_split(data, a.split_file)
    fit_mask = (np.asarray(data.split).astype(str) != "test" if manifest.get("fit_population") == "development_train_val"
                else np.asarray(data.split).astype(str) == "train")
    programme_manifest = manifest.get("programme_targets", {})
    if programme_manifest.get("source") == "hallmark":
        if a.programme_targets:
            raise ValueError("Hallmark-supervised checkpoint must not be exported with a PBS target artifact")
        attach_v2_targets(data, fit_mask)
    else:
        if not a.programme_targets:
            raise ValueError("PBS-supervised checkpoint export requires --programme-targets")
        attached = attach_external_programme_targets(data, a.programme_targets, fit_mask)
        if int(attached["target_dimension"]) != int(programme_manifest.get("target_dimension", -1)):
            raise ValueError("PBS target width differs from the training checkpoint")
        expected_pbs = programme_manifest.get("manifest")
        if not isinstance(expected_pbs, dict):
            raise ValueError("PBS-supervised checkpoint lacks immutable target provenance")
        actual_pbs = attached["manifest"]
        immutable = ("target_kind", "patient_id_digest", "split_digest", "overlap_gene_digest",
                     "fit_patient_id_digest", "scores_sha256", "singular_values_sha256",
                     "gene_basis_sha256", "atom_coordinates_sha256", "atom_id_digest")
        changed = [key for key in immutable if expected_pbs.get(key) != actual_pbs.get(key)]
        if changed:
            raise ValueError(f"PBS export target provenance differs from training checkpoint: {changed}")
    trained = tuple(manifest.get("trained_states", ("wsi_identity", "rna_identity")))
    semantic_dim = int(manifest.get("semantic_dim", a.semantic_dim)); dims = manifest.get("modal_dims", {})
    use_clinical = int(dims.get("clinical", 0)) > 0
    if use_clinical:
        _standardize_clinical(data)
    for name, path in (("snv", a.snv_features), ("cnv", a.cnv_features)):
        expected = int(dims.get(name, 0))
        if expected and not path:
            raise ValueError(f"checkpoint uses {name}; provide --{name}-features for export")
        if path:
            _attach_numeric_table(data, path, name)
            if expected and getattr(data, f"_v2_{name}").shape[1] != expected:
                raise ValueError(f"{name} feature width differs from checkpoint")
    if bool(manifest.get("mlp_clip_anchor", False)):
        if not a.mlp_clip_anchor:
            raise ValueError("anchored checkpoint export requires --mlp-clip-anchor")
        attach_mlp_clip_anchor(data, a.mlp_clip_anchor)
    if manifest.get("model_config"):
        cfg = V2ModelConfig(**manifest["model_config"])
    else:
        cfg = V2ModelConfig(rna_dim=int(data.rna.shape[1]), hidden_dim=a.hidden_dim, layers=a.layers, heads=a.heads,
                             semantic_dim=semantic_dim, use_mlp_clip_anchor=bool(manifest.get("mlp_clip_anchor", False)))
    model = TumorStateV2(cfg, clinical_dim=int(dims.get("clinical", data.clinical.shape[1])),
                         snv_dim=int(dims.get("snv", 0)) or None, cnv_dim=int(dims.get("cnv", 0)) or None,
                         programme_dim=_programme_head_dim(manifest, data)).to(a.device)
    model.load_state_dict(checkpoint["model"]); model.eval()
    widths = {"wsi_identity": 256, "rna_identity": 256, "wsi_biology": 256, "rna_biology": 256,
              "full_identity": 256, "full_biology": 256, "full_patient": 256}
    if semantic_dim:
        widths["semantic_wsi"] = semantic_dim
    n = len(data.patient_ids); arrays = {key: np.zeros((n, width), np.float32) for key, width in widths.items() if key in trained}
    diagnostics: dict[str, list[float]] = {}
    loader = UncappedHoptimusBatches(data, np.arange(n), a.token_budget, 917,
                                     include_clinical=use_clinical)
    with torch.no_grad():
        for batch in loader:
            indices = batch.pop("indices").numpy()
            batch = {k: v.to(a.device, non_blocking=True) for k,v in batch.items()}
            wsi, rna, full = model(batch,"wsi"), model(batch,"rna"), model(batch,"full")
            values = {"wsi_identity":wsi["z_identity"], "rna_identity":rna["z_identity"],
                      "wsi_biology":wsi["z_biology"], "rna_biology":rna["z_biology"],
                      "full_identity":full["z_identity"], "full_biology":full["z_biology"],
                      "full_patient":full["z_patient"]}
            if "z_semantic" in wsi:
                values["semantic_wsi"] = wsi["z_semantic"]
            for key,value in values.items():
                if key in arrays:
                    arrays[key][indices] = value.cpu().float().numpy()
            for key in ("anchor_residual_scale", "anchor_gate_mean", "anchor_correction_norm",
                        "local_slot_entropy", "coordinate_present_fraction"):
                if key in wsi:
                    diagnostics.setdefault(key, []).append(float(wsi[key].detach().float().mean().cpu()))
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.output, patient_ids=np.asarray(data.patient_ids), split=np.asarray(data.split), cancers=np.asarray(data.cancers),
                        trained_states=np.asarray(trainable_state_names(*trained)),
                        artifact_version=np.asarray(int(manifest.get("artifact_version", 4))),
                        manifest_json=np.asarray(json.dumps(manifest, sort_keys=True)), **arrays)
    sidecar = Path(a.diagnostics_output) if a.diagnostics_output else Path(f"{a.output}.diagnostics.json")
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(json.dumps({
        "artifact": str(Path(a.output).resolve()),
        "artifact_version": int(manifest.get("artifact_version", 4)),
        "source_manifest": manifest.get("source_manifest"),
        "teacher_mode": manifest.get("teacher_mode", "legacy_unknown"),
        "diagnostics_mean": {key: float(np.mean(values)) for key, values in diagnostics.items()},
    }, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__": main()
