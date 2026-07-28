"""Export fair raw H-Optimus, Ridge, and CCA baselines from one paired cohort.

Example:
``python -m morpheus.v2.export_baselines --data-config ... --split-file ...
--output-dir ...``.  The command opens the canonical H-Optimus token store,
uses every patch for each patient (no cap), and writes artifacts in exactly the
same patient/cancer/split order used by V2.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Protocol, Sequence

import numpy as np

from morpheus.src.training.train_bio_query_former import load_bio_query_data

from .baseline_exports import (export_cca_alignment_baseline,
                               export_raw_hoptimus_baseline,
                               export_ridge_alignment_baseline)
from .preflight import sha256_json, validate_runtime_split
from .provenance import source_manifest


class PatchStore(Protocol):
    def load_many_patient_tokens(self, patient_ids: Sequence[str], max_tokens: int | None, seed: int,
                                 slide_balanced: bool = True): ...


def uncapped_patient_mean_std(store: PatchStore, patient_ids: Sequence[str], *, read_batch: int = 8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute exact all-token per-patient mean/std with bounded host memory."""
    if read_batch < 1:
        raise ValueError("read_batch must be positive")
    means: list[np.ndarray] = []; stds: list[np.ndarray] = []; counts: list[int] = []
    for start in range(0, len(patient_ids), read_batch):
        identifiers = [str(value) for value in patient_ids[start:start + read_batch]]
        loaded = store.load_many_patient_tokens(identifiers, max_tokens=None, seed=0, slide_balanced=True)
        if len(loaded) != len(identifiers):
            raise RuntimeError("patch store returned a different number of patient bags")
        for patient, (tokens, _) in zip(identifiers, loaded):
            values = np.asarray(tokens, dtype=np.float32)
            if values.ndim != 2 or len(values) == 0 or not np.isfinite(values).all():
                raise RuntimeError(f"canonical H-Optimus store returned invalid uncapped tokens for {patient}")
            means.append(values.mean(axis=0, dtype=np.float64).astype(np.float32))
            stds.append(values.std(axis=0, dtype=np.float64).astype(np.float32))
            counts.append(int(len(values)))
    return np.stack(means), np.stack(stds), np.asarray(counts, dtype=np.int64)


def export_all_baselines(*, data_config: str, split_file: str, output_dir: str, components: int = 256,
                         cca_components: int = 128, patient_read_batch: int = 8,
                         expected_development_cancers: int = 11, expected_heldout_cancers: int = 22,
                         fit_development: bool = False) -> dict[str, object]:
    """Load exact paired data once and write all three fair baseline artifacts."""
    data = load_bio_query_data(data_config, split_file, wsi_mode="hoptimus_patch")
    split_manifest = validate_runtime_split(split_file, data.patient_ids, data.cancers, data.split,
                                            expected_development_cancers, expected_heldout_cancers)
    if data.hoptimus_store is None:
        raise RuntimeError("hoptimus_patch data load did not provide canonical patch store")
    store_info = data.hoptimus_store.validate()
    mean, std, counts = uncapped_patient_mean_std(data.hoptimus_store, data.patient_ids, read_batch=patient_read_batch)
    if np.any(counts < 1):
        raise RuntimeError("paired baseline cohort contains a patient without H-Optimus patches")
    destination = Path(output_dir); destination.mkdir(parents=True, exist_ok=True)
    provenance = {"data_config": str(Path(data_config)), "split_file": str(Path(split_file)), "split_digest": split_manifest["split_digest"],
                  "cohort_digest": split_manifest["cohort_digest"], "hoptimus_store": store_info,
                  "patch_count_digest": sha256_json(counts.tolist()), "patch_cap": None, "rna_source": "paired_bulkformer_embedding",
                  "fit_population": "development_train_val" if fit_development else "train_only_inner_fit"}
    common = {"patient_ids": data.patient_ids, "cancers": data.cancers, "split": data.split, "source_provenance": provenance}
    paths = {
        "raw_hoptimus": export_raw_hoptimus_baseline(destination / "raw_hoptimus_meanstd.npz", pooled_mean=mean, pooled_std=std, **common),
        "ridge": export_ridge_alignment_baseline(destination / "ridge_alignment.npz", wsi=np.concatenate([mean, std], axis=1), rna=data.rna, components=components, fit_development=fit_development, **common),
        "cca": export_cca_alignment_baseline(destination / "cca_alignment.npz", wsi=np.concatenate([mean, std], axis=1), rna=data.rna, components=cca_components, fit_development=fit_development, **common),
    }
    # Reuse exact complete-bag sufficient statistics for every stochastic
    # projection baseline.  This avoids six redundant reads of the ragged WSI
    # store and makes feature identity auditable across MLP-CLIP variants.
    pooled_cache = destination / "uncapped_all_patch_meanstd_cache.npz"
    np.savez_compressed(pooled_cache, patient_ids=np.asarray(data.patient_ids, dtype=str),
                        cancers=np.asarray(data.cancers, dtype=str), split=np.asarray(data.split, dtype=str),
                        pooled_mean=mean, pooled_std=std, patch_counts=counts,
                        source_provenance_json=np.asarray(json.dumps(provenance, sort_keys=True)))
    manifest = {"format_version": 1, "split": split_manifest, "store": store_info, "source_provenance": provenance,
                "n_patients": len(data.patient_ids), "token_count": {"min": int(counts.min()), "max": int(counts.max()), "total": int(counts.sum())},
                "artifacts": {name: str(path) for name, path in paths.items()},
                "uncapped_feature_cache": str(pooled_cache),
                "source": source_manifest(configuration={"data_config": data_config, "split_file": split_file,
                                                         "components": components, "cca_components": cca_components})}
    (destination / "baseline_export_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-config", required=True); parser.add_argument("--split-file", required=True); parser.add_argument("--output-dir", required=True)
    parser.add_argument("--components", type=int, default=256); parser.add_argument("--cca-components", type=int, default=128)
    parser.add_argument("--patient-read-batch", type=int, default=8)
    parser.add_argument("--fit-development", action="store_true", help="final refit on train+validation patients")
    parser.add_argument("--expected-development-cancers", type=int, default=11); parser.add_argument("--expected-heldout-cancers", type=int, default=22)
    args = parser.parse_args()
    print(json.dumps(export_all_baselines(data_config=args.data_config, split_file=args.split_file, output_dir=args.output_dir,
                                           components=args.components, cca_components=args.cca_components, patient_read_batch=args.patient_read_batch,
                                           expected_development_cancers=args.expected_development_cancers,
                                           expected_heldout_cancers=args.expected_heldout_cancers,
                                           fit_development=args.fit_development), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
