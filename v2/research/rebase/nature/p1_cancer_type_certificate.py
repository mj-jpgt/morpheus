"""P1's cancer-type confound certificate, on the Phase 1 cohort, written to an artifact.

Why this file exists
--------------------
``PHASE1_RESULT.md`` published a cancer-type balanced accuracy pair --
``0.463 raw -> 0.035 adjusted``, chance 0.048, n = 2,530 -- in prose, from a session
probe that never persisted. There is no run output anywhere on persistent storage
recording it, ``run_calibra`` has never emitted such a metric, and the canonical
estimator (``confound_certificate.lda_oof_balanced_accuracy``) was written on
2026-08-02, two days *after* the pair was committed. So the pair could be neither
traced nor reproduced, and was withdrawn.

This module is the repair: the same check, run from the same canonical functions the
rest of the project uses, against a hash-pinned artifact, writing its answer to disk.
A number that exists only in prose is the defect; a number with a path and a digest is
the fix.

What it measures
----------------
On the cohort ``run_calibra`` itself selects (``split == "test"`` in the artifact,
intersected with the frozen RNA target table), for every declared representation state:

* **raw** -- joint out-of-fold balanced accuracy for cancer type over all axes at once;
* **adjusted** -- the same statistic after ``cross_fitted_residuals`` against the
  99-column ``cancer`` + pooled-``tss`` design, i.e. the identical adjustment CALIBRA
  applies before it measures the channel;
* **chance** -- ``1 / n_classes``.

Every statistic is imported. ``lda_oof_balanced_accuracy``, ``nearest_class_mean_oof``,
``_encode`` and ``_stratified_folds`` come from ``v2/calibra/confound_certificate.py``;
``confound_design``, ``cross_fitted_residuals`` and ``pooled_tissue_source_site`` from
``v2/calibra/residualise.py``; the nonlinear probes from
``v2/calibra/nonlinear_confound_probe.py``. Nothing is reimplemented here. Re-deriving a
statistic inline where an import belonged is the single most repeated defect on this
project, and it is what produced the unattributable pair in the first place.

The nonlinear probes under ``--secondary`` are **diagnostic, not certification**. They
are reported because a mean-based certificate is a statement about the first moment
only, and because they bound what estimator the vanished July probe could have been.
They are not evidence of provenance: an estimator that happens to land near a published
number does not thereby become the estimator that produced it.

Example
-------
    PYTHONPATH=$WS python -m morpheus.v2.research.rebase.nature.p1_cancer_type_certificate \
        --artifact  .../runs_misc/calibra_run/artifacts/diagnostic_full_seed42.npz \
        --targets   .../runs_misc/calibra_run/artifacts/frozen_rna_targets.npz \
        --output    .../P1_CANCER_TYPE_CERTIFICATE.json --secondary
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))

from morpheus.v2.calibra.confound_certificate import (_encode, _stratified_folds,  # noqa: E402
                                                      lda_oof_balanced_accuracy,
                                                      nearest_class_mean_oof)
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,  # noqa: E402
                                             pooled_tissue_source_site)

__all__ = ["sha256_file", "select_cohort", "certificate"]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select_cohort(artifact, targets, partition: str = "test") -> np.ndarray:
    """``run_calibra.main``'s own cohort rule: the partition, aligned to the targets.

    Kept as one function so the cohort this certificate runs on cannot drift away from
    the cohort the channel measurement runs on. The n = 2,530 figure that made the
    published pair traceable at all is this mask's size.
    """
    patient_ids = np.asarray([str(p) for p in artifact["patient_ids"]])
    split = np.asarray([str(s) for s in artifact["split"]])
    target_index = {str(p): i for i, p in enumerate(targets["patient_ids"])}
    aligned = np.asarray([target_index.get(p, -1) for p in patient_ids])
    mask = np.ones(len(patient_ids), dtype=bool) if partition == "all" else (split == partition)
    return mask & (aligned >= 0)


def _standardise(block: np.ndarray) -> np.ndarray:
    """Per-axis standardisation, exactly as ``certify_axes`` does it before scoring."""
    return (block - block.mean(axis=0)) / np.maximum(block.std(axis=0), 1e-12)


def certificate(artifact_path: str, targets_path: str, *, partition: str = "test",
                min_site_count: int = 10, n_splits: int = 5, seed: int = 42,
                secondary: bool = False, n_jobs: int = 1) -> dict:
    artifact = np.load(artifact_path, allow_pickle=True)
    targets = np.load(targets_path, allow_pickle=True)

    cancers = np.asarray([str(c) for c in artifact["cancers"]])
    patient_ids = np.asarray([str(p) for p in artifact["patient_ids"]])
    split = np.asarray([str(s) for s in artifact["split"]])
    mask = select_cohort(artifact, targets, partition)

    site, frequent = pooled_tissue_source_site(patient_ids[mask], min_site_count=min_site_count)
    design = confound_design(pd.DataFrame({"cancer": cancers[mask], "tss": site}), ["cancer", "tss"])
    classes, class_index = _encode(cancers[mask])
    n_classes = int(len(classes))
    folds = _stratified_folds(class_index, n_splits, seed)

    report = {
        "artifact": {"path": str(artifact_path), "sha256": sha256_file(artifact_path)},
        "targets": {"path": str(targets_path), "sha256": sha256_file(targets_path)},
        "cohort": {
            "partition": partition,
            "n_patients": int(mask.sum()),
            "n_classes": n_classes,
            "chance_rate": 1.0 / n_classes,
            "n_confound_columns": int(design.shape[1]),
            "n_sites_kept": int(len(frequent)),
            "min_site_count": int(min_site_count),
            "split_counts": {k: int((split == k).sum()) for k in sorted(set(split.tolist()))},
        },
        "estimator": {
            "primary": "confound_certificate.lda_oof_balanced_accuracy",
            "per_axis": "confound_certificate.nearest_class_mean_oof",
            "n_splits": int(n_splits), "seed": int(seed), "shrinkage": 0.1,
            "standardisation": "per-axis, as certify_axes does it",
            "adjustment": ("residualise.cross_fitted_residuals against a cancer + pooled-tss "
                           "design, ridge alpha=1.0, 5 folds"),
        },
        "states": {},
    }

    for state in sorted(str(s) for s in artifact["trained_states"]):
        x = np.asarray(artifact[state], dtype=np.float64)[mask]
        adjusted = cross_fitted_residuals(x, design, seed=seed)
        raw_scaled, adjusted_scaled = _standardise(x), _standardise(adjusted)
        per_axis_raw = nearest_class_mean_oof(raw_scaled, class_index, n_classes,
                                              n_splits=n_splits, seed=seed, folds=folds)
        per_axis_adjusted = nearest_class_mean_oof(adjusted_scaled, class_index, n_classes,
                                                   n_splits=n_splits, seed=seed, folds=folds)
        row = {
            "n_axes": int(x.shape[1]),
            "raw_joint_lda": lda_oof_balanced_accuracy(raw_scaled, class_index, n_classes,
                                                       n_splits=n_splits, seed=seed, folds=folds),
            "adjusted_joint_lda": lda_oof_balanced_accuracy(adjusted_scaled, class_index, n_classes,
                                                            n_splits=n_splits, seed=seed, folds=folds),
            "raw_joint_lda_unstandardised": lda_oof_balanced_accuracy(
                x, class_index, n_classes, n_splits=n_splits, seed=seed, folds=folds),
            "adjusted_joint_lda_unstandardised": lda_oof_balanced_accuracy(
                adjusted, class_index, n_classes, n_splits=n_splits, seed=seed, folds=folds),
            "raw_per_axis_max": float(np.max(per_axis_raw)),
            "raw_per_axis_median": float(np.median(per_axis_raw)),
            "adjusted_per_axis_max": float(np.max(per_axis_adjusted)),
            "adjusted_per_axis_median": float(np.median(per_axis_adjusted)),
        }
        if secondary:
            from morpheus.v2.calibra.nonlinear_confound_probe import (
                forest_balanced_accuracy_oof, knn_balanced_accuracy_oof,
                rbf_svm_balanced_accuracy_oof)
            for label, block in (("raw", raw_scaled), ("adjusted", adjusted_scaled)):
                row[f"{label}_knn15"] = knn_balanced_accuracy_oof(
                    block, class_index, n_classes, k=15, n_splits=n_splits, seed=seed, folds=folds)
                row[f"{label}_knn15_prior_corrected"] = knn_balanced_accuracy_oof(
                    block, class_index, n_classes, k=15, n_splits=n_splits, seed=seed, folds=folds,
                    prior_correction=True)
                row[f"{label}_forest"] = forest_balanced_accuracy_oof(
                    block, class_index, n_classes, n_splits=n_splits, seed=seed, folds=folds,
                    n_jobs=n_jobs)
                row[f"{label}_rbf_svm"] = rbf_svm_balanced_accuracy_oof(
                    block, class_index, n_classes, n_splits=n_splits, seed=seed, folds=folds)
        report["states"][state] = row
        print(f"[{state:14s}] raw={row['raw_joint_lda']:.6f} adjusted={row['adjusted_joint_lda']:.6f} "
              f"chance={1.0 / n_classes:.6f}", flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--targets", required=True, help="frozen_rna_targets.npz")
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition", default="test", choices=("test", "val", "all"))
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--secondary", action="store_true",
                        help="also run the k-NN / forest / RBF-SVM probes. Diagnostic only: these "
                             "are not the certificate and cannot establish provenance.")
    args = parser.parse_args()

    report = certificate(args.artifact, args.targets, partition=args.partition,
                         min_site_count=args.min_site_count, n_splits=args.n_splits,
                         seed=args.seed, secondary=args.secondary, n_jobs=args.n_jobs)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(f"[done] -> {output}")


if __name__ == "__main__":
    main()
