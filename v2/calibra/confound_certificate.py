"""Confound certificate — must-FAIL control 1: axes may not predict site/scanner/batch.

This is Track 1's T1.3 and, verbatim, condition 4 of the five-point certification
rule in ``MULTIMODAL_EXPANSION.md`` §1: *an axis may only be exposed if it FAILS
to predict tissue source site*. It closes G3.5, which
``e0_basis_transfer.py`` currently records as ``unavailable_no_site_labels`` — a
gap that was never a data problem, because a TCGA barcode carries its TSS code in
the second field and ``residualise.pooled_tissue_source_site`` already derives it.

Why the null permutes labels WITHIN CANCER
------------------------------------------
Tissue source sites are not exchangeable across cancer types: a site contributes
patients to a handful of cancers, so cancer type alone predicts site well above
chance. A globally permuted label null destroys that association and therefore
sets an unreasonably low bar — every axis carrying ordinary lineage information
would "predict site" and be refused certification, and the certificate would
refuse everything without distinguishing a site code from a biology axis.

Permuting the site label *within cancer type* holds the legitimate cancer↔site
association fixed and asks the only question certification cares about: does this
axis carry site information **beyond what cancer type already explains**? That is
the same convention ``calibration.permutation_null`` uses for the pairing null, so
the two nulls in this instrument agree about what "chance" means.

Why balanced accuracy, and why nearest-class-mean
-------------------------------------------------
Site classes are wildly imbalanced, so plain accuracy is maximised by always
predicting the largest site and would certify a constant axis. Balanced accuracy
weights every site equally and has an exactly known chance rate, ``1/n_classes``.

The per-axis classifier is a nearest-class-mean rule on the single standardised
axis, fit out-of-fold. For a one-dimensional feature under uniform class priors
this *is* linear discriminant analysis, it has no hyperparameters to tune (so
there is no protocol in which a defender could weaken it), and it is cheap enough
to support the >=1,000-permutation null the pass criterion requires.

The joint test is not optional
------------------------------
A per-axis test is the certification unit, but it is also the weakest possible
version of the question: site information can be spread across many coordinates
and be invisible in every one of them individually. ``certify_axes`` therefore
also runs a *joint* LDA over all axes at once. A representation that passes
per-axis and fails jointly has not been shown to be site-free — it has been shown
that nobody looked with enough resolution — and the joint row is reported with
equal prominence for exactly that reason.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site

__all__ = ["balanced_accuracy", "nearest_class_mean_oof", "lda_oof_balanced_accuracy",
           "within_stratum_permutations", "certify_axes"]


def _encode(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    classes, index = np.unique(np.asarray(labels).astype(str), return_inverse=True)
    return classes, index.astype(np.int64)


def balanced_accuracy(true_index: np.ndarray, predicted_index: np.ndarray, n_classes: int) -> float:
    """Mean per-class recall over the classes actually present in ``true_index``."""
    true_index = np.asarray(true_index, dtype=np.int64)
    predicted_index = np.asarray(predicted_index, dtype=np.int64)
    correct = np.bincount(true_index, weights=(predicted_index == true_index).astype(np.float64),
                          minlength=n_classes)
    total = np.bincount(true_index, minlength=n_classes).astype(np.float64)
    present = total > 0
    if not present.any():
        return float("nan")
    return float(np.mean(correct[present] / total[present]))


def _stratified_folds(class_index: np.ndarray, n_splits: int, seed: int) -> list[np.ndarray]:
    """Deterministic class-stratified fold assignment.

    Rolled by hand rather than taken from sklearn because the permutation null
    re-folds thousands of times and must do so on permuted labels *identically*;
    keeping the rule in one visible place is worth more here than the import.
    """
    rng = np.random.default_rng(seed)
    fold = np.empty(len(class_index), dtype=np.int64)
    for value in np.unique(class_index):
        rows = np.flatnonzero(class_index == value)
        rows = rng.permutation(rows)
        fold[rows] = np.arange(len(rows)) % n_splits
    return [np.flatnonzero(fold == f) for f in range(n_splits)]


def nearest_class_mean_oof(features: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                           n_splits: int = 5, seed: int = 42,
                           folds: list[np.ndarray] | None = None) -> np.ndarray:
    """Out-of-fold balanced accuracy for EVERY column of ``features``, vectorised.

    ``features`` is ``[patient, axis]``. Each axis is treated as its own
    one-dimensional classifier; the return value is one balanced accuracy per
    axis. Vectorising over axes is what makes a per-axis 1,000-permutation null
    affordable.
    """
    features = np.asarray(features, dtype=np.float64)
    class_index = np.asarray(class_index, dtype=np.int64)
    n, n_axes = features.shape
    folds = folds if folds is not None else _stratified_folds(class_index, n_splits, seed)
    predicted = np.empty((n, n_axes), dtype=np.int64)
    for test_rows in folds:
        train_rows = np.setdiff1d(np.arange(n), test_rows, assume_unique=False)
        onehot = np.zeros((n_classes, len(train_rows)))
        onehot[class_index[train_rows], np.arange(len(train_rows))] = 1.0
        counts = onehot.sum(axis=1)
        # Class means per axis. Empty training classes are pushed to +inf so the
        # nearest-mean rule can never assign a class it has not seen.
        means = np.divide(onehot @ features[train_rows], np.maximum(counts, 1)[:, None])
        means[counts == 0] = np.inf
        # |x - m| over [test, class, axis]; chunked so the broadcast stays bounded.
        block = max(1, int(4_000_000 // max(n_classes * n_axes, 1)))
        for start in range(0, len(test_rows), block):
            rows = test_rows[start:start + block]
            distance = np.abs(features[rows][:, None, :] - means[None, :, :])
            predicted[rows] = np.argmin(distance, axis=1)
    return np.asarray([balanced_accuracy(class_index, predicted[:, a], n_classes)
                       for a in range(n_axes)], dtype=np.float64)


def lda_oof_balanced_accuracy(features: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                              n_splits: int = 5, seed: int = 42, shrinkage: float = 0.1,
                              folds: list[np.ndarray] | None = None) -> float:
    """JOINT out-of-fold balanced accuracy: all axes at once, shrunk-covariance LDA.

    Uniform class priors, because the statistic being reported is balanced
    accuracy; empirical priors would optimise a different loss and quietly make
    the test easier. Shrinkage is fixed, not tuned, so the null and the observed
    statistic are produced by exactly the same estimator.
    """
    features = np.asarray(features, dtype=np.float64)
    class_index = np.asarray(class_index, dtype=np.int64)
    n, d = features.shape
    folds = folds if folds is not None else _stratified_folds(class_index, n_splits, seed)
    predicted = np.empty(n, dtype=np.int64)
    for test_rows in folds:
        train_rows = np.setdiff1d(np.arange(n), test_rows, assume_unique=False)
        onehot = np.zeros((n_classes, len(train_rows)))
        onehot[class_index[train_rows], np.arange(len(train_rows))] = 1.0
        counts = onehot.sum(axis=1)
        means = np.divide(onehot @ features[train_rows], np.maximum(counts, 1)[:, None])
        centred = features[train_rows] - means[class_index[train_rows]]
        covariance = (centred.T @ centred) / max(len(train_rows) - int((counts > 0).sum()), 1)
        covariance.flat[:: d + 1] += shrinkage * float(np.trace(covariance)) / d + 1e-8
        weights = np.linalg.solve(covariance, means.T).T                    # [class, dim]
        bias = -0.5 * np.einsum("kd,kd->k", means, weights)
        score = features[test_rows] @ weights.T + bias
        score[:, counts == 0] = -np.inf
        predicted[test_rows] = np.argmax(score, axis=1)
    return balanced_accuracy(class_index, predicted, n_classes)


def within_stratum_permutations(class_index: np.ndarray, strata: np.ndarray, *,
                                n_permutations: int, seed: int) -> list[np.ndarray]:
    """Permute site labels WITHIN cancer type — see the module docstring."""
    class_index = np.asarray(class_index, dtype=np.int64)
    strata = np.asarray(strata).astype(str)
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(strata == level) for level in np.unique(strata)]
    permutations = []
    for _ in range(n_permutations):
        permuted = class_index.copy()
        for rows in groups:
            permuted[rows] = class_index[rng.permutation(rows)]
        permutations.append(permuted)
    return permutations


def _bootstrap_ci(features_axis: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                  n_boot: int, seed: int, n_splits: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(class_index)
    values = []
    for _ in range(n_boot):
        rows = rng.integers(0, n, size=n)
        if len(np.unique(class_index[rows])) < 2:
            continue
        classes_here, local = np.unique(class_index[rows], return_inverse=True)
        values.append(float(nearest_class_mean_oof(features_axis[rows][:, None], local,
                                                   len(classes_here), n_splits=n_splits,
                                                   seed=seed)[0]))
    if not values:
        return float("nan"), float("nan")
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def certify_axes(features: np.ndarray, patient_ids: np.ndarray, cancers: np.ndarray, *,
                 min_site_count: int = 10, n_splits: int = 5, n_permutations: int = 1000,
                 seed: int = 42, n_boot: int = 200, n_boot_axes: int = 8,
                 residualise: bool = False, n_jobs: int = 1) -> dict:
    """Run the confound certificate for one representation state.

    Pass criterion, stated in advance so "it failed" is falsifiable:
      * per axis -- out-of-fold balanced accuracy for pooled TSS must NOT exceed
        the 95th percentile of the within-cancer label-permutation null, and
      * the axis bootstrap CI must include the chance rate ``1/n_classes``.
    An axis breaching either bar is **not certified** and is named in
    ``breaching_axes``; it is never silently dropped.
    """
    features = np.asarray(features, dtype=np.float64)
    site, frequent = pooled_tissue_source_site(np.asarray(patient_ids).astype(str),
                                               min_site_count=min_site_count)
    classes, class_index = _encode(site)
    n_classes = int(len(classes))
    cancers = np.asarray(cancers).astype(str)
    if n_classes < 2 or len(features) != len(class_index):
        return {"status": "unavailable_insufficient_site_labels", "n_classes": n_classes}
    if residualise:
        # The SAME adjustment CALIBRA applies before it measures the channel:
        # cancer + pooled TSS, cross-fitted. Certifying the raw state answers
        # P4's question ("may this axis be exposed?"); certifying the adjusted
        # state answers P1's ("does our own adjustment actually discharge the
        # site confound, or only appear to?"). They are different questions and
        # a battery that reports only one of them is not a battery.
        import pandas as pd
        design = confound_design(pd.DataFrame({"cancer": cancers, "tss": site}), ["cancer", "tss"])
        features = cross_fitted_residuals(features, design, seed=seed)
    # Standardising per axis is what makes the nearest-class-mean rule comparable
    # across axes of different scale; it cannot change any single axis's accuracy.
    scaled = (features - features.mean(axis=0)) / np.maximum(features.std(axis=0), 1e-12)
    folds = _stratified_folds(class_index, n_splits, seed)
    observed = nearest_class_mean_oof(scaled, class_index, n_classes, n_splits=n_splits,
                                      seed=seed, folds=folds)
    joint_observed = lda_oof_balanced_accuracy(scaled, class_index, n_classes,
                                               n_splits=n_splits, seed=seed, folds=folds)
    permutations = within_stratum_permutations(class_index, cancers,
                                               n_permutations=n_permutations, seed=seed)

    from .calibration import _map

    def _one_axis_null(permuted):
        return nearest_class_mean_oof(scaled, permuted, n_classes, n_splits=n_splits,
                                      seed=seed, folds=folds)

    def _one_joint_null(permuted):
        return lda_oof_balanced_accuracy(scaled, permuted, n_classes, n_splits=n_splits,
                                         seed=seed, folds=folds)

    null = np.asarray(_map(_one_axis_null, [(p,) for p in permutations], n_jobs))   # [perm, axis]
    joint_null = np.asarray(_map(_one_joint_null, [(p,) for p in permutations], n_jobs))
    null_p95 = np.percentile(null, 95, axis=0)
    null_median = np.median(null, axis=0)
    # Floored at 1/(n+1) exactly as G4.5 requires; the resolution is reported so a
    # p at the floor can never be read as a smaller number than the design supports.
    axis_p = (np.sum(null >= observed[None, :], axis=0) + 1.0) / (n_permutations + 1.0)
    chance = 1.0 / n_classes
    exceeds_null = observed > null_p95

    # The bootstrap is the expensive half, so it is spent where it decides
    # something: every axis that already breached the permutation bar, plus the
    # strongest few that did not, so a near-miss cannot hide behind a missing CI.
    order = np.argsort(-observed)
    boot_axes = sorted(set(np.flatnonzero(exceeds_null).tolist())
                       | set(order[:max(0, int(n_boot_axes))].tolist()))
    ci_low = np.full(len(observed), np.nan)
    ci_high = np.full(len(observed), np.nan)
    for axis in boot_axes:
        low, high = _bootstrap_ci(scaled[:, axis], class_index, n_classes,
                                  n_boot=n_boot, seed=seed + axis, n_splits=n_splits)
        ci_low[axis], ci_high[axis] = low, high
    ci_excludes_chance = np.zeros(len(observed), dtype=bool)
    measured = np.isfinite(ci_low) & np.isfinite(ci_high)
    ci_excludes_chance[measured] = ~((ci_low[measured] <= chance) & (chance <= ci_high[measured]))

    breaching = np.flatnonzero(exceeds_null & (ci_excludes_chance | ~measured))
    joint_p = float((np.sum(joint_null >= joint_observed) + 1.0) / (n_permutations + 1.0))
    return {
        "status": "scored",
        "residualised": bool(residualise),
        "adjustment": "cancer+pooled_tss_cross_fitted" if residualise else "none_raw_state",
        "n_patients": int(len(features)), "n_axes": int(features.shape[1]),
        "n_classes": n_classes, "n_sites_kept": int(len(frequent)),
        "min_site_count": int(min_site_count), "chance_rate": float(chance),
        "n_permutations": int(n_permutations), "permutation_resolution": 1.0 / (n_permutations + 1.0),
        "permutation_strata": "cancer_type",
        "observed_max": float(np.max(observed)), "observed_median": float(np.median(observed)),
        "null_p95_median": float(np.median(null_p95)), "null_median_median": float(np.median(null_median)),
        "n_axes_exceeding_null_p95": int(exceeds_null.sum()),
        "fraction_axes_exceeding_null_p95": float(exceeds_null.mean()),
        "min_axis_permutation_p": float(axis_p.min()),
        "n_axes_p_below_0p05": int((axis_p < 0.05).sum()),
        "n_breaching_axes": int(len(breaching)),
        "breaching_axes": [{"axis": int(a), "balanced_accuracy": float(observed[a]),
                            "null_p95": float(null_p95[a]), "null_median": float(null_median[a]),
                            "permutation_p": float(axis_p[a]), "ci95": [float(ci_low[a]), float(ci_high[a])],
                            "chance_rate": float(chance)}
                           for a in breaching[np.argsort(-observed[breaching])][:25].tolist()],
        "certified": bool(len(breaching) == 0),
        "joint_lda_balanced_accuracy": float(joint_observed),
        "joint_null_median": float(np.median(joint_null)), "joint_null_p95": float(np.percentile(joint_null, 95)),
        "joint_permutation_p": joint_p,
        "joint_certified": bool(joint_observed <= float(np.percentile(joint_null, 95))),
        "per_axis_balanced_accuracy": observed.tolist(),
        "per_axis_null_p95": null_p95.tolist(),
        "per_axis_permutation_p": axis_p.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition", default="test", choices=("test", "val", "all"))
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--n-boot", type=int, default=200)
    parser.add_argument("--n-boot-axes", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--residualise", action="store_true",
                        help="certify the state AFTER CALIBRA's own cancer+pooled-TSS adjustment")
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict] = {}
    rows: list[dict] = []
    for artifact_path in args.artifacts:
        path = Path(artifact_path)
        method = path.stem
        raw = np.load(path, allow_pickle=True)
        declared = {str(s) for s in raw["trained_states"]} if "trained_states" in raw.files else set()
        patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
        cancers = np.asarray([str(c) for c in raw["cancers"]])
        split = np.asarray([str(s) for s in raw["split"]])
        mask = np.ones(len(patient_ids), dtype=bool) if args.partition == "all" else (split == args.partition)
        for state in sorted(declared):
            key = f"{method}::{state}"
            if state not in raw.files:
                results[key] = {"status": "unavailable_state_declared_but_absent"}
                continue
            result = certify_axes(np.asarray(raw[state], dtype=np.float64)[mask],
                                  patient_ids[mask], cancers[mask],
                                  min_site_count=args.min_site_count, n_splits=args.n_splits,
                                  n_permutations=args.n_permutations, seed=args.seed,
                                  n_boot=args.n_boot, n_boot_axes=args.n_boot_axes,
                                  residualise=args.residualise, n_jobs=args.n_jobs)
            results[key] = result
            for metric in ("observed_max", "observed_median", "chance_rate", "null_p95_median",
                           "n_axes_exceeding_null_p95", "fraction_axes_exceeding_null_p95",
                           "n_breaching_axes", "joint_lda_balanced_accuracy", "joint_null_p95",
                           "joint_permutation_p"):
                rows.append({"method": method, "representation_state": state, "task": "confound_certificate",
                             "target": "tissue_source_site", "metric": metric,
                             "value": float(result.get(metric, float("nan"))), "note": result["status"]})
            print(f"[{key}] axes={result.get('n_axes')} breach={result.get('n_breaching_axes')} "
                  f"per_axis_max={result.get('observed_max', float('nan')):.4f} "
                  f"chance={result.get('chance_rate', float('nan')):.4f} "
                  f"joint={result.get('joint_lda_balanced_accuracy', float('nan')):.4f} "
                  f"joint_null_p95={result.get('joint_null_p95', float('nan')):.4f} "
                  f"{'CERTIFIED' if result.get('certified') else 'NOT-CERTIFIED'}"
                  f"/{'joint-OK' if result.get('joint_certified') else 'JOINT-FAIL'}", flush=True)
    (output / "confound_certificate.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    import pandas as pd
    pd.DataFrame(rows).to_csv(output / "task_rows.csv", index=False)
    print(f"[done] {len(rows)} rows -> {output}", flush=True)


if __name__ == "__main__":
    main()
