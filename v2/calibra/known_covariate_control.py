"""Positive control (b): a known-legible covariate must come back at its KNOWN strength.

``HANDOFF_GATES.md``'s governing rule is that a negative result is reportable
only if the positive control passed in the same run, on the same data, through
the same code path. CALIBRA already ships two of the three positive controls
§9 asks for: ``run_calibra --require-rna-positive-control`` (RNA-input states
predicting RNA-derived targets, circular by construction and therefore obliged
to be strong) and ``spike_recovery_curve`` (a synthetic spike above the floor,
recovered). This module is the third and the only one that is not self-referential:
a covariate whose morphological legibility was established by *other people, on
other cohorts, before we looked*.

The thing that makes it a control rather than a demo
----------------------------------------------------
"We recovered something" is not a pass. The expected strength is written down
**before the run**, from the published literature, with citations, and the
control passes only if the measured confidence interval overlaps it. A CI lying
entirely below the interval is under-recovery — the pipeline is destroying signal
it should transmit. A CI lying entirely *above* it is not a triumph, it is a leak:
recovering TP53 from H&E far better than four published pan-cancer studies is
evidence of circularity or of the covariate having entered the features, and it
fails for that reason.

Which AUROC is the comparable one
---------------------------------
Pooled pan-cancer AUROC is inflated whenever the covariate's prevalence varies by
cancer type and the representation encodes cancer type — both of which are true
here — because the probe can score by guessing lineage. Published numbers are
per-cancer, so ``within_cancer_auroc`` (computed inside each cancer with both
classes present, then weighted by cancer size) is the primary statistic and the
one graded against the pre-registered interval. The pooled figure is reported
beside it precisely so the gap between the two is visible.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site

__all__ = ["out_of_fold_scores", "within_cancer_auroc", "evaluate_known_covariate"]

_MIN_PER_CANCER = 20
_MIN_MINORITY = 5


def _auroc(labels: np.ndarray, score: np.ndarray) -> float:
    """Rank-based AUROC with explicit tie handling; NaN when a class is absent."""
    labels = np.asarray(labels).astype(bool)
    score = np.asarray(score, dtype=np.float64)
    positive, negative = int(labels.sum()), int((~labels).sum())
    if positive == 0 or negative == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    ranks = np.empty(len(score), dtype=np.float64)
    ranks[order] = np.arange(1, len(score) + 1, dtype=np.float64)
    # average ranks within ties so a constant score scores exactly 0.5
    sorted_score = score[order]
    start = 0
    for end in range(1, len(sorted_score) + 1):
        if end == len(sorted_score) or sorted_score[end] != sorted_score[start]:
            ranks[order[start:end]] = np.mean(ranks[order[start:end]])
            start = end
    return float((ranks[labels].sum() - positive * (positive + 1) / 2.0) / (positive * negative))


def out_of_fold_scores(features: np.ndarray, labels: np.ndarray, *, n_splits: int = 5,
                       alpha: float = 1.0, seed: int = 42) -> np.ndarray:
    """Out-of-fold ridge probe scores, folds stratified on the label.

    Ridge rather than logistic regression on purpose: it has one fixed
    hyperparameter, no convergence behaviour that could differ between the
    observed run and the permutation null, and AUROC only reads the ordering it
    produces. Nothing about the pass decision should depend on an optimiser.
    """
    from sklearn.linear_model import Ridge

    features = np.asarray(features, dtype=np.float64)
    labels = np.asarray(labels).astype(float)
    rng = np.random.default_rng(seed)
    fold = np.empty(len(labels), dtype=np.int64)
    for value in np.unique(labels):
        rows = rng.permutation(np.flatnonzero(labels == value))
        fold[rows] = np.arange(len(rows)) % n_splits
    score = np.empty(len(labels), dtype=np.float64)
    for f in range(n_splits):
        test = np.flatnonzero(fold == f)
        train = np.flatnonzero(fold != f)
        if len(test) == 0 or len(np.unique(labels[train])) < 2:
            score[test] = 0.0
            continue
        score[test] = Ridge(alpha=alpha).fit(features[train], labels[train]).predict(features[test])
    return score


def within_cancer_auroc(labels: np.ndarray, score: np.ndarray, cancers: np.ndarray) -> tuple[float, dict]:
    """Size-weighted mean of the per-cancer AUROCs — the published comparator."""
    labels = np.asarray(labels).astype(bool)
    cancers = np.asarray(cancers).astype(str)
    per_cancer: dict[str, dict] = {}
    values, weights = [], []
    for cancer in sorted(set(cancers)):
        rows = np.flatnonzero(cancers == cancer)
        minority = min(int(labels[rows].sum()), int((~labels[rows]).sum()))
        if len(rows) < _MIN_PER_CANCER or minority < _MIN_MINORITY:
            per_cancer[cancer] = {"status": "unavailable_insufficient_cases", "n": int(len(rows)),
                                  "n_minority": minority}
            continue
        value = _auroc(labels[rows], score[rows])
        per_cancer[cancer] = {"status": "scored", "n": int(len(rows)), "n_positive": int(labels[rows].sum()),
                              "auroc": value}
        values.append(value)
        weights.append(float(len(rows)))
    if not values:
        return float("nan"), per_cancer
    return float(np.average(values, weights=weights)), per_cancer


def evaluate_known_covariate(features: np.ndarray, labels: np.ndarray, cancers: np.ndarray,
                             patient_ids: np.ndarray, *, expected_low: float, expected_high: float,
                             n_splits: int = 5, seed: int = 42, n_boot: int = 1000,
                             n_permutations: int = 1000, min_site_count: int = 10,
                             residualise: bool = False, covariate_name: str = "") -> dict:
    """Score one state against one pre-registered covariate expectation."""
    features = np.asarray(features, dtype=np.float64)
    labels = np.asarray(labels).astype(bool)
    cancers = np.asarray(cancers).astype(str)
    if labels.sum() < 30 or (~labels).sum() < 30:
        return {"status": "unavailable_insufficient_labelled_cases", "n_positive": int(labels.sum()),
                "n_negative": int((~labels).sum())}
    adjustment = "none_raw_state"
    if residualise:
        import pandas as pd
        site, _ = pooled_tissue_source_site(np.asarray(patient_ids).astype(str), min_site_count=min_site_count)
        design = confound_design(pd.DataFrame({"cancer": cancers, "tss": site}), ["cancer", "tss"])
        features = cross_fitted_residuals(features, design, seed=seed)
        adjustment = "cancer+pooled_tss_cross_fitted"

    score = out_of_fold_scores(features, labels, n_splits=n_splits, seed=seed)
    within, per_cancer = within_cancer_auroc(labels, score, cancers)
    pooled = _auroc(labels, score)

    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(n_boot):
        rows = rng.integers(0, len(labels), size=len(labels))
        value, _ = within_cancer_auroc(labels[rows], score[rows], cancers[rows])
        if np.isfinite(value):
            boot.append(value)
    ci = ([float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
          if boot else [float("nan"), float("nan")])

    # Within-cancer label permutation: the chance level for a within-cancer AUROC,
    # measured rather than assumed to be 0.5.
    null = []
    for _ in range(n_permutations):
        permuted = labels.copy()
        for cancer in np.unique(cancers):
            rows = np.flatnonzero(cancers == cancer)
            permuted[rows] = labels[rng.permutation(rows)]
        value, _ = within_cancer_auroc(permuted, score, cancers)
        if np.isfinite(value):
            null.append(value)
    null_p95 = float(np.percentile(null, 95)) if null else float("nan")
    permutation_p = float((sum(v >= within for v in null) + 1.0) / (len(null) + 1.0)) if null else float("nan")

    overlaps = bool(np.isfinite(ci[0]) and np.isfinite(ci[1])
                    and ci[0] <= float(expected_high) and ci[1] >= float(expected_low))
    return {
        "status": "scored", "covariate": covariate_name, "adjustment": adjustment,
        "n_patients": int(len(labels)), "n_positive": int(labels.sum()),
        "prevalence": float(labels.mean()),
        "within_cancer_auroc": within, "within_cancer_auroc_ci95": ci,
        "pooled_auroc": pooled,
        "expected_interval": [float(expected_low), float(expected_high)],
        "ci_overlaps_expected": overlaps,
        "below_expected": bool(np.isfinite(ci[1]) and ci[1] < float(expected_low)),
        "above_expected": bool(np.isfinite(ci[0]) and ci[0] > float(expected_high)),
        "null_p95": null_p95, "permutation_p": permutation_p,
        "permutation_resolution": 1.0 / (len(null) + 1.0) if null else float("nan"),
        "n_cancers_scored": int(sum(1 for v in per_cancer.values() if v.get("status") == "scored")),
        "per_cancer": per_cancer,
        "passed": overlaps,
    }


def load_covariate_table(path: str, column: str) -> dict[str, bool]:
    """Read one binary covariate per canonical TCGA patient."""
    import pandas as pd
    table = (pd.read_parquet(path) if str(path).lower().endswith(".parquet")
             else pd.read_csv(path, sep=None, engine="python"))
    patient = next((c for c in ("patient_id", "bcr_patient_barcode", "Patient ID", "sample_id")
                    if c in table.columns), None)
    if patient is None or column not in table.columns:
        raise ValueError(f"covariate table needs a patient ID column and a {column!r} column")
    values = pd.to_numeric(table[column], errors="coerce")
    keep = np.isfinite(values)
    ids = table.loc[keep, patient].astype(str).map(lambda v: "-".join(str(v).split("-")[:3]))
    return dict(zip(ids, values[keep].astype(bool)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--covariate-table", required=True)
    parser.add_argument("--covariate-column", required=True)
    parser.add_argument("--expected-low", type=float, required=True)
    parser.add_argument("--expected-high", type=float, required=True)
    parser.add_argument("--preregistration", default="", help="path to the JSON written BEFORE the run")
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition", default="test", choices=("test", "val", "all"))
    parser.add_argument("--states", nargs="*", default=[])
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.preregistration and not Path(args.preregistration).is_file():
        raise SystemExit("the pre-registered expectation must exist BEFORE this run; refusing to grade without it")
    covariate = load_covariate_table(args.covariate_table, args.covariate_column)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict] = {}
    rows: list[dict] = []
    for artifact_path in args.artifacts:
        path = Path(artifact_path)
        method = path.stem
        raw = np.load(path, allow_pickle=True)
        declared = {str(s) for s in raw["trained_states"]} if "trained_states" in raw.files else set()
        if args.states:
            declared &= set(args.states)
        patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
        cancers = np.asarray([str(c) for c in raw["cancers"]])
        split = np.asarray([str(s) for s in raw["split"]])
        mask = np.ones(len(patient_ids), dtype=bool) if args.partition == "all" else (split == args.partition)
        mask &= np.asarray([pid in covariate for pid in patient_ids])
        labels = np.asarray([bool(covariate.get(pid, False)) for pid in patient_ids])[mask]
        for state in sorted(declared):
            for residualise in (False, True):
                key = f"{method}::{state}::{'adjusted' if residualise else 'raw'}"
                result = evaluate_known_covariate(
                    np.asarray(raw[state], dtype=np.float64)[mask], labels, cancers[mask],
                    patient_ids[mask], expected_low=args.expected_low, expected_high=args.expected_high,
                    seed=args.seed, n_boot=args.n_boot, n_permutations=args.n_permutations,
                    residualise=residualise, covariate_name=args.covariate_column)
                results[key] = result
                for metric in ("within_cancer_auroc", "pooled_auroc", "null_p95", "permutation_p",
                               "n_positive", "prevalence"):
                    rows.append({"method": method, "representation_state": state,
                                 "task": f"known_covariate_{'adjusted' if residualise else 'raw'}",
                                 "target": args.covariate_column, "metric": metric,
                                 "value": float(result.get(metric, float("nan"))), "note": result["status"]})
                ci = result.get("within_cancer_auroc_ci95", [float("nan")] * 2)
                print(f"[{key}] within_cancer_auroc={result.get('within_cancer_auroc', float('nan')):.4f} "
                      f"CI=[{ci[0]:.4f},{ci[1]:.4f}] expected=[{args.expected_low},{args.expected_high}] "
                      f"pooled={result.get('pooled_auroc', float('nan')):.4f} "
                      f"null_p95={result.get('null_p95', float('nan')):.4f} "
                      f"{'PASS' if result.get('passed') else 'FAIL'}", flush=True)
    (output / "known_covariate_control.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    import pandas as pd
    pd.DataFrame(rows).to_csv(output / "task_rows.csv", index=False)
    print(f"[done] -> {output}", flush=True)


if __name__ == "__main__":
    main()
