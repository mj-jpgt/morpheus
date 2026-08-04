"""Positive control for a newly built modality block: it must recover lineage.

Every negative result on this project is reportable only if a positive control
passed in the same run, on the same data, through the same code path
(``HANDOFF_GATES.md``).  When a *new modality* is wired in, the specific way it
can go wrong is not subtle: a mangled identifier join, an inverted mask, a
transposed matrix or a silently-empty column produces a block that is shaped
correctly, hashes fine, and carries nothing.  Every downstream "no channel"
result would then be a builder bug wearing the costume of a finding.

The control this module implements is the cheapest one that can distinguish
those two: **a molecular block must recover cancer type.**  Somatic mutation,
copy number and protein abundance are all strongly lineage-patterned --- BRAF in
SKCM, APC and KRAS in colorectum, VHL in KIRC, arm-level gains that are near
diagnostic of a tissue.  A block from which lineage is unreadable is broken.

Note what this control does *not* license.  Passing it says the block carries
information, not that the modality is legible from morphology; those are
different questions and the second one is measured elsewhere, by
``run_calibra``, against a floor.  A block can pass here and still show no
cross-modal channel at all, and that combination is the interesting one --- it
is exactly the case in which a null is a statement about the representation
rather than about the pipeline.

Design choices, each for the same reason ``known_covariate_control`` gives:
nothing about the pass decision should depend on an optimiser.  The probe is
ridge regression onto one-hot lineage with one fixed alpha and an argmax
readout, so the observed run and the permutation null go through identical
arithmetic.  Chance is **measured** by permuting the lineage labels, never
assumed to be ``1/n_cancers``: class sizes here are wildly unequal and an
accuracy-shaped statistic under imbalance does not sit where intuition puts it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

__all__ = ["balanced_accuracy", "out_of_fold_lineage_prediction", "evaluate_modality_block"]


def balanced_accuracy(labels: np.ndarray, predicted: np.ndarray) -> tuple[float, dict[str, float]]:
    """Mean per-class recall, plus the per-class recalls it averaged.

    Unweighted on purpose: a weighted score on this cohort is dominated by the
    three largest cancers, and a block that recovered only those would read as a
    pass.
    """
    labels = np.asarray(labels).astype(str)
    predicted = np.asarray(predicted).astype(str)
    per_class: dict[str, float] = {}
    for value in sorted(set(labels.tolist())):
        rows = labels == value
        per_class[value] = float((predicted[rows] == value).mean())
    return float(np.mean(list(per_class.values()))) if per_class else float("nan"), per_class


def out_of_fold_lineage_prediction(features: np.ndarray, cancers: np.ndarray, *, n_splits: int = 5,
                                   alpha: float = 1.0, seed: int = 42) -> np.ndarray:
    """Out-of-fold argmax-ridge lineage prediction, folds stratified on lineage."""
    from sklearn.linear_model import Ridge

    features = np.asarray(features, dtype=np.float64)
    cancers = np.asarray(cancers).astype(str)
    classes = np.asarray(sorted(set(cancers.tolist())))
    index = {value: position for position, value in enumerate(classes)}
    onehot = np.zeros((len(cancers), len(classes)), dtype=np.float64)
    onehot[np.arange(len(cancers)), [index[value] for value in cancers]] = 1.0

    rng = np.random.default_rng(seed)
    fold = np.empty(len(cancers), dtype=np.int64)
    for value in classes:
        rows = rng.permutation(np.flatnonzero(cancers == value))
        fold[rows] = np.arange(len(rows)) % n_splits

    predicted = np.empty(len(cancers), dtype=object)
    for f in range(n_splits):
        test = np.flatnonzero(fold == f)
        train = np.flatnonzero(fold != f)
        if len(test) == 0 or len(train) == 0:
            continue
        model = Ridge(alpha=alpha).fit(features[train], onehot[train])
        predicted[test] = classes[np.argmax(model.predict(features[test]), axis=1)]
    return np.asarray([("" if value is None else str(value)) for value in predicted])


def evaluate_modality_block(features: np.ndarray, cancers: np.ndarray, *, n_splits: int = 5,
                            alpha: float = 1.0, seed: int = 42, n_permutations: int = 100,
                            minimum_balanced_accuracy: float = 0.25, block_name: str = "") -> dict:
    """Score one modality block against the pre-registered lineage-recovery bar."""
    features = np.asarray(features, dtype=np.float64)
    cancers = np.asarray(cancers).astype(str)
    if features.ndim != 2 or len(features) != len(cancers):
        raise ValueError("modality block control needs aligned [patient, feature] values and lineages")
    if not np.isfinite(features).all():
        raise ValueError("modality block contains non-finite values")

    predicted = out_of_fold_lineage_prediction(features, cancers, n_splits=n_splits, alpha=alpha, seed=seed)
    observed, per_class = balanced_accuracy(cancers, predicted)

    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_permutations):
        shuffled = cancers[rng.permutation(len(cancers))]
        null_predicted = out_of_fold_lineage_prediction(features, shuffled, n_splits=n_splits,
                                                        alpha=alpha, seed=seed)
        value, _ = balanced_accuracy(shuffled, null_predicted)
        if np.isfinite(value):
            null.append(value)
    null_p95 = float(np.percentile(null, 95)) if null else float("nan")
    null_median = float(np.median(null)) if null else float("nan")
    permutation_p = (float((sum(v >= observed for v in null) + 1.0) / (len(null) + 1.0))
                     if null else float("nan"))

    above_null = bool(np.isfinite(null_p95) and observed > null_p95)
    above_absolute = bool(observed >= minimum_balanced_accuracy)
    return {
        "status": "scored", "block": block_name,
        "n_patients": int(len(cancers)), "n_features": int(features.shape[1]),
        "n_cancers": int(len(set(cancers.tolist()))),
        "balanced_accuracy": observed, "per_cancer_recall": per_class,
        "null_balanced_accuracy_median": null_median, "null_balanced_accuracy_p95": null_p95,
        "n_permutations": len(null),
        "permutation_p": permutation_p,
        "permutation_resolution": 1.0 / (len(null) + 1.0) if null else float("nan"),
        "minimum_balanced_accuracy": float(minimum_balanced_accuracy),
        "above_measured_null": above_null, "above_absolute_bar": above_absolute,
        "passed": bool(above_null and above_absolute),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", required=True, help="target artifact in the frozen-target schema")
    parser.add_argument("--exclude-groups", nargs="*", default=["random_control"],
                        help="target groups to leave out of the probe")
    parser.add_argument("--control-groups", nargs="*", default=["random_control"],
                        help="target groups to score SEPARATELY as a must-fail comparison")
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition", default="all", choices=("train", "val", "test", "all"))
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-permutations", type=int, default=100)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--minimum-balanced-accuracy", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    raw = np.load(args.targets, allow_pickle=True)
    groups = np.asarray([str(g) for g in raw["target_groups"]])
    cancers = np.asarray([str(c) for c in raw["cancers"]])
    split = np.asarray([str(s) for s in raw["split"]])
    scores = np.asarray(raw["scores"], dtype=np.float64)
    rows = np.ones(len(cancers), dtype=bool) if args.partition == "all" else (split == args.partition)

    results = {}
    real = ~np.isin(groups, args.exclude_groups)
    results["real"] = evaluate_modality_block(
        scores[np.ix_(rows, real)], cancers[rows], n_splits=args.n_splits, alpha=args.alpha,
        seed=args.seed, n_permutations=args.n_permutations,
        minimum_balanced_accuracy=args.minimum_balanced_accuracy,
        block_name=f"{Path(args.targets).stem}::real")
    control = np.isin(groups, args.control_groups)
    if control.any():
        results["control"] = evaluate_modality_block(
            scores[np.ix_(rows, control)], cancers[rows], n_splits=args.n_splits, alpha=args.alpha,
            seed=args.seed, n_permutations=args.n_permutations,
            minimum_balanced_accuracy=args.minimum_balanced_accuracy,
            block_name=f"{Path(args.targets).stem}::random_control")

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "modality_block_control.json").write_text(
        json.dumps({"partition": args.partition, "targets": str(Path(args.targets).resolve()),
                    "results": results}, indent=2), encoding="utf-8")
    for key, result in results.items():
        print(f"[{key}] balanced_accuracy={result['balanced_accuracy']:.4f} "
              f"null_p95={result['null_balanced_accuracy_p95']:.4f} "
              f"null_median={result['null_balanced_accuracy_median']:.4f} "
              f"p={result['permutation_p']:.4g} n={result['n_patients']} d={result['n_features']} "
              f"{'PASS' if result['passed'] else 'FAIL'}", flush=True)
    print(f"[done] -> {output}", flush=True)


if __name__ == "__main__":
    main()
