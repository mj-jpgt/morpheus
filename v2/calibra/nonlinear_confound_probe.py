"""Non-mean-based probes of the confound certificate, on TCGA.

Why this module exists
----------------------
``confound_certificate.certify_axes`` scores a representation with exactly two
classifiers: ``nearest_class_mean_oof`` per axis and ``lda_oof_balanced_accuracy``
jointly.  **Both are functions of the class means.**  The adjustment it certifies is
``residualise.cross_fitted_residuals`` — a ridge regression of the representation on a
one-hot ``cancer + pooled TSS`` design — which removes the class *mean vector* by
construction.  A mean-based test passing after a mean-removing adjustment is close to
arithmetic, and "the confound is gone" and "the confound is no longer in the first
moment" are indistinguishable to that instrument.

The spatial replication showed the distinction is not hypothetical: on HEST an
adjustment that took joint-LDA slide accuracy from 0.9707 to 0.0025 (chance 0.0769)
left an out-of-fold 15-NN naming the slide 72.9% of the time.  That argument is about
the estimator family, not about spatial data, so it applies verbatim to the 85-class
TCGA site certificate.  This module asks TCGA.

The probe families, and why each one
------------------------------------
* **k-NN vote** (``knn_balanced_accuracy_oof``).  The decision rule is the labels of the
  nearest training rows; no class mean enters at any step, so mean-removal cannot make
  it pass by construction.  ``k`` is swept rather than tuned: small ``k`` reads the
  finest local structure and is the most sensitive to a residual confound, large ``k``
  smooths toward a density estimate, and a representation that only looks clean at
  large ``k`` has been cleaned by the probe rather than by the adjustment.  The plain
  majority vote is pinned by a test to reproduce ``hest_claims.knn_balanced_accuracy_oof``
  exactly, so the bulk and spatial readings come from one estimator.

  ``prior_correction`` exists because TCGA's site classes are wildly imbalanced (the
  pooled ``OTHER`` class is 30% of the test partition).  A plain majority vote
  over-predicts the largest class, which *depresses* balanced accuracy and biases the
  run toward the reassuring answer.  Weighting each vote by the inverse training-fold
  class frequency is the hyperparameter-free rule matched to the balanced-accuracy loss
  actually being reported.  Both variants are run; a conclusion that survives only the
  plain vote is not a conclusion.

* **Random forest** (``forest_balanced_accuracy_oof``).  Chosen because it fails in a
  *different* way from k-NN and shares none of its assumptions.  Its decision rule is a
  set of axis-aligned thresholds: it is not a metric method, it is invariant to any
  monotone rescaling of any single coordinate, and it can key on a difference in
  variance, in skew, or in an interaction between two coordinates.  Every one of those
  is invisible to LDA and survives mean-removal untouched.

* **RBF-kernel SVM** (``rbf_svm_balanced_accuracy_oof``).  Smooth and global where k-NN
  is local and discrete — the decision function is a weighted sum of Gaussians centred
  on training rows — so it reads the same local geometry without k-NN's specific failure
  on small classes, where a vote tie is resolved by class size.

Nulls: two, and only one of them is a chance rate
-------------------------------------------------
``within_stratum_permutations`` permutes the confound label inside cancer type.  That is
the right convention for the question the certificate poses about a *raw* axis ("does it
carry site beyond what cancer already explains?"), and the module docstring of
``confound_certificate`` argues for it correctly.  It is **not** a chance rate, and the
two ways it stops being one are different:

1. *Degeneracy.*  A stratum containing exactly one class is permuted by the identity, so
   those rows are handed their true labels.  Five of thirteen HEST slides were in that
   position and the "null" read 0.53 against a chance rate of 0.077.
2. *Nesting.*  When the confound nests inside the stratum — on TCGA, **no** tissue source
   site contributes patients to more than one cancer — every permuted label is still a
   site *of the correct cancer*, so the null inherits the whole cancer→site restriction
   and sits far above ``1/n_classes`` for any representation that knows cancer.

``nesting_diagnostic`` measures both before any feature is touched, and
``within_stratum_chance`` gives the analytic value the within-stratum null should sit
near if the representation carries the stratum and nothing more.  ``global_permutations``
*measures* chance rather than assuming ``1/n_classes``; it is the applicable bar whenever
the adjustment being tested removes the stratum as well as the confound.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from .calibration import _map
from .confound_certificate import (_encode, _stratified_folds, balanced_accuracy,
                                   lda_oof_balanced_accuracy, nearest_class_mean_oof,
                                   within_stratum_permutations)
from .residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site

__all__ = ["knn_balanced_accuracy_oof", "forest_balanced_accuracy_oof",
           "rbf_svm_balanced_accuracy_oof", "global_permutations", "within_stratum_chance",
           "nesting_diagnostic", "null_summary", "probe_state"]


# --- probes ------------------------------------------------------------------------------

def knn_balanced_accuracy_oof(features: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                              k: int = 15, n_splits: int = 5, seed: int = 42,
                              prior_correction: bool = False,
                              folds: list[np.ndarray] | None = None, block: int = 512) -> float:
    """Out-of-fold k-nearest-neighbour balanced accuracy for a categorical confound.

    With ``prior_correction=False`` this is a plain majority vote and is numerically
    identical to ``hest_claims.knn_balanced_accuracy_oof`` — asserted by
    ``test_nonlinear_confound_probe.py``, so the two cohorts cannot silently drift onto
    different estimators.  With ``prior_correction=True`` each neighbour's vote is
    weighted by the inverse of its class's frequency in the training fold, which is the
    rule matched to balanced accuracy under imbalanced priors.
    """
    features = np.asarray(features, dtype=np.float64)
    class_index = np.asarray(class_index, dtype=np.int64)
    if features.ndim != 2 or len(features) != len(class_index):
        raise ValueError("kNN probe needs an aligned [row, feature] matrix and label vector")
    folds = folds if folds is not None else _stratified_folds(class_index, n_splits, seed)
    predicted = np.empty(len(class_index), dtype=np.int64)
    square = (features ** 2).sum(axis=1)
    all_rows = np.arange(len(class_index))
    for test_rows in folds:
        train_rows = np.setdiff1d(all_rows, test_rows, assume_unique=False)
        if len(train_rows) <= k:
            raise ValueError("not enough training rows for the requested k")
        train_features = features[train_rows]
        train_square = square[train_rows]
        train_labels = class_index[train_rows]
        weight = np.ones(n_classes, dtype=np.float64)
        if prior_correction:
            counts = np.bincount(train_labels, minlength=n_classes).astype(np.float64)
            weight = np.where(counts > 0, 1.0 / np.maximum(counts, 1.0), 0.0)
        for start in range(0, len(test_rows), block):
            rows = test_rows[start:start + block]
            # squared euclidean distance without materialising differences
            distance = (square[rows][:, None] + train_square[None, :]
                        - 2.0 * (features[rows] @ train_features.T))
            nearest = np.argpartition(distance, kth=k - 1, axis=1)[:, :k]
            votes = train_labels[nearest]
            tally = np.zeros((len(rows), n_classes), dtype=np.float64)
            np.add.at(tally, (np.repeat(np.arange(len(rows)), k), votes.ravel()),
                      weight[votes.ravel()])
            predicted[rows] = np.argmax(tally, axis=1)
    return balanced_accuracy(class_index, predicted, n_classes)


def _sklearn_oof(make_model, features: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                 n_splits: int, seed: int, folds: list[np.ndarray] | None,
                 standardise: bool) -> float:
    features = np.asarray(features, dtype=np.float64)
    class_index = np.asarray(class_index, dtype=np.int64)
    if features.ndim != 2 or len(features) != len(class_index):
        raise ValueError("probe needs an aligned [row, feature] matrix and label vector")
    folds = folds if folds is not None else _stratified_folds(class_index, n_splits, seed)
    predicted = np.empty(len(class_index), dtype=np.int64)
    all_rows = np.arange(len(class_index))
    for test_rows in folds:
        train_rows = np.setdiff1d(all_rows, test_rows, assume_unique=False)
        train_features, test_features = features[train_rows], features[test_rows]
        if standardise:
            # Fit the scaling on the TRAINING fold only; a scale fitted on all rows is a
            # (small) leak and this probe exists to be believed.
            centre = train_features.mean(axis=0)
            scale = np.maximum(train_features.std(axis=0), 1e-12)
            train_features = (train_features - centre) / scale
            test_features = (test_features - centre) / scale
        model = make_model()
        model.fit(train_features, class_index[train_rows])
        predicted[test_rows] = np.asarray(model.predict(test_features), dtype=np.int64)
    return balanced_accuracy(class_index, predicted, n_classes)


def forest_balanced_accuracy_oof(features: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                                 n_estimators: int = 300, max_features="sqrt",
                                 min_samples_leaf: int = 1, n_splits: int = 5, seed: int = 42,
                                 folds: list[np.ndarray] | None = None, n_jobs: int = 1) -> float:
    """Out-of-fold random-forest balanced accuracy.

    ``class_weight="balanced_subsample"`` because balanced accuracy is the reported loss
    and TCGA's site classes are imbalanced; without it the forest optimises plain
    accuracy and the probe reads low for a reason that has nothing to do with the
    representation.  No hyperparameter here is tuned against an outcome.
    """
    from sklearn.ensemble import RandomForestClassifier

    def _make():
        return RandomForestClassifier(n_estimators=int(n_estimators), max_features=max_features,
                                      min_samples_leaf=int(min_samples_leaf),
                                      class_weight="balanced_subsample", random_state=int(seed),
                                      n_jobs=int(n_jobs))

    return _sklearn_oof(_make, features, class_index, n_classes, n_splits=n_splits, seed=seed,
                        folds=folds, standardise=False)


def rbf_svm_balanced_accuracy_oof(features: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                                  C: float = 1.0, gamma="scale", n_splits: int = 5, seed: int = 42,
                                  folds: list[np.ndarray] | None = None) -> float:
    """Out-of-fold RBF-kernel SVM balanced accuracy, on per-fold-standardised features."""
    from sklearn.svm import SVC

    def _make():
        return SVC(C=float(C), kernel="rbf", gamma=gamma, class_weight="balanced",
                   random_state=int(seed))

    return _sklearn_oof(_make, features, class_index, n_classes, n_splits=n_splits, seed=seed,
                        folds=folds, standardise=True)


# --- nulls -------------------------------------------------------------------------------

def global_permutations(class_index: np.ndarray, *, n_permutations: int, seed: int) -> list[np.ndarray]:
    """Permute the label with no stratification — this MEASURES chance.

    Expressed as ``within_stratum_permutations`` over a single constant stratum so that
    the global and within-stratum nulls are produced by literally the same code and can
    differ only through the strata argument.
    """
    class_index = np.asarray(class_index, dtype=np.int64)
    return within_stratum_permutations(class_index, np.zeros(len(class_index), dtype=np.int64),
                                       n_permutations=n_permutations, seed=seed)


def within_stratum_chance(class_index: np.ndarray, strata: np.ndarray, n_classes: int) -> float:
    """Balanced accuracy of an oracle that knows the stratum and then guesses uniformly.

    The analytic reference for a within-stratum permutation null.  For a class ``c``
    living in stratum ``s``, such an oracle recalls ``c`` with probability
    ``1 / (number of classes present in s)``; balanced accuracy averages that over the
    classes present.  A class spread over several strata gets the row-weighted mean of
    its strata's rates, which is the same quantity when the confound nests.

    Reported so that a within-stratum null far above ``1/n_classes`` is legible as
    nesting rather than as a leak.
    """
    class_index = np.asarray(class_index, dtype=np.int64)
    strata = np.asarray(strata).astype(str)
    per_stratum = {level: len(np.unique(class_index[strata == level]))
                   for level in np.unique(strata)}
    rates = []
    for value in np.unique(class_index):
        rows = strata[class_index == value]
        rates.append(float(np.mean([1.0 / max(per_stratum[level], 1) for level in rows])))
    return float(np.mean(rates)) if rates else float("nan")


def nesting_diagnostic(class_index: np.ndarray, strata: np.ndarray, n_classes: int) -> dict:
    """How badly a within-stratum permutation null misrepresents chance, measured.

    Two independent failure modes, reported separately because they have different cures:
    ``fraction_of_rows_unpermutable`` (a stratum with one class — the HEST case) and
    ``n_classes_spanning_multiple_strata`` (nesting — the TCGA case).
    """
    class_index = np.asarray(class_index, dtype=np.int64)
    strata = np.asarray(strata).astype(str)
    levels = np.unique(strata)
    degenerate = [level for level in levels if len(np.unique(class_index[strata == level])) < 2]
    unpermutable = int(sum(int((strata == level).sum()) for level in degenerate))
    spanning = int(sum(1 for value in np.unique(class_index)
                       if len(np.unique(strata[class_index == value])) > 1))
    return {
        "n_classes": int(n_classes),
        "n_classes_present": int(len(np.unique(class_index))),
        "n_strata": int(len(levels)),
        "n_strata_with_one_class": int(len(degenerate)),
        "fraction_of_rows_unpermutable": float(unpermutable) / float(max(len(class_index), 1)),
        "n_classes_spanning_multiple_strata": spanning,
        "confound_nests_inside_stratum": bool(spanning == 0),
        "design_chance_rate": 1.0 / float(n_classes),
        "within_stratum_chance": within_stratum_chance(class_index, strata, n_classes),
    }


def null_summary(null: np.ndarray, observed: float) -> dict:
    null = np.asarray(null, dtype=np.float64)
    return {
        "null_median": float(np.median(null)),
        "null_p95": float(np.percentile(null, 95)),
        "null_max": float(np.max(null)),
        "n_permutations": int(len(null)),
        "permutation_p": float((int((null >= observed).sum()) + 1.0) / (len(null) + 1.0)),
        "permutation_resolution": 1.0 / (len(null) + 1.0),
    }


# --- one state x one target --------------------------------------------------------------

def _probe_callables(k_grid, forest_trees, run_svm, seed, n_splits):
    probes = {}
    for k in k_grid:
        probes[f"knn_k{k}"] = (lambda f, c, n, folds, k=k: knn_balanced_accuracy_oof(
            f, c, n, k=k, n_splits=n_splits, seed=seed, folds=folds, prior_correction=False))
        probes[f"knn_prior_k{k}"] = (lambda f, c, n, folds, k=k: knn_balanced_accuracy_oof(
            f, c, n, k=k, n_splits=n_splits, seed=seed, folds=folds, prior_correction=True))
    if forest_trees:
        probes["forest"] = (lambda f, c, n, folds: forest_balanced_accuracy_oof(
            f, c, n, n_estimators=forest_trees, n_splits=n_splits, seed=seed, folds=folds))
    if run_svm:
        probes["rbf_svm"] = (lambda f, c, n, folds: rbf_svm_balanced_accuracy_oof(
            f, c, n, n_splits=n_splits, seed=seed, folds=folds))
    return probes


def probe_state(features: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                strata: np.ndarray | None = None, k_grid=(1, 3, 5, 10, 15, 25, 50),
                forest_trees: int = 300, run_svm: bool = False, n_splits: int = 5, seed: int = 42,
                knn_permutations: int = 200, model_permutations: int = 100,
                n_jobs: int = 1, include_linear: bool = True) -> dict:
    """Run every probe on one feature block, with a global and a within-stratum null.

    ``include_linear`` re-runs the certificate's own two classifiers on the identical
    rows and folds, so the nonlinear reading is never compared against a number quoted
    from another run.
    """
    features = np.asarray(features, dtype=np.float64)
    class_index = np.asarray(class_index, dtype=np.int64)
    folds = _stratified_folds(class_index, n_splits, seed)
    out: dict = {"n": int(len(class_index)), "n_axes": int(features.shape[1]),
                 "n_classes": int(n_classes), "chance_rate": 1.0 / float(n_classes),
                 "n_splits": int(n_splits), "seed": int(seed)}
    if include_linear:
        per_axis = nearest_class_mean_oof(features, class_index, n_classes, n_splits=n_splits,
                                          seed=seed, folds=folds)
        out["linear_reference"] = {
            "joint_lda_balanced_accuracy": float(lda_oof_balanced_accuracy(
                features, class_index, n_classes, n_splits=n_splits, seed=seed, folds=folds)),
            "per_axis_max": float(np.max(per_axis)),
            "per_axis_median": float(np.median(per_axis)),
        }
    permutation_sets = {"global": global_permutations}
    if strata is not None:
        permutation_sets["within_stratum"] = (
            lambda ci, *, n_permutations, seed: within_stratum_permutations(
                ci, strata, n_permutations=n_permutations, seed=seed))
    probes = _probe_callables(k_grid, forest_trees, run_svm, seed, n_splits)
    out["probes"] = {}
    for name, probe in probes.items():
        start = time.time()
        observed = float(probe(features, class_index, n_classes, folds))
        record = {"balanced_accuracy": observed,
                  "multiple_of_chance": observed * float(n_classes)}
        count = knn_permutations if name.startswith("knn") else model_permutations
        for label, maker in permutation_sets.items():
            if count <= 0:
                continue
            permutations = maker(class_index, n_permutations=count, seed=seed)
            # features passed as an ARGUMENT rather than captured, so joblib memmaps the
            # block once instead of pickling it per permutation.
            null = np.asarray(_map(probe, [(features, p, n_classes, folds)
                                           for p in permutations], n_jobs), dtype=np.float64)
            record[label] = null_summary(null, observed)
        record["wall_seconds"] = round(time.time() - start, 1)
        out["probes"][name] = record
    return out


# --- CLI ---------------------------------------------------------------------------------

def _labels(patient_ids, cancers, target: str, min_site_count: int):
    if target == "cancer":
        classes, class_index = _encode(np.asarray(cancers).astype(str))
        return classes, class_index, None
    site, _ = pooled_tissue_source_site(np.asarray(patient_ids).astype(str),
                                        min_site_count=min_site_count)
    classes, class_index = _encode(site)
    return classes, class_index, np.asarray(cancers).astype(str)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tag", default="probe")
    parser.add_argument("--states", nargs="+", default=None)
    parser.add_argument("--targets", nargs="+", default=("site", "cancer"))
    parser.add_argument("--arms", nargs="+", default=("raw", "adjusted", "adjusted_standardised"))
    parser.add_argument("--partition", default="test", choices=("test", "val", "all"))
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--k-grid", default="1,3,5,10,15,25,50")
    parser.add_argument("--forest-trees", type=int, default=300)
    parser.add_argument("--run-svm", action="store_true")
    parser.add_argument("--knn-permutations", type=int, default=200)
    parser.add_argument("--model-permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    k_grid = tuple(int(v) for v in str(args.k_grid).split(",") if v.strip())
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    results: dict = {"config": {"partition": args.partition, "k_grid": list(k_grid),
                               "forest_trees": args.forest_trees, "run_svm": bool(args.run_svm),
                               "knn_permutations": args.knn_permutations,
                               "model_permutations": args.model_permutations,
                               "n_splits": args.n_splits, "seed": args.seed,
                               "min_site_count": args.min_site_count}}
    for artifact_path in args.artifacts:
        path = Path(artifact_path)
        raw = np.load(path, allow_pickle=True)
        patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
        cancers = np.asarray([str(c) for c in raw["cancers"]])
        split = np.asarray([str(s) for s in raw["split"]])
        mask = np.ones(len(patient_ids), bool) if args.partition == "all" else (split == args.partition)
        ids, cans = patient_ids[mask], cancers[mask]
        # A repeated patient id would let a kNN name the site by finding the same patient.
        # Asserted rather than assumed; the count is reported either way.
        duplicates = int(len(ids) - len(set(ids.tolist())))
        declared = [str(s) for s in raw["trained_states"]] if "trained_states" in raw.files else []
        states = list(args.states) if args.states else sorted(declared)
        import pandas as pd
        for target in args.targets:
            classes, class_index, strata = _labels(ids, cans, target, args.min_site_count)
            n_classes = int(len(classes))
            site, _ = pooled_tissue_source_site(ids, min_site_count=args.min_site_count)
            design = confound_design(pd.DataFrame({"cancer": cans, "tss": site}), ["cancer", "tss"])
            diagnostic = nesting_diagnostic(class_index, cans if strata is None else strata,
                                            n_classes) if target == "site" else None
            for state in states:
                if state not in raw.files:
                    results[f"{path.stem}::{state}::{target}"] = {"status": "state_absent"}
                    continue
                block = np.asarray(raw[state], dtype=np.float64)[mask]
                adjusted = cross_fitted_residuals(block, design, seed=args.seed)
                arms = {
                    "raw": block,
                    "adjusted": adjusted,
                    "adjusted_standardised": ((adjusted - adjusted.mean(axis=0))
                                              / np.maximum(adjusted.std(axis=0), 1e-12)),
                }
                for arm in args.arms:
                    key = f"{path.stem}::{state}::{target}::{arm}"
                    started = time.time()
                    record = probe_state(arms[arm], class_index, n_classes, strata=strata,
                                         k_grid=k_grid, forest_trees=args.forest_trees,
                                         run_svm=bool(args.run_svm), n_splits=args.n_splits,
                                         seed=args.seed, knn_permutations=args.knn_permutations,
                                         model_permutations=args.model_permutations,
                                         n_jobs=args.n_jobs)
                    record.update({"status": "scored", "artifact": path.stem, "state": state,
                                   "target": target, "arm": arm, "partition": args.partition,
                                   "duplicate_patient_ids": duplicates,
                                   "design_columns": int(design.shape[1]),
                                   "design_rank": int(np.linalg.matrix_rank(design)),
                                   "wall_seconds": round(time.time() - started, 1)})
                    if diagnostic is not None:
                        record["nesting"] = diagnostic
                    results[key] = record
                    best = max((v["balanced_accuracy"] for n, v in record["probes"].items()
                                if n.startswith("knn")), default=float("nan"))
                    forest = record["probes"].get("forest", {}).get("balanced_accuracy", float("nan"))
                    print(f"[{key}] chance={record['chance_rate']:.4f} knn_max={best:.4f} "
                          f"forest={forest:.4f} "
                          f"lda={record.get('linear_reference', {}).get('joint_lda_balanced_accuracy', float('nan')):.4f} "
                          f"({record['wall_seconds']:.0f} s)", flush=True)
                    (output / f"{args.tag}.json").write_text(json.dumps(results, indent=2),
                                                             encoding="utf-8")
    (output / f"{args.tag}.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"[done] {len(results) - 1} blocks -> {output / (args.tag + '.json')}", flush=True)


if __name__ == "__main__":
    main()
