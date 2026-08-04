"""Nonlinear confound adjustment, and the channel measured through it.

Why this module exists
----------------------
``NOTEBOOK_ENTRIES/tcga_nonlinear_confound_probe_result_20260804T2100Z.md`` established
that after ``residualise.cross_fitted_residuals`` -- the ridge-on-one-hot adjustment
CALIBRA applies before every channel number in P1 -- an out-of-fold k-NN still recovers
tissue source site at 4.3-4.9x chance and cancer type at 3.4-4.9x chance, at the
permutation floor.  The residual is real.  The open question is whether it is large
enough to *produce* the channel, and answering it needs an adjustment the residual
cannot survive, plus a check that the adjustment actually worked.

A structural fact that shapes every arm here
--------------------------------------------
``confound_design(frame, ["cancer", "tss"])`` is **purely one-hot and purely additive**.
For a categorical predictor set the linear model on a one-hot design already spans the
conditional-mean model, so the only thing a nonlinear learner can add on top is the
``cancer x tss`` **interaction** -- and ``nonlinear_confound_probe.nesting_diagnostic``
measured that TCGA site **nests inside** cancer (0 of 84 kept sites contributes patients
to two cancers).  A site indicator therefore already determines the cell for all but the
pooled ``OTHER`` rows, and the interaction is nearly degenerate by construction.

So "a nonlinear estimate of ``E[X | cancer, site]``" is predicted in advance to be close
to a no-op on this design.  :func:`saturated_cell_residuals` measures the *limit* of that
family rather than arguing about it: it removes the exact conditional mean of every
observed ``(cancer, site)`` cell, which no kernel, forest or boosting estimate of the same
conditional mean can beat.  If the k-NN recovery survives **that**, the surviving
structure is provably not a first-moment effect, and the only adjustment that can touch
it is one that removes a higher moment -- :func:`cross_fitted_location_scale`, which is
offered as a *different and stronger* operation, never as "the same adjustment made
nonlinear".

The adjuster contract
---------------------
Every adjuster in this module is a callable ``matrix -> residual`` with the *identical*
fold contract as ``cross_fitted_residuals``: ``KFold(n_splits=min(n_splits, len(matrix)),
shuffle=True, random_state=seed)``, nuisance fit on the training rows only, and a final
column re-centring.  That is what makes the arms comparable to each other and to every
published number.  :func:`make_adjuster` at ``model="ridge"`` returns
``cross_fitted_residuals`` itself, so the incumbent arm is not a re-implementation --
``v2/tests/test_nonlinear_adjustment.py`` pins both that identity and the identity of
:func:`channel_under_adjustment` with ``calibration.permutation_null`` under it.

Why kernel ridge is the primary nonlinear arm
---------------------------------------------
Three reasons, all stated in the predeclaration before any number existed:

1. On a one-hot design the RBF kernel is a monotone function of the Hamming distance
   between two ``(cancer, site)`` rows, so its hypothesis space is exactly "smoothed cell
   means with arbitrary interaction" -- a strict superset of the additive ridge's span.
2. It is multi-output in one solve: all 256 image axes and all molecular targets come
   from a single factorisation.
3. The kernel and its per-fold Cholesky factors depend on the **design only, never on
   Y**, so :class:`KernelRidgeAdjuster` factorises once and every permutation of the
   pairing null is a triangular solve.  A reduced permutation count was the one declared
   item that defeated the previous run (its forest null was cut from 100 to 50); this
   choice removes that failure mode instead of repeating it.

A multi-output random forest is the second family, because it fails differently: axis-
aligned thresholds on indicator columns, arbitrary depth-driven interaction, no metric
assumption.  It cannot be precomputed, so its null is affordable only at a reduced count,
which is declared rather than discovered.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from .calibration import _map, permutation_null
from .residualise import (apply_pooled_tissue_source_site, confound_design,
                          cross_fitted_residuals, pooled_tissue_source_site)
from .spectral import (effective_rank, heldout_top_cca, top_canonical_correlation)

__all__ = ["cell_codes", "cell_design", "saturated_cell_residuals", "kernel_ridge_residuals",
           "KernelRidgeAdjuster", "forest_residuals", "cross_fitted_location_scale",
           "make_adjuster", "adjuster_agreement", "cross_fitted_r2",
           "channel_under_adjustment", "labels_only_ceiling", "ADJUSTER_MODELS",
           "retention_of_excess", "corrected_multiple",
           "in_sample_residuals", "row_shuffled", "cross_fitting_offset_energy",
           "regenerated_adjustment_null"]

#: Every adjustment arm this module can build, by name.  ``ridge`` is the incumbent and
#: is ``cross_fitted_residuals`` itself, not a copy.
ADJUSTER_MODELS = ("none", "ridge", "saturated", "kernel_ridge", "forest", "location_scale",
                   "in_sample_additive", "in_sample_saturated")


def _folds(n_rows: int, n_splits: int, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """The residualiser's own fold contract, in one place.

    ``cross_fitted_residuals`` builds ``KFold(n_splits=min(n_splits, len(matrix)),
    shuffle=True, random_state=seed)``.  Every adjuster here uses this, so a difference
    between two arms is never a difference in how the rows were split.
    """
    from sklearn.model_selection import KFold

    splitter = KFold(n_splits=min(n_splits, n_rows), shuffle=True, random_state=seed)
    return [(train, test) for train, test in splitter.split(np.zeros((n_rows, 1)))]


# --- cells -------------------------------------------------------------------------------

def cell_codes(*label_arrays) -> tuple[np.ndarray, np.ndarray]:
    """Integer code per row for the joint level of the supplied categorical arrays.

    Returns ``(codes, levels)`` where ``levels`` are the joined string labels in the
    order of the codes.  Used to build the saturated ``cancer x site`` design and the
    cells of the location-scale adjustment from one definition, so the two cannot drift.
    """
    columns = [np.asarray(a).astype(str) for a in label_arrays]
    if not columns or any(len(c) != len(columns[0]) for c in columns):
        raise ValueError("cell_codes needs at least one array and all must be the same length")
    joined = np.asarray(["\x1f".join(parts) for parts in zip(*columns)])
    levels, codes = np.unique(joined, return_inverse=True)
    return codes.astype(np.int64), levels


def cell_design(codes: np.ndarray, n_cells: int | None = None) -> np.ndarray:
    """One-hot indicator matrix of the saturated cell structure.

    The column space of this design is every function of ``(cancer, site)`` there is, so
    a linear model on it is the *most flexible conditional-mean model that exists* for a
    categorical confound.  That is exactly why it bounds the whole
    kernel/forest/boosting family from above on this question.
    """
    codes = np.asarray(codes, dtype=np.int64)
    width = int(n_cells if n_cells is not None else (codes.max() + 1 if codes.size else 0))
    design = np.zeros((len(codes), width), dtype=np.float64)
    if width:
        design[np.arange(len(codes)), codes] = 1.0
    return design


def saturated_cell_residuals(matrix: np.ndarray, codes: np.ndarray, *, n_splits: int = 5,
                             alpha: float = 1e-6, seed: int = 42) -> np.ndarray:
    """Cross-fitted removal of the exact conditional mean of every observed cell.

    ``alpha`` is nominal rather than a shrinkage choice: the cell design is orthogonal
    by construction (one indicator per row), so the ridge solution is the cell mean
    scaled by ``n_cell / (n_cell + alpha)`` and a tiny alpha only regularises cells a
    training fold never saw.  Passing the cell design to ``cross_fitted_residuals``
    rather than reimplementing the fit is deliberate -- the fold contract, the intercept
    and the final re-centring then come from the incumbent code path.
    """
    return cross_fitted_residuals(matrix, cell_design(codes), n_splits=n_splits,
                                  alpha=float(alpha), seed=seed)


# --- kernel ridge ------------------------------------------------------------------------

def _rbf_kernel(design: np.ndarray, gamma: float) -> np.ndarray:
    """``exp(-gamma * ||a - b||^2)`` for every row pair of ``design``.

    On a one-hot confound design ``||a - b||^2`` is the Hamming distance between two
    ``(cancer, site)`` rows and takes very few values (0, 2 or 4 here), so ``gamma``
    sweeps kernel width from near-identity to near-constant across the grid the
    predeclaration fixed.
    """
    design = np.asarray(design, dtype=np.float64)
    square = (design ** 2).sum(axis=1)
    distance = square[:, None] + square[None, :] - 2.0 * (design @ design.T)
    return np.exp(-float(gamma) * np.maximum(distance, 0.0))


class KernelRidgeAdjuster:
    """Cross-fitted RBF kernel-ridge residualiser that factorises once.

    The whole point: the kernel and its per-fold Cholesky factors are functions of the
    **design**, and the design does not change when the pairing null permutes ``y``.
    Building them in ``__init__`` turns every permutation into a triangular solve, which
    is what lets the null run at the project's own 2,000 permutations instead of being
    quietly cut.

    No intercept term is fitted; the training-fold column mean is removed before the
    solve and added back to the prediction, which is the kernel-ridge equivalent of
    ``Ridge(fit_intercept=True)`` and keeps the arm comparable to the incumbent.
    """

    def __init__(self, design: np.ndarray, *, alpha: float = 1.0, gamma: float = 0.5,
                 n_splits: int = 5, seed: int = 42):
        from scipy.linalg import cho_factor

        self.design = np.asarray(design, dtype=np.float64)
        self.alpha, self.gamma = float(alpha), float(gamma)
        self.n_splits, self.seed = int(n_splits), int(seed)
        self.n_rows = len(self.design)
        self.kernel = _rbf_kernel(self.design, self.gamma)
        self.folds = _folds(self.n_rows, self.n_splits, self.seed)
        self._factors = []
        for train, _ in self.folds:
            block = self.kernel[np.ix_(train, train)]
            block = block + self.alpha * np.eye(len(train))
            self._factors.append(cho_factor(block, lower=True, check_finite=False))

    def __call__(self, matrix: np.ndarray) -> np.ndarray:
        from scipy.linalg import cho_solve

        matrix = np.asarray(matrix, dtype=np.float64)
        if len(matrix) != self.n_rows:
            raise ValueError("kernel adjuster was factorised for a different number of rows")
        residual = np.empty_like(matrix)
        for (train, test), factor in zip(self.folds, self._factors):
            centre = matrix[train].mean(axis=0, keepdims=True)
            dual = cho_solve(factor, matrix[train] - centre, check_finite=False)
            prediction = self.kernel[np.ix_(test, train)] @ dual + centre
            residual[test] = matrix[test] - prediction
        return residual - residual.mean(axis=0, keepdims=True)


def kernel_ridge_residuals(matrix: np.ndarray, design: np.ndarray, *, alpha: float = 1.0,
                           gamma: float = 0.5, n_splits: int = 5, seed: int = 42) -> np.ndarray:
    """One-shot form of :class:`KernelRidgeAdjuster`, for callers scoring a single matrix."""
    return KernelRidgeAdjuster(design, alpha=alpha, gamma=gamma, n_splits=n_splits,
                               seed=seed)(matrix)


# --- forest ------------------------------------------------------------------------------

def forest_residuals(matrix: np.ndarray, design: np.ndarray, *, n_estimators: int = 300,
                     min_samples_leaf: int = 1, n_splits: int = 5, seed: int = 42,
                     n_jobs: int = 1) -> np.ndarray:
    """Cross-fitted multi-output random-forest residuals against the same design.

    Multi-output so that one fit per fold serves every column; a per-column forest would
    make the permutation null unaffordable and would not change the estimand.  Trees on
    a one-hot design split on "is this level or not", which is the clean encoding for a
    categorical predictor, and depth supplies the interaction the additive design cannot
    express.
    """
    from sklearn.ensemble import RandomForestRegressor

    matrix = np.asarray(matrix, dtype=np.float64)
    design = np.asarray(design, dtype=np.float64)
    residual = np.empty_like(matrix)
    for train, test in _folds(len(matrix), n_splits, seed):
        model = RandomForestRegressor(n_estimators=int(n_estimators),
                                      min_samples_leaf=int(min_samples_leaf),
                                      random_state=int(seed), n_jobs=int(n_jobs))
        model.fit(design[train], matrix[train])
        residual[test] = matrix[test] - np.asarray(model.predict(design[test])).reshape(
            matrix[test].shape)
    return residual - residual.mean(axis=0, keepdims=True)


# --- second moment -----------------------------------------------------------------------

def cross_fitted_location_scale(matrix: np.ndarray, codes: np.ndarray, *, n_splits: int = 5,
                                seed: int = 42, prior_count: float = 10.0) -> np.ndarray:
    """Cross-fitted per-cell mean removal **and** per-cell per-axis scale normalisation.

    **This is not a residualisation and is not offered as "the same adjustment made
    nonlinear".**  A residual is a first-moment operation whatever the flexibility of the
    nuisance model: ``X - E[X | C]`` still carries ``Var(X | C)``, the conditional skew
    and every conditional interaction between axes.  Those are precisely what a k-NN, a
    forest and an RBF-SVM read and what a class-mean certificate cannot see.  This
    function removes the second of them.

    The per-cell variance is shrunk toward the pooled variance at ``prior_count``
    pseudo-observations, ``(n_c * v_c + m * v_pool) / (n_c + m)``, because TCGA cells hold
    ~33 patients against 256 axes and an unshrunk per-cell SD would be dominated by its
    own sampling error -- which would inject noise rather than remove structure, and
    would make a channel collapse uninterpretable.  Cells absent from a training fold
    fall back to the pooled mean and variance, which is the same rule evaluated at
    ``n_c = 0``.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    codes = np.asarray(codes, dtype=np.int64)
    if len(codes) != len(matrix):
        raise ValueError("location-scale adjustment needs one cell code per row")
    n_cells = int(codes.max()) + 1 if codes.size else 0
    prior = float(prior_count)
    adjusted = np.empty_like(matrix)
    for train, test in _folds(len(matrix), n_splits, seed):
        train_matrix, train_codes = matrix[train], codes[train]
        pooled_mean = train_matrix.mean(axis=0)
        pooled_var = train_matrix.var(axis=0)
        counts = np.bincount(train_codes, minlength=n_cells).astype(np.float64)
        totals = np.zeros((n_cells, matrix.shape[1]))
        np.add.at(totals, train_codes, train_matrix)
        safe = np.maximum(counts, 1.0)[:, None]
        means = np.where(counts[:, None] > 0, totals / safe, pooled_mean[None, :])
        squares = np.zeros((n_cells, matrix.shape[1]))
        np.add.at(squares, train_codes, (train_matrix - means[train_codes]) ** 2)
        variances = np.where(counts[:, None] > 0, squares / safe, pooled_var[None, :])
        shrunk = ((counts[:, None] * variances + prior * pooled_var[None, :])
                  / (counts[:, None] + prior))
        scale = np.sqrt(np.maximum(shrunk, 1e-24))
        adjusted[test] = (matrix[test] - means[codes[test]]) / scale[codes[test]]
    return adjusted - adjusted.mean(axis=0, keepdims=True)


# --- the cross-fitting artefact, and the two controls that isolate it --------------------

def in_sample_residuals(matrix: np.ndarray, design: np.ndarray) -> np.ndarray:
    """Residuals of ``matrix`` on ``design`` fitted on **all** rows, no cross-fitting.

    Not an adjustment anyone should use for a channel: fitting the nuisance model on the
    rows whose residuals it produces removes more than the confound, which is the
    over-residualisation failure ``residualise.cross_fitted_residuals`` exists to avoid.
    It is here as a **diagnostic**, because it has one property no cross-fitted residual
    has: the residual mean of every design cell is **exactly zero** by the normal
    equations.

    That makes it the clean contrast for the question "is the confound recovery a
    first-moment effect?".  A cross-fitted residual leaves each (cell x fold) group
    displaced by that fold's estimation error ``mu_c - mu_hat_c^(-f)``, an offset of
    order ``sigma / sqrt(n_c)`` **shared by every patient in the group** -- a group mean
    shift that no cross-fitted adjustment of any flexibility can remove, because it is
    created by the cross-fitting itself.  If a k-NN's recovery is that offset, it
    survives the saturated arm and vanishes here.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    design = np.asarray(design, dtype=np.float64)
    if design.shape[1] == 0:
        return matrix - matrix.mean(axis=0, keepdims=True)
    fitted = design @ np.linalg.lstsq(design, matrix, rcond=None)[0]
    residual = matrix - fitted
    return residual - residual.mean(axis=0, keepdims=True)


def row_shuffled(matrix: np.ndarray, *, seed: int = 42) -> np.ndarray:
    """``matrix`` with its rows permuted -- the negative control for an adjustment artefact.

    Shuffling rows destroys every association between the block and the confound labels
    while preserving the block's marginal geometry exactly: the same rows, the same
    covariance, the same effective rank.  Push the shuffled block through the *same*
    design and the *same* adjuster and probe it with the *same* k-NN, and any recovery
    that remains cannot be surviving confound structure, because there is none left to
    survive.  It can only be structure the adjustment introduced.

    This is the control the previous run's defence needed and did not have.
    ``tcga_nonlinear_confound_probe_result_20260804T2100Z.md`` item 10 argues that
    "any structure the residualiser itself introduced is inside the null as well as
    inside the observed".  That is not so: the global null permutes the **labels** of the
    already-residualised features, so a cell-tied artefact is broken in the null and
    intact in the observed, and would be scored as signal.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    order = np.random.default_rng(int(seed)).permutation(len(matrix))
    return matrix[order]


def cross_fitting_offset_energy(residual: np.ndarray, codes: np.ndarray,
                                folds: list[tuple[np.ndarray, np.ndarray]] | None = None, *,
                                n_splits: int = 5, seed: int = 42) -> dict:
    """How much of a residual block's energy sits in per-cell and per-(cell x fold) means.

    ``cell_mean_energy_fraction`` is ``0`` for :func:`in_sample_residuals` by
    construction and strictly positive for any cross-fitted residual;
    ``cell_fold_mean_energy_fraction`` isolates the part that is tied to the residualiser's
    own folds and is therefore unambiguously an artefact of cross-fitting rather than a
    property of the representation.  Reported so that "the adjustment leaves a group mean
    offset" is a measured quantity and not an argument.
    """
    residual = np.asarray(residual, dtype=np.float64)
    codes = np.asarray(codes, dtype=np.int64)
    folds = folds if folds is not None else _folds(len(residual), n_splits, seed)
    total = float((residual ** 2).sum())
    n_cells = int(codes.max()) + 1 if codes.size else 0

    def _group_energy(group: np.ndarray) -> float:
        levels = np.unique(group)
        energy = 0.0
        for level in levels:
            rows = np.flatnonzero(group == level)
            energy += float(len(rows)) * float((residual[rows].mean(axis=0) ** 2).sum())
        return energy

    fold_of = np.empty(len(residual), dtype=np.int64)
    for index, (_, test) in enumerate(folds):
        fold_of[test] = index
    cell_fold = codes.astype(np.int64) * (len(folds) + 1) + fold_of
    # For G groups of a block with no group structure, the expected share of energy in the
    # group means is G / n -- each group's mean carries p * sigma^2 / n_g of energy over n_g
    # rows.  Quoting the raw fraction without it says nothing, because a block cut into 525
    # groups puts 19% of its energy in group means by arithmetic alone.
    expected = {"cell": float(len(np.unique(codes))) / float(len(residual)),
                "cell_fold": float(len(np.unique(cell_fold))) / float(len(residual))}
    observed = {"cell": _group_energy(codes) / total if total > 0 else float("nan"),
                "cell_fold": _group_energy(cell_fold) / total if total > 0 else float("nan")}
    return {
        "total_energy": total,
        "cell_mean_energy_fraction": observed["cell"],
        "cell_fold_mean_energy_fraction": observed["cell_fold"],
        "cell_mean_energy_expected_under_random_grouping": expected["cell"],
        "cell_fold_mean_energy_expected_under_random_grouping": expected["cell_fold"],
        "cell_mean_energy_ratio_to_expected": observed["cell"] / expected["cell"],
        "cell_fold_mean_energy_ratio_to_expected": observed["cell_fold"] / expected["cell_fold"],
        "n_cells": n_cells,
        "n_folds": int(len(folds)),
    }


def regenerated_adjustment_null(x: np.ndarray, class_index: np.ndarray, n_classes: int, adjust, *,
                                k_grid=(1, 3, 5, 10, 15, 25, 50), n_splits: int = 5,
                                seed: int = 42, n_permutations: int = 200, n_jobs: int = 1) -> dict:
    """The confound probe against a null that **regenerates the adjustment inside itself**.

    This is the null the k-NN confound probe needs and does not have, and the asymmetry is
    the whole point:

    * ``calibration.permutation_null`` -- the channel's null -- permutes the pairing and
      then **re-residualises** ``y`` on every permutation (``calibration.py:158``).  Any
      correlation the shared residualisation *induces* between two blocks is therefore
      regenerated inside the null, which is why P1 4.6's induced-correlation floor does
      not inflate the channel's ``permutation_p``.
    * ``nonlinear_confound_probe.probe_state`` permutes the **labels** of an
      already-adjusted feature matrix.  A structure that the adjustment tied to the *true*
      cells is broken in the null and intact in the observed, so it is scored as surviving
      confound rather than as chance.

    Cross-fitted residualisation leaves exactly such a structure: each (cell x fold) group
    is displaced by that fold's estimation error, an offset shared by every patient in the
    group.  Here the null permutes the **rows of the block before the adjustment**, so
    every draw carries an adjustment artefact of its own size and the observed value is
    graded against it.  The labels, the design, the cell structure and the probe folds are
    untouched; only the association between the block and the labels is destroyed.

    The reading is the maximum over ``k_grid`` and over both vote rules, computed the same
    way in the observed and in every permutation, so the null is a null *for the statistic
    being quoted* rather than for one arbitrary member of the grid.
    """
    from .confound_certificate import _stratified_folds
    from .nonlinear_confound_probe import knn_balanced_accuracy_oof, null_summary

    x = np.asarray(x, dtype=np.float64)
    class_index = np.asarray(class_index, dtype=np.int64)
    folds = _stratified_folds(class_index, n_splits, seed)
    grid = [(int(k), bool(prior)) for k in k_grid for prior in (False, True)]

    def _reading(block: np.ndarray) -> float:
        return max(knn_balanced_accuracy_oof(block, class_index, n_classes, k=k,
                                             n_splits=n_splits, seed=seed, folds=folds,
                                             prior_correction=prior)
                   for k, prior in grid)

    adjusted = adjust(x)
    observed = float(_reading(adjusted))
    per_probe = {f"knn_{'prior_' if prior else ''}k{k}":
                 float(knn_balanced_accuracy_oof(adjusted, class_index, n_classes, k=k,
                                                 n_splits=n_splits, seed=seed, folds=folds,
                                                 prior_correction=prior))
                 for k, prior in grid}

    rng = np.random.default_rng(seed)
    orders = [(rng.permutation(len(x)),) for _ in range(n_permutations)]

    def _one(order):
        return _reading(adjust(x[order]))

    record = {"observed": observed, "multiple_of_chance": observed * float(n_classes),
              "chance_rate": 1.0 / float(n_classes), "n_classes": int(n_classes),
              "per_probe_observed": per_probe, "k_grid": [int(k) for k in k_grid]}
    if n_permutations > 0:
        null = np.asarray(_map(_one, orders, n_jobs), dtype=np.float64)
        record["regenerated_null"] = null_summary(null, observed)
        record["regenerated_null"]["null_median_multiple_of_chance"] = (
            float(np.median(null)) * float(n_classes))
        record["excess_multiple"] = (observed - float(np.median(null))) * float(n_classes)
    return record


# --- arm construction --------------------------------------------------------------------

def make_adjuster(model: str, *, design: np.ndarray | None = None,
                  codes: np.ndarray | None = None, n_splits: int = 5, seed: int = 42,
                  alpha: float = 1.0, gamma: float = 0.5, n_estimators: int = 300,
                  prior_count: float = 10.0, n_jobs: int = 1):
    """Build one adjustment arm as a ``matrix -> residual`` callable.

    ``model="ridge"`` returns a closure over ``cross_fitted_residuals`` itself rather than
    a reimplementation, so the incumbent arm is byte-identical to every published number
    by construction rather than by assertion (and ``test_nonlinear_adjustment.py`` asserts
    it anyway).
    """
    if model == "none":
        return lambda matrix: (np.asarray(matrix, dtype=np.float64)
                               - np.asarray(matrix, dtype=np.float64).mean(axis=0, keepdims=True))
    if model == "ridge":
        if design is None:
            raise ValueError("ridge adjustment needs a design")
        return lambda matrix: cross_fitted_residuals(matrix, design, n_splits=n_splits,
                                                     alpha=alpha, seed=seed)
    if model == "saturated":
        if codes is None:
            raise ValueError("saturated adjustment needs cell codes")
        return lambda matrix: saturated_cell_residuals(matrix, codes, n_splits=n_splits,
                                                       alpha=alpha, seed=seed)
    if model == "kernel_ridge":
        if design is None:
            raise ValueError("kernel ridge adjustment needs a design")
        return KernelRidgeAdjuster(design, alpha=alpha, gamma=gamma, n_splits=n_splits, seed=seed)
    if model == "forest":
        if design is None:
            raise ValueError("forest adjustment needs a design")
        return lambda matrix: forest_residuals(matrix, design, n_estimators=n_estimators,
                                               n_splits=n_splits, seed=seed, n_jobs=n_jobs)
    if model == "location_scale":
        if codes is None:
            raise ValueError("location-scale adjustment needs cell codes")
        return lambda matrix: cross_fitted_location_scale(matrix, codes, n_splits=n_splits,
                                                          seed=seed, prior_count=prior_count)
    if model == "in_sample_additive":
        if design is None:
            raise ValueError("in-sample additive adjustment needs a design")
        return lambda matrix: in_sample_residuals(matrix, design)
    if model == "in_sample_saturated":
        if codes is None:
            raise ValueError("in-sample saturated adjustment needs cell codes")
        return lambda matrix: in_sample_residuals(matrix, cell_design(codes))
    raise ValueError(f"unknown adjustment model {model!r}; expected one of {ADJUSTER_MODELS}")


def adjuster_agreement(a: np.ndarray, b: np.ndarray) -> dict:
    """How much a candidate adjustment differs from the incumbent, per axis.

    Predeclared distrust item 2: an arm whose residuals correlate with the incumbent's at
    >= 0.99 on every axis is a relabelled incumbent, and its channel reading says nothing
    about nonlinearity whatever it is.  Reported for every arm so the reader does not
    have to take "this arm is different" on trust.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError("adjuster_agreement compares two adjustments of the same block")
    per_axis = np.full(a.shape[1], np.nan)
    for j in range(a.shape[1]):
        if np.std(a[:, j]) > 1e-12 and np.std(b[:, j]) > 1e-12:
            per_axis[j] = float(np.corrcoef(a[:, j], b[:, j])[0, 1])
    frobenius = float(np.linalg.norm(a - b) / max(np.linalg.norm(a), 1e-12))
    return {"per_axis_correlation_median": float(np.nanmedian(per_axis)),
            "per_axis_correlation_min": float(np.nanmin(per_axis)) if np.isfinite(per_axis).any() else float("nan"),
            "per_axis_correlation_p05": float(np.nanpercentile(per_axis, 5)) if np.isfinite(per_axis).any() else float("nan"),
            "relative_frobenius_difference": frobenius,
            "is_relabelled_incumbent": bool(np.nanmedian(per_axis) >= 0.99)}


def cross_fitted_r2(matrix: np.ndarray, adjust) -> dict:
    """Fraction of a block's centred variance an adjustment removes, out of fold.

    ``1 - ||residual||^2 / ||centred||^2``, pooled over columns and per column.  Used for
    the labels-only variance accounting of the predeclaration's step 2: it says directly
    how much of each side of the channel *is* the confound.  Only meaningful for
    mean-removing arms; a scale-changing arm (``location_scale``) can report a negative
    value and that is not a defect, so the value is emitted without clipping.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    centred = matrix - matrix.mean(axis=0, keepdims=True)
    residual = np.asarray(adjust(matrix), dtype=np.float64)
    total = float((centred ** 2).sum())
    per_column_total = (centred ** 2).sum(axis=0)
    per_column = np.where(per_column_total > 1e-24,
                          1.0 - (residual ** 2).sum(axis=0) / np.maximum(per_column_total, 1e-24),
                          np.nan)
    return {"pooled_r2": float(1.0 - (residual ** 2).sum() / total) if total > 1e-24 else float("nan"),
            "per_column_r2_median": float(np.nanmedian(per_column)),
            "per_column_r2_max": float(np.nanmax(per_column)) if np.isfinite(per_column).any() else float("nan")}


# --- the channel, measured through an arbitrary adjustment -------------------------------

def channel_under_adjustment(x: np.ndarray, y: np.ndarray, adjust, *, strata=None,
                             n_permutations: int = 2000, n_components: int = 16, seed: int = 42,
                             n_jobs: int = 1, with_heldout: bool = True) -> dict:
    """P1's channel statistic and its pairing null, under any adjustment.

    Generalises ``calibration.permutation_null`` in exactly one respect: the fixed
    ``cross_fitted_residuals(., design)`` becomes a caller-supplied ``adjust``.  Every
    other step -- which block is residualised once, which is re-residualised per
    permutation, how the permutation orders are drawn, the p-value convention -- is
    unchanged, and ``test_nonlinear_adjustment.py`` asserts that with the incumbent
    adjuster this function reproduces ``permutation_null``'s dictionary exactly.  That
    identity is the whole basis on which an arm's number is comparable to P1 4.4's.

    ``strata=None`` gives the **global** pairing null; passing cancer type gives the
    project's within-cancer convention.  Note that the nesting degeneracy which makes a
    within-cancer *label* null uninterpretable on TCGA does not arise here: this
    permutes the patient **pairing** between two blocks, not a label inside a stratum.
    Both are computed by the CLI and the within-cancer one is quoted, as P1 4.4 does.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    rng = np.random.default_rng(seed)
    x_adjusted = adjust(x)
    strata = np.zeros(len(x), dtype=int) if strata is None else np.asarray(strata)

    y_adjusted = adjust(y)
    observed = top_canonical_correlation(x_adjusted, y_adjusted, n_components=n_components)
    orders = []
    for _ in range(n_permutations):
        order = np.arange(len(y))
        for level in np.unique(strata):
            idx = np.flatnonzero(strata == level)
            order[idx] = rng.permutation(idx)
        orders.append((order,))

    def _one_permutation(order):
        return top_canonical_correlation(x_adjusted, adjust(y[order]), n_components=n_components)

    record = {"observed_top_cca": float(observed), "n_permutations": int(n_permutations),
              "n_components": int(n_components)}
    if n_permutations > 0:
        null = np.asarray(_map(_one_permutation, orders, n_jobs), dtype=np.float64)
        exceed = int(np.sum(null >= observed))
        record.update({
            "null_median": float(np.median(null)),
            "null_p95": float(np.percentile(null, 95)),
            "null_max": float(np.max(null)),
            "excess_over_null_median": float(observed - np.median(null)),
            "permutation_p": float((exceed + 1) / (n_permutations + 1)),
            "permutation_resolution": 1.0 / (n_permutations + 1),
        })
    if with_heldout:
        record["heldout_top_cca"] = float(heldout_top_cca(x_adjusted, y_adjusted,
                                                          n_components=n_components, seed=seed))
    record.update({
        "effective_rank_x_adjusted": float(effective_rank(x_adjusted)),
        "effective_rank_y_adjusted": float(effective_rank(y_adjusted)),
        "unadjusted_top_cca": float(top_canonical_correlation(x, y, n_components=n_components)),
    })
    return record


def retention_of_excess(arm: dict, reference: dict) -> float:
    """``excess(arm) / excess(reference)`` -- how much of the channel an arm leaves standing.

    The predeclared grading quantity, and it is a ratio of **excesses over each arm's own
    null**, never of raw ``S1`` values.  A stronger adjustment can move the capacity floor
    as well as the signal, so comparing a raw 0.60 against a null measured under a
    different adjustment is exactly the error P1 4.4 warns about.  Kept here rather than
    written out in a notebook entry so that every quoted retention comes from one
    definition.
    """
    denominator = float(reference["excess_over_null_median"])
    if abs(denominator) < 1e-12:
        return float("nan")
    return float(arm["excess_over_null_median"]) / denominator


def corrected_multiple(record: dict) -> float:
    """A probe reading as a multiple of the **measured** chance rate under its adjustment.

    ``multiple_of_chance`` divides by the design rate ``1 / n_classes``.  That is the right
    denominator for a raw block and the wrong one for an adjusted one, because a
    cross-fitted adjustment leaves cell-tied structure that a k-NN reads and the
    label-permutation null cannot see.  This divides by the median of
    :func:`regenerated_adjustment_null` instead, which measures the rate rather than
    assuming it.
    """
    median = float(record["regenerated_null"]["null_median"])
    if median <= 0:
        return float("nan")
    return float(record["observed"]) / median


def labels_only_ceiling(design: np.ndarray, y: np.ndarray, *, n_components: int = 16,
                        seed: int = 42, adjust=None, strata=None, n_permutations: int = 0,
                        n_jobs: int = 1) -> dict:
    """The channel obtainable from the confound labels alone, with no image features.

    The bound from the other side.  ``design`` stands in for the image representation --
    either the additive one-hot design or the saturated cell design -- and the same
    16-component statistic is read, so capacity is matched to the real measurement rather
    than assumed comparable.  With ``adjust`` supplied the labels are pushed through the
    identical pipeline as well, which answers the sharper question: what does a purely
    confound "representation" score *after* CALIBRA's adjustment?
    """
    design = np.asarray(design, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    record = {
        "raw_top_cca": float(top_canonical_correlation(design, y, n_components=n_components)),
        "raw_heldout_top_cca": float(heldout_top_cca(design, y, n_components=n_components,
                                                     seed=seed)),
        "n_design_columns": int(design.shape[1]),
        "design_rank": int(np.linalg.matrix_rank(design)),
    }
    if adjust is not None:
        record["adjusted"] = channel_under_adjustment(
            design, y, adjust, strata=strata, n_permutations=n_permutations,
            n_components=n_components, seed=seed, n_jobs=n_jobs)
    return record


# --- CLI ---------------------------------------------------------------------------------

def _load_block(artifact: Path, targets: Path, state: str, partition: str, min_site_count: int):
    """Rebuild ``run_calibra``'s x / y / design / strata for one artifact and state.

    Deliberately mirrors ``run_calibra.main`` line for line -- the same target artifact
    validation, the same ``RANDOM_CONTROL__`` exclusion, the same partition mask, the same
    alignment to the molecular patient index, the same pooling -- because a channel number
    measured on a differently-assembled cohort is not comparable to P1 4.4's.
    """
    import pandas as pd

    raw = np.load(artifact, allow_pickle=True)
    targets_raw = np.load(targets, allow_pickle=True)
    target_ids = np.asarray([str(p) for p in targets_raw["patient_ids"]])
    target_names = np.asarray([str(t) for t in targets_raw["target_names"]])
    scores = np.asarray(targets_raw["scores"], dtype=np.float64)
    keep = ~np.char.startswith(target_names, "RANDOM_CONTROL__")
    scores, target_names = scores[:, keep], target_names[keep]
    target_index = {pid: i for i, pid in enumerate(target_ids)}

    patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
    cancers = np.asarray([str(c) for c in raw["cancers"]])
    split = np.asarray([str(s) for s in raw["split"]])
    mask = np.ones(len(patient_ids), bool) if partition == "all" else (split == partition)
    aligned = np.asarray([target_index.get(pid, -1) for pid in patient_ids])
    mask &= aligned >= 0

    y = scores[aligned[mask]]
    x = np.asarray(raw[state], dtype=np.float64)[mask]
    tss, frequent = pooled_tissue_source_site(patient_ids[mask], min_site_count=min_site_count)
    frame = pd.DataFrame({"cancer": cancers[mask], "tss": tss})
    design = confound_design(frame, ["cancer", "tss"])
    codes, levels = cell_codes(cancers[mask], tss)
    return {"x": x, "y": y, "design": design, "codes": codes, "cell_levels": levels,
            "cancers": cancers[mask], "tss": tss, "patient_ids": patient_ids[mask],
            "n": int(mask.sum()), "n_targets": int(y.shape[1]), "n_sites_kept": int(len(frequent))}


def _arm_specs(arms: list[str], kernel_grid: str, forest_trees: int, prior_count: float):
    """Expand the arm list into concrete (name, model, kwargs) triples."""
    specs = []
    for arm in arms:
        if arm == "kernel_ridge":
            for cell in str(kernel_grid).split(","):
                if not cell.strip():
                    continue
                alpha, gamma = (float(v) for v in cell.split(":"))
                specs.append((f"kernel_ridge_a{alpha:g}_g{gamma:g}", "kernel_ridge",
                              {"alpha": alpha, "gamma": gamma}))
        elif arm == "forest":
            specs.append(("forest", "forest", {"n_estimators": int(forest_trees)}))
        elif arm == "saturated":
            specs.append(("saturated", "saturated", {"alpha": 1e-6}))
        elif arm == "location_scale":
            specs.append(("location_scale", "location_scale", {"prior_count": float(prior_count)}))
        elif arm in ("none", "ridge", "in_sample_additive", "in_sample_saturated"):
            specs.append((arm, arm, {}))
        else:
            raise ValueError(f"unknown arm {arm!r}")
    return specs


def _run_probe(features: np.ndarray, block: dict, k_grid, args) -> dict:
    """The k-NN confound probe on one already-adjusted block, for both targets.

    Imported wholesale from ``nonlinear_confound_probe`` -- same estimator, same fold
    rule, same nulls as the run that produced the 4.3-4.9x readings this run is testing.
    The **global** null is the applicable bar because TCGA site nests inside cancer; the
    within-cancer null is carried beside the site rows for decomposition only.
    """
    from .confound_certificate import _encode
    from .nonlinear_confound_probe import nesting_diagnostic, probe_state

    probes = {}
    for target, labels, strata in (("site", block["tss"], block["cancers"]),
                                   ("cancer", block["cancers"], None)):
        classes, class_index = _encode(np.asarray(labels).astype(str))
        probe = probe_state(features, class_index, int(len(classes)), strata=strata,
                            k_grid=k_grid, forest_trees=0, run_svm=False,
                            n_splits=args.n_splits, seed=args.seed,
                            knn_permutations=args.probe_permutations, model_permutations=0,
                            n_jobs=args.n_jobs,
                            nulls="both" if strata is not None else "global")
        if target == "site":
            probe["nesting"] = nesting_diagnostic(class_index, strata, int(len(classes)))
        probes[target] = probe
    return probes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tag", default="nonlinear_adjustment")
    parser.add_argument("--states", nargs="+", default=("wsi_biology",))
    parser.add_argument("--partition", default="test", choices=("test", "val", "all"))
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--arms", nargs="+",
                        default=("none", "ridge", "saturated", "kernel_ridge", "location_scale"))
    parser.add_argument("--kernel-grid", default="1.0:0.25,1.0:0.5,1.0:1.0,0.1:0.25,0.1:0.5,0.1:1.0")
    parser.add_argument("--forest-trees", type=int, default=300)
    parser.add_argument("--prior-count", type=float, default=10.0)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-permutations", type=int, default=2000,
                        help="pairing permutations for the channel null")
    parser.add_argument("--global-null", action="store_true",
                        help="ALSO compute the unstratified pairing null beside the within-cancer one")
    parser.add_argument("--probe", action="store_true", help="step 3: re-run the k-NN confound probe")
    parser.add_argument("--probe-permutations", type=int, default=200)
    parser.add_argument("--probe-k-grid", default="1,3,5,10,15,25,50")
    parser.add_argument("--ceiling", action="store_true", help="step 2: the labels-only ceiling")
    parser.add_argument("--regenerated-null", action="store_true",
                        help="grade the confound probe against a null that REGENERATES the "
                             "adjustment inside each permutation, as the channel's null does")
    parser.add_argument("--shuffle-control", action="store_true",
                        help="also probe the ROW-SHUFFLED block through the same adjuster: the "
                             "negative control for structure the adjustment itself introduces")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    specs = _arm_specs(list(args.arms), args.kernel_grid, args.forest_trees, args.prior_count)
    results: dict = {"config": {k: (list(v) if isinstance(v, (list, tuple)) else v)
                                for k, v in vars(args).items()},
                     "arms": [name for name, _, _ in specs]}

    def _flush():
        (output / f"{args.tag}.json").write_text(json.dumps(results, indent=2), encoding="utf-8")

    for artifact_path in args.artifacts:
        path = Path(artifact_path)
        for state in args.states:
            block = _load_block(path, Path(args.targets), state, args.partition, args.min_site_count)
            key_prefix = f"{path.stem}::{state}"
            results[f"{key_prefix}::cohort"] = {
                "n": block["n"], "n_targets": block["n_targets"],
                "n_design_columns": int(block["design"].shape[1]),
                "design_rank": int(np.linalg.matrix_rank(block["design"])),
                "n_cells": int(len(block["cell_levels"])),
                "n_sites_kept": block["n_sites_kept"],
                "duplicate_patient_ids": int(len(block["patient_ids"])
                                             - len(set(block["patient_ids"].tolist()))),
                "n_axes": int(block["x"].shape[1]),
            }
            _flush()

            incumbent = make_adjuster("ridge", design=block["design"], n_splits=args.n_splits,
                                      seed=args.seed)
            incumbent_x = incumbent(block["x"])

            if args.ceiling:
                for label, matrix in (("additive_design", block["design"]),
                                      ("saturated_cell_design", cell_design(block["codes"]))):
                    started = time.time()
                    record = labels_only_ceiling(
                        matrix, block["y"], n_components=args.n_components, seed=args.seed,
                        adjust=incumbent, strata=block["cancers"],
                        n_permutations=args.n_permutations, n_jobs=args.n_jobs)
                    record["r2_of_targets_on_labels"] = cross_fitted_r2(
                        block["y"], make_adjuster("saturated", codes=block["codes"],
                                                  n_splits=args.n_splits, seed=args.seed,
                                                  alpha=1e-6))
                    record["r2_of_image_on_labels"] = cross_fitted_r2(
                        block["x"], make_adjuster("saturated", codes=block["codes"],
                                                  n_splits=args.n_splits, seed=args.seed,
                                                  alpha=1e-6))
                    record["wall_seconds"] = round(time.time() - started, 1)
                    results[f"{key_prefix}::ceiling::{label}"] = record
                    print(f"[ceiling {key_prefix} {label}] raw_top_cca="
                          f"{record['raw_top_cca']:.4f} heldout={record['raw_heldout_top_cca']:.4f} "
                          f"adjusted={record['adjusted']['observed_top_cca']:.4f}", flush=True)
                    _flush()

            for name, model, kwargs in specs:
                started = time.time()
                adjust = make_adjuster(model, design=block["design"], codes=block["codes"],
                                       n_splits=args.n_splits, seed=args.seed,
                                       n_jobs=args.n_jobs, **kwargs)
                record = channel_under_adjustment(
                    block["x"], block["y"], adjust, strata=block["cancers"],
                    n_permutations=args.n_permutations, n_components=args.n_components,
                    seed=args.seed, n_jobs=args.n_jobs)
                record["null_strata"] = "cancer"
                if args.global_null:
                    record["global_pairing_null"] = channel_under_adjustment(
                        block["x"], block["y"], adjust, strata=None,
                        n_permutations=args.n_permutations, n_components=args.n_components,
                        seed=args.seed, n_jobs=args.n_jobs, with_heldout=False)
                adjusted_x = adjust(block["x"])
                record["agreement_with_incumbent"] = adjuster_agreement(adjusted_x, incumbent_x)
                record["r2_of_image_removed"] = cross_fitted_r2(block["x"], adjust)
                record["r2_of_targets_removed"] = cross_fitted_r2(block["y"], adjust)
                record["arm"] = name
                record["model"] = model
                record["params"] = kwargs
                record["wall_seconds"] = round(time.time() - started, 1)

                record["offset_energy"] = cross_fitting_offset_energy(
                    adjusted_x, block["codes"], n_splits=args.n_splits, seed=args.seed)

                if args.regenerated_null:
                    from .confound_certificate import _encode
                    k_grid = tuple(int(v) for v in str(args.probe_k_grid).split(",") if v.strip())
                    regenerated = {}
                    for target, labels in (("site", block["tss"]), ("cancer", block["cancers"])):
                        classes, class_index = _encode(np.asarray(labels).astype(str))
                        regenerated[target] = regenerated_adjustment_null(
                            block["x"], class_index, int(len(classes)), adjust, k_grid=k_grid,
                            n_splits=args.n_splits, seed=args.seed,
                            n_permutations=args.probe_permutations, n_jobs=args.n_jobs)
                        r = regenerated[target]["regenerated_null"]
                        print("      [regen %s::%s %s] obs=%.4f (%.2fx) null_med=%.4f (%.2fx) "
                              "p95=%.4f p=%.4f" % (
                                  key_prefix, name, target, regenerated[target]["observed"],
                                  regenerated[target]["multiple_of_chance"], r["null_median"],
                                  r["null_median_multiple_of_chance"], r["null_p95"],
                                  r["permutation_p"]), flush=True)
                    record["regenerated_adjustment_null"] = regenerated

                if args.probe:
                    k_grid = tuple(int(v) for v in str(args.probe_k_grid).split(",") if v.strip())
                    record["probe"] = _run_probe(adjusted_x, block, k_grid, args)
                    if args.shuffle_control:
                        # The block's rows permuted BEFORE the adjustment: same geometry,
                        # no confound association left. Anything the probe still reads
                        # here was put there by the adjustment.
                        shuffled = row_shuffled(block["x"], seed=args.seed)
                        record["shuffle_control"] = {
                            "adjusted": _run_probe(adjust(shuffled), block, k_grid, args),
                            "raw": _run_probe(shuffled - shuffled.mean(axis=0, keepdims=True),
                                              block, k_grid, args),
                        }

                results[f"{key_prefix}::{name}"] = record
                observed = record["observed_top_cca"]
                median = record.get("null_median", float("nan"))
                print(f"[{key_prefix}::{name}] top_cca={observed:.4f} null_median={median:.4f} "
                      f"excess={observed - median:.4f} heldout={record['heldout_top_cca']:.4f} "
                      f"rank={record['effective_rank_x_adjusted']:.2f} "
                      f"({record['wall_seconds']:.0f} s)", flush=True)
                _flush()
    _flush()
    print(f"[done] -> {output / (args.tag + '.json')}", flush=True)


if __name__ == "__main__":
    main()
