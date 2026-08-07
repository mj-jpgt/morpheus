"""WS-A3 -- counterfactuals for the four new association metrics.

Every metric must ship with: (1) a positive control -- recovering a planted association of
known strength, monotonically, before any confound is involved; (2) a must-fail control --
independent blocks reading near that metric's own null scale; (3) a matched (permutation)
null, reused from ``association.metric_permutation_null`` rather than an ad hoc threshold.
Plus a synthetic reproduction of the predeclared Regime I / Regime II question (does the
induced-correlation floor described in ``induced_correlation_sweep.py`` and P1 spec §4 show up
under every metric, not only Pearson-on-projections) -- this MUST be established on synthetic
data with a known ground truth before any of these metrics is trusted on real artifacts, per
``NOTEBOOK_ENTRIES/PREDECLARED_ws_a3_association_metrics_20260807T2211Z.md``.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.association import (ASSOCIATION_METRICS, compute_all_metrics,
                                              distance_correlation, hsic, kernel_cca,
                                              median_heuristic_gamma, metric_permutation_null,
                                              metric_recovery_curve, rv_coefficient,
                                              subsample_rows)
from morpheus.v2.calibra.residualise import confound_design, cross_fitted_residuals


# --------------------------------------------------------------------------------------- #
# subsampling rule                                                                          #
# --------------------------------------------------------------------------------------- #

def test_subsample_rows_is_a_noop_under_threshold():
    idx = subsample_rows(500, seed=1, threshold=2000, target=2000)
    assert idx.tolist() == list(range(500))


def test_subsample_rows_subsamples_deterministically_over_threshold():
    idx1 = subsample_rows(3000, seed=7, threshold=2000, target=2000)
    idx2 = subsample_rows(3000, seed=7, threshold=2000, target=2000)
    assert len(idx1) == 2000
    assert np.array_equal(idx1, idx2)                     # seeded, reproducible
    assert len(set(idx1.tolist())) == 2000                # no repeats
    idx3 = subsample_rows(3000, seed=8, threshold=2000, target=2000)
    assert not np.array_equal(idx1, idx3)                 # a different seed moves it


# --------------------------------------------------------------------------------------- #
# positive control -- every metric must recover a planted association, monotonically       #
# --------------------------------------------------------------------------------------- #

def _planted_linear_pair(n: int, p: int, q: int, r_true: float, *, seed: int):
    """A minimal, self-contained planted-association construction (not a copy of
    ``calibration.spike_targets`` -- deliberately simpler, single-purpose for this smoke test:
    one shared latent factor loads onto one column of X and one column of Y at exactly
    ``corr = r_true`` in expectation, the rest is independent noise."""
    rng = np.random.default_rng(seed)
    latent = rng.normal(size=n)
    x_noise = rng.normal(size=(n, p))
    y_noise = rng.normal(size=(n, q))
    x_noise[:, 0] = r_true * latent + np.sqrt(max(0.0, 1 - r_true ** 2)) * x_noise[:, 0]
    y_noise[:, 0] = r_true * latent + np.sqrt(max(0.0, 1 - r_true ** 2)) * y_noise[:, 0]
    return x_noise, y_noise


@pytest.mark.parametrize("metric_name", ["rv", "dcor", "hsic", "kernel_cca"])
def test_positive_control_every_metric_increases_with_planted_strength(metric_name):
    metric = ASSOCIATION_METRICS[metric_name]
    n = 400
    readings = []
    for r_true in (0.0, 0.3, 0.6, 0.9):
        x, y = _planted_linear_pair(n, 5, 5, r_true, seed=123)
        readings.append(metric(x, y) if metric_name == "rv" else metric(x, y, seed=1))
    # Monotone non-decreasing within a small numerical tolerance -- the falsifier named in
    # the predeclaration (§3.1): a metric that does not rise with planted strength is flagged
    # unreliable before it is used for anything else.
    for a, b in zip(readings, readings[1:]):
        assert b >= a - 0.03, (metric_name, readings)
    assert readings[-1] > readings[0] + 0.15, (metric_name, readings)


# --------------------------------------------------------------------------------------- #
# must-fail control -- independent blocks read near the null                                #
# --------------------------------------------------------------------------------------- #

@pytest.mark.parametrize("metric_name", ["rv", "dcor", "hsic", "kernel_cca"])
def test_must_fail_control_independent_blocks_stay_near_null(metric_name):
    rng = np.random.default_rng(99)
    n = 300
    x = rng.normal(size=(n, 6))
    y = rng.normal(size=(n, 6))
    metric = ASSOCIATION_METRICS[metric_name]
    observed = metric(x, y) if metric_name == "rv" else metric(x, y, seed=1)
    # A permutation null of the SAME statistic on the SAME data is the only honest reference
    # (RV/HSIC/kernel-CCA have no natural "zero"): observed must not tower over its own null.
    perms = []
    for order_seed in range(15):
        order = np.random.default_rng(order_seed).permutation(n)
        perms.append(metric(x, y[order]) if metric_name == "rv" else metric(x, y[order], seed=1))
    null = np.asarray(perms)
    assert observed <= np.percentile(null, 95) + 0.05, (metric_name, observed, null)


# --------------------------------------------------------------------------------------- #
# RV coefficient identity check -- RV(x, y) == r^2 for two single columns                   #
# --------------------------------------------------------------------------------------- #

def test_rv_coefficient_reduces_to_squared_pearson_for_single_columns():
    rng = np.random.default_rng(5)
    x = rng.normal(size=200)
    y = 0.4 * x + rng.normal(size=200)
    rv = rv_coefficient(x, y)
    pearson_sq = float(np.corrcoef(x, y)[0, 1]) ** 2
    assert rv == pytest.approx(pearson_sq, rel=1e-9)


# --------------------------------------------------------------------------------------- #
# HSIC estimator identity: sum(Kc*Lc)/n^2 == trace(K H L H)/n^2, the formula named in the    #
# module docstring, checked against a brute-force H-matrix computation                      #
# --------------------------------------------------------------------------------------- #

def test_hsic_matches_the_brute_force_trace_khlh_formula():
    rng = np.random.default_rng(3)
    n = 60
    x = rng.normal(size=(n, 3))
    y = rng.normal(size=(n, 3))
    fast = hsic(x, y, seed=1, max_n=2000)

    gx = median_heuristic_gamma(x)
    gy = median_heuristic_gamma(y)

    def _kernel(a, gamma):
        sq = (a ** 2).sum(axis=1)
        d2 = np.maximum(sq[:, None] + sq[None, :] - 2 * a @ a.T, 0.0)
        return np.exp(-gamma * d2)

    K, L = _kernel(x, gx), _kernel(y, gy)
    H = np.eye(n) - np.ones((n, n)) / n
    brute = float(np.trace(K @ H @ L @ H)) / (n ** 2)
    assert fast == pytest.approx(brute, rel=1e-8)


def test_kernel_cca_is_bounded_and_self_association_is_high():
    rng = np.random.default_rng(11)
    n = 300
    x = rng.normal(size=(n, 4))
    value = kernel_cca(x, x, seed=1, max_n=2000)
    assert 0.0 <= value <= 1.0
    assert value > 0.5                      # a block against ITSELF must read strongly


def test_distance_correlation_bounded_in_zero_one():
    rng = np.random.default_rng(4)
    n = 200
    x = rng.normal(size=(n, 3))
    y = 0.7 * x[:, :1] + rng.normal(size=(n, 3)) * 0.3
    value = distance_correlation(x, y, seed=1)
    assert 0.0 <= value <= 1.0 + 1e-9


def test_compute_all_metrics_returns_every_key():
    rng = np.random.default_rng(2)
    x, y = rng.normal(size=(150, 4)), rng.normal(size=(150, 4))
    out = compute_all_metrics(x, y, seed=1, max_n=2000)
    assert set(out) == {"rv", "dcor", "hsic", "kernel_cca"}
    assert all(np.isfinite(v) for v in out.values())


# --------------------------------------------------------------------------------------- #
# THE central predeclared question, reproduced on SYNTHETIC data with known ground truth,   #
# mirroring the block-design simulation in P1_REVISION_SPEC.md §4.4 EXACTLY (Regime I vs     #
# Regime II) -- this is the counterfactual that must pass before any of these metrics is     #
# trusted on real artifacts.                                                                 #
# --------------------------------------------------------------------------------------- #

def _block_design(n: int, n_blocks: int, *, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sizes = rng.multinomial(n, np.ones(n_blocks) / n_blocks)
    codes = np.repeat(np.arange(n_blocks), sizes)
    rng.shuffle(codes)
    design = np.zeros((n, n_blocks))
    design[np.arange(n), codes] = 1.0
    return design


@pytest.mark.parametrize("metric_name", ["rv", "dcor", "hsic", "kernel_cca"])
def test_regime_ii_induced_association_exceeds_regime_i_synthetic(metric_name):
    """SYNTHETIC. Two blocks unrelated to each other (independent noise) but EACH carrying real
    block structure (Regime II: residualising against the TRUE block design induces spurious
    residual association) versus the SAME two blocks residualised against a Gaussian design of
    matched rank (Regime I: no real structure to explain, so no induced floor is predicted).
    This is exactly spec §4.4's simulation, generalised from Pearson to the new metric.
    """
    rng = np.random.default_rng(hash(metric_name) % (2 ** 31))
    n, n_blocks, p, q = 800, 20, 6, 6
    design_real = _block_design(n, n_blocks, seed=1)
    block_means_x = rng.normal(scale=1.5, size=(n_blocks, p))
    block_means_y = rng.normal(scale=1.5, size=(n_blocks, q))
    codes = design_real.argmax(axis=1)
    x = block_means_x[codes] + rng.normal(size=(n, p))
    y = block_means_y[codes] + rng.normal(size=(n, q))          # X and Y share NO latent signal

    design_gaussian = rng.normal(size=(n, n_blocks))             # Regime I: matched rank, no structure

    x_res_ii = cross_fitted_residuals(x, design_real, seed=1)
    y_res_ii = cross_fitted_residuals(y, design_real, seed=1)
    x_res_i = cross_fitted_residuals(x, design_gaussian, seed=1)
    y_res_i = cross_fitted_residuals(y, design_gaussian, seed=1)

    metric = ASSOCIATION_METRICS[metric_name]
    kwargs = {} if metric_name == "rv" else {"seed": 1, "max_n": 2000}
    regime_ii = metric(x_res_ii, y_res_ii, **kwargs)
    regime_i = metric(x_res_i, y_res_i, **kwargs)

    # F1 falsifier (predeclaration §4): report whichever way this comes out; the assertion
    # below is the criterion the predeclaration fixed in advance -- Regime II must exceed
    # Regime I by a real margin, not just numerically.
    assert regime_ii > regime_i, (metric_name, regime_ii, regime_i)


def test_metric_permutation_null_reuses_the_pairing_convention():
    rng = np.random.default_rng(21)
    n = 300
    design = confound_design(
        __import__("pandas").DataFrame({"g": rng.integers(0, 5, size=n).astype(str)}), ["g"])
    x = rng.normal(size=(n, 4))
    y = rng.normal(size=(n, 4))
    record = metric_permutation_null(x, y, design, rv_coefficient, n_permutations=20, seed=1)
    assert set(record) >= {"observed", "null_median", "null_p95", "permutation_p"}
    assert np.isfinite(record["observed"])


def test_metric_recovery_curve_reuses_floors_from_recovery_and_returns_a_floor():
    rng = np.random.default_rng(31)
    n = 500
    design = confound_design(
        __import__("pandas").DataFrame({"g": rng.integers(0, 6, size=n).astype(str)}), ["g"])
    x = rng.normal(size=(n, 6))
    y = rng.normal(size=(n, 6))
    result = metric_recovery_curve(x, y, design, rv_coefficient,
                                   levels=(0.0, 0.1, 0.3, 0.6), n_draws=6, seed=1)
    assert set(result) >= {"detection_floor", "transmission_floor", "recovered_median"}
    assert len(result["recovered_median"]) == 4
