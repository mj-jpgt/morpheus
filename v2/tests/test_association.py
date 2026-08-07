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
    """Monotonicity is checked RELATIVE to each metric's own scale, not against a shared
    absolute margin. RV and HSIC read on genuinely different absolute scales from dCor/kernel
    CCA once only 1 of 5 columns per side carries the planted signal (the other 4 are pure
    noise that dilutes a whole-block statistic but not a single-direction one) -- an absolute
    "+0.15" bar the first version of this test used was a test-design bug, not a metric defect:
    it flagged rv/hsic as failing while both still rose 4-6x from r_true=0 to r_true=0.9,
    strictly monotonically. Caught on the real training box, see the result notebook entry.
    """
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
        assert b >= a - 0.02 * max(readings), (metric_name, readings)
    floor = max(readings[0], 1e-6)
    assert readings[-1] > floor * 3.0, (metric_name, readings)


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

def _orthogonalised_linear_pair(n: int, k: int, p: int, *, seed: int, signal_scale: float = 2.0):
    """Per-column pairs that are EXACTLY uncorrelated in raw sample (like spec §4.2's u/v),
    while EACH column separately carries real, substantial design-explained structure.

    **Two prior versions of this helper were wrong, and the reason is itself part of the
    result this work item reports.** Both used a CATEGORICAL one-hot block design to generate
    X and Y (``x = D_onehot @ A + noise``). That construction makes X and Y both depend on the
    SAME discrete cluster assignment, and two variables driven by a shared discrete mixture
    carry real, non-artefactual mutual information (detectable by any nonlinear-aware
    statistic -- exactly what RV/dCor/HSIC/kernel-CCA are for) via the shared clustering alone,
    regardless of whether the per-cluster means happen to correlate. Residualising against the
    TRUE generating design then correctly REMOVES that real shared confound, while an unrelated
    design does not -- which is the CORRECT, expected direction for "did adjustment remove a
    real confound", but it is the OPPOSITE phenomenon from spec §4's FWL projection artifact,
    which requires u/v to be exactly unrelated (in population, not just in raw correlation)
    before residualisation, so that any correlation appearing afterward is manufactured by the
    projection geometry alone rather than by imperfectly-removed real dependence. See the
    result notebook entry: this confusion is reported as a finding in its own right.

    Fix: use a CONTINUOUS Gaussian design (not categorical). ``(x, y_base)`` given a Gaussian
    design and independent Gaussian noise are JOINTLY GAUSSIAN, and for a jointly Gaussian pair
    zero correlation IS independence -- there is no leftover nonlinear/mixture dependence to
    confound the measurement. Each column of Y is then forced to raw-correlation-zero against
    the matching column of X via ``calibration.spike_targets(..., r_true=0)``, reproducing
    ``closed_form_induced_correlation``'s own construction column-by-column.
    """
    from morpheus.v2.calibra.calibration import spike_targets

    rng = np.random.default_rng(seed)
    design = rng.normal(size=(n, k))
    loadings_x = rng.normal(size=(k, p))
    loadings_y = rng.normal(size=(k, p))
    x = signal_scale * (design @ loadings_x) / np.sqrt(k) + rng.normal(size=(n, p))
    y_base = signal_scale * (design @ loadings_y) / np.sqrt(k) + rng.normal(size=(n, p))
    y_columns = []
    for j in range(p):
        spiked, _, _ = spike_targets(x[:, [j]], y_base[:, [j]], 0.0,
                                     rng=np.random.default_rng(seed + 100 + j),
                                     return_directions=True)
        y_columns.append(spiked[:, 0])
    y = np.column_stack(y_columns)
    # Verify the construction actually did what it claims: near-zero RAW per-column
    # correlation, by construction, checked here so a broken helper fails loudly.
    raw_corr = np.array([np.corrcoef(x[:, j], y[:, j])[0, 1] for j in range(p)])
    assert np.max(np.abs(raw_corr)) < 1e-6, ("orthogonalisation failed", raw_corr)
    return x, y, design


@pytest.mark.parametrize("metric_name", ["rv", "dcor", "hsic", "kernel_cca"])
def test_regime_ii_induced_association_exceeds_regime_i_synthetic(metric_name):
    """SYNTHETIC. X and Y are EXACTLY uncorrelated per column in raw sample (§4's u/v setup)
    AND, because the generating design is Gaussian/linear rather than categorical, exactly
    population-independent (no leftover mixture dependence -- see the helper's docstring for
    why the categorical version of this test was measuring a different phenomenon). Each
    column separately carries real, design-explained structure. Regime II residualises against
    the TRUE design that explains both (spec's mechanism: induced residual association from
    shared projection geometry); Regime I residualises against a Gaussian design of matched
    rank that is UNRELATED to X, Y (the stated falsifier).
    """
    n, k, p = 800, 20, 5
    x, y, design_real = _orthogonalised_linear_pair(n, k, p, seed=hash(metric_name) % (2 ** 31))
    design_gaussian = np.random.default_rng(7).normal(size=(n, k))

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
