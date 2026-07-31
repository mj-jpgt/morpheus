"""Tests for CALIBRA. The spike-recovery tests ARE the instrument's self-test:
if these fail, no CALIBRA number means anything."""
from __future__ import annotations

import numpy as np
import pandas as pd

from morpheus.v2.calibra import (cca_spectrum, confound_design, cross_fitted_residuals,
                                 effective_rank, spike_recovery_curve, spike_targets,
                                 top_canonical_correlation)


def _paired(n=300, p=24, q=16, shared=0.0, seed=0):
    """Paired modalities with a controllable shared component."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, p))
    y = rng.normal(size=(n, q))
    if shared:
        z = rng.normal(size=n)
        x[:, 0] += shared * z
        y[:, 0] += shared * z
    return x, y


# --- spectral -------------------------------------------------------------

def test_effective_rank_matches_known_rank():
    rng = np.random.default_rng(0)
    low = rng.normal(size=(200, 3)) @ rng.normal(size=(3, 40))
    assert effective_rank(low) < 4.0
    assert effective_rank(rng.normal(size=(200, 40))) > 20.0
    assert effective_rank(np.ones((50, 8))) == 0.0


def test_cca_spectrum_is_bounded_and_ordered():
    x, y = _paired(shared=1.5)
    spectrum = cca_spectrum(x, y, n_components=8)
    assert spectrum.size == 8
    assert np.all((spectrum >= 0) & (spectrum <= 1.0 + 1e-9))
    assert np.all(np.diff(spectrum) <= 1e-9)  # descending


# --- spike construction ---------------------------------------------------

def test_spike_targets_installs_the_requested_correlation():
    """The spike must be exactly the strength we asked for, pre-residualisation."""
    rng = np.random.default_rng(3)
    x, y = _paired(n=800, seed=3)
    for r_true in (0.0, 0.05, 0.2, 0.5):
        gen = np.random.default_rng(11)
        u = gen.normal(size=x.shape[1]); u /= np.linalg.norm(u)
        v = gen.normal(size=y.shape[1]); v /= np.linalg.norm(v)
        spiked = spike_targets(x, y, r_true, rng=np.random.default_rng(11))
        s = x @ u; a = spiked @ v
        observed = abs(np.corrcoef(s, a)[0, 1])
        assert abs(observed - r_true) < 0.05, f"asked {r_true}, installed {observed}"


def test_spike_targets_rejects_invalid_strength():
    x, y = _paired()
    for bad in (-0.1, 1.0, 1.5):
        try:
            spike_targets(x, y, bad, rng=np.random.default_rng(0))
        except ValueError:
            continue
        raise AssertionError(f"r_true={bad} should have raised")


# --- residualisation ------------------------------------------------------

def test_cross_fitted_residuals_remove_the_confound():
    rng = np.random.default_rng(5)
    site = rng.integers(0, 4, size=400)
    x = rng.normal(size=(400, 10)) + 5.0 * np.eye(4)[site] @ rng.normal(size=(4, 10))
    frame = pd.DataFrame({"site": site.astype(str)})
    design = confound_design(frame, ["site"])
    residual = cross_fitted_residuals(x, design)
    # Site should no longer be linearly recoverable from the residual.
    before = top_canonical_correlation(np.eye(4)[site], x, n_components=3)
    after = top_canonical_correlation(np.eye(4)[site], residual, n_components=3)
    assert after < before - 0.2, (before, after)


def test_confound_design_handles_missing_numerics():
    frame = pd.DataFrame({"purity": [0.5, np.nan, 0.8, 0.2], "site": ["A", "B", "A", "B"]})
    design = confound_design(frame, ["purity", "site"])
    assert np.isfinite(design).all()
    assert design.shape[0] == 4
    # a missingness indicator column must exist for the NaN
    assert design.shape[1] >= 4


# --- THE INSTRUMENT SELF-TEST --------------------------------------------

def test_recovery_curve_is_monotone_and_floor_is_sane():
    """Stronger spikes must be recovered more strongly, and the floor must sit
    above the level-0 null but below the largest spike."""
    x, y = _paired(n=400, seed=7)
    design = np.zeros((len(x), 0))
    result = spike_recovery_curve(x, y, design, levels=(0.0, 0.05, 0.2, 0.5),
                                  n_draws=6, n_components=8, seed=7)
    median = np.nanmedian(result.recovered, axis=1)
    assert median[-1] > median[0], f"strongest spike not recovered above null: {median}"
    assert result.attenuation_slope > 0
    assert np.isfinite(result.detection_floor)
    assert 0.0 < result.detection_floor <= 0.5


def test_residualising_the_signal_away_raises_the_floor():
    """The instrument must DETECT over-residualisation: if we adjust away the very
    axis the spike lives on, the detection floor must get worse (or vanish).
    This is the check that makes a CALIBRA null interpretable."""
    rng = np.random.default_rng(9)
    n = 400
    z = rng.normal(size=n)                      # the axis the spike will use
    x = rng.normal(size=(n, 12)); x[:, 0] += 3.0 * z
    y = rng.normal(size=(n, 10)); y[:, 0] += 3.0 * z
    benign = np.zeros((n, 0))
    destructive = z[:, None]                    # residualise out the shared axis
    lenient = spike_recovery_curve(x, y, benign, levels=(0.0, 0.05, 0.2),
                                   n_draws=6, n_components=6, seed=9)
    harsh = spike_recovery_curve(x, y, destructive, levels=(0.0, 0.05, 0.2),
                                 n_draws=6, n_components=6, seed=9)
    lenient_floor = lenient.detection_floor if np.isfinite(lenient.detection_floor) else 1.0
    harsh_floor = harsh.detection_floor if np.isfinite(harsh.detection_floor) else 1.0
    assert harsh_floor >= lenient_floor, (lenient_floor, harsh_floor)


def test_summary_reports_observed_against_floor():
    x, y = _paired(n=300, shared=2.0, seed=13)
    result = spike_recovery_curve(x, y, np.zeros((len(x), 0)), levels=(0.0, 0.1),
                                  n_draws=4, n_components=6, seed=13)
    summary = result.summary()
    for key in ("detection_floor", "observed", "observed_above_floor", "attenuation_slope",
                "recovered_median", "n_patients"):
        assert key in summary
    assert isinstance(summary["observed_above_floor"], bool)


# --- held-out CCA: the fix for in-sample capacity inflation -----------------

def test_heldout_cca_is_lower_than_insample_on_noise():
    """In-sample CCA is a multivariate maximum and is upward-biased at finite n.
    On pure noise the held-out estimate must not inherit that bias."""
    from morpheus.v2.calibra.spectral import heldout_top_cca
    rng = np.random.default_rng(0)
    x = rng.normal(size=(300, 40))
    y = rng.normal(size=(300, 40))
    insample = top_canonical_correlation(x, y, n_components=16)
    held = heldout_top_cca(x, y, n_components=16, seed=0)
    assert insample > 0.3, insample                      # capacity inflation is real
    assert held < insample - 0.15, (insample, held)


# --- REGRESSION: the targeted-readout defect -------------------------------
# The first implementation scored recovery with top-CCA, a maximum over all
# components, while the spike lives on one known direction. On real data the
# ambient top-CCA sits near 0.97, so every detection floor came back NaN. These
# three tests reproduce that failure mode on synthetic data and would have caught
# it before it reached a results file.

def _ambient(n=500, p=24, q=20, n_shared=6, strength=4.0, seed=0):
    """Paired modalities with STRONG pre-existing multivariate structure, i.e. the
    regime real data lives in (adjusted top-CCA ~0.9)."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, p))
    y = rng.normal(size=(n, q))
    for k in range(n_shared):                       # many shared axes -> high top-CCA
        z = rng.normal(size=n)
        x[:, k] += strength * z
        y[:, k] += strength * z
    return x, y


def test_level_zero_baseline_is_null_like_despite_strong_ambient_signal():
    """THE regression test. With a heavily structured pair, a max-based readout
    reports ~0.97 at r_true=0; a targeted readout must report ~0."""
    x, y = _ambient(seed=21)
    assert top_canonical_correlation(x, y, n_components=8) > 0.85, "fixture is not ambient-heavy"
    result = spike_recovery_curve(x, y, np.zeros((len(x), 0)), levels=(0.0, 0.05, 0.2),
                                  n_draws=8, n_components=8, seed=21)
    summary = result.summary()
    baseline = summary["baseline_recovered_median"]
    assert baseline < 0.15, f"level-0 readout {baseline:.3f} is picking up ambient structure"
    assert summary["baseline_is_null_like"] is True


def test_floor_is_finite_when_ambient_structure_is_strong():
    """The defect's visible symptom was NaN floors on every real state."""
    x, y = _ambient(seed=22)
    result = spike_recovery_curve(x, y, np.zeros((len(x), 0)), levels=(0.0, 0.02, 0.05, 0.1, 0.2),
                                  n_draws=8, n_components=8, seed=22)
    assert np.isfinite(result.detection_floor), "floor is NaN under ambient structure (the old bug)"
    assert 0.0 < result.detection_floor <= 0.2


def test_recovery_tracks_the_injected_strength():
    """Recovered correlation must rise with r_true and stay in its neighbourhood.
    The old readout was flat at ~0.97 and even went DOWN at the largest spike."""
    x, y = _ambient(seed=23)
    levels = (0.0, 0.1, 0.3, 0.6)
    result = spike_recovery_curve(x, y, np.zeros((len(x), 0)), levels=levels,
                                  n_draws=8, n_components=8, seed=23)
    median = np.nanmedian(result.recovered, axis=1)
    assert np.all(np.diff(median) > 0), f"not monotone in r_true: {median}"
    assert abs(median[-1] - 0.6) < 0.15, f"r_true=0.6 recovered as {median[-1]:.3f}"
    assert 0.5 < result.attenuation_slope <= 1.2, result.attenuation_slope


def test_n_jobs_does_not_change_the_numbers():
    """Parallelism must be a wall-clock decision only."""
    x, y = _ambient(n=240, seed=24)
    kwargs = dict(levels=(0.0, 0.1, 0.3), n_draws=4, n_components=6, seed=24)
    serial = spike_recovery_curve(x, y, np.zeros((len(x), 0)), n_jobs=1, **kwargs)
    parallel = spike_recovery_curve(x, y, np.zeros((len(x), 0)), n_jobs=2, **kwargs)
    np.testing.assert_allclose(serial.recovered, parallel.recovered, rtol=1e-10, atol=1e-12)


def test_heldout_cca_recovers_a_planted_signal():
    from morpheus.v2.calibra.spectral import heldout_top_cca
    rng = np.random.default_rng(1)
    z = rng.normal(size=600)
    x = rng.normal(size=(600, 20)); x[:, 0] += 4 * z
    y = rng.normal(size=(600, 20)); y[:, 0] += 4 * z
    assert heldout_top_cca(x, y, n_components=8, seed=1) > 0.6
