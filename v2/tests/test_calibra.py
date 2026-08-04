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


def test_cross_fitted_residuals_handle_a_single_column():
    """One column must residualise, and must agree with the same column inside a block.

    sklearn ravels a single-column Ridge target, so ``model.predict`` returns ``(n,)``
    where the matrix slice is ``(n, 1)``; the subtraction then broadcast to ``(n, n)``
    and raised. Per-AXIS residualisation is one column at a time, and it is the unit
    the P4 certification rule works in, so the one-column path is not a corner case --
    it is the case. The second assertion pins that the fix changes nothing for k >= 2:
    column 0 of a two-column call must equal the one-column call exactly.
    """
    rng = np.random.default_rng(11)
    site = rng.integers(0, 4, size=200)
    onehot = np.eye(4)[site]
    first = rng.normal(size=200) + 3.0 * onehot @ rng.normal(size=4)
    second = rng.normal(size=200) + 3.0 * onehot @ rng.normal(size=4)
    design = confound_design(pd.DataFrame({"site": site.astype(str)}), ["site"])

    single = cross_fitted_residuals(first[:, None], design)
    assert single.shape == (200, 1)
    assert np.isfinite(single).all()

    block = cross_fitted_residuals(np.column_stack([first, second]), design)
    np.testing.assert_allclose(single[:, 0], block[:, 0], rtol=0, atol=1e-12)

    # And it actually removes the confound rather than merely running.
    before = top_canonical_correlation(onehot, first[:, None], n_components=1)
    after = top_canonical_correlation(onehot, single, n_components=1)
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
                "recovered_median", "n_patients", "summary_schema_version",
                "observed_above_floor_status", "channel_statistic"):
        assert key in summary
    assert summary["summary_schema_version"] == 2
    # SCHEMA 2 CONTRACT: with no channel statistic supplied there is nothing to
    # grade, and the flag says so instead of inventing a verdict from a random
    # direction. `None` is deliberate -- a silent False would be the schema-1 bug.
    assert summary["observed_above_floor"] is None
    assert summary["observed_above_floor_status"] == "ungraded_no_channel_statistic"
    result.channel_statistic = 0.9
    result.channel_statistic_name = "test"
    graded = result.summary()
    assert isinstance(graded["observed_above_floor"], bool)
    assert graded["observed_above_floor_status"] == "graded"


def test_the_flag_is_a_magnitude_and_cannot_be_flipped_by_the_sign_of_the_channel():
    """B1 of PREDECLARED_observed_above_floor_20260804T1843Z, as a regression test.

    A detection floor is a threshold on the SIZE of an association. At schema 1 the
    flag compared a signed statistic against a positive floor, so two datasets
    differing only in the sign of the planted association returned opposite
    verdicts for identical evidence. Nothing about the fixed flag may depend on
    the sign."""
    from morpheus.v2.calibra.calibration import channel_clears_floor

    for magnitude in (0.05, 0.2, 0.6):
        positive = channel_clears_floor(magnitude, 0.1)
        negative = channel_clears_floor(-magnitude, 0.1)
        assert positive == negative, (magnitude, positive, negative)
    assert channel_clears_floor(0.6, 0.1)[0] is True
    assert channel_clears_floor(-0.6, 0.1)[0] is True
    assert channel_clears_floor(0.05, 0.1)[0] is False
    # Ungradable is None, never False: "we could not grade this" and "it is below
    # the floor" are opposite statements.
    assert channel_clears_floor(float("nan"), 0.1) == (None, "ungraded_no_channel_statistic")
    assert channel_clears_floor(0.6, float("nan")) == (None, "ungraded_floor_did_not_resolve")


def test_floors_from_recovery_reproduces_the_curve_s_own_floors():
    """The floor rule is lifted out so a different readout can use the SAME rule.
    If the extracted function ever drifts from the instrument, every
    direction-matched / like-for-like floor becomes a re-implementation."""
    from morpheus.v2.calibra.calibration import floors_from_recovery

    x, y = _paired(n=260, shared=1.5, seed=21)
    result = spike_recovery_curve(x, y, np.zeros((len(x), 0)), levels=(0.0, 0.05, 0.2, 0.4),
                                  n_draws=8, n_components=6, seed=21)
    floors = floors_from_recovery(result.levels, result.recovered,
                                  recovery_fraction=result.recovery_fraction)
    summary = result.summary()
    for key in ("detection_floor", "transmission_floor", "null_reference_p90"):
        assert np.allclose(floors[key], summary[key], equal_nan=True), key
    assert np.allclose(floors["unpaired_hit_rate"], summary["unpaired_hit_rate"], equal_nan=True)
    assert np.allclose(floors["paired_hit_rate"], summary["paired_hit_rate"], equal_nan=True)


def test_spike_targets_accepts_a_supplied_image_direction():
    """Without this the floor can only ever be measured on a direction pair that is
    random on the image side, while the channel it grades is fitted on both."""
    rng = np.random.default_rng(5)
    x, y = _paired(n=600, seed=5)
    u = np.zeros(x.shape[1]); u[3] = 1.0
    v = np.zeros(y.shape[1]); v[2] = 1.0
    for r_true in (0.0, 0.3, 0.6):
        spiked, got_u, got_v = spike_targets(x, y, r_true, rng=np.random.default_rng(0),
                                             image_direction=u, molecular_direction=v,
                                             return_directions=True)
        assert np.allclose(got_u, u) and np.allclose(got_v, v)
        observed = abs(np.corrcoef(x @ u, spiked @ v)[0, 1])
        assert abs(observed - r_true) < 0.05, (r_true, observed)
    # A supplied direction consumes no draw, exactly as molecular_direction already did.
    both = spike_targets(x, y, 0.2, rng=np.random.default_rng(0), image_direction=u,
                         molecular_direction=v)
    again = spike_targets(x, y, 0.2, rng=np.random.default_rng(77), image_direction=u,
                          molecular_direction=v)
    assert np.array_equal(both, again)
    for bad in (np.zeros(x.shape[1]),):
        try:
            spike_targets(x, y, 0.2, rng=np.random.default_rng(0), image_direction=bad)
        except ValueError:
            continue
        raise AssertionError("a zero image_direction must raise")


def test_direction_pairs_pin_the_planted_axis_across_draws():
    """A direction-matched floor is a floor measured on a SUPPLIED pair, one per draw."""
    from morpheus.v2.calibra.spectral import heldout_cca_directions

    x, y = _paired(n=400, shared=2.0, seed=31)
    order = np.random.default_rng(0).permutation(len(x))
    pairs = [heldout_cca_directions(x, y, order[:200], n_components=6),
             heldout_cca_directions(x, y, order[200:], n_components=6)]
    assert all(len(u) == x.shape[1] and len(v) == y.shape[1] for u, v in pairs)
    result = spike_recovery_curve(x, y, np.zeros((len(x), 0)), levels=(0.0, 0.1, 0.4),
                                  n_draws=4, n_components=6, seed=31, direction_pairs=pairs)
    assert result.summary()["direction_source"] == "fitted_pairs_supplied"
    # Pairs cycle, so draws 0 and 2 share a pair and must be bit-identical.
    assert np.array_equal(result.recovered[:, 0], result.recovered[:, 2])
    assert not np.array_equal(result.recovered[:, 0], result.recovered[:, 1])
    try:
        spike_recovery_curve(x, y, np.zeros((len(x), 0)), levels=(0.0, 0.1), n_draws=2,
                             direction_pairs=pairs, molecular_directions=np.eye(y.shape[1])[:1])
    except ValueError:
        pass
    else:
        raise AssertionError("direction_pairs and molecular_directions must be exclusive")


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


def _confounded(n, strength, seed, n_sites=8):
    rng = np.random.default_rng(seed)
    site = rng.integers(0, n_sites, size=n)
    onehot = np.eye(n_sites)[site]
    x = rng.normal(size=(n, 20)) + strength * onehot @ rng.normal(size=(n_sites, 20))
    y = rng.normal(size=(n, 18)) + strength * onehot @ rng.normal(size=(n_sites, 18))
    design = confound_design(pd.DataFrame({"site": site.astype(str)}), ["site"])
    return x, y, design


def test_floor_survives_a_confound_induced_baseline():
    """On real data, residualising through a wide confound design induces a non-zero
    level-0 baseline even though the spike is built orthogonal to the image score.
    The paired floor must still resolve; an unpaired threshold against the level-0
    p90 returned NaN here, which is what it did on every real state."""
    x, y, design = _confounded(n=800, strength=1.0, seed=31)
    result = spike_recovery_curve(x, y, design, levels=(0.0, 0.05, 0.1, 0.2, 0.4),
                                  n_draws=12, n_components=8, seed=31)
    summary = result.summary()
    assert summary["confound_induced_baseline"] > 0.01, "fixture induces no baseline; test is vacuous"
    assert np.isfinite(result.detection_floor), (
        f"paired floor failed at induced baseline {summary['confound_induced_baseline']:.3f}")
    assert summary["baseline_is_null_like"] is True, summary["confound_induced_baseline"]


def test_gate_warns_when_the_induced_baseline_rivals_the_real_signal():
    """The complement of the test above: the gate must not be decorative. Under an
    aggressive design (few patients, many strong site levels) the induced baseline
    grows until it is a substantial fraction of the measurable channel, and at that
    point a floor is not trustworthy and the run must be flagged."""
    x, y, design = _confounded(n=500, strength=3.0, seed=32)
    summary = spike_recovery_curve(x, y, design, levels=(0.0, 0.1, 0.4),
                                   n_draws=12, n_components=8, seed=32).summary()
    assert summary["confound_induced_baseline"] > summary["baseline_gate_threshold"]
    assert summary["baseline_is_null_like"] is False, "gate slept through a pathological baseline"


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
