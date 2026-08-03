"""Tests for the induced-correlation sweep (Track 2 / phase-gate C).

The headline claim these guard is that residualising two ORTHOGONAL signals
through a shared confound design manufactures correlation between them. That is
a methodological claim about every confound-adjusted cross-modal analysis, so
the tests are written to attack it: they check the derived identity against the
pipeline, they check that a design of matched rank with no relationship to the
patients induces nothing (the ledger's stated falsifier), and they check the
scaling behaviour that separates the derived law from the plan's original guess.

All fixtures are small and synthetic; the whole module runs in seconds on CPU.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from morpheus.v2.calibra.calibration import spike_targets
from morpheus.v2.calibra.induced_correlation_sweep import (
    DesignSpec, build_design, closed_form_induced_correlation, design_participation_ratio,
    fit_scaling_law, planted_directions, predict_induced_correlation, shared_design_geometry,
    stratified_subsample, sweep_cell)
from morpheus.v2.calibra.residualise import confound_design


def _confounded(n, strength, seed, n_sites=8, p=20, q=18):
    """Two modalities that share a site effect — the regime real cohorts live in."""
    rng = np.random.default_rng(seed)
    site = rng.integers(0, n_sites, size=n)
    onehot = np.eye(n_sites)[site]
    x = rng.normal(size=(n, p)) + strength * onehot @ rng.normal(size=(n_sites, p))
    y = rng.normal(size=(n, q)) + strength * onehot @ rng.normal(size=(n_sites, q))
    design = confound_design(pd.DataFrame({"site": site.astype(str)}), ["site"])
    return x, y, design


def _induced(x, y, design, **kwargs):
    return sweep_cell(x, y, design, n_draws=kwargs.pop("n_draws", 12),
                      n_components=kwargs.pop("n_components", 8), **kwargs)


# --- the mechanism ---------------------------------------------------------

def test_closed_form_reproduces_the_pipeline_draw_for_draw():
    """THE mechanism test. ``r_induced = -corr(s_hat, a_hat) * kappa`` is an
    identity, not a fit: the induced correlation IS minus the correlation of the
    two design-explained parts. If this ever stops holding to ~1e-12, the written
    mechanism no longer describes what the instrument does and the claim must be
    withdrawn before anything else is debugged."""
    for strength, n_sites in ((1.0, 8), (3.0, 8), (2.0, 30)):
        x, y, design = _confounded(n=700, strength=strength, seed=31, n_sites=n_sites)
        row = _induced(x, y, design, seed=31)
        assert row["closed_form_max_abs_error"] < 1e-10, (strength, n_sites, row["closed_form_max_abs_error"])
        assert row["closed_form_pearson"] > 1 - 1e-9
        assert abs(row["induced_correlation_median"] - row["closed_form_induced_median"]) < 1e-10


def test_planted_directions_match_the_spike_the_pipeline_plants():
    """The closed form is only a valid comparator if it scores the SAME (u, v)."""
    x, y, _ = _confounded(n=300, strength=1.0, seed=5)
    pairs = planted_directions(x.shape[1], y.shape[1], n_draws=3, seed=17)
    rng = np.random.default_rng(17)
    draw_seeds = [int(rng.integers(1 << 31)) for _ in range(3)]
    for (u, v), draw_seed in zip(pairs, draw_seeds):
        _, u_pipeline, v_pipeline = spike_targets(x, y, 0.0, rng=np.random.default_rng(draw_seed),
                                                  return_directions=True)
        np.testing.assert_allclose(u, u_pipeline, rtol=0, atol=0)
        np.testing.assert_allclose(v, v_pipeline, rtol=0, atol=0)


# --- the falsifier: matched rank, no structure -----------------------------

def test_matched_rank_design_without_structure_induces_nothing():
    """The ledger's stated falsifier, run deliberately.

    A Gaussian design and a row-permuted real design have the SAME rank and the
    SAME n as the real design, and differ only in whether they explain anything
    about the patients. If induced correlation were an artefact of design width
    or of our residualiser, all three would agree. They must not."""
    x, y, design = _confounded(n=800, strength=2.0, seed=41, n_sites=20)
    rng = np.random.default_rng(41)
    real = _induced(x, y, design, seed=41)["induced_correlation_median"]
    gaussian = _induced(x, y, rng.normal(size=(len(x), design.shape[1])), seed=41)["induced_correlation_median"]
    permuted = _induced(x, y, design[rng.permutation(len(x))], seed=41)["induced_correlation_median"]
    assert real > 0.05, f"fixture induces nothing; the test would be vacuous ({real})"
    assert gaussian < real / 4, (real, gaussian)
    assert permuted < real / 4, (real, permuted)


def test_prediction_is_zero_when_the_design_explains_nothing():
    """Both declared predictors must return ~0 for a structureless design rather
    than NaN — a NaN would silently drop the falsifier arm out of every fit."""
    rng = np.random.default_rng(3)
    x, y = rng.normal(size=(500, 20)), rng.normal(size=(500, 18))
    design = rng.normal(size=(500, 40))
    closed = closed_form_induced_correlation(x, y, design, n_draws=8, seed=3)
    predicted = predict_induced_correlation(closed, design_participation_ratio(design))
    assert predicted["predicted_induced_correlation_p1"] < 0.05
    assert predicted["predicted_induced_correlation_p2"] < 0.05


# --- scaling ---------------------------------------------------------------

def test_induced_correlation_does_not_vanish_as_n_grows():
    """Separates the derived law from the plan's ``k/n`` guess.

    The plan expected the induced correlation to scale with design rank over n,
    i.e. to be a small-sample artefact that more patients remove. The derivation
    says it is a *bias*: R_s, R_a and k_eff converge, so it converges to a
    non-zero constant. Quadrupling n must therefore NOT quarter it."""
    small = _induced(*_confounded(n=400, strength=2.0, seed=7), seed=7)["induced_correlation_median"]
    large = _induced(*_confounded(n=1600, strength=2.0, seed=7), seed=7)["induced_correlation_median"]
    assert small > 0.02 and large > 0.02, (small, large)
    assert large > 0.4 * small, f"k/n scaling would give ~{small / 4:.4f}; measured {large:.4f}"


def test_induced_correlation_falls_as_the_design_gains_rank():
    """The counter-intuitive half of the law: at comparable explained variance a
    WIDER design induces a SMALLER per-direction correlation, because the two
    design-explained parts have more room to be orthogonal (1/sqrt(k_eff))."""
    narrow = _induced(*_confounded(n=2000, strength=2.0, seed=11, n_sites=4), seed=11)
    wide = _induced(*_confounded(n=2000, strength=2.0, seed=11, n_sites=64), seed=11)
    assert wide["k_eff"] > 4 * narrow["k_eff"], (narrow["k_eff"], wide["k_eff"])
    assert wide["induced_correlation_median"] < narrow["induced_correlation_median"], (
        narrow["induced_correlation_median"], wide["induced_correlation_median"])


def test_design_participation_ratio_tracks_a_balanced_one_hot_rank():
    """``k_eff`` must be the design's usable dimension, not its column count: a
    one-hot block is rank deficient and ridge shrinks weak levels away, so
    quoting the column count would misstate the law's denominator."""
    rng = np.random.default_rng(2)
    site = rng.integers(0, 10, size=2000)
    design = confound_design(pd.DataFrame({"site": site.astype(str)}), ["site"])
    geometry = design_participation_ratio(design)
    assert design.shape[1] == 11                      # 10 levels + dummy_na
    assert geometry["design_rank"] == 9               # collinear with the intercept
    assert 8.0 < geometry["k_eff"] < 10.0, geometry
    assert design_participation_ratio(np.zeros((50, 0)))["k_eff"] == 0.0


def test_shared_geometry_ignores_rank_that_neither_modality_loads_on():
    """Why the law's denominator is ``k_eff_shared`` and not the design's rank.

    Appending design columns that explain nothing about either modality raises
    the rank and the plain participation ratio, but must NOT change the
    dimension over which the two design-explained parts can be orthogonal. This
    is the reason the measured induced correlation is flat from k=33 to k=501 on
    the real cohort while the rank-based predictors decay as 1/sqrt(k)."""
    rng = np.random.default_rng(19)
    n, n_sites = 1200, 6
    site = rng.integers(0, n_sites, size=n)
    onehot = np.eye(n_sites)[site]
    x = rng.normal(size=(n, 15)) + 3.0 * onehot @ rng.normal(size=(n_sites, 15))
    y = rng.normal(size=(n, 12)) + 3.0 * onehot @ rng.normal(size=(n_sites, 12))
    narrow = confound_design(pd.DataFrame({"site": site.astype(str)}), ["site"])
    inert = rng.integers(0, 60, size=n)                 # a nuisance neither modality knows about
    wide = np.hstack([narrow, confound_design(pd.DataFrame({"z": inert.astype(str)}), ["z"])])

    narrow_geometry = shared_design_geometry(narrow, x, y)
    wide_geometry = shared_design_geometry(wide, x, y)
    assert design_participation_ratio(wide)["k_eff"] > 5 * design_participation_ratio(narrow)["k_eff"]
    assert wide_geometry["k_eff_shared"] < 2.0 * narrow_geometry["k_eff_shared"], (
        narrow_geometry, wide_geometry)
    # ...and the measured induced correlation follows k_eff_shared, not the rank.
    narrow_induced = _induced(x, y, narrow, seed=19)["induced_correlation_median"]
    wide_induced = _induced(x, y, wide, seed=19)["induced_correlation_median"]
    assert wide_induced > 0.4 * narrow_induced, (narrow_induced, wide_induced)


def test_shared_geometry_separates_alignment_from_explained_variance():
    """The reason (4) uses the CROSS term ``A'B`` and not the two variance
    profiles separately.

    Both cases below have the same design, the same rank, and near-identical R^2
    for each modality. They differ only in WHETHER the design explains the two
    modalities along the same directions. The induced correlation differs ~5x,
    and ``k_eff_shared`` must track that; a variance-only denominator cannot see
    it and would predict the same number for both."""
    rng = np.random.default_rng(23)
    n = 1500
    factor_a, factor_b = rng.integers(0, 10, size=n), rng.integers(0, 10, size=n)
    onehot_a, onehot_b = np.eye(10)[factor_a], np.eye(10)[factor_b]
    x = 3.0 * onehot_a @ rng.normal(size=(10, 12)) + rng.normal(size=(n, 12))
    y_shared = 3.0 * onehot_a @ rng.normal(size=(10, 12)) + rng.normal(size=(n, 12))
    y_disjoint = 3.0 * onehot_b @ rng.normal(size=(10, 12)) + rng.normal(size=(n, 12))
    design = confound_design(pd.DataFrame({"a": factor_a.astype(str), "b": factor_b.astype(str)}),
                             ["a", "b"])
    shared = shared_design_geometry(design, x, y_shared)
    disjoint = shared_design_geometry(design, x, y_disjoint)
    assert abs(shared["predicted_r2_y"] - disjoint["predicted_r2_y"]) < 0.05, (shared, disjoint)
    assert disjoint["k_eff_shared"] > 10 * shared["k_eff_shared"], (shared, disjoint)
    induced_shared = _induced(x, y_shared, design, seed=23)["induced_correlation_median"]
    induced_disjoint = _induced(x, y_disjoint, design, seed=23)["induced_correlation_median"]
    assert induced_disjoint < 0.3 * induced_shared, (induced_shared, induced_disjoint)


def test_fit_scaling_law_recovers_a_planted_exponent():
    rows = []
    for k_eff in (5.0, 20.0, 80.0, 320.0):
        for n in (500.0, 2000.0, 8000.0):
            rows.append({"k_eff": k_eff, "n_patients": n,
                         "induced_correlation_median": 0.3 * k_eff ** -0.5 * n ** 0.0})
    fit = fit_scaling_law(rows)
    assert abs(fit["exponent_k_eff"] + 0.5) < 1e-6, fit
    assert abs(fit["exponent_n"]) < 1e-6, fit
    assert fit["residual_rms_log"] < 1e-9


# --- T2.5: the estimator is not the cause ----------------------------------

def test_fold_count_does_not_move_the_induced_correlation():
    """Cross-fitting is a deliberate choice (in-sample residualisation removes
    more than the confound). If the induced correlation were an artefact of the
    fold count it would move with it."""
    x, y, design = _confounded(n=800, strength=2.0, seed=13, n_sites=16)
    values = [_induced(x, y, design, seed=13, n_splits=k)["induced_correlation_median"]
              for k in (2, 3, 5, 10, 20)]
    assert min(values) > 0.05, values
    assert max(values) / min(values) < 1.15, values


def test_shrinkage_moves_the_effect_only_through_how_much_it_residualises():
    """The honest reading of the alpha knob.

    Ridge alpha does eventually change the number — but only by making the design
    stop residualising at all. Equation (1) says the induced correlation is a
    function of the design's R^2 for the two scores, so a shrinkage so heavy that
    R^2 collapses MUST take the induced correlation with it. That is the
    mechanism working, not an estimator artefact. The estimator-artefact
    hypothesis would instead show the effect moving while R^2 stayed put."""
    x, y, design = _confounded(n=800, strength=2.0, seed=13, n_sites=16)
    rows = {alpha: _induced(x, y, design, seed=13, alpha=alpha)
            for alpha in (0.01, 0.1, 1.0, 10.0, 1000.0)}
    mild = [rows[a]["induced_correlation_median"] for a in (0.01, 0.1, 1.0, 10.0)]
    assert max(mild) / min(mild) < 1.3, mild            # four decades of alpha, ~flat
    crushed = rows[1000.0]
    assert crushed["design_r2_x_median"] < 0.2 * rows[1.0]["design_r2_x_median"]
    assert crushed["induced_correlation_median"] < 0.2 * rows[1.0]["induced_correlation_median"]
    for row in rows.values():                            # the identity never breaks
        assert row["closed_form_max_abs_error"] < 1e-10, row["residualiser_alpha"]


# --- cohort control --------------------------------------------------------

def test_stratified_subsample_holds_composition_and_hits_the_target_size():
    """An unstratified subsample would move cohort composition together with n,
    and composition is exactly what the design residualises — the two effects
    would then be inseparable in the n-sweep."""
    rng = np.random.default_rng(0)
    strata = np.repeat(np.array(["A", "B", "C", "D"]), [1000, 500, 200, 50])
    for n_target in (100, 875, 1500):
        idx = stratified_subsample(strata, n_target, seed=1)
        assert len(idx) == n_target
        assert len(set(idx.tolist())) == n_target
        for level in ("A", "B", "C", "D"):
            share_full = float(np.mean(strata == level))
            share_sub = float(np.mean(strata[idx] == level))
            assert abs(share_sub - share_full) < 0.02, (n_target, level, share_full, share_sub)
    assert len(stratified_subsample(strata, 10_000, seed=1)) == len(strata)


def test_design_spec_modes_and_missing_covariates_are_explicit():
    ids = np.array([f"TCGA-{s:02d}-{i:04d}" for i in range(200) for s in (1, 2)][:200])
    cancers = np.array(["BRCA", "LUAD"] * 100)
    real, meta = build_design(DesignSpec("d", ("cancer", "tss"), min_site_count=1), ids, cancers, seed=0)
    assert meta["n_confound_columns"] == real.shape[1] > 0
    permuted, _ = build_design(DesignSpec("p", ("cancer", "tss"), min_site_count=1, mode="permuted"),
                               ids, cancers, seed=0)
    assert permuted.shape == real.shape
    np.testing.assert_allclose(np.sort(permuted.sum(axis=0)), np.sort(real.sum(axis=0)))
    gaussian, meta = build_design(DesignSpec("g", mode="gaussian", width=7), ids, cancers, seed=0)
    assert gaussian.shape == (200, 7) and meta["n_confound_columns"] == 7
    for bad in (dict(mode="nonsense"), dict(mode="gaussian")):
        try:
            DesignSpec("bad", ("cancer",), **bad)
        except ValueError:
            continue
        raise AssertionError(f"DesignSpec accepted {bad}")
