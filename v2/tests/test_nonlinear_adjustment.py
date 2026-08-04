"""Tests for `v2/calibra/nonlinear_adjustment.py`.

Two of these are contracts rather than checks, and they are the reason any number this
module produces is comparable to P1 4.4's:

* the generic adjuster at ``model="ridge"`` returns **exactly** what
  ``residualise.cross_fitted_residuals`` returns, and
* ``channel_under_adjustment`` under that adjuster reproduces
  ``calibration.permutation_null``'s dictionary **exactly**.

The rest establish that the arms do what their docstrings claim, in particular that a
first-moment adjustment of *any* flexibility leaves a variance-only confound completely
legible to a k-NN while the location-scale arm removes it.  That construction is the
whole question this module exists to answer, so it is a test and not an argument.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.calibration import permutation_null
from morpheus.v2.calibra.nonlinear_adjustment import (ADJUSTER_MODELS, KernelRidgeAdjuster,
                                                      adjuster_agreement, cell_codes, cell_design,
                                                      channel_under_adjustment,
                                                      cross_fitted_location_scale, cross_fitted_r2,
                                                      cross_fitting_offset_energy, forest_residuals,
                                                      in_sample_residuals, kernel_ridge_residuals,
                                                      labels_only_ceiling, make_adjuster,
                                                      row_shuffled, saturated_cell_residuals)
from morpheus.v2.calibra.nonlinear_confound_probe import knn_balanced_accuracy_oof
from morpheus.v2.calibra.residualise import confound_design, cross_fitted_residuals


def _toy(seed: int = 0, n: int = 240, p: int = 6, q: int = 5, n_cancer: int = 4, n_site: int = 8):
    """A small cohort with site nested inside cancer, as TCGA's is."""
    import pandas as pd

    rng = np.random.default_rng(seed)
    site = rng.integers(0, n_site, n)
    cancer = site % n_cancer                      # nested: site determines cancer
    frame = pd.DataFrame({"cancer": [f"C{c}" for c in cancer], "tss": [f"S{s}" for s in site]})
    design = confound_design(frame, ["cancer", "tss"])
    codes, _ = cell_codes(frame["cancer"].to_numpy(), frame["tss"].to_numpy())
    x = rng.normal(size=(n, p)) + np.eye(n_site)[site] @ rng.normal(size=(n_site, p))
    y = rng.normal(size=(n, q)) + x[:, :1] * 0.8 + np.eye(n_site)[site] @ rng.normal(size=(n_site, q))
    return x, y, design, codes, frame["cancer"].to_numpy()


# --- the two contracts -------------------------------------------------------------------

def test_ridge_adjuster_is_the_incumbent_bit_for_bit():
    x, _, design, _, _ = _toy()
    adjust = make_adjuster("ridge", design=design, n_splits=5, seed=42)
    assert np.array_equal(adjust(x), cross_fitted_residuals(x, design, n_splits=5, alpha=1.0,
                                                            seed=42))


def test_channel_under_incumbent_reproduces_permutation_null_exactly():
    x, y, design, _, cancers = _toy()
    adjust = make_adjuster("ridge", design=design, n_splits=5, seed=42)
    mine = channel_under_adjustment(x, y, adjust, strata=cancers, n_permutations=25,
                                    n_components=4, seed=42)
    theirs = permutation_null(x, y, design, strata=cancers, n_permutations=25, n_components=4,
                              seed=42)
    for key in ("observed_top_cca", "null_median", "null_p95", "null_max",
                "excess_over_null_median", "permutation_p", "n_permutations"):
        assert mine[key] == theirs[key], key


# --- cells -------------------------------------------------------------------------------

def test_cell_codes_and_design_are_the_saturated_structure():
    a = np.array(["x", "x", "y", "y"])
    b = np.array(["1", "2", "1", "1"])
    codes, levels = cell_codes(a, b)
    assert len(levels) == 3 and codes[2] == codes[3] and codes[0] != codes[1]
    design = cell_design(codes)
    assert design.shape == (4, 3)
    assert np.array_equal(design.sum(axis=1), np.ones(4))


def test_cell_codes_rejects_ragged_input():
    with pytest.raises(ValueError):
        cell_codes(np.array(["a", "b"]), np.array(["a"]))


def test_saturated_arm_removes_a_pure_cell_mean_signal():
    """A block that IS its cell means residualises to noise; the additive arm cannot."""
    rng = np.random.default_rng(3)
    n, n_cancer, n_site = 600, 5, 10
    site = rng.integers(0, n_site, n)
    cancer = rng.integers(0, n_cancer, n)         # crossed, so the interaction is estimable
    codes, levels = cell_codes([f"C{c}" for c in cancer], [f"S{s}" for s in site])
    interaction = rng.normal(size=(len(levels), 4)) * 3.0
    x = interaction[codes] + rng.normal(size=(n, 4)) * 0.1
    import pandas as pd
    design = confound_design(pd.DataFrame({"cancer": [f"C{c}" for c in cancer],
                                           "tss": [f"S{s}" for s in site]}), ["cancer", "tss"])
    additive = cross_fitted_residuals(x, design, seed=42)
    saturated = saturated_cell_residuals(x, codes, seed=42)
    assert np.linalg.norm(saturated) < 0.5 * np.linalg.norm(additive)


# --- kernel ridge ------------------------------------------------------------------------

def test_kernel_ridge_is_deterministic_centred_and_shaped():
    x, _, design, _, _ = _toy()
    first = kernel_ridge_residuals(x, design, alpha=1.0, gamma=0.5, seed=42)
    second = kernel_ridge_residuals(x, design, alpha=1.0, gamma=0.5, seed=42)
    assert first.shape == x.shape
    assert np.array_equal(first, second)
    assert np.allclose(first.mean(axis=0), 0.0, atol=1e-9)


def test_kernel_ridge_factorises_once_and_reuses_it():
    """The adjuster is reusable across matrices -- the property the null depends on."""
    x, y, design, _, _ = _toy()
    adjuster = KernelRidgeAdjuster(design, alpha=1.0, gamma=0.5, seed=42)
    assert np.array_equal(adjuster(x), kernel_ridge_residuals(x, design, alpha=1.0, gamma=0.5,
                                                              seed=42))
    assert np.array_equal(adjuster(y), kernel_ridge_residuals(y, design, alpha=1.0, gamma=0.5,
                                                              seed=42))


def test_kernel_ridge_rejects_a_block_of_the_wrong_height():
    x, _, design, _, _ = _toy()
    adjuster = KernelRidgeAdjuster(design, alpha=1.0, gamma=0.5, seed=42)
    with pytest.raises(ValueError):
        adjuster(x[:10])


def test_kernel_ridge_removes_a_confound_mean_at_least_as_well_as_the_additive_ridge():
    x, _, design, codes, _ = _toy()
    additive = cross_fitted_residuals(x, design, seed=42)
    kernel = kernel_ridge_residuals(x, design, alpha=0.1, gamma=0.5, seed=42)
    assert np.linalg.norm(kernel) <= np.linalg.norm(additive) * 1.05


# --- forest ------------------------------------------------------------------------------

def test_forest_arm_runs_multi_output_and_differs_from_the_incumbent():
    x, _, design, _, _ = _toy()
    forest = forest_residuals(x, design, n_estimators=20, seed=42)
    assert forest.shape == x.shape
    assert np.allclose(forest.mean(axis=0), 0.0, atol=1e-9)
    assert not np.allclose(forest, cross_fitted_residuals(x, design, seed=42))


# --- the question this module exists for --------------------------------------------------

def _variance_only_confound(seed: int = 11, n: int = 900, p: int = 8, n_classes: int = 3):
    """Classes with IDENTICAL means and different variances.

    No conditional-mean adjustment of any flexibility can touch this, because there is no
    conditional mean to remove.  It is the situation the whole question is about.
    """
    rng = np.random.default_rng(seed)
    labels = rng.integers(0, n_classes, n)
    scale = np.array([1.0, 2.5, 5.0])[labels][:, None]
    x = rng.normal(size=(n, p)) * scale
    x = x - x.mean(axis=0, keepdims=True)
    return x, labels


def test_first_moment_adjustment_leaves_a_variance_confound_fully_legible():
    x, labels = _variance_only_confound()
    design = np.eye(3)[labels]
    chance = 1.0 / 3.0
    for name, adjusted in (("ridge", cross_fitted_residuals(x, design, seed=42)),
                           ("saturated", saturated_cell_residuals(x, labels, seed=42)),
                           ("kernel_ridge", kernel_ridge_residuals(x, design, alpha=1.0, gamma=0.5,
                                                                   seed=42))):
        recovered = knn_balanced_accuracy_oof(adjusted, labels, 3, k=15, seed=42)
        assert recovered > chance + 0.15, f"{name} unexpectedly removed a variance-only confound"


def test_location_scale_arm_removes_what_no_residualiser_can():
    x, labels = _variance_only_confound()
    design = np.eye(3)[labels]
    chance = 1.0 / 3.0
    ridge = knn_balanced_accuracy_oof(cross_fitted_residuals(x, design, seed=42), labels, 3,
                                      k=15, seed=42)
    scaled = knn_balanced_accuracy_oof(cross_fitted_location_scale(x, labels, seed=42), labels, 3,
                                       k=15, seed=42)
    assert scaled < ridge - 0.10
    assert scaled < chance + 0.10


def test_location_scale_rejects_a_mismatched_cell_vector():
    x, labels = _variance_only_confound(n=120)
    with pytest.raises(ValueError):
        cross_fitted_location_scale(x, labels[:10], seed=42)


def test_location_scale_shrinkage_falls_back_to_pooled_for_an_unseen_cell():
    """A cell absent from a training fold must not produce a NaN or an infinity."""
    rng = np.random.default_rng(5)
    x = rng.normal(size=(60, 3))
    codes = np.zeros(60, dtype=np.int64)
    codes[:2] = 1                                  # a cell that some folds will not see
    adjusted = cross_fitted_location_scale(x, codes, n_splits=5, seed=42)
    assert np.isfinite(adjusted).all()


# --- the cross-fitting artefact -----------------------------------------------------------

def test_in_sample_residuals_zero_every_cell_mean_and_cross_fitted_ones_do_not():
    """The property the diagnostic turns on, asserted rather than argued."""
    rng = np.random.default_rng(7)
    n, n_cells = 400, 20
    codes = rng.integers(0, n_cells, n)
    x = rng.normal(size=(n, 5)) + np.eye(n_cells)[codes] @ rng.normal(size=(n_cells, 5)) * 2.0
    design = cell_design(codes)
    in_sample = in_sample_residuals(x, design)
    cross_fitted = saturated_cell_residuals(x, codes, seed=42)
    in_energy = cross_fitting_offset_energy(in_sample, codes, seed=42)
    cross_energy = cross_fitting_offset_energy(cross_fitted, codes, seed=42)
    assert in_energy["cell_mean_energy_fraction"] < 1e-12
    assert cross_energy["cell_mean_energy_fraction"] > 1e-4
    assert cross_energy["cell_fold_mean_energy_fraction"] > in_energy["cell_fold_mean_energy_fraction"]


def test_row_shuffle_preserves_geometry_and_destroys_the_association():
    x, _, design, codes, _ = _toy()
    shuffled = row_shuffled(x, seed=42)
    assert sorted(map(tuple, shuffled.tolist())) == sorted(map(tuple, x.tolist()))
    assert not np.array_equal(shuffled, x)


def test_a_pure_noise_block_can_still_be_probed_after_cross_fitted_adjustment():
    """Guard rail for the negative control: it must run and return a finite reading.

    The control's *value* on real data is a measurement and is not asserted here; what is
    asserted is that the path exists and produces a number, so that a run reporting it
    cannot be reporting a silently-skipped step.
    """
    rng = np.random.default_rng(2)
    n, n_cells = 500, 10
    codes = rng.integers(0, n_cells, n)
    noise = rng.normal(size=(n, 12))
    adjusted = saturated_cell_residuals(noise, codes, seed=42)
    reading = knn_balanced_accuracy_oof(adjusted, codes, n_cells, k=5, seed=42)
    assert np.isfinite(reading) and 0.0 <= reading <= 1.0


# --- reporting helpers --------------------------------------------------------------------

def test_adjuster_agreement_flags_a_relabelled_incumbent():
    x, _, design, _, _ = _toy()
    incumbent = cross_fitted_residuals(x, design, seed=42)
    nudged = incumbent + 1e-6 * np.random.default_rng(0).normal(size=incumbent.shape)
    assert adjuster_agreement(nudged, incumbent)["is_relabelled_incumbent"]
    different = cross_fitted_location_scale(x, np.arange(len(x)) % 8, seed=42)
    assert not adjuster_agreement(different, incumbent)["is_relabelled_incumbent"]


def test_adjuster_agreement_rejects_mismatched_shapes():
    x, _, design, _, _ = _toy()
    incumbent = cross_fitted_residuals(x, design, seed=42)
    with pytest.raises(ValueError):
        adjuster_agreement(incumbent, incumbent[:, :2])


def test_cross_fitted_r2_is_the_variance_an_arm_removes():
    x, _, design, codes, _ = _toy()
    saturated = cross_fitted_r2(x, make_adjuster("saturated", codes=codes, seed=42, alpha=1e-6))
    none = cross_fitted_r2(x, make_adjuster("none"))
    assert none["pooled_r2"] == pytest.approx(0.0, abs=1e-9)
    assert saturated["pooled_r2"] > 0.05


def test_labels_only_ceiling_reports_both_readouts():
    x, y, design, codes, cancers = _toy()
    record = labels_only_ceiling(design, y, n_components=4, seed=42,
                                 adjust=make_adjuster("ridge", design=design, seed=42),
                                 strata=cancers, n_permutations=10)
    assert 0.0 <= record["raw_top_cca"] <= 1.0
    assert np.isfinite(record["raw_heldout_top_cca"])
    assert record["adjusted"]["observed_top_cca"] < record["raw_top_cca"]


def test_make_adjuster_covers_every_declared_model_and_rejects_others():
    x, _, design, codes, _ = _toy()
    for model in ADJUSTER_MODELS:
        kwargs = {"n_estimators": 5} if model == "forest" else {}
        adjust = make_adjuster(model, design=design, codes=codes, seed=42, **kwargs)
        assert adjust(x).shape == x.shape
    with pytest.raises(ValueError):
        make_adjuster("gradient_boosting", design=design)
    with pytest.raises(ValueError):
        make_adjuster("ridge")
    with pytest.raises(ValueError):
        make_adjuster("saturated")
