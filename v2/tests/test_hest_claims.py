"""Tests for the spot-level claim-replication harness.

Two of these encode findings rather than mechanics.  ``test_knn_sees_what_lda_cannot``
constructs the exact situation claim 1c was predeclared to detect -- a confound that lives
in the class covariance rather than the class mean -- and asserts that the certificate's
mean-based classifier is blind to it while the kNN probe is not.  ``test_capacity_floor_...``
pins the 2*sqrt(k/n) law the spatial permutation null is graded against, so the prediction
and the measurement cannot drift apart in separate files.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.confound_certificate import lda_oof_balanced_accuracy
from morpheus.v2.calibra.hest import HestAdapterError, spot_key
from morpheus.v2.calibra.hest_claims import (between_slide_variance_share,
                                             capacity_floor_prediction,
                                             knn_balanced_accuracy_oof, matched_random_gene_sets,
                                             slide_bootstrap_interval, slide_stratified_subsample,
                                             spatial_confound_design)
from morpheus.v2.calibra.residualise import cross_fitted_residuals, pooled_tissue_source_site


def _slides(sizes):
    return np.concatenate([np.full(n, f"S{i}") for i, n in enumerate(sizes)])


# --- the declared subsample -------------------------------------------------------------

def test_subsample_is_slide_balanced_and_capped_by_the_smallest_slide():
    slides = _slides([50, 30, 12])
    mask = np.ones(len(slides), dtype=bool)
    rows = slide_stratified_subsample(slides, mask, 20, seed=42)
    counts = dict(zip(*np.unique(slides[rows], return_counts=True)))
    assert counts == {"S0": 20, "S1": 20, "S2": 12}
    assert (np.sort(rows) == rows).all()


def test_subsample_is_seed_reproducible_and_respects_the_base_mask():
    slides = _slides([40, 40])
    mask = np.zeros(len(slides), dtype=bool)
    mask[:20] = True
    mask[40:60] = True
    a = slide_stratified_subsample(slides, mask, 10, seed=7)
    b = slide_stratified_subsample(slides, mask, 10, seed=7)
    assert np.array_equal(a, b)
    assert np.array_equal(a, np.sort(a)) and mask[a].all()
    assert not np.array_equal(a, slide_stratified_subsample(slides, mask, 10, seed=8))


def test_subsample_zero_means_every_row_not_zero_rows():
    slides = _slides([5, 7])
    mask = np.ones(len(slides), dtype=bool)
    assert len(slide_stratified_subsample(slides, mask, 0)) == len(slides)


def test_subsample_refuses_misaligned_inputs():
    with pytest.raises(HestAdapterError):
        slide_stratified_subsample(_slides([4]), np.ones(3, dtype=bool), 2)


# --- the design really does resolve slide -----------------------------------------------

def test_spatial_design_resolves_the_slide_as_the_site_term():
    ids = np.asarray([spot_key(f"SLIDE{i}", f"BC{j}-1") for i in range(3) for j in range(12)])
    cancers = np.asarray(["COAD"] * 36)
    design, site, frequent = spatial_confound_design(ids, cancers, min_site_count=10)
    assert sorted(set(site.tolist())) == ["SLIDE0", "SLIDE1", "SLIDE2"]
    assert frequent == {"SLIDE0", "SLIDE1", "SLIDE2"}
    assert np.linalg.matrix_rank(design) == 3


def test_naive_spot_ids_would_have_collapsed_the_site_term():
    """The bug the key layout exists to prevent, pinned as an executable statement."""
    naive = np.asarray([f"SLIDE{i}__BC{j}-1" for i in range(3) for j in range(12)])
    site, _ = pooled_tissue_source_site(naive)
    assert len(set(site.tolist())) == 1          # every spot lands in one constant "site"


# --- claim 1c: the probe that is not mean-based -----------------------------------------

def test_knn_recovers_a_balanced_confound_and_matches_chance_on_noise():
    rng = np.random.default_rng(0)
    labels = np.repeat(np.arange(4), 60)
    features = rng.normal(size=(240, 6)) + 6.0 * np.eye(4)[labels] @ rng.normal(size=(4, 6))
    assert knn_balanced_accuracy_oof(features, labels, 4, k=5, seed=1) > 0.9
    noise = rng.normal(size=(240, 6))
    assert knn_balanced_accuracy_oof(noise, labels, 4, k=15, seed=1) < 0.45


def test_knn_sees_what_the_mean_based_certificate_cannot():
    """Claim 1c's premise, constructed: a confound carried by scale, not by location.

    Each class has the same mean and a different variance.  Removing class means -- which is
    exactly what ``cross_fitted_residuals`` on a one-hot design does -- leaves the confound
    fully intact.  LDA reads class means and shrunk pooled covariance, so it is near chance;
    a neighbourhood vote is not.
    """
    rng = np.random.default_rng(3)
    labels = np.repeat(np.arange(3), 200)
    scale = np.asarray([0.2, 1.0, 5.0])[labels][:, None]
    features = rng.normal(size=(600, 8)) * scale
    design = np.eye(3)[labels]
    adjusted = cross_fitted_residuals(features, design, seed=0)
    joint = lda_oof_balanced_accuracy(adjusted, labels, 3, seed=0)
    knn = knn_balanced_accuracy_oof(adjusted, labels, 3, k=15, seed=0)
    # measured on this construction: joint 0.417, kNN 0.653, chance 0.333. The bars are set
    # around those rather than at round numbers, because the point of the test is the GAP.
    assert joint < 0.45, joint
    assert knn > 0.58, knn
    assert knn - joint > 0.20, (knn, joint)


def test_knn_refuses_a_k_it_cannot_honour():
    rng = np.random.default_rng(0)
    with pytest.raises(HestAdapterError):
        knn_balanced_accuracy_oof(rng.normal(size=(10, 3)), np.repeat(np.arange(2), 5), 2, k=15)


# --- claim 2: the predeclared law -------------------------------------------------------

def test_capacity_floor_prediction_matches_the_predeclared_table():
    assert capacity_floor_prediction(2766, 16) == pytest.approx(0.1521, abs=5e-4)
    assert capacity_floor_prediction(53217, 16) == pytest.approx(0.0347, abs=5e-4)
    assert capacity_floor_prediction(144162, 16) == pytest.approx(0.0211, abs=5e-4)


def test_capacity_floor_shrinks_with_n_and_grows_with_components():
    assert capacity_floor_prediction(10000, 16) < capacity_floor_prediction(1000, 16)
    assert capacity_floor_prediction(1000, 32) > capacity_floor_prediction(1000, 16)
    # the design eats degrees of freedom, so a larger design predicts a HIGHER floor
    assert capacity_floor_prediction(1000, 16, design_rank=200) > capacity_floor_prediction(1000, 16)


def test_capacity_floor_is_close_to_a_measured_null_on_independent_blocks():
    """The law is not decoration: an actual top-CCA on independent data must land near it."""
    from morpheus.v2.calibra.spectral import top_canonical_correlation

    rng = np.random.default_rng(11)
    values = [top_canonical_correlation(rng.normal(size=(4000, 40)), rng.normal(size=(4000, 40)),
                                        n_components=16) for _ in range(5)]
    predicted = capacity_floor_prediction(4000, 16)
    assert abs(float(np.median(values)) - predicted) < 0.2 * predicted


# --- claim 3: matched random gene sets --------------------------------------------------

def test_matched_random_sets_avoid_the_panel_are_disjoint_and_track_the_panel_statistics():
    rng = np.random.default_rng(5)
    names = np.asarray([f"G{i}" for i in range(600)])
    mean = rng.uniform(0.0, 4.0, size=600)
    logvar = 0.8 * mean + rng.normal(scale=0.1, size=600)
    panel = names[np.argsort(-logvar)[:50]]
    sets = matched_random_gene_sets(mean, logvar, names, panel, n_sets=3, seed=101)
    panel_rows = {int(np.flatnonzero(names == g)[0]) for g in panel}
    for chosen in sets:
        assert len(chosen) == 50 and len(set(chosen.tolist())) == 50
        assert not (set(chosen.tolist()) & panel_rows)
        assert abs(mean[chosen].mean() - mean[list(panel_rows)].mean()) < 0.5
    assert not np.array_equal(sets[0], sets[1])


def test_matched_random_sets_reject_a_panel_gene_that_has_no_statistics():
    names = np.asarray([f"G{i}" for i in range(100)])
    with pytest.raises(HestAdapterError):
        matched_random_gene_sets(np.zeros(100), np.zeros(100), names, np.asarray(["MISSING"]))


# --- claim 5: intervals and the variance decomposition ----------------------------------

def test_slide_bootstrap_interval_brackets_the_point_and_counts_slides():
    values = {f"S{i}": np.asarray([0.1 * i, 0.1 * i + 0.02]) for i in range(13)}
    out = slide_bootstrap_interval(values, n_boot=500, seed=1)
    assert out["n_slides"] == 13
    assert out["ci95"][0] < out["point"] < out["ci95"][1]


def test_pooled_r_of_the_slide_mean_baseline_is_the_root_of_the_between_slide_share():
    """The identity that makes claim 5 quotable rather than rhetorical.

    A per-slide-mean predictor's pooled Pearson r against the truth is exactly
    ``sqrt(between-slide variance share)``.  Correlation is the square root of a variance
    ratio, so a slide effect explaining only a third of the variance still yields a pooled r
    near 0.6.  Measured on the real test partition: share 0.3568, pooled r 0.5973, and
    sqrt(0.3568) = 0.5973.
    """
    from morpheus.v2.calibra.hest import per_slide_mean_baseline, pooled_r

    rng = np.random.default_rng(9)
    slides = _slides([80, 120, 60])
    offsets = np.repeat(rng.normal(scale=1.3, size=(3, 4)), [80, 120, 60], axis=0)
    targets = offsets + rng.normal(size=(260, 4))
    mask = np.ones(len(slides), dtype=bool)
    prediction = per_slide_mean_baseline(targets, slides, mask)
    share = between_slide_variance_share(targets, slides, mask)
    np.testing.assert_allclose(pooled_r(targets, prediction, mask), np.sqrt(share), atol=1e-12)


def test_between_slide_variance_share_is_one_for_a_pure_slide_offset_and_zero_without_one():
    slides = _slides([50, 50, 50])
    offsets = np.repeat(np.asarray([[-2.0], [0.0], [2.0]]), 50, axis=0)
    assert between_slide_variance_share(offsets, slides)[0] == pytest.approx(1.0)
    rng = np.random.default_rng(2)
    noise = rng.normal(size=(150, 4))
    assert float(np.nanmax(between_slide_variance_share(noise, slides))) < 0.15
