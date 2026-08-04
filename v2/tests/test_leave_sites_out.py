"""Tests for the leave-sites-out generalisation test.

A generalisation test is only evidence if it can return "collapses" when the
channel really does not transfer, and "survives" when it does. Each test below
plants one of those two situations and asserts the instrument reacts. The fold
builder gets its own tests because a site that straddles a fold would leave the
same scanner on both sides of the split and quietly make the test trivial.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.leave_sites_out import (evaluate_fold, fold_composition,
                                                 leave_sites_out_channel,
                                                 matched_random_folds, site_folds)
from morpheus.v2.calibra.residualise import pooled_tissue_source_site
from morpheus.v2.calibra.spectral import (heldout_cca_projection, heldout_top_cca,
                                          heldout_top_cca_indexed)


def _cohort(n_per_site: int = 24, n_sites: int = 8, n_cancers: int = 3, seed: int = 0):
    """TCGA-shaped barcodes whose TSS field is real, with site NESTED IN CANCER."""
    rng = np.random.default_rng(seed)
    ids, cancers, sites = [], [], []
    site_number = 0
    for c in range(n_cancers):
        for _ in range(n_sites):
            code = f"{site_number:02d}"
            site_number += 1
            for _ in range(n_per_site):
                ids.append(f"TCGA-{code}-{rng.integers(1000, 9999)}{len(ids)}")
                cancers.append(f"CANCER{c}")
                sites.append(code)
    return np.asarray(ids), np.asarray(cancers), np.asarray(sites)


# --------------------------------------------------------------------------- folds

def test_site_folds_never_straddles_a_site():
    """The defining property: a site is wholly held out or wholly trained on."""
    _, cancers, sites = _cohort()
    fold = site_folds(sites, cancers, n_folds=5)
    for site in np.unique(sites):
        assert len(set(fold[sites == site].tolist())) == 1, f"site {site} straddles folds"


def test_site_folds_is_deterministic_and_order_invariant():
    """Seed-free: no split a caller could re-roll until it liked the answer."""
    _, cancers, sites = _cohort()
    first = site_folds(sites, cancers, n_folds=5)
    assert np.array_equal(first, site_folds(sites, cancers, n_folds=5))
    order = np.random.default_rng(1).permutation(len(sites))
    shuffled = site_folds(sites[order], cancers[order], n_folds=5)
    # Same site -> same fold, whatever order the rows arrived in.
    by_site = {s: int(first[sites == s][0]) for s in np.unique(sites)}
    for site in np.unique(sites):
        assert len(set(shuffled[sites[order] == site].tolist())) == 1
        assert int(shuffled[sites[order] == site][0]) == by_site[site]


def test_site_folds_keeps_every_cancer_on_both_sides():
    """Balanced WITHIN cancer, so this stays a site test and not a cancer test."""
    _, cancers, sites = _cohort(n_sites=8, n_cancers=3)
    fold = site_folds(sites, cancers, n_folds=5)
    for f in range(5):
        held = fold == f
        assert held.any()
        assert set(cancers[held].tolist()) == set(cancers[~held].tolist())


def test_site_folds_balances_patients_across_folds():
    _, cancers, sites = _cohort(n_per_site=24, n_sites=10, n_cancers=2)
    fold = site_folds(sites, cancers, n_folds=5)
    sizes = np.asarray([(fold == f).sum() for f in range(5)])
    assert sizes.max() - sizes.min() <= 0.5 * sizes.mean()


def test_site_folds_rejects_bad_input():
    _, cancers, sites = _cohort()
    with pytest.raises(ValueError):
        site_folds(sites, cancers[:-1], n_folds=5)
    with pytest.raises(ValueError):
        site_folds(sites, cancers, n_folds=1)


def test_matched_random_folds_preserve_cancer_and_fold_sizes():
    """The comparator must differ from the site split ONLY in respecting sites."""
    _, cancers, sites = _cohort()
    fold = site_folds(sites, cancers, n_folds=5)
    matched = matched_random_folds(cancers, fold, seed=3)
    for cancer in np.unique(cancers):
        for f in range(5):
            assert (((cancers == cancer) & (fold == f)).sum()
                    == ((cancers == cancer) & (matched == f)).sum())


def test_matched_random_folds_do_break_sites():
    """If the comparator respected sites it would not be a comparator."""
    _, cancers, sites = _cohort()
    fold = site_folds(sites, cancers, n_folds=5)
    matched = matched_random_folds(cancers, fold, seed=5)
    straddling = sum(1 for s in np.unique(sites) if len(set(matched[sites == s].tolist())) > 1)
    assert straddling > 0


def test_fold_composition_counts_add_up():
    _, cancers, sites = _cohort()
    fold = site_folds(sites, cancers, n_folds=5)
    rows = fold_composition(fold, sites, cancers)
    assert len(rows) == 5
    assert sum(r["n_heldout_patients"] for r in rows) == len(sites)
    for row in rows:
        assert row["n_heldout_patients"] + row["n_training_patients"] == len(sites)
        assert row["n_heldout_sites"] >= 1
        assert 0.0 < row["largest_heldout_site_share"] <= 1.0


# ------------------------------------------------------------------- the estimator

def test_heldout_top_cca_still_equals_its_indexed_core():
    """The refactor must not move the published statistic by a single bit."""
    rng = np.random.default_rng(0)
    x = rng.normal(size=(300, 12))
    y = x @ rng.normal(size=(12, 8)) + 0.5 * rng.normal(size=(300, 8))
    order = np.random.default_rng(42).permutation(300)
    cut = int(300 * 0.5)
    assert heldout_top_cca(x, y, n_components=6, seed=42) == heldout_top_cca_indexed(
        x, y, order[:cut], order[cut:], n_components=6)


def test_heldout_projection_is_the_statistic():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(200, 10))
    y = x @ rng.normal(size=(10, 6)) + 0.3 * rng.normal(size=(200, 6))
    train, test = np.arange(120), np.arange(120, 200)
    px, py = heldout_cca_projection(x, y, train, test, n_components=5)
    assert px.shape == py.shape == (80,)
    assert abs(np.corrcoef(px, py)[0, 1]) == pytest.approx(
        heldout_top_cca_indexed(x, y, train, test, n_components=5))


def test_indexed_cca_is_near_zero_on_unpaired_noise():
    """Held-out, so capacity cannot manufacture a correlation the way in-sample does."""
    rng = np.random.default_rng(2)
    x, y = rng.normal(size=(400, 20)), rng.normal(size=(400, 20))
    value = heldout_top_cca_indexed(x, y, np.arange(200), np.arange(200, 400), n_components=10)
    assert value < 0.3


def test_indexed_cca_recovers_a_planted_signal():
    rng = np.random.default_rng(3)
    shared = rng.normal(size=(400, 1))
    x = np.hstack([shared, rng.normal(size=(400, 9))])
    y = np.hstack([shared, rng.normal(size=(400, 5))])
    value = heldout_top_cca_indexed(x, y, np.arange(200), np.arange(200, 400), n_components=5)
    assert value > 0.8


def test_too_small_a_split_is_unavailable_not_wrong():
    rng = np.random.default_rng(4)
    x, y = rng.normal(size=(12, 4)), rng.normal(size=(12, 4))
    assert np.isnan(heldout_top_cca_indexed(x, y, np.arange(6), np.arange(6, 12), n_components=2))


# ----------------------------------------------------------------- the whole test

def test_channel_collapses_on_unpaired_noise():
    """Must-return-collapses: no pairing exists, so no fold may survive."""
    ids, cancers, sites = _cohort(n_per_site=20, n_sites=6, n_cancers=2, seed=7)
    rng = np.random.default_rng(8)
    x, y = rng.normal(size=(len(ids), 24)), rng.normal(size=(len(ids), 12))
    result = leave_sites_out_channel(x, y, ids, cancers, n_folds=4, n_components=6,
                                     n_permutations=60, n_boot=60, seed=11)
    assert result["verdict"] == "collapses"
    assert result["n_folds_surviving"] == 0


def test_channel_survives_a_site_independent_planted_signal():
    """Must-return-survives: pairing that does not depend on site must transfer."""
    ids, cancers, sites = _cohort(n_per_site=40, n_sites=6, n_cancers=2, seed=9)
    rng = np.random.default_rng(10)
    n = len(ids)
    shared = rng.normal(size=(n, 1))
    x = np.hstack([shared, rng.normal(size=(n, 15))])
    y = np.hstack([shared, rng.normal(size=(n, 7))])
    result = leave_sites_out_channel(x, y, ids, cancers, n_folds=4, n_components=5,
                                     n_permutations=60, n_boot=60, seed=12)
    assert result["verdict"] in {"survives", "attenuated_but_present"}
    assert result["n_folds_surviving"] >= 3


def test_a_purely_site_specific_signal_does_not_transfer():
    """The situation this whole test exists to catch.

    Each site gets its OWN pairing direction, so the channel is real and strong
    in-sample and worth nothing on a site never seen. A test that scored this as
    survival would be measuring nothing.
    """
    ids, cancers, sites = _cohort(n_per_site=30, n_sites=6, n_cancers=2, seed=13)
    rng = np.random.default_rng(14)
    n, d = len(ids), 12
    x = rng.normal(size=(n, d))
    y = np.zeros((n, d))
    for site in np.unique(sites):
        rows = np.flatnonzero(sites == site)
        rotation = np.linalg.qr(rng.normal(size=(d, d)))[0]
        y[rows] = x[rows] @ rotation
    result = leave_sites_out_channel(x, y, ids, cancers, n_folds=4, n_components=5,
                                     n_permutations=60, n_boot=60, seed=15)
    matched = [f["observed"] for f in result["matched_random_folds"] if f.get("status") == "scored"]
    # The comparator splits patients at random, so most sites appear on BOTH sides
    # and the per-site rotations are learnable: it must score higher than the
    # site-respecting split. That gap is exactly what the comparator is for.
    assert np.nanmedian(matched) > result["median_heldout_site_channel"]
    assert result["verdict"] == "collapses"


def test_adjusted_arm_runs_and_is_labelled():
    ids, cancers, sites = _cohort(n_per_site=30, n_sites=5, n_cancers=2, seed=16)
    rng = np.random.default_rng(17)
    n = len(ids)
    shared = rng.normal(size=(n, 1))
    x = np.hstack([shared, rng.normal(size=(n, 11))])
    y = np.hstack([shared, rng.normal(size=(n, 5))])
    result = leave_sites_out_channel(x, y, ids, cancers, n_folds=4, n_components=5,
                                     adjust=True, min_site_count=1,
                                     n_permutations=40, n_boot=40, seed=18)
    assert result["adjusted"] is True
    assert result["adjustment"] == "cancer+pooled_tss_cross_fitted_within_block"
    assert result["status"] == "scored"


def test_the_test_never_claims_to_discharge_the_external_cohort_blocker():
    """The one thing this result must never be read as buying."""
    ids, cancers, _ = _cohort(n_per_site=20, n_sites=5, n_cancers=2, seed=19)
    rng = np.random.default_rng(20)
    x, y = rng.normal(size=(len(ids), 10)), rng.normal(size=(len(ids), 6))
    result = leave_sites_out_channel(x, y, ids, cancers, n_folds=4, n_components=4,
                                     n_permutations=30, n_boot=30, seed=21)
    assert result["discharges_no_external_cohort"] is False
    assert result["test_name"] == "held-out-site"
    assert result["not_this_test"] == "held-out-cancer"


def test_fold_result_reports_null_and_both_intervals():
    ids, cancers, sites = _cohort(n_per_site=25, n_sites=5, n_cancers=2, seed=22)
    rng = np.random.default_rng(23)
    n = len(ids)
    shared = rng.normal(size=(n, 1))
    x = np.hstack([shared, rng.normal(size=(n, 9))])
    y = np.hstack([shared, rng.normal(size=(n, 5))])
    pooled, _ = pooled_tissue_source_site(ids, min_site_count=1)
    fold = site_folds(sites, cancers, n_folds=4)
    result = evaluate_fold(x, y, cancers, sites, pooled, fold == 0, n_components=4,
                           adjust=False, n_permutations=40, n_boot=40, seed=24)
    assert result["status"] == "scored"
    assert result["permutation_resolution"] == pytest.approx(1 / 41)
    assert 0.0 <= result["permutation_p"] <= 1.0
    assert len(result["site_cluster_ci95"]) == 2
    assert len(result["patient_ci95"]) == 2
    assert result["site_cluster_ci95"][0] <= result["site_cluster_ci95"][1]


def test_tiny_fold_is_reported_unavailable_not_silently_dropped():
    ids, cancers, sites = _cohort(n_per_site=3, n_sites=2, n_cancers=1, seed=25)
    rng = np.random.default_rng(26)
    x, y = rng.normal(size=(len(ids), 6)), rng.normal(size=(len(ids), 4))
    pooled, _ = pooled_tissue_source_site(ids, min_site_count=1)
    fold = site_folds(sites, cancers, n_folds=2)
    result = evaluate_fold(x, y, cancers, sites, pooled, fold == 0, n_components=3,
                           adjust=False, n_permutations=10, n_boot=10, seed=27)
    assert result["status"].startswith("unavailable_")
