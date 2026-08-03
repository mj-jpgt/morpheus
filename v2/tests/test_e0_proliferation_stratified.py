"""Guards for E0's proliferation stratification.

The arm whose survival IS the result is `responsive_nonprolif`. Every way of
accidentally inflating it is a way of manufacturing a discharge, so those are
what these tests target.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.calibra.e0_proliferation_stratified import (ARM_SEED_OFFSET,
                                                             STRATIFIED_SEED_OFFSET,
                                                             proliferation_gene_set,
                                                             stratify_responsive)


def test_every_arm_has_a_distinct_seed_offset():
    """Caught a real bug: a `.get(name, default[name])` whose default was evaluated
    eagerly and raised for the names that WERE present. Shared offsets would also
    correlate two arms' Haar draws and bootstrap resamples, which would make an
    arm comparison look more stable than it is."""
    required = {"responsive_matched", "nonresponsive", *STRATIFIED_SEED_OFFSET}
    assert required <= set(ARM_SEED_OFFSET), sorted(required - set(ARM_SEED_OFFSET))
    assert len(set(ARM_SEED_OFFSET.values())) == len(ARM_SEED_OFFSET)


def test_unknown_targets_are_excluded_from_both_strata_not_banked_as_non_proliferation():
    """Banking unparsed targets as non-proliferation would inflate exactly the arm
    whose survival is the result."""
    targets = ["MKI67", None, "ALB", "CCNB1", None, "GAPDH"]
    rows = np.arange(6)
    strata, meta = stratify_responsive(rows, targets, {"MKI67", "CCNB1"}, seed=1)
    assert sorted(strata["responsive_nonprolif"].tolist()) == [2, 5]
    assert sorted(strata["responsive_prolif_only"].tolist()) == [0, 3]
    assert meta["n_target_unknown"] == 2
    # The two unknowns appear in neither stratum.
    assert len(strata["responsive_nonprolif"]) + len(strata["responsive_prolif_only"]) == 4


def test_placebo_is_size_matched_to_the_stratified_arm_and_is_not_the_same_rows():
    """Without this the 'placebo' would not control for the n-match being redrawn."""
    rng = np.random.default_rng(0)
    targets = [f"G{i}" for i in range(200)]
    proliferation = {f"G{i}" for i in range(40)}          # 40 of 200 are proliferation
    rows = np.arange(200)
    strata, meta = stratify_responsive(rows, targets, proliferation, seed=7)
    assert len(strata["responsive_nonprolif"]) == 160
    assert len(strata["responsive_placebo"]) == 160
    assert meta["placebo_is_size_matched"] is True
    assert meta["n_dropped_for_placebo"] == 40
    # The placebo must not coincide with the real stratum, or it controls nothing.
    assert set(strata["responsive_placebo"].tolist()) != set(strata["responsive_nonprolif"].tolist())
    # ...and it must draw from the WHOLE arm, so it should include proliferation rows.
    assert len(set(strata["responsive_placebo"].tolist()) & set(range(40))) > 0


def test_missing_targets_yield_no_strata_rather_than_a_silent_proxy():
    strata, meta = stratify_responsive(np.arange(5), None, {"MKI67"}, seed=1)
    assert strata == {}
    assert meta["stratification_status"] == "unavailable_no_parsed_targets"
    strata, meta = stratify_responsive(np.arange(5), [], {"MKI67"}, seed=1)
    assert strata == {} and meta["stratification_status"] == "unavailable_no_parsed_targets"


def test_stratification_is_deterministic_for_a_fixed_seed():
    targets = [f"G{i}" for i in range(100)]
    a, _ = stratify_responsive(np.arange(100), targets, {"G1", "G2"}, seed=99)
    b, _ = stratify_responsive(np.arange(100), targets, {"G1", "G2"}, seed=99)
    assert np.array_equal(a["responsive_placebo"], b["responsive_placebo"])


def test_proliferation_gene_set_reads_the_flag_and_uppercases(tmp_path):
    path = tmp_path / "ann.parquet"
    pd.DataFrame({"gene": ["mki67", "ALB", "ccnb1"],
                  "proliferation_loading": [1.0, 0.0, 1.0],
                  "essentiality_loading": [0.1, 0.2, 0.3]}).to_parquet(path)
    assert proliferation_gene_set(path) == {"MKI67", "CCNB1"}


def test_an_empty_proliferation_set_leaves_the_arm_whole_which_must_be_visible():
    """A misconfigured annotation file would make the 'stratified' arm identical to
    the responsive arm and the discharge trivially true. The counts expose it."""
    targets = [f"G{i}" for i in range(50)]
    strata, meta = stratify_responsive(np.arange(50), targets, set(), seed=3)
    assert len(strata["responsive_nonprolif"]) == 50
    assert meta["n_proliferation_targets"] == 0
    assert meta["proliferation_target_fraction"] == pytest.approx(0.0)
    assert meta["n_dropped_for_placebo"] == 0
