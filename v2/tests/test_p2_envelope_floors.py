"""The measured floors must agree with the two numbers §4.1 already published, and
the absent ones must stay named.

`v2/research/rebase/p2/p2_envelope_floors.py` measures the same-seed retraining
floor for every statistic, view and block recoverable from the five exported
repeats. Three things have to hold or the measurement is not usable:

* **it reproduces §4.1.** The R1 floors it computes from the artifacts must be the
  3.295× / 3.111× §4.1 parses out of `d1_envelope_readout.py`'s printed log. Two
  independent paths to one number is the point; if they part, neither may be used;
* **no statistic is reimplemented.** R1/R2/R3 resolve to the one `effective_rank`
  object, and the published alternatives resolve to the ones §4.6 scores;
* **the absent floors stay absent and stay named.** A block with no floor is a
  recorded result with a cost attached, not a silence.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2 import p2_competing_metrics as C
from morpheus.v2.research.rebase.p2 import p2_envelope_floors as F
from morpheus.v2.research.rebase.p2 import p2_floor_audit as A

MEASURED = A.DATA / "ws_floor" / "out" / "P2_ENVELOPE_FLOORS.json"


@pytest.fixture(scope="module")
def measured():
    return json.loads(MEASURED.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------
# It reproduces the floor the paper already publishes
# --------------------------------------------------------------------------
def test_the_r1_floors_reproduce_the_ones_section_4_1_publishes(measured):
    """Two independent paths to 3.295× and 3.111×.

    §4.1's numbers are parsed out of the box readout log by the figure
    extractor. These are recomputed from the five `.npz` artifacts on a
    workspace verified file-by-file against `git ls-tree`. A disagreement is
    reported, never resolved by preferring one of them -- which is why this is a
    test and not a comment.
    """
    floors = measured["floors"]["wsi_biology"]
    assert floors["residualised"]["R1"]["fold"] == pytest.approx(3.295, abs=0.0015)
    assert floors["raw"]["R1"]["fold"] == pytest.approx(3.111, abs=0.0015)
    assert measured["floors"]["channel"]["residualised"]["top_cca_16"]["fold"] == pytest.approx(
        1.055, abs=0.0015)

    printed = json.loads((A.DATA / "extracted" / "F1_RETRAINING_REPEAT.json")
                         .read_text(encoding="utf-8"))["printed_spread"]
    for key, block in (("rank_residualised", "residualised"), ("rank_raw", "raw")):
        for end in ("min", "max"):
            assert floors[block]["R1"][end] == pytest.approx(printed[key][end], abs=5e-4), (
                f"the readout log and the recomputation disagree on {key} {end} -- "
                "STOP and report both, do not use either")


def test_the_five_repeats_are_the_five_the_paper_reports(measured):
    assert sorted(measured["reps"]) == ["rep1", "rep2", "rep3", "rep4", "rep5"]
    for rep in measured["reps"].values():
        assert rep["n_test"] == 2766
        assert len(rep["sha256"]) == 64


# --------------------------------------------------------------------------
# Nothing is reimplemented
# --------------------------------------------------------------------------
def test_every_statistic_is_the_imported_one():
    """Object identity where it is available, value identity where it is not.

    Four inline-formula substitutions have been caught in this paper and the
    tell was identical each time: an arithmetic expression where an import
    belonged.
    """
    assert F.effective_rank is effective_rank
    assert F.participation_ratio is C.participation_ratio
    assert F.participation_ratio_rownorm is C.participation_ratio_rownorm
    assert F.rankme is C.rankme
    assert F.stable_rank is C.stable_rank
    assert F.lidar is C.lidar
    assert F.alpha_req is C.alpha_req

    rng = np.random.default_rng(11)
    sample = rng.normal(size=(120, 24))
    for name in ("R1", "R2", "R3"):
        assert F.STATISTICS[name](sample) == pytest.approx(
            effective_rank(sample, variant=RANK_VARIANTS[name]), rel=1e-12)
    assert F.STATISTICS["hard_rank"](sample) == float(np.linalg.matrix_rank(sample))


def test_alpha_and_lidar_are_measured_at_the_settings_t1_is_scored_with(measured):
    """A floor measured at other settings would not be a floor for T1's numbers."""
    metrics = json.loads((A.DATA / "ws_p2" / "out" / "P2_METRICS_D2.json")
                         .read_text(encoding="utf-8"))["_config"]
    assert list(F.ALPHA_INDEX_RANGE) == metrics["alpha_index_range"]
    assert F.LIDAR_DELTA == metrics["lidar_delta"]
    assert measured["config"]["alpha_index_range"] == metrics["alpha_index_range"]
    assert measured["config"]["lidar_delta"] == metrics["lidar_delta"]


# --------------------------------------------------------------------------
# The shape rule
# --------------------------------------------------------------------------
def test_the_shape_rule_reproduces_section_4_1s_own_description():
    """"Four repeats agree to within 2% and one lands at a third of them.\""""
    shape = F._shape({"rep1": 28.3202, "rep2": 8.8340, "rep3": 28.3482,
                      "rep4": 29.1057, "rep5": 28.9588})
    assert shape["outlier"] == "rep2" and shape["outlier_is_low"]
    assert shape["rest_fold"] == pytest.approx(1.0277, abs=5e-4)
    assert shape["bimodal"] is True


def test_a_smooth_spread_is_not_called_bimodal():
    """The rule has to be able to say no, or reporting it per statistic is theatre."""
    assert F._shape({f"rep{i}": 1.0 + 0.1 * i for i in range(5)})["bimodal"] is False


def test_a_non_finite_value_yields_no_fold_rather_than_a_number():
    entry = F._fold({"rep1": 1.0, "rep2": float("nan"), "rep3": 2.0, "rep4": 2.0, "rep5": 2.0})
    assert entry["fold"] is None and entry["min"] is None
    assert entry["shape"]["defined"] is False


# --------------------------------------------------------------------------
# The findings the draft leans on, asserted against the measurement
# --------------------------------------------------------------------------
def test_the_floor_depends_on_the_statistic(measured):
    """§4.3's withdrawn claim was that it does not. On one block it spans 3.3×."""
    block = measured["floors"]["wsi_biology"]["residualised"]
    assert block["hard_rank"]["fold"] == 1.0
    assert block["R1"]["fold"] / block["hard_rank"]["fold"] > 3.0
    assert block["R1"]["fold"] > block["R3"]["fold"] > block["PR"]["fold"] \
        > block["stable_rank"]["fold"]


def test_the_floor_depends_on_the_view_and_the_divergent_run_is_wsi_only(measured):
    """The catastrophic one-in-five is a property of that run's WSI encoder.

    On the same five artifacts the R1 floor is 3.295× on `wsi_biology`, 1.019× on
    `rna_biology` and 1.020× on `full_biology`, and repeat 2 -- whose WSI-view
    rank is a third of its siblings' -- is not separated from them on either
    other view.
    """
    floors = measured["floors"]
    assert floors["wsi_biology"]["residualised"]["R1"]["fold"] > 3.0
    for view in ("rna_biology", "full_biology"):
        entry = floors[view]["residualised"]["R1"]
        assert entry["fold"] < 1.05, view
        assert entry["shape"]["bimodal"] is False, view
    assert floors["wsi_biology"]["residualised"]["R1"]["shape"]["bimodal"] is True


def test_published_rankme_is_more_reproducible_than_our_centred_statistic(measured):
    """A result against our own instrument, on our own artifacts.

    RankMe's uncentred normalisation retains the mean-offset direction, and every
    exported row has L2 norm 1, so that direction is large and stable. On the
    residualised block, where the mean is gone, the two statistics coincide.
    """
    raw = measured["floors"]["wsi_biology"]["raw"]
    residualised = measured["floors"]["wsi_biology"]["residualised"]
    assert raw["RankMe"]["fold"] < 0.7 * raw["R1"]["fold"]
    assert residualised["RankMe"]["fold"] == pytest.approx(
        residualised["R1"]["fold"], abs=0.0015)


def test_the_exported_rows_are_unit_norm_which_is_why_r2_and_r3_coincide_raw(measured):
    """R3 is R2 on L2-normalised rows, so on a unit-norm block they are one number.

    Recorded because it is the mechanism behind the RankMe result above, and
    because two statistics agreeing exactly is otherwise the signature of a bug.
    """
    for view in ("wsi_biology", "rna_biology", "full_biology"):
        raw = measured["floors"][view]["raw"]
        for a, b in (("R2", "R3"), ("PR", "PR_rownorm")):
            for rep, value in raw[a]["values"].items():
                # Agreement to float noise, not to the bit: the row-normalisation
                # still divides by a norm, it just divides by 1.
                assert value == pytest.approx(raw[b]["values"][rep], rel=1e-7), (view, a, rep)
        residualised = measured["floors"][view]["residualised"]
        assert residualised["R2"]["fold"] != residualised["R3"]["fold"], view


def test_the_divergent_run_is_visible_to_every_statistic_that_moves(measured):
    """The shape is statistic-dependent; the divergence is not hidden from any of them.

    A statistic in which the divergent run was NOT the outlier would mean the
    divergence is invisible to some measures, which would be a different and
    larger finding. It is not what the data say: every statistic on the WSI block
    that moves at all puts repeat 2 at the extreme, in the degradation direction.
    """
    for block in ("raw", "residualised"):
        for name, entry in measured["floors"]["wsi_biology"][block].items():
            if name == "hard_rank":          # does not move at all; floor 1.000x
                assert entry["fold"] == 1.0
                continue
            assert entry["shape"]["outlier"] == "rep2", (block, name)
            # lower rank, or a STEEPER alpha-ReQ decay exponent: both are the
            # degradation direction for the statistic in question.
            expect_low = not name.startswith("alpha_req")
            assert entry["shape"]["outlier_is_low"] is expect_low, (block, name)


# --------------------------------------------------------------------------
# The floors that do NOT exist
# --------------------------------------------------------------------------
def test_the_absent_blocks_are_named_with_what_each_would_cost(measured):
    """A named absent floor is a result. It has to survive in a file, not in prose."""
    absent = measured["absent_blocks"]
    assert set(absent) == set(F.ABSENT_BLOCKS)
    assert "fixed held-out probe" in absent
    for block, why in absent.items():
        assert len(why) > 80, block
        assert "GPU" in why or "re-created" in why, block
    assert "§5" in absent["fixed held-out probe"], (
        "the block every rank number in draft §5 sits on must say so")


def test_the_audit_and_the_measurement_agree_on_which_blocks_have_no_floor(measured):
    """One list, two files, asserted equal rather than maintained twice."""
    audit = A.load()
    assert set(audit["blocks_with_no_measured_floor"]) == set(measured["absent_blocks"])
    assert audit["statistics_with_no_measured_floor"] == [], (
        "every statistic T1 scores now has a floor on the exported block")
