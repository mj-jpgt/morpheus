"""§5.4 limit 2's pass must survive being pushed on, and the pushing must be honest.

`p2_limit2_stress.py` exists to break, or fail to break, the one selection this
paper makes at the hyperparameter the project actually ships — `m = 0.999` over
`m = 0.99` on the fixed held-out probe at step 600. Three things have to hold or
the pushing is theatre:

* **no statistic is reimplemented.** Every entry of `STATS` resolves to the one
  `effective_rank` object or to the `p2_competing_metrics` callable §4.6 scores,
  and the extra variants are keys of `RANK_VARIANTS` rather than inventions;
* **duplicates are not counted as corroboration.** `z_biology` is L2-normalised
  at the model's output, so R2 = R3 and R1 = R1_rownorm and PR = PR_rownorm on
  this block. "Six statistics agree" would be a lie if three of them are the
  other three;
* **the direction of the evidence is labelled.** Both sides of Test B move
  against the pass as repeats are added, so a failure at larger n is arithmetic
  and a survival is evidence. The module has to say which it is reporting.

The last test is the one that matters: it asserts the monotonicity claim the
whole write-up rests on, directly against the measured n = 5 and n = 10 files.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from morpheus.v2.calibra.spectral import CANONICAL, RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2 import p2_envelope_floors as F
from morpheus.v2.research.rebase.p2 import p2_floor_audit as A
from morpheus.v2.research.rebase.p2 import p2_limit2_stress as S

N5 = A.DATA / "e0_run" / "d1_probefloor600" / "out" / "P2_LIMIT2_STRESS_N5.json"
N10 = A.DATA / "e0_run" / "d1_probefloor600" / "out" / "P2_LIMIT2_STRESS_N10.json"
LATE5 = A.DATA / "e0_run" / "d1_probefloor600" / "out" / "P2_LIMIT2_STRESS_LATE5.json"
GRID = A.DATA / "e0_run" / "d1_probefloor600" / "out" / "P2_MOMENTUM_GRID.json"


def _load(path):
    if not path.is_file():
        pytest.skip(f"{path.name} not vendored")
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def n5():
    return _load(N5)


@pytest.fixture(scope="module")
def n10():
    return _load(N10)


# --------------------------------------------------------------------------
# Nothing is reimplemented
# --------------------------------------------------------------------------
def test_every_statistic_resolves_to_an_existing_definition():
    """`STATS` is `STATISTICS` plus three keys of `RANK_VARIANTS`, and nothing else.

    The failure this guards is the one that has been caught four times in this
    paper: an arithmetic expression where an import belonged. A fourth statistic
    once lived under the R2/R3 names in an analysis script, changed a published
    count, and survived review, the tests and its own author.
    """
    assert set(S.STATS) == set(F.STATISTICS) | set(S.EXTRA_VARIANTS)
    for name in S.EXTRA_VARIANTS:
        assert name in RANK_VARIANTS, f"{name} is not a key of RANK_VARIANTS -- it was invented"
    # And they really are the one function: a random matrix scored through STATS
    # must equal the same variant scored directly.
    x = np.random.default_rng(0).normal(size=(64, 32))
    for name in S.EXTRA_VARIANTS:
        assert S.STATS[name](x) == effective_rank(x, variant=RANK_VARIANTS[name])
    assert RANK_VARIANTS["R1"] == CANONICAL


def test_the_module_computes_no_rank_inline():
    """The tokens come from the tree-wide guard, not from a second list of them.

    `RANK_SHAPED` is imported rather than restated for two reasons. It is the one
    place this repository defines what an inline spectrum looks like, so a token
    added there is enforced here for free. And writing the literals out would put
    them in THIS file's source, which `test_no_second_definition_exists_in_the_tree`
    scans -- a check for inline SVD that trips the check for inline SVD. That is
    not hypothetical: it is what happened when this test was first written.
    """
    from morpheus.v2.tests.test_effective_rank_canonical import RANK_SHAPED

    src = (A.REPO / "v2" / "research" / "rebase" / "p2" / "p2_limit2_stress.py").read_text(
        encoding="utf-8")
    assert "from morpheus.v2.calibra.spectral import" in src
    for token in RANK_SHAPED:
        assert token not in src, (
            f"`{token}` in p2_limit2_stress.py -- a spectrum computed here would be a "
            "second definition of effective rank")


# --------------------------------------------------------------------------
# Duplicates are detected, not assumed
# --------------------------------------------------------------------------
def test_row_normalised_variants_duplicate_their_bases_on_row_unit_data():
    """On L2-normalised rows, `normalise_rows` is a no-op -- checked, not asserted.

    This is the property that makes R2 and R3 the same statistic on the probe
    block, which §3.1 has to record because §3.1's whole point is that three
    statistics travel under one name.
    """
    x = np.random.default_rng(1).normal(size=(128, 40))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    for base, rownorm in (("R1", "R1_rownorm"), ("R2", "R3")):
        assert S.STATS[rownorm](x) == pytest.approx(S.STATS[base](x), rel=S.DUPLICATE_REL)
    # ... and NOT a no-op when the rows carry genuine norm variation.
    y = np.random.default_rng(2).normal(size=(128, 40))
    y[:20] *= 50.0
    assert S.STATS["R3"](y) != pytest.approx(S.STATS["R2"](y), rel=1e-3)


def test_the_measured_duplicates_are_the_ones_the_block_predicts(n5):
    dup = n5["duplicates_on_this_block"]["duplicate_of"]
    assert dup["R3"] == "R2" and dup["R1_rownorm"] == "R1" and dup["PR_rownorm"] == "PR"
    # RankMe is exp(Shannon entropy) of the UNCENTRED spectrum with an eps the
    # RankMe paper adds outside the normalisation, so it is R1_uncentred plus a
    # perturbation -- near, not equal, and it must not be counted as a second
    # statistic agreeing with R1_uncentred.
    assert dup["R1_uncentred"] == "RankMe"
    near = n5["duplicates_on_this_block"]["near_duplicates"]["R1_uncentred ~= RankMe"]
    assert S.DUPLICATE_REL < near < S.NEAR_DUPLICATE_REL


# --------------------------------------------------------------------------
# Test A exists only under R3, and it agrees with the audit it re-tests
# --------------------------------------------------------------------------
def test_test_A_is_R3_only_and_matches_the_audit_row():
    """The row can only be judged under the statistic it FAILS the strong test on.

    `long_m0.999.log` and `long_m0.99.log` predate both the canonical column and
    `export_dir`, so no other statistic has a published pair for this row. That
    is not a gap in this module; it is a property of the evidence.
    """
    assert set(S.PUBLISHED_TEST_A) == {"R3"}
    row = next(c for c in A.load()["comparisons"] if c["id"] == "5.4-m0999-over-m099")
    assert row["statistic"] == "R3"
    assert S.PUBLISHED_TEST_A["R3"]["high"] == row["a"]["value"]
    assert S.PUBLISHED_TEST_A["R3"]["low"] == row["b"]["value"]
    assert S.PUBLISHED_TEST_A["R3"]["ratio"] == pytest.approx(row["ratio"], abs=5e-4)


def test_a_statistic_pinned_in_every_repeat_is_not_reported_as_a_failure():
    """`hard_rank` reads 256 in every repeat of both arms.

    `sep > floor` is then `1.0 > 1.0`, which is False, and reporting that as a
    FAILURE would be wrong in the direction that flatters this paper's thesis:
    the statistic did not reject the comparison, it cannot rule on it at all.
    """
    v = S.verdict(1.0, 1.0, monotone_in_n=True)
    assert v["passes"] is None and "NO DISCRIMINATION" in v["why"]
    assert S.verdict(1.5, 1.2, monotone_in_n=True)["passes"] is True
    assert S.verdict(1.1, 1.2, monotone_in_n=True)["passes"] is False


def test_the_direction_of_the_evidence_is_labelled():
    assert "PASSING is informative" in S.verdict(1.5, 1.2, monotone_in_n=True)["informative"]
    assert "EITHER outcome" in S.verdict(1.5, 1.2, monotone_in_n=False)["informative"]


# --------------------------------------------------------------------------
# Test C, and the sign convention that would otherwise fake a disagreement
# --------------------------------------------------------------------------
def test_the_direction_table_agrees_with_the_selection_rule():
    """`alpha-ReQ |alpha - 1|` is smaller-is-better here because it is there.

    A separate direction table would be a second definition of "better", and the
    failure mode is specific and ugly: alpha-ReQ's arms are PERFECTLY ordered on
    this block, and grading them without the sign would print INVERTED -- a
    statistic reported as contradicting the others when it agrees with them.
    """
    from morpheus.v2.research.rebase.p2 import p2_selection_rule as SR
    by_name = {n: d for n, _, d in SR.METRICS}
    assert by_name["alpha_ReQ_raw |alpha-1|"] == -1
    assert S.DIRECTION["alpha_req_abs_dev"] == -1
    for name in ("rankme_raw (as published)", "participation_ratio_raw", "stable_rank_raw"):
        assert by_name[name] == +1
    for name in ("RankMe", "PR", "stable_rank", "R1", "R3"):
        assert S.DIRECTION[name] == +1
    # The raw exponent has no agreed direction anywhere on this project, so it
    # gets none here either rather than being given one.
    assert S.DIRECTION["alpha_req"] == 0
    assert set(S.DIRECTION) == set(S.STATS)


def test_complete_separation_respects_the_direction():
    high, low = [1.0, 1.1, 1.2], [2.0, 2.1, 2.2]
    assert S.ordering(high, low, +1)["separated"] is False
    assert S.ordering(high, low, +1)["inverted"] is True
    # Same numbers, smaller-is-better: now it is the high arm that wins.
    assert S.ordering(high, low, -1)["separated"] is True
    assert S.ordering(high, low, -1)["inverted"] is False
    assert S.ordering(high, low, 0)["separated"] is None


def test_complete_separation_p_is_the_exact_permutation_probability():
    """1/C(10,5) at five per arm, 1/C(20,10) at ten -- it STRENGTHENS with n.

    This is what makes Test C worth running beside Test B: Test B can only get
    harder to pass as repeats accumulate, and Test C can only get more surprising
    when it does pass. A result supported by both is supported for two reasons
    that do not share a failure mode.
    """
    five = S.ordering([10.0, 11.0, 12.0, 13.0, 14.0], [1.0, 2.0, 3.0, 4.0, 5.0])
    assert five["separated"] is True
    assert five["exact_one_sided_permutation_p"] == pytest.approx(1 / 252)
    ten = S.ordering([10.0 + i for i in range(10)], [1.0 + i * 0.1 for i in range(10)])
    assert ten["exact_one_sided_permutation_p"] == pytest.approx(1 / 184756)
    assert ten["exact_one_sided_permutation_p"] < five["exact_one_sided_permutation_p"]
    # And the scope note travels with it: same-seed repeats are not seed replicates.
    assert "NOT a p-value for the momentum effect" in five["scope"]


def test_a_constant_statistic_is_not_reported_as_an_overlap():
    """`hard_rank` reads 256 in every repeat of both arms.

    "Overlap" would read as a statistic that looked and disagreed. It did not
    look. §4.9's "16/16" instance is the same statistic pinned at a batch size.
    """
    c = S.ordering([256.0] * 5, [256.0] * 5)
    assert c["separated"] is None and c["constant"] is True


def test_every_hill_variant_orders_the_arms_even_where_the_ratio_test_fails(n5):
    """The finding Test B alone hides.

    Under EVERY variant of `RANK_VARIANTS` -- both Hill orders, centred and
    uncentred, row-normalised and not -- every m = 0.999 repeat is above every
    m = 0.99 repeat. The statistics do not disagree about the ORDER of the two
    arms; they disagree only about whether the gap is large enough relative to
    the within-arm spread. That is a much narrower disagreement than "R1 says yes
    and R3 says no", and it is the honest way to state it.
    """
    t = n5["tests"]
    for name in ("R1", "R2", "R3", "R1_uncentred", "R1_rownorm", "R2_uncentred", "RankMe"):
        c = t[name]["test_C_complete_separation"]
        assert c["separated"] is True, f"{name} does not order the arms"
        assert c["exact_one_sided_permutation_p"] == pytest.approx(1 / 252)
    # alpha-ReQ orders them too, once its sign convention is applied.
    assert t["alpha_req_abs_dev"]["test_C_complete_separation"]["separated"] is True
    # The published alternatives that are NOT effective ranks do not order them.
    for name in ("PR", "PR_rownorm", "stable_rank"):
        c = t[name]["test_C_complete_separation"]
        assert c["separated"] is False and c["inverted"] is False, f"{name} unexpectedly ordered"
    assert t["hard_rank"]["test_C_complete_separation"]["constant"] is True


# --------------------------------------------------------------------------
# The measured result, and the arithmetic the write-up leans on
# --------------------------------------------------------------------------
def test_the_n5_numbers_reproduce_the_floors_already_published(n5):
    """Two independent paths to the same two floors, as §4.1's test insists on.

    `p2_probe_floors.py` published R1 1.1547× and R3 1.1947× for this pair from
    the same ten states. This module reaches them through its own arm-fold path,
    and a disagreement would mean one of the two is wrong.
    """
    t = n5["tests"]
    assert t["R1"]["floor"] == pytest.approx(1.1547, abs=5e-4)
    assert t["R3"]["floor"] == pytest.approx(1.1947, abs=5e-4)
    assert t["R1"]["test_B_worst_case_over_repeats"]["separation"] == pytest.approx(1.4533, abs=5e-4)
    assert t["R3"]["test_B_worst_case_over_repeats"]["separation"] == pytest.approx(1.1384, abs=5e-4)
    assert t["R1"]["test_B_worst_case_over_repeats"]["passes"] is True
    assert t["R3"]["test_B_worst_case_over_repeats"]["passes"] is False
    assert t["R3"]["test_A_published_single_draw"]["passes"] is True


def test_the_split_is_the_hill_order_and_not_the_preprocessing(n5):
    """Every order-1 statistic passes Test B; every order-2 statistic fails it.

    This is the finding that makes "R1 passes, R3 fails" more than a clash of two
    labels: it survives switching centring on and off and row normalisation on
    and off, and it tracks exactly one thing -- the Hill order of the same
    spectrum. `stable_rank` and `PR` are the order-2 family's other members and
    they fail hardest.
    """
    t = n5["tests"]
    order1 = ["R1", "R1_rownorm", "R1_uncentred", "RankMe"]
    order2 = ["R2", "R3", "R2_uncentred", "PR", "PR_rownorm", "stable_rank"]
    for name in order1:
        assert t[name]["test_B_worst_case_over_repeats"]["passes"] is True, name
    for name in order2:
        assert t[name]["test_B_worst_case_over_repeats"]["passes"] is False, name
    # and the order-1/order-2 pairs that differ ONLY in the order disagree:
    for one, two in (("R1", "R2"), ("R1_rownorm", "R3"), ("R1_uncentred", "R2_uncentred")):
        assert (t[one]["test_B_worst_case_over_repeats"]["passes"]
                is not t[two]["test_B_worst_case_over_repeats"]["passes"]), (one, two)


def test_the_ratio_verdict_flips_with_the_reading_step_and_the_order_verdict_does_not(n5):
    """The sharpest thing these ten runs say, and it is not in the paper.

    §5.4 limit 2 is read at step 600. The probe states were exported at every
    step, so the SAME two arms and the SAME repeats can be given the SAME tests
    at 100, 200, 300, 400 and 500 for no GPU at all -- and the ratio test's
    verdict is not stable across them. Under R3 it reads
    FAIL, FAIL, PASS, FAIL, PASS, FAIL: it passes at **two of six** readings and
    **changes its answer four times in six consecutive steps**. Under RankMe it
    changes four times as well (PASS, FAIL, PASS, FAIL, FAIL, PASS). Only under
    canonical R1 does it settle -- one transition, and stable from step 300 on.

    Test C, meanwhile, is SEP at every step under every one of them. So the
    statistics never disagree about which arm is higher -- they disagree about
    whether the gap clears the noise, and that disagreement is as much with the
    reading step as with each other.
    """
    by_step = n5["by_step"]
    assert sorted(by_step, key=int) == ["100", "200", "300", "400", "500", "600"]

    def bs(stat):
        return [by_step[s][stat]["test_B_worst_case_over_repeats"]["passes"]
                for s in sorted(by_step, key=int)]

    assert bs("R3") == [False, False, True, False, True, False]
    assert bs("R1") == [False, False, True, True, True, True]
    assert bs("RankMe") == [True, False, True, False, False, True]
    # R3's ratio verdict is a coin flip over the reading step; R1's settles. The
    # count that matters is how often it CHANGES, not how often it passes.
    def flips(v):
        return sum(a is not b for a, b in zip(v[:-1], v[1:]))

    assert sum(bs("R3")) == 2 and flips(bs("R3")) == 4
    assert flips(bs("RankMe")) == 4
    assert flips(bs("R1")) == 1
    assert all(bs("R1")[2:]), "R1 does not settle from step 300"

    # And the ORDER verdict is unanimous at every step, for every Hill variant.
    for step in by_step:
        for stat in ("R1", "R2", "R3", "R1_uncentred", "R1_rownorm", "R2_uncentred", "RankMe"):
            c = by_step[step][stat]["test_C_complete_separation"]
            assert c["separated"] is True, f"{stat} fails to order the arms at step {step}"

    # Test A -- the audit's OWN rule -- would have read FAIL at step 200. The row
    # is read at 600 so this does not overturn it, but the pass is step-contingent
    # and the draft does not say so.
    a = {s: by_step[s]["R3"]["test_A_published_single_draw"]["passes"] for s in by_step}
    assert a["200"] is False and a["600"] is True


def test_at_ten_repeats_the_row_stops_clearing_its_own_floor(n10):
    """The predeclared primary falsifier fired. §5.4 limit 2 does not clear at n = 10.

    Doubling the repeats per arm moves the R3 floor from **1.195x to 1.326x**,
    past the row's own published ratio of **1.262x**, so `p2_floor_audit.check`'s
    rule -- the paper's rule, not a new one -- now returns FAIL for the one
    selection this paper makes at the hyperparameter the project ships.

    The whole increase is on the m = 0.99 arm, and it is one run: `m0.99 rep6`
    reads R1 5.520 / R3 4.465 where the other nine span 6.462-7.124 / 5.142-5.922.
    It is not excludable. It trained (`biology_contrastive` 7.575 at step 600
    against chance ln 80 = 4.382), its R1 is nowhere near the collapsed arm's
    ~1.5-2.7, and the project's own bimodality rule does not flag it
    (`concordant` is False, so `bimodal` is False). It is what same-seed GPU
    non-determinism does at m = 0.99 roughly one run in ten, and five repeats did
    not see it.
    """
    t = n10["tests"]
    assert t["R3"]["floor"] == pytest.approx(1.3263, abs=5e-4)
    assert t["R3"]["test_A_published_single_draw"]["passes"] is False
    assert t["R3"]["floor"] > S.PUBLISHED_TEST_A["R3"]["ratio"]
    # The floor changed hands: at n = 5 the m = 0.999 arm carried it, at n = 10 the
    # m = 0.99 arm does. The 2026-08-05 00:00 entry's "carried by m = 0.999" and
    # "the smallest of the sixteen probe floors" are both overtaken.
    assert t["R3"]["floor_arm"] == "m099"
    assert t["R1"]["floor_arm"] == "m099"
    # R1's ratio test still survives -- it was the one thing that could have
    # broken and did not -- but it is now alone in doing so.
    assert t["R1"]["test_B_worst_case_over_repeats"]["passes"] is True
    for name in ("R1_uncentred", "RankMe", "R3", "R2", "R2_uncentred"):
        assert t[name]["test_B_worst_case_over_repeats"]["passes"] is False, name


def test_the_independent_second_five_reaches_the_same_verdict(n10):
    """Not one unlucky draw: repeats 6-10 ALONE also put the floor above 1.262x.

    This is the check that decides whether the n = 10 floor is a real widening or
    an artifact of pooling two batches. Repeats 6-10 are an independent n = 5 --
    same arms, same seed, same workspace, same 10-way concurrency, run 2h46m
    after the first five. On their own they give an R3 floor of **1.279x**, which
    is also above the row's 1.262x. Both halves of the twenty runs agree; the
    ORIGINAL five were the favourable draw.

    And there is no batch effect to blame it on: at m = 0.99 repeats 6-10 span
    5.520-7.029 against repeats 1-5's 6.462-7.124, which OVERLAP, so the two
    batches are samples of one spread rather than two levels.
    """
    late = _load(LATE5)
    assert late["config"]["rep_offset"] == 5 and late["config"]["reps"] == 5
    assert late["tests"]["R3"]["floor"] == pytest.approx(1.279, abs=1e-3)
    assert late["tests"]["R3"]["test_A_published_single_draw"]["passes"] is False
    assert late["tests"]["R1"]["test_B_worst_case_over_repeats"]["passes"] is True

    early = {r: v["views"]["wsi_biology"]["R1"]
             for r, v in n10["reps"]["m099"].items() if int(r[3:]) <= 5}
    later = {r: v["views"]["wsi_biology"]["R1"]
             for r, v in n10["reps"]["m099"].items() if int(r[3:]) > 5}
    assert min(later.values()) < min(early.values())      # the later five reach lower
    assert max(later.values()) > min(early.values())      # ... and still overlap
    assert min(early.values()) < max(later.values())


def test_the_arms_stay_completely_separated_at_ten_repeats(n10):
    """What survives, and it survives more strongly than before.

    Test C is the one test here that gets MORE surprising with more repeats, and
    it holds: all ten m = 0.999 repeats above all ten m = 0.99 repeats, exact
    one-sided permutation p = 1/184,756. So the honest residue of this row is an
    ORDERING claim, not a ratio claim.
    """
    for name in ("R1", "R2", "R3", "R1_rownorm", "R1_uncentred", "RankMe"):
        c = n10["tests"][name]["test_C_complete_separation"]
        assert c["separated"] is True, name
        assert c["exact_one_sided_permutation_p"] == pytest.approx(1 / 184756)
    # It is not universal, and the exception is recorded rather than dropped:
    # the uncentred order-2 variant loses the separation at ten repeats.
    assert n10["tests"]["R2_uncentred"]["test_C_complete_separation"]["separated"] is False


def test_the_momentum_grid_is_a_smooth_monotone_trend_not_a_step():
    """The shipped comparison is two rungs of a ladder, not a discontinuity.

    Five same-seed repeats at each of m = 0, 0.98, 0.99, 0.995, 0.999 at step 600.
    Under canonical R1 the ranges are

        m0     1.564- 2.738
        m0.98  4.783- 6.008
        m0.99  6.462- 7.124
        m0.995 7.710- 8.249
        m0.999 10.353-11.955

    -- strictly increasing, with EVERY adjacent pair completely separated and no
    overlap anywhere on the grid. So the predeclared reading is **SMOOTH**: the
    two-point m = 0.999 vs m = 0.99 contrast is a coarse read of a monotone trend,
    and the honest claim is about the trend rather than about 0.999 beating its
    neighbour.

    And the individual rungs are small. Only ONE of the four adjacent increments
    clears its own floor under R1 (0.995 -> 0.999, 1.255x against 1.155x) and NONE
    clears under R3. The m = 0.999-over-m = 0.99 gap is two increments each of
    which is individually inside the floor.
    """
    grid = _load(GRID)
    assert grid["config"]["arms"] == ["0", "0.98", "0.99", "0.995", "0.999"]
    order = ["m0", "m098", "m099", "m0995", "m0999"]
    for stat in ("R1", "R3"):
        mins = [grid["arm_floor"][a][stat]["min"] for a in order]
        assert mins == sorted(mins), f"{stat} is not monotone in m across the grid"
    # Under R1 the grid never overlaps; under R3 m = 0.995 and m = 0.99 do.
    r1 = {a: grid["arm_floor"][a]["R1"] for a in order}
    for lo, hi in zip(order[:-1], order[1:]):
        assert r1[hi]["min"] > r1[lo]["max"], f"R1 overlaps between {lo} and {hi}"
    r3 = {a: grid["arm_floor"][a]["R3"] for a in order}
    assert r3["m0995"]["min"] < r3["m099"]["max"], "R3 no longer overlaps at 0.99/0.995"

    # The file is written with `sort_keys=True`, so compare as a set: the ORDER
    # of the four adjacent pairs is alphabetical on disk, not the grid's order.
    steps = grid["adjacent_steps"]
    assert set(steps) == {"m098_over_m0", "m099_over_m098",
                          "m0995_over_m099", "m0999_over_m0995"}
    passing = [k for k, v in steps.items()
               if v["R1"]["test_B_worst_case_over_repeats"]["passes"]]
    assert passing == ["m0999_over_m0995"], passing
    assert not any(v["R3"]["test_B_worst_case_over_repeats"]["passes"]
                   for v in steps.values()), "an R3 adjacent step unexpectedly clears"
    # Every adjacent pair is still completely ORDERED under R1 -- the trend is
    # real even though each rung is inside its own floor.
    for k, v in steps.items():
        assert v["R1"]["test_C_complete_separation"]["separated"] is True, k


def test_more_repeats_can_only_move_test_B_against_the_pass(n5, n10):
    """The arithmetic the honesty of this whole exercise rests on.

    The floor is `max/min` over an arm's repeats and so is non-decreasing in the
    repeat count; the Test B separation is `min(high)/max(low)` and so is
    non-increasing in it. Therefore a statistic that fails Test B at n = 10 has
    told us nothing it had not already told us at n = 5, and a statistic that
    PASSES at n = 10 has survived a test that could only have got harder. If this
    assertion ever fails, the two files are not nested samples of one experiment
    and no comparison between them is valid.
    """
    for name, e5 in n5["tests"].items():
        e10 = n10["tests"][name]
        if e5["floor"] is None or e10["floor"] is None:
            continue
        assert e10["floor"] >= e5["floor"] - 1e-9, f"{name}: floor shrank with more repeats"
        s5 = e5["test_B_worst_case_over_repeats"].get("separation")
        s10 = e10["test_B_worst_case_over_repeats"].get("separation")
        if s5 is None or s10 is None:      # `alpha_req` has no agreed direction
            continue
        assert s10 <= s5 + 1e-9, f"{name}: worst-case separation grew with more repeats"
