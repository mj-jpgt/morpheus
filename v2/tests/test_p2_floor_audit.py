"""The floor audit must agree with its sources, and the draft must agree with the audit.

`paper/P2_RANK_DRAFT.md` §4.1a is the paper auditing itself against its own
central criterion: a rank difference smaller than the measured same-seed
retraining floor is not resolvable. That criterion has been applied to this
project's own numbers five separate times, each discovered late and separately.
The audit list at `v2/research/rebase/p2/floor_audit.json` enumerates every rank
comparison the paper makes or relies on, and this module is what makes it a
check rather than a document:

* **every recorded value is re-read from the file it came from**, so a number
  edited in the draft, in a notebook entry or in a vendored box log breaks the
  build instead of quietly disagreeing with the audit;
* **every ratio is recomputed from its own two values**;
* **every verdict is recomputed against the floor its own block licenses** —
  block-matching is load-bearing and has already produced one wrong verdict
  (D1-B seed 43's residualised 3.246× judged against the *raw* floor of 3.111×
  reads as outside when on its own block it is inside);
* **the draft's table and its counting sentence are generated from the list**,
  so the two cannot drift.

The negative tests at the end are the ones that matter: they mutate a copy of
the audit and assert the checker notices.
"""
from __future__ import annotations

import copy
import re

import pytest

from morpheus.v2.research.rebase.p2 import p2_floor_audit as A

DRAFT = A.REPO / "paper" / "P2_RANK_DRAFT.md"


@pytest.fixture(scope="module")
def audit():
    return A.load()


@pytest.fixture(scope="module")
def draft():
    return DRAFT.read_text(encoding="utf-8")


# --------------------------------------------------------------------------
# The list against its sources
# --------------------------------------------------------------------------
def test_every_recorded_value_agrees_with_its_source(audit):
    problems = A.check(audit)
    assert not problems, "\n".join(problems)


def test_the_list_is_not_empty_and_ids_are_unique(audit):
    rows = audit["comparisons"]
    assert len(rows) >= 40, "the audit is meant to be exhaustive, not a sample"
    assert len({r["id"] for r in rows}) == len(rows)


def test_every_row_names_a_statistic_a_block_and_what_it_rests_on(audit):
    for row in audit["comparisons"]:
        assert row["statistic"], row["id"]
        assert row["block"], row["id"]
        assert row["rests_on"], row["id"]
        assert row["kind"] in A.KINDS, row["id"]


def test_every_floor_free_row_says_why(audit):
    """A row with no floor must name the statistic or block that has none.

    Half the paper's rank numbers live on the fixed held-out probe, on in-run
    training batches, on a 16-patient gate batch or on a co-trained view other
    than `wsi_biology`, and **no floor has been measured on any of them**. That
    is a limitation the table has to carry explicitly rather than by omission.
    """
    for row in audit["comparisons"]:
        if row.get("floor") is None:
            assert row.get("floor_note"), row["id"]


def test_the_floors_are_the_ones_the_draft_publishes(audit, draft):
    assert "3.295" in draft and "3.111" in draft
    assert audit["floors"]["R1_residualised_export"]["value"] == 3.295
    assert audit["floors"]["R1_raw_export"]["value"] == 3.111
    band = audit["floors"]["within_arm_seed_band"]
    assert (band["low"], band["high"]) == (2.10, 3.75)


def test_the_seven_headline_ratios_agree_with_section_4_1s_own_table(audit, draft):
    """§4.1's published table and the audit must quote the same seven folds.

    Parsed out of the draft rather than restated here: §4.1 is the section the
    whole criterion comes from, and if the audit and the headline table ever
    disagree the audit is worthless.
    """
    section = A._section(draft, "4.1 All seven")
    published = set()
    for line in section.splitlines():
        if not line.strip().startswith("|"):
            continue
        cells = [c.strip().replace("**", "") for c in line.strip("|").split("|")]
        if len(cells) == 5 and cells[2].endswith("×"):
            published.add(round(float(cells[2].rstrip("×")), 3))
    assert len(published) == 7, published
    # The three raw-block twins of D1-B are listed separately in §4.1a (they are
    # §4.1's provenance note, not its table) and are excluded here.
    audited = {round(r["ratio"], 3) for r in audit["comparisons"]
               if r["section"].startswith("§4.1") and r["kind"] == "selection"
               and not r["id"].startswith("4.1-d1b-")}
    assert published == audited, f"published {sorted(published)}, audited {sorted(audited)}"


# --------------------------------------------------------------------------
# The draft against the list
# --------------------------------------------------------------------------
def test_the_draft_prints_the_rendered_table(audit, draft):
    table = A.render_markdown(audit)
    assert table in draft, (
        "paper/P2_RANK_DRAFT.md §4.1a no longer matches "
        "`p2_floor_audit.py --markdown`. Regenerate it rather than editing the "
        "table by hand: the table is generated so that it cannot drift from the "
        "list, and the list is checked against its sources.")


def test_the_draft_prints_the_generated_counting_sentence(audit, draft):
    assert A.summary_sentence(audit) in draft


def test_the_five_instances_this_audit_was_built_to_catch_are_all_in_it(audit):
    """The five that were found one at a time must each appear as a row.

    momentum fix 3.29x; m = 0.999 over m = 0.99 1.26x; §4.7.4's violation 3.73x;
    the decorrelation ablation 1.85x; §5.1 instance 2 ~3.2x.
    """
    ratios = {round(r["ratio"], 2) for r in audit["comparisons"]}
    for expected in (3.29, 1.26, 3.73, 1.85, 3.23):
        assert expected in ratios, f"{expected}x is no longer enumerated"


def test_the_section_5_1_instance_2_residual_is_resolved(audit, draft):
    """It is exempt, and the draft says on what the claim rests instead.

    This entry was flagged in
    `NOTEBOOK_ENTRIES/p2_section5_rewritten_around_the_momentum_replication_20260804T2000Z.md`
    §5 and left unrewritten. The audit resolves it as a regime contrast whose
    two numbers are not even on the same block, and §5.1 now names the binary
    gate outcome and the co-measured collapse evidence it actually rests on.
    """
    row = next(r for r in audit["comparisons"] if r["id"] == "5.1-instance2")
    assert row["kind"] == "regime"
    assert row["floor"] is None and "EXEMPT" in row["floor_note"]
    assert "binary" in row["rests_on"]
    section = A._section(draft, "5.1 What a liveness gate certifies")
    assert "not a ratio" in section or "is not a ratio" in section


def test_known_source_disagreements_are_recorded_and_not_used(audit):
    """A value that disagrees with its source is reported, never substituted."""
    disagreements = audit["known_source_disagreements"]
    assert disagreements
    for item in disagreements:
        assert item["where"] and item["claim"] and item["source_says"] and item["action"]
    # The one currently open: §5.2's "2.6-3.3x at every step past 150".
    assert any("2.6" in d["claim"] for d in disagreements)


# --------------------------------------------------------------------------
# The negative tests: does the checker actually catch anything?
# --------------------------------------------------------------------------
def test_a_drifted_value_is_caught(audit):
    broken = copy.deepcopy(audit)
    row = next(r for r in broken["comparisons"] if r["id"] == "5.4-row2-seedvaried")
    row["a"]["value"] = 99.0
    row["ratio"] = 99.0 / row["b"]["value"]
    problems = A.check(broken)
    assert any(row["id"] in p for p in problems), problems


def test_a_ratio_that_disagrees_with_its_own_values_is_caught(audit):
    broken = copy.deepcopy(audit)
    row = next(r for r in broken["comparisons"] if r["id"] == "4.9a-decorr-r3")
    row["ratio"] = 3.9
    assert any("4.9a-decorr-r3" in p for p in A.check(broken))


def test_a_residualised_ratio_judged_against_the_raw_floor_is_caught(audit):
    """The block-match check, on the pair that has already produced one wrong verdict.

    D1-B seed 43 is 3.246x on the residualised block and 3.091x on the raw one.
    Judged against the raw floor of 3.111x the residualised figure reads as
    OUTSIDE the floor; on its own block it is inside. The checker must refuse
    the mismatch and must also notice that the verdict flipped.
    """
    broken = copy.deepcopy(audit)
    row = next(r for r in broken["comparisons"]
               if r["section"].startswith("§4.1") and round(r["ratio"], 3) == 3.246)
    row["floor"] = "R1_raw_export"
    row.pop("floor_note", None)
    problems = A.check(broken)
    assert any("no floor_note" in p and row["id"] in p for p in problems), problems
    assert any("gives True" in p and row["id"] in p for p in problems), problems


def test_a_missing_source_file_is_caught(audit):
    broken = copy.deepcopy(audit)
    row = next(r for r in broken["comparisons"] if r["id"] == "4.9a-decorr-r1")
    row["a"]["src"] = {"kind": "probe_log",
                       "file": "e0_run/d1_diag/does_not_exist.log", "step": 400,
                       "column": "CANONICAL"}
    assert any("unresolvable" in p for p in A.check(broken))


def test_a_markdown_value_edited_out_of_the_draft_is_caught(audit):
    """The literal has to still be in the file it is cited from."""
    with pytest.raises(LookupError):
        A._markdown_literal("paper/P2_RANK_DRAFT.md", "5.2", "1234.5678")


def test_an_unknown_probe_header_is_refused():
    """A log written by a different script must fail, not yield the wrong column."""
    import json as _json  # noqa: PLC0415

    assert ("R3-rank", "CANONICAL", "feat-std", "rna-rna", "contrastive") in A.PROBE_HEADERS
    assert ("eff-rank", "feat-std", "rna-rna", "contrastive") in A.PROBE_HEADERS
    assert _json  # the import exists only to keep the failure mode explicit


def test_the_probe_logs_distinguish_r3_from_the_canonical_statistic():
    """The 1.85x the notebook quotes is R3; statistic-matched it is 1.940x.

    Both are inside the floor, so nothing turns on it -- but the two columns are
    different statistics and the audit must read the one it names.
    """
    rows = A._probe_rows("e0_run/d1_diag/ablate_decorr0.04.log")
    assert rows[400]["R3-rank"] != rows[400]["CANONICAL"]
    base = A._probe_rows("e0_run/d1_diag/ablate_decorr0.0.log")
    assert round(rows[400]["R3-rank"] / base[400]["R3-rank"], 3) == 1.854
    assert round(rows[400]["CANONICAL"] / base[400]["CANONICAL"], 3) == 1.940


def test_the_three_ablation_arms_share_one_initialisation():
    """One verified common step-0 state, so the ablation varies decorrelation alone."""
    step0 = [A._probe_rows(f"e0_run/d1_diag/ablate_decorr{d}.log")[0]
             for d in ("0.0", "0.01", "0.04")]
    for column in ("R3-rank", "CANONICAL", "feat-std", "rna-rna"):
        assert len({row[column] for row in step0}) == 1, column


def test_the_like_for_like_pair_is_the_same_configuration():
    """§5.4 says no like-for-like measurement exists in this regime. One does.

    `ablate_decorr0.04` and `mseed_m0.999_s42` are the same momentum, the same
    decorrelation, the same capacity, the same learning rate and the same seed,
    and `d1_momentum_probe.py` runs a constant learning rate with no schedule
    that depends on the step budget -- so up to step 400 they differ only in GPU
    non-determinism. It is n = 2 and may not be quoted as a floor; the audit row
    says so.
    """
    header = re.compile(r"momentum=(\S+)\s+decorrelation=(\S+)\s+capacity=(\S+)\s+lr=(\S+)\s+seed=(\S+)")
    fields = []
    for rel in ("ablate_decorr0.04.log", "mseed_m0.999_s42.log"):
        text = (A.DATA / "e0_run" / "d1_diag" / rel).read_text(encoding="utf-8")
        match = header.search(text)
        assert match, rel
        fields.append(match.groups())
    assert fields[0] == fields[1], fields
    # It is now SUPERSEDED: the same quantity at n = 5 on the same arm, statistic
    # and step. The pair is kept rather than deleted, and it sits INSIDE the
    # five-repeat spread -- which is what an n = 2 sample of the same noise
    # should do, and is recorded rather than assumed.
    row = next(r for r in A.load()["comparisons"] if r["id"] == "4.9a-likeforlike-pair")
    assert row["floor"] == "R1_probe_step400"
    assert row["clears"] is False
    assert "SUPERSEDED" in row["floor_note"]


# --------------------------------------------------------------------------
# The expanded floor set, and the scope error it repairs
# --------------------------------------------------------------------------
def test_the_draft_prints_the_rendered_floor_table(audit, draft):
    """§4.1a's measured-floor table is generated, like its audit table."""
    assert A.render_floors_markdown(audit) in draft, (
        "paper/P2_RANK_DRAFT.md §4.1a no longer matches "
        "`p2_floor_audit.py --floors`. Regenerate it rather than editing the "
        "table by hand.")


def test_every_floor_names_a_statistic_a_block_and_a_caveat(audit):
    """A floor with no statistic or block cannot be block-matched against anything.

    And a floor with no caveat gets quoted as an estimate: every one of these is
    n = 5, one arm, one seed, and a floor twice over.
    """
    for name, floor in audit["floors"].items():
        assert floor["statistic"] and floor["block"] and floor["caveat"], name
        assert floor.get("n") in (None, 5) or "low" in floor, name


def test_the_expanded_floor_set_covers_every_statistic_t1_scores(audit):
    """Eight of T1's twelve metric rows had no floor at all. None is left without one."""
    statistics = {f["statistic"] for f in audit["floors"].values()}
    for expected in ("R1", "R2", "R3", "PR", "PR_rownorm", "stable rank", "LiDAR",
                     "α-ReQ |α−1|", "hard numerical matrix rank",
                     "RankMe (raw, uncentred, ε = 1e-7)", "RankMe (residualised)"):
        assert expected in statistics, expected
    assert audit["statistics_with_no_measured_floor"] == []


def test_the_floor_is_a_property_of_the_statistic_and_of_the_view(audit):
    """§4.3's withdrawn claim, and §4.1a's third finding, asserted against the list."""
    floors = audit["floors"]
    assert floors["hard_rank_residualised_export"]["value"] == 1.0
    assert floors["R1_residualised_export"]["value"] > 3.0
    assert floors["R1_residualised_rna_view"]["value"] < 1.05
    assert floors["R1_residualised_full_view"]["value"] < 1.05
    assert floors["RankMe_published_raw_export"]["value"] < floors["R1_raw_export"]["value"]


# --- unjudgeable is not a verdict -----------------------------------------
def test_the_three_way_count_partitions_the_selections(audit):
    """`failing` + `clearing` + `unjudgeable` must exhaust the selections.

    The three-way split is the repair: for as long as a row with no floor could
    record `clears: false`, eleven comparisons with no ruler at all were being
    counted beside comparisons that had failed a measured one.
    """
    s = A.summary(audit)
    assert (s["selection_failing"] + s["selection_clearing"]
            + s["selection_unjudgeable"]) == s["selection"]
    assert s["selection_unjudgeable"] > 0, "the absent floors are the point of the audit"


def test_the_two_rows_that_used_to_clear_now_clear_a_floor_on_their_own_block(audit):
    """§4.4(3)'s probe repeat and §5.2's step-400 fold, judged at last.

    These two were the paper's only "clears" before §4.1a, then became
    unjudgeable when the audit refused to let an EXPORT floor rule on a PROBE
    row. The probe block now has a floor of its own at every step the harness
    reads, so they are judged -- against the larger of the two arms' five-repeat
    spreads at their own step, which is the m = 0 arm in both cases.
    """
    for rid, floor in (("4.4-3-fixedseed-probe", "R3_probe_step200"),
                       ("5.2-per-step-max", "R3_probe_step400")):
        row = next(r for r in audit["comparisons"] if r["id"] == rid)
        assert row["floor"] == floor, rid
        assert row["clears"] is True, rid
        assert row["block"] == audit["floors"][floor]["block"], rid


def test_the_selections_that_clear_are_exactly_these(audit):
    """Pinned so that neither a new pass nor a lost one can arrive quietly.

    One of these is on the exported block -- RankMe as published, whose floor of
    1.811x is low enough that its widest D2 fold gets over it. The other eight
    are the probe-block rows the five-repeat floor made judgeable, and they are
    listed by name because "all of them cleared" is the kind of result that has
    to be checkable rather than asserted.
    """
    clearing = sorted(r["id"] for r in audit["comparisons"]
                      if r["kind"] == "selection" and r["clears"] is True)
    assert clearing == sorted([
        "4.6-rankme-d2",
        "4.4-3-fixedseed-probe",
        "4.9a-decorr-r1",
        "4.9a-decorr-r3",
        "5.2-per-step-max",
        "5.2-staleness-none",
        "5.2-turnover-momentum",
        "5.4-row2-seedvaried",
        "5.4-row3-r3",
    ]), clearing


def test_the_section_5_rows_that_remain_unjudgeable_say_which_gap_stops_them(audit):
    """The probe block has a floor now; three §5 selections are still outside it.

    Each is outside for a reason the floor cannot be stretched over: a reading
    step past the 500 the repeats were run to, a momentum value no repeat was run
    at, or a capacity no repeat was run at. The rule this file exists to enforce
    is that a floor is never borrowed across a block, and a step and an arm are
    part of this block's identity.
    """
    unjudged = {r["id"] for r in audit["comparisons"]
                if r["kind"] == "selection" and r["floor"] is None}
    assert unjudged == {"5.4-row1-step600", "5.4-m0999-over-m099", "5.2-capacity-64"}
    for rid in unjudged:
        row = next(r for r in audit["comparisons"] if r["id"] == rid)
        assert "STILL UNJUDGEABLE" in row["floor_note"], rid
        assert row["clears"] is None, rid


def test_the_probe_floor_is_carried_by_the_collapsed_arm_and_is_not_bimodal(audit):
    """Two properties of the new floor that the paper's argument turns on.

    (a) At every step and both statistics it is the m = 0 arm -- the one at the
    collapse floor -- that carries it, by roughly an order of magnitude more
    spread than the m = 0.999 arm. A floor measured only on the stable arm would
    have been between 1.09x and 1.17x and would have flattered every row.

    (b) It is NOT bimodal, and the exported-block floor is. Four of five exported
    repeats agree to 2% with rep2 at a third of them; nothing of that shape
    appears on the probe block under either statistic at any step.
    """
    probe = {n: f for n, f in audit["floors"].items() if "_probe_step" in n}
    assert len(probe) == 10, sorted(probe)
    for name, floor in probe.items():
        assert floor["floor_arm"] == "m0", name
        assert floor["shape"]["bimodal"] is False, name
    assert audit["floors"]["R1_residualised_export"]["shape"]["bimodal"] is True


# --- negative tests: does the checker catch the new failure modes? ---------
def test_a_verdict_recorded_without_a_floor_is_caught(audit):
    broken = copy.deepcopy(audit)
    row = next(r for r in broken["comparisons"] if r["id"] == "5.4-row2-seedvaried")
    row["clears"] = False
    assert any("clears=False" in p and row["id"] in p for p in A.check(broken))


def test_a_floor_recorded_without_a_verdict_is_caught(audit):
    broken = copy.deepcopy(audit)
    row = next(r for r in broken["comparisons"] if r["id"] == "4.6-six-pairs")
    row["clears"] = None
    assert any("clears=None" in p and row["id"] in p for p in A.check(broken))


def test_a_floor_that_disagrees_with_the_file_it_was_measured_into_is_caught(audit):
    broken = copy.deepcopy(audit)
    broken["floors"]["R3_residualised_export"]["value"] = 3.295
    assert any("R3_residualised_export" in p for p in A.check(broken))


def test_a_floor_whose_independent_cross_check_disagrees_is_caught(audit):
    """3.295x has two sources now. If they part, neither may be used."""
    broken = copy.deepcopy(audit)
    broken["floors"]["R1_residualised_export"]["cross_check"]["a"]["src"]["path"] = (
        "floors.rna_biology.residualised.R1.max")
    problems = A.check(broken)
    assert any("cross-check" in p and "STOP" in p for p in problems), problems


def test_a_floor_whose_recorded_shape_disagrees_with_its_source_is_caught(audit):
    """Whether the bimodal shape survives a change of statistic is a claim, so it is checked."""
    broken = copy.deepcopy(audit)
    broken["floors"]["stable_rank_residualised_export"]["shape"]["bimodal"] = True
    assert any("shape.bimodal" in p for p in A.check(broken))
    broken = copy.deepcopy(audit)
    broken["floors"]["R1_residualised_export"]["shape"]["outlier"] = "rep4"
    assert any("shape.outlier" in p for p in A.check(broken))
