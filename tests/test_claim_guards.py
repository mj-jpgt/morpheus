"""The caveats, as tests. If these pass vacuously the guards are theatre.

Each caveat encoded here shares one property: **it yields good numbers while the claim is
wrong.** These tests exist so a claim cannot be published without the specific control
that separates the real finding from its look-alike.

Deliberately NOT tested: completeness. Partial coverage of biology is normal and is
extended by adding scales; it is not a defect. What is guarded is the narrower case where
a *causal attribution* silently becomes wrong.
"""
from __future__ import annotations

import pytest

from morpheus.v2.calibra.claim_guards import CAVEATS, validate_claim


def _fully_evidenced(kind: str) -> dict:
    """A claim carrying every discharge. Guards must let this through, or they are noise."""
    return {
        "kind": kind,
        "composition_control": {"cell_of_origin_label": "tumour_autonomous",
                                "beats_composition_baseline": True},
        "purity_in_adjustment_set": True,
        "statistic_is_signed": True,
        "proliferation_controlled": True,
        "platforms": ["perturb_seq_crispri", "orfeome_overexpression"],
        "external_cohorts": ["CPTAC"],
    }


@pytest.mark.parametrize("kind", ["legible_axis", "gene_attribution", "transfer", "direction", "cross_platform"])
def test_fully_evidenced_claims_are_admissible(kind):
    """Guards that block everything are as useless as guards that block nothing."""
    verdict = validate_claim(_fully_evidenced(kind))
    assert verdict.admissible, [b.code for b in verdict.blockers]


def test_composition_blocks_a_legibility_claim():
    """THE caveat. The dictionary is built from pure cell lines -- no immune cells, stroma
    or architecture -- yet those dominate an H&E slide. Because tumour-intrinsic programmes
    correlate with composition in patients, the fit succeeds numerically while the
    coefficients absorb composition, and the catalogue misattributes the morphology to the
    gene. Misattribution, not missing coverage."""
    claim = _fully_evidenced("legible_axis")
    claim["composition_control"] = None
    verdict = validate_claim(claim)
    assert not verdict.admissible
    assert "composition_attribution" in [b.code for b in verdict.blockers]

    # An unlabelled control that failed its baseline does not discharge it either.
    claim["composition_control"] = {"cell_of_origin_label": "", "beats_composition_baseline": False}
    assert not validate_claim(claim).admissible

    # Either route alone is enough.
    claim["composition_control"] = {"cell_of_origin_label": "immune_mediated"}
    assert validate_claim(claim).admissible
    claim["composition_control"] = {"beats_composition_baseline": True}
    assert validate_claim(claim).admissible


def test_purity_blocks_morphology_molecular_claims():
    """Bulk RNA is a 30-90% tumour MIXTURE; dictionary atoms are pure populations, so the
    coefficients absorb purity -- which is also plainly visible on a slide. Open since
    Phase 1 and never closed."""
    for kind in ("legible_axis", "gene_attribution"):
        claim = _fully_evidenced(kind)
        claim["purity_in_adjustment_set"] = False
        verdict = validate_claim(claim)
        assert not verdict.admissible
        assert "purity_confound" in [b.code for b in verdict.blockers]


def test_sign_blind_statistic_blocks_a_directional_claim():
    """svdvals(Va^T Vb) is invariant to response sign, so an ANTI-aligned effect scores
    identically to an aligned one -- and CRISPRi only measures loss of function while
    tumours are largely gain of function. A direction claim can come out inverted."""
    claim = _fully_evidenced("direction")
    claim["statistic_is_signed"] = False
    verdict = validate_claim(claim)
    assert not verdict.admissible
    assert [b.code for b in verdict.blockers] == ["sign_blind"]


def test_proliferation_blocks_a_transfer_claim():
    """The responsive arm is selected on HAVING an effect, enriching for essential /
    ribosome / cell-cycle genes. Proliferation-matching-proliferation looks identical to a
    real transfer result and deflates it to the most generic axis in cancer biology."""
    claim = _fully_evidenced("transfer")
    claim["proliferation_controlled"] = False
    verdict = validate_claim(claim)
    assert not verdict.admissible
    assert "proliferation_deflation" in [b.code for b in verdict.blockers]


def test_two_cell_lines_on_one_protocol_is_not_platform_replication():
    """K562 and RPE1 are different lineages on the SAME Perturb-seq protocol. A shared
    platform artifact replicates across them exactly as readily as shared biology. This is
    the caveat most likely to be quietly overstated in a write-up."""
    claim = _fully_evidenced("transfer")
    claim["platforms"] = ["perturb_seq_crispri", "perturb_seq_crispri"]   # K562 + RPE1
    verdict = validate_claim(claim)
    assert not verdict.admissible, "same protocol twice must not count as two platforms"
    assert "single_platform" in [b.code for b in verdict.blockers]

    claim["platforms"] = ["perturb_seq_crispri"]
    assert not validate_claim(claim).admissible


def test_external_cohort_required_for_certification():
    """Every morphology result so far is TCGA, whose site/scanner effects are well
    documented. Confound removal was verified for cancer type, not across cohorts."""
    claim = _fully_evidenced("legible_axis")
    claim["external_cohorts"] = []
    verdict = validate_claim(claim)
    assert not verdict.admissible
    assert "no_external_cohort" in [b.code for b in verdict.blockers]


def test_unknown_claim_kind_is_inadmissible_by_default():
    """A claim shape nobody has thought about has not been checked. Defaulting to
    permissive is exactly how an unreviewed claim ships."""
    verdict = validate_claim({"kind": "vibes"})
    assert not verdict.admissible
    assert verdict.notes and "unknown claim kind" in verdict.notes[0]


def test_current_e0_result_is_not_yet_an_admissible_transfer_claim():
    """Pin the ACTUAL state of the project as of the E0 run (commit 24d1bff).

    E0 returned 'supported' for K562 at ~10% of ceiling. That is a real result, but it is
    not yet a publishable transfer claim: proliferation was never controlled, and RPE1 --
    even once decidable -- is the same platform. If someone discharges these, this test
    fails and must be updated deliberately rather than drifting."""
    e0 = {"kind": "transfer", "proliferation_controlled": False,
          "platforms": ["perturb_seq_crispri", "perturb_seq_crispri"]}
    verdict = validate_claim(e0)
    assert not verdict.admissible
    assert {b.code for b in verdict.blockers} == {"proliferation_deflation", "single_platform"}


def test_every_caveat_explains_its_silent_failure_mode():
    """A guard whose rationale is not written down gets deleted by the next person who
    finds it inconvenient. Each caveat must say HOW it fails quietly and WHAT clears it."""
    for code, blocker in CAVEATS.items():
        assert blocker.code == code
        assert len(blocker.mechanism) > 120, f"{code}: mechanism must explain the silent failure"
        assert len(blocker.discharge) > 40, f"{code}: must say what discharges it"


def test_inadmissible_claims_emit_visible_status_rows():
    """Never drop a blocked claim silently -- it must appear in the task rows."""
    verdict = validate_claim({"kind": "transfer", "proliferation_controlled": False,
                              "platforms": ["perturb_seq_crispri"]})
    rows = verdict.as_rows(method="e0", state="k562")
    assert rows and all(r["metric"] == "status" for r in rows)
    assert all(r["note"].startswith("inadmissible_") for r in rows)
