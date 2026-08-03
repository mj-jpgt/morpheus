"""The caveats, as tests. If these pass vacuously the guards are theatre.

Each caveat encoded here shares one property: **it yields good numbers while the claim is
wrong.** These tests exist so a claim cannot be published without the specific control
that separates the real finding from its look-alike.

Deliberately NOT tested: completeness. Partial coverage of biology is normal and is
extended by adding scales; it is not a defect. What is guarded is the narrower case where
a *causal attribution* silently becomes wrong.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from morpheus.v2.calibra.claim_guards import (CAVEATS, CLAIM_EVIDENCE_PATH, evidence_digest,
                                              load_claim_evidence, validate_claim,
                                              validate_recorded_claim)

_REPO_ROOT = Path(__file__).resolve().parents[1]
_EVIDENCE = _REPO_ROOT / CLAIM_EVIDENCE_PATH


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


def test_the_recorded_e0_claim_is_still_inadmissible():
    """The project's REAL claim state, read from the evidence file rather than pinned here.

    This replaces a hardcoded fixture that encoded one verdict. The guard could not
    catch anything then: nothing in production built a claim dict, so 'discharging a
    blocker' meant editing this test to say so. Now the state is data with provenance
    and the guard reads it.

    E0 remains inadmissible because `single_platform` stands -- K562 and RPE1 are two
    lineages on ONE Perturb-seq protocol. If someone adds a genuinely different
    platform this fails and must be updated deliberately."""
    verdict = validate_recorded_claim("E0_basis_transfer_K562", _EVIDENCE, _REPO_ROOT)
    assert not verdict.admissible, "the guard stopped biting; this refactor must not loosen it"
    assert "single_platform" in {b.code for b in verdict.blockers}


def test_missing_or_unreadable_evidence_is_inadmissible_never_permissive(tmp_path):
    """The same default that already governs an unknown claim kind. An absent record
    must not read as 'nothing blocking'."""
    verdict = validate_recorded_claim("E0_basis_transfer_K562", tmp_path / "nope.json", tmp_path)
    assert not verdict.admissible and verdict.blockers
    claims, notes = load_claim_evidence(tmp_path / "nope.json", tmp_path)
    assert claims == {} and notes and "unreadable" in notes[0]


def test_a_silently_edited_value_is_caught_by_the_digest(tmp_path):
    """Not tamper-proofing -- tamper EVIDENCE. Changing a value without recomputing
    the digest voids the whole record, so the change cannot be silent."""
    entry = tmp_path / "entry.md"; entry.write_text("x")
    evidence = {"proliferation_controlled": {"value": True, "run": "r", "entry": "entry.md", "commit": "abc1234"}}
    evidence["sha256"] = evidence_digest(evidence)
    path = tmp_path / "e.json"
    document = {"claims": {"C": {"kind": "transfer", "evidence": evidence}}}
    path.write_text(json.dumps(document))
    claims, notes = load_claim_evidence(path, tmp_path)
    assert claims["C"]["proliferation_controlled"] is True and not notes

    document["claims"]["C"]["evidence"]["proliferation_controlled"]["value"] = "tampered"
    path.write_text(json.dumps(document))
    claims, notes = load_claim_evidence(path, tmp_path)
    assert "proliferation_controlled" not in claims["C"]
    assert any("digest mismatch" in note for note in notes)


@pytest.mark.parametrize("broken,reason", [
    ({"value": True, "entry": "entry.md", "commit": "abc1234"}, "no run"),
    ({"value": True, "run": "r", "commit": "abc1234"}, "does not exist"),
    ({"value": True, "run": "r", "entry": "absent.md", "commit": "abc1234"}, "does not exist"),
    ({"value": True, "run": "r", "entry": "entry.md"}, "not a git hash"),
    ({"value": True, "run": "r", "entry": "entry.md", "commit": "not-a-hash"}, "not a git hash"),
    (True, "not an evidence record"),
])
def test_evidence_without_resolvable_provenance_is_treated_as_absent(tmp_path, broken, reason):
    """'Somebody typed True' must not discharge anything. This is the whole point."""
    (tmp_path / "entry.md").write_text("x")
    evidence = {"proliferation_controlled": broken}
    evidence["sha256"] = evidence_digest(evidence)
    path = tmp_path / "e.json"
    path.write_text(json.dumps({"claims": {"C": {"kind": "transfer", "evidence": evidence}}}))
    claims, notes = load_claim_evidence(path, tmp_path)
    assert "proliferation_controlled" not in claims["C"]
    assert any(reason in note for note in notes), notes
    assert "proliferation_deflation" in {b.code for b in validate_claim(claims["C"]).blockers}


def test_REFACTOR_FALSIFIER_a_claim_cannot_be_discharged_by_editing_values_alone(tmp_path):
    """The falsifier pre-declared for this refactor itself.

    If a claim can be made admissible by writing the desired answers into the
    evidence file with no analysis and no provenance, the evidence file is just the
    old hardcoded fixture in a new location and the refactor has FAILED.
    """
    evidence = {
        "proliferation_controlled": {"value": True},
        "platforms": {"value": ["a", "b"]},
        "purity_in_adjustment_set": {"value": True},
        "composition_control": {"value": {"beats_composition_baseline": True}},
        "statistic_is_signed": {"value": True},
        "external_cohorts": {"value": ["CPTAC"]},
    }
    evidence["sha256"] = evidence_digest(evidence)          # digest honestly recomputed
    path = tmp_path / "e.json"
    path.write_text(json.dumps({"claims": {"C": {"kind": "transfer", "evidence": evidence}}}))
    verdict = validate_recorded_claim("C", path, tmp_path)
    assert not verdict.admissible, "REFACTOR FAILED: values alone discharged a claim"
    assert {b.code for b in verdict.blockers} == {"proliferation_deflation", "single_platform"}


def test_the_shipped_evidence_file_is_internally_consistent():
    """Every record's digest matches and every claimed entry file exists, or CI says so."""
    claims, notes = load_claim_evidence(_EVIDENCE, _REPO_ROOT)
    assert claims, "the shipped evidence file produced no claims"
    assert not notes, notes


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
