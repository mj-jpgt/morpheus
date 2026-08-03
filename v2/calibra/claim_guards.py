"""Claim admissibility guards — caveats as executable policy, not prose.

Every limitation encoded here has one property in common: **it produces good-looking
numbers while the claim is wrong.** None of them announce themselves. A prose caveat in
a handoff gets skimmed; a guard that refuses to mark a claim admissible does not.

The guards are deliberately narrow. They do NOT enforce completeness -- partial coverage
is normal in computational biology and is extended by adding scales. They enforce the
much smaller set of conditions under which a *causal attribution* silently becomes wrong.
That distinction matters: a predictive model may use a confounded feature freely, because
nobody asks why the prediction works. A certified causal catalogue may not, because
misattribution is its only real failure mode.

Usage::

    verdict = validate_claim({"kind": "legible_axis", ...})
    if not verdict.admissible:
        # emit the repo status convention, never a silent drop
        row = {"metric": "status", "value": float("nan"),
               "note": "inadmissible_" + verdict.blockers[0].code}
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path

__all__ = ["Blocker", "ClaimVerdict", "validate_claim", "CAVEATS", "CLAIM_EVIDENCE_PATH",
           "evidence_digest", "load_claim_evidence", "validate_recorded_claim"]

# The project's recorded claim state. Before this file existed, `validate_claim`
# read evidence from NOWHERE: nothing in production built a claim dict, and the
# record of E0's admissibility lived only as a hardcoded fixture in
# tests/test_claim_guards.py. A guard that can only restate a decision somebody
# already made is not a guard -- it is documentation that looks like a check.
CLAIM_EVIDENCE_PATH = "v2/research/rebase/nature/claim_evidence.json"

_COMMIT_RE = re.compile(r"^[0-9a-f]{7,40}$")


@dataclass(frozen=True)
class Blocker:
    code: str
    why: str
    mechanism: str          # HOW the claim goes wrong while the numbers look fine
    discharge: str          # what evidence would clear it


CAVEATS: dict[str, Blocker] = {
    "composition_attribution": Blocker(
        "composition_attribution",
        "A morphology-legibility claim needs a cell-of-origin / composition control.",
        "The intervention dictionary is built from pure cell-line populations, which contain no "
        "immune cells, stroma, vasculature or tissue architecture -- much of what an H&E slide "
        "actually shows. In patients, tumour-intrinsic programmes CORRELATE with composition "
        "(proliferative tumours skew immune-cold, mesenchymal ones stroma-rich), so a fit onto the "
        "dictionary still succeeds numerically while the coefficients quietly absorb the "
        "composition signal the dictionary cannot represent. The catalogue then reports 'gene g is "
        "morphologically legible' when the truth is 'g correlates with an infiltration pattern, and "
        "the pattern is what is visible'. This is misattribution, not missing coverage, and adding "
        "further modalities multiplies it rather than resolving it.",
        "Label the axis tumour-autonomous / immune-mediated / stromal / composition-driven using "
        "spatial or deconvolution evidence, OR beat a capacity-matched cell-composition baseline.",
    ),
    "purity_confound": Blocker(
        "purity_confound",
        "Tumour purity must be in the adjustment set before any morphology<->molecular claim.",
        "TCGA bulk RNA is a MIXTURE (roughly 30-90% tumour); the dictionary atoms are PURE "
        "populations. Fitting a mixture with pure atoms lets the coefficients absorb purity. "
        "Purity is simultaneously one of the most visually obvious features on a slide, so "
        "'morphology predicts molecular state' has a trivial explanation that reproduces the "
        "result exactly. Flagged as open since Phase 1 and still not closed.",
        "Include purity/mRNAsi (or an equivalent tumour-fraction estimate) in the confound design "
        "and show the effect survives it.",
    ),
    "sign_blind": Blocker(
        "sign_blind",
        "A directional claim cannot rest on a sign-blind statistic.",
        "Subspace overlap svdvals(Va^T Vb) is invariant to the sign of a response: an ANTI-aligned "
        "perturbation effect scores identically to an aligned one. Separately, CRISPRi measures "
        "loss of function, while tumours are substantially driven by gain of function "
        "(amplification, activating mutation). So the dictionary is one-directional AND the "
        "statistic cannot read direction -- a discovery claim can come out inverted.",
        "Use a signed statistic (e.g. signed canonical correlation on matched directions) and state "
        "the knockdown-vs-gain-of-function asymmetry explicitly.",
    ),
    "proliferation_deflation": Blocker(
        "proliferation_deflation",
        "A transfer claim must be shown not to reduce to proliferation.",
        "The responsive arm is selected on HAVING a detectable transcriptional effect, which "
        "enriches for essential, core-machinery, ribosome and cell-cycle genes. If the measured "
        "alignment is proliferation-in-cell-lines matching proliferation-in-tumours, the result is "
        "the most generic axis in cancer biology and deflates to near-nothing -- while looking "
        "identical to a real finding.",
        "Re-run with proliferation/cell-cycle programme regressed out, or with the responsive arm "
        "stratified by proliferation loading, and show the gap survives.",
    ),
    "single_platform": Blocker(
        "single_platform",
        "Replication across cell lines from one assay is not platform-independent replication.",
        "K562 and RPE1 are different lineages but the SAME Perturb-seq protocol (CRISPRi + scRNA-seq, "
        "same processing). A shared platform artifact replicates across them exactly as readily as "
        "shared biology. Cross-lineage agreement rules out 'a K562 quirk'; it does not rule out "
        "'a Perturb-seq quirk'.",
        "Replicate on a perturbation resource from a different platform, or scope the claim to "
        "'replicates across lineages within Perturb-seq'.",
    ),
    "no_external_cohort": Blocker(
        "no_external_cohort",
        "Certification requires at least one cohort outside the discovery cohort.",
        "Every morphology result so far is TCGA, which carries well-documented site and scanner "
        "effects. Confound removal was verified for cancer type, not across an independent cohort, "
        "so a site-specific artifact would survive the current checks intact.",
        "Replicate in an external cohort (CPTAC, HEST-1k, or equivalent).",
    ),
}

# Which blockers apply to which kind of claim.
_REQUIREMENTS: dict[str, tuple[str, ...]] = {
    "legible_axis": ("composition_attribution", "purity_confound", "no_external_cohort"),
    "gene_attribution": ("composition_attribution", "purity_confound", "sign_blind", "no_external_cohort"),
    "transfer": ("proliferation_deflation", "single_platform"),
    "direction": ("sign_blind",),
    "cross_platform": ("single_platform",),
}

# Evidence field that discharges each blocker, and the value that counts as discharged.
_DISCHARGED_BY: dict[str, str] = {
    "composition_attribution": "composition_control",
    "purity_confound": "purity_in_adjustment_set",
    "sign_blind": "statistic_is_signed",
    "proliferation_deflation": "proliferation_controlled",
    "single_platform": "platforms",
    "no_external_cohort": "external_cohorts",
}


@dataclass
class ClaimVerdict:
    admissible: bool
    blockers: list[Blocker] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def as_rows(self, method: str = "", state: str = "") -> list[dict]:
        """Repo-standard status rows -- an inadmissible claim is VISIBLE, never dropped."""
        return [{"method": method, "representation_state": state, "task": "claim_guard",
                 "target": b.code, "metric": "status", "value": float("nan"),
                 "note": f"inadmissible_{b.code}"} for b in self.blockers]


def _is_discharged(code: str, claim: dict) -> bool:
    field_name = _DISCHARGED_BY[code]
    value = claim.get(field_name)
    if code == "single_platform":
        # Distinct PLATFORMS, not distinct cell lines. Two lineages on one protocol is one platform.
        return len({str(v) for v in (value or [])}) >= 2
    if code == "no_external_cohort":
        return len(value or []) >= 1
    if code == "composition_attribution":
        # Either the axis is labelled by cell-of-origin, or it beat a composition baseline.
        if isinstance(value, dict):
            return bool(value.get("cell_of_origin_label")) or bool(value.get("beats_composition_baseline"))
        return False
    return value is True


def evidence_digest(evidence: dict) -> str:
    """Tamper-EVIDENCE digest over a claim's evidence block, excluding the digest.

    This is not tamper-proofing and must not be described as such: anyone editing
    the file can recompute it. What it buys is that a value cannot be changed
    *silently* -- the digest must be updated too, which is a visible act in a diff
    and is what review can catch.
    """
    payload = {k: v for k, v in evidence.items() if k != "sha256"}
    return sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _resolved_value(field_name: str, record: object, repo_root: Path) -> tuple[object, str]:
    """Extract an evidence value only if its provenance resolves.

    **Evidence with no provenance is not evidence.** A field set to ``true`` with
    no run, no notebook entry and no commit is exactly the failure mode this whole
    file exists to prevent -- somebody typed the answer they wanted. Such a field
    is reported as ABSENT so its blocker fires, rather than quietly discharging.
    """
    if not isinstance(record, dict) or "value" not in record:
        return None, f"{field_name}: not an evidence record with a 'value'"
    run, entry, commit = (str(record.get(k, "")).strip() for k in ("run", "entry", "commit"))
    if not run:
        return None, f"{field_name}: no run recorded"
    if not entry or not (repo_root / entry).exists():
        return None, f"{field_name}: notebook entry {entry!r} does not exist"
    if not _COMMIT_RE.match(commit):
        return None, f"{field_name}: commit {commit!r} is not a git hash"
    return record["value"], ""


def load_claim_evidence(path: str | Path, repo_root: str | Path) -> tuple[dict[str, dict], list[str]]:
    """Read the recorded claim state into the flat dicts ``validate_claim`` takes.

    Returns ``(claims, notes)``. Anything unreadable yields an EMPTY mapping and a
    note -- never a permissive default, matching the policy already applied to an
    unknown claim ``kind``.
    """
    repo_root = Path(repo_root)
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8"))
        records = document["claims"]
    except Exception as exc:
        return {}, [f"claim evidence unreadable at {path}: {exc!r}; every claim is inadmissible"]
    claims, notes = {}, []
    for name, record in records.items():
        evidence = record.get("evidence", {})
        expected, actual = str(evidence.get("sha256", "")), evidence_digest(evidence)
        if expected != actual:
            notes.append(f"{name}: evidence digest mismatch (recorded {expected[:12] or 'none'}, "
                         f"computed {actual[:12]}); treated as no evidence")
            claims[name] = {"kind": record.get("kind", "")}
            continue
        flat: dict[str, object] = {"kind": record.get("kind", "")}
        for field_name, item in evidence.items():
            if field_name == "sha256":
                continue
            value, why = _resolved_value(field_name, item, repo_root)
            if why:
                notes.append(f"{name}: {why}; treated as absent")
                continue
            flat[field_name] = value
        claims[name] = flat
    return claims, notes


def validate_recorded_claim(name: str, path: str | Path, repo_root: str | Path) -> ClaimVerdict:
    """Validate a claim from the recorded evidence file rather than a literal dict."""
    claims, notes = load_claim_evidence(path, repo_root)
    if name not in claims:
        kind_blockers = [CAVEATS[code] for code in _REQUIREMENTS.get("transfer", ())]
        return ClaimVerdict(False, kind_blockers,
                            [f"no recorded evidence for claim {name!r}", *notes])
    verdict = validate_claim(claims[name])
    relevant = [note for note in notes if note.startswith(f"{name}:")]
    return ClaimVerdict(verdict.admissible, verdict.blockers, [*verdict.notes, *relevant])


def validate_claim(claim: dict) -> ClaimVerdict:
    """Return whether a claim may be published, and precisely what is missing.

    Unknown ``kind`` is inadmissible by default: a claim shape nobody has thought about
    has not been checked, and defaulting to permissive is how unreviewed claims ship.
    """
    kind = str(claim.get("kind", ""))
    if kind not in _REQUIREMENTS:
        return ClaimVerdict(False, [], [f"unknown claim kind {kind!r}; add it to _REQUIREMENTS "
                                        f"(known: {sorted(_REQUIREMENTS)})"])
    blockers = [CAVEATS[code] for code in _REQUIREMENTS[kind] if not _is_discharged(code, claim)]
    return ClaimVerdict(not blockers, blockers)
