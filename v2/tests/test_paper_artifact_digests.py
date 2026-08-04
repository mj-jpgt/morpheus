"""A quoted number must be traceable to an artifact by **content hash**, not by filename.

`test_paper_paths_resolve.py` checks that a cited path *exists*. That is a weaker property than it
looks, and the gap is not hypothetical. P1 §4.2's headline table cites its output directory and its
method, and names no artifact at all. Three distinct files called `d2_h_seed42.npz` sit on persistent
storage, written from three different commits, and the same estimator returns **0.3633**, **0.1782**
and **0.3785** from them. Only one of the three reproduces the published table. A filename resolved;
the number was still unattributable.

So this file checks the next thing along: wherever a draft quotes a number that came out of an
artifact, the draft must identify **which bytes**. Concretely:

* every artifact basename appearing in a draft must have a SHA-256 near it (`test_every_artifact_
  basename_in_a_draft_is_accompanied_by_a_hash`);
* every pinned number must appear in its draft together with the full digest of the artifact that
  produced it (`test_every_pinned_number_is_hash_identified_in_its_paper`);
* a number we **cannot** attribute must be marked as such rather than quietly attributed to whichever
  artifact happens to match (`test_an_unidentified_number_is_declared_unidentified`). This is the rule
  that keeps the registry honest: the failure mode being guarded against is not an absent hash, it is a
  *confident wrong* hash.

The digests below were taken on the box and are recorded here so the check runs in a checkout, where
the artifacts themselves are not present. `test_recorded_digests_match_the_box_when_reachable`
re-verifies them against the real files whenever this suite happens to run somewhere they are mounted.
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Root of the persistent NFS mount the artifacts live on. Absent in a checkout, present on the box.
BOX_ROOT = Path("/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d")

#: Artifacts that published numbers were computed from, by content hash.
#:
#: `d2_h_seed42.npz` deliberately appears three times. Two of the three are NOT the published artifact
#: and are pinned precisely so that a future edit citing them is caught rather than accepted.
ARTIFACTS = {
    "d2_final::d2_h_seed42": {
        "path": "runs/d2_final/artifacts/d2_h_seed42.npz",
        "sha256": "4a18b94f1017b85dd576f30ee8e3caf92d7897630a7054efb70166191cbe69e3",
        "note": "the artifact P1 4.2 was computed from; named by p1_evidence/run_track1.sh:6",
    },
    "d2_final::d2_i_seed42": {
        "path": "runs/d2_final/artifacts/d2_i_seed42.npz",
        "sha256": "028e8635465dd3c6d3dbead25a8c204ca1ae0cee4aabb20e5412847fb147b665",
        "note": "the d2_i rows of P1 4.2, same run",
    },
    "d2_v3_recovered::d2_h_seed42": {
        "path": "e0_run/d2_v3/recovered_artifacts/d2_h_seed42.npz",
        "sha256": "053490d685bf0dc47f2094831048db2bb884fe99f7ada3b57508ca23b561b899",
        "note": "August re-run, NOT published; raw joint LDA 0.1782",
    },
    "d2_v3_s42::d2_h_seed42": {
        "path": "e0_run/d2_v3/d2_v3_s42/artifacts/d2_h_seed42.npz",
        "sha256": "e81f4496f82c503a0dd1833e77cde2ea383cf79b0a6a7423a11c977c7f8f2625",
        "note": "August re-run, NOT published; raw joint LDA 0.3785",
    },
    "d2_v3_s42::d2_i_seed42": {
        "path": "e0_run/d2_v3/d2_v3_s42/artifacts/d2_i_seed42.npz",
        "sha256": "b49dc3efaf0a25dd9aacbd4792b843bc8cefff3840f15c64002b9ca45df8c4e9",
        "note": "August re-run, NOT published",
    },
}

#: Numbers a draft quotes, and the artifact each came from.
#:
#: ``artifact`` is a key of ``ARTIFACTS``, or ``None`` when the number cannot be attributed. ``None``
#: is not a placeholder to be filled in later with a guess — it is a claim in its own right, and
#: ``test_an_unidentified_number_is_declared_unidentified`` requires the draft to say so in the text.
PINNED_NUMBERS = [
    {
        "value": "0.3633",
        "what": "d2_h::wsi_biology raw joint LDA, P1 4.2",
        "artifact": "d2_final::d2_h_seed42",
        "papers": ["paper/P1_CALIBRA_DRAFT.md", "paper/P1_FIGURES.md"],
    },
    {
        "value": "0.2348",
        "what": "d2_i::wsi_biology raw joint LDA, P1 4.2",
        "artifact": "d2_final::d2_i_seed42",
        "papers": ["paper/P1_CALIBRA_DRAFT.md", "paper/P1_FIGURES.md"],
    },
    {
        "value": "0.463",
        "what": "cancer-type raw balanced accuracy, P1 4.2",
        "artifact": None,
        "unidentified_because":
            "PHASE1_RESULT.md states it in prose with no artifact path and no hash; no run output "
            "under p1_evidence/, p1_out/ or e0_run/ on the box records a cancer-type balanced "
            "accuracy of 0.463 or 0.035, and the cohort (n = 2,530) differs from the site arm "
            "(n = 2,766) so it cannot be inherited from it",
        "papers": ["paper/P1_CALIBRA_DRAFT.md", "paper/P1_FIGURES.md"],
    },
]

#: Phrases that count as declaring a number unattributable. One must appear near the number.
UNIDENTIFIED_MARKERS = ("artifact not identified", "provenance unresolved", "no identified artifact")

#: Basenames that name an artifact whose identity matters. Citing one without a hash nearby is the
#: defect this file exists to prevent.
HASH_REQUIRED_BASENAMES = ("d2_h_seed42.npz", "d2_i_seed42.npz")

#: How far from a citation a hash may sit and still count as attached to it.
WINDOW = 1600

_SHA256 = re.compile(r"\b[0-9a-f]{64}\b")
#: An abbreviated digest, as drafts write it inline: `4a18b94f1017b85d…`.
_SHA256_PREFIX = re.compile(r"\b[0-9a-f]{16}\b")


def _text(paper: str) -> str:
    return (REPO_ROOT / paper).read_text(encoding="utf-8")


def _papers_under_test() -> set[str]:
    return {paper for pinned in PINNED_NUMBERS for paper in pinned["papers"]}


def _occurrences(text: str, needle: str) -> list[int]:
    """Occurrences of ``needle`` that are not the prefix of a longer number.

    Without the trailing-digit guard, `0.463` matches `0.4637` in an unrelated table and the check
    reports a defect that is not there -- the same false-positive class `test_paper_paths_resolve`
    records having hit once already.
    """
    out, start = [], 0
    while True:
        index = text.find(needle, start)
        if index < 0:
            return out
        after = text[index + len(needle): index + len(needle) + 1]
        if not after.isdigit():
            out.append(index)
        start = index + 1


def _window(text: str, index: int) -> str:
    return text[max(0, index - WINDOW): index + WINDOW]


# --- the registry has to be well formed before it can guard anything --------------------

def test_every_recorded_digest_is_a_well_formed_sha256() -> None:
    for key, artifact in ARTIFACTS.items():
        assert _SHA256.fullmatch(artifact["sha256"]), f"{key}: not a 64-hex SHA-256"
        assert artifact["note"], f"{key}: needs a note saying what it is"


def test_every_pinned_number_names_a_known_artifact_or_none() -> None:
    for pinned in PINNED_NUMBERS:
        key = pinned["artifact"]
        assert key is None or key in ARTIFACTS, f"{pinned['value']}: unknown artifact {key!r}"
        if key is None:
            reason = pinned.get("unidentified_because", "")
            assert len(reason) > 40, (
                f"{pinned['value']}: an unidentified number must record WHY it could not be "
                "attributed, so the gap is a finding rather than an omission")
        for paper in pinned["papers"]:
            assert (REPO_ROOT / paper).exists(), f"{pinned['value']}: {paper} does not exist"


def test_the_registry_is_not_vacuous() -> None:
    """Guard against the checks passing because nothing is pinned."""
    assert len(PINNED_NUMBERS) >= 3
    assert any(p["artifact"] is None for p in PINNED_NUMBERS), (
        "at least one number is known to be unattributable; if that is ever fixed, remove this "
        "assertion in the same commit that fixes it")
    assert len({a["sha256"] for a in ARTIFACTS.values()}) == len(ARTIFACTS)


# --- the actual guarantee ---------------------------------------------------------------

@pytest.mark.parametrize("paper", sorted(_papers_under_test()))
def test_every_artifact_basename_in_a_draft_is_accompanied_by_a_hash(paper: str) -> None:
    """A bare filename is not an identifier. Three files share `d2_h_seed42.npz`.

    Every mention of one of these basenames must have a digest -- full or 16-hex abbreviated --
    within ``WINDOW`` characters, so a reader can tell *which* file is meant.
    """
    text = _text(paper)
    offenders: list[str] = []
    for basename in HASH_REQUIRED_BASENAMES:
        for index in _occurrences(text, basename):
            window = _window(text, index)
            if _SHA256.search(window) or _SHA256_PREFIX.search(window):
                continue
            line = text.count("\n", 0, index) + 1
            offenders.append(f"{paper}:{line} cites {basename} with no SHA-256 within "
                             f"{WINDOW} characters")
    assert not offenders, (
        "\n  ".join(["artifact filenames cited without a content hash:"] + offenders)
        + "\n\nThree distinct files are called d2_h_seed42.npz and they give raw joint LDA 0.3633, "
          "0.1782 and 0.3785. Name the SHA-256 beside the citation.")


@pytest.mark.parametrize("paper", sorted(_papers_under_test()))
def test_every_pinned_number_is_hash_identified_in_its_paper(paper: str) -> None:
    """Wherever a draft quotes a pinned number, the producing artifact's digest must be present.

    The full digest is required *somewhere* in the document -- so the paper commits to the bytes --
    and a digest (full or abbreviated) must sit within ``WINDOW`` of every occurrence of the number,
    so a reader at that point in the text can resolve it without hunting.
    """
    text = _text(paper)
    offenders: list[str] = []
    for pinned in PINNED_NUMBERS:
        if paper not in pinned["papers"] or pinned["artifact"] is None:
            continue
        digest = ARTIFACTS[pinned["artifact"]]["sha256"]
        occurrences = _occurrences(text, pinned["value"])
        if not occurrences:
            continue
        if digest not in text:
            offenders.append(f"{paper} quotes {pinned['value']} ({pinned['what']}) but never states "
                             f"the full SHA-256 {digest} of the artifact it came from")
            continue
        for index in occurrences:
            window = _window(text, index)
            if _SHA256.search(window) or _SHA256_PREFIX.search(window):
                continue
            line = text.count("\n", 0, index) + 1
            offenders.append(f"{paper}:{line} quotes {pinned['value']} ({pinned['what']}) with no "
                             f"digest within {WINDOW} characters")
    assert not offenders, "\n  ".join(["numbers quoted without hash identification:"] + offenders)


@pytest.mark.parametrize("paper", sorted(_papers_under_test()))
def test_an_unidentified_number_is_declared_unidentified(paper: str) -> None:
    """A number we cannot attribute must say so, near where it is quoted.

    This is the check that stops the registry being completed by guesswork: the tempting failure is
    not an absent hash, it is attributing 0.463 to whichever artifact happens to reproduce it.
    """
    text = _text(paper).lower()
    offenders: list[str] = []
    for pinned in PINNED_NUMBERS:
        if paper not in pinned["papers"] or pinned["artifact"] is not None:
            continue
        for index in _occurrences(text, pinned["value"]):
            if any(marker in _window(text, index) for marker in UNIDENTIFIED_MARKERS):
                continue
            line = text.count("\n", 0, index) + 1
            offenders.append(f"{paper}:{line} quotes {pinned['value']} ({pinned['what']}), which has "
                             "no identified artifact, without saying so")
    assert not offenders, (
        "\n  ".join(["unattributable numbers quoted as though attributed:"] + offenders)
        + "\n\nSay plainly that the artifact is not identified. Do not attribute it to whichever "
          "artifact reproduces it.")


# --- and the digests themselves, when the bytes are reachable ---------------------------

def test_recorded_digests_match_the_box_when_reachable() -> None:
    """On the box, re-hash. In a checkout, skip -- but say how many were skipped and why."""
    checked, mismatched = 0, []
    for key, artifact in ARTIFACTS.items():
        path = BOX_ROOT / artifact["path"]
        if not path.exists():
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        checked += 1
        if digest.hexdigest() != artifact["sha256"]:
            mismatched.append(f"{key}: recorded {artifact['sha256']}, actual {digest.hexdigest()}")
    if checked == 0:
        pytest.skip(f"none of the {len(ARTIFACTS)} artifacts are mounted here; "
                    f"expected {BOX_ROOT} (present only on the GPU box)")
    assert not mismatched, "\n  ".join(["recorded digests are stale:"] + mismatched)
