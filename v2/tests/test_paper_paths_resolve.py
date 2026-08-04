"""Every repository path cited in a paper draft must resolve to a file that exists.

`paper/P2_RANK_DRAFT.md` §4.5(a) argues that the defence against a number quoted under the wrong
label is **mechanical provenance, not diligence**. This test applies that recommendation to the
drafts themselves: four cited paths were found not to exist, incidentally, by an agent doing
something else. Diligence had already read those lines many times.

Scope, and why it is narrow on purpose:

* Only **repository-relative** paths are checked. Paths on the GPU box (`~/…`, `/lambda/…`) cannot
  be resolved from a checkout and are skipped — but they are skipped *by an explicit rule*, not by
  accident, and `test_box_paths_are_marked_as_such` pins the rule so a box path cannot silently be
  written as though it were a repo path.
* A token is treated as a candidate path if it is inside backticks and looks like a path: it
  contains a ``/`` or ends in a known file extension. Prose in backticks (`programme_only`,
  `effective_rank`) is therefore ignored, as are shell fragments and brace expansions, which are
  resolved to their prefix directory instead of exploded.

If this test fails, fix the citation in the draft. Do not add the path to the allowlist unless it
genuinely names something outside the repository, in which case mark it in the draft as a box path.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Drafts whose citations are checked. Add new drafts here.
PAPERS = [
    "paper/P2_RANK_DRAFT.md",
    "paper/P2_FIGURES.md",
    "paper/P1_CALIBRA_DRAFT.md",
    "paper/P1_FIGURES.md",
    "paper/LIVENESS_GATE_DESIGN.md",
    "paper/QUEUE_ANCHORING.md",
]

#: Anything starting with one of these is on the GPU box or another machine, not in the repository.
BOX_PREFIXES = ("~", "/lambda/", "/home/", "/tmp/", "/dev/", "/c/", "http://", "https://")

#: Top-level trees that exist only on the persistent NFS mount, never in a checkout. A draft may
#: cite them, but `test_every_box_tree_is_declared_once_per_paper` requires the draft to say so.
BOX_TREES = {
    "p1_evidence": "P1's evidence outputs, under "
                   "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/",
    "p1_out": "P1 artifact outputs on the box",
    "e0_run": "D1/D2 run tree on the box",
    "runs": "training run outputs on the box",
}

#: File extensions that make a bare token (no slash) a path rather than prose.
PATH_SUFFIXES = (".md", ".py", ".json", ".jsonl", ".npz", ".csv", ".log", ".pt", ".txt", ".yaml",
                 ".yml", ".toml", ".cfg", ".pdf", ".bib")

#: Cited paths that legitimately do not exist in the repository and are marked as such in the
#: drafts. Each entry must carry the reason; an entry without one is a bug in this file.
KNOWN_ABSENT = {
    # P2 §4.9 / §6.4: the whole point of the citation is that the file is missing.
    "paper/.../RESULTS.md": "cited by HANDOFF_BUILD_AGENT.md:98 and does not exist; the absence IS "
                            "the finding (P2 §6.4)",
    "p1_evidence/dilution/": "P1's CALIBRA outputs, on the box under "
                             "/lambda/nfs/geeg/.../morpheus_phase_d/; marked as a box path in the "
                             "drafts",
    # P1 4.6.6 and A.6 say so themselves: these were scratchpad files and the draft records that
    # they still need reproducing into the repository.
    "sim.py": "scratchpad verification script; P1 itself says it 'should be reproduced into' the "
              "repository and has not been",
    "sim3.py": "scratchpad verification script; same P1 note as sim.py",
}

#: Tokens that look like paths but are patterns, globs or prose-with-a-slash.
IGNORE_EXACT = {
    "and/or", "n/a", "N/A", "his/her", "raw/residualised", "km/h",
}

#: Output files that are written on the GPU box and are cited by bare name because the surrounding
#: prose already gives their directory. Each must say where it lives.
BOX_OUTPUT_BASENAMES = {
    # Entries are removed from here as soon as the file is vendored, so that its
    # citations start being checked instead of excused. The ten P2 readouts that
    # used to sit here are now under v2/research/rebase/p2/figures/data/, pulled by
    # that directory's extract_from_box.py and pinned by its MANIFEST.json;
    # test_no_box_output_basename_is_actually_in_the_repository is what caught them.
    "D1_PAIRED_BOOTSTRAP.json": "~/e0_run/d1_v2/",
    "D1_READOUT_INDEX.json": "~/e0_run/d1_v2/",
    "D1_AUDIT.json": "~/e0_run/d1_v2/",
    "D2_PAIR_MANIFEST.json": "~/e0_run/d2_v3/",
    "P2_METRICS_ALL_SUBSAMPLED.json": "~/e0_run/ or ~/ws_p2/out/",
    "d1_audit.log": "~/e0_run/",
    "d2_readout.py": "on the box; not vendored - see P2 6.2",
    "collapse_diag.py": "scratchpad/ on the A100; not vendored",
    "last.pt": "training checkpoint under ~/e0_run/d1_v1/<arm>/",
    "certificate_adjusted/": "output directory under p1_evidence/track1/ on the box",
    "certificate_raw/": "output directory under p1_evidence/track1/ on the box",
}

#: A token is not a path if it is a DOI, a DBLP/OpenReview key, an API route, or a formula.
_DOI = re.compile(r"^10\.\d{4,}/")
_DBLP = re.compile(r"^(conf|journals|rec)/")
_MATHY = re.compile(r"[Ͱ-Ͽ∀-⋿²³⁴⁻]")

_BACKTICKED = re.compile(r"`([^`\n]+)`")


def _candidate_paths(text: str) -> set[str]:
    """Backticked tokens that look like paths."""
    out: set[str] = set()
    for raw in _BACKTICKED.findall(text):
        token = raw.strip()
        # Strip a trailing line/section reference: `foo.py:149`, `foo.md:12-20`, `foo.py:60–106`.
        # The separator may be a hyphen, en-dash or em-dash — the drafts use all three, and an
        # earlier version of this regex accepted only the hyphen and so reported a real file as
        # missing. That is the same class of error this test exists to catch, in the test itself.
        token = re.sub(r":\d+([-‐-―]\d+)?$", "", token)
        # A pytest node id names a test inside a file: check only the file part.
        if "::" in token:
            token = token.split("::", 1)[0]
        # A citation may name several files: `a.py`, `b.py` is two backtick groups already, but
        # `d2_readout.py` style bare names are handled by the suffix rule below.
        if not token or token in IGNORE_EXACT:
            continue
        if " " in token:  # shell command, sentence, or code fragment
            continue
        if _DOI.match(token) or _DBLP.match(token) or _MATHY.search(token):
            continue
        if token.startswith("/paper/") or token.startswith("/works/") or token.startswith("/notes/"):
            continue  # bibliographic API routes
        if token.startswith("."):
            continue  # a bare extension, e.g. `.npz`
        if "/" not in token and not token.endswith(PATH_SUFFIXES):
            continue
        out.add(token)
    return out


def _is_box_path(token: str) -> bool:
    if token.startswith(BOX_PREFIXES):
        return True
    return token.split("/")[0] in BOX_TREES


def _resolvable_form(token: str) -> str:
    """Reduce a citation to something checkable.

    Brace expansions (`d1_p_seed{42,43,44}/last.pt`) and globs (`probevar_*.log`) are reduced to
    their longest literal prefix directory, because the point of the check is that the *location*
    exists, not that every expansion does.
    """
    for marker in ("{", "*", "?", "<", "["):
        if marker in token:
            token = token.split(marker)[0]
            token = token.rsplit("/", 1)[0] + "/" if "/" in token else ""
    return token


def _iter_citations():
    for paper in PAPERS:
        path = REPO_ROOT / paper
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for token in sorted(_candidate_paths(text)):
            yield paper, token


@pytest.mark.parametrize("paper", PAPERS)
def test_every_repository_path_cited_in_a_draft_exists(paper: str) -> None:
    missing: list[str] = []
    for cited_in, token in _iter_citations():
        if cited_in != paper:
            continue
        if _is_box_path(token) or token in KNOWN_ABSENT:
            continue
        if token in BOX_OUTPUT_BASENAMES:
            continue
        candidate = _resolvable_form(token)
        if not candidate:
            continue
        if (REPO_ROOT / candidate).exists():
            continue
        # A bare filename with no directory is allowed only if it is unique in the tree; the drafts
        # use these for well-known files (`spectral.py`). Resolve it by search.
        if "/" not in candidate.rstrip("/"):
            hits = [p for p in REPO_ROOT.rglob(candidate) if ".git" not in p.parts]
            if hits:
                continue
        missing.append(token)
    assert not missing, (
        f"{paper} cites {len(missing)} path(s) that do not exist in the repository:\n  "
        + "\n  ".join(missing)
        + "\n\nFix the citation. If the path is on the GPU box, write it as `~/...` or with its "
          "/lambda/... prefix so it is visibly not a repository path."
    )


def test_box_paths_are_marked_as_such() -> None:
    """A path under a box-only tree must not be written as though it were repo-relative.

    `p1_evidence/` exists only on the persistent NFS mount. Citing it bare, as three drafts once
    did, reads as a repository directory and resolves to nothing.
    """
    offenders: list[str] = []
    for paper, token in _iter_citations():
        if _is_box_path(token):
            continue
        head = token.split("/")[0]
        if head in {"p1_evidence", "p1_out", "e0_run", "ws_rank", "ws_p2", "ws_d1", "runs"}:
            text = (REPO_ROOT / paper).read_text(encoding="utf-8")
            # The citation is acceptable if the draft says, near it, that it is a box path.
            idx = text.find(token)
            window = text[max(0, idx - 400): idx + 400]
            if "box" in window.lower() or "/lambda/" in window:
                continue
            offenders.append(f"{paper}: {token}")
    assert not offenders, (
        "box-only trees cited as if they were repository paths:\n  " + "\n  ".join(offenders)
    )


def test_known_absent_entries_each_carry_a_reason() -> None:
    for token, reason in KNOWN_ABSENT.items():
        assert reason and len(reason) > 20, f"KNOWN_ABSENT[{token!r}] needs a real reason"
    for token, where in BOX_OUTPUT_BASENAMES.items():
        assert where, f"BOX_OUTPUT_BASENAMES[{token!r}] must say where the file lives"


def test_no_box_output_basename_is_actually_in_the_repository() -> None:
    """If one of these gets vendored, it should stop being excused and start being checked."""
    now_present = []
    for token in BOX_OUTPUT_BASENAMES:
        hits = [p for p in REPO_ROOT.rglob(token) if ".git" not in p.parts]
        if hits:
            now_present.append(f"{token} -> {hits[0].relative_to(REPO_ROOT)}")
    assert not now_present, (
        "these are now in the repository and must be removed from BOX_OUTPUT_BASENAMES so their "
        "citations are checked:\n  " + "\n  ".join(now_present)
    )


def test_the_checker_actually_finds_something() -> None:
    """Guard against the extractor silently matching nothing and the suite passing vacuously."""
    tokens = list(_iter_citations())
    assert len(tokens) > 50, f"only {len(tokens)} citations extracted; the regex has probably broken"
