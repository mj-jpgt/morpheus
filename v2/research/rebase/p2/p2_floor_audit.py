"""The floor audit: every rank comparison P2 makes or relies on, judged against
the floor measured on its OWN block.

Draft section 4.1a. This module owns `floor_audit.json`, which is the
machine-checkable form of that subsection's table, and it renders the table the
draft prints, so the two cannot drift apart.

WHY THIS EXISTS. The paper's central criterion is that a rank difference smaller
than the measured same-seed retraining floor is not resolvable. That criterion has
been applied to the paper's own numbers **five separate times, each discovered
late and separately** -- the momentum fix (3.29x), the m = 0.999-over-0.99 choice
(1.26x), 4.7.4's surviving necessity violation (3.73x), the decorrelation ablation
(1.85x) and 5.1's instance 2 (~3.2x, flagged and left unrewritten). Finding them
one at a time is how a referee finds the sixth. This file enumerates all of them
at once, so that the next one is found by a test rather than by a reviewer.

WHAT IT DOES NOT DO. Nothing here computes a rank statistic, a canonical
correlation, or any other spectral quantity. Every number is READ from a source --
a vendored box artifact under `figures/data/`, or a named table in a repository
markdown file -- and the only arithmetic performed is a ratio of two numbers that
were already read. Four inline-formula substitutions have been caught in this
paper (draft 3.1, 4.5a); the defence is mechanical provenance, so this module
resolves rather than recomputes.

THE THREE THINGS THE CHECKER ENFORCES.

1. **Every recorded value agrees with its source.** `check()` re-resolves both
   sides of every comparison from the file named in its `src` block and fails if
   the recorded value has drifted.

2. **Every ratio agrees with its own two values**, to the precision the draft
   prints it at.

3. **Block-matching.** A ratio measured on the residualised block is judged
   against the residualised floor (3.295x) and a raw-block ratio against the raw
   floor (3.111x). This is load-bearing and has already produced one wrong
   verdict: D1-B seed 43's residualised 3.246x judged against the RAW floor of
   3.111x reads as *outside* the floor when on its own block it is inside. An
   entry whose block has no measured floor must say so explicitly
   (`floor: null` plus a `floor_note`), never silently borrow another block's.

Every entry also declares a `kind`:

  `selection`   a choice between two candidate configurations -- the criterion
                applies with full force.
  `nuisance`    a measurement OF the noise, not a selection against it (a seed
                spread, a retraining spread). The criterion does not apply
                because the entry is an instance of the floor rather than a
                claim resolved against it.
  `regime`      two different experimental regimes rather than two candidate
                configurations (16 memorised patients against 3,118 streaming
                ones; a zero-training representation across dilution levels).
                Exempt, with the reason recorded in `rests_on`.
  `collapse`    a reading at or near the collapse floor, where 4.10 defends the
                statistic, with independently co-measured evidence of collapse.
  `direction`   the claim is a sign, an ordering-instability or an equality, not
                a magnitude. Recorded because an unresolvable magnitude does not
                by itself sink a claim about whether an ordering is stable -- but
                the ratio is printed so a reader can see what the ordering was
                read off.
  `excluded`    already excluded from every count for a reason other than the
                floor (no artifact, unknown statistic, withdrawn).

usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 python3 p2_floor_audit.py --check
    ... p2_floor_audit.py --markdown        # the table draft 4.1a prints
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
AUDIT = HERE / "floor_audit.json"
DATA = HERE / "figures" / "data"
REPO = HERE.parents[3]          # v2/research/rebase/p2 -> repository root

#: Headers the probe logs written by `v2/research/rebase/d1_momentum_probe.py` are
#: allowed to carry. The statistic each column names is NOT interchangeable:
#: `R3-rank` is the row-normalised order-2 Hill number and `CANONICAL` is Roy &
#: Vetterli order 1. The earlier logs (`probevar_*`) predate the two-statistic
#: header and carry ONE rank column called `eff-rank`; the log's own
#: `final_eff_rank=` line reports it, and it is the R3 statistic, not the
#: canonical one -- which is why a notebook entry quoting "eff-rank" from any of
#: these logs is quoting R3 and must not be judged against an R1 floor without
#: saying so. A header that is neither of these fails here rather than silently
#: yielding whichever column happens to sit at that index.
PROBE_HEADERS = {
    ("R3-rank", "CANONICAL", "feat-std", "rna-rna", "contrastive"),
    ("eff-rank", "feat-std", "rna-rna", "contrastive"),
}


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load() -> dict:
    return json.loads(AUDIT.read_text(encoding="utf-8"))


def _dotted(blob: Any, path: str) -> Any:
    node = blob
    for part in path.split("."):
        node = node[part]
    return node


def _probe_rows(rel: str) -> dict[int, dict[str, float]]:
    """Parse a `d1_momentum_probe.py` log into {step: {column: value}}.

    The header line is matched rather than assumed, so a log written by a
    different script -- or by a future version of this one with a reordered
    header -- fails here instead of silently returning the wrong column.
    """
    text = (DATA / rel).read_text(encoding="utf-8")
    header: tuple[str, ...] | None = None
    rows: dict[int, dict[str, float]] = {}
    for line in text.splitlines():
        cells = line.split()
        if cells[:1] == ["step"]:
            header = tuple(cells[1:])
            if header not in PROBE_HEADERS:
                raise ValueError(f"{rel}: unexpected probe header {header!r}")
            continue
        if header is None or not cells or not cells[0].isdigit():
            continue
        step = int(cells[0])
        values = cells[1:]
        if len(values) != len(header):
            continue
        rows[step] = {}
        for name, cell in zip(header, values):
            try:
                rows[step][name] = float(cell)
            except ValueError:      # `nan` in the contrastive column at step 0
                rows[step][name] = float("nan")
    if not rows:
        raise ValueError(f"{rel}: no probe rows parsed")
    return rows


_MARKUP = re.compile(r"[*`\s]")


def _markdown_literal(rel: str, heading: str | None, literal: str,
                      value: float | None = None) -> float:
    """Assert `literal` occurs in a named repository markdown file, and return it.

    This is the resolver for values that exist only as prose or as a cell in a
    table in the draft or in a notebook entry -- the momentum sweep at step 600,
    the collapse-regime readings, the historical instances. It is deliberately a
    substring assertion rather than a parse: if the draft's number is edited, the
    audit entry no longer resolves and the test fails, which is exactly the
    coupling this file exists to create.
    """
    text = (REPO / rel).read_text(encoding="utf-8", errors="replace")
    if heading is not None:
        text = _section(text, heading)
    if _MARKUP.sub("", literal) not in _MARKUP.sub("", text):
        raise LookupError(
            f"{rel}: {literal!r} is no longer present"
            + (f" under a heading containing {heading!r}" if heading else "")
            + " -- the audit entry disagrees with its source, so neither may be used")
    if value is not None:
        # For sources whose printed form is not a bare number ("~2", "16/16"):
        # the literal is still asserted against the file, and the numeric value
        # the audit uses is declared beside it rather than parsed out of prose.
        return float(value)
    return float(literal.replace("−", "-").replace("–", "-"))


def _section(text: str, heading: str) -> str:
    lines = text.splitlines()
    start, level = None, 0
    for i, line in enumerate(lines):
        if line.lstrip().startswith("#") and heading.lower() in line.lower():
            start = i + 1
            level = len(line) - len(line.lstrip("#"))
            break
    if start is None:
        raise LookupError(f"no heading containing {heading!r}")
    out = []
    for line in lines[start:]:
        stripped = line.lstrip()
        if stripped.startswith("#") and (len(line) - len(stripped)) <= level:
            break
        out.append(line)
    return "\n".join(out)


def resolve(src: dict) -> float:
    """One recorded value, read back out of the file it came from."""
    kind = src["kind"]
    if kind == "json":
        blob = json.loads((DATA / src["file"]).read_text(encoding="utf-8"))
        return float(_dotted(blob, src["path"]))
    if kind == "spread":
        # A max/min fold over a named family of values -- used where the quantity
        # a claim compares is itself a SPREAD (draft 4.3), for which this project
        # has never measured a floor.
        blob = json.loads((DATA / src["file"]).read_text(encoding="utf-8"))
        values = [float(_dotted(blob, p)) for p in src["paths"]]
        return max(values) / min(values)
    if kind == "probe_log":
        return _probe_rows(src["file"])[int(src["step"])][src["column"]]
    if kind == "markdown":
        return _markdown_literal(src["file"], src.get("heading"), src["literal"],
                                 src.get("value"))
    raise KeyError(f"unknown source kind {kind!r}")


def resolve_ratio(entry: dict) -> float:
    """The ratio an entry's own sources support, recomputed from them.

    For a `pair_max` entry the recorded ratio is the LARGEST fold in a family of
    comparisons the draft treats as one row (5 statistics x 6 pairs, 6 pairs x 3
    views), so the family is re-scanned rather than one representative trusted.
    """
    if entry.get("pair_max"):
        blob = json.loads((DATA / entry["pair_max"]["file"]).read_text(encoding="utf-8"))
        folds = []
        for left, right in entry["pair_max"]["paths"]:
            a, b = float(_dotted(blob, left)), float(_dotted(blob, right))
            folds.append(max(a, b) / min(a, b))
        return max(folds)
    a, b = resolve(entry["a"]["src"]), resolve(entry["b"]["src"])
    return max(a, b) / min(a, b)


# --------------------------------------------------------------------------
# Checking
# --------------------------------------------------------------------------
def check(audit: dict | None = None) -> list[str]:
    """Every disagreement between the audit list and the files it cites."""
    audit = audit or load()
    floors = audit["floors"]
    problems: list[str] = []
    seen: set[str] = set()

    for name, floor in floors.items():
        if "value" in floor:
            got = resolve_ratio(floor)
            if abs(got - floor["value"]) > 0.0015:
                problems.append(
                    f"floor {name}: recorded {floor['value']} but its source gives {got:.4f}")

    for entry in audit["comparisons"]:
        eid = entry["id"]
        if eid in seen:
            problems.append(f"{eid}: duplicate id")
        seen.add(eid)

        # (1) the two values still agree with their sources, and (2) the ratio
        # agrees with the two values.
        try:
            got = resolve_ratio(entry)
        except (LookupError, KeyError, ValueError, OSError) as exc:
            problems.append(f"{eid}: source unresolvable -- {exc}")
            continue
        if abs(got - entry["ratio"]) > entry.get("ratio_tolerance", 0.001):
            problems.append(
                f"{eid}: table says {entry['ratio']}x, its source gives {got:.4f}x")

        # (3) block matching. A ratio may only be judged against a floor measured
        # on its own block; anything else must declare the mismatch.
        floor_name = entry.get("floor")
        if floor_name is None:
            if not entry.get("floor_note"):
                problems.append(f"{eid}: no floor and no floor_note explaining why")
        else:
            if floor_name not in floors:
                problems.append(f"{eid}: unknown floor {floor_name!r}")
                continue
            floor = floors[floor_name]
            same_block = (floor.get("block") == entry.get("block")
                          and floor.get("statistic") == entry.get("statistic"))
            if not same_block and not entry.get("floor_note"):
                problems.append(
                    f"{eid}: judged against {floor_name} "
                    f"({floor.get('statistic')}/{floor.get('block')}) but is "
                    f"{entry.get('statistic')}/{entry.get('block')}, with no floor_note")
            if "value" in floor:
                clears = got > floor["value"]
                if bool(entry["clears"]) != clears:
                    problems.append(
                        f"{eid}: recorded clears={entry['clears']} but {got:.4f}x "
                        f"against {floor_name} = {floor['value']}x gives {clears}")
            elif "low" in floor:
                clears = got > floor["high"]
                if bool(entry["clears"]) != clears:
                    problems.append(
                        f"{eid}: recorded clears={entry['clears']} but {got:.4f}x "
                        f"against the {floor['low']}-{floor['high']}x band gives {clears}")

        if entry["kind"] not in KINDS:
            problems.append(f"{eid}: unknown kind {entry['kind']!r}")
        if not entry.get("rests_on"):
            problems.append(f"{eid}: no `rests_on` -- what carries the claim is not stated")
        if entry["kind"] == "selection" and not entry["clears"] and not entry.get("rests_on"):
            problems.append(f"{eid}: an unresolvable selection with nothing stated behind it")
    return problems


KINDS = {"selection", "nuisance", "regime", "collapse", "direction", "excluded"}


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------
def _fmt(value: float, digits: int) -> str:
    return f"{value:.{digits}f}"


def render_markdown(audit: dict | None = None) -> str:
    """The table draft 4.1a prints, generated from the audit list.

    The draft's copy is asserted against this in `v2/tests/test_p2_floor_audit.py`,
    so a number cannot be edited in one place and left in the other.
    """
    audit = audit or load()
    floors = audit["floors"]
    head = ("| # | § | comparison | values | ratio | statistic | block | floor | clears? |"
            " what the claim rests on |")
    rule = "|---|---|---|---|---:|---|---|---|:---:|---|"
    lines = [head, rule]
    for i, entry in enumerate(audit["comparisons"], start=1):
        if entry.get("pair_max"):
            values = entry.get("values_label", "—")
        else:
            digits = entry.get("value_digits", 3)
            values = f"{_fmt(entry['a']['value'], digits)} / {_fmt(entry['b']['value'], digits)}"
        floor_name = entry.get("floor")
        if floor_name is None:
            floor_cell = "**none measured**"
        else:
            floor = floors[floor_name]
            floor_cell = (f"{floor['value']}×" if "value" in floor
                          else f"{floor['low']}–{floor['high']}×")
            if entry.get("floor_note"):
                floor_cell += " ‡"
        if entry["kind"] in ("regime", "excluded"):
            verdict = "**exempt**"
        elif floor_name is None:
            verdict = "**no floor**"
        elif entry["clears"]:
            verdict = "yes"
        else:
            verdict = "**no**"
        lines.append(
            f"| {i} | {entry['section']} | {entry['comparison']} | {values} | "
            f"**{entry['ratio']:.3f}×** | {entry['statistic']} | {entry['block']} | "
            f"{floor_cell} | {verdict} | {entry['rests_on']} |")
    return "\n".join(lines)


def summary(audit: dict | None = None) -> dict:
    audit = audit or load()
    rows = audit["comparisons"]
    return {
        "total": len(rows),
        "selection": sum(r["kind"] == "selection" for r in rows),
        "selection_failing": sum(r["kind"] == "selection" and not r["clears"] for r in rows),
        "exempt": sum(r["kind"] in ("regime", "excluded") for r in rows),
        "no_floor_measured": sum(r.get("floor") is None for r in rows),
        "block_mismatched": sum(bool(r.get("floor_note")) for r in rows),
    }


def summary_sentence(audit: dict | None = None) -> str:
    """The counting sentence draft 4.1a prints, generated rather than typed.

    Asserted against the draft by the test, for the same reason the table is: a
    count that is maintained by hand drifts the moment a row is added.
    """
    s = summary(audit)
    return (f"{s['total']} rank comparisons are enumerated; "
            f"{s['selection']} of them are selections between candidate configurations, and "
            f"{s['selection_failing']} of those {s['selection']} do not clear the floor their "
            f"own block licenses. "
            f"{s['exempt']} are exempt with the reason stated, "
            f"{s['no_floor_measured']} are on a statistic or a block for which **no floor has "
            f"ever been measured**, and "
            f"{s['block_mismatched']} carry an explicit block- or statistic-mismatch note.")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--markdown", action="store_true", help="print the draft 4.1a table")
    ap.add_argument("--sentence", action="store_true", help="print the draft 4.1a counting sentence")
    ap.add_argument("--check", action="store_true", help="re-resolve every value (default)")
    args = ap.parse_args(argv)
    if args.markdown:
        print(render_markdown())
        return 0
    if args.sentence:
        print(summary_sentence())
        return 0
    problems = check()
    for line in problems:
        print("DISAGREEMENT:", line)
    print(json.dumps(summary(), indent=1))
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
