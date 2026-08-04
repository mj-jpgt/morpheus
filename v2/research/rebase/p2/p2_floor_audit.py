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
#: The last of these is the CAPACITY-SWEEP header, from the older
#: `decorr_causal.py` harness behind §5.2 measurement 3. It carries the same four
#: columns plus a trailing `decorr` column reporting the decorrelation LOSS TERM,
#: which is not a rank and is never read here. It was previously refused, and
#: naming it rather than guessing is what recovered §5.2 measurement 3's reading
#: step: both of its arms print a row at **step 150**, which §6.2 recorded as
#: never having been written down.
PROBE_HEADERS = {
    ("R3-rank", "CANONICAL", "feat-std", "rna-rna", "contrastive"),
    ("eff-rank", "feat-std", "rna-rna", "contrastive"),
    ("eff-rank", "feat-std", "rna-rna", "contrastive", "decorr"),
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


def resolve_shape(src: dict) -> dict:
    """A floor's bimodality descriptor, read back out of the file it was measured into."""
    blob = json.loads((DATA / src["file"]).read_text(encoding="utf-8"))
    return _dotted(blob, src["path"])


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
        # A SECOND, independent source for the same floor where one exists. The
        # two floors the paper's whole criterion rests on are parsed from
        # `d1_envelope_readout.py`'s printed log; the cross-check recomputes the
        # same statistic on the same five artifacts from a workspace verified
        # file-by-file against `git ls-tree`. If the two ever part, neither may
        # be used -- which is the rule, not a preference for one of them.
        if "cross_check" in floor:
            other = resolve_ratio(floor["cross_check"])
            if abs(other - floor["value"]) > 0.0015:
                problems.append(
                    f"floor {name}: recorded {floor['value']} but its independent cross-check "
                    f"({floor['cross_check']['what']}) gives {other:.4f} -- STOP, and report both")
        # The shape is a claim about the paper's thesis, so it is checked like a
        # value: whether the divergent run is an outlier under THIS statistic, and
        # whether the four-run agreement survives, are read back out of the file
        # they were measured into.
        if "shape" in floor and "shape_src" in floor:
            source = resolve_shape(floor["shape_src"])
            for key in ("outlier", "outlier_is_low", "bimodal"):
                if floor["shape"][key] != source[key]:
                    problems.append(
                        f"floor {name}: shape.{key} recorded {floor['shape'][key]!r} but its "
                        f"source gives {source[key]!r}")
            if abs(floor["shape"]["rest_fold"] - source["rest_fold"]) > 0.0002:
                problems.append(
                    f"floor {name}: shape.rest_fold recorded {floor['shape']['rest_fold']} but "
                    f"its source gives {source['rest_fold']:.4f}")

    for entry in audit["comparisons"]:
        eid = entry["id"]
        if eid in seen:
            problems.append(f"{eid}: duplicate id")
        seen.add(eid)

        # (0) UNJUDGEABLE IS NOT A VERDICT. A comparison whose statistic-and-block
        # has no measured floor may record neither a pass nor a failure, and the
        # coupling is enforced rather than left to convention: for as long as
        # `clears: false` was allowed to sit beside `floor: null`, thirteen rows
        # with no ruler at all were being counted alongside rows that had failed
        # a measured one.
        if (entry.get("floor") is None) != (entry.get("clears") is None):
            problems.append(
                f"{eid}: floor={entry.get('floor')!r} but clears={entry.get('clears')!r} — a row "
                f"with no floor must record `clears: null`, and a row with a floor must record a "
                f"boolean")

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

    # The floors must each name a statistic and a block, or block-matching cannot
    # be checked at all, and each must carry the caveat that travels with it.
    for name, floor in floors.items():
        for field in ("statistic", "block", "caveat"):
            if not floor.get(field):
                problems.append(f"floor {name}: no {field}")
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
            # Not a pass and not a failure: no floor has ever been measured for
            # this row's statistic on this row's block, so the criterion cannot
            # be applied to it in either direction.
            verdict = "**unjudgeable**"
        elif entry["clears"]:
            verdict = "yes"
        else:
            verdict = "**no**"
        lines.append(
            f"| {i} | {entry['section']} | {entry['comparison']} | {values} | "
            f"**{entry['ratio']:.3f}×** | {entry['statistic']} | {entry['block']} | "
            f"{floor_cell} | {verdict} | {entry['rests_on']} |")
    return "\n".join(lines)


#: Display order for the floor table: the block 4.1's floor is measured on
#: first, then the paired block LiDAR needs, then the two views nothing had ever
#: been measured on. Within a group, largest floor first, so the 3.3x spread
#: BETWEEN statistics on one block is the first thing the column shows.
_FLOOR_GROUPS = [
    ("the exported `wsi_biology` block — the one §4.1 measures, every statistic T1 scores",
     lambda name: (not name.startswith("LiDAR") and "_probe_step" not in name
                   and not name.endswith(("_rna_view", "_full_view")))),
    ("the `wsi_biology` + `rna_biology` positive pair — LiDAR's block",
     lambda name: name.startswith("LiDAR")),
    ("the `rna_biology` view, canonical R1 — same five runs, same statistic, other view",
     lambda name: name.endswith("_rna_view")),
    ("the `full_biology` view, canonical R1 — same five runs, same statistic, other view",
     lambda name: name.endswith("_full_view")),
    # The block draft §5's every rank number is measured on. It is listed last
    # because it is the newest and its rows read differently from the four above:
    # they are per STEP, not per view, and the `divergent run` column is a repeat
    # of whichever ARM carried the larger fold rather than of a single arm.
    ("the fixed held-out probe — the block every §5 rank number sits on, per reading step",
     lambda name: "_probe_step" in name),
]


def render_floors_markdown(audit: dict | None = None) -> str:
    """The measured-floor table draft 4.1a prints, generated from the audit list.

    Min and max are re-resolved from `P2_ENVELOPE_FLOORS.json` rather than
    restated, so the table cannot drift from the file the floors were measured
    into. The last three columns are the SHAPE, not the magnitude: whether the
    four-run agreement and the factor-scale gap that make 4.1's floor bimodal
    survive a change of statistic is a claim about this paper's thesis, not a
    detail of presentation.
    """
    audit = audit or load()
    floors = audit["floors"]
    lines = ["| statistic | block | min | max | **floor** | other four agree to |"
             " divergent run | bimodal? |",
             "|---|---|---:|---:|---:|---:|---|:---:|"]
    placed: set[str] = set()
    for label, predicate in _FLOOR_GROUPS:
        group = [(n, f) for n, f in floors.items()
                 if n not in placed and predicate(n) and "value" in f]
        group.sort(key=lambda item: (-item[1]["value"], item[0]))
        if not group:
            continue
        lines.append(f"| **{label}** | | | | | | | |")
        for name, floor in group:
            placed.add(name)
            lo, hi = resolve(floor["b"]["src"]), resolve(floor["a"]["src"])
            shape = floor.get("shape", {})
            block = floor["block"].split(" (")[0].replace(", `", " — `")
            lines.append(
                f"| {floor['statistic'].replace('|', chr(92) + '|')} | {block} | "
                f"{lo:.4f} | {hi:.4f} | **{floor['value']:.3f}×** | "
                f"{shape['rest_fold']:.3f}× | {shape['outlier']} "
                f"({'low' if shape['outlier_is_low'] else 'high'}) | "
                f"{'**yes**' if shape['bimodal'] else 'no'} |")
    return "\n".join(lines)


def summary(audit: dict | None = None) -> dict:
    """The three-way count. `failing` + `clearing` + `unjudgeable` = `selection`.

    The three-way split is the point. A row whose block has no measured floor is
    NOT a failure, and counting it as one flatters the paper: it says the
    criterion was applied and the comparison lost, when in fact the criterion
    could not be applied at all.
    """
    audit = audit or load()
    rows = audit["comparisons"]
    sel = [r for r in rows if r["kind"] == "selection"]
    return {
        "total": len(rows),
        "selection": len(sel),
        "selection_failing": sum(r["clears"] is False for r in sel),
        "selection_clearing": sum(r["clears"] is True for r in sel),
        "selection_unjudgeable": sum(r.get("floor") is None for r in sel),
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
    # "clears"/"clear" and the reason clause both used to be hard-coded for a
    # count of one and for a block with no floor at all. Neither is true now: the
    # fixed held-out probe HAS a floor, and what stops the remaining rows is that
    # it was not measured at their reading step or on their arm.
    verb = "clears" if s["selection_clearing"] == 1 else "clear"
    return (f"{s['total']} rank comparisons are enumerated; "
            f"{s['selection']} of them are selections between candidate configurations. "
            f"**{s['selection_failing']} of those {s['selection']} do not clear the floor their own "
            f"statistic and block license, {s['selection_clearing']} {verb} it, and "
            f"{s['selection_unjudgeable']} cannot be judged at all** because no floor has been "
            f"measured for the statistic, block, reading step and arm they sit on. "
            f"{s['exempt']} rows are exempt with the reason stated, "
            f"{s['no_floor_measured']} of the {s['total']} are unjudgeable for want of a floor, and "
            f"{s['block_mismatched']} carry an explicit statistic-, block- or kind-mismatch note.")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--markdown", action="store_true", help="print the draft 4.1a table")
    ap.add_argument("--floors", action="store_true", help="print the draft 4.1a floor table")
    ap.add_argument("--sentence", action="store_true", help="print the draft 4.1a counting sentence")
    ap.add_argument("--check", action="store_true", help="re-resolve every value (default)")
    args = ap.parse_args(argv)
    if args.markdown:
        print(render_markdown())
        return 0
    if args.floors:
        print(render_floors_markdown())
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
