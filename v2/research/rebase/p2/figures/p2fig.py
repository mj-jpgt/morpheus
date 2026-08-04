"""Shared style, data access and provenance helpers for the P2 figures.

Matplotlib only. No seaborn. Nothing here imports `morpheus`, so every figure
script runs standalone from this directory with numpy + matplotlib alone.

Three things this module exists to enforce.

**No hand-entered numbers.** Every accessor here reads a file — either a
vendored box artifact under `data/` (see `extract_from_box.py`) or a markdown
source in the repository. Figures call these accessors; they do not carry
literals. Where a value survives only as prose in a notebook entry, the parser
for that entry lives here and names the file and the line it read.

**The statistic and the block are never implicit.** Four rank statistics have
been called `effective_rank` in this repository and a fifth lives under two
names (`P2_FIGURES.md` binding constraint 1). `STAT` holds the labels; axis
helpers require one. `axis_label()` refuses to build a rank axis label without
both a statistic and a block.

**A pending measurement is loud.** `pending_panel()` draws a hatched placeholder
naming the blocking measurement and the path its data will arrive at, and
records the block so `make_all.py` can report it. Under `--strict` (what the
test suite uses) it raises `PendingMeasurement` instead. It never silently
plots nothing.

Palette: Okabe & Ito, colourblind-safe, print-targeted. Validated with the
dataviz skill's `validate_palette.js` on 2026-08-04:

    Palette (light, surface #fcfcfb, categorical): 6 slots
      [PASS] Lightness band     all 6 inside L 0.43-0.77
      [PASS] Chroma floor       all 6 >= 0.1
      [WARN] CVD separation     worst all-pairs #CC79A7<->#009E73 dE 7.6 (deutan)
      [PASS] Normal-vision floor worst all-pairs #E69F00<->#D55E00 dE 15.6
      [WARN] Contrast vs surface below 3:1: #E69F00, #CC79A7, #56B4E9
      -> ALL CHECKS PASS

Both WARNs are relieved by design and not by assertion: every series in every
panel carries a distinct marker shape *and* a distinct linestyle or hatch, so
identity never rests on hue (this is also what makes the figures legible in
greyscale, which the print brief requires); and every series is directly
labelled with its value, which is the relief the contrast WARN obligates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "out"
REPO = HERE.parents[4]  # v2/research/rebase/p2/figures -> repository root


# --------------------------------------------------------------------------
# Palette. Fixed order, never cycled.
# --------------------------------------------------------------------------
class C:
    """Okabe & Ito, in fixed assignment order."""

    BLUE = "#0072B2"        # slot 1 - rank / the statistic under test
    VERMILLION = "#D55E00"  # slot 2 - the information channel
    GREEN = "#009E73"       # slot 3
    ORANGE = "#E69F00"      # slot 4
    PURPLE = "#CC79A7"      # slot 5
    SKY = "#56B4E9"         # slot 6

    INK = "#1a1a1a"
    MUTED = "#5a5a5a"
    FAINT = "#9a9a9a"
    RULE = "#d4d4d4"
    ENVELOPE = "#b9b9b9"    # the n=1 envelope: grey, never a series hue
    SURFACE = "#ffffff"


ORDER = [C.BLUE, C.VERMILLION, C.GREEN, C.ORANGE, C.PURPLE, C.SKY]

#: Markers and linestyles carry identity alongside hue, so the figures survive
#: greyscale printing and the CVD floor-band warning above.
MARKERS = ["o", "s", "^", "D", "v", "P"]
LINESTYLES = ["-", (0, (5, 2)), (0, (1, 1.4)), (0, (7, 2, 1, 2)), (0, (3, 1, 1, 1, 1, 1)), (0, (4, 1))]
HATCHES = ["", "///", "...", "\\\\\\", "xxx", "+++"]


# --------------------------------------------------------------------------
# The statistics, by name. Binding constraint 1.
# --------------------------------------------------------------------------
STAT = {
    "R1": "R1 (Roy & Vetterli order-1 effective rank, centred; v2/calibra/spectral.py CANONICAL)",
    "R1_short": "R1 (canonical effective rank, centred)",
    "R2": "R2 (order-2 Hill number of the centred singular values, $(\\Sigma\\sigma)^2/\\Sigma\\sigma^2$)",
    "R2_short": "R2 (order-2 Hill number, centred)",
    "R3": "R3 (order-2 Hill number after L2-normalising rows)",
    "R3_short": "R3 (order-2 Hill, row-normalised)",
    "PR": "PR (eigenvalue participation ratio $(\\Sigma\\sigma^2)^2/\\Sigma\\sigma^4$)",
    "PR_rownorm": "PR$_{\\mathrm{rownorm}}$ (eigenvalue participation ratio, row-normalised)",
    "RANKME": "RankMe as published (uncentred, $\\varepsilon=10^{-7}$)",
    "HARD": "hard numerical matrix rank (torch.linalg.matrix_rank) - NOT an effective rank",
    "STABLE": "stable rank $\\|A\\|_F^2/\\|A\\|_2^2$",
    "LIDAR": "LiDAR (adapted, q = 2, modalities as views)",
    "ALPHA": "$\\alpha$-ReQ $|\\alpha-1|$ (fastssl estimator)",
}

BLOCK = {
    "raw": "raw exported block",
    "resid": "confound-residualised block (cancer + TSS, 84 sites)",
    "train": "training batch, in-run",
    "probe": "fixed held-out probe",
}


def axis_label(stat: str, block: str, *, quantity: str = "effective rank") -> str:
    """Build a rank axis label. Both the statistic and the block are mandatory.

    This project has four statistics under three names and a raw-versus-
    residualised distinction that flips arm orderings (P2_FIGURES.md F5(b)), so
    an unqualified "effective rank" axis is not a well-defined object.
    """
    if stat not in STAT:
        raise KeyError(f"unknown statistic {stat!r}; one of {sorted(STAT)}")
    if block not in BLOCK:
        raise KeyError(f"unknown block {block!r}; one of {sorted(BLOCK)}")
    return f"{quantity}\n{STAT[stat]}\n{BLOCK[block]}"


# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------
def apply_style() -> None:
    plt.rcParams.update({
        "figure.facecolor": C.SURFACE,
        "axes.facecolor": C.SURFACE,
        "savefig.facecolor": C.SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 7.5,
        "axes.labelsize": 7.5,
        "axes.titlesize": 8.5,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelcolor": C.INK,
        "axes.edgecolor": C.RULE,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "xtick.color": C.MUTED,
        "ytick.color": C.MUTED,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "text.color": C.INK,
        "legend.fontsize": 6.8,
        "legend.frameon": False,
        "lines.linewidth": 1.6,
        "lines.markersize": 5.0,
        "pdf.fonttype": 42,     # embed as TrueType, not Type 3: editable vector text
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # keep text as text in the SVG
        "figure.dpi": 110,
    })


def grid(ax, axis: str = "both") -> None:
    ax.set_axisbelow(True)
    ax.grid(True, axis=axis, color=C.RULE, linewidth=0.5, alpha=0.9)


def caption(fig, text: str, *, y: float = 0.012, size: float = 6.0) -> None:
    """Print the binding caption inside the artwork.

    The plan requires several statements to travel with the figure rather than
    with a manuscript caption a reader may not have (F8's withdrawal, T1's
    underpowering, F6's negative). Rendering them into the figure makes that
    unconditional.
    """
    import textwrap
    width = int(fig.get_size_inches()[0] * 72 / (size * 0.505))
    body = "\n".join(textwrap.fill(" ".join(p.split()), width)
                     for p in text.strip().split("\n\n"))
    fig.text(0.012, y, body, ha="left", va="bottom", fontsize=size,
             color=C.MUTED, linespacing=1.55)


# --------------------------------------------------------------------------
# Pending measurements
# --------------------------------------------------------------------------
class PendingMeasurement(RuntimeError):
    """Raised when a panel's measurement has not landed and --strict is set."""


PENDING_LOG: list[str] = []


def pending_panel(ax, *, title: str, marker: str, blocked_on: str,
                  arrives_at: str, predeclaration: str = "", strict: bool = False) -> None:
    """Draw a loud placeholder for a measurement that has not reported.

    Never silently empty: the panel is hatched, carries the marker text the plan
    requires (e.g. ``[MOMENTUM SEED REPLICATION PENDING]``), names the
    measurement it is blocked on and the path its data will arrive at, and is
    recorded in ``PENDING_LOG`` for ``make_all.py`` to surface. With
    ``strict=True`` it raises instead, so the test suite and any automated
    regeneration fail rather than emit a figure with a hole in it.
    """
    msg = f"{marker} {title}: blocked on {blocked_on}; arrives at {arrives_at}"
    PENDING_LOG.append(msg)
    if strict:
        raise PendingMeasurement(msg)
    print(f"PENDING: {msg}", file=sys.stderr)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(True)
        side.set_color(C.VERMILLION)
        side.set_linewidth(1.2)
        side.set_linestyle((0, (4, 2)))
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor="none",
                               edgecolor=C.ENVELOPE, hatch="///", linewidth=0))
    ax.set_title(title)
    body = (f"{marker}\n\n"
            f"blocked on: {blocked_on}\n\n"
            f"data arrives at:\n{arrives_at}")
    if predeclaration:
        body += f"\n\npredeclared in:\n{predeclaration}"
    ax.text(0.5, 0.5, body, transform=ax.transAxes, ha="center", va="center",
            fontsize=6.4, color=C.INK, linespacing=1.7)


# --------------------------------------------------------------------------
# Data access. Everything a figure reads goes through here.
# --------------------------------------------------------------------------
def data_path(rel: str) -> Path:
    p = DATA / rel
    if not p.exists():
        raise FileNotFoundError(
            f"{p} is missing. Run `python extract_from_box.py` from "
            f"{HERE} to refresh the vendored box data.")
    return p


def load_json(rel: str) -> Any:
    return json.loads(data_path(rel).read_text(encoding="utf-8"))


def load_text(rel: str) -> str:
    return data_path(rel).read_text(encoding="utf-8")


def repo_text(rel: str) -> str:
    p = REPO / rel
    if not p.exists():
        raise FileNotFoundError(f"{p} is missing (repository source for a figure value)")
    return p.read_text(encoding="utf-8", errors="replace")


def manifest() -> dict:
    return load_json("MANIFEST.json")


def sha256(rel: str) -> str:
    return hashlib.sha256(data_path(rel).read_bytes()).hexdigest()


# ---- markdown table parsing (repository sources) --------------------------
_NUM = re.compile(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?")


def md_section(text: str, heading: str) -> str:
    """The body of one markdown section, up to the next heading at the same level.

    Several sources in this repository carry more than one table with the same
    shape (`D2_RESULT.md` has an unrestricted, a stratified and a
    contamination table, all `| seed | ... |`), and the figures must not read
    the wrong one. Scoping the parse to a named section makes the choice
    explicit and makes it fail if the section is renamed.
    """
    lines = text.splitlines()
    start = None
    level = 0
    for i, line in enumerate(lines):
        if line.lstrip().startswith("#") and heading.lower() in line.lower():
            start = i + 1
            level = len(line) - len(line.lstrip("#"))
            break
    if start is None:
        raise LookupError(f"no markdown heading containing {heading!r}")
    out = []
    for line in lines[start:]:
        s = line.lstrip()
        if s.startswith("#") and (len(line) - len(line.lstrip("#"))) <= level:
            break
        out.append(line)
    return "\n".join(out)


def md_rows(text: str, *, contains: str = "", n_cols: int | None = None) -> list[list[str]]:
    """Pipe-table rows from a markdown source, filtered by a substring.

    Used so that values quoted in `D2_RESULT.md` etc. are read from the file
    rather than retyped. Returns cells with markdown emphasis stripped. Header
    and `---` separator rows are dropped.
    """
    rows = []
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith("|") or contains not in s:
            continue
        cells = [c.strip().replace("**", "").replace("`", "") for c in s.strip("|").split("|")]
        if all(set(c) <= set("-: ") for c in cells):
            continue
        if n_cols is not None and len(cells) != n_cols:
            continue
        rows.append(cells)
    return rows


def nums(cell: str) -> list[float]:
    """Numbers in a cell. U+2212 MINUS SIGN and en dash are normalised to '-'.

    This matters: the markdown sources are typeset with a real minus sign, and a
    naive parse silently returns the ABSOLUTE value of every negative number in
    this paper - which would have flipped the sign of F1(c)'s channel
    differences and reversed the result.
    """
    return [float(m) for m in _NUM.findall(cell.replace("−", "-").replace("–", "-"))]


def one_num(cell: str) -> float:
    v = nums(cell)
    if len(v) != 1:
        raise ValueError(f"expected exactly one number in {cell!r}, got {v}")
    return v[0]


def ci(cell: str) -> tuple[float, float]:
    """Parse a `[lo, hi]` interval cell, tolerating the en dash / minus sign."""
    v = nums(cell)
    if len(v) != 2:
        raise ValueError(f"expected an interval in {cell!r}, got {v}")
    return v[0], v[1]


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
def save(fig, stem: str) -> list[Path]:
    """Write PDF and SVG (vector) plus a 300 dpi PNG."""
    OUT.mkdir(parents=True, exist_ok=True)
    written = []
    for ext, kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 300})):
        p = OUT / f"{stem}.{ext}"
        fig.savefig(p, bbox_inches="tight", pad_inches=0.06, **kw)
        written.append(p)
    plt.close(fig)
    print(f"wrote {', '.join(p.name for p in written)} -> {OUT}")
    return written


def cli(description: str) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--strict", action="store_true",
                    help="raise PendingMeasurement instead of drawing a pending placeholder")
    ap.add_argument("--out", type=Path, default=None, help="override the output directory")
    args = ap.parse_args()
    if args.out is not None:
        globals()["OUT"] = args.out
    apply_style()
    return args


__all__ = [
    "C", "ORDER", "MARKERS", "LINESTYLES", "HATCHES", "STAT", "BLOCK",
    "axis_label", "apply_style", "grid", "caption", "save", "cli",
    "PendingMeasurement", "pending_panel", "PENDING_LOG",
    "DATA", "OUT", "REPO", "HERE", "data_path", "load_json", "load_text",
    "repo_text", "manifest", "sha256", "md_rows", "nums", "one_num", "ci", "np", "plt",
]
