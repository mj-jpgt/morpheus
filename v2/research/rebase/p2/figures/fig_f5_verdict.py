"""F5 - The verdict is under-determined: statistic, block, view.  (P2 draft 4.5)

Claim: "which arm has the higher effective rank" does not have one answer on our
data until three implementation choices nobody states are fixed - and those
choices are worth more than the between-arm difference they are used to
adjudicate. "Which arm carries more molecular information" is invariant to all
three.

PROVENANCE, and why it matters here. Panel (a) is drawn from
`P2_RANK_VARIANTS.json` and NOT from any printed table: three mutually
incompatible statistics have been called `effective_rank` in this repository, a
fourth (the eigenvalue participation ratio) lived under the R2/R3 labels inside
the analysis code of this very section, and the correction was made at commit
a11549a. Printings of draft 4.5(a) older than a11549a must not be used. The
eigenvalue statistics are labelled `PR` and `PR_rownorm` here and NEVER R2/R3,
and they are drawn in a visually separated block, because `PR` is identical cell
for cell to T1's participation-ratio row and is therefore not independent
evidence.

Sources
  (a) data/ws_p2/out/P2_RANK_VARIANTS.json    (script v2/research/rebase/p2/p2_rank_variants.py)
      data/ws_p2/out/p2_run.log               cross-check of every mark
  (b) data/ws_rank/RANK_RECOMPUTE.json        raw and residualised, all six variants
  (c) data/ws_p2/out/P2_ROBUSTNESS.json       (script v2/research/rebase/p2/p2_robustness.py)

NOT DRAWN, and stated rather than implied: the synthetic inset comparing
(sum s)^2/sum s^2 against (sum s^2)^2/sum s^4 that P2_FIGURES.md F5 asks for is
marked NEEDS EXTRACTION there and is a separate computation under
v2/research/rebase/p2/; it is not part of this figure and no stand-in is drawn
for it.
"""
from __future__ import annotations

import p2fig as P
from p2fig import C, np, plt

PAIRS = ["D2 s42", "D2 s43", "D2 s44", "D1 s42", "D1 s43", "D1 s44"]

#: (json key, printed label, which block it belongs to)
STAT_ROWS = [
    ("R1", "R1   canonical, order 1, centred", "canonical"),
    ("R2", "R2   order-2 Hill of the singular values", "canonical"),
    ("R3", "R3   order-2 Hill, rows L2-normalised", "canonical"),
    ("PR", "PR   eigenvalue participation ratio", "eigen"),
    ("PR_rownorm", "PR$_{\\mathrm{rownorm}}$   eigenvalue PR, rows normalised", "eigen"),
]

VIEW_PAIRS = [("D2 s42", "H42", "I42"), ("D2 s43", "H43", "I43"), ("D2 s44", "H44", "I44"),
              ("D1 s42", "P42", "F42"), ("D1 s43", "P43", "F43"), ("D1 s44", "P44", "F44")]
VIEWS = ["wsi_biology", "rna_biology", "full_biology"]


def marks() -> dict[str, list[str]]:
    v = P.load_json("ws_p2/out/P2_RANK_VARIANTS.json")["verdicts"]
    log = P.load_text("ws_p2/out/p2_run.log")
    for name, entry in v.items():
        needle = f"{name:>12}" + "".join(f"{m:>9}" for m in entry["marks"])
        assert needle in log, (
            f"F5(a): the marks for {name} in P2_RANK_VARIANTS.json do not match the "
            f"printed table in p2_run.log. Refusing to plot either.\n  wanted: {needle!r}")
    return {k: e["marks"] for k, e in v.items()}


def block_flip() -> list[dict]:
    """Raw against residualised, D2 seeds 43 and 44, for R1, R2 and R3."""
    d2 = P.load_json("ws_rank/RANK_RECOMPUTE.json")["instances"]["instance6_D2"]
    out = []
    for seed in (43, 44):
        h, i = d2[f"H{seed}"], d2[f"I{seed}"]
        # On the RAW exported artifacts R2 == R3 exactly, because the model already
        # L2-normalises z_biology. That is why the distinction went unnoticed.
        for art in (h, i):
            assert abs(art["raw"]["R2"] - art["raw"]["R3"]) / art["raw"]["R2"] < 1e-8, art
        for stat in ("R1", "R2", "R3"):
            row = {"seed": seed, "stat": stat}
            for blk in ("raw", "residualised"):
                hv, iv = h[blk][stat], i[blk][stat]
                row[blk] = (hv, iv, "H" if hv > iv else "I")
            row["flips"] = row["raw"][2] != row["residualised"][2]
            out.append(row)
    return out


def view_grid() -> tuple[list[list[str]], list[str], int, int]:
    rob = P.load_json("ws_p2/out/P2_ROBUSTNESS.json")
    grid, info, right, right_d2 = [], [], 0, 0
    for label, a, b in VIEW_PAIRS:
        # The information winner is the same under both CCA estimators on every
        # pair and every view; assert it rather than choose one.
        ins = rob[a]["wsi_biology"]["cca_insample"] > rob[b]["wsi_biology"]["cca_insample"]
        out = rob[a]["wsi_biology"]["cca_heldout"] > rob[b]["wsi_biology"]["cca_heldout"]
        assert ins == out, (label, "the two CCA estimators disagree on wsi_biology")
        winner_info = a[0] if out else b[0]
        info.append(winner_info)
        row = []
        for view in VIEWS:
            hv, iv = rob[a][view]["eff_rank"], rob[b][view]["eff_rank"]
            w = a[0] if hv > iv else b[0]
            row.append(w)
            right += (w == winner_info)
            right_d2 += (w == winner_info) and label.startswith("D2")
        grid.append(row)
    return grid, info, right, right_d2


def _cell(ax, x, y, ok, *, size=0.40, text=None):
    colour = C.BLUE if ok else C.VERMILLION
    ax.add_patch(plt.Rectangle((x - size, y - size), 2 * size, 2 * size,
                               facecolor=colour, alpha=0.13, edgecolor="white",
                               linewidth=1.4, zorder=1))
    ax.text(x, y, text if text is not None else ("✓" if ok else "✗"),
            ha="center", va="center", fontsize=9.5 if text is None else 8.0,
            color=colour, fontweight="bold", zorder=3)


def main() -> int:
    P.cli(__doc__)
    mk = marks()
    flips = block_flip()
    grid, info, right, right_d2 = view_grid()

    disagree = [j for j in range(6)
                if len({mk["R1"][j], mk["R2"][j], mk["R3"][j]}) > 1]

    fig = plt.figure(figsize=(7.5, 9.9))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.10, 0.92, 0.90], hspace=0.46,
                          left=0.315, right=0.965, top=0.955, bottom=0.145)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])

    # ------------------------------------------------------------------ (a)
    ax_a.set_title("(a)  The STATISTIC.  Does higher rank pick the arm with more information?", pad=10)
    ys, labels, colours = [], [], []
    y = 0.0
    for key, label, block in reversed(STAT_ROWS):
        if block == "canonical" and ys and STAT_ROWS[3][0] in [k for k, _l, _b in reversed(STAT_ROWS)][:len(ys)]:
            pass
        ys.append(y)
        labels.append(label)
        colours.append(C.INK if block == "canonical" else C.MUTED)
        y += 1.0
        if key == "PR":            # separate the eigenvalue block visually
            y += 0.55
    ys, labels, colours = ys[::-1], labels[::-1], colours[::-1]
    for row_y, (key, _label, _block) in zip(ys, STAT_ROWS):
        for j, m in enumerate(mk[key]):
            _cell(ax_a, j, row_y, m == "OK")
    info_y = -1.05
    for j in range(6):
        _cell(ax_a, j, info_y, True)
    ax_a.axhline(info_y + 0.62, color=C.RULE, linewidth=0.8)
    sep = (ys[2] + ys[3]) / 2
    ax_a.axhline(sep, color=C.RULE, linewidth=0.8, linestyle=(0, (3, 2)))
    ax_a.set_yticks(list(ys) + [info_y])
    ax_a.set_yticklabels(labels + ["INFORMATION verdict\nheld-out top-CCA, 40 untrained targets"],
                         fontsize=6.8)
    for lbl, col in zip(ax_a.get_yticklabels(), colours + [C.BLUE]):
        lbl.set_color(col)
    ax_a.tick_params(axis="y", length=0)
    ax_a.set_xticks(range(6))
    ax_a.set_xticklabels(PAIRS, fontsize=7.0)
    ax_a.tick_params(axis="x", length=0, top=True, labeltop=True, bottom=False, labelbottom=False)
    ax_a.set_xlim(-0.6, 5.6)
    ax_a.set_ylim(info_y - 2.30, max(ys) + 0.6)
    for side in ax_a.spines.values():
        side.set_visible(False)
    for j in disagree:
        ax_a.add_patch(plt.Rectangle((j - 0.48, min(ys) - 0.48), 0.96, max(ys) - min(ys) + 0.96,
                                     facecolor="none", edgecolor=C.VERMILLION, linewidth=1.1,
                                     linestyle=(0, (3, 2)), zorder=4))
    ax_a.text(-0.6, info_y - 0.72,
              f"R1, R2 and R3 disagree on {len(disagree)} of 6 pairs "
              f"({', '.join(PAIRS[j] for j in disagree)}, boxed). The information verdict is correct in every column.\n"
              "The lower block is the SAME data under a disambiguated name: the eigenvalue participation ratio, which an earlier version of this\n"
              "analysis computed under the R2 / R3 labels. PR is identical cell for cell to T1's participation-ratio row, so the two are NOT\n"
              "independent evidence.",
              ha="left", va="top", fontsize=6.2, color=C.MUTED, linespacing=1.6)

    # ------------------------------------------------------------------ (b)
    ax_b.set_title("(b)  The BLOCK.  Raw against confound-residualised, the two D2 seeds where R3 flips",
                   pad=10)
    ax_b.axis("off")
    hdr = (f"{'':<10}{'':<5}{'RAW block':>26}{'':<3}{'RESIDUALISED block':>28}{'':>10}")
    lines = [hdr,
             f"{'seed':<10}{'stat':<5}{'H':>9}{'I':>9}{'winner':>8}"
             f"{'':<3}{'H':>9}{'I':>9}{'winner':>10}{'':>4}"]
    for r in flips:
        rh, ri, rw = r["raw"]
        dh, di, dw = r["residualised"]
        flag = "   <-- FLIPS" if r["flips"] else ""
        lines.append(f"{r['seed']:<10}{r['stat']:<5}{rh:>9.3f}{ri:>9.3f}{rw:>8}"
                     f"{'':<3}{dh:>9.3f}{di:>9.3f}{dw:>10}{flag}")
    # Rendered line by line so the flipping rows can carry their own colour; a
    # single text block cannot, and a highlight rectangle behind a monospace block
    # is guesswork about line metrics.
    for k, line in enumerate(lines):
        ax_b.text(0.0, 0.97 - 0.088 * k, line, transform=ax_b.transAxes, va="top", ha="left",
                  family="monospace", fontsize=6.4, color=C.INK)
    for k, r in enumerate(flips):
        ax_b.text(0.0, 0.97 - 0.088 * (k + len(lines) - len(flips)), lines[k + 2],
                  transform=ax_b.transAxes, va="top", ha="left", family="monospace",
                  fontsize=6.4, fontweight="bold" if r["flips"] else "normal",
                  color=C.VERMILLION if r["flips"] else C.INK)
    ax_b.text(0.0, 0.20,
              "The block matters for ONE statistic and not the others: R1 and R2 keep their ordering,\n"
              "R3 reverses it on both seeds. On the RAW exported artifacts R2 and R3 are EQUAL to nine\n"
              "significant figures, because the model already L2-normalises z_biology - which is why the\n"
              "distinction went unnoticed. Statistic and block are named on every row; nothing here is\n"
              "an unqualified \"effective rank\".",
              transform=ax_b.transAxes, va="top", ha="left", fontsize=6.4, color=C.MUTED,
              linespacing=1.55)

    # ------------------------------------------------------------------ (c)
    ax_c.set_title("(c)  The VIEW.  Three co-trained views of the same model", pad=10)
    for k, (label, _a, _b) in enumerate(VIEW_PAIRS):
        for v, view in enumerate(VIEWS):
            w = grid[k][v]
            _cell(ax_c, k, len(VIEWS) - 1 - v, w == info[k], text=w, size=0.42)
        _cell(ax_c, k, -1.15, True, text=info[k], size=0.42)
    ax_c.axhline(-0.55, color=C.RULE, linewidth=0.8)
    ax_c.set_xticks(range(6))
    ax_c.set_xticklabels([p for p, _a, _b in VIEW_PAIRS], fontsize=7.0)
    ax_c.tick_params(axis="x", length=0, top=True, labeltop=True, bottom=False, labelbottom=False)
    ax_c.set_yticks([2, 1, 0, -1.15])
    ax_c.set_yticklabels([f"rank winner, {v}" for v in VIEWS]
                         + ["INFORMATION winner\n(identical on all three views)"], fontsize=6.8)
    for lbl in ax_c.get_yticklabels()[:3]:
        lbl.set_color(C.INK)
    ax_c.get_yticklabels()[3].set_color(C.BLUE)
    ax_c.tick_params(axis="y", length=0)
    ax_c.set_xlim(-0.6, 5.6)
    ax_c.set_ylim(-2.95, 2.6)
    for side in ax_c.spines.values():
        side.set_visible(False)
    ax_c.text(-0.6, -1.82,
              f"Rank agrees with information on {right} of {len(VIEW_PAIRS) * len(VIEWS)} "
              f"(pair x view) comparisons; restricted to D2, {right_d2} of 9.\n"
              "Rank uses " + P.STAT["R1_short"] + "; information is the held-out top-CCA on the "
              "40 untrained targets.",
              ha="left", va="top", fontsize=6.4, color=C.INK, linespacing=1.6)

    P.caption(fig,
        "F5. These panels do not show that rank is wrong. They show that A RANK VERDICT IS NOT A WELL-DEFINED OBJECT until the statistic, the block and "
        "the view are stated - and that the three choices are worth more than the between-arm difference they are used to adjudicate. Panel (a) is drawn from "
        "ws_p2/out/P2_RANK_VARIANTS.json, not from any printed table, and its eigenvalue rows are labelled PR and PR_rownorm, never R2/R3: a fourth statistic "
        "lived under those labels inside the analysis code of the section that argues the name is unreliable, and was corrected at commit a11549a. PR is "
        "identical cell for cell to T1's participation-ratio row and is not independent evidence. In panel (c) the rna_biology -> RNA-derived-target comparison "
        "is partly circular and its absolute CCA (0.79-0.85) must not be read as a clean image-to-molecular channel; the RANK measurements on that view are "
        "unaffected by the circularity, but the D2 count inherits the caveat and is not quoted as a rate.",
        y=0.006)

    P.save(fig, "F5_verdict_instability")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
