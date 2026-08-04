"""T1 - Rank as a selection rule against the published alternatives, AND IT IS
UNDERPOWERED, AND ITS GROUND TRUTH IS BLOCK-DEPENDENT.  (P2 draft 4.6, 4.6a)

Twelve metric rows x six matched pairs. The design can detect a perfect rule and
nothing weaker, and THE UNDERPOWERING IS PART OF THE TABLE: it is printed inside
the table body, as a banded row above the metrics, not relegated to a caption a
reader may not have. D2 s44's column is marked unresolvable in every row, because
F4(a) puts that pair's between-arm gap at 1.4 sampling standard deviations - so
effective rank's honest D2 record reads "1 clear hit, 1 clear miss, 1 pair it
cannot resolve".

A SECOND BAND, below the table, carries the caveat that arrived on 2026-08-04 and
is larger than the first. Every mark above is scored against ONE ground truth: the
held-out channel onto 40 gene-set targets. That channel's arm contrast turns out to
be a property of the coordinate system the targets are written in
(NOTEBOOK_ENTRIES/d2_coordinate_system_result_20260804T0800Z.md 1a): the -0.12
appears on the gene sets and on none of the five other target blocks on disk. Re-run
the same rule with each block as truth in turn and seed 43's arm ordering flips, so
every metric row's D2 count moves - and the ordering between the two rows this table
is quoted for REVERSES. Drawing the table without that band would present a
coordinate choice as a fact.

The exact two-sided binomial p is RECOMPUTED here from math.comb rather than
copied, and the marks are recomputed from the per-artifact metrics and asserted
against the printed table before anything is drawn. The band's counts are recomputed
the same way, from the same metrics, against the per-block contrasts.

Rendered as a figure (PDF / SVG / PNG) so that the table travels with its own
banded caveats rather than depending on a typesetting step.

Sources
  data/ws_p2/out/p2_run.log            the printed selection-rule table
  data/ws_p2/out/P2_METRICS_D2.json    the per-artifact metrics the marks come from
  data/ws_p2/out/P2_METRICS_D1.json
  v2/research/rebase/nature/d2_coordinate_system/out/EXAM_PANEL.json
                                       the arm contrast on each of the six blocks
  scripts v2/research/rebase/p2/p2_competing_metrics.py, p2_selection_rule.py,
          p2_selection_rule_blocks.py
"""
from __future__ import annotations

import json
import math
import re

import p2fig as P
from p2fig import C, np, plt

PAIRS = [("D2 s42", "d2", "H42", "I42"), ("D2 s43", "d2", "H43", "I43"),
         ("D2 s44", "d2", "H44", "I44"), ("D1 s42", "d1", "P42", "F42"),
         ("D1 s43", "d1", "P43", "F43"), ("D1 s44", "d1", "P44", "F44")]

#: printed row label -> (metric key in the JSON, higher-is-better?)
#: alpha-ReQ and LiDAR rows are scored by the run script from quantities the
#: metrics JSON does not expose in the same form, so they are read from the
#: printed table and are not independently recomputed here; that is stated in
#: the row note rather than glossed over.
RECOMPUTABLE = {
    "effective_rank_raw": ("effective_rank_raw", True),
    "effective_rank_residualised": ("effective_rank_residualised", True),
    "rankme_raw (as published)": ("rankme_raw", True),
    "rankme_residualised": ("rankme_residualised", True),
    "participation_ratio_raw": ("participation_ratio_raw", True),
    "participation_ratio_residualised": ("participation_ratio_residualised", True),
    "stable_rank_raw": ("stable_rank_raw", True),
    "stable_rank_residualised": ("stable_rank_residualised", True),
}

#: how each row must be described, per the plan's "Must carry".
ROW_NOTE = {
    "alpha_ReQ_raw |alpha-1|": "authors' released fastssl estimator, not the paper text",
    "alpha_ReQ_resid |alpha-1|": "|alpha-1| is OUR operationalisation of their \"Goldilocks zone\"",
    "LiDAR_raw": "ADAPTED: q = 2, the two modalities as views",
    "LiDAR_residualised": "outside the authors' tested q range; ordering invariant over delta 1e-8 to 1e0",
}


#: The six target blocks the D2 arms were re-scored on, shortest label first.
#: Order follows §1a of the coordinate-system entry: the published block, then the
#: four other 128-column expression-derived blocks, then the negative control.
EXAM_BLOCKS = [
    ("geneset_untrained40", "gene sets 40\n(PUBLISHED)"),
    ("pbs_codes128_ARM_I_OWN", "PBS codes\n128"),
    ("pca_basis128", "PCA basis\n128"),
    ("pbs_shuffled_s1", "shuffled PBS\n128"),
    ("geneset_random_control", "rand control\n90 (NEG CTL)"),
    ("random_dictionary128", "rand dict\n128"),
]

#: The two rows §4.6 quotes against one another. If their ordering is not stable
#: across blocks then neither is the sentence the section is built on.
BAND_ROWS = [("effective_rank_residualised", "canonical effective rank (resid.)"),
             ("rankme_raw (as published)", "RankMe (raw, as published)")]


def per_block_counts() -> tuple[dict, dict, list[int]]:
    """D2 count per metric per block, recomputed; plus each block's winner pattern.

    The metrics are held fixed and only the truth is swapped, which is the whole
    content of the counterfactual. Δ in the panel is PBS − Hallmark, so Δ < 0 means
    the Hallmark arm carries more; that sign is asserted against the published
    gene-set row rather than assumed.
    """
    panel = json.loads(P.repo_text(
        "v2/research/rebase/nature/d2_coordinate_system/out/EXAM_PANEL.json"))
    truth: dict[str, dict[int, float]] = {}
    for row in panel:
        assert row["block"] == "residualised", row
        truth.setdefault(row["exam"], {})[int(row["seed"])] = float(row["pbs_minus_hallmark"])

    d2 = P.load_json("ws_p2/out/P2_METRICS_D2.json")
    # The published truth must agree with the metrics JSON's own untrained-40
    # channel, or the panel and the table are describing different runs.
    for seed in (42, 43, 44):
        published = (d2[f"I{seed}"]["points"]["untrained40"]["top_cca"]
                     - d2[f"H{seed}"]["points"]["untrained40"]["top_cca"])
        got = truth["geneset_untrained40"][seed]
        assert abs(published - got) < 5e-4, (
            f"T1: EXAM_PANEL.json gives {got:+.4f} for the gene-set arm contrast on seed "
            f"{seed}, but P2_METRICS_D2.json gives {published:+.4f}. Refusing to draw either.")

    counts: dict[str, dict[str, int]] = {}
    patterns: dict[str, str] = {}
    for key, _label in EXAM_BLOCKS:
        patterns[key] = "".join("H" if truth[key][s] < 0 else "I" for s in (42, 43, 44))
        counts[key] = {}
        for metric, (json_key, higher_better) in RECOMPUTABLE.items():
            ok = 0
            for seed in (42, 43, 44):
                ma = d2[f"H{seed}"]["metrics"][json_key]
                mb = d2[f"I{seed}"]["metrics"][json_key]
                picked_a = (ma > mb) if higher_better else (ma < mb)
                ok += picked_a == (truth[key][seed] < 0)
            counts[key][metric] = ok
    flips = [s for s in (42, 43, 44)
             if len({truth[k][s] < 0 for k, _l in EXAM_BLOCKS}) > 1]
    return counts, patterns, flips


def exact_two_sided_binomial(k: int, n: int) -> float:
    """Exact two-sided binomial p against p0 = 1/2, by the point-probability rule."""
    pmf = [math.comb(n, i) / 2 ** n for i in range(n + 1)]
    return sum(p for p in pmf if p <= pmf[k] + 1e-12)


def printed_table() -> list[dict]:
    log = P.load_text("ws_p2/out/p2_run.log")
    block = log.split("SELECTION-RULE SCORE")[1].split("WITHIN-ARM SEED VARIABILITY")[0]
    rows = []
    for m in re.finditer(
            r"^(\S.*?)\s{2,}((?:(?:OK|MISS)\s+){6})(\d)/3\s+(\d)/3\s+(\d)/6\s+([\d.]+)\s*$",
            block, re.M):
        rows.append({"metric": m.group(1).strip(),
                     "marks": m.group(2).split(),
                     "d2": int(m.group(3)), "d1": int(m.group(4)),
                     "all": int(m.group(5)), "p": float(m.group(6))})
    assert len(rows) == 12, [r["metric"] for r in rows]
    return rows


def recompute(rows: list[dict]) -> None:
    """Recompute every mark the metrics JSON can support, and assert it."""
    src = {"d2": P.load_json("ws_p2/out/P2_METRICS_D2.json"),
           "d1": P.load_json("ws_p2/out/P2_METRICS_D1.json")}
    for r in rows:
        spec = RECOMPUTABLE.get(r["metric"])
        if spec is None:
            continue
        key, higher_better = spec
        for j, (label, which, a, b) in enumerate(PAIRS):
            ma, mb = src[which][a]["metrics"][key], src[which][b]["metrics"][key]
            ca = src[which][a]["points"]["untrained40"]["top_cca"]
            cb = src[which][b]["points"]["untrained40"]["top_cca"]
            picked_a = (ma > mb) if higher_better else (ma < mb)
            got = "OK" if picked_a == (ca > cb) else "MISS"
            assert got == r["marks"][j], (
                f"T1: recomputing {r['metric']} on {label} from the metrics JSON gives "
                f"{got}, but data/ws_p2/out/p2_run.log prints {r['marks'][j]}. "
                "Refusing to draw either.")
        # And the exact binomial p, recomputed rather than copied.
        want = exact_two_sided_binomial(r["all"], 6)
        assert abs(want - r["p"]) < 6e-4, (
            f"T1: exact two-sided binomial for {r['all']}/6 is {want:.4f}, but the "
            f"run log prints {r['p']:.3f} for {r['metric']}.")


def main() -> int:
    P.cli(__doc__)
    rows = printed_table()
    recompute(rows)
    block_counts, patterns, flips = per_block_counts()

    # F4(a): the pair whose between-arm gap is smaller than three sampling sd.
    unresolvable = 2   # index of D2 s44 in PAIRS
    d2 = P.load_json("ws_p2/out/P2_METRICS_D2.json")
    sa = d2["H44"]["subsample"]["effective_rank_residualised"]
    sb = d2["I44"]["subsample"]["effective_rank_residualised"]
    gap_sd = abs(sa["mean"] - sb["mean"]) / math.hypot(sa["sd"], sb["sd"])
    assert gap_sd < 3.0, gap_sd

    n_rows = len(rows)
    # The table keeps the geometry it was designed with (7.6 x 6.7 in); the figure
    # grows underneath it so the block-dependence band has its own space rather
    # than being squeezed into the footnotes.
    fig = plt.figure(figsize=(7.6, 8.75))
    ax = fig.add_axes([0.0, 0.235, 1.0, 0.765])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    x_metric = 0.015
    x_pair = [0.345 + 0.070 * i for i in range(6)]
    x_tot = [0.770, 0.822, 0.874]
    x_p = 0.985
    top = 0.905
    row_h = 0.0455

    ax.text(x_metric, 0.975,
            "T1   Rank as a selection rule against the published alternatives - "
            "and it is UNDERPOWERED, and its GROUND TRUTH MOVES",
            ha="left", va="bottom", fontsize=9.0, fontweight="bold", color=C.INK)
    ax.text(x_metric, 0.955,
            "Does the metric pick the arm that carries more molecular information? "
            "Ground truth: held-out top-CCA on the 40 GENE-SET targets neither arm trained "
            "on, confound-residualised block - one of the six blocks banded below.",
            ha="left", va="bottom", fontsize=6.6, color=C.MUTED)

    # ---- the underpowering, INSIDE the table body -------------------------
    band_h = 0.075
    ax.add_patch(plt.Rectangle((x_metric - 0.008, top - band_h + 0.006), 1.0 - x_metric,
                               band_h, facecolor=C.VERMILLION, alpha=0.10,
                               edgecolor=C.VERMILLION, linewidth=1.0, zorder=1))
    p6, p5, p4 = (exact_two_sided_binomial(k, 6) for k in (6, 5, 4))
    ax.text(x_metric + 0.004, top - 0.006,
            f"n = 6.  Exact two-sided binomial:  6/6 -> p = {p6:.3f};   5/6 -> {p5:.3f};   "
            f"4/6 -> {p4:.3f}.\n"
            "THIS DESIGN CAN DETECT A PERFECT RULE AND NOTHING WEAKER.  "
            "No comparison between two rows of this table is evidence in either direction.",
            ha="left", va="top", fontsize=7.0, color=C.INK, fontweight="bold",
            linespacing=1.7, zorder=3)

    head = top - band_h - 0.018
    for label, x in zip([p[0] for p in PAIRS], x_pair):
        ax.text(x, head, label, ha="center", va="bottom", fontsize=6.8, color=C.INK,
                fontweight="bold")
    for label, x in zip(["D2", "D1", "ALL"], x_tot):
        ax.text(x, head, label, ha="center", va="bottom", fontsize=6.8, color=C.INK,
                fontweight="bold")
    ax.text(x_p, head, "exact\ntwo-sided $p$", ha="right", va="bottom", fontsize=6.4,
            color=C.INK, fontweight="bold", linespacing=1.4)
    ax.text(x_metric, head, "metric", ha="left", va="bottom", fontsize=6.8,
            color=C.INK, fontweight="bold")

    # The unresolvable column, marked in EVERY row.
    ax.add_patch(plt.Rectangle((x_pair[unresolvable] - 0.031, head - row_h * n_rows - 0.012),
                               0.062, row_h * n_rows + 0.020, facecolor="none",
                               edgecolor=C.VERMILLION, linewidth=1.1, linestyle=(0, (3, 2)),
                               hatch="///", zorder=0, alpha=0.55))

    for k, r in enumerate(rows):
        y = head - 0.020 - row_h * k
        if k % 2 == 0:
            ax.add_patch(plt.Rectangle((x_metric - 0.008, y - row_h * 0.32), 1.0 - x_metric,
                                       row_h * 0.92, facecolor=C.RULE, alpha=0.28,
                                       edgecolor="none", zorder=0))
        emph = r["metric"].startswith("effective_rank")
        ax.text(x_metric, y, r["metric"], ha="left", va="center", fontsize=6.8,
                color=C.INK, fontweight="bold" if emph else "normal", zorder=3)
        note = ROW_NOTE.get(r["metric"])
        if note:
            ax.text(x_metric + 0.005, y - row_h * 0.34, note, ha="left", va="center",
                    fontsize=5.3, color=C.MUTED, style="italic", zorder=3)
        for j, (mark, x) in enumerate(zip(r["marks"], x_pair)):
            if j == unresolvable:
                ax.text(x, y, "unres.", ha="center", va="center", fontsize=5.8,
                        color=C.VERMILLION, fontweight="bold", zorder=3)
                continue
            ok = mark == "OK"
            ax.text(x, y, "✓" if ok else "✗", ha="center", va="center", fontsize=9.0,
                    color=C.BLUE if ok else C.VERMILLION, fontweight="bold", zorder=3)
        for val, x in zip((f"{r['d2']}/3", f"{r['d1']}/3", f"{r['all']}/6"), x_tot):
            ax.text(x, y, val, ha="center", va="center", fontsize=6.8, color=C.INK,
                    fontweight="bold" if x == x_tot[2] else "normal", zorder=3)
        ax.text(x_p, y, f"{r['p']:.3f}", ha="right", va="center", fontsize=6.8,
                color=C.INK if r["p"] < 0.05 else C.MUTED, zorder=3)

    foot = head - 0.020 - row_h * n_rows - 0.006
    ax.text(x_metric, foot,
            f"D2 s44 IS MARKED UNRESOLVABLE IN EVERY ROW. F4(a) puts that pair's between-arm gap at "
            f"{gap_sd:.1f} sampling standard deviations under 40 draws of 80% of\n"
            "patients, so no row of this table can score it. Effective rank's honest D2 record therefore "
            "reads: 1 clear hit, 1 clear miss, 1 pair it cannot resolve.\n"
            "The counts in the D2, D1 and ALL columns are the run's own, which score that cell; they are "
            "printed as the run produced them and are NOT corrected here.",
            ha="left", va="top", fontsize=6.3, color=C.VERMILLION, linespacing=1.7)

    # ---- the block-dependence band, in its own axes ------------------------
    axb = fig.add_axes([0.0, 0.0, 1.0, 0.235])
    axb.set_xlim(0, 1)
    axb.set_ylim(0, 1)
    axb.axis("off")
    axb.add_patch(plt.Rectangle((x_metric - 0.008, 0.020), 1.0 - x_metric, 0.945,
                                facecolor=C.ORANGE, alpha=0.11, edgecolor=C.ORANGE,
                                linewidth=1.0, zorder=1))
    axb.text(x_metric + 0.004, 0.925,
             "AND EVERY MARK ABOVE IS SCORED AGAINST ONE COORDINATE SYSTEM OUT OF SIX.",
             ha="left", va="top", fontsize=7.4, color=C.INK, fontweight="bold", zorder=3)
    axb.text(x_metric + 0.004, 0.860,
             "The ground truth is the held-out channel onto 40 GENE-SET targets. That arm "
             "contrast exists on the gene sets and on no other target block on disk "
             "(draft 4.6a). Re-scoring the same two arms on each block\nin turn flips seed "
             f"{flips[0] if flips else '-'}, and with it every metric row's D2 count. "
             "The D1 half cannot be re-scored -- those arms were never scored on any block "
             "but the gene sets -- so only the D2 count moves here.",
             ha="left", va="top", fontsize=6.3, color=C.INK, linespacing=1.6, zorder=3)

    x_blk = [0.345 + 0.108 * i for i in range(len(EXAM_BLOCKS))]
    y_hdr = 0.660
    for (key, label), x in zip(EXAM_BLOCKS, x_blk):
        axb.text(x, y_hdr, label, ha="center", va="top", fontsize=5.9, color=C.INK,
                 fontweight="bold", linespacing=1.35, zorder=3)
    band_rows = [("arm ordering, seeds 42 / 43 / 44",
                  [patterns[k] for k, _l in EXAM_BLOCKS], None)]
    for metric, label in BAND_ROWS:
        d1 = next(r["d1"] for r in rows if r["metric"] == metric)
        band_rows.append((f"{label}   D2 count", [block_counts[k][metric] for k, _l in EXAM_BLOCKS],
                          d1))
    for j, (label, values, d1) in enumerate(band_rows):
        y = 0.505 - 0.098 * j
        axb.text(x_metric + 0.004, y, label, ha="left", va="center", fontsize=6.4,
                 color=C.INK, zorder=3)
        for v, x in zip(values, x_blk):
            if d1 is None:
                axb.text(x, y, str(v), ha="center", va="center", fontsize=7.0,
                         color=C.VERMILLION if "I" in str(v) else C.INK,
                         fontweight="bold", zorder=3)
            else:
                axb.text(x, y, f"{v}/3", ha="center", va="center", fontsize=7.0,
                         color=C.INK, fontweight="bold", zorder=3)
                axb.text(x, y - 0.042, f"({v + d1}/6, p={exact_two_sided_binomial(v + d1, 6):.3f})",
                         ha="center", va="center", fontsize=5.4,
                         color=C.VERMILLION if v + d1 == 6 else C.MUTED, zorder=3)
    er_key, rm_key = BAND_ROWS[0][0], BAND_ROWS[1][0]
    reversals = sum(
        (block_counts[k][er_key] > block_counts[k][rm_key])
        != (block_counts["geneset_untrained40"][er_key]
            > block_counts["geneset_untrained40"][rm_key])
        for k, _l in EXAM_BLOCKS)
    axb.text(x_metric + 0.004, 0.220,
             f"THE ORDERING BETWEEN THESE TWO ROWS REVERSES ON {reversals} OF "
             f"{len(EXAM_BLOCKS)} BLOCKS. On the published gene sets RankMe is ahead of "
             "canonical effective rank on D2; on the dictionary's own\ncodes and on a plain "
             "PCA basis the two swap, and canonical effective rank reaches 6/6 overall - "
             "p = 0.031, 'significant' by this table's own bar, produced by nothing but the "
             "choice\nof which coordinate system the exam is written in. No count in this "
             "table may be quoted without its target block.",
             ha="left", va="top", fontsize=6.3, color=C.INK, linespacing=1.6, zorder=3)

    ax.text(x_metric, foot - 0.095,
            "LiDAR is ADAPTED - q = 2 with the two modalities as views, licensed by the paper's own "
            "footnote 4 but outside the authors' tested q range - and its\n"
            "ordering is invariant across delta from 1e-8 to 1e0. alpha-ReQ follows the authors' "
            "released fastssl estimator rather than the paper text, and |alpha - 1| is\n"
            "OUR operationalisation of their \"Goldilocks zone\", not theirs: all twelve artifacts sit "
            "at alpha between 2.6 and 4.8, far outside it. RankMe as published is\n"
            "uncentred on the RAW block, per its own definition; no value in this table is compared "
            "against a PUBLISHED RankMe number.",
            ha="left", va="top", fontsize=6.3, color=C.MUTED, linespacing=1.7)

    P.save(fig, "T1_selection_rule")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
