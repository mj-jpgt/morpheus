"""F9 - Rank rises while a co-measured collapse measure rises with it.  (P2 draft 4.9a)

Claim: raising the covariance-decorrelation weight raises effective rank
monotonically across three levels AND raises the RNA-view patient-to-patient
mutual cosine monotonically at the same time - on the SAME three runs, printed
on the same log lines. Rank reports more occupied directions; a direct
measurement of the condition rank exists to detect reports the patients' states
converging on one vector. That is a stronger dissociation than anything in draft
4.9, for two reasons about the SHAPE of the evidence rather than its size: it is
monotone across three levels rather than a single contrast, and the
contradicting quantity is CO-MEASURED rather than inferred from a downstream
readout in another table.

WHAT THE PANEL MUST NOT BE READ AS. The rank change is 1.854x under R3 and
1.940x under the canonical statistic. Both are inside draft 4.1's 3.295x figure -
AND 3.295x NEVER LICENSED THEM, because it is canonical R1 on the residualised
EXPORTED wsi_biology block of a DIFFERENT ARM at 40 epochs. The band drawn on (a)
and (b) is therefore this block's own: five same-seed repeats per arm on the
fixed held-out probe at STEP 400 give 1.4489x under R3 and 1.5702x under
canonical R1, and against those the sweep CLEARS. Clearing is not the claim:
at ONE SEED PER LEVEL a magnitude that clears an n = 5 floor is not evidence of
a dose-response. The monotonicity and the co-measured cosine carry this
observation; the magnitude of the rank change does not, before or after the floor
existed. Both statements are annotated inside the artwork, not left to a caption
a reader may not have.

CHANGED 2026-08-05. This script previously drew 3.295x and printed "BOTH INSIDE",
which is the block mismatch draft 4.1a exists to catch, applied to ourselves. See
NOTEBOOK_ENTRIES/the_probe_block_has_a_floor_at_last_20260804T1620Z.md and draft
4.1a rows 54-55.

TWO RANK STATISTICS, TWO PANELS. Binding constraint 1 of paper/P2_FIGURES.md
forbids putting two rank statistics on one axis, so (a) carries R3 - the column
these logs' own `final_eff_rank=` line reports, and therefore the column a
notebook entry quoting "eff-rank" is quoting - and (b) carries the canonical
Roy & Vetterli order-1 statistic, which is the one draft 4.1's floor is measured
in. The cosine series is identical in both; only the rank statistic changes, and
the direction survives it.

THE FLOOR IS ON A THIRD BLOCK AGAIN. These runs are read on a fixed 256-patient
held-out probe. Draft 4.1's floor is canonical R1 on the residualised EXPORTED
`wsi_biology` block, measured on `programme_only` at 40 epochs. The band drawn
here is therefore indicative and is labelled as such - no floor has been
measured on the probe. See draft 4.1a rows 48-50.

Sources
  data/e0_run/d1_diag/ablate_decorr{0.0,0.01,0.04}.log
      one seed per level, m = 0.999, capacity 4096, lr 2e-4, 400 steps, produced
      by v2/research/rebase/d1_momentum_probe.py, which imports BOTH rank
      statistics from v2/calibra and computes neither inline.
  data/e0_run/d1_diag/mseed_m0.999_s42.log
      an INDEPENDENT SAME-SEED REPEAT of the decorrelation = 0.04 arm: identical
      momentum, decorrelation, capacity, learning rate, seed and step-0 state,
      and the probe runs a constant learning rate with no schedule keyed to the
      step budget. n = 2 is a PAIR, not a floor (draft 4.1 says so in terms), and
      the panel labels it that way.
  data/e0_run/d1_probefloor/out/P2_PROBE_FLOORS.json
      THIS BLOCK'S OWN floor, per statistic, at the reading step these runs are
      read at (400): ten runs, five identical repeats of each of two arms, on the
      same fixed held-out probe. Measured by
      v2/research/rebase/p2/p2_probe_floors.py; nothing recomputed here.
"""
from __future__ import annotations

import p2fig as P
from p2fig import C, np, plt

DIAG = "e0_run/d1_diag"
LEVELS = ["0.0", "0.01", "0.04"]
#: The two headers `d1_momentum_probe.py` has ever written. A log carrying any
#: other one fails here rather than yielding whichever column sits at that index.
HEADERS = {
    ("R3-rank", "CANONICAL", "feat-std", "rna-rna", "contrastive"),
    ("eff-rank", "feat-std", "rna-rna", "contrastive"),
}
READ_AT = 400

LEVEL_STYLE = [
    (C.SKY, "o", P.LINESTYLES[0]),
    (C.BLUE, "s", P.LINESTYLES[1]),
    (C.PURPLE, "^", P.LINESTYLES[2]),
]
COSINE_LABEL = ("RNA-view patient-to-patient mutual cosine\n"
                "mean off-diagonal cosine of L2-normalised rows\n"
                "block: " + P.BLOCK["probe"] + "     1.0 = total collapse")


def probe_log(rel: str) -> tuple[dict, dict[int, dict[str, float]]]:
    """Parse one probe log into its header fields and its per-step rows.

    Nothing is computed here: every value is a column of the log, which was
    written by a script that imports its rank statistics from `v2/calibra`.
    """
    text = P.load_text(f"{DIAG}/{rel}")
    meta: dict[str, str] = {}
    header: tuple[str, ...] | None = None
    rows: dict[int, dict[str, float]] = {}
    for line in text.splitlines():
        if line.startswith("momentum="):
            for token in line.split():
                key, _, value = token.partition("=")
                meta[key] = value
            continue
        cells = line.split()
        if cells[:1] == ["step"]:
            header = tuple(cells[1:])
            assert header in HEADERS, f"{rel}: unexpected header {header!r}"
            continue
        if header is None or not cells or not cells[0].isdigit():
            continue
        if len(cells) - 1 != len(header):
            continue
        row = {}
        for name, cell in zip(header, cells[1:]):
            try:
                row[name] = float(cell)
            except ValueError:                    # `nan` at step 0
                row[name] = float("nan")
        rows[int(cells[0])] = row
    assert rows, rel
    return meta, rows


def ablation() -> dict:
    """The three levels, with the guards that make them one ablation.

    The three arms must share a single verified initialisation, or the sweep
    varies more than the decorrelation weight and the monotonicity means
    nothing. They must also share momentum, capacity, learning rate and seed.
    """
    out = {}
    for level in LEVELS:
        meta, rows = probe_log(f"ablate_decorr{level}.log")
        assert meta["decorrelation"] == level, (level, meta)
        out[level] = {"meta": meta, "rows": rows}
    ref = out[LEVELS[0]]
    for level in LEVELS[1:]:
        for key in ("momentum", "capacity", "lr", "seed", "steps", "probe"):
            assert out[level]["meta"][key] == ref["meta"][key], (level, key)
        for column in ("R3-rank", "CANONICAL", "feat-std", "rna-rna"):
            assert out[level]["rows"][0][column] == ref["rows"][0][column], (
                f"F9: the decorrelation = {level} arm does not start from the same state as the "
                f"{LEVELS[0]} arm on {column}. Refusing to draw a sweep that varies more than "
                "the decorrelation weight.")
    return out


def repeat_of_the_top_level() -> tuple[dict, dict[int, dict[str, float]]]:
    """The independent same-seed repeat of the decorrelation = 0.04 arm."""
    meta, rows = probe_log("mseed_m0.999_s42.log")
    ref = probe_log("ablate_decorr0.04.log")[0]
    for key in ("momentum", "decorrelation", "capacity", "lr", "seed", "probe"):
        assert meta[key] == ref[key], (
            f"F9: mseed_m0.999_s42 differs from ablate_decorr0.04 on {key} "
            f"({meta[key]!r} against {ref[key]!r}) and is therefore not a same-seed repeat.")
    return meta, rows


#: The reading step these three runs are quoted at, and therefore the only step
#: whose probe floor may be drawn beside them. Written here rather than passed in
#: so that a change of reading step has to change the floor too.
FLOOR_STEP = "400"
FLOOR_VIEW = "wsi_biology"


def floor(statistic: str) -> float:
    """THIS block's own same-seed retraining floor, at THIS reading step.

    Not draft 4.1's 3.295x. That is canonical R1 on the residualised EXPORTED
    `wsi_biology` block of a different arm at 40 epochs, and these runs are
    `programme_free` on the fixed held-out probe - a different block, so it never
    licensed this comparison. Drawing it here is exactly the mismatch draft 4.1a
    was built to find, and this script used to commit it.

    The floor is read per statistic, because draft 4.1a's own measurement is that
    the floor is a property of the statistic (1.000x to 3.295x on one block), so
    judging the R3 panel against R1's floor would be the second half of the same
    error.
    """
    blob = P.load_json("e0_run/d1_probefloor/out/P2_PROBE_FLOORS.json")
    cell = blob["floors"][FLOOR_STEP][FLOOR_VIEW][statistic]
    fold = cell["floor"]
    recomputed = cell["floor_max"] / cell["floor_min"]
    assert abs(recomputed - fold) < 0.002, (
        f"F9: P2_PROBE_FLOORS records a {statistic} floor of {fold} at step {FLOOR_STEP} but its "
        f"own min and max give {recomputed:.4f}. Refusing to draw either.")
    return fold


def _dissociation_panel(ax, arms, repeat, column, stat_key, floor_fold, ticks, *, letter, title):
    """One twin-axis panel: a rank statistic on the left, the cosine on the right.

    Every long sentence lives in the text band beneath the panels; only the two
    series, their values and the floor band are inside the axes, because this is
    the panel whose whole job is to be legible at a glance.
    """
    x = np.arange(len(LEVELS), dtype=float)
    rank = [arms[l]["rows"][READ_AT][column] for l in LEVELS]
    cosine = [arms[l]["rows"][READ_AT]["rna-rna"] for l in LEVELS]
    fold = max(rank) / min(rank)

    # The floor, anchored at this panel's own decorrelation = 0 value. Drawn
    # first, and as a BAND rather than a line, because what it bounds is a ratio.
    ax.axhspan(rank[0], rank[0] * floor_fold, color=C.ENVELOPE, alpha=0.38,
               linewidth=0, zorder=0)
    verdict = "clears it" if fold > floor_fold else "is inside it"
    ax.annotate(f"THIS BLOCK'S OWN same-seed retraining floor, ×{floor_fold:.3f} — five repeats per "
                f"arm on the\nfixed held-out probe at step {READ_AT}, anchored at this panel's own "
                f"decorrelation = 0 value:\n{rank[0]:.2f} → {rank[0] * floor_fold:.2f}."
                f"   The sweep's ×{fold:.3f} {verdict}   —   and at one seed per level "
                f"that carries nothing.",
                (len(LEVELS) - 1 + 0.44, rank[0]), xytext=(0, 4),
                textcoords="offset points", ha="right", va="bottom", fontsize=6.3,
                color=C.MUTED, linespacing=1.5, zorder=1)

    line_rank, = ax.plot(x, rank, color=C.BLUE, marker="o", linestyle=P.LINESTYLES[0],
                         markersize=7.5, markeredgecolor="white", markeredgewidth=1.1,
                         zorder=4, label="effective rank  —  rises,  reads as “better”")
    for xi, v in zip(x, rank):
        ax.annotate(f"{v:.2f}", (xi, v), xytext=(-8, 9), textcoords="offset points",
                    ha="right", va="bottom", fontsize=7.4, color=C.BLUE, fontweight="bold")

    ax.set_yscale("log")
    # The top of the axis has to clear BOTH the band and the series. When the
    # floor was 3.295x the band always won; on this block's own floor the sweep
    # rises above it, which is the point, so the limit is taken over both.
    lo, hi = min(rank) * 0.66, max(rank[0] * floor_fold, max(rank)) * 1.34
    ax.set_ylim(lo, hi)
    ax.set_ylabel(P.axis_label(stat_key, "probe"), fontsize=6.2, color=C.BLUE)
    ax.tick_params(axis="y", colors=C.BLUE)
    ax.minorticks_off()
    shown = [t for t in ticks if lo <= t <= hi]
    ax.set_yticks(shown)
    ax.set_yticklabels([str(t) for t in shown])
    ax.set_xticks(x)
    ax.set_xticklabels(LEVELS, fontsize=8.2)
    ax.set_xlim(-0.46, len(LEVELS) - 1 + 0.46)
    ax.set_xlabel("`feature_decorrelation` weight      ONE SEED PER LEVEL", fontsize=7.0,
                  labelpad=5)
    P.grid(ax, axis="x")

    twin = ax.twinx()
    twin.spines["right"].set_visible(True)
    twin.spines["right"].set_color(C.RULE)
    line_cos, = twin.plot(x, cosine, color=C.VERMILLION, marker="s",
                          linestyle=P.LINESTYLES[1], markersize=7.5,
                          markeredgecolor="white", markeredgewidth=1.1, zorder=4,
                          label="RNA-view mutual cosine  —  rises,  means WORSE")
    for xi, v in zip(x, cosine):
        twin.annotate(f"{v:.4f}", (xi, v), xytext=(9, -6), textcoords="offset points",
                      ha="left", va="top", fontsize=7.4, color=C.VERMILLION,
                      fontweight="bold")
    twin.set_ylim(0.30, 1.06)
    twin.set_ylabel(COSINE_LABEL, fontsize=6.2, color=C.VERMILLION)
    twin.tick_params(axis="y", colors=C.VERMILLION)

    # The same-seed repeat of the top level: an uncertainty mark on a sweep that
    # has one seed per level. n = 2 is a pair, and the legend says so.
    r_rank = repeat[1][READ_AT][column]
    r_cos = repeat[1][READ_AT]["rna-rna"]
    mark_rank, = ax.plot([x[-1]], [r_rank], marker="o", markersize=8.5,
                         markerfacecolor="none", markeredgecolor=C.BLUE,
                         markeredgewidth=1.4, linestyle="none", zorder=5)
    twin.plot([x[-1]], [r_cos], marker="s", markersize=8.5, markerfacecolor="none",
              markeredgecolor=C.VERMILLION, markeredgewidth=1.4, linestyle="none", zorder=5)

    ax.set_title(f"({letter})  {title}", pad=11)
    ax.legend([line_rank, line_cos, mark_rank],
              [line_rank.get_label(), line_cos.get_label(),
               "open markers: independent SAME-SEED repeat of the 0.04 arm  —  "
               "n = 2 is a pair, NOT a floor"],
              loc="upper left", handlelength=2.6, fontsize=6.5)
    return fold, r_rank, r_cos


def _track_panel(ax, arms, repeat, column, *, letter, title, ylabel, ticks, legend_loc):
    steps = sorted(arms[LEVELS[0]]["rows"])
    for level, (colour, marker, linestyle) in zip(LEVELS, LEVEL_STYLE):
        ax.plot(steps, [arms[level]["rows"][s][column] for s in steps], color=colour,
                marker=marker, linestyle=linestyle, markersize=4.2,
                markeredgecolor="white", markeredgewidth=0.7,
                label=f"decorrelation = {level}", zorder=3)
    shared = [s for s in steps if s in repeat[1]]
    ax.plot(shared, [repeat[1][s][column] for s in shared], color=C.MUTED, marker="v",
            linestyle=(0, (1, 1.4)), markersize=3.6, linewidth=1.0,
            markerfacecolor="none", zorder=2, label="same-seed repeat of the 0.04 arm")
    ax.set_xlabel("training step", fontsize=6.8, labelpad=4)
    ax.set_ylabel(ylabel, fontsize=6.0)
    if ticks is not None:
        ax.set_yscale("log")
        ax.minorticks_off()
        ax.set_yticks(ticks)
        ax.set_yticklabels([str(t) for t in ticks])
    P.grid(ax)
    ax.set_title(f"({letter})  {title}", pad=9, fontsize=7.6)
    ax.legend(loc=legend_loc, handlelength=2.4, fontsize=6.0)


def main() -> int:
    P.cli(__doc__)
    arms = ablation()
    repeat = repeat_of_the_top_level()
    floor_r3 = floor("R3")
    floor_r1 = floor("R1")

    fig = plt.figure(figsize=(7.9, 11.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.78], hspace=0.60, wspace=0.42,
                          left=0.132, right=0.845, top=0.958, bottom=0.212)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])

    fold_r3, _, _ = _dissociation_panel(
        ax_a, arms, repeat, "R3-rank", "R3_short", floor_r3, [3, 4, 5, 6, 8, 10, 14, 18],
        letter="a",
        title=("Both quantities rise, monotonically, on the SAME three runs\n"
               "      rank statistic R3 — the column these logs' own `final_eff_rank=` line reports"))
    fold_r1, _, _ = _dissociation_panel(
        ax_b, arms, repeat, "CANONICAL", "R1_short", floor_r1, [5, 7, 10, 14, 20, 26],
        letter="b",
        title=("The direction survives the choice of rank statistic\n"
               "      canonical R1 — the statistic §4.1's exported-block floor is measured in"))

    _track_panel(ax_c, arms, repeat, "R3-rank", letter="c",
                 title="the rank track",
                 ylabel="effective rank\n" + P.STAT["R3_short"]
                        + "\nblock: " + P.BLOCK["probe"],
                 ticks=[2, 5, 10, 20, 40, 70], legend_loc="upper right")
    _track_panel(ax_d, arms, repeat, "rna-rna", letter="d",
                 title="the collapse measure",
                 ylabel="RNA-view mutual cosine\nblock: " + P.BLOCK["probe"]
                        + "\n1.0 = total collapse",
                 ticks=None, legend_loc="lower right")
    ax_d.set_ylim(0.18, 1.04)
    ax_d.axhline(1.0, color=C.MUTED, linewidth=0.9, linestyle=(0, (4, 2)), zorder=1)

    fig.text(0.012, 0.176,
             f"ONE SEED PER LEVEL.  The rank change is ×{fold_r3:.3f} under R3 and "
             f"×{fold_r1:.3f} under the canonical statistic — both CLEAR this block's own step-{READ_AT} "
             f"floors of ×{floor_r3:.3f} and ×{floor_r1:.3f},\nthe grey bands in (a) and (b). "
             f"§4.1's ×3.295 is a DIFFERENT BLOCK and a different arm and is not drawn here.\n"
             "THE MONOTONICITY AND THE CO-MEASURED COSINE CARRY THIS RESULT, NOT THE MAGNITUDE OF "
             "THE RANK CHANGE, and no sentence about it may be quoted as though they did — "
             "least of all the pass.",
             ha="left", va="top", fontsize=7.0, color=C.INK, linespacing=1.7,
             fontweight="bold")

    P.caption(fig,
        "F9. Draft section 4.9a. Raising the covariance-decorrelation weight raises effective rank monotonically across three levels while the RNA-view "
        "patient-to-patient mutual cosine - a direct measurement of the condition rank exists to detect - rises monotonically WITH it, co-measured on the "
        "identical runs and printed on the same log lines. At 0.8696 the RNA-view states of different patients are nearly the same vector, and rank moves "
        "the wrong way as that worsens. This is a stronger dissociation than draft section 4.9's instances for two reasons about the SHAPE of the evidence "
        "rather than its size: it is monotone across three levels rather than a single contrast, and the contradicting quantity is co-measured rather than "
        "inferred from a downstream readout in another table. "
        "The floor drawn in (a) and (b) is THIS BLOCK'S OWN, per statistic, at the reading step these runs are quoted at: five same-seed repeats of each of "
        "two arms on the fixed held-out probe at step 400 give 1.4489x under R3 and 1.5702x under canonical R1, and the sweep clears both. Draft section "
        "4.1's 3.295x is canonical R1 on the residualised EXPORTED wsi_biology block of a DIFFERENT ARM at 40 epochs; it never licensed this comparison and "
        "is deliberately not drawn - an earlier version of this figure drew it and printed BOTH INSIDE, which is the block mismatch draft section 4.1a "
        "exists to catch. Clearing is not what the figure claims: at one seed per level a magnitude that clears an n = 5 floor is not evidence of a "
        "dose-response, and the floor itself is n = 5 per arm at one seed on one stack. Draft section 4.1a rows 54-55. "
        "It also carries a correction: \"feature_decorrelation is defective\" was CONDITIONAL ON A QUERY-WRITTEN QUEUE. Without a momentum key encoder the "
        "same term aggravated the collapse (1.59 against 2.17 at step 250, draft section 5.1 instance 3); with one it raises rank. Every claim this project "
        "makes about that term needs \"in the absence of a momentum key encoder\" attached. "
        "Statistic and block are named on every axis, and R3 and the canonical statistic occupy separate panels because binding constraint 1 forbids two "
        "rank statistics on one axis.",
        y=0.006)

    P.save(fig, "F9_decorrelation_dissociation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
