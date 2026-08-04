"""F7 - The dose-response: rank under-reports the loss.  (P2 draft 4.8)

Claim: rank is miscalibrated in magnitude against the information channel by
between about 2x and about 22x depending on implementation choices, over a
monotone seven-level sweep on a representation with ZERO FITTED PARAMETERS - so
retraining noise does not apply, and both quantities are read from the same
representation through the same instrument.

ONE AXIS, NOT TWO. Panel (a) puts both rank curves and the channel on a SINGLE
axis of percentage change from each series' own d = 0 value. That is what the
plan asks for ("or the panel understates the divergence") and it also avoids the
dual-scale reading that a twin axis invites. Both rank curves are drawn: the raw
block is where the published -17.8% came from and the residualised block is the
one the channel is measured on, and the correction this paper carries is
precisely that the two were mixed.

Sources
  v2/research/rebase/nature/DILUTION_LOWER_BOUND.md 2   the seven-level table
  data/ws_rank/RANK_RECOMPUTE.json  (`instance4_dilution`)  all six rank variants,
      raw and residualised, on ~/p1_out/dilution/dilution_foreign_tumour_pca256.npz
  Build scripts: v2/research/dilution/
"""
from __future__ import annotations

import re

import p2fig as P
from p2fig import C, np, plt

DILUTION_MD = "v2/research/rebase/nature/DILUTION_LOWER_BOUND.md"
LEVELS = ["000", "010", "020", "030", "040", "060", "080"]

#: (variant key, block, printed label)
MISCAL = [
    ("R1", "raw", "R1   raw block"),
    ("R2", "raw", "R2   raw block"),
    ("R3", "raw", "R3   raw block"),
    ("R1_uncentred", "raw", "R1 uncentred   raw block"),
    ("R1", "residualised", "R1   RESIDUALISED block"),
    ("R2", "residualised", "R2   residualised block"),
    ("R3", "residualised", "R3   residualised block"),
]


def sweep() -> dict:
    """The seven-level table, read from its own source file."""
    md = P.repo_text(DILUTION_MD)
    rows = [r for r in P.md_rows(md, n_cols=10)
            if re.fullmatch(r"[\d.]+", r[0]) and len(P.nums(r[1])) == 1]
    assert len(rows) == 7, [r[0] for r in rows]
    d = {"achieved": [], "cca_heldout": [], "ratio_raw": [], "ratio_null_corrected": [],
         "attenuation": [], "rank_printed": []}
    for r in rows:
        d["achieved"].append(P.one_num(r[1]))
        d["cca_heldout"].append(P.one_num(r[3]))
        d["ratio_raw"].append(P.one_num(r[4]))
        d["ratio_null_corrected"].append(P.one_num(r[5]))
        d["attenuation"].append(P.one_num(r[7]))
        d["rank_printed"].append(P.one_num(r[8]))
    # NB: the source writes the null as an EN-DASHED RANGE, "0.145-0.147". It must
    # not go through nums(), which normalises dashes to minus signs and would
    # return [0.145, -0.147]; the two endpoints are matched positively instead.
    null = re.findall(r"\d+\.\d+",
                      re.search(r"permutation null median is (\S+) at every level",
                                md).group(1))
    assert len(null) == 2, null
    d["null_lo"], d["null_hi"] = float(null[0]), float(null[1])
    assert 0.0 < d["null_lo"] <= d["null_hi"] < 1.0, null
    return d


def ranks() -> dict:
    rr = P.load_json("ws_rank/RANK_RECOMPUTE.json")["instances"]["instance4_dilution"]
    return {lvl: rr[f"wsi_dilution_{lvl}"] for lvl in LEVELS}


def pct_change(series: list[float]) -> list[float]:
    return [100.0 * (v / series[0] - 1.0) for v in series]


def main() -> int:
    P.cli(__doc__)
    s = sweep()
    rk = ranks()
    d = s["achieved"]

    raw_r1 = [rk[l]["raw"]["R1"] for l in LEVELS]
    res_r1 = [rk[l]["residualised"]["R1"] for l in LEVELS]
    # The printed rank column of the source table must be the raw R1 curve.
    for got, want in zip(raw_r1, s["rank_printed"]):
        assert abs(got - want) < 0.06, (
            f"F7: the recomputed raw R1 dilution curve is {raw_r1!r} but "
            f"{DILUTION_MD} 2 prints {s['rank_printed']!r}. Refusing to plot either.")

    chan_pct = pct_change(s["ratio_null_corrected"])
    curves = [
        ("R1, RAW block", pct_change(raw_r1), C.BLUE, "o", P.LINESTYLES[0]),
        ("R1, RESIDUALISED block", pct_change(res_r1), C.SKY, "^", P.LINESTYLES[1]),
        ("null-corrected information channel", chan_pct, C.VERMILLION, "s", P.LINESTYLES[2]),
    ]

    fig = plt.figure(figsize=(7.5, 9.9))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.0, 0.85], hspace=0.52,
                          left=0.155, right=0.975, top=0.955, bottom=0.185)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])

    # ------------------------------------------------------------------ (a)
    ax_a.set_title("(a)  Both rank curves and the channel, on ONE axis of percentage change\n"
                   "      from each series' own $d$ = 0 value", pad=10)
    ax_a.axhline(0.0, color=C.RULE, linewidth=0.9)
    for name, vals, colour, mk, ls in curves:
        ax_a.plot(d, vals, color=colour, marker=mk, linestyle=ls, markersize=5.0,
                  markeredgecolor="white", markeredgewidth=0.8, label=name, zorder=3)
        ax_a.annotate(f"{vals[-1]:+.1f}%", (d[-1], vals[-1]), xytext=(7, 0),
                      textcoords="offset points", ha="left", va="center",
                      fontsize=7.0, color=colour, fontweight="bold")
    ax_a.set_xlim(-0.03, 0.95)
    ax_a.set_ylim(-72, 8)
    ax_a.set_xticks(d)
    ax_a.set_xticklabels([f"{v:.3f}" for v in d], fontsize=6.6)
    ax_a.set_xlabel("achieved dilution $d$  (test median fraction of foreign-donor patches)",
                    fontsize=6.8, labelpad=6)
    ax_a.set_ylabel("percentage change from that series' own $d$ = 0 value\n"
                    "rank: " + P.STAT["R1_short"] + "\n"
                    "channel: null-corrected held-out top-CCA ratio", fontsize=6.4)
    P.grid(ax_a)
    ax_a.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), handlelength=2.6)
    ax_a.text(0.985, 0.62,
              "Single seed (42), single draw of donor assignments,\n"
              "NO error bar on any level-to-level difference.\n"
              "The ratio curve is not monotone: it dips to "
              f"{s['ratio_null_corrected'][1]:.3f} at $d$ = {d[1]:.3f},\n"
              "a single-seed wobble on a level where the channel moved by 0.001.",
              transform=ax_a.transAxes, ha="right", va="top", fontsize=6.1,
              color=C.MUTED, linespacing=1.55)

    # ------------------------------------------------------------------ (b)
    ax_b.set_title("(b)  The miscalibration, as a factor:  how many times more the channel moved\n"
                   "      than the rank statistic did, over the same seven levels", pad=10)
    chan_drop = abs(chan_pct[-1])
    bars = []
    for key, blk, label in MISCAL:
        series = [rk[l][blk][key] for l in LEVELS]
        drop = abs(pct_change(series)[-1])
        bars.append((label, chan_drop / drop, drop, blk))
    bars.sort(key=lambda b: b[1])
    for k, (label, factor, drop, blk) in enumerate(bars):
        colour = C.BLUE if blk == "raw" else C.SKY
        ax_b.barh(k, factor, height=0.6, color=colour, edgecolor="white", linewidth=1.5,
                  hatch="" if blk == "raw" else "///", zorder=2)
        ax_b.annotate(f"{factor:.2f}$\\times$    (rank {-drop:.2f}%)", (factor, k),
                      xytext=(6, 0), textcoords="offset points", ha="left", va="center",
                      fontsize=6.8, color=C.INK, fontweight="bold")
    ax_b.set_yticks(range(len(bars)))
    ax_b.set_yticklabels([b[0] for b in bars], fontsize=6.8)
    ax_b.tick_params(axis="y", length=0)
    ax_b.set_ylim(-0.65, len(bars) - 0.35)
    ax_b.set_xlim(0, 27)
    ax_b.set_xlabel(f"channel change ({-chan_drop:.1f}%) divided by the rank statistic's own change\n"
                    "over the same seven dilution levels", fontsize=6.6, labelpad=6)
    P.grid(ax_b, axis="x")
    ax_b.text(0.99, 0.06,
              "One message: the DIRECTION survives every implementation\n"
              "choice and the MAGNITUDE does not.",
              transform=ax_b.transAxes, ha="right", va="bottom", fontsize=6.4,
              color=C.INK, linespacing=1.55)

    # ------------------------------------------------------------------ (c)
    ax_c.set_title("(c)  The instrument's own controls over the same levels", pad=10)
    ax_c.axhline(1.0, color=C.MUTED, linewidth=1.0, linestyle=(0, (4, 2)), zorder=1)
    ax_c.plot(d, s["attenuation"], color=C.GREEN, marker="D", linestyle=P.LINESTYLES[0],
              markersize=5.0, markeredgecolor="white", markeredgewidth=0.8, zorder=3,
              label="confound-adjustment attenuation (reference = 1)")
    ax_c.plot(d, s["cca_heldout"], color=C.VERMILLION, marker="s",
              linestyle=P.LINESTYLES[2], markersize=5.0, markeredgecolor="white",
              markeredgewidth=0.8, zorder=3, label="held-out top-CCA, absolute")
    ax_c.fill_between([-0.03, 0.95], s["null_lo"], s["null_hi"], color=C.ENVELOPE,
                      alpha=0.55, linewidth=0, zorder=0)
    ax_c.annotate(f"permutation null median, {s['null_lo']:.3f}-{s['null_hi']:.3f}\n"
                  "at EVERY level (within-cancer, this sweep's own null)",
                  (0.02, s["null_hi"]), xytext=(0, 6), textcoords="offset points",
                  ha="left", va="bottom", fontsize=6.2, color=C.MUTED)
    ax_c.set_xlim(-0.03, 0.95)
    ax_c.set_ylim(0, 1.32)
    ax_c.set_xticks(d)
    ax_c.set_xticklabels([f"{v:.3f}" for v in d], fontsize=6.6)
    ax_c.set_xlabel("achieved dilution $d$", fontsize=6.8, labelpad=6)
    ax_c.set_ylabel("attenuation (unitless) and\nabsolute held-out top-CCA", fontsize=6.5)
    P.grid(ax_c)
    ax_c.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), handlelength=2.6)
    ax_c.text(0.0, -0.245,
              "The readout did not degrade: attenuation stays near 1 at every level, and every level "
              "clears this sweep's OWN\nshuffled-pairing null. That null (0.145-0.147) belongs to this "
              "experiment and to no other; it must never appear on\nan axis with D2's 0.140, which "
              "differs in n, in component count and in the permutation procedure itself.",
              transform=ax_c.transAxes, ha="left", va="top", fontsize=6.2, color=C.INK,
              linespacing=1.6)

    P.caption(fig,
        "F7. THIS INSTANCE DOES NOT CONTRADICT RANKME AS STATED. High rank with degraded information is exactly the necessary-not-sufficient case "
        "RankMe reserves for itself, and LiDAR's own text says so. What is added here is the MAGNITUDE of the miscalibration, on a representation with "
        "zero fitted parameters - so section 3.5's retraining noise does not apply, and rank and channel are read from the same representation through "
        "the same instrument. Single seed (42), single draw of donor assignments, no error bar on any level-to-level difference; the detection floor is "
        "censored at >= 0.40 from d = 0.09 onward and the transmission floor at 0.05 from below; the curve is a property of unweighted mean pooling, not "
        "of the modality. The source file's own section 4 WITHDRAWS the phrase \"lower bound\": the measured quantity is the cost of preparation-matched, "
        "information-free contamination. Do not describe (a)'s ratio curve as monotone.",
        y=0.006)

    P.save(fig, "F7_dilution_dose_response")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
