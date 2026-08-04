"""F1 - The retraining envelope, measured, and the seven arm differences inside it.  (P2 draft 4.1)

Claim: every between-arm rank difference this project has ever measured is
smaller than the spread of the same statistic when ONE configuration is retrained
with the SAME seed.

THE ENVELOPE IS NO LONGER n = 1. Until 2026-08-04 it rested on a single retraining
pair (2.69x) and panel (d) was a hatched placeholder for the controlled repeat that
would replace it. That repeat has reported: five identical `programme_only` runs at
seed 42, differing only in GPU non-determinism, spreading 3.295x in rank and 1.055x
in the channel those same five runs carry.

Two things the drawing must not smooth over.

**The distribution is bimodal, not a spread.** Four repeats agree to within 2% and
one lands at a third of them. Panel (a) plots the five values individually and
never a mean or a band, because a band would invite the reader to imagine a
distribution the data does not have.

**The envelope is a FLOOR.** `programme_only` is the stable arm and same-seed
repeats exclude seed variation entirely, so 3.295x understates what a practitioner
faces. Every panel that draws it labels it a floor.

The floor is drawn PER BLOCK. Six of the seven arm comparisons are on the
residualised block (floor 3.295x) and Phase 1b is on the raw block (floor 3.111x);
comparing a residualised ratio against a raw floor would flip D1-B seed 43, which
is exactly the raw/residualised confusion draft 4.5(b) is about.

Sources
  (a) data/extracted/F1_RETRAINING_REPEAT.json          the five repeats, parsed from
      data/e0_run/d1_envelope_readout.log               the readout log, byte for byte
  (b) data/ws_rank/RANK_RECOMPUTE.json                  D2 + D1-B, R1, residualised
      data/ws_rank/RANK_RECOMPUTE_P1B.json              Phase 1b, R1, RAW
  (c) v2/research/rebase/nature/D2_RESULT.md 2          stratified 40-target readout
  (d) data/e0_run/d2_v3/RECOVERED_SEED42_READOUT.json   the superseded n = 1 estimate
      data/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json
"""
from __future__ import annotations

import p2fig as P
from p2fig import C, np, plt


def repeats() -> dict:
    """The five same-seed retrains, and the floor they establish.

    Nothing is recomputed. The per-repeat values were written by
    `v2/research/rebase/d1_envelope_readout.py`, which imports every statistic from
    `calibra`; this reads them out of that run's log. The spreads ARE recomputed
    from the per-repeat values and asserted against the spreads the log itself
    printed, so a partial re-extraction that mixed two runs would fail here rather
    than be drawn.
    """
    rec = P.load_json("extracted/F1_RETRAINING_REPEAT.json")
    if not rec.get("complete"):
        raise SystemExit(
            "F1: the controlled retraining repeat has NOT reported, but this script has "
            "been rewritten around it. Restore the pending placeholder rather than drawing "
            "a partial envelope; see NOTEBOOK_ENTRIES/"
            "PREDECLARED_retraining_envelope_20260804T0330Z.md.")
    reps = [(name, rec["reps"][name]) for name in sorted(rec["reps"])]
    floors = {}
    for key, tag in (("rank_residualised", "resid"), ("rank_raw", "raw"), ("channel", "channel")):
        values = [r[key] for _n, r in reps]
        got = max(values) / min(values)
        printed = rec["printed_spread"][key]
        assert abs(got - printed["spread"]) < 1e-3, (
            f"F1: the five repeats give a {key} spread of {got:.4f}x, but "
            f"{rec['readout_log']} printed {printed['spread']}x. Refusing to draw either.")
        assert abs(max(values) - printed["max"]) < 5e-4, (key, max(values), printed)
        assert abs(min(values) - printed["min"]) < 5e-4, (key, min(values), printed)
        floors[tag] = got
    return {"reps": reps, "floor": floors, "meta": rec}


def superseded_n1() -> dict:
    """The n = 1 estimate the measured floor replaces: re-export against retrain."""
    rec = P.load_json("e0_run/d2_v3/RECOVERED_SEED42_READOUT.json")["H42_recovered"]
    ret = P.load_json("e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json")["H42"]
    lo = rec["effective_rank_residualised"]
    hi = ret["effective_rank_residualised"]
    return {"reexport_rank": lo, "retrained_rank": hi, "fold": hi / lo,
            "reexport_channel": rec["points"]["unrestricted"]["top_cca"],
            "retrained_channel": ret["points"]["unrestricted"]["top_cca"]}


def arm_ratios() -> list[tuple[str, float, str]]:
    """The seven between-arm rank ratios, as (label, fold, statistic/block tag)."""
    rr = P.load_json("ws_rank/RANK_RECOMPUTE.json")
    d2 = rr["instances"]["instance6_D2"]
    d1 = rr["instances"]["instance5_D1B"]
    out = []
    for seed in (42, 43, 44):
        h, i = d2[f"H{seed}"]["residualised"]["R1"], d2[f"I{seed}"]["residualised"]["R1"]
        hi = "H" if h > i else "I"
        out.append((f"D2 s{seed} ({hi} hi)", max(h, i) / min(h, i), "resid"))
    for seed in (42, 43, 44):
        p = d1[f"d1_p_seed{seed}"]["residualised"]["R1"]
        f = d1[f"d1_f_seed{seed}"]["residualised"]["R1"]
        hi = "prog_only" if p > f else "prog_free"
        out.append((f"D1-B s{seed} ({hi} hi)", max(p, f) / min(p, f), "resid"))
    # Phase 1b is the one RAW-block point and its arms were never verified
    # matched (draft 3.4). It is tagged so panel (b) compares it against the RAW
    # floor and marks it.
    p1b = P.load_json("ws_rank/RANK_RECOMPUTE_P1B.json")["v21_release_20260720_retry3_resume_safe"]
    full = p1b["diagnostic_full_seed42"]["wsi_biology"]["raw"]["R1"]
    prog = p1b["diagnostic_programme_only_seed42"]["wsi_biology"]["raw"]["R1"]
    out.append(("Phase 1b s42 (full hi)", max(full, prog) / min(full, prog), "raw"))
    return sorted(out, key=lambda r: r[1])


def d2_paired_channel() -> list[dict]:
    """Paired within-run channel differences, PBS - Hallmark, with both bootstraps.

    Scoped to the STRATIFIED section: D2_RESULT.md carries three tables of this
    shape and the figure must read the 40 targets neither arm trained on, not the
    90 unrestricted ones (50 of which are Hallmark's own supervision).
    """
    md = P.repo_text("v2/research/rebase/nature/D2_RESULT.md")
    sec = P.md_section(md, "THE STRATIFIED READOUT")
    rows = [r for r in P.md_rows(sec, n_cols=8) if r[0] in ("42", "43", "44")]
    assert len(rows) == 3, rows
    out = []
    for r in rows:
        d = {
            "seed": r[0],
            "hallmark": P.one_num(r[1]),
            "pbs": P.one_num(r[2]),
            "delta": P.one_num(r[3]),
            "patient": P.ci(r[4]),
            "cancer": P.ci(r[6]),
        }
        # The table quotes PBS - Hallmark; guard the sign so a source edit that
        # reversed the direction would fail here rather than flip the figure.
        assert d["delta"] < 0 < d["hallmark"] - d["pbs"], d
        out.append(d)
    return out


def main() -> int:
    P.cli(__doc__)
    env = repeats()
    n1 = superseded_n1()
    ratios = arm_ratios()
    pairs = d2_paired_channel()
    floor_res, floor_raw = env["floor"]["resid"], env["floor"]["raw"]
    ch_fold = env["floor"]["channel"]

    fig = plt.figure(figsize=(7.6, 10.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=1.02, wspace=0.80,
                          left=0.135, right=0.985, top=0.958, bottom=0.300)
    gs_a = gs[0, 0].subgridspec(2, 1, height_ratios=[1.45, 1.0], hspace=0.30)
    ax_a = fig.add_subplot(gs_a[0])
    ax_a2 = fig.add_subplot(gs_a[1], sharex=ax_a)
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # ---------------------------------------------------------------- (a)
    # Five identical runs, rank above and channel below on a shared x-axis. The
    # two are stacked rather than sharing an axis: they are different quantities
    # in different units, and the comparison the paper makes is between their
    # SPREADS, which is what stacking shows without implying a common scale.
    ax_a.set_title("(a)  The envelope, measured: 5 identical retrains", pad=9)
    x = np.arange(1, len(env["reps"]) + 1)
    rank = np.array([r["rank_residualised"] for _n, r in env["reps"]])
    chan = np.array([r["channel"] for _n, r in env["reps"]])

    ax_a.set_yscale("log")
    ax_a.axhspan(rank.min(), rank.max(), color=C.ENVELOPE, alpha=0.30, zorder=0, linewidth=0)
    ax_a.plot(x, rank, marker="o", linestyle="none", color=C.BLUE, markersize=7.5,
              markeredgecolor="white", markeredgewidth=1.0, zorder=3)
    for xi, v in zip(x, rank):
        ax_a.annotate(f"{v:.3f}", (xi, v), xytext=(0, 9), textcoords="offset points",
                      ha="center", fontsize=6.3, color=C.BLUE, fontweight="bold")
    ax_a.set_ylim(rank.min() / 1.6, rank.max() * 1.9)
    ax_a.set_xlim(0.45, len(x) + 0.55)
    ax_a.set_ylabel("R1 effective rank\ncentred, residualised block",
                    fontsize=6.2, color=C.BLUE)
    ax_a.tick_params(axis="y", colors=C.BLUE, labelsize=6.2)
    ax_a.tick_params(axis="x", labelbottom=False)
    P.grid(ax_a, axis="y")
    ax_a.text(0.985, 0.04, f"spread $\\times${floor_res:.3f}", transform=ax_a.transAxes,
              ha="right", va="bottom", fontsize=7.0, color=C.BLUE, fontweight="bold")

    ax_a2.axhspan(chan.min(), chan.max(), color=C.ENVELOPE, alpha=0.30, zorder=0, linewidth=0)
    ax_a2.plot(x, chan, marker="s", linestyle="none", color=C.VERMILLION, markersize=6.0,
               markeredgecolor="white", markeredgewidth=0.9, zorder=3)
    for xi, v in zip(x, chan):
        ax_a2.annotate(f"{v:.4f}", (xi, v), xytext=(0, 8), textcoords="offset points",
                       ha="center", fontsize=6.1, color=C.VERMILLION, fontweight="bold")
    ax_a2.set_ylim(chan.min() - 0.012, chan.max() + 0.022)
    ax_a2.set_ylabel("held-out top-CCA\n40 untrained targets", fontsize=6.2, color=C.VERMILLION)
    ax_a2.tick_params(axis="y", colors=C.VERMILLION, labelsize=6.2)
    ax_a2.set_xticks(x)
    ax_a2.set_xticklabels([n.replace("rep", "") for n, _r in env["reps"]], fontsize=6.8)
    ax_a2.set_xlabel("repeat  (same seed 42, same configuration;\n"
                     "GPU non-determinism is the only source of variation)",
                     fontsize=6.3, labelpad=5)
    P.grid(ax_a2, axis="y")
    ax_a2.text(0.985, 0.03, f"spread $\\times${ch_fold:.3f}", transform=ax_a2.transAxes,
               ha="right", va="bottom", fontsize=7.0, color=C.VERMILLION, fontweight="bold")

    odd = int(np.argmin(rank))
    ax_a.annotate(f"rep {x[odd]}", (x[odd], rank[odd]), xytext=(16, 0),
                  textcoords="offset points", fontsize=6.4, color=C.INK, va="center",
                  fontweight="bold",
                  arrowprops=dict(arrowstyle="-", color=C.MUTED, linewidth=0.7))
    ax_a2.text(0.0, -0.72,
               f"rank $\\times${floor_res:.3f}  against  channel $\\times${ch_fold:.3f}, "
               "on identical inputs.\n"
               f"BIMODAL, not a spread: four repeats agree to 2%, rep {x[odd]} lands at a third\n"
               f"(rank $\\div${rank.max() / rank.min():.2f}, channel $-$"
               f"{(1 - chan.min() / chan.max()) * 100:.0f}%). Reproducible ~80% of the time and\n"
               "catastrophically not ~20% of the time. A FLOOR: programme_only is the\n"
               "stable arm and the seed is held fixed, so both understate the truth.",
               transform=ax_a2.transAxes, ha="left", va="top", fontsize=6.2, color=C.INK,
               linespacing=1.5)

    # ---------------------------------------------------------------- (b)
    ax_b.set_title("(b)  The seven between-arm differences, against it", pad=10)
    ax_b.axvspan(1.0, floor_res, color=C.ENVELOPE, alpha=0.32, zorder=0, linewidth=0)
    ax_b.axvline(floor_res, color=C.MUTED, linewidth=1.0, linestyle=(0, (4, 2)), zorder=1)
    ax_b.axvline(floor_raw, color=C.MUTED, linewidth=0.8, linestyle=(0, (1, 1.6)), zorder=1)
    ax_b.axvline(1.0, color=C.RULE, linewidth=0.8, zorder=1)
    inside = 0
    for k, (lab, val, tag) in enumerate(ratios):
        # Each ratio is judged against the floor measured on ITS OWN block.
        within = val <= (floor_raw if tag == "raw" else floor_res)
        inside += within
        ax_b.plot([val], [k], marker="o" if tag == "resid" else "^",
                  color=C.BLUE if within else C.VERMILLION, markersize=6.5,
                  markeredgecolor="white", markeredgewidth=0.9, zorder=3)
        ax_b.annotate(f"{val:.3f}$\\times$", (val, k), xytext=(0, 7), textcoords="offset points",
                      ha="center", fontsize=6.4,
                      color=C.BLUE if within else C.VERMILLION, fontweight="bold")
    ax_b.set_yticks(range(len(ratios)))
    ax_b.set_yticklabels([r[0] for r in ratios], fontsize=6.2)
    ax_b.tick_params(axis="y", length=0)
    ax_b.set_xscale("log")
    ax_b.set_xlim(0.95, 4.3)
    ax_b.set_xticks([1.0, 1.5, 2.0, round(floor_res, 2), 4.0])
    ax_b.set_xticklabels(["1.0", "1.5", "2.0", f"{floor_res:.2f}", "4.0"])
    ax_b.minorticks_off()
    ax_b.set_ylim(-0.6, len(ratios) - 0.3)
    ax_b.set_xlabel("between-arm rank ratio, higher / lower  (log scale)\n"
                    "R1 (canonical effective rank, centred);  circle = residualised block,\n"
                    "triangle = Phase 1b, RAW block, arms never verified matched", fontsize=6.4,
                    labelpad=6)
    P.grid(ax_b, axis="x")
    ax_b.text(0.0, -0.52,
              f"shaded: the measured retraining FLOOR, 1.00-{floor_res:.3f}$\\times$ "
              f"(residualised).\n"
              f"Dotted line: the RAW-block floor, {floor_raw:.3f}$\\times$ - the Phase 1b point\n"
              f"is judged against that one, not against the residualised floor.\n"
              "Five repeats, one arm, one seed: a floor, not a bound.",
              transform=ax_b.transAxes, fontsize=6.0, color=C.MUTED, va="top", ha="left",
              linespacing=1.5)
    ax_b.text(0.03, 0.975, f"{inside} of {len(ratios)} inside the floor",
              transform=ax_b.transAxes, ha="left", va="top",
              fontsize=7.2, color=C.BLUE, fontweight="bold")

    # ---------------------------------------------------------------- (c)
    ax_c.set_title("(c)  The asymmetry, beside it", pad=10)
    ax_c.axvline(0.0, color=C.INK, linewidth=0.9, zorder=1)
    for k, p in enumerate(pairs):
        y = len(pairs) - 1 - k
        for name, key, colour, mk, off in [
                ("patient bootstrap", "patient", C.BLUE, "o", +0.17),
                ("cancer-cluster bootstrap", "cancer", C.VERMILLION, "s", -0.17)]:
            clo, chi = p[key]
            ax_c.plot([clo, chi], [y + off, y + off], color=colour, linewidth=2.0,
                      solid_capstyle="butt", zorder=2, label=name if k == 0 else None)
            ax_c.plot([p["delta"]], [y + off], marker=mk, color=colour, markersize=5.0,
                      markeredgecolor="white", markeredgewidth=0.9, zorder=3)
        ax_c.annotate(f"{p['delta']:+.4f}", (p["delta"], y + 0.40), ha="center", va="bottom",
                      fontsize=6.5, color=C.INK, fontweight="bold")
    ax_c.set_yticks(range(len(pairs)))
    ax_c.set_yticklabels([f"D2 s{p['seed']}" for p in reversed(pairs)], fontsize=6.8)
    ax_c.tick_params(axis="y", length=0)
    ax_c.set_ylim(-0.75, len(pairs) - 0.15)
    ax_c.set_xlim(-0.205, 0.215)
    ax_c.set_xlabel("paired channel difference, PBS $-$ Hallmark\n"
                    "held-out top-CCA, 40 targets neither arm trained on\n"
                    "confound-residualised block;  CI$_{95}$, 2,000 resamples", fontsize=6.4,
                    labelpad=6)
    P.grid(ax_c, axis="x")
    ax_c.legend(loc="upper right", bbox_to_anchor=(1.0, 1.02), ncol=1, handlelength=1.8)
    same_sign = sum(p["delta"] < 0 for p in pairs)
    excl = sum(p["patient"][1] < 0 and p["cancer"][1] < 0 for p in pairs)
    ax_c.text(0.985, 0.40,
              f"same sign {same_sign}/3\nboth CIs exclude zero {excl}/3\n\n"
              "the channel is a PAIRED within-run\ndifference; rank is not, and\n"
              "there is no paired form of\n\"this run has higher rank\"",
              transform=ax_c.transAxes, ha="right", va="center", fontsize=6.4,
              color=C.INK, linespacing=1.55)

    # ---------------------------------------------------------------- (d)
    # What the measured floor replaced, and what replacing it cost and bought.
    ax_d.set_title("(d)  What the floor replaced, and what it costs us", pad=10)
    ax_d.set_xscale("log")
    ax_d.axvspan(1.0, floor_res, color=C.ENVELOPE, alpha=0.32, zorder=0, linewidth=0)
    ax_d.axvline(1.0, color=C.RULE, linewidth=0.8, zorder=1)
    rows_d = [
        (f"MEASURED floor\n5 repeats, seed 42", floor_res, C.BLUE, "o"),
        (f"superseded $n$ = 1 estimate\n{n1['reexport_rank']:.3f} $\\to$ "
         f"{n1['retrained_rank']:.3f}", n1["fold"], C.MUTED, "D"),
    ]
    d1_ratios = [(lab, val) for lab, val, _t in ratios if lab.startswith("D1-B")]
    for k, (lab, val, colour, mk) in enumerate(rows_d):
        y = len(rows_d) + len(d1_ratios) - 1 - k
        ax_d.plot([val], [y], marker=mk, color=colour, markersize=8,
                  markeredgecolor="white", markeredgewidth=1.0, zorder=3)
        ax_d.annotate(f"{val:.3f}$\\times$", (val, y), xytext=(0, 9), textcoords="offset points",
                      ha="center", fontsize=6.5, color=colour, fontweight="bold")
    d1_inside = 0
    for k, (lab, val) in enumerate(sorted(d1_ratios, key=lambda r: -r[1])):
        y = len(d1_ratios) - 1 - k
        within = val <= floor_res
        d1_inside += within
        ax_d.plot([val], [y], marker="s", color=C.VERMILLION if within else C.GREEN,
                  markersize=6.0, markeredgecolor="white", markeredgewidth=0.9, zorder=3)
        ax_d.annotate(f"{val:.3f}$\\times$", (val, y), xytext=(0, 8), textcoords="offset points",
                      ha="center", fontsize=6.2,
                      color=C.VERMILLION if within else C.GREEN, fontweight="bold")
    labels_d = [r[0] for r in rows_d] + [lab.split(" (")[0] + "\n(D1's necessity test)"
                                         for lab, _v in sorted(d1_ratios, key=lambda r: -r[1])]
    ax_d.set_yticks(range(len(labels_d)))
    ax_d.set_yticklabels(list(reversed(labels_d)), fontsize=6.0)
    ax_d.tick_params(axis="y", length=0)
    ax_d.set_ylim(-0.55, len(labels_d) - 0.25)
    ax_d.set_xlim(0.95, 4.3)
    ax_d.set_xticks([1.0, 1.5, 2.0, round(floor_res, 2), 4.0])
    ax_d.set_xticklabels(["1.0", "1.5", "2.0", f"{floor_res:.2f}", "4.0"])
    ax_d.minorticks_off()
    ax_d.set_xlabel("rank ratio, higher / lower  (log scale)\n"
                    "R1, residualised block", fontsize=6.4, labelpad=6)
    P.grid(ax_d, axis="x")
    # Stated, not asserted: if a future recomputation moved a D1 ratio outside the
    # floor the figure must SAY so rather than fail, because that outcome is one
    # of the four the envelope was predeclared to distinguish.
    ax_d.text(0.0, -0.40,
              f"{d1_inside} of {len(d1_ratios)} D1 rank ratios are inside the floor, so D1 is "
              "UNINFORMATIVE about\n"
              "rank in either direction. The necessity result is NOT refuted: programme_only\n"
              "still wins the channel 3/3 with patient CIs excluding zero (F6), and it stays at\n"
              "full prominence flagged as unresolvable rather than deleted. The count moved\n"
              "6/7 to 7/7 IN OUR FAVOUR - and rests on one measurement, five repeats, one\n"
              "arm, one seed. Reported with the scepticism a result against us would get.",
              transform=ax_d.transAxes, fontsize=6.0, color=C.MUTED, va="top", ha="left",
              linespacing=1.5)

    P.caption(fig,
        "F1. The envelope of panel (a) is FIVE identical retrains of ONE configuration at ONE seed on ONE stack "
        "(D1 `programme_only`, seed 42, GPU non-determinism the only source of variation), read out by "
        "v2/research/rebase/d1_envelope_readout.py, which imports every statistic from calibra. It is a FLOOR twice over: "
        "programme_only is this project's stable arm, and same-seed repeats exclude seed variation entirely. It cannot separate "
        "rank-specific variance from stack non-determinism, architecture or schedule. Statistic R1 (Roy & Vetterli order-1, centred), "
        "block confound-residualised for every D2 and D1-B point and RAW for the Phase 1b point, which is judged against the RAW-block "
        f"floor of {floor_raw:.3f}x rather than the residualised {floor_res:.3f}x - and the Phase 1b arms were never verified matched (draft 3.4). "
        "Panel (b)'s x-axis is a ratio of two R1 values within one comparison; no two rank statistics share an axis anywhere in this figure. "
        "The claim does not rest on this figure alone: F2 reaches the same conclusion from 8 within-arm degrees of freedom without using the "
        "floor at all. Panels (b) and (c) sit side by side because the channel is quoted as a paired within-run difference and rank is not - "
        "there is no paired form of \"this run has higher rank\".",
        y=0.006)

    P.save(fig, "F1_envelope")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
