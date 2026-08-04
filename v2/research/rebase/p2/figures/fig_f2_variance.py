"""F2 - Where rank's variance lives: arm against training seed.  (P2 draft 4.2)

THE PAPER'S MOST IMPORTANT DISPLAY ITEM. It is the contribution that does not
depend on a sign count, and it reaches F1's conclusion from 8 within-arm degrees
of freedom without using the 2.69x number.

Claim: two-thirds of the variation in effective rank across twelve matched
artifacts is training-seed nuisance; two percent of the variation in the
information channel is. The arm effect on rank is not significant; the arm
effect on information is overwhelming.

Design: 4 arms (D2-H, D2-I, D1-P, D1-F) x 3 seeds (42/43/44), one frozen
artifact each. The seed changes nothing about objective, architecture, data,
split or schedule - that is what makes it the nuisance factor.

The decomposition is RECOMPUTED here from the per-artifact metrics and then
asserted against the printed decomposition in the verified-workspace run log, so
a divergence between the two stops the figure instead of being drawn.

Sources
  data/ws_p2/out/P2_METRICS_D2.json      12 artifacts, all metrics
  data/ws_p2/out/P2_METRICS_D1.json
  data/ws_p2/out/p2_run.log              the printed decomposition, cross-check
"""
from __future__ import annotations

import math
import re

import p2fig as P
from p2fig import C, np, plt

ARMS = [
    ("D2-H", "Hallmark", "H", "d2"),
    ("D2-I", "PBS", "I", "d2"),
    ("D1-P", "programme_only", "P", "d1"),
    ("D1-F", "programme_free", "F", "d1"),
]
SEEDS = (42, 43, 44)


def artifacts() -> dict[str, list[dict]]:
    d2 = P.load_json("ws_p2/out/P2_METRICS_D2.json")
    d1 = P.load_json("ws_p2/out/P2_METRICS_D1.json")
    src = {"d2": d2, "d1": d1}
    return {tag: [src[which][f"{key}{s}"] for s in SEEDS] for tag, _, key, which in ARMS}


def decompose(values: dict[str, list[float]], *, log: bool) -> dict:
    """One-way ANOVA with the ARM as the only factor; the seed is the residual.

    Rank-type metrics are decomposed on the LOG scale because they are
    multiplicative and span 6.4-34.1; the CCA on the raw scale. F is
    F(3, 8) = (SS_arm/3) / (SS_seed/8).
    """
    v = {a: [math.log(x) if log else x for x in xs] for a, xs in values.items()}
    flat = [x for xs in v.values() for x in xs]
    grand = sum(flat) / len(flat)
    ss_arm = sum(len(xs) * (sum(xs) / len(xs) - grand) ** 2 for xs in v.values())
    ss_seed = sum((x - sum(xs) / len(xs)) ** 2 for xs in v.values() for x in xs)
    return {
        "ss_arm": ss_arm, "ss_seed": ss_seed,
        "arm_share": 100.0 * ss_arm / (ss_arm + ss_seed),
        "F": (ss_arm / 3.0) / (ss_seed / 8.0),
        "scale": "log" if log else "raw",
    }


def printed_decomposition() -> list[tuple[str, str, float, float, float, float]]:
    """The same table as printed by p2_necessity_and_variance.py, for cross-check."""
    log = P.load_text("ws_p2/out/p2_run.log")
    rows = []
    for line in log.splitlines():
        m = re.match(r"^(.*?)\s{2,}(log|raw)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)%\s+([\d.]+)\s*$",
                     line.strip())
        if m:
            rows.append((m.group(1).strip(), m.group(2), float(m.group(3)), float(m.group(4)),
                         float(m.group(5)), float(m.group(6))))
    assert len(rows) == 3, f"expected 3 decomposition rows in p2_run.log, got {len(rows)}"
    return rows


def main() -> int:
    args = P.cli(__doc__)
    art = artifacts()

    quantities = [
        ("canonical effective rank",
         P.STAT["R1_short"] + "\nresidualised block, log scale",
         lambda a: a["metrics"]["effective_rank_residualised"], True),
        ("RankMe as published",
         P.STAT["RANKME"] + "\nRAW block, log scale",
         lambda a: a["metrics"]["rankme_raw"], True),
        ("held-out top-CCA  (GROUND TRUTH)",
         "top canonical correlation, 40 untrained targets\nresidualised block, raw scale",
         lambda a: a["points"]["untrained40"]["top_cca"], False),
    ]

    dec = []
    for name, sub, get, log in quantities:
        dec.append((name, sub, decompose({t: [get(x) for x in art[t]] for t in art}, log=log)))

    # Cross-check against the printed table. A disagreement is a stop, not a warning.
    for (name, _sub, d), row in zip(dec, printed_decomposition()):
        _label, scale, ss_a, ss_s, share, f = row
        assert scale == d["scale"], (name, scale, d["scale"])
        for got, want, tol, what in ((d["ss_arm"], ss_a, 5e-5, "SS_arm"),
                                     (d["ss_seed"], ss_s, 5e-5, "SS_seed"),
                                     (d["arm_share"], share, 5e-2, "arm share"),
                                     (d["F"], f, 5e-3, "F")):
            assert abs(got - want) <= tol, (
                f"F2: recomputed {what} for {name} is {got!r} but "
                f"data/ws_p2/out/p2_run.log prints {want!r}. Refusing to plot either.")

    fig = plt.figure(figsize=(7.5, 9.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.80, 1.15, 0.72], hspace=0.62, wspace=0.16,
                          left=0.175, right=0.985, top=0.955, bottom=0.135)
    ax_a = fig.add_subplot(gs[0, :])
    ax_br = fig.add_subplot(gs[1, 0])
    ax_bc = fig.add_subplot(gs[1, 1])
    ax_c = fig.add_subplot(gs[2, :])

    # ------------------------------------------------------------------ (a)
    ax_a.set_title("(a)  The decomposition: how much of each quantity is the ARM,\n"
                   "      and how much is the training SEED", pad=8)
    ys = [2.8, 1.4, 0.0]
    for y, (name, sub, d) in zip(ys, dec):
        arm, seed = d["arm_share"], 100.0 - d["arm_share"]
        # 2px surface gap between the two segments.
        ax_a.barh(y, arm, height=0.50, color=C.BLUE, edgecolor="white", linewidth=1.6, zorder=2)
        ax_a.barh(y, seed, left=arm, height=0.50, color=C.ORANGE, hatch="///",
                  edgecolor="white", linewidth=1.6, zorder=2)
        ax_a.text(arm / 2, y, f"arm\n{arm:.1f}%", ha="center", va="center",
                  fontsize=7.4, color="white", fontweight="bold")
        if seed >= 12:
            ax_a.text(arm + seed / 2, y, f"training seed\n{seed:.1f}%", ha="center", va="center",
                      fontsize=7.4, color=C.INK, fontweight="bold")
        else:
            ax_a.annotate(f"training seed {seed:.1f}%", (100, y), xytext=(7, 0),
                          textcoords="offset points", ha="left", va="center",
                          fontsize=7.0, color=C.INK, fontweight="bold")
        sig = "n.s." if d["F"] < 4.07 else "p < 0.05"   # F(3,8) 5% critical value = 4.066
        ax_a.text(-1.5, y + 0.74, name, ha="right", va="top", fontsize=7.2, color=C.INK,
                  fontweight="bold")
        ax_a.text(-1.5, y + 0.52, sub, ha="right", va="top", fontsize=6.0, color=C.MUTED,
                  linespacing=1.4)
        ax_a.text(-1.5, y - 0.40, f"$F(3,8)$ = {d['F']:.2f}   ({sig})", ha="right", va="bottom",
                  fontsize=7.2, color=C.BLUE if d["F"] >= 4.07 else C.INK, fontweight="bold")
    ax_a.set_xlim(0, 100)
    ax_a.set_ylim(-0.85, 3.75)
    ax_a.set_yticks([])
    ax_a.spines["left"].set_visible(False)
    ax_a.set_xticks([0, 25, 50, 75, 100])
    ax_a.set_xticklabels(["0", "25", "50", "75", "100%"])
    ax_a.set_xlabel("share of the total sum of squares over the 12 matched artifacts\n"
                    "rank-type rows decomposed on the LOG scale (they are multiplicative, "
                    "spanning 6.4-34.1); the CCA on the RAW scale", fontsize=6.5, labelpad=6)
    P.grid(ax_a, axis="x")
    ax_a.text(0.0, 1.0, "", transform=ax_a.transAxes)

    # ------------------------------------------------------------------ (b)
    for ax, title, get, logscale, xlab in (
            (ax_br, "(b)  The twelve artifacts, per arm: rank",
             lambda a: a["metrics"]["effective_rank_residualised"], True,
             P.axis_label("R1_short", "resid")),
            (ax_bc, "        ... and information, same rows",
             lambda a: a["points"]["untrained40"]["top_cca"], False,
             "held-out top-CCA, 40 untrained targets\n"
             "confound-residualised block\nlinear scale")):
        ax.set_title(title, pad=8)
        for k, (tag, arm_name, key, _which) in enumerate(ARMS):
            y = len(ARMS) - 1 - k
            vals = [get(a) for a in art[tag]]
            ax.plot([min(vals), max(vals)], [y, y], color=C.RULE, linewidth=6,
                    solid_capstyle="round", zorder=1)
            for j, (s, v) in enumerate(zip(SEEDS, vals)):
                ax.plot([v], [y], marker=P.MARKERS[j], color=P.ORDER[j], markersize=6.5,
                        markeredgecolor="white", markeredgewidth=0.9, zorder=3,
                        label=f"seed {s}" if k == 0 else None)
            fold = max(vals) / min(vals)
            ax.annotate(f"{fold:.3f}$\\times$", (max(vals), y), xytext=(6, 0),
                        textcoords="offset points", ha="left", va="center",
                        fontsize=6.8, color=C.INK, fontweight="bold")
        ax.set_yticks(range(len(ARMS)))
        ax.set_yticklabels([f"{t}\n{n}" for t, n, _, _ in reversed(ARMS)], fontsize=6.5)
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.6, len(ARMS) - 0.25)
        if logscale:
            ax.set_xscale("log")
            ax.set_xlim(5.0, 60.0)
            ax.set_xticks([6, 10, 20, 34])
            ax.set_xticklabels(["6", "10", "20", "34"])
            ax.minorticks_off()
        else:
            ax.set_xlim(0.44, 0.70)
            ax.set_xticks([0.45, 0.50, 0.55, 0.60, 0.65])
        ax.set_xlabel(xlab, fontsize=6.4, labelpad=6)
        P.grid(ax, axis="x")
    ax_bc.set_yticklabels([])
    ax_br.legend(loc="lower left", bbox_to_anchor=(0.0, -0.02), ncol=3, handletextpad=0.3,
                 columnspacing=1.0)
    ax_br.text(0.0, -0.36,
               "The two columns side by side are the argument: within an arm, nothing about the "
               "objective, data, split,\nschedule or architecture changed - only the training seed.",
               transform=ax_br.transAxes, ha="left", va="top", fontsize=6.4, color=C.MUTED)

    # ------------------------------------------------------------------ (c)
    h44, h43 = P.load_json("ws_p2/out/P2_METRICS_D2.json")["H44"], \
        P.load_json("ws_p2/out/P2_METRICS_D2.json")["H43"]
    r44 = h44["metrics"]["effective_rank_residualised"]
    r43 = h43["metrics"]["effective_rank_residualised"]
    c44 = h44["points"]["untrained40"]["top_cca"]
    c43 = h43["points"]["untrained40"]["top_cca"]
    ax_c.set_title("(c)  The single cleanest object, enlarged: D2 arm H, seed 44 against seed 43\n"
                   "      one arm, no arm contrast to argue about, inside RankMe's reserved scope "
                   "by RankMe's own sentence", pad=8)
    ax_c.axis("off")

    sub_l = ax_c.inset_axes([0.02, 0.10, 0.42, 0.62])
    sub_r = ax_c.inset_axes([0.56, 0.10, 0.42, 0.62])
    for sub, vals, labels, xlab, fmt, logscale in (
            (sub_l, (r44, r43), ("seed 44", "seed 43"), P.axis_label("R1_short", "resid"),
             "{:.3f}", True),
            (sub_r, (c44, c43), ("seed 44", "seed 43"),
             "held-out top-CCA, 40 untrained targets\nconfound-residualised block", "{:.4f}", False)):
        lower = min(range(2), key=lambda i: vals[i])
        for j, (v, lab) in enumerate(zip(vals, labels)):
            side = -1 if j == lower else +1
            ha = "right" if side < 0 else "left"
            sub.plot([v], [0], marker=P.MARKERS[2 - j], color=P.ORDER[2 - j], markersize=13,
                     markeredgecolor="white", markeredgewidth=1.2, zorder=3)
            sub.annotate(fmt.format(v), (v, 0), xytext=(10 * side, 13),
                         textcoords="offset points", ha=ha, fontsize=8.0,
                         color=P.ORDER[2 - j], fontweight="bold")
            sub.annotate(lab, (v, 0), xytext=(10 * side, -14), textcoords="offset points",
                         ha=ha, va="top", fontsize=6.8, color=C.INK)
        sub.plot(list(vals), [0, 0], color=C.RULE, linewidth=4, zorder=1, solid_capstyle="round")
        sub.set_yticks([])
        sub.set_ylim(-0.5, 0.9)
        sub.spines["left"].set_visible(False)
        sub.set_xlabel(xlab, fontsize=6.3, labelpad=4)
        if logscale:
            sub.set_xscale("log")
            sub.set_xlim(7.5, 35)
            sub.set_xticks([8, 10, 15, 20, 30])
            sub.set_xticklabels(["8", "10", "15", "20", "30"])
            sub.minorticks_off()
        else:
            sub.set_xlim(0.5885, 0.6075)
            sub.set_xticks([0.590, 0.594, 0.598, 0.602, 0.606])
        P.grid(sub, axis="x")
    sub_l.text(0.5, 1.10, f"{r43 / r44:.2f}$\\times$ apart in rank", transform=sub_l.transAxes,
               ha="center", va="bottom", fontsize=8.5, color=C.INK, fontweight="bold")
    sub_r.text(0.5, 1.10, f"{abs(c44 - c43):.4f} apart in channel -\n"
                          "with the LOWER-rank run marginally ahead",
               transform=sub_r.transAxes, ha="center", va="bottom", fontsize=7.6,
               color=C.INK, fontweight="bold")

    P.caption(fig,
        "F2. Four arms and three seeds is 8 within-arm degrees of freedom, not a large design. "
        f"F(3,8) = {dec[0][2]['F']:.2f} for canonical effective rank is a FAILURE TO REJECT: the arm effect on rank is not resolvable here, "
        "which is not the same as there being none. The seed changes nothing about objective, architecture, data, split or schedule - "
        "that is what makes it the nuisance factor. Statistic R1 (Roy & Vetterli order-1, centred), block confound-residualised; "
        "the RankMe row is uncentred and RAW, per its own definition, and is included because it is the answer to \"you evaluated a centred variant\" - "
        "it is WORSE than ours, not better. Rank-type quantities are decomposed on the log scale and the CCA on the raw scale, "
        "as the axis label states; no two rank statistics share an axis. Panel (c) is inside RankMe's own reserved scope: three seeds of one arm "
        "are \"different runs of a given method\". Values recomputed here from the per-artifact metrics and asserted equal to the printed "
        "decomposition in data/ws_p2/out/p2_run.log.",
        y=0.006)

    P.save(fig, "F2_variance_decomposition")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
