"""F3 - The floor is a property of the arm, not of the statistic.  (P2 draft 4.3)

Claim: there is no one reproducibility envelope to calibrate and reuse. At the
same global training step, on the same configuration, with only the seed
differing, one arm reproduces to a fifth and its sibling spans a factor of six.

STATISTIC AND BLOCK. Statistic R3 - the in-run tripwire, centred and
L2-row-normalised, order-2 - NOT R1, and the block is a TRAINING BATCH at
global step 200 (epoch 11), not a held-out set. These are in-training
measurements whose states were never saved, so they are [NOT RECOMPUTABLE]
under R1 and no number in this figure may appear on an axis with any R1 value
in this paper. The key is `train_rank_tripwire_observed`, WITH the `train_`
prefix: querying it without the prefix returns [] and [] reads as a confident
negative (NOTEBOOK_ENTRIES/operational_shared_box_rules_20260804T0730Z.md).

n = 5, NOT n = 3. P2_FIGURES.md quotes this panel at n = 3 with
`programme_only` spanning 1.003x. That printing predates seeds 45 and 46, which
landed in ~/e0_run/d1_seeds4546/. At n = 5 the stable arm spans 1.18x and the
unstable one 6.05x. The n = 3 subset is printed on the panel as superseded, so
a reader holding the older table can see what changed and why.

Source
  data/extracted/F3_TRIPWIRE_STEP200_R3_n5.json
    distilled by extract_from_box.py from
      ~/e0_run/d1_v2/d1_{p,f}_seed{42,43,44}/train_metrics.jsonl
      ~/e0_run/d1_seeds4546/d1_{p,f}_seed{45,46}/train_metrics.jsonl
    each input's SHA-256 is recorded inside that file.
  data/ws_rank/RANK_RECOMPUTE.json  (`instability_tripwire_step200_R3`)
    the independently written n = 3 table, used here as a cross-check.
"""
from __future__ import annotations

import p2fig as P
from p2fig import C, np, plt

ROWS = [
    ("programme_only", C.BLUE, "o", "the stable arm"),
    ("programme_free", C.VERMILLION, "s", "the interesting arm"),
]


def tripwire() -> dict:
    d = P.load_json("extracted/F3_TRIPWIRE_STEP200_R3_n5.json")
    assert d["statistic"] == "R3", d["statistic"]
    assert d["key"] == "train_rank_tripwire_observed", d["key"]

    # Cross-check the three seeds that the canonicalisation run also wrote.
    ref = P.load_json("ws_rank/RANK_RECOMPUTE.json")["instances"]["instability_tripwire_step200_R3"]
    for name, entry in ref.items():
        got = d["runs"][name]["value"]
        assert abs(got - entry["R3_at_step_200"]) < 1e-9, (
            f"F3: {name} is {got!r} in the n=5 extraction but "
            f"{entry['R3_at_step_200']!r} in ws_rank/RANK_RECOMPUTE.json. "
            "Refusing to plot either.")
        assert d["runs"][name]["epoch"] == entry["epoch"]
    return d


def main() -> int:
    P.cli(__doc__)
    d = tripwire()
    by_arm = {}
    for name, r in sorted(d["runs"].items(), key=lambda kv: kv[1]["seed"]):
        by_arm.setdefault(r["arm"], []).append((r["seed"], r["value"]))

    fig = plt.figure(figsize=(7.4, 4.55))
    ax = fig.add_axes([0.235, 0.60, 0.62, 0.245])

    folds = {}
    for k, (arm, colour, marker, gloss) in enumerate(ROWS):
        y = len(ROWS) - 1 - k
        pts = by_arm[arm]
        vals = [v for _s, v in pts]
        lo, hi = min(vals), max(vals)
        folds[arm] = hi / lo

        # The spread, drawn as a bracket with the fold printed. The two rows are
        # deliberately NOT connected: they are separate arms, not a series.
        ax.plot([lo, hi], [y, y], color=colour, linewidth=2.0, solid_capstyle="butt", zorder=2)
        for edge in (lo, hi):
            ax.plot([edge, edge], [y - 0.13, y + 0.13], color=colour, linewidth=2.0, zorder=2)
        for _seed, v in pts:
            ax.plot([v], [y], marker=marker, color=colour, markersize=7.5,
                    markeredgecolor="white", markeredgewidth=1.1, zorder=4)
        # Seed labels: points closer than 3% cannot be labelled separately at this
        # scale - and their being on top of one another IS the panel's point - so
        # they get one grouped label rather than four unreadable ones.
        groups: list[list[tuple[int, float]]] = []
        for seed, v in sorted(pts, key=lambda sv: sv[1]):
            if groups and v / groups[-1][-1][1] < 1.03:
                groups[-1].append((seed, v))
            else:
                groups.append([(seed, v)])
        for j, g in enumerate(groups):
            xs = [v for _s, v in g]
            mid = (min(xs) * max(xs)) ** 0.5
            ax.annotate(",".join(str(s) for s, _v in g), (mid, y),
                        xytext=(0, 13 + 9 * (j % 2)), textcoords="offset points",
                        ha="center", va="bottom", fontsize=5.8, color=C.MUTED,
                        arrowprops=dict(arrowstyle="-", color=C.RULE, linewidth=0.5,
                                        shrinkA=0, shrinkB=6))
        ax.annotate(f"spread {folds[arm]:.2f}$\\times$", (hi, y), xytext=(13, 0),
                    textcoords="offset points", ha="left", va="center",
                    fontsize=9.0, color=colour, fontweight="bold")

    ax.set_yticks([1, 0])
    ax.set_yticklabels([f"{a}\n{g}" for a, _c, _m, g in ROWS], fontsize=7.2)
    for lbl, (_a, colour, _m, _g) in zip(ax.get_yticklabels(), ROWS):
        lbl.set_color(colour)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.62, 1.62)
    ax.set_xscale("log")
    ax.set_xlim(5.5, 190)
    ax.set_xticks([6, 10, 20, 50, 100, 150])
    ax.set_xticklabels(["6", "10", "20", "50", "100", "150"])
    ax.minorticks_off()
    ax.set_xlabel(
        "in-run rank tripwire at global step 200 (epoch 11)  -  log scale\n"
        + P.STAT["R3"] + "\n"
        + "block: " + P.BLOCK["train"]
        + ";  key `train_rank_tripwire_observed`", fontsize=6.8, labelpad=8)
    P.grid(ax, axis="x")
    ax.set_title("F3   The reproducibility floor is a property of the ARM, not of the statistic\n"
                 "      same configuration, same global step, only the training seed differs",
                 pad=10)

    # The values, and the superseded n = 3 subset, as a printed table.
    seeds = sorted({s for pts in by_arm.values() for s, _v in pts})
    head = f"{'arm':<16}" + "".join(f"{'seed ' + str(s):>11}" for s in seeds) + f"{'spread':>11}"
    lines = [head]
    for arm, _c, _m, _g in ROWS:
        vals = dict(by_arm[arm])
        lines.append(f"{arm:<16}" + "".join(f"{vals[s]:>11.3f}" for s in seeds)
                     + f"{folds[arm]:>10.2f}x")
    sub = {}
    for arm, pts in by_arm.items():
        v3 = [v for s, v in pts if s in (42, 43, 44)]
        sub[arm] = max(v3) / min(v3)
    fig.text(0.028, 0.315, "\n".join(lines), fontsize=6.4, color=C.INK,
             family="monospace", linespacing=1.8, va="top")
    fig.text(0.028, 0.185,
             "n = 3 (seeds 42/43/44 only), as printed in paper/P2_FIGURES.md F3:  "
             f"programme_only {sub['programme_only']:.3f}x,  programme_free {sub['programme_free']:.2f}x.\n"
             "SUPERSEDED. Seeds 45 and 46 landed afterwards in ~/e0_run/d1_seeds4546/ and are included "
             f"above; the stable arm's floor is {folds['programme_only']:.2f}x, not "
             f"{sub['programme_only']:.3f}x.",
             fontsize=6.4, color=C.INK, linespacing=1.7, va="top")

    P.caption(fig,
        "F3. Statistic R3 - the in-run tripwire statistic, NOT R1 - and the block is a TRAINING batch at global step 200, not a held-out set. "
        "The states behind these numbers were never saved, so they are [NOT RECOMPUTABLE] under R1, and no number in this figure may share an axis "
        "with any R1 value in this paper. The key is `train_rank_tripwire_observed` with the `train_` prefix, because querying it without the prefix "
        "returns [] and an empty list reads as a confident negative. The two rows are not connected and share only the axis: they are two arms, not a series. "
        "The reading is that rank's reproducibility is worst exactly where the arm is interesting - which is the configuration a practitioner is looking at "
        "when they reach for rank. F3 and F4(a) are not in conflict: F4(a) holds the model fixed and varies the patient sample; F3 holds the step fixed and "
        "varies the seed.",
        y=0.008)

    P.save(fig, "F3_reproducibility_floor")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
