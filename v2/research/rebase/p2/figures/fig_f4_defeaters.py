"""F4 - Defeater check: the instability is in training, not in estimation. (P2 draft 4.4)

The Leavitt & Morcos (ICLR 2021) 4.2 analogue, which a referee will ask for.

Claim: the result is not an artifact of a poor rank estimator. Four independent
measurements place the variance in training.

STATISTICS DO NOT MIX. Panel (a) and (d) are R1 on held-out artifacts; panel (c)
is R3 on a training probe. They are drawn on separate axes with separate labels
and must not be read across.

Sources
  (a) data/ws_p2/out/P2_METRICS_{D2,D1}.json  (`subsample` block: 40 draws, 80% of patients)
      cross-checked against the printed table in data/ws_p2/out/p2_run.log
  (b) v2/research/rebase/nature/D2_RESULT.md 4  +  the two D2 readout JSONs
  (c) data/e0_run/d1_diag/probevar_m{0.999,0}_{1,2,3}.log
  (d) data/ws_p2/out/P2_METRICS_D2.json  (faithful RankMe on our own artifacts)
      NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md
        sections 2 and 4 - the 68 artifact x block ratio table and the tolerance-change
        result exist ONLY as that entry's prose and table, and are parsed from it here.
"""
from __future__ import annotations

import math
import re

import p2fig as P
from p2fig import C, np, plt

RANK_ENTRY = ("NOTEBOOK_ENTRIES/"
              "effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md")

PAIRS = [
    ("D2 s42", "d2", "H42", "I42", "Hallmark", "PBS"),
    ("D2 s43", "d2", "H43", "I43", "Hallmark", "PBS"),
    ("D2 s44", "d2", "H44", "I44", "Hallmark", "PBS"),
    ("D1 s42", "d1", "P42", "F42", "programme_only", "programme_free"),
    ("D1 s43", "d1", "P43", "F43", "programme_only", "programme_free"),
    ("D1 s44", "d1", "P44", "F44", "programme_only", "programme_free"),
]


def subsample_pairs() -> list[dict]:
    src = {"d2": P.load_json("ws_p2/out/P2_METRICS_D2.json"),
           "d1": P.load_json("ws_p2/out/P2_METRICS_D1.json")}
    log = P.load_text("ws_p2/out/p2_run.log")
    printed = {}
    for m in re.finditer(
            r"^\s*(D[12] s4[234]):\s+\S+=\s*([\d.]+)\+-\s*([\d.]+)\s+\S+=\s*([\d.]+)"
            r"\+-\s*([\d.]+)\s+gap=\s*([\d.]+)\s+gap/sd=\s*([\d.]+)\s*$", log, re.M):
        printed.setdefault(m.group(1), []).append(
            tuple(float(x) for x in m.groups()[1:]))

    out = []
    for label, which, a, b, name_a, name_b in PAIRS:
        sa = src[which][a]["subsample"]["effective_rank_residualised"]
        sb = src[which][b]["subsample"]["effective_rank_residualised"]
        gap = abs(sa["mean"] - sb["mean"])
        sd = math.hypot(sa["sd"], sb["sd"])          # the printed table's definition
        want = printed[label][0]                      # first block = effective_rank_residualised
        for got, exp, tol, what in ((sa["mean"], want[0], 6e-4, f"{a} mean"),
                                    (sa["sd"], want[1], 6e-4, f"{a} sd"),
                                    (sb["mean"], want[2], 6e-4, f"{b} mean"),
                                    (sb["sd"], want[3], 6e-4, f"{b} sd"),
                                    (gap, want[4], 6e-4, f"{label} gap"),
                                    (gap / sd, want[5], 6e-3, f"{label} gap/sd")):
            assert abs(got - exp) <= tol, (
                f"F4(a): {what} is {got!r} in the metrics JSON but {exp!r} in "
                "data/ws_p2/out/p2_run.log. Refusing to plot either.")
        out.append({"label": label, "a": sa, "b": sb, "name_a": name_a, "name_b": name_b,
                    "gap": gap, "sd": sd, "gap_sd": gap / sd, "n": sa["n"], "frac": sa["frac"]})
    return out


def probe_repeats() -> dict[str, list[float]]:
    """Final R3 of each 200-step controlled probe repeat, from its own log."""
    out: dict[str, list[float]] = {}
    for m in ("0.999", "0"):
        vals = []
        for rep in (1, 2, 3):
            text = P.load_text(f"e0_run/d1_diag/probevar_m{m}_{rep}.log")
            hit = re.search(r"^DONE momentum=(\S+) decorr=\S+ final_eff_rank=([\d.]+)",
                            text, re.M)
            assert hit, f"no DONE line in probevar_m{m}_{rep}.log"
            assert float(hit.group(1)) == float(m), (m, hit.group(1))
            vals.append(float(hit.group(2)))
        out[m] = vals
    return out


def ratio_table() -> dict[str, tuple[float, float, float]]:
    """The 68-artifact x block ratio table. Lives only in the notebook entry."""
    text = P.repo_text(RANK_ENTRY)
    rows = {}
    for r in P.md_rows(text, contains="/ R1", n_cols=4):
        rows[r[0]] = (P.one_num(r[1]), P.one_num(r[2]), P.one_num(r[3]))
    assert {"R2 / R1", "R3 / R1", "R1 uncentred / R1"} <= set(rows), rows
    m = re.search(r"over all (\d+) recomputed", text)
    assert m, "cannot find the combination count in the notebook entry"
    rows["_n"] = (float(m.group(1)), 0.0, 0.0)
    tol = re.search(r"the new relative one, over all \d+ recomputed\s*\nartifact . block "
                    r"combinations, is \*\*([\d.e+-]+)\*\*", text)
    assert tol, "cannot find the tolerance-change result in the notebook entry"
    rows["_tol"] = (float(tol.group(1)), 0.0, 0.0)
    return rows


def main() -> int:
    P.cli(__doc__)
    subs = subsample_pairs()
    probes = probe_repeats()
    ratios = ratio_table()
    n_comb = int(ratios["_n"][0])

    d2 = P.load_json("ws_p2/out/P2_METRICS_D2.json")
    rec = P.load_json("e0_run/d2_v3/RECOVERED_SEED42_READOUT.json")["H42_recovered"]
    md = P.repo_text("v2/research/rebase/nature/D2_RESULT.md")
    recorded = P.one_num(P.md_rows(md, contains="recorded in the lost run", n_cols=3)[0][1])
    reexport = rec["points"]["unrestricted"]["top_cca"]

    fig = plt.figure(figsize=(7.5, 8.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.30, 0.66], hspace=0.50, wspace=0.42,
                          left=0.135, right=0.985, top=0.935, bottom=0.255)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[0, 1])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # ------------------------------------------------------------------ (a)
    ax_a.set_title("(a)  Sampling noise with the trained model held FIXED\n"
                   f"      {subs[0]['n']} draws of {subs[0]['frac']:.0%} of patients per artifact", pad=9)
    for k, s in enumerate(subs):
        y = len(subs) - 1 - k
        unresolvable = s["gap_sd"] < 3.0
        for j, (side, colour, mk) in enumerate(((s["a"], C.BLUE, "o"), (s["b"], C.VERMILLION, "s"))):
            off = 0.15 if j == 0 else -0.15
            ax_a.errorbar([side["mean"]], [y + off], xerr=[side["sd"]], fmt=mk, color=colour,
                          markersize=5.0, elinewidth=1.4, capsize=2.4,
                          markeredgecolor="white", markeredgewidth=0.8, zorder=3,
                          label=(s["name_a"] if j == 0 else s["name_b"]) if k == 0 else None)
        ax_a.annotate(f"gap/sd  {s['gap_sd']:.1f}" + ("   UNRESOLVABLE" if unresolvable else ""),
                      (max(s["a"]["mean"], s["b"]["mean"]), y), xytext=(9, 0),
                      textcoords="offset points", ha="left", va="center", fontsize=6.4,
                      color=C.VERMILLION if unresolvable else C.INK,
                      fontweight="bold" if unresolvable else "normal")
    ax_a.set_yticks(range(len(subs)))
    ax_a.set_yticklabels([s["label"] for s in reversed(subs)], fontsize=7.0)
    ax_a.tick_params(axis="y", length=0)
    ax_a.set_ylim(-0.7, len(subs) - 0.3)
    ax_a.set_xscale("log")
    ax_a.set_xlim(5.0, 95)
    ax_a.set_xticks([6, 10, 20, 30])
    ax_a.set_xticklabels(["6", "10", "20", "30"])
    ax_a.minorticks_off()
    ax_a.set_xlabel(P.axis_label("R1_short", "resid") + "\nlog scale;  bars are +-1 sd over the draws",
                    fontsize=6.4, labelpad=6)
    P.grid(ax_a, axis="x")
    ax_a.legend(loc="lower right", bbox_to_anchor=(1.0, -0.03), ncol=1, handletextpad=0.4)
    # The D2 s44 annotation is drawn beside its own row, in the same colour as the
    # marker pair it qualifies; the sentence that explains what it costs the paper
    # is in the caption, where it cannot collide with a neighbouring panel.

    # ------------------------------------------------------------------ (c)
    ax_c.set_title("(c)  Fixed seed, fixed 200-step horizon, identical inputs\n"
                   "      three repeats of a controlled probe", pad=9)
    lo_gap, hi_gap = max(probes["0"]), min(probes["0.999"])
    ax_c.axhspan(lo_gap, hi_gap, color=C.ENVELOPE, alpha=0.30, linewidth=0, zorder=0)
    for j, (m, colour, mk, name) in enumerate((("0.999", C.BLUE, "o", "momentum $m$ = 0.999"),
                                               ("0", C.VERMILLION, "s", "momentum $m$ = 0"))):
        vals = probes[m]
        xs = [1, 2, 3]
        ax_c.plot(xs, vals, linestyle="none", marker=mk, color=colour, markersize=7,
                  markeredgecolor="white", markeredgewidth=1.0, zorder=3, label=name)
        spread = (max(vals) - min(vals)) / (sum(vals) / len(vals))
        ax_c.annotate(f"spread {spread:.0%}", (3, max(vals)), xytext=(10, 0),
                      textcoords="offset points", ha="left", va="center",
                      fontsize=7.0, color=colour, fontweight="bold")
        for x, v in zip(xs, vals):
            ax_c.annotate(f"{v:.2f}", (x, v), xytext=(0, 8), textcoords="offset points",
                          ha="center", fontsize=6.3, color=colour)
    ax_c.text(0.5, (lo_gap + hi_gap) / 2, f"empty band\n{lo_gap:.2f} to {hi_gap:.2f}",
              ha="left", va="center", fontsize=6.4, color=C.MUTED)
    ax_c.set_xticks([1, 2, 3])
    ax_c.set_xticklabels(["repeat 1", "repeat 2", "repeat 3"], fontsize=6.8)
    ax_c.set_xlim(0.45, 3.9)
    ax_c.set_ylim(0, 9.2)
    ax_c.set_ylabel("rank at step 200\n" + P.STAT["R3_short"] + "\n" + P.BLOCK["probe"]
                    + ", 256 patients", fontsize=6.4)
    P.grid(ax_c, axis="y")
    ax_c.legend(loc="upper left", bbox_to_anchor=(0.0, 0.99), handletextpad=0.4)

    # ------------------------------------------------------------------ (b)
    ax_b.set_title("(b)  Readout determinism", pad=9)
    for j, (v, lab, colour, mk) in enumerate((
            (recorded, "recorded in the lost run", C.MUTED, "x"),
            (reexport, "re-export of the surviving\nseed-42 checkpoint", C.BLUE, "o"))):
        ax_b.plot([v], [j], marker=mk, color=colour, markersize=10,
                  markeredgecolor="white", markeredgewidth=1.0, zorder=3)
        ax_b.annotate(f"{v:.5f}".rstrip("0"), (v, j), xytext=(0, 11), textcoords="offset points",
                      ha="center", fontsize=7.6, color=colour, fontweight="bold")
        ax_b.text(0.58635, j, lab, ha="left", va="center", fontsize=6.5, color=C.INK)
    ax_b.set_yticks([])
    ax_b.spines["left"].set_visible(False)
    ax_b.set_ylim(-0.7, 1.8)
    ax_b.set_xlim(0.5854, 0.5872)
    ax_b.set_xticks([0.5855, 0.5860, 0.5865, 0.5870])
    ax_b.set_xlabel("unrestricted held-out top-CCA, 90 targets\n"
                    "same checkpoint, same arguments, re-exported", fontsize=6.4, labelpad=6)
    P.grid(ax_b, axis="x")
    ax_b.text(0.0, -0.42, f"deterministic to five significant figures\n"
                          f"({reexport:.5f} against {recorded:.4f} recorded).\n"
                          "The export and readout path is not the source of the variance.",
              transform=ax_b.transAxes, ha="left", va="top", fontsize=6.4, color=C.INK)

    # ------------------------------------------------------------------ (d)
    ax_d.set_title("(d)  Definitional insensitivity", pad=9)
    r1 = d2["H42"]["metrics"]["effective_rank_residualised"]
    rm = d2["H42"]["metrics"]["rankme_residualised"]
    for j, (v, lab, colour, mk) in enumerate(((r1, "canonical R1", C.BLUE, "o"),
                                              (rm, "faithful RankMe\n(uncentred, $\\varepsilon = 10^{-7}$)",
                                               C.GREEN, "^"))):
        ax_d.plot([v], [j], marker=mk, color=colour, markersize=10,
                  markeredgecolor="white", markeredgewidth=1.0, zorder=3)
        ax_d.annotate(f"{v:.3f}", (v, j), xytext=(0, 11), textcoords="offset points",
                      ha="center", fontsize=7.6, color=colour, fontweight="bold")
        ax_d.annotate(lab, (v, j), xytext=(0, -12), textcoords="offset points",
                      ha="center", va="top", fontsize=6.5, color=C.INK)
    ax_d.set_yticks([])
    ax_d.spines["left"].set_visible(False)
    ax_d.set_ylim(-1.05, 1.75)
    ax_d.set_xlim(23.3838, 23.3942)
    ax_d.set_xticks([23.385, 23.387, 23.389, 23.391, 23.393])
    ax_d.set_xlabel("effective rank on H42, confound-residualised block\n"
                    "two definitions, one artifact", fontsize=6.4, labelpad=6)
    P.grid(ax_d, axis="x")
    unc = ratios["R1 uncentred / R1"]
    ax_d.text(0.0, -0.42,
              f"Over all {n_comb} recomputed artifact x block combinations:\n"
              f"  absolute -> relative tolerance change moved NO value "
              f"(max relative difference {ratios['_tol'][0]:.3e});\n"
              f"  uncentred R1 / centred R1 has median {unc[1]:.3f} "
              f"(min {unc[0]:.3f}, max {unc[2]:.3f}).\n"
              "Source: NOTEBOOK_ENTRIES/effective_rank_canonicalised_..._20260804T0005Z.md, 2 and 4.",
              transform=ax_d.transAxes, ha="left", va="top", fontsize=6.3, color=C.INK)

    P.caption(fig,
        "F4. In panel (a), D2 s44 is the ONE unresolvable pair - the between-arm gap is 1.4 sampling standard deviations - and section 4.6's D2 "
        "\"2 of 3\" depends on it. It is not a hit; it is a pair the metric cannot separate, which is why T1 marks that cell unresolvable in every row. "
        "Panel (c) constrains the TYPICAL spread and essentially not the tail: three repetitions cannot rule out the liveness gate's 25% divergence rate "
        "(P(0 in 3 | p = 0.25) = 0.42; the exact upper 95% bound from 0 of 3 is p <= 0.63), and the design was cut from ten repetitions to three because a "
        "ten-way launch exhausted GPU memory. Panel (c)'s statistic is R3 on a training probe and must not share an axis with (a)'s and (d)'s R1 on held-out "
        "artifacts; they are drawn on separate axes for that reason. Panel (a) and F3 are not in conflict: (a) holds the model fixed and varies the patient "
        "sample, F3 holds the step fixed and varies the seed. Panel (b) shares its data with F1(a). The faithful RankMe of panel (d) is computed on OUR OWN "
        "artifacts and is an internal comparison; no number in this paper is plotted against a published RankMe value, because RankMe's epsilon sits outside "
        "the division and its p_k do not sum to 1.",
        y=0.006)

    P.save(fig, "F4_defeater_check")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
