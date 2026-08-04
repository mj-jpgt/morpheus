"""F6 - The necessity test, which went against us.  (P2 draft 4.7)

PLACED BEFORE F7 AND F8, mirroring the draft's order, because this is the result
that falsified the previous framing and the draft reports it before the instances
that favour the paper. It is not supplementary and it is not drawn smaller than
F1 or F2.

Claim: RankMe's necessary-not-sufficient hedge is NOT violated by our
best-matched three-seed experiment. It is confirmed, 3 of 3 - and the gaps at
which it is confirmed are themselves inside rank's own nuisance band.

BOTH ESTIMATORS ARE DRAWN AND THE CONSERVATIVE ONE IS WEIGHTED. The patient
bootstrap is decisive 3 of 3; the cancer-cluster bootstrap is 2 of 3, and seed
43's cancer interval grazes zero at +0.0006. A panel showing only the patient
interval would be the selective quotation section 4.6 exists to refuse, so the
cancer interval is drawn thicker, listed first in the legend, and its grazing
endpoint is annotated.

D1_PAIRED_BOOTSTRAP.json (unsuffixed) MUST NOT BE USED - it scores all 90
non-control targets, 50 of which are programme_only's own supervision. It is not
vendored into this repository for that reason; only the STRATIFIED (40 targets)
and RANDOM_CONTROL (90 targets) files are.

Sources
  data/ws_p2/out/P2_METRICS_D1.json                        rank and channel points
  data/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json    40 targets, 2,000 resamples
  data/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json 90 random_control targets
  data/e0_run/d1_v2/D1_PAIR_MANIFEST.json                  matching + the preregistered escalation
  data/ws_p2/out/p2_run.log                                the 66-pair violation scan
  v2/research/rebase/p2/p2_necessity_and_variance.py        the predeclared thresholds
  NOTEBOOK_ENTRIES/PREDECLARED_D1_necessity_test_20260803T2300Z.md   the outcome predeclaration
  NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md  where the
      fold / dCCA criterion is declared "fixed before inspecting the pairs"
"""
from __future__ import annotations

import re

import p2fig as P
from p2fig import C, np, plt

PREDECL = "NOTEBOOK_ENTRIES/PREDECLARED_D1_necessity_test_20260803T2300Z.md"
CRITERION = ("NOTEBOOK_ENTRIES/"
             "p2_competing_metrics_and_necessity_test_20260803T2326Z.md")
SCRIPT = "v2/research/rebase/p2/p2_necessity_and_variance.py"
SEEDS = (42, 43, 44)


def thresholds() -> tuple[float, float]:
    """The pre-declared violation criterion, read from the script that declares it."""
    src = P.repo_text(SCRIPT)
    fold = float(re.search(r"^RANK_FOLD\s*=\s*([\d.]+)", src, re.M).group(1))
    dcca = float(re.search(r"^CCA_DELTA\s*=\s*([\d.]+)", src, re.M).group(1))
    # The same criterion is stated in the notebook entry as fixed before the pairs
    # were inspected; confirm the two agree before either is drawn.
    entry = P.repo_text(CRITERION)
    quoted = re.search(r"a pair \(lo, hi\) is a violation iff.*?`\s*\*?\*?and\*?\*?\s*\n"
                       r">\s*`?CCA\(lo\).{1,3} CCA\(hi\).{1,3}\s*([\d.]+)`",
                       entry, re.S)
    assert quoted, "cannot find the quoted violation criterion in the notebook entry"
    fold_quoted = float(re.search(r"eff_rank\(hi\)/eff_rank\(lo\)\D*([\d.]+)",
                                  entry).group(1))
    assert fold_quoted == fold and float(quoted.group(1)) == dcca, (
        f"F6: the criterion declared in {SCRIPT} is fold >= {fold}, dCCA >= {dcca}, but "
        f"{CRITERION} quotes fold >= {fold_quoted}, dCCA >= {quoted.group(1)}. "
        "Refusing to draw either.")
    return fold, dcca


def d1_points() -> list[dict]:
    d1 = P.load_json("ws_p2/out/P2_METRICS_D1.json")
    strat = P.load_json("e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json")
    ctrl = P.load_json("e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json")
    assert strat["n_targets"] == 40 and ctrl["n_targets"] == 90, (strat, ctrl)
    assert ctrl["target_groups"] == ["random_control"], ctrl["target_groups"]

    out = []
    for k, seed in enumerate(SEEDS):
        p, f = d1[f"P{seed}"], d1[f"F{seed}"]
        sp, sc = strat["pairs"][k], ctrl["pairs"][k]
        assert f"seed{seed}" in sp["programme_only"] and f"seed{seed}" in sc["programme_only"]
        # The bootstrap's own point estimates must equal the metric file's.
        for got, want in ((sp["point_programme_only"], p["points"]["untrained40"]["top_cca"]),
                          (sp["point_programme_free"], f["points"]["untrained40"]["top_cca"])):
            assert abs(got - want) < 1e-9, (
                f"F6: seed {seed} channel is {want!r} in P2_METRICS_D1.json but {got!r} in "
                "D1_PAIRED_BOOTSTRAP_STRATIFIED.json. Refusing to plot either.")
        d = sp["programme_free_minus_programme_only"]
        c = sc["programme_free_minus_programme_only"]
        out.append({
            "seed": seed,
            "rank_p": p["metrics"]["effective_rank_residualised"],
            "rank_f": f["metrics"]["effective_rank_residualised"],
            "cca_p": sp["point_programme_only"], "cca_f": sp["point_programme_free"],
            "delta": d["patient"]["point_delta"],                 # signed free - only
            "patient": (d["patient"]["ci95_low"], d["patient"]["ci95_high"]),
            "cancer": (d["cancer"]["ci95_low"], d["cancer"]["ci95_high"]),
            "ctrl_delta": c["patient"]["point_delta"],
            "ctrl_patient": (c["patient"]["ci95_low"], c["patient"]["ci95_high"]),
            "ctrl_cancer": (c["cancer"]["ci95_low"], c["cancer"]["ci95_high"]),
            "repeats": d["patient"]["repeats"],
        })
    return out


def nuisance_band() -> tuple[float, float]:
    """Section 4.2's own within-arm seed folds: the band F6(b) draws behind its bars."""
    d2 = P.load_json("ws_p2/out/P2_METRICS_D2.json")
    d1 = P.load_json("ws_p2/out/P2_METRICS_D1.json")
    folds = []
    for src, keys in ((d2, "HI"), (d1, "PF")):
        for key in keys:
            vals = [src[f"{key}{s}"]["metrics"]["effective_rank_residualised"] for s in SEEDS]
            folds.append(max(vals) / min(vals))
    return min(folds), max(folds)


def violations() -> list[dict]:
    """The 66-pair scan's hits, parsed from the run log that produced them."""
    log = P.load_text("ws_p2/out/p2_run.log")
    n = int(re.search(r"(\d+) violating pair\(s\) of (\d+) examined", log).group(2))
    hits = []
    for m in re.finditer(
            r"^\s*([A-Z]\d\d)\s+([\d.]+)\s+([\d.]+)\s+\|\s+([A-Z]\d\d)\s+([\d.]+)\s+([\d.]+)"
            r"\s+([\d.]+)x\s+\+([\d.]+)\s+(\S.*?)\s*$", log, re.M):
        hits.append({"lo": m.group(1), "lo_rank": float(m.group(2)), "lo_cca": float(m.group(3)),
                     "hi": m.group(4), "hi_rank": float(m.group(5)), "hi_cca": float(m.group(6)),
                     "fold": float(m.group(7)), "dcca": float(m.group(8)), "scope": m.group(9)})
    assert len(hits) == 2, hits
    return hits, n


def main() -> int:
    P.cli(__doc__)
    fold_thr, dcca_thr = thresholds()
    pts = d1_points()
    band_lo, band_hi = nuisance_band()
    hits, n_examined = violations()
    manifest = P.load_json("e0_run/d1_v2/D1_PAIR_MANIFEST.json")
    assert manifest["objective_only_difference"] is True, manifest

    fig = plt.figure(figsize=(7.5, 10.9))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.05, 0.70], hspace=0.78, wspace=0.34,
                          left=0.135, right=0.985, top=0.960, bottom=0.185)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    ax_d = fig.add_subplot(gs[2, 0])
    ax_leg = fig.add_subplot(gs[2, 1])

    # ------------------------------------------------------------------ (a)
    ax_a.set_title("(a)  The predeclaration, drawn", pad=9)
    xmax, ymin, ymax = 4.4, -0.13, 0.16
    ax_a.add_patch(plt.Rectangle((fold_thr, dcca_thr), xmax - fold_thr, ymax - dcca_thr,
                                 facecolor=C.VERMILLION, alpha=0.13, edgecolor=C.VERMILLION,
                                 linewidth=1.0, linestyle=(0, (4, 2)), zorder=1))
    ax_a.text((fold_thr + xmax) / 2, (dcca_thr + ymax) / 2,
              "VIOLATION\nthe only configuration\nthe hedge cannot absorb\n"
              f"fold $\\geq$ {fold_thr}  AND  $\\Delta$CCA $\\geq$ +{dcca_thr}",
              ha="center", va="center", fontsize=6.3, color=C.VERMILLION,
              fontweight="bold", linespacing=1.6, zorder=2)
    ax_a.axhline(0.0, color=C.INK, linewidth=0.8, zorder=1)
    for k, p in enumerate(pts):
        fold = p["rank_p"] / p["rank_f"]
        ax_a.plot([fold], [p["delta"]], marker=P.MARKERS[k], color=C.BLUE, markersize=8,
                  markeredgecolor="white", markeredgewidth=1.0, zorder=4,
                  label=f"D1 seed {p['seed']}")
        ax_a.annotate(f"s{p['seed']}", (fold, p["delta"]), xytext=(0, -13),
                      textcoords="offset points", ha="center", va="top",
                      fontsize=6.5, color=C.BLUE)
    ax_a.set_xlim(1.0, xmax)
    ax_a.set_ylim(ymin, ymax)
    ax_a.set_xlabel("rank fold, programme_only / programme_free\n"
                    + P.STAT["R1_short"] + ", residualised block", fontsize=6.4, labelpad=6)
    ax_a.set_ylabel("channel difference\nprogramme_free $-$ programme_only\n"
                    "held-out top-CCA, 40 untrained targets", fontsize=6.4)
    P.grid(ax_a)
    ax_a.text(0.0, -0.34,
              "All three land in the CONFIRMING quadrant: lower rank, LESS information.\n"
              "Both thresholds, and outcomes O1-O4, were pre-declared: see the caption.",
              transform=ax_a.transAxes, ha="left", va="top", fontsize=6.2,
              color=C.INK, linespacing=1.6)

    # ------------------------------------------------------------------ (b)
    ax_b.set_title("(b)  The three seeds, in rank", pad=9)
    ax_b.axhspan(band_lo, band_hi, color=C.ENVELOPE, alpha=0.35, linewidth=0, zorder=0)
    ax_b.text(0.99, 0.995,
              f"shaded: section 4.2's own WITHIN-ARM seed band,\n"
              f"{band_lo:.2f}-{band_hi:.2f}$\\times$ (F2(b), four arms)",
              transform=ax_b.transAxes, ha="right", va="top", fontsize=6.2,
              color=C.MUTED, linespacing=1.5)
    xs = np.arange(len(pts))
    for k, p in enumerate(pts):
        fold = p["rank_p"] / p["rank_f"]
        inside = band_lo <= fold <= band_hi
        ax_b.bar(k, fold, width=0.5, color=C.BLUE if inside else C.VERMILLION,
                 edgecolor="white", linewidth=1.6, zorder=2,
                 hatch="" if inside else "///")
        ax_b.annotate(f"{fold:.2f}$\\times$", (k, fold), xytext=(0, 5),
                      textcoords="offset points", ha="center", fontsize=7.2,
                      color=C.INK, fontweight="bold")
        ax_b.annotate(f"{p['rank_p']:.3f}\n{p['rank_f']:.3f}", (k, 0.10),
                      xytext=(0, 0), textcoords="offset points", ha="center", va="bottom",
                      fontsize=6.0, color=C.INK, zorder=5,
                      bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.92))
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels([f"seed {p['seed']}" for p in pts], fontsize=7.0)
    ax_b.set_xlim(-0.6, len(pts) - 0.4)
    ax_b.set_ylim(0, 4.9)
    ax_b.set_ylabel("rank fold, programme_only / programme_free\n"
                    + P.STAT["R1_short"] + ", residualised block", fontsize=6.4)
    P.grid(ax_b, axis="y")
    n_inside = sum(band_lo <= p["rank_p"] / p["rank_f"] <= band_hi for p in pts)
    ax_b.text(0.0, -0.34,
              f"{n_inside} of 3 sit inside the nuisance band; all three are below the largest\n"
              "within-arm seed fold this project has measured. Bar labels are the two R1 values.",
              transform=ax_b.transAxes, ha="left", va="top", fontsize=6.2, color=C.INK,
              linespacing=1.6)

    # -------------------------------------------------------------- (b)+(c)
    ax_c.set_title("(b)  ... and in the information channel, with BOTH bootstrap estimators - "
                   "and (c)  the negative control on the same axes", pad=9)
    ax_c.axvline(0.0, color=C.INK, linewidth=1.0, zorder=2)
    rows = []
    for p in pts:
        rows.append((f"seed {p['seed']}\n40 real targets", p["delta"], p["patient"], p["cancer"],
                     False))
    for p in pts:
        rows.append((f"seed {p['seed']}\n90 random_control", p["ctrl_delta"], p["ctrl_patient"],
                     p["ctrl_cancer"], True))
    for k, (lab, delta, pat, can, is_ctrl) in enumerate(rows):
        y = len(rows) - 1 - k
        base = C.FAINT if is_ctrl else None
        for name, (clo, chi), colour, mk, off, lw in (
                ("cancer-cluster bootstrap (conservative)", can,
                 base or C.VERMILLION, "s", -0.17, 3.4),
                ("patient bootstrap", pat, base or C.BLUE, "o", +0.17, 1.7)):
            ax_c.plot([clo, chi], [y + off, y + off], color=colour, linewidth=lw,
                      solid_capstyle="butt", zorder=3,
                      label=name if k == 0 else None)
            ax_c.plot([delta], [y + off], marker=mk, color=colour, markersize=5.0,
                      markeredgecolor="white", markeredgewidth=0.9, zorder=4)
            if not is_ctrl and chi > 0:
                ax_c.annotate(f"grazes zero at {chi:+.4f}", (chi, y + off), xytext=(7, 0),
                              textcoords="offset points", ha="left", va="center",
                              fontsize=6.2, color=C.VERMILLION, fontweight="bold")
        ax_c.annotate(f"{delta:+.4f}", (delta, y + 0.40), ha="center", va="bottom",
                      fontsize=6.6, color=C.MUTED if is_ctrl else C.INK, fontweight="bold")
    ax_c.axhline(2.5, color=C.RULE, linewidth=0.8)
    ax_c.set_yticks(range(len(rows)))
    ax_c.set_yticklabels([r[0] for r in reversed(rows)], fontsize=6.6)
    ax_c.tick_params(axis="y", length=0)
    ax_c.set_ylim(-0.75, len(rows) - 0.25)
    ax_c.set_xlim(-0.175, 0.075)
    ax_c.set_xlabel("paired channel difference, programme_free $-$ programme_only\n"
                    "held-out top-CCA;  CI$_{95}$ from "
                    f"{pts[0]['repeats']:,} resamples;  confound-residualised block",
                    fontsize=6.4, labelpad=6)
    P.grid(ax_c, axis="x")
    ax_c.legend(loc="upper left", bbox_to_anchor=(0.0, -0.155), ncol=2, handlelength=2.0)
    dec_pat = sum(p["patient"][1] < 0 for p in pts)
    dec_can = sum(p["cancer"][1] < 0 for p in pts)
    ax_c.text(0.985, 0.985,
              f"patient bootstrap decisive {dec_pat} of 3;  "
              f"cancer-cluster bootstrap {dec_can} of 3.\n"
              "The conservative estimator is drawn thicker and named first.",
              transform=ax_c.transAxes, ha="right", va="top", fontsize=6.3,
              color=C.INK, linespacing=1.6)

    # ------------------------------------------------------------------ (d)
    ax_d.set_title("(d)  The one violation, small and last", pad=9)
    ax_d.axis("off")
    lines = [f"{'lower-rank':<8}{'rank':>9}{'CCA':>8}   {'higher-rank':<12}{'rank':>9}{'CCA':>8}"
             f"{'fold':>8}{'dCCA':>9}",
             ""]
    for h in hits:
        lines.append(f"{h['lo']:<8}{h['lo_rank']:>9.3f}{h['lo_cca']:>8.4f}   "
                     f"{h['hi']:<12}{h['hi_rank']:>9.3f}{h['hi_cca']:>8.4f}"
                     f"{h['fold']:>7.2f}x{h['dcca']:>+9.4f}")
        lines.append(f"{'':<8}{h['scope']}")
    ax_d.text(0.0, 0.98, "\n".join(lines), transform=ax_d.transAxes, va="top", ha="left",
              family="monospace", fontsize=5.3, color=C.INK, linespacing=1.9)
    ax_d.text(0.0, 0.10,
              f"{len(hits)} of {n_examined} ordered pairs. Drawn deliberately small: both are\n"
              "cross-arm, one is cross-experiment, and RankMe reserves itself to\n"
              "\"different runs of a given method\". Partially pre-empted by Aldeneh\n"
              "et al. (ICASSP 2025), who published that lower-ranked layers can\n"
              "outperform higher-ranked ones.",
              transform=ax_d.transAxes, va="top", ha="left", fontsize=6.1,
              color=C.MUTED, linespacing=1.6)

    # The preregistered escalation, carried with the figure rather than around it.
    ax_leg.axis("off")
    ax_leg.text(0.0, 0.98,
                "THE PREREGISTERED ESCALATION THIS FIGURE TRIPS\n\n"
                "D1_PAIR_MANIFEST.json records, before the run:\n"
                f"\"{manifest['preregistered_prediction']}\"\n\n"
                "programme_only won on all three seeds. This paper flags that\n"
                "and does not resolve it.\n\n"
                "The arms differ in OBJECTIVE, so D1 sits closer to a between-\n"
                "method comparison than D2 and lands further outside RankMe's\n"
                "stated scope. Section 4.7 is therefore a CONFIRMATION of\n"
                "RankMe, not a refutation of it, and this figure carries the\n"
                "paper's own negative.",
                transform=ax_leg.transAxes, va="top", ha="left", fontsize=6.1,
                color=C.INK, linespacing=1.7)

    P.caption(fig,
        "F6. This is the paper's own negative and it is reported before the instances that favour the paper. "
        "D1_PAIRED_BOOTSTRAP.json (unsuffixed) MUST NOT be used and is not vendored here: it scores all 90 non-control targets, of which 50 are "
        "programme_only's own supervision. The intervals drawn are the 40-target STRATIFIED bootstrap and the 90-target RANDOM_CONTROL bootstrap. "
        "Both estimators are drawn and the conservative one is weighted: the patient bootstrap is decisive 3 of 3, the cancer-cluster bootstrap 2 of 3, "
        "and seed 43's cancer interval grazes zero - the panel shows that it does. The arms are matched on everything but the objective "
        "(D1_PAIR_MANIFEST.json, \"objective_only_difference\": true). Statistic R1, block confound-residualised, on the held-out wsi_biology block. "
        f"PANEL (a)'S PREDECLARATION: the violation region is fold >= {fold_thr} AND dCCA >= +{dcca_thr}, both fixed before the pair list was inspected. "
        f"They are declared in {SCRIPT} and quoted as pre-fixed in {CRITERION}; the four outcomes O1-O4 and how each would be read were pre-declared "
        f"2026-08-03 23:00 UTC in {PREDECL}, before any channel value had been seen.",
        y=0.006)

    P.save(fig, "F6_necessity_test")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
