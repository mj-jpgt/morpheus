"""F8 - The use that survives, its boundary, and the withdrawal.  (P2 draft 4.9, 4.10)

This figure carries the paper's self-correction and must not be drawn without it.

Claim: effective rank near its floor (about 1-2, with patient-to-patient mutual
cosine near 1) is reliable evidence of TOTAL COLLAPSE. Anywhere above that -
including at 3.6% of nominal dimensionality - it is uninformative about the
channel in both directions. And a HARD matrix rank is worthless for the
surviving use.

PANEL (b) IS THE WITHDRAWAL, AND IT IS NOT OPTIONAL. A version of this panel
showing only the 16/16 pinning is the single most misleading figure this project
could publish, and this project has already described this instance incorrectly
once. The panel therefore shows, stacked: the hard `z_biology` matrix rank pinned
at 16/16 with its STRUCTURAL ceiling drawn; the CENTRED effective rank of the
SAME objective falling 12.88 -> 1.00 by step 50; and the collapse evidence for
the same arm, including a retrieval accuracy that ends BELOW its own chance line.

NEEDS EXTRACTION, discharged and reported honestly:
  * the centred effective-rank track IS a per-step array and is drawn as one -
    steps 0/25/50/100/200/400 from `~/e0_run/d1_diag/diag_d.log`.
  * the 16/16 track and the cosine / retrieval evidence are recorded as ENDPOINT
    PAIRS ONLY in `~/e0_run/collapse_diag.log`, which probes before and after
    rather than on a schedule. There is no array to extract. They are drawn as
    before/after paired markers, labelled "endpoint values as recorded; per-step
    array not retained", and NOTHING IS INTERPOLATED. The source script for the
    original 16/16 measurement is named in NOTEBOOK.md as a scratchpad file on
    the A100 and is not in this repository.

Sources
  data/e0_run/d1_diag/diag_d.log        per-step track, steps 0-400
  data/e0_run/collapse_diag.log         the "16/16" instance, endpoint pairs only
  data/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json   the boundary marker and its null
  NOTEBOOK_ENTRIES/d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md
  NOTEBOOK_ENTRIES/d1_programme_free_collapsing_in_training_20260803T1930Z.md
  NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md
"""
from __future__ import annotations

import re

import p2fig as P
from p2fig import C, np, plt

D1B = ("NOTEBOOK_ENTRIES/"
       "d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md")
D1_MID = ("NOTEBOOK_ENTRIES/"
          "d1_programme_free_collapsing_in_training_20260803T1930Z.md")
D1A_END = ("NOTEBOOK_ENTRIES/"
           "d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md")


def diag_d_track() -> dict:
    """The per-step diagnostic trace. This one IS an array."""
    text = P.load_text("e0_run/d1_diag/diag_d.log")
    rows = []
    for m in re.finditer(
            r"^\s*(\d+) loss ([\d.]+) acc ([\d.]+) pos ([\d.]+) worst-neg ([\d.]+)"
            r" min-margin (-?[\d.]+) wsi-wsi ([\d.]+) eff-rank ([\d.]+) std ([\d.]+)\s*$",
            text, re.M):
        rows.append({"step": int(m.group(1)), "loss": float(m.group(2)),
                     "acc": float(m.group(3)), "pos": float(m.group(4)),
                     "worst_neg": float(m.group(5)), "margin": float(m.group(6)),
                     "wsi_wsi": float(m.group(7)), "rank": float(m.group(8)),
                     "std": float(m.group(9))})
    assert len(rows) >= 5, rows
    return {"rows": rows,
            "chance": float(re.search(r"cosine margin required[^\n]*?:\s*([\d.]+)", text).group(1))}


def collapse_endpoints() -> dict:
    """Arm A of the '16/16' diagnostic. ENDPOINT PAIRS ONLY - no array exists."""
    text = P.load_text("e0_run/collapse_diag.log")
    arm = text.split("--- A.")[1].split("--- B.")[0]
    def pair(label):
        m = re.search(rf"{label}\s+([\d.]+)\s*->\s*([\d.]+)", arm)
        assert m, f"{label!r} not found in collapse_diag.log arm A"
        return float(m.group(1)), float(m.group(2))
    hard = re.search(r"z_biology matrix rank (\d+) -> (\d+)\s+\(max (\d+)\)", arm)
    assert hard, "the 16/16 track is not in collapse_diag.log arm A"
    return {
        "infonce": pair("in-batch InfoNCE"),
        "acc": pair(r"retrieval acc@1"),
        "chance_acc": float(re.search(r"retrieval acc@1[^\n]*\(chance ([\d.]+)\)", arm).group(1)),
        "cross_pos": pair("cross pos cos"),
        "cross_neg": pair("cross neg cos"),
        "within": pair("WSI within-modality offdiag cos"),
        "hard": (int(hard.group(1)), int(hard.group(2))),
        "hard_max": int(hard.group(3)),
        "chance_infonce": float(re.search(r"in-batch InfoNCE[^\n]*\(chance ([\d.]+)", arm).group(1)),
    }


def collapse_rows(track: dict) -> list[dict]:
    """The three collapse instances of panel (a), each with its own statistic named."""
    rows = track["rows"]
    at50 = next(r for r in rows if r["step"] == 50)
    start = rows[0]

    d1b = P.repo_text(D1B)
    init = float(re.search(r"effective rank ([\d.]+) at step 0", d1b).group(1))
    # Columns: decorrelation | biology_full_consistency | step 50 | 100 | 150 | 200 | 250.
    # The header row also has seven cells and its step-150 cell reads "150", so rows
    # are kept only when the first cell is a bare number.
    tbl = [r for r in P.md_rows(d1b, n_cols=7)
           if re.fullmatch(r"[\d.]+", r[0].strip())]
    assert len(tbl) == 5, tbl
    step150 = [P.one_num(r[4]) for r in tbl if P.nums(r[4])]
    assert len(step150) == 4 and max(step150) < 10, step150
    rna_cos = float(re.search(r"RNA-view mutual cosine is ([\d.]+) with", d1b).group(1))

    mid = P.md_rows(P.repo_text(D1_MID), contains="programme_free", n_cols=6)
    end = P.md_rows(P.repo_text(D1A_END), contains="programme_free", n_cols=6)
    assert len(mid) == 1 and len(end) == 1, (mid, end)
    mid_epoch, mid_rank, mid_hard, mid_cos = (int(P.one_num(mid[0][1])), P.one_num(mid[0][2]),
                                              int(P.one_num(mid[0][3])), P.one_num(mid[0][4]))
    end_epoch, end_rank, end_hard, end_cos = (int(P.one_num(end[0][1])), P.one_num(end[0][2]),
                                              int(P.one_num(end[0][3])), P.one_num(end[0][4]))

    return [
        {"title": "clean in-batch InfoNCE,\nno queue",
         "stat": "diagnostic-script centred effective rank\n(the diagnostic's own definition; NOT\n"
                 "asserted equal to R1, R2 or R3 here)",
         "rank": f"{start['rank']:.2f} -> {at50['rank']:.2f} by step {at50['step']}",
         "evidence": [f"positive / worst-negative cosine {at50['pos']:.4f} / {at50['worst_neg']:.4f}",
                      f"minimum margin {start['margin']:+.3f} -> {at50['margin']:+.4f}"],
         "cos": at50["pos"]},
        {"title": "five programme_free configurations,\none common initialisation",
         "stat": P.STAT["R3_short"] + "\ncentred, fixed 256-patient held-out probe",
         "rank": f"{init:.2f} -> about 2 by step 150 "
                 f"({min(step150):.2f}-{max(step150):.2f} across five arms)",
         "evidence": [f"RNA-view mutual cosine {rna_cos:.4f} at step 200",
                      "no loss weighting prevents it, including both terms at zero"],
         "cos": rna_cos},
        {"title": "programme_free seed 42, real D1\ntraining, 282 held-out patients",
         "stat": P.STAT["R3_short"] + "\nheld-out probe on the live checkpoints",
         "rank": f"{mid_rank:.2f} at epoch {mid_epoch} and {end_rank:.2f} at epoch {end_epoch}",
         "evidence": [f"RNA-RNA mutual cosine {mid_cos:.3f} / {end_cos:.3f}",
                      f"hard rank {mid_hard} / {end_hard} of 256"],
         "cos": max(mid_cos, end_cos)},
    ]


def main() -> int:
    P.cli(__doc__)
    track = diag_d_track()
    ends = collapse_endpoints()
    rows = collapse_rows(track)
    d2 = P.load_json("e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json")

    fig = plt.figure(figsize=(7.6, 11.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.72, 1.15, 0.78], hspace=0.50, wspace=0.30,
                          left=0.075, right=0.975, top=0.960, bottom=0.225)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b1 = fig.add_subplot(gs[1, 0])
    ax_b2 = fig.add_subplot(gs[1, 1])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])

    # ------------------------------------------------------------------ (a)
    ax_a.set_title("(a)  The collapse regime, where rank works - and where we did NOT falsify it", pad=10)
    ax_a.axis("off")
    for k, r in enumerate(rows):
        y = 0.93 - 0.325 * k
        ax_a.text(0.0, y, r["title"], transform=ax_a.transAxes, ha="left", va="top",
                  fontsize=7.0, color=C.INK, fontweight="bold")
        ax_a.text(0.0, y - 0.115, r["stat"], transform=ax_a.transAxes, ha="left", va="top",
                  fontsize=5.9, color=C.MUTED, linespacing=1.5)
        ax_a.text(0.34, y, r["rank"], transform=ax_a.transAxes, ha="left", va="top",
                  fontsize=7.6, color=C.BLUE, fontweight="bold")
        ax_a.text(0.34, y - 0.06, "\n".join(r["evidence"]), transform=ax_a.transAxes,
                  ha="left", va="top", fontsize=6.2, color=C.INK, linespacing=1.6)
    ax_a.text(0.34, 0.99, "rank value and co-measured collapse evidence",
              transform=ax_a.transAxes, ha="left", va="bottom", fontsize=6.0, color=C.MUTED)

    # ----------------------------------------------------------------- (b1)
    ax_b1.set_title("(b)  The withdrawal: the same objective, two \"ranks\"", pad=10)
    steps = [r["step"] for r in track["rows"]]
    ranks = [r["rank"] for r in track["rows"]]
    hard0, hard1 = ends["hard"]
    hmax = ends["hard_max"]

    ax_b1.axhline(hmax, color=C.MUTED, linewidth=1.2, linestyle=(0, (5, 2)), zorder=1)
    ax_b1.annotate(f"structural maximum = batch size {hmax}", (steps[-1], hmax),
                   xytext=(-4, 5), textcoords="offset points", ha="right", va="bottom",
                   fontsize=6.4, color=C.MUTED, fontweight="bold")
    # ENDPOINT PAIR ONLY - two markers, no line between them that could be read
    # as an interpolation of a track that was never recorded.
    for x, v in ((steps[0], hard0), (steps[-1], hard1)):
        ax_b1.plot([x], [v], marker="X", color=C.VERMILLION, markersize=10,
                   markeredgecolor="white", markeredgewidth=1.1, zorder=4)
    ax_b1.annotate(f"hard `z_biology` matrix rank\n{hard0}/{hmax} -> {hard1}/{hmax}\n"
                   "ENDPOINT VALUES AS RECORDED;\nper-step array not retained",
                   (steps[-1], hard1), xytext=(-6, -8), textcoords="offset points",
                   ha="right", va="top", fontsize=6.2, color=C.VERMILLION,
                   fontweight="bold", linespacing=1.5)
    ax_b1.plot(steps, ranks, color=C.BLUE, marker="o", markersize=5.0,
               markeredgecolor="white", markeredgewidth=0.8, zorder=3,
               label="centred effective rank,\nSAME objective, PER-STEP")
    for x, v in zip(steps, ranks):
        if x in (0, 25, 50, 400):
            ax_b1.annotate(f"{v:.2f}", (x, v), xytext=(0, 7), textcoords="offset points",
                           ha="center", fontsize=6.4, color=C.BLUE, fontweight="bold")
    ax_b1.set_xlim(-18, 430)
    ax_b1.set_ylim(0, 19)
    ax_b1.set_xlabel("optimisation step", fontsize=6.8, labelpad=6)
    ax_b1.set_ylabel("rank\nTWO DIFFERENT STATISTICS, drawn together only because\n"
                     "the withdrawal is about mistaking one for the other", fontsize=6.2)
    P.grid(ax_b1)
    ax_b1.legend(loc="lower left", bbox_to_anchor=(0.10, 0.14), handlelength=1.8)

    # ----------------------------------------------------------------- (b2)
    ax_b2.set_title("      ... and the collapse evidence for the same arm", pad=10)
    ev = [("within-modality off-diagonal cosine", ends["within"], C.BLUE, "o"),
          ("cross-modal POSITIVE pair cosine", ends["cross_pos"], C.GREEN, "^"),
          ("cross-modal NEGATIVE pair cosine", ends["cross_neg"], C.VERMILLION, "s"),
          ("retrieval accuracy@1", ends["acc"], C.PURPLE, "D")]
    for k, (name, (a, b), colour, mk) in enumerate(ev):
        y = len(ev) - 1 - k
        ax_b2.annotate("", xy=(b, y), xytext=(a, y),
                       arrowprops=dict(arrowstyle="-|>", color=colour, linewidth=1.4,
                                       shrinkA=5, shrinkB=5))
        ax_b2.plot([a], [y], marker="o", color="white", markersize=7,
                   markeredgecolor=colour, markeredgewidth=1.6, zorder=3)
        ax_b2.plot([b], [y], marker=mk, color=colour, markersize=7,
                   markeredgecolor="white", markeredgewidth=0.9, zorder=3)
        ax_b2.annotate(f"{a:.4f}", (a, y), xytext=(0, -10), textcoords="offset points",
                       ha="center", va="top", fontsize=6.2, color=C.MUTED)
        ax_b2.annotate(f"{b:.4f}", (b, y), xytext=(0, 8), textcoords="offset points",
                       ha="center", fontsize=6.6, color=colour, fontweight="bold")
    ax_b2.axvline(ends["chance_acc"], color=C.INK, linewidth=1.0,
                  linestyle=(0, (3, 2)), zorder=1)
    ax_b2.annotate(f"chance = {ends['chance_acc']:.3f}\nretrieval ENDS BELOW IT",
                   (ends["chance_acc"], -0.62), xytext=(4, 0), textcoords="offset points",
                   ha="left", va="bottom", fontsize=6.3, color=C.INK, fontweight="bold")
    ax_b2.set_yticks(range(len(ev)))
    ax_b2.set_yticklabels([e[0] for e in reversed(ev)], fontsize=6.5)
    ax_b2.tick_params(axis="y", length=0)
    ax_b2.set_ylim(-0.85, len(ev) - 0.35)
    ax_b2.set_xlim(-0.06, 1.10)
    ax_b2.set_xlabel("value (all four are unitless and bounded)\n"
                     "open circle = before, filled = after; ENDPOINT PAIRS AS RECORDED,\n"
                     "nothing between them is interpolated", fontsize=6.4, labelpad=6)
    P.grid(ax_b2, axis="x")
    # The "positives and negatives are now indistinguishable" reading is left to
    # the caption: the two endpoint values are already printed on their own rows,
    # and a sentence here can only be placed on top of one of them.

    # ------------------------------------------------------------------ (c)
    ax_c.set_title("(c)  The boundary: where rank stops being informative", pad=10)
    nominal = d2["H44"]["n_test"] and 256
    pts = []
    for key, arm in (("H44", "Hallmark"), ("I44", "PBS")):
        e = d2[key]
        pts.append((arm, e["effective_rank_residualised"],
                    e["points"]["untrained40"]["top_cca"],
                    e["points"]["untrained40"]["null_mean"]))
    null = sum(p[3] for p in pts) / len(pts)
    ax_c.axvline(null, color=C.INK, linewidth=1.1, linestyle=(0, (3, 2)), zorder=1)
    ax_c.annotate(f"permutation null {null:.3f}\n(D2's own 200-draw null;\n"
                  "never to be drawn with F7's)",
                  (null, 5.6), xytext=(5, 0), textcoords="offset points", ha="left",
                  va="center", fontsize=6.3, color=C.INK)
    for k, (arm, rank, cca, _n) in enumerate(pts):
        ax_c.plot([cca], [rank], marker="o" if k == 0 else "s", color=P.ORDER[k],
                  markersize=17, markeredgecolor="white", markeredgewidth=1.6, zorder=4)
        ax_c.annotate(f"D2 s44 {arm}\nrank {rank:.2f} of a nominal {nominal}  "
                      f"({100 * rank / nominal:.1f}% of ambient)\nchannel {cca:.4f}",
                      (cca, rank), xytext=(0, 24 if k == 0 else -24),
                      textcoords="offset points",
                      ha="center", va="bottom" if k == 0 else "top", fontsize=6.4,
                      color=P.ORDER[k], fontweight="bold", linespacing=1.5)
    ax_c.set_xlim(0.10, 0.72)
    ax_c.set_ylim(3.9, 12.6)
    ax_c.set_xlabel("held-out top-CCA, 40 untrained targets\nconfound-residualised block",
                    fontsize=6.6, labelpad=6)
    ax_c.set_ylabel(P.axis_label("R1_short", "resid"), fontsize=6.3)
    P.grid(ax_c)
    ax_c.text(0.0, -0.40,
              "Two artifacts whose rank agrees to within the sampling noise of the estimator\n"
              f"({pts[0][1]:.2f} against {pts[1][1]:.2f}; F4(a) puts the gap at 1.4 sampling sd)\n"
              f"and whose channels differ by {abs(pts[0][2] - pts[1][2]):.4f}. At 3.6% of ambient "
              "dimensionality,\nrank is uninformative about the channel in BOTH directions.",
              transform=ax_c.transAxes, ha="left", va="top", fontsize=6.2, color=C.INK,
              linespacing=1.6)

    # ------------------------------------------------------------------ (d)
    ax_d.set_title("(d)  The cheaper alarm", pad=10)
    for k, r in enumerate(rows):
        y = len(rows) - 1 - k
        ax_d.barh(y, r["cos"], height=0.5, color=C.GREEN, edgecolor="white",
                  linewidth=1.5, zorder=2)
        ax_d.annotate(f"{r['cos']:.4f}", (r["cos"], y), xytext=(6, 0),
                      textcoords="offset points", ha="left", va="center",
                      fontsize=7.0, color=C.INK, fontweight="bold")
    ax_d.axvline(1.0, color=C.MUTED, linewidth=1.1, linestyle=(0, (5, 2)), zorder=1)
    ax_d.annotate("natural maximum = 1", (1.0, len(rows) - 0.42), xytext=(-4, 0),
                  textcoords="offset points", ha="right", va="bottom", fontsize=6.3,
                  color=C.MUTED, fontweight="bold")
    ax_d.set_yticks(range(len(rows)))
    ax_d.set_yticklabels([f"instance {len(rows) - k}" for k in range(len(rows))], fontsize=6.8)
    ax_d.tick_params(axis="y", length=0)
    ax_d.set_ylim(-0.6, len(rows) - 0.2)
    ax_d.set_xlim(0, 1.18)
    ax_d.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_d.set_xlabel("patient-to-patient mutual cosine, the SAME arms as panel (a)\n"
                    "one matrix product, no SVD", fontsize=6.6, labelpad=6)
    P.grid(ax_d, axis="x")
    lo = min(r["cos"] for r in rows)
    hi = max(r["cos"] for r in rows)
    ax_d.text(0.0, -0.40,
              f"{lo:.3f}-{hi:.4f} in every case, saturating at a natural maximum of 1.\n"
              "This is why the draft recommends the cosine over rank even for the\n"
              "one use of rank that survives.",
              transform=ax_d.transAxes, ha="left", va="top", fontsize=6.2, color=C.INK,
              linespacing=1.6)

    P.caption(fig,
        "F8. THE WITHDRAWAL, in full, because without (iv) and (v) this figure must not be published. "
        "(i) The 16/16 column is a HARD NUMERICAL RANK (torch.linalg.matrix_rank), not R1, R2 or R3, and a 16 x 256 float matrix has full row rank under "
        "essentially any perturbation; (ii) its maximum is 16 because the batch is 16; (iii) it is a TRAIN batch of 16, not held out; "
        "(iv) the CENTRED EFFECTIVE RANK OF THE SAME OBJECTIVE FALLS TO 1.00 by step 50, so this instance is evidence FOR the collapse-diagnostic use, "
        "not against it; (v) this project previously listed it among its two strongest instances and THAT DESCRIPTION IS WITHDRAWN, here and in P1. "
        "Cross-modal POSITIVE and NEGATIVE pairs end indistinguishable, 0.9959 against 0.9960 - the negatives marginally HIGHER - and in-batch InfoNCE ends "
        "at 2.7734 against its own chance of 2.7726. Panel (b)'s 16/16 track and all four collapse-evidence quantities are ENDPOINT PAIRS as "
        "recorded in ~/e0_run/collapse_diag.log, which probes "
        "before and after rather than on a schedule; no per-step array exists and none is interpolated. The source script for the original 16/16 "
        "measurement is a scratchpad file on the A100 and is not in this repository. The centred effective-rank track IS per-step and is drawn as one. "
        "Panel (a) reads \"we did not falsify the collapse-diagnostic use\", never \"we verified it\": we have not found a case of total collapse that "
        "effective rank missed, and the figure for that absence would imply a symmetry the data do not support. Panel (c)'s null is D2's own 200-draw "
        "null and must never share an axis with F7's within-cancer 0.145-0.147; they differ in n, in component count and in the permutation procedure.",
        y=0.006)

    P.save(fig, "F8_collapse_boundary_withdrawal")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
