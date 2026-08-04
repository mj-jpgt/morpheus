"""F1 - The envelope, and the seven arm differences inside it.  (P2 draft 4.1)

Claim: six of the seven between-arm rank differences this project has ever
measured are smaller than the spread of the same statistic when one
configuration is retrained with the same seed.

THE ENVELOPE IS n = 1. It is one retraining pair on one configuration on one
stack. It is therefore drawn as a single observation - two markers and the
interval between them, in grey, with the count printed - and never as a shaded
confidence region, because there is no interval to shade. Panel (d) is the
controlled repeat that would give it one; it has not reported.

Sources
  (a) data/e0_run/d2_v3/RECOVERED_SEED42_READOUT.json   re-export of the surviving checkpoint
      data/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json    the retrained arm
      v2/research/rebase/nature/D2_RESULT.md 4           the recorded original (rank not recorded)
  (b) data/ws_rank/RANK_RECOMPUTE.json                  D2 + D1-B, R1, residualised
      data/ws_rank/RANK_RECOMPUTE_P1B.json              Phase 1b, R1, RAW
  (c) v2/research/rebase/nature/D2_RESULT.md 2          stratified 40-target readout
  (d) data/extracted/F1_RETRAINING_REPEAT.json          state of ~/e0_run/d1_envelope
"""
from __future__ import annotations

import p2fig as P
from p2fig import C, np, plt


def envelope() -> dict:
    """The 2.69x retraining envelope, as three observations of one configuration."""
    rec = P.load_json("e0_run/d2_v3/RECOVERED_SEED42_READOUT.json")["H42_recovered"]
    ret = P.load_json("e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json")["H42"]
    # The recorded original: its channel survives in D2_RESULT.md 4; its rank was
    # never recorded, which is why that row has no rank marker in panel (a).
    md = P.repo_text("v2/research/rebase/nature/D2_RESULT.md")
    row = P.md_rows(md, contains="recorded in the lost run", n_cols=3)
    assert len(row) == 1, row
    recorded_channel = P.one_num(row[0][1])
    return {
        "recorded_channel": recorded_channel,
        "reexport_rank": rec["effective_rank_residualised"],
        "reexport_channel": rec["points"]["unrestricted"]["top_cca"],
        "retrained_rank": ret["effective_rank_residualised"],
        "retrained_channel": ret["points"]["unrestricted"]["top_cca"],
    }


def arm_ratios() -> list[tuple[str, float, str]]:
    """The seven between-arm rank ratios, as (label, fold, statistic/block tag)."""
    rr = P.load_json("ws_rank/RANK_RECOMPUTE.json")
    d2 = rr["instances"]["instance6_D2"]
    d1 = rr["instances"]["instance5_D1B"]
    out = []
    for seed in (42, 43, 44):
        h, i = d2[f"H{seed}"]["residualised"]["R1"], d2[f"I{seed}"]["residualised"]["R1"]
        hi = "H" if h > i else "I"
        out.append((f"D2 s{seed}  ({hi} higher)", max(h, i) / min(h, i), "resid"))
    for seed in (42, 43, 44):
        p = d1[f"d1_p_seed{seed}"]["residualised"]["R1"]
        f = d1[f"d1_f_seed{seed}"]["residualised"]["R1"]
        hi = "programme_only" if p > f else "programme_free"
        out.append((f"D1-B s{seed}  ({hi} higher)", max(p, f) / min(p, f), "resid"))
    # Phase 1b is the one RAW-block point and its arms were never verified
    # matched (draft 3.4). It is tagged so panel (b) can mark it.
    p1b = P.load_json("ws_rank/RANK_RECOMPUTE_P1B.json")["v21_release_20260720_retry3_resume_safe"]
    full = p1b["diagnostic_full_seed42"]["wsi_biology"]["raw"]["R1"]
    prog = p1b["diagnostic_programme_only_seed42"]["wsi_biology"]["raw"]["R1"]
    out.append(("Phase 1b s42  (full higher)", max(full, prog) / min(full, prog), "raw"))
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
    args = P.cli(__doc__)
    env = envelope()
    fold = env["retrained_rank"] / env["reexport_rank"]
    ratios = arm_ratios()
    pairs = d2_paired_channel()
    repeat = P.load_json("extracted/F1_RETRAINING_REPEAT.json")

    fig = plt.figure(figsize=(7.5, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.72, wspace=0.50,
                          left=0.155, right=0.985, top=0.945, bottom=0.235)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # ---------------------------------------------------------------- (a)
    ax_a.set_title("(a)  The envelope, established: $n$ = 1", pad=10)
    rows_a = [
        ("recorded original\n(lost run)", None, env["recorded_channel"], "x"),
        ("re-export of the\nsurviving checkpoint", env["reexport_rank"], env["reexport_channel"], "o"),
        ("retrained,\nsame seed 42", env["retrained_rank"], env["retrained_channel"], "D"),
    ]
    # The envelope drawn as what it is: the interval between TWO observations.
    ax_a.plot([env["reexport_rank"], env["retrained_rank"]], [1, 0], color=C.ENVELOPE,
              linewidth=7, solid_capstyle="round", zorder=1, alpha=0.6)
    for y, (lab, r, ch, mk) in zip((2, 1, 0), rows_a):
        if r is None:
            ax_a.text(6.6, y, "rank not recorded", ha="left", va="center",
                      fontsize=6.2, color=C.MUTED, style="italic")
        else:
            ax_a.plot([r], [y], marker=mk, color=C.BLUE, markersize=8,
                      markeredgecolor="white", markeredgewidth=1.0, zorder=3)
            ax_a.annotate(f"{r:.4f}", (r, y), xytext=(0, 9), textcoords="offset points",
                          ha="center", fontsize=6.8, color=C.BLUE, fontweight="bold")
        ax_a.annotate(f"channel {ch:.5f}".rstrip("0").rstrip("."), (26.6, y),
                      xytext=(0, -11), textcoords="offset points",
                      ha="right", va="center", fontsize=6.3, color=C.VERMILLION)
    ax_a.set_xlim(6.0, 27.0)
    ax_a.set_ylim(-0.65, 2.65)
    ax_a.set_yticks([2, 1, 0])
    ax_a.set_yticklabels([r[0] for r in rows_a], fontsize=6.4)
    ax_a.tick_params(axis="y", length=0)
    ax_a.spines["left"].set_visible(False)
    ax_a.set_xlabel(P.axis_label("R1_short", "resid"), fontsize=6.5, labelpad=6)
    P.grid(ax_a, axis="x")
    ax_a.text(0.0, -0.40,
              f"retraining: rank $\\times${fold:.2f},  channel "
              f"{env['retrained_channel'] - env['recorded_channel']:+.3f}\n"
              f"re-export deterministic to 5 s.f. "
              f"({env['reexport_channel']:.5f} against {env['recorded_channel']:.4f} recorded)",
              transform=ax_a.transAxes, ha="left", va="top", fontsize=6.4, color=C.INK)

    # ---------------------------------------------------------------- (b)
    ax_b.set_title("(b)  The seven between-arm differences, against it", pad=10)
    hi = fold
    ax_b.axvspan(1.0, hi, color=C.ENVELOPE, alpha=0.32, zorder=0, linewidth=0)
    ax_b.axvline(hi, color=C.MUTED, linewidth=1.0, linestyle=(0, (4, 2)), zorder=1)
    ax_b.axvline(1.0, color=C.RULE, linewidth=0.8, zorder=1)
    inside = 0
    for k, (lab, val, tag) in enumerate(ratios):
        within = val <= hi
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
    ax_b.set_xticks([1.0, 1.5, 2.0, round(hi, 2), 4.0])
    ax_b.set_xticklabels(["1.0", "1.5", "2.0", f"{hi:.2f}", "4.0"])
    ax_b.minorticks_off()
    ax_b.set_ylim(-1.35, len(ratios) - 0.35)
    ax_b.set_xlabel("between-arm rank ratio, higher / lower  (log scale)\n"
                    "R1 (canonical effective rank, centred);  circle = residualised block,\n"
                    "triangle = Phase 1b, RAW block, arms never verified matched", fontsize=6.4,
                    labelpad=6)
    P.grid(ax_b, axis="x")
    ax_b.text(1.02, -1.28,
              f"shaded: the $n$ = 1 retraining envelope, 1.00-{hi:.2f}$\\times$.\n"
              "One pair, one configuration. Not an interval, so no band edge is a bound.",
              fontsize=6.1, color=C.MUTED, va="bottom", ha="left")
    ax_b.text(0.03, 0.975, f"{inside} of {len(ratios)} inside the envelope",
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
    ax_c.set_xlim(-0.205, 0.145)
    ax_c.set_xlabel("paired channel difference, PBS $-$ Hallmark\n"
                    "held-out top-CCA, 40 targets neither arm trained on\n"
                    "confound-residualised block;  CI$_{95}$, 2,000 resamples", fontsize=6.4,
                    labelpad=6)
    P.grid(ax_c, axis="x")
    ax_c.legend(loc="upper right", bbox_to_anchor=(1.0, 1.02), ncol=1, handlelength=1.8)
    same_sign = sum(p["delta"] < 0 for p in pairs)
    excl = sum(p["patient"][1] < 0 and p["cancer"][1] < 0 for p in pairs)
    ax_c.text(0.99, 0.50,
              f"same sign {same_sign}/3\nboth CIs exclude zero {excl}/3\n\n"
              "the channel is a PAIRED within-run\ndifference; rank is not, and\n"
              "there is no paired form of\n\"this run has higher rank\"",
              transform=ax_c.transAxes, ha="right", va="center", fontsize=6.4,
              color=C.INK, linespacing=1.55)

    # ---------------------------------------------------------------- (d)
    n_rep = len(repeat["reps"])
    ep = sorted({v["epochs_logged"] for v in repeat["reps"].values()})
    progress = (f"{n_rep}/{n_rep} logging, "
                f"{ep[0] if len(ep) == 1 else f'{ep[0]}-{ep[-1]}'} epochs so far, "
                f"0/{n_rep} exported")
    if repeat.get("complete"):
        raise SystemExit(
            "F1(d): the controlled retraining repeat has REPORTED. This panel must now be "
            "redrawn as a real per-repeat rank/channel display, per "
            "NOTEBOOK_ENTRIES/PREDECLARED_retraining_envelope_20260804T0330Z.md. Refusing to "
            "emit the pending placeholder over live data.")
    P.pending_panel(
        ax_d,
        title="(d)  The controlled repeat that would make (a) an interval",
        marker="[RETRAINING ENVELOPE PENDING]",
        blocked_on=f"{n_rep} same-seed retraining repeats\n({progress})",
        arrives_at=f"{repeat['root']}\n{repeat['expected']}",
        predeclaration=repeat["predeclaration"],
        strict=args.strict)

    P.caption(fig,
        "F1. The envelope of panel (a) is ONE retraining pair on ONE configuration on ONE stack (D2 arm H, seed 42, identical arguments). "
        "It cannot separate rank-specific variance from stack non-determinism, architecture or schedule, and it is drawn as a single observation, "
        "not as an estimated distribution: no confidence band is shown because none has been measured. Statistic R1 (Roy & Vetterli order-1, centred), "
        "block confound-residualised for every D2 and D1-B point and RAW for the Phase 1b point - and the Phase 1b arms were never verified matched (draft 3.4). "
        "Panel (b)'s x-axis is a ratio of two R1 values within one comparison; no two rank statistics share an axis anywhere in this figure. "
        "The claim does not rest on this figure alone: F2 reaches the same conclusion from 8 within-arm degrees of freedom without using the "
        f"{fold:.2f}x number. Panels (b) and (c) sit side by side because the channel is quoted as a paired within-run difference and rank is not - "
        "there is no paired form of \"this run has higher rank\".",
        y=0.006)

    P.save(fig, "F1_envelope")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
