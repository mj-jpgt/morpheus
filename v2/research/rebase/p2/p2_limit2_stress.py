"""§5.4 limit 2 under stress: more repeats, every statistic, and a momentum grid.

WHAT THIS IS FOR. `p2_probe_floors.py` measured a floor for the comparison this
project actually ships — `m = 0.999` over `m = 0.99` on the fixed held-out probe
at step 600 — from five same-seed repeats per arm, and the row cleared it by
5.6%. The same ten runs then undercut the pass under one statistic and supported
it under another, and the disagreement was left standing:

    Test A (the audit rule, `p2_floor_audit.check`): the PUBLISHED single-draw
        ratio must exceed the floor.  R3: 1.262x against 1.195x -- clears.
    Test B (worst case over the repeats): min(high arm) / max(low arm) must
        exceed the floor.  R3: 1.138x against 1.195x -- FAILS.
                           R1: 1.453x against 1.155x -- clears.

This module computes both tests, for every statistic, at any number of repeats,
so that the disagreement can be pushed on rather than restated. Predeclared in
`NOTEBOOK_ENTRIES/PREDECLARED_probe_floor_n10_and_momentum_grid_20260805T0200Z.md`.

THE ARITHMETIC THAT MAKES HALF OF IT UNINFORMATIVE, SAID UP FRONT. The floor is
`max/min` over an arm's repeats, so it is NON-DECREASING in the repeat count. The
Test B separation is `min(high)/max(low)`, so it is NON-INCREASING in it. Adding
repeats can therefore only ever make Test B harder. A statistic that fails Test B
at n = 5 and fails it at n = 10 has learned nothing; a statistic that PASSES it
at n = 10 has survived a test that had every opportunity to break it. Test A is
the one that can genuinely flip, because the published ratio is fixed while the
floor grows. `verdict()` labels each accordingly rather than leaving a reader to
work it out.

NOTHING IS COMPUTED INLINE. The statistic table, the fold and the sha are
`p2_envelope_floors`'s; the rank variants are
`effective_rank(..., variant=RANK_VARIANTS[...])`; the extra alternatives are
`p2_competing_metrics`'s. Four inline-formula substitutions have been caught in
this paper and each one was an arithmetic expression where an import belonged.

THE EXTRA VARIANTS, AND WHY THEY ARE NOT NEW STATISTICS. `STATISTICS` scores
`R1`, `R2` and `R3`. `RANK_VARIANTS` carries three more -- `R1_uncentred`,
`R1_rownorm`, `R2_uncentred` -- which are the same one implementation at other
settings of its two documented preprocessing switches. They are added here so
that "R1 passes and R3 fails" can be decomposed along the two axes that actually
separate them (Hill ORDER, and row NORMALISATION) instead of being left as a
clash between two labels. Nothing is invented: every one of the six is a key of
`RANK_VARIANTS` already.

AND DUPLICATES ARE DETECTED, NOT ASSUMED. `z_biology` is L2-normalised at the
model's output, so on this block `normalise_rows` is a no-op and R2 = R3 and
R1 = R1_rownorm to float precision. That was established for R2/R3 on 2026-08-04
and is CHECKED here numerically for every pair, because a statistic that
duplicates another must not be counted as independent corroboration of it.

usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 python3 p2_limit2_stress.py \
        --root ~/e0_run/d1_probefloor600 --high 0.999 --low 0.99 \
        --reps 10 --step 600 --output .../P2_LIMIT2_STRESS_N10.json

    ... --reps 5 --rep-offset 5    # repeats 6-10 alone, an independent n = 5
    ... --grid 0,0.98,0.99,0.995,0.999 --output .../P2_MOMENTUM_GRID.json
"""
from __future__ import annotations

import argparse
import json
from math import comb
from pathlib import Path

import numpy as np

from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2.p2_envelope_floors import STATISTICS, _fold, _sha256
from morpheus.v2.research.rebase.p2.p2_probe_floors import VIEWS, arm_label

#: The three `RANK_VARIANTS` keys `STATISTICS` does not already carry. Same
#: function, other settings of its two documented switches; see the module
#: docstring for why they are here and why they are not new statistics.
EXTRA_VARIANTS = ("R1_uncentred", "R1_rownorm", "R2_uncentred")

#: Everything scored, in one table. `STATISTICS` first so its keys keep their
#: meaning, then the extra variants.
STATS = dict(STATISTICS)
STATS.update({name: (lambda x, _v=RANK_VARIANTS[name]: effective_rank(x, variant=_v))
              for name in EXTRA_VARIANTS})

#: The published single-run values §5.4 limit 2's Test A is computed from, and the
#: ONLY statistic for which that test exists on this row. `long_m0.999.log` and
#: `long_m0.99.log` were written by the version of `d1_momentum_probe.py` that
#: printed FOUR columns -- `eff-rank` (R3), `feat-std`, `rna-rna`, `contrastive`.
#: There is no `CANONICAL` column in either file and the states were never
#: exported, so no R1 ratio and no PR / RankMe / stable-rank / alpha-ReQ ratio can
#: be recovered for this row at all. Test A is therefore an R3-only test, which is
#: itself part of the finding: the statistic under which the row FAILS Test B is
#: the only statistic under which it can be given Test A.
PUBLISHED_TEST_A = {
    "R3": {
        "high": 7.42, "low": 5.88, "ratio": 7.42 / 5.88,
        "src": ("e0_run/d1_diag/long_m0.999.log and long_m0.99.log, step 600, "
                "column `eff-rank` -- the two runs `floor_audit.json`'s "
                "`5.4-m0999-over-m099` resolves its `a` and `b` from"),
    },
}

#: Which way "better" points, per statistic. `+1` means larger is better quality.
#: The one exception is the same one `p2_selection_rule.METRICS` records and for
#: the same reason: alpha-ReQ is scored as PROXIMITY TO 1 (`|alpha - 1|`, smaller
#: is better), which is THIS PROJECT'S operationalisation of the paper's
#: "Goldilocks zone" and not the authors' rule. Without this table alpha-ReQ's
#: perfectly ordered arms would be printed as INVERTED -- a sign convention
#: reported as a contradiction, which is the kind of error that makes a table look
#: like a disagreement between statistics when it is a bug.
#:
#: `alpha_req` itself (the raw exponent, not the deviation) has NO agreed
#: direction on this project -- `p2_selection_rule` scores only the deviation --
#: so it is left at 0 and its ordering is printed as a raw fact with no verdict
#: attached rather than being given a direction it does not have.
DIRECTION = {"R1": +1, "R2": +1, "R3": +1, "R1_uncentred": +1, "R1_rownorm": +1,
             "R2_uncentred": +1, "PR": +1, "PR_rownorm": +1, "RankMe": +1,
             "stable_rank": +1, "hard_rank": +1, "alpha_req_abs_dev": -1, "alpha_req": 0}

#: What no amount of repeating fixes, recorded for the same reason
#: `p2_probe_floors.ABSENT` is: a silence reads as an oversight.
ABSENT = {
    "Test A under any statistic but R3": (
        "The two runs the row is a claim ABOUT were logged before the harness printed "
        "a canonical column and before `export_dir` existed. Their states are gone. "
        "R1's support for this row is therefore Test B only, and R3's Test A pass "
        "cannot be corroborated by a second statistic on the same two runs."),
    "a seed-varied floor": (
        "Every repeat here is seed 42. The floor excludes seed variation entirely and "
        "is smaller than one a reader varying the seed would measure -- a floor twice "
        "over, exactly as `p2_probe_floors.py` says of its own."),
    "an interval on any number here": (
        "n = 5 or n = 10 same-seed repeats on one stack. These are ranges, not "
        "estimated distributions, and no confidence statement is made from them."),
}


# --------------------------------------------------------------------------
def score(path: str) -> dict:
    """Every statistic in `STATS` on every probe view, for one state file."""
    raw = np.load(path, allow_pickle=False)
    blocks = {v: np.asarray(raw[v], dtype=np.float64) for v in VIEWS}
    return {"path": str(Path(path).resolve()), "sha256": _sha256(path),
            "views": {v: {n: float(fn(blocks[v])) for n, fn in STATS.items()} for v in VIEWS}}


#: A pair whose values agree to better than this RELATIVE difference on every
#: repeat of every arm is the same statistic on this block. 1e-6 rather than
#: something tighter because the known case is not exact: R2 and R3 coincide
#: because `z_biology` is L2-normalised at output, and the 2026-08-04 measurement
#: of them on one state read 6.9711779953 against 6.9711779832 -- a relative
#: 1.7e-9 that a 1e-9 test misses. The number that decides this is REPORTED, not
#: only the verdict, so the threshold cannot hide a near-miss.
DUPLICATE_REL = 1e-6

#: A looser second tier, because one pair on this block sits between the two and
#: would otherwise be counted as independent corroboration when it is not:
#: `RankMe` and `R1_uncentred` are both exp(Shannon entropy) of the UNCENTRED
#: singular-value distribution and differ only by RankMe's ``eps``, which the
#: RankMe paper adds OUTSIDE the normalisation. They read 3.013207 against
#: 3.012465 on one state -- a relative 2.5e-4, far above `DUPLICATE_REL` and far
#: below anything that makes them two statistics. Reported in its own tier with
#: the number attached, rather than by loosening the exact test until it catches.
NEAR_DUPLICATE_REL = 1e-3


def duplicates(per_arm_per_rep: dict[str, dict[str, dict[str, float]]]) -> dict:
    """Which statistics are the SAME statistic on this block, checked not assumed.

    Two statistics agreeing to `DUPLICATE_REL` on EVERY repeat of EVERY arm are
    not two pieces of evidence, and counting them as two would inflate any
    "N statistics agree" statement. The worst (largest) relative difference seen
    over all repeats is reported for every duplicate pair, so a reader can see how
    close "the same" was.
    """
    rows = [v for arm in per_arm_per_rep.values() for v in arm.values()]
    names = [n for n in STATS if all(np.isfinite(r.get(n, np.nan)) for r in rows)]
    pairs, near, dup_of = {}, {}, {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            worst = max(abs(r[a] - r[b]) / abs(r[b]) for r in rows if r[b] != 0)
            if worst <= DUPLICATE_REL:
                pairs[f"{b} == {a}"] = worst
                dup_of.setdefault(b, a)
            elif worst <= NEAR_DUPLICATE_REL:
                near[f"{b} ~= {a}"] = worst
                dup_of.setdefault(b, a)
    return {"duplicate_of": dup_of, "worst_relative_difference": pairs,
            "near_duplicates": near,
            "thresholds": {"duplicate": DUPLICATE_REL, "near_duplicate": NEAR_DUPLICATE_REL},
            "why": ("`z_biology` is L2-normalised at the model's output, so on this block "
                    "`normalise_rows` is a no-op. A duplicate is NOT independent "
                    "corroboration and is excluded from any count of agreeing statistics; "
                    "a near-duplicate is excluded too and its relative difference is "
                    "printed so the exclusion can be argued with.")}


def arm_values(root: Path, arm: str, step: int, reps: int, offset: int, view: str) -> dict:
    """{rep: {statistic: value}} for one arm at one step, plus the state records."""
    per_rep, records = {}, {}
    for r in range(offset + 1, offset + reps + 1):
        f = root / f"m{arm}_rep{r}" / f"probe_step{step}.npz"
        if not f.is_file():
            raise SystemExit(f"missing {f}")
        rec = score(str(f))
        records[f"rep{r}"] = rec
        per_rep[f"rep{r}"] = rec["views"][view]
    return {"per_rep": per_rep, "records": records}


def ordering(high: list[float], low: list[float], direction: int = +1) -> dict:
    """Test C: are the two arms COMPLETELY SEPARATED, and how surprising is that?

    Test A and Test B are both statements about a RATIO -- how big the gap is
    against how big the within-arm spread is -- and Test B's two sides both move
    against the pass as repeats are added. This asks the other question, which no
    number in this paper has yet asked of this row: **is every repeat of the high
    arm above every repeat of the low arm?** That is a statement about ORDER, not
    about magnitude, and unlike Test B it gets STRONGER with more repeats, so the
    two tests fail and succeed for different reasons and neither substitutes for
    the other.

    `p` is the EXACT one-sided permutation probability of complete separation
    under exchangeability of the `n_high + n_low` repeats -- one arrangement out
    of `C(n_high + n_low, n_high)` -- which is the smallest p any rank test can
    return at these sample sizes. At n = 5 per arm it is 1/252; at n = 10 per arm
    it is 1/184,756.

    **WHAT IT IS NOT.** The repeats differ only in GPU non-determinism: same seed,
    same workspace, same stack, one card. Exchangeability holds within THAT noise
    source and nothing else, so this p is a statement about run-to-run
    reproducibility and is NOT a p-value for the momentum effect, which would need
    the seed varied. §4.2 measures the seed as the dominant term. Recorded here in
    the same words it will be reported in.
    """
    if not high or not low or not all(np.isfinite(high + low)):
        return {"separated": None, "why": "a non-finite value"}
    if direction == 0:
        return {"separated": None, "direction": 0,
                "high_range": [min(high), max(high)], "low_range": [min(low), max(low)],
                "why": ("no agreed direction for this statistic on this project -- "
                        "`p2_selection_rule.METRICS` scores only the |alpha - 1| deviation, "
                        "so the raw exponent is reported and not judged")}
    # A statistic whose "better" is SMALLER is separated when the high arm sits
    # entirely BELOW the low arm. Comparing the raw values without this would
    # print alpha-ReQ's perfectly ordered arms as a contradiction.
    if min(high + low) == max(high + low):
        return {"separated": None, "direction": direction, "constant": True,
                "value": high[0], "n_high": len(high), "n_low": len(low),
                "why": ("the statistic is PINNED at one value in every repeat of both arms, "
                        "so it neither orders nor overlaps them -- it cannot rule. Printing "
                        "this as an overlap would read as a statistic that disagreed, which "
                        "is the opposite of what it did.")}
    h, l = ([v * direction for v in high], [v * direction for v in low])
    sep = min(h) > max(l)
    inverted = max(h) < min(l)
    n = len(high) + len(low)
    return {
        "separated": bool(sep), "inverted": bool(inverted), "direction": direction,
        "high_range": [min(high), max(high)], "low_range": [min(low), max(low)],
        "gap": min(h) - max(l),
        "exact_one_sided_permutation_p": (1.0 / comb(n, len(high))) if sep else None,
        "n_high": len(high), "n_low": len(low),
        "informative": ("EITHER outcome is informative, and this test gets HARDER TO PASS "
                        "and MORE surprising when passed as repeats are added -- the "
                        "opposite of Test B."),
        "scope": ("exchangeability holds over GPU non-determinism only: same seed, same "
                  "workspace, one card. NOT a p-value for the momentum effect."),
    }


def verdict(sep: float | None, floor: float | None, monotone_in_n: bool) -> dict:
    """`sep > floor`, with what the comparison is WORTH attached to it."""
    if sep is None or floor is None:
        return {"passes": None, "why": "a non-finite value on one side"}
    # `hard_rank` is 256 in every repeat of every arm, so separation and floor are
    # both exactly 1 and `sep > floor` is False. Reporting that as a FAILURE would
    # be wrong in the direction that flatters the paper's thesis: the statistic did
    # not reject the comparison, it is incapable of ruling on it at all. §4.9's
    # "16/16" instance is the same statistic doing the same thing at batch size.
    if sep == floor == 1.0:
        return {"passes": None, "separation": sep, "floor": floor, "margin": 1.0,
                "why": ("NO DISCRIMINATION -- the statistic is pinned at the same value in "
                        "every repeat of both arms, so its floor and its separation are "
                        "both exactly 1. It cannot rule either way and is not counted as "
                        "a failure."),
                "informative": "neither outcome is available from this statistic"}
    return {
        "passes": bool(sep > floor), "separation": sep, "floor": floor,
        "margin": sep / floor,
        "informative": (
            "PASSING is informative, FAILING is not: both sides of this test move "
            "against the pass as repeats are added, so a failure at larger n is "
            "arithmetically expected and a survival is not."
            if monotone_in_n else
            "EITHER outcome is informative: the ratio is fixed (it is a property of "
            "two specific published runs) while the floor grows with repeats, so a "
            "crossing is a real change of verdict."),
    }


# --------------------------------------------------------------------------
def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="/home/ubuntu/e0_run/d1_probefloor600")
    ap.add_argument("--high", default="0.999", help="the arm the row claims is higher")
    ap.add_argument("--low", default="0.99", help="the arm it is claimed to be higher than")
    ap.add_argument("--grid", default="", help="score these arms side by side instead of testing a pair")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--rep-offset", type=int, default=0,
                    help="skip this many repeats first, so repeats 6-10 can be scored alone")
    ap.add_argument("--step", type=int, default=600)
    ap.add_argument("--view", default="wsi_biology")
    ap.add_argument("--output", required=True)
    args = ap.parse_args(argv)

    root = Path(args.root).expanduser()
    arms = args.grid.split(",") if args.grid else [args.high, args.low]
    data = {a: arm_values(root, a, args.step, args.reps, args.rep_offset, args.view) for a in arms}

    # Per arm, per statistic: the fold across its own repeats -- its arm floor.
    per_arm = {}
    for a in arms:
        names = sorted({s for r in data[a]["per_rep"].values() for s in r})
        per_arm[arm_label(a)] = {
            s: _fold({rep: v[s] for rep, v in data[a]["per_rep"].items()}) for s in names}

    out = {
        "_": (f"§5.4 limit 2 under stress. {args.reps} same-seed repeats per arm "
              f"(offset {args.rep_offset}) at step {args.step} on `{args.view}`, arms "
              f"{arms}. Ranges, not distributions: seed 42 throughout, one stack."),
        "config": {"root": str(root), "arms": arms, "reps": args.reps,
                   "rep_offset": args.rep_offset, "step": args.step, "view": args.view,
                   "statistic_sites": {
                       "STATISTICS/_fold/_sha256":
                           "v2/research/rebase/p2/p2_envelope_floors.py -- imported, not restated",
                       "R1/R2/R3/R1_uncentred/R1_rownorm/R2_uncentred":
                           "v2/calibra/spectral.py RANK_VARIANTS -- one implementation",
                       "PR/PR_rownorm/RankMe/stable_rank/alpha_req/LiDAR":
                           "v2/research/rebase/p2/p2_competing_metrics.py",
                       "hard_rank": "numpy.linalg.matrix_rank -- NOT an effective rank"}},
        "absent": ABSENT,
        "arm_floor": per_arm,
        "duplicates_on_this_block": duplicates({a: data[a]["per_rep"] for a in arms}),
        "reps": {arm_label(a): {rep: rec for rep, rec in data[a]["records"].items()}
                 for a in arms},
    }

    if args.grid:
        # STEP or SMOOTH? Every ADJACENT pair on the grid gets the same three
        # readings the shipped pair gets, so "m = 0.999 beats m = 0.99" can be
        # placed against its neighbours instead of being read as a two-point
        # contrast. A step function and a smooth trend support very different
        # claims and the draft currently asserts monotonicity while its own §5.2
        # table has m = 0.9 BELOW m = 0 at this step.
        steps = {}
        for a, b in zip(arms[:-1], arms[1:]):
            key = f"{arm_label(b)}_over_{arm_label(a)}"
            names = sorted(per_arm[arm_label(a)])
            steps[key] = {}
            for s in names:
                d = DIRECTION[s]
                fa, fb = per_arm[arm_label(a)][s]["fold"], per_arm[arm_label(b)][s]["fold"]
                floor = max([f for f in (fa, fb) if f is not None], default=None)
                if d > 0:
                    himin, lomax = per_arm[arm_label(b)][s]["min"], per_arm[arm_label(a)][s]["max"]
                elif d < 0:
                    himin, lomax = per_arm[arm_label(a)][s]["min"], per_arm[arm_label(b)][s]["max"]
                else:
                    himin = lomax = None
                sep = (himin / lomax) if (himin is not None and lomax and lomax > 0) else None
                steps[key][s] = {
                    "floor": floor, "direction": d,
                    "test_B_worst_case_over_repeats": verdict(sep, floor, monotone_in_n=True),
                    "test_C_complete_separation": ordering(
                        [v[s] for v in data[b]["per_rep"].values()],
                        [v[s] for v in data[a]["per_rep"].values()], d)}
        out["adjacent_steps"] = steps
        out["reading_rule"] = (
            "STEP: the grid is flat within its own floor across the lower arms and rises "
            "only at the top one -- the two-point comparison is then well posed and the "
            "closeness of the neighbour is not evidence against it. SMOOTH: every adjacent "
            "pair separates by roughly the same amount, so the shipped comparison is a "
            "coarse read of a trend and the honest claim is about the trend. NEITHER: "
            "non-monotone in m, which undercuts the momentum story more broadly than "
            "§5.4 limit 2 does. Predeclared at "
            "NOTEBOOK_ENTRIES/PREDECLARED_probe_floor_n10_and_momentum_grid_20260805T0200Z.md §6.")
    else:
        hi, lo = arm_label(args.high), arm_label(args.low)
        names = sorted(per_arm[hi])
        tests = {}
        for s in names:
            folds = [per_arm[x][s]["fold"] for x in (hi, lo)]
            floor = max([f for f in folds if f is not None], default=None)
            floor_arm = None
            if floor is not None:
                floor_arm = hi if per_arm[hi][s]["fold"] == floor else lo
            # The worst-case ratio has to be taken in the direction "better"
            # actually points, or a smaller-is-better statistic is graded upside
            # down: for alpha-ReQ's |alpha - 1| the adverse case is the high arm's
            # LARGEST deviation against the low arm's SMALLEST.
            d = DIRECTION[s]
            if d > 0:
                himin, lomax = per_arm[hi][s]["min"], per_arm[lo][s]["max"]
            elif d < 0:
                himin, lomax = per_arm[lo][s]["min"], per_arm[hi][s]["max"]
            else:
                himin = lomax = None
            sep = (himin / lomax) if (himin is not None and lomax and lomax > 0) else None
            entry = {"floor": floor, "floor_arm": floor_arm,
                     "arm_fold": {x: per_arm[x][s]["fold"] for x in (hi, lo)},
                     "high_min": himin, "low_max": lomax,
                     "test_B_worst_case_over_repeats": verdict(sep, floor, monotone_in_n=True),
                     "direction": d,
                     "test_C_complete_separation": ordering(
                         [v[s] for v in data[args.high]["per_rep"].values()],
                         [v[s] for v in data[args.low]["per_rep"].values()], d)}
            pub = PUBLISHED_TEST_A.get(s)
            entry["test_A_published_single_draw"] = (
                dict(verdict(pub["ratio"], floor, monotone_in_n=False), **{
                    "published_high": pub["high"], "published_low": pub["low"],
                    "src": pub["src"]}) if pub else
                {"passes": None, "why": ABSENT["Test A under any statistic but R3"]})
            tests[s] = entry
        out["tests"] = tests

    path = Path(args.output).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")

    # ---- printed table ----
    dup = out["duplicates_on_this_block"]["duplicate_of"]
    print(f"\nn = {args.reps} per arm (offset {args.rep_offset}), step {args.step}, "
          f"`{args.view}`, arms {arms}")
    for tier, key in (("DUPLICATE", "worst_relative_difference"), ("NEAR-DUPLICATE", "near_duplicates")):
        for pair, worst in out["duplicates_on_this_block"][key].items():
            print(f"  {tier} ON THIS BLOCK (not independent evidence): {pair}  "
                  f"worst relative difference {worst:.2e}")
    if args.grid:
        print(f"\n{'statistic':>18} " + " ".join(f"{arm_label(a):>22}" for a in arms))
        for s in sorted(per_arm[arm_label(arms[0])]):
            cells = []
            for a in arms:
                e = per_arm[arm_label(a)][s]
                cells.append("n/a".rjust(22) if e["min"] is None else
                             f"{e['min']:8.3f}-{e['max']:<8.3f}{e['fold']:5.2f}x")
            print(f"{s:>18} " + " ".join(cells))
        print("\nadjacent steps on the grid -- is the shipped pair a STEP or part of a trend?")
        for key, per_stat in out["adjacent_steps"].items():
            row = []
            for s in ("R1", "R3", "RankMe"):
                b, c = (per_stat[s]["test_B_worst_case_over_repeats"],
                        per_stat[s]["test_C_complete_separation"])
                mark = ("SEP" if c.get("separated") else
                        "INV" if c.get("inverted") else "ovl")
                sep = "n/a" if b.get("separation") is None else f"{b['separation']:.3f}x"
                v = {True: "PASS", False: "FAIL", None: "-"}[b.get("passes")]
                row.append(f"{s}: {sep} vs {per_stat[s]['floor']:.3f}x {v:<4} {mark}")
            print(f"  {key:>22}  " + " | ".join(row))
    else:
        print(f"\n{'statistic':>18} {'floor':>8} {'by':>7} {'TestB sep':>10} {'B':>6} "
              f"{'TestA ratio':>12} {'A':>6} {'C order':>9} {'p':>10}  dup")
        for s, e in out["tests"].items():
            b, a = e["test_B_worst_case_over_repeats"], e["test_A_published_single_draw"]
            c = e["test_C_complete_separation"]
            fl = "n/a" if e["floor"] is None else f"{e['floor']:.3f}x"
            bs = "n/a" if b.get("separation") is None else f"{b['separation']:.3f}x"
            bv = {True: "PASS", False: "FAIL", None: "none"}[b.get("passes")]
            as_ = "n/a" if a.get("separation") is None else f"{a['separation']:.3f}x"
            av = {True: "PASS", False: "FAIL", None: "-"}[a.get("passes")]
            cv = ("SEP" if c.get("separated") else
                  "INVERTED" if c.get("inverted") else
                  "constant" if c.get("constant") else
                  "no-dir" if c.get("direction") == 0 else "overlap")
            cp = "-" if c.get("exact_one_sided_permutation_p") is None else \
                 f"{c['exact_one_sided_permutation_p']:.2e}"
            print(f"{s:>18} {fl:>8} {str(e['floor_arm']):>7} {bs:>10} {bv:>6} "
                  f"{as_:>12} {av:>6} {cv:>9} {cp:>10}  {dup.get(s, '')}")
    print("\nper repeat, never a mean")
    for a in arms:
        for rep in sorted(data[a]["per_rep"], key=lambda r: int(r[3:])):
            v = data[a]["per_rep"][rep]
            print(f"  {arm_label(a):<8} {rep:<7} R1={v['R1']:8.3f}  R3={v['R3']:8.3f}  "
                  f"RankMe={v['RankMe']:7.3f}  PR={v['PR']:8.3f}  stable={v['stable_rank']:6.3f}")
    print("\nwrote", path)
    return out


if __name__ == "__main__":
    main()
