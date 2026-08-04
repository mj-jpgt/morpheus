"""The retraining floor on the FIXED HELD-OUT PROBE — the block every §5 number sits on.

`p2_envelope_floors.py` gave every statistic and every exported view a floor of
its own and then had to name four blocks it could not reach. The first of those
is the fixed held-out probe, and it is not a minor omission: **every rank number
in draft §5 is measured there**, so until this module's inputs existed all of §5
was `unjudgeable` by the paper's own criterion — neither passing it nor failing
it. Its `ABSENT_BLOCKS` entry states the cost exactly:

    "A floor here needs the five repeats re-run with `d1_momentum_probe.py`
     attached — GPU, five runs."

That is the run this module reads. Draft §6.2 predeclared the same thing in the
same words ("five same-seed repeats of the `programme_free` / 500-step
configuration with `d1_momentum_probe.py` attached, read at a fixed step").

WHAT WAS RUN, AND WHAT IT IS A FLOOR OF. Ten runs, not five: **five identical
repeats of each of the two arms** §5's comparisons are drawn between —

    d1_momentum_probe.py {0.999,0.0} 0.04 500 4096 2e-4 42 <export_dir>

— same momentum, same decorrelation, same capacity, same learning rate, **same
seed 42**, same 500-step budget, same workspace, one after another on one A100.
The only difference between the five repeats of an arm is GPU non-determinism.
It is therefore a FLOOR TWICE OVER in exactly the sense §4.1 uses of the
exported-block floor: same-seed repeats exclude seed variation entirely, so the
spread measured here is strictly smaller than the spread a reader who varies the
seed would see. n = 5 per arm, no interval, one seed, one stack. A table of these
numbers is not a table of estimated distributions and must not be read as one.

BOTH ARMS, BECAUSE A COMPARISON HAS TWO SIDES. §4.1's floor was measured on one
arm (`programme_only`, the stable one) and applied to comparisons that span two.
That was defensible there and would not have been here: the low arm of every §5
comparison sits at or near the collapse floor, where a fold is the ratio of two
small numbers and can be enormous for reasons that have nothing to do with the
claim. `arm_floor` reports each arm separately and `floor` reports
`max(arm floors)` — the conservative one — with `floor_arm` naming which arm
carried it. Nothing here silently averages the two.

READ AT A FIXED STEP — AT EVERY STEP, ACTUALLY. The unjudgeable rows do not all
read at the same step: §4.4(3) reads at 200, §4.9a at 400, §5.4 rows 2 and 3 at
500, §5.2a at 200. A floor measured at step 500 does not license a verdict at
step 200 any more than a floor on one block licenses one on another, so the
floor is measured at **every step the harness reads** (0, 50, ... 500) and a
comparison is judged against the floor at its own step. §5.2's step-600 numbers
are outside this run's budget and are named as such rather than judged against
step 500.

NOTHING IS COMPUTED INLINE, AND NOTHING IS REDEFINED. The statistic table, the
fold, the bimodality descriptor, α-ReQ's index range and LiDAR's ridge are all
IMPORTED FROM `p2_envelope_floors.py`, so a floor measured here is measured under
exactly the definitions the exported-block floors were measured under and the two
tables are comparable by construction. Four inline-formula substitutions have been
caught in this paper and the tell was identical each time: an arithmetic
expression where an import belonged.

THE RAW BLOCK, AND WHY THERE IS NO RESIDUALISED ONE. `d1_momentum_probe.py`
reads its two statistics on the probe states as the model emits them; the rank
statistics centre internally and no confound residualisation is applied. That is
the block draft §5 quotes, so it is the block given a floor here. A residualised
probe block would need the probe rows' cancer and pooled-TSS labels, and — more
to the point — **no number in this paper is measured on it**, so a floor for it
would be a floor for nothing. `ABSENT` below says so rather than leaving it
ambiguous.

THE CONFIGURATION IS AN ARGUMENT, BECAUSE A FLOOR EXISTS ONLY WHERE IT WAS RUN.
The first run of this module left three audit rows unjudgeable and each named the
gap that stopped it: a reading at **step 600** past the 500 the repeats were run
to, an arm at **m = 0.99** that no repeat was run at, and an arm at **capacity 64**
that no repeat was run at. Closing them is the same measurement at other settings,
so `--steps`, `--arms`, `--arm-kind`, `--capacity` and `--momentum` are arguments
and they travel into the output's `config`, from which `p2_floor_audit`'s block
strings are built. Nothing is borrowed across any of them.

`--arm-kind capacity` exists for the third: §5.2 measurement 3's two arms differ
in QUEUE CAPACITY at fixed momentum rather than in momentum at fixed capacity, so
the arms are named `cap64`/`cap4096` and their repeats live in `cap<n>_rep<r>`.

AND A THREE-ARM RUN IS SCORED ONCE PER PAIR. §5.4's two step-600 rows are
m = 0.999 against m = 0 and m = 0.999 against m = 0.99. `combine()` takes the max
over the arms it is given, so scoring all three at once would judge the second row
against the collapsed arm's spread, which is not one of its two sides. Each pair
is therefore a separate invocation over the same repeat directories, and each
writes its own output file.

usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 python3 p2_probe_floors.py \
        --root /home/ubuntu/e0_run/d1_probefloor \
        --output ~/e0_run/d1_probefloor/out/P2_PROBE_FLOORS.json

    ... --root ~/e0_run/d1_probefloor600 --arms 0.999,0 --steps 600 \
        --read-step 600 --output .../P2_PROBE_FLOORS_S600_m0999_m0.json

    ... --root ~/e0_run/d1_capfloor --arm-kind capacity --arms 64,4096 \
        --momentum 0 --steps 200 --read-step 150 \
        --output .../P2_PROBE_FLOORS_CAP.json
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np

# Every definition is the exported-block module's, imported rather than restated:
# same statistics, same α-ReQ index range, same LiDAR ridge, same fold, same
# bimodality rule. If those change there, they change here.
from morpheus.v2.research.rebase.p2.p2_envelope_floors import (ALPHA_INDEX_RANGE, LIDAR_DELTA,
                                                               STATISTICS, _fold, _sha256)
from morpheus.v2.research.rebase.p2.p2_competing_metrics import lidar

#: The two views `d1_momentum_probe.py`'s `geometry()` reads. There is no
#: `full_biology` here: the probe forward pass takes `view="wsi"` and
#: `view="rna"` and nothing else, so unlike the exported artifact this block has
#: no third co-trained view to give a floor to.
VIEWS = ("wsi_biology", "rna_biology")

#: LiDAR's pseudo-view, named exactly as `p2_envelope_floors.py` names it so a
#: paired-block floor can never be matched to a single-view comparison.
LIDAR_VIEW = "wsi_rna_paired"

#: What a floor still cannot be measured for after this run, and why. Recorded
#: for the same reason `p2_envelope_floors.ABSENT_BLOCKS` is: a silence reads as
#: an oversight, a named absence reads as a result.
ABSENT = {
    "the residualised probe block": (
        "NOT A GAP. The probe states are read as the model emits them and no confound "
        "residualisation is applied at any point in `d1_momentum_probe.py`; every rank "
        "number draft §5 quotes is on the raw block. A residualised probe floor would be "
        "a floor for no number in this paper."),
    "step 700 and beyond, and every budget between the ones measured": (
        "This module is run once per step budget and a floor exists only at the budgets "
        "actually run — 500 (`--steps 500`) and 600 (`--steps 600`). A floor at one budget "
        "licenses no verdict at another any more than one block's floor licenses another's, "
        "and the reading step is written into the block string so that block-matching "
        "enforces it. Nothing in this paper reads past step 600."),
    "momentum values other than 0, 0.99 and 0.999, and capacities other than 64 and 4,096": (
        "The arms are an argument (`--arms`), and a floor exists only for the arms run. "
        "§5.2's turnover sweep also visits m = 0.9 and m = 0.95 and capacities 512, 2,048 "
        "and 8,192; no repeat has been run at any of them."),
    "learning rates other than 2e-4": (
        "Every repeat is at the project's training rate. §5.2a's own result is that the "
        "LEARNING RATE is the variable that moves rank most, so a floor measured at one "
        "rate is the last quantity that may be borrowed for another — which is why §5.2a's "
        "four contrasts stay unjudgeable."),
    "the `full_biology` view of the probe": (
        "`geometry()` performs two forward passes, `view='wsi'` and `view='rna'`. There is "
        "no third view on this block to measure."),
}

_STEP = re.compile(r"probe_step(\d+)\.npz$")

#: How an arm is named and where its repeats live, per kind of sweep. The arm
#: names an output KEY, and `p2_floor_audit.resolve` addresses this file with a
#: DOT-SEPARATED path, so a key containing a `.` would be silently split in two
#: and resolve to nothing. Momentum values are the one place in this project
#: where that is a live hazard, and capacities are integers so they are safe.
#:
#: `momentum` is the original behaviour and is byte-identical to it: the same
#: labels, the same directory names, so `P2_PROBE_FLOORS.json` regenerates
#: unchanged. `capacity` exists because §5.2 measurement 3's two arms differ in
#: QUEUE CAPACITY rather than in momentum, and a floor measured at one capacity
#: may not be borrowed for another.
ARM_KINDS = {
    "momentum": {
        "label": lambda arm: "m0" if float(arm) == 0 else "m" + str(arm).replace(".", ""),
        "directory": lambda arm, rep: f"m{arm}_rep{rep}",
        "varies": "momentum",
    },
    "capacity": {
        "label": lambda arm: f"cap{int(float(arm))}",
        "directory": lambda arm, rep: f"cap{arm}_rep{rep}",
        "varies": "queue capacity",
    },
}


def arm_label(arm: str, kind: str = "momentum") -> str:
    """`0.999` -> `m0999`, `64` -> `cap64`. See `ARM_KINDS`."""
    return ARM_KINDS[kind]["label"](arm)


# --------------------------------------------------------------------------
def score_state(path: str) -> dict:
    """Every statistic on every probe view, for one repeat at one step."""
    raw = np.load(path, allow_pickle=False)
    blocks = {v: np.asarray(raw[v], dtype=np.float64) for v in VIEWS}
    record = {
        "path": str(Path(path).resolve()),
        "sha256": _sha256(path),
        "n_probe": int(blocks["wsi_biology"].shape[0]),
        "views": {v: {name: float(fn(blocks[v])) for name, fn in STATISTICS.items()}
                  for v in VIEWS},
    }
    record["views"][LIDAR_VIEW] = {
        "LiDAR": float(lidar([blocks["wsi_biology"], blocks["rna_biology"]],
                             delta=LIDAR_DELTA)["lidar"])}
    return record


def collect_arm(reps: dict[str, dict[int, dict]], steps: list[int]) -> dict:
    """{step: {view: {statistic: fold entry}}} across one arm's repeats."""
    out: dict[str, dict] = {}
    for step in steps:
        out[str(step)] = {}
        for view in list(VIEWS) + [LIDAR_VIEW]:
            names = sorted({s for r in reps.values() for s in r[step]["views"][view]})
            out[str(step)][view] = {
                stat: _fold({rep: r[step]["views"][view][stat] for rep, r in reps.items()})
                for stat in names}
    return out


def combine(arms: dict[str, dict], steps: list[int]) -> dict:
    """The conservative floor: the larger of the two arms', with the arm named.

    A §5 comparison is a ratio between a run of one arm and a run of the other,
    so the noise it has to clear is the noise of whichever arm is noisier. Taking
    the max is a choice and it is recorded as one; `arm_floor` keeps both so a
    reader can see what was taken the max of.
    """
    out: dict[str, dict] = {}
    for step in steps:
        key = str(step)
        out[key] = {}
        for view in list(VIEWS) + [LIDAR_VIEW]:
            names = sorted({s for arm in arms.values() for s in arm[key][view]})
            out[key][view] = {}
            for stat in names:
                per_arm = {a: arms[a][key][view][stat] for a in arms}
                folds = {a: e["fold"] for a, e in per_arm.items()}
                usable = {a: f for a, f in folds.items() if f is not None}
                best = max(usable, key=usable.get) if usable else None
                out[key][view][stat] = {
                    "floor": (usable[best] if best else None),
                    "floor_arm": best,
                    # The two numbers the floor IS, exposed at the top level so
                    # `p2_floor_audit.resolve_ratio` can re-derive it from its own
                    # two sources rather than trust the ratio recorded beside them.
                    "floor_min": (per_arm[best]["min"] if best else None),
                    "floor_max": (per_arm[best]["max"] if best else None),
                    "why": None if best else "no arm yielded a defined fold",
                    "arm_floor": {a: {"fold": e["fold"], "min": e["min"], "max": e["max"],
                                      "values": e["values"], "shape": e["shape"]}
                                  for a, e in per_arm.items()},
                }
    return out


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="/home/ubuntu/e0_run/d1_probefloor")
    ap.add_argument("--arms", default="0.999,0.0")
    ap.add_argument("--arm-kind", default="momentum", choices=sorted(ARM_KINDS),
                    help="what the arms vary — `momentum` (the original) or `capacity`")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--read-step", type=int, default=500,
                    help="the step whose floor is printed in the summary table")
    # The three settings held fixed ACROSS the arms. They are arguments rather
    # than constants because a floor exists only at the configuration it was run
    # at, and the configuration has to travel with the numbers into the output —
    # `p2_floor_audit`'s block strings are built from it.
    ap.add_argument("--steps", type=int, default=500, help="the run's step budget")
    ap.add_argument("--capacity", default="4096", help="queue capacity, when it is not the arm")
    ap.add_argument("--momentum", default=None, help="momentum, when it is not the arm")
    ap.add_argument("--lr", default="2e-4")
    ap.add_argument("--decorrelation", default="0.04")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output", required=True)
    args = ap.parse_args(argv)

    kind = ARM_KINDS[args.arm_kind]
    if args.arm_kind == "capacity" and args.momentum is None:
        raise SystemExit("--arm-kind capacity needs --momentum: it is held fixed across the arms")
    root = Path(args.root)
    arms_raw: dict[str, dict[str, dict[int, dict]]] = {}
    arm_argv: dict[str, str] = {}
    steps_seen: set[int] | None = None
    for arm in args.arms.split(","):
        label = arm_label(arm, args.arm_kind)
        mom = arm if args.arm_kind == "momentum" else args.momentum
        cap = arm if args.arm_kind == "capacity" else args.capacity
        arm_argv[label] = (f"d1_momentum_probe.py {mom} {args.decorrelation} {args.steps} "
                           f"{cap} {args.lr} {args.seed} <export_dir>")
        reps: dict[str, dict[int, dict]] = {}
        for r in range(1, args.reps + 1):
            d = root / kind["directory"](arm, r)
            files = sorted(glob.glob(str(d / "probe_step*.npz")))
            if not files:
                raise SystemExit(f"no probe states under {d}")
            per_step = {}
            for f in files:
                m = _STEP.search(Path(f).name)
                per_step[int(m.group(1))] = score_state(f)
            reps[f"rep{r}"] = per_step
            got = set(per_step)
            steps_seen = got if steps_seen is None else (steps_seen & got)
            wsi = per_step[max(per_step)]["views"]["wsi_biology"]
            print(f"  {label:<7} rep{r}  step{max(per_step):<4} "
                  f"R1={wsi['R1']:8.3f}  R3={wsi['R3']:8.3f}  "
                  f"RankMe={wsi['RankMe']:7.3f}  hard={wsi['hard_rank']:6.0f}", flush=True)
        arms_raw[label] = reps

    steps = sorted(steps_seen or [])
    arms = {arm: collect_arm(reps, steps) for arm, reps in arms_raw.items()}
    floors = combine(arms, steps)

    out = {
        "_": (f"The retraining floor on the FIXED HELD-OUT PROBE — the block every rank "
              f"number in draft §5 is measured on, and the block `p2_envelope_floors.py` "
              f"named as unreachable from the exported repeats. {args.reps} identical "
              f"same-seed repeats of EACH of the {len(arms)} arms compared, which here vary "
              f"in {kind['varies']}. A FLOOR TWICE OVER (same-seed repeats exclude seed "
              f"variation entirely), n = {args.reps} per arm, no interval, one seed, one "
              f"stack. Not a distribution."),
        "config": {
            "harness": "v2/research/rebase/d1_momentum_probe.py",
            "argv": arm_argv,
            "objective_profile": "programme_free",
            "arms": list(arms), "repeats_per_arm": args.reps,
            "arm_kind": args.arm_kind, "arms_vary": kind["varies"],
            "seed": args.seed, "steps": args.steps,
            "capacity": (None if args.arm_kind == "capacity" else int(args.capacity)),
            "momentum": (None if args.arm_kind == "momentum" else float(args.momentum)),
            "lr": float(args.lr), "decorrelation": float(args.decorrelation),
            "block": "fixed held-out probe, raw (256 held-out patients, as the model emits them)",
            "steps_measured": steps,
            "alpha_index_range": list(ALPHA_INDEX_RANGE), "lidar_delta": LIDAR_DELTA,
            "statistic_sites": {
                "STATISTICS/_fold/_shape/alpha range/LiDAR delta":
                    "v2/research/rebase/p2/p2_envelope_floors.py — imported, not restated",
                "R1/R2/R3": "v2/calibra/spectral.py RANK_VARIANTS",
                "PR/PR_rownorm/RankMe/stable_rank/alpha_req/LiDAR":
                    "v2/research/rebase/p2/p2_competing_metrics.py",
                "hard_rank": "numpy.linalg.matrix_rank — NOT an effective rank"},
        },
        "shape_rule": ("outlier = the run whose removal minimises the remaining four's fold; "
                       "concordant = that fold <= 1.05; bimodal = concordant and "
                       "full fold / rest fold >= 2.0 — imported from p2_envelope_floors.py"),
        "floor_rule": (f"floor = max over the {len(arms)} arms of that arm's own max/min "
                       f"across its {args.reps} same-seed repeats, at the SAME step and the "
                       f"SAME statistic; `floor_arm` names which arm carried it and "
                       f"`arm_floor` keeps both. A comparison is judged against the floor "
                       f"built from ITS OWN two arms, so a three-arm run is scored once per "
                       f"pair rather than once overall."),
        "absent": ABSENT,
        "reps": {arm: {rep: {str(s): rec for s, rec in per_step.items()}
                       for rep, per_step in reps.items()} for arm, reps in arms_raw.items()},
        "arms": arms,
        "floors": floors,
    }
    path = Path(args.output).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")

    read = str(args.read_step)
    print(f"\nfloor on the fixed held-out probe at step {read} "
          f"(max/min over five same-seed repeats, per arm)")
    print(f"{'view':>16} {'statistic':>18} {'arm':>7} {'min':>10} {'max':>10} {'fold':>8}  shape")
    for view, stats in floors[read].items():
        for stat, entry in stats.items():
            for arm, a in entry["arm_floor"].items():
                fold = "n/a" if a["fold"] is None else f"{a['fold']:.3f}x"
                shape = a["shape"]
                tag = ("undefined" if not shape.get("defined") else
                       (f"BIMODAL outlier={shape['outlier']}" if shape["bimodal"]
                        else f"not bimodal (outlier={shape['outlier']}, "
                             f"rest={shape['rest_fold']:.3f}x, sep={shape['separation']:.2f})"))
                lo = "n/a" if a["min"] is None else f"{a['min']:10.4f}"
                hi = "n/a" if a["max"] is None else f"{a['max']:10.4f}"
                print(f"{view:>16} {stat:>18} {arm:>7} {lo:>10} {hi:>10} {fold:>8}  {tag}")
            combined = "n/a" if entry["floor"] is None else "%.3fx" % entry["floor"]
            print(f"{'':>16} {stat:>18} {'FLOOR':>7} {'':>10} {'':>10} "
                  f"{combined:>8}  carried by {entry['floor_arm']}")

    print("\nper-repeat, per-step canonical R1 and R3 on `wsi_biology` — never a mean")
    for arm, reps in arms_raw.items():
        for rep, per_step in sorted(reps.items()):
            row = " ".join(f"{s}:{per_step[s]['views']['wsi_biology']['R1']:.2f}/"
                           f"{per_step[s]['views']['wsi_biology']['R3']:.2f}" for s in steps)
            print(f"  {arm:<8} {rep}  {row}")
    print("\nwrote", path)
    return out


if __name__ == "__main__":
    main()
