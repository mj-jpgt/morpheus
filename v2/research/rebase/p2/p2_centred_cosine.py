"""Does centring explain §5.2a's dissociation? The mutual cosine, centred and uncentred, on the same states.

THE QUESTION. §5.2a records that at `lr = 1e-3` the three momentum arms hold
centred effective rank flat — 1.06 / 1.05 / 1.05, a spread of 1.01× — while the
co-measured RNA-view mutual cosine on the same three runs falls 0.9946 → 0.9257
→ 0.5207, a factor of 1.91. §4.10 and §6.2 both refuse to adjudicate between two
accounts:

  (A) the dissociation is REAL — rank is insensitive to a genuine difference in
      how degenerate these representations are, at the reading where §4.10 says
      the collapse diagnostic is reliable; or
  (B) the difference is a MEAN OFFSET — our rank is column-centred and the mutual
      cosine is not, so a difference confined to the mean-offset direction
      produces exactly this pattern, and rank is right to ignore it.

The measurement that separates them is the cosine recomputed on the CENTRED
representation, on the same states. §6.2 names it and prices it: it needs those
runs' activations, which the six vendored `lr_L*` logs do not carry.

WHAT THIS MODULE READS. The probe states written by `d1_momentum_probe.py`'s
purely additive `export_dir` argument, for

    d1_momentum_probe.py {0, 0.9, 0.999} 0.04 200 4096 1e-3 42 <export_dir>

— the three `lr = 1e-3` arms of §5.2a (L3, L1, L5) at their own decorrelation,
capacity, seed and step budget — with THREE same-seed repeats of each arm, so
that "moves" and "flat" are verdicts against a measured same-seed spread rather
than against an eyeball. n = 3 per arm, one seed, one stack: a floor twice over
in §4.1's sense, and not a distribution.

THE THIRD ACCOUNT, WHICH THE DRAFT DOES NOT NAME. `geometry()` takes two forward
passes: the rank columns are computed on ``view="wsi"`` and the `rna-rna` cosine
on ``view="rna"``. §5.2a and §4.10 place the 1.01× rank spread and the 1.91×
cosine movement side by side as two instruments on one block; they are two
instruments on two VIEWS. So

  (C) the difference is real and lives in the RNA view, which the quoted rank
      number does not look at.

is live, and is separated by scoring R1 and R3 on the `rna_biology` states as
well as the `wsi_biology` ones. Every table below therefore carries both views.

NOTHING IS COMPUTED INLINE. R1/R2/R3 come from `v2/calibra/spectral.py` through
`p2_envelope_floors.STATISTICS`, and the fold and the bimodality descriptor are
that module's too — imported, not restated. Four inline-formula substitutions
have been caught in this paper and the tell was identical each time: an
arithmetic expression where an import belonged. The one quantity that has no
existing implementation is the mutual cosine itself; it is defined ONCE here,
with a `centre` flag, so that the centred and uncentred forms cannot drift apart.

THE BIT-LEVEL GUARD, AND IT IS A STOPPING CONDITION. The uncentred cosine
recomputed from the saved states must reproduce the harness's own printed
`rna-rna` column to four decimal places at every step of every run. If it does
not, the saved states are not the states the printed column was read from and
nothing computed here may be used. `--verify-logs` performs the check and
`check_against_logs` raises rather than warns.

usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 python3 p2_centred_cosine.py \
        --root /home/ubuntu/e0_run/d1_lrcentre \
        --logs /home/ubuntu/e0_run/d1_lrcentre \
        --output ~/e0_run/d1_lrcentre/out/P2_CENTRED_COSINE.json
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np

# Imported, not restated: the same statistic table, the same fold and the same
# bimodality rule the exported-block and probe-block floors were measured under.
from morpheus.v2.research.rebase.p2.p2_envelope_floors import STATISTICS, _fold, _sha256

#: The two views `d1_momentum_probe.py`'s `geometry()` reads. The rank columns it
#: PRINTS are the `wsi_biology` ones; the `rna-rna` cosine it prints is the
#: `rna_biology` one. Both are scored here under every statistic, because the
#: draft pairs a WSI-view rank with an RNA-view cosine and that pairing is itself
#: one of the three accounts under test.
VIEWS = ("wsi_biology", "rna_biology")

#: The view the mutual cosine is read on — the one `geometry()`'s `rna-rna`
#: column is computed from, named rather than positional.
COSINE_VIEW = "rna_biology"

#: The harness column the uncentred cosine must reproduce, and the precision at
#: which the log prints it. `d1_momentum_probe.py` formats it `%9.4f`, so four
#: decimal places is the most the log can carry and the tolerance is set to half
#: a unit in the last place plus a margin for the float32→float64 path.
LOG_COSINE_COLUMN = "rna-rna"
LOG_COSINE_TOLERANCE = 1e-4

_STEP = re.compile(r"probe_step(\d+)\.npz$")
_ARM = re.compile(r"^m([0-9.]+)_rep([0-9]+)$")


# --------------------------------------------------------------------------
# The one quantity with no existing implementation
# --------------------------------------------------------------------------
def mutual_cosine(x, *, centre: bool) -> float:
    """Mean off-diagonal cosine between the rows of ``x``.

    ``centre=False`` is exactly what `d1_momentum_probe.geometry()` prints in its
    `rna-rna` column: rows L2-normalised, Gram matrix formed, diagonal dropped,
    mean taken. ``centre=True`` subtracts the COLUMN mean first and then
    normalises — the same preprocessing `spectral.CANONICAL` applies before its
    SVD, so that the cosine and the canonical rank are finally reading the same
    matrix.

    One function with a flag rather than two functions, because the whole point
    of the measurement is that the two differ in exactly one preprocessing step
    and nothing else.
    """
    value = np.asarray(x, dtype=np.float64)
    if value.ndim != 2:
        raise ValueError(f"mutual_cosine expects a 2-D matrix, got shape {value.shape}")
    if centre:
        value = value - value.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(value, axis=-1, keepdims=True)
    unit = value / np.maximum(norms, 1e-300)
    gram = unit @ unit.T
    off = ~np.eye(len(gram), dtype=bool)
    return float(gram[off].mean())


def mean_offset_ratio(x) -> float:
    """‖column mean‖ / RMS row norm of the column-centred matrix.

    The mean-offset account's direct observable, and the secondary predeclared in
    `NOTEBOOK_ENTRIES/PREDECLARED_centred_cosine_20260804T1700Z.md` §5: account
    (B) predicts this falls sharply across the three arms while the centred
    statistics do not move. It EXPLAINS the primary verdict; it cannot overturn
    it, and no verdict in this module is taken on it.
    """
    value = np.asarray(x, dtype=np.float64)
    centre = value.mean(axis=0, keepdims=True)
    deviation = value - centre
    rms = float(np.sqrt((deviation ** 2).sum(axis=1).mean()))
    return float(np.linalg.norm(centre) / rms) if rms > 0 else float("nan")


#: Label -> callable on one (n, d) matrix, for the two quantities this module adds.
#: Kept in a table of the same shape as `p2_envelope_floors.STATISTICS` so that
#: `_fold` treats them identically and a reader cannot mistake one for a rank.
COSINE_STATISTICS = {
    "mutual_cosine_uncentred": lambda x: mutual_cosine(x, centre=False),
    "mutual_cosine_centred": lambda x: mutual_cosine(x, centre=True),
    "mean_offset_ratio": mean_offset_ratio,
}


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------
def score_state(path: str) -> dict:
    """Every rank statistic on both views, plus the two cosines, for one state."""
    raw = np.load(path, allow_pickle=False)
    blocks = {v: np.asarray(raw[v], dtype=np.float64) for v in VIEWS}
    record = {
        "path": str(Path(path).resolve()),
        "sha256": _sha256(path),
        "n_probe": int(blocks[COSINE_VIEW].shape[0]),
        "views": {v: {name: float(fn(blocks[v])) for name, fn in STATISTICS.items()}
                  for v in VIEWS},
    }
    for name, fn in COSINE_STATISTICS.items():
        record["views"][COSINE_VIEW][name] = float(fn(blocks[COSINE_VIEW]))
    return record


def probe_log_cosines(path: str) -> dict[int, float]:
    """{step: `rna-rna`} out of a `d1_momentum_probe.py` log.

    The header is matched rather than assumed, exactly as `p2_floor_audit`
    matches it, so a log written by a different script fails here instead of
    silently yielding whichever column happens to sit at that index.
    """
    header: tuple[str, ...] | None = None
    out: dict[int, float] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        cells = line.split()
        if cells[:1] == ["step"]:
            header = tuple(cells[1:])
            if LOG_COSINE_COLUMN not in header:
                raise ValueError(f"{path}: no {LOG_COSINE_COLUMN!r} column in header {header!r}")
            continue
        if header is None or not cells or not cells[0].isdigit() or len(cells) - 1 != len(header):
            continue
        out[int(cells[0])] = float(cells[1:][header.index(LOG_COSINE_COLUMN)])
    if not out:
        raise ValueError(f"{path}: no probe rows parsed")
    return out


def check_against_logs(reps: dict[str, dict[str, dict[int, dict]]],
                       logs: dict[str, dict[str, str]]) -> list[dict]:
    """Recomputed uncentred cosine against the harness's own printed column.

    RAISES on any disagreement. The states either are the states the printed
    column was read off or they are not, and if they are not then no number in
    this module means what it says.
    """
    checked: list[dict] = []
    problems: list[str] = []
    for arm, per_rep in reps.items():
        for rep, per_step in per_rep.items():
            printed = probe_log_cosines(logs[arm][rep])
            for step, record in sorted(per_step.items()):
                got = record["views"][COSINE_VIEW]["mutual_cosine_uncentred"]
                if step not in printed:
                    problems.append(f"{arm}/{rep}: step {step} has a state but no log row")
                    continue
                delta = abs(got - printed[step])
                checked.append({"arm": arm, "rep": rep, "step": step,
                                "recomputed": got, "printed": printed[step], "delta": delta})
                if delta > LOG_COSINE_TOLERANCE:
                    problems.append(
                        f"{arm}/{rep} step {step}: recomputed uncentred cosine {got:.6f} but the "
                        f"log prints {printed[step]:.4f} — STOP, and report both")
    if problems:
        raise SystemExit("STATE/LOG DISAGREEMENT — nothing here may be used:\n  "
                         + "\n  ".join(problems))
    return checked


# --------------------------------------------------------------------------
# Folding
# --------------------------------------------------------------------------
def within_arm(reps: dict[str, dict[int, dict]], steps: list[int]) -> dict:
    """{step: {view: {statistic: fold across this arm's repeats}}} — the floor."""
    out: dict[str, dict] = {}
    for step in steps:
        out[str(step)] = {}
        for view in VIEWS:
            names = sorted({s for r in reps.values() for s in r[step]["views"][view]})
            out[str(step)][view] = {
                stat: _fold({rep: r[step]["views"][view][stat] for rep, r in reps.items()})
                for stat in names}
    return out


def across_arms(arms: dict[str, dict[str, dict[int, dict]]], steps: list[int],
                repeat: str) -> dict:
    """The between-arm spread, one repeat index at a time.

    Read per repeat and never as a mean: §5.2a is one seed per cell, and the
    question is what the difference BETWEEN the arms is, so each repeat index
    yields its own three-arm comparison and all of them are printed.
    """
    out: dict[str, dict] = {}
    for step in steps:
        out[str(step)] = {}
        for view in VIEWS:
            names = sorted({s for arm in arms.values()
                            for s in arm[repeat][step]["views"][view]})
            out[str(step)][view] = {}
            for stat in names:
                values = {arm: arms[arm][repeat][step]["views"][view][stat] for arm in arms}
                entry = _fold(values)
                entry["spread"] = float(max(values.values()) - min(values.values()))
                out[str(step)][view][stat] = entry
    return out


def _spread(entry: dict) -> float:
    return float(entry["max"] - entry["min"]) if entry.get("max") is not None else float("nan")


def verdict(within: dict[str, dict], across: dict[str, dict], step: int) -> dict:
    """The reading predeclared at `PREDECLARED_centred_cosine_20260804T1700Z.md` §3.

    Taken on the ABSOLUTE spread across the three arms (max − min), because a
    cosine is bounded in [−1, 1] and may sit near zero, where a fold is not a
    meaningful scale. Computed here rather than asserted in prose, and nothing
    about the direction was predeclared because nothing about it was ours to
    choose.
    """
    key = str(step)
    out: dict = {"step": step, "rule": "PREDECLARED_centred_cosine_20260804T1700Z.md §3",
                 "moves_threshold": 0.20, "flat_threshold": 0.10}
    for stat in ("mutual_cosine_centred", "mutual_cosine_uncentred"):
        # The floor: the largest within-arm spread over the three repeats, taken
        # over the arms, exactly as `p2_probe_floors.combine` takes the max over
        # the two arms rather than averaging them.
        floor_by_arm = {arm: _spread(w[key][COSINE_VIEW][stat])
                        for arm, w in within.items()}
        floor_arm = max(floor_by_arm, key=floor_by_arm.get)
        # The across-arm spread: the largest over the repeat indices, so that the
        # verdict is not carried by the most flattering draw.
        by_repeat = {rep: a[key][COSINE_VIEW][stat]["spread"] for rep, a in across.items()}
        out[stat] = {
            "across_arm_spread_per_repeat": by_repeat,
            "across_arm_spread_max": max(by_repeat.values()),
            "across_arm_spread_min": min(by_repeat.values()),
            "within_arm_spread": floor_by_arm,
            "floor": floor_by_arm[floor_arm],
            "floor_arm": floor_arm,
            # Each branch is tested against the repeat draw that is LEAST
            # favourable to it: "moves" must hold for the smallest of the three
            # across-arm spreads and clear the floor there, "flat" must hold for
            # the largest and stay inside the floor there. The predeclaration
            # fixes the thresholds; taking the worst draw rather than the mean is
            # this paper's standing rule and is applied to both branches alike.
            "clears_own_floor": min(by_repeat.values()) > floor_by_arm[floor_arm],
            "inside_own_floor": max(by_repeat.values()) <= floor_by_arm[floor_arm],
        }
    centred, uncentred = out["mutual_cosine_centred"], out["mutual_cosine_uncentred"]
    c_spread = centred["across_arm_spread_min"]
    u_spread = uncentred["across_arm_spread_min"]
    if c_spread >= 0.20 and centred["clears_own_floor"]:
        reading = ("A — the dissociation is REAL: the centred cosine moves too, so rank is "
                   "missing a change a centred co-measure can see")
    elif (centred["across_arm_spread_max"] <= 0.10 and centred["inside_own_floor"]
          and u_spread >= 0.20 and uncentred["clears_own_floor"]):
        reading = ("B — the difference is a MEAN OFFSET: only the uncentred cosine moves, rank "
                   "is right to ignore it, and the dissociation DISSOLVES")
    else:
        reading = ("neither — report the magnitudes and do not adjudicate "
                   "(predeclared third branch)")
    out["reading"] = reading
    return out


# --------------------------------------------------------------------------
def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="/home/ubuntu/e0_run/d1_lrcentre",
                    help="directory holding `m<arm>_rep<n>/probe_step*.npz`")
    ap.add_argument("--logs", default=None,
                    help="directory holding `lrc_m<arm>_rep<n>.log`; defaults to --root")
    ap.add_argument("--arms", default="0,0.9,0.999")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--read-step", type=int, default=200,
                    help="the step §5.2a reads at, and the step the verdict is taken at")
    ap.add_argument("--output", required=True)
    args = ap.parse_args(argv)

    root, logdir = Path(args.root), Path(args.logs or args.root)
    arms: dict[str, dict[str, dict[int, dict]]] = {}
    logs: dict[str, dict[str, str]] = {}
    steps_seen: set[int] | None = None
    for arm in args.arms.split(","):
        label = "m" + str(arm).replace(".", "")
        arms[label], logs[label] = {}, {}
        for r in range(1, args.reps + 1):
            d = root / f"m{arm}_rep{r}"
            files = sorted(glob.glob(str(d / "probe_step*.npz")))
            if not files:
                raise SystemExit(f"no probe states under {d}")
            per_step = {}
            for f in files:
                m = _STEP.search(Path(f).name)
                per_step[int(m.group(1))] = score_state(f)
            arms[label][f"rep{r}"] = per_step
            logs[label][f"rep{r}"] = str(logdir / f"lrc_m{arm}_rep{r}.log")
            got = set(per_step)
            steps_seen = got if steps_seen is None else (steps_seen & got)
            last = per_step[max(per_step)]["views"]
            print(f"  m={arm:<6} rep{r}  step{max(per_step):<4} "
                  f"wsi R1={last['wsi_biology']['R1']:7.3f} R3={last['wsi_biology']['R3']:6.3f}  "
                  f"rna R1={last['rna_biology']['R1']:7.3f} R3={last['rna_biology']['R3']:6.3f}  "
                  f"cos={last['rna_biology']['mutual_cosine_uncentred']:7.4f} "
                  f"centred={last['rna_biology']['mutual_cosine_centred']:7.4f}", flush=True)
        steps = sorted(steps_seen or [])

    steps = sorted(steps_seen or [])
    checked = check_against_logs(arms, logs)
    print(f"\nstate/log guard: {len(checked)} rows, largest |recomputed − printed| = "
          f"{max(c['delta'] for c in checked):.2e}  (tolerance {LOG_COSINE_TOLERANCE:.0e})")

    within = {arm: within_arm(reps, steps) for arm, reps in arms.items()}
    across = {f"rep{r}": across_arms(arms, steps, f"rep{r}") for r in range(1, args.reps + 1)}
    decided = verdict(within, across, args.read_step)

    out = {
        "_": ("The RNA-view mutual cosine, CENTRED and UNCENTRED, on the same probe states as "
              "the canonical rank — the measurement draft §6.2 named as the one that would "
              "settle whether §4.10's surviving use is under strain. Three same-seed repeats of "
              "each of §5.2a's three `lr = 1e-3` arms. n = 3 per arm, one seed, one stack. Not a "
              "distribution."),
        "predeclaration": "NOTEBOOK_ENTRIES/PREDECLARED_centred_cosine_20260804T1700Z.md",
        "config": {
            "harness": "v2/research/rebase/d1_momentum_probe.py",
            "argv": {("m" + str(a).replace(".", "")):
                     f"d1_momentum_probe.py {a} 0.04 200 4096 1e-3 42 <export_dir>"
                     for a in args.arms.split(",")},
            "objective_profile": "programme_free",
            "arms": list(arms), "repeats_per_arm": args.reps,
            "seed": 42, "steps": 200, "capacity": 4096, "lr": 1e-3, "decorrelation": 0.04,
            "block": "fixed held-out probe, raw (256 held-out patients, as the model emits them)",
            "steps_measured": steps,
            "cosine_view": COSINE_VIEW,
            "statistic_sites": {
                "STATISTICS/_fold/_shape": "v2/research/rebase/p2/p2_envelope_floors.py — imported",
                "R1/R2/R3": "v2/calibra/spectral.py RANK_VARIANTS",
                "mutual cosine, both forms": "this module, one function with a `centre` flag",
            },
            "views_note": ("`geometry()` prints its rank columns on `view='wsi'` and its "
                           "`rna-rna` cosine on `view='rna'`. §5.2a and §4.10 pair the two as "
                           "though they were one block. Both views are scored here under every "
                           "statistic so that the view mismatch is separable from the centring "
                           "question — account (C) of the predeclaration."),
        },
        "state_log_guard": {
            "rule": (f"the recomputed uncentred cosine must equal the harness's printed "
                     f"`{LOG_COSINE_COLUMN}` column to {LOG_COSINE_TOLERANCE:.0e}; "
                     f"`check_against_logs` raises rather than warns"),
            "rows": len(checked),
            "max_abs_delta": max(c["delta"] for c in checked),
            "checked": checked,
        },
        "verdict": decided,
        "reps": {arm: {rep: {str(s): rec for s, rec in per_step.items()}
                       for rep, per_step in reps.items()} for arm, reps in arms.items()},
        "within_arm_floor": within,
        "across_arm": across,
    }
    path = Path(args.output).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")

    key = str(args.read_step)
    print(f"\nper repeat, never a mean — step {args.read_step}, "
          f"fixed held-out probe, lr = 1e-3, seed 42")
    print(f"{'arm':>8} {'rep':>5} | {'wsi R1':>8} {'wsi R3':>8} | {'rna R1':>8} {'rna R3':>8} | "
          f"{'cos':>8} {'cos(ctr)':>9} {'offset':>8}")
    for arm, reps in arms.items():
        for rep, per_step in sorted(reps.items()):
            w = per_step[args.read_step]["views"]["wsi_biology"]
            r = per_step[args.read_step]["views"]["rna_biology"]
            print(f"{arm:>8} {rep:>5} | {w['R1']:8.3f} {w['R3']:8.3f} | "
                  f"{r['R1']:8.3f} {r['R3']:8.3f} | "
                  f"{r['mutual_cosine_uncentred']:8.4f} {r['mutual_cosine_centred']:9.4f} "
                  f"{r['mean_offset_ratio']:8.3f}")

    print(f"\nacross-arm spread per repeat, and the within-arm floor — step {args.read_step}")
    for stat in ("mutual_cosine_uncentred", "mutual_cosine_centred", "mean_offset_ratio"):
        line = "  ".join(f"{rep}:{a[key][COSINE_VIEW][stat]['spread']:.4f}"
                         for rep, a in sorted(across.items()))
        floors = "  ".join(f"{arm}:{_spread(w[key][COSINE_VIEW][stat]):.4f}"
                           for arm, w in within.items())
        print(f"  {stat:>26}  across {line}   |  within {floors}")

    print("\nverdict:", decided["reading"])
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
