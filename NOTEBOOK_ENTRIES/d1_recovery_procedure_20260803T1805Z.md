## 2026-08-03 18:05 UTC — How to recover the D1 run if the instance stops (read this before restarting anything)

**Logged:** 2026-08-03 18:05 UTC. **How obtained:** inspection of `V2Trainer.save_checkpoint` /
`load_checkpoint`, `run_d1`, and the live run state on the A100 (`150.136.45.194`). Written while D1
is still training, because this project has already lost a completed six-hour sweep to an instance
stop and the recovery path should not have to be reconstructed under time pressure.

### Technical

**What survives an instance stop.** `/home/ubuntu` is wiped; `~/e0_run` is a symlink to
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/e0_run` and persists. The D1 run root
is `~/e0_run/d1_v1`, so every checkpoint and metric file survives. Each run writes `last.pt` (~292 MB)
and appends to `train_metrics.jsonl` **every epoch**, so at most one epoch of work is ever at risk.

**What `--resume` restores.** `V2Trainer.load_checkpoint` restores the model, the optimiser, the
scheduler, the torch and numpy RNG states, `programme_memory`, **and `biology_memory`** — the paired
WSI/RNA queue. That last one matters: `programme_free`'s entire contrastive objective is defined
against that queue, so a resume that silently reset it would restart the arm with an empty negative
set. It does not.

**The trap.** `run_d1` refuses to reuse a non-empty output directory —

```
raise RuntimeError(f"refusing stale D1 output directory {run_dir}; use a new run root")
```

— so **re-running `phase_d d1 --execute` after a crash will not resume; it will refuse, and pointing
it at a fresh run root would discard all completed work.** That guard is correct (it exists to stop
half-finished arms being compared), but it means recovery is manual.

**Recovery procedure.**

1. Read `~/e0_run/d1_v1/D1_LAUNCH_PLAN.json`. It contains the exact `argv` for all six runs, along
   with each one's arm, profile and seed. This file is the authority — do not retype the command.
2. For each run whose directory has a `last.pt` but no `TRAIN_SUCCESS.json`, re-launch its `argv`
   directly, appending `--resume <run_dir>/last.pt`. Runs that already have `TRAIN_SUCCESS.json` are
   complete; leave them alone.
3. Cap threads, and prefer `setsid` so the job outlives the ssh session:

```
cd ~/ws_d1 && PYTHONPATH=$HOME/ws_d1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  setsid nohup ~/venv/bin/python -m morpheus.v2.runner <argv from the plan> \
  --resume <run_dir>/last.pt > <run_dir>/resume.log 2>&1 < /dev/null &
```

4. G2.6 re-runs on resume. That is fine and is a feature — it re-proves the memorisation gate against
   the restored state — but budget ~40 minutes per run for it at 2400 steps.
5. The stages `run_d1` would have done after training must then be run by hand, in order, from
   `run_d1`'s own definitions: per-run `TRAIN_SUCCESS.json`, `morpheus.v2.export` per arm/seed,
   `morpheus.v2.calibra.run_calibra` with `--require-rna-positive-control --require-channel-gates`,
   then `d2_compare`. Do not skip the CALIBRA controls; without them no D1 number is valid.

**Launch hygiene that made the current run recoverable**, worth repeating on the next one: `setsid
nohup` meant the whole pipeline survived an ssh disconnect that killed the launching shell — the
`phase_d` parent was reparented to init and kept its three runners as children, so the export /
CALIBRA / bootstrap stages it owns are still queued rather than orphaned. Verified live:
`phase_d` PID 126670 alive, PPID re-parented, three runners as its children.

**One thing to do differently next time:** put the four thread-cap variables in `phase_d`'s own
environment at launch. The training subprocesses and, more importantly, the CALIBRA and bootstrap
stages inherit it, and those are joblib workloads on a shared box.

### In plain terms

If the machine stops, nothing is lost except at most the epoch in progress — the saved state lives on
network storage that survives, and it includes the memory queue the programme-free arm needs, not
just the model weights.

The one thing not to do is re-run the original launch command. It is deliberately built to refuse to
touch a directory that already has work in it, so it will either stop or, if pointed somewhere fresh,
quietly start over from nothing. Recovery means re-issuing the individual training commands from the
saved plan file with a resume flag, then running the measurement steps by hand.

### Meaning for the claim

Nothing scientific. It protects six GPU-hours of the run that supplies P2's primary evidence and the
biology-head-without-programme-supervision arm, and it removes a plausible way to accidentally
destroy that run while trying to save it.

### Files / commits

- `~/e0_run/d1_v1/D1_LAUNCH_PLAN.json`, `D1_PAIR_MANIFEST.json`, `~/e0_run/d1_v1/d1_*/last.pt`
- `v2/training.py` — `save_checkpoint` / `load_checkpoint`
- `v2/research/rebase/phase_d.py` — `run_d1`, the stale-directory guard and the post-training stages
