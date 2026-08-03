## 2026-08-03 01:50 UTC — D2 v3 relaunch: environment reconstituted, 5 of 6 arms must be retrained

**Logged:** 2026-08-03 01:50 UTC. **How obtained:** Lambda A100 `150.136.45.194` (80 GB, 30 cores,
216 GB RAM). `morpheus.v2.research.rebase.phase_d d2 --execute`, three concurrent run roots under
`~/e0_run/d2_v3/` (persistent NFS).

### Technical
Two environment problems had to be fixed before any GPU work, both consequences of the
`/home/ubuntu` wipe:

1. **`~/morpheus-rebase` had lost its `.git`.** `phase_d._require_clean_worktree` runs
   `git status --porcelain` and fail-closes G0.2, so execution was impossible. Rather than trust the
   tree blindly, I hashed all 176 `.py`/`.yaml`/`.json` files under `v2/ src/ configs/ tests/` on
   both the Lambda box and the local repo at HEAD: **byte-identical, 0 differing, 0 only-on-one-side.**
   Only then did I `git init` + commit the tree (`9a4f307`) to restore the gate. The gate is
   therefore satisfied by a tree provably equal to local HEAD `7600179`, not by a rubber stamp.
2. `~/e0_run` is confirmed a symlink to
   `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/e0_run`, verified persistent by
   round-tripping a file through the absolute NFS path. **Every D2 v3 output is written under
   `~/e0_run/d2_v3/`.** Nothing lands on `/home/ubuntu`.

Run plan. Seed-42 Hallmark survives and re-exports exactly (separate entry), but the seed-42 PBS
checkpoint is a partial epoch-13 snapshot, so **5 arms genuinely need training**. I am nonetheless
retraining all 6 (three `phase_d` run roots, `--seeds 42`, `43`, `44`, one per process) because
`phase_d` trains both arms per root, and the streams run concurrently — so retraining H42 costs
**zero additional wall clock** (stream s42 would otherwise idle while s43/s44 each run two arms) and
buys a fully gated, self-consistent seed-42 pair with its own `TRAIN_SUCCESS`, liveness and CALIBRA.
The recovered H42 export is kept alongside as an independent reproducibility control.

Concurrency measurement, which is the interesting part. The handoff expected these runs to be
data-loading bound at ~6 GB of an 80 GB card, so concurrency should have been nearly free. It is
not, on this box, right now:

- 3 concurrent runners: ~20 GB GPU total, 89% reported utilisation, **~3.5 min/epoch each**
  → aggregate **0.86 epochs/min**.
- Previous serial run (from `d2_i_seed42`'s surviving metrics: 14 epochs in ~27 min): ~1.9 min/epoch
  → aggregate **0.53 epochs/min**.

So concurrency is buying about **1.6×**, not 3×. Cause is CPU, not GPU: load average peaked at 67 on
30 cores. Two contributors, one mine and one not. Mine: my own `d2_readout.py` was silently taking
**16 cores** of BLAS threads for its SVDs — killed, and it will be re-run niced and thread-capped
once training is done. Not mine: another agent is concurrently running
`morpheus.v2.calibra.induced_correlation_sweep` (8 loky workers), a second 6-worker joblib pool, and
`morpheus.v2.build_pbs_targets` on the same box. Each of my three runners is only getting ~80% of one
core's worth of scheduling for its dataloader. Reported GPU "utilisation" of 89% is time-any-kernel-
resident across three interleaved processes and does not mean the card is saturated.

Also: the 800-step G2.6 `_overfit_programme_only_actual` liveness gate runs before epoch 0, costing
~15 min of pre-epoch time per arm. That is per-arm, so it is paid 6 times, ~1.5 h of the total.

### In plain terms
The machine had been reset, so the code checkout had lost its identity papers and the safety gate
that refuses to run experiments from unversioned code was blocking. Instead of switching the gate
off, I proved file-by-file that the code on the machine is exactly the code in the repository, then
restored its papers. Running three trainings at once helps, but only about half as much as hoped,
because the bottleneck is the processor feeding data, not the graphics card — and another agent is
using the same processors.

### Meaning for the claim
No effect on the claim yet; this is provenance and scheduling. The one thing worth carrying forward
is that **G0.2 was satisfied by verified code identity rather than bypassed** — if anyone asks
whether these runs came from committed code, the answer is yes and it is checkable.

### Files / commits
- Lambda repo commit `9a4f307` (tree verified == local `7600179`), `3853d40` (d2_compare fix)
- `~/e0_run/d2_v3/d2_v3_s{42,43,44}/` — run roots, each with `D2_PAIR_MANIFEST.json`
- `~/e0_run/d2_v3/logs/` — launch and export logs
