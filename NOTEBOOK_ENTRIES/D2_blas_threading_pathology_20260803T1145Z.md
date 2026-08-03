## 2026-08-03 11:45 UTC — the D2 readout was 23x slower than it needed to be: multithreaded BLAS is pathological on this box

**Logged:** 2026-08-03 11:45 UTC. **How obtained:** Lambda A100 `150.136.45.194`, `~/venv` numpy
2.2.6. Direct timing of `np.linalg.svd(a, full_matrices=False)` on a random 2766x256 matrix at
different `OMP_NUM_THREADS`, after `d2_compare` runs failed to finish in 3.5 hours.

### Technical
`d2_compare`'s bootstrap calls `top_canonical_correlation` 8,000 times per seed-pair
(2 modes x 2000 repeats x 2 arms), and each call PCA-whitens a 2766x256 representation and a
2766xk target block. I timed the real thing and got **7.5-8.7 seconds per canonical correlation**,
which puts a single-pair bootstrap at ~16.6 wall-hours and a 3-pair readout at ~50. That is absurd
for this matrix size, so I timed the SVD alone:

| `OMP_NUM_THREADS` | SVD of 2766x256 | SVD of 2766x90 |
|---:|---:|---:|
| 1 | **0.205 s** | 0.051 s |
| 4 | **4.80 s** | 0.644 s |

**Four threads are 23x slower than one.** This is not contention arithmetic — contention would cost
a small constant factor, not 23x. It is threading pathology in the installed OpenBLAS on a box that
is heavily oversubscribed: the LAPACK driver spawns and synchronises threads per call, and for a
matrix this small the barrier cost swamps the 1.8e8-flop factorisation, especially when the runnable
queue is 3x the core count.

Consequences, and both of my earlier decisions were wrong in the same direction:
- I "protected" training by capping analysis at `OMP_NUM_THREADS=6`. That cap **caused** the
  slowness rather than limiting it.
- `phase_d`'s own bootstrap subprocesses inherit no thread cap at all, so each of the three burned
  13+ CPU-hours and would not have finished today. I killed them; the three streams therefore raise
  `RuntimeError: D2 paired bootstrap failed` and write no `SUCCESS.json`. **Everything upstream of
  that step is complete and intact and was verified after the kill: 6/6 `TRAIN_SUCCESS`, all six
  `.npz` artifacts, and all three CALIBRA `calibra_gates.json` / `calibra_summary.json` /
  `task_rows.csv`.** Only the orchestrator's final convenience bootstrap was sacrificed, and it is
  being recomputed correctly.

Relaunched with `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1` and parallelism moved to **processes**:
4 target-group readouts x 3 seed-pairs = 12 concurrent single-threaded `d2_compare` jobs plus the
per-artifact readout. Predicted ~0.26 s/CCA -> ~40 min per job. Load fell from 162 to 34 immediately.

Note this also means the per-pair invocation sets `pair_index = 0` for every job, so the bootstrap
RNG seed is `42` for each seed-pair rather than `42/43/44` as in a single 3-pair invocation. The
point estimates are unaffected (deterministic); only CI endpoints move in the third decimal.

### In plain terms
The statistics step was taking eight seconds per repetition when it should take a quarter of a
second. The cause was letting the linear-algebra library use several processor cores per calculation:
the calculations are small, so the cost of organising the cores dwarfed the work, and four cores came
out twenty-three times slower than one. Telling every job to use exactly one core, and instead
running thirteen jobs side by side, fixed it.

### Meaning for the claim
No effect on any number — this is purely how fast the same arithmetic runs. It does change what is
practical: a fully stratified D2 readout is a ~40-minute job done right and a multi-day job done
wrong, and the previous D2 write-up's "it is CPU-only and takes minutes" was optimistic in a way
that would have blocked this experiment again.

### Files / commits
- `~/e0_run/d2_v3/run_d2_analysis2.sh` — single-threaded, process-parallel driver
- `~/e0_run/d2_v3/bootstrap/D2_<group>_seed<S>.json` — 12 outputs, pending
