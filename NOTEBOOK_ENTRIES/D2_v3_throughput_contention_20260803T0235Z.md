## 2026-08-03 02:35 UTC — D2 v3 throughput is co-tenant-bound, not GPU-bound (and a niced-analysis dead end)

**Logged:** 2026-08-03 02:35 UTC. **How obtained:** Lambda A100 `150.136.45.194`;
`~/e0_run/d2_v3/progress.sh` (5-minute sampler writing `logs/progress.log`), `ps`/`uptime`/`nvidia-smi`.

### Technical
The handoff's expectation was that D2 arms are data-loading bound at ~6 GB of an 80 GB card, so
running seeds concurrently should be close to free. Measured steady state with three concurrent
`phase_d` streams:

| | epochs/min per stream | aggregate | GPU mem | GPU util | load avg |
|---|---:|---:|---:|---:|---:|
| previous serial run (from surviving `d2_i_seed42` metrics) | 0.53 | 0.53 | ~6 GB | — | — |
| D2 v3, 3 concurrent | ~0.18 | **~0.55** | ~20 GB | 82–89% | **87–93** |

So concurrency is currently buying **essentially nothing** — about 5.5 min/epoch per stream against
1.9 min/epoch serial. The card is not the constraint (20 of 80 GB; "89% utilisation" is
time-any-kernel-resident across three interleaved processes, not saturation).

The constraint is CPU, and most of it is not mine. A co-tenant agent is running **31–38 loky worker
processes** (`morpheus.v2.calibra.induced_correlation_sweep` at 8 workers, a second 6-worker joblib
pool, `morpheus.v2.build_pbs_targets`), each at 150–200% CPU, on a 30-core box. My three runners sit
in `R` state at ~64–70% of a single core each — they are runnable and starved, not blocked on I/O
(checked: threads are `Sl`, not `D`, so this is not NFS stalling on the H-Optimus patch store).

Two of my own mistakes contributed and are fixed:
1. My first `d2_readout.py` invocation silently took **16 cores** of BLAS threads for its SVDs.
   Killed; all subsequent analysis runs pin `OMP/OPENBLAS/MKL_NUM_THREADS`.
2. **Dead end:** I then re-ran it at `nice -n 19` with 1 thread, reasoning it would be free. Under
   load ~90 it got **3% of a core** and had produced nothing after 43 minutes. Killed. The
   permutation nulls are folded into the end-of-run analysis instead, at normal priority with 4
   threads. Lesson: `nice 19` on a box at 3× oversubscription is not "low priority", it is "never".

I am **not** restarting the streams to run all six arms at once. It would halve the sequential depth
but requires hand-rolling `phase_d`'s post-training export/CALIBRA/TRAIN_SUCCESS path, and the
brief's priority is a defensible result over a fast one. Nor am I renicing the co-tenant's jobs.
Projected completion at the current rate: H arms ~05:10 UTC, I arms ~08:45 UTC.

### In plain terms
Running three trainings side by side was supposed to be nearly free because each one barely touches
the graphics card. It isn't, because another agent is running about thirty-five processes on the
same thirty processors, so all three of my trainings are queueing for a slice of CPU to feed the
card. The graphics card is mostly idle and waiting to be fed. Nothing is broken; it is simply a
shared machine, and the honest wall-clock number will reflect that rather than the hardware.

### Meaning for the claim
No effect on the claim — this is scheduling, and every scientific gate is unchanged. It does mean
the wall-clock figure I report is a co-tenancy figure, not a hardware figure, and should not be
quoted as "what a D2 sweep costs on an A100".

### Files / commits
- `~/e0_run/d2_v3/logs/progress.log` — 5-minute throughput sampler (persistent)
- `~/e0_run/d2_v3/run_d2_analysis.sh` — end-of-run driver, thread-pinned
