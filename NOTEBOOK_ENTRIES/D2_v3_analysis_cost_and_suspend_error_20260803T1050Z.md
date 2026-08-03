## 2026-08-03 10:50 UTC — D2 v3 training complete (6/6); the paired bootstrap is 96 CPU-hours, and I stalled it for two hours by misreading `ps`

**Logged:** 2026-08-03 10:50 UTC. **How obtained:** Lambda A100 `150.136.45.194`. `phase_d` stream
logs, `TRAIN_SUCCESS.json`/`artifacts/` under `~/e0_run/d2_v3/`, and `ps -o pid,stat,time` on the
analysis processes.

### Technical
**Training is done: 6/6 arms.** All three Hallmark and all three PBS arms reached 40/40 epochs and
passed the runner's G2 liveness gate and `phase_d`'s independent assertion. All six representation
artifacts exported (~18.4 MB each) under `~/e0_run/d2_v3/d2_v3_s{42,43,44}/artifacts/`. CALIBRA G4
controls ran per stream and every channel gate reported `OK`, e.g.
`[d2_i_seed42::rna_biology] adj_cca=0.8534 floor=0.4 atten=0.962 OK`,
`[d2_h_seed44::wsi_biology] adj_cca=0.6075 floor=0.4 atten=0.944 OK`.

Wall clock, three concurrent streams on a heavily shared box: **launch 01:13 → 6/6 TRAIN_SUCCESS
07:55 = 6 h 42 min**, exports complete 08:04.

**The paired bootstrap is far more expensive than the previous entry in `D2_RESULT.md` implied**
("CPU-only and takes minutes"). Measured from `phase_d`'s own single-pair runs: ~8 CPU-hours for
1 pair x 2 modes x 2000 repeats x 2 metric evaluations = 8,000 canonical correlations, i.e.
**~3.6 CPU-seconds per CCA** (each one whitens a 2766x256 and a 2766xk matrix). A 3-pair readout is
24,000 CCAs ~ 24 CPU-hours; four of them is **~96 CPU-hours**, which is >3 wall-hours even with all
30 cores. That number should be written down because it will bite the next person who plans a D2
readout as an afterthought.

**Dead end, and an expensive one.** At 08:47, with `phase_d`'s CALIBRA and my four bootstraps
fighting at load 162, I `kill -STOP`ped my five analysis processes to let CALIBRA through, intending
to resume them. At 09:07 I checked `ps -eo pid,%cpu` and saw them at ~91%, concluded they had been
resumed by the terminal cleanup when their parent ssh died, and left them alone. **They had not
been resumed.** `%CPU` in `ps` is the average over the process's whole lifetime, not its current
usage, so a stopped process keeps reporting whatever it averaged before it stopped. They sat in `T`
state for **1 h 57 min** doing nothing; I only caught it by noticing their cumulative `TIME` was
frozen at exactly `00:56:56` across two polls 20 minutes apart. Resumed at 10:45.
**The check for "is this process running" is `stat` (`R`/`S` vs `T`) or a delta in `TIME` — never
`%CPU`.**

Given the real cost, I dropped two of the four bootstraps rather than run everything at 2000
repeats:
- **kept** `untrained40` (the decisive fair contrast) and `random_control` (the negative control) —
  neither is obtainable any other way;
- **dropped** `unrestricted`, because `phase_d` independently runs exactly that comparison per
  stream at the same 2000 repeats through the same code path (only the bootstrap RNG seed differs:
  `seed + pair_index` vs `seed`), so the three per-seed `D2_PAIRED_BOOTSTRAP.json` files give it;
- **dropped** `hallmark_in_training`'s bootstrap, keeping its point estimates from the per-artifact
  readout, since it is a mechanism check rather than a claim.

### In plain terms
All six models trained and passed their health checks, and the expensive part is now the statistics,
not the training: the resampling test has to redo a heavy linear-algebra step ninety-six thousand
times. I also lost two hours to my own error — I paused my analysis jobs to let another step
through, then misread a progress column and believed they had restarted when they were still frozen.
The column I trusted reports a lifetime average, so a paused job looks busy forever.

### Meaning for the claim
No effect on the claim. The two readouts that decide it — the 40 untrained targets and the random
controls — are the two still running at full 2000 repeats. The unrestricted number, which only has
to reproduce the old headline, comes from the orchestrator's own canonical output.

### Files / commits
- `~/e0_run/d2_v3/d2_v3_s{42,43,44}/artifacts/d2_{h,i}_seed*.npz` — six frozen artifacts
- `~/e0_run/d2_v3/d2_v3_s*/calibra/` — G4 controls; `d2_v3_s*/D2_PAIRED_BOOTSTRAP.json` — unrestricted per seed
- `~/e0_run/d2_v3/D2_PAIRED_BOOTSTRAP_{untrained40,random_control}.json` — pending
