## 2026-08-04 01:30 UTC — PREDECLARATION: the critical momentum should track queue turnover

**Logged:** 2026-08-04 01:30 UTC. **Written and committed BEFORE the runs are launched**, so the
predictions cannot be fitted to the results. **How obtained:** arithmetic from measured quantities;
the experiment it predicts is described below and had not been started at commit time.

### The criterion under test

An EMA key encoder with momentum `m` has a time constant `τ = 1/(1−m)` steps. A queue of `capacity`
slots fed `batch` patients per step turns over completely every `T = capacity / batch` steps.

**Claim: the key encoder must change more slowly than the queue refreshes — `τ > T` — or there is no
independent reference frame and collapse is free.**

If `τ < T`, the key encoder has caught up to the query encoder before the queue has been replaced, so
every key in the queue was written by a key encoder that already agrees with the current query
encoder, and a transformation applied to all patients at once moves queries and keys together again.

The claim is **dimensionless**: what matters is the ratio `τ/T`, not `m` and not `capacity`
separately. Measured batch size is **214 patients** (from 4309.4 effective negatives against a 4,096
queue), and batch size is set by the token budget, independent of capacity.

### Already-measured points, restated in the ratio

| capacity | T | m | τ | τ/T | observed (eff-rank @300) |
|---:|---:|---:|---:|---:|---|
| 4096 | 19.1 | 0 | 0 | 0 | **fails** (3.32) |
| 4096 | 19.1 | 0.9 | 10 | **0.52** | **fails** (2.70) |
| 4096 | 19.1 | 0.99 | 100 | **5.2** | **works** (6.01) |
| 4096 | 19.1 | 0.999 | 1000 | 52 | **works** (7.33) |

The boundary lies between τ/T = 0.52 and 5.2. The criterion predicts it is at **τ/T = 1**.

### Predictions, made before running

| # | capacity | T | m | τ | τ/T | **predicted** |
|---|---:|---:|---:|---:|---:|---|
| P1 | 2048 | 9.6 | 0.9 | 10 | **1.04** | marginal — weakly works |
| P2 | 2048 | 9.6 | 0.95 | 20 | **2.1** | works |
| P3 | 4096 | 19.1 | 0.95 | 20 | **1.05** | marginal |
| P4 | 8192 | 38.3 | 0.95 | 20 | **0.52** | **fails** |
| P5 | 8192 | 38.3 | 0.99 | 100 | **2.6** | works |

**The sharpest single prediction is P4.** `m = 0.95` at capacity 8192 has τ/T = 0.52 — the *same
ratio* as `m = 0.9` at capacity 4096, which failed at 2.70. If the criterion holds, P4 must fail at
roughly that level, **even though m = 0.95 is a larger momentum than a value that works at a smaller
capacity**. That is a prediction no account based on "more momentum is better", on staleness, or on
"decoupling at all" would make: all three predict P4 works.

**The second sharpest is the pair P2 vs P4**, which hold `m = 0.95` fixed and vary only capacity —
predicted to work at 2048 and fail at 8192, i.e. the *same momentum* changing verdict purely because
the turnover moved.

### Falsification conditions, also predeclared

- If P4 **works**, the time-constant account is wrong and τ/T is not the controlling variable.
- If the threshold does not move with capacity — e.g. m = 0.95 behaves the same at 2048, 4096 and
  8192 — the account is wrong.
- If outcomes track `m` alone regardless of capacity, the account is wrong and "more momentum is
  better" is the better description.

Any of these returns the mechanism to an open question, which is an acceptable outcome and preferable
to a third explanation that later has to be withdrawn. This is the third mechanism proposed for this
collapse; the first two (regulariser weighting, MoCo staleness) were both falsified by measurement,
and this predeclaration exists so the third is held to the same standard rather than a lower one
because it is ours and it is prettier.

### Reading rule, predeclared

"Works" = centred effective rank ≥ 5 at step 300 on the 256-patient held-out probe; "fails" = ≤ 3.5;
between those is "marginal". Those thresholds are taken from the existing 4096 series, where the
working arms sit at 6.0–7.3 and the failing arms at 2.7–3.3, and are fixed here so the verdict is not
chosen after seeing the numbers.

### Files / commits

- To be produced: `~/e0_run/d1_diag/turn_cap{2048,4096,8192}_m*.log`
- Harness: `v2/research/rebase/d1_collapse_causal_test.py` / `~/ws_d1/momentum_test.py`
