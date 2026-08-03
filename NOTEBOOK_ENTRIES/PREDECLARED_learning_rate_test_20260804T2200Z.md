## 2026-08-04 22:00 UTC — PREDECLARATION: does the momentum threshold live in steps, or in parameter-space drift?

**Logged:** 2026-08-04 22:00 UTC. **Written and committed BEFORE the runs**, as with the turnover
criterion — which this same discipline falsified. **How obtained:** design only; no run started at
commit time.

### The question

A momentum key encoder rescues the biology representation, monotonically in `m`, with a threshold
between `τ = 20` (fails) and `τ = 100` (works) where `τ = 1/(1−m)` steps. That threshold sat in the
**same place** at capacities 2048, 4096 and 8192, which killed the ratio account (`τ/T`). What is left
is an absolute threshold in `τ` — and two incompatible readings of it:

**Account A — steps.** What matters is the *number of steps* of lag. `τ` is measured in steps and the
learning rate does not change it, so **the critical `m` is unchanged by learning rate.**

**Account B — parameter-space drift.** What matters is how far the key encoder lags in *parameter
space*, which is roughly `τ × (drift per step)`, and drift per step scales with learning rate. So
**`τ_critical ∝ 1/lr`**: raise the learning rate and less lag suffices; lower it and more is needed.

The learning rate is the one manipulation that separates them, because it moves per-step drift without
touching either `τ` or `T`.

### Predictions, made before running

Baseline (already measured, lr 2e-4): `m = 0.9` (τ=10) **fails**, `m = 0.99` (τ=100) **works**.
Test at 5× and 1/5× that rate.

| # | lr | m | τ | **Account A predicts** | **Account B predicts** |
|---|---:|---:|---:|---|---|
| L1 | 1e-3 (5×) | 0.9 | 10 | **fails** (τ still < threshold) | **works** (τ_crit falls ~5× to 4–20) |
| L2 | 4e-5 (1/5×) | 0.99 | 100 | **works** (τ still > threshold) | **fails** (τ_crit rises ~5× to 100–500) |

**L1 and L2 are opposite-signed discriminators**, which is the point: no single-direction artefact can
satisfy both accounts. If L1 works *and* L2 fails, Account B. If L1 fails *and* L2 works, Account A.
Any other combination falsifies both, and that is a real possible outcome.

**Two controls, because changing the learning rate changes more than the lag.**

| # | lr | m | purpose | both accounts predict |
|---|---:|---:|---|---|
| L3 | 1e-3 | 0 | is the high rate healthy at all without momentum? | fails |
| L4 | 4e-5 | 0.999 | is the low rate merely too slow to learn in 200 steps? | works |

**L4 is the load-bearing control.** If L2 fails *and* L4 also fails, the low-rate arm simply has not
trained yet and L2 says nothing about lag — the test is confounded and must be rerun at more steps.
L2 only discriminates if L4 works.

### Reading rule, fixed in advance

Centred effective rank on the held-out probe at step 200, as before. **The primary read is the
contrast against the learning-rate-matched control, not the absolute value**, because changing the
learning rate may shift the whole scale. Concretely:

- L1 counts as "works" if it lands **closer to the working band than to L3**, its matched control.
- L2 counts as "fails" if it lands **closer to L3-like collapse than to L4**, its matched control.

Where the arms are unambiguous the absolute bands from the earlier sweep still apply: ≥ 5 works,
≤ 3.5 fails.

### What I expect, stated so it can be wrong

I weakly favour **Account B**. A fixed number of *steps* is a strange thing for the dynamics to care
about; what should matter is how far apart the two encoders are, and that is a distance in parameter
space, which the learning rate scales directly. But the earlier data is genuinely more consistent with
A on its face — the threshold sat at the same `τ` across a 4× range of capacity — so I am not
confident, and this is exactly the situation where writing the prediction down first is worth the
five minutes.

**This is the fourth account proposed for this collapse.** Regulariser weighting, MoCo staleness and
the `τ/T` ratio were each falsified by measurement, the last by an experiment built to test it. If
both L1 and L2 come out against their accounts, the honest report is that we have four falsified
explanations and a robust, unexplained empirical fix — which is a stronger methods section than one
unfalsified story, and considerably stronger than a story we stopped testing once it survived.

### Files / commits

- To be produced: `~/e0_run/d1_diag/lr_lr{1e-3,4e-5}_m*.log`
- Harness: `~/ws_d1/momentum_test.py` (needs an `lr` argument added)
- Prior falsification under the same discipline: `PREDECLARED_turnover_criterion_20260804T0130Z.md`,
  `turnover_criterion_FALSIFIED_20260804T0330Z.md`
