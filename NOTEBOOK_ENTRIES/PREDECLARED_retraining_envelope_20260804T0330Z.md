## 2026-08-04 03:30 UTC — PREDECLARATION: how the retraining envelope will be read, including the outcome that would flatter us

**Logged:** 2026-08-04 03:30 UTC. **Committed before the five retraining repeats report.** The runs
are queued behind the seed extension and the momentum replication; no envelope number exists yet.

### Why this needs predeclaring

P2 §4.1's headline is that six of seven between-arm rank differences fall inside a **2.69×**
retraining envelope. That envelope rests on **one retraining pair** — seed 42, retrained once. The
central claim has the same n=1 defect as the momentum sweep, applied to the load-bearing number
rather than a supporting one.

Measuring it properly is decisive because D1's own rank ratios are **2.02×, 3.09×, 1.68×**, so two of
three already sit inside the current envelope. **One of the possible outcomes is suspiciously
convenient for us**, and the reading is fixed here so it cannot be chosen afterwards.

### The convenient outcome, and why it is not a rescue

If the envelope swallows D1's rank differences, the tempting write-up is *"the 3/3 selection success
dissolves, so the counter-evidence goes away."* **That is wrong and is foreclosed here.**

A comparison inside the noise floor cannot be evidence **for** rank's reliability. It equally cannot
be evidence **against** it. The honest statement is that **D1 does not resolve the question in either
direction**, and §4.7 must say so while keeping the necessity result reported at its current
prominence, flagged as *not resolvable* rather than deleted.

**What would make a wide envelope a genuine demonstration is the asymmetry, not the width.** On the
same three seeds the channel differences have patient CIs excluding zero 3/3 (−0.0705, −0.0863,
−0.0961). So if the envelope swallows the rank ratios *while the channel intervals still exclude
zero*, the identical comparison shows **the channel resolvable and rank not** — which is the paper's
thesis measured on one pair of arms. That framing holds **only if both sides are reported**; quoting
the rank half alone converts a demonstration into a convenience.

### The four outcomes, fixed in advance

| envelope | reading |
|---|---|
| **< 1.68×** | All three D1 ratios clear it. Rank's 3/3 stands as a real selection success; §4.7 keeps it at full strength. **Worst outcome for the paper's claim** and to be reported as such. |
| **1.68× – 3.09×** | Some seeds resolvable, some not. **Report per seed, do not pool.** Only the seed(s) clearing the envelope carry weight. |
| **> 3.09×** | No D1 rank difference is resolvable. D1 is **uninformative about rank**; the channel remains resolvable; the asymmetry is the finding. Explicitly **do not** claim the necessity result is refuted. |
| **> the six-of-seven threshold in §4.1** | Our own headline count changes, in our favour. To be reported with **the same scepticism we would apply to a result going the other way**, including that it rests on this one measurement. |

### The number will be a FLOOR, twice over

1. **`programme_only` is the stable arm.** Its channel varies 1.018–1.026× across seeds where
   `programme_free` varies 1.056×, and its step-200 tripwire rank spans 1.003× against
   `programme_free`'s 6.05×. Measuring the envelope on the stable arm understates it.
2. **Same-seed repeats exclude seed variation entirely.** The only source of variation is GPU
   non-determinism. The seed extension measures the other axis separately and both belong in the paper.

So a floor that already swallows two of three arm separations is a stronger statement than a point
estimate that happens to. It must be reported as a floor, not as *the* envelope.

### Channel measured alongside rank, on the same runs

Each repeat is exported anyway, so both quantities are computed per repeat: canonical effective rank
and top-CCA on the untrained-40 block, residualised identically to `d2_compare`. That gives **rank
spread against channel spread on identical configurations, in one table**, rather than requiring a
reader to compare across tables — the asymmetry read directly off the measurement that produces it.

### Statistic and block, stated because a fourth definition was recently found hiding under the others

- **Statistic:** canonical Roy & Vetterli order 1 (centred, rows at own norms), imported from
  `v2.calibra.spectral`. Nothing computed inline.
- **Block:** exported artifact `wsi_biology`, held-out test partition, cancer + pooled-TSS
  residualised, top-CCA at 16 components.

### Files / commits

- To be produced: `~/e0_run/d1_envelope/rep{1..5}/`, `rep{1..5}.npz`, readout log
- Chain: `~/chain_retrain_envelope.sh`
- Prior predeclarations under the same discipline: `PREDECLARED_turnover_criterion` (falsified),
  `PREDECLARED_D1_necessity_test` (outcome O2 fired, against us)
