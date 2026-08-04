## 2026-08-04 03:00 UTC — D1's A3 passes on arm difference and only on arm difference; and the effect we measure is half the instrument's own floor

**Logged:** 2026-08-04 03:00 UTC. **How obtained:** `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_*.json`,
produced by `d1_audit.py` on `~/ws_d1` (the one workspace verified equal to HEAD: 0 code files
differing, 0 missing). Statistic and block are stated for every number below.

### A3 — the negative control

| seed | `programme_only` | `programme_free` | Δ | patient CI | cancer CI |
|---|---:|---:|---:|---|---|
| 42 | 0.4728 | 0.4504 | −0.0224 | [−0.0594, +0.0068] | [−0.0986, +0.0895] |
| 43 | 0.4810 | 0.4738 | −0.0072 | [−0.0480, +0.0278] | [−0.0822, +0.1137] |
| 44 | 0.4747 | 0.4425 | −0.0322 | [−0.0807, +0.0077] | [−0.1055, +0.0768] |

*(90 `random_control` targets, top canonical correlation at 16 components, cancer + pooled-TSS
residualised, 2,766 held-out patients.)*

**Verdict, in the only form it supports: the instrument does not manufacture an arm difference out of
noise; the absolute level on random controls is high and is separately explained.**

Not "controls are at chance", which is false and is how an unqualified pass would be read.

**Why the absolute level is not chance.** I had offered two readings and said the data did not choose
between them. It does. The within-cancer permutation null on this cohort is **0.140**
(`D2_stratified_result`, 200 draws) — that is the 16-component capacity floor at n=2,766. Random gene
sets scoring 0.44–0.48 sit at **~3.2× that floor**, so they carry real signal and the
"it is just CCA bias" reading is out.

*(Note on provenance: the traced figure in the notebook is **0.140**; 0.147 was quoted to me. The
difference does not affect the conclusion — 0.44/0.147 = 3.0× versus 0.44/0.140 = 3.2× — but the
value used here is the one that traces to a file.)*

**And this is not a new problem, which is the reassuring part.** T1.4 measured covariate-matched
random gene sets reproducing **76–82%** of the per-target channel, invariant to dilution, and
concluded it is a property of the readout rather than of image quality. D1's random controls sit at
**77–85%** of the real-target level (0.4504/0.5412 … 0.4747/0.6087). Two independent estimates, two
different axes, same quantity: **76–82% and 77–85%.** The agreement is worth stating explicitly,
because a reader meeting either number alone will suspect a one-off.

### The effect is half the instrument's floor — and this belongs in P2

The number that matters most here is a comparison P2 has not yet made:

| quantity | value |
|---|---:|
| real-target channel, `programme_only` | 0.6117 |
| random-control channel, `programme_only` | 0.4728 |
| **real-versus-noise margin** | **0.1389** |
| **arm difference on real targets (D1 seed 42)** | **0.0705** |
| ratio | **0.51** |

**The effect D1 is built to detect is roughly half the margin between real targets and random ones.**

That is the same shape as P2 §4.1's envelope argument, arriving from the *channel* side rather than
the *rank* side. §4.1 says rank differences are smaller than rank's own retraining noise; this says
the *information difference we are reporting* is small relative to the readout's own floor.

Stating it strengthens the paper rather than weakening it, and it pre-empts the obvious objection —
that we apply a strict standard to RankMe and a loose one to ourselves. We should be seen to apply the
envelope test to our own headline and report where it lands, not only to the metric we are
criticising.

### D1's paired bootstrap exists — quote both bootstraps

The draft's `[D1 PAIRED BOOTSTRAP PENDING]` is stale.

| seed | Δ (free − only) | patient CI | cancer-cluster CI |
|---|---:|---|---|
| 42 | −0.0705 | [−0.0938, −0.0444] | [−0.0957, −0.0180] |
| 43 | −0.0863 | [−0.1186, −0.0522] | [−0.1386, **+0.0006**] |
| 44 | −0.0961 | [−0.1314, −0.0618] | [−0.1535, −0.0016] |

**3/3 decisive on the patient bootstrap; 2/3 on the cancer-cluster bootstrap**, seed 43 crossing zero
at +0.0006. The cluster bootstrap is the conservative one — it resamples cancers rather than patients
and so respects the fact that patients within a cancer are not independent. Quoting only the patient
result would be the exact kind of selective reporting this project has spent two days catching.

This is predeclared outcome **O2** (`PREDECLARED_D1_necessity_test`): *rank vindicated on this pair,
reported as a limitation with the same prominence a confirmation would have had.*

### The single-seed momentum sweep was a defect, not a decision

P2 §5's momentum choice rests on 2.81 vs 7.42 (R3, held-out, 500 steps) — a 2.64× ratio from **one
seed per momentum value**, which §4.1's own 2.69× retraining envelope would disqualify.

**That was not a design decision. `momentum_test.py` had `SEED = 42` hardcoded**, with no way to vary
it; every run in the sweep was necessarily the same seed. "We only ran one seed" and "the harness
could only run one seed" are different admissions and the second is the useful one, because it says
the limitation was invisible to whoever read the sweep rather than accepted by them.

Fixed: the harness now takes a seed argument and reports the **canonical** Roy & Vetterli order-1
statistic beside the R3 participation ratio on the same states, so the replication is comparable to
what P2 quotes rather than to what the sweep happened to use. Six runs — m ∈ {0, 0.999} × seeds
{42, 43, 44}, 500 steps, the sweep's own regime — are armed behind the seed extension.

### In plain terms

The negative control does what it is there for: the two arms score the same on nonsense targets, so
the difference we report between them is not an artefact of the measuring device. But the nonsense
targets themselves score far above the floor — three times it — so "random controls pass" must not be
written as "random controls are at chance". A different experiment already found the same thing from
another direction, and the two agree closely.

Separately, the difference between our two arms is about half the gap between real targets and random
ones. We are criticising a metric for having differences smaller than its own noise; intellectual
honesty requires us to say where our own effect sits against our own noise, and it sits at half.

### Files / commits

- `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json`, `..._RANDOM_CONTROL.json`
- Null: `D2_stratified_result_20260803T1210Z.md` (0.140, 200 draws)
- Agreement: `t12_t14_t16_t17_calibra_ledger_20260803T0230Z.md` (76–82%)
- Predeclaration: `PREDECLARED_D1_necessity_test_20260803T2300Z.md` (outcome O2)
- Harness: `~/ws_d1/momentum_test.py` — seed argument and canonical statistic added
