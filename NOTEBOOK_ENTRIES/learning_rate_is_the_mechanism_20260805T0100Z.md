## 2026-08-05 01:00 UTC — The two missing cells report: it is the learning rate. Momentum is neither necessary nor sufficient, and this is the first of four accounts to survive a predeclared test

**Logged:** 2026-08-05 01:00 UTC. **How obtained:** `~/e0_run/d1_diag/lr_L{5,6}.log`, read against
`~/e0_run/d1_diag/lr_L{1,2,3,4}.log`. 400 steps, `programme_free`, fixed held-out probe, capacity
4,096, one verified initialisation (**67.55**). Predeclaration
`PREDECLARED_learning_rate_test_20260804T2200Z.md` (`f68a7ac`), committed before any arm ran; the two
missing cells were predeclared a second time, with their discriminating predictions, in
`lr_test_and_decorrelation_reversal_20260804T1130Z.md` before L5 and L6 reported.

### 1. All six arms

| arm | lr | m | τ | **final eff-rank (R3)** | **rna_rna (mutual cosine)** |
|---|---:|---:|---:|---:|---:|
| L3 | 1e-3 | 0 | — | **1.06** | 0.9946 |
| L1 | 1e-3 | 0.9 | 10 | **1.05** | 0.9257 |
| **L5** | **1e-3** | **0.999** | **1000** | **1.05** | **0.5207** |
| **L6** | **4e-5** | **0** | **—** | **12.30** | **0.8199** |
| L2 | 4e-5 | 0.99 | 100 | **27.88** | 0.5436 |
| L4 | 4e-5 | 0.999 | 1000 | **35.24** | 0.3807 |

L1–L4 reproduce the earlier entry's table exactly. The `eff-rank` column of these logs is **R3** (the
row-normalised order-2 Hill number), per the header contract asserted in
`v2/research/rebase/p2/p2_floor_audit.py`; it is not canonical R1 and is not compared with any R1
number.

### 2. The predeclared reading, applied

The second predeclaration said: *"If L5 works and L6 fails, Account A survives a real test. If L5
fails and L6 works, the momentum threshold is an artefact of the learning rates it was measured at."*

**L5 fails at maximal momentum. L6 works at zero momentum.** The momentum-threshold account is
falsified on its own predeclared terms. What replaces it is the account that predeclaration named as
the third explanation the original design had omitted: **a pure learning-rate effect**.

Read as folds, which is the form that shows how completely the two knobs separate:

| | **rank (R3) moves by** | **mutual cosine moves by** |
|---|---:|---:|
| learning rate, at m = 0 (L3 → L6) | **11.60×** | 1.213× |
| learning rate, at m = 0.999 (L5 → L4) | **33.56×** | 1.368× |
| momentum 0 → 0.999, at lr 1e-3 (L3 → L5) | **1.01×** | **1.910×** |
| momentum 0 → 0.999, at lr 4e-5 (L6 → L4) | 2.87× | **2.154×** |

### 3. The honest full statement, which is what goes in the draft

The collapse is **primarily a learning-rate phenomenon**. At the training rate actually used, `2e-4`
(`d1_momentum_probe.py`'s default), momentum **does** rescue it and the seed-varied replication is
unambiguous — canonical R1 11.26 / 10.45 / 10.55 at m = 0.999 against 3.18 / 1.13 / 2.36 at m = 0,
every seed of one arm above every seed of the other
(`retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §3). But momentum is:

- **not necessary** — L6 reaches 12.30 at `4e-5` with no momentum encoder at all; and
- **not sufficient** — L5 sits at 1.05 at `1e-3` with `m = 0.999`.

**Lowering the learning rate would have solved the original problem more simply, and we did not try
it.** That sentence is the useful part of this entry and it belongs at the front of the write-up, not
in a limitation.

### 4. The sequence, which is part of the contribution

This is the **fourth** account proposed for this collapse and the **first to survive a predeclared
test**.

| # | account | how it died |
|---|---|---|
| 1 | regulariser weighting | five-arm sweep: all five arms collapse from 67.55 to 1.59–3.43, both regularisers at zero included |
| 2 | MoCo key **staleness** | three measurements: 19-step turnover; agreement does not predict rank; the freshest queue is the healthiest |
| 3 | the `τ/T` turnover criterion | predeclared, five predictions, four wrong and the discriminating one produced **no effect** (`turnover_criterion_FALSIFIED_20260804T0330Z.md`) |
| 4 | **momentum above a threshold in τ** | **this entry** — L5 fails at m = 0.999, L6 works at m = 0 |
| 5 | **the learning rate** | **survives**, on a predeclared pair of opposite-signed predictions |

Three of the four dead accounts were killed by experiments built to test them, two of those
predeclared. That is a stronger methods story than one unfalsified survivor, and the draft reports the
sequence rather than only the survivor.

### 5. What this does NOT establish, stated because the same standard applies here

- **One seed per cell.** Six arms, six runs.
- **Every arm is on the fixed held-out probe, which has no measured retraining floor** and cannot get
  one from the five exported repeats (`P2_ENVELOPE_FLOORS.json` `absent_blocks`). All four contrasts
  are therefore **unjudgeable** by this project's own criterion, and they are entered in
  `floor_audit.json` as `direction` rows with `clears: null` (ids `5.2a-*`, rows 57–60). **What
  carries the result is the predeclared sign, not a magnitude**: the rank readings fall in two bands
  (1.05–1.06 against 12.30–35.24) with nothing between them, and the contrast the falsified account
  required to be large is **1.01×**, which is an equality.
- **It does not make the anchoring story wrong; it makes it incomplete.** An account of this collapse
  now has to say why decoupling the key encoder matters at `2e-4`, is unnecessary at `4e-5` and is
  useless at `1e-3`. We do not have one.
- **It does not touch §5.4's conclusion**, which rests on a binary training outcome
  (`programme_free`: 0 of 3 seeds completing uncollapsed before the fix, 3 of 3 after) and not on why
  the fix works.

### 6. The dissociation this produces, and the alternative account we cannot exclude

At `lr = 1e-3` the three momenta are **indistinguishable in rank** — 1.06 / 1.05 / 1.05, a spread of
1.01× — while the RNA-view mutual cosine on those same three runs falls **0.9946 → 0.9257 → 0.5207**,
a factor of 1.91. Rank says three identical collapsed runs; the cosine says one of them is half as
degenerate. Across all six arms the two instruments order the two knobs **oppositely**: rank is moved
11.6–33.6× by the learning rate and at most 2.9× by momentum; the cosine is moved 1.2–1.4× by the
learning rate and 1.9–2.2× by momentum.

**One alternative account is not excluded and must travel with this observation.** The rank is
**centred** and the cosine is **uncentred**, so a difference living entirely in the mean-offset
direction would produce exactly this pattern — the same asymmetry that makes RankMe as published more
reproducible than our centred statistic on the exported artifacts (1.811× against 3.111×,
`retraining_floors_for_every_statistic_view_and_block_20260804T1220Z.md`). The measurement that would
settle it is the mutual cosine recomputed on the centred representation for those three runs; it needs
their activations, which the logs do not carry. It is recorded as a missing measurement in draft §6.2
rather than resolved in our favour.

### 7. A provenance gap, stated at the prominence it deserves

**`lr_L{1..6}.log` are not vendored into this repository.** `ablate_decorr*` and `mseed_*` were copied
into `v2/research/rebase/p2/figures/data/e0_run/d1_diag/` and every value read from them is re-parsed
from the copy by `v2/tests/test_p2_floor_audit.py`. These six were not, so the four audit rows for
this result resolve against **the draft's own table** — the weakest of the three source kinds the
audit supports, and the kind this project has twice found insufficient. **This is the paper's only
established mechanism result and it currently has its weakest provenance.** Closing it is a file copy;
nothing needs re-running. Recorded in draft §6.4 and in `P2_FIGURES.md` S9.

### 8. What changed in the paper

- **§5.2a is new** and carries the six arms, the fold table, the honest full statement and the
  four-account sequence.
- **§5.2** no longer claims a momentum threshold as a mechanism; the threshold between m = 0.9 and
  m = 0.99 is scoped to `lr = 2e-4`, and the τ/T falsification is restated as *a predicted effect
  absent* rather than *an inversion* (the ordering was read off a 1.04× difference).
- **§5.4** limit 4 is new: the binary before/after outcome supports momentum versus none **at one
  learning rate**.
- **§4.10** carries the collapse-floor dissociation with its alternative account attached; **§6.2**
  carries the missing centred-cosine measurement; **§6.4** carries the vendoring gap.
- **`floor_audit.json`** gains rows 57–60; the audit is now 60 rows. The 25-selection split is
  unchanged at **13 fail / 11 unjudgeable / 1 clears** because all four new rows are `direction`.

### Files / commits

- `~/e0_run/d1_diag/lr_L{1..6}.log` — **not vendored**, see §7.
- Predeclarations: `PREDECLARED_learning_rate_test_20260804T2200Z.md` (`f68a7ac`);
  `lr_test_and_decorrelation_reversal_20260804T1130Z.md` §1 (the two missing cells).
- Prior falsifications: `turnover_criterion_FALSIFIED_20260804T0330Z.md`,
  `d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`.
- Harness: `v2/research/rebase/d1_momentum_probe.py` (`LR` default `2e-4`).
