## 2026-08-04 17:00 UTC — PREDECLARED: does the §5.2a dissociation survive centring? The reading rule, fixed before the runs

**Logged:** 2026-08-04 17:00 UTC, **before any arm was launched and before any state was scored.**
**Status at logging:** nothing run, nothing looked at. The three `lr = 1e-3` arms exist only as the six
vendored `lr_L*` logs, which carry no activations.

### 0. The question

§5.2a records that at `lr = 1e-3` the three momentum arms hold centred effective rank flat —
**1.06 / 1.05 / 1.05**, a spread of **1.01×** — while the co-measured RNA-view mutual cosine on the
same three runs falls **0.9946 → 0.9257 → 0.5207**, a factor of **1.91**. §4.10 and §6.2 both refuse to
adjudicate between two accounts of that:

* **(A) the dissociation is real** — rank is insensitive to a genuine difference in how degenerate
  these three representations are, at the exact reading where §4.10 says the collapse diagnostic is
  reliable; or
* **(B) the difference is a mean offset** — our rank is **column-centred** and the mutual cosine is
  **not**, so a difference confined to the mean-offset direction produces exactly this pattern, rank
  is right to ignore it, and the dissociation dissolves.

Centring has been the hidden variable twice already: it is also the asymmetry behind RankMe's 1.811×
floor against our centred statistic's 3.111× on identical runs (§4.1b).

**(B) is the outcome favourable to us** — it removes a claim from the paper and relieves the one place
§4.10's surviving use is put under strain. That is stated first, so that everything below is read as a
guard against wanting it.

### 1. What will be run

```
d1_momentum_probe.py {0, 0.9, 0.999} 0.04 200 4096 1e-3 42 <export_dir>
```

— the three `lr = 1e-3` arms of §5.2a (L3, L1, L5), at their own decorrelation 0.04, capacity 4,096,
seed 42 and 200-step budget, with the harness's purely additive `export_dir` argument attached so the
probe states `geometry()` already computes are saved rather than thrown away.

**Three same-seed repeats of each arm, not one.** §5.2a is one seed per cell, and this paper's own rule
is that a difference which has not been measured against a floor carries nothing. Three repeats per arm
give the cosine — centred and uncentred — a same-seed retraining floor of its own on this block, at
this learning rate, so that "moves" and "flat" are verdicts against a measured spread rather than
against an eyeball. n = 3 per arm, one seed, one stack; a floor twice over in §4.1's sense. Per repeat,
never a mean.

### 2. The statistics, and where each is defined

| quantity | definition | where it lives |
|---|---|---|
| canonical R1, R3 | `effective_rank(x, variant=RANK_VARIANTS[...])` | `v2/calibra/spectral.py` — **imported, never restated** |
| **uncentred** RNA mutual cosine | mean of the off-diagonal of `G = R̂ R̂ᵀ`, `R̂` = L2-row-normalised `rna_biology` | one new function, defined once |
| **centred** RNA mutual cosine | the same, after subtracting the **column mean** of `rna_biology` and *then* row-normalising | the same function, `centre=True` |

Both cosines come from one function with a `centre` flag, so the two cannot drift apart. No rank
statistic is computed inline anywhere; four inline-formula substitutions have been caught in this paper
and the tell was identical each time.

**Bit-level guard, declared as a stopping condition.** The uncentred cosine recomputed from the saved
states must reproduce the harness's own printed `rna-rna` column **to four decimal places at every
step of every run**. If it does not, the saved states are not the states the printed column was read
from, and **nothing in this entry may be used** — I stop and report both numbers rather than pick one.

### 3. The reading, fixed now

The uncentred cosine's published movement is **0.9946 → 0.5207**, an absolute change of **0.474**.
Because a cosine is bounded in [−1, 1] and may sit near zero, the decision is taken on the **absolute
spread across the three arms** (max − min), not on a fold; the fold is reported beside it.

| centred-cosine spread across the three arms | reading I will take |
|---|---|
| **≥ 0.20** *and* larger than the largest within-arm spread over three repeats | **(A). The dissociation is real.** Rank is missing a change a centred co-measure can see. §4.10's surviving use is under strain and the draft must say so; §6.2's missing-measurement row closes against us. |
| **≤ 0.10** *and* not larger than the largest within-arm spread, **while** the uncentred spread is ≥ 0.20 and does clear its own floor | **(B). The difference is a mean offset. The dissociation dissolves.** Rank is right to ignore it, §4.10's strain is relieved, and the claim comes out of the paper. |
| anything else | **report the magnitudes and do not adjudicate.** No third reading is invented after the fact. |

0.20 is ~40% of the published movement and 0.10 is ~20% of it; both are fixed here rather than chosen
against the numbers.

### 4. A third account this entry adds, because it is not in the draft and it must not be discovered afterwards

**The two instruments are not read on the same view.** `geometry()` takes two forward passes: the rank
columns (`R3-rank`, `CANONICAL`) are computed on **`view="wsi"`** and the `rna-rna` cosine on
**`view="rna"`**. §5.2a and §4.10 place the 1.01× rank spread and the 1.91× cosine movement side by
side as though they were two instruments on one block; they are two instruments on **two views**.

So a third account is live: **(C) the difference is real and lives in the RNA view, which the quoted
rank number does not look at.** It is distinguished by measuring canonical R1 and R3 on the
**`rna_biology`** states as well as the `wsi_biology` ones. Predeclared: if RNA-view centred rank moves
across the three arms by more than the WSI-view rank does *and* more than its own repeat floor, the
finding is a **view mismatch** and is reported as such — neither (A) nor (B) — and the draft's pairing
of a WSI-view rank with an RNA-view cosine is a prose defect to be flagged, not a result.

### 5. A named secondary, which explains whichever way the primary goes

The mean-offset account has a direct observable: **ρ = ‖column mean of `rna_biology`‖ ⁄ RMS row norm of
the column-centred `rna_biology`** — the size of the shared offset relative to the between-patient
variation around it. Account (B) predicts ρ falls sharply across the three arms while the centred
statistics do not move. It is reported for all nine runs whichever reading the primary rule selects. It
is a **secondary**: it cannot overturn §3's verdict, only explain it.

### 6. What would make me distrust a result favourable to us — i.e. distrust (B)

Written now, because the incentive runs one way.

1. **The premise must reproduce first.** If the re-run's *uncentred* RNA cosine does not reproduce the
   published ordering — m = 0.999 clearly below m = 0, by more than the three-repeat spread — then
   there was never a stable dissociation to dissolve, and the correct report is *"the §5.2a cosine
   movement does not replicate"*, which is a different and worse finding. Likewise if the WSI-view rank
   is **not** flat at ~1.01× across the three arms, §5.2a's premise is not reproduced and no reading of
   the centred cosine is licensed.
2. **Flatness at a degenerate value proves nothing.** If the centred R3 is ~1.05 in all three arms, the
   centred representation is essentially one-dimensional, and for a one-dimensional family the mean
   off-diagonal centred cosine is a function of the sign structure of a single coefficient — it is near
   0 by construction and its flatness is **partly entailed by the rank number itself**. If that is what
   I find, the honest statement is that the test was **partly circular** and (B) is *consistent with*
   rather than *established by* it. I will say so in those words rather than bank the favourable
   reading.
3. **No floor, no verdict.** If the within-arm three-repeat spread of the centred cosine is of the same
   order as its across-arm spread, "flat" is unmeasured, not measured, and the row stays open.
4. **Centring must not be doing the work alone.** If the centred *cosine* is flat but centred *rank on
   the RNA view* moves, then centring did not dissolve the difference — this one statistic did — and
   (B) is refused.
5. **The uncentred cosine must clear its own floor before its movement counts as a movement.** The
   0.474 change is quoted from one seed per cell; if three repeats put it inside its own spread, then
   §5.2a's headline dissociation was inside the noise all along and both (A) and (B) are moot.

### 7. What this cannot settle whatever it returns

Three repeats, one seed, one stack, one learning rate, one 200-step budget, one architecture. It settles
what the difference between those three specific runs consists of. It does not establish anything about
centred rank's sensitivity in general, and no sentence written from it may say that it does.

### 8. Commit

This entry is committed **before** the workspace is built and before any run is launched. Its commit
hash is the predeclaration.
