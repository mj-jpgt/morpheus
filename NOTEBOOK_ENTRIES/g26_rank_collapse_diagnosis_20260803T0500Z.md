## 2026-08-03 05:00 UTC — G2.6's third defect is rank-1 collapse of the biology head under the contrastive objective itself, not under any regulariser

**Logged:** 2026-08-03 05:00 UTC. **How obtained:** `~/ws_d1/diag_{a,b,c,d}.py` on the A100
(`150.136.45.194`), real cohort, `paired_split_maximal`, seed 42, the same fixed 16-patient batch the
gate uses. Logs under `~/e0_run/d1_diag/`.

### Technical

The brief attributed `programme_free`'s G2.6 failure to `feature_decorrelation` (whose global minimum
is total collapse) with `variance_floor` as an insufficient counterweight. **That attribution is
wrong.** Four measurements, in increasing order of how much was stripped away:

| configuration | decorrelation | variance | memory queue | full-consistency | outcome |
|---|---|---|---|---|---|
| shipped `programme_free` gate | 0.04 | 0.01 | 64 frozen keys | 1.0 | collapse, pinned ln(16) |
| decorrelation and variance BOTH zeroed | 0 | 0 | 64 frozen keys | 1.0 | **collapse by step 10**, wsi-wsi cos 1.0000 |
| memory queue removed as well | 0 | 0 | none | 1.0 | collapse by step 100 |
| full-consistency removed as well (clean in-batch InfoNCE) | 0 | 0 | none | 0 | **effective rank 12.88 → 1.00 by step 50** |

The last row is the finding. A plain symmetric in-batch InfoNCE on `z_biology`, with every other
term in the objective deleted, drives the biology head to effective rank 1.00 within 50 steps and
then stops moving. From `diag_d`:

```
    0 loss 3.0762 acc 0.188 pos 0.1235 worst-neg 0.2186 min-margin -0.2190 eff-rank 12.88
   25 loss 2.4322 acc 0.125 pos 0.9344 worst-neg 0.9872 min-margin -0.1103 eff-rank  1.21
   50 loss 2.3958 acc 0.188 pos 0.9993 worst-neg 0.9993 min-margin -0.0001 eff-rank  1.00
  400 loss 2.3958 acc 0.125 pos 0.9999 worst-neg 0.9999 min-margin -0.0000 eff-rank  1.00
```

**Why collapse is a descent direction here.** At initialisation the cross-modal InfoNCE is
**3.0762, which is ABOVE chance** (ln 16 = 2.7726). The positives are systematically *worse* than the
negatives: mean positive cosine 0.1235 against a worst-negative of 0.2186, minimum margin −0.219.
The WSI biology states are 0.80 mutually collinear at init and the RNA states 0.62, so the
modality-specific mean dominates both views and the residual patient-specific signal is uncorrelated
across modalities — the cross-modal pairing carries no usable signal at step 0. From a
worse-than-chance start, **erasing every distinction is a genuine downhill move**: it takes the loss
from 3.0762 to chance. The optimiser takes it, and the destination is the permutation-symmetric
configuration where every `z` is the same vector — at which the InfoNCE gradient is *exactly zero* by
symmetry. It is a saddle, not a minimum, but it is an absorbing one for an encoder that has just
destroyed its own input sensitivity.

**The memory queue makes it worse, and that explains defect 1's partial success.** At step 0,
72–82% of the softmax mass sits on the 64 frozen queue keys, which are themselves 0.79 (WSI) and
0.56 (RNA) mutually collinear — effectively a *single* hard-negative direction. "Move away from one
shared direction" is a common-mode gradient, identical for all 16 patients, so it translates the
whole batch together. With the queue, collapse arrives by step 10; without it, by step 100.

**Neither declared anti-collapse term is positioned or weighted to stop this.**

1. `feature_decorrelation` is applied to `out_wsi["z_biology"]`, but **`variance_floor` is applied to
   `output["z_biology"]` — the FUSED view only** (`v2/training.py`, `if trains_biology:` branch). The
   two states that the D1 contrastive objective actually consumes, and that actually collapse, are
   the WSI-only and RNA-only biology states. The RNA biology view has no variance floor at all.
   The anti-collapse term and the collapse pressure act on different tensors.
2. Its weight is 0.01 with `target_std = d**-0.5 = 0.0625` at d=256, so the term contributes at most
   `0.01 × 0.0625 = 6.25e-4` against an InfoNCE of ~2.8 — three orders of magnitude too small to be a
   force. For reference VICReg, the method this pairing is borrowed from, uses variance weight equal
   to the invariance weight and 25× the covariance weight; here variance is 1/100 of the contrastive
   term and 1/4 of the covariance term. **The variance:covariance ratio is inverted by ~100×.**

**Consequences for two claims already in the notebook.** (a) The "decorrelation-driven collapse"
diagnosis is superseded — decorrelation aggravates but does not cause. (b) The reference number
"with decorrelation deleted, 800 steps at lr 1e-3 gave 2.0789 / retrieval 0.188 / patient cos 0.4946",
treated as a healthy comparator, is **itself a collapsed state**: per-row loss is exactly 2.08 for all
16 rows and 2.0794 = ln 8, i.e. the batch fell into two mutually orthogonal clusters of 8 with
within-cluster cosine 1.0. Retrieval 0.125 = 1/8 is chance *within a cluster*. Nothing measured so far
on this gate has been non-collapsed.

**Not a data problem.** The fixed batch is clean: 16 unique patient ids, 20–30 valid patches each,
one slide each, all RNA present, no duplicate RNA rows (max off-diagonal RNA cosine 0.9479), raw
patch-mean off-diagonal cosine 0.3265. The information is present; the encoder destroys it.

**Arithmetic of the criterion.** With 15 in-batch negatives at temperature 0.07, InfoNCE ≤ 0.10
requires a cosine margin of 0.3472 between the positive and every negative. That is a modest ask for
16 memorised patients and is not the blocker.

### In plain terms

The model was blamed for collapsing because of a regulariser. It is not the regulariser. Strip the
objective down to nothing but the contrastive loss and it still collapses — faster, in fact, than the
gate's own trajectory suggested. The reason is that the model starts out *worse than guessing*: at
step zero, a patient's tissue image looks slightly less like its own RNA profile than like a random
other patient's, because both representations are dominated by a large "average patient" direction.
When you are worse than guessing, one of the easy ways to improve is to stop making any distinctions
at all and just guess — and that is exactly what the optimiser discovers within fifty steps. Once
every patient is represented by literally the same vector, the loss has no idea which way to push, so
the model sits there forever.

The two safeguards that were supposed to prevent this are pointed at the wrong thing and turned down
too low. The variance floor is applied to the fused image+RNA state, while the states that actually
collapse are the image-only and RNA-only states. And its weight makes it 1000× smaller than the force
it is meant to oppose.

### Meaning for the claim

G2.6 has been failing for a real and previously unidentified reason, and the two earlier fixes
addressed real but secondary defects. The gate's `<= 0.10` threshold is not implicated: the required
cosine margin is 0.35 and the inputs are separable. No goalpost should move. The open question is
whether an anti-collapse term that is correctly *positioned* (on the per-view biology states the
contrastive loss consumes) and correctly *weighted* is sufficient, or whether the worse-than-chance
initialisation has to be addressed directly. That sweep is running.

Also note for D1's eventual interpretation: **if the biology head collapses to rank 1 under
`programme_free`'s own objective at gate scale, the same pressure exists in full training**, and any
D1 rank comparison between arms must be read with that in mind.

### Files / commits

- `~/ws_d1/diag_a.py`, `diag_b.py`, `diag_c.py`, `diag_d.py`, `diag_e.py` (A100)
- `~/e0_run/d1_diag/diag_{a,b,c,d}.log`, `diag_e_*.log`
- `v2/training.py` — `variance_floor` applied to `output["z_biology"]` (fused) not the per-view states
- `v2/losses.py` — `variance_floor`, `feature_decorrelation`
- `v2/runner.py` — `_overfit_programme_free_contrastive`, `_require_programme_free_overfit`
