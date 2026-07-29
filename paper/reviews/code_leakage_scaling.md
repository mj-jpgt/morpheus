# Code review: LEAKAGE and SCALING fixes (F-R2 decorrelation, F-R4 honest metrics)

Scope: `v2/training.py`, `v2/losses.py`, `v2/runner.py` (plus `v2/data.py`,
`v2/honest_metrics.py` for the two questions asked). Protocol assumed:
leakage-controlled 11-dev / 21-heldout, uncapped H-Optimus patches, A100-40GB.

VERDICT: code lane -- severity of worst finding = HIGH (the F-R2 decorrelation
fix is silently inert or degenerate on real uncapped batches; the reported
synthetic rank recovery 10.3->21.2 will not reproduce on the A100 real-data run
as currently wired).

---

## Q1. Does `feature_decorrelation` / `honest_metrics` compute any statistic across the held-out split or leak test info?

**No leakage found in either function. Both are batch-local / group-local.**

### `feature_decorrelation` (losses.py:46-70)
- Operates only on the current batch tensor `state` (z_biology). Mean/var are
  computed over `dim=0` of that batch. No dataset-, split-, or test-level
  statistic is referenced.
- In training it fires on **train** batches only (gradient path,
  training.py:250). In `evaluate_epoch` it also runs on the **val** loader, but
  purely as a logged metric under `@torch.no_grad()` (training.py:341-351) --
  no backward, no cross-batch aggregation, and val != test. The held-out test
  split is never loaded anywhere in `runner.run` (only `train`/`val` indices are
  materialised into `UncappedHoptimusBatches`, runner.py:426-432). Test is
  untouched by the trainer.

### `honest_metrics` (honest_metrics.py)
- `macro_group_pearson` / `control_adjusted_specificity` compute *within-group*
  (per-cancer) Pearson on whatever prediction/target/groups arrays the evaluator
  passes. That is the intended held-out reporting computation, not leakage: no
  train statistic crosses into it and no held-out statistic crosses back into
  training.
- Important robustness note (in the fix's favour): within-group Pearson is
  invariant to per-group affine transforms. So even if the evaluator
  residualises/standardises targets using per-cancer statistics estimated on the
  held-out set, the honest number is numerically unaffected -- a per-cancer mean
  subtraction and std division cancel in the correlation. The *pooled*
  cross-cancer Pearson is not affine-invariant, but that is exactly the
  confounded number F-R4 is de-emphasising. Good.

### Upstream leakage controls (verified clean)
- `residualise_programmes` (runner.py:33-40) and `attach_v2_targets`
  (runner.py:167-188) fit programme mean/scale on the `train` (or, under
  `--fit-development`, train+val) mask only; held-out cancers fall back to the
  **train** global mean/scale -- correct.
- Positive graph `_v2_positive` / `_v2_neighbour_indices` is built from
  `train_rows` only; positives can only reference train patients (runner.py:181-187).
- `TrainOnlyStandardizer`, `_standardize_clinical`, `_attach_numeric_table` all
  fit on `train & present` only.
- `_validation_selection` (runner.py:209-227) runs on **val** only.

**Minor / non-leakage flags**
- Protocol says 21 held-out cancers but `--expected-heldout-cancers` defaults to
  **22** (runner.py:486) and the Lambda scripts pass it through. If the true
  protocol is 11/21, the run must explicitly pass `--expected-heldout-cancers 21`
  or `validate_runtime_split` will assert against the wrong count (either a hard
  preflight failure, or -- worse if the guard is permissive -- a silent
  off-by-one in the split contract). Confirm the flag matches the paper's 21.

---

## Q2. Does the min-batch guard interact badly with dynamic token batching (small ragged B)? -- YES (HIGH)

`feature_decorrelation` skips the batch when `len(state) < min_batch` (=8)
(losses.py:64). "len(state)" here is **B = number of patients in the batch**, not
patches. The dynamic sampler makes B systematically tiny.

### Why B is small and often < 8
`DynamicTokenBatchSampler.batches` (data.py:60-80) greedily packs whole,
**uncapped** patient bags until `used + count > token_budget`. So
`B ~= floor(token_budget / mean_tokens_per_patient)`.
- `token_budget` = 32768 (runner default), 16384 in `run_v22_a10`, and
  `run_v21_recovery` OOM-retries by *lowering* the budget (min 16384).
- Uncapped H-Optimus (224px tiles) on TCGA slides routinely yields ~5k-50k
  tokens/patient. At budget 16384-32768 that gives **B ~= 1-3 per batch** for
  typical/large slides, and B=1 (solo batch) for any patient above the budget.
- Decorrelation needs B>=8, i.e. mean tokens/patient <= ~2048-4096. That is the
  small-WSI tail only.

**Consequence 1 (correctness of the fix): F-R2 is largely inert on real data.**
The Barlow off-diagonal term is the central LEAKAGE/collapse fix (rank 10.3->21.2
on synthetic *dense* batches). On the uncapped A100 run it will return exactly
0.0 for most/all batches -- no gradient -- so the biology head is free to collapse
back to eff-rank ~5-6 exactly as the paper diagnoses. The synthetic validation
used batch sizes that never occur in the real loader. This is the worst finding:
the headline fix may not actually be exercised in the queued real-data ablation.

**Consequence 2 (scaling coupling): the memory knob silently disables the fix.**
OOM handling lowers `token_budget`, which lowers B, which pushes more batches
below `min_batch` -- so the standard A100-40GB response to OOM directly turns off
F-R2 with no warning. The decorrelation-weight ablation arm
(`--decorrelation-weight 0.0`) and the "on" arm may be indistinguishable on real
data because the "on" arm rarely fires.

**Consequence 3 (biased subsample): when it does fire, it fires on the wrong
patients.** `batches()` sorts by token_count and shuffles only within
contiguous 64-wide chunks (data.py:61-69), so batches are size-homogeneous.
B>=8 is reached only in the small-token strata. Decorrelation is therefore
applied almost exclusively to patients with few patches (smaller sections/tumors)
-- a cohort-correlated confound, ironic for a loss whose purpose is
de-confounding.

**Consequence 4 (degenerate estimator even when active).** z_biology is 256-D.
The batch correlation is a rank-(B-1) estimator (losses.py:68). At B=8 you
estimate a 256x256 correlation from 7 effective samples; the off-diagonal is
dominated by sampling noise, so the penalty injects noise gradients rather than a
reliable decorrelation signal. Barlow-Twins-style penalties assume N well above
the feature dimension. `min_batch=8` is far too permissive for dim=256.

**Consequence 5 (same pathology hits the other batch-guarded losses).**
`symmetric_infonce` returns 0 for B<2 (losses.py:14) and has only B-1 negatives;
solo/2-patient batches give the identity contrastive little or no signal --
and retrieval R@10 is the selection gate (runner.py:439-445), so selection is
driven by whichever epochs happened to pack larger batches.
`programme_neighbourhood_loss`/`supervised_programme_contrastive` (B<3) similarly
go quiet. Note the authors already built `ProgrammeMemoryBank` (training.py:77)
to rescue supcon cross-batch positives from exactly this ragged-batch starvation
-- but built no analogous rescue for decorrelation or identity InfoNCE.

**Consequence 6 (silent).** The skip returns zeros with no metric. The logged
`decorrelation` value is just averaged over rows that recorded it; a batch that
skipped still runs `step`, and `metrics["decorrelation"]` is only set inside the
`if weights["decorrelation"]:` block (training.py:247-253) which *does* execute
(the guard is inside the loss fn), so the logged mean is dominated by 0.0s with
no visibility into how many batches were genuinely active.

### Mitigations (in priority order)
1. **Decouple decorrelation from ragged B via a feature bank / accumulation.**
   Maintain a detached running buffer of z_biology (mirror `ProgrammeMemoryBank`,
   or an EMA feature covariance) and compute the Barlow off-diagonal over the
   accumulated N each step, so the estimate is stable and token_budget-invariant.
   Alternatively accumulate the correlation across micro-batches within one
   optimizer step. This is the fix that actually makes F-R2 real on A100.
2. **Guarantee a patient floor in the sampler.** Add a min-patients-per-batch
   target: for the decorrelation view, cap tokens-per-patient (subsample patches
   for *this loss only*, keeping the uncapped bag for identity/programme) so B
   reliably reaches a Barlow-adequate size; or split very large patients across
   the token budget while grouping enough distinct patients.
3. **Raise/justify `min_batch` relative to dim, and shrink the projector.**
   256-D from B~8 is degenerate. Either compute decorrelation on a lower-dim
   projection of z_biology, or require N >> dim via (1). Make `min_batch` a
   config surfaced in the manifest, not a buried default.
4. **Instrument, don't silence.** Emit `decorrelation_active_fraction` (share of
   train batches with B>=min_batch) and mean-B per epoch. If active_fraction is
   near 0 on the real run, the ablation is invalid and must be re-run with (1)/(2).
   Currently there is no way to tell from the logs.
5. **Mix sizes across batches.** Replace within-64-chunk shuffling with
   size-stratified or token-balanced batching so decorrelation isn't confined to
   the small-WSI stratum.
6. **Pin the ablation's token_budget and disable OOM-driven budget lowering for
   the F-R2 arm**, or the on/off arms are confounded by how often the loss fired.
7. **Fix the 21-vs-22 held-out flag** before the real-data run so preflight
   asserts the intended contract.

---

## Bottom line for the paper
The leakage controls are sound: neither `feature_decorrelation` nor
`honest_metrics` touches the held-out split, and the within-cancer honest metric
is even affine-robust to eval-time residualisation. The SCALING side is where the
risk concentrates: as wired, the F-R2 decorrelation fix is effectively disabled
or degenerate on real uncapped A100 batches, so the queued real-data ablation
could show "no rank recovery" for an implementation reason, not a scientific one.
Land mitigations (1)+(4) at minimum before trusting the A100 T4-ablation result.
