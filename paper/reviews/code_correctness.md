# Code Correctness Review — F-R1 / F-R2 / F-R4 fix diffs

Adversarial review of the collapse-fix diffs: `feature_decorrelation`
(`v2/losses.py`), biology normalize F-R1 (`v2/model.py`), decorrelation wiring +
profile gating (`v2/training.py`), `--decorrelation-weight` (`v2/runner.py`),
`v2/honest_metrics.py`, and the two test files. Findings ranked by severity.

---

## F1 (HIGH) — `min_batch=8` can silently no-op the entire F-R2 fix on real data

`feature_decorrelation(state, min_batch=8, ...)` returns `state.new_zeros(())`
whenever the batch has fewer than 8 patients (`losses.py:64-65`). In
`training.py:247-250` the guarded call is `if weights["decorrelation"]:
feature_decorrelation(output["z_biology"])` — no floor is enforced upstream.

The batch size in *patients* is set by `DynamicTokenBatchSampler`
(`data.py:70-80`): it packs complete, uncapped patient bags until
`token_budget` is exceeded. With the default `--token-budget 32768` and
uncapped H-Optimus tile bags (TCGA WSIs routinely yield thousands to tens of
thousands of tiles per patient), the per-batch patient count is frequently well
below 8 — a single high-token patient can even form a batch of 1
(`data.py:74`, only splits `if batch and used+count > budget`).

Consequences:
- On the real-data T4 ablation, the decorrelation term (the whole F-R2
  mechanism) may fire on only a minority of batches, or effectively never,
  depending on the tile-count distribution. The measured "fix arm" could be a
  near-no-op vs the `--decorrelation-weight 0.0` baseline, quietly nulling the
  headline result.
- The collapse-*inducing* structural losses have far lower floors:
  `programme_neighbourhood_loss` needs `len>=3`, `supervised_programme_contrastive`
  needs `len>=3`, Gaussian NLL has no floor. So on 3–7-patient batches the
  low-rank programme-manifold pull is applied while the anti-collapse
  decorrelation is skipped — the exact regime the fix targets is the one where
  it is silent.
- No test exercises the boundary: every `_big_batch` fixture uses `n=16`
  (`test_stress_collapse.py:45,138`), always above the guard, so the silent
  no-op path is uncovered.

Recommend: confirm the real batch-size distribution at
`token_budget=32768`; if batches are commonly `<8`, either lower `min_batch`
(a rank-`B-1` estimator is still informative at B=4–6), accumulate the
correlation across micro-batches, or raise the token budget. At minimum, log
the fraction of batches where decorrelation was skipped so a silent no-op is
visible in `train_metrics.jsonl`.

## F2 (HIGH/MEDIUM) — decorrelation is applied to `full_biology` only; the view carrying the collapse pressure (`wsi_biology`) is never decorrelated

`training.py:250` calls `feature_decorrelation(output["z_biology"])` where
`output` is the `view="full"` forward. It is never applied to `out_wsi` or
`out_rna`. But the programme loop (`training.py:222-226`) gives the *structural*
collapse-inducing losses (neighbourhood + supcon) to the WSI view in the `full`
profile:
- `full` profile: `wsi_biology` gets NLL + neighbourhood + supcon (collapse
  pressure) but **no decorrelation**; `full_biology` gets NLL +
  decorrelation but **no structure**.
- `programme_only` profile: all three views get structure
  (`structure_for_all_biology_views=True`, `training.py:221`), yet decorrelation
  still acts on `full` only — so `wsi_biology` and `rna_biology` receive
  structural pull with no counter-pressure.

WSI→molecular prompting is the primary reported path (see the comment at
`training.py:217-219` and `_validation_selection` using `out_wsi`
programme/identity in `runner.py:218`). The exported biology fingerprint states
include `wsi_biology`/`rna_biology`/`full_biology` (`runner.py:315`,
`_trained_states_for_profile`). If the reported eff-rank fingerprint is measured
on `wsi_biology`, the fix acts on a *different* state than the one measured, and
only benefits it indirectly through the shared `self.biology` weights and query
stack — not through matched batch statistics. The synthetic 10.3→21.2 claim was
presumably measured on the decorrelated view; that improvement need not transfer
to `wsi_biology`.

Recommend: apply `feature_decorrelation` to whichever biology view(s) the paper
reports rank on (most likely `out_wsi["z_biology"]`), or state explicitly that
the fingerprint is `full_biology`. As written the anti-collapse term and the
collapse-inducing term are on different views.

## F3 (MEDIUM) — the F-R2 collapse tests do not actually test that the fix works

`test_biology_head_effective_rank_is_not_degenerate` asserts
`effective_rank(out["z_biology"]) > 1.5` (`test_stress_collapse.py:145`). A
threshold of 1.5 only rules out near-total collapse; the fix claim is
10.3→21.2. The test's own comment concedes it "documents the instrument". There
is:
- no A/B assertion comparing `--decorrelation-weight 0.0` vs `0.04` (or
  `weights["decorrelation"]` on/off) showing the term raises rank;
- no assertion that `metrics["decorrelation"]` is present and non-zero, so a
  silent no-op (F1) or a mis-wired weight would still pass;
- `test_full_training_step_metrics_are_finite` only checks finiteness, and
  `test_every_active_head_receives_gradient` checks `biology` gets *some*
  gradient — but biology also receives programme/neighbourhood/supcon gradient,
  so decorrelation being dead is invisible to it.

Net: no test in the suite would fail if F-R2 were accidentally disabled. Add a
test that runs the same fixture with decorrelation weight 0 vs 0.04 and asserts
the trained eff-rank is materially higher with it on (this also directly guards
F1 and F2).

## F4 (LOW-MEDIUM) — `variance_floor(z_biology)` is near-unsatisfiable and partly fights F-R1 normalization

`variance_floor` uses `target_std=1.0` (`losses.py:29`). After F-R1,
`z_biology` is per-row L2-normalized in 256-D (`model.py:289`), so per-feature
std across the batch is bounded near `1/sqrt(256) ≈ 0.06` and can never approach
1.0. `variance_floor(output["z_biology"])` is therefore a nearly constant ~0.94
penalty whose gradient pushes to concentrate per-feature variance — mildly
opposing the unit-norm geometry the fix relies on. Weight is only 0.01
(`separation_after_warmup`/`variance_after_warmup`) and it is zeroed entirely in
`programme_only`, so impact is small, but the target_std is mis-scaled for a
unit-normed state. Consider a target consistent with the normalized geometry
(e.g. `~1/sqrt(dim)`) or dropping the floor on `z_biology` now that
decorrelation is the rank tool.

## F5 (LOW) — biased/unbiased mismatch inside `feature_decorrelation`

Standardization divides by `sqrt(var(unbiased=False))` (N) but the correlation
matrix divides by `len(state)-1` (N-1) (`losses.py:66-68`). Off-diagonals are
thus scaled by `N/(N-1)` relative to true Pearson correlations (≈+7% at N=8,
+2% at N=16). Harmless to the objective (a constant reweighting absorbed by the
loss weight) but the diagonal is not exactly 1 and the code is not a textbook
Pearson matrix; worth a one-line note if the term's magnitude is quoted.

## F6 (LOW) — `honest_metrics` edge paths lightly tested

`control_adjusted_specificity` correctly NaN-propagates when either the real or
control macro is NaN (`honest_metrics.py:77-78`), but no test covers that path
(all-small-groups → NaN) for the *adjusted* function; only
`macro_group_pearson`'s NaN path is tested. `test_honest_metrics.py:56-57`
uses a loose bound (`abs(null_spec) < 0.4`) that would tolerate a fairly large
spurious specificity. Logic itself is correct: `_pearson` guards `len<2` and
`denom<=0` (`honest_metrics.py:28-34`), `macro_group_pearson` drops sub-`min_group`
and NaN groups (`:55-59`). No correctness bug found here.

---

## Wiring items verified correct
- `--decorrelation-weight` default 0.04 → `V2LossSchedule(decorrelation_after_warmup=...)`
  → `weights["decorrelation"]` → guarded add (`runner.py:380`, `training.py:59,247`).
  Setting `0.0` cleanly disables the term (`if weights["decorrelation"]:` false),
  so the ablation baseline arm is correctly wired.
- `programme_only` keeps `decorrelation` in its active-key set
  (`training.py:69-73`); `identity_only` drops it. Correct per intent.
- Decorrelation is 0 during warmup (`epoch < warmup_epochs`, `training.py:51,59`);
  fine given epochs=40/warmup=4 defaults.
- F-R1: programme heads now consume the normalized `biology_state`
  (`model.py:289,296`); gradient path programme→`biology_state`→`self.biology`
  intact.

---

### VERDICT
- severity of worst finding: **HIGH** — F1 (guard may silently no-op F-R2 on
  real uncapped batches) and F2 (decorrelation acts on `full_biology` while the
  collapse pressure and likely-reported fingerprint sit on `wsi_biology`), both
  compounded by F3 (no test would catch either). These directly threaten the
  validity of the queued real-data T4 rank ablation.
- Recommended gate before trusting T4 numbers: log per-batch patient count +
  decorrelation-skip fraction, confirm the term fires, and align the
  decorrelated view with the measured fingerprint view.
