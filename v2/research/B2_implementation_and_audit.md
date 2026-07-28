# Branch B — Implementation & Audit Synthesis

Master synthesis over the 7 B1 ledgers and the adversarial critiques. Scope: fix the observed
biology-head rank collapse (effective rank ~5-6 of 256; identity ~84) and the confounded
molecular-prompting metric, under the leakage-controlled 11-dev/21-heldout protocol, uncapped
H-Optimus patches, Lambda A10+.

The single most important synthesis finding: **there are three competing root-cause theories for
the rank collapse, and the critics disagree on which is primary.** They are NOT mutually exclusive
and the cheapest correct plan resolves the disagreement *empirically* with an effective-rank test,
not by fiat. The three candidate drivers:

- **T1 (Lane 7 / F6-A):** the programme Gaussian-NLL (weight 1.0) regresses the *unnormalized*
  `biology_state` while every rank-spreader acts on the *normalized* `z_biology` — one strong
  low-rank pull on a different geometry than the four weak spreaders.
- **T2 (Lane 2 / F2):** `programme_neighbourhood_loss` KL-distills `z_biology`'s Gram matrix onto
  the ~50-D (effective rank ~5-6) Hallmark target Gram every step, at weight 0.20 — a distillation
  collapse that re-imposes low rank regardless of geometry.
- **T3 (Lane 1 / F1-covariance):** there is simply no cross-dimension decorrelation term; the
  per-dim `variance_floor` provably cannot raise rank.

**The gating test that adjudicates all three is the same: an effective-rank assertion on
`z_biology`. It does not exist today. Build it first.**

---

## 0. Cross-cutting facts every fix must respect (verified in ledgers + critiques)

- **No test catches the failure.** The 31 CPU tests check per-head gradient *existence* and
  `isfinite` only. Nothing catches rank collapse, NaN in KL/NLL, `log_variance`->clamp mode
  collapse, no-op loss, or anchor-gate saturation. Any fix that ships without a rank/finiteness
  guard can silently no-op and pass CI. This is the load-bearing gap.
- **Profile gates zero terms.** `weights()` in `programme_only` zeros `separation`, `variance`,
  and anything gated with them (training.py:63-65). A new anti-collapse term MUST be added to the
  `full` AND `programme_only` gate paths or it is dead exactly where biology geometry is trained.
- **Small, ragged batches.** DynamicTokenBatchSampler + uncapped patches make patients-per-batch
  `B` small and variable (can be 2-3). Any batch-covariance term is a rank-`(B-1)` estimator and
  needs a min-`B` guard; any whitening/std normalization needs a NaN guard.
- **The metric fix is orthogonal to the model fix.** The within-cancer honest metric (F6-B) can and
  should land independently; it changes reporting/gating, not training. Expect it to *lower* the
  headline (honest ~+0.07, baseline-matched) — that is the correct outcome, not a failure.

---

## 1. RANKED, vetted fix list (payoff × minimality)

Ranking key: payoff = probability × magnitude of moving biology effective rank (or metric honesty);
minimality = lines changed × blast radius (checkpoint/export/protocol). Every fix names an exact
`file:line`, the concrete change, the FIT/FAIL/SCALING verdicts distilled from the critiques, and the
ONE gating stress test.

### R0 (PREREQUISITE, land before any model fix) — Effective-rank + finiteness guard test
- **What / where:** New file `morpheus/v2/tests/test_stress_collapse.py`. Add `effective_rank()`
  (Roy-Vetterli: `exp(entropy of L1-normalized svdvals)`) and the six stress tests specified in
  §3. Reuse `_batch()` (tests/test_v21_model.py:19) and `_model_and_batch()`
  (morpheus/tests/test_v2_recovery_contracts.py:20).
- **Verdicts:** FIT good (Lane 6, all critics converge that this is the missing instrument).
  SCALING trivial (<5 s CPU). Payoff: this is what makes every downstream fix *attributable and
  non-silent*.
- **Gating test:** itself — it must FAIL on a synthetic rank-deficient `z_biology` and PASS on a
  well-spread one before it is trusted to gate anything.

### R1 (TRY FIRST) — F6-A: NLL heads consume the normalized biology state
- **What / where:** model.py:283 (and the dict at model.py:283-290). Change
  `biology_state = self.biology(biology)` → `biology_state = nn.functional.normalize(self.biology(biology), dim=-1)`,
  and let `z_biology`, `programme_mean`, `programme_log_variance` all read that single normalized
  tensor. ~2 lines, zero added FLOPs (normalize already computed).
- **Verdicts:**
  - FIT: contested. Lane 7 / one F6 critic = **good_fit, primary driver (T1)**. Two F6 critics =
    **bad_fit / no-op**: they argue `self.biology` already gets NLL gradient (test_v21_model.py:108
    passes), so this connects no *new* gradient and does not by itself raise rank; and that
    normalizing a magnitude regressor (Gaussian-NLL) away from its radial DOF risks pushing
    `log_variance` toward the +8 clamp (mode collapse).
  - FAIL modes: (a) no-op on rank if T2 (neighbour-KL) is the real driver — the KL still pins the
    Gram to rank ~5-6 after the change; (b) `log_variance` saturation / NLL miscalibration because
    unit-norm input bounds `(target-mean)^2`; (c) checkpoint non-equivalence — the head must re-fit
    scale, so retrain from warmup, do not hot-swap into a late checkpoint.
  - SCALING: unanimous **scales** — O(1), reuses one tensor, no export/protocol change.
- **Gating test:** `test_biology_effective_rank_does_not_collapse_below_identity` (§3a) must rise
  after N steps, AND `test_programme_log_variance_not_saturated_and_mean_not_constant` (§3c) —
  `log_variance.max() < 7.9`. Both must pass. If rank stays low, T1 is falsified → proceed to R2.
- **Why first:** smallest diff, zero blast radius, and it is the clean *experiment* that
  adjudicates T1 vs T2. Its result (rank moves or not) tells you whether to stop or escalate to R2.

### R2 (PRIMARY MECHANISM FIX if R1 under-delivers) — F1-covariance / Lane-1 feature decorrelation
- **What / where:** losses.py — add `feature_decorrelation(state)` = mean-squared off-diagonal of
  the centered batch covariance (Lane 1 code, losses.py:61-76 spec, mirrors
  `centered_cross_covariance` at losses.py:21-26). training.py:32 — add
  `decorrelation_after_warmup: float = 0.04` to V2LossSchedule and include it in `weights()` for the
  **full AND programme_only** gate paths (training.py:63-65). training.py:232-237 — after the
  variance-floor block, add the weighted term on the **normalized** `z_biology`. Guard: skip when
  `B < 8`.
- **Verdicts:**
  - FIT: **good_fit** across critics as the mechanism-matched anti-rank-collapse term (the missing
    off-diagonal decorrelation; `variance_floor` is per-dim and cannot raise rank). Leakage-safe by
    construction (within-batch statistic, fits no transform).
  - FAIL modes: (a) can be dominated/no-op if left too weak against the ~0.40 combined
    neighbour-KL+supcon collapse pressure — one critic argues the KL driver simply wins unless it is
    also reduced (see R3); (b) rank-deficient at small `B` (<<256): decorrelates only `B-1`
    directions/step; (c) whitening/std eps blow-up at near-constant dims → needs the min-`B` guard +
    a NaN assert.
  - SCALING: **scales** (O(B·256²), negligible vs patch attention) but **batch-composition
    sensitive** — one critic downgrades to "concern" on small ragged `B`. Mitigate with the min-`B`
    guard and a modest weight (0.04, do not exceed 0.1).
- **Gating test:** §3a effective-rank ratio gate must clear `>= 0.4 × rank(z_identity)` after N
  steps AND §3b finiteness (`test_kl_and_nll_are_finite_under_degenerate_inputs`) must stay finite
  with the new term active.

### R3 (ESCALATION if R1+R2 still collapse) — F2 (Lane 2): replace neighbour-KL with RNA-paired biology InfoNCE
- **What / where:** training.py:157-160 — replace the `programme_neighbourhood_loss(...)` call with a
  new `biology_paired_contrastive(z_biology_wsi, z_biology_rna, positive_mask)` reusing
  `symmetric_infonce` (losses.py:13-18) semantics; RNA view `out_rna` is already computed
  (training.py:183), positives = existing top-8 programme mask (runner.py:141). Keep the
  `neighbourhood_after_warmup=0.20` weight, retarget its role. Keep Gaussian-NLL unchanged (Lane 2
  argues it is benign). Bundle with the R2 covariance term for batches with sparse positives.
- **Verdicts:**
  - FIT: **good_fit for the diagnosed T2 driver** — gives biology the same rank-expanding paired
    InfoNCE identity already has (why identity is rank ~84). No new forward pass, no export/protocol
    change.
  - FAIL modes: (a) InfoNCE rank depends on in-batch negatives → small ragged `B` weakens it; the
    detached memory bank (training.py:104) can supply *negatives only*, positives must stay in-batch;
    (b) batch composition risk — an all-same-cancer batch makes the paired positive trivial and can
    re-collapse; (c) **do NOT keep any nonzero neighbour-KL weight** — Lane 2 shows any nonzero weight
    re-imposes the low-rank Gram.
  - SCALING: **scales** (reuses out_rna, O(256²) extra) but shares the small-`B` sensitivity of all
    contrastive terms.
- **Gating test:** §3a effective-rank ratio gate + a NaN guard on the new contrastive when a batch
  has zero positives (mirror losses.py:48,57 guards). Larger change than R1/R2 → land last.

### R4 (INDEPENDENT, land in parallel) — F6-B: within-cancer honest metric as primary/gated
- **What / where:** paired_bootstrap.py (~line 14) — add `macro_group_metric(base, groups)` closure.
  v21_evaluation.py:553 — pass macro metric to the paired-vs-teacher bootstrap. v21_evaluation.py:439
  — difference biology-vs-matched-control on `macro_cancer_pearson`. select_v21_profile.py:46-53 —
  change the promotion-gate filter from `metric == "pearson"` to `metric == "macro_cancer_pearson"`.
  Report BOTH pooled and macro everywhere; never drop pooled.
- **Verdicts:** FIT **good_fit** (machinery already exists: `_macro_cancer_pearson`
  v21_evaluation.py:188, `build_matched_random_controls` discovery_targets.py:178, cancer-cluster
  bootstrap). SCALING **scales** (numpy/pandas, ~1 unit test). Grounded in Howard 2021 / HESCAPE.
  FAIL/expectation: the gate gets harder and the headline drops to the honest ~+0.07 baseline-matched
  number — intended, per Lane 5.
- **Gating test:** unit test that `macro_group_metric` over 2 groups = mean of within-group pearsons
  and returns NaN when every group has <3 pairs (§3, honest-metric guard).

### R5 (DIAGNOSTIC ONLY, optional) — F5 measurement half: fix the gradient-conflict metric
- **What / where:** training.py:273-277 — replace the shared-param set (which wrongly uses pre-trunk
  `wsi.patch[-1]`/`rna.projection[-1]`) with the real post-attention trunk
  `self.model.queries + self.model.blocks.parameters() + self.model.norm.parameters()`.
  training.py:255-261 — un-collapse `_last_loss_components` into `programme_nll`,
  `programme_neighbour`, `programme_supcon` so intra-biology conflict is visible. Run ONE epoch at
  `--gradient-diagnostics-every 1`.
- **Verdicts:** FIT **good_fit but "measure, don't add"** — the logging already exists
  (training.py:264-305); the task is to *fix* it, not add it. SCALING **concern**: `retain_graph=True`
  raises peak memory on diagnostic steps under uncapped patches → possible every-25th-step OOM on
  A10; keep the cadence coarse (e.g. 100) and compute conflict before `loss.backward()`.
- **Gating test:** a smoke test that the diagnostic path runs finite at `gradient_diagnostics_every=1`
  on `_batch()`. **PCGrad is NOT approved** (see §5).

---

## 2. Dead-wiring audit table

Column key: **Verdict** = bug vs intentional (per Lane 7 + critics). **Minimal safe change** = the
smallest change that is safe; "none" means leave as-is.

| # | Wiring item | file:line | Trained? | Exported? | Verdict | Minimal safe change |
|---|---|---|---|---|---|---|
| 1 | `z_context` (Linear 512→128, NOT normalized) | model.py:287; export.py:54-55 (absent) | No | No | **Intentional-but-wasteful.** Consumes cross-attn capacity for nothing. | Leave wired-off; freeze the contract with a test (§3d `test_dead_heads_are_not_silently_exported`). Reclaim slots only in a separate capacity PR. |
| 2 | `z_uncertainty` (2 slots) | model.py:291; export absent | No | No | **Intentional-but-wasteful.** | Same as #1 — freeze non-export, do not resurrect now. |
| 3 | 20 residual slots (5×4: wsi/rna/clinical/snv/cnv) | model.py:257-259, 291 | No | No | **Intentional-but-wasteful.** 75% of the 32-slot bank feeds dead heads; residuals feed `z_patient` via `query.mean(1)` (model.py:288). | Reclaiming requires bumping `query_slots` (model.py:44) → changes `self.queries` shape + manifest → **breaks checkpoints/export**; retrain-from-scratch. **Do NOT bundle with a collapse fix.** Keep as separate hygiene PR. |
| 4 | mean-pool over 4 slots (`take()`) | model.py:254 | n/a | n/a | **Intentional, minor — NOT the collapse cause.** Decisive counterexample: identity uses the *identical* 4-slot mean-pool yet is rank ~84. | **None.** Do not attention-pool/concat as a rank fix (F4 sunk — §5). Concat also breaks `self.biology` width (512→2048) → state_dict break. |
| 5 | Normalize inconsistency: NLL reads unnormalized `biology_state`; spreaders read normalized `z_biology` | model.py:283,286,290; losses.py:46-66 | Yes | Yes (z_biology) | **BUG-ish, contested (T1).** Candidate primary collapse driver per Lane 7; disputed as no-op by two F6 critics. | R1: normalize once, feed everywhere (model.py:283). ~2 lines. Ship with §3a rank + §3c log_variance guards; retrain from warmup. |
| 6 | Detached memory bank (`states.detach()`) | training.py:104; read side training.py:95 | Current-batch anchors carry gradient; bank entries do not | via z_biology | **CORRECT as designed.** Standard MoCo queue; back-prop into stale activations is impossible/undesirable. Not a bug. | **None.** May be reused as *negatives* for R3's contrastive; positives must stay in-batch/attached. |
| 7 | Coordinate path gated off (`coordinate_present_fraction=0`) | model.py:117,132-139; runner.py:43-50 | No | No | **Intentional leakage hygiene** (guards against manufactured zero-coords). | **None.** Do not force on; not collapse-relevant. |
| 8 | Semantic / PLIP head (`semantic_dim=0` default) | model.py:293-295 | No | No | **Intentional, capability-gated** (needs PLIP cache). | **None.** |

---

## 3. Fast stress-test suite spec

New file `morpheus/v2/tests/test_stress_collapse.py`. Shared helpers `effective_rank()` and
`_train_n()` per Lane 6 (svdvals in fp32; `gradient_diagnostics_every=0`; seed each test). Runtime
target <5 s CPU. Small fixtures (batch=3, hidden=16) cannot show absolute rank-256 floors, so assert
**relative/mechanistic** properties on a locally-built `B>=8` k-cluster **separable** synthetic batch.

| Test name | Assertion | Fixture | Fails on |
|---|---|---|---|
| **(a)** `test_biology_effective_rank_does_not_collapse_below_identity` | `effective_rank(out["z_biology"]) >= 2.0` AND `>= 0.4 * effective_rank(out["z_identity"])` after 20 steps, profile `full` | `_batch()` w/ B≥8 separable clusters | **Rank collapse** (real failure ratio ≈ 5.5/84 ≈ 0.065 ≪ 0.4) |
| **(b)** `test_kl_and_nll_are_finite_under_degenerate_inputs` + `test_full_step_metrics_all_finite` | `isfinite(programme_neighbourhood_loss(state,targets))` with an all-zero row (÷~0 in `_norm`); `isfinite(gaussian_nll(mean, log_variance=40, target-with-masked-NaN))`; every value in one `V2Trainer.step` `metrics` finite | `_batch()` | **NaN** in the `-1e4`→softmax→kl_div chain (losses.py:51-52) and `exp(-log_variance)` overflow (losses.py:77) |
| **(c)** `test_programme_log_variance_not_saturated_and_mean_not_constant` | `out["programme_log_variance"].max() < 7.9` (not pinned at +8 clamp) AND `out["programme_mean"].std(0).mean() > 1e-3` after 20 steps, profile `programme_only` | `_batch()` | **Mode collapse** — NLL explaining everything as noise / constant mean |
| **(d)** `test_every_declared_active_head_receives_gradient` + `test_dead_heads_are_not_silently_exported` | per profile, `assert_module_used(active_head, loss)` (reuses contracts.py:76 helper) and grad None/zero for declared-inactive heads; `export.widths` contains none of `z_context/z_uncertainty/*_residual` | `_batch()`, `_model_and_batch()` | **No-op loss** / accidental export of an untrained (dead) head |
| **(e)** `test_anchor_gate_not_saturated_and_residual_can_move` | `out["anchor_gate_mean"] ∈ (0.02, 0.98)`; `out["anchor_correction_norm"] > 1e-4` after 20 steps; `anchor_residual_scale` moved off 0 OR gate variance > 1e-4 | `_model_and_batch()` (anchor=True), separable batch | **Anchor saturation** (observed residual ≈ 0; residual never learns to move) |
| **(f)** `test_take_consumes_every_query_slot_exactly_once` | `identity+biology+context+5*residual+uncertainty == config.query_slots`; instrumented `take` cursor ends at `query.shape[1]`; residual dict has exactly `{wsi,rna,clinical,snv,cnv}` | `_batch()` | Slot-accounting desync (guards any slot-count edit, e.g. R3 dead-slot reclamation) |
| **(g)** `test_macro_group_metric` (honest metric) | macro over 2 groups == mean of within-group pearsons; returns NaN when every group has <3 pairs; `_control_comparison_rows` emits a macro-based row | tests/test_data_evaluation.py synthetic arrays | Silent pooled-metric confound / NaN mishandling for R4 |

---

## 4. Recommended ORDER of operations, and the FIRST change

1. **R0 — build `test_stress_collapse.py` (effective-rank + finiteness + log_variance guards).**
   THIS IS THE FIRST CHANGE. Rationale: every candidate model fix is disputed as a possible no-op,
   and no existing test can tell success from silent failure. R0 turns each subsequent step into an
   attributable experiment and is zero-risk (<5 s CPU, no training change). Verify it FAILS on
   synthetic rank-deficient `z_biology` and PASSES on well-spread before proceeding.
2. **R1 — normalize the biology state (model.py:283).** Smallest diff, zero blast radius. Run R0.
   - If §3a rank rises and §3c log_variance stays < 7.9 → T1 confirmed, largely done; add R2 only if
     rank is still short of ~30-50.
   - If rank stays low → T1 falsified; escalate.
3. **R2 — add the feature-decorrelation term (losses.py + training.py:232-237), gated in full AND
   programme_only, min-B≥8 guard, weight 0.04.** Re-run R0. This is the mechanism-matched anti-rank
   term and should be the standing fix.
4. **R3 — only if R1+R2 still collapse: replace neighbour-KL with RNA-paired biology InfoNCE
   (training.py:157-160), neighbour-KL weight → 0.** Largest change, land last; keep R2 covariance as
   a floor for sparse-positive batches.
5. **R4 — land the honest within-cancer metric in parallel** (independent of the model track; only
   touches eval/gating). Expect and communicate a lower, baseline-matched headline.
6. **R5 (optional) — fix the gradient-conflict *measurement* and run one diagnostic epoch** to check
   whether the collapse is directional conflict (likely NOT) before ever considering surgery.

Run R1→R2→R3 as a *sequence of single-mechanism experiments*, each gated by R0, so the rank
movement is attributable to exactly one change.

---

## 5. Explicit "do NOT do" (fixes a critic sank)

- **Do NOT ship the F1 `separation` bump (0.01→0.1) as a rank fix.** `whitened_cross_covariance`
  (losses.py:37-43) decorrelates identity↔biology (BETWEEN blocks); a rank-6 biology can be perfectly
  decorrelated from identity. It cannot raise biology's *own* rank — category error (Lane 2:28,
  Lane 1:111). It is also zeroed in `programme_only`. Ship the intra-block decorrelation term (R2)
  instead. If you still want the bump, ramp it and ablate it separately behind a rank metric — never
  as a hard 10× jump.
- **Do NOT do F2's "256-D contrastive" sub-variant.** `z_biology` is ALREADY a 256-D L2-normalized
  head (model.py:203,286); widening a head that is already 256-D is a pure no-op. Only the RNA-paired
  InfoNCE (R3) or per-gene targets change anything.
- **Do NOT scope F2 per-gene targets as a `losses.py:46-66` edit.** It changes `programme_dim`
  (model.py:189, runner.py:360), the Gaussian-NLL head widths (model.py:211-212), the export widths
  dict + manifest (export.py:52,54), and re-fits the residualizer. Naive edit crashes the
  losses.py:80-81 shape check. Do it only as a full target-pipeline PR, and only if higher-rank
  *regression* is independently wanted — it is NOT needed for the mechanism.
- **Do NOT ship F3 SigLIP as a fix for the biology collapse.** SigLIP swaps the *identity* loss
  (training.py:194, z_identity only); biology gets no paired contrastive and is untouched. It also
  hits the same anchor bottleneck (residual ≈ 0 at init, 0.25-bounded) that already caps InfoNCE, so
  the expected delta is ~+0.005. SigLIP is a legitimate *identity* experiment (register `t`,`b` as
  `nn.Parameter` on the model, use `.mean()` reduction, add a grad-finite test) but must not be sold
  as the rank remedy.
- **Do NOT do F4 (attention-pool/concat the biology slots; reclaim 20 slots) as a rank fix.**
  Decisive: identity and biology share the identical 4-slot mean-pool (model.py:254) yet are rank ~84
  vs ~5-6, so pooling is provably not the biology-specific cause. Concat breaks `self.biology` width
  (512→2048) and every checkpoint; slot reclamation bumps `query_slots` and breaks manifest/export and
  perturbs `z_patient` (model.py:288). Slot reclamation is a fine *separate* capacity PR, not a
  collapse fix.
- **Do NOT enable PCGrad (F5 mitigation half) now.** Identity and biology are disjoint heads sharing
  only the trunk; trunk-level surgery never edits the biology head, so it cannot raise biology rank.
  It is also premature until the *fixed* measurement (R5) shows sustained NEGATIVE directional cosine
  (likely a null result — the collapse is low-rank supervision, not conflict). Every-step PCGrad also
  keeps `retain_graph=True` → OOM risk on uncapped patches.
- **Do NOT drop the F6-A normalization change into a late checkpoint.** The biology head's input scale
  changes; it is not weight-equivalent. Retrain from scratch/warmup.
- **Do NOT keep any nonzero `programme_neighbourhood_loss` weight if you go to R3.** Any nonzero weight
  re-imposes the rank-5-6 target Gram (Lane 2:111-114).

---

## Summary (top-3 ranked fixes + first change)

The observed failure is biology effective rank ~5-6 of 256 (identity ~84) and a confounded metric.
Three theories compete — normalize inconsistency (T1), neighbour-KL low-rank distillation (T2), and a
missing decorrelation term (T3) — and every candidate fix is disputed as a possible silent no-op. The
adjudicating instrument is identical in all cases: an effective-rank test on `z_biology`, which does
not exist. **First change: build `morpheus/v2/tests/test_stress_collapse.py` (effective-rank +
NaN/log_variance guards), verified to fail on synthetic rank-deficient input.** Then, ranked by
payoff×minimality: **(1)** normalize the biology state so NLL and the rank-spreaders share one
geometry — model.py:283, ~2 lines, zero blast radius, gated by the new rank + log_variance<7.9 tests
(retrain from warmup, do not hot-swap). **(2)** add the VICReg off-diagonal feature-decorrelation term
on normalized `z_biology` — losses.py + training.py:232-237, weight 0.04, gated in full AND
programme_only with a min-B≥8 guard; the mechanism-matched anti-rank term. **(3)** if collapse
persists, replace neighbour-KL with an RNA-paired biology InfoNCE (training.py:157-160, neighbour-KL
weight→0), reusing the already-computed RNA view. Land the honest within-cancer metric (R4) in
parallel. Explicitly reject the separation bump, 256-D contrastive, SigLIP-as-rank-fix, slot pooling,
and PCGrad.
