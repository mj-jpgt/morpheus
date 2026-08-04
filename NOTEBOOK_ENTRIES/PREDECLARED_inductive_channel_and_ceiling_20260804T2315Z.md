# PREDECLARED — does P1's morphology→molecular channel survive an INDUCTIVE confound adjustment, or only a transductive one?

**Written:** 2026-08-04 23:15 UTC, at commit `742fc3a`, **before the wiring was written and before any
number in the result entry existed.** Only *setup composition* has been read: that the local tree has
no artifacts, that `runs/d2_final/artifacts/d2_{h,i}_seed42.npz` and `data/frozen_rna_targets.npz`
exist on the A100 box (`150.136.45.194`, 30 cores), and that `v2/calibra/nonlinear_adjustment.py`'s
`channel_under_adjustment` calls its `adjust` on the same rows it scores. No channel number, no null,
no ceiling has been computed.

---

## 1. The question, and why it is now open

`NOTEBOOK_ENTRIES/p4_inductive_adjustment_measured_20260804T2300Z.md` measured that **every CALIBRA
number computed on an "adjusted" state on this project was computed transductively**: the
confound-removal nuisance model (5-fold ridge on a one-hot cancer + pooled-TSS design) is fitted on
the very rows it scores. On TCGA site recoverability the difference is a factor of **24** — joint LDA
balanced accuracy 0.0109 in sample against **0.2643** out of sample, on the identical 1,382 patients.

`NOTEBOOK_ENTRIES/nonlinear_adjustment_channel_result_20260804T2130Z.md` is P1's flagship result:
`d2_h::wsi_biology` reads top-CCA **0.6052 → 0.6051** under a saturated cancer×site cell-mean design
that provably upper-bounds any conditional-mean adjustment, retention of excess 0.987–1.007 across
eleven adjusted arms. Its second leg is a **labels-only ceiling**: a representation that is *nothing
but* the confound labels, pushed through the identical pipeline, scores **0.1237** (additive
108-column design) or **0.0903** (saturated 105-column cell design) after adjustment — *below* the
real channel's own within-cancer pairing-null median of **0.1483**, i.e. **11.2%** and **6.0%** of the
channel's excess over that null.

**Both legs were computed transductively.** Verified against the code path, not assumed:
`nonlinear_adjustment.channel_under_adjustment` computes `x_adjusted = adjust(x)`,
`y_adjusted = adjust(y)` and `adjust(y[order])` per permutation, and every `adjust` built by
`make_adjuster` — `ridge` (which *is* `residualise.cross_fitted_residuals`), `saturated`,
`kernel_ridge`, `forest`, `location_scale` — fits its nuisance model inside the call on the rows
handed to it. `labels_only_ceiling` routes through the same function. There is no inductive path in
that module.

**So the open question is:** does the channel still stand above what confound-label information alone
can explain, when the nuisance model is fitted on a *separate* discovery fold and applied out of
sample? If site information comes back at 24× under an honest inductive adjustment, the confound-only
representation's post-adjustment score is free to rise too — and *that* number, not the channel's own
retention, is what decides whether P1's headline sentence needs the same narrowing the
site-recoverability claim just got.

## 2. Protocol, fixed here

**Cohort and split — reused, not rebuilt.** `p4_certify.exposure_split` and
`p4_certify.prepare_state` are called unchanged: patients of the `test` partition permuted **within
cancer type** at `seed = 42`, discovery fraction **0.5** (n_D = 1,384, n_E = 1,382), which is the
identical split the P4 inductive entry measured 0.2643 on. D2's own train+val/test split cannot serve,
because it shares **0 cancer types and 0 TSS codes** with `test` (P4 predeclaration §2b), so the
nuisance operator needs its own discovery fold carved out of the certification cohort.

**Operators — reused, not rebuilt.** `inductive_adjustment.ConfoundAdjustmentOperator.fit` on the D
rows (5-fold ridge, `alpha = 1.0`, one-hot cancer + pooled TSS at `min_site_count = 10`, frozen
levels, frozen site pooling, `on_unseen_level="refuse"`), one operator for the 256-column image block
and a second, fitted identically, for the 90-column target block. **No second operator is written.**

**Provenance assertions, checked in code and reported, run void if they fail.**

* 0 patients in both D and E.
* The image operator's `reference_digest` equals the P4 entry's recorded
  `2060a635fa83756a1c3b7aa8506b7b19fcc4431f5d1a303da39b3cb2bf9d62ce`. If it does not, this run is not
  on P4's state and the comparison to 0.2643 is withdrawn.
* The adjuster closures this run builds reproduce `p4_certify.prepare_state(...)`'s
  `adjusted_features` and `adjusted_targets` **bit-for-bit** (`np.array_equal`).

**Statistic — the canonical one, no inline arithmetic.**
`nonlinear_adjustment.channel_under_adjustment`, which is pinned by
`test_channel_under_incumbent_reproduces_permutation_null_exactly` to reproduce
`calibration.permutation_null` under the incumbent adjuster. It reads
`spectral.top_canonical_correlation(·, ·, n_components=16)` (S1, the §4.4 headline),
`spectral.heldout_top_cca` (S2), `spectral.effective_rank`, and the within-cancer **pairing**
permutation null at **2,000** permutations (p floor 1/2001 = 0.0005). The global pairing null is
computed beside it. Retention is `nonlinear_adjustment.retention_of_excess` — a ratio of excess over
each arm's own null, never of raw S1. The ceiling is `nonlinear_adjustment.labels_only_ceiling`.

**One code change to a canonical module, declared here.** `channel_under_adjustment` and
`labels_only_ceiling` take a new optional `adjust_y=None` which defaults to `adjust`, because an
inductive operator for a 256-column block and one for a 90-column block are necessarily different
fitted objects while the transductive `adjust` is column-count agnostic. Default behaviour is
unchanged and the existing identity test is left in place; a new test asserts
`adjust_y=None` is byte-identical to passing `adjust_y=adjust`.

**Arms.** All on `d2_h_seed42::wsi_biology`, `test` partition, targets `frozen_rna_targets.npz` with
the 90 `RANDOM_CONTROL__` columns excluded exactly as `run_calibra.py` does. `d2_i_seed42` if the box
has CPU to spare.

| arm | n | x adjusted by | y adjusted by |
|---|---:|---|---|
| `transductive_full` (**reproduction gate**) | 2,766 | `cross_fitted_residuals` on all test rows | same |
| `none_exposure` | 1,382 | column centring only | same |
| `transductive_exposure` (**matched-n control**) | 1,382 | `cross_fitted_residuals` on the E rows, E's own design | same |
| `inductive_exposure` | 1,382 | D-fitted operator applied to E | D-fitted target operator applied to E |

**Ceiling arms.** The labels-only "representation" is the **additive 108-column
`confound_design`** and the **saturated 105-column `cell_design`**, both built on the whole test
partition so D and E rows are encoded identically and so the numbers are comparable to the published
108/105-column figures. Each is scored under `transductive_exposure` and under `inductive_exposure`
(a third `ConfoundAdjustmentOperator` fitted on the D rows of the labels block), against the same
inductively adjusted targets, with the same 2,000-permutation pairing null.

**Diagnostics reported for every arm, whatever the direction:** the `p4_certify._adjustment_audit`
per-axis raw-vs-adjusted correlation and residual variance ratio; `cross_fitted_r2` of both blocks on
the labels; `effective_rank` of both adjusted blocks; `adjuster_agreement` against the matched
transductive control; S1 and S2; both nulls.

CPU only, thread caps `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1`, process parallelism at `--n-jobs 6`
on a 30-core box carrying other agents' load. If a genuine GPU need appears the run stops and says so.

## 3. What would make me say the P1 claim SURVIVES

Both of these, on the inductively adjusted state:

* **S1.** The channel's `permutation_p` is at the 0.0005 floor **and**
  `retention_of_excess(inductive_exposure, transductive_exposure) ≥ 0.8` — the channel is not
  materially attenuated by making the adjustment honest.
* **S2 — the deciding leg.** The labels-only ceiling under the *identical inductive treatment* stays
  **below the channel's own within-cancer pairing-null median**, and its excess over its own null is
  **≤ 20% of the channel's excess**. Transductively these read 0.0903/0.1237 against a null median of
  0.1483, i.e. 6.0% and 11.2%; I allow a near-doubling of that headroom before calling it a change.

## 4. What would make me say it needs the SAME NARROWING as the site-recoverability claim

Either of:

* **N1 (the one I expect to decide it).** The labels-only ceiling under inductive adjustment rises
  **above the channel's own null median**, or its excess exceeds **20%** of the channel's excess. Then
  a representation that is *nothing but* the confound labels becomes a material share of the measured
  channel once the nuisance model is not allowed to see the rows it scores, and the supportable
  sentence narrows to: *the channel exceeds the labels-only ceiling when the confound is removed from
  the first moment of the rows the nuisance model was fitted on.*
* **N2.** `retention_of_excess(inductive, transductive_exposure) < 0.8`, i.e. the channel itself is
  materially attenuated out of sample.

**Prior.** I put **0.7** on N1 firing. Mechanism, stated before measuring: only 43.8% of E rows come
from a site frequent in D (P4 predeclaration §4), so for 56.2% of rows the operator has no site column
at all and the labels block's site indicators for those rows pass through essentially untouched. A
labels-only representation is the *most* exposed thing there is to that failure. **0.35** on N2 —
the inductive residual is noisier, but noise attenuates both the channel and its null.

## 5. What would make me distrust a FAVOURABLE result

A favourable result here is "retention ≥ 0.8 and the ceiling still tiny". I will not report it as one
if any of the following holds, and every one is checked and reported:

1. **The inductive pairing null cannot regenerate what the transductive one does, and I predict this
   inflates the inductive arm's excess.** `calibration.permutation_null` re-adjusts `y` on every
   permutation, so correlation *induced* by shared residualisation is regenerated inside the null
   (P1 §4.6). Under a transductive adjuster the refit happens on the permuted rows and the induced
   floor is reproduced. Under an **inductive** adjuster `adjust(y[order])` applies a **fixed** map
   keyed to the *exposure rows' own* design, so a permuted patient's un-removed site effect is no
   longer aligned with the design at that position and the shared un-removed confound is **not**
   regenerated in the null. The inductive null should therefore sit *lower* and retention should read
   *higher* for a reason that has nothing to do with the channel. **Check:** if the inductive arm's
   null median falls more than 30% below the matched transductive control's while its raw S1 rises,
   retention is declared uninterpretable in that direction and the verdict is taken from §3's S2 leg
   (the ceiling) alone, which is immune to this because ceiling and channel are read against nulls
   built the same way.
2. **The adjustment is doing nothing.** Any axis with raw-vs-adjusted correlation > 0.99, or a median
   residual variance ratio near 1, and a "the channel survives" reading is a statement about a null
   adjustment. `_adjustment_audit` is reported for every arm. (The P4 entry measured 0.7536 median
   correlation and 0.5896 variance ratio for this operator, so I expect this not to fire.)
3. **The reproduction gate misses.** `transductive_full` must read S1 **0.6052**, null median
   **0.1483**, excess **0.4569** on `d2_h::wsi_biology` to four decimals. If it does not, nothing
   downstream is comparable to §4.4 and the run is reported as a failed gate, not as a result.
4. **Provenance.** Any D∩E overlap, any digest mismatch against P4's operator, or any failure of the
   bit-for-bit agreement with `p4_certify.prepare_state`.
5. **The two artifacts disagree.** If `d2_i` runs and points the other way, both are reported and no
   direction is claimed.
6. **Capacity.** If the inductively adjusted block's `effective_rank` collapses toward the
   16-component budget, S1 is propped up by rank collapse and is not quotable. Reported for every arm.

## 6. What would make me distrust an UNFAVOURABLE result

1. **Over-removal masquerading as a confound finding.** If the channel collapses (retention < 0.5), the
   inductive adjustment's attenuation of a *known planted* signal is measured through the identical
   operator with `calibration.spike_recovery_curve` before the collapse may be called a confound
   result. Declared now so its absence later is a consequence, not an omission.
2. **A ceiling that rises only because its own null rose.** The ceiling is graded on excess over its
   own null as well as on raw S1 against the channel's null median, and both are reported.
3. **n.** Every exposure arm sits at n = 1,382 against the published 2,766. The matched transductive
   control exists precisely so that "inductive vs transductive" is never confounded with "n halved",
   and no inductive number is compared to a published n = 2,766 number without it.

## 7. Recorded scope limits

* One cohort, one partition, `wsi_biology` only, one discovery fraction (0.5) unless CPU allows 0.7.
* `min_site_count` stays at the project default of 10; no sensitivity sweep.
* The **detection-floor machinery and the confound certificate are not re-run here** — this run
  measures the channel and the ceiling. P4's entry already carries the certificate on this exact
  state.
* `claim_guards.py`, `claim_evidence.json` and `paper/P1_CALIBRA_DRAFT.md` are **not touched**; any
  prose correction is written into the result entry for the owners of those files to apply.
