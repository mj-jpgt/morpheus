# PREDECLARED — replacing P4's counterfactual "28 of 90" with a measured number from a real inductive adjustment operator

**Written:** 2026-08-04 22:45 UTC, at commit `9025045`, **before any Test B or Test C number below
has been computed** and before the wiring in `p4_certify.py` was written. Only *setup composition*
has been read off the artifact so far (partition sizes, cancer levels, TSS counts, fold sizes) — no
certificate, no channel correlation, no floor.

Supersedes nothing; it is the follow-up to
`NOTEBOOK_ENTRIES/p4_certification_end_to_end_20260804T2000Z.md`, which reported **0 of 90** queries
answerable and a **counterfactual 28 of 90** obtained by *relaxing* condition 3 rather than by
running an operator.

---

## 1. The question

The 20:00 entry named one buildable gap: `cross_fitted_residuals` is transductive, so the adjusted
state that passes certification cannot be produced for a patient outside the cohort it was fitted on.
It then computed what a certified interface *would* answer if that operator existed — 28 of 90 — by
dropping condition 3 from the policy. That number is a projection. This work replaces it with a
measurement.

## 2. Two things established before predeclaring, because they dictate the design

**(a) The operator already exists and P4 does not need a second one.**
`v2/calibra/inductive_adjustment.ConfoundAdjustmentOperator` is exactly the required map: a frozen
`DesignSpec` (one-hot cancer + pooled TSS, levels and any numeric moments frozen against the
reference), a frozen `SitePooling` (unseen site → `OTHER`, which is the existing pooling rule
evaluated at a count of zero), and the K persisted ridge models plus the reference centre. Its
`adjust_reference` is bit-identical to `cross_fitted_residuals` and its `adjust` applies the same
fitted map to rows the operator never saw. Fourteen instances of it are already fitted and persisted
against TCGA (`runs_misc/tcga_operators/`, entry `tcga_operator_fitted_and_persisted_20260805T0045Z`).
**No new operator will be written.** The work is a split, a wiring, and a measurement.

**(b) P4's condition 1 split cannot serve as the operator's discovery fold. Measured, not assumed.**
On `runs/d2_final/artifacts/d2_h_seed42.npz`:

| | train+val | test |
|---|---:|---:|
| patients | 3,661 | 2,766 |
| cancer types | 11 (BRCA, GBM, HNSC, KIRC, LGG, LUAD, LUSC, OV, PRAD, THCA, UCEC) | 21 (ACC … UVM) |
| raw TSS codes | 258 | 352 |
| **cancer types shared with the other side** | **0** | **0** |
| **TSS codes shared with the other side** | **0** | **0** |

D2's split holds out whole **cancers**, and TCGA assigns a TSS code per site *per disease*, so the
site sets are disjoint too. An operator fitted on the representation's own discovery fold would see
every test cancer as an unseen level and every test site as `OTHER`: its adjustment would degenerate
to subtracting a constant. **So conditions 1 and 3 of the five-point rule cannot be discharged by the
same split**, and the honest inductive question for P4 has to be asked *inside* the certification
cohort: fit the nuisance operator on a discovery half of the test partition, apply it unchanged to
the other half, and certify and query only that other half.

## 3. Protocol, fixed here

**Cohort.** `runs/d2_final/artifacts/d2_h_seed42.npz`, state `wsi_biology`, `split == "test"`,
n = 2,766. Targets `data/frozen_rna_targets.npz`, the same **90** non-control targets
(50 `hallmark_in_training`, 24 `heldout_pathway`, 8 `immune_tme`, 8 `tumour_state`).

**Exposure split.** Patients permuted *within cancer type* at `seed = 42` and dealt into a discovery
fold **D** and an exposure fold **E**. Primary fraction **0.5** (n_D = 1,384, n_E = 1,382; all 21
cancers present on both sides). Sensitivity arm at **0.7** (n_D = 1,937, n_E = 829) if the primary
completes with CPU to spare.

**Operator.** `ConfoundAdjustmentOperator.fit(features[D], frame={"cancer": cancers[D]},
columns=["cancer", "tss"], patient_ids=ids[D], site_column="tss", min_site_count=10, n_splits=5,
alpha=1.0, seed=42)` — i.e. the *published* 5-fold ridge on the one-hot cancer + pooled-TSS design,
unchanged. A second operator, fitted identically on the D rows of the 90-column target block, adjusts
the targets. Both are applied to E with `.adjust(...)`, and the adjusted E block is what everything
downstream sees.

**Three arms, each scored end to end.**

| arm | adjusted state | what it isolates |
|---|---|---|
| `transductive_full` | `cross_fitted_residuals` on all 2,766 test rows | the published 20:00 run; the **28** baseline. Not re-run — read from the vendored `out/C_competitor_d2h_seed42_wsi.json` |
| `transductive_exposure` | `cross_fitted_residuals` on the **1,382 E rows only** | the matched-n control. Same rows, same n, same design width as the inductive arm, but fitted *on the rows it scores* |
| `inductive_exposure` | operator fitted on D, applied to E | the real thing |

The matched control is not optional: without it a drop in the inductive arm confounds *inductive vs
transductive* with *n halved*, and n halving alone raises every detection floor.

**Detection floors are computed identically in both exposure arms** — `spike_recovery_curve` on the
E rows against the E design, 25 draws, level grid (0, .01, .02, .05, .10, .20, .40). They are
therefore *the same numbers* in both arms by construction and cannot explain any difference between
them. That the floor machinery is itself transductive is a stated limitation, not a hidden one: the
floor is a property of the audit, not of the state a patient is shown.

**Condition 3 in the inductive arm** is `certify_axes(adjusted_E, ids[E], cancers[E],
residualise=False, min_site_count=10, n_permutations=1000, n_boot=200)` — the certificate run on the
state that would actually be exposed, with no further in-sample residualisation. Site classes for the
certificate are defined on **E**, by E's own `min_site_count = 10` pooling; the *adjustment* uses D's
frozen pooling. That asymmetry is deliberate and is the deployed situation: the operator's site list
is frozen at fit time, the auditor defines site classes on the cohort in front of them.

**Policy C, unchanged from the 20:00 entry:** a query is answered iff its supporting axis passes the
confound certificate on the exposed state **and** `|heldout_single_direction_correlation|` exceeds
that (axis, target) pair's own detection floor **and** it exceeds the 95th percentile of a 200-draw
within-cancer row-permutation null. Both the full count and the count with condition 3 relaxed will
be reported, the latter being the direct successor to the 28.

## 4. The number I expect the operator to be limited by, stated in advance

The TCGA test partition has **352 raw TSS codes over 2,766 patients; median site size 3**. Only 84
sites reach `min_site_count = 10` on the whole partition, and on a half-sized discovery fold only
**31** do. Measured on the frozen split: **606 of the 1,382 exposure rows (43.8%)** come from a site
that is frequent in D and therefore receive a genuine site adjustment; the other **776 (56.2%) are
pooled to `OTHER` and receive a cancer adjustment only.** At 0.7 the coverage rises only to 54.8% on
829 rows.

This is not an artefact of my split. It is TCGA's site-size distribution, and it is the structural
reason an inductive site adjustment is intrinsically weaker than a transductive one: the transductive
path gets to define its site levels on the very rows it is adjusting, so its coverage is 100% by
construction.

## 5. Predictions

**Primary prediction: the inductively adjusted state FAILS the joint row of the confound
certificate, and the measured answerable count is 0 of 90.** Probability I would put on it: **0.65**.
Mechanism: over half the exposure rows get no site adjustment at all, and the ones that do get a site
mean estimated from ~10–30 discovery patients, so a joint LDA over 256 axes should still find site
structure well above a within-cancer permutation null.

**Point predictions, to be graded verbatim:**

| quantity | prediction |
|---|---|
| joint LDA balanced accuracy, `inductive_exposure` | **0.06–0.20** (raw arm reads 0.3633 at n = 2,766; transductive-adjusted reads 0.0118) |
| joint LDA balanced accuracy, `transductive_exposure` | **≤ 0.03**, certified |
| axes breaching per-axis, `inductive_exposure` | **5–40** of 256 |
| **`n_C` full policy, `inductive_exposure`** | **0** |
| **`n_C` with condition 3 relaxed, `inductive_exposure`** | **8–20** of 90 (vs 28 at n = 2,766) |
| **`n_C` with condition 3 relaxed, `transductive_exposure`** | **10–22** of 90 |
| category breakdown of whatever clears the floor, either exposure arm | ≥ 60% `hallmark_in_training`; **0 or 1** `heldout_pathway` |
| targets with a NaN (unresolvable) detection floor on E | **25–40** of 90 (23 of 90 at n = 2,766; halving n can only make it worse) |

**Conditional branch.** If the inductive arm *does* certify, the reportable number is its full-policy
`n_C`, and I predict it lands **below** the transductive-exposure relaxed count — i.e. the real
operator costs queries relative to the counterfactual, it does not gain them.

## 6. What would make me distrust a favourable result

A favourable result here is "the inductively adjusted state certifies and answers a decent number of
queries". I will not report that as a pass if any of the following holds, and each is checked:

1. **The adjustment is doing nothing.** If the inductive arm's adjusted E block is nearly the raw E
   block (correlation with the raw block > 0.99 per axis, or if the fitted design has collapsed to
   ~22 cancer columns because almost no site cleared the threshold), then a "PASS" on the certificate
   would be a statement about a null adjustment. The design width and the per-axis raw-vs-adjusted
   correlation will be reported for both exposure arms.
2. **The adjustment is doing too much.** Test A′ of the 20:00 entry already showed that the *adjusted*
   certificate is near-inert against a planted site code at any realistic strength — a blanket PASS on
   an adjusted state is a statement about the adjustment and the classifier family, not evidence that
   an axis is site-free. **The planted-axis ladder will be re-run on the inductive arm** (site code,
   cancer code, noise, SNR 0.5/1.0/2.0/∞). If a noiseless planted site code does **not** breach on the
   inductively adjusted state, the certificate is inert there and no PASS from it may be quoted.
3. **Leakage through the split.** D and E must share no patient. If any patient id appears in both, or
   if the operator's `reference_digest` does not match the D rows, the run is void. Asserted in code.
4. **A count that beats the counterfactual.** If the measured relaxed count exceeds 28 despite half
   the n, that is a red flag, not a win: it would most likely mean the inductive residual has *more*
   variance left in it (less signal removed), inflating correlations. In that case the per-axis
   raw-vs-adjusted correlations and the residual variance ratio decide, and the finding is reported
   as an artefact unless they exonerate it.
5. **The `heldout_pathway` count moving up.** 1 of 24 genuinely untrained pathways survived at
   n = 2,766. If more than 3 survive at n = 1,382 without an accompanying rise in the trained
   categories, I will treat it as a fold-composition accident and say so.

## 7. Recorded failure conditions

* If E carries a cancer level absent from D, `DesignSpec.transform` raises `UnseenLevelError` and the
  run stops. It is not silenced with `on_unseen_level="zero"`; the stratified-within-cancer split
  exists precisely so this cannot arise, and if it does the split is wrong.
* If the GPU box is loaded such that the CPU arms cannot complete, the completed arms are reported
  and the missing ones are named as not run — never inferred.
* `claim_guards.py` and `claim_evidence.json` are not touched by this work.
