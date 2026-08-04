> **SUPERSEDED IN PART — see `PHASE1B_TARGETED_READOUT.md`.**
> Every `detection_floor` below is NaN because the spike readout was a top-CCA *maximum* while the
> spike lived on one direction; that is fixed and floors now resolve. **F2 (below) is WITHDRAWN as an
> objective claim**: the identity head is numerically invariant to molecular supervision
> (max|diff| 2.6e-04 between `full` and `identity_only`), so it is the frozen MLP-CLIP teacher and
> F2 restates "the teacher beats our biology head". F1 and the confound-removal check stand.

# CALIBRA Phase 1 — first calibrated measurement

**Cohort:** 2,530 held-out-cancer TCGA patients (21 test cancers, disjoint from the **11** dev
cancers — this file said "14" until 2026-08-05; the artifact has 11, which is the documented
11-train / 21-test partition, and the correction is recorded rather than silently applied).
**Artifact:** `runs_misc/calibra_run/artifacts/diagnostic_full_seed42.npz` on persistent storage,
SHA-256 `72dcefcf05482288e4a353f7697678b9f82f7648078e223345eb3f6444b82c71`, against
`frozen_rna_targets.npz` SHA-256 `d526a7adc7456ac4f0e5e3ff71c0ef2bac96dc8488435ea714ba9840d8b51fb2`;
split totals train 3,124 / val 538 / test 2,530 = 6,192, the **pre-rebuild** cohort. Run output:
`runs/calibra_v2_local/` (commit `4c7166b`).
**Adjustment:** cancer type + tissue-source-site (sites with <10 patients pooled to OTHER → 75 sites
kept, 99 confound columns), cross-fitted ridge.
**Statistic:** top canonical correlation, 16 whitened components/side, directions fit on one half and
scored on the held-out half. **Null:** 60 permutations of patient pairing *within cancer strata*.

| state | in-sample | **held-out** | null median | excess | perm p |
|---|---:|---:|---:|---:|---:|
| **wsi_biology** | 0.520 | **0.477** | 0.151 | +0.369 | 0.016 |
| **wsi_identity** | 0.598 | **0.539** | 0.154 | +0.444 | 0.016 |
| full_biology | 0.891 | 0.876 | 0.155 | +0.736 | 0.016 |
| full_patient | 0.847 | 0.824 | 0.155 | +0.692 | 0.016 |
| rna_biology * | 0.905 | 0.898 | 0.156 | +0.750 | 0.016 |
| rna_identity * | 0.845 | 0.846 | 0.158 | +0.687 | 0.016 |

\* RNA-input states predicting RNA-derived targets are circular; reference only, excluded from claims.
perm p = 0.016 is the resolution floor at 60 permutations (no permutation exceeded the observed).

## Two findings

**F1 — the morphology↔molecular channel survives confound adjustment and generalises.**
Held-out 0.477 vs 0.151 chance for WSI-only features. The held-out estimate is close to in-sample
(0.520), so this is **not** capacity inflation. "Full" would overstate it: the adjustment is
conditional-mean only. The stronger and now-measured statement is that the channel survives the
*upper bound* on conditional-mean adjustment — a saturated cancer × site cell design moves it
0.6052 → 0.6051 — while the confound labels alone reach only 6.0–11.2% of its excess over null.
See "Validity checks passed".

**F2 — the head trained for biology is WORSE at biology than the head trained for identity.**
wsi_identity 0.539 > wsi_biology 0.477, and the gap *widens* held-out. The retrieval-trained head
carries more molecular signal than the head explicitly supervised on molecular programmes. This is the
functional counterpart of the earlier geometric finding (biology effective rank ~38 vs identity ~191).

## Validity checks passed

> **`0.463 → 0.035` IS WITHDRAWN — it could not be reproduced, and the measured values replace it.**
> This file was the origin of the pair P1 §4.2 and P1_FIGURES panel (d) quoted, stated in prose with
> no artifact path and no hash.
>
> The **cohort** has since been identified, by structure rather than by matching the number:
> `runs_misc/calibra_run/artifacts/diagnostic_full_seed42.npz`, SHA-256
> `72dcefcf05482288e4a353f7697678b9f82f7648078e223345eb3f6444b82c71`, under `run_calibra`'s own
> selection rule (`split == "test"` intersected with the frozen RNA target table) gives exactly
> **2,530** patients, **21** cancer types, **99** confound columns and **75** TSS sites kept — six of
> six against this file's header. Its split totals 6,192, i.e. the **pre-rebuild** cohort, and that
> alone explains the n = 2,530 vs n = 2,766 gap against the site arm, which runs on the rebuilt
> 6,427-patient split.
>
> The **numbers** were not reproduced. Re-run on that cohort with the project's own canonical
> estimator (`confound_certificate.lda_oof_balanced_accuracy`, standardised, 5 stratified folds,
> seed 42) the `wsi_biology` state gives **0.734 → 0.031**, not 0.463 → 0.035. None of the 56
> readings on this cohort (7 states x 8 estimator variants) has both published endpoints; the published "before" lies in the *nonlinear* band and
> its "after" in the *linear* band (table below). Only the chance rate reproduces, exactly:
> 1/21 = 0.047619. The published pair is withdrawn rather than rounded toward, and it is not
> attributed to whichever reading happens to sit nearest it.
>
> *Why it cannot be traced:* the pair was committed on 2026-07-30 (`4c7166b`, in the commit message).
> `v2/calibra/confound_certificate.py`, which defines the only balanced-accuracy functions in the
> repository, was written on 2026-08-02. `run_calibra` has never emitted such a metric, the committed
> run for this file (`runs/calibra_v2_local/`) contains no balanced accuracy in any of its 122 data rows,
> and no script computing one exists in git history, live or deleted. It came from a session probe
> that never persisted, against code that no longer exists.
>
> Detail: `NOTEBOOK_ENTRIES/p1_cancer_type_pair_withdrawn_20260805T0230Z.md`. Artifact:
> `v2/research/rebase/nature/p1_cancer_type/out/P1_CANCER_TYPE_CERTIFICATE.json`, produced by
> `v2/research/rebase/nature/p1_cancer_type_certificate.py`. Pinned in
> `v2/tests/test_paper_artifact_digests.py`.

- **Confound removal verified, not assumed — and bounded, not absolute:** on the `wsi_biology` state,
  cancer-type balanced accuracy from the residualised representation drops to **0.031** (chance
  0.048) from **0.734** raw, a **23.8×** drop to below chance. That is a statement about the **first
  moment**: the adjustment removes cancer from the class means, which is what a mean-based scorer
  such as LDA can see. It is *not* the same as "cancer is gone", and the stronger sentence is
  **refuted** — and the like-for-like nonlinear reading on this same cohort, state and folds is much
  weaker than the mean-based one: k-NN falls only **0.445 → 0.177**, i.e. 2.5×, leaving the adjusted
  state at **3.7× chance**; the random forest 0.524 → 0.268 (5.6× chance) and the RBF-SVM
  0.616 → 0.296 (6.2× chance).

  | estimator (`wsi_biology`, n = 2,530, 21 classes, chance 0.0476) | raw | adjusted | drop |
  |---|---:|---:|---:|
  | **joint LDA, standardised** (the certificate) | **0.7339** | **0.0308** | **23.8×** |
  | joint LDA, unstandardised | 0.7283 | 0.0333 | 21.9× |
  | k-NN, k = 15 | 0.4447 | 0.1766 | 2.5× |
  | k-NN, k = 15, prior-corrected | 0.5083 | 0.2284 | 2.2× |
  | random forest, 300 trees | 0.5240 | 0.2681 | 2.0× |
  | RBF-SVM | 0.6157 | 0.2958 | 2.1× |
  | per-axis nearest-class-mean, max over 256 axes | 0.1650 | 0.0516 | 3.2× |
  | *withdrawn, published 2026-07-30* | *0.463* | *0.035* | *13.2×* |

  The withdrawn row is the reason the replacement matters rather than being a decimal correction:
  its 13.2× ratio is not a ratio any single estimator on this cohort returns. The nonlinear ratios in
  the paragraph above are **raw ratios to chance on this cohort, not netted against a null that
  regenerates the adjustment**, so they are not the same quantity as the 3.45× below and must not be
  quoted interchangeably with it. On the separately measured n = 2,766 arm a nonlinear probe recovers cancer
  from the adjusted state at **3.45× chance** and site at **3.15× chance**, netted against a null
  that regenerates the adjustment inside every permutation, *p* at the resolution floor, with
  k-NN, random forest and RBF-SVM agreeing. An equal-means synthetic confirms the mechanism:
  where classes differ only in conditional variance, LDA reads 0.231 against a chance of 0.250
  while k-NN reads 0.554, the forest 0.625 and the SVM 0.658 — a mean-based certificate cannot
  see this by construction.
  **The residual is real, and it is bounded.** A saturated cancer × site cell design — which spans
  every function of the confound labels and therefore upper-bounds *any* conditional-mean
  adjustment — moves the channel by 0.0001, **0.6052 → 0.6051**; and the confound labels on their
  own account for only **6.0–11.2%** of the channel's excess over its own null. The surviving
  confound cannot explain the finding.
  *See `NOTEBOOK_ENTRIES/tcga_nonlinear_confound_probe_result_20260804T2100Z.md` and
  `NOTEBOOK_ENTRIES/nonlinear_adjustment_channel_result_20260804T2130Z.md`.*
- **Held-out canonical directions** — removes the in-sample maximisation bias.
- **Within-cancer permutation null** — preserves cohort structure, destroys only patient pairing.
- **Rare-site pooling** — prevents singleton-site dummies acting as per-patient indicators.

## Open / next
- perm p is at its resolution floor; more permutations for a precise p.
- The 0.477 is a *multivariate maximum*; it is not comparable to the ~+0.07 per-target univariate
  within-cancer specificity. Both are true — the multivariate view finds a channel that per-target
  analysis dilutes. Reconciling the two quantitatively is the next analysis.
- Purity / morphological-footprint (mRNAsi) not yet in the adjustment set.
- Second target modality (spatial) and external cohorts not yet run.
