> **SUPERSEDED IN PART — see `PHASE1B_TARGETED_READOUT.md`.**
> Every `detection_floor` below is NaN because the spike readout was a top-CCA *maximum* while the
> spike lived on one direction; that is fixed and floors now resolve. **F2 (below) is WITHDRAWN as an
> objective claim**: the identity head is numerically invariant to molecular supervision
> (max|diff| 2.6e-04 between `full` and `identity_only`), so it is the frozen MLP-CLIP teacher and
> F2 restates "the teacher beats our biology head". F1 and the confound-removal check stand.

# CALIBRA Phase 1 — first calibrated measurement

**Cohort:** 2,530 held-out-cancer TCGA patients (21 test cancers, disjoint from the 14 dev cancers).
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

**F1 — the morphology↔molecular channel survives full confound adjustment and generalises.**
Held-out 0.477 vs 0.151 chance for WSI-only features. The held-out estimate is close to in-sample
(0.520), so this is **not** capacity inflation.

**F2 — the head trained for biology is WORSE at biology than the head trained for identity.**
wsi_identity 0.539 > wsi_biology 0.477, and the gap *widens* held-out. The retrieval-trained head
carries more molecular signal than the head explicitly supervised on molecular programmes. This is the
functional counterpart of the earlier geometric finding (biology effective rank ~38 vs identity ~191).

## Validity checks passed
- **Confound removal verified, not assumed:** cancer-type balanced accuracy from the residualised
  representation drops to **0.035** (chance 0.048) from **0.463** raw. Cancer is gone.
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
