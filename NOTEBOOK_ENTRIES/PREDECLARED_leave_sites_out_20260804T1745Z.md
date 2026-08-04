# PREDECLARED — leave-sites-out generalisation of the channel

**UTC** 2026-08-04T17:45Z
**Commit at predeclaration** `80076b4`
**Workspace** `~/ws_lso/morpheus` on `150.136.45.194`, from `git -c core.autocrlf=false archive HEAD`,
586/586 files verified against `git ls-tree -r HEAD` blob SHA-1, zero mismatches.

Nothing in this entry has been run. The reading below is fixed before the first statistic.

---

## 1. What this test is, and what it is not

`claim_guards.no_external_cohort` blocks `legible_axis` and `gene_attribution`. Its stated
mechanism is: *"a site-specific artifact would survive the current checks intact."* A true
external cohort is weeks away. The *purpose* of an external cohort is to ask whether the
result survives different acquisition conditions. That question can be asked today by
holding out whole tissue source sites.

**This does not discharge `no_external_cohort` and `claim_guards.py` is not edited.** Same
institution network, same protocols, same era, same extraction pipeline, same encoder.
It is evidence about acquisition-condition generalisation, which is the main — not the
only — thing an external cohort tests.

**Two different generalisation tests, named distinctly everywhere:**

| Name | Held out | Where it already exists |
|---|---|---|
| **held-out-cancer** | whole cancer types | the existing `tumor_state_heldout_cancer` split; the `test` partition (2,766 patients, 21 cancers) |
| **held-out-site** (this test) | whole tissue source sites | new, this entry |

## 2. The structural fact that determines the design

Measured on `d2_h_seed42.npz` (n=6,427), via `residualise.pooled_tissue_source_site`:

- 610 raw TSS codes overall; 186 pooled classes at `min_site_count=10`.
- `train` 3,118 patients / 246 sites / 11 cancers; `val` 543 / 159 / 11; `test` 2,766 / 352 / 21.
- **Sites in `test` that also appear in `train`+`val`: 0 of 352.**
- **Sites contributing more than one cancer: 0 of 352.**

TCGA assigns a TSS code per submitting site *per disease*, so **site is nested within
cancer**. Two consequences, both load-bearing:

1. Because the existing split holds out cancers, the encoder has **already** never seen any
   `test` site. The existing channel number is therefore already measured on
   encoder-unseen sites — but confounded with unseen cancers, so it cannot be read as a
   site result. This test separates the two.
2. A naive "hold out whole sites" over the pooled cohort would remove whole *cancers* and
   silently re-run the held-out-cancer test. Folds must therefore be built
   **within cancer**.

All 21 test cancers have >= 2 sites (minimum 3, ACC), so within-cancer site folds exist for
every patient: 2,766 of 2,766 are in a cancer with >= 2 sites.

## 3. Design

**Cohort.** The `test` partition of `d2_h_seed42.npz`, n=2,766, intersected with
`frozen_rna_targets.npz`. Encoder saw none of these patients, none of these sites, none of
these cancers. Targets: the `all_non_control` block, `RANDOM_CONTROL__` columns excluded
exactly as `run_calibra.main` excludes them.

**States.** Primary **`wsi_biology`** — the morphology->molecular channel, the only one that
bears on the claim. `full_biology` and `rna_biology` reported alongside; `rna_biology` is
the circular positive control and is expected to survive trivially. If it does not, the
harness is broken and no other number in this entry may be read.

**Folds.** K=5, **whole sites, stratified within cancer**. Deterministic and seed-free: for
each cancer, sites are sorted by descending patient count (ties by site code) and each site
is assigned to the fold currently holding fewest patients *of that cancer*. Whole sites move
together. Per-fold site and patient counts, and the per-cancer breakdown, are reported.

**Statistic — the channel on held-out sites.** `heldout_top_cca_indexed` at
`n_components=16` (the paper's channel budget): PCA-whitening maps and the top canonical
direction pair are fit on the **training-site** patients only and applied to the
**held-out-site** patients; the statistic is `|corr|` of the projected pair on held-out
sites. Directions are never fit on the rows they are graded on.

**Arms, both run.**
- **Unadjusted** — no residualisation. If the channel survives here on unseen sites that is
  much stronger than surviving after adjustment, because site information is free to be
  used and simply does not transfer to sites never seen.
- **Adjusted** — `cancer + pooled TSS` cross-fitted residualisation, CALIBRA's own design,
  applied **separately within the training-site block and within the held-out-site block**.
  A held-out site has no column in a design fit on training sites, so it cannot be adjusted
  from there; adjusting each block by its own site structure is exactly what one does with
  a real external cohort. Because site is nested in cancer, the TSS dummies span the cancer
  dummies, so this arm removes cancer as well and measures within-site, within-cancer
  covariation only. That is the most stringent reading and is stated as such.

**Null.** Within-cancer permutation of the cross-modal pairing among the **held-out-site**
patients, rescored through the same fitted directions. `n_permutations=1000`, p floored at
1/1001 and the resolution reported. Same convention as `calibration.permutation_null` and
`confound_certificate.within_stratum_permutations`.

**Intervals.**
- **Site-cluster bootstrap (primary)** — resample held-out *sites* with replacement, each
  bringing all its patients. This is the interval that answers "would another draw of sites
  give this?".
- Patient bootstrap (secondary), and the across-fold spread of the 5 fold values.
- 1,000 resamples each, percentile 2.5/97.5.

**Matched comparator, predeclared.** The same estimator on a **random patient split matched
to the site folds** — identical per-(cancer, fold) counts, patients permuted within cancer,
site ignored. This isolates the cost of *site shift* from the cost of fitting the read-out
out-of-sample. Without it a drop cannot be attributed to sites.

## 4. The reading, fixed now

Per arm, on `wsi_biology`, pooled across the 5 folds and per fold:

- **SURVIVES** — observed > permutation null p95 **and** the site-cluster bootstrap 95% CI
  lower bound > null p95, in **>= 3 of 5** folds.
- **COLLAPSES** — observed <= null p95, or the CI covers null p95, in **>= 3 of 5** folds.
- **ATTENUATED BUT PRESENT** — survives as above, but the ratio of the held-out-site channel
  to the matched random-split channel is **< 0.5**. A ratio near 1 means site shift cost
  essentially nothing; a ratio well below 1 is a real acquisition-condition penalty that
  must be quoted even though the channel cleared its null.

A survival in the **unadjusted** arm is the strong result. A survival only in the
**adjusted** arm is weaker and will be reported as such, not merged with it.

Both arms are reported whatever they show. If the channel collapses across sites, that is
the headline and it is reported first.

## 5. Threats this test does not answer

- Same institution network, same protocols, same era. `no_external_cohort` stands.
- The encoder is frozen and is not retrained per fold; this is a **read-out** transfer test.
  The encoder having never seen any test site (section 2) is what keeps that honest, but a
  retrained-per-fold version would be a stronger test and is not run here.
- TCGA-to-elsewhere domain shift is far larger than TCGA site-to-site shift. The project's
  own measurement: a cohort classifier separates TCGA from HEST at AUC 0.99999
  (`spatial_baselines_20260803T0620Z.md`). Surviving site shift inside TCGA does not
  predict surviving that.
