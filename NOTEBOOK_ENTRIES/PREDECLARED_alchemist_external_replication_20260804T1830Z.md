# PREDECLARED — what an ALCHEMIST result would have to look like to count

**UTC** 2026-08-04T18:30Z
**Written before any ALCHEMIST channel number exists.** At the time of writing, patch
extraction is running (42 of 1,106 slides staged) and no expression file has been scored.
The cohort-classifier AUC and the channel are both unmeasured.

**Predecessor** `NOTEBOOK_ENTRIES/external_cohort_options_verified_20260804T1900Z.md`, which
selected ALCHEMIST and recorded the tumour-polygon problem as an irreducible deviation.

---

## 0. Bad news that is already certain, stated before the bar

Three facts are known now, are not results, and each of them weakens what this experiment
can conclude. They are recorded here so that a later good number cannot be read as if they
had not applied.

1. **The cohort costs 1.789 TB, not "close to zero effort".** `api.gdc.cancer.gov/files`,
   program ALCHEMIST, `data_type = Slide Image`: 1,349 files, **2.076 TB** total
   (mean 1.54 GB, median 1.35 GB, max 4.27 GB). Restricted to one slide per paired patient:
   1,106 files, **1.788832467192 TB**. The scout entry gave no size and the plan called the
   effort delta "close to zero"; on the download alone it is ~3 h at the 170 MB/s this link
   sustains across 6 parallel GDC streams (28.3 MB/s each, measured).

2. **ALCHEMIST publishes no tissue source site at all.** `api.gdc.cancer.gov/cases` faceted
   on `tissue_source_site.code` and `tissue_source_site.name` returns `_missing` for **all
   1,176 cases**. `pooled_tissue_source_site` takes field 2 of the barcode, which for
   `ALCH-B0CW` is a per-case code, so with `min_site_count=10` every ALCHEMIST patient
   collapses to a single `OTHER` level — a constant column, not a site adjustment. **The
   site half of "cancer + pooled site" cannot be applied to ALCHEMIST.** The only batch-like
   variable that exists is the Aperio `ScanScope ID` in the slide header, which is being
   recorded per slide. It is a scanner identity, not a tissue source site, and it will not
   be renamed into one.

3. **No WSI tumour polygons exist for ALCHEMIST**, or for any non-TCGA GDC project. TCGA
   patches are restricted to pathologist-drawn tumour polygons; ALCHEMIST sampling is
   whole-slide tissue mask. This is the declared deviation from
   `external_cohort_options_verified`, and §3 below is the pre-stated penalty for it.
   The fabrication vector recorded in that entry stands: IDC's `cptac_*_tumor_annotations`
   are `Modality: RTSTRUCT` radiology contours on CT/MR/PT and are not pathology
   annotations. Nothing in this experiment will cite them.

**Two corrections to documents this task was briefed from**, reported rather than worked
around:

- The brief and `v2/research/dilution/extract_normal_patches.py:9` both point at
  **`DILUTION_CURVE.md`. That file does not exist** in the repo or anywhere in its git
  history. The seven retained-channel numbers live in
  `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2. The docstring reference is stale.
- The scout entry says ALCHEMIST histology "must come from the 1,176 clinical supplement
  files, which were not audited". It does not: **`diagnoses.primary_diagnosis` is populated
  on `/cases` for all 1,176**, giving 16 histologies (adenocarcinoma NOS 853, squamous cell
  carcinoma NOS 229, adenosquamous 24, ...). Stage genuinely is `_missing`. The cancer
  stratum for the confound design therefore comes from the API, not from an unaudited file.

---

## 1. The instrument, pinned before it is pointed at anything

Identical on both cohorts. Every number below is produced by this and nothing else.

**Patches.** 256x256 px from a fixed 128 um x 128 um field of view (crop `round(128/mpp)`
native pixels, resample to 256), JPEG quality 75, 4:2:0, no colour normalisation. Verified
empirically against the actual TCGA store rather than taken on trust: 40 randomly read
TCGA-UT patches are all 256x256, all carry JPEG component layout `(1,2,2,0),(2,1,1,1),(3,1,1,1)`
= 4:2:0, and all carry the quantisation table `[8,6,5,8,12,20,...]`, which is the IJG table at
quality 75. ALCHEMIST slides read `MPP = 0.2476-0.2527`, so `crop_px` = 507-517 native pixels
for the same 128 um. The renderer is **imported** from
`v2/research/dilution/extract_normal_patches.py`, not restated, so ALCHEMIST sits inside the
guard `v2/tests/test_hest.py::test_adapter_constants_match_the_tcga_extractor` already
applies to those constants.

**H-Optimus-0 `pretrained_cfg` centre-crops to 87.5%, so the analysed field is 112 um, not
128.** Identical on both cohorts, so comparability holds. Stated, not corrected for.

**Patch budget.** 30 patches per slide, one slide per patient. This matches the TCGA store,
whose per-slide count is 30 at the median *and* at the 10th and 90th percentiles
(271,710 patches / 8,736 slides / 7,175 patients). TCGA patients are subsampled to exactly
30 tokens (seed 42) so the two sides have the same aggregation budget; the uncapped TCGA
number is reported alongside as a sensitivity.

**Representation.** `concat(mean, std)` over a patient's 30 raw 1536-d H-Optimus-0 tokens
= 3072-d, then PCA-256 fitted within cohort. **Zero fitted parameters that see any label.**
This is deliberately the same representation the dilution curve was measured on
(`build_dilution_artifact.py`), because the dilution curve is the only calibration available
for the deviation in §3 and a penalty measured on one representation does not transfer to
another. The trained `wsi_biology` head is *not* used: `d2_compare._load` demands it, but it
only comes out of a TCGA-trained V2 checkpoint, and running that on ALCHEMIST would add a
transfer question on top of the cohort question and confound the two.

**Targets.** The 90 non-control signatures of `frozen_rna_targets.npz`
(50 `hallmark_in_training`, 24 `heldout_pathway`, 8 `immune_tme`, 8 `tumour_state`).
Scoring is `within_sample_gene_rank` with `fit_population: none` and
`cohort_fit_free_target_scoring: true` — recorded in that artifact's own `metadata_json` —
so it transfers to a new cohort without refitting anything on TCGA. The 90 random controls
are **not** used: their gene draws are not recoverable from any artifact on disk, so an
ALCHEMIST "control" would not be the same control. The permutation null in §2 is the
negative control instead.

**Adjustment.** `confound_design` on `["cancer"]`, cross-fitted residuals
(Ridge alpha=1.0, KFold n_splits=5, shuffle, seed 42), applied identically to the targets and
to the representation. `cancer` is LUAD/LUSC on the TCGA side and `primary_diagnosis` on the
ALCHEMIST side, with levels of fewer than 10 patients pooled to `OTHER` — the same
`min_site_count=10` rule `pooled_tissue_source_site` uses. **Site is omitted on both sides**,
because per §0.2 it does not exist for ALCHEMIST and an adjustment applied to only one arm is
not the same instrument. TCGA is additionally reported with `["cancer", "tss"]` so the cost
of dropping the site term is visible rather than assumed.

**Metric.** `top_canonical_correlation(residual_x, residual_y, n_components=16)`. In-sample,
therefore capacity-inflated, therefore always reported null-corrected.

**Null.** `permutation_null`, rows of the target matrix permuted **within cancer strata**,
300 permutations, seed 42.

**Comparator, and this is the load-bearing design choice.** The headline TCGA number
(pan-cancer, 6,427 patients, 32 cancers, adjusted top-CCA 0.5573, null median ~0.146) is
**not** the right bar for ALCHEMIST, and using it would guarantee a spurious failure.
ALCHEMIST is one disease. A pan-cancer channel is measured across 32 tumour types whose
morphology and expression differ enormously; residualising cancer type does not remove
disease-breadth from the variance that CCA has to work with. So the matched comparator is
**TCGA LUAD + LUSC, tumour-polygon-restricted: 846 patients** (420 LUAD, 426 LUSC), all 846
present in both the frozen target table and the patch store. n = 846 against ALCHEMIST's
1,106 is a close match, so the capacity inflation of in-sample CCA is comparable on both
sides as well. The pan-cancer number is reported for context and is explicitly not the bar.

---

## 2. Validation gates that must pass before any channel number is quoted

These are not results. If any fails, the experiment is **uninterpretable** and no channel
number is reported at all.

| Gate | Requirement |
|---|---|
| **G1 expression reimplementation** | The `within_sample_gene_rank` scorer, run on TCGA's own raw RNA table, must reproduce the frozen `frozen_rna_targets.npz` columns at Pearson r >= 0.999. Targets that do not reproduce are **dropped from both cohorts** and named. If fewer than 60 of the 90 survive, stop. |
| **G2 gene coverage** | Each retained signature must reach >= 0.95 gene coverage in the ALCHEMIST expression matrix — the `minimum_required_coverage` the frozen artifact used. |
| **G3 cohort-classifier control** | The within-TCGA control AUC (`cohort_classifier_auc(tcga[:n], tcga[n:2n])`) must land in [0.45, 0.55]. HEST measured 0.5012. A control away from 0.5 means the sampling is broken, as it was once before when h5py sorted-index reads produced disjoint patient sets and a bogus 0.903. |
| **G4 extraction integrity** | >= 95% of the 1,106 slides must yield 30 patches; every slide's md5 must match GDC's published md5; every embedding block finite and (n, 1536). |

---

## 3. The predeclared penalty, from the dilution curve

The deviation is that ALCHEMIST patches are tissue-level and TCGA's are tumour-restricted.
The dilution experiment measured exactly what that costs. From
`v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2, null-corrected retained channel
against foreign-tissue fraction *d*:

| d | 0.00 | 0.10 | 0.20 | 0.30 | 0.40 | 0.60 | 0.80 |
|---|---|---|---|---|---|---|---|
| retained | 1.000 | 0.999 | 0.968 | 0.905 | 0.804 | 0.607 | 0.333 |

Half-loss at d ~ 0.68.

**Predeclared expectation.** A whole-slide tissue mask over an FFPE NSCLC resection admits
adjacent normal lung, stroma, vessels and immune aggregates. I commit in advance to a
plausible band of **d in [0.20, 0.50]**, giving a predicted retained ratio of
**0.97 down to ~0.71** (0.71 by linear interpolation between the measured 0.804 at d=0.40 and
0.607 at d=0.60). The lower edge of that band, 0.71, is the number the bar in §4 is built
from.

**What this calibration does not cover, and must not be used to excuse.** The dilution curve
was measured *within* TCGA — same scanners, same stain distribution, same fixation protocol,
same barcode space. It prices contamination and nothing else. It does **not** price scanner
and staining shift between two institutions' slide pipelines, nor the ALCHEMIST-specific
sampling of one slide per patient. If the observed ratio falls below the band, "contamination"
is not an available explanation for the excess; the excess is unexplained and will be
reported as unexplained.

---

## 4. The bar

Let

```
R = (ALCH_observed - ALCH_null_median) / (TCGA_NSCLC_observed - TCGA_NSCLC_null_median)
```

with both terms measured by the §1 instrument, and let *p* be the ALCHEMIST permutation
p-value (resolution 1/301 = 0.0033 at 300 permutations).

| Verdict | Condition |
|---|---|
| **REPLICATES** | p < 0.01 **and** R >= 0.60 |
| **ATTENUATED BUT PRESENT** | p < 0.01 **and** 0.30 <= R < 0.60 |
| **FAILS TO REPLICATE** | p >= 0.01 **or** R < 0.30 |
| **UNINTERPRETABLE** | any gate in §2 fails |

**Why 0.60 and not 0.71.** 0.71 is the dilution floor for contamination alone. The extra
margin down to 0.60 is the allowance for cross-institution scanner/stain shift, which is real,
is not calibrated by anything this project has measured, and would otherwise be charged
against a channel that did in fact transfer. Setting the bar at 0.71 would let an unpriced
nuisance term produce a "failure"; setting it below 0.30 would let almost anything produce a
"replication". Both edges are chosen now, with no ALCHEMIST number in hand.

**ATTENUATED BUT PRESENT is a real verdict, not a hedge.** It means the channel is
statistically present but loses more than contamination predicts, and the excess loss is then
an open question — not something to be retro-fitted to a larger assumed *d*. If it lands
there, I will say the excess is unattributed.

---

## 5. The mandatory control, and what it can void

Before any channel number is quoted: `cohort_classifier_auc` (logistic regression, C=1.0,
StandardScaler fit per fold, StratifiedKFold(5, shuffle, seed 0), OOF `roc_auc_score`) on
row-L2-normalised raw 1536-d patch embeddings, 20,000 per cohort, TCGA vs ALCHEMIST, plus the
within-TCGA control. This is `v2/calibra/hest.py:cohort_classifier_auc` unmodified — the same
instrument that measured **TCGA vs HEST at AUC 0.99999**.

Predeclared reading:

- **AUC >= 0.99** — residual batch signal is available to every downstream result. This does
  not by itself invalidate the channel, because the channel is measured *within* each cohort
  separately and never pools them, so a between-cohort direction cannot be used by either
  measurement. It does mean that no claim of the form "the representation is cohort-invariant"
  survives, and it must be printed next to the channel number every time the channel number is
  printed.
- **AUC in [0.90, 0.99)** — same, stated less strongly.
- **AUC < 0.90** — worth remarking on, given HEST's 0.99999.

The within-TCGA control is gate G3. A within-cohort control that is not ~0.5 voids the
cross-cohort number too.

---

## 6. What is not being attempted

- No tumour segmentation model to synthesise polygons. It would put a second model's errors
  inside the result.
- No recomputation of the published TCGA numbers. The pan-cancer 0.5573 stands as published;
  this entry adds a TCGA-NSCLC measurement beside it and does not touch it.
- No edit to `NOTEBOOK.md`, the paper drafts, or `claim_guards.py`. Acquiring the cohort does
  not by itself discharge a blocker; whether it does is a separate, deliberate decision made
  after these numbers exist.

---

## 7. Provenance

- Cohort manifest: `alchemist_paired.csv`, sha256
  `b40c1909a323c12afef50a3358ee4af5a0cc8aa207d50f9694ba6f78b10c1773`;
  expression file-id list sha256
  `be88aafe90cd5e32dab64a80acd7597625c0ac45667b7c27dcfdc4f110a9c5bd`.
- GDC counts, re-derived independently of the scout entry and agreeing with it: slide files
  1,349 / 1,175 cases; STAR-Counts expression files 1,138 / 1,107 cases; **intersection 1,106**.
- Selection rule: one slide per paired case, lexicographically first `file_name` — independent
  of file size and of tissue area. 947 of 1,106 patients have exactly one slide regardless.
- Code: `v2/research/external/build_alchemist_manifest.py`,
  `v2/research/external/extract_alchemist_patches.py`, commit `24e3466`.
