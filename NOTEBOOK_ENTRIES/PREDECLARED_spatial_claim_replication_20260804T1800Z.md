## 2026-08-04 18:00 UTC — PREDECLARATION: re-running five project claims at spot-level spatial resolution, with the reading of each written down before any run

**Logged:** 2026-08-04 18:00 UTC, **before** any statistic on this page has been computed. **Cohort:**
HEST-1k spatial artifact `hest_spatial_hoptimus.npz` / `hest_spatial_targets.npz`,
`/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/out/artifacts/`. **Nothing in this file is a
result.** The numbered readings below are falsifiers; the run that follows either meets them or does
not, and either way the outcome is reported against *this* text.

### 0. What is fixed before the run

**Unit, partition, state.** A row is a Visium spot. Analysis partition is `test`: **13 slides,
53,217 spots**, slide sizes 1,653–9,079, cancer labels COAD 21,134 / ORGAN_PROSTATE 11,436 /
COADREAD 9,079 / HGSOC 4,669 / PRAD 3,043 / READ 2,203 / IDC 1,653. Representation state is
`wsi_identity` — row-L2-normalised frozen H-Optimus-0, 1536-d. Targets are the 50 train-selected
variable genes, log1p CPM(1e4), plus the 16 shipped `RANDOM_CONTROL__` permuted columns.

**The slide-key fix is verified, not assumed.** `residualise.pooled_tissue_source_site` derives the
site confound as `identifier.split("-")[1]`. On the built artifact, that field equals `slide_ids`
for **all 144,162 rows**, yielding **44 distinct sites** (13 in the test partition). Had the keys
been naive spot ids the field would have been the constant `"1"` from the barcode suffix and the
site term would have silently vanished from a design still reported as slide-adjusted. Checked
first, because every claim below rests on it.

**Subsample rule — declared, not silent.** CALIBRA's permutation null re-whitens a 53,217×1536 block
once per permutation; the 20-minute abort logged on 2026-08-03 was that cost. Subsampling is
therefore **slide-stratified**: for a per-slide cap `m`, take a seed-42 random draw of
`min(m, n_slide)` spots from each of the 13 test slides. Slide balance is preserved because slide is
the confound under test. The grid is

| m (spots/slide) | 40 | 100 | **213** | 400 | 800 | 1600 | ALL |
|---|---:|---:|---:|---:|---:|---:|---:|
| n | 520 | 1,300 | **2,769** | 5,200 | 10,400 | 20,800 | 53,217 |

**n = 2,769 (m = 213) is the anchor**, chosen to match TCGA's n = 2,766 as closely as 13-way slide
balance allows, so that the matched-n comparison isolates the *data* rather than the sample size.
Every reported number names its own n. Nothing is truncated without appearing in this table.

**Confound design.** `confound_design(DataFrame({"cancer", "tss"}), ["cancer", "tss"])` with `tss`
from `pooled_tissue_source_site(patient_ids, min_site_count=10)` — i.e. **slide**. All 13 test slides
clear the pooling threshold. Slide nests inside cancer, so the design has ~13 free columns against
TCGA's 99–108 at n = 2,530: **a much lower design-rank-to-n ratio, which is itself predictive** (see
claim 1's second reading).

**Nothing is recomputed inline.** Every rank, channel, null, floor and certificate statistic is the
imported `v2/calibra/` function: `calibration.permutation_null`, `calibration.spike_recovery_curve`,
`spectral.{cca_spectrum, top_canonical_correlation, heldout_top_cca, effective_rank}`,
`confound_certificate.certify_axes`, `run_calibra.{score_target_block_per_column,
random_direction_column_correlation, grade_random_controls}`, `residualise.*`, `hest.{pooled_r,
within_slide_r, per_slide_mean_baseline, normalise_expression}`. The one new statistic (claim 1c) is
added to `v2/calibra/` with tests, not written into a run script.

### 1. Confound adjustment works two-sided — *the spatial analogue of "site" is "slide"*

TCGA reference: cancer-type balanced accuracy **0.463 → 0.035** (chance 0.048, n = 2,530); joint-LDA
site accuracy **0.3633 → 0.0118** (chance 1/85 = 0.0118, n = 2,766); spike-recovery attenuation slope
**0.974–1.039**.

Spatial instrument: `certify_axes(..., residualise=False)` and `(..., residualise=True)` on the
anchor subsample, 13 slide classes, **chance = 1/13 = 0.0769**, within-cancer label permutations,
`n_permutations = 200` (reduced from 1,000; the cost is the joint LDA at d = 1536 and the reduction
is declared here, not discovered later). Attenuation from `spike_recovery_curve` on the same rows.

**1a — removal.** *Holds* if adjusted joint-LDA slide accuracy falls to ≤ its own permutation null
p95 and within 1.5× of chance (≤ 0.115). *Weakens* if it lands between 1.5× and 3× chance.
*Reverses* if the adjusted state still recovers slide at > 3× chance (> 0.23), because that would
mean the adjustment CALIBRA applies before every channel number on this project does not discharge
the dominant batch factor when that factor is a slide.

**1b — cost to the channel.** *Holds* if `attenuation_slope` ∈ [0.90, 1.10] and the level-0
confound-induced baseline (`baseline_recovered_median`) is **below** TCGA's 0.08–0.13.
Predeclared direction and reason: the induced correlation scales with how much of a random (u, v)
pair lies in the design span, and this design spans ~13 dimensions of 2,769 rows against TCGA's ~99
of 2,530. **The prediction is that the induced floor collapses at spot scale**, ≤ 0.03 at the anchor
and smaller still at full n. If instead it is *larger*, the mechanism written into
`calibration.py`'s docstring is wrong and that must be said.

**1c — the classifier family, declared as a limitation before it is measured.** `certify_axes` uses
nearest-class-mean per axis and shrunk-covariance LDA jointly. Both are **mean-based**, and
cross-fitted residualisation on a one-hot slide design removes slide *means* by construction. A pass
is therefore close to arithmetic and is weak evidence that "the batch signal is gone" — the same
objection applies verbatim to the TCGA site certificate this claim comes from. A second,
non-mean-based probe is therefore run on the identical rows: **out-of-fold 15-nearest-neighbour
balanced accuracy for slide**, raw and adjusted, against the same 1/13 chance rate. *Reading:* if
kNN recovers slide at ≫ chance from the adjusted residuals, then "adjustment discharges the
confound" is a statement about linear-mean classifiers and not about the representation, and both
this claim and its TCGA parent must be re-scoped in the paper. This probe is being added because the
weakness was noticed while reading the code, before the run, and recording it afterwards would have
looked like a rescue.

### 2. Chance is not zero — the spatial permutation null *(the transferable number)*

TCGA reference: 16-component top-CCA permutation null median **0.140–0.147** at n = 2,766 (three
procedures; `p2_rewritten_around_the_surviving_claim` records that 0.140 is a row-shuffle and
0.145–0.147 a within-cancer permutation, and the two must not be quoted interchangeably).

Spatial instrument: `calibration.permutation_null(x, y, design, strata=…, n_permutations=100,
n_components=16, seed=42)` at every n on the grid, in **two strata arms** — `strata=cancer` (the
literal port of the TCGA convention) and `strata=slide` (the spatial analogue of the site stratum).
Both are reported; neither is designated primary after the fact.

**Predeclared quantitative expectation, written as a formula so it can be wrong.** The top canonical
correlation between two independent PCA-whitened k-dimensional blocks at sample size n concentrates
near the Wachter/Marchenko–Pastur edge **2·√(k/n)**. At k = 16 this predicts

| n | 520 | 1,300 | 2,769 | 5,200 | 10,400 | 20,800 | 53,217 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2·√(16/n) | 0.351 | 0.222 | **0.152** | 0.111 | 0.078 | 0.055 | **0.035** |

At the anchor this predicts **0.152** against TCGA's measured 0.140–0.147 — within 4–8%. *Holds* if
the measured spatial null at n = 2,769 lies in **[0.12, 0.18]** and the sweep tracks 2·√(k/n) to
within 20% relative at every n. *Weakens* if the anchor is in range but the scaling departs by more
than 20% at the extremes. *Reverses* if the anchor null falls outside [0.12, 0.18], because then the
capacity floor is a property of the data's covariance structure and not of (n, k), and no floor can
be carried between cohorts at all.

**The consequence is stated in advance in both directions.** If the law holds, then the floor at the
full 144,162-spot cohort is **≈0.021** and at the 53,217-spot test partition **≈0.035** — roughly
**4–7× smaller than the 0.14 this project quotes from bulk TCGA**. Any spatial channel graded
against a borrowed 0.14 would be judged against a floor several times too high, and a genuine small
channel would be discarded as noise. If the law fails, the opposite warning applies and the number
is not transferable at all.

### 3. Random gene sets reproduce most of the per-target signal

TCGA reference: covariate-matched random gene sets read at **76–82%** of the level real curated sets
read at, in every state, invariant to dilution (fitted-direction per-column statistic).

Two arms, because the shipped control is **not** the same object as TCGA's:

**3a — permuted-column arm (what the artifact ships).** The 16 `RANDOM_CONTROL__` columns are
row-permutations of real target columns: marginals preserved, spot correspondence destroyed. This is
a *pairing* null, not a random gene set. **Predeclared: this arm should read ≈ 0** (fitted-direction
median within ±0.02 of zero) and a non-zero reading would indicate leakage in the residualisation or
the fold structure, not a gene-set effect. It is reported so that the two controls are not confused.

**3b — genuine random-gene-set arm.** Draw **5 independent 50-gene sets** from the 17,197-gene
cross-slide intersection **excluding the 50 panel genes**, each gene matched to a panel gene by
nearest neighbour in (train-spot mean, train-spot log-variance), seeds 101–105. Expression rebuilt
from the raw `st/*.h5ad` for the 13 test slides through `hest.normalise_expression` — the same
transform the panel got. Score both blocks with `run_calibra.score_target_block_per_column` on the
identical residualised rows, and report the ratio of median fitted-direction correlation.

*Holds* if the random-set / real-set ratio is **≥ 0.65** and its spread across the 5 draws is within
±0.10. *Weakens* if 0.40–0.65. *Reverses* if **< 0.40**, i.e. localising the target to the tissue
recovers a pathway-specific signal that bulk averaging had destroyed — which is the single most
interesting way this claim could fail, and the reason the arm is worth building.

### 4. Per-axis certification is insufficient — the leak is smeared

TCGA reference: raw per-axis maximum **0.0506–0.0551** against chance 0.0118 (4.3–4.7× chance) while
joint LDA reaches **0.2348–0.3633** (20–31× chance).

Spatial instrument: the same `certify_axes` run as claim 1, raw arm, chance 1/13 = 0.0769. Reported
as **multiples of chance** so that the 13-class and 85-class certificates are comparable at all.

*Holds* if joint/chance exceeds per-axis-max/chance by **≥ 3×**. *Weakens* at 1.5–3×. *Reverses* if
the best single axis reaches ≥ 0.8 of the joint accuracy, i.e. the slide leak is concentrated on one
coordinate and per-axis certification would in fact have caught it. Prior expectation, recorded:
TCGA↔HEST cohort AUC is 0.99999, so slide/batch structure in this representation is large; the joint
accuracy is expected near ceiling, and the informative quantity is therefore the **per-axis
maximum**, not the joint.

### 5. The zero-parameter baseline, with intervals

Reference (point estimates only, no intervals): per-slide mean pooled r **0.5706** / within-slide
**0.0000**; ridge on H-Optimus-0 **0.3565** / **0.1506**; global mean 0 / 0.

Instrument: `hest.per_slide_mean_baseline`, `hest.pooled_r`, `hest.within_slide_r`, ridge refit
train→test at the same alpha = 1000. Intervals from a **2,000-replicate bootstrap over the 13 test
slides** — the slide is the resampling unit because spots inside a slide are neither independent nor
non-overlapping. Plus a variance decomposition: per-gene between-slide share of total test-spot
variance.

*Holds* if the per-slide-mean pooled 95% CI lies entirely above the ridge pooled 95% CI, the ridge
within-slide CI lies entirely above 0, and the between-slide variance share is > 0.5 for the median
gene. *Weakens* if the pooled CIs overlap (13 slides is a small bootstrap universe and this is a
real possibility). *Reverses* if the ridge within-slide CI includes 0, in which case there is no
demonstrated morphology→expression channel at spot level at all and everything above is moot.

### Constraints that ride with every number below, without exception

* **44 slides, 18 of them bowel; 13 test slides over 7 labels, 4 of them colorectal.** This cohort
  can measure a channel. It cannot support any cross-cancer generality claim.
* **The image window is not the assay.** A 128 µm field is 6.90× the tissue a 55 µm Visium spot
  assayed; after H-Optimus-0's own 87.5% centre crop the encoder sees 112 µm, **5.28×**. The target
  is not localised to the window. Every spatial claim carries this.
* **Licence.** CC BY-NC-SA 4.0 collection; the 41 CC BY-NC-ND samples are excluded and must stay
  excluded from any released derivative.
* **Cohorts do not mix.** TCGA vs HEST embedding AUC 0.99999 against a within-TCGA control of
  0.5012. No number here is transferred to or from TCGA except as an explicitly labelled comparison.

### Cost, declared in advance

The full-scale run is not attempted and its cost is estimated instead. `permutation_null` calls
`cca_spectrum`, which re-whitens the image block **inside** the permutation loop — an SVD of
n×1536 per permutation. That is the multi-hour term, and it is O(n·1536²) per permutation. The
measured per-permutation cost at each n on the grid is reported, and the full-cohort (144,162-spot,
1,000-permutation) figure extrapolated from it rather than guessed.
