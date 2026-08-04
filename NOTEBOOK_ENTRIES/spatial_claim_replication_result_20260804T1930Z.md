## 2026-08-04 19:30 UTC — Five claims re-run at spot level: two survive intact, one strengthens past its own bar, and the confound adjustment is shown to remove slide from the first moment only — a kNN still names the slide 72.9% of the time on the adjusted state

**Logged:** 2026-08-04 19:30 UTC. **Predeclared:**
`NOTEBOOK_ENTRIES/PREDECLARED_spatial_claim_replication_20260804T1800Z.md`, committed `80076b4`,
before any number below existed. **How obtained:** `v2/calibra/hest_claims.py` on the A100 box
(150.136.45.194), workspaces `~/ws_claims{,2,3,4}` deployed by `git archive HEAD` and verified
file-by-file against `git ls-tree -r HEAD` blob SHA-1 (611–623 files, 0 mismatches, every time).
Outputs `/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/out/claims/`.

### Verdicts, worst first

| # | claim | TCGA reference | spot-level | verdict vs the predeclared bar |
|---|---|---|---|---|
| **1a** | adjustment removes the confound | site 0.3633 → 0.0118 | slide 0.9707 → **0.0025** | **holds for the certificate — and the certificate is not enough** |
| **1c** | *(added before the run)* a non-mean probe | not measured | kNN **0.8837 → 0.7291**, global null p95 **0.0856** | **REVERSES the reassurance.** Slide survives adjustment |
| **1b** | adjustment costs the channel nothing | attenuation 0.974–1.039 | **1.030** | holds |
| **1b'** | induced level-0 floor should collapse | 0.08–0.13 at 99 design columns | **0.1299 at 13** | **my prediction FAILS**; the mechanism does not |
| **2** | chance is not zero | 0.140–0.147 at n = 2,766 | **0.1452 / 0.1464** at n = 2,769 | **holds, and the scaling law holds to ≤5% at every n** |
| **3** | random gene sets reproduce the signal | 76–82% | **90–99%**, and **128–147%** across slides | **holds, then overshoots: random sets BEAT the panel** |
| **4** | per-axis certification is insufficient | best axis 0.055, joint 0.235–0.363 | best axis **0.3687**, joint **0.9707** | **weakens on the letter, holds on the substance** |
| **5** | the zero-parameter baseline | 0.5706 vs 0.3565 pooled | CIs **disjoint**; within-slide CI excludes 0 | **holds; my variance bar fails and the reason is instructive** |

### 0. What the run was

Test partition, 13 slides, 53,217 spots, state `wsi_identity` (frozen H-Optimus-0, 1536-d), targets the
50 train-selected genes. **The slide-key fix is verified, not assumed:** `identifier.split("-")[1]`
equals `slide_ids` for all 144,162 rows and yields 44 distinct sites. Subsampling is the predeclared
slide-stratified draw; the anchor is **m = 213 spots/slide, n = 2,769**, matched to TCGA's n = 2,766.
Confound design is `cancer + pooled TSS(=slide)`: **22 columns, rank 13**, against TCGA's 99–108 at
n = 2,530.

### 1. Claim 2 — the spatial permutation null. Holds, and the law is now measured.

`calibration.permutation_null`, k = 16, seed 42, two strata arms.

| n | design rank | perms | predicted 2√(k/n_eff) | null median (strata = cancer) | null p95 | null median (strata = slide) | measured / predicted |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 520 | 13 | 100 | 0.3553 | **0.3423** | 0.3862 | 0.3441 | 0.963 |
| 1,300 | 13 | 100 | 0.2230 | **0.2133** | 0.2403 | 0.2180 | 0.956 |
| **2,769** | 13 | 100 | 0.1524 | **0.1452** | 0.1649 | 0.1464 | 0.953 |
| 5,200 | 13 | 100 | 0.1111 | **0.1057** | 0.1199 | 0.1084 | 0.951 |
| 10,400 | 13 | 100 | 0.0785 | **0.0745** | 0.0865 | 0.0762 | 0.949 |
| 20,800 | 13 | 50 | 0.0555 | **0.0538** | 0.0608 | 0.0535 | 0.969 |
| **53,217** | 13 | 100 | 0.0347 | **0.0333** | 0.0364 | 0.0334 | 0.960 |

At the matched anchor the spatial null is **0.1452**, inside TCGA's measured 0.140–0.147 band. The
predeclared law **2·√(k/n_eff)** tracks the measurement to within **3.1–5.1% at every n over a 100-fold
range**, always slightly above it. Strata choice moves the null by ≤ 0.003 — the cancer/slide
distinction is immaterial for this statistic here.

**This is the transferable number and it transfers as a formula, not as a value.** The capacity floor
is a function of (n, k) and essentially nothing else. Consequences, stated in both directions as
predeclared:

* At the 53,217-spot test partition the floor is **0.033**; extrapolating the same law to the full
  144,162-spot cohort gives **≈0.021**. A spatial channel graded against a 0.14 borrowed from bulk
  TCGA would be judged against a floor **4–7× too high**, and a real channel of r ≈ 0.05 would be
  thrown away as noise.
* Conversely, any spot-level paper quoting a top-CCA-like statistic at small n and calling it
  "significant" against zero is quoting capacity. At n = 520 chance alone is **0.34**.
* The observed adjusted top-CCA is 0.81–0.84 at every n, i.e. 5–25× the null. That is not in dispute;
  what claims 1c, 3 and 4 dispute is what it is *made of*.

### 2. Claim 1 — the adjustment removes slide from the mean, and only from the mean

**1a, the certificate's own verdict (`confound_certificate.certify_axes`, n = 2,769, 13 slide classes,
chance 1/13 = 0.0769, 200 within-cancer permutations, n_boot 100 — the permutation count is reduced
from the 1,000 the TCGA run used and that reduction was declared in advance):**

| arm | per-axis max | per-axis median | per-axis null p95 (median) | axes breaching | joint LDA | joint null p95 | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| raw | 0.3687 | 0.1603 | 0.1452 | **1,212 / 1,536** | **0.9707** | 0.5465 | NOT CERTIFIED |
| adjusted | 0.0744 | 0.0639 | 0.0845 | **0 / 1,536** | **0.0025** | 0.0470 | CERTIFIED |

Joint slide accuracy falls **0.9707 → 0.0025**, a factor of 384, to *far below* the 0.0769 chance rate —
the same below-chance behaviour cross-fitted residualisation produced on TCGA. By the instrument this
project uses, the slide confound is discharged completely. **1a holds.**

**1c, the probe added before the run because the classifier family is mean-based and the adjustment is
a mean.** Out-of-fold 15-NN balanced accuracy for slide, same rows, same design, same folds:

| arm | kNN balanced accuracy | × chance | global-permutation null p95 | permutation p |
|---|---:|---:|---:|---:|
| raw | 0.8837 | 11.5× | 0.0878 | 0.005 (floor) |
| **adjusted** | **0.7291** | **9.5×** | **0.0856** | **0.005 (floor)** |
| adjusted, per-axis standardised (the certificate's own scaling) | 0.7403 | 9.6× | 0.0863 | 0.005 (floor) |

**After the exact adjustment CALIBRA applies before every channel number on this project, a
nearest-neighbour vote still names which of thirteen slides a spot came from 72.9% of the time.** The
global-permutation null for the same classifier on the same features is 0.086, so this is not a
property of kNN being powerful; it is slide information that survived. The mechanism is not subtle and
is not a bug: cross-fitted ridge on a one-hot slide design removes the slide *mean vector*. Spots from
one slide remain each other's nearest neighbours in every other respect — scanner, fixation, tissue
architecture, and the fact that a 128 µm window on a 100 µm pitch physically overlaps its neighbours.

**What this costs, precisely.** It does not invalidate the adjusted top-CCA numbers as *measurements*.
It invalidates the sentence "the site signal is gone", which
`t13_adjusted_certificate_and_p6_20260803T0300Z.md` states and `P1_CALIBRA_DRAFT.md` §4.2 carries. The
supportable sentence is: *the confound is removed from the first moment, and a mean-based certificate
therefore certifies; a classifier that reads anything else still recovers it.* The TCGA site
certificate has never been asked this question, and it should be — the argument is about the
classifier family, not about spatial data, and it applies verbatim to the 85-class TSS result.

**A second, separable defect found while grading 1c: the certificate's null convention is degenerate
in this cohort.** `within_stratum_permutations` permutes slide labels inside cancer type. Five of the
thirteen test slides (TENX70/COADREAD, TENX65/HGSOC, TENX68/IDC, TENX46/PRAD, ZEN49/READ) are the only
slide of their oncotree label, so for those spots the permutation is the identity and the "null" is
handed the true labels. That is why the raw joint null median is **0.533** and the raw kNN
within-cancer null p95 is **0.541**, against a chance rate of 0.0769. The certificate still decides
correctly here because the observed values are far from the null on both sides, but **the within-cancer
null must not be quoted as a chance rate whenever the confound is nested inside the stratum**, and the
global-permutation null (0.086, ≈ 1/13) is the one that means what it says. `hest_claims` now computes
both and reports the unpermuted fraction (5/13 = 0.385) beside them.

**1b — what the adjustment costs the channel.** `spike_recovery_curve`, 25 draws, k = 16, n = 2,769:

* `attenuation_slope` = **1.0300**, inside the predeclared [0.90, 1.10] and inside TCGA's 0.974–1.039.
  **Holds.** The adjustment does not destroy transmitted signal.
* `transmission_floor` = 0.01 (finest level on the grid), `detection_floor` = 0.40 — TCGA read 0.01 and
  0.30/0.40. Unchanged.
* Unadjusted top-CCA 0.9642 → adjusted 0.8105, held-out 0.7935.
* `observed_matched_direction` = 0.0246, far below the 0.40 detection floor, so `observed_above_floor`
  is False — exactly as on TCGA, and for the same reason: the channel is concentrated in particular
  directions, and a random direction pair sees nothing.

**1b' — my own predeclared prediction failed, and the failure is worth more than the prediction.** I
predicted the confound-induced level-0 baseline would **collapse to ≤ 0.03** at spot scale, reasoning
from `calibration.py`'s docstring that it depends on how much of a random (u, v) pair lies in the
design span, and this design spans 13 dimensions of 2,769 rows against TCGA's 99 of 2,530. Measured:
**0.1299** — the top of TCGA's 0.08–0.13 range, at one seventh the design rank. **The prediction as
written is falsified.** The mechanism survives and the correct statement is the one Track 2's alpha
sweep already reached: the induced correlation scales with *how much the nuisance model actually
removes*, not with its rank. Slide removes a great deal here (see claim 5), so a rank-13 design induces
what a rank-99 design induced on TCGA. Rank is not a usable proxy and should not be used as one in the
paper.

### 3. Claim 4 — per-axis certification. Weakens on the letter, holds on the substance.

| cohort | chance | best single axis | × chance | joint | × chance | joint ÷ per-axis (in × chance) | best axis ÷ joint |
|---|---:|---:|---:|---:|---:|---:|---:|
| TCGA d2_h wsi_biology | 0.0118 | 0.0532 | 4.5× | 0.3633 | 30.8× | **6.8×** | 0.15 |
| TCGA d2_i wsi_biology | 0.0118 | 0.0511 | 4.3× | 0.2348 | 19.9× | **4.6×** | 0.22 |
| **HEST spatial, slide** | 0.0769 | **0.3687** | **4.8×** | **0.9707** | **12.6×** | **2.6×** | **0.38** |

The predeclared bar was "joint/chance exceeds per-axis-max/chance by ≥ 3×". It is **2.6×**, so on the
letter the claim **weakens**. But the bar is the wrong instrument here and saying so is not a rescue:
joint accuracy is **0.9707 against a ceiling of 1.0**, so the ratio is compressed by saturation, not by
the leak being concentrated. The predeclared *reversal* condition — best axis ≥ 0.8 of joint — is
nowhere near met (0.38). And 1,212 of 1,536 axes individually breach their own null. The substantive
claim is intact: **no single axis carries the leak, and per-axis certification of these 1,536 axes
would have passed 324 of them while a joint test reads 0.97.** The transferable lesson is that
"multiples of chance" is not a safe way to compare certificates across class counts once the joint test
saturates; report the raw pair and the ceiling.

### 4. Claim 3 — random gene sets. Holds, then overshoots.

Three arms, all on the same residualised rows, all scored by
`run_calibra.score_target_block_per_column`.

**(a) The shipped `RANDOM_CONTROL__` columns are a pairing null, not a random gene set.** Fitted-direction
median **−0.0021** against the real panel's 0.4816. Predeclared to read ≈ 0; it does. Reported so the
two controls are never conflated.

**(b) Genuine matched random gene sets, within-fold statistic (the literal port of TCGA's T1.4).** Five
independent 50-gene sets drawn from the 17,197-gene cross-slide intersection, each gene matched to a
panel gene on train-spot (mean, log-variance), zero overlap with the panel:

| | real panel | set 0 | set 1 | set 2 | set 3 | set 4 |
|---|---:|---:|---:|---:|---:|---:|
| fitted-direction median | 0.4816 | 0.4765 | 0.4575 | 0.4335 | 0.4606 | 0.4612 |
| ratio to real | — | **0.990** | **0.950** | **0.900** | **0.957** | **0.958** |

Median ratio **0.957** against TCGA's 76–82%. Predeclared "holds if ≥ 0.65"; it holds with room.

**(c) The arm that actually generalises, and the one that reverses the reading.** `score_target_block_per_column`
folds spots at random, so its train and test folds share slides and — at Visium density — share
overlapping 128 µm windows. So the same question was asked where no fold can share tissue: ridge
(α = 1000) fit on the **22 training slides**, scored by `within_slide_r` on the **13 test slides**.

| | real panel | set 0 | set 1 | set 2 | set 3 | set 4 |
|---|---:|---:|---:|---:|---:|---:|
| within-slide r | **0.1506** | 0.2219 | 0.1984 | 0.1929 | 0.1936 | 0.2020 |
| ratio to real | — | **1.473** | **1.317** | **1.281** | **1.285** | **1.341** |
| pooled r | 0.3565 | 0.4521 | 0.4311 | 0.4117 | 0.4398 | 0.4296 |

**Random matched gene sets are predicted 1.28–1.47× better than the curated 50-gene panel, on the only
statistic in this project that survives a slide boundary.** Not 76–82% of the signal — 128–147% of it.
Variance matching is good (mean per-column train variance 5.63–5.96 for the random sets against 5.73
for the panel, ≤ 4%), so this is not a variance artefact; the random sets do sit at higher mean
expression after their own-set normalisation (3.18–3.24 vs 2.32), because the panel is the
top-variance/high-zero tail, and that is the honest residual caveat on the direction of the excess.

**The reading.** At spot level there is no evidence of *any* gene-set-specific component in the
morphology→expression channel. Whatever H-Optimus-0 reads off a 112 µm field predicts a random
variance-matched set of 50 genes at least as well as the 50 most variable ones. TCGA's 76–82% was
already the strongest single argument that this channel is non-specific; spatially it is 96% by the
same statistic and >100% by a stricter one. **Bulk averaging was, if anything, hiding how
non-specific the channel is.**

### 5. Claim 5 — the zero-parameter baseline, with intervals. Holds; my variance bar fails; the
identity behind it is exact.

Bootstrap unit is the **slide** (13 of them), 2,000 replicates, because spots inside a slide are
neither independent nor non-overlapping.

| predictor | parameters | pooled r | pooled 95% CI | within-slide r | within-slide 95% CI |
|---|---:|---:|---|---:|---|
| per-slide mean | **0** | **0.5706** | **[0.4652, 0.6032]** | −7e−18 | [−3e−17, 1e−17] |
| global (train) mean | 0 | 0.0000 | [−3e−17, 3e−17] | 7e−18 | [−9e−18, 2e−17] |
| ridge on H-Optimus-0 | 1536×50 | 0.3565 | **[0.2695, 0.4208]** | **0.1506** | **[0.1211, 0.1831]** |

* The per-slide-mean pooled CI lies **entirely above** the ridge pooled CI (0.4652 > 0.4208): the
  zero-parameter baseline's win on the metric the HEST leaderboard reports is not a point estimate,
  it survives a slide bootstrap. **Holds.**
* The ridge within-slide CI **excludes 0** (0.1211–0.1831). There is a real, slide-identity-free
  morphology→expression channel. **Holds** — and this is the condition whose failure would have made
  everything above moot.
* **My predeclared variance bar fails.** Between-slide share of total test-spot target variance is
  **median 0.357** (p10 0.100, p90 0.621), not > 0.5.

**The reason the bar was wrong is the most quotable thing in claim 5.** The pooled Pearson r of a
per-slide-mean predictor is *exactly* the square root of the between-slide variance share — measured:
median-gene share 0.3568, median-gene pooled r 0.5973, and √0.3568 = 0.5973 (now pinned as a test).
Correlation is the square root of a variance ratio, so **a slide effect explaining barely a third of
the variance already buys a pooled r of 0.60**. "Pooled correlation is dominated by between-slide
variance" is true, but the honest phrasing is stronger and less obvious: *a minority of the variance
being between-slide is sufficient to dominate the pooled correlation, because r hides the square.*
Any paper reporting pooled spot-level r should report the between-slide variance share beside it.

### What a full-scale run costs (measured, not guessed)

`permutation_null` re-whitens the image block inside its permutation loop — an SVD of n×1536 per
permutation — which is the whole cost.

| n | s / permutation | workers |
|---:|---:|---:|
| 2,769 | 0.41 | 12 |
| 10,400 | 0.69 | 12 |
| 20,800 | 1.43 | 12 |
| 53,217 | 6.48 | 6 (≈3.2 at 12) |

Linear in n above ~5,000. Extrapolated to the **full 144,162-spot cohort at 12 workers: ≈8.8 s per
permutation, so 1,000 permutations ≈ 2.4 h** — and this was measured on a box carrying 13 co-tenant
processes at load 30+, so a quiet box is faster. The 2026-08-03 note that a spot-level CALIBRA sweep is
"a deliberate multi-hour job" is correct in magnitude but the blocker is not the null: it is
**`certify_axes`**, whose joint LDA is O(n·d²) per permutation at d = 1536. At n = 2,769 the raw arm
cost 368 s for 200 permutations; at n = 53,217 the same arm would cost **≈2 h at 200 permutations and
≈10 h at the 1,000 the TCGA certificate used**. A full-cohort certificate is the thing to budget for,
and a numerically identical hoist of the x-side whitening out of `permutation_null`'s loop would cut
the null's cost by roughly the permutation count — worth doing before anyone runs this at 144k.

### Honest constraints on every number above

* **44 slides, 18 bowel; 13 test slides over 7 labels, 5 of them colorectal.** This measures a channel.
  It cannot support a cross-cancer generality claim, and the slide/cancer nesting is what made the
  certificate's null degenerate.
* **The image window is not the assay.** 128 µm cut, **112 µm** seen after H-Optimus-0's own 87.5%
  centre crop, against a 55 µm Visium spot: **5.28× more tissue than the transcriptome came from**. The
  target is not localised to the window. This dilutes every association reported here rather than
  manufacturing one.
* **Licence** CC BY-NC-SA 4.0, NoDerivatives samples excluded and staying excluded.
* **Cohorts do not mix**: TCGA vs HEST embedding AUC 0.99999 against a within-TCGA control of 0.5012.
  No number here is transferred to TCGA except as a labelled comparison.
* **Subsampling is declared**: the anchor is 213 spots/slide, n = 2,769; the null sweep spans
  520–53,217; the certificate and the channel arms are anchor-only, and the reasons are cost figures
  in the table above.

### Prose implications, flagged rather than edited

Not touching `NOTEBOOK.md` or the drafts. Three things in the current text are now wrong or
unsupportable as written:

1. **`P1_CALIBRA_DRAFT.md` §4.2 and `t13_..._20260803T0300Z.md`: "the site signal is gone."** Supported
   only for mean-based classifiers. Spatially, the same adjustment leaves a 9.5×-chance kNN signal. The
   TCGA site certificate needs the same probe before that sentence can stand.
2. **`calibration.py`'s docstring attributes the induced level-0 floor to how much of (u, v) lies in
   the design span**, which invites reading it as rank-driven. A rank-13 design induced 0.1299 where a
   rank-99 design induced 0.08–0.13. Explanatory power, not rank.
3. **The 76–82% random-gene-set figure should not be quoted as an upper bound on non-specificity.**
   Spatially it is 96% within-fold and 128–147% across slides.

One addition rather than a correction: **the capacity floor should be quoted as 2·√(k/n), with the
measured 0.95–0.97 calibration factor**, not as the scalar 0.140/0.147. The scalar is only correct at
n ≈ 2,766, k = 16, and every spatial paper in this literature runs at 10–100× that n.

### Files / commits

* `v2/calibra/hest_claims.py`, `v2/tests/test_hest_claims.py` (18 tests) — commits `9785dbb`,
  `99ac306`-fix, and the four that follow on `research/rebase-vision`.
* Predeclaration `80076b4`.
* Results: `/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/out/claims/`
  `{claim1_claim4_certificate,claim1b_channel,claim1c_knn_probe,claim3_gene_sets,claim5_baselines}.json`,
  `null_small/claim2_null_sweep.json`, `null_full/claim2_null_sweep.json`,
  `genesets_cross/claim3_gene_sets.json`, `random_gene_sets_{test,all}.npz`.
* Logs: `/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/logs/claims_*.log`.

### Suite status

**At the commit this entry lands on (`be0ce4e`): `morpheus/v2/tests` is 408 passed, 0 failed, 27
errors**, the 27 being `test_p2_figures.py` needing matplotlib, absent from `~/venv` by policy. My 18
new tests are inside that 408.

Two transient failures were observed on the way and are recorded because they were real when seen,
both in other people's files and both since resolved by their owners:

* `test_leave_sites_out.py::test_indexed_cca_recovers_a_planted_signal` — asserted
  `heldout_top_cca_indexed(...) > 0.8`, measured **0.7923**, at the commit I started from. Passes now.
* `test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree` — flagged
  `v2/tests/test_p2_centring_amplification.py` (added by another agent at `88e37ff`) as SVD-based rank
  outside `calibra/spectral.py`. Fixed by that agent at `7c435dd`.

Neither was touched by me. The guard test was also re-run with my two files removed from the tree, to
confirm the offender list never contained them: 13 passed.

**Note on the shared checkout.** This repository is being edited concurrently by other agents; the
tracked file count moved 611 → 637 during this session and `git ls-tree -r HEAD` grew under commits
that are not mine. Every workspace deploy was a fresh `git archive HEAD` verified file-by-file against
that tree's blob SHA-1 (611 / 618 / 621 / 623 / 637 files, zero mismatches), so every number above was
produced from a tree provably identical to a named commit.
