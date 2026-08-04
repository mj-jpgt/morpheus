## 2026-08-04 18:40 UTC — PREDECLARATION: asking the TCGA confound certificate the question the spatial replication asked, with a classifier that is not a function of class means

**Logged:** 2026-08-04 18:40 UTC, **before** any probe statistic in this file exists. **Nothing below
is a result.** Every number that appears here is either (a) a fixed cohort fact measured before the
predeclaration and labelled as such, or (b) a threshold I am committing to in advance. The run that
follows is graded against *this* text and against nothing written afterwards.

---

### 0. Why this run exists

`P1_CALIBRA_DRAFT.md` §4.2 says the confound adjustment is *verified, not assumed*, on two
confounders: cancer-type balanced accuracy **0.463 → 0.035** against chance **0.048**, and joint-LDA
pooled-tissue-source-site accuracy **0.3633 → 0.0118** against chance **0.0118**, with **zero**
breaching axes in all six artifact × state blocks. `t13_adjusted_certificate_and_p6_20260803T0300Z.md`
states the conclusion as *"the site signal is gone."*

`v2/calibra/confound_certificate.py` scores that claim with exactly two classifiers:
`nearest_class_mean_oof` (per axis) and `lda_oof_balanced_accuracy` (jointly). **Both are functions of
the class means.** The adjustment being certified is `residualise.cross_fitted_residuals` — a ridge
regression of the representation on a one-hot `cancer + pooled TSS` design, which by construction
removes the class *mean vector*. Passing a mean-based test after a mean-removing adjustment is close
to arithmetic.

The spatial replication (`spatial_claim_replication_result_20260804T1930Z.md`, claim 1c) showed this
is not a hypothetical. On HEST, an adjustment that took joint-LDA slide accuracy from 0.9707 to
**0.0025** (chance 0.0769) left an out-of-fold 15-NN naming the slide **72.9%** of the time, against a
global-permutation null p95 of 0.086 and a permutation *p* at the 1/201 floor.

**The argument there was about the estimator family, not about spatial data. TCGA has never been
asked.** This run asks it.

---

### 1. Cohort facts, measured before this predeclaration was written

These are structural properties of the artifacts, not outcomes. They are recorded here because the
thresholds below depend on them and because the null convention turns on the nesting.

Artifacts `~/e0_run/d2_v3/d2_v3_s4{2,3,4}/artifacts/d2_{h,i}_seed4{2,3,4}.npz` and D1-B's
`~/e0_run/d1_v2/artifacts/d1_{p,f}_seed4{2,3,4}.npz`, all 6,427 patients × 256 dims, states
`{wsi_biology, rna_biology, full_biology}`, splits train 3,118 / val 543 / test 2,766.

| partition | n | cancer classes | raw TSS codes | pooled TSS classes (`min_site_count=10`) |
|---|---:|---:|---:|---:|
| **test** (the §4.2 partition) | **2,766** | **21** (chance 1/21 = **0.047619**) | 352 | **85** (84 kept + `OTHER`; chance 1/85 = **0.011765**) |
| all | 6,427 | 32 (chance 0.031250) | 610 | 186 |

**The nesting is total, and the direction of the trap is not the one HEST had.** On the test
partition **0 of the 84 kept sites contributes patients to more than one cancer type** — a TSS code
determines the cancer. Only the pooled `OTHER` class (829 of 2,766 rows, 30.0%) spans all 21 cancers.

Consequently:

* The HEST failure mode — *a stratum containing a single class, so within-stratum permutation is the
  identity* — is nearly absent here: **1 of 21 test cancers has only one pooled site**, covering
  **28 of 2,766 rows (1.0%)**. Not 5 of 13 slides / 38.5% as on HEST.
* A **different** consequence of nesting is live and is the one that matters. Because site nests
  inside cancer, a within-cancer permutation hands the null the entire cancer→site restriction: a
  permuted label is always a site *of the correct cancer*. The within-cancer null therefore measures
  *"site beyond cancer"* and is structurally far above 1/85 for any representation that knows cancer.
  **It is not a chance rate and will not be quoted as one.**

I therefore predeclare a third reference alongside the two permutation nulls: **`within_stratum_chance`**,
the balanced accuracy of an oracle that knows the cancer exactly and then guesses uniformly among the
pooled sites present in that cancer — mean over classes of 1/(sites in that class's cancer). It is
computed from the label/stratum tables alone, before any feature is touched, and it is the analytic
value the within-cancer permutation null should sit near if the representation carries cancer and
nothing more.

**Which null applies, declared now.** The claim under test is that the adjusted representation is not
reading site — and the adjustment removes cancer *and* site. A state from which both have been
regressed out should not be able to name the site at all. **The applicable bar for the adjusted arms
is the global-permutation null, i.e. measured chance.** The within-cancer null is the correct bar only
for the narrower question the certificate poses about a *raw* axis ("does it carry site beyond what
cancer already explains?"). Both are computed for every arm; the verdict is read against the global
null; neither is assumed and both are measured by permutation rather than taken as 1/n_classes.

---

### 2. What will be run

**Instrument.** New module `v2/calibra/nonlinear_confound_probe.py`, imported — nothing computed
inline in a shell. Folds, label encoding, balanced accuracy and the permutation machinery are imported
from `confound_certificate.py` so the probe and the certificate share fold assignment and scoring
exactly. The plain k-NN vote is pinned by a test to reproduce
`hest_claims.knn_balanced_accuracy_oof` bit-for-bit on a fixture, so the spatial and bulk readings are
produced by the same estimator.

**Anchor grid** (the blocks §4.2 quotes): `d2_h_seed42` and `d2_i_seed42`, state `wsi_biology`,
partition `test`, n = 2,766.

**Targets:** pooled TSS (85 classes) and cancer type (21 classes). Cancer type is also run on
partition `all` (32 classes) so the 32-class figure quoted in the project's framing is covered.

**Feature arms**, all three, every time:

1. `raw` — the state as stored.
2. `adjusted` — `cross_fitted_residuals(state, confound_design(cancer + pooled TSS))`, seed 42,
   `n_splits=5`, `alpha=1.0`: the *identical* call `certify_axes(..., residualise=True)` makes, so this
   is the same adjusted state the published numbers come from.
3. `adjusted_standardised` — arm 2 with the certificate's own per-axis standardisation, because a
   distance-based probe is scale-sensitive and the certificate standardises before scoring.

**Probe families and why each was chosen.**

* **k-NN, out-of-fold, k ∈ {1, 3, 5, 10, 15, 25, 50}.** The decision rule is the labels of the nearest
  training points. It is a function of local neighbourhood structure and involves no class mean at any
  step, so removing class means cannot make it pass by construction. Sweeping k is not a
  hyperparameter search for the best number: small k reads the finest local structure and is the most
  sensitive to a residual confound, large k smooths toward a density estimate. **The reading is taken
  at the maximum over k**, declared now, so that a favourable number at one k cannot be selected after
  the fact — and so that a probe which only looks clean because k was large is caught.
  A **prior-corrected** variant is run alongside the plain vote: TCGA's site classes are wildly
  imbalanced (`OTHER` is 30% of rows), and a plain majority vote over-predicts the largest class,
  which *depresses* balanced accuracy and biases the run toward the favourable conclusion. Votes are
  weighted by the inverse training-fold class frequency — hyperparameter-free, and the correct rule
  for the balanced-accuracy loss actually being reported.
* **Random forest, out-of-fold, 300 trees, `max_features="sqrt"`, `class_weight="balanced_subsample"`,
  seed 42.** Chosen because it fails in a *different* way from k-NN and shares none of its
  assumptions. Its decision rule is a set of axis-aligned thresholds; it is invariant to any monotone
  rescaling of any single axis, it is not a metric method, and it can key on a difference in variance,
  in skew, or in an interaction between two coordinates — every one of which is invisible to LDA and
  survives mean-removal untouched. If a confound persists as a second-moment or interaction effect,
  the forest is the family that should see it and the certificate's two classifiers are the family
  that cannot.
* **RBF-kernel SVM, out-of-fold, `C=1`, `gamma="scale"`, `class_weight="balanced"`, on standardised
  features.** A third family, smooth and global where k-NN is local and discrete: its decision
  function is a weighted sum of Gaussians centred on training points, so it reads the same local
  geometry as k-NN but fits it globally and is not vulnerable to k-NN's specific failure on small
  classes (vote ties resolved by class size). It is declared **anchor-only**; if measured wall cost
  makes it infeasible it will be reported as not run, and its absence will be a stated cost decision
  rather than a choice made after seeing a number.

**Nulls.** For every (block, target, arm, probe): a **global** permutation null (labels permuted with
a single constant stratum — this *measures* chance and is not assumed to be 1/n_classes) and a
**within-cancer** permutation null (`within_stratum_permutations`, the certificate's own convention).
Both via the same imported function. For the cancer target the within-cancer null does not exist and
only the global null is reported. Null median, p95 and the permutation *p* floored at 1/(P+1) are
reported for each, together with the resolution, so a *p* at the floor is never read as smaller than
the design supports.

**Permutation counts, declared now:** P = 200 for k-NN (matching the spatial run), P = 100 for the
random forest and the SVM (cost). If measured wall cost puts the anchor grid over 2 hours the forest
count may be cut, **to no fewer than 50**, and any cut is reported next to the number it produced. No
count is ever increased after seeing a result.

**Secondary blocks, no thresholds attached, reported for breadth:** the remaining d2_v3 artifacts
(seeds 43, 44; both `d2_h` and `d2_i`), the other two states (`full_biology`, `rna_biology`), and the
six D1-B artifacts `~/e0_run/d1_v2/artifacts/d1_{p,f}_seed4{2,3,4}.npz`. These run k-NN at k ∈ {5, 15}
and the forest, global null only. Every block that does not run is named in the result entry.

---

### 3. The readings — written down before the run

All thresholds refer to the **adjusted** arm, **site** target, **maximum over k**, at the anchor, and
are stated as multiples of the measured global-null median (with the 1/n_classes value in brackets as
the design chance).

| # | condition | reading |
|---|---|---|
| **A — SOUND** | k-NN ≤ **2× chance** (≤ 0.0235 for site, ≤ 0.0952 for cancer) at **every** k, **and** ≤ the global-null p95, **and** the forest likewise, in **both** artifacts | **The adjustment removes the confound in a way that survives a nonlinear reader.** P1 §4.2 stands as written and is now much better evidenced. The spatial 1c finding is a property of the spatial cohort (13 slides, physically overlapping windows), not of the estimator. |
| **B — FIRST MOMENT ONLY** | k-NN ≥ **5× chance** (≥ 0.0588 site, ≥ 0.2381 cancer) at any k, with permutation *p* at the floor against the **global** null | **"The site signal is gone" must become "removed from the first moment".** Correction to `P1_CALIBRA_DRAFT.md` §4.2, `t13_adjusted_certificate_and_p6_20260803T0300Z.md`, and every downstream restatement, at the prominence of the original claim. Every adjusted number on this project inherits the caveat. |
| **C — BETWEEN** | 2× < k-NN < 5× chance | Report the magnitudes for every block and **do not adjudicate**. Name what would resolve it. |
| **D — UNINFORMATIVE** | the **raw** arm's k-NN is itself below 2× chance, **or** the global-null p95 exceeds 2× chance | The probe has no power at this n, d and class count. Nothing is concluded about the adjusted state in either direction, and the run is reported as a failed measurement rather than as a favourable one. |

The same four readings apply to the cancer target with its own chance rate substituted.

---

### 4. What would make me distrust a **favourable** result

This is the section that matters, because outcome A is the comfortable one and the one I have the
least incentive to interrogate. Every item below is a check I commit to running and reporting **even
if the headline comes out clean**.

1. **No power ⇒ no finding.** If the *raw* arm's k-NN is not clearly above chance, the adjusted arm's
   being at chance demonstrates nothing. The raw arm is graded first and reported first, and a clean
   adjusted number sitting behind a weak raw number is reported as reading D, not A.
2. **Below-chance is a symptom, not a reassurance.** Cross-fitted residualisation is known on this
   project to push mean-based accuracy *below* chance (joint LDA 0.0118 on TCGA, 0.0025 on HEST). If
   the adjusted k-NN lands materially below the global-null **median**, that is evidence the residual
   has been anti-correlated with the label, which is a statement about the first moment and not about
   neighbourhoods. The null median is reported beside every observed value for exactly this reason,
   and a far-below-median result will be flagged as suspicious rather than banked as clean.
3. **A clean k=15 with a dirty k=1 is a dirty result.** Hence the sweep, and hence the max-over-k
   rule fixed above. Smoothing by a large k is the probe removing the signal, not the adjustment.
4. **Disagreement between families resolves upward.** If k-NN reads chance and the forest reads
   signal (or the reverse), the **higher** reading is the finding. A probe that finds structure is
   positive evidence; a probe that does not is not evidence of absence.
5. **Scale sensitivity is a probe artefact, not a property of the representation.** If `adjusted` and
   `adjusted_standardised` disagree, neither is quoted alone; both are reported and the disagreement
   is the finding.
6. **Prior correction must not change the sign of the conclusion.** With `OTHER` at 30% of rows, a
   plain vote can look clean purely by collapsing onto the majority class. Both the plain and the
   prior-corrected vote are reported; if only the plain one is clean, the reading is B or C, not A.
7. **Capacity check.** 2,766 patients over 85 classes is ~33 per class in 256 dimensions. If the
   measured global-null p95 for k-NN is far above 1/85, the probe is capacity-bound and the correct
   statement is "this probe is too weak here", not "the representation is clean".
8. **A favourable TCGA result does not retract the spatial finding.** The two cohorts are not
   exchangeable (HEST spots physically overlap; TCGA patients do not). If TCGA reads A, the supportable
   conclusion is that the certificate's pass generalises *on TCGA*, and the spatial 1c caveat stays
   attached to the spatial numbers.

### 5. What would make me distrust an **unfavourable** result

Stated for symmetry, so that outcome B is not accepted more cheaply than outcome A.

* **Leakage through the patient/site key.** If any patient identifier appears twice, a k-NN could name
  the site by finding the same patient. Uniqueness of `patient_ids` within the analysis partition is
  asserted before any probe runs, and the assertion is reported.
* **Cancer standing in for site.** Because site nests inside cancer, a probe that has merely recovered
  *cancer* would score far above 1/85 on site while carrying no site information beyond cancer. This
  is why the within-cancer null and the analytic `within_stratum_chance` are both computed: an
  unfavourable site reading that sits *below* the within-cancer null is a cancer finding, not a site
  finding, and will be reported as such.
* **Residual folds versus probe folds.** `cross_fitted_residuals` uses its own `KFold(seed=42)` and the
  probe uses `_stratified_folds(seed=42)`. These are different partitions, so a probe test row can be
  a residual train row. This cannot manufacture the confound (the residual is a *removal*), but it is
  recorded as a known non-independence.

---

### 6. Discipline

Thread caps `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1`, process parallelism, CPU only (the GPU is on an
ALCHEMIST download). Workspace deployed by `git -c core.autocrlf=false archive HEAD` and verified
file-by-file against `git ls-tree -r HEAD` blob SHA-1 before anything runs. Outputs to persistent NFS.
`NOTEBOOK.md`, the paper drafts and `claim_guards.py` are not edited; prose corrections are **flagged**
in the result entry, not applied.
