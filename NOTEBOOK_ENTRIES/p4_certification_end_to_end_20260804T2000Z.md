## 2026-08-04 20:00 UTC — The five-point certification rule, run end to end for the first time: 3 of 5 conditions pass on the best axis, the exposable state is answerable from ZERO axes, and the certificate refuses a planted site code only above a measured strength

**Logged:** 2026-08-04 20:00 UTC. **Predeclared in**
`NOTEBOOK_ENTRIES/PREDECLARED_p4_certification_system_tests_20260804T1750Z.md`, committed (`dfefc9b`)
before the harness that produces any number below was written, with one amendment appended and
committed before the full runs (`§AMENDMENT`, on evidence from a reduced-setting smoke run).

**How obtained.** Workspace `~/ws_p4/morpheus` on the A100 (`150.136.45.194`), built from
`git -c core.autocrlf=false archive HEAD` and verified **618/618 files by git blob SHA-1** (0
mismatched, 0 missing, 0 extra). CPU only, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`; the GPU was at 97% utilisation for other agents throughout
and was not touched. Driver `v2/research/rebase/nature/p4_certification/p4_certify.py`, which
**defines no statistic** — every number comes from `calibra.confound_certificate.certify_axes`,
`calibra.calibration.spike_recovery_curve`, `calibra.spectral.heldout_single_direction_correlation`
and `paired_absolute_correlation`, `calibra.leave_sites_out.site_folds` and `evaluate_fold`,
`calibra.residualise.*`, and `calibra.confound_certificate.within_stratum_permutations`.
`claim_guards.py` and `claim_evidence.json` are **not edited**.

Artifacts: `runs/d2_final/artifacts/d2_{h,i}_seed42.npz`, `test` partition, n = 2,766, 85 pooled TSS
classes at `min_site_count = 10`, chance rate 1/85 = 0.011765, 1,000 within-cancer label
permutations (resolution 1/1001), 200-draw axis bootstraps.

---

### Bad news first

**On the state a user would actually be shown, the fraction of axes a query layer is allowed to
answer from is 0.000, in 6 of 6 artifact × state cells.** Not small — zero. Under the predeclared
reading, **P4 is not a system today, it is a refusal.**

**The one condition of five that fails is condition 3, and it fails on an axis that is individually
innocent.** The most legible axis in the whole representation predicts tissue source site *below* its
own permutation null. It is refused because the representation it belongs to leaks site jointly, at
31× chance, in a direction no per-axis screen looks along.

**Condition 4b cannot be evaluated at all**, because no external cohort with paired morphology is on
disk. That is absence of data, not a measured failure, and the two are kept apart throughout.

---

## Test A — the abstention rule

**Definition, fixed in the predeclaration.** An axis is *answerable* iff (A1) its out-of-fold
balanced accuracy for pooled TSS does not exceed the 95th percentile of the within-cancer permutation
null, **and** (A2) it is not named in `breaching_axes`, **and** (A3) the state's joint LDA test passes.
A3 is a property of the state and gates every axis in it at once, because a state that passes per-axis
and fails jointly has not been shown to be site-free — it has been shown that nobody looked with
enough resolution.

**Result.** All twelve cells (2 artifacts × 3 states × raw/adjusted) were **re-derived from scratch in
this workspace** at 1,000 permutations, and every breaching count and every joint accuracy
reproduces the published run of 2026-08-03 (`t13_confound_certificate_20260803T0152Z.md`,
`t13_adjusted_certificate_and_p6_20260803T0300Z.md`) to the digits those entries quote. The
answerable fractions are computed from that run under the definition above.

| artifact | state | arm | breaching axes / 256 | joint LDA | joint null p95 | joint OK? | **answerable, per-axis only** | **answerable, strict** |
|---|---|---|---:|---:|---:|:---:|---:|---:|
| d2_h | wsi_biology | **raw** | 17 | 0.3633 | 0.1539 | ✗ | 0.934 | **0.000** |
| d2_h | full_biology | raw | 60 | 0.2630 | 0.1742 | ✗ | 0.766 | **0.000** |
| d2_h | rna_biology | raw | 58 | 0.2563 | 0.1744 | ✗ | 0.773 | **0.000** |
| d2_i | wsi_biology | **raw** | 43 | 0.2348 | 0.1224 | ✗ | 0.832 | **0.000** |
| d2_i | full_biology | raw | 61 | 0.2689 | 0.1856 | ✗ | 0.762 | **0.000** |
| d2_i | rna_biology | raw | 48 | 0.2744 | 0.1861 | ✗ | 0.813 | **0.000** |
| d2_h | wsi_biology | **adjusted** | 0 | 0.0118 | 0.0528 | ✓ | 1.000 | **1.000** |
| d2_h | full_biology | adjusted | 0 | 0.0101 | 0.0668 | ✓ | 1.000 | 1.000 |
| d2_h | rna_biology | adjusted | 0 | 0.0074 | 0.0654 | ✓ | 1.000 | 1.000 |
| d2_i | wsi_biology | **adjusted** | 0 | 0.0052 | 0.0418 | ✓ | 1.000 | 1.000 |
| d2_i | full_biology | adjusted | 0 | 0.0085 | 0.0758 | ✓ | 1.000 | 1.000 |
| d2_i | rna_biology | adjusted | 0 | 0.0079 | 0.0732 | ✓ | 1.000 | 1.000 |

**Reading, against the predeclared rule.** The first predeclared branch is taken and it is met in its
strict form: the strict answerable fraction on the raw state is **0 in 6 of 6 cells**, not 5 of 6.
The second branch is also taken — the adjusted fraction is 1.000, which is above the 0.90 threshold
at which the predeclaration forbids reading it as a pass until the discrimination test returns.

**The gap between the two readings is the finding.** A per-axis-only certificate would answer from
76–93% of axes on the raw state. The joint row takes that to zero. The difference between
"a certificate" and "a certificate with a joint row" is, on these artifacts, the difference between a
working interface and no interface at all.

---

## Test A′ — does the certificate discriminate, or does it pass everything?

Three axes of known character appended to `d2_h_seed42::wsi_biology` before certification: a site
code, a cancer code, and pure Gaussian noise, over a signal-to-noise ladder. SNR 1.0 is the
predeclared arm; the ladder is a declared addition, recorded in the harness docstring and in the
commit that introduced it, because a single-axis nearest-class-mean rule at 85 classes is a weak
detector by construction and a binary verdict at one SNR would not say *how strong* a site code has
to be before it is refused.

**Predeclared pass pattern: the site code must breach, the cancer code and the noise must not.**
Chance rate 0.011765; per-axis nulls are within-cancer permutations, 1,000 draws; 259 axes scored
(256 real + 3 planted).

**Raw state — the certificate discriminates, at every strength tested.**

| SNR | PLANT_site | its null p95 | perm p | breaches? | PLANT_cancer | its null p95 | breaches? | PLANT_noise | breaches? | state joint LDA |
|---|---:|---:|---:|:---:|---:|---:|:---:|---:|:---:|---:|
| 0.5 | 0.0245 | 0.0205 | 0.0040 | **YES** | 0.0159 | 0.0253 | no | 0.0083 | no | 0.3795 |
| **1.0 (predeclared)** | 0.0287 | 0.0229 | 0.0010 | **YES** | 0.0280 | 0.0362 | no | 0.0083 | no | 0.4292 |
| 2.0 | 0.0705 | 0.0264 | 0.0010 | **YES** | 0.0441 | 0.0554 | no | 0.0083 | no | 0.5264 |
| ∞ (noiseless) | **1.0000** | 0.0397 | 0.0010 | **YES** | 0.2372 | 0.2527 | no | 0.0083 | no | 0.7281 |

**The predeclared pattern holds at 4 of 4 strengths, including one weaker than the predeclared arm.**
Three things make this a real control rather than a formality:

* A **noiseless site code reads balanced accuracy 1.0000** and is refused at the resolution floor. The
  instrument is not blunt.
* A **noiseless cancer code reads 0.2372 and is NOT refused**, because its own within-cancer null sits
  at 0.2527. That is the null doing exactly the job its docstring claims: cancer type predicts site
  well above chance, and the certificate correctly declines to charge an axis for lineage information
  it did not add. Had the null been a global permutation, this axis would have been refused and every
  biology axis with it.
* The **state's joint LDA tracks the injected leak monotonically** — 0.3795 → 0.4292 → 0.5264 →
  0.7281 as the planted code strengthens, from a baseline of 0.3633 with nothing planted. The joint
  row is sensitive to a single added axis out of 257.

**Adjusted state — the certificate is near-inert against a site code at every realistic strength.**

| SNR | PLANT_site | its null p95 | perm p | breaches? | PLANT_cancer breaches? | state joint LDA | joint OK? | axes breaching |
|---|---:|---:|---:|:---:|:---:|---:|:---:|---:|
| 0.5 | 0.0028 | 0.0183 | 0.9990 | no | no | 0.0067 | ✓ | 0 |
| **1.0 (predeclared)** | 0.0090 | 0.0183 | 0.8601 | no | no | 0.0067 | ✓ | 0 |
| 2.0 | 0.0100 | 0.0180 | 0.7872 | no | no | 0.0081 | ✓ | 0 |
| ∞ (noiseless) | **0.2751** | 0.0236 | 0.0010 | **YES** | **YES** | 0.0387 | ✓ | 2 |

**Reading, against the predeclared rule, and the first branch is taken.** `PLANT_site` **stops
breaching after adjustment** at SNR 0.5, 1.0 and 2.0 — its accuracy falls to 0.0028–0.0100, at or
below the chance rate of 0.0118. Per the predeclaration:

> the adjusted certificate cannot tell a site code from a biology axis at any realistic strength,
> because the adjustment removes site by construction. **"The adjusted state certifies" is a
> statement about the adjustment, not evidence that any axis is site-free**, and the adjusted arm may
> not be quoted as certification for exposure.

The boundary is measured rather than asserted: the adjusted certificate retains power only against a
**perfectly noiseless** site code, somewhere above SNR 2.0. It is not inert; it is close to it over
the range any real axis occupies.

**And it errs in the other direction too, which is new.** At SNR = ∞ the adjusted arm refuses
`PLANT_cancer` as well (0.2459, p at the resolution floor). A pure cancer code carries no site
information beyond cancer and must not breach — the raw arm correctly declines to refuse it. It
breaches after adjustment because the design **annihilates** it, leaving numerical residue that the
certificate's per-axis standardisation then rescales to unit variance before scoring. **So on
degenerate or near-null axes the adjusted certificate produces false refusals as well as false
passes.** That is a caveat that belongs with the instrument and is recorded here rather than in the
place where it would be convenient to omit it.

---

## Test B — the five conditions, end to end, on one real axis

**Axis, selected before any condition was scored**, as the argmax of
`|heldout_single_direction_correlation|` over 256 axes × 90 non-control targets on the adjusted block
of `d2_h_seed42::wsi_biology` — 23,040 pairs:

> **axis 46**, against **HALLMARK_ALLOGRAFT_REJECTION** (`hallmark_in_training`), out-of-fold
> single-direction correlation **0.47035**.

This is the most favourable case available, so a failure here decides the gate in the negative for
every weaker axis.

| # | condition (`MULTIMODAL_EXPANSION.md` §1, verbatim) | verdict | the number |
|---|---|:---:|---|
| 1 | the operator is estimated **on a discovery fold only** | **PASS** | 3,118 train / 543 val / 2,766 test on the artifact's own `split`; **0** patients in both; the readout direction is fit out of fold by construction |
| 2 | its axes clear the **CALIBRA detection floor** | **PASS** | `detection_floor` **0.05**, `transmission_floor` 0.01, level-0 confound-induced baseline 0.022; the axis reads **0.4703**, i.e. **9.4×** the floor |
| 3 | they pass the **confound certificate** | **FAIL on the exposable state** | the axis is clean per-axis on both arms (raw balanced accuracy **0.0272** vs its null p95 **0.0358**; adjusted 0.0071 vs 0.0200) — but the **raw state's joint LDA is 0.3633** against a null p95 of 0.1539 and a chance rate of 0.0118, permutation p at the 1/1001 floor. Adjusted state: joint 0.0118, ≤ chance, **PASS** |
| 4a | they **replicate** in untouched patients | **PASS** | whole sites held out, 5 folds, 453–864 patients and 32–89 sites each. Out-of-site 0.4334 / 0.4753 / 0.5253 / 0.4959 / 0.5167, median **0.4959** against a within-cohort **0.4765** → **retained fraction 1.041** against a predeclared bar of 0.50. **5/5 folds clear their own permutation null (p95 ≈ 0.09–0.10) and their site-cluster CI clears it too.** Sign agrees 5/5 |
| 4b | … and in **≥1 external cohort** | **UNEVALUABLE** | the only external material on the box is `external/cptac_gdc_rna` — 2,724 STAR-counts files, **RNA only, no slides**. Verified: `api.gdc.cancer.gov` returns **0** CPTAC slide images. HEST-1k exists but is a different modality against a different target space and carries no D2 axis |
| 5 | failures are **recorded and exposed** alongside successes | **PASS (procedural)** | `certify_axes` names `breaching_axes` rather than dropping them; `ClaimVerdict.as_rows` emits **3** visible NaN status rows for a `legible_axis` claim (`composition_attribution`, `purity_confound`, `no_external_cohort`) |

**Gate: not met.** Conditions 1, 2 and 5 pass; condition 3 fails on the state that would be exposed;
condition 4 passes in its untouched-patients half and is unevaluable in its external-cohort half. An
UNEVALUABLE condition counts as *not passed* — you cannot certify on evidence that does not exist —
and is kept distinct from a measured failure throughout.

### The three things this run establishes that were not known before

1. **A per-axis-only certificate would have certified this axis.** Its own site accuracy is *below*
   its own null on both arms. It is refused only by the joint row. The joint test is therefore a
   required field of P4's certificate schema, demonstrated on the axis a query layer would most want
   to expose, not argued from principle.
2. **Condition 2 is comfortably met and had never been checked.** The channel on the single best axis
   is 9.4× the smallest spike the same analysis can resolve, with a confound-induced level-0 baseline
   of 0.022 — small, so the floor is not being manufactured by the adjustment.
3. **Condition 4a is met at full strength on the adjusted state and fails on the raw one.** Holding
   out whole tissue source sites costs this axis nothing (retained fraction 1.041). A concurrent
   agent's block-level run on the same artifact and targets
   (`leave_sites_out_result_20260804T1830Z.md`) reports the same for the multivariate channel —
   adjusted survives 5/5 at a ratio of **1.010** to a matched random split — and reports that the
   **unadjusted arm collapses at 2 of 5 folds**. So the raw state fails conditions 3 *and* 4a, and the
   adjusted state passes both.

### The structural obstacle this makes unavoidable

Everything above reduces to one fact that had not been named: **the state that certifies cannot be
exposed.** `cross_fitted_residuals` fits its nuisance model on folds of the very rows it is scoring,
against a design built from the cohort's own cancer and tissue-source-site composition. There is no
operator on this project that takes a *new* slide and returns its adjusted coordinates. "Expose the
adjusted axis" is therefore not an available option at query time — it is an option only for patients
already inside the cohort the adjustment was fitted on.

That makes condition 3 a **build item, not a measurement item**: an *inductive* adjustment operator —
nuisance model fitted once on the discovery fold, applied unchanged to a held-out or external patient
— must exist before any axis can be both certified and exposed. It does not exist, it is not large,
and nothing else on the project has needed it, which is why nobody has written it.

---

## Test C — the competitor gap, and the first number P4's contribution has ever had

CellWhisperer (Schaefer et al., *Nature Biotechnology*, online 2025-11-11, DOI
10.1038/s41587-025-02857-9) ships a chat box over CELLxGENE whose abstract contains no uncertainty,
abstention, calibration or refusal language. P4's contribution is exactly the refusals such a system
does not make, and until now that contribution had no number.

**Query set.** The 90 non-control targets of `frozen_rna_targets.npz` (50 `hallmark_in_training`,
24 `heldout_pathway`, 8 `immune_tme`, 8 `tumour_state`). Query *t* = "report target *t* for this
patient from the slide." State `d2_h_seed42::wsi_biology`, test partition, adjusted block, one
supporting axis per target — the axis with the largest `|heldout_single_direction_correlation|`.

**Policy N (competitor-style).** Answer every query. `n = 90` by construction; that is what "no
abstention" means.
**Policy C (ours).** Answer only if the supporting axis passes the confound certificate on the
exposed state *and* the correlation exceeds that (axis, target) pair's own CALIBRA detection floor
*and* it exceeds the 95th percentile of a 200-draw within-cancer row-permutation null.

| | count |
|---|---:|
| queries | 90 |
| **answered by Policy N (competitor-style)** | **90** |
| **answered by Policy C (certified)** | **0** |
| **gap — confidently-worded answers our rule refuses** | **90** |
| refused: supporting axis / state fails the site certificate | 90 |
| refused: below the CALIBRA detection floor | **62** |
| refused: inside the within-cancer permutation null | 1 |

**The predeclared `n_C = 0` branch is taken, so the honest statement is the one the predeclaration
fixed in advance: our system answers nothing, and P4's contribution today is entirely the refusal.**

**The number that matters more than the gap, because it survives the one buildable fix.** The site
condition refuses all 90 through a single state-level failure (the joint LDA row), so the gap is
dominated by one fact. Removing that condition — i.e. asking what a certified interface would answer
*after* the inductive adjustment operator of §"structural obstacle" exists — leaves:

> **28 of 90 queries answerable**, and **19 of those 28 are `hallmark_in_training`** — targets the
> representation was supervised on. Of the 24 genuinely untrained `heldout_pathway` targets, **1**
> survives. The rest: 5 of 8 `immune_tme`, 3 of 8 `tumour_state`.

So even in the most favourable counterfactual, a certified interface over this representation answers
under a third of the queries a CellWhisperer-style one answers, and two thirds of what it does answer
are pathways it was trained to reproduce.

**Why 62 of 90 fall below the floor, and it is not a small effect.** The per-(axis, target) detection
floors, from `spike_recovery_curve` on the level grid (0, .01, .02, .05, .10, .20, .40):

| detection floor | 0.02 | 0.05 | 0.10 | 0.20 | 0.40 | **unresolvable (NaN)** |
|---|---:|---:|---:|---:|---:|---:|
| targets | 1 | 8 | 11 | 20 | 27 | **23** |

against observed `|correlation|` running 0.0568 to 0.4703, median 0.2392 (p05 0.0909, p25 0.1761,
p75 0.3137, p95 0.4338). For 23 of 90 targets the floor is **NaN** — not even a planted spike of
r = 0.40 is recovered above the level-0 upper tail in 80% of draws — which is scored FAIL, because an
unresolvable floor is the absence of a floor, not a pass. Those 23 are queries a competitor-style
interface answers in fluent prose about an analysis that demonstrably could not have detected the
effect it is describing.

**Supporting axes.** 90 queries are supported by only **46 distinct axes**; 10 queries are supported
by an axis that itself breaches the certificate (3 distinct axes). The strongest queries cluster
tightly: `HALLMARK_ALLOGRAFT_REJECTION` 0.470 and `HALLMARK_INTERFERON_GAMMA_RESPONSE` 0.447 both on
axis 46; `immune_t_cell_inflammation` 0.470, `immune_cytolytic_activity` 0.461 and `immune_ifng`
0.446 all on axis 79. The weakest are KEGG_MEDICUS environmental-factor signalling pathways at
0.057–0.086 against floors of 0.10–0.40. **This is P4's figure**: the same panel, sorted, with the
floor drawn per target and the refusal reason coloured.

---

## Test D — the spatial modality

**Not run here, deliberately.** A concurrent agent predeclared the spatial replication in
`PREDECLARED_spatial_claim_replication_20260804T1800Z.md` (commit `80076b4`), whose Claim 1 and
Claim 4 are exactly `certify_axes` on the HEST-1k spot artifact with **slide** in place of tissue
source site — 44 slides, 13 in the test partition, chance 1/13 = 0.0769, verified on the built
artifact by checking that `residualise.pooled_tissue_source_site` returns `slide_ids` for all 144,162
rows. Duplicating it would produce a second number under the same name.

What this work contributes instead:

* **D2 — which of the five conditions the spatial modality can even be asked.** From the artifact
  schema rather than from hope: conditions **1, 3 and 5 are askable** (slide-grouped split, slide
  labels present, the same reporting machinery); condition **2 is askable but expensive** —
  `run_calibra` on 53,217 spots was killed after 20 minutes at ~1,300% CPU; condition **4b is not
  askable**, because reading a *D2 axis* on HEST requires transporting the TCGA-trained
  representation onto HEST spots, and TCGA vs HEST embeddings separate at **AUC 0.99999** against a
  within-TCGA control of 0.5012. The transport is not free and nothing on the project has done it.
* **D3 — HEST does not discharge `no_external_cohort` as a matter of fact.** The guard's discharge
  condition is mechanically trivial (`len(external_cohorts) >= 1`) and its own remedy text names
  HEST-1k. But `claim_evidence.json` has no `external_cohorts` key, no spatial replication result
  exists, and both spatial entries state that no existing scientific claim has been re-run on spatial
  data. **No edit to `claim_evidence.json` is made by this work.**
* **D1 — their result landed while this entry was being written**
  (`spatial_claim_replication_result_20260804T1930Z.md`), and it answers the transfer question in both
  directions at once.

  **The pattern transfers.** n = 2,769 spots, 13 slide classes, chance 1/13 = 0.0769: raw **1,212 of
  1,536** axes breaching with joint LDA **0.9707** against a null p95 of 0.5465 — NOT CERTIFIED;
  adjusted **0** breaching with joint LDA **0.0025** against a null p95 of 0.0470 — CERTIFIED, a
  factor of 384 down and far below chance. Their predeclared Claim 1a bar (≤ 1.5 × chance = 0.115)
  clears with room. So the raw-FAIL / adjusted-PASS behaviour measured here at 85 site classes is not
  a TCGA artefact: it reproduces in a second modality at a seventh of the design rank. **That is
  evidence for condition 3 in a second modality, and it is not evidence for condition 4b.**

  **And it makes the adjusted PASS mean much less than it appears to.** Their Claim 1c — a probe added
  before the run because the certificate's classifier family is mean-based and the adjustment is a
  mean — runs an out-of-fold **15-NN** vote for slide on the *adjusted* residual. It reads **0.7291
  balanced accuracy, 9.5× chance**, against a **global**-permutation null of 0.0856 for the same
  classifier on the same features. Cross-fitted ridge on a one-hot design removes the confound's
  **mean vector** and nothing else.

  **This lands directly on Test A′ above and sharpens it.** My planted-axis ladder shows the adjusted
  certificate is near-inert against a *linear* site code at any realistic strength. Their kNN probe
  shows that what survives the adjustment is not merely undetectable by the certificate — it is
  recoverable by an ordinary non-linear classifier. Two independent routes to the same conclusion:
  **the adjusted state's blanket PASS is a statement about the adjustment and the classifier family,
  not evidence that any axis is confound-free.**

  **The argument is about the classifier family, not about spatial data, and it applies verbatim to
  the 85-class TCGA site certificate — which has never been asked this question.** Running the same
  kNN probe on TCGA is hours of CPU and is now the highest-value open measurement on P4. Until it is
  run, the sentence "the site signal is gone" in
  `t13_adjusted_certificate_and_p6_20260803T0300Z.md` and `paper/P1_CALIBRA_DRAFT.md` §4.2 is not
  supported; the supportable version is *the confound is removed from the first moment, and a
  mean-based certificate therefore certifies*.

* **A defect in the null convention that touches every number in Test A above.**
  `within_stratum_permutations` permutes the confound label inside cancer type, so wherever the
  confound is **nested** in the stratum the permutation is the identity and the null is handed the
  true labels. Five of thirteen HEST test slides are the only slide of their oncotree label, which is
  why their raw joint null median reads 0.533 against a chance rate of 0.0769. On TCGA, tissue source
  site is partially nested in cancer, so a milder version of the same inflation is present in the
  null p95 values quoted in Test A. **It makes the certificate more permissive, so it cannot rescue a
  failure** — the raw joint accuracies of 0.235–0.363 exceed even an inflated null p95 of 0.12–0.19 —
  but any *pass* read against a within-cancer null should be quoted with the global-permutation null
  beside it, and the adjusted arm's passes here are not.

---

## Suite

Run on this workspace at the commit the tests were launched from:
`pytest morpheus/v2/tests morpheus/tests --ignore=morpheus/v2/tests/test_p2_figures.py -q` →
**448 passed, 0 failed**; re-verified at the final commit on a *fresh* workspace built the same way
(**643/643 files by git blob SHA-1**, 0 mismatched, 0 missing, 0 extra) → **476 passed, 0 failed in
56.40 s**, the growth being other agents' tests landing between the two runs. `test_p2_figures.py` (24 tests) still cannot run in `~/venv`
(`ModuleNotFoundError: No module named 'matplotlib'`); nothing was installed into that environment.

## One library defect found and fixed, with a test

`residualise.cross_fitted_residuals` **crashed on a single-column matrix**. sklearn ravels a
one-column Ridge target, so `model.predict` returns `(n,)` against an `(n, 1)` slice and the
subtraction broadcast to `(n, n)` and raised. Per-axis residualisation is one column at a time and it
is the unit the certification rule works in, so this was not a corner case — it made the entire
per-axis path unreachable. Fixed by reshaping the prediction to the slice's shape, which is a strict
no-op for every `k >= 2` call site; `test_cross_fitted_residuals_handle_a_single_column` pins that by
asserting column 0 of a two-column call equals the one-column call to 1e-12. **No existing number on
the project moves.**

## One defect in a shipped flag, found and reported rather than fixed

`SpikeRecoveryResult.summary()["observed_above_floor"]` compares a **signed** correlation along a
randomly drawn direction pair against the detection floor. With a single-column `x` and `y` that sign
is the sign of a scalar draw and carries nothing: on the smoke run the flag read **False** for an
association of magnitude 0.4765 against a floor of 0.05. For a multi-column `x` the same
median-over-random-directions is centred near zero, so the flag is close to always-`False` there too.
The criterion was amended before the full runs, in a committed amendment to the predeclaration, to
the pairing the library's own docstring prescribes — magnitude of the single-direction statistic
against the floor. **The library is left unchanged**; the amendment applies at reporting time, which
is where the module docstring says magnitudes belong.

## Files / provenance

Harness `v2/research/rebase/nature/p4_certification/p4_certify.py`.
Outputs `~/ws_p4/out/{A_abstention_seed42.json, A_plant_ladder.json,
B_endtoend_d2h_seed42_wsi.json, C_competitor_d2h_seed42_wsi.json}`, vendored into
`v2/research/rebase/nature/p4_certification/`.
