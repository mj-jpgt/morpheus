# KILL REPORT — Feasibility, Thesis T3

**Verdict: KILLED AS SPECIFIED.**
**Date:** 2026-07-29 · **Role:** feasibility assassin (no wet lab, open-access data only, one A100)
**Method:** direct queries to the GDC REST API (`api.gdc.cancer.gov`, plain-URL, no auth), NCBI E-utilities,
Semantic Scholar, and inspection of the MORPHEUS repo on disk. WebSearch not used (exhausted).
**Every number below is a live API return or a file on disk, reproduced verbatim. No citation is invented.**

---

## 0. What is being killed, precisely

T3 is not killed because the data does not exist, and it is **not** killed on compute. It is killed because the
two properties the thesis needs simultaneously — **(A) full confound control** and **(B) statistical power to
resolve a 0.07 effect** — are **mutually exclusive on this cohort**, and the arithmetic proving that is
one line long. That mutual exclusivity destroys the *insurance policy* ("publishable even when k = 0") which is
the sole stated justification for running T3 at all.

The thesis's own framing is what makes this lethal:

> "This makes T3 worth running only if it is designed as a MEASUREMENT with a high floor rather than a HUNT
> with a low hit rate."

Correct. And the floor is not high. The floor is **not computable**.

---

## 1. THE FATAL OBSTACLE — the confound-complete analysis cannot detect its own target effect

### 1.1 Power arithmetic (Fisher z, α = 0.05 two-sided, 80% power)

Required n to detect a correlation r:  `n = ((1.96 + 0.8416) / atanh(r))² + 3`
Minimum detectable r at a given n:     `r_min = tanh(2.8016 / sqrt(n − 3))`

| r | n required |
|---|---|
| **0.07** (MORPHEUS's established control-adjusted signal) | **1,600** |
| 0.10 | 783 |
| 0.15 | 347 |
| 0.20 | 194 |

### 1.2 Apply it to the actual analysis strata

| Analysis stratum | n | df burned by confounds | effective n | **r_min** | vs. target 0.07 |
|---|---|---|---|---|---|
| Full TCGA paired cohort, cancer-type + site + purity + composition only | 6,192 | ~260 | 5,930 | **0.036** | powered (1.9×) |
| Cases with **both** grade and stage (§2.1) | 2,077 | ~160 | 1,917 | **0.064** | marginal (1.09×) |
| Held-out-cancer share of the grade+stage cohort | ~1,300 | ~100 | ~1,200 | **0.081** | **FLOOR EXCEEDS SIGNAL** |
| CPTAC, best single cancer (CCRCC, 262 TCIA subjects) | ≤262 | — | ~250 | **0.176** | 2.5× too coarse |
| CPTAC, typical cancer (~150) | ~150 | — | ~147 | **0.227** | 3.2× too coarse |
| CPTAC, smallest usable (SAR, 88) | 88 | — | 85 | **0.295** | 4.2× too coarse |

**Read row 3.** The paper's headline estimand — a confound-certified component evaluated in the
leakage-controlled held-out cancers with a stage+grade baseline — has a **detection floor of r ≈ 0.081 against a
target effect of r ≈ 0.07**. The experiment cannot see the thing it is built to measure.

### 1.3 Why this specifically kills the k = 0 insurance policy

The null claim is:

> "we measured the spectrum under full confound control and it is flat at r ~ 0.07 — the +0.07 is not a floor
> imposed by weak methods, it IS the channel"

**A null observed in an analysis whose detection floor exceeds the target effect size is not evidence of
absence.** It is evidence of insufficient n. The claim "this is the channel, not our floor" is precisely the
claim the arithmetic forbids: in row 3, the floor *is* the floor, and it is 0.081.

You can have power (row 1) or you can have the confound battery that makes the null citable (row 3). You cannot
have both. **The high-floor measurement paper does not exist.** What exists is either an underpowered null or a
powered-but-uncertified correlation — and the latter is already published (Fu et al., *Nature Cancer* 2020,
17,355 WSIs, which attributed the associations to tumour composition).

---

## 2. Supporting kills — the confound battery is not constructible

Every covariate in the pre-registered residualisation set was checked against GDC. Three of them do not exist at
the required coverage.

### 2.1 "stage, grade" — exists for 2,077 cases across 8 cancer types, not 6,192 across 32

`GET /cases?facets=diagnoses.tumor_grade&filters=[program=TCGA]` (11,428 TCGA cases):

```
g3 2,037 · g2 1,653 · high grade 398 · g1 319 · gx 109 · g4 104 · low grade 21
not reported 8 · unknown 3 · gb 2 · _missing 6,774
```
→ interpretable grade on **4,532 / 11,428 = 40%**.

`GET /cases?facets=diagnoses.ajcc_pathologic_stage&filters=[program=TCGA]`:
→ `_missing 4,423`, i.e. stage on **~61%**.

Intersection, queried directly (`grade ∈ {G1..G4} AND ajcc_pathologic_stage ∈ {Stage I..IVC}`):

> **pagination.total = 2,077**, distributed as
> KIRC 526 · HNSC 435 · STAD 408 · LIHC 351 · PAAD 179 · ESCA 129 · CHOL 48 · CESC 1
> (sums exactly to 2,077 — the list is complete: **8 projects**)

**Consequences:**
- BRCA, LUAD, LUSC, COAD, PRAD, THCA, SKCM, BLCA, GBM, LGG, OV — none of the large cohorts — appear. The
  grade+stage-complete slice is **8 of 32 cancer types and ~34% of the paired cohort**.
- Missingness is **perfectly confounded with cancer type**: GBM/LGG/LAML have no AJCC stage *by definition*;
  BRCA has no recorded grade. Once you residualise on cancer type (which T3 does), stage and grade become
  within-stratum covariates that are **structurally absent from whole strata**. There is no imputation that is
  honest here — you would be imputing a variable that does not exist for that disease.
- The prognostic claim (v) requires a baseline of "stage + grade + subtype + known signatures" **in the
  held-out cancers**. That baseline is unbuildable for most of the held-out set.

### 2.2 "RNA-quality proxies" / RIN — missing for 86% of TCGA

The GDC analyte record *does* carry the field (confirmed: `rna_integrity_number` present on RNA analytes, e.g.
`TCGA-B0-5695-01A-11R` = 9.0). Coverage:

`GET /cases?filters=[program=TCGA AND NOT rna_integrity_number]` → **pagination.total = 9,800**

→ **9,800 / 11,428 = 85.8% of TCGA cases have no RIN on any analyte.** RIN is available for ~1,628 cases (14%).

This is the worst possible configuration, because RIN is simultaneously:
1. **the confound with the largest demonstrated morphology channel** — PathQC (*Bioengineering* 13 (2026),
   doi:10.3390/bioengineering13010005; preprint bioRxiv 10.1101/2025.09.29.679347) predicts RIN from H&E alone at
   **mean r = 0.47** (autolysis r = 0.45) on GTEx, i.e. **~7× the entire effect T3 is chasing**; and
2. **unmeasurable in 86% of the discovery cohort.**

And on the 14% where it exists, TCGA's own biospecimen QC applied an **RIN inclusion threshold** before
sequencing — the surviving RIN distribution is range-restricted (the sampled record is 9.0). Range restriction
mechanically attenuates any correlation with RIN. So validity certificate item **(ii) — "FAILS to predict …
RIN above chance" — passes for a trivial, artefactual reason on a non-random 14% subsample.** A certificate that
passes because the variable was truncated certifies nothing, and a competent reviewer will say exactly that.

### 2.3 "submitting site/batch" — 241 levels, and site is the reason the design has no headroom

`GET /cases?facets=tissue_source_site.code&filters=[program=TCGA]` → **241 distinct TSS codes** on 11,428 cases
(34 codes with <20 cases; top site AB = 200). Mean ≈ 47 cases/site.

Site *must* be controlled — Howard et al. (*Nat Commun* 12 (2021), doi:10.1038/s41467-021-24698-1, 294 citations)
showed TCGA submitting site is trivially decodable from H&E and survives stain normalisation. But on the 6,192
paired subset this is ~200 dummy levels, heavily nested within cancer type. In the small held-out cancers
(CHOL n=48 in GDC; ACC, UVM, DLBC, UCS, KICH, MESO all <100) the site df approaches a double-digit fraction of n.
This is where the ~100–260 df in the §1.2 table comes from, and it is not optional.

---

## 3. Supporting kills — the replication legs (iii) and (iv) are arithmetically closed

### 3.1 CPTAC is not a replication cohort for a 0.07 effect; and the paired n is not even queryable in one place

`GET /cases?filters=[project=CPTAC-3]` → **1,683 cases**.
`GET /cases?filters=[project=CPTAC-3 AND files.data_type="Slide Image" AND files.experimental_strategy="RNA-Seq"]`
→ **pagination.total = 0.**

CPTAC WSIs are **not in GDC** — they live on TCIA. So the paired WSI∩RNA∩proteomics cohort requires a manual
**three-portal join (TCIA imaging ↔ GDC RNA-Seq ↔ PDC proteomics)** across three different ID conventions, and
the resulting n is *unverified* (the scout doc was right to refuse to assert it; TCIA collection subject counts
of 88–262 are an **upper bound** on each leg, and the triple intersection is strictly smaller).

Even taking the upper bound as real, §1.2 rows 4–6 apply: **r_min = 0.176 to 0.295 per cancer, against a target
of 0.07. CPTAC is 2.5–4× underpowered per cancer, and no engineering fixes that.**

**The self-contradiction that closes the loop:** the only way to get CPTAC over the power line is to *pool across
cancer types* (n ≈ 1,200 → r_min = 0.081, still marginal). But pooling across cancer types reintroduces exactly
the cross-cancer cohort structure that MORPHEUS has itself already established as a **46–49% artefact** and which
T3's residualisation exists to remove. **The only configuration with adequate power is the one the thesis
definitionally forbids.** Item (iii) "loadings replicate when refit independently in CPTAC" and item (iv)
"protein shadow in CPTAC proteomics/phospho" both die here — (iv) harder, because phospho coverage is a subset of
proteomics coverage and phospho is acutely ischemia-sensitive.

### 3.2 GTEx replication is an estimand swap, not a replication

GTEx is large (25,306 samples / 970 donors), fully open, and — uniquely — its H&E and RNA come from the **same
sample**. But it is **non-cancer with no outcomes**, and PathQC (above) establishes that its leading
morphology-legible molecular axes are **RIN (r = 0.47) and autolysis (r = 0.45)**. Jones/Gundersen/Engelhardt
(bioRxiv 10.1101/2022.06.10.495669, n = 13,360) ran the joint analysis and recovered, verbatim, *"cell-type
heterogeneity, sample ischemic time, and donor health and demographics"* — including mechanical ventilation.

"Loadings of a tumour morphology↔transcriptome channel replicate when refit in GTEx" is not a replication of the
same estimand. If the loadings *do* replicate in post-mortem normal tissue, the most parsimonious reading is that
you have rediscovered a tissue-composition/preanalytic axis — which is failure mode #1 in the thesis's own list.

### 3.3 The ΔC-index claim (v) has too few events

`GET /cases?filters=[project ∈ {KIRC,HNSC,STAD,LIHC,PAAD,ESCA,CHOL} AND vital_status="Dead"]`
→ **pagination.total = 907** deaths across *all* cases of the only 7 projects that carry both grade and stage
(per-project, from the global facet: HNSC 224 · KIRC 177 · STAD 175 · LIHC 132 · PAAD 100 · ESCA+CHOL = remainder).

Restricting to the held-out-cancer subset leaves roughly **400–500 events**. A bootstrap CI on ΔC-index at that
event count has a half-width on the order of 0.03; the increment from a covariate with marginal r ≈ 0.07 over an
already-strong stage+grade+subtype baseline is realistically **<0.01**. "Bootstrap CI excluding 0" is not
attainable. This is consistent with the independent ST benchmark's own survival result — *"borderline
statistically significant for three methods"* (*Nat Commun* 16 (2025), doi:10.1038/s41467-025-56618-y; verified
live via Semantic Scholar as a real paper, "Benchmarking the translational potential of spatial gene expression
prediction from histology").

---

## 4. The uncalibratable-instrument problem (why the null is not merely underpowered but *uninterpretable*)

This is independent of n and cannot be fixed with more data.

To claim "r ≈ 0.07 is the channel, not a floor imposed by weak methods," you must exclude the
over-residualisation explanation. The only instrument for that is a **positive control**: a signal of known,
non-trivial strength that survives the *identical* residualisation and is still recovered.

**Every known-strong morphology↔transcriptome axis is inside the residualisation set by construction:**

| Known-strong axis | Status in T3 |
|---|---|
| tumour purity | residualised out |
| cell-type composition / deconvolution proportions | residualised out |
| immune infiltration | residualised out (it *is* the deconvolution proportions) |
| cancer type | residualised out |
| grade | residualised out — **and grade is literally a morphological variable** |
| stage | residualised out |

And validity certificate **(ii) actively requires the surviving component to FAIL at predicting them.**

**Therefore, after residualisation, there is by construction no signal of known ground-truth strength left with
which to demonstrate that the pipeline still works.** A flat spectrum is observationally equivalent to:
(a) there is no channel; (b) the ~260-df residualisation removed the biology along with the nuisance (composition
*is* partly the biology in a tumour); (c) attenuation from image/RNA tissue mismatch (§5).

The reviewer's one-sentence rejection writes itself: *"The authors have not demonstrated that their adjustment
preserves signal known to be present; their null is therefore uninterpretable."* There is no experiment in the
proposed design that answers it.

---

## 5. An unmeasured attenuation term that makes "we report the spectrum" unclaimable

`GET /files?facets=experimental_strategy&filters=[program=TCGA AND data_type="Slide Image"]`
→ **30,326 slide images = "Tissue Slide" 18,425 + "Diagnostic Slide" 11,901.**

These are **not equivalent for this study**:
- **Diagnostic (DX)** slides are FFPE sections from the *diagnostic block* — a **different tissue block** from the
  one the RNA aliquot was extracted from.
- **Tissue (TS/BS)** slides are the frozen top/bottom sections that *bracket* the piece used for molecular
  extraction — genuinely adjacent to the RNA.

Grepping the repo (`grep -rniE "diagnostic|\bDX\b|slide_type|frozen|ffpe"` over `*.py|*.md|*.yaml|*.json`)
returns **no slide-type field anywhere** — not in `src/data/hoptimus_patch_store.py`
(`REQUIRED_METADATA_COLUMNS`), not in `src/data/wsi_patch_bags.py`, not in the split builders. Slide provenance is
**uncontrolled and unrecorded** in the 6,192-case asset.

This injects an unknown, cancer-type- and site-correlated attenuation into the image↔RNA pairing. It is
survivable for a *prediction* paper. It is fatal for a **measurement** paper, because the whole product being
sold is a calibrated number ("the spectrum, reported for the first time"), and you cannot convert an observed
correlation into a statement about the underlying channel without the attenuation constant. That constant is
estimable only in a same-sample cohort — i.e. GTEx — which is normal tissue dominated by ischemia and RIN, so it
does not transfer.

---

## 6. What is NOT a valid kill (stated so the surviving kills are credible)

- **Compute is a non-issue.** Patch features are already extracted and frozen. Residualisation, cross-fitted
  sparse CCA, and a 1,000-permutation null on a 6,192 × d matrix are minutes-to-hours on one A100-40GB. The
  ">10,000 MSigDB/Reactome set" screen is a single matrix multiply against precomputed ssGSEA/mean-z scores.
  **Compute is a rounding error and must not be offered as the reason T3 fails.**
- **The discovery cohort exists.** 6,192 paired TCGA cases with H-Optimus features + BulkFormer embeddings are on
  disk. The pairing is real.
- **Purity and composition are obtainable.** ESTIMATE runs locally from expression; xCell/quanTIseq/EPIC run
  locally; published pan-cancer consensus purity tables are open.
- **The "no head-to-head ST-vs-bulk" caveat in the thesis is correct** and correctly stated. That part of the
  thesis is intellectually honest and I did not find grounds to dispute it.
- **The cell-line→tumour gap is irrelevant to T3.** T3 is a human-tumour observational measurement; CCLE/GDSC/
  Perturb-seq are not on its critical path. That standard attack does not apply here.

---

## 7. Secondary defects worth recording

1. **The split spec in the thesis does not match the repo.** The thesis says "leakage-controlled held-out-cancer
   split (14 dev / 21 test cancers)" — that sums to **35**, against a stated 32 cancers. The code
   (`v2/build_paired_split.py`) defaults to `--expected-development-cancers 11` and
   `--expected-heldout-cancers 22` (= 33). A pre-registered measurement paper cannot ship with its central split
   miscounted in the abstract.
2. **The |r| < 0.8 novelty threshold is far too permissive.** A component correlating **0.79** with a
   proliferation or EMT meta-signature would be reported as "not reducible to any catalogued gene set." Screening
   10,000+ sets and accepting the null at 0.8 is an *accept-the-null* procedure with no power analysis attached.
   Reviewers will demand ~0.5, at which threshold novelty almost certainly fails. This is cheap to fix and should
   be fixed regardless of T3's fate.
3. **Prior art occupies more of the design than the thesis concedes** (per the existing scout doc, whose
   citations I did not re-verify except where noted): HPL (*Nat Commun* 2024, doi:10.1038/s41467-024-48666-7) and
   its colon follow-up (*Nat Commun* 2025, independent trial N = 1,213) already execute
   discover→replicate→prognostic; Fu et al. (*Nature Cancer* 2020) already ran the pan-cancer sweep at 17,355
   WSIs and attributed it to composition. This is a novelty problem, not a feasibility one, but it removes the
   fallback of "even if the effect is small, the framing is new."

---

## 8. What survives (narrow, and it is not a paper)

One thing here is cheap, sound, and worth doing **as an internal diagnostic**:

> On the full 6,192, residualise both modalities on **only the covariates that actually exist at full coverage**
> — cancer type, submitting site, purity/ESTIMATE, deconvolution proportions — and report the cross-fitted
> canonical-correlation spectrum against a permutation null. This stratum *is* powered (r_min = 0.036, §1.2 row 1)
> and costs hours.

Its honest title is *"how much of our +0.07 is composition"* — a **repo-internal diagnostic that informs T1/T2**,
not a Nature-tier deliverable. It **cannot** carry the words "confound-certified" (no grade, no stage, no RIN),
cannot carry certificate (ii) (RIN untestable), cannot carry (iii)/(iv) (CPTAC underpowered), and cannot carry
(v) (too few events). Run it as a 1-day sanity check on the existing finding. Do not build a thesis on it.

---

## 9. Citation hygiene

**Verified live this session:** all GDC facet/count queries above (`api.gdc.cancer.gov/cases`, `/files`, public,
no auth — re-runnable verbatim from the URLs in §2–§5); doi:10.1038/s41467-025-56618-y confirmed real via
Semantic Scholar Graph API as *"Benchmarking the translational potential of spatial gene expression prediction
from histology"* (bioRxiv 2023 → Nat Commun 2025, 51 citations on the preprint record).

**Carried from `thesis_T3_denovo_discovery.md` without independent re-verification this session** (that document
logs its own sources; flagged here so nothing is laundered): PathQC doi:10.3390/bioengineering13010005; Howard
et al. doi:10.1038/s41467-021-24698-1; Fu et al. doi:10.1038/s43018-020-0085-8; Quiros et al.
doi:10.1038/s41467-024-48666-7; Jones/Gundersen/Engelhardt bioRxiv 10.1101/2022.06.10.495669; CPTAC TCIA
per-collection subject counts (88–262).

**COULD-NOT-VERIFY:** the CPTAC WSI∩RNA∩proteomics paired sample count — GDC returns **0** slide images for
CPTAC-3, confirming imaging is TCIA-only; the true triple-intersection n is not obtainable from any single API
and remains unasserted. Semantic Scholar returned HTTP 429 on one of two calls (the HEST leaderboard query); the
HEST numbers quoted in the thesis statement were **not** re-verified here and nothing in this kill report depends
on them.
