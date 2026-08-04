# 2026-08-04 19:30 UTC — E1: clinical endpoints from morphology on held-out cancers. **Nothing clears its pre-registered band.** TP53 is real but under-recovers; grade is the only strong readout and it has no band; four actionable drivers and stage are at chance after adjustment

**Pre-registration:** `NOTEBOOK_ENTRIES/PREDECLARED_E1_clinical_endpoints_20260804T1745Z.md` +
`configs/PREREG_E1_clinical_endpoints.json`, committed as **f505fdc** at 17:45 UTC. First endpoint
statistic computed at 17:52 UTC. The commit precedes the run and is on `origin/research/rebase-vision`.

**How obtained:** `v2/calibra/known_covariate_control.py`, unmodified, blob SHA-1
`7ca73815cce98221b95e91c659d13ece0bb66b56` verified identical on the box after
`git -c core.autocrlf=false archive HEAD`. No statistic computed inline. Workspace
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/e1_endpoints/`, results under
`results/`, labels `inputs/e1_endpoint_labels.parquet`.

---

## Headline, bad news first

**Zero endpoints cleared the pre-registered band in the sense that matters.** Of six endpoints with
literature bands, **four (BRAF, KRAS, APC, PIK3CA) and dMMR fail outright** — their intervals cover
*measured* chance after cancer+site adjustment. **Stage fails too.** TP53 is the single endpoint that
is reproducibly above measured chance and survives adjustment, but its point estimate sits **below**
the published band in three of four arms; it clears the band only because the upper tail of its
bootstrap interval reaches 0.60. That is a **replication, not an application contribution**, and a
weak one.

The one strong readout — **histologic grade, 0.678 [0.642, 0.715] adjusted** — is exploratory, has
**no pre-registered band** because I could not retrieve one I was willing to cite, and is
morphology→morphology (a pathologist assigns grade by looking at the slide), which I declared
near-circular *before* the run.

**Per the pre-registered success criterion, this is the "zero clears" outcome**, and the finding is
about how little morphology carries once adjustment is done properly.

---

## 1. Three corrections to the record, found by re-checking rather than inheriting

The handoff warned that a claimed absence had previously turned out to be a lookup failure. It had
happened again.

**(a) `t17b_known_covariate_20260803T0155Z.md` says "there is no mutation table ... on disk". This is
false.** The MC3 public pan-cancer MAF is at
`/lambda/nfs/geeg/biorag3_persistent_20260711/raw_tcga_v1/mc3.v0.2.8.PUBLIC (1).maf.gz`, 753 MB.
Parsed: **3,600,963 variant rows, 10,224 sequenced patients**. The filename contains a space and a
`(1)`, which is the most plausible reason a glob missed it. This single file supplies **TP53 status
for 2,686 of the 2,766 held-out patients across 14 scoreable cancers** — by a wide margin the
best-covered molecular endpoint available to this project, and it was recorded as unavailable.

Label extraction validated against published TCGA frequencies — the per-cancer TP53 rates are
textbook-correct, which is strong evidence the parse is right:

| UVM | THCA | KIRC | KIRP | PCPG | CESC | LIHC | BLCA | STAD | COAD | HNSC | LUSC | ESCA | UCS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.000 | 0.003 | 0.023 | 0.024 | 0.023 | 0.049 | 0.308 | 0.493 | 0.486 | 0.611 | 0.749 | 0.847 | 0.918 | 0.918 |

**(b) MSI — partly false.** The claim that `microsatellite_instability` carries "no MSI-H/MSI-L/MSS
calls whatsoever" is **confirmed** (test partition: 7 YES / 74 NO, an assay-performed flag). But the
adjacent column **`loss_expression_of_mismatch_repair_proteins_by_ihc` carries genuine dMMR IHC calls
for 147 held-out patients, 42 positive**, and `..._by_ihc_result` gives per-protein MLH1/MSH2/PMS2/MSH6
detail for 44. dMMR by IHC is the companion diagnostic for MSI-H. A usable, if underpowered, label
existed and was missed.

**(c) No consensus-subtype table — confirmed.** No PAM50, no CMS, nothing under `/lambda/nfs/geeg`.

---

## 2. The structural constraint the handoff asked me to check — it is worse than stated, and it is decisive

The handoff noted ER/PR covers 690 patients, all BRCA, 0 in test. The underlying cause is more
general: **the D2 split is leave-cancer-types-out.** Every cancer is wholly train/val or wholly test.

* **train/val (11):** BRCA, LGG, LUSC, LUAD, KIRC, THCA, HNSC, PRAD, UCEC, GBM, OV
* **test (21):** SKCM 306, BLCA 283, LIHC 268, STAD 244, KIRP 211, SARC 210, CESC 198, COAD 197,
  PAAD 123, TGCT 114, ESCA 110, THYM 82, MESO 60, READ 59, UVM 54, KICH 50, UCS 49, ACC 48, PCPG 43,
  CHOL 29, DLBC 28

So **ER, PR, PAM50 (BRCA) and IDH (LGG/GBM) have zero held-out patients and cannot be measured at
all** — not by me, not by anyone, without a re-split. The 0.8781 ER precedent is not replicable as a
held-out result on this split. This is a property of the split, not a negative finding.

**Endpoints dropped for insufficient held-out coverage, before any modelling:**

| endpoint | why dropped |
|---|---|
| ER, PR, PAM50 | BRCA is train-only — 0 held-out patients |
| IDH1/IDH2 status | LGG and GBM both train-only — 0 held-out patients |
| HPV | **0 labelled patients on all 2,766 held-out**, in all three HPV columns, despite CESC (n=198) being a test cancer |
| HRD | no HRD/LOH score table on disk; deriving one from the raw SNP6 segments is a new method with no pre-registered basis |
| COAD CMS, molecular subtype | no subtype table on disk |

**The compensating advantage:** because the test partition spans 21 lineages, cancer-type adjustment
is a real test here, in a way the ER control explicitly could not exercise (ER was one cancer, so
pooled and within-cancer were identical by construction).

---

## 3. Results — image-only (`wsi_biology`), held-out test partition, d2_h_seed42

1,000-draw bootstrap CI; chance **measured** by 1,000 within-cancer label permutations (`null_p95`).
"Excludes chance" = CI lower bound > `null_p95`.

### 3a. Endpoints with a pre-registered band

| endpoint | n / n+ | adj | within-cancer AUROC | CI95 | **pooled** | measured chance | band | verdict |
|---|---|---|---:|---|---:|---:|---|---|
| **TP53** | 2686 / 893 | raw | 0.5791 | [0.5437, 0.6134] | **0.6883** | 0.5306 | [0.60,0.80] | above chance |
| **TP53** | 2686 / 893 | **adj** | **0.5912** | **[0.5618, 0.6274]** | 0.5777 | 0.5310 | [0.60,0.80] | **PASS band (marginal — point below band)** |
| BRAF | 2686 / 211 | raw | 0.5608 | [0.4611, 0.6244] | **0.7324** | 0.5675 | [0.64,0.82] | **FAIL — covers chance** |
| BRAF | 2686 / 211 | adj | 0.4748 | [0.3976, 0.5736] | 0.4808 | 0.5660 | [0.64,0.82] | **FAIL** |
| KRAS | 2686 / 262 | raw | 0.4824 | [0.4279, 0.5539] | **0.7529** | 0.5518 | [0.60,0.78] | **FAIL — covers chance** |
| KRAS | 2686 / 262 | adj | 0.5049 | [0.4538, 0.5706] | 0.5289 | 0.5505 | [0.60,0.78] | **FAIL** |
| APC | 2686 / 306 | raw | 0.4684 | [0.3965, 0.5313] | **0.7285** | 0.5486 | [0.60,0.78] | **FAIL — covers chance** |
| APC | 2686 / 306 | adj | 0.4564 | [0.4004, 0.5237] | 0.4584 | 0.5521 | [0.60,0.78] | **FAIL** |
| PIK3CA | 2686 / 276 | raw | 0.5369 | [0.4675, 0.5894] | 0.6871 | 0.5469 | [0.60,0.78] | **FAIL — covers chance** |
| PIK3CA | 2686 / 276 | adj | 0.5551 | [0.5017, 0.6127] | 0.5693 | 0.5479 | [0.60,0.78] | **FAIL — covers chance** |
| dMMR | 147 / 42 | raw | 0.6587 | [0.5116, 0.7527] | 0.6315 | 0.5961 | [0.66,0.80] | **FAIL — covers chance** |
| dMMR | 147 / 42 | adj | 0.5050 | [0.3557, 0.6148] | 0.4828 | 0.5904 | [0.66,0.80] | **FAIL** |

**Note the pooled column.** For BRAF, KRAS and APC the *pooled* AUROC is 0.73–0.75 while the
within-cancer AUROC is 0.47–0.56 — at or below chance. Those pooled numbers are **lineage
identification**, nothing else: BRAF is 154/303 in SKCM, APC and KRAS concentrate in colorectum. A
paper reporting "BRAF predicted from H&E at 0.73 pan-cancer" off this representation would be
reporting melanoma detection. **This is precisely the artifact the handoff asked to catch, and it is
the single most reusable result in this entry.** Cancer+TSS adjustment collapses every one of those
pooled figures to ~0.46–0.53.

### 3b. TP53 across arms and seeds — the one reproducible signal

| artifact | raw | raw CI | adj | **adj CI** | chance |
|---|---:|---|---:|---|---:|
| d2_h_seed42 | 0.5791 | [0.5437, 0.6134] | **0.5912** | **[0.5618, 0.6274]** | 0.5310 |
| d2_i_seed42 | 0.5398 | [0.5106, 0.5834] | **0.5827** | **[0.5566, 0.6307]** | 0.5303 |
| d2_h_seed43 | 0.5937 | [0.5540, 0.6273] | **0.6017** | **[0.5687, 0.6358]** | 0.5318 |
| d2_h_seed44 | 0.5570 | [0.5280, 0.5964] | **0.5847** | **[0.5550, 0.6210]** | 0.5305 |

**All four adjusted intervals exclude measured chance.** Point estimates 0.583–0.602. Adjustment
*raises* the within-cancer figure every time while collapsing pooled from ~0.69 to ~0.58 — removing
the cancer-type direction removes a nuisance axis that was actively hurting within-cancer separation.

**Honest grading.** My pre-registered rule was "PASS if the CI overlaps the band". These CIs do
overlap [0.60, 0.80], so by the letter of the rule TP53 passes. But the **point estimate is below the
band in three of four arms**, and it overlaps only via its upper tail. The correct reading is
**under-recovery**: morphology carries real TP53 information here, at roughly half the effect size the
published per-cancer literature reports (Kather 0.60–0.78; Noorbakhsh 0.65–0.80). I am not claiming
TP53 as a cleared endpoint.

**TP53 per cancer, adjusted, d2_h_seed42** (14 scored; 7 cancers skipped for minority < 5):

| SKCM | BLCA | LIHC | STAD | KIRP | SARC | COAD | CESC | PAAD | ESCA | MESO | READ | KICH | ACC |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.633 | 0.678 | 0.675 | 0.596 | 0.504 | 0.499 | 0.605 | 0.422 | 0.647 | 0.524 | 0.493 | 0.594 | 0.688 | 0.705 |

Above 0.60 in 7 of 14 (SKCM, BLCA, LIHC, COAD, PAAD, KICH, ACC); at or below chance in KIRP, SARC,
CESC, MESO. The four weakest have tiny minority classes (KIRP 5 positives, CESC 9, MESO 8) and should
not be over-read. **The signal is not uniform across lineages and must not be reported as pan-cancer.**

**The RNA states are not the ceiling-breaker I expected.** Pre-declared as circular and therefore not
evidence about morphology: `rna_biology` 0.6134 and `full_biology` 0.6139 adjusted — only ~0.02 above
image-only 0.5912. Even *with* RNA input, TP53 within-cancer prediction is ~0.61. The limit here is
the 256-d representation, not the imaging modality.

### 3c. Exploratory — no pre-registered band, graded against measured chance only

| endpoint | n / n+ | adj | within-cancer | CI95 | pooled | chance | verdict |
|---|---|---|---:|---|---:|---:|---|
| **grade_high** | 1209 / 672 | raw | 0.7056 | [0.6741, 0.7377] | 0.6524 | 0.5364 | above chance |
| **grade_high** | 1209 / 672 | **adj** | **0.6777** | **[0.6422, 0.7145]** | 0.6285 | 0.5364 | **above chance, survives adjustment** |
| stage_late | 2061 / 855 | raw | 0.5571 | [0.5237, 0.5844] | 0.6223 | 0.5241 | **covers chance** |
| stage_late | 2061 / 855 | adj | 0.5469 | [0.5177, 0.5733] | 0.5348 | 0.5239 | **FAIL — covers chance** |

**Grade is the strongest readout in this entry and the only one with a comfortable margin.** Unlike
the mutation endpoints, its *within-cancer* value (0.706) **exceeds** its pooled value (0.652), which
is the signature of a genuinely within-cancer property rather than a lineage proxy. Per cancer,
adjusted: **ESCA 0.811, BLCA 0.704, LIHC 0.683, STAD 0.667, PAAD 0.633, CESC 0.628, CHOL 0.495** —
consistent across 6 of 7 cancers rather than driven by one.

**Three reasons grade is not an application contribution as it stands**, all declared before the run:
1. **No pre-registered band.** It cannot pass a band it does not have, and I will not retrofit one.
2. **Morphology→morphology.** Grade is assigned by a pathologist reading the slide. Recovering it is
   close to circular and is *not* evidence the representation carries molecular information.
3. BLCA's label is 267/282 positive (95% high-grade); a near-constant label inflates confidence.

**Stage is a clean negative.** Adjusted 0.5469 [0.5177, 0.5733] against measured chance 0.5239 — the
interval covers chance. Per cancer it scatters 0.393–0.721 with six of fourteen *below* 0.5, which is
what noise looks like. Note the raw pooled 0.6223 vs within-cancer 0.5571: the same lineage inflation
again. **Anatomic stage is not legible in the morphology of the primary, in this representation.**

### 3d. Exploratory driver genes — twelve, no band, adjusted, graded against measured chance only

| gene | n+ | within-cancer (adj) | CI95 | pooled | measured chance | cancers | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| ATM | 165 | 0.6135 | [0.5563, 0.6659] | 0.6125 | 0.5468 | 10 | above chance |
| KMT2D | 317 | 0.5997 | [0.5516, 0.6380] | 0.5913 | 0.5334 | 13 | above chance |
| ARID1A | 243 | 0.5981 | [0.5485, 0.6480] | 0.5854 | 0.5412 | 9 | above chance |
| SMAD4 | 113 | 0.6544 | [0.5492, 0.6967] | 0.5943 | 0.5708 | 7 | covers chance |
| NRAS | 127 | 0.6171 | [0.4694, 0.6841] | 0.5292 | 0.5854 | 4 | covers chance |
| FBXW7 | 140 | 0.6030 | [0.5349, 0.6665] | 0.5952 | 0.5483 | 7 | covers chance |
| RB1 | 138 | 0.5916 | [0.5427, 0.6476] | 0.6186 | 0.5534 | 6 | covers chance |
| BRCA2 | 130 | 0.5800 | [0.5143, 0.6555] | 0.5802 | 0.5518 | 6 | covers chance |
| NF1 | 143 | 0.5701 | [0.5022, 0.6308] | 0.5758 | 0.5530 | 7 | covers chance |
| CTNNB1 | 145 | 0.5679 | [0.5074, 0.6514] | 0.5957 | 0.5619 | 6 | covers chance |
| PTEN | 116 | 0.5396 | [0.4822, 0.6181] | 0.5543 | 0.5547 | 10 | covers chance |
| CDKN2A | 112 | 0.4768 | [0.4130, 0.5517] | 0.4614 | 0.5682 | 6 | covers chance |

**Nine of twelve cover measured chance.** Three (ATM, KMT2D, ARID1A) clear it. They are **not
promoted** — the pre-registration forbids moving an endpoint from exploratory to confirmatory after
the fact, they have no band, and with fourteen exploratory endpoints tested and no multiplicity
correction, three clearing at this margin is close to what one should expect anyway.

One observation offered as a **hypothesis, not a claim**: the four genes that clear chance across many
cancers — TP53, ATM, KMT2D, ARID1A — are all large genes whose mutation correlates with tumour
mutational burden, and TMB plausibly tracks nuclear atypia, which is visible. The same explanation
would cover why *grade* is the strongest readout in the entry. If that is what is happening, this
representation carries **one** morphological axis of roughly grade/atypia strength, and every
"mutation prediction" above is that axis re-expressed — not gene-specific information. Testing that
would mean checking whether the probes for TP53, ATM, KMT2D, ARID1A and grade are mutually
predictable and collapse onto a single direction. **Not tested here; recorded as the obvious next
experiment.**

---

## 4. What this means

* **The pooled-vs-within-cancer gap is the transferable finding.** Across BRAF, KRAS, APC, PIK3CA,
  TP53 and stage, the pooled pan-cancer AUROC overstates the within-cancer AUROC by 0.07–0.27. On a
  leave-cancer-types-out split with 21 test lineages, any pan-cancer "molecular readout from H&E"
  claim that does not adjust for cancer type is measuring lineage. This generalises the artifact
  caught in the spatial work to a new statistic and a new task family.
* **Chance is not 0.5.** Measured within-cancer chance ran 0.524–0.596 depending on endpoint and class
  balance. dMMR's was **0.596** — grading it against 0.5 would have turned a null into a "hit".
  Four of the five failures fail *specifically* because their interval covers measured chance, and
  two of them (PIK3CA raw/adj, BRAF raw) would have looked like passes against an assumed 0.5.
* **The split is the binding constraint on this project's clinical claims.** The best-precedented
  endpoints (ER, PR, PAM50, IDH) are unmeasurable held-out by construction. Any future clinical claim
  needs a patient-wise or site-wise split, not this one.
* **A negative that is worth as much as the positive:** with 2,686 held-out patients, correct labels,
  a pre-registered protocol and a positive control that passed on the same code path, four clinically
  actionable driver mutations and stage are **at chance** from morphology once cancer type is
  adjusted for.

## 5. Threats to this entry's own conclusions

* **The representation is 256-d and was never trained for these endpoints.** These are floors for
  *this* representation, not for H&E. A dedicated tile-level model would likely do better; the
  literature bands come from such models. The honest statement is "this representation does not carry
  it", not "morphology does not carry it".
* **TP53's marginal pass depends on the overlap rule.** Under a stricter rule (point estimate inside
  the band) it would fail. The rule was fixed in advance and I am reporting against it, with the
  marginality stated rather than buried.
* **Multiplicity.** Six banded endpoints plus fourteen exploratory ones. Grade, ATM, KMT2D and ARID1A
  clear chance in the exploratory set and are *not* promoted; per the pre-registration, no endpoint
  moves from exploratory to confirmatory after the fact. Four of fourteen clearing, uncorrected, is
  close to what chance would give at this margin.
* **The endpoints that clear may all be one axis.** See §3d: TP53, ATM, KMT2D, ARID1A and grade could
  plausibly be a single atypia/TMB-like direction rather than five readouts. If so, the honest count
  of distinct capabilities here is **one**, not five. Untested.
* **dMMR is uninformative, not negative.** 147 patients in one cancer, 42 positive, measured chance
  0.596. It had no power to detect the published 0.74. This was declared underpowered in advance.

## 6. Files

Pre-registration `NOTEBOOK_ENTRIES/PREDECLARED_E1_clinical_endpoints_20260804T1745Z.md`,
`configs/PREREG_E1_clinical_endpoints.json` (commit f505fdc, before the run).
Labels `.../e1_endpoints/inputs/e1_endpoint_labels.parquet` (MC3 + PanCan clinical mirror);
MC3 calls `.../inputs/mc3_calls.pkl`; results `.../e1_endpoints/results/<endpoint>/`;
logs `.../run_e1.log`, `.../run_tp53_rep.log`, `.../run_tp53_h42.log`.
Suite at the run commit: **407 passed, 27 errors**, all 27 in `v2/tests/test_p2_figures.py` from
absent matplotlib, as documented. No code was modified by this work.
