# PRE-REGISTRATION — E1: clinically meaningful endpoints predicted from morphology alone, held-out cancers

**Written:** 2026-08-04 17:45 UTC. **Nothing in this file was informed by any endpoint run.**
Label coverage was audited first (that is a counting operation, not a statistic); no probe, no AUROC,
no permutation has been computed on any endpoint at the time of writing. The run command is given at
the bottom and will be executed *after* this file is committed.

## 0. The structural fact that governs everything below

The D2 split is **leave-cancer-types-out**, not patient-wise. Every cancer type sits wholly in
train/val or wholly in test:

* **train/val (11 cancers, 3,661 patients):** BRCA, LGG, LUSC, LUAD, KIRC, THCA, HNSC, PRAD, UCEC, GBM, OV
* **test (21 cancers, 2,766 patients):** SKCM 306, BLCA 283, LIHC 268, STAD 244, KIRP 211, SARC 210,
  CESC 198, COAD 197, PAAD 123, TGCT 114, ESCA 110, THYM 82, MESO 60, READ 59, UVM 54, KICH 50,
  UCS 49, ACC 48, PCPG 43, CHOL 29, DLBC 28

**This is confirmed and it does constrain the endpoint list, as the handoff asked me to check.** The
consequence is stronger than "ER/PR has 0 test patients":

* **ER/PR and PAM50** are BRCA-only → BRCA is train → **no held-out measurement is possible.**
* **IDH status** is LGG/GBM → both train → **no held-out measurement is possible.**
* **HPV in HNSC** → train. HPV in CESC is test but see §2.

So the ER precedent at 0.8781 cannot be replicated as a held-out result on this split, by anyone,
ever, without a re-split. Endpoints below were chosen *because* they live in the test cancers.

The compensating advantage: because the test partition spans 21 lineages, the within-cancer statistic
is exercised properly here in a way the ER control explicitly could not exercise it (ER was one
cancer, so pooled and within-cancer were identical by construction). Cancer-type adjustment is a real
test on this partition rather than a formality.

## 1. Verification of the earlier "unusable on disk" claim — it was partly a lookup failure

`t17b_known_covariate_20260803T0155Z.md` states: *"there is no mutation table and no consensus-subtype
table on disk"* and that MSI carries *"no MSI-H/MSI-L/MSS calls whatsoever"*. I re-checked rather than
inherited, per instruction. Findings, all three reported as discrepancies:

1. **FALSE — the mutation table is on disk.**
   `/lambda/nfs/geeg/biorag3_persistent_20260711/raw_tcga_v1/mc3.v0.2.8.PUBLIC (1).maf.gz`, 753 MB,
   the MC3 public pan-cancer MAF. Parsed: 3,600,963 variant rows, **10,224 sequenced patients**.
   The filename contains a space and a `(1)`, which is the most likely reason a glob missed it.
   **2,686 of the 2,766 held-out patients have MC3 calls.** TP53 is mutated in 893 of them.
   The claim that mutation status is unavailable is wrong and is corrected here.

2. **PARTLY FALSE — MSI.** The `microsatellite_instability` column is indeed an assay-performed flag
   (test partition: 7 YES / 74 NO) and carries no MSI-H/MSI-L/MSS calls — that part of the claim is
   **confirmed**. But the adjacent column
   `loss_expression_of_mismatch_repair_proteins_by_ihc` carries genuine **mismatch-repair-deficiency
   IHC calls for 147 held-out patients (42 dMMR-positive)**, and
   `loss_expression_of_mismatch_repair_proteins_by_ihc_result` carries per-protein MLH1/MSH2/PMS2/MSH6
   detail for 44. dMMR by IHC is the clinical companion diagnostic for MSI-H and is the same
   actionable readout. So a usable, if small, MSI-family label does exist and was missed.

3. **CONFIRMED — no consensus-subtype table.** No PAM50, no CMS, no TCGA molecular-subtype table found
   anywhere under `/lambda/nfs/geeg`. This claim stands. (And PAM50/CMS would in any case be
   unmeasurable-or-BRCA-only per §0.)

Also confirmed absent for the test partition: **HPV** — `hpv_status_by_p16_testing`,
`hpv_status_by_ish_testing` and `hpv_test` are **all zero-coverage on all 2,766 held-out patients**
despite CESC (n=198) being a test cancer. HPV is therefore dropped for insufficient coverage, not for
a negative result. **HRD** — no HRD/LOH/telomeric-imbalance score table on disk; the raw Gistic2 and
SNP6 segment files are present but deriving an HRD score from them is a new method, not a label
lookup, and would have no pre-registered basis. HRD is dropped.

## 2. Endpoint slate, with held-out coverage measured before any modelling

Coverage counted on the D2 test partition. "cancers scoreable" = cancers with n>=20 and minority
class >=5, the thresholds already hard-coded in `v2/calibra/known_covariate_control.py`
(`_MIN_PER_CANCER=20`, `_MIN_MINORITY=5`).

| endpoint | n held-out labelled | n positive | prevalence | cancers scoreable |
|---|---:|---:|---:|---|
| `mut_TP53` | 2,686 | 893 | 0.332 | **14** |
| `mut_KRAS` | 2,686 | 262 | 0.098 | 9 |
| `mut_BRAF` | 2,686 | 211 | 0.079 | 4 |
| `mut_PIK3CA` | 2,686 | 276 | 0.103 | 10 |
| `mut_APC` | 2,686 | 306 | 0.114 | 8 |
| `mmr_deficient` | 147 | 42 | 0.286 | 1 (COAD) |
| `grade_high` | 1,209 | 672 | 0.556 | 7 |
| `stage_late` | 2,061 | 855 | 0.415 | 14 |

Eleven further genes (ARID1A, KMT2D, PTEN, CTNNB1, NRAS, SMAD4, FBXW7, RB1, ATM, BRCA2, CDKN2A, NF1)
have coverage and are declared **exploratory** in §4.

## 3. Pre-registered bands, from literature retrieved and read today

Every number below was read off a page I fetched. The URL is given. No citation appears here that I
did not open.

### 3a. Primary confirmatory endpoints

**`mut_TP53` — band [0.60, 0.80], point estimate 0.70.**
* Kather JN, Heij LR, Grabsch HI, et al. *Pan-cancer image-based detection of clinically actionable
  genetic alterations.* **Nature Cancer** 2020;1:789–799. PMID 33763651.
  Fetched: <https://pmc.ncbi.nlm.nih.gov/articles/PMC7610412/>
  Quoted: *"the mean cross-validated area under the receiver operating curve (AUROC) for the top eight
  mutations ranged from 0.60 to 0.78 in lung adenocarcinoma; from 0.65 to 0.76 in colorectal cancer;
  from 0.62 to 0.78 in breast cancer and from 0.66 to 0.78 in gastric cancer."* TP53 is stated to be
  *"significantly detected in all four of these cancer types"*.
* Noorbakhsh J, Farahmand S, Foroughi pour A, et al. *Deep learning-based cross-classifications reveal
  conserved spatial behaviors within tumor histological images.* **Nature Communications**
  2020;11:6367. PMID 33311458. Fetched: <https://pmc.ncbi.nlm.nih.gov/articles/PMC7733499/>
  Quoted, TP53: *"Self-cohort predictions (diagonal values) have AUC values ranging from 0.65–0.80 for
  per-slide"*, with LUAD 0.80 and STAD 0.65; cross-tissue TP53 *"AUCs 0.62–0.72 for slides"*.
* Band is the union of the two sources' per-cancer ranges, 0.60 to 0.80.

**`mmr_deficient` (dMMR by IHC) — band [0.66, 0.80], point estimate 0.74.**
* Echle A, Grabsch HI, Quirke P, et al. *Clinical-Grade Detection of Microsatellite Instability in
  Colorectal Tumors by Deep Learning.* **Gastroenterology** 2020;159(4):1406–1416. PMID 32562722.
  Fetched: <https://pmc.ncbi.nlm.nih.gov/articles/PMC7578071/>
  Quoted: *"training the deep learning system on individual cohorts yielded an intra-cohort AUROC of
  0.74 [0.66, 0.80] in the TCGA cohort (n=426)"*.
* Band is exactly that published TCGA interval. Declared in advance as **underpowered**: 147 patients
  in one cancer, 42 positive. The module's own guard requires >=30 per class, which this passes only
  just. A wide CI is expected and a non-informative result here is uninformative, not negative.

### 3b. Secondary endpoints, same band source, same grading rule

**`mut_BRAF` — band [0.64, 0.82], point 0.77.** Kather 2020 (PMC7610412), external validation in
colorectal cancer: *"0.77 (0.64 – 0.82, p<10−5)"*.

**`mut_KRAS`, `mut_PIK3CA`, `mut_APC` — band [0.60, 0.78], point 0.69.** Kather 2020 (PMC7610412),
the stated top-eight-mutation range across the four cancer types. These are members of that
actionable set; no gene-specific figure is printed in the fetched text, so the general range is used
and this is a weaker basis than TP53's, declared as such now.

### 3c. Exploratory — NO pre-registered band

`grade_high` and `stage_late`, plus the eleven further genes in §2.

I searched for a literature band for histologic grade and for stage and **could not retrieve one I was
willing to cite**. Rather than invent a band, these are graded **only against measured chance** and are
reported as exploratory. They cannot "pass a pre-registered band" because they have none; that is
stated now, before the run, so that a good number on grade cannot later be presented as a confirmatory
result.

Two further cautions declared in advance:
* **Grade is not a molecular inference.** Histologic grade is assigned by a pathologist looking at the
  slide. Recovering it from morphology is morphology→morphology and is close to circular. A high value
  is expected and is *not* evidence that the representation carries molecular information.
* **Stage is anatomic extent**, only partly encoded in the primary tumour's morphology.

## 4. Grading rule — fixed now

Inherited unchanged from the ER protocol so the two are comparable:

1. Primary statistic is the **within-cancer AUROC** (size-weighted mean of per-cancer AUROCs), on the
   **test partition**, `--partition test`. Pooled AUROC is reported beside it *specifically so the gap
   is visible*: pooled >> within-cancer means the probe is guessing lineage, which is the artifact
   this project exists to catch.
2. **PASS** requires the 1,000-draw bootstrap CI of the within-cancer AUROC to **overlap the band**.
   CI entirely **below** = under-recovery, **FAIL**. CI entirely **above** = leak/circularity,
   **FAIL** — not a triumph.
3. **Additionally**, and independently, the CI lower bound must **exceed the measured chance level**
   (`null_p95` from 1,000 within-cancer label permutations). Chance is **measured, not assumed to be
   0.5** — the ER run measured 0.542–0.546. An endpoint whose interval covers measured chance is a
   **negative** and is reported as one.
4. **Both raw and cancer+pooled-TSS-adjusted are reported for every endpoint.** The adjusted number is
   the one that counts. An endpoint that passes raw and fails adjusted is a cancer-type-identification
   artifact and is reported as a failure.
5. **Per-cancer AUROCs are reported for every endpoint**, never pooled only.
6. **Image-only means `wsi_biology`.** `rna_biology` and `full_biology` take RNA as input; mutation
   status is correlated with expression, so those states are circular for these endpoints by
   construction. They are run for completeness and **pre-declared as not evidence about morphology**.
   Only `wsi_biology` can pass.
7. **Multiplicity.** Six endpoints carry bands (§3a, §3b). The exploratory set in §3c is not corrected
   and is not confirmatory. No endpoint will be promoted from exploratory to confirmatory after the
   fact.

## 5. Success criterion, stated before the run

Per the handoff: two or three endpoints above their pre-registered band, surviving cancer+site
adjustment, with an interval excluding measured chance, on held-out patients = an application
contribution. One = a replication. **Zero = a finding about how little morphology carries once
adjustment is done properly, and will be reported at the same prominence.**

My prior, recorded now: TP53 is the only primary endpoint I expect to have a real chance, and I
expect it to land at the **bottom** of its band or below, because the D2 representation is 256-d and
was not trained for mutation prediction, and because these are entirely unseen lineages. dMMR I expect
to be uninformative through lack of power.

## 6. Exact run command

```
python -m morpheus.v2.calibra.known_covariate_control \
  --artifacts ~/e0_run/d2_v3/d2_v3_s42/artifacts/d2_h_seed42.npz \
              ~/e0_run/d2_v3/d2_v3_s42/artifacts/d2_i_seed42.npz \
  --covariate-table .../e1_endpoints/inputs/e1_endpoint_labels.parquet \
  --covariate-column <endpoint> --expected-low <lo> --expected-high <hi> \
  --preregistration <this file's JSON sibling> \
  --partition test --n-boot 1000 --n-permutations 1000 --seed 42
```

Statistics come from `v2/calibra/known_covariate_control.py`, unmodified. Nothing is computed inline.
Label table: `e1_endpoint_labels.parquet`, built from MC3 + the PanCan clinical mirror, described in
§1–2 above.
