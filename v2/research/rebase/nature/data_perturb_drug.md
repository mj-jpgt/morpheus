# Perturbation & drug-response data scout — validating an interventional claim without wet lab

**Scout date:** 2026-07-29 · **Tools:** WebFetch + direct `curl` (WebSearch exhausted).
**Verification policy:** every count below is either (a) *measured by me this session* by downloading
the file and counting, marked **[MEASURED]**, or (b) resolved to a PMID/DOI/URL that returned 200,
marked **[VERIFIED]**. Anything I could not confirm is marked **COULD-NOT-VERIFY** and is not
presented as fact. Nothing here is cited from memory.

---

## 0. READ THIS FIRST — three access routes changed in the last 12 months

These will silently break any pipeline written against older documentation.

| Resource | What changed | Status today (2026-07-29) | Workaround |
|---|---|---|---|
| **GDSC — `cancerrxgene.org`** | **Site retired.** Every path returns **HTTP 410 Gone** [MEASURED], serving a notice page: *"The Genomics of Drug Sensitivity in Cancer datasets and features have now been integrated into the Sanger DepMap where users can explore them using enhanced functionality… Download the GDSC1 & GDSC2 datasets from the Cell Model Passports."* | Data itself is **fine and open** | Two live routes, both tested 200 OK: `https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/` (rel 8.4) and `https://cog.sanger.ac.uk/cmp/download/<filename>` (rel 8.5 metadata + 27Oct23 dose-response) |
| **CMap / LINCS — `clue.io`** | **Retired 2026-01-31.** Verbatim banner still served: *"Effective January 31, 2026, the clue.io site, including its tools, will be retired due to financial constraints… All publicly available data will still be accessible at GEO."* [MEASURED] | The **S3 bucket is still up and serving** — I downloaded 465 MB from it today | **Mirror `s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/` NOW.** The LINCS-2020 rebuild (1.2 M signatures) does **not** appear to have a GEO accession — GEO only has the older GSE92742/GSE70138 builds. This is a genuine data-loss risk. |
| **DepMap — `depmap.org/portal`** | Portal is now behind a **Cloudflare Turnstile** challenge. `curl` and WebFetch both get a `DepMap — Verification` HTML page (HTTP 200, 5,175 bytes, loads `challenges.cloudflare.com/turnstile/v0/api.js`) [MEASURED]. Portal text: *"Need DepMap data in bulk? … please don't scrape the portal."* | **No scriptable bulk API.** Downloads require a human browser session. | Two scriptable back doors verified below: **figshare** (Chronos gene-effect matrix, CC BY 4.0) and a third-party **Zenodo** freeze of 25Q3. |

---

## 1. Tier-1: the datasets that carry the argument

### 1.1 DepMap — CRISPR gene dependency (Broad)

* **Exact name:** DepMap Public **26Q1** (announced 2026-04-01 on `forum.depmap.org`, topic 4606) [VERIFIED]. Prior release 25Q3.
* **Size [MEASURED]:** I downloaded the 26Q1 Chronos `gene_effect.csv` (431,100,903 bytes) and counted:
  **1,208 cell-line models × 18,531 genes.** That is the interventional matrix — a Chronos gene-effect
  score (≈0 = no effect, ≈−1 = median common-essential) for every gene × every screened line.
* **Model registry [MEASURED]** (from 25Q3 `Model.csv`, 2,132 rows — the registry is a *superset* of the
  1,208 CRISPR-screened lines): **35 OncotreeLineages, 96 primary diseases.**
  Lymphoid 263 · Lung 260 · Skin 150 · CNS/Brain 127 · Esophagus/Stomach 104 · Bowel 99 · Breast 96 ·
  Head&Neck 94 · Soft Tissue 93 · Bone 93 · Myeloid 88 · Ovary/FT 75 · Kidney 74 · Pancreas 68 ·
  PNS 60 · Biliary 45 · Fibroblast 45 · Uterus 44 · Bladder 39 · Pleura 36 · Liver 29 · Cervix 26 ·
  Thyroid 25 · Eye 24 · **Prostate 15** · Testis 12 · Ampulla 5 · Vulva/Vagina 5 · (151 models are
  annotated `Non-Cancerous`). PrimaryOrMetastasis: 1,135 primary / 630 metastatic / 346 blank.
* **Modalities & PAIRING:** all modalities are keyed on a single `ModelID` (`ACH-xxxxxx`), so
  **CRISPR gene effect ⟷ RNA-seq TPM ⟷ WES/WGS mutations ⟷ CN ⟷ RPPA ⟷ MS proteomics ⟷ Olink
  proteomics ⟷ PRISM/CTRP drug AUC are all matched on the same physical cell line.** This is the
  single most valuable property of DepMap for your purpose: it is a *fully paired* multi-omic +
  interventional panel. `Model.csv` also carries `SangerModelID` and `COSMICID`, which is your join
  key to GDSC and Project Score, and `RRID`/`CCLEName` for everything else.
  26Q1 added: Olink Explore HT proteomics on **161 lines across 24 lineages**, and the harmonised
  Sanger MS proteomics matrix (Gonçalves et al. 2022, re-hosted) [VERIFIED from release notes].
* **Licence / access:**
  * Portal (human browser only, Turnstile): `https://depmap.org/portal/data_page/?tab=allData`
  * **Scriptable:** Chronos parameters + `gene_effect.csv` on figshare, **CC BY 4.0**,
    DOI `10.6084/m9.figshare.31660582` → direct file
    `https://ndownloader.figshare.com/files/62677015` (431 MB) [MEASURED, 200 OK].
    Sibling files: `t0_offset.csv`, `library_effect.csv`, `guide_efficacy.csv`, `replicate_efficacy.csv`.
  * **Scriptable (omics + model table):** third-party frozen mirror, **CC BY 4.0**,
    DOI `10.5281/zenodo.20355477` "DepMap immutable 25Q3 files" (M. Kovaliov, TAU) →
    `OmicsExpressionTPMLogp1HumanAllGenesStranded.csv` (1.11 GB), `OmicsProfiles.csv`, `Model.csv`
    [MEASURED, 200 OK]. Useful precisely because DepMap silently rewrites past releases.
  * Bioconductor `depmap` package ships up to **24Q2** [VERIFIED] — i.e. two years stale; don't use it
    if you need current lineage coverage.
* **Biggest limitation:** **the dependency is measured in 2-D monoculture with no immune, stromal, or
  vascular compartment, and the lineage census is badly skewed** — 15 prostate lines vs 260 lung.
  Any claim of the form "H&E → predicted dependency" is untestable in DepMap for prostate, thyroid,
  uterus, cervix or liver at meaningful n. Concretely: **your validation cohort is effectively
  lung / skin / CNS / bowel / breast / haematologic, and nothing else.**

### 1.2 GDSC1 + GDSC2 — drug dose-response (Sanger)

* **Exact name:** GDSC release **8.4** (`current_release`, dated 2022-07-26); a **8.5** metadata refresh
  exists on the Cell Model Passports bucket (`screened_compounds_rel_8.5.csv`, last-modified 2026-01-14).
* **Size [MEASURED]** — I downloaded and counted both fitted-dose-response tables:
  | | rows (fitted curves) | cell models | drug IDs | drug names | TCGA descriptors |
  |---|---|---|---|---|---|
  | **GDSC1** | 333,161 | 970 | 402 | 378 | 32 |
  | **GDSC2** | 242,036 | 969 | 295 | 286 | 33 |
  Compound annotation table (8.4 **and** 8.5, byte-different but identical in content shape):
  **621 DRUG_IDs / 542 unique drug names**, each with `TARGET` and `TARGET_PATHWAY`.
* **⚠️ The "newer" release is not newer [MEASURED]:** `GDSC2_fitted_dose_response_27Oct23.xlsx`
  (22.1 MB, re-uploaded 2026-02-27) contains **exactly 242,036 rows / 969 models / 295 drugs** — the
  same numbers as the July-2022 8.4 CSV, with columns renamed (`CANCER_TYPE` for `TCGA_DESC`).
  **GDSC drug-response data has not grown since July 2022.**
* **Cancers — GDSC2 unique cell lines by TCGA type [MEASURED]:** UNCLASSIFIED 181 · LUAD 62 · SCLC 59 ·
  SKCM 54 · BRCA 51 · COREAD 46 · HNSC 39 · ESCA 35 · GBM 34 · OV 34 · DLBC 34 · NB 32 · KIRC 32 ·
  PAAD 29 · LAML 26 · ALL 26 · STAD 24 · MESO 21 · BLCA 18 · MM 17 · LGG 17 · THCA 16 · LIHC 15 ·
  CESC 14 · LUSC 14 · LCML 10 · UCEC 9 · **PRAD 6** · MB 4 · CLL 2 · ACC 1.
* **Modalities & PAIRING:** GDSC ships **only** `LN_IC50`, `AUC`, `RMSE`, `Z_SCORE` per (model, drug).
  It carries **no omics of its own.** Pairing is via `SANGER_MODEL_ID` / `COSMIC_ID` → Cell Model
  Passports → DepMap `Model.csv`. Drug → gene target is via the `PUTATIVE_TARGET` / `TARGET` column,
  which is how you close the loop from a *genetic* perturbation prediction to a *pharmacologic* readout.
* **Licence / access:** open, no registration. Verified live:
  * `https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/` — 13 files, incl.
    `GDSC1_fitted_dose_response_24Jul22.csv` (52 M), `GDSC2_…` (38 M), raw data (743 M / 2.0 G),
    `Cell_Lines_Details.xlsx`, ANOVA result workbooks, `GDSCtools_mobems.zip`.
  * `https://cog.sanger.ac.uk/cmp/download/GDSC2_fitted_dose_response_27Oct23.xlsx` etc.
    (bucket **listing** is 403, but **direct object GETs are 200** — you must know the filename).
* **Biggest limitation:** **frozen since 2022, and the drug→target map is "putative".** More
  damagingly for you: `PUTATIVE_TARGET` is a *nominal* annotation, so "predicted dependency on gene G
  → sensitivity to drug D targeting G" inherits every polypharmacology error in that column. And with
  6 prostate and 9 endometrial lines, per-lineage validation is statistically empty outside the big five.

### 1.3 PRISM Repurposing — 4,518 drugs × 578 lines

* **Exact name:** PRISM Repurposing secondary screen. **Corsello SM … Golub TR,
  *Nature Cancer* 2020;1(2):235-248, DOI 10.1038/s43018-019-0018-6, PMID 32613204** [VERIFIED].
* **Size [VERIFIED, verbatim abstract]:** *"growth inhibitory activity of **4,518 drugs** tested across
  **578 human cancer cell lines**"* — barcoded pooled viability screening.
* **Modalities & PAIRING:** viability AUC/logfold only, keyed on DepMap `ModelID` → therefore paired to
  the full DepMap omics + CRISPR stack. This is the **widest chemical space paired to CRISPR
  dependency** that exists open-access.
* **Licence / access:** distributed inside DepMap downloads (`https://depmap.org/repurposing`,
  cited in the paper). **Subject to the Turnstile problem in §0** — no scriptable route found this
  session. **COULD-NOT-VERIFY** a direct machine-fetchable URL for the current PRISM release.
* **Biggest limitation:** pooled barcoded viability means slow-growing lines and lines with unusual
  media requirements drop out; and it is a *viability-only* readout — no transcriptional response,
  so it cannot arbitrate mechanism, only outcome.

### 1.4 CTRPv2 — 481 compounds, connectivity-oriented

* **Exact name:** Cancer Therapeutics Response Portal v2. **Seashore-Ludlow B … Schreiber SL,
  *Cancer Discovery* 2015;5(11):1210-23, DOI 10.1158/2159-8290.CD-15-0235, PMID 26482930** [VERIFIED];
  original **Basu A … Schreiber SL, *Cell* 2013;154(5):1151-1161, DOI 10.1016/j.cell.2013.08.003,
  PMID 23993102** [VERIFIED].
* **Size:** the abstracts do not state the compound/line counts verbatim, so I will not assert the
  commonly-quoted "481 compounds × 860 lines" — **COULD-NOT-VERIFY this session.**
* **Pairing:** CTRP AUC is re-distributed **inside the DepMap portal** (confirmed by an active
  `forum.depmap.org` thread "CTRP AUC values interpretation", 2026-04-17) [VERIFIED], hence keyed on
  DepMap `ModelID` and paired to CRISPR + omics.
* **Biggest limitation:** CTRP's chemical library is heavily weighted to probe compounds and
  tool molecules, not clinical agents — poor for anything you want to argue is clinically actionable.

### 1.5 Sanger Project Score — the *second, independent* CRISPR dependency map

* **Exact name:** Project Score (Sanger DepMap). **Behan FM … Garnett MJ, *Nature* 2019;568:511-516,
  DOI 10.1038/s41586-019-1103-9, PMID 30971826** [VERIFIED]; database paper **Dwane L … Garnett MJ,
  *NAR* 2021;49(D1):D1365-D1372, DOI 10.1093/nar/gkaa882, PMID 33068406** [VERIFIED], which states
  *"the fitness effect of **18,009 genes** tested across **323 cancer cell models**."*
* **Size [MEASURED]:** I downloaded `essentiality_matrices.zip` (241,487,762 bytes) and counted —
  **17,995 genes × 325 cell lines**, shipped as **seven aligned matrices**:
  `00_logFCs.tsv`, `01_corrected_logFCs.tsv` (CRISPRcleanR), `01a_qnorm_corrected_logFCs.tsv`,
  `02a_BayesianFactors.tsv` (BAGEL), `02b_MageckFDRs.tsv`, `03_scaledBayesianFactors.tsv`,
  `04_binaryDepScores.tsv`.
* **Licence / access:** open, direct, **no Turnstile** —
  `https://cog.sanger.ac.uk/cmp/download/essentiality_matrices.zip` [MEASURED, 200 OK].
  Portal `https://score.depmap.sanger.ac.uk/` is a JS SPA (unusable headlessly).
* **Why this matters more than its size suggests:** it is a **different lab, different library,
  different scoring pipeline, overlapping cell lines**. That makes it the *only* clean way to
  measure how much of your predicted-dependency signal is real biology versus Broad-pipeline
  batch structure — a reviewer will ask this, and DepMap alone cannot answer it.
* **Biggest limitation:** 325 lines is a quarter of DepMap's, and only genes present in the KY
  library are covered; cross-map comparison requires restricting to the intersection.

### 1.6 LINCS L1000 / CMap — the largest *genetic-and-chemical* perturbation transcriptome panel

I parsed the full CMap-LINCS-2020 signature metadata (465 MB) this session. All numbers **[MEASURED]**:

* **Level-5 signatures: 1,201,944**, of which **423,422 are exemplar signatures** (`is_exemplar_sig=1`)
  — the exemplar subset is what you actually want.
* **83,588 unique `pert_id`s across 248 cell lines.**
* **Perturbation classes:**
  | class | signatures | unique perturbagens | cell lines |
  |---|---|---|---|
  | `trt_cp` (compound) | 720,216 | 34,418 compounds | 230 |
  | `trt_sh*` (shRNA KD) | 238,351 | 8,767 genes | 22 |
  | **`trt_xpr` (CRISPR KO)** | **140,945** | **5,158 genes** | **27** |
  | `trt_oe` (overexpression) | 34,171 | — | — |
  | `trt_lig` (ligand) | 7,546 | — | — |
  | controls (`ctl_*`) | 59,978 | — | — |
* **CRISPR-KO by line:** ES2 14,567 · A549 14,556 · U251MG 14,246 · A375 11,344 · AGS 10,643 ·
  HT29 10,643 · YAPC 10,643 · BICR6 10,643 · PC3 10,641 · MCF7 8,861 · HS944T 3,603 · then a tail of
  ~2,514-signature lines (KELLY, SNGM, HCC1806, KYSE30…).
* **Compound annotation:** `compoundinfo_beta.txt` — 39,321 rows, **34,419 `pert_id`s, 891 distinct
  targets, 658 distinct MoAs**, with canonical SMILES and InChIKey.
* **Cell annotation:** `cellinfo_beta.txt` — 240 entries: **188 tumour, 23 normal, 29 pooled**;
  21 lineages but **82 entries have `cell_lineage = unknown`**. Lung 32 · haematopoietic 21 ·
  large intestine 19 · ovary 15 · breast 11 · endometrium 9 · skin 9 · CNS 7 · **prostate 6** ·
  urinary tract 5 · bone/kidney/soft tissue/liver 4 each · stomach 3 · pancreas/cervix/placenta 1 each.
  `ccle_name` column gives the join to DepMap.
* **Modalities & PAIRING:** L1000 measures **978 landmark genes directly**, inferring ~12,328 more.
  Pairing to DepMap is via `ccle_name`/`cell_iname`. Crucially, **`trt_xpr` + `trt_cp` in the same
  cell line lets you connect a genetic perturbation signature to a chemical one within a single
  context** — this is the "genetic → pharmacologic" bridge, measured, not modelled.
* **Licence / access:** public, no registration.
  * **CMap LINCS 2020 (recommended, richest):** `https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/`
    — `siginfo_beta.txt` (465 MB), `cellinfo_beta.txt` (38 KB), `compoundinfo_beta.txt` (4.6 MB),
    plus the level-5 GCTX matrices. **All 200 OK today. See §0 — mirror it.**
  * **GEO (durable):** GSE92742 (Phase I, 1,319,138 samples, GPL20573; Level 5 COMPZ 19.9 Gb)
    and GSE70138 (Phase II, 354,123 samples; ~307 GB total) [VERIFIED].
* **Biggest limitation:** **brutal cell-line concentration.** The top 9 lines (MCF7, PC3, A549, A375,
  HT29, VCAP, HA1E, HEPG2, HCC515) account for **59.4 % of all 1.2 M signatures** [MEASURED], and the
  CRISPR-KO arm exists in only **27** lines. So LINCS gives you *chemical* breadth across 230 lines but
  *genetic* breadth across essentially a dozen. Compounding this: L1000 directly measures only 978
  genes; the remaining ~94 % of the "transcriptome" is inferred, which makes it a weak ground truth
  for anything sensitive to genes outside the landmark set.

### 1.7 Tahoe-100M — the largest single-cell perturbation atlas, and the biggest single win here

* **Exact name:** **Tahoe-100M** (Vevo Therapeutics / Arc). Preprint: *"Tahoe-100M: A Giga-Scale
  Single-Cell Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling"*,
  Zhang J, … Goodarzi H, **Yu J**. bioRxiv **10.1101/2025.02.20.639398**, v3 dated 2025-05-10,
  category Genomics. **Still a preprint — not journal-published as of today** [VERIFIED via bioRxiv API].
* **Size [MEASURED from the HF dataset config]:** **95,624,334 cells**; 1,693,653,078,843 bytes
  (**1.69 TB**) uncompressed, **337.6 GB** download; 14 plates × 96 wells = **1,344 samples**.
* **⚠️ Two numbers in the headline are softer than they look [MEASURED]:**
  * *"1,100 small-molecule perturbations"* → the shipped `sample_metadata.parquet` has
    **1,138 unique `drugname_drugconc` combinations but only ~379-380 unique drug names**.
    So it is **~379 distinct chemicals at ~3 doses each**, not 1,100 chemicals. Quote it that way.
  * *"50 cancer cell lines"* (card + preprint) → the shipped `cell_line_metadata.parquet` enumerates
    **102 unique `cell_name` / 99 unique DepMap IDs**. The metadata table appears to be a superset of
    the released panel. **Confirm against `obs_metadata` before quoting a line count.**
* **Organ distribution of the metadata panel [MEASURED]:** Lung 27 · Bowel 25 · Pancreas 11 · Skin 10 ·
  Breast 7 · Esophagus/Stomach 6 · CNS/Brain 4 · Liver 3 · Uterus 3 · and 1 each for
  Vulva/Vagina, Bladder, Ovary/FT, Cervix, PNS, Kidney. **Solid-tumour dominant** — this is the
  single best answer to "solid-tumour perturbation atlas".
* **Modalities & PAIRING — this is the important part:**
  * Cells are **pooled across cell lines in the same well** and demultiplexed by SNP, so drug effect
    and cell-line identity are **measured in the same physical well** — batch confounding between
    context and perturbation is structurally eliminated. Each cell row carries
    `cell_line_id`, `drug`, `moa-fine`, `canonical_smiles`, `pubchem_cid`, `plate`, `sample`.
  * `cell_line_metadata.parquet` carries **`Cell_ID_DepMap` (ACH-…) and `Cell_ID_Cellosaur`** →
    **direct, keyed join to DepMap CRISPR gene effect, CCLE omics, GDSC and PRISM.** [MEASURED: I read
    the parquet; e.g. A549 = ACH-000681 = CVCL_0023.]
  * It also ships per-line **driver-gene annotation** (`Driver_Gene_Symbol`, `Driver_VarZyg`,
    `Driver_VarType`, `Driver_ProtEffect_or_CdnaEffect`, `Driver_Mech_InferDM` LoF/GoF,
    `Driver_GeneType_DM` Oncogene/Suppressor) — e.g. A549: CDKN2A Hom Deletion LoF, KRAS p.G12S GoF.
  * **`metadata/pseudobulk_differential_expression/` is the practical entry point** — 1,026 parquet
    shards. I opened shard 0 [MEASURED]: 3,986,181 rows, columns
    `gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, plate, n_cells_trt, n_cells_ctrl,
    Cell_ID_Cellosaur, Cell_ID_DepMap, drug, concentration, concentration_unit, Cell_Name_Vevo`;
    that shard = **A549 / ACH-000681 / plate 1 / 64 drugs**. **You can get a full DESeq2-style
    (cell line × drug × dose × gene) log2FC table without touching the 337 GB cell-level data.**
* **Licence / access:** **CC0-1.0 — public domain, no restrictions at all.**
  `https://huggingface.co/datasets/tahoebio/Tahoe-100M` (41,823 downloads).
  Direct file example (tested 200 OK):
  `https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/cell_line_metadata.parquet`.
  Companions: `tahoebio/Tahoe-x1-embeddings` (Apache-2.0), `tahoebio/tahoe-de-rhaister`
  (CC0, DE summary statistics), plus community SLAF-format mirrors [VERIFIED via HF API].
* **Biggest limitation:** **the perturbations are chemical, not genetic.** If your response space is
  built from genetic (CRISPRi) Perturb-seq, Tahoe cannot supply held-out *genetic* perturbations —
  it can only test whether the *context* dimension (cell line → response) transfers. Secondarily,
  cells are grown as **mixed-cell-line spheroids in 96-well plates** at a single 24 h-ish timepoint;
  mean transcript counts are low (~1,900-2,200 tscp/cell, ~1,200-1,400 genes/cell in the samples I
  inspected [MEASURED]), so per-cell resolution is shallow and pseudobulking is essentially mandatory.

### 1.8 Genome-scale Perturb-seq (the response space itself)

* **Replogle JM … Weissman JS, *Cell* 2022;185(14):2559-2575.e28, DOI 10.1016/j.cell.2022.05.013,
  PMID 35688146** [VERIFIED, verbatim abstract]: *"genome-scale Perturb-seq targeting all expressed
  genes with CRISPR interference (CRISPRi) across **>2.5 million human cells**."*
  Cell lines are **K562 (CML) and RPE1 (immortalised retinal epithelium)**.
  **Biggest limitation for you: neither is a solid tumour.** K562 is a leukaemia line; RPE1 is not
  cancer at all. Any "solid-tumour patient → this response space" claim crosses a lineage gap that
  is *itself* the thing under test.
* **Companion analysis resource:** **Nadig A … O'Connor LJ, "Transcriptome-wide analysis of
  differential expression in perturbation atlases", *Nature Genetics* 2025;57(5):1228-1237,
  DOI 10.1038/s41588-025-02169-3, PMID 40259084** [VERIFIED] — reprocessed DE statistics over
  perturbation atlases; use this rather than re-deriving DE yourself.

### 1.9 Solid-tumour Perturb-seq — what actually exists (all [VERIFIED] with verbatim abstracts)

| Dataset | Model | Scale | Why it matters | Limitation |
|---|---|---|---|---|
| **Ursu O … Boehm JS, *Nat Biotechnol* 2022;40:896-905, DOI 10.1038/s41587-021-01160-7, PMID 35058622** | **Lung cancer cells** (solid) | *"**200 TP53 and KRAS variants** … in over **300,000 single lung cancer cells**"* | The only large **variant-level** Perturb-seq in a solid tumour line; gives GoF/LoF/dominant-negative continua rather than binary KO | Two genes only; variant impact ≠ dependency |
| **Liu SJ … Gilbert LA, *Genome Biology* 2024;25, DOI 10.1186/s13059-024-03404-6, PMID 39375777** | **In vivo GBM** (mouse, intracranial CED) | multiplex in vivo CRISPRi perturb-seq across malignant **and** microenvironment cells, ± radiotherapy | The **only** in-vivo cancer Perturb-seq I could verify with a TME compartment and a therapy interaction; explicitly shows *"radiotherapy rewires transcriptional responses to genetic perturbations in an in vivo-dependent manner"* | Mouse model; not patient tissue; modest gene-set size |
| **Spisak S … Sethi NS, *Nat Commun* 2024;15, DOI 10.1038/s41467-024-46285-w, PMID 38472198** | Colorectal | dual endogenous reporter + perturb-seq for stem/differentiation regulators | Solid-tumour, differentiation-state phenotype | Targeted library, single phenotype axis |
| **Hou J … Peng W, *NAR Cancer* 2022;4, DOI 10.1093/narcan/zcac038, PMID 36518525** | Tumour-intrinsic immune factors | single-cell CRISPR immune screens | Links perturbation to immune phenotype | Small, immune-axis only |

**Key negative finding:** I found **no genome-scale (all-expressed-genes) Perturb-seq in a solid-tumour
cell line** in this search. The Liu 2024 GBM paper explicitly frames in-vivo cancer perturb-seq as
newly enabled. **This gap is real and is the strongest honest novelty argument available** — but it
also means you cannot validate a solid-tumour genetic-response space against a like-for-like ground
truth, and a reviewer will notice.

### 1.10 Harmonised aggregators (use these to avoid re-inventing preprocessing)

* **scPerturb — Peidli S … Sander C, *Nature Methods* 2024;21(3):531-540,
  DOI 10.1038/s41592-023-02144-y, PMID 38279009** [VERIFIED, verbatim]: *"**44 publicly available
  single-cell perturbation-response datasets** with molecular readouts, including transcriptomics,
  proteomics and epigenomics… uniform quality control pipelines and harmonize feature annotations."*
  Also defines **E-distance** (energy statistics) as a perturbation-effect metric — a defensible,
  already-peer-reviewed effect-size measure you can adopt instead of inventing one.
  Access: `scperturb.org`. **Limitation:** harmonisation is at the annotation level, not batch-corrected;
  the 44 datasets are overwhelmingly non-cancer or leukaemia.
* **PerturBase — Wei Z … Liu Q, *NAR* 2025, DOI 10.1093/nar/gkae858, PMID 39377396** [VERIFIED] —
  database for single-cell perturbation data analysis and visualisation. **COULD-NOT-VERIFY** its
  dataset/cell counts this session.

### 1.11 Preclinical models beyond cell lines (the transfer bridge)

* **Novartis PDXE — Gao H … Sellers WR, *Nature Medicine* 2015;21(11):1318-25, DOI 10.1038/nm.3954,
  PMID 26479923** [VERIFIED]: *"High-throughput screening using patient-derived tumor xenografts to
  predict clinical trial drug response"* — the "1×1×1" population design. This is the canonical
  **in-vivo, patient-derived, drug-response** dataset with matched omics.
  **COULD-NOT-VERIFY** exact PDX/treatment counts or a live download URL this session.
  **Limitation:** human stroma is replaced by mouse within a few passages; immune-deficient host, so
  no immune-mediated response.
* **Cell Model Passports — van der Meer D … Garnett MJ, *NAR* 2019;47(D1):D923-D929,
  DOI 10.1093/nar/gky872, PMID 30260411** [VERIFIED]. The Sanger-side model registry and the
  designated new home for GDSC. Machine-fetchable objects at `https://cog.sanger.ac.uk/cmp/download/`
  (e.g. `model_list_20241120.csv`, 925 KB, verified 200 OK); the web UI is a JS SPA.

---

## 2. The actual question: what can validate cell-line → tumour transfer **without wet lab**

This is the part of your claim that is load-bearing, so here is the verified state of the art,
including the result that cuts against you.

### 2.1 The alignment layer — makes cell-line and tumour expression comparable

**Celligner — Warren A … McFarland JM, *Nature Communications* 2021;12:22,
DOI 10.1038/s41467-020-20294-x, PMID 33397959** [VERIFIED, verbatim abstract].
*"An unsupervised alignment method (Celligner)… integrate several large-scale cell line and tumor
RNA-Seq datasets. Although our method aligns the majority of cell lines with tumor samples of the same
cancer type, it also reveals **large differences in tumor similarity across cell lines**. Using this
approach, we identify **several hundred cell lines** from diverse lineages that present a more
**mesenchymal and undifferentiated** transcriptional state and that **exhibit distinct chemical and
genetic dependencies**."*

Read that last clause carefully: **cell lines that fail to align to tumours have systematically
different dependencies.** That is a directly measurable, publishable *confounder check* for you —
and it is the mechanism by which a DepMap-validated prediction can be right in vitro and wrong in
the patient. DepMap ships Celligner distances (an active forum thread "Celligner median distances",
2026-04-30, confirms they are portal-available) [VERIFIED].

### 2.2 The transfer layer — three peer-reviewed, wet-lab-free frameworks

| Method | Citation | What it transfers | Validation cohort — **this is the template you should copy** |
|---|---|---|---|
| **PRECISE** | Mourragui S … Wessels LFA, *Bioinformatics* 2019;35(14):i510-i519, DOI 10.1093/bioinformatics/btz372, PMID 31510654 [VERIFIED] | Linear domain adaptation, preclinical → tumour | TCGA |
| **TRANSACT** | Mourragui SMC … Wessels LFA, *PNAS* 2021;118(49), DOI 10.1073/pnas.2106682118, PMID 34873056 [VERIFIED] | Non-linear consensus space across cell lines, PDX and human tumours | *"**23 drug prediction challenges on The Cancer Genome Atlas and 226 metastatic tumors from the Hartwig Medical Foundation**"*; robust for *"platinum-based chemotherapies, gemcitabine, and paclitaxel"*; recovers known targeted-therapy biomarkers |
| **CellHit** | Carli F … Raimondi F, *NAR* 2025;53(W1):W143-W150, DOI 10.1093/nar/gkaf414, PMID 40377071 [VERIFIED] | GDSC1/2 + PRISM models → patient transcriptomes | Ships **precomputed predictions across all TCGA samples**, batch correction + *"enhanced Celligner methodology"* + Parametric UMAP + SHAP. `https://cellhit.bioinfolab.sns.it/` |

**CellHit is your fastest baseline.** It already publishes per-TCGA-sample drug-sensitivity predictions
derived from GDSC+PRISM+CCLE. If your image-derived predictions do not beat CellHit's
transcriptome-derived ones on the same TCGA samples, the whole H1 hop (image → molecular state) is
adding nothing, and you will have found that out in a week without generating any data.

### 2.3 The result that cuts against the claim — cite it before a reviewer does

**Zhang et al., "Characteristics and prognostic value of potential dependency genes in clear cell
renal cell carcinoma based on… DepMap", *Int J Med Sci* 2021;18(10), DOI 10.7150/ijms.51703,
PMID 33850477** [VERIFIED]. They took 16 DepMap-nominated ccRCC dependency genes into real patient
cohorts; only **GET4 and CRB3** survived as independent prognostic factors, and the paper concludes
that **"a dependency gene validated in cell lines didn't directly represent its role in corresponding
patients" of the same cancer type.**

This is a small paper in a modest journal, but it is the cleanest published statement of your
central risk: **DepMap-validated ≠ patient-relevant, even within-lineage.** Your design must
pre-register how you handle the case where DepMap agreement is high and survival association is null.

### 2.4 Adjacent recent work worth knowing

* **FORGE — Bhattacharjee N … Sengupta D, *Nature Communications* 2026;17(1),
  DOI 10.1038/s41467-026-73977-2, PMID 42259811** [VERIFIED]: *"Gene dependency-informed inference of
  response to targeted cancer therapies"* — joint modelling of drug response **and** gene dependency,
  validated in PDX. This is the closest published thing to "dependency prediction → drug response",
  from omics rather than images. **You must position against it.**
* **STATE — Adduri A, … Roohani YH, bioRxiv 10.1101/2025.06.26.661135** (v2, 2025-07-10) [VERIFIED]:
  transformer trained on *"gene expression data from over 100 million perturbed cells"*, predicting
  perturbation response across contexts, with a cell-type-specific evaluation. This is the incumbent
  model of the response space you are proposing to map patients into.
* **Arc Virtual Cell Challenge benchmark** is built on **H1 human embryonic stem cells**
  [VERIFIED from three independent arXiv abstracts: 2604.13986, 2603.25240, 2511.16954] —
  **not cancer**. Do not let the field's headline benchmark stand in for cancer-context transfer.
* **Evaluation caution:** Marrakchi Y, D'Ascenzo D, Cultrera di Montesano S, arXiv **2607.04595**
  (2026-07-06) [VERIFIED] — *"per-cell accuracy and raw-pseudobulk scores should be used with caution"*;
  proposes a population-level Classifier Discrimination Score. Also Liu Q et al., arXiv **2511.16954**:
  the Perturbation Discrimination Score *"is highly sensitive to the choice of similarity or distance
  measure."* **Your metric choice is itself contestable — declare and justify it.**
* **MIX-Seq — McFarland JM … Tsherniak A, *Nature Communications* 2020;11:4296,
  DOI 10.1038/s41467-020-17440-w, PMID 32855387** [VERIFIED, verbatim]: multiplexed post-perturbation
  scRNA-seq *"across pools of **100 or more cancer cell lines**"*, chemical **and genetic**
  perturbations, SNP-demultiplexed, and — critically — *"enable prediction of long-term cell viability
  from short-term transcriptional responses to treatment."* This is the **published precedent that a
  transcriptional response predicts a viability endpoint across contexts** — i.e. the H3→H4 link.
  It is Broad-authored and DepMap-keyed, so it joins straight to your CRISPR matrix.

---

## 3. A concrete no-wet-lab validation ladder

Ordered by cost. Every input below is open-access and I have verified a live URL for all except
PRISM and PDXE.

1. **Consistency (free, 1 day).** Predicted dependency vs **DepMap 26Q1 Chronos** (1,208 × 18,531),
   stratified by `OncotreeLineage`. Report per-lineage, never pooled — pooled correlation is inflated
   by the common-essential gene block.
2. **Pipeline-independence (free, 1 day).** Repeat against **Sanger Project Score**
   (17,995 × 325, `essentiality_matrices.zip`). Restrict to the model intersection via
   `Model.csv:SangerModelID`. Any signal that does not survive both is a Broad artefact.
3. **Genetic → pharmacologic closure (2 days).** Map predicted dependency on gene *G* to
   **GDSC2** `LN_IC50`/`AUC` for drugs whose `PUTATIVE_TARGET` contains *G* (295 drugs, 969 models),
   then to **PRISM** (4,518 × 578) for breadth. Pre-register the target-mapping rule.
4. **Genetic → transcriptional closure (2 days).** For the **27 lines with LINCS `trt_xpr`**,
   check that your predicted response to KO of *G* matches the measured L1000 CRISPR-KO signature —
   and that the corresponding `trt_cp` compound signature in the *same* line agrees. This is the only
   place you can test genetic-vs-chemical concordance in a matched context, measured.
5. **Context transfer (1 week).** Use **Tahoe-100M pseudobulk DE** (CC0, per cell-line × drug × dose
   log2FC, no 337 GB download needed) to test whether your model's *context* dimension predicts
   which cell line responds — held out by lineage, not at random.
6. **Alignment confound (2 days).** Compute or pull **Celligner** tumour-similarity per cell line;
   show your accuracy does **not** collapse on the "mesenchymal/undifferentiated" poorly-aligned
   subset. If it does, say so — Warren et al. predicts it will.
7. **Patient endpoint (1-2 weeks).** TCGA survival stratified by predicted dependency, **plus**
   head-to-head against **CellHit**'s precomputed TCGA predictions. Report the PMID 33850477 failure
   mode explicitly.
8. **External clinical cohort (weeks + application).** **Hartwig Medical Foundation**, 226 metastatic
   tumours as used by TRANSACT — this is the strongest available human-outcome validation, but access
   is by data-request application, **not** open download. **COULD-NOT-VERIFY** current Hartwig
   application terms this session.

---

## 4. Availability verdict

**Fully open, scriptable, verified live today:**
Tahoe-100M (CC0) · Sanger Project Score essentiality matrices · GDSC 8.4/8.5 via Sanger FTP and
`cog.sanger.ac.uk` · CMap LINCS 2020 via S3 · LINCS GSE92742/GSE70138 via GEO · DepMap 26Q1 Chronos
gene-effect via figshare (CC BY 4.0) · DepMap 25Q3 omics/Model via Zenodo (CC BY 4.0) · scPerturb ·
all Perturb-seq papers' GEO deposits (per-paper).

**Open but browser-only (Cloudflare Turnstile — no API):** the DepMap portal itself, and therefore
**PRISM Repurposing** and **CTRP** in their current releases.

**Retired sites, data survived:** `cancerrxgene.org` (410 Gone) · `clue.io` (retired 2026-01-31).
Neither is paywalled; both simply moved or stopped.

**Controlled access (application required, not paywalled):** Hartwig Medical Foundation.

**COULD-NOT-VERIFY this session:** CTRPv2 exact compound/cell-line counts · a scriptable PRISM
download URL · PDXE model counts and download URL · PerturBase dataset counts · current Hartwig
access terms · any GEO accession for the CMap **LINCS-2020** rebuild (only the older Phase I/II
builds are in GEO).

---

## 5. Provenance of measured numbers

Files downloaded and counted during this session (scratchpad
`…/183cb9cc-8cb0-4a9d-ad6f-8c70cbac5491/scratchpad/`):
`GDSC1_fitted_dose_response_24Jul22.csv` (54,084,519 B) ·
`GDSC2_fitted_dose_response_24Jul22.csv` (39,834,719 B) ·
`GDSC2_fitted_dose_response_27Oct23.xlsx` (22,112,537 B) ·
`screened_compounds_rel_8.4.csv` / `rel_8.5.csv` ·
DepMap 26Q1 Chronos `gene_effect.csv` (431,100,903 B) ·
DepMap 25Q3 `Model.csv` (699,474 B) ·
`essentiality_matrices.zip` (241,487,762 B) ·
`siginfo_beta.txt` (465,242,319 B) · `cellinfo_beta.txt` · `compoundinfo_beta.txt` ·
Tahoe `cell_line_metadata.parquet`, `drug_metadata.parquet`, `sample_metadata.parquet`,
`pseudobulk_differential_expression/train-00000-of-01026.parquet` (91,442,763 B).
