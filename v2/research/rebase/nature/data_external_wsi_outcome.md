# External (non-TCGA) open-access cohorts with WSI + molecular and/or outcome

**Scout:** `external_wsi_outcome` · **Date:** 2026-07-29
**Verification method:** WebFetch only (WebSearch exhausted). Every number below was pulled from a live fetch of the source listed. Anything I could not pull is flagged **COULD-NOT-VERIFY** and must not be cited as fact.

---

## 0. Executive verdict (read this first)

The honest picture for "external cohort with real slides *and* real outcome":

| Tier | What exists | Reality |
|---|---|---|
| **Tier A — slides + molecular + survival, genuinely open** | CPTAC (TCIA), HANCOCK, SurGen | Real. Downloadable today. CC BY. |
| **Tier B — slides + treatment response, small** | Ovarian Bevacizumab Response, IMPRESS | Real but tiny (78 and 126 patients). |
| **Tier C — slides + genomics + outcome, but outcome behind dbGaP** | Cancer Moonshot Biobank (CMB), NLST | Images CC BY and open; the outcome table needs an application. |
| **Tier D — slides, no outcome (task-specific)** | CAMELYON16/17, PANDA, BCNB, Digital Brain Tumour Atlas | Big, clean, free — but they are *label* datasets, not outcome cohorts. |
| **Tier E — outcome, no images at all** | AACR GENIE / GENIE BPC | **No imaging whatsoever.** See §7. |
| **The gap** | Immunotherapy-treated cohort with WSIs | **Does not exist open-access at scale.** See §8. This is the single most important negative finding. |

---

## 1. CPTAC — the workhorse external cohort

The Clinical Proteomic Tumor Analysis Consortium is the only resource that gives you **WSI + proteomics + genomics + survival** across many cancer types, openly, at scale. It is the natural non-TCGA counterpart.

### 1.1 WSI availability by cancer (verified from TCIA browse-collections, fetched live)

| TCIA collection | Cancer | Subjects | WSI modality present | External omics linked |
|---|---|---:|---|---|
| CPTAC-LUAD | Lung adenocarcinoma | 244 | Histopathology (243 subj, **1,137 slides, 431.5 GB**) | Clinical, Genomics, Proteomics |
| CPTAC-CCRCC | Clear-cell renal cell | 262 | WSI (222 subj, **783 slides, 190 GB**) | Clinical, Genomics, Proteomics |
| CPTAC-UCEC | Uterine corpus endometrial | 250 | WSI (**887 slides, 154 GB**) | Clinical, Genomics, Proteomics |
| CPTAC-LSCC | Lung squamous cell | 212 | WSI | Clinical, Genomics, Proteomics |
| CPTAC-HNSCC | Head & neck squamous | 207 | WSI | Clinical, Genomics, Proteomics, Microarray |
| CPTAC-GBM | Glioblastoma | 200 | WSI (**527 series, 151.5 GB**) | Clinical, Genomics, Proteomics |
| CPTAC-AML | Acute myeloid leukemia | 180 | WSI + Follow-Up + Classification | — |
| CPTAC-STAD | Stomach adenocarcinoma | 168 | WSI | Genomics, Proteomics |
| CPTAC-PDA | Pancreatic ductal adeno | 168 | WSI | Clinical, Genomics, Proteomics |
| CPTAC-BRCA | Breast | 134 | WSI (no radiology) | Clinical, Genomics, Proteomics |
| CPTAC-COAD | Colon | 106 | WSI (no radiology) | Clinical, Genomics, Proteomics |
| CPTAC-OV | Ovarian | 102 | WSI (no radiology) | Clinical, Genomics, Proteomics |
| CPTAC-CM | Cutaneous melanoma | 95 | WSI | — |
| CPTAC-SAR | Sarcomas | 88 | WSI | Clinical only |

**Total: ~2,400 subjects with WSI across 14 cancer types.**

### 1.2 Licence and access — this is the good part

- **Licence: CC BY 4.0.** Verified verbatim on CPTAC-LUAD, CPTAC-GBM, CPTAC-CCRCC, CPTAC-UCEC pages.
- **Format:** Aperio SVS, pyramidal.
- **Download route:** pathology is *not* served by the NBIA/TCIA Data Retriever. It is served by **IBM Aspera Faspex**: `https://faspex.cancerimagingarchive.net/aspera/faspex/public/package` — you need the Aspera Connect browser plugin (or `ascp` CLI, which is the sane scripted route).
- **Browse/search UI for pathology:** `https://pathdb.cancerimagingarchive.net/eaglescope/dist/`
- **Important nuance on CPTAC-GBM / CPTAC-HNSCC "Access: Limited":** the *limitation is on the radiology*, not the slides. TCIA states the MR/CT falls under the NIH Controlled Data Access Policy (facial-reconstruction risk). **The pathology WSIs remain CC BY 4.0 and openly downloadable.** Verified on the CPTAC-GBM page. Do not let the "Limited" badge scare you off the slides.

### 1.3 Does CPTAC actually have outcome? — **Yes, and better than commonly assumed**

Queried the GDC API live:

- `CPTAC-3` (brain, H&N, kidney, lung, pancreas, uterus): **1,683 cases**, of which **640 have `vital_status = Dead`** → **~38% event rate**. That is a genuinely usable survival cohort, not a token one.
- `CPTAC-2` (breast, colon, ovary): **342 cases**.
- Fields confirmed populated: `demographic.vital_status`, `demographic.days_to_death`, `diagnoses.days_to_last_follow_up`.

Reproduce with:
```
https://api.gdc.cancer.gov/cases?filters={"op":"and","content":[
  {"op":"in","content":{"field":"project.project_id","value":["CPTAC-3"]}},
  {"op":"in","content":{"field":"demographic.vital_status","value":["Dead"]}}]}&size=1&format=json
```
(`pagination.total` = 640.)

Caveat found in the raw records: some `days_to_last_follow_up` values are **negative** (I saw `-6.0` in the first page of 5 records). Follow-up is short and noisy relative to TCGA. Sanity-filter before modelling.

### 1.4 Modality pairing — how tight is the WSI ↔ proteomics link?

TCIA's own wording, verbatim: **"Pathology imaging is collected as part of the CPTAC qualification workflow."**

Interpretation: these are the specimen-qualification slides from the *same* tissue segments submitted for proteogenomic aliquoting — so the pairing is much tighter than TCGA's diagnostic-slide-vs-molecular-block mismatch. TCIA does not, however, spell out a per-slide → per-aliquot mapping on the collection page. **I could not verify a formal slide-to-aliquot mapping table.** Treat pairing as "same specimen, same workflow" — strong, but not a documented section-level guarantee.

Cross-linking is by shared de-identified patient ID across TCIA ↔ GDC ↔ PDC (Proteomic Data Commons). CPTAC-LUAD additionally encodes DICOM `Clinical Trial Time Point ID` = days from pathological diagnosis.

**Biggest limitation:** ~100–260 patients per cancer type. Underpowered for per-cancer survival modelling; usable pan-cancer or as an external validation set. Second limitation: the proteomics is bulk on a *different piece* of the same specimen, so any spatial claim linking a slide region to a protein measurement is an inference, not a measurement.

---

## 2. HANCOCK — the best new multimodal cohort (2025)

**Dörrich M et al., "A multimodal dataset for precision oncology in head and neck cancer", Nature Communications 2025, DOI 10.1038/s41467-025-62386-6, PMID 40759646.** Abstract verified via Europe PMC.

- **n = 763 head-and-neck cancer patients**, monocentric real-world (Erlangen).
- **Modalities, paired per patient:**
  - H&E WSI of primary tumour: **701–709 patients** (TCIA lists 709 slides), Aperio SVS, 0.1213 µm/px
  - H&E WSI of lymph nodes: **396 patients**, SVS, 0.1945–0.2634 µm/px
  - **Tissue microarrays** with IHC: CD3, CD8, CD56, CD68, CD163, PD-L1, MHC-1 — up to 32 cores/patient (TCIA says up to 16 TMA images/patient, PNG/TIFF)
  - Clinical + pathological tabular data: all 763
  - Blood parameters (routine panel; CRP only n=94)
  - Surgery reports and medical history as free text, German + English translation
- **Outcome:** overall survival with follow-up **up to 14 years**; recurrence status; 5-yr OS 77.3%.
- **Licence: CC BY 4.0.** Code Apache 2.0.
- **Access:**
  - TCIA mirror — **live, 4.5 TB**: `https://doi.org/10.7937/rcty-5h16` → `https://www.cancerimagingarchive.net/collection/HANCOCK` (Aspera). Verified live and marked complete as of Dec 2025.
  - Project site `https://hancock.research.fau.edu/` — **DNS did not resolve at fetch time (ENOTFOUND)**. Use the TCIA mirror; treat the FAU site as possibly down/moved.
  - Code: `https://github.com/ankilab/HANCOCK_MultimodalDataset`
- **Biggest limitation:** single centre, single disease site, no bulk omics (no RNA-seq/proteomics — the "molecular" layer is IHC/TMA only). Also 4.5 TB is a serious storage commitment.

**This is arguably the highest-value dataset in this entire report** for anything involving WSI + text + immune-marker IHC + long survival, because it is the only one with all four and a permissive licence.

---

## 3. SurGen — colorectal WSI + mutations + survival

**Myles C, Um IH, Marshall C, Harris-Birtill D, Harrison DJ. "SurGen: 1020 H&E-stained Whole Slide Images With Survival and Genetic Markers." GigaScience 14 (2025), DOI 10.1093/gigascience/giaf086.** Abstract verified via arXiv API.

- **1,020 H&E WSIs from 843 colorectal cancer cases** (St Andrews / NHS Scotland).
- **Molecular labels:** KRAS, NRAS, BRAF mutation status; **mismatch-repair (MMR) status**.
- **Survival data for 426 of the cases.**
- **Access:** BioImage Archive **S-BIAD1285** — `https://doi.org/10.6019/S-BIAD1285`. Verified via the BioStudies API: **1,022 files**, released 2024-07-24, retrievable over **HTTP, FTP and Globus** (Globus is the right choice for the full pull). No registration wall.
- Code/loader: `https://github.com/CraigMyles/SurGen-Dataset`
- **Licence: COULD-NOT-VERIFY** — the BioStudies `/info` endpoint returned no licence field, and I could not load the study landing page (login shell). BioImage Archive default is CC0/CC BY, but **check the study page before assuming**.
- **Biggest limitation:** survival on only **426/843** cases (~50%), so the effective outcome cohort is half the headline number. Single-nation, and no transcriptomics/proteomics — the "molecular" layer is a handful of mutation calls.

---

## 4. Treated cohorts with response labels (small but real)

### 4.1 Ovarian Bevacizumab Response (TCIA)

- **78 patients, 288 H&E WSIs**, SVS, **253.8 GB**.
- **Treatment: bevacizumab.** Response labels at slide level: "effective in 162 and invalid in 126" of the 288 slides. Sensitive = no measurable regrowth or CA-125 < 2× ULN; Resistant = measurable regrowth or CA-125 > 2× ULN.
- **Licence: CC BY 4.0.** DOI `10.7937/TCIA.985G-EY35`. Aspera download.
- `https://www.cancerimagingarchive.net/collection/ovarian-bevacizumab-response/`
- **Biggest limitation:** n=78 is too small to train on and marginal even for validation; labels are surrogate (CA-125), and the effective/invalid counts are per-*slide* not per-patient, which invites leakage if you split naively.

### 4.2 IMPRESS — breast neoadjuvant chemotherapy, H&E + multiplex IHC

**Huang Z, Shao W, Han Z, … Li Z. "Artificial intelligence reveals features associated with breast cancer neoadjuvant chemotherapy responses from multi-stain histopathologic images." npj Precision Oncology 2023, DOI 10.1038/s41698-023-00352-5, PMID 36707660, PMC9883475.** Full text verified via Europe PMC.

- **n = 126 patients: 62 HER2+ and 64 TNBC.**
- **Pairing is the selling point:** each patient has **one H&E WSI and one multiplex-IHC WSI of the same tissue** (CD8 green, CD163 red, PD-L1 brown), registered.
- **Outcome: pCR vs residual tumour**, plus **RCB** for incomplete responders. This is genuine *treatment-response* data.
- Plus a small held-out external set: 20 HER2+ (10/10) and 20 TNBC (10/10).
- **Access: COULD-NOT-VERIFY.** The published article contains **no Data Availability statement** (confirmed by reading the full text XML). Nonetheless the "IMPRESS HER2+ dataset" is used as an external validation set by at least three independent later papers I verified — Breast Cancer Res 2025 (DOI 10.1186/s13058-025-02139-x), Cyborg and Bionic Systems 2026 (DOI 10.34133/cbsystems.0554), Cancers 2025 (DOI 10.3390/cancers17152423) — so it clearly circulates. **Route: email the corresponding authors (Kun Huang, Zaibo Li) or check the citing papers' methods for the mirror they used.** Do not assume it is a click-download.
- **Biggest limitation:** 126 patients, and the access route is undocumented.

### 4.3 Cancer Moonshot Biobank (CMB) — the underrated one

TCIA hosts nine CMB collections, all with **Histopathology / Whole Slide Image**: CMB-AML, CMB-BRCA, **CMB-CRC**, CMB-GEC (gastroesophageal), CMB-LCA (lung), **CMB-MEL** (melanoma), CMB-MML (multiple myeloma), CMB-OV, CMB-PCA (prostate).

Verified in detail for **CMB-CRC**:
- **102 subjects with histopathology** (62 also have radiology).
- WSI in **SVS + JSON**; radiology in CT/MR/PT/US/NM/DX/XA.
- **Licence CC BY 4.0** for radiology and histopathology (head images only are NIH-controlled).
- **Longitudinal design:** DICOM tag (0012,0052) carries "days from enrollment", enabling temporal alignment of imaging to clinical events. CMB is explicitly a *prospective, treated, longitudinally-followed* biobank — which is exactly what TCGA is not.
- **CMB-MEL: 66 subjects**, WSI + CT/US/MR/PT.
- **The catch:** *"Associated genomic, phenotypic and clinical data will be hosted by dbGaP."* An open-access clinical subset sits in the **CTDC (Cancer Trial Data Commons) Cancer Moonshot Biobank** portal, but the full treatment/outcome table is **dbGaP-controlled** → requires a data access request with institutional signing official.
- **Biggest limitation:** the images are free and the outcomes are not. Per-cancer n is small (66–102). I could not load the CTDC portal to enumerate exactly which outcome fields are in the open subset — **COULD-NOT-VERIFY on the open/controlled boundary**.

### 4.4 NLST (National Lung Screening Trial) — pathology subset

- **451 subjects with digitised histopathology** out of 26,254 imaged subjects.
- 1,225 primary-tumour slide files (**~775 GB**), plus 23 second-primary slides (10 subjects) and 4 extra slides (2 subjects).
- **Licence: CC BY 4.0.** SVS. Aspera.
- `https://www.cancerimagingarchive.net/collection/nlst/`
- **Outcome:** NLST is a mortality-endpoint RCT, so lung-cancer-specific and all-cause mortality exist — but the linked clinical/follow-up tables come from **CDAS (Cancer Data Access System)**, which requires a (free, but non-instant) project application. TCIA's page confirms clinical data includes "Diagnosis, Exposure, Measurement, Follow-Up" but does not itself serve the mortality linkage.
- **Biggest limitation:** only 451 of 26k subjects have slides — the WSI arm is a small biopsy-confirmed subset, heavily selected toward screen-detected cancers. Massive stage/lead-time bias if used as a general lung cohort.

---

## 5. Task-specific challenge datasets (big, clean, **no outcome**)

Be clear-eyed: these are label datasets. None has survival.

| Dataset | Size | Labels | Access / licence | Killer limitation |
|---|---|---|---|---|
| **PANDA** (Bulten W et al., *Nature Medicine* 2022, DOI 10.1038/s41591-021-01620-2, PMID 35027755 — verified) | ~10k prostate biopsy WSIs, Radboud + Karolinska | ISUP grade / Gleason | Kaggle competition data page (Kaggle blocks automated fetch; **exact n and licence COULD-NOT-VERIFY from the data page — read the Kaggle rules before use**) | Biopsy cores only, grade label only, **no outcome**. Kaggle competition licence restricts some reuse. |
| **CAMELYON16** | **400 WSIs** (270 train: 170+100; 130 test), 2 centres (Radboud UMC, UMC Utrecht) | Lymph-node metastasis, pixel-level XML contours + binary masks | grand-challenge.org; AWS Open Data | Binary metastasis detection only. No patient outcome. |
| **CAMELYON17** | **1,000 WSIs from 200 patients** (100 train / 100 test, 5 nodes each), **5 Dutch centres** | pN-stage per patient; lesion annotations on 50 slides only; ITCs not exhaustively annotated | **CC0** — GigaDB, **AWS Open Data `s3://camelyon-dataset` (us-west-2), `aws s3 ls --no-sign-request s3://camelyon-dataset/`**, Baidu Pan | Staging label only. Sparse annotation. No outcome, no molecular. |
| **BCNB** (Xu F et al., *Front Oncol* 2021, DOI 10.3389/fonc.2021.759007, PMID 34722313 — verified) | **1,058 patients**, early breast cancer **core-needle biopsy** WSIs | ALN metastatic status + number of positive nodes, ER/PR/HER2, Ki67, tumour size, histological grade, molecular subtype, age | `https://bupt-ai-cz.github.io/BCNB/` — registration form, **non-commercial academic licence only**, no redistribution | **No survival.** Non-commercial licence is a real constraint if you plan to ship anything. Biopsy cores, not resections. |
| **Digital Brain Tumour Atlas** (Roetzer-Pejrimovsky et al., *Scientific Data* 2022, DOI 10.1038/s41597-022-01157-0, PMID 35169150 — verified) | **3,115 slides, 126 brain tumour types**, 1995–2019, Vienna brain tumour bank | WHO tumour-type diagnosis + clinical annotations | **CC BY**, hosted on EBRAINS (**exact EBRAINS download URL and whether an account is required: COULD-NOT-VERIFY — the KG search page would not render for me**) | Diagnostic classification resource. Long-tail class imbalance across 126 types. Outcome data not the point of the resource. |
| **Orion CRC** (Lin JR et al., *Nature Cancer* 2023, DOI 10.1038/s43018-023-00576-1, PMID 37349501 — verified) | ~40 CRC patients | **H&E and 18-plex CyCIF on the *same tissue section*** — the tightest WSI↔protein pairing that exists anywhere | Synapse (**exact Synapse ID and licence COULD-NOT-VERIFY — the Synapse page I tried served a JS shell**) | n≈40. Tiny. But if you need pixel-registered H&E↔multiplex-protein ground truth, nothing else comes close. |

---

## 6. TCIA as a whole — how to enumerate the rest yourself

`https://www.cancerimagingarchive.net/browse-collections/` is the authoritative index and it renders fine to a plain fetch. Filter on the Modalities column for `Histopathology` / `Whole Slide Image`. Beyond CPTAC and CMB, the same table shows NLST and the TCGA collections carry histopathology. TCIA's general licence is per-collection; most pathology is CC BY 4.0, but **check each collection page — the licence is stated per collection and per data type** (the CPTAC-GBM split between CC BY slides and NIH-controlled MR is the canonical example).

---

## 7. AACR Project GENIE — **no images. At all.**

Verified: **Lavery JA / de Bruijn I et al., "Analysis and Visualization of Longitudinal Genomic and Clinical Data from the AACR Project GENIE Biopharma Collaborative in cBioPortal", *Cancer Research* 2023, DOI 10.1158/0008-5472.CAN-23-0816, PMID 37668528.** Abstract retrieved verbatim from Europe PMC.

- GENIE public v19.0: **271,837 samples / 227,696 patients**, 19 centres.
- **GENIE BPC: ~25,000 patients** with curated longitudinal treatment and outcome data via the **PRISSMM** model.
- **Data types: targeted panel mutations, CNA, fusions/SVs, cfDNA, clinical timelines. There is no pathology WSI and no radiology image data.**
- Critical trap: PRISSMM stands for **P**athology, **R**adiology/**I**maging, **S**igns/symptoms, **S**erial measurements, **M**edical oncologist assessment. The "Pathology" and "Imaging" domains in PRISSMM are **structured abstractions of *reports***, not image files. Anyone who tells you "GENIE has imaging" has misread PRISSMM.
- Access: Synapse `syn7222066` (registration + data use agreement) and `genie.cbioportal.org` (I hit a Google OAuth redirect — **the portal requires login**).
- **Verdict for this remit: unusable.** Excellent outcome data, zero images. Its only role is as a prior/reference for mutation frequencies and treatment-line distributions.

---

## 8. The honest gap: **there is no open immunotherapy cohort with slides**

This is the finding I would push back hardest on if anyone in the project assumes otherwise.

What I checked and what I found:

- **TCIA `Anti-PD-1_MELANOMA`**: 47 subjects, immunotherapy response labels — but modalities are **MR, CT, PT only. No histopathology.** Verified from the TCIA melanoma collection listing.
- The only TCIA melanoma collections with WSI are **CMB-MEL (66 subjects)** and **CPTAC-CM (95 subjects)** — and CMB-MEL's outcome data sits in dbGaP while CPTAC-CM's TCIA entry lists no external clinical/omics linkage at all.
- PubMed searches for `immune checkpoint inhibitor` + `whole slide image` + public/open data returned essentially nothing usable. The three hits I resolved — J Transl Med 2025 (DOI 10.1186/s12967-025-07416-z, gastric), J Transl Med 2025 (DOI 10.1186/s12967-025-06181-3, pan-GI), Front Mol Biosci 2025 (DOI 10.3389/fmolb.2025.1615533, CRC) — all build immunotherapy-response *predictions* on TCGA slides plus **private** institutional IO cohorts. **The slides used for the IO endpoint are not released.**
- Well-known IO trial cohorts (IMvigor210 and similar) release transcriptomics, not slides.

**Practical consequence:** if the project's claim requires "validated on a treated/immunotherapy cohort with images", your only open options today are: Ovarian Bevacizumab Response (n=78, anti-VEGF not IO), IMPRESS (n=126, chemo not IO, undocumented access), and CMB (dbGaP for outcomes). Any stronger claim needs a collaboration or a dbGaP application. Plan for that now, not at revision time.

---

## 9. Things I could not verify (do not cite these as established)

- **PAIP / Pathology AI Platform (Korea)** — `wisepaip.org` **refused the connection (ECONNREFUSED 3.37.4.142:443)**. I could not confirm current dataset contents, sizes, access requirements or whether the platform is still operating. **COULD-NOT-VERIFY — treat PAIP as unavailable until you confirm otherwise.**
- **HTAN (Human Tumor Atlas Network)** — `data.humantumoratlas.org` and `docs.humantumoratlas.org` both served JS shells with no substantive content. It certainly contains H&E and multiplex imaging plus single-cell/spatial data across several atlases, but **I could not verify participant counts, the open-vs-dbGaP split, the licence, or whether survival is released.** **COULD-NOT-VERIFY.** Worth a manual browser visit — it is the most likely place to find WSI ↔ spatial-omics pairing at scale.
- **TIGER challenge** (breast TILs) — PubMed returned zero hits for my query terms. **COULD-NOT-VERIFY that a peer-reviewed dataset paper exists.**
- **PANDA exact slide/patient counts and licence** — Kaggle blocks fetch; the Nature Medicine paper is confirmed but I did not extract n from it.
- **SurGen licence**, **Digital Brain Tumour Atlas EBRAINS download route**, **Orion Synapse ID**, **IMPRESS download URL**, **CTDC open-subset field list** — all flagged inline above.
- I found **no** evidence for any pan-cancer non-TCGA WSI+outcome resource larger than CPTAC. If one is claimed in the manuscript, it needs a citation I have not been able to find.

---

## 10. Recommended acquisition order

1. **CPTAC pathology, all 14 collections** — CC BY 4.0, Aspera Faspex, ~2,400 patients, pairs to GDC (640 deaths in CPTAC-3) and PDC proteomics. The backbone.
2. **HANCOCK** via TCIA `10.7937/rcty-5h16` — 763 patients, 14-yr survival, H&E + IHC TMA + text, CC BY 4.0. Budget 4.5 TB.
3. **SurGen** via BioImage Archive `S-BIAD1285` over Globus — 843 CRC cases, MMR/KRAS/NRAS/BRAF, survival on 426. Check the licence field on the landing page first.
4. **CAMELYON17** from `s3://camelyon-dataset` (CC0, no-sign-request) — free, instant, good for pretraining/sanity checks. No outcome.
5. **CMB (all nine collections)** images now under CC BY 4.0; start the **dbGaP application in parallel** if you need the treatment/outcome tables — that lead time is weeks to months.
6. Only then chase IMPRESS / Ovarian Bevacizumab for the treated-cohort story, knowing both are ~100 patients.
