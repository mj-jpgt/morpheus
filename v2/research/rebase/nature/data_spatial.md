# Open-Access Data Scout — Spatial Transcriptomics Paired with H&E

**Scout date:** 2026-07-29
**Remit:** Spatial transcriptomics datasets with paired H&E histology. Exact sizes, cancers covered, resolution, licence, download route.

**Verification method (important):** WebSearch was exhausted for this session. Every claim below was verified by direct WebFetch against arXiv API, PubMed E-utilities, HuggingFace dataset pages, GitHub repos, or vendor dataset pages. Anything I could not confirm from a primary source is explicitly marked **COULD-NOT-VERIFY**. No citation here was written from memory.

---

## Quick comparison table

| Dataset | n slides/samples | Cancer coverage | Platform(s) | Licence | Gated? |
|---|---|---|---|---|---|
| **HEST-1k** | 1,229 (paper) / 1,276 (current HF) | 367 cancer samples, 25 subtypes | ST, Visium, Visium HD, Xenium | CC BY-NC-SA 4.0 | Yes (HF token + accept terms) |
| **HEST-Benchmark** | 259 rows (187 train / 72 test), 42.2 GB | 9 tasks, 9 cancer types | subset of HEST-1k | CC BY-NC-SA 4.0 | Yes |
| **STimage-1K4M** | 1,149 slides / 4,293,195 spot-image pairs | 456 slides (39.7%) cancer-related | ST (151), Visium (994), Visium HD (4) | MIT (repo card) | No |
| **10x public datasets** | catalogue, per-dataset download | many single tumours (CRC, breast, RCC, lung…) | Visium, Visium HD, Xenium | CC BY 4.0 | Registration form typical |
| **HESCAPE** | 54 donors, 6 gene panels | breast, lung, colon, IO panels | 10x Xenium | MIT | No |
| **SpaRED** | 26 curated public datasets | mixed | Visium-family | MIT | No (Google Drive) |
| **HER2ST (Andersson)** | 36 sections | HER2+ breast cancer | legacy ST (100 µm) | not stated | Password-protected 7z |
| **Mendeley 29ntw7sh4r** | 23 patients, triplicate sections | breast cancer (LumA/LumB/TNBC/HER2+) | legacy ST | CC BY 4.0 | No |
| **HumanST-1k** | 1.8M H&E–transcriptome pairs | pan-organ, oncology-focused | multi-platform | **COULD-NOT-VERIFY** | **COULD-NOT-VERIFY** |

---

## 1. HEST-1k — the anchor dataset

**Exact name:** HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis
**Citation (verified via arXiv API):** Jaume G, Doucet P, Song AH, Lu MY, Almagro-Pérez C, Wagner SJ, Vaidya AJ, Chen RJ, Williamson DFK, Kim A, Mahmood F. arXiv:2406.16192v2 (submitted 2024-06-23, updated 2024-11-02). NeurIPS 2024 Spotlight.

### Size (exact, as published)
Verbatim from the abstract: *"HEST-1k, a collection of 1,229 spatial transcriptomic profiles, each linked to a WSI and extensive metadata. HEST-1k was assembled from 153 public and internal cohorts encompassing 26 organs, two species (Homo Sapiens and Mus Musculus), and 367 cancer samples from 25 cancer types."*

Derived assets stated in the paper: **2.1 million expression–morphology pairs** and **76 million detected nuclei**.

### Size (current, as of this scout)
The live HuggingFace card reports **1,276 samples** and **2.01 TB** total file size, with **>1.5 million expression/morphology pairs and >76 million nuclei**. The card notes incremental additions: **27 Visium HD samples added January 2026** and **18 Xenium samples added February 2026** (transcripts + segmentation now public). So the repo has grown past the paper.

> **COULD-NOT-VERIFY:** a clean per-technology sample breakdown. The arXiv HTML does not contain a single consolidated technology table — it is spread across appendix tables A2–A10. Directional statement that *is* supported: Visium dominates, with substantial Xenium, and smaller counts of legacy ST and Visium HD. Do **not** quote a numeric split without opening the HF metadata yourself (`HEST_v1_1_0.csv` in the repo root carries per-sample technology fields).

### Modalities and how they are PAIRED
This is the key property and the reason HEST-1k is the anchor:

- Each sample = **one WSI (H&E) + one spatially-resolved transcriptome**, physically the same tissue section, registered in the same pixel coordinate frame.
- Pairing is at the **spot level**: for each spot, HEST extracts a **224×224 px patch at 20× magnification**, which the paper states is equivalent to a **112×112 µm H&E region**.
- Quality gate (verbatim): *"All images with a pixel size higher than 1.15 μm/px were discarded to ensure an acceptable image quality."*
- Xenium harmonisation (verbatim): *"we generated 'pseudo-Visium' spots by pooling transcripts on 55×55-μm patches without spacing."* — i.e. Xenium subcellular data is downsampled to Visium-like spots so all technologies share one interface. If you need true subcellular Xenium resolution, HEST is **not** the right entry point; go to the raw 10x/GEO release.
- Extra layers shipped: nuclei segmentation (76M nuclei), tissue segmentation masks, and rich metadata (organ, oncotree code, species, cohort).

### Cancers
367 cancer samples across 25 cancer types (26 organs total including non-cancer). The benchmark tasks name the well-represented ones: breast IDC, prostate PRAD, pancreas PAAD, skin SKCM, colon COAD, rectum READ, kidney ccRCC, lung LUAD, and axillary lymph node IDC metastases.

### Licence and download route
**Licence:** CC BY-NC-SA 4.0 (Attribution-NonCommercial-ShareAlike). Verbatim from the paper: *"HEST-1k, HEST-Benchmark, and HEST library are released under the Attribution-NonCommercial-ShareAlike 4.0 International license (CC BY-NC-SA 4.0)."*

**NonCommercial + ShareAlike is a real constraint** — it forecloses commercial licensing of anything derived, and ShareAlike propagates to derivative datasets. If the plan is a commercial product or a permissively licensed model release, HEST-1k contaminates it.

**Download:**
- Dataset: https://huggingface.co/datasets/MahmoodLab/hest — **gated**. You must create a HuggingFace account, accept the terms of use on the card, and generate an auth token. (Confirmed empirically: fetching `.../raw/main/README.md` returns **HTTP 401**, so the gate is live.)
- Library / code: https://github.com/mahmoodlab/HEST
  ```bash
  git clone https://github.com/mahmoodlab/HEST.git
  cd HEST
  conda create -n hest python=3.11 && conda activate hest
  pip install -e .
  ```
- Subsetting: `huggingface_hub.snapshot_download()` with `allow_patterns` filtered by sample ID, organ, species, or oncotree code. **Do this** — the full pull is 2 TB.

### Biggest limitation
**CC BY-NC-SA 4.0.** Non-commercial and viral. Secondary limitation: the Xenium samples are pseudo-Visium-binned inside HEST, so the dataset's effective resolution ceiling is ~55 µm spots / 112 µm patches — it cannot support single-cell-resolution claims without going back to source.

---

## 2. HEST-Benchmark — the evaluation harness

**Exact name:** HEST-Benchmark (a.k.a. HEST-bench)
**Composition:** 9 gene-expression-prediction tasks, each predicting **50 genes** from an H&E patch. Table extracted from arXiv:2406.16192v2:

| # | Task | Organ | Patients | Samples | Genes |
|---|---|---|---|---|---|
| 1 | IDC | Breast | 4 | 4 | 50 |
| 2 | PRAD | Prostate | 2 | 23 | 50 |
| 3 | PAAD | Pancreas | 3 | 3 | 50 |
| 4 | SKCM | Skin | 2 | 2 | 50 |
| 5 | COAD | Colon | 2 | 4 | 50 |
| 6 | READ | Rectum | 2 | 4 | 50 |
| 7 | ccRCC | Kidney | 24 | 24 | 50 |
| 8 | LUAD | Lung | 2 | 2 | 50 |
| 9 | IDC | Axillary lymph node | 4 | 4 | 50 |

**Current HF artifact:** https://huggingface.co/datasets/MahmoodLab/hest-bench — **259 rows** (187 train / 72 test), **42.2 GB**, sample IDs `INT1`–`INT24`, patches as `.h5`, expression as `.h5ad`. Image + text modalities. From "AI for Pathology Image Analysis Lab @ HMS/BWH".

**Drift note:** the GitHub README lists the tasks as IDC, PRAD, PAAD, SKCM, COAD, READ, CCRCC, LUNG, LYMPH_IDC — and the HF card additionally mentions *"a new colorectal adenocarcinoma task based on 4 Xenium samples"* in the updated leaderboard. So the live benchmark has ≥10 tasks. **Treat the 9-task table above as the paper-of-record version and re-derive counts from the current repo before quoting.**

**Biggest limitation:** brutally small patient counts. Seven of nine tasks have **2–4 patients**. PRAD has 23 samples from **2 patients**; ccRCC (24/24) is the only task with a defensible patient-level n. Any HEST-Bench improvement is a claim about a handful of donors, and patient-level cross-validation is barely meaningful. This is the single most attackable point in any paper that leans on HEST-Bench as primary evidence.

---

## 3. STimage-1K4M — the largest spot-level pairing

**Exact name:** STimage-1K4M: A histopathology image-gene expression dataset for spatial transcriptomics
**Citation (verified via arXiv API):** Chen J, Zhou M, Wu W, Zhang J, Li Y, Li D. arXiv:2406.06393v2 (2024-06-10).

### Size (exact, from arXiv HTML)
- **1,149 slides** total.
- **4,293,195 spot–image pairs** — the headline number and the reason to care.
- Slide-level technology split: **legacy ST 151 slides (13.1%)**, **Visium 994 slides (86.5%)**, **Visium HD 4 slides (0.3%)**.
- Spot-level split: ST 1.4%, Visium 54.4%, **Visium HD 44.2%** — i.e. 4 Visium HD slides supply nearly half of all spots. Effective slide diversity is far lower than the spot count implies.
- **50 distinct tissue types**; brain 251 slides (21.8%), breast 205 slides (17.8%).
- **456 slides (39.7%) cancer-related.**
- **10 species**, predominantly human and mouse.

### Resolution (exact, from the paper's technology table)
| Technology | Spot diameter | Centre-to-centre | Genes measured |
|---|---|---|---|
| Legacy ST | 100 µm | 200 µm | ~15k–30k |
| Visium | 55 µm | 100 µm | ~15k–30k |
| Visium HD | — | 8 µm × 8 µm bins | ~15k–30k |

### Pairing
Each slide is tiled into sub-tiles centred on spatial spots; **each tile is paired with a 15,000–30,000-dimensional gene expression vector** (full transcriptome, not a panel). Spot coordinates and per-slide spot radius are shipped alongside so you can re-tile at your own patch size. This is the distinguishing feature versus HEST: full transcriptome per spot rather than a 50-gene benchmark target.

### Licence and download route
**Licence:** **MIT** on both the GitHub repo (*"All code is licensed under the MIT License"*) and the HuggingFace dataset card. This is the permissive alternative to HEST-1k.

> **Caveat you must check yourself:** MIT is stated for the *repo/card*. The underlying slides are aggregated from GEO accessions (confirmed on the HF card: GSE144239, GSE148612, GSE153424, GSE175540, GSE197023, GSE210616, GSE211895, GSE212323, among others) and from 10x releases. Upstream per-cohort terms are not necessarily MIT and the aggregators cannot relicense them. For a commercial release, audit source-by-source.

**Download:**
- https://huggingface.co/datasets/jiawennnn/STimage-1K4M — **ungated, public**.
- Top-level structure confirmed: folders `ST/`, `Visium/`, `VisiumHD/`, `annotation/`, `aux/`, `meta/`, plus `README.md`.
- Code: https://github.com/JiawenChenn/STimage-1K4M
- Docs site https://jiawenchenn.github.io/STimage-1K4M/ exists but its `/docs/01-download` path returned **404** at scout time — use the HF repo tree directly.
- **COULD-NOT-VERIFY:** total repository size in GB/TB. HF reports only the size category `100B < n < 1T`.

### Biggest limitation
**The HuggingFace dataset viewer is currently broken** — schema-mismatch error: 17,138 gene columns present but inconsistent across files, and `xaxis`/`yaxis`/`r` coordinate columns missing from some files. You cannot rely on `load_dataset()` streaming; you must `snapshot_download` the folders and parse per-slide yourself, handling heterogeneous gene panels. Secondary limitation: only 39.7% is cancer, and brain (21.8%) is the largest tissue — the composition is skewed toward neuroscience, not oncology.

---

## 4. 10x Genomics public releases — Visium / Visium HD / Xenium

**Route:** https://www.10xgenomics.com/datasets

**Licence — verified on individual dataset pages:** *"This dataset is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license."* Confirmed independently on both a Visium HD page and a Xenium page. **CC BY 4.0 is the most permissive licence of anything in this scout** — commercial use allowed with attribution.

**Verified example datasets (each carries H&E):**

1. **Visium HD Spatial Gene Expression Library, Human Colorectal Cancer (FFPE)** — sigmoid colon, male, age 60, 5 µm FFPE section, Space Ranger 3.0.0, H&E stained and imaged, 11 mm capture area, 2 µm bins (binned to 8 µm/16 µm). URL path `/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc`. Also surfaced as `/datasets/visium-hd-cytassist-11mm-human-colon-cancer-HE`.
2. **FFPE Human Breast with Custom Add-on Panel (Xenium)** — Xenium Human Breast Gene Expression v1 panel **plus 100 custom genes**; two tissue samples with **574,527** and **365,604** cells detected; **registered post-Xenium H&E image provided**. The page carries an explicit warning that the H&E *"may or may not be suitable for more refined exact sub-cellular resolution correspondence due to minute differences in the optics of the microscopy systems"* — read that as: registration is good to ~cell-neighbourhood, not reliably to sub-cellular.
3. **Xenium In Situ Gene and Protein Expression, Human Renal Cell Carcinoma (FFPE)** — path `/datasets/xenium-protein-ffpe-human-renal-carcinoma`.

**Associated peer-reviewed anchor (verified via PubMed):** Janesick A, …, Taylor SEB. *High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis.* Nat Commun. 2023;14(1):8353. DOI 10.1038/s41467-023-43458-x. PMID 38114474. This is the breast-cancer study behind the widely used matched scRNA-seq + Visium + Xenium FFPE breast sections.

**Biggest limitation:** **n = 1 per dataset.** These are technology demonstrations, not cohorts. There is no patient metadata, no outcome data, no consistent processing, and no programmatic index — the catalogue page is JavaScript-rendered and returned "Showing 0 datasets" to a plain fetch, so you cannot enumerate it via API. Downloads generally require submitting an email/registration form even though the licence is open.

> **COULD-NOT-VERIFY:** the complete count of public 10x spatial datasets, and whether registration is mandatory for every file (the CRC page did not state a registration requirement, but the catalogue is JS-gated so this could not be confirmed at scale).

---

## 5. HESCAPE — pan-organ Xenium benchmark (2025)

**Exact name:** HESCAPE — *A Large-Scale Benchmark of Cross-Modal Learning for Histology and Gene Expression in Spatial Transcriptomics*
**Citation (verified via arXiv API):** Gindra RH, Palla G, Nguyen M, Wagner SJ, Tran M, Theis FJ, Saur D, Crawford L, Peng T. arXiv:2508.01490 (2025-08-02, updated 2025-08-27).

- **Size:** curated pan-organ dataset spanning **6 gene panels and 54 donors**; 10x **Xenium**. Panels cover human breast, lung, colon, immuno-oncology, plus multi-tissue combinations.
- **Pairing:** H&E image patches ↔ Xenium panel expression; built for cross-modal *contrastive* pretraining (image↔gene retrieval, Recall@5), plus two downstream tasks: **gene mutation classification** and gene expression prediction.
- **Licence:** MIT (repo).
- **Download:** code https://github.com/peng-lab/hescape ; data https://huggingface.co/datasets/Peng-AI/hescape-pyarrow — ungated, loadable directly:
  ```python
  ds = load_dataset("Peng-AI/hescape-pyarrow", name="human-breast-panel", split="train")
  ```
  Pretrained gene-model weights at https://huggingface.co/Peng-AI/hescape-drvi. Vision encoder weights must be obtained separately from their own repos.
- **Why it matters strategically:** this is the only source here that reports a *negative* result you can cite defensively — verbatim: *"while contrastive pretraining consistently improves gene mutation classification performance, it degrades direct gene expression prediction compared to baseline encoders trained without cross-modal objectives."* They identify **batch effects** as the interfering factor. If your thesis involves H&E↔ST contrastive alignment, this paper is either your strongest prior-art threat or your strongest motivation, depending on framing.
- **Biggest limitation:** Xenium panels only — a few hundred targeted genes, not the transcriptome. Cross-panel generalisation is exactly what it shows is hard. 54 donors is small for "pan-organ".

---

## 6. SpaRED / SpaCKLE — curated benchmark of 26 public ST datasets

**Citations (both verified via arXiv API — note these are two related papers, cite the right one):**
- Ruiz D, Cárdenas P, Manrique L, Vega D, Mejia GM, Arbeláez P. *Completing Spatial Transcriptomics Data for Gene Expression Prediction Benchmarking.* arXiv:2505.02980 (2025-05-05, updated 2025-09-04).
- Mejia G, Ruiz D, Cárdenas P, Manrique L, Vega D, Arbeláez P. *SpaRED benchmark: Enhancing Gene Expression Prediction from Histology Images with Spatial Transcriptomics Completion.* arXiv:2407.13027 (2024-07-17, updated 2024-09-27).

- **Size:** **26 public datasets** systematically curated and reprocessed — described as an **8.6× increase over previous works**. Visium-family, histology-paired.
- **Contribution beyond curation:** **SpaCKLE**, a transformer-based gene-expression *completion* model addressing dropout; reported **>82.5% MSE reduction** vs prior approaches. The benchmark evaluates **8 prediction models** on both raw and SpaCKLE-completed data.
- **Licence:** MIT (repo).
- **Download:** https://github.com/BCV-Uniandes/SpaRED ; datasets distributed via a **Google Drive folder** linked from the README; installable as a PyPI package (`spared`).
- **Biggest limitation:** distribution via a personal Google Drive link is fragile and unversioned — no DOI, no checksums, no guarantee of persistence. Also, the imputed/completed values are model outputs; training on SpaCKLE-completed data and then evaluating gene-expression prediction risks circularity that must be handled explicitly.

---

## 7. HER2ST (Andersson et al.) — the classic breast cancer ST cohort

**Citation (verified via PubMed):** Andersson A, …, Lundeberg J. *Spatial deconvolution of HER2-positive breast cancer delineates tumor-associated cell type interactions.* Nat Commun. 2021;12(1):6012. DOI 10.1038/s41467-021-26271-2. PMID 34650042.

- **Size:** **36 breast cancer sections** (HER2-positive). Exact patient count not stated on the repo page — commonly cited as 8 patients (A–H) but **COULD-NOT-VERIFY** from a primary source in this session.
- **Modalities and pairing:** count matrices + **H&E images** + **pathologist annotations** per section, plus spot-selection tables. The pathologist annotations are the rare asset — most ST datasets have no expert region labels.
- **Resolution:** legacy ST platform — **100 µm spots, 200 µm centre-to-centre**. Coarse. Roughly 10–40 cells per spot.
- **Licence:** **no explicit licence stated** on the repository. Treat as all-rights-reserved until you confirm with the authors.
- **Download route:** https://github.com/almaan/her2st points to **https://zenodo.org/record/3957257**. Files are **7z-encrypted with passwords published in the repo README**: count matrices and images `zNLXkYk3Q9znUseS`; metadata and spot selection `yUx44SzG6NdB32gY`. Contact listed for access questions.
- **Biggest limitation:** the password-protection scheme is a deliberate friction/consent mechanism, and combined with the absent licence statement it makes redistribution legally murky. Plus 100 µm resolution is now two generations behind.

---

## 8. Mendeley 29ntw7sh4r — 23-patient breast cancer ST ("the STNet data")

**Dataset:** *Human breast cancer in situ capturing transcriptomics.* Stenbeck L, Bergenstråhle L, Lundeberg J, Borg Å. Mendeley Data, DOI **10.17632/29ntw7sh4r.5**, published 2021-11-30.

- **Size:** **triplicate sections from 23 breast cancer patients**, spanning **luminal A, luminal B, triple-negative, and HER2-positive** subtypes. Sample metadata, replicate structure, and **tumour annotations** provided.
- **Licence:** **CC BY 4.0** — permissive, commercial use permitted with attribution. Rare and valuable for a clinical ST cohort.
- **Download:** https://data.mendeley.com/datasets/29ntw7sh4r — public, no gate.
- **Related model paper (verified via PubMed):** He B, …, Zou J. *Integrating spatial gene expression and breast tumour morphology via deep learning.* Nat Biomed Eng. 2020;4(8):827–834. DOI 10.1038/s41551-020-0578-x. PMID 32572199. This is ST-Net.
  > **Flagged uncertainty:** the Mendeley record itself lists its related article as Ståhl et al. Science 2016 (10.1126/science.aaf2403), **not** He et al. The 23-patient breast cohort is universally referred to as "the ST-Net dataset" in the literature, but the direct record-level linkage between this Mendeley DOI and He et al. 2020 is **COULD-NOT-VERIFY** from primary sources this session. Confirm from the He et al. Data Availability statement before citing the DOI as ST-Net's training data.
- **Biggest limitation:** legacy ST resolution (100 µm spots) and 2016-era chemistry — low genes-per-spot and high dropout relative to Visium. Also fully superseded in size by HEST-1k/STimage, which both ingest it.

---

## 9. HumanST-1k (via STAMP) — largest claimed pairing, but no verified public route

**Citation (verified via arXiv API):** Zhou F, Xu Y, Zhang Z, Wang Y, Guo Z, Liang L, Ma J, Jin C, Liu Z, Zhou H, Wang H, Cai D, Zhao C, Wang X, Yang C, Wang Y, Li W, Gao F, Wang Z, Li Z, Zhang X, Liang L, Chen H. *Spatial Transcriptomics-Guided Alignment Enhances Molecular Profiling in Pathology Foundation Model.* arXiv:**2606.03644v1** (2026-05-29).

- **Claimed size:** *"HumanST-1k, a human ST dataset spanning diverse anatomical organs and sequencing platforms. This atlas yields **1.8 million pairs of H&E patches and corresponding transcriptomic profiles**."*
- **Purpose:** used to fine-tune pathology foundation models via pathway-informed alignment (STAMP), aggregating raw transcripts into functional pathways before alignment — an explicit denoising move.
- **Status:** **COULD-NOT-VERIFY** licence, download route, or whether HumanST-1k has been released at all. The paper is very recent (May 2026). Treat as *announced, not obtained*. Do not cite it as an available resource.
- **If released, it would be the largest H&E↔ST pairing available** — worth setting a watch on the arXiv listing / the authors' GitHub (Hao Chen's lab, HKUST).

---

## 10. Aggregator databases (secondary — use for discovery, not as primary corpora)

| Resource | Scale | URL | Citation (verified) |
|---|---|---|---|
| **STOmicsDB** | 218 manually curated datasets, 17 species; includes archiving standards + submission system | https://db.cngb.org/stomics/ | Xu Z, …, Wei X. *Nucleic Acids Res.* 2024;52(D1):D1053–D1061. DOI 10.1093/nar/gkad933. PMID 37953328 |
| **CROST** | 182 ST datasets, **1,033 sub-datasets**, **48,043 tumour-related spatially variable genes**; explicitly tumour-focused | https://ngdc.cncb.ac.cn/crost | Wang G, …, Bao Y. *Nucleic Acids Res.* 2024;52(D1):D882–D890. DOI 10.1093/nar/gkad782. PMID 37791883 |
| **SODB** | **>2,400 experiments** from **>25 spatial omics technologies** in unified format; SOView visualisation | (see paper) | Yuan Z, Pan W, Zhao X, Zhao F, Xu Z, Li X, Zhao Y, Zhang MQ, Yao J. *Nat Methods.* 2023. DOI 10.1038/s41592-023-01773-7 (PMID 36797409); publisher correction DOI 10.1038/s41592-023-01844-9 (PMID 36932185) |

**CROST is the most relevant of the three** for an oncology thesis — it is the only one built around tumour-associated spatially variable genes, and it integrates transcriptome + epigenome + genome to relate SVGs to cancer progression and prognosis.

**Shared limitation of all three:** they index and visualise, but do **not** guarantee paired, registered, quality-controlled H&E for every entry. None of them ships a HEST-style harmonised image–expression interface. Use them to find cohorts, then fetch from source.

---

## 11. Human Tumor Atlas Network (HTAN) — high potential, unverifiable this session

**Status: COULD-NOT-VERIFY.** Both https://humantumoratlas.org/explore and https://docs.humantumoratlas.org/ are JavaScript-rendered single-page apps; plain fetches returned only navigation chrome. `https://docs.humantumoratlas.org/data_access/` returned **404**. A PubMed search for an HTAN data-release paper returned only irrelevant hits.

I therefore make **no claims** about HTAN assay counts, spatial platform coverage, H&E availability, cancer types, or open-vs-controlled (dbGaP) access tiers. HTAN is well known to host multimodal tumour atlases including imaging, and is very likely the highest-value unexplored source on this list for an oncology thesis — but it needs a browser session or the `htan` CLI / Synapse API to characterise properly.

**Action item:** re-scout HTAN using the Synapse API or `claude-in-chrome`, not WebFetch.

---

## What is unavailable, gated, or paywalled

1. **HEST-1k is gated** — HF account + accept terms + auth token. Confirmed: `raw/main/README.md` returns HTTP 401. Not a blocker, but it is a click-through licence acceptance, and the licence is **CC BY-NC-SA 4.0** (non-commercial, share-alike). Same for `hest-bench`.
2. **HER2ST has no stated licence** and is distributed as password-protected 7z archives on Zenodo. Legally ambiguous for redistribution.
3. **STimage-1K4M's HF viewer is broken** (schema mismatch, 17,138 inconsistent gene columns, missing coordinate columns). Data is downloadable but not streamable; you must write your own loader.
4. **10x catalogue is not machine-enumerable** — JS-rendered, returns "Showing 0 datasets" to a fetch. Individual dataset pages work and are CC BY 4.0. Downloads typically behind an email/registration form.
5. **HTAN could not be characterised at all** — see §11.
6. **HumanST-1k (1.8M pairs) has no verified public route** — announced in arXiv:2606.03644 (May 2026) only.
7. **SpaRED data lives on a Google Drive link** — no DOI, no versioning, no checksums.
8. **No paywall blocked anything material.** Every paper cited here was verified through open APIs (arXiv, PubMed E-utilities). Nature Biomedical Engineering and Nature Methods full texts are paywalled, but PubMed abstracts and metadata were sufficient for all claims made.

---

## Recommended acquisition order

1. **STimage-1K4M first** (MIT, ungated, 4.29M pairs, full transcriptome per spot) — write the loader once, it's the largest permissive corpus.
2. **HEST-1k second** (accept the NC-SA terms; subset by oncotree code rather than pulling 2 TB) — it is the field's de-facto standard and reviewers will expect HEST-Bench numbers.
3. **HESCAPE third** if the work touches cross-modal alignment — you need its negative result in your related-work section either way.
4. **10x CC BY 4.0 releases** for any component that must be commercially unencumbered, and for anything requiring true Visium HD (8 µm) or Xenium subcellular resolution.
5. **CROST** for tumour-SVG discovery; **HTAN** as the highest-value re-scout target.

## Citations verified in this session (all confirmed against primary APIs)

- arXiv:2406.16192v2 — HEST-1k (arXiv API)
- arXiv:2406.06393v2 — STimage-1K4M (arXiv API)
- arXiv:2508.01490 — HESCAPE (arXiv API)
- arXiv:2505.02980 and arXiv:2407.13027 — SpaRED / SpaCKLE (arXiv API)
- arXiv:2506.05361 — STFlow (arXiv API; method paper, uses HEST-1k + STimage-1K4M)
- arXiv:2606.03644v1 — STAMP / HumanST-1k (arXiv API)
- PMID 38114474 — Janesick et al., Nat Commun 2023, DOI 10.1038/s41467-023-43458-x (PubMed)
- PMID 34650042 — Andersson et al., Nat Commun 2021, DOI 10.1038/s41467-021-26271-2 (PubMed)
- PMID 32572199 — He et al., Nat Biomed Eng 2020, DOI 10.1038/s41551-020-0578-x (PubMed)
- PMID 37953328 — STOmicsDB, Nucleic Acids Res 2024, DOI 10.1093/nar/gkad933 (PubMed)
- PMID 37791883 — CROST, Nucleic Acids Res 2024, DOI 10.1093/nar/gkad782 (PubMed)
- PMID 36797409 / 36932185 — SODB, Nat Methods 2023, DOI 10.1038/s41592-023-01773-7 (PubMed)
- Mendeley Data DOI 10.17632/29ntw7sh4r.5 (Mendeley record page)
