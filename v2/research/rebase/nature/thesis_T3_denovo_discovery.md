# Thesis T3 — De novo discovery of morphology-encoded molecular axes

**Date:** 2026-07-29
**Scout:** feasibility + prior-art sweep, no WebSearch (exhausted). Sources actually queried: Europe PMC REST
(`/search`, `/fullTextXML`), OpenAlex `/works`, NCBI E-utilities (esearch/efetch), WebFetch on arXiv abs +
raw.githubusercontent + cancerimagingarchive.net. arXiv API and Semantic Scholar API both returned **HTTP 429 for
the whole session** and could not be used; every claim below is sourced from an index that answered.
**No citation in this document is unverified.** Where a number could not be pinned down it says so.

---

## 0. The claim under test

> Stop predicting the curated canon (saturated, confounded, +0.07 ceiling). Instead **discover de novo molecular
> axes that are morphologically encoded**, and validate them by (a) cross-cohort replication in open data and
> (b) independent prognostic value — not by wet lab.

Verdict up front: **the discovery *mechanism* is thoroughly published; the discovery *epistemics* are not.**
The surviving contribution is not "run an unsupervised joint model on WSI+RNA" — three groups have done that at
larger scale than MORPHEUS would. It is the **certification protocol** that separates a discovered axis from
composition, batch, site, stage and RNA-quality — and the honest measurement of how much signal is left after
that subtraction. That reframing is defensible. The original framing is not.

---

## 1. Has this been done? The closest existing works

### 1.1 Jones, Gundersen, Engelhardt — "Linking histology and molecular state across human tissues"
**bioRxiv 10.1101/2022.06.10.495669 (2022). Still a preprint (Europe PMC PPR506074; no
commentCorrection/published-link recorded; 1–3 citations). Confidence: READ-ABSTRACT-ONLY (verbatim from
Europe PMC core record).**

This is **the single closest prior art to T3's core method** and it was not on the repo ledger.

What they did, verbatim from the abstract: *"we jointly analyzed 13,360 human tissue samples with paired bulk gene
expression profiles and histology images across 935 donors from the Genotype and Tissue Expression (GTEx)
Consortium v8 study. This analysis reveals relationships among gene expression and cellular morphology through
shared sources of expression and morphological heterogeneity both within and between tissue types."*

**The load-bearing sentence — read it twice:** *"We describe shared sources of variation including **cell-type
heterogeneity, sample ischemic time, and donor health and demographics**. We find specific correlated effects in
both morphology and transcription linked to specific donor characteristics, such as **their use of mechanical
ventilation**."*

So: the largest paired-modality unsupervised joint analysis ever run recovered, as its shared morphology↔molecular
axes, **cell composition + a preanalytic variable + demographics + an ICU treatment flag.** Not one novel
biological programme. This is simultaneously the closest prior art and the strongest available statement of T3's
failure mode.

What they did **not** do, and where T3's delta lives: no cancer cohort; no external replication of the recovered
axes in a second cohort; no prognostic evaluation; no confound *certificate* (they name the confounds
descriptively rather than requiring the axis to fail at predicting them); no test of whether an axis is already in
the curated canon.

### 1.2 Quiros et al. — Histomorphological Phenotype Learning (HPL)
**"Mapping the landscape of histomorphological cancer phenotypes using self-supervised learning on unannotated
pathology slides," Nature Communications 15 (2024), doi:10.1038/s41467-024-48666-7, 46 citations.
Confidence: READ-ABSTRACT-ONLY (PubMed 38862472 + Europe PMC).**

Self-supervised, label-free discovery of tile clusters → an "HP-Atlas" of histomorphological phenotypes, showing
benign→malignant trajectories through inflammatory/reactive states. Abstract: clusters *"have distinct features
which can be identified using orthogonal methods, **linking histologic, molecular and clinical phenotypes**"*; in
lung cancer they align with survival, with recognised tumour types and growth patterns, and with **transcriptomic
measures of immunophenotype**; *"These properties are maintained in a multi-cancer study."*

This is de novo morphology-side discovery, explicitly linked to molecular and clinical phenotypes, published in a
Nature-family journal. **Any T3 pitch phrased as "unsupervised morphology clusters that turn out to be
molecularly meaningful" is pre-empted outright.**

Two follow-ups extend it in exactly the directions T3 would want to claim:
- **Colon.** "Self-supervised learning reveals clinically relevant histomorphological patterns for therapeutic
  strategies in colon cancer," Nature Communications (2025), PubMed 40057490. Barlow Twins on 435 TCGA COAD WSIs,
  Leiden → 47 HPCs; *"HPC reproducibility and predictive ability for overall survival are confirmed in an
  **independent clinical trial (N = 1213 WSIs)**"*; gene-set enrichment and immune-landscape analysis per HPC.
  **This already executes the "discover de novo, replicate in an independent cohort, show prognostic value"
  template.**
- **Kidney.** "Unbiased self supervised learning of kidney histology reveals phenotypic and prognostic insights,"
  Scientific Reports (2025), PubMed 41062686. HPCs transferred to five internal/external validation sets
  (1,421 WSIs), reproducibly associated with pathologist quantifications (interstitial fibrosis AUC = 0.83), and
  specific HPCs predicted longitudinal kidney-function decline.
- **Colorectal, multicentre.** "Multicenter self-supervised computational pathology identifies prognostic
  histomorphological phenotypes in colorectal cancer," bioRxiv PPR1282959 (2026), doi 10.64898/2026.07.15.738753.
  HPL-PanColon; *"reduced institution- and dataset-specific batch effects"* vs general-purpose pathology FMs;
  leave-one-institution-out on 1,024 patients; derives the Colon Histomorphology Prognostic Score (CHiPS).
  **Note the batch-effect framing — the confound axis T3 wants to own is already being claimed here.**

### 1.3 Fu et al. — Pan-cancer computational histopathology
**Nature Cancer 1 (2020), doi:10.1038/s43018-020-0085-8, 445 citations. Confidence: READ-ABSTRACT-ONLY.**

17,355 H&E WSIs, 28 cancer types, transfer-learned features correlated against matched genomic, transcriptomic and
survival data. Recovers whole-genome duplication signatures, aneuploidies, focal amps/dels, driver mutations.

**The sentence that most directly threatens T3:** *"There are widespread associations between bulk gene expression
levels and histopathology, **which reflect tumor composition** and enable the localization of transcriptomically
defined tumor-infiltrating lymphocytes."*

The pan-cancer morphology↔bulk-expression association sweep is **done**, at 2.8× MORPHEUS's slide count, and its
own authors attribute the associations to composition. T3 cannot claim the sweep. It can only claim what survives
composition adjustment — and it must beat the prior that nothing does.

### 1.4 HistoGWAS — automated discovery of morphology-linked molecular variation
**"HistoGWAS: an AI-enabled framework for automated genetic analysis of tissue phenotypes in histology cohorts,"
Genome Biology (2026), doi:10.1186/s13059-026-04031-z, Europe PMC 41918129. Confidence: READ-ABSTRACT-ONLY.**

*"foundation models for automated trait definition, variance component models for efficient association testing,
and generative models for variant effect interpretation. Applied to 11 tissues from the Genotype-Tissue Expression
project, HistoGWAS identifies four genome-wide significant loci associated with tissue histology — tissueQTLs —
which we link to molecular changes and complex traits."*

This is a **published, general-purpose framework for de novo discovery of molecular correlates of FM-derived
histology traits**, with a formal significance framework and a variance-component null. The molecular side is
germline genotype rather than a transcriptomic axis — that is T3's only structural gap against it — but the
"automated trait definition → association → molecular interpretation" pipeline is occupied, on GTEx, with a
statistical null already worked out. **A reviewer will ask why T3 is not just HistoGWAS with expression on the
right-hand side.** Have an answer.

### 1.5 The supervised-target saturation evidence (context, not collision)
- **Diao et al., "Human-interpretable image features derived from densely mapped cancer pathology slides predict
  diverse molecular phenotypes," Nature Communications 12 (2021), doi:10.1038/s41467-021-21896-9, 176 citations.**
  >1.6M pathologist annotations across >5,700 samples → 607 HIFs across five cancer types, predicting diverse
  molecular signatures at **AUROC 0.601–0.864**, *"with performance comparable to 'black-box' methods."*
  The interpretable-morphology→known-molecular-target map is saturated.
- **Winter et al., arXiv:2606.29949** (already adjudicated in `NEAR_COLLISIONS.md` §c, READ-FULL-TEXT there):
  frozen WSI+RNA adapter, per-hallmark R² spectrum (G2M 0.78, IFN-γ 0.75 → OxPhos 0.22, FA metabolism 0.20),
  cross-cohort R@10 collapse to 10.8%/4.0%. This is the "which curated programmes are morphologically visible"
  curve, already published.
- **SPADE**, Nucleic Acids Research 49 (2021), doi:10.1093/nar/gkab095, 42 citations — discovers genes correlated
  with image latent features at ST spot level ("Discovery of molecular features underlying the morphological
  landscape…"). The *spot-level* version of T3's discovery step exists and is five years old.
- **Way et al., "Morphology and gene expression profiling provide complementary information for mapping cell
  state," Cell Systems 13 (2022), doi:10.1016/j.cels.2022.10.001, 123 citations.** Cell Painting + L1000 on A549
  across 1,327 compounds × 6 doses: the two modalities are *"partially shared but also complementary."*
  This is the best published *support* for T3's premise — morphology encodes state that expression signatures
  miss — but it is in cell lines under perturbation, not tissue.

### 1.6 Verdict on novelty
| Component of T3 | Status |
|---|---|
| Unsupervised joint factor/axis discovery on paired WSI+bulk RNA | **DONE** — Jones/Engelhardt on GTEx, n=13,360 |
| De novo morphology phenotypes linked to molecular + survival, multi-cancer | **DONE** — HPL, Nat Commun 2024 |
| Discover → replicate in independent cohort → prognostic value | **DONE** — HPL-colon, Nat Commun 2025 (N=1,213) |
| Pan-cancer morphology↔bulk-expression association sweep | **DONE** — Fu, Nature Cancer 2020 (17,355 WSIs) |
| Automated FM-trait → molecular association with a formal null | **DONE** — HistoGWAS, Genome Biol 2026 |
| Spot-level image-latent → gene discovery | **DONE** — SPADE, NAR 2021 |
| **Requiring the axis to FAIL at predicting composition / site / stage / RIN (certificate)** | **NOT FOUND** |
| **Screening a discovered axis against the curated canon for restatement** | **NOT FOUND** |
| **Reporting the full spectrum of the channel, not just the top component** | **NOT FOUND** |

Targeted negative searches (Europe PMC, all HITS: 0): `pathology foundation model AND gene expression AND
confound`; `morphology AND cell type composition AND deconvolution AND whole-slide images`; `discover AND novel
AND morphology AND transcriptomic axis`; `multimodal AND factor analysis AND histopathology AND RNA`.
These are weak negatives (abstract-field only) but consistent.

---

## 2. Is the core claim true? — and the T1 sub-question

### 2.1 The T1 question, answered with numbers
> *Is there published evidence that spot-level ST targets yield materially higher morphology-predictability than
> bulk targets?*

**No — not for the same estimand. The apparent gap is an estimand swap, and it collapses under transfer.**

**Evidence A — HEST-Benchmark leaderboard** (Jaume et al., "HEST-1k: A Dataset for Spatial Transcriptomics and
Histology Image Analysis," arXiv:2406.16192; abstract verified via arXiv abs page; leaderboard numbers verified
via `raw.githubusercontent.com/mahmoodlab/HEST/main/README.md`). Dataset: 1,229 ST profiles (1,276 in the current
library) each paired to a WSI, 153 cohorts, 26 organs, 367 cancer samples from 25 cancer types, 2.1M
expression–morphology pairs. Benchmark: 25 pathology FMs, **ridge regression on PCA-reduced patch embeddings**,
nine cancer/organ tasks, Pearson correlation. Top results: **H-Optimus-1 0.4229**, GenBio-PathFM 0.4197,
**H-Optimus-0 0.4150**, UNI2-h 0.4141, Virchow 0.4061.

Note carefully: MORPHEUS's own patch encoder (H-Optimus) sits at the top of that board, and the spread across
25 foundation models is **0.4229 → 0.4061 for the top five** — i.e. HEST reproduces MORPHEUS's own
**method-invariance** finding in the spatial setting. Whatever is being measured, the encoder barely matters.

**Evidence B — the independent benchmark** ("Benchmarking the translational potential of spatial gene expression
prediction from histology," Nature Communications 16 (2025), doi:10.1038/s41467-025-56618-y, PMC11814321,
24 citations. **Confidence: READ-FULL-TEXT** via Europe PMC fullTextXML). Eleven methods, five SRT datasets,
external validation on TCGA. Verbatim numbers:

| Setting | Metric | Value |
|---|---|---|
| Within-image spot-level, 4-fold CV, HER2+ (785 genes) & cSCC (997 genes) | best method EGNv2, mean PCC over genes | **0.28** (MI 0.06, SSIM 0.22, AUC 0.65) |
| Cross-subtype transfer (train Visium-Hercep-Test2+ → test Visium-HER2+) | best method EGNv1, avg correlation | **0.11** (HVGs) / **0.13** (SVGs) |
| TCGA-BRCA, predicted pseudobulk vs true bulk, **patient-level** (n=671) | mean over methods | 0.53 (HisToGene/EGNv2), 0.52, 0.49, 0.44, 0.42, 0.39, 0.37, **0.33** (ST-Net) |
| Survival risk-group separation from predicted pseudobulk | log-rank | *"borderline statistically significant for three methods"* (DeepSpaCE, Hist2ST, ST-Net) |
| Image-QC dependence | corr(pred. performance, H&E QC PC1) | DeepPT **r = 0.53** (undesirable); HisToGene 0.05, Hist2ST 0.19 |

The paper's own words on generalisation: *"models trained on a specific breast cancer subtype struggle to
generalise to others, as the biological characteristics differ between subtypes."*

**Why the 0.28–0.42 numbers do not refute MORPHEUS's +0.07.** Three separate estimand problems:

1. **Spot-level Pearson pools within-slide variance.** Within a slide, the dominant axis of expression variation is
   tumour vs stroma vs immune vs necrosis — i.e. *tissue composition*, which H&E is a direct optical readout of.
   That is not molecular-state prediction; it is segmentation with a transcriptomic label. Fu et al. said this in
   2020 (associations *"reflect tumor composition"*), and no ST benchmark controls for it.
2. **The "patient-level correlation" of 0.33–0.53 in the Nature Communications benchmark is computed across genes
   within a patient.** Read the methods: *"we calculated the correlation between predicted GE pseudobulk and bulk
   GE for each image"*, restricted to *"log-transformed gene values that were greater than 5"*. That is the
   correlation of the **profile shape** — the average expression profile — and it is exactly the
   ~46–49% cohort-structure artefact MORPHEUS has already characterised. A constant per-gene mean predictor scores
   well on it. It is **not** a cross-patient control-adjusted statistic and must never be compared to +0.07.
3. **Under a genuine distribution shift, spot-level performance drops to 0.11–0.13** — the same order as the bulk
   control-adjusted signal, from the same benchmark, in the same modality.

**Conclusion for T1:** the literature contains **no** head-to-head evidence that spot-level ST targets carry more
*cross-patient, composition-adjusted* morphology-predictable molecular signal than bulk targets. What it shows is
that spot-level tasks have a large, easy, within-slide composition component that inflates the headline number and
does not survive transfer. **If T1 is pitched on "ST targets have a higher ceiling," it is pitched on an
artefact.** The defensible version of T1 is narrower: *ST resolves where in the tissue an axis lives*, which is a
**localisation** argument, not a **predictability** argument.

### 2.2 Is T3's own claim plausible?

**The premise is half-right.** "The curated-canon prediction task is saturated and confounded" is well supported
(Diao AUROC 0.601–0.864 on known signatures; Winter's hallmark R² spectrum; Fu's composition attribution;
MORPHEUS's own method-invariance at +0.07). **"Therefore de novo axes will be richer" does not follow**, and the
strongest evidence available cuts against it:

- **The +0.07 is an average over targets; discovery seeks the argmax.** This is the honest steelman and it is the
  only reason T3 is not dead. If the morphology↔expression channel's control-adjusted canonical-correlation
  spectrum is flat, +0.07 is also the ceiling for the best de novo direction and T3 collapses. If the spectrum has
  a heavy head, the leading direction can be several times +0.07 while the mean over 50 hallmarks stays at 0.07.
  **This is a measurable, falsifiable, and — critically — currently unpublished quantity.** Nobody has reported
  the spectrum. Measuring it is the single highest-value thing T3 can do, and it is cheap.
- **Against:** the two largest paired studies both report that the leading shared axes are composition and
  nuisance. Jones/Engelhardt: cell-type heterogeneity, ischemic time, demographics, mechanical ventilation.
  Fu: tumour composition.
- **Against, hard:** **PathQC** (Bioengineering 13 (2026), doi:10.3390/bioengineering13010005; preprint bioRxiv
  10.1101/2025.09.29.679347) predicts, from H&E alone using UNI features on GTEx (25,306 samples, 29 tissues,
  970 donors), **RNA Integrity Number at mean r = 0.47** and **autolysis at r = 0.45** (adrenal R = 0.82).
  RNA quality is *directly legible in morphology*, and RNA quality *directly distorts the transcriptome*.
  Any morphology↔expression axis discovered in GTEx must be shown not to be RIN. Nobody has done this. It is also
  a real hazard in TCGA (FFPE/frozen, ischemia, section quality).
- **Against, structural:** Howard et al., "The impact of site-specific digital histology signatures on deep
  learning model accuracy and bias," Nature Communications 12 (2021), doi:10.1038/s41467-021-24698-1,
  294 citations — TCGA submitting site is trivially decodable from H&E, survives stain normalisation, and biases
  survival/mutation/stage predictions. Site is a live confound for every TCGA-discovered axis.
- **For:** Way et al. (Cell Systems 2022) establishes in a controlled perturbation setting that morphology and
  expression are *complementary*, not redundant. This is the cleanest published existence proof that
  morphology-encoded state exists outside the expression canon. It is in A549 cells, not tissue.

**Net:** the claim is *plausible but the prior is unfavourable*, and the unfavourable prior is quantified and
published. T3 is a legitimate high-risk bet only if it is run as a **measurement** ("how much non-confound signal
is in this channel, and does the top of it replicate?") rather than as a **hunt** ("we will find a new axis").
Run as a hunt, the modal outcome is a beautifully validated rediscovery of tumour purity.

---

## 3. Open-access data: what exists, what it really gives you

### 3.1 Discovery cohort
**TCGA** — on disk: 6,192 patients, paired WSI (H-Optimus patch features, uncapped) + bulk RNA (BulkFormer),
32 cancers, leakage-controlled 14-dev/21-test held-out-cancer split. Adequate. Caveats: slide and RNA aliquot are
adjacent pieces, not the same tissue (composition mismatch is baked in); site confound per Howard et al.

### 3.2 Replication cohorts that actually exist and are actually open

| Cohort | Paired content | Real size | Access | Hard limitations |
|---|---|---|---|---|
| **CPTAC** (via TCIA for WSIs, GDC for RNA-seq, PDC for proteomics) | H&E WSI + RNA-seq + proteomics/phospho | **~2,416 subjects with public WSIs on TCIA**, summed across 14 collections: CCRCC 262, UCEC 250, LUAD 244, LSCC 212, HNSCC 207, GBM 200, AML 180, PDA 168, STAD 168, BRCA 134, COAD 106, OV 102, CM 95, SAR 88 | Public, except **CPTAC-HNSCC and CPTAC-GBM are "Limited" access** on TCIA | Per-cancer n≈100–260, i.e. **~10× smaller than TCGA per cancer**; survival follow-up immature; **the exact WSI∩RNA-seq intersection is NOT the subject count and must be computed — do not quote 2,416 as the paired n** |
| **GTEx v8/v10** | H&E WSI + bulk RNA from the *same* tissue sample | **25,306 samples, 29 tissues, 970 donors** (PathQC figure); ~30,000 specimens from ~1,000 donors (Arch Pathol Lab Med 2025, doi:10.5858/arpa.2023-0467-oa); Jones/Engelhardt used 13,360 samples / 935 donors; RNAPath used 838 donors / 23 tissues | Fully open (GTEx portal expression matrices; GTEx Histology Viewer images) | **Non-cancer, no outcomes.** Post-mortem: ischemic time, autolysis and RIN are first-order axes (PathQC r≈0.45–0.47). Its value is as a **negative control / confound stress test**, not as a prognostic replication set |
| **HEST-1k** | H&E WSI + spatial transcriptomics | 1,229 profiles (1,276 in current library), 153 cohorts, 26 organs, 367 cancer samples / 25 cancer types, 2.1M expression–morphology pairs | HuggingFace, **CC-BY-NC-SA 4.0 (non-commercial)**; full dataset >2TB, queryable subsets | Few patients per cohort; heterogeneous platforms; NC licence. Use for **spatial localisation of an axis**, not for statistical replication |
| **Molecular-pole-only replication** (no images needed) | bulk RNA + outcome | METABRIC (~2,000, cBioPortal), ICGC/PCAWG, recount3/GEO | Open | Tests only half the axis — but it is the *cheapest and statistically strongest* half, and reviewers accept it |
| **HTAN** | imaging + scRNA/ST | open-access tier varies by atlas | Open tier | Paired *bulk* RNA + WSI + outcome is thin; verify per atlas before relying on it |
| **MOSAIC** (Owkin, used by MISO, doi:10.1038/s41467-025-66691-y, 348 samples / 5 indications) | WSI + spTx | — | **NOT open** | Do not plan around it |
| **POSEIDON** (used by Winter et al., n=265 NSCLC) | — | — | **NOT open** | Do not plan around it |

### 3.3 The data verdict, stated bluntly
**There is essentially one substantial open external cohort with paired WSI + bulk RNA + cancer outcomes, and it
is CPTAC, and it is small and its follow-up is immature.** Everything else is either non-cancer (GTEx), tiny and
non-commercially licensed (HEST), or images-free (METABRIC/ICGC/recount3).

The only replication design that actually closes is a **split replication**:
1. **Morphology pole** → refit the axis independently in CPTAC and in GTEx; test that the *loading* replicates
   (Procrustes / cross-fitted correlation of loadings), not just that some correlation exists.
2. **Molecular pole** → score the axis's gene weights in METABRIC / ICGC / recount3 and test prognostic value at
   n in the thousands, where you actually have power.
3. **Prognostic independence** → TCGA held-out cancers (leakage-controlled split already built) + CPTAC where
   follow-up permits, with ΔC-index over a stage+grade+subtype+known-signature base model, bootstrap CI.
4. **Orthogonal-modality corroboration** → CPTAC proteomics/phospho (already inventoried in the repo). If a
   morphology-encoded axis is real biology, it should have a protein shadow. **This is the strongest non-wet-lab
   validation available to this project and it is under-used.**
5. **Spatial corroboration** → HEST-1k: does the axis's high-scoring morphology localise to a coherent tissue
   compartment? This is where T1's ST assets earn their keep — as *localisation*, per §2.1.

---

## 4. Highest honest claim, and what falsifies it

### 4.1 The claim
> **The morphology↔transcriptome channel in human tumours has a measurable, confound-certified spectrum. We report
> it for the first time. After removing cancer type, submitting site/batch, stage/grade, tumour purity and
> cell-type composition, and RNA quality, the channel retains *k* identifiable components; the leading component
> is not reducible to any catalogued gene set, its loadings replicate when refit independently in CPTAC/GTEx, it
> has a protein shadow in CPTAC proteomics, and it adds prognostic information beyond stage, grade, subtype and
> known signatures in held-out cancers.**

Note the shape: **the paper is publishable even if k = 0.** "We measured the spectrum of the morphology–molecular
channel under full confound control and it is flat at r ≈ 0.07 — the +0.07 is not a floor imposed by weak methods,
it is the channel" is a genuine, citable, field-correcting result that kills a large amount of ongoing work. That
is the insurance policy, and it is the reason to run T3 as a measurement. It is worth Nature Communications /
Nature Cancer on its own. Do not design the study such that a null is a failure.

### 4.2 What would falsify it — ranked by probability of actually happening
1. **The axis is composition.** Regress the axis on CIBERSORTx/xCell/ESTIMATE/purity; if R² > ~0.6, it is
   composition in new coordinates. *This is the modal outcome.* Precedents: Fu 2020 ("reflect tumor composition"),
   Jones/Engelhardt ("cell-type heterogeneity").
2. **The spectrum is flat.** Cross-fitted, confound-adjusted canonical correlations show no head. Then the argmax
   direction ≈ the mean over curated targets ≈ +0.07, and the entire "discovery beats prediction" premise dies.
   Directly implied by MORPHEUS's own method-invariance result. **Measure this first, before building anything.**
3. **The certificate fails.** The axis predicts submitting site / scanner / batch above chance (Howard et al.
   shows site is trivially decodable) or predicts RIN/autolysis above chance (PathQC shows morphology carries
   r≈0.47 of RIN). Either kills it.
4. **Restatement, not discovery.** |r| > 0.8 against any MSigDB Hallmark / Reactome score, consensus molecular
   subtype, proliferation/immune/EMT/hypoxia meta-signature, or the first few PCs of the expression matrix.
   Pre-register the screen and the threshold; run it against ≥10,000 gene sets, not 50.
5. **No replication.** Loadings refit in CPTAC/GTEx correlate with the TCGA loadings at chance. Note this is
   *underpowered* by design given CPTAC's per-cancer n — a null here is ambiguous, so power-analyse before
   claiming it either way.
6. **No prognostic increment.** ΔC-index bootstrap CI over stage+grade+subtype+known-signatures includes 0.
   Given Winter's transfer collapse and the ST benchmark's *"borderline statistically significant"* survival
   result, expect a small increment at best.
7. **Reviewer collision.** "This is HistoGWAS with expression instead of genotype" / "this is HPL with a
   continuous score" / "Jones & Engelhardt did the joint GTEx analysis in 2022." All three are fair. The answer
   must be the certificate + the spectrum + the canon-screen, and it must be in the abstract.

### 4.3 What would make it credible without wet lab
From the precedents that succeeded, in order of how much they buy:
1. **Independent-cohort replication with a pre-registered analysis** — the standard set by HPL-colon
   (Nat Commun 2025, N=1,213 independent trial WSIs) and by Beck et al., "Systematic analysis of breast cancer
   morphology uncovers stromal features associated with survival," Science Translational Medicine 3 (2011),
   doi:10.1126/scitranslmed.3002564, 457 citations — the canonical precedent for *"a de novo morphology feature is
   a real discovery because it replicated and it was prognostic,"* with no wet lab.
2. **Orthogonal modality.** A protein shadow in CPTAC. Morphology→RNA→protein agreement is much harder to explain
   as an artefact than morphology→RNA alone.
3. **Adversarial nulls.** The confound certificate (must *fail* at site/RIN/purity), a permutation null, and a
   text-prior null (is the axis recoverable from gene *names* alone, GenePT-style).
4. **REMARK compliance** for the prognostic claim (doi:10.1093/jnci/djy088, JNCI 2018; doi:10.1371/journal.pmed.1001216).
   Reviewers of a prognostic-biomarker claim will apply it whether or not you cite it.
5. **Spatial localisation** in HEST-1k — showing *where* the axis lives converts a number into a picture a
   pathologist can adjudicate, which is what makes the "new biology" claim land.

---

## 5. Ratings

| Axis | Score | Reasoning |
|---|---|---|
| **Novelty** | **3 / 5** | The discovery *mechanism* is fully occupied: Jones/Engelhardt (joint GTEx factor analysis, n=13,360), HPL (Nat Commun 2024 + colon/kidney/CRC follow-ups with independent-cohort replication and prognostic scores), Fu (pan-cancer sweep, 17,355 WSIs), HistoGWAS (published discovery framework on GTEx), SPADE (spot-level, 2021). What is genuinely unoccupied is the **epistemics**: a confound *certificate* (axis must FAIL at site/purity/composition/RIN), a **canon-restatement screen**, and **reporting the whole spectrum rather than the top component**. Real, defensible, but it is an evaluation contribution wearing a discovery costume. Claiming otherwise gets it desk-rejected. |
| **Feasibility (open data)** | **3 / 5** | Discovery is trivial — assets are on disk and the compute is a rounding error. The bottleneck is replication: **one** substantial open cohort with paired WSI+RNA+outcome (CPTAC, ~100–260 subjects/cancer, immature follow-up, and the WSI∩RNA intersection is smaller than the subject count); GTEx is large and open but non-cancer and outcome-free; HEST is small and non-commercially licensed. The split-replication design (morphology pole in CPTAC/GTEx, molecular pole in METABRIC/ICGC/recount3, protein shadow in CPTAC-PDC, prognostic increment in held-out TCGA cancers) closes the loop but is fiddly and each leg is individually underpowered. |
| **Ceiling** | **4 / 5** | *Conditional on a hit*: a replicated, composition-adjusted, prognostically-independent, previously-uncatalogued morphology–molecular axis with a protein shadow is Nature Cancer / Nature Medicine material — the Beck-2011 template at pan-cancer scale with modern instrumentation. Not Nature-main without wet-lab or trial-level validation. *Unconditionally*: the null result ("the channel is flat at 0.07 under full confound control") is itself a strong Nature Communications / Nature Cancer paper that invalidates a research programme, so the **floor is unusually high** — which is the best argument for running this thesis at all. |

---

## 6. Concrete next step (cheapest decisive experiment)

Before any architecture work: **estimate the confound-adjusted canonical-correlation spectrum of the
morphology↔transcriptome channel** on the 6,192 TCGA cases already on disk.

- Residualise both modalities on: cancer type, submitting site, stage, grade, purity/ESTIMATE, cell-type
  deconvolution proportions, and available RNA-quality proxies.
- Cross-fitted sparse CCA / joint factor model, held-out-cancer split (the existing 14/21 split).
- Report the **whole eigenvalue spectrum with a permutation null**, not the top component.
- Decision rule, pre-registered: if the top confound-adjusted, cross-fitted component does not exceed ~2–3×
  the +0.07 baseline, **T3 becomes the null/measurement paper** and stops being a discovery hunt.

This costs days, not months, and it converts an unfalsifiable ambition into a decision. It also produces the
paper's central figure either way.

---

## 7. Citation hygiene log

**Verified this session** (index that answered in parentheses):
Jones/Gundersen/Engelhardt bioRxiv 10.1101/2022.06.10.495669 (Europe PMC PPR506074, core record incl. author list
and full abstract) · Quiros et al. Nat Commun 2024 doi:10.1038/s41467-024-48666-7 (PubMed 38862472) · colon HPL
Nat Commun 2025 (PubMed 40057490) · kidney HPL Sci Rep 2025 (PubMed 41062686) · HPL-PanColon bioRxiv
10.64898/2026.07.15.738753 (Europe PMC PPR1282959) · Fu et al. Nature Cancer 2020 doi:10.1038/s43018-020-0085-8
(Europe PMC 35122049) · HistoGWAS Genome Biol 2026 doi:10.1186/s13059-026-04031-z (Europe PMC 41918129) · Diao
et al. Nat Commun 2021 doi:10.1038/s41467-021-21896-9 (Europe PMC 33712588) · SPADE NAR 2021 doi:10.1093/nar/gkab095
(Europe PMC 33619564) · Way et al. Cell Systems 2022 doi:10.1016/j.cels.2022.10.001 (Europe PMC 36395727) ·
HE2RNA Nat Commun 2020 doi:10.1038/s41467-020-17678-4 (Europe PMC 32747659) · Tangle CVPR 2024 /
arXiv:2405.11618 (OpenAlex, doi:10.1109/cvpr52733.2024.00920) · HEST-1k arXiv:2406.16192 (arXiv abs page +
mahmoodlab/HEST README leaderboard) · ST-prediction benchmark Nat Commun 2025 doi:10.1038/s41467-025-56618-y
(**full text**, PMC11814321) · RNAPath Nat Commun 2024 doi:10.1038/s41467-024-50317-w (Europe PMC 39003292) ·
PathQC Bioengineering 2026 doi:10.3390/bioengineering13010005 + bioRxiv 10.1101/2025.09.29.679347 (Europe PMC
41595938 / PPR1095001) · GTEx histology QC Arch Pathol Lab Med 2025 doi:10.5858/arpa.2023-0467-oa (Europe PMC
38797720) · Howard et al. Nat Commun 2021 doi:10.1038/s41467-021-24698-1 (OpenAlex) · Beck et al. Sci Transl Med
2011 doi:10.1126/scitranslmed.3002564 (Europe PMC 22072638) · REMARK JNCI 2018 doi:10.1093/jnci/djy088 and PLoS
Med 2012 doi:10.1371/journal.pmed.1001216 (Europe PMC 29873743 / 22675273) · SEQUOIA Nat Commun 2024
doi:10.1038/s41467-024-54182-5 · MISO Nat Commun 2025 doi:10.1038/s41467-025-66691-y · STPath npj Digit Med 2025
doi:10.1038/s41746-025-02020-3 · CPTAC TCIA collection subject counts (cancerimagingarchive.net/browse-collections).

**Carried from the repo, adjudicated previously in `NEAR_COLLISIONS.md`, not re-verified here:**
Winter et al. arXiv:2606.29949 (READ-FULL-TEXT there); SurvPath arXiv:2304.06819.

**COULD-NOT-VERIFY / not attempted:** arXiv API and Semantic Scholar API were HTTP-429 for the entire session;
no claim in this document depends on either. bioRxiv `api.biorxiv.org` did not resolve (DNS); bioRxiv records were
obtained through Europe PMC instead. Publication status of Jones/Gundersen/Engelhardt beyond the 2022 preprint
**could not be confirmed** — Europe PMC records no comment/correction link. Do not assert it is unpublished; assert
only that no published version was found.

**Numbers deliberately NOT asserted:** the CPTAC WSI∩RNA-seq paired sample count (subject counts are TCIA
collection sizes, which are an upper bound); exact HTAN open-tier paired-modality counts.
