# PRIOR-ART ASSASSINATION ATTEMPT — Thesis T3

**Target:** "The morphology<->transcriptome channel in human tumours has a measurable, confound-certified spectrum, and we report it for the first time."

**Verdict: NOT KILLED (survives, but wounded on two sub-claims).**

**Date of sweep:** 2026-07-29/30. **Search constraint:** WebSearch exhausted; all evidence below obtained via WebFetch against arXiv API, PubMed E-utilities, Europe PMC REST, and OpenAlex. Semantic Scholar returned HTTP 429 on every attempt (COULD-NOT-USE).

---

## 1. What I searched

| Axis of T3 | Queries run (APIs) |
|---|---|
| Confound-residualised morphology<->expression measurement | PubMed (`confound*`/`pitfall*`/`limitations` in title x `gene expression` x `histology`); Europe PMC (`residualized` + `whole slide` + `gene expression`; `deconfound*`/`confounder-aware`/`confound-controlled` x histology x transcriptome); arXiv (`histology` AND `transcriptome` AND `confounder` -> **0 results**) |
| Spectrum / k-components / permutation null | arXiv (`canonical correlation`/`shared latent`/`permutation null`/`identifiable` x histology x transcriptomic); Europe PMC (`canonical correlation` x whole slide x transcriptom* x cancer); OpenAlex |
| Predictability ceiling / information content | arXiv (`upper bound`/`information content`/`how much information`/`ceiling` x H&E x gene expression) |
| Tissue-composition / purity reduction | arXiv & PubMed & Europe PMC (`cell type composition`, `tumor purity`, `tissue composition` x gene expression prediction) |
| Site/batch/scanner confounding in pathology FMs | arXiv (`pathology foundation model` x `batch effect`/`confound`/`site`) |
| Within-cancer vs cross-cancer inflation | Europe PMC (`within-cancer`; `pan-cancer` + `cancer type` + confounded/inflated) |
| Bulk WSI->RNA on TCGA | arXiv (`bulk RNA` + `whole slide` + `TCGA`); Europe PMC/OpenAlex for SEQUOIA, HE2RNA, Commun Med pan-cancer study |
| Non-tumour analogues (GTEx, Cell Painting) | Europe PMC/OpenAlex for Ash et al., imageQTL, Way et al. |

Roughly 20 distinct API queries, ~450 titles screened.

---

## 2. Verified prior art found (all citations below were retrieved from a live API this session)

### Tier A — genuinely threatening (overlaps a named sub-claim of T3)

**A1. DECAT — "When Are Multimodal Predictions Biologically Supported? A Diagnostic Evaluation Framework."**
Steiner D, Arango-Argoty G, Sun G, Jacob E. arXiv:2605.31504v1, 2026-05-29. VERIFIED (arXiv API + abs page).
Verbatim from abstract: *"accurate prediction does not reveal whether the model has learned biology that is shared across modalities, biology confined to one modality, or spurious correlations that reflect confounders rather than genuine biology. We introduce DECAT, a model-agnostic post-hoc evaluation framework that classifies multimodal representations into four diagnostic scenarios ... using five null-referenced metrics and a rule-based decision procedure ... validated ... on real data from 8,979 TCGA patients, evaluating both multimodal embeddings and five pretrained pathology foundation models ... DECAT detects confounding invisible to AUROC without requiring the confounder labels."*

**Why it hurts:** this is the single closest published thing to T3's framing. Same cohort family (TCGA, n=8,979 vs MORPHEUS n=6,192), same modality pair (pathology FM embeddings + paired RNA), same core question ("is the cross-modal signal shared biology or confounding?"), and it is explicitly **null-referenced** — i.e. it already instantiates the "measure against a null, and be publishable when the answer is 'confounded'" move. It also already reports the field-correcting negative that entangled multimodal models "falsely claim shared biology in the majority of cases where it is absent."

**What it does NOT do:** it is a *diagnostic classifier over model representations for a given downstream task*, not a *measurement of the channel itself*. It does not (i) residualise both modalities on a stated confound battery (site/batch/stage/grade/purity/ESTIMATE/deconvolution/RIN), (ii) estimate a **dimension** k of surviving shared structure, (iii) screen leading components against >10,000 MSigDB/Reactome sets for irreducibility, (iv) replicate loadings by independent refit in CPTAC and GTEx, (v) test a protein/phospho shadow, or (vi) test prognostic increment over stage+grade+subtype. It also does not report a scalar effect size for the residual channel.

**Damage assessment:** kills the phrase *"for the first time"* as applied to **confound-certification of cross-modal claims in TCGA**. T3 must now cite DECAT and position against it. It does **not** kill the measurement.

**A2. HESCAPE — "A Large-Scale Benchmark of Cross-Modal Learning for Histology and Gene Expression in Spatial Transcriptomics."**
Gindra RH, Palla G, Nguyen M, Wagner SJ, Tran M, Theis FJ, Saur D, Crawford L, Peng T. arXiv:2508.01490, 2025-08-02. VERIFIED (arXiv API).
Finds that contrastive pretraining improves mutation classification but **degrades** direct gene-expression prediction, and identifies **batch effects as the primary interfering factor**, calling for batch-robust multimodal methods. 6 gene panels, 54 donors, spatial transcriptomics.

**Damage:** pre-empts a weak version of "the channel is confounded by batch" — but at ST spot level, in a small donor cohort, with no residualisation, no null spectrum, no k. It is a benchmark, not a measurement of a confound-adjusted estimand.

**A3. Wang C, Chan AS, Fu X, Ghazanfar S, Kim J, Patrick E, Yang JYH. "Benchmarking the translational potential of spatial gene expression prediction from histology." Nat Commun 2025; 10.1038/s41467-025-56618-y.** VERIFIED (PubMed esummary, PMID 39934114). This is the paper T3 already cites; the citation is real and the DOI is correct. It is prior art for "benchmark under transfer," not for "confound-certified spectrum."

### Tier B — adjacent, does part of the machinery but in a different setting

**B1. Howard FM et al. "The impact of site-specific digital histology signatures on deep learning model accuracy and bias." Nat Commun 2021; 10.1038/s41467-021-24698-1.** VERIFIED (Europe PMC).
Establishes that TCGA submitting-site signatures are recoverable from H&E and confound downstream molecular predictions. This is prior art for **the existence of the site confounder** and partially for T3's "validity certificate" idea (a model should not be able to read site). It does not measure the residual channel.

**B2. Ash JT, Darnell G, Munro D, Engelhardt BE. "Joint analysis of expression levels and histological images identifies genes associated with tissue morphology." Nat Commun 2021; 10.1038/s41467-021-21727-x** (bioRxiv 10.1101/458711, 2018). VERIFIED (Europe PMC + OpenAlex).
GTEx, normal tissue: joint image-expression association with covariate adjustment. This is the closest published *statistical* precedent for "associate morphology with expression while adjusting for nuisance." Normal tissue, gene-level associations, no tumour confound battery, no k, no external replication design.

**B3. "Machine-learning models based on histological images from healthy donors identify imageQTLs and predict chronological age." PNAS 2025; 10.1073/pnas.2423469122.** VERIFIED (Europe PMC).
**B4. "HistoGWAS: an AI-enabled framework for automated genetic analysis of tissue phenotypes in histology cohorts." Genome Biol 2026; 10.1186/s13059-026-04031-z.** VERIFIED (Europe PMC).
Both do image-phenotype <-> molecular/genetic association with nuisance control in **non-tumour** cohorts. They matter because they occupy the "GTEx replication" ground T3 wants for sub-claim (iii) — T3 must cite them rather than present GTEx morphology-molecular linkage as untrodden.

**B5. "Histology image analysis of 13 healthy tissues reveals molecular-histological correlations." Sci Rep 2025; 10.1038/s41598-025-11853-7.** VERIFIED (Europe PMC). Same category as B3/B4.

**B6. Way GP et al. "Morphology and gene expression profiling provide complementary information for mapping cell state." Cell Systems 2022; 10.1016/j.cels.2022.10.001** (bioRxiv 10.1101/2021.10.21.465335). VERIFIED (Europe PMC).
The canonical *measurement* of the morphology<->transcriptome channel — but in Cell Painting / L1000 cell-line perturbation space, not human tumours. It is the strongest conceptual template for "we measured the channel and it is partly complementary rather than redundant," and a reviewer will ask why T3 is different. Answer: different substrate (human tumour tissue), different confound structure, different estimand.

### Tier C — the prediction/hunt literature T3 explicitly is NOT (but must beat on framing)

- Schmauch B et al. "A deep learning model to predict RNA-Seq expression of tumours from whole slide images." Nat Commun 2020; 10.1038/s41467-020-17678-4 (HE2RNA). VERIFIED.
- Pizurica M et al. "Digital profiling of gene expression from histology images with linearized attention" (SEQUOIA). Nat Commun 2024; **10.1038/s41467-024-54182-5**. VERIFIED (Europe PMC core). 7,584 TCGA tumours / 16 cancer types, generalisation validated on two independent cohorts (1,368 tumours), plus breast-cancer recurrence-risk stratification. **This is the closest asset-overlap paper**: TCGA -> external cohort validation + prognostic stratification. It is a HUNT with no confound residualisation, no null, no k.
- "A systematic pan-cancer study on deep learning-based prediction of multi-omic biomarkers from routine pathology images." Commun Med 2024; 10.1038/s43856-024-00471-5. VERIFIED (OpenAlex). ~12,000 models, ~4,000 biomarkers, 32 cancer types. Reports that predictive accuracy "was not significantly affected by tumor purity." **This is the nearest published statement about purity and predictability** — and it is a crude, per-biomarker, non-residualised claim that T3 can directly contradict/refine. Useful as the strawman rather than the assassin.
- Weitz P et al. "Transcriptome-wide prediction of prostate cancer gene expression from histopathology images using co-expression-based CNNs" (arXiv 2104.09310, 2021) and "Evaluation and Prognostic Validation of Deep Regression Models for WSI-Based Gene-Expression Prediction" (arXiv, 2024-10). VERIFIED via arXiv listing (titles/dates only).
- Pathology-FM robustness line: "Current Pathology Foundation Models are Unrobust to Medical Center Differences" (arXiv, 2025-01-29); "Scanner-Induced Domain Shifts Undermine the Robustness of Pathology Foundation Models" (arXiv, 2026-01-07); "Beyond Counts: A Distributional Robustness Margin For Pathology Foundation Models" (arXiv, 2026-07-28); "Mitigating Batch Effects in Histopathology via Language-Mediated Robust Embedding Generation" (arXiv, 2026-06-27). VERIFIED as arXiv listings.

---

## 3. Explicit negatives (searches that returned nothing)

- arXiv `histology AND transcriptome AND confounder`: **0 results**.
- PubMed `histology + deep learning + gene expression prediction + tumor purity + confounding`: **0 results**.
- PubMed title-restricted `(pitfall|caveat|confound|critical assessment|limitations|benchmark) x (gene expression|transcriptom) x (histology|pathology|whole slide|image)`: **1 result** — Wang et al. 2025 (A3). Nothing else.
- Europe PMC `"residualized" AND "whole slide" AND "gene expression"`: no on-topic hits.
- No paper found anywhere that reports a **number of identifiable components (k)** in the tumour morphology<->transcriptome relationship against a permutation null.
- No paper found that reports a **confound-adjusted scalar effect size** for the tumour morphology<->transcriptome channel (T3's "r ~ 0.07" estimand).
- No paper found that pairs a **negative validity certificate** (component must FAIL to predict site/scanner/purity/RIN) with a **positive replication certificate** (loadings refit in CPTAC and GTEx) for this modality pair.

---

## 4. Why T3 survives

The composite estimand is unoccupied. Specifically, no verified publication does **all or even most** of:
1. residualise **both** modalities on the full battery (cancer type, submitting site/batch, stage, grade, purity/ESTIMATE, deconvolution proportions, RNA-quality proxies) in TCGA at n>6,000;
2. estimate the **dimension** k of surviving cross-modal structure with cross-fitting and a permutation null;
3. pre-registered irreducibility screen of the leading component against >10,000 gene sets + consensus subtypes;
4. a **negative** validity certificate (fails to predict site/scanner/purity/RIN);
5. independent-refit loading replication in CPTAC **and** GTEx;
6. a CPTAC protein/phospho shadow;
7. leakage-controlled held-out-cancer prognostic increment over stage+grade+subtype+known signatures.

Items 1, 2, 3, 5, 6, 7 have no verified precedent for this modality pair in human tumours. Item 4 is partially anticipated by Howard 2021 (site) and, at the framework level, by DECAT.

## 5. Required defensive moves (do these or a reviewer kills you)

1. **Cite and differentiate DECAT (arXiv:2605.31504).** State plainly: DECAT diagnoses whether a *model's representation* for a *given task* is confounded; T3 measures the *channel itself* and reports its dimension and effect size. Consider running DECAT as a baseline diagnostic on the MORPHEUS embeddings — it is model-agnostic and post-hoc, so it is cheap, and passing it strengthens the certificate.
2. **Drop unqualified "for the first time."** Use: "first confound-residualised estimate of the *dimension and magnitude* of the tumour morphology<->transcriptome channel."
3. **Cite HESCAPE** as the ST-level batch-effect precedent, and be precise that T3's estimand is bulk, within-cancer, confound-residualised — not comparable to HESCAPE's spot-level PCCs. (This is consistent with the thesis's own "explicitly NOT claimable" clause about ST vs bulk head-to-head; keep that clause.)
4. **Cite Commun Med 2024 (10.1038/s43856-024-00471-5)** and directly engage its claim that predictability "was not significantly affected by tumor purity." T3's residualisation result is the sharpened correction to that claim — this is a strength, not a threat.
5. **Cite Way et al. Cell Systems 2022** as the cell-line-level precedent for the "channel measurement" genre, and Ash 2021 / imageQTL PNAS 2025 / HistoGWAS Genome Biol 2026 as the GTEx-side precedents. Sub-claim (iii) (GTEx replication) is the weakest novelty leg — it sits on occupied ground. Consider reframing GTEx as a *negative control / specificity check* (loadings should replicate only for tumour-intrinsic components, not composition components) rather than as a discovery claim.
6. **SEQUOIA (10.1038/s41467-024-54182-5)** already does TCGA -> external validation + prognostic stratification. T3's prognostic sub-claim (v) is only novel if it is an *increment over stage+grade+subtype+known signatures with a bootstrap CI on delta C-index in held-out cancers*. Keep that framing verbatim; it is what SEQUOIA does not do.

## 6. Caveats on this sweep

- **Semantic Scholar was unavailable (HTTP 429 on all attempts).** A fraction of CS-venue work (MICCAI/NeurIPS workshop papers) indexed there but not in arXiv/PubMed/Europe PMC/OpenAlex may be missed.
- Europe PMC's relevance ranking is fuzzy-OR; several searches returned noise. Absence of a hit in those runs is weaker evidence than the exact-zero arXiv and PubMed results reported in section 3.
- Titles listed in section 2 Tier C from arXiv listing output have verified titles and dates but I did not open each abstract page; treat their detailed content as unconfirmed.
- No citation in this document was written from memory. Anything I could not retrieve live is not listed.
