# Prior-Art Kill Attempt — T2_virtual_perturbation

**Date:** 2026-07-29
**Verdict:** **NOT KILLED** (survives, with three mandatory narrowings and one citation correction)
**Method:** WebSearch exhausted. All hits verified via Europe PMC REST, NCBI E-utilities (esearch/esummary/efetch), and the arXiv API. Semantic Scholar returned HTTP 429 for every call and was not used. Every citation below was retrieved live; anything unverified is flagged COULD-NOT-VERIFY.

---

## 0. CITATION CORRECTION (important — the thesis mis-cites its own strongest prior art)

The thesis says the phrase "image-only virtual perturbation" is "already taken by **Coladan, Genome Medicine 2026**."

**"Coladan" is not an author. It is the METHOD NAME.** PubMed author search for `Coladan[Author]` returns **0 results** — "Coladan" is explicitly listed by PubMed among phrases not found. The real citation is:

> **Wang Z, Yang C, Tang X, Yin E, Yao Y, Luo Y, He J, Sun N.** "Trimodal, uncertainty-guided whole-slide framework for genome-scale spatial expression and image-only virtual perturbation in cancer cohorts." *Genome Medicine* 2026. PMID **42449400**. DOI **10.1186/s13073-026-01713-y**.
> Method names: **Coladan** (framework) and **Coladan-human3K** (~3,000-profile ST resource).

Given that three fabricated citations already contaminated this project's earlier sweep, flag this: the thesis was one step away from producing a fourth ("Coladan et al."). The paper is REAL; the attribution was wrong.

**Verified abstract facts (fetched verbatim from PubMed):** 32 Visium datasets; Pearson **0.230 → 0.431 (~1.9×)**; zero-shot transfer to VisiumHD and spot-level Xenium; "**CLS embedding-only perturbation performs on par with expression-based baselines, enabling image-only virtual perturbation without measured expression**"; illustrated on normal and cancer prostate sections for in-situ hypothesis generation.

The thesis's numbers for Coladan are correct. Its author attribution is not.

---

## 1. What I tried to kill it with, and what I actually found

### 1.1 Direct hit search: H&E → genome-scale CRISPR dependency in patients — **NOTHING FOUND**

Queries run (all returned zero on-target):

| Source | Query | On-target hits |
|---|---|---|
| PubMed | `("whole slide"[tiab] OR histopathology[tiab] OR "histology images"[tiab]) AND ("gene dependency"[tiab] OR "gene essentiality"[tiab] OR DepMap[tiab] OR "CRISPR screen"[tiab])` | 5 total, **0 on-target** (GenAR, STPath, an influenza CRISPR screen, Perturb-DBiT, Dhainaut spatial CRISPR) |
| Europe PMC | `"H&E" AND (DepMap OR Perturb-seq OR "CRISPR screen") AND predict*` | **0 on-target** |
| Europe PMC | `(histolog* OR "whole-slide" OR "H&E" OR pathology) AND ("dependency map" OR "gene effect" OR CERES OR Chronos) AND ("deep learning" OR "machine learning") AND TCGA` | **0 on-target** |
| arXiv | `all:"histopathology" AND all:"CRISPR"` | **0 total results** |
| arXiv | `all:"whole slide" AND all:"perturbation" AND all:"dependency"` | 3 results, all MIL-interpretability, **0 on-target** |
| arXiv | `all:"in silico perturbation" AND all:"pathology"` | **0 total results** |
| arXiv | `all:"virtual perturbation" AND all:"histology"` | **0 total results** |
| PubMed | `"virtual perturbation" AND (image OR histology OR morphology)` | 8 total; only **Coladan (42449400)** on-target |

**Conclusion: no published work predicts DepMap gene-effect scores or Perturb-seq responses for patient tumors from H&E.** The core object of T2 is unoccupied.

### 1.2 The nearest neighbours that DO exist (and must be cited)

**(a) Per-patient dependency prediction from OMICS (not images) — the estimand is taken; the modality is not.**

- **Chiu YC, Zheng S, Wang LJ, et al.** "Predicting and characterizing a cancer dependency map of tumors with deep learning." *Science Advances* 2021. **PMID 34417181**. — **DeepDEP**. Unsupervised pretraining on unlabeled tumor genomics; applied to ~8,000 TCGA tumors; "first pan-cancer synthetic dependency map with clinical relevance."
- **Shi X, Gekas C, Verduzco D, ... Flister MJ, Dezso Z.** "Building a translational cancer dependency map for The Cancer Genome Atlas." *Nature Cancer* 2024. **PMID 39009815**, DOI 10.1038/s43018-024-00789-y. — ML dependency maps for patient tumors; predicts drug response and outcome; PAPSS1/PAPSS2 and CNOT7/CNOT78 synthetic lethals validated in vitro and in vivo; also maps gene tolerability in healthy tissue for therapeutic-window prioritization; web app released.
- **Shi Y, Xu W, Hu P.** "Deep Unsupervised Domain Adaptation for Translating Cancer Dependency Maps From Cell Lines to Breast Cancer Tumor Genomics." *Genetic Epidemiology* 2026. **PMID 42339999**.
- **Yu G, Gong Y, Fan B.** "Context-aware gene dependency modeling via graph attention networks for precision oncology" (GATDep). *J Transl Med* 2025. **PMID 41291726**. Sample-specific gene dependency from transcriptomics.
- **Lin CH, Lichtarge O.** BioVNN. *Bioinformatics* 2021. **PMID 34042953**.
- **Gross B, Dauvin A, Cabeli V, et al.** "Robust evaluation of deep learning-based representation methods for survival and gene essentiality prediction on bulk RNA-seq data." *Sci Rep* 2024. **PMID 39048590**. Notably reports that "baseline methods achieve comparable or superior performance compared to more complex models" on the survival task.

**Implication for T2:** the "place a tumor on a dependency manifold" idea is ~5 years old in the RNA modality (DeepDEP, 2021) and has a Nature Cancer flagship (2024). **T2's contribution cannot be "we do it for patients" — it must be "we do it from pixels, and we measure how much is real."** The Nature Cancer 2024 paper in particular makes falsification test #4 (WSI-only vs RNA-only vs WSI+RNA) not optional but *obligatory*: a reviewer will ask why an image is needed when Shi et al. already did this from omics.

**(b) The "beat the mean / linear baseline" move — published, in single-cell, three times over. Thesis correctly concedes this.**

- **Ahlmann-Eltze C, Huber W, Anders S.** "Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines." *Nature Methods* 2025. **PMID 40759747**, DOI 10.1038/s41592-025-02772-6. (bioRxiv preprint 10.1101/2024.09.16.613342.)
- **Wenteler A, Occhetta M, Branson N, et al.** "PertEval-scFM: Benchmarking Single-Cell Foundation Models for Perturbation Effect Prediction." *bioRxiv* 2024. DOI 10.1101/2024.10.02.616248.
- **Wong DR, Hill AS, Moccia R.** "Simple controls exceed best deep learning algorithms and reveal foundation model effectiveness for predicting genetic perturbations." *Bioinformatics* 2025. **PMID 40407144**, DOI 10.1093/bioinformatics/btaf317.

These three own the *rhetorical move*. None of them touches images, tissue, lineage decomposition, or patient-level residuals. The move is a template, not a claim — templates are reusable, but only if the new instantiation carries its own finding.

**(c) The decomposition itself — I could not find it in ANY modality in the form T2 proposes.**

Searched for: within-lineage vs lineage-mean dependency evaluation; "lineage mean" / "lineage-matched" + dependency + CRISPR (Europe PMC, 25 results, 0 on-target — closest was PRODE, *Genome Biology* 2025, **PMID 40022167**, which is neighborhood-informed essentiality scoring, not an evaluation decomposition); DepMap + lineage + baseline/confounder (0 on-target); arXiv `all:"dependency prediction" AND all:"cancer" AND all:"baseline"` (0 results).

**No published three-way decomposition of predicted dependency into cancer-type-mean / covariate-explained / patient-residual, reported per perturbation, exists in the dependency-prediction literature.** This is the genuinely open slot and it is the thesis's stated deliverable. Good.

**(d) The WSI→molecular "it's all cohort structure" finding — also not published in the form the project has it.**

Searched Europe PMC and arXiv for within-cancer vs across-cancer stratified reporting of WSI→expression performance. Retrieved SEQUOIA (*Nat Commun* 2024), Path2Omics (*Cancer Research* 2025, PMID 41166699), BiSCALE (*Adv Sci* 2026, PMID 41736695), HiST (*Adv Sci* 2026, PMID 41487073), STimage (*Nat Commun* 2026, PMID 41545411), sCellST (*Nat Commun* 2026), HistoGWAS (*Genome Biology* 2026, PMID 41918129), OmiCLIP/Loki (*Nat Methods* 2025), MISO (*Nat Commun* 2025), the *J Ovarian Res* 2026 genome-wide r=0.36 paper (PMID 42277915). **None of them separately reports within-cancer-type control-adjusted performance vs. across-type performance.** Path2Omics (30 TCGA cancer types) and HiST (five cancer types) both aggregate. The project's established finding (i) is, as far as this sweep can tell, unpublished — which is an asset, not a liability, but it means T2 rests on an unpublished premise that reviewers will demand be shown.

**(e) Closest methodological threat to falsification test #4 (modality necessity) — a real, recent, and adjacent preprint.**

> **Richter T, Zimmermann E, Hall J, Theis FJ, Raghavan S, Winter PS, Amini AP, Crawford L.** "Beyond alignment: synergistic integration is required for multimodal cell foundation models." *bioRxiv* 2026. DOI **10.64898/2026.02.23.707420**, Europe PMC ID **PPR1221797**.

Introduces the **Synergistic Information Score (SIS)**, grounded in **partial information decomposition**, to test when cross-modal fusion adds genuine complementary information vs. redundancy. Benchmarks ten fusion approaches on spatial transcriptomics. Key finding: "standard alignment-based fusion objectives on frozen encoders inherently collapse to detecting linear redundancies"; **"tasks dominated by linear redundancies are sufficiently served by unimodal baselines."**

This is the single most dangerous paper found. It does not touch dependency or perturbation prediction, but it **pre-empts the framing and the formalism** of T2's modality-necessity test. If T2 runs WSI-only vs RNA-only vs WSI+RNA as a bespoke ablation without engaging PID/SIS, a reviewer who knows this preprint will say the analysis was done better elsewhere. Mitigation: adopt SIS as the modality-necessity metric and cite it. That converts a threat into a rigor credential.

**(f) Batch/alignment critique in the ST modality — partially occupies the T1 side-finding.**

> **Gindra RH, Palla G, Nguyen M, Wagner SJ, Tran M, Theis FJ, Saur D, Crawford L, Peng T.** "HESCAPE: A Large-Scale Benchmark of Cross-Modal Learning for Histology and Gene Expression in Spatial Transcriptomics." arXiv, Aug 2025 (v2 2025-08-27).

Pan-organ, 6 gene panels, 54 donors. Finds contrastive pretraining improves mutation classification but **degrades** expression prediction, and identifies **batch effects as the key interfering factor**. This is adjacent to — but not identical with — the T1 side-finding. HESCAPE does not do the mean-spot-subtracted, held-out-slide, like-for-like comparison against bulk. So the side-finding's "COULD-NOT-VERIFY / nobody has done it" status **holds after an independent sweep** (I searched arXiv `all:"mean expression baseline"` → 0 results; `all:"spatial gene expression" AND all:"baseline" AND all:"evaluation"` → 1 irrelevant result; Europe PMC ST+"critical assessment"/"pitfalls"/"overestimate" → 0 on-target). Two independent searches failing to find it is weak positive evidence of absence, not proof. Note also the *J Ovarian Res* 2026 datapoint: genome-wide mean Pearson **r = 0.36** for bulk-level H&E→expression in a single cancer type — a published number worth having when arguing about estimands.

**(g) Cell Painting / morphological profiling → perturbation response — a real conceptual precedent in a different tissue regime.**

- **Chandrasekaran SN, Alix E, Arevalo J, et al.** "Morphological map of under- and overexpression of genes in human cells." *Nature Methods* 2025. **PMID 40775081**. JUMP-CP: 15,243 genes perturbed in U-2 OS via ORF and CRISPR KO; morphological profiles reveal gene clusters and functional relationships.
- **Peng R, Liu Z, Gao Y, Wang J.** "A generative framework for predicting cellular morphological and transcriptomic perturbation responses" (MultiVCDiff). *Cell Reports Methods* 2026. **PMID 42167225**. Jointly predicts morphology images and expression from perturbations; 1,118 chemical + 130 genetic perturbations; zero-shot screening.
- **Seal S, Dee W, Shah A, et al.** "Counting cells can accurately predict small-molecule bioactivity benchmarks." *Nature Communications* 2026. **PMID 41651839**. — **Directly relevant to T2's positive sub-claim.** This is the Cell-Painting analogue of "the only thing your image predicts is proliferation": it shows a cell-count baseline is a serious competitor on bioactivity benchmarks and recommends filtering benchmarks to properly assess phenotypic signatures. If T2's surviving positive is "proliferation-coupled core essentiality," Seal et al. is the paper that already made that argument one modality over, and it must be cited as the intellectual precedent — including its conclusion that such benchmarks need filtering, not celebrating.
- **Kalinin AA, Arevalo J, Serrano E, et al.** *Nat Commun* 2025, **PMID 40467541** (mAP framework for profile strength) — relevant evaluation machinery.

**Implication:** morphology→perturbation-manifold is an established idea in cell culture (JUMP-CP, Cell Painting). T2 must not present "images encode perturbation-relevant state" as novel. What is novel is doing it on **patient H&E**, where the perturbation was never applied and the manifold must be transported from cell lines.

**(h) Miscellaneous checked and cleared**

Counterfactual/synthetic-modality work in pathology exists (*PLoS Comput Biol* 2026 PMID 41990107, cross-modal diffusion + counterfactual analysis on TCGA; *IEEE JBHI* 2025 PMID 40042950, counterfactual histology-genomic risk stratification) but none is interventional in the perturbation sense. HistoGWAS (*Genome Biology* 2026, PMID 41918129) uses variance-component methods on histology across 11 tissues — worth reading for the variance-decomposition machinery, but it is a GWAS framework, not a perturbation benchmark.

---

## 2. Why the kill fails

Three independent conditions all hold:

1. **The object is unoccupied.** No paper predicts genome-scale genetic dependency from patient H&E. Confirmed across PubMed, Europe PMC (incl. preprints), and arXiv with eight distinct query formulations.
2. **The deliverable is unoccupied.** No paper — in any modality — decomposes predicted dependency into cancer-type-mean / lineage-explained / patient-residual and reports the residual per perturbation. The Ahlmann-Eltze/PertEval/Wong template exists; the lineage-stratified instantiation does not.
3. **The null-survivability is real and is what protects it.** Because Shi et al. 2024 (Nature Cancer) and Chiu et al. 2021 already established the *positive* framing from omics, a well-powered image-modality **null** is a genuine addition to the record rather than a non-result. Nothing in the literature pre-empts that null.

Coladan (Wang et al., *Genome Medicine* 2026) owns the *phrase* and a spot-level in-silico perturbation demonstration, but it perturbs a CLS embedding to generate hypotheses on two prostate sections. It has no lineage decomposition, no patient-residual estimand, no held-out-perturbation-class evaluation, no genome-scale CRISPR ground truth, and no cross-cancer leakage control. It is a name collision, not a scoop.

---

## 3. Mandatory narrowings (the price of survival)

1. **Never write "we predict patient-specific dependencies."** Shi et al. 2024 and Chiu et al. 2021 own that sentence in omics. Write "we measure what fraction of an image-derived dependency prediction is patient-specific."
2. **Falsification test #4 must use partial information decomposition / SIS (Richter et al. 2026, DOI 10.64898/2026.02.23.707420), not a bespoke ablation.** Otherwise the modality-necessity analysis is a weaker version of an existing preprint.
3. **The proliferation/core-essentiality positive must be framed against Seal et al. 2026 (PMID 41651839).** "The image predicts the growth-rate-coupled part" is exactly the cell-count-baseline result, one modality over. Present it as cross-modality replication of a known confound, not as discovery.
4. **Add an RNA-only arm using the Shi/Chiu protocol as the reference ceiling**, not as a strawman. The honest headline may well be "H&E adds nothing over RNA for dependency placement" — which is publishable given the assets, but only if the RNA arm is a faithful reimplementation of the state of the art.
5. **Fix the citation.** Wang Z et al., *Genome Medicine* 2026, PMID 42449400 — not "Coladan et al."

## 4. Title that survives this sweep

> *"How much of an image-derived cancer dependency prediction is about the patient? A lineage-decomposed benchmark on 6,192 tumors and 11,258 genome-scale perturbations."*

Every noun in that title was checked against the literature above and none of them is claimed.

---

## 5. Search audit trail

- **WebSearch:** not used (exhausted per constraint).
- **Semantic Scholar API:** attempted 3×, HTTP 429 every time. **No S2 result is cited here.**
- **arXiv API:** 9 queries. Note: `abs:"phrase" AND abs:"phrase"` returns 0 spuriously; `all:` form works. All arXiv results above used the `all:` form.
- **NCBI E-utilities:** 8 esearch + 3 esummary + 1 efetch (full XML for PMID 42449400).
- **Europe PMC REST:** 11 queries, `resultType=core`, includes bioRxiv/medRxiv preprints.
- **COULD-NOT-VERIFY:** (i) the non-existence of a like-for-like mean-spot-subtracted held-out-slide ST-vs-bulk comparison — absence of evidence only, across two independent search strategies; (ii) whether Shi et al. 2024 or Chiu et al. 2021 internally report any lineage-stratified control analysis (abstracts do not say; full texts not retrieved in this sweep — **this must be checked before writing, as it is the one thing that could still kill the decomposition claim**).
