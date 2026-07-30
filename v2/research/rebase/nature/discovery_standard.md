# What Counts as a Credible Biological Discovery Without Wet-Lab Validation

**Scope:** computational-only (or near-computational-only) papers, 2022–2026, whose *biological* discovery claims were accepted in high-tier venues.
**Purpose:** derive the concrete evidentiary bar Morpheus must clear, given a hard no-wet-lab constraint.

---

## 0. Verification protocol for this document

- Every citation below was resolved through **PubMed E-utilities** (`esearch` → `esummary` → `efetch`) this session. PMID, journal, year and DOI are machine-returned, not recalled.
- **Semantic Scholar API returned HTTP 429 on every attempt** this session; nothing here depends on it.
- Where an abstract could not be retrieved in full, the entry is marked **PARTIAL-VERIFY** and the claim about its validation design is stated conservatively.
- Nothing in this file is cited from memory. If it is not in the ledger, it is not asserted.

### Citation ledger (all verified)

| # | Paper | First author | Venue | Year | DOI | PMID |
|---|---|---|---|---|---|---|
| 1 | Accurate proteome-wide missense variant effect prediction with AlphaMissense | Cheng J | Science | 2023 | 10.1126/science.adg7492 | 37733863 |
| 2 | Genome-wide prediction of disease variant effects with a deep protein language model | Brandes N | Nature Genetics | 2023 | 10.1038/s41588-023-01465-0 | 37563329 |
| 3 | Identifying disease-critical cell types and cellular processes by integrating scRNA-seq and human genetics (sc-linker) | Jagadeesh KA | Nature Genetics | 2022 | 10.1038/s41588-022-01187-9 | 36175791 |
| 4 | Plasma proteomic associations with genetics and health in the UK Biobank (UKB-PPP) | Sun BB | Nature | 2023 | 10.1038/s41586-023-06592-6 | 37794186 |
| 5 | Rare variant associations with plasma protein levels in the UK Biobank | Dhindsa RS | Nature | 2023 | 10.1038/s41586-023-06547-x | 37794183 |
| 6 | Clustering predicted structures at the scale of the known protein universe | Barrio-Hernandez I | Nature | 2023 | 10.1038/s41586-023-06510-w | 37704730 |
| 7 | Machine learning-based penetrance of genetic variants | Forrest IS | Science | 2025 | 10.1126/science.adm7066 | 40875860 |
| 8 | Advancing regulatory variant effect prediction with AlphaGenome | Avsec Ž | Nature | 2026 | 10.1038/s41586-025-10014-0 | 41606153 |
| 9 | A plasma proteomics-based candidate biomarker panel predictive of ALS | Chia R | Nature Medicine | 2025 | 10.1038/s41591-025-03890-6 | 40830661 |
| 10 | Predicting transcriptional outcomes of novel multigene perturbations with GEARS | Roohani Y | Nature Biotechnology | 2024 | 10.1038/s41587-023-01905-6 | 37592036 |
| 11 | Pan-cancer integrative histology-genomic analysis via multimodal deep learning | Chen RJ | Cancer Cell | 2022 | 10.1016/j.ccell.2022.07.004 | 35944502 |
| H1 | Uncovering new families and folds in the natural protein universe | Durairaj J | Nature | 2023 | 10.1038/s41586-023-06622-3 | 37704037 |
| H2 | Using artificial intelligence to document the hidden RNA virosphere (LucaProt) | Hou X | Cell | 2024 | 10.1016/j.cell.2024.09.027 | 39389057 |
| C1 | The impact of site-specific digital histology signatures on deep learning model accuracy and bias | Howard FM | Nature Communications | 2021 | 10.1038/s41467-021-24698-1 | 34285218 |
| C2 | Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines | Ahlmann-Eltze C | Nature Methods | 2025 | 10.1038/s41592-025-02772-6 | 40759747 |
| C3 | A comparison of computational methods for expression forecasting | Kernfeld E | Genome Biology | 2025 | 10.1186/s13059-025-03840-y | 41250104 |
| X1 | Transfer learning enables predictions in network biology (Geneformer) | Theodoris CV | Nature | 2023 | 10.1038/s41586-023-06139-9 | 37258680 |

`H` = hybrid (see §2, important negative result for the no-wet-lab premise). `C` = counter-evidence. `X` = contrast case.

---

## 1. The short answer

A computational-only biological discovery is believed when the authors have **manufactured an independent test that the model could not have optimised against, and that a confound-only explanation would fail.** Wet-lab validation is not special because it is wet; it is special because it is *causal, orthogonal in measurement, and generated after the prediction*. Anything that reproduces two or three of those three properties can substitute.

Across the eleven accepted computational-only papers below, **every single one used at least two independent lines of the following six**, and none was accepted on within-dataset performance alone:

| Substitute | Reproduces | Strength |
|---|---|---|
| **S1. Held-out / external cohort replication** | post-hoc-ness | medium |
| **S2. Orthogonal data modality** (different assay measuring the same latent) | measurement independence | high |
| **S3. Known-biology recovery + novel extension** (positive controls that were never trained on) | calibration of the discovery channel | medium-high |
| **S4. Perturbation / dependency-atlas concordance** (Perturb-seq, DepMap, DMS) | *causality* | **highest — this is the true wet-lab surrogate** |
| **S5. Clinical or organismal outcome association** | biological consequence at the phenotype level | high |
| **S6. Prospective/temporal separation** (pre-symptomatic samples, later-collected data) | post-hoc-ness, strongly | high |

**The single most load-bearing observation:** S4 is the only substitute that carries causal content. Papers that make *mechanistic* claims without S4 (or wet lab) are graded down to "resource" or "association" papers regardless of venue. Papers that make *predictive/annotation* claims can survive on S1+S2+S3.

---

## 2. Negative result that reshapes the premise (read first)

Two of the flagship "AI-made-a-discovery" papers I set out to use as pure-computational exemplars **are not pure**. Verbatim from their abstracts:

- **Durairaj 2023 (Nature, H1):** after uncovering the β-flower fold computationally, they "experimentally demonstrated that one of these belongs to a new superfamily of translation-targeting toxin-antitoxin systems, TumE-TumA."
- **Hou 2024 (Cell, H2):** "A subset of these novel RNA viruses was confirmed by RT-PCR and RNA/DNA sequencing."
- **Theodoris 2023 (Nature, X1):** Geneformer's *in silico* deletion predictions were carried into iPSC-cardiomyocyte experiments.

**Lesson:** in the *discovery-of-new-entities* genre (new folds, new viruses, new therapeutic targets), Nature/Cell/Science top tier has effectively required a token wet-lab anchor — one confirmed instance out of thousands predicted. The one anchor is doing enormous work: it converts "our classifier fires on these sequences" into "at least one of these is real."

**Consequence for Morpheus:** the no-wet-lab constraint largely **forecloses the "we discovered N new biological entities" genre**. It does not foreclose the genres in §3 that succeeded without any wet lab: *variant/effect annotation at scale*, *genetic-architecture mapping*, *cell-type–disease attribution*, *penetrance/risk quantification*, and *evaluation-paradigm papers that overturn a field's belief*. Morpheus should target one of those, not entity discovery.

---

## 3. The eleven computational-only cases

### 3.1 AlphaMissense — Cheng 2023, Science (PMID 37733863)
**Discovered:** pathogenicity classification for all ~71M possible human missense variants; 89% classified as likely benign or likely pathogenic; plus the emergent finding that **gene-average pathogenicity predicts cell essentiality**, identifying short essential genes that statistical approaches are underpowered to detect.
**Stood in for wet lab:**
- **S3 (dominant):** "achieves state-of-the-art results across a wide range of genetic and experimental benchmarks, **all without explicitly training on such data**." The refusal to train on the evaluation channel is the entire credibility argument.
- **S2:** deep mutational scanning (a real experimental modality, just other people's) as an orthogonal readout.
- **S3-extension:** the essentiality result is the novel extension — a *different biological quantity* recovered by a model trained on population frequency, which no confound in the training data explains.
**Why believed:** the discovery channel (structure + conservation) is demonstrably calibrated on things we already know, then produces a claim about something it was never shown.

### 3.2 ESM1b variant effects — Brandes 2023, Nature Genetics (PMID 37563329)
**Discovered:** ~450M missense effect predictions, and the specific biological claim that **~2 million variants are damaging only in specific protein isoforms** — i.e. isoform context materially changes variant interpretation.
**Stood in for wet lab:**
- **S2:** 28 deep mutational scan datasets (orthogonal experimental modality).
- **S1/S3:** ~150,000 ClinVar/HGMD variants, held out.
- **Generalisation stress test:** extended to in-frame indels and stop-gains — a domain shift the model was not fit to.
**Note:** the isoform claim is a *reinterpretation* claim, not an entity-discovery claim. That genre survives without wet lab.

### 3.3 sc-linker — Jagadeesh 2022, Nature Genetics (PMID 36175791)
**Discovered:** cell types and cellular *programs* through which GWAS variants act, including GABAergic neurons in major depressive disorder, a **disease-dependent M-cell program in ulcerative colitis**, and a **disease-specific complement cascade process in multiple sclerosis**; plus the structural claim that in autoimmune disease, epithelial programs are *disease-dependent only* (response, not initiation).
**Stood in for wet lab:**
- **S3 (textbook execution):** "The inferred disease enrichments **recapitulated known biology** and highlighted notable cell–disease relationships." Known-biology recovery is stated first, novel claims second. This ordering is the rhetorical signature of the genre.
- **S2:** integration of three independent data types — scRNA-seq, epigenomic SNP-to-gene maps, GWAS summary statistics — none of which shares a technical confound with the others.
- **Mechanistic asymmetry as evidence:** the healthy-vs-disease-dependent program contrast is an *internal control* that a generic confound would not produce.
**This is the closest structural template for a Morpheus discovery paper.**

### 3.4 / 3.5 UKB-PPP — Sun 2023 & Dhindsa 2023, Nature (PMIDs 37794186, 37794183)
**Discovered:** 14,287 primary pQTLs across 2,923 proteins, **81% previously undescribed**; long-range epistasis of ABO × FUT2 secretor status on GI-tissue-enriched proteins; extension of PCSK9 genetically-proxied effects to new endpoints; disentangling of COVID-19 susceptibility loci.
**Stood in for wet lab:**
- **S1:** ancestry-stratified replication (non-European pQTL mapping) — replication across a population that differs in LD structure kills a large class of artifacts.
- **S3:** "technical and biological validations" as an explicit, named section of the paper. PCSK9 is the positive control; the new endpoints are the extension.
- **S5:** disease endpoint association via genetic proxying (Mendelian-randomisation logic), which additionally supplies a *causal direction* from the randomisation of alleles.
- **Coherence constraint:** the ABO/FUT2 result is believable because it lands on tissue-enriched expression that was predicted independently.
**Note the substitute for causality here is not perturbation but *Mendelian randomisation* — germline genotype as a natural experiment.** Morpheus has no germline randomisation; it must get causality from Perturb-seq instead.

### 3.6 AFDB structural clustering — Barrio-Hernandez 2023, Nature (PMID 37704730)
**Discovered:** 2.30M non-singleton structural clusters, 31% unannotated; ~4% species-specific; remote structural homology including **human immune-related proteins with putative prokaryotic homologs**.
**Stood in for wet lab:**
- **S3:** evolutionary analysis showing most clusters are ancient — an independent, non-circular prior on which clusters are real.
- **G10 executed explicitly (rare and admirable):** the species-specific 4% is reported as "**lower-quality predictions or examples of de novo gene birth**." They name their own artifact hypothesis in the abstract and refuse to resolve it. This is why the rest of the paper is trusted.
**Lesson:** publicly conceding the fraction of your result that is probably artifact buys credibility for the fraction that is not.

### 3.7 ML penetrance — Forrest 2025, Science (PMID 40875860)
**Discovered:** quantitative penetrance for 1,648 rare variants in 31 autosomal-dominant genes; penetrance varies by variant class; refined interpretation of VUS and LoF variants via clinical trajectories over time.
**Stood in for wet lab (this is the most complete example in the ledger):**
- **S1:** models trained on 1,347,298 EHR participants, **"then applied them to an independent cohort with linked exome data."** Train and test are different cohorts *and* different data types.
- **S5:** "was associated with **clinical outcomes**."
- **S2/S4-lite:** "and **functional data**" — orthogonal assay concordance.
- **S6:** "delineating clinical trajectories **over time**" — temporal structure the model cannot fake.
- **Baseline comparison stated:** "Compared with conventional case-versus-control approaches, ML penetrance provided refined quantitative estimates." They name the incumbent and beat it.
**Score: 4 of 6 substitutes, plus an explicit baseline. This is what a no-wet-lab Science paper looks like in 2025.**

### 3.8 AlphaGenome — Avsec 2026, Nature (PMID 41606153) — **PARTIAL-VERIFY**
Metadata verified (Nature, 2026, DOI 10.1038/s41586-025-10014-0). Full abstract not retrieved; the verified fragment concerns the sequence-length/resolution trade-off in sequence-to-function models. **I will not assert its validation design.** Included because it establishes that regulatory-variant-effect prediction remains a live top-tier venue in 2026, and because a same-family bioRxiv item surfaced in the same search — "Evaluating sequence-to-function deep learning models for ancestry-stratified regulatory variant effect prediction using multi-ancestry blood eQTLs" (Sun X, bioRxiv 2026, DOI 10.64898/2026.06.22.730889, PMID 42395544) — indicating that **ancestry-stratified external eQTL replication is the current accepted validation currency in this subfield.**

### 3.9 ALS plasma proteomics — Chia 2025, Nature Medicine (PMID 40830661) — **PARTIAL-VERIFY**
Full abstract blocked by a quoting limit; the retrieved summary states: 33 differentially abundant plasma proteins in ALS vs controls, **replication in an independent cohort**, an ML model with high reported diagnostic accuracy, and analysis of **pre-symptomatic samples** indicating the process begins years before symptom onset.
**Stood in for wet lab:** **S1** (independent cohort) + **S6** (pre-symptomatic/temporal separation) + biological coherence (muscle, nerve, energy metabolism — the expected ALS axes, i.e. S3).
**The pre-symptomatic analysis is the key move:** it is a natural prospective experiment. Samples collected before the outcome existed cannot be confounded by the outcome.

### 3.10 GEARS — Roohani 2024, Nature Biotechnology (PMID 37592036)
**Discovered:** predicted transcriptional responses to multigene perturbations "consisting of genes that were **never experimentally perturbed**"; 40% higher precision than existing approaches across four genetic-interaction subtypes.
**Stood in for wet lab:**
- **S4:** held-out perturbations within a combinatorial Perturb-seq screen — perturbation-atlas concordance, the causal surrogate.
- **Recovery of a structured biological taxonomy** (four distinct GI subtypes), not just aggregate correlation.
**But see §4 — this claim did not fully survive.** GEARS is in the ledger as an accepted paper *and* as the object of the field's most important recent correction.

### 3.11 PORPOISE — Chen 2022, Cancer Cell (PMID 35944502) — **PARTIAL-VERIFY**
Verified premise: joint image-omic prognostic modelling across TCGA, motivated by the observation that "most prognostic models are either based on histology or genomics alone." Full abstract not retrieved; **I will not assert whether an external non-TCGA cohort was used.**
**Relevance to Morpheus:** this is the direct predecessor of the WSI+bulk-RNA setup, in a strong-but-not-Nature venue (Cancer Cell). The venue gap between this and §3.7 is informative: TCGA-internal cross-validation with interpretability gets Cancer Cell; independent-cohort + clinical-outcome + functional-concordance gets Science.

---

## 4. The counter-evidence — what these substitutes fail to catch

These three papers define the traps. Morpheus must be built to survive them, and can potentially *weaponise* them.

**C1 — Howard 2021, Nature Communications (PMID 34285218), "The impact of site-specific digital histology signatures on deep learning model accuracy and bias."**
TCGA slides carry submitting-site signatures. Models predicting molecular or outcome labels from H&E can be reading *where the slide was made*, because site correlates with patient ancestry, stage, and processing. **Any Morpheus WSI claim that is not site-stratified is dead on arrival with a competent reviewer.** This is the single most likely reason a WSI→molecular discovery claim gets rejected.

**C2 — Ahlmann-Eltze 2025, Nature Methods (PMID 40759747), "Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines."**
Directly negates the accepted-in-2024 GEARS-class claim. The perturbation-prediction field spent two years believing a result that a linear baseline dissolves.

**C3 — Kernfeld 2025, Genome Biology (PMID 41250104), "A comparison of computational methods for expression forecasting."**
Eleven large-scale perturbation datasets, unified engine: "**it is uncommon for expression forecasting methods to outperform simple baselines.**"

**Synthesis of C1–C3:** the dominant failure mode in computational biology 2022–2026 is **not fabrication, it is baseline omission plus confound omission.** Two of the three corrections above are literally "you did not run the linear model."

**Direct read on Morpheus's own established findings:**
- Finding (i) — WSI→molecular signal is ~+0.07 Pearson within-cancer, control-adjusted, and **METHOD-INVARIANT across every method including baselines** — is a C2/C3-shaped result *that Morpheus already discovered internally*. Method-invariance means the model is not the finding; the *data relationship* is the finding, and it is small. Reported as a method win, it dies. Reported as a **field-level correction with the invariance as the central evidence**, it is exactly the C2/C3 genre, which published in Nature Methods and Genome Biology in 2025.
- Finding (ii) — effective-rank collapse (~40 vs ~180), recovered +53 (2.1×, 3 seeds) with **no change in benchmark score** — is currently a *model-internal* observation, not a biological discovery. Under the checklist below it fails G0 (claim typing): "a representation has higher rank" is not a proposition about biology. It becomes one only if the recovered rank directions are shown to carry biological content that the collapsed representation loses — which requires S4 (Perturb-seq concordance) or S5 (outcome association) on the *recovered dimensions specifically*.

---

## 5. THE CHECKLIST

A computational discovery must satisfy **all of G0–G4** and **at least two of G5–G9**, with **G10–G12** as tier-raisers. This threshold is calibrated to the ledger: every accepted paper in §3 clears G0–G4 and ≥2 of G5–G9; Forrest 2025 clears nine of thirteen.

### Mandatory gates

**G0 — Claim typing.** State the discovery as a falsifiable proposition **about biology**, in one sentence, with no model noun in it. If the sentence requires naming your architecture, you have a methods result, not a discovery.
*Test:* would the sentence still be interesting if a competitor's model produced it?

**G1 — Confound-only null.** Build the best model you can that uses **only** the confounders (site, batch, sequencing platform, ancestry, sex, age, stage, tumour purity, cohort membership) and no biology. Report its number next to yours. The discovery is the residual.
*For Morpheus:* submitting-site (C1), cancer type, purity. Held-out-cancer splits control cross-cancer leakage but **do not** control within-cancer site.

**G2 — Baseline dominance with the baseline's number printed.** Linear/ridge model, nearest-neighbour on raw features, mean-of-training-set, and the field-standard incumbent. State the delta and its CI. C2 and C3 exist because this gate was skipped.

**G3 — Method-invariance disambiguation.** Run ≥3 structurally different methods. Then interpret honestly:
- Methods disagree, yours wins → the claim is *about the method*; you owe G2 hard.
- **All methods agree** → the claim is *about the data*, and must be re-typed as a data-level or field-level claim. Method-invariance is not a null result; it is a positive result about where the signal lives. It is also the strongest possible evidence that the signal is not an artifact of your architecture.

**G4 — Stability under nuisance variation.** The discovery survives ≥3 seeds, ≥2 preprocessing pipelines, and reasonable hyperparameter perturbation. Report the discovery's effect size distribution across these, not the best run.

### Independent-evidence gates (need ≥2)

**G5 — External cohort replication.** A cohort used for **zero** decisions: no tuning, no threshold selection, no early stopping, no "we also checked." Ideally differing in population, assay platform, and collection site. Preregister the analysis before touching it.
*Morpheus:* any open-access non-TCGA WSI+outcome cohort. This is the highest-value missing asset.

**G6 — Orthogonal modality corroboration.** An independent measurement of the same latent quantity, sharing no technical confound. Bulk RNA and WSI from the same block are *not* orthogonal (shared patient, shared site, shared processing). CPTAC proteomics against TCGA transcriptomics is closer to orthogonal; imputed-from-A vs measured-in-B is orthogonal.
*Morpheus:* CPTAC is inventoried and unwired. Wiring it is the cheapest available G6.

**G7 — Known-biology recovery + registered novel extension.** Two halves, both required. (a) Recover a set of positive controls **specified before the analysis**, at a stated sensitivity. (b) Emit ≥1 novel claim of the same type through the *same* channel, with an explicit statement of why the channel's calibration in (a) licenses (b). Recovery alone is a sanity check, not a discovery; novelty alone is unfalsifiable.

**G8 — Perturbational / causal-atlas concordance.** The strongest substitute and the only one with causal content. The prediction must be **directional and pre-specified**: "perturbing X moves module M in direction d." Evaluate against Perturb-seq / DepMap / GDSC held out from every stage of model fitting. Report the null distribution over random gene sets, not just the hit.
*Morpheus:* K562 GWPS (11,258 × 8,248), K562-essential (2,285), RPE1 (2,679), plus DepMap and GDSC. This is Morpheus's *single strongest and least-exploited* credibility asset — it is the closest thing to a wet lab that exists on disk. The Squires ~1-intervention/node regime supports single-node causal identifiability claims; the absence of doubles forecloses Varici-style interaction identifiability, and the paper must say so.

**G9 — Clinical / organismal outcome association.** Association with survival, progression, response, or diagnosis, adjusted for the standard prognostic covariates of that disease (stage, grade, age, treatment), with the covariate-only model's C-index printed alongside. Unadjusted survival curves are not evidence.

### Tier-raisers

**G10 — Named artifact hypotheses.** A section that enumerates every way the result could be fake, tests each, and **reports the fraction of the result that is probably artifact**. Barrio-Hernandez 2023 states in its abstract that ~4% of clusters may be "lower-quality predictions." Conceding the bad fraction is what makes the good fraction credible.

**G11 — Temporal or prospective separation.** Data collected before the outcome existed (pre-symptomatic samples: Chia 2025), or a model frozen before a new dataset was released. Removes the entire class of post-hoc objections.

**G12 — A dated, public, falsifiable prediction.** Name a specific experiment, in a specific system, whose result would kill the claim; commit to it in the paper. Costs nothing without a wet lab and converts a static claim into a standing bet. This is the honest substitute for "we validated it": *we did not validate it, here is exactly how you can refute it.*

### Disqualifiers (any one sinks the claim)

- Reporting only aggregate correlation/AUC when the claim is mechanistic.
- Interpreting attention/attribution maps as biology without a perturbational or orthogonal check.
- Any test-set contact during model selection.
- Novelty claimed by absence from a literature search rather than by positive evidence.
- Effect size within the range induced by seed variation (G4).
- "Recapitulates known biology" as the *entire* validation (that is G7a without G7b — a calibration check dressed as a finding).

### Fast scoring rubric

| Gates cleared | Verdict |
|---|---|
| G0–G4 only | Methods/benchmark paper. Not a discovery. |
| G0–G4 + one of G5–G9 | Association paper. Nature Communications / specialist tier. |
| G0–G4 + two of G5–G9 | Credible discovery. Nature Genetics / Cancer Cell tier. |
| G0–G4 + G8 + two others + G10 | Mechanistic discovery without wet lab. Nature / Science tier. |
| Anything, minus G1 or G2 | Rejected, or published and later corrected (C1, C2, C3). |

---

## 6. Application to Morpheus

**The no-wet-lab constraint rules out** the entity-discovery genre (new folds, new viruses, new targets), which §2 shows has required a token wet-lab anchor even at the top tier.

**It leaves open, in descending order of fit to assets on disk:**

1. **Field-level correction / new evaluation paradigm.** Finding (i) — a genuine within-cancer control-adjusted WSI→molecular signal of ~+0.07 Pearson that is *method-invariant across every method including baselines*, with ~46–49% of apparent cross-cancer performance attributable to cohort structure — is a C2/C3-class result with a cleaner confound decomposition than either. G1, G2, G3 are already satisfied; this is rare. **Adding G5 (one external non-TCGA cohort) and G9 (outcome association) converts it from Nature Methods-tier to a serious Nature-tier evaluation paper.** The published precedents (C2 in Nature Methods, C3 in Genome Biology, C1 in Nature Communications) show the genre travels.
2. **Perturbation-anchored causal claim (G8).** Genome-scale Perturb-seq is the only asset that supplies causal content. Any mechanistic Morpheus claim should be routed through it with a pre-specified direction and a random-gene-set null. State the Squires-yes / Varici-no identifiability boundary explicitly — reviewers reward that and punish its absence.
3. **Finding (ii) needs re-typing.** Effective-rank collapse and its +53 recovery currently fails G0. Route the recovered dimensions through G8 (do they predict Perturb-seq responses the collapsed representation cannot?) or G9 (do they carry outcome information?). "No change in benchmark score" is not a weakness here — it is the *point*, and it is the same shape as the method-invariance finding: **the benchmark is not measuring the biology.** Two independent results converging on "the benchmark is blind" is a stronger thesis than either alone.

**Cheapest credibility purchases, ranked by evidence gained per unit effort:** (a) wire CPTAC → G6; (b) acquire one open-access external WSI+outcome cohort and preregister → G5 + G11; (c) route any biological claim through Perturb-seq with a pre-specified direction → G8; (d) write the named-artifact section with a stated artifact fraction → G10; (e) add site-stratified confound-only nulls to every WSI result → hardens G1 against C1.
