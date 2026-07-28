# Lane b_batch_confound — Batch-effect / cohort-confound critiques in comp-path & multimodal genomics

Remit: Find critiques showing apparent performance in computational pathology / multimodal (image+omics) models is driven by site/batch/cohort/cancer-type shortcuts rather than genuine biology. Assess whether "molecular prediction from WSI is mostly cancer-type identity" and our specific MORPHEUS V2 findings are already reported, partially known, or novel.

## Queries run
1. Howard batch effects TCGA digital histology site-specific signatures deep learning
2. computational pathology models learn site batch shortcuts not biology critique
3. gene expression prediction from histology mostly cancer type identity confounder within-cancer
4. Dehkharghanian batch effect deep learning pan-cancer TCGA whole slide image
5. "Biased data, biased AI" Dehkharghanian predict acquisition site TCGA images findings
6. molecular subtype prediction histology confounded by cancer type random gene set control specificity
7. SEQUOIA gene expression prediction whole slide image cancer type baseline permutation control critique
8. predicting gene expression from histology performance driven by tissue type not cell-level baseline mean expression
9. spatial transcriptomics prediction histology models beaten by mean expression baseline critique benchmark HEST
10. pan-cancer histology gene expression prediction inflated correlation cancer type as feature within-cancer drop
11. multimodal image omics survival prediction confounded cancer type identity leakage critique 2024
12. "predicting the mean" gene expression histology benchmark correlation inflated cross-sample structure
13. histology foundation model gene expression prediction correlation collapses when tissue type removed patient-level confound HEST-benchmark critique
14. CLIP contrastive pathology molecular alignment representation collapse rank effective dimensionality VICReg variance
15. hallmark gene set pathway prediction from histology tissue-of-origin baseline confound within-cancer-type Pearson

## Sources
- **Howard et al. 2021, Nat Commun** (https://pmc.ncbi.nlm.nih.gov/articles/PMC8292530/) — "The impact of site-specific digital histology signatures on deep learning model accuracy and bias." The canonical site-preserved-split critique with hard numbers: mean AUROC drop 0.069 (up to 0.291) under preserved-site CV; 91.1% of features degrade; 35.7% no longer detectable. Directly analogous to our cross-cohort inflation of Pearson.
- **Dehkharghanian et al. 2023, Diagnostic Pathology** (https://pmc.ncbi.nlm.nih.gov/articles/PMC10189924/) — "Biased data, biased AI: deep networks predict the acquisition site of TCGA images." Shows a model trained for cancer-type classification (KimiaNet) reveals acquisition site at >86% accuracy; DenseNet 70%. Proves cancer-type classifiers exploit medically-irrelevant site patterns — the shortcut mechanism.
- **Jaume/Mahmood-lab benchmarking, 2025, Nat Commun** (https://pmc.ncbi.nlm.nih.gov/articles/PMC11814321/) — "Benchmarking the translational potential of spatial gene expression prediction from histology." Mean PCC ~0.28 at best; simpler CNN (DeepPT) beats transformer/GNN models; correlations inflated on HVG/SVG subsets vs all genes; no random-gene baseline; downstream survival C-index barely reaches RNA-seq baseline (0.55 vs 0.58). Establishes "simple baseline / inflated-metric" critique for WSI->expression.
- **Gindra et al. HESCAPE, 2025, ICCV-W / arXiv 2508.01490** (https://arxiv.org/html/2508.01490v1) — Large-scale cross-modal histology/gene-expression benchmark. Baseline image encoder WITHOUT gene-expression pretraining matches or beats aligned counterparts (Gigapath baseline PCC 0.338 vs DRVI-Gigapath 0.277). Batch effects in the expression modality skew image-encoder representations. Baseline PCC 0.338 is strikingly close to our global 0.348.
- **CHRep, 2026, arXiv 2604.21573** (https://arxiv.org/html/2604.21573v1) — Under leave-one-slide-out evaluation, existing models suffer "slide-level appearance shifts and regression-driven over-smoothing that suppress biologically meaningful variation." Names the exact mechanism we observe (regression to a low-D manifold suppressing within-sample variance).
- **Deep-learning perturbation prediction does NOT beat linear baselines, 2025, Nat Methods** (https://www.nature.com/articles/s41592-025-02772-6) — Five foundation + two DL models fail to beat deliberately simple baselines for perturbation-effect prediction. Genomics-wide precedent for "complex model = baseline once you control properly."
- **Dimensional collapse in contrastive learning** (https://arxiv.org/pdf/2403.18699 ; VICReg https://www.abhik.ai/papers/vicreg) — Documents dimensional collapse, feature suppression, and modality gap as generic contrastive failure modes; VICReg variance term meant to prevent collapse. Context for our rank-collapse mechanism.
- **PEaRL, 2025, arXiv 2510.03455** (https://arxiv.org/html/2510.03455) — Recent Hallmark-pathway-from-histology representation-learning method. Shows the Hallmark-from-WSI task is an active target, but does NOT report a random-gene-set / within-cancer confound decomposition.

## Findings
- Site/batch is a well-established confounder in comp-path. **Howard 2021**: preserved-site CV drops AUROC by mean 0.069 (max 0.291); 91% of predictable features degrade, 36% vanish; site remains predictable at >0.85 AUROC even after stain normalization. [Howard PMC8292530]
- Cancer-type classifiers demonstrably encode site, not just biology: **Dehkharghanian 2023** recovers acquisition site at >86% accuracy from a cancer-type-trained network's features. [PMC10189924]
- WSI->gene-expression correlations are low and inflated by variable-gene selection and lack random controls: **benchmarking paper** mean PCC ~0.28, correlations rise to ~0.45 only on HVGs, no white-noise/random-gene baseline, and clinical downstream barely beats RNA-seq. [PMC11814321]
- Cross-modal alignment does NOT add biology over a plain image encoder: **HESCAPE** baseline Gigapath PCC 0.338 > aligned 0.277; batch effects in expression corrupt the image encoder. The 0.338 baseline ~ our 0.348 global. [arXiv 2508.01490]
- The specific mechanism — regression/alignment to a low-D expression manifold that over-smooths and suppresses genuine within-sample variance — is named by **CHRep** ("regression-driven over-smoothing suppresses biologically meaningful variation") and generically by the **dimensional-collapse** literature. [arXiv 2604.21573; arXiv 2403.18699]
- "Complex model collapses to a simple baseline once confounds are controlled" is an established cross-domain pattern (**perturbation prediction, Nat Methods 2025**). [s41592-025-02772-6]
- No source found that (a) quantifies the fraction of WSI->molecular-prompting Pearson attributable to cross-cancer cohort structure (~46-49%), (b) uses a random-gene-set control to define residual within-cancer specificity (~+0.07 Pearson for every method including baseline), or (c) reports effective-rank collapse of a biology head (rank ~5-6/256) with a healthy identity head (~84) as the internal signature of this confound.

## Novelty verdict
**Partially known — the direction is well-trodden; our specific quantitative decomposition and the internal-mechanism signature appear novel.**

Reasoning:
- The *high-level claim* — "apparent WSI->molecular performance is largely cancer-type / site / cohort shortcut, not biology" — is NOT novel. It is the explicit thesis of Howard 2021 (site) and Dehkharghanian 2023 (cancer-type classifiers encode site), and is strongly implied by the benchmarking and HESCAPE papers (baselines match complex models; metrics inflated). Anyone citing our work must credit this prior art; framing it as a wholly new critique would be wrong.
- The *within-cancer + random-gene-set-adjusted specificity decomposition* appears novel in this precise form. Prior art uses preserved-SITE splits (Howard) or patient-level splits (benchmarks); we use within-CANCER-TYPE stratification AND a random-gene-set null to isolate residual specificity, yielding the sharp result that genuine, control-adjusted signal is only ~+0.07 Pearson and is essentially identical across all methods and the MLP-CLIP baseline. The benchmarking paper notably LACKS a random-gene baseline (it even flags that white-noise genes should show zero correlation but does not build the control), so our random-control-adjusted number fills an acknowledged gap.
- The ~46-49% "fraction of Pearson that is cross-cancer cohort structure" (global 0.348 -> within-cancer 0.188) is a novel quantitative attribution. The HESCAPE baseline of 0.338 independently corroborates our global 0.348, strengthening (not undermining) our claim, but no source decomposes the global number into cohort vs within-cancor components.
- The *internal mechanistic signature* — biology head effective rank collapsing to ~5-6/256 while the identity head stays at ~84, anchor residual ~0 (identity == frozen teacher), driven by neighbour-KL + supcon alignment to a ~50-D Hallmark manifold with an ineffective per-dimension VICReg variance floor — appears novel as applied here. The generic failure modes (dimensional collapse, over-smoothing) are known (CHRep, VICReg literature), but no source reports this dual-head rank asymmetry as the fingerprint of the cohort-confound in a WSI->molecular model. This is our strongest novel contribution.
- SigLIP marginally beating MLP-CLIP (+0.005 within-cancer, wins 62% of targets) is novel but weak/incremental and consistent with the benchmarks' finding that method choice barely matters once confounds dominate.

Net: publish as a *sharper, better-controlled quantification and mechanistic diagnosis* of an already-recognized confound, explicitly citing Howard, Dehkharghanian, the Nat Commun benchmark, and HESCAPE as prior art. Do not claim to have discovered that WSI->molecular is cancer-type-driven; claim the random-control-adjusted specificity floor and the dual-head rank-collapse signature.
