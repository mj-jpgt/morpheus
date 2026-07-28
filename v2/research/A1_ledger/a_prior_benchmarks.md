# Lane A1 — Prior WSI→molecular/pathway prediction benchmarks

Remit: find existing benchmarks/papers that predict gene expression / pathway / Hallmark / molecular
signatures from H&E WSIs, and determine whether our specific confound (~50% of the WSI→Hallmark
Pearson signal is cross-cancer cohort identity; random-gene-set control ~0.30; genuine within-cancer
random-adjusted specificity only ~+0.07 for every method incl. baseline) is already documented, along
with whether prior work uses within-cohort controls or random-gene-set null models.

## Queries run
- HEST-1k benchmark gene expression prediction H&E whole slide image cross-cancer confound
- predicting spatial transcriptomics from histology within-cohort control random gene set null model
- gene expression prediction histology mean expression baseline outperforms spatial variation correlation inflated
- TITAN CONCH slide foundation model molecular signature prediction pan-cancer cohort confound benchmark
- "within-cancer" OR "within-cohort" gene expression prediction correlation drops evaluation stratified by cancer type pathology foundation model
- molecular subtype prediction histology confounded by tumor type proxy AUROC decomposition pan-cancer vs single cancer
- systematic pan-cancer deep learning multi-omic biomarkers pathology images "pan-cancer" confounder tissue type removing correlation Kather
- comprehensive benchmark image foundation models predict gene expression whole slide images biorxiv 2026 baseline mean expression cohort
- predicting gene expression from histology confounded by slide-level batch effect correlation across spots not biology critique
- HEST-bench highly variable genes Pearson correlation foundation model UNI gene prediction results table
- Gindra Palla large-scale benchmark cross-modal histology gene expression spatial transcriptomics baseline expression predictor findings
- gene expression prediction histology "random genes" OR "housekeeping genes" achieve high correlation control negative baseline warning
- HEST benchmark evaluation Pearson computed across patients confound tissue heterogeneity criticism per-slide correlation

## Sources
- HEST-1k / HEST-Benchmark (Jaume et al., NeurIPS 2024) — https://github.com/mahmoodlab/hest and https://openreview.net/forum?id=mlhFJE7PKo — The canonical WSI→gene-expression benchmark. Crucially, it is structured as **nine separate per-cancer tasks** (top-50 HVGs, one organ/cancer at a time, patient-stratified CV, Pearson across patients). By construction it evaluates WITHIN a single cancer type, so it does not itself expose or quantify a cross-cancer confound — but it is the framework everyone else's confound would live inside.
- "Benchmarking the translational potential of spatial gene expression prediction from histology" (Nat Commun 2025, PMC11814321) — https://pmc.ncbi.nlm.nih.gov/articles/PMC11814321/ — 11 methods, leave-one-slide-out + TCGA external validation. Documents generalization failure and that predicted-GE survival models barely match/underperform bulk RNA-seq baselines (C-index ~0.52–0.58 vs baseline ~0.57–0.58). Explicitly notes it does NOT run mean-expression, random-gene, or permutation null controls — a gap.
- Systematic pan-cancer multi-omic biomarker study (Arslan, Kather et al., Commun Med 2024, PMC10942985) — https://pmc.ncbi.nlm.nih.gov/articles/PMC10942985/ — 12,093 models, 4,031 biomarkers, 32 cancer types. Median AUC 0.644. Reports per-cancer AUC spread (0.585–0.768) but does NOT isolate how much of pan-cancer performance is cancer-type-as-proxy.
- MDPI Genes 2026 "Pan-Cancer Prediction of Genomic Alterations from H&E in a Real-World Clinical Cohort" — https://www.mdpi.com/2073-4425/17/4/371 — Directly warns that "a substantial proportion of apparently high-performing models were capturing tumor type, histological subtype, or staging-related features as a proxy for mutation status." This is the closest qualitative statement to our confound.
- TITAN multimodal WSI foundation model (Ding/Mahmood et al., Nat Med 2025, PMC12618242) — https://pmc.ncbi.nlm.nih.gov/articles/PMC12618242/ — 39 molecular-classification tasks; reports beating a mean-pooling baseline, i.e. uses a trivial baseline as control, but does not decompose cross-cancer structure for molecular tasks.
- HESCAPE large-scale cross-modal benchmark (Gindra, Palla et al., ICCV-W 2025, arXiv:2508.01490) — https://arxiv.org/abs/2508.01490 — Finds contrastive cross-modal pretraining *degrades* direct gene-expression prediction vs a plain regression baseline, and names **batch effects** as the dominant interfering factor. Establishes that non-biological structure dominates alignment.
- Multimodal-contrastive / linearized-attention line (e.g. Nat Commun 2024 s41467-024-54182-5) — https://www.nature.com/articles/s41467-024-54182-5 — Representative of the standard control in this subfield: comparing correlation to a **random untrained model of the same architecture** (an architecture null, NOT a random-gene-set / cohort-identity null).

## Findings
- **The cross-cancer / tumor-type-as-proxy confound is qualitatively documented, but not quantified in our decomposed form.** The MDPI 2026 real-world cohort paper and the pan-cancer molecular-subtype literature (biorxiv 333914; PMC5932236 "molecular classes transcending lineage") explicitly warn that pan-cancer molecular predictors ride on tumor-type/lineage signal, and recommend "computationally subtracting the molecular differences between cancer types." Nobody reports a clean number like "≈46–49% of Pearson is cross-cancer structure, collapsing global 0.348 → within-cancer 0.188."
- **HEST-Bench sidesteps rather than measures the confound.** By running one cancer type per task, HEST implicitly acknowledges cross-cancer structure would inflate a pooled metric, but it never states the magnitude or reports the pooled-vs-stratified gap. So the field's default benchmark neither documents nor quantifies our effect.
- **Random-gene-set / cohort-identity null models are essentially absent.** The standard "control" in the WSI→expression subfield is a random *untrained network* of the same architecture, or a mean-expression baseline (TITAN). No source in this lane runs a **random-gene-set null** (predicting an arbitrary matched-size gene panel) to show a baseline method already scores ~0.30 Pearson. Our "random-control-adjusted specificity ~+0.07 for every method including the baseline" style of analysis is not present in prior art.
- **Within-cohort/within-cancer controls exist but are used for generalization, not for confound decomposition.** The translational-benchmark (PMC11814321) uses leave-one-slide-out and cross-subtype transfer to show models fail to generalize; the Kather pan-cancer study reports per-cancer AUCs. Neither frames within-cancer evaluation as a *control that removes a cross-cancer confound and quantifies how much signal is lost* — which is exactly our contribution.
- **The "non-biological structure dominates" intuition is established (batch/slide effects).** HESCAPE names batch effects as the primary interference; the batch-effect literature (arXiv:2503.07173 batch-agnostic encoder; BLEEP) says apparent performance may reflect batch-level rather than biological correlation. This is adjacent to our finding but is about *slide/technical* batch, not *cross-cancer cohort identity as ~half the metric*.
- **Representation-collapse of a biology head (effective rank ~5–6 of 256) is not a documented phenomenon in WSI→pathway benchmarks.** This lane found no prior report tying benchmark inflation to a collapsed regression head / VICReg per-dimension variance floor failing to prevent rank collapse. That mechanism appears entirely outside the prior-benchmark literature.
- **SigLIP-beats-CLIP on molecular prompting** has no direct prior-art comparison in these benchmarks; MedSigLIP appears as one encoder among several (biorxiv 2026.03.02.709012) but no source isolates a SigLIP-vs-MLP-CLIP molecular-prompting margin.

## Novelty verdict
**Appears novel in its quantified, decomposed form; partially known in qualitative spirit.**

Reasoning:
1. Prior art *qualitatively* warns that pan-cancer molecular/mutation prediction from H&E rides on tumor-type-as-proxy (MDPI 2026; pan-cancer subtype literature; Kather pan-cancer study), and that non-biological (batch/slide) structure dominates cross-modal alignment (HESCAPE). So the *existence* of a cross-cohort confound is known — our result is not claiming a brand-new phenomenon class.
2. However, the **specific, falsifiable decomposition** — that ≈46–49% of the WSI→Hallmark Pearson is cross-cancer cohort structure (global 0.348 → within-cancer 0.188), that a **random-gene-set null** already scores ~0.30–0.32 globally, and that the genuine within-cancer, random-adjusted specificity is only ~+0.07 for *every* method including the baseline — is **not reported anywhere in the surveyed benchmarks**. HEST avoids the confound by design without measuring it; translational and pan-cancer benchmarks report generalization gaps and per-cancer spreads but never subtract a random-gene null or attribute a percentage of the metric to cohort identity.
3. The **random-gene-set null model** as a specificity control is, per this lane, absent from the WSI→expression/pathway benchmark literature (which uses random-untrained-network or mean-expression baselines instead). This makes our control methodology itself a novel contribution to the benchmark's critique.
4. The **mechanistic link** (256-D biology head collapsing to rank ~5–6 against a ~50-D Hallmark manifold; per-dimension VICReg floor failing to prevent it) is orthogonal to and unaddressed by prior benchmark papers.

Net: treat the *phenomenon* (cross-cohort confound in H&E→molecular prediction) as **partially known / previously flagged qualitatively**, but treat our **quantitative ~50% decomposition + random-gene-set-null-adjusted within-cancer specificity (~+0.07 for all methods) + the head-collapse mechanism** as **appears novel** and publishable as the concrete, numeric confound audit the field has gestured at but not executed.
