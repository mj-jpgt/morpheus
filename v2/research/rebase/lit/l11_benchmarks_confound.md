## Benchmarks & confound-aware evaluation

*Lane l11. Remit: WSI->molecular/omics benchmarks (HEST/HEST-bench, SpaRED/SpatialBenchmark), cohort/site/batch confounding in comp-path (Howard, Dawood/Buyer-Beware, batch-effect studies), confound-aware and generalization evaluation, ORBIT-like held-out-cancer protocols. Every entry maps to MORPHEUS axes A1-A5. The recurring thesis of this lane: **current comp-path headline metrics are structurally confounded, so any MORPHEUS novelty claim about a "promptable / identified / interventional" representation is only credible if it ships an evaluation that the field does not yet have.** This directly reinforces the session finding in REBASE_THESIS.md §2 that the biology head's rank collapse is invisible to the confounded benchmark score.*

---

### A. WSI -> molecular / omics benchmarks

**1. HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis** (Jaume, Doucet, Song, Lu, ... Mahmood; NeurIPS 2024 Datasets & Benchmarks, Spotlight) — arXiv:2406.16192
- *Takeaway:* The reference benchmark that turns "predict expression from morphology" into a standardized, foundation-model-comparable task across 9 tumor types.
- *Technical summary:* 1,229 spatial-transcriptomic profiles paired to WSIs from 153 cohorts (26 organs, 2 species, 367 cancer samples/25 types), yielding 2.1M expression-morphology pairs and 76M nuclei. HEST-Benchmark freezes a patch encoder and regresses the top-50 highly-variable genes per organ (Ridge/XGBoost on frozen features), scoring Pearson correlation. It is the de facto probe for whether a pathology FM's frozen embedding carries molecular signal, and rankings differ markedly from ImageNet-style classification leaderboards.
- *Plain-English:* A big, clean library of matched tissue-picture + gene-readout pairs, plus a standard test for "how much of the gene activity can a model read off the image."
- *Applicability:* A1/A4 — HEST is the ready-made external harness for MORPHEUS's frozen-trunk claim: if TumorStateV2's WSI+RNA latent is genuinely molecular, its frozen features should be competitive on HEST-Benchmark *without* task-specific retraining, which is exactly the A1 "new question = new query, not new model" promise. A4: HEST is where "encode RNA vs. retrieve it" can be measured head-to-head.
- *Novelty implication:* Pre-empts any MORPHEUS claim to be "the first WSI->molecular representation" — that ground is taken. MORPHEUS must instead claim *promptability/identifiability on top of* a HEST-competitive trunk, and should report HEST numbers to earn credibility.

**2. Completing Spatial Transcriptomics Data for Gene Expression Prediction Benchmarking (SpaRED / SpaCKLE)** (Ruiz, Cárdenas, Manrique, Vega, Mejia, Arbeláez; Medical Image Analysis vol.106, 2025) — arXiv:2505.02980
- *Takeaway:* A 26-dataset curated benchmark (SpaRED) plus a transformer that *completes* dropout-corrupted ST data before benchmarking, cutting MSE >82.5%.
- *Technical summary:* SpaRED standardizes 26 public ST datasets; SpaCKLE is a transformer gene-expression completion model that imputes missing/zero-inflated spots. Evaluating 8 SOTA histology->expression predictors on raw vs. completed data shows completion "substantially improves results across all models," implying much prior benchmarking was measuring noise, not method quality.
- *Plain-English:* Spatial gene data is full of holes; if you patch the holes fairly first, the leaderboard of image->gene predictors reshuffles.
- *Applicability:* A3/A4 — the confound here is *label noise in the target*, a distinct axis from site confounding. For MORPHEUS's emergence/elicitation benchmark (A3), it warns that a molecular-knowledge readout must control target-side dropout or it will mis-rank representations exactly as headline metrics do.
- *Novelty implication:* Reframes "measure molecular knowledge" as needing target-denoising; strengthens the case that MORPHEUS's contribution should be an *evaluation protocol*, not just a predictor.

**3. Benchmarking the translational potential of spatial gene expression prediction from histology** (Wang C., Chan A.S., Fu X., et al.; Nature Communications, 2025) — https://www.nature.com/articles/s41467-025-56618-y (PMC11814321)
- *Takeaway:* Eleven image->ST methods judged not on within-slide accuracy but on *cross-study generalizability and downstream survival* — where simple CNNs (DeepPT) beat fancy architectures and no method wins everywhere.
- *Technical summary:* Five SRT datasets + external TCGA-BRCA validation, scored across five hierarchical categories: within-image accuracy, cross-study generalizability, clinical/survival translational impact, usability, and compute. Headline within-image correlation is a poor predictor of cross-study or survival utility; architectural complexity does not track translational value.
- *Plain-English:* A model that predicts genes well on the same slide often fails on a new cohort or at predicting patient outcome — so "accuracy" is the wrong yardstick.
- *Applicability:* A1/A3 (guardrail) — this is the canonical "held-out cohort, not held-out spot" argument, and it is precisely MORPHEUS's cancer-held-out guardrail (§4). It supplies a ready five-tier evaluation template MORPHEUS can adopt so its promptable representation is judged on transfer, not memorization.
- *Novelty implication:* Strengthens MORPHEUS's guardrail framing but pre-empts "cross-study evaluation of image->expression" as itself novel; MORPHEUS must layer promptability/identifiability on top.

**4. HiST: A Hierarchical Sparse Transformer for Cross-Modal Spatial Transcriptomics Modeling** (Wu, Xu, Diao, Li, Wei, Andersson, Gui; ICML 2026) — arXiv:2606.14251
- *Takeaway:* A gigapixel-efficient ST predictor whose key design is a *slide-calibration token* that absorbs acquisition/batch variation as a first-class modeling target.
- *Technical summary:* Sparse window attention over the measured-gene lattice plus resolution-changing operators; a per-slide calibration token handles acquisition drift so cost scales with observed spots, not slide area. Improves multi-organ benchmark accuracy while cutting runtime/memory.
- *Plain-English:* Instead of pretending every slide was made the same way, it gives each slide a small "correction knob," then predicts genes cheaply on huge images.
- *Applicability:* A2/A4 — the calibration-token idea is a concrete mechanism for MORPHEUS to make batch/site an *addressable, factored-out slot* rather than a leak into the biology head; directly relevant to the covariance-term rank-recovery result in §2.
- *Novelty implication:* Partially pre-empts "model site as an explicit latent"; MORPHEUS should cite it and differentiate by making the site-slot *identified and promptable-around*, not just a nuisance token.

**5. A deep learning framework for efficient pathology image analysis (EAGLE)** (Neidlinger, Lenz, Foersch, +24; 2025) — arXiv:2502.13027
- *Takeaway:* A frozen two-FM pipeline (CHIEF tile-selection + Virchow2 features) evaluated on a broad 43-task battery spanning morphology, biomarker, treatment response, and prognosis.
- *Technical summary:* Selectively analyzes informative regions rather than all tiles, reaching up to +23% over patch-aggregation baselines with >99% compute reduction (~2.27s/slide) across 43 tasks / 9 cancers. Demonstrates that a *frozen*-encoder + light task head is competitive across a wide, heterogeneous task set.
- *Plain-English:* Pick the few tissue regions that matter, run a frozen model on them, and you match or beat slower methods on dozens of clinical tasks.
- *Applicability:* A1/A4 — a direct precedent for MORPHEUS's frozen-trunk-plus-many-tasks thesis, and a strong "how many tasks with one frozen backbone" baseline. A1 must beat "43 hard-coded heads" by *inferring* the task from NL, not enumerating heads.
- *Novelty implication:* Pre-empts "one frozen backbone, many tasks" as novel per se; sharpens that MORPHEUS's delta is NL task auto-detection/routing, not multi-task frozen probing.

**6. THItoGene: predicting spatial transcriptomics from histology** (Jia et al.; Briefings in Bioinformatics, 2024) — PMC10749789
- *Takeaway:* A representative image->ST method (hybrid CNN + Vision-Transformer + graph attention) of the kind repeatedly used as a benchmark baseline in HEST/SpaRED.
- *Technical summary:* Uses dynamic convolutional and capsule/graph modules to capture spot neighborhoods and predict spatial gene expression from H&E. Serves as one of the SOTA baselines whose within-slide performance the benchmarking papers (#2,#3) show does not transfer.
- *Plain-English:* A specialized image-to-gene predictor that looks good on its home dataset but is exactly the kind of model the benchmark critiques show fails to generalize.
- *Applicability:* A4 — a concrete "encode-only, no retrieval, no site-control" baseline against which MORPHEUS's encode-vs-retrieve study can quantify the gain from RAG context and site handling.
- *Novelty implication:* Neutral/baseline; useful as a contrast, not a threat.

---

### B. Cohort / site / batch confounding in computational pathology

**7. The impact of site-specific digital histology signatures on deep learning model accuracy and bias** (Howard F.M., Dolezal, Kochanny, ... Pearson; Nature Communications 12:4423, 2021) — https://www.nature.com/articles/s41467-021-24698-1
- *Takeaway:* The foundational confounding result: TCGA submitting-site is deep-learning-detectable and inflates/biases predictions of survival, mutations, stage — and even ethnicity — surviving color normalization.
- *Technical summary:* Across ~3,000 patients / 6 cancer subtypes, CNNs identify the tissue-submitting site from H&E despite stain normalization and augmentation; site correlates with clinical/genomic labels, so measured accuracy for survival, driver mutations, and stage is partly a site-detector in disguise. They quantify the image characteristics constituting the "site signature."
- *Plain-English:* Hospitals leave an invisible fingerprint on their slides; models "cheat" by reading the fingerprint, so their reported skill at predicting genes or survival is partly fake.
- *Applicability:* Guardrail/A3 — the single most important citation for MORPHEUS's confound-aware evaluation. Any A3 "emergent biological knowledge" claim must show the signal is not the site signature; mandates site-stratified / leave-site-out splits.
- *Novelty implication:* Reframes the whole lane: MORPHEUS cannot claim "reads molecular biology from morphology" without a Howard-style site-confound audit. Strengthens the guardrail; raises the evidentiary bar.

**8. Buyer Beware: confounding factors and biases abound when predicting omics-based biomarkers from histological images** (Dawood, et al.; bioRxiv 2024.06.23.600257, published Nature Biomedical Engineering, 2026) — https://www.biorxiv.org/content/10.1101/2024.06.23.600257v1
- *Takeaway:* Omics-biomarker-from-WSI models capture *co-dependencies among biomarkers* rather than isolating each, so single-biomarker AUCs are confounded aggregates.
- *Technical summary:* On n=8,221 patients, statistical co-dependence testing shows biomarker statuses are strongly interdependent; models trained per-biomarker actually predict the correlated bundle. Recommends co-dependence-aware variable definition, diverse cohorts, disentangling architectures, and stringent stratification testing.
- *Plain-English:* When a model says it predicts biomarker X, it's often just riding X's correlation with Y and Z — the "single-biomarker" score is borrowed skill.
- *Applicability:* A2/A3 — this is the identifiability argument in evaluation form: without *disentangled, per-programme addressable slots* (A2), a biomarker "prompt" returns a confounded bundle. Directly motivates MORPHEUS's per-pathway addressability as the fix, and a co-dependence-controlled elicitation benchmark (A3).
- *Novelty implication:* Strongly *strengthens* A2's rationale (identifiability improves prompt reliability) while pre-empting naive "predict biomarker by prompt" claims. MORPHEUS should test whether identified slots reduce the co-dependence leakage Dawood measures.

**9. Current Pathology Foundation Models are unrobust to Medical Center Differences** (de Jong, Marcus, Teuwen; 2025) — arXiv:2501.18055
- *Takeaway:* Introduces the *Robustness Index*; across 10 public pathology FMs, medical-center identity is encoded more strongly than tissue or cancer type, and only one model exceeds RI>1.
- *Technical summary:* RI = degree to which biological neighbors dominate same-center neighbors in embedding space. For 9/10 FMs RI<1 (center dominates biology); center-of-origin is predicted more accurately than cancer type, and cancer-type errors are systematically same-center confounders.
- *Plain-English:* Today's big pathology models organize their memory by "which hospital" more than "which cancer," so their errors cluster by hospital.
- *Applicability:* A2/guardrail — gives MORPHEUS a single scalar (RI) to report for its trunk; a promptable representation with RI<1 would make prompts route on site, not biology. A2 identifiability should be validated by an RI improvement.
- *Novelty implication:* Supplies a concrete confound metric MORPHEUS can *adopt and beat*; if MORPHEUS's identified latent doesn't raise RI, the identifiability claim is empty. Reframes success criterion.

**10. Do Histopathological Foundation Models Eliminate Batch Effects? A Comparative Study** (Kömen, Marienwald, Dippel, Hense; AIM-FM Workshop @ NeurIPS 2024) — arXiv:2411.05489
- *Takeaway:* Large-scale SSL pretraining does *not* remove hospital signatures; they persist in FM embeddings, dominate feature-space distances, and span multiple PCs.
- *Technical summary:* Comparative probing of pathology FM embeddings shows distinct hospital clusters remain after stain normalization; the batch signal is a leading source of embedding variance, so downstream heads can silently learn it.
- *Plain-English:* Even the newest, biggest pathology models still smell of the hospital that made the slide.
- *Applicability:* A4/guardrail — kills the assumption that "just use a big frozen FM" solves confounding; MORPHEUS's frozen-trunk plug-in (A4) inherits the batch signal and must actively factor it out.
- *Novelty implication:* Pre-empts "scale removes batch"; strengthens MORPHEUS's need for an explicit site-slot / IRM guardrail rather than relying on the backbone.

**11. Comparing Computational Pathology Foundation Models using Representational Similarity Analysis** (Mishra, Lotter; ML4H 2025) — arXiv:2509.15482
- *Takeaway:* Borrowing RSA from neuroscience, six pathology FMs show pronounced *slide-dependence but minimal disease-dependence* in their internal geometry.
- *Technical summary:* RSA over TCGA H&E patches: UNI2/Virchow2 most distinct, Prov-GigaPath most central; all models' representational geometry is organized by slide, weakly by disease. Stain normalization reduces (not removes) slide artifacts.
- *Plain-English:* Compare the "mental maps" of six models and they mostly agree on which slide a patch came from, not which disease it is.
- *Applicability:* A2/A3 — RSA is a candidate tool for MORPHEUS's *representation-quality* evaluation (the missing eval in §2): compare identified-slot geometry to a biology-derived reference matrix, not just downstream accuracy.
- *Novelty implication:* Offers a method (representational geometry vs. biological reference) that could become MORPHEUS's A3 emergence metric — a way to measure biological knowledge that "isn't just downstream accuracy." Directly answers the A3 novelty question with a borrowable technique.

**12. Pathology Foundation Models are Scanner Sensitive: Benchmark and Mitigation with Contrastive ScanGen Loss** (Carloni, Brattoli, Keum, Park, Lee, Ahn, Pereira; MedAGI @ MICCAI 2025, Oral) — arXiv:2507.22092
- *Takeaway:* Scanner hardware is a distinct confounder from site; a contrastive ScanGen loss at fine-tune time restores cross-scanner generalization on EGFR-mutation prediction.
- *Technical summary:* Multi-scanner benchmark quantifies scanner-driven performance swings in FM features; ScanGen adds a scanner-invariance contrastive term during task fine-tuning, improving cross-scanner transfer while preserving EGFR-mutation AUC.
- *Plain-English:* Two scanners of the same slide can flip a model's answer; a special training penalty makes it ignore the scanner.
- *Applicability:* A4/A5 — for MORPHEUS's molecular-prompt tasks (EGFR is a "prompt"), scanner invariance must be part of the frozen-trunk contract; ScanGen is a concrete recipe.
- *Novelty implication:* Neutral-to-strengthening; a mitigation MORPHEUS can incorporate but not claim.

**13. Hospital-Specific Bias in Patch-Based Pathology Models** (Zhang M.; 2025) — arXiv:2508.14779
- *Takeaway:* A lightweight adversarial adaptor strips hospital-domain information from latent features while preserving disease classification.
- *Technical summary:* Adversarial adaptor separates disease-relevant from institution-specific factors in the embedding; t-SNE and cross-hospital benchmarks show reduced hospital clustering at negligible disease-accuracy cost.
- *Plain-English:* Bolt on a small module that erases the hospital fingerprint but keeps the diagnosis.
- *Applicability:* A2/A4 — an adversarial route to the *factored* latent MORPHEUS wants; complements the covariance-term approach in §2. Suggests MORPHEUS could compare adversarial vs. covariance vs. IRM for site-invariance.
- *Novelty implication:* Pre-empts "adversarial site removal" as the novel bit; MORPHEUS's differentiator must be that removal is *promptable-around and identifiable*, not just invariant.

**14. Mitigating Batch Effects in Histopathology via Language-Mediated Robust Embedding Generation (GLMP)** (Zhang Y., Wu S., Zhang Z., ... Zhu, Wu, Zhang D.; 2026) — arXiv:2606.28697
- *Takeaway:* Route image patches through an *intermediate textual description* (via a multimodal LLM) before embedding, so the numeric vector is anchored to biology-language rather than pixel-level batch artifacts.
- *Technical summary:* GLMP generates a text description of histological features, then embeds that, claiming the first pathology model to use text as an intermediate batch-robust representation; improves cross-institution generalization by prioritizing describable biological signal over institution-specific pixels.
- *Plain-English:* Instead of turning the image straight into numbers, first say in words what's in it — words carry the biology, not the hospital's staining quirks.
- *Applicability:* A3/A4 — this is *NL<->biology grounding used as a confound defense*, exactly MORPHEUS's A3 territory. Evidence that a language bottleneck both grounds and de-confounds, informing whether MORPHEUS's NL interface should sit in the representation path (not just the output).
- *Novelty implication:* Partially pre-empts "language grounding improves robustness"; MORPHEUS must differentiate by making the language interface *task-promptable and identified*, and by measuring emergence, which GLMP does not.

**15. Hidden Variables in Deep Learning Digital Pathology and Their Potential to Cause Batch Effects** (Schmitt et al., 17 authors; Journal of Medical Internet Research, 2021) — PMC7886613
- *Takeaway:* CNNs readily learn non-biological metadata — scanner type (100% balanced acc), institution, patient age, slide-prep date (56%) — any of which can become a batch confounder.
- *Technical summary:* On melanoma WSIs from five institutes, networks recovered all four hidden variables at above-chance-to-perfect accuracy, demonstrating that arbitrary acquisition metadata is learnable and can correlate with labels.
- *Plain-English:* Models can guess the scanner, the hospital, even the slide's prep date from the image alone — so any of those can secretly drive predictions.
- *Applicability:* Guardrail — broadens the confounder list beyond "site" to scanner/date/age; MORPHEUS's confound audit should probe several, and its per-sample missingness handling (§4) intersects with acquisition metadata.
- *Novelty implication:* Strengthens guardrail scope; no direct novelty threat.

**16. Deep feature batch correction using ComBat for machine learning applications in computational pathology** (Murchan, Ó Broin, Baird, Sheils, Finn; Journal of Pathology Informatics, 2024) — PMC11470259
- *Takeaway:* ComBat on raw patch embeddings drops site-predictability from AUROC>0.95 to ~0.5 while *preserving* MSI predictability — a working post-hoc de-confounder.
- *Technical summary:* Applied to frozen deep features across multiple tissue-source sites; outperforms Macenko color normalization (which barely helps), and external-cohort validation confirms retained clinical signal (MSI) after harmonization.
- *Plain-English:* A statistical scrub of the model's features removes the hospital signal but keeps the cancer signal — and stain-color tricks don't.
- *Applicability:* A4/guardrail — a concrete baseline de-confounder MORPHEUS should compare against; establishes that feature-space correction beats image-space normalization, informing where in the pipeline MORPHEUS factors out site.
- *Novelty implication:* Baseline; MORPHEUS's in-model identified-slot approach must beat post-hoc ComBat to justify complexity.

**17. Examining Batch Effect in Histopathology as a Distributionally Robust Optimization Problem** (bioRxiv preprint, 2021) — https://www.biorxiv.org/content/10.1101/2021.09.14.460365
- *Takeaway:* Frames batch effect as distribution shift and applies group-DRO so the model optimizes worst-site rather than average performance.
- *Technical summary:* Treats each site/batch as a group and minimizes worst-group risk (DRO), improving cross-site robustness relative to ERM that silently overfits the majority batch.
- *Plain-English:* Instead of doing well on average across hospitals, force the model to do well on its *worst* hospital.
- *Applicability:* Guardrail/A5 — DRO is a principled alternative/complement to IRM in MORPHEUS's "environment-balanced IRM or drop it" guardrail; relevant to held-out-cancer worst-group evaluation.
- *Novelty implication:* Supports the guardrail; positions IRM-vs-DRO as a design choice MORPHEUS must justify, not a novelty.

**18. A Survey of Pathology Foundation Model: Progress and Future Directions** (Xiong, Chen, Sung; IJCAI 2025 Survey Track) — arXiv:2504.04045
- *Takeaway:* Taxonomy of pathology FMs by scope/pretraining/architecture and by evaluation tier (slide, patch, multimodal, biological), naming generalization and confounding as open challenges.
- *Technical summary:* Organizes the field's models and their evaluation tasks, flags data-vs-capacity scaling, end-to-end training, and robustness as unsolved, and situates biological-task evaluation as under-developed relative to classification.
- *Plain-English:* A map of all the big pathology models and how they're tested, pointing out that "does it know biology / does it generalize" is the weak spot.
- *Applicability:* A1/A3 — positions MORPHEUS's promptable + biological-emergence evaluation as filling a named gap; useful for the related-work framing.
- *Novelty implication:* Confirms the *evaluation gap* MORPHEUS targets is acknowledged-but-open — strengthens the "novel evaluation" pitch.

---

### C. Confound-aware & held-out generalization evaluation (ORBIT-like protocols)

**19. Monotherapy cancer drug-blind response prediction is limited to intraclass generalization** (Herbert, Chia, Jensen, Walther-Antonio; PLOS Computational Biology, 2026) — https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013232
- *Takeaway:* "Drug-blind" generalization only works *within a mechanism-of-action class*; the strongest predictor of performance is the drug itself for 97% of drugs.
- *Technical summary:* Elastic-net analysis over 1,641 models of drug-to-drug information transfer: predictive coefficients cluster strictly by MoA (RAF/MEK/ERK, PI3K/AKT/mTOR); training per-mechanism beats global multi-drug training with ~10x less data. Cross-mechanism transfer is minimal.
- *Plain-English:* A model can predict response to a new drug only if that drug works like ones it already knows; across different mechanisms, it can't.
- *Applicability:* A5 — the direct evidence bound for MORPHEUS's interventional-query axis: a frozen representation answering "what if we give drug D" will only generalize within-MoA unless it encodes mechanism. Sets the honest ceiling for A5 counterfactual claims and defines the held-out protocol (held-out drug/MoA).
- *Novelty implication:* Pre-empts over-strong A5 claims (drug-blind counterfactuals across mechanisms); reframes A5 success as *within-mechanism* interventional generalization measured against a mean/MoA baseline.

**20. Simple controls exceed best deep learning algorithms and reveal foundation model effectiveness for predicting genetic perturbations** (Wong, Hill, Moccia; Bioinformatics 41(6), 2025) — https://academic.oup.com/bioinformatics/article/41/6/btaf317/8142305
- *Takeaway:* A biology-informed *mean* baseline beats scGPT and GEARS at held-out perturbation prediction; fine-tuned scGPT is statistically indistinguishable from random-init.
- *Technical summary:* CRISPR-informed mean (0 for CRISPRi target, 2x for CRISPRa) outperforms GEARS (ΔPearson 0.08, p=9.3e-4) and scGPT (0.11, p=8.1e-6) across three datasets, with 4.7-213.9x better rank scores; ablating pretrained weights and attention barely changes scGPT.
- *Plain-English:* For predicting what happens when you knock out a gene, a dead-simple rule beats the fancy foundation models — and the models' "pretraining" was doing almost nothing.
- *Applicability:* A5/guardrail — the mandatory baseline discipline for MORPHEUS's interventional queries: any counterfactual head must beat a mechanism-informed mean, or the "causal geometry" claim is unearned. Mirrors the §2 finding that headline scores can be blind to representation quality.
- *Novelty implication:* Strongly *reframes* A5 — MORPHEUS must pre-register simple controls; a counterfactual query that doesn't beat the CRISPR-mean is not novel. High-value adversarial prior-art for the causal axis.

**21. Benchmarking foundation cell models for post-perturbation RNA-seq prediction** (BMC Genomics, 2025) — https://link.springer.com/article/10.1186/s12864-025-11600-2
- *Takeaway:* Independent confirmation that scGPT/scFoundation are beaten by a train-set mean baseline on held-out-perturbation RNA-seq.
- *Technical summary:* Benchmarks single-cell FMs against simple baselines under held-out-perturbation splits; the mean-of-training-examples predictor matches or exceeds the FMs, underscoring that evaluation must use held-out perturbations (not held-out cells) and strong baselines.
- *Plain-English:* Another team finds that averaging beats the big single-cell models at guessing perturbation effects.
- *Applicability:* A5/guardrail — reinforces #20; the held-out-perturbation (not held-out-cell) split is the ORBIT-like protocol MORPHEUS must use for A5.
- *Novelty implication:* Strengthens the baseline-discipline guardrail; further raises the bar for any generative/interventional claim.

**22. Fourier Asymmetric Attention on Domain Generalization for Pan-Cancer Drug Response Prediction (FourierDrug)** (Song, Bai, Liu; 2025) — arXiv:2502.04034
- *Takeaway:* Treat each cancer type as a domain and learn cancer-type-invariant features so drug-response prediction transfers to *unseen* cancer types with no target data.
- *Technical summary:* Fourier transform + asymmetric attention clusters sensitive samples compactly and disperses resistant ones in frequency space; trained only on in-vitro cell lines, it matches/exceeds SOTA on single-cell and patient-level unseen-cancer prediction.
- *Plain-English:* Learn what "drug-sensitive" looks like in a way that doesn't depend on the cancer type, so you can predict for cancers you never trained on.
- *Applicability:* A5/guardrail — a concrete held-out-cancer domain-generalization recipe and the "each cancer = a domain" framing that operationalizes MORPHEUS's cancer-held-out guardrail for interventional queries.
- *Novelty implication:* Pre-empts "cancer-type domain generalization for drug response"; MORPHEUS's differentiator is answering it as a *query on a frozen promptable representation*, not a bespoke DG model.

**23. Generalizable AI predicts immunotherapy outcomes across cancers and treatments (COMPASS)** (medRxiv 2025.05.01.25326820, 2025) — https://www.medrxiv.org/content/10.1101/2025.05.01.25326820v3
- *Takeaway:* A concept-bottleneck/pathway model reports cross-indication generalization — identifying responders in cancer types *excluded from training* at ~76.5% on held-out cohorts.
- *Technical summary:* COMPASS(-PFT) maps expression to interpretable pathway/concept features and predicts immunotherapy response; the marquee test is cross-indication (train on some cancers, predict on held-out indications), an ORBIT-style protocol.
- *Plain-English:* A model that explains itself via biological pathways can still spot immunotherapy responders in cancer types it never saw.
- *Applicability:* A2/A3/A5 — a *pathway-addressable, concept-bottleneck* precedent that couples interpretability (A2/A3) with held-out-cancer transfer (A5); strong template for MORPHEUS's identified-slot + cross-indication evaluation.
- *Novelty implication:* Partially pre-empts "pathway-slot representation generalizes across cancers"; MORPHEUS must add NL-promptability and interventional queries, and show identifiability (not just interpretable concepts).

**24. Learning and actioning general principles of cancer cell drug sensitivity** (Nature Communications, 2025) — https://www.nature.com/articles/s41467-025-56827-5
- *Takeaway:* Learns transferable, mechanism-level principles of drug sensitivity that generalize across cell lines/contexts and can be "actioned" for hypothesis generation.
- *Technical summary:* Models drug sensitivity to extract generalizable determinants rather than dataset-specific correlations, with held-out evaluation across cell-line contexts and downstream actionable predictions.
- *Plain-English:* Instead of memorizing which cell dies to which drug, learn the underlying rules that carry over to new cells.
- *Applicability:* A5 — supports the feasibility of *transferable* (not memorized) sensitivity structure that a MORPHEUS interventional query could exploit; a comparator for "principle-level vs. correlational" generalization.
- *Novelty implication:* Neutral-to-strengthening for A5's "beat correlational baselines" question; a benchmark/comparator, not a blocker.

**25. Beyond the Failures: Rethinking Foundation Models in Pathology** (Tizhoosh; 2025) — arXiv:2510.23807
- *Takeaway:* Argues pathology FMs fail structurally because *dense embeddings cannot represent the combinatorial richness of tissue* — a direct architectural critique aligned with MORPHEUS's slot/compositionality thesis.
- *Technical summary:* Contends low accuracy/instability stem from inherited natural-image SSL, patch design, and dense-vector bottlenecks fragile to noise; calls for domain-specific, compositional architectures over general ViT adaptations.
- *Plain-English:* One dense vector per tissue can't capture how many things are going on at once, which is why these models are shaky.
- *Applicability:* A2 — an authority-level argument *for* MORPHEUS's move away from a single fused latent toward identified, addressable, compositional slots (VSA binding thesis in §1).
- *Novelty implication:* *Strengthens* A2's premise (dense fusion is the wrong substrate) while warning MORPHEUS not to reintroduce a single dense vector; the retracted single-fused-latent (§2) is exactly this failure mode.

**26. A General-Purpose Self-Supervised Model for Computational Pathology (UNI)** (Chen R.J., Ding, Lu, ... Mahmood; Nature Medicine, 2024) — arXiv:2308.15474
- *Takeaway:* The widely used pathology FM whose frozen features anchor most benchmarks — and which the confounding studies (#9-#11) show is site/slide-organized.
- *Technical summary:* ViT-L pretrained (DINOv2) on ~100M patches / 100k+ WSIs; evaluated on 34 clinical tasks including morphological and some molecular-adjacent tasks, establishing the frozen-linear-probe evaluation convention this lane critiques.
- *Plain-English:* A foundational, heavily benchmarked pathology model — the backbone others test for hidden hospital bias.
- *Applicability:* A1/A4 — a candidate frozen trunk / baseline for MORPHEUS and the reference for the "linear-probe on frozen features" evaluation MORPHEUS must both use and transcend (probe -> prompt).
- *Novelty implication:* Baseline/context; defines the "hard-coded probe" evaluation paradigm A1 aims to replace with NL task auto-detection.

**27. PLUTO-4: Frontier Pathology Foundation Models** (2025) — arXiv:2511.02826
- *Takeaway:* A current frontier pathology FM release reporting broad task-suite evaluation, representative of the leaderboard MORPHEUS's frozen trunk will be measured against.
- *Technical summary:* Scales pretraining and reports patch/slide-level benchmark performance across a wide task battery; extends the frozen-feature evaluation convention to the frontier scale.
- *Plain-English:* One of the newest, biggest pathology models, tested on many tasks at once.
- *Applicability:* A1/A4 — a moving-target baseline: MORPHEUS's promptable-representation claim must be framed relative to frontier FMs, not superseded models, and its delta (NL routing, identifiability, interventional queries) must be orthogonal to raw scale.
- *Novelty implication:* Pre-empts "bigger frozen backbone" as the contribution; keeps MORPHEUS honest that scale alone is not its axis.

---

### Cross-cutting synthesis for MORPHEUS

- **The lane's unifying verdict:** headline comp-path metrics (HEST correlations, biomarker AUCs, perturbation scores) are pervasively confounded by site/scanner/batch (Howard #7, de Jong #9, Kömen #10, Schmitt #15), by biomarker co-dependence (Dawood #8), and by weak baselines (Wong #20, BMC #21). This is the empirical backbone of REBASE_THESIS §2's claim that the current MORPHEUS eval is "structurally blind to representation quality."
- **Adopt-and-beat metrics MORPHEUS should ship:** Robustness Index (#9), ComBat site-AUROC drop (#16), RSA-to-biological-reference geometry (#11), leave-site-out and leave-cancer-out splits (#3, #22, #23), CRISPR/MoA-mean baselines for interventional queries (#20, #19).
- **Where the lane pre-empts MORPHEUS:** frozen-trunk multi-task probing (#5, #26, #27), image->molecular benchmarking (#1, #2, #3), site-invariance via adversarial/ComBat/language routes (#13, #16, #14), pathway-concept cross-cancer transfer (#23). MORPHEUS's defensible delta is the *combination* — NL task auto-detection (A1) over *identified, co-dependence-disentangled slots* (A2) evaluated with a *confound-aware emergence benchmark* (A3) and *baseline-beating interventional queries* (A5).
