# Decision-support methods (retrospective only)

Lane l10. Remit: uncertainty quantification (UQ), selective prediction / abstention, next-best-test / active feature acquisition (AFA), and prognosis / treatment-response prediction from multimodal cancer data. Framed strictly as **methods + retrospective evaluation**, not deployed CDS. Every entry is mapped to MORPHEUS rebase axes A1-A5.

Axis legend: A1 promptable unified representation + NL task auto-detection; A2 identified, pathway-addressable slots; A3 NL<->biology grounding + emergent-knowledge elicitation and its evaluation; A4 multimodal prompting (encode-vs-RAG, frozen-trunk plug-in); A5 interventional/causal-geometry queries.

---

### Uncertainty quantification, calibration, selective prediction

1. **Good Rankings, Wrong Probabilities: A Calibration Audit of Multimodal Cancer Survival Models** (S. Ghawami, arXiv preprint, 2026) — https://arxiv.org/abs/2604.04239
- Takeaway: multimodal survival models rank patients well (high C-index) yet their predicted survival *probabilities* are systematically miscalibrated.
- Technical summary: Audited 3 models on TCGA-BRCA and 11 architectures across 5 TCGA cancer types with fold-level calibration tests and Benjamini-Hochberg FDR correction. 166 of 290 fold-level tests rejected correct calibration at median event time; MCAT reached C-index 0.817 on GBMLGG yet failed calibration on all five folds. Gating-based fusion calibrated better than bilinear/concatenation, and post-hoc Platt scaling fixed miscalibration without hurting discrimination.
- Plain-English: A model can correctly say who is sicker than whom while being badly wrong about the actual survival odds it reports — and ranking metrics hide this.
- Applicability: A2, A4. MORPHEUS decision-support outputs need a calibration contract, not just a C-index; fusion mechanism (gating vs concatenation) is a calibration lever, informing how modality slots are combined. Argues for reporting calibration on any prompted prognosis task.
- Novelty implication: Pre-empts a naive MORPHEUS claim that multimodal fusion "improves prognosis." Reframes the win as calibrated, addressable probabilities — strengthening the case for identified slots + explicit UQ over raw ranking gains.

2. **Deep learning uncertainty quantification for clinical text classification** (A. Chen et al., *Journal of Biomedical Informatics*, 2023) — https://www.sciencedirect.com/science/article/pii/S1532046423002976 (PMC11467893)
- Takeaway: compares UQ methods (softmax, MC-dropout, ensembles, evidential) for selective clinical text classification and shows abstention lifts accuracy at a rejection cost.
- Technical summary: Benchmarks predictive-uncertainty estimators on clinical NLP tasks (mortality, ICD coding, outpatient classification) and uses them to drive selective prediction; reports the accuracy-vs-abstention trade-off curves. Ensemble/evidential approaches give better selective-classification behaviour than raw softmax confidence.
- Plain-English: Teaching a clinical text model to say "I'm not sure" on the hard cases raises its accuracy on the ones it does answer.
- Applicability: A1, A3. MORPHEUS auto-detects an NL task then answers; this supports attaching a per-task uncertainty head so the promptable system can abstain rather than hallucinate a label.
- Novelty implication: Supports (does not pre-empt) a MORPHEUS abstention capability; reframes selective prediction as a per-prompt property rather than a fixed-probe property.

3. **A review of uncertainty quantification in medical image analysis: probabilistic and non-probabilistic methods** (L. Huang et al., arXiv, 2023) — https://arxiv.org/abs/2310.06873
- Takeaway: taxonomy of UQ methods (Bayesian, ensembles, evidential, conformal, test-time augmentation) for medical imaging and their evaluation metrics.
- Technical summary: Surveys aleatoric vs epistemic decomposition, calibration metrics (ECE, reliability diagrams), and downstream uses (segmentation quality control, OOD detection). Emphasizes that UQ method choice must match the decision it supports.
- Plain-English: A field guide to the many ways an imaging model can express doubt, and which are trustworthy.
- Applicability: A4. Informs which modality-encoders MORPHEUS should equip with UQ vs treat as RAG context; a modality worth encoding is one whose uncertainty is decision-relevant.
- Novelty implication: Background/reframing; positions MORPHEUS's UQ choices within an established taxonomy so novelty must be the *promptable* delivery, not UQ per se.

4. **Trustworthy clinical AI solutions: a unified review of uncertainty quantification in deep learning models for medical image analysis** (B. Lambert et al., *Artificial Intelligence in Medicine*, 2024) — https://www.sciencedirect.com/science/article/pii/S0933365724000721
- Takeaway: unifies UQ methods around the clinical-trust use case (failure detection, referral, abstention).
- Technical summary: Reviews single-pass (evidential, deterministic-distance), sampling (MC-dropout, ensembles), and post-hoc (conformal, temperature scaling) UQ, mapped to clinical workflows like selective referral. Highlights that time-critical settings favour single-forward-pass UQ.
- Plain-English: If a hospital AI is going to be trusted, it must flag its own likely mistakes; this catalogs how.
- Applicability: A1, A4. Argues MORPHEUS's frozen-trunk plug-in should carry a cheap single-pass UQ so any prompted task inherits abstention without extra sampling cost.
- Novelty implication: Reframes UQ as a trust-delivery mechanism; MORPHEUS's contribution should be *task-agnostic* UQ under a single representation, not a new estimator.

5. **Uncertainty Quantification for Machine Learning in Healthcare: A Survey** (L.J. Kimmel et al., arXiv, 2025) — https://arxiv.org/abs/2505.02874
- Takeaway: broad survey linking UQ methods to healthcare decision types and evaluation protocols.
- Technical summary: Organizes UQ by source (data, model, distribution shift) and by task (classification, regression, survival), and stresses retrospective evaluation pitfalls (leakage, miscalibration under shift). Recommends decision-curve / net-benefit style evaluation over raw calibration alone.
- Plain-English: A map of how to measure and report a healthcare model's confidence so the numbers mean something clinically.
- Applicability: A1, A2. Supports per-slot UQ reporting; informs the retrospective-eval protocol MORPHEUS should adopt for decision-support claims.
- Novelty implication: Background; sets the evaluation bar MORPHEUS must clear.

6. **A Comprehensive Survey on Evidential Deep Learning and Its Applications** (J. Gao et al., arXiv, 2024) — https://arxiv.org/abs/2409.04720
- Takeaway: consolidates evidential deep learning (Dirichlet-over-classes) as a single-pass UQ paradigm spanning classification, regression, segmentation.
- Technical summary: Formalizes EDL replacing softmax with a Dirichlet whose concentration encodes evidence, giving aleatoric+epistemic uncertainty in one forward pass; surveys extensions to regression, meta/active learning, few-shot. Notes calibration and evidence-collapse failure modes.
- Plain-English: A way to make one network output both its answer and a principled measure of how much evidence it has, without running it many times.
- Applicability: A4, A5. Cheap single-pass uncertainty suits a frozen multimodal trunk; evidential heads could gate whether a modality is encoded or deferred to RAG, and could quantify confidence in counterfactual queries.
- Novelty implication: Provides an off-the-shelf UQ head; MORPHEUS novelty cannot be "evidential UQ" but can be its use for per-programme addressable confidence.

7. **Deep evidential fusion with uncertainty quantification and reliability learning for multimodal medical image segmentation** (L. Huang et al., arXiv/*Information Fusion*, 2023) — https://arxiv.org/abs/2309.05919
- Takeaway: fuses modalities by weighting each by its learned reliability/uncertainty (evidential Dempster-Shafer combination).
- Technical summary: Each modality produces evidential mass functions; a reliability coefficient per modality is learned and combined via Dempster's rule, so unreliable modalities are down-weighted at fusion time. Improves segmentation robustness when a modality is noisy or degraded.
- Plain-English: When combining scans, the model trusts each source only as much as it has earned, instead of averaging blindly.
- Applicability: A4. Directly informs MORPHEUS's encode-vs-RAG decision: reliability-weighted fusion is a principled criterion for how much a proteomics/CNV/SNV slot should influence the answer.
- Novelty implication: Pre-empts a claim of "novel reliability-aware fusion"; MORPHEUS should differentiate via promptable, identified per-pathway reliability rather than per-modality masses.

8. **Uncertainty-aware abstention in medical diagnosis based on medical texts** (H. Alaa et al., arXiv, 2025) — https://arxiv.org/abs/2502.18050
- Takeaway: selective-prediction framework for text-based diagnosis that abstains under high predictive uncertainty across binary, multi-label, and multi-class tasks.
- Technical summary: Models predictive uncertainty on MIMIC-III mortality, MIMIC-IV ICD-10 multi-label, and outpatient multi-class tasks; sets abstention thresholds to hit target accuracy while minimizing rejection rate. Shows accuracy/rejection frontiers across task types.
- Plain-English: A diagnostic text model that declines the cases it's least sure about, tuned to reach a chosen accuracy.
- Applicability: A1. Reinforces that a single promptable model spanning multiple NL task shapes can share one abstention mechanism — aligned with MORPHEUS's NL task auto-detection.
- Novelty implication: Supports MORPHEUS's promptable-abstention direction; the differentiator remains biology grounding, not text-only diagnosis.

9. **Uncertainty-Aware Deep Learning Classification of Adamantinomatous Craniopharyngioma from Preoperative MRI** (E. Prince et al., *Cancers/Diagnostics*, 2023) — https://pmc.ncbi.nlm.nih.gov/articles/PMC10047069/
- Takeaway: an abstention policy converts a modest classifier into a high-accuracy selective one on a rare tumor.
- Technical summary: Uses predictive-uncertainty thresholding on an MRI CNN; reported accuracy improves from ~80.8% to ~95.5% at a ~34% abstention rate. Demonstrates the accuracy-coverage trade-off on a small, rare-disease cohort.
- Plain-English: By skipping its one-in-three least-confident cases, a tumor classifier becomes far more accurate on the rest.
- Applicability: A2. Concrete evidence that abstention is most valuable in low-data / rare-programme regimes — the exact regime where per-pathway addressable slots matter.
- Novelty implication: Supporting datapoint for MORPHEUS's abstention value proposition in rare settings.

---

### Conformal prediction & distribution-free UQ for survival / risk

10. **Conformal Risk Control** (A.N. Angelopoulos, S. Bates, A. Fisch, L. Lei, T. Schuster, arXiv 2208.02814, ICLR 2024) — https://arxiv.org/abs/2208.02814
- Takeaway: extends split conformal prediction from coverage to controlling the expectation of any monotone loss, with a distribution-free guarantee.
- Technical summary: Generalizes conformal prediction so a tunable threshold controls E[loss] (e.g., false-negative rate, miscoverage) at a user-set level; subsumes standard coverage guarantees. Model-agnostic post-hoc wrapper.
- Plain-English: A statistical wrapper that lets you promise "this system's error rate stays below X" for many error types, regardless of the model inside.
- Applicability: A1, A5. Provides the formal risk guarantee MORPHEUS decision-support answers could carry per prompted task; loss-agnostic control suits heterogeneous NL-detected tasks and could bound risk on counterfactual queries.
- Novelty implication: Reframes MORPHEUS UQ claims: rigor should be stated as conformal risk control, not ad hoc confidence. Strengthens if adopted, pre-empts if ignored.

11. **Conformal predictive intervals in survival analysis: a resampling approach** (2025) — https://pmc.ncbi.nlm.nih.gov/articles/PMC12104816/
- Takeaway: builds distribution-free predictive intervals for survival time via resampling, valid under censoring.
- Technical summary: Adapts conformal machinery to right-censored survival with a resampling scheme to restore exchangeability; produces prediction intervals with finite-sample coverage. Evaluated on clinical survival datasets.
- Plain-English: Honest "the event will likely happen within this window" intervals for survival predictions, even when many patients' outcomes are censored.
- Applicability: A2, A4. Turns a prognosis slot's point estimate into a guaranteed interval — the kind of calibrated output MORPHEUS should expose per addressable programme.
- Novelty implication: Supports interval-valued prognosis outputs; MORPHEUS novelty is promptable delivery, not the conformal method.

12. **Weighted Conformal Prediction for Survival Analysis under Covariate Shift** (arXiv 2512.03738, 2025) — https://arxiv.org/abs/2512.03738
- Takeaway: keeps conformal survival coverage valid when the test population differs from training (covariate shift).
- Technical summary: Uses likelihood-ratio weighting to re-weight calibration scores under covariate shift and censoring, restoring coverage guarantees off-distribution. Tested on shifted survival cohorts.
- Plain-English: Makes survival-prediction guarantees hold up even when the new patients don't look like the training patients.
- Applicability: A4. Directly relevant to MORPHEUS deployment across cohorts/modalities with shift; motivates shift-aware UQ on encoded modalities.
- Novelty implication: Pre-empts any implicit MORPHEUS assumption of in-distribution validity; strengthens claims only if shift is handled explicitly.

13. **Two-sided conformalized survival analysis** (arXiv 2410.24136, 2024) — https://arxiv.org/abs/2410.24136
- Takeaway: constructs both lower and upper conformal bounds on survival time under censoring.
- Technical summary: Extends one-sided conformalized survival (lower-bound) methods to two-sided intervals with coverage guarantees despite Type-I censoring, giving bounded rather than half-open predictive regions.
- Plain-English: Instead of only "you'll survive at least this long," it gives a calibrated upper end too.
- Applicability: A2. Bounded prognosis intervals per slot; informs how MORPHEUS should present survival answers.
- Novelty implication: Background method; supports calibrated interval outputs.

---

### Active feature acquisition / next-best-test

14. **Evaluation of Active Feature Acquisition Methods for Time-varying Feature Settings** (H. von Kleist, A. Zamanian, I. Shpitser, N. Ahmidi, *JMLR*, 2025; arXiv 2312.01530) — https://arxiv.org/abs/2312.01530
- Takeaway: formalizes how to *retrospectively* evaluate a next-best-test agent from logged data without deploying it (AFAPE).
- Technical summary: Frames active feature acquisition performance evaluation under two assumptions (no direct effect of acquisition on feature values; no unobserved confounding), enabling missing-data and offline-RL estimators. Introduces a semi-offline RL framework with direct-method, IPW, and double-RL estimators requiring only a weaker positivity assumption.
- Plain-English: A rigorous way to ask "would ordering tests differently have helped?" using only records you already have, without experimenting on patients.
- Applicability: A5. The causal/offline-evaluation backbone MORPHEUS needs if it claims next-best-test value; connects AFA to counterfactual "what if we had measured X" queries.
- Novelty implication: Strongly reframes MORPHEUS's decision-support evaluation — retrospective AFA claims must use this kind of causal estimator, not naive replay. Pre-empts sloppy eval.

15. **Evaluation of Active Feature Acquisition Methods for Static Feature Settings** (H. von Kleist, A. Zamanian, I. Shpitser, N. Ahmidi, arXiv 2312.03619, 2023) — https://arxiv.org/abs/2312.03619
- Takeaway: extends AFAPE to static (time-invariant) features with IPW / direct-method / double-RL estimators handling MAR and MNAR missingness.
- Technical summary: Adapts the semi-offline RL evaluation framework to one-shot feature panels; derives three estimators and shows improved data efficiency of semi-offline estimators on synthetic and real data under MAR/MNAR.
- Plain-English: The same "would different testing have helped" question, for one-time panels of tests rather than repeated measurements.
- Applicability: A5. Applies to molecular-panel ordering (e.g., which omics assay to run) — a MORPHEUS next-best-assay use case.
- Novelty implication: Pre-empts naive next-best-assay claims; supplies the correct estimator family.

16. **Learning-To-Measure: In-context Active Feature Acquisition** (Y. Kobayashi, Z. Jing, J. Yao, H. Namkoong, S. Joshi, arXiv 2510.12624, 2024) — https://arxiv.org/abs/2510.12624
- Takeaway: a transformer decides the next feature to acquire *in-context*, generalizing across tasks without retraining.
- Technical summary: Trains on diverse tasks so a transformer conditions on observed features + a few in-context examples to select the next measurement, replacing per-step retraining and hand-built information-theoretic scores. Beats EMC, info-theoretic, and random baselines while avoiding retraining overhead.
- Plain-English: One model that learns a general "what should I measure next" instinct and applies it to new problems on the fly.
- Applicability: A1, A5. Closest external analogue to MORPHEUS's promptable, no-retrain next-best-test routing; validates that acquisition policy can be an in-context capability of a unified model rather than a hard-coded probe.
- Novelty implication: **Candidate pre-emption risk** for A1/A5 — in-context AFA already exists. MORPHEUS must differentiate by grounding the acquisition query in identified biological programmes and NL task auto-detection, not tabular features.

17. **NOCTA: Non-Greedy Objective Cost-Tradeoff Acquisition for Longitudinal Data** (arXiv 2507.12412, 2025) — https://arxiv.org/abs/2507.12412
- Takeaway: non-greedy next-best-test policy that trades information gain against acquisition cost over a longitudinal horizon.
- Technical summary: Plans acquisitions to optimize a long-horizon objective rather than greedy per-step utility, explicitly modeling cost vs benefit across time. Evaluated on longitudinal clinical data.
- Plain-English: Chooses tests looking several steps ahead and weighing each test's price against its expected value.
- Applicability: A5. Cost-aware, horizon-aware acquisition is the form a clinically credible MORPHEUS next-best-test query should take.
- Novelty implication: Pre-empts "cost-aware acquisition" novelty; MORPHEUS differentiator is biological addressability + NL prompting.

18. **Towards Dynamic Feature Acquisition on Medical Time Series by Maximizing Conditional Mutual Information** (arXiv 2407.13429, 2024) — https://arxiv.org/html/2407.13429
- Takeaway: greedily acquires the measurement that maximizes conditional mutual information with the outcome on clinical time series.
- Technical summary: Estimates CMI between candidate measurements and the target given already-observed data and acquires the maximizer; handles irregular sampling in medical time series. Reports cost-accuracy trade-offs vs baselines.
- Plain-English: At each step, order the test expected to tell you the most about the patient's outcome.
- Applicability: A5. Information-theoretic acquisition criterion MORPHEUS could expose as an interpretable next-test justification.
- Novelty implication: Background method; MORPHEUS should route such criteria through identified slots.

19. **Distribution Guided Active Feature Acquisition** (arXiv 2410.03915, 2024) — https://arxiv.org/abs/2410.03915
- Takeaway: uses a learned generative model of the feature distribution to guide which feature to acquire next.
- Technical summary: A generative surrogate imputes/estimates unobserved features and their value-of-information, guiding acquisition without exhaustive rollouts. Improves sample efficiency over model-free RL AFA.
- Plain-English: Uses a model of "what typical patients look like" to guess which missing test would be most worth doing.
- Applicability: A4, A5. Generative surrogate over modalities connects to encode-vs-RAG: a modality worth encoding is one the surrogate can't cheaply impute.
- Novelty implication: Supports a value-of-information framing for MORPHEUS modality decisions.

20. **Deep Sensing: Active Sensing using Multi-directional Recurrent Neural Networks** (J. Yoon, W.R. Zame, M. van der Schaar, ICLR 2018) — https://openreview.net/pdf?id=r1SnX5xCb
- Takeaway: foundational active-sensing model that jointly learns what and when to measure under measurement cost.
- Technical summary: An M-RNN imputes missing streams (using lagged and advanced timing) and an active-sensing policy selects the most informative measurements at test time to maximize prediction performance per unit cost. Evaluated on ICU prediction tasks.
- Plain-English: An early system that learns both to fill in missing vitals and to pick which vital to actually measure next.
- Applicability: A5. Historical anchor for next-best-test; situates MORPHEUS's contribution as promptable, biology-grounded active sensing.
- Novelty implication: Establishes prior art; MORPHEUS novelty must be above generic active sensing.

21. **Dynamic Measurement Scheduling for Event Forecasting using Deep Reinforcement Learning** (C. Chang et al., arXiv 1901.09699, 2019) — https://arxiv.org/abs/1901.09699
- Takeaway: deep-RL scheduler that decides which lab panels to order and when, accounting for panel structure and cost.
- Technical summary: Trains an RL policy on ICU data to schedule measurements (tests come in panels of differing cost) so as to maximize event-forecasting utility net of cost. Demonstrates cost savings at matched forecasting performance.
- Plain-English: Learns a test-ordering schedule that forecasts deterioration while avoiding unnecessary labs.
- Applicability: A5. Panel-structured, cost-aware scheduling is the realistic next-best-test setting; MORPHEUS assays also come in panels.
- Novelty implication: Prior art for cost/panel-aware scheduling; reinforces that MORPHEUS's edge is biological addressability + NL query, not scheduling per se.

22. **A reinforcement learning guided adaptive cost-sensitive feature acquisition method** (*Applied Soft Computing*, 2022) — https://www.sciencedirect.com/science/article/abs/pii/S1568494622000163
- Takeaway: RL policy that adaptively balances misclassification cost against feature-acquisition cost.
- Technical summary: Jointly minimizes acquisition and misclassification costs via an RL agent that acquires features until expected net benefit turns negative; compared against static cost-sensitive baselines. Reports improved cost-accuracy trade-offs.
- Plain-English: Keeps ordering tests only while each new test is worth more than it costs.
- Applicability: A5. Net-benefit stopping rule MORPHEUS could surface as the rationale for stopping acquisition.
- Novelty implication: Background; supports value-of-information stopping in decision support.

---

### Multimodal prognosis / survival from cancer data (encode-vs-RAG, slots, fusion)

23. **Modeling Dense Multimodal Interactions Between Biological Pathways and Histology for Survival Prediction (SurvPath)** (G. Jaume, A. Vaidya, R.J. Chen, D. Williamson, P.P. Liang, F. Mahmood, CVPR 2024; arXiv 2304.06819) — https://arxiv.org/abs/2304.06819
- Takeaway: tokenizes transcriptomics into *biological pathway tokens* and fuses them with histology patch tokens via a memory-efficient transformer, giving pathway-level interpretability.
- Technical summary: Transcriptomics is grouped into pathway tokens (each encoding a cellular function) rather than raw genes; histology WSIs give patch tokens; a memory-efficient multimodal transformer models dense pathway<->patch interactions. Achieves SOTA over unimodal and multimodal baselines on five TCGA datasets and surfaces genotype-phenotype interaction as prognostic factors.
- Plain-English: Instead of feeding thousands of raw genes, it bundles them into biological "pathways" and lets the model reason about how each pathway interacts with what the tumor looks like.
- Applicability: A2, A3, A4. This is the strongest external validation of MORPHEUS's *pathway-addressable slots*: pathway tokens are exactly the identified, per-programme units MORPHEUS proposes, and cross-modal attention grounds them in histology.
- Novelty implication: **Strong pre-emption risk for A2.** Pathway tokenization + cross-modal attention already exists and is interpretable. MORPHEUS must claim novelty in *promptability* (A1), NL grounding/emergent-knowledge elicitation (A3), and interventional queries (A5) — not merely pathway tokens for prognosis.

24. **Multimodal prototyping for cancer survival prediction (MMP)** (A.H. Song, R.J. Chen, G. Jaume, et al., ICML 2024) — https://github.com/mahmoodlab/MMP
- Takeaway: compresses each modality (WSI, transcriptomics) into a small set of morphological/pathway *prototypes* before fusion, improving efficiency and interpretability.
- Technical summary: Learns prototype tokens per modality via optimal-transport-style assignment, reducing gigapixel WSIs and high-dim omics to compact prototype sets that a fusion module combines for survival prediction. Matches or beats dense-token baselines with far less compute.
- Plain-English: Summarizes a whole slide and a genome into a handful of representative "concepts," then predicts survival from those.
- Applicability: A2, A4. Prototypes are a candidate mechanism for identified slots; supports the idea that a frozen trunk can expose compact, addressable units.
- Novelty implication: Pre-empts "compact interpretable slot" novelty; MORPHEUS differentiator remains prompt-time addressability and NL routing.

25. **Pathology-and-Genomics Multimodal Transformer for Survival Outcome Prediction (PathOmics)** (K. Ding et al., MICCAI 2023) — https://conferences.miccai.org/2023/papers/485-Paper1847.html
- Takeaway: pretrains a multimodal transformer to align pathology and genomics, then fine-tunes for survival, tolerating missing genomics at inference.
- Technical summary: Uses unsupervised multimodal pretraining to capture pathology-genomics interactions, then supervised finetuning; the aligned pathology branch can predict survival even when genomics is unavailable at test time. Evaluated on TCGA colorectal/gastric cohorts.
- Plain-English: Learns how tissue images and gene data relate, so it can still make good survival calls when the gene data is missing.
- Applicability: A4. Directly informs MORPHEUS's encode-vs-RAG and missing-modality behaviour: an aligned trunk degrades gracefully when a modality is absent.
- Novelty implication: Pre-empts "graceful missing-modality multimodal survival"; MORPHEUS must add promptability/causal queries on top.

26. **Machine learning-based multimodal prognostic models integrating pathology images and high-throughput omic data for overall survival prediction in cancer: a systematic review** (arXiv 2507.16876, 2025) — https://arxiv.org/abs/2507.16876
- Takeaway: systematic review establishing that multimodal (histology+omics) models generally beat unimodal for OS prediction, with recurring methodological gaps.
- Technical summary: Reviews fusion strategies (early/late/intermediate, attention, graph), datasets (mostly TCGA), and evaluation; flags small cohorts, weak external validation, and inconsistent calibration/UQ reporting as field-wide weaknesses.
- Plain-English: A survey confirming that combining tissue images with molecular data usually predicts survival better, but the field's evaluation is often shaky.
- Applicability: A4. Baseline evidence and a gap list (external validation, UQ, calibration) MORPHEUS can target.
- Novelty implication: Reframes MORPHEUS opportunity — the open gaps (calibration, external validity, promptability) are where novelty should sit, not fusion accuracy.

27. **Multimodal deep learning for cancer prognosis prediction with clinical information prompts integration** (npj Digital Medicine, 2025) — https://www.nature.com/articles/s41746-025-02257-y
- Takeaway: integrates clinical information as *prompts* into a multimodal prognosis model.
- Technical summary: Injects structured clinical variables as prompt tokens conditioning a histology/omics prognosis model, improving survival prediction over non-prompted fusion. Evaluated on multi-cancer cohorts.
- Plain-English: Feeds a patient's clinical facts to the model as a prompt that steers its survival prediction.
- Applicability: A1, A4. Demonstrates prompt-conditioned multimodal prognosis — adjacent to MORPHEUS's promptable representation, but with clinical (not NL-task) prompts.
- Novelty implication: **Candidate pre-emption for A1/A4** (prompts already used in multimodal prognosis). MORPHEUS must show NL *task auto-detection* and biology grounding beyond clinical-variable prompting.

28. **Continually Evolved Multimodal Foundation Models for Cancer Prognosis** (J. Peng, S. Zhou, ..., T. Chen, arXiv 2501.18170, 2025) — https://arxiv.org/abs/2501.18170
- Takeaway: a foundation-model framework that incorporates newly arriving modalities/distributions without full retraining.
- Technical summary: Targets continual incorporation of new data distributions into a multimodal prognosis model, moving beyond concatenation to model cross-modal interdependencies; validated on TCGA. Emphasizes adaptive integration over task-specific pipelines.
- Plain-English: A prognosis model that keeps absorbing new kinds of data over time instead of being rebuilt from scratch.
- Applicability: A4. Aligns with MORPHEUS's frozen-trunk plug-in ambition (add a modality without retraining the trunk).
- Novelty implication: Pre-emption risk for the "add-modality-without-retrain" framing; MORPHEUS must anchor its version in identified slots + prompting.

29. **BioFusionNet: Survival Risk Stratification in ER+ Breast Cancer through Multifeature and Multimodal Data Fusion** (R. Al-Thani et al., arXiv 2402.10717, *IEEE JBHI*, 2024) — https://arxiv.org/abs/2402.10717
- Takeaway: attention-weighted fusion of histology, genomics, and clinical data for breast-cancer survival risk stratification.
- Technical summary: Combines self-supervised histology features, gene expression, and clinical variables via a co-dual cross-attention fusion and a weighted Cox loss for imbalance; reports improved C-index and risk-group separation on ER+ breast cancer.
- Plain-English: Merges tissue, gene, and clinical signals with attention to sort ER+ breast-cancer patients into risk groups.
- Applicability: A4. Concrete fusion-architecture datapoint and a Cox-loss handling of imbalance.
- Novelty implication: Background baseline; MORPHEUS differentiator is promptability/UQ, not another fusion net.

30. **Deep learning-driven survival prediction in pan-cancer studies by integrating multimodal histology-genomic data** (2025) — https://pmc.ncbi.nlm.nih.gov/articles/PMC11926983/
- Takeaway: pan-cancer histology+genomics survival model showing cross-cancer generalization of multimodal fusion.
- Technical summary: Integrates WSI and genomic features across many TCGA cancer types with an attention fusion, evaluating pan-cancer and per-cancer survival; multimodal beats unimodal broadly.
- Plain-English: A single multimodal survival model trained across many cancer types, confirming the image+gene combo travels across cancers.
- Applicability: A4. Evidence that a unified multimodal trunk can span cancer types — a prerequisite for MORPHEUS's unified representation claim.
- Novelty implication: Supports feasibility of a pan-cancer unified trunk; novelty must be promptability atop it.

---

### Missing-modality robustness & prompt-based modality handling

31. **Distilled Prompt Learning for Incomplete Multimodal Survival Prediction (DisPro)** (Y. Xu, F. Zhou, C. Zhao, Y. Wang, C. Yang, H. Chen, CVPR 2025; arXiv 2503.01653) — https://arxiv.org/abs/2503.01653
- Takeaway: uses available modalities as *prompts to an LLM* to infer missing modalities, with knowledge distillation restoring modality-specific detail.
- Technical summary: Two-stage prompting — UniPro distills each modality's knowledge distribution; MultiPro uses available modalities as prompts to an LLM to infer the missing modality (common cross-modal info) while injecting distilled unimodal knowledge. Outperforms baselines across missing-modality scenarios.
- Plain-English: When gene data is missing, it prompts a language model with the tissue data to reconstruct the gist of the missing gene signal.
- Applicability: A1, A3, A4. Direct precedent for MORPHEUS's encode-vs-RAG-vs-prompt spectrum: modalities function as prompts, and an LLM supplies missing-modality context — the RAG-context idea in action.
- Novelty implication: **Strong pre-emption risk for A4.** "Modalities as prompts, LLM infers missing modality" already exists. MORPHEUS must claim identifiability/addressability (A2) and interventional queries (A5) that DisPro does not offer.

32. **Deep Multimodal Learning with Missing Modality: A Survey** (R. Wu et al., arXiv 2409.07825, 2024) — https://arxiv.org/abs/2409.07825
- Takeaway: taxonomy of missing-modality strategies (imputation, dropout, distillation, prompt-based) and their trade-offs.
- Technical summary: Categorizes methods by whether they reconstruct, discard, or adapt to missing modalities, and by training-time vs test-time missingness; discusses evaluation protocols for incomplete multimodal learning.
- Plain-English: A map of every way to cope when one data type is missing at prediction time.
- Applicability: A4. Frames the design space for MORPHEUS's encode-vs-RAG decision and graceful degradation.
- Novelty implication: Background; positions MORPHEUS's approach within known strategies.

33. **Robust Multimodal Survival Prediction with the Latent Differentiation Conditional Variational AutoEncoder** (arXiv 2503.09496, 2025) — https://arxiv.org/abs/2503.09496
- Takeaway: a conditional VAE learns modality-shared vs modality-specific latents so survival prediction is robust to a missing modality.
- Technical summary: Disentangles a shared latent (recoverable from any modality) from modality-specific latents; at inference, missing modalities' specific latents are sampled from the conditional prior, preserving the shared signal. Improves robustness under incomplete inputs on TCGA survival tasks.
- Plain-English: Splits what's common across data types from what's unique to each, so losing one type keeps the common signal intact.
- Applicability: A2, A4. Shared/specific disentanglement is a candidate identifiability mechanism for MORPHEUS slots.
- Novelty implication: Pre-empts "disentangled robust survival" novelty; MORPHEUS must tie disentanglement to *addressable biological programmes*, not abstract latents.

34. **PRIME: Prototype-Driven Multimodal Pretraining for Cancer Prognosis with Missing Modalities** (arXiv 2604.04999, 2026) — https://arxiv.org/abs/2604.04999
- Takeaway: prototype-based multimodal pretraining that stays robust when modalities are missing at prognosis time.
- Technical summary: Learns cross-modal prototypes during pretraining that anchor a modality-agnostic representation, so downstream prognosis tolerates absent modalities. Evaluated on TCGA multimodal cohorts.
- Plain-English: Pretrains reusable "concept anchors" shared by all data types, so prognosis holds up when some data is absent.
- Applicability: A2, A4. Prototypes-as-anchors again suggests a slot mechanism robust to missingness.
- Novelty implication: Reinforces prototype/slot prior art; differentiator is promptable addressability.

---

### Treatment-response prediction & causal/counterfactual decision support

35. **BITES: balanced individual treatment effect for survival data** (S. Schrod, ..., M. Altenbuchinger, *Bioinformatics*, 2022) — https://academic.oup.com/bioinformatics/article/38/Supplement_1/i60/6617509
- Takeaway: deep counterfactual survival model that recommends individualized treatment by estimating balanced individual treatment effects.
- Technical summary: Treatment-specific Cox-DL outcome heads with a Sinkhorn/IPM penalty balancing treated vs control latent distributions to correct selection bias. On Rotterdam-trained / GBSG-validated node-positive breast cancer, it recommends hormone therapy for only 83.4% of patients and yields the strongest survival separation between recommended vs anti-recommended groups (P=1.6e-5); SHAP shows menopausal status and grade modulate treatment benefit.
- Plain-English: Predicts, per patient, whether hormone therapy will actually help, and its recommendations separate responders from non-responders in an independent trial.
- Applicability: A5. This is the counterfactual "treatment as a query" pattern MORPHEUS wants: predict outcome under treat vs no-treat, not a retrained classifier. Provides an evaluation template (recommended vs anti-recommended survival separation).
- Novelty implication: **Anchors A5.** Counterfactual survival ITE already exists in oncology. MORPHEUS must frame drug/perturbation as a *prompt* over identified programmes and generalize beyond binary treatment, else A5 is pre-empted.

36. **Enabling counterfactual survival analysis with balanced representations** (P. Chapfuwa et al., ACM CHIL 2021) — https://www.researchgate.net/publication/350736076
- Takeaway: balanced-representation neural network predicts factual and counterfactual survival curves, yielding CATE under censoring.
- Technical summary: Learns a representation invariant to treatment assignment (IPM balancing) and predicts treatment-specific survival distributions, giving counterfactual survival and CATE despite right-censoring. Evaluated on clinical trial and observational survival data.
- Plain-English: Estimates each patient's survival both if treated and if not, to reveal who benefits.
- Applicability: A5. Method precedent for counterfactual prognosis under censoring; informs how MORPHEUS could answer perturbation queries.
- Novelty implication: Prior art for counterfactual survival; MORPHEUS differentiator is promptable, biology-addressable counterfactuals.

37. **Estimating individual treatment effect: generalization bounds and algorithms (CFR-Net / TARNet)** (U. Shalit, F. Johansson, D. Sontag, ICML 2017) — https://proceedings.mlr.press/v70/shalit17a/shalit17a.pdf
- Takeaway: foundational representation-balancing framework (TARNet/CFR) with generalization bounds for individual treatment-effect estimation.
- Technical summary: Introduces treatment-specific heads over a shared representation with an IPM (Wasserstein/MMD) penalty bounding the counterfactual generalization error; establishes the theory most oncology ITE methods build on.
- Plain-English: The theoretical blueprint for predicting "what would happen under the other treatment" from observational data.
- Applicability: A5. The causal-geometry backbone (balanced representations) MORPHEUS's interventional queries would rest on.
- Novelty implication: Establishes A5 prior art at the theory level; MORPHEUS novelty must be architectural/promptable, not the ITE principle.

38. **Informing immunotherapy with multi-omics driven machine learning** (N. Kang et al., *npj Digital Medicine*, 2024) — https://www.nature.com/articles/s41746-024-01043-6
- Takeaway: multi-omics ML for immunotherapy response prediction and biomarker elicitation.
- Technical summary: Integrates multiple omics layers to predict immune-checkpoint response and surfaces candidate predictive biomarkers; positions ML as both predictor and hypothesis generator for immunotherapy. (Access-gated; details from listing.)
- Plain-English: Uses many molecular data types together to predict who responds to immunotherapy and to flag why.
- Applicability: A3, A4. Connects multimodal encoding to *emergent biomarker knowledge* — relevant to MORPHEUS's emergent-knowledge elicitation and its evaluation.
- Novelty implication: Supports A3's biomarker-elicitation framing; MORPHEUS must make elicitation *promptable and evaluated*, not a post hoc feature-importance readout.

39. **Cancer immunotherapy response prediction from multi-modal clinical and image data using semi-supervised deep learning** (X. Wang et al., *Radiotherapy and Oncology*, 2023; PMID 37414254) — https://www.sciencedirect.com/science/article/abs/pii/S0167814023003316
- Takeaway: semi-supervised multimodal (imaging+clinical) model for immunotherapy response, improved by a tissue biomarker.
- Technical summary: Uses a large unlabeled cohort via semi-supervised learning to improve generalization; combining the model with a clinically approved tissue biomarker further raises accuracy. Reports response-prediction AUCs on internal/external cohorts.
- Plain-English: Learns from many unlabeled scans to better predict immunotherapy response, and gets even better paired with a lab biomarker.
- Applicability: A4. Shows an ML prediction complementing (not replacing) an established biomarker — a template for MORPHEUS treating a biomarker as RAG context.
- Novelty implication: Supports encode-vs-RAG framing where a validated biomarker is context, not an encoded modality.

40. **Non-invasive multimodal CT deep learning biomarker to predict pathological complete response of NSCLC following neoadjuvant immunochemotherapy: a multicenter study** (2024) — https://pmc.ncbi.nlm.nih.gov/articles/PMC11409329/
- Takeaway: multimodal CT-derived deep-learning biomarker predicts pathologic complete response before therapy, externally validated.
- Technical summary: Combines multi-phase/region CT deep features (and clinical data) to predict pCR to neoadjuvant immunochemotherapy across multiple centers, reporting AUCs with external validation — a rigor level often missing in the field.
- Plain-English: Reads pretreatment CT scans to forecast which lung-cancer patients' tumors will fully respond, tested across hospitals.
- Applicability: A4. Multicenter external validation is the retrospective-eval standard MORPHEUS treatment-response claims should meet.
- Novelty implication: Sets external-validation bar; supports encode-imaging-for-response but not promptability.

---

### Foundation-model / broader methodological context

41. **Multimodal deep learning-based prognostication in glioma patients: a systematic review** (M. Nafe et al., 2023) — https://pmc.ncbi.nlm.nih.gov/articles/PMC9856816/
- Takeaway: disease-focused review confirming multimodal > unimodal for glioma prognosis, with data/validation caveats.
- Technical summary: Synthesizes imaging+molecular+clinical prognostic models for glioma, cataloguing modality combinations (e.g., MGMT status, pathology) and their incremental value, and flags reproducibility gaps.
- Plain-English: For brain tumors, combining scan, molecular, and clinical data predicts outcome better — but studies vary in quality.
- Applicability: A4. Concrete per-disease evidence for which modality combinations add prognostic value.
- Novelty implication: Background; supports multimodal value, novelty must be elsewhere.

42. **Pan-cancer survival prediction using a deep learning architecture with multimodal representation and integration (MultiSurv-style)** (L. Vale-Silva & K. Rohr region; *Bioinformatics Advances*, 2023) — https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad006/6998218
- Takeaway: a unified multimodal architecture predicts survival across many cancer types from combined clinical/omics/imaging modalities.
- Technical summary: Modality-specific submodels feed a fusion layer producing conditional survival distributions across 33 cancer types; supports flexible modality subsets at inference. Reports strong pan-cancer C-index.
- Plain-English: One survival model spanning dozens of cancers that works with whatever data types are available.
- Applicability: A4. Demonstrates flexible-modality-subset inference in a unified model — a MORPHEUS prerequisite.
- Novelty implication: Prior art for flexible-modality unified survival; MORPHEUS edge is promptability/UQ/causal queries.

43. **RadFlag: A Black-Box Hallucination Detection Method for Medical Vision-Language Models** (S. Zhang et al., arXiv 2411.00299, 2024) — https://arxiv.org/abs/2411.00299
- Takeaway: black-box, sampling-consistency method flags likely-hallucinated findings in medical VLM outputs for selective abstention.
- Technical summary: Samples multiple generations at varying temperature and flags low-consistency claims as hallucinations, enabling report-level abstention/referral without model internals. Improves precision on radiology report generation by withholding flagged findings.
- Plain-English: Detects when a medical AI is likely making things up by checking whether it says the same thing across repeated tries.
- Applicability: A1, A3. If MORPHEUS elicits NL/biology statements, a black-box consistency check gives abstention for generated biological claims — relevant to evaluating emergent-knowledge reliability.
- Novelty implication: Supports A3 evaluation tooling (flagging unreliable elicited knowledge); MORPHEUS should pair emergent-knowledge claims with such consistency checks.

44. **A Meta-Learner Framework to Estimate Individualized Treatment Effects for Survival Outcomes** (2025) — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12440118/
- Takeaway: meta-learner (T-/X-learner) framework for individualized treatment effects on censored survival outcomes.
- Technical summary: Adapts meta-learner CATE estimation to survival by pairing pseudo-outcomes/base learners (random survival forest, Bayesian AFT, survival NNs) with T- and X-learner constructions; discusses censoring-aware estimation and time-dependence limitations.
- Plain-English: A flexible recipe for estimating each patient's treatment benefit on time-to-event outcomes using off-the-shelf survival learners.
- Applicability: A5. Provides model-agnostic ITE machinery MORPHEUS could wrap around a frozen trunk for counterfactual queries.
- Novelty implication: Prior art for survival CATE; MORPHEUS differentiator is promptable, representation-level counterfactuals.

45. **Multimodal deep learning-based prognostication via clinical-prompt integration and other prompt-conditioned prognosis** — see entry 27; cross-listed for A1. (npj Digital Medicine, 2025) — https://www.nature.com/articles/s41746-025-02257-y
- Takeaway: (cross-reference) prompt-conditioned multimodal prognosis is an active, near-neighbor line to MORPHEUS A1.
- Applicability: A1, A4. Explicitly flagged as the closest "prompt + multimodal prognosis" prior art.
- Novelty implication: Reiterated pre-emption watch for A1 — MORPHEUS must show NL *task auto-detection/routing*, not just clinical-variable prompting.

---

## Cross-cutting notes for MORPHEUS

- **Strongest slot precedent (A2):** SurvPath's pathway tokens (#23) and prototype methods (#24, #34) show identified, interpretable units are already used for prognosis. MORPHEUS's identifiability/addressability claim must go beyond "pathway tokens exist" to *prompt-time per-programme addressability* with a stated identifiability argument.
- **Prompting + multimodal already exists (A1/A4):** #27 (clinical prompts), #31 (modalities-as-prompts + LLM infers missing modality) are the sharpest pre-emption risks. MORPHEUS's defensible novelty is NL task auto-detection/routing and biology grounding, not "prompts in a multimodal model."
- **Counterfactual treatment queries already exist (A5):** #35-#37 establish counterfactual survival/ITE as prior art. MORPHEUS must recast drug/perturbation as a query over a frozen representation's causal geometry (not a retrained ITE head) to be novel.
- **UQ/eval is the moat, if claimed rigorously:** #1 (calibration audit), #10 (conformal risk control), #14-#15 (retrospective AFA evaluation via causal estimators) collectively define the evaluation bar. Adopting conformal risk control + AFAPE-style offline evaluation + calibration auditing would let MORPHEUS make decision-support claims that most cited papers cannot back.
