## Adversarial prior-art / steelman-the-null

Lane id: `l15_steelman_prior_art`. Remit: actively hunt papers that ALREADY do what MORPHEUS hopes to claim across A1-A5 — existing promptable multimodal cancer models, existing NL-queryable tumor/cell representations, existing encode-vs-retrieve studies, existing emergent-biology evaluations, existing perturbation/intervention-as-query systems. This lane is deliberately hostile: its job is to surface scoops, pre-emptions, and the strongest version of "MORPHEUS's headline is already done."

**Cross-cutting takeaway for MORPHEUS.** The individual *mechanisms* MORPHEUS wants to claim are, almost without exception, already published: (a) NL-promptable, task-inferring multimodal medical models exist and are strong (Med-PaLM M, BiomedGPT, RadFM, PathChat, BiomedParse); (b) NL-queryable *cell/tumor* representations exist (Cell2Sentence, GenePT, scGPT, Geneformer); (c) perturbation/drug-as-a-query from a learned latent exists and is mature (scGen latent arithmetic, CPA, chemCPA, GEARS, CellOT, Biolord); (d) pathway-*addressable* multimodal tumor tokens exist (SurvPath); (e) the "emergent biological knowledge" of these models has *already been adversarially evaluated and often found wanting* (Kedzierska, Boiarsky, Ahlmann-Eltze — linear baselines match or beat FMs). The honest novelty gap for MORPHEUS is therefore NOT any single axis primitive but the **synthesis**: a single *frozen, identified* oncology trunk that (i) auto-routes NL tasks, (ii) exposes *biologically identified* pathway slots (not just attention tokens named after pathways), (iii) accepts perturbation queries *without retraining*, and (iv) is evaluated with an emergence/elicitation benchmark that the prior-art critiques show the field currently lacks. Every claim below is a real, findable paper; the ones that most threaten a naive MORPHEUS pitch are flagged in "Novelty implication."

---

### Group A — Promptable / instruction-following multimodal medical & cancer models (A1, A4)

**1. Towards Generalist Biomedical AI (Med-PaLM M)** (Tu, Azizi, Driess, Schaekermann et al., arXiv 2023) — https://arxiv.org/abs/2307.14334
- *Takeaway:* One generative model that "flexibly encodes and interprets" language, imaging, and genomics, task-specified purely by the text prompt, across 14 tasks (MultiMedBench).
- *Technical summary:* A single PaLM-E-derived multimodal model is instruction-tuned to route 14 heterogeneous biomedical tasks (VQA, report generation/summarization, genomic variant classification, image classification) from a unified NL interface; it matches or beats task-specific SOTA on all 14 and produces chest-X-ray reports clinicians preferred to radiologists' in up to 40.5% of cases. Task identity is *inferred from the prompt*, not selected by a hard-coded head.
- *Plain-English:* You tell one AI, in plain English, which medical question you want answered — about an image, a report, or a genome — and it figures out the task and answers it.
- *Applicability:* A1 (this is the reference implementation of NL-task-specified multimodal medical inference), A4 (it *encodes* genomics as a modality alongside imaging). Design implication: MORPHEUS's "NL task auto-detection" is not new as a capability; the defensible delta is a *frozen, identified oncology-specific trunk* with pathway-addressable slots, not a generalist instruction-tuned decoder.
- *Novelty implication:* **Strongly pre-empts A1 as stated.** "Promptable unified multimodal biomedical representation with task auto-detection" is Med-PaLM M. MORPHEUS must reframe A1 around *identifiability + frozen-trunk oncology specialization + causal query*, not around "you can prompt it."

**2. BiomedGPT: A Generalist Vision-Language Foundation Model for Diverse Biomedical Tasks** (Zhang et al., Nature Medicine 2024; arXiv 2305.17100) — https://arxiv.org/abs/2305.17100
- *Takeaway:* A lightweight open-source generalist that handles classification, VQA, captioning, and report generation across imaging+text from one prompt-driven interface, SOTA on 16/25 benchmarks.
- *Technical summary:* A single sequence-to-sequence transformer is multi-task pretrained then instruction-tuned across many biomedical modalities; the task is specified by the textual instruction, and one set of weights serves all tasks. Human evaluation shows low error rates in radiology QA (3.8%) and report generation (8.3%).
- *Plain-English:* A small, open version of a "do-everything" medical AI where the plain-English request selects the job.
- *Applicability:* A1 (prompt-selected multitask), A4 (multimodal encoding). Design implication: an open generalist baseline MORPHEUS's promptability must beat *on oncology-specific, identifiability-sensitive* metrics, since it cannot win on generalist breadth.
- *Novelty implication:* Pre-empts A1's "one model, many tasks, NL-selected." Reinforces that promptable-generalist is commodity; the oncology + identified-slot + causal angle is where novelty must live.

**3. LLaVA-Med: Training a Large Language-and-Vision Assistant for Biomedicine in One Day** (Li et al., NeurIPS 2023; arXiv 2306.00890) — https://arxiv.org/abs/2306.00890
- *Takeaway:* Cheap curriculum-style visual instruction tuning yields an open biomedical chat model that answers open-ended NL questions about medical images.
- *Technical summary:* PubMed Central figure-caption data plus GPT-4-generated instruction-following dialogues train a vision-language assistant with a two-stage curriculum (concept alignment then instruction following), reaching SOTA on biomedical VQA. Demonstrates that NL grounding to images is achievable with modest compute.
- *Plain-English:* A ChatGPT-for-medical-images built in under a day, that you converse with in natural language.
- *Applicability:* A1, A3 (NL⇄image grounding via instruction data). Design implication: the *grounding data recipe* (GPT-4-authored instructions over captioned images) is a known template; MORPHEUS's A3 novelty cannot be "we instruction-tune on captions."
- *Novelty implication:* Pre-empts the generic NL-grounding recipe. Neutral-to-cautionary for A3.

**4. Towards Generalist Foundation Model for Radiology (RadFM)** (Wu et al., arXiv 2023; publ. Nat. Commun. 2025) — https://arxiv.org/abs/2308.02463
- *Takeaway:* A radiology generalist over 2D+3D scans that interleaves image and text and is driven by NL task descriptions, beating GPT-4V and Med-Flamingo on RadBench.
- *Technical summary:* Visually-conditioned generative pretraining on 16M 2D/3D scans with paired text enables diverse radiologic tasks (VQA, report gen, diagnosis, rationale) from a unified prompt interface. Introduces RadBench for holistic NL-task evaluation.
- *Plain-English:* A single radiology AI you query in words, working on both flat and volumetric scans.
- *Applicability:* A1 (NL-task routing over volumetric medical data), A4. Design implication: 3D/volumetric prompting is already solved in imaging; MORPHEUS should not claim novelty for "multimodal + promptable" per se.
- *Novelty implication:* Pre-empts A1 in the imaging setting.

**5. A Multimodal Generative AI Copilot for Human Pathology (PathChat)** (Lu, Chen, Williamson et al., Nature 2024; arXiv 2312.07814) — https://arxiv.org/abs/2312.07814
- *Takeaway:* A pathology-specific vision-language assistant: histology image + NL clinical context in, free-form diagnostic reasoning out; 87% MCQ accuracy with context, preferred over GPT-4V.
- *Technical summary:* A histology vision encoder (100M images) is coupled to an LLM and instruction-tuned on 250k+ vision-language instructions, yielding an interactive "copilot" that answers open-ended pathology questions and integrates clinical text as *context*. Human pathologists rated its answers more accurate than GPT-4V.
- *Plain-English:* A pathologist's chat assistant that looks at the slide, reads the case notes, and reasons out loud.
- *Applicability:* A1 (NL-promptable cancer-image model), A4 (clinical text used as RAG-like *context* rather than encoded), A3 (NL diagnostic reasoning). Design implication: the *cancer-specific* promptable model already exists; MORPHEUS's differentiator is that PathChat is a chat wrapper over an encoder, with **no identified latent and no counterfactual query** — exactly MORPHEUS's intended delta.
- *Novelty implication:* **The single most threatening A1 pre-emption for oncology.** "Promptable multimodal cancer model" = PathChat. MORPHEUS must claim the *identified/causal/frozen-query* properties PathChat lacks, or it is a marginal PathChat variant.

**6. Towards a Visual-Language Foundation Model for Computational Pathology (CONCH)** (Lu et al., Nature Medicine 2024; arXiv 2307.12914) — https://arxiv.org/abs/2307.12914
- *Takeaway:* A contrastive vision-language pathology model (1.17M image-caption pairs) enabling zero-shot, *text-prompted* classification/segmentation/retrieval without fine-tuning.
- *Technical summary:* CoCa-style contrastive + captioning pretraining aligns histology tiles with pathology language, so downstream tasks are performed by *NL label prompts* (CLIP-style prototypes) with minimal/no supervised tuning, SOTA on 13 benchmarks. This is precisely the "text-prototype prompting" MORPHEUS's current stack already imitates.
- *Plain-English:* A pathology CLIP: name the thing you want to find in words and it finds it, no retraining.
- *Applicability:* A1 (zero-shot NL-prototype prompting), A3 (histology⇄language alignment). Design implication: MORPHEUS's existing "CLIP label prototypes" *are* CONCH-style; that is not the promptable-representation novelty. The gap is *task inference/routing* beyond fixed label prototypes.
- *Novelty implication:* Pre-empts "text-prompt the pathology model." Confirms the drift diagnosis: MORPHEUS's current prompting = CONCH-level, not TQI-level.

**7. BiomedParse: a foundation model for image parsing of everything everywhere all at once** (Zhao et al., Nature Methods 2025; arXiv 2405.12971) — https://arxiv.org/abs/2405.12971
- *Takeaway:* A single model that segments/detects/recognizes 82 object types across 9 modalities driven entirely by a **text prompt**, replacing per-task heads and manual boxes.
- *Technical summary:* 6M+ image-mask-text triples (GPT-4-aligned to biomedical ontologies) train a joint segmentation/detection/recognition model where the NL prompt names the target; it beats box-prompted baselines and enables "segment all instances of X" from language alone.
- *Plain-English:* Say what to outline in an image, in words, and it outlines it — across many kinds of medical images.
- *Applicability:* A1 (NL prompt → task/target selection), A3 (ontology-grounded language). Design implication: text-prompt-as-task-selector for biomedical images is done; MORPHEUS's routing must operate over a *representation/latent*, not pixels, to differ.
- *Novelty implication:* Pre-empts "NL prompt selects the task" in the image-parsing setting.

**8. A visual-language foundation model for pathology image analysis using medical Twitter (PLIP)** (Huang, Bianchi, Yuksekgonul, Montine, Zou; Nature Medicine 2023) — https://www.nature.com/articles/s41591-023-02504-3
- *Takeaway:* A pathology CLIP trained on social-media image-text pairs, enabling zero-shot NL classification and text-to-image retrieval.
- *Technical summary:* OpenPath (208k pathology image-text pairs from Twitter) fine-tunes CLIP into PLIP, giving zero-shot classification via text prompts and cross-modal retrieval, outperforming vanilla CLIP on pathology.
- *Plain-English:* A pathology search engine and classifier you drive with plain-English descriptions, learned from doctors' tweets.
- *Applicability:* A1, A3. Design implication: another confirmation that NL-prototype prompting of cancer images is a solved, published capability.
- *Novelty implication:* Pre-empts the NL-prompt-classify primitive for pathology.

**9. Segment Anything in Medical Images (MedSAM)** (Ma, He, Li, Han, You, Wang; Nature Communications 2024) — https://www.nature.com/articles/s41467-024-44824-z
- *Takeaway:* A promptable universal medical segmentation model — prompt-conditioned inference across modalities from a single frozen backbone.
- *Technical summary:* SAM adapted on 1.5M medical image-mask pairs yields a promptable segmentation FM where a *geometric* prompt (box) selects the target at inference on a frozen trunk; establishes the "one frozen backbone, prompt selects the task" pattern in medicine.
- *Plain-English:* One medical segmentation model that you point at what you want, no retraining per organ.
- *Applicability:* A1 (frozen trunk + prompt-conditioned inference — MORPHEUS's exact frozen-trunk ambition, but with box prompts not NL), A4 (frozen plug-in). Design implication: "frozen trunk, prompt at inference" is established; MORPHEUS's novelty is the *NL* prompt + *task inference*, not the frozen-prompt architecture.
- *Novelty implication:* Pre-empts the frozen-trunk-promptable *architecture*; leaves the NL-routing + identifiability delta intact.

**10. Capabilities of Gemini Models in Medicine (Med-Gemini)** (Saab, Tu, Weng et al., arXiv 2024) — https://arxiv.org/abs/2404.18416
- *Takeaway:* A multimodal medical generalist with web-search-augmented reasoning, prompt-specified across imaging, text, genomics, and long clinical context.
- *Technical summary:* Fine-tuned Gemini variants combine multimodal encoding with *retrieval* (uncertainty-guided web search) for medical QA/reasoning, setting SOTA on MedQA and multimodal benchmarks; explicitly blends *encoded* modalities with *retrieved* context.
- *Plain-English:* A frontier multimodal medical AI that also looks things up when unsure.
- *Applicability:* A1, A4 (an explicit *encode-vs-retrieve* production system — modalities encoded, facts retrieved). Design implication: the encode/retrieve split MORPHEUS wants to formalize is already an engineering pattern here; MORPHEUS must formalize *which molecular modality* to encode vs retrieve, a question Med-Gemini does not pose.
- *Novelty implication:* Partially pre-empts A4's spirit (hybrid encode+retrieve exists) but leaves the *molecular-modality-specific* encode-vs-RAG question open.

---

### Group B — General-purpose frozen oncology trunks (the "unified representation" substrate) (A1, A2, A4)

**11. A General-Purpose Self-Supervised Model for Computational Pathology (UNI)** (Chen et al., Nature Medicine 2024; arXiv 2308.15474) — https://arxiv.org/abs/2308.15474
- *Takeaway:* A frozen self-supervised WSI encoder (100M patches) that transfers to 33 tasks including few-shot class-prototype slide classification and 108-way subtyping.
- *Technical summary:* DINOv2 pretraining on 100k WSIs yields a frozen embedding used with lightweight probes/prototypes; strong few-shot and resolution-agnostic transfer establish the frozen-trunk + probe paradigm MORPHEUS builds on.
- *Plain-English:* A general histology feature extractor others build task heads on top of.
- *Applicability:* A1/A4 (the frozen unified trunk), A2 (but its latent is *not* identified/addressable). Design implication: UNI is the "solid unified representation" analog to MORPHEUS's TumorStateV2; the delta is *identifiability + prompting*, not the trunk.
- *Novelty implication:* Neutral — establishes the substrate, not the promptable/identified claim.

**12. A whole-slide foundation model for digital pathology from real-world data (Prov-GigaPath)** (Xu, Usuyama, Bagga et al., Nature 2024) — https://www.nature.com/articles/s41586-024-07441-w
- *Takeaway:* A slide-level (not just tile) frozen foundation model with a vision-language extension, trained on 1.3B tiles from real-world data.
- *Technical summary:* A tile encoder plus a LongNet slide aggregator produces whole-slide embeddings; a CONCH-style vision-language head enables zero-shot NL tasks. Slide-level context is the level MORPHEUS's tumor-state operates at.
- *Plain-English:* A foundation model that understands a whole biopsy slide at once, with a language interface bolted on.
- *Applicability:* A1, A2 (slide-level unified embedding), A3. Design implication: slide-level unified + NL is published; MORPHEUS's addressable-slot structure is the missing piece here too.
- *Novelty implication:* Pre-empts "slide-level unified representation with NL"; leaves identified-slots open.

**13. A foundation model for clinical-grade computational pathology and rare cancers (Virchow)** (Vorontsov, Bozkurt, Casson et al., Nature Medicine 2024) — https://www.nature.com/articles/s41591-024-03141-0
- *Takeaway:* A 632M-parameter frozen pathology FM (1.5M WSIs) with strong rare-cancer generalization from a single trunk.
- *Technical summary:* Large-scale ViT pretraining yields tile embeddings supporting biomarker prediction and pan-cancer detection via lightweight heads, emphasizing generalization to unseen/rare cancers — MORPHEUS's held-out-cancer guardrail concern.
- *Plain-English:* A very large histology model that still works on cancers it barely saw in training.
- *Applicability:* A1/A4 (frozen trunk), guardrail (rare/held-out cancer generalization). Design implication: strong generalization from a frozen trunk is achievable at scale; MORPHEUS's small-scale identifiability claims must be shown *not* to sacrifice this.
- *Novelty implication:* Neutral substrate; sets a generalization bar.

**14. Foundation model for cancer imaging biomarkers** (Pai, Bontempi, Hadzic et al., Nature Machine Intelligence 2024) — https://www.nature.com/articles/s42256-024-00807-9
- *Takeaway:* A frozen self-supervised CT foundation model whose embeddings serve multiple cancer biomarker tasks with better stability/efficiency than task-specific training.
- *Technical summary:* Contrastive pretraining on unlabeled CT lesions yields a frozen representation reused across biomarker-discovery tasks (malignancy, prognosis, molecular status), showing frozen-embedding reuse for *oncology biomarkers* specifically.
- *Plain-English:* A reusable cancer-imaging feature space that many biomarker predictors can share.
- *Applicability:* A4 (frozen-trunk reuse across oncology tasks), A2. Design implication: reinforces frozen-embedding reuse as standard; novelty must be in *how* tasks are specified (NL) and *what* the latent guarantees (identifiability).
- *Novelty implication:* Neutral substrate.

---

### Group C — NL-queryable cell/gene representations & "language of biology" (A1, A2, A3)

**15. Transfer learning enables predictions in network biology (Geneformer)** (Theodoris, Xiao, Chopra et al., Nature 2023) — https://www.nature.com/articles/s41586-023-06139-9
- *Takeaway:* A transformer pretrained on ~30M single-cell transcriptomes whose *attention* encodes gene-network structure and supports in-silico deletion (a perturbation query).
- *Technical summary:* Rank-value gene encoding + masked pretraining yields context-aware gene/cell embeddings; *in-silico* gene deletion shifts embeddings to predict disease-driving genes and candidate therapeutic targets — a counterfactual query executed on a learned representation.
- *Plain-English:* A model that learned how genes work together well enough that deleting a gene in-silico predicts real disease drivers.
- *Applicability:* A2 (gene-addressable structure), A5 (in-silico deletion = intervention-as-query), A3 (emergent network knowledge). Design implication: "delete/perturb a gene as a query on a frozen embedding" is *already* Geneformer's headline. MORPHEUS's A5 must beat it or add causal-geometry guarantees Geneformer lacks.
- *Novelty implication:* **Pre-empts A5's "perturbation as a query, not a retrained classifier."** Geneformer does exactly this for gene deletion.

**16. scGPT: toward building a foundation model for single-cell multi-omics using generative AI** (Cui, Wang, Maan et al., Nature Methods 2024) — https://www.nature.com/articles/s41592-024-02201-0
- *Takeaway:* A single-cell generative FM supporting cell-type annotation, batch integration, multi-omic integration, and *perturbation prediction* from one pretrained trunk with prompt-like conditioning.
- *Technical summary:* Transformer over (gene, expression) tokens with specialized attention; fine-tunes to diverse tasks including genetic-perturbation response prediction, and uses condition tokens as lightweight task/state prompts.
- *Plain-English:* A GPT for single cells that can annotate, integrate, and predict how cells respond to genetic changes.
- *Applicability:* A1 (condition-token prompting), A2 (gene tokens), A4 (multi-omic *encoding*), A5 (perturbation prediction). Design implication: scGPT already spans A1/A2/A4/A5 at the molecular level; MORPHEUS's addition of *histology + identifiability + NL task inference* is the differentiator, not the multi-omic promptable trunk.
- *Novelty implication:* **Broadly pre-empts the molecular half of MORPHEUS across four axes.** Strongest single "we already did the omics version" threat.

**17. Large-scale foundation model on single-cell transcriptomics (scFoundation)** (Hao, Gong, Wang et al., Nature Methods 2024) — https://www.nature.com/articles/s41592-024-02305-7
- *Takeaway:* A 100M-parameter scRNA FM whose frozen embeddings boost downstream tasks including *drug-response* and perturbation prediction.
- *Technical summary:* Read-depth-aware masked pretraining on ~50M cells yields embeddings that, used as frozen context, improve perturbation and bulk drug-response models — a frozen-trunk-plus-context pattern directly analogous to MORPHEUS A4.
- *Plain-English:* A big cell model whose features make downstream drug-response predictors better.
- *Applicability:* A4 (frozen-trunk context reuse), A5 (drug-response). Design implication: frozen molecular embeddings feeding drug-response is published; MORPHEUS's encode-vs-retrieve novelty must be a *principled decision rule*, not the reuse itself.
- *Novelty implication:* Pre-empts frozen-molecular-context reuse for drug response.

**18. Cell2Sentence: Teaching Large Language Models the Language of Biology** (Levine, Rizvi, Lévy et al., ICML 2024; bioRxiv 2023.09.11.557287) — https://www.biorxiv.org/content/10.1101/2023.09.11.557287
- *Takeaway:* Represent a cell as a *sentence* of ranked gene names so a general LLM can read/generate/query cells in natural language.
- *Technical summary:* Cells are serialized to gene-name sentences; fine-tuning an LLM on these enables cell-type prediction, cell generation from text, and NL interrogation of transcriptomic state — literally an NL⇄cell-state interface. C2S-Scale later shows scaling laws and emergent capabilities.
- *Plain-English:* Turn each cell into a sentence and let a language model talk about, and generate, cells.
- *Applicability:* A1/A3 (NL⇄cell-state grounding is the whole method), A2. Design implication: "NL-queryable tumor/cell representation" in the literal sense already exists. MORPHEUS's grounding must add *identifiability + multimodal (histology) + causal query* that a gene-name sentence cannot express.
- *Novelty implication:* **Directly pre-empts "NL-queryable tumor representation."** MORPHEUS cannot claim first NL interface to cell state.

**19. GenePT: A Simple But Hard-to-Beat Foundation Model for Genes and Cells Built from ChatGPT** (Chen & Zou, bioRxiv 2023/2024) — https://www.biorxiv.org/content/10.1101/2023.10.16.562533
- *Takeaway:* Embedding genes via their NCBI text descriptions (GPT-3.5 embeddings) and cells as weighted gene-text averages rivals specialized single-cell FMs — with essentially no training.
- *Technical summary:* Gene = LLM embedding of its text description; cell = expression-weighted average of gene embeddings; matches or beats Geneformer/scGPT on several tasks. A devastating "steelman-the-null" for expensive omics FMs.
- *Plain-English:* Just using ChatGPT's text embeddings of gene descriptions works about as well as billion-cell models.
- *Applicability:* A3 (NL grounding *is* the representation), A4 (a "retrieve the text, don't encode the omics" existence proof). Design implication: this is a live null hypothesis MORPHEUS must beat — if a text-embedding baseline matches its encoded representation, the "encode" side of A4 is unjustified for that modality.
- *Novelty implication:* **Steelman-the-null flagship.** Forces MORPHEUS to demonstrate its encoded latent beats a cheap NL/RAG baseline; otherwise A4's "encode" recommendation is unsupported.

**20. Large language models encode clinical knowledge (Med-PaLM)** (Singhal, Azizi, Tu et al., Nature 2023) — https://www.nature.com/articles/s41586-023-06291-2
- *Takeaway:* Prompting/instruction-tuning elicits latent clinical knowledge from an LLM, evaluated with a purpose-built human-axes framework (MultiMedQA + rubric).
- *Technical summary:* Introduces instruction prompt-tuning and a multi-axis human evaluation (factuality, harm, bias) to *measure elicited* medical knowledge, reaching passing-level MedQA. The evaluation *methodology* — measuring elicited knowledge, not downstream accuracy — is what MORPHEUS A3 wants for biology.
- *Plain-English:* A framework for measuring how much real medical knowledge you can coax out of a language model by asking well.
- *Applicability:* A3 (emergent/elicited-knowledge *evaluation* methodology). Design implication: borrow the multi-axis elicitation-evaluation design, but for *biological* (pathway/mechanism) knowledge — an under-served niche MORPHEUS can own.
- *Novelty implication:* Pre-empts the *elicitation-evaluation* concept in medicine; MORPHEUS's A3 novelty is porting it to *emergent biological/mechanistic* knowledge, which is not yet standardized.

---

### Group D — Adversarial evaluations / steelman-the-null on emergent biology (A3, A2)

**21. Assessing the limits of zero-shot foundation models in single-cell biology** (Kedzierska, Crawford, Amini, Lu; bioRxiv 2023) — https://www.biorxiv.org/content/10.1101/2023.10.16.561085
- *Takeaway:* Zero-shot embeddings from Geneformer and scGPT often *fail to beat* simple baselines on clustering and batch integration.
- *Technical summary:* Systematic zero-shot evaluation shows the pretrained embeddings underperform highly-variable-gene + PCA and other simple pipelines on core tasks, questioning claims of emergent biological structure without fine-tuning.
- *Plain-English:* Out of the box, these big cell models often lose to a basic PCA.
- *Applicability:* A3 (directly measures "emergent biological knowledge" and finds it lacking), A2. Design implication: MORPHEUS *must* pre-register a zero-shot emergence benchmark with strong classical baselines, or its representation-quality claims will be dismissed exactly as these were.
- *Novelty implication:* **Reframes A3 as an open problem** — the field lacks a credible emergence metric and current FMs fail naive ones. This is MORPHEUS's opening *and* its bar.

**22. A Deep Dive into Single-Cell RNA Sequencing Foundation Models** (Boiarsky, Singh, Buendia, Getz, Sontag; bioRxiv 2023) — https://www.biorxiv.org/content/10.1101/2023.10.19.563100
- *Takeaway:* Simple logistic regression on standard features matches scGPT/Geneformer on cell-type annotation, challenging the value-add of pretraining.
- *Technical summary:* Controlled comparison shows a well-tuned linear classifier is competitive with FM embeddings for annotation, and pretraining benefits are marginal for the tested tasks.
- *Plain-English:* For labeling cell types, plain logistic regression basically ties the fancy models.
- *Applicability:* A2/A3 (null baseline for "identified/addressable" adding value). Design implication: MORPHEUS must show identifiability/addressability yields *measurable transfer/reliability gains over linear baselines* — the whole point of A2 — not just competitive accuracy.
- *Novelty implication:* Steelman-the-null for A2. Sets the "beat logistic regression meaningfully" bar.

**23. Deep learning-based predictions of gene perturbation effects do not yet outperform simple linear baselines** (Ahlmann-Eltze, Huber, Anders; bioRxiv 2024) — https://www.biorxiv.org/content/10.1101/2024.09.16.613342
- *Takeaway:* State-of-the-art perturbation predictors (incl. GEARS-class and FM-based) do *not* beat a simple additive/linear model at predicting unseen perturbation effects.
- *Technical summary:* Rigorous re-benchmark with proper baselines and metrics shows deep perturbation-response models fail to exceed a linear "mean-shift" baseline on held-out perturbations, exposing evaluation leakage/optimism in the sub-field.
- *Plain-English:* For predicting what a genetic tweak does, today's deep models don't beat a straight-line guess.
- *Applicability:* A5 (the causal/interventional-query axis's null hypothesis). Design implication: MORPHEUS's interventional-query claims require an evaluation that a linear mean-shift *cannot* pass — otherwise A5 is unfalsifiable hype.
- *Novelty implication:* **Steelman-the-null flagship for A5.** The entire perturbation-prediction premise is under credible attack; MORPHEUS's causal-geometry framing must clear this bar to be worth anything.

**24. Evaluating the Utilities of Foundation Models in Single-cell Data Analysis (scEval)** (Liu, Chen, Wang et al., bioRxiv 2023) — https://www.biorxiv.org/content/10.1101/2023.09.08.555192
- *Takeaway:* A multi-task benchmark quantifying where single-cell FMs help vs. where they don't (integration, annotation, perturbation, imputation).
- *Technical summary:* Systematic evaluation across scGPT/Geneformer/scFoundation and tasks establishes task-conditional utility and failure modes, a template for honest FM assessment.
- *Plain-English:* A scoreboard showing which jobs cell foundation models are actually good at.
- *Applicability:* A3 (evaluation infrastructure), A2. Design implication: reuse/extend this harness; do not invent an unvetted benchmark of MORPHEUS's own that lacks external baselines.
- *Novelty implication:* Reframes A3 — the evaluation scaffolding partly exists; MORPHEUS should add the *emergence/elicitation* dimension it omits.

---

### Group E — Perturbation / drug / intervention-as-a-query systems (A5, A2)

**25. scGen predicts single-cell perturbation responses** (Lotfollahi, Wolf, Theis; Nature Methods 2019) — https://www.nature.com/articles/s41592-019-0494-8
- *Takeaway:* Latent-space *vector arithmetic* predicts a cell's response to a perturbation it never saw — an intervention applied *after* encoding.
- *Technical summary:* A VAE learns a latent where a perturbation is a difference vector; adding it to control cells generates predicted perturbed states across unseen cell types. This is precisely "drug/perturbation = a spec applied post-encoding," MORPHEUS's A5 design.
- *Plain-English:* Learn the "direction" a treatment moves cells in latent space, then apply that arrow to new cells to predict their response.
- *Applicability:* A5 (intervention-as-latent-operation, MORPHEUS's exact framing), A2. Design implication: the "apply the intervention after encoding" mechanism is 2019 prior art; MORPHEUS must claim *causal-geometry validity* (geodesic ≈ causal distance) that scGen's Euclidean arithmetic lacks.
- *Novelty implication:* **Pre-empts the mechanism of A5.** Post-encoding intervention arithmetic is established; MORPHEUS's only opening is the *geometry/identifiability guarantee*.

**26. Predicting cellular responses to complex perturbations (Compositional Perturbation Autoencoder, CPA)** (Lotfollahi, Klimovskaia Susmelj, De Donno et al., Molecular Systems Biology 2023) — https://www.embopress.org/doi/full/10.15252/msb.202211517
- *Takeaway:* A disentangled latent with *addressable*, composable perturbation/dose/covariate embeddings enabling counterfactual response prediction to unseen drug/dose/combination.
- *Technical summary:* CPA factorizes cell state into basal + additive perturbation + covariate latents; because perturbations are *separately addressed* and linearly composable, it answers counterfactual "what if drug A at dose d in cell type c" queries and unseen combinations.
- *Plain-English:* Splits a cell's state into "baseline" plus "drug effect" plus "context," so you can mix and match to predict new treatment scenarios.
- *Applicability:* A2 (addressable, composable slots), A5 (counterfactual/compositional intervention queries). Design implication: MORPHEUS's "identified, pathway-addressable, compositional slots for counterfactual queries" is *substantially* CPA's design — MORPHEUS must move addressability from *perturbation labels* to *biological pathway programmes* and add histology/NL.
- *Novelty implication:* **Pre-empts A2+A5 jointly** (addressable + compositional + counterfactual). One of the two most dangerous prior-art papers for MORPHEUS's core thesis.

**27. Predicting cellular responses to novel drug perturbations at single-cell resolution (chemCPA)** (Hetzel, Böhm, Kilbertus, Günnemann, Lotfollahi, Theis; NeurIPS 2022; arXiv 2204.13545) — https://arxiv.org/abs/2204.13545
- *Takeaway:* Extends CPA to *unseen drugs* via molecular-structure conditioning + transfer from bulk data.
- *Technical summary:* An encoder-decoder with a chemical-structure encoder and architecture-surgery transfer predicts single-cell responses to drugs absent from training, enabling in-silico screening as a query.
- *Plain-English:* Predict how cells react to a brand-new drug from its chemical structure, without ever having tested it.
- *Applicability:* A5 (drug-as-query generalizing to unseen molecules). Design implication: generalization of intervention queries to unseen drugs is solved for transcriptomics; MORPHEUS's contribution must be the *multimodal/frozen-trunk* setting, not the query capability.
- *Novelty implication:* Pre-empts "unseen-drug intervention query."

**28. Predicting transcriptional outcomes of novel multigene perturbations with GEARS** (Roohani, Huang, Leskovec; Nature Biotechnology 2024) — https://www.nature.com/articles/s41587-023-01905-6
- *Takeaway:* A GNN over a gene-gene knowledge graph predicts responses to *unseen* single and *combinatorial* genetic perturbations, including non-additive genetic interactions.
- *Technical summary:* GEARS couples perturbation embeddings with a prior gene-relationship graph to extrapolate to novel multi-gene perturbations and flag epistatic (synergistic/suppressive) interactions — interventional queries over combinatorial space.
- *Plain-English:* Predicts what happens when you knock out gene combinations you've never tested, including surprising interaction effects.
- *Applicability:* A5 (combinatorial intervention-as-query), A2 (graph-structured gene addressing). Design implication: combinatorial counterfactual queries with interaction structure exist; MORPHEUS's "causal geometry" must add something beyond a knowledge-graph prior.
- *Novelty implication:* Pre-empts combinatorial interventional queries (though note: Ahlmann-Eltze #23 disputes GEARS-class gains over linear baselines — cite both).

**29. Learning single-cell perturbation responses using neural optimal transport (CellOT)** (Bunne, Stark, Gut et al.; Nature Methods 2023) — https://www.nature.com/articles/s41592-023-01969-x
- *Takeaway:* Neural optimal transport learns the *map* from control to perturbed cell distributions, generalizing to held-out patients/cells.
- *Technical summary:* Learns an OT map (a transport plan) representing a perturbation's distributional effect, predicting unseen single-cell responses and out-of-sample patients — a principled, geometry-flavored intervention operator.
- *Plain-English:* Learns how a treatment "transports" the whole cloud of cells from untreated to treated, then applies it to new samples.
- *Applicability:* A5 (intervention as a learned transport map — close cousin of MORPHEUS's "geodesic causal" idea), causal-geometry. Design implication: OT gives a *distributional/geometric* intervention operator already; MORPHEUS's Riemannian-geodesic framing overlaps and must be differentiated (e.g., identifiability-conditioned geodesics).
- *Novelty implication:* **Pre-empts the "geometric intervention operator" spirit of A5.** Together with CPA and scGen, the causal-geometry axis is crowded.

**30. Disentanglement of single-cell data with biolord** (Piran, Cohen, Hoshen, Nitzan; Nature Biotechnology 2024) — https://www.nature.com/articles/s41587-023-02079-x
- *Takeaway:* Explicitly *disentangles* known + unknown attributes so each is separately controllable, enabling counterfactual generation over held-out attribute combinations.
- *Technical summary:* biolord partitions the latent into labeled biological attributes and residual unknowns with decorrelation, giving *addressable* factors you can set independently to generate counterfactual cell states (e.g., unseen cell-type × perturbation × time).
- *Plain-English:* Separates a cell's traits into independent dials you can turn one at a time to imagine new combinations.
- *Applicability:* A2 (identified, separately-addressable attribute slots — MORPHEUS's exact A2 goal), A5 (counterfactual generation). Design implication: attribute-addressable disentanglement for counterfactuals *in cancer-relevant single-cell data* is published; MORPHEUS's differentiator is *pathway-programme* addressing + multimodal + NL, not the disentangled-slot idea.
- *Novelty implication:* **Pre-empts A2's "identified, addressable slots for counterfactual queries."** With CPA, the strongest A2 pre-emption.

---

### Group F — Multimodal oncology integration & encode-vs-retrieve (A4, A2)

**31. Modeling Dense Multimodal Interactions Between Biological Pathways and Histology for Survival Prediction (SurvPath)** (Jaume, Vaidya, Chen, Williamson, Mahmood; CVPR 2024; arXiv 2304.06819) — https://arxiv.org/abs/2304.06819
- *Takeaway:* Tokenizes transcriptomics into **biological pathway tokens** and models their dense interactions with histology tokens in a multimodal transformer — pathway-*addressable* multimodal tumor tokens already exist.
- *Technical summary:* Transcriptomics is grouped into pathway tokens (named, biologically-defined units); a memory-efficient transformer computes pathway×histology co-attention for survival, and the learned attention exposes which pathways interact with which morphology.
- *Plain-English:* Groups genes into named biological pathways as tokens, then lets the model see how each pathway "talks to" the tumor's appearance.
- *Applicability:* A2 (pathway-addressable multimodal slots — the literal thing MORPHEUS's unwired TQI wants), A4 (transcriptomics *encoded* as pathway tokens fused with histology). Design implication: MORPHEUS's "pathway slots" are *architecturally* pre-empted by SurvPath; the delta is that SurvPath's pathway tokens are *defined a priori by gene sets*, NOT *identified* (iVAE-style) and NOT NL-promptable/counterfactual. MORPHEUS must show *identifiability* buys reliability SurvPath's fixed gene-set tokens don't.
- *Novelty implication:* **Pre-empts the *architecture* of A2** (pathway-addressable multimodal tokens). MORPHEUS's only defensible A2 novelty is *identified* (vs. gene-set-defined) slots + prompting/causal use.

**32. Pan-cancer integrative histology-genomic analysis via multimodal deep learning (PORPOISE)** (Chen, Lu, Williamson et al.; Cancer Cell 2022) — https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00317-8
- *Takeaway:* A pan-cancer multimodal (WSI + molecular) survival model with interpretability linking morphology and molecular features across 14 cancer types.
- *Technical summary:* Fuses WSI and genomic/transcriptomic features for prognosis across the TCGA pan-cancer cohort and extracts multimodal prognostic markers — the canonical encode-both-modalities oncology baseline.
- *Plain-English:* Combines slide images and genetics across many cancers to predict outcome and explain which features matter.
- *Applicability:* A4 (both modalities *encoded and fused* — the default MORPHEUS compares against), A2. Design implication: full-encode fusion is the standard; A4's contribution is a *decision rule* for when NOT to encode (RAG instead), which PORPOISE never asks.
- *Novelty implication:* Neutral baseline; defines the "encode everything" default A4 must improve on.

**33. Harnessing multimodal data integration to advance precision oncology** (Boehm, Khosravi, Vanguri, Gao, Shah; Nature Reviews Cancer 2022) — https://www.nature.com/articles/s41568-021-00408-3
- *Takeaway:* A structured review of *how* to combine imaging, pathology, molecular, and clinical data (early/late/joint fusion), i.e., the encode-vs-integrate design space.
- *Technical summary:* Surveys fusion strategies, missing-modality handling, and evaluation pitfalls for multimodal oncology, framing the trade-offs MORPHEUS A4 wants to formalize.
- *Plain-English:* A field guide to the ways you can merge different cancer data types, and their pros and cons.
- *Applicability:* A4 (the encode/integrate taxonomy), guardrail (missing-modality). Design implication: use its taxonomy to *position* the encode-vs-retrieve question; note that "retrieve as context (RAG)" is largely *absent* from this 2022 taxonomy — a genuine gap MORPHEUS can occupy.
- *Novelty implication:* Reframes A4 — the encode/fuse space is mapped, but *RAG-as-a-modality-choice* is not, leaving MORPHEUS's specific formalization open.

**34. Almanac — Retrieval-Augmented Language Models for Clinical Medicine** (Zakka, Shad, Chaurasia et al.; NEJM AI 2024) — https://ai.nejm.org/doi/full/10.1056/AIoa2300068
- *Takeaway:* Retrieval augmentation improves factuality/safety of clinical LLM answers vs. parametric-only models — the "treat knowledge as retrieved context" side of encode-vs-retrieve.
- *Technical summary:* Grounds clinical LLM responses in retrieved guideline/reference text, improving factuality, completeness, and safety over non-retrieval baselines in physician evaluation.
- *Plain-English:* A medical chatbot that looks up trusted references before answering, and is safer for it.
- *Applicability:* A4 (RAG-as-context existence proof and safety argument — mirrors MORPHEUS's closed-RAG card), A3. Design implication: supports keeping *sparse, safety-sensitive, verbal* knowledge as RAG; the open question MORPHEUS uniquely poses is when a *quantitative molecular modality* is better retrieved than encoded — not addressed here.
- *Novelty implication:* Pre-empts RAG-as-context in the *textual/clinical* setting; leaves the *molecular-modality* encode-vs-retrieve question genuinely open.

---

### Group G — Causal / interventional representation learning (A5, A2)

**35. Learning Causal Representations of Single Cells via Sparse Mechanism Shift Modeling (sVAE+)** (Lopez, Tagasovska, Ra, Cho, Pritchard, Regev; CLeaR/arXiv 2022) — https://arxiv.org/abs/2211.03553
- *Takeaway:* Treats each genetic perturbation as a *sparse mechanism shift*, giving an *identifiable* causal latent whose dimensions correspond to perturbation-affected mechanisms.
- *Technical summary:* Under sparse-mechanism-shift assumptions, interventions (Perturb-seq targets) identify latent factors up to permutation; the model recovers which latent dimensions each perturbation acts on — perturbation-conditioned identifiability, the exact iVAE-flavored claim MORPHEUS A2 rests on.
- *Plain-English:* Because each genetic knockout changes only a few hidden mechanisms, you can pin down what those mechanisms are and which ones each knockout hits.
- *Applicability:* A2 (perturbation-conditioned identifiability of addressable mechanisms), A5. Design implication: MORPHEUS's "perturbation-conditioned iVAE → latent dims = real pathway states" is *this paper's* thesis in single-cell. MORPHEUS must claim the *multimodal + histology + NL-prompting* extension, not the identification-by-sparse-shift principle.
- *Novelty implication:* **Pre-empts A2's identification mechanism** in exactly the perturbation-conditioned form MORPHEUS proposed.

**36. Combinatorial prediction of therapeutic perturbations using causally-inspired neural networks (PDGrapher)** (Gonzalez, Ricci-Tam, Zitnik et al.; bioRxiv 2024) — https://www.biorxiv.org/content/10.1101/2024.01.03.573985
- *Takeaway:* Predicts *which* interventions (gene targets) shift a diseased cell to a healthy state — inverts the usual response-prediction into a causal *target-discovery* query.
- *Technical summary:* A causally-inspired GNN over protein-interaction/gene-regulatory graphs predicts combinatorial perturbagens driving a desired state transition, framed explicitly as an intervention-planning (not correlation) task.
- *Plain-English:* Instead of asking "what does this drug do," it asks "what should we perturb to fix this cell," and answers with gene combinations.
- *Applicability:* A5 (interventional *target-discovery* query — arguably a harder counterfactual than MORPHEUS's), A2. Design implication: the "intervention as a query to plan, not classify" framing is published; MORPHEUS's causal-geometry must offer something PDGrapher's graph approach doesn't (e.g., frozen-trunk, NL-specified goals).
- *Novelty implication:* Pre-empts "intervention planning as a query." Strengthens the case that A5 as a *concept* is not novel; only a specific instantiation could be.

---

## Synthesis for the fleet

**Where MORPHEUS is most pre-empted (kill-or-reframe candidates):**
- **A1 (promptable, NL-task-inferring multimodal cancer model):** essentially done — Med-PaLM M (#1), PathChat (#5, oncology), BiomedParse (#7), scGPT condition-prompting (#16). Reframe A1 around *frozen + identified + causal*, not "you can prompt it."
- **A2 (identified, pathway-addressable slots):** the *architecture* is SurvPath (#31, pathway tokens) and the *identification/disentanglement* is CPA (#26), biolord (#30), sVAE+ (#35). MORPHEUS's only opening: **identified** pathway slots (vs. a priori gene-set tokens) whose identifiability *measurably improves prompt reliability/transfer* — and no prior work has shown that link.
- **A5 (intervention/drug-as-query):** crowded and mature — scGen latent arithmetic (#25), CPA/chemCPA (#26/#27), GEARS (#28), CellOT geometric operator (#29), Geneformer in-silico deletion (#15). And #23 argues the whole sub-field doesn't beat linear baselines. MORPHEUS must clear the *linear-mean-shift null* and show *causal-geometry* (geodesic≈causal) adds calibrated value — a claim no cited paper validates.

**Where a genuine gap survives (defensible novelty):**
- **A3 emergence/elicitation evaluation:** the critiques (#21 Kedzierska, #22 Boiarsky, #23 Ahlmann-Eltze) show the field *lacks* a credible emergent-biological-knowledge metric and current FMs fail naive ones; the elicitation-evaluation methodology exists only for clinical text (#20 Med-PaLM). A rigorous *emergent-mechanistic-knowledge* benchmark for a multimodal cancer model is genuinely unoccupied.
- **A4 molecular encode-vs-retrieve decision rule:** the fusion taxonomy (#33) and RAG-as-context (#34, GenePT #19) both exist, but *no* cited work formalizes *when a quantitative molecular modality (proteomics/phospho/CNV/SNV) should be encoded vs. retrieved as context* on a frozen trunk. GenePT (#19) makes this urgent by showing a text/RAG baseline can match encoded omics.
- **The synthesis itself:** no single cited system combines a *frozen, identified* oncology trunk + NL task routing + pathway-addressable identified slots + post-encoding causal-geometry query + an emergence benchmark. Each piece is pre-empted; the *integrated, identifiability-load-bearing* whole is not.
