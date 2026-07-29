## Promptable / task-general unified interfaces + task auto-detection

Lane remit: generalist/unified models with a single promptable interface that *infers* the task (Gato, Painter/SegGPT, Unified-IO, in-context task inference, prompt-conditioned routing / mixture-of-experts). This is the A1 core: one representation + task auto-detection instead of hard-coded probes. Entries also tagged where they inform A2-A5.

MORPHEUS rebase axes: **A1** promptable unified rep + NL task auto-detection; **A2** identified, pathway-addressable slots making prompting reliable; **A3** NL<->biology grounding + emergent-knowledge elicitation & its evaluation; **A4** multimodal prompting (encode vs RAG; frozen-trunk plug-in); **A5** interventional/causal-geometry queries (perturbation/drug as a query).

---

### 1. A Generalist Agent (Gato)
**Reed et al., DeepMind — TMLR 2022.** arXiv:2205.06175
- **Takeaway:** One transformer with one set of weights performs 604 tasks across text, images, Atari, and robot control by inferring the task from context.
- **Technical summary:** All modalities/tasks are serialized into a flat token sequence and modeled autoregressively by a single 1.18B-param transformer. At inference Gato decides from the prompt/context whether to emit text, button presses, or joint torques — the task is never explicitly selected. It exceeds expert performance on ~450 of 604 tasks with no per-task head.
- **Plain-English:** A single AI brain plays games, chats, captions images, and drives a robot arm, working out which job to do just from what it's shown.
- **Applicability (A1):** Canonical proof that a shared token stream + context-conditioned decoding can replace task-specific heads — the exact substrate MORPHEUS wants for auto-detecting a biological query. Design implication: commit to a single serialized interface where the *prompt selects the behavior*, not a routing classifier.
- **Novelty implication:** **Pre-empts** any generic claim that "one promptable model spanning many task types" is itself novel — MORPHEUS must locate novelty in the *biological* axes (A2 identifiability, A5 causal queries), not in the generalist framing.

### 2. Images Speak in Images: A Generalist Painter for In-Context Visual Learning
**Wang et al., BAAI — CVPR 2023.** arXiv:2212.02499
- **Takeaway:** Vision tasks are redefined as image-to-image, and the task is specified by an input/output *example pair* — a visual prompt — rather than a task label.
- **Technical summary:** Trained by masked image modeling on stitched input-output pairs; at inference a demonstration pair conditions which task to perform (depth, segmentation, denoising, keypoints). It matches task-specific models on 7 tasks and generalizes to out-of-domain tasks (open-category keypoints) never seen in training.
- **Plain-English:** Show the model one worked example ("here's the input, here's the answer") and it repeats that operation on new images — no retraining, no task menu.
- **Applicability (A1, A3):** Demonstrates task auto-detection purely from a demonstration, the tightest form of "the prompt is the task specification." For MORPHEUS: a canonical (perturbed-state, response) example pair could itself define a query in-context.
- **Novelty implication:** **Reframes** — pushes MORPHEUS to distinguish *NL* task auto-detection (its A1) from *demonstration-pair* task specification; NL grounding (A3) is the differentiator, since Painter has no language channel.

### 3. SegGPT: Segmenting Everything In Context
**Wang et al., BAAI — ICCV 2023.** arXiv:2304.03284
- **Takeaway:** Unifies all segmentation variants (instance, part, contour, text, video) into one in-context "coloring" problem solved from a prompt example.
- **Technical summary:** Formulated as in-context coloring with random per-sample color mapping so the model must infer the intended grouping from context rather than memorizing fixed colors. One model does few-shot semantic seg, video object seg, panoptic seg via in-context inference, in- and out-of-domain.
- **Plain-English:** Give it an example of what to outline and it outlines the same kind of thing anywhere, even in videos it wasn't trained on.
- **Applicability (A1, A2):** The random-color-map trick forces the model to bind meaning to *context*, not to a fixed output slot — relevant to A2's goal of addressable-but-not-hardcoded programme slots.
- **Novelty implication:** **Strengthens** the feasibility of "same trunk, many task granularities via prompt," but again shows the vision community owns the generalist framing; MORPHEUS novelty is the biology-specific addressability.

### 4. Unified-IO: A Unified Model for Vision, Language, and Multi-Modal Tasks
**Lu et al., AI2 — ICLR 2023.** arXiv:2206.08916
- **Takeaway:** Homogenizes every input/output (images, masks, boxes, depth, text) into one discrete token vocabulary so a single seq2seq transformer does 90+ tasks with no task-specific heads.
- **Technical summary:** Per-pixel maps, boxes, and language are all discretized into a shared vocabulary; one encoder-decoder is jointly trained across 90+ datasets. First model to solve all 7 GRIT tasks and competitive on 16 benchmarks with no task-specific fine-tuning.
- **Plain-English:** By turning pictures, outlines, and words into the same kind of "tokens," one model can be asked to do almost any vision-or-language job.
- **Applicability (A1, A4):** The tokenization-of-heterogeneous-modalities strategy directly informs A4 (how to *encode* a modality into the shared stream vs keep it external). Design implication: decide MORPHEUS's discretization boundary — which omics get tokenized into the trunk vs supplied as RAG context.
- **Novelty implication:** **Pre-empts** "unified heterogeneous I/O" as a novelty; sharpens MORPHEUS's A4 claim toward the *encode-vs-retrieve decision rule*, which Unified-IO does not address (it encodes everything).

### 5. Unified-IO 2: Scaling Autoregressive Multimodal Models with Vision, Language, Audio, and Action
**Lu et al., AI2 — CVPR 2024.** arXiv:2312.17172
- **Takeaway:** Extends the unified-token idea to text+image+audio+action generation in one autoregressive model trained from scratch with a multimodal mixture-of-denoisers.
- **Technical summary:** Tokenizes images, text, audio, action, and boxes into a shared semantic space processed by one encoder-decoder; instruction-tuned on 120 datasets. SOTA on GRIT and strong across 30+ benchmarks, including generation and understanding, from a single set of weights.
- **Plain-English:** A bigger single model that can both understand and *produce* pictures, sound, text, and robot actions, choosing what to output from the instruction.
- **Applicability (A1, A4):** Shows instruction-conditioned generation across many modalities is trainable at scale — supports a MORPHEUS trunk that both reads and *writes* biological states (e.g., predicted expression) under NL control.
- **Novelty implication:** **Strengthens** feasibility of generation-as-query but again **pre-empts** generic multimodal-generalist novelty.

### 6. OFA: Unifying Architectures, Tasks, and Modalities Through a Simple Seq2Seq Framework
**Wang et al., Alibaba DAMO — ICML 2022.** arXiv:2202.03052
- **Takeaway:** Uses natural-language *instructions* (not task IDs) inside a single seq2seq model to specify each of a dozen cross-modal tasks, in both pretraining and fine-tuning.
- **Technical summary:** Handcrafted instruction templates ("what does the image describe?", "which region does X refer to?") route behavior through one Transformer encoder-decoder over a shared vocabulary. Pretrained on only 20M pairs, it sets SOTA on multiple cross-modal tasks with no task-specific layers.
- **Plain-English:** You tell one model what to do in plain words, and it captions, detects, generates, or grounds — no separate networks.
- **Applicability (A1):** This is the closest classical analog to MORPHEUS's A1: *NL instruction = task selector* inside a unified trunk. Design implication: adopt an instruction schema as the front door, but MORPHEUS must add *auto-detection* (OFA still relies on fixed template phrasing).
- **Novelty implication:** **Pre-empts** "NL instruction selects the task in a unified model." MORPHEUS's defensible delta is inferring an *unspecified/ambiguous* biological task and routing to identified programmes (A2), not template matching.

### 7. Flamingo: a Visual Language Model for Few-Shot Learning
**Alayrac et al., DeepMind — NeurIPS 2022.** arXiv:2204.14198
- **Takeaway:** Bridges frozen vision and frozen language backbones with trainable cross-attention, enabling few-shot task adaptation from interleaved image-text prompts.
- **Technical summary:** Gated cross-attention layers ("Perceiver Resampler" + adapters) let a frozen LLM attend to visual features; trained on interleaved web corpora. New tasks (VQA, captioning) are learned in-context from a handful of examples with no weight updates.
- **Plain-English:** Keep two big pretrained models frozen, glue them with a small trainable bridge, and the combined system learns new picture-and-text tasks from a few examples.
- **Applicability (A4, A1):** The **frozen-trunk plug-in** pattern is exactly MORPHEUS's A4 thesis: freeze a pretrained biological trunk, attach lightweight adapters to bring in a new modality. Design implication: prefer gated cross-attention adapters over full retraining when adding proteomics/phospho.
- **Novelty implication:** **Strengthens** A4's frozen-trunk plug-in as established practice — lowers MORPHEUS's technical risk but means the plug-in mechanism itself is prior art; novelty is *which biological modality to encode vs retrieve*.

### 8. Segment Anything (SAM)
**Kirillov et al., Meta — ICCV 2023.** arXiv:2304.02643
- **Takeaway:** A promptable segmentation foundation model: points, boxes, or masks are prompts, and the model returns valid masks with zero-shot generalization.
- **Technical summary:** ViT image encoder computes one embedding; a lightweight prompt encoder + mask decoder resolve arbitrary prompts in real time. Trained on SA-1B (1B masks/11M images) via a promptable-segmentation objective that handles ambiguity by emitting multiple valid masks.
- **Plain-English:** Click or draw a box on anything in an image and it cuts out the object, even categories it never explicitly learned.
- **Applicability (A1, A2):** SAM's explicit **ambiguity-aware** promptable design (return several valid interpretations) is a template for MORPHEUS handling under-specified biological queries. Design implication: model prompt ambiguity as multiple outputs rather than forcing a single answer.
- **Novelty implication:** **Reframes** — establishes "promptable foundation model" as a design pattern MORPHEUS inherits; the biological analog (promptable pathway addressing) is the open, novel territory.

### 9. SAM 2: Segment Anything in Images and Videos
**Ravi et al., Meta — 2024.** arXiv:2408.00714
- **Takeaway:** Extends promptable segmentation to video with a streaming memory, so a prompt propagates through time.
- **Technical summary:** Adds a memory bank/attention over prior frames to the SAM architecture, enabling promptable video object segmentation with real-time interaction and correction. Trained on a large video mask dataset collected in the loop.
- **Plain-English:** Point at an object once and it tracks and segments it through a whole video.
- **Applicability (A1, A5):** Memory-conditioned prompting over a temporal axis is analogous to conditioning a biological query on a trajectory (e.g., time-course after perturbation). Informs A5 where a query implies a temporal/counterfactual rollout.
- **Novelty implication:** **Strengthens** the idea that a single prompt can drive a *dynamical* prediction, supporting MORPHEUS's A5 framing of perturbation-as-query-over-a-trajectory.

### 10. Sequential Modeling Enables Scalable Learning for Large Vision Models (LVM)
**Bai et al., UC Berkeley/JHU — CVPR 2024.** arXiv:2312.00785
- **Takeaway:** Trains a vision model purely by next-token prediction over "visual sentences," with test-time behavior selected by visual prompts — no language at all.
- **Technical summary:** Converts images, videos, and annotations into 420B visual tokens; a single next-token objective yields a model whose task is set by a carefully constructed visual prompt at inference. Demonstrates scaling with data diversity and prompt-defined task solving.
- **Plain-English:** Feed a model billions of "sentences made of image pieces," and afterward you steer it to new tasks just by arranging example images as a prompt.
- **Applicability (A1):** Pure-sequence, language-free task specification is an existence proof that the *prompt structure* alone can carry the task — relevant if MORPHEUS supports non-NL (data-only) query modes alongside NL.
- **Novelty implication:** **Reframes** MORPHEUS's A1: NL is one prompt channel among several; the claim should be that MORPHEUS *fuses* NL task auto-detection with data-shaped prompts, which LVM (no NL) does not.

### 11. Pix2seq: A Language Modeling Framework for Object Detection
**Chen et al., Google — ICLR 2022.** arXiv:2109.10852
- **Takeaway:** Recasts a structured perception task (detection) as plain sequence generation, showing "read out what you know as tokens" beats bespoke heads.
- **Technical summary:** Bounding boxes + classes are discretized into a token sequence and produced by a standard encoder-decoder with only data augmentation as task-specific machinery. Competitive with specialized detectors on COCO.
- **Plain-English:** Instead of a custom detection module, just have the model *describe* the boxes as a sentence of numbers.
- **Applicability (A1, A2):** Foundational argument that specialized probes/heads can be replaced by generation over an addressable token grammar — directly supports MORPHEUS replacing hard-coded biological probes with a generated, addressable output grammar (A2).
- **Novelty implication:** **Strengthens** the "no hard-coded probes" thesis with an early, clean precedent; MORPHEUS should cite it as the lineage for probe-free task read-out.

### 12. Finetuned Language Models Are Zero-Shot Learners (FLAN)
**Wei et al., Google — ICLR 2022.** arXiv:2109.01652
- **Takeaway:** Instruction-tuning across many tasks phrased as NL instructions makes a model generalize zero-shot to unseen task *types*.
- **Technical summary:** 137B LM fine-tuned on 60+ NLP tasks verbalized as instruction templates; outperforms zero-shot (and some few-shot) GPT-3 on 20/25 datasets. Ablations show NL instructions and task *diversity* are essential — removing instructions collapses generalization.
- **Plain-English:** Teaching a model to follow many worded instructions makes it able to follow brand-new instructions it never saw.
- **Applicability (A1, A3):** Establishes that NL instructions are the mechanism enabling task auto-detection to *unseen* tasks — the theoretical backbone of MORPHEUS A1. Design implication: breadth and NL phrasing of the training task suite governs zero-shot biological task coverage.
- **Novelty implication:** **Pre-empts** the general claim "NL instructions unlock unseen-task generalization"; MORPHEUS must show this transfers to *biological* task inference and that emergent biological knowledge (A3) is measurable, not assumed.

### 13. Multitask Prompted Training Enables Zero-Shot Task Generalization (T0)
**Sanh et al. (BigScience) — ICLR 2022.** arXiv:2110.08207
- **Takeaway:** Explicit multitask training on prompted datasets induces zero-shot generalization in a model 16x smaller than GPT-3.
- **Technical summary:** Converts many supervised datasets into multiple prompt templates and fine-tunes a T5 encoder-decoder; the resulting T0 generalizes to held-out task types, often beating far larger models on unseen tasks and BIG-bench.
- **Plain-English:** Train on lots of tasks written as prompts and a small model learns to handle new tasks it was never trained on.
- **Applicability (A1):** Confirms zero-shot task routing emerges from *explicit* prompt-format multitask training, not just scale — encouraging for a compute-modest biological trunk.
- **Novelty implication:** **Strengthens** feasibility at moderate scale; **pre-empts** "prompt-format multitask training" novelty.

### 14. Cross-Task Generalization via Natural Language Crowdsourcing Instructions / Super-NaturalInstructions
**Wang, Mishra et al. — EMNLP 2022.** arXiv:2204.07705
- **Takeaway:** A 1,616-task benchmark with declarative NL instructions; a model trained on it (Tk-Instruct) beats larger instruction models on unseen tasks.
- **Technical summary:** Tasks are specified by plain-language definitions plus optional k-shot examples across 76 categories; generalization scales with number of tasks and model size. Tk-Instruct outperforms InstructGPT by >9% while much smaller.
- **Plain-English:** A giant catalog of tasks each described in words lets a model learn to follow instructions for tasks it has never seen.
- **Applicability (A1, A3):** The declarative-instruction benchmark methodology is a template for building a *biological* task-instruction suite to train and evaluate MORPHEUS's auto-detection (A1) and to measure held-out biological task generalization (A3 evaluation).
- **Novelty implication:** **Reframes** — points to a concrete gap MORPHEUS can fill: there is no equivalent broad, declarative *biology-task* instruction benchmark; building one is defensible novelty.

### 15. Language Models are Few-Shot Learners (GPT-3)
**Brown et al., OpenAI — NeurIPS 2020.** arXiv:2005.14165
- **Takeaway:** Establishes in-context learning: a frozen LM infers and performs a task from a few prompt examples with no gradient updates.
- **Technical summary:** A 175B autoregressive LM performs many NLP tasks in zero/one/few-shot regimes purely by conditioning on a natural-language prompt; performance scales with model size. The task is never selected explicitly — it is inferred from the prompt.
- **Plain-English:** A big model does new tasks just from a worded prompt and a couple of examples, without any retraining.
- **Applicability (A1):** The origin point for "the prompt infers the task." Everything in MORPHEUS A1 descends from this; design implication: the trunk must be large/diverse enough for in-context task inference to *emerge*.
- **Novelty implication:** **Pre-empts** ownership of "in-context task inference" as a concept — MORPHEUS's contribution is domain grounding + identifiability, not the phenomenon.

### 16. An Explanation of In-Context Learning as Implicit Bayesian Inference
**Xie, Raghunathan, Liang, Ma — ICLR 2022.** arXiv:2111.02080
- **Takeaway:** Formalizes task auto-detection: the model infers a latent *concept/task* from the prompt via implicit Bayesian inference over pretraining structure.
- **Technical summary:** Using a mixture-of-HMMs generative model, they prove when latent-concept inference from a prompt succeeds despite train/prompt distribution mismatch, and validate on the synthetic GINC dataset where both Transformers and LSTMs show ICL.
- **Plain-English:** In-context learning works because the model quietly guesses which "topic" the prompt belongs to, then acts accordingly.
- **Applicability (A1, A2):** Gives a principled account of *how* a task is auto-detected — a latent-concept posterior. This is the theoretical scaffold for A2: if biological programmes are the latent concepts, identifiability determines whether the posterior can address them reliably.
- **Novelty implication:** **Strengthens** and **reframes** A2 — it justifies framing "reliable prompting" as "identifiable latent-concept inference," a sharper and more novel claim than generic prompting.

### 17. What Can Transformers Learn In-Context? A Case Study of Simple Function Classes
**Garg, Tsipras, Liang, Valiant — NeurIPS 2022.** arXiv:2208.01066
- **Takeaway:** Transformers trained from scratch can in-context learn whole function *classes* (linear, sparse, 2-layer nets, decision trees) at near-optimal accuracy.
- **Technical summary:** Given in-context (x, f(x)) pairs, the model infers an unseen f and predicts on new x, matching least-squares/optimal estimators and remaining robust to distribution shift. Establishes ICL as genuine algorithm learning, not memorization.
- **Plain-English:** Show a model example input-output pairs of an unknown rule and it figures out the rule and applies it — like statistics done inside the forward pass.
- **Applicability (A5, A1):** Directly relevant to A5: a (perturbation, response) example set could let MORPHEUS *infer a response function in-context* rather than retraining a classifier per perturbation.
- **Novelty implication:** **Strengthens** A5's "perturbation as a query, not a retrained model" by showing in-context function inference is real; MORPHEUS's delta is doing this over causal biological geometry, not synthetic functions.

### 18. Why Can GPT Learn In-Context? Language Models Implicitly Perform Gradient Descent as Meta-Optimizers
**Dai et al. — ACL 2023 Findings.** arXiv:2212.10559
- **Takeaway:** Interprets in-context learning as the forward pass implicitly running gradient-descent-like updates ("meta-optimization").
- **Technical summary:** Draws a dual form between attention over demonstrations and gradient updates, arguing ICL produces "meta-gradients" applied via attention; supported by similarity to explicit fine-tuning on several tasks.
- **Plain-English:** When a model learns from prompt examples, it's as if it secretly does a quick training step in its head.
- **Applicability (A1):** Mechanistic support that prompting can substitute for retraining — the premise behind MORPHEUS treating queries as prompts rather than new training runs.
- **Novelty implication:** **Strengthens** the "prompt not retrain" stance (A5/A1); a mechanistic citation for why MORPHEUS need not train a classifier per query.

### 19. In-Context Learning Creates Task Vectors
**Hendel, Geva, Globerson — EMNLP 2023 Findings.** arXiv:2310.15916
- **Takeaway:** ICL compresses a set of demonstrations into a single internal "task vector" that then modulates the model to execute the task.
- **Technical summary:** Shows the demonstration set maps to a compact vector at a specific layer; injecting that vector alone (without demos) reproduces the task behavior. Task and query are thus separable, with the task represented explicitly and causally.
- **Plain-English:** A prompt's instructions get squeezed into one internal "task setting" knob the model then turns.
- **Applicability (A2, A1):** Evidence that the inferred task lives in an *identifiable, manipulable* internal vector — a direct precedent for A2's addressable programme slots. Design implication: MORPHEUS could expose/steer a biological "task vector" as the addressing mechanism.
- **Novelty implication:** **Strengthens** A2 feasibility (task representations are extractable/addressable) but partly **pre-empts** "the task lives in an addressable latent"; MORPHEUS must show *biological programme* addressability with identifiability guarantees, beyond a single opaque vector.

### 20. Function Vectors in Large Language Models
**Todd, Li, Sen Sharma, Mueller, Wallace, Bau — ICLR 2024.** arXiv:2310.15213
- **Takeaway:** Specific attention heads carry compact, causal "function vectors" that transport a task abstraction and can be composed by addition.
- **Technical summary:** Via causal mediation analysis they isolate FVs that trigger zero-shot task execution across contexts and are partially *composable* (summing FVs yields novel composite tasks). Demonstrates internal task representations are localized and steerable.
- **Plain-English:** Inside the model, small bundles of neurons encode "do this operation," and you can add two together to get a combined operation.
- **Applicability (A2, A5):** Compositionality of function vectors is a strong precedent for A2 (per-programme addressable slots) and A5 (composing perturbations). Design implication: aim for programme representations that *add* to yield combinatorial perturbations, echoing GEARS-style compositional prediction but at the representation level.
- **Novelty implication:** **Strengthens and partly pre-empts** A2 — identifiable, composable task representations already exist in LLMs; MORPHEUS's novel burden is *guaranteed identifiability tied to named biological pathways*, not merely discovered-post-hoc vectors.

### 21. Learning Task Representations from In-Context Learning
**(2025).** arXiv:2502.05390
- **Takeaway:** Proposes learning explicit, robust task representations from ICL rather than reading them out post-hoc, improving controllability of task inference.
- **Technical summary:** Introduces a method to distill/optimize task representations from in-context demonstrations (e.g., attention-weighted aggregation across heads) that are more faithful and transferable than single extracted vectors, improving zero-shot transfer and steerability.
- **Plain-English:** Instead of hoping the model's internal "task knob" is good, deliberately train a clean version of it.
- **Applicability (A2):** Moves from *discovering* task vectors to *engineering* them — the trajectory MORPHEUS needs for reliable, addressable programme slots. Design implication: treat programme addressability as a trained objective, not an emergent accident.
- **Novelty implication:** **Reframes** A2 as an active design problem and **pre-empts** the generic version; MORPHEUS's identifiability-with-biological-semantics remains open.

### 22. Data Distributional Properties Drive Emergent In-Context Learning in Transformers
**Chan et al., DeepMind — NeurIPS 2022.** arXiv:2205.05055
- **Takeaway:** In-context (vs in-weight) task inference emerges only under specific data properties: burstiness, many rare classes, and dynamic meaning.
- **Technical summary:** Controlled experiments show ICL appears when training data is bursty and has a Zipfian/long-tailed class distribution; otherwise models default to in-weights memorization. Identifies distributional preconditions for the auto-detection capability.
- **Plain-English:** Whether a model learns to figure tasks out on the fly (rather than memorizing) depends on the *shape* of its training data.
- **Applicability (A1, A3):** Prescriptive for MORPHEUS's pretraining corpus: to get *biological* task auto-detection to emerge, the single-cell/omics corpus must be bursty and long-tailed (many rare cell states/perturbations). Design implication: curate for distributional properties, not just volume.
- **Novelty implication:** **Reframes** A1/A3 — success is contingent on data statistics; a MORPHEUS claim about emergent biological knowledge (A3) should be *conditioned on and tested against* these distributional preconditions, a rarely-made and defensible point.

### 23. Visual Prompting via Image Inpainting
**Bar, Gandelsman, Darrell, Globerson, Efros — NeurIPS 2022.** arXiv:2209.00647
- **Takeaway:** Task specification by a demonstration pair solved as inpainting — the earliest clean "visual in-context prompting" result, trained on arXiv figures.
- **Technical summary:** Concatenate an example (input, output) pair with a new input, mask the output region, and let a masked autoencoder inpaint it; trained on 88k unlabeled academic-figure images. Handles segmentation, detection, colorization, edge detection with no task-specific fine-tuning.
- **Plain-English:** Lay out "example in / example out / new in / blank" and the model fills the blank with the right answer.
- **Applicability (A1):** Minimalist proof that task auto-detection needs no labels or task IDs — only a well-posed prompt layout. Informs a data-only query mode for MORPHEUS.
- **Novelty implication:** **Strengthens** the demonstration-as-task idea; underscores that MORPHEUS's *NL grounding* (A3) is what separates it from pure visual-analogy prompting.

### 24. Mixture-of-Experts with Expert Choice Routing
**Zhou et al., Google — NeurIPS 2022.** arXiv:2202.09368
- **Takeaway:** Inverts routing — experts choose tokens (fixed bucket sizes) instead of tokens choosing top-k experts — fixing load imbalance and improving specialization.
- **Technical summary:** Each expert selects its top tokens, giving variable experts-per-token and balanced utilization; improves training efficiency and downstream quality over token-choice top-k routing.
- **Plain-English:** Rather than each input picking its specialist, each specialist picks the inputs it's best for, so no specialist is starved.
- **Applicability (A1):** Core reference for *prompt/content-conditioned routing* — the mechanism by which a detected task is dispatched to specialized sub-networks. Informs whether MORPHEUS routes biological queries to programme-specialized experts vs one dense trunk.
- **Novelty implication:** **Reframes** the routing-vs-dense design choice; MORPHEUS should justify its routing story (per-pathway experts?) against expert-choice as a known-good baseline rather than claim routing itself as novel.

### 25. Symbolic Mixture-of-Experts: Adaptive Skill-based Routing for Heterogeneous Reasoning
**(2025).** arXiv:2503.05641
- **Takeaway:** Routes queries to experts by *symbolic skill descriptions* (natural-language skill tags), enabling instance-level, interpretable routing across heterogeneous tasks.
- **Technical summary:** Rather than learned opaque gating, it selects pretrained expert models/adapters based on textual skill requirements inferred from the query, mixing them per instance; improves heterogeneous reasoning over fixed or task-label routing.
- **Plain-English:** The system reads what *kind of skill* a question needs, in words, and calls the matching specialists.
- **Applicability (A1, A2, A4):** The tightest analog to MORPHEUS's A1+A2: NL-inferred routing to *named, addressable* skill experts. For A4, symbolic routing could decide whether to invoke an encoded-modality expert or a RAG-context expert.
- **Novelty implication:** **Pre-empts** "NL-inferred routing to named skill modules" in the general case — MORPHEUS must ground the named modules in *identifiable biological pathways* (A2/A3) to stay novel; otherwise this is close prior art for the routing claim.

### 26. InstructSeq: Unifying Vision Tasks with Instruction-conditioned Multi-modal Sequence Generation
**Fang et al. — 2023.** arXiv:2311.18835
- **Takeaway:** Free-form NL instructions (LLM-generated) condition a single autoregressive model to perform diverse vision tasks without task-specific fine-tuning.
- **Technical summary:** A visual encoder + text encoder feed an autoregressive transformer that generates the task-appropriate output sequence; trained with diverse LLM-written instructions. Strong on semantic/referring segmentation, comprehension, and captioning under one model.
- **Plain-English:** Describe the vision task however you like in words and one model produces the answer.
- **Applicability (A1, A3):** Shows *open-vocabulary NL* (not fixed templates) can drive task auto-detection — the free-form-instruction target for MORPHEUS A1, and a recipe (LLM-generated instructions) for building A3 training data.
- **Novelty implication:** **Pre-empts** free-form NL task specification in unified models; strengthens the case that MORPHEUS's novelty is biological grounding and identifiability, not the NL front-end.

### 27. Towards Generalist Biomedical AI (Med-PaLM M)
**Tu et al., Google — 2023 (npj Digital Medicine 2024).** arXiv:2307.14334
- **Takeaway:** A single generalist model handles 14 biomedical tasks over text, images, and *genomics* with one set of weights, via instructions.
- **Technical summary:** Med-PaLM M is evaluated on MultiMedBench (14 tasks incl. genomic variant calling, report generation, med QA) using shared weights; reaches competitive/SOTA results, with clinicians sometimes preferring its chest-X-ray reports over radiologists'.
- **Plain-English:** One medical AI answers questions, reads scans, and interprets genomes, told what to do in plain language.
- **Applicability (A1, A3, A4):** The clearest existing "generalist biomedical model with NL task interface incorporating genomics" — the nearest neighbor to MORPHEUS in domain. Design implication and warning: MORPHEUS must differentiate on A2 (identifiable pathway slots) and A5 (causal/interventional queries), which Med-PaLM M does not attempt.
- **Novelty implication:** **Strongest pre-emption risk in-domain** — a reviewer will ask "how is this not Med-PaLM M for cells?" MORPHEUS's answer must be identifiability (A2), interventional causal queries (A5), and measured emergent-knowledge elicitation (A3), none of which Med-PaLM M claims.

### 28. Medical Vision Generalist: Unifying Medical Imaging Tasks in Context
**Ren et al., JHU — 2024.** arXiv:2406.05565
- **Takeaway:** An in-context, image-to-image generalist for medical imaging across CT/MRI/X-ray/ultrasound, task set by a demonstration.
- **Technical summary:** Combines masked image modeling + autoregressive training to do conditional image-to-image generation; a context example specifies the task. Beats prior vision generalists across 13 datasets/4 modalities and adapts to unseen datasets.
- **Plain-English:** Show a medical example of the transformation you want and it applies it to a new scan, across many imaging types.
- **Applicability (A1, A4):** Demonstrates Painter/SegGPT-style in-context generalism *inside biomedicine* — proof the paradigm transfers to clinical data; informs A4 handling of multiple imaging modalities in one trunk.
- **Novelty implication:** **Pre-empts** "in-context medical generalist"; reinforces that MORPHEUS's molecular/causal axes (A2/A5), not the generalist-in-medicine framing, must carry novelty.

### 29. scGPT: Toward Building a Foundation Model for Single-Cell Multi-Omics Using Generative AI
**Cui, Wang, Maan et al. — Nature Methods 2024.** DOI:10.1038/s41592-024-02201-0 (bioRxiv 2023.04.30.538439)
- **Takeaway:** A generative transformer pretrained on 33M cells that, via fine-tuning/prompting, serves cell-type annotation, batch integration, multi-omic integration, perturbation prediction, and gene-network inference from one backbone.
- **Technical summary:** Uses a gene/expression token scheme with attention over gene sets and a generative masked-prediction objective; downstream tasks are adapted with task-specific heads/fine-tuning. Reports competitive-to-SOTA performance across the five task families, and its attention maps recover gene programmes.
- **Plain-English:** A "GPT for cells" pretrained on tens of millions of cells that can be pointed at many single-cell analysis jobs.
- **Applicability (A1, A2, A4, A5):** The central in-domain reference. Its *per-task fine-tuning* is exactly the hard-coded-probe pattern MORPHEUS A1 wants to replace with prompt-based auto-detection; its gene-programme attention is a starting point for A2; its perturbation head is a non-promptable A5 baseline.
- **Novelty implication:** **Defines the primary contrast for MORPHEUS.** scGPT still uses task-specific heads/fine-tuning rather than a unified promptable interface — so MORPHEUS's A1 (NL auto-detection, no per-task head) and A5 (perturbation as a *query*, not a retrained head) are genuine deltas. Must be argued explicitly; this is the paper reviewers will compare against.

### 30. Geneformer: Transfer Learning Enables Predictions in Network Biology
**Theodoris et al. — Nature 2023.** DOI:10.1038/s41586-023-06139-9
- **Takeaway:** A transformer pretrained on ~30M single-cell transcriptomes enables transfer to many network-biology tasks and supports *in silico* perturbation to predict disease-relevant genes.
- **Technical summary:** Rank-value gene encoding + self-attention pretraining yields context-aware gene embeddings; fine-tuning transfers to tasks with limited data. *In silico* deletion/activation of genes in the embedding space predicts downstream network shifts and candidate therapeutic targets, validated experimentally.
- **Plain-English:** Pretraining on millions of cells lets the model do many biology tasks and simulate "what if we knock out this gene?" in the computer.
- **Applicability (A5, A2, A3):** Geneformer's *in silico perturbation* is a concrete precedent for A5's "perturbation as a query" — but it is a manual embedding-space edit, not an NL-promptable causal query. Its context-aware embeddings inform A2/A3.
- **Novelty implication:** **Partly pre-empts A5** (perturbation-in-latent-space exists) but **leaves MORPHEUS's core A5 open**: making the counterfactual a *promptable NL query with causal geometry* rather than a hand-specified embedding edit, plus measuring the elicited emergent knowledge (A3).

### 31. Predicting Transcriptional Outcomes of Novel Multigene Perturbations with GEARS
**Roohani, Huang, Leskovec — Nature Biotechnology 2024.** DOI:10.1038/s41587-023-01905-6 (bioRxiv 2022.07.12.499735)
- **Takeaway:** A graph-based model predicts post-perturbation expression for *unseen* single and *combinatorial* genetic perturbations by using a gene-relationship knowledge graph.
- **Technical summary:** Couples a GNN over a gene co-expression/GO knowledge graph with perturbation embeddings to generalize to perturbations absent from training, including 2-gene combinations, capturing genetic-interaction subtypes (synergy/suppression).
- **Plain-English:** It forecasts how cells will respond to gene knockouts it has never seen, even combinations, by leveraging known gene relationships.
- **Applicability (A5, A4):** The canonical "perturbation prediction" system — a strong A5 baseline and a template for combinatorial/counterfactual queries. Its knowledge-graph use informs A4's RAG-context question (relationships as external structured context vs encoded).
- **Novelty implication:** **Pre-empts** "predict unseen combinatorial perturbations" as a capability, but does so as a *dedicated retrained predictor*. MORPHEUS's A5 novelty is reframing this as one *prompt* to a unified promptable trunk (no bespoke model per query family), and unifying it with A1/A3 — a positioning, not a capability, claim, so it must be made carefully.

### 32. Exploring the Limits of Transfer Learning with a Unified Text-to-Text Transformer (T5)
**Raffel et al., Google — JMLR 2020.** arXiv:1910.10683
- **Takeaway:** Casts *every* NLP task as text-to-text so one model, format, and loss handle translation, QA, classification, summarization — the template for unified interfaces.
- **Technical summary:** Introduces the text-to-text framing with task-prefix strings ("translate English to German: ...") and studies the design space at scale (C4 corpus). Achieves SOTA across many benchmarks with a single architecture and objective.
- **Plain-English:** Turn every language task into "text in, text out" and one model with worded task prefixes does them all.
- **Applicability (A1):** The conceptual ancestor of promptable unified models: a *task-prefix string selects behavior*. MORPHEUS's NL task tag is a biological descendant of the T5 prefix.
- **Novelty implication:** **Pre-empts** the unified-interface-via-text-prefix idea at its root; confirms MORPHEUS should not claim the unified promptable *format* as novel, only its biological identifiability, grounding, and causal-query axes.

---

### Cross-cutting synthesis for MORPHEUS positioning

- **A1 (promptable unified + NL auto-detection) is mature prior art** across vision, language, and multimodal (Gato, OFA, Unified-IO 1/2, FLAN/T0, GPT-3, T5, InstructSeq). MORPHEUS cannot claim the *paradigm*; it can claim the first rigorous *biological* instantiation with measured task auto-detection.
- **Task auto-detection has a mechanistic account** (Xie implicit-Bayes; Garg function classes; ICL task vectors; function vectors) that MORPHEUS can borrow to make A2 precise: reliable prompting = identifiable latent-concept/programme inference. This is the strongest bridge from the ML literature into MORPHEUS's novel A2 claim.
- **Identifiable, composable task representations already exist** (Task Vectors, Function Vectors, Learning Task Representations). This both *supports* A2 feasibility and *raises the bar*: MORPHEUS must add biological-pathway semantics + identifiability guarantees, not just "a steerable vector."
- **In-domain pre-emption is concentrated in Med-PaLM M, scGPT, Geneformer, GEARS.** These collectively cover generalist biomedical NL interfaces, single-cell foundation trunks, and perturbation/in-silico queries — but *none* unify a promptable NL front-end with identifiable pathway slots and causal/counterfactual queries. That intersection is the defensible MORPHEUS core.
- **A4 (encode-vs-RAG decision; frozen-trunk plug-in):** Flamingo establishes frozen-trunk adapters; Unified-IO establishes encode-everything tokenization; GEARS uses a knowledge graph as external structure. No lane paper articulates a *decision rule* for when to encode a modality vs treat it as retrieval context — a genuine gap.
- **A5 (perturbation/drug as a query, not a retrained classifier):** Garg (in-context function inference), Geneformer (in-silico perturbation), and GEARS (combinatorial perturbation) each realize *pieces*, but always via bespoke heads or manual latent edits. Reframing the counterfactual as a *prompt to a unified trunk* is the least-pre-empted MORPHEUS axis.
