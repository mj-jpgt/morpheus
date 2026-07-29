## Compositional & instruction-following multimodal

Lane id: `l14_compositional_instruction`. Remit: compositional generalization, instruction-following multimodal models, in-context learning (ICL) for new tasks, compositional binding / VSA / hyperdimensional computing, rare-combination generalization. Every entry maps to one or more MORPHEUS rebase axes:
- **A1** promptable unified representation + NL task auto-detection
- **A2** identified, pathway-addressable slots (identifiability, per-programme addressability)
- **A3** NL<->biology grounding + emergent-knowledge elicitation AND its evaluation
- **A4** multimodal prompting: when to ENCODE a modality vs treat it as RAG context; frozen-trunk plug-in
- **A5** interventional/causal-geometry queries (perturbation/drug as a query, not a retrained classifier)

---

### Instruction-following & prompted task generalization (A1)

1. **Finetuned Language Models Are Zero-Shot Learners (FLAN)** — Wei, Bosma, Zhao, Guu, Yu, Lester, Du, Dai, Le (Google), ICLR 2022. arXiv:2109.01652.
   - *Takeaway:* Fine-tuning on a diverse collection of tasks phrased as **natural-language instructions** unlocks zero-shot generalization to held-out task types.
   - *Technical summary:* Instruction-tuning a 137B LM on ~60 NLP datasets grouped into task clusters, evaluated by holding out whole clusters. FLAN beats zero-shot GPT-3 (175B) on 20/25 datasets and few-shot GPT-3 on 6. Ablations show the instruction phrasing itself (not just multitask exposure) drives generalization, and the benefit emerges only above ~100B parameters.
   - *Plain-English:* If you train a model to follow written task descriptions across many tasks, it can then do brand-new tasks it was only told about in words.
   - *Applicability:* **A1.** Direct template for MORPHEUS's "infer the requested task, route it" claim: task auto-detection is a learned instruction-following behavior, not a set of hard-coded probe heads. Design implication — MORPHEUS should train on a broad *mixture of biological tasks framed as NL instructions* and hold out task families to prove genuine routing rather than memorized heads.
   - *Novelty implication:* Pre-empts any claim that "NL task specification for a foundation model" is itself novel — it is table stakes in NLP. MORPHEUS novelty must live in the *biological grounding* of the instruction space and in identifiable slots, not in the mere existence of NL prompting.

2. **Multitask Prompted Training Enables Zero-Shot Task Generalization (T0)** — Sanh, Webson, Raffel, et al. (BigScience), ICLR 2022. arXiv:2110.08207.
   - *Takeaway:* Explicit **multi-prompt** training (many human-authored prompt templates per task) yields robust zero-shot task generalization at 11B, far smaller than GPT-3.
   - *Technical summary:* T0 fine-tunes T5 on a large hand-curated set of prompted datasets (PromptSource), with multiple natural-language templates per dataset. It matches or beats GPT-3 zero-shot on many held-out tasks despite being 16x smaller, and prompt diversity improves robustness to prompt wording.
   - *Plain-English:* Showing a model the same task worded many different ways teaches it to respond to the *intent* rather than the exact phrasing.
   - *Applicability:* **A1.** For MORPHEUS, argues that reliable task auto-detection needs *many paraphrases* of each biological ask (e.g., "which pathway is active?", "score TGF-beta signaling", "is EMT engaged?") so routing is intent-based and prompt-robust.
   - *Novelty implication:* Reframes reliability of prompting as a training-data-design problem; MORPHEUS should not claim reliability from architecture alone.

3. **ZeroPrompt: Scaling Prompt-Based Pretraining to 1,000 Tasks Improves Zero-Shot Generalization** — Xu, Wang, et al., EMNLP Findings 2022. arXiv:2201.06910.
   - *Takeaway:* Scaling the *number* of prompted training tasks (to ~1000) improves zero-shot generalization more than scaling model size in the small-model regime.
   - *Technical summary:* Introduces a genome-of-tasks style curriculum plus a "task-scaling" study showing near-log-linear zero-shot gains with task count, and a proposal-based method to search prompts. Demonstrates that task *breadth* is a first-class scaling axis.
   - *Plain-English:* The more different jobs you teach a model to follow instructions on, the better it improvises on jobs it has never seen.
   - *Applicability:* **A1/A3.** Motivates a MORPHEUS "task census": enumerate hundreds of distinct biological queries as instruction tasks. The number of distinct addressable biological asks is itself a scaling lever for emergent routing.
   - *Novelty implication:* Strengthens the case that MORPHEUS's differentiator is *task breadth in biology*, an axis the NLP world has already validated but that is under-explored for omics foundation models.

---

### Multimodal instruction-following & unified representation (A1, A4)

4. **Flamingo: a Visual Language Model for Few-Shot Learning** — Alayrac, Donahue, Luc, et al. (DeepMind), NeurIPS 2022. arXiv:2204.14198.
   - *Takeaway:* Interleaving frozen vision and frozen language backbones with lightweight cross-attention ("Perceiver Resampler" + gated xattn) yields a single model doing many multimodal tasks by **in-context prompting**.
   - *Technical summary:* Frozen pretrained vision encoder and frozen LM are bridged by trainable cross-attention layers; the model ingests arbitrarily interleaved image/video-text sequences. A single Flamingo sets SOTA on numerous benchmarks via few-shot prompting, beating models fine-tuned on far more task data.
   - *Plain-English:* Keep two strong pretrained models frozen, add a small connector, and you get a system you can teach new visual tasks just by showing a few examples in the prompt.
   - *Applicability:* **A4 (frozen-trunk plug-in) + A1.** Canonical blueprint for MORPHEUS's "when to ENCODE a modality vs RAG": encode a new omics modality by training only a lightweight adapter into a *frozen* biological trunk, preserving the base representation. Supports "plug-in modality" over full retraining.
   - *Novelty implication:* Pre-empts any claim that frozen-trunk multimodal plug-in is novel per se. MORPHEUS must justify *which* biological modalities warrant an encoder adapter vs cheap RAG context, and quantify the tradeoff — that decision rule is the novel contribution, not the frozen-trunk mechanism.

5. **Visual Instruction Tuning (LLaVA)** — Liu, Li, Wu, Lee, NeurIPS 2023 (oral). arXiv:2304.08485.
   - *Takeaway:* First to extend NL instruction tuning to the image-text space using GPT-generated multimodal instruction data, connecting a frozen vision encoder to an LLM via a simple projection.
   - *Technical summary:* Uses language-only GPT-4 to synthesize 150K image-grounded instruction-following samples; trains a projection + LLM on a CLIP vision encoder's features. Achieves strong general-purpose visual chat and SOTA on ScienceQA when combined.
   - *Plain-English:* Auto-generate "look at this and answer" training examples, and a chat model learns to follow open-ended visual instructions.
   - *Applicability:* **A1/A4.** For MORPHEUS, argues that instruction-following over a new modality can be bootstrapped with *synthetically generated* NL-omics instruction pairs, lowering the annotation barrier for A3 grounding.
   - *Novelty implication:* Reframes data cost: MORPHEUS's grounding corpus can be partly LLM-synthesized, but that raises an evaluation burden (A3) to prove the elicited knowledge is real, not hallucinated.

6. **InstructBLIP: Towards General-purpose Vision-Language Models with Instruction Tuning** — Dai, Li, Li, Tiong, Zhao, Wang, Li, Fung, Hoi, NeurIPS 2023. arXiv:2305.06500.
   - *Takeaway:* An **instruction-aware Query Transformer** extracts features *conditioned on the instruction*, giving SOTA zero-shot on 13 held-out datasets.
   - *Technical summary:* Curates 26 datasets into instruction format across held-in/held-out splits; the Q-Former receives the instruction text so the visual features it pulls are query-dependent rather than fixed. Beats BLIP-2 and Flamingo zero-shot.
   - *Plain-English:* Let the question steer which parts of the image the model looks at, and it answers unseen question types better.
   - *Applicability:* **A1/A2.** Strong architectural argument for MORPHEUS A2: the *prompt should condition the readout* of the shared representation (instruction-aware pooling over pathway slots), not just post-hoc probe a fixed embedding. Query-conditioned addressing is the mechanism that makes prompting reliable.
   - *Novelty implication:* Pre-empts "prompt-conditioned readout" as novel; MORPHEUS's contribution is making that readout land on *identifiable biological programmes* (A2), a stronger guarantee than instruction-aware attention alone.

7. **ImageBind: One Embedding Space To Bind Them All** — Girdhar, El-Nouby, Liu, Singh, Alwala, Joulin, Misra (Meta), CVPR 2023 (highlight). arXiv:2305.05665.
   - *Takeaway:* Six modalities can be bound into one embedding space using **only image-paired data** — emergent cross-modal alignment without exhaustive pairwise data.
   - *Technical summary:* Each modality (image, text, audio, depth, thermal, IMU) is aligned to the image embedding via contrastive learning on image-paired data; modalities never trained together still align ("emergent" zero-shot cross-modal retrieval, arithmetic composition, and generation).
   - *Plain-English:* Anchor everything to images and all the other data types line up with each other for free, even pairs you never trained on.
   - *Applicability:* **A4/A1.** Key design lesson for MORPHEUS multimodal prompting: a single "anchor" modality (e.g., transcriptome) can bind proteomics/phospho/CNV/SNV without needing every pairwise dataset — supports encode-via-anchor for sparse omics pairings. Emergent cross-modal transfer is exactly the property MORPHEUS wants for rarely co-measured assays.
   - *Novelty implication:* Reframes A4: the question "encode vs RAG" is partly answered by data availability — anchor-binding lets you *encode* even sparsely-paired modalities. MORPHEUS should test whether omics modalities bind emergently through a shared anchor.

8. **MetaMorph: Multimodal Understanding and Generation via Instruction Tuning** — Tong, Fan, et al., 2024. arXiv:2412.14164.
   - *Takeaway:* A single instruction-tuned model that both *understands* and *generates* across modalities, showing understanding and generation reinforce each other.
   - *Technical summary:* Unifies visual understanding and generation in one autoregressive instruction-tuned model; finds that co-training the two objectives yields mutual benefit and that visual generation can emerge from predominantly understanding-oriented data.
   - *Plain-English:* Teaching one model to both read and draw makes it better at both than training them separately.
   - *Applicability:* **A1/A5.** Supports a MORPHEUS design where the same trunk both *reads* omics states and *generates* counterfactual states (A5 perturbation-as-query) — generation and understanding as one instruction interface.
   - *Novelty implication:* Strengthens the A5 claim that a generative/counterfactual head can share the trunk with the readout, rather than being a bolted-on classifier.

9. **SVIT: Scaling up Visual Instruction Tuning** — Zhao, Liu, Zhang, et al., 2023. arXiv:2307.04087.
   - *Takeaway:* Scaling visual-instruction *data* (to millions of GPT-generated QA/description/reasoning items) improves multimodal instruction following.
   - *Technical summary:* Constructs a 4.2M-sample multimodal instruction corpus (conversation, complex reasoning, detailed description, referring QA) and shows consistent gains, isolating data scale as a lever.
   - *Plain-English:* More and richer "describe/answer/reason about this image" examples make a multimodal assistant noticeably sharper.
   - *Applicability:* **A1/A3.** Quantitative precedent for how much NL-omics instruction data MORPHEUS may need; and a reminder that data scale, not just architecture, gates instruction-following quality.
   - *Novelty implication:* Neutral-to-pre-empting; underscores that MORPHEUS's grounding-corpus scale is a decisive, reportable variable.

10. **Otter: A Multi-Modal Model with In-Context Instruction Tuning** — Li, Zhang, Chen, Gao, et al., 2023. arXiv:2305.03726.
    - *Takeaway:* Adds explicit **in-context instruction tuning** on interleaved examples to an OpenFlamingo backbone, improving few-shot instruction following.
    - *Technical summary:* Builds MIMIC-IT (multimodal in-context instruction tuning data) with in-context example triplets; the resulting model follows instructions and demonstrations provided in-context better than the base VLM.
    - *Plain-English:* Train the model on prompts that include worked examples, and it gets better at learning new visual tasks from examples at inference.
    - *Applicability:* **A1/A4.** For MORPHEUS, supports a demonstration-based prompting mode ("here are 3 labeled samples, now score the 4th") layered on the NL interface — useful when a biological task has no dedicated head.
    - *Novelty implication:* Pre-empts "few-shot demonstration prompting for a bio model" as inherently novel; the novelty is demonstrations over *addressable biological slots*.

---

### Compositional generalization: theory & measurement (A2, A3)

11. **When Does Compositional Structure Yield Compositional Generalization? A Kernel Theory** — Lippl, Stachenfeld, ICLR 2025. arXiv:2405.16391.
    - *Takeaway:* Even with perfectly compositional representations, kernel/readout models are restricted to **conjunction-wise additive** functions and fail on tasks needing non-additive composition (e.g., transitive relations).
    - *Technical summary:* Analyzes fixed compositionally-structured representations under kernel regression; proves the learnable function class is limited to summing per-component values, and identifies failure modes — *memorization leakage* and *shortcut bias* from biased training data. Predictions match CNNs, ResNets, and ViTs.
    - *Plain-English:* Having the right building blocks in your representation is not enough — a simple readout can only *add up* the pieces, so it still botches tasks where the pieces interact.
    - *Applicability:* **A2/A3.** Critical caution for MORPHEUS: even if pathway slots are cleanly identified (A2), a *linear/additive probe* on them cannot answer non-additive biology (epistasis, synergy). Implies MORPHEUS needs an expressive, possibly interventional readout (A5), and that training data bias will leak into apparent "knowledge."
    - *Novelty implication:* **Reframes** the A2 claim — identifiable slots are necessary but not sufficient; MORPHEUS must show its readout can express non-additive combinations, else it over-claims compositional generalization.

12. **A Theoretical Analysis of Compositional Generalization in Neural Networks: A Necessary and Sufficient Condition** — Yuanpeng Li, 2025. arXiv:2505.02627.
    - *Takeaway:* States a necessary-and-sufficient condition: the computational graph must **align with the true compositional structure** AND components must encode exactly the right amount of information.
    - *Technical summary:* Proves that compositional generalization requires structural alignment between model and task plus information-tight component encodings, unifying architecture, regularization, and data properties. Offers an a-priori test of whether a model *can* compositionally generalize before training.
    - *Plain-English:* A model generalizes compositionally only if its internal wiring matches how the task actually decomposes and each part carries just the right information — no more, no less.
    - *Applicability:* **A2.** Gives MORPHEUS a principled target for slot identifiability: slots should encode *exactly* one biological programme's information (not entangled, not lossy). "Right amount of information" is a design spec for A2 addressability.
    - *Novelty implication:* Strengthens A2 by supplying the theoretical yardstick MORPHEUS can cite for why identifiable, information-tight slots enable reliable prompting.

13. **Measuring Compositional Generalization: A Comprehensive Method on Realistic Data (CFQ)** — Keysers, Schärli, Scales, et al. (Google), ICLR 2020. arXiv:1912.09713.
    - *Takeaway:* Introduces **distribution-based compositional splits** (maximize compound divergence, minimize atom divergence — "MCD") as a rigorous way to measure compositional generalization.
    - *Technical summary:* Builds the CFQ semantic-parsing benchmark and the DBCA method: train/test share atoms but differ maximally in atom *combinations*. Standard seq2seq models degrade sharply as compound divergence rises.
    - *Plain-English:* The fair test of "can it recombine known parts" is to make sure the test uses familiar parts in unfamiliar arrangements — and models struggle exactly there.
    - *Applicability:* **A3 (evaluation).** Directly usable for MORPHEUS's emergent-knowledge evaluation: construct omics splits where individual genes/pathways are seen but their *co-perturbation combinations* are held out (rare-combination generalization). MCD is a ready-made metric for A3/A5 claims.
    - *Novelty implication:* Provides the evaluation rigor MORPHEUS needs; without an MCD-style split, claims of "emergent combinatorial knowledge" are not credible.

14. **Compositionality Decomposed: How Do Neural Networks Generalise?** — Hupkes, Dankers, Mul, Bruni, JAIR 2020. arXiv:1908.08351.
    - *Takeaway:* Decomposes "compositionality" into five testable behaviors (systematicity, productivity, substitutivity, localism, overgeneralization) with concrete tests.
    - *Technical summary:* Proposes a taxonomy and PCFG-SET test suite; shows different architectures pass different aspects, arguing the field must specify *which* compositional property is claimed.
    - *Plain-English:* "Compositional" means several different things; you have to say which one you're testing and measure each separately.
    - *Applicability:* **A3 (evaluation).** MORPHEUS should specify which compositional property it claims for biology (e.g., systematicity = new gene x known context; productivity = longer perturbation combos) and test each, not use a blanket "generalizes" claim.
    - *Novelty implication:* **Reframes** how MORPHEUS should phrase compositional claims — decomposed and measured, reducing over-claim risk.

15. **Human-like Systematic Generalization Through a Meta-Learning Neural Network (MLC)** — Lake, Baroni, Nature 2023, 623:115-121. doi:10.1038/s41586-023-06668-3.
    - *Takeaway:* A standard seq2seq network **meta-trained on a stream of compositional tasks** reaches human-level systematic generalization, rebutting Fodor-Pylyshyn.
    - *Technical summary:* MLC optimizes a network over a dynamic curriculum of few-shot compositional episodes (learn new "words," use them in rule-like combinations). It matches or exceeds humans on systematic generalization benchmarks without built-in symbolic machinery.
    - *Plain-English:* Instead of hand-coding rules, keep giving a network fresh little compositional puzzles, and it learns to recombine new concepts the way people do.
    - *Applicability:* **A1/A3.** Suggests MORPHEUS could *meta-learn* biological task routing/composition by training over a stream of synthetic biological composition episodes (new pathway defined in-context, then queried), improving rare-combination generalization without symbolic scaffolding.
    - *Novelty implication:* Strengthens the feasibility of emergent compositional biology from ordinary architectures + the right training regime — supports MORPHEUS's "emergent, not hard-coded" thesis while warning the *training curriculum* is the real innovation.

16. **From Frege to ChatGPT: Compositionality in Language, Cognition, and Deep Neural Networks** — Baroni, Pavlick, et al., 2024 (review). arXiv:2405.15164.
    - *Takeaway:* Synthesizes decades of compositionality debate and argues modern LLMs exhibit partial, graded compositionality rather than the all-or-nothing symbolic kind.
    - *Technical summary:* Reviews formal, cognitive, and empirical notions of compositionality; catalogs where neural nets succeed/fail and how evaluation choices shape conclusions.
    - *Plain-English:* Compositionality in neural nets is real but partial and messy — how you test decides what you conclude.
    - *Applicability:* **A3.** Framing reference for MORPHEUS's positioning of "emergent biological knowledge" as graded and evaluation-dependent.
    - *Novelty implication:* Reframes claims toward graded/measured language; discourages binary "MORPHEUS understands biology" statements.

17. **Consistency of Compositional Generalization Across Multiple Levels** — Sun, Li, et al., AAAI 2025. arXiv:2412.13636.
    - *Takeaway:* Models can generalize at one compositional level while failing at others; proposes training for *multi-level* consistency.
    - *Technical summary:* Defines a hierarchy of compositional levels and a meta-learning method enforcing consistency across them on VQA/semantic-parsing; shows single-level success overstates true compositional ability.
    - *Plain-English:* A model may recombine small pieces fine but fall apart on bigger combinations — you have to check every level.
    - *Applicability:* **A3.** MORPHEUS should evaluate rare-combination generalization at multiple scales (pairwise, triple, higher-order perturbations) rather than reporting one level.
    - *Novelty implication:* Pre-empts single-level "combinatorial generalization" claims; raises the evaluation bar for A5.

18. **Systematic Generalization in Language Models Scales with Information Entropy** — 2025. arXiv:2505.13089.
    - *Takeaway:* Systematic-generalization performance correlates with the **information entropy of the training distribution's compositional structure**.
    - *Technical summary:* Shows that as the entropy (diversity) of component combinations in training rises, systematic generalization improves predictably, giving a quantitative handle on data design.
    - *Plain-English:* The more varied the combinations you train on, the better the model recombines — and you can measure this with entropy.
    - *Applicability:* **A3/A5.** Prescriptive for MORPHEUS training-set design: maximize entropy over perturbation/context combinations to buy combinatorial generalization; entropy becomes a reportable design metric.
    - *Novelty implication:* Strengthens a data-centric route to MORPHEUS's combinatorial claims and gives a metric to defend them.

19. **Swing-by Dynamics in Concept Learning and Compositional Generalization** — Park, Lubana, Tanaka, et al., 2024. arXiv:2410.08309.
    - *Takeaway:* Compositional abilities emerge via non-monotonic ("swing-by") learning dynamics — abrupt, structured phase transitions during training.
    - *Technical summary:* Analyzes diffusion/generative concept learning and finds a two-stage trajectory where compositional generalization appears suddenly after the model organizes concept factors; connects to feature-learning geometry.
    - *Plain-English:* The ability to combine concepts doesn't creep in smoothly — it snaps into place partway through training.
    - *Applicability:* **A3.** Warns MORPHEUS that emergent compositional/biological knowledge may appear only after a training phase transition — early checkpoints will look like failure; evaluate across training, not just at the end.
    - *Novelty implication:* Adds nuance to "emergent knowledge" claims — emergence is dynamical, so MORPHEUS should document the trajectory, not a single snapshot.

20. **Representational Homomorphism Predicts and Improves Compositional Generalization in Transformer Language Models** — An, Du, 2026. arXiv:2601.18858.
    - *Takeaway:* A **Homomorphism Error (HE)** metric — misalignment between syntactic composition rules and the model's hidden-state combination rules — predicts OOD compositional accuracy (R^2=0.73) and, when minimized, improves it.
    - *Technical summary:* On adapted SCAN with decoder-only transformers, HE measures whether "combine in representation space" mirrors "combine in symbol space." Regularizing to reduce HE gives a statistically significant OOD gain (p=0.023).
    - *Plain-English:* If the way a model blends internal representations mirrors the grammar's way of blending words, it generalizes better — and you can measure and enforce that.
    - *Applicability:* **A2/A3.** Offers MORPHEUS a concrete, trainable objective for A2: enforce that combining slot activations in representation space is *homomorphic* to combining biological programmes — a measurable identifiability/addressability criterion.
    - *Novelty implication:* Strengthens A2 with an operational metric; MORPHEUS could adapt HE to "biological homomorphism" as an evaluation of addressable slots.

---

### In-context learning as task inference (A1, A5)

21. **In-Context Learning Creates Task Vectors** — Hendel, Geva, Globerson, EMNLP Findings 2023. arXiv:2310.15916.
    - *Takeaway:* ICL compresses the demonstration set into a **single task vector** that modulates the transformer — ICL ≈ (infer task vector) then (apply to query).
    - *Technical summary:* Shows a mid-layer activation captures the demonstrated task; patching just this vector into a fresh forward pass reproduces ICL behavior across many tasks/models, evidencing an implicit "learn a rule, then apply it" decomposition.
    - *Plain-English:* When you show a model examples, it quietly boils them down to one internal "task setting" and then just runs that setting on your query.
    - *Applicability:* **A1/A2.** Direct mechanistic support for MORPHEUS's "infer the requested task and route it": task routing can be a *vector in activation space*, and an explicit, addressable "task slot" (A2) is biologically plausible. Suggests MORPHEUS could expose/steer a task vector.
    - *Novelty implication:* Pre-empts "task inference is magic"; it's a locatable vector. MORPHEUS's novelty is tying that vector to *named biological programmes*, making it inspectable and steerable.

22. **Function Vectors in Large Language Models** — Todd, Li, Sharma, Mueller, Wallace, Bau, ICLR 2024. arXiv:2310.15213.
    - *Takeaway:* A few attention heads transport a compact **function vector** encoding the demonstrated input-output mapping; FVs trigger the task zero-shot and **compose** via vector arithmetic.
    - *Technical summary:* Causal mediation analysis localizes FVs to specific mid-layer heads; adding an FV to unrelated contexts executes the task, and FVs for sub-behaviors can be combined to produce new complex behaviors.
    - *Plain-English:* The "what to do" of a task lives in a small handful of internal directions you can copy, paste, and even add together to build new tasks.
    - *Applicability:* **A2/A5.** Powerful evidence for MORPHEUS A2 (addressable, composable task/programme vectors) and A5 (compose perturbation queries by *vector arithmetic* rather than retraining a classifier). Implies interventional queries can be posed as algebra over identified directions.
    - *Novelty implication:* **Strengthens** A2 and A5 substantially — but also *pre-empts* the mechanism; MORPHEUS must show biological function vectors are identifiable, composable, and causally faithful in omics space, not just claim the analogy.

23. **Emergent World Representations: Exploring a Sequence Model Trained on a Synthetic Task (Othello-GPT)** — Li, Hopkins, Bau, Viégas, Pfister, Wattenberg, ICLR 2023 (oral). arXiv:2210.13382.
    - *Takeaway:* A GPT trained only to predict legal moves develops an **emergent, causally-manipulable internal representation** of board state; interventions on it steer predictions.
    - *Technical summary:* Probes recover a nonlinear board-state representation; latent interventions causally change legal-move outputs, and "latent saliency maps" interpret decisions — evidence of a genuine world model, not surface statistics.
    - *Plain-English:* A model taught only to predict moves secretly builds a picture of the board — and if you edit that picture, its predictions change accordingly.
    - *Applicability:* **A3/A5.** Template for MORPHEUS's emergent-knowledge *evaluation*: probe for an emergent biological state representation (e.g., pathway-activity manifold) and prove it is causal by intervening on it and observing consistent downstream shifts. This is the gold standard for A5's "causal geometry."
    - *Novelty implication:* **Strengthens** A3/A5 methodology — MORPHEUS can adopt probe-then-intervene as its emergent-knowledge test. The intervention step is what separates a real claim from a correlational probe.

24. **Emergent Abilities of Large Language Models** — Wei, Tay, Bommasani, et al., TMLR 2022. arXiv:2206.07682.
    - *Takeaway:* Some capabilities appear abruptly (non-linearly) only above a scale threshold — "emergent abilities."
    - *Technical summary:* Documents tasks where performance is near-random until a critical model scale, then jumps; frames emergence as unpredictable from small-scale behavior. (Later work contests whether some emergence is a metric artifact — see #25.)
    - *Plain-English:* Certain skills just don't exist in small models and then suddenly show up once the model is big enough.
    - *Applicability:* **A3.** Motivates MORPHEUS's emergent-biological-knowledge framing — but demands scale-sweep evidence to claim emergence rather than assert it.
    - *Novelty implication:* Supports the emergence narrative while flagging the burden of proof.

25. **Are Emergent Abilities in Large Language Models Just In-Context Learning?** — Lu, Bansal, et al., ACL 2024. arXiv:2309.01809.
    - *Takeaway:* Many "emergent" abilities are explained by in-context learning + instruction tuning + memory, not a genuinely new capability.
    - *Technical summary:* Controlled experiments separate ICL from claimed emergence; finds most purported emergent skills reduce to competent ICL and prompt-format effects, cautioning against over-attribution.
    - *Plain-English:* A lot of the "surprising new skills" are really the model cleverly using examples and instructions, not brand-new understanding.
    - *Applicability:* **A3 (evaluation caution).** Direct warning for MORPHEUS: distinguish *emergent biological knowledge* from clever prompt-following/retrieval. A3 evaluation must control for ICL and RAG leakage.
    - *Novelty implication:* **Pre-empts/threatens** naive emergence claims — MORPHEUS must design ablations showing knowledge is in the weights, not supplied by the prompt/context (bears directly on A4's encode-vs-RAG boundary).

---

### Compositional binding, VSA & hyperdimensional computing (A2)

26. **A Survey on Hyperdimensional Computing aka Vector Symbolic Architectures, Part I: Models and Data Transformations** — Kleyko, Rachkovskij, Osipov, Rahimi, ACM Computing Surveys 2023. doi:10.1145/3538531 (arXiv:2111.06077).
    - *Takeaway:* Comprehensive taxonomy of HDC/VSA models where **binding** (outer-product-like) and **bundling** (superposition) operations solve the variable-binding problem with high-dim distributed vectors.
    - *Technical summary:* Reviews TPR, HRR, MAP, Binary Spatter Codes, Sparse Binary Distributed Representations; shows how role-filler binding, permutation, and unbinding support structured symbolic computation in vector space with graceful degradation.
    - *Plain-English:* A math toolkit that lets you attach "this value goes in that slot" inside a single big vector, then pull items back out — structured symbols living in distributed vectors.
    - *Applicability:* **A2.** The formal backbone for MORPHEUS A2: identified, *addressable slots* = role-filler binding. VSA gives operations to bind "programme = TGF-beta" to "activity = high," store many bindings in superposition, and *unbind* to query — a principled mechanism for per-programme addressability and compositional prompts.
    - *Novelty implication:* **Reframes/strengthens** A2 by offering a rigorous vocabulary. If MORPHEUS's slots behave like VSA role-filler bindings, addressability claims gain formal grounding; if they don't, VSA exposes what's missing.

27. **A Survey on HDC/VSA, Part II: Applications, Cognitive Models, and Challenges** — Kleyko, Rachkovskij, Osipov, Rahimi, ACM Computing Surveys 2023. doi:10.1145/3558000 (arXiv:2112.15424).
    - *Takeaway:* Catalogs VSA applications (analogical reasoning, classification, sequence/graph encoding) and open challenges (capacity, learning bindings from data).
    - *Technical summary:* Surveys where role-filler binding + superposition have delivered, plus the capacity limits and the difficulty of *learning* rather than hand-designing bindings — the crux for integrating VSA with deep nets.
    - *Plain-English:* Part II shows what these bind/bundle tricks are good for and where they break — especially learning the bindings automatically.
    - *Applicability:* **A2.** Flags the key risk for MORPHEUS: learned slot-bindings must respect capacity limits; too many superposed programmes causes crosstalk. Informs how many pathway slots can be reliably addressed at once.
    - *Novelty implication:* Warns that A2 addressability has a *capacity ceiling* — MORPHEUS should quantify how many programmes it can bind before interference degrades prompting.

28. **Capacity Analysis of Vector Symbolic Architectures** — Clarkson, Ubaru, Yang (IBM), 2023. arXiv:2301.10352.
    - *Takeaway:* Rigorous bounds on the dimension needed for reliable membership/set operations, connecting VSAs to sketching and Bloom filters.
    - *Technical summary:* Analyzes MAP-I, MAP-B, and two sparse-binary VSAs; derives dimensional requirements for membership testing and set-intersection estimation, and studies a Hopfield variant — establishing VSA-sketching-Bloom-filter links.
    - *Plain-English:* Tells you exactly how big your vectors must be to reliably store and query a set of bound items.
    - *Applicability:* **A2.** Gives MORPHEUS quantitative guidance: representation width vs number of simultaneously addressable biological programmes. If slots are VSA-like, this bounds reliable prompting capacity.
    - *Novelty implication:* Strengthens A2 with hard numbers; MORPHEUS can cite capacity theory to justify how many slots it exposes.

29. **Vector Symbolic Architectures as a Computing Framework for Emerging Hardware** — Kleyko, Davies, Frady, Kanerva, Sommer, Rahimi, et al., Proceedings of the IEEE 2022. arXiv:2106.05268.
    - *Takeaway:* Positions VSA/HDC as a general programming model for structured, compositional computation, including neuromorphic implementations.
    - *Technical summary:* Reviews the algebra (bind/bundle/permute/cleanup), error resilience, and how VSA enables computing-in-superposition and factorization of composite representations ("resonator networks") — recovering constituents of a bound vector.
    - *Plain-English:* A blueprint for building compositional, symbol-like computation directly in vectors, robust to noise.
    - *Applicability:* **A2/A5.** "Factorization / resonator" is directly relevant: given a composite omics embedding, *factor it back* into constituent programme activations — an unbinding operation MORPHEUS could use to make slots addressable and to decompose combinatorial perturbations (A5).
    - *Novelty implication:* Offers a concrete mechanism (resonator factorization) MORPHEUS could adopt or must differentiate from — decomposing a state into its programmes is a VSA-solved problem, so MORPHEUS's contribution is doing it in *learned biological* space.

---

### Rare-combination / compositional zero-shot generalization (A2, A5)

30. **Learning Conditional Attributes for Compositional Zero-Shot Learning** — Wang, Chen, et al., CVPR 2023. arXiv:2305.17940.
    - *Takeaway:* Attributes are **context-conditioned** ("wet" in "wet apple" ≠ "wet cat"); learning conditional attribute embeddings improves recognition of unseen attribute-object pairs.
    - *Technical summary:* Proposes an attribute hyper-learner + base-learner producing object-conditioned attribute representations, improving CZSL on unseen compositions by modeling attribute-object interaction rather than treating them independently.
    - *Plain-English:* The meaning of a property depends on what it's attached to; modeling that dependence helps recognize new property-object combos.
    - *Applicability:* **A2/A5.** Direct analogue to biology: a pathway's "activation" means something different in different cell contexts. MORPHEUS slots should be *context-conditioned* (programme x cell-state), enabling rare-combination (perturbation x context) generalization — a core A5 requirement.
    - *Novelty implication:* Reframes A2: slots are not context-free labels but conditional embeddings; strengthens MORPHEUS's rare-combination claim if implemented, weakens it if slots are treated as context-independent.

31. **Graph-guided Cross-composition Feature Disentanglement for Compositional Zero-Shot Learning** — 2024. arXiv:2408.09786.
    - *Takeaway:* Uses a graph over compositions to **disentangle** attribute and object features across compositions, improving unseen-pair generalization.
    - *Technical summary:* Builds a composition graph and cross-composition disentanglement objective so attribute/object factors are shared consistently across pairs, reducing entanglement that hurts rare compositions.
    - *Plain-English:* By linking related combinations in a graph, the model cleanly separates "the property" from "the thing," so it handles new combinations better.
    - *Applicability:* **A2.** Suggests MORPHEUS could use a *pathway/interaction graph* to regularize slot disentanglement across contexts — graph-guided identifiability.
    - *Novelty implication:* Supports graph-regularized A2; a candidate mechanism MORPHEUS should consider or distinguish from prior-knowledge-graph approaches.

32. **Hybrid Discriminative Attribute-Object Embedding Network for Compositional Zero-Shot Learning** — 2024. arXiv:2412.00121.
    - *Takeaway:* Combines discriminative and compositional embeddings to better separate fine-grained attribute-object pairs for unseen combinations.
    - *Technical summary:* A hybrid network fuses holistic composition embeddings with disentangled attribute/object branches, improving fine-grained discrimination of rare compositions.
    - *Plain-English:* Blend a "whole combo" view with a "separate parts" view to tell apart subtle new combinations.
    - *Applicability:* **A2/A5.** Argues MORPHEUS may need *both* a holistic state embedding and disentangled programme slots — relevant to whether a modality is encoded holistically vs decomposed.
    - *Novelty implication:* Neutral; informs the holistic-vs-disentangled design tension.

33. **ConceptMix: A Compositional Image Generation Benchmark with Controllable Difficulty** — Wu, Zhu, Xie, et al., NeurIPS 2024 (Datasets & Benchmarks). arXiv:2408.14339.
    - *Takeaway:* A benchmark that **scales compositional difficulty** by the number of bound attributes (k objects x attributes) and auto-grades with VQA, exposing how generation degrades with compositional load.
    - *Technical summary:* Procedurally composes prompts with controllable k concepts, grades each concept's presence via a VQA model; SOTA text-to-image models drop sharply as k grows, revealing binding failures.
    - *Plain-English:* Ask a generator to combine more and more properties at once and grade each one — quality collapses as you pile on requirements.
    - *Applicability:* **A3/A5 (evaluation).** Template for a MORPHEUS combinatorial-load benchmark: prompt for k simultaneous perturbations/programmes and auto-score each — measuring where binding breaks under compositional load.
    - *Novelty implication:* Provides an evaluation paradigm MORPHEUS should adopt to substantiate combinatorial-generalization claims and to find its binding ceiling.

34. **Multi-Sourced Compositional Generalization in Visual Question Answering** — 2025. arXiv:2505.23045.
    - *Takeaway:* Studies compositional generalization when primitives come from **multiple modalities/sources**, showing cross-source recombination is harder than single-source.
    - *Technical summary:* Defines multi-sourced compositional splits in VQA (primitives from vision vs language) and shows models fail to recombine primitives that originate in different modalities, proposing alignment methods to close the gap.
    - *Plain-English:* Combining an idea learned from pictures with one learned from words is harder than combining two ideas from the same place.
    - *Applicability:* **A4/A5.** Highly relevant to MORPHEUS multimodal prompting: recombining a programme learned from transcriptomics with one from proteomics/phospho may fail without explicit cross-modal alignment — a concrete risk for A4 "encode multiple modalities" and A5 cross-modal counterfactuals.
    - *Novelty implication:* **Reframes** A4/A5 — cross-modal compositional generalization is a known hard case; MORPHEUS must explicitly evaluate cross-source recombination, not assume it transfers.

35. **Identifiability Guarantees for Causal Disentanglement from Soft Interventions** — Zhang, Squires, Greenewald, Srivastava, Shanmugam, Uhler, NeurIPS 2023. arXiv:2307.06250.
    - *Takeaway:* Proves latent causal variables are **identifiable** (up to an equivalence class) from unpaired observational + soft-interventional data, and predicts effects of **unseen intervention combinations** — demonstrated on combinatorial genomic perturbations.
    - *Technical summary:* Under a generalized faithfulness condition, an autoencoding VAE recovers the latent causal graph and enables prediction of novel combinations of interventions in the infinite-data limit; validated on predicting combinatorial perturbation effects in genomics.
    - *Plain-English:* From a mix of "no intervention" and "gentle nudge" data you can pin down the hidden causal factors and forecast what unseen combinations of nudges will do — shown for combining genetic perturbations.
    - *Applicability:* **A2 + A5 (both, tightly).** This is the closest prior art to MORPHEUS's core thesis: *identifiable* latent factors (A2) that support *counterfactual combinatorial intervention prediction* (A5) in biology, without retraining a classifier per query. The soft-intervention formalism is exactly "perturbation-as-query."
    - *Novelty implication:* **Strongest pre-emption risk in this lane.** MORPHEUS's A2+A5 claim (identifiable, pathway-addressable slots enabling counterfactual perturbation queries) overlaps heavily with this work in the genomics setting. MORPHEUS must differentiate on: (a) NL-promptable interface (A1) over the identified factors, (b) multimodal scope beyond single-omics perturbations (A4), and (c) scale/foundation-model framing rather than a task-specific causal VAE. Failure to distinguish would let a reviewer say "identifiable causal factors for combinatorial perturbation is already solved."
