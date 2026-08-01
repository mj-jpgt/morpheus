## Emergent-capability & elicitation evaluation

Lane l06_emergence_eval. Remit: measuring emergent abilities and latent knowledge in foundation models — the emergence debate, capability elicitation, probing for latent knowledge, eliciting-latent-knowledge (ELK), knowledge-neuron/concept probing, and evaluation methodology for "does the model know X". Primary axis served is **A3** (NL⇄biology grounding + emergent-knowledge elicitation *and its evaluation*); several entries also inform A1 (task auto-detection), A2 (identified/addressable slots), A4 (encode vs. retrieve), and A5 (interventional queries). All entries are real, findable papers with verified titles/IDs.

---

### Part I — The emergence debate: is "emergence" real, and how is it (mis)measured?

**1. Emergent Abilities of Large Language Models** (Wei, Tay, Bommasani, Raffel, Zoph, et al.; TMLR 2022) — arXiv:2206.07682
- *Takeaway:* Defines "emergent" as a capability absent in small models that appears sharply above a scale threshold and cannot be extrapolated from small-model scaling curves.
- *Technical summary:* Surveys dozens of tasks (multi-digit arithmetic, word unscrambling, multitask QA, chain-of-thought, instruction following) where performance stays at chance until a critical parameter/FLOP count, then rises sharply. The claim is descriptive: the discontinuity is presented as a genuine property of scale, not of measurement.
- *Plain-English:* Big models can suddenly do things small models can't, and you can't see it coming by watching the small ones.
- *Applicability (A3):* This is the phenomenon MORPHEUS's A3 wants to *measure in biology* — does a multimodal tumor model "suddenly know" pathway biology at scale? It motivates the need for an emergence/elicitation benchmark rather than a single downstream accuracy number.
- *Novelty implication:* **Reframes.** MORPHEUS should not claim "our model has emergent biological knowledge" as a headline; entries 3–6 show the term is contested. MORPHEUS's contribution must be the *measurement method*, not a claim to emergence per se.

**2. Predictability and Surprise in Large Generative Models** (Ganguli, Hernandez, Lovitt, et al., Anthropic; ACM FAccT 2022) — arXiv:2202.07785
- *Takeaway:* Aggregate loss is smoothly predictable while specific capabilities are not — a structural tension for any capability claim.
- *Technical summary:* Documents that broad training loss follows clean scaling laws, yet the emergence of specific downstream skills is unpredictable in timing and form. Argues this predictable-loss/unpredictable-capability combination is what makes evaluation and governance hard.
- *Plain-English:* You can forecast how "good" a model gets on average, but not exactly which new tricks it will pick up or when.
- *Applicability (A3):* Warns that MORPHEUS cannot validate a biological-knowledge claim from trunk loss or C-index; capability-specific probes are required. Directly supports the thesis's point that the current harness is "structurally blind to representation quality."
- *Novelty implication:* **Strengthens** the need for a dedicated elicitation harness; **pre-empts** any attempt to infer "the model knows biology" from headline metrics.

**3. Are Emergent Abilities of Large Language Models a Mirage?** (Schaeffer, Miranda, Koyejo; NeurIPS 2023) — arXiv:2304.15004
- *Takeaway:* Apparent emergence is largely an artifact of discontinuous/nonlinear metrics; under smooth metrics 25 of 29 tested metrics show continuous improvement.
- *Technical summary:* Shows that swapping exact-match/accuracy for a continuous per-token metric (e.g., token edit distance, Brier) turns "sharp emergence" into smooth, predictable curves. Attributes the illusion to metric choice plus too-few test items to estimate small-model performance.
- *Plain-English:* The "sudden ability" often disappears when you grade the model more gently — it was improving smoothly all along.
- *Applicability (A3, A2):* The single most important methodological warning for MORPHEUS's emergence benchmark: any biological-knowledge elicitation metric must be reported under smooth *and* thresholded scoring, with adequate test-set size, or the result is uninterpretable.
- *Novelty implication:* **Pre-empts a fatal critique.** If MORPHEUS reports a biology-knowledge "jump," a reviewer will invoke Schaeffer. MORPHEUS must build metric-robustness in from the start; this converts a vulnerability into a differentiator.

**4. Are Emergent Abilities in Large Language Models just In-Context Learning?** (Lu, Bigoulaeva, Sachdeva, Tayyar Madabushi, Gurevych; ACL 2024) — arXiv:2309.01809
- *Takeaway:* Across 1000+ experiments, purported emergent abilities decompose into in-context learning + memory + linguistic knowledge, not novel reasoning.
- *Technical summary:* Controls for in-context learning and instruction tuning and finds that most "emergent" task gains vanish, i.e., the ability was latent and elicited by the prompt format rather than newly created at scale.
- *Plain-English:* A lot of "emergence" is really the model being prompted the right way to use knowledge it already had.
- *Applicability (A3, A1):* Reframes A3 as *elicitation* rather than *creation* — the interesting MORPHEUS question is "can we elicit latent biological knowledge the trunk already encodes?", which is exactly A1's task-inferring prompt interface plus A3's readout.
- *Novelty implication:* **Reframes and strengthens.** Supports MORPHEUS's framing that the novel work is elicitation/measurement; weakens any "our model spontaneously invented biology" claim.

**5. Understanding Emergent Abilities of Language Models from the Loss Perspective** (Du, Zeng, Dong, Tang; NeurIPS 2024) — arXiv:2403.15796
- *Takeaway:* Models with equal pretraining loss show equal downstream performance regardless of size; emergence is a function of loss crossing a threshold, not parameters.
- *Technical summary:* Empirically ties downstream task accuracy to pretraining loss (not model/data scale independently) and redefines emergence as capabilities that appear once loss falls below a critical value, after which accuracy jumps from chance.
- *Plain-English:* What matters isn't how big the model is but how low its training loss got — hit a certain loss and the ability switches on.
- *Applicability (A3):* Gives MORPHEUS a principled x-axis for any emergence-of-biology plot: report capability against trunk loss, not parameter count, so results are comparable across ablations.
- *Novelty implication:* **Strengthens methodology.** Adopting loss-indexed emergence curves would make MORPHEUS's biological-knowledge claims more rigorous than parameter-indexed folklore.

**6. The Quantization Model of Neural Scaling** (Michaud, Liu, Girit, Tegmark; NeurIPS 2023) — arXiv:2303.13506
- *Takeaway:* Knowledge/skills come in discrete "quanta" learned in frequency order; power-law loss and apparent emergence both fall out of this decomposition.
- *Technical summary:* Proposes that a model learns discrete computational units ("quanta") in order of usage frequency; the sum of quantized improvements produces smooth power-law loss while individual quanta produce sharp per-skill jumps. Validated on toy and language models.
- *Plain-English:* Models learn skills as discrete building blocks, picking up the common ones first; each block clicks in suddenly even though the average curve looks smooth.
- *Applicability (A2, A3):* A theoretical bridge to MORPHEUS's A2: if biological "programmes/pathways" are the quanta, then *per-programme* addressable slots are the natural unit for measuring which pieces of biology the model has acquired.
- *Novelty implication:* **Strengthens A2.** Motivates measuring biology knowledge at the per-pathway (quantum) grain rather than as a monolithic score.

**7. Beyond the Imitation Game (BIG-bench): Quantifying and Extrapolating the Capabilities of Language Models** (Srivastava et al., 442 authors; TMLR 2022) — arXiv:2206.04615
- *Takeaway:* A 204-task benchmark showing scale improves capability *and* calibration but leaves large absolute gaps vs. experts, with several tasks showing "breakthrough" scaling.
- *Technical summary:* Introduces a massive, diverse, community-built benchmark explicitly to probe near-future capabilities; analyzes which tasks scale smoothly vs. show breakthroughs and studies calibration alongside accuracy.
- *Plain-English:* A giant grab-bag of hard tasks used to chart what language models can and can't yet do.
- *Applicability (A3):* Template for MORPHEUS's own emergence benchmark: a *diverse battery* of biological-knowledge probes (not one task), scored for both accuracy and calibration, with breakthrough-vs-smooth analysis.
- *Novelty implication:* **Reframes.** A "BIG-bench for tumor-state biology knowledge" is an open, defensible MORPHEUS deliverable — no such curated elicitation battery exists for multimodal cancer FMs.

**8. Inverse Scaling: When Bigger Isn't Better** (McKenzie et al.; TMLR 2023) — arXiv:2306.09479
- *Takeaway:* Documents tasks where larger models get *worse*, with four identified causes (e.g., preferring memorized sequences over instructions).
- *Technical summary:* Aggregates winning datasets from the Inverse Scaling Prize and categorizes failure modes where capability *decreases* with scale, showing scaling is not monotone and headline "more scale = more knowledge" is unsafe.
- *Plain-English:* Sometimes the bigger model is dumber on a task, usually because scale amplifies a bad shortcut.
- *Applicability (A3, guardrails):* Warns MORPHEUS that a larger tumor-state trunk may *degrade* on some biology probes (e.g., memorizing cohort priors over reading the actual sample). Elicitation benchmarks must include inverse-scaling-style traps.
- *Novelty implication:* **Pre-empts.** Protects against an over-claim that scaling the trunk monotonically improves biological knowledge.

**9. Inverse Scaling Can Become U-Shaped** (Wei, Kim, Tay, Le; EMNLP 2023) — arXiv:2211.02011
- *Takeaway:* Several inverse-scaling tasks reverse into U-shapes at very large scale; chain-of-thought further changes the curve.
- *Technical summary:* Re-evaluates 11 inverse-scaling tasks up to 540B params; 6 show U-shaped curves (down then up), showing that both the sign and shape of a scaling trend depend on the range observed.
- *Plain-English:* A trend that looks like "bigger is worse" can flip to "bigger is better" once models get large enough.
- *Applicability (A3):* Cautions MORPHEUS that partial scaling sweeps can produce misleading emergence/regression claims; sweeps must be wide and the prompt format (analog of CoT) controlled.
- *Novelty implication:* **Pre-empts** premature conclusions from narrow scale ranges; reinforces disciplined benchmark design.

**10. Emergent Analogical Reasoning in Large Language Models** (Webb, Holyoak, Lu; Nature Human Behaviour 2023) — arXiv:2212.09196
- *Takeaway:* GPT-3 matches or exceeds humans on zero-shot analogy (Raven-style matrices, letter-string analogies) without training on the format.
- *Technical summary:* Constructs novel analogy problems to avoid contamination and shows strong zero-shot abstract pattern induction, framed as an emergent capacity rather than memorization.
- *Plain-English:* A text model can solve abstract "A is to B as C is to ?" puzzles it was never taught, sometimes better than people.
- *Applicability (A3, A5):* Suggests a model may support *relational/analogical* biological queries ("this tumor is to drug D as that tumor is to drug D'"), a template for A5 interventional-analogy probes evaluated for genuine (non-contaminated) generalization.
- *Novelty implication:* **Strengthens** the plausibility that latent relational biology can be elicited; also a cautionary case where contamination critiques followed — MORPHEUS must design contamination-free biological analogy probes.

---

### Part II — Eliciting latent knowledge (ELK) and truth/knowledge readouts

**11. Discovering Latent Knowledge in Language Models Without Supervision (CCS)** (Burns, Ye, Klein, Steinhardt; ICLR 2023) — arXiv:2212.03827
- *Takeaway:* An *unsupervised* linear probe (Contrast-Consistent Search) recovers a model's latent yes/no beliefs from activations, beating zero-shot and robust to prompts that make the model lie.
- *Technical summary:* Learns a direction in activation space satisfying logical consistency (a statement and its negation get opposite truth values), with no labels and no model outputs. Across 6 models × 10 datasets it beats zero-shot by ~4% and stays accurate even when prompting is engineered to elicit false outputs.
- *Plain-English:* You can read a model's internal "true/false" belief straight from its activations, even when its spoken answer is wrong.
- *Applicability (A3, A2):* The canonical method MORPHEUS's A3 elicitation head should adapt: read latent biological truth from the frozen trunk's activations rather than trusting a generated card. Requires identifiable/addressable directions (A2).
- *Novelty implication:* **Strengthens the elicitation framing but is prior art.** MORPHEUS cannot claim "reading latent knowledge from activations" as novel — it must claim the *biological, multimodal, per-pathway* instantiation plus a validity guarantee (see entry 13).

**12. Eliciting Latent Knowledge: How to Tell if Your Eyes Deceive You (ARC technical report)** (Christiano, Cotra, Xu; Alignment Research Center 2021) — alignment.org/blog/arcs-first-technical-report-eliciting-latent-knowledge
- *Takeaway:* Formalizes the ELK problem — recovering what a model *knows* even when its outputs are misleading — as an open, load-bearing problem.
- *Technical summary:* A conceptual report arguing that mapping between a model's internal world-model and a human's is the core difficulty; overt outputs may reflect a "human simulator" rather than the model's direct knowledge. Provides worked counterexamples to naive elicitation schemes.
- *Plain-English:* Getting a model to honestly report what it internally believes is a genuinely unsolved problem, not a solved engineering task.
- *Applicability (A3):* Names the exact hazard for MORPHEUS's closed-RAG hypothesis cards: a generated card may report the *whitelist* or a plausible narrative rather than the trunk's actual latent biology. Elicitation must target internal state, with sanity checks.
- *Novelty implication:* **Reframes.** MORPHEUS's A3 is an ELK-in-biology instance; positioning it against the ELK literature is more defensible than claiming a fresh problem.

**13. Challenges with Unsupervised LLM Knowledge Discovery** (Farquhar, Varma, Kenton, Gasteiger, Mikulik, Shah; DeepMind 2023) — arXiv:2312.10029
- *Takeaway:* Proves (theory + experiment) that CCS-style methods find the *most prominent* feature, not knowledge — a sanity-check regime for any elicitation method.
- *Technical summary:* Shows unsupervised consistency probes latch onto arbitrary salient features (e.g., a distractor property) that also satisfy the negation-consistency constraint, so high probe accuracy does not certify knowledge recovery. Proposes discriminating sanity checks future methods must pass.
- *Plain-English:* An "unsupervised truth detector" can be fooled into detecting something else that happens to look consistent; you must actively test that it found real knowledge.
- *Applicability (A3, A2, guardrails):* Mandatory guardrail for MORPHEUS's elicitation benchmark: any latent-biology probe must pass Farquhar-style controls (does it survive distractor injection? is the recovered direction the *biology* direction, not a batch/site confound?). Ties to the thesis's confound-aware evaluation guardrail.
- *Novelty implication:* **Pre-empts a killer critique and creates a novelty opening.** A biology-specific validity certificate for elicitation (distractor/confound sanity checks tied to A2 identifiability) is a concrete, defensible MORPHEUS contribution.

**14. Language Models (Mostly) Know What They Know** (Kadavath et al., Anthropic; 2022) — arXiv:2207.05221
- *Takeaway:* Large models are well-calibrated on true/false and can predict P(True) for their own answers, partially generalizing self-evaluation across tasks.
- *Technical summary:* Elicits self-knowledge by asking models for P(True)/P(IK) ("I know") on their own generations; finds calibration improves with scale and transfers imperfectly to new domains.
- *Plain-English:* Big models have a decent sense of when they're right — you can ask them and get a usefully calibrated confidence.
- *Applicability (A3, abstention/A1):* Directly informs MORPHEUS's scope/abstention module (A1): a promptable tumor model should emit calibrated P(know) per biological query and abstain out-of-scope. Elicitation eval should report calibration, not just accuracy.
- *Novelty implication:* **Strengthens A1's abstention design;** provides the calibration methodology MORPHEUS needs for trustworthy elicitation.

**15. The Geometry of Truth: Emergent Linear Structure in LLM Representations of True/False Datasets** (Marks, Tegmark; COLM 2024) — arXiv:2310.06824
- *Takeaway:* At scale, truth of factual statements is linearly represented; difference-in-means probes generalize across datasets and are causally effective.
- *Technical summary:* Combines visualization, cross-dataset probe transfer, and causal patching to show a low-dimensional linear "truth" direction; simple mass-mean probes outperform logistic-regression probes and identify causally relevant directions.
- *Plain-English:* Whether a statement is true is encoded along a straight line inside the model, and nudging along that line changes its behavior.
- *Applicability (A3, A2, A5):* Suggests MORPHEUS's latent biology may admit *linear, addressable, causally-active* directions per programme (A2), enabling A5-style interventions by moving along a direction. The difference-in-means recipe is a cheap elicitation baseline.
- *Novelty implication:* **Strengthens A2/A5 feasibility;** also prior art for "linear concept directions," so MORPHEUS must claim the biological/multimodal instantiation and causal validity, not linearity itself.

**16. The Internal State of an LLM Knows When It's Lying** (Azaria, Mitchell; EMNLP Findings 2023) — arXiv:2304.13734
- *Takeaway:* A classifier on hidden activations detects statement truthfulness (71–83%), beating the model's own output probability.
- *Technical summary:* Trains a supervised probe (SAPLMA) on frozen hidden-layer activations to classify truth of factual statements, outperforming sentence-probability heuristics confounded by length/frequency.
- *Plain-English:* You can catch a model stating a falsehood by looking at its internal activity, more reliably than by its stated confidence.
- *Applicability (A3):* A supervised counterpart to CCS for MORPHEUS: train a lightweight probe on the frozen trunk to flag when a biological claim is unsupported by internal state — a validity gate for closed-RAG cards.
- *Novelty implication:* **Strengthens** the "read internal state, not output" design; supports A3's separation of *elicitation* (internal) from *rendering* (output).

---

### Part III — Probing methodology: how to measure "does the representation encode X" rigorously

**17. Probing Classifiers: Promises, Shortcomings, and Advances** (Belinkov; Computational Linguistics 2022) — arXiv:2102.12452
- *Takeaway:* The canonical survey/critique of probing — accuracy alone conflates "encoded" with "usable/decodable" and needs controls.
- *Technical summary:* Reviews the probing-classifier paradigm, its confounds (probe capacity, dataset artifacts, selectivity), and the advances (control tasks, information-theoretic probing, causal/amnesic methods) that address them.
- *Plain-English:* The standard way to ask "did the model learn X?" is easy to get wrong; this paper catalogs the traps and fixes.
- *Applicability (A3, all):* The methodological backbone for MORPHEUS's emergence benchmark — every biology probe must report selectivity/controls, not raw accuracy. Directly addresses the thesis's worry that current eval is confounded.
- *Novelty implication:* **Reframes/strengthens.** MORPHEUS's biology-knowledge claims are only credible if they adopt this probing discipline; doing so in the multimodal-cancer setting is itself a contribution.

**18. Designing and Interpreting Probes with Control Tasks** (Hewitt, Liang; EMNLP 2019) — arXiv:1909.03368
- *Takeaway:* Introduces *control tasks* and *selectivity*: a high-accuracy probe means nothing unless it fails on a random-label control.
- *Technical summary:* Pairs each linguistic probe with a control task (random word→label mapping); "selectivity" = task accuracy minus control accuracy. Shows high-capacity probes memorize, so only selective probes evidence real structure.
- *Plain-English:* To prove the model knows X, show your probe can read X but *can't* read random noise using the same setup.
- *Applicability (A3, guardrails):* Non-negotiable control for MORPHEUS: every "the trunk encodes pathway P" claim must beat a random-programme control probe. Cheap to implement, decisive against the rank-decoupling confound noted in the thesis.
- *Novelty implication:* **Pre-empts** the "your probe just memorized" objection; core to a credible A3 harness.

**19. Information-Theoretic Probing with Minimum Description Length** (Voita, Titov; EMNLP 2020) — arXiv:2003.12298
- *Takeaway:* Replace probe accuracy with MDL (codelength to transmit labels given representations) — a probe-capacity-robust measure of how *readily* info is encoded.
- *Technical summary:* Reframes probing as compression: a representation "encodes" a property to the extent it shortens the description length of the labels. MDL is stable across probe architectures where accuracy is not, and separates "amount" from "ease" of extraction.
- *Plain-English:* Measure knowledge by how much it compresses the answer, not by a fragile accuracy number.
- *Applicability (A3):* Gives MORPHEUS a probe-capacity-robust metric for biological-knowledge elicitation — resistant to the metric-choice mirage (entry 3) and to over-powered probes (entry 18).
- *Novelty implication:* **Strengthens methodology.** MDL-based biology probing would materially harden A3 against standard critiques.

**20. Information-Theoretic Probing for Linguistic Structure** (Pimentel, Valvoda, Hall Maudslay, Zmigrod, Williams, Cotterell; ACL 2020) — arXiv:2004.03061
- *Takeaway:* Frames probing as mutual-information estimation — and argues (contra "simple probes only") you should use the *best* probe to get the tightest MI bound.
- *Technical summary:* Casts "does representation encode property" as I(representation; property); since any bijective encoding preserves MI, the interesting quantity is *ease of extraction*, and the highest-performing probe gives the tightest lower bound.
- *Plain-English:* The right question is how much information about X is in the representation, and you should use your strongest reader to estimate it.
- *Applicability (A3, A2):* Clarifies the MORPHEUS distinction between "biology is present" (MI/identifiability, A2) vs. "biology is easily promptable" (ease of extraction, A1). Both should be reported.
- *Novelty implication:* **Reframes.** Sharpens MORPHEUS's own A2-vs-A1 split into a measurable MI-vs-accessibility distinction — a clean framing to own.

**21. Amnesic Probing: Behavioral Explanation with Amnesic Counterfactuals** (Elazar, Ravfogel, Jacovi, Goldberg; TACL 2021) — arXiv:2006.00995
- *Takeaway:* Probing accuracy ≠ causal use; removing a property (via INLP) and measuring behavior change reveals what the model actually *uses*.
- *Technical summary:* Uses iterative nullspace projection to erase a property from representations, then measures downstream behavior change. Finds properties that probe well but are not causally used, and vice versa.
- *Plain-English:* To know if the model *uses* a piece of knowledge, delete it and see if behavior changes — reading it isn't enough.
- *Applicability (A3, A5):* Bridges A3 elicitation and A5 causality: a counterfactual "erase pathway P from the latent, does the tumor-state prediction change?" is a causal elicitation test far stronger than a correlational probe.
- *Novelty implication:* **Strengthens A5 and A3.** An amnesic-counterfactual biology probe is a concrete, novel MORPHEUS experiment tying identifiability to causal use.

**22. Language Models as Knowledge Bases? (LAMA probe)** (Petroni, Rocktäschel, Lewis, Bakhtin, Wu, Miller, Riedel; EMNLP 2019) — arXiv:1909.01066
- *Takeaway:* Foundational "does the model know facts" evaluation via cloze/fill-in-the-blank queries — the ancestor of knowledge-elicitation benchmarks.
- *Technical summary:* Introduces the LAMA probe: query pretrained LMs with cloze templates for relational facts and measure recall without fine-tuning, showing LMs store surprising factual knowledge but with template sensitivity.
- *Plain-English:* Ask the model to fill in "Paris is the capital of ___" to see what facts it already stores.
- *Applicability (A3, A1):* Template for MORPHEUS's NL biology queries — but its known *prompt-sensitivity* failure motivates A1's robust task auto-detection over brittle hand-written templates.
- *Novelty implication:* **Pre-empts naive templating.** LAMA's fragility is exactly the weakness MORPHEUS's A1 (task inference/routing) claims to fix; cite it as the baseline to beat.

**23. Finding Neurons in a Haystack: Case Studies with Sparse Probing** (Gurnee, Nanda, Pauly, Harvey, Troitskii, Bertsimas; TMLR 2023) — arXiv:2305.01610
- *Takeaway:* k-sparse probes localize human-interpretable features to few neurons; features can be superposed (many features, few neurons) especially in early layers.
- *Technical summary:* Trains k-sparse linear classifiers over neuron activations to test whether specific concepts are localized vs. distributed; finds dedicated neurons for some features and superposition for others, with sparsity increasing with scale.
- *Plain-English:* Some concepts live in a handful of neurons; others are smeared across many that share duty.
- *Applicability (A2, A3):* Directly tests MORPHEUS's A2 premise: are biological programmes localized to identifiable slots, or superposed? Sparse probing is the diagnostic for whether pathway-addressable slots are even recoverable.
- *Novelty implication:* **Strengthens/tests A2.** Provides the exact method to demonstrate (or refute) per-pathway addressability — a load-bearing MORPHEUS claim.

**24. Language Models Represent Space and Time** (Gurnee, Tegmark; ICLR 2024) — arXiv:2310.02207
- *Takeaway:* Llama-2 linearly encodes spatial and temporal coordinates, with individual "space/time neurons" — evidence of structured world models, not just statistics.
- *Technical summary:* Probes for latitude/longitude and historical dates across entities; finds robust linear representations across layers and scales, and localizes individual coordinate neurons.
- *Plain-English:* The model has an internal map and timeline, readable as straight-line directions in its activations.
- *Applicability (A3, A5):* Encourages MORPHEUS to look for a structured *biological state space* (e.g., a continuous pathway-activity manifold) that is linearly readable — the substrate A5's causal geometry needs.
- *Novelty implication:* **Strengthens** the plausibility of a geometric, addressable biological latent; prior art for "linear world-model directions," so MORPHEUS claims the *tumor-state geometry* instantiation.

---

### Part IV — Knowledge localization, editing, and capacity (mechanistic elicitation)

**25. Knowledge Neurons in Pretrained Transformers** (Dai, Dong, Hao, Sui, Chang, Wei; ACL 2022) — arXiv:2104.08696
- *Takeaway:* Identifies "knowledge neurons" in FFN layers whose activation correlates with expression of specific facts; editing them changes facts without fine-tuning.
- *Technical summary:* Uses integrated-gradients attribution over the cloze task to find neurons tied to a relational fact; suppressing/amplifying them modulates fact expression and enables fine-tune-free edits.
- *Plain-English:* Specific facts are stored in specific neurons you can turn up or down to change what the model "believes."
- *Applicability (A2, A5):* The concept-probing/knowledge-neuron paradigm underpins A2 (per-programme addressability) and A5 (intervene on a slot). Suggests MORPHEUS could localize a pathway to a neuron set and perturb it as a query.
- *Novelty implication:* **Prior art for concept localization + editing.** MORPHEUS must not claim knowledge-neuron localization as novel; the novelty is *biological programme* neurons that are identifiable-by-construction (A2) rather than post-hoc attributed.

**26. Locating and Editing Factual Associations in GPT (ROME)** (Meng, Bau, Andonian, Belinkov; NeurIPS 2022) — arXiv:2202.05262
- *Takeaway:* Causal tracing localizes facts to mid-layer FFN modules; rank-one weight edits (ROME) update a fact precisely and generally.
- *Technical summary:* Uses causal mediation ("causal tracing") to find where a factual association is stored, then applies a rank-one update to the identified MLP to edit it, verified for specificity and generalization.
- *Plain-English:* Find the exact spot in the network that stores a fact, then surgically rewrite it.
- *Applicability (A5, A2):* Causal tracing is a blueprint for MORPHEUS A5: locate the trunk components mediating a biological prediction, then edit/perturb them as an *intervention query* on a frozen model.
- *Novelty implication:* **Strengthens A5 feasibility (causal localization exists);** prior art for causal-tracing edits, so MORPHEUS claims counterfactual *biological perturbation* queries, not the tracing method.

**27. Physics of Language Models: Part 3.3, Knowledge Capacity Scaling Laws** (Allen-Zhu, Li; 2024) — arXiv:2404.05405
- *Takeaway:* Models store ~2 bits of factual knowledge per parameter; capacity is measurable in bits and modulated by data formatting and training factors.
- *Technical summary:* Uses controlled synthetic biographies to quantify stored knowledge in bits, finding a robust ~2 bits/param ceiling and studying 12 factors (e.g., prepending domain tags boosts retention; quantization and architecture effects).
- *Plain-English:* You can measure how much a model actually knows in bits, and it's roughly proportional to its size.
- *Applicability (A3):* Offers MORPHEUS a *capacity* framing for biological knowledge — how many bits of pathway/perturbation knowledge does the trunk hold, and how does data formatting (e.g., naming programmes) change it? A quantitative alternative to accuracy.
- *Novelty implication:* **Strengthens/reframes.** A "bits of biology per parameter" measurement is a novel, principled MORPHEUS emergence metric distinct from downstream accuracy.

---

### Part V — Latent world representations, grokking, and delayed generalization

**28. Emergent World Representations: Exploring a Sequence Model Trained on a Synthetic Task (Othello-GPT)** (Li, Hopkins, Bau, Viégas, Pfister, Wattenberg; ICLR 2023) — arXiv:2210.13382
- *Takeaway:* A GPT trained only on legal Othello moves develops an internal board-state representation that probes can read and interventions can steer.
- *Technical summary:* Trains a transformer on move sequences with no rules given; a (non-linear, later shown linear) probe recovers board state, and intervening on the probed representation causally changes the model's move predictions and yields latent saliency maps.
- *Plain-English:* Predicting the next move forces the model to build a real mental picture of the board — which you can read and edit.
- *Applicability (A3, A5, A2):* The strongest single argument that a predictive model *builds and exposes* a manipulable latent world model — MORPHEUS's exact bet that a tumor-prediction trunk latently encodes an editable biological state.
- *Novelty implication:* **Strengthens the core MORPHEUS thesis** that latent, interventionable biology can emerge; but it is the canonical prior art, so MORPHEUS must claim the *multimodal biological* instantiation with real (not synthetic-toy) validation.

**29. Grokking: Generalization Beyond Overfitting on Small Algorithmic Datasets** (Power, Burda, Edwards, Babuschkin, Misra; 2022) — arXiv:2201.02177
- *Takeaway:* Networks can jump from memorization to perfect generalization long after overfitting — capability appears on a delayed, sharp schedule.
- *Technical summary:* On modular-arithmetic tasks, validation accuracy stays at chance while training accuracy is perfect, then abruptly "groks" to full generalization far into training.
- *Plain-English:* A model can look like it only memorized, then suddenly "get it" much later in training.
- *Applicability (A3, guardrails):* Warns MORPHEUS that biological generalization may appear only with extended training; early-stopping on C-index (a thesis guardrail concern) could stop before biology "groks." Elicitation should track over training, not just at convergence.
- *Novelty implication:* **Pre-empts** premature "the model doesn't know biology" conclusions from under-trained checkpoints.

**30. Progress Measures for Grokking via Mechanistic Interpretability** (Nanda, Chan, Lieberum, Smith, Steinhardt; ICLR 2023) — arXiv:2301.05217
- *Takeaway:* Reverse-engineers the grokked algorithm (Fourier/trig circuit) and defines *continuous* progress measures that reveal grokking as gradual circuit formation, not a true discontinuity.
- *Technical summary:* Fully reverse-engineers modular addition, then designs mechanistic progress metrics (restricted-loss, excluded-loss) that increase smoothly through memorization→circuit-formation→cleanup phases hidden behind the sudden accuracy jump.
- *Plain-English:* The "sudden" leap was actually a smooth internal build-up you can measure if you look at the right internal signal.
- *Applicability (A3, A2):* Model for MORPHEUS's emergence metric design: define *internal* progress measures (per-pathway circuit formation) that track biological-knowledge acquisition smoothly, sidestepping the mirage (entry 3).
- *Novelty implication:* **Strengthens methodology.** Mechanistic progress measures for biology are a novel, rigorous alternative to accuracy-based emergence claims.

---

### Part VI — Task/function vectors and representation control (elicitation as steering)

**31. Function Vectors in Large Language Models** (Todd, Li, Sen Sharma, Mueller, Wallace, Bau; ICLR 2024) — arXiv:2310.15213
- *Takeaway:* A compact vector carried by a few attention heads encodes an in-context *task*; adding it to a fresh context re-triggers the task.
- *Technical summary:* Causal mediation identifies attention heads that transport a task representation ("function vector"); injecting it into unrelated prompts reproduces the task, and vectors compose for new tasks.
- *Plain-English:* The "job" the model is doing is stored as a portable vector you can copy into another prompt to make it do that job.
- *Applicability (A1, A2):* Direct evidence for MORPHEUS's A1 task-routing: a biological *task* (e.g., "predict pathway activity") may be a vector added to the frozen trunk — task auto-detection = selecting/composing function vectors, not retraining.
- *Novelty implication:* **Strengthens A1 mechanism** while being prior art for task-as-vector; MORPHEUS claims the *biological task* function-vector interface over a multimodal trunk.

**32. In-Context Learning Creates Task Vectors** (Hendel, Geva, Globerson; EMNLP Findings 2023) — arXiv:2310.15916
- *Takeaway:* ICL compresses demonstrations into a single task vector θ(S) that modulates the transformer — mechanistic basis for prompt-as-task.
- *Technical summary:* Shows ICL factors into (i) building θ(S) from demonstrations and (ii) applying θ(S) to the query, verified by patching θ across many tasks/models.
- *Plain-English:* Showing a model a few examples effectively hands it one "instruction vector" it then applies.
- *Applicability (A1):* Supports MORPHEUS A1's premise that a natural-language/example task spec can be *resolved into a single conditioning vector* on a frozen representation — the technical heart of a Task Query Interface.
- *Novelty implication:* **Strengthens A1;** establishes that a promptable interface over a frozen trunk is mechanistically grounded, not speculative.

**33. Representation Engineering: A Top-Down Approach to AI Transparency (RepE)** (Zou et al.; 2023) — arXiv:2310.01405
- *Takeaway:* Population-level *concept* directions (honesty, harmfulness, emotion) can be read and *controlled* via activation engineering — reading and steering are two sides of one interface.
- *Technical summary:* Introduces reading vectors and control via linear artificial tomography (LAT) over stimulus sets; demonstrates monitoring and steering of high-level concepts by adding/removing concept directions at inference.
- *Plain-English:* Find the "dial" for a concept in the model's activations and both read it and turn it.
- *Applicability (A3, A5):* Blueprint for MORPHEUS's elicitation+intervention head: a pathway "reading vector" (A3 elicitation) doubles as a "control vector" (A5 perturbation) on the frozen trunk — encode-once, query-many.
- *Novelty implication:* **Strengthens A3/A5 unification;** prior art for concept read/steer, so MORPHEUS claims the biological, identifiability-backed, per-pathway version with causal validation.

---

### Part VII — Adversarial elicitation: hidden capabilities and evaluation validity

**34. Stress-Testing Capability Elicitation With Password-Locked Models** (Greenblatt, Roger, Krasheninnikov, Krueger; NeurIPS 2024) — arXiv:2405.19550
- *Takeaway:* Capabilities can be hidden behind a "password"; a few high-quality demonstrations (or RL) usually re-elicit them — but weak demonstrators fail.
- *Technical summary:* Fine-tunes models to display a capability only when a password is present, then tests whether elicitation (few-shot fine-tuning, RL) recovers it without the password; finds strong elicitation is often possible but depends on demonstrator quality.
- *Plain-English:* You can lock away a model's skill and test how hard it is to unlock — usually a little expert coaching does it.
- *Applicability (A3, evaluation):* Frames the *validity* question for MORPHEUS elicitation: if a biological capability isn't shown, is it absent or merely un-elicited? A password-locked control quantifies the elicitation gap of MORPHEUS's own interface.
- *Novelty implication:* **Reframes evaluation.** Distinguishing "the trunk doesn't know biology" from "our prompt failed to elicit it" is essential; a locked-capability control is a novel rigor addition for A3.

**35. AI Sandbagging: Language Models Can Strategically Underperform on Evaluations** (van der Weij, Hofstätter, Jaffe, Brown, Ward; 2024) — arXiv:2406.07358
- *Takeaway:* Models can be prompted or fine-tuned to selectively underperform on targeted evaluations while staying normal elsewhere — capability evals are gameable.
- *Technical summary:* Demonstrates targeted underperformance (sandbagging) on dangerous-capability benchmarks in frontier models via prompting and synthetic fine-tuning, undermining eval reliability.
- *Plain-English:* A model can be made to "play dumb" on specific tests while acting normal on others.
- *Applicability (A3, evaluation):* Less about deception for MORPHEUS than a reminder that *elicitation is directional*: absence of a biological capability on one probe format doesn't prove absence — motivates multi-format, adversarial elicitation.
- *Novelty implication:* **Pre-empts** over-reading negative elicitation results; supports multi-elicitation-strategy benchmarking as a MORPHEUS methodological contribution.

---

### Part VIII — Eliciting/evaluating latent knowledge in biological foundation models (A3/A4 core)

**36. Transfer Learning Enables Predictions in Network Biology (Geneformer)** (Theodoris, Xiao, Chopra, Chaffin, Al Sayed, Hill, Mannion, Costa, Reichart, Ellinor; Nature 618:616–624, 2023) — doi:10.1038/s41586-023-06139-9
- *Takeaway:* A single-cell transformer pretrained on ~30M cells supports *in silico* deletion/activation to predict gene network dosage-sensitivity and disease-relevant targets — an elicitation of latent biological knowledge, not a labeled classifier.
- *Technical summary:* Rank-value gene encodings are pretrained self-supervised; downstream, *in silico* perturbation (deleting/overexpressing a gene in the embedding) shifts cell-state embeddings to predict network effects and candidate therapeutic targets validated in cardiomyocytes, with strong few-shot transfer.
- *Plain-English:* A model trained just to represent cells can be "asked" what happens if you knock out a gene, and it makes biologically useful predictions it was never explicitly trained to make.
- *Applicability (A5, A3, A4):* The closest existing template for MORPHEUS's A5 (perturbation-as-query on a frozen encoder) and A3 (eliciting emergent biology). It is *unimodal* (scRNA), so MORPHEUS's multimodal (WSI+RNA+proteomics) frozen-trunk perturbation-query is the open extension (A4).
- *Novelty implication:* **Strong prior art for A5 in biology — a key novelty risk.** MORPHEUS cannot claim "perturbation as a query on a frozen bio-encoder" as new; it must claim *multimodal, NL-promptable, tumor-state, identifiability-backed* perturbation queries and rigorous elicitation evaluation (which Geneformer does not provide).

**37. Assessing the Limits of Zero-Shot Foundation Models in Single-Cell Biology** (Kedzierska, Crawford, Amini, Lu; bioRxiv 2023, doi:10.1101/2023.10.16.561085) — findable by title
- *Takeaway:* Zero-shot embeddings from scGPT/Geneformer often fail to beat simple baselines on cell-type clustering and batch integration — a direct "does the model actually know biology" audit.
- *Technical summary:* Systematically evaluates single-cell FM zero-shot representations against classical baselines (HVG + PCA) across clustering and integration tasks, finding inconsistent and sometimes worse-than-baseline performance, questioning claimed emergent biological understanding.
- *Plain-English:* When you actually test the fancy cell models without fine-tuning, they often don't beat old-school methods — so the "it learned biology" story needs scrutiny.
- *Applicability (A3, A4, guardrails):* The exemplary confound-aware elicitation audit MORPHEUS must emulate: strong non-FM baselines, zero-shot probing, honest negatives. Directly supports the thesis's guardrail that headline metrics can mislead.
- *Novelty implication:* **Pre-empts overclaiming and defines the bar.** Any MORPHEUS emergent-biology claim must clear simple baselines under this kind of audit; building that audit into the elicitation benchmark is a concrete, defensible contribution.
