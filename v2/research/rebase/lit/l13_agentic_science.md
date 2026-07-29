## Agentic LLM scientific workflows & tool-use

Lane l13. Remit: LLM agents that plan/execute scientific tasks; tool-use over scientific models; LLM-as-orchestrator for biology; model-as-tool APIs; function-calling to specialist predictors. Primary MORPHEUS relevance: **A1** (NL task inference/routing) and **A3** (NL<->biology grounding + emergent-knowledge elicitation & its evaluation), with secondary touches on A4/A5.

**Cross-cutting novelty tension for MORPHEUS (read first).** The dominant paradigm in this lane is *orchestrator + external tools*: an LLM parses a natural-language request, infers the task, and **routes it to a frozen specialist predictor called as a tool** (Gorilla, ToolLLM, HuggingGPT, TxAgent, Biomni, CRISPR-GPT). This directly pre-empts any naive MORPHEUS claim that "NL-driven task auto-detection over biological models" is itself novel — agents already do exactly that at the *system* level. MORPHEUS's defensible wedge is therefore **internalizing** the task routing inside a single promptable trunk (one model, many programmes, addressable slots) rather than gluing many models together with LLM control code — i.e., novelty lives in A2 (identifiable, pathway-addressable slots) and the frozen-trunk plug-in of A4, NOT in "an LLM picks the task." Every entry below is scored against that tension.

---

### 1. ReAct: Synergizing Reasoning and Acting in Language Models
Yao, Zhao, Yu, Du, Shafran, Narasimhan, Cao — ICLR 2023. arXiv:2210.03629.
**Takeaway:** Interleaving chain-of-thought reasoning traces with tool/environment actions is the base recipe that makes LLM task-routing reliable and inspectable.
**Technical summary:** ReAct prompts an LLM to emit alternating "thought" and "action" tokens, where actions call an external API (e.g., Wikipedia) and observations feed back into context. On HotpotQA/Fever it curbs CoT hallucination via grounding, and on ALFWorld/WebShop it beats imitation/RL baselines by 34%/10% absolute with 1–2 in-context examples. The reasoning trace is the mechanism that decides *which* tool to call and *how* to fold results back.
**Plain-English:** Let the model "think out loud," take an action in the world, look at the result, and think again — this loop is far more reliable than either pure reasoning or pure action.
**Applicability:** A1 (the canonical think→route→act loop for task inference), A3 (interpretable NL traces ground each decision). Design implication: any MORPHEUS "task auto-detection" front-end should expose a ReAct-style trace so routing decisions are auditable, not a black-box classifier head.
**Novelty implication:** Pre-empts "LLM infers the task then acts" as a generic contribution; MORPHEUS must claim the *unified-trunk* internalization, not the reasoning loop.

### 2. Toolformer: Language Models Can Teach Themselves to Use Tools
Schick, Dwivedi-Yu, Dessì, et al. — NeurIPS 2023. arXiv:2302.04761.
**Takeaway:** A model can self-supervise *when/which/how* to call tools by filtering its own API-call insertions on downstream perplexity.
**Technical summary:** Toolformer samples candidate API calls (calculator, QA, search, translate, calendar) inside training text, keeps only calls that reduce next-token loss, and fine-tunes on the filtered corpus. The result decides autonomously which API to call, with what arguments, and how to splice the return — improving zero-shot downstream tasks without task-specific demonstrations. Tool-use becomes a learned property of the weights, not a prompt scaffold.
**Plain-English:** The model learns on its own that "here I should call a calculator" by checking whether calling it would have helped it predict the text.
**Applicability:** A1 (learned, not hard-coded, tool/task selection), A4 (a template for *learning* when a modality-specialist should be invoked vs. handled in-context). Design implication: supports MORPHEUS's "learned encode-vs-RAG gating" — the decision to encode a modality can be trained by held-out utility, mirroring Toolformer's perplexity filter.
**Novelty implication:** Reframes A4: "when to encode vs. retrieve" is learnable self-supervision, strengthening MORPHEUS if it makes that gate *internal* rather than an external agent policy.

### 3. Gorilla: Large Language Model Connected with Massive APIs
Patil, Zhang, Wang, Gonzalez — NeurIPS 2024. arXiv:2305.15334.
**Takeaway:** The clearest "model-as-tool" precedent: NL instruction → correct API call into ML model zoos (HuggingFace/TorchHub/TensorHub), with retrieval to stay current.
**Technical summary:** Gorilla fine-tunes LLaMA with retriever-aware training over APIBench, mapping instructions like "build me a classifier for medical images" to the exact model API + dependencies. Retrieval-aware fine-tuning lets it adapt to changing/versioned APIs and cuts hallucinated calls, surpassing GPT-4 on API-writing accuracy. It formalizes calling *specialist models* as callable tools selected from documentation.
**Plain-English:** You ask in English for a capability, and Gorilla writes the exact code to call the right pretrained model for the job.
**Applicability:** A1 (NL→specialist-model routing is the whole point), A4 (models-as-tools is precisely the "treat a predictor as a callable" pattern). Design implication: this is the strongest steelman that "NL routes to biological predictors" is solved externally; MORPHEUS must show that *one addressable trunk* beats a retriever-over-model-zoo on identifiability/consistency.
**Novelty implication:** **Pre-empts** the router-level A1 claim; pushes MORPHEUS novelty toward A2 addressable slots and cross-task weight sharing.

### 4. ToolLLM: Facilitating LLMs to Master 16000+ Real-world APIs
Qin, Liang, Ye, et al. (incl. Gerstein) — ICLR 2024. arXiv:2307.16789.
**Takeaway:** Scales tool-use to 16k+ real APIs with a DFS decision-tree over reasoning traces; open ToolLLaMA matches ChatGPT on complex tool orchestration.
**Technical summary:** ToolBench auto-generates instructions and multi-tool solution paths across 16,464 RapidAPI endpoints; a DFSDT search expands multiple reasoning traces to find valid tool chains. ToolLLaMA (LLaMA + neural API retriever) generalizes zero-shot to unseen APIs and out-of-distribution APIBench. Demonstrates that *breadth* of callable tools is tractable with retrieval + search.
**Plain-English:** Trains an open model to correctly chain thousands of web APIs to satisfy a request, rivaling ChatGPT.
**Applicability:** A1 (large-scale task→tool routing), A3 (instruction grounding). Design implication: at 16k tools, retrieval + search dominate — a caution that a MORPHEUS "many biological programmes" interface must keep the *addressable set identifiable* to avoid the same search-cost blow-up.
**Novelty implication:** Reinforces the router paradigm; MORPHEUS differentiates by bounded, named, per-pathway slots vs. an open retrieval pool.

### 5. Reflexion: Language Agents with Verbal Reinforcement Learning
Shinn, Cassano, Berman, Gopinath, Narasimhan, Yao — NeurIPS 2023. arXiv:2303.11366.
**Takeaway:** Agents improve by writing natural-language self-critiques into episodic memory instead of updating weights.
**Technical summary:** After a failed trajectory the agent verbally reflects on feedback, stores the reflection, and conditions future attempts on it, lifting HumanEval to 91% (vs. 80% GPT-4 baseline) and improving decision/reasoning tasks. Learning is linguistic and test-time, requiring no gradient updates. It shows self-verification/critique as a reusable agent primitive.
**Plain-English:** The agent writes itself a note about what went wrong and reads it before trying again, getting better without retraining.
**Applicability:** A3 (self-critique as an elicitation/verification mechanism), A1 (iterative task refinement). Design implication: MORPHEUS's emergent-knowledge *evaluation* can borrow Reflexion-style verbal verification loops to probe whether the trunk "knows" a biological fact under challenge.
**Novelty implication:** Neutral/strengthens A3 evaluation tooling; not a biology-specific claim.

### 6. Voyager: An Open-Ended Embodied Agent with LLMs
Wang, Xie, Jiang, et al. — 2023. arXiv:2305.16291.
**Takeaway:** An LLM agent that builds a reusable *skill library* of executable code, composing prior skills into new tasks — a template for accumulating callable capabilities.
**Technical summary:** In Minecraft, Voyager uses an automatic curriculum, a growing library of verified code "skills," and iterative prompting with environment feedback to acquire and recombine abilities lifelong. It outperforms prior agents on exploration and skill transfer without fine-tuning. Skills are stored as callable functions and retrieved for novel goals.
**Plain-English:** The agent writes little programs for things it learns to do, saves them, and snaps them together to solve new problems.
**Applicability:** A1 (compositional task solving), A4 (skill library ≈ callable specialist inventory). Design implication: argues for a *composable* MORPHEUS interface where pathway-addressable slots act like reusable skills, not one-shot heads.
**Novelty implication:** Reframes A2 positively — addressable programmes as a persistent skill inventory is a coherent, differentiated framing vs. ad hoc probes.

### 7. Autonomous Chemical Research with Large Language Models (Coscientist)
Boiko, MacKnight, Kline, Gomes — Nature 2023, 624:570. doi:10.1038/s41586-023-06792-0.
**Takeaway:** GPT-4 orchestrating web search, docs, code, and lab automation autonomously designs and runs real chemistry (Pd cross-couplings).
**Technical summary:** Coscientist binds multiple LLM modules (planner, web searcher, doc reader, code executor, automation controller) to plan and physically execute experiments, including reaction optimization on a cloud lab. It demonstrates (semi-)autonomous design→execute→analyze cycles across six tasks. Tool-use bridges language reasoning to wet-lab actuation.
**Plain-English:** An AI reads the literature, writes the code, and drives real lab robots to run and optimize chemical reactions largely on its own.
**Applicability:** A1 (end-to-end task inference→execution), A5 (proposes and *runs* interventions — experiments as queries). Design implication: the canonical proof that "intervention-as-a-query" can be operationalized by an orchestrator; MORPHEUS's A5 counterfactual-query claim must beat this by not needing a retrained/rerun loop.
**Novelty implication:** **Pre-empts** "agent proposes interventions" at the system level; MORPHEUS must claim *in-model* counterfactual geometry vs. external experiment loops.

### 8. Emergent Autonomous Scientific Research Capabilities of LLMs
Boiko, MacKnight, Gomes — 2023. arXiv:2304.05332.
**Takeaway:** The precursor study documenting that a multi-tool LLM agent exhibits emergent planning/execution of chemical research.
**Technical summary:** A GPT-4 agent with web-search, code, and documentation tools plans multi-step syntheses and reasons about experimental design, showing capabilities that appear only when tools + reasoning are combined. It frames these as *emergent* research skills of tool-augmented LLMs. Establishes the empirical basis later hardened into Coscientist.
**Plain-English:** When you give a strong LLM the right tools, it starts doing genuine multi-step research planning that it can't do from text alone.
**Applicability:** A3 (emergent capability measurement), A1. Design implication: motivates MORPHEUS's A3 emphasis on *measuring* emergent knowledge — but warns that "emergence" here is a system property of tools+LLM, so MORPHEUS must isolate emergence *of the trunk itself*.
**Novelty implication:** Cautions A3 evaluation: attribute emergence to the model, not the scaffold; strengthens the case for careful ablation.

### 9. Augmenting Large Language Models with Chemistry Tools (ChemCrow)
Bran, Cox, Schilter, Baldassari, White, Schwaller — Nature Machine Intelligence 2024, 6:525. arXiv:2304.05376.
**Takeaway:** GPT-4 + 18 expert chemistry tools autonomously plans/executes syntheses and guides discovery of a novel chromophore.
**Technical summary:** ChemCrow wraps 18 curated cheminformatics/synthesis tools behind a ReAct-style agent that plans organic synthesis, drug/material tasks, and reasons over tool outputs. Expert + LLM evaluation shows it automates diverse chemical tasks and physically synthesizes an insect repellent and organocatalysts. Tool encapsulation with strict I/O is what suppresses chemistry hallucination.
**Plain-English:** Bolting real chemistry software onto GPT-4 lets it plan and actually make molecules, rather than bluffing.
**Applicability:** A1 (NL→chemistry-task routing), A4 (specialist tools as callables). Design implication: the "curated tool with strict schema stops hallucination" lesson maps to MORPHEUS's need for *identifiable* slots — vague interfaces hallucinate, constrained ones don't.
**Novelty implication:** Pre-empts NL-routing-to-chemistry-predictors; supports A2 (constraint = reliability) as the real lever.

### 10. Empowering Biomedical Discovery with AI Agents
Gao, Fang, Huang, Giunchiglia, Noori, Schwarz, Ektefaie, Kondic, Zitnik — Cell 2024 (arXiv:2404.02831).
**Takeaway:** The reference vision paper framing biomedical "AI scientists" as LLM-reasoning + biomedical tools + experimental platforms with structured memory.
**Technical summary:** Positions collaborative agents that combine ML predictors, databases, and lab platforms under an LLM controller, using persistent memory for continual learning across virtual-cell simulation, circuit design, and therapy development. Emphasizes human-in-the-loop augmentation and levels of agent autonomy. Provides the taxonomy the field now cites.
**Plain-English:** A blueprint for AI "lab partners" in biomedicine that reason, call specialist tools, remember, and help run experiments.
**Applicability:** A1, A3, A4, A5 (covers the full stack conceptually). Design implication: use its autonomy-levels taxonomy to position MORPHEUS as the *unified perception/representation layer* these agents would call, rather than another orchestrator.
**Novelty implication:** **Reframes** MORPHEUS's role: it can slot in as the callable "virtual-cell" trunk beneath such agents — complementary, not competing, which is a cleaner novelty story.

### 11. Biomni: A General-Purpose Biomedical AI Agent
Huang, Zhong, ... Leskovec — bioRxiv 2025.05.30.656746 (PMC12157518).
**Takeaway:** A generalist agent (LLM reasoning + retrieval-augmented planning + code execution over ~150 tools/databases) that generalizes across many biomedical tasks *without task-specific tuning*.
**Technical summary:** Biomni auto-mines tools/software/databases into a unified action space and uses code-based execution to solve genomics, drug-repurposing, rare-disease diagnosis, microbiome, and cloning tasks. It reports strong zero-prompt-tuning generalization and a no-code web interface used by 15k+ scientists. Demonstrates one agent spanning heterogeneous omics tasks.
**Plain-English:** One AI research assistant that can analyze genomes, design protocols, and prioritize genes across many biology tasks without being reprogrammed each time.
**Applicability:** A1 (general task inference across biology — direct competitor framing), A3, A4. Design implication: Biomni is the closest *general* prior art to MORPHEUS's "one interface, many biological tasks" pitch — MORPHEUS must show that a single representation trunk gives identifiability/consistency guarantees Biomni's tool-glue cannot.
**Novelty implication:** **Strongest pre-emption risk in the lane** for the "generalist biological task router." Steelman it explicitly; differentiate on A2/A4 (encoded unified trunk vs. code-orchestrated toolbox).

### 12. TxAgent: Therapeutic Reasoning Across a Universe of Tools
Gao, Wang, ... Zitnik — 2025. arXiv:2503.10970. (ToolUniverse, 211 tools.)
**Takeaway:** White-box multi-step therapeutic reasoning agent that dynamically selects among 211 verified biomedical tools (openFDA, Open Targets, Monarch/HPO).
**Technical summary:** TxAgent interleaves reasoning with real-time tool calls for drug interactions, contraindications, and patient-specific therapy, hitting 92.1% on open-ended drug reasoning and beating GPT-4o by up to 25.8% and DeepSeek-R1-671B on structured multi-step tasks. ToolUniverse standardizes tool schemas for reliable selection. Real-time retrieval keeps answers grounded in verified sources.
**Plain-English:** A drug-reasoning AI that looks things up in trusted medical databases step by step before recommending treatments.
**Applicability:** A1 (dynamic tool selection), A3 (grounding in verified biology), A4 (databases as RAG context vs. encoded). Design implication: TxAgent's "verified-source retrieval" is the A4 *RAG side* — MORPHEUS should adopt its schema discipline for the modalities it chooses to retrieve rather than encode.
**Novelty implication:** Reframes A4: some biomedical knowledge is better as governed RAG than encoded — MORPHEUS's contribution is *deciding* the split, not encoding everything.

### 13. CRISPR-GPT: Agentic Automation of Gene-Editing Experiments
Huang, Qu, Cheng, ... Cong, Zou, et al. — Nature Biomedical Engineering 2025. arXiv:2404.18021.
**Takeaway:** Multi-agent LLM copilot (User Proxy, Planner, Task Executor, Tool Provider) that designs full CRISPR experiments across four editing modalities and is wet-lab validated.
**Technical summary:** CRISPR-GPT selects CRISPR systems, designs gRNAs, recommends delivery, drafts protocols, and plans validation for knockout/base/prime/epigenetic editing, with Meta/Auto/QA interaction modes. It was validated by knocking out four genes (Cas12a) and epigenetically activating two (dCas9) in human cell lines. Domain tools + planner suppress design errors.
**Plain-English:** An AI copilot that takes a gene-editing goal and produces the complete, validated experimental design.
**Applicability:** A1 (task decomposition/routing), A2 (per-modality "slots" ≈ editing modalities), A5 (designs interventions). Design implication: its modality-typed sub-agents mirror MORPHEUS's pathway-addressable slots — evidence that *typed, named* sub-capabilities improve reliability.
**Novelty implication:** Supports A2 framing (typed addressability works); pre-empts "agent designs interventions" for gene editing specifically.

### 14. BioDiscoveryAgent: An AI Agent for Designing Genetic Perturbation Experiments
Roohani, Lee, Huang, ... Leskovec, Zou — ICLR 2025. arXiv:2405.17631.
**Takeaway:** An LLM-only closed-loop agent (with lit/gene-search + AI critic) that beats trained Bayesian-optimization baselines at choosing which genes to perturb next.
**Technical summary:** Without training an ML surrogate or acquisition function, the agent builds prompts from task + prior-round results to pick perturbation targets, using PubMed/Reactome tools and a self-critique. With Claude 3.5 Sonnet it improves relevant-hit prediction by ~21% (and ~46% on hard non-essential genes) over Bayesian optimization across six datasets. Reasoning over biology literature substitutes for a learned acquisition model.
**Plain-English:** Instead of a specialized optimizer, a well-prompted LLM that reads papers picks better gene-knockout experiments round by round.
**Applicability:** A5 (perturbation selection as an interventional query — central), A1, A3. Design implication: shows an LLM's *prior biological knowledge* can drive interventional search — MORPHEUS's A5 counterfactual queries should be benchmarked against this LLM-only baseline, not just against BO.
**Novelty implication:** **Pre-empts** "LLM knowledge guides perturbation choice"; MORPHEUS must show its causal-geometry queries add calibrated, in-representation counterfactuals beyond prompt-driven literature reasoning.

### 15. The Virtual Lab: AI Agents Design New SARS-CoV-2 Nanobodies
Swanson, Wu, Bulaong, Pak, Zou — Nature 2025 (bioRxiv 2024.11.11.623004).
**Takeaway:** A "virtual lab" of role-specialized LLM agents (PI + domain scientists + critic) runs interdisciplinary research meetings and designs nanobodies later validated experimentally.
**Technical summary:** A human sets a high-level goal; an AI PI agent convenes specialist agents (immunologist, ML expert, computational biologist) plus a scientific-critic agent through structured "team meetings," producing a computational pipeline (ESM/AlphaFold-Multimer/Rosetta) that proposes nanobody designs. A subset bound recent SARS-CoV-2 variants in wet-lab assays. Agent-team structure with a critic improves scientific output quality.
**Plain-English:** A simulated team of AI scientists holds meetings, argues, and designs real antibody-like proteins that actually worked in the lab.
**Applicability:** A1 (goal→multi-role decomposition), A3, A5 (design-as-intervention). Design implication: the *critic agent* is the reusable idea — MORPHEUS's A3 emergent-knowledge evaluation benefits from an adversarial critic that stress-tests elicited biology.
**Novelty implication:** Pre-empts "agent team designs biomolecules"; orthogonal to MORPHEUS's representation claim but a strong prior-art anchor for A5.

### 16. BioinformaticsAgent (BIA): LLMs to Reshape Bioinformatics Workflows
Xin, et al. — bioRxiv 2024.05.22.595240.
**Takeaway:** An LLM agent that encapsulates bioinformatics functions (with semantic descriptions + strict parameter constraints) to map analytic intent to correct tool calls without hallucinating parameters.
**Technical summary:** BIA wraps toolkit functions (e.g., Scanpy/Squidpy/Monocle-style) as schema-constrained callables and uses an LLM planner to compose multi-step analysis pipelines. Strict I/O specs prevent invented parameters and let high-level intent become concrete function calls. Targets reproducible, agent-driven omics analysis.
**Plain-English:** An AI that turns "analyze this single-cell dataset" into the exact, valid sequence of bioinformatics commands.
**Applicability:** A1 (intent→pipeline routing), A4 (tools-as-callables). Design implication: reinforces that **schema constraint = reliability**; MORPHEUS's addressable slots should carry explicit type/parameter contracts (A2).
**Novelty implication:** Supports A2; another instance pre-empting generic "NL→omics-task routing."

### 17. AutoBA: Fully Automated Multi-omic Analysis via an LLM Agent
Zhou, Zhang, Chen, Li, Xu, Chen, Gao — Advanced Science 2024 (arXiv:2309.03242).
**Takeaway:** An LLM agent that self-designs end-to-end omics workflows (WGS, RNA-seq, scRNA-seq, ChIP-seq, spatial) from minimal input, locally deployable for privacy.
**Technical summary:** AutoBA generates detailed step-by-step analysis plans and code from a description of the data and goal, adapting workflows to the input's characteristics and to emerging tools. It is robust across sequencing modalities and runs locally. Demonstrates one agent covering many omics pipelines.
**Plain-English:** Tell it what data you have and what you want to learn, and it writes and runs the whole bioinformatics pipeline.
**Applicability:** A1 (auto workflow inference), A4 (multi-modality handling). Design implication: AutoBA's per-modality workflow branching is the *external* analog of MORPHEUS's internal multimodal prompting — useful contrast for the encode-vs-orchestrate argument.
**Novelty implication:** Pre-empts "auto-detect the omics task"; MORPHEUS's edge is shared representation across modalities vs. separate per-modality pipelines.

### 18. GeneAgent: Self-Verification Language Agent for Gene-Set Knowledge Discovery
Wang, Jin, Wei, Tian, Lai, Zhu, Day, Ross, Lu — 2024. arXiv:2405.16205.
**Takeaway:** A gene-set-annotation agent that autonomously queries domain databases to *self-verify* claims, cutting hallucination vs. GPT-4.
**Technical summary:** GeneAgent generates candidate gene-set functions then verifies each against curated biological databases before finalizing, evaluated on 1,106 gene sets where it consistently beats standard GPT-4. Manual review confirms the self-verification module improves factual reliability of the narratives. Grounding-by-verification is the core mechanism.
**Plain-English:** An AI that names what a group of genes does, then double-checks each claim against trusted databases before answering.
**Applicability:** A3 (NL<->biology grounding + hallucination control — central), A1. Design implication: MORPHEUS's A3 emergent-knowledge claims need exactly this *database-grounded verification* to be credible; adopt GeneAgent-style self-checks in the eval harness.
**Novelty implication:** **Strengthens A3 methodology** but warns that ungrounded "emergent knowledge" is easily contaminated by LLM priors — MORPHEUS must show elicited knowledge exceeds a verify-against-DB agent.

### 19. Language Agents Achieve Superhuman Synthesis of Scientific Knowledge (PaperQA2)
Skarlinski, Cox, Laurent, Braza, Hinks, Hammerling, Ponnapati, Rodriques, White — 2024. arXiv:2409.13740.
**Takeaway:** An optimized RAG agent that beats human experts at cited literature summarization and contradiction detection in biology.
**Technical summary:** PaperQA2 runs retrieval, evidence gathering, and iterative answer refinement over full-text literature, producing Wikipedia-grade cited summaries and detecting contradictions across papers (with 70% of detected contradictions validated by experts). Rigorous human-vs-AI protocols show superhuman synthesis on defined tasks. Establishes agentic literature grounding as measurable.
**Plain-English:** An AI that reads the scientific literature and writes better-cited, more accurate summaries than expert scientists — and spots where papers disagree.
**Applicability:** A3 (grounding + evaluation methodology — central), A4 (literature as RAG). Design implication: its contradiction-detection eval is a template for MORPHEUS's A3 emergent-knowledge *evaluation* — test whether the trunk's beliefs contradict curated literature.
**Novelty implication:** **Reframes A4/A3**: for text-expressible biological knowledge, RAG already achieves superhuman synthesis — so MORPHEUS's A3 novelty must be knowledge *not* recoverable from literature RAG (i.e., latent structure only the encoder sees).

### 20. Accelerating Scientific Discovery with the AI Co-Scientist
Gottweis, Weng, Daryin, Tu, et al. (Google) — 2025. arXiv:2502.18864.
**Takeaway:** A tournament-based multi-agent system that evolves and ranks hypotheses via test-time compute scaling, validated on drug repurposing / AMR with experimental confirmation (novel AML candidates).
**Technical summary:** An asynchronous multi-agent architecture generates, debates, and iteratively refines hypotheses through a self-improving tournament (generation, reflection, ranking, evolution agents) that scales compute at inference. It surfaced experimentally confirmed drug-repurposing candidates for AML and mechanisms for AMR. Emphasizes hypothesis *evolution* rather than single-shot generation.
**Plain-English:** An AI that generates many scientific hypotheses, has agents critique and compete them like a tournament, and proposes the survivors — some of which panned out in the lab.
**Applicability:** A3 (hypothesis generation/grounding), A5 (drug-repurposing = intervention proposal). Design implication: a large-scale steelman that hypothesis/intervention generation is an *agentic search* problem; MORPHEUS's A5 must argue geometry-based counterfactuals beat tournament search on efficiency/calibration.
**Novelty implication:** **Pre-empts** "AI proposes and ranks novel therapeutic hypotheses"; MORPHEUS repositions A5 as in-representation counterfactual *querying* vs. external evolutionary search.

### 21. SciAgents: Automating Scientific Discovery via Multi-Agent Graph Reasoning
Ghafarollahi, Buehler — Advanced Materials 2024. arXiv:2409.05556.
**Takeaway:** Couples ontological knowledge graphs with multi-agent LLMs to traverse a domain graph and generate/refine interdisciplinary hypotheses.
**Technical summary:** SciAgents samples paths through an ontological knowledge graph and assigns specialist agents to expand, critique, and ground hypotheses with retrieved data, applied to bio-inspired materials. The graph structure gives a "swarm intelligence" that finds non-obvious cross-domain links. Knowledge-graph grounding constrains agent hypotheses.
**Plain-English:** AI agents walk a big map of scientific concepts and connect distant ideas into new, testable hypotheses.
**Applicability:** A3 (structured NL<->knowledge grounding), A1. Design implication: the KG-as-scaffold idea supports treating MORPHEUS's pathway graph as an *addressable ontology* — slots indexed by a biological knowledge graph (A2).
**Novelty implication:** Reframes A2: pathway-addressability can be presented as ontology-indexed slots, a design SciAgents validates at the multi-agent level.

### 22. The AI Scientist: Towards Fully Automated Open-Ended Scientific Discovery
Lu, Lu, Lange, Foerster, Clune, Ha (Sakana) — 2024. arXiv:2408.06292.
**Takeaway:** End-to-end pipeline where an LLM ideates, codes/runs ML experiments, writes the paper, and an automated reviewer scores it — full research loop under ~$15/paper.
**Technical summary:** The system chains idea generation, experiment implementation via code execution, result visualization, paper writing, and an LLM reviewer approximating human review, demonstrated across ML subfields. It closes the loop from hypothesis to reviewed manuscript autonomously. Highlights automated *evaluation* as part of the discovery loop.
**Plain-English:** An AI that comes up with a research idea, runs the experiments, writes the paper, and reviews it — mostly by itself, cheaply.
**Applicability:** A1 (full task autonomy), A3 (automated review = capability evaluation). Design implication: its LLM-reviewer is a caution — automated evaluation can be gamed; MORPHEUS's A3 emergent-knowledge eval needs grounded, non-self-referential metrics.
**Novelty implication:** Neutral to MORPHEUS's representation claim; relevant as a warning about self-evaluation validity.

### 23. Agent Laboratory: Using LLM Agents as Research Assistants
Schmidgall, Su, Wang, et al. (AMD/Johns Hopkins) — 2025. arXiv:2501.04227.
**Takeaway:** A human-guided pipeline (literature review → experimentation → report writing) where specialized LLM agents assist rather than fully automate, with cost/quality analysis.
**Technical summary:** Agent Laboratory decomposes research into stages handled by role agents (PhD/postdoc/ML-engineer personas) using code execution and retrieval, and studies how human feedback at each stage improves outputs and lowers cost. It positions agents as accelerators of human research with measurable quality trade-offs. Emphasizes human-in-the-loop staging.
**Plain-English:** A pipeline of AI research assistants that do the literature review, run experiments, and draft the report, with a human steering at key steps.
**Applicability:** A1 (staged task decomposition), A3. Design implication: supports a human-steerable MORPHEUS front-end where NL prompts specify *stage-level* intent; reinforces that autonomy should be dialable, not maximal.
**Novelty implication:** Neutral; contributes autonomy-level framing for positioning MORPHEUS as a callable component.

### 24. CACTUS: Chemistry Agent Connecting Tool-Usage to Science
McNaughton, Ramalaxmi, Kruel, Knutson, Varikoti, Kumar — 2024. arXiv:2405.00972.
**Takeaway:** An open-model chemistry agent integrating cheminformatics tools for property prediction/similarity/drug-likeness, runnable on consumer hardware.
**Technical summary:** CACTUS wraps RDKit-style tools behind small open LLMs (Gemma/Mistral/Llama2/Falcon/MPT) via a ReAct loop and benchmarks them on chemistry QA, with Gemma-7b and Mistral-7b outperforming baselines. It shows capable tool-use agents need not be frontier-scale. Demonstrates local, reproducible molecular reasoning.
**Plain-English:** A lightweight chemistry AI that uses molecule software to answer questions and runs on an ordinary machine.
**Applicability:** A1, A4 (specialist tools as callables). Design implication: evidence that a modest controller suffices when tools are strong — supports MORPHEUS keeping the *trunk* strong and the routing thin.
**Novelty implication:** Minor pre-emption of NL→chemistry routing; strengthens "strong tools > big controller."

### 25. LAB-Bench: Measuring Capabilities of Language Models for Biology Research
Laurent, Janizek, Ruzo, Hinks, Hammerling, Narayanan, Ponnapati, White, Rodriques (FutureHouse) — 2024. arXiv:2407.10362.
**Takeaway:** 2,400+ MCQs across literature recall, figure interpretation, database navigation, and DNA/protein-sequence manipulation, with human-expert baselines.
**Technical summary:** LAB-Bench evaluates practical biology-research skills (cloning, protocol, sequence tasks) beyond textbook QA, benchmarking frontier LLMs vs. expert biologists. It exposes where agents are/aren't reliable enough to serve as research assistants. Provides task-grounded capability measurement for biology agents.
**Plain-English:** A biology "driving test" for AI, checking whether it can really navigate databases, read figures, and handle sequences like a working scientist.
**Applicability:** A3 (capability/emergence evaluation — central), A1. Design implication: MORPHEUS should report A3 emergent-knowledge claims against LAB-Bench-style task-grounded items, not free-form generation, to be credible.
**Novelty implication:** **Strengthens A3 evaluation rigor**; sets the bar MORPHEUS's emergent-knowledge claims must clear.

### 26. ScienceAgentBench: Rigorous Assessment of Language Agents for Data-Driven Discovery
Chen, Chen, Ning, et al., Ning, Gao, Su, Sun (OSU) — ICLR 2025. arXiv:2410.05080.
**Takeaway:** 102 expert-validated tasks from 44 papers across four sciences; best agent solves only 32.4% (o1-preview 42.2%) — current agents are far from reliable.
**Technical summary:** Each task requires generating a Python program to solve a real data-driven discovery problem, scored on output/program quality and cost. Across five models the ceiling is low, revealing brittleness in scientific code generation. A sobering, leakage-controlled benchmark.
**Plain-English:** A tough, realistic test of AI agents doing actual data-analysis science — and they mostly fail.
**Applicability:** A1 (task inference/execution eval), A3. Design implication: use as an external-validity check — if MORPHEUS is pitched as reducing agent brittleness on such tasks, this is the yardstick.
**Novelty implication:** Neutral; provides the confound-aware bar and argues headroom exists that a better representation could fill.

### 27. DiscoveryWorld: A Virtual Environment for Automated Scientific Discovery Agents
Jansen, Côté, Khot, Bransom, Dalvi Mishra, Majumder, Tafjord, Clark (AI2) — NeurIPS 2024. arXiv:2406.06769.
**Takeaway:** A simulated multi-task world (120 challenges, 8 topics incl. proteomics) requiring full hypothesize→experiment→analyze→act cycles; baseline agents struggle.
**Technical summary:** Tasks have difficulty tiers and parametric variants, with three automatic metrics (task completion, task-relevant actions, knowledge discovered) to score the *discovery process*, not just the answer. It isolates discovery reasoning from domain-specific tooling. Reveals current agents' weakness at genuine experimental reasoning.
**Plain-English:** A video-game-like lab world where AI must run experiments to figure out hidden rules — and today's agents mostly can't.
**Applicability:** A3 (measuring the discovery/emergence process), A5 (experiment-as-action). Design implication: its "knowledge discovered" metric is a model for evaluating whether MORPHEUS's counterfactual queries actually yield *new* knowledge vs. restating priors.
**Novelty implication:** Strengthens A3/A5 evaluation design; process-level metrics reduce the risk of overclaiming emergence.

### 28. ResearchAgent: Iterative Research Idea Generation over Scientific Literature
Baek, Jauhar, Cucerzan, Hwang (Microsoft/KAIST) — 2024. arXiv:2404.07738.
**Takeaway:** LLM + academic knowledge graph generate problems/methods/experiment designs, refined by multiple LLM ReviewingAgents aligned to human preferences.
**Technical summary:** ResearchAgent grounds ideation in an academic graph + entity knowledge store and iterates via review agents whose criteria are derived from human judgments. Cross-discipline evaluation shows more novel, valid proposals. Grounds idea generation in structured literature rather than free association.
**Plain-English:** An AI that reads connected research papers and brainstorms new, reviewed research ideas across fields.
**Applicability:** A3 (grounded hypothesis generation), A1. Design implication: the human-preference-aligned reviewer criteria are reusable for MORPHEUS's A3 evaluation of whether elicited biology is *novel and valid*, not just fluent.
**Novelty implication:** Neutral/strengthens A3 eval; another literature-grounded ideation precedent to distinguish from latent-representation elicitation.

### 29. ProtAgents: Protein Discovery via LLM Multi-Agent Collaborations
Ghafarollahi, Buehler (MIT) — Digital Discovery 2024. arXiv:2402.04268.
**Takeaway:** Multiple role-specialized LLM agents combine physics-based simulation, ML predictors, and retrieval to do de novo protein design.
**Technical summary:** Agents with distinct capabilities (structure analysis, physics simulation, knowledge retrieval) collaborate to design proteins with targeted mechanical properties, integrating first-principles simulation with ML in one autonomous loop. Shows multi-objective biomolecular design via agent division of labor. Physics tools ground the design choices.
**Plain-English:** A team of AI specialists — one simulates physics, one predicts, one looks things up — designs new proteins together.
**Applicability:** A4 (physics/ML predictors as callable specialists), A5 (design-as-intervention). Design implication: illustrates *when to call an external simulator* vs. encode — a concrete A4 encode-vs-RAG decision point (expensive physics → call it, don't encode).
**Novelty implication:** Supports A4 decision framing; pre-empts "agents design proteins" but is orthogonal to a unified-representation claim.

### 30. MLAgentBench: Evaluating Language Agents on Machine-Learning Experimentation
Huang, Vora, Liang, Leskovec (Stanford) — ICML 2024. arXiv:2310.03302.
**Takeaway:** 13 ML-experimentation tasks; best (Claude-Opus) agent averages 37.5% success with high variance and long-horizon/hallucination failures.
**Technical summary:** Agents (ReAct-based) must edit code, run experiments, and improve models; success ranges from strong on established datasets to 0% on recent Kaggle tasks. Exposes weaknesses in long-term planning and hallucination control. A rigorous experimentation benchmark.
**Plain-English:** A test of whether AI agents can actually run and improve machine-learning experiments — they manage only some of the time.
**Applicability:** A1 (experimentation-as-task-execution eval). Design implication: long-horizon planning fragility argues MORPHEUS should keep *task routing* short-horizon (one prompt → one addressed programme) rather than long agentic chains.
**Novelty implication:** Neutral; supports the "internalize routing, avoid long chains" design bet behind A1/A2.

### 31. Large Language Model Agents for Biological Intelligence across Genomics, Proteomics, Spatial Biology and Biomedicine
(Review) — Briefings in Bioinformatics 2026, 27(2):bbag110. Oxford Academic.
**Takeaway:** Survey codifying the five agent capabilities for biology — tool-use, multi-step planning, memory/self-reflection, knowledge grounding, multimodal integration.
**Technical summary:** Reviews LLM-agent systems across genomics/proteomics/spatial/biomedicine and abstracts the recurring capability stack, emphasizing ontology/database grounding and multimodal integration (sequences, structures, expression, images). Provides a capability taxonomy and open challenges. Consolidates the field's design patterns.
**Plain-English:** A map of how AI agents are being used across biology and the five core skills they all need.
**Applicability:** A1, A3, A4 (its "knowledge grounding" + "multimodal integration" axes map directly to MORPHEUS A3/A4). Design implication: adopt its five-capability checklist to position MORPHEUS as supplying the *multimodal-integration + grounding* substrate agents currently bolt on externally.
**Novelty implication:** **Reframes** MORPHEUS as the missing "encoded multimodal representation" layer under the agent stack — a cleaner complementarity story than competing on orchestration.

### 32. From AI for Science to Agentic Science: A Survey on Autonomous Scientific Discovery
(Survey) — 2025. arXiv:2508.14111.
**Takeaway:** Recent, broad survey charting the shift from tool-AI to autonomous agentic science, with taxonomy of autonomy levels and open problems.
**Technical summary:** Synthesizes agentic-science systems across domains, categorizing by autonomy, orchestration pattern (single vs. multi-agent), tool integration, and evaluation, and enumerates reliability/grounding gaps. Situates chemistry/biology exemplars (Coscientist, ChemCrow, Biomni, co-scientist) in one framework. A current field map.
**Plain-English:** An up-to-date overview of how AI is moving from "a tool scientists use" to "an agent that does science."
**Applicability:** A1, A3, A5 (frames the whole lane). Design implication: use its autonomy-level and orchestration taxonomy to explicitly place MORPHEUS (unified promptable trunk) *below* the orchestration layer, sharpening the "not-another-agent" positioning.
**Novelty implication:** **Reframes** the entire novelty argument: it documents that orchestration is crowded, implying MORPHEUS's defensible novelty is in the *representation/interface* primitive (A2/A4), not the agent loop.

---

### Lane synthesis for MORPHEUS

- **A1 is largely pre-empted at the system level.** Gorilla, ToolLLM, HuggingGPT, TxAgent, Biomni, CRISPR-GPT, BIA, AutoBA all do "NL request → infer task → call the right specialist predictor." MORPHEUS cannot claim NL task-routing as novel; it can claim *internalizing* routing in one frozen, promptable trunk with shared weights across programmes.
- **Biomni is the sharpest general prior-art collision** (one biomedical agent, many tasks, no per-task tuning). Steelman it explicitly; differentiate on identifiable/addressable slots (A2) and encoded-multimodal representation (A4) vs. code-orchestrated toolbox.
- **A5 has strong agentic precedents** (BioDiscoveryAgent, co-scientist, Coscientist, Virtual Lab) that already treat perturbations/experiments as proposable interventions — but via external search/experiment loops or LLM literature priors. MORPHEUS's wedge: counterfactual queries answered *inside* the representation's geometry, calibrated and one-shot, benchmarked against the LLM-only BioDiscoveryAgent baseline.
- **A3 evaluation is the most transferable win.** GeneAgent (DB self-verification), PaperQA2 (contradiction detection), LAB-Bench, ScienceAgentBench, DiscoveryWorld give ready-made, grounded, leakage-aware protocols. Reuse them so emergent-knowledge claims aren't attributable to LLM priors or scaffolding.
- **A4 gains a clean decision rule from the lane.** TxAgent/ProtAgents show verified databases and expensive physics are better *called as tools/RAG* than encoded; MORPHEUS's contribution is *learning the encode-vs-retrieve gate* (Toolformer-style utility filtering) rather than encoding everything.
