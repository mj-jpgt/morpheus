## AI for biological discovery / hypothesis generation

Lane l07_ai_discovery. Remit: model-driven biological discovery — hypothesis generation, AI-scientist systems, in-silico perturbation for discovery, LLM-driven biomedical hypothesis, mechanistic discovery from omics/imaging, and the epistemics of "real discovery vs restatement." Every entry maps to one or more MORPHEUS rebase axes:
- **A1** promptable unified representation + NL task auto-detection
- **A2** identified, pathway-addressable slots for reliable prompting (identifiability, per-programme addressability)
- **A3** NL<->biology grounding + emergent-knowledge elicitation AND its evaluation
- **A4** multimodal prompting: when to ENCODE a modality vs treat it as RAG context; frozen-trunk plug-in
- **A5** interventional/causal-geometry queries (counterfactual perturbation/drug as a query, not a retrained classifier)

---

### 1. Autonomous chemical research with large language models (Coscientist)
**Boiko, MacKnight, Kline, Gomes — Nature 623 (2023)** / preprint "Emergent autonomous scientific research capabilities of large language models," arXiv:2304.05332.
- **Takeaway:** An LLM orchestrator that plans, codes, and physically executes wet-lab chemistry (incl. Suzuki/Sonogashira cross-couplings) closes the ideation-to-execution loop autonomously.
- **Technical summary:** Coscientist combines a planner LLM with web-search, documentation-reading, and code-execution modules that drive robotic liquid handlers. It reasoned over reaction conditions, wrote control code, and optimized a real catalyzed cross-coupling from a natural-language goal, demonstrating tool-augmented autonomy rather than text-only suggestion.
- **Plain-English:** You tell an AI "run this kind of reaction," and it figures out the steps and actually operates the lab robot to do it.
- **Applicability:** A1 (NL goal -> task routing across heterogeneous tools), A5 (executes interventions, not just classifies). Design implication: MORPHEUS's "task auto-detection" should be validated end-to-end against an execution proxy, not only on held-out labels — the discovery value lives in the loop, not the classifier.
- **Novelty implication:** PRE-EMPTS any bare claim that "LLM proposes biology" is novel. MORPHEUS must differentiate on grounded, identifiable *representation-level* interventions, not agentic orchestration alone.

### 2. Towards an AI co-scientist (Google AI co-scientist)
**Gottweis, Weng, Daryin, Tu, Palepu, Vashisht, Natarajan et al. — arXiv:2502.18864 (2025); Nature (2026).**
- **Takeaway:** A multi-agent Gemini system that generates, debates (Elo tournament), and evolves hypotheses produced experimentally validated leads in liver fibrosis, drug repurposing (AML), and antimicrobial resistance.
- **Technical summary:** Specialized agents (Generation, Reflection/peer-review, Ranking via pairwise scientific-debate tournaments, Evolution, Meta-review) iterate on a research goal grounded in literature. Reported wet-lab validations include repurposed drugs inhibiting AML cell viability and a novel gene-transfer mechanism for cAMR later matched to an unpublished experimental finding.
- **Plain-English:** A committee of AI agents argues its way to research ideas, ranks them like a chess ladder, and the best ones held up in real experiments.
- **Applicability:** A1 (goal -> orchestrated task decomposition), A3 (elicits and ranks emergent hypotheses; the tournament *is* an evaluation of novelty), A5 (proposes interventions). Design implication: MORPHEUS should borrow the tournament-as-evaluation idea to *measure* emergent knowledge (A3) rather than assert it.
- **Novelty implication:** REFRAMES the bar for "discovery." Restatement-vs-discovery is now adjudicated by tournament + external validation; MORPHEUS's A3 evaluation must be at least this rigorous or it will read as weaker.

### 3. Robin: a multi-agent system for automating scientific discovery
**Ghareeb, Rodriques et al. (FutureHouse) — arXiv:2505.13400 (2025); Nature (2026), s41586-026-10652-y.**
- **Takeaway:** End-to-end multi-agent system that autonomously generated the hypothesis, chose assays, analyzed its own RNA-seq, and nominated ripasudil (a ROCK inhibitor) as a novel dry-AMD candidate.
- **Technical summary:** Crow/Falcon (shallow/deep literature agents) plus Finch (data-analysis agent over flow-cytometry and RNA-seq) coordinate a full loop. Robin proposed enhancing RPE-cell phagocytosis, then a circadian-modulation mechanism (KL001), and via follow-up RNA-seq surfaced ABCA1 (lipid efflux) as a candidate mechanism — humans only executed the bench steps.
- **Plain-English:** The AI ran the whole thinking part of a drug-discovery project — idea, experiment choice, data crunching — and found a plausible new blindness treatment.
- **Applicability:** A1 (autonomous task routing over modalities), A5 (drug candidate as the output of an interventional loop). Design implication: MORPHEUS's differentiator vs Robin is representation identifiability (A2) — Robin's mechanism claims come from LLM reasoning over external data, not from an addressable, causal latent.
- **Novelty implication:** PRE-EMPTS "AI proposes a validated drug candidate" as a MORPHEUS headline. MORPHEUS must claim the *how* (identified pathway slots + counterfactual geometry), not the *that*.

### 4. BioDiscoveryAgent: an AI agent for designing genetic perturbation experiments
**Roohani, Lee, Huang, Vora, Steinhart, K. Huang, Marson, Liang, Leskovec — arXiv:2405.17631 (2024/25).**
- **Takeaway:** An LLM agent (Claude 3.5 Sonnet) designs the next round of CRISPR perturbations in a closed loop, beating Bayesian optimization by 21% (46% on non-essential genes) across six datasets.
- **Technical summary:** The agent integrates literature search, dataset code-execution, and a second peer-review agent to iteratively select genes/gene-pairs to perturb, framing experiment design as sequential decision-making. It predicts productive combinations >2x better than random, addressing closed-loop combinatorial design without a trained surrogate.
- **Plain-English:** Instead of a statistical optimizer, an LLM reads the literature and the data to decide which genes to knock out next — and it picks better.
- **Applicability:** A5 (intervention selection as query), A1 (NL-driven experiment routing). Design implication: MORPHEUS could position its interventional queries as *complementary* to agentic design — the model supplies the causal-geometry prior the agent lacks.
- **Novelty implication:** STRENGTHENS the need for MORPHEUS's A5 causal geometry: agentic selection works but is data-hungry per loop; an identifiable interventional representation would cut loop count.

### 5. Empowering biomedical discovery with AI agents
**Gao, Bryson, ... Zitnik — Cell 187 (2024); arXiv:2404.02831.**
- **Takeaway:** A conceptual framework casting AI agents as "collaborators" combining LLM reasoning, ML tools, and experimental platforms across a spectrum of human-vs-AI autonomy.
- **Technical summary:** Defines agent roles (hypothesis generation, experiment planning, analysis), a taxonomy of autonomy levels, and the tool/skill/memory scaffolding needed for virtual-cell and therapeutic-design tasks. Argues discovery agents must integrate multi-scale biological knowledge and feedback from real experiments.
- **Plain-English:** A blueprint for how AI "lab partners" should be built to actually help discover biology, not just chat about it.
- **Applicability:** A1 (unified agent interface), A4 (multi-tool/multimodal integration). Design implication: use its autonomy taxonomy to place MORPHEUS — a *representation* that answers interventional queries sits below full agentic autonomy but supplies the missing grounded substrate.
- **Novelty implication:** REFRAMES MORPHEUS as infrastructure ("the virtual instrument") within an agent stack rather than a competing agent — sharpens the positioning.

### 6. Language agents achieve superhuman synthesis of scientific knowledge (PaperQA2)
**Skarlinski, Cox, Laurent, Braza, Hinks, Hammerling, Ponnapati, Rodriques, White — arXiv:2409.13740 (2024).**
- **Takeaway:** A retrieval agent matches/exceeds subject-matter experts on literature search, summarization, and contradiction detection; found ~2.34 contradictions per biology paper (70% expert-validated).
- **Technical summary:** PaperQA2 uses full-text retrieval, iterative evidence gathering, and citation-grounded generation, benchmarked on LitQA2. It produces Wikipedia-style summaries more accurate than human-written ones and detects factual contradictions across the corpus.
- **Plain-English:** An AI that reads papers better than experts for finding facts and spotting where papers disagree.
- **Applicability:** A4 (literature as high-quality RAG context — exactly the "treat-as-RAG" side of the encode-vs-RAG decision), A3 (grounding claims in citations). Design implication: MORPHEUS's A4 should adopt citation-grounded RAG for modalities that are cheap-to-retrieve and hard-to-encode (rare CNV/SNV annotations), reserving encoding for dense continuous modalities.
- **Novelty implication:** STRENGTHENS A4's encode-vs-RAG axis with evidence that RAG-grade literature synthesis is now superhuman — MORPHEUS shouldn't try to encode what RAG already does better.

### 7. The AI Scientist-v2: workshop-level automated scientific discovery via agentic tree search
**Yamada, Lange, Lu, Foerster, Ha, Clune et al. (Sakana AI) — arXiv:2504.08066 (2025).**
- **Takeaway:** Fully automated pipeline (idea -> experiment -> paper) whose output passed peer review at a workshop — the first AI-generated peer-reviewed paper.
- **Technical summary:** Replaces human-authored code templates with agentic tree search over experiment plans and a VLM-based reviewer feedback loop; generates hypotheses, runs ML experiments, and writes manuscripts autonomously. Highlights both capability and the fragility/hallucination of the current generation.
- **Plain-English:** An AI ran an entire ML research project start to finish, and one of its papers got accepted at a workshop.
- **Applicability:** A1 (task auto-detection across the research workflow), A3 (self-generated hypotheses and their evaluation). Design implication: its documented failure modes (fabricated results, shallow novelty) define exactly the "restatement" trap MORPHEUS's A3 evaluation must screen out.
- **Novelty implication:** REFRAMES "novelty" as adversarial — MORPHEUS needs an explicit restatement detector, since automated systems reliably produce fluent non-discoveries.

### 8. Large language models as biomedical hypothesis generators: a comprehensive evaluation
**Qi, Yang, Jiang et al. — arXiv:2407.08940 (2024); OpenReview.**
- **Takeaway:** Systematic eval showing LLMs generate novel, literature-consistent hypotheses even on background-knowledge unseen during training, with multi-agent/tool use improving uncertainty exploration.
- **Technical summary:** Builds a background-hypothesis dataset with a strict temporal cutoff to test genuine generalization vs memorization, evaluating zero-shot, few-shot, fine-tuning, and multi-agent settings. Finds real generative capacity but sensitivity to prompt/agent structure.
- **Plain-English:** A careful test of whether LLMs can invent (not just recall) biomedical hypotheses — they can, somewhat, if you test on truly new literature.
- **Applicability:** A3 (the temporal-cutoff design is a template for measuring emergent knowledge without leakage). Design implication: MORPHEUS's A3 emergent-knowledge claims must use a strict train/eval temporal split or they are unfalsifiable.
- **Novelty implication:** PROVIDES the evaluation methodology MORPHEUS needs for A3; failing to adopt a leakage-controlled protocol would be a novelty-claim liability.

### 9. Scientific hypothesis generation by LLMs: laboratory validation in breast cancer treatment
**Abdel-Rehim, Zenil, Orhobor, ..., Soldatova, King — J. R. Soc. Interface (2025); PMC12134935.**
- **Takeaway:** GPT-4 proposed synergistic drug pairs selectively killing MCF7 (cancer) over MCF10A (normal); multiple pairs (e.g., disulfiram+simvastatin, disulfiram+quinacrine) validated with above-control synergy in the lab.
- **Technical summary:** Two-round human-in-the-loop pipeline: GPT-4 nominates combinations, wet-lab viability/SynergyFinder assays test them, and results feed a second prompting round. 3/12 then 3/4 combinations validated, demonstrating iterative model-guided experimental design with real hit rates.
- **Plain-English:** GPT-4 guessed which drug pairs would team up against breast-cancer cells while sparing normal cells, and the lab confirmed several.
- **Applicability:** A5 (drug-combination as an interventional query), A3 (validated emergent hypotheses). Design implication: MORPHEUS's counterfactual drug queries should be reported with hit-rate-style metrics against wet-lab ground truth, matching this concreteness.
- **Novelty implication:** PRE-EMPTS "LLM finds validated drug synergies." MORPHEUS must claim mechanism-level, pathway-addressable prediction — not just ranked combinations.

### 10. Scientific hypothesis generation and validation: methods, datasets, and future directions (survey)
**Xiong et al. — arXiv:2505.04651 (2025).**
- **Takeaway:** Survey structuring the field into generation methods (LLM prompting, KG-augmented, multi-agent) and validation infrastructure (LabBench, AgentClinic, simulation), foregrounding the generation-vs-validation gap.
- **Technical summary:** Catalogs datasets and benchmarks and argues that hypothesis *validation* is the field's bottleneck, since fluent generation vastly outpaces reliable, reproducible confirmation. Maps where automated validation harnesses do and don't exist.
- **Plain-English:** A map of how AI makes scientific guesses and — the harder part — how we check whether they're actually right.
- **Applicability:** A3 (defines the validation landscape MORPHEUS must plug into). Design implication: pick an existing validation harness/benchmark rather than inventing one, so A3 claims are comparable.
- **Novelty implication:** REFRAMES the contribution surface — validation, not generation, is scarce; MORPHEUS's defensible novelty should sit on the validation/identifiability side.

### 11. Transfer learning enables predictions in network biology (Geneformer)
**Theodoris, Xiao, ..., Ellinor — Nature 618 (2023); s41586-023-06139-9.**
- **Takeaway:** A transformer pretrained on ~30M single-cell transcriptomes enables in-silico gene deletion to nominate therapeutic targets in data-limited settings.
- **Technical summary:** Rank-value gene encoding + masked pretraining yields context-aware embeddings; in-silico deletion perturbs the attention-derived network and reads out state shifts (e.g., moving diseased cardiomyocytes toward healthy), used to prioritize targets later validated experimentally.
- **Plain-English:** Train a language-model-like network on millions of cells, then "delete" genes on the computer to see which ones a disease cell needs — those become drug-target candidates.
- **Applicability:** A5 (in-silico deletion as interventional query), A2 (gene-level tokens as addressable units). Design implication: Geneformer's deletion is a discrete token-masking hack, not an identified causal latent — MORPHEUS's A2/A5 advantage is doing this on a *pathway-addressable, identifiable* representation.
- **Novelty implication:** PRE-EMPTS "foundation model does in-silico perturbation for discovery." MORPHEUS must differentiate on identifiability (A2) and true counterfactual geometry (A5), since Geneformer already owns the naive version.

### 12. scGPT: toward building a foundation model for single-cell multi-omics
**Cui, Wang, ..., Wang — Nature Methods 21 (2024); s41592-024-02201-0.**
- **Takeaway:** Generative transformer over 33M cells with a single backbone fine-tuned for annotation, integration, multi-omic fusion, perturbation prediction, and GRN inference.
- **Technical summary:** Uses gene/expression tokenization and specialized attention masking; one pretrained trunk is adapted per task via fine-tuning. Perturbation prediction and network inference are downstream heads rather than emergent zero-shot capabilities.
- **Plain-English:** A GPT for cells that, after light retraining, handles many single-cell tasks including predicting how cells respond to genetic changes.
- **Applicability:** A1 (unified trunk, many tasks), A2 (gene tokens), A4 (multi-omic integration). Design implication: scGPT still *fine-tunes per task* — MORPHEUS's A1 claim (prompt-time task auto-detection on a frozen trunk) is precisely the delta; make the frozen-vs-finetuned contrast explicit.
- **Novelty implication:** PRE-EMPTS "unified single-cell foundation model." MORPHEUS's novelty is *promptability without retraining* (A1) + *addressability* (A2), not the unified backbone per se.

### 13. Predicting transcriptional outcomes of novel multigene perturbations with GEARS
**Roohani, Huang, Leskovec — Nature Biotechnology 42 (2023); s41587-023-01905-6.**
- **Takeaway:** A GNN over a GO/knowledge graph predicts single- and combinatorial-perturbation transcriptional responses, including gene pairs never experimentally perturbed.
- **Technical summary:** Combines perturbation embeddings with gene embeddings on the GO graph via a compositional module to capture genetic interactions; achieves ~40% higher precision on interaction-subtype prediction and identifies strongest interactions 2x better than baselines.
- **Plain-English:** Using a graph of what genes do, it predicts what happens when you knock out gene combos you've never tested.
- **Applicability:** A5 (perturbation-response as query), A2 (knowledge-graph priors give pathway structure), A4 (GO graph as structured prior). Design implication: GEARS injects pathway structure via an *external* graph; MORPHEUS's A2 aims to make pathway structure *intrinsic and addressable* in the latent — contrast the two.
- **Novelty implication:** STRENGTHENS the case that pathway priors help; but PRE-EMPTS "predict unseen combos" — MORPHEUS must claim identifiability, not just combinatorial generalization.

### 14. Predicting cellular responses to complex perturbations in high-throughput screens (CPA)
**Lotfollahi, Klimovskaia Susmelj, ..., Theis — Molecular Systems Biology 19 (2023).**
- **Takeaway:** An adversarial autoencoder disentangles basal cell state from perturbation/dose/covariate embeddings, enabling compositional counterfactual prediction for unseen dose/cell-type/time.
- **Technical summary:** CPA factorizes the latent into a basal state and additively composed perturbation/covariate latents, adversarially removing perturbation info from the basal code. It predicts out-of-distribution responses by recombining learned embeddings — an explicitly compositional, interpretable interventional model.
- **Plain-English:** It separates "what kind of cell this is" from "what the drug did," so you can mix and match to predict responses to conditions you never measured.
- **Applicability:** A5 (counterfactual generation), A2 (disentangled, addressable perturbation slots). Design implication: CPA is the closest prior art to MORPHEUS's identifiable-slots idea — MORPHEUS must show it goes beyond additive compositionality (e.g., nonlinear pathway addressability, NL promptability A1) or it risks being "CPA + a prompt."
- **Novelty implication:** PRE-EMPTS naive A2. This is the sharpest novelty risk: MORPHEUS's identifiable, pathway-addressable slots overlap heavily with CPA's disentanglement — differentiate on grounding (A3) and promptability (A1).

### 15. Causal identification of single-cell perturbation effects with CINEMA-OT
**Dong, Wang, ..., Kluger — Nature Methods 20 (2023).**
- **Takeaway:** Optimal-transport matching separates confounding variation from true perturbation effects to build counterfactual cell pairs.
- **Technical summary:** CINEMA-OT decomposes variation into confounder and treatment factors, then uses OT to match treated cells to their counterfactual controls, yielding individualized treatment-effect estimates and synergy analysis at single-cell resolution.
- **Plain-English:** It pairs each treated cell with the untreated version it "would have been," isolating what the perturbation actually did.
- **Applicability:** A5 (counterfactual pairing is the core of interventional geometry). Design implication: MORPHEUS's counterfactual queries need an identifiability guarantee like CINEMA-OT's confounder separation, or its "causal" claims are correlational.
- **Novelty implication:** REFRAMES what "causal" must mean. If MORPHEUS calls its queries causal, it inherits CINEMA-OT's bar for confounder control — cite it as the standard.

### 16. Causal disentanglement for single-cell representations and controllable counterfactual generation (CausCell)
**(Nature Communications, 2025); s41467-025-62008-1.**
- **Takeaway:** A diffusion-based generative framework enforcing causal structure among latent concepts for controllable single-cell counterfactual generation.
- **Technical summary:** CausCell combines concept disentanglement with a structural causal prior and diffusion decoding, letting users intervene on individual (supervised or unsupervised) concepts and generate coherent counterfactual cells while preserving downstream structure.
- **Plain-English:** It learns cause-and-effect "knobs" inside a cell model, so you can turn one knob and generate a realistic cell showing the consequence.
- **Applicability:** A2 (concept-level addressable slots), A5 (controllable counterfactuals). Design implication: CausCell demonstrates concept-addressable intervention with generative fidelity — MORPHEUS's A2 must exceed it on *biological identifiability* (concepts = real pathways) and NL addressability (A1).
- **Novelty implication:** PRE-EMPTS "controllable counterfactual generation with disentangled concepts." MORPHEUS's edge must be pathway-grounded, promptable concepts, not disentanglement alone.

### 17. CRADLE-VAE: counterfactual reasoning-based artifact disentanglement for perturbation modeling
**arXiv:2409.05484 (2024).**
- **Takeaway:** A VAE that disentangles technical artifacts from true perturbation biology via counterfactual reasoning, improving perturbation-effect estimation quality.
- **Technical summary:** Learns separate latent subspaces for biological perturbation response and technical/quality artifacts, using counterfactual objectives so generated responses are artifact-free. Improves robustness of single-cell perturbation prediction on noisy screens.
- **Plain-English:** It tells apart "real biology" from "measurement junk" when modeling how cells react, giving cleaner predictions.
- **Applicability:** A2 (artifact vs signal slot separation), A5 (counterfactual objective). Design implication: MORPHEUS's identifiable slots should include an explicit technical-artifact channel, or batch effects will contaminate "pathway" slots.
- **Novelty implication:** STRENGTHENS the argument that identifiability requires modeling nuisance factors — a design requirement, not a novelty threat.

### 18. Scouter: predicting transcriptional responses to genetic perturbations with LLM embeddings
**(PMC12855003, 2025).**
- **Takeaway:** Uses LLM-derived gene/perturbation text embeddings to predict perturbation responses, bridging NL knowledge and single-cell prediction.
- **Technical summary:** Encodes genes via language-model embeddings of their descriptions and couples them to an expression predictor, letting textual prior knowledge inform quantitative perturbation-response forecasting, including for less-characterized genes.
- **Plain-English:** It reads what's known about a gene in plain text and uses that to predict how perturbing it changes the cell.
- **Applicability:** A3 (NL<->biology grounding operationalized), A4 (text embeddings as RAG-style context feeding a quantitative head). Design implication: this is a concrete A3/A4 pattern — inject NL knowledge as an *encoder input*, exactly MORPHEUS's grounding thesis; cite as feasibility evidence.
- **Novelty implication:** PRE-EMPTS "use LLM knowledge to ground perturbation prediction." MORPHEUS must go further: bidirectional grounding + emergent-knowledge *elicitation and measurement* (A3), not one-directional text->prediction.

### 19. How to build the virtual cell with artificial intelligence: priorities and opportunities
**Bunne, Roohani, ..., Leskovec, Theis, Regev et al. — Cell (2024); arXiv:2409.11654.**
- **Takeaway:** Community roadmap defining the AI Virtual Cell around universal representations, predictive perturbation simulation, interpretability, and "virtual instruments" for in-silico experiments.
- **Technical summary:** Frames three capability pillars — universal multi-scale representation, perturbation prediction, and virtual-instrument querying — and argues interpretability (knowing *why*) is as essential as predictive accuracy. Sets community standards to avoid mistaking artifacts for biology.
- **Plain-English:** The field's shared plan for an AI model of a cell you can run experiments on inside a computer.
- **Applicability:** A1 (universal representation), A4 (multimodal/multi-scale), A5 (virtual-instrument = interventional query). Design implication: MORPHEUS should explicitly claim which AIVC pillars it addresses; the "virtual instrument" language is the canonical framing for A5.
- **Novelty implication:** REFRAMES MORPHEUS within an accepted vocabulary — position A1/A2/A5 as concrete instantiations of AIVC pillars the roadmap left open (esp. identifiability, which the roadmap flags but does not solve).

### 20. Interpretable deep learning in single-cell omics
**Molho, Ding, ..., Tang — Bioinformatics 40 (2024), btae374; arXiv:2401.06823.**
- **Takeaway:** Review taxonomizing interpretability methods for single-cell DL and their use for mechanistic discovery (identity genes, gene sets, gene programs).
- **Technical summary:** Distinguishes post-hoc attribution from architecturally interpretable models (e.g., biologically-informed/gene-program layers) and connects each to discovery tasks — arguing black-box embeddings impede mechanistic claims.
- **Plain-English:** A survey of how to make cell-AI models explain themselves so their outputs count as biological insight, not just numbers.
- **Applicability:** A3 (interpretability as prerequisite for "real discovery"), A2 (program-layer architectures = addressable slots). Design implication: MORPHEUS should adopt biologically-informed (gene-program) layers so slots are inherently interpretable, supporting the discovery-vs-restatement distinction.
- **Novelty implication:** STRENGTHENS A2's rationale (interpretable-by-construction slots) and REFRAMES the restatement problem as an interpretability problem.

### 21. Assessing the limits of zero-shot foundation models in single-cell biology
**Kedzierska, Jiang, Eraslan, ..., Pinello — bioRxiv 2023.10.16.561085 (2023).**
- **Takeaway:** scGPT and Geneformer often fail to beat simple baselines (HVG + PCA/clustering) in zero-shot settings, questioning claimed emergent capability.
- **Technical summary:** Systematic zero-shot evaluation on clustering and batch integration shows foundation-model embeddings frequently underperform classical pipelines without task-specific fine-tuning, implying much "capability" is fine-tuning, not pretraining knowledge.
- **Plain-English:** The big cell "foundation models" don't actually beat basic methods until you retrain them — casting doubt on their out-of-the-box smarts.
- **Applicability:** A1 (directly challenges frozen-trunk promptability), A3 (emergent-knowledge claims need this baseline). Design implication: MORPHEUS's A1 promptability MUST be benchmarked against HVG+PCA baselines zero-shot, or reviewers will apply this critique.
- **Novelty implication:** PRE-EMPTS naive A1/A3 claims. This is a key novelty *risk*: MORPHEUS must show its frozen-trunk promptability beats trivial baselines, which prior FMs did not.

### 22. SciKnowEval: evaluating multi-level scientific knowledge of large language models
**Feng, Wei, ..., Chen et al. — arXiv:2406.09098 (2024).**
- **Takeaway:** Benchmark grading LLMs across five cognitive levels (memory, comprehension, reasoning, discernment, application) in biology/chemistry/physics/materials.
- **Technical summary:** Assembles a large multi-level dataset (incl. biology/chemistry) to separate rote recall from genuine reasoning/application; reveals models strong on memory but weaker on discernment/application — the levels most relevant to discovery.
- **Plain-English:** A test that distinguishes an AI that memorized science facts from one that can actually reason and apply them.
- **Applicability:** A3 (multi-level framing operationalizes "emergent knowledge" vs "restatement"). Design implication: MORPHEUS's A3 evaluation should report the cognitive *level* of the knowledge elicited, not a single accuracy number.
- **Novelty implication:** PROVIDES vocabulary for A3 — MORPHEUS can claim it elicits application-level (not memory-level) biological knowledge, a sharper, testable novelty claim.

### 23. Evaluating large language models in scientific discovery (SDE framework)
**arXiv:2512.15567 (2025/26).**
- **Takeaway:** Argues QA benchmarks are poor proxies for discovery and proposes a Scientific Discovery Evaluation tying questions to modular research scenarios.
- **Technical summary:** Constructs evaluation where each item maps to a step in a realistic discovery workflow (hypothesis, design, interpretation), measuring contribution to discovery rather than answer-lookup accuracy.
- **Plain-English:** A better way to grade AI on whether it helps *discover* things, not whether it aces a science quiz.
- **Applicability:** A3 (the core "discovery vs restatement" measurement problem). Design implication: adopt a scenario-anchored evaluation for MORPHEUS's discovery claims so a hypothesis is scored by its role in a workflow, not by literature match.
- **Novelty implication:** REFRAMES the entire lane's evaluation. MORPHEUS's strongest defensible novelty may be *how it is evaluated* (scenario-anchored discovery) as much as what it predicts.

### 24. Knowledge graphs for drug repurposing: a review of databases and methods
**Gu, ..., (Briefings in Bioinformatics 25, 2024); bbae461.**
- **Takeaway:** Surveys KG-based repurposing (embedding, path-based, GNN) — the dominant non-LLM route to interventional hypotheses (drug->disease).
- **Technical summary:** Catalogs biomedical KGs and completion methods (TransE/RotatE/ComplEx, GNNs, path reasoning) for predicting drug-disease links, emphasizing explainability via reasoning paths.
- **Plain-English:** A tour of how connecting biomedical facts into a graph lets you predict which existing drugs might treat which diseases.
- **Applicability:** A4 (structured KG as RAG-style context vs encoding), A5 (drug-repurposing as interventional query). Design implication: for CNV/SNV/annotation modalities, KG-as-context (A4 RAG side) is a proven pattern; MORPHEUS should route sparse symbolic knowledge through a KG rather than encode it.
- **Novelty implication:** PRE-EMPTS the KG-completion route to interventional hypotheses. MORPHEUS's A5 must offer something KG completion cannot — quantitative, cell-state-specific counterfactual *geometry*, not link prediction.

### 25. Multi-agent LLMs for biomedical hypothesis generation in drug-combination discovery
**(iScience, 2025); S2589-0042(25)02245-X.**
- **Takeaway:** A collaboration-inspired multi-agent LLM framework (Researcher/Reviewers/Moderator) generated and validated an Alzheimer's drug combination that reduced amyloid aggregation in vitro (external validation accuracy 0.82).
- **Technical summary:** Specialized agents iterate hypotheses via in-context learning and mutual critique, outperforming knowledge-based baselines; a nominated combination was experimentally shown to reduce amyloid aggregation, closing generation-to-validation.
- **Plain-English:** A team of AI agents debated their way to a drug combo for Alzheimer's that actually cut amyloid clumping in a dish.
- **Applicability:** A3 (multi-agent critique as emergent-knowledge elicitation + eval), A5 (drug combination as intervention). Design implication: multi-agent critique is an alternative to representation-level identifiability for reliability — MORPHEUS should argue why an identified representation (A2) is more sample-efficient than agent debate.
- **Novelty implication:** PRE-EMPTS "multi-agent LLM -> validated drug combination." MORPHEUS differentiates by making the *representation* the reasoner, not an ensemble of prompted agents.

### 26. Augmenting large language models with chemistry tools (ChemCrow)
**M. Bran, Cox, Schilter, Baldassari, White, Schwaller — Nature Machine Intelligence 6 (2024); arXiv:2304.05376.**
- **Takeaway:** GPT-4 + 18 expert tools autonomously plans syntheses and drives discovery tasks, showing tool-grounding sharply reduces hallucination vs bare LLM.
- **Technical summary:** A ReAct-style agent selects among cheminformatics/RxN/lab tools to plan and (with hardware) execute synthesis of small molecules and materials; expert evaluation shows tool-augmented outputs are far more reliable than the base model's.
- **Plain-English:** Give an LLM real chemistry software as tools and it can plan (and run) actual molecule synthesis instead of making things up.
- **Applicability:** A1 (NL goal -> tool routing), A4 (tools/knowledge as external context vs internal encoding). Design implication: reinforces MORPHEUS's encode-vs-RAG axis — grounding via external tools/knowledge beats forcing everything into weights.
- **Novelty implication:** STRENGTHENS A4; the "tool-grounding cuts hallucination" result supports MORPHEUS's choice to RAG some modalities rather than encode them.

### 27. PyTDC: a multimodal ML platform for biomedical foundation models
**Velez-Arce, ..., Zitnik — arXiv:2505.05577 (2025).**
- **Takeaway:** Training/evaluation/inference platform standardizing multimodal biomedical FM benchmarking, including single-cell and perturbation tasks.
- **Technical summary:** Extends Therapeutics Data Commons with streaming multimodal datasets and model-agnostic APIs for cell-type/perturbation/therapeutic tasks, enabling reproducible cross-model comparison.
- **Plain-English:** A standardized toolbox and leaderboard so different biomedical AI models can be compared fairly on the same tasks.
- **Applicability:** A4 (multimodal task harness), A3 (reproducible eval). Design implication: MORPHEUS should report on a shared harness like PyTDC/TDC so its multimodal and interventional claims are comparable, not bespoke.
- **Novelty implication:** Neutral/STRENGTHENS — provides the comparison substrate; using it makes MORPHEUS's novelty claims credible rather than self-graded.

### 28. Predicting cellular responses to perturbation across diverse contexts with STATE
**Arc Institute (Adduri, Kim, Goodarzi, Hsu et al.) — bioRxiv (2025).**
- **Takeaway:** A large perturbation-prediction model separating a cell-set "State Transition" model from a cell-"State Embedding," trained across many contexts to predict responses to unseen perturbations/contexts.
- **Technical summary:** STATE models perturbation as a transition over sets of cells (capturing distributional, not just mean, shifts) atop a self-supervised cell embedding, trained on large observational + perturbational atlases; targets cross-context generalization where prior FMs underperformed.
- **Plain-English:** A big model that predicts how whole populations of cells shift when you perturb them, even in cell types it wasn't trained on.
- **Applicability:** A5 (perturbation-response across contexts = interventional query at scale), A1 (shared trunk across contexts). Design implication: STATE's set-level distributional transition is a strong A5 baseline — MORPHEUS's counterfactual geometry must add identifiability/addressability (A2) that STATE's embeddings leave implicit.
- **Novelty implication:** PRE-EMPTS "scale + context-general perturbation prediction." MORPHEUS's novelty must be the *addressable/identifiable* structure of the transition, not scale or context-generalization.

### 29. Self-driven biological discovery through automated hypothesis generation and experimental validation
**bioRxiv 2025.06.24.661378 (2025).**
- **Takeaway:** A closed-loop system that autonomously generates biological hypotheses and drives their experimental validation, reported as a self-driving discovery pipeline.
- **Technical summary:** Couples LLM-based hypothesis generation with automated experiment selection and analysis in an iterative loop, aiming to demonstrate genuine (not restated) discovery via prospective wet-lab confirmation.
- **Plain-English:** An AI system that comes up with biological ideas and then runs the experiments to check them, largely on its own.
- **Applicability:** A1 (autonomous task routing), A3 (validated emergent hypotheses), A5 (interventions in the loop). Design implication: another end-to-end competitor — reinforces that MORPHEUS should not compete on loop autonomy but on representation identifiability feeding such loops.
- **Novelty implication:** PRE-EMPTS end-to-end autonomous discovery framing; MORPHEUS is a *component* (the grounded interventional representation) inside these loops, not a rival loop.

---

### Cross-cutting synthesis for MORPHEUS

- **Biggest novelty risk (A2/A5):** CPA (#14), CausCell (#16), and GEARS (#13) already deliver disentangled/addressable/compositional interventional prediction. MORPHEUS's identifiable, pathway-addressable slots overlap heavily; the defensible delta must be (a) *biological* identifiability (slots = real pathways, validated), (b) NL promptability (A1) on a frozen trunk, and (c) grounding/elicitation with a measurable emergent-knowledge protocol (A3).
- **A1 is the sharpest exposed flank:** Kedzierska (#21) shows FMs fail zero-shot vs trivial baselines; scGPT (#12) still fine-tunes per task. MORPHEUS's "prompt-time task auto-detection on frozen trunk" is genuinely novel *only if* benchmarked to beat HVG+PCA and non-fine-tuned baselines.
- **A3 needs an evaluation, not an assertion:** Qi et al (#8, temporal cutoff), SciKnowEval (#22, cognitive levels), SDE (#23, scenario-anchored) collectively define how to measure emergent knowledge; MORPHEUS should commit to one to avoid the "restatement" critique that AI-Scientist-v2 (#7) exemplifies.
- **A4 encode-vs-RAG has strong support:** PaperQA2 (#6), ChemCrow (#26), KG-repurposing (#24), Scouter (#18) all show symbolic/textual knowledge is best RAG'd, while dense continuous modalities are best encoded — MORPHEUS's routing rule is well-founded.
- **Discovery-vs-restatement is the epistemic core:** the AI-scientist systems (#1-3, #5, #7, #29) prove agentic autonomy is commoditizing; MORPHEUS's durable claim is *representation-level* grounded interventional reasoning that these agents currently lack.
