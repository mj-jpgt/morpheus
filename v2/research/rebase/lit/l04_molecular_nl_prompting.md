## Molecular prompting & NL-prompting of scientific models

*Lane l04. Remit: predicting molecular/pathway state from images or other modalities via "prompting"; instruction-tuned scientific LLMs; natural-language interfaces to biological predictors; text-conditioned biomedical prediction. Every entry maps to MORPHEUS rebase axes A1-A5. NOT covered here: unified task-routing interfaces (l05), agentic tool-use (l13), pure fusion architectures (l01), omics FMs per se (l03).*

Axis legend: **A1** promptable rep + NL task auto-detection · **A2** identified/pathway-addressable slots · **A3** NL⇄biology grounding + emergent-knowledge elicitation & its evaluation · **A4** multimodal prompting: encode-vs-RAG, frozen-trunk plug-in · **A5** interventional/causal queries.

---

### Group I — Molecule ⇄ natural-language models (the founding "molecular prompting" line)

**1. Text2Mol: Cross-Modal Molecule Retrieval with Natural Language Queries** (Edwards, Zhai, Ji; EMNLP 2021). https://aclanthology.org/2021.emnlp-main.47/
- *Takeaway:* the first system to retrieve molecules from free-text descriptions by learning a shared molecule–language embedding.
- *Technical summary:* Learns an aligned semantic embedding space between SMILES/graph molecule encoders and a text encoder (SciBERT), trained with a contrastive/paired objective on ~33k molecule-description pairs from ChEBI, plus a cross-modal attention reranker for explainability. Establishes the retrieval framing ("find the molecule the sentence is about") that all later molecule-text work builds on.
- *Plain-English:* You type a sentence describing what a chemical does and it finds the matching molecule, by putting words and molecules in the same "meaning space."
- *Applicability:* **A3** (foundational NL⇄biology grounding via shared embedding) and **A1** (text-as-query over a molecular representation). Design implication for MORPHEUS: the shared-space + reranker pattern is a minimal, evaluable prototype for the TQI's "task-text → query" conditioning, and its retrieval metric (Hits@k, MRR) is a template for measuring grounding quality without downstream accuracy.
- *Novelty implication:* Pre-empts any claim that "text-as-a-query-over-a-biological-representation" is itself novel — MORPHEUS must locate novelty in *task inference/routing over a unified multimodal tumor state*, not in text-conditioned retrieval.

**2. MolT5 — Translation between Molecules and Natural Language** (Edwards et al.; EMNLP 2022). https://arxiv.org/abs/2204.11817
- *Takeaway:* casts molecule captioning and text-to-molecule generation as a single T5-style translation problem.
- *Technical summary:* Self-supervised denoising pretraining on unpaired SMILES + natural-language text, then fine-tuning for (a) molecule→caption and (b) description→SMILES generation; released small/base/large checkpoints. Overcomes chemistry data scarcity by single-modal pretraining before cross-modal fine-tuning.
- *Plain-English:* Treats "describe this molecule" and "draw the molecule from this description" like translating between two languages.
- *Applicability:* **A1/A3.** Implication: the encoder-decoder "translate the biology into language and back" formulation is a candidate for MORPHEUS's NL output head (A3), but MolT5 is generative-only with no task inference — reinforces that generation ≠ the hard part.
- *Novelty implication:* Reframes MORPHEUS's A3 claim: emergent-knowledge *elicitation and its measurement* is the delta, since captioning/generation is a solved, commoditized capability.

**3. MoMu — A Molecular Multimodal Foundation Model Associating Molecule Graphs with Natural Language** (Su et al.; 2022). https://arxiv.org/abs/2209.05481
- *Takeaway:* contrastively ties a molecule *graph* encoder to a text encoder for zero-shot cross-modal tasks.
- *Technical summary:* Weakly-supervised contrastive pretraining pairing molecular graphs (GIN) with sentences mined from scientific literature; enables zero-shot text-to-graph retrieval, molecule captioning, and text-guided generation. Extends the Text2Mol idea from fingerprints/SMILES to graph structure.
- *Plain-English:* Links the 2-D structure diagram of a molecule to the words scientists use about it, learned from paper text.
- *Applicability:* **A3/A4.** Implication: shows literature text as a cheap supervision signal for grounding — relevant to MORPHEUS's closed-RAG hypothesis-cards and to A4's "when is text best as context."
- *Novelty implication:* Strengthens the case that contrastive text-grounding of a structured biological object is well-trodden; MORPHEUS's structured object (tumor state) must add identifiability (A2) to differentiate.

**4. MoleculeSTM — Multi-modal Molecule Structure–Text Model for Text-based Retrieval and Editing** (Liu et al.; Nature Machine Intelligence 2023). https://arxiv.org/abs/2212.10789
- *Takeaway:* open-vocabulary molecule-text model whose headline feat is *editing a molecule toward a natural-language goal*.
- *Technical summary:* Contrastive pretraining on 280k structure-text pairs (PubChemSTM); demonstrates zero-shot structure-text retrieval and, crucially, *text-guided molecule editing* ("make this more soluble") via latent optimization against the text encoder. Reports SOTA generalization to novel biochemical concepts.
- *Plain-English:* You describe a property in words and it nudges a molecule's structure to acquire that property.
- *Applicability:* **A5** (edit-toward-a-described-goal is an *interventional query expressed in NL*) and **A3.** Implication: MoleculeSTM is the closest existing analogue to MORPHEUS's "drug/perturbation as a query" — but the intervention lives in *chemical* space, not a *causal tumor-state* space.
- *Novelty implication:* **Candidate novelty risk.** Pre-empts a naïve "NL intervention as a query is new" claim. MORPHEUS's A5 novelty must be the *causal-geometry* grounding (geodesic ≈ causal distance) applied *after encoding* on a frozen tumor state, not latent-space editing of the input object.

**5. BioT5 — Enriching Cross-modal Integration in Biology with Chemical Knowledge and Natural Language Associations** (Pei et al.; EMNLP 2023). https://arxiv.org/abs/2310.07276
- *Takeaway:* unifies molecules (SELFIES), proteins (FASTA), and literature text in one T5 pretraining with robust molecular tokenization.
- *Technical summary:* Uses SELFIES for 100%-valid molecule strings and mines bio-entity context from unstructured literature; trained across molecule/protein/text with structured wrapping tokens, improving molecule property prediction, captioning, and protein QA.
- *Plain-English:* One model that reads chemistry, protein sequences, and English about them together, with a molecule encoding that never produces "invalid" molecules.
- *Applicability:* **A1/A3.** Implication: the SELFIES-style *always-valid* token choice is a lesson for MORPHEUS's slot vocabulary — a promptable representation needs a token/slot format that cannot express ill-formed states.
- *Novelty implication:* Strengthens A1 prior art volume; a multi-biomodal text model is not itself novel.

**6. BioT5+ — Towards Generalized Biological Understanding with IUPAC Integration and Multi-task Tuning** (Pei et al.; 2024). https://arxiv.org/abs/2402.17810
- *Takeaway:* adds IUPAC-name grounding and broad multi-task instruction tuning on top of BioT5.
- *Technical summary:* Integrates IUPAC nomenclature (the names appearing in text) with molecular structure and scales multi-task tuning across dozens of molecule/protein/text tasks, closing gaps to specialist models on captioning and property tasks.
- *Plain-English:* Teaches the model the *human names* of molecules so its text and structure knowledge line up better.
- *Applicability:* **A1/A3.** Implication: name↔structure alignment is a grounding trick; for MORPHEUS, pathway *names* are the natural "IUPAC" that should index addressable slots (A2).
- *Novelty implication:* Neutral; reinforces that multi-task instruction tuning is standard.

**7. MolCA — Molecular Graph-Language Modeling with Cross-Modal Projector and Uni-Modal Adapter** (Liu et al.; EMNLP 2023). https://arxiv.org/abs/2310.12798
- *Takeaway:* a Q-Former projector bridges a frozen graph encoder into an LLM's text space, with LoRA adapters for tasks.
- *Technical summary:* Connects a 2-D graph encoder to a language model via a Q-Former cross-modal projector, plus a uni-modal LoRA adapter for downstream tuning; strong on captioning and molecule-text retrieval/QA. Architecturally mirrors BLIP-2's frozen-encoder + Q-Former recipe.
- *Plain-English:* Plugs a molecule "reader" into a chatbot using a small translator module, keeping the reader frozen.
- *Applicability:* **A4** (frozen-trunk plug-in via a projector is exactly MORPHEUS's frozen-representation ambition) and **A2** (query-tokens as addressable slots). Implication: the Q-Former's learned query tokens are a concrete template for the `(batch, n_pathways, D)` slots the TQI needs but the current encoder never exposes.
- *Novelty implication:* Pre-empts "projector-into-LLM over a frozen biological encoder" as novel — MORPHEUS should cite this as the mechanism and locate novelty in *what the slots identify* (pathways) and *how they are queried* (causal).

**8. GIT-Mol — A Multi-modal Large Language Model for Molecular Science with Graph, Image, and Text** (P. Liu et al.; Computers in Biology and Medicine 2024). https://arxiv.org/abs/2308.06911
- *Takeaway:* aligns molecule graph + structure image + text into a unified latent via a "GIT-Former."
- *Technical summary:* Introduces GIT-Former to map three modalities (graph, 2-D image, text) into one latent space with any-to-language generation; improves property prediction and captioning over single-modality baselines.
- *Plain-English:* One model that looks at a molecule as a graph, as a picture, and as words, all at once.
- *Applicability:* **A4** (multi-encoder alignment; when to add an image modality vs text). Implication: informs MORPHEUS's encode-vs-RAG study — GIT-Mol *encodes* every modality; MORPHEUS's open question is when that beats retrieval.
- *Novelty implication:* Neutral prior art for A4's "encode everything" pole.

**9. Unifying Molecular and Textual Representations via Multi-task Language Modelling (Text+Chem T5)** (Christofidellis et al.; ICML 2023). https://arxiv.org/abs/2301.12586
- *Takeaway:* a *single* multi-domain LM handles chemical and natural-language tasks with shared weights and no task-specific heads.
- *Technical summary:* Multi-task language modeling over natural + chemical language (forward/retro synthesis, captioning, text→molecule) with one set of weights; performs single- and cross-domain tasks without per-task fine-tuning, an early demonstration of *implicit task handling* from the prompt.
- *Plain-English:* One model that, depending on what you ask, does chemistry translation, description, or synthesis — no switching models.
- *Applicability:* **A1** (closest molecular-domain analogue to "infer the task from the prompt"). Implication: directly relevant to MORPHEUS A1's task auto-detection — but note the task is signaled by an explicit prefix, not *inferred*.
- *Novelty implication:* **Candidate novelty risk for A1.** A one-model-many-molecular-tasks system already exists; MORPHEUS's A1 delta must be genuine *auto-detection/routing with abstention* over a *multimodal tumor state*, not prompt-prefixed task selection.

**10. nach0 — Multimodal Natural and Chemical Languages Foundation Model** (Livne et al.; Chemical Science 2024). https://arxiv.org/abs/2311.12410
- *Takeaway:* an instruction-tuned encoder-decoder spanning literature text, patents, and molecule strings across QA/NER/generation.
- *Technical summary:* Pretrained on unlabeled scientific text, patents, and SMILES, then instruction-tuned; solves biomedical QA, entity recognition, and molecular generation in a unified seq2seq format.
- *Plain-English:* A science-and-chemistry model you instruct in plain language to answer questions or invent molecules.
- *Applicability:* **A1/A3.** Implication: another instance of the instruction-tuned scientific-LLM template MORPHEUS's NL interface will be compared against.
- *Novelty implication:* Adds to A1 prior-art density.

**11. 3D-MolT5 — Leveraging Discrete Structural Information for Molecule-Text Modeling** (Pei et al.; ICLR 2025). https://arxiv.org/abs/2406.05797
- *Takeaway:* injects tokenized 3-D geometry into a T5 molecule-text model.
- *Technical summary:* Discretizes 3-D structural information into tokens fused with 1-D sequence and text, improving property prediction and 3-D-dependent tasks over 1-D-only molecule-text models.
- *Plain-English:* Adds the molecule's 3-D shape as extra "words" so the model reasons about geometry, not just formula.
- *Applicability:* **A4** (adding a geometric modality by *encoding* it as tokens). Implication: supports MORPHEUS's "encode as slots" option for structured modalities where geometry matters.
- *Novelty implication:* Neutral.

**12. ChemLML — Chemical Language Model Linker: blending text and molecules with modular adapters** (Zhang et al.; 2024). https://arxiv.org/abs/2410.20182
- *Takeaway:* a modular *adapter* that links a frozen text encoder to a frozen molecule generator — mix-and-match trunks.
- *Technical summary:* Trains only a lightweight linker between pretrained text and molecule models to enable text-conditioned molecule generation, avoiding full retraining; systematically compares different frozen text/molecule backbones.
- *Plain-English:* A small connector that lets an off-the-shelf language model drive an off-the-shelf molecule generator.
- *Applicability:* **A4** (frozen-trunk plug-in; the cleanest "adapter, not retrain" evidence). Implication: strong precedent for MORPHEUS's frozen-representation + light adapter design and for the *marquee proteomics plug-in test*.
- *Novelty implication:* Pre-empts "adapters over frozen trunks is novel"; MORPHEUS's novelty is what the frozen trunk *is* (identified tumor state), not the adapter pattern.

**13. Mol-Instructions — A Large-Scale Biomolecular Instruction Dataset for Large Language Models** (Fang et al.; ICLR 2024). https://arxiv.org/abs/2306.08018
- *Takeaway:* the standard 2M+ instruction dataset covering molecule, protein, and biomolecular-text tasks.
- *Technical summary:* 148k molecule-oriented, 505k protein-oriented, and 53k biomolecular-text instructions across 17 subtasks and 11+ properties; used to instruction-tune general LLMs into biomolecular assistants.
- *Plain-English:* A giant question-and-answer workbook that teaches general chatbots to handle molecules and proteins.
- *Applicability:* **A1/A3.** Implication: template for how MORPHEUS would build an *instruction/prompt corpus* for tumor-state tasks; also a warning that instruction data alone yields text competence, not identified representations.
- *Novelty implication:* Establishes instruction-tuning-for-biomolecules as commodity; MORPHEUS's evaluation (A3) must go beyond such instruction-following accuracy.

**14. InstructMol — Multi-Modal Integration for a Versatile and Reliable Molecular Assistant in Drug Discovery** (Cao et al.; COLING 2025). https://arxiv.org/abs/2311.16208
- *Takeaway:* two-stage instruction tuning aligns molecule graphs to an LLM, yielding *plugin* checkpoints.
- *Technical summary:* Stage 1 aligns graph-text caption pairs; stage 2 task-specific tuning for property prediction, captioning, reagent/reaction tasks; produces lightweight LoRA "plugin" checkpoints that can be loaded/combined while retaining general dialogue.
- *Plain-English:* A molecule chat-assistant you extend by snapping in small skill "plugins."
- *Applicability:* **A1/A4.** Implication: the composable-plugin design supports MORPHEUS's per-task adapters over a frozen trunk (A4) and NL task surface (A1).
- *Novelty implication:* Adds to A1/A4 prior art; the plugin-over-frozen-LLM idea is established.

**15. Property-Enhanced Instruction Tuning for Multi-Task Molecule Generation with LLMs** (Xu et al.; 2024). https://arxiv.org/abs/2412.18084
- *Takeaway:* conditions molecule generation on explicit numeric property targets inside instructions.
- *Technical summary:* Injects computed molecular properties into instruction tuning so an LLM generates molecules meeting quantitative property constraints across multiple tasks, improving controllability over caption-only conditioning.
- *Plain-English:* Tell the model the exact numbers you want (e.g., a solubility value) and it makes molecules hitting them.
- *Applicability:* **A5** (numeric-target conditioning ≈ a constrained interventional query) and **A4** (numbers-as-context vs numbers-encoded). Implication: relevant to MORPHEUS's A4 question of *encoding numeric modalities (phospho/CNV) vs passing them as prompt context*.
- *Novelty implication:* Informs A4/A5 framing; shows numeric conditioning of a generator is feasible but not causal.

**16. ChatGPT-powered Conversational Drug Editing Using Retrieval and Domain Feedback (ChatDrug)** (S. Liu et al.; ICLR 2024). https://arxiv.org/abs/2305.18090
- *Takeaway:* iterative NL drug editing that *retrieves* exemplars and *feeds back* domain checks around a frozen LLM.
- *Technical summary:* A prompt module + retrieval-and-domain-feedback loop + conversation module edits small molecules, peptides, and proteins toward described properties, with retrieved neighbors as in-context guidance and domain tools as validators.
- *Plain-English:* A back-and-forth chat that redesigns a drug toward what you ask, looking up similar known drugs to help.
- *Applicability:* **A4** (RAG-as-context around a frozen model — a direct encode-vs-retrieve data point) and **A5** (editing-as-query). Implication: strongest lane evidence for MORPHEUS's A4 "retrieve rather than encode" pole and for RAG-augmented interventional queries.
- *Novelty implication:* **Candidate novelty risk for A4/A5** — a retrieval-augmented, feedback-guided NL editing loop already exists; MORPHEUS must differentiate on *causal geometry* and *identified slots*, not the RAG loop itself.

**17. Emerging Opportunities of Using LLMs for Translation Between Drug Molecules and Indications** (Meyers et al.; 2024). https://arxiv.org/abs/2402.09588
- *Takeaway:* frames drug↔indication mapping as a bidirectional translation task for LLMs.
- *Technical summary:* Benchmarks LLMs on generating a molecule from a therapeutic indication and vice versa, arguing this clinically-grounded translation is a distinct, underexplored capability from captioning.
- *Plain-English:* Ask "what molecule treats this disease?" (or the reverse) and have a language model answer.
- *Applicability:* **A3** (clinical-language ⇄ molecule grounding closest to MORPHEUS's clinical-question framing). Implication: an indication-conditioned query is a clinical-flavored precedent for the TQI.
- *Novelty implication:* Reframes A3 toward clinically-actionable grounding, where MORPHEUS's multimodal tumor context could add value beyond a text-only LLM.

---

### Group II — Cell / gene expression ⇄ natural language (prompting omics readouts)

**18. GenePT — Simple and effective embedding model for single-cell biology built from ChatGPT** (Chen & Zou; Nature Biomedical Engineering 2024; preprint bioRxiv 2023.10.16.562533). https://www.nature.com/articles/s41551-024-01284-6
- *Takeaway:* LLM text-embeddings of gene descriptions rival expression-pretrained single-cell FMs.
- *Technical summary:* Embeds NCBI gene summaries with GPT-3.5 to get per-gene vectors; builds cell embeddings by expression-weighted averaging or by embedding a rank-ordered "cell sentence." Matches/beats Geneformer/scGPT on gene-property and cell-type tasks with no omics pretraining.
- *Plain-English:* Just using ChatGPT's understanding of what each gene does — read from the literature — is enough to represent cells surprisingly well.
- *Applicability:* **A3** (literature-derived, *addressable per-gene* text knowledge) and **A2** (gene-indexed slots). Implication: strong evidence that MORPHEUS's pathway slots could be *initialized/grounded from LLM literature embeddings*, and that emergent biological knowledge is measurable via such probes.
- *Novelty implication:* **Candidate novelty risk for A3/A2** — LLM literature knowledge already encodes substantial "biology"; MORPHEUS must show its multimodal state adds knowledge *not present in text priors* (a decisive emergence experiment: beat GenePT-style text-only baselines).

**19. Cell2Sentence — Teaching Large Language Models the Language of Biology** (Levine et al.; ICML 2024; bioRxiv 2023.09.11.557287). https://pmc.ncbi.nlm.nih.gov/articles/PMC11565894/
- *Takeaway:* turns a cell into a "sentence" of rank-ordered gene names so a plain LLM can be fine-tuned on it.
- *Technical summary:* Rank-orders expressed genes into space-separated "cell sentences"; fine-tunes GPT-2-scale LLMs for cell-type annotation, conditional cell generation, and cell↔text tasks, with an invertible transform back to expression.
- *Plain-English:* Write each cell as a list of its most-active genes, then teach a chatbot to read and write those lists.
- *Applicability:* **A1/A3.** Implication: the sequence-ization trick is a route to make MORPHEUS's transcriptomic modality *natively promptable* by an LM; but rank-ordering discards magnitude — a caution for numeric-modality encoding (A4).
- *Novelty implication:* Establishes "expression-as-text-prompt" prior art; MORPHEUS's imaging+RNA joint state is the differentiator.

**20. Scaling LLMs for Next-Generation Single-Cell Analysis (C2S-Scale)** (Rizvi et al.; bioRxiv 2025.04.14.648850). https://www.biorxiv.org/content/10.1101/2025.04.14.648850v2.full
- *Takeaway:* scales Cell2Sentence to 27B params over 1B+ tokens of transcriptomic + text + metadata.
- *Technical summary:* Trains Gemma-based LLMs on a mixed corpus of cell sentences, biological text, and metadata; reports consistent scaling gains on predictive and generative single-cell tasks and emergent multi-task behavior.
- *Plain-English:* A much bigger "cells-as-sentences" model that gets better and more general as it grows.
- *Applicability:* **A1/A3** (scale-driven emergence of biological capability). Implication: relevant to MORPHEUS's A3 emergence-evaluation — provides a scaling baseline against which "emergent knowledge" claims must be measured, not assumed.
- *Novelty implication:* Raises the bar: emergence-from-scale is expected, so MORPHEUS's A3 novelty must isolate emergence *attributable to the identified/multimodal design*, not to scale.

**21. LangCell — Language-Cell Pre-training for Cell Identity Understanding** (Zhao et al.; ICML 2024). https://arxiv.org/abs/2405.06708
- *Takeaway:* the only single-cell model that does *zero-shot* cell-identity classification via joint cell-text pretraining.
- *Technical summary:* Jointly pretrains a cell encoder and text encoder on cell-data paired with identity-rich descriptions (contrastive + matching + masked-gene objectives); enables zero-shot and strong few-shot cell-type/identity recognition where transcriptome-only FMs need fine-tuning.
- *Plain-English:* Because it learned cells alongside descriptions of what they are, it can name a new cell type it was never explicitly trained to label.
- *Applicability:* **A1/A3** (zero-shot NL labeling = task specified in language) and **A2** (identity-grounded representation). Implication: LangCell is the tightest omics-side analogue to MORPHEUS's promptable state — the design lesson is that *text-grounded pretraining* is what buys zero-shot promptability.
- *Novelty implication:* **Candidate novelty risk for A1** — a zero-shot, NL-promptable *single-cell* representation exists; MORPHEUS's novelty must rest on the *WSI+RNA+clinical multimodal tumor* scope and causal/identifiable structure.

**22. ChatCell — Facilitating Single-Cell Analysis with Natural Language** (Fang et al.; 2024; *preprint later withdrawn*). https://arxiv.org/abs/2402.08303
- *Takeaway:* a conversational front-end for single-cell tasks via vocabulary adaptation + unified sequence generation.
- *Technical summary:* Adapts an LLM vocabulary to gene tokens and casts pseudo-cell generation, cell-type annotation, and drug-response tasks as unified text generation, lowering the coding barrier to single-cell analysis. (Note: authors withdrew the arXiv version; cite with that caveat.)
- *Plain-English:* Do single-cell analysis by chatting instead of coding.
- *Applicability:* **A1.** Implication: an existence proof of a NL front-end to omics prediction — informs the TQI's UX, with the withdrawal a reminder to validate rigor.
- *Novelty implication:* Adds to A1 prior-art breadth (with reliability caveat).

**23. scReader — Prompting Large Language Models to Interpret scRNA-seq Data** (Li et al.; 2024). https://arxiv.org/abs/2412.18156
- *Takeaway:* prompts LLMs with hybrid gene-embedding + text to interpret cross-species single-cell data.
- *Technical summary:* Combines a gene-expression embedding module with LLM prompting to enable cell-type interpretation and cross-species generalization without full omics-FM pretraining; uses natural-language prompts as the task interface.
- *Plain-English:* Feed a language model a compact numeric summary of a cell plus a question, and it interprets the cell.
- *Applicability:* **A4** (numeric embedding *injected into a prompt* — an encode-into-context hybrid) and **A1.** Implication: concrete design point for MORPHEUS's "embed the modality as a soft prompt token" option in the encode-vs-RAG spectrum.
- *Novelty implication:* Informs A4; the soft-prompt-injection pattern is established.

---

### Group III — Biological sequence ⇄ natural language (DNA/RNA/protein prompting)

**24. ChatNT — A Multimodal Conversational Agent for DNA, RNA and Protein Tasks** (de Almeida et al.; Nature Machine Intelligence 2025; bioRxiv 2024.04.30.591835). https://www.biorxiv.org/content/10.1101/2024.04.30.591835v2
- *Takeaway:* one English-conversation model solves 18+ genomics/transcriptomics/proteomics tasks with a single architecture.
- *Technical summary:* Couples a frozen Nucleotide-Transformer DNA encoder to an LLM via a projection so users pose classification/regression tasks in English; sets SOTA on the NT benchmark (avg MCC 0.77, +8 pts over NTv2-500M) solving all tasks simultaneously, and is extensible to longer encoders/other modalities.
- *Plain-English:* Ask questions about DNA/RNA/protein in plain English and one model answers dozens of different scientific questions.
- *Applicability:* **A1** (the flagship "many tasks, one NL interface, single frozen encoder" system) and **A4** (frozen sequence-encoder plug-in). Implication: ChatNT is the single most important comparator for MORPHEUS's A1 thesis — it proves the *frozen-encoder + NL-task interface* is achievable and benchmark-competitive.
- *Novelty implication:* **Major candidate novelty risk for A1.** A frozen-trunk, NL-promptable, multi-omics-task agent already exists and wins benchmarks. MORPHEUS's A1/A4 novelty must be sharply relocated to (i) *WSI+molecular tumor-state* scope, (ii) *task auto-detection with abstention* (ChatNT is told the task in the question), and (iii) *identifiable pathway slots* — otherwise A1 is largely pre-empted.

**25. ProtST — Multi-Modality Learning of Protein Sequences and Biomedical Texts** (Xu et al.; ICML 2023, oral). https://arxiv.org/abs/2301.12040
- *Takeaway:* aligns protein-sequence encoders with biomedical text for supervised and *zero-shot* protein classification.
- *Technical summary:* Builds ProtDescribe (protein-sequence↔text-description pairs) and trains with multimodal contrastive + masked objectives, enabling zero-shot protein function prediction from text prompts and boosting supervised property tasks.
- *Plain-English:* Pair protein sequences with descriptions so the model can classify a new protein's function from a text prompt alone.
- *Applicability:* **A3/A1.** Implication: template for building MORPHEUS's paired (tumor-state ↔ description) corpus and for zero-shot functional prompting; the ProtDescribe construction is a reusable recipe.
- *Novelty implication:* Reinforces that text-grounded zero-shot functional prompting is established across modalities.

**26. InstructProtein — Aligning Human and Protein Language via Knowledge Instruction** (Wang et al.; ACL 2024; arXiv 2023). https://arxiv.org/abs/2310.03269
- *Takeaway:* uses a knowledge-graph-driven instruction generator to align protein and human language bidirectionally.
- *Technical summary:* Pretrains on protein + NL corpora, then instruction-tunes with KG-derived instructions to reduce annotation bias, enabling protein→text (function description) and text→protein (design) generation.
- *Plain-English:* Teaches a model to talk about proteins both ways — describe one, or design one from a description — using a knowledge graph to make the lessons.
- *Applicability:* **A2/A3** (KG-structured instructions ≈ pathway-addressable grounding). Implication: MORPHEUS could use pathway knowledge graphs to *generate* prompt/instruction data addressing specific programme slots (A2).
- *Novelty implication:* Informs A2/A3; KG-grounded instruction is a known bias-reduction tool.

**27. Prot2Text — Multimodal Protein's Function Generation with GNNs and Transformers** (Abdine et al.; AAAI 2024). https://arxiv.org/abs/2307.14367
- *Takeaway:* generates free-text protein function descriptions from sequence+structure via GNN-encoder → text-decoder.
- *Technical summary:* Encoder-decoder fusing sequence, ESM features, and structural graph, decoding to natural-language function paragraphs; introduces a generative (not classification) protein-function benchmark.
- *Plain-English:* Give it a protein and it writes a paragraph describing what the protein does.
- *Applicability:* **A3.** Implication: a concrete "biology→NL description" head design and an evaluation of *generated* biological knowledge — relevant to MORPHEUS's flexible NL output head and to A3's measurement problem.
- *Novelty implication:* Reframes A3 evaluation: free-text generation needs generative metrics (and the risk of fluent-but-wrong), sharpening MORPHEUS's emergence-eval design.

---

### Group IV — Instruction-tuned scientific LLMs & multimodal biomedical assistants (context/comparators)

**28. BioMedGPT — Open Multimodal Generative Pre-trained Transformer for BioMedicine** (Luo et al.; 2023). https://arxiv.org/abs/2308.09442
- *Takeaway:* aligns molecule and protein encoders to an LLM through natural language into one biomedical assistant.
- *Technical summary:* Bridges 2-D molecule and protein FMs to an LLM via cross-modal feature aligners trained on text, enabling free-text QA over molecules/proteins and competitive results on biomedical QA benchmarks.
- *Plain-English:* A single chatbot that reasons across molecules, proteins, and the biomedical literature.
- *Applicability:* **A1/A4** (align frozen domain encoders into an LLM — the encode-and-project pattern). Implication: direct architectural precedent for MORPHEUS attaching a tumor-state encoder to an NL interface.
- *Novelty implication:* Adds to A1/A4 prior art; "attach FM encoders to an LLM for a biomedical assistant" is now standard.

**29. BioInstruct — Instruction Tuning of LLMs for Biomedical Natural Language Processing** (Tran et al.; JAMIA 2024). https://arxiv.org/abs/2310.19975
- *Takeaway:* a 25k-instruction corpus (GPT-4-bootstrapped) that lifts biomedical QA/IE/generation.
- *Technical summary:* Self-generated biomedical instructions from 3-seed prompting of GPT-4; instruction-tunes LLaMA-1/2, yielding +17.3% QA, +5.7% IE, large generation gains over untuned baselines.
- *Plain-English:* A workbook of biomedical tasks that makes general chatbots much better at medical text.
- *Applicability:* **A3** (instruction data for biomedical grounding). Implication: recipe for MORPHEUS's instruction corpus; also a baseline showing text-only tuning's ceiling.
- *Novelty implication:* Establishes biomedical instruction tuning as commodity — not a novelty vector for MORPHEUS.

**30. Biology-Instructions — A Dataset and Benchmark for Multi-Omics Sequence Understanding of LLMs** (He et al.; 2024). https://arxiv.org/abs/2412.19191
- *Takeaway:* a benchmark exposing that current LLMs are weak at multi-omics sequence tasks even with prompting.
- *Technical summary:* Curates multi-omics (DNA/RNA/protein) instruction tasks and shows off-the-shelf and naively-tuned LLMs struggle without dedicated encoders, motivating a two-stage encoder-plus-LLM recipe.
- *Plain-English:* A test set proving that chatbots alone can't read raw omics sequences well — they need a specialist encoder.
- *Applicability:* **A4/A1** (evidence that some modalities must be *encoded*, not prompted as raw text). Implication: empirical support for MORPHEUS's A4 thesis that structured omics should be *encoded* into slots rather than passed as text/RAG.
- *Novelty implication:* **Strengthens A4** — provides a citable result that raw-sequence-as-text prompting fails, legitimizing the encode-vs-retrieve question as a real open problem.

---

### Group V — Predicting molecular/pathway state from images via "prompting" foundation models

**31. SEQUOIA — Digital profiling of gene expression from histology images with linearized attention** (Zhu et al.; Nature Communications 2024). https://www.nature.com/articles/s41467-024-54182-5
- *Takeaway:* predicts bulk/spatial transcriptomic profiles from H&E whole-slide images atop a pathology FM.
- *Technical summary:* Uses UNI histology-FM features with a linearized-attention transformer to model whole-slide context and regress gene-expression profiles; generalizes across tissue types and recovers expression-based prognostic signal.
- *Plain-English:* Reads a tissue slide and estimates which genes are switched on, without doing the sequencing.
- *Applicability:* **A4** (image→molecular *prediction* — the imaging modality carries molecular state) and **A1** (a hard-coded probe, not a prompt). Implication: SEQUOIA is exactly the *current-MORPHEUS-style* fixed molecular probe; the rebase delta is making such queries *promptable* rather than a bespoke regressor.
- *Novelty implication:* Reframes MORPHEUS's contribution: image→expression is solved as a fixed model; novelty is turning it into an *addressable, promptable pathway query* (A1+A2), not the prediction itself.

**32. FmH2ST — Foundation model-based spatial transcriptomics generation from histological images** (2025; Nucleic Acids Research, gkaf865). https://academic.oup.com/nar/article/53/17/gkaf865/8249850
- *Takeaway:* dual-branch FM + spot-image model to generate spatial gene expression from histology.
- *Technical summary:* Integrates prior knowledge from a histology foundation model with fine-grained spot-image detail in a dual-branch architecture to predict spatially-resolved expression, improving generalization across platforms.
- *Plain-English:* Turns a tissue image into a map of where each gene is active.
- *Applicability:* **A4/A2** (spatially-addressable molecular readout from imaging). Implication: informs whether MORPHEUS's molecular slots can be made *spatial* and image-derived.
- *Novelty implication:* Adds to the image→molecular prior art; the *fixed-task* framing again motivates the promptable delta.

**33. PEKA — Teaching pathology foundation models to accurately predict gene expression with parameter-efficient knowledge transfer** (2025). https://arxiv.org/abs/2504.07061
- *Takeaway:* a parameter-efficient adapter (Block-Affine Adaptation + distillation) transfers a frozen pathology FM to expression prediction.
- *Technical summary:* Freezes the pathology FM and trains Block-Affine adapters with knowledge-distillation and structure-alignment losses for cross-modal (image→expression) transfer, beating full fine-tuning at a fraction of parameters.
- *Plain-English:* Cheaply teach an existing tissue-image model to also predict gene activity, without retraining it.
- *Applicability:* **A4** (frozen-trunk + parameter-efficient plug-in for a *new molecular task*). Implication: strong precedent for MORPHEUS's frozen-representation-plus-adapter strategy in the imaging→molecular direction.
- *Novelty implication:* Pre-empts "adapter transfers a frozen pathology FM to molecular readouts is new"; MORPHEUS must own the *promptable/identified* framing, not the transfer mechanics.

**34. Diffusion Generative Modeling for Spatially Resolved Gene Expression Inference from Histology Images (Stem)** (2025). https://arxiv.org/abs/2501.15598
- *Takeaway:* a conditional diffusion model captures the multimodal *distribution* of expression given an image patch.
- *Technical summary:* Conditions a diffusion generator on histology features to sample spatial expression, modeling uncertainty/multimodality rather than a single point estimate; improves calibration over regression baselines.
- *Plain-English:* Instead of one guess of gene activity per tissue spot, it produces a realistic range of possibilities.
- *Applicability:* **A5/A4** (generative, distributional molecular prediction ≈ substrate for counterfactual sampling). Implication: relevant to MORPHEUS A5 — a *generative* molecular head could support "what-if" sampling rather than a point classifier.
- *Novelty implication:* Informs A5's generative-query framing; distribution modeling is established, so novelty is the *interventional/causal* conditioning.

**35. CellSymphony — Deciphering molecular and phenotypic orchestration of cells with single-cell pathomics** (2025). https://arxiv.org/abs/2508.10232
- *Takeaway:* fuses single-cell spatial-omics with pathology-FM image embeddings for cell-level molecular+phenotype states.
- *Technical summary:* Aligns single-cell expression with morphology embeddings from a pathology FM to characterize each cell's molecular and phenotypic program in spatial context.
- *Plain-English:* Matches what each cell looks like to what it is doing molecularly, cell by cell.
- *Applicability:* **A2/A4** (cell-level, programme-addressable multimodal state). Implication: closest to MORPHEUS's per-programme slot ambition at single-cell resolution; a design reference for addressable biological units.
- *Novelty implication:* Adds to A2/A4 prior art on morphology↔molecular fusion; MORPHEUS's promptable/causal layer remains the differentiator.

---

### Group VI — Surveys / landscape (orientation only)

**36. Leveraging Biomolecule and Natural Language through Multi-Modal Learning: A Survey** (Pei et al.; 2024). https://arxiv.org/abs/2403.01528
- *Takeaway:* the reference map of molecule/protein/genome ⇄ language methods, taxonomizing alignment strategies.
- *Technical summary:* Surveys contrastive, generative, and instruction-tuned biomolecule-text models across tasks and modalities; organizes the field into alignment paradigms and identifies open gaps (multi-omics, evaluation).
- *Plain-English:* A big-picture review of every way people connect biology to language models.
- *Applicability:* **A1-A4** (landscape). Implication: use to bound MORPHEUS's related-work claims and ensure no prior system is missed.
- *Novelty implication:* Meta-evidence that the molecule/protein-text space is crowded; MORPHEUS's whitespace is the *multimodal tumor-state + causal/identified* combination the survey does not cover.

**37. LLM4Cell — A Survey of Large Language and Agentic Models for Single-Cell Biology** (2025). https://arxiv.org/abs/2510.07793
- *Takeaway:* recent survey of LLM/agentic approaches to single-cell biology, including cell-text and prompting methods.
- *Technical summary:* Catalogs cell-language pretraining, cell-sentence methods, prompting/instruction approaches, and agentic pipelines for single-cell tasks, with evaluation-gap discussion.
- *Plain-English:* A current review of using language models (and agents) for single-cell data.
- *Applicability:* **A1/A3** (landscape for the omics-NL half of the lane). Implication: identifies which cell-NL capabilities are saturated vs open, informing where MORPHEUS should not re-claim novelty.
- *Novelty implication:* Confirms cell-NL prompting is an active, crowded area — sharpens the emergence-evaluation (A3) as the under-served problem.

---

### Cross-cutting synthesis for MORPHEUS

- **A1 is the most pre-empted axis.** ChatNT (#24), Text+Chem T5 (#9), and LangCell (#21) already deliver frozen-trunk, NL-promptable, multi-task biological models — several told-the-task, but ChatNT + C2S-Scale approach implicit multi-task competence. MORPHEUS's defensible A1 novelty narrows to **task auto-detection with abstention over a WSI+molecular tumor state**, not "an NL interface to a biological model."
- **A4 has a genuine, citable open problem.** Biology-Instructions (#30) shows raw-sequence-as-text prompting fails; ChatDrug (#16)/ChemLML (#12)/PEKA (#33) span the encode↔retrieve↔adapter spectrum. No lane paper formalizes *which modality to encode vs retrieve* — this is the cleanest whitespace.
- **A3 emergence-evaluation is under-served.** GenePT (#18) is the key adversary: LLM literature priors already encode much "biology," so MORPHEUS must prove emergent knowledge *beyond text priors and beyond scale* (C2S-Scale #20) — a decisive experiment, not an assertion.
- **A5's NL-intervention framing exists in chemical space** (MoleculeSTM #4, ChatDrug #16, Stem #34) but **not as post-encoding causal-geometry queries on a frozen tumor state** — MORPHEUS's strongest remaining novelty foothold in this lane.
- **A2 (identifiable pathway-addressable slots) is the least-covered here**; MolCA's Q-Former query-tokens (#7) and CellSymphony (#35) hint at addressability but none claim *identifiability*. Cross-lane with l09/l12 for the theory.
