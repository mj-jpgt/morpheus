## Multimodal RAG vs. encoded-modality integration

*Lane l08 — Retrieval-augmented multimodal/biomedical models; when retrieval beats parametric
encoding; the kNN-LM/RETRO lineage; retrieval-augmented biology; the encode-vs-retrieve tradeoff for a
modality. Maps primarily to MORPHEUS axis **A4** (which modalities to ENCODE vs. treat as RAG context,
frozen-trunk plug-in), with secondary hits on A1 (when-to-retrieve routing / task inference), A2
(addressable slots), A3 (grounding/elicitation), and A5 (retrieval-as-query / intervention).*

Scope note: entries are ordered roughly foundational-lineage → encode-vs-retrieve theory →
biomedical/multimodal RAG → retrieval-augmented biology. Every citation was fetched and verified; no
fabricated references.

---

### Part A — The kNN-LM / RETRO / RALM lineage (the retrieval-vs-parametric backbone)

**1. Generalization through Memorization: Nearest Neighbor Language Models (kNN-LM)**
Khandelwal, Levy, Jurafsky, Zettlemoyer, Lewis. ICLR 2020. arXiv:1911.00172.
- *Takeaway:* You can improve a frozen LM at inference with zero retraining by interpolating its
  softmax with a k-NN lookup over a datastore of hidden-state→next-token pairs.
- *Technical summary:* A single forward pass over a corpus builds a (key=context embedding,
  value=next token) datastore queried at test time via FAISS; the final distribution linearly mixes
  the parametric softmax with the kNN distribution. On WikiText-103 it set SOTA perplexity (15.79,
  −2.9) with no gradient updates, and enabled domain adaptation by simply swapping the datastore.
  Gains concentrate on rare/long-tail and factual tokens.
- *Plain-English:* Instead of memorizing everything in its weights, the model keeps a searchable
  "notebook" of things it has seen and looks up close matches when unsure — which helps most on rare
  facts.
- *Applicability (A4, A1):* The founding evidence that a **frozen trunk + external datastore** beats
  pushing everything into parameters, *especially for the long tail*. For MORPHEUS this is the
  theoretical license for A4's frozen-trunk plug-in: a rarely-measured modality (phospho, CPTAC
  proteomics) with few training samples is exactly the "rare pattern" regime where retrieval, not
  encoding, wins. Datastore-swap = per-cohort adaptation without retraining the trunk.
- *Novelty implication:* **Strengthens** A4's framing but **pre-empts** any claim that "retrieve vs.
  encode" is itself novel — the tradeoff is 6 years old in NLP. MORPHEUS's novelty must be the
  *biological modality-selection rule*, not the mechanism.

**2. Improving Language Models by Retrieving from Trillions of Tokens (RETRO)**
Borgeaud, Mensch, Hoffmann, ... Sifre (DeepMind). ICML 2022. arXiv:2112.04426.
- *Takeaway:* Chunked cross-attention over retrieved neighbors lets a 7.5B model match GPT-3 (175B)
  on the Pile with 25× fewer parameters.
- *Technical summary:* A frozen BERT retriever fetches nearest-neighbor chunks from a 2T-token
  database; a differentiable encoder + interleaved chunked cross-attention conditions generation on
  them. Performance improves log-linearly with database size up to 2T tokens, decoupling knowledge
  capacity from parameter count.
- *Plain-English:* A small model that looks things up in a giant library can be as knowledgeable as a
  huge model that memorized everything.
- *Applicability (A4):* Direct architectural template for "encode the reasoning core small, retrieve
  the knowledge." Motivates keeping MORPHEUS's WSI+RNA trunk compact while a proteomic/pathway
  datastore supplies rare-programme evidence via cross-attention rather than adapters.
- *Novelty implication:* **Pre-empts** "retrieval scales knowledge without retraining" as a MORPHEUS
  contribution; **reframes** A4 toward the harder, unclaimed question — *which biological modalities
  are knowledge-like (retrieve) vs. reasoning-like (encode)?*

**3. REALM: Retrieval-Augmented Language Model Pre-Training**
Guu, Lee, Tung, Pasupat, Chang. ICML 2020. arXiv:2002.08909.
- *Takeaway:* First to make the retriever *end-to-end trainable during pre-training* via a latent
  knowledge-retrieval step, then freeze it for downstream QA.
- *Technical summary:* A masked-LM objective backpropagates through a maximum-inner-product-search
  retriever so the model learns *what* to retrieve; fine-tuning on Open-QA keeps the retriever frozen.
  Establishes the parametric-memory (weights) vs. non-parametric-memory (corpus) decomposition.
- *Plain-English:* The model is trained not just to answer, but to learn which documents are worth
  fetching in the first place.
- *Applicability (A1, A4):* The learned-retriever precedent for MORPHEUS's router: deciding *whether
  and what* to retrieve is trainable, not hand-set. Relevant to A1's task-inference routing.
- *Novelty implication:* **Pre-empts** "trainable retrieval into a frozen core" novelty; supports a
  design where MORPHEUS's scope/abstention router is learned end-to-end.

**4. Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks (RAG)**
Lewis, Perez, Piktus, ... Kiela. NeurIPS 2020. arXiv:2005.11401.
- *Takeaway:* The canonical parametric+non-parametric hybrid — a seq2seq generator conditioned on
  dense-retrieved passages — that named the "RAG" paradigm.
- *Technical summary:* DPR retrieves Wikipedia passages; BART marginalizes over them (RAG-Sequence /
  RAG-Token). Beats parametric-only seq2seq and task-specific retrieve-extract pipelines on open-domain
  QA, and generates more factual, specific, diverse text. The knowledge index is hot-swappable.
- *Plain-English:* Bolt a search engine onto a generator and it becomes more factual and updatable
  without retraining.
- *Applicability (A4, A3):* The reference definition of "modality as RAG context." MORPHEUS's
  safety-gated closed-RAG card renderer is already a constrained RAG; A4 asks which molecular
  modalities *should* enter as retrieved context vs. encoded features. A3: hot-swap index = updatable
  grounding.
- *Novelty implication:* **Pre-empts** RAG-as-mechanism novelty entirely; MORPHEUS cannot claim RAG,
  only a novel *biological application/selection criterion*.

**5. Atlas: Few-shot Learning with Retrieval Augmented Language Models**
Izacard, Lewis, Lomeli, ... Grave. JMLR 2023 / arXiv:2208.03299.
- *Takeaway:* Retrieval-augmented models can be strongly *few-shot* — 42% on NaturalQuestions with 64
  examples, beating a 540B model with 50× fewer parameters.
- *Technical summary:* Co-trains a Contriever retriever + Fusion-in-Decoder generator with
  retrieval-in-the-loop pretraining; studies index content and shows the datastore can be edited to
  update knowledge. Demonstrates sample-efficiency, not just capacity, from non-parametric memory.
- *Plain-English:* Looking things up doesn't just add facts — it lets a model learn new tasks from a
  handful of examples.
- *Applicability (A4, A2):* Directly relevant to MORPHEUS's data-scarce modalities (CPTAC proteomics
  is inventory-only): retrieval buys *few-shot* competence where an encoder would overfit. A2:
  editable index ≈ addressable knowledge slots.
- *Novelty implication:* **Strengthens** the A4 thesis that scarce modalities favor retrieval; the
  *few-shot* argument is the sharpest quantitative case MORPHEUS can borrow.

**6. Memorizing Transformers**
Wu, Rabe, Hutchins, Szegedy. ICLR 2022 (Spotlight). arXiv:2203.08913.
- *Takeaway:* An internal kNN memory of past (key,value) activations, read via attention, lets a model
  use newly defined facts/functions at test time without weight updates.
- *Technical summary:* Extends attention with a non-differentiable kNN lookup into a large memory (up
  to 262K tokens); benefits grow with memory size across code, math, and formal-theorem domains, where
  the model exploits functions/theorems defined only at inference.
- *Plain-English:* The model keeps a scratchpad of what it just saw and can attend back to it, so it
  "learns" new definitions on the fly.
- *Applicability (A4, A5):* Template for treating a modality as *memory read via attention* rather than
  a fused input — the middle ground between encode and RAG. Relevant to A5: newly specified
  perturbations/definitions usable at query time.
- *Novelty implication:* **Reframes** A4 as a spectrum (fuse ↔ attend-to-memory ↔ retrieve-text), not
  a binary; MORPHEUS can claim novelty in *placing biological modalities on this spectrum*.

**7. In-Context Retrieval-Augmented Language Models (In-Context RALM)**
Ram, Levine, Dalmedigos, ... Shoham. TACL 2023. arXiv:2302.00083.
- *Takeaway:* Simply prepending retrieved documents to a frozen LM's context — no architecture change,
  no training — yields large gains across model sizes.
- *Technical summary:* Off-the-shelf retrievers + document prepending; a lightweight rerank tuned to
  the LM setting further helps. Establishes the cheapest possible RAG and a strong training-free
  baseline.
- *Plain-English:* Often you don't need to modify the model at all — just paste the relevant documents
  into the prompt.
- *Applicability (A4, A1):* The zero-surgery baseline any MORPHEUS "encode the modality" proposal must
  beat: if prepending retrieved proteomic evidence to the closed-RAG renderer matches an encoded
  adapter, the adapter isn't earning its place (guardrail: earn-your-place modality gating).
- *Novelty implication:* **Pre-empts** claims for elaborate fusion by setting a hard, cheap baseline —
  a useful adversarial control for A4's encode-vs-retrieve study.

---

### Part B — The encode-vs-retrieve tradeoff and *why* retrieval works (A4 core theory)

**8. Nonparametric Masked Language Modeling (NPM)**
Min, Shi, Lewis, Chen, Yih, Hajishirzi, Zettlemoyer. ACL 2023 Findings. arXiv:2212.01349.
- *Takeaway:* Replace the softmax over a fixed vocabulary with a *nonparametric* distribution over
  every phrase in a reference corpus — retrieval literally *replaces* parametric prediction.
- *Technical summary:* Contrastively trained encoder predicts masked spans by retrieving phrases from
  a corpus (no output vocabulary matrix). Across 16 zero-shot tasks it beats significantly larger
  parametric models, with the largest wins on rare words, rare senses, rare facts, and non-Latin
  scripts.
- *Plain-English:* Instead of memorizing a fixed dictionary, the model points at the exact phrase it
  needs in a library — which shines on rare things big models fumble.
- *Applicability (A4):* The strongest single argument for A4's biological version: for a
  *sparsely-observed, heavy-tailed* modality (mutation combinations, rare phospho-sites), a
  nonparametric/retrieval head can out-generalize an encoder that must allocate parameters to the tail.
  Suggests a concrete MORPHEUS experiment: encoded proteomic adapter vs. nonparametric retrieval head,
  scored on rare-programme recall.
- *Novelty implication:* **Reframes** A4 from "add RAG context" to "retrieve *instead of* encode a
  modality" — a sharper, more defensible MORPHEUS claim, but one NPM already owns in NLP; the
  biological-modality instantiation is where novelty must live.

**9. Why do Nearest Neighbor Language Models Work?**
Xu, Alon, Neubig. arXiv:2301.02828 (2023).
- *Takeaway:* kNN-LM gains come from three mechanistic factors (a different input representation for
  the final prediction, the softmax-temperature, and ensembling via approximate search) — and can
  largely be *internalized into the parametric model*.
- *Technical summary:* Ablations attribute the improvement to using a distinct representation for the
  next-token distribution, temperature calibration, and stochastic approximate-kNN ensembling; folding
  these into standard LM training recovers much of the benefit *without* an explicit datastore.
- *Plain-English:* The magic of "looking it up" is partly a trick you can bake back into the model —
  retrieval isn't always fundamentally necessary.
- *Applicability (A4):* Critical devil's-advocate for A4: before MORPHEUS declares a modality
  "retrieve-only," it must rule out that a better-trained encoder captures the same signal. Motivates a
  control arm where the retrieval benefit is distilled into the trunk.
- *Novelty implication:* **Reframes / partially pre-empts** the encode-vs-retrieve dichotomy —
  retrieval advantage may be an artifact of under-trained encoders. Raises the evidentiary bar for any
  MORPHEUS "must-retrieve" claim; a genuine novelty needs a modality where the gap *survives*
  distillation.

**10. Fine-Tuning or Retrieval? Comparing Knowledge Injection in LLMs**
Ovadia, Brief, Mishaeli, Elisha. arXiv:2312.05934 (2023/24); EMNLP 2024.
- *Takeaway:* For injecting *new* knowledge, RAG consistently beats unsupervised fine-tuning — models
  struggle to absorb new facts via fine-tuning unless massively paraphrased.
- *Technical summary:* Controlled study on knowledge-intensive tasks contrasts continued fine-tuning
  vs. RAG for both seen and unseen facts; RAG wins across the board and avoids catastrophic forgetting,
  while fine-tuning needs many paraphrase variants to learn a single fact.
- *Plain-English:* If you want a model to *know* something new, giving it a lookup usually beats
  retraining it.
- *Applicability (A4):* Direct decision rule for MORPHEUS: a new modality/knowledge source added
  *after* the trunk is frozen is better served by retrieval than by re-tuning — supporting A4's
  frozen-trunk-plug-in stance over adapter fine-tuning for late-arriving data (e.g., a new CPTAC
  cohort).
- *Novelty implication:* **Strengthens** the "freeze the trunk, retrieve the new modality" design and
  gives it empirical backing; **pre-empts** any claim that this preference is MORPHEUS-original.

**11. Retrieval-Enhanced Machine Learning (REML)**
Zamani, Diaz, Dehghani, Metzler, Bendersky. SIGIR 2022. arXiv:2205.01230.
- *Takeaway:* Generalizes retrieval-augmentation beyond text: any ML model can be an "end user" of an
  IR system, gaining generalization, scalability, robustness, and interpretability.
- *Technical summary:* Proposes a framework (indexing/representation/retrieval/ranking for models, not
  humans) that subsumes kNN-LM/RAG as special cases and articulates open problems in optimizing
  retrieval *for downstream model consumption*.
- *Plain-English:* The idea of "look it up instead of memorizing" isn't just for chatbots — it's a
  general design pattern for any model.
- *Applicability (A4, A3):* Gives MORPHEUS the vocabulary to argue A4 as a principled, formalizable
  design axis rather than an engineering hack — a modality is a "retrieval-enhanced input." A3:
  interpretability from surfacing retrieved evidence.
- *Novelty implication:* **Reframes** A4 as an instance of a known general framework; MORPHEUS's
  contribution is the *biomedical-multimodal specialization + a modality-selection criterion*, which
  REML explicitly lists as open.

**12. Self-RAG: Learning to Retrieve, Generate, and Critique through Self-Reflection**
Asai, Wu, Wang, Sil, Hajishirzi. arXiv:2310.11511 (2023); ICLR 2024.
- *Takeaway:* Train the model to decide *on demand* whether to retrieve, and to critique retrieved
  passages and its own output, via reflection tokens.
- *Technical summary:* Special tokens control adaptive retrieval, relevance/critique judgments, and
  supportedness at inference; 7B/13B Self-RAG beats ChatGPT and retrieval-augmented Llama2-chat on QA,
  reasoning, and long-form factuality/citation.
- *Plain-English:* The model learns when it actually needs to look something up, and to check whether
  what it found is trustworthy.
- *Applicability (A1, A4):* The clearest template for A1/A4's *routing* decision — MORPHEUS's
  scope/abstention router as a learned "retrieve-or-not / encode-or-retrieve" gate, with self-critique
  feeding the safety gate on closed-RAG cards.
- *Novelty implication:* **Pre-empts** "adaptive when-to-retrieve" as novel in general; **strengthens**
  the feasibility of MORPHEUS's router, but MORPHEUS must specialize the decision to *modalities*, not
  passages.

**13. Retrieval-Augmented Generation for Large Language Models: A Survey**
Gao, Xiong, Gao, ... Wang. arXiv:2312.10997 (2023/24).
- *Takeaway:* Canonical taxonomy — Naive / Advanced / Modular RAG over a retrieval–generation–
  augmentation tripartite — plus evaluation frameworks and open problems.
- *Technical summary:* Systematizes chunking, indexing, query transformation, reranking, fusion, and
  iterative/adaptive retrieval; catalogs benchmarks and failure modes (hallucination under irrelevant
  retrieval, lost-in-the-middle).
- *Plain-English:* A map of the whole RAG design space and how to evaluate it.
- *Applicability (A4, A3):* Design checklist for MORPHEUS's molecular-RAG: where to place a modality on
  the Naive→Modular axis, and which evaluation pitfalls to pre-register in the encode-vs-retrieve study.
- *Novelty implication:* **Pre-empts** generic RAG-architecture claims; useful mainly to position
  MORPHEUS *against* known design points and avoid reinventing them.

---

### Part C — Multimodal & biomedical RAG (the marquee A4 application surface)

**14. Retrieval-Augmented Multimodal Language Modeling (RA-CM3)**
Yasunaga, Aghajanyan, Shi, ... Yih. ICML 2023. arXiv:2211.12561.
- *Takeaway:* First retrieval-augmented model that both *retrieves and generates* across text+images,
  letting a base multimodal model refer to external multimodal memory instead of memorizing it.
- *Technical summary:* CLIP retriever + CM3 Transformer generator over LAION; retrieves relevant
  text/images into context. Gains 12 FID / 17 CIDEr on MS-COCO at lower training compute, plus
  multimodal in-context learning and faithful generation.
- *Plain-English:* An image-and-text model that fetches relevant pictures and captions instead of
  cramming all of them into its weights.
- *Applicability (A4):* The direct multimodal precedent for MORPHEUS: a modality (image *or* text) can
  be retrieved into a shared model rather than encoded. Supports treating WSI-adjacent evidence
  (reference tiles, matched reports) as retrieved multimodal context.
- *Novelty implication:* **Pre-empts** "multimodal retrieval" as novel; MORPHEUS's edge is *molecular*
  modalities (proteomics/phospho/CNV) with no natural text/image analog — retrieval over structured
  biological memory is the unclaimed territory.

**15. Benchmarking Retrieval-Augmented Generation for Medicine (MedRAG / MIRAGE)**
Xiong, Jin, Lu, Zhang. arXiv:2402.13178 (2024).
- *Takeaway:* Systematic 41-way benchmark (corpora × retrievers × LLMs) on 7,663 medical questions;
  RAG lifts GPT-3.5/Mixtral to GPT-4 level (+up to 18% over CoT).
- *Technical summary:* MIRAGE benchmark + MedRAG toolkit; finds combining heterogeneous medical corpora
  and retrievers is best, documents log-linear corpus scaling and a medical "lost-in-the-middle"
  effect.
- *Plain-English:* Carefully tested recipe book showing that medical lookups reliably make models more
  accurate — up to a point, and with pitfalls.
- *Applicability (A3, A4):* The evaluation-design reference for MORPHEUS's grounding/elicitation
  benchmark (A3) and for how to run a rigorous encode-vs-retrieve comparison (A4) with confound-aware
  controls (the lane's guardrail).
- *Novelty implication:* **Pre-empts** naive "RAG helps medicine" claims; **strengthens** the need for
  MORPHEUS's confound-aware evaluation harness by exposing corpus/position confounds.

**16. MMed-RAG: Versatile Multimodal RAG for Medical Vision-Language Models**
Xia, Zhu, Li, ... Yao. arXiv:2410.13085 (2024); ICLR 2025.
- *Takeaway:* Domain-aware retrieval + adaptive context selection + RAG-theory-grounded preference
  tuning cut Med-LVLM factual errors by 43.8% across radiology/ophthalmology/pathology.
- *Technical summary:* Three components — a specialty-aware retriever, an adaptive selector for *how
  many* contexts to use, and preference fine-tuning that provably balances parametric vs. retrieved
  reliance; evaluated on 5 datasets and VQA + report generation.
- *Plain-English:* A medical image assistant that knows which specialty's references to pull, how many,
  and how much to trust them.
- *Applicability (A1, A4):* "How much to trust retrieved vs. internal knowledge" is exactly MORPHEUS's
  encode-vs-retrieve *arbitration*; the adaptive count-selection maps to A1's routing and the safety
  gate.
- *Novelty implication:* **Pre-empts** adaptive multimodal-medical RAG as novel; MORPHEUS must move the
  contribution from *imaging+text* to *molecular* modality arbitration.

**17. RULE: Reliable Multimodal RAG for Factuality in Medical Vision-Language Models**
Xia, Zhu, Li, ... Yao. EMNLP 2024. arXiv:2407.05131.
- *Takeaway:* Calibrated selection of the *number* of retrieved contexts + preference tuning against
  over-reliance yields +47.4% factual accuracy.
- *Technical summary:* Statistically calibrates retrieval-count to control factuality risk; builds a
  preference set from over-reliance failures to fine-tune the balance between inherent and retrieved
  knowledge.
- *Plain-English:* Retrieving *more* isn't always better — this calibrates the right amount and teaches
  the model not to blindly follow it.
- *Applicability (A4):* Directly informs the failure mode MORPHEUS's guardrails target — a modality's
  retrieved context can *hurt* if over-weighted; supports a calibrated gate rather than always-on
  fusion.
- *Novelty implication:* **Strengthens** the "retrieval can degrade if unmanaged" argument central to
  A4's honest tradeoff; the calibration mechanism is prior art to build on, not claim.

**18. HeteroRAG: Heterogeneous Retrieval-Augmented Generation for Medical Vision-Language Tasks**
Chen, Liao, Zhu, ... Wang. ACL 2026 Findings. arXiv:2508.12778.
- *Takeaway:* Retrieve from *heterogeneous* sources (multimodal reports + varied text corpora) via
  modality-specific CLIPs and a multi-corpora query generator, with preference tuning to align sources.
- *Technical summary:* Introduces MedAtlas (mixed report/text repositories), modality-specific
  retrievers, a query generator adapting queries per source, and Heterogeneous Knowledge Preference
  Tuning; SOTA over 11 datasets / 3 imaging modalities.
- *Plain-English:* A medical model that pulls from many different kinds of references at once and
  reconciles them.
- *Applicability (A4, A2):* The multi-source retrieval design MORPHEUS needs if different molecular
  modalities live in *different* datastores (proteomic vs. mutational vs. literature) — per-source
  retrievers ≈ A2's pathway-addressable slots realized as addressable stores.
- *Novelty implication:* **Pre-empts** heterogeneous-source medical RAG; **reframes** A2 toward
  "addressable *retrieval* slots," a bridge MORPHEUS can exploit between A2 and A4.

**19. FactMM-RAG: Fact-Aware Multimodal Retrieval Augmentation for Radiology Report Generation**
Sun, Zhao, Han, Xiong. NAACL 2025. arXiv:2407.15268.
- *Takeaway:* Mine factually-grounded report pairs (via RadGraph) to train a *fact-aware* multimodal
  retriever; factual supervision transfers to the generator without diagnostic labels.
- *Technical summary:* Two stages — RadGraph-based factual pair mining, then universal multimodal
  retriever training; +6.5% F1CheXbert / +2% F1RadGraph, showing factual retrieval supervision needs no
  explicit labels.
- *Plain-English:* Teach the "look-up" step to prize factual correctness, and the whole system writes
  more accurate reports.
- *Applicability (A3, A4):* A3's grounding/elicitation: the retriever is trained on *biological fact
  structure* (graph relations), not just similarity — a template for MORPHEUS retrieving
  pathway-consistent evidence and evaluating factual grounding.
- *Novelty implication:* **Strengthens** the idea that retrieval keys should encode *biological fact
  structure*; a MORPHEUS pathway-fact retriever would be a novel instantiation, not a novel mechanism.

**20. RA-RRG: Multimodal Retrieval-Augmented Radiology Report Generation with Key Phrase Extraction**
Park, Yoon, Kim, Choi. ACL 2026 Findings. arXiv:2504.07415.
- *Takeaway:* Retrieve *clinical key phrases* (not whole documents) conditioned on the image, reducing
  hallucination at far lower compute/data than a multimodal LLM.
- *Technical summary:* LLM extracts essential key phrases from reports; image-conditioned retrieval of
  phrases conditions the generator; SOTA CheXbert with competitive RadGraph F1 on MIMIC-CXR/IU X-Ray,
  and natural multi-view aggregation.
- *Plain-English:* Instead of pasting whole reports, it fetches the few clinically critical phrases an
  image warrants.
- *Applicability (A4, A2):* Fine-grained retrieval units ≈ A2 addressability: MORPHEUS could retrieve
  *pathway-level phrases/claims* rather than whole cards — a lighter frozen-trunk plug-in that beats a
  heavy encoder.
- *Novelty implication:* **Reframes** the retrieval unit granularity question; supports MORPHEUS's
  closed-RAG "whitelisted claims" design as a principled (not ad hoc) choice.

**21. AMANDA: Agentic Medical Knowledge Augmentation for Data-Efficient Medical VQA**
Wang, Mao, Wen, Luo, Ding. EMNLP 2025 Findings. arXiv:2510.02328.
- *Takeaway:* Training-free agentic framework: coarse-to-fine question decomposition (intrinsic) +
  biomedical *knowledge-graph* retrieval (extrinsic) improves zero/few-shot Med-VQA on 8 benchmarks.
- *Technical summary:* Decomposes the clinical question, then grounds each sub-question in KG retrieval;
  no fine-tuning, strong in data-scarce regimes.
- *Plain-English:* An agent breaks a hard medical question into parts and looks each part up in a
  medical knowledge graph.
- *Applicability (A1, A3):* A1's NL task decomposition/routing + A3's structured-knowledge grounding
  without retraining — a plausible outer loop for MORPHEUS's TQI that decomposes an NL task into
  slot/pathway sub-queries.
- *Novelty implication:* **Pre-empts** agentic KG-grounded medical QA; MORPHEUS's differentiator is
  grounding in a *learned tumor-state representation*, not just an external KG.

**22. Retrieval-Augmented Generation in Biomedicine: A Survey of Technologies, Datasets, and Clinical
Applications**
He, Zhang, Rouhizadeh, ... Teodoro. arXiv:2505.01146 (2025).
- *Takeaway:* Frames a "biomedical RAG trilemma" — reasoning depth vs. inference latency vs. data
  privacy — that constrains clinical deployment, and flags multimodal RAG misalignment as unsolved.
- *Technical summary:* Surveys 2020–2025 biomedical RAG across naive/advanced/modular paradigms; argues
  agentic workflows deepen reasoning but add latency, and cloud-vs-local is a privacy tradeoff; calls
  for self-correcting, verifiable clinical agents.
- *Plain-English:* A map of medical RAG that names the three-way tension every clinical system must
  navigate.
- *Applicability (A4, A3):* Positions MORPHEUS's frozen-trunk closed-RAG as a privacy/latency-favorable
  point in the trilemma; explicitly notes *multimodal* misalignment — the exact gap A4 targets.
- *Novelty implication:* **Reframes** A4 as addressing a named open problem (multimodal biomedical RAG
  misalignment), giving MORPHEUS a citable gap to claim novelty against.

---

### Part D — Retrieval-augmented biology (the encode-vs-retrieve tradeoff *for a biological modality*)

**23. BioBridge: Bridging Biomedical Foundation Models via Knowledge Graphs**
Wang, Wang, Srinivasan, Ioannidis, Rangwala, Anubhai. ICLR 2024. arXiv:2310.03320.
- *Takeaway:* Parameter-efficient KG "bridges" connect independent unimodal biomedical FMs (protein,
  molecule, text) *without fine-tuning them*, enabling cross-modal retrieval and extrapolation to
  unseen modalities/relations.
- *Technical summary:* Learns transformations between frozen unimodal embedding spaces guided by a
  biomedical KG (~76% avg. improvement over KG-embedding baselines); serves as a general-purpose
  cross-modal retriever and boosts multimodal QA and guided drug discovery.
- *Plain-English:* Rather than build one giant model for all of biology, it teaches small connectors so
  separate expert models can talk to each other via a knowledge graph.
- *Applicability (A4, A2):* The strongest biology-side template for A4's frozen-trunk plug-in: keep the
  WSI+RNA trunk frozen and *bridge* proteomics/phospho via a KG-guided connector + retrieval, rather
  than fusing an adapter. A2: cross-modal bridges are addressable per relation/pathway.
- *Novelty implication:* **Pre-empts** "connect frozen biomedical FMs without retraining" as novel;
  **strengthens** A4's frozen-trunk feasibility. MORPHEUS must differentiate via *identified/causal*
  slots (A2/A5) rather than KG-embedding bridges alone.

**24. Retrieval-based Controllable Molecule Generation (RetMol)**
Wang, Nie, Qiao, Xiao, Baraniuk, Anandkumar. ICLR 2023. arXiv:2208.11126.
- *Takeaway:* Steer a *frozen* generative model with a small set of retrieved exemplar molecules — no
  task-specific fine-tuning — and generalize beyond the retrieval database.
- *Technical summary:* Retrieval fuses exemplar + input molecules; self-supervised nearest-neighbor
  objective + iterative refinement dynamically updates outputs and the retrieval set; solves real tasks
  (e.g., SARS-CoV-2 main-protease binders) across multiple base generators.
- *Plain-English:* Give the model a few example molecules close to what you want, and it composes new
  ones — without retraining.
- *Applicability (A5, A4):* A5's "drug/perturbation as a query": a design criterion (exemplars) applied
  *after* encoding, over a frozen model — precisely MORPHEUS's intervention-as-query framing but in
  chemistry. A4: retrieval steering vs. fine-tuning.
- *Novelty implication:* **Pre-empts** "retrieval steering of a frozen generator as intervention" in
  the molecule domain; **reframes** A5 — MORPHEUS's novelty is doing this over an *identified tumor-
  state* with causal geometry, not a molecular generator.

**25. scRAG: Hybrid Retrieval-Augmented Generation for LLM-based Cross-Tissue Single-Cell Annotation**
Yu, Zheng, Chen, Hua, Luo. ACL 2025 Findings. aclanthology 2025.findings-acl.53.
- *Takeaway:* Cross-tissue cell-type annotation by retrieving *both* structured KG triples *and*
  similar cells from a reference database, then refining with marker genes — beating tissue-specific
  methods and trained classifiers.
- *Technical summary:* Two-stage hybrid retrieval (KG triples + reference-cell kNN) generates candidate
  types; marker genes from candidates and neighbors refine the LLM's prediction; strong cross-tissue
  generalization.
- *Plain-English:* To label a cell, it looks up both biological facts and the most similar known cells,
  then double-checks with marker genes.
- *Applicability (A4, A1, A3):* Direct biology instance of "retrieve a modality (reference cells)
  rather than encode it" for generalization to *unseen tissues* — MORPHEUS's cancer-held-out guardrail
  regime. Hybrid structured+unstructured retrieval mirrors A3 grounding + A1 NL routing.
- *Novelty implication:* **Pre-empts** retrieval-augmented cell annotation; MORPHEUS's edge is that its
  retrieval keys come from an *identified, causal* tumor-state (A2/A5), not raw expression similarity.

**26. GenePT: A Simple but Effective Foundation Model for Genes and Cells Built From ChatGPT
Embeddings** *(retrieval/embedding-vs-encoding baseline)*
Chen, Zou. bioRxiv 2023.10.16.562533 (2023/24); Genome Biology-track.
- *Takeaway:* Representing genes/cells by *retrieved* LLM text embeddings of gene descriptions (no
  transcriptomic pretraining) rivals or beats heavily-trained single-cell FMs (Geneformer/scGPT) on
  several tasks.
- *Technical summary:* Each gene → GPT text-embedding of its NCBI description; a cell → weighted average
  of its genes' text embeddings. Competitive on gene-property and cell-type tasks despite *no*
  parametric training on expression data — knowledge is retrieved from the LLM/text side, not encoded.
- *Plain-English:* Just describing each gene in words and embedding that description works about as well
  as training a giant model on millions of cells.
- *Applicability (A4):* A pointed encode-vs-retrieve datapoint *inside biology*: a retrieved/textual
  representation can match an expensively *encoded* one — evidence that some modalities are better
  represented via retrieved prior knowledge than learned encoders. Motivates a MORPHEUS ablation:
  encoded RNA trunk vs. text-embedding-retrieval baseline.
- *Novelty implication:* **Reframes / partially pre-empts** the premise that a modality must be encoded
  to be useful; a strong adversarial baseline A4's study must include. (Note: preprint — verify final
  venue before citing as peer-reviewed.)

---

### Cross-cutting synthesis for MORPHEUS

- **The encode-vs-retrieve tradeoff is well-established prior art (entries 1–13).** MORPHEUS *cannot*
  claim the mechanism, the frozen-trunk-plus-datastore pattern, adaptive when-to-retrieve, or
  "retrieval beats fine-tuning for new knowledge." These are settled. A4's defensible novelty must be a
  *modality-selection criterion specific to molecular biology* — a principled, testable rule for which
  of {proteomics, phospho, CNV, SNV} is encode-like vs. retrieve-like, ideally tied to the
  long-tail/identifiability arguments (entries 1, 5, 8) and *surviving distillation* (entry 9).
- **The multimodal-biomedical RAG surface is crowded but imaging/text-centric (entries 14–22).** No
  entry addresses retrieval over *structured molecular* modalities with no natural image/text form.
  That is the open lane the survey (22) explicitly names as "multimodal RAG misalignment."
- **Biology-side retrieval exists (23–26) but keys on raw similarity or external KGs, not an identified
  causal tumor-state.** MORPHEUS's differentiator is retrieving against A2-identified,
  A5-causal-geometric slots — a claim none of these pre-empt, but one that *depends* on A2/A5 actually
  delivering identifiability (unproven in the current build).
