# Axis A1 — Promptable unified representation + NL task auto-detection

**Axis remit:** one representation + a natural-language front-end that *infers* the requested
task and routes it, in place of hard-coded per-task probes/heads. The rebase question is not
"can you build an NL interface to a biological model" (answered: yes, many times) but "can a
single frozen trunk *auto-detect* an under-specified biological task and route it to an
identified programme, rather than being *told* the task via a fixed prefix/instruction."

---

## Key findings (with citations)

### F1. The promptable-unified paradigm is saturated prior art — MORPHEUS cannot own the framing
- **T5** (Raffel, JMLR 2020, arXiv:1910.10683) established task-prefix strings selecting behavior
  in one text-to-text model — the root ancestor of every promptable interface.
- **GPT-3** (Brown, NeurIPS 2020, arXiv:2005.14165) is the origin of "the prompt infers the task"
  (in-context task inference, no explicit task selection).
- **Gato** (Reed, TMLR 2022, arXiv:2205.06175): one weight set, 604 tasks across modalities, task
  inferred from context — canonical proof a shared token stream + context-conditioned decoding
  replaces task-specific heads.
- **OFA** (Wang, ICML 2022, arXiv:2202.03052) and **Unified-IO 1/2** (Lu, ICLR 2023 / CVPR 2024):
  NL instructions (not task IDs) route behavior through one seq2seq trunk over a shared vocabulary.
- **FLAN / T0 / Super-NaturalInstructions** (Wei, Sanh, Wang; ICLR/EMNLP 2022): instruction-tuning
  across many NL-verbalized tasks unlocks zero-shot generalization to *held-out task types*, and
  ablations show the NL phrasing itself (not just multitask exposure) drives it.
- **Implication:** A1's *paradigm* is table stakes. Novelty must live in the biological
  instantiation — specifically auto-detection + abstention + routing to *identified* slots (A2).

### F2. Task auto-detection has a mechanistic, *addressable* account MORPHEUS can borrow
- **In-context learning as implicit Bayesian inference** (Xie, ICLR 2022, arXiv:2111.02080):
  the model infers a latent *concept/task* from the prompt via a posterior over pretraining
  structure. This is the exact scaffold for framing "reliable prompting = identifiable
  latent-concept inference," bridging A1→A2.
- **What Can Transformers Learn In-Context** (Garg, NeurIPS 2022, arXiv:2208.01066): transformers
  in-context infer whole *function classes* — real algorithm learning, not memorization.
- **ICL Creates Task Vectors** (Hendel, EMNLP 2023 Findings, arXiv:2310.15916) and **Function
  Vectors** (Todd, ICLR 2024, arXiv:2310.15213): the inferred task lives in a compact, causal,
  *composable* internal vector that can be extracted and injected. Task and query are separable;
  the "task" is a manipulable slot.
- **Implication:** the detected task is a locatable, steerable vector — a direct precedent for
  exposing a biological "task/programme vector" as the routing mechanism, not a black-box head.

### F3. Auto-detection *emergence* is contingent on data statistics, not just scale
- **Data Distributional Properties Drive Emergent ICL** (Chan, NeurIPS 2022, arXiv:2205.05055):
  in-context (vs in-weight) task inference emerges *only* under burstiness + long-tailed/Zipfian
  class distributions; otherwise the model defaults to in-weight memorization.
- **ZeroPrompt** (Xu, EMNLP Findings 2022, arXiv:2201.06910): scaling the *number* of prompted
  tasks (~1000) improves zero-shot routing more than model size in the small-model regime — task
  *breadth* is a first-class scaling axis.
- **Implication:** for biological task auto-detection to *emerge*, the omics/WSI corpus must be
  curated for burstiness + long-tail (many rare cell states/perturbations) and for task breadth —
  a prescriptive, testable precondition no omics FM currently reports.

### F4. Nearly every in-domain system is TOLD the task — genuine auto-detection is the gap
- **ChatNT** (de Almeida, Nat. Mach. Intell. 2025, bioRxiv 2024.04.30.591835): one English-conversation
  model, 18+ genomics/transcriptomics/proteomics tasks, single frozen Nucleotide-Transformer
  encoder — but the **task is stated in the question**.
- **Med-PaLM M** (Tu, arXiv:2307.14334): generalist biomedical model, 14 tasks incl. genomics,
  task inferred from a prompt that **names the task**.
- **Text+Chem T5** (Christofidellis, ICML 2023, arXiv:2301.12586) and **BioT5/nach0/Mol-Instructions**:
  one-model-many-molecular-tasks, but signaled by an explicit **task prefix/instruction**.
- **LangCell** (Zhao, ICML 2024, arXiv:2405.06708): zero-shot NL cell-identity labeling — closest
  omics analog to promptable state, but a fixed classification framing.
- **PathChat** (Lu, Nature 2024, arXiv:2312.07814) and **CONCH/PLIP**: promptable *oncology* image
  models, but chat-wrappers / fixed NL-prototype classification over an encoder — no task inference
  beyond label prototypes, no identified latent, no counterfactual.
- **Implication:** the defensible A1 delta narrows to **task auto-detection with abstention over an
  ambiguous NL query on a multimodal WSI+molecular tumor state** — inferring *which* biological task
  is being asked (and declining when none is addressable), not executing a stated one.

### F5. Agentic systems already do NL→task routing — but *externally*, by gluing tools
- **Gorilla / ToolLLM / Biomni / TxAgent / CRISPR-GPT / BIA / AutoBA** (l13): "NL request → infer
  task → call the right frozen specialist as a tool." **Biomni** (bioRxiv 2025.05.30.656746) is the
  sharpest collision: one biomedical agent, many tasks, no per-task tuning.
- **ToolLLM** (Qin, ICLR 2024, arXiv:2307.16789): at 16k tools, retrieval + tree-search dominate —
  a caution that an open routing pool blows up in search cost.
- **MLAgentBench / ScienceAgentBench** (arXiv:2310.03302 / 2410.05080): long-horizon agent chains
  are brittle (best agents ~32–38% success), motivating short-horizon internalized routing.
- **Symbolic-MoE** (2025, arXiv:2503.05641): routes by *NL skill descriptions* to *named* expert
  modules — the tightest analog to A1+A2, but generic (not biologically identified).
- **Implication:** MORPHEUS's wedge is **internalizing** routing inside one frozen promptable trunk
  with shared weights across programmes (bounded, named, per-pathway slots), avoiding the
  search-cost blow-up and long-horizon fragility of external tool-glue — and giving
  identifiability/consistency guarantees a code-orchestrated toolbox cannot.

### F6. Ambiguity-aware promptable design is a template for under-specified biological queries
- **SAM** (Kirillov, ICCV 2023, arXiv:2304.02643): a promptable foundation model that handles
  ambiguity by emitting *multiple valid masks* rather than forcing one answer.
- **InstructBLIP** (Dai, NeurIPS 2023, arXiv:2305.06500): an instruction-aware Q-Former pulls
  *query-conditioned* features — the prompt should condition the *readout* of the representation,
  not post-hoc probe a fixed embedding.
- **Implication:** an under-specified biological query ("is this tumor immune-active?") should
  resolve to multiple valid programme interpretations, and the readout should be prompt-conditioned
  over addressable slots — not a single hard-coded probe.

### F7. The "hard-coded probe" paradigm A1 replaces — and its external evaluation harness
- **EAGLE** (arXiv:2502.13027): frozen two-FM pipeline, **43 hard-coded task heads** across 9 cancers
  — exactly the enumerate-a-head paradigm A1 must beat by *inferring* the task instead.
- **UNI / SEQUOIA / PEKA** (arXiv:2308.15474 / Nat. Commun. 2024 / arXiv:2504.07061): frozen trunk +
  linear probe / bespoke molecular regressor = the fixed-probe pattern.
- **HEST-bench** (Jaume, NeurIPS 2024 D&B, arXiv:2406.16192): the external harness — a promptable
  MORPHEUS trunk should be HEST-competitive *without task-specific retraining* to earn the A1 claim.
- **Held-out-task-family protocol** (FLAN / Super-NaturalInstructions): the credibility test for
  "genuine routing vs memorized heads" — hold out whole task clusters, not just held-out examples.

---

## What this implies for MORPHEUS design (build-on vs new-design)

**Build ON (established, do not re-claim):**
- The text-prefix / instruction interface (T5, OFA) and instruction-tuning recipe (FLAN, T0,
  ZeroPrompt) — including *many paraphrases per biological ask* for intent-based, prompt-robust
  routing (T0's multi-prompt lesson).
- Frozen-trunk + lightweight adapter / query-conditioned readout (Flamingo gated x-attn; MolCA /
  InstructBLIP Q-Former learned query tokens) as the mechanism to attach the NL front-end.
- Task-as-addressable-vector (Hendel task vectors, Todd function vectors) as the routing substrate.
- ReAct-style auditable routing trace (arXiv:2210.03629) so task detection is inspectable, not a
  black-box classifier.

**NEW-DESIGN (the defensible A1 deltas):**
1. **Auto-detection + abstention**, not stated-task execution. Every in-domain system (F4) is told
   the task; MORPHEUS infers it from an ambiguous NL query over a multimodal tumor state and
   abstains when no identified programme addresses it.
2. **Internalized routing** to *identified, named pathway slots* (F5 wedge → A2), replacing external
   tool-orchestration and hard-coded heads (EAGLE's 43 heads, F7).
3. **Emergence conditioned on data statistics** (F3): the auto-detection claim is tied to measured
   burstiness/long-tail + task breadth of the corpus, and validated by held-out-task-family
   generalization (F7) — a discipline no omics FM (scGPT uses per-task heads; ChatNT is told the
   task) currently applies.

---

## Candidate research directions

### D1. Auto-detecting and abstaining on under-specified biological queries
**Claim:** Build a router that infers the latent biological task from an *ambiguous* NL query over
a frozen WSI+molecular tumor state, returns *multiple valid programme interpretations* when the
query is under-specified (SAM-style), and *abstains* when no identified pathway slot addresses it —
formalized as a posterior over identified programmes (Xie's implicit-Bayes latent-concept
inference), evaluated with a held-out-task-family split (FLAN / Super-NaturalInstructions).
**Why novel:** Every in-domain promptable system (ChatNT, Med-PaLM M, Text+Chem T5, LangCell,
PathChat) is *told* the task via an explicit instruction/prefix/question and none abstains; combining
auto-detection + ambiguity-aware multi-output + abstention over *biologically identified* slots is
unoccupied. It converts A1 from "an NL interface" into a measurable inference-under-ambiguity claim.

### D2. Internalized prompt-conditioned routing to identified pathway slots (vs external agent glue)
**Claim:** Replace external tool-orchestration (Gorilla/Biomni/TxAgent) with routing *inside* one
frozen trunk: adapt Symbolic-MoE's NL-skill routing so the skill tags are *named, identified*
biological programme slots, with the detected task materialized as a Hendel/Todd-style function
vector tied to a specific pathway. Route short-horizon (one prompt → one addressed programme),
avoiding ToolLLM search-cost blow-up and MLAgentBench long-horizon fragility.
**Why novel:** Prior routing is either generic (Symbolic-MoE, no biological grounding) or external
tool-glue with no shared representation and no identifiability (Biomni). Internalized routing with
shared weights across *identified* pathway slots — benchmarked on routing consistency/identifiability
against an agent-orchestrated toolbox — is the un-pre-empted intersection of A1 and A2, and directly
answers the reviewer's "how is this not Biomni / just another agent" objection.

### D3. A declarative biological task-instruction benchmark + emergence-conditioned corpus
**Claim:** Ship (a) a broad NL *biological task census* (T0/ZeroPrompt-style: hundreds of distinct
tumor-state asks, many paraphrases each, organized into held-out task families) and (b) a
pretraining-corpus curation targeting Chan et al.'s distributional preconditions (burstiness +
long-tailed rare cell states/perturbations), then report auto-detection performance *as a function
of* those measured statistics and held-out-task-family generalization.
**Why novel:** There is no biology equivalent of Super-NaturalInstructions with held-out task
families, and no omics FM conditions its promptability claim on measured burstiness/long-tail
preconditions for emergent task inference. It makes "task auto-detection emerges in biology" a
falsifiable, pre-registered claim rather than an assertion — and gives the field the missing harness
that distinguishes genuine routing from memorized per-task heads (EAGLE, scGPT).
