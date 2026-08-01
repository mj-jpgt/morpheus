# MORPHEUS Rebase — Recovered Vision & Five-Axis Formalization

*The up-front conceptual frame that targets the literature search. The fleet's findings will
sharpen, challenge, and re-rank these axes; treat this as a hypothesis, not a conclusion.*

---

## 1. The founding thesis (recovered from PRISM_ROADMAP / prismarchitecture / dev-plan)

A single **identified, compositional, causally-geometric tumor-state representation** such that
*"every clinical question becomes a query on that representation, and new questions require only a new
query, not a new model."* The intended surface was a **Task Query Interface (TQI)**: natural-language
task specs resolved against a **frozen** representation, extensible to new **modalities without
retraining the trunk** (proteomics as the marquee test), and capable of **counterfactual /
interventional** answers (intervention simulator; geodesic causal similarity).

Original novelty properties (declared necessary, not optional): **identification** (perturbation-
conditioned iVAE → latent dims = real pathway states), **compositional fidelity** (VSA binding → rare
combinations lossless), **causal geometry** (metric tensor → geodesic ≈ causal distance),
**intervention simulation** (drug = intervention spec applied *after* encoding), **modality-agnostic
pathway slots**, and **NL querying (TQI)**.

## 2. The drift (from the code inventory + root notes)

The live `v2/TumorStateV2` is a solid **unified typed representation** (WSI+RNA aligned; anchorable
identity/biology/context/patient heads; clean `DenseAdapter`+availability-mask fusion; leakage-safe
targets; safety-gated closed-RAG cards) — **but the promptable/identified/causal ambition is gone**:

- **No NL prompt to the model; no task auto-detection.** Every task is a hard-coded frozen-state probe
  (`comprehensive_evaluation.py:_rows_for`). "Prompting" = fixed CLIP label prototypes + numeric
  soft-kNN "molecular prompting" + a closed-RAG Qwen renderer that only *reorders whitelisted claims*.
- **No per-pathway slot structure.** PRISM's `TaskQueryInterface` (`PRISM/models/tqi.py`) is an
  **unwired scaffold** needing `(batch, n_pathways, D)` slots the model never exposes.
- **Proteomics / phospho / CPTAC are inventory-only** (a `NaN` eval placeholder); SNV/CNV/clinical are
  optional fused adapters; the encoder core is really just WSI+RNA.
- Root notes retracted the single-fused-latent and demoted identification / causal-geometry /
  "disentangled" — recovering C-index at the cost of the original novelty.
- **This session's finding reinforces the stakes:** the biology head collapses to low effective rank,
  a covariance term reverses it (+53, ~2.1×, 3 seeds), yet the confounded benchmark score is
  *unchanged* — i.e. the current evaluation is structurally blind to representation quality. A
  promptable/identified representation needs evaluation the current harness cannot provide.

## 3. The five axes (search themes = candidate-paper seeds)

Each axis: recovered-vision root → build-on (exists) vs. new-design (missing) → why it could be novel.

**Axis 1 — Promptable unified representation + NL task auto-detection.** *The core delta.* A single
latent with an NL interface that **infers** the requested task and routes it, vs. hard-coded probes.
Build on: `QueryBlock` typed-slot attention; `text_prototypes` embedding infra; PRISM `tqi.py`
(ScopeDetector + head registry) as reference. New: task-text→query conditioning + a real router +
scope/abstention. *Novelty question the fleet must answer:* has anyone built a **task-inferring,
NL-promptable multimodal cancer representation** (vs. per-task probes or a chat wrapper)?

**Axis 2 — Identified, pathway-addressable slots that make prompting meaningful.** Prompting is only
reliable if the latent is **identifiable and addressable** (per-programme/pathway). Build on: biology
head + programme regression; perturbation-conditioned-iVAE identification thesis. New: expose
`(batch, n_programme, D)` slots; recover identification as the load-bearing claim. *Novelty question:*
does identifiability/addressability measurably improve prompt reliability & transfer?

**Axis 3 — NL⇄biology grounding, emergent knowledge & its evaluation.** Move between prompt-language
and biological "language"; **elicit and rigorously measure emergent/latent biological knowledge** —
the hard, novel part is the *elicitation + emergence evaluation*, not generation. Build on:
hypothesis-card + closed-RAG. New: flexible NL output head + an emergence/elicitation benchmark.
*Novelty question:* is there a principled way to measure emergent biological knowledge in a
multimodal cancer model that isn't just downstream accuracy?

**Axis 4 — Multimodal prompting: encode vs. retrieve.** A principled account (and system) of when a
modality (proteomics, phospho, CNV, SNV) is best an **encoded input** vs. **RAG context** — the user's
own open question, and the marquee frozen-trunk generalization test. Build on: `DenseAdapter`+mask
extensibility; CPTAC inventory. New: proteomics/phospho adapters + a RAG-vs-encode study. *Novelty
question:* is "which modalities to encode vs. retrieve" an open, formalizable problem in multimodal
biomedical FMs?

**Axis 5 — Interventional / causal-geometry queries.** Counterfactual "what if we perturb X / give
drug D" as a **query**, not a retrained classifier. Build on: PRISM intervention-simulator +
Riemannian-metric design; DepMap/Perturb-seq (ORBIT Axis 5). New: wire it to the promptable interface.
*Novelty question:* can interventional queries be answered from a frozen, promptable representation
better than correlational baselines, and is that framing novel?

## 4. Guardrails to keep (honest engineering discipline, not reasons to abandon the vision)
Cancer-held-out validation (not within-seen-cancer early stopping); per-sample missingness as a
first-class training signal; environment-balanced IRM (or drop it); earn-your-place modality gating;
confound-aware evaluation (this session's rank-decoupling result shows why headline metrics mislead).

## 5. What "success" looks like
≥3 concrete, novelty-vetted paper directions, each mapped to one or more axes and to real MORPHEUS
assets, with an explicit build-on/new-design split and the experiments required. The fleet must be
willing to **kill** an axis as non-novel (the adversarial prior-art lane) — the goal is truth, not
confirmation.
