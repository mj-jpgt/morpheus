# MORPHEUS Rebase — Novelty Ledger

*Each candidate direction with its adversarial 3-vote outcome and the prior art it must survive.
Verdict source: the fleet's adversarial novelty JSON (3 independent prior-art votes per candidate).
**Outcome across all 8 adjudicated candidates: SURVIVED (0/3 refuted votes each).** No candidate was
refuted — the adversarial lane found a *nearest-neighbour* for every one but no *scoop*. The honest
reading is therefore not "everything is clear" but "every candidate has a reframe-or-die near-collision;"
those collisions are named per row as the load-bearing differentiator.*

Legend — **Verdict:** SURVIVED (≥2/3 survive votes) / REFUTED. **Kill-risk:** the single closest prior
art that would refute a *naively-stated* version of the claim (the "CLOSEST"/nearest-neighbour). Axes: A1
promptable+auto-detect · A2 identified/addressable slots · A3 grounding+emergence-eval · A4 encode-vs-RAG ·
A5 interventional/causal-geometry.

---

## §1. Adjudicated candidates (all SURVIVED, 0/3 refuted)

### C1 — Auto-detecting & abstaining on under-specified biological queries · **A1** · **SURVIVED (0/3 refuted)**
**Claim.** A router that infers the latent biological task from an *ambiguous* NL query over a frozen
WSI+molecular tumor state, returns *multiple valid programme interpretations* when under-specified (SAM-style),
and *abstains* when no identified pathway slot addresses it — formalized as a posterior over identified
programmes (Xie implicit-Bayes) and evaluated on a held-out-task-family split.
**Prior art it must survive:** general LLM abstention (AbstentionBench 2025); ambiguous clinical text-to-SQL
two-stage interpretations (CLUES 2026); uncertainty-alignment abstention (SUA-TR 2026); abstention on a fixed
imaging task (stroke agentic AI 2026); task ambiguity in LMs (Tamkin 2022, adjacent); SAM ambiguity-as-multiple-masks
(adjacent); Tree-of-Clarifications / CLAM / Reasoning-about-Intent (open-domain multi-interpretation); ELISA/scPilot
(NL task routing over scFMs, no ambiguity/abstention); promptable biomedical multimodal models that are *told* the
task and do not abstain (Med-PaLM M, PathChat, ChatNT, LangCell).
**Kill-risk (CLOSEST):** *"Answer, Clarify, or Abstain: Fine-Grained Selective Prediction for Medical VLMs"*
(Nguyen et al., ICML 2026 workshop, non-archival). Shares the ANSWER/CLARIFY/ABSTAIN triage + a frozen-hidden-state
latent-mode auto-detector — **but** the latent variable is *input-quality condition*, not *biological task*;
radiology/CT/fundus, not WSI+molecular; no multi-interpretation; MoE+RankNet, not implicit-Bayes; cross-condition,
not held-out-task-family. **Differentiator that keeps it alive:** the latent being detected is the *biological
task/programme*, over a WSI+molecular tumor state, with multi-interpretation and held-out-task-family evaluation.

### C2 — Internalized prompt-conditioned routing to identified pathway slots (vs external agent glue) · **A1×A2** · **SURVIVED (0/3 refuted)**
**Claim.** Replace external tool-orchestration (Gorilla/Biomni/TxAgent) with routing *inside* one frozen trunk:
adapt Symbolic/Skill-MoE NL-skill routing so skill tags are named, identified pathway slots, materializing the
detected task as a Hendel/Todd function vector tied to a specific pathway, routed short-horizon (one prompt → one
addressed programme).
**Prior art it must survive:** Function Vectors (Todd 2023, generic); Skill/Symbolic-MoE (2503.05641, generic
skills, MedMCQA only); ProtoPathway (Reactome GNN, not FV routing); standard MoE over pathology+genomics;
SpatialFusion (pathway-activity fusion, no internalized routing); MoPE-MOI (pathway-guided sparse MoE over MSigDB
Hallmark experts); Task-Conditioned Routing Signatures (routing identifiability benchmark, generic LLM); TIGER
(NL-task-conditioned internalized routing in a vision FM); ICL-Router; GC-MoE (named cell-type expert routing, not
pathway slots, not frozen prompt-conditioned); Discovering Interpretable Biological Concepts in scRNA FMs (post-hoc
sparse dictionary probing, not internalized routing).
**Kill-risk (CLOSEST):** MoPE-MOI (named-MSigDB-Hallmark experts with patient-level traceability) + Skill-MoE
(NL-skill routing to named modules). **Differentiator:** routing is *internalized in one frozen prompt-conditioned
trunk*, the named modules are *identified* pathway slots (not a-priori gene-set experts), and the task materializes
as a *function vector tied to a pathway* — the union none of them hold.

### C3 — Declarative biological task-instruction benchmark + emergence-conditioned corpus · **A1×A3** · **SURVIVED (0/3 refuted)**
**Claim.** Ship (a) a broad NL biological task census (T0/ZeroPrompt-style: hundreds of tumor-state asks, many
paraphrases each, held-out task families) and (b) a corpus curated for Chan et al.'s distributional preconditions
(burstiness + long-tailed rare states/perturbations), then report auto-detection *as a function of* those measured
statistics and of held-out-task-family generalization.
**Prior art it must survive:** InstructCell (NL-instruction/omics copilot; no held-out task-family census, no
distributional conditioning); Cell-o1/CellPuzzles (held-out tissues, reasoning puzzles not declarative
auto-detection); SC-Arena (closest NL scFM benchmark, no held-out families or emergence-precondition conditioning);
BioLLM; Teddy; Stack (ICL without burstiness/long-tail conditioning); Biology-Instructions (multi-omics sequence
understanding); Chan et al. 2022 (the distributional-precondition theory being ported).
**Kill-risk (CLOSEST):** SC-Arena (NL single-cell reasoning benchmark) + Biology-Instructions (multi-omics
instruction benchmark). **Differentiator:** *held-out task families* + auto-detection reported *as a measured
function of burstiness/long-tail preconditions* — the emergence-conditioning that no biological benchmark applies.

### C4 — Identifiability-as-reliability: prove the missing empirical link · **A2** · **SURVIVED (0/3 refuted)**
**Claim.** Demonstrate that identified pathway slots (iVAE-style conditional prior + sparse-mechanism-shift on
perturbation data) yield measurably MORE reliable and transferable NL prompting than a-priori gene-set pathway
tokens (SurvPath) and unsupervised/linear factors (MOFA+/PCA), on leave-site-out and leave-cancer-out splits under
a confound-aware protocol; ship the reliability instrument (Dawood co-dependence-leakage reduction, de Jong
Robustness-Index gain, Roeder cross-run linear-identifiability stability).
**Prior art it must survive:** Lopez sVAE+ / Bereket SAMS-VAE (identifiability machinery, evaluate disentanglement/accuracy
not prompt reliability/transfer); Jiang "What Makes a Representation Good…" 2026 (identifiability for prediction
accuracy); PAMT (a-priori pathway tokens — the baseline class); Hou npj Digital Medicine (clinical-prompt integration,
no identifiability-vs-gene-set test); general-ML prompt-disentanglement for DG; SurvPath (nearest-neighbour, claims
accuracy/interpretability); iVAE / mechanism-sparsity (components only); FLARE / ProtoPathway (multi-site robustness /
interpretability, no identifiability→reliability claim).
**Kill-risk (CLOSEST):** SurvPath (pathway tokens, but a-priori/fixed and claims interpretability, *not*
identifiability-driven reliability/transfer). **Differentiator:** the causal test that *identifiability itself*
(not interpretability) buys prompt reliability/transfer — a link **no cited work demonstrates**. This is the least
pre-empted claim in the entire corpus.

### C5 — Frozen-trunk, partially-identified pathway slots with honest block-level guarantees · **A2×A4×A5** · **SURVIVED (0/3 refuted)**
**Claim.** Attach a lightweight adapter to a FROZEN multimodal trunk (WSI+RNA+optional proteomics/phospho) that
identifies pathway-GROUP (block) latent slots using perturbation sparse-mechanism-shift (sVAE+/SAMS) + a known
pathway-membership grouping (Morioka) + a known pathway causal graph (CauCA) as the named Locatello bias, projects
them onto named pathways (SENA-style), and exposes them for free-form NL addressing — claiming only
partial/programme-group identifiability against the JMLR-2024 ceiling.
**Prior art it must survive:** SAMS-VAE; sparse-mechanism-shift causal discovery (2206.02013); identifiability from
soft interventions (2307.06250); identifiability from purely observational data (JMLR 2024, 2410.22038); SENA layer
(Cell Systems 2025); Cradle-VAE; ProbVLM (probabilistic adapter for frozen VLMs); Language-Guided CBMs; SpatialFusion;
SurvConvMixer; sVAE+; SAMS-VAE; Disentangled Multimodal WSI+RNA (2508.16479, tumor/TME split only); PAST (multimodal FM).
**Kill-risk (CLOSEST):** Winter et al., *"Data-Efficient Multimodal Alignment for Histopathology-based Molecular
Prediction"* (arXiv:2606.29949, Jun 2026) — frozen WSI+RNA FMs + lightweight alignment adapter + open-vocabulary NL
gene-set/pathway prompting + an honest "graduated predictability spectrum." **Anticipates the system-level combination
but NOT the identifiability guarantee or latent-slot identification.** Also SENA-discrepancy-VAE (identified +
pathway-named + causal, but single-omic, non-frozen, no adapter/NL/block-identifiability). **Differentiator:**
the *identifiability guarantee* (block/programme-group, JMLR-2024-bounded) on identified *latent* slots — Winter has
the plumbing, not the guarantee; SENA has the guarantee, not the frozen-multimodal-NL-adapter substrate.

### C6 — Context-conditioned compositional slots with a biological-homomorphism objective + capacity budget · **A2×A5** · **SURVIVED (0/3 refuted)**
**Claim.** Represent slots as programme×cell-state conditional embeddings (CZSL) trained under a biological
Homomorphism-Error objective (adapt An & Du 2026) so slot arithmetic approximates programme composition, paired with
a non-additive/interventional readout (satisfying Lippl) and a VSA-derived capacity bound (Clarkson 2023) reporting
how many programmes can be co-addressed before crosstalk; deliver a ConceptMix-style combinatorial-load benchmark that
finds the binding ceiling.
**Prior art it must survive:** VSA capacity theory (Kleyko PIEEE 2022, general); compositional-generalization
diagnostics (Hupkes JAIR 2020); zero-shot scFM limitations (2025, critique only); interpretable cellular/gene-signature
embeddings (2021); "From modality-specific to compositional foundation models for cell" (Cell Systems 2026, perspective
only); brain-inspired HDC for DNA methylation (HDC-as-classifier only); An & Du Homomorphism-Error (the metric being
adapted, general-language); Clarkson VSA capacity (the theory being applied); Leveraging VL embeddings for zero-shot
histopathology (VLM embeddings, not CZSL slot arithmetic); "What Do Biological Foundation Models Compute? Sparse
Autoencoders" (interpretability, no VSA capacity/addressability metric).
**Kill-risk (CLOSEST):** the An & Du Homomorphism-Error metric + Clarkson VSA capacity theory (both general/non-biological)
and the Cell-Systems compositional-FM perspective (modular unimodal combination, no capacity bound). **Differentiator:**
first *biological* instantiation — HE objective on programme×cell-state slots + a *measured binding ceiling* on a
biological FM, which no cited work reports.

### C7 — BioELK-Bench: confound-aware emergent-biology elicitation benchmark with a validity certificate · **A3** · **SURVIVED (0/3 refuted)**
**Claim.** A curated probe battery measuring how much pathway/mechanism knowledge a frozen multimodal tumor trunk
encodes and how readily NL prompting elicits it, scored to survive four killer critiques: (i) MDL/MI + thresholded
accuracy under smooth AND hard metrics; (ii) Hewitt–Liang selectivity vs a random-programme control; (iii) a
Farquhar distractor/site sanity check requiring the recovered direction to FAIL at predicting submitting-site/scanner
while succeeding on biology; (iv) a decisive delta vs the GenePT text-prior and Kedzierska PCA baselines; (v) a
password-locked control separating "trunk lacks knowledge" from "prompt failed." Emergence plotted loss-indexed, at
per-pathway quantum grain, across multiple elicitation formats.
**Prior art it must survive:** Kedzierska/SCMBench SC-FM evaluations (no confound-certified elicitation battery);
"Deeper evaluation of a SC-FM" (NMI 2024); "Biology-driven insights…" (Genome Biology 2025); "Systematic Evaluation
of a SC-FM" (2602.17532); "Intermediate Layers Encode Optimal Biological Representations"; Hewitt–Liang / Voita–Titov /
Greenblatt / Elicitation-Game (borrowed method ingredients); Nature-2024 precision-oncology VLM (no elicitation
certificate); "Probing, Fusion, and Trustworthiness" (2606.17115 — closest "probe a multimodal cancer FM", lacks the
certificate); Exhaustive Circuit Mapping of a SC-FM (unimodal Geneformer); Pathology-FMs-are-Scanner-Sensitive /
fmMAP (confound ingredient only); GenePT (the text-prior null to beat); Knowledge Elicitation for cancer staging
(NL over clinical text); VCBench (closest integrated benchmark — single-cell not multimodal-cancer, not an
NL-elicitation/validity-certificate battery); Med-PaLM (elicitation rubric, clinical text only).
**Kill-risk (CLOSEST):** *"Probing, Fusion, and Trustworthiness: A Systematic Evaluation of Foundation Model
Representations for Multimodal Cancer Analysis"* (arXiv:2606.17115, 2026) — probes a multimodal cancer FM but **lacks
the confound certificate, selectivity, MDL, NL elicitation, text-prior baseline, password-lock, and emergence curves**.
And VCBench (pre-registered baselines + confound control, but single-cell, not NL-elicitation). **Differentiator:**
the *validity certificate* (recovered direction is biology, not the most-prominent site/distractor feature) assembled
for a *multimodal cancer* trunk — the unfilled union that ties A3 evaluation to A2 identifiability.

### C8 — Amnesic-counterfactual causal elicitation: elicited ≠ used, as a validity gate on hypothesis cards · **A3×A5** · **SURVIVED (0/3 refuted)**
**Claim.** Move A3 from correlational probing to causal elicitation: locate a pathway's latent direction, erase it
(INLP/amnesic) and measure whether the tumor-state/survival/drug-response prediction causally changes (Elazar), and
dually steer along it (RepE, Geometry-of-Truth) to confirm sufficiency. Use this as a validity gate on the closed-RAG
hypothesis card — the card's biological claims are accepted only if they match directions the trunk causally uses —
directly defusing the ELK "human-simulator" hazard where the card reports the whitelist rather than the latent state.
**Prior art it must survive:** Amnesic Probing (Elazar TACL 2021, NLP only); LEACE (Belrose NeurIPS 2023);
RepE (Zou 2023, LLM only); ELK human-simulator framing (Christiano 2021, alignment only); Geneformer in-silico
deletion (GENE space, not activation-space causal-use gate on a generated explanation); Multimodal prototyping for
survival (no causal-elicitation gate); causal dictionary learning in genomic LMs; XpertCausal (radiologist-guided
causal CBM for chest X-ray); CausalMixNet (medical imaging causal intervention); Decode-gLM (SAE interpret+steer of
genomic LMs, no generated-card validity gate); multimodal-FM-in-CRC trustworthiness review (aspirational, no gate);
biomedical-LM-interpretability scoping review.
**Kill-risk (CLOSEST):** Decode-gLM (SAE interpret+steer of Nucleotide-Transformer) + Elazar amnesic probing.
**Differentiator:** a *causal-use validity gate on a generated biological explanation* (the closed-RAG hypothesis
card) for a *multimodal tumor trunk* — converting the card from a fluency risk into a certified readout; no cited
work applies amnesic-counterfactual erasure as a gate on generated cards.

---

## §2. Ledger summary

| # | Candidate | Axes | Votes | Verdict | Closest kill-risk | Reframe-or-die differentiator |
|---|-----------|------|-------|---------|-------------------|-------------------------------|
| C1 | Auto-detect + abstain on ambiguous bio queries | A1 | 3/3 survive | **SURVIVED** | Nguyen "Answer/Clarify/Abstain" med-VLM | latent = *biological task*, WSI+molecular, multi-interp, held-out-task-family |
| C2 | Internalized routing to identified pathway slots | A1×A2 | 3/3 | **SURVIVED** | MoPE-MOI + Skill-MoE | internalized frozen trunk + *identified* slots + pathway function vector |
| C3 | Declarative bio task census + emergence-cond. corpus | A1×A3 | 3/3 | **SURVIVED** | SC-Arena + Biology-Instructions | held-out task families + burstiness/long-tail conditioning |
| C4 | **Identifiability-as-reliability (the missing link)** | A2 | 3/3 | **SURVIVED** | SurvPath | identifiability *buys* reliability/transfer — proven by no one |
| C5 | Frozen-trunk block-identified pathway slots | A2×A4×A5 | 3/3 | **SURVIVED** | Winter 2606.29949 + SENA | identifiability *guarantee* on identified latent slots |
| C6 | Context-cond. compositional slots + capacity budget | A2×A5 | 3/3 | **SURVIVED** | An&Du HE + Clarkson VSA | first biological HE objective + measured binding ceiling |
| C7 | **BioELK-Bench elicitation benchmark + certificate** | A3 | 3/3 | **SURVIVED** | Probing/Fusion 2606.17115 + VCBench | confound *validity certificate* on a multimodal cancer trunk |
| C8 | Amnesic causal-use validity gate on hypothesis cards | A3×A5 | 3/3 | **SURVIVED** | Decode-gLM + Elazar amnesic | causal-use gate on a *generated card*, multimodal tumor |

**No candidate was refuted.** The adversarial process found a nearest-neighbour for all eight but a *scoop* for
none — consistent with the corpus-wide finding that every mechanism is individually published while the
identified/multimodal/evaluated *synthesis* is not.

---

## §3. Un-adjudicated axis directions (surfaced but NOT put to the adversarial vote)

The applicability syntheses proposed additional directions on A4 and A5 that did **not** enter the 8-candidate
adversarial JSON. They are recorded here as *provisionally novel, pending a formal prior-art vote* — do not treat
them as vetted:

- **A4-D1 · Modality Encodability Score (MES):** a pre-registered, computable encode-vs-retrieve *decision rule*
  for {proteomics, phospho, CNV, SNV, bulk-RNA} with In-Context-RALM and GenePT baselines + a survives-distillation
  gate. *Status: strong whitespace (survey-named "multimodal RAG misalignment"), un-voted.*
- **A4-D2 · Pathway-slot-keyed molecular memory:** structured molecular modalities integrated as attention-readable
  memory keyed by *identified* A2 slots (fuse↔attend↔retrieve spectrum). *Status: bridges A2×A4, un-voted; depends
  on A2 delivering identifiability.*
- **A4-D3 · Mosaic frozen-trunk plug-in with adaptive encode/retrieve/marginalize gate.** *Status: un-voted.*
- **A5-D1 · Measured geodesic-coherence causal geometry** on the tumor manifold (Shesha coherence + TopOMetry
  fidelity). *Status: un-voted; highest-risk (must beat linear null).*
- **A5-D2 · Drug/perturbation as NL prompt over identified slots** with auditable chem→gene translation + calibrated
  abstention (forward + inverse/PDGrapher). *Status: partially covered by C5/C6, dedicated form un-voted.*
- **A5-D3 · "Identification-buys-intervention" benchmark:** show block-identified slots beat the linear mean-shift
  null and a matched non-identified trunk, reporting coverage against Squires(1/node)/Varıcı(2/node). *Status: the
  A5 analogue of C4; un-voted but high-value.*

**Recommendation:** route A4-D1 (MES) and A5-D3 (identification-buys-intervention) through the same 3-vote
adversarial process before committing — both are high-leverage and land in named-open gaps, but neither has been
formally screened for scoops.

---

## §4. Retired framings (dead as headline claims — do NOT re-claim)

These are *not* refuted candidates (none was submitted); they are the **bare axis framings** the corpus proves are
saturated. Retire them explicitly so no deliverable accidentally leads with them:

- **"An NL-promptable / task-inferring multimodal cancer model"** — dead (Med-PaLM M, PathChat, ChatNT, BiomedParse).
- **"Pathway-addressable multimodal tumor tokens"** — dead as *architecture* (SurvPath); survives only as *identified*
  (vs gene-set) slots whose identifiability is shown to buy reliability (→ C4).
- **"Perturbation/drug as a query on a learned latent"** — dead (scGen, CPA, GEARS, biolord, STATE, Geneformer);
  survives only on the frozen/identified/NL/multimodal substrate with *measured* geometry (→ C5, A5-D1/D3).
- **"Encode-vs-retrieve / frozen-trunk-plus-datastore mechanism"** — dead (kNN-LM, RETRO, RAG, In-Context RALM);
  survives only as a *molecular-modality decision rule* (→ A4-D1, un-voted).
- **"Our model has emergent biological knowledge"** — dead as a *claim* (Schaeffer mirage; Kedzierska/Boiarsky/
  Ahlmann-Eltze nulls); survives only as a *measurement method with a validity certificate* (→ C7, C8).
