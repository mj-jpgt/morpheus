# MORPHEUS Rebase — Literature Review (deduped, 5-axis synthesis)

*Master synthesis across 15 harvested lanes (`l01`–`l15`), organized by the five rebase axes.
Every mechanism MORPHEUS wants to claim is, individually, prior art; the defensible ground is the
**synthesis** and the **evaluation**. This review names the strongest papers per axis, the
pre-emptions to differentiate against, and the surviving whitespace.*

---

## 0. Corpus scope, provenance & quality caveats

**Volume.** ~500 lane-entries were harvested across 15 lanes (l01 multimodal-repr 34 · l02 pathology-FM 30 ·
l03 omics-FM 33 · l04 molecular-NL-prompting 37 · l05 promptable-unified 32 · l06 emergence-eval 37 ·
l07 AI-discovery 29 · l08 multimodal-RAG 26 · l09 alignment-identifiability 33 · l10 decision-support 44 ·
l11 benchmarks-confound 27 · l12 interventional-causal 34 · l13 agentic-science 32 · l14 compositional-instruction 35 ·
l15 steelman-prior-art 36). After cross-lane deduplication (papers such as Geneformer, scGPT, CPA, GEARS,
SurvPath, GenePT, Function Vectors, SENA-discrepancy-VAE, sVAE+ each recur in 3–6 lanes) the corpus is
**~300–350 distinct papers**.

**Quality.** All 15 lanes passed an adversarial referee pass (audit date 2026-07-29). Headline verdict across
the audits: **0 fabricated citations**; every spot-check (≈100 papers, deliberately weighted toward the
highest-risk 2026 arXiv IDs) resolved to a real, accurately-summarized work; all entries structurally complete.

**Caveats to carry forward:**
- **One withdrawn preprint** — ChatCell (arXiv:2402.08303, l04 #22): cite as existence-proof only.
- **One low-confidence entry** — "Self-driven biological discovery" (bioRxiv 2025.06.24.661378, l07 #29): no
  authors/DOI; argumentative role already covered by Coscientist/Robin/AI-co-scientist. Drop or backfill.
- **Metadata slips (papers real, details off):** LangCell author list (l03 #31), BulkFormer author parse
  (l03 #12), Symbolic-MoE→Skill-MoE rename (l05 #25), ConceptMix authors (l14 #33), GEARS off-by-a-year
  (Nat Biotech 2024 not 2023), SENA title abbreviates "…for Biological Data", several l10 date mislabels
  (L2M 2025 not 2024; SurvPGC 2026 not 2025).
- **Very recent 2026 preprints** (verified real but lightly cited): treat as provisional — e.g. Shesha
  coherence, MapPFN, Chem2Gen-Bench, CITE-VAE, scBIG (l12); An & Du Homomorphism-Error (l14); Winter et al.
  data-efficient alignment (l09/verdicts).
- **Remit stretches (defensible, flagged):** l04 Group V (SEQUOIA, FmH2ST, Stem, CellSymphony) are
  *pre-promptable fixed probes*, not prompting systems; l11 Section C drug/perturbation papers are
  *methodological analogues*, not WSI→molecular; some l11 venue/honor tags are lower-confidence.

---

## Axis A1 — Promptable unified representation + NL task auto-detection

**The paradigm is saturated prior art.** The text-prefix/instruction interface (T5, Raffel JMLR 2020),
in-context task inference (GPT-3, Brown NeurIPS 2020), the generalist-with-context (Gato, Reed TMLR 2022),
and instruction-tuned zero-shot-to-unseen-task-*families* (FLAN, Wei ICLR 2022; T0, Sanh ICLR 2022;
Super-NaturalInstructions, Wang EMNLP 2022) collectively own "one promptable model, NL selects the task."
Unified-IO 1/2 and OFA own "unified heterogeneous I/O via NL instructions." **MORPHEUS cannot claim the framing.**

**In-domain, the strongest neighbours are TOLD the task.** ChatNT (de Almeida, Nat Mach Intell 2025) —
one English model, 18+ genomics/transcriptomics/proteomics tasks over a frozen Nucleotide-Transformer, but
the task is *stated in the question* — is the single most important A1 comparator. Med-PaLM M (Tu 2023,
14 tasks incl. genomics), PathChat (Lu, Nature 2024 — the sharpest *oncology* pre-emption: promptable
histology+clinical-text copilot), CONCH/PLIP (CLIP-style NL-prototype pathology prompting, which the current
MORPHEUS "text_prototypes" stack already imitates), BiomedParse, Cell2Sentence/C2S-Scale, LangCell (zero-shot
NL cell-identity), and scMulan (attribute-conditioned) all deliver NL-conditioned biological inference — none
performs genuine **task auto-detection with abstention**, and none routes to *identified* slots.

**The mechanistic account MORPHEUS can borrow (bridges A1→A2).** Task auto-detection is implicit Bayesian
latent-concept inference (Xie, ICLR 2022); transformers in-context infer whole function *classes* (Garg,
NeurIPS 2022); the inferred task is a compact, causal, *composable* internal vector (Hendel task vectors,
EMNLP 2023; Todd function vectors, ICLR 2024). Ambiguity-as-multiple-valid-outputs is SAM's design
(Kirillov, ICCV 2023); instruction-conditioned readout is InstructBLIP's Q-Former (Dai, NeurIPS 2023).

**Emergence is contingent on data statistics, not scale** (Chan et al., NeurIPS 2022): in-context task
inference emerges only under burstiness + long-tailed/Zipfian classes. Task *breadth* is itself a scaling
axis (ZeroPrompt, EMNLP 2022). No omics FM conditions its promptability claim on these preconditions.

**Agentic systems already do NL→task routing externally** (Gorilla, ToolLLM, Biomni, TxAgent, CRISPR-GPT —
l13). Biomni is the sharpest "generalist biological task router" collision. Their weakness — long-horizon
brittleness (MLAgentBench ~37%, ScienceAgentBench ~32%) and search-cost blow-up at 16k tools — is MORPHEUS's
opening: *internalized*, short-horizon routing with shared weights across named slots.

**Strongest A1 papers to cite:** Xie 2022 (mechanism), FLAN/T0/Super-NaturalInstructions (held-out-task-family
protocol), ChatNT + PathChat + Med-PaLM M (in-domain pre-emptions to beat), Chan 2022 (emergence
preconditions), Hendel/Todd (task-as-addressable-vector), SAM + InstructBLIP (ambiguity + query-conditioned readout).

**A1 whitespace:** *task auto-detection + ambiguity-aware multi-output + abstention over biologically
identified slots*, validated on a held-out-task-family split — unoccupied.

---

## Axis A2 — Identified, pathway-addressable slots that make prompting reliable

**The identifiability price is fixed and known; name the inductive bias.** Locatello et al. (ICML 2019, best
paper) prove unsupervised disentanglement is impossible without inductive bias — the load-bearing guardrail.
The canonical recipe is a conditional prior p(z|u) (iVAE, Khemakhem AISTATS 2020) with metadata as the
auxiliary variable (Hyvärinen–Morioka 2016). Expect residual *linear* indeterminacy across runs and measure
it (Roeder, ICML 2021). Perturbation/sparse-mechanism-shift identifiers (Lachapelle CLeaR 2022 / JMLR 2024)
buy only **block/programme-group** identifiability under realistic soft perturbations — the honest ceiling.
Intervention-free fallbacks exist: known variable *grouping* (Morioka–Hyvärinen ICML 2024), mixture priors
(Kivva NeurIPS 2022), known causal *graph* (CauCA, Wendong NeurIPS 2023 — biology's pathways are causally
dependent, so drop the independence assumption).

**Bio-specific addressable-slot precedents largely pre-empt the architecture and disentanglement.**
SurvPath (Jaume, CVPR 2024) already tokenizes transcriptomics into **named pathway tokens** co-attended with
histology — *but defined a priori by gene sets, not identified, not promptable*. CPA (Lotfollahi, MSB 2023)
and biolord (Piran, Nat Biotech 2024) give per-attribute/perturbation addressable subspaces for counterfactuals.
**SENA-discrepancy-VAE** (de la Fuente, ICLR 2025) is the single sharpest A2 risk: identified + pathway-named +
causal + perturbation-predictive in one model. sVAE+ (Lopez 2023) and SAMS-VAE (Bereket–Karaletsos NeurIPS 2023)
own perturbation-conditioned identifiability; FCR (Mao 2024) adds identifiable context×treatment interaction.
Program-level is the right granularity (scBIG 2026; COMPASS 2025 couples a pathway bottleneck to ~76.5%
held-out-cancer transfer). The classical null is linear factor analysis (MOFA+, Argelaguet 2020) and logistic
regression (Boiarsky 2023) — identifiability must beat these *meaningfully*, not merely tie.

**Task/function vectors raise the bar** (Hendel, Todd): tasks are localized, steerable, additively composable
directions — so "the task lives in an addressable latent" is prior art; MORPHEUS must add pathway semantics +
identifiability guarantees. Move from *discovering* to *engineering* addressability ("Learning Task
Representations from ICL", 2025; instruction-conditioned readout, InstructBLIP).

**Operational addressability criteria (measurable, not asserted):** biological Homomorphism-Error (An & Du
2026 — R²=0.73 predicts OOD compositional accuracy); information-tightness spec (Yuanpeng Li 2025); the warning
that additive readouts can't express epistasis/synergy so slots need a non-additive/interventional readout
(Lippl–Stachenfeld, ICLR 2025); context-conditioned slots (conditional-attribute CZSL, Wang CVPR 2023);
VSA/HDC capacity ceilings on simultaneously-addressable items (Kleyko surveys 2023; Clarkson 2023).

**The confound lane *is* the A2 argument in disguise:** WSI→biomarker models predict a co-dependent bundle,
not an isolated biomarker (Dawood "Buyer Beware", Nat BME 2026); pathology FMs organize by medical center over
cancer type (de Jong Robustness Index, 2025; Mishra–Lotter RSA, 2025); dense fusion is argued the wrong
substrate (Tizhoosh 2025). These hand MORPHEUS concrete adopt-and-beat instruments (RI gain, co-dependence-leakage
reduction).

**Strongest A2 papers to cite:** Locatello 2019 (guardrail), iVAE (Khemakhem 2020), Lachapelle JMLR 2024
(block-identifiability ceiling), SENA-discrepancy-VAE (closest bio pre-emption), SurvPath (architecture
pre-emption), Dawood + de Jong (evaluation-form argument), An & Du 2026 (measurable addressability).

**A2 whitespace (MORPHEUS's defensible core):** the **identifiability→reliability link is unproven by anyone** —
no work shows identifiability (vs mere interpretability) *measurably improves prompt reliability and
cross-cohort/cross-cancer transfer*.

---

## Axis A3 — NL⇄biology grounding + emergent-knowledge elicitation and its evaluation

**"Emergence" is contested and easy to fake — the contribution is the metric, not the claim.** Apparent
emergence is often a metric mirage (Schaeffer, NeurIPS 2023: 25/29 metrics smooth out under continuous
scoring); aggregate loss is predictable while capabilities are not (Ganguli, FAccT 2022 — direct backing for
the thesis's "harness blind to representation quality"); most "emergence" is elicitation of latent knowledge,
not creation (Lu, ACL 2024). Index emergence by *loss* not parameters (Du, NeurIPS 2024) at a *per-skill quantum*
grain (Michaud, NeurIPS 2023), with mechanistic progress measures (Nanda, ICLR 2023) and awareness that
generalization may grok late (Power 2022).

**Reading latent knowledge (ELK) is open with mandatory sanity checks.** Unsupervised (CCS, Burns ICLR 2023),
linear-and-causal (Geometry-of-Truth, Marks–Tegmark COLM 2024), and supervised (SAPLMA, Azaria–Mitchell 2023)
latent-truth probes exist — **but** unsupervised probes provably latch onto the most *prominent* feature, not
knowledge (Farquhar, DeepMind 2023): a "latent-biology" probe on a tumor trunk could be reading a site direction.
ELK is genuinely unsolved (Christiano ARC 2021) — the closed-RAG hypothesis card may report the *whitelist /
human-simulator narrative*, not the trunk's latent state.

**Probing must be disciplined or it proves nothing:** control tasks + selectivity vs a random-programme
control (Hewitt–Liang EMNLP 2019); MDL / information-theoretic probing for probe-capacity robustness (Voita–Titov
2020; Pimentel 2020); **amnesic probing** — encoded ≠ causally used; erase (INLP) and measure behavior change
(Elazar TACL 2021), bridging A3→A5; password-locked controls to separate "trunk lacks knowledge" from "prompt
failed" (Greenblatt NeurIPS 2024); multi-format elicitation against sandbagging.

**The biology-specific nulls any A3 claim must clear:** GenePT (Chen–Zou, Nat Biomed Eng 2024) — GPT-3.5
text-embeddings of gene descriptions rival Geneformer/scGPT with *no* omics pretraining — the flagship A3
adversary; zero-shot single-cell FMs often lose to PCA (Kedzierska 2023), logistic regression ties them
(Boiarsky 2023), deep perturbation predictors don't beat a linear mean-shift (Ahlmann-Eltze 2024).
Emergence-from-scale is expected (C2S-Scale), so novelty must isolate emergence attributable to the
*identified/multimodal design*. The elicitation-*evaluation* methodology exists only for clinical **text**
(Med-PaLM multi-axis rubric, Singhal Nature 2023) — porting it to *mechanistic/pathway* knowledge is
unoccupied.

**Grounding is also a confound defense, but the "unified space" is fractured:** a language bottleneck
de-confounds (GLMP 2026); site/scanner/batch confounds pervade and survive normalization (Howard 2021; de Jong
RI 2025; Kömen 2024); grounding quality is measurable via alignment/uniformity (Wang–Isola 2020) and the
modality gap warns the biology-vs-language shared space is geometrically fractured by construction (Liang,
NeurIPS 2022).

**Strongest A3 papers to cite:** Schaeffer 2023 (metric mirage), Farquhar 2023 (distractor certificate),
Hewitt–Liang + Voita–Titov (probing discipline), Elazar 2021 (amnesic causal probing), GenePT + Kedzierska +
Ahlmann-Eltze (the nulls), Med-PaLM (elicitation rubric to port), Othello-GPT (Li ICLR 2023 — probe-then-intervene
gold standard).

**A3 whitespace (the clearest of all five):** no **confound-aware, validity-certified emergent-biological-knowledge
benchmark for a multimodal cancer FM** exists that separates elicited biology from text priors, scale, and site
confounds under probe-capacity-robust metrics.

---

## Axis A4 — Multimodal prompting: ENCODE vs RETRIEVE, and the frozen-trunk plug-in

**The tradeoff is settled prior art.** Frozen-trunk + external datastore beats parametric encoding, gains
concentrating on the rare/long-tail (kNN-LM, Khandelwal ICLR 2020; RETRO, Borgeaud ICML 2022); RAG (Lewis 2020)
and In-Context RALM (Ram TACL 2023 — the zero-surgery baseline any "encode" proposal must beat) define
modality-as-context; RAG beats unsupervised fine-tuning for new knowledge (Ovadia EMNLP 2024). **The mechanism
is not the contribution.**

**Long-tail/data-scarce modalities are exactly where retrieval wins** (NPM, Min ACL 2023 — retrieve *instead of*
encode; Atlas, Izacard JMLR 2023 — few-shot competence) — mapping onto CPTAC proteomics (inventory-only), rare
phospho-sites, rare mutation combinations. **But the advantage may be an under-trained encoder artifact** (Xu–Alon–Neubig
2023 — kNN gains can be distilled back): a genuine "must-retrieve" claim needs a gap that *survives distillation*.

**Frozen-trunk plug-in mechanisms are mature — pick, don't invent:** Q-Former soft-prompt tokens (BLIP-2, Li ICML
2023 — the concrete template for the `(batch, n_pathways, D)` slots the unwired TQI needs); soft-prompt-for-frozen-LM
(Frozen, Tsimpoukelli 2021); locked-encoder tuning (LiT 2022); gated cross-attention (Flamingo 2022);
tokenizer-only-per-modality (Meta-Transformer 2023); modality-expert FFNs (BEiT-3, ONE-PEACE). Subset-marginalizing
fusion for mosaic/missing modalities (PoE/MMVAE/MoPoE; MultiVI/totalVI/scGLUE) plus modality-dropout training
(Ma et al. CVPR 2022) handle the clinical reality that samples rarely carry every omic. Modality-specific molecular
encoders are mature (ESM-2, DeepPhospho, Phosformer, Prosit, PINNACLE) — the "encode" path can lean on frozen
existing encoders.

**The biomedical-multimodal-RAG surface is crowded but imaging/text-centric** (RA-CM3, MedRAG/MIRAGE, MMed-RAG,
RULE — calibrated retrieval count +47.4% factuality, HeteroRAG, FactMM-RAG). The survey (He 2025) explicitly names
**"multimodal RAG misalignment"** as unsolved. **No cited work retrieves over *structured molecular* modalities
with no natural text/image form.** Biology-side retrieval (scRAG, BioBridge, RetMol, GenePT) keys on raw similarity
or generic KGs, never an *identified causal tumor-state*. It is a spectrum — fuse ↔ attend-to-memory (Memorizing
Transformers) ↔ retrieve-text — not a binary.

**In-domain pre-emptions:** SurvPath (encode-everything a-priori pathway tokens), PORPOISE (encode-both default),
Med-Gemini (already blends encoded + retrieved but never asks *which molecular modality*); the Boehm et al.
(Nat Rev Cancer 2022) fusion taxonomy omits RAG-as-a-modality-choice entirely. GenePT is the live null: if a
text/RAG baseline matches an encoded latent, "encode" is unjustified for that modality.

**Strongest A4 papers to cite:** kNN-LM + NPM + Atlas (long-tail-favors-retrieval theory), Xu–Alon–Neubig
(survives-distillation bar), In-Context RALM (zero-surgery baseline), BLIP-2 (soft-prompt slots), MoPoE + Ma et al.
(mosaic/dropout), GenePT (the null), He 2025 survey (names the open gap).

**A4 whitespace:** a **pre-registered, quantitative molecular-modality-selection rule** (encode vs retrieve vs
marginalize per modality) with GenePT/In-Context-RALM/survives-distillation baselines built in — and retrieval
keyed by *identified* A2 slots over structured molecular memory. *(These A4 directions were surfaced by the
applicability synthesis but were not put through the adversarial novelty vote — see NOVELTY_LEDGER §3.)*

---

## Axis A5 — Interventional / causal-geometry queries

**Perturbation-as-a-query is thoroughly pre-empted — zero novelty for the capability itself.** scGen
(Lotfollahi Nat Methods 2019 — perturbation = post-encoding difference vector), CPA/chemCPA (additive composable
counterfactuals, incl. unseen drugs), GEARS (unseen combinatorial perturbations over a GO graph), biolord
(recombine attribute subspaces), STATE (Arc Institute, >100M perturbed cells), Geneformer in-silico deletion,
MapPFN (in-context zero-shot perturbation — closest "new question = new query, not a new model"), PertAdapt
(frozen-trunk + condition adapter beats GEARS).

**Identifiability-from-interventions is owned by the causal-representation literature** (Lachapelle
mechanism-sparsity; sVAE+; SAMS-VAE; soft-intervention identifiability, Zhang NeurIPS 2023; Ahuja/Squires/Varıcı
give an explicit **perturbation-coverage budget** — 1/node necessary+sufficient, 2/node assumption-light; FCR
identifies context×treatment). Claim only **block/pathway-group** identifiability (JMLR 2024). SENA-discrepancy-VAE
is the sharpest A5+A2 collision.

**Causal *geometry* is real but must be measured, not asserted:** directional *coherence*, not magnitude,
encodes regulatory architecture and predicts stress (Shesha S_p, Raju 2026); intervention-as-optimal-transport-map
(CellOT; scalable W1 solver) is geometry-native and now cheap; the pullback-metric geodesic machinery exists
(Chen ICANN 2019), with a Finslerian correction for stochastic decoders (Pouplin 2022) and a *train-the-latent-flat*
alternative (FlatVI 2025); manifold-fidelity is evaluable (TopOMetry, eLife 2024, which warns VAE latents distort
geometry).

**The field's evaluation is a live wound:** deep perturbation models frequently fail to beat a linear/mean-shift
baseline (Ahlmann-Eltze 2024; PerturBench; PertEval-scFM; and l11's Wong et al. — a CRISPR-informed mean beats
scGPT/GEARS). Drug-blind generalization is limited to *within-mechanism-of-action* class (Herbert 2026). Do **not**
assume drug = clean target knockdown — the chem→gene map is context-dependent (Chem2Gen-Bench 2026). The Virtual
Cell Challenge (Cell 2025) is the held-out-perturbation, cross-context protocol to adopt.

**Clinical counterfactuals ("treatment as a query") set the rigor bar:** BITES, counterfactual-survival-with-balanced-representations,
CFR-Net/TARNet (the IPM-balancing causal-geometry backbone); CINEMA-OT (confounder-separated counterfactual pairs
— the standard for calling a query "causal"); PDGrapher (the *inverse* query: which perturbation reaches a target
state); GPerturb + Conformal Risk Control (answers must carry sparse attributions + calibrated, abstaining UQ).

**Strongest A5 papers to cite:** scGen/CPA/GEARS/biolord/STATE (the pre-empted primitives), sVAE+ + soft-intervention
identifiability + JMLR-2024 ceiling (identification theory + honest bound), Shesha + TopOMetry (measure the geometry),
Ahlmann-Eltze + PerturBench + Virtual Cell Challenge (the null + protocol), PDGrapher (inverse query),
Conformal Risk Control (calibrated abstention).

**A5 whitespace:** every primitive is pre-empted; the only defensible novelty is the **synthesis on a substrate no
cited system uses** — a frozen, identified, pathway-named, NL-promptable, multimodal tumor trunk answering
forward *and* inverse counterfactuals whose causal geometry is *measured*, whose chem→gene translation is
*auditable*, and whose answers are *calibrated/abstaining*, benchmarked to beat the linear mean-shift null.
*(Dedicated A5 directions were surfaced by the applicability synthesis; only their A2-coupled forms entered the
adversarial vote — see NOVELTY_LEDGER §3.)*

---

## Cross-cutting verdict

1. **No single-axis mechanism is novel.** A1 (Med-PaLM M/PathChat/ChatNT), A2 (SurvPath/SENA/CPA/biolord),
   A4 (kNN-LM/RAG/GenePT), A5 (scGen/CPA/GEARS/Geneformer) are each thoroughly occupied.
2. **The defensible novelty is (a) the integrated substrate** — a *frozen, identified* oncology trunk that
   auto-routes NL tasks, exposes *identified* (not gene-set) pathway slots, accepts perturbation queries
   without retraining, evaluated confound-aware — **and (b) the evaluation/measurement contributions**, which
   are the hardest to scoop and directly repair the diagnosed "harness blind to representation quality."
3. **Two hard guardrails recur across lanes:** name the Locatello inductive bias explicitly (A2), and pair
   identified slots with a non-additive/interventional readout (Lippl). Two nulls recur: GenePT text-prior (A3/A4)
   and the linear mean-shift (A5). Any headline claim must clear them.
4. **Ranked whitespace (strongest→weakest):** A3 emergence/elicitation *evaluation* (genuinely unoccupied) →
   A2 identifiability→reliability *link* (unproven by anyone) → A4 molecular encode-vs-retrieve *decision rule*
   (named-open by the survey) → A1 auto-detection+abstention (sharp narrowed delta) → A5 measured causal geometry
   (highest-risk; must beat the linear null).
