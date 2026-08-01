# Axis A5 — Interventional / causal-geometry queries: applicability synthesis

**Axis:** A5 — counterfactual perturbation / drug as a *query* on a frozen identified representation,
not a retrained classifier; causal-geometry (geodesic / metric-tensor / optimal-transport) queries on
a tumor-state manifold.

**Method:** read across all 15 harvested literature lanes; the A5-load-bearing material concentrates in
`l12_interventional_causal`, `l15_steelman_prior_art`, `l09_alignment_identifiability`,
`l10_decision_support`, and `l07_ai_discovery`. Findings below are the ones that most constrain or
enable A5, with the design verdict (build-on vs new-design) and three concrete candidate directions.

---

## 1. Key findings (with citations)

### 1a. "Perturbation-as-a-query on a learned latent" is thoroughly pre-empted — MORPHEUS gets zero novelty for the capability itself
- **scGen** (Lotfollahi, Wolf, Theis, *Nat. Methods* 2019) established the exact A5 primitive: a
  perturbation is a *difference vector applied after encoding*; add it to control cells to predict
  responses in unseen cell types. Post-encoding intervention arithmetic is 2019 prior art.
- **CPA / chemCPA** (Lotfollahi et al., *MSB* 2023; Hetzel et al., *NeurIPS* 2022) — encode to a
  disentangled *basal* state, then compose additive perturbation/dose/covariate embeddings for OOD
  counterfactuals; chemCPA reads molecular structure so the query generalizes to structurally novel
  drugs. This *is* the "drug = intervention spec applied after encoding" contract.
- **GEARS** (Roohani, Huang, Leskovec, *Nat. Biotech.* 2023) — GNN over a GO/co-expression graph
  predicts unseen single- and multi-gene perturbations (additive-embedding conditioning = a query).
- **biolord** (Piran et al., *Nat. Biotech.* 2024) — per-attribute latent subspaces recombined to
  generate unseen cell-type × perturbation × dose counterfactuals.
- **STATE** (Arc Institute, 2025) — industrial-scale (>100M perturbed cells, 70 contexts) State-Embedding
  trunk + State-Transition head taking (state, perturbation) as an explicit query pair.
- **Geneformer** (Theodoris et al., *Nature* 2023) — *in-silico* gene deletion on a frozen embedding
  already IS "perturbation as a query, not a retrained classifier," for discovery.
- **MapPFN** (Sextro et al., 2026) — prior-fitted network doing *in-context, zero-shot* perturbation
  prediction: a new context is a prompt, not a retrain — the closest architectural analog to
  "new question = new query, not a new model."
- **PertAdapt** (bioRxiv 2025) — a lightweight condition-sensitive adapter makes a *frozen* scFM trunk
  beat GEARS, validating the frozen-trunk + plug-in-head pattern for interventions.

### 1b. Identifiability-from-interventions (the theory MORPHEUS leans on for "slots = real pathways") is owned by the causal-representation literature; only *partial/block* identifiability is honest
- **Mechanism-sparsity nonlinear ICA** (Lachapelle et al., *CLeaR* 2022; *JMLR* 2024) — sparse actions
  identify latents up to permutation; the 2024 nonparametric result says realistic sparse/soft
  perturbations buy **block-level (programme-group)**, not per-dimension, identifiability. This bounds
  what MORPHEUS can honestly claim.
- **sVAE+** (Lopez et al., *CLeaR* 2023) and **SAMS-VAE** (Bereket & Karaletsos, *NeurIPS* 2023) — single-cell
  sparse-mechanism-shift models where each perturbation edits a sparse, inferable subset of latents;
  this is MORPHEUS's "which pathway slots does this query touch" mechanism, already published.
- **Soft-intervention identifiability** (Zhang et al., *NeurIPS* 2023, discrepancy-VAE base) — latent
  causal variables + graph are identifiable from *soft* (non-atomic) interventions, the regime CRISPRi/drug
  perturbations actually live in; also predicts unseen perturbation *combinations*.
- **Interventional CRL** (Ahuja et al., *ICML* 2023), **Linear causal disentanglement via interventions**
  (Squires et al., *ICML* 2023 — one perfect intervention per node necessary+sufficient), **General
  identifiability** (Varıcı et al., *AISTATS* 2024 — two uncoupled interventions/node, assumption-light):
  these give an explicit **perturbation-coverage budget** MORPHEUS must report against.
- **FCR** (Mao et al., 2024) — proven block-/component-wise identifiability of covariate, treatment, and
  covariate×treatment *interaction* blocks: the "same drug, different tumor state → different answer"
  factor MORPHEUS needs, already identified.
- **SENA-discrepancy-VAE** (de la Fuente et al., *ICLR* 2025) — **the single sharpest A5/A2 novelty risk**:
  identified + pathway-named + causal + perturbation-predictive in one model, at no accuracy cost.

### 1c. Causal *geometry* is real and biologically meaningful — but it must be **measured, not asserted**
- **Shesha directional coherence S_p** (Raju, 2026) — across 2,200+ perturbations, the *direction*
  (not magnitude) of response vectors encodes regulatory architecture and predicts UPR/stress
  (p<10⁻¹⁸); pleiotropic regulators pay a "geometric tax" of large-but-incoherent shifts. Direct
  evidence that interventional geometry carries causal signal, and that *coherence/direction*, not
  geodesic length alone, is the loaded quantity.
- **CellOT** (Bunne et al., *Nat. Methods* 2023) and the **scalable W1 neural-OT solver** (Uscidda,
  Cuturi et al., 2025) — a perturbation is an optimal-transport map between control and treated
  distributions; geometry-native and now scalable, so "intervention = geometric map" is not novel.
- **Pullback-metric geodesics** (Chen et al., *ICANN* 2019) — the canonical decoder-Jacobian Riemannian
  machinery MORPHEUS's "metric tensor → geodesic ≈ causal distance" would instantiate; establishes both
  correctness and the query-time compute cost.
- **Finslerian latent geometry** (Pouplin, Hauberg et al., 2022/23) — for *stochastic* decoders (tumor-state
  decoders are stochastic) the correct metric is Finslerian, not Riemannian.
- **FlatVI / Euclidean-latent VAE** (Palma, Theis et al., 2025) — the design fork: *train the latent flat*
  so straight lines = geodesics (cheap) vs *compute* geodesics at query time (expensive).
- **TopOMetry** (Sidarta-Oliveira et al., *eLife* 2024) — tooling to *evaluate* how faithfully an embedding
  preserves manifold geometry; warns default VAE latents distort it. Closes the "eval blind to geometry" gap.

### 1d. The field's evaluation is a live wound: deep perturbation models frequently fail to beat linear/mean-shift baselines
- **Ahlmann-Eltze, Huber, Anders** (bioRxiv 2024) — SOTA perturbation predictors (GEARS-class, FM-based)
  do **not** beat a simple additive/linear mean-shift on held-out perturbations. The A5 null hypothesis.
- **PerturBench** (GSK.ai, 2024) and **PertEval-scFM** (2024) — model rankings rank-shuffle across
  metrics/splits; frozen scFM embeddings give little advantage over mean baselines under a fixed probe,
  gains concentrated in easy (weak-effect) genes.
- **Virtual Cell Challenge** (Bunne, Roohani et al., *Cell* 2025) — community "Turing test": held-out
  perturbation, cross-context evaluation. The protocol to adopt.
- **Chem2Gen-Bench** (Lin & Chen, 2026) — chemical↔genetic perturbation equivalence is "measurable but
  heterogeneous"; **do not assume drug = clean target knockdown** — the chem→gene map is context-dependent.

### 1e. Clinical counterfactuals ("treatment as a query") already exist and set the rigor bar
- **BITES** (Schrod et al., *Bioinformatics* 2022), **counterfactual survival with balanced
  representations** (Chapfuwa et al., *CHIL* 2021), **CFR-Net/TARNet** (Shalit et al., *ICML* 2017):
  individualized treatment-effect / counterfactual-survival is established oncology prior art; the
  balanced-representation (IPM) trick is the causal-geometry backbone.
- **CINEMA-OT** (Dong et al., *Nat. Methods* 2023) — OT-matched counterfactual cell pairs with explicit
  confounder separation: the bar for calling any A5 query "causal" rather than correlational.
- **AFAPE** (von Kleist et al., *JMLR* 2025) — how to *retrospectively* evaluate next-best-test /
  next-best-assay counterfactuals from logged data without deploying; **Learning-to-Measure** (Kobayashi
  et al., 2024) does in-context AFA (a pre-emption risk for "acquisition as a query").
- **PDGrapher** (Gonzalez, Zitnik et al., *Nat. Biomed. Eng.* 2025) — the *inverse* query: "which
  perturbation combination moves this diseased state to a target state?" — MORPHEUS should support both
  forward and inverse modes.
- **GPerturb** (Wu et al., *Nat. Commun.* 2025) and **Conformal Risk Control** (Angelopoulos et al.,
  *ICLR* 2024) — interventional answers should carry sparse attributions + calibrated, distribution-free
  uncertainty and be able to abstain.

---

## 2. What this implies for MORPHEUS design (build-on vs new-design)

**Build-on (do not reinvent):**
- **Intervention-after-encoding contract** (scGen/CPA/biolord) — keep the "drug/perturbation = spec applied
  *after* a frozen encode" interface; **replace additive linear composition** with identified, possibly
  geodesic, latent moves.
- **Frozen-trunk + causal plug-in head** — CITRIS (normalizing-flow causal head on a frozen encoder,
  *l09*) and PertAdapt (condition-sensitive adapter, *l12*) show the frozen-trunk + intervention-head
  pattern works; adopt it rather than fine-tuning the trunk.
- **Geometry machinery** — pullback-metric geodesics (Chen 2019), FlatVI (train-flat), OT maps
  (CellOT / W1 solver), Finsler correction for stochastic decoders. Pick a design; don't invent the tooling.
- **Identification objective** — sparse-mechanism-shift (Lachapelle / sVAE+ / SAMS-VAE) under soft
  interventions (Zhang 2023); **claim only block/pathway-group identifiability** (JMLR 2024).
- **Evaluation** — Virtual Cell Challenge held-out/cross-context protocol, PerturBench distributional
  metrics (energy distance / MMD), PertEval frozen-probe methodology, AFAPE offline causal estimators,
  conformal risk control + calibration auditing.

**New-design (where the only defensible A5 novelty survives):** every A5 *primitive* is pre-empted, so
novelty is the **synthesis on a substrate no cited system uses** — a *frozen, identified, pathway-named,
NL-promptable, multimodal tumor (WSI + multi-omic) trunk* answering *forward and inverse* counterfactuals
whose **causal geometry is measured (coherence + geodesic fidelity), whose chem→gene translation is
auditable, and whose answers are calibrated and abstaining**. Three concrete instantiations follow.

---

## 3. Candidate research directions

### D1 — Measured geodesic-coherence causal geometry on a frozen multimodal tumor manifold
**Claim:** On MORPHEUS's frozen multimodal tumor-state latent, define a *learned metric tensor whose
geodesic length approximates interventional (causal) distance*, and evaluate the "geodesic ≈ causal"
hypothesis directly — using **directional coherence** (Shesha S_p) of predicted response vectors as the
primary geometric signal and **TopOMetry-style manifold-fidelity scoring** as the yardstick — rather than
asserting Euclidean latent arithmetic is causal.
**Why novel:** scGen/CPA use *Euclidean* latent arithmetic; CellOT/W1 give OT geometry; FlatVI flattens;
pullback/Finsler give the metric machinery — but **no cited paper (a) ties geodesic length to causal
distance, (b) uses coherence/direction rather than magnitude as the loaded quantity, and (c) validates
geometric fidelity on a *multimodal clinical (WSI+omics) tumor* manifold instead of scRNA.** Shesha shows
direction beats magnitude biologically; TopOMetry shows VAE latents distort geometry — so a *measured*
causal-geometry claim on a tumor trunk is genuinely unoccupied.

### D2 — Drug/perturbation as an NL prompt over identified pathway slots, with auditable chem→gene translation and calibrated abstention
**Claim:** Expose forward ("effect of drug/knockdown X on this tumor state") *and* inverse ("which
perturbation moves this tumor toward remission", à la PDGrapher) counterfactuals as **natural-language
queries routed to sparse, block-identified pathway slots**, where the drug→target step is an explicit,
*evaluated* chem→gene translation (not an identity assumption) and every answer carries sparse
attributions + conformal/GPerturb-style uncertainty and can abstain.
**Why novel:** SENA-discrepancy-VAE already gives identified pathway-named causal factors; CPA/GEARS/STATE
give the prediction; PDGrapher gives inverse queries; Chem2Gen-Bench warns drug≠clean knockdown;
GPerturb/conformal give UQ — but **no cited system combines NL task-inference, pathway-slot addressability,
an auditable chem→gene map, and calibrated abstention into a single frozen-trunk query interface.** The
delta over SENA/CPA is the *NL-promptable, hedged, forward+inverse* interface; the delta over MapPFN/STATE
is *pathway-named identified slots + multimodal tumor grounding* instead of numeric prompts on scRNA.

### D3 — An "identification-buys-intervention" benchmark: prove block-identified pathway slots beat the linear mean-shift null
**Claim:** Build a held-out-perturbation, cross-context evaluation (Virtual Cell Challenge + PerturBench
metrics) that **isolates the causal contribution of *identifiability itself*** — showing MORPHEUS's
block-identified, pathway-addressable slots beat (i) a linear/additive mean-shift baseline (the
Ahlmann-Eltze null), (ii) a matched *non-identified* frozen trunk under an identical probe (PertEval
methodology), while (iii) **reporting perturbation coverage against the Squires (1/node) and Varıcı
(2/node) identifiability budgets** so the "identified" claim is falsifiable, not decorative.
**Why novel:** PerturBench/PertEval/Ahlmann-Eltze establish that current models — including frozen scFMs —
*fail* to beat linear baselines, but **none tests whether an explicitly *identified* representation buys
interventional accuracy over a non-identified one, nor ties the claim to intervention-coverage theory.**
This turns MORPHEUS's identifiability thesis from an assumption into a measured, adversarially-baselined
result — and directly repairs the project's own diagnosed "harness blind to representation quality" gap.
