# Axis A2 — Identified, pathway-addressable slots that make prompting reliable

**Axis definition.** A2 is the claim that MORPHEUS exposes *identified* latent slots (recoverable up to a benign transform across runs) that are *per-programme addressable* (each slot corresponds to a named biological programme/pathway), so that an NL prompt lands on the same biology every time. It is the reliability substrate under A1 (routing is only reliable if slots are stable) and the structural precondition for A5 (you can only "perturb programme X" if X is an addressable slot).

This file synthesizes the A2-relevant findings across all 15 harvested lanes (`l01`, `l03`, `l05`, `l09`, `l11`, `l12`, `l14`, `l15` carry the load) into: key findings, what they imply for build-on vs new-design, and concrete candidate directions.

---

## Key findings (with citations)

### The identifiability price is fixed and known — MORPHEUS must *name* its inductive bias, not claim the primitive

- **Unsupervised disentanglement is provably impossible without inductive bias.** Locatello et al. (ICML 2019, best paper) prove infinitely many entangled generative models fit any factorized latent; a 12,000-model study finds no unsupervised metric reliably picks disentangled models (`l09` #8). This is the **load-bearing guardrail for A2**: MORPHEUS cannot claim identified pathway slots from unsupervised scRNA/WSI alone. Every addressability claim must name its bias (pathway priors, perturbations, multimodal pairing, conditioning).
- **The canonical recipe that buys identifiability is a conditional prior.** iVAE (Khemakhem, Kingma, Monti & Hyvärinen, AISTATS 2020) proves deep latent-variable models are identifiable up to permutation + pointwise transform *iff* the prior is factorized and conditioned on an auxiliary variable `p(z|u)` (`l09` #3). Metadata (batch, tissue, patient, timepoint) is the **identifiability driver**, not a nuisance to normalize away (Hyvärinen & Morioka 2016, `l09` #1).
- **Expect a residual *linear* indeterminacy, and measure it.** Roeder, Metz & Kingma (ICML 2021) show softmax/contrastive models are identifiable only up to a single linear transform across runs/architectures (`l09` #4). Implication: slot addressability is a learned linear map; "linear identifiability" is the correct *cross-run stability metric* to report, and MORPHEUS should quotient the indeterminacy, not claim exact recovery.
- **Interventions and sparse mechanism shift are the biology-native identifiers.** Lachapelle et al. (CLeaR 2022; JMLR 2024) prove factors are identifiable when mechanisms depend *sparsely* on actions — but realistic soft/sparse perturbations buy only **block/programme-group identifiability, not per-dimension** (`l12` #31/#32, `l09` #9). This sets the *honest ceiling* for what MORPHEUS may claim.
- **Intervention-free fallbacks exist.** Morioka & Hyvärinen (ICML 2024) show a *known grouping* of observed variables gives causal-representation identifiability with **no** temporal structure, interventions, or supervision (`l09` #28). Gene-module / pathway-membership groupings are exactly this bias — a valuable hedge when perturbation data is absent. Kivva et al. (NeurIPS 2022) similarly get identifiability from a *mixture* prior (cell-state clusters) with no labels (`l09` #6); CellPLM already uses a Gaussian-mixture latent prior (`l03` #7).
- **Biology's programmes are causally *dependent*, so drop the independence assumption.** CauCA (Wendong et al., NeurIPS 2023) generalizes ICA to causally-linked latents given a *known causal graph* and interventions (`l09` #27). Pathways interact — MORPHEUS should adopt "known-graph causal programmes," a biology-honest reframe away from "independent slots."

### Bio-specific addressable-slot precedents already exist — the architecture and disentanglement are largely pre-empted

- **Pathway-addressable multimodal tokens already exist architecturally.** SurvPath (Jaume et al., CVPR 2024) tokenizes transcriptomics into **named biological-pathway tokens** and co-attends them with histology (`l15` #31, `l12` context). Critically, its pathway tokens are **defined a priori by gene sets — not *identified*** (iVAE-style) and not promptable/counterfactual. This is the sharpest **architecture pre-emption** of A2, and it pinpoints MORPHEUS's only opening: *identified* slots vs. fixed gene-set tokens.
- **Per-attribute addressable subspaces are published.** biolord (Piran et al., Nat Biotech 2024) gives each known attribute (type, time, perturbation, dose) its own latent subspace you can set independently for counterfactuals (`l12` #5, `l15` #30) — the closest published analog to per-programme addressable slots. CPA (Lotfollahi et al., MSB 2023) disentangles basal + additive perturbation + covariate slots (`l03` #17, `l12` #3, `l15` #26). Both address slots by **perturbation/metadata label**, not by *pathway programme*.
- **Identified latent dims ↔ named pathways is nearly solved by one paper.** SENA-discrepancy-VAE (de la Fuente et al., ICLR 2025) makes each latent *causal* factor an interpretable combination of learned **pathway-activity scores**, at no cost to predictive accuracy (`l12` #10). This is **the single biggest novelty risk for A2/A3** — the "identified + pathway-named + causal + perturbation-predictive" combination in one model.
- **Perturbation-conditioned identifiability is the mechanism MORPHEUS rests on.** sVAE+ (Lopez et al., 2023) treats perturbations as sparse unknown-target interventions and proves identifiable causal latents whose dims map to perturbation-affected mechanisms (`l12` #8, `l15` #35). SAMS-VAE (Bereket & Karaletsos, NeurIPS 2023) learns a *sparse mask* over which latent slots each perturbation touches (`l12` #7). FCR (Mao et al., 2024) adds identifiable **context×treatment interaction** blocks (`l12` #9).
- **Program-level is the right granularity.** scBIG (2026) shows programme/module-level representations beat gene-level for unseen/combinatorial perturbation (`l12` #16), corroborating a `(batch, n_programme, D)` slot exposure. COMPASS (2025) couples a pathway/concept bottleneck with **held-out-cancer transfer** at ~76.5% (`l11` #23).
- **The classical baseline is linear factor analysis.** MOFA+ (Argelaguet et al., 2020) already yields sparse, interpretable, pathway-annotatable factors from multi-omics (`l03` #20). MORPHEUS must show its slots capture *nonlinear, prompt-conditioned, causal* structure MOFA+ cannot — and, per Boiarsky et al. (2023), that identifiability yields **measurable transfer/reliability gains over logistic regression / PCA**, not just competitive accuracy (`l15` #22).

### Task/function vectors show slots are addressable and composable — but this raises the bar

- **The inferred task lives in an addressable, composable vector.** ICL "task vectors" (Hendel et al., EMNLP 2023) and "function vectors" (Todd et al., ICLR 2024) show tasks are localized to steerable mid-layer directions that **compose by addition** (`l05` #19/#20, `l14` #21/#22). This supports A2 feasibility *and* pre-empts "the task lives in an addressable latent" — MORPHEUS must add **biological-pathway semantics + identifiability guarantees**, not just a steerable vector.
- **Move from *discovering* to *engineering* addressability.** "Learning Task Representations from ICL" (2025) treats task representations as a *trained objective* rather than a post-hoc extraction (`l05` #21). InstructBLIP (Dai et al., 2023) makes the readout **instruction-conditioned** so the prompt steers which features are pooled (`l14` #6) — query-conditioned addressing is the mechanism that makes prompting reliable.
- **Reliable prompting = identifiable latent-concept inference.** Xie et al. (ICLR 2022) formalize task auto-detection as implicit Bayesian inference over a latent concept (`l05` #16). If biological programmes are the latent concepts, **identifiability determines whether the posterior can address them reliably** — the tightest theoretical bridge from A1 into A2.

### Operational criteria for addressability — measurable, not asserted

- **Homomorphism error is a trainable addressability metric.** An & Du (2026) show a *Homomorphism Error* (misalignment between symbol-space composition and hidden-state composition) predicts OOD compositional accuracy (R²=0.73) and improves it when minimized (`l14` #20). MORPHEUS can adapt this to "**biological homomorphism**": combining slot activations should be homomorphic to combining programmes.
- **Information-tightness is the identifiability spec.** Li (2025) proves compositional generalization needs components that encode *exactly* one factor's information — no more, no less (`l14` #12). A concrete design target for slots.
- **Identified slots + additive readout is insufficient for non-additive biology.** Lippl & Stachenfeld (ICLR 2025) prove kernel/linear readouts on perfectly compositional representations can only sum per-component values, failing on epistasis/synergy (`l14` #11). A2 must be paired with an **expressive/interventional readout** (A5), or it over-claims.
- **Slots must be context-conditioned.** Conditional-attribute CZSL (Wang et al., CVPR 2023) shows an attribute's meaning depends on its object (`l14` #30) — a pathway's "activation" means something different per cell state, so slots should be `programme × cell-state` conditional embeddings.
- **Addressability has a capacity ceiling.** VSA/HDC theory (Kleyko et al. 2023; Clarkson et al. 2023) gives hard bounds on representation width vs. number of simultaneously bound/addressable items before crosstalk (`l14` #26–29). MORPHEUS should *quantify* how many programmes it can reliably address at once.

### Evaluation lane: the confound literature *is* the A2 argument in disguise

- **"Buyer Beware" is the identifiability argument in evaluation form.** Dawood et al. (Nat BME 2026) show WSI→biomarker models predict a **co-dependent bundle**, not an isolated biomarker (`l11` #8). Without disentangled, per-programme addressable slots, a biomarker *prompt* returns a confounded aggregate. This directly motivates A2 and hands MORPHEUS a test: **do identified slots reduce co-dependence leakage?**
- **Center/site is encoded more strongly than biology.** de Jong et al. (2025) Robustness Index finds 9/10 pathology FMs organize by medical center over cancer type (`l11` #9); RSA studies concur (`l11` #11). A promptable representation with RI<1 routes prompts on *site*, not biology — so **A2 identifiability should be validated by an RI improvement** and by modeling site as an addressable, factored-out slot (HiST calibration token, `l11` #4).
- **Dense fusion is argued to be the wrong substrate.** Tizhoosh (2025) contends single dense embeddings cannot represent tissue's combinatorial richness (`l11` #25) — an authority-level argument *for* MORPHEUS's move to identified, addressable, compositional slots.

---

## What this implies for MORPHEUS design (build-on vs new-design)

**Build on (established machinery — cite, don't claim):**
- **Conditional-prior identification** (iVAE) with **metadata as the auxiliary `u`** — the recipe for "slots that mean the same thing every run." Not a linear/orthogonality penalty.
- **Sparse-mechanism-shift / perturbation-conditioned identification** (sVAE+, SAMS-VAE, Lachapelle) as the *identifying signal* — but claim only **block/pathway-group** identifiability (JMLR 2024 ceiling).
- **Pathway projection of latent factors** (SENA) and **pathway-token architecture** (SurvPath) for the naming/addressing layer.
- **Query-conditioned readout** (InstructBLIP) and **slot substrates** (Perceiver latents, attention bottlenecks, `l01` #3/#14) as the addressing mechanism.
- **Intervention-free biases** (known pathway grouping — Morioka; mixture prior — Kivva/CellPLM; known causal graph — CauCA) for the label-scarce/clinical regime.

**New design (the genuinely open ground):**
1. **The identifiability→reliability link is unproven by anyone.** No cited work shows that *identifiability* (vs. mere interpretability) *measurably improves prompt reliability and cross-cohort/cross-cancer transfer*. SurvPath's tokens are named-but-fixed; SENA's are identified-but-not-prompted; biolord/CPA address *labels* not *pathway programmes*. This link is MORPHEUS's defensible core.
2. **Frozen-trunk exposure of identified, NL-addressable slots.** Every bio precedent trains a bespoke model. Doing perturbation-conditioned/grouping-based identification via *light adapters on a frozen multimodal trunk*, then exposing the slots for free-form NL prompting, is unoccupied.
3. **An addressability *measurement* protocol.** Adopt-and-adapt: biological Homomorphism Error (An & Du), co-dependence-leakage reduction (Dawood), Robustness Index gain (de Jong), and a VSA-style capacity budget — turning "our slots are addressable" from an assertion into a reported number.

**Two hard guardrails MORPHEUS must foreground:** (a) name the Locatello inductive bias explicitly; (b) pair identified slots with a non-additive/interventional readout (Lippl), or the compositional claim collapses.

---

## Candidate research directions

### Direction 1 — Identifiability-as-reliability: the missing empirical link
**Claim.** Prove that *identified* pathway slots (conditioned via iVAE-style priors + sparse-mechanism-shift on perturbation data) produce measurably **more reliable and transferable NL prompting** than (i) a-priori gene-set pathway tokens (SurvPath) and (ii) unsupervised/linear factors (MOFA+, PCA) — on leave-site-out and leave-cancer-out splits, under a confound-aware protocol. Ship the reliability instrument: reduction in Dawood-style biomarker co-dependence leakage, Robustness-Index gain (de Jong), and cross-run linear-identifiability stability (Roeder).
**Why novel.** Every prior (SurvPath, SENA, biolord, CPA, COMPASS) claims *interpretability* or *accuracy*; **none demonstrates that identifiability itself buys prompt reliability/transfer**. Boiarsky and the confound lane show the field lacks exactly this measurement, making it both the open gap and the credibility bar.

### Direction 2 — Frozen-trunk, partially-identified pathway slots with honest block-level guarantees
**Claim.** Attach a lightweight adapter to a frozen multimodal (WSI + RNA + optional proteomics/phospho) trunk that, using perturbation data as sparse mechanism-shift signal (sVAE+/SAMS) *plus* a known pathway-membership grouping (Morioka) *plus* a known pathway causal graph (CauCA) as the named Locatello bias, identifies **pathway-group (block) latent slots** projected onto named pathways (SENA-style) and exposed for NL addressing. Claim only *partial/programme-group* identifiability (JMLR 2024), reported against the block-identifiability ceiling.
**Why novel.** sVAE+/SAMS/SENA each do a *piece* on a bespoke, single-omic, non-frozen model; **no cited system combines frozen-trunk + adapter-based identification + intervention-free grouping fallback + pathway naming + NL addressability**, and none states the honest block-level guarantee rather than over-claiming per-dimension disentanglement.

### Direction 3 — Context-conditioned compositional slots with a biological-homomorphism objective and capacity budget
**Claim.** Represent slots as `programme × cell-state` conditional embeddings (CZSL conditional attributes), trained under a **biological Homomorphism-Error objective** (adapt An & Du 2026) so that slot arithmetic ≈ programme composition, paired with a non-additive/interventional readout (to satisfy Lippl) and a VSA-derived capacity bound (Clarkson 2023) that reports how many programmes can be co-addressed before crosstalk. Deliver a ConceptMix-style combinatorial-load benchmark that finds MORPHEUS's binding ceiling.
**Why novel.** Moves A2 from *discovered* task/function vectors (Hendel, Todd) and *asserted* pathway tokens to **engineered, measured, capacity-bounded addressability** — supplying the operational addressability metric and capacity theory that no biological FM currently reports, and directly answering the compositional-generalization cautions (Lippl, Hupkes) the field raises.
