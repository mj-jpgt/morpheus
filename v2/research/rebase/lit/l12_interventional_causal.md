## Interventional / causal representation & perturbation modeling

Lane l12 — the **A5 core**: counterfactual/interventional prediction, perturbation-response
modeling (Perturb-seq, GEARS, scGPT-perturb, DepMap), causal representation learning &
identifiability, drug-response from representations, and Riemannian/geodesic latent geometry.
Every entry is mapped to MORPHEUS's five rebase axes (A1 promptable NL task routing; A2 identified,
pathway-addressable slots; A3 NL⇄biology grounding & emergent-knowledge eval; A4 multimodal
prompting / encode-vs-RAG / frozen trunk; A5 interventional/causal-geometry queries).

Convention for **Novelty implication**: *strengthens* = supports a MORPHEUS claim as feasible/valuable;
*pre-empts* = someone already did (part of) the claimed novelty, so MORPHEUS must differentiate;
*reframes* = shifts how the claim should be stated or evaluated.

---

### 1. GEARS: Predicting transcriptional outcomes of novel multigene perturbations
**Authors/venue:** Roohani, Huang, Leskovec — *Nature Biotechnology* 2023. https://www.nature.com/articles/s41587-023-01905-6 (bioRxiv 2022.07.12.499735)
**Takeaway:** A GNN over a gene–gene knowledge graph predicts Perturb-seq responses to gene combinations that were *never* individually perturbed.
**Technical summary:** GEARS represents each gene and each perturbation with embeddings tied together by a Gene Ontology / co-expression knowledge graph, giving a graph inductive bias that lets a GNN extrapolate to unseen single and multi-gene perturbations. On two Perturb-seq screens (1,543 RPE-1 and 1,092 K562 perturbations) it predicted four genetic-interaction subtypes with ~40% higher precision than prior methods and found the strongest interactions ~2× better. The perturbation is injected as an additive embedding conditioning the decoder — a query, not a retrained classifier.
**Plain-English:** Feed it a combination of genes to knock out and it forecasts how the whole transcriptome will shift, even for combinations no one has ever tested, by leaning on a map of which genes are related.
**Applicability:** A5 (the canonical "perturbation as a query on a representation" baseline MORPHEUS's intervention simulator must beat or subsume); A2 (its knowledge-graph structure is an external prior — MORPHEUS's claim is to have *internal* identified pathway slots instead). Design implication: GEARS is the SOTA MORPHEUS's counterfactual head is measured against; its additive-embedding conditioning is the minimal interface to exceed.
**Novelty implication:** *pre-empts* the bare "predict perturbation response" claim — MORPHEUS cannot claim novelty for counterfactual transcriptomics per se, only for doing it on a *frozen, identified, NL-promptable* multimodal tumor representation.

### 2. scGPT: toward building a foundation model for single-cell multi-omics
**Authors/venue:** Cui, Wang, Maan, Pang, Luo, Wang — *Nature Methods* 2024. https://doi.org/10.1038/s41592-024-02201-0
**Takeaway:** A 33M-cell generative transformer fine-tunable to perturbation prediction, but only marginally beating GEARS.
**Technical summary:** scGPT is a masked-generative transformer over gene tokens pretrained on >33M cells; for perturbation it adds a perturbation token to the perturbed gene and predicts the post-perturbation profile. Downstream benchmarks show it matches or slightly exceeds GEARS after fine-tuning, and its zero-shot in-silico perturbation *lags* GEARS. It demonstrates the "one pretrained trunk, many task heads" pattern MORPHEUS wants.
**Plain-English:** A GPT-style model for cells; you can adapt it to predict perturbation effects, but on that specific task it barely edges out a much smaller specialized model.
**Applicability:** A1 (foundation-model trunk + task heads, the promptable-representation template); A4 (frozen-ish trunk + plug-in adaptation); A5 (in-silico perturbation as fine-tuned head). Design implication: evidence that a generic pretrained trunk does *not* automatically win at intervention — MORPHEUS must justify identification/geometry as the differentiator, not scale alone.
**Novelty implication:** *reframes* — foundation-model scale is not sufficient for causal/perturbation gains; strengthens MORPHEUS's bet that *structure* (identifiability, geometry) is the load-bearing novelty.

### 3. CPA: Compositional Perturbation Autoencoder
**Authors/venue:** Lotfollahi, Klein, … Theis — *Molecular Systems Biology* 2023 (bioRxiv 2021.04.14.439903). https://www.biorxiv.org/content/10.1101/2021.04.14.439903v1.full
**Takeaway:** An adversarial autoencoder that additively composes a "basal" cell state with drug/dose/covariate embeddings to answer counterfactuals.
**Technical summary:** CPA uses an adversarial classifier to scrub perturbation information out of the basal latent, so perturbation, dose, and covariate embeddings can be *added back* at inference to generate out-of-distribution counterfactuals (unseen drug combinations, dose-response curves) with uncertainty estimates. This "encode to a disentangled basal state, then apply intervention as additive latent shift" is precisely the MORPHEUS "drug = intervention spec applied *after* encoding" pattern. Learned drug embeddings are interpretable and support dose extrapolation.
**Plain-English:** It strips a cell's latent code of any drug signature, so you can then paste in a drug's effect — including combinations and doses never measured — to predict the treated cell.
**Applicability:** A5 (the archetype of intervention-after-encoding and additive-latent counterfactuals); A2 (adversarial disentanglement of basal vs. perturbation — a weak form of the addressable-slot idea). Design implication: MORPHEUS's intervention simulator should adopt CPA's compose-after-encode contract but replace additive linear composition with identified, possibly geodesic, latent moves.
**Novelty implication:** *pre-empts* additive counterfactual composition; MORPHEUS must show that *identified/geodesic* composition beats CPA's additive latent arithmetic and works from a frozen tumor trunk.

### 4. chemCPA: Predicting cellular responses to novel drug perturbations at single-cell resolution
**Authors/venue:** Hetzel, Böhm, Kilbertus, Günnemann, Lotfollahi, Theis — *NeurIPS* 2022. https://arxiv.org/abs/2204.13545
**Takeaway:** Extends CPA with a molecular (chemical-structure) encoder so it generalizes to *structurally novel* drugs, not just held-out combinations of seen drugs.
**Technical summary:** chemCPA replaces CPA's per-drug lookup embedding with a learned function of molecular structure (fingerprints / GNN), and transfers from bulk (L1000) to single-cell via architecture-surgery fine-tuning. This lets counterfactual prediction extend to drugs absent from training, a genuine OOD generalization over CPA. It is the drug-structure→latent-intervention bridge.
**Plain-English:** Instead of memorizing each drug, it reads the molecule's chemistry, so it can guess the effect of a brand-new compound.
**Applicability:** A5 (OOD drug counterfactuals); A4 (chemical structure as an *encoded* modality feeding an intervention spec — informs the encode-vs-RAG decision for drug identity). Design implication: for MORPHEUS, a drug's chemistry is a candidate *encoded* intervention conditioner (like chemCPA) rather than RAG context; supports treating "the drug" as a structured query.
**Novelty implication:** *pre-empts* structure-conditioned drug counterfactuals; MORPHEUS differentiates by grounding the *response* in an identified tumor-state manifold rather than a generic basal AE.

### 5. Biolord: Disentanglement of single-cell data
**Authors/venue:** Piran, Cohen, … Nitzan — *Nature Biotechnology* 2024. https://www.nature.com/articles/s41587-023-02079-x
**Takeaway:** Encodes each known attribute (cell type, time, perturbation, dose) into its *own* latent subspace plus one unknown-attribute block, then recombines subspaces to generate unseen states.
**Technical summary:** biolord enforces disentanglement not adversarially but by construction: separate encodings per labeled attribute define a decomposed latent that the generator recombines, so "virtually shifting" a cell across attribute values produces experimentally inaccessible counterfactuals. It outperforms SOTA on unseen-drug and unseen-genetic-perturbation prediction. Attribute subspaces are directly addressable — the closest published analog to MORPHEUS's "per-programme addressable slots."
**Plain-English:** It gives each biological factor (type, age, drug) its own private set of latent dials, so you can turn one dial to a value never seen and read out the resulting cell.
**Applicability:** A2 (per-attribute addressable subspaces — the strongest existing instance of addressability); A5 (recombination = counterfactual query); A1 (attribute selection is a primitive form of task routing). Design implication: MORPHEUS should treat biolord as the reference for *addressable slots*, extending from known metadata attributes to *identified pathway programmes*.
**Novelty implication:** *pre-empts* attribute-addressable latent slots most directly — MORPHEUS's A2 claim must be sharpened to *pathway/programme* addressability with identifiability guarantees, and to slots that a *frozen trunk* exposes for NL prompting, which biolord does not do.

### 6. CellOT: Learning single-cell perturbation responses using neural optimal transport
**Authors/venue:** Bunne, … Krause, Rätsch, Alvarez-Melis, Cuturi — *Nature Methods* 2023. https://www.nature.com/articles/s41592-023-01969-x
**Takeaway:** Models a perturbation as an *optimal-transport map* between unpaired control and treated cell distributions, capturing heterogeneous subpopulation responses.
**Technical summary:** CellOT parameterizes each perturbation's response as an OT map via a pair of dual potentials (input-convex neural nets), learning couplings between control and perturbed states without pairing. It outperforms prior methods on scRNA-seq and multiplexed protein imaging, and generalizes to held-out lupus patients under IFN-β. OT gives a *geometry-native* notion of "how a cell moves under intervention."
**Plain-English:** It learns the cheapest way to morph a population of untreated cells into treated ones, respecting that different cells react differently.
**Applicability:** A5 (intervention as a transport map — a geometric counterfactual); geometry link to MORPHEUS's geodesic causal-distance ambition. Design implication: OT maps are a concrete instantiation of "intervention = movement on a manifold"; MORPHEUS's geodesic causal-similarity should be compared to / built on OT displacement geometry.
**Novelty implication:** *reframes* MORPHEUS's causal-geometry claim as OT-adjacent; strengthens feasibility of geometric interventions but means "intervention as geometric map" is not itself novel — the novelty is a *learned metric tensor whose geodesics equal causal distance* on a frozen tumor latent.

### 7. SAMS-VAE: Modeling cellular perturbations with the Sparse Additive Mechanism Shift VAE
**Authors/venue:** Bereket, Karaletsos (Insitro) — *NeurIPS* 2023. https://arxiv.org/abs/2311.02794
**Takeaway:** Factorizes each perturbation's effect as a *sparse, additive* latent shift with a learned mask over which latent dimensions the perturbation touches.
**Technical summary:** SAMS-VAE models a perturbed cell as basal latent + Σ (mask ⊙ perturbation-latent), with Bayesian sparsity so each intervention edits only a few latent coordinates — an explicit "which slots does this perturbation address" inference. This improves compositional generalization and yields interpretable, per-perturbation latent targets. It operationalizes the sparse-mechanism-shift hypothesis inside a generative perturbation model.
**Plain-English:** Every drug or knockout is assumed to nudge just a handful of the cell's latent knobs, and the model figures out which ones — making effects composable and interpretable.
**Applicability:** A2 (sparse mask = learned per-perturbation slot addressing — a direct mechanism for MORPHEUS's addressable slots); A5 (additive sparse counterfactuals). Design implication: MORPHEUS can use a SAMS-style sparse mask to *identify which pathway slots a query touches*, tying A2 addressability to A5 interventions.
**Novelty implication:** *pre-empts* "perturbation addresses a sparse identified set of latent programmes." MORPHEUS must add: (a) grounding those slots in *named pathways* (A3), (b) NL prompting over them (A1), and (c) a frozen multimodal tumor trunk (A4) — none present in SAMS-VAE.

### 8. sVAE+: Learning causal representations of single cells via sparse mechanism shift modeling
**Authors/venue:** Lopez, Tagasovska, Ra, Cho, Pritchard, Regev — *AISTATS/PMLR* 2023 (arXiv 2211.03553). https://arxiv.org/pdf/2211.03553
**Takeaway:** Treats genetic perturbations as sparse *unknown-target* interventions on latent factors and proves this yields identifiable causal latents.
**Technical summary:** Building on mechanism-sparsity identifiability, sVAE+ models each perturbation as shifting a sparse subset of latent variables, giving a Bayesian, low-tuning way to recover latent factors up to permutation from interventional single-cell data. It empirically disentangles biological factors better than non-causal VAEs on Perturb-seq. This is the single-cell instantiation of the identifiability-from-interventions program.
**Plain-English:** By assuming each genetic perturbation only jostles a few hidden factors, the model can pin down what those hidden factors actually are.
**Applicability:** A2 (identifiability via interventions — the theoretical backbone of MORPHEUS's "latent dims = real pathway states" claim); A5 (interventions are the *identifying signal*, not just the query). Design implication: MORPHEUS's perturbation-conditioned iVAE identification thesis is essentially sVAE+; adopt its sparse-shift objective as the identification mechanism.
**Novelty implication:** *pre-empts* the core identification claim (A2). MORPHEUS's defensible novelty is applying it to a *frozen multimodal tumor trunk with NL-addressable named slots and geodesic queries*, not the identification result itself.

### 9. FCR: Learning Identifiable Factorized Causal Representations of Cellular Responses
**Authors/venue:** Mao, Lopez, Liu, Hütter, Richmond, Benos, Qiu — arXiv 2024. https://arxiv.org/abs/2410.22472
**Takeaway:** Decomposes each cellular response into covariate-specific, treatment-specific, and *interaction* blocks with proven block-/component-wise identifiability.
**Technical summary:** FCR uses nonlinear-ICA theory to guarantee that the covariate, treatment, and covariate×treatment-interaction factors are identifiable (block-wise for main effects, component-wise for interactions). This explicitly models *context-dependence* — how a drug's effect changes with cell context — which is central to tumor-heterogeneity counterfactuals. It beats prior methods on single-cell perturbation datasets and surfaces context-dependent targets.
**Plain-English:** It cleanly separates "what the cell already was," "what the drug does," and "what the drug does *only in this kind of cell*," with math guaranteeing those pieces are recoverable.
**Applicability:** A2 (identifiability with explicit interaction terms); A5 (context-conditioned counterfactuals — exactly the tumor-subtype-specific drug-response question). Design implication: MORPHEUS should carry an explicit interaction/context factor so that "same drug, different tumor state → different answer" is identified, not confounded.
**Novelty implication:** *pre-empts* identified context×treatment factorization. MORPHEUS differentiates via pathway-named slots + NL routing + frozen multimodal trunk; FCR is the identifiability baseline to cite and beat on interaction recovery.

### 10. SENA-discrepancy-VAE: Interpretable Causal Representation Learning in the Pathway Space
**Authors/venue:** de la Fuente, Lehmann, … Lagani, Hernaez — *ICLR* 2025 (arXiv 2506.12439). https://arxiv.org/pdf/2506.12439
**Takeaway:** Makes each latent causal factor an interpretable combination of *learned biological-pathway activities*, at no cost to predictive accuracy.
**Technical summary:** SENA-discrepancy-VAE augments the discrepancy-VAE with a SENA layer that maps latent causal factors to (linear) combinations of pathway-activity scores, so each recovered causal factor is explicable in pathway terms. It matches the non-interpretable baseline on predicting unseen genetic and drug perturbations while yielding biologically meaningful causal latents. This is the closest published realization of "identified latent dims ↔ named pathways."
**Plain-English:** Its hidden causal knobs are each written as a recipe of known biological pathways, so you can read what each knob means without losing prediction quality.
**Applicability:** A2 (pathway-addressable identified slots — nearly the MORPHEUS A2 spec); A3 (latent↔biology grounding via pathway names); A5 (perturbation prediction on those factors). Design implication: MORPHEUS's pathway-slot layer can borrow SENA's pathway-projection; the remaining delta is *NL prompting*, *frozen multimodal trunk*, and *geodesic causal queries*.
**Novelty implication:** *pre-empts* the "identified + pathway-named + causal + perturbation-predictive" combination most completely of any single paper — **the single biggest novelty risk for A2/A3**. MORPHEUS must own A1 (NL task-inference) + A4 (multimodal/frozen-trunk) + A5-geometry to remain distinct.

### 11. Identifiability Guarantees for Causal Disentanglement from Soft Interventions
**Authors/venue:** Zhang, Greenewald, Squires, Srivastava, Shanmugam, Uhler — *NeurIPS* 2023. https://arxiv.org/abs/2307.06250
**Takeaway:** Proves latent causal variables and their causal graph are identifiable from *soft* (non-atomic) interventions — the regime real perturbations live in.
**Technical summary:** The paper (the base of discrepancy-VAE) shows that with one soft intervention per latent node, the latent causal model is identifiable up to a permutation/scaling, using a discrepancy objective between observational and interventional distributions. It formalizes when Perturb-seq-style soft knockdowns suffice to recover latent causal structure. This underwrites the whole "genetic perturbations identify pathway latents" premise.
**Plain-English:** Even when a genetic perturbation only partially nudges a hidden factor (not a clean on/off), you can still mathematically recover the hidden factors and their causal wiring.
**Applicability:** A2 (identifiability theory for the addressable-slot claim under realistic soft perturbations); A5 (interventions as identifying evidence). Design implication: MORPHEUS's identification claims should cite soft-intervention identifiability so they survive the fact that CRISPRi/drug perturbations are soft, not atomic.
**Novelty implication:** *strengthens* the theoretical legitimacy of MORPHEUS's identification axis while *pre-empting* the underlying theorem — the theory is not MORPHEUS's to claim, only to apply.

### 12. STATE: Arc Institute's first virtual cell model
**Authors/venue:** Arc Institute virtual-cell team (Adduri, Roohani, Hsu, et al.) — 2025. https://arcinstitute.org/news/virtual-cell-model-state (code: github.com/ArcInstitute/state)
**Takeaway:** Two-module virtual cell: a State-Embedding trunk (167M observational cells) + a State-Transition transformer (100M+ perturbed cells, 70 contexts) that predicts perturbation responses across contexts.
**Technical summary:** SE learns a noise-robust cell embedding from ~167M cells; ST is a bidirectional self-attention transformer over *cell sets* that predicts how a population transitions given a starting transcriptome + a perturbation, modeling heterogeneity without distributional assumptions. Trained on >100M perturbed cells spanning 70 cell contexts, it is the largest cross-context perturbation-response model to date and takes (state, perturbation) as an explicit query pair.
**Plain-English:** A giant model trained on hundreds of millions of cells that, given a cell and a perturbation, predicts how a whole population of such cells will change.
**Applicability:** A5 (industrial-scale perturbation-as-query across contexts); A1/A4 (separate embedding trunk vs. transition head = frozen-trunk + task module pattern). Design implication: STATE is the scale frontier MORPHEUS must position against — MORPHEUS's argument is *identifiability + NL promptability + tumor/clinical multimodality*, not out-scaling Arc.
**Novelty implication:** *pre-empts* "perturbation response across many contexts from a shared embedding" at massive scale; MORPHEUS cannot win on scale, so must win on identified/promptable/geodesic structure and clinical (WSI+multiomic) grounding.

### 13. Virtual Cell Challenge: Toward a Turing test for the virtual cell
**Authors/venue:** Bunne, Roohani, … (Arc Institute) — *Cell* 2025. https://www.cell.com/cell/fulltext/S0092-8674(25)00675-0
**Takeaway:** Defines the evaluation frame and community benchmark for perturbation-prediction ("virtual cell") models — a Turing-test-style held-out perturbation task.
**Technical summary:** The paper articulates what a virtual cell must do (predict unseen perturbations, generalize across contexts) and operationalizes it as a challenge (1,200+ teams, 70+ contexts), standardizing metrics for perturbation-response prediction. It is the field's current statement of *how interventional models should be evaluated*. Directly relevant to MORPHEUS's finding that its current harness is "structurally blind to representation quality."
**Plain-English:** A community-agreed exam for models that simulate cells, focused on predicting the effects of perturbations they were never trained on.
**Applicability:** A5 (evaluation of interventional queries); A3 (measuring what a representation actually "knows"). Design implication: MORPHEUS should adopt held-out-perturbation, cross-context evaluation from this challenge to replace confounded C-index-style scoring that is blind to representation quality.
**Novelty implication:** *reframes* MORPHEUS's evaluation gap as a solved-elsewhere problem — adopt these protocols rather than inventing new ones; the novelty is the *representation*, not the benchmark.

### 14. PerturBench: Benchmarking ML models for cellular perturbation analysis
**Authors/venue:** (GSK.ai / collaborators) — arXiv 2408.10609, 2024. https://arxiv.org/html/2408.10609v2
**Takeaway:** A unified benchmark showing many perturbation models barely beat simple baselines and rank-shuffle across metrics/splits.
**Technical summary:** PerturBench standardizes datasets, splits, and metrics for single-cell perturbation prediction and finds that model rankings are highly sensitive to the metric and the train/test split, with strong "predict the mean" baselines that many deep models fail to beat. It exposes evaluation fragility in exactly the counterfactual regime MORPHEUS targets. It provides leaderboard infrastructure and distributional metrics (e.g., energy distance, MMD).
**Plain-English:** When you carefully test perturbation-prediction models on equal footing, many of them barely outperform "just predict the average response."
**Applicability:** A5 (interventional evaluation rigor); A3 (representation-quality measurement). Design implication: MORPHEUS must beat strong mean/additive baselines under distributional metrics and multiple splits, or its counterfactual claims are unfalsifiable.
**Novelty implication:** *reframes* — sets the evidentiary bar; strengthens MORPHEUS's own critique that current evals are confounded, and warns that "we predict perturbations" is worthless without beating trivial baselines.

### 15. PertEval-scFM: Benchmarking single-cell foundation models for perturbation effect prediction
**Authors/venue:** (Wanjala, Zappia, et al.) — bioRxiv 2024.10.02.616248. https://www.biorxiv.org/content/10.1101/2024.10.02.616248
**Takeaway:** Foundation-model embeddings give little to no advantage over simple baselines for perturbation-effect prediction, especially for strongly perturbed genes.
**Technical summary:** PertEval-scFM probes frozen scFM embeddings (Geneformer, scGPT, etc.) with a standardized MLP "probe" for perturbation effects and finds they rarely beat a mean baseline, with gains concentrated in easy (weak-effect) genes. It isolates the *representation's* contribution by holding the probe fixed. This is a direct test of the premise that a good trunk enables good intervention prediction.
**Plain-English:** Big pretrained cell models don't actually carry much extra signal for predicting perturbation effects once you test them fairly.
**Applicability:** A4 (frozen-trunk plug-in probing — the exact methodology MORPHEUS uses); A5; A3 (measuring emergent knowledge). Design implication: MORPHEUS's frozen-trunk-probe evaluation should mirror PertEval; and it warns that a frozen trunk alone will *not* deliver intervention gains without identification/geometry.
**Novelty implication:** *reframes/strengthens* — validates MORPHEUS's diagnosis that scale-only trunks are "structurally blind," and raises the bar: MORPHEUS must show its *identified* trunk beats generic scFM trunks under identical probes.

### 16. scBIG: Beyond Independent Genes — Module-Inductive Representations for Perturbation Prediction
**Authors/venue:** Ruan, Quan, Xu, Yang, Yang — arXiv 2026 (q-bio.GN 2602.04901). https://arxiv.org/pdf/2602.04901
**Takeaway:** Predicts perturbation responses at the level of *coordinated gene programmes/modules* rather than independent genes, improving unseen/combinatorial generalization.
**Technical summary:** scBIG clusters genes into coherent programmes (Gene-Relation Clustering), encodes inter-programme interactions (Gene-Cluster-Aware Encoder), and adds structure-aware alignment to preserve modular coordination, yielding ~6.7% average improvement over the strongest baselines, concentrated on unseen and combinatorial perturbations. It argues perturbation effects are inherently program-level — the same premise behind MORPHEUS's pathway slots. Modules are learned, not imposed from an ontology.
**Plain-English:** Instead of predicting each gene separately, it groups genes into functional teams and predicts how the teams move, which helps most for perturbations it hasn't seen.
**Applicability:** A2 (programme-level representation = addressable slots); A5 (program-level counterfactuals). Design implication: MORPHEUS's `(batch, n_programme, D)` slot exposure is corroborated as the right granularity for intervention prediction; borrow module-coordination alignment losses.
**Novelty implication:** *pre-empts* the "program-level beats gene-level for perturbation" claim; MORPHEUS's differentiation is *identifiability of* those programmes + NL addressability, not merely modular grouping.

### 17. Geometric coherence of single-cell CRISPR perturbations (Shesha perturbation stability)
**Authors/venue:** Raju — arXiv 2026 (q-bio.QM 2604.16642). https://arxiv.org/pdf/2604.16642
**Takeaway:** A geometry-of-response metric (directional coherence S_p) predicts which perturbations cause reproducible, stress-inducing shifts vs. large-but-incoherent ones.
**Technical summary:** Shesha perturbation stability S_p measures the directional consistency of transcriptomic response vectors across cells under a perturbation; across 2,200+ perturbations in five datasets, pleiotropic regulators (e.g., CEBPA) show large-magnitude but incoherent shifts (a "geometric tax"), while lineage factors give coherent responses, and S_p predicts UPR/stress activation (p<10⁻¹⁸) better than magnitude. It shows response *geometry*, not just magnitude, encodes regulatory architecture. This is direct evidence that the geometry of interventional shifts is biologically meaningful.
**Plain-English:** Whether all cells respond to a perturbation by moving the *same direction* (not just moving a lot) tells you about the gene's role and whether it triggers cellular stress.
**Applicability:** A5 (causal-geometry of interventions — the geodesic/metric-tensor thesis); A2 (coherence as a property of identified programmes). Design implication: MORPHEUS's geodesic causal-distance should incorporate *directional coherence* of response vectors as a validated geometric signal, not just Euclidean magnitude.
**Novelty implication:** *strengthens* MORPHEUS's "geometry ≈ causal structure" claim with fresh empirical support, while *reframing* it: coherence/direction, not just geodesic length, is the biologically-loaded geometric quantity.

### 18. MapPFN: Learning Causal Perturbation Maps in Context
**Authors/venue:** Sextro, Kłos, Dernbach — arXiv 2026 (cs.LG 2601.21092). https://arxiv.org/abs/2601.21092
**Takeaway:** A prior-fitted network that does *in-context* (amortized, zero-shot) perturbation-effect estimation, adapting to a new biological context at inference with no fine-tuning.
**Technical summary:** MapPFN is pretrained on a synthetic biological prior with explicit causal interventions, then performs in-context learning: given a few observations from a new context, it maps sequential experimental evidence to post-perturbation distributions without dataset-specific retraining, matching supervised baselines zero-shot. It reframes perturbation prediction as *prompted inference* over a frozen model — architecturally the closest analog to MORPHEUS's "new question = new query, not a new model." Treatment effects emerge from the prompt/context, not from weights.
**Plain-English:** A model you "prompt" with a handful of measurements from a new cell context, and it immediately predicts perturbation effects there without being retrained.
**Applicability:** A1 (query-not-retrain; in-context task specification); A5 (counterfactuals as amortized inference); A4 (frozen model + contextual prompting). Design implication: strongly validates MORPHEUS's TQI premise — perturbation queries can be answered by a frozen model via context/prompt; MapPFN is a concrete blueprint.
**Novelty implication:** *pre-empts* "perturbation as a query on a frozen model, no retraining" (A1/A5 core). MORPHEUS must differentiate via NL (not few-shot numeric) prompting, identified pathway slots, and multimodal tumor grounding — MapPFN is numeric/synthetic-prior and single-modality.

### 19. Chem2Gen-Bench: Benchmarking Chemical-to-Genetic Translation in Perturbation Response Space
**Authors/venue:** Lin, Chen — arXiv 2026 (cs.LG 2606.21109). https://arxiv.org/pdf/2606.21109
**Takeaway:** Tests whether chemical and genetic perturbations produce *interchangeable* response signatures — the basis for translating drug effects into target-gene effects.
**Technical summary:** The benchmark aligns 260,084 chemical and 1,099,045 genetic perturbation profiles by cell-target context and evaluates pairwise alignment, retrieval, and foundation-model embeddings, finding translation fidelity is "measurable but heterogeneous," and that scFM embeddings often fail to beat simple gene-delta baselines in target-matched settings. It provides an auditable test of the "drug ≈ its target's knockdown" assumption central to intervention simulation. Background adjustment materially improves alignment.
**Plain-English:** It checks how often a drug's effect really looks like knocking out its target gene — sometimes yes, often only partly.
**Applicability:** A5 (chemical vs. genetic intervention equivalence — validity of drug-as-target-perturbation queries); A4 (when chemical structure adds signal over gene-delta). Design implication: MORPHEUS must *not* assume drug = clean target knockdown; the chemical→genetic map is context-dependent and should be an evaluated, hedged translation.
**Novelty implication:** *reframes* — cautions against a naive "drug = intervention on its target slot" shortcut; MORPHEUS's intervention simulator needs an explicit, auditable chem→gene translation rather than an identity assumption.

### 20. CITE-VAE: Learning Latent Dynamical Causal Processes for Single-Cell Perturbation Prediction
**Authors/venue:** Jiang, Liu, Gao, Abbasnejad, Yao, Shi — *SIGKDD* 2026 (AI4Science) (arXiv 2605.25581). https://arxiv.org/pdf/2605.25581
**Takeaway:** Adds *temporal dynamics* to identifiable perturbation modeling — latent programmes + perturbation mechanisms evolve over time with recoverability guarantees.
**Technical summary:** The model integrates latent cellular programmes, perturbation-conditioned mechanisms, and temporal dynamics, with identifiability analysis for when latent causal variables are recoverable up to standard equivalence classes; CITE-VAE learns these from sequencing data and improves generalization to unseen perturbations on synthetic and CRISPR data. It extends static identifiable-perturbation models into a dynamical/causal-process regime. Time is the extra identifying signal.
**Plain-English:** It models not just the end state of a perturbed cell but the causal process over time, and proves the hidden factors driving it are recoverable.
**Applicability:** A2 (identifiability, now dynamical); A5 (counterfactual *trajectories*, not just endpoints). Design implication: if MORPHEUS ever needs response-over-time or dose-time counterfactuals, dynamical identifiability is the framework; otherwise a static special case.
**Novelty implication:** *strengthens/pre-empts* the identifiable-perturbation-latents claim in the temporal direction; MORPHEUS's static tumor-state framing is a subset — differentiate by multimodality + NL, not by dynamics.

### 21. PDGrapher: Combinatorial prediction of therapeutic perturbations via causally-inspired neural networks
**Authors/venue:** (Zitnik lab; Gonzalez et al.) — *Nature Biomedical Engineering* 2025 (bioRxiv 2024.01.03.573985). https://www.nature.com/articles/s41551-025-01481-x
**Takeaway:** Solves the *inverse* intervention problem — given a diseased and desired state, predict which combination of targets to perturb to get there.
**Technical summary:** PDGrapher casts genes as nodes in a structural causal model approximated by PPI/GRN graphs, with a perturbagen-discovery GNN proposing target sets and a response GNN evaluating them, assuming no unobserved confounders. Across 19 datasets (11 cancer types, chemical + genetic interventions) it outperforms prior methods at identifying phenotype-reversing perturbagens. It is the "what intervention achieves outcome Y?" query — the inverse of GEARS/CPA. This is exactly a clinical "what should we do?" counterfactual.
**Plain-English:** Instead of predicting a drug's effect, it works backwards: given a sick cell and a target healthy state, it proposes which genes/drugs to hit to get there.
**Applicability:** A5 (inverse/prescriptive counterfactual — "which intervention?" as a query); A3 (grounding to actionable targets). Design implication: MORPHEUS's intervention interface should support *both* forward ("effect of this drug") and inverse ("which drug moves this tumor toward remission") queries — PDGrapher is the reference for the inverse mode.
**Novelty implication:** *pre-empts* prescriptive/inverse intervention prediction; MORPHEUS differentiates by doing it as an NL query on a *frozen identified* representation rather than a bespoke graph model with an assumed causal graph.

### 22. GPerturb: Gaussian-process modeling of single-cell perturbation data
**Authors/venue:** (Wu et al.) — *Nature Communications* 2025. https://www.nature.com/articles/s41467-025-61165-7
**Takeaway:** A structured GP that separates perturbation-induced from baseline variation and gives *sparse, interpretable, uncertainty-quantified* gene-level effects.
**Technical summary:** GPerturb decomposes expression into a baseline component and a sparse additive perturbation effect via a Bayesian GP, identifying which genes each perturbation significantly changes with calibrated uncertainty. It offers interpretability and UQ that neural counterfactual models often lack, on both CRISPR and drug data. Sparsity again localizes effects to few dimensions.
**Plain-English:** A statistical model that says, with error bars, exactly which genes a perturbation actually changes versus normal cell-to-cell noise.
**Applicability:** A5 (counterfactual effects with uncertainty); A2 (sparse effect localization); A3 (interpretability). Design implication: MORPHEUS's intervention outputs should carry calibrated uncertainty and sparse attributions — an abstention/UQ requirement for trustworthy prompting.
**Novelty implication:** *strengthens* the case that intervention answers need UQ + sparsity; *reframes* MORPHEUS toward reporting calibrated, abstaining counterfactuals rather than point predictions.

### 23. Enforcing Latent Euclidean Geometry in Single-Cell VAEs for Manifold Interpolation
**Authors/venue:** (Palma, Theis, et al.) — arXiv 2507.11789, 2025. https://arxiv.org/pdf/2507.11789
**Takeaway:** Regularizes a single-cell VAE so that *straight lines in latent space* are biologically valid interpolations, improving trajectory/interpolation fidelity.
**Technical summary:** The method (FlatVI / Euclidean-latent VAE) trains the decoder so the pullback metric is approximately Euclidean, making linear latent interpolation match geodesics of the data manifold and improving latent + gene-wise trajectory reconstruction in population-dynamics modeling. It directly engages the metric-tensor / geodesic question: rather than compute geodesics, it *flattens* the latent so Euclidean = geodesic. This is a concrete design choice on MORPHEUS's causal-geometry axis.
**Plain-English:** It reshapes the cell latent space so that drawing a straight line between two cells actually traces a realistic biological path.
**Applicability:** A5 (latent geometry for interventional interpolation — the geodesic-causal-distance axis). Design implication: MORPHEUS faces a fork — *compute* geodesics under a learned metric (expensive) vs. *train* the latent flat so Euclidean moves are causal (this paper); the latter may make interventions cheap and geodesic-faithful.
**Novelty implication:** *reframes* MORPHEUS's metric-tensor plan: geometry can be *baked into training* rather than computed at query time; a serious alternative design MORPHEUS must weigh and cite.

### 24. TopOMetry: systematically learning and evaluating the latent geometry of single-cell data
**Authors/venue:** (Sidarta-Oliveira et al.) — *eLife* reviewed preprint 2024 (bioRxiv 2022.03.14.484134). https://elifesciences.org/reviewed-preprints/100361
**Takeaway:** A framework that models cell state as a Riemannian manifold via Laplacian-type operators and *evaluates* how well an embedding preserves that geometry.
**Technical summary:** TopOMetry builds multiple candidate latent geometries from Laplacian-type operators and scores them for how faithfully they capture the data manifold, providing systematic geometry evaluation rather than assuming a single embedding. It formalizes "cell landscapes are Riemannian manifolds" and gives tooling to *measure* geometric fidelity — relevant to evaluating whether MORPHEUS's latent respects causal geometry. It exposes that default embeddings (incl. VAE) can distort manifold structure.
**Plain-English:** A toolkit that treats the space of cell states as a curved surface and checks which of many possible maps of it least distorts the true shape.
**Applicability:** A5 (geometry as a first-class, evaluated property); A3 (evaluating representation geometry). Design implication: MORPHEUS should *measure* geodesic/manifold fidelity of its latent (à la TopOMetry) rather than assert it — closes the "evaluation blind to representation quality" gap on the geometry axis.
**Novelty implication:** *reframes* MORPHEUS's causal-geometry claim into a measurable one and warns that VAE latents distort geometry — MORPHEUS must demonstrate, not assume, geodesic≈causal.

### 25. Fast Approximate Geodesics for Deep Generative Models
**Authors/venue:** Chen, Klushyn, Kurle, Jiang, Bayer, van der Smagt — *ICANN* 2019 (arXiv 1812.08284). https://arxiv.org/pdf/1812.08284
**Takeaway:** Foundational method for computing geodesics on a generative model's latent manifold via the decoder-Jacobian pullback metric, made tractable with a fast approximation.
**Technical summary:** It endows a VAE latent with a Riemannian metric = pullback of the ambient metric through the decoder Jacobian, then approximates geodesics efficiently (avoiding expensive ODE solves) so that latent distances reflect data geometry, not raw Euclidean coordinates. This is the canonical machinery MORPHEUS's "metric tensor → geodesic ≈ causal distance" would instantiate. It establishes both the correctness and the cost problem of query-time geodesics.
**Plain-English:** A way to measure "real" distances on the curved surface a generative model has learned, cheaply enough to be practical.
**Applicability:** A5 (the literal metric-tensor/geodesic tooling for causal-geometry queries). Design implication: gives MORPHEUS the concrete pullback-metric + fast-geodesic recipe, and quantifies the compute cost that motivates alternatives like flattened latents (entry 23).
**Novelty implication:** *pre-empts* the *mechanism* of geodesic latent distances (not novel); MORPHEUS's novelty must be the *causal interpretation* (geodesic length ≈ causal/interventional distance on a tumor manifold), which this paper does not claim.

### 26. Identifying latent distances with Finslerian geometry
**Authors/venue:** Pouplin, Eklund, Hauberg, et al. — arXiv 2212.10010, 2022/2023. https://arxiv.org/pdf/2212.10010
**Takeaway:** Argues the correct latent metric for stochastic generative models is *Finslerian*, not Riemannian, giving better-behaved geodesics under decoder uncertainty.
**Technical summary:** For generative models with stochastic decoders, the induced geometry is more faithfully Finslerian; the paper derives a Finsler metric and geodesics that account for decoder variance, improving distance estimates over the naive Riemannian pullback. It refines the theory MORPHEUS's causal-geometry axis rests on. Relevant because tumor-state decoders are stochastic.
**Plain-English:** A more careful math for measuring distances in uncertain generative latent spaces, improving on the standard curved-distance approach.
**Applicability:** A5 (latent-geometry rigor for causal distance). Design implication: if MORPHEUS uses geodesic causal distance with a stochastic decoder, Finslerian (not Riemannian) geometry may be the technically correct choice.
**Novelty implication:** *reframes* the geometric-distance machinery; a technical refinement MORPHEUS should be aware of so its geometry claim is stated correctly for stochastic decoders.

### 27. CODE-AE: Context-aware deconfounding autoencoder for personalized drug-response prediction
**Authors/venue:** He, Liu, Zhang, Xie — *Nature Machine Intelligence* 2022. https://www.nature.com/articles/s42256-022-00541-0
**Takeaway:** Disentangles biological signal from confounders across domains (cell lines → patients) to transfer drug-response predictions to real tumors.
**Technical summary:** CODE-AE learns shared vs. private latent factors and adversarially removes context/confounder signal so a model trained on cell-line compound screens transfers to patient tumors, improving personalized drug-response prediction. It directly addresses the cell-line→patient generalization gap central to clinical intervention queries. Deconfounding = separating causal drug effect from domain artifacts.
**Plain-English:** It cleans out the differences between lab cell lines and real patients so that drug-response predictions learned in the lab actually work on patients.
**Applicability:** A5 (drug-response counterfactuals that transfer to patients); A4 (multi-domain/multimodal integration); A2 (shared/private factor separation). Design implication: MORPHEUS's drug-response queries on tumor states must deconfound domain/context, or cell-line-trained interventions won't hold clinically — build a shared/private split.
**Novelty implication:** *pre-empts* deconfounded cell-line→patient drug transfer; MORPHEUS adds the frozen multimodal tumor trunk + NL query, but must cite CODE-AE as the transfer baseline.

### 28. Enhancing drug and cell-line representations via contrastive learning for anti-cancer drug prioritization
**Authors/venue:** (Nguyen et al.) — *npj Precision Oncology* 2024. https://www.nature.com/articles/s41698-024-00589-8
**Takeaway:** Contrastive learning on DepMap expression + drug structure yields representations that improve drug prioritization and transfer.
**Technical summary:** The method contrastively aligns cell-line (DepMap gene-expression) and drug representations so that response-relevant structure is captured, improving anti-cancer drug ranking over supervised baselines. It shows representation *quality* (via contrastive objectives) drives drug-response performance — a representation-first stance MORPHEUS shares. Uses expression as the information-rich, widely available modality.
**Plain-English:** By teaching the model to pull matching drug–cell pairs together in latent space, it builds better representations for ranking which drugs will work.
**Applicability:** A4 (which modality to encode — expression as primary; drug as encoded partner); A5 (drug prioritization as a representation query). Design implication: supports MORPHEUS's frozen-representation + contrastive-alignment approach for drug queries, and the choice of expression as the load-bearing molecular modality.
**Novelty implication:** *strengthens* the representation-first bet for drug response; not directly competing on causality/geometry, so low novelty risk.

### 29. DIPK: Improving drug-response prediction by integrating gene relationships with deep learning
**Authors/venue:** (Li et al.) — *Briefings in Bioinformatics* 2024 (PMC11006795). https://pmc.ncbi.nlm.nih.gov/articles/PMC11006795/
**Takeaway:** Self-supervised integration of gene-interaction networks, expression, and molecular topology for robust cell-line drug-response prediction.
**Technical summary:** DIPK fuses prior gene-relationship knowledge (via self-supervised graph pretraining), expression profiles, and drug molecular graphs to predict IC50-style responses with improved accuracy and robustness on DepMap/GDSC. It exemplifies knowledge-guided multimodal fusion for drug response. Prior gene-relationship structure is injected rather than identified.
**Plain-English:** It combines what we know about gene interactions with expression and drug chemistry to more reliably predict how a cell line responds to a drug.
**Applicability:** A4 (multimodal fusion: expression + drug graph + prior network); A5 (drug-response prediction). Design implication: a fusion baseline; contrasts with MORPHEUS's aim to *identify* structure internally rather than inject an external gene network.
**Novelty implication:** *pre-empts* multimodal knowledge-guided drug-response prediction (a crowded area); reinforces that MORPHEUS's differentiation is identification/promptability, not multimodal fusion accuracy alone.

### 30. Fast, scalable Wasserstein-1 neural optimal-transport solver for single-cell perturbation prediction
**Authors/venue:** (Uscidda, Cuturi, et al.) — 2025 (PMC12261427). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12261427/
**Takeaway:** A scalable W1 neural-OT solver that makes optimal-transport perturbation prediction (à la CellOT) fast enough for large screens.
**Technical summary:** It solves the Wasserstein-1 OT problem with a neural potential that scales to large single-cell perturbation datasets, predicting control→perturbed distribution maps more efficiently than input-convex W2 solvers. It keeps the geometry-native (transport-map) view of interventions while removing the scalability bottleneck. Interventions remain distributional maps, not point predictions.
**Plain-English:** A faster engine for the "morph untreated cells into treated cells" optimal-transport approach, so it works on big datasets.
**Applicability:** A5 (scalable geometric/transport counterfactuals). Design implication: if MORPHEUS models interventions as transport/geodesic maps, this is the scalable solver template; reinforces OT as a practical intervention geometry.
**Novelty implication:** *reframes* — geometric interventions are now scalable, so MORPHEUS can't claim OT-geometry counterfactuals are impractical; must show its *identified/named-slot* geometry adds interpretability OT lacks.

### 31. Disentanglement via Mechanism Sparsity Regularization: a new principle for nonlinear ICA
**Authors/venue:** Lachapelle, Rodríguez López, … Lacoste-Julien — *CLeaR* 2022 (arXiv 2107.10098). https://arxiv.org/pdf/2107.10098
**Takeaway:** Foundational theory: latent factors are identifiable (up to permutation) if their mechanisms depend *sparsely* on actions/auxiliary variables and a graphical criterion holds.
**Technical summary:** The paper proves nonlinear-ICA identifiability under mechanism sparsity — when interventions/actions touch few latent factors, those factors are recoverable — and shows unknown-target interventions can disentangle latents. This is the theoretical parent of sVAE+ and SAMS-VAE and thus of MORPHEUS's identification-from-perturbations claim. It ties identifiability to the *sparsity of interventions*, exactly the Perturb-seq regime.
**Plain-English:** If each action only affects a few hidden factors, you can mathematically recover those factors — a general recipe behind identifiable perturbation models.
**Applicability:** A2 (the identifiability theorem underpinning addressable, identified slots). Design implication: MORPHEUS's identification objective should invoke mechanism-sparsity conditions to justify that its latent dims = real pathway states.
**Novelty implication:** *pre-empts* (owns) the identifiability principle itself; MORPHEUS applies, not invents, it — novelty must lie downstream (naming, prompting, geometry, multimodality).

### 32. Nonparametric Partial Disentanglement via Mechanism Sparsity
**Authors/venue:** Lachapelle, Mahajan, Mitliagkas, Lacoste-Julien — *JMLR* 2024. https://www.jmlr.org/papers/volume27/24-0771/24-0771.pdf
**Takeaway:** Extends mechanism-sparsity identifiability to nonparametric settings with sparse actions, interventions, and temporal dependencies — *partial* (block) disentanglement guarantees.
**Technical summary:** It generalizes the 2022 result to nonparametric mechanisms and characterizes exactly which groups of latents become identifiable ("partial disentanglement") given the sparsity graph of interventions/temporal links. This tells you *how much* identifiability realistic sparse-intervention data buys — often block-level, not full component-level. Directly bounds what MORPHEUS can honestly claim about slot identifiability.
**Plain-English:** A more general and realistic version of the sparsity-identifiability theory that tells you which *groups* of hidden factors you can pin down, not just whether you can.
**Applicability:** A2 (realistic, partial identifiability of slot groups). Design implication: MORPHEUS should claim *block/pathway-group* identifiability (what soft, sparse tumor perturbations support), not per-dimension identification — this paper sets the honest ceiling.
**Novelty implication:** *reframes* MORPHEUS's A2 claim toward *partial* (programme-block) identifiability; over-claiming full disentanglement would be refuted by this theory.

### 33. PertAdapt: Unlocking single-cell foundation models for perturbation via condition-sensitive adaptation
**Authors/venue:** (bioRxiv 2025.11.21.689655), 2025. https://www.biorxiv.org/content/10.1101/2025.11.21.689655
**Takeaway:** A lightweight, condition-sensitive adapter that finally lets frozen scFM trunks beat GEARS on perturbation prediction.
**Technical summary:** PertAdapt inserts a condition-sensitive adaptation module on top of a frozen single-cell foundation model so the trunk's representations become perturbation-informative without full fine-tuning, reportedly closing/reversing the gap where scFMs previously lagged GEARS. It is a direct instance of the frozen-trunk + plug-in-adapter pattern for interventions. The trunk stays frozen; only a small conditioning module trains.
**Plain-English:** A small add-on that makes a big frozen cell model good at predicting perturbation effects, without retraining the whole thing.
**Applicability:** A4 (frozen-trunk plug-in for a new task — precisely MORPHEUS's adapter/head design); A5; A1 (task-specific conditioning module ≈ routed head). Design implication: strong evidence MORPHEUS's frozen-trunk + intervention-head architecture is viable *if* the trunk is good; adapter conditioning is the mechanism.
**Novelty implication:** *pre-empts* "frozen-trunk + adapter for perturbation prediction"; MORPHEUS differentiates by the trunk being *identified/multimodal* and the adapter being *NL-promptable*, not a fixed perturbation head.

### 34. Scalable single-cell gene-expression generation with latent diffusion models
**Authors/venue:** (Palma, Theis, et al.) — arXiv 2511.02986, 2025. https://arxiv.org/pdf/2511.02986
**Takeaway:** Latent-diffusion generation of single-cell expression, including conditional/counterfactual generation under covariates and perturbations.
**Technical summary:** The model performs diffusion in a learned single-cell latent to generate realistic expression profiles conditioned on covariates/perturbations, offering high-fidelity conditional sampling and counterfactual generation at scale. It represents the diffusion-based alternative to VAE/OT for intervention simulation, with strong sample fidelity. Conditioning enters via classifier-free guidance rather than additive latent shifts.
**Plain-English:** A diffusion (image-generator-style) model for cells that can synthesize what cells would look like under given conditions or perturbations.
**Applicability:** A5 (counterfactual generation via conditional diffusion — an alternative intervention engine); A4 (conditioning as multimodal control). Design implication: diffusion guidance is a candidate mechanism for MORPHEUS's intervention simulator, distinct from additive (CPA) or transport (CellOT) composition.
**Novelty implication:** *reframes* the design space of intervention simulators (VAE-additive vs. OT-map vs. diffusion-guidance); MORPHEUS should position its choice against all three and justify identifiability/geometry, which diffusion alone does not provide.

---

### Lane synthesis — where the A5 novelty actually lives

- **Perturbation-as-query is thoroughly pre-empted** (GEARS, CPA/chemCPA, biolord, STATE, MapPFN, PertAdapt). MORPHEUS gets *no* novelty for "predict a perturbation/drug response from a representation." Its A5 claim must be the *combination*: frozen + identified + pathway-named + NL-promptable + geodesic — on a **multimodal tumor** trunk, not scRNA alone.
- **Identifiability-from-interventions is owned by the causal-representation literature** (Lachapelle mechanism-sparsity; sVAE+; SAMS-VAE; FCR; soft-intervention identifiability). MORPHEUS applies these; it must claim only *partial/block (pathway-group)* identifiability to stay honest (JMLR 2024).
- **The single sharpest novelty risk is SENA-discrepancy-VAE (entry 10)**: identified + pathway-named + causal + perturbation-predictive in one model. MORPHEUS's remaining defensible ground is A1 (NL task auto-detection/routing) + A4 (multimodal encode-vs-RAG, frozen trunk) + geodesic causal geometry — not identifiability or pathway-naming alone.
- **Geometry is real but must be measured, not asserted** (Shesha coherence; TopOMetry; pullback-metric geodesics; Finsler; flattened-latent VAEs). "Geodesic ≈ causal distance" is a *measurable, contestable* hypothesis; MORPHEUS must evaluate it (TopOMetry-style) and pick a design (compute geodesics vs. flatten the latent).
- **Evaluation is the field's live wound** (PerturBench, PertEval-scFM, Virtual Cell Challenge): strong baselines beat many models; foundation trunks add little. This *validates* MORPHEUS's own diagnosis that its harness is blind to representation quality — and sets the bar: beat mean/additive baselines under distributional metrics and held-out perturbations, with a frozen-trunk probe.
