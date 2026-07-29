## Cross-modal alignment, identifiability & disentanglement theory

Lane id: `l09_alignment_identifiability`. Remit: theory that tells us when a learned latent is *identified* (recoverable up to a benign transform), when cross-modal alignment recovers *shared* structure, and when disentanglement/addressability is even possible. This is the theoretical spine under MORPHEUS axis A2 (identified, pathway-addressable slots that make prompting reliable) and A5 (interventional/causal-geometry queries), with strong bearing on A4 (multimodal encode-vs-context) and A1 (a promptable unified representation only "routes reliably" if its slots are stable across runs).

Cross-cutting takeaway for MORPHEUS: **unconstrained nonlinear representation learning is non-identifiable** (Locatello, Hyvärinen). Identifiability is *bought* with structure — auxiliary conditioning, multiple views/modalities, temporal/sparsity mechanisms, or interventions. Every one of these "prices" maps onto a MORPHEUS design lever (conditioning tokens, multimodal pairing, perturbation data). The novelty risk cuts both ways: much of what MORPHEUS wants to claim as "addressable programme slots" is *already theoretically characterized*, so MORPHEUS must claim the *engineering synthesis + biological grounding + prompting interface*, not the identifiability primitive itself.

---

### Group A — Foundational nonlinear ICA & identifiability of learned representations

**1. Unsupervised Feature Extraction by Time-Contrastive Learning and Nonlinear ICA** (Hyvärinen & Morioka, NIPS 2016) — https://papers.nips.cc/paper/2016/hash/d305281faf947ca7acade9ad5c8c818c-Abstract.html
- *Takeaway:* The first constructive proof that nonlinear ICA is identifiable if you exploit non-stationarity (time segments as an auxiliary label).
- *Technical summary:* Trains a feature extractor by discriminating time segments (multinomial logistic regression on segment index); shows the optimal features recover the independent sources up to a linear transformation when sources are non-stationary across segments. Establishes the template "auxiliary variable + contrastive discrimination ⇒ identifiability" that all later results generalize.
- *Plain-English:* If different chunks of your data have different statistics, a classifier trained to tell the chunks apart is secretly recovering the true hidden factors.
- *Applicability:* A2 — "conditioning context" (batch, tissue, patient, timepoint) is not a nuisance to be normalized away but the *auxiliary variable that buys identifiability*. Design implication: MORPHEUS should treat metadata tokens as identifiability drivers, not just covariates.
- *Novelty implication:* Pre-empts any claim that "our slots are identifiable because we condition on context" — that mechanism is 2016. MORPHEUS must show the biological *addressability* is new, not the conditioning trick.

**2. Nonlinear ICA Using Auxiliary Variables and Generalized Contrastive Learning** (Hyvärinen, Sasaki & Turner, AISTATS 2019) — https://arxiv.org/abs/1805.08651
- *Takeaway:* Unifies time-, non-stationarity-, and label-conditioned nonlinear ICA into one auxiliary-variable framework with a full identifiability + estimation-consistency proof.
- *Technical summary:* Sources are conditionally independent given an auxiliary variable u (time, history, class, or any side info); a discriminator separating real (x,u) from randomized (x,u*) provably recovers the sources. Provides the general theorem later specialized by iVAE.
- *Plain-English:* Almost any extra observed signal — a label, a timestamp, another sensor — can be used to pin down the true hidden variables that generated your data.
- *Applicability:* A2, A4 — any co-observed modality (proteomics, CNV, clinical) can serve as the auxiliary variable u for RNA sources. Directly motivates A4's "encode a modality to sharpen identifiability" over "treat as RAG context."
- *Novelty implication:* Reframes MORPHEUS's multimodal encoding as *identifiability-by-auxiliary-conditioning*, a known principle — strengthens the theoretical footing but caps the novelty of the mechanism.

**3. Variational Autoencoders and Nonlinear ICA: A Unifying Framework (iVAE)** (Khemakhem, Kingma, Monti & Hyvärinen, AISTATS 2020) — https://arxiv.org/abs/1907.04809
- *Takeaway:* Deep latent-variable models are identifiable up to simple (permutation + pointwise) transforms *iff* the prior is factorized and conditioned on an auxiliary variable.
- *Technical summary:* Introduces iVAE: a VAE with a conditionally-factorized exponential-family prior p(z|u). Proves the true joint is recovered up to an equivalence class, giving disentanglement with universal-approximator decoders — including "vanilla" VAEs as a special (non-identifiable) case, exposing why standard VAEs entangle.
- *Plain-English:* A standard autoencoder can scramble the meaning of its latent knobs; adding a conditioning signal and the right prior makes each knob correspond to a real underlying factor.
- *Applicability:* A2 — the canonical recipe for "slots that mean the same thing every run." MORPHEUS's per-programme slots need a conditional prior p(slot | context) to be identified, not just an L1/orthogonality penalty.
- *Novelty implication:* This *is* the identifiability primitive MORPHEUS depends on. Claim novelty in biological instantiation + prompting, not the theorem.

**4. On Linear Identifiability of Learned Representations** (Roeder, Metz & Kingma, ICML 2021) — https://arxiv.org/abs/2007.00810
- *Takeaway:* A broad family of discriminative/softmax models (incl. contrastive, InfoNCE) is identifiable *up to a single linear transform* across runs and architectures.
- *Technical summary:* Shows that despite overparameterization, models trained with softmax-type objectives converge to representations related by a linear map — so a learned linear probe transfers across seeds/architectures. Connects contrastive learning to the iVAE identifiability class.
- *Plain-English:* Two networks trained the same way learn representations that differ only by a rotation/stretch — so a simple linear readout is enough to line them up.
- *Applicability:* A2, A1 — justifies MORPHEUS reading out pathway programmes with *linear* probes/prompts and expecting cross-run stability. Design implication: slot addressability can be a learned linear map, and "linear identifiability" is the right stability metric to report.
- *Novelty implication:* Strengthens A2's feasibility but pre-empts "our slots are uniquely learned" — expect a *linear* indeterminacy, so MORPHEUS should measure and quotient it out, not claim exact recovery.

**5. Contrastive Learning Inverts the Data Generating Process** (Zimmermann, Sharma, Schneider, Bethge & Brendel, ICML 2021) — https://arxiv.org/abs/2102.08850
- *Takeaway:* InfoNCE-trained encoders implicitly invert the generative model, recovering latents up to affine/orthogonal transforms under stated assumptions.
- *Technical summary:* Formalizes the loss's optimum as the inverse of the ground-truth mixing under a von Mises–Fisher/uniform latent model; validates on the 3DIdent benchmark. Ties contrastive learning to nonlinear ICA and prescribes loss variants that match different latent geometries.
- *Plain-English:* Contrastive learning doesn't just cluster data — under the right conditions it reconstructs the actual hidden factors that generated it.
- *Applicability:* A2, A4 — if MORPHEUS uses contrastive multimodal pairing, its shared latent can be argued to *recover* biology, not just correlate. The latent geometry (sphere vs box) must match the assumed prior — a concrete design choice.
- *Novelty implication:* Supports A4's "encode to recover shared biology" claim with theory; the geometry-matching requirement is a design constraint MORPHEUS should state explicitly.

**6. Identifiability of Deep Generative Models Without Auxiliary Information** (Kivva, Rajendran, Ravikumar & Aragam, NeurIPS 2022) — https://arxiv.org/abs/2206.10044
- *Takeaway:* With a *mixture* (e.g., clustered) prior and piecewise-linear (ReLU) decoders, deep generative models are identifiable *without* any auxiliary label.
- *Technical summary:* Proves an identifiability hierarchy for VAE decoders with Gaussian-mixture / non-factorized priors and ReLU networks, generalizing iVAE beyond the conditioning requirement. Covers practical architectures (VaDE, MFC-VAE).
- *Plain-English:* If your data naturally falls into clusters, you can recover the true generative factors even without extra labels.
- *Applicability:* A2, A4 — cell-type / cell-state clustering in scRNA is exactly the mixture structure this needs. Implies MORPHEUS can claim identifiable slots even for unlabeled data, provided the latent is modeled as a mixture (cell states) rather than a single Gaussian.
- *Novelty implication:* Reframes — MORPHEUS can drop the "must condition on labels" caveat where cell-state mixtures exist, a useful hedge for A4's frozen-trunk/label-scarce regime.

**7. Towards Nonlinear Disentanglement in Natural Data with Temporal Sparse Coding (SlowVAE)** (Klindt, Schott, Sharma, Ustyuzhaninov, Brendel, Bethge & Paiton, ICLR 2021) — https://arxiv.org/abs/2007.10930
- *Takeaway:* Temporal *sparsity* of factor transitions (Laplacian-slow priors) yields identifiable disentanglement in natural sequences.
- *Technical summary:* SlowVAE places a sparse (heavy-tailed) transition prior over consecutive frames so that few factors change at a time; proves/empirically shows recovery of ground-truth factors on natural video.
- *Plain-English:* In the real world only a few things change from moment to moment; assuming that sparsity lets a model discover what those things are.
- *Applicability:* A2, A5 — trajectories / pseudotime / drug time-courses satisfy "few programmes change per step." Motivates a sparse-transition prior for MORPHEUS's temporal or perturbation-response slots.
- *Novelty implication:* Pre-empts a naive "sparsity gives us disentangled programmes" claim; the principle exists — MORPHEUS should cite it and extend to biological interventions.

---

### Group B — Disentanglement: limits & the structures that break them

**8. Challenging Common Assumptions in the Unsupervised Learning of Disentangled Representations** (Locatello, Bauer, Lucic, Rätsch, Gelly, Schölkopf & Bachem, ICML 2019, best paper) — https://arxiv.org/abs/1811.12359
- *Takeaway:* Fully unsupervised disentanglement is *provably impossible* without inductive bias on model and data; 12,000-model study finds no unsupervised metric reliably picks disentangled models.
- *Technical summary:* Gives an impossibility theorem (infinitely many equally-fitting entangled generative models exist for any factorized latent) and a large reproducible benchmark showing random seed / hyperparameters act as implicit supervisors. Concludes future work must be explicit about inductive biases and supervision.
- *Plain-English:* You cannot get meaningful, human-aligned "knobs" from data alone — some assumption or supervision must be baked in, or you're just picking one of infinitely many arbitrary answers.
- *Applicability:* A2 (the load-bearing caveat) — MORPHEUS *cannot* claim identifiable pathway slots from unsupervised scRNA alone. Every addressability claim must name its inductive bias (pathway priors, perturbations, multimodal pairing, conditioning).
- *Novelty implication:* The central novelty *guardrail* for the whole rebase. If A2's slots are "just learned," reviewers will invoke this. MORPHEUS must foreground the specific bias (biological priors/interventions) that makes its slots identified.

**9. Disentanglement via Mechanism Sparsity Regularization: A New Principle for Nonlinear ICA** (Lachapelle, Rodríguez, Sharma, Everett, Le Priol, Lacoste & Lacoste-Julien, CLeaR 2022) — https://arxiv.org/abs/2107.10098
- *Takeaway:* When latent factors depend *sparsely* on past factors / actions, jointly learning the factors and their sparse causal graph yields identifiability up to permutation.
- *Technical summary:* VAE with learned binary masks over a latent transition/mechanism graph; an L0-style sparsity penalty on the mechanism recovers factors provably. Extended (JMLR 2025) to sparse actions, interventions, and temporal dependencies.
- *Plain-English:* If each cause touches only a few effects, enforcing that sparseness reveals what the individual causes are.
- *Applicability:* A2, A5 — perturbations (drug, CRISPR) are *sparse actions* on latent programmes; MORPHEUS can regularize the perturbation-to-slot mechanism to be sparse for identified, addressable slots.
- *Novelty implication:* This is arguably the closest prior to MORPHEUS's "pathway-addressable slots under perturbation." Strong novelty risk — MORPHEUS must differentiate via NL-prompting + emergent-knowledge grounding, not the sparsity principle.

**10. Synergies Between Disentanglement and Sparsity: Generalization and Identifiability in Multi-Task Learning** (Lachapelle, Deleu, Mahajan, Mitliagkas, Bengio, Lacoste-Julien & Bertrand, ICML 2023) — https://arxiv.org/abs/2211.14666
- *Takeaway:* A sparse task-specific readout on top of a disentangled representation both improves generalization *and* is itself an identifiability condition.
- *Technical summary:* Proves that if downstream tasks use sparse subsets of latents, requiring sparse task-predictors provably disentangles the representation (permutation-identifiability), connecting multi-task learning to nonlinear ICA.
- *Plain-English:* If each task only needs a few of the hidden factors, forcing the task-heads to be sparse forces the factors to be clean and separable.
- *Applicability:* A1, A2 — MORPHEUS's promptable tasks each engage a few programmes; enforcing sparse task-routing is simultaneously a routing mechanism (A1) and a slot-identification mechanism (A2). Elegant dual justification for prompt-routing.
- *Novelty implication:* Strengthens A1↔A2 coupling with theory, but shows "sparse task heads disentangle" is known — MORPHEUS's contribution is the NL task auto-detection layer, not the sparsity-identifiability link.

**11. Independent Mechanism Analysis, a New Concept?** (Gresele, von Kügelgen, Stimper, Schölkopf & Besserve, NeurIPS 2021) — https://arxiv.org/abs/2106.05200
- *Takeaway:* Adds a causal "independent-influences" constraint (orthogonality of the mixing Jacobian columns) that resolves nonlinear-ICA nonidentifiability where statistical independence alone fails.
- *Technical summary:* Shows statistical independence of sources is insufficient for nonlinear identifiability; proposes IMA — each source influences the observation via a geometrically independent mechanism — as an extra principle, with local isometry-style conditions.
- *Plain-English:* Just assuming the hidden factors are statistically unrelated isn't enough; you also need each factor to affect the output in its own independent way.
- *Applicability:* A2, A5 — motivates a Jacobian-orthogonality regularizer so each MORPHEUS programme perturbs expression through a distinct mechanism, sharpening interventional (A5) interpretability.
- *Novelty implication:* Offers a *reframe* MORPHEUS could adopt (mechanism-independence as the biological prior) — a differentiator from plain sparsity approaches.

---

### Group C — Cross-modal alignment, the modality gap & multimodal identifiability

**12. Understanding Contrastive Representation Learning through Alignment and Uniformity on the Hypersphere** (Wang & Isola, ICML 2020) — https://arxiv.org/abs/2005.10242
- *Takeaway:* Contrastive loss asymptotically optimizes two measurable properties — alignment (positive pairs close) and uniformity (features spread on the sphere) — that predict downstream performance.
- *Technical summary:* Decomposes InfoNCE into alignment + uniformity terms, shows each is separately optimizable, and that directly optimizing both matches or beats contrastive learning. Provides two closed-form metrics.
- *Plain-English:* Good representations do two things: matching items land near each other, and everything else is spread out evenly so the space is used fully.
- *Applicability:* A3, A4 — alignment/uniformity are ready-made *evaluation metrics* for MORPHEUS's NL↔biology grounding quality and for diagnosing modality collapse when encoding proteomics/CNV.
- *Novelty implication:* Supplies A3 evaluation instrumentation (measure grounding, don't just assert it) — strengthens the "and its evaluation" clause of A3.

**13. Mind the Gap: Understanding the Modality Gap in Multi-modal Contrastive Representation Learning** (Liang, Zhang, Kwon, Yeung & Zou, NeurIPS 2022) — https://arxiv.org/abs/2203.02053
- *Takeaway:* CLIP-style models leave the two modalities in *separate cones* of the shared space (the "modality gap"), driven by init + temperature, and the gap size affects downstream performance and fairness.
- *Technical summary:* Shows narrow-cone effect from random init persists through contrastive training; the temperature-controlled repulsion holds modalities apart. Manipulating the gap changes zero-shot accuracy and fairness.
- *Plain-English:* Even after training to match them, images and text live in two different neighborhoods of the shared map, not mixed together — and that separation matters.
- *Applicability:* A4 (critical) — if MORPHEUS encodes RNA + proteomics + NL text via contrastive pairing, expect a *biology-vs-language gap* that will silently degrade prompting. Design implication: measure and possibly close the gap (or exploit it) rather than assume a truly shared space.
- *Novelty implication:* Pre-empts a naive "unified shared representation" claim (A1/A4): the shared space is geometrically fractured by construction. MORPHEUS must address the gap explicitly — a concrete unaddressed risk.

**14. Identifiability Results for Multimodal Contrastive Learning** (Daunhawer, Bizeul, Palumbo, Marx & Vogt, ICLR 2023) — https://arxiv.org/abs/2303.09166
- *Takeaway:* Multimodal contrastive learning *block-identifies* the content (shared) latents even with modality-specific nuisances and nontrivial factor dependencies — a strictly more general setting than multi-view.
- *Technical summary:* Distinguishes multi-view (one generative mechanism) from multimodal (distinct mechanisms per modality); proves the shared block is recovered up to invertible transform even when only a subset of content is shared and factors are statistically dependent.
- *Plain-English:* Pairing two different kinds of data (say gene expression and text) provably recovers the information they share, even though each also carries its own private, irrelevant details.
- *Applicability:* A4 (core theoretical license) — justifies MORPHEUS *encoding* a modality specifically to recover the biology it *shares* with RNA, while private modality nuisances are quarantined. Decision rule: encode when a modality carries shared content; RAG-context it when it's mostly modality-specific.
- *Novelty implication:* This is the strongest existing theory *for* A4's shared-content claim — strengthens it but also means the identifiability is not MORPHEUS-novel. Novelty must be the biological encode-vs-context decision policy.

**15. Self-Supervised Learning with Data Augmentations Provably Isolates Content from Style** (von Kügelgen, Sharma, Gresele, Brendel, Schölkopf, Besserve & Locatello, NeurIPS 2021) — https://arxiv.org/abs/2106.04619
- *Takeaway:* Augmentation-based SSL provably identifies the *invariant content* partition up to invertible map, even under causal dependence between content and style.
- *Technical summary:* Models augmentations as soft interventions on "style"; proves content (augmentation-invariant) latents are block-identified in both generative and discriminative settings. Formal content/style split with dependency allowed.
- *Plain-English:* Training a model to ignore superficial changes (crops, color) makes it lock onto the stable "what it is" content — provably.
- *Applicability:* A4, A2 — batch/technical variation = "style"; biological programme = "content." Motivates augmentation/invariance objectives so MORPHEUS slots capture biology, not batch. Encode-vs-context: encode modalities that share *content*, treat batch-like signals as style to be invariant to.
- *Novelty implication:* Directly underwrites "identified biological content, batch-invariant" — a claim MORPHEUS wants for A2/A4. Cite as the mechanism; novelty is the biological content/style instantiation.

**16. Multi-View Causal Representation Learning with Partial Observability** (Yao, Xu, Lachapelle, Magliacane, Taslakian, Martius, von Kügelgen & Locatello, ICLR 2024) — https://arxiv.org/abs/2311.04056
- *Takeaway:* With one encoder per view and contrastive learning, the information shared across *any subset* of views is identifiable up to smooth bijection — with a graphical "identifiability algebra" saying exactly which latents are recoverable.
- *Technical summary:* Unifies multi-view nonlinear ICA, disentanglement, and CRL; each view is a nonlinear mixture of a *subset* of (possibly causally related) latents. Provides rules for which shared/partially-shared variables are identified given the view-overlap structure.
- *Plain-English:* When several partial "views" of a system overlap, you can work out precisely which hidden factors each combination of views lets you recover.
- *Applicability:* A4 (decision framework) — MORPHEUS's modalities are exactly partial views over shared biology; the identifiability algebra tells you *which programmes become addressable* when you add proteomics vs phospho vs CNV. Turns "which modality to encode" into a graphical calculation.
- *Novelty implication:* Provides a principled selection rule MORPHEUS can *adopt* for A4 — strengthens the framework; the encode/RAG policy can be derived from, not merely asserted against, this theory.

**17. Rethinking Minimal Sufficient Representation in Contrastive Learning** (Wang, Lin, Zhu, Zhou, Chen & Wang, CVPR 2022) — https://arxiv.org/abs/2203.07004
- *Takeaway:* Contrastive learning tends to a *minimal-sufficient* shared representation that discards non-shared info — risky when downstream tasks need that discarded info; add reconstruction/regularization to retain it.
- *Technical summary:* Shows the InfoNCE optimum keeps only inter-view shared information (minimal sufficient for the pretext task), causing overfitting to shared features; proposes regularizers to preserve task-relevant non-shared info.
- *Plain-English:* Matching two views can throw away everything that isn't shared — sometimes the very thing you later need — so you should deliberately keep some extra information.
- *Applicability:* A4 (caution) — if MORPHEUS aligns modalities purely contrastively, *modality-private biology* (e.g., a CNV signal absent from RNA) is discarded. Design implication: pair contrastive alignment with reconstruction so private-but-informative signals survive — informs when to *encode* vs *RAG-context*.
- *Novelty implication:* Pre-empts an over-strong "shared representation captures everything" claim; MORPHEUS should explicitly retain modality-private information — a design nuance reviewers will probe.

**18. Identifiable Multimodal Causal Representation Learning under Partial Latent Sharing** (Benhamza, Clausel & Tami, 2026) — https://arxiv.org/abs/2605.19135
- *Takeaway:* Extends multimodal identifiability to the causal setting where modalities share only a *subset* of latents, with guarantees that need no parametric latent distribution.
- *Technical summary:* Proves recovery of shared causal latent structure under partial sharing; introduces a Wasserstein-based alignment module and reports gains over prior multimodal CRL baselines.
- *Plain-English:* Even when two data types overlap only partly, you can still recover the shared *causal* factors linking them, without assuming a specific statistical form.
- *Applicability:* A4, A5 — the most current template for MORPHEUS's partial-sharing multimodal biology with causal semantics; supports counterfactual queries over shared causal slots (A5).
- *Novelty implication:* Very recent and close — a live novelty risk for any "identifiable multimodal causal slots" claim. MORPHEUS must position against this explicitly (biological grounding + prompting, not the identifiability result).

**19. Hierarchical Contrastive Learning for Multimodal Data** (2026) — https://arxiv.org/abs/2604.05462
- *Takeaway:* Learns *globally-shared, partially-shared, and modality-specific* representations jointly with proven identifiability of the hierarchical decomposition.
- *Technical summary:* A hierarchical contrastive objective factorizes latents into nested sharing levels; proves subspace-identifiability of each level and recovery of loading matrices, so the shared subspace driving cross-modal dependence is recoverable even when individual components are not.
- *Plain-English:* Instead of one shared space, it builds a layered map — what all modalities share, what pairs share, and what's unique to each — and proves each layer is recoverable.
- *Applicability:* A4 (architecture blueprint) — directly maps to MORPHEUS's need to separate pan-modality biology from modality-private signal; the hierarchy *is* the encode-vs-RAG decision surface made explicit.
- *Novelty implication:* Offers a concrete architecture MORPHEUS could build on; reframes A4 from binary (encode/RAG) to a *hierarchy of sharing* — a stronger, defensible framing but not MORPHEUS-original.

**20. MVEB: Self-Supervised Learning with Multi-View Entropy Bottleneck** (Wen, Chen, Li, Zhu et al., IEEE TPAMI 2024) — https://arxiv.org/abs/2403.19078
- *Takeaway:* Operationalizes "minimal sufficient" via a multi-view entropy bottleneck: maximize view agreement while maximizing embedding entropy to drop superfluous info.
- *Technical summary:* Objective = agreement between two views + differential entropy of the embedding (uniformity-like) with a tractable entropy estimator; yields minimal-sufficient representations with strong linear-probe results.
- *Plain-English:* Keep exactly what two views agree on and nothing more, by rewarding agreement and penalizing wasted capacity.
- *Applicability:* A4, A3 — a practical loss for MORPHEUS's shared latent that controls how much modality-private info is kept; the entropy term doubles as a uniformity diagnostic (A3 evaluation).
- *Novelty implication:* Tooling for A4/A3; supports the minimal-sufficient framing while (per #17) flagging its risks — MORPHEUS should tune the bottleneck to biology, not vision defaults.

---

### Group D — Causal representation learning, interventions & counterfactual geometry (A5-heavy)

**21. Towards Causal Representation Learning** (Schölkopf, Locatello, Bauer, Ke, Kalchbrenner, Goyal & Bengio, Proc. IEEE 2021) — https://arxiv.org/abs/2102.11107
- *Takeaway:* Position paper defining causal representation learning — recovering high-level causal variables from raw data — as the route to transfer, intervention, and counterfactual reasoning.
- *Technical summary:* Frames the ICM (independent causal mechanisms) principle, SCM-based representations, and why disentanglement is a special case of causal factorization; motivates interventions/multi-environment data as the identifiability source.
- *Plain-English:* To generalize and answer "what if" questions, models should discover the actual causal building blocks of the world, not just correlated features.
- *Applicability:* A5 (foundational), A2 — the manifesto behind "drug/perturbation as a query, not a retrained classifier." Grounds MORPHEUS's causal-geometry ambition in an authoritative reference.
- *Novelty implication:* Provides vocabulary and legitimacy for A5; MORPHEUS's novelty is *executing* CRL on multimodal omics with a prompting interface, which this paper only motivates.

**22. Interventional Causal Representation Learning** (Ahuja, Mahajan, Wang & Bengio, ICML 2023) — https://arxiv.org/abs/2209.11924
- *Takeaway:* Access to *interventional* data provably identifies latent causal factors, with guarantees under weak assumptions and even unknown intervention targets in some regimes.
- *Technical summary:* Shows single-node interventions enable block-affine/permutation identification of latents without distributional assumptions on the factors; characterizes what observational vs interventional data each buy.
- *Plain-English:* If you can poke the system (intervene) and watch what changes, you can pin down the true underlying causes — something you can't do from passive observation alone.
- *Applicability:* A5 (core) — CRISPR/drug perturbations are exactly the interventions this theory needs; grounds MORPHEUS's claim that perturbation data makes slots *causally* addressable and queryable.
- *Novelty implication:* Strengthens A5's theoretical basis; the identifiability-from-intervention result is established, so MORPHEUS claims the *biological perturbation-as-prompt interface*, not the theorem.

**23. Linear Causal Disentanglement via Interventions** (Squires, Seigal, Bhate & Uhler, ICML 2023) — https://arxiv.org/abs/2211.16467
- *Takeaway:* For linear latent SCMs under linear mixing, *one perfect intervention per latent node* is necessary and sufficient for identifiability.
- *Technical summary:* Proves interventions are necessary (a missed node ⇒ indistinguishable models) and that a single intervention per node suffices, via a partial-order-generalized RQ matrix decomposition.
- *Plain-English:* To untangle a chain of causes, you need to have poked each cause at least once; do that and the whole causal structure becomes recoverable.
- *Applicability:* A5 (design budget) — tells MORPHEUS the *perturbation coverage* required to claim identified causal programmes: roughly one intervention per programme. A concrete experimental-design prescription.
- *Novelty implication:* Pre-empts "we identify causal programmes" claims made without adequate perturbation coverage — MORPHEUS must report coverage against this necessity result.

**24. Learning Linear Causal Representations from Interventions under General Nonlinear Mixing** (Buchholz, Rajendran, Rosenfeld, Aragam, Schölkopf & Ravikumar, NeurIPS 2023, oral) — https://arxiv.org/abs/2306.02235
- *Takeaway:* Identifiability of latent causal variables holds even under *nonlinear* mixing (Gaussian latents), using a contrastive algorithm exploiting precision-matrix geometry — no paired interventions needed.
- *Technical summary:* Extends interventional CRL beyond linear decoders; recovers latents from single-node interventions under general nonlinear observation maps via geometric structure in per-environment precision matrices.
- *Plain-English:* Even when the mapping from causes to data is complicated and nonlinear, poking one cause at a time is enough to recover them.
- *Applicability:* A5, A4 — realistic since omics readouts are nonlinear functions of latent programmes; supports MORPHEUS using nonlinear encoders while retaining causal identifiability from perturbations.
- *Novelty implication:* Closes the "but our decoder is nonlinear" objection to A5 — strengthens feasibility; again the theorem is prior art.

**25. Identifiability Guarantees for Causal Disentanglement from Soft Interventions** (Zhang, Squires, Greenewald, Srivastava, Shanmugam & Uhler, NeurIPS 2023) — https://arxiv.org/abs/2307.06250
- *Takeaway:* *Soft* (imperfect) interventions still give identifiability under a generalized faithfulness condition, even with unobserved causal variables; applied to combinatorial perturbation prediction in genomics.
- *Technical summary:* Variational-Bayes method recovering latent SCM from soft interventions; demonstrates predicting effects of unseen perturbation *combinations* on real genomic data.
- *Plain-English:* You don't need to fully clamp a gene — even partial perturbations let you recover the causal structure and predict combinations you've never seen.
- *Applicability:* A5 (core), A3 — most real drug/CRISPR perturbations are soft; directly supports MORPHEUS's counterfactual/combinatorial "what if I perturb X" queries and their evaluation on held-out combinations.
- *Novelty implication:* Very close prior for "predict unseen perturbation combinations from an identified causal latent" — a live novelty risk for A5's combinatorial-counterfactual claim. Differentiate via the NL-query interface and cross-modal grounding.

**26. General Identifiability and Achievability for Causal Representation Learning** (Varıcı, Acartürk, Shanmugam & Tajer, AISTATS 2024, oral) — https://arxiv.org/abs/2310.15450
- *Takeaway:* Two hard uncoupled interventions per node give full identifiability *without* the faithfulness assumptions prior work needed; provides a constructive score-based algorithm.
- *Technical summary:* Uses score-function variation across environments to invert the mixing and recover the latent DAG + variables; sharpens the sufficient-intervention count and removes faithfulness when observational data is available.
- *Plain-English:* With a couple of clean pokes per hidden cause, you can recover the whole causal picture, and you don't need extra "niceness" assumptions.
- *Applicability:* A5 — refines the perturbation budget (two interventions/node for the assumption-light regime) MORPHEUS should target for strong identifiability guarantees.
- *Novelty implication:* Tightens the experimental-design story for A5; MORPHEUS should quote its coverage against these bounds rather than claim a new guarantee.

**27. Causal Component Analysis (CauCA)** (Wendong, Kekić, von Kügelgen, Buchholz, Besserve, Gresele & Schölkopf, NeurIPS 2023) — https://arxiv.org/abs/2305.17225
- *Takeaway:* Generalizes ICA to *causally dependent* latents given a known causal graph; interventional datasets identify both the unmixing and the mechanisms.
- *Technical summary:* Assumes the latent causal graph is known, drops the independence assumption of ICA, and proves identifiability from interventions via a normalizing-flow estimator of unmixing + mechanisms.
- *Plain-English:* Classic ICA assumes hidden factors are independent; this handles factors that cause each other, as long as you know the wiring diagram and can intervene.
- *Applicability:* A5, A2 — pathway knowledge gives MORPHEUS a *prior causal graph* over programmes; CauCA shows how to exploit it to identify causally-linked (not independent) slots — realistic for biology where pathways interact.
- *Novelty implication:* Reframes A2 away from "independent programmes" (biologically false) toward "known-graph causal programmes" — a differentiating, biology-honest stance MORPHEUS should adopt.

**28. Causal Representation Learning Made Identifiable by Grouping of Observational Variables** (Morioka & Hyvärinen, ICML 2024) — https://arxiv.org/abs/2310.15709
- *Takeaway:* A *known grouping* of observed variables gives CRL identifiability with **no** temporal structure, interventions, or weak supervision required.
- *Technical summary:* Exploits an assumed grouping (which observables belong together) as the structural bias; proves identifiability of the causal latents and gives a self-supervised estimation method.
- *Plain-English:* If you know which measurements naturally cluster together, that alone is enough to recover the causal factors — no perturbations needed.
- *Applicability:* A4, A2 — gene-module / pathway-membership groupings are exactly this bias; lets MORPHEUS claim identified slots even *without* perturbation data by using pathway groupings as the structural prior.
- *Novelty implication:* Gives A2/A4 an intervention-free fallback — a valuable hedge when perturbation data is absent; the grouping principle is prior art but its biological instantiation is open.

**29. Weakly Supervised Causal Representation Learning** (Brehmer, de Haan, Lippe & Cohen, NeurIPS 2022) — https://arxiv.org/abs/2203.16437
- *Takeaway:* *Paired* pre-/post-intervention samples (without labels) identify causal representations and the causal graph.
- *Technical summary:* Introduces implicit latent causal models trained with VAEs on before/after-intervention image pairs; proves identifiability of causal variables and structure from the pairing alone.
- *Plain-English:* Show the model the same scene before and after a change, and it learns both the hidden causes and how they connect — no annotations needed.
- *Applicability:* A5 — matched control/perturbed cell pairs (or pseudo-pairing via optimal transport) are the biological analog; supports learning MORPHEUS's causal slots from perturbation *contrasts*.
- *Novelty implication:* Supports A5's "perturbation contrast as supervision" but pre-empts the pairing mechanism; MORPHEUS's novelty is handling *unpaired* single-cell perturbations (cells are destroyed on measurement).

**30. CITRIS: Causal Identifiability from Temporal Intervened Sequences** (Lippe, Magliacane, Löwe, Asano, Cohen & Gavves, ICML 2022) — https://arxiv.org/abs/2202.03169
- *Takeaway:* Temporal sequences with *known intervention targets* identify scalar and multidimensional causal factors, extensible to frozen pretrained encoders via normalizing flows.
- *Technical summary:* VAE/flow framework exploiting time + intervention-target info to disentangle causal variables in video; handles multidimensional factors and plugs onto pretrained autoencoders.
- *Plain-English:* Watching a system over time while knowing what was changed at each step lets you recover the causal variables — even bolting onto an existing pretrained model.
- *Applicability:* A5, A4 (frozen-trunk) — the normalizing-flow-on-frozen-encoder trick is exactly MORPHEUS's "frozen-trunk plug-in": add a causal head onto a frozen foundation model to get interventional slots without retraining.
- *Novelty implication:* Provides a concrete mechanism for A4's frozen-trunk plug-in + A5's interventional head — MORPHEUS can build on this; novelty is the biological + NL-prompt wrapper.

**31. Learning Causal Representations of Single Cells via Sparse Mechanism Shift Modeling** (Lopez, Tagasovska, Ra, Cho, Pritchard & Regev, CLeaR 2023) — https://arxiv.org/abs/2211.03553
- *Takeaway:* Treats genetic/chemical perturbations as *stochastic interventions on a sparse unknown subset of latents*, giving a biology-native CRL model that transfers to unseen perturbations.
- *Technical summary:* Deep generative model with sparse mechanism-shift priors over latent variables; identifies which latents a perturbation targets and outperforms baselines on real perturbation transfer tasks.
- *Plain-English:* Each drug or gene knockout nudges just a few hidden cell-state factors; modeling that sparsity lets the model generalize to perturbations it hasn't seen.
- *Applicability:* A5 (direct biological precedent), A2 — the closest existing instantiation of MORPHEUS's "perturbation ⇒ sparse programme shift ⇒ addressable slot"; validates the whole A5 premise on real Perturb-seq-style data.
- *Novelty implication:* **Highest-overlap prior for MORPHEUS's core A5 story.** MORPHEUS must clearly differentiate: NL-promptable interface, cross-modal grounding, emergent-knowledge elicitation — not the sparse-mechanism-shift single-cell model itself.

**32. Self-Supervised Contrastive Learning Performs Non-Linear System Identification** (González Laiz, Schmidt & Schneider, ICLR 2025) — https://arxiv.org/abs/2410.14673
- *Takeaway:* Contrastive SSL with a dynamics prior ("dynamics contrastive learning") recovers latent linear/switching/nonlinear dynamical systems under nonlinear observation maps, with identifiability guarantees.
- *Technical summary:* Connects SSL to control-theoretic system identification; proves recovery of latent dynamics and validates across system classes, extending contrastive-ICA theory to temporal/dynamical latents.
- *Plain-English:* Contrastive learning on time series can reverse-engineer the hidden dynamical system generating the data, not just static factors.
- *Applicability:* A5, A2 — for trajectory/pseudotime data, MORPHEUS's latent programmes can be treated as a recoverable *dynamical system*, enabling counterfactual "roll the dynamics forward under perturbation" queries.
- *Novelty implication:* Extends identifiability to dynamics — supports A5 counterfactual-trajectory queries; a reframe MORPHEUS can use to go beyond static-slot claims.

**33. Unifying Causal Representation Learning with the Invariance Principle** (Yao, Chen, von Kügelgen, Locatello et al., 2024/2025) — https://arxiv.org/abs/2409.02772
- *Takeaway:* Recasts many CRL identifiability results as instances of a single *invariance* principle — align representations to invariances implied by the data structure (views, time, interventions).
- *Technical summary:* Shows multi-view, temporal, and interventional CRL identifiability all follow from enforcing the right invariances; unifies the zoo of prior results and clarifies which data assumption yields which guarantee.
- *Plain-English:* All the different ways to recover hidden causes boil down to one idea — force the representation to respect the invariances baked into how the data was collected.
- *Applicability:* A1, A2 — gives MORPHEUS a *single organizing principle* for choosing which structural bias (views vs time vs perturbation) to enforce per task, unifying the routing story.
- *Novelty implication:* Reframes the whole lane — MORPHEUS can present its design choices as "which invariance to enforce," a clean narrative, but must acknowledge the unification is prior art.

---

### Synthesis for the rebase

- **The identifiability ladder MORPHEUS is climbing is well-mapped.** Auxiliary conditioning (iVAE), multi-view/multimodal pairing (Daunhawer, Yao), content/style split (von Kügelgen), sparsity (Lachapelle), and interventions (Ahuja, Squires, Lopez) each *buy* identifiability at a known price. MORPHEUS's defensible novelty is not any single primitive but (a) the NL-promptable interface over identified slots (A1), (b) biological grounding + emergent-knowledge elicitation and its measurement (A3), and (c) the encode-vs-RAG decision policy derived from partial-observability theory (A4).
- **Three papers are direct novelty threats to A5/A2:** Lopez et al. 2023 (single-cell sparse mechanism shift), Zhang et al. 2023 (soft-intervention combinatorial genomics), and Benhamza et al. 2026 (identifiable multimodal causal, partial sharing). MORPHEUS must cite and differentiate against these explicitly.
- **Two hard guardrails:** Locatello (no identifiability without inductive bias — name yours) and the modality gap (Liang — the "unified" space is geometrically fractured). Both must be addressed head-on, not assumed away.
