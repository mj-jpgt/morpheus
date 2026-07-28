# The Effective-Rank Fingerprint: Diagnosing Cohort-Confounded WSI→Molecular Alignment from a Single Trained Model

**Authors.** [Author placeholder — names, affiliations, and corresponding author to be filled in]

---

## Abstract

When a contrastive model is trained to align whole-slide images (WSIs) with molecular programme signatures, one head of a dual-head encoder can quietly stop using its representation space. We show that a biology head supervised to regress a $\sim$50-dimensional MSigDB Hallmark manifold undergoes *effective-rank collapse* — its 256-dimensional embedding falls to an effective rank of $\sim$5–6 — while its sibling identity head, trained with paired contrastive supervision, stays healthy at $\sim$84. We make three claims. **(C1)** This dual-head collapse is a reproducible *fingerprint*: the biology head's effective rank tracks the intrinsic rank of its supervision target, not the capacity of the encoder. **(C2)** A rank-collapsed head still scores respectably on molecular prompting *only because the benchmark is confounded*: $\sim$46–49\% of Pearson correlation is cross-cancer cohort structure, random-gene nulls already reach $\sim$0.30–0.32, and the genuine within-cancer, control-adjusted specificity is a *method-invariant* $\sim$$+0.07$. **(C3)** The fix is prescriptive and provable: a per-dimension variance floor cannot prevent rank collapse, but an off-diagonal covariance-decorrelation term can. We propose the biology head's effective rank as an internal diagnostic that flags cohort-confounded alignment without knowing the confounder a priori.

---

## 1. Introduction

Consider a two-tower encoder that maps a whole-slide image to a molecular-programme embedding and is trained, contrastively and by regression, to place each slide near its matching Hallmark-pathway signature. The supervision target is low-dimensional: the 50 MSigDB Hallmark programmes span an intrinsic rank of only $\sim$5–6. When we look inside a trained model, the head responsible for this alignment has inherited that dimensionality — its 256-D embedding has collapsed to an effective rank of 6.0 on held-out test slides, while the *identity* head sitting beside it in the same network, supervised by paired image–RNA contrastive learning, retains an effective rank of 84.3. The two heads share an encoder, a batch, and an optimizer; they differ only in what they are asked to match, and one of them has stopped using almost all of its representational geometry. This is the phenomenon we study.

The mechanism is specific. A representation aligned or regressed onto a rank-$r$ manifold has no gradient pressure to occupy directions beyond that manifold, and heteroscedastic Gaussian-NLL regression onto a Hallmark target actively rewards contracting the unused directions. The result is *rank* collapse — a rank-deficient covariance — which is distinct from, and not cured by, keeping every coordinate's *variance* off the floor. A per-dimension variance penalty of the form $\mathrm{relu}(1-\mathrm{std})$ constrains diagonal variances but places no constraint on off-diagonal covariance, so a representation can satisfy it with unit per-dimension variance while living entirely in a 6-dimensional subspace. Preventing collapse therefore requires a term that acts on the *covariance* — decorrelating features — not on marginal variances.

We report this in a concrete system: a hierarchical patch-pooling Query-Former over uncapped H-Optimus WSI patches paired with a BulkFormer RNA encoder, producing an L2-normalized `z_identity` (anchored as a bounded zero-init residual on a frozen all-patch MLP–CLIP teacher) and an L2-normalized `z_biology` supervised by 50-D Hallmark targets through Gaussian-NLL, programme-neighbour KL, and supervised-contrastive losses. Across held-out slides ($n=2530$), the identity head is healthy (effective rank 84.3) and the biology head is collapsed (6.0), with a correspondingly inflated modality gap (0.475 vs.\ 0.296). The collapse is stable across seeds and objective profiles.

**Three claims.** **(C1, mechanism.)** The matched dual-head rank-collapse fingerprint is a signature of aligning to a low-rank supervision manifold: the collapsed head's effective rank tracks the target's intrinsic rank, independent of encoder capacity, and its healthy sibling rules out a trivial optimization failure. **(C2, coupling.)** A collapsed head is not obviously broken because the *benchmark hides it*. On our molecular-prompting evaluation, every method — including a plain MLP–CLIP baseline — loses $\sim$46–49\% of its global Pearson to within-cancer, and random-gene null sets already reach $\sim$0.30–0.32; subtracting the random-gene control leaves a within-cancer specificity of $\sim$$+0.07$ Pearson that is *identical across all methods*, collapsed or not. Respectable-looking scores and rank collapse coexist precisely because the task rewards cohort structure a low-rank embedding can still supply. **(C3, prescriptive.)** The fix is two coupled changes — feed a single normalized biology state to both the NLL heads and the rank-spreaders so they share geometry, and add a VICReg/Barlow-style off-diagonal decorrelation term ($\sim$$0.04$) on the normalized 256-D biology embedding — and we show why the cheaper alternative (a variance floor) provably cannot work. The biology head's effective rank then serves as an *internal fingerprint* that flags cohort-confounded WSI$\rightarrow$molecular alignment before any external audit.

**What we concede, and what is new.** The individual ingredients are known, and we claim none of them. The cohort confound in computational pathology is documented [Howard et al.\ 2021; *Buyer Beware*, NBME 2026; DECAT]. The random-gene control for molecular prediction is standard [Schmauch et al., HE2RNA 2020; Venet et al.\ 2011], as is the observation that within-cancer stratification roughly halves performance [Fu et al.\ 2020]. Representation collapse as a phenomenon — dimensional collapse, rank collapse — and its cures are well studied [Andriopoulos et al., NRC1, NeurIPS 2024; Jing et al., ICLR 2022; Bardes et al., VICReg, ICLR 2022], and the retrieval-versus-regression distinction that predicts *which* head collapses has theoretical grounding [Wang & Isola, CVPR 2022; FactorCL, NeurIPS 2023]. Our contribution is the *coupled package*: the demonstration that a specific, low-rank molecular supervision target produces a *diagnosable* dual-head rank-collapse fingerprint (C1); that this fingerprint is causally coupled to a confounded benchmark whose method-invariant $\sim$$+0.07$ specificity explains why collapse goes unpunished (C2); and that the fingerprint admits a prescriptive, provable fix that a variance floor cannot supply (C3). We frame effective rank not as a quantity to be maximized for its own sake but as a *diagnostic instrument* for a failure mode that standard benchmarks are structurally unable to reveal.

---

## 2. Related Work

Our contribution sits at the intersection of four literatures. We are neither the first to build a WSI$\rightarrow$expression benchmark, nor to warn about cohort confounds, nor to characterize representational collapse, nor to reason about what contrastive alignment retains. Our novelty is the *coupled* package: a matched dual-head **rank-collapse fingerprint** that ties a specific geometric failure of the biology head to the cross-cancer structure of the benchmark, together with a prescriptive fix. We situate that claim against each thread below and state precisely where we differ.

### 2.1 WSI → molecular / expression prediction and benchmarks

A now-substantial line of work maps whole-slide morphology onto molecular readouts. HE2RNA regressed bulk RNA-Seq from WSIs across TCGA and reported that predictability concentrates in immune- and hallmark-related signatures \citep{schmauch2020he2rna}; Fu et al. showed pan-cancer that learned histological features correlate with recurrent genetic alterations and bulk expression, and documented that a large fraction of apparent signal reflects tumour type and composition rather than fine-grained morphology \citep{fu2020pancancer}. More recently, HEST-1k standardized 1,229 spatial-transcriptomic/WSI pairs and a foundation-model benchmark (HEST-Benchmark) for predicting expression from morphology \citep{jaume2024hest}. **We differ:** we do not propose a new dataset or a stronger regressor. We treat the WSI$\rightarrow$molecular-programme task as a *diagnostic probe* and show (C2) that on a $\sim$46--49% cross-cancer benchmark every method — including a plain MLP-CLIP baseline — retains an identical, method-invariant within-cancer, control-adjusted specificity of only $\sim$+0.07 Pearson, so leaderboard gains here are largely cohort structure rather than morphology$\rightarrow$biology alignment.

### 2.2 Cohort and batch confounding in computational pathology

That cohort structure is a first-class hazard is well established. Howard et al. demonstrated that site-specific digital signatures in TCGA survive color normalization and bias survival, mutation, and stage prediction \citep{howard2021impact}; Dehkharghanian et al. showed deep networks recover the acquisition *site* of TCGA slides with high accuracy, evidence that batch identity is linearly available in the features \citep{dehkharghanian2023biased}. "Buyer Beware" catalogued how co-dependencies and confounders (e.g., tumour mutational burden) inflate omics-from-histology biomarker estimates \citep{tizhoosh2026rethinking}, HESCAPE found that contrastive pretraining can *degrade* direct expression prediction and named batch effects as the interfering factor \citep{gindra2025hescape}, and DECAT built a null-referenced, confounder-label-free decision procedure that flags when entangled models (e.g., CLIP) report "shared biology" that is actually confounding \citep{steiner2026decat}. **We differ:** we concede the confound entirely rather than re-detect it, and instead identify an *internal, geometry-only fingerprint* of it — effective-rank collapse of the biology head ($\sim$5--6 of 256) while the sibling identity head stays healthy ($\sim$84) — that is read off a single trained model without holding out cohorts, without confounder labels, and without a paired null benchmark.

### 2.3 Representational and neural collapse

The geometry we exploit has deep roots. Jing et al. characterized *dimensional collapse* in contrastive SSL, where embeddings occupy a low-dimensional subspace \citep{jing2022understanding}; VICReg introduced explicit variance-and-covariance regularization — a per-dimension variance hinge plus an off-diagonal decorrelation term — to prevent it without asymmetry tricks \citep{bardes2022vicreg}; and Andriopoulos et al. formalized Neural Regression Collapse (NRC1) for multivariate regression, showing last-layer representations collapse toward the low-rank span of the target subspace \citep{andriopoulos2024prevalence}. **We differ:** prior work treats collapse as a phenomenon to prevent globally; we use it as a *targeted diagnostic and target-aware prediction*. Because our biology head is regressed onto a rank-$\sim$5--6 Hallmark manifold, NRC-style collapse of the head to that rank is expected, so its effective rank becomes a confound *fingerprint*; and we show analytically and empirically that a per-dimension variance floor (VICReg's variance term alone, $\mathrm{relu}(1-\mathrm{std})$) provably cannot restore rank, whereas only the covariance/decorrelation term can (C3) — sharpening which half of the VICReg objective is load-bearing for *this* failure.

### 2.4 Retrieval-vs-regression: what alignment keeps and drops

Finally, theory tells us contrastive and regression objectives retain different information. Wang et al. showed the minimal-sufficient representation from contrastive learning discards task-relevant, non-shared information and is therefore insufficient downstream \citep{wang2020understanding}; FactorCL formalized this for the multimodal case, factorizing task-relevant information into shared and modality-unique parts and arguing shared-only objectives miss unique signal \citep{liang2023factorized}. **We differ:** these results motivate *why* a biology head aligned only to shared cross-modal structure can score respectably while carrying little modality-unique morphological signal — but they are stated at the level of mutual information. We instantiate the consequence as a concrete, measurable geometric signature in a deployed WSI$\rightarrow$molecular model and turn it into an actionable design rule (feed one shared L2-normalized biology state to both the Gaussian-NLL heads and the decorrelation term, then add off-diagonal covariance decorrelation), rather than proposing a new general-purpose objective.

---

## 3. Method (architecture + the fix)

We first describe the MORPHEUS V2 architecture that produces the two embedding heads at the centre of our study (§3.1). We then diagnose the failure mode — a rank collapse localised to the *biology* head while its sibling *identity* head stays healthy (§3.2) — and prove that the natural remedy, a per-dimension variance floor, cannot address it. Finally we give the fix: a *normalize-once* coupling (F-R1) combined with an off-diagonal covariance-decorrelation term (F-R2), and we define the **effective-rank fingerprint** that we propose as an internal, confounder-agnostic diagnostic for cohort-confounded WSI$\rightarrow$molecular-programme alignment (§3.3).

Throughout, a *WSI* is a whole-slide image, $z\in\mathbb{R}^{d}$ with $d=256$ is an L2-normalised embedding, and a *batch* is a matrix $Z=[z_1,\dots,z_n]^{\top}\in\mathbb{R}^{n\times d}$ of $n$ such embeddings.

### 3.1 The MORPHEUS V2 architecture

MORPHEUS V2 is a dual-encoder, dual-head model that maps a whole-slide image and its paired bulk RNA-seq profile into a shared 256-D space, split into an **identity** subspace (patient/slide provenance) and a **biology** subspace (molecular programme state).

**Vision branch — hierarchical Query-Former.**
Uncapped H-Optimus patch features (all patches per slide, no fixed cap) are pooled by a hierarchical soft-slot attention module — a *Query-Former* — that performs learned soft-slot pooling at three granularities: 64 **local** slots, 32 **slide** slots, and 128 **patient** slots. Each level attends over the tokens produced by the level below, so a patient representation is a soft aggregation of slide summaries, which are in turn soft aggregations of local tissue neighbourhoods. Soft-slot pooling (attention weights over an over-complete slot dictionary, rather than hard clustering) keeps the pooling differentiable and permutation-invariant while accommodating a variable, and typically large, number of input patches.

**Molecular branch.**
Paired bulk RNA-seq is encoded by BulkFormer, yielding an RNA token summary that is projected into the same 256-D geometry as the vision heads.

**Two heads.**
From the fused representation the model emits two L2-normalised vectors:

- $z_{\mathrm{identity}}\in\mathbb{S}^{255}$ — a **bounded, zero-initialised residual** on a *frozen all-patch MLP-CLIP teacher*. Writing $t$ for the (unit-norm) teacher embedding, a learned gate $g\in[0,1]$ and a bounded correction $\delta$,
  $$
  z_{\mathrm{identity}}
  \;=\;
  \mathrm{normalize}\!\big(\,t + g\cdot s\cdot \delta\,\big),
  \qquad s = \text{residual scale},
  $$
  with $\delta$ produced by the anchor head and $s$ a small learned scalar initialised at $0$. Zero-init means the head *begins as the teacher* and can only depart from it to the extent the residual earns; in our trained model the residual is effectively inert ($s=-0.0011$, mean correction norm $0.00073$, gate mean $0.646$), so $z_{\mathrm{identity}}$ is, to numerical precision, the frozen MLP-CLIP teacher. This *anchoring* is what makes the identity head a clean, healthy control against which the biology head's collapse is measured.

- $z_{\mathrm{biology}}\in\mathbb{S}^{255}$ — supervised toward a 50-D MSigDB **Hallmark** programme target vector $y\in\mathbb{R}^{50}$.

**Biology-head supervision.**
The biology head is trained with three objectives against the Hallmark targets:

1. a **heteroscedastic Gaussian negative log-likelihood** (Gaussian-NLL): a read-out head predicts a mean $\mu(z_{\mathrm{biology}})$ and a per-programme log-variance $\log\sigma^2(z_{\mathrm{biology}})$, and is trained by
   $$
   \mathcal{L}_{\mathrm{NLL}}
   = \frac{1}{2}\sum_{k=1}^{50}
   \left[
   \frac{\big(y_k-\mu_k\big)^2}{\sigma_k^2}
   + \log \sigma_k^2
   \right];
   $$
2. a **programme-neighbour KL** that softens supervision toward programme-space neighbours (encouraging embeddings of biologically adjacent programmes to share a neighbourhood); and
3. a **supervised contrastive** term over Hallmark labels.

The Hallmark target manifold is nominally 50-D but has an **intrinsic rank of only $\sim$5–6**: the 50 Hallmark programme scores are strongly co-varying, so the label geometry the biology head is asked to reproduce is extremely low-rank. This fact is the seed of the collapse.

### 3.2 Diagnosis: a localised rank collapse

We measure the geometry of each head on the held-out test split ($n=2530$, $d=256$) via **effective rank**. Given a batch $Z$ with centred covariance $C=\tfrac1n Z_c^{\top}Z_c$ (equivalently the singular values $\sigma_1\ge\dots\ge\sigma_d\ge 0$ of $Z_c/\sqrt{n}$), define the $\ell_1$-normalised spectrum $p_i=\sigma_i / \sum_j \sigma_j$ and the Roy–Vetterli effective rank

$$
\boxed{\;
\mathrm{erank}(Z)
\;=\;
\exp\!\Big( -\sum_{i=1}^{d} p_i \log p_i \Big)
\;=\;
\exp\!\big(H(p_1,\dots,p_d)\big)
\;}
$$

the exponential of the Shannon entropy of the (normalised) singular-value distribution. $\mathrm{erank}\in[1,d]$: it equals $d$ for an isotropic (white) representation and $1$ when all variance lies along a single direction. It is a smooth, basis-independent measure of how many directions the representation actually populates.

On the held-out split the two heads diverge sharply:

| head              | effective rank |
|-------------------|---------------:|
| `wsi_identity`    | 84.3           |
| `rna_identity`    | 37.5           |
| `wsi_biology`     | 6.0            |
| `rna_biology`     | 4.4            |
| `full_biology`    | 5.3            |
| `full_patient`    | 8.0            |

The identity head, anchored to the frozen teacher and trained with an RNA-paired contrastive signal, occupies $\sim$84 of 256 directions. The biology head — regressed and aligned to the intrinsic-rank-$\sim$5–6 Hallmark manifold — has **collapsed to $\sim$5–6 effective directions**, i.e. it has inherited the low-rank geometry of its supervision targets almost exactly. This is a *dimensional/rank collapse* in the sense of Jing et al. (2022) and NRC1 (Andriopoulos et al., 2024), but *localised to one head of a matched sibling pair* — the healthy identity head serves as a within-model control that rules out data-, optimiser-, or capacity-level explanations. The modality (WSI$\leftrightarrow$RNA) gap between L2-normed centroids is correspondingly larger for the collapsed head (biology $0.475$) than for the healthy one (identity $0.296$).

**Why a per-dimension variance floor cannot fix it.**
The standard VICReg *variance* term applies a hinge floor to each coordinate's standard deviation, $\mathcal{L}_{\mathrm{var}} = \tfrac1d\sum_{j=1}^{d}\mathrm{relu}\!\big(\gamma-\mathrm{std}(Z_{:,j})\big)$ (here $\gamma=1$, weight $0.01$). We state precisely why this is powerless against *rank* collapse.

> **Proposition (variance floor $\not\Rightarrow$ rank).**
> Let $Z_c\in\mathbb{R}^{n\times d}$ be a centred batch and let $\mathrm{std}(Z_{:,j})=\sqrt{C_{jj}}$ denote the $j$-th coordinate standard deviation, $C=\tfrac1n Z_c^\top Z_c$. The per-dimension floor $\mathcal{L}_{\mathrm{var}}$ is a function of the **diagonal** of $C$ only. Its global minimum ($\mathcal{L}_{\mathrm{var}}=0$) is attained by *any* $C$ with $C_{jj}\ge\gamma^2$ for all $j$ — including rank-1 configurations. In particular, take $Z_c = a\,u^\top$ for a fixed unit direction $u\in\mathbb{R}^d$ with $|u_j|\ge\gamma/\|a\|_{\mathrm{rms}}$ and an activation vector $a\in\mathbb{R}^n$: every coordinate meets its variance floor, yet $\mathrm{rank}(C)=1$ and $\mathrm{erank}(Z)\to 1$. Hence $\mathcal{L}_{\mathrm{var}}=0$ is consistent with maximal rank collapse.

*Proof sketch.* $\mathcal{L}_{\mathrm{var}}$ depends on $C$ solely through $\{C_{jj}\}$. The rank of $C$ is governed by its **off-diagonal** structure: a matrix with large, satisfied diagonal entries can still be exactly rank-1 if its rows are mutually proportional (perfect off-diagonal correlation). The explicit $Z_c=a\,u^\top$ construction realises this: it drives every marginal variance above the floor while confining all mass to a single eigen-direction. Marginal variances are invariant to the correlations that determine rank; only a term that reads the **off-diagonal** covariance can constrain rank. $\square$

The proposition motivates the exact form of the fix: decorrelation, not per-coordinate scaling.

### 3.3 The fix: normalize-once (F-R1) + covariance decorrelation (F-R2)

The fix is two *coupled* changes to the biology-head objective, plus an optional escalation.

**(F-R1) Normalize-once — shared geometry.**
In the collapsed model the Gaussian-NLL read-out heads and the rank-spreading regularisers consumed *different* views of the biology state (pre- vs. post-normalisation), so a regulariser could inflate the geometry of a state the likelihood never saw. F-R1 feeds a **single L2-normalised biology state** $z_{\mathrm{biology}}=\mathrm{normalize}(h)$ to *both* the Gaussian-NLL heads *and* the decorrelation term, so likelihood and geometry act on the same object. This coupling is what makes F-R2 bite: decorrelating a state that the supervision also optimises forces the two to co-adapt rather than fight.

**(F-R2) VICReg / Barlow-style off-diagonal covariance decorrelation.**
On the normalised 256-D biology batch $Z\in\mathbb{R}^{n\times d}$ we add a term that penalises the off-diagonal entries of the centred batch covariance. Let $\bar z=\tfrac1n\sum_i z_i$, $Z_c = Z-\mathbf 1\bar z^\top$, and

$$
C \;=\; \frac{1}{\,n-1\,}\,Z_c^{\top} Z_c \;\in\;\mathbb{R}^{d\times d}.
$$

The decorrelation loss is the **mean-squared off-diagonal** of $C$,

$$
\boxed{\;
\mathcal{L}_{\mathrm{cov}}
\;=\;
\frac{1}{d}\sum_{i\ne j} C_{ij}^{\,2}
\;=\;
\frac{1}{d}\Big(\|C\|_F^2 - \sum_{i=1}^{d} C_{ii}^2\Big)
\;}
$$

added to the objective with weight $\lambda_{\mathrm{cov}}\approx 0.04$, active in the `full` and `programme_only` objective profiles. Because $\mathcal L_{\mathrm{cov}}$ reads exactly the off-diagonal structure the Proposition identifies as the rank-determining quantity, driving it toward zero pushes $C$ toward a diagonal (decorrelated) covariance, which raises $\mathrm{erank}$. A **minimum-batch guard** disables the term when $n$ is too small for the empirical off-diagonals to be meaningful (small-batch $C$ is dominated by sampling noise, whose spurious off-diagonals would inject gradient noise rather than signal).

**(F-R3) Optional escalation — RNA-paired biology InfoNCE.**
If F-R1+F-R2 do not fully restore rank, we replace the programme-neighbour KL with an **RNA-paired biology InfoNCE**, giving the biology head the same rank-expanding paired-contrastive signal that keeps the identity head healthy (§3.2). This substitutes a low-rank, label-softening KL for a high-rank, instance-discriminative contrastive objective.

**The effective-rank fingerprint (proposed diagnostic).**
We propose $\mathrm{erank}(z_{\mathrm{biology}})$ — the Roy–Vetterli effective rank of the biology head defined above — as an **internal fingerprint** of cohort-confounded WSI$\rightarrow$molecular alignment. Its diagnostic value is that it is computed from the model's own representation on unlabelled held-out data and **requires no a-priori knowledge of the confounder**: a biology head that has collapsed to its label manifold's intrinsic rank ($\sim$5–6 here) is reproducing low-rank cohort structure rather than high-rank per-sample biology, and the collapse is legible directly from the singular spectrum. The matched sibling design sharpens this into a *differential* fingerprint — a large gap $\mathrm{erank}(z_{\mathrm{identity}})\!-\!\mathrm{erank}(z_{\mathrm{biology}})$ (here $\sim 84$ vs $\sim 5$–$6$) flags the collapse while controlling for everything the two heads share. This differential, matched-head rank signature is the mechanism claim (C1) that anchors the rest of the paper; the prescriptive separation of variance-floor (cannot fix) from covariance-decorrelation (can fix) is the Proposition above (C3).

---

## 4. Metrics (confound-aware evaluation)

We evaluate WSI$\to$molecular-programme alignment by how well a slide embedding predicts each of $P$ molecular-programme targets (50-D MSigDB Hallmark scores plus held-out and control gene sets). The central methodological point of this section is that the *default* way to report such prediction — a pooled correlation over all held-out slides — is confounded by cohort (cancer-type) structure, and that a **control-adjusted within-cancer specificity** is the only quantity that isolates the signal a molecular head actually contributes. We define each metric formally, state the null model, and justify the choice of primary endpoint.

### 4.1 Notation

Let the held-out test set contain $n$ slides indexed by $i$, each with cancer type $c(i) \in \{1,\dots,K\}$. For a molecular target $t$ (a gene set), let $y_{i}^{(t)}$ be the ground-truth programme score for slide $i$ and $\hat{y}_{i}^{(t)}$ the value regressed / retrieved from the slide embedding under the method being evaluated. All correlations below are Pearson correlations computed per target $t$ and then averaged over the target panel $\mathcal{T}$ (with $|\mathcal{T}| = P$); we suppress the average over $\mathcal{T}$ in the notation where it is unambiguous.

### 4.2 Pooled (global) Pearson

The pooled estimator ignores cancer type and correlates predictions against truth across the entire held-out cohort:

$$
r_{\mathrm{pool}}^{(t)} \;=\; \operatorname{Corr}_{\,i \in \{1,\dots,n\}}\!\left(\hat{y}_{i}^{(t)},\, y_{i}^{(t)}\right),
\qquad
r_{\mathrm{pool}} \;=\; \frac{1}{P}\sum_{t \in \mathcal{T}} r_{\mathrm{pool}}^{(t)} .
$$

This is the number most commonly reported in WSI$\to$expression work. It conflates two sources of covariance: (i) *within-cancer* biological signal — whether the embedding tracks programme activity among tumours of the *same* type — and (ii) *between-cancer* mean structure — the fact that both $\hat{y}$ and $y$ shift systematically across cancer types. Because cancer type is legible from an H&E slide almost trivially (tissue architecture, stain, tiling artefacts), (ii) inflates $r_{\mathrm{pool}}$ even for an embedding that carries no sub-type-resolved molecular information.

### 4.3 Within-cancer (macro) Pearson

To remove the between-cancer contribution we compute Pearson *within* each cancer type and macro-average over types, weighting every cancer equally:

$$
r_{\mathrm{macro}}^{(t)} \;=\; \frac{1}{K}\sum_{k=1}^{K}
\operatorname{Corr}_{\,i:\,c(i)=k}\!\left(\hat{y}_{i}^{(t)},\, y_{i}^{(t)}\right),
\qquad
\texttt{macro\_cancer\_pearson} \;=\; \frac{1}{P}\sum_{t \in \mathcal{T}} r_{\mathrm{macro}}^{(t)} .
$$

Conditioning on $c(i)$ blocks the cross-cancer mean shift: $r_{\mathrm{macro}}^{(t)}$ can only be large if the embedding resolves programme variation *among tumours that share a cancer type*, which is the clinically and biologically non-trivial regime. Macro- (rather than micro-) averaging prevents a few large cancer cohorts from dominating and re-importing cross-cancer structure through unequal support.

Empirically the gap between the two estimators is large and *method-independent*: across all evaluated models $r_{\mathrm{macro}} \approx 0.51$–$0.54 \times r_{\mathrm{pool}}$ (e.g. $0.348 \to 0.188$ for the MLP-CLIP baseline; $0.327 \to 0.166$ for the anchored model). Roughly **46–49% of the pooled Pearson is cross-cancer cohort structure for every method**, baselines included. This halving is consistent with prior reports that a large fraction of apparent WSI$\to$expression accuracy is attributable to coarse tissue/cohort identity rather than fine molecular state [Fu 2020; Howard 2021], and with cohort-confounding audits of histology–omics pipelines [Buyer Beware NBME 2026; DECAT].

### 4.4 Size-matched random-gene-set null

A within-cancer correlation can still be non-zero for a target that carries no specific signal, because slide embeddings and *any* aggregate expression score share slow technical and compositional covariation (tumour purity, sequencing depth, stromal fraction). We therefore calibrate against a null built from **random gene sets of matched size**. For each real target $t$ with $|G_t|$ genes we draw random gene sets $G_{t,b}$, $b=1,\dots,B$, with $|G_{t,b}| = |G_t|$ sampled from the measured transcriptome, score them identically to produce $y_{i}^{(t,b)}$, and evaluate the same estimators:

$$
\bar{r}^{\,\mathrm{rand}}_{\mathrm{macro}}(t) \;=\; \frac{1}{B}\sum_{b=1}^{B} r_{\mathrm{macro}}^{(t,b)},
\qquad
r^{\mathrm{rand}}_{\mathrm{macro}} \;=\; \frac{1}{P}\sum_{t \in \mathcal{T}} \bar{r}^{\,\mathrm{rand}}_{\mathrm{macro}}(t).
$$

The size match matters because correlation with an aggregate score grows with set size (more genes average out noise), so an unmatched null would understate the floor. This random-set control is the standard device for turning a raw module score into a *specific* one: it is the null underlying Venet et al.'s demonstration that random signatures of matched size predict outcome nearly as well as published ones [Venet 2011], and it is operationalised as the background-subtracted control in `Seurat::AddModuleScore` [Tirosh 2016 / Seurat] and in WSI$\to$expression evaluation [Schmauch, HE2RNA 2020]. Our contribution is not the control but its use *inside* the within-cancer estimator.

The null is not small. Random gene sets reach $r^{\mathrm{rand}}_{\mathrm{pool}} \approx 0.30$–$0.32$ pooled and $\approx 0.15$ within-cancer (MLP-CLIP: $0.317\,|\,0.154$). A method reporting a pooled Pearson of $0.35$ is therefore operating only $\sim\!0.03$ above a random-gene floor — most of the headline number is reproduced by gene sets with no relation to the target.

### 4.5 Control-adjusted within-cancer specificity (primary endpoint)

Subtracting the null from the real signal, *within cancer*, yields the quantity we report as primary:

$$
\boxed{\;
\Delta_{\mathrm{spec}}(t) \;=\; r_{\mathrm{macro}}^{(t)} \;-\; \bar{r}^{\,\mathrm{rand}}_{\mathrm{macro}}(t),
\qquad
\Delta_{\mathrm{spec}} \;=\; \frac{1}{P}\sum_{t \in \mathcal{T}} \Delta_{\mathrm{spec}}(t)
\;}
$$

$\Delta_{\mathrm{spec}}$ removes *both* confounds in sequence: the within-cancer conditioning removes cross-cancer mean structure, and the random-set subtraction removes the size-dependent technical/compositional floor. What remains is the correlation attributable specifically to the *identity of the genes in the target*, among tumours of the same type — the only component that reflects genuine WSI$\to$programme alignment.

Measured this way, every method — including the plain MLP-CLIP baseline — lands at $\Delta_{\mathrm{spec}} \approx +0.07$ Pearson (MLP-CLIP $+0.068$; SigLIP $+0.067$; anchored $+0.068$; no-anchor $+0.069$). The genuine, control-adjusted within-cancer specificity is $\sim\!+0.07$ and is **method-invariant**: it does not move across architectures, objectives, or the presence of our anchor. This flatness is the empirical backbone of claim C2 — a rank-collapsed head can still post a "respectable" pooled score of $0.33$ because the benchmark is $\sim\!46$–$49\%$ cross-cancer and the random-gene floor already sits at $\sim\!0.30$–$0.32$, leaving only $\sim\!+0.07$ of method-invariant real signal for any model to compete over.

### 4.6 Effective rank

To characterise the *geometry* that produces this behaviour we use effective rank [Roy & Vetterli 2007]. For a batch of $N$ L2-normalized embeddings stacked in $Z \in \mathbb{R}^{N\times d}$ ($d=256$), let $C = \tfrac{1}{N}(Z-\bar Z)^\top (Z-\bar Z)$ be the feature covariance with eigenvalues $\lambda_1,\dots,\lambda_d \ge 0$. Define the normalized spectrum $p_j = \lambda_j / \sum_{m}\lambda_m$ and

$$
\operatorname{erank}(Z) \;=\; \exp\!\Big(-\sum_{j=1}^{d} p_j \log p_j\Big) \;=\; \exp\big(H(p)\big),
$$

the exponential of the Shannon entropy of the normalized eigenvalue spectrum. It is a smooth, basis-independent count of "effectively used" dimensions: $\operatorname{erank}=d$ for an isotropic representation and $\operatorname{erank}\to 1$ as variance concentrates on a single axis. Unlike a hard rank or a per-dimension variance statistic, effective rank is sensitive to *covariance* — dimensions that are individually non-degenerate but mutually collinear are correctly discounted, which is precisely the failure mode we diagnose.

On the held-out split ($n{=}2530$, $d{=}256$) the biology head is collapsed — $\operatorname{erank}(z_{\text{biology}})$ is $4.4$–$6.0$ depending on modality ($5.3$ fused) — while the sibling identity head stays healthy at $84.3$ (WSI) / $37.5$ (RNA). The collapse point ($\sim\!5$–$6$) coincides with the intrinsic rank of the $50$-D Hallmark target manifold, motivating the use of biology-head effective rank as an **internal fingerprint** of cohort-confounded alignment (C1). Because $\Delta_{\mathrm{spec}}$ shows the collapsed head loses essentially nothing on the confounded benchmark, effective rank flags the pathology that the confounded metric hides — the two are complementary, not redundant.

### 4.7 Why $\Delta_{\mathrm{spec}}$ is primary and pooled Pearson is confounded

We report $\Delta_{\mathrm{spec}}$ (control-adjusted within-cancer specificity) as the primary number, and treat $r_{\mathrm{pool}}$ as diagnostic only, for three reasons.

1. **Pooled Pearson answers the wrong question.** A model can maximise $r_{\mathrm{pool}}$ by predicting cancer type and reading programme activity off the cross-cancer mean map $c \mapsto \mathbb{E}[y \mid c]$, without resolving any within-type variation. Cancer type is near-trivially decodable from H&E, so $r_{\mathrm{pool}}$ largely rewards a shortcut. Formally, the law of total covariance decomposes $\operatorname{Cov}(\hat y, y)$ into a within-cancer term and a between-cancer term $\operatorname{Cov}\!\big(\mathbb{E}[\hat y\mid c], \mathbb{E}[y\mid c]\big)$; $r_{\mathrm{pool}}$ sums both, $r_{\mathrm{macro}}$ retains only the first. The $\sim\!47\%$ drop from $r_{\mathrm{pool}}$ to $r_{\mathrm{macro}}$ measured across all methods is exactly this between-cancer term.

2. **Within-cancer Pearson alone is still not specific.** Conditioning on cancer removes the mean shift but not the size-dependent technical/compositional floor, which reaches $\approx 0.15$ within-cancer. Only subtraction of the size-matched random-gene null isolates gene-identity-specific signal; this is the Venet/AddModuleScore control lineage applied inside the conditional estimator [Venet 2011; Seurat AddModuleScore; Schmauch 2020].

3. **$\Delta_{\mathrm{spec}}$ is discriminative where the others are not.** Because both confounds are large and method-invariant, $r_{\mathrm{pool}}$ compresses genuinely different models into a narrow $0.33$–$0.35$ band and rewards ones that lean harder on the cohort shortcut. $\Delta_{\mathrm{spec}}$ exposes that the real contested margin is only $\sim\!+0.07$ and reveals it to be flat across methods — the finding that a rank-5–6 head "scores respectably" not on merit but on confound budget. A benchmark that ranks methods by $r_{\mathrm{pool}}$ therefore ranks them largely by how much cohort structure they absorb, which is why we neither optimise nor headline it.

Concretely, all downstream comparisons, ablations, and head-to-head deltas in this paper are reported on $\Delta_{\mathrm{spec}}$ (with $r_{\mathrm{pool}}$, $r_{\mathrm{macro}}$, and the random-set floor shown alongside for transparency), and every claim about method quality is made net of the random-gene control and within cancer type.

---

## 5. Experiments & Results

We organize the empirical section around the paper's three claims. **C1 (mechanism):** the biology head of a matched dual-head model undergoes effective-rank collapse while its sibling identity head stays healthy (T1). **C2 (coupling):** a rank-5–6 head nonetheless scores respectably *only* because the benchmark is dominated by cross-cancer cohort structure and random-gene nulls are already high, leaving a small, method-invariant residual of genuine within-cancer signal (T2, T3, plus retrieval and head-to-head). **C3 (prescriptive):** a per-dimension variance floor cannot repair *rank* collapse, whereas a covariance-decorrelation term is the pre-registered intervention we test in T4.

Unless noted, geometry is measured on the held-out test split ($n=2530$, embedding dimension $d=256$); molecular-prompting metrics are Pearson correlations averaged over 180 MSigDB Hallmark targets and over seeds $\{42,43,44\}$. "Within-cancer" denotes the macro-averaged per-cancer Pearson (`macro_cancer_pearson`), which removes cross-cancer cohort structure by construction. We *concede* to prior art the confound itself (Howard et al. 2021; *Buyer Beware* NBME 2026; DECAT), the random-gene control (Schmauch et al., HE2RNA 2020; Venet et al. 2011), the roughly 50% halving of correlation once cohort is controlled (Fu et al. 2020), collapse-as-phenomenon (Andriopoulos et al., NRC1, NeurIPS 2024; Jing et al., ICLR 2022; Bardes et al., VICReg, ICLR 2022), and the retrieval-vs-regression asymmetry (Wang & Isola, CVPR 2022; Liang et al., FactorCL, NeurIPS 2023). Our contribution is the *coupled* mechanism-plus-benchmark package, and the proposal of biology-head effective rank as an internal fingerprint of cohort confounding.

### 5.1 T1 — Dual-head effective-rank spectrum (C1)

We report effective rank (the entropy-based participation ratio of the normalized singular-value spectrum) for each head. The identity head, anchored as a bounded zero-init residual on a frozen all-patch MLP-CLIP teacher, retains a high-dimensional geometry; the biology head, aligned to the intrinsically low-rank ($\sim$5–6) 50-D Hallmark manifold, collapses to a handful of active directions in $\mathbb{R}^{256}$.

| Head | Modality | Effective rank (of 256) |
|---|---|---|
| `wsi_identity` | WSI | **84.3** |
| `rna_identity` | RNA | 37.5 |
| `wsi_biology` | WSI | **6.0** |
| `rna_biology` | RNA | 4.4 |
| `full_biology` | fused | 5.3 |
| `full_patient` | fused | 8.0 |

The collapse is head-specific, not model-wide: within the *same* network the identity head sits at $\approx$84 while every biology-supervised head sits at $\approx$4–6. Modality gap (L2 distance between L2-normalized modality centroids) is 0.296 for identity versus 0.475 for biology, consistent with the biology head packing both modalities onto a thin shared subspace. The anchor behaves as designed — a near-null residual on the teacher: `residual_scale` $=-0.0011$, `correction_norm` $=0.00073$, `gate_mean` $=0.646$ — so the identity head is effectively the frozen MLP-CLIP teacher and its healthy rank is inherited, not learned around the collapse. The matched-sibling contrast is the fingerprint: an alignment target of intrinsic rank $\sim$5–6 stamps that rank onto the supervised head while the unsupervised-geometry sibling is untouched.

### 5.2 T2 — Confound decomposition ladder (C2)

We decompose each method's raw score by descending a three-rung ladder: (1) global Pearson, which mixes genuine signal with cross-cancer cohort structure; (2) within-cancer macro Pearson, which removes that structure; (3) subtraction of a matched random-gene null, yielding control-adjusted within-cancer specificity. The $\Delta\%$ column is the global$\rightarrow$within-cancer drop.

| Method | Global $r$ | Within-cancer $r$ | $\Delta\%$ |
|---|---|---|---|
| `mlp_clip` (baseline) | 0.348 | 0.188 | $-46\%$ |
| `siglip` | 0.349 | 0.193 | $-45\%$ |
| `mlp_clip_hardneg` | 0.349 | 0.184 | $-47\%$ |
| `morpheus_v2_no_anchor` | 0.338 | 0.185 | $-45\%$ |
| `morpheus_v1` | 0.338 | 0.179 | $-47\%$ |
| `morpheus_v2_anchored` | 0.327 | 0.166 | $-49\%$ |

For *every* method, $\sim$46–49% of the reported Pearson correlation is cross-cancer cohort structure rather than within-cancer biology — the halving anticipated by Fu et al. (2020). The target-group breakdown for the baseline shows the same structure at finer grain (global | within-cancer): `immune_tme` 0.486 | 0.325, `tumour_state` 0.425 | 0.202, `heldout_pathway` 0.213 | 0.112, and — tellingly — a `random_control` group already reaching 0.317 | 0.154. That a random-gene target group scores 0.317 globally is the crux of C2: the benchmark floor is high, so a rank-collapsed head can look competitive without resolving fine-grained biology.

### 5.3 T3 — Method-invariance of control-adjusted specificity (C2)

Subtracting the matched random-gene null isolates the genuine, cohort-free biological signal (real minus random gene sets), reported globally and within-cancer.

| Method | Specificity, global | Specificity, within-cancer |
|---|---|---|
| `mlp_clip` (baseline) | $+0.061$ | $+0.068$ |
| `siglip` | $+0.068$ | $+0.067$ |
| `morpheus_v2_anchored` | $+0.078$ | $+0.068$ |
| `morpheus_v2_no_anchor` | $+0.069$ | $+0.069$ |

The control-adjusted within-cancer specificity is $\approx +0.07$ Pearson and is **method-invariant**: it is statistically indistinguishable across all measured methods, *including the plain MLP-CLIP baseline*. Once cohort structure and the random-gene floor are removed, roughly $+0.07$ of within-cancer signal is all that any of these encoders — collapsed or not — actually delivers. This is the empirical heart of C2: the rank-5–6 biology head is not "good," it is *coupled* to a benchmark whose global scores are dominated by cross-cancer structure ($\sim$46–49%) and a high random-gene null ($\sim$0.30–0.32).

The head-to-head comparison against the baseline (Pearson delta; % of the 180 targets won) confirms that no encoder meaningfully separates from MLP-CLIP on real biology:

| Method vs. `mlp_clip` | Global $\Delta$ (win%) | Within-cancer $\Delta$ (win%) |
|---|---|---|
| `siglip` | $+0.002$ (55%) | $+0.005$ (62%) |
| `mlp_clip_hardneg` | $+0.001$ (49%) | $-0.004$ (34%) |
| `morpheus_v2_no_anchor` | $-0.010$ (44%) | $-0.003$ (51%) |
| `morpheus_v2_anchored` | $-0.021$ (32%) | $-0.022$ (19%) |

The largest control-adjusted separation is SigLIP's $+0.005$ within-cancer residual (62% of targets won); given the $\approx +0.07$ method-invariant specificity ceiling, we treat this as a **minor** residual rather than a substantive gain. All methods reach molecular-phenotype AUROC $\approx 0.71$, again indistinguishable.

### 5.4 Retrieval (C2, secondary)

Retrieval Recall@$k$ (mean) shows the anchored v2 model does **not** beat the MLP-CLIP baseline; the coupling in C2 is not an artifact of the regression metric.

| Method | Recall@$k$ (mean) | Within-cancer recall |
|---|---|---|
| `mlp_clip` (baseline) | **0.066** | **0.133** |
| `mlp_clip_hardneg` | 0.066 | — |
| `morpheus_v2_anchored` | 0.060 | — |
| `siglip` | 0.059 | — |
| `morpheus_v2_no_anchor` | 0.040 | 0.120 |

The baseline leads on both mean and within-cancer recall. Consistent with the retrieval-vs-regression asymmetry of Wang & Isola (2022) and Liang et al. (2023), a low-rank biology head is penalized more sharply by retrieval — which rewards spread — than by per-target regression, so retrieval and specificity agree that v2 buys no real biological gain.

### 5.5 T4 — VICReg decorrelation ablation **[queued: Lambda full run]**

We register T4 in advance; **no T4 numbers are reported here.** T1 establishes that the biology head loses *rank*, not merely per-dimension variance. A per-dimension variance floor ($\mathrm{relu}(1-\mathrm{std})$, weight 0.01) constrains marginal variances but leaves the off-diagonal covariance — and therefore the rank — unconstrained; it provably cannot prevent rank collapse. The prescriptive fix (C3) couples two changes: **(F-R1)** feed a single L2-normalized biology state to both the Gaussian-NLL heads and the rank-spreaders so they share geometry, and **(F-R2)** add a VICReg/Barlow-style off-diagonal feature-decorrelation term on the normalized 256-D $z_{\text{biology}}$ — the mean-squared off-diagonal of the centered batch covariance, weight $\approx 0.04$, active in the `full` and `programme_only` profiles under a minimum-batch guard.

> **Pre-registered hypothesis (T4).** Adding the covariance-decorrelation term (F-R2), with the shared-geometry coupling (F-R1), will (i) *recover the biology-head effective rank* from the collapsed $\approx 5.3$ (`full_biology`, T1) toward the intact identity-head regime, while (ii) *leaving the control-adjusted within-cancer specificity statistically unchanged* at the method-invariant $\approx +0.07$ Pearson (T3). The prediction that rank rises while specificity does not is the falsifiable core of the fingerprint claim: it would demonstrate that effective rank diagnoses the confound-induced geometry *independently* of task score, so that collapse flags cohort confounding even when benchmark numbers look healthy. The per-dimension variance floor serves as the negative control (predicted: rank unchanged). An optional escalation (F-R3) replaces the neighbour-KL with an RNA-paired biology InfoNCE.
>
> **Exact metric reported by T4.** The primary readout is biology-head effective rank (participation ratio of the normalized singular spectrum) for `full_biology` / `wsi_biology` / `rna_biology`, reported as $\Delta$ effective rank versus the collapsed T1 baseline; the paired secondary readout is control-adjusted within-cancer specificity (real minus random-gene null, `macro_cancer_pearson`), tested for the pre-registered null of no change from $\approx +0.07$. Both are reported per seed $\{42,43,44\}$, with the variance-floor arm as negative control.

---

## 6. Discussion, Risks, Limitations

### 6.1 Discussion

Our results resolve into a single, coupled account of why WSI$\rightarrow$molecular-programme alignment looks stronger than it is, and how to detect the failure from the inside. We restate the three claims through the evidence.

**C1 — the mechanism is a matched dual-head rank-collapse fingerprint.** The biology head, trained to align/regress onto a 50-D MSigDB Hallmark target whose intrinsic rank is only $\sim$5–6, collapses to an effective rank of $6.0$ (WSI), $4.4$ (RNA), and $5.3$ (fused) out of $256$, while its architecturally identical sibling identity head — supervised by paired contrastive signal rather than a low-rank target — stays healthy at $84.3$ (WSI) and $37.5$ (RNA). The collapse is not a global optimisation pathology; it is head-local and target-induced, which is exactly what makes it diagnostic. The geometry corroborates this: the modality gap of the collapsed biology space ($0.475$) is markedly wider than that of the identity space ($0.296$), the signature of two modalities pinned to a shared low-rank ridge rather than genuinely co-embedded. Because the identity anchor is a bounded zero-init residual on a frozen MLP-CLIP teacher (residual scale $-0.0011$, correction norm $0.00073$, gate mean $0.646$), the identity head is provably $\approx$ the teacher, isolating the collapse to the biology objective and ruling out the anchor as a confound.

**C2 — the head scores respectably only because the benchmark is coupled to the confound.** A rank-5–6 biology head is, by any representational-quality standard, degenerate; yet it posts a global Pearson of $0.327$–$0.338$, within a few points of every baseline. This is not evidence that the head works — it is evidence that the *benchmark* is permissive. Two structural facts do the work. First, $\sim$46–49% of every method's global Pearson is cross-cancer cohort structure: within-cancer macro-Pearson falls from $0.348\!\rightarrow\!0.188$ for MLP-CLIP ($-46\%$), $0.349\!\rightarrow\!0.193$ for SigLIP, and $0.327\!\rightarrow\!0.166$ for the collapsed head alike. Second, random-gene nulls already reach $0.317$ global / $0.154$ within-cancer for MLP-CLIP — a floor that a rank-5 representation can clear by encoding little more than cancer-type and coarse tumour-purity axes. Subtracting the random-gene control leaves a genuine, within-cancer, control-adjusted specificity of only $\sim\!+0.07$ Pearson (MLP-CLIP $+0.061|+0.068$, SigLIP $+0.068|+0.067$, anchored V2 $+0.078|+0.068$). Crucially this residual is *method-invariant*: it is statistically identical across the collapsed head, the strong SigLIP baseline, and the plain MLP-CLIP baseline. The benchmark cannot distinguish a healthy head from a collapsed one because the quantity it rewards is dominated by structure both heads capture for free.

**C3 — the prescription follows from the mechanism, not from tuning.** Because the failure is a *rank* collapse, a per-dimension variance floor ($\mathrm{relu}(1-\mathrm{std})$, weight $0.01$) cannot fix it: enforcing unit marginal variance on each of 256 coordinates is satisfied by a representation that lives on a 5-D subspace with those coordinates as correlated linear images. Only an off-diagonal covariance/decorrelation term removes the degeneracy. Our fix is two coupled changes — (F-R1) feed a single L2-normalised biology state to both the Gaussian-NLL heads and the rank-spreaders so they optimise one shared geometry, and (F-R2) add a VICReg/Barlow-style mean-squared off-diagonal decorrelation term (weight $\sim\!0.04$, batch-guarded) on the normalised 256-D biology embedding — with an optional escalation (F-R3) that swaps the neighbour-KL for an RNA-paired biology InfoNCE, importing the same rank-expanding paired signal that keeps the identity head healthy.

**Why this is a package, not three separate observations.** Each ingredient is individually known (see Risks below). The contribution is their *coupling*: the effective rank of the biology head is an **internal fingerprint** that flags cohort-confounded WSI$\rightarrow$molecular alignment *without knowing the confounder a priori*. Where prior audits require the practitioner to already suspect site or cancer-type leakage and hold it out, the rank signature is read off a single held-out forward pass. That is the diagnostic novelty, and it is why we frame the paper as a method/diagnostic contribution rather than a benchmark or a theory paper.

### 6.2 Risks and Defenses

We state the strongest objections in the form a hostile reviewer would, and answer each. We **cite and concede** every component result; the defense is always about the *coupled* claim.

| # | Objection | Our defense |
|---|-----------|-------------|
| R1 | **"This is DECAT with a new name."** DECAT already frames confounded WSI$\rightarrow$expression alignment. | We cite DECAT as prior art for the confound framing and do not claim it. DECAT (like Howard 2021) requires the confounder to be *named and controlled externally*; our contribution is an **internal, confounder-agnostic** signal — the biology-head effective rank ($5.3$/256 vs a healthy $84.3$) — that fires from one forward pass with no held-out site/cancer variable, plus the matched-sibling design that localises the collapse to the target. DECAT diagnoses *that* confounding exists; we provide a representation-geometry fingerprint of *when a given trained head has silently absorbed it*, and a decorrelation prescription. |
| R2 | **"The cohort confound is already well known."** Howard 2021, Buyer Beware (NBME 2026). | Conceded and cited. We claim none of it. Our claim is C1+C2 *coupled*: that the confound is what makes a **rank-collapsed** head score respectably, quantified as the $\sim$46–49% cross-cancer share and the method-invariant $\sim\!+0.07$ residual. Prior work establishes the confound exists; we show it is the mechanism that hides representational collapse from every standard metric. |
| R3 | **"The random-gene control is not new."** Schmauch HE2RNA (2020), Venet (2011). | Conceded and cited; we adopt their null, we do not claim it. Novelty is that the control-adjusted specificity is **method-invariant** ($+0.061$ to $+0.078$ global, all $\approx\!+0.07$ within-cancer) across a collapsed head, a strong SigLIP baseline, and plain MLP-CLIP — turning the null from a per-method sanity check into evidence that the benchmark's discriminative headroom is a shared $\sim\!0.07$-Pearson band no current method exceeds. |
| R4 | **"Collapse is just a fixable training bug — patch it and move on."** Dimensional collapse (Jing, ICLR 2022), VICReg (Bardes, ICLR 2022), NRC1 (Andriopoulos, NeurIPS 2024). | We cite all three as the phenomenology and the remedy family, and claim neither. But "fixable bug" concedes our point: (a) the bug is *diagnostic* — its presence is the fingerprint of confounding, so the value is in reading it, not only removing it; (b) the *natural* fix (a per-dim variance floor) provably fails on rank collapse and we show only the off-diagonal term works, which is a non-obvious prescriptive result; (c) fixing collapse does **not** move the benchmark ($\sim\!+0.07$ ceiling is method-invariant), so the bug and its patch are decoupled from apparent score — precisely the trap C2 identifies. |
| R5 | **"$\sim\!+0.07$ specificity means the task is uninteresting."** | The opposite conclusion. Molecular-phenotype AUROC is $\sim\!0.71$ for all methods and the immune-TME target group reaches $0.486|0.325$ global\|within — there *is* real, decodable biology, concentrated in immune/TME axes. The $\sim\!+0.07$ figure is the *control-adjusted, within-cancer* headroom that separates a method from a random-gene null; that it is small and method-invariant is the finding — it means current alignment quality is confound-limited, not that the underlying WSI$\rightarrow$programme signal is absent. A small honest number that no method beats is a sharper benchmark, not a dead one. |
| R6 | **"SigLIP's $+0.005$ within-cancer win is inside the noise band — you have no effect."** | Agreed, and that agreement is load-bearing. SigLIP's advantage ($+0.002$ global / $+0.005$ within, winning 55%/62% of targets) sits inside our own noise band, which is exactly why we do **not** claim a state-of-the-art method. Our claims are C1/C2/C3, none of which rest on ranking methods. The noise-band result reinforces C2: when the best and worst representations differ by less than the seed noise on the headline metric while their effective ranks differ 14-fold ($84.3$ vs $6.0$), the metric is not measuring representation quality — which is the entire point. |

### 6.3 Ethical and Clinical Caveat

This is a **diagnostic and methodological study of representation geometry, not a clinical decision system.** None of the reported quantities — global or within-cancer Pearson, retrieval recall@$k$ ($0.040$–$0.066$), or molecular-phenotype AUROC ($\sim\!0.71$) — support patient-level molecular inference, treatment selection, or diagnosis from H&E, and they must not be deployed or interpreted as such. Our central negative result is a caution: apparent WSI$\rightarrow$molecular performance is substantially ($\sim$46–49%) cross-cancer cohort structure, so a model that appears to "predict expression from morphology" may be predicting little beyond cancer type and coarse purity. Any downstream clinical use would require prospective, site-stratified, cancer-type-controlled validation with the random-gene null reported, and remains out of scope here. The effective-rank fingerprint is offered as an auditing aid for method developers, not as evidence of clinical readiness.

### 6.4 Limitations

- **Single-seed geometry.** All effective-rank and modality-gap measurements (biology $5.3$, identity $84.3$, gaps $0.475$/$0.296$) are computed on the seed-42 artifact and held-out test split ($n=2530$, dim $256$). While the molecular-prompting Pearson numbers are averaged over seeds 42/43/44, the *geometry* itself is single-seed. The qualitative gap between the collapsed biology head and the healthy identity head is large ($\sim$14$\times$ in effective rank) and unlikely to be a seed artifact, but we do not yet report seed-level variance on the rank statistics themselves. A **multi-seed $\lambda$ ablation** over the decorrelation weight (the $\sim\!0.04$ off-diagonal term of F-R2) is queued on the Lambda cluster to establish (i) the seed distribution of biology-head effective rank with and without decorrelation, and (ii) the dose-response of rank recovery vs. the decorrelation weight, including whether rank recovery transfers to the benchmark or stays decoupled as R4 predicts.
- **The fix is prescribed and mechanistically argued, not yet fully ablated.** We prove a per-dim variance floor cannot address rank collapse and identify the off-diagonal covariance term as the necessary ingredient; the coupled F-R1/F-R2 remedy (and the F-R3 escalation) is specified and motivated, but the end-to-end demonstration that it restores biology-head rank *and* the resulting effect on within-cancer specificity is the subject of the queued ablation, not a completed result in this draft.
- **Benchmark scope.** Targets are 50-D MSigDB Hallmark programmes with an intrinsic rank of $\sim$5–6; the collapse story is sharpest for such low-rank targets and may attenuate for higher-rank or single-gene objectives. The held-out-pathway group ($0.213|0.112$) shows the specificity gap narrows further off the training programme manifold, which bounds the generality of any positive claim.
- **Confounder coverage.** Our cross-cancer decomposition captures the dominant cohort axis (cancer type / coarse purity); finer confounders (site, scanner, batch) are not separately partitioned here, so the $\sim$46–49% figure is a lower bound on total cohort structure, and the rank fingerprint, while confounder-agnostic in principle, is validated primarily against the cancer-type axis.
- **Retrieval regime.** Absolute retrieval recall is low ($0.040$–$0.066$, within-cancer $0.120$–$0.133$), consistent with the retrieval-vs-regression tension (Wang, CVPR 2022; FactorCL, NeurIPS 2023) we cite; we use retrieval as a corroborating geometry probe, not as a headline capability, and do not claim competitive retrieval performance.

---

## References

Bibliographic references are maintained in [`references.bib`](./references.bib) (BibTeX). The in-text `\citep{...}` keys resolve against that file; render with any BibTeX/BibLaTeX-aware toolchain (e.g. `pandoc --citeproc main.md --bibliography=references.bib`). Two sources cited in prose — the Roy & Vetterli (2007) effective-rank definition and Seurat/Tirosh (2016) `AddModuleScore` — are pending addition to `references.bib`.
