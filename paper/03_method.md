## Method (architecture + the fix)

We first describe the MORPHEUS V2 architecture that produces the two embedding
heads at the centre of our study (§M.1). We then diagnose the failure mode —
a rank collapse localised to the *biology* head while its sibling *identity*
head stays healthy (§M.2) — and prove that the natural remedy, a per-dimension
variance floor, cannot address it. Finally we give the fix: a *normalize-once*
coupling (F-R1) combined with an off-diagonal covariance-decorrelation term
(F-R2), and we define the **effective-rank fingerprint** that we propose as an
internal, confounder-agnostic diagnostic for cohort-confounded
WSI$\rightarrow$molecular-programme alignment (§M.3).

Throughout, a *WSI* is a whole-slide image, $z\in\mathbb{R}^{d}$ with $d=256$
is an L2-normalised embedding, and a *batch* is a matrix
$Z=[z_1,\dots,z_n]^{\top}\in\mathbb{R}^{n\times d}$ of $n$ such embeddings.

---

### M.1 The MORPHEUS V2 architecture

MORPHEUS V2 is a dual-encoder, dual-head model that maps a whole-slide image
and its paired bulk RNA-seq profile into a shared 256-D space, split into an
**identity** subspace (patient/slide provenance) and a **biology** subspace
(molecular programme state).

**Vision branch — hierarchical Query-Former.**
Uncapped H-Optimus patch features (all patches per slide, no fixed cap) are
pooled by a hierarchical soft-slot attention module — a *Query-Former* — that
performs learned soft-slot pooling at three granularities: 64 **local** slots,
32 **slide** slots, and 128 **patient** slots. Each level attends over the
tokens produced by the level below, so a patient representation is a soft
aggregation of slide summaries, which are in turn soft aggregations of local
tissue neighbourhoods. Soft-slot pooling (attention weights over an
over-complete slot dictionary, rather than hard clustering) keeps the pooling
differentiable and permutation-invariant while accommodating a variable, and
typically large, number of input patches.

**Molecular branch.**
Paired bulk RNA-seq is encoded by BulkFormer, yielding an RNA token summary
that is projected into the same 256-D geometry as the vision heads.

**Two heads.**
From the fused representation the model emits two L2-normalised vectors:

- $z_{\mathrm{identity}}\in\mathbb{S}^{255}$ — a **bounded, zero-initialised
  residual** on a *frozen all-patch MLP-CLIP teacher*. Writing $t$ for the
  (unit-norm) teacher embedding, a learned gate $g\in[0,1]$ and a bounded
  correction $\delta$,
  $$
  z_{\mathrm{identity}}
  \;=\;
  \mathrm{normalize}\!\big(\,t + g\cdot s\cdot \delta\,\big),
  \qquad s = \text{residual scale},
  $$
  with $\delta$ produced by the anchor head and $s$ a small learned scalar
  initialised at $0$. Zero-init means the head *begins as the teacher* and can
  only depart from it to the extent the residual earns; in our trained model
  the residual is effectively inert ($s=-0.0011$, mean correction norm
  $0.00073$, gate mean $0.646$), so $z_{\mathrm{identity}}$ is, to numerical
  precision, the frozen MLP-CLIP teacher. This *anchoring* is what makes the
  identity head a clean, healthy control against which the biology head's
  collapse is measured.

- $z_{\mathrm{biology}}\in\mathbb{S}^{255}$ — supervised toward a 50-D MSigDB
  **Hallmark** programme target vector $y\in\mathbb{R}^{50}$.

**Biology-head supervision.**
The biology head is trained with three objectives against the Hallmark
targets:

1. a **heteroscedastic Gaussian negative log-likelihood** (Gaussian-NLL): a
   read-out head predicts a mean $\mu(z_{\mathrm{biology}})$ and a
   per-programme log-variance $\log\sigma^2(z_{\mathrm{biology}})$, and is
   trained by
   $$
   \mathcal{L}_{\mathrm{NLL}}
   = \frac{1}{2}\sum_{k=1}^{50}
   \left[
   \frac{\big(y_k-\mu_k\big)^2}{\sigma_k^2}
   + \log \sigma_k^2
   \right];
   $$
2. a **programme-neighbour KL** that softens supervision toward
   programme-space neighbours (encouraging embeddings of biologically adjacent
   programmes to share a neighbourhood); and
3. a **supervised contrastive** term over Hallmark labels.

The Hallmark target manifold is nominally 50-D but has an **intrinsic rank of
only $\sim$5–6**: the 50 Hallmark programme scores are strongly co-varying, so
the label geometry the biology head is asked to reproduce is extremely
low-rank. This fact is the seed of the collapse.

---

### M.2 Diagnosis: a localised rank collapse

We measure the geometry of each head on the held-out test split
($n=2530$, $d=256$) via **effective rank**. Given a batch $Z$ with centred
covariance $C=\tfrac1n Z_c^{\top}Z_c$ (equivalently the singular values
$\sigma_1\ge\dots\ge\sigma_d\ge 0$ of $Z_c/\sqrt{n}$), define the
$\ell_1$-normalised spectrum $p_i=\sigma_i / \sum_j \sigma_j$ and the
Roy–Vetterli effective rank

$$
\boxed{\;
\mathrm{erank}(Z)
\;=\;
\exp\!\Big( -\sum_{i=1}^{d} p_i \log p_i \Big)
\;=\;
\exp\!\big(H(p_1,\dots,p_d)\big)
\;}
$$

the exponential of the Shannon entropy of the (normalised) singular-value
distribution. $\mathrm{erank}\in[1,d]$: it equals $d$ for an isotropic
(white) representation and $1$ when all variance lies along a single
direction. It is a smooth, basis-independent measure of how many directions
the representation actually populates.

On the held-out split the two heads diverge sharply:

| head              | effective rank |
|-------------------|---------------:|
| `wsi_identity`    | 84.3           |
| `rna_identity`    | 37.5           |
| `wsi_biology`     | 6.0            |
| `rna_biology`     | 4.4            |
| `full_biology`    | 5.3            |
| `full_patient`    | 8.0            |

The identity head, anchored to the frozen teacher and trained with an
RNA-paired contrastive signal, occupies $\sim$84 of 256 directions. The
biology head — regressed and aligned to the intrinsic-rank-$\sim$5–6 Hallmark
manifold — has **collapsed to $\sim$5–6 effective directions**, i.e. it has
inherited the low-rank geometry of its supervision targets almost exactly.
This is a *dimensional/rank collapse* in the sense of Jing et al. (2022) and
NRC1 (Andriopoulos et al., 2024), but *localised to one head of a matched
sibling pair* — the healthy identity head serves as a within-model control
that rules out data-, optimiser-, or capacity-level explanations. The modality
(WSI$\leftrightarrow$RNA) gap between L2-normed centroids is correspondingly
larger for the collapsed head (biology $0.475$) than for the healthy one
(identity $0.296$).

**Why a per-dimension variance floor cannot fix it.**
The standard VICReg *variance* term applies a hinge floor to each coordinate's
standard deviation,
$\mathcal{L}_{\mathrm{var}}
= \tfrac1d\sum_{j=1}^{d}\mathrm{relu}\!\big(\gamma-\mathrm{std}(Z_{:,j})\big)$
(here $\gamma=1$, weight $0.01$). We state precisely why this is powerless
against *rank* collapse.

> **Proposition (variance floor $\not\Rightarrow$ rank).**
> Let $Z_c\in\mathbb{R}^{n\times d}$ be a centred batch and let
> $\mathrm{std}(Z_{:,j})=\sqrt{C_{jj}}$ denote the $j$-th coordinate standard
> deviation, $C=\tfrac1n Z_c^\top Z_c$. The per-dimension floor
> $\mathcal{L}_{\mathrm{var}}$ is a function of the **diagonal** of $C$ only.
> Its global minimum ($\mathcal{L}_{\mathrm{var}}=0$) is attained by *any* $C$
> with $C_{jj}\ge\gamma^2$ for all $j$ — including rank-1 configurations. In
> particular, take $Z_c = a\,u^\top$ for a fixed unit direction
> $u\in\mathbb{R}^d$ with $|u_j|\ge\gamma/\|a\|_{\mathrm{rms}}$ and an
> activation vector $a\in\mathbb{R}^n$: every coordinate meets its variance
> floor, yet $\mathrm{rank}(C)=1$ and $\mathrm{erank}(Z)\to 1$. Hence
> $\mathcal{L}_{\mathrm{var}}=0$ is consistent with maximal rank collapse.

*Proof sketch.* $\mathcal{L}_{\mathrm{var}}$ depends on $C$ solely through
$\{C_{jj}\}$. The rank of $C$ is governed by its **off-diagonal** structure:
a matrix with large, satisfied diagonal entries can still be exactly rank-1 if
its rows are mutually proportional (perfect off-diagonal correlation). The
explicit $Z_c=a\,u^\top$ construction realises this: it drives every marginal
variance above the floor while confining all mass to a single eigen-direction.
Marginal variances are invariant to the correlations that determine rank; only
a term that reads the **off-diagonal** covariance can constrain rank. $\square$

The proposition motivates the exact form of the fix: decorrelation, not
per-coordinate scaling.

---

### M.3 The fix: normalize-once (F-R1) + covariance decorrelation (F-R2)

The fix is two *coupled* changes to the biology-head objective, plus an
optional escalation.

**(F-R1) Normalize-once — shared geometry.**
In the collapsed model the Gaussian-NLL read-out heads and the
rank-spreading regularisers consumed *different* views of the biology state
(pre- vs. post-normalisation), so a regulariser could inflate the geometry of
a state the likelihood never saw. F-R1 feeds a **single L2-normalised
biology state** $z_{\mathrm{biology}}=\mathrm{normalize}(h)$ to *both* the
Gaussian-NLL heads *and* the decorrelation term, so likelihood and geometry
act on the same object. This coupling is what makes F-R2 bite: decorrelating
a state that the supervision also optimises forces the two to co-adapt rather
than fight.

**(F-R2) VICReg / Barlow-style off-diagonal covariance decorrelation.**
On the normalised 256-D biology batch $Z\in\mathbb{R}^{n\times d}$ we add a
term that penalises the off-diagonal entries of the centred batch covariance.
Let $\bar z=\tfrac1n\sum_i z_i$, $Z_c = Z-\mathbf 1\bar z^\top$, and

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

added to the objective with weight $\lambda_{\mathrm{cov}}\approx 0.04$, active
in the `full` and `programme_only` objective profiles. Because $\mathcal L_
{\mathrm{cov}}$ reads exactly the off-diagonal structure the Proposition
identifies as the rank-determining quantity, driving it toward zero pushes
$C$ toward a diagonal (decorrelated) covariance, which raises $\mathrm{erank}$.
A **minimum-batch guard** disables the term when $n$ is too small for the
empirical off-diagonals to be meaningful (small-batch $C$ is dominated by
sampling noise, whose spurious off-diagonals would inject gradient noise
rather than signal).

**(F-R3) Optional escalation — RNA-paired biology InfoNCE.**
If F-R1+F-R2 do not fully restore rank, we replace the programme-neighbour KL
with an **RNA-paired biology InfoNCE**, giving the biology head the same
rank-expanding paired-contrastive signal that keeps the identity head healthy
(§M.2). This substitutes a low-rank, label-softening KL for a high-rank,
instance-discriminative contrastive objective.

**The effective-rank fingerprint (proposed diagnostic).**
We propose $\mathrm{erank}(z_{\mathrm{biology}})$ — the Roy–Vetterli effective
rank of the biology head defined above — as an **internal fingerprint** of
cohort-confounded WSI$\rightarrow$molecular alignment. Its diagnostic value is
that it is computed from the model's own representation on unlabelled held-out
data and **requires no a-priori knowledge of the confounder**: a biology head
that has collapsed to its label manifold's intrinsic rank ($\sim$5–6 here) is
reproducing low-rank cohort structure rather than high-rank per-sample
biology, and the collapse is legible directly from the singular spectrum. The
matched sibling design sharpens this into a *differential* fingerprint — a
large gap $\mathrm{erank}(z_{\mathrm{identity}})\!-\!\mathrm{erank}
(z_{\mathrm{biology}})$ (here $\sim 84$ vs $\sim 5$–$6$) flags the collapse
while controlling for everything the two heads share. This differential,
matched-head rank signature is the mechanism claim (C1) that anchors the rest
of the paper; the prescriptive separation of variance-floor (cannot fix) from
covariance-decorrelation (can fix) is the Proposition above (C3).
