> **Note:** `main.md` is the authoritative assembled paper. This section draft predates the P4 metric correction (effective rank restated as Roy–Vetterli / singular-value: biology ~33–47, identity ~142–191); any `~5–6`/`84` figures here are superseded by main.md.

## Metrics (confound-aware evaluation)

We evaluate WSI$\to$molecular-programme alignment by how well a slide embedding predicts each of $P$ molecular-programme targets (50-D MSigDB Hallmark scores plus held-out and control gene sets). The central methodological point of this section is that the *default* way to report such prediction — a pooled correlation over all held-out slides — is confounded by cohort (cancer-type) structure, and that a **control-adjusted within-cancer specificity** is the only quantity that isolates the signal a molecular head actually contributes. We define each metric formally, state the null model, and justify the choice of primary endpoint.

### Notation

Let the held-out test set contain $n$ slides indexed by $i$, each with cancer type $c(i) \in \{1,\dots,K\}$. For a molecular target $t$ (a gene set), let $y_{i}^{(t)}$ be the ground-truth programme score for slide $i$ and $\hat{y}_{i}^{(t)}$ the value regressed / retrieved from the slide embedding under the method being evaluated. All correlations below are Pearson correlations computed per target $t$ and then averaged over the target panel $\mathcal{T}$ (with $|\mathcal{T}| = P$); we suppress the average over $\mathcal{T}$ in the notation where it is unambiguous.

### Pooled (global) Pearson

The pooled estimator ignores cancer type and correlates predictions against truth across the entire held-out cohort:

$$
r_{\mathrm{pool}}^{(t)} \;=\; \operatorname{Corr}_{\,i \in \{1,\dots,n\}}\!\left(\hat{y}_{i}^{(t)},\, y_{i}^{(t)}\right),
\qquad
r_{\mathrm{pool}} \;=\; \frac{1}{P}\sum_{t \in \mathcal{T}} r_{\mathrm{pool}}^{(t)} .
$$

This is the number most commonly reported in WSI$\to$expression work. It conflates two sources of covariance: (i) *within-cancer* biological signal — whether the embedding tracks programme activity among tumours of the *same* type — and (ii) *between-cancer* mean structure — the fact that both $\hat{y}$ and $y$ shift systematically across cancer types. Because cancer type is legible from an H&E slide almost trivially (tissue architecture, stain, tiling artefacts), (ii) inflates $r_{\mathrm{pool}}$ even for an embedding that carries no sub-type-resolved molecular information.

### Within-cancer (macro) Pearson

To remove the between-cancer contribution we compute Pearson *within* each cancer type and macro-average over types, weighting every cancer equally:

$$
r_{\mathrm{macro}}^{(t)} \;=\; \frac{1}{K}\sum_{k=1}^{K}
\operatorname{Corr}_{\,i:\,c(i)=k}\!\left(\hat{y}_{i}^{(t)},\, y_{i}^{(t)}\right),
\qquad
\texttt{macro\_cancer\_pearson} \;=\; \frac{1}{P}\sum_{t \in \mathcal{T}} r_{\mathrm{macro}}^{(t)} .
$$

Conditioning on $c(i)$ blocks the cross-cancer mean shift: $r_{\mathrm{macro}}^{(t)}$ can only be large if the embedding resolves programme variation *among tumours that share a cancer type*, which is the clinically and biologically non-trivial regime. Macro- (rather than micro-) averaging prevents a few large cancer cohorts from dominating and re-importing cross-cancer structure through unequal support.

Empirically the gap between the two estimators is large and *method-independent*: across all evaluated models $r_{\mathrm{macro}} \approx 0.51$–$0.54 \times r_{\mathrm{pool}}$ (e.g. $0.348 \to 0.188$ for the MLP-CLIP baseline; $0.327 \to 0.166$ for the anchored model). Roughly **46–49% of the pooled Pearson is cross-cancer cohort structure for every method**, baselines included. This halving is consistent with prior reports that a large fraction of apparent WSI$\to$expression accuracy is attributable to coarse tissue/cohort identity rather than fine molecular state [Fu 2020; Howard 2021], and with cohort-confounding audits of histology–omics pipelines [Buyer Beware (Dawood et al. 2024); DECAT].

### Size-matched random-gene-set null

A within-cancer correlation can still be non-zero for a target that carries no specific signal, because slide embeddings and *any* aggregate expression score share slow technical and compositional covariation (tumour purity, sequencing depth, stromal fraction). We therefore calibrate against a null built from **random gene sets of matched size**. For each real target $t$ with $|G_t|$ genes we draw random gene sets $G_{t,b}$, $b=1,\dots,B$, with $|G_{t,b}| = |G_t|$ sampled from the measured transcriptome, score them identically to produce $y_{i}^{(t,b)}$, and evaluate the same estimators:

$$
\bar{r}^{\,\mathrm{rand}}_{\mathrm{macro}}(t) \;=\; \frac{1}{B}\sum_{b=1}^{B} r_{\mathrm{macro}}^{(t,b)},
\qquad
r^{\mathrm{rand}}_{\mathrm{macro}} \;=\; \frac{1}{P}\sum_{t \in \mathcal{T}} \bar{r}^{\,\mathrm{rand}}_{\mathrm{macro}}(t).
$$

The size match matters because correlation with an aggregate score grows with set size (more genes average out noise), so an unmatched null would understate the floor. This random-set control is the standard device for turning a raw module score into a *specific* one: it is the null underlying Venet et al.'s demonstration that random signatures of matched size predict outcome nearly as well as published ones [Venet 2011], and it is operationalised as the background-subtracted control in `Seurat::AddModuleScore` [Tirosh 2016 / Seurat] and in WSI$\to$expression evaluation [Schmauch, HE2RNA 2020]. Our contribution is not the control but its use *inside* the within-cancer estimator.

The null is not small. Random gene sets reach $r^{\mathrm{rand}}_{\mathrm{pool}} \approx 0.30$–$0.32$ pooled and $\approx 0.15$ within-cancer (MLP-CLIP: $0.317\,|\,0.154$). A method reporting a pooled Pearson of $0.35$ is therefore operating only $\sim\!0.03$ above a random-gene floor — most of the headline number is reproduced by gene sets with no relation to the target.

### Control-adjusted within-cancer specificity (primary endpoint)

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

### Effective rank

To characterise the *geometry* that produces this behaviour we use effective rank [Roy & Vetterli 2007]. For a batch of $N$ L2-normalized embeddings stacked in $Z \in \mathbb{R}^{N\times d}$ ($d=256$), let $C = \tfrac{1}{N}(Z-\bar Z)^\top (Z-\bar Z)$ be the feature covariance with eigenvalues $\lambda_1,\dots,\lambda_d \ge 0$. Define the normalized spectrum $p_j = \lambda_j / \sum_{m}\lambda_m$ and

$$
\operatorname{erank}(Z) \;=\; \exp\!\Big(-\sum_{j=1}^{d} p_j \log p_j\Big) \;=\; \exp\big(H(p)\big),
$$

the exponential of the Shannon entropy of the normalized eigenvalue spectrum. It is a smooth, basis-independent count of "effectively used" dimensions: $\operatorname{erank}=d$ for an isotropic representation and $\operatorname{erank}\to 1$ as variance concentrates on a single axis. Unlike a hard rank or a per-dimension variance statistic, effective rank is sensitive to *covariance* — dimensions that are individually non-degenerate but mutually collinear are correctly discounted, which is precisely the failure mode we diagnose.

On the held-out split ($n{=}2530$, $d{=}256$) the biology head is collapsed — $\operatorname{erank}(z_{\text{biology}})$ is $4.4$–$6.0$ depending on modality ($5.3$ fused) — while the sibling identity head stays healthy at $84.3$ (WSI) / $37.5$ (RNA). The collapse point ($\sim\!5$–$6$) coincides with the intrinsic rank of the $50$-D Hallmark target manifold, motivating the use of biology-head effective rank as an **internal fingerprint** of cohort-confounded alignment (C1). Because $\Delta_{\mathrm{spec}}$ shows the collapsed head loses essentially nothing on the confounded benchmark, effective rank flags the pathology that the confounded metric hides — the two are complementary, not redundant.

### Why $\Delta_{\mathrm{spec}}$ is primary and pooled Pearson is confounded

We report $\Delta_{\mathrm{spec}}$ (control-adjusted within-cancer specificity) as the primary number, and treat $r_{\mathrm{pool}}$ as diagnostic only, for three reasons.

1. **Pooled Pearson answers the wrong question.** A model can maximise $r_{\mathrm{pool}}$ by predicting cancer type and reading programme activity off the cross-cancer mean map $c \mapsto \mathbb{E}[y \mid c]$, without resolving any within-type variation. Cancer type is near-trivially decodable from H&E, so $r_{\mathrm{pool}}$ largely rewards a shortcut. Formally, the law of total covariance decomposes $\operatorname{Cov}(\hat y, y)$ into a within-cancer term and a between-cancer term $\operatorname{Cov}\!\big(\mathbb{E}[\hat y\mid c], \mathbb{E}[y\mid c]\big)$; $r_{\mathrm{pool}}$ sums both, $r_{\mathrm{macro}}$ retains only the first. The $\sim\!47\%$ drop from $r_{\mathrm{pool}}$ to $r_{\mathrm{macro}}$ measured across all methods is exactly this between-cancer term.

2. **Within-cancer Pearson alone is still not specific.** Conditioning on cancer removes the mean shift but not the size-dependent technical/compositional floor, which reaches $\approx 0.15$ within-cancer. Only subtraction of the size-matched random-gene null isolates gene-identity-specific signal; this is the Venet/AddModuleScore control lineage applied inside the conditional estimator [Venet 2011; Seurat AddModuleScore; Schmauch 2020].

3. **$\Delta_{\mathrm{spec}}$ is discriminative where the others are not.** Because both confounds are large and method-invariant, $r_{\mathrm{pool}}$ compresses genuinely different models into a narrow $0.33$–$0.35$ band and rewards ones that lean harder on the cohort shortcut. $\Delta_{\mathrm{spec}}$ exposes that the real contested margin is only $\sim\!+0.07$ and reveals it to be flat across methods — the finding that a rank-5–6 head "scores respectably" not on merit but on confound budget. A benchmark that ranks methods by $r_{\mathrm{pool}}$ therefore ranks them largely by how much cohort structure they absorb, which is why we neither optimise nor headline it.

Concretely, all downstream comparisons, ablations, and head-to-head deltas in this paper are reported on $\Delta_{\mathrm{spec}}$ (with $r_{\mathrm{pool}}$, $r_{\mathrm{macro}}$, and the random-set floor shown alongside for transparency), and every claim about method quality is made net of the random-gene control and within cancer type.
