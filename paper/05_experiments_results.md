## Experiments & Results

We organize the empirical section around the paper's three claims. **C1 (mechanism):** the biology head of a matched dual-head model undergoes effective-rank collapse while its sibling identity head stays healthy (T1). **C2 (coupling):** a rank-5–6 head nonetheless scores respectably *only* because the benchmark is dominated by cross-cancer cohort structure and random-gene nulls are already high, leaving a small, method-invariant residual of genuine within-cancer signal (T2, T3, plus retrieval and head-to-head). **C3 (prescriptive):** a per-dimension variance floor cannot repair *rank* collapse, whereas a covariance-decorrelation term is the pre-registered intervention we test in T4.

Unless noted, geometry is measured on the held-out test split ($n=2530$, embedding dimension $d=256$); molecular-prompting metrics are Pearson correlations averaged over 180 MSigDB Hallmark targets and over seeds $\{42,43,44\}$. "Within-cancer" denotes the macro-averaged per-cancer Pearson (`macro_cancer_pearson`), which removes cross-cancer cohort structure by construction. We *concede* to prior art the confound itself (Howard et al. 2021; *Buyer Beware* (Dawood et al. 2024); DECAT), the random-gene control (Schmauch et al., HE2RNA 2020; Venet et al. 2011), the roughly 50% halving of correlation once cohort is controlled (Fu et al. 2020), collapse-as-phenomenon (Andriopoulos et al., NRC1, NeurIPS 2024; Jing et al., ICLR 2022; Bardes et al., VICReg, ICLR 2022), and the retrieval-vs-regression asymmetry (Wang et al., CVPR 2022; Liang et al., FactorCL, NeurIPS 2023). Our contribution is the *coupled* mechanism-plus-benchmark package, and the proposal of biology-head effective rank as an internal fingerprint of cohort confounding.

### T1 — Dual-head effective-rank spectrum (C1)

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

### T2 — Confound decomposition ladder (C2)

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

### T3 — Method-invariance of control-adjusted specificity (C2)

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

### Retrieval (C2, secondary)

Retrieval Recall@$k$ (mean) shows the anchored v2 model does **not** beat the MLP-CLIP baseline; the coupling in C2 is not an artifact of the regression metric.

| Method | Recall@$k$ (mean) | Within-cancer recall |
|---|---|---|
| `mlp_clip` (baseline) | **0.066** | **0.133** |
| `mlp_clip_hardneg` | 0.066 | — |
| `morpheus_v2_anchored` | 0.060 | — |
| `siglip` | 0.059 | — |
| `morpheus_v2_no_anchor` | 0.040 | 0.120 |

The baseline leads on both mean and within-cancer recall. Consistent with the retrieval-vs-regression asymmetry of Wang et al. (2022) and Liang et al. (2023), a low-rank biology head is penalized more sharply by retrieval — which rewards spread — than by per-target regression, so retrieval and specificity agree that v2 buys no real biological gain.

### T4 — VICReg decorrelation ablation **[queued: Lambda full run]**

We register T4 in advance; **no T4 numbers are reported here.** T1 establishes that the biology head loses *rank*, not merely per-dimension variance. A per-dimension variance floor ($\mathrm{relu}(1-\mathrm{std})$, weight 0.01) constrains marginal variances but leaves the off-diagonal covariance — and therefore the rank — unconstrained; it provably cannot prevent rank collapse. The prescriptive fix (C3) couples two changes: **(F-R1)** feed a single L2-normalized biology state to both the Gaussian-NLL heads and the rank-spreaders so they share geometry, and **(F-R2)** add a VICReg/Barlow-style off-diagonal feature-decorrelation term on the normalized 256-D $z_{\text{biology}}$ — the mean-squared off-diagonal of the centered batch covariance, weight $\approx 0.04$, active in the `full` and `programme_only` profiles under a minimum-batch guard.

> **Pre-registered hypothesis (T4).** Adding the covariance-decorrelation term (F-R2), with the shared-geometry coupling (F-R1), will (i) *recover the biology-head effective rank* from the collapsed $\approx 5.3$ (`full_biology`, T1) toward the intact identity-head regime, while (ii) *leaving the control-adjusted within-cancer specificity statistically unchanged* at the method-invariant $\approx +0.07$ Pearson (T3). The prediction that rank rises while specificity does not is the falsifiable core of the fingerprint claim: it would demonstrate that effective rank diagnoses the confound-induced geometry *independently* of task score, so that collapse flags cohort confounding even when benchmark numbers look healthy. The per-dimension variance floor serves as the negative control (predicted: rank unchanged). An optional escalation (F-R3) replaces the neighbour-KL with an RNA-paired biology InfoNCE.
>
> **Exact metric reported by T4.** The primary readout is biology-head effective rank (participation ratio of the normalized singular spectrum) for `full_biology` / `wsi_biology` / `rna_biology`, reported as $\Delta$ effective rank versus the collapsed T1 baseline; the paired secondary readout is control-adjusted within-cancer specificity (real minus random-gene null, `macro_cancer_pearson`), tested for the pre-registered null of no change from $\approx +0.07$. Both are reported per seed $\{42,43,44\}$, with the variance-floor arm as negative control.
