# SOTA / Novelty Experiment Review — moving MORPHEUS from "diagnostic" to a genuine contribution

**Reviewer lane:** research (deep-dive)
**Date:** 2026-07-28
**Scope:** What is the strongest, feasible additional experiment or reframing that converts the effective-rank-fingerprint diagnostic into a novel + potentially SOTA contribution on WSI→molecular / multimodal cancer representation?

---

## 0. What the paper currently is (and its central vulnerability)

The paper is a clean **negative/diagnostic** result:

- **C1:** biology head collapses to eff-rank ~5–6/256 (identity sibling ~84) — a matched dual-head fingerprint of aligning to a low-rank (intrinsic ~5–6) 50-D Hallmark manifold.
- **C2:** the collapsed head still "scores respectably" only because the benchmark is ~46–49% cross-cancer cohort structure and the random-gene null is already ~0.30–0.32, leaving a **method-invariant ~+0.07** within-cancer, control-adjusted specificity.
- **C3:** variance floor provably can't fix rank; off-diagonal decorrelation can (synthetic 10.3→21.2). T4 real-data ablation queued on A100.

**The vulnerability that blocks acceptance as anything but a diagnostic:** the pre-registered T4 hypothesis *predicts* that fixing collapse (rank ↑) leaves specificity unchanged at +0.07 — i.e. the fix is **decoupled from any payoff**. A reviewer will say: "You diagnosed a geometric artifact that, by your own pre-registration, does not matter for the task. Why should I care?" Every proposal below is designed to answer exactly that: **show that fixing collapse (or exploiting the fingerprint) buys something real that current SOTA does not have.**

A second structural point that the proposals exploit: the paper's critique is scoped to the **cross-cancer bulk** regime. The regimes where morphology→molecular signal is genuinely real (within-slide spot-level, and confound-controlled downstream clinical endpoints) are *untouched* by the critique — and are exactly where a SOTA claim can live without contradicting C2.

---

## 1. Direct answers to the two sub-questions

### 1a. "Does fixing collapse actually improve a real downstream task beyond the confound?"

**Unknown, and this is the single most important gap.** The literature is split and that split is the opportunity:

- Multiple 2024–25 results argue dimensional/rank collapse *is* a downstream bottleneck and effective rank is a good proxy for downstream quality (LDReg, ICLR 2024; WERank 2024; "Preventing Dimensional Collapse in SSL," NeurIPS 2024; Matrix-Information-Theory SSL). This predicts a **positive** payoff for the fix.
- But the paper's own T4 pre-registration predicts **no** payoff on the confounded prompting metric — which is consistent, because that metric is confound-saturated and cannot see representation quality (their own C2).

These are only contradictory if you measure payoff on the *wrong* endpoint. The prompting metric is the wrong endpoint. The right endpoints are **confound-controlled downstream clinical tasks** (survival, subtype, mutation) and **high-rank molecular targets** — neither of which is confound-saturated. **No experiment in the current draft measures either.** Proposal 1 and 2 close this.

### 1b. "What recent (2024–2026) methods could we beat or combine with?"

| Method | Type | Relevance | Beat or combine |
|---|---|---|---|
| **PLUTO-4G** (r=0.427 HEST), **H-Optimus-1** (0.422), **UNI2-H** (0.413), **Virchow-2** (0.396) | Pathology FMs, HEST expression leaderboard | The public leaderboard to beat for morphology→expression | **Beat** (Prop 2) by fixing target rank + decorrelation on top of H-Optimus features |
| **SEQUOIA** (Nat Commun 2024, linearized attention, UNI features) | Bulk RNA from WSI | Direct competitor for WSI→bulk regression | **Beat/combine** — add decorrelated head |
| **MINT** (2026, spatial-transcriptomics-supervised FM) / **STPath** (npj Digital Med 2025) | Molecular-supervised pathology FMs | These *do exactly* the molecular-supervision the paper warns collapses — perfect foil | **Combine/critique**: show their heads collapse too, then fix |
| **POMP** (IJCAI 2025), **CPathomic**, **MurreNet** (MICCAI 2025), **M4Survive** (C-index 81.27), **Multimodal Prototyping** (MICCAI 2024), **THREADS/TANGLE** | Multimodal path+omics survival, contrastive alignment | SOTA survival with contrastive alignment — none report effective-rank or confound-robustness | **Beat on confound-robust C-index** (Prop 1); **audit** their alignment heads (Prop 3) |
| **HESCAPE** (2025), **DECAT** (2026), **Buyer Beware** (NBME 2026) | Confound audits | Prior art on the confound — must differentiate | **Differentiate**: erank is internal, label-free, per-model (Prop 3) |
| **VICReg / Barlow Twins / LDReg / WERank** | Decorrelation regularizers | The fix family | **Combine** — LDReg (local dim reg) and WERank (weight-space) are stronger F-R2 variants worth benchmarking against plain off-diag covariance |

---

## 2. Three ranked experiment proposals

### 🥇 Proposal 1 — Confound-controlled downstream causal test: does decorrelation transfer to survival / subtype / mutation? (TOP)

**One-line thesis:** Convert the fingerprint from a *diagnostic* into a *training principle with payoff* by showing the rank-restored biology embedding beats the collapsed head — and beats frozen UNI/Virchow/H-Optimus features — on real clinical endpoints under **site-stratified** evaluation, where the +0.07 confound ceiling does not apply.

**Why this is #1:** It directly answers the reviewer's killer question (1a), it reuses the exact infra already on disk (TCGA-CDR survival, MC3 MAF mutations, GISTIC CNV, PAM50/molecular subtypes, H-Optimus patches, the trained V2 heads), and it de-risks the whole paper: *whatever the outcome is publishable.*

**Design:**
1. Freeze three embeddings from the same model: `z_biology` (collapsed, erank ~5), `z_biology_decorrelated` (post-F-R1+F-R2, erank restored), `z_identity` (healthy control ~84). Add external baselines: frozen UNI2, Virchow2, H-Optimus-1 slide embeddings (mean/ABMIL pooled).
2. Three downstream endpoints, each with a linear/Cox probe on frozen features:
   - **Survival** — Cox C-index, TCGA-CDR OS/PFI, per-cancer.
   - **Molecular subtype** — PAM50 (BRCA), TCGA subtypes (multiple), macro-AUROC.
   - **Driver mutation** — top-N recurrent mutations from MC3 (e.g. TP53, PIK3CA, KRAS, BRAF), AUROC.
3. **The load-bearing design choice:** every probe is evaluated under **site-preserved / site-stratified splits** (Howard 2021 protocol) *and* a naive-random split, reported side by side. The naive split re-imports the confound; the site-stratified split isolates real signal. This is the axis on which erank should matter.
4. **Primary hypothesis (the novel, falsifiable claim):** decorrelation-restored rank *does not* help on naive splits (confound-saturated, matches T4) but *does* help on site-stratified splits (confound-controlled) — i.e. **fixing collapse buys confound-robustness specifically.** That is a genuinely new, mechanistically-grounded claim connecting representation geometry → confound-robust generalization.

**What to beat/combine:** frozen UNI2/Virchow2/H-Optimus (the FM features everyone uses); as a stretch, wire the decorrelated head into a survival head and compare C-index to POMP / MurreNet / M4Survive numbers on shared TCGA cohorts.

**Expected payoff:**
- *Best case:* decorrelation delivers a consistent C-index/AUROC gain under site-stratified eval that collapsed heads and frozen FMs do not — reframes the paper as "a label-free geometric criterion that predicts and produces confound-robust multimodal representations," a SOTA-relevant training principle. High-impact.
- *Worst case:* no gain even on site-stratified splits → **strengthens** C2/C3 into a sharper claim ("collapse is genuinely orthogonal to downstream utility; the fix is cosmetic") — still a clean, defensible paper.
- Either way the "so what?" objection dies.

**Feasibility:** High. All data present on disk (`TCGA-CDR-SupplementalTableS1.xlsx`, `mc3...maf.gz`, GISTIC CNV, clinical followup TSV, H-Optimus patch store). Only new code is linear/Cox probes + site-stratified splitter. Fits the queued A100 window; no retraining strictly required for the frozen-probe version (decorrelated head needs the T4 run that is already queued).

**Risk:** Medium. If decorrelation doesn't help even under site-stratification, the SOTA angle weakens — but the diagnostic strengthens, so downside is bounded.

---

### 🥈 Proposal 2 — Break collapse at the source with high-rank targets, then compete on the HEST leaderboard

**One-line thesis:** C1 says erank tracks *target* rank. Turn that mechanism into a design law: swap the intrinsic-rank-~5 Hallmark target for a **high-rank per-gene target** (highly-variable genes), which should *not* collapse — then benchmark on the public **HEST-Benchmark** (spot-level, within-slide), the regime the paper's cross-cancer critique explicitly does not cover, and try to pass **PLUTO-4G's r=0.427 SOTA**.

**Why #2:** Highest SOTA ceiling (a real public leaderboard with a live number to beat) and it validates C1 as a *prescriptive law* rather than a bug report. It's #2 only because it needs a new dataset (HEST spatial transcriptomics) and a retrain, i.e. more work and more risk than Prop 1.

**Design:**
1. **Target-rank sweep (mechanism → law):** retrain the biology head against targets of increasing intrinsic rank — 50 Hallmark (rank ~5) → 50 HVG (HEST protocol) → 1000 HVG → per-gene. Plot biology-head erank vs target intrinsic rank. Prediction (from C1/NRC1): erank rises monotonically with target rank; collapse disappears when the target is high-rank. This is a strong, clean confirmation of the mechanism and it is currently *asserted but not swept*.
2. **HEST-Benchmark head-to-head:** on the 9 HEST tasks (50 HVG, spot-level, ridge-after-PCA protocol), evaluate H-Optimus features + (a) plain regression head, (b) +F-R2 decorrelation, (c) high-rank-target-trained biology head. Compare mean Pearson to PLUTO-4G (0.427), H-Optimus-1 (0.422), UNI2-H (0.413), Virchow-2 (0.396).
3. **Consistency with the thesis:** frame HEST (within-slide, spot-level) as the *un-confounded* complement to the cross-cancer bulk critique — beating SOTA here does not contradict C2, it completes the story ("morphology→molecular is real where cohort structure is removed by construction").

**What to beat/combine:** PLUTO-4G / H-Optimus-1 on HEST; combine decorrelation (F-R2) or LDReg/WERank on top of the strongest frozen FM.

**Expected payoff:**
- *Best case:* a SOTA or near-SOTA HEST number *and* a clean target-rank→erank scaling law → the paper now has both a mechanism (C1 as law) and a leaderboard win. Strong, publishable-at-a-good-venue upgrade.
- *Likely case:* competitive-but-not-SOTA Pearson + a decisive erank-vs-target-rank curve → mechanism confirmed, SOTA claim softened to "matches SOTA while explaining why."

**Feasibility:** Medium. HEST-1k is public and standardized (mahmoodlab, NeurIPS 2024) with a fixed protocol, so evaluation is turnkey. Cost is downloading HEST + retraining heads on the new targets. Fits A100 but is a bigger job than Prop 1.

**Risk:** Medium-high. Beating PLUTO-4G (a frontier FM) with an H-Optimus-based head is not guaranteed; the *scaling-law* result is nearly certain even if the leaderboard win isn't, so partial payoff is protected.

---

### 🥉 Proposal 3 — The fingerprint as a label-free, cross-model confound-robustness predictor (reframing)

**One-line thesis:** Elevate the differential erank fingerprint from "a thing that happened in our one model" to a **general, unsupervised model-selection criterion**: across a panel of public pathology FMs and multimodal alignment methods, does biology-head / alignment-head effective rank *predict* confound-robust downstream generalization (the site-stratified gap from Prop 1) — with no labels and no held-out confounder?

**Why #3:** It's the cheapest reframing and it maximally differentiates from DECAT/HESCAPE/Buyer-Beware (all of which need the confounder named/held-out; erank needs neither). It's #3 because it is a *methodological* contribution rather than a SOTA number, so its ceiling is "useful new metric" not "beats leaderboard."

**Design:**
1. Assemble a panel: UNI2, Virchow2, H-Optimus, DINOv2, plus 2–3 molecular-supervised/aligned methods that expose an alignment head (MINT, STPath, POMP, CPathomic, TANGLE/THREADS if reproducible).
2. For each, compute the erank of the molecular/alignment head (and the differential vs an identity/instance head where available) — purely from a forward pass on unlabeled held-out slides.
3. Correlate erank (and the differential) against the **site-stratified minus naive downstream gap** from Prop 1's endpoints. Claim: low biology-head erank predicts a large confound-robustness gap (i.e. flags models that will over-rely on cohort structure) — a label-free early-warning metric.
4. Bonus: erank as an **early-stopping / model-selection** signal during training (matrix-info-theory SSL shows erank rises through training — does the erank peak predict the confound-robust optimum?).

**What to beat/combine:** DECAT / HESCAPE as the differentiator baseline — show erank recovers their confound verdicts without their required inputs (named confounder, paired null, held-out site).

**Expected payoff:**
- A validated, unsupervised, per-model metric that predicts confound-robustness → a genuinely novel *tool* the community can adopt, and a much stronger C1 (fingerprint generalizes beyond MORPHEUS). Moderate-high methodological impact; low SOTA impact.

**Feasibility:** Medium. Depends on reproducing/loading several external models; erank computation itself is trivial. Some alignment methods may be hard to run.

**Risk:** Medium. If erank does *not* correlate with the confound-robustness gap across models, C1's generality is undercut — but that itself is informative and bounds the fingerprint's scope.

---

## 3. Recommended sequencing

1. **Do Prop 1 first** — it reuses on-disk data, rides the queued A100/T4 run, and is the fastest path to killing the "so what?" objection. It also *produces* the site-stratified downstream gap that Prop 3 needs as its target variable.
2. **Then Prop 2** if the target-rank mechanism holds and A100 time remains — it's the SOTA upside.
3. **Prop 3** as the integrative reframing once Prop 1's gaps exist to correlate against.

Props 1→3 compose into a single arc: *diagnose (C1) → fix (C3) → show the fix buys confound-robustness (Prop 1) → generalize the fingerprint into a label-free selection metric (Prop 3)*, with Prop 2 as the SOTA leaderboard flourish.

---

## Sources

- [Preventing Dimensional Collapse in SSL (NeurIPS 2024)](https://proceedings.neurips.cc/paper_files/paper/2024/file/ad7922fd4650f8aba5d8b067e622ca84-Paper-Conference.pdf)
- [LDReg: Local Dimensionality Regularized SSL (ICLR 2024)](https://arxiv.org/pdf/2401.10474)
- [WERank: Rank Degradation Prevention via Weight Regularization](https://arxiv.org/pdf/2402.09586)
- [Matrix Information Theory for Self-Supervised Learning](https://arxiv.org/pdf/2305.17326)
- [HEST-1k dataset + HEST-Benchmark (NeurIPS 2024)](https://github.com/mahmoodlab/hest)
- [A Large-Scale Benchmark of Cross-Modal Learning for Histology and Gene Expression](https://arxiv.org/html/2508.01490v1)
- [PLUTO-4: Frontier Pathology Foundation Models (HEST SOTA r=0.427)](https://arxiv.org/pdf/2511.02826)
- [MINT: Molecularly Informed Training with Spatial Transcriptomics Supervision](https://arxiv.org/pdf/2603.07895)
- [STPath generative FM for spatial transcriptomics + WSI (npj Digital Medicine 2025)](https://www.nature.com/articles/s41746-025-02020-3)
- [SEQUOIA — digital gene-expression profiling with linearized attention (Nat Commun 2024)](https://www.nature.com/articles/s41467-024-54182-5)
- [POMP: Pathology-omics Multimodal Pre-training (IJCAI 2025)](https://www.ijcai.org/proceedings/2025/869)
- [MurreNet (MICCAI 2025)](https://papers.miccai.org/miccai-2025/0624-Paper0057.html)
- [M4Survive: Multi-Modal Mamba survival (C-index 81.27)](https://arxiv.org/html/2503.10057v1)
- [Multimodal Survival Modeling in the Age of Foundation Models](https://arxiv.org/html/2505.07683v2)
- [BatMan: Mitigating Batch Effects via Stratification for Survival](https://arxiv.org/pdf/2209.03902)
- [Comprehensive benchmark of image FMs for gene expression from WSI (bioRxiv 2026)](https://www.biorxiv.org/content/10.64898/2026.03.02.709012v1.full)
