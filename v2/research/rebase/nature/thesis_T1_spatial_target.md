# Thesis T1 — "The target is wrong: spatial transcriptomics as correctly-localized ground truth"

**Scout verdict: the thesis as stated is empirically dead on arrival. An inverted version of it is
the strongest, most executable result in the set — but it is a benchmark-audit contribution, not a
Nature paper.**

Date: 2026-07-29. Search constraint: WebSearch exhausted; all citations below verified this session
via WebFetch against arXiv abstract pages, NCBI E-utilities (esearch/esummary/efetch), GitHub, and
HuggingFace dataset pages. Anything not verified is explicitly marked COULD-NOT-VERIFY and is
**not** counted as prior art.

---

## 0. The claim under test

> WSI→molecular prediction is stuck (~+0.07 genuine signal, method-invariant) because BULK RNA IS A
> MIS-SPECIFIED TARGET — it averages the whole tumor while morphology is local. A model
> trained/evaluated against spatially-localized targets breaks the ceiling that every bulk-target
> method ties at; the bulk ceiling is a property of the TASK, not the models.

The claim decomposes into two sub-claims that must be judged separately:

- **(a) "Bulk averages; morphology is local."** TRUE, uncontroversial, and *already the founding
  motivation of the entire image→ST subfield*. This is the first paragraph of HEST-1k, STimage-1K4M,
  ST-Net's descendants, and every paper below. It is not a novel observation and cannot carry a paper.
- **(b) "Spatial targets break the ceiling."** **NOT SUPPORTED. The published evidence points the
  other way.** See §2. This is the load-bearing half of the thesis and it fails.

---

## 1. Has this been done? — the 5 closest works

### 1.1 HEST-1k / HEST-Benchmark — *the proposed benchmark already exists and is two years old*
Jaume, Doucet, Song, Lu, Almagro-Pérez, Wagner, Vaidya, Chen, Williamson, Kim, Mahmood.
**NeurIPS 2024 Spotlight.** arXiv:2406.16192. VERIFIED (arXiv abs page; GitHub `mahmoodlab/HEST`;
HuggingFace `MahmoodLab/hest`).

What it did, precisely: assembled 1,229 spatial-transcriptomic profiles (1,276 in the current HF
release) each linked to a WSI, from 153 public+internal cohorts (180 per the current HF card), 26
organs, 2 species, 367 cancer samples across 25 cancer types; yielded **2.1M expression-morphology
pairs and >76M nuclei**. HEST-Benchmark formalizes the task as **nine per-organ tasks** (IDC, PRAD,
PAAD, SKCM, COAD, READ, CCRCC, LUNG, LYMPH_IDC), each predicting the **top-50 highly-variable genes
from a 112×112 μm H&E patch centered on each ST spot**, as multivariate regression, with Ridge (+PCA
256) or Random Forest on frozen patch-encoder features, scored by Pearson.

**This is the deliverable T1 proposes.** "Spot-level ground truth aligned to morphology, open access,
standardized, pan-cancer" is HEST. The reframing is not available.

### 1.2 The HEST-Benchmark leaderboard — *the ceiling reproduces at spot level*
VERIFIED (GitHub `mahmoodlab/HEST` README, fetched this session).

| Encoder | Mean Pearson, 9 tasks |
|---|---|
| H-Optimus-1 | **0.4229** |
| GenBio-PathFM | 0.4197 |
| H-Optimus-0 | 0.4150 |
| ... (25 FMs evaluated) | ... |
| **ResNet50 (ImageNet, 2015)** | **0.3252** |

Read this table the way MORPHEUS reads its own: **spread among the top pathology foundation models
≈ 0.008; best-FM-minus-ResNet50 = +0.098.** That is the *same shape* as MORPHEUS's bulk finding —
a small, method-invariant delta over a naive baseline — one decimal place up. Twenty-five foundation
models, trained on 100M+ pathology patches, tie within a hundredth of each other and beat a 2015
ImageNet CNN by ~0.1.

Corroborating: **"From histology to spatial transcriptomics: establishing a lightweight single-patch
baseline"** (BMC Bioinformatics 2026; PMID 42151778; DOI 10.1186/s12859-026-06447-7; VERIFIED via
efetch) reports **EfficientNet-B0, 5.3M parameters, max PCC 0.310** on liver and 50 genes at
PCC ≥ 0.30 — i.e. a model ~1/60th the size of a pathology FM sits essentially at the HEST ResNet50
floor. The authors further note accuracy "correlates with spatial organization patterns rather than
transcript abundance," which is the composition hypothesis in §2.3 stated by someone else.

### 1.3 HESCAPE — *spatial supervision does not rescue expression prediction; the confound just moves*
Gindra, Palla, Nguyen, Wagner, Tran, Theis, Saur, Crawford, Peng. arXiv:2508.01490 (ICCV-W 2025).
VERIFIED (arXiv abs page).

Curated pan-organ ST benchmark, **6 gene panels, 54 donors**. Findings: (i) the *gene-expression*
encoder, not the image encoder, drives cross-modal alignment; (ii) contrastive cross-modal
pretraining **improves gene-mutation classification but degrades direct gene-expression prediction
relative to baseline encoders**; (iii) **batch effects are named as the central obstacle** to
cross-modal integration.

Why this matters to T1: it is the direct empirical refutation of "move to spatial and the signal
appears." Moving to spatial targets does not remove the non-biological structure; it swaps
*cross-cancer cohort identity* for *per-slide / per-donor batch*. HiST (arXiv:2606.14251, ICML 2026;
previously logged in `lit/l11_benchmarks_confound.md`) independently treats this as first-class by
adding a **per-slide calibration token** — i.e. the field already concedes slide identity is a
dominant nuisance in the spot-level task.

### 1.4 "Benchmarking the translational potential of spatial gene expression prediction from histology"
Wang C., Chan A.S., Fu X., et al. **Nature Communications 2025**; PMID 39934114;
DOI 10.1038/s41467-025-56618-y. VERIFIED (esummary + efetch).

Eleven image→ST methods evaluated on five spatially-resolved-transcriptomics datasets plus **external
TCGA validation**, across five tiers: within-image accuracy, cross-study generalizability,
translational/survival impact, usability, compute. Findings: within-image spot accuracy is a **poor
predictor** of cross-study transfer or survival utility; architectural complexity does not track
translational value; simple CNNs (DeepPT-class) match or beat elaborate architectures.

Why this matters to T1: this is the closest published test of "does the spatial target buy real
payoff," and the answer is *no downstream*. Per the earlier MORPHEUS lane note
(`A1_ledger/a_prior_benchmarks.md`), predicted-expression survival models land at C-index ~0.52–0.58
against a bulk-RNA baseline of ~0.57–0.58. The spatially-localized target does not beat bulk where it
counts.

### 1.5 SEPAL (+ SpaRED / SpaCKLE) — *the central methodological move is already published*
SEPAL: Mejia, Cárdenas, Ruiz, Castillo, Arbeláez. arXiv:2309.01036 (v3, Jan 2024). VERIFIED.
SpaRED/SpaCKLE: arXiv:2407.13027 and arXiv:2505.02980; **Medical Image Analysis 2025**, PMID 40885036,
DOI 10.1016/j.media.2025.103754. VERIFIED.

SEPAL **directly supervises relative differences with respect to mean expression**, and **proposes a
refined benchmark that restricts the prediction variables to only those genes with clear spatial
patterns**. That is exactly the "subtract the trivially-predictable mean, evaluate only what is
genuinely spatial" control T1 would need to invent. It exists, on two human breast-cancer datasets,
since 2023. SpaRED curates 26 ST datasets and shows that **target-side dropout** is severe enough
that completing it (SpaCKLE) cuts MSE >82.5% and **reshuffles the leaderboard** — i.e. much prior
image→ST benchmarking was measuring target noise.

### 1.6 HistoPrism — *the sharpest collision: pan-cancer + spatial + pathway-level, already at ICLR*
Hu, Zeng, Bhasker, Kather, Speidel. arXiv:2601.21560, **accepted ICLR 2026**. VERIFIED (arXiv abs).

Transformer predicting **spatial gene expression from H&E pan-cancer**, whose stated key contribution
is **"a pathway-level benchmark, shifting assessment from isolated gene-level variance to coherent
functional pathways,"** with claimed **strong pan-cancer generalization**.

This is the closest thing in the literature to "reframed task + spatial target + pathway/Hallmark
readout + pan-cancer." **Read the full paper before committing to T1 in any form** — if HistoPrism
also reports variance/mean decomposition, the remaining whitespace in §4 shrinks materially.

### Also relevant (verified, secondary)
- **STimage-1K4M** (Chen, Zhou, Wu, Zhang, Li, Li; arXiv:2406.06393): 1,149 ST slides, **4,293,195
  sub-tile/expression pairs**, 15,000–30,000-dim expression per tile.
- **STFlow** (Huang, Liu, Babadi, Jin, Ying; arXiv:2506.05361, **ICML 2025**): whole-slide flow
  matching modeling the joint spot distribution; **>18% relative improvement over pathology FMs** on
  HEST-1k and STimage-1K4M. Note: +18% *relative* on ~0.42 is ~+0.075 absolute — again the same order.
- **DiffBulk** (IEEE TMI 2026; PMID 42048193; DOI 10.1109/TMI.2026.3688322): diffusion-based ST
  prediction on Xenium.
- Bulk-target reference points for the comparison: **HE2RNA** (Nat Commun 2020, PMID 32747659),
  **tRNAsformer** (Commun Biol 2023, PMID 36949169), **SEQUOIA** (Nat Commun 2024, PMID 39543087;
  7,584 tumors / 16 cancer types, validated on 1,368 tumors in two independent cohorts),
  **ENLIGHT-DeepPT** (Nature Cancer 2024, PMID 38961276; 16 TCGA cohorts).

### COULD-NOT-VERIFY
- **TANGLE / "Transcriptomics-guided Slide Representation Learning" (Jaume et al., CVPR 2024).**
  Recalled but **not verifiable this session** (arXiv export API returned persistent 429; CVF
  OpenAccess returned 403; Semantic Scholar API returned persistent 429). If it exists as recalled it
  is the bulk-RNA-as-slide-supervision precedent and directly relevant. **Do not cite until verified.**
- **Semantic Scholar search was entirely unavailable** (429 on every call). A second sweep should redo
  the S2 queries; the arXiv full-text search and PubMed coverage below may therefore under-sample
  CS-venue preprints.

---

## 2. Is the core claim true? — evidence on "spot-level > bulk in morphology-predictability"

### 2.1 Direct answer to the posed question
**No published work found this session reports a controlled comparison of bulk vs spot-level targets
on matched cohorts.** The premise is *untested*, not established. That is a genuine whitespace — but
it is whitespace because the comparison is hard to make meaningful (§2.2), not because nobody thought
of it.

### 2.2 The numbers do not support the thesis, and the naive comparison is invalid
The tempting comparison is: HEST spot-level **0.42** vs MORPHEUS within-cancer bulk **~0.19**.
Therefore spatial wins by 2×. **This is an apples-to-oranges comparison of variance denominators, and
it will be the first thing a reviewer destroys.**

- **Spot-level Pearson is computed over spots spanning tumor / stroma / necrosis / immune infiltrate
  within one slide.** The dominant variance component is **cell-type composition**, which H&E reports
  almost by definition. Predicting "this spot is lymphoid, that one is stromal" is closer to a
  segmentation task than a molecular-state task.
- **Bulk within-cancer Pearson is computed over patients, after composition has been averaged out.**
  What remains is precisely the hard residual.

So the higher spot-level number plausibly reflects an *easier* target, not a *correctly specified*
one. The thesis assumes correct localization; the data are equally consistent with easier variance.

### 2.3 The method-invariance signature reproduces at spot level — this is the decisive fact
| | Bulk (MORPHEUS, established) | Spot-level (HEST-Bench, published) |
|---|---|---|
| Best-vs-naive-baseline delta | **~+0.07** (random-gene-null-adjusted, within-cancer) | **~+0.098** (H-Optimus-1 0.4229 − ResNet50 0.3252) |
| Spread among strong methods | method-invariant | **~0.008** across the top 25 pathology FMs |
| Tiny/old model performance | baselines tie | EfficientNet-B0 @ 5.3M params → 0.310 |
| Named dominant nuisance | cross-cancer cohort identity (~46–49%) | per-slide/donor batch (HESCAPE); slide-calibration token (HiST) |

**Switching from bulk to spot-level does not break the tie. It reproduces it.** The ceiling survives
the change of target modality. That is the finding — and it is the opposite of the stated hypothesis.

### 2.4 What is genuinely true in the neighbourhood
- The *absolute* signal is higher spot-level (0.42 vs 0.19). Something real is there.
- Nobody has decomposed it. HEST-Bench documents **no mean-expression / constant predictor baseline**
  and no random-gene-panel null (confirmed: the benchmark tutorial specifies only Ridge+PCA and RF).
  SEPAL's mean-subtraction and spatially-patterned-gene restriction is the only published gesture in
  this direction, and only on two breast datasets.
- Therefore: **the decomposition of the spot-level number into (per-slide mean) + (cell composition) +
  (composition-independent molecular state) is unpublished.** That is the real opportunity.

---

## 3. Open-access data: what exists, real sizes, real limitations

| Dataset | Access route | Real size | License | Contents |
|---|---|---|---|---|
| **HEST-Bench** | HF `MahmoodLab/hest-bench` (gated) | **42.2 GB**; 259 samples (187 train / 72 test) | CC BY-NC-SA 4.0 | 9 tasks; **2–24 samples per task**; 112 μm patches + top-50 HVG `.h5ad` |
| **HEST-1k (full)** | HF `MahmoodLab/hest` (gated: account + accept terms) | **2.01 TB** total (~1 TB core) | CC BY-NC-SA 4.0 | 1,276 ST profiles, 180 cohorts, 26 organs, 2 species, 25 cancer types; pyramidal TIFF WSIs, `.h5ad` expression, tissue masks, **224×224 patches around spots**, **CellViT nuclei segmentation**, Xenium transcript parquet |
| **STimage-1K4M** | HF `jiawennnn/STimage-1K4M`; GitHub `JiawenChenn/STimage-1K4M` | 1,149 slides; **4,293,195** tile/expression pairs | code MIT; **per-source data licenses must be checked individually** | ST / Visium / VisiumHD; human + mouse; 15k–30k genes per spot |
| **SpaRED (+SpaCKLE)** | arXiv:2407.13027 / 2505.02980 repos | 26 curated ST datasets | per-source | dropout-completed ST for fair benchmarking |
| **TCGA bulk + WSI** | **already on disk** | 6,192 patients, 32 cancers | — | H-Optimus patch features (uncapped) + BulkFormer bulk RNA, held-out-cancer split |

**Feasibility headline: the entry point is 42 GB, not 2 TB.** HEST-Bench is small, gated-but-free, and
its leaderboard is computed on **H-Optimus** features — the exact encoder MORPHEUS already has
extracted for 6,192 TCGA patients. CellViT nuclei segmentations **ship with HEST-1k**, so the
composition-adjustment covariate in §4 requires no new model.

**Real limitations, stated honestly:**
1. **Tiny per-task n.** 2–24 samples per HEST task. Patient-level CV over ≤24 samples gives wide
   confidence intervals. *Before promising any "these methods tie" claim, compute the CI width — a tie
   you cannot distinguish from noise is not a finding.* This is the top practical risk.
2. **No matched bulk+ST cohort at scale.** There is no large public cohort with **both** bulk RNA and
   spot-level ST on the same patients. The clean controlled bulk-vs-spot comparison the thesis needs
   cannot be run directly; it must be done as a cross-dataset comparison of *control-adjusted deltas*,
   which is weaker and must be argued carefully.
3. **Platform heterogeneity is severe.** Legacy ST (100 μm spots), Visium (55 μm, ~1–10 cells),
   VisiumHD (2 μm bins), Xenium (**targeted 300–500 gene panel, not transcriptome-wide**). "Spot-level
   ground truth" is not one measurement — HESCAPE's "6 gene panels" is the same problem.
4. **Spots are not cells.** A 55 μm Visium spot is itself a mini-bulk average. The thesis's own logic
   ("averaging is the problem") applies recursively; the difference from bulk is one of degree.
5. **Target-side dropout / zero-inflation** is severe enough to reshuffle leaderboards (SpaRED).
6. **CC BY-NC-SA 4.0** on HEST — non-commercial and share-alike. Fine for academic publication;
   note it if any industry collaboration is contemplated.

---

## 4. Highest honest claim, and what falsifies it

### 4.1 The claim to abandon
> "Spatial targets break the bulk ceiling."

Abandon this. A reviewer who opens the HEST README sees 25 FMs within 0.008 of each other and +0.098
over a 2015 ResNet50, and the paper is over in one sentence.

### 4.2 The claim to make instead — **target-invariance of the ceiling**

> **The ~+0.07–0.10 method-invariant ceiling in H&E→molecular prediction is TARGET-INVARIANT.** It
> holds identically whether the target is bulk RNA over 6,192 TCGA patients across 32 cancers, or
> spot-level spatial transcriptomics over 1,276 HEST samples across 26 organs. The apparently larger
> spot-level signal (r ≈ 0.42) is not better-localized molecular state: it is dominated by per-slide
> mean and cell-type composition — quantities a 5.3M-parameter CNN and a 2015 ImageNet ResNet50
> already capture — and it collapses to the same +0.07–0.10 once you subtract the per-slide mean,
> restrict to spatially-patterned genes, and adjust against a random-gene-panel null. **The bottleneck
> is therefore not target mis-specification. It is the information content of H&E morphology about
> molecular state beyond cell composition — and the field's leaderboards are measuring composition.**

This inverts the thesis and is *stronger* than it, because:
- it is a law-like statement tested on two independent target modalities and two independent cohorts;
- it *explains* both the bulk tie and the HEST tie with one mechanism;
- it ships a concrete diagnostic the field demonstrably lacks (HEST has no mean/null baseline);
- it converts a negative result into an evaluation-paradigm contribution — the axis `REBASE_THESIS`
  §5 and `RESEARCH_PATHS` already rank as most defensible ("an evaluation contribution is hardest to
  scoop").

### 4.3 Deliverable ladder
- **Tier 1 — near-certain (weeks).** Port MORPHEUS's exact control-adjusted decomposition
  (random-gene-panel null, per-slide/per-cohort mean subtraction, held-out-cancer protocol) onto
  HEST-Bench using the H-Optimus features already on disk. Publish the **control-adjusted HEST
  leaderboard**. Expected outcome: the leaderboard compresses or reorders.
- **Tier 2 — likely (the damning number).** Regress spot expression on **CellViT nuclei counts and
  types shipped with HEST-1k** as a composition-only baseline. Report the pathology-FM delta *over
  composition*. If that delta is ≈0, the field's premier molecular benchmark is a nuclei-counting
  benchmark. That is a real, quotable, previously-unreported result.
- **Tier 3 — speculative (the only route to a high-tier venue).** Identify the gene/pathway subset
  carrying **composition-independent** morphological signal; test whether it transfers to held-out
  cancers *and* to bulk TCGA. If that subset is small, stable, and biologically coherent
  (proliferation / hypoxia / EMT / stress), the claim upgrades from "the benchmark is broken" to
  "morphology encodes exactly these molecular programs and no others" — a positive biological
  finding. **Caveat: the Nat Commun 2025 translational benchmark already suggests such transfer is
  weak, so price Tier 3 as unlikely.**

### 4.4 Falsifiers (pre-register these)
| # | Falsifies | Condition |
|---|---|---|
| F1 | **Target-invariance claim** | After mean-subtraction + spatially-patterned-gene restriction + random-panel null + composition adjustment, HEST-Bench retains a delta materially >+0.10 **and** that delta separates methods (top FM clearly > ResNet50 / EfficientNet-B0). Then spatial targets really do carry more, and the original T1 is right. |
| F2 | **Composition explanation** | A CellViT-nuclei-composition regressor scores far below the pathology FMs on the adjusted metric — i.e. the FMs know something beyond composition. |
| F3 | **Novelty** | HistoPrism (ICLR 2026, arXiv:2601.21560) or a HEST v2 already reports the composition-adjusted / mean-decomposed leaderboard. HistoPrism explicitly claims to address "prior variance-based assessment limitations." **Read it in full before committing.** |
| F4 | **Statistical viability** | Per-task n = 2–24. If bootstrap CIs on the per-task Pearson are wider than the ~0.008 inter-method spread, "the methods tie" is unprovable from HEST alone, and Tier 1 must be re-scoped (pool tasks, or add STimage-1K4M / SpaRED for power). **Check this first — it is cheap and it gates everything.** |
| F5 | **Whole thesis** | Someone publishes the matched bulk+ST controlled comparison on a real cohort with a clean result either way. No such cohort is currently public at scale, so this is low near-term risk. |

---

## 5. Ratings

| Axis | Score | Reasoning |
|---|---|---|
| **Novelty** | **2 / 5** | As literally stated: **1–2**. The dataset (HEST, 2024), the pan-cancer spatial task (HistoPrism, ICLR 2026), the mean-subtraction move (SEPAL, 2023), the spatially-patterned-gene restriction (SEPAL), the target-noise critique (SpaRED, 2024–25), the batch-confound critique (HESCAPE, 2025), and the translational-transfer critique (Nat Commun 2025) are **all published**. The one unclaimed piece is the *control-adjusted cross-target audit* (§4.2), which lifts the reframe to ~3. Scored **2** for the thesis as posed. |
| **Feasibility (open data)** | **5 / 5** | Very likely the most executable thesis in the set. HEST-Bench is **42 GB**, free, gated-only. Its leaderboard already uses **H-Optimus**, which MORPHEUS has extracted for 6,192 patients. **CellViT nuclei ship with the dataset**, so the composition baseline needs no new model. MORPHEUS already owns the exact control machinery (random-gene null, cohort decomposition, held-out-cancer). A100-40GB is ample. Timeline: weeks. |
| **Ceiling** | **2 / 5** | Best realistic outcome is a rigorous benchmark audit: NeurIPS D&B, Nature Methods correspondence/analysis, or Nature Communications at the very top. Tier 3 could reach 3–4 but depends on a composition-independent gene subset that transfers, which existing evidence (Nat Commun 2025) suggests is weak. Not Nature. |

---

## 6. Recommendation

**Do not run T1 as posed.** Run the inverted version, and run it as **the control experiment for
whichever thesis you do pick** rather than as a standalone paper. Its true value to MORPHEUS is that
it is cheap (42 GB, features already computed) and it converts the existing +0.07 bulk finding from
"a fact about our pipeline" into "a target-invariant law about the task" — which materially hardens
the confound-audit contribution that `A1_ledger/a_prior_benchmarks.md` already rates as the project's
most defensible novelty. Do F4 (the power check) first; it costs an afternoon and gates everything
else.

**Immediate next actions:**
1. Read HistoPrism (arXiv:2601.21560) in full — it is the F3 falsifier.
2. Bootstrap the per-task CI width on HEST-Bench — it is the F4 gate.
3. Re-run the Semantic Scholar sweep (unavailable this session) and verify TANGLE.
