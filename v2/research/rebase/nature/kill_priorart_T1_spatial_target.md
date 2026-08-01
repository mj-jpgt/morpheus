# Prior-Art Kill Attempt — T1_spatial_target (inverted: "target-invariant ceiling")

**Date:** 2026-07-29
**Assassin verdict:** **NOT KILLED** on novelty — but the *framing* is partially pre-occupied, and the venue ceiling is lower than the thesis assumes.
**Search constraint honoured:** WebSearch not used (exhausted). All hits below retrieved via arXiv API, Europe PMC REST, Semantic Scholar graph API, OpenAlex, and direct publisher fetch. Every citation was fetched in-session; anything I could not fetch is explicitly marked COULD-NOT-VERIFY.

---

## 1. What would have killed it

A published paper that does any of:

- **K1.** Reports a **mean-expression / constant-predictor baseline** on HEST-Bench (or an equivalent spot-level ST-from-H&E benchmark) and shows the deep/FM models' margin over it is small.
- **K2.** Performs a **per-slide-mean-subtracted (centered) re-evaluation** of spot-level ST prediction and reports the collapse in Pearson r.
- **K3.** Shows a **cell-composition-only baseline** (nuclei counts/types, e.g. CellViT/HoVer-Net) matches pathology-foundation-model features for spot expression prediction.
- **K4.** Uses a **random-gene-panel null** to adjust ST-prediction leaderboards.
- **K5.** Explicitly compares the **bulk-RNA-from-WSI ceiling to the spatial-ST-from-H&E ceiling** and argues they are the same quantity ("target invariance").

**None of K1–K5 was found in published form.** K5 in particular returned zero hits across every API.

---

## 2. Closest prior art (verified)

### 2.1 HESCAPE — the single biggest novelty threat
**"A Large-Scale Benchmark of Cross-Modal Learning for Histology and Gene Expression in Spatial Transcriptomics"**, Gindra, Palla, Nguyen, Wagner, Tran, Theis, Saur, Crawford, Peng. **arXiv:2508.01490v2** (2025-08-02, v2 2025-08-27). Abstract retrieved verbatim.

Verbatim, load-bearing sentences:
> "...downstream task evaluation reveals a striking contradiction: while contrastive pretraining consistently improves gene mutation classification performance, it degrades direct gene expression prediction compared to baseline encoders trained without cross-modal objectives. **We identify batch effects as a key factor that interferes with effective cross-modal alignment.**"

Why it matters: "batch effect" in a 6-panel / 54-donor pan-organ ST benchmark **is** the per-slide/per-donor mean. A reviewer will say HESCAPE already named the confound MORPHEUS proposes to subtract. What HESCAPE does *not* do, from the abstract: no constant/mean-expression baseline, no centered (mean-removed) leaderboard, no composition-only baseline, no random-gene null, no bulk comparison. It diagnoses batch effect as an obstacle to *alignment*, not as the thing the *leaderboard metric is actually measuring*. So it erodes the "the field has never noticed" rhetoric without doing the decomposition.

### 2.2 CHRep — names the exact failure mode, as a method paper
**"CHRep: Cross-modal Histology Representation and Post-hoc Calibration for Spatial Gene Expression Prediction"**, **arXiv:2604.21573v1** (2026-04-23).
States that under realistic **leave-one-slide-out** evaluation, existing models struggle with **"slide-level appearance shifts and regression-driven over-smoothing."** "Regression-driven over-smoothing" is precisely "the model predicts the slide mean." CHRep treats this as a problem to *fix* (post-hoc calibration), never as a measurement artifact to *quantify*. No mean baseline, no decomposition.

### 2.3 COAST and HEXST — implicit concessions that absolute expression is mean-dominated
- **COAST**, **arXiv:2607.09166v1** (2026-07-10), abstract verbatim: *"Existing context-aware methods mainly supervise absolute expression, while relative expression relationships between spots are rarely used explicitly... trained with a joint objective that combines absolute expression regression with **signed differential regression between the target and context spots**."*
- **HEXST**, **arXiv:2605.04682v1** (2026-05-06): existing regression methods produce **"over-smoothed gene expression profiles"**; adds a **"contrast-sensitive objective."**
- **HiST** (2026-06-12, arXiv id COULD-NOT-VERIFY): adds a **"slide calibration token"** for acquisition variation.

These three are the field independently converging on "the within-slide, mean-removed component is the hard part." They do it as *training tricks*, not as an *evaluation correction*, and none reports what the metric looks like after the mean is removed. This is the strongest evidence that MORPHEUS's Tier 1 number is (a) not published and (b) something the field is circling.

### 2.4 HistoPrism — the "the metric measures the wrong thing" move, already made once
**"HistoPrism: Unlocking Functional Pathway Analysis from Pan-Cancer Histology via Gene Expression Prediction"**, **arXiv:2601.21560v3** (2026-01-29).
Argues prior work is limited to "per-cancer settings and variance-based evaluation" and shifts assessment "from isolated gene-level variance to coherent functional pathways." Same *genre* of contribution (re-specify the ST-prediction metric), different axis (pathway vs. gene, not composition vs. residual). Reduces the novelty of "we ship a new evaluation paradigm for ST prediction."

### 2.5 sCellST — closest thing to Tier 2, but framed positively
**"sCellST predicts single-cell gene expression from H&E images"**, Nature Communications 2026, **DOI 10.1038/s41467-025-67965-1**. Abstract verbatim: a segmented-cell MIL model *"predicts single-cell gene expression from morphology, **matching patch-based methods on spot level prediction tasks**."*

This is a nuclei-derived model matching patch-FM models at spot level — i.e. the empirical content of Tier 2 is arguably already on the record. But: (i) it is sCellST's *positive* selling point, not a debunk; (ii) it is cell-*morphology* MIL, not composition-only counts/types; (iii) it is two cancer datasets, not HEST-Bench; (iv) no FM-minus-composition delta is reported. Tier 2 survives, but the parent should expect "sCellST already showed this" as a reviewer reflex, and must run composition-*only* (counts/types, no morphology embedding) to differentiate.

### 2.6 Composition-from-H&E is a solved, well-populated channel
- **Hist2Cell**, Cell Genomics 2026, **DOI 10.1016/j.xgen.2025.101137** — fine-grained cellular architecture from histology.
- **"Integrating Pathology Foundation Models and Spatial Transcriptomics for Cellular Decomposition from Histology Images"**, arXiv 2025-07-09 (arXiv id COULD-NOT-VERIFY) — a *lightweight regressor* on frozen FM features predicts cell-type composition accurately.
- **VirtualST**, **DOI 10.18063/cbr.v7i2.1919** (2026) — conditions ST prediction jointly on FM features **and** "local cell-type composition derived from nuclei segmentation," i.e. already treats composition as a distinct, additive input channel.
- **SHEST**, bioRxiv **10.1101/2025.11.19.689364**; **STHELAR**, Sci Data 2026 **10.1038/s41597-026-06937-6**.

Implication: the *premise* of Tier 2 (H&E → composition is easy) is thoroughly established. The *unmeasured* quantity is the FM-over-composition **delta** on HEST-Bench. That delta is the novel object, and nobody has reported it.

### 2.7 Negative evidence — the strongest asset for the thesis
**"A comprehensive survey of computer vision methods for spatial transcriptomics"**, Briefings in Bioinformatics 2026, **DOI 10.1093/bib/bbag255** (fetched via OUP). Catalogues 46+ H&E→ST models, 2020–2025. Fetched and interrogated directly on all four controls:
- mean-expression / constant-predictor baselines: **not mentioned**
- slide-level mean / batch confound inflating Pearson r: **not addressed**
- cell-type composition as the primary performance driver: **not discussed**
- random-gene-panel nulls: **not mentioned**

The only adjacent caution it raises is about adjacent-section leakage:
> "because ST data often come from adjacent layers or regions of the same sample, these sections usually share similar tissue morphology and molecular patterns, which can inflate reported performance"

A 2026 survey of the whole subfield missing all four controls is close to proof that the control-adjusted leaderboard does not exist.

### 2.8 The leaderboard is still being chased on raw r
Verified 2026 HEST-Bench SOTA claims with no control adjustment anywhere: **MINT** (2026-03-09, "best overall performance on HEST-Bench, mean Pearson r = 0.440"; arXiv id COULD-NOT-VERIFY), **SEAL** (2026-02-15), **MoLF** (2026-02-02), **TICON** (2025-12-24), **STPath** (npj Digit Med 2025, **10.1038/s41746-025-02020-3**). The r≈0.42–0.44 regime the thesis targets is confirmed as the live, uncorrected number.

### 2.9 Other checked-and-cleared
- **HEST-1k**, Jaume et al., **arXiv:2406.16192** — abstract verified (1,229 profiles, 153 cohorts, 26 organs, 76M nuclei). The nuclei needed for Tier 2 are confirmed to ship with the dataset. The full baseline list could not be extracted from the abstract page — **COULD-NOT-VERIFY** the thesis's specific claim that HEST-Bench ships *only* Ridge+PCA and RF, though nothing contradicts it.
- **SpaRED/SpaCKLE** (arXiv 2025-05-05, id COULD-NOT-VERIFY) — curated benchmark fixing dataset/preprocessing inconsistency; a *fair-comparison* fix, not a confound decomposition.
- **MV_Hybrid** (arXiv 2025-08-01, id COULD-NOT-VERIFY) — LOSO robustness across backbones; no mean baseline.
- **"Generalization of deep learning models for predicting spatial gene expression profiles using histology images"**, bioRxiv **10.1101/2023.09.20.558624** — full text returned HTTP 403; **COULD-NOT-VERIFY** whether it contains a mean baseline. Title/scope suggest generalization, not decomposition. *Residual risk: low but non-zero. Recheck before writing.*
- **SEQUOIA**, Nat Commun 2024 **10.1038/s41467-024-54182-5** (bulk from WSI) — no cross-cancer-structure confound analysis found.
- Zero hits on OpenAlex/Europe PMC/arXiv for any bulk-vs-spatial ceiling comparison (K5).

---

## 3. Verdict and the honest damage report

**killed = false.** No paper does the control-adjusted HEST-Bench leaderboard, the composition-only delta, the random-gene null, or the bulk/spatial target-invariance claim.

But three things about the thesis do not survive contact with the literature:

1. **"The field demonstrably lacks the diagnostic" is too strong.** HESCAPE (2.1) named batch effects as the confound in a large ST benchmark in Aug 2025; CHRep named "slide-level appearance shifts and regression-driven over-smoothing" in Apr 2026; COAST/HEXST/HiST all ship training-time fixes for exactly this. The correct claim is narrower and still true: *nobody has quantified it as an evaluation correction.* Write it that way or a reviewer will kill it for you.

2. **Tier 2 is partially pre-empted by sCellST.** A nuclei-based model already "matches patch-based methods on spot level prediction" in Nature Communications. To keep Tier 2 damning, the baseline must be composition-**only** — CellViT counts/types per spot, no morphology embedding, no learned cell representation — and the deliverable must be the FM-minus-composition delta on HEST-Bench, which is genuinely unreported.

3. **Venue ceiling.** The literature shape says this is a benchmark-audit paper. The comparable prior work (HESCAPE, SpaRED, HistoPrism) sits at MICCAI/NeurIPS-D&B/Brief Bioinform, not Nature. Tier 1+2 as scoped is a strong NeurIPS Datasets & Benchmarks paper or a *Nature Methods* Brief Communication / Matters Arising at best. Only Tier 3 (composition-independent gene set that transfers to held-out cancers *and* to bulk TCGA) reaches Nature-tier, and the thesis itself prices Tier 3 as unlikely.

**Recommendation to the parent:** consistent with the thesis's own closing line — do not run this as the headline thesis. Run Tier 1 + the composition-only Tier 2 as the ~2–4 week control experiment that hardens whichever thesis is chosen, and cite HESCAPE/CHRep/COAST up front as independent corroboration that the confound is real rather than pretending the field is blind to it. Gate on the per-task bootstrap-CI power check first.

**Recheck before writing:** bioRxiv 10.1101/2023.09.20.558624 full text (403 this session).
