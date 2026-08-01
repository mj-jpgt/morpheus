# T2 — Virtual perturbation: reading interventional response off histology

**Scout date:** 2026-07-29 · **Search tools:** WebFetch only (PubMed E-utilities, bioRxiv API;
Semantic Scholar and arXiv APIs were rate-limited/HTTP-429 for most of the session — noted where
this forced a COULD-NOT-VERIFY).

**Verification policy for this document:** every citation below was resolved to a PMID/DOI in this
session unless explicitly marked `COULD-NOT-VERIFY (this session)`. Nothing is cited from memory.

---

## 0. The claim under test

> Learn a perturbation-response space from genome-scale Perturb-seq (Replogle: 11,258 K562 KDs +
> RPE1/K562-essential), map TCGA patients (WSI ± bulk RNA) into it, and read out **per-patient
> predicted response to unseen genetic perturbations** — an *interventional* quantity — from an
> H&E image. Validate without wet lab via DepMap (lineage-matched), GDSC/CCLE (drug→target), and
> survival stratified by predicted dependency.

Decomposed, this is a chain of four inferences:

| # | Hop | Status in literature |
|---|-----|----------------------|
| H1 | H&E morphology → molecular state of *this* tumor | done many times; **effect size small within cancer type** |
| H2 | molecular state → position in a genetic-perturbation response space | done (GEARS/CPA/STATE/biolord lineage) |
| H3 | that position → response to an **unseen** perturbation | done in-vitro; **cross-context generalization is the open problem** |
| H4 | predicted response → a validatable clinical/functional quantity (DepMap, GDSC, survival) | done from **omics**, not from images |

The composition H1∘H2∘H3∘H4 — image-in, unseen-genetic-perturbation-response-out, at patient scale,
validated against DepMap — **I could not find published**. Every individual hop, and three of the
four adjacent pairs, **are** published. That is the honest novelty position.

---

## 1. Has this been done? — the five closest works

### 1.1 Coladan — *"image-only virtual perturbation"* — the sharpest collision
**Wang Z, Yang C, Tang X, Yin E, Yao Y, Luo Y, He J, Sun N.** "Trimodal, uncertainty-guided
whole-slide framework for genome-scale spatial expression and image-only virtual perturbation in
cancer cohorts." *Genome Medicine* 2026. **DOI 10.1186/s13073-026-01713-y** (PMID 42449400). ✅ verified

*What they did (from the verbatim abstract):* released **Coladan-human3K**, ~3,000 human spatial
transcriptomics profiles ("the largest human ST resource"); trained **Coladan**, a trimodal
(image · language · spatial-gene) whole-slide model predicting **genome-wide expression per spot with
calibrated uncertainty** while preserving foundation-model representations. Across **32 Visium
datasets** they improve spot Pearson **0.230 → 0.431 (~1.9×)**, show pathway-enrichment consistency,
and transfer **zero-shot to VisiumHD and spot-level Xenium**. Critically: *"Classification token (CLS)
embedding-only perturbation performs on par with expression-based baselines, enabling **image-only
virtual perturbation** without measured expression, illustrated on normal and cancer prostate
sections for in-situ hypothesis generation."*

*Why it matters to T2:* the **name and the headline framing are taken**, in a Genome Medicine paper,
this year. Anyone reviewing T2 will have this in hand.

*Where T2 is still distinct:* (i) their perturbation is an **embedding-space nudge** used for *in-situ
hypothesis generation*, not a prediction of a specific gene's knockdown transcriptional response;
(ii) no genome-scale **Perturb-seq** response space — no Replogle, no unseen-perturbation
generalization; (iii) **no external functional validation** — no DepMap, no GDSC, no survival;
(iv) illustrated on prostate sections, not a 32-cancer / 6,192-patient cohort; (v) spot-level, not
patient-level. So T2's differentiator must be *the validation apparatus and the unseen-perturbation
generalization*, **not** the phrase "virtual perturbation from images."

### 1.2 DeepDEP — omics → per-tumor cancer dependency map (the H2∘H4 precedent)
**Chiu YC, Zheng S, Wang LJ, Iskra BS, Rao MK, Houghton PJ, Huang Y, Chen Y.** "Predicting and
characterizing a cancer dependency map of tumors with deep learning." *Science Advances* 2021.
**DOI 10.1126/sciadv.abh1275** (PMID 34417181). ✅ verified

*What they did:* trained on CCLE/DepMap cell-line multi-omics (mutation, expression, CNV,
methylation) → CRISPR gene-effect, then **transferred the model to TCGA tumors** to produce a
predicted, per-patient dependency map, and characterized the resulting tumor dependencies.

*Why it matters:* the exact output T2 promises — *per-patient predicted genetic dependency for real
tumors, validated the DepMap way* — has existed for five years **from omics**. T2 gets no novelty
for the *output*; only for the *input modality* (H&E) and for the *perturbation-response-space*
formulation (predicting a full transcriptional response rather than a scalar gene-effect score).

### 1.3 Cell-line→tumor dependency transfer, modernized
"Deep Unsupervised Domain Adaptation for Translating Cancer Dependency Maps From Cell Lines to
Breast Cancer Tumor Genomics." *Genetic Epidemiology* 2026. **DOI 10.1002/gepi.70044**
(PMID 42339999). ✅ verified (title/venue/DOI; abstract not fetched)

And: **Bhattacharjee et al.**, "Gene dependency-informed inference of response to targeted cancer
therapies" (**FORGE**). *Nature Communications* 2026. **DOI 10.1038/s41467-026-73977-2**
(PMID 42259811). ✅ verified — jointly models drug response *and* gene essentiality from **basal gene
expression**, derives a "Benefit Score", reports erlotinib dependency concordance 0.69 / IC50
concordance 0.62, joint > single-task (p = 0.039), and validates on **PDX** and **Tahoe-100M**.

*Why they matter:* the cell-line→patient **domain-shift problem** for dependency is an active,
publishing subfield with strong baselines. T2's DepMap validation will be benchmarked against these,
and they use richer inputs (full omics) than T2's image.

### 1.4 Dawood et al. — drug sensitivity straight off H&E (the H1∘H4 precedent)
**Dawood M, Vu QD, Young LS, Branson K, Jones L, Rajpoot N, Minhas FUAA.** "Cancer drug sensitivity
prediction from routine histology images." *npj Precision Oncology* 2024.
**DOI 10.1038/s41698-023-00491-9** (PMID 38184744). ✅ verified

*What they did:* a proof-of-concept in **breast cancer** — patient-specific drug sensitivities were
**inferred from cell-line data via gene-expression mapping** (i.e. pseudo-labels, exactly T2's
label-generation trick), then a deep model was trained to predict them from H&E WSIs, across
multiple approved and experimental drugs.

*Why it matters:* T2's proposed pipeline is structurally **the same recipe**, swapping "drugs" for
"genetic perturbations" and "breast" for "32 cancers." The authors themselves frame it as
*proof-of-concept*, which is a fair signal for how far this recipe got.

Adjacent: **Liu H, Xie X, Wang B**, "Deep learning infers clinically relevant protein levels and drug
response in breast cancer from unannotated pathology images" (*npj Breast Cancer* 2024,
**DOI 10.1038/s41523-024-00620-y**, PMID 38413598) — WSI → 223 RPPA biomarkers on TCGA-BRCA +
CPTAC-BRCA, 0.79 AUC for trastuzumab response. ✅ verified

### 1.5 The morphology↔perturbation prior — and its ceiling
**Ramezani M, et al.** "A genome-wide atlas of human cell morphology." *Nature Methods* 2025.
**DOI 10.1038/s41592-024-02537-7** (PMID 39870862). ✅ verified — PERISCOPE: the first unbiased
**morphology-based genome-wide perturbation atlas**, CRISPR-Cas9 knockouts of **>20,000 genes** in
**>30 million cells**, three genome-wide genotype→phenotype maps; recovered TMEM251/LYSET biology.

**Way GP, Natoli T, Subramanian A, … Carpenter AE.** "Morphology and gene expression profiling
provide complementary information for mapping cell state." *Cell Systems* 2022.
**DOI 10.1016/j.cels.2022.10.001** (PMID 36395727). ✅ verified — A549 cells, **1,327 compounds × 6
doses**, L1000 vs Cell Painting head-to-head: the two assays capture *"both shared and complementary
information"*; Cell Painting profiles are **more reproducible and more diverse but measure fewer
distinct feature groups**.

*Why they matter — both directions:* PERISCOPE is the strongest existing evidence that **image
morphology genuinely encodes genetic-perturbation identity at genome scale**, which is T2's founding
premise. Way et al. is the caveat: morphology and expression are **partially disjoint** views, so a
morphology→expression-response map is *lossy by construction*, not merely noisy.

Also verified and directly relevant: **Wu B, et al.**, "Concordant transcriptional and morphological
remodeling revealed by in vivo Perturb-CLEAR", bioRxiv 2026, **DOI 10.64898/2026.04.06.716787** —
pooled CRISPR + whole-mount imaging **paired with Perturb-seq**, explicitly linking structural
phenotype to transcriptomic change *in vivo* (mouse cortex, NDD risk genes). This is the closest
existing "morphology ↔ Perturb-seq" joint measurement, and it is **not cancer, not H&E, not human**.

### 1.6 The evaluation warning shot (must be cited, cannot be dodged)
**Ahlmann-Eltze C, Huber W, Anders S.** "Deep-learning-based gene perturbation effect prediction does
not yet outperform simple linear baselines." *Nature Methods* 2025.
**DOI 10.1038/s41592-025-02772-6** (PMID 40759747). ✅ verified — *"We compared five foundation models
and two other deep learning models against deliberately simple baselines for predicting transcriptome
changes after single or double perturbations. **None outperformed the baselines.**"*

Combined with the on-disk lane-l12 entries (PerturBench; PertEval-scFM), the field's consensus is
that **perturbation-response prediction is a graveyard of unfalsified claims**. T2 inherits that
evidentiary bar *plus* an extra lossy modality hop.

### 1.7 What I searched and did **not** find (the open gap)
- PubMed: `(histology OR "whole slide" OR pathology) AND ("gene dependency" OR "CRISPR screen" OR
  "DepMap")` → 16 hits, **zero** relevant (all wet-lab papers where "pathology" means disease).
- PubMed: `("Perturb-seq" OR "perturbation atlas") AND (histology OR "whole slide" OR pathology OR
  "H&E" OR morphology)` → 15 hits; the only morphology↔Perturb-seq link is Perturb-CLEAR (mouse
  brain, above). **No H&E, no tumor, no patient.**
- PubMed: `("genome-scale" OR "genome-wide") AND "Perturb-seq"` → 41 hits; genome-scale contexts are
  **K562, RPE1, hiPSC/pluripotent, CD4+ T cells, mouse brain**. **No genome-scale Perturb-seq in a
  solid-tumor line exists** (see §3.3).

**Verdict on Q1:** *Not done as specified.* The nearest thing (Coladan) owns the phrase but not the
substance; the nearest substance (DeepDEP / FORGE / domain-adapted dependency transfer) owns the
output but not the modality; the nearest recipe (Dawood) owns the pipeline but for drugs, one cancer,
and self-describes as proof-of-concept. **The composition is open. The gap is real but narrow, and it
is narrow in the specific direction that makes it hard to defend as *important* rather than merely
*unattempted*.**

---

## 2. Is the core claim TRUE/plausible?

### 2.1 The arithmetic of the chain — this is the crux

Take the project's own established, real finding: **within-cancer, control-adjusted WSI→bulk-molecular
signal is ≈ +0.07 Pearson, and it is METHOD-INVARIANT** (every method, including baselines, lands in
the same place). Roughly 46–49% of the apparent cross-cancer performance is cohort/lineage structure.

Now propagate:

- **H1 contributes ≈ +0.07** of genuine patient-specific molecular information (within cancer type).
- **H2/H3**: the mapping from expression to perturbation response is itself not better than linear
  baselines (Nat Methods 2025), and cross-context transfer is the field's stated open problem.
- **H4**: DepMap gene-effect is dominated by **lineage** and by **expression of the target gene or a
  paralog** — meaning a "predicted dependency" that correlates with truth *across* cancer types is
  almost entirely a **re-encoding of predicted lineage**.

The composition of a +0.07 signal with two additional lossy, lineage-confounded hops does **not**
plausibly yield a per-patient interventional readout. **The claim as literally stated ("per-patient
predicted response to unseen genetic perturbations, from an image") is not plausible at a
clinically or biologically meaningful effect size.** The claim that *will* survive is a much weaker,
much more carefully specified one (§4).

What *is* plausible, and is supported by PERISCOPE + Perturb-CLEAR: **morphology carries
perturbation-relevant information in principle.** The failure mode is not "images are uninformative";
it is that **H&E of a fixed archival tumor section, at TCGA quality, aggregated to patient level, is a
vastly weaker morphological readout than Cell Painting of a live perturbed cell**, and the
patient-level variance available within a cancer type is small.

### 2.2 The specific T1 sub-question: do spot-level ST targets beat bulk targets for
morphology-predictability?

**Yes in reported absolute numbers, but the comparison as usually made is not valid, and I could not
find a like-for-like published head-to-head.**

Verified evidence:
- **Coladan (Genome Medicine 2026)**: spot-level genome-wide prediction across 32 Visium datasets,
  Pearson **0.230 (prior SOTA) → 0.431**. These are far above the ~+0.07 within-cancer
  control-adjusted bulk number.
- **"Benchmarking the translational potential of spatial gene expression prediction from histology"**,
  *Nature Communications* 2025, **DOI 10.1038/s41467-025-56618-y** (PMID 39934114) ✅ verified —
  eleven methods, five ST datasets, **external validation via TCGA**, evaluating "performance of
  predicted gene expression, model generalisability, translational potential, usability and
  computational efficiency." (The abstract does not expose the baseline-adjusted numbers; I did not
  retrieve full text — `COULD-NOT-VERIFY (this session)` for the specific Δ-over-baseline figures.)
- **SpaRED / SpaCKLE**, *Medical Image Analysis* 2025, **DOI 10.1016/j.media.2025.103754**
  (PMID 40885036) ✅ verified — 26 curated public ST datasets, 8 prediction models benchmarked; a
  transformer completion model cuts MSE by >82.5%, i.e. **a large fraction of the apparent difficulty
  is dropout/data quality, not biology**.
- **Hallinan C, Lucas CG, Fan J**, "Impact of Data Quality on Deep Learning Prediction of Spatial
  Transcriptomics from Histology Images", bioRxiv 2025, **DOI 10.1101/2025.09.04.674228**
  (PMID 40964396) ✅ verified — ablations show **sparsity and noise in the molecular data drive
  predictive performance**; imputation gives only marginal, non-generalizing gains; lower image
  resolution degrades accuracy *and* interpretability.

**The honest reading.** Spot-level numbers are higher, but they are a **different estimand**. Spot
Pearson is computed *across spots within a slide*, where the dominant variance is **cell-type
composition and tissue architecture** — stroma vs. tumor vs. immune vs. necrosis — which is exactly
what H&E trivially encodes. Bulk within-cancer control-adjusted Pearson removes the cohort mean and
asks a much harder question: *does this patient's tumor differ from the average tumor of this cancer
type in a way the image predicts?* Almost no spot-level ST paper reports the analogous
mean-spot-subtracted, cross-slide-held-out number.

**So: spot-level ST targets are materially easier, but the easiness is largely composition, not new
patient-specific biology.** If T2 (or T1) leans on ST to get a bigger number, a competent reviewer
will ask for the mean-baseline-subtracted, held-out-slide version — and the Nat Commun 2025 benchmark
and the Hallinan preprint are exactly the papers they will cite when asking.

*Caveat honestly stated:* I found **no published direct comparison** of "spot-level ST target vs.
bulk target, same slides, same encoder, same baseline adjustment." **COULD-NOT-VERIFY.** That absence
is itself a small, publishable methodological result — and it is a cleaner, safer thesis than T2.

### 2.3 The K562/RPE1 lineage gap — fatal or bridgeable?

**Verified state of the world:** genome-scale Perturb-seq exists only in **K562 (CML), RPE1 (retinal
epithelium, hTERT-immortalized, near-normal), hiPSC/pluripotent lines, primary CD4+ T cells, and
mouse brain** (PubMed sweep, 41 hits, §1.7). **There is no genome-scale solid-tumor Perturb-seq.**
Replogle et al., *Cell* 2022, **DOI 10.1016/j.cell.2022.05.013** (PMID 35688146) ✅ verified —
genome-scale CRISPRi across **>2.5 million human cells**, all expressed genes.

**Is it fatal for T2 as stated? Yes, for the strong version.** A K562 CRISPRi response manifold
encodes: erythroid/myeloid differentiation programs, BCR-ABL-driven signalling, and — the parts that
*do* transfer — **core-essential machinery** (ribosome biogenesis, proteostasis, mitochondrial
respiration, spliceosome, transcription). Replogle's own headline discoveries are exactly of this
core-essential type (CCDC86/ZNF236/SPATA5L1 ribosome biogenesis; C7orf26 transcription; TMEM242
mitochondrial respiration). Those are **pan-lineage housekeeping axes** — precisely the ones that are
*least* patient-discriminative in a tumor cohort, because every tumor needs them.

The perturbations that would make T2 clinically interesting (lineage-specific oncogene dependencies:
*EGFR* in lung, *ESR1* in breast, *AR* in prostate, *MYCN* in neuroblastoma) are the ones a
K562/RPE1 manifold has **no information about**. This is a structural, not statistical, gap.

**Bridges, ranked by how much they actually help:**

1. **Tahoe-100M (strongest bridge, verified).** Zhang J, Ubas AA, de Borja R, Svensson V, … Goodarzi H,
   Yu J. "Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent Gene Function
   and Cellular Modeling." bioRxiv, v1 2025-02-24 / v3 2025-05-10, **DOI 10.1101/2025.02.20.639398**
   ✅ verified — **100 million transcriptomic profiles**, **1,100 small-molecule perturbations**,
   **50 cancer cell lines**, Mosaic platform with reduced batch effects. This is *chemical*, not
   genetic, perturbation — but it spans **50 lineages**, which is the exact axis Replogle lacks. It is
   already being used as a validation set in the peer-reviewed literature (FORGE, Nat Commun 2026).
   **Recommendation: if T2 proceeds, Tahoe-100M should be the primary response space and Replogle the
   secondary, not the reverse.** That single swap changes T2's risk profile more than any modelling
   choice.
2. **LINCS L1000 (broad, shallow, verified).** Subramanian A, Narayan R, Corsello SM, Peck DD, et al.
   "A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles." *Cell* 2017,
   **DOI 10.1016/j.cell.2017.10.049** (PMID 29195078) ✅ verified — **1.3M profiles**, 978 measured
   landmark genes with inference for 81% of unmeasured transcripts, includes **shRNA/cDNA genetic**
   perturbagens as well as compounds, across many lineages. Available at **clue.io**. Advantage:
   genetic *and* chemical, multi-lineage, and directly comparable to the classic Connectivity-Map
   signature-reversal paradigm that already maps patient tumors into a perturbation space.
   Disadvantage: 978 landmarks, bulk, noisy, well-mined.
3. **Chemical→genetic translation is *not* free.** The on-disk lane-l12 entry Chem2Gen-Bench reports
   translation fidelity is "measurable but heterogeneous" — i.e. **drug ≠ clean knockdown of its
   target**. T2's GDSC validation via drug→target mapping inherits this error term and must report it,
   not assume it away. *(Chem2Gen-Bench itself: `COULD-NOT-VERIFY (this session)` — arXiv API was
   429-blocked; it is a prior-session on-disk entry, treat as unverified until re-checked.)*
4. **Perturb-CLEAR / spatial CRISPR screens** (verified, §1.5; and "Simultaneous CRISPR screening and
   spatial transcriptomics reveal intracellular, intercellular, and functional transcriptional
   circuits", *Cell* 2025, **DOI 10.1016/j.cell.2025.02.012**, PMID 40081369 ✅ verified) — the future
   right answer, but neither is human tumor H&E at cohort scale today.

---

## 3. Open-access data: what exists, real sizes, real limits

### 3.1 Already on disk (no acquisition risk)
| Asset | Scale | Limitation for T2 |
|---|---|---|
| TCGA WSI + bulk RNA | 6,192 patients, 32 cancers, H-Optimus features + BulkFormer, held-out-cancer split | archival FFPE H&E; one/few slides per patient; **patient-level variance within cancer type is the scarce resource** |
| Replogle GWPS | K562 11,258×8,248; K562-essential 2,285; RPE1 2,679 | **two lineages, neither a solid tumor**; CRISPRi (soft) not KO |
| CCLE (436 MB) + GDSC (131 MB) | ~1,000 lines × ~500 compounds (GDSC1+2) | cell lines ≠ tumors; IC50 batch/assay effects |
| CPTAC proteomics/phospho | inventoried, unwired | small n per cancer; would be the strongest *orthogonal* validation if wired |

### 3.2 Open, needed, and obtainable
| Dataset | Access route | Real scale | Real limitation |
|---|---|---|---|
| **DepMap CRISPR gene effect (Chronos)** | depmap.org quarterly public release (`CRISPRGeneEffect.csv`), scripts already on disk | ~1,100+ lines × ~18,000 genes *(exact release counts `COULD-NOT-VERIFY (this session)`)* | dependency is **lineage-dominated**; most selective dependencies are predicted by target/paralog expression — a lineage-confounded validation target |
| **Tahoe-100M** | bioRxiv 10.1101/2025.02.20.639398; released publicly (HuggingFace `tahoebio/Tahoe-100M` — *release channel not re-verified this session*) | 100M cells, 1,100 compounds, **50 cancer cell lines** | chemical not genetic; cell lines; enormous (needs streaming/subsampling on a single A100-40GB) |
| **LINCS L1000** | clue.io; GEO series *(accessions not re-verified this session)* | 1.3M profiles ✅ | 978 landmarks; bulk; heavily mined |
| **Coladan-human3K** | with Genome Medicine 2026 paper (10.1186/s13073-026-01713-y) | ~3,000 human ST profiles ✅ | brand-new; the authors are direct competitors on this exact framing |
| **SpaRED** | with Med Image Anal 2025 (10.1016/j.media.2025.103754) | 26 curated public ST datasets, 8 benchmarked models ✅ | Visium dropout; standardization is the point, scale is modest |
| **Replogle GWPS raw** | authors' public release *(portal/accession `COULD-NOT-VERIFY (this session)`)* | already on disk | — |

### 3.3 What does **not** exist and cannot be bought with an open download
- **Genome-scale Perturb-seq in a solid-tumor line.** Verified absent (§1.7). No open-access route.
- **Paired H&E + genetic-perturbation response for the same tumor.** Does not exist for human tumors.
  Perturb-CLEAR is the closest and is mouse brain.
- **Any ground truth for "this patient's tumor's response to knocking down gene X."** This is the
  quantity T2 claims to predict, and **there is no dataset in which it is measured.** Every proposed
  validation (DepMap lineage-matched, GDSC drug→target, survival) is a **proxy**, and every proxy is
  lineage-confounded. This is not a data-acquisition problem; it is the thesis's epistemic ceiling.

---

## 4. The highest-ceiling **honest** claim, and what falsifies it

### 4.1 What must be abandoned
- ❌ "Per-patient interventional readout from an image." Not supportable at +0.07 within-cancer signal.
- ❌ "Virtual perturbation from histology" as a novelty phrase. Taken (Coladan, Genome Medicine 2026).
- ❌ "Predicted dependency stratifies survival" as evidence of interventional content. Survival is
  predicted by grade/stage/subtype, all visible in H&E, all lineage-linked. This validation **cannot
  fail** and therefore proves nothing.

### 4.2 The claim that survives
> **An image-only estimator can place a tumor on a genome-scale perturbation-response manifold, and
> the fraction of that placement which is patient-specific rather than lineage-explained is
> measurable, small, and concentrated in a nameable class of perturbations.**
>
> Concretely: build WSI → perturbation-response-space mapping; then **decompose** predicted
> dependency into (a) cancer-type-mean, (b) predicted-lineage-explained, (c) patient-residual. Report
> the residual honestly, per perturbation. The deliverable is *the decomposition and the benchmark*,
> not the predictor.

This reframes T2 from *"we built a virtual perturbation oracle"* (unsupportable, and pre-empted) to
*"we established how much interventional information is actually in an H&E slide, and separated it
from the lineage shortcut that makes every prior result look better than it is."* That is the same
intellectual move as Ahlmann-Eltze/Huber/Anders (Nat Methods 2025) and PertEval-scFM, applied to a
modality where nobody has done it — and it is **methodologically bulletproof precisely because it
survives a null result.** A rigorous, well-powered null here is publishable; a weak positive is not.

Secondary, and the only place a *positive* claim is defensible: **core-essential / proteostasis /
ribosome-biogenesis axes are pan-lineage** (Replogle's own discoveries are of this type), so if H&E
predicts anything interventional, it is most likely proliferation-coupled essentiality. That is a
narrow, testable, biologically coherent hypothesis with a real mechanism (mitotic figures, nuclear
morphology, and cellularity are visible in H&E and are the morphological correlates of ribosome-biogenesis
load). **It is also the least clinically actionable class of dependency** — an honest paper must say so.

### 4.3 Falsification tests (pre-register all four; run #1 first, it is cheap and decisive)

1. **The lineage-ablation test (the kill switch).** Within each cancer type, rank per-patient
   predicted dependencies and compare Spearman-to-DepMap-lineage-mean against a **cancer-type-mean
   predictor that ignores the image entirely**. If Δ ≈ 0 — as the +0.07 method-invariant finding
   predicts — **T2 is dead** and no architecture rescues it. *Run this in week 1, on existing on-disk
   assets, before writing a line of model code.*
2. **The permutation test.** Shuffle patient↔image within cancer type. If performance barely drops,
   the model is reading lineage, not the patient.
3. **The unseen-perturbation test.** Hold out entire perturbation classes (not random perturbations)
   and require beating (i) the mean-response baseline and (ii) a linear baseline, per
   Ahlmann-Eltze et al. If it does not, the "unseen perturbation" claim is void.
4. **The modality-necessity test.** WSI-only vs. RNA-only vs. WSI+RNA. Given the project's own
   established finding, RNA-only will likely dominate and WSI will add ~nothing. If so, the paper's
   honest title contains the words "does not."

### 4.4 Recommended de-risking, in order
1. **Run falsification test #1 before anything else.** It costs a day and settles the thesis.
2. **Swap the response space: Tahoe-100M (50 cancer lineages) primary, Replogle secondary.** This is
   the single highest-value change and it directly addresses the stated risk.
3. **Drop survival as a validation.** It cannot fail; it therefore carries no evidence.
4. **Make the decomposition the deliverable.** Reframe as an evaluation-paradigm paper (the project's
   demonstrated strength: leakage control, control-adjustment, method-invariance detection). That
   plays to the established +0.07 / effective-rank findings rather than fighting them.
5. **If a positive result appears, immediately check it is not Coladan's result restated** at patient
   rather than spot level.

---

## 5. Ratings

| Axis | Score | Justification |
|---|---|---|
| **Novelty** | **3 / 5** | The *composition* (genome-scale Perturb-seq response space + WSI → unseen-perturbation response + DepMap/GDSC validation, 32 cancers) is unpublished — verified by three targeted PubMed sweeps returning zero. But the *phrase and framing* are taken (Coladan, Genome Med 2026), the *output* is taken (DeepDEP 2021, and 2026 domain-adapted successors), the *recipe* is taken (Dawood 2024), and the *premise* is taken (PERISCOPE 2025). Novelty is real but compositional and thin — a reviewer can name four papers that each cover a face of it. |
| **Feasibility (open data)** | **3 / 5** | Data and compute are fully in hand — TCGA/WSI/RNA, Replogle, CCLE/GDSC on disk; DepMap, Tahoe-100M, LINCS all openly obtainable; A100-40GB suffices with streaming. Engineering risk ≈ 0. But feasibility of reaching a *defensible* result is middling: the chain has three independently-documented-weak hops, and **no dataset anywhere measures the ground truth being predicted**. Executable: 5. Executable-to-a-positive-result: 2. Blended honest: 3. |
| **Ceiling** | **3 / 5** | The strong version, if true, is Nature Medicine / Nature Cancer / Nature BME tier — but it is not true at the effect size the project's own measurements imply. The honest version (lineage-vs-patient decomposition of image-derived interventional predictions; first calibrated estimate of interventional information in H&E) is a strong *Nature Methods / Nature Communications* paper and a genuinely useful field correction. It is not Nature, because the headline is a negative/bounding result. Ceiling 3, with 4 reachable **only** if the Tahoe-100M swap unlocks a lineage-specific dependency signal that survives falsification test #1. |

---

## 6. Verification ledger

**Verified this session (PMID and/or DOI resolved via PubMed E-utilities or bioRxiv API):**
Coladan (PMID 42449400, 10.1186/s13073-026-01713-y) · DeepDEP (PMID 34417181, 10.1126/sciadv.abh1275) ·
Dawood npj Precis Oncol (PMID 38184744, 10.1038/s41698-023-00491-9) · wsi2rppa (PMID 38413598,
10.1038/s41523-024-00620-y) · DUDA cell-line→breast-tumor dependency (PMID 42339999, 10.1002/gepi.70044) ·
FORGE (PMID 42259811, 10.1038/s41467-026-73977-2) · Ahlmann-Eltze/Huber/Anders (PMID 40759747,
10.1038/s41592-025-02772-6) · PERISCOPE morphology atlas (PMID 39870862, 10.1038/s41592-024-02537-7) ·
Way et al. Cell Systems (PMID 36395727, 10.1016/j.cels.2022.10.001) · Perturb-CLEAR (PMID 41993505,
10.64898/2026.04.06.716787) · CRISPR+spatial Cell 2025 (PMID 40081369, 10.1016/j.cell.2025.02.012) ·
Replogle GWPS (PMID 35688146, 10.1016/j.cell.2022.05.013) · Tahoe-100M (bioRxiv
10.1101/2025.02.20.639398) · LINCS L1000 (PMID 29195078, 10.1016/j.cell.2017.10.049) ·
Nat Commun ST benchmark (PMID 39934114, 10.1038/s41467-025-56618-y) · SpaRED/SpaCKLE (PMID 40885036,
10.1016/j.media.2025.103754) · Hallinan data-quality preprint (PMID 40964396, 10.1101/2025.09.04.674228).

**COULD-NOT-VERIFY (this session):** exact DepMap public-release line/gene counts; Replogle raw-data
portal/accession; LINCS GEO series accessions; Tahoe-100M HuggingFace release channel; HEST-1k
(arXiv API 429 throughout); Chem2Gen-Bench (on-disk prior-session entry, arXiv API 429);
baseline-adjusted Δ figures inside the Nat Commun 2025 ST benchmark (full text not retrieved).
**Semantic Scholar and arXiv APIs returned HTTP 429 for essentially the whole session; PubMed and
bioRxiv carried the sweep.** No citation in this document was written from memory.
