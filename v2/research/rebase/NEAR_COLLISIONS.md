# NEAR_COLLISIONS — Adjudicated by reading the actual papers

**Date:** 2026-07-29
**Method:** Direct WebFetch on arXiv abs/HTML, ar5iv full text, arXiv API, OpenReview API, Semantic Scholar API.
WebSearch was **not** used (budget exhausted). Where a paper could not be located, this is stated explicitly
rather than guessed.

> **Headline warning before anything else.** Two of the five "threats" carried in
> `NOVELTY_LEDGER.md` — **"Nguyen et al., *Answer, Clarify, or Abstain*"** and **"Decode-gLM"** — **could not be
> found to exist** in any index (arXiv, OpenReview, Semantic Scholar). A third, **"VCBench (single-cell)"**, also
> could not be found; the only `VCBench` on arXiv is a venture-capital benchmark. These appear to be artefacts of
> the earlier sweeps. Real, *verified* nearest neighbours are substituted below and are in two cases **closer
> threats than the phantom citations were**.

---

## (a) Auto-detect the latent BIOLOGICAL task; multi-interpretation; abstain when no pathway slot addresses it

### 1. What the paper actually does

**The named threat — "Nguyen et al., *Answer, Clarify, or Abstain: Fine-Grained Selective Prediction for Medical
VLMs*, ICML 2026 workshop" — could not be located.** Negative results from four independent indexes:

| Query | Index | Result |
|---|---|---|
| `all:"Answer, Clarify, or Abstain"` | arXiv API | 1 hit — **Baan et al.**, not Nguyen |
| `ti:"Fine-Grained Selective Prediction"` | arXiv API | **0 results** |
| `au:Nguyen AND abs:"abstain"` | arXiv API | 7 hits, **none** med-VLM triage |
| `"Answer Clarify or Abstain Medical VLMs"` | Semantic Scholar | **no match** |
| `"Answer Clarify or Abstain ... Selective Prediction Medical"` | OpenReview API | **no match** (returns unrelated fine-grained/selective papers) |

A genuinely non-archival workshop paper can be invisible to all five, so I cannot prove non-existence — but I can
say the citation as written is **unverifiable**, and every specific attribute the ledger asserts about it
("MoE+RankNet", "frozen-hidden-state latent-mode auto-detector", "input-quality condition", "radiology/CT/fundus")
is currently **unsourced**. Do not cite it.

**Verified real nearest neighbours** (these are the actual prior art):

- **Baan, Aziz, Plank, Fernández — "Clarify, Abstain or Answer? Strategising in Conversation with
  Belief-Augmented Generation" (arXiv:2605.25831, May 2026).** *This is the real three-way-triage paper.*
  Belief-Augmented Generation (BAG): sample multiple outputs from an LLM, prompt the model to reason over its own
  sample set, then choose answer / clarify / decline. Six models, multi-turn **ambiguous** QA. Reports improved QA
  accuracy and strategy choices better aligned to actual model uncertainty; explicitly notes that
  **distinguishing clarification from abstention remains difficult**. Text-only, no vision, no biology.
- **Khanmohammadi, Thind, Ghassemi — "Calibrated Triage, Not Autonomy: Confidence Estimation for Medical
  Vision-Language Models" (arXiv:2606.15910, Jun 2026).** Five open-weight LVLMs, three medical VQA datasets
  (broad clinical imaging, radiology, pathology). Latent = **answer trustworthiness / calibration**, motivated by
  "a VLM can answer a question about a medical image fluently and confidently while barely using the image."
  **Two-way** automate-vs-defer, not three-way. Best estimators cut high-confidence errors from 41–45% to 1–4%.
  Confidence estimators trained on natural images transfer without adaptation (a cross-domain, not
  cross-task-family, generalisation test).
- **Presacan et al. — "Ask Before You Diagnose: Safe-Psych" (arXiv:2607.13036, 2026).** A sequential benchmark
  with literal **DIAGNOSE / CLARIFY / ABSTAIN** actions for LLMs in psychiatry under incomplete clinical
  information. Text-only, no imaging, no molecular data.

### 2. The precise overlap

The **answer/clarify/abstain triage primitive is thoroughly published** — independently, at least three times in
2026 alone, across conversational LLMs (Baan), medical VLMs (Khanmohammadi), and clinical benchmarks (Presacan).
Baan additionally overlaps on *ambiguity as the trigger*: the choice is driven by the model's latent uncertainty
over an under-specified query, which is structurally the same control flow MORPHEUS proposes. Khanmohammadi
overlaps on *medical multimodal* + selective prediction. Presacan overlaps on the *exact three-way label set* in a
clinical setting.

**The ledger's stated differentiator was right for the wrong reason.** It claimed the differentiator holds because
Nguyen's latent is "input-quality condition, not biological task." That reasoning is sound and *does* hold for
Khanmohammadi (latent = calibration/answerability) — but Baan's latent is *semantic ambiguity of the query*, which
is a closer relative and was not on the ledger at all.

### 3. The defensible remaining delta

**DELTA SURVIVES — but the triage mechanism itself must be dropped as a novelty claim.** Frame it as adopted prior
art (cite Baan, Khanmohammadi, Presacan), and claim only:

- The latent variable is the **biological task/programme**, inferred over a **frozen WSI + molecular tumour
  representation** — no verified prior work detects a *biological* latent task; all three detect uncertainty,
  calibration, or information sufficiency.
- **Abstention is grounded in representation coverage, not confidence**: "no identified pathway slot addresses
  this question" is a *structural* abstention criterion tied to the slot inventory. Every paper found abstains on
  a *scalar uncertainty threshold*. This is the sharpest and most defensible difference.
- **Returning multiple valid interpretations** when under-specified (SAM-style). Baan's model *asks* for
  clarification; it does not *enumerate* the competing readings. Note Baan explicitly reports the
  clarify-vs-abstain boundary is hard — cite this as motivation, and beat it.
- **Held-out task-family** evaluation. Khanmohammadi's cross-domain transfer is the nearest analogue and is
  weaker (same task, new imaging domain).

⚠️ Do not claim novelty on "abstention in a medical multimodal model." That is now firmly occupied.

### 4. Confidence

**COULD-NOT-ACCESS** for the named Nguyen threat (unverifiable; likely a fabricated citation).
**READ-ABSTRACT-ONLY** for Baan, Khanmohammadi, Presacan (abstracts + metadata verified via arXiv API/abs pages).

---

## (b) IDENTIFIABILITY (not interpretability) buys prompt reliability and cross-cohort transfer

**Threat: SurvPath — Jaume, Vaidya, Chen, Williamson, Liang, Mahmood; CVPR 2024; arXiv:2304.06819.**

### 1. What the paper actually does

Multimodal survival prediction fusing whole-slide images with bulk transcriptomics. Two stated contributions:
(i) tokenize transcriptomics "in a semantically meaningful and interpretable way"; (ii) capture dense multimodal
interactions between the two modalities.

**Pathway tokenization (the load-bearing detail).** Verified from ar5iv full text:

- Pathways come from **Reactome (1,281 available)** and **MSigDB Hallmarks (50 available)**.
- After filtering for ≥90% transcriptomic coverage they retain **331 pathways over 4,999 genes** (281 Reactome +
  50 Hallmarks).
- Each pathway token is produced by a **sparse MLP**, `x_i^(P) = φ_i(g_p_i)` — a *per-pathway* learned encoder
  whose input support is **fixed by the curated gene set membership**.
- **Therefore: the pathway inventory is fixed a priori by curated gene sets; only the within-pathway MLP weights
  are learned.** The abstract's phrase "learn biological pathway tokens" refers to learning the *embedding*, not
  discovering the *pathway set*. This confirms the repo's existing note in `applicability/axis2.md`.

**Evaluation.** Five TCGA cohorts: BLCA (n=359), BRCA (n=869), STAD (n=317), COADREAD (n=296), HNSC (n=392).
**No external / non-TCGA validation cohort.** Claims SOTA c-index plus an "interpretability framework" for
identifying prognostic factors.

**The word "identifiability" / "identifiable" does not appear anywhere in the paper.**

### 2. The precise overlap

Overlap is **architectural and nomenclatural only**. SurvPath establishes that (i) a multimodal tumour model can
carry units *named after biological pathways*, (ii) those units can be co-attended with histology, and (iii) they
can be read out post hoc for interpretation. If MORPHEUS's pitch is "we have pathway-named multimodal tokens," it
is pre-empted outright. That is the whole of the collision.

### 3. The defensible remaining delta

**DELTA SURVIVES — cleanly and with a wide margin.** SurvPath does not touch any of the four things the claim
actually asserts:

- **Not identified.** Fixed curated gene sets with learned within-set weights are the textbook *non*-identified
  case. No latent-variable identification, no interventional/perturbation signal, no sparse-mechanism-shift, no
  iVAE-style conditioning. There is no theory in the paper at all.
- **Never claims identifiability.** Zero occurrences. The paper's own framing is *interpretability*.
- **Never tests prompt reliability.** SurvPath has **no prompting interface whatsoever** — it is a fixed-head
  survival regressor. There is no query, so there is nothing to be reliable about.
- **Never tests transfer.** Five TCGA cohorts, within-cohort cross-validation, **no external cohort, no
  leave-site-out, no leave-cancer-out.** The transfer axis is entirely untouched.

So the specific proposition — *identifiability, as distinct from interpretability, measurably improves prompt
reliability and cross-cohort transfer* — is not merely unproven by SurvPath; SurvPath does not instrument either
dependent variable. It is the ideal **non-identified baseline arm** for the experiment, not a scoop.

⚠️ Two cautions. (1) The experiment must be run *against* SurvPath's 331-token construction as the control arm,
or the claim is untested rhetoric. (2) SurvPath's absence of external validation means you cannot cite its numbers
as the transfer baseline — you must re-run it under your own leave-site-out / leave-cancer-out protocol.

### 4. Confidence

**READ-FULL-TEXT** (ar5iv full text for the tokenization and dataset sections; arXiv abs for the abstract).

---

## (c) Frozen multimodal trunk + pathway slots with a BLOCK-IDENTIFIABILITY guarantee, NL-addressable

**Threat: Winter, Vonficht, Le Bescond, Gebbe, Rosati, Chen, Schick, Stewart, Brieu — "Data-Efficient Multimodal
Alignment for Histopathology-based Molecular Prediction," arXiv:2606.29949, submitted 29 Jun 2026.**

### 1. What the paper actually does

Trains a **lightweight alignment module on top of frozen histopathology and RNA-seq foundation models** to enable
"open-vocabulary molecular prompting" — querying H&E slides with gene-set signatures to predict pathway activity
without sequencing or end-to-end retraining.

- **Encoders (all frozen).** Histology: **THREADS** image encoder over CONCH ViT features aggregated by attention
  MIL (d_v = 1024). RNA: two variants tested — scGPT + THREADS RNA-head (d_r = 1024), and **BulkFormer**
  (d_r = 512, purpose-built for bulk RNA-seq). Best configuration: **BulkFormer MLP-CLIP with soft-kNN querying**.
- **Training.** Contrastive (CLIP-style) alignment on a multi-cancer cohort, **N = 1,720**.
- **Prompting mechanism — important.** Queries are **NOT free natural language.** Gene sets act as
  "open-vocabulary molecular prompts" by computing **ssGSEA scores against a reference RNA library**, then
  performing **soft-kNN retrieval in the aligned latent space**. Open-vocabulary in the *gene-set* sense, not the
  *natural-language* sense.
- **Headline result.** 25-fold improvement in retrieval over baselines.
- **Graduated predictability spectrum.** Per-hallmark R² over the 50 MSigDB hallmarks under 5-fold CV, organised
  by "morphological grounding": G2M Checkpoint R² = 0.78, IFN-γ Response R² = 0.75 (R² > 0.5 = morphologically
  visible); Oxidative Phosphorylation R² = 0.22, Fatty Acid Metabolism R² = 0.20 ("no detectable signal in H&E").
- **Clinical + cross-cohort.** POSEIDON trial (n = 265, NSCLC, inference-only): H&E-predicted squamous scores
  recapitulate NSCLC subtype; predicted IFN-γ mirrors PD-L1 tumour-cell expression groups. Domain adaptation on
  TCGA-LUAD (n = 223) and TCGA-BRCA (n = 1,042): **transfer degrades sharply — R@10 falls to 10.8% and 4.0%
  respectively** — and fine-tuning on 5% of target data recovers much of it.
- **Theory.** Essentially none. Verified absent from the full text: **"identifiab", "disentangle", "causal
  representation", "block-identifiab".** No latent slot structure (a single shared CLIP latent; all gene sets
  query the same space), no abstention, no task routing, no multi-interpretation.

### 2. The precise overlap

**This is the most serious of the five, and it is more serious than the ledger recorded.** Winter et al. already
build, end to end and with clinical validation, the *entire system substrate* of the claim:

1. frozen multimodal (WSI + RNA) foundation-model trunk ✔
2. lightweight adapter rather than end-to-end retraining ✔
3. pathway-level addressing at query time, open-vocabulary over gene sets ✔
4. inference-time pathway activity prediction from H&E alone ✔
5. cross-cohort transfer *measured*, with an honest failure report ✔
6. data-efficient domain adaptation ✔
7. an honest per-pathway capability spectrum ✔

Beyond the plumbing, two of their *findings* pre-empt likely MORPHEUS results. (i) The **graduated predictability
spectrum** already publishes the "which pathways are readable from a frozen trunk" curve. (ii) Their organising
principle — predictability tracks **morphological footprint in H&E** — is a strong confounding explanation for any
MORPHEUS slot-quality result. If identified slots read out better, morphological grounding is the null hypothesis
you must control for, and Winter has already named it.

### 3. The defensible remaining delta

**DELTA SURVIVES, but substantially narrowed. Claim only the guarantee and the latent structure.**

- **No identifiability claim of any kind.** Verified absent from the full text. Winter's latent is a *single
  shared contrastive embedding space*; pathways are **queries into** that space, not **slots within** it. The
  block-identifiability guarantee — identified pathway-*group* latent slots, bounded against the JMLR-2024
  ceiling, obtained via perturbation sparse-mechanism-shift + known pathway-membership grouping — is untouched.
  Winter has the plumbing; MORPHEUS must supply the guarantee. That is the whole of the surviving delta on C5.
- **Prompting is gene-sets-via-ssGSEA, not natural language.** Every query must be expressible as a gene set
  scored against a reference RNA library. Free-form NL addressing of a slot ("is this tumour hypoxic?") is not
  supported. This is a real but *modest* delta — an NL→gene-set front-end is an obvious extension a reviewer will
  see immediately, so do not lean on it.
- **No perturbation/interventional data.** Purely observational contrastive alignment. The sparse-mechanism-shift
  identification route is unavailable to them by construction.
- **No abstention, no routing, no multi-interpretation.**

⚠️ Strategic consequences, stated plainly:
1. **Winter is now the mandatory baseline for C5**, not merely related work. A CLIP-aligned frozen-trunk adapter
   must be an arm in every experiment.
2. **Retire any claim to "the system-level combination."** Winter owns it, with clinical-trial validation.
3. **The R@10 → 10.8% / 4.0% transfer collapse is a gift**: it is a published, quantified failure of the
   *non*-identified frozen adapter on exactly the cross-cohort axis C4 targets. It is the strongest available
   motivation for the identifiability→transfer hypothesis — provided you reproduce it and then beat it.
4. **Control for morphological grounding.** Stratify results by Winter's spectrum, or a reviewer will attribute
   any gain to morphology rather than identifiability.

### 4. Confidence

**READ-FULL-TEXT** (arXiv HTML full text; abs page cross-checked).

---

## (d) Confound-aware elicitation benchmark with a VALIDITY CERTIFICATE, selectivity/MDL, and text-prior nulls

**Threat A: Hu, Tripodi, Naidoo, McGough, Chakraborti — "Probing, Fusion, and Trustworthiness: A Systematic
Evaluation of Foundation Model Representations for Multimodal Cancer Analysis," arXiv:2606.17115, submitted
15 Jun 2026.**
**Threat B: "VCBench" (single-cell) — could not be found to exist.**

### 1. What the papers actually do

**2606.17115.** A systematic *evaluation study* (explicitly **not** a benchmark release) of frozen FM
representations under distribution shift, on **two real-world commercial in-house oncology cohorts, IH-BC
(breast) and IH-NSCLC**. Two modalities: WSI and transcriptomics.

- **5 foundation models:** image — CONCH, UNI, Virchow, MUSK; omics — UCE (plus scVI and PCA baselines).
- **8 downstream tasks:** IH-BC — LOH, biomarker PR, PIK3CA, **Biopsy Site**, Subtype; IH-NSCLC — **Biopsy Site**,
  **Tumor Site**, TMB.
- **Method:** frozen FMs produce case-level representations; unimodal and multimodal heads (linear layers + MLPs)
  are trained on top. Three image–omics fusion strategies compared.
- **Trustworthiness:** **split conformal prediction** — prediction sets with coverage guarantees, plus a
  "rescue rate" measuring how often a wrong top-1 prediction still contains the true label in its set.
  "In the majority of cases where a point prediction fails, the true diagnosis remains recoverable within the
  prediction set."
- **Findings:** image and omics representations carry complementary signal; FM representations are competitive
  OOD; multimodal fusion helps mainly when no single modality dominates.

Verified **absent** from the full text: Hewitt–Liang selectivity / control tasks; MDL / Voita–Titov probing;
any site/scanner/batch **confound** analysis; any natural-language querying or elicitation; any text-prior
baseline; the word "emergent"; any benchmark or password-lock contribution. CONCH is image–text pretrained, but
only its **frozen visual** embeddings are used — no language path is exercised.

**The single most important finding for this collision:** 2606.17115 treats **Biopsy Site and Tumor Site as
legitimate downstream prediction tasks**, sitting in the main results table alongside PIK3CA and TMB —
*"information about cancer subtypes, biopsy sites, and biomarkers is used as the downstream prediction tasks."*
There is **no** discussion of site/scanner/batch as nuisance, and **nowhere** is a representation required to
*fail* at predicting a confound.

**VCBench.** Searched arXiv for `all:"VCBench"` — 4 hits, all irrelevant: a **venture-capital** LLM benchmark
(2509.14448) and a **multimodal-mathematics** benchmark (2504.18589), plus two VC follow-ups. A broad sweep of
`abs:"virtual cell" AND abs:"benchmark"` (25 results) returned **no artefact named VCBench**. The real virtual-cell
benchmark landscape is: **AssayBench** (2605.10876), **"Benchmarking virtual cell models for in-the-wild
perturbation response"** (2604.27646), the **CZI Virtual Cells Workshop** recommendations (2507.10502), and —
closest to the NL angle and *not previously on the ledger* — **SC-Arena** (2602.23199), "A Natural Language
Benchmark for Single-Cell Reasoning with Knowledge-Augmented Evaluation," five biological reasoning tasks.
**The ledger's "VCBench (single-cell, pre-registered baselines + confound control)" is unverifiable; SC-Arena is
the citation that should replace it.**

### 2. The precise overlap

With 2606.17115: both take a **multimodal cancer trunk**, freeze it, train **linear/shallow probes** on its
representations, and evaluate **under distribution shift** with **uncertainty quantification**. It is the closest
published instance of "probe a multimodal cancer foundation model and report what it knows." If MORPHEUS's pitch
reduces to that sentence, it is pre-empted.

With SC-Arena: overlap on **natural-language evaluation of a biological foundation model**, but single-cell
transcriptomics only, and knowledge-augmented LLM-judge scoring rather than a confound certificate.

### 3. The defensible remaining delta

**DELTA SURVIVES — the widest margin of the five.** Every distinguishing component of the claim is verified absent
from 2606.17115:

- **The validity certificate is not merely absent — it is inverted.** MORPHEUS requires a recovered direction to
  **FAIL** at predicting site/scanner while succeeding on biology. 2606.17115 **rewards** site predictability by
  scoring Biopsy Site and Tumor Site as tasks. This is a quotable, unambiguous contrast and the strongest single
  sentence available for the related-work section. (Their site results are also useful *evidence*: they
  demonstrate site is richly encoded in these FM representations, which is precisely why the certificate is
  needed.)
- **No probe selectivity / control tasks (Hewitt–Liang).** Their probes are unregularised task heads; probe
  capacity is uncontrolled, so "the representation knows X" is not separated from "the probe learned X."
- **No MDL / description-length probing (Voita–Titov).**
- **No natural-language elicitation.** Their pipeline never queries the model in language.
- **No text-prior null.** No GenePT-style baseline testing whether a pathway is recoverable from the *name* alone.
- **No emergence framing, no loss-indexed curves, no password-locked control.**
- **Not a benchmark.** Self-described as an evaluation study on two proprietary in-house cohorts — meaning it is
  also **not reproducible by others**, which strengthens the case for a public benchmark artefact.

⚠️ One caution: conformal prediction is now the established trustworthiness instrument in this exact setting.
Either adopt it or justify not doing so.

### 4. Confidence

**READ-FULL-TEXT** for 2606.17115 (arXiv HTML, two independent targeted passes).
**COULD-NOT-ACCESS / LIKELY-NONEXISTENT** for the alleged single-cell VCBench.
**READ-ABSTRACT-ONLY** (metadata level) for AssayBench, 2604.27646, and SC-Arena.

---

## (e) Amnesic-counterfactual CAUSAL-USE gate on GENERATED hypothesis cards

**Threat A: Elazar, Ravfogel, Jacovi, Goldberg — "Amnesic Probing: Behavioral Explanation with Amnesic
Counterfactuals," TACL 2021, arXiv:2006.00995.**
**Threat B: "Decode-gLM" — could not be found to exist.**

### 1. What the papers actually do

**Elazar et al.** Argues that probing cannot support behavioural conclusions, and offers an alternative that asks
*how information is used* rather than *what is encoded*. Amnesic Probing removes a property from a representation
via a causal intervention (**INLP**, iterative nullspace projection) and measures the effect on the model's task
behaviour. Applied to **BERT**, asking e.g. "is part-of-speech information important for word prediction?"
Headline finding: **conventional probing accuracy is not correlated with task importance**, and the authors call
for scrutiny of causal/behavioural claims drawn from probes.

The intervention target is the **masked-language-model's own predictive distribution** — i.e. an analysis of an
already-trained model's behaviour. **It is a post-hoc interpretability method, not a gate or filter on output.**
No generated artefact is accepted or rejected. NLP only; no biology, no multimodality.

**"Decode-gLM."** Verified absent: `all:"Decode-gLM"` → **0 results** on arXiv; `"Nucleotide Transformer" AND
"steering"` → **0 results**. The name does not correspond to any indexed paper. The real paper the ledger was
probably reaching for:

- **Sarwan Ali — "Causal dictionary learning reveals and validates transcription-factor binding features in
  genomic language models," arXiv:2607.19618, 21 Jul 2026.** Trains **sparse autoencoders on Nucleotide
  Transformer and DNABERT-2**, recovers thousands of monosemantic features mapping to TF sequence motifs, and
  validates them causally by *"ablating individual dictionary directions during the model's forward pass and
  measuring the induced shift in the model's own predictive distribution."* Includes validation protocols for
  **GC-composition confounds**. **DNA sequence only.** No INLP, no LEACE. **No gating or filtering of generated
  natural-language output** — it operates entirely at the activation-feature level.

### 2. The precise overlap

Elazar owns the **primitive**: erase a direction, measure whether behaviour changes, conclude encoded-vs-used.
MORPHEUS must cite it as the method source and can claim **zero** novelty on the mechanism itself. 2607.19618
extends the same erase-and-measure logic into a *biological* foundation model, which closes the "nobody has done
amnesic-style intervention in biology" gap — the sequence-genomics half of it, at least.

### 3. The defensible remaining delta

**DELTA SURVIVES — cleanly, on the application and the control flow, not the mechanism.**

- **Elazar is diagnostic; MORPHEUS is a gate.** Elazar produces a *finding about a model*. MORPHEUS proposes an
  **accept/reject rule on a generated artefact** — a hypothesis card is emitted only if the programme it names
  survives the amnesic test. Turning a post-hoc analysis into an inline production-time admission criterion is a
  genuine and clearly statable difference. No verified work does this.
- **The object gated is generated natural language.** Neither Elazar (BERT token distributions) nor 2607.19618
  (TF-motif dictionary features) gates a generated explanation. This directly addresses the ELK
  "human-simulator" hazard — the card reporting the whitelist rather than the trunk's latent state — which
  neither paper engages.
- **Modality and substrate.** Elazar: text. 2607.19618: DNA sequence, unimodal. MORPHEUS: a **multimodal
  (WSI + molecular) tumour trunk**, with the downstream target being survival / drug-response prediction.
- **Method stack.** LEACE (Belrose 2023) alongside INLP, plus the dual steering direction for sufficiency —
  neither comparison paper uses LEACE or pairs erasure with steering.

⚠️ Two notes. (1) 2607.19618 is *very recent* (Jul 2026) and should be added to the ledger — it is a closer
neighbour than the phantom "Decode-gLM" and it establishes the causal-ablation-in-bio-FM precedent. (2) Its
**GC-composition confound protocol** is a small but real precedent for collision (d)'s certificate idea in a
different domain; cite it there too rather than letting a reviewer surface it.

### 4. Confidence

**READ-ABSTRACT-ONLY** for Elazar (abstract verbatim + method characterisation; full text not retrieved — though
the post-hoc-vs-gate distinction is unambiguous from the abstract's own framing).
**COULD-NOT-ACCESS / LIKELY-NONEXISTENT** for "Decode-gLM".
**READ-ABSTRACT-ONLY** for arXiv:2607.19618 (abstract + targeted Q&A on the abs page).

---

## Summary table

| # | Collision | Verdict | One-line surviving delta |
|---|---|---|---|
| **(a)** | Biological-task auto-detect + multi-interpretation + abstain **vs** "Nguyen" Answer/Clarify/Abstain | **COULD-NOT-VERIFY** (threat unfindable) → **DELTA-SURVIVES vs real prior art** | Latent is the *biological task*, and abstention is **structural** ("no identified slot covers this"), not a confidence threshold — but the triage primitive itself is now published ≥3× (Baan 2605.25831, Khanmohammadi 2606.15910, Presacan 2607.13036) and must be conceded. |
| **(b)** | Identifiability→reliability/transfer **vs** SurvPath (2304.06819) | **DELTA-SURVIVES** (wide) | SurvPath's 331 pathway tokens are fixed a priori (281 Reactome + 50 Hallmarks); the word "identifiability" never appears; it has **no prompting interface** and **no external cohort** — it instruments neither dependent variable, so it is the control arm, not a scoop. |
| **(c)** | Frozen trunk + block-identifiable NL-addressable slots **vs** Winter (2606.29949) | **DELTA-SURVIVES — but badly narrowed** | Winter already ships the frozen WSI+RNA adapter, gene-set pathway prompting, cross-cohort transfer, and the predictability spectrum, with **zero** identifiability language: the *only* surviving claim is the **block-identifiability guarantee on latent slots**. Winter is now a mandatory baseline; the system-level-combination claim must be retired. |
| **(d)** | Confound validity certificate + selectivity/MDL + text-prior nulls **vs** 2606.17115 / "VCBench" | **DELTA-SURVIVES** (widest) | 2606.17115 scores **Biopsy Site and Tumor Site as downstream tasks** — the exact inverse of the certificate — and has no selectivity, no MDL, no NL elicitation, no text-prior null, no emergence framing, and is not a benchmark; the "single-cell VCBench" does not exist (substitute **SC-Arena 2602.23199**). |
| **(e)** | Amnesic causal-use gate on generated cards **vs** Elazar (2006.00995) / "Decode-gLM" | **DELTA-SURVIVES** | Elazar is post-hoc explanation of BERT behaviour, never a gate; "Decode-gLM" does not exist — the real neighbour, **arXiv:2607.19618** (SAE + causal ablation on Nucleotide Transformer), is DNA-only, uses no INLP/LEACE, and gates nothing. The delta is **accept/reject on a generated NL artefact** over a multimodal tumour trunk. |

### Ledger hygiene actions implied

1. **Delete** "Nguyen et al., *Answer, Clarify, or Abstain*" (C1 kill-risk) — unverifiable. **Replace** with Baan
   2605.25831 (closer threat), Khanmohammadi 2606.15910, Presacan 2607.13036.
2. **Delete** "Decode-gLM" (C8 kill-risk) — does not exist. **Replace** with arXiv:2607.19618.
3. **Delete** "VCBench (single-cell)" (C7 prior art) — does not exist. **Replace** with SC-Arena 2602.23199,
   AssayBench 2605.10876, and 2604.27646.
4. **Upgrade** Winter 2606.29949 from "kill-risk" to **required baseline**; downgrade C5's scope to the
   identifiability guarantee alone; add "control for morphological grounding" as an experimental requirement.
5. **Add** Winter's R@10 → 10.8%/4.0% cross-cohort collapse as the quantified motivation for C4.
