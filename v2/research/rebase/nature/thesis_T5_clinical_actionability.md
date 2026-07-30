# T5 — Clinical Actionability: Multimodal Prediction of Treatment Response with External Open-Cohort Validation

**Scout date:** 2026-07-29 · **Search mode:** WebFetch only (WebSearch exhausted) · **Verification policy:** every citation below was resolved this session via PubMed E-utilities, the GDC API, TCIA, arXiv abs pages, or Semantic Scholar. Anything not resolvable is explicitly marked **COULD-NOT-VERIFY** or **BACKGROUND-ONLY (unverified)**. No citation in this document was written from memory alone without a flag.

---

## 0. Verdict (read this first)

| Axis | Score | One-line |
|---|---|---|
| Novelty | **2 / 5** | "Multimodal → immunotherapy response" is already a *Nature 2025* headline result (MUSK) and a *Nature Cancer 2022* result (Vanguri/DyAM). The plain claim is occupied. |
| Feasibility with open data | **2 / 5** | **The paired-WSI ICI cohort does not exist in open access.** Everything runnable is either chemo-era TCGA or n = 64–126 external sets. |
| Ceiling as proposed | **2 / 5** | Would land as a solid npj Precision Oncology / Medical Image Analysis paper, not Nature. |
| Ceiling if reframed as an audit | **3 / 5** | A pre-registered, adequately-powered *falsification* of clinical-actionability claims is the only Nature-family-shaped object here (Nature Methods / Nature Medicine "reality check"), and it is a hard sell. |

**KILLER RISK (confirmed, not speculative):** Open-access immune-checkpoint-inhibitor cohorts (Hugo, Riaz, Liu, Van Allen, IMvigor210, Braun, Samstein/MSK-IMPACT, GENIE BPC) supply **RNA and/or panel DNA and clinical outcome — never whole-slide images**. Conversely, every open cohort that *does* have WSI + per-patient treatment + response is either (a) TCGA, which is chemo-era, or (b) n = 64–126. TCGA's own ICI content, queried live from the GDC API today: **ipilimumab 20 cases, pembrolizumab 4, nivolumab 1, durvalumab 1.** There is no bridge.

---

## 1. Has this been done? — the 5 closest works (all verified)

### 1.1 MUSK — *A vision-language foundation model for precision oncology*
Xiang J, Wang X, Zhang X, et al. **Nature** 2025;638(8051):769–778. doi:10.1038/s41586-024-08378-w. PMID 39779851. *(verified: PubMed efetch)*

**What it did:** Pretrained on 50M pathology images from 11,577 patients plus ~1B pathology text tokens using unified masked modelling on *unpaired* image and text, then contrastive alignment. Evaluated on 23 patch- and slide-level benchmarks. Abstract states verbatim: *"MUSK showed strong performance in outcome prediction, including melanoma relapse prediction, pan-cancer prognosis prediction and immunotherapy response prediction in lung and gastro-oesophageal cancers."*

**Why it hurts T5:** This is the exact sentence T5 wants to write, already in Nature, already at scale, already framed as "the clinically consequential endpoint." **The specific ICI cohort sizes and identities were not obtainable from the abstract this session — flagged for full-text follow-up**, but the headline claim is pre-empted regardless of whose cohorts they were. MUSK's response prediction is a supervised classifier on a frozen trunk, i.e. structurally identical to what T5 proposes.

### 1.2 Vanguri et al. — *Multimodal integration of radiology, pathology and genomics for prediction of response to PD-(L)1 blockade in NSCLC*
Vanguri RS, Luo J, Aukerman AT, Egger JV, Fong CJ, Horvat N, et al. **Nature Cancer** 2022;3:1151–1164. doi:10.1038/s43018-022-00416-8. PMID 36038778. *(verified: PubMed esummary)*

**What it did:** The DyAM (dynamic attention with masking) model fusing CT radiomics + H&E pathology + targeted genomics to predict PD-(L)1 blockade response in NSCLC at MSKCC (~247 patients, **n BACKGROUND-ONLY, unverified this session**). Handles per-patient missing modalities via attention masking — the same problem MORPHEUS solves with `DenseAdapter` + availability mask.

**Why it hurts T5:** This is the canonical "multimodal ICI response with paired images" paper, and it establishes both the architecture pattern *and* the reason it is hard to reproduce: the cohort is institutional. It is the proof that the T5 idea works **and** the proof that you cannot do it on open data.

### 1.3 ViGNet — *A clinical data-supported deep learning approach for NSCLC immunotherapy response prediction in digital pathology*
**Medical Image Analysis** 2026. doi:10.1016/j.media.2026.104154. PMID 42269196. *(verified: PubMed efetch)*

**What it did:** Integrates WSI features with **gene expression profiles** and specific clinical modalities; reports 82.55% discrimination for immunotherapy response; validated on *one internal and two independent external cohorts*. (Per-cohort n not stated in abstract.)

**Why it hurts T5:** This is T5's exact modality triple (WSI + bulk RNA + clinical), exact endpoint (ICI response), exact validation design (external cohorts), published in a top venue two months before this scout. T5 as literally stated is **scooped**.

### 1.4 ENLIGHT-DP — *Prediction of clinical outcomes of advanced cutaneous squamous cell carcinoma to PD-1 inhibition directly from histopathology slides using inferred transcriptomics*
**Frontiers in Immunology** 2026. doi:10.3389/fimmu.2026.1822422. PMID 42212120. *(verified: PubMed efetch)*

**What it did:** n = 38 advanced cSCC treated with cemiplimab. Predicts response **from H&E via inferred transcriptomics** — i.e. WSI → predicted expression → response. PPV 84.2%; PFS stratification HR = 0.22.

**Why it matters:** This is precisely MORPHEUS's WSI→bulk-molecular pipeline pointed at a treatment endpoint. It is also a cautionary tale: n = 38, single-arm, retrospective, with an HR of 0.22 that is not credibly estimable at that sample size. It shows the *literature bar in this niche is low* — which cuts both ways: easy to publish, impossible to publish at Nature tier.

### 1.5 HEiST — *Histopathology-inferred spatial transcriptomics characterizes the tumor microenvironment in 1,500 head and neck tumors and predicts clinical outcomes*
**bioRxiv** 2026. doi:10.64898/2026.05.16.725687. PMID 42388794. *(verified: PubMed efetch)*

**What it did:** Infers spatial transcriptomics from H&E across 1,500 HNSC tumors from two public datasets plus a CCRT-treated cohort and an **immunotherapy-treated cohort**; validated on two independent external ST cohorts; claims immunotherapy-response predictors that *"markedly surpass FDA-approved biomarkers."*

**Why it matters:** This is the T1 (spatial) + T5 (response) fusion, at n = 1,500, six months old, on public data. If MORPHEUS pursues the ST-mediated route to response prediction, HEiST is the incumbent.

### Also-relevant (verified, second ring)
- **PORPOISE** — Chen RJ, Lu MY, Williamson DFK, et al. *Pan-cancer integrative histology-genomic analysis via multimodal deep learning.* **Cancer Cell** 2022;40(8):865–878.e6. PMID 35944502. 14 cancer types, WSI + molecular fusion. **Endpoint is survival, not response** — this is the architectural incumbent T5 must beat, and the reason "pan-cancer WSI+omics fusion" earns zero novelty.
- **Anagnostou V, et al.** *Multimodal genomic features predict outcome of immune checkpoint blockade in NSCLC.* **Nature Cancer** 2020. doi:10.1038/s43018-019-0008-8. PMID 32984843. Genomics-only multimodal ICI outcome.
- **Farahmand S, Fernandez AI, Ahmed FS, Rimm DL, Chuang JH, Reisenbichler E, Zarringhalam K.** *Deep learning trained on H&E tumor ROI predicts HER2 status and trastuzumab treatment response in HER2+ breast cancer.* **Mod Pathol** 2022;35(1):44–51. PMID 34493825. 188 WSIs for HER2 (AUC 0.90 CV, 0.81 on TCGA); **187 pre-treatment HER2+ trastuzumab-treated samples, AUC 0.80 in 5-fold CV.** The public release is smaller than the paper (see §3).
- **IMPRESS** — *Artificial intelligence reveals features associated with breast cancer neoadjuvant chemotherapy responses from multi-stain histopathologic images.* **npj Precision Oncology** 2023;7:14. doi:10.1038/s41698-023-00352-5. PMID 36707660. H&E + multiplex IHC WSI → NAC response; AUC 0.8975 (HER2+), 0.7674 (TNBC). *(cohort composition n=126, 62 HER2+/64 TNBC is **BACKGROUND-ONLY, unverified** — the abstract retrieved did not state n.)*
- **COMPASS** — pathway-concept-bottleneck cross-indication ICI response, reported ~76.5% on held-out indications. medRxiv 2025.05.01.25326820. *Carried from `lit/l11_benchmarks_confound.md` #23; **not re-verified this session** — treat as preprint-grade.*

**Verdict on Task 1:** T5's literal proposal is **occupied in triplicate** (MUSK Nature 2025 / Vanguri Nat Cancer 2022 / ViGNet MedIA 2026). There is no version of "we predict treatment response from multimodal tumor state and validate externally" that is novel in 2026.

---

## 2. Is the core claim true / plausible?

### 2.1 "The clinically consequential endpoint is treatment response, not pathway correlation" — TRUE, and it is the right instinct
This is correct as a *values* statement and it is the correct diagnosis of what is wrong with the current MORPHEUS evaluation. Nobody should dispute it. It is also, unfortunately, not a research contribution — it is a target selection, and the target is crowded.

### 2.2 "Multimodal tumor state predicts response better than unimodal" — PLAUSIBLE but likely SMALL, and MORPHEUS's own data argues against a large effect
Your own established finding (i) is the strongest available prior on this: **within-cancer, control-adjusted WSI→bulk-molecular signal is ~+0.07 Pearson and method-invariant across every method including baselines.** If morphology adds only +0.07 to *molecular state* — the easier, denser, more directly-observable target — then the honest prior for morphology adding to *treatment response* (a sparser, noisier, more clinically-confounded binary label mediated through molecular state) is **smaller, not larger**.

Concretely: to detect a ΔAUC of 0.03 over a clinical baseline at 80% power with α=0.05 in a balanced binary-response cohort, you need roughly **n ≳ 1,500–2,500 patients** (DeLong-based, correlated ROCs, moderate baseline AUC ≈ 0.65). The largest external validation cohort available to you is **n = 126**. That cohort can detect ΔAUC ≈ 0.15 at best. **The study is underpowered by an order of magnitude for the effect size your own pilot implies.** This is the quantitative form of the killer risk and it should be the first line of any go/no-go memo.

### 2.3 The T1 sub-question: do spot-level ST targets yield materially higher morphology-predictability than bulk targets?

**Answer: NO PUBLISHED HEAD-TO-HEAD EVIDENCE — and the comparison as usually reported is not well-posed.** I could not verify any paper that runs both target types under a common protocol.

- **HEST-1k** (Jaume et al., NeurIPS 2024 D&B Spotlight, arXiv:2406.16192 — *abstract verified this session*) is the reference spot-level harness: 1,229 ST profiles from 153 cohorts, 26 organs, 367 cancer samples across 25 cancer types, 2.1M expression–morphology pairs. **The arXiv abstract reports no Pearson values**; the per-gene correlations circulating in secondary sources were **NOT verifiable this session** and must not be quoted.
- **SEQUOIA** (Zhu et al., Nat Commun 2024, `lit/l04` #31) is the bulk-target incumbent. Its correlations are computed *across patients*.
- **The denominators differ, so the numbers are not commensurable.** Spot-level Pearson is computed across spots *within a slide*: a large fraction of that variance is tumour-vs-stroma-vs-immune compartment identity, which is trivially readable from morphology. Bulk Pearson is computed across patients, where between-patient variance is dominated by factors morphology cannot see (germline, clonal architecture, sampling, batch). **A higher spot-level number does not mean spot-level targets carry more morphology-recoverable biology; it usually means the easy axis of variation is inside the frame.** The known failure mode — a within-slide mean-expression predictor being competitive with trained models — is the spot-level analogue of your own control-adjustment finding. I could **NOT** locate a citable benchmark paper making this argument explicitly (**COULD-NOT-VERIFY**; PubMed returned 0 and arXiv API 429'd/timed out repeatedly this session).

**Consequence for T5:** routing response prediction through inferred spot-level ST (the HEiST strategy) does **not** get you a free accuracy uplift. It gets you a different, prettier intermediate representation with the same patient-level information ceiling.

---

## 3. Open-access data audit — the part that decides this

All counts below are **live queries run 2026-07-29**, not recalled.

### 3.1 GDC / TCGA — the only large open WSI + treatment + response resource

`GET api.gdc.cancer.gov/cases?facets=diagnoses.treatments.treatment_outcome` filtered to program = TCGA:

| treatment_outcome | cases |
|---|---|
| complete response | 1,643 |
| progressive disease | 598 |
| stable disease | 327 |
| partial response | 225 |
| treatment ongoing | 928 |
| unknown | 548 |
| not reported | 300 |
| no measurable disease | 177 |
| persistent disease | 22 |
| no response | 12 |
| treatment stopped due to toxicity | 1 |
| **_missing** | **7,509** |

Filtered to **RECIST-4 outcome AND `files.data_type = "Slide Image"`** → **2,898 cases total**, of which **TCGA 2,581**, CGCI 251, HCMI 38, CCG 28.

**Intersection with MORPHEUS's 6,192 paired WSI+RNA patients is NOT yet computed.** Naive estimate 2,581 × (6,192 / ~11,000) ≈ **1,400–1,500 patients**. *This join must be run on disk before any planning decision — it is a 20-minute job and it is the single most decision-relevant number in this document.*

**Which drugs?** Top therapeutic agents among TCGA cases with a RECIST-4 outcome:

| agent | cases | | agent | cases |
|---|---|---|---|---|
| cisplatin | 494 | | capecitabine | 107 |
| fluorouracil | 312 | | bevacizumab | 90 |
| carboplatin | 290 | | doxorubicin | 85 |
| paclitaxel | 282 | | **ipilimumab** | **20** |
| tamoxifen | 235 | | **pembrolizumab** | **4** |
| cyclophosphamide | 221 | | **nivolumab** | **1** |
| gemcitabine | 184 | | **durvalumab** | **1** |
| temozolomide | 170 | | | |

**This is the killer, in one table.** TCGA accrual predates the checkpoint era. TCGA can support a *chemotherapy* response study (platinum n≈494 pre-WSI/RNA-join, realistically ~300–350 after) and nothing else.

**Additional TCGA caveats that are not cosmetic:**
- "Complete response" (1,643 — 59% of all labelled cases) in TCGA frequently encodes *post-surgical disease status*, not drug response. The label is not RECIST-clean.
- Treatment records are free-text-derived, institution-heterogeneous, and missing for 7,509 / ~11,000 cases — and **missingness is site-correlated**, which is exactly the confound Howard et al. and de Jong et al. (in `lit/l11`) show pathology FMs exploit.
- Response is confounded with stage, resectability, and cancer type — the same structural confound that produced your 46–49% "cross-cancer cohort structure" artefact.

### 3.2 CPTAC — dead for response
`program = CPTAC`, 2,025 cases: complete response 493, persistent disease 248, unknown 75, **stable disease 33, partial response 9, progressive disease 6**, _missing 1,205. With PR + PD + SD = 48 cases, CPTAC cannot support a response endpoint under any framing. Its WSI + proteomics pairing remains valuable for T4-type work; for T5 it is not usable.

### 3.3 The complete inventory of open cohorts with **WSI + per-patient treatment + response**

Queried against the TCIA collections index (verified 2026-07-29):

| Resource | n subjects | Modalities | Endpoint | Access | Usable? |
|---|---|---|---|---|---|
| **TCGA (GDC)** | ~1,400–1,500 est. after WSI+RNA join | WSI, bulk RNA, clinical | RECIST-4, chemo era | Open | Only large option; chemo only; label noisy |
| **HER2-TUMOR-ROIS** (Yale) | **273 total; Yale trastuzumab-response sub-cohort = 85** (36 pCR responders / 49 non-responders) | H&E WSI (.svs) + tumour ROI XML, 40 GB | pCR to trastuzumab | **CC BY 4.0**, DOI 10.7937/E65C-AM96 | Yes — but n=85 |
| **Ovarian Bevacizumab Response** (TCIA) | **78** | Histopathology / WSI, Follow-Up, Measurement, Treatment | Bevacizumab response | Open | Yes — n=78 |
| **IMPRESS** (npj Prec Onc 2023) | ~126 *(unverified)* | H&E + multiplex IHC WSI | pCR to NAC (HER2+, TNBC) | Public | Yes — n≈126 |
| **Post-NAT-BRCA** (TCIA) | **64** | Histopathology / WSI, Treatment, Follow-Up | Post-neoadjuvant status | Open | Marginal |
| **AURORA-Metastatic-Breast-Multiomics** (TCIA) | **55** | WSI, IF, pathology detail | Metastatic breast multi-omics | Open | No response endpoint confirmed |

**Total open external-validation capacity for treatment response with paired images: ~350 patients, spread across four disease/treatment contexts that share nothing with each other.**

Note: the TCIA sweep was a single index query, so this is a **strong but not exhaustive** negative. A manual pass over TCIA's full collection list and over Zenodo/Synapse WSI releases is cheap insurance before killing the thesis outright.

### 3.4 The ICI cohorts everyone names — what they actually contain

**All entries in this block are BACKGROUND-ONLY and were NOT verified this session** (GEO/dbGaP/EGA/Synapse were not queried). Treat the *sizes* as approximate and the *"no WSI"* column as the load-bearing claim, which is corroborated by the TCIA negative in §3.3.

| Cohort | Modality actually released | WSI? | Access |
|---|---|---|---|
| Hugo 2016 melanoma anti-PD-1 (GSE78220) | bulk RNA-seq, ~28–38 pts | **No** | GEO, open |
| Riaz 2017 melanoma nivolumab (GSE91061) | bulk RNA-seq pre/on-tx, ~68 pts | **No** | GEO, open |
| Liu 2019 melanoma anti-PD-1 | WES + RNA, ~144 pts | **No** | dbGaP — **controlled, not open** |
| Van Allen 2015 melanoma ipilimumab | WES, ~110 pts | **No** | dbGaP — controlled |
| IMvigor210 (Mariathasan 2018) urothelial atezolizumab | bulk RNA + clinical, ~348 pts via `IMvigor210CoreBiologies` R package | **No** | R package open; raw via EGA controlled |
| Braun 2020 ccRCC (CheckMate) | WES + RNA, ~300 pts | **No** | dbGaP — controlled |
| Samstein 2019 MSK-IMPACT ICI | panel DNA (TMB) + OS, ~1,600 pts | **No** | cBioPortal, open |
| **AACR GENIE BPC** | panel DNA + PRISSMM-curated treatment/PFS/OS. **NSCLC cohort n = 1,846** (Choudhury NJ, et al. *Clin Cancer Res* 2023;29(17):3418–3428, PMID 37223888 — **verified**). 10 solid-tumour cohorts overall (Acebedo A, et al. *ESMO Real World Data Digit Oncol* 2025;7:100097 — **verified**) | **No** | Synapse, registration + DUA |
| ORIEN / AVATAR | WES + RNA + real-world outcomes | Some, institutionally | **Consortium membership required — NOT open** |

**The pattern is structural, not accidental.** ICI cohorts came from trials and trial correlative science, where the deliverable was a sequencing assay on a core biopsy. Slides stayed with the pathology department under institutional control. This will not change on your timeline.

---

## 4. Highest honest claim, and what falsifies it

### 4.1 What T5 CANNOT honestly claim
- ❌ "First multimodal prediction of immunotherapy response" — MUSK, Vanguri, ViGNet.
- ❌ "Externally validated on open ICI cohorts" — those cohorts have no images.
- ❌ Any clinical-utility claim from n ≤ 126 external validation.
- ❌ Any ΔAUC < 0.10 claim on the available external cohorts — it is not detectable.

### 4.2 The highest honest claim (as a positive result)
> *On the complete universe of open-access cohorts with paired H&E, molecular profiling, per-patient therapy and per-patient response, a pan-cancer multimodal tumour-state representation transfers to therapy-response prediction in cancer types held out of training, and does so with a control-adjusted margin over stage/subtype/cancer-type clinical baselines of Δ — reported with its confidence interval and its detectable-effect floor.*

This is honest, and it is a **Medical Image Analysis / npj Precision Oncology / Nature Communications** paper at best. It is not Nature.

### 4.3 The highest-ceiling honest claim (as an audit — the only Nature-family shape here)
> *Clinical-actionability claims for multimodal tumour-state representations are not currently testable on open data, and where they are testable they do not survive control adjustment. We (a) enumerate the entire open-access universe of WSI + treatment + response — ~1,500 TCGA chemo-era patients plus ~350 external patients across four disparate contexts; (b) show that on this universe the control-adjusted multimodal margin over clinical baselines is indistinguishable from zero and **method-invariant** — identical for a foundation-model trunk, a linear probe, and a stage-and-subtype-only model; (c) show this is the same method-invariance we previously demonstrated for WSI→molecular prediction (+0.07 Pearson), i.e. a single phenomenon, not two; and (d) derive the sample sizes and cohort designs that would be required to detect a real effect, establishing a power floor the field is currently ~10× below.*

**Why this could reach Nature-family tier:** it is a *new evaluation paradigm* claim (one of the four things the user wants), it is backed by two independent method-invariance results (molecular target + clinical target), and it converts a null into a prescription. Precedent for this genre landing well exists in the field (site-signature confounding, single-cell FM baseline failures — `lit/l11` #7, #9, #20, #21). **Why it probably does not:** null results about model performance are much harder to place than null results about *data* confounding, and reviewers will correctly note that the negative may reflect the poverty of open data rather than a property of multimodal representations. That ambiguity is intrinsic and cannot be engineered away.

### 4.4 What would falsify each
| Claim | Falsifier |
|---|---|
| Multimodal adds to response prediction | Δ over a stage + subtype + cancer-type logistic baseline whose 95% CI includes 0, on the held-out-cancer split. **Pre-register the baseline before looking.** |
| The effect is real, not cohort structure | Δ collapses to ~0 after site/institution-stratified adjustment (TCGA TSS code), or a site-only model matches the multimodal model. |
| Method-invariance generalises from molecular to clinical targets | A frozen H-Optimus + BulkFormer trunk **beats** a mean/linear/clinical baseline by a margin outside CI on ≥2 of the 4 external cohorts. This would kill the audit framing and revive the positive framing. |
| Open data suffices | Any open cohort with ≥500 ICI-treated patients and paired WSI surfaces. *(Actively watch for this; it would change everything.)* |
| Chemo response is a usable proxy for actionability | Platinum-response models trained on TCGA fail to transfer to the ovarian-bevacizumab and breast-NAC external cohorts at any margin. |

---

## 5. Ratings

| Dimension | Score | Justification |
|---|---|---|
| **Novelty** | **2 / 5** | MUSK (Nature 2025), Vanguri (Nat Cancer 2022), ViGNet (MedIA 2026), HEiST (2026), ENLIGHT-DP (2026), COMPASS. The literal claim is scooped three times over. The +1 above floor is for the held-out-cancer × treatment-response × control-adjusted combination, which I did not find published (PubMed returned 0 for that conjunction) — but absence in one query is weak evidence. |
| **Feasibility with open data** | **2 / 5** | Runnable tomorrow on TCGA (~1,400–1,500 patients, chemo). **Not** creditable as clinical actionability: external validation caps at ~350 patients across four unrelated contexts, and the study is ~10× underpowered for the effect size MORPHEUS's own pilot implies. +1 above floor only because the TCGA join genuinely exists and the four external cohorts are genuinely open and CC-licensed. |
| **Ceiling** | **2 / 5** as proposed; **3 / 5** if reframed as the §4.3 audit | Positive framing lands mid-tier. Audit framing has a Nature-family shape and rides two real prior results of yours, but null-result placement risk is high and the "poverty of data vs. property of models" confound is unresolvable. |

---

## 6. Recommendation

**Do not run T5 as a treatment-response prediction thesis.** The data does not exist and the claim is scooped.

**Two salvage routes, in preference order:**

1. **Fold T5 into the audit thesis as its clinical arm.** Your two established findings (method-invariant +0.07 on molecular targets; rank collapse decoupled from benchmark score) plus a third — *method-invariance also holds for the clinically consequential endpoint, on the entire open universe* — is a coherent three-legged evaluation paper. T5's real contribution is the **data-universe enumeration in §3**, which is genuinely useful and which nobody has published as a systematic audit. That table is worth writing up regardless of what happens to the rest.

2. **Keep response as a downstream *probe* of a promptable trunk (A1/A5), not as the thesis.** "Response prediction is one query among many on a frozen representation, and it degrades gracefully to held-out cancers" is defensible, unscooped in that framing, and does not require the claim to carry the paper.

**Immediate next actions (cheap, high information):**
- [ ] Run the TCGA join on disk: how many of the 6,192 paired WSI+RNA patients have a RECIST-4 outcome and a named agent? Break down by agent and by dev/test cancer split. **This one number decides the salvage route.**
- [ ] Download and inventory the four external cohorts (HER2-TUMOR-ROIS 10.7937/E65C-AM96; Ovarian Bevacizumab Response; IMPRESS; Post-NAT-BRCA) — confirm actual label availability, not just modality listings.
- [ ] Retrieve MUSK full text and extract its immunotherapy cohort identities and sizes. If MUSK used public data, that route is closed; if institutional, the §3 audit gains force.
- [ ] Retrieve ViGNet full text: which two external cohorts? If they are open, they belong in §3.3 and this document's central negative weakens.
- [ ] Exhaustive TCIA + Zenodo + Synapse sweep for WSI + treatment-response releases, to convert §3.3's strong negative into an exhaustive one.

---

## 7. Verification log

**Verified live this session:**
GDC API facet queries (TCGA treatment_outcome; TCGA outcome ∩ Slide Image by program; TCGA therapeutic_agents; CPTAC treatment_outcome) · TCIA collections index · TCIA HER2-TUMOR-ROIS collection page · PubMed efetch/esummary for PMIDs 39779851, 36038778, 32984843, 34493825, 35944502, 37223888, 41647353, 42212120, 42269196, 42388794, 42027950 · Semantic Scholar for IMPRESS (PMID 36707660) · arXiv abs 2406.16192 (HEST-1k).

**COULD-NOT-VERIFY this session:**
- CHIMERA challenge (multimodal bladder / BCG response) — grand-challenge URL 404, PubMed 0 hits, arXiv API 429. Existence unconfirmed; **do not cite**.
- Any benchmark paper explicitly critiquing spot-level ST evaluation vs. mean-expression baselines — PubMed 0 hits, arXiv API timeouts. The §2.3 argument stands on reasoning, not citation.
- IMPRESS cohort composition (n=126 / 62 HER2+ / 64 TNBC) — abstract retrieved did not state n.
- MUSK and ViGNet immunotherapy cohort identities and sizes — abstracts insufficient.
- Vanguri n ≈ 247 — not stated in retrieved esummary.
- COMPASS (medRxiv 2025.05.01.25326820) — carried from prior lit sweep, not re-verified.

**BACKGROUND-ONLY (asserted from model knowledge, explicitly unverified):** all entries in §3.4 except the two GENIE BPC citations. The individual accession numbers and patient counts for Hugo/Riaz/Liu/Van Allen/IMvigor210/Braun/Samstein should be checked against GEO/dbGaP/EGA before appearing in any manuscript.

**Tooling note:** arXiv API and Semantic Scholar both rate-limited (429) or timed out repeatedly; PubMed E-utilities and the GDC API were reliable. Prefer those two for future sweeps in this session.
