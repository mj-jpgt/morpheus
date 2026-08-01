# FEASIBILITY KILL — T2 "virtual perturbation" (decomposition/benchmark reframe)

**Date:** 2026-07-29 · **Assessor role:** feasibility assassin (no wet lab, open-access data only, one A100-40GB)
**Search constraint:** WebSearch exhausted. All external claims below verified via **WebFetch** against
PubMed E-utilities and the HuggingFace datasets API. Unverifiable items are marked **COULD-NOT-VERIFY**.
All on-disk numbers were **recomputed from the actual files this session**, not taken from the thesis doc.

**VERDICT: KILLED.**

Not killed because it is hard. Killed because the central deliverable — *"report the patient-residual
per perturbation"* — is (1) **unscoreable against any ground truth that exists or can be downloaded**,
and (2) **arithmetically a restatement of a number the project already has**. The "it survives a null
result" defence is inverted: the apparatus produces a null **by construction**, which is vacuous, not
rigorous.

---

## 0. What I verified on disk (recomputed, not quoted)

| Asset | Verified reality | Path |
|---|---|---|
| TCGA paired WSI(H-Optimus) + bulk RNA(BulkFormer) | **6,443 patients, 32 cancer types** (thesis says 6,192 — consistent post-QC) | `data/processed/rna/bulkformer/tcga_bulkformer_embedding_metadata.parquet` ∩ `data/processed/wsi/tcga_ut_hoptimus0_patch_embeddings_metadata.parquet` |
| Per-cancer-type n | **median 160**, min 28 (DLBC), max 733 (BRCA) | same |
| Replogle GWPS | `X` = **(11,258 perturbations × 8,248 genes)** — **pseudobulk, one row per perturbation, one cell line** | `PRISM/data/perturbseq/K562_gwps_normalized_bulk_01.h5ad` |
| Replogle K562-essential / RPE1 | 80 MB / 95 MB h5ad, same pseudobulk form | same dir |
| DepMap CRISPR gene effect | **1,178 cell lines × 17,916 genes** (429 MB) | `PRISM/data/ccle/raw/CRISPRGeneEffect.csv` |
| DepMap lineage depth | **median 34 cell lines/lineage**; Prostate **10**, Thyroid **11**, Liver 24, Cervix 18, Eye 15 | `PRISM/data/ccle/raw/Model.csv` (OncotreeLineage, restricted to CRISPR-screened lines) |
| CCLE / GDSC processed | 437 MB / 131 MB pickles present | `PRISM/data/{ccle,drug_response}` |
| **Free disk on C:** | **49 GB of 927 GB (95% full)** | `df -h` |

Engineering/compute is **not** the kill. H-Optimus patch features are precomputed; Replogle is a
375 MB pseudobulk matrix; DepMap is a 429 MB CSV. Every model in this thesis trains in hours on one
A100-40GB. I looked for a compute kill and there isn't one. The kill is epistemic and statistical.

---

## 1. FATAL OBSTACLE #1 — the residual has no ground truth, and its DepMap score is zero *by construction*

Write the decomposition exactly as the thesis specifies, for patient *i* in cancer type *c*, perturbation *p*:

```
ŷ_ip  =  μ_c(p)        [cancer-type mean]
       + λ_ic(p)       [predicted-lineage-explained]
       + ε_ip          [PATIENT RESIDUAL  ← the entire deliverable]
```

To report ε as an *information* quantity (not merely as prediction variance) you must score it against
a measured per-patient response `y*_ip`. Inventory of every candidate instrument:

| Candidate ground truth | Entity measured | Has H&E? | Varies across patients *within* a cancer type? | Can score ε? |
|---|---|---|---|---|
| DepMap CRISPR (1,178 lines × 17,916 genes) | cell line | **no** | **no** — collapses to a lineage constant | **No** |
| GDSC / CCLE drug response | cell line | **no** | **no** | **No** |
| Replogle GWPS (on disk) | K562 / RPE1 | **no** | **no** (2 contexts total) | **No** |
| Tahoe-100M (50 lines, 1,100 drugs) | cell line | **no** | **no** | **No** |
| CPTAC proteomics | patient | yes (some) | yes — but it is **not a perturbation response** | No (tests H1 again) |
| Survival | patient | yes | yes — but grade/stage-driven | Dropped by the thesis itself, correctly |
| Organoid/PDX CRISPR on the same tumor | patient | yes | yes | **Does not exist** (§4) — and needs wet lab |

The decisive line is row 1. Mapped to patients, DepMap supplies exactly one number per *(lineage, gene)*.
That target is **constant in *i* within *c***. The residual ε_ip is, by definition, mean-zero in *i* within *c*.
**The correlation of a mean-zero-in-*i* quantity with a constant-in-*i* target is identically zero — for
every perturbation, at every sample size, whether or not the underlying biology is real.**

So the proposed validation apparatus is **mathematically orthogonal to the proposed deliverable**.
DepMap can score components (a) and (b) and *only* (a) and (b). Component (c) — the paper — is
unscoreable, not noisily, but structurally.

**This detonates the thesis's own load-bearing defence.** "It is bulletproof precisely because it
survives a null result" is only true when the measurement *could have come out positive*. Here the
DepMap-scored residual is pinned to zero by the algebra of the estimand. That is a **vacuous null**,
and a reviewer at any of the target venues will see it in one reading. A rigorous null requires a
demonstrated-sensitive instrument; there is no instrument.

**Entity disjointness is the root cause and it is not fixable by download:** every entity in the world
that has a *measured perturbation response* (DepMap/GDSC/Replogle/Tahoe cell lines) has **no histology**,
and every entity that has *histology* (TCGA patients) has **no measured perturbation response**. Cell
lines have no H&E — they are cultured, and generating H&E-comparable sections of them is wet lab.
There is no open-access dataset in which one entity has both. This is the epistemic ceiling the scout
doc itself named in §3.3 and then walked past when it kept the benchmark as the deliverable.

---

## 2. FATAL OBSTACLE #2 — the only computable surrogate makes the study information-free (data-processing inequality)

The single escape from §1 is RNA-derived pseudo-labels: fit `g: RNA → dependency` on DepMap cell lines
(the DeepDEP recipe, verified: PMID 34417181), apply `g` to each patient's bulk RNA to manufacture a
per-patient target `ŷ*_ip = g(RNA_i)`, then ask whether WSI recovers its within-cancer residual.

This is computable on the assets in hand. It is also **guaranteed to produce no new information**, and
I can show that from the shape of the files on disk:

- `g` is a **fixed function with no patient-specific parameters**. Verified: `K562_gwps_normalized_bulk_01.h5ad`
  holds `X = (11,258 × 8,248)` — **one response vector per perturbation, in one cell line**. The
  "perturbation-response manifold" has **n = 1 in the context dimension** (n = 2 distinct lines across
  all three on-disk files). A manifold with no context axis cannot supply patient-specific interventional
  content; "placing a tumor on it" is a fixed linear/kernel readout of the tumor's predicted transcriptome.
- Therefore `WSI → predicted molecular state → predicted dependency` is a **Markov chain**, and
  `I(WSI ; dependency | cancer type) ≤ I(WSI ; molecular state | cancer type)`.
- The project has **already measured the right-hand side**: within-cancer, control-adjusted WSI→bulk-molecular
  Pearson ≈ **+0.07, method-invariant across every method including baselines**.

**So the headline number of the T2 paper is already in hand, and the entire Perturb-seq/DepMap apparatus
is a coordinate change applied to it.** The decomposition of predicted *dependency* into (a)/(b)/(c) is
isomorphic to the decomposition of predicted *expression* into (a)/(b)/(c), which the project has done.
Eleven thousand perturbations is presentational width, not evidential depth.

A reviewer's one-line objection: *"You have measured your own WSI→RNA correlation a second time and
passed it through a fixed map that can only shrink it."* There is no rebuttal available, because the
shrinkage direction is a theorem, not an empirical question.

---

## 3. FATAL OBSTACLE #3 — the per-perturbation deliverable is unpowered by ~5×, and unfixable with open data

Even granting the surrogate in §2, the stated deliverable ("report the residual **per perturbation**")
fails a power calculation on the real cohort. Using SE(r) ≈ 1/√(n−3) on the **verified n = 6,443**
(SE = 0.01246), 80% power, two-sided:

| Estimand | Multiple-testing burden | Min detectable r | n needed at r = 0.03 | n needed at r = 0.07 *(impossible best case)* |
|---|---|---|---|---|
| One pooled scalar | 1 | **0.035** | 8,711 | 1,600 |
| ~50 perturbation classes | 50 | **0.052** | 18,957 | 3,482 |
| **11,258 Replogle perturbations** | 11,258 | **0.068** | **32,756** | 6,016 |
| **17,916 DepMap genes** | 17,916 | **0.069** | **33,925** | 6,231 |
| **Per cancer type** (median n = 160) | 1 | **0.223** | 8,711 *per cancer type* | 1,600 *per cancer type* |

Read the r = 0.07 column carefully, because it is the charitable case and it still fails: **0.07 is the
project's *un-attenuated* WSI→RNA ceiling**, achievable only if the RNA→dependency map `g` loses
*nothing* — which §2 shows is impossible. At that physically unreachable upper bound the genome-wide
design lands at r_min = 0.069 vs. effect 0.07, i.e. **~50% power at best**. At any realistic attenuation
(r ≈ 0.02–0.04) genome-wide power is **≈ 0**.

- The **per-perturbation** deliverable needs **~34,000 paired-WSI patients**. TCGA is the largest open
  paired WSI + bulk-RNA cohort in existence; it yields **6,443**. Adding all of CPTAC (~1,000) gets to
  ~7,400 — still **4.6× short**. There is no open-access route to 34,000. This is not a funding or
  effort problem; the cohort does not exist.
- The **perturbation-class** fallback ("concentrated in a nameable class") needs **~19,000**. Also
  unreachable. Class-pooling reduces the multiple-testing burden but **not** the patient-side SE — the
  patients are the same 6,443.
- The **per-cancer-type** resolution that any clinical framing would require is hopeless by two orders
  of magnitude: median n = 160 → nothing below r = 0.223 is detectable.

**What the assets can actually deliver is one pooled scalar**, marginally powered, with no ground truth
behind it (§1) and already known in a different coordinate system (§2). One unscoreable scalar is not a
benchmark and is not a Nature-tier deliverable.

---

## 4. FATAL OBSTACLE #4 — the validation a reviewer will demand cannot be produced, and I checked

Any claim of the form "this *patient's* tumor is dependent on gene X" invites exactly one demand:
show a functional readout on patient-derived material. Verified searches (PubMed E-utilities, this session):

- `"patient-derived organoid" AND "CRISPR screen" AND cancer` → **1 result total**, PMID **40838977**,
  Venkadakrishnan et al., *"Epigenetic Derepression of PROX1 Promotes Neuroendocrine Prostate Cancer
  Progression"*, **Cancer Research 2025** — a single-gene mechanism paper, **not a paired resource**.
- `"patient-derived xenograft" AND "deep learning" AND (histology OR "whole slide") AND (drug response
  OR treatment response)` → **2 results**: PMID **40595000** (*Sci Rep* 2025) and PMID **36960342**
  (*Front Med* 2023, "Data augmentation and multimodal learning for predicting drug response in
  patient-derived xenografts from gene expressions and histology images"). Both are **drug** response,
  neither is a genetic-perturbation resource, and neither is at cohort scale.

There is **no open-access corpus** pairing human tumor histology with measured genetic-perturbation
response, at any n. With wet lab excluded by constraint, this reviewer demand is **unanswerable**. The
scout doc's §3.3 states this; it is a kill, not a caveat.

**Compounding — the positive-control circularity.** A null is only publishable if you can show the
pipeline detects signal where signal must exist. The only positive control available here is
proliferation/mitotic activity. But proliferation-coupled core-essentiality is *also the thesis's only
claimed positive finding* (§4.2 of the scout doc). **You cannot use the same axis as both the control
that legitimizes the null and the finding the null is contrasted against.** And if it is used as the
finding, the honest restatement is *"H&E predicts tumor grade, and grade predicts core essentiality"* —
both halves already known.

---

## 5. Secondary obstacles (real, but not individually fatal)

- **Cell-line→tumor gap, quantified.** Median **34** CRISPR-screened cell lines per lineage. The two
  TCGA cancer types with the largest patient counts after BRCA/LGG/LUSC/LUAD are **THCA (368 patients →
  11 thyroid lines)** and **PRAD (301 → 10 prostate lines)**. The lineage reference profile for the
  cancers where you have the most patients is estimated from ~10 lines. This degrades component (b),
  the *one* component that is scoreable. Not fatal on its own; fatal in combination with §1.
- **The recommended de-risk (Tahoe-100M swap) does not fit on this machine.** Verified via HF API:
  `tahoebio/Tahoe-100M`, **usedStorage 1.69 TB**, download size **338 GB**, 3,389 parquet files,
  CC0-1.0, 41,823 downloads. The `pseudobulk_differential_expression` subset is **1,026 shards at
  ~75–104 MB ≈ 92 GB**; `obs_metadata.parquet` alone is **2.29 GB**. Free space on C: is **49 GB**.
  So the single highest-value recommendation in the scout doc is **not executable without new storage**.
  Surmountable (buy a drive, or stream a subset) — flagging it as amber, not the kill. Note also that
  Tahoe is **chemical**, 50 cell lines, **still no histology** — it widens the lineage axis from 2 to 50
  but does nothing at all about §1.
- **Coladan-human3K availability.** The competitor's ~3,000-profile ST resource release channel was
  **COULD-NOT-VERIFY** this session.

---

## 6. What survives — the T1 side-finding is feasible and should absorb the effort

I checked whether the ST side-finding is data-feasible, and unlike T2 it is:

- **HEST-1k** (`MahmoodLab/hest`, HF API, verified this session): **1,276 spatial transcriptomic profiles,
  each linked and aligned to a Whole Slide Image at < 1.15 µm/px**, 26 organs, **398 cancer samples
  across 25 cancer types**, >1.5M expression/morphology pairs, 77,766 downloads. Openly downloadable.
- The pairing T2 needs does not exist; **the pairing T1 needs already exists, aligned, at usable
  resolution.** The mean-spot-subtracted, held-out-slide analogue of spot Pearson is a straightforward
  computation on one A100 over precomputed patch features.
- The scout doc's own finding stands: **no published like-for-like head-to-head** of spot-level vs. bulk
  targets under matched baseline adjustment (**COULD-NOT-VERIFY** that one exists). That absence is a
  real, cheap, scoreable, adequately-powered result — the estimand is defined, the ground truth is
  *measured* (real ST counts, not a pseudo-label), and the null is informative because the instrument is
  demonstrably sensitive (spot Pearson 0.230→0.431 proves the pipeline detects signal).

Everything §1–§4 denies T2, HEST-1k grants T1: a measured target on the same entity as the image.

## 7. Cheapest disconfirmation of *this* kill

If you want to overturn this verdict rather than accept it, one week settles it, and it is the scout
doc's own falsification test #1 — with the estimand fixed:

1. Build `g: RNA → DepMap gene effect` on the 1,178 lines (on disk).
2. Apply to all 6,443 patients' BulkFormer RNA → per-patient pseudo-dependency.
3. Residualize within cancer type; correlate the **WSI-predicted** residual against the **RNA-derived**
   residual, pooled, one pre-specified scalar.
4. **Kill confirmed if that pooled r < 0.05** (which §2 says it must be, since it is upper-bounded by
   the +0.07 that is already measured, times the attenuation of `g`).

If it lands materially above 0.07 the data-processing argument in §2 is empirically wrong and this kill
should be revisited. It will not.

---

## 8. Verification ledger for this document

**Verified this session (WebFetch):** HuggingFace `tahoebio/Tahoe-100M` — 1.69 TB usedStorage,
338 GB download, 3,389 files, 50 cell lines, 1,100 compounds, CC0-1.0; its `metadata/` and
`pseudobulk_differential_expression/` file listing and sizes. HuggingFace `MahmoodLab/hest` — 1,276
ST profiles aligned to WSI at <1.15 µm/px, 26 organs, 398 cancer samples / 25 cancer types, 77,766
downloads. PubMed esearch `"patient-derived organoid" AND "CRISPR screen" AND cancer` → count 1
(PMID 40838977). PubMed esummary PMID 40838977 → Venkadakrishnan et al., *Cancer Research* 2025.
PubMed esearch PDX+deep learning+histology+drug response → count 2 (PMIDs 40595000, 36960342);
esummary confirms titles/venues/years.

**Recomputed from disk this session:** all of §0 (cohort sizes, per-cancer-type counts, Replogle matrix
shape, DepMap matrix shape, DepMap per-lineage cell-line counts, free disk).

**Carried from the scout doc without re-verification** (each already carried a PMID/DOI there):
DeepDEP PMID 34417181; Ahlmann-Eltze PMID 40759747; Coladan PMID 42449400; Replogle PMID 35688146.

**COULD-NOT-VERIFY (this session):** Coladan-human3K release channel; exact provenance of the
6,192-vs-6,443 discrepancy in the paired-cohort count (both are internally consistent; likely a QC filter).

**Not verified because not needed:** nothing in this kill depends on a citation I could not resolve.
The two fatal obstacles (§1, §2) are algebraic and rest on file shapes I read off disk myself.
