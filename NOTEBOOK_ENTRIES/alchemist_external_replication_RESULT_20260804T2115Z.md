# ALCHEMIST external replication — RESULT

**UTC** 2026-08-04T21:15Z
**Predeclared in** `PREDECLARED_alchemist_external_replication_20260804T1830Z.md`
plus addenda `..._ADDENDUM_dilution_arm_20260804T1905Z.md` (the dilution curve does not price
this contaminant) and `..._ADDENDUM_cancer_labels_20260804T2015Z.md` (the cancer-label choice).
All three were committed before any channel number existed.

---

## The number that has to come first

**Cohort classifier, TCGA vs ALCHEMIST: AUC = 0.99906.**
Within-TCGA control: **0.50016** (gate G3 requires [0.45, 0.55] — passes).
12,615 patches per cohort, row-L2-normalised raw 1536-d H-Optimus-0 embeddings,
`v2/calibra/hest.py:cohort_classifier_auc` unmodified — the instrument that measured
TCGA vs HEST at 0.99999.

**ALCHEMIST separates from TCGA essentially as cleanly as HEST does.** Residual batch signal is
available to every downstream result and this line must travel with every channel number below.
What it does *not* license: the channel is measured **within** each cohort separately, never
pooled, so a between-cohort direction is not available to either measurement. What it kills
outright: any claim that this representation is cohort-invariant.

The batch structure is not mysterious. ALCHEMIST carries 10 distinct Aperio `ScanScope ID`s
(dominant one 48% of slides) and scan dates spanning 2015–2024, against TCGA's ~2010–2013.

## Verdict: **REPLICATES**, on the predeclared bar, with margin

`R = (ALCH_observed − ALCH_null_median) / (TCGA_NSCLC_observed − TCGA_NSCLC_null_median)`.
Bar was: REPLICATES at p < 0.01 and R ≥ 0.60.

| Target block | TCGA-NSCLC obs / null | ALCHEMIST obs / null | **R** | p |
|---|---|---|---:|---:|
| **59 targets** (G1 + coverage, primary) | 0.7660 / 0.2541 | 0.8201 / 0.2239 | **1.165** | 0.0033 |
| 74 targets (G1 only, no coverage cut) | 0.7717 / 0.2552 | 0.8062 / 0.2229 | **1.129** | 0.0033 |
| **59 targets, ALCHEMIST subsampled to n = 841** | 0.7660 / 0.2541 | 0.8207 / 0.2523 | **1.110** | 0.0033 |

p = 0.0033 is the 1/301 resolution floor at 300 permutations — the finest the design supports,
on every arm.

**R > 1: the channel does not attenuate, it comes out slightly larger on ALCHEMIST.** That is
the surprising direction and it is not being sold as a win. Two things about it:

- **Part of it was sample size, and that part is now removed.** In-sample top-CCA is a maximum
  over 16 directions, so it is capacity-inflated, and the inflation depends on n. ALCHEMIST's
  1,106 patients gave it a *lower* null (0.2239) than TCGA-NSCLC's 841 (0.2541), flattering R
  for a reason unrelated to either cohort. Subsampling ALCHEMIST to 841 raises its null to
  0.2523 — within 0.002 of TCGA's, confirming the mechanism — and R falls 1.165 → **1.110**.
  **1.110 is the number to quote.**
- **The residual ~11% is not explained here.** The two arms still differ in confound-design
  richness (8 design columns vs 3), expression platform (STAR-Counts TPM vs EBPlusPlus RSEM,
  neutralised in principle by within-sample ranking but not proven so), and everything the
  0.99906 cohort AUC represents. I am not attributing the excess.

## What did not move the answer

| Sensitivity | R (59 targets) |
|---|---:|
| primary: ALCHEMIST's own pooled ICD-O histology, 7 levels | 1.165 |
| declared `LUAD` / `LUSC` / `OTHER` mapping, 3 levels (822 / 209 / 75) | 1.169 |
| TCGA with `cancer + pooled TSS` instead of cancer only | excess 0.5042 vs 0.5118 |
| ALCHEMIST with `cancer + scanner` instead of cancer only | excess 0.5966 vs 0.5961 |

**The cancer-label decision is not load-bearing** — it moves R by 0.004. This was the specific
risk flagged against this experiment (a basis change moved the whole D2 effect by +0.118–0.120),
and it does not materialise here. Adding TCGA's site term costs 1.5% of its excess; adding
ALCHEMIST's scanner term changes nothing, which is worth noting given the scanner heterogeneity
above — the scanner identity is not what the channel is riding on.

## What this result is **not** allowed to say

**"Tissue-level sampling costs little."** The dilution curve cannot price this contaminant. Its
measured arm is `foreign_tumour` — other patients' same-cancer tumour — and the three
normal-tissue arms (`pooled`, `matched`, `dx_normal`) were never run; `DILUTION_LOWER_BOUND.md`
says so itself and forbids the lower-bound reading. ALCHEMIST's contamination is adjacent normal
lung and stroma from the *same* patient, the untested case. There is no measured number behind
"the missing polygons didn't matter", and this entry does not supply one. **The clearly indicated
next experiment is the `dx_normal` arm**, whose blocker was GPU re-embedding of TCGA normal
slides — machinery that is now deployed and working.

**"The channel transfers."** This is a *replication*: the same instrument, applied independently
inside a second cohort, finds a channel of comparable size. No TCGA-fitted model was run on
ALCHEMIST. Transfer is a different question and is untouched here.

## Gates

| Gate | Result |
|---|---|
| **G1** expression reimplementation reproduces the frozen block at r ≥ 0.999 | **PASS — 74/74 MSigDB targets at r = 1.000000**, max abs diff 0.0000 on 400 duplicate-free TCGA samples. All 15 curated mechanism programmes FAIL (r 0.76–0.99) and were dropped from **both** cohorts. |
| **G2** ≥ 0.95 gene coverage in ALCHEMIST | Applied in the shared 17,407-symbol universe; removes 15 further targets, leaving 59. The 74-target block is reported alongside and agrees (R 1.129 vs 1.165). |
| **G3** within-TCGA control AUC in [0.45, 0.55] | **PASS — 0.50016** |
| **G4** ≥ 95% of slides yield 30 patches, md5 match, finite blocks | **PASS — 1,106/1,106 slides, every one exactly 30 patches, 0 bad blocks, every md5 matched GDC** |

G1's stop-rule ("if fewer than 60 of the 90 survive, stop") was written against G1 and G1
returned 74. The intersection with G2 leaves 59, one below that figure; both blocks are
therefore reported and they agree.

## Cohort and cost

- **1,106 paired patients**, GDC intersection re-derived independently and matching the scout
  entry. Slides 1,349 files / 1,175 cases; STAR-Counts expression 1,138 files / 1,107 cases.
- **1.789 TB downloaded** (one slide per patient; all 1,349 slides would be 2.076 TB). ~3 h at
  ~170 MB/s across 6 parallel GDC streams. Peak disk 22 GB via download → patch → embed → delete.
- Manifest sha256 `b40c1909a323c12afef50a3358ee4af5a0cc8aa207d50f9694ba6f78b10c1773`;
  expression file-id list sha256 `be88aafe90cd5e32dab64a80acd7597625c0ac45667b7c27dcfdc4f110a9c5bd`.
- TCGA comparator: **841** of 846 LUAD+LUSC patients (5 have fewer than 30 tokens in the store).

## Protocol notes worth keeping

- **The 128 µm field of view is encoded in the TCGA source, not merely inferred.** TCGA-UT patch
  filenames end in a third field that is constant within each slide and equals `round(128/mpp)`
  native pixels: 507 → 0.2525 µm/px, 506 → 0.2530, 255 → 0.502 (a 20× slide), 1157 → 0.1106.
  ALCHEMIST slides read MPP 0.2462–0.2634 and crop to 486–518 px for the same 128 µm.
- The renderer constants were checked against real TCGA patches rather than trusted: 40 sampled
  patches are all 256×256, all JPEG 4:2:0, all carrying the IJG quality-75 quantisation table.
- H-Optimus-0's `pretrained_cfg` centre-crops to 87.5%, so the analysed field is **112 µm**, not
  128 — identically on both cohorts, so comparability holds. Stated, not corrected for.

## Discrepancies found in source documents

1. **`DILUTION_CURVE.md` does not exist** anywhere in the repo or its git history. The brief and
   `extract_normal_patches.py:9` both point at it. The numbers live in
   `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2.
2. **The dilution curve's applicability was overstated** — see the addendum; measured arm is
   `foreign_tumour`, normal-tissue arms never run.
3. **ALCHEMIST histology is *not* confined to unaudited clinical supplements.** The scout entry
   says it is; `diagnoses.primary_diagnosis` is populated on `/cases` for all 1,176. Stage and
   `tissue_source_site` genuinely are `_missing`.
4. **The persisted `cancer_only` TCGA operator cannot adjust a lung comparison at all.** Its
   provenance lists 22 design columns — 21 cancers plus `cancer_nan` — with **no `cancer_LUAD`
   and no `cancer_LUSC`**, because the TCGA test split it was fitted on holds 0 of each (all 846
   NSCLC patients are in train/val under the 11v21 holdout). TCGA's own NSCLC rows raise
   `UnseenLevelError` against it exactly as ALCHEMIST does.
5. **The frozen block did not enforce its own recorded `minimum_required_coverage: 0.95`** as a
   per-signature drop: its manifest retains `HALLMARK_ADIPOGENESIS` at 189/200 = 0.945.

## A bug found and fixed mid-run, recorded because it nearly cost 20 slides silently

Three extraction shards shared one `--output-dir` and therefore one scratch `slides/` directory.
`shutil.rmtree` at the end of the function meant the earliest-finishing shard deleted `.svs`
files the other two were still reading: 17 `Unsupported or missing image file`, 2 `Cannot read
raw tile`, 1 `No such file or directory` — **all after md5 had already verified**. Bytes that
verify and then cannot be opened are absent, not corrupt, which fixed the diagnosis. Scratch is
now per shard; the 20 were re-run to completion. Final state 1,106/1,106, 0 failures.

## Artifacts

Under `~/e0_run/ext/alchemist/` on the Lambda box (NFS-backed, survives instance stop):

- `RESULT_strict/`, `RESULT_nocoveragecut/`, `RESULT_matchedN/` — `alchemist_channel.json`,
  `cohort_control.json`, and both `*_wsi_identity.npz` in the `v2/contracts.py` artifact shape
  (`patient_ids`, `cancers`, `split`, `trained_states=['wsi_identity']`, `artifact_version=4`,
  `manifest_json`, `wsi_identity` (n, 3072)). `split` is deterministic and **nominal** — nothing
  is fitted on it, since the representation has zero fitted parameters and target scoring is
  cohort-fit-free.
- `targets/` — `{alchemist,tcga_nsclc}_targets{,_nocoveragecut}.npz` and `G1_validation.json`.
- `run/staging/` — 1,106 per-slide `.npz`, each 30 × 1536, with full per-slide provenance
  (GDC file id, md5, mpp, crop_px, Aperio header including scanner id).
- Code: `v2/research/external/` — `build_alchemist_manifest.py`, `extract_alchemist_patches.py`,
  `rank_target_scoring.py`, `build_alchemist_targets.py`, `alchemist_channel.py`.

Suite: **517 passed, 0 failures** locally. On the box 490 pass and 27 error, all 27 in
`test_p2_figures.py` on a missing `matplotlib` — an environment gap, unrelated to these changes.
