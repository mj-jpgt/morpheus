# External cohort options for paired H&E WSI + molecular readout, verified against live APIs

**UTC** 2026-08-04T19:00Z
**Remit** non-TCGA, paired H&E whole-slide imaging + a molecular readout. Ranked by
effort-to-first-result.
**Standard** every number carries the endpoint that produced it. Anything not confirmed is
marked COULD NOT VERIFY and left as a gap. Three fabricated citations have contaminated this
project; nothing here is estimated and presented as measured.

Supersedes the WSI+*outcome* scout `v2/research/rebase/nature/data_external_wsi_outcome.md`
(2026-07-29) for the *molecular* question. That report's CPTAC section remains broadly right
but three of its premises are now corrected below.

---

## The result that changes the plan

**ALCHEMIST is a larger, more open, lower-effort paired cohort than CPTAC, and it is already
in a portal this project knows how to use.** It did not appear in the previous scout at all.

## Ranked by effort-to-first-result

| # | Cohort | Paired patients (H&E WSI + expression) | Licence / gating | Format | Tumour polygons | Route |
|---|---|---:|---|---|---|---|
| 1 | **ALCHEMIST-ALCH** (GDC) | **1,106** | open, no dbGaP | SVS, 40x, ~0.25 um/px | no | GDC API |
| 2 | **CPTAC** (IDC + GDC) | **1,580** (1,554 also proteomics) | CC BY 4.0, open | DICOM VL WSM, 20x, ~0.494 um/px | **no** | IDC v3 API + s5cmd |
| 3 | CGCI-BLGSP | 252 | open | SVS, 40x (MPP 0.2527) | no | GDC API |
| 4 | HCMI-CMDC | 378 (351 tissue-derived ceiling) | open | SVS, mostly **frozen** | no | GDC API |
| 5 | CGCI-HTMCP-CC | 123 | open | SVS, 40x | no | GDC API |
| 6 | HTAN | **314** | good half **dbGaP-gated** | SVS plurality, ~1/3 non-pyramidal | no | Synapse/CRDC |
| 7 | CDDP_EAGLE-1 / CGCI-HTMCP-DLBCL / -LC / CCG-CUPP | 47 / 34 / 26 / 29 | open | SVS | no | GDC API |
| — | GTEx | 13,133 image-to-RNA-seq sample links, 937 donors | bespoke ToU | Aperio | no | **wrong-target control only** |

Non-TCGA GDC total with both slides and expression: **1,995 cases**. TCGA reference: 10,361.

## Verified in detail

### 1. ALCHEMIST-ALCH — 1,106 paired, fully open (top pick)

Verified twice, independently, by me and by a second agent, and by two different join keys
(`case_id` and `submitter_id`), all giving 1,106.

- `api.gdc.cancer.gov/cases`, program ALCHEMIST: slides **1,175** cases, expression **1,107**
  cases, **intersection 1,106**.
- `api.gdc.cancer.gov/files` facets, slide images: total **1,349**, `access` = `[('open',
  1349)]`, `data_format` = `[('svs', 1349)]`.
- Expression: total **1,138**, `access` = `[('open', 1138)]`, workflow `STAR - Counts`.
- Preservation **100% FFPE**, `tissue_type = tumor` for all 1,349 — this matters, because our
  TCGA patches are FFPE tumour and HCMI's are mostly frozen.
- Resolution read from actual SVS TIFF headers: `AppMag = 40, MPP = 0.2476` (78368x83007),
  and two further slides at MPP 0.2525 / 0.2476. Consistent 40x, ~0.25 um/px.
- Open access proven empirically, not from metadata: unauthenticated ranged GET on
  `api.gdc.cancer.gov/data/850506b2-abbe-4dfb-969f-c7b393c653b6` returned **HTTP 206**,
  `Content-Range: 0-199/1408981844`, magic bytes `II*\0` (little-endian TIFF).
- Cohort: resected stage IB-IIIA NSCLC, screening trial NCT02194738.
- Caveat: GDC reports `primary_site`/`disease_type` as `Not Reported`; histology and stage
  must come from the 1,176 clinical supplement files, which were **not** audited.

**Why this is rank 1 and not rank 2 despite CPTAC being larger:** same API, same auth (none),
same file format as the TCGA slides already in the pipeline, one disease, FFPE tumour. The
effort delta to first result is close to zero.

### 2. CPTAC — 1,580 paired, and the current plan's numbers need correcting

Pinned to **IDC data version 24.0** (`/v2/versions` -> `{"active":true,
"date_active":"2026-04-27","idc_data_version":"24.0"}`).

- **"GDC has no CPTAC slides" — CONFIRMED, verified by me directly.**
  `api.gdc.cancer.gov/files` with program CPTAC + `data_type=Slide Image` -> `total: 0`.
  Same for `data_format=SVS` -> 0. Facet over all 109,259 CPTAC files returns 19 `data_type`
  buckets and Slide Image is not among them. Sanity check: the identical query shape returns
  30,326 for TCGA.
- IDC SM-filtered counts (`POST /v3/cohort/counts`, `Modality: SM`): **2,107 unique slide
  patients, 7,377 series, 2.569 TB** across 13 collections.
- **1,580 patients have an H&E WSI in IDC and open bulk RNA-Seq in GDC**; **1,554** also have
  PDC proteomics. Computed by literal intersection of the three patient-ID sets, not estimated.
- All 7,377 series are H&E, verified **exhaustively** rather than sampled: the complete
  manifest's `SeriesDescription` values are all `HE *` (3,770 `HE tumor`, 1,643 `HE normal`,
  1,230 `HE`, ...).
- Format: 100% of IDC SM is `1.2.840.10008.5.1.4.1.1.77.1.6` (VL Whole Slide Microscopy
  Image Storage). Header read from a real slide: 20x, PixelSpacing 0.0004942 mm =
  **0.494 um/px**, JPEG baseline tiles, 43823x36780.
- Licence **CC BY 4.0** for all 7,377 series. No dbGaP, no login — proven by an anonymous
  successful S3 GET of a 66.9 MB `.dcm` from `idc-open-data`.

**Three corrections to the plan as briefed:**

| Briefed | Verified |
|---|---|
| ~1,286 paired cases | **1,580** with RNA-Seq (1,554 with RNA-Seq + proteomics) |
| ~200 GB | **2.569 TB** for all 13 CPTAC SM collections |
| CC BY 3.0 | **CC BY 4.0** (CC BY 3.0 exists in IDC but covers no CPTAC pathology) |

Also: the IDC **v1 API is dead** (HTTP 410, *"IDC v1 API endpoints have been deprecated"*).
Anything in the project citing `/v1/collections` must be re-derived from v2/v3.

Two traps that would produce wrong numbers:
- **`cptac_stad` has no pathology at all** (`Modality=SM` -> 0 patients; `image_types` is
  `CT, US`). Do not count it.
- `subject_count` in `/v2/collections` is all-modality, not SM (cptac_ucec 254 vs 250 SM).
- `cptac_cm` (95) and `cptac_sar` (88) have slides but **zero** expression in GDC or PDC.

**The annotation problem is confirmed and is worse than "no polygons exist" — it is a
fabrication vector.** IDC exposes `analysis_result_id` values literally named
`cptac_ccrcc_tumor_annotations` (636), `cptac_ucec_tumor_annotations` (617),
`cptac_pda_tumor_annotations` (536). These are **`Modality: RTSTRUCT` radiology contours on
CT/MR/PT**, with SeriesDescriptions like `'Pre-Dose, RIGHT KIDNEY - 1'`. TCIA's own page
lists their source images as `CT, MR, PT`. **If anyone cites these as WSI tumour
annotations, that is exactly the kind of citation that has already contaminated this
project three times.**

Genuine WSI annotations in IDC (`Modality=ANN`, Microscopy Bulk Simple Annotations, 7,102
series) exist in exactly 15 of 176 collections: 14 `tcga_*` plus
`bonemarrowwsi_pediatricleukemia`. **Zero CPTAC.** Cross-checked externally: Zenodo
`"CPTAC" AND "tumor annotation"` -> 0; GitHub `CPTAC+whole+slide` -> `total_count: 0`;
HuggingFace `search=CPTAC` -> 11 datasets, none annotations.

### 3. HTAN — the honest answer to "does it reach usable scale?" is **no**

This was the genuine open question. Computed from
`humantumoratlas.org/processed_syn_data.json.gz` (221,281 file records, 14 atlases, Release
7.0), cross-validated byte-for-byte against the portal's live ClickHouse backend
(`SELECT count() FROM files` = 221281).

| Cohort definition | Participants |
|---|---:|
| any file with `assayName == "H&E"` | 1,138 |
| H&E ∩ any RNA | 724 |
| H&E restricted to >=4 mm slide-scale | 537 |
| **WSI-scale H&E ∩ any RNA** | **314** |

The collapse from 1,138 to 537 has a single decisive cause: **HTAN Duke's 694 largest-N
"H&E" files are TMA cores, not slides** — all exactly 11264x11264 px at 0.161 um/px =
**1.81 mm square**, `Pyramid: "No"`.

**And the gating kills most of what survives.** All 2,474 ImagingLevel2 H&E files are open,
but bulk RNA-seq Levels 1 and 2 are `downloadSource: "dbGaP"` and only Level 3 is on
Synapse. Duke has 1,064 L1 + 953 L2 versus **8** L3 files — so Duke's 216-patient pairing,
the largest single block, is almost entirely behind dbGaP. BU (525 L3 files) is the
practical one.

H&E is **~9%** of HTAN imaging by file count (2,502 H&E against 110,398 electron microscopy,
21,156 ExSEQ, 18,626 MERFISH, 3,852 CyCIF, 2,921 Visium, 2,839 CODEX) — so the "multiplex
rather than H&E" suspicion was right. Five atlases have zero H&E: MSK, DFCI, CHOP, OHSU,
TNP-TMA. Format is SVS plurality (1,411), and **~1/3 of H&E files are non-pyramidal**, which
will break tiled WSI readers. Tumour annotations: **zero `.geojson` across all 221,281
records**; specimen-level `PercentTumorCells` exists, per-slide region masks do not.

COULD NOT VERIFY for HTAN: dbGaP accession, licence terms, whether Release 7.0 is current.
The >=4 mm WSI cut is a derived threshold, not an HTAN field, and 105 H&E files lack the
metadata to classify.

### 4. Non-TCGA GDC — a decisive negative

**Only six programs in all of GDC have any slide images**, faceted over all 35,940:
TCGA 30,326, CGCI 3,115, ALCHEMIST 1,349, HCMI 810, CCG 291, CDDP_EAGLE 49. That is the
complete list. Verified by me directly.

**TARGET, BEATAML1.0-COHORT and MMRF-COMMPASS have zero slide image files** — so the
smear-versus-tissue question about them is moot, there is nothing to evaluate. Same for
ORGANOID-PANCREATIC, EXCEPTIONAL_RESPONDERS, REBC-THYR, MP2PRT, APOLLO-LUAD, WCDT-MCRPC,
FM-AD, TRIO-CRU, CTSP-DLBCL1, NCICCR-DLBCL, CMI-MBC.

Method note worth keeping: the obvious query
`and(files.data_type="Slide Image", files.data_type="Gene Expression Quantification")` on
`/cases` returns **0**, because GDC applies file filters per file and no single file is both
types. Case-ID lists must be pulled separately and intersected. This is a trap that would
make a real cohort look empty.

HCMI caveat: of its 1,153 expression files, **666 are `next generation cancer model`**
(organoid/cell line), not tissue. Restricting to tissue-derived gives ~351 as a ceiling, and
its slides are mostly frozen rather than FFPE.

### 5. GTEx — counted, and correctly excluded from the total

- **25,713 histology images / 980 donors / 41 tissues**
  (`gtexportal.org/api/v2/histology/image` -> `totalNumberOfItems: 25713`).
- **13,133 images carry a sampleId that is literally a GTEx v10 RNA-seq sample** (937 donors)
  — same-aliquot linkage, not inferred.
- **Normal tissue, so it counts as a wrong-target control and nothing else.** It is not
  added to any paired total above. Two caveats if used even for that: only 617 of 25,713
  images are annotated `no_abnormalities` (fibrosis 1,524, congestion 1,412, steatosis 263),
  and the tissue is postmortem with pathologist-scored autolysis — so it is morphologically
  distinguishable from surgical tissue for reasons that have nothing to do with the target.
- The BRD image route in the brief is **dead**: `brd.nci.nih.gov/brd/image-search/searchhome`
  -> 301 -> `brd.cancer.gov`, and the path 404s on the new host. Bulk access is now AnVIL
  `AnVIL_GTEx_public_data` on Terra, ~8 TB, **requester-pays**.

### 6. Others checked

- **KPMP** — real and better than expected: live GraphQL at `atlas.kpmp.org/graphql` reports
  **WSI 625 participants / 5,321 files**, spatial transcriptomics 203 participants, and H&E
  specifically 1,268 images (each biopsy serially stained H&E/PAS/Trichrome/Jones). But bulk
  RNA is only 163 files, and it is kidney-only. Its public webpage counts are stale versus
  its own API.
- **HEST-1k** — already in this project. Now **1,276** samples (HF card / 1,276 `st/*.h5ad`),
  not the paper's 1,229; do not mix the two number sets. **`cc-by-nc-sa-4.0` and gated**
  (unauthenticated fetch -> HTTP 401). >2 TB.
- **STimage-1K4M** — MIT licence, `"gated": false`, 1,149 slides, but only **672 human**
  (419 mouse). Spot-level, not bulk. Overlaps HEST heavily; deduplicate by source accession.
- **AURORA** (TCIA) — 55 patients, images CC BY 4.0, but expression is in dbGaP/GEO.
- **Verified negatives, no expression data:** CAMELYON17 (CC0, labels only), PANDA, PAIP 2019
  (*"Additional clinicopathologic data are not provided"*), ANHIR, BCNB (IHC receptor status
  only, and `.jpg` not pyramidal WSI), Human Protein Atlas (IHC not H&E, and the 45-tissue
  IHC panel and 51-tissue RNA panel are **different tissue sets**, so no per-image expression
  vector exists), UK Biobank (~1,100 participants, *"not yet available for general
  research"*), Genomics England (`/digital-pathology` and `/imaging/` both 404 — put no
  number on it).
- **TCIA's API cannot enumerate pathology at all**: `getModalityValues` returns 21 modalities
  and **no `SM`**; `getSeries?Modality=SM` returns zero bytes. TCIA WSI are non-DICOM
  supplementary downloads. Use IDC for CPTAC, not TCIA.

---

## Does CPTAC remain the best option? Plainly: no, not for first result.

**ALCHEMIST should go first.** 1,106 paired open-access FFPE tumour cases in SVS at 40x,
through the GDC API this project already uses. CPTAC is larger (1,580) and adds proteomics
(1,554), but costs a new API, a new file format (DICOM VL WSM rather than SVS), and 2.569 TB.

**The tumour-polygon problem is real, unfixable, and identical for both.** No tumour region
annotations exist for CPTAC *or* ALCHEMIST *or* any non-TCGA GDC project. Our TCGA patches
are tumour-restricted; sampling either cohort would be tissue-level. That is a **declared
deviation, not a fixable one** — the alternatives are all worse:

- generating polygons with a segmentation model introduces a second model into the
  comparison and its errors become part of the result;
- comparing tissue-level external patches to tumour-restricted TCGA patches confounds
  cohort with sampling region;
- the only clean fix is to *also* recompute the TCGA side tissue-level, which changes the
  published numbers.

The honest option is the third one, run as a matched pair: tissue-level on both sides, with
the tumour-restricted TCGA number reported beside it as the deviation it is.

**And the ceiling on all of this is not the cohort.** This project's own measurement is that
a classifier separates TCGA from HEST at AUC 0.99999
(`spatial_baselines_20260803T0620Z.md`). Today's held-out-site result
(`leave_sites_out_result_20260804T1830Z.md`) shows the channel survives site shift *within*
TCGA at ratio 1.010 — but TCGA-to-elsewhere shift is a different order of magnitude, and
none of these cohorts removes that.

## Not verified — do not cite as established

Per-study PDC proteomics case counts (GraphQL timed out at 270 s); same-tissue-block
correspondence between CPTAC slides and sequenced aliquots (patient-level match only, so
1,580 is a **ceiling**); MPP across the full CPTAC slide set (one header read, not a survey);
why `cptac_lscc` has 212 slide patients but 108 with RNA-Seq; GTEx slide-to-sample
intersection at tissue level; HTAN dbGaP accession and licence; `CPTAC-Glioblastoma-CODEX`
contents; ALCHEMIST histology/stage labels; HCMI magnification. IMPRESS, Post-NAT-BRCA and
EGA-hosted cohorts were **not investigated at all** and no claim is made about them.
