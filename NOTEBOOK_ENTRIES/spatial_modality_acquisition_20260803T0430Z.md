## 2026-08-03 04:30 UTC — HEST-1k acquired as the second modality; it is CC BY-NC-SA, not CC BY, and 44 of 1,276 samples survive a clean-licence + protocol-match filter

**Logged:** 2026-08-03 04:30 UTC. **How obtained:** `huggingface_hub` listing/download of `MahmoodLab/hest` (metadata `HEST_v1_3_0.csv`) and `v2/calibra/hest_build.py scan`, run on the A100 box (150.136.45.194) in `~/ws_spatial`, staged to `/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/`.

### Technical

**Licence — the brief was wrong, and the error is in the unsafe direction.** HEST-1k is released under
**CC BY-NC-SA 4.0**, not CC BY: "HEST-1k, HEST-Library, and HEST-Benchmark are released under the
Attribution-NonCommercial-ShareAlike 4.0 International license." The collection is also **HF-gated**
("You need to agree to share your contact information to access this dataset"). Repo *listing* is public,
but every `resolve` returns 401 `GatedRepoError` without an approved token. The token already on the box
(user `mj0jpgg`) is approved, so nothing was blocked.

The per-sample `license` column matters more than the collection licence. Across the 386 human-cancer
samples: 242 permissive CC BY variants, **41 CC BY-NC-ND 4.0**, 29 `Internal`, 23 CC BY-NC 3.0, 2 CC BY-NC
4.0, 49 blank. The **NoDerivatives** 41 are the ones that genuinely bite — tiling and embedding is
plausibly a derivative work — and `Internal`/blank are unresolvable without contacting depositors.

**Cohort selection.** 1,276 total → 689 human → 386 human cancer → 179 human cancer **Visium** → 54 also
permissively licensed → 45 also with a usable pixel size → **44** after dropping one targeted panel
(`nb_genes` 1,056; the rest are whole-transcriptome). Visium was chosen over the 144 legacy
"Spatial Transcriptomics" samples (100–150 µm spots on a 200 µm pitch — too coarse to localise anything at
our field of view) and over Xenium/CosMx (imaging-based, subcellular, a different target semantics).

Final cohort: **44 slides, 144,205 usable spots**, 17,197 genes in the cross-slide intersection, organs
Bowel 18 / Brain 7 / Breast 6 / Prostate 6 / Skin 4 / Lung 2 / Ovary 1. Download was 88 GB in 226 s with
zero failures. Slide-grouped split: 22 train / 9 val / 13 test slides = 58,556 / 32,418 / 53,231 spots.

**Protocol match — the load-bearing decision.** Tiles are cut to the **same fixed 128 µm × 128 µm field**
as all 271,710 existing TCGA patches (crop `128/mpp` native px, BICUBIC to 256×256, JPEG q75 4:2:0, no
colour normalisation), centred on each spot centroid. The renderer is *imported* from
`v2/research/dilution/extract_normal_patches.py` rather than reimplemented, and `test_hest.py` asserts the
constants on both sides still agree, so the two stores cannot drift apart silently.

Selecting on pixel size was necessary, not cosmetic: HEST `pixel_size_um_estimated` runs 0.137–3.24 µm/px,
and at the coarse end a 128 µm field is only ~40 native px, so reaching 256×256 would be a 6× upsample of
nothing. Requiring ≥224 native px per field keeps every retained slide in **downsampling** territory
(observed mpp 0.172–0.460, crop 278–746 px) — the same regime as TCGA, never invented detail.

**What the spot→window map costs.** A Visium spot is 55 µm across on a 100 µm pitch (2,376 µm² assayed).
A 128 µm window is 16,384 µm² — **6.90× more tissue than the transcriptome came from** — and its corners
reach 90.5 µm, so it also contains tissue assayed by the six neighbouring spots. The image therefore sees
strictly more than the target does. This ratio is computed by `window_area_ratio()`, asserted in tests, and
written into every artifact manifest. `--fov-microns` makes the spot-matched (55 µm) ablation a one-flag
change, but that ablation would *break* comparability with the TCGA store, which is the entire reason this
cohort is usable — so 128 µm is the primary and the cost is declared rather than removed.

**Unit of analysis.** A row is a **spot**, not a patient: `patient_ids` holds `<sample_id>__<barcode>`.
HEST's `patient` field is unpopulated for most public 10x samples (6 distinct labels across 44 slides), so
slides are the smallest safe grouping unit. `write_spatial_artifact` **refuses** any artifact in which a
slide straddles two partitions — spot-level random splitting would put neighbouring, often physically
overlapping, tissue on both sides of the split.

### In plain terms

We now have a second kind of molecular measurement, and it is the kind that fixes a real problem. Until
now every molecular target was an average over a whole slide, while the picture we fed the model came from
a hand-picked tumour region — so we could never tell whether a result was about morphology predicting
biology, or just about that averaging mismatch. Spatial transcriptomics measures gene expression at
thousands of individual spots on the slide, so each picture can be paired with the biology of the exact
tissue in it.

Two caveats. First, the dataset is not as freely licensed as we were told: it is non-commercial and
share-alike, and a chunk of it forbids derivative works outright. We used only the cleanly-licensed
portion, which cost us most of the data — 44 slides out of 1,276. Second, the pairing is not perfect: to
keep the images comparable to our existing TCGA images we must cut them to a fixed physical size, and that
size happens to be about seven times larger than the patch of tissue each expression measurement came
from. The image sees more than the measurement does. We measured that factor exactly rather than hiding it.

### Meaning for the claim

This satisfies P4's phase-gate requirement of at least one modality beyond bulk RNA, and it is the
prerequisite for testing whether our existing results are a property of the morphology→molecular channel
or an artifact of bulk averaging. It does **not** yet test that — no scientific claim has been re-run on
spatial data, deliberately.

Two things constrain how far this cohort can be pushed. The clean-licence subset is 44 slides across 7
organs, dominated by bowel (18); that is enough to measure a channel but not enough to claim it
generalises across cancer types. And the 6.9× window/assay area mismatch means any measured
morphology→expression association is an association between a 128 µm image and the expression of the
~15% of it that was actually assayed — a dilution, and one that biases *against* finding a channel rather
than manufacturing one.

The licence finding should propagate to the paper: HEST-derived results cannot be released under a
permissive licence without checking the ShareAlike obligation, and the 41 NoDerivatives samples must stay
out of any released derivative regardless.

### Files / commits

- `v2/calibra/hest.py` — geometry, expression normalisation, slide-grouped splits, baselines, artifact/target writers
- `v2/calibra/hest_build.py` — `scan` / `tile` / `assemble` / `baselines` stages
- `v2/tests/test_hest.py` — 23 CPU-fast tests; protocol-constant guard, split-straddle guard, panel-leakage guard
- Cohort + plan: `/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/{hest_meta/selected_visium.csv,out/plan.json}`
- Raw data (persistent): `/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/hest_data/` (88 GB)
