## 2026-08-03 06:20 UTC — The zero-parameter per-slide mean beats H-Optimus ridge on the metric the field reports (pooled r 0.571 vs 0.357), and TCGA/HEST separate at AUC 0.99999

**Logged:** 2026-08-03 06:20 UTC. **How obtained:** `hest_build tile/assemble/baselines/cohort-control` on the A100 box (150.136.45.194), workspace `~/ws_spatial`, outputs in `/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/out/`.

### Technical

**Cohort built.** All 44 selected slides tiled and embedded, **zero failures**: 144,162 spots
at 1536-d H-Optimus-0. Of 144,205 usable spots, 43 were dropped because a full 128 µm window
would have run off the slide edge (dropped rather than clipped — a clipped window is no
longer a 128 µm field). Slide-grouped split: 58,538 train / 32,407 val / 53,217 test spots
over 22 / 9 / 13 slides. Target panel: 50 most-variable genes chosen on **train spots only**,
from the 17,197-gene cross-slide intersection, log1p CPM(1e4), plus 16 `RANDOM_CONTROL__`
permuted columns. `validate_artifact` passes; no slide straddles a partition.

**Test-set baselines, 50 genes, mean Pearson r:**

| predictor | parameters | pooled r | within-slide r |
|---|---|---|---|
| per-slide mean | **0** | **0.5706** | **-0.0000** |
| global (train) mean | 0 | 0.0000 | 0.0000 |
| ridge on H-Optimus-0 | 1536×50 | 0.3565 | 0.1506 |

The per-slide mean predicts every spot by the mean expression of its own slide. It is
constant within a slide, so its within-slide correlation is 0 by construction. **It
nonetheless beats the H-Optimus ridge model by 0.21 on pooled r — the number the HEST
leaderboard and the 46+ models in the 2026 *Brief Bioinform* survey report.** Pooled
correlation over spots drawn from many slides is dominated by between-slide variation, and
slide identity alone explains that. The ridge model is the only predictor here that carries
any within-slide signal at all (0.151), which is real but roughly a third of what its own
pooled number advertises.

**Cohort classifier control.** Logistic regression on row-L2-normalised H-Optimus embeddings,
20,000 patches per cohort, 5-fold cross-validated: **TCGA vs HEST AUC = 0.99999**. A
within-TCGA split of the same size and classifier gives **0.5012**, i.e. chance. The two
cohorts are perfectly linearly separable in embedding space.

That null took two attempts and the first one was wrong in an instructive way. `h5py` fancy
indexing requires sorted indices, so I sorted the sampled row numbers — but the TCGA patch
store is ordered by patient and slide, so the two halves of the "within-TCGA" comparison were
disjoint sets of *patients*. It returned 0.903, which reads as "the classifier overfits, so
the cross-cohort AUC is meaningless" — the opposite of the truth. Shuffling after the read
fixes it.

**A protocol correction.** H-Optimus-0's timm `pretrained_cfg` is 224 px at `crop_pct` 0.875,
`crop_mode` center, so `create_transform` centre-crops our 256 px patch to 224. The encoder
therefore sees **112 µm, not the 128 µm we cut**, and the window/assay area ratio against a
55 µm Visium spot is **5.28×**, not 6.90×. This transform is resolved from the same cfg for
all 271,710 TCGA patches, so comparability is untouched — but "128 µm" overstates the
analysed field by 14%, and both numbers now ride in every artifact manifest.

**Throughput.** Two bottlenecks, neither the one expected. Serial `read_region` off NFS is
155 ms per 744×744 window (~6 spots/s, six hours). After threading the reads, py-spy caught
the main thread inside torchvision `to_tensor`/`normalize` in 5 of 6 samples with all twelve
reader threads idle — the **CPU transform**, not the GPU, was the limit. H-Optimus-0 itself
benchmarks at 103 img/s even while three D2 trainings share the A100. Moving render+transform
into the workers brought the sweep to ~50 minutes.

### In plain terms

We can now ask the question this modality was acquired for, and the first answer is a
warning about how the field measures itself.

The standard way to report "can we predict gene expression from an image?" is to pool every
spot from every slide together and correlate prediction against truth. We built a predictor
that ignores the image completely and just says "every spot on this slide has this slide's
average expression". It has no parameters and cannot possibly know anything about local
tissue. On the pooled metric it scored 0.57. A real model, using H-Optimus-0's view of the
actual image, scored 0.36. The do-nothing baseline won.

The reason is that slides differ from each other much more than spots differ within a slide,
so pooling lets a model score well by recognising which slide it is looking at. When we
remove that — correlate only within each slide — the do-nothing baseline drops to exactly
zero, as it must, and the image model retains 0.15. So the image genuinely carries local
information. It is just far less than the usual headline number implies.

Separately, a simple classifier can tell our TCGA images from these spatial images perfectly
— 99.999% — while being unable to distinguish two random halves of TCGA at all. The two datasets
are trivially distinguishable to the model, which matters for anything that mixes them.

### Meaning for the claim

Three consequences, in order of how much they should change behaviour.

**1. Pooled spot-level correlation is not evidence of a morphology→molecular channel, and we
should not report it as one.** A zero-parameter predictor beats a foundation-model ridge on
it. Any spatial result we publish must lead with the within-slide number, and should carry
the per-slide-mean row explicitly so a reader can see the gap. This baseline is absent from
the HEST leaderboard and from the survey's 46+ models; reporting it is the most useful thing
this cohort has produced so far.

**2. The channel is real but small.** Within-slide r = 0.151 for a plain ridge on frozen
embeddings is a genuine, slide-identity-free morphology→expression signal — the first such
measurement on this project that is not confounded by bulk averaging. It is a floor, not a
ceiling: no tuning, one alpha, 50 genes, and the image window covers 5.28× the tissue the
transcriptome came from, which dilutes the association rather than inventing it.

**3. AUC 0.99999 between TCGA and HEST is a hard constraint on cross-cohort work.** Residual
batch signal is fully available to any model that sees both cohorts, so no result that pools
or transfers between them can be attributed to biology without explicit correction, and the
number must be declared wherever such a transfer appears. The within-TCGA null at 0.5012
establishes that this is a property of the cohorts, not of the classifier.

Scope note: no existing scientific claim has been re-run on spatial data. That is a later
stage and depends on decisions not yet made.

### Files / commits

- `v2/calibra/hest.py`, `v2/calibra/hest_build.py`, `v2/tests/test_hest.py` (25 tests)
- Commits `b3ad33f`, `c86cc76`, `3cdeb93`, `7cc3a20`, `4d9b146`, `8110c87`, `0789d3d` on `research/rebase-vision`
- Artifacts: `/lambda/nfs/geeg/biorag3_persistent_20260711/spatial/out/artifacts/hest_spatial_{hoptimus,targets}.npz`
- Results: `.../spatial/out/{baselines.json,cohort_control.json}`; plan `.../spatial/out/plan.json`
