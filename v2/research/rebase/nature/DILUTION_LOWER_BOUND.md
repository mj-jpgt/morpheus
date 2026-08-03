# Dilution — what non-informative patches cost the channel

Task **D** for P1. No GPU, no new data: the `foreign_tumour` arm of
`v2/research/dilution/build_dilution_artifact.py` contaminates each patient's bag with **same-cancer,
different-patient** tumour patches drawn from the existing H-Optimus store, then pushes every level
through CALIBRA.

**All four predictions were written down before the run** (`v2/research/rebase/nature/P1_PREDECLARATION.md`,
commit 1c4b4b5). **Two of them are falsified, both in the same direction: the channel is far more robust
to contamination than I predicted.**

---

## 1. Construction

Dilution level `d` is the fraction of the **final** bag that is not the patient's own tumour, so a
patient with `T` tumour patches receives `round(d/(1−d)·T)` foreign patches. Levels are **nested** (the
patches used at 10% are a prefix of those used at 80%), so the curve is not confounded by independent
draw noise between levels. Foreign patches are drawn **donor-slide-first** — a patient takes as many
patches as it needs from one randomly assigned donor slide before moving to the next — which mimics a
whole-slide bag containing a contiguous stretch of the wrong tissue rather than an unrealistic mosaic.
Donors are same-cancer, never the patient itself.

Cohort: 6,427 patients, 238,610 tumour patches, 7,644 tumour slides; 2,766 evaluated on the held-out
`test` partition. Representation: `concat(mean, std)` over 1,536-d H-Optimus tokens — the unweighted
global skip path, **no fitted parameters** — then reduced to 256 dimensions by PCA **refit per level on
train rows only** (`v2/research/dilution/reduce_dilution.py`). The reduction is not cosmetic and is
justified twice: the permutation null re-whitens the block once per permutation and a 3,072-wide SVD a
thousand times over for seven levels does not finish; and 256 is capacity-matched to the D2 states this
is read beside. Refitting per level rather than freezing the level-0 basis is the honest choice —
a frozen basis would measure how far the diluted bags drift off the level-0 axes, which builds the
answer into the transform. Retained variance 0.879–0.923 at every level.

Instrument: identical to Track 1 — `--partition test`, 108-column cancer + pooled-TSS design, seed 42,
16 components, 20 draws, **300 permutations** (resolution 1/301 = 0.0033), levels
0.0, 0.05, 0.10, 0.20, 0.40.

## 2. The curve

| requested d | achieved d (test median) | adjusted top-CCA | held-out top-CCA | **ratio to d = 0** | **null-corrected ratio** | detection floor | attenuation | effective rank | perm p |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.000 | **0.5573** | 0.4932 | 1.000 | 1.000 | 0.20 | 1.130 | 196.2 | 0.0033 |
| 0.10 | 0.091 | 0.5571 | 0.5017 | 0.9996 | 0.999 | 0.40 | 0.985 | 194.1 | 0.0033 |
| 0.20 | 0.211 | 0.5447 | 0.5129 | 0.977 | 0.968 | 0.40 | 1.003 | 190.5 | 0.0033 |
| 0.30 | 0.302 | 0.5190 | 0.4986 | 0.931 | 0.905 | 0.40 | 1.057 | 187.5 | 0.0033 |
| 0.40 | 0.400 | 0.4774 | 0.4619 | 0.857 | 0.804 | 0.40 | 1.014 | 184.7 | 0.0033 |
| 0.60 | 0.600 | 0.3971 | 0.3680 | 0.713 | 0.607 | 0.40 | 0.855 | 176.5 | 0.0033 |
| 0.80 | 0.800 | **0.2844** | 0.1922 | **0.510** | **0.333** | 0.40 | 0.863 | 161.2 | 0.0033 |

"Null-corrected ratio" = `(observed − permutation null median) / (observed₀ − null median₀)`. The
permutation null median is 0.145–0.147 at every level — the capacity floor of a 16-component top-CCA on
2,766 patients — so the raw ratio flatters the surviving channel and the null-corrected column is the
one to quote. Every level clears its own within-cancer shuffled-pairing null at the 1/301 resolution
floor; `channel_gate_failures` is empty at every level.

## 3. The predictions, graded

| # | prediction | result | verdict |
|---|---|---|---|
| **D1** | channel declines monotonically in d | 0.5573 → 0.5571 → 0.5447 → 0.5190 → 0.4774 → 0.3971 → 0.2844 | **PASS** (flat then strictly monotone; the d = 0 → 0.09 step is within noise) |
| **D2** | `channel(d)/channel(0) ≈ (1−d)` to within ±0.15 for d ≤ 0.40 | measured 1.000 / 0.977 / 0.931 / 0.857 against predicted 0.909 / 0.789 / 0.698 / 0.600 | **FALSIFIED for d ≥ 0.21** — the decline is far slower than proportional (null-corrected: 0.999 / 0.968 / 0.905 / 0.804, still outside the band) |
| **D3** | detection floor non-decreasing in d | 0.20 → 0.40, then flat | **PASS, but the measurement is censored** — see §5 |
| **D4** | at d = 0.60 the channel retains > 0 but < 55% of its d = 0 value | 71.3% raw, **60.7% null-corrected** | **FALSIFIED** — it retains more than predicted, on either statistic |

**Both falsifications point the same way: mean-pooled bags are markedly more robust to non-informative
patches than a proportional-dilution model predicts.** At 40% contamination — two foreign patches for
every three of the patient's own — 80% of the null-corrected channel survives. At 80% contamination —
four foreign patches for every one of the patient's own — a third of it still survives, and it is still
distinguishable from a shuffled pairing at the finest resolution the run supports.

The half-loss point on the null-corrected curve is **d ≈ 0.68**: you must replace roughly two thirds of
the bag with same-cancer tumour from other patients before you lose half the channel.

## 4. The "lower bound" claim — what is measured and what is assumed

**Measured:** the cost of contaminating a bag with same-cancer, different-patient tumour patches drawn
from the same store, at matched stain, matched preparation (100% FFPE diagnostic on both sides) and
matched pipeline. There is no domain shift; the only thing added is patches carrying no information
about *this* patient.

**Assumed, and I could not test it:** that this is the *most benign* possible non-informative
contaminant, hence a floor on what any other contaminant costs. Declared in the predeclaration before
the run, and repeated here because it did not become true by being measured:

> There is a mechanism arguing the other way. Normal tissue is off-manifold in a way that is similar
> across **all** patients, so it adds a near-constant offset to a mean-pooled bag, which damages
> between-patient variation less than adding a *patient-specific random* tumour shift does. If that
> mechanism dominates, `foreign_tumour` is an **upper** bound, not a lower one.

The four-arm design in `build_dilution_artifact.py` exists precisely to settle this (`pooled`,
`matched`, `dx_normal` are the normal-tissue arms) and all three need GPU re-embedding of normal
slides, which was out of scope here. **Until those arms run, the correct phrasing is "the cost of
preparation-matched, information-free contamination", not "a lower bound".** The number should not be
quoted with the word "lower bound" attached unless this paragraph travels with it.

One thing the arm *does* establish unconditionally: whatever a normal-tissue arm eventually shows, the
component of the loss attributable to **domain shift** can be isolated, because this arm has none.

## 5. Three caveats that limit what the curve supports

1. **The detection floor is censored.** The level grid tops out at 0.40, and the floor reads 0.40 from
   d = 0.09 onward. It is therefore "≥ 0.40", not "= 0.40", and D3's monotonicity claim is supported
   only over the one step it could resolve (0.20 → 0.40). The transmission floor reads 0.05 everywhere,
   which is the *finest* level in this grid, so it is likewise censored from below. A finer grid at
   both ends would be needed to publish either floor as a function of d.
2. **The whole curve is one representation.** `raw_hoptimus_meanstd` has no fitted parameters, which is
   the right choice for isolating the effect of the patches, but a trained aggregator with attention
   could plausibly down-weight foreign patches and would be more robust still. The number is a property
   of unweighted mean pooling, not of the modality.
3. **Single seed (42) and a single draw of donor assignments.** The nesting makes the curve internally
   consistent but gives no error bar on the level-to-level differences.

## 6. Two things the curve says that were not asked for

**Effective rank badly under-reports the information loss.** Across the full sweep, effective rank
falls **196.2 → 161.2, i.e. −18%**, while the null-corrected channel falls **1.000 → 0.333, i.e. −67%**.
Rank is nearly preserved while two thirds of the cross-modal information is destroyed. This is a fourth
independent instance of the project's "effective rank does not track information content" observation,
and it is the *opposite* limb from the G2.6 in-batch collapse diagnosis (where rank stayed pinned at
16/16 while the representation collapsed to a single direction) — here rank drifts down gently while
information collapses. Both directions of the dissociation now have evidence.

**Non-specificity is invariant to dilution.** The matched random gene-set controls track the real
targets at a fitted-direction ratio of **0.815, 0.810, 0.799, 0.792, 0.798, 0.777, 0.727** across
d = 0 → 0.80. Contamination removes real and random-control signal at essentially the same rate, so the
~76–82% non-specificity found in T1.4 is a property of the readout rather than of bag quality, and
cleaning up the patches would not fix it. In the floor's own random-direction units, 0/90 controls clear
the detection floor at every level.

---

**Logged:** 2026-08-03, 01:27–03:55 UTC.

**How obtained:** Lambda A100 box `ubuntu@150.136.45.194`, workspace `~/ws_p1`, `~/venv/bin/python`,
CPU only. (1) `python -m morpheus.v2.research.dilution.build_dilution_artifact --cohort-artifact
d2_h_seed42.npz --embeddings-h5 hoptimus_patch_embeddings.h5 --metadata-parquet
hoptimus_patch_metadata.parquet --levels 0.0,0.10,0.20,0.30,0.40,0.60,0.80 --arms foreign_tumour
--seed 42`; (2) `python -m morpheus.v2.research.dilution.reduce_dilution --n-components 256`;
(3) `python -m morpheus.v2.calibra.run_calibra --partition test --levels 0.0,0.05,0.10,0.20,0.40
--n-draws 20 --n-components 16 --n-permutations 300 --seed 42 --n-jobs 10 --score-random-controls`.
Outputs under `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/dilution/`
(persistent). The first attempt at (3) ran on the un-reduced 3,072-wide states and was abandoned after
48 minutes for 2 of 7 levels; that dead end is why (2) exists.

### Technical

Seven nested dilution levels, donor-slide-first, same-cancer different-patient tumour patches from the
canonical H-Optimus store; 2,766 held-out patients; 108-column cancer + pooled-TSS design; 300
within-cancer permutations. Adjusted top-CCA falls 0.5573 → 0.2844 from d = 0 to d = 0.80; corrected for
the permutation null median of 0.145–0.147 the surviving fraction is 1.000, 0.999, 0.968, 0.905, 0.804,
0.607, 0.333. Half-loss at d ≈ 0.68. Every level clears its own shuffled-pairing null at p = 1/301.
Attenuation stays at 0.99–1.13 up to d = 0.40 and falls to 0.855/0.863 at d = 0.60/0.80. Effective rank
196.2 → 161.2. Predeclared D2 (proportional decline ±0.15) and D4 (< 55% surviving at d = 0.60) are both
falsified in the robust direction; D1 (monotonicity) passes; D3 (non-decreasing detection floor) passes
but is censored by the level grid.

### In plain terms

We took each patient's slide, deliberately mixed in patches from *other* patients with the same cancer,
and asked how much of the image-to-molecular signal survived. We wrote down beforehand that we expected
the signal to fall roughly in proportion to how much of the bag we had replaced — replace 40%, lose
about 40%.

That is not what happens. Replace 40% and you lose about 20%. Replace 80% — four foreign patches for
every one real one — and a third of the signal is still there, still clearly distinguishable from
random. Averaging over a bag of patches turns out to be far more forgiving of useless patches than we
assumed. Two of our four written-down predictions were wrong, both because we were too pessimistic.

The caveat we flagged in advance and still cannot remove: we called this a "lower bound" on the cost of
useless patches, on the assumption that patches from the same cancer are the gentlest possible kind of
contamination. We cannot check that without the normal-tissue arms, which need GPU work we did not do,
and there is a plausible argument the assumption is backwards. So the honest label is "the cost of
preparation-matched, information-free contamination", not "a lower bound".

### Meaning for the claim

* **For P1:** the instrument produces a clean, monotone dose–response on a representation with zero
  fitted parameters, with every level passing its own permutation null and its own random gene-set
  control. That is a usable demonstration that CALIBRA measures a channel and not a fitting artefact.
* **For any patch-selection or tissue-segmentation claim on this project:** the returns are much smaller
  than assumed. Removing a *modest* amount of uninformative tissue buys almost nothing — the curve is
  flat to 20% contamination. Effort spent on patch curation below the ~40% level is not recoverable in
  channel terms.
* **For the "effective rank" thread:** rank falls 18% while information falls 67%. Rank must not be used
  as a proxy for information content, and this is now the fourth independent demonstration.
* **For T1.4/T1.5:** the ~78% random-control ratio is invariant to dilution, so the non-specificity of
  the channel is a property of the readout and cannot be blamed on patch quality.
* **Not done:** the `pooled`, `matched` and `dx_normal` normal-tissue arms (GPU re-embedding required),
  which are what would turn the assumed direction of the bound into a measured one; finer level grids at
  both floor extremes; a second seed.
