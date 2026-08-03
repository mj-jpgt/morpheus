## 2026-08-03 03:55 UTC — Dilution: two of four predeclared predictions FALSIFIED, both because I was too pessimistic. Replacing 40% of a bag with other patients' tumour costs only 20% of the channel

**Logged:** 2026-08-03 03:55 UTC. **How obtained:** `build_dilution_artifact --arms foreign_tumour --levels 0.0,0.10,0.20,0.30,0.40,0.60,0.80 --seed 42` → `reduce_dilution --n-components 256` → `run_calibra --partition test --levels 0.0,0.05,0.10,0.20,0.40 --n-draws 20 --n-components 16 --n-permutations 300 --seed 42 --n-jobs 10 --score-random-controls`, Lambda box `~/ws_p1`, CPU only.

### Technical

Seven **nested** dilution levels (the patches used at 10% are a prefix of those used at 80%), foreign
patches drawn **donor-slide-first** from same-cancer, different-patient tumour slides in the canonical
H-Optimus store. 6,427 patients, 238,610 tumour patches, 7,644 slides; 2,766 evaluated on `test`.
Representation `concat(mean, std)` — zero fitted parameters — reduced to 256 dims by PCA refit per level
on train rows only (retained variance 0.879–0.923). Identical instrument to Track 1: 108-column
cancer + pooled-TSS design, seed 42, 300 within-cancer permutations (resolution 1/301).

| achieved d | adjusted top-CCA | held-out | ratio to d=0 | **null-corrected ratio** | detection floor | attenuation | eff. rank |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.000 | 0.5573 | 0.4932 | 1.000 | **1.000** | 0.20 | 1.130 | 196.2 |
| 0.091 | 0.5571 | 0.5017 | 0.9996 | **0.999** | 0.40 | 0.985 | 194.1 |
| 0.211 | 0.5447 | 0.5129 | 0.977 | **0.968** | 0.40 | 1.003 | 190.5 |
| 0.302 | 0.5190 | 0.4986 | 0.931 | **0.905** | 0.40 | 1.057 | 187.5 |
| 0.400 | 0.4774 | 0.4619 | 0.857 | **0.804** | 0.40 | 1.014 | 184.7 |
| 0.600 | 0.3971 | 0.3680 | 0.713 | **0.607** | 0.40 | 0.855 | 176.5 |
| 0.800 | 0.2844 | 0.1922 | 0.510 | **0.333** | 0.40 | 0.863 | 161.2 |

Null-corrected ratio = `(observed − null median) / (observed₀ − null median₀)`; the permutation null
median is 0.145–0.147 at every level (the 16-component capacity floor at n = 2,766), so the raw ratio
flatters the surviving channel. Every level clears its own shuffled-pairing null at p = 1/301;
`channel_gate_failures` empty at every level.

**Grading the predeclaration (`P1_PREDECLARATION.md`, commit 1c4b4b5):**

* **D1 (monotone decline) — PASS.** Flat over the first step (0.5573 → 0.5571, within noise) then
  strictly monotone.
* **D2 (`channel(d)/channel(0) ≈ 1−d` within ±0.15 for d ≤ 0.40) — FALSIFIED for d ≥ 0.21.** Measured
  1.000 / 0.977 / 0.931 / 0.857 against predicted 0.909 / 0.789 / 0.698 / 0.600. Null-corrected
  (0.999 / 0.968 / 0.905 / 0.804) is still outside the band from d = 0.21 on.
* **D3 (detection floor non-decreasing) — PASS but CENSORED.** 0.20 → 0.40 then flat, and 0.40 is the
  top of the level grid, so the floor is "≥ 0.40" and only one step was resolvable. The transmission
  floor reads 0.05 everywhere, which is the *finest* level in this grid, so it is censored from below
  too. Neither floor can be published as a function of d without a finer grid at both ends.
* **D4 (< 55% of the channel surviving at d = 0.60) — FALSIFIED.** 71.3% raw, **60.7% null-corrected**.

Half-loss on the null-corrected curve is at **d ≈ 0.68**.

**The "lower bound" label is not claimed, and the predeclaration says why.** What is measured is the
cost of *preparation-matched, information-free* contamination (100% FFPE diagnostic on both sides, same
store, same pipeline, no domain shift). Whether it is a *floor* on the cost of any contaminant depends
on an assumption I could not test: that same-cancer tumour is the most benign contaminant. The
counter-mechanism, written down before the run, is that normal tissue is off-manifold in a way similar
across *all* patients, so it adds a near-constant offset to a mean-pooled bag and may damage
between-patient variation **less** than a patient-specific random tumour shift does — in which case
this arm is an upper bound. Settling it needs the `pooled`/`matched`/`dx_normal` arms, which require
GPU re-embedding of normal slides and were out of scope.

**Two findings that were not asked for.**

1. **Effective rank badly under-reports the information loss.** Rank falls 196.2 → 161.2 (**−18%**)
   while the null-corrected channel falls 1.000 → 0.333 (**−67%**). This is a fourth independent
   instance of "effective rank does not track information content", and it is the *opposite* limb from
   today's G2.6 in-batch collapse diagnosis (rank pinned at 16/16 while the representation collapsed to
   cosine 0.9999). Both directions of the dissociation now have evidence.
2. **Non-specificity is invariant to dilution.** The matched random gene-set controls track the real
   targets at fitted-direction ratios of 0.815, 0.810, 0.799, 0.792, 0.798, 0.777, 0.727 across
   d = 0 → 0.80. Contamination removes real and control signal at the same rate, so the ~76–82%
   non-specificity found in T1.4 is a property of the readout, not of bag quality. In the floor's own
   random-direction units, 0/90 controls clear the floor at every level.

**Dead end, recorded.** The first attempt ran CALIBRA on the un-reduced 3,072-wide states and was
abandoned after 48 minutes having completed 2 of 7 levels: the permutation null re-whitens the
representation block once per permutation, and a 2,766 × 3,072 SVD a thousand times over, seven times,
does not finish. Diagnosis showed most pool workers idle — the cost is in the serial whitening, not in
the parallel draws. `v2/research/dilution/reduce_dilution.py` (PCA-256, refit per level on train rows
only) exists because of that, and it is also the fairer comparison: 256 is capacity-matched to the D2
states the curve is read beside. Refitting per level rather than freezing the level-0 basis is
deliberate — a frozen basis would measure drift off the level-0 axes and build the answer into the
transform.

### In plain terms

We took each patient's slide, deliberately mixed in patches taken from *other* patients with the same
cancer, and measured how much of the image-to-molecular signal survived. We had written down in advance
that we expected the signal to fall roughly in step with how much of the bag we had replaced: swap out
40%, lose about 40%.

That is not what happens. Swap out 40% and you lose about 20%. Swap out 80% — four foreign patches for
every real one — and a third of the signal is still there and still clearly distinguishable from
randomly paired data. You have to replace about two thirds of the bag before you lose half the signal.
Averaging over a bag of patches is far more forgiving of useless patches than we assumed, and two of our
four written-down predictions were wrong for that reason.

We also cannot honestly call this a "lower bound" yet, which is what we set out to measure. That label
rests on same-cancer tumour being the gentlest possible contaminant, and there is a reasonable argument
that normal tissue would be gentler still because it looks alike for everybody. Testing it needs GPU
work we did not do.

### Meaning for the claim

* **P1:** a clean monotone dose–response on a representation with zero fitted parameters, every level
  passing its own permutation null and its own random gene-set control. A usable demonstration that
  CALIBRA is measuring a channel rather than a fitting artefact.
* **Any patch-curation or tissue-segmentation claim on this project must be re-scoped.** The curve is
  flat to 20% contamination and has lost only 20% at 40%. Effort spent removing modest amounts of
  uninformative tissue is not recoverable in channel terms.
* **The rank thread gains a fourth data point**, and the *opposite* limb from the collapse diagnosis.
* **T1.4/T1.5 are reinforced:** the channel's non-specificity survives dilution unchanged, so it cannot
  be blamed on patch quality and must be a property of the readout.
* **Not done:** the three normal-tissue arms (GPU), which are what would turn the assumed direction of
  the bound into a measured one; finer level grids at both floor extremes; a second seed.

### Files / commits

`v2/research/dilution/reduce_dilution.py` (new), `v2/research/dilution/build_dilution_artifact.py`
(pre-existing, `--normal-staging` already optional).
Write-up: `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md`.
Data: `p1_evidence/dilution/{dilution_foreign_tumour.npz,dilution_foreign_tumour_pca256.npz,calibra_pca256/}`
under `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/`.
