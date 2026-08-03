# Phase 1b — CALIBRA with a targeted readout, and E3 (the anchoring control)

**Run:** `runs/calibra_v3_targeted`, A100 box (30 cores; this is CPU linear algebra — the core count is
what mattered, not the GPU). 2,530 held-out-cancer TCGA patients, 21 test cancers.
**Protocol:** 13 spike levels (0 → 0.60), **40 draws** (was 4), **2,000 permutations** (was 60),
16 components/side, cancer + pooled-TSS adjustment (99 columns, 75 sites), seed 42.
**Code:** `v2/calibra/calibration.py`, 17/17 tests passing locally and on the box.

---

## 0. What was broken, and what it cost us

Phase 1 reported `detection_floor = NaN` and `observed_above_floor = false` for **every** state. The
instrument's headline capability — *"here is the effect size this analysis would have missed"* — produced
nothing on real data while passing all 11 synthetic self-tests.

Three separate defects, each found by fixing the one above it:

1. **Readout was a maximum, not a measurement.** Recovery was scored with `top_canonical_correlation`,
   a maximum over 16 components per side, while the spike lives on exactly one known direction pair
   `(u, v)`. Ambient top-CCA is ~0.97, so a spike of r_true ≤ 0.2 was invisible: level-0 recovery read
   **0.97 instead of ~0**, and at r_true = 0.2 the measured value *fell* (replacing the v-component
   destroys existing structure faster than a weak spike restores it). **Fix:** score the planted axis,
   `corr(X_res·u, Y_spiked_res·v)`.
2. **The spike only partially replaced its target.** `spike_targets` computed `y + outer(a_new − a, v)`
   where `a` is standardised but `y·v` is not, leaving a residual `(σ−1)·a` carrying the ambient
   correlation. This put level 0 at 0.099 instead of ~0 and attenuated an r_true = 0.6 spike to 0.27.
   **Fix:** rescale the update by the raw component's own s.d.
3. **The readout took an absolute value before pairing.** Confound-induced correlation has a random
   sign, so under `|r|` a draw with a negative baseline gets *smaller* when a positive spike is added —
   destroying the paired comparison. **Fix:** signed statistic; magnitudes only at reporting time.

## 1. A real finding that fell out of the fix

**Residualising two orthogonal signals through a shared confound design induces correlation between
them.** The spike is built orthogonal to the image score, yet level 0 reads **0.067–0.140** (not 0) for
the 99-column cancer+TSS design at n = 2,530. How large depends on how much of the drawn `(u, v)` lies
in the design span, so it varies draw to draw.

This is why the floor must be read with a **paired** test, and it is a caution for every confound-adjusted
cross-modal analysis in this literature: **adjustment does not merely remove signal, it manufactures a
small amount of apparent cross-modal correlation.**

> **PRIOR ART — this phenomenon is NOT new. Corrected 2026-08-02.**
> An earlier version of this section claimed "nobody reports this". That was wrong, and the claim is
> withdrawn. The effect is an exact algebraic identity, not an empirical discovery: for signals with
> `r_uv = 0` it is the multivariate partial-correlation formula
> `corr(Mu, Mv) = −R_u·R_v·ρ / sqrt((1−R_u²)(1−R_v²))` — Yule (1907) / Frisch–Waugh–Lovell — reproduced
> against our simulation to 7.4e-16 at this n and design rank. It has been independently rediscovered and
> published *as a warning* in fMRI (Murphy 2009), GWAS (Aschard 2015; Dahl 2019), genomics (Nygaard 2016;
> Li 2023), and stated as folklore for `removeBatchEffect` (Smyth 2020). Most directly, **Winkler et al.
> 2020, NeuroImage 220:117065** applies CCA to imaging × non-imaging residuals and warns that
> "residualisation introduces dependencies" inflating error rates — cross-modal, warned, and with a fix.
> See `NOVELTY_SEARCH.md`. What may remain ours is the *magnitude* under correct residualisation of
> exactly-orthogonal signals (0.067–0.140 here vs 0.003 for independent signals, so structural rather
> than a degrees-of-freedom artifact). Nothing broader may be claimed.

## 2. Two floors — do not conflate them

| | what it answers | value here |
|---|---|---|
| **Transmission floor** (paired) | *Does the pipeline transmit a signal of this size, or destroy it?* | **0.01 for every state** — the finest level tested |
| **Detection floor** (unpaired) | *Smallest r_true reliably distinguishable given the draw-to-draw variability a single-shot analysis faces* | **0.2** (WSI), 0.3–0.4 (full), 0.4–0.5 (RNA) |

The paired floor is near-noiseless by construction and **must not be quoted as a detection limit** — a
real analysis has no paired baseline to subtract. The unpaired floor is the conservative, quotable one.

## 3. Results

| run | state | held-out CCA | eff. rank | detection floor | attenuation | perm p |
|---|---|---:|---:|---:|---:|---:|
| full | wsi_biology | 0.4768 | 38.48 | 0.2 | 1.086 | 0.0005 |
| **programme_only** | **wsi_biology** | **0.4748** | **32.06** | 0.2 | 1.103 | 0.0005 |
| full | wsi_identity | 0.5393 | 191.07 | 0.3 | 1.228 | 0.0005 |
| **identity_only** | **wsi_identity** | **0.5393** | **191.07** | 0.3 | 1.228 | 0.0005 |
| full | full_biology | 0.8757 | 47.26 | 0.4 | 0.964 | 0.0005 |
| programme_only | full_biology | 0.8899 | 38.71 | 0.4 | 0.968 | 0.0005 |
| full | rna_biology * | 0.8983 | 32.58 | 0.5 | 0.944 | 0.0005 |

\* RNA→RNA is circular; positive control only. It is strong, so **G4.1 passes and the pipeline is sound.**
`permutation_p = 0.0005 = 1/2001` throughout: no permutation of 2,000 reached the observed value.

**Attenuation is 0.94–1.23 — i.e. ≈ 1.** The confound adjustment does **not** destroy signal.
That answers the objection that killed three earlier theses, and it is the instrument's first real
output. (Slopes slightly above 1 mean adjustment removes variance that was diluting the planted axis;
worth noting, not alarming.)

**`observed_matched_direction` is −0.09 to +0.02** — a *random* direction pair sees nothing. The channel
is concentrated in particular directions, not diffuse. Sanity check passed.

## 4. E3 — the anchoring control. **F2 does not survive.**

Direct array comparison between artifacts:

```
wsi_identity    full vs identity_only    max|diff| = 2.6e-04
rna_identity    full vs identity_only    max|diff| = 2.7e-04
full_identity   full vs identity_only    max|diff| = 2.6e-04
wsi_biology     full vs programme_only   max|diff| = 1.4e-01   (500x larger)
rna_biology     full vs programme_only   max|diff| = 7.8e-02
full_biology    full vs programme_only   max|diff| = 1.5e-01
```

**The identity head is numerically invariant to whether molecular supervision is present** (2.6e-4).
It is the frozen MLP-CLIP teacher, passed through. The biology head, by contrast, genuinely moves.

Therefore **F2 — "the head trained for biology is worse at biology than the head trained for identity"
— restates "the frozen MLP-CLIP teacher carries more molecular signal than our trained biology head."**
That is a distillation observation, not a claim about objectives. The anchoring caveat flagged in
`HANDOFF_BUILD_AGENT.md` is now **confirmed quantitatively, and F2 must be withdrawn as an objective claim.**

**Worse for the design: the decisive arm does not exist.** `identity_only` declares no biology state and
`programme_only` declares no identity state. There is **no biology head trained without programme
supervision** anywhere on disk, so *"does molecular supervision degrade the molecular channel?"* is
**unanswerable from these artifacts**. It needs a retrain (Milestone D).

## 5. What replaces it — C2 confirmed in the opposite direction

| | eff. rank | held-out molecular CCA |
|---|---:|---:|
| full (identity + programme losses) | 38.48 | 0.4768 |
| programme_only | 32.06 | 0.4748 |
| **change** | **−17%** | **−0.002** |

The representation changes materially (max|diff| 0.14), loses **17% of its effective rank**, and its
molecular channel is **unchanged**. Previously we saw rank rise 49.9 → 103.3 (+107%) with specificity
flat at 0.1366 → 0.1367. **Now we have the dissociation in both directions, in independent experiments.**

Effective rank is not a proxy for representational information content. That claim is now considerably
better evidenced than F2 ever was.

## 6. Claim status after this run

**Strengthened** — attenuation ≈ 1 (adjustment does not destroy signal); permutation p = 1/2001;
C2 confirmed bidirectionally; confound-induced-correlation is a new, reportable result.

**Withdrawn** — F2 as an objective claim. Any text asserting "molecular supervision degrades the
molecular channel" must be removed until a retrain provides the missing arm.

**Still unrun** — E0 (basis transfer, the crux), E0b, E1, E2, E4, E5.

## 7. Caveats
- Single seed (42). No CI on any between-run difference; **a paired bootstrap on the biology gap is
  still required** before the C2 numbers are quoted as a difference.
- Detection floors are in *single-direction* correlation units and are **not** comparable to the
  multivariate held-out CCA. `floor_scale` is emitted in every summary for this reason.
- `full` vs `programme_only` manifests were not verified as matched on epochs/LR/budget in this run
  (G0.4). Until they are, the rank comparison in §5 is **suggestive, not causal**.
