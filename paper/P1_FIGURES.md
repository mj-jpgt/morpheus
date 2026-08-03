# P1 — figure and table plan

Companion to `paper/P1_CALIBRA_DRAFT.md`. One row per display item. For each: the **exact data**
it must be drawn from, the **single claim** it carries, and its **status** — `PLOTTABLE` (data on
disk, no new computation), `NEEDS EXTRACTION` (data exists but must be pulled from a run output that
is not yet in a plot-ready file), or `NOT MEASURED` (the figure cannot be drawn and the paper says so
in text).

Nothing in this file may be drawn from a number that is not in the cited source. If a panel needs a
value that does not exist, the row says `NOT MEASURED` and the draft states the absence in prose
rather than the figure implying it.

Box paths are relative to `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/` on
`ubuntu@150.136.45.194` (persistent NFS). Repo paths are relative to the repository root.

---

## Main figures

### F1 — The instrument's own failure

**Claim.** A test suite that a broken instrument passes is not a control; a spike pushed through the
real pipeline is.

**Panels.**
- (a) Pre-fix recovery curve on real data: recovered value against `r_true` over the 13-level grid,
  showing level-0 reading ~0.97 (ambient top-CCA) and the value *falling* at `r_true = 0.2`, with
  `detection_floor = NaN` annotated for every state.
- (b) The same curve after the three fixes, level 0 near the induced baseline and a monotone
  recovery with slope ≈ 1.
- (c) A three-row inset naming the three nested defects and the readout change each one required.

**Data.** `runs/calibra_v3_targeted` summaries (post-fix); pre-fix values as recorded in
`v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §0 and `PHASE1_RESULT.md` (all
`detection_floor = NaN`, `observed_above_floor = false`).

**Status.** Panel (b) `PLOTTABLE`. Panel (a) **`NEEDS EXTRACTION`** — the pre-fix per-level curve was
described in prose, not preserved as a level-by-level array; only the endpoints (level-0 ≈ 0.97, the
fall at 0.2, NaN floors) are recorded. If the array cannot be recovered, draw (a) as an annotated
schematic explicitly labelled "values as recorded; per-level array not retained" or drop the panel and
keep the text.

---

### F2 — The adjustment is verified, not assumed

**Claim.** The confound adjustment removes what it claims to remove, on two different confounders,
and a per-axis certificate would have missed the leak.

**Panels.**
- (a) Grouped bars, six state × artifact combinations: joint LDA balanced accuracy raw vs adjusted,
  with the chance line at 0.0118 and each bar's own joint null p95 marked. Raw 0.2348–0.3633; adjusted
  0.0052–0.0118.
- (b) The shape of the leak: per-axis balanced-accuracy distribution over 256 axes for
  `d2_h::wsi_biology`, with the per-axis null p95 line, the per-axis maximum (0.0532), the median
  (0.0334, *below* its own null p95 of 0.0392), and the joint value (0.3633) as a separate marker far
  to the right. This panel is the figure that carries the certification-rule finding.
- (c) Breaching-axis counts raw → adjusted (17→0, 60→0, 58→0, 43→0, 61→0, 48→0).
- (d) Small companion bar: cancer-type balanced accuracy 0.463 → 0.035 against chance 0.048 (different
  cohort, n = 2,530 — label it as such, do not merge with (a)).

**Data.** `p1_evidence/track1/certificate_raw/{confound_certificate.json,task_rows.csv}` and
`certificate_adjusted/` (same filenames); tabulated in
`v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md` §T1.3 and
`NOTEBOOK_ENTRIES/t13_confound_certificate_20260803T0152Z.md`. Panel (d) from
`v2/research/rebase/nature/PHASE1_RESULT.md`.

**Status.** (a), (c), (d) `PLOTTABLE`. (b) `NEEDS EXTRACTION` — needs the per-axis vector from
`task_rows.csv` on the box.

---

### F3 — Attenuation ≈ 1: the adjustment does not destroy signal

**Claim.** Confound adjustment costs the injected signal essentially nothing — the measurement that
retires the objection that killed three earlier theses.

**Panels.**
- (a) Recovered vs injected correlation, one line per Track 1 state (six lines), with the identity
  line. Slopes 0.974–1.039.
- (b) Attenuation slope as a strip plot across all measurement contexts, one point per cell: Track 1
  six states (0.974–1.039), Track 2 twelve real design × n cells (1.07–1.12), Track 2 Gaussian control
  three cells (1.000–1.001), dilution seven levels (0.855–1.130), Phase 1b seven states
  (0.944–1.228). Colour by run; **do not pool** — the point is that every context lands near 1, not
  that they are one population.

**Data.** `p1_evidence/track1/` CALIBRA summaries (`attenuation` field per state);
`p1_evidence/track2/` sweep outputs; `p1_evidence/dilution/`. Tabulated in
`TRACK1_NEGATIVE_CONTROLS.md` §T1.7(c), `TRACK2_INDUCED_CORRELATION.md` §8,
`DILUTION_LOWER_BOUND.md` §2, `PHASE1B_TARGETED_READOUT.md` §3.

**Status.** `PLOTTABLE`.

---

### F4 — The two floors, and the scale on which they are not comparable

**Claim.** The paired transmission floor and the unpaired detection floor are different quantities;
only the unpaired one is a detection limit; and neither is on the same scale as the reported channel.

**Panels.**
- (a) One draw's recovery curve with both floors marked: transmission floor at the finest grid level
  (≤ 0.01, drawn as a left-censored marker) and detection floor at 0.30, with the level-0 upper tail
  across draws shaded.
- (b) Detection floor as a strip plot across the cells where it was measured, showing it is a property
  of (representation × target block × design × n): 0.20/0.30/0.40 across five target blocks on two
  artifacts; 0.25–0.30 across twelve real design × n cells; 0.050/0.015/0.010 for the structureless
  Gaussian control.
- (c) **The scale panel.** Three markers on one axis, explicitly labelled with their units:
  detection floor 0.30 (single random direction), `observed_matched_direction` −0.028 to +0.036 (real
  channel through a *random* direction pair), held-out top-CCA 0.4703–0.6052 (16-component
  multivariate maximum, chance 0.147). The panel exists to make it visually obvious that the first and
  third may not be compared, which was the original defect of §4.1.

**Data.** `p1_evidence/track1/` CALIBRA summaries (`transmission_floor`, `detection_floor`,
`floor_scale`, `observed_matched_direction`, `observed_above_floor`);
`v2/research/rebase/nature/GATE_LOG.md` rows `T1.7c_spike_recovery::*` and `T1.2_baseline_block::*`
(per-block `detection_floor=`); `p1_evidence/track2/` for the design × n cells.

**Status.** (b), (c) `PLOTTABLE` from `GATE_LOG.md` alone. (a) `NEEDS EXTRACTION` — needs one draw's
per-level array from the run output.

---

### F5 — Induced correlation: identity, and the falsifier that gives it content

**Claim.** The mechanism is an exact identity (conceded, not claimed), and the magnitude is
structural rather than the classical degrees-of-freedom term.

**Panels.**
- (a) Measured induced correlation against the closed-form Yule/FWL prediction, one point per draw
  over all 270 sweep cells, with the identity line. Annotate max |disagreement| = 8.6 × 10⁻¹⁶, median
  3.1 × 10⁻¹⁶, per-draw Pearson r = 1.000000. This panel is the concession, drawn as a figure.
- (b) The P5 falsifier: induced |r| against n for three arms — real cancer+TSS (0.0866 → 0.0718),
  the same design row-permuted (0.0087 → 0.0020), and Gaussian k = 99 (0.0164 → 0.0016) — on a log y
  axis, with the ratio annotated at each n (9.9× → 35.4×) and the predeclared bar (≤ 0.025 at
  n = 2,530) drawn as a horizontal line.
- (c) The n ladder alone with the sampling expectation √(2530/n) overlaid, showing the measured ratio
  0.960 against a sampling expectation of 0.627.

**Data.** `p1_evidence/track2/` sweep outputs (main tag), 270 cells, seeds 42/43/44;
`TRACK2_INDUCED_CORRELATION.md` §1–§3. Predeclaration `P1_PREDECLARATION.md` §B (commit `1c4b4b5`),
grading by `p1_evidence/grade_t2.py`.

**Status.** `PLOTTABLE`.

---

### F6 — Design rank is the wrong axis

**Claim.** The induced correlation is set by how much of *both* modalities the confound design
explains, not by how many columns it has.

**Panels.**
- (a) Induced |r| against `k_eff` across the 11-design rank ladder at n = 2,530 (k_eff 0 → 451), flat
  apart from the two low-R² designs — with `cancer` (k_eff 31, |r| 0.0844) and `cancer + tss_pool10`
  (k_eff 104, |r| 0.0748) labelled to make the point that adding 478 columns *lowers* it slightly.
- (b) The same values against `R_x·R_y/√k_eff_shared`, collapsing onto a line. Include `tss_pool10`
  (k_eff 74 but |r| 0.0379) as the labelled counter-example to a rank explanation.
- (c) Predictor ratio (measured / predicted) for the four candidate forms as box plots: P1 median 1.98,
  P2 median 2.98, P3 median 0.886 (p10–p90 0.76–1.07). **Panel caption must state that P3 is post hoc
  on the rank axis** and out of sample only on n and seed.

**Data.** `p1_evidence/track2/` rank-ladder cells; `TRACK2_INDUCED_CORRELATION.md` §4, §5.

**Status.** `PLOTTABLE`.

---

### F7 — More patients do not buy sensitivity

**Claim.** The detection floor is set by the induced correlation, not by sample size — the paper's
most operationally consequential result.

**Panels.**
- (a) Detection floor against n (1,000 / 2,530 / 6,427), one line per design: four real designs flat
  at 0.25–0.30, `gaussian_k99` falling 0.050 → 0.015 → 0.010. Log y axis.
- (b) The same cells plotted as detection floor against induced correlation, showing the floor tracks
  the induced baseline across a 60× range of the latter.

**Data.** `p1_evidence/track2/` floor sweep (tag with the 0.0…0.50 level grid, 40 draws, 2 seeds);
`TRACK2_INDUCED_CORRELATION.md` §8.

**Status.** `PLOTTABLE`.

---

### F8 — The negative-control battery, including the losses

**Claim.** A battery that only ever passes proves nothing; this one caught a real leak and two
uncomfortable results.

**Panels.**
- (a) Modality-shuffled pairing: observed adjusted top-CCA against the within-cancer permutation null
  distribution, six states, with the null median (0.1465–0.1483) and p95 marked and
  `p = 1/2001` annotated with its resolution. **The null median line is the panel's real content.**
- (b) Random gene sets, two statistics side by side: (left) 90 controls in single-random-direction
  units against the detection floor, 0/90 exceedances in every state; (right) the fitted-direction
  ratio control/real, 0.759–0.819 in every state on both artifacts. The two halves answering
  differently is the point.
- (c) Gene-label shuffle: paired difference (true − shuffled) with CI95 for three shuffle seeds, all
  three CIs crossing zero, with the true value's own CI95 [0.4874, 0.5962] as a shaded band for
  contrast.

**Data.** `p1_evidence/track1/`; `GATE_LOG.md` rows `T1.6_modality_shuffled_pairing::*`,
`T1.4_random_gene_sets::*`, `T1.4_random_vs_real_fitted_direction::*`,
`T1.5_shuffled_minus_true::*`; tabulated in `TRACK1_NEGATIVE_CONTROLS.md` §T1.4–T1.6.

**Status.** (b), (c) `PLOTTABLE`. (a) `NEEDS EXTRACTION` — the full 2,000-value null distribution is
in the run output; only median and p95 are in the markdown.

---

### F9 — Positive controls, and the measured chance level

**Claim.** The instrument recovers a covariate at the strength four published papers had already
measured, and the chance level had to be measured rather than assumed.

**Panels.**
- (a) ER-status within-cancer AUROC with CI95 for the six rows measured, the pre-registered band
  [0.78, 0.92] as a shaded region, the pre-registered point estimate 0.86 as a line, and the
  **measured** null p95 (0.542–0.546) marked against a dashed 0.5 reference. The four-point gap
  between assumed and measured chance is the panel's second message.
- (b) PR status as the second anchor against band [0.70, 0.85].

**Data.** `GATE_LOG.md` rows `T1.7b_known_covariate::covariate_{er,pr}::*`; pre-registration
`p1_evidence/inputs/PREREG_known_covariate.json`; `TRACK1_NEGATIVE_CONTROLS.md` §T1.7(b).

**Status.** `PLOTTABLE`. Caption must record that BRCA is a *development* cancer, so this ran on
`--partition all` and is in-distribution.

---

### F10 — Dose–response under information-free contamination

**Claim.** The instrument produces a clean monotone dose–response on a representation with zero
fitted parameters, and mean-pooled bags are far more robust than a proportional model predicts.

**Panels.**
- (a) Null-corrected surviving channel against achieved dilution d (1.000, 0.999, 0.968, 0.905, 0.804,
  0.607, 0.333), with the predeclared proportional model `1 − d` and its ±0.15 band overlaid so the
  falsification of D2 is visible, and the half-loss point d ≈ 0.68 marked. Plot the raw ratio as a
  faint second line to show how much the null correction matters.
- (b) Twin axis: effective rank (196.2 → 161.2, −18%) against null-corrected channel (−67%) over the
  same levels — this panel does double duty for F11.
- (c) Random-control ratio against d (0.815 → 0.727), flat, establishing that non-specificity is a
  property of the readout and not of patch quality.

**Data.** `p1_evidence/dilution/` CALIBRA outputs per level; `DILUTION_LOWER_BOUND.md` §2, §6;
`NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`.

**Status.** `PLOTTABLE`. Caption must carry the "cost of preparation-matched, information-free
contamination — **not** a lower bound" phrasing, since the normal-tissue arms were not run.

---

### F11 — Effective rank does not track information content

**Claim.** Four independent instances, in both directions.

**Panel.** A single 2 × 2 panel, one quadrant per instance, each showing Δrank against Δinformation
with the units named: (1) +107% rank / specificity flat 0.1366 → 0.1367; (2) −17% rank / channel
−0.002; (3) rank pinned 16/16 while retrieval falls 0.062 → 0.000 (below chance) and within-modality
cosine rises 0.7089 → 0.9999; (4) −18% rank / −67% channel.

**Data.** (1) `v2/research/rebase/ENGINE_CLD.md`, `HANDOFF_BUILD_AGENT.md`; (2)
`PHASE1B_TARGETED_READOUT.md` §5; (3) `NOTEBOOK.md` entry 2026-08-02 01:20 UTC, source
`scratchpad/collapse_diag.py`; (4) `DILUTION_LOWER_BOUND.md` §6.

**Status.** `NEEDS EXTRACTION` for instances (1) and (3) — both are recorded as endpoint numbers in
markdown, and instance (1) comes from an earlier codebase generation with a different benchmark
statistic. The caption must say so, and must record that instance (2)'s two arms were not verified
matched on epochs/LR/step budget. Draw as a labelled comparison chart of recorded values, not as a
regression.

---

## Main tables

### T1 — The floors, side by side

Transmission (paired, ≤ 0.01, censored) vs detection (unpaired, 0.20–0.40) with the question each
answers, the units, and an explicit "not quotable as a detection limit" row for the paired floor.

**Provenance.** `TRACK1_NEGATIVE_CONTROLS.md` §T1.7(c); `GATE_LOG.md` `T1.7c_spike_recovery::*`.
**Status.** `PLOTTABLE`.

### T2 — Negative- and positive-control battery verdict table

The ten controls with direction required and verdict, exactly as §4.8, including the two that came
back the wrong way.

**Provenance.** `TRACK1_NEGATIVE_CONTROLS.md` headline table; all rows in `GATE_LOG.md`.
**Status.** `PLOTTABLE`.

### T3 — Predeclared predictions and grades

Appendix A of the draft, promoted to a main table if the venue allows: 11 predictions, 3 failed, 1
holding with a provenance caveat.

**Provenance.** `P1_PREDECLARATION.md` (commit `1c4b4b5`); grading `p1_evidence/grade_t2.py`;
results in `TRACK2_INDUCED_CORRELATION.md` and `DILUTION_LOWER_BOUND.md`.
**Status.** `PLOTTABLE`.

### T4 — Target blocks scored through one instrument

The 16-row must-beat table of §4.13, with the corrected column header. Carries the PCA loss.

**Provenance.** `TRACK1_NEGATIVE_CONTROLS.md` §T1.1/T1.2; `GATE_LOG.md` `T1.2_baseline_block::*`.
**Status.** `PLOTTABLE`. Note the header defect in the source table (prints "baseline" twice); the
column order was verified against the `heldout=` field of the ledger rows.

---

## Supplementary

### S1 — Prior-art map for the induced correlation

A one-page table: literature, reference, verbatim quote, what it covers, what it does not. Exists to
make the concession unmissable and to pre-empt "this is Yule 1907".

**Provenance.** `NOVELTY_SEARCH.md` §Q1; `TRACK2_INDUCED_CORRELATION.md` §0.
**Status.** `PLOTTABLE` (text table).

### S2 — Winkler unit audit

The five-pass full-text audit of Winkler et al. 2020: every table and figure, and the units each
reports (pcer %, FWER %, power %, p-values), establishing that no correlation-unit magnitude appears.

**Provenance.** `TRACK2_INDUCED_CORRELATION.md` §0 (full text from PMC7573815).
**Status.** `PLOTTABLE` (text table).

### S3 — Estimator sweep

The 8-cell `n_splits` × α grid of §4.6.6, with the predeclared 25% bar drawn and the failure at
α = 100 isolated. Must be labelled **FAILED**, not "robust to fold count".

**Provenance.** `TRACK2_INDUCED_CORRELATION.md` §7.
**Status.** `PLOTTABLE`.

### S4 — Cross-fitting versus the classical term

Measured / classical `0.6745/√(n − R)` over 24 structureless cells, median 0.379 (p10 0.238, p90
0.699), with `gaussian_k600` at n = 500 marked as the one regime where cross-fitting cannot help
(0.055 measured vs 0.617 classical).

**Provenance.** `TRACK2_INDUCED_CORRELATION.md` §6.
**Status.** `PLOTTABLE`.

### S5 — Induced baseline by target block and artifact

The 5-block × 2-artifact table of §4.6.7 (0.0604–0.1062), showing the magnitude replicating on a
second artifact and varying with how much of the block the design explains.

**Provenance.** `GATE_LOG.md`, `induced_baseline=` field of `T1.2_baseline_block::*` rows.
**Status.** `PLOTTABLE`.

### S6 — Supervision-target ablation

The three tables of §4.12 (untrained-40 readout, negative control, effective rank), with the
"paired differences only" reproducibility note.

**Provenance.** `D2_RESULT.md`.
**Status.** `PLOTTABLE`.

### S7 — Gate ledger summary

101 rows: 62 gates, 39 observations, 7 failed gates, with the test that plants a baseline beating us
0.95 to 0.40 and asserts the verdict does not move.

**Provenance.** `GATE_LOG.md`;
`v2/tests/test_track1_battery_ledger.py::test_a_losing_baseline_is_an_observation_and_cannot_move_the_verdict`.
**Status.** `PLOTTABLE`.

### S8 — Field-of-view correction

The centre-crop arithmetic: 256 px cut at 128 µm → `crop_pct` 0.875 → 224 px → **112 µm**;
window/assay area ratio 5.28× rather than 6.90×.

**Provenance.** `v2/calibra/hest.py:60–106`;
`v2/tests/test_hest.py::test_effective_field_accounts_for_the_encoder_centre_crop`;
`NOTEBOOK_ENTRIES/spatial_baselines_20260803T0620Z.md`.
**Status.** `PLOTTABLE`.

---

## Figures the paper does NOT have, and says so

| would-be figure | why it cannot be drawn | where the draft says so |
|---|---|---|
| Both floors on an **external cohort** | No external cohort has been through the instrument. Deliberate scope decision. | §1.3, §5.1 |
| Induced correlation at a **second design rank on a second artifact** | The `d2_i` Track 2 sweep was queued and did not complete; the rank/n ladders are one artifact, one state. The magnitude alone replicates on `d2_i` via S5. | §4.6.7, §5.9 |
| Detection floor as a **continuous function of dilution** | Censored: the dilution level grid tops out at 0.40 and the floor reads 0.40 from d = 0.09 onward. | §4.10, §5.8 |
| Transmission floor as a **function of anything** | Censored from below at the finest grid level in every cell (≤ 0.01, ≤ 0.05 in the dilution grid). | §4.5, §5.8 |
| **Normal-tissue** contamination arms | `pooled`, `matched`, `dx_normal` need GPU re-embedding of normal slides; not run. This is why "lower bound" is withdrawn. | §4.10, §5.11 |
| A **purity** rank point in the induced-correlation ladder | No TCGA purity table on either machine; an expression-derived surrogate was rejected because it is computed from the RNA targets forming Y and would manufacture the effect under study. | §5.10 |
| **Text-prior** and **cell-composition** baseline blocks | Need external resources (gene text-embedding table; deconvolution signature matrix) on neither machine. `claim_guards.composition_attribution` stays undischarged. | §4.13, §5.10 |
| The real channel **above the floor in the floor's own units** | `observed_above_floor = 0` for every state, correctly: the floor is single-direction and the channel is a multivariate maximum. No such measurement exists. | §4.5, §5.7 |
| A **seed-level reproducibility** figure for trained artifacts | Training is not seed-reproducible on this stack; only paired within-run differences are quotable. | §3.10, §5.6 |
