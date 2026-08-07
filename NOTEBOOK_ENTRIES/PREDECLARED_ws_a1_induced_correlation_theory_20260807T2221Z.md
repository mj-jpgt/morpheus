## PREDECLARED — WS-A1: formalizing the induced-correlation predictor (complaint #3)

**Logged:** 2026-08-07 22:21 UTC, by agent W1-A. **Status: predeclared, not yet measured.** This file is
committed before `v2/calibra/induced_correlation_theory.py` is run against any real artifact.

---

## 0. What is being built and why

`paper/P1_REVISION_SPEC.md` §4 derives two regimes for the induced correlation a shared confound
design manufactures between two orthogonal signals under cross-fitted residualisation:

- **Regime I (structureless design):** `RMS|r_induced| ≈ sqrt(R) / (n − R)`, a classical
  degrees-of-freedom term that vanishes as `n → ∞`.
- **Regime II (a design that genuinely explains both modalities):** `RMS|r_induced| →
  ρ_u ρ_a E|cosθ| / sqrt((1−ρ_u²)(1−ρ_a²))`, a population constant, invariant in `n`.

The deliverable is `v2/calibra/induced_correlation_theory.py`:
- `isotropic_induced_correlation(n, design_rank)` — the Regime I law, `sqrt(R)/(n−R)`.
- `plugin_induced_correlation(R_s, R_a, cos_theta)` — the zero-free-parameter Regime II identity,
  equation (1) of the spec, evaluated on supplied (not fitted) `R_s, R_a, cos_theta`.
- `decompose_induced_correlation(x, y, design)` — imports the existing
  `induced_correlation_sweep.closed_form_induced_correlation` (PROJECT_GUIDE §2.5: do not re-derive
  the identity) and reports `(R_s, R_a, cos_theta, r_predicted, r_isotropic, excess_ratio)`.

This spec's own §4.4 is a **synthetic** simulation (labelled as such there) verifying the two regimes
against a one-hot block-design fixture, not the real 270-cell TCGA sweep. My job is to re-run the
falsifiers against the REAL sweep (`track2/main_rows.csv`, `floors_rows.csv`, `knobs_rows.csv`,
`d2i_rows.csv` on the training box, now reachable), not to re-derive or re-trust the synthetic table.

## 1. Falsifiers, verbatim from spec §4.6

- **F1.** On `permuted` and `gaussian` designs, measured `|r_induced|` must agree with `sqrt(R)/(n−R)`
  within a factor of 1.5 across the whole n-ladder. If it does not, Regime I is wrong.
- **F2.** On real designs, `R_s²` and `R_a²` must be approximately constant in `n` (regression slope of
  `R_s²` on `log n` not distinguishable from zero). If `R_s²` declines like `R/n`, the entire account is
  wrong and P3 stays post hoc.
- **F3.** The plug-in predictor using **measured** `R_s`, `R_a`, `cos θ` and no fitted constant must beat
  P3 out of sample on the rank ladder — the ladder P3 was written after seeing. If it does not, report
  that plainly; it is still an improvement in *interpretation* but not in *prediction*, and must be
  labelled that way.
- **F4.** `cos θ` on real designs must fall with rank more slowly than `1/sqrt(R)` wherever the paper
  reports rank-invariance. If `cos θ` tracks `1/sqrt(R)` exactly on real designs too, then §4.5's
  explanation of the invariance is wrong.

## 2. What "measured" means here, precisely, and a known wrinkle flagged in advance

`floors_rows.csv`/`knobs_rows.csv`/`main_rows.csv` (real, already-computed, on `main`/`floors`/`knobs`
tags from `run_track2.sh`) carry, per cell: `design_rank` (the ACTUAL matrix rank, via
`np.linalg.matrix_rank` inside `design_participation_ratio`, not the raw column count), `k_eff`,
`k_eff_shared`, `design_r2_x_median`/`design_r2_y_median` (= `R_s²`, `R_a²`), `fitted_part_cos_rms`
(= RMS of `cos θ` over the 40 planted draws), `induced_correlation_median` (measured `|r|`), and the
existing `predicted_induced_correlation_p1/p2/p3`. `R_s = sqrt(design_r2_x_median)`,
`R_a = sqrt(design_r2_y_median)`. I test F1–F4 directly against these real fields — no rerun of the
sweep is needed or planned.

**Wrinkle, flagged before running anything:** `v2/research/rebase/nature/TRACK2_INDUCED_CORRELATION.md`
§6 already reports, on this exact real data, "Measured / classical across 24 structureless cells:
median 0.379 (p10 0.238, p90 0.699)" against **`0.6745/√(n−R)`** — a DIFFERENT classical formula from
the spec's own `sqrt(R)/(n−R)` (they differ by a factor of `sqrt(R)*sqrt(n-R)/(n-R)` relative to each
other's normalisation and are not interchangeable; the first is the median of |correlation between two
INDEPENDENT residualised vectors|, the second is this spec's derived Regime-I RMS induced-correlation
bias). F1 as stated must be tested against the spec's own `sqrt(R)/(n−R)`, not against the notebook's
`0.6745/√(n−R)`. Both numbers will be reported, labelled, so the discrepancy is visible rather than
picked between silently.

## 3. What would make me distrust a FAVOURABLE result

- If F1–F4 all pass cleanly, that itself is mildly suspicious given the synthetic-vs-real distinction
  above and the pre-existing knowledge that P1/P2 already failed on this exact rank ladder — a
  suspiciously clean pass invites re-checking that `R_s`, `R_a`, `design_rank` were not silently
  swapped between real and structureless design rows, or that the F1 comparison was not accidentally
  run against `k_eff` instead of the true `design_rank`.
- If F3 passes only because the plug-in predictor is being scored on the SAME rows P3 was fitted on
  (the anchor design, n=2,530, seed 42) rather than genuinely out-of-sample cells, that is not a pass —
  it must be scored on the full ladder exactly as P3's own out-of-sample claim was scored in
  `TRACK2_INDUCED_CORRELATION.md` §5.
- If the excess ratio (`r_predicted / r_isotropic`) for real designs is reported without also reporting
  it for `permuted`/`gaussian` rows (where it must be ≈1, since Regime I *is* the isotropic law there),
  that is an incomplete comparison and must be flagged, not shipped as headline evidence alone.

## 4. Falsifier bar for F3, made concrete

"Beat P3" means: lower `log10_rms_error` (same metric `fit_scaling_law` already reports for P1/P2/P3)
of `measured / plugin_prediction` than of `measured / predicted_induced_correlation_p3`, computed over
the identical set of real-design rows in `main_rows.csv` (the ladder rows, `design_mode == "real"`).
Reported both ways; whichever is smaller wins, and a tie within 10% relative is reported as a tie, not
rounded to a win.

## 5. Scope note

This predeclaration covers WS-A1 only. WS-A2 (capacity law) and WS-A6 (exchangeability) are
predeclared separately, in sibling files committed alongside this one.
