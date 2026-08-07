# PREDECLARED — WS-A4, mis-specification and omitted-variable bias

**Written before `v2/calibra/misspecification.py` exists and before any number has been measured.**
Agent W1-B, executing `paper/P1_REVISION_SPEC.md` §5 WS-A4 (complaints #4, #8). Read in full:
`v2/calibra/nonlinear_adjustment.py`, `nonlinear_confound_probe.py`, `residualise.py`. Confirmed per
spec §3: the interaction design (`nonlinear_adjustment.cell_design`, the saturated cancer x site
design) and the nonlinear confound probes already exist and are **not rebuilt here**. What is
genuinely absent, confirmed by grep, is any omitted-variable-bias (OVB) test: nothing in the repo
currently drops a declared covariate and re-measures drift, or plants a confound deliberately absent
from the design and measures leakage.

## 0. Compute status

**CORRECTION, same day, before any measurement was taken.** Per
`PREDECLARED_ws_a3_association_metrics_20260807T2211Z.md` §0, `ssh p1box` is now reachable (the
coordinator reprovisioned it) and real artifacts are confirmed present with real SHA-256 hashes. **This
predeclaration now covers a real-data OVB measurement**, not a synthetic-only fallback. Synthetic
counterfactuals (§3) remain mandatory and run first. Real-data protocol:

* **Cohort/targets/design base:** identical to WS-A3 §0 — `d2_h_seed42.npz`/`d2_i_seed42.npz`, test
  partition, `frozen_rna_targets.npz`'s 90 real targets, `cancer`+`tss` design via
  `residualise.confound_design`.
* **`leave_one_covariate_out`** on the real design: drop `tss` (keep `cancer` alone) and drop `cancer`
  (keep `tss` alone), plus a real *available* third covariate — `dx_year` from
  `p1_evidence/inputs/tcga_clinical_covariates.parquet`, present on the training box (confirmed) and
  already used by `induced_correlation_sweep._load_clinical_covariates` — added to the design and then
  dropped, so the OVB test has at least one real covariate genuinely absent from the paper's shipped
  design to omit.
* **`inject_unobserved_confound`**: the planted confound itself is, by construction, a synthetic
  quantity added on top of the real `x`/`y` blocks (exactly as `calibration.spike_targets` already
  plants a synthetic spike on real targets everywhere else in this project) — this is not a "synthetic
  data" fallback, it is the standard planted-signal instrument applied to real data, and is reported as
  such (not prefixed `[SYNTHETIC]`, but the planted component's construction is stated explicitly).
* **Nonlinear-encoded arm**: built from the real `dx_year` covariate (continuous), nonlinearly encoded
  (declared before running: a decade-boundary indicator, `1[year mod 10 < 2]`, chosen because it is
  visibly non-monotonic and not recoverable by any linear/additive encoding of `dx_year`), planted as a
  confound of both real blocks and tested against the additive vs. saturated-cell designs.

## 1. The three arms, and what each is FOR

1. **`leave_one_covariate_out(frame, columns, x, y, ...)`** — for each declared covariate column,
   rebuild the design with that column dropped, re-run the channel (`channel_under_adjustment` or
   `calibration.permutation_null`, reused not reimplemented) and the detection floor
   (`calibration.spike_recovery_curve`, reused), and report the drift relative to the full design. This
   is the direct OVB test the spec asks for: "what happens to the channel/floor if a covariate we
   *do* have was left out."
2. **`inject_unobserved_confound(x, y, design, strength, ...)`** — plant a confound of known strength
   correlated with both blocks (following the same construction contract as
   `calibration.spike_targets`, imported not copied) that is **deliberately never added to the
   design**, and measure how much of it leaks into the observed channel. This answers "what if we
   missed one entirely" rather than "what if we drop one we already have."
3. **Nonlinearly-encoded-confound arm** — plant a confound whose *true* functional form is nonlinear
   in an existing design column (e.g. confound = f(numeric_column) for a non-monotonic f, or an XOR-
   style interaction of two categorical levels not captured by the additive design), include the
   **linear/additive** encoding of the same raw column in the design, and measure whether the linear
   design catches it. Reuses `residualise.confound_design` for the additive arm and
   `nonlinear_adjustment.cell_design`/`saturated_cell_residuals` for the saturated-cell upper bound,
   per spec's instruction to "report against the existing saturated-cell upper bound."

## 2. Predeclared expectations and falsifiers

**A. `leave_one_covariate_out`.** For a covariate that is genuinely unrelated to either block (planted
independent noise column, included in the design for this control only), dropping it should **not**
materially move the channel or the floor (predeclared null expectation: drift within the run's own
draw-to-draw noise band, no systematic direction).
*Falsifier A1.* If dropping a covariate that is provably independent of both X and Y (verified by
construction) moves the floor or channel by more than the noise band established by re-running the
full design at a different seed, that is evidence the OVB machinery itself introduces an artifact
(e.g. a rank-count effect unrelated to confounding), and this must be reported as a defect in the
instrument before any real-covariate result is trusted.
For a covariate that IS a real confound of both blocks (by construction), dropping it should inflate
the observed channel/floor relative to the full design — this is the informative, expected-favourable
direction, and is reported as such only after A1's null check passes.

**B. `inject_unobserved_confound`.** Predeclared expectation: the leakage into the channel scales
**monotonically** with the planted confound's strength (its R² against both blocks), and at the
strongest planted level the channel with the confound omitted should read detectably higher than the
channel with the same confound included in the design (the latter recovers the pre-injection channel,
which is the built-in check that the injected confound is not simply unremovable noise).
*Falsifier B1.* If leakage does **not** increase with planted strength (non-monotonic or flat), the
injection construction is not doing what it claims and must be diagnosed before any leakage number is
reported as a "how much a missed confound costs you" statement.
*Falsifier B2.* If including the planted confound in the design does not reduce the channel back
toward the pre-injection baseline (within the run's noise band), the design-fitting side of the
machinery, not the injection side, is broken.

**C. Nonlinear-encoded arm.** Spec's own predeclared expectation, restated and adopted here: **the
purely additive/linear design does not fully catch a nonlinearly-encoded confound** — i.e. a
detectable residual channel/leakage remains after the additive design's adjustment that is **reduced**
(not necessarily eliminated) by the saturated-cell design. This is stated as the expected direction,
not assumed as the conclusion — if the additive design catches the nonlinear confound just as well as
the saturated design, that is the more interesting (favourable-to-the-paper) result and must be
reported plainly as such, not chosen after the fact.
*Falsifier C1.* If the saturated-cell design does **not** reduce the leaked channel relative to the
additive design (i.e. the saturated design, which spans strictly more functions, performs no better),
that contradicts the design-nesting logic already established in `nonlinear_adjustment.py`'s own
docstring and must be treated as a bug in this new code, not a finding, until re-checked.

## 3. Counterfactuals (rule §2.6)

* **Positive control:** an OVB run where the dropped/omitted covariate is a strong, real (by
  construction) confound of both blocks must show a clear, large channel/floor inflation — if it does
  not, the whole OVB machinery is not sensitive enough to trust on a real, smaller effect.
* **Must-fail control:** dropping/omitting a covariate that is pure independent noise must show no
  systematic inflation (Falsifier A1 above is this control, stated as a falsifier because failing it
  is informative either way).
* **Matched null:** every reported drift is computed against a same-seed, same-`n`, same-draw-count
  re-run of the *full* design (not a naive zero baseline), so that instrument-level draw noise is not
  mistaken for OVB signal.

## 4. Reporting rules

* Every arm is reported whether the result is favourable or not, bad news first.
* Every synthetic number is labelled `[SYNTHETIC]`.
* No inline reimplementation: `cross_fitted_residuals`, `confound_design`, `cell_design`,
  `saturated_cell_residuals`, `spike_targets`, `spike_recovery_curve`, `permutation_null` are all
  imported from their existing modules; `misspecification.py` contains only the three new functions
  named in spec §5 WS-A4 plus thin orchestration.
