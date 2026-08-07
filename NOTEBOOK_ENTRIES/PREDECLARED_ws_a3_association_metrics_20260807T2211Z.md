# PREDECLARED — WS-A3, alternative association metrics (RV, dCor, HSIC, kernel CCA)

**Written before `v2/calibra/association.py` exists and before any number has been measured.**
Agent W1-B, executing `paper/P1_REVISION_SPEC.md` §5 WS-A3 (complaint #5). At the time of writing:
`v2/calibra/spectral.py`, `calibration.py`, `residualise.py`, `induced_correlation_sweep.py` have
been read in full. A repo-wide grep for HSIC, dCor, RV coefficient and kernel CCA returns zero hits,
confirmed independently of the spec's own claim. No association metric other than Pearson-on-CCA-
projections exists in this repository as of this entry.

## 0. Compute status — stated before anything else, because it constrains what follows

**CORRECTION, same day, before any measurement was taken.** `ssh p1box` initially refused
authentication from this checkout (public-key rejected); the coordinator subsequently reprovisioned
the alias and it is now reachable, verified by `ssh p1box hostname` and a GPU occupancy check
(idle A100, 0 MiB used). Real artifacts are present under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/` on the training box: the D2-final
representation artifacts (`runs/d2_final/artifacts/d2_h_seed42.npz`, `d2_i_seed42.npz`) and the frozen
RNA targets (`data/frozen_rna_targets.npz`), traced to real files with real SHA-256 hashes per
PROJECT_GUIDE §2.7 (recorded in the result entry). **This predeclaration now covers a REAL-DATA
measurement, not a synthetic-only fallback.** The synthetic verification in §3 (positive control,
must-fail control) remains mandatory and is run FIRST, before any metric is trusted on real data — that
sequencing is unchanged. The real-data protocol:

* **Cohort.** `d2_h_seed42.npz` / `d2_i_seed42.npz`, `test` partition (n≈2,530, matching every other
  P1 number quoted at that size), states `wsi_biology`, `rna_biology`, `full_biology`.
* **Targets.** `frozen_rna_targets.npz`'s 90 non-control (real) target columns, the same ones
  `run_calibra.py` scores everywhere else in this project.
* **Design.** Reused, not reimplemented: `induced_correlation_sweep.build_design` for the
  `cancer_tss_pool10` (Regime II / real) design and its `permuted_cancer_tss_pool10` and
  `gaussian_k99` falsifier arms (Regime I), exactly the arms `induced_correlation_sweep.py` already
  ships. `induced_correlation_sweep.py` is W1-A's file and is **imported**, never edited or copied.
* **Subsampling.** At n≈2,530 (test partition) no subsampling is needed for the O(n²) metrics under
  the declared `n>2000` rule; the rule itself is still implemented and unit-tested against a synthetic
  n=3,000 draw so it is exercised at least once this session.

## 1. The question, restated precisely

`induced_correlation_sweep.closed_form_induced_correlation` and spec §4 establish that residualising
two *unrelated* blocks against a shared confound design induces a nonzero Pearson correlation between
their residuals whenever the design explains a real, non-vanishing fraction of both blocks (Regime
II), and this induced correlation is n-invariant. The paper's account (Frisch–Waugh–Lovell) is a
statement about the **residual subspaces**, not about the statistic used to read association off
them — so the prediction is that RV, distance correlation, HSIC and kernel CCA, none of which is
Pearson-on-a-single-direction, should **all** show the same qualitative pattern: elevated,
n-stable/rank-explained residual association under Regime II, and small, rank-decaying residual
association under Regime I (structureless design of matched rank — Gaussian, or the real design with
rows permuted).

## 2. Kernel and bandwidth rule — declared BEFORE running, not fitted after

* **HSIC** (Gretton et al. 2005): RBF kernel on both sides, `k(a,b) = exp(-gamma * ||a-b||^2)`.
  **Median heuristic**: `gamma = 1 / (2 * median(pairwise squared Euclidean distance))`, computed
  separately for the X-side and Y-side kernels, each from its own pairwise-distance median. Biased
  empirical estimator `HSIC = (1/n^2) * trace(K H L H)`, `H = I - (1/n) * 11^T` (Gretton et al. 2005,
  eq. 4), because the biased estimator is the one the paper's target readers will expect by name and
  its bias is a constant offset under the null this analysis is designed to detect, not a source of
  false structure.
* **Kernel CCA** (Bach & Jordan 2002): RBF kernel on both sides, same median-heuristic gamma as HSIC.
  Regularised generalised eigenvalue problem solved via a symmetric (`scipy.linalg.eigh`) formulation
  with ridge `kappa = 1e-3` added to each centred Gram matrix as `kappa * n * I` (Bach & Jordan's own
  scaling convention, §4). Implementation deliberately avoids `numpy.linalg.svd`/`svdvals` (the tokens
  `test_effective_rank_canonical.py`'s AST scan polices) by using `eigh` on the symmetric generalised
  problem — this is a legitimate alternative solution method for kernel CCA, not a workaround, and is
  stated here so the choice cannot look retrofitted after the fact.
* **Distance correlation** (Székely, Rizzo & Bakirov 2007): unbiased/biased **V-statistic** form
  (double-centred Euclidean distance matrices), the standard definition; no kernel choice needed.
* **RV coefficient** (Josse & Holmes 2016 review of Escoufier 1973): no kernel; computed from the
  cross-product Gram matrices in closed form (`trace(Sx Sy) / sqrt(trace(Sx^2) trace(Sy^2))`).

**Subsampling rule, declared before any n=6,427-scale run is attempted** (HSIC, dCor and kernel CCA
are O(n^2) memory/O(n^3) compute; RV coefficient is not): if `n > 2000`, subsample rows via
`np.random.default_rng(seed).choice(n, 2000, replace=False)` before computing dCor/HSIC/kernel-CCA,
seed fixed and logged, never silently truncated. This rule is exercised at synthetic n up to 2,000 in
this predeclaration's own measurement (no real n=6,427 draw exists to subsample from this session) and
is implemented in `association.py` so it is ready for the real run.

## 3. Counterfactuals every metric must ship (rule §2.6 / spec §7.5)

1. **Positive control.** A planted linear association of known strength (`r_true` in
   `{0.0, 0.3, 0.6, 0.9}`, built exactly as `calibration.spike_targets` builds its spike) must be
   recovered by every metric in a monotonically increasing reading, before any confound design is
   applied. Failure of monotonicity for any metric is reported and that metric is flagged unreliable
   before it is used for anything else.
2. **Must-fail control.** Two independent (unrelated) Gaussian blocks, no shared design, must read
   near that metric's own null scale under every metric (checked against a permutation null of the
   same blocks, not an arbitrary absolute threshold, since HSIC/RV/dCor have no natural "zero").
3. **Matched null.** A label/pairing permutation null computed *per metric*, exactly analogous to
   `calibration.permutation_null`, reused via `residualise.cross_fitted_residuals` (never
   reimplemented) with the metric substituted for `top_canonical_correlation`.

## 4. Predeclared expectation and falsifier — stated before measuring

**Expected (favourable) result:** all four metrics show the Regime II / Regime I separation — a
residual association under the real-like structured design that is materially larger than, and does
not decay with `n` at the rate of, the same statistic under a rank-matched structureless design. This
is the paper's account holding up under a metric change.

**F1 (falsifier, per-metric).** For a given metric, if the Regime II reading is **not** distinguishably
larger than the Regime I reading at matched `n` and matched design rank (e.g. Regime II falls inside
the Regime I bootstrap/permutation interval in a majority of synthetic draws), that metric does **not**
detect the projection artifact predicted by FWL, **and this must lead the results section for that
metric, not be buried under the three that behave as predicted.** Per the task brief: "if HSIC or dCor
shows no floor, the FWL account is incomplete and that is a major finding."
**F2.** If a metric's Regime II reading visibly decays toward zero as `n` grows over the synthetic grid
`n in {200, 500, 1000, 2000}` (rather than staying flat/slowly-decaying the way the closed-form Pearson
account predicts), that is evidence the metric's induced-association floor is a finite-sample artifact
for that statistic rather than the population-parameter effect spec §4.3 describes, and must be
reported as a point of disagreement with the Pearson account, not smoothed into agreement.

## 5. What would make me distrust a favourable outcome

* All four metrics could trivially "pass" F1 if the synthetic Regime II design is made strong enough
  that any reasonable statistic detects it. Guard: the synthetic block design uses the **same**
  unbalanced one-hot block parameters as spec §4.4's own verified simulation (R=20 columns, matched to
  its n-grid where feasible), not a design hand-tuned to be maximally detectable.
* A metric passing only because its own null scale is inflated (e.g. HSIC's biased estimator has a
  strictly positive expectation under the null) is not evidence of the FWL effect specifically —
  guarded by always reporting the Regime II reading **relative to its own Regime-I-matched or
  permutation null**, never as a raw magnitude against zero.

## 6. Reporting rules

* Every metric is reported whether or not it confirms the paper's account.
* Every number in this line of work is prefixed `[SYNTHETIC]` in the result notebook entry and never
  merged into any real-data table.
* Nothing is computed inline that duplicates `calibra`'s existing statistics: only `rv_coefficient`,
  `distance_correlation`, `hsic`, `kernel_cca` are new; `cross_fitted_residuals`, `confound_design` and
  `spike_targets`/`_standardise`/`_correlation` (via import, not copy) come from existing modules.
