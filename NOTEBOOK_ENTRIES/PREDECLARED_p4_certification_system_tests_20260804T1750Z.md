## 2026-08-04 17:50 UTC — PREDECLARED: four system-level tests of P4's certification rule, run before any interface exists

**Logged:** 2026-08-04 17:50 UTC. Nothing in this file has been computed. It is committed before
the harness is written so that the pass/fail readings below cannot be chosen after seeing a number.

P4's claim is *you cannot prompt what you cannot certify*, and its governing artifact is the
five-point rule of `v2/research/rebase/MULTIMODAL_EXPANSION.md` §1. That rule has never been
executed end to end on anything. It does not need an interface to be executed: the interface is a
query layer over a set of certified axes, and **the set of certified axes is computable today**.
These four tests compute it.

Standing constraints for every run below: CPU only, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` (other agents hold the GPU); workspace built with
`git -c core.autocrlf=false archive HEAD` and verified per-file by git blob SHA-1; **no statistic
computed inline** — every number comes from `v2/calibra/` (`confound_certificate.certify_axes`,
`calibration.spike_recovery_curve`, `spectral.heldout_single_direction_correlation`,
`residualise.confound_design` / `cross_fitted_residuals`, `confound_certificate.within_stratum_permutations`).
`v2/calibra/claim_guards.py` and `v2/research/rebase/nature/claim_evidence.json` are **not edited**:
discharging a blocker is a deliberate act, not a side effect of an analysis run.

---

### Fixed setup, shared by tests A–C

* **Artifacts.** `runs/d2_final/artifacts/d2_{h,i}_seed42.npz` and, where stated,
  `e0_run/d2_v3/d2_v3_s{43,44}/artifacts/d2_{h,i}_seed{43,44}.npz`, under
  `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/`.
* **Partition.** `split == "test"`, n = 2,766. The representation never saw these patients.
* **States.** `wsi_biology`, `rna_biology`, `full_biology`. Where one state must be named,
  it is **`wsi_biology`**: a promptable interface over a slide can only answer from the image-only
  channel, so `full_biology` and `rna_biology` are not exposable objects at query time.
* **Adjustment.** `confound_design(cancer + pooled TSS)` via `residualise`, `min_site_count = 10`,
  cross-fitted — the exact adjustment CALIBRA applies before it measures any channel.
* **Targets.** `data/frozen_rna_targets.npz`, 180 columns:
  `hallmark_in_training` 50, `heldout_pathway` 24, `immune_tme` 8, `tumour_state` 8,
  `random_control` 90. The **90 non-control** targets are the union of the first four groups.

---

## Test A — the abstention rule: what fraction of axes may a query layer answer from?

**Question.** Under condition 3 of the five-point rule alone, how many of the 256 axes of a D2 state
would a query layer be *allowed* to answer from?

**Definition, fixed now.** An axis is **answerable** iff all three hold on the state that would be
exposed:

* **A1 (per-axis)** its out-of-fold balanced accuracy for pooled TSS does not exceed the 95th
  percentile of the within-cancer label-permutation null (`n_permutations = 1000`);
* **A2 (per-axis CI)** it is not named in `breaching_axes`, i.e. it does not additionally carry a
  bootstrap CI excluding chance (`certify_axes` already ANDs A1 with A2 and treats a breaching axis
  with an unmeasured CI as breaching; that convention is kept, not softened);
* **A3 (joint)** the state's joint LDA test passes (`joint_certified`). A state that passes per-axis
  and fails jointly has not been shown to be site-free — it has been shown that nobody looked with
  enough resolution (`confound_certificate` module docstring; `t13_confound_certificate` entry).
  A3 is a property of the state, so it gates every axis in that state at once.

**Reported.** `answerable_fraction = |{axes satisfying A1∧A2∧A3}| / 256`, per artifact × state ×
seed, separately for the **raw** state (what a user would be shown) and the **adjusted** state
(what CALIBRA measures on). Both the strict reading (A1∧A2∧A3) and the per-axis-only reading
(A1∧A2) are reported, because they differ and the difference is the finding.

**Predeclared readings.**

* If the strict answerable fraction on the **raw** state is **0 in ≥5 of 6 artifact×state cells** →
  **P4 is not a system, it is a refusal**, under condition 3 alone, and must be written that way.
* If the strict answerable fraction on the **adjusted** state is **> 0.90** → that is *not* read as
  a pass until Test A′ below returns. A certificate that passes everything certifies nothing.
* Any intermediate fraction (0.1–0.9) is the interesting outcome and is reported as measured.

### Test A′ — does the certificate discriminate, or does it pass everything?

Three axes of known character are **appended** to each state before certification and scored by the
identical machinery. All three are built from labels already on the artifact, with a fixed seed:

1. `PLANT_site` — the pooled-TSS class index mapped through one fixed random real value per site,
   plus Gaussian noise at a signal-to-noise ratio of 1.0 (so it is a *noisy* site code, not a
   degenerate one). **Must breach.**
2. `PLANT_cancer` — the identical construction on cancer type instead of site. **Must not breach.**
   The null permutes site *within cancer*, so a pure lineage axis carries no site information beyond
   what cancer already explains; if this breaches, the null is mis-specified and the certificate is
   refusing ordinary biology.
3. `PLANT_noise` — `N(0,1)`, independent of everything. **Must not breach.**

**Pass = exactly the pattern (1 breaches, 2 and 3 do not) on the raw state.** Any other pattern
means the certificate does not discriminate and every "certified" verdict on this project is
uninformative in a direction that must then be named.

**The near-circularity check, predeclared with both outcomes.** The same three planted axes are
also certified with `--residualise`, i.e. after the cancer+pooled-TSS adjustment.

* If `PLANT_site` **stops breaching** after adjustment → the adjusted certificate cannot tell a pure
  site code from a biology axis, because the adjustment removes site by construction. "The adjusted
  state certifies" is then a statement about the adjustment, **not** evidence that any axis is
  site-free, and the adjusted arm may not be quoted as certification for exposure.
* If `PLANT_site` **still breaches** after adjustment → the adjusted certificate retains
  discriminating power and the adjusted pass recorded on 2026-08-03 is a substantive pass.

Both outcomes are informative and both will be reported. Neither is a defect of the instrument.

---

## Test B — the five conditions, end to end, on one real axis

The core gate of P4, executed start to finish for the first time.

**Axis selection, fixed before any result is seen.** State `wsi_biology`, artifact `d2_h_seed42`,
test partition, **adjusted** block. For every one of the 256 axes and every one of the 90
non-control targets, compute `spectral.heldout_single_direction_correlation(axis[:, None], target)`
(out-of-fold ridge, `n_splits = 5`, `alpha = 1.0`, `seed = 42`). The end-to-end axis is the axis
attaining the largest **absolute** value over that 256 × 90 grid, and its target is the
corresponding target. Rationale: this is the single most favourable case available. **If the best
case fails a condition, no weaker axis passes it**, so one axis suffices to decide the gate in the
negative; only a pass would require extending the run.

Each of the five conditions is scored **PASS / FAIL / UNEVALUABLE**, and an UNEVALUABLE condition
counts as *not passed* for the gate — you cannot certify on evidence that does not exist. The
distinction between a measured failure and absent data is preserved in the report.

1. **Operator estimated on a discovery fold only.** PASS iff (i) no `split == "test"` patient
   appears in `train`/`val`, verified on the artifact's own `split` vector, and (ii) the readout
   direction is fit out of fold — which `heldout_single_direction_correlation` enforces by
   construction. FAIL on any leak.
2. **Clears the CALIBRA detection floor.** `calibration.spike_recovery_curve` with
   `x` = the single axis (2,766 × 1), `y` = the single target column (2,766 × 1), `design` = the
   cancer + pooled-TSS design, default levels `(0, .01, .02, .05, .10, .20, .40)`, `n_draws = 25`,
   `seed = 42`. PASS iff `summary()["observed_above_floor"]` is `True` — the scale-safe flag, which
   compares the targeted single-direction observation against the unpaired `detection_floor` and is
   the only comparison the module permits. **FAIL if the floor is NaN**: an unresolvable floor is
   not a pass, it is the absence of one.
3. **Confound certificate.** From Test A, for this axis, on the raw state (the exposable object) and
   on the adjusted state, each with its state's joint row. PASS iff A1∧A2∧A3 hold on the exposed
   state.
4. **Replicates in untouched patients and ≥1 external cohort.** Two sub-conditions, scored
   separately, and **both** are required by the rule as written.
   * **4a untouched patients.** Beyond the test partition, the stricter reading the rule invites:
     leave-sites-out. Fold the test patients by pooled TSS into 5 groups of whole sites; fit the
     readout direction on 4 groups, score on the held-out sites, using
     `spectral.heldout_single_direction_correlation` on the site-grouped index split via the
     repository's leave-sites-out helper. PASS iff the out-of-site value has the **same sign** as
     the within-cohort value **and** retains **≥ 50%** of its magnitude. The 50% threshold is
     declared here and is not moved afterwards.
   * **4b external cohort.** `claim_guards.no_external_cohort` is undischarged for every morphology
     result on the project. The only external material on the box is `external/cptac_gdc_rna`,
     which is RNA-only STAR counts; HEST-1k spatial data exists but no D2 axis has been read on it.
     Predeclared **UNEVALUABLE**, conditional on verifying that no paired external H&E + RNA cohort
     exists on the box. If one is found, this is re-scored.
5. **Failures recorded and exposed alongside successes.** Procedural, scored by inspection, and
   labelled as procedural in the report. Met iff (i) `certify_axes` names `breaching_axes` rather
   than dropping them, (ii) `ClaimVerdict.as_rows` emits a visible status row for an inadmissible
   claim rather than silence, and (iii) every failure of this run is logged in the result entry.

**Gate.** The axis is certifiable iff 1–5 all PASS. Predeclared: **the report leads with the
conditions that fail**, and the count of failing conditions is the headline number, not the count
of passing ones.

---

## Test C — the competitor gap: what would a CellWhisperer-style interface answer that we refuse?

CellWhisperer (*Nature Biotechnology*, Nov 2025) ships a chat interface over CELLxGENE with no
uncertainty or abstention language. P4's contribution is precisely the refusals such a system does
not make, and that contribution currently has **no number**. This test produces one.

**Query set.** The 90 non-control targets of `frozen_rna_targets.npz`. Query *t* = "report target
*t* for this patient from the slide." State `wsi_biology`, `d2_h_seed42`, test partition, adjusted
block primary and raw block secondary.

**Policy N (competitor-style, no certificate).** Answer every query, from the axis with the largest
absolute out-of-fold correlation for that target. `n_answered = 90` by construction — that is what
"no abstention" means.

**Policy C (ours).** Answer query *t* only if there exists an axis satisfying **all** of:
* the confound certificate on the exposed state (A1∧A2∧A3 of Test A);
* `spike_recovery_curve` for that (axis, target) pair returns `observed_above_floor == True`;
* the axis's out-of-fold correlation for *t* exceeds the 95th percentile of a **within-cancer row
  permutation null** of the target column, built with
  `confound_certificate.within_stratum_permutations` (an imported permuter) and scored with the same
  imported correlation statistic, `n_permutations = 200`.

**Reported numbers.**
* `n_N = 90`, `n_C`, and the gap `n_N − n_C` — the count of confidently-worded answers our rule
  refuses. **That gap is P4's contribution expressed as a number.**
* Among Policy N's 90 answers, the count refused for each reason, reported separately and not
  collapsed: site-breaching support / below the detection floor / inside the permutation null.
  An answer may fail more than one; the report gives both the per-reason counts and the union.
* The distribution of the supporting correlation over the 90 answers, so the reader can see what
  magnitude of evidence the competitor policy would be asserting on.

**Predeclared readings.**
* If `n_C = 0` → the honest statement is that **our system answers nothing**, and P4's contribution
  is entirely the refusal. That is publishable as a negative-result system paper and must not be
  dressed as a working interface.
* If `n_C` is between 1 and 89 → the gap is the contribution and is quotable directly.
* If `n_C = 90` → our rule refuses nothing the competitor answers, the certificate is inert on this
  query set, and P4 has no measured contribution over CellWhisperer. This outcome is reported as
  such.

---

## Test D — does certification survive the spatial modality?

**Not run here.** A concurrent agent predeclared the spatial replication in
`NOTEBOOK_ENTRIES/PREDECLARED_spatial_claim_replication_20260804T1800Z.md` (committed `80076b4`),
whose Claim 1 and Claim 4 are exactly `certify_axes` on the HEST-1k spot artifact with **slide** in
place of tissue source site (13 test slides, chance 1/13 = 0.0769). Duplicating it would be waste
and would produce a second number under the same name.

What this task contributes instead, and it is the part their predeclaration does not cover:

* **D1.** Read their result when it lands and record whether the certificate *transfers* — i.e.
  whether the raw-FAIL / adjusted-PASS pattern found on TCGA at 85 site classes reproduces on HEST
  at 13 slide classes. Their Claim 1a fixes the bar (adjusted joint-LDA ≤ 1.5 × chance = 0.115).
* **D2.** State, from the artifact schema rather than from hope, which of the five conditions can
  even be *asked* of the spatial modality. Predeclared expectation, to be checked and corrected:
  conditions 1, 3 and 5 are askable; condition 2 is askable but expensive; condition 4b is **not**
  askable, because reading a *D2 axis* on HEST requires transporting the TCGA-trained representation
  onto HEST spots, which nothing on the project has done — and TCGA vs HEST embeddings separate at
  AUC 0.99999, so the transport is not free.
* **D3.** Record whether HEST discharges `claim_guards.no_external_cohort` **as a matter of fact
  rather than of mechanism**. The guard's discharge condition is mechanically trivial
  (`len(external_cohorts) >= 1`), and the guard's own remedy text names HEST-1k. That is not the
  same as having replicated a claim there. No edit to `claim_evidence.json` is made by this task.

---

## What this predeclaration does not do

It does not test P3. P3's hypothesis was refuted by its own predeclared test and the surviving
claim is narrow; the P3 work in this task is a written re-plan against what is now true, in
`paper/P3_P4_PLAN.md`, not a new experiment. No number in that plan will be new.

---

## AMENDMENT, 2026-08-04 18:35 UTC — condition 2's pass flag is sign-dependent in the single-axis case

Recorded **before the full run**, on evidence from a smoke run at reduced settings
(`--n-draws 4 --n-permutations 20`), and committed before any full-setting number was produced.

The predeclaration said condition 2 is PASS iff `SpikeRecoveryResult.summary()["observed_above_floor"]`
is `True`. That flag is `observed_matched_direction > detection_floor`, and
`observed_matched_direction` is a **signed** correlation along the randomly drawn direction pair
`(u, v)` the spike was planted on. With a single-column `x` and a single-column `y` — which is
exactly the per-axis certification unit — `u` and `v` are scalars, so the statistic is
`sign(u·v) · |corr(axis, target)|` and its **sign is a coin flip with no content**. On the smoke run
it returned **−0.4765** against a `detection_floor` of 0.05, i.e. the flag said FAIL for an
association whose magnitude is nearly ten times the floor.

**Amended criterion, fixed now.** Condition 2 is PASS iff

> `abs(spectral.heldout_single_direction_correlation(axis, target)) > SpikeRecoveryResult.detection_floor`

which is the pairing `heldout_single_direction_correlation`'s own docstring prescribes: *"the CALIBRA
`detection_floor` is expressed in single-direction correlation units … grading a per-target control
against that floor therefore requires a per-target statistic on the same scale."* A NaN floor remains
a FAIL. The shipped `observed_above_floor` flag and the signed `observed_matched_direction` are still
reported next to the amended verdict so the discrepancy is visible rather than papered over.

The same substitution applies to `refused_floor` in Test C, for the same reason.

**This is a defect in the shipped flag, not only in my use of it**, and it is recorded as one: for a
multi-column `x` the same median-over-random-directions is centred near zero, so
`observed_above_floor` is close to always-`False` there too. No library code is changed by this task
— the flag is left exactly as it is and the amended criterion is applied at reporting time, which is
where the module docstring says magnitudes belong ("Magnitudes are taken at reporting time, never
before the pairing").
