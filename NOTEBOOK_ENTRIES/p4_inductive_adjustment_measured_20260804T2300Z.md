## 2026-08-04 23:00 UTC — The inductive operator was built and wired, and the answer is still ZERO of 90: out of sample the confound adjustment leaves 24× more site signal than in sample, and the counterfactual's premise was the part that was wrong

**Logged:** 2026-08-04 23:00 UTC. **Predeclared in**
`NOTEBOOK_ENTRIES/PREDECLARED_p4_inductive_adjustment_20260804T2245Z.md`, committed (`af1c1a9`)
before the wiring was written and before any number below existed.

**How obtained.** Workspace `~/ws_p4i/morpheus` on the A100 (`150.136.45.194`), built from
`git -c core.autocrlf=false archive HEAD` at commit `5ac6d91` and verified **685/685 tracked files**
against an LF-normalised md5 manifest generated from that commit: **0 code files differ, 0 missing,
0 extra** (the 20 files reported differing are all `.pdf`/`.png` figures, whose manifest digest is an
artefact of applying LF-normalisation to binary content; the `phase_d --workspace-manifest` rule
checks `v2/ tests/ src/ configs/` × `.py|.json|.yaml`, of which **0** differ). CPU only,
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`; the GPU was at
99% for other agents throughout and was not touched. numpy 2.2.6, scikit-learn 1.7.2, scipy 1.15.3.

Driver `v2/research/rebase/nature/p4_certification/p4_certify.py`, which **defines no statistic**.
The adjustment operator is `v2/calibra/inductive_adjustment.ConfoundAdjustmentOperator`, used
unchanged. `claim_guards.py` and `claim_evidence.json` are **not edited**.

Artifacts `runs/d2_final/artifacts/d2_h_seed42.npz`, state `wsi_biology`, `test` partition;
targets `data/frozen_rna_targets.npz`, the same **90** non-control targets.

---

## Bad news first

**1. The buildable gap named at 20:00 has been built, and it does not open the interface. The
measured answerable count is 0 of 90 — the same as before.** Condition 3 still fails, but it now
fails on the *exposable* state itself rather than being unaskable.

**2. The counterfactual was right about the arithmetic and wrong about its premise.** "28 of 90 would
clear the floor if the operator existed" assumed that the inductively adjusted state would certify.
It does not. On the same 1,382 patients, the same 5-fold ridge, the same one-hot cancer + pooled-TSS
design and the same certificate:

| adjusted state on the exposure fold (n = 1,382, 34 pooled site classes, chance 0.02941) | joint LDA balanced accuracy | its within-cancer null p95 | axes breaching / 256 | certified? |
|---|---:|---:|---:|:---:|
| **raw** (no adjustment) | 0.4855 | 0.3195 | 46 | ✗ |
| **transductive** — nuisance model fitted on the rows it scores | **0.0109** | 0.0913 | **0** | **✓** |
| **inductive** — nuisance model fitted on the discovery fold only | **0.2643** | 0.1495 | **29** | **✗** |

**The in-sample adjustment reads 0.0109, below the chance rate. The out-of-sample adjustment reads
0.2643, nine times chance and 1.8× its own permutation null. That is a factor of 24.** The site
confound is not removed by this adjustment; it is removed from the rows the nuisance model was fitted
on. Everything the project has said about the adjusted state has been said about the first case.

**3. A larger discovery fold does not rescue it.** At a 0.7 discovery fraction — 1,937 fitting rows,
site coverage up from 43.8% to 54.8% — the inductively adjusted state reads joint LDA **0.2502**
against a null p95 of 0.1839, still not certified. The failure is not a discovery-fold-size artefact.

**4. Conditions 1 and 3 of the five-point rule cannot be discharged by the same split, as a matter of
fact about this artifact.** `train+val` carries 11 cancer types and 258 TSS codes, `test` carries 21
and 352, and the two sets have **0 cancer types and 0 TSS codes in common** — D2 holds out whole
cancers and TCGA assigns a TSS code per submitting site *per disease*. An operator
fitted on the representation's own discovery fold sees every test cancer as an unseen level and every
test site as `OTHER`, so its adjustment degenerates to subtracting a constant. The nuisance operator
therefore needs its own discovery fold, carved out of the certification cohort. This was measured
before predeclaring and is why the protocol is what it is.

**5. Nothing new was built that already existed.** `ConfoundAdjustmentOperator` *is* the operator P4
needed — same nuisance family (5-fold ridge, `alpha = 1.0`), same design, frozen levels, frozen site
pooling, bit-identical to `cross_fitted_residuals` on its fitting cohort. No second operator was
written. What was missing was a split and a wiring, and one measurement nobody had made.

---

## 1. What was built

`--adjustment {transductive, transductive_exposure, inductive}` on `p4_certify.py`.

* `transductive` (default) is the published 20:00 path, unmoved. `prepare_state` in that mode is
  asserted equal to `cross_fitted_residuals` with `np.array_equal`
  (`test_transductive_mode_still_equals_cross_fitted_residuals`).
* `transductive_exposure` scores the exposure fold and adjusts it **in sample** — the matched-n
  control. Without it, a drop in the inductive arm would confound *inductive vs transductive* with
  *n halved*.
* `inductive` fits a `ConfoundAdjustmentOperator` on the discovery fold and applies it unchanged.
  Only in this mode is the adjusted state marked **exposable**, and only then does the adjusted arm
  decide condition 3. A transductively adjusted coordinate cannot be produced for a new patient, so
  a user is still shown a raw axis there — the semantics of the published run are preserved.

The split permutes patients **within cancer type** (seed 42), so every cancer appears on both sides;
`DesignSpec.transform` is left on `on_unseen_level="refuse"` so a split that broke that rule would
stop the run rather than silently adjust those rows as the design's baseline level.

**The published path is verified unmoved, not asserted unmoved.** `--adjustment transductive` was
re-run end to end on all 2,766 test rows and its output compared with the vendored 20:00 artifact
`out/C_competitor_d2h_seed42_wsi.json` field by field: **0 mismatches across 90 rows × 15 graded
fields** (`axis`, `correlation`, `abs_correlation`, `detection_floor`, `null_p95`, `null_median`,
`axis_balanced_accuracy`, `axis_null_p95`, `observed_matched_direction`, `clears_detection_floor`,
`axis_breaches_certificate`, `refused_site`, `refused_floor`, `refused_null`,
`answered_by_policy_C`), and identical `joint_lda_balanced_accuracy` 0.36326266849728855, identical
refusal counts (90 / 62 / 1), identical correlation quantiles, identical 23 NaN floors. The wiring
did not move a published number.

**Two tests carry the claim, and they are a pair.**
`test_an_exposure_row_is_unmoved_by_the_other_exposure_rows` replaces every exposure row but the
first with different data and asserts row 0's adjusted coordinates are **bit-identical**
(`array_equal`). `test_the_matched_transductive_control_does_NOT_have_that_property` runs the
identical manipulation on the control and asserts it moves. The first test cannot be passing for a
trivial reason.

---

## 2. The measurement — Test C, re-run for real

Query set, policy and statistic identical to the 20:00 entry. Policy C answers a query iff its
supporting axis passes the confound certificate **on the exposed state**, *and*
`|heldout_single_direction_correlation|` exceeds that (axis, target) pair's own CALIBRA detection
floor, *and* it exceeds the 95th percentile of a 200-draw within-cancer row-permutation null.

| | 20:00 run<br>transductive, n = 2,766 | matched control<br>transductive, n = 1,382 | **inductive**<br>**n = 1,382** | inductive<br>n = 829 (f = 0.7) |
|---|---:|---:|---:|---:|
| exposed state | raw | raw | **adjusted (inductive)** | adjusted (inductive) |
| queries | 90 | 90 | **90** | 90 |
| answered by Policy N (competitor-style) | 90 | 90 | **90** | 90 |
| **answered by Policy C — MEASURED** | 0 | 0 | **0** | 0 |
| gap | 90 | 90 | **90** | 90 |
| refused: site certificate | 90 | 90 | **90** | 90 |
| refused: below the detection floor | 62 | 57 | **61** | 64 |
| refused: inside the permutation null | 1 | 7 | **13** | 26 |
| **with condition 3 relaxed** (the 20:00 counterfactual) | **28** | 31 | **28** | 23 |

**`n_C = 0`. That is the number this work was commissioned to replace the projection with, and the
projection's headline arithmetic survives while its premise does not.**

### The category breakdown, reported the way the original entry did

Counts are of queries that clear **floor and null** — i.e. condition 3 relaxed, which is the only
level at which any query survives in any arm. Under the full policy every cell is **0**.

| target category | queries | 20:00 counterfactual (n = 2,766) | matched transductive (n = 1,382) | **inductive (n = 1,382)** | inductive (n = 829) |
|---|---:|---:|---:|---:|---:|
| `hallmark_in_training` | 50 | 19 | 20 | **18** | 15 |
| `heldout_pathway` | 24 | 1 | 4 | **2** | 2 |
| `immune_tme` | 8 | 5 | 5 | **5** | 5 |
| `tumour_state` | 8 | 3 | 2 | **3** | 1 |
| **total** | **90** | **28** | **31** | **28** | **23** |

**The 20:00 entry's qualitative finding holds under a real operator: what clears the floor is mostly
what the model was supervised on.** 18 of 28 are `hallmark_in_training`; **2 of the 24 genuinely
untrained `heldout_pathway` targets survive**, against 1 in the projection. The immune panel is the
one place the representation does better than its training set predicts — 5 of 8 `immune_tme`
targets clear in every arm.

The identity of the 28 is not the identity of the projection's 28: the argmax axis/target moves from
axis 46 / `HALLMARK_ALLOGRAFT_REJECTION` (0.4703 at n = 2,766) to **axis 133 /
`immune_cytolytic_activity`, 0.4683** on the inductive exposure fold, the observed-correlation median
falls from 0.2392 to 0.2058, and the number of targets with an unresolvable (NaN) detection floor
falls from 23 to **16**. The total coinciding at 28 is a coincidence, and is reported as one.

### Why 61 of 90 fall below the floor

Per-(axis, target) detection floors on the inductive exposure fold, from `spike_recovery_curve` on
the level grid (0, .01, .02, .05, .10, .20, .40), 25 draws:

| detection floor | 0.01 | 0.02 | 0.05 | 0.10 | 0.20 | 0.40 | **unresolvable (NaN)** |
|---|---:|---:|---:|---:|---:|---:|---:|
| targets | 4 | 4 | 4 | 9 | 25 | 28 | **16** |

against observed `|correlation|` running 0.0596 to 0.4683, median 0.2058 (p05 0.0756, p25 0.1135,
p75 0.2851, p95 0.4200). The floors are computed **from the same rows and the same design in both
exposure arms**, so they are near-identical between them and cannot explain the difference between
the inductive and the transductive arm. That the floor machinery is itself transductive is a stated
limitation of the audit, not of the exposed state.

90 queries are supported by **35 distinct axes**; **12 queries are supported by an axis that itself
breaches the certificate** on the exposed state (5 distinct axes). The comparable figure for the
matched control, whose exposed state is raw, is 8 queries on 5 axes; the 20:00 run's 10 on 3 axes is
also a raw-state figure, so the three are not like-for-like and no trend is claimed. The strongest
queries cluster as before: `immune_cytolytic_activity` 0.468 and `immune_ifng` 0.436 on axis 133,
`immune_t_cell_inflammation` 0.460 on axis 173, `HALLMARK_ALLOGRAFT_REJECTION` 0.459 and
`HALLMARK_INTERFERON_GAMMA_RESPONSE` 0.442 on axis 46.

---

## 3. Test B — the five conditions on the exposure fold

Axis selected before any condition was scored, as the argmax of
`|heldout_single_direction_correlation|` over 256 axes × 90 targets on the **inductively adjusted**
block: **axis 133** against **`immune_cytolytic_activity`** (`immune_tme`), out-of-fold
single-direction correlation **0.46831**.

| # | condition | inductive arm | matched transductive control |
|---|---|---|---|
| 1 | operator estimated on a discovery fold only | **PASS** — 1,384 discovery / 1,382 exposure, **0 patients in both**, asserted in code; representation split 3,118/543/2,766 with 0 patients in both | PASS |
| 2 | clears the CALIBRA detection floor | **PASS** — `detection_floor` **0.01**, observed **0.4683**, i.e. 47× the floor | PASS (floor 0.01, observed 0.4633) |
| 3 | passes the confound certificate | **FAIL on the exposed (inductively adjusted) state** — the axis is individually innocent (0.04714 against its own null p95 of 0.04730, a margin of 0.00016) but the **state's joint LDA is 0.2643 against a null p95 of 0.1495**, with 29 of 256 axes breaching | FAIL — exposed state is raw (joint 0.4855 vs 0.3195, 46 breaching). Its *adjusted* state certifies (0.0109, 0 breaching) but is not exposable |
| 4a | replicates in untouched patients (whole sites held out) | **PASS** — within-cohort 0.4773, median out-of-site 0.4950, **retained fraction 1.037**, 5/5 folds clear their own null and their site-cluster CI clears it too | PASS — 0.4795 / 0.4950 / **1.032**, 5/5 and 5/5 |
| 4b | ≥ 1 external cohort | **UNEVALUABLE** — unchanged; no paired external H&E+RNA cohort on the box | UNEVALUABLE |
| 5 | failures recorded and exposed | **PASS (procedural)** — unchanged | PASS |

**Gate: not met, 4 passed / 1 failed / 1 unevaluable, in both arms.** The condition that fails is the
same one, and it fails for a *new* reason in the inductive arm: not "the state that certifies cannot
be exposed", but "the state that can be exposed does not certify".

**Note the near-miss on the axis itself.** Axis 133's per-axis site accuracy on the inductively
adjusted state is 0.04714 against a null p95 of 0.04730. It clears by 1.6e-4. A per-axis-only
certificate would have certified it by a margin no one should be willing to defend, and the joint row
refuses it outright — the same lesson as the 20:00 entry, on a different axis and a different state.

---

## 4. Where the residual site signal lives, and the checks the predeclaration required

### §6.1 — is the adjustment doing anything at all?

| arm | median per-axis corr(raw, adjusted) | axes with corr > 0.99 | median residual variance ratio | design width |
|---|---:|---:|---:|---:|
| transductive, n = 2,766 | 0.7443 | **0** | 0.5918 | 108 |
| transductive_exposure, n = 1,382 | 0.7524 | **0** | 0.5931 | 57 |
| **inductive, n = 1,382** | **0.7536** | **0** | **0.5896** | 55 (operator) / 57 (exposure) |
| inductive, n = 829 | 0.7436 | 0 | 0.5742 | — |

**The inductive and the transductive adjustment remove the same amount of variance** — 41.0% vs
40.7% of the median axis — and neither leaves any axis untouched. So the inductive arm's failure is
**not** a null adjustment: it removes as much as the transductive one does, and removes *the wrong
part*. That is the finding, and it is the one the audit was written to be able to distinguish.

### Coverage is not the whole story — the covered rows fail too

The operator can only site-adjust an exposure row whose site was frequent **in the discovery fold**:
606 of 1,382 rows (43.8%) at `min_site_count = 10`. The obvious hypothesis is that the failure is
entirely those 776 uncovered rows. It is not. `site_coverage_probe.py` splits the exposure fold on
exactly that line and runs the *same* certificate on each side (1,000 permutations, 200 bootstrap
draws; site classes and chance rate are re-derived within each stratum, so the strata are comparable
to their own nulls and not to each other):

| stratum | n | site classes | chance | arm | joint LDA | its null p95 | certified? | breaching / 256 |
|---|---:|---:|---:|---|---:|---:|:---:|---:|
| all | 1,382 | 34 | 0.0294 | raw | 0.4855 | 0.3195 | ✗ | 46 |
| all | 1,382 | 34 | 0.0294 | transductive | 0.0109 | 0.0913 | ✓ | 0 |
| all | 1,382 | 34 | 0.0294 | **inductive** | **0.2643** | 0.1495 | **✗** | 29 |
| **site-adjustable** | **606** | 26 | 0.0385 | raw | 0.5540 | 0.4035 | ✗ | 39 |
| **site-adjustable** | **606** | 26 | 0.0385 | transductive | 0.0237 | 0.0788 | ✓ | 0 |
| **site-adjustable** | **606** | 26 | 0.0385 | **inductive** | **0.2836** | 0.2004 | **✗** | 30 |
| pooled to `OTHER` | 776 | 9 | 0.1111 | raw | 0.7026 | 0.6236 | ✗ | 15 |
| pooled to `OTHER` | 776 | 9 | 0.1111 | transductive | 0.0081 | 0.2208 | ✓ | 0 |
| pooled to `OTHER` | 776 | 9 | 0.1111 | **inductive** | **0.3972** | 0.3238 | **✗** | 8 |

**The rows the operator *could* site-adjust also fail**: 0.2836 against their own null p95 of 0.2004,
7.4× their chance rate — while the in-sample adjustment on the *same 606 rows* reads 0.0237 and
certifies with 0 breaching axes. So the residual is not merely the uncovered rows: it is also
estimation error on the covered ones, where a site's mean comes from a handful of discovery patients
instead of from the rows being scored. That is why a larger discovery fold (§ bad news 3) does not
rescue it either. The `all`-stratum rows reproduce the Test B certificate to every quoted digit,
which is the probe's internal consistency check.

### §6.2 — is the certificate on the inductively adjusted state inert, like the transductive one?

**No, and this is the check that makes the FAIL in §3 credible.** The 20:00 entry showed that the
*transductively* adjusted certificate cannot tell a planted site code from a biology axis at any
realistic strength — it stops refusing `PLANT_site` at SNR 0.5, 1.0 and 2.0 — so a PASS from it says
nothing. The same ladder, re-run on the **inductively adjusted** state (259 axes = 256 real + a site
code, a cancer code and pure noise; 1,000 within-cancer permutations; n = 1,382; chance 0.02941):

| SNR | arm | `PLANT_site` acc | its null p95 | perm p | refuses site? | `PLANT_cancer` refused? | `PLANT_noise` refused? | state joint LDA |
|---|---|---:|---:|---:|:---:|:---:|:---:|---:|
| 0.5 | raw | 0.0576 | 0.0493 | 0.0110 | exceeds null, CI wide | no | no | 0.5124 |
| 0.5 | **inductively adjusted** | 0.0391 | 0.0468 | 0.2418 | **no** | no | no | 0.2463 |
| **1.0 (predeclared)** | raw | 0.0872 | 0.0576 | 0.0010 | exceeds null, CI wide | no | no | 0.5706 |
| **1.0 (predeclared)** | **inductively adjusted** | 0.0599 | 0.0489 | 0.0050 | **YES** | no | no | 0.2542 |
| 2.0 | raw | 0.1248 | 0.0709 | 0.0010 | **YES** | no | no | 0.6618 |
| 2.0 | **inductively adjusted** | 0.0865 | 0.0518 | 0.0010 | **YES** | no | no | 0.2761 |
| ∞ | raw | **0.9706** | 0.1029 | 0.0010 | **YES** | no | no | 0.7183 |
| ∞ | **inductively adjusted** | **0.9706** | 0.0749 | 0.0010 | **YES** | **YES** (false refusal) | no | 0.3351 |

**Three things follow.**

* **The instrument is live on the exposed state.** At the predeclared SNR of 1.0 the inductively
  adjusted certificate **refuses a planted site code** and accepts a cancer code and pure noise. The
  transductively adjusted certificate at the same strength does not (0.0090 against a null p95 of
  0.0183, p = 0.86, per the 20:00 ladder). So the FAIL reported above comes from a certificate that
  can still discriminate, not from one that refuses everything.
* **A noiseless site code passes through the inductive adjustment essentially untouched: 0.9706 raw
  → 0.9706 adjusted.** The transductive adjustment took the same construction from 1.0000 to 0.2751.
  That is the §2 finding again, at maximum contrast: a nuisance model that must be fitted on other
  patients cannot remove a site effect for the 56% of rows whose site it never had a column for.
* **The false-refusal caveat reproduces.** At SNR = ∞ the inductively adjusted arm refuses
  `PLANT_cancer` (0.9118, p at the resolution floor). Cancer *is* fully covered by the operator's
  design, so the axis is annihilated and the certificate's per-axis standardisation rescales the
  numerical residue to unit variance. Same mechanism as the 20:00 entry recorded for the transductive
  arm; it is a property of adjusting-then-standardising, not of this operator.

**One honest complication in the per-axis flag at this n.** `certify_axes` marks an axis
`breaching` only if it exceeds its null p95 **and** its bootstrap CI excludes chance. At n = 1,382
the CI is wide enough that the raw arm's `PLANT_site` exceeds its null at p = 0.0010 and is still not
flagged breaching at SNR 0.5 and 1.0. `exceeds_null_p95` and the permutation p are the informative
columns at this cohort size, and both are reported above.

### §6.3 — the split

0 patients appear in both folds; asserted in `prepare_state` and re-asserted by
`test_the_inductive_operator_never_saw_the_rows_it_scores`.

### §6.4 — did the count beat the counterfactual?

No: 28 against 28 at half the n, and 28 against the matched control's 31 on the identical rows. The
inductive arm is uniformly the more conservative of the two (61 vs 57 below floor, 13 vs 7 inside the
null), which is the direction a noisier out-of-sample residual predicts. No red flag.

### §6.5 — did `heldout_pathway` move up?

From 1 to 2 of 24, inside the "≤ 3" band the predeclaration fixed as unremarkable, and accompanied by
no rise in the trained categories (19 → 18). No accident to report.

---

## 5. How the predictions did

Graded verbatim against §5 of the predeclaration.

| predicted | measured | |
|---|---|:---:|
| the inductive state FAILS the joint certificate (p = 0.65) | **FAILS** | ✓ |
| `n_C` full policy, inductive = **0** | **0** | ✓ |
| joint LDA, inductive: 0.06–0.20 | **0.2643** | ✗ **worse than predicted** |
| joint LDA, transductive_exposure ≤ 0.03, certified | **0.0109**, certified | ✓ |
| axes breaching, inductive: 5–40 | **29** | ✓ |
| `n_C` relaxed, inductive: 8–20 | **28** | ✗ (above the band) |
| `n_C` relaxed, transductive_exposure: 10–22 | **31** | ✗ (above the band) |
| ≥ 60% of survivors `hallmark_in_training` | 18/28 = **64%** | ✓ |
| 0 or 1 `heldout_pathway` | **2** | ✗ (just outside) |
| targets with a NaN floor: 25–40 | **16** | ✗ **better than predicted** |

**Four of ten point predictions missed, and they miss in a consistent direction:** I over-predicted
how much halving n would cost the channel (the floors got *better*, not worse, and more queries
cleared them than I allowed for) and under-predicted how much site signal an out-of-sample adjustment
leaves behind. The one prediction the conclusion rests on — that the exposable state would not
certify — was right, and it was right by more than I expected.

---

## 6. What this changes

**For P4.** The interface is still a refusal, and the reason has moved from "not yet built" to
"measured and failing". The 20:00 entry could say condition 3 was a *build item*. It is now a
measurement item with a measured answer: **the exposed adjusted state does not certify, at either
discovery fraction tested.** A certified promptable interface over this representation requires an
adjustment that survives out of sample, and cross-fitted ridge on a one-hot cancer + site design is
not it.

**For P1, and this is the part that reaches beyond P4.** Every CALIBRA number computed on an
"adjusted" state was computed transductively. On these artifacts the transductive adjustment reads
0.0109 for joint site recoverability and the out-of-sample one reads 0.2643 on the same patients.
The supportable sentence is now narrower than the one in
`t13_adjusted_certificate_and_p6_20260803T0300Z.md` and `paper/P1_CALIBRA_DRAFT.md` §4.2, and
narrower even than the correction already recorded in
`prose_corrections_first_moment_20260804T2145Z.md`:

> the confound is removed from the first moment **of the rows the nuisance model was fitted on**, and
> a mean-based certificate evaluated on those same rows therefore certifies.

This sits alongside, and is independent of, the two limits already on record — that the adjustment
removes only the first moment (the spatial 15-NN probe reads 0.7291 on an adjusted state) and that
the adjusted certificate is near-inert against a planted linear site code. **Three separate routes
now say the adjusted state's blanket PASS is a statement about the adjustment, the classifier family
and the fitting rows, and not evidence that any axis is confound-free.** No edit is made here to
`claim_evidence.json` or to the paper; the correction is recorded and left for the owners of those
files.

---

## 7. Suite

Run on this workspace at the commit the runs were launched from (`5ac6d91`):
`pytest morpheus/v2/tests morpheus/tests --ignore=morpheus/v2/tests/test_p2_figures.py -q` →
**565 passed, 0 failed in 64.05 s**. `test_p2_figures.py` run separately reads **1 passed, 27 errors
in 2.59 s**, every error `ModuleNotFoundError: No module named 'matplotlib'` — the known condition of
`~/venv`. **Nothing was installed into that environment.**

The 15 new tests are the whole delta: `test_p4_inductive_wiring.py` alone reads **15 passed**, and
`test_inductive_adjustment.py` (the operator's own 20 tests, untouched) still reads **20 passed** —
35 together.

## 8. Files / provenance

Harness `v2/research/rebase/nature/p4_certification/p4_certify.py` (commit `5ac6d91`).
Tests `v2/tests/test_p4_inductive_wiring.py`, 15 tests.
Outputs `~/ws_p4i/out/{A_plant_inductive_f50, B_inductive_f50, B_transductive_exposure_f50,
C_inductive_f50, C_inductive_f70, C_transductive_exposure_f50, C_transductive_REGRESSION,
site_coverage_probe}.json`, vendored into `v2/research/rebase/nature/p4_certification/out/`.
Probe driver `site_coverage_probe.py`, vendored beside them.
Operator `v2/calibra/inductive_adjustment.py`, used unchanged. Reference digest of the f = 0.5
discovery-fold operator, as fitted for Tests B and C:
`2060a635fa83756a1c3b7aa8506b7b19fcc4431f5d1a303da39b3cb2bf9d62ce`. The planted-axis ladder fits its
own operator on the 259-column augmented matrix and therefore carries a different digest
(`6b3a9fdc…`) — correctly, since the digest covers the matrix shape as well as the design and the
patient identifiers.
