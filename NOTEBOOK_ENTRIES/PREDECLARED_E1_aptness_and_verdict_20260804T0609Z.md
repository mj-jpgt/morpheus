# 2026-08-04 06:09 UTC — E1 is not apt for the surviving claim, and is not runnable on any artifact that exists. Predeclaration first, then the verdict.

**Logged:** 2026-08-04 06:09 UTC (real UTC clock on this machine *and* on `150.136.45.194`; note that
several entries in this directory carry same-day timestamps later than the current clock — that
inconsistency predates this entry and is not corrected here).
**How obtained:** read of `v2/calibra/e1_rank_information.py`, `v2/calibra/aggregate_e1.py`,
`v2/research/rebase/run_e1_training.py`, `v2/tests/test_e1_rank_information.py`,
`v2/tests/test_aggregate_e1.py`; manifest inspection of the six D2 and six D1-B frozen artifacts on
`ubuntu@150.136.45.194`. CPU only; the A100 was at 99% on another chain throughout and was not touched.

---

## 0. Predeclaration — written before any E1 output was produced or inspected

This is recorded first because the project's standard (`PREDECLARED_D1_necessity_test_20260803T2300Z.md`,
outcome **O2**, reported at abstract prominence when it went against us) is that the reading is fixed
before the number is seen. **No E1 output was inspected in the course of this work, because none
exists**: the only `e1_arrays.npz` files anywhere on the box are inside `.pytest-e1*/` temporary
directories from the unit test, on synthetic data, with `gates_pass == False` asserted by the test
itself (`test_e1_rank_information.py:103`). Nothing was run. Nothing was looked at.

**The claim E1 would have had to discharge**, verbatim from `paper/P2_RANK_DRAFT.md` §1.3:

> Effective rank is unusable as a selection signal because its between-arm differences are smaller
> than its own within-arm reproducibility floor — inside the regime its proponents reserve for it.

**What would have supported it.** A measurement in which, for a fixed arm, the *same configuration
re-instantiated* (new training seed, or a retrain at the same seed) produces a spread in the canonical
statistic `spectral.CANONICAL` — Roy & Vetterli order 1, column-centred, rows at own norms — that
equals or exceeds the between-arm difference in that same statistic, on the same block, with both
quantities computed by the same imported function.

**What would have refuted it.** The between-arm difference exceeding the within-arm spread with a
margin, on the canonical statistic, in the majority of pairs; that is, effective rank resolving its
own arms.

**What would have been uninformative.** (i) Any comparison whose uncertainty band comes from
resampling *patients* rather than re-instantiating *training* — §4.1(iv) has already established that
patient resampling gives an SD of ≈0.1 on a rank of 25 while retraining spans 2.69×, so a
patient-bootstrap interval is a measurement of the wrong variance component and cannot bear on this
claim in either direction. (ii) Any result in which the rank statistic and the comparator are not both
stated with their block (raw vs residualised) and their variant name.

**What would have made me distrust a favourable result.** Predeclared, per instruction, in advance:
(a) a favourable verdict that depended on the *choice* of block or variant — i.e. that reversed under
`RANK_VARIANTS["R2"]`/`["R3"]` or under residualisation, both of which are known to flip arm orderings
on two D2 seeds; (b) a within-arm floor estimated from `n = 1` retraining pair, which is exactly the
limitation already conceded for the 2.69× envelope in §4.1(i); (c) any number produced by a formula
written inline rather than imported from `v2/calibra/spectral.py`.

---

## 1. Verdict: **E1 is not apt for this claim. It was not run.**

E1 asks a different question, and its decision rule cannot represent the current claim's supporting
outcome.

**1.1 E1's question is the falsified one.** The module docstring states it exactly: *"test whether
decorrelation-created rank carries molecular information"* (`e1_rank_information.py:1`). Its
aggregated verdict variable is
`supports_rank_without_detectable_information` (`aggregate_e1.py:48`). That is the
*rank-does-not-track-information* claim — the claim the D1 necessity test falsified and which §1 of
the draft records as removed. The surviving claim is not about whether rank tracks information; it is
about whether rank can **resolve its own arms above its own re-measurement noise**. Those are
different measurements of different quantities, and discharging one with the other is precisely the
substitution this project caught at §3.1 and again at §4.5(a).

**1.2 E1 contains no within-arm reproducibility term at all.** The claim is a ratio: between-arm
difference over within-arm floor. E1 measures the numerator (`delta_effective_rank`, `aggregate_e1.py:12`)
and never the denominator. Its only uncertainty bands are `paired_bootstrap_deltas`
(`e1_rank_information.py:218`), which resample **patients** and **cancer clusters** — the evaluation
cohort — with training held fixed. Under the predeclaration above that is the *uninformative* case,
and it is uninformative for a reason the paper has already measured: §4.1(iv) puts patient-resampling
SD at ≈0.1 against a retraining envelope of 2.69×. E1's three seeds are three *before/after pairs*,
scored for sign agreement; at no point does it compare an arm to a re-instantiation of itself.

**1.3 E1's aggregator hard-codes the sign the claim disputes.** `aggregate_e1.py:38` computes
`rank_positive = (frame.delta_effective_rank > 0).all()` and line 48 makes it a **necessary
conjunct** of the reported verdict. Under the surviving claim, a between-arm rank delta that fails to
hold its sign across seeds is the *supporting* observation. E1 would score that as
`all_seeds_rank_increase = False` and report `supports_rank_without_detectable_information = False` —
i.e. the design's "claim not supported" output is produced by both the supporting and the refuting
outcome of the current claim. A decision rule that maps two opposite worlds onto the same verdict
cannot discharge a claim about either. This is not a defect in E1; it is a design built for a
different hypothesis, and the draft already says so at §3.1 while nevertheless calling it *"the
experiment this paper should have rested on"*.

**Per the instruction governing this work, I stopped here rather than running it and presenting the
output as evidence for the surviving claim.**

---

## 2. Independent of aptness, E1 is not runnable on any artifact that exists

Three hard stops, each of which fires before any number is computed. Verified by reading the
manifests of the six D2 (`~/e0_run/d2_v3/d2_v3_s4{2,3,4}/artifacts/`) and six D1-B
(`~/e0_run/d1_v2/artifacts/`) frozen artifacts on the box.

**2.1 The preregistered intervention does not exist in any pair.** `_validate_intervention`
(`e1_rank_information.py:73-90`) requires the declared field to be `0.0` in the before arm and `> 0.0`
in the after arm, and raises otherwise. Measured:

| pair | `run_configuration.decorrelation_weight` before → after | E1 outcome |
|---|---|---|
| D2 arm H vs arm I, seed 42 | **0.04 → 0.04** | `ValueError: E1 intervention ... must be 0.0 before and >0.0 after` |
| D1-B `d1_p` vs `d1_f`, seed 42 | **0.04 → 0.04** | same |

D2's arms differ in `programme_targets` (`hallmark` vs a perturbation basis); D1-B's differ in
`objective_profile` (`programme_only` vs `programme_free`). Neither is the decorrelation intervention.

**2.2 D2's pair additionally fails the matched-manifest guard.** `validate_matched_artifacts`
(`e1_rank_information.py:113-150`) requires exact manifest equality after removing only the
decorrelation field. Arm I's manifest carries **60+ keys absent from arm H** (the whole
`programme_targets.manifest.*` block) and **22 differing shared keys**. The pair is refused with
`matched artifacts have non-intervention manifest differences`. D1-B's pair is closer — it differs in
four shared keys — but still fails, on `objective_profile`.

**2.3 There is no frozen targets file of the shape `_load_targets` requires.** `_load_targets`
(`e1_rank_information.py:260`) needs an `.npz` with `patient_ids`, `scores`, `target_names`. No such
file exists under `~/e0_run/`; the P2 scripts use `~/e0_run/data/frozen_rna_targets.npz`, which is a
different contract and is consumed through `p2_competing_metrics.load_targets`.

**Loosening any of these to admit the available artifacts would be rewriting a preregistration to fit
the data**, and would in any case convert E1 into the D1 necessity test — same arms, same seeds — which
has already been run and is reported at §4.7 with a better-suited estimator (a stratified paired
bootstrap). No such loosening was made; `e1_rank_information.py` and `aggregate_e1.py` are unmodified.

---

## 3. Discrepancies found against documents, reported rather than worked around

**3.1 "It is CPU work" is false.** `NOTEBOOK_ENTRIES/p2_rank_draft_20260803T2134Z.md:278` says of E1:
*"It is built, preregistered, three-seed, equivalence-margin, gate-enforced, and it is CPU work."*
E1's analysis stage is CPU work, but its **inputs do not exist**. `v2/research/rebase/run_e1_training.py`
is the only driver that can produce them — its docstring says so: *"only this driver can produce the
matched three-seed artifacts required for an E1 claim"* — and it launches
`morpheus.v2.runner` twice per seed with `--decorrelation-weight` 0 and >0
(`run_e1_training.py:60, 121-122`). Running E1 therefore costs **six GPU trainings**, not a CPU
afternoon. That materially changes the cost/benefit sentence in the entry and in `paper/P2_FIGURES.md`
S6, and it should be corrected wherever the "CPU work" framing is repeated.

**3.2 "Never run" is confirmed, and is the one claim here that checks out.** No `E1_*` rows in
`v2/research/rebase/nature/GATE_LOG.md`; no E1 outputs under `runs/`; on the box, every `e1_*.npz` is
inside a `.pytest-e1*/` tmpdir. `paper/P2_RANK_DRAFT.md` §3.1 and `paper/P2_FIGURES.md` S6 are accurate
on this point.

---

## 4. What an apt experiment would look like — stated, not run

For the record, so that the gap is named rather than filled by a substitute: the apt design compares,
in `spectral.CANONICAL`, on a stated block, (a) the between-arm difference for a pair against (b) the
spread of the same statistic across re-instantiations of *one* of those arms. §4.1(ii)'s variance
decomposition over the 12 artifacts (34.5% arm / 65.5% training-seed nuisance, F(3,8) = 1.41) is
already that measurement in its best-powered form, and §4.1(i)'s 2.69× envelope is its weakest
(`n = 1` retraining pair). The measurement that would actually strengthen the claim is **more
retraining repeats of a fixed configuration**, which is GPU work and is already armed as the
retraining envelope (`PREDECLARED_retraining_envelope_20260804T0330Z.md`). E1 is not that experiment
and cannot be made into it by rerunning it on more seeds.

---

## Files / commits

Read-only with respect to `v2/calibra/`: **no source file was modified**. `NOTEBOOK.md` and
`paper/P2_RANK_DRAFT.md` untouched. Suite unchanged at this point: **336 collected, 335 passed**, the
single failure being `test_paper_paths_resolve.py::test_no_box_output_basename_is_actually_in_the_repository`,
caused by an untracked `v2/research/rebase/p2/figures/data/` tree that appeared in this working copy
during the session and belongs to concurrent work; it is not touched here.
