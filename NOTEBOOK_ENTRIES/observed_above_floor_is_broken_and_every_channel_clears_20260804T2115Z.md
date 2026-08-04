# `observed_above_floor` is broken, and every channel clears its own floor once it is fixed

**Predeclaration:** `NOTEBOOK_ENTRIES/PREDECLARED_observed_above_floor_20260804T1843Z.md`, committed
before a single diagnostic case was constructed.
**Code:** `v2/research/rebase/nature/floor_flag_diagnostic.py` (constructed cases),
`v2/research/rebase/nature/floor_flag_recompute.py` (real artifacts).
**Evidence:** `runs/floor_flag_audit/floor_flag_diagnostic.json`,
`runs/floor_flag_audit/floor_flag_recompute.json`.
**Workspaces:** local `git -c core.autocrlf=false archive HEAD`, 648–654 files verified by git blob
SHA-1, 0 mismatches; box workspace at
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/runs_misc/floorflag_ws`, 650 files,
0 mismatches. Threads capped to 1 throughout. Every statistic imported from `v2/calibra/`.
**Suite:** `pytest morpheus/v2/tests morpheus/tests -q` on the verified workspace →
**534 passed**, 0 failed, 107s. Collection at the last pre-existing commit (`9b9d070`) is 517, so the
delta is exactly the 17 tests added by `v2/tests/test_inductive_adjustment.py` and no existing test
changed. (517 − 28 `test_p2_figures` tests = 489 on the matplotlib-less box venv; the brief's 476
baseline predates the tests added by the last few commits on this branch.)

---

## 0. Verdict

**`observed_above_floor` is broken, in two independent ways, and all 13 states of
`runs/calibra_v3_targeted` clear their own detection floor by 1.80–2.74× on a corrected comparator
that was already present, unread, in the shipped artifact.**

The project's standing explanation — that `observed_above_floor = 0` is *correct*, a units mismatch
between a single-random-direction floor and a 16-component multivariate maximum — does not survive
contact with cases whose answers are known in advance. The flag is not conservatively reporting a
comparison it cannot make. It is a constant.

---

## 1. What the flag actually computes

`calibration.py:238`:

```python
"observed_above_floor": bool(np.isfinite(matched) and np.isfinite(self.detection_floor)
                             and matched > self.detection_floor),
```

with `matched = meta["observed_matched_direction"]`, built in `spike_recovery_curve` as

```python
matched = _correlation(x_residual @ u, y_residual @ v)     # per draw, i == 0 branch
observed_matched = float(np.nanmedian([m for _, m in results]))
```

`(u, v)` is **the random direction pair drawn for that draw's spike**, and `_correlation` is
**signed**. The right-hand side, `detection_floor`, is always a positive grid level when finite.

### Defect 1 — a signed statistic against a positive threshold

Documented as deliberate at `calibration._correlation` (lines 147–162). The argument given there is
sound *for the paired spike comparison*, where a draw is compared against its own level-0 value and
the induced baseline's random sign would break the pairing under `|r|`. It does not transfer to the
**unpaired** comparison against a positive floor, where a magnitude is the only coherent statistic.

This is the same defect class this project already fixed once, mirrored.
`v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §0 defect 3 records a readout that took an
absolute value *before* pairing and destroyed the paired comparison. The repair correctly made the
paired statistic signed — and then carried the signed statistic into an unpaired comparison, where it
is wrong in the opposite direction. Both defects are the same underlying error: the sign convention
was chosen for one comparison and applied to the other.

### Defect 2 — a random direction is not the channel

For any multi-column representation, `corr(X_res u, Y_res v)` along a *random* pair carries
essentially none of the channel: the overlap of a random direction with the channel's direction falls
as `1/sqrt(d)` on each side. The comparator is therefore ≈ 0 whatever the channel's strength, and it
**cannot exceed a positive floor by construction**. It is worse than that: the `detection_floor` is
*defined* as the smallest level whose readout clears the level-0 p90, and `observed_matched_direction`
is drawn from approximately that same level-0 distribution. The flag compares a quantity to a
threshold that was calibrated to sit above it.

This is not an inference. It is measured in §2.

---

## 2. Constructed cases with known answers

`floor_flag_diagnostic.py`. The "real" association is planted with `spike_targets` **itself**, i.e.
with the identical construction the floor is calibrated against, so the data carries a
single-direction correlation of exactly the strength the instrument measures floors in and the units
argument cannot apply. `truth_planted_direction_abs` is `paired_absolute_correlation` along the
planted pair after the identical residualisation — the ground truth, on the floor's own scale.

All four predeclared criteria for **broken** are met.

### B1 — sign defect (the paired sign test)

Same cohort, same seed, same planted magnitude; only the sign of the association is reflected
(an exact reflection in the `v0` coordinate). One column per side, which is the case P4 unit
certification runs in.

| seed | floor | truth | truth/floor | shipped `matched` (pos) | flag (pos) | shipped `matched` (neg) | flag (neg) |
|---|---|---|---|---|---|---|---|
| 42 | 0.20 | 0.6955 | **3.48×** | −0.6955 | **False** | +0.6955 | True |
| 43 | 0.30 | 0.5556 | 1.85× | +0.5556 | True | −0.5556 | **False** |
| 44 | 0.60 | 0.4693 | 0.78× | +0.4693 | False | −0.4693 | False |
| 45 | 0.075 | 0.6940 | **9.25×** | +0.6940 | True | −0.6940 | **False** |

Three of the four pairs return **different verdicts for identical evidence**. The fourth agrees only
because it is below the floor in both. Nothing differs between the two rows of a pair except a sign
that no correlation-magnitude readout in the repository can even see —
`top_canonical_correlation`, `heldout_top_cca` and `paired_absolute_correlation` all take `abs`.

### B2 — power defect at 10× the floor

Predeclared: "a genuine association whose true single-direction correlation is ≥ 10× the
`detection_floor` the same call reports returns `False`."

| case | floor | truth | truth/floor | shipped flag |
|---|---|---|---|---|
| `J_mild_confound_rho0.8` (16×16) | 0.075 | **0.7897** | **10.53×** | **False** |
| `I_signpair_seed45_neg` (1×1) | 0.075 | 0.6940 | 9.25× | **False** |

Met exactly as predeclared.

### B3 — seed lottery

Eight one-column cases at a fixed planted magnitude: the sign of `observed_matched_direction` is
`±truth` at random (with `p = q = 1`, `u` and `v` are literally `±1` from a normal draw). Two of
eight seeds (42, 48) land negative, i.e. 25% of seeds, above the predeclared 20% bar. The flag is a
coin flip on data that has not changed.

### B4 — irrelevance

The 16×16 sweep, holding everything fixed but the planted strength:

| planted ρ | truth (floor units) | held-out top-CCA | shipped \|`matched`\| | shipped flag |
|---|---|---|---|---|
| 0.0 | 0.0271 | 0.0366 | 0.0017 | False |
| 0.2 | 0.1716 | 0.0497 | 0.0007 | False |
| 0.4 | 0.3739 | 0.2635 | 0.0022 | False |
| 0.6 | 0.5802 | 0.5493 | 0.0039 | False |
| 0.8 | 0.7897 | 0.7833 | 0.0013 | False |

The channel moves through a factor of 29 and the shipped comparator does not move at all — it sits at
0.001–0.004 throughout, three orders of magnitude below the thing it is supposed to measure. The
first `G_` family (strong confound) is the same picture: truth 0.012 → 0.627, comparator 0.0001 →
0.003.

### Guards on the corrected comparator, on the constructed cases

`G1` (pairing destroyed by permuting `y` within cancer strata): held-out top-CCA falls to
0.005–0.028 and never clears the floor, in every case, including those where the real channel reads
0.78. `G3`: the corrected comparator (0.26–0.78 at ρ ≥ 0.4) sits far above the
`confound_induced_baseline` (0.020). The corrected comparator is also *conservative* — at
`J rho=0.2` the truth (0.172) clears the floor and the held-out statistic (0.050) does not, because it
must find the direction from data.

---

## 3. Recomputation on the shipped artifacts — all 13 states

`floor_flag_recompute.py`, rebuilding `x`, `y` and the 99-column cancer+pooled-TSS design exactly as
`run_calibra.py` does (2,530 test patients, 90 non-control targets, 75 sites kept, `min_site_count`
10), against `runs/calibra_v3_targeted/task_rows.csv`.

**Reconstruction check first.** The rebuilt pipeline reproduces the shipped `heldout_top_cca` to a
maximum absolute difference of **9.4e-16** across all 13 states. The reconstruction is exact; nothing
below rests on a re-derivation that drifted.

**The corrected comparator was already in the artifact.** `heldout_top_cca` — canonical directions fit
on one half of the patients, correlation scored on the other half, `n_components = 16`, the same
budget the floor uses — has been emitted by `run_calibra.py` for every state since before this run.
It is a single-direction correlation between one image score and one molecular score after the
identical residualisation, which is the floor's own scale.

| method | state | floor | held-out top-CCA | **ratio** | corrected | sign-fix-only | G1 null max | G1 p | G4 min over 3 seeds |
|---|---|---|---|---|---|---|---|---|---|
| full | full_biology | 0.40 | 0.8757 | 2.19× | **PASS** | 0.2114 | 0.0747 | 0.0099 | 0.8757 |
| full | full_identity | 0.30 | 0.7792 | 2.60× | **PASS** | 0.0720 | 0.0732 | 0.0099 | 0.7792 |
| full | full_patient | 0.30 | 0.8235 | 2.74× | **PASS** | 0.1134 | 0.0875 | 0.0099 | 0.8235 |
| full | rna_biology | 0.50 | 0.8983 | 1.80× | **PASS** | 0.2337 | 0.0711 | 0.0099 | 0.8969 |
| full | rna_identity | 0.40 | 0.8460 | 2.11× | **PASS** | 0.0940 | 0.0914 | 0.0099 | 0.8310 |
| full | wsi_biology | 0.20 | 0.4768 | 2.38× | **PASS** | 0.0953 | 0.0894 | 0.0099 | 0.4768 |
| full | wsi_identity | 0.30 | 0.5393 | 1.80× | **PASS** | 0.0569 | 0.0730 | 0.0099 | 0.5393 |
| identity_only | full_identity | 0.30 | 0.7792 | 2.60× | **PASS** | 0.0720 | 0.0731 | 0.0099 | 0.7792 |
| identity_only | rna_identity | 0.40 | 0.8460 | 2.11× | **PASS** | 0.0940 | 0.0915 | 0.0099 | 0.8310 |
| identity_only | wsi_identity | 0.30 | 0.5393 | 1.80× | **PASS** | 0.0568 | 0.0730 | 0.0099 | 0.5393 |
| programme_only | full_biology | 0.40 | 0.8899 | 2.22× | **PASS** | 0.2305 | 0.0665 | 0.0099 | 0.8899 |
| programme_only | rna_biology | 0.50 | 0.9039 | 1.81× | **PASS** | 0.2506 | 0.0614 | 0.0099 | 0.9008 |
| programme_only | wsi_biology | 0.20 | 0.4748 | 2.37× | **PASS** | 0.1097 | 0.0680 | 0.0099 | 0.4748 |

**13 of 13 states clear their own floor. No state fails.**

### The two defects are separable, and both matter

The `sign-fix-only` column is the median `|corr|` along 40 random direction pairs — the shipped
comparator with its sign defect repaired and nothing else. It reads 0.057–0.251 and **clears the floor
in zero of 13 states**. Repairing the sign alone would not have changed a single verdict. Both defects
had to be found; either one alone is sufficient to make the flag useless on a multivariate channel.

### Predeclared guards, all cleared

* **G1 (destroyed pairing).** 100 permutations of `y` within cancer strata, per state. The permuted
  held-out top-CCA maxes out at **0.0614–0.0915**, below every floor (0.20–0.50) in every state; the
  observed value exceeds all 100 permutations in all 13 states (`p = 1/101 = 0.0099`, the minimum
  attainable). **This is the quantitative refutation of the capacity objection.** The in-sample
  `adjusted_top_cca` genuinely *is* capacity-inflated — its own permutation null p95 is 0.171 — which
  is why comparing *it* to the floor was rightly forbidden. The **held-out** statistic's null is
  0.06–0.09. The objection is real for one statistic and empirically false for the other.
* **G2 (no in-sample maximum).** `heldout_top_cca` fits directions on train rows and scores on rows it
  never saw. The in-sample `observed` is reported alongside but never used for the verdict.
* **G3 (above the induced baseline).** `confound_induced_baseline` is 0.065–0.140; the corrected
  comparator is 0.475–0.904, clearing it by 3.4–13×.
* **G4 (seed stability).** Three seeds moving both the residualiser folds and the CCA split: the
  minimum over seeds still clears in all 13 states, and the spread is ≤ 0.015.

One further point in the correction's favour, which was not predeclared and is noted because it
strengthens the case: the floor's own noise reference (`null_reference_p90`, 0.185 for `full_biology`)
is **larger** than the corrected comparator's permutation maximum (0.075). The floor is calibrated
against a noisier baseline than the comparator carries, so the comparison is conservative in the
direction that matters.

---

## 4. What would make me distrust this, stated as promised

The result runs in our favour, so the caveats are given the same prominence as the numbers.

1. **The floor is measured on random direction pairs; the channel's direction is fitted.** The
   induced level-0 correlation depends on how much of `(u, v)` lies in the design span, and a
   *fitted* direction may lie in the design span differently from a random one. I did **not** measure
   a direction-matched floor. The measurement that would close this is
   `spike_recovery_curve(..., molecular_directions=V)` with `V` the held-out CCA's molecular
   direction — though note `spike_targets` still draws `u` at random, so even that is only
   half-matched, and closing it properly needs an image-direction argument the function does not
   currently take. **This is the largest open hole in the favourable result and it should be filled
   before the corrected claim is load-bearing in a paper.**
2. **The ratio is suspiciously stable.** 1.80–2.74× across 13 states spanning channels of 0.47 to
   0.90. Floors and channels co-vary: `baseline_recovered_median` runs 0.067 (wsi_biology, channel
   0.475, floor 0.20) to 0.140 (rna_biology, channel 0.898, floor 0.50). A representation with more
   structure in the design span has both a larger induced baseline and a larger channel, so the floor
   and the thing it grades are not independent. The verdict is robust to this — a 2× margin is a 2×
   margin — but "the channel is 2.2× its floor" should not be read as an effect size.
3. **`heldout_top_cca` splits patients at random, not by site.** A held-out patient can share a
   tissue source site with training patients. The spike readout uses all rows too, so the two sides of
   the comparison are on the same footing, but a leave-sites-out version (the machinery exists at
   `v2/calibra/leave_sites_out.py`) would be stricter and is not run here.
4. **`n_draws = 40` and a 13-level grid** mean the floor resolves only to a grid level. A floor of
   0.40 could be anywhere in (0.30, 0.40].

None of these overturns the finding that the *flag* is broken, which rests on constructed cases with
known answers and does not depend on any real-data judgement.

---

## 5. A source-versus-source disagreement, reported as required

Two files in `v2/calibra/` state opposite things about whether a fitted-direction statistic is on the
floor's scale.

* `spectral.heldout_single_direction_correlation` docstring: *"The CALIBRA `detection_floor` is
  expressed in single-direction correlation units… Grading a per-target negative control against that
  floor therefore requires a per-target statistic on the same scale,"* and proceeds to define a
  fitted, out-of-fold direction as that statistic.
* `run_calibra.random_direction_column_correlation` docstring: *"Grading a fitted-direction readout
  against a random-direction floor is not a like-for-like comparison, and any per-target claim that
  does so is reading a floor that was never measured for it,"* and `run_calibra.py` tags every
  fitted-direction row `NOT_floor_units`.

The measurements in §2–3 favour the first: the capacity objection that motivates `NOT_floor_units` is
true of an *in-sample* maximum and empirically false of a *held-out* one (G1: null max 0.06–0.09
against floors of 0.20–0.50). But the second is not simply wrong either — its concern is caveat 1
above, which remains open. Both are recorded; neither file is edited here.

A third disagreement: `run_calibra.random_direction_column_correlation` is *also* signed
(`np.median` of signed per-draw products), so the per-target random-control grading in
`grade_random_controls` inherits defect 1. That grading compares those signed values against a
positive `detection_floor`, so its "controls sit below the floor" verdict is partly structural. **Not
investigated here — flagged for whoever owns Track 1.**

---

## 6. Prose corrections flagged, not made

Drafts are out of scope for edits by instruction. These passages are now contradicted by the
measurements above and by the shipped artifact's own `heldout_top_cca` column.

1. **`paper/P1_CALIBRA_DRAFT.md` §4.5 (≈ line 760)** — *"`observed_above_floor = 0` for every state,
   and that is the correct answer rather than a failure… **This repository therefore contains no
   measurement of the real channel, in the floor's own units, that exceeds the floor.**"* The second
   sentence is false: `heldout_top_cca` is such a measurement, it is in the artifact, and it exceeds
   the floor in 13 of 13 states. The first is false for the reason given in §1–2.
2. **`paper/P1_CALIBRA_DRAFT.md` §4.5** — *"Comparing the 0.47–0.61 channel against the 0.30–0.40
   floor directly was the original defect of §4.1 and must not be done."* This conflates two
   statistics. The original defect (`PHASE1B_TARGETED_READOUT.md` §0 defect 1) was a *spike readout*
   that was a maximum. Comparing the **in-sample** `adjusted_top_cca` to the floor is indeed
   forbidden (permutation null p95 0.171). Comparing the **held-out** statistic is not, and G1
   measures the difference.
3. **`paper/P1_CALIBRA_DRAFT.md` limitation 7 (≈ line 1368)** — *"There is no measurement in this
   repository of the real channel, in the floor's own units, that exceeds the floor. The paper
   therefore certifies pipeline sensitivity, not the significance of the channel."* **This is the
   disclaimer that reverses.** The corrected statement is that the instrument certifies pipeline
   sensitivity *and* places the observed channel 1.80–2.74× above that sensitivity limit, subject to
   caveat 1 of §4 (the floor is measured on random directions, the channel's is fitted). Significance
   still rests on the permutation null; what changes is that the floor now also speaks to the channel.
4. **`paper/P1_FIGURES.md` line 399** — the evidence-table row *"`observed_above_floor = 0` for every
   state, correctly… No such measurement exists"* — same correction; the row should carry
   `heldout_top_cca` vs `detection_floor` with the G1 null.
5. **`paper/P1_FIGURES.md` lines 110, 116** — `observed_matched_direction` −0.028 to +0.036 is quoted
   as evidence that "a random pair sees nothing, because the channel is concentrated in particular
   directions". The observation is right; the inference drawn from it (that no comparison is
   therefore possible) is not. The number is also a *signed median of random-sign draws*, so it
   measures approximately nothing at all — the magnitude version is 0.057–0.251 (§3).
6. **`v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md` line 203** — same "`observed_above_floor`
   = 0 for every state and that is the correct answer" claim.
7. **`paper/P2_RANK_DRAFT.md`** — the task brief states P2 carries this explanation. A grep of
   `paper/P2_RANK_DRAFT.md` for `observed_above_floor`, `single random direction`, `floor_scale` and
   `floor's own units` returns **no hits**, so no P2 passage is identified here. If the explanation is
   in P2 under different wording, it has not been located and should be searched again by whoever
   owns that draft.

Historical notebook entries (`t12_t14_t16_t17_calibra_ledger_20260803T0230Z.md`,
`p1_submission_draft_20260803T1230Z.md`) also carry the claim. They are the record of what was
believed at the time and are left as written; this entry supersedes them.

## 7. What was NOT changed

`v2/calibra/calibration.py` is **unmodified**. The flag is still broken in the library. Fixing it
changes the value of a key that appears in shipped artifacts and in `claim_guards.py`'s reach, and
`v2/tests/test_calibra.py:167` asserts only its type. The fix wanted is a `floor_scale`-aware
comparator that takes a magnitude and accepts a channel statistic from the caller rather than
inventing one from a random direction — that is a library change with an artifact-compatibility
question attached, and it is left for a decision rather than made unilaterally here.
