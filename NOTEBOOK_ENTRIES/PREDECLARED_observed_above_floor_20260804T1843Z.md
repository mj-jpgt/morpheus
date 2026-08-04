# PREDECLARED — is `observed_above_floor` broken?

**Written before running a single diagnostic case.** Source read (`v2/calibra/calibration.py`,
`v2/calibra/spectral.py`, `v2/calibra/residualise.py`) is complete; no synthetic case has been
constructed and no artifact has been recomputed at the time of writing. The criteria below are
fixed from here.

## 0. What the source says the flag is

`SpikeRecoveryResult.summary()` (calibration.py:238):

```python
"observed_above_floor": bool(np.isfinite(matched) and np.isfinite(self.detection_floor)
                             and matched > self.detection_floor),
```

where `matched = meta["observed_matched_direction"]`, produced in `spike_recovery_curve` as

```python
matched = _correlation(x_residual @ u, y_residual @ v)      # i == 0 branch, per draw
observed_matched = float(np.nanmedian([m for _, m in results]))
```

and `(u, v)` is the **random direction pair drawn for that draw's spike** — `u ~ N(0, I_p)` normalised,
`v ~ N(0, I_q)` normalised (or a supplied programme loading when `molecular_directions` is passed).
`_correlation` is **signed** by deliberate design (calibration.py:147-162), for a reason that applies
to the *paired* spike comparison and is documented there.

`detection_floor` is by construction a **non-negative level from the `levels` grid** (smallest
`r_true > 0` clearing the level-0 p90 in ≥ `recovery_fraction` of draws), so the right-hand side of the
comparison is always positive when finite.

## 1. Predeclared criteria for BROKEN

The flag is broken if **any** of the following holds. Each is a property a floor comparison must have
independent of the units argument.

* **B1 — sign defect.** Two datasets identical except for the sign of the planted association
  (`+rho` vs `-rho`, same magnitude, same everything else) return different values of the flag.
  A detection floor is a magnitude threshold; a real channel is not less real for being negatively
  oriented, and neither `top_canonical_correlation` nor `heldout_top_cca` (both of which take
  `abs`, spectral.py:`paired_absolute_correlation`) can distinguish the two.
* **B2 — power defect.** A dataset carrying a genuine association whose true single-direction
  correlation is **≥ 10× the `detection_floor` the same call reports** returns `False`.
* **B3 — seed lottery.** For fixed data with a fixed genuine association well above the floor, the
  flag changes when only `seed` (the direction draw) changes, in ≥ 20% of seeds. A certificate that
  is a coin flip is not a certificate.
* **B4 — irrelevance.** `observed_matched_direction` does not increase with the true strength of the
  association in the unspiked data: i.e. sweeping the real channel's strength from 0 to strong leaves
  the comparator statistically flat. That would mean the comparator measures something other than the
  channel and no threshold on it can ever be informative.

## 2. Predeclared criteria for NOT BROKEN

`observed_above_floor` is sound if, across both signs and ≥ 5 seeds:
`True` whenever the true single-direction association magnitude exceeds the reported
`detection_floor`, `False` whenever it is below, with disagreement confined to a boundary band of
±1 grid level. In that case the "units mismatch" explanation in `P1_CALIBRA_DRAFT` and
`P2_RANK_DRAFT` stands as written and nothing is rewritten.

## 3. If it is broken: predeclared corrected statistic

The floor is in **single-direction correlation units**, so the comparator must be a single-direction
correlation on the same data, and it must be honest about direction selection. The corrected
comparator is `spectral.heldout_single_direction_correlation` / `spectral.heldout_top_cca_indexed`
(magnitude), both of which are already in the library — **no statistic will be written inline**. The
corrected flag is `abs(comparator) > detection_floor`.

## 4. Predeclared guards — what would make me distrust a FAVOURABLE result

This correction, if it lands, runs **in our favour**: it would retire P1's "the instrument certifies
pipeline sensitivity but never the channel" limitation, one of the paper's most damaging self-
assessments. It is therefore held to the standard of a result going against us. A PASS is reported
only if all four guards clear, and each guard's value is reported whether or not it clears.

* **G1 — destroyed-pairing control.** The corrected comparator must NOT clear the floor when the
  cross-modal pairing is destroyed by permuting rows of `y` within cancer strata
  (`calibration.permutation_null`'s stratification). If a permuted cohort also PASSes, the
  comparator is capacity-inflated and the PASS is an artifact of direction fitting, not a channel.
* **G2 — no in-sample maximum.** The comparator must be fit out of fold / on held-out rows. If a
  PASS requires the in-sample multivariate `top_canonical_correlation` (`observed`), that is the
  ORIGINAL scale error re-committed under a new name and the PASS is withdrawn.
* **G3 — above the induced baseline.** The corrected comparator must exceed
  `meta["confound_induced_baseline"]` (the level-0 magnitude, 0.067–0.140 on the 99-column design).
  A comparator below the correlation that residualisation *manufactures* between orthogonal signals
  is not evidence of a channel, even if it exceeds a grid level.
* **G4 — seed and direction stability.** The verdict must hold at ≥ 3 seeds and must not depend on a
  post-hoc choice of which direction to report.

A PASS that clears the grid level but fails G1 or G3 will be reported as **NO PASS**, with the
numbers shown.

## 5. Predeclared reporting rule

Whatever the outcome:

* If a source disagrees with a document, both are reported.
* If the flag is broken, the recomputation is run on **every** state present in the existing
  artifacts, and states that still do not clear are reported with the same prominence as states that
  do. No state is dropped for being unfavourable.
* Prose corrections to `paper/P1_CALIBRA_DRAFT.md`, `paper/P1_FIGURES.md`,
  `paper/P2_RANK_DRAFT.md` and `v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md` are
  **flagged, not made** (drafts are out of scope for edits by instruction).
