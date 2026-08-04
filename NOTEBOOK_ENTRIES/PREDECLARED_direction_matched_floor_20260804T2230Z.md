# PREDECLARED — does the 13-of-13 reversal survive a DIRECTION-MATCHED floor?

**Written before any direction-matched floor has been measured.** At the time of writing:
`v2/calibra/calibration.py`, `spectral.py`, `residualise.py`, `run_calibra.py` and
`v2/research/rebase/nature/floor_flag_recompute.py` have been read in full; `spike_targets` has
**not** yet been given an image-direction argument; no spike has ever been planted along a fitted
direction anywhere in this repository. The criteria below are fixed from here and are not revised
after seeing a number.

Antecedent: `NOTEBOOK_ENTRIES/observed_above_floor_is_broken_and_every_channel_clears_20260804T2115Z.md`
(and its own predeclaration `PREDECLARED_observed_above_floor_20260804T1843Z.md`). That entry reports
13 of 13 states clearing their detection floor by 1.80–2.74× on `heldout_top_cca`, and names as its
**own largest open hole** the fact that the floor is measured on **random** direction pairs while the
channel's direction is **fitted**. This predeclaration governs the measurement that closes, or fails
to close, that hole. **This is the second favourable result in a row on this thread, so the bar is
set here, in writing, before the instrument exists.**

---

## 0. The hole, stated precisely

`spike_targets(x, y, r_true, rng=..., molecular_direction=None)` draws the image direction `u` from
`rng` and never accepts one from the caller. `spike_recovery_curve` can therefore pin the *molecular*
direction (via `molecular_directions`) but never the *image* direction. Consequently every
`detection_floor` in every shipped artifact is measured on a direction pair that is random on at
least one side, and in the shipped runs (`molecular_directions` unset) random on **both** sides.

The channel it is graded against, `heldout_top_cca`, is a correlation between an image direction and
a molecular direction that were both **fitted** to the data (out of fold, on train rows, scored on
held-out rows).

Two distinct mismatches are bundled inside "random vs fitted", and they must be separated because
they push in **opposite** directions:

* **M1 — direction geometry.** The confound-induced level-0 correlation depends on how much of
  `(u, v)` lies in the span of the 99-column design. A fitted direction need not have the same
  overlap as a random one, so its floor need not be the random-direction floor. Sign of the effect:
  **unknown in advance**, and it is exactly what is measured here.
* **M2 — oracle vs estimated readout.** The floor's readout *knows* `(u, v)`; the channel's readout
  must *estimate* them from data and pays a direction-estimation cost. An oracle detects a smaller
  planted signal than an estimator does, so an oracle floor is a **lower bound** on the floor that
  applies to the channel's own statistic. M2 therefore runs **anti-conservatively** — in our favour —
  and fixing M1 alone does **not** close the hole.

Fixing M1 alone would be the easy, self-flattering move. It is declared here in advance as
**insufficient**.

## 1. Is planting along a fitted direction ill-posed? — declared position, before measuring

It is **well-posed** as a geometry question and **incompletely posed** as a like-for-like question,
for three reasons that are stated now so that they cannot be discovered conveniently later:

1. **Circularity.** A direction fitted on the same rows the floor is then measured on is not an
   arbitrary axis: it was selected *because* it carries the channel. Mitigation declared in advance:
   the pair is fitted on a random half of patients (`heldout_cca_projection`'s train split), and the
   primary floors are measured on a cohort whose cross-modal pairing has been **destroyed** (rows of
   `y` permuted within cancer strata, exactly `permutation_null`'s stratification), which preserves
   the design geometry and the marginal covariance of both modalities while removing the channel that
   selected the direction.
2. **The replacement construction is orthogonalised in RAW space, not residual space.**
   `spike_targets` builds `a_new = r*s + sqrt(1-r^2)*a_perp` with `a_perp` orthogonal to
   `s = standardise(x @ u)` in raw coordinates. On a *fitted* pair the raw and residualised
   correlations differ, so at `r_true = 0` a residue of the pre-existing channel can survive along
   the planted axis and inflate the level-0 baseline. This is a property of the instrument, not of
   the channel. Declared discriminator: the same floor measured on the pairing-destroyed cohort,
   where there is no pre-existing channel to leak.
3. **M2 is not a direction-matching problem at all.** No choice of `(u, v)` makes an oracle readout
   into an estimated one.

**Nearest well-posed comparison, declared in advance:** plant the spike into a cohort whose
cross-modal pairing has been destroyed, along the state's own fitted direction pair, and score it
with **the channel's own statistic** (`spectral.heldout_top_cca`, `n_components = 16`, directions
refit on train rows and scored on held-out rows). Both sides of the comparison then use the same
statistic, the same component budget, the same residualisation and a matched direction. That is the
`F_lfl_match` cell below and it is the **primary** result of this entry.

## 2. The 2×2 (plus two references) that will be measured

Every cell uses the shipped protocol: the 13-level grid
`0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6`, `n_draws = 40`,
`n_components = 16`, `recovery_fraction = 0.8`, `seed = 42`, the 2,530 test patients, 90 non-control
targets, the 99-column cancer + pooled-TSS design (`min_site_count = 10`) — i.e. exactly
`runs/calibra_v3_targeted/calibra_protocol.json`. All 13 states are run. **Floors are read with the
UNPAIRED rule** (`detection_floor`: smallest level clearing the level-0 p90 in ≥ 80% of draws), which
is the quotable one; `transmission_floor` is recorded but never quoted as a detection limit.

| cell | direction pair | cohort | readout |
|---|---|---|---|
| `F_rand_real` | random `(u, v)` per draw | real | oracle: `corr(X_res u, Y_spiked_res v)` |
| `F_match_real` | fitted `(u, v)` per draw | real | oracle |
| `F_rand_perm` | random `(u, v)` per draw | pairing destroyed | oracle |
| `F_match_perm` | fitted `(u, v)` per draw | pairing destroyed | oracle |
| `F_lfl_rand` | random `(u, v)` per draw | pairing destroyed | `heldout_top_cca` (fitted, k=16) |
| `F_lfl_match` | fitted `(u, v)` per draw | pairing destroyed | `heldout_top_cca` (fitted, k=16) |

Draw-to-draw variability, which the unpaired floor rule requires: draw *k* uses split seed
`42 + k` for the CCA fit that produces its direction pair and for its within-cancer permutation
order. A pinned direction with no other randomness would give a degenerate, noiseless floor; that
degeneracy is itself a declared failure mode (H2 below).

Graded quantity: the shipped `heldout_top_cca` per state (0.4748–0.9039), which the antecedent entry
reproduced from raw data to 9.4e-16.

## 3. Predeclared criteria

Let `channel(s)` be the state's `heldout_top_cca` and `ratio_C(s) = channel(s) / F_C(s)`.

**THE REVERSAL SURVIVES** iff **both** of:

* **S1** — `F_match_perm` is finite for ≥ 11 of 13 states and `channel(s) > F_match_perm(s)` in
  **13 of 13**, with `min_s ratio_match_perm(s) >= 1.25`; and
* **S2** — `F_lfl_match` is finite for ≥ 11 of 13 states and `channel(s) > F_lfl_match(s)` in
  **13 of 13**, with `min_s ratio_lfl_match(s) >= 1.25`.

S2 is the load-bearing one. S1 without S2 is **PARTIAL**, not a survival.

**IT WAS AN ARTIFACT OF RANDOM-vs-FITTED** iff either:

* **A1** — `F_match_perm(s) >= channel(s)` in ≥ 4 of 13 states, or the median of
  `ratio_match_perm` falls below 1.0; or
* **A2** — `F_lfl_match(s) >= channel(s)` in ≥ 4 of 13 states, or the median of `ratio_lfl_match`
  falls below 1.0.

In that case the correct reading is that the fitted direction sits where the pipeline's floor is
higher, the 1.80–2.74× margins were bought against a floor measured somewhere easier, and P1's
limitation 7 should **not** be rewritten as the antecedent entry proposes.

**PARTIAL** — anything between: some states clear and some do not. Then every state is reported with
equal prominence, the claim is restricted to the states that clear on `F_lfl_match`, and the
antecedent's "13 of 13" headline is withdrawn in favour of the restricted count.

**UNINFORMATIVE** if any of:

* **U1** — a floor fails to resolve on the grid (NaN: even `r_true = 0.6` does not clear its own
  level-0 p90 in 80% of draws) in ≥ 4 of 13 states, for the cell in question;
* **U2** — the level-0 baseline of a cell exceeds `0.5 * channel(s)` in ≥ 4 states, meaning the spike
  construction cannot null the pre-existing structure along that axis and the cell measures leakage
  rather than a floor;
* **U3** — a cell's floor sits at the **finest non-zero grid level** (0.01) in ≥ 11 of 13 states. A
  floor pinned to the edge of the grid is a bound, not a measurement, and no ratio computed against
  it is quotable.

## 4. What would make me distrust a FAVOURABLE outcome — declared in advance

This is the second favourable result in a row on this thread. Each of these is reported with its
value whether or not it fires.

* **H1 — a floor that conveniently *drops*.** If `F_match_perm` is lower than `F_rand_perm` by more
  than one grid level in ≥ 7 of 13 states, that is a *larger* margin obtained by changing the
  measurement, and it is treated as suspect until the mechanism is shown. Declared discriminator: the
  matched-pair level-0 baseline (`baseline_recovered_median`) must be correspondingly smaller than
  the random-pair one, i.e. the fitted direction must be shown to lie *further out* of the design
  span. If the floor drops without the baseline dropping, the drop is an artifact and the result is
  reported as **NOT ESTABLISHED**.
* **H2 — degeneracy from pinning.** If any cell's floor equals the finest non-zero grid level in
  ≥ 11 of 13 states, the readout has lost its noise (see U3). A 50× margin against a 0.01 floor is
  not evidence; it is a grid artifact.
* **H3 — real-cohort vs destroyed-cohort disagreement.** If `F_match_real` and `F_match_perm` differ
  by more than one grid level in ≥ 7 states, both are reported and neither is quoted alone; the
  difference is attributed to the raw-space orthogonalisation residue of §1.2 and that is stated as
  an open instrument defect.
* **H4 — S1 without S2.** Clearing the oracle-readout matched floor while failing, or not measuring,
  the like-for-like floor is **PARTIAL**. M2 stays open and the hole is not closed. Declared now so
  that a cheap PASS cannot be presented as the full one.
* **H5 — the ratio remains suspiciously stable.** Carried forward unresolved from the antecedent
  entry: floors and channels co-vary through ambient structure, so a band of ratios across channels
  of 0.47–0.90 is **not** an effect size and will not be upgraded to one, however favourable. If the
  new ratios again land in a narrow band, that is restated as a caveat, not as corroboration.
* **H6 — refactor drift.** The extension adds an `image_direction` argument to `spike_targets` and a
  `direction_pairs` argument to `spike_recovery_curve`, and lifts the floor rule into a reusable
  function. **Declared blocking control:** re-running `F_rand_real` through the extended code must
  reproduce the shipped `detection_floor` (0.20–0.50) and `transmission_floor` for **all 13 states
  exactly**. If any shipped floor moves, the refactor changed the instrument and every number in this
  entry is void; I stop and report the drift instead of the result.

## 5. What this settles about the two docstrings

`spectral.heldout_single_direction_correlation` says a fitted out-of-fold direction *is* a statistic
on the floor's scale. `run_calibra.random_direction_column_correlation` says grading a fitted
direction against a random-direction floor "is not a like-for-like comparison" and tags such rows
`NOT_floor_units`. Declared reading rule, fixed before the numbers:

* If `F_lfl_match` ≈ `F_rand_real` (within one grid level in ≥ 11 of 13 states), the random-direction
  floor **was** on the fitted statistic's scale after all, and `heldout_single_direction_correlation`
  is right while `NOT_floor_units` is over-cautious.
* If `F_lfl_match` differs materially from `F_rand_real` but the channel still clears it,
  `random_direction_column_correlation` is right that the shipped comparison was not like-for-like,
  and wrong that no comparison is possible; the like-for-like floor is the one to quote and the
  shipped floor should not be.
* If `F_lfl_match` exceeds the channel, `random_direction_column_correlation` is right on both counts
  and `heldout_single_direction_correlation`'s docstring must be corrected.

Whatever the outcome, both docstrings are quoted in the result entry and the losing one is corrected
in source (both are library files, not drafts).

## 6. Reporting rules

* Every state is reported. No state is dropped for being unfavourable.
* Every guard H1–H6 is reported with its value, fired or not.
* Both sides of any source-vs-source disagreement are quoted.
* Nothing is computed inline: every statistic comes from `v2/calibra/`. The floor rule itself is
  lifted into `calibration.floors_from_recovery` so that the like-for-like cell uses the **same**
  floor rule as the shipped instrument rather than a re-implementation.
* Workspace verified per-file by git blob SHA-1 from `git -c core.autocrlf=false archive HEAD`;
  threads capped to 1; persistent storage only.
