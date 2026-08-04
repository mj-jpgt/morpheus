# The direction-matched floor: the reversal is NOT confirmed, because the floor cannot be measured that way

**Predeclaration:** `NOTEBOOK_ENTRIES/PREDECLARED_direction_matched_floor_20260804T2230Z.md`, committed
(`acf4c3c`) before `spike_targets` had an `image_direction` argument and before any spike had ever
been planted along a fitted direction in this repository.
**Antecedent:** `NOTEBOOK_ENTRIES/observed_above_floor_is_broken_and_every_channel_clears_20260804T2115Z.md`,
which reported 13 of 13 states clearing their floor by 1.80–2.74× and named this as its own largest
open hole.
**Code:** `v2/research/rebase/nature/direction_matched_floor.py` (the 2×2),
`direction_matched_floor_controls.py` (H6 + mechanism), `track1_random_control_sign.py`.
**Evidence:** `runs_misc/dmfloor_out/dm_*.json` (13 states × 6 cells),
`runs_misc/h6_out/h6_*.json` (H6), `runs_misc/h6_out/mech_*.json` (mechanism),
`runs_misc/track1_sign/track1_random_control_sign.json`, all under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/`.
**Workspaces:** four, each `git -c core.autocrlf=false archive HEAD` verified per file by git blob
SHA-1 — 659, 662, 663 and 663 files, **0 mismatches** each. Threads capped to 1 throughout. Every
statistic imported from `v2/calibra/`; the floor rule is `calibration.floors_from_recovery`, the
instrument's own.
**Suite:** 540 passed, 0 failed (535 at the parent commit `7482a38`; the delta is exactly the 5 tests
added).

---

## 0. Verdict

**The 13-of-13 reversal is NOT confirmed by a direction-matched floor — and it is not refuted
either. A direction-matched floor cannot be measured with the current spike construction, and this
run identifies and quantifies why.** By the criteria fixed in advance, the primary cell is
**UNINFORMATIVE**: `U1` fires (the like-for-like floor fails to resolve on the grid in 6 of 13 states,
bar was ≥ 4) and `U2` fires hard (its level-0 baseline exceeds half the channel in 12 of 13 states,
bar was ≥ 4).

Neither predeclared survival criterion is met. `S1` fails (12 of 13 clear the direction-matched
oracle floor, not 13; minimum ratio 0.95, bar was ≥ 1.25). `S2` fails (the like-for-like floor
resolves in 7 of 13 states, bar was ≥ 11). Neither artifact criterion is met either (`A1`: 1 state of
13, bar ≥ 4; `A2`: 0 of the 7 resolved states, bar ≥ 4).

**The hole named by the antecedent entry stays open.** What has changed is that it now has a measured
cause and a named repair, and that the headline "13 of 13 at 1.80–2.74×" cannot be quoted as if the
floor it is measured against were the right one.

---

## 1. H6 — the blocking control, cleared

Declared blocking: if the extension moved any shipped floor, everything here is void. Re-running
`spike_recovery_curve` with **no direction arguments**, at the shipped protocol (13-level grid,
`n_draws = 40`, `k = 16`, `seed = 42`, 2,530 test patients, 99-column design):

**13 of 13 states reproduce the shipped `detection_floor` and `transmission_floor` exactly, and the
shipped signed `observed_matched_direction` to `atol = 1e-9`.** 0 failures. The refactor —
`image_direction`, `direction_pairs`, `floors_from_recovery`, and the shared `_heldout_cca_fit` —
did not change the instrument.

## 2. Is planting along a fitted direction ill-posed? Yes, and here is the number

Predeclared §1.2 said the risk in advance: `spike_targets` builds
`a_perp = standardise(a - rho*s)` with `rho = corr(s, a)` taken in **raw** space, so on a fitted pair
a residue can survive into residual space and inflate the level-0 baseline. Measured, per state
(median over 8 pairs):

| cohort | `rho_raw` fitted | `rho_raw` random | level-0 readout \|r\| fitted | level-0 readout \|r\| random |
|---|---|---|---|---|
| real | **0.590 – 0.889** | 0.089 – 0.337 | 0.082 – 0.302 | 0.066 – 0.126 |
| pairing destroyed | **0.189 – 0.348** | 0.072 – 0.151 | 0.190 – 0.405 | 0.063 – 0.151 |

A fitted pair's raw-space alignment is 2–9× a random pair's, so the quantity the construction
subtracts is 2–9× larger, and the residue it leaves behind is correspondingly larger. The signature
is visible in the recovery curves themselves: on a matched pair they are **non-monotone**. For
`full_biology`, the like-for-like curve runs

`0.450, 0.438, 0.425, 0.411, 0.381, 0.343, 0.301, 0.208, 0.100, 0.020, 0.181, 0.424, 0.615`

across `r_true = 0 … 0.6` — it *falls* by a factor of 22 before it rises. A stronger planted signal
reads *lower*, because the replacement is destroying pre-existing structure along the axis faster
than the spike restores it. That is the same pathology `PHASE1B_TARGETED_READOUT.md` §0 defect 1
recorded for a max readout, reappearing on a matched direction. **A curve that is not monotone in
`r_true` has no detection floor**, and the floor rule returns NaN or a meaningless grid level for it.

Crucially the effect survives destroying the pairing (`rho_raw` 0.19–0.35, level-0 0.19–0.41 vs the
random pair's 0.06–0.15). Since a within-cancer permutation preserves cancer structure on both sides,
the surviving alignment is **cancer-mediated, not channel-mediated** — so it is a property of the
instrument and the design, not of the thing being graded. Measured, not asserted.

**The nearest well-posed comparison and what it costs.** The two cells that come closest —
`F_match_perm` (matched pair, pairing destroyed, oracle readout) and `F_lfl_match` (matched pair,
pairing destroyed, the channel's own `heldout_top_cca` readout) — are the right *design*, and they
are still contaminated by the above in 4 of 13 and 12 of 13 states respectively. Their numbers are
reported in full below, and neither is quotable as *the* direction-matched floor.

**The repair, named not made.** Orthogonalise the planted component in **residual** space rather than
raw space (or plant additively on residuals), so that `corr(x_res u, resid(spiked) v) = 0` at
`r_true = 0` by construction for any pair. That changes every floor this project has ever measured
and is a separate, declared decision; it is not made here.

## 3. The 2×2, all 13 states, nothing dropped

Channel = the shipped `heldout_top_cca` (reproduced from raw data to ≤ 9.4e-16). Floors are the
UNPAIRED `detection_floor`. `P`/`F` = channel above/below that cell's floor.

| method | state | channel | shipped floor | `F_rand_real` (re-draw) | `F_match_real` | `F_rand_perm` | **`F_match_perm`** | `F_lfl_rand` | **`F_lfl_match`** |
|---|---|---|---|---|---|---|---|---|---|
| full | full_biology | 0.8757 | 0.40 | 0.40 | 0.075 | 0.50 | **0.40 (2.19× P)** | NaN | **0.60 (1.46× P)** |
| full | full_identity | 0.7792 | 0.30 | 0.30 | 0.30 | 0.40 | **0.40 (1.95× P)** | NaN | **NaN** |
| full | full_patient | 0.8235 | 0.30 | 0.30 | 0.075 | 0.40 | **0.30 (2.74× P)** | NaN | **0.40 (2.06× P)** |
| full | rna_biology | 0.8983 | 0.50 | 0.50 | 0.15 | 0.50 | **0.40 (2.25× P)** | NaN | **NaN** |
| full | rna_identity | 0.8460 | 0.40 | 0.40 | 0.15 | 0.40 | **0.30 (2.82× P)** | NaN | **0.60 (1.41× P)** |
| full | wsi_biology | 0.4768 | 0.20 | 0.40 | 0.30 | 0.40 | **0.50 (0.95× F)** | 0.60 (0.79× F) | **NaN** |
| full | wsi_identity | 0.5393 | 0.30 | 0.30 | 0.20 | 0.30 | **0.30 (1.80× P)** | NaN | **0.50 (1.08× P)** |
| identity_only | full_identity | 0.7792 | 0.30 | 0.30 | 0.30 | 0.40 | **0.40 (1.95× P)** | NaN | **NaN** |
| identity_only | rna_identity | 0.8460 | 0.40 | 0.40 | 0.15 | 0.40 | **0.30 (2.82× P)** | NaN | **0.60 (1.41× P)** |
| identity_only | wsi_identity | 0.5393 | 0.30 | 0.30 | 0.20 | 0.30 | **0.30 (1.80× P)** | NaN | **0.50 (1.08× P)** |
| programme_only | full_biology | 0.8899 | 0.40 | 0.40 | 0.075 | 0.50 | **0.40 (2.22× P)** | NaN | **0.60 (1.48× P)** |
| programme_only | rna_biology | 0.9039 | 0.50 | 0.40 | 0.15 | 0.50 | **0.40 (2.26× P)** | NaN | **NaN** |
| programme_only | wsi_biology | 0.4748 | 0.20 | 0.30 | 0.30 | 0.30 | **0.40 (1.19× P)** | 0.60 (0.79× F) | **NaN** |

* **`F_match_perm`** (matched, geometry-clean, oracle readout): **12 of 13 clear**, ratios 0.95–2.82,
  median 2.19. `full::wsi_biology` **fails** at 0.95×, and `programme_only::wsi_biology` is marginal
  at 1.19×. Stable to the permutation: two independent within-cancer permutations reproduce every
  one of the 13 floors exactly.
* **`F_lfl_match`** (the primary cell, the channel's own statistic on both sides): resolves in **7 of
  13** states, all 7 clearing at **1.08–2.06×, median 1.41** — roughly **half** the antecedent's
  1.80–2.74×. It is NaN in the other 6. Permutation-stable where it resolves (one state moves
  0.5↔0.6).
* **`F_lfl_rand`** (random pair, the channel's readout): **NaN in 11 of 13**, and in the two states
  where it resolves the channel *fails* it. This is the sharpest single result in the run and it is
  reported in §5.
* **`F_match_real`** (matched pair, real cohort): floors as low as 0.075, giving ratios up to
  **11.9×**. This is the number that would have been reported had the pairing-destroyed control not
  been predeclared. It is **withdrawn** — see H3/H5.

## 4. Every guard, fired or not

| guard | outcome | value |
|---|---|---|
| **H1** — a floor that conveniently drops | **does not fire** | `F_match_perm` is never more than one grid level below `F_rand_perm` in any state |
| **H2** — degeneracy from pinning | **does not fire** | no cell's detection floor sits at the finest grid level; floors span 0.075–0.6. (`transmission_floor` is 0.01 in every oracle cell, the known paired degeneracy — never quoted as a detection limit) |
| **H3** — real vs destroyed cohort disagree | **FIRES** | the two differ by more than one grid level in **8 of 13** states (e.g. 0.075 vs 0.40). Both reported; neither quoted alone; attributed to the raw-space residue of §2, now measured |
| **H4** — S1 without S2 is only PARTIAL | n/a | neither S1 nor S2 holds, so not even PARTIAL is claimable |
| **H5** — favourable only on the circular arm | **FIRES** | the matched floor collapses to 0.075–0.30 on the cohort whose channel selected the direction, and returns to 0.30–0.50 once the pairing is destroyed. The 11.9× is withdrawn as declared |
| **H6** — refactor drift (blocking) | **CLEARS** | 13/13 exact on `detection_floor`, `transmission_floor` and `observed_matched_direction` |
| **U1** — floor unresolved on the grid | **FIRES** | `F_lfl_match` NaN in 6/13; `F_lfl_rand` NaN in 11/13 |
| **U2** — level-0 baseline > 0.5 × channel | **FIRES** | `F_lfl_match` in 12/13 (0.45–0.60 against channels of 0.47–0.90); `F_match_perm` in 4/13 |
| **U3** — floor pinned at the grid edge | **does not fire** | — |

### Carried forward unresolved, as required

**H5' — the ratio is still suspiciously stable, and it is still not an effect size.** `F_match_perm`
gives 0.95–2.82× across channels spanning 0.47–0.90; `F_lfl_match` gives 1.08–2.06×. Floors and
channels co-vary through ambient structure — a representation with more structure in the design span
has both a larger induced baseline and a larger channel — so the band is not evidence of a common
effect size and is not upgraded to one.

**New, and against us: the shipped floor is itself unstable to the direction draw.** `F_rand_real`
is the *same* protocol as the shipped floor with the random pairs re-drawn (H6 shows the code path is
identical). It moves by one to two grid levels in **3 of 13** states: `full::wsi_biology` 0.20 → 0.40,
`programme_only::wsi_biology` 0.20 → 0.30, `programme_only::rna_biology` 0.50 → 0.40. So "the channel
is 1.80–2.74× its floor" was already resting on a denominator carrying about ±1 grid level of draw
noise, before any of the direction-matching question is reached.

## 5. Which docstring was right

Both are quoted, as required.

* `spectral.heldout_single_direction_correlation`: *"The CALIBRA `detection_floor` is expressed in
  single-direction correlation units… Grading a per-target negative control against that floor
  therefore requires a per-target statistic on the same scale"* — and it defines a fitted, out-of-fold
  direction as that statistic.
* `run_calibra.random_direction_column_correlation`: *"Grading a fitted-direction readout against a
  random-direction floor is not a like-for-like comparison, and any per-target claim that does so is
  reading a floor that was never measured for it"* — and it tags such rows `NOT_floor_units`.

The predeclared reading rule fixed three cases. Case 1 (the two floors agree within one grid level in
≥ 11 states) is **not** met: where `F_lfl_match` resolves it sits **one to two grid levels above** the
shipped random-direction floor in 7 of 7. Case 3 (the like-for-like floor exceeds the channel) is not
met either. **Case 2 holds: `run_calibra.random_direction_column_correlation` is right that the
shipped comparison is not like-for-like, and its `NOT_floor_units` tag is correct caution.**

The sharpest evidence is `F_lfl_rand`. Planting a spike on a **random** pair and reading it with the
**channel's** statistic gives no detectable floor at all in **11 of 13** states, even at `r_true = 0.6`
— and where it does resolve, at 0.6, the channel does not clear it. A random-direction floor is
therefore not merely a different number from the fitted-direction floor; on the channel's own
readout it is largely **not measurable**. `heldout_single_direction_correlation`'s claim that a
fitted out-of-fold direction *is* on the floor's scale is **not supported** by this run.

What `random_direction_column_correlation` is *wrong* about is the implied conclusion that no
comparison is possible. Where a like-for-like floor resolves, the channel does clear it, at
1.08–2.06×. The correct position is: **the comparison is possible in principle, is not the one the
shipped artifacts make, and cannot be completed until the spike construction is repaired.**

**Both docstrings are therefore corrected in source** (they are library files, not drafts) rather than
one being declared the winner: `heldout_single_direction_correlation` gains the caveat that the
shipped floor is measured on a random pair and is not that statistic's floor;
`random_direction_column_correlation` gains a pointer to this measurement and drops the implication
that the comparison is impossible.

## 6. Track 1's random control: the sign fix does not change the verdict, and the control is weak

The antecedent flagged this and left it. Measured on all 13 states, both conventions, re-graded with
`grade_random_controls` (90 control columns, ceiling 5%):

| convention | control median | control p95 | control max | exceedance fraction | T1.4 passes |
|---|---|---|---|---|---|
| schema 1 (signed) | −0.150 … +0.032 | 0.012 – 0.072 | 0.025 – 0.147 | 0.000 in 13/13 | 13/13 |
| schema 2 (abs) | 0.036 – 0.220 | 0.049 – 0.399 | 0.053 – 0.439 | 0.000 in 12/13, 0.022 in 1 | 13/13 |

**0 of 13 verdicts change.** The fix raises the controls by up to 6× (`full_identity` 0.006 → 0.058;
`rna_biology` −0.150 → 0.220) and one state gains exceedances (`programme_only::full_biology`, 2.2%,
under the 5% ceiling). Under schema 1 the control median was **negative** in 10 of 13 states, so it
cleared a positive floor without reference to anything the controls contained — the pass was
structural, exactly as suspected.

**The more important finding is that the control is uninformative either way.** The *real* target
block, scored through the identical statistic, reads essentially the same as the controls: 0.218 vs
0.191 (`full_biology`), 0.228 vs 0.220 (`rna_biology`), 0.040 vs 0.036 (`wsi_identity`). A statistic
that cannot separate 90 real RNA programme targets from 90 matched random gene sets is not testing
the controls; it is reporting the weakness of a random direction. "Controls sit below the floor"
should not be quoted as a negative control passing until it is shown that the same statistic puts
real targets somewhere else.

## 7. Consequences for the antecedent's claims

The antecedent's arithmetic is unchanged and its finding that the *flag* was broken stands — that
rests on constructed cases with known answers and is independent of everything here. What does not
survive is the **grading**:

1. **"13 of 13 states clear their own floor by 1.80–2.74×" must not be quoted as a like-for-like
   result.** The floor in that sentence is a random-direction floor, the channel's direction is
   fitted, and §5 shows those are not the same scale.
2. The best-supported restatements available today are: **12 of 13 at 0.95–2.82× against a
   direction-matched oracle floor on a pairing-destroyed cohort** (one state, `full::wsi_biology`,
   fails), and **7 of 13 at 1.08–2.06× where a like-for-like floor resolves at all**. Both carry the
   U2 contamination of §2.
3. The antecedent's proposed rewrite of `P1_CALIBRA_DRAFT.md` limitation 7 — that the instrument
   *"places the observed channel 1.80–2.74× above that sensitivity limit"* — **should not be made in
   that form.** The defensible replacement is that the flag was broken and is fixed, that the
   repository does now contain a fitted-direction channel statistic exceeding the random-direction
   floor in 13 of 13 states, and that a **direction-matched floor cannot yet be measured**, so
   whether the channel clears the floor that actually applies to it remains open. The other five
   prose corrections the antecedent flagged are unaffected.

## 8. Prose corrections flagged, not made

Drafts are out of scope by instruction.

* `paper/P1_CALIBRA_DRAFT.md` limitation 7 — see §7.3. Do **not** adopt the antecedent's proposed
  wording; adopt the restricted version.
* Anywhere the 1.80–2.74× range appears, it must carry (a) that the floor is measured on random
  direction pairs, (b) that a re-draw of those pairs moves it by one to two grid levels in 3 of 13
  states, and (c) that where a like-for-like floor resolves the margin roughly halves.
* Anywhere `random_control` T1.4 is described as a passing negative control, it must carry §6: the
  statistic does not separate real targets from matched random ones.

## 9. What was NOT done

* The spike construction was **not** repaired. Residual-space orthogonalisation is the named fix and
  it would move every floor in the project; that is a declared decision, not a side effect.
* No leave-sites-out version of any floor was run (`v2/calibra/leave_sites_out.py` exists); the
  held-out CCA split is still uniformly random, so a held-out patient can share a site with training
  patients. Carried forward from the antecedent, still open.
* `run_calibra.py` was **not** re-run over the CALIBRA cohort; no shipped artifact was rewritten.
* `NOTEBOOK.md`, the paper drafts and `claim_guards.py` are untouched by instruction.
