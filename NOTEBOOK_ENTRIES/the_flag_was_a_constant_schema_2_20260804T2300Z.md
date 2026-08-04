# Schema 2: `observed_above_floor` was a constant, not a measurement — fixed, with provenance

**Scope:** the library fix deliberately deferred by
`NOTEBOOK_ENTRIES/observed_above_floor_is_broken_and_every_channel_clears_20260804T2115Z.md` §7,
which left `v2/calibra/calibration.py` unmodified because *"fixing it changes the value of a key that
appears in shipped artifacts"*. It is fixed here, and the artifact-compatibility question is answered
explicitly rather than absorbed.

**Suite:** `pytest morpheus/v2/tests morpheus/tests -q` on a workspace built by
`git -c core.autocrlf=false archive HEAD` → **540 passed, 0 failed**, 108 s. Collection at the parent
commit (`7482a38`) is **535**, so the delta is exactly the **5 tests added here** and no existing test
was removed or weakened. (`test_p2_figures` needs matplotlib, absent from the box venv, so the box
count is lower by those 28; nothing was installed into `~/venv`.)

---

## 0. What was wrong, in one line each

Two keys in the emitted schema were not measurements.

1. **`observed_above_floor`** compared the **signed** `observed_matched_direction` — a correlation
   along a **random** direction pair — against a strictly positive `detection_floor`. On any
   multi-column representation a random pair carries essentially none of the channel, so the
   left-hand side sits near zero however strong the channel is; the random sign then makes the
   verdict a coin flip on top of that. The flag was `0` for every state of every shipped run because
   it could not be anything else.
2. **`run_calibra.random_direction_column_correlation`** returned `np.median(signed_scores, axis=0)`
   and was likewise graded against a positive floor by `grade_random_controls`. With `u` random the
   per-draw sign is random, so the signed median collapses towards zero for **every** column whatever
   that column carries — and a control that always reads ~0 always passes. Track 1's T1.4 negative
   control was therefore passing partly by construction.

Both are the same error, and it is the *mirror* of one this project already fixed once
(`v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §0 defect 3): a sign convention chosen for
the **paired** within-draw spike comparison, where it is correct and documented at
`calibration._correlation`, carried into an **unpaired** comparison against a positive threshold,
where only a magnitude is coherent.

## 1. The fix

`v2/calibra/calibration.py`:

* new module constant **`SUMMARY_SCHEMA_VERSION = 2`**, documented in place with what schema 1 emitted
  and why it must be ignored;
* new function **`channel_clears_floor(channel_statistic, detection_floor) -> (verdict, status)`**.
  It grades on **magnitude**, and it returns **`None`** — never a silent `False` — when the comparison
  cannot be made. "We could not grade this" and "the channel is below its floor" are opposite
  scientific statements and schema 1 emitted the same value for both;
* `SpikeRecoveryResult` gains **`channel_statistic`** and `channel_statistic_name`, supplied by the
  **caller**. This module has no way to measure a channel — its own matched readout is a random
  direction — and inventing the comparator internally is precisely what made the flag a constant;
* `summary()` emits `summary_schema_version`, `observed_above_floor` (bool or `None`),
  `observed_above_floor_status`, `channel_statistic`, `channel_statistic_name`, and
  `observed_matched_direction_abs` beside the retained signed `observed_matched_direction`;
* the module docstring's **"Scale warning"** is corrected. The prohibition on comparing the
  *in-sample* `top_canonical_correlation` to the floor stands unchanged (its within-strata
  permutation null reaches p95 = 0.171). What was wrong was the following sentence, which offered
  `observed_matched_direction` as "the same-units comparator".

`v2/calibra/run_calibra.py`:

* `_channel_measurement` computes `heldout_top_cca` **before** `summary()` and hands it to the result
  as the channel statistic, so `observed_above_floor` is now a real verdict on the fitted-direction,
  out-of-fold statistic. It is emitted as `NaN` when ungradable, never `0.0`;
* `random_direction_column_correlation` now returns `np.median(np.abs(scores), axis=0)`;
* `calibra_protocol.json` gains **`summary_schema_version`** and a `schema_2_changes` list naming
  every key whose *meaning* changed.

`v2/calibra/spectral.py` gains `heldout_cca_directions` (used by the direction-matched floor work;
`heldout_cca_projection` is refactored onto a shared private fit and keeps its exact association
order so no published held-out number moves by a float ulp).

## 2. The consequence, stated rather than absorbed

**No shipped artifact is rewritten.** `runs/calibra_v3_targeted/task_rows.csv` and every artifact
produced before this commit carry `observed_above_floor = 0` for all 13 states, and:

> that value is a **constant, not a measurement**. It is not evidence that any channel is below its
> floor. Any downstream reading of it — in a draft, a figure, or a ledger — is void.

Those artifacts carry no `summary_schema_version` key at all, and **absence of the key means schema
1**. That is the rule a future reader needs and it is recorded here and in
`calibration.SUMMARY_SCHEMA_VERSION`. Three keys change meaning across the boundary and are
therefore **not comparable** between a pre-fix and a post-fix run:

| key | schema 1 | schema 2 |
|---|---|---|
| `observed_above_floor` | signed random-direction statistic vs positive floor, i.e. a constant `0` | `abs(channel_statistic) > detection_floor`, `NaN` when ungradable |
| `observed_matched_direction` (per-column rows) | signed median over random directions | median absolute correlation |
| `channel_statistic` | absent | `heldout_top_cca(k, seed)` |

Re-running `run_calibra.py` on the same inputs will therefore change these columns. That is the
intended effect and it is why the version field exists.

**One call site changes verdict-bearing behaviour and it is the correct direction.**
`v2/research/rebase/nature/p4_certification/p4_certify.py` had already *amended* its criterion by
hand to `abs(heldout_single_direction_correlation) > detection_floor`, explicitly because the shipped
flag was a coin flip on a one-column `x`. That amendment is now the library's own rule, so the
certification's `shipped_flag_observed_above_floor` and its `clears_detection_floor` agree by
construction instead of disagreeing. The P4 verdicts themselves are unchanged — the amended criterion
was already what they used.

`v2/research/rebase/nature/floor_flag_diagnostic.py` is the record of the schema-1 defect and
**cannot be re-run to reproduce it**: at schema 2 the library refuses to invent a comparator, so a
call supplying no channel statistic returns `None` / `ungraded_no_channel_statistic`. The script's
docstring now says so and the JSON in `runs/floor_flag_audit/` is the record. A re-run is visibly a
different measurement rather than a silently changed one.

## 3. The sign defect in Track 1's random control — investigated, not just flagged

The antecedent entry flagged this and left it: *"Not investigated here — flagged for whoever owns
Track 1."* `v2/research/rebase/nature/track1_random_control_sign.py` measures it on the 13 shipped
states, scoring the control block and the real-target block through **both** conventions and
re-grading T1.4 with `grade_random_controls` each way. The schema-1 convention is reproduced inside
that script and **nowhere in the library**, so the defective statistic cannot be called by accident
while the comparison that condemns it stays reproducible.

The fix moves the controls **up**, i.e. it makes our own negative control *harder* to pass. Numbers
are reported in the direction-matched-floor result entry alongside the rest of this run.

## 4. New tests (5)

| test | what it pins |
|---|---|
| `test_summary_reports_observed_against_floor` (extended) | schema is 2; ungraded is `None` plus a status, never a silent `False`; graded when a channel statistic is supplied |
| `test_the_flag_is_a_magnitude_and_cannot_be_flipped_by_the_sign_of_the_channel` | B1 of the original predeclaration as a regression test: the verdict for `+r` and for `−r` must be equal at every `r` |
| `test_floors_from_recovery_reproduces_the_curve_s_own_floors` | the extracted floor rule has not drifted from the instrument — otherwise every like-for-like floor is a re-implementation |
| `test_spike_targets_accepts_a_supplied_image_direction` | the planted correlation is exact on a supplied pair; a supplied direction consumes no rng draw; a zero direction raises |
| `test_random_direction_statistic_is_a_magnitude_not_a_signed_median` | the statistic is non-negative, does **not** sit at zero on a column a random direction genuinely sees, and still reads low on noise |

## 5. Prose corrections flagged, not made

Drafts are out of scope for edits by instruction. In addition to the six passages already flagged by
the antecedent entry, this commit makes one more statement false wherever it appears:

* any text asserting that `observed_above_floor = 0` is *"the correct answer"* is now contradicted by
  the library itself, which refuses to emit that comparison at all. The flag's schema-1 value should
  be described as **a constant carried by pre-schema-2 artifacts**, and the sentence rewritten around
  `channel_statistic` vs `detection_floor`.

## 6. What was NOT done

* Shipped artifacts were **not** rewritten and `run_calibra.py` was **not** re-run over the CALIBRA
  cohort. Regenerating them is a separate, declared decision: it changes three columns of
  `task_rows.csv` and every downstream ledger that reads them.
* `v2/calibra/claim_guards.py` is untouched by instruction. It does not reference
  `observed_above_floor` (checked by grep), so no guard silently depends on the old meaning.
* `NOTEBOOK.md` and the paper drafts are untouched by instruction.
