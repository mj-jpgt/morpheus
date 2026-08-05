## 2026-08-05 04:30 UTC — The n = 10 reversal of §5.4 limit 2 applied to the paper: five flagged locations, all five done, and the ordering preserved everywhere

**Logged:** 2026-08-05 04:30 UTC. **Nothing was measured.** Every number written here was taken from
`NOTEBOOK_ENTRIES/limit2_breaks_at_ten_repeats_20260805T0330Z.md` and re-checked against the vendored
files that entry cites (`P2_LIMIT2_STRESS_{N5,N10,LATE5,N10_RNA}.json` under
`v2/research/rebase/p2/figures/data/e0_run/d1_probefloor600/out/`) before being written into the
draft. The source entry's §8 said the paper was being edited concurrently and flagged seven exact
locations rather than editing them; this is that list worked through.

### 0. What the correction is, in one paragraph, because it is easy to overstate in both directions

§5.4 limit 2 — `m = 0.999` over `m = 0.99` on the fixed held-out probe at step 600, the value this
project ships — **cleared** a floor measured from five same-seed repeats per arm (1.262× against
1.195×) and **does not clear** the same floor measured from ten (1.262× against **1.3263×**). The row
now fails `ratio > floor`, the paper's own rule, with no threshold moved and no statistic swapped.
**What survives untouched is the ordering**: all ten `m = 0.999` repeats sit above all ten `m = 0.99`
repeats under seven statistics, exact one-sided permutation probability **1/184,756**, and the
five-point momentum grid {0, 0.98, 0.99, 0.995, 0.999} is strictly monotone with complete separation
at every adjacent pair under canonical R1. **The trend and the ordering hold; the claim that this
particular gap clears this block's measured noise does not.** Every edit below was written to carry
both halves, because "the momentum result is dead" would be a larger claim than anything measured.

### 1. The numbers, re-checked against the vendored files rather than copied from the entry

| quantity | value | where it was read from |
|---|---:|---|
| the row's ratio (two specific published runs, fixed) | **1.2619×** | `floor_audit.json` row `5.4-m0999-over-m099`, unchanged |
| R3 floor, n = 5, carried by m = 0.999 | 1.1947× | `P2_LIMIT2_STRESS_N5.json` `arm_floor.m0999.R3.fold` |
| **R3 floor, n = 10, carried by m = 0.99** | **1.3263×** | `P2_LIMIT2_STRESS_N10.json` `arm_floor.m099.R3.fold` (max 5.9216 / min 4.4649) |
| R3 floor, repeats 6–10 alone (independent n = 5) | 1.2791× | `P2_LIMIT2_STRESS_LATE5.json` `arm_floor.m099.R3.fold` |
| R3 floor, `rna_biology`, n = 10 | 1.2764× | `P2_LIMIT2_STRESS_N10_RNA.json` `arm_floor.m099.R3.fold` |
| R1 floor, n = 10, carried by m = 0.99 | 1.2906× | `P2_LIMIT2_STRESS_N10.json` `arm_floor.m099.R1.fold` |
| R1 worst-case separation at n = 10 | 1.453× | source entry §2 (B2) and `tests.R1` |
| the shape of the m = 0.99 arm at n = 10 | `outlier: rep6`, `outlier_is_low: true`, `bimodal: false`, `rest_fold: 1.15155` | `arm_floor.m099.R3.shape` |

**Three of the prompt's restated figures were checked and are right as restated** (1.326×, 1.262×,
1.279×). The one that needed care is the R1 floor: at n = 5 it is 1.155× and at n = 10 it is 1.291×,
and the draft previously quoted "1.453× against a 1.155× floor". Quoting a ten-repeat separation
against a five-repeat floor would have been the same error in the opposite direction, so every R1
sentence now quotes **1.453× against 1.291×**.

### 2. The five locations, and what each became

1. **`v2/research/rebase/p2/floor_audit.json`.** Floor `R3_probe_step600_m0999_vs_m099`: `value`
   1.1947 → **1.3263**, `floor_arm` `m0999` → **`m099`**, `n` 5 → **10**, both `src` blocks and the
   `shape_src` repointed from `P2_PROBE_FLOORS_S600_m0999_m099.json` to `P2_LIMIT2_STRESS_N10.json`,
   and the caveat rewritten to record the re-measurement, the one repeat that carries it and the four
   exclusion routes that were checked and did not apply. **Its R1 twin `R1_probe_step600_m0999_vs_m099`
   was moved with it** (1.1547 → **1.2906**, same arm change, same file) — the two are the same twenty
   runs and leaving one at n = 5 beside the other at n = 10 is exactly the drift this audit exists to
   catch. Row `5.4-m0999-over-m099`: `clears` `true` → **`false`**, `ratio` unchanged at 1.2619, and
   `rests_on` replaced with the §0–§1 and §5 reading of the source entry — the failure, the two
   independent confirmations, the non-excludable repeat, **and the ordering that survives**.
   **Then regenerated**, not hand-edited: `p2_floor_audit.py --check` reports no disagreements, and
   `--markdown` / `--floors` / `--sentence` were written into the draft by a script, so the two
   draft-versus-generator tests still pass.
2. **§4.1a's counts.** Verified against the regenerated output rather than computed by hand:
   `summary()` returns `selection_failing: 14`, `selection_clearing: 11`, `selection_unjudgeable: 0`.
   The generated counting sentence moved from **13 fail / 12 clear** to **14 fail / 11 clear**, and
   every prose restatement of it moved with it — §4.1a finding 1, §4.1a's "What this costs and what it
   buys", §4.1b's "**Eleven** further selections clear" (now ten), §5.4's "13 of the 25 selections",
   and Appendix C's audit bullet. **A new finding 2a** in §4.1a states the reversal in full. The
   counting history gains its first step in which a row leaves the *clearing* column — every previous
   step moved rows out of *unjudgeable* — and the line "the thirteen failures have never moved" was
   rewritten to "no failure has ever moved into the clearing column; one has now moved the other way",
   which is the same fact stated so that it stays true.
3. **§5.4 limit 2's prose.** The "clears it by 5.6%" / "the rule now licenses it" framing is retracted
   in place with the retraction stated rather than the text quietly swapped. Limit 2 now runs: the
   failure and the predeclared falsifier that called it; **(a)** both independent halves and the second
   view agreeing; **(b)** the one repeat that moves it and the four exclusion routes checked against
   it; **(c)** the ordering, exact across 10 × 10 under seven statistics, with the explicit note that
   1/184,756 is a statement about GPU non-determinism at a fixed seed and **not** a p-value for the
   momentum effect; **(d)** the five-point grid, monotone and completely separated at every adjacent
   pair under R1, with the fact that only one of the four adjacent rungs clears its own floor — so
   this row's gap is the sum of two rungs each inside it. The limit specific to this row is now
   stated: **it can be judged only under R3**, because the two runs it quotes predate the canonical
   column and the state export, so the one statistic that can rule on it is the one it fails.
   §6.2's closed-row entry and Appendix C's momentum bullet were rewritten the same way.
4. **The Status block.** Item 12 was **not** edited — it records what was true when it was written.
   **Item 14 was appended** (13 was the highest), recording the reversal, the predeclaration that
   called it, the two independent confirmations, and the ordering and grid that survive, in the voice
   of the existing items and citing the source entry.
5. **The 00:00 UTC entry.** `three_floors_close_the_last_three_unjudgeable_rows_20260805T0000Z.md`
   gained a **`CORRECTION APPENDED`** block in the style of `t13_adjusted_certificate_and_p6`;
   nothing above it was altered. It withdraws seven located passages by line, all of them about
   §5.4 limit 2 only — the other two rows that entry closed are untouched — and it states the thing
   that entry got most wrong: its §2 offered this row as the example showing that *"measure both
   sides"* returns the stable arm. **The rule is vindicated and the example inverts**: at ten repeats
   the m = 0.99 arm is the noisier one and carries the floor, and it does so on the strength of one
   repeat in ten, which is the sharpest argument in the project against reading a five-repeat floor as
   a property of an arm.

### 3. Tests

`v2/tests/test_p2_floor_audit.py` pinned the old verdict in five places and each was updated to pin
the new one rather than relaxed:

* `test_every_floor_names_a_statistic_a_block_and_a_caveat` — `n` may now be 10, **only** for the
  `_probe_step600_m0999_vs_m099` pair and only with the re-measurement recorded in its caveat.
* `test_the_selections_that_clear_are_exactly_these` — `5.4-m0999-over-m099` removed from the list.
* `test_the_last_three_unjudgeable_selections_were_closed_by_their_own_floors` — now asserts per row
  what its verdict and repeat count are; the point it still enforces is that limit 2 is **judged**,
  and that the verdict against it is a failure rather than a silence.
* `test_the_value_this_project_runs_passes_narrowly_and_the_row_says_so` → renamed
  `..._does_not_clear_at_ten_repeats_and_the_row_says_so`, and it asserts **both halves**: the three
  numbers that make the failure not a one-draw artifact (1.3263, 1.2791, 1.2764) **and** the three
  that make the ordering survive (1/184,756, "monotone", 1.453). A correction that dropped the second
  half would overshoot, so the test refuses it.
* `test_the_probe_floor_is_carried_by_the_collapsed_arm_and_is_not_bimodal` — the arm carrying this
  floor is now pinned as `m099` at `n = 10`, which is what would catch it inverting again.

Nothing in `claim_guards.py`, `claim_evidence.json`, `p2_limit2_stress.py`,
`test_p2_limit2_stress.py`, or any `PREDECLARED_*` file was touched.

### 4. Suite state, verbatim

Run as `python -m pytest v2/tests -q -p no:randomly --basetemp=<scratch>` with the repository
reachable on `PYTHONPATH` under the package name `morpheus`. The **before** run is on a detached
worktree at `38e81c2`, so it is the same suite on the same machine rather than a figure quoted from
an earlier entry:

* **Before (worktree at `38e81c2`): `610 passed, 1 skipped, 445 warnings in 131.77s`.**
* **After all edits: `610 passed, 1 skipped, 445 warnings in 136.59s`.**

Zero failures either side and no change in the number of tests: the five assertions above were
rewritten, not added or removed. `test_the_draft_prints_the_rendered_table` and
`test_the_draft_prints_the_rendered_floor_table` both pass, which is the check that §4.1a's two
tables came from the generator and not from a hand edit.

### 5. What was deliberately not done

* **§5.2's `m = 0.9` dip is still open.** The source entry §9 records that §5.2's table reads
  `m = 0.9` at 2.23 against `m = 0`'s 2.81, one seed each, against prose saying the effect is monotone
  in `m`, and that the five-point grid does not cover that interval. No repeat was run at `m = 0.9`.
  It is not one of the flagged locations and it is not touched here; closing it costs five runs.
* **Status items 6, 10, 12 and 13 were left exactly as written.** They are the record of what was
  believed when each was written, and this project appends rather than rewrites.
* **The abstract was not edited.** It carries no audit count, so the reversal does not reach it; its
  §5.4 sentence is about the momentum fix (row 2), which still clears by 2.4×.
