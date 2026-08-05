## 2026-08-05 09:00 UTC — Composing two same-day corrections to §4.1a: the honest total is **15 fail / 10 clear / 0 unjudgeable**, not either correction's own arithmetic

**Logged:** 2026-08-05 09:00 UTC. **Outcome: RESULT (a composition, not a new measurement).**
No new statistic was computed in this entry. Every number below is read back out of
`v2/research/rebase/p2/floor_audit.json` as it exists on disk, before and after the edits this entry
describes, via `python v2/research/rebase/p2/p2_floor_audit.py`.

---

## 0. The awkward finding, stated first

**Two corrections to `paper/P2_RANK_DRAFT.md` §4.1a landed on 2026-08-05, each internally correct, and
their own stated totals do not compose — not because either is wrong, but because the second one's
recomputation was performed against a workspace frozen *before* the first one landed.**

1. `4961a0d` ("The paper records that §5.4 limit 2 fails...") applied the ten-repeat momentum-floor
   result to `floor_audit.json` and the draft, moving the audit from **13 fail / 12 clear** to
   **14 fail / 11 clear** of 25 selections. Committed 2026-08-05 02:04 EDT (06:04 UTC).
2. `NOTEBOOK_ENTRIES/unstable_arm_exported_floor_measured_20260805T0755Z.md` (commit `a181ed1`,
   2026-08-05 03:57 EDT / 07:57 UTC) measured the exported floor's `programme_free` arm, confirmed
   3.295× stands, and found that the same measurement reverses the paper's RankMe self-criticism —
   row 30 (`4.6-rankme-d2`) flips clear → fail — and shrinks the `rna_biology`/`full_biology` view
   counts from "12 of 12" to "7 of 12" (§4, §5 of that entry). It states the resulting count as
   **"Selections: 12 clear / 13 fail → 11 clear / 14 fail of 25, 0 unjudgeable"**.

**Those two headline totals — item 1's `14 fail / 11 clear` and item 2's stated `11 clear / 14 fail` —
are numerically identical, and that identity is a coincidence, not a confirmation.** Item 2's workspace
was verified byte-equal to commit `a392c0a` — the predeclaration commit, timestamped 2026-08-05 00:45
UTC, which is *before* item 1 (06:04 UTC) landed. Item 2's own recomputation in its §4 ("Recomputed
read-only from `floor_audit.json`'s recorded ratios... nothing was written to `floor_audit.json`")
therefore started from the **pre-item-1** baseline of 12 clear / 13 fail, not the actual on-disk state
at the time item 2 was written. Applying row 30's flip to *that* baseline gives 11 clear / 14 fail —
which happens to equal item 1's own already-published total, purely because both corrections move
exactly one selection from clear to fail and the two baselines differ by exactly one selection in the
same direction. **Composing both corrections against the true current state of `floor_audit.json`
gives a different, larger total: 15 fail / 10 clear.** Nobody wrote 15/10 down before this entry.

This is reported first, per this project's standing rule (`PROJECT_GUIDE.md` §2 rule 2), because it is
exactly the failure mode rule 16 exists to catch — two correct, independently-measured deltas that do
not compose the way either one's own restatement implies, discovered only by re-deriving from the file
on disk rather than trusting either summary's arithmetic.

## 1. What was actually on disk, verified before touching anything

```
$ python v2/research/rebase/p2/p2_floor_audit.py
{
 "total": 62, "selection": 25,
 "selection_failing": 14, "selection_clearing": 11, "selection_unjudgeable": 0,
 "exempt": 5, "no_floor_measured": 16, "block_mismatched": 38
}
EXIT: 0   # zero DISAGREEMENT lines
```

`floor_audit.json` on disk, before this entry's edits, already carried item 1's momentum-floor
correction (`5.4-m0999-over-m099` reads `clears: false`, floor `R3_probe_step600_m0999_vs_m099` at
`n=10`, value `1.3263`) and did **not** yet carry item 2's floor changes — confirmed by inspecting the
seven floors item 2 names (`RankMe_published_raw_export` = 1.811, `LiDAR_residualised_export` = 1.06,
`LiDAR_raw_export` = 1.034, `R1_residualised_rna_view` = 1.019, `R1_raw_rna_view` = 1.023,
`R1_residualised_full_view` = 1.02, `R1_raw_full_view` = 1.014 — all still the `programme_only`-only
values item 2's own §7 says it deliberately did not write). This confirms item 2's entry text
("`floor_audit.json` deliberately not edited here") and independently confirms the baseline
discrepancy described in §0: the file's actual pre-this-entry state is 14 fail / 11 clear, not the
12/13 item 2's own arithmetic assumed.

## 2. Applying item 2's Class E changes to the actual current file

Seven floor values updated, sourced from
`v2/research/rebase/p2/figures/data/e0_run/d1_envelope_pf/out/P2_ENVELOPE_FLOORS_PF.json` (the
`programme_free` arm, vendored in commit `a181ed1`) and cross-checked against
`P2_ENVELOPE_FLOORS_PO_RECHECK.json` (the `programme_only` re-check, confirming it reproduces the
already-published values exactly):

| floor | old (`programme_only` alone) | new (both arms, `max`) | carried by |
|---|---:|---:|---|
| `RankMe_published_raw_export` | 1.811 | **3.5484** | `programme_free` |
| `LiDAR_residualised_export` | 1.060 | **2.5588** | `programme_free` |
| `LiDAR_raw_export` | 1.034 | **2.1708** | `programme_free` |
| `R1_residualised_rna_view` | 1.019 | **1.2305** | `programme_free` |
| `R1_raw_rna_view` | 1.023 | **1.2376** | `programme_free` |
| `R1_residualised_full_view` | 1.020 | **1.4220** | `programme_free` |
| `R1_raw_full_view` | 1.014 | **1.4355** | `programme_free` |

Each floor's `a`/`b` sources were repointed at the `P2_ENVELOPE_FLOORS_PF.json` max/min for the
relevant statistic/block/view, its `shape` block updated to that file's shape descriptor (all seven are
non-bimodal, unlike the `wsi_biology` residualised/raw floors), and a `floor_arm: "programme_free"`
field added — matching the existing convention used by the `*_probe_step*` floors, which already carry
`floor_arm` for the same reason (a floor measured as the max of two arms names which arm carries it).
`R1_residualised_export` and `R1_raw_export` (the `wsi_biology` floor itself, 3.295×/3.111×) were **not**
touched: `programme_free`'s reading there (1.4254×/1.4395×) is smaller than `programme_only`'s, so
`programme_only` continues to carry it and the value is unchanged — this is the "3.295× stands" half of
item 2's finding, and it required no edit.

Row `4.6-rankme-d2`'s `clears` flipped `true → false` and its `rests_on` rewritten to state the
reversal. Rows `4.5c-rna` and `4.5c-full` (kind `direction`, not `selection` — they do not enter the
25-selection count, but their prose was wrong under the old floors) had their `rests_on` rewritten with
the new floors and the per-pair breakdown below.

## 3. The composed result, verified against the checker and `summary()`

```
$ python v2/research/rebase/p2/p2_floor_audit.py
{
 "total": 62, "selection": 25,
 "selection_failing": 15, "selection_clearing": 10, "selection_unjudgeable": 0,
 "exempt": 5, "no_floor_measured": 16, "block_mismatched": 38
}
EXIT: 0   # zero DISAGREEMENT lines — every source re-resolves, every block-match holds
```

**15 fail / 10 clear / 0 unjudgeable of 25 selections.** Not 13/12 (the prompt's naive guess, which
double-composes neither correction). Not 11/14 as item 2's own text states (which silently drops item
1's momentum correction, per §0 above). `check()` reports zero disagreements: every recorded value
re-resolves from its named source, every ratio agrees with its two values, and every block-match holds.

Zero selections now clear a floor on the exported artifact block. Ten selections clear, all ten on the
fixed held-out probe (§5's block) — the same ten that were already clearing before this entry; row 30
was the only one clearing on the exported block, and it is the row that flipped.

## 4. The view counts and the D2 arm-floor gap, verified per pair

Per-pair recomputation, direct from `v2/research/rebase/p2/figures/data/ws_p2/out/P2_ROBUSTNESS.json`
(rank) and `P2_METRICS_D2.json` (RankMe), against the new floors:

**`rna_biology`, floor 1.2305×:**

| pair | fold | clears? |
|---|---:|:---:|
| D2 s42 (H42/I42) | 1.2054 | no |
| D2 s43 (H43/I43) | 1.1160 | no |
| D2 s44 (H44/I44) | 1.2380 | **yes (by 0.61%)** |
| D1 s42 (P42/F42) | 1.7659 | yes |
| D1 s43 (P43/F43) | 2.8522 | yes |
| D1 s44 (P44/F44) | 3.0137 | yes |

**4 of 6.**

**`full_biology`, floor 1.4220×:**

| pair | fold | clears? |
|---|---:|:---:|
| D2 s42 | 1.2345 | no |
| D2 s43 | 1.0424 | no |
| D2 s44 | 1.1404 | no |
| D1 s42 | 2.2484 | yes |
| D1 s43 | 3.6057 | yes |
| D1 s44 | 5.2498 | yes |

**3 of 6.** Total **7 of 12**, matching item 2's own §5 table exactly (confirming that half of item
2's arithmetic, unlike its selection-count summary, was already correct against the pre-item-1
baseline and remains correct against the composed one, because the view-count computation in item 2
did not depend on the momentum row at all).

**RankMe D2, floor 3.5484×:**

| pair | fold | clears? |
|---|---:|:---:|
| H42/I42 | 1.6771 | no |
| H43/I43 | 3.3817 | no |
| H44/I44 | 1.2481 | no |

**0 of 3.** All three fail; s43 (3.3817×) is the row that used to clear the old 1.811× floor.

**Two open items this composition makes explicit, both already true before this entry but sharper
now:**

- **No D2 arm (Hallmark or PBS) has ever had its own retraining floor measured, on either training
  arm.** Every both-arm correction in this paper — the view floors, the RankMe floor — reaches only
  the D1-licensed pairs (`programme_only`/`programme_free`) directly. The five losses from 12/12 to
  7/12 are a floor *transferred* from D1's arms, not a measurement made on D2's own arms. Closing this
  is five same-seed retrains per D2 arm (Hallmark and PBS), ten runs, GPU — not currently scheduled.
- **The one D2 survivor, `rna_biology` seed 44, clears by 0.61%** (1.2380× against 1.2305×) — a margin
  this paper does not otherwise treat as safe to lean on (compare §5.4 limit 2's 5.6%, which broke at
  n=10). It is flagged rather than quoted as a clean pass.

## 5. What changed in the draft and the audit, and what a test enforces now

`floor_audit.json`: seven floor values, their sources, shapes and caveats; row `4.6-rankme-d2`'s
verdict; rows `4.5c-rna`/`4.5c-full`'s prose. `paper/P2_RANK_DRAFT.md` §4.1a's table, floor table and
counting sentence were **regenerated** via `p2_floor_audit.py --markdown --floors --sentence` and
pasted in verbatim (never hand-edited) — `v2/tests/test_p2_floor_audit.py::test_the_draft_prints_*`
assert this. Every other prose location the abstract, §1.3, §1.4, §4.1a's "five things" and "what the
audit found", §4.1b, §4.3, §4.5(c), §6.1, §6.2, the Conclusion, Appendix A and Appendix C carries was
walked and updated to the composed state — a `grep` sweep for `1.811`, `12 of 12`, `only selection`,
`1.019`, `1.020×`, `3.2× between` and `fifty times` across the draft, resolved location by location,
leaving only the Status block's historical items (1–14) untouched per this project's append-only
convention (`PROJECT_GUIDE.md` §2 rule 9) and a new item 15 appended describing this composition.

Two tests were updated **as a decision, not a repair**, per item 2's own §7 prediction that this would
be necessary:

- `test_the_selections_that_clear_are_exactly_these` — `4.6-rankme-d2` removed from the pinned
  clearing-selection list (ten remain, all on the probe block).
- `test_the_floor_is_a_property_of_the_statistic_and_of_the_view` — the `rna_biology`/`full_biology`
  floor assertions widened from `< 1.05` (the single-arm bound) to `> 1.05` and `< 1.30`/`< 1.50` (the
  both-arm bound, still well below `wsi_biology`'s), and the RankMe-vs-R1 raw-floor inequality
  **reversed** (RankMe's floor is now asserted **greater than** R1's, not less).

No other test's expectations needed to change; `check()` reports zero disagreements against the
composed file.

## 6. Suite

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python -m pytest v2/tests tests -q --basetemp=./pytmp_composed
```

Result: **704 passed, 1 skipped, 447 warnings** (one failure fixed mid-session: a path-resolution test
caught two citations in the draft that did not yet exist — this entry's own path, before it was
written, and a brace-expansion path glob that is not a literal file. Both fixed; see commit.)

## 7. Files

- `v2/research/rebase/p2/floor_audit.json` — seven floors, one selection verdict, two `direction` rows'
  prose.
- `paper/P2_RANK_DRAFT.md` — §4.1a's three generated blocks; abstract, §1.3, §1.4, §4.1a prose,
  §4.1b, §4.3, §4.5(c), §6.1, §6.2, Conclusion, Appendix A, Appendix C; Status item 15 appended.
- `v2/tests/test_p2_floor_audit.py` — two assertions updated as a decision.
- **Not touched:** `v2/calibra/claim_guards.py`, `v2/research/rebase/nature/claim_evidence.json`, any
  other agent's `PREDECLARED_*` file, `figures/data/e0_run/d1_envelope_pf/out/` (read only).
