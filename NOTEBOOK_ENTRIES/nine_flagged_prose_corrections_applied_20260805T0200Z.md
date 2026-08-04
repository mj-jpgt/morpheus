## 2026-08-05 02:00 UTC — The nine prose claims two measuring agents flagged and deliberately did not edit, applied. No measurement, no code, and four more stale sentences found and left flagged rather than guessed at

**Logged:** 2026-08-05 02:00 UTC. **How obtained:** nothing was measured. Every number written into
the draft in this entry is quoted from one of two notebook entries that measured it, and each is
attributed at the edit site. This is a prose pass and its only claim is fidelity to those two entries.

**Sources, read in full before any edit:**

* `NOTEBOOK_ENTRIES/the_dissociation_does_not_survive_its_own_floor_20260804T1800Z.md` §6, "Three
  sentences in the draft are now wrong and are flagged rather than edited".
* `NOTEBOOK_ENTRIES/three_floors_close_the_last_three_unjudgeable_rows_20260805T0000Z.md` §8, "Prose
  the other agent owns, flagged not edited".

### 1. The nine, all applied

| # | site | what was wrong | what it now says |
|---|---|---|---|
| 1 | §5.2a | *"Rank says the three high-rate arms are the same run; a co-measured collapse statistic says one of them is half as degenerate as the others"* | The dissociation does not survive a floor on the cosine: across-arm spread **0.223** against a within-arm spread of **0.250** at n = 3, so no arm difference may be read off that statistic there. The m = 0.999 arm's own four same-seed runs span **0.5207–0.9292**. The paragraph's conclusion — momentum does nothing **that rank can see** — is unchanged, because that is the only thing it ever claimed. |
| 2 | §4.10 | *"The two statistics are ordering the three runs differently at the exact reading where this section says rank is reliable"*, and the two competing accounts (A)/(B) it could not choose between | Neither statistic orders those runs at all. The competing-accounts framing is replaced by the measured fact, with the mean-offset secondary (varies 5× between identical retrains, separates no arms) as a parenthetical. **The surviving use is not under strain from this observation.** |
| 3 | §4.10 | a **new cost**, not a deletion: the mutual cosine had never had a floor | New paragraph after the "partial remedy" sentence: the cosine's first retraining floor on this project is **0.25 in absolute cosine units** (fixed held-out probe, collapse floor, `lr = 1e-3`, three same-seed retrains) — wide enough to swallow the 0.474 movement §5.2a read off it. The cosine is offered as *more legible* than rank and **not** as more reproducible, and *"and the seed spread of both"* is named as the load-bearing half of the recommendation. |
| 4 | §6.2, capacity-sweep row | *"The step was never recorded and the sweep's own logs are not vendored"*, and *"would not make the row judgeable"* | Row struck through and closed. The step **was** recorded — 150, in `qsweep_d0.04_cap64.log` and `decorr_causal_0.04.log` — and a five-column header, not missing data, kept it out of the audit. The second half of the old row is kept as what it was right about: recovering the step did not make the row judgeable; the GPU did. Clears by **1.68×**, with the cross-harness note (the sweep's 6.17 sits 0.8% above the top of its own five-run range) on the row. |
| 5 | §6.2, step-600 row | *"not measured, and it is what keeps §5.4 row 1 and §5.4 limit 2 unjudgeable"* | Struck through and closed at n = 5 per arm, fifteen runs, with the arms now part of the block string. Row 1 clears by **1.51×**, limit 2 by **1.06×**. |
| 6 | §5.4 limit 2 | *"a 1.26× difference, smaller than every floor this paper has measured on any statistic, view or block, and on a block where no floor has been measured at all"* | Both clauses withdrawn as measured-false. It clears its own two arms' floor (1.262× against **1.195×**) by **1.06×** — and the fragility is stated at equal length, not smoothed: the same ten runs separate the arms by only **1.138×** worst case under R3, *inside* that floor, while under canonical R1 they separate cleanly at **1.453×** against a 1.155× floor. Presented as §4.5's statistic-conditionality landing on the one hyperparameter value this project ships. |
| 7 | §5.4 row 1's table cell, §5.4's *"Three remain **unjudgeable**, each naming the specific run that would settle it"* | three pending | Cell reads **clears, by 1.51×** against its own 1.749× floor. The prose records that the three named runs were run and returned three passes (1.51× / 1.06× / 1.68×) and that no selection in this paper is unjudgeable any longer. |
| 8 | §4.1a, "What this costs and what it buys" | *"Nine selections do clear, eight of them §5's"* | **Twelve** clear, **eleven** of them §5's. |
| 9 | §3.1 | R2 ≡ R3 was stated only for the raw exported artifacts | Extended to the fixed held-out probe, **checked rather than assumed**: R2 = 6.9711779953 against R3 = 6.9711779832, a floating-point-level difference, because `z_biology` leaves the model L2-normalised so R3's row normalisation is a no-op. Written as a property of those blocks, **not** as a general identity — after residualisation the two separate again. |

**The "eleven of them §5's" sub-count of correction 8 was verified against the source rather than
inferred.** The instruction anticipated that the entry might not state it and asked for it to be
flagged as unverified if so. It does state it, in those words: *"It is twelve and eleven"*
(`three_floors_close_the_last_three_unjudgeable_rows` §8, item 5). It is also arithmetically
consistent with the audit's own counts — 12 clear, exactly one of them (row 30, RankMe as published)
on the exported artifact block.

### 2. Three further edits made because the nine would otherwise have orphaned them

These are not a tenth correction; each is a sentence that my own edit made false or dangling, in the
same paragraph or table as one of the nine.

1. **§5.4's table intro** — *"Two of the three quantifications … The third is still unjudgeable"* —
   became false the moment the third row's cell changed. Now "all three".
2. **§5.4's second limit on the floor** — *"The comparison clears at step 500 and cannot be judged at
   step 600"*. Rewritten as the history it now is: it *was* unjudgeable for a revision, it said what
   would close it, that was run, and the point survives the closing because it is the reason the row
   sat unjudged rather than being waved through against a floor that was never its own.
3. **§5.2 measurement 3** — *"This third measurement is the one that is still unjudgeable, and the
   reason is ours: §5.2 records no reading step"*. Both halves are now measured-false and the
   sentence is the §5-side statement of correction 4. Rewritten, with §5.2's own capacity confound
   explicitly **not** repaired by the pass.

### 3. Five stale claims found, **not** edited, and why

The instruction was nine corrections and nothing else, and to say so rather than guess where the
draft's wording did not match what the entries described. These four are contradicted by the same two
entries but are not among the nine and are not orphaned by any edit above (the fifth is the closest
to being an exception, and is called out as such). **They are left for
whoever owns them, and they are listed here so the next agent does not have to find them again:**

1. **Status block, item 8**: *"Eight further selections clear, and all eight are §5's."* It is
   eleven, by the same audit that produced correction 8.
2. **Status block, item 10**: *"Three selections remain unjudgeable, each naming the run that would
   settle it (§6.2)."* None does. Also the italic restatement in item 9 — *"Since the probe block was
   measured (item 10) they are 13 fail, 9 clear, 3 unjudgeable"* — which is now 13 / 12 / 0. The
   status block reads as a running revision log, so whether these are live claims or dated records of
   what each revision did is a call for its owner, not for a prose-correction pass.
3. **§4.1a findings 1 and 2**: *"13 fail a floor their own statistic and block license, 9 clear it and
   3 cannot be judged at all"*, *"eight of the nine passes"*, *"**Three remain**, each naming its own
   gap"*. Same counts, different sentences from the one correction 8 names. Note that the *generated*
   counting sentence at the head of the same subsection is already correct at 13 / 12 / 0, so §4.1a
   currently disagrees with itself.
4. **§5.4's "What clearing does and does not buy"**: *"This one returns nine, eight of them on a block
   we had to build the floor for, and it still fails thirteen."* Twelve and eleven.
5. **§6.2's centred-cosine row** still reads *"Not measured, and it is the one measurement that would
   settle whether §4.10's surviving use is under strain"* — it was measured, and
   `the_dissociation_does_not_survive_its_own_floor` §6 says so in as many words. Correction 2 removed
   §4.10's pointer to it, so §4.10 and §6.2 now disagree about whether that measurement exists. **This
   is the one omission that a reader is most likely to trip over**, and it was left only because it
   was outside the nine.

### 4. The generated blocks were not touched, and did not need to be

§4.1a's audit table, its measured-floor table and its counting sentence are rendered by
`v2/research/rebase/p2/p2_floor_audit.py`, and `v2/tests/test_p2_floor_audit.py` asserts the draft
prints exactly what the module renders. **They were already correct.** The measuring agent regenerated
them when it closed the three rows, so rows 43 (§5.4 row 1), 46 (§5.4 limit 2) and 50 (§5.2
measurement 3) already carried `yes` and their floors before this pass began, and the counting
sentence already read *"13 … do not clear …, 12 clear it, and none is unjudgeable"*. Nothing rendered
was hand-edited; `test_the_draft_prints_the_rendered_table` and
`test_the_draft_prints_the_rendered_floor_table` pass unchanged, which is the check that this is true
rather than an assurance that it is.

The hand-maintained tables that *do* need editing — §5.4's three-row verdict table and §6.2's
would-be-measurement table — were edited by hand, because nothing generates them.

### 5. One citation had to be reworded to pass a test, and the test was right

`test_paper_paths_resolve.py::test_every_repository_path_cited_in_a_draft_exists` failed on the first
draft of the §6.2 capacity row, on `decorr_causal.py`: the sweep's **logs** are vendored
(`v2/research/rebase/p2/figures/data/e0_run/d1_diag/`) but its **harness** is not. Reworded to "the
`decorr_causal` sweep harness (on the box; only its logs are vendored)", which is both what the test
asks for and more accurate than the source entry's bare filename. Recorded because it is exactly the
class of error that test exists to catch and it caught one on the first try.

### 6. Suite

Run as the project's convention has it — the repository as `morpheus/` in a workspace, from the
parent:

```
python -m pytest morpheus/v2/tests morpheus/tests -q
```

| | before | after |
|---|---|---|
| result | **638 passed, 1 skipped** | **638 passed, 1 skipped** |

No test changed, none was added, and no test needed to change: this pass moved prose only. The
draft-versus-generator tests passing in both states is the load-bearing part of that.

### 7. Files

* `paper/P2_RANK_DRAFT.md` — the only file edited. §3.1, §4.1a, §4.10, §5.2, §5.2a, §5.4, §6.2.
* Untouched, as instructed and as they should be: `claim_guards.py`, `claim_evidence.json`, every
  `PREDECLARED_*` file, and every block rendered by `p2_floor_audit.py`.
