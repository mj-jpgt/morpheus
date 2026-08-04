## 2026-08-04 19:00 UTC — §5 rewritten around the probe-block floor: our own fix clears, four prose errors are corrected at source, and §4.1's bimodality is scoped to the block it was measured on

**Logged:** 2026-08-04 19:00 UTC. **How obtained:** prose only. CPU only, thread-capped; no run
launched, no statistic recomputed. Every number below is read from
`v2/research/rebase/p2/figures/data/e0_run/d1_probefloor/out/P2_PROBE_FLOORS.json`,
`v2/research/rebase/p2/floor_audit.json`, or the eighteen logs vendored under
`v2/research/rebase/p2/figures/data/e0_run/d1_diag/`. Written against
`the_probe_block_has_a_floor_at_last_20260804T1620Z.md` and
`lr_logs_vendored_and_the_step_budget_they_contradict_20260804T1400Z.md`, both of which flagged rather
than changed the prose because §5 is this agent's.

### 1. The structural change: §5.4's headline reverses, and the reason is the whole point

**§5.4 row 2 — the momentum fix, the paper's own intervention — now clears.** It has held three
different verdicts in three revisions, and §5.4 now prints all three because the sequence is the
argument:

| revision | judged against | verdict |
|---|---|---|
| before block-matching | 3.295× — canonical R1, **exported residualised `wsi_biology`**, `programme_only`, 40 epochs | **fails by 0.3%** |
| after block-matching, before the probe floor | nothing — the probe block had no measured floor | **unjudgeable** |
| now | **1.367×** — canonical R1, **fixed held-out probe, step 500**, ten same-seed repeats of the two arms compared | **clears by 2.4×** |

**No threshold in this paper was changed at any point.** 3.295× was never this comparison's floor —
different block, different arm, different duration — and §4.1a's block-matching rule said so in
writing *before* the right floor existed, which is why the middle state is `unjudgeable` rather than a
pass. §5.4 states this explicitly, at length, because the distinction between measuring the right
floor and moving a threshold is the difference between an audit and a rationalisation, and a reader
who sees a paper's self-criticism reverse will want it.

The same reversal applies to §4.4(3)'s fixed-seed probe repeat (3.495× against R3 probe step-200
2.041×) and to §4.9a's decorrelation ablation (1.854× / 1.940× against the step-400 floors of 1.449× /
1.570×). **§4.9a's clearing does not change what §4.9a rests on** — one seed per level means the
magnitude carries nothing and the monotonicity plus the co-measured cosine carry everything — and the
draft says so rather than upgrading the claim.

### 2. The bimodality is scoped to the exported block, throughout

**Twenty probe floors — two statistics × five reading steps × two arms — are none of them bimodal**,
with remaining-four folds of 1.17×–1.64× against the 1.05× a four-run agreement requires. Divergence
on the probe is **graded across all five runs**; divergence on the exported artifact is **one run
falling off a cliff**. These are not the same phenomenon, and "the retraining envelope is bimodal" was
being stated unscoped in the draft, in this project's summaries, and in `P2_FIGURES.md`.

Scoped in eleven places: the status block, both abstracts, §1.3's practitioner rule, §1.4's
contribution 1, §4.1's heading and its "80% / 20%" sentence and its "several repeats" rule, §6.2, the
conclusion, Appendix C, and F1's caption and the figure plan's pending-dependency row. The
eighty/twenty sentence in particular may not travel without "on the exported artifact block".

### 3. The probe floor is carried by the collapsed arm, and that is a limitation of §4.1's floor too

| | across the five reading steps |
|---|---|
| **m = 0 (collapsed) alone** | **1.333×–2.057×** |
| m = 0.999 (stable) alone | **1.089×–1.165×** |

**Measuring one arm — the healthy one, the natural choice — would have published about 1.1× and
flattered every row in the audit by a factor of two.** That is exactly what §4.1's exported-block
floor does: five repeats of `programme_only`, this project's *stable* arm. §4.1 already called itself
a floor twice over for that reason; this is the first **quantification** of what that costs, and it is
now stated beside both floors rather than only beside the new one.

### 4. Four prose errors, corrected at source, all in §5

Each was reported in `floor_audit.json`'s `known_source_disagreements` before being corrected, and
none was ever substituted silently.

1. **§5.2a said 400 steps; the logs say `steps=200`.** Corrected in §5.2a's table intro and provenance
   line, the status block, Appendix A, Appendix C, `QUEUE_ANCHORING.md` and `P2_FIGURES.md` S9. Every
   rank value was unaffected — 200 is the last row each log has, which is what makes the quoted values
   the "final eff-rank" — and the predeclaration had fixed step 200 in advance, so the logs agreed
   with the predeclaration and the prose did not. A correction note is appended to
   `learning_rate_is_the_mechanism_20260805T0100Z.md`.
2. **§5.2's table and its "measurement 2" are two run families at one configuration, presented as
   one.** They are now named **family A** (`long_m*.log`, 1,500-step budget, read to 600) and
   **family B** (`mom_*_d0.04.log`, 300-step budget), and every subsequent sentence says which. See §5
   below for what the split buys.
3. **§5.2's "40 epochs = 583 steps"** is the epoch equivalence of the training this objective is used
   for, not the harness budget, which the logs give as **1,500**. Both are now stated, in both the
   draft and `QUEUE_ANCHORING.md`.
4. **§6.4 said "six rows"** where §5.2a has **four** (audit rows 57–60). Corrected, and the whole
   bullet is rewritten from an open provenance gap into a closed one.

### 5. The silver lining, used rather than buried

Family A against family B is an **n = 2 same-configuration retraining spread on the fixed held-out
probe at step 100** — the exact block and step the new floor covers — from two launches that were
never intended as repeats:

| m | family A | family B | fold |
|---:|---:|---:|---:|
| **0** | 1.62 | 2.58 | **1.593×** |
| 0.99 | 6.49 | 6.65 | 1.025× |
| 0.999 | 7.03 | 6.89 | 1.020× |

The five-repeat floor measured at that step, on that arm, under the same statistic, is **1.494×**.
**Two entirely separate sets of runs agree that this block's step-100 floor is about
one-and-a-half-fold, and neither was designed to measure it.** Stated honestly: the n = 2 figure sits
*above* the n = 5 range rather than inside it, and the residue is the one thing the two families do not
share — a 1,500-step budget against a 300-step one. **An n = 2 spread is not a floor** (§4.1's own
argument), and both rows are in the audit as `nuisance` — measurements *of* the noise, not selections
against it (rows 61–62).

### 6. The three that remain unjudgeable, and whether any can be closed by reading a log

**None of the three can.** Each is recorded in §6.2 as its own row rather than folded into one.

* **§5.4 row 1** (2.641×, step 600) — the *reading* exists in the family-A logs; the *floor* stops at
  step 500, and a step-500 floor may not be applied to a step-600 reading any more than a raw floor
  may be applied to a residualised ratio. Needs five repeats of each arm at a 600-step budget. GPU.
* **§5.4 limit 2** (1.262×, step 600, m = 0.99) — the same, plus a third arm at m = 0.99, which is
  neither arm the repeats cover. GPU.
* **§5.2 measurement 3** (2.857×, capacity 64) — **§5.2 records no reading step for the capacity
  sweep at all**, and the sweep's own logs are not vendored. Vendoring them would recover the step —
  a file copy — but **would not make the row judgeable**, because every probe repeat is at capacity
  4,096 and a floor measured at one capacity may not be borrowed for another. This is the one of the
  three whose gap is partly ours: we did not record the step.

§5.2a's four `direction` rows also stay unjudgeable, and the reason is now stated rather than
inherited: the probe floor exists but was measured at **`lr = 2e-4`**, and those arms are at `1e-3`
and `4e-5` — and §5.2a's own result is that the learning rate is the variable that moves rank most, so
a floor measured at one rate is the last quantity that may be borrowed for another.

### 7. What clearing does not buy, said in the draft before a reader has to ask

**13 of the 25 selections still fail** a floor their own statistic and block license, including all
seven of §4.1's between-arm differences, and the one selection that clears on the *exported* block is
still **RankMe as published** rather than ours. The paper's criterion now returns all three of its
possible answers — 13 fail, 9 clear, 3 unjudgeable — which is the shape a criterion should have. A
criterion that never returned a pass would have been a rhetorical device.

### 8. Suite

**427 passing**, unchanged. `p2_floor_audit.py --check` reports no disagreements at 62 rows.

One mechanical note worth recording because it cost a debugging cycle: `p2_floor_audit._section`
terminates a section at the next heading of **any** level, so a new `####` heading inside §5.2 silently
truncated the section and broke row 50's markdown source resolution. The family-split block is a bold
lead-in paragraph rather than a heading for that reason.

### Files / commits

- `v2/research/rebase/p2/figures/data/e0_run/d1_probefloor/out/P2_PROBE_FLOORS.json`;
  `v2/research/rebase/p2/floor_audit.json` (62 rows).
- `paper/P2_RANK_DRAFT.md` §4.1, §4.1a, §4.9a, §5.2, §5.2a, §5.4, §6.2, §6.4, §7, Appendices A and C,
  and the status block (new items 10 and 11).
- `paper/P2_FIGURES.md` — new **T10**; F1, F9, S4, S9 and two pending-dependency rows updated.
- `paper/QUEUE_ANCHORING.md` — header verdict, the family-A budget, the 200-step correction.
- Sources: `the_probe_block_has_a_floor_at_last_20260804T1620Z.md`,
  `lr_logs_vendored_and_the_step_budget_they_contradict_20260804T1400Z.md`.
