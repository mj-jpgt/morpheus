## 2026-08-04 20:00 UTC — P2 §5 rewritten around the momentum seed replication: the fix's own rank difference is inside our floor, and the fix never rested on it

**Logged:** 2026-08-04 20:00 UTC. **How obtained:** writing pass on `paper/P2_RANK_DRAFT.md`,
`paper/QUEUE_ANCHORING.md` and `paper/P2_FIGURES.md`, CPU only, no job touched, no GPU used, no number
recomputed. Every value below is quoted from an existing source with its statistic and block; nothing
was computed inline and no rank statistic was evaluated in this pass.

### 1. What landed, and what it costs

The armed `m ∈ {0, 0.999} × seeds {42, 43, 44}` sweep reported
(`retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §3). Canonical Roy & Vetterli order 1,
fixed held-out probe, 500 steps:

| m | s42 | s43 | s44 | within-arm spread |
|---|---:|---:|---:|---:|
| 0.999 | 11.26 | 10.45 | 10.55 | 1.08× |
| 0 | 3.18 | 1.13 | 2.36 | 2.81× |

Every m = 0.999 seed exceeds every m = 0 seed, on canonical R1 and on R3. The single-seed defect that
made §5 strain against §4.1 is **closed** — and the single seed was a **hardcoded harness parameter**,
i.e. a defect rather than a design choice, which the draft now says in §5.2 rather than leaving it as a
"Limits" bullet that reads like a judgement.

**And the number that must not be buried.** Worst-case separation is `min(m = 0.999) / max(m = 0)` =
**10.45 / 3.18 = 3.29×**, against §4.1's measured floor of **3.295×**. By the paper's own criterion the
momentum fix's *rank* difference is **not resolvable either**. On R3 it is 6.85 / 2.81 = 2.44×, also
inside. The original single-seed step-600 comparison (7.42 / 2.81 = 2.64×) was inside it too.

So **every rank comparison this project has ever made now fails this paper's own test** — the seven
between-arm differences of §4.1, and the hyperparameter choice of §5.2.

### 2. How §5 now reads

- **§5.1** unchanged (the liveness gate).
- **§5.2** retitled *"…found in rank, established on a binary outcome"*. It gains the seed-replication
  table with its statistic and block stated (canonical R1 and R3 on the **fixed held-out probe**, not
  §4.1's residualised exported block), and the hardcoded-seed defect moved up from §5.3 into it.
- **§5.3** is now only *"Why §5 is a demonstration and not a contradiction"*. The two structural
  defences that survive — monotone across four m values, readings at the collapse floor — stay there,
  with an explicit line that neither makes the ratio resolvable.
- **§5.4 is new and is a subsection, not a footnote**: *"The paper's own standard, applied to the
  paper: our own fix's rank difference is inside our own floor."* Bad news first, in a table of all
  three quantifications against the floor; then the three qualifications (arm, duration, block) stated
  and then refused, because all three point toward a **larger** floor in this regime rather than a
  smaller one; then the resolution.

**The resolution, which is the section's argument.** The momentum fix was never justified by a rank
ratio. Rank is how the problem was *noticed* — 67.55 → ~2 within 150 steps under every regulariser
setting, invisible to a gate that freezes the queue — which is the collapse-diagnostic use §4.10
defends. What established the repair is a **binary outcome**:

| | before the fix (D1-A) | after the fix (D1-B, `--biology-key-momentum 0.999`) |
|---|---|---|
| `programme_free` reaching epoch 40 | 1 of 3 seeds; the other two refused by the gate at 0.50883 and 2.14122 | **3 of 3** |
| state of the arm that did | **collapsed** — R3 1.71 at epoch 39, RNA-view mutual cosine 0.986, hard rank 11 | **not collapsed** — canonical R1 13.418 / 7.600 / 6.394 |
| exports, CALIBRA, bootstrap | **none** (`run_d1` raises on the first non-zero return code) | all six runs complete; **three paired bootstraps with CI₉₅ per seed** |
| held-out channel | not measurable | 0.5412 / 0.5336 / 0.5126, above its own `random_control` |

**Uncollapsed completions: 0 of 3 before, 3 of 3 after.** That is a change of kind, not of degree, and
it is not a ratio inside anyone's floor.

Three limits are stated with it: it is a before/after across two launches rather than a matched
ablation; what it supports is **momentum versus none**, not `m = 0.999` over `m = 0.99` (1.26× at step
600, deep inside every floor here, so **the specific value we run is selected by a comparison our own
rule disqualifies**); and the replication is 500 steps, so the 40-epoch evidence is D1-B's completion
rather than the rank curve.

### 3. The predeclared disjunction was under-specified, and the draft says so

§5.3 had committed, before the result: *separation ⇒ §5.2's fix clears §4.1's bar and the tension
disappears; overlap ⇒ the momentum choice is a rank comparison this paper's rule disqualifies and §5.2
must be rewritten to rest on downstream behaviour.* **The distributions separated and the ratio still
did not clear the bar.** The disjunction conflated *do the arms separate* with *does the separation
exceed the floor*, and was written before the floor's value existed. §5.4 records that it was
under-specified rather than pretending it resolved cleanly, and takes the branch it assigned to
overlap.

### 4. Also done

- **Block-matching made load-bearing in §5.4.** D1's ratios are quoted in both forms — raw
  2.02 / 3.09 / 1.68 and residualised 2.190 / 3.246 / 1.738 — with the observation that seed 43's
  residualised 3.246× judged against the **raw** floor of 3.111× would read as *outside* when on its own
  block it is inside. The momentum numbers sit on a third block with no floor of its own.
- **§4.4(3)'s fixed-seed repeat reported honestly rather than dropped.** It is the one rank measurement
  on this fix that clears the floor (3.5×, empty band 1.98–6.92) — by about 6%, and it holds the seed
  fixed, which is the term §4.2 measures as dominant. The seed-varied version is the check §1.3 asks
  for and it fails.
- **Floor quoted at three decimals (3.295×) everywhere it appears in this pass.** The stale 3.30× in the
  status block and the superseded n = 1 estimate of 2.69× in `QUEUE_ANCHORING.md`'s header are gone;
  the latter is now explicitly flagged there as not-the-floor.
- **Framing preserved.** §5's preamble restates that this is a **worked example of the metric used in
  the regime the paper says it works, not a separate contribution**, and that it is not to be split out.
- Updated: draft status block item (c), the long abstract, §3.6 rule 6, §4.7.3's deflation bullet, §6.2's
  momentum row, §7, Appendix A (two new provenance rows) and Appendix C's momentum bullet;
  `paper/QUEUE_ANCHORING.md`'s superseding header; `paper/P2_FIGURES.md` S4 status and caption, the
  "figures the paper does not have" table and the pending-dependencies table.

### 5. One source disagreement, reported rather than quoted

**The assignment's framing that `programme_free` "had never trained to 40 epochs in this project's
history" does not match the source, and the draft does not say it.**
`d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md` records `d1_f_seed42` as
"40 epochs ✓" in D1-A — pre-fix — with that arm sitting at R3 **1.71** and RNA-view mutual cosine
**0.986**, i.e. it *ran* to 40 epochs and was collapsed, and the same entry records that `run_d1`
produced no exports, CALIBRA or bootstrap for that launch at all.

The binary claim survives in the accurate form and is written that way throughout: **0 of 3 seeds
completed 40 epochs uncollapsed before the fix, 3 of 3 after**, with the one pre-fix completion
identified as collapsed and with the absence of any exported readout stated. Nothing was quoted in the
stronger form.

*(A second, unresolved oddity noticed while tracing this and deliberately not relied on: the notebook
timeline for `~/e0_run/d1_v2` is not internally consistent — `p2_rank_draft_20260803T2134Z.md` records
D1-B arms in flight at 21:34 on 08-03, while `d1b_blocked_gate_does_not_exercise_the_fix_20260804T0500Z.md`
describes D1-B being launched with `--biology-key-momentum 0.999` and failing its gate at 05:00 on
08-04. §5.4's before/after is written as a before/after across two launches and explicitly not as a
matched ablation, so nothing in the draft turns on the ordering — but if D1-B's provenance is ever
quoted as a controlled contrast, this needs resolving first.)*

### In plain terms

Our own fix is graded by the same ruler as everyone else's, and it fails: the rank gap behind it,
3.29×, is smaller than the 3.295× the same number moves when you retrain the same model. The reason we
still believe the fix is that the objective went from never once finishing training in a usable state
to finishing on all three seeds with a measurable, interval-backed signal. That is a yes/no, not a
ratio, and no noise floor touches it. Rank told us where to look; it did not tell us the repair worked.

### Files / commits

- `paper/P2_RANK_DRAFT.md` §5 (new §5.4), status block, abstract, §3.6, §4.7.3, §6.2, §7, Appendices A and C
- `paper/QUEUE_ANCHORING.md` superseding header
- `paper/P2_FIGURES.md` S4, "figures the paper does not have", pending dependencies
- Sources: `NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §1, §3;
  `d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`;
  `d1b_blocked_gate_does_not_exercise_the_fix_20260804T0500Z.md`;
  `PREDECLARED_retraining_envelope_20260804T0330Z.md`
