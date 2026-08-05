## 2026-08-05 00:00 UTC — The last three unjudgeable selections have floors on their own settings, and all three clear. One of them passes by 5.6% and the same ten runs undercut it

**Logged:** 2026-08-05 00:00 UTC. **How obtained:** twenty-five GPU runs on the A100
(`150.136.45.194`) from `~/ws_j` / `~/ws_j2`, workspaces built with
`git -c core.autocrlf=false archive HEAD` and verified **file by file against `git ls-tree`** —
584 and 585 files, 0 mismatches, and `d1_momentum_probe.py` byte-identical in both. Scored by
`v2/research/rebase/p2/p2_probe_floors.py`, which imports its statistic table, its fold, its
bimodality rule, α-ReQ's index range and LiDAR's ridge from `p2_envelope_floors.py` rather than
restating any of them. Outputs on persistent NFS at `~/e0_run/d1_probefloor600/` and
`~/e0_run/d1_capfloor/`, vendored under `v2/research/rebase/p2/figures/data/`.

### 0. What was missing, and what closed it

Three of §4.1a's sixty-two rows could not be judged, and each named the specific gap that stopped it
rather than repeating "no floor exists on this block" — which stopped being true when the probe block
got a floor. None of the three could be closed by stretching that floor:

| row | what stopped it |
|---|---|
| §5.4 row 1 | read at **step 600**, past the 500 the repeats were run to |
| §5.4 limit 2 | step 600 **and** an arm at **m = 0.99**, which is neither arm the repeats cover |
| §5.2 measurement 3 | an arm at **capacity 64**, which no repeat was run at |

**Each was closed by running the same measurement at the setting the row sits on**, five same-seed
repeats per arm, both arms of every comparison, seed 42 throughout — the variation being measured is
GPU non-determinism, not seed. **And the arms are now written into the block string**, exactly as the
reading step already was, so block-matching enforces them: a floor built from m = 0.999 and m = 0 may
not rule on a comparison between m = 0.999 and m = 0.99.

```
d1_momentum_probe.py {0.999, 0, 0.99} 0.04 600 4096 2e-4 42 <export_dir>     15 runs
d1_momentum_probe.py 0 0.04 200 {64, 4096} 2e-4 42 <export_dir>              10 runs
```

**A three-arm run is scored once per pair.** `combine()` takes the max over the arms it is given, so
scoring m = 0.999 / m = 0.99 / m = 0 in one go would judge §5.4 limit 2 against the collapsed arm's
spread, which is not one of its two sides. Two invocations over the same fifteen repeat directories,
two output files.

### 1. The three verdicts

| row | ratio | statistic | its own floor | carried by | verdict |
|---|---:|---|---:|---|---|
| §5.4 row 1 — m = 0.999 vs m = 0, one seed, **step 600** | **2.641×** | R3 | **1.749×** | m = 0 | **clears by 1.51×** |
| §5.4 limit 2 — **m = 0.999 over m = 0.99**, step 600 | **1.262×** | R3 | **1.195×** | m = 0.999 | **clears by 1.06×** |
| §5.2 measurement 3 — capacity 64 vs 4,096, **step 150** | **2.857×** | R3 | **1.705×** | capacity 4,096 | **clears by 1.68×** |

**Nothing about the direction was predeclared and no threshold moved.** Every verdict is
`ratio > floor` computed by `p2_floor_audit.check`, which re-resolves both sides of every comparison
and both ends of every floor from the files they were measured into.

### 2. Canonical R1 and R3 side by side, per arm, with the statistic and block named

Block: **fixed held-out probe, `wsi_biology` view, raw** — the block §5 quotes.

| floor | step | **R1** | R1 by arm | **R3** | R3 by arm |
|---|---:|---:|---|---:|---|
| m = 0.999 / m = 0 | 600 | **1.751×** | m0: 1.751× · m0.999: 1.155× | **1.749×** | m0: 1.749× · m0.999: 1.195× |
| m = 0.999 / m = 0.99 | 600 | **1.155×** | m0.99: 1.102× · m0.999: 1.155× | **1.195×** | m0.99: 1.152× · m0.999: 1.195× |
| capacity 64 / 4,096 | 150 | **1.662×** | cap4096: 1.662× · cap64: 1.112× | **1.705×** | cap4096: 1.705× · cap64: 1.224× |

**Both arms, and take the larger — the rule earns its keep twice more and fails to once.** In the two
comparisons that have a collapsed side, the collapsed arm carries the floor by roughly a factor of
1.5: **m = 0 gives 1.749× where m = 0.999 gives 1.195×**, and **capacity 4,096 gives 1.705× where
capacity 64 gives 1.224×**. Measuring the healthy arm alone would have published floors of 1.20× and
1.22× and flattered both rows. **The third comparison is the exception that shows what the rule is
actually about**: m = 0.999 against m = 0.99 has **no collapsed arm** — both train — so its floor is
the smallest of the sixteen probe floors at 1.195×, and it is carried by m = 0.999. The rule is not
"the collapsed arm is noisier". It is "measure both sides", and here that returns the stable one.

### 3. Per repeat, never a mean — canonical R1 / R3 on `wsi_biology`

**Step 600, capacity 4,096, lr 2e-4, decorrelation 0.04, seed 42.**

| rep | m = 0.999 | m = 0.99 | m = 0 |
|---|---|---|---|
| 1 | 10.939 / 6.971 | 6.935 / 5.631 | 2.738 / 2.265 |
| 2 | 11.955 / 8.054 | 6.913 / 5.407 | 2.456 / 2.188 |
| 3 | 10.696 / 7.205 | 6.462 / 5.142 | 2.129 / 1.956 |
| 4 | 10.353 / 6.741 | 7.124 / 5.922 | 1.564 / 1.295 |
| 5 | 10.973 / 7.348 | 6.805 / 5.532 | 2.246 / 1.895 |

**Step 150, m = 0, lr 2e-4, decorrelation 0.04, seed 42.**

| rep | capacity 64 | capacity 4,096 |
|---|---|---|
| 1 | 8.804 / 5.003 | 3.123 / 2.088 |
| 2 | 9.470 / 5.841 | 3.753 / 2.869 |
| 3 | 9.734 / 6.037 | 3.645 / 2.718 |
| 4 | 9.787 / 6.123 | 3.800 / 2.906 |
| 5 | 9.228 / 5.436 | 2.286 / 1.704 |

### 4. The middle row passes narrowly, and the same ten runs undercut it. Both are recorded

**§5.4 limit 2 is the value this project actually runs**, `m = 0.999` over `m = 0.99`, and §5.4 says
of it: *"The specific value this project runs was selected by a rank comparison our own rule cannot
license."* The rule now licenses it — **by 5.6%**. Three things have to be said about that pass and
they are all in the audit row's `rests_on`:

1. **It clears**: 1.262× against 1.195×, computed by the checker.
2. **The ten runs that measured the floor separate the two arms by only 1.138× worst case under
   R3** — m = 0.999's lowest repeat is 6.741 and m = 0.99's highest is 5.922 — **which is inside that
   floor.** The row passes on the particular single-seed draw §5.2's table happens to record and
   would not pass on the worst draw of five. That is not a contradiction: the row is a claim about
   two specific runs and the floor is the noise those runs are drawn from. It *is* a warning about
   how much the pass is worth.
3. **Under canonical R1 the same ten runs separate them cleanly**: worst case **1.453×** against a
   **1.155×** floor, and every m = 0.999 repeat (10.353–11.955) is above every m = 0.99 repeat
   (6.462–7.124). **So the arms are cleanly ordered under R1 and marginally so under R3** — which is
   §4.5's statistic-conditionality landing on the one hyperparameter value this project ships.

**What it does not change**: the binary training outcome supports *momentum against none*, not
`0.999` over `0.99`, and §5.4's third limit says so on grounds that have nothing to do with the
floor.

### 5. §5.2 measurement 3 — the reading step was recoverable, and the two harnesses agree

§6.2 said the capacity sweep's step *"was never recorded"* and its logs *"are not vendored"*. **Both
halves were wrong and the fix cost a file copy.** `qsweep_d0.04_cap64.log` prints **6.17 at step 150**
and `decorr_causal_0.04.log` prints **2.16 at step 150**, both at decorrelation 0.04 — the two numbers
§5.2 quotes, at one step, in the same row. What kept them out of the audit was the **header**: the
older `decorr_causal.py` harness prints a fifth column carrying the decorrelation *loss term*, and
`PROBE_HEADERS` refused it by design rather than guess which column was the rank. Naming the header
recovered the step, and row 50's two values now resolve from bytes instead of from the draft's own
prose. **Recorded in `floor_audit.json`'s `known_source_disagreements` as a prose correction for
§6.2's owner.** The second half of that §6.2 row was right and stands: knowing the step does not make
the row judgeable, because a floor at one capacity may not be borrowed for another. That needed the
GPU.

**The floor is measured with a different harness from the row, and the note travels with the row.**
The sweep was run with `decorr_causal.py`; the floor is `d1_momentum_probe.py` at the same momentum,
decorrelation, learning rate, seed, capacity and reading step. **Where the two overlap they agree**:
the sweep's capacity-4,096 arm reads R3 **2.16** at step 150 against a same-seed span of
**1.704–2.906** here. **The capacity-64 arm does not quite**: the sweep reads **6.17** where the five
repeats span **5.003–6.123**, so the published value sits **0.8% above the top of its own five-run
range.** It does not move the verdict (2.857× against 1.705×) and it is recorded rather than smoothed
over — a `floor_note` on the row says the harnesses are not the same file.

### 6. Shape, and every other statistic

**None of the sixteen probe floors is bimodal**, under any statistic, at any step, on any arm — which
extends the finding of 2026-08-04 16:20 from ten floors to sixteen and keeps §4.1's bimodality scoped
to the **exported** block, where four of five repeats agree to 2% and one lands at a third of them.

At the same step and block, the other statistics:

| | m0.999/m0, step 600 | m0.999/m0.99, step 600 | cap64/cap4096, step 150 |
|---|---:|---:|---:|
| canonical R1 | **1.751×** | **1.155×** | **1.662×** |
| R2 / R3 | **1.749×** | **1.195×** | **1.705×** |
| PR / PR_rownorm | 1.805× | 1.432× | 1.758× |
| RankMe | 1.263× | 1.155× | 1.530× |
| stable rank | 1.537× | 1.371× | 1.576× |
| α-ReQ \|α−1\| | 1.711× | 1.145× | **2.144×** |
| **hard numerical rank** | **1.000×** | **1.000×** | **1.000×** |
| LiDAR (`wsi_rna_paired`) | 1.320× | 1.187× | 1.515× |
| canonical R1, `rna_biology` | 1.937× | 1.104× | 1.930× |

**R2 and R3 coincide on this block, and the reason is a property of the block rather than of the
statistics.** `z_biology` is L2-normalised at the model's output, so every probe row already has unit
norm and R3's row normalisation is a no-op: measured directly on one state, R2 = 6.9711779953 and
R3 = 6.9711779832. The same is true of PR and PR_rownorm. **On the fixed held-out probe, R2 and R3 are
the same statistic** — which is worth stating in §3.1, because §3.1's whole point is that three
statistics travel under one name and here two of them collapse into one.

### 7. The audit, re-run

`p2_floor_audit.py --check` reports **no disagreements**. Six floors were added (R1 and R3 for each of
the three comparisons) and the three rows were closed.

| | before | after |
|---|---:|---:|
| selections | 25 | 25 |
| — do not clear | 13 | **13** |
| — clear | 9 | **12** |
| — **unjudgeable** | **3** | **0** |
| rows total | 62 | 62 |
| rows unjudgeable for want of a floor | 19 | **16** |

**Every selection between candidate configurations in this paper now has a floor measured on its own
statistic, block, reading step and arms.** Thirteen of the twenty-five still fail one. **Sixteen
non-selection rows are still unjudgeable**, and the audit has not stopped saying "we cannot tell": the
in-run training batch, the 16-patient gate batch and the 282-patient live checkpoint have no floor and
cannot get one without re-creating runs that were never exported, and §5.2a's four `direction` rows
sit at learning rates `1e-3` and `4e-5` where nothing has been measured — and where, as of today's
other entry, the co-measured statistic they are read against has a retraining spread wider than the
difference read off it.

**The generated blocks were regenerated, not edited**: §4.1a's audit table, its measured-floor table
and its counting sentence all come from `p2_floor_audit.py` and a test asserts the draft prints
exactly what the module renders. `summary_sentence` now has a zero branch, because "0 cannot be
judged at all because no floor has been measured" reads as though a floor were missing.

### 8. Prose the other agent owns, flagged not edited

1. **§6.2**, the capacity-sweep row: *"The step was never recorded and the sweep's own logs are not
   vendored"* — both false now; and the row's conclusion *"would not make the row judgeable"* has been
   overtaken by the measurement. The row should close.
2. **§6.2**, the step-600 row: *"A probe floor at step 600, and at momentum values other than 0 and
   0.999 — not measured, and it is what keeps §5.4 row 1 and §5.4 limit 2 unjudgeable."* Measured;
   both rows close.
3. **§5.4 limit 2**: *"a **1.26×** difference, smaller than every floor this paper has measured on any
   statistic, view or block, and on a block where no floor has been measured at all"* — no longer
   true in either clause. It is larger than the floor measured on its own two arms at its own step
   (1.195×), and §4 now carries a floor smaller than 1.26× on more than one block. The honest
   replacement is §4 of this entry, including the 1.138× worst-case separation that cuts the other
   way.
4. **§5.4 row 1's table cell** and the §4.1a prose *"Three remain **unjudgeable**, each naming the
   specific run that would settle it"* — the runs were named, run, and returned three passes.
5. **§4.1a's "What this costs and what it buys"** still says *"Nine selections do clear, eight of them
   §5's"*. It is twelve and eleven.
6. **§3.1** should record that R2 and R3 are the same statistic on the probe block, for the reason in
   §6 above.

### 9. Files

* Floors: `~/e0_run/d1_probefloor600/out/P2_PROBE_FLOORS_S600_m0999_m0.json`,
  `..._S600_m0999_m099.json`, `~/e0_run/d1_capfloor/out/P2_PROBE_FLOORS_CAP.json` — each carrying the
  sha256 of every probe state it read, every per-repeat value at every step, its `config` (the
  settings a floor exists only at) and its own `absent` list.
* The twenty-five run logs, vendored beside them.
* `v2/research/rebase/p2/p2_probe_floors.py` — configuration is now arguments, `--arm-kind capacity`
  added; the default path is unchanged and `P2_PROBE_FLOORS.json` regenerates as it stands.
* `v2/research/rebase/p2/floor_audit.json`, `v2/tests/test_p2_floor_audit.py`.

---

## CORRECTION APPENDED 2026-08-05 04:10 UTC — the middle row's pass does not survive ten repeats per arm, and this entry's illustrating example inverts

This entry is append-only, so nothing above has been altered. Read this block as binding wherever it
contradicts the text above. **Source:**
`NOTEBOOK_ENTRIES/limit2_breaks_at_ten_repeats_20260805T0330Z.md`, predeclared in full at
`NOTEBOOK_ENTRIES/PREDECLARED_probe_floor_n10_and_momentum_grid_20260805T0200Z.md` before any of the
new runs started.

**What is withdrawn.** Every clause below concerns **§5.4 limit 2 only**. §5.4 row 1 (2.641× against
1.749×) and §5.2 measurement 3 (2.857× against 1.705×) are untouched, as is everything this entry
says about the capacity floor, the reading step, and the two harnesses agreeing.

| where | withdrawn text |
|---|---|
| 1 (title) | "One of them passes by 5.6%" — it does not pass at ten repeats per arm |
| 45 (§1 table) | "§5.4 limit 2 … **1.195×** … carried by m = 0.999 … **clears by 1.06×**" |
| 59 (§2 table) | the R3 entry "**1.195×** · m0.99: 1.152× · m0.999: 1.195×" for the m = 0.999 / m = 0.99 floor |
| 66–69 (§2) | "its floor is the smallest of the sixteen probe floors at 1.195×, and it is carried by m = 0.999. The rule is not 'the collapsed arm is noisier'. It is 'measure both sides', and here that returns **the stable one**." |
| 93–110 (§4, heading and items 1–2) | "The middle row passes narrowly" / "**It clears**: 1.262× against 1.195×" / "The rule now licenses it — **by 5.6%**" |
| 8.3 | "It is larger than the floor measured on its own two arms at its own step (1.195×)" |
| 8.5 | "It is twelve and eleven" — it is **eleven and ten** |

**What replaces them.** Both arms were taken from five same-seed repeats to **ten**, everything else
identical (seed 42, `~/ws_j2`, one A100, ten concurrent). At n = 10:

| | n = 5 (this entry) | **n = 10** |
|---|---:|---:|
| the row's ratio — two specific published runs, fixed | 1.262× | 1.262× |
| **R3 floor, its own two arms, its own step** | 1.195× | **1.326×** |
| which arm carries it | m = 0.999 | **m = 0.99** |
| verdict under `ratio > floor`, the paper's own rule | clears by 5.6% | **DOES NOT CLEAR** |

**It is not one unlucky draw.** Repeats 6–10 scored alone, as an independent n = 5 never pooled with
the first five, give **1.2791×**, also above 1.262×; on the `rna_biology` view the n = 10 floor is
**1.2764×**, also above it. The entire increase is one m = 0.99 repeat reading R1 5.520 / R3 4.465
where the other nine span 6.462–7.124 / 5.142–5.922, and every predeclared exclusion route was
checked and none applies: it trained (`biology_contrastive` 7.575 against chance ln 80 = 4.382), it
is not collapsed (R1 5.520 against the m = 0 arm's 1.564–2.738), `p2_envelope_floors._shape` does not
call the arm bimodal, and repeats 6–10 overlap repeats 1–5, so it is not a batch effect. **Same-seed
GPU non-determinism at m = 0.99 produces a run about 20% low roughly one time in ten, and five
repeats did not see it.**

**§2's rule is vindicated; the example this entry chose to illustrate it with inverts.** "Measure
both sides" is exactly right, and this row is now the sharpest argument for it in the project — but
not for the reason §2 gives. §2 says the comparison with no collapsed arm "returns the stable one".
At ten repeats the **m = 0.99** arm is the noisier of the two and carries the floor, and it does so
on the strength of a single repeat in ten. **The lesson is not that measuring both sides returns the
stable arm. It is that which arm is "the stable one" is not a property you can read off five
repeats.** The floor is `max/min` and therefore non-decreasing in the repeat count, so its growth is
not itself a finding; that it grew *past the ratio of the row it governs* is.

**What this entry got right and is not withdrawn — the ordering, which got stronger.** §4 item 3's
"every m = 0.999 repeat above every m = 0.99 repeat" holds at ten per arm and under seven statistics
(R1, R2, R3, R1_rownorm, R1_uncentred, RankMe, and α-ReQ |α−1| with its sign applied): exact
one-sided permutation probability **1/184,756 = 5.4e-06**, up from 1/252 at five per arm. It remains
a statement about GPU non-determinism at a fixed seed and **not** a p-value for the momentum effect.
Canonical R1's worst-case separation of **1.453×** clears that statistic's own ten-repeat floor of
**1.291×** — the one test that had ten repeats per arm to break it and did not. And a five-point
momentum grid {0, 0.98, 0.99, 0.995, 0.999}, run afterwards, is **strictly monotone with complete
separation at every adjacent pair** under canonical R1, so the shipped comparison is two rungs of a
ladder rather than a step. **The trend and the ordering hold. The floor-clearing claim for exactly
this pair does not.**

**Applied to the paper.** `floor_audit.json`'s floor `R3_probe_step600_m0999_vs_m099` now reads
1.3263× from the ten-repeat file (its R1 twin reads 1.2906× from the same runs), row
`5.4-m0999-over-m099` records `clears: false`, §4.1a's generated table and counting sentence were
regenerated to **14 fail / 11 clear / 0 unjudgeable** of 25 selections, §5.4 limit 2 and Appendix C
were rewritten around §0–§1 and §5 of the source entry, and the draft's Status log gained **item 14**
recording the reversal rather than editing item 12.
