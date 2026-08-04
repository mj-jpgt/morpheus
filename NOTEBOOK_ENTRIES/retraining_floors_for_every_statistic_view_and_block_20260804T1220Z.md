## 2026-08-04 13:00 UTC — The criterion had a floor for one statistic on one block. Measuring the rest: it is a property of the statistic (1.000×–3.295×), of the view (1.019× vs 3.295×), and RankMe as published beats us on it

**Logged:** 2026-08-04 12:20 UTC. **How obtained:** CPU only, thread-capped
(`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`). The GPU was at
0% and was not touched. Every rank statistic is **imported** — R1/R2/R3 and the channel from
`v2/calibra/spectral.py`, RankMe / participation ratios / stable rank / α-ReQ / LiDAR from
`v2/research/rebase/p2/p2_competing_metrics.py`, the hard rank from `numpy.linalg.matrix_rank`;
nothing computed inline. Run from `~/ws_floor/morpheus`, built with
`git -c core.autocrlf=false -c core.eol=lf archive HEAD` at commit `146d9d6` and verified per file by
git blob SHA-1 against `git ls-tree -r HEAD` **before** execution: **543/543 files, 0 missing, 0
extra, 0 differing**. One tarball, no per-file `scp`. No other agent's workspace was read or written.

### 1. Why — the scope error, stated plainly

`P2_RANK_DRAFT.md` §4.1's criterion is that a rank difference smaller than the measured same-seed
retraining floor is not resolvable. **The floor had been measured for canonical R1 on the exported
`wsi_biology` block and for nothing else** — 3.295× residualised / 3.111× raw, five identical
`programme_only` retrains at seed 42, bimodal.

§4.1a then applied that one number to fifty comparisons. Thirteen sat on a statistic (R2, R3, PR,
RankMe, participation ratio, stable rank, α-ReQ, LiDAR, hard rank) or a block (the fixed held-out
probe, in-run training batches, the 16-patient gate batch, the `rna_biology`/`full_biology` views)
with **no floor at all** — and were nonetheless printed with a verdict of "**no**". Eleven more were
held-out-probe measurements judged against an exported-artifact floor. Eight of T1's twelve metric
rows and every rank number in §5 were in that position.

**A comparison with no measured floor has not failed the criterion; the criterion has not been
applied to it.** Printing those two outcomes in one column flattered us: it made an unmeasured ruler
look like a failed test. This is the same objection the paper makes to RankMe — a criterion applied
outside the scope of the measurement that licenses it — made to us, and it had been true of us for as
long as §4.1a existed.

### 2. What was measured, and what it cost (nothing)

The five repeats were **already exported** and each `.npz` carries all three co-trained views. So
every missing floor except the block-level ones is a re-derivation, not a re-run.
`v2/research/rebase/p2/p2_envelope_floors.py`: **ten statistics × three views × {raw, residualised}**,
plus LiDAR over the `wsi`/`rna` positive pair and the top-CCA channel on the same five runs.
Output vendored to
`v2/research/rebase/p2/figures/data/ws_floor/out/P2_ENVELOPE_FLOORS.json` through
`extract_from_box.py`, the only sanctioned path; **every previously vendored file re-hashed
byte-identical**.

**No source disagreement.** The recomputation from the artifacts reproduces the readout log the
paper's numbers are parsed from, to four decimal places at both extremes: rank_raw
8.0326 / 24.9895 (3.111×), rank_residualised 8.8340 / 29.1057 (3.295×), channel 0.5859 / 0.6182
(1.055×), and every per-repeat value in §4.1's table. **3.295× and 3.111× now have two independent
sources**, and `p2_floor_audit.check()` fails if they ever part rather than preferring one.

#### 2.1 The floors — five values, min, max, fold, and whether the bimodal shape holds

**Caveat, carried per statistic and not negotiable:** each of these is **n = 5**, one arm, one seed,
one stack, no interval, and it is a **FLOOR TWICE OVER** — `programme_only` is this project's stable
arm, and same-seed repeats exclude seed variation entirely (§4.2 measures that other axis as the
larger term). These are not estimated distributions. The right sentence is "a floor of 2.290×
measured on five same-seed repeats", never "R3 varies 2.29×".

**Bimodality rule, fixed in the module and not chosen per statistic:** the *divergent run* is the one
whose removal minimises the remaining four's fold; *rest* is that four's own fold; **bimodal** =
rest ≤ 1.05 **and** full fold ≥ 2 × rest. Under this rule §4.1's own floor reads exactly as §4.1
describes it (rep2, four others within 2.8%, separation 3.21×).

**`wsi_biology`, residualised — the block §4.1 measures**

| statistic | rep1 | rep2 | rep3 | rep4 | rep5 | min | max | fold | divergent | rest | bimodal |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|:--:|
| R1 | 28.3202 | 8.8340 | 28.3482 | 29.1057 | 28.9588 | 8.8340 | 29.1057 | **3.295×** | rep2 | 1.028× | **YES** |
| RankMe | 28.3248 | 8.8359 | 28.3529 | 29.1105 | 28.9635 | 8.8359 | 29.1105 | **3.295×** | rep2 | 1.028× | **YES** |
| α-ReQ \|α−1\| | 1.2800 | 2.6247 | 1.2399 | 1.1419 | 1.1935 | 1.1419 | 2.6247 | **2.299×** | rep2 | 1.121× | no |
| R3 | 14.1768 | 6.3663 | 14.2608 | 14.5795 | 14.3560 | 6.3663 | 14.5795 | **2.290×** | rep2 | 1.028× | **YES** |
| R2 | 12.2858 | 5.5972 | 12.1480 | 12.4506 | 12.2726 | 5.5972 | 12.4506 | **2.224×** | rep2 | 1.025× | **YES** |
| α-ReQ α | 2.2800 | 3.6247 | 2.2399 | 2.1419 | 2.1935 | 2.1419 | 3.6247 | **1.692×** | rep2 | 1.064× | no |
| PR_rownorm | 5.4474 | 3.7152 | 5.3433 | 5.4407 | 5.3176 | 3.7152 | 5.4474 | **1.466×** | rep2 | 1.024× | no |
| PR | 4.4913 | 3.1643 | 4.2728 | 4.4225 | 4.2798 | 3.1643 | 4.4913 | **1.419×** | rep2 | 1.051× | no |
| stable rank | 2.6898 | 2.1977 | 2.5520 | 2.6419 | 2.5584 | 2.1977 | 2.6898 | **1.224×** | rep2 | 1.054× | no |
| hard rank | 256 | 256 | 256 | 256 | 256 | 256 | 256 | **1.000×** | — | 1.000× | no |
| *channel (top-CCA 16)* | 0.6182 | 0.5859 | 0.6123 | 0.6110 | 0.6098 | 0.5859 | 0.6182 | **1.055×** | rep2 | 1.014× | no |

**`wsi_biology`, raw**

| statistic | rep1 | rep2 | rep3 | rep4 | rep5 | min | max | fold | divergent | rest | bimodal |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|:--:|
| R1 | 24.4806 | 8.0326 | 24.5040 | 24.9895 | 24.9116 | 8.0326 | 24.9895 | **3.111×** | rep2 | 1.021× | **YES** |
| R2 | 10.8603 | 5.1004 | 10.7676 | 10.9157 | 10.8102 | 5.1004 | 10.9157 | **2.140×** | rep2 | 1.014× | **YES** |
| R3 | 10.8603 | 5.1004 | 10.7676 | 10.9157 | 10.8102 | 5.1004 | 10.9157 | **2.140×** | rep2 | 1.014× | **YES** |
| α-ReQ \|α−1\| | 1.3972 | 2.6304 | 1.3685 | 1.2291 | 1.2707 | 1.2291 | 2.6304 | **2.140×** | rep2 | 1.137× | no |
| RankMe (published) | 3.5732 | 1.9883 | 3.5318 | 3.6011 | 3.5401 | 1.9883 | 3.6011 | **1.811×** | rep2 | 1.020× | no |
| α-ReQ α | 2.3972 | 3.6304 | 2.3685 | 2.2291 | 2.2707 | 2.2291 | 3.6304 | **1.629×** | rep2 | 1.075× | no |
| PR = PR_rownorm | 3.7480 | 2.5913 | 3.6464 | 3.6578 | 3.5796 | 2.5913 | 3.7480 | **1.446×** | rep2 | 1.047× | no |
| stable rank | 2.2082 | 1.7729 | 2.1729 | 2.1692 | 2.1338 | 1.7729 | 2.2082 | **1.246×** | rep2 | 1.035× | no |
| hard rank | 256 | 256 | 256 | 256 | 256 | 256 | 256 | **1.000×** | — | 1.000× | no |

**`rna_biology`, residualised** — nothing is bimodal and nothing exceeds 1.032×

| statistic | rep1 | rep2 | rep3 | rep4 | rep5 | min | max | fold | divergent | rest | bimodal |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|:--:|
| R1 | 27.3131 | 27.2245 | 27.5654 | 27.7497 | 27.7267 | 27.2245 | 27.7497 | **1.019×** | rep2 | 1.016× | no |
| RankMe | 27.3178 | 27.2292 | 27.5701 | 27.7545 | 27.7315 | 27.2292 | 27.7545 | **1.019×** | rep2 | 1.016× | no |
| R2 | 15.3483 | 15.2687 | 15.4966 | 15.5923 | 15.6221 | 15.2687 | 15.6221 | **1.023×** | rep2 | 1.018× | no |
| R3 | 16.2409 | 16.1913 | 16.3968 | 16.5703 | 16.5713 | 16.1913 | 16.5713 | **1.023×** | rep2 | 1.020× | no |
| α-ReQ \|α−1\| | 3.6959 | 3.6705 | 3.6362 | 3.6774 | 3.7186 | 3.6362 | 3.7186 | **1.023×** | rep3 | 1.013× | no |
| stable rank | 3.9529 | 3.8939 | 3.9985 | 3.9487 | 3.9667 | 3.8939 | 3.9985 | **1.027×** | rep2 | 1.013× | no |
| PR | 6.8638 | 6.7747 | 6.9503 | 6.9395 | 6.9832 | 6.7747 | 6.9832 | **1.031×** | rep2 | 1.017× | no |
| PR_rownorm | 7.4323 | 7.3563 | 7.5288 | 7.5512 | 7.5928 | 7.3563 | 7.5928 | **1.032×** | rep2 | 1.022× | no |
| hard rank | 256 | 256 | 256 | 256 | 256 | 256 | 256 | **1.000×** | — | 1.000× | no |

**`rna_biology`, raw:** R1 **1.023×**, R2 = R3 1.028×, RankMe 1.030×, PR = PR_rownorm 1.032×,
α-ReQ |α−1| 1.026×, stable rank 1.022×, hard rank 1.000×. Nothing bimodal.

**`full_biology`, residualised:** R1 **1.020×**, RankMe 1.020×, R2 1.006×, R3 1.019×,
α-ReQ |α−1| 1.058×, stable rank 1.017×, PR 1.006×, PR_rownorm 1.009×, hard rank 1.000×.
**`full_biology`, raw:** R1 **1.014×**, R2 = R3 1.006×, RankMe 1.005×, PR = PR_rownorm 1.018×,
α-ReQ |α−1| 1.058×, stable rank 1.023×, hard rank 1.000×. Nothing bimodal on either.

**LiDAR** (the `wsi`/`rna` positive pair, q = 2, δ = 1e-4): residualised 38.7877 / 41.1002 =
**1.060×** (rep2 low, rest 1.024×, not bimodal); raw 39.6212 / 40.9816 = **1.034×** (rep2 low, rest
1.026×, not bimodal).

### 3. What the floors say — three results, two of them against us

**(a) The floor is a property of the STATISTIC, and the spread between statistics is the size of the
effect the paper is about.** On one block — exported `wsi_biology`, residualised, the *same five
runs* — it runs from **1.000×** (hard numerical rank, which does not move at all) through 1.224×
(stable rank), 1.419× (PR) and 2.224× / 2.290× (R2 / R3) to **3.295×** (R1). §4.3's heading said the
floor is a property of the arm *and not of the statistic*; that half is now measurably wrong and is
**withdrawn** in §4.3. Concretely: judging an R3 comparison against R1's floor — which §4.1a did on
fourteen rows — is **1.4× too strict**.

**(b) RankMe as published is the more reproducible statistic on our own artifacts.** Raw block:
**1.811×** against canonical R1's **3.111×**, same five runs, and it is not bimodal (separation 1.78
against R1's 3.05). The mechanism is the one §4.6 already names for RankMe's D2 advantage: the
uncentred normalisation retains the mean-offset direction, and **every row of every exported view has
L2 norm exactly 1.000** (verified directly on `rep1.npz`: min = max = 1.0 for all three views), so
that direction is both large and stable. On the *residualised* block the column mean is gone and the
two statistics coincide — floors of 3.295× against 3.295×, levels agreeing to 0.02%. **This is a
result against our own instrument and it is now the one selection in the whole audit that clears a
floor its own statistic and block license.**

A related structural check that falls out of the same table: **R2 and R3 agree to float noise on the raw block (≤ 1e-8 relative), and so do
PR and PR_rownorm**, for all three views, because R3 is R2 on L2-normalised rows and the rows
are already unit-norm. They separate on the residualised block, where residualisation destroys the
norm. Two statistics agreeing exactly is otherwise the signature of a bug, so it is asserted in the
test with its reason.

**(c) The floor is a property of the VIEW, and the divergent run is divergent in only one of them.**
Same five artifacts, same statistic, same runs: **3.295×** on `wsi_biology`, **1.019×** on
`rna_biology`, **1.020×** on `full_biology`. Repeat 2 — whose WSI-view rank is a third of its
siblings' — sits within **1.6%** of them on the RNA view and is the *high* member on the full view.
**The catastrophic one-in-five is not a property of the run. It is a property of that run's WSI
encoder.** Consequence for §4.5(c): every one of its twelve `rna_biology` / `full_biology` arm
comparisons is resolvable against its own view's floor, and none of the six on `wsi_biology` is.

### 4. Is the bimodal shape statistic-dependent? Yes — but the divergence is not hidden from anything

This was asked as a direct test of the paper's thesis, so both halves of the answer matter.

**Every statistic that moves at all identifies repeat 2 as the outlier, in the degradation
direction** — lower rank for the eight rank-type statistics, *higher* α for α-ReQ (a steeper
eigenspectrum decay, which is α-ReQ's own bad direction). On the WSI block that is 10 of 10 in both
raw and residualised, hard rank excepted because it does not move. **There is no statistic in which
the divergent run looks ordinary**, which would have been the larger and more disruptive finding.

**What is statistic-dependent is the shape and the magnitude.** The four-concordant-plus-a-factor
signature §4.1 calls bimodal survives under **R1, R2, R3 and residualised RankMe** and under nothing
else. Raw RankMe misses it narrowly (rest 1.020×, separation 1.776 against a 2.0 bar) and is reported
with its numbers rather than on the boolean. The same run is 1.22× away under stable rank, 1.06×
under LiDAR, 1.00× under the hard rank.

**The reading.** The entropy-based effective ranks (R1, R2, R3, RankMe) are functions of the *whole*
normalised spectrum and are dominated by the tail; stable rank, PR and LiDAR are dominated by the top
of the spectrum. That the first family moves by a factor while the second moves by percent says the
divergence is a **redistribution of spectral mass in the tail, not a change of the dominant
subspace**. That is a narrower and better-supported claim than "rank is unreliable", and it is the
one the draft now carries.

### 5. What is NOT recoverable, named — and what §5 would need

The five repeats were **exported, not probed**. `~/e0_run/d1_envelope/` holds `rep{n}.npz` and a
per-run `train_metrics.jsonl`; neither carries a probe forward pass, a training batch's activations,
or a gate batch. So **no floor exists, and none can be recovered from these files, for**:

| block | why not, and what it would cost |
|---|---|
| **fixed held-out probe** | Needs the five repeats re-run with `d1_momentum_probe.py` attached — GPU, five runs. **Every rank number in draft §5 is on this block**, as are §4.9a's decorrelation ablation and §4.4(3)'s probe repeat. |
| **training batch, in-run** | The tripwire reads R3 on the live batch; those activations are never saved (`F3_TRIPWIRE_STEP200_R3_n5.json`: "states never saved, so [NOT RECOMPUTABLE]"). Needs the repeats re-run with the tripwire logged at a fixed step — GPU. |
| **16-patient gate batch** | Constructed inside the gate run, never exported. Needs five identical gate runs — GPU. |
| **282 held-out patients, live checkpoint** | Read off a D1-A checkpoint that was never exported; §4.9 records it `[NOT RECOMPUTED]` under R1. Needs the checkpoint re-created. |

This list lives in `P2_ENVELOPE_FLOORS.json`'s `absent_blocks` with the cost attached, and
`floor_audit.json`'s own list is asserted equal to it by a test, so the absence is a **recorded
result** rather than a silence.

**What §5 would need in order to become judgeable at all is one specific run: five same-seed repeats
of the `programme_free` / 500-step configuration with `d1_momentum_probe.py` attached, read at a
fixed step.** §6.2 now names it. Until it exists §5's rank comparisons are neither confirmed nor
refuted by this paper's criterion — they are outside its reach.

### 6. The audit re-run — the counts, bad news first

`p2_floor_audit.py --check`, zero disagreements. **56 rows** (§4.5(a)'s thirty and §4.5(c)'s eighteen
were single rows spanning five statistics and three views respectively; now that each has a floor of
its own they are split into five and three rows, judged separately).

| | before | after |
|---|---:|---:|
| selections between candidate configurations | 25 | 25 |
| **fail a floor their own statistic and block license** | **23** | **13** |
| **clear it** | **2** | **1** |
| **unjudgeable — no floor on this block** | **0** | **11** |
| rows with no measured floor (of all rows) | 13 / 50 | 25 / 56 |

**The failure count fell from 23 to 13, and not one of the ten that left the failing column moved
into the clearing one.** They moved to *unjudgeable*, which is the honest place for them: they are
probe-block, in-run-batch, gate-batch and live-checkpoint measurements, and the criterion cannot
reach them in either direction.

**And the two rows the paper reported as clearing were both on the block with no floor.** §4.4(3)'s
fixed-seed probe repeat (3.495×, quoted as "clears by 6%") and §5.2's step-400 fold (3.596×) are both
R3 or R1 on the *fixed held-out probe*, judged against an exported-artifact floor. They cleared a
floor that does not license them. **The one selection in this paper that clears a floor its own
statistic and block license is RankMe as published** — and only 1 of its 3 D2 pairs does (s43 at
3.382× against a 1.811× floor; s42 is 1.677× and s44 is 1.248×), so its 3/3 count still rests on two
unresolvable orderings.

**One more thing the split rows show: the resolvability verdict is itself under-determined.** Scored
against each statistic's own floor, the number of §4.5(a)'s six arm pairs that are resolvable is
**0 under R1, 1 under R2, 1 under R3, 2 under PR, 2 under PR_rownorm**; scored against each view's own
floor it is **0 on `wsi_biology`, 6 on `rna_biology`, 6 on `full_biology`**. §4.5's thesis is that the
*verdict* flips with the statistic, the block and the view. So does the question of whether there is a
verdict to be had.

### 7. What was enforced, so this cannot recur

- **`clears: null` ⟺ `floor: null`**, checked. For as long as a floor-less row could record
  `clears: false`, rows with no ruler were counted beside rows that had failed a measured one.
  `summary()` is now a three-way partition and `render_markdown` prints **unjudgeable**, not
  "**no floor**" in a column of verdicts.
- **Every floor is re-resolved from `P2_ENVELOPE_FLOORS.json`** — the recorded value against
  `max/min` from the file it was measured into, at tolerance 0.0015.
- **The two floors the paper's criterion rests on carry a second, independent source**, and a
  disagreement is reported (`STOP, and report both`) rather than resolved.
- **The recorded bimodality of each floor is checked against its source**, because whether the shape
  survives a change of statistic is a claim, not presentation.
- Negative tests for each of the above, plus the pre-existing block-mismatch one.

### 8. §5 — not touched, and one thing that needs correcting there

**§5 is another agent's section and is byte-identical to `cf1b0cf`** (asserted, not assumed). Two
items for whoever owns it:

1. **§5.4 says the momentum fix's rank difference "is not resolvable" by this paper's own criterion.
   That is stronger than the evidence now supports.** The comparison is on the fixed held-out probe,
   which has no measured floor and cannot get one from the exports, so the accurate word is
   **unjudgeable** — the criterion has not been applied to it, in either direction. §5.4's conclusion
   is unaffected because it rests on a binary training outcome, not on the ratio. Flagged in §6.2.
2. The same applies to every other rank comparison in §5 (§5.1's instance 3 sweep, §5.2's five
   measurements, §5.4's four rows): all are probe-block or gate-batch and all are now recorded as
   unjudgeable in `floor_audit.json`.

### 9. Suite

**426 passed** (baseline 399), thread-capped, `--basetemp` redirected. New:
`v2/tests/test_p2_envelope_floors.py` (14) and 13 added to `v2/tests/test_p2_floor_audit.py`.
`test_p2_figures.py` needs matplotlib, present locally; **nothing was installed into `~/venv`**.

### 10. In plain terms

We had one ruler — one way of measuring rank, on one kind of matrix — and we had been holding
everything in the paper up against it, including things it was never made to measure. The five
retrained models we already had on disk let us build the missing rulers for free, and they say three
things. Different ways of measuring rank disagree about how noisy rank is by a factor of three. The
published metric we spend the paper criticising is *more* reproducible on our own data than ours is.
And the one model out of five that came out wrong only came out wrong in its image half — its RNA half
was fine — so the failure we built the whole argument on is a failure of one encoder, not of the run.
For eleven of our own comparisons there is still no ruler at all and we cannot make one without the
GPU; those are now marked as unmeasured rather than as failed, which is a smaller-sounding claim and a
truer one.

### Files / commits

- New: `v2/research/rebase/p2/p2_envelope_floors.py`, `v2/tests/test_p2_envelope_floors.py`,
  `v2/research/rebase/p2/figures/data/ws_floor/out/{P2_ENVELOPE_FLOORS.json,floors_run.log}`
- Edited: `v2/research/rebase/p2/floor_audit.json` (26 floors, 56 rows),
  `p2_floor_audit.py` (three-way count, `--floors`, the coupling and cross-check enforcement),
  `v2/tests/test_p2_floor_audit.py`, `figures/extract_from_box.py` (two vendored files; a manifest
  bug that listed `.gitattributes` as box evidence and rewrote the file as CRLF on every Windows run)
- Draft: `paper/P2_RANK_DRAFT.md` §4.1a (both tables generated), §4.3 (heading claim withdrawn),
  §6.2 (one row closed, one row opened), status-block item 6, Appendix C. **§5 untouched.**
  `paper/P2_FIGURES.md` T8.
- Box: `~/ws_floor/morpheus` at `146d9d6`, 543/543 verified; `~/ws_floor/out/`.
- Sources: `NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §1;
  `p2_floor_audit_and_decorrelation_dissociation_20260804T2330Z.md`;
  `PREDECLARED_retraining_envelope_20260804T0330Z.md`.
