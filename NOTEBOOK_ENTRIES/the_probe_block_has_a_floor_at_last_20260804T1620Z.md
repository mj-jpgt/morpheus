## 2026-08-04 16:20 UTC — The fixed held-out probe has a retraining floor. It is carried by the collapsed arm, it is not bimodal, and eight of §5's eleven unjudgeable selections clear it

**Logged:** 2026-08-04 16:20 UTC. **How obtained:** ten GPU runs on the A100 (`150.136.45.194`) from
`~/ws_pf`, a workspace built with `git -c core.autocrlf=false archive HEAD` and verified **file by
file against `git ls-tree`** — 568 files, 0 mismatches. Scored by
`v2/research/rebase/p2/p2_probe_floors.py`, which imports its statistic table, its fold, its
bimodality rule, α-ReQ's index range and LiDAR's ridge from `p2_envelope_floors.py` rather than
restating any of them. Outputs on persistent NFS at `~/e0_run/d1_probefloor/`, vendored to
`v2/research/rebase/p2/figures/data/e0_run/d1_probefloor/`.

### 0. Why this run and not another

Draft §6.2 named it, in these words:

> **What would close it is one specific run: five same-seed repeats of the `programme_free` /
> 500-step configuration with `d1_momentum_probe.py` attached, read at a fixed step.**

Every rank number in §5 is measured on the fixed held-out probe, and no floor had ever been measured
for that block — so all of §5 was `unjudgeable` by the paper's own criterion, which is neither
passing nor failing. `p2_envelope_floors.ABSENT_BLOCKS` priced it at "GPU, five runs".

### 1. What was run

```
d1_momentum_probe.py {0.999,0.0} 0.04 500 4096 2e-4 42 <export_dir>
```

**Ten runs, not five: five identical repeats of *each* of the two arms §5 compares.** Same momentum,
decorrelation, capacity, learning rate, **same seed 42**, same 500-step budget, same workspace, one
A100. The only difference between the five repeats of an arm is GPU non-determinism — all ten begin
at the same verified initialisation, canonical R1 **101.38** and R3 **67.55**, which is also the
initialisation `mseed_m0.999_s42.log` and every `ablate_*` log begin at.

It is a **floor twice over** in exactly the sense §4.1 uses of the exported-block floor: same-seed
repeats exclude seed variation entirely, and §4.2 measures seed as the dominant term. n = 5 per arm,
no interval, one seed, one stack. Not a distribution.

`export_dir` is a new, purely additive seventh argument to the harness (commit `0688980`): the probe
states `geometry()` already computes are written beside the log at every step it reads. The printed
columns are untouched, so these logs are comparable to the six `mseed_*` logs written without it.
Without it only R1 and R3 could ever have had a floor here, because the states were being thrown
away.

### 2. Per repeat, never a mean — canonical R1 / R3 on `wsi_biology`, fixed held-out probe

**m = 0.999**

| rep | step 100 | 200 | 250 | 400 | 500 |
|---|---|---|---|---|---|
| 1 | 13.95 / 7.35 | 12.49 / 7.81 | 11.74 / 7.15 | 12.14 / 7.80 | **11.70 / 7.74** |
| 2 | 14.00 / 7.47 | 12.10 / 7.15 | 12.95 / 7.97 | 12.38 / 7.95 | **11.69 / 7.64** |
| 3 | 13.28 / 7.07 | 12.15 / 7.09 | 11.11 / 6.93 | 11.25 / 7.20 | **11.05 / 7.45** |
| 4 | 13.65 / 7.20 | 12.80 / 7.69 | 11.92 / 7.30 | 11.80 / 7.61 | **11.18 / 7.53** |
| 5 | 12.59 / 6.67 | 11.19 / 6.87 | 11.42 / 6.91 | 11.90 / 7.91 | **10.74 / 7.02** |

**m = 0**

| rep | step 100 | 200 | 250 | 400 | 500 |
|---|---|---|---|---|---|
| 1 | 4.00 / 2.83 | 3.73 / 2.99 | 3.73 / 3.00 | 2.84 / 2.16 | **2.30 / 2.02** |
| 2 | 3.36 / 2.14 | 2.20 / 1.69 | 1.93 / 1.47 | 2.57 / 1.98 | **2.35 / 1.93** |
| 3 | 3.00 / 1.96 | 2.10 / 1.51 | 2.23 / 1.74 | 2.16 / 1.85 | **1.84 / 1.52** |
| 4 | 3.55 / 2.52 | 2.37 / 1.97 | 2.24 / 2.00 | 2.51 / 2.06 | **1.84 / 1.53** |
| 5 | 3.98 / 2.93 | 1.81 / 1.46 | 3.17 / 2.38 | 1.81 / 1.49 | **1.72 / 1.33** |

### 3. The floor, canonical R1 and R3 side by side, per reading step

Block: **fixed held-out probe, `wsi_biology` view, raw** — the block §5 quotes, the one on which the
rank statistics centre internally and no confound residualisation is applied.

| step | **canonical R1** floor | m = 0.999 alone | m = 0 alone | **R3** floor | m = 0.999 alone | m = 0 alone |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | **1.333×** | 1.112× | 1.333× | **1.494×** | 1.119× | 1.494× |
| 200 | **2.057×** | 1.143× | 2.057× | **2.041×** | 1.137× | 2.041× |
| 250 | **1.933×** | 1.165× | 1.933× | **2.035×** | 1.153× | 2.035× |
| 400 | **1.570×** | 1.101× | 1.570× | **1.449×** | 1.105× | 1.449× |
| 500 | **1.367×** | 1.089× | 1.367× | **1.516×** | 1.103× | 1.516× |

**The collapsed arm carries the floor at every step and under both statistics, by roughly a factor of
two.** Had this been measured on one arm the way §4.1's exported floor was, it would have read
1.09×–1.17× and would have flattered every row in the audit. Measuring both arms is not thoroughness;
it is the difference between a floor of 1.1× and one of 2.1×.

### 4. Every other statistic, at step 500, same block

| statistic | `wsi_biology` | `rna_biology` |
|---|---:|---:|
| canonical R1 | **1.367×** | 1.271× |
| R2 / R3 | **1.516×** | 1.428× |
| PR / PR_rownorm | 1.594× | 1.601× |
| RankMe | 1.300× | 1.292× |
| stable rank | 1.339× | 1.352× |
| α-ReQ (α) | 1.372× | 1.727× |
| α-ReQ \|α−1\| | 1.490× | **2.270×** |
| **hard numerical rank** | **1.000×** | **1.000×** |

LiDAR is a per-pair quantity rather than a per-view one and is recorded under its own pseudo-view
`wsi_rna_paired`, so that a paired-block floor can never be matched to a single-view comparison by
accident: **1.447×** (m = 0 arm; 1.163× on the m = 0.999 arm).

The hard numerical rank is pinned at 256 in all ten runs, which is the probe size — the same
degenerate behaviour it shows on the exported block, and the reason §4.9's "16/16" instance is not
evidence of anything. Every other statistic on this block spreads by 1.3×–2.3× between identical
same-seed retrains.

### 5. Is it bimodal? **No — and that is a result**

§4.1's exported-block floor is bimodal: four of five repeats agree to 2% and rep2 lands at a third of
them. **Nothing of that shape appears on the probe block**, under either statistic, at any of the ten
steps, on either arm. Under the rule fixed in `p2_envelope_floors._shape` (`concordant` = remaining
four within 1.05×; `bimodal` = concordant and separation ≥ 2.0) every one of the twenty probe floors
is `not bimodal`, with remaining-four folds of 1.17×–1.64× — nowhere near the 1.05× that a
four-run agreement requires. The divergence here is **graded across all five runs**, not one run
falling off a cliff.

The two shapes are therefore not the same phenomenon, and a reader who was told the retraining
envelope is bimodal should be told that this holds on the exported artifact and not on the probe.

### 6. The audit, re-run

`p2_floor_audit.py --check` reports no disagreements. Ten floors were added (R1 and R3 at steps 100,
200, 250, 400, 500), and the reading **step is written into the block string** so that
block-matching enforces it — a floor at step 500 may not be applied to a reading at step 400 any more
than a raw floor may be applied to a residualised ratio.

| | before | after |
|---|---:|---:|
| selections | 25 | 25 |
| — do not clear | 13 | **13** |
| — clear | 1 | **9** |
| — **unjudgeable** | **11** | **3** |
| rows total | 60 | 62 |
| rows unjudgeable for want of a floor | 29 | **19** |

**Eight of the eleven became judgeable, and all eight clear.** Nothing about the direction was
predeclared; every verdict is `ratio > floor` computed by the checker.

| row | comparison | ratio | floor on its own block and step | verdict |
|---|---|---:|---|---|
| §4.4(3) | m = 0.999 vs m = 0, seed fixed, step 200 | 3.495× | R3, probe, step 200 — 2.041× | **clears** |
| §5.4 row 2 | m = 0.999 vs m = 0, worst of three seeds, step 500 | 3.286× | **R1**, probe, step 500 — 1.367× | **clears** |
| §5.4 row 3 | the same under the tripwire statistic | 2.438× | R3, probe, step 500 — 1.516× | **clears** |
| §5.2 prose | the widest per-step fold, step 400 | 3.596× | R3, probe, step 400 — 1.449× | **clears** |
| §5.2 meas. 2 | m = 0.999 vs no momentum, step 100 | 2.671× | R3, probe, step 100 — 1.494× | **clears** |
| §5.2 turnover | m = 0.95 → 0.999 at capacity 4,096, step 250 | 2.687× | R3, probe, step 250 — 2.035× | **clears** |
| §4.9a | decorrelation 0.0 → 0.04, step 400 (R3) | 1.854× | R3, probe, step 400 — 1.449× | **clears** |
| §4.9a | the same three runs under canonical R1 | 1.940× | R1, probe, step 400 — 1.570× | **clears** |

**The one that matters most is §5.4 row 2.** It is the momentum fix, the paper's own intervention. It
failed the exported-block R1 floor of 3.295× by 0.3% — which §5.4 records at length as the paper
applying its criterion to itself and losing. On the block it is actually measured on, at the step it
is actually read at, it clears by a factor of 2.4. Both statements are true and the second is the one
that is statistic-, block- and step-matched.

Two rows moved the other way and are recorded rather than buried. §5.2's "the best-agreeing arm does
not have the best rank" (1.036×) is **inside** the step-100 floor of 1.494×, which is what §4.1a
already said of it — an equality claim may be read off it, an ordering may not. And the n = 2
like-for-like pair §5.4 said did not exist (1.066×) is inside the n = 5 spread it is now superseded
by, which is what an n = 2 sample of the same noise should do.

### 7. The three that remain unjudgeable, and exactly what stops each

The block has a floor now; these three sit outside it, and each says why rather than repeating "no
floor exists on this block", which is no longer true.

* **§5.4 row 1** (2.641×) — read at **step 600**, past the 500 the repeats were run to. Five more
  runs at a 600-step budget close it.
* **§5.4 limit 2** (1.262×) — step 600 **and** one arm at **m = 0.99**, which is neither arm.
* **§5.2 measurement 3** (2.857×) — one arm at **capacity 64**, and §5.2 records **no reading step**
  for the sweep at all. Its capacity-4,096 comparator's own log carries a six-column header this
  audit's parser refuses by design, so it is named rather than guessed at.

§5.2a's four `direction` rows also stay unjudgeable, and for a reason worth stating: they are at
learning rates **1e-3 and 4e-5** against the repeats' 2e-4, and §5.2a's own result is that the
learning rate is the variable that moves rank most — so a floor measured at one rate is the last
thing that may be borrowed for another.

### 8. What this floor still does not cover

Recorded in `absent` of `P2_PROBE_FLOORS.json` rather than left as silence: step 600 and beyond;
momentum values other than 0 and 0.999; capacities other than 4,096; and the **residualised** probe
block — which is not a gap, because no number in this paper is measured on it. Three blocks still
have no floor at all (the in-run training batch, the 16-patient gate batch, the 282-patient live
checkpoint) for the same reason as before: the `d1_envelope` repeats were exported, not probed.

### 9. Suite

427 passing. Four tests that asserted the old state as protective invariants — "the two rows that
used to clear are on a block with no floor", "the only selection that clears is the published
metric", "every rank number in §5 is unjudgeable", and the like-for-like pair's `floor is None` —
now assert the new one, by name and by value. A fifth was added pinning the two properties §3 and §5
above turn on: that the floor is carried by the m = 0 arm at every step, and that it is not bimodal
anywhere while the exported-block floor is.
