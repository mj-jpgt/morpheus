## 2026-08-04 18:00 UTC — §5.2a's dissociation does not survive a floor on the *cosine*. The centring test could not be run, because the effect it was to explain is inside the co-measured statistic's own retraining spread

**Logged:** 2026-08-04 18:00 UTC. **Predeclared at commit `8635de0`**
(`NOTEBOOK_ENTRIES/PREDECLARED_centred_cosine_20260804T1700Z.md`), written and pushed **before the
workspace was built and before any arm was launched**. **How obtained:** nine GPU runs on the A100
(`150.136.45.194`) from `~/ws_j` / `~/ws_j2`, workspaces built with
`git -c core.autocrlf=false archive HEAD` and verified **file by file against `git ls-tree`** —
584 and 585 files, 0 mismatches. Scored by `v2/research/rebase/p2/p2_centred_cosine.py`, which imports
every rank statistic from `p2_envelope_floors.py` (and so from `v2/calibra/spectral.py`) and defines
the mutual cosine **once**, with a `centre` flag, so the two forms cannot drift apart. Outputs on
persistent NFS at `~/e0_run/d1_lrcentre/`.

### 0. The headline, and it is not either of the two answers the question offered

**The centring question could not be decided, and the reason is worse for the paper than either answer
would have been: the effect it was supposed to explain is not there.** At `lr = 1e-3`, on three
same-seed repeats of each of §5.2a's three arms, **the RNA-view mutual cosine's difference between the
arms is inside the cosine's own same-seed retraining spread.** The largest across-arm spread of any
repeat is **0.223**; the largest within-arm spread over three identical retrains is **0.250**. By this
paper's own criterion — the one §4.1 is built on and §4.1a applies to sixty-two rows — **no arm
difference may be read off that statistic on this block at this learning rate.**

**§5.2a's dissociation therefore dissolves, but not by the route §4.10 names.** It does not dissolve
because the difference lives in the mean-offset direction that centring removes. It dissolves because
**one seed per cell was never enough to establish it**, and the co-measured statistic the paper offers
as *more legible than rank* has a retraining floor of its own that had never been measured and is
wider than the difference read off it.

**Both instruments now agree**, which is the opposite of a dissociation: at step 200 the WSI-view R3
across the three arms spans **1.016×–1.204×** while the same three arms' own repeats span
**1.080×–1.270×**, and the cosine across the three arms spans **0.114–0.223** while their own repeats
span **0.059–0.250**. Rank says the three arms are indistinguishable. The cosine, once it is given a
floor, says the same thing.

### 1. What was run

```
d1_momentum_probe.py {0, 0.9, 0.999} 0.04 200 4096 1e-3 42 <export_dir>
```

— §5.2a's L3, L1 and L5 at their own decorrelation 0.04, capacity 4,096, seed 42 and 200-step budget,
with the harness's purely additive `export_dir` argument attached, **three same-seed repeats of each
arm**. Nine runs, one A100, concurrent, all from the same verified initialisation of canonical R1
**101.38** / R3 **67.55** — the same step-0 state every `lr_L*`, `mseed_*`, `ablate_*` and
`probefloor_*` log in this project begins at. n = 3 per arm, one seed, one stack: a **floor twice
over** in §4.1's sense, and not a distribution.

**The state/log guard passed.** The uncentred cosine recomputed from the saved states reproduces the
harness's own printed `rna-rna` column on **45 of 45 rows**, largest absolute difference
**4.97 × 10⁻⁵** against a tolerance of 10⁻⁴ — which is the rounding of the log's own `%9.4f`. The
states are the states the printed column was read off. `check_against_logs` raises rather than warns,
and the predeclaration made that a stopping condition.

### 2. The RNA-view mutual cosine, per repeat, never a mean — step 200

Uncentred, i.e. exactly the column §5.2a quotes.

| arm | rep1 | rep2 | rep3 | **within-arm spread (n = 3)** | §5.2a's own run |
|---|---:|---:|---:|---:|---:|
| m = 0 | 0.7144 | 0.9640 | 0.8723 | **0.2495** | 0.9946 |
| m = 0.9 | 0.9270 | 0.9789 | 0.9863 | 0.0593 | 0.9257 |
| m = 0.999 | 0.9292 | 0.7561 | 0.8872 | 0.1731 | **0.5207** |

**The across-arm spread, computed within each repeat index so that it is a like-for-like comparison of
three runs made under identical conditions:**

| | rep1 | rep2 | rep3 | floor (largest within-arm spread) | resolvable? |
|---|---:|---:|---:|---:|:---:|
| uncentred mutual cosine | 0.2148 | 0.2229 | 0.1140 | **0.2495** (m = 0 arm) | **no — every draw is inside** |

**§5.2a's headline movement is 0.9946 → 0.5207, a spread of 0.4739. Not one of the three repeats
reproduces it, and the arm it is read off does not reproduce it either.** With §5.2a's own run
included as a fourth same-seed run of the same configuration, the **m = 0.999 arm alone spans
0.5207–0.9292, a spread of 0.4085** — which is 86% of the across-arm movement the section reads as the
result. The m = 0 arm spans 0.7144–0.9946, 0.2802.

*(The four-run figure combines runs from two workspaces and two launches, which is why it is quoted
second and the n = 3 same-launch figure carries the verdict. The two are consistent: the n = 3 spread
is 0.2495 and adding a fourth run widens it, as adding a draw to a max/min must.)*

### 3. The centred cosine — and why this test could not decide anything

| arm | rep1 | rep2 | rep3 | within-arm spread |
|---|---:|---:|---:|---:|
| m = 0 | −0.0017 | 0.2264 | 0.0037 | 0.2281 |
| m = 0.9 | 0.0370 | 0.6362 | 0.0060 | **0.6302** |
| m = 0.999 | 0.0081 | −0.0037 | 0.0285 | 0.0322 |

Across-arm spread per repeat: **0.0388 / 0.6398 / 0.0249**, against a floor of **0.6302**. The
predeclared rule (§3 of the predeclaration) returns its **third branch — "report the magnitudes and do
not adjudicate"** — computed by `p2_centred_cosine.verdict`, not asserted here.

**Two of the five distrust conditions written into the predeclaration fired, and they are the reason
this is reported as a non-result rather than as the favourable answer.**

* **Condition 2, degeneracy.** Seven of the nine centred readings sit between **−0.004 and +0.037**.
  That is what a nearly one-dimensional centred representation must give: centred R3 is 1.02–1.29 in
  every run, so after the column mean is removed there is essentially one direction left and the mean
  off-diagonal cosine of a one-dimensional family is a function of the sign structure of a single
  coefficient. **A flat centred cosine here is partly entailed by the rank number itself**, so it
  could never have established account (B) on its own — it could only have been consistent with it.
  The predeclaration says this in those words, and it is why the entry does not bank the reading.
* **Condition 3, no floor no verdict.** The centred cosine's own within-arm spread (up to 0.630) is
  *larger* than its across-arm spread in two of three repeats. "Flat" would be unmeasured, not
  measured.

**Condition 5 is the one that decides the entry**: the uncentred movement had to be real before either
account of it could be tested, and it is not. Both (A) and (B) are moot.

### 4. The third account, which the draft does not name — and it does not fire either

`geometry()` takes two forward passes. **The rank columns are computed on `view="wsi"` and the
`rna-rna` cosine on `view="rna"`.** §5.2a and §4.10 place the 1.01× rank spread and the 1.91× cosine
movement side by side as two instruments on one block; they are two instruments on **two views**. The
predeclaration named the possibility that the difference is real and lives in the RNA view, which the
quoted rank number does not look at, and fixed the test: score rank on the RNA states too.

**It does not fire.** At step 200 the RNA-view R3 across the three arms spans **1.063×–1.247×** by
repeat, against within-arm spreads of **1.053×–1.398×**. The RNA-view rank is as flat as the WSI-view
rank, and equally unresolvable against its own floor. **The view mismatch is real and should be
corrected in the prose, but it is not carrying the dissociation.**

### 5. The named secondary: the mean-offset ratio

‖column mean‖ ⁄ RMS row norm of the centred `rna_biology` block, at step 200:

| arm | rep1 | rep2 | rep3 |
|---|---:|---:|---:|
| m = 0 | 1.586 | 5.184 | 2.620 |
| m = 0.9 | 3.572 | 6.832 | 8.504 |
| m = 0.999 | 3.630 | 1.765 | 2.811 |

Across-arm spread 2.04 / 5.07 / 5.88 against within-arm spreads of 3.60 / 4.93 / 1.87. **The offset
does vary enormously — by a factor of five between identical retrains — and it does not separate the
arms either.** It explains *why* the uncentred cosine is so unstable (the cosine is dominated by the
shared offset, and the offset is not reproducible), and it is consistent with the mean-offset account
of §4.10 as a description of *what the uncentred cosine is measuring*. It is not evidence that the
arms differ in it. It is a secondary and it does not overturn §3.

### 6. What this costs the paper, and what it buys

**It closes §6.2's last open measurement, in the direction of removing a claim.** The row read: *"the
RNA-view mutual cosine recomputed on the CENTRED representation, for §5.2a's three
high-learning-rate arms — not measured, and it is the one measurement that would settle whether
§4.10's surviving use is under strain."* It has been measured. **§4.10's surviving use is not under
strain from this observation, because the observation does not replicate.**

**Three sentences in the draft are now wrong and are flagged rather than edited** (prose is another
agent's):

1. **§5.2a**: *"Rank says the three high-rate arms are the same run; a co-measured collapse statistic
   says one of them is half as degenerate as the others."* The co-measured statistic says no such
   thing once it is given a floor: at n = 3 the arms are indistinguishable on it, and the m = 0.999
   arm's own four same-seed runs span 0.5207–0.9292. The paragraph's conclusion — *"momentum does
   nothing that rank can see"* — survives; its evidence for a **disagreement between the two
   instruments** does not.
2. **§4.10**: *"The two statistics are ordering the three runs differently at the exact reading where
   this section says rank is reliable."* They are not ordering them differently; neither of them
   orders them at all. The paragraph's two competing accounts (A) and (B) should be replaced by the
   measured fact that the cosine's difference is inside its own retraining spread.
3. **§4.10 again**, and this one is a **new cost rather than a correction**: *"if you report effective
   rank, at least also report the patient-to-patient mutual cosine and the seed spread of both."* The
   recommendation is sound and this measurement is the first time the mutual cosine has been given a
   floor on any block in this project — **and the floor is wide**: 0.25 in absolute cosine units, on
   the fixed held-out probe, at the collapse floor, at `lr = 1e-3`, between three identical same-seed
   retrains. §4.10 offers the cosine as *more legible* than rank; it should now also say that it is
   **not more reproducible**, at least here, and that the second half of its own recommendation —
   *"and the seed spread of both"* — is the load-bearing half.

**Scope, stated because the result is convenient.** Nine runs, three per arm, one seed, one stack, one
learning rate (`1e-3`, not the project's `2e-4`), one 200-step budget, one architecture, one block.
It establishes what the difference between these particular runs consists of and nothing more. In
particular it says **nothing** about the mutual cosine's reproducibility at `2e-4`, where §4.10's
collapsed-arm readings of 0.977–0.999 live and where the representations are not at the collapse floor
in the same way.

### 7. Per repeat, per step — the whole surface, because the step matters

Uncentred cosine / centred cosine, RNA view:

| arm | rep | step 50 | 100 | 150 | 200 |
|---|---|---|---|---|---|
| m = 0 | 1 | 0.644 / 0.040 | 0.026 / 0.003 | 0.111 / 0.095 | 0.714 / −0.002 |
| m = 0 | 2 | 0.932 / 0.024 | 0.959 / 0.160 | 0.918 / 0.054 | 0.964 / 0.226 |
| m = 0 | 3 | 0.790 / 0.020 | 0.935 / 0.009 | 0.866 / 0.004 | 0.872 / 0.004 |
| m = 0.9 | 1 | 0.717 / 0.209 | 0.991 / 0.062 | 0.968 / 0.054 | 0.927 / 0.037 |
| m = 0.9 | 2 | 0.430 / 0.126 | 0.984 / 0.801 | 0.978 / 0.618 | 0.979 / 0.636 |
| m = 0.9 | 3 | 0.872 / 0.103 | 0.897 / 0.059 | 0.969 / 0.003 | 0.986 / 0.006 |
| m = 0.999 | 1 | 0.435 / 0.009 | 0.995 / 0.018 | 0.987 / 0.008 | 0.929 / 0.008 |
| m = 0.999 | 2 | 0.568 / 0.059 | 0.781 / −0.004 | 0.569 / 0.003 | 0.756 / −0.004 |
| m = 0.999 | 3 | 0.885 / 0.112 | 0.955 / 0.036 | 0.960 / 0.020 | 0.887 / 0.029 |

**One run's cosine moves by more across four reading steps than the three arms move across each
other**, which is visible in §5.2a's own vendored logs before any of this ran: `lr_L1` reads
0.7057 → 0.9527 → 0.9812 → 0.9257 across steps 50–200 (a within-run range of 0.276), and the
across-arm spread at the four steps is 0.286 / 0.090 / 0.431 / 0.474. **The 0.474 the section quotes
is the largest of four, at the last step, from one seed.** That was on record and readable without a
GPU; it should have been read.

**And the rank is not flat at every step either.** At step 50 the three `lr_L*` arms read R3 5.83 /
1.07 / 1.06 — a fold of **5.500×**, not 1.01×. §5.2a's "flat at 1.01×" is a statement about step 200
specifically and is quoted as though it were about the arms.

### 8. Files

* Predeclaration: `NOTEBOOK_ENTRIES/PREDECLARED_centred_cosine_20260804T1700Z.md`, commit `8635de0`.
* Scorer: `v2/research/rebase/p2/p2_centred_cosine.py`; tests `v2/tests/test_p2_centred_cosine.py`.
* Measurement: `~/e0_run/d1_lrcentre/out/P2_CENTRED_COSINE.json` (carries the sha256 of every state
  read, every per-repeat value at every step, and the state/log guard's 45 per-row deltas), the nine
  `lrc_m*_rep*.log`, vendored to
  `v2/research/rebase/p2/figures/data/e0_run/d1_lrcentre/`.
* Harness: `v2/research/rebase/d1_momentum_probe.py`, unchanged.
