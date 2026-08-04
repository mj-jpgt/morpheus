## 2026-08-04 02:55 UTC — P2's analysis scripts are in the repository, and re-running them on a verified tree moved one table: §4.5(a)'s R2 and R3 rows were not R2 and R3

**Logged:** 2026-08-04 02:55 UTC. **How obtained:** the five scripts at `~/e0_run/p2_*.py` on
`ubuntu@150.136.45.194` vendored to `v2/research/rebase/p2/` (commit `7b37dce`), then re-run on a
freshly created workspace `~/ws_p2/morpheus` verified byte-equal to HEAD before anything was
executed. CPU only, `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=1`.
GPU untouched (a training chain was live on it at 25.7 GB throughout). Log: `~/ws_p2/out/p2_run.log`;
outputs `~/ws_p2/out/P2_METRICS_D2.json`, `P2_METRICS_D1.json`, `P2_RANK_VARIANTS.json`,
`P2_ROBUSTNESS.json`.

---

### 0. Bad news first

**`paper/P2_RANK_DRAFT.md` §4.5(a) reports two rows under the wrong statistic.** The rows labelled
**R2** and **R3** were not computed with `RANK_VARIANTS["R2"]` and `RANK_VARIANTS["R3"]`. They were
computed with the order-2 Hill number of the **eigenvalue** distribution, `(Σσ²)²/Σσ⁴`, where
`d1_audit.py`'s R2 — and therefore the canonical R2 — is the order-2 Hill number of the
**singular-value** distribution, `(Σσ)²/Σσ²`. `p2_rank_variants.py` defined both inline, under the
labels `R2` and `R3`, and both were wrong by that substitution. Three consequences, all against
the paper:

| draft statement | published | measured under the canonical statistics |
|---|---|---|
| §4.5(a) R2 row | `OK MISS`… → **2/3 D2, 2/3 D1, 4/6** | `OK OK MISS OK OK OK` → **2/3 D2, 3/3 D1, 5/6** |
| §4.5(a) R3 row | `OK OK MISS OK OK OK` | `OK MISS OK OK OK OK` — the two D2 verdicts swap |
| §4.5(a) headline: *"On 3 of the 6 pairs the three functions return different verdicts (D2 s43, D2 s44, D1 s42)"* | 3 of 6 | **2 of 6** (D2 s43, D2 s44). D1 s42 no longer disagrees. |
| §4.7.3: *"D1 scores 3/3 under the canonical R1 and under R3, but **2/3 under R2**… the only qualification that survives"* | 2/3 under R2 | **3/3 under canonical R2.** The qualification does not survive. |
| §1.4 contribution 5: *"statistic (3/6 pairs)"* | 3/6 | 2/6 |

The values that *were* published as R2 and R3 are not wrong numbers — they are correct values of a
different statistic, and `p2_competing_metrics.py` already reports that statistic under its own name.
The draft's §4.6 row **"participation ratio (raw / resid.)" is exactly §4.5(a)'s "R2" row**, cell for
cell, 4/6 — which is the arithmetic signature of the substitution, and reproduces here.

**A second, smaller error in the same section.** §4.5's provenance note explains the discrepancy
between the two source entries' R3 rows as *"one is computed on the raw block and the other on the
residualised block"*. That is not the reason. Both are on the residualised block; the difference is
the statistic. The R3 **levels** in §4.5(b) (H43 14.746 / I43 15.915; H44 6.564 / I44 6.302) are the
canonical R3 and reproduce **exactly** — it is §4.5(a)'s *verdicts* that came from `PR_rownorm`.

**A third thing the draft should know:** `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json`
**now exists** (written 2026-08-04 00:08 UTC), together with `D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json`
(02:40 UTC). §4.7.2's six `[D1 PAIRED BOOTSTRAP PENDING]` cells and §6.2's "blocked by a one-line
stale-path bug" row are out of date. Values in §4 below. I have not edited the draft — its owner should.

---

### 1. The workspace, and why a new one was made

The audit of 2026-08-03 23:59 UTC found all seven GPU workspaces drifted from HEAD, `~/ws` — where
these scripts ran — worst of the seven, with a `spectral.py` predating the rank canonicalisation
entirely. `p2_rank_variants.py` began `sys.path.insert(0, "/home/ubuntu/ws")`, so it imported from
exactly that tree.

`~/ws_p2/morpheus` was created from `git archive` of HEAD and **verified before use**, per file, by
git blob SHA-1 against `git ls-tree -r <commit>`:

```
commit 7b37dce   files 402/402   missing 0   extra 0   differ 0
```

Two notes on the verification, because both were nearly missed:

- **The first attempt failed and looked like it had succeeded.** `git archive` on Windows honours
  `core.autocrlf=true`, so every file arrived CRLF and all 402 blob hashes differed. Had I checked
  only file *count* and a couple of spot sizes, a workspace with every line ending altered would have
  passed. `git -c core.autocrlf=false -c core.eol=lf archive` is what produced the clean tree.
- The commit HEAD pointed at moved during the work (another agent pushed `c3484bb`). `c3484bb`
  differs from `7b37dce` by **one added `NOTEBOOK_ENTRIES/` file and no code**, so the workspace is
  code-identical to the current tree; the run is attributed to `7b37dce` rather than to "HEAD".

**No other agent's workspace was read from, written to or imported.**

---

### 2. What reproduced, exactly

Every value below is from `~/ws_p2/out/p2_run.log` against the draft as written.

**§4.2 — the variance decomposition, the paper's most important display item. Exact, all six cells.**

| quantity | SS_arm | SS_seed | arm share | F(3,8) | |
|---|---:|---:|---:|---:|---|
| canonical effective rank (residualised) | 1.3047 | 2.4762 | 34.5% | 1.41 | ✓ |
| RankMe as published (raw, ε = 1e-7) | 0.9772 | 2.3817 | 29.1% | 1.09 | ✓ |
| ground truth: held-out top-CCA, 40 untrained targets | 0.0353 | 0.0007 | 98.0% | 128.20 | ✓ |

Per-arm seed spreads exact: H 23.387/28.771/9.143 (3.15×), I 14.868/34.117/9.105 (3.75×),
P 29.381/24.673/11.115 (2.64×), F 13.418/7.600/6.394 (2.10×); channel folds 1.026/1.026/1.018/1.056×.
The paper's "cleanest object" — H44 against H43, 3.15× apart in rank and 0.0012 apart in channel with
the lower-rank run ahead — reproduces exactly.

*One presentational nit, not a moved number.* §4.2's prose says the channel "moves 1.8–5.6%" across
seeds. 5.6% is `fold − 1` for arm F; the script's own `rel.range` column, `(max−min)/mean`, prints
**5.4%** for the same arm. Both are correct readings of 0.5412/0.5336/0.5126; the draft should say
which it means.

**§4.4(1) — patient subsampling, 80%, 40 draws. Exact, all six rows including `gap/sd`.**
D2 s42 23.325±0.072 / 14.834±0.037, gap 8.491, 105.3 sd. D2 s43 28.657±0.115 / 33.912±0.111, 33.0 sd.
**D2 s44 9.132±0.021 / 9.093±0.019, gap 0.039, 1.39 sd** — the "hit that is not a hit". D1 150.3 /
238.2 / 145.1 sd.

**§4.5(b) — R1/R2/R3 levels on the residualised block. Exact.** H43 R2 13.227 / I43 12.972;
H43 R3 14.746 / I43 15.915; H44 R2 5.733 / I44 5.815; H44 R3 6.564 / I44 6.302.

**§4.5(c) — the view table. Exact.** Information winner H/H/H and P/P/P on all three views for all
six pairs; rank winner H/I/I, I/I/I, H/I/I on D2 and P/P/P throughout D1. Aggregate **11/18**;
D2-only **2/9**.

**§4.6 — the selection-rule table. Exact, all twelve metric rows and all six ground-truth values.**
Effective rank 5/6 (p = 0.219) raw and residualised; RankMe raw 3/3 D2 and 1/3 D1; LiDAR raw
**0/3 on D2**; stable rank 3/6; α-ReQ 4/6.

**§3.2's held-out check. Exact.** D2 s42 `wsi_biology` +0.1325 in-sample → **+0.1541** held-out, and
every arm ordering survives the held-out estimator on all three views.

**§4.7.4 — the necessity scan. Exact.** 2 violating pairs of 66; the usable one is
**H44 against I43, 3.73× lower rank (9.143 against 34.117) carrying +0.1101 more channel**,
scope "within-experiment, across arms". The other is P44 against I43, 3.07×, +0.1206, cross-experiment.

**§4.1's rank ratios.** Every D2 and D1-B ratio in §4.1's seven-row table recomputes from the values
above to the digits published (D2 s44 1.004×, s43 1.186×, s42 1.573×; D1-B s44 1.738×, s42 2.190×,
s43 3.246×).

---

### 3. What did not reproduce: §4.5(a)

`p2_rank_variants.py` now reports five statistics rather than three, because the two extra ones are
what the original was actually computing. All on the residualised held-out `wsi_biology` block:

| artifact | R1 | R2 | R3 | PR | PR_rownorm | CCA40 |
|---|---:|---:|---:|---:|---:|---:|
| H42 | 23.387 | 10.733 | 12.614 | 4.151 | 4.930 | 0.6126 |
| I42 | 14.868 | 8.187 | 9.141 | 4.101 | 4.773 | 0.4800 |
| H43 | 28.771 | 13.227 | 14.746 | 5.179 | 5.748 | 0.5970 |
| I43 | 34.117 | 12.972 | 15.915 | 4.189 | 4.825 | 0.4882 |
| H44 | 9.143 | 5.733 | 6.564 | 3.131 | 3.767 | 0.5983 |
| I44 | 9.105 | 5.815 | 6.302 | 3.543 | 3.842 | 0.4757 |
| P42 | 29.381 | 12.555 | 14.762 | 4.416 | 5.482 | 0.6117 |
| F42 | 13.418 | 8.671 | 9.118 | 4.968 | 5.286 | 0.5412 |
| P43 | 24.673 | 11.569 | 13.173 | 4.760 | 5.500 | 0.6198 |
| F43 | 7.600 | 4.858 | 5.320 | 3.239 | 3.329 | 0.5336 |
| P44 | 11.115 | 7.015 | 7.622 | 4.045 | 4.566 | 0.6087 |
| F44 | 6.394 | 4.191 | 4.726 | 2.687 | 2.928 | 0.5126 |

| statistic | D2 s42 | D2 s43 | D2 s44 | D1 s42 | D1 s43 | D1 s44 | D2 | D1 | ALL |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **R1** canonical | OK | MISS | OK | OK | OK | OK | 2/3 | 3/3 | 5/6 |
| **R2** canonical, `(Σσ)²/Σσ²` | OK | OK | MISS | OK | OK | OK | 2/3 | **3/3** | **5/6** |
| **R3** canonical, R2 on normalised rows | OK | MISS | OK | OK | OK | OK | 2/3 | 3/3 | 5/6 |
| `PR` = `(Σσ²)²/Σσ⁴` — **published as "R2"** | OK | OK | MISS | MISS | OK | OK | 2/3 | 2/3 | 4/6 |
| `PR_rownorm` — **published as "R3"** | OK | OK | MISS | OK | OK | OK | 2/3 | 3/3 | 5/6 |

The last two rows are the draft's §4.5(a) R2 and R3 rows, cell for cell.

**Which direction this cuts.** Against the paper, twice. The statistic-instability count drops from
3/6 pairs to 2/6, weakening §4.5(a) and §1.4's fifth contribution. And §4.7.3's *"partly true, and it
is the only qualification that survives"* — that D1's confirmation of necessity is 2/3 rather than
3/3 under the statistic an earlier draft nominated — is gone: under the canonical R2 D1 is a clean
3/3, so the negative result of §4.7 is **cleaner than the draft admits**, not messier. The 2/3 belongs
to `PR`, which is not one of the three statistics §3.1 enumerates and has no call site in this
repository.

**What is unaffected.** §4.5(b), which reports R1/R2/R3 *levels*, was computed by
`rank_recompute_all_instances.py` through the canonical function and reproduces exactly. Nothing in
§4.1, §4.2, §4.3, §4.4, §4.6, §4.7.2 or §4.8 depends on §4.5(a).

---

### 4. The D1 paired bootstrap has landed

`~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json`, 2026-08-04 00:08 UTC, 40 targets
(`heldout_pathway` + `immune_tme` + `tumour_state`) — i.e. the pre-restricted readout, not the
90-target one §4.7 forbids. Seed 42, `programme_free − programme_only`: point −0.0705, patient CI₉₅
**[−0.0938, −0.0444]** with `p_improve` 0.0000, cancer CI₉₅ **[−0.0957, −0.0180]** with `p_improve`
0.0045; `n_valid` 2,000 of 2,000 repeats, 84 distinct sites retained. The other two seeds and the
`D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json` negative control are in the same files.

I am not editing `paper/P2_RANK_DRAFT.md` — another agent owns it. Flagging that §4.7.2's six
`[D1 PAIRED BOOTSTRAP PENDING]` cells and §6.2's row are now fillable, and that at seed 42 the
interval excludes zero under both clusterings, so the deflation §4.7.3 records as *"Not resolved,
and marked as such"* is resolvable.

---

**§4.4(4)'s definitional check. Exact.** Canonical R1 on H42-residualised **23.3868** against a
faithful RankMe (uncentred, ε = 1e-7) **23.3909** — the same number to four significant figures, which
is the numerical answer to "you evaluated a centred variant and the published metric does better".

---

### 5. Three stale pointers found while checking every path

- `paper/P2_RANK_DRAFT.md` Appendix B cites the vendored rank-recomputation scripts as
  `v2/research/rebase/recompute_rank.py`, `recompute_p1b.py`. **Neither exists.** The files are
  `rank_recompute_all_instances.py` and `rank_recompute_phase1b.py`.
- `paper/P2_FIGURES.md` F3 cited the probe-repeat logs as `~/ws_d1/probevar_*.log`. **That path does
  not exist.** They are at `~/e0_run/d1_diag/probevar_{m0,m0.999}_{1,2,3}.log`, which is what draft
  §4.4(3) says. Corrected in the figure-plan rewrite.
- `paper/P2_RANK_DRAFT.md` §4.9 cites the surviving D1-A checkpoints as
  `~/e0_run/d1_v1/d1_{p42,p43,p44,f42}/last.pt`. **Those paths do not resolve**; the directories are
  `d1_p_seed42`, `d1_p_seed43`, `d1_p_seed44`, `d1_f_seed42`, and the four `last.pt` files are there.
  The draft's substantive statement is right — `d1_f_seed43` and `d1_f_seed44` exist as directories
  but contain no checkpoint, because the gate refused those arms.

---

### In plain terms

The scripts that produced four sections of the paper lived only on the GPU box, where nothing tested
them and nothing noticed that the copy of the code they imported was out of date. They are now in the
repository with a test that runs all five on synthetic data.

Re-running them on a tree checked file by file against the repository reproduced every published
number except one table. That table lists three functions the repository calls "effective rank" and
compares their verdicts. Two of the three had been computed with the wrong formula — a closely
related statistic with the same informal name. Correcting it makes the disagreement between the
functions smaller than the paper claims, and removes the one caveat the paper attached to the result
that went against it. Both changes are unflattering and both are stated.

Separately: the D1 confidence intervals the paper marks as pending were written to disk at
00:08 this morning.

### Files / commits

- `v2/research/rebase/p2/` — the five vendored scripts and their README (commit `7b37dce`)
- `v2/tests/test_p2_analysis_scripts.py` — 9 tests; suite 317 → 326, 53 s capped
- `v2/tests/test_effective_rank_canonical.py` — `p2_competing_metrics.py` added to `SVD_ALLOWLIST`
  with the reason (it implements RankMe, α-ReQ, LiDAR, stable rank and the eigenvalue PR, none of
  which is effective rank; effective rank itself is imported, never reimplemented)
- Box: `~/ws_p2/morpheus` (verified tree), `~/ws_p2/out/` (log and JSON)
- Prior: `WORKSPACE_DRIFT_AUDIT_ALL_20260803T2359Z.md`,
  `p2_competing_metrics_and_necessity_test_20260803T2326Z.md`,
  `effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`
