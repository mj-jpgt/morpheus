# P2: competing label-free metrics as selection rules, and the necessity test

UTC 2026-08-03T23:26Z. Branch `research/rebase-vision`.

Scope: Tasks 2 and 3 of the P2 assignment — compute the published competing metrics on our own
frozen artifacts and score them **as selection rules**, and run the necessity test that RankMe's
"necessary but not sufficient" hedge is vulnerable to. Task 1 (prior-art citation-graph sweep) is
logged separately.

**This entry does not edit `paper/P2_RANK_DRAFT.md`, `NOTEBOOK.md` or `v2/`.** Another agent owns the
draft and a third unified `effective_rank` (commit `85c0fa8`); I touched none of their files. Their
canonical definition — Roy & Vetterli 2007 Definition 1, order 1, column-centred — is exactly the one
used for every "effective rank" number below, and §4.3 uses their R1/R2/R3 taxonomy.

---

## 0. Bad news first

Four results here cut against the paper as currently framed. They are stated before the supporting
material so they cannot be buried.

1. **D1 does not violate necessity — it confirms rank.** The necessity test was set up expecting
   `programme_free` (low rank) to match or beat `programme_only` (high rank) on the molecular
   channel. It does not. `programme_only` wins the channel on **3/3 seeds**, in the same direction
   as its rank advantage. Under §4.8.4 of the draft this is the *first* row of the preregistered
   outcome table — "large rank gap, large channel gap in the same direction" — which that table
   says must be **reported as a negative for this paper's generality, in §4.1 and the abstract, at
   the same prominence as the instances that dissociate.** Qualified by §4.3: this is 3/3 under the
   canonical R1 and under R3, but **2/3 under R2**, the statistic draft §4.8.3 nominates for the D1
   rank column.
2. **D1's result also trips a preregistered escalation elsewhere.** `d1_v2/D1_PAIR_MANIFEST.json`
   records `"preregistered_prediction": "programme_free >= programme_only on the held-out molecular
   channel; if programme_only wins, the collapse story is wrong -- escalate, do not proceed to D2"`.
   `programme_only` wins 3/3. This is flagged, not resolved, here.
3. **On the canonical readout view, plain effective rank is the *best* of the metrics tested, not
   the worst.** Over the 6 matched pairs it scores 5/6; RankMe-as-published scores 4/6, LiDAR 3–4/6,
   α-ReQ 4/6. "Why did you not just use the better metric" currently has the answer "there was not
   one", which is useful, but it is not the answer the framing anticipated.
4. **No metric, including ours, is statistically distinguishable from a coin flip on 6 pairs.**
   Best exact two-sided binomial p = 0.219. The selection-rule experiment as designed is
   underpowered and must not be quoted as if it settled anything.

The claim that **does** survive, and survives in a stronger and better-powered form than the
selection-rule count, is in §4 and §5: effective rank's variation across these artifacts is
**two-thirds training-seed nuisance**, against a ground-truth channel that is **98% arm effect**;
and the rank verdict on a pair **flips depending on which co-trained view of the same model you
measure it on**, while the information verdict does not.

---

## 1. What was run

Scripts (on the box, `~/e0_run/`, all single-threaded, no GPU):

| script | purpose | output |
|---|---|---|
| `p2_competing_metrics.py` | computes every metric + ground truth per artifact | `P2_METRICS_D2.json`, `P2_METRICS_D1.json`, `P2_METRICS_ALL_SUBSAMPLED.json` |
| `p2_selection_rule.py` | scores each metric as a selection rule; noise floor | stdout (reproduced below) |
| `p2_necessity_and_variance.py` | arm-vs-seed variance decomposition; necessity scan | stdout |
| `p2_robustness.py` | held-out CCA; all three representation views | `P2_ROBUSTNESS.json` |

`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=1` throughout; GPU
untouched (a training chain is live on it). Whole sweep runs in ~10 s per 6 artifacts.

Artifacts, all frozen:

- D2 — `~/e0_run/d2_v3/d2_v3_s{42,43,44}/artifacts/d2_{h,i}_seed{42,43,44}.npz`. `H` = Hallmark,
  `I` = PBS.
- D1 — `~/e0_run/d1_v2/artifacts/d1_{f,p}_seed{42,43,44}.npz`. Per `d1_v2/D1_PAIR_MANIFEST.json`,
  `F` = `programme_free`, `P` = `programme_only`; `"objective_only_difference": true`.
  **All six D1 runs are complete and live**: `~/e0_run/d1_audit.log` records
  `[PASS] A1_all_six_runs_complete` with `overfit_present: True` and `all_grads_positive: True` for
  all six. The audit's subsequent readout step failed only on a stale absolute path
  (`/lambda/nfs/.../d1_v2/artifacts/...` vs the real `/home/ubuntu/e0_run/...`), not on the data.
  That is a **one-line path bug in the D1 audit chain and it should be fixed** — as of this entry
  the official `D1_PAIRED_BOOTSTRAP_STRATIFIED.json` still does not exist.

Ground truth = top canonical correlation (16 components) between the confound-residualised
representation and the residualised **40 targets neither arm trained on**
(`heldout_pathway` + `immune_tme` + `tumour_state`), reproducing `d2_readout.py`'s residualisation,
seed and component budget exactly. Verified against `D2_PER_ARTIFACT_READOUT.json`: my D2 numbers
reproduce it to 4 d.p., and the recovered Δ values (−0.1325 / −0.1089 / −0.1226) match the
assignment's stated ground truth exactly.

---

## 2. The metrics, and what is and is not faithful

| metric | source | fidelity |
|---|---|---|
| effective rank | Roy & Vetterli, EUSIPCO 2007; this repo's `v2/calibra/spectral.py` | canonical; column-centred, exp of entropy of L1-normalised singular values |
| **RankMe** | Garrido, Balestriero, Najman & LeCun, ICML 2023, arXiv:2210.02885 | **faithful.** `p_k = σ_k(Z)/‖σ(Z)‖₁ + ε`, ε = 1e-7, on the **uncentred** embedding matrix, no zero-pruning — per the verbatim definition already verified at full text in draft §2.1 |
| participation ratio | standard second-moment dimensionality statistic | not attributed to any paper; included so the claim does not rest on one functional form |
| stable rank | ‖X‖²_F/‖X‖²₂ | as above |
| **α-ReQ** | Agrawal, Mondal, Ghosh & Richards, NeurIPS 2022 | **faithful to the authors' released code** — see §6.2 |
| **LiDAR** | Thilak et al., **ICLR 2024** | **adapted (licensed by the paper's own footnote 4); δ swept** — see §6.1 |

### Reference corrections the draft needs (all verified this pass — provenance in §6.4)

1. **LiDAR has a peer-reviewed venue: ICLR 2024.** Draft §2.2 says `[UNVERIFIED: peer-reviewed
   venue … Cite as an arXiv preprint.]` That is now wrong. DBLP record
   `https://dblp.org/rec/conf/iclr/Thilak0SDGNSL24` gives *The Twelfth International Conference on
   Learning Representations, ICLR 2024, Vienna, Austria, May 7–11, 2024*, OpenReview.net,
   `https://openreview.net/forum?id=f3g5XpL9Kb`. (The acceptance *track* could not be retrieved —
   OpenReview returns a bot challenge — so do not state poster/spotlight/oral.)
2. **α-ReQ's NeurIPS 2022 venue is confirmed; draft §2.1's `[COULD-NOT-VERIFY]` can be cleared.**
   Four independent retrieved records agree: DBLP `conf/nips/AgrawalMGR22`; the NeurIPS proceedings
   page and its official BibTeX (vol. 35, **pp. 17626–17638**, `doi:10.52202/068431-1281`);
   OpenAlex W7133259885; and the camera-ready PDF footer "36th Conference on Neural Information
   Processing Systems (NeurIPS 2022)". Two title caveats: the PDF's own title page omits "in
   Self-Supervised Learning" (present only in the proceedings metadata), and the NeurIPS BibTeX
   renders the first author "Agrawal, Kumar K" where DBLP has "Kumar Krishna Agrawal".
3. **Roy & Vetterli full citation**, for §2.1: Olivier Roy & Martin Vetterli, "The effective rank: A
   measure of effective dimensionality", *15th European Signal Processing Conference (EUSIPCO
   2007)*, Poznan, Poland, 3–7 September 2007, IEEE, **pp. 606–610**; IEEE Xplore document
   **7098875**; DOI **`10.5281/zenodo.40328`** (a Zenodo/EURASIP archive DOI — **no IEEE
   `10.1109/…` DOI was found**, so do not invent one). The original has **no ε** and normalises
   **singular values**; the ε in both RankMe and LiDAR is a later addition, not Roy & Vetterli's.

### Two verbatim admissions in LiDAR that the paper should use

Both from arXiv:2312.04000v1 and both help us:

- §5: *"with VICReg, LiDAR achieves optimal results when applied to the representation rather than
  the embedding, as embedding-based evaluations result in dramatic performance degradation, a
  phenomenon aligning with the non-monotonic relationship between rank and performance reported by
  (Garrido et al., 2022). **This illustrates that high rank is a necessary but not a sufficient
  condition for high performance.**"*
- Limitations: *"we have observed instances where **the LiDAR metric exhibits a negative
  correlation with probe accuracy**, particularly pronounced in scenarios like VICReg when dealing
  with higher dimensional embeddings."*

The second is directly on point for our D2 result: LiDAR's own authors report anti-correlation with
downstream performance in some regimes, so our 0/3 on D2 is a documented failure mode of the metric
rather than an artifact of our adaptation alone.

### A finding about RankMe that the paper needs

`rankme_residualised` = 23.391 where `effective_rank_residualised` = 23.387 (H42), and equally close
on all 12 artifacts. **The ε and the RankMe-vs-Roy-Vetterli distinction are numerically irrelevant.**
Everything that separates RankMe's score from ours is the single preprocessing choice of **column
centring**: RankMe as published is computed on the uncentred `Z`, so its leading singular value
absorbs the mean offset of the embedding cloud.

This matters because it is a referee's best line of attack — "you evaluated a centred variant; the
published metric does better" — and it is true on D2: `rankme_raw` scores 3/3 on D2 where centred
effective rank scores 2/3. It is equally true that this reverses on D1, where `rankme_raw` scores
1/3 against centred rank's 3/3. **RankMe's advantage on D2 comes from the mean-offset direction, not
from rank**, and that mechanism is not stable across experiments. Both halves must be reported.

---

## 3. Selection-rule scores

Ground truth (top-CCA, 40 untrained targets, residualised, `wsi_biology`):

| pair | arm A | arm B | A | B | Δ(A−B) | winner |
|---|---|---|---:|---:|---:|:---:|
| D2 s42 | H42 | I42 | 0.6126 | 0.4800 | +0.1325 | H |
| D2 s43 | H43 | I43 | 0.5970 | 0.4882 | +0.1089 | H |
| D2 s44 | H44 | I44 | 0.5983 | 0.4757 | +0.1226 | H |
| D1 s42 | P42 | F42 | 0.6117 | 0.5412 | +0.0705 | P |
| D1 s43 | P43 | F43 | 0.6198 | 0.5336 | +0.0863 | P |
| D1 s44 | P44 | F44 | 0.6087 | 0.5126 | +0.0961 | P |

Does the metric pick the arm that actually carries more molecular information?

| metric | D2 s42 | D2 s43 | D2 s44 | D1 s42 | D1 s43 | D1 s44 | D2 | D1 | ALL | p (2-sided binomial) |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|---:|
| effective rank (raw) | OK | MISS | OK | OK | OK | OK | 2/3 | 3/3 | **5/6** | 0.219 |
| effective rank (residualised) | OK | MISS | OK | OK | OK | OK | 2/3 | 3/3 | **5/6** | 0.219 |
| **RankMe (raw, as published)** | OK | OK | OK | MISS | OK | MISS | **3/3** | 1/3 | 4/6 | 0.688 |
| RankMe (residualised) | OK | MISS | OK | OK | OK | OK | 2/3 | 3/3 | 5/6 | 0.219 |
| participation ratio (raw) | OK | OK | MISS | MISS | OK | OK | 2/3 | 2/3 | 4/6 | 0.688 |
| participation ratio (resid.) | OK | OK | MISS | MISS | OK | OK | 2/3 | 2/3 | 4/6 | 0.688 |
| stable rank (raw) | OK | OK | MISS | MISS | OK | MISS | 2/3 | 1/3 | 3/6 | 1.000 |
| stable rank (resid.) | MISS | OK | MISS | MISS | OK | OK | 1/3 | 2/3 | 3/6 | 1.000 |
| α-ReQ \|α−1\| (raw) | OK | MISS | OK | OK | OK | MISS | 2/3 | 2/3 | 4/6 | 0.688 |
| α-ReQ \|α−1\| (resid.) | OK | MISS | OK | OK | OK | MISS | 2/3 | 2/3 | 4/6 | 0.688 |
| **LiDAR (raw)** | MISS | MISS | MISS | OK | OK | OK | **0/3** | 3/3 | 3/6 | 1.000 |
| LiDAR (residualised) | MISS | OK | MISS | OK | OK | OK | 1/3 | 3/3 | 4/6 | 0.688 |

**Answer to "why did you not use the better metric": on our data there is not one.** LiDAR — the
strongest published alternative, and the one the draft flags as most constraining — scores **0/3 on
D2**, the in-scope matched-arm contrast, choosing the information-poor arm every time, and it does so
**at every value of its unspecified δ across eight orders of magnitude** (§6.1). RankMe wins D2 and
loses D1. α-ReQ's 4/6 is squeezed out of α differences its own authors would call meaningless, since
all 12 artifacts sit far outside its "Goldilocks zone" (§6.2). Nothing here is significant at n = 6,
and that limitation belongs in the paper.

### The design cannot establish that a selection rule works — only that it is perfect

Exact two-sided binomial at n = 6: a **perfect 6/6 gives p = 0.031**, 5/6 gives 0.219, 4/6 gives
0.688. So the matched-pair design (2 experiments × 3 seeds) has just enough power to detect a
*flawless* rule and none at all to detect a merely good one. **No conclusion of the form "metric X is
better than metric Y" can be drawn from these six pairs**, in either direction, and the paper must
say so rather than letting 5/6 vs 4/6 do rhetorical work. This is the main reason §4's variance
decomposition — which uses magnitudes and estimates the nuisance term from 8 within-arm degrees of
freedom — should carry the argument instead.

### The D2 s44 "hit" is not a hit

Patient-subsampling (80%, 40 draws) gives the measurement noise of each metric **given a fixed
trained model**:

| pair | eff. rank A | eff. rank B | gap | gap/sd |
|---|---:|---:|---:|---:|
| D2 s42 | 23.325 ± 0.072 | 14.834 ± 0.037 | 8.491 | 105.3 |
| D2 s43 | 28.657 ± 0.115 | 33.912 ± 0.111 | 5.255 | 33.0 |
| **D2 s44** | **9.132 ± 0.021** | **9.093 ± 0.019** | **0.039** | **1.39** |
| D1 s42 | 29.244 ± 0.102 | 13.393 ± 0.027 | 15.850 | 150.3 |
| D1 s43 | 24.584 ± 0.069 | 7.586 ± 0.016 | 16.997 | 238.2 |
| D1 s44 | 11.098 ± 0.031 | 6.387 ± 0.011 | 4.710 | 145.1 |

At D2 s44 the between-arm rank gap is **1.4 sampling standard deviations** — an unresolvable tie,
against a ground-truth gap of +0.1226 that the paired bootstrap resolves with a CI excluding zero.
Effective rank's honest D2 record is therefore **1 clear hit, 1 clear miss, 1 pair it cannot
resolve at all**, not 2/3.

Note the direction of this: **measurement noise is negligible** (SD ≈ 0.1 on a rank of 25; the
metric is a *precise* measurement). The instability is entirely in **training**. That distinction
should be made explicitly in the paper, because §4.7's "reproducibility floor" can otherwise be
misread as an estimator problem.

---

## 4. The result that survives, and is better powered than any of the above

### 4.1 Variance decomposition: arm (signal) vs training seed (nuisance)

12 artifacts = 4 arms (D2-H, D2-I, D1-P, D1-F) × 3 seeds. The seed changes **nothing** about
objective, architecture, data, split or schedule. Rank-type metrics decomposed on the log scale
(they are multiplicative and span 6.4–34.1); the CCA ground truth on the raw scale.

| quantity | SS_arm | SS_seed | **arm share** | F(3,8) |
|---|---:|---:|---:|---:|
| effective rank (residualised) | 1.3047 | 2.4762 | **34.5%** | 1.41 |
| RankMe (raw, as published) | 0.9772 | 2.3817 | **29.1%** | 1.09 |
| **ground truth: untrained40 top-CCA** | 0.0353 | 0.0007 | **98.0%** | **128.20** |

**Two-thirds of the variation in effective rank across these artifacts is training-seed nuisance.
Two percent of the variation in the information channel is.** The arm effect on rank is not even
significant (F = 1.41 on 3,8 df); the arm effect on information is overwhelming (F = 128).

This is the paper's strongest statement and it does not depend on a 6-way sign count. It says the
proxy's dynamic range is dominated by exactly the factor that carries no information.

Per-arm, holding everything but the seed fixed:

| arm | effective rank across 3 seeds | fold | untrained40 CCA across 3 seeds | fold |
|---|---|---:|---|---:|
| D2 H (Hallmark) | 23.387 / 28.771 / 9.143 | **3.15×** | 0.6126 / 0.5970 / 0.5983 | 1.026× |
| D2 I (PBS) | 14.868 / 34.117 / 9.105 | **3.75×** | 0.4800 / 0.4882 / 0.4757 | 1.026× |
| D1 P (`programme_only`) | 29.381 / 24.673 / 11.115 | **2.64×** | 0.6117 / 0.6198 / 0.6087 | 1.018× |
| D1 F (`programme_free`) | 13.418 / 7.600 / 6.394 | **2.10×** | 0.5412 / 0.5336 / 0.5126 | 1.056× |

**This is inside RankMe's own reserved scope.** RankMe states verbatim that it *"should however only
be used to compare different runs of a given method"*. Three seeds of one arm **are** different runs
of a given method. Across them the metric moves 2.1–3.75× while the quantity it is a proxy for moves
1.8–5.6%. The between-arm gaps it is being asked to resolve are 12–24% relative — i.e. **the
nuisance-induced range of the proxy exceeds the signal it must resolve by roughly an order of
magnitude, in the one regime its authors reserve for it.**

Sharpest single instance, strictly within one arm: **H44 vs H43 — 3.15× apart in effective rank
(9.143 vs 28.771), 0.0012 apart in the molecular channel (0.5983 vs 0.5970, the lower-rank run
marginally ahead).** Same objective, same data, same schedule, same architecture; only the seed
differs.

### 4.2 The rank verdict is view-dependent; the information verdict is not

Each artifact stores three co-trained views of the same model: `wsi_biology` (the canonical readout
view), `rna_biology`, `full_biology`. Recomputing both quantities on each:

| pair | info winner (wsi / rna / full) | rank winner (wsi / rna / full) | rank verdict stable? |
|---|---|---|:---:|
| D2 s42 | H / H / H | **H / I / I** | **NO** |
| D2 s43 | H / H / H | I / I / I | yes (wrong) |
| D2 s44 | H / H / H | **H / I / I** | **NO** |
| D1 s42 | P / P / P | P / P / P | yes |
| D1 s43 | P / P / P | P / P / P | yes |
| D1 s44 | P / P / P | P / P / P | yes |

**The information ordering is identical under all three views for all six pairs. The rank ordering
is not: two of six pairs give opposite verdicts depending on which head of the same model you
measure.** Aggregated over all 18 (pair × view) comparisons effective rank is right 11/18
(p ≈ 0.48); restricted to D2, the in-scope contrast, it is right **2/9** — worse than chance, in a
contrast where the information ordering is 9/9 consistent.

*Caveat, stated plainly:* `rna_biology`→RNA-derived pathway targets is partly circular and its
absolute CCA (0.79–0.85) should not be read as a clean image→molecular channel. That is why the
canonical readout is `wsi_biology`. The **rank** measurements on the other views are not affected by
this circularity, and the instability of the rank *verdict* across views stands on its own.

### 4.3 The three statistics this repository calls "effective rank" disagree as selection rules

Commit `85c0fa8` ("One effective_rank, ten call sites") establishes that three mutually incompatible
statistics were in use here under one name: **R1** = Roy–Vetterli, order 1, column-centred (canonical,
`spectral.py`); **R2** = order-2 participation ratio (`d1_audit.py`); **R3** = R2 on L2-normalised rows
(`d1_geometry_probe.py`, `training.py`, `runner.py`). All three are computed on the same 12 artifacts,
same residualisation, same ground truth:

| statistic | D2 s42 | D2 s43 | D2 s44 | D1 s42 | D1 s43 | D1 s44 | D2 | D1 | ALL |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **R1** canonical (Roy–Vetterli) | OK | **MISS** | **OK** | **OK** | OK | OK | 2/3 | 3/3 | 5/6 |
| **R2** `d1_audit.py` | OK | **OK** | **MISS** | **MISS** | OK | OK | 2/3 | 2/3 | 4/6 |
| **R3** `d1_geometry_probe.py` | OK | **OK** | **MISS** | **OK** | OK | OK | 2/3 | 3/3 | 5/6 |

**On 3 of the 6 pairs the three statistics return different verdicts** (D2 s43, D2 s44, D1 s42). So
"does effective rank track information content" does not have one answer even holding the data,
the arms and the ground truth fixed — it depends on which of three functions, all of them called
effective rank in this repository and all of them in use in the published literature, is chosen.

**This directly qualifies §0's bad news item 1.** The draft's §4.8.3 states that D1's rank column
"will be **statistic R2**". **Under R2, D1 scores 2/3, not 3/3.** The finding that "D1 confirms
necessity 3/3" holds under R1 and R3 and *not* under the statistic the draft nominated for that
table. The D1 rows must state which statistic they use, and the honest summary of D1 is "rank
ordering agrees with the information ordering on 2–3 of 3 seeds depending on the rank statistic",
not a clean 3/3.

### 4.4 Held-out ground truth

A referee will note the headline CCA is in-sample at a 16-component budget and therefore
upward-biased. Re-run with `heldout_top_cca` (directions fit on half the patients, scored on the
other half): **every arm ordering survives on all three views, with equal or larger deltas** (e.g.
D2 s42 wsi: +0.1325 in-sample → +0.1541 held-out). The ground truth is not an in-sample artifact.

---

## 5. The necessity test

### 5.1 What was pre-declared as a violation

RankMe's defence is that high rank is *"a necessary (but not sufficient) condition for good
downstream performance"*. Under that hedge **high rank + low information is predicted and is not a
counterexample**; only **low rank + high information** breaks it.

Criterion, fixed before inspecting the pair list, in `p2_necessity_and_variance.py`:

> a pair (lo, hi) is a violation iff `eff_rank(hi)/eff_rank(lo) ≥ 2.0` **and**
> `CCA(lo) − CCA(hi) ≥ 0.0705`.

Both thresholds come from quantities established independently of this analysis. 2.0× is *below* the
2.10–3.75× that the seed alone produces within a single arm, so it is a conservative floor. 0.0705 is
the **smallest** between-arm channel gap this project has accepted as real (D1 seed 42).

### 5.2 Result: the D1 candidate FAILS to violate necessity

The assignment's candidate was D1 `programme_free` (low rank) vs `programme_only` (high rank).
Measured on `d1_v2`:

| seed | rank, `programme_only` | rank, `programme_free` | ratio | channel, `programme_only` | channel, `programme_free` | Δ | rank ordering correct? |
|---|---:|---:|---:|---:|---:|---:|:---:|
| 42 | 29.381 | 13.418 | 2.19× | 0.6117 | 0.5412 | +0.0705 | **yes** |
| 43 | 24.673 | 7.600 | 3.25× | 0.6198 | 0.5336 | +0.0863 | **yes** |
| 44 | 11.115 | 6.394 | 1.74× | 0.6087 | 0.5126 | +0.0961 | **yes** |

`programme_free` has both the lower rank **and** the lower channel, 3/3. **Necessity is not
violated. It is confirmed.** This is the draft's §4.8.4 row 1 and must be reported at abstract
prominence per that table's own instruction.

Two notes on interpretation:

- `programme_free` is **not collapsed** this time — ranks 6.4–13.4, not the ~1.7 of D1-A, and its
  channel (0.513–0.541) sits clearly above its own `random_control` (0.443–0.474). So this is the
  "both arms train" branch, not the collapse branch. The queue fix worked.
- The rank ratios (1.74–3.25×) are **not larger than the within-arm seed range** (2.10–3.75×). So
  even in the instance where rank gets the answer right, the size of the gap it is reading is
  inside its own nuisance band.
- The assignment's stated figures (~12 vs ~111, ≈9×) do not correspond to `d1_v2`. The ~111 figure
  is not reproduced by any statistic here; the earlier `d1_v1` numbers in draft §4.8.2 are statistic
  **R3** (9.81 / 10.47 / 1.71), a different statistic on a different, collapsed run. The two must
  not be compared — §4.8.3 already says this.

### 5.3 A violation does exist, but not from D1's contrast

Scanning all 66 ordered pairs of the 12 artifacts against the pre-declared criterion, **2 pairs
violate**:

| lower-rank artifact | rank | CCA | higher-rank artifact | rank | CCA | fold | ΔCCA | scope |
|---|---:|---:|---|---:|---:|---:|---:|---|
| H44 | 9.143 | 0.5983 | I43 | 34.117 | 0.4882 | **3.73×** | **+0.1101** | within-experiment (D2), across arms |
| P44 | 11.115 | 0.6087 | I43 | 34.117 | 0.4882 | 3.07× | +0.1206 | cross-experiment |

The first is the usable one: **H44 sits at 3.73× lower effective rank than I43 and carries +0.110
more molecular channel** — a gap comparable to the headline D2 arm effect. Both are the same
architecture, cohort, schedule and modality pair; they differ in arm and seed.

**How hard to push this.** It is a genuine low-rank/high-information instance and it is the only
configuration RankMe's hedge cannot absorb. But it is a cross-arm, cross-seed comparison, and RankMe
reserves itself to *"different runs of a given method"* — a referee will argue H44 and I43 are not
that. My recommendation is to lead with §4.1's variance decomposition and the H44-vs-H43 within-arm
instance (strictly in scope, kills usefulness, does not formally break necessity), and to present
the H44/I43 pair as a **supporting** necessity counterexample with its scope stated, rather than as
the load-bearing one. Claiming it as an in-scope necessity violation is the one place in this
analysis where the paper could be caught overreaching.

**And it is partially pre-empted.** The companion prior-art sweep
(`NOTEBOOK_ENTRIES/p2_prior_art_citation_graph_sweep_20260803T2326Z.md`, §2 A1) found that Aldeneh,
Thilak, Higuchi, Theobald & Likhomanenko, *"Towards Automatic Assessment of Self-Supervised Speech
Models using Rank"*, **ICASSP 2025** (DOI `10.1109/ICASSP49660.2025.10889651`, arXiv:2409.10787) already
reports, verbatim, that *"rank does not reliably predict the best-performing layer for specific
downstream tasks, **as lower-ranked layers can outperform higher-ranked ones**"* — a published
low-rank/high-information instance, in a peer-reviewed venue, co-authored by LiDAR's first author. It
varies **layer depth within one trained encoder** rather than comparing separately trained runs, so
it is not the same regime, but **this section must be written as corroborating that result in a new
regime, not as discovering the configuration.** A version of §5.3 that presents low-rank/high-
information as novel is not submittable.

### 5.4 Standing instruction for future artifacts

The comparison is set up and re-runnable in ~10 s per 6 artifacts. For any new arm pair, a
**violation of necessity** is declared iff, on the canonical `wsi_biology` view with the same
residualisation and the 40 untrained targets:

`eff_rank(hi)/eff_rank(lo) ≥ 2.0` **and** `CCA(lo) − CCA(hi) ≥ 0.0705`, with the CCA gap's paired
bootstrap CI₉₅ excluding zero, **and** the ordering surviving `heldout_top_cca`.

Anything weaker than that is a dissociation of *usefulness*, not of *necessity*, and must be labelled
as such.

---

## 6. Fidelity — what is now solid, and what remains a debt

The source definitions were retrieved and the two provisional implementations were corrected. The
selection-rule verdicts did **not** change under either correction, which is itself worth stating.

### 6.1 LiDAR — definition confirmed, adaptation licensed, δ swept

Definition retrieved verbatim from arXiv:2312.04000v1 §3, Eqs. (1)–(4). My implementation matches
it exactly: `Σ_b = E[(µ_x − µ)(µ_x − µ)ᵀ]`, `Σ_w = E E[(e(x̃) − µ_x)(e(x̃) − µ_x)ᵀ] + δI_p`,
`Σ_lidar = Σ_w^{−1/2} Σ_b Σ_w^{−1/2}`, then the Roy–Vetterli smooth rank on its eigenvalues with
`p_i = λ_i/‖λ‖₁ + ε`.

- **The modality-as-view adaptation is explicitly licensed by the paper.** Footnote 4 defines the
  within-class samples as *"Data augmentations, **or otherwise data points which are treated as
  positive samples**"*. Our objective treats a patient's WSI-derived and RNA-derived embeddings as
  exactly that positive pair. This is stronger footing than the previous draft of this entry
  assumed and the adaptation can be defended, not merely disclosed.
- **Their `n > p` recommendation is met**: n = 2,766 patients > p = 256 dimensions. (Verbatim:
  *"we recommend choosing a value for n that is greater than the length of feature vectors (p)"*.)
- **The remaining, genuine gap is q.** They use q = 50 (recommended) or q = 10, sweep q only over
  10–100, and **never discuss q = 2 or state a minimum q**. We have q = 2, the bare minimum at which
  Σ_w is estimable. This must be disclosed; it is the one respect in which our LiDAR is outside the
  authors' tested regime.
- **δ and ε are given NO numeric value anywhere in the paper** — only Proposition 1's symbolic
  constraints (`ε < 1 − ‖λ‖_∞/‖λ‖₁`, `δ < (e⁻¹ − ε)‖λ‖₁`) and the remark "for δ ≪ ‖λ‖₁ …". Since
  the constant is unspecified, I **swept it over eight orders of magnitude**:

| δ | LiDAR (raw): D2 / D1 / all | LiDAR (residualised): D2 / D1 / all |
|---|---|---|
| 1e-8 | 0/3 · 3/3 · 3/6 | 0/3 · 3/3 · 3/6 |
| 1e-6 | 0/3 · 3/3 · 3/6 | 0/3 · 3/3 · 3/6 |
| 1e-4 | 0/3 · 3/3 · 3/6 | 1/3 · 3/3 · 4/6 |
| 1e-3 | 0/3 · 3/3 · 3/6 | 1/3 · 3/3 · 4/6 |
| 1e-2 | 0/3 · 3/3 · 3/6 | 1/3 · 3/3 · 4/6 |
| 1e-1 | 0/3 · 3/3 · 3/6 | 0/3 · 3/3 · 3/6 |
| 1e+0 | 0/3 · 3/3 · 3/6 | 0/3 · 3/3 · 3/6 |

  **LiDAR scores 0/3 on D2 at every δ tested, and 3/3 on D1 at every δ tested.** The unspecified
  constant does not rescue it and cannot be blamed for the result. The absolute LiDAR value moves a
  lot with δ (H42 raw: 176.1 at δ=1e-8 down to 7.0 at δ=1) but the *ordering* — which is all a
  selection rule uses — is invariant.

### 6.2 α-ReQ — re-implemented from the authors' code, verdict unchanged

The **paper text is not sufficient to reproduce a number.** It says only *"we compute the full set of
numerical eigenvalues, and estimate α by fitting a linear model on the eigenspectrum in log-log
scale"* and gives **no index range and no weighting**. The index range exists only in the authors'
released `fastssl` repository (`fastssl/utils/powerlaw.py::stringer_get_powerlaw`, invoked from
`scripts/train_model.py` as `trange=np.arange(3, 100)`): centre → PCA → `explained_variance_ratio_`
→ **weighted** least squares of log λ on −log(rank) over **eigenvalue ranks 4–100**, weights
**w = 1/rank**. I replaced my provisional unweighted 11–50 OLS fit with this procedure. **The
selection-rule verdict is byte-identical (4/6, same per-pair OK/MISS pattern)**, so the earlier
result was not an artifact of the fit range.

**Flag for the draft:** the paper says "full set", the code fits ranks ~4–100. The index range must
**not** be attributed to the paper text.

Measured α (authors' procedure, residualised): H 3.14 / 3.06 / 4.18, I 3.75 / 2.65 / 4.34,
P 2.97 / 3.28 / 4.38, F 4.66 / 3.96 / 4.16 (fit R² 0.87–0.97).

**A caveat on how I scored it.** α-ReQ does **not** state a `|α − 1|` rule anywhere. It states a
*"Goldilocks zone"*: *"the best representations are those where α is in a range that is close to 1"*
and *"either too high or too low an α value implies poor generalization"*; its Algorithm 1 uses an
adaptively-updated interval, initialised to `[0, 1]`. `|α − 1|` is **our** operationalisation and
must be labelled as such. Note also that **every one of our 12 artifacts has α between 2.6 and 4.8**,
i.e. all far outside the zone α-ReQ calls good — so on α-ReQ's own terms it declines to distinguish
them, and the 4/6 above is squeezed out of differences the metric's authors would call meaningless.

Note further that α-ReQ makes **the same hedge as RankMe**, verbatim: *"a task-agnostic measure like
α is a **necessary but not sufficient condition**"*. §5's necessity framing therefore applies to
α-ReQ too, and the paper can say so.

### 6.3 Remaining debts

1. **RankMe is faithful** and needs no further verification — the definition was verified at full
   text in draft §2.1, and my implementation reproduces the canonical `effective_rank` to 4 d.p. on
   a centred input, an independent check on both.
2. The D1 numbers here are **point estimates without a paired bootstrap**. The official D1 readout
   has not run (path bug, §1). Direction is consistent 3/3 and survives held-out estimation, but the
   D1 rows should carry CIs before publication.
3. LiDAR's **ICLR 2024 camera-ready could not be retrieved** (OpenReview bot challenge); everything
   quoted is from arXiv v1. If the camera-ready changed constants, that is unchecked.
4. The LiDAR paper contains an **internal inconsistency** in its own hyperparameters (main text says
   n = 10k for VICReg; Appendix 11.1 says n = 5000 for SimCLR *and* VICReg; the VICReg appendix says
   10,000 images). Irrelevant to our numbers, but do not cite a single n for VICReg as if the paper
   were self-consistent.

### 6.4 Retrieval provenance for §2's reference corrections

- LiDAR venue: DBLP API `https://dblp.org/search/publ/api?q=LiDAR+Sensing+Linear+Probing+Performance&format=json`
  → BibTeX `https://dblp.org/rec/conf/iclr/Thilak0SDGNSL24.bib?param=1` (booktitle, publisher, year,
  OpenReview URL verified from this record). arXiv metadata (authors, "Technical report" comment,
  v1-only submission history) from `http://export.arxiv.org/api/query?id_list=2312.04000` and
  `https://arxiv.org/abs/2312.04000`. Formulas from `https://arxiv.org/pdf/2312.04000v1`.
- α-ReQ venue: DBLP `https://dblp.org/rec/conf/nips/AgrawalMGR22`; NeurIPS proceedings page
  `http://papers.nips.cc/paper_files/paper/2022/hash/70596d70542c51c8d9b4e423f4bf2736-Abstract-Conference.html`;
  official BibTeX `https://papers.nips.cc/paper_files/paper/19418-/bibtex` (pages and DOI verified
  from this record); OpenAlex W7133259885; camera-ready PDF
  `https://proceedings.neurips.cc/paper_files/paper/2022/file/70596d70542c51c8d9b4e423f4bf2736-Paper-Conference.pdf`
  (all α-ReQ quotes above are from this PDF). Estimator code from
  `https://raw.githubusercontent.com/kumarkrishna/fastssl/main/fastssl/utils/powerlaw.py` and
  `.../scripts/train_model.py`.
- Roy & Vetterli: DBLP `https://dblp.org/rec/conf/eusipco/RoyV07.bib?param=1` (venue, pages, IEEE
  Xplore URL); OpenAlex W2102991522 (DOI, pages); definition verbatim from
  `https://infoscience.epfl.ch/record/110188/files/RoyV07.pdf`.
- **Could not retrieve:** Semantic Scholar (HTTP 429 throughout); OpenReview forum pages for both
  LiDAR and α-ReQ (bot challenge / `ChallengeRequiredError`); numeric δ and ε for LiDAR (absent from
  the paper); any IEEE `10.1109/…` DOI for Roy & Vetterli (none exists on the records checked).

---

## 7. Provenance

- `~/e0_run/P2_METRICS_D2.json`, `P2_METRICS_D1.json`, `P2_METRICS_ALL_SUBSAMPLED.json`,
  `P2_ROBUSTNESS.json`
- `~/e0_run/p2_competing_metrics.py`, `p2_selection_rule.py`, `p2_necessity_and_variance.py`,
  `p2_robustness.py`
- cross-check against `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json` (reproduces to 4 d.p.)
- `~/e0_run/d1_v2/D1_PAIR_MANIFEST.json`, `~/e0_run/d1_audit.log`
