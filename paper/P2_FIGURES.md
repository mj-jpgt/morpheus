# P2 — figure and table plan

Companion to `paper/P2_RANK_DRAFT.md`. **Rewritten 2026-08-04** against the claim that survived the
necessity test; the previous version of this file was written for the falsified claim ("effective
rank does not track information content") and carried a STALE banner. Nothing below is inherited
from it unverified: every path in this file was checked to exist on 2026-08-04, and the ones that
did not are named in "Paths corrected in this rewrite" at the end.

**The claim every figure here serves** (`P2_RANK_DRAFT.md` §1.3):

> Effective rank is unusable as a selection signal because its between-arm differences are smaller
> than its own within-arm reproducibility floor — inside the regime its proponents explicitly
> reserve for it. Its dynamic range is dominated by the one factor that carries no information: the
> training seed.

One row per display item. For each: what it shows, the **exact data file** it draws from, the
**single claim** it carries, and the **draft section** it belongs to.

**Status vocabulary.** `PLOTTABLE` — the numbers are on disk or in a cited markdown table, no new
computation. `NEEDS EXTRACTION` — the data exist but must be pulled from a run output that is not
yet plot-ready. `NOT MEASURED` — the figure cannot be drawn and the draft says so in text.
`PENDING` — blocked on a measurement that has not been made.

Nothing in this file may be drawn from a number that is not in the cited source. If a panel needs a
value that does not exist, the row says so and the draft states the absence in prose rather than the
figure implying it.

**Path conventions.** Repo paths are from the repository root. `~/…` is `ubuntu@150.136.45.194`;
that box's `~/e0_run` is the same tree as
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/e0_run` (persistent NFS), and JSON
files written by the D1/D2 chains record the `/lambda/nfs/…` form of the artifact paths internally.
`~/ws_p2/out/` is the verified-workspace reproduction of 2026-08-04
(`NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260804T0255Z.md`) and is the preferred source for
anything in §4.2–§4.7, because the scripts that wrote it are in the repository at
`v2/research/rebase/p2/`.

---

## Binding constraints on every figure in this file

Not stylistic preferences. Each exists because the evidence would be misrepresented without it.

1. **No panel may place values from two different rank statistics on one axis, and every axis label
   and legend entry must name which one.** Three mutually incompatible statistics have been called
   `effective_rank` in this repository (draft §3.1): **R1** = Roy & Vetterli order 1, centred
   (`v2/calibra/spectral.py`, `CANONICAL`); **R2** = order-2 Hill number of the centred singular
   values, `(Σσ)²/Σσ²`; **R3** = R2 after L2-normalising rows. The "16/16" instance is a **hard
   numerical rank**, a fourth thing again. And a fifth statistic, the order-2 Hill number of the
   *eigenvalues* `(Σσ²)²/Σσ⁴`, appears in the draft under two names — as "participation ratio" in
   §4.6 and, until commit `a11549a`, as "R2" in §4.5(a). **Any panel derived from §4.5(a) must be
   drawn from `~/ws_p2/out/P2_RANK_VARIANTS.json`, and no printing of that table older than
   `a11549a` may be used.**
2. **Every panel showing a rank *level* must show a reproducibility envelope**, as a shaded band or
   an annotated bracket: the **3.295× measured same-seed retraining floor** of §4.1 (residualised
   block; **3.111×** on the raw block — use the one matching the panel's own block), or the arm's own
   seed spread from §4.2's table. A level comparison drawn without it is the practice this paper
   exists to criticise. The superseded n = 1 estimate was 2.69× and may appear only where it is
   labelled as superseded (F1(d)).
3. **No number in this paper may be plotted against a published RankMe value.** RankMe normalises
   with `p_k = σ_k/‖σ‖₁ + ε` (ε outside the division, so the `p_k` do not sum to 1); Roy & Vetterli
   use `0 log 0 = 0`; ours uses a relative LAPACK cut. Draft §2.6, §3.1. The faithful RankMe
   computed *on our own artifacts* for §4.2 and §4.6 is an internal comparison and is fine.
4. **Every figure carrying an instance that does NOT contradict RankMe as stated must say so in its
   caption** — that is F7, F8 and part of F6. RankMe restricts itself to same-method comparisons and
   to a necessary-not-sufficient reading.
5. **No pooled scatter of Δrank against Δinformation across instances.** Different statistics,
   cohorts, information measures and units. Where several instances appear together they appear as
   small multiples with their own axes and units named.
6. **Any panel drawn from arms admitted by the liveness gate must state which arms were admitted**,
   because admission is a stochastic filter: 6/8 pass rate over eight identical runs, 650× value
   spread (draft §5.1).
7. **Panels touching the historical instances must carry their withdrawals**: "16/16" is withdrawn
   as a rank instance and the decorrelation instance is `[NOT RECOMPUTABLE]`. See F8 and T7.
8. **A permutation null drawn on a panel must be the null of that panel's own experiment.** At least
   three are in play (0.140, 0.145–0.147, 0.151–0.158) and they differ in *n*, component count and
   procedure; §4.7.3 mixed two of them and was corrected at commit `9fee55b`. See T5.

---

## Main figures

### F1 — The measured retraining floor, and the seven arm differences inside it

**Draft section.** §4.1. **This is the paper's headline figure.**

**REWRITTEN 2026-08-04.** The controlled repeat that panel (d) used to be a hatched
`[RETRAINING ENVELOPE PENDING]` placeholder for **has reported**, and the script was written to
refuse the placeholder once it did. The envelope is no longer an n = 1 estimate of 2.69×; it is a
measured **3.295×** floor over five identical same-seed retrains, and the figure is rebuilt around
it. The panel letters below are the new ones.

**Claim.** **All seven** of the between-arm rank differences this project has ever measured are
smaller than the spread of the same statistic when one configuration is retrained, at the same seed,
five times — **on the `wsi_biology` view, under canonical R1**. Both conditions are load-bearing and
**T9 is the figure that says so**: on `rna_biology` and `full_biology` the same five retrains spread
1.019× and 1.020× and every between-arm difference clears. F1 must not be captioned as though the
floor were a property of the metric.

**Panels.**

- **(a) The envelope, measured.** Five identical `programme_only` retrains at seed 42 — GPU
  non-determinism the only source of variation — as two stacked strips on a shared repeat axis:
  **rank (R1, residualised, log)** 28.320 / **8.834** / 28.348 / 29.106 / 28.959, spread
  **×3.295**; **channel (top-CCA, 40 untrained targets, linear)** 0.6182 / **0.5859** / 0.6123 /
  0.6110 / 0.6098, spread **×1.055**. The five values are plotted individually and **never as a mean
  or a band**: the distribution is bimodal, and a band would invite the reader to imagine a
  distribution the data does not have. Repeat 2 is annotated. The strips are stacked rather than on a
  twin axis, because the two are different quantities in different units and the comparison being made
  is between their *spreads*.
- **(b) The seven comparisons against it.** A single horizontal axis of *rank ratio*, log-scaled,
  with the floor drawn as a shaded region from 1.0 to **3.295** and each comparison as a labelled
  point: D2 s44 **1.004×**, D2 s43 **1.186×**, Phase 1b **1.200×**, D2 s42 **1.573×**,
  D1-B s44 **1.738×**, D1-B s42 **2.190×**, D1-B s43 **3.246×**. **Every ratio is judged against the
  floor measured on its own block** — 3.295× residualised, drawn dashed; **3.111× raw**, drawn dotted,
  and the Phase 1b triangle is judged against that one. Comparing a residualised ratio against the raw
  floor would put D1-B s43 outside, which is the raw/residualised confusion §4.5 is about. **The
  visual message is all seven points inside the shading.**
- **(c) The asymmetry, beside it.** The same three D2 seeds on a channel axis: paired differences
  **−0.1325 / −0.1089 / −0.1226** with both patient CI₉₅ ([−0.1605,−0.0993] / [−0.1460,−0.0749] /
  [−0.1502,−0.0866]) and cancer CI₉₅ ([−0.1792,−0.0632] / [−0.1623,−0.0118] / [−0.1653,−0.0411]),
  and a zero line. Same sign 3/3, both CIs excluding zero 3/3. **The point of putting (b) and (c)
  side by side is that the channel is quoted as a paired within-run difference and rank is not, and
  there is no paired form of "this run has higher rank".**
- **(d) What the floor replaced, and what it costs us.** The measured floor (**3.295×**) against the
  superseded n = 1 estimate (8.681 → 23.387, **2.694×**), and beneath them D1's three necessity-test
  rank ratios (**3.246× / 2.190× / 1.738×**) on the same axis. All three are inside the floor, so
  **D1 is uninformative about rank in either direction** — and the panel states in text that the
  necessity result is **not refuted** (the channel still separates 3/3 with patient CIs excluding
  zero) and that the count moving 6/7 → 7/7 runs **in our favour** and is reported with the
  scepticism that requires. Whether each D1 point is inside is **computed and stated, never
  asserted**: if a future recomputation moved one outside, the figure must say so rather than fail,
  because that is one of the four outcomes the envelope was predeclared to distinguish.

**Data.** (a) `v2/research/rebase/p2/figures/data/extracted/F1_RETRAINING_REPEAT.json`, parsed by
`extract_from_box.py` out of `~/e0_run/d1_envelope_readout.log` (vendored verbatim at
`v2/research/rebase/p2/figures/data/e0_run/d1_envelope_readout.log`), written by `v2/research/rebase/d1_envelope_readout.py`, which
imports **every** statistic from `v2/calibra`. The extractor records the readout module's git blob
SHA-1 and each `rep{1..5}.npz`'s SHA-256; the figure recomputes the three spreads from the per-repeat
values and **asserts them against the spreads the log itself printed** before drawing anything.
(b) rank ratios from `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`
§5–§6 and `~/ws_rank/RANK_RECOMPUTE.json`; independently reproduced for the ten D2/D1-B values in
`~/ws_p2/out/P2_RANK_VARIANTS.json`; Phase 1b from
`v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §3, §5, §7 and `~/ws_rank/RANK_RECOMPUTE_P1B.json`.
(c) `v2/research/rebase/nature/D2_RESULT.md` §2 and §7; bootstrap outputs
`~/e0_run/d2_v3/bootstrap/`, `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json`.
(d) `v2/research/rebase/p2/figures/data/e0_run/d2_v3/RECOVERED_SEED42_READOUT.json`, `v2/research/rebase/p2/figures/data/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json`.

**Status.** `PLOTTABLE`.

**Caption must carry.** The floor is **five identical retrains of one configuration at one seed on one
stack**, and it is a **floor twice over**: `programme_only` is this project's stable arm, and same-seed
repeats exclude seed variation entirely. It cannot separate rank-specific variance from stack
non-determinism, architecture or schedule. Statistic **R1**, block **residualised** for every D2 and
D1-B point and **raw** for the Phase 1b point, which is judged against the raw-block floor of 3.111×
and whose arms were never verified matched (§3.4). The seventh point clears by **1.5%**. The claim does
not rest on this figure alone: F2 reaches the same conclusion from 8 within-arm degrees of freedom
without using the floor at all. **And the view and the statistic must be named in the caption as
conditions, not as details**: this floor is `wsi_biology` under canonical R1; the same five runs give
1.019× on `rna_biology`, 1.020× on `full_biology`, 2.290× under R3, 1.224× under stable rank and
1.000× under the hard numerical rank, and **1.811× under RankMe as published on the raw block against
our own 3.111×** (T9).

---

### F2 — Where rank's variance lives: arm against training seed

**Draft section.** §4.2. **This is the paper's most important display item and the previous version
of this file gave it no row at all.** It is the contribution that does not depend on a sign count.

**Claim.** Two-thirds of the variation in effective rank across twelve matched artifacts is
training-seed nuisance; two percent of the variation in the information channel is. The arm effect
on rank is not significant; the arm effect on information is overwhelming.

**Panels.**

- **(a) The decomposition itself.** Three stacked 100% bars — canonical effective rank
  (residualised), RankMe as published (raw, ε = 1e-7), and the ground-truth held-out top-CCA —
  each split into **arm** and **seed** shares: **34.5 / 65.5**, **29.1 / 70.9**, **98.0 / 2.0**.
  Print `F(3,8)` beside each bar: **1.41 (n.s.)**, **1.09 (n.s.)**, **128.20**. Rank-type metrics
  are decomposed on the log scale (they are multiplicative, spanning 6.4–34.1), the CCA on the raw
  scale; **the axis label must say so**. Including the RankMe bar is not optional: it is the answer
  to "you evaluated a centred variant", and it is *worse* than ours, not better.
- **(b) The twelve artifacts, per arm.** Four rows (D2-H, D2-I, D1-P, D1-F), each with its three
  seeds plotted on a shared log rank axis and, on a second panel at the same row positions, its
  three channel values on a linear axis. Rank folds **3.15× / 3.75× / 2.64× / 2.10×** against
  channel folds **1.026× / 1.026× / 1.018× / 1.056×**. **The two columns side by side are the
  argument**: within an arm, nothing about the objective, data, split, schedule or architecture
  changed.
- **(c) The single cleanest object, enlarged.** D2 arm H, seed 44 against seed 43: **9.143 against
  28.771 in rank (3.15×)** and **0.5983 against 0.5970 in channel (0.0012 apart)** — with the
  *lower*-rank run marginally ahead. One arm, no arm contrast to argue about, and inside RankMe's
  reserved scope by RankMe's own sentence.

**Data.** `~/ws_p2/out/P2_METRICS_D2.json`, `~/ws_p2/out/P2_METRICS_D1.json`, and the printed
decomposition in `~/ws_p2/out/p2_run.log`; scripts
`v2/research/rebase/p2/p2_competing_metrics.py` and
`v2/research/rebase/p2/p2_necessity_and_variance.py`. Originals
`~/e0_run/P2_METRICS_D2.json`, `~/e0_run/P2_METRICS_D1.json`, tabulated in
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §4.1 and reproduced
exactly in `NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260804T0255Z.md` §2.

**Status.** `PLOTTABLE`. Every value has been computed three times — the original run, the
drift-audit recomputation against `~/ws_d1`, and the 2026-08-04 run from the vendored scripts on a
per-file-verified tree — and agrees to the digit each time.

**Caption must carry.** Four arms and three seeds is 8 within-arm degrees of freedom, not a large
design; `F(3,8) = 1.41` is a *failure to reject*, not a demonstration of no arm effect, and the
caption must say "the arm effect on rank is not resolvable here" rather than "there is none". The
seed changes nothing about objective, architecture, data, split or schedule — that is what makes it
the nuisance factor. Statistic **R1**, block **residualised**; RankMe row is uncentred and raw, per
its own definition.

---

### F3 — The floor is a property of the arm — and, T9 adds, of the statistic and of the view

**Draft section.** §4.3. **Title corrected 2026-08-05:** this figure was specified as *"a property of
the arm, **not** of the statistic"*. That is measurably wrong — §4.1a measures the floor for ten
statistics on one block and finds it running 1.000×–3.295× — and the claim is withdrawn from the
title, from §4.3 and from §1.4's contribution list. F3 keeps the **arm** axis; **T9 carries the
statistic and view axes**, and the two must be read as one result in three parts.

**Claim.** There is no one reproducibility envelope to calibrate and reuse **on any of its three
axes**. This panel shows the arm axis: at the same global
training step, on the same configuration, with only the seed differing, one arm reproduces to three
parts in a thousand and its sibling spans a factor of six.

**Panel.** One panel, two rows on a shared log rank axis at **global step 200 (epoch 11)**:

- `programme_only`: **110.765 / 110.879 / 111.078** — spread **1.003×**, three markers visually on
  top of one another.
- `programme_free`: **7.545 / 45.646 / 12.194** — spread **6.05×**.

Draw each row's spread as a bracket with the fold printed. Do not connect the two rows; do not put
them on the same tick labels beyond the shared axis.

**Data.** `~/e0_run/d1_v2/d1_p_seed{42,43,44}/train_metrics.jsonl` and
`d1_f_seed{42,43,44}/train_metrics.jsonl`, key **`train_rank_tripwire_observed`**, epoch 11.
Tabulated in `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`
§6 and draft §4.3.

**Status.** `PLOTTABLE`.

**Caption must carry.** **Statistic R3** — the in-run tripwire statistic, not R1 — and the key is
`train_rank_tripwire_observed`, with the `train_` prefix, because querying it without the prefix
returns `[]` and `[]` reads as a confident negative
(`NOTEBOOK_ENTRIES/operational_shared_box_rules_20260804T0730Z.md`). These are in-training
measurements whose states were never saved, so they are `[NOT RECOMPUTABLE]` under R1 and **no
number in this figure may appear on an axis with any R1 value in this paper.** The reading is that
rank's reproducibility is worst exactly where the arm is interesting — which is the configuration a
practitioner is looking at when they reach for rank.

---

### F4 — Defeater check: the instability is in training, not in estimation

**Draft section.** §4.4. This is the Leavitt & Morcos (ICLR 2021) §4.2 analogue and a referee will
ask for it; the previous version of this file had no row for it.

**Claim.** The result is not an artifact of a poor rank estimator. Four independent measurements
place the variance in training.

**Panels.**

- **(a) Sampling noise given a fixed trained model.** Patient subsampling at 80%, 40 draws, per
  artifact: six paired points with error bars — 23.325±0.072 / 14.834±0.037; 28.657±0.115 /
  33.912±0.111; **9.132±0.021 / 9.093±0.019**; 29.244±0.102 / 13.393±0.027; 24.584±0.069 /
  7.586±0.016; 11.098±0.031 / 6.387±0.011 — with **gap / sd** printed against each pair
  (105.3, 33.0, **1.39**, 150.3, 238.2, 145.1). **Annotate D2 s44 as the one unresolvable pair**,
  because §4.6's D2 "2/3" depends on it and it is not a hit.
- **(b) Readout determinism.** Two markers: recorded channel 0.5861 against re-export 0.58612;
  "deterministic to five significant figures". Shares its data with F1(a) and may be drawn as an
  inset there instead, but must appear in one of the two.
- **(c) Fixed seed, fixed short horizon.** Three repeats of a controlled 200-step probe with
  identical inputs: m = 0.999 → **7.15 / 6.92 / 7.25** (spread 4.7%); m = 0 → **1.80 / 1.46 / 1.98**
  (30%), with the empty band from 1.98 to 6.92 shaded.
- **(d) Definitional insensitivity.** A small paired plot: canonical R1 **23.387** against faithful
  RankMe (uncentred, ε = 1e-7) **23.391** on H42-residualised, plus the summary that the
  absolute→relative tolerance change moved **no** value (max relative difference `0.000e+00` over 68
  artifact × block combinations) and that uncentred R1 / centred R1 has median **0.995** over those
  same 68.

**Data.** (a) `~/ws_p2/out/P2_METRICS_D2.json`, `~/ws_p2/out/P2_METRICS_D1.json` (`subsample` block),
printed in `~/ws_p2/out/p2_run.log`; original `~/e0_run/P2_METRICS_ALL_SUBSAMPLED.json`, tabulated in
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §3.
(b) `v2/research/rebase/nature/D2_RESULT.md` §4.
(c) `NOTEBOOK_ENTRIES/rank_probe_repeat_variance_20260804T0900Z.md`; logs
**`~/e0_run/d1_diag/probevar_m0.999_{1,2,3}.log`** and **`~/e0_run/d1_diag/probevar_m0_{1,2,3}.log`**.
(d) `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`
§4; `~/ws_rank/RANK_RECOMPUTE.json`; `~/ws_p2/out/P2_METRICS_D2.json` for the RankMe value.

**Status.** `PLOTTABLE`.

**Caption must carry.** Panel (c) constrains the *typical* spread and essentially not the *tail*:
three repetitions cannot rule out the gate's 25% divergence rate
(`P(0 in 3 | p = 0.25) = 0.42`; exact upper 95% bound from 0/3 is `p ≤ 0.63`), and the design was cut
from ten repetitions to three because a ten-way launch exhausted GPU memory. Panel (c)'s statistic is
**R3** and must not share an axis with (a)'s **R1**. Panel (a) and F3 are not in conflict: (a) holds
the model fixed and varies the patient sample, F3 holds the step fixed and varies the seed.

---

### F5 — The verdict is under-determined: statistic, block, view

**Draft section.** §4.5.

**Claim.** "Which arm has the higher effective rank" does not have one answer on our data until
three implementation choices nobody states are fixed — and those choices are worth more than the
between-arm difference they are used to adjudicate. "Which arm carries more molecular information"
is invariant to all three.

**Panels.** Three small multiples, one per degree of freedom, each a 6-column grid of ✓/✗ cells with
the information verdict drawn beneath as a **flat reference row that is correct in every column**.

- **(a) Statistic.** Five rows, columns D2 s42/s43/s44, D1 s42/s43/s44. Canonical
  R1 `OK MISS OK OK OK OK`; canonical **R2 `OK OK MISS OK OK OK`**; canonical
  **R3 `OK MISS OK OK OK OK`** — the three disagree on **2 of 6 pairs (D2 s43, D2 s44)**. Then, drawn
  in a visually separated block below them and **labelled `PR` and `PR_rownorm`, never R2 and R3**:
  `PR` `OK OK MISS MISS OK OK` and `PR_rownorm` `OK OK MISS OK OK OK` — the eigenvalue participation
  ratio `(Σσ²)²/Σσ⁴` and its row-normalised form, which are what an earlier version of this analysis
  computed under the R2/R3 labels. `PR` is identical cell for cell to T1's "participation ratio" row,
  and the panel must say so, because the two are not independent evidence. Draft §4.5(a) was
  corrected to this at commit `a11549a`; earlier printings of that table must not be used.
- **(b) Block.** Raw against confound-residualised, for the two D2 seeds where it flips under R3:
  seed 43 R3 raw H 11.720 / I 11.111 (H) → resid H 14.746 / I 15.915 (**I**); seed 44 R3 raw
  H 5.385 / I 5.449 (I) → resid H 6.564 / I 6.302 (**H**). Show R1 and R2 in the same panel as the
  non-flipping comparison, so the reader sees that the block matters for one statistic and not the
  others. Note in the panel that on the **raw** exported artifacts R2 ≡ R3 exactly, because the
  model already L2-normalises `z_biology` — which is why the distinction went unnoticed.
- **(c) View.** The three co-trained views of the same model — `wsi_biology`, `rna_biology`,
  `full_biology`. Rank winner H/I/I, I/I/I, H/I/I on D2 and P/P/P on all three D1 seeds; information
  winner H/H/H and P/P/P for all six pairs on all three views. Aggregate over the 18 (pair × view)
  comparisons: rank right **11/18**; restricted to D2, **2/9**.

**Data.** (a) `~/ws_p2/out/P2_RANK_VARIANTS.json` — **this file, not the draft table**; script
`v2/research/rebase/p2/p2_rank_variants.py`. (b)
`NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md` §4;
`~/ws_rank/RANK_RECOMPUTE.json`. (c) `~/ws_p2/out/P2_ROBUSTNESS.json`; script
`v2/research/rebase/p2/p2_robustness.py`; original `~/e0_run/P2_ROBUSTNESS.json`, tabulated in
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §4.2.

**Status.** `PLOTTABLE`. (a) must be drawn **from `P2_RANK_VARIANTS.json`**, which is the source
draft §4.5(a) now cites; §4.5's provenance note previously attributed the two source entries'
differing R3 rows to raw-versus-residualised block, which was wrong — both are residualised and the
difference is the statistic — and that too was corrected at `a11549a`.

**(e) Inset — the size of the substitution, on a synthetic spectrum family. `PLOTTABLE`.** The
draft's own framing is that the correction to (a) *is* the paper's result: a fourth statistic living
under a disambiguated name, inside the analysis code of the section that argues the name is
unreliable. This inset makes the size of it concrete rather than asserted. Power-law spectra
`σ_k = k^−a`, 64 components on 256 rows, `a` from 0 to 2 in steps of 0.1, each realised **exactly**
(the matrix is built as an orthonormal basis of the mean-zero subspace scaled by the requested
spectrum, so centring is a no-op and the centred singular values are the requested ones):

| `a` | R1 | **R2** `(Σσ)²/Σσ²` | **PR** `(Σσ²)²/Σσ⁴` | R2/PR |
|---:|---:|---:|---:|---:|
| 0.0 | 64.000 | 64.000 | 64.000 | 1.000 |
| 0.5 | 54.990 | 44.946 | 13.811 | **3.254** |
| 0.9 | 34.517 | 17.922 | 3.026 | **5.923** (max) |
| 1.0 | 29.118 | 13.811 | 2.453 | 5.630 |
| 2.0 | 4.668 | 2.453 | 1.167 | 2.103 |

Draw R2 and PR as two curves against `a`, with the R2/PR ratio on a secondary axis; mark `a = 0`,
where the two agree exactly, so the reader sees that they coincide **only** on a flat spectrum.
Annotate the hand-checkable anchor `σ ∝ (2,1,1)`: **R1 = 2√2 = 2.8284271247**, **R2 = 8/3**,
**PR = 2** — three different numbers on one spectrum. Annotate also the closed-form fingerprint
visible in the table: `PR` at decay `a` equals `R2` at decay `2a` exactly (PR(0.5) = 13.811 = R2(1.0)),
because PR is the order-2 Hill number of the *squared* spectrum and squaring a power law doubles its
exponent.

**Data.** `v2/research/rebase/p2/p2_hill_order_inset.py`, deterministic (seed 20260804), CPU-only,
thread-capped; test at `v2/tests/test_p2_hill_order_inset.py`. Neither statistic is implemented in
that script: `R1`/`R2` come from `spectral.effective_rank` under the named `RANK_VARIANTS`, `PR` from
`p2_competing_metrics.participation_ratio`, and the test asserts the script contains no
decomposition of its own — which is the discipline whose absence produced the error the inset is
about.

**Caption must carry.** This is a synthetic illustration of the *arithmetic*, not evidence about our
artifacts. The empirical size of the substitution is (a)'s table, where the count moved from 3 of 6
to 2 of 6.

**Caption must carry.** These panels do not show that rank is wrong; they show that **a rank verdict
is not a well-defined object until the statistic, the block and the view are stated**. The
`rna_biology` → RNA-derived-target comparison is partly circular and its absolute CCA (0.79–0.85)
must not be read as a clean image→molecular channel; the **rank** measurements on that view are
unaffected by the circularity, but the 2/9 count inherits the caveat and is not quoted as a rate.

---

### F6 — The necessity test, which went against us

**Draft section.** §4.7. **Placed here, before F7 and F8, mirroring the draft's order**, because
§4.7 is the result that falsified the previous framing and the draft reports it before the instances
that favour the paper. It must not be moved to supplementary and it must not be drawn smaller than
F1 or F2.

**Claim.** RankMe's necessary-not-sufficient hedge is not violated by our best-matched three-seed
experiment. It is confirmed, 3/3 — and the gaps at which it is confirmed are themselves inside
rank's own nuisance band.

**Panels.**

- **(a) The predeclaration, drawn.** A two-axis quadrant plot: rank fold on one axis, channel
  difference on the other, with the pre-declared violation region (`fold ≥ 2.0` **and**
  `ΔCCA ≥ 0.0705`) shaded and labelled "the only configuration the hedge cannot absorb". Both
  thresholds came from quantities established independently of this analysis. Plot the three D1-B
  pairs; **all three land in the confirming quadrant**. The predeclaration's date and file must be
  printed in the panel, not only in the caption.
- **(b) The three seeds.** Paired bars: rank `programme_only` / `programme_free` = 29.381/13.418
  (2.19×), 24.673/7.600 (3.25×), 11.115/6.394 (1.74×), with **§4.2's own within-arm seed band
  (2.10–3.75×) shaded behind them** — the whole point being that two of the three ratios are inside
  the retraining envelope of §4.1 and all three are inside the nuisance band. Beside them, channel
  0.6117/0.5412, 0.6198/0.5336, 0.6087/0.5126 with Δ = +0.0705 / +0.0863 / +0.0961 and the paired
  bootstrap CIs: patient CI₉₅ [−0.0938,−0.0444] / [−0.1186,−0.0522] / [−0.1314,−0.0618] and cancer
  CI₉₅ [−0.0957,−0.0180] / [−0.1386,**+0.0006**] / [−0.1535,−0.0016] (signed `programme_free −
  programme_only`). **Seed 43's cancer CI grazes zero and the panel must show that it does.**
- **(c) The negative control, greyed, on the same axes.** `random_control`, 90 targets: Δ = −0.0224 /
  −0.0072 / −0.0322, patient CI₉₅ [−0.0594,+0.0068] / [−0.0480,+0.0278] / [−0.0807,+0.0077], cancer
  CI₉₅ covering zero in 3/3. It separates from (b) by 3–13×, which is what makes (b) readable.
- **(d) The one violation, small and last.** H44 against I43: 3.73× lower rank (9.143 against
  34.117) carrying **+0.1101** more channel — 1 of 66 ordered pairs, plus P44 against I43 (3.07×,
  +0.1206, cross-experiment). Drawn deliberately small, because it is cross-arm and cross-seed and
  RankMe reserves itself to "different runs of a given method".

**Data.** Rank — `~/ws_p2/out/P2_RANK_VARIANTS.json` and `~/e0_run/d1_v2/artifacts/d1_{p,f}_seed{42,43,44}.npz`,
canonical R1 on the residualised held-out `wsi_biology` block. Channel points and the violation scan —
`~/ws_p2/out/P2_METRICS_D1.json` and `~/ws_p2/out/p2_run.log` (necessity scan section); script
`v2/research/rebase/p2/p2_necessity_and_variance.py`. Intervals —
**`~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json`** (40 targets, written 2026-08-04 00:08 UTC)
and **`~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json`** (90 `random_control` targets, 02:40
UTC). Predeclaration — `NOTEBOOK_ENTRIES/PREDECLARED_D1_necessity_test_20260803T2300Z.md`; readout
restriction — `NOTEBOOK_ENTRIES/d1_readout_preregistration_20260803T1700Z.md`. Matching —
`~/e0_run/d1_v2/D1_PAIR_MANIFEST.json` (`"objective_only_difference": true`).

**Status.** `PLOTTABLE`. The bootstrap draft §4.7.2 previously marked `[D1 PAIRED BOOTSTRAP PENDING]`
in six cells **existed all along** — it was hidden by the audit chain's stale absolute path — and was
folded into the draft at commit `a11549a`. The values above are read directly from the two JSON files.
**Both estimators must be drawn, and the conservative one weighted**: the patient bootstrap is
decisive 3/3, the cancer-cluster bootstrap 2/3. A panel showing only the patient interval would be
the selective quotation §4.6 exists to refuse.

**Caption must carry.** `D1_PAIRED_BOOTSTRAP.json` (unsuffixed) **must not be used** — it scores all
90 non-control targets, of which 50 are `programme_only`'s own supervision. The arms differ in
*objective*, so D1 sits closer to a between-method comparison than D2 and lands further outside
RankMe's stated scope; §4.7 is therefore a confirmation of RankMe rather than a refutation of it, and
this figure carries the paper's own negative. The result also trips a preregistered escalation
recorded in `D1_PAIR_MANIFEST.json` — *"if programme_only wins, the collapse story is wrong —
escalate"* — which this paper flags and does not resolve. Panel (d) is **partially pre-empted**:
Aldeneh et al. (ICASSP 2025) have already published that lower-ranked layers can outperform
higher-ranked ones (§2.2).

---

### F7 — The dose–response: rank under-reports the loss, by a factor that depends on preprocessing

**Draft section.** §4.8.

**Claim.** Rank is miscalibrated in magnitude against the information channel by between 1.95× and
21.5× depending on implementation choices, over a monotone seven-level sweep on a representation
with **zero fitted parameters** — so §3.5's retraining noise does not apply and both quantities are
read from the same representation through the same instrument.

**Panels.**

- **(a) Both curves against achieved dilution *d*** (0.000, 0.091, 0.211, 0.302, 0.400, 0.600,
  0.800): **R1 raw block** 196.2 / 194.1 / 190.5 / 187.5 / 184.7 / 176.5 / 161.2 and **R1
  residualised block** 210.2 / 210.2 / 209.5 / 208.4 / 207.0 / 205.1 / 203.7 on the left axis, and
  null-corrected channel ratio 1.000 / 0.999 / 0.968 / 0.905 / 0.804 / 0.607 / 0.333 on the right.
  **Both axes scaled as percentage change from their own d = 0 value**, or the panel understates the
  divergence. **Both rank curves must be shown**: the raw block is where the published numbers came
  from (−17.8%) and the residualised block is the one the channel is measured on (−3.10%), and the
  correction the paper carries is precisely that the two were mixed.
- **(b) The miscalibration table as a figure.** Seven rows — raw R1 3.74×, raw R2 1.95×, raw R3
  2.22×, raw R1-uncentred 3.34×, **residualised R1 21.53×**, residualised R2 11.70×, residualised R3
  14.82× — as a horizontal bar chart against the channel's −66.7%. One message: the direction
  survives every implementation choice and the magnitude does not.
- **(c) The instrument's own controls over the same levels**, to close off "the readout degraded":
  attenuation 1.130 / 0.985 / 1.003 / 1.057 / 1.014 / 0.855 / 0.863 with a reference line at 1, and
  the permutation null median (0.145–0.147 at every level) drawn as a flat line under the raw channel
  ratio.

**Data.** `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2, §6;
`NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`; rank recomputation
`NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md` §5
and `~/ws_rank/RANK_RECOMPUTE.json`; artifact `~/p1_out/dilution/dilution_foreign_tumour_pca256.npz`.
Build scripts in the repository: `v2/research/dilution/`.

**Status.** `PLOTTABLE`.

**Caption must carry.** **This instance does not contradict RankMe as stated** — high rank with
degraded information is the necessary-not-sufficient case RankMe reserves, and LiDAR's own text says
so. What is added is the magnitude. Single seed (42), single draw of donor assignments, **no error
bar on level-to-level differences**; the detection floor is censored at ≥ 0.40 from d = 0.09 onward
and the transmission floor at 0.05 from below; the curve is a property of unweighted mean pooling,
not of the modality; and the source file's own §4 **withdraws the phrase "lower bound"** — the
measured quantity is "the cost of preparation-matched, information-free contamination". Do not
describe (a)'s ratio curve as monotone: the dip to 0.990 at d = 0.091 is a single-seed wobble on a
level where the channel moved by 0.001.

---

### F8 — The use that survives, its boundary, and the withdrawal

**Draft section.** §4.9, §4.10. This figure carries the paper's self-correction and must not be
drawn without it.

**Claim.** Effective rank near its floor (≈ 1–2, with patient-to-patient mutual cosine ≈ 1) is
reliable evidence of total collapse. Anywhere above that — including at 3.6% of nominal
dimensionality — it is uninformative about the channel in both directions. And a **hard** matrix
rank is worthless for the surviving use.

**Panels.**

- **(a) The collapse regime, where rank works.** Three rows, each a rank value with its co-measured
  collapse evidence printed beside it, **statistic labelled per row**: 12.88 → **1.00** at step 50
  (**R3** — see the note below; pos/worst-neg cosine 0.9993/0.9993, min margin
  −0.219 → −0.0001); 67.55 → ~2 at step 150 (**R3**; RNA-view mutual cosine 0.9813); 1.76 at epoch 21
  and 1.71 at epoch 39 on 282 held-out patients (**R3**; RNA–RNA mutual cosine 0.977 / 0.986, hard
  rank 9 / 11).
- **(b) The withdrawal.** The "16/16" instance drawn honestly: a track pinned flat at **16/16** with
  the ceiling drawn as a dashed line labelled "**structural maximum = batch size 16**", stacked
  directly above the **R3** rank of the same objective falling **12.88 → 1.00 by step
  50**. Beside them, the collapse evidence for the same arm: within-modality off-diagonal cosine
  0.7089 → 0.9999, cross-modal positive 0.0538 → 0.9959 and negative 0.0816 → **0.9960** (the
  negatives marginally *higher*), retrieval acc@1 0.062 → **0.000** with its chance line at 0.062
  marked so the below-chance endpoint is visible. **A version of this panel showing only the 16/16
  pinning is the single most misleading figure this project could publish, and this project has
  already described this instance incorrectly once.**
- **(c) The boundary.** One large annotated marker: D2 seed 44, R1 effective rank **9.11 and 9.14 of
  a nominal 256** — 3.6% of ambient — with held-out channels **0.5983 and 0.4757** against a
  permutation null of **0.140**, drawn on a channel axis with the null line marked. This marker is
  the panel's whole argument and should be its largest element.
- **(d) The cheaper alarm.** Patient-to-patient mutual cosine for every collapsed arm in (a), beside
  its rank: 0.977–0.9999 in every case, saturating at a natural maximum of 1, one matrix product, no
  SVD. Supports the draft's recommendation of cosine over rank even for the surviving use.

**Data.** (a), (b) `NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md` (`~/e0_run/d1_diag/`,
`diag_d` trace at steps 0 / 25 / 50 / 400); `NOTEBOOK_ENTRIES/g26_variance_floor_fix_20260803T0210Z.md`;
`NOTEBOOK.md` entry 2026-08-02 01:20 UTC; `NOTEBOOK_ENTRIES/d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`;
`NOTEBOOK_ENTRIES/d1_programme_free_collapsing_in_training_20260803T1930Z.md`;
`NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`.
(c) `v2/research/rebase/nature/D2_RESULT.md` §2, §4; `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json`.
(d) as (a).

**Status.** (a), (c), (d) `PLOTTABLE`. (b) **`PLOTTABLE`, and the extraction resolved it in two
different directions — draw each half accordingly.**

* **The rank track HAS per-step values.** Six recorded steps — 0, 25, 50, 100, 200, 400 — with rank,
  loss, retrieval acc@1, positive cosine, worst-negative cosine, minimum margin and WSI-WSI cosine on
  each. Draw it as a curve through the recorded steps. Source vendored at
  `v2/research/rebase/p2/collapse_tracks/diag_d.log`, parsed by
  `v2/research/rebase/p2/p2_f8b_tracks.py`, pinned by `v2/tests/test_p2_f8b_tracks.py`.
* **The collapse-evidence quantities have ONLY endpoints.** The script that produced the 16/16
  instance probes before and after training, not on a schedule, so no array exists to recover. The
  original fallback stands for this half: **before/after paired markers labelled "endpoint values as
  recorded; per-step array not retained"**. Do not interpolate. Source vendored at
  `v2/research/rebase/p2/collapse_tracks/collapse_diag.log`; every value the panel quotes for arm A
  is asserted against it by the test.

**The statistic, corrected.** The `12.88 → 1.00` column is **R3** — `(Σσ)²/Σσ²` on the column-centred,
**row-L2-normalised** matrix — not `spectral.CANONICAL`. It is computed by an inline formula at
`v2/research/rebase/p2/collapse_tracks/diag_d.py:50-51`, which is vendored verbatim (and allowlisted
in `v2/tests/test_effective_rank_canonical.py` for that reason) so the identification can be read off
the source. Earlier text in this file, and in `paper/P1_CALIBRA_DRAFT.md` §4.11 and `paper/P1_FIGURES.md`,
calls it "the centred effective rank", which in this repository's vocabulary reads as the canonical
order-1 statistic. **The canonical value cannot be recovered** — the diagnostic kept no checkpoint, so
it is `[NOT RECOMPUTED — needs a GPU re-run]`. At the collapsed endpoint the label barely matters (a
rank-1 matrix scores 1.00 under every variant); at 12.88 it is unconstrained.

**Caption must carry.** (i) The 16/16 column is a **hard numerical rank**, not R1/R2/R3, and a
16 × 256 float matrix has full row rank under essentially any perturbation; (ii) its maximum is 16
because the batch is 16; (iii) it is a **train** batch of 16, not held-out; (iv) the centred
effective rank of the same objective **falls to 1.00**, so this instance is evidence *for* the
collapse-diagnostic use, not against it; (v) this project previously listed it among its two
strongest instances and **that description is withdrawn**, here and in P1. Without (iv) and (v) this
figure must not be published. Also: we have not found a case of total collapse that effective rank
missed, so panel (a) reads *"we did not falsify the collapse-diagnostic use"*, never *"we verified
it"*.

---

### F9 — Rank rises while a co-measured collapse measure rises with it

**Draft section.** §4.9a. **NEW 2026-08-04.** The dissociation this figure carries is **stronger than
anything in §4.9**, and until now it had no figure at all.

**Claim.** Raising the covariance-decorrelation weight raises effective rank **monotonically across
three levels** while the RNA-view patient-to-patient mutual cosine — a *direct* measurement of the
condition rank exists to detect — rises **monotonically with it**, co-measured on the identical runs
and printed on the same log lines. Rank reports more occupied directions; the cosine reports the
patients' states converging on one vector.

**Why it is stronger than §4.9's instances**, and both reasons are about the *shape* of the evidence
rather than its size: it is **monotone across three levels** rather than a single contrast, and the
contradicting quantity is **co-measured** rather than inferred from a downstream readout in another
table. §4.9's instances can each be answered with "different instrument, different table"; this one
cannot.

**Panels.**

- **(a) The dissociation, twin axes, statistic R3.** `feature_decorrelation` ∈ {0.0, 0.01, 0.04} on
  the x axis, one point per level. Left axis: **R3** — 4.32 / 6.22 / 8.01. Right axis: **RNA-view
  mutual cosine** — 0.4774 / 0.7657 / 0.8696. Both series rise. R3 is the column these logs' own
  `final_eff_rank=` line reports, so it is the column a notebook entry quoting "eff-rank" is quoting.
- **(b) The same, under canonical R1** — 6.29 / 9.32 / 12.20 against the identical cosine series.
  **Two panels rather than one because binding constraint 1 forbids two rank statistics on one axis**;
  the point of (b) is that the direction survives the choice of statistic, and that R1 is the
  statistic §4.1's floor is measured in.
- **(c), (d) The tracks over training**, steps 0–400, all three arms plus the same-seed repeat —
  rank in (c), cosine in (d) — from **one verified common initialisation** (R3 67.55, canonical R1
  101.38, cosine 0.3650 at step 0 in every arm), so the endpoint is not an artefact of where the
  reading was taken.

**Required annotations, drawn inside the artwork and not left to a caption.**

1. **ONE SEED PER LEVEL.** Printed on the x-axis label of (a) and (b) and again in the bold band
   beneath the panels.
2. **The rank change is ×1.854 (R3) / ×1.940 (canonical R1) and both are INSIDE §4.1's ×3.295
   floor** — drawn as a shaded band anchored at each panel's own decorrelation = 0 value
   (4.32 → 14.23 in (a), 6.29 → 20.73 in (b)), so the reader sees the whole sweep sitting inside it.
   **The monotonicity and the co-measured cosine carry this result, not the magnitude of the rank
   change**, and the bold band says exactly that.
3. The **same-seed repeat** of the 0.04 arm is drawn as open markers, labelled in the legend as
   **n = 2 — a pair, NOT a floor**.

**Data.** `v2/research/rebase/p2/figures/data/e0_run/d1_diag/ablate_decorr{0.0,0.01,0.04}.log` and
`mseed_m0.999_s42.log`, vendored and hashed through `extract_from_box.py` like every other figure's
data (`v2/research/rebase/p2/figures/data/MANIFEST.json`); produced on the box by `v2/research/rebase/d1_momentum_probe.py`, which
imports **both** rank statistics from `v2/calibra` and computes neither inline. The floor is read from
`v2/research/rebase/p2/figures/data/extracted/F1_RETRAINING_REPEAT.json` and the figure **asserts the printed spread against the
log's own min and max** before drawing it. Script `v2/research/rebase/p2/figures/fig_f9_decorrelation.py`;
tests in `v2/tests/test_p2_figures.py`; reported in
`NOTEBOOK_ENTRIES/lr_test_and_decorrelation_reversal_20260804T1130Z.md` §2.

**Status.** `PLOTTABLE`.

**Caption must carry.** (i) **One seed per level, 400 steps, one objective, one cohort.** (ii) The
rank change is **inside §4.1's floor**, and the floor drawn is on a **different arm, duration and
block** — canonical R1 on the residualised exported `wsi_biology` block, `programme_only`, 40 epochs
— while these runs are `programme_free` at 400 steps on a fixed held-out probe **for which no floor
has been measured**, so the band is indicative (§4.1a rows 48–50). (iii) The correction the figure
carries: ***"`feature_decorrelation` is defective" was conditional on a query-written queue.***
Without a momentum key encoder the same term *aggravated* the collapse (1.59 against 2.17 at step
250, §5.1 instance 3); with one it raises rank. **Every claim this project makes about that term
needs *"in the absence of a momentum key encoder"* attached** — §2.4, §4.9, §4.9a and Appendix C now
do. (iv) The `mseed_m0.999_s42` repeat is `n = 2` and **may not be quoted as a floor** for this
regime; §4.1's own argument is that a pair drawn from concordant repeats can license everything.

---

## Main tables

### T8 — The floor audit: every rank comparison the paper makes or relies on

**Draft section.** §4.1a. **NEW 2026-08-04.**

60 rows, one per rank comparison, each with its two values, its fold, its **statistic**, its
**block**, the floor **its own statistic and block license**, whether it clears, and what the claim
rests on if it does not. Of the 25 selections between candidate configurations, **13 fail, 11 are
unjudgeable because no floor has ever been measured for the block they sit on, and 1 clears** — and
the one that clears is **RankMe as published**, whose floor (1.811× raw) is nearly half of canonical
R1's (3.111× raw) on the same five retrains. 29 of the 60 rows are unjudgeable; 5 are exempt with the
reason stated. **Unjudgeable is not a pass and not a failure**, and the audit enforces the
distinction: a row with no floor records `clears: null`.

**A second generated table sits above it in §4.1a**: the floors themselves, one row per statistic per
block per view, measured from the five exported same-seed repeats by
`v2/research/rebase/p2/p2_envelope_floors.py` and rendered by `p2_floor_audit.py --floors`. It
carries min, max, the floor, the other four repeats' own agreement, which repeat diverges and whether
the bimodal shape §4.1 describes survives that change of statistic. It does not: the shape holds
under R1, R2, R3 and residualised RankMe and under nothing else, though every statistic that moves at
all puts the same repeat at the extreme.

**This table is generated, not typed.** `v2/research/rebase/p2/floor_audit.json` is the list;
`v2/research/rebase/p2/p2_floor_audit.py --markdown` renders exactly the table the draft prints, and
`v2/tests/test_p2_floor_audit.py` **re-reads every value out of the file it came from** — a vendored
box log, a JSON readout, or a named section of the draft — and fails if a ratio disagrees with its
source, if a verdict disagrees with its floor, or if the draft's copy of the table has drifted from
the list. A ‡ marks a row whose statistic, block or kind does not match the floor it is judged
against.

**Provenance.** `floor_audit.json` cites, per row, a file under
`v2/research/rebase/p2/figures/data/` or a section of `paper/P2_RANK_DRAFT.md`,
`paper/QUEUE_ANCHORING.md` or a notebook entry. Nothing in it is recomputed.

**Status.** `PLOTTABLE` (text table).

**Must carry.** Block-matching is load-bearing: D1-B seed 43 is **3.246×** residualised and
**3.091×** raw, and its residualised figure against the *raw* floor of 3.111× reads as **outside**
the floor when on its own block it is inside. The table also records one source
disagreement, reported and never substituted, and **now corrected at source**: §5.2's prose said the
momentum effect was *"2.6–3.3× at every step past 150"*, where the section's own table gives
**2.208×–3.596×** over steps 200–600 — **both ends of the quoted range were wrong** — and **4.340×**
at step 100. §5.2 and `paper/QUEUE_ANCHORING.md` now print the corrected range. *An earlier version of
this note and of `floor_audit.json` gave the step-100 fold as 4.343×; the only source in this
repository is §5.2's table, 7.03 / 1.62 = 4.3395, so 4.340× is what may be quoted and 4.343× is
withdrawn.*

**Four rows are new (57–60): §5.2a's learning-rate arms.** They are `direction` rows, not selections,
so the 13 / 11 / 1 split of the 25 selections is unchanged; all four are on the fixed held-out probe
and are therefore **unjudgeable**. They are in the table because §5.2a is the paper's only established
mechanism result and its rank comparisons must be audited like everyone else's — and because the one
that matters most is an **equality**: momentum from 0 to 0.999 moves rank **1.01×** at `lr = 1e-3`,
which is the shape a 1.01× difference can carry and the shape the momentum-threshold account cannot
survive.

### T9 — The floor is a property of the VIEW, and of the STATISTIC — **NEW 2026-08-05**

**Draft section.** §4.1b, and §4.1a's floor table it is drawn from. **This result had no display item
and it is now the claim §4 is organised around**, which is why it gets one.

**Claim.** The reproducibility floor is not a property of effective rank. It is a property of the
**view** the rank is read from and the **statistic** it is read with, on the same five same-seed
retrains, from the same artifacts. **The difference between "unusable" and "usable" is a choice of
co-trained view that no paper we have read reports making.**

**Panel A — the view axis.** Canonical R1, residualised, the five `programme_only` repeats, one row
per co-trained view, with the between-arm resolvability count printed beside each:

| view | min | max | **floor** | divergent repeat | **§4.5(c) pairs resolvable** |
|---|---:|---:|---:|---|:---:|
| `wsi_biology` | 8.8340 | 29.1057 | **3.295×** | rep2, low, **bimodal** | **0 of 6** |
| `full_biology` | 29.8042 | 30.3991 | **1.020×** | rep2, **high** | **6 of 6** |
| `rna_biology` | 27.2245 | 27.7497 | **1.019×** | rep2, low, by **1.9%** | **6 of 6** |

**Panel B — the statistic axis**, one block (exported `wsi_biology`, residualised), same five runs:
**1.000×** hard numerical rank · **1.224×** stable rank · **1.419×** PR · **1.466×** PR_rownorm ·
**2.224×** R2 · **2.290×** R3 · **2.299×** α-ReQ |α−1| · **3.295×** canonical R1 · **3.295×**
residualised RankMe. And on the **raw** block, the row that costs the paper most: **1.811×** RankMe as
published against **3.111×** canonical R1.

**Data.** `v2/research/rebase/p2/figures/data/ws_floor/out/P2_ENVELOPE_FLOORS.json`, measured by
`v2/research/rebase/p2/p2_envelope_floors.py` (CPU only, thread-capped, workspace verified 543/543
files by git blob SHA-1). Every number above is already rendered by
`p2_floor_audit.py --floors`, which generates §4.1a's floor table; the resolvability counts are rows
26–28 of `floor_audit.json`. **Nothing here is recomputed and nothing is typed by hand** — if this
table is drawn rather than printed, it must be drawn from that JSON.

**Status.** `PLOTTABLE` (text table; the two panels are already generated as one table by
`p2_floor_audit.py --floors` and may be printed as-is if the figure budget does not stretch to a
drawn version).

**Caption must carry.** Four things, in this order. **(1)** The floor is a floor **twice over** on
every row — the stable arm, and same-seed repeats that exclude seed variation entirely — n = 5, one
arm, one seed, one stack, no interval; these are not estimated distributions. **(2)** The divergent
repeat is identified by **ten of the ten statistics that move at all**, so nothing hides it; what
changes is the magnitude, and the four-concordant-plus-a-factor shape survives **only** under R1, R2,
R3 and residualised RankMe. That reads as a **redistribution of spectral mass in the tail**, which
entropy-based statistics see by construction and dominant-subspace statistics do not. **(3)** The
`rna_biology` ground truth is partly circular (§4.5(c)) and this table is **not** a recommendation to
switch views: on that view every difference is resolvable and the rank ordering is wrong on **3 of 6**
pairs, against **1 of 6** on `wsi_biology`. **(4)** RankMe's 1.811× against our 3.111× is **on the raw
block only**; on the residualised block, where the mean offset is already gone, the two coincide at
3.295×, which is the evidence for the mechanism rather than a qualification of it.

### T1 — Rank as a selection rule against the published alternatives — **and it is underpowered**

**Draft section.** §4.6.

Twelve metric rows × six matched pairs, with D2, D1, ALL and the exact two-sided binomial *p*:
canonical effective rank raw and residualised 5/6 (p = 0.219); RankMe as published **3/3 on D2 and
1/3 on D1**, 4/6; RankMe residualised 5/6; participation ratio 4/6; stable rank 3/6; α-ReQ |α−1| 4/6;
**LiDAR raw 0/3 on D2**, 3/6; LiDAR residualised 4/6.

**The underpowering is part of the table and must be printed inside it, not relegated to the
caption.** A header note in the table body: *"n = 6. Exact two-sided binomial: 6/6 → p = 0.031;
5/6 → 0.219; 4/6 → 0.688. This design can detect a perfect rule and nothing weaker. No comparison
between two rows of this table is evidence in either direction."* Additionally mark **D2 s44's cell
in every row as unresolvable** (F4(a): 1.4 sampling sd), so effective rank's honest D2 record reads
"1 clear hit, 1 clear miss, 1 pair it cannot resolve".

**A SECOND BAND, added 2026-08-04, and it is larger than the first: the ground truth is a coordinate
choice.** Draft §4.6a. Every mark in the table is scored against the held-out channel onto 40
**gene-set** targets, and that arm contrast exists on the gene sets and on **none** of the five other
molecular target blocks on disk
(`NOTEBOOK_ENTRIES/d2_coordinate_system_result_20260804T0800Z.md` §1a). The band prints, for six
blocks as columns: the **arm ordering across seeds 42/43/44** (`HHH` on gene sets, shuffled PBS,
`random_control` and random dictionary; **`HIH`** on PBS codes and PCA basis — seed 43 is the seed the
coordinate system flips), and the **D2 count** for the two rows §4.6 quotes against one another,
with each cell's ALL/6 and exact binomial beneath it. The band's closing sentence is the finding:
*the ordering between those two rows reverses on 2 of 6 blocks, and canonical effective rank reaches
6/6 — p = 0.031, "significant" by this table's own bar — produced by nothing but the choice of which
coordinate system the exam is written in. No count in this table may be quoted without its target
block.* The band also states that the **D1 half cannot be re-scored** — those arms were never scored
against any block but the gene sets — so only the D2 count moves, and that is an absent measurement
rather than evidence of block-stability.

**Provenance.** `~/ws_p2/out/P2_METRICS_D2.json`, `~/ws_p2/out/P2_METRICS_D1.json`, printed table in
`~/ws_p2/out/p2_run.log`; scripts `v2/research/rebase/p2/p2_competing_metrics.py`,
`v2/research/rebase/p2/p2_selection_rule.py` and — for the band —
`v2/research/rebase/p2/p2_selection_rule_blocks.py` over
`v2/research/rebase/nature/d2_coordinate_system/out/EXAM_PANEL.json`; originals
`~/e0_run/P2_METRICS_D2.json`, `P2_METRICS_D1.json`, tabulated in
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §3. The band's counts are
recomputed in the figure from the same per-artifact metrics the table's marks are recomputed from, and
the exam panel's gene-set column is asserted against the metrics JSON's own untrained-40 contrast
before anything is drawn.
**Status.** `PLOTTABLE`.

**Must carry.** LiDAR is **adapted** — q = 2 with the two modalities as views, licensed by the
paper's own footnote 4 but outside the authors' tested q range — and its ordering is invariant across
δ from 1e-8 to 1e0. α-ReQ follows the authors' released `fastssl` estimator, not the paper text, and
`|α − 1|` is our operationalisation of their "Goldilocks zone", not theirs; all 12 artifacts sit at
α between 2.6 and 4.8, far outside it. **And the ground-truth band**: no count in the table is
quotable without naming its target block.

### T2 — What RankMe claims and what it restricts

Two columns of verbatim quotation from arXiv:2210.02885v3: the claims ("indicative of";
"RankMe Consistently Predicts Downstream performances From Representations"; "a predictor of
representations' performance") beside the restrictions ("a necessary (but not sufficient)
condition"; "RankMe should however only be used to compare different runs of a given method";
"there is no inherent reason for the rank of embeddings to transfer in a monotonic way"). **This
table is what makes the paper's scope claim checkable by a referee** and belongs in the main body.
It is also what licenses F2 and F3: three seeds of one arm *are* "different runs of a given method".

**Provenance.** Full-text PDF of arXiv:2210.02885v3, verified 2026-08-04; quoted in draft §2.1.
**Status.** `PLOTTABLE` (text table).

### T3 — Prior negative results, and what each already establishes

Row per prior negative with a verbatim claim and the regime it does *not* test: LiDAR
(arXiv:2312.04000v1, **ICLR 2024** — "RankMe correlates poorly with downstream performance for most
models"; Spearman 0.3174 / Kendall 0.2056 on VICReg at 100 epochs); Otero, Mateus & Balestriero
(arXiv:2410.04289v1); Kulkarni et al. (arXiv:2602.20433v2, LLM domain); Cheng (arXiv:2607.13432v1,
plasticity domain); **Aldeneh et al., ICASSP 2025** — "lower-ranked layers can outperform
higher-ranked ones", which partially pre-empts F6(d) and must appear here rather than only in prose.

**Exists to make the concession unmissable and to pre-empt "this is Thilak et al.".** Directly
analogous to P1's S1 prior-art map. **Provenance.** Full-text PDF for LiDAR; abstracts only for the
others (draft §2.6 records the verification level of each); venue corrections in
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §2; census in
`NOTEBOOK_ENTRIES/p2_prior_art_citation_graph_sweep_20260803T2326Z.md` (453 de-duplicated works).
**Status.** `PLOTTABLE` (text table).

### T4 — The statistics this repository has called `effective_rank`

Draft §3.1's table: definition, implementation site, range, which instances used it, whether it
equals the published definition — **plus a fifth row for `(Σσ²)²/Σσ⁴`**, the eigenvalue
participation ratio, which is what §4.6 calls "participation ratio" and what §4.5(a) mislabelled as
R2. Include the measured ratios over 68 artifact × block combinations: R2/R1 median 0.629 (min
0.338), R3/R1 median 0.655, R1-uncentred/R1 median 0.995.

**Provenance.** `v2/calibra/spectral.py` (`CANONICAL`, `RANK_VARIANTS`);
`v2/research/rebase/d1_audit.py`; `v2/research/rebase/d1_geometry_probe.py`;
`v2/research/rebase/p2/p2_competing_metrics.py` (`participation_ratio`, and the docstring stating
why it is not R2); `v2/tests/test_effective_rank_canonical.py`; Roy & Vetterli Definition 1.
**Status.** `PLOTTABLE`.

### T5 — Information measures and their measured chance levels

Draft §3.2's table: measure, definition, chance level, code location — so no reader reads a channel
number against an assumed null of zero. Must include the four distinct InfoNCE chance levels
(ln 16 = 2.7726, ln 80 = 4.38, ln 2576 = 7.854, ln 4310 = 8.369) with the warning that they belong to
different configurations.

**And it must carry the permutation nulls the same way, because this paper has already mixed two of
them.** At least three appear on the project — **0.140** (200-draw row shuffle of the residualised
target matrix, D2), **0.145–0.147** (within-cancer, dilution sweep), and 0.151–0.158 elsewhere — and
they differ in *n*, component count and, for 0.140, in the permutation procedure itself. §4.7.3
compared D2's random controls against 0.147 and was corrected to 0.140 at commit `9fee55b`. **No
figure may draw a null line from a different experiment than the values it sits under**, and F8(c)'s
0.140 and F7(c)'s 0.145–0.147 must never appear on one axis. Draft §3.2 also records an unresolved
label conflict in our own notes — one entry calls 0.140 a within-cancer null, another a row shuffle;
the value is agreed and the procedure label is not — and the table should carry that as stated rather
than pick one.

**Provenance.** `v2/calibra/spectral.py`, `v2/calibra/run_calibra.py`, `v2/losses.py`,
`NOTEBOOK.md:1554`, `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2,
`v2/research/rebase/nature/D2_RESULT.md` §3. **Status.** `PLOTTABLE`.

### T6 — Reference verification status

Draft §2.6's table, as a main table if the venue allows, because three fabricated citations have
contaminated this project and a reader is entitled to see the status stated rather than inferred.
Must retain the `INCOMPLETE` row for the prior-art census and record that LiDAR's ICLR 2024 venue and
α-ReQ's NeurIPS 2022 venue were resolved on 2026-08-03.

**Provenance.** draft §2.6; `NOTEBOOK_ENTRIES/winkler_prior_art_20260803T0120Z.md` (the protocol);
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §2.
**Status.** `PLOTTABLE` (text table).

### T7 — The historical instances, recomputed: what survives, what is qualified, what is withdrawn

Draft §4.9's table: instance, manipulation, **statistic and block**, rank, information, verdict now —
five rows (D2, dilution, Phase 1b, "16/16", decorrelation) with the last two marked **WITHDRAWN** and
**`[NOT RECOMPUTABLE — artifact never existed]`**, plus D1-A's 9.81/1.71 marked `[NOT RECOMPUTED]`.
This is the table F8(b) and S5 exist to support.

**Provenance.** `v2/research/rebase/nature/D2_RESULT.md`;
`v2/research/rebase/nature/DILUTION_LOWER_BOUND.md`;
`v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §3, §5, §7; readout `runs/calibra_v3_targeted`;
`NOTEBOOK.md` 2026-08-02 01:20 UTC; `v2/research/rebase/ENGINE_CLD.md` §1 and
`HANDOFF_BUILD_AGENT.md` §1–2; recomputation
`NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md` §5.
**Status.** `PLOTTABLE`.

---

## Supplementary

### S1 — The five-arm collapse sweep: no loss weighting prevents it

R3 effective rank at steps 50/100/150/200/250 for five `programme_free` configurations from one
verified common initialisation (67.55 at step 0), spanning `decorrelation ∈ {0, 0.04}` ×
`biology_full_consistency ∈ {0, 0.1, 1.0}`: (0.04, 1.0) 4.08/1.95/2.16/1.68/1.59; (0, 1.0)
2.62/2.16/2.47/1.94/2.17; (0.04, 0.1) 2.99/3.43/—/—/—; (0, 0.1) 2.97/2.00/2.50/—/—; (0, 0)
2.98/1.98/1.86/—/—. Carries a claim the main paper only touches: a regulariser family introduced to
raise rank does not prevent this collapse at any setting tested, **including both terms at zero**.

**Provenance.** `NOTEBOOK_ENTRIES/d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`; logs
`~/e0_run/d1_diag/`. **Statistic R3**, centred, fixed 256-patient held-out probe. **Status.**
`PLOTTABLE`. **Caption:** about our implementation, not about VICReg or Barlow Twins.

### S2 — The decorrelation term's own minimum is collapse

The term-isolation ladder, raw graded contrastive at `(consistency, decorrelation, variance)`:
(0,0,0) 0.00340; (0.01,0,0) 0.04613; (0.1,0,0) 0.08343; (1.0,0,0) 1.84745; (0,0,0.01) 0.53165;
(0,0.004,0.01) 0.13706; (0,0.001,0.01) 2.77288; (0,0.04,0.01) 2.60579; (1.0,0.04,0.01) 2.63086. Plus
the term's value on a healthy batch (38.97) against an all-identical one (1.19e-17), and its
self-extinction 20.74 → 0.00 within 25 steps at every weight 0.001–4.0.

**Provenance.** `NOTEBOOK_ENTRIES/g26_term_isolation_20260803T0930Z.md`;
`NOTEBOOK_ENTRIES/g26_centring_fix_20260803T0730Z.md`. **Status.** `PLOTTABLE`.

### S3 — Per-feature spread fails too, and in the opposite direction

`programme_free` at epoch 21 has **higher** mean per-feature standard deviation than `programme_only`
(0.0137 against 0.0044) and **lower** effective rank (1.76 against 7.38), because the collapse is to
the family `zᵢ = m + aᵢ·u` rather than to a point. At epoch 39: 0.0156 against 0.0056 and 1.71
against 9.81 / 10.47. Isotropic per-feature std for d = 256 is 0.0625. Exists so the paper cannot be
read as recommending per-feature spread as the replacement scalar.

**Provenance.** `NOTEBOOK_ENTRIES/d1_programme_free_collapsing_in_training_20260803T1930Z.md`;
`NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`.
**Statistic R3.** **Status.** `PLOTTABLE`.

### S4 — The worked example: momentum, graded in rank — and the seed replication that does not exist

**Draft section.** §5.2, and §5.4's admission.

Centred R3 on a fixed held-out probe, capacity held at 4,096 in every arm, one verified common
initialisation, steps 0→600 (40 epochs of this objective is 583 steps): m = 0 gives
67.55/4.10/1.62/1.62/2.26/3.32/2.18/2.43/2.81; m = 0.9 gives 67.55/3.88/3.51/2.15/1.65/2.70/2.31/
2.34/2.23; m = 0.99 gives 67.55/8.65/6.49/4.56/5.70/6.01/5.50/5.50/5.88; m = 0.999 gives
67.55/9.35/7.03/6.99/7.60/7.33/7.84/7.61/7.42. Draw all four as curves; annotate that m = 0.9 tracks
the no-momentum arm, so the threshold lies between 0.9 and 0.99, and that both working arms are flat
from step 200 to 600.

**Provenance.** `NOTEBOOK_ENTRIES/queue_size_implicates_the_key_set_20260803T2200Z.md`;
`NOTEBOOK_ENTRIES/momentum_rescues_rank_but_staleness_is_not_the_mechanism_20260803T2330Z.md`; logs
`~/e0_run/d1_diag/`.

**Status.** `PLOTTABLE` as drawn — **and the seed replication has now landed, so the panel's marker
must change rather than be dropped. Draft §5.2 and the new §5.4 are rewritten around it; this panel is
not yet redrawn.** The original sweep is **one seed per momentum value**, and §4.2
measures the seed term as dominant for exactly this statistic on this stack. The single seed was a
**defect, not a design choice** — the momentum harness had its seed hardcoded, so the sweep could not
have varied seeds had we asked it to — and the panel must say which of the two it was.

**The replication (2026-08-04, `~/e0_run/d1_diag/mseed_*`, three seeds per momentum, 500 steps):**
canonical R1 **11.26 / 10.45 / 10.55** at m = 0.999 against **3.18 / 1.13 / 2.36** at m = 0; R3
7.40 / 6.85 / 7.15 against 2.81 / 1.05 / 2.06. **Every m = 0.999 seed exceeds every m = 0 seed on both
statistics**, so the predeclared disjunction (now recorded in §5.4) resolves in favour of separation and the single-seed defect is
closed. Replace `[MOMENTUM SEED REPLICATION PENDING]` with a **per-arm seed band drawn the way F2(b)
draws one**.

**But the panel must now carry the awkward number instead.** The worst-case separation is
**min(m = 0.999) / max(m = 0) = 10.45 / 3.18 = 3.29×, against F1's measured retraining floor of
3.295×**. **That floor is on a different block**, and the fixed held-out probe these runs are read
on has never had one measured, so the correct word is **unjudgeable**, not "not resolvable": the
criterion cannot reach the comparison in either direction, and the panel must say that rather than
print a failure it cannot support. On R3 it is 6.85 / 2.81 = 2.44×, likewise unjudgeable. **The panel
must also point at S9**: the mechanism behind this fix is not momentum but the learning rate, and the
threshold in `m` this panel draws is a property of the one learning rate it was measured at. Three caveats make this indicative rather than a like-for-like
disqualification, and none of them rescues it: the floor was measured on `programme_only` at 40 epochs
on an exported artifact, and these runs are `programme_free` at 500 steps on a held-out probe —
different arm, duration and block, and no like-for-like floor for this regime has been measured — though one same-seed **pair** in it now has been, and is concordant (1.066× at step 400; §4.1a row 56). **n = 2 is a pair, not a floor**, and the panel must not draw it as one.
**What does not depend on rank at all is the reason the fix was adopted** (the unfixed configuration
collapses and the fixed one does not, visible in retrieval, in the contrastive loss, and in whether
`programme_free` reaches 40 epochs without the tripwire firing), and the panel must not be drawn as
though it did. **Draft §5.4 states that binary outcome in the form the panel should annotate**:
`programme_free` completed 40 epochs uncollapsed **0 of 3** seeds before the fix — the one arm that
reached epoch 40 sat at R3 **1.71** with RNA-view mutual cosine **0.986**, and `run_d1` produced no
exports, CALIBRA or bootstrap at all — and **3 of 3** after it, at canonical R1 13.418 / 7.600 / 6.394
with a held-out channel of 0.5412 / 0.5336 / 0.5126 above its own `random_control` and a paired
bootstrap interval per seed. **The choice of `m = 0.999` over `m = 0.99` (1.26×) is supported by
nothing that clears the floor; only momentum against none is.** *Numbers:
`NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §3;
`d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`; draft §3.3, §4.7.2, §4.7.3.
Draft §5.2 and §5.4 **are** now rewritten around this.*

**Caption must carry.** MoCo's staleness account is ruled out three ways, and the anchoring
alternative is stated as unconfirmed and now also as **incomplete** (S9 shows this collapse is
primarily a learning-rate phenomenon, and an anchoring account has to explain why the same decoupling
is unnecessary at `4e-5` and useless at `1e-3`): the queue turns over completely every 19 steps; key-to-encoder cosine at step
100 is 0.427 / **0.908** / 0.441 for m = 0 / 0.99 / 0.999 against ranks 2.58 / 6.65 / **6.89**, so the
best-agreeing arm does not have the best rank; and at fixed key encoder the strongest arm is capacity
**64**, entirely overwritten every step (rank 6.17 against 2.16 at capacity 4,096). MoCo advances its
account twice **as a hypothesis** (arXiv:1911.05722, verified at full text) and a falsification must
say so.

### S5 — The gate is a stochastic filter

Eight runs with identical inputs, identical seed and an identical 2,400-step budget gave
`final_biology_contrastive` values spanning **650×** — 0.00859, 0.01076, 0.01770, 0.02019, 0.02407,
0.03266, 0.38009, 5.58511 — a **6/8** pass rate against an unchanged 0.10 threshold, and a bimodal
shape: six clustered within 4×, two divergent. At a 75% per-arm pass rate the probability of all
three contrastive arms clearing is 0.42. Plus the harness-versus-runner inversion: seeds 42/43/44 gave
0.01871 / 0.01206 / 0.05666 standalone and passed / **0.50883 ✗** / **2.14122 ✗** inside the runner
(chance is ln 16 = 2.7726).

Exists because binding constraint 6 requires it, and because the replacement rank probe of F4(c) is
better behaved on *spread* but rests on a thinner basis on *divergence rate*.

**Provenance.** `NOTEBOOK_ENTRIES/g26_is_not_reproducible_20260804T0700Z.md`;
`NOTEBOOK_ENTRIES/d1_relaunch_20260803T1530Z.md`;
`NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`.
**Status.** `PLOTTABLE`.

### S6 — E1: the experiment this paper should have been built on

A one-page description of `v2/calibra/e1_rank_information.py` and `v2/calibra/aggregate_e1.py`: a
preregistered, gate-enforced, three-seed, equivalence-margin design (margin 0.10) whose aggregated
endpoints are `delta_effective_rank`, `delta_direction_count_above_floor` and
`delta_information_density`, with a paired spike-calibrated detection floor per arm — **never run**.
Exists because it is the honest answer to "why didn't you run the obvious experiment", and because
`information_density = direction_count_above_floor / effective_rank` is itself an example of the
practice under criticism.

**Two corrections to this panel, from the 2026-08-04 aptness audit** (full account in
`NOTEBOOK_ENTRIES/PREDECLARED_E1_aptness_and_verdict_20260804T0609Z.md`):

1. **It is not the experiment this paper should have been built on, and the panel must stop saying so.**
   E1's question is *"does decorrelation-created rank carry molecular information"* — the claim that was
   falsified and removed. The surviving claim is about **resolvability against a within-arm
   reproducibility floor**, and E1 contains no within-arm term at all: its only intervals resample
   patients and cancer clusters, the variance component §4.1(iv) measures at SD ≈ 0.1 against a
   retraining floor of 3.295×. Worse, `aggregate_e1.py:38,48` makes `delta_effective_rank > 0` a
   necessary conjunct of its verdict, so a sign-unstable delta — the *supporting* observation for the
   surviving claim — is reported as "claim not supported". The panel keeps its point (we built the
   practice we criticise) and loses the counterfactual.
2. **"It is CPU work" is wrong.** E1's inputs do not exist. `v2/research/rebase/run_e1_training.py`
   is the only driver that can produce them and it runs `morpheus.v2.runner` twice per seed at
   `--decorrelation-weight` 0 and >0. Running E1 costs **six GPU trainings**. Every available
   artifact pair carries `decorrelation_weight = 0.04` in *both* arms, so E1's own
   `_validate_intervention` refuses all of them.

**Provenance.** `v2/calibra/e1_rank_information.py`; `v2/calibra/aggregate_e1.py`;
`v2/research/rebase/run_e1_training.py`; absence verified
against `v2/research/rebase/nature/GATE_LOG.md` and `runs/`. **Status.** `PLOTTABLE`
(text/schematic).

### S7 — Instance 1's provenance gap, stated in full

The chain: `HANDOFF_BUILD_AGENT.md:98` cites `paper/.../RESULTS.md`; **no such file exists in the
repository**; the numbers (rank 49.9 → 103.3 against "within-cancer specificity" 0.1366 → 0.1367)
survive only as prose in `HANDOFF_BUILD_AGENT.md` §1–2 and `v2/research/rebase/ENGINE_CLD.md` §1;
"within-cancer specificity" is defined in no file now present; the rank statistic predates the
`spectral.py` consolidation and cannot be assigned to R1, R2 or R3.

**Status.** `PLOTTABLE` (text table). **This instance appears in no main figure**, is excluded from
every count in the paper, and belongs in a history paragraph.

### S8 — One implementation, and the test that keeps it

Schematic of the consolidation: ten call sites, three statistics, one function; the hand-computed
pin (`σ ∝ (2,1,1)` ⟹ `erank = 2√2 = 2.8284271247461903`, order-2 = 8/3); object identity asserted
across every importable call site; and the AST + SVD scan that fails the build if a second definition
or an unallowlisted SVD-based rank reappears. Add the P2 scripts' own entry: they are now in
`v2/research/rebase/p2/` with `v2/tests/test_p2_analysis_scripts.py` running all five on synthetic
input.

**Provenance.** `v2/calibra/spectral.py`; `v2/tests/test_effective_rank_canonical.py`;
`v2/tests/test_p2_analysis_scripts.py`; `v2/research/rebase/p2/README.md`.
**Status.** `PLOTTABLE` (schematic).

---

### S9 — The learning-rate test: the mechanism, and three dead accounts before it — **NEW 2026-08-05**

**Draft section.** §5.2a.

**Claim.** The collapse §5.2 fixes with a momentum key encoder is primarily a **learning-rate**
phenomenon. Learning rate separates the outcomes perfectly; momentum separates none of them at the
high rate. Momentum is **neither necessary nor sufficient**.

**Panel.** Six arms as a 2 × 2 (+2) grid, all 400 steps from the same verified initialisation
**67.55**, centred **R3** on the fixed held-out probe, with the RNA-view mutual cosine printed in each
cell:

| | m = 0 | m = 0.9 | m = 0.99 | m = 0.999 |
|---|---|---|---|---|
| **lr 1e-3** | L3 — **1.06** (cos 0.9946) | L1 — **1.05** (cos 0.9257) | — | **L5 — 1.05** (cos 0.5207) |
| **lr 4e-5** | **L6 — 12.30** (cos 0.8199) | — | L2 — **27.88** (cos 0.5436) | L4 — **35.24** (cos 0.3807) |

**L5 and L6 are the load-bearing cells and must be marked as such**: L5 fails at maximal momentum and
L6 succeeds at zero momentum, which is the pair the predeclaration named as decisive. **L1–L4 were run
and read first, against a predeclaration that discriminated a different pair of accounts; the 2 × 2
above is what showed that design could not separate momentum from learning rate, and the two empty
cells were predeclared before they ran.** The panel must show the empty cells as having been empty.

**Provenance.** `~/e0_run/d1_diag/lr_L{1..6}.log`, `v2/research/rebase/d1_momentum_probe.py` with an
`lr` argument (module default `2e-4`, the project's training rate). Predeclaration
`NOTEBOOK_ENTRIES/PREDECLARED_learning_rate_test_20260804T2200Z.md` (`f68a7ac`); L1–L4 in
`lr_test_and_decorrelation_reversal_20260804T1130Z.md`; L5/L6 in
`NOTEBOOK_ENTRIES/learning_rate_is_the_mechanism_20260805T0100Z.md`.

**Status.** `PENDING VENDORING` — the six logs are **not** in the vendored `d1_diag` directory, unlike
the `ablate_decorr*` and `mseed_*` logs beside them, so no value in this panel is re-parsed from a copy
in this repository. It is `PLOTTABLE` from the draft's table and must not be drawn as though its
provenance were equal to F9's. Closing it is a file copy; draft §6.4 states the gap.

**Caption must carry.** **(1)** This is the **fourth** account proposed for this collapse and the
**first to survive a predeclared test**; regulariser weighting, MoCo staleness and the `τ/T` turnover
criterion were each falsified, the last two by experiments built to test them. **(2)** At the training
rate actually used, `2e-4`, momentum **does** rescue the objective, and the seed-varied replication of
that (S4) is unambiguous — **lowering the learning rate would have solved the original problem more
simply and we did not try it**. **(3)** One seed per cell, on a block with **no measured floor**: all
four contrasts are **unjudgeable** (§4.1a rows 57–60) and the result rests on a predeclared sign, not
on a magnitude. **(4)** The cosine column is the reason the panel prints two numbers per cell: rank is
moved 11.6–33.6× by the learning rate and at most 2.9× by momentum, while the cosine on the same six
runs is moved 1.2–1.4× by the learning rate and **1.9–2.2× by momentum** — the two instruments order
the two knobs oppositely. The cosine is **uncentred** and the rank is **centred**, and that alternative
account is not excluded (§4.10, §6.2).

## Figures the paper does NOT have, and says so

| would-be figure | why it cannot be drawn | where the draft says so |
|---|---|---|
| ~~**A controlled repeat design for the envelope** — N retrainings of one configuration with rank and channel measured on each~~ | **CLOSED at N = 5 (2026-08-04), and this row was stale.** F1 panel (a) is that design: five identical `programme_only` retrains at seed 42, rank and channel measured on each, giving a **3.295×** floor against a **1.055×** channel spread. It replaced the n = 1 estimate of 2.69×. What it still cannot do is attribute retraining variance to the metric rather than to this stack, and it is measured on the stable arm at a fixed seed, so it is a floor rather than an envelope. F2 remains the panel that does not depend on it. | §4.1, §4.4, §6.2 |
| ~~A seed replication of S4's momentum sweep~~ | **REPORTED 2026-08-04, and it does not clear F1's floor.** Three seeds per momentum, 500 steps; every m = 0.999 seed exceeds every m = 0 seed, closing a single-seed defect that was **hardcoded in the harness**, not a design choice. But the worst-case separation is **3.29×** against a **3.295×** floor measured on a *different block*, and the block these runs are read on has no floor at all — so the comparison is **unjudgeable**, not failing. S4 and §5.4 now carry that word, and the draft rests the fix on a binary training outcome (§5.4) rather than on the ratio. | §5.2, §5.4, §6.2 |
| **The six `lr_L*.log` learning-rate logs, vendored** | **Not vendored, and this is the weakest provenance in the paper's only established mechanism result.** `ablate_decorr*` and `mseed_*` were copied into the vendored `d1_diag` directory under `v2/research/rebase/p2/figures/data/` and every value read from them is re-parsed from the copy by `v2/tests/test_p2_floor_audit.py`. `lr_L{1..6}.log` were not, so §5.2a's four audit rows resolve against **the draft's own table** — the weakest of the three source kinds the audit supports. **Closing it is a file copy; nothing needs re-running.** | §5.2a, §6.4, S9 |
| **The RNA-view mutual cosine on the CENTRED representation, for the three `lr = 1e-3` arms** | **Not measured.** Rank reads 1.06 / 1.05 / 1.05 across m = 0 / 0.9 / 0.999 while the uncentred cosine reads 0.9946 / 0.9257 / 0.5207. Either rank is insensitive at the collapse floor or the whole difference is in the mean-offset direction centring removes — the same asymmetry that makes RankMe more reproducible than our statistic (T9, §4.1b). Needs the runs' activations, which the logs do not carry. | §4.10, §6.2, S9 |
| **A like-for-like retraining floor for the momentum regime** | **Not measured — but one same-seed PAIR in that regime now has been, and it is concordant** (`ablate_decorr0.04` against `mseed_m0.999_s42`: 1.066× at step 400, ≤ 1.128× over the shared steps; draft §4.1a row 56). **n = 2 is a pair, not a floor**, and §4.1's own argument is that a pair drawn from concordant repeats can license everything, so it is recorded and not used. F1's floor is `programme_only`, 40 epochs, residualised exported block; S4's runs are `programme_free`, 500 steps, fixed held-out probe. Every quotation of the 3.29× / 3.295× comparison must say so — and must also say that all three mismatches point toward a *larger* floor in this regime, so the mismatch is not a reprieve. | §5.4, §6.2 |
| **A labelled linear probe on every artifact** | Not run. It is the reference standard RankMe and LiDAR were validated against; ours is a held-out canonical correlation against unsupervised molecular targets. | §3.2, §6.2, §6.3 |
| Error bars on any dilution rank or channel value | Single seed, single donor draw; the source states there is no error bar on level-to-level differences. | §4.8, §6.2 |
| An equivalence test on Phase 1b's channel difference | The paired bootstrap its own source says "is still required" was never run; "unchanged" means the point estimates differ by 0.002 and nothing more. | §4.9, §6.2 |
| Instance 1 with per-seed values or a run artifact | The cited source file does not exist in the repository; only summary prose survives. See S7. | §4.9, §6.4 |
| Per-step curves for F8(b)'s **collapse-evidence** quantities (cosines, retrieval) | **Confirmed unrecoverable.** `collapse_diag.log` records endpoint pairs only, because the script probes before and after training rather than on a schedule. F8(b)'s **rank** track, by contrast, *does* have per-step values and is now extracted. | F8 status |
| A rank-versus-channel figure from D1-A | D1-A's `programme_free` arm never trained; its source entry forbids concluding anything about supervision from it. | §4.9 |
| D1-A's 9.81 / 1.71 recomputed under R1 | Needs a GPU forward pass from surviving checkpoints (`~/e0_run/d1_v1/d1_{p_seed42,p_seed43,p_seed44,f_seed42}/last.pt`; the `f_seed43` and `f_seed44` directories exist but hold **no `last.pt`** — the gate refused those arms) and the GPU was in use. `[NOT RECOMPUTED]`. Draft §4.9 writes these paths as `d1_{p42,…}`, which do not resolve. | §4.9 |
| Any instance on a second architecture, cohort or modality pair | Not measured. Every number is TCGA, one architecture family, morphology → bulk expression. `claim_guards.no_external_cohort` is undischarged. | §5.2, §6.2 |
| A case where effective rank **missed** a total collapse | We have not found one; the figure would imply a symmetry the data do not support. | §4.10, §6.2 |
| A pooled scatter of Δrank against Δinformation across instances | Different statistics, cohorts, measures and units. Meaningless, and would look authoritative. | §3.1, binding constraint 5 |
| Any of our rank values against a published RankMe value | Three normalisation conventions that differ on near-collapsed spectra. | §2.6, §3.1 |
| E1's three-seed equivalence result | Built and never run — **and it should not be run for this paper's surviving claim.** E1 tests whether decorrelation-created rank carries information; the surviving claim is about resolvability against a within-arm floor, which E1 measures nowhere. Its aggregator makes `delta_effective_rank > 0` a necessary conjunct of its verdict, so the supporting outcome of the current claim is scored as "claim not supported". It is also unrunnable: no artifact pair has the preregistered decorrelation 0 → >0 intervention. See S6 and `NOTEBOOK_ENTRIES/PREDECLARED_E1_aptness_and_verdict_20260804T0609Z.md`. | §3.1, §6.2 |

---

## Pending dependencies, collected

| item | figure | state as of 2026-08-04 03:10 UTC |
|---|---|---|
| D1 paired bootstrap, 40 targets | **F6(b), F6(c)** | **RESOLVED.** `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json` (00:08 UTC) and `D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json` (02:40 UTC); folded into draft §4.7.2 at commit `a11549a`. Draw **both** estimators; weight the cancer-cluster one. |
| §4.5(a)'s statistic labels | **F5(a)** | **RESOLVED.** Corrected in the draft at `a11549a`; plot from `~/ws_p2/out/P2_RANK_VARIANTS.json`. |
| **Momentum seed replication** | **S4** | **REPORTED 2026-08-04. The draft is now rewritten around it (§5.2, §5.4); S4 itself is still to be redrawn.** Three seeds per momentum at 500 steps (`~/e0_run/d1_diag/mseed_*`): canonical R1 11.26 / 10.45 / 10.55 at m = 0.999 against 3.18 / 1.13 / 2.36 at m = 0, so every m = 0.999 seed exceeds every m = 0 seed on both statistics and the single-seed defect is closed. **But the worst-case separation is 3.29× against F1's 3.295× floor**, so S4 must carry that the fix's *rank* difference is not resolvable by §4.1's own criterion — different arm, duration and block, so indicative rather than a like-for-like disqualification, and all three mismatches point toward a larger floor in this regime. **S4 must also carry what the fix does rest on**, per new draft §5.4: `programme_free` completing 40 epochs uncollapsed **0 of 3** seeds before the fix and **3 of 3** after, with a channel and paired bootstrap intervals where no export existed at all. Numbers in `NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §3. |
| **Controlled retraining repeat design** | **F1(a), F1(b), F1(d), F4, F6(b)** | **RESOLVED at N = 5, 2026-08-04.** `~/e0_run/d1_envelope/rep{1..5}.npz`, readout `~/e0_run/d1_envelope_readout.log`, vendored and extracted to `v2/research/rebase/p2/figures/data/extracted/F1_RETRAINING_REPEAT.json`. Rank floor **3.295×** residualised / **3.111×** raw against a **1.055×** channel spread. Do **not** draw it as a fitted distribution or a confidence band: it is **bimodal** (four repeats within 2%, one at a third) and it is a **floor**, measured on the stable arm at a fixed seed. F1(a) plots the five values individually for that reason. |
| **Per-block ground truth for the selection rule** | **T1** | **RESOLVED for the D2 half, 2026-08-04.** `v2/research/rebase/nature/d2_coordinate_system/out/EXAM_PANEL.json` re-scores both D2 arms on six target blocks; `v2/research/rebase/p2/p2_selection_rule_blocks.py` re-runs the rule with each as truth. Every metric row's D2 count moves; the ordering between canonical effective rank and RankMe reverses on two blocks. **Unresolved for the D1 half** — those arms were never scored on any block but the gene sets, so T1's band holds D1 fixed and says so. |
| Per-step arrays for the collapse tracks | F8(b) | **RESOLVED, split.** The rank track has per-step values and is extracted (`v2/research/rebase/p2/collapse_tracks/diag_d.log`); the collapse-evidence quantities have endpoints only and keep the paired-marker treatment. Statistic corrected to **R3**. |
| Synthetic `(Σσ)²/Σσ²` against `(Σσ²)²/Σσ⁴` inset | F5 | **RESOLVED.** `v2/research/rebase/p2/p2_hill_order_inset.py` + `v2/tests/test_p2_hill_order_inset.py`; R2/PR spans 1.000 to 5.923 over `a ∈ [0, 2]`. |
| Phase 1b equivalence test | T7 | Not run; the row stays qualified, and "unchanged" means the point estimates differ by 0.002. |
| Labelled linear probe on every artifact | — | Not run. No figure depends on it; §6.2 and §6.3 state the absence. |

---

## Paths corrected in this rewrite

Checked 2026-08-04. Every path in this file exists unless the row says it does not.

- `~/ws_d1/probevar_*.log` (previous F3) — **does not exist**. The probe-repeat logs are at
  `~/e0_run/d1_diag/probevar_m0.999_{1,2,3}.log` and `~/e0_run/d1_diag/probevar_m0_{1,2,3}.log`,
  which is what draft §4.4(3) says. Now cited at F4(c).
- `DILUTION_LOWER_BOUND.md`, `PHASE1B_TARGETED_READOUT.md` (previous F1, F4, T1) — bare filenames
  that do not resolve from the repository root. They are
  `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` and
  `v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md`.
- `p1_evidence/dilution/` (previous F1, F4) — **not in this repository**. The dilution CALIBRA
  outputs are P1's; the artifact this paper reads is
  `~/p1_out/dilution/dilution_foreign_tumour_pca256.npz` and the build scripts are
  `v2/research/dilution/`.
Two stale pointers **in `paper/P2_RANK_DRAFT.md`**, not fixed here because the draft is another
agent's, and flagged in
`NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260804T0255Z.md` §5:

- Draft **line 1699** (§4.9) cites the surviving D1-A checkpoints as
  `~/e0_run/d1_v1/d1_{p42,p43,p44,f42}/last.pt`. Those paths do not resolve: the directories are
  `d1_p_seed42`, `d1_p_seed43`, `d1_p_seed44`, `d1_f_seed42`. The draft's substantive statement is
  right — `d1_f_seed43` and `d1_f_seed44` exist as directories and hold no checkpoint.
- Draft **line 2234** (Appendix B) cites the vendored rank-recomputation scripts as
  `v2/research/rebase/rank_recompute_all_instances.py` and `v2/research/rebase/rank_recompute_phase1b.py`. **Neither exists**; the files are
  `v2/research/rebase/rank_recompute_all_instances.py` and
  `v2/research/rebase/rank_recompute_phase1b.py`.

---

## Cross-paper deconfliction with P1 — **EXECUTED 2026-08-04**

All four edits were made in the 2026-08-04 rewrite pass and are retained here for the record.
Verification: `paper/P1_FIGURES.md` F11 is replaced by a `DELETED` stub with the reason;
`paper/P1_CALIBRA_DRAFT.md` §4.11 is now a two-paragraph pointer with no table and no rank numbers;
§2.6 carries the verified RankMe and Roy & Vetterli citations, drops the Jing et al. mis-grouping, and
its `[CITATION NEEDED]` is closed in §2.7; `paper/P1_FIGURES.md` F10(b) now reports both curves with
the block named and makes no rank-versus-information claim.

The four edits, as originally stated:

1. **Delete P1 F11** and reduce P1 §4.11 to the two sentences it already ends with — that a geometric
   quality metric is computed on the representation rather than through the analysis pipeline whose
   null is in question, and therefore cannot substitute for a sensitivity statement — with a pointer
   to this paper. The table and the F11 panel move here in full.
2. **Correct P1 §4.11's description of instance 3.** It called it one of the two "strongest"
   instances and described it as "rank pinned while information collapses". That column is a hard
   numerical rank at a structural ceiling of 16, and the **R3** rank of the same objective
   falls 12.88 → 1.00. See F8(b) — P1's own wording, "the centred effective rank", needs the same
   correction, in `P1_CALIBRA_DRAFT.md` §4.11 and `P1_FIGURES.md`.
3. **Correct P1 §2.6.** It grouped Jing et al. (ICLR 2022) among proposals of "geometric proxies for
   representation quality"; that paper contains no rank→performance claim. Its
   `[CITATION NEEDED: RankMe…]` is closed: Garrido, Balestriero, Najman & LeCun, ICML 2023,
   arXiv:2210.02885; Roy & Vetterli, EUSIPCO 2007, pp. 606–610, DOI 10.5281/zenodo.40328, who make
   **no** quality claim.
4. **P1 §4.10 panel F10(b)** keeps its twin axis, because the dilution result is P1's own, but its
   caption drops the rank-versus-information *claim* and simply reports both curves. The claim is made
   here, in F7.

**P1 §4.12(iv)** ("a capacity explanation is contradicted: effective rank is H 23.39/28.77/9.14
against PBS 14.87/34.12/9.11") **may stay in P1**, because there it discharges an objection to P1's
own ablation rather than making a claim about rank. F1 of this paper and that sentence of P1 must not
both be described as the paper's finding.
