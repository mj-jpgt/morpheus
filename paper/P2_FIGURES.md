# P2 — figure and table plan

> **STALE AGAINST THE 2026-08-04 REWRITE — read this first.** The draft was reorganised around a new
> claim (`P2_RANK_DRAFT.md` §1.3): *effective rank is unusable as a selection signal because its
> between-arm differences are smaller than its own within-arm reproducibility floor.* The old claim
> ("effective rank does not track information content") was falsified by the necessity test and is
> gone. Consequences for this file, to be executed before any figure is drawn:
>
> - **The headline figure is no longer F2.** It is now the reproducibility-floor figure (F3 below) plus
>   a new panel carrying draft §4.1's seven-ratio table against the 2.69× envelope, and a new panel
>   carrying draft §4.2's variance decomposition (34.5% / 65.5% against 98.0%, F 1.41 against 128.2).
>   The variance decomposition is the single most important display item in the paper and currently has
>   no row in this file.
> - **A defeater-check figure is required** (draft §4.4): subsampling SDs, re-export determinism, the
>   fixed-seed probe repeats, and the centring/tolerance insensitivity, on one panel. It is the
>   Leavitt & Morcos §4.2 analogue and a referee will ask for it.
> - **A verdict-instability figure is required** (draft §4.5): statistic × block × view, three small
>   multiples, with the information verdict drawn beside it as a flat reference.
> - **F2 is demoted.** The seed-43 inversion is now reported as implementation-dependent (draft §4.9);
>   any panel showing it must show all four statistic × block combinations or not be drawn.
> - **A figure for the necessity test is required** (draft §4.7) and must be placed *before* the
>   figures that favour the paper, mirroring the draft's order.
> - **F8 is no longer `PENDING` on the rank side.** D1-B's rank column is filled; only the paired
>   bootstrap is pending.
> - **§"Cross-paper deconfliction" is now EXECUTED** — see the end of this file.
>
> Everything below the binding constraints is retained because the per-panel data provenance is still
> correct; only the claim each figure carries and their ordering change.

Companion to `paper/P2_RANK_DRAFT.md`. One row per display item. For each: the **exact data** it must
be drawn from, the **single claim** it carries, and its **status** — `PLOTTABLE` (data on disk or in a
cited markdown table, no new computation), `NEEDS EXTRACTION` (data exists but must be pulled from a
run output that is not yet in a plot-ready file), `NOT MEASURED` (the figure cannot be drawn and the
paper says so in text), or `PENDING` (blocked on D1-B).

Nothing in this file may be drawn from a number that is not in the cited source. If a panel needs a
value that does not exist, the row says `NOT MEASURED` and the draft states the absence in prose
rather than the figure implying it.

Box paths are relative to `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/` on
`ubuntu@150.136.45.194` (persistent NFS). Repo paths are relative to the repository root.

---

## Binding constraints on every figure in this file

These are not stylistic preferences. Each exists because the evidence would be misrepresented without
it.

1. **No panel may place values from two different rank statistics on one axis.** Three mutually
   incompatible functions named `effective_rank` exist in this repository (draft §3.1): **R1** =
   Roy & Vetterli exactly, `exp(−Σ pᵢ ln pᵢ)` on L1-normalised singular values
   (`v2/calibra/spectral.py:14`); **R2** participation ratio of centred singular values
   (`v2/research/rebase/d1_audit.py:149`); **R3** participation ratio of centred *L2-row-normalised*
   singular values (`v2/research/rebase/d1_geometry_probe.py:50`). **Every axis label and legend entry
   must name which one.** Instance 3 (F5) is a *hard numerical* rank — a fourth thing again.
2. **No panel may plot rank change against information change as a scatter pooled across instances.**
   The instances differ in rank statistic, cohort, information measure and units. The overview figure
   F1 is a **small multiple**, one panel per instance with its own axes and units named — never a
   pooled regression.
3. **No number in this paper may be plotted against a published RankMe value.** RankMe normalises with
   `p_k = σ_k/‖σ‖₁ + ε`; Roy & Vetterli use the `0 log 0 = 0` convention; our implementation filters
   at `> 1e-12`. On near-collapsed spectra the three differ measurably (draft §2.6).
4. **Every figure carrying an instance that does NOT contradict RankMe as stated must say so in its
   caption.** That is F2, F5, F6 and part of F1. RankMe restricts itself to same-method comparisons
   and to a necessary-not-sufficient reading; a figure implying otherwise misrepresents the source.
5. **F5's caption must state that the 16/16 column is a hard matrix rank at a structural ceiling of
   16, and the panel must show the counter-measurement (12.88 → 1.00).** A version of F5 showing only
   the 16/16 pinning is the single most misleading figure this project could publish, and this project
   has already described that instance incorrectly once.
6. **Any panel showing a rank *level* must show the reproducibility envelope** from draft §4.7 (3.7×
   across seeds; 2.7× across retrainings at one seed) as a shaded band or an annotated bracket, so no
   reader takes a level comparison at face value.
7. **Every panel carrying instance 1 must be visually de-emphasised** and its caption must state that
   the cited source file does not exist in this repository (draft §4.5, §5.4).

---

## Main figures

### F1 — Six dissociations, ranked by evidential strength, with scope marked

**Claim.** Rank and measured information move independently; the instances are not of equal weight;
and only one of them falls inside the regime RankMe reserves for itself.

**Panels.** A 2 × 3 grid of small multiples **ordered by evidential strength, not by instance
number**, with the strength rank printed in each panel corner and a **scope badge** — "in scope" or
"out of scope" — in the opposite corner. Each panel plots its rank quantity and its information
quantity on **twin axes with both units named**.

- (a) **Instance 6, strength 1, IN SCOPE** — D2 three seeds. Thumbnail of F3.
- (b) **Instance 4, strength 2, out of scope** — dilution dose–response. R1 rank 196.2 → 161.2 against
  null-corrected channel 1.000 → 0.333, both normalised to their d = 0 values so the 18%-versus-67%
  divergence is the panel's content. Thumbnail of F2.
- (c) **Instance 2, strength 3, out of scope** — Phase 1b. Two points: R1 rank 38.48 → 32.06 against
  held-out top-CCA 0.4768 → 0.4748, with the `wsi_identity` internal control (191.07 / 0.5393 in both
  arms) as a second, greyed pair. Caption must state both differences are inside the stack's
  retraining noise.
- (d) **Instance 3, strength 4, out of scope** — hard rank pinned. Thumbnail of F5, drawn with the
  counter-measurement already visible.
- (e) **Instance 1, strength 5, out of scope** — greyed. Rank 49.9 → 103.3 against benchmark
  0.1366 → 0.1367. Caption band: "primary source file does not exist in this repository; benchmark
  statistic undefined; rank statistic unknown."
- (f) **Instance 5** — `[D1 RESULTS PENDING]`. Empty framed panel with labelled axes and "pending"
  printed, so the figure is complete-looking without implying a result. If D1-B does not land before
  submission, replace with draft §4.8.4's outcome table as a text panel.

**Data.** (a) `v2/research/rebase/nature/D2_RESULT.md` §2, §4; `~/e0_run/d2_v3/bootstrap/`,
`~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json`. (b) `DILUTION_LOWER_BOUND.md` §2, §6 and
`NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`; `p1_evidence/dilution/`.
(c) `PHASE1B_TARGETED_READOUT.md` §3, §5; `runs/calibra_v3_targeted`. (d) `NOTEBOOK.md` 2026-08-02
01:20 UTC; `NOTEBOOK_ENTRIES/g26_variance_floor_fix_20260803T0210Z.md`;
`NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md`. (e) `ENGINE_CLD.md` §1;
`HANDOFF_BUILD_AGENT.md` §1–2. (f) pending.

**Status.** (a)–(e) `PLOTTABLE`. (f) `PENDING`.

---

### F2 — The in-scope failure: the rank ordering inverts, with intervals

**Claim.** In a within-method, single-architecture, in-distribution comparison at non-degenerate ranks
— the regime RankMe reserves for itself — a rank-based selection rule picks the right arm once, the
wrong arm once, and cannot choose once, while the performance ordering is stable in all three.

**This is the paper's headline figure.**

**Panels.**
- (a) Three seed groups on a shared axis. **Upper half:** R1 effective rank of arms H and I as paired
  bars (23.39/14.87; 28.77/**34.12**; 9.14/9.11), with the seed-43 inversion and the seed-44 tie
  annotated, and a ✅/❌/❌ "ordering correct?" strip above. **Lower half:** Δ(I − H) on the untrained-40
  readout as a point with **both** the patient CI₉₅ and the cancer CI₉₅ ([−0.1605,−0.0993] /
  [−0.1792,−0.0632]; [−0.1460,−0.0749] / [−0.1623,−0.0118]; [−0.1502,−0.0866] / [−0.1653,−0.0411]),
  and a zero line. **The visual message is that the lower half is flat and the upper half is not.**
- (b) The negative control on the same axes, greyed: Δ = −0.0099 / −0.0280 / −0.0268 with its CIs, so
  the 4–13× separation is visible and the "10–20% may be generic quality" caveat is legible from the
  figure.
- (c) Seed 44 alone, enlarged: two arms equal in rank to two decimals (9.11 vs 9.14) and differing by
  −0.1226 in channel with both CIs excluding zero. **This single matched-rank pair is the cleanest
  refutation in the paper and deserves its own panel.**

**Data.** `v2/research/rebase/nature/D2_RESULT.md` §2, §4; run `d2_v3`; outputs
`~/e0_run/d2_v3/bootstrap/`, `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json`.

**Status.** `PLOTTABLE` from the cited tables.

**Caption must carry.** Both arms use `--objective-profile programme_only` and differ **only** in the
supervision target table; matching is by construction via `D2_PAIR_MANIFEST.json` (`pair_manifest_sha256`
`ce1352e0…`). The rank column is audit check A5, recorded with the instruction "reported, **not
interpreted**", added for an unrelated reason (initialisation collinearity) — i.e. it was not collected
to support this conclusion. Point estimates are not reproducible from the seed alone (draft §3.5); only
the paired within-run difference is quoted. Three seeds is a small number: the finding is that rank
ordering is unreliable, not a quantified error rate.

---

### F3 — The reproducibility floor: rank cannot resolve its own re-measurement

**Claim.** On this stack effective rank moves further between identical retrainings than the
between-configuration differences it would be used to select on.

**Panels.**
- (a) **Strip plot of every available rank level for one experiment.** The six D2 levels (H
  23.39/28.77/9.14, I 14.87/34.12/9.11) plus the re-export/retrain pair (8.68 vs 23.39 at one seed),
  on one axis, with the 3.7× seed spread and the 2.7× retraining spread drawn as brackets. Beside it,
  on a second axis, the three paired channel differences (−0.1325 / −0.1089 / −0.1226) with their CIs.
  **The two axes side by side are the whole argument.**
- (b) **The retraining pair, isolated.** Three markers: recorded original (channel 0.5861), re-export
  of the surviving checkpoint (0.58612, rank 8.68), retrained identical configuration (0.6214, rank
  23.39). Annotate "re-export deterministic to 5 s.f." and "retraining: rank ×2.7, channel +0.035".
  This panel establishes that the variance is in training, not in the readout.
- (c) **The counter-evidence panel, required.** The short-horizon controlled probe: three repeats at
  200 steps, m = 0.999 giving 7.15 / 6.92 / 7.25 (relative spread 4.7%) and m = 0 giving 1.80 / 1.46 /
  1.98 (30%). Rank *can* be tight under a controlled short protocol. The caption must state that this
  shows the variance is accumulated over training rather than intrinsic to the measurement, and that
  three repetitions cannot rule out the gate's 25% divergence rate (`P(0 in 3 | p = 0.25) = 0.42`;
  exact upper 95% bound from 0/3 is `p ≤ 0.63`).

**Data.** (a), (b) `v2/research/rebase/nature/D2_RESULT.md` §4. (c)
`NOTEBOOK_ENTRIES/rank_probe_repeat_variance_20260804T0900Z.md`; logs `~/ws_d1/probevar_*.log`.

**Status.** `PLOTTABLE`.

**Caption must carry.** Two configurations and one stack; no controlled repeat design has been run
(draft §5.2), so this cannot distinguish rank-specific variance from stack non-determinism,
architecture or schedule. The channel readout is not exempt either (0.5861 vs 0.6214 is a 6% move) —
the difference is that the channel is quoted as a paired within-run difference and rank is not.

---

### F4 — The dose–response: rank under-reports the loss by 3.7×

**Claim.** Rank is not merely uncorrelated with the channel; it is miscalibrated in magnitude against
it by a factor that grows with the damage, over a monotone seven-level sweep on a representation with
**zero fitted parameters**.

**Panels.**
- (a) Both curves against achieved dilution *d* (0.000, 0.091, 0.211, 0.302, 0.400, 0.600, 0.800):
  R1 effective rank (196.2, 194.1, 190.5, 187.5, 184.7, 176.5, 161.2) on the left axis and
  null-corrected channel ratio (1.000, 0.999, 0.968, 0.905, 0.804, 0.607, 0.333) on the right. **Both
  axes must be scaled as percentage change from their own d = 0 value**, or the panel understates the
  divergence. Mark the half-loss point at d ≈ 0.68.
- (b) The same data as a single ratio curve — (fraction of rank retained) / (fraction of channel
  retained) — **1.000, 0.990, 1.003, 1.056, 1.171, 1.482, 2.467** across the seven levels. One line,
  one message: the ratio sits at ~1 while almost nothing is happening and rises steeply exactly when
  the damage becomes real. Do **not** describe the curve as monotone — the 0.990 at d = 0.091 is a
  single-seed wobble on a level where the channel moved by 0.001.
- (c) The instrument's own controls over the same levels, to close off "the readout degraded":
  attenuation (1.130, 0.985, 1.003, 1.057, 1.014, 0.855, 0.863) with a reference line at 1, and the
  raw channel ratio as a faint second line to show how much the null correction matters (permutation
  null median 0.145–0.147 at every level).

**Data.** `p1_evidence/dilution/` CALIBRA outputs per level; tabulated in `DILUTION_LOWER_BOUND.md`
§2, §6 and `NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`.

**Status.** `PLOTTABLE`.

**Caption must carry.** **This instance does not contradict RankMe as stated** — high rank with
degraded information is the necessary-not-sufficient case RankMe reserves, and LiDAR's noise-dimension
argument already covers the existence of it. What is added is the *magnitude*. Also: single seed (42),
single draw of donor assignments, **no error bar on level-to-level differences**; detection floor
censored at ≥ 0.40 from d = 0.09 onward; the curve is a property of unweighted mean pooling, not of
the modality; and "lower bound" is withdrawn by the source file's own §4 — the measured quantity is
"the cost of preparation-matched, information-free contamination".

---

### F5 — Instance 3: a withdrawal, drawn honestly

**Claim.** A hard matrix rank at its structural ceiling is insensitive to total collapse; the centred
effective rank of the same objective is not. This figure carries the paper's self-correction.

**Panel.** A single panel, three horizontal tracks sharing a step axis, with the two rank tracks
**stacked so the contrast is unavoidable**.

- **Track 1 — collapse evidence,** arm A (full `programme_free` schedule, 16 fixed patients, 800
  steps): within-modality off-diagonal cosine 0.7089 → **0.9999**; cross-modal positive 0.0538 →
  0.9959 and negative 0.0816 → 0.9960 drawn as two lines converging; retrieval acc@1 0.062 → **0.000**
  with its chance line at 0.062 marked so the *below-chance* endpoint is visible.
- **Track 2 — hard `z_biology` matrix rank: flat at 16/16** across all three arms, with the ceiling
  drawn as a dashed line labelled "structural maximum = batch size 16".
- **Track 3 — centred effective rank on the same objective: 12.88 → 1.00 by step 50**, from `diag_d`
  (clean in-batch InfoNCE), with pos/worst-neg cosine 0.1235 / 0.2186 → 0.9993 / 0.9993 and minimum
  margin −0.219 → −0.0001 annotated at the same steps.

**Data.** Tracks 1–2: `NOTEBOOK.md` entry 2026-08-02 01:20 UTC (source `scratchpad/collapse_diag.py`
on the A100) and `NOTEBOOK_ENTRIES/g26_variance_floor_fix_20260803T0210Z.md`. Track 3:
`NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md` (`~/e0_run/d1_diag/`, `diag_d` trace
at steps 0 / 25 / 50 / 400).

**Status.** Tracks 1–2 `PLOTTABLE` as endpoint markers; track 3 `PLOTTABLE` at four steps.

**`NEEDS EXTRACTION` note.** Tracks 1–2 are recorded as **endpoint pairs** ("0.7089 → 0.9999"), not
per-step arrays. Either (i) extract the per-step arrays from the logs under `~/e0_run/d1_diag/` and
`scratchpad/collapse_diag.py`'s output if retained, or (ii) **draw tracks 1–2 as before/after paired
markers rather than curves, labelled "endpoint values as recorded; per-step array not retained"**. Do
not interpolate.

**Caption must carry.** (i) The 16/16 column is a **hard numerical rank**, not R1/R2/R3; (ii) its
maximum is the batch size of 16; (iii) this is a train batch of 16, not held-out; (iv) the centred
effective rank of the same objective **falls to 1.00**, so this instance is evidence *for* the
collapse-diagnostic use, not against it; (v) this project previously described this instance as one
of its two strongest and that description is **withdrawn**. Without (iv) and (v) in the caption this
figure must not be published.

---

### F6 — The boundary of the collapse diagnostic

**Claim.** Rank near its floor is reliable evidence of total collapse; at 3.6% of nominal
dimensionality it is not evidence of anything.

**Panels.**
- (a) **The collapse regime, where rank works.** Four rows on one axis, each a rank value with its
  co-measured collapse evidence printed beside it: 12.88 → 1.00 (pos/worst-neg cosine 0.9993/0.9993);
  67.55 → ~2 at step 150 (RNA-view mutual cosine 0.9813); 1.76 at epoch 21 and 1.71 at epoch 39 on 282
  held-out patients (RNA–RNA mutual cosine 0.977 / 0.986, hard rank 9 / 11). **Rank statistic labelled
  per row** — the first is a diagnostic-script centred effective rank, the rest are R3.
- (b) **The boundary.** A single large annotated marker: D2 seed 44, R1 effective rank **9.11 and 9.14
  of a nominal 256** (3.6% of ambient), with held-out channels **0.5983 and 0.4757** against a
  permutation null of **0.140**, drawn on a channel axis with the null line marked. This one marker is
  the panel's whole argument and should be its largest element.
- (c) **The cheaper alarm.** Patient-to-patient mutual cosine for every collapsed arm in (a) beside its
  rank, showing that cosine signalled collapse in every case, saturates at a natural maximum of 1, and
  requires no SVD. Supports the draft's recommendation of cosine over rank even for the surviving use.

**Data.** (a) `NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md`;
`NOTEBOOK_ENTRIES/d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`;
`NOTEBOOK_ENTRIES/d1_programme_free_collapsing_in_training_20260803T1930Z.md`;
`NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`.
(b) `v2/research/rebase/nature/D2_RESULT.md` §2, §4. (c) as (a).

**Status.** `PLOTTABLE`.

**Caption must carry.** We have not found a case of total collapse that effective rank missed, so
panel (a) reads "we did not falsify the collapse-diagnostic use", not "we verified it".

---

### F7 — Three statistics, one name

**Claim.** The scalar's definition is not stable across call sites in the codebase that produced these
results, and the historical instances were not all measured with the same one.

**Panels.**
- (a) A comparison on **one synthetic matrix family** — a 282 × 256 matrix interpolated from isotropic
  to rank-1 along the `zᵢ = m + aᵢ·u` family the real collapse actually follows — showing R1, R2, R3
  and hard rank as four curves against the interpolation parameter. **The one panel in this paper that
  requires new computation**; it is cheap, CPU-only and deterministic, and it makes the
  incomparability concrete instead of asserting it. Add RankMe's ε variant of R1 as a fifth curve if
  the ε value can be read from the paper, to show the near-collapse divergence noted in draft §2.6.
- (b) A table-as-figure: which instance used which statistic, with the code location. Instance 6 → R1;
  instance 4 → R1; instance 2 → R1; instance 3 → hard rank; instance 1 → unknown; instance 5 → R2
  (pending), with D1-A's control values in R3.

**Data.** (a) new computation from `v2/calibra/spectral.py:14`, `v2/research/rebase/d1_audit.py:149`,
`v2/research/rebase/d1_geometry_probe.py:50` re-implemented side by side on synthetic input.
(b) draft §3.1.

**Status.** (a) **`NEEDS EXTRACTION`** — new CPU computation. Must be written into the repo as a
script under `v2/research/rebase/` with a test, not a scratchpad file. Run with thread caps
(`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`). (b) `PLOTTABLE`.

**Caption must carry.** Panel (a) is synthetic and illustrative; it is not evidence about any
representation in this paper. Its purpose is to show that the statistics disagree by construction,
which is why draft §3.1 forbids cross-statistic comparison.

---

### F8 — `[D1 RESULTS PENDING]`

**Claim (to be confirmed or withdrawn).** With arms matched by construction, three seeds, and a paired
bootstrap on **both** the rank and the channel, does the rank gap track the channel gap?

**Panels (planned).**
- (a) Per seed: R2 effective rank for `programme_only` and `programme_free` as paired bars, with the
  §4.7 reproducibility envelope shaded behind them.
- (b) Per seed: Δ channel on the 40 untrained targets with patient and cancer CI₉₅, plus the 90
  `random_control` negative control in grey.
- (c) The admission record: which arms cleared the in-runner liveness gate and at what value, against
  the 0.10 threshold and the chance value ln 16 = 2.7726 — because gate admission is a stochastic
  filter (6/8 pass rate over eight identical runs, 650× value spread) and arms that failed it are not
  a random sample.

**Data (planned).** `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json` (headline),
`D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json` (negative control), `D1_AUDIT.json` check A5 (rank column,
**statistic R2**), `D1_READOUT_INDEX.json`. **`D1_PAIRED_BOOTSTRAP.json` must not be used for panel
(b)** — it scores all 90 non-control targets, of which 50 are `programme_only`'s own supervision.

**Status.** `PENDING`. Do not draw any part of this figure from D1-A: D1-A's `programme_free` arm never
trained and its source entry states that nothing about programme supervision may be concluded from it.

**Caption must carry, when drawn.** The two arms differ in *objective*, so D1 sits closer to a
between-method comparison than F2 does and therefore lands further outside RankMe's stated scope. F2
remains the in-scope result whatever D1 shows.

---

## Main tables

### T1 — The six dissociations, ranked, with scope

Draft §4.1's overview table: manipulation, rank statistic, rank change, information change, strength
rank, and **"contradicts RankMe as stated?"**. The last two columns are the table's point and must not
be dropped for space.

**Provenance.** `DILUTION_LOWER_BOUND.md` §2/§6; `D2_RESULT.md` §2/§4; `PHASE1B_TARGETED_READOUT.md`
§3/§5/§7; `NOTEBOOK.md` 2026-08-02 01:20 UTC; `ENGINE_CLD.md` §1 + `HANDOFF_BUILD_AGENT.md` §1–2.
**Status.** `PLOTTABLE`.

### T2 — What RankMe claims and what it restricts

A two-column table of verbatim quotations from arXiv:2210.02885v3: the claims (abstract "indicative
of"; body "RankMe Consistently Predicts Downstream performances From Representations", "a predictor of
representations' performance") beside the restrictions ("a necessary (but not sufficient) condition";
"RankMe should however only be used to compare different runs of a given method"; "there is no
inherent reason for the rank of embeddings to transfer in a monotonic way"; "Except for some
degenerate solutions at full-rank…").

**This table is what makes the paper's scope claim checkable by a referee** and it should be a main
table, not supplementary.

**Provenance.** Full-text PDF of arXiv:2210.02885v3, verified 2026-08-04; quoted in draft §2.1.
**Status.** `PLOTTABLE` (text table).

### T3 — Prior negative results, and what each already establishes

Row per prior negative with a verbatim claim and what it leaves open: LiDAR (arXiv:2312.04000v1 —
"RankMe correlates poorly with downstream performance for most models"; Spearman 0.3174 / Kendall
0.2056 on VICReg at 100 epochs; "a high rank does not guarantee superior performance"); Otero, Mateus
& Balestriero (arXiv:2410.04289v1 — "current methods like RankMe fail to adequately evaluate
representation quality"); Kulkarni et al. (arXiv:2602.20433v2, LLM domain); Cheng (arXiv:2607.13432v1,
plasticity domain). Final column: the regime each does *not* test.

**Exists to make the concession unmissable and to pre-empt "this is Thilak et al. 2023".** Directly
analogous to P1's S1 prior-art map.

**Provenance.** Full-text PDF for LiDAR; abstracts only for the other three (draft §2.6 records the
verification level of each).
**Status.** `PLOTTABLE` (text table). **Row for any additional prior negative found when the §2.2
census is completed must be added before submission.**

### T4 — Rank statistics implemented in this repository

Draft §3.1's table: definition, implementation site, range/maximum, which instances used it, and
whether it equals the published definition.

**Provenance.** `v2/calibra/spectral.py:14-29`; `v2/research/rebase/d1_audit.py:149-153`;
`v2/research/rebase/d1_geometry_probe.py:50-53`; duplicates at `v2/run_rank_ablation.py:35-42` and
`v2/tests/test_stress_collapse.py:23-35`; Roy & Vetterli Definition 1.
**Status.** `PLOTTABLE`.

### T5 — Information measures and their measured chance levels

Draft §3.2's table: measure, definition, chance level, code location. Exists so no reader reads a
channel number against an assumed null of zero. Must include the four distinct InfoNCE chance levels
(ln 16 = 2.7726, ln 80 = 4.38, ln 2576 = 7.854, ln 4310 = 8.369) with the warning that they belong to
different configurations and must not be mixed.

**Provenance.** `v2/calibra/spectral.py`, `v2/calibra/run_calibra.py`, `v2/losses.py:13`,
`NOTEBOOK.md:1554`, `DILUTION_LOWER_BOUND.md` §2, `D2_RESULT.md` §3.
**Status.** `PLOTTABLE`.

### T6 — Reference verification status

Draft §2.6's table, reproduced as a main table if the venue allows, because three fabricated citations
have contaminated this project and a reader is entitled to see the status stated rather than inferred.
Must retain the `[COULD-NOT-VERIFY]` row for α-ReQ's venue and the `INCOMPLETE` row for the prior-art
census.

**Provenance.** draft §2.6; `HANDOFF_BUILD_AGENT.md:156`;
`NOTEBOOK_ENTRIES/winkler_prior_art_20260803T0120Z.md` (the protocol).
**Status.** `PLOTTABLE` (text table).

---

## Supplementary

### S1 — The five-arm collapse sweep: no loss weighting prevents it

R3 effective rank at steps 50/100/150/200/250 for five `programme_free` configurations from one
verified common initialisation (67.55 at step 0), spanning `decorrelation ∈ {0, 0.04}` ×
`biology_full_consistency ∈ {0, 0.1, 1.0}`: (0.04, 1.0) 4.08/1.95/2.16/1.68/1.59; (0, 1.0)
2.62/2.16/2.47/1.94/2.17; (0.04, 0.1) 2.99/3.43/—/—/—; (0, 0.1) 2.97/2.00/2.50/—/—; (0, 0)
2.98/1.98/1.86/—/—.

Carries a claim the main paper only touches: a regulariser family introduced to raise rank does not
prevent this collapse at any setting tested, **including both terms at zero**. Caption must state this
is about our implementation, not about VICReg or Barlow Twins.

**Provenance.** `NOTEBOOK_ENTRIES/d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`; logs
`~/e0_run/d1_diag/`. Centred R3 on a fixed 256-patient held-out probe.
**Status.** `PLOTTABLE`.

### S2 — The decorrelation term's own minimum is collapse

The term-isolation ladder, raw graded contrastive at `(consistency, decorrelation, variance)`:
(0,0,0) 0.00340; (0.01,0,0) 0.04613; (0.1,0,0) 0.08343; (1.0,0,0) 1.84745; (0,0,0.01) 0.53165;
(0,0.004,0.01) 0.13706; (0,0.001,0.01) 2.77288; (0,0.04,0.01) 2.60579; (1.0,0.04,0.01) 2.63086.
Plus the term's value on a healthy batch (38.97) against an all-identical one (1.19e-17), and its
self-extinction 20.74 → 0.00 within 25 steps at every weight 0.001–4.0.

**Provenance.** `NOTEBOOK_ENTRIES/g26_term_isolation_20260803T0930Z.md`;
`NOTEBOOK_ENTRIES/g26_centring_fix_20260803T0730Z.md`; `NOTEBOOK.md` 2026-08-03 decision.
**Status.** `PLOTTABLE`.

### S3 — Per-feature spread fails too, and in the opposite direction

`programme_free` at epoch 21 has **higher** mean per-feature standard deviation than `programme_only`
(0.0137 vs 0.0044) and **lower** effective rank (1.76 vs 7.38), because the collapse is to the family
`zᵢ = m + aᵢ·u` rather than to a point. At epoch 39: 0.0156 vs 0.0056 and 1.71 vs 9.81/10.47.
Isotropic per-feature std for d = 256 is 0.0625.

Exists so the paper cannot be read as recommending per-feature spread as the replacement scalar.

**Provenance.** `NOTEBOOK_ENTRIES/d1_programme_free_collapsing_in_training_20260803T1930Z.md`;
`NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`.
**Status.** `PLOTTABLE`.

### S4 — E1: the experiment this paper should have been built on

A one-page description of `v2/calibra/e1_rank_information.py` and `v2/calibra/aggregate_e1.py`: a
preregistered, gate-enforced, three-seed, equivalence-margin design (margin 0.10) whose aggregated
endpoints are `delta_effective_rank`, `delta_direction_count_above_floor` and
`delta_information_density`, with a paired spike-calibrated detection floor per arm — **never run**.

Exists because it is the honest answer to "why didn't you run the obvious experiment", and because
`information_density = direction_count_above_floor / effective_rank` is itself an example of the
practice under criticism.

**Provenance.** `v2/calibra/e1_rank_information.py` (docstring, `:248`, `:295-298`, `:383-404`);
`v2/calibra/aggregate_e1.py:12,38`; absence verified against `v2/research/rebase/nature/GATE_LOG.md`
and `runs/`.
**Status.** `PLOTTABLE` (text/schematic).

### S5 — Instance 1's provenance gap, stated in full

The chain: `HANDOFF_BUILD_AGENT.md:98` cites `paper/.../RESULTS.md`; no such file exists in the
repository; the numbers survive only as prose in `HANDOFF_BUILD_AGENT.md` §1/§2 and `ENGINE_CLD.md`
§1; "within-cancer specificity" is defined in no file now present; the rank statistic predates the
`spectral.py` consolidation and cannot be assigned to R1/R2/R3.

**Provenance.** as above.
**Status.** `PLOTTABLE` (text table).

---

## Figures the paper does NOT have, and says so

| would-be figure | why it cannot be drawn | where the draft says so |
|---|---|---|
| **LiDAR, or any published alternative criterion, computed on our artifacts** | Identified during this draft's prior-art sweep; requires a labelled probe on every artifact; not run. The most valuable missing measurement in the paper. | §2.5, §5.2 |
| A controlled repeat design for the reproducibility floor | Not run. F3 rests on three seeds plus one accidental retrain and cannot separate rank-specific variance from stack non-determinism, architecture or schedule. | §4.7, §5.2 |
| A pooled scatter of Δrank against Δinformation across instances | Different rank statistics, cohorts, information measures and units. Meaningless, and would look authoritative. | §3.1, §3.6 rule 5 |
| Any of our rank values plotted against a published RankMe value | Three different normalisation conventions (ε inside, `0 log 0`, `> 1e-12` filter) that differ measurably on near-collapsed spectra. | §2.6 |
| Error bars on any dilution rank or channel value | Single seed, single donor draw; the source states there is no error bar on level-to-level differences. | §4.2, §5.2 |
| An equivalence test on instance 2's channel difference | The paired bootstrap the source says "is still required" was never run; "unchanged" means the point estimates differ by 0.002 and nothing more. | §4.4, §5.2 |
| Instance 1 with per-seed values or a run artifact | The cited source file does not exist in the repository; only summary prose survives. | §4.5, §5.4 |
| Per-step curves for F5's tracks 1–2 | Recorded as endpoint pairs, not arrays. Drawn as before/after markers unless recovered from `~/e0_run/d1_diag/`. | F5 `NEEDS EXTRACTION` note |
| Any instance on a second architecture, cohort or modality pair | Not measured. Every number is TCGA, one architecture family, morphology → bulk expression. | §5.2 |
| A case where effective rank **missed** a total collapse | We have not found one; the figure would imply a symmetry the data do not support. | §4.9, §5.2 |
| A rank-versus-channel figure from D1-A | D1-A's `programme_free` arm never trained; its source entry forbids concluding anything about supervision from it. | §4.8.2 |
| E1's three-seed equivalence result | Built and never run. | §3.1, §5.2, S4 |

---

## Cross-paper deconfliction with P1 — **EXECUTED 2026-08-04**

All four edits below were made in this pass. Verification: `paper/P1_FIGURES.md` F11 is replaced by a
`DELETED` stub with the reason; `paper/P1_CALIBRA_DRAFT.md` §4.11 is now a two-paragraph pointer with no
table and no rank numbers; §2.6 carries the verified RankMe and Roy & Vetterli citations, drops the Jing
et al. mis-grouping, and its `[CITATION NEEDED]` is closed in §2.7; `paper/P1_FIGURES.md` F10(b) now
reports both curves with the block named and makes no rank-versus-information claim. **P1 §4.12(iv)
retains its one sentence citing D2's rank values, per the note at the end of this section.** The original
statement of the four edits is preserved below for the record.

### Original statement (for the record)

`paper/P1_FIGURES.md` currently contains **F11 — "Effective rank does not track information content"**,
a single 2 × 2 panel over the same four instances, and `paper/P1_CALIBRA_DRAFT.md` §4.11 contains the
same four-row table. P1's version now also contains **two statements this paper withdraws**.

Four required edits:

1. **Delete P1 F11** and reduce P1 §4.11 to the two sentences it already ends with — that a geometric
   quality metric is computed on the representation rather than through the analysis pipeline whose
   null is in question, and therefore cannot substitute for a sensitivity statement — with a pointer
   to this paper. The table and the F11 panel move here in full.
2. **Correct P1 §4.11's description of instance 3.** It calls instance 3 one of the two "strongest"
   instances and describes it as "rank pinned while information collapses". That column is a hard
   numerical rank at a structural ceiling of 16, and the centred effective rank of the same objective
   falls 12.88 → 1.00 (draft §4.6). The description is not sustainable.
3. **Correct P1 §2.6.** It groups Jing et al. (ICLR 2022) among proposals of "geometric proxies for
   representation quality"; the paper contains no rank→performance claim (draft §2.3). It also carries
   `[CITATION NEEDED: RankMe…]` — **that gap is now closed**: RankMe is Garrido, Balestriero, Najman &
   LeCun, ICML 2023, arXiv:2210.02885, and Roy & Vetterli is EUSIPCO 2007, pp. 606–610, DOI
   10.5281/zenodo.40328. Propagate both, and add that Roy & Vetterli make **no** quality claim.
4. **P1 §4.10 panel F10(b)** ("twin axis: effective rank against null-corrected channel over the
   dilution levels") currently does double duty for F11. It should keep the twin axis, because the
   dilution result is P1's own, but its caption must drop the rank-versus-information *claim* and
   simply report both curves. The claim is made here, in F4.

**P1 §4.12(iv)** ("a capacity explanation is contradicted: effective rank is H 23.39/28.77/9.14 against
PBS 14.87/34.12/9.11") **may stay in P1**, because there it discharges an objection to P1's own
ablation rather than making a claim about rank. F2 of this paper and that sentence of P1 must not both
be described as the paper's finding.

Until edit (1) is executed, **§4 of this paper and §4.11 of P1 must not both be submitted.**
