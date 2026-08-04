## 2026-08-04 23:30 UTC — The floor audit: 23 of 25 selections in P2 are inside a floor, 13 rows have no floor at all, and the sixth late-found instance was P44/I43

**Logged:** 2026-08-04 23:30 UTC. **How obtained:** CPU only, thread-capped
(`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`); the GPU was
running L5/L6 and was not touched. **No rank statistic was computed in this pass.** Every number below
is *read* from a vendored box log, from a JSON readout already in `figures/data/`, or from a named
section of a repository markdown file — the audit module resolves, it does not recompute. Nine new box
logs were vendored through `extract_from_box.py`, the only sanctioned path, and every previously
vendored file re-hashed byte-identical.

### 1. Why

`P2_RANK_DRAFT.md` §4.1's criterion — a rank difference smaller than the measured same-seed retraining
floor is not resolvable — had been applied to this paper's own numbers **five times, each discovered
separately and late**: the momentum fix (3.29×), `m = 0.999` over `m = 0.99` (1.26×), §4.7.4's
surviving necessity violation (3.73×), the decorrelation ablation (1.85×) and §5.1's instance 2
(≈3.2×, flagged in `p2_section5_rewritten_around_the_momentum_replication_20260804T2000Z.md` §5 and
left unrewritten). Finding them one at a time is how a referee finds the sixth. So: enumerate all of
them at once, and make the enumeration a test rather than a document.

### 2. What was built

- **`v2/research/rebase/p2/floor_audit.json`** — 50 comparisons. Per row: both values with a `src`
  block naming the file each came from, the fold, the statistic, the **block**, which floor applies,
  whether it clears, what the claim rests on if it does not, and a `kind`
  (`selection` / `nuisance` / `regime` / `collapse` / `direction` / `excluded`).
- **`v2/research/rebase/p2/p2_floor_audit.py`** — resolves every value back out of its source
  (`json` path, `probe_log` step + column, `spread` over a family, `markdown` literal in a named
  section), recomputes every ratio, recomputes every verdict against the floor **its own block
  licenses**, and renders the table §4.1a prints. It computes no spectral statistic of any kind.
- **`v2/tests/test_p2_floor_audit.py`** — 20 tests. The load-bearing ones are negative: a drifted
  value, a ratio that disagrees with its own two values, a **residualised ratio judged against the raw
  floor**, a missing source file and a markdown literal edited out of the draft each have to be caught.
  It also asserts the draft's copy of the table *is* the rendered table and that the counting sentence
  is the generated one, so the two cannot drift.
- **Draft §4.1a**, table plus prose, and **Appendix A / Appendix C** rows.

### 3. The counts, bad news first

**23 of the 25 selections between candidate configurations do not clear the floor their own block
licenses.** The two that do: §4.4(3)'s fixed-seed probe repeat (3.495×, clears by 6%, and it holds the
**seed** fixed, which §4.2 measures as the dominant term), and one step of a single-seed sweep
(§5.2 step 400, 3.596×). Five rows are exempt with the reason stated; **13 sit on a statistic or a
block for which no floor has ever been measured**; 38 carry an explicit block- or statistic-mismatch
note.

### 4. Six things the audit found that the draft did not already say

1. **No floor exists for most of the paper's rank numbers.** The 3.295× / 3.111× floor is canonical R1
   on the exported `wsi_biology` block and nothing else. R2, R3, PR, RankMe, participation ratio,
   stable rank, α-ReQ, LiDAR and hard rank have none; nor do the fixed held-out probe, in-run training
   batches, the 16-patient gate batch, or the `rna_biology` / `full_biology` views. **Eight of T1's
   twelve metric rows, and every rank number in §5, are in that position.** New §6.2 row.
2. **§4.6's counts have a fourth defect and it is ours.** Beyond n = 6, D2 s44's 1.4 sampling SDs and
   §4.6a's coordinate choice: **not one of the six pairs is resolvable** under the one metric that has
   a floor.
3. **§4.3's headline is the one claim the criterion cannot be applied to.** The compared quantity is a
   *spread* (6.05× against 1.18×); this project has measured a floor on rank **levels** and never on a
   spread.
4. **THE SIXTH INSTANCE: §4.7.4's second violation, P44 against I43, 3.070×.** It is inside the 3.295×
   floor *and* inside §4.2's 2.10–3.75× band, and the draft names it **only** in `P2_FIGURES.md` F6(d)
   with no floor verdict attached. Related: the predeclared violation threshold was **2.0×** — 1.65×
   *smaller* than the floor and below the seed band — so the 66-pair scan could only ever return
   unresolvable violations, which both of the two it returned are.
5. **Two §5.2 claims read an ordering off a 1.04× rank difference.** *"The best-agreeing arm does not
   have the best rank"* is 6.89 against 6.65; *"the discriminating one inverted"* is 3.67 against 3.53
   against a predeclared *"fails ≤ 3.5"* bar. **Both claims survive when restated as what they are** —
   agreement varies 2.06× while rank is *indistinguishable*, and the predicted effect was **absent**,
   not inverted. Rows 45 and 46 are the same two numbers used two ways, which is the clearest
   demonstration in the table that an unresolvable difference cannot carry an ordering and can carry a
   null.
6. **A like-for-like same-seed pair in §5.4's own regime exists**, where §5.4 and §6.2 both say no
   like-for-like measurement does. `~/e0_run/d1_diag/ablate_decorr0.04.log` and `mseed_m0.999_s42.log`
   are the same momentum, decorrelation, capacity, learning rate and seed, from an identical step-0
   state (67.55 / 101.38 / 0.0342 / 0.3650); `d1_momentum_probe.py` runs a **constant** learning rate
   with no schedule keyed to the step budget, so up to step 400 they differ only in GPU
   non-determinism. They span **1.066×** at step 400 and at most **1.128×** over the eight shared
   logged steps. **n = 2 is a pair, not a floor** — §4.1's own argument is that a pair drawn from four
   concordant repeats spans at most 1.028× and would license everything — so it is **recorded and not
   used**, and §5.4, §6.2, F9 and `P2_FIGURES.md`'s pending table now say the accurate sentence
   instead of the merely-true one.

### 5. §5.1's instance-2 residual: resolved as exempt, and §5.1 rewritten to say what it rests on

**§5 is another agent's section and this is a deliberate touch, noted here.** The residual is exempt,
and on a stronger ground than "different regimes": **5.81 is read on a 16-patient *train* batch against
a frozen 64-key queue and ~1.8 on 282 *held-out* patients at cohort scale**, so §3.1's own rule forbids
forming the ratio at all — the ≈3.2× is not a quantity this paper may quote. §5.1 now states that and
names what the instance actually rests on: the gate's **binary** pass (contrastive 0.012–0.057 against
a predeclared ≤ 0.10 bar, retrieval 16/16, three seeds) and the cohort-scale arm being independently
collapsed (RNA-view mutual cosine 0.977 / 0.986, hard rank 9 / 11, sibling arm at 7.38 / 7.35). The
claim would stand unchanged if neither rank number had been recorded.

### 6. One source disagreement — reported, not substituted

**`P2_RANK_DRAFT.md` §5.2 and `paper/QUEUE_ANCHORING.md` say the momentum effect is *"2.6–3.3× at every
step past 150"*. That disagrees with the section's own table.** The per-step m = 0.999 / m = 0 folds
are **3.363× (200), 2.208× (300), 3.596× (400), 3.132× (500), 2.641× (600)** — a range of
**2.208–3.596×** — and **4.343× at step 100**. Both ends of the quoted range are wrong.

Per the rule, **neither number is used**: the audit's step-400 row quotes the table's own two cells,
and the disagreement is recorded in `floor_audit.json`'s `known_source_disagreements` (asserted by a
test) and flagged in §4.1a. **§5.2's prose is not rewritten here** — §5 is another agent's and the
correct range is a judgement about what claim that section wants to make (note that at step 400 the
single-seed fold *does* exceed the floor, while the seed-varied worst case does not, so the choice of
step is not neutral).

### 7. F9 — the dissociation figure that did not exist

`v2/research/rebase/p2/figures/fig_f9_decorrelation.py`, PDF + SVG + 300 dpi PNG, matplotlib only,
Okabe & Ito, statistic and block named on every axis, no value hardcoded — every number is parsed from
a vendored log and the floor is asserted against the envelope log's own printed min and max before
being drawn.

At `m = 0.999`, 400 steps, one seed per level, one verified common initialisation:

| decorrelation | R3 | canonical R1 | RNA-view mutual cosine |
|---|---:|---:|---:|
| 0.0 | 4.32 | 6.29 | 0.4774 |
| 0.01 | 6.22 | 9.32 | 0.7657 |
| 0.04 | **8.01** | **12.20** | **0.8696** |

**Rank rises while a direct measure of the collapse it exists to detect rises with it, monotonically,
co-measured on the same runs and printed on the same log lines.** Stronger than §4.9's material for two
reasons about the *shape* of the evidence: monotone across three levels rather than a single contrast,
and the contradicting quantity is co-measured rather than inferred downstream.

**A statistic correction worth recording.** The 4.32 / 6.22 / 8.01 the earlier entry quotes as
"eff-rank" is the **R3** column — it is what the log's own `final_eff_rank=` line reports. Matched to
the statistic §4.1's floor is measured in, the same three runs read 6.29 / 9.32 / 12.20 and the fold is
**1.940×** rather than 1.854×. Both are inside the floor, so nothing turns on it, but the figure and
the audit name the statistic per axis and per row.

Panels: (a) R3 with the cosine on a twin axis; (b) the same under canonical R1 — **two panels because
binding constraint 1 forbids two rank statistics on one axis**; (c) and (d) the two tracks over
training from the common initialisation. Required annotations drawn inside the artwork: **ONE SEED PER
LEVEL**, and **the ×1.854 / ×1.940 rank change is inside §4.1's ×3.295 floor — the monotonicity and the
co-measured cosine carry this, not the magnitude**, with the floor drawn as a band anchored at each
panel's own decorrelation = 0 value. The same-seed repeat is drawn as open markers labelled *n = 2 is a
pair, NOT a floor*.

**And the conditionality.** *"`feature_decorrelation` is defective"* was **conditional on a
query-written queue**. Without momentum the term aggravated the collapse (1.59 against 2.17 at step
250); with it, it raises rank. §2.4, §4.9, §4.9a, Appendix C and F9's caption now attach *"in the
absence of a momentum key encoder"* to every claim about it.

### 8. Suite

**399 passed**, thread-capped, run from a workspace symlinked to HEAD
(`pytest morpheus/tests morpheus/v2/tests`; `--basetemp` redirected because the default temp root is
not writable on this box). Baseline before this pass was 377. `test_p2_figures.py` needs matplotlib,
which is present locally and **absent from the box venv — nothing was installed into `~/venv`**.

### 9. In plain terms

We wrote down every rank comparison the paper makes and checked each against the same ruler the paper
holds RankMe to. Almost none of them clears it — twenty-three of the twenty-five real comparisons are
smaller than the amount the number moves when you retrain the same model — and for thirteen of the
rows we have never measured the ruler at all, because the floor was only ever measured for one
statistic on one kind of matrix. The audit found one comparison the draft had never judged, two places
where a claim was read off a 4% difference, and two identical runs sitting in the log directory that
we had said did not exist. None of that changes the paper's conclusion; all of it is the kind of thing
a referee finds, and it is much cheaper for us to find it.

### Files / commits

- `v2/research/rebase/p2/floor_audit.json`, `p2_floor_audit.py`; `v2/tests/test_p2_floor_audit.py`
- `v2/research/rebase/p2/figures/fig_f9_decorrelation.py`, `make_all.py`, `extract_from_box.py`
  (nine logs added to `VENDORED`), `figures/data/MANIFEST.json`, `figures/out/F9_*.{pdf,svg,png}`
- `paper/P2_RANK_DRAFT.md` — new §4.1a and §4.9a; status-block items 6 and 7; §2.4, §4.9, §5.1,
  §5.4, §6.2, Appendices A, B and C
- `paper/P2_FIGURES.md` — new F9 and T8 rows; S4 and the pending-dependencies table amended
- `v2/tests/test_p2_figures.py` — F9 added to `MODULES` and to the synthetic corpus
- Sources: `NOTEBOOK_ENTRIES/lr_test_and_decorrelation_reversal_20260804T1130Z.md` §2;
  `retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §1, §3;
  `p2_section5_rewritten_around_the_momentum_replication_20260804T2000Z.md` §5;
  `turnover_criterion_FALSIFIED_20260804T0330Z.md`
