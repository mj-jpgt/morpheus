# PREDECLARATION — the exported-artifact retraining floor on the UNSTABLE arm (`programme_free`, 40 epochs)

**Logged: 2026-08-05 00:45 UTC.** Committed **before** the five `programme_free` retrains are
launched and before any number below exists. Nothing here is measured; every value quoted as existing
is re-read from a named file that predates this entry.

---

## 0. Why this exists, and why it is the paper's own objection turned on its headline

`P2_RANK_DRAFT.md` §4.1's floor — **3.295×** canonical R1 on the residualised exported `wsi_biology`
block, and **3.111×** raw — is five identical **`programme_only`** retrains at seed 42, 40 epochs
(`~/e0_run/d1_envelope/rep{1..5}.npz`, `chain_retrain_envelope.sh`). `programme_only` is this
project's **stable** arm, and the draft says so itself, in the same paragraph that quotes the number:
*"Measuring the retraining spread on the stable arm understates it."*

Since then the paper has measured what "understates" is worth. On the fixed held-out probe, ten
same-seed repeats of **both** arms
(`NOTEBOOK_ENTRIES/the_probe_block_has_a_floor_at_last_20260804T1620Z.md` §3) give:

| step | floor, both arms | stable (m = 0.999) arm alone | ratio |
|---:|---:|---:|---:|
| 100 | 1.333× | 1.112× | 1.20 |
| 200 | 2.057× | 1.143× | 1.80 |
| 250 | 1.933× | 1.165× | 1.66 |
| 400 | 1.570× | 1.101× | 1.43 |
| 500 | 1.367× | 1.089× | 1.26 |

*"The collapsed arm carries the floor at every step and under both statistics, by roughly a factor of
two. Had this been measured on one arm the way §4.1's exported floor was, it would have read
1.09×–1.17× and would have flattered every row in the audit."*

The completeness audit of 2026-08-04 20:30 filed this as open item 3, in these words: *"Every
exported-block floor is `programme_only`. The probe measurement showed the collapsed arm carries the
floor by ~2×, so §4.1's 3.295× is very likely an underestimate by about that factor — and the
direction of that error is against the paper's own headline count."*

**So this is not a robustness check. It is the paper applying to its own most-quoted number the
scope rule it applies to everybody else's.**

---

## 1. Exactly what will be run, and what will not change

**Protocol: `~/chain_retrain_envelope.sh`, unchanged except for one flag.**

```
python -m morpheus.v2.runner \
  --data-config  ~/e0_run/data/v1_abs_hallmark.json \
  --split-file   ~/e0_run/data/paired_split_maximal.json \
  --output-dir   ~/e0_run/d1_envelope_pf/rep$r \
  --objective-profile programme_free            <-- THE ONLY CHANGE (was programme_only)
  --epochs 40 --token-budget 8192 --hidden-dim 512 --layers 4 --heads 8 \
  --learning-rate 2e-4 --weight-decay 1e-2 --decorrelation-weight 0.04 \
  --variance-weight 0.01 --separation-weight 0.01 --loss-warmup-epochs 4 \
  --programme-head-dim 256 --restrict-to-split --seed 42 --device cuda \
  --fit-development --fixed-final-epoch \
  --biology-key-momentum 0.999 --rank-tripwire-step 200 --rank-tripwire-minimum 4.0 \
  --gate-repeats 0 --rank-probe-repeats 0 \
  --expected-development-cancers 11 --expected-heldout-cancers 21
```

then `morpheus.v2.export` on each `last.pt`, byte-for-byte the same export invocation, to
`~/e0_run/d1_envelope_pf/rep$r.npz`.

Held fixed and stated so a later reader can check them: **5 repeats** (matching n = 5, not more);
**seed 42** in every repeat, so the only source of variation is GPU non-determinism, exactly as on the
stable arm; **40 epochs**; the same data config, split file, architecture, optimiser, schedule,
momentum and tripwire; the same export path and the same `wsi_biology`/`rna_biology`/`full_biology`
views; **five runs concurrent on one A100**, as the stable arm's five were.

**Readout, also unchanged.** `v2/research/rebase/p2/p2_envelope_floors.py --reps
'.../d1_envelope_pf/rep*.npz'` — the module that produced every published exported-block floor. Every
statistic is imported (`calibra.spectral` for R1/R2/R3 and the top-CCA channel;
`p2_competing_metrics` for RankMe / PR / stable rank / α-ReQ / LiDAR; `numpy.linalg.matrix_rank` for
the hard rank). **No statistic is written inline and no new module defines one.**
`v2/research/rebase/d1_envelope_readout.py` is run as the second, independent source for R1 and the
channel, exactly as it is for 3.295×.

**Combining the arms.** The convention already fixed in `p2_probe_floors.combine()` is used and not
re-invented: the combined floor is **`max` of the two arms' folds, with the arm that carried it
named** — never a pool and never an average. Per-arm folds are reported beside it.

**What this will still be.** A floor **once** over, not twice: same-seed repeats still exclude seed
variation entirely (§4.2 measures that axis and finds it larger). n = 5 per arm, one seed, one stack,
no interval. It is not a distribution and will not be written as one.

---

## 2. Prediction, written before the runs

| quantity | prediction | interval I would be surprised to fall outside |
|---|---|---|
| `wsi_biology` R1 residualised, `programme_free` alone | **≈ 5×**, larger than 3.295× | 2.5× – 12× |
| — its shape | **not bimodal**; graded across the five, as on the probe block | — |
| `wsi_biology` R1 raw, `programme_free` | ≈ 4.5× | 2.5× – 12× |
| **`rna_biology` R1 residualised, `programme_free`** | **≈ 1.10×** | 1.02× – 1.8× |
| **`full_biology` R1 residualised, `programme_free`** | **≈ 1.15×** | 1.02× – 2.0× |
| channel (top-CCA 16), `programme_free` | ≈ 1.10× | 1.02× – 1.4× |
| combined `wsi_biology` R1 residualised floor (max of arms) | **≈ 5×** — i.e. the headline roughly **doubles** | 3.295× – 12× |

Reasoning, so the prediction is checkable rather than a hedge: `programme_free` is the arm that failed
to complete 40 epochs uncollapsed in 0 of 3 seeds before the momentum fix; its step-200 tripwire rank
spans **6.05×** across five seeds against `programme_only`'s **1.18×** (§4.3); and on the probe block
the collapsed arm carried the floor by 1.2–1.8×. Against that, its *seed*-varied exported
`wsi_biology` fold is **2.099×** against the stable arm's **2.643×** (§4.1a rows 14–15), which cuts the
other way, and is why the interval is wide rather than narrow.

---

## 3. What each outcome MEANS — fixed now so it cannot be chosen afterwards

The paper's central claim is in its own title: **effective rank resolves its own re-measurement on one
co-trained view and not on another.** Concretely, §4.1a rows 26–28: against each view's own floor,
**0 of 6** between-arm pairs are resolvable on `wsi_biology` (floor 3.295×), **6 of 6** on
`rna_biology` (floor 1.019×) and **6 of 6** on `full_biology` (floor 1.020×).

**So the `wsi_biology` half of the claim cannot be hurt by this measurement and the `rna_biology` /
`full_biology` half is the entire exposure.** That asymmetry is stated here, before the numbers, because
it is the opposite of what "the headline floor may double" sounds like.

The six per-pair ratios are already on disk (`~/e0_run/P2_ROBUSTNESS.json`, the twelve D2/D1
artifacts), so the count that each possible floor produces is arithmetic and is tabulated **now**:

| view | the six pair ratios, sorted | resolvable count as a function of the floor `f` |
|---|---|---|
| `wsi_biology` | 1.004, 1.186, 1.573, 1.738, 2.190, **3.246** | `f ≥ 3.246` → **0/6**. Already 0/6 at 3.295×; **any larger floor leaves it 0/6.** |
| `rna_biology` | **1.116**, 1.205, 1.238, 1.766, 2.852, **3.014** | ≤1.116 → 6/6; ≤1.205 → 5/6; ≤1.238 → 4/6; ≤1.766 → 3/6; ≤2.852 → 2/6; ≤3.014 → 1/6; >3.014 → **0/6** |
| `full_biology` | **1.042**, 1.140, 1.234, 2.248, 3.606, **5.250** | ≤1.042 → 6/6; ≤1.140 → 5/6; ≤1.234 → 4/6; ≤2.248 → 3/6; ≤3.606 → 2/6; ≤5.250 → 1/6; >5.250 → **0/6** |

### The four readings, fixed in advance

**(A) `rna_biology` floor ≤ 1.116× AND `full_biology` ≤ 1.042×.**
The claim **survives untouched** and is strengthened on the `wsi_biology` side: 6/6 and 6/6 stand, the
view-conditionality stands, and the only edits needed are numerical (3.295× → the combined value
wherever it is quoted, plus the "3.2× between views" ratio recomputed). *The paper's central claim
survives.*

**(B) `rna_biology` floor in (1.116×, 1.766×] or `full_biology` in (1.042×, 2.248×].**
The claim **survives in substance and its arithmetic changes**: the "usable" view still resolves the
majority of pairs, but 6/6 becomes 5/6, 4/6 or 3/6, and every sentence in the draft that says
"**6 of 6**", "twelve resolvable comparisons" or "separates twelve resolvable comparisons from six
unresolvable ones" must be rewritten to the measured count. *The central claim survives; the counts do
not.* This is the outcome I expect.

**(C) `rna_biology` floor > 3.014× or `full_biology` > 5.250×** (either one).
The view on which rank was said to be usable resolves **nothing**. **The title claim dies**: there is
no view on which rank resolves its own re-measurement, and the paper is not "view-conditional
reproducibility", it is "rank does not resolve its own re-measurement anywhere we looked". That is a
different paper, and the draft needs a real rewrite of §1, §4.5(c), §4.1b, the abstract and the title —
not a number swap. **I will say so plainly and will not soften it into "the counts changed".**

**(D) Any outcome in which the `programme_free` `wsi_biology` floor comes back *smaller* than
3.295×.** Then the stable/unstable framing of the probe block does **not** generalise to the exported
block, the audit's open item 3 is answered in the paper's favour, and the correct write-up is that the
"~2× flattering" inference from the probe was **wrong on this block** — reported as a failed
prediction of mine, with the probe's own numbers left standing. **`programme_only`'s 3.295× would then
remain the arm-conservative reading and 3.295× stays.**

### One row that flips against the paper regardless of view

§4.1a **row 13** (§4.2's within-arm seed fold for D2-I, **3.747×**) is currently the *only*
`wsi_biology` row recorded as **clearing** the exported floor. If the combined floor exceeds 3.747× it
stops clearing, and §4.2's "the seed alone moves rank further than five same-seed retrains do" is
**withdrawn**. That withdrawal is committed to here, in advance, because it is the sentence a wider
floor most tempts one to keep.

Similarly, §4.1a row 6 (D1-B seed 43, 3.246× against 3.295×, *"the closest of the seven, clearing by
1.5% in the wrong direction"*) becomes less close, and the *"seven of seven fall inside a floor of
3.295×"* quotation becomes seven of seven inside a larger floor — a change of number, not of verdict.

### And one that would be a genuine problem for the instrument, not for the claim

If `programme_free`'s five repeats are **not all trainable** — i.e. one or more trips the rank
tripwire, fails a gate, or does not reach 40 epochs — then the five are not five repeats of one
configuration and **no floor may be computed from the survivors**. Selecting the runs that completed
is exactly the selection effect that would flatter the floor. In that case the outcome reported is
"the unstable arm does not admit a same-seed floor at 40 epochs because it does not reliably reach 40
epochs", the completion count is reported as the result, and **3.295× stays in the paper with that
sentence attached.** This is stated in advance because it is a plausible outcome and it is not a null
result.

---

## 4. What this measurement cannot do

- It does not vary the seed; §4.2 owns that axis and finds it larger.
- It does not touch a second architecture, cohort or modality. `no_external_cohort` stays undischarged.
- It says nothing about the probe block, the in-run training batch, the 16-patient gate batch or the
  282-patient live checkpoint. Those remain in `absent_blocks`.
- n = 5 per arm is n = 5. No interval will be quoted, and "the floor is 5×" will not be written as
  "rank varies 5×".

## 5. Files

- To be produced: `~/e0_run/d1_envelope_pf/rep{1..5}/`, `rep{1..5}.npz`,
  `~/e0_run/d1_envelope_pf/out/P2_ENVELOPE_FLOORS_PF.json`,
  `.../P2_ENVELOPE_FLOORS_BOTH_ARMS.json`, `d1_envelope_pf_readout.log`
- Reused unchanged: `v2/research/rebase/p2/p2_envelope_floors.py`,
  `v2/research/rebase/d1_envelope_readout.py`, `~/chain_retrain_envelope.sh` (arm flag only)
- Sources: `NOTEBOOK_ENTRIES/PREDECLARED_retraining_envelope_20260804T0330Z.md`;
  `retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §1;
  `retraining_floors_for_every_statistic_view_and_block_20260804T1220Z.md`;
  `the_probe_block_has_a_floor_at_last_20260804T1620Z.md` §3;
  `p2_completeness_pass_stale_open_states_and_the_figure_that_drew_the_wrong_floor_20260804T2030Z.md` §7 item 3
- **Not touched:** `claim_guards.py`, `claim_evidence.json`, any other agent's `PREDECLARED_*`,
  `paper/P2_RANK_DRAFT.md`.
