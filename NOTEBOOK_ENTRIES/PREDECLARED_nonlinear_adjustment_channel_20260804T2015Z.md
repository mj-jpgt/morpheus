# PREDECLARED — does the residual nonlinear confound structure explain P1's channel?

**Written:** 2026-08-04 20:15 UTC (box clock, `date -u` on 150.136.45.194 = `2026-08-04T20:08:44Z`;
the local repo clock agrees). Note that several `NOTEBOOK_ENTRIES/` stamps run ahead of the wall
clock — `tcga_nonlinear_confound_probe_result_20260804T2100Z.md` and
`direction_matched_floor_20260805T0030Z.md` are both already committed. **This file is written after
both**, at commit `ff9a6f9`, and **before any number in it exists**.

**Question.** `NOTEBOOK_ENTRIES/tcga_nonlinear_confound_probe_result_20260804T2100Z.md` (`7482a38`)
established that after the exact `cancer + pooled TSS` cross-fitted ridge adjustment CALIBRA applies
before every channel number in P1, an out-of-fold k-NN still recovers tissue source site at 4.3–4.9×
chance and cancer type at 3.4–4.9× chance, at the permutation floor. The residual is real. **Is it
large enough to produce the channel?**

---

## 1. What "the channel" is here, stated exactly, because two statistics carry the name

The brief describes P1's headline as "a morphology→molecular channel of 0.58–0.62 held-out against a
permutation null of ~0.147". Read against the source that is **P1 §4.4**, and the two halves come
from two different places, so both are pinned here:

| | statistic | value quoted in P1 | null |
|---|---|---|---|
| **S1 (headline)** | `spectral.top_canonical_correlation(adjust(x), adjust(y), n_components=16)` — an in-sample maximum over 16 whitened directions per side, computed on the **held-out patient partition** (`--partition test`, n = 2,766) | §4.4: `d2_h::wsi_biology` **0.6052**, `d2_i::wsi_biology` **0.4703** | §4.4: within-cancer pairing permutation, median **0.1483 / 0.1472**, p = 0.0005 (1/2001) |
| **S2** | `spectral.heldout_top_cca(adjust(x), adjust(y), n_components=16, seed=42)` — directions fitted on one half of the held-out patients, correlation scored on the other | `wsi_biology` 0.47–0.54 | its own destroyed-pairing null 0.06–0.09 (`calibration.py:74`) |

"0.58–0.62" is the **S1** family (§4.4's 0.6052; §3.10's 0.5861/0.6214; §4.12's 0.5970–0.6126) and
"0.147" is **S1's** null. **Pairing an S2 value with 0.147 would mix two nulls**, which
`paper/P2_RANK_DRAFT.md:1049-1058` already warns about for 0.140 vs 0.147. **Both S1 and S2 are
reported for every arm below**, and every number carries its statistic name and its block.

**Anchor block, fixed now.** `runs/d2_final/artifacts/d2_h_seed42.npz` and `d2_i_seed42.npz` on
persistent NFS — the copies the probe entry §6 showed reproduce §4.2 to four decimals, in preference
to `~/e0_run/d2_v3/*`; state `wsi_biology`; `--partition test`; targets
`morpheus_phase_d/data/frozen_rna_targets.npz` with the `RANDOM_CONTROL__` columns excluded exactly
as `run_calibra.py:364-371` does; design `confound_design({cancer, pooled tss(min_site_count=10)},
["cancer","tss"])`; `n_components=16`; `seed=42`; `KFold(5, shuffle=True, random_state=42)`.

**Gate before anything else is believed:** arm A1 below must reproduce §4.4's `0.6052 / 0.4703` and
null median `0.1483 / 0.1472` on `wsi_biology`. If it does not, the artifact provenance defect of the
probe entry §6 has bitten again and **no verdict is issued** — the discrepancy is reported instead.

---

## 2. A structural fact that must be stated before the arms, because it constrains the answer

`confound_design` on `["cancer", "tss"]` produces a **purely one-hot, purely additive** design: 21
cancer indicators + 84 pooled-site indicators + `OTHER` + two `dummy_na` columns, ≈ 108 columns. Two
consequences:

1. For a **categorical** predictor set, the linear model on a one-hot design already spans the
   conditional-mean model. The only thing a nonlinear learner can add on top of an *additive* one-hot
   design is the **cancer × site interaction**.
2. The probe entry's `nesting_diagnostic` measured that **site nests totally inside cancer on the test
   partition** — 0 of 84 kept sites contributes patients to two cancers, only the pooled `OTHER` class
   spans them. So the site indicator already *determines* the cell for all but the `OTHER` rows, and
   the cancer × site interaction is nearly degenerate by construction.

**Therefore a "nonlinear adjustment against the same design", read literally as a nonlinear estimate of
E[X | cancer, site], is predicted in advance to be near a no-op**, and if it is, step 3 will show the
k-NN recovery unmoved and the exercise will have tested nothing. That prediction is written here so
that it is a prediction and not a post-hoc excuse. It is why the arm list below runs **past** the
literal request: A2 is the exact upper limit of any first-moment adjustment on this design, and A5
removes structure no residualiser removes at all.

---

## 3. The arms (all cross-fitted, all applied identically to X and to Y, all `KFold(5, seed=42)`)

| arm | adjustment | why |
|---|---|---|
| **A0** | none (column-centred only) | the raw block, for the ×-of-raw denominators |
| **A1** | `residualise.cross_fitted_residuals`, Ridge α=1.0 on the 108-column additive one-hot | **the incumbent.** Every channel number in P1 is measured here. Reproduction gate. |
| **A2** | cross-fitted ridge, α = 1e-6, on the **saturated cell** one-hot (one column per observed (cancer, site) cell) | the **exact limit** of first-moment adjustment on a categorical design. If A2 ≈ A1 and A2 leaves the k-NN at 4.8×, the surviving structure is provably *not* a first-moment nonlinearity, and no kernel/forest/boosting estimate of E[X\|C] can do better. |
| **A3** | cross-fitted **RBF kernel ridge** on the 108-column one-hot | the primary nonlinear arm. Chosen over gradient boosting and random forest for three stated reasons: (i) on a one-hot design the RBF kernel is a monotone function of Hamming distance between (cancer, site) pairs, so its hypothesis space is exactly "smoothed cell means with arbitrary interaction" — a strict superset of A1's span; (ii) it is multi-output in **one** solve, so all 256 image axes and all molecular targets come from a single factorisation; (iii) crucially, the kernel and its per-fold Cholesky factors depend on the **design only**, never on Y, so they are computed once and the permutation null runs at the project's own **2,000** permutations. A reduced permutation count was the one declared item that defeated the previous run (its forest null was cut to 50); this choice removes that failure mode rather than repeating it. |
| **A4** | cross-fitted **multi-output random forest**, 300 trees, on the 108-column one-hot | the second family, which fails differently: axis-aligned thresholds on indicator columns, arbitrary depth-driven interaction, no metric assumption. Multi-output so one fit serves every column. Its null is affordable only at **200** permutations (p floor 0.0050), declared now. |
| **A5** | cross-fitted **conditional location–scale**: subtract the per-cell mean and divide by the per-cell per-axis SD, SD shrunk toward the pooled SD at prior count m = 10, cells = observed (cancer, site) | **This is NOT "the same adjustment made nonlinear" and is not offered as one.** It removes second-moment structure, which no residualiser of any flexibility removes, because a residual is by construction only a first-moment operation. It is included because step 3 requires that *something* actually reduce the k-NN recovery, or the whole exercise is circular. Its result is reported under a different heading and with the over-removal guard of §6. |

Every arm's hyperparameters are fixed here. A3 is run over the grid α ∈ {1.0, 0.1}, γ ∈ {0.25, 0.5,
1.0} (squared Hamming distance between two distinct one-hot confound rows is 2 or 4, so this spans
kernel widths from near-identity to near-constant). **The quoted A3 cell is the one that minimises the
step-3 k-NN recovery — the adjustment-validity criterion — never the one that minimises the channel.**
The whole grid is tabulated so the selection can be checked.

---

## 4. The three readings, with bands, fixed before any number exists

Grading is on **retention of excess over the arm's own null median**:

```
excess(arm) = S1(arm) - null_median(arm)          # same statistic, same permutation procedure
retention   = excess(arm) / excess(A1)
```

Grading on `excess` rather than on `S1` is deliberate and is fixed here: a stronger adjustment can
move the capacity floor as well as the signal, and comparing a raw 0.60 against a null measured under
a different adjustment is exactly the error §4.4 of P1 warns about.

* **Reading 1 — CHANNEL SURVIVES.** Under the strongest arm that **passes the step-3 validity check of
  §5**, `retention ≥ 0.80` on `wsi_biology`, on both `d2_h` and `d2_i`, with S2 agreeing in direction.
  Conclusion: the residual is real but is not what the channel reads; P1's claim stands with
  *"removed from the first moment"* replacing *"gone"*.
* **Reading 2 — CHANNEL COLLAPSES.** `retention ≤ 0.20` under a step-3-valid arm, **and** the
  over-removal guard of §6 passes. Conclusion: the headline was substantially confound. Reported as a
  project-level finding.
* **Reading 3 — PARTIAL.** `0.20 < retention < 0.80`. **Report the fraction, do not adjudicate.**
* **Reading 4 — NO TEST.** No arm passes step 3. Then this run has not tested the question, and the
  honest output is (a) that fact, (b) the structural reason of §2 measured rather than asserted, and
  (c) A5's numbers reported separately as a *different* experiment. **Reading 4 is not a favourable
  result and will not be written as one.**

Disagreement between arms **resolves downward for the channel** (the arm showing the largest channel
loss is quoted as the finding), mirroring the probe entry's "families resolve upward" rule for the
confound: in both cases the rule points away from the comfortable answer.

---

## 5. Step 3 — the validity check that stops this being circular

For every arm, re-run `nonlinear_confound_probe.probe_state` on the arm-adjusted `wsi_biology`, test
partition, targets `site` (85 pooled classes, design chance 0.011765) and `cancer` (21 classes,
0.047619), k ∈ {1,3,5,10,15,25,50}, plain **and** inverse-training-frequency-weighted vote, reading =
**max over k and over both vote rules** (as in the probe entry), **200 global permutations**, p floor
1/201 = 0.0050.

**Which null.** The **global** null, and the reason is the probe entry's: TCGA site nests inside
cancer, so a within-cancer permutation of *site labels* returns them essentially unchanged and its
oracle rate is 0.1739 against a design chance of 0.0118. The within-cancer label null is reported
beside the site rows for decomposition only. Chance is **measured** by the global null, not assumed.

*(The channel's own null in §4 is a different object and the degeneracy does not apply to it: it
permutes patient **pairing** between two blocks, not a label within a stratum, so within-cancer is
non-degenerate there and is the project's convention. Both the within-cancer and the global pairing
null are computed for the channel; the within-cancer one is quoted.)*

Bar, fixed now, against A1's measured multiples (site 4.80×, cancer 4.67× on `d2_h`, from the probe
entry's `adjusted, standardised` rows — recomputed here on the identical rows rather than quoted):

* **verified** — arm's k-NN max ≤ **2.0× chance** for **both** targets;
* **partially verified** — reduced by ≥ **50%** relative to A1's multiple on both targets, but above 2.0×;
* **not verified** — otherwise. **A "not verified" arm cannot support Reading 1.**

---

## 6. Over-removal guard, for the unfavourable direction

A collapse under an arm that also destroys real signal is not evidence of confound. If any arm reaches
Reading 2 or the lower half of Reading 3, its **attenuation is measured through the identical
adjustment** with `calibration.spike_recovery_curve` (planted spikes of known strength pushed through
the same adjuster; `attenuation_slope`, `transmission_floor`, `detection_floor` read by
`calibration.floors_from_recovery`, the instrument's own rule). An arm whose attenuation slope is far
below A1's has over-adjusted and its collapse is **not** reported as a confound finding. If that
measurement does not run, any collapse is reported as **ambiguous between confound and
over-adjustment**, explicitly.

---

## 7. Step 2 — the labels-only ceiling (the bound from the other side)

How much channel is obtainable from `cancer + pooled TSS` **alone**, with no image features. Three
readouts, all with the same 16-component budget so capacity is matched:

* **C1, raw ceiling.** `top_canonical_correlation(D, y_raw, n_components=16)` and the held-out
  `heldout_top_cca_indexed` equivalent, for `D` = the saturated cell one-hot **and** the additive
  108-column design. This is the largest channel confound labels can carry before any adjustment.
* **C2, pipeline ceiling.** Push `D_cell` through the identical pipeline **as if it were the image
  representation** — adjust both sides under each arm, S1 at 16 components, with its own within-cancer
  pairing null. This is what a purely-confound "representation" scores after CALIBRA's adjustment.
* **C3, variance accounting.** Cross-fitted R² of the molecular block on the labels, and of
  `wsi_biology` on the labels, per block and pooled. A direct statement of how much of each side is
  confound.

C1/C2 are the ceiling the brief asks for. They are cheap and are run regardless of the outcome of §4.

---

## 8. Eight ways a FAVOURABLE result (Reading 1) would be untrustworthy — all reported either way

Two favourable results landed today already; the standard has to hold in this direction too.

1. **Step 3 fails.** If the arm did not reduce the k-NN recovery, "the channel survives adjustment" is
   the linear adjustment under a new name and says nothing. Reading 4, not Reading 1.
2. **The adjuster is degenerate on this design.** Per-axis correlation between the arm's residuals and
   A1's residuals is reported for every arm. At ≥ 0.99 the arm is a relabelled A1 and is reported as
   such, whatever its channel reads.
3. **The null moved.** If a nonlinear arm raises its own null median, a flat `S1` is not a surviving
   channel. Grading is on `excess` over the arm's **own** null for exactly this reason, and both the
   raw `S1` and the null are tabulated so the reader can check the arithmetic.
4. **The ceiling is large.** If §7's labels-only ceiling is comparable to the observed channel, a
   surviving `S1` is not evidence of a morphology→molecular channel. The ceiling is printed beside the
   channel in the headline table, not in an appendix.
5. **Capacity.** `S1` is a maximum over 16 whitened directions at n = 2,766. If an arm reduces the
   effective rank of the adjusted block, `S1` can *rise* for capacity reasons alone.
   `spectral.effective_rank` (CANONICAL variant) of every adjusted block is reported.
6. **S1 is in-sample across directions.** If the channel survives on S1 but not on S2 (fitted
   directions, held-out rows), that is an in-sample-maximum artefact. Both are reported for every arm;
   a disagreement between them blocks Reading 1.
7. **Artifact disagreement.** If `d2_h` and `d2_i` disagree in direction, no verdict is issued.
8. **Folds.** The adjuster's `KFold(5, seed=42)` and the probe's `_stratified_folds(seed=42)` are
   different partitions, so a probe test row can be an adjuster train row — the same non-independence
   the probe entry recorded at its §7 item 10. It is recorded here, not waved away; it is inside the
   permutation nulls as well as inside the observed, because every null is recomputed through the same
   adjuster.

## 8b. Three ways an UNFAVOURABLE result (Reading 2) would be untrustworthy

1. **Over-removal**, handled by §6.
2. **A5 is not the same adjustment.** If the collapse appears only under A5, the finding is "a
   second-moment adjustment removes the channel", which is a different and weaker sentence than "the
   nonlinear adjustment removes the channel". It will be written with that distinction intact.
3. **Loss of rank.** If an arm drives the adjusted block's effective rank near the 16-component
   budget, `S1` falls for a capacity reason. Item 5 above covers it; it cuts both ways.

---

## 9. What will not be run, declared in advance

* Cohorts other than TCGA; partitions other than `test`; `min_site_count` other than 10.
* States other than `wsi_biology` for any verdict (`full_biology` and `rna_biology` are reported for
  breadth only, and both are RNA-derived, so their channel to RNA targets is near-circular at ~0.89
  and is not a morphology→molecular measurement).
* Gradient boosting, for the cost reason given in §3 (single-output, so 5 folds × every column per
  permutation). A3 and A4 are the two families run.
* Kernel/forest hyperparameters beyond the grid fixed in §3. **No permutation count will be increased
  after seeing a result**, and none of the counts declared here will be quietly cut without saying so.

---

## 10. Where the code and the numbers go

New module `v2/calibra/nonlinear_adjustment.py` with tests in
`v2/tests/test_nonlinear_adjustment.py`. **No statistic is computed inline**: `top_canonical_correlation`,
`heldout_top_cca`, `effective_rank`, `floors_from_recovery`, `spike_recovery_curve`,
`nonlinear_confound_probe.probe_state`, `global_permutations`, `null_summary` are imported. Two
equivalence tests are the contract that keeps this comparable to the published run:

* the generic adjuster at `model="ridge"` returns **exactly** `cross_fitted_residuals`' output;
* the generic channel measurement under that adjuster reproduces `calibration.permutation_null`'s
  dictionary **exactly**.

Outputs under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/nonlinear_adjustment/`.
Workspace deployed by `git -c core.autocrlf=false archive HEAD` and verified per file against
`git ls-tree -r HEAD` blob SHA-1. Threads capped
(`OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1`), process parallelism, CPU only.
