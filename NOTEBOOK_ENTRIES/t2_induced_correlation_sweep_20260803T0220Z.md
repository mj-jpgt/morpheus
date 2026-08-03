## 2026-08-03 02:20 UTC — Track 2: the induced correlation is structural, not degrees-of-freedom. It is flat in n over a 13× range, flat in design rank over a 15× range, and 20× larger than a matched-rank structureless design

**Logged:** 2026-08-03 02:20 UTC. **How obtained:** `python -m morpheus.v2.calibra.induced_correlation_sweep --artifact d2_h_seed42.npz --state wsi_biology --partition all --n-grid 500,1000,2000,2530,4000,6427 --seeds 42,43,44 --n-draws 40 --n-jobs 8 --clinical-covariates tcga_clinical_covariates.parquet`, 270 cells, Lambda box `~/ws_p1`; graded by `p1_evidence/grade_t2.py` against the predictions in `v2/research/rebase/nature/P1_PREDECLARATION.md`, which was written and hashed before this run.

### Technical

**The closed form is an identity, and it is verified draw for draw.** Across all 270 cells the maximum
absolute disagreement between the pipeline's measured induced correlation and the Yule (1907) /
Frisch–Waugh–Lovell partial-correlation formula, evaluated on the *same* planted (u, v) pairs, is
**8.6 × 10⁻¹⁶** (median 3.1 × 10⁻¹⁶); Pearson correlation between the two, per draw, is
**1.000000**. There is no empirical question about the mechanism. The only open questions are how its
two factors scale, and whether the magnitude is larger than the classical degrees-of-freedom term.

**P5 — the falsifier that decides whether Track 2 has content. PASSES.**
Matched design rank, matched n, structure removed:

| n | k_eff | real cancer+TSS | rows-permuted, same design | Gaussian, k=99 | real / permuted |
|---:|---:|---:|---:|---:|---:|
| 500 | 32.7 | 0.0866 | 0.0087 | 0.0164 | 9.9× |
| 1,000 | 45.7 | 0.0804 | 0.0065 | 0.0101 | 12.4× |
| 2,000 | 87.4 | 0.0809 | 0.0040 | 0.0057 | 20.2× |
| **2,530** | **104.3** | **0.0748** | **0.0037** | **0.0035** | **20.4×** |
| 4,000 | 144.4 | 0.0789 | 0.0031 | 0.0026 | 25.7× |
| 6,427 | 215.1 | 0.0718 | 0.0020 | 0.0016 | 35.4× |

Predeclared bar: structureless arms ≤ 0.025 at n = 2,530 and ≥ 3× smaller than real. Measured 0.0037
and 0.0035, ratio 20.4×. And the separation *widens* with n, because the structureless arms decay like
a sampling term while the real design does not.

**P4 — bias, not sampling fluctuation. PASSES.** At the anchor design, `|r|(6,427)/|r|(2,530) = 0.960`
against a predeclared bar of ≥ 0.60 and a pure-sampling expectation of √(2530/6427) = 0.627. Across a
13× change in n the induced correlation moves from 0.0866 to 0.0718 — a 17% decline, not the 72%
decline sampling noise would give.

**P0 — the plan's own guess, `|r| ~ k/n`. FALSIFIED.** Fitted log–log exponents on real designs:
`b_k_eff = +0.288` (predicted +1.0), `b_n = −0.180` (predicted −1.0), R² = 0.398.
The module's own derived Eq. (3) prediction of `b_k = −0.5, b_n = 0` is also wrong on the rank axis.

**The rank ladder is what explains both.** At n = 2,530, holding the cohort fixed and moving design
rank over 15×:

| design | k | k_eff | k_eff_shared | R²_x | R²_y | measured | P1 | P2 | P3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| none | 0 | 0 | 0 | 0 | 0 | 0.0000 | — | — | — |
| tss_pool50 | 2 | 0.0 | 0.0 | 0.00 | 0.00 | 0.0003 | — | — | — |
| cancer | 33 | 31.0 | 9.04 | 0.297 | 0.473 | **0.0844** | 0.0706 | 0.0453 | 0.0859 |
| cancer+tss_pool50 | 35 | 31.0 | 9.05 | 0.297 | 0.474 | 0.0844 | 0.0703 | 0.0448 | 0.0860 |
| tss_pool10 (no cancer) | 76 | 74.0 | 9.73 | 0.166 | 0.245 | **0.0379** | 0.0190 | 0.0162 | 0.0492 |
| **cancer+tss_pool10 (anchor)** | 109 | 104.3 | 10.26 | 0.321 | 0.470 | **0.0748** | 0.0407 | 0.0258 | 0.0896 |
| cancer+tss+dx_year | 111 | 106.3 | 10.27 | 0.321 | 0.470 | 0.0744 | 0.0402 | 0.0255 | 0.0897 |
| cancer+tss+dx_year+age | 113 | 108.3 | 10.29 | 0.320 | 0.469 | 0.0732 | 0.0397 | 0.0252 | 0.0897 |
| cancer+tss+dx_year_cat | 141 | 132.3 | 10.45 | 0.313 | 0.464 | 0.0743 | 0.0354 | 0.0227 | 0.0900 |
| cancer+tss_pool3 | 283 | 269.5 | 11.28 | 0.309 | 0.450 | 0.0729 | 0.0232 | 0.0154 | 0.0924 |
| cancer+tss_pool1 | 511 | 450.9 | 11.85 | 0.305 | 0.445 | 0.0727 | 0.0176 | 0.0118 | 0.0936 |

Three readings, all reportable:

1. **Design rank is almost irrelevant.** From k_eff = 31 to k_eff = 451 — a 14.5× increase, 478 extra
   columns — the induced correlation moves 0.0844 → 0.0727, a **−14% change**. Widening a nuisance
   model does not manufacture more induced correlation.
2. **Cancer type alone produces essentially the whole effect** (0.0844 of an anchor 0.0748 — the
   anchor is in fact *slightly lower*). The 478 tissue-source-site columns add nothing, because they
   add rank without adding explanatory power over *both* modalities.
3. **What does track it is how much the design explains of BOTH sides.** `tss_pool10` without cancer
   has *more* rank than `cancer` (74 vs 31) but half the induced correlation (0.0379 vs 0.0844),
   because its R² over the two modalities is 0.166/0.245 instead of 0.297/0.473. `k_eff_shared` — the
   effective dimension of the jointly explained subspace — moves only 9.0 → 11.9 across the entire
   ladder, which is why the effect is flat.

**Predictor accuracy across all 270 cells** (ratio measured / predicted; a perfect predictor is 1.0):

| predictor | cells | median ratio | p10–p90 | log₁₀ RMS error |
|---|---:|---:|---|---:|
| P1 `0.6745·κ/√k_eff` | 167 | 1.98 | 1.15 – 3.98 | 0.365 |
| P2 `0.6745·R_s R_a/√k_eff` | 167 | 2.98 | 1.65 – 6.41 | 0.526 |
| **P3 `0.6745·R_x R_y/√k_eff_shared`** (real designs only) | 169 | **0.886** | **0.76 – 1.07** | **0.079** |
| P3, all cells incl. structureless | 241 | 0.839 | 0.39 – 1.05 | 0.204 |

P1 and P2 under-predict by 2–3× and drift. P3 predicts within roughly ±20% across the whole grid from
the (design, X, Y) spectra alone — no draws, no spike, no recovery curve. **P3 remains post hoc on the
rank axis** (it was written after P1/P2 failed at n = 2,530, seed 42, and the predeclaration says so);
what is genuinely out of sample here is every other n and every other seed, and it holds there.

**Cross-fitting suppresses the classical term, which is why the structural term is visible.** For the
structureless arms the classical N−R sampling scale is `0.6745/√(n−R)`. Measured / classical across
24 structureless cells: **median 0.379** (p10 0.238, p90 0.699). Cross-fitted residualisation delivers
only about a third of the in-sample degrees-of-freedom inflation that Winkler et al. (2020) and the
Muirhead/Anderson N−R result describe. The one place it cannot help is when the design nearly spans
the sample: `gaussian_k600` at n = 500 (k_eff = 499) gives 0.055, and the classical formula there
predicts 0.617, so even suppressed it is the largest structureless value in the grid.

`baseline_is_null_like` is True for 100% of cells in all three design modes, so nothing here is the
readout defect the gate exists to catch.

### In plain terms

When you remove the same background variables from two unrelated measurements, you accidentally make
the leftovers look related. We already knew the formula for this — it is a 1907 identity, and our
pipeline reproduces it to fifteen decimal places. The question was whether the size we see, about 0.08,
is just the well-known statistical bookkeeping cost of using up degrees of freedom, or something else.

It is something else. If you take the exact same nuisance model and shuffle its rows so it still has
the same size but no longer has anything to do with the patients, the effect drops from 0.075 to
0.0037 — twenty times smaller. If you make the nuisance model fifteen times bigger, the effect barely
moves. If you collect thirteen times more patients, the effect barely moves. What does move it is how
much of *both* the image side and the molecular side the nuisance model can explain — and for us that
is almost entirely cancer type, not the 478 hospital columns.

So the practical warning is not "keep your nuisance model small". It is "the induced correlation is
about the size of the thing you are adjusting for, not the number of columns you spend on it — and
more data will not make it go away."

### Meaning for the claim

* **Track 2's magnitude framing survives.** Combined with the Winkler verdict (no magnitude published),
  the reportable claim is: *the induced correlation under correctly cross-fitted residualisation on
  real cross-modal data is ~0.07–0.09, is a bias rather than a sampling fluctuation, and is ~20×
  larger than the residual degrees-of-freedom term that the classical N−R result accounts for.*
* **The predicted-floor claim is now portable.** P3 lets another cohort compute its expected induced
  correlation from its own design and data spectra before running the instrument. That is the
  certificate field P4 needs and the sizing rule Track 3 needs.
* **The ledger row "confound adjustment does not destroy signal (attenuation 0.94–1.23)" is hardened.**
  Its stated falsifier was "attenuation far from 1 under a differently constructed confound design at
  comparable rank"; eleven designs spanning 15× in rank all sit at the same induced baseline and every
  cell passes the null-like gate.
* **A correction to our own earlier reasoning must be published with it.** The plan predicted `k/n` and
  the module's first two derivations predicted a `1/√k_eff` decay. All three are wrong. Design rank is
  the wrong axis; jointly-explained variance is the right one.
* **P3 (the biology paper) inherits a caution:** D3 asks whether the channel survives purity entering
  the adjustment set. On this evidence, adding columns to the design changes the induced baseline by
  ~1%, so any channel change D3 sees when purity is added is attributable to purity, not to rank. That
  is the disambiguation T2.2 was for, and it comes out in D3's favour.
* **Not done:** no TCGA purity table exists on either machine, so the purity rank point named in the
  plan was replaced by `dx_year`/`age` and by the TSS pooling threshold. Recorded, not silently
  substituted — the module docstring and the sweep protocol both carry it.

### Files / commits

`v2/calibra/induced_correlation_sweep.py` (inherited, audited), `p1_evidence/grade_t2.py`.
Data: `p1_evidence/track2/main_rows.csv` (270 cells), `main_law.json`, `graded_predictions.json`, under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/`.
Predeclaration: `v2/research/rebase/nature/P1_PREDECLARATION.md`, commit 1c4b4b5.
Estimator-robustness (P6) and the floor sweep (T2.6) are still running.
