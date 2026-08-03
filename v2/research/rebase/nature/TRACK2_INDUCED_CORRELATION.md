# Track 2 — induced correlation under shared residualisation

Phase-gate item **C** for P1. Every prediction below was written into
`v2/research/rebase/nature/P1_PREDECLARATION.md` and committed (1c4b4b5) **before** the sweep ran.

---

## 0. The prior-art question, settled first — because it could have voided the whole track

**Winkler AM, Renaud O, Smith SM, Nichols TE. "Permutation inference for canonical correlation
analysis." *NeuroImage* 2020;220:117065.** Full text retrieved from PMC7573815 (open access), five
targeted passes over Theory §2.6, Simulations, Results §4.3/4.5, Discussion §5.3, every figure caption
and every table header.

**VERDICT: NO. It does not report an induced-correlation magnitude.**

Everything it reports is in error-rate, power or p-value units: Table 6 per-comparison error rate (%),
Table 7 pcer % (simple residualisation 83.85% [82.17–85.40] vs ~5% for Huh–Jhun/Theil), Table 8 pcer %,
Table 9 FWER % against N = 100…1000, Table 10 power %. Figure 2 plots canonical correlation on the
*horizontal* axis and p-value on the vertical, for a scenario with no nuisance variables. **No figure or
table anywhere reports a correlation value obtained from null data after residualisation, and no
numeric correlation coefficient appears in the narrative at all.**

The closest statement, §2.6 verbatim:

> "While **Y** occupies an N-dimensional space, **Ỹ** occupies a smaller one; its dimensions are, at
> most, of a size given by the rank of R_Z, which is N−R…"

> "With fewer effective observations determined by this lower space after residualisation, and the same
> number of variables, the sample canonical correlations in the unpermuted case are **stochastically
> larger** than in the permuted, which in turn leads to an excess of spuriously small p-values."

A stochastic-ordering statement, relative to the permutation distribution, with no magnitude. Discussion
§5.3 confirms the framing: "…leads to inflated error rates and an invalid test… particularly if the
number of nuisance variables is relatively large compared to the sample size". The R/N dependence is
asserted qualitatively and shown only through error rates. The only place R enters a formula is the
Bartlett/Wilks df correction (`C = R` for partial CCA) — used to *undo* the effect for a parametric
p-value, never inverted into a magnitude. Winkler et al. 2014 (PMC4010955) takes the same posture.

**What we must therefore cite, so that no reviewer can present it as an omission:** (a) Yule (1907) /
Frisch–Waugh–Lovell for the identity; (b) Winkler et al. 2020 for the inferential consequence and the
fix; (c) **Muirhead (1982) / Anderson (2003)** for the classical result that rank-R residualisation
leaves an effective sample size of N−R, which Winkler's own correction uses.

**The claim that remains ours** is the magnitude *in correlation units*, on real cross-modal data,
under **cross-fitted** residualisation, together with the demonstration below that it is structural
rather than the classical degrees-of-freedom term. Nothing broader.

## 1. The mechanism is an identity, and it is verified draw for draw

For `corr(u, v) = 0` the residualised correlation is the multivariate partial-correlation identity

```
r_induced  =  − R_s R_a cos(θ) / sqrt((1 − R_s²)(1 − R_a²))
```

with `R_s², R_a²` the design's cross-fitted R² for the two scores and `cos(θ) = corr(ŝ, â)`.

Across **all 270 sweep cells**, the maximum absolute disagreement between the pipeline's measured
value and this formula evaluated on the *same* planted (u, v) pairs is **8.6 × 10⁻¹⁶** (median
3.1 × 10⁻¹⁶); the per-draw Pearson correlation between the two is **1.000000**. There is no empirical
question about the mechanism. The open questions are how its factors scale, and whether the magnitude
exceeds the classical term.

## 2. The falsifier that decides whether Track 2 has content — P5 — PASSES

If a structureless design of matched rank at matched n reproduced the effect, the phenomenon would be
degrees-of-freedom bookkeeping, fully covered by the classical N−R result, and we would have nothing.

| n | k_eff | real cancer+TSS | same design, **rows permuted** | Gaussian, k = 99 | real / permuted |
|---:|---:|---:|---:|---:|---:|
| 500 | 32.7 | 0.0866 | 0.0087 | 0.0164 | 9.9× |
| 1,000 | 45.7 | 0.0804 | 0.0065 | 0.0101 | 12.4× |
| 2,000 | 87.4 | 0.0809 | 0.0040 | 0.0057 | 20.2× |
| **2,530** | **104.3** | **0.0748** | **0.0037** | **0.0035** | **20.4×** |
| 4,000 | 144.4 | 0.0789 | 0.0031 | 0.0026 | 25.7× |
| 6,427 | 215.1 | 0.0718 | 0.0020 | 0.0016 | 35.4× |

Predeclared bar: structureless arms ≤ 0.025 at n = 2,530 **and** ≥ 3× smaller than real. Measured
0.0037 and 0.0035, ratio **20.4×**. The row-permuted design has *identical rank and identical column
marginals* and differs only in having no relationship to any patient.

The separation widens with n because the structureless arms decay like a sampling term while the real
design does not.

## 3. It is a bias, not a sampling fluctuation — P4 PASSES, P0 is FALSIFIED

At the anchor design, across a **13× change in n**:

| n | 500 | 1,000 | 2,000 | 2,530 | 4,000 | 6,427 |
|---|---:|---:|---:|---:|---:|---:|
| induced (median of 3 seeds) | 0.0866 | 0.0804 | 0.0809 | 0.0748 | 0.0789 | 0.0718 |

`|r|(6,427) / |r|(2,530) = 0.960`, against a predeclared bar of ≥ 0.60 and a pure-sampling expectation
of √(2530/6427) = 0.627. A 17% decline over 13× more patients, not the 72% decline a sampling term
would give. **More data does not remove it.**

**P0 — the plan's own guess, `|r| ~ k/n` — is falsified.** Fitted log–log exponents on real designs:
`b_k_eff = +0.288` (predicted +1.0), `b_n = −0.180` (predicted −1.0), R² = 0.398. The module's own
derived Eq. (3) form (`b_k = −0.5`, `b_n = 0`) is also wrong on the rank axis. We were wrong twice and
both are reported.

## 4. Design rank is the wrong axis; jointly-explained variance is the right one

Rank ladder at n = 2,530, cohort held fixed, rank moved over 15×:

| design | k | k_eff | k_eff_shared | R²_x | R²_y | **measured** | P1 | P2 | P3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| none | 0 | 0 | 0 | 0 | 0 | 0.0000 | — | — | — |
| tss_pool50 | 2 | 0.0 | 0.0 | 0.00 | 0.00 | 0.0003 | — | — | — |
| **cancer** | 33 | 31.0 | 9.04 | 0.297 | 0.473 | **0.0844** | 0.0706 | 0.0453 | 0.0859 |
| cancer + tss_pool50 | 35 | 31.0 | 9.05 | 0.297 | 0.474 | 0.0844 | 0.0703 | 0.0448 | 0.0860 |
| **tss_pool10 (no cancer)** | 76 | 74.0 | 9.73 | 0.166 | 0.245 | **0.0379** | 0.0190 | 0.0162 | 0.0492 |
| **cancer + tss_pool10 (anchor)** | 109 | 104.3 | 10.26 | 0.321 | 0.470 | **0.0748** | 0.0407 | 0.0258 | 0.0896 |
| cancer + tss + dx_year | 111 | 106.3 | 10.27 | 0.321 | 0.470 | 0.0744 | 0.0402 | 0.0255 | 0.0897 |
| cancer + tss + dx_year + age | 113 | 108.3 | 10.29 | 0.320 | 0.469 | 0.0732 | 0.0397 | 0.0252 | 0.0897 |
| cancer + tss + dx_year_cat | 141 | 132.3 | 10.45 | 0.313 | 0.464 | 0.0743 | 0.0354 | 0.0227 | 0.0900 |
| cancer + tss_pool3 | 283 | 269.5 | 11.28 | 0.309 | 0.450 | 0.0729 | 0.0232 | 0.0154 | 0.0924 |
| cancer + tss_pool1 | 511 | 450.9 | 11.85 | 0.305 | 0.445 | 0.0727 | 0.0176 | 0.0118 | 0.0936 |

1. **Rank is nearly irrelevant.** k_eff 31 → 451 (14.5×, 478 extra columns) moves the effect −14%.
2. **Cancer type alone produces essentially all of it** (0.0844 vs an anchor of 0.0748 — the anchor is
   in fact slightly *lower*). The 478 tissue-source-site columns add nothing.
3. **What tracks it is how much the design explains of BOTH sides.** `tss_pool10` has more than twice
   the rank of `cancer` but half the induced correlation, because its R² over the two modalities is
   0.166/0.245 instead of 0.297/0.473. `k_eff_shared` — the effective dimension of the *jointly*
   explained subspace — moves only 9.0 → 11.9 across the whole ladder, which is why the effect is flat.

**The practical warning is therefore not "keep your nuisance model small". It is: the induced
correlation is set by how much of both modalities the confound explains, not by how many columns you
spend on it, and it does not shrink with n.**

## 5. Predictors

Ratio measured / predicted; a perfect predictor is 1.0.

| predictor | cells | median ratio | p10–p90 | log₁₀ RMS error |
|---|---:|---:|---|---:|
| P1 `0.6745·κ/√k_eff` | 167 | 1.98 | 1.15 – 3.98 | 0.365 |
| P2 `0.6745·R_s R_a/√k_eff` | 167 | 2.98 | 1.65 – 6.41 | 0.526 |
| **P3 `0.6745·R_x R_y/√k_eff_shared`** (real designs) | 169 | **0.886** | **0.76 – 1.07** | **0.079** |
| P3, all cells including structureless | 241 | 0.839 | 0.39 – 1.05 | 0.204 |

P1 and P2 under-predict by 2–3× and drift. P3 predicts within ~±20% across the whole grid from the
(design, X, Y) spectra alone — no draws, no spike, no recovery curve — which makes an induced-correlation
floor computable for a new cohort *before* the instrument is run.

**Provenance of P3, stated and not hidden:** it was written after P1 and P2 failed on the rank ladder at
n = 2,530, seed 42, and is therefore **post hoc on the rank axis**. What is genuinely out of sample here
is every other n (500 – 6,427) and every other seed (43, 44), and it holds there. The predeclaration
records this; it is not rehabilitated by the fit.

## 6. Cross-fitting suppresses the classical term — which is why the structural term is visible

For structureless designs the classical N−R sampling scale is `0.6745/√(n − R)`. Measured / classical
across 24 structureless cells: **median 0.379** (p10 0.238, p90 0.699).

Cross-fitted residualisation delivers only about **a third** of the in-sample degrees-of-freedom
inflation the Winkler/Muirhead literature describes. This matters for the framing: the structural effect
we report is ~20× what remains of the classical mechanism *because* we residualise correctly. An
in-sample analysis would see the two mixed together.

The one regime where cross-fitting cannot help is a design that nearly spans the sample:
`gaussian_k600` at n = 500 (k_eff = 499) gives 0.055 — the largest structureless value in the grid,
though still far below the classical prediction of 0.617 there.

## 7. Estimator robustness (T2.5) — the predeclared prediction P6 FAILED

Anchor design, induced correlation:

| n | n_splits | α = 0.01 | α = 1.0 | α = 100 |
|---:|---:|---:|---:|---:|
| 2,530 | 2 | 0.0828 | 0.0807 | 0.0380 |
| 2,530 | 5 | 0.0825 | 0.0817 | 0.0531 |
| 2,530 | 10 | 0.0805 | 0.0802 | 0.0573 |
| 2,530 | 20 | 0.0815 | 0.0809 | 0.0592 |
| 6,427 | 2 | 0.0922 | 0.0910 | 0.0632 |
| 6,427 | 5 | 0.0901 | 0.0904 | 0.0725 |
| 6,427 | 10 | 0.0905 | 0.0907 | 0.0760 |
| 6,427 | 20 | 0.0893 | 0.0894 | 0.0772 |

Predeclared: < 25% relative movement across the full grid. **Measured 55.7% at n = 2,530. FAILED.**
The falsifier's own ">2× movement" threshold is also marginally breached (max/min = 2.18 at n = 2,530);
that is reported rather than rounded down.

Decomposed:
* **fold count is irrelevant** — ≤ 2.4% spread over `n_splits` ∈ {2, 5, 10, 20} at fixed α;
* **shrinkage over the usable range is irrelevant** — α 0.01 vs 1.0 differ by ≤ 3%;
* **the entire failure is α = 100**, a 30–53% reduction.

The mechanism is not an artefact and must not be sold as robustness: α = 100 on a one-hot design at
n = 2,530 heavily under-fits the nuisance model, the design then explains less of *both* modalities
(R_s, R_a fall), and Equation (1) requires the induced correlation to fall with them. **Under-adjusting
reduces the induced correlation and leaves the confound in — that is a trade-off, not robustness.**

The quotable statement is: **the induced correlation is invariant to the cross-fitting scheme and to
shrinkage in the range anyone would use (0.07–0.09 for α ∈ [0.01, 1] at any fold count), and it scales
with how much the nuisance model actually removes.**

For the structureless `gaussian_k99` arm the relative spread looks enormous (1.70) on absolute values of
0.0010–0.0079; relative spread is meaningless there and is reported only so its meaninglessness is
explicit.

## 8. T2.6 — both floors as functions of (design rank × n), and the result changes cohort sizing

Full level grid (0.0 … 0.50, 40 draws, 2 seeds), same driver, so the floors cannot come from a second
implementation.

| design | n | k_eff | induced | **detection floor** | transmission floor | attenuation | ambient top-CCA |
|---|---:|---:|---:|---:|---:|---:|---:|
| cancer | 1,000 | 30.9 | 0.0847 | **0.25** | 0.01 | 1.124 | 0.684 |
| cancer | 2,530 | 31.0 | 0.0857 | **0.30** | 0.01 | 1.109 | 0.668 |
| cancer | 6,427 | 31.0 | 0.0844 | **0.30** | 0.01 | 1.099 | 0.663 |
| cancer + tss_pool50 | 1,000 | 30.9 | 0.0847 | 0.25 | 0.01 | 1.124 | 0.684 |
| cancer + tss_pool50 | 2,530 | 31.0 | 0.0857 | 0.30 | 0.01 | 1.109 | 0.668 |
| cancer + tss_pool50 | 6,427 | 50.0 | 0.0842 | 0.30 | 0.01 | 1.101 | 0.663 |
| cancer + tss_pool10 | 1,000 | 45.3 | 0.0799 | 0.25 | 0.01 | 1.112 | 0.683 |
| cancer + tss_pool10 | 2,530 | 102.3 | 0.0774 | 0.25 | 0.01 | 1.091 | 0.664 |
| cancer + tss_pool10 | 6,427 | 215.1 | 0.0811 | 0.25 | 0.01 | 1.086 | 0.660 |
| cancer + tss_pool1 | 1,000 | 316.6 | 0.0806 | 0.25 | 0.01 | 1.070 | 0.677 |
| cancer + tss_pool1 | 2,530 | 446.2 | 0.0790 | 0.30 | 0.01 | 1.086 | 0.661 |
| cancer + tss_pool1 | 6,427 | 579.9 | 0.0831 | 0.30 | 0.01 | 1.080 | 0.654 |
| **gaussian_k99** | 1,000 | 99.0 | 0.0097 | **0.050** | 0.01 | 1.000 | 0.877 |
| **gaussian_k99** | 2,530 | 99.0 | 0.0040 | **0.015** | 0.01 | 1.000 | 0.860 |
| **gaussian_k99** | 6,427 | 99.0 | 0.0014 | **0.010** | 0.01 | 1.001 | 0.863 |

**The detection floor is set by the induced correlation, not by n.** For every real design it sits at
0.25–0.30 and *does not improve* as n goes from 1,000 to 6,427 — a 6.4× increase in sample size buys
nothing. For the structureless Gaussian design of comparable rank it falls 0.050 → 0.015 → 0.010, i.e.
roughly like the sampling scale, exactly as it should when there is nothing structural to hit.

This is the most operationally consequential result in Track 2. **An external cohort cannot buy a lower
floor by recruiting more patients**, if its confound design explains both modalities the way cancer type
does here. What lowers the floor is reducing `R_x·R_y/√k_eff_shared` — i.e. a cohort whose nuisance
structure is less predictive of both modalities — not a bigger n. Track 3's cohort sizing must be
written against the induced-correlation prediction (P3), not against a power calculation in n.

The transmission floor is 0.01 — the finest level on the grid — in every cell, so it is censored from
below and should be reported as "≤ 0.01", not "= 0.01". Attenuation is 1.07–1.12 for real designs and
1.000 for the Gaussian control across all n.

## 9. Recorded substitutions and what is not done

* **No TCGA purity table exists on either machine.** The `["cancer","tss","purity"]` rank point named in
  the plan could not be built. An expression-derived surrogate was rejected because it is computed from
  the very RNA targets that form Y and would inflate R_a by construction — i.e. it would manufacture the
  effect under study. Replaced by `dx_year` (numeric and categorical) and `age` from the TCGA PanCan
  clinical mirror, plus the TSS pooling threshold, which moves rank over an order of magnitude while
  holding the *kind* of confound fixed. Recorded, not silently substituted.
* **Second representation (d2_i)** sweep queued, not complete — the whole grid is one artifact
  (`d2_h_seed42`) and one state (`wsi_biology`), with three seeds per cell.
* Both floors are **censored** at the grid edges: transmission at 0.01 (finest level tested) and, in the
  dilution run, detection at the grid's top. Report as inequalities.
* Single artifact (`d2_h_seed42`), single state (`wsi_biology`) for the main grid; three seeds per cell.

---

**Logged:** 2026-08-03, 01:22–03:30 UTC.

**How obtained:** Lambda A100 box `ubuntu@150.136.45.194`, workspace `~/ws_p1`, `~/venv/bin/python`,
CPU only, 270 + 72 cells. `python -m morpheus.v2.calibra.induced_correlation_sweep --artifact
d2_h_seed42.npz --targets frozen_rna_targets.npz --state wsi_biology --partition all --n-grid
500,1000,2000,2530,4000,6427 --seeds 42,43,44 --n-draws 40 --n-jobs 8 --clinical-covariates
tcga_clinical_covariates.parquet --tag main`, then `--n-splits-grid 2,5,10,20 --alpha-grid
0.01,1.0,100.0 --tag knobs`; graded by `p1_evidence/grade_t2.py` against
`v2/research/rebase/nature/P1_PREDECLARATION.md` (commit 1c4b4b5). Winkler full text from
[PMC7573815](https://pmc.ncbi.nlm.nih.gov/articles/PMC7573815/). Outputs under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/track2/`.

### Technical

Level-0 spike recovery over a 15-design × 6-n × 3-seed grid on 6,427 TCGA patients. The Yule/FWL
identity reproduces the pipeline draw for draw to 8.6 × 10⁻¹⁶ over all 270 cells. At matched rank and
n = 2,530 the real cancer+TSS design induces 0.0748 while its row-permuted twin induces 0.0037 and a
Gaussian design of the same width induces 0.0035 (20.4×, widening to 35.4× at n = 6,427). Induced
correlation is 0.0866 → 0.0718 across n = 500 → 6,427 (ratio 0.960 vs a sampling expectation of 0.627)
and 0.0844 → 0.0727 across k_eff = 31 → 451. Fitted log–log exponents +0.288 in k_eff and −0.180 in n,
against a predeclared +1 and −1. `0.6745·R_x·R_y/√k_eff_shared` predicts the measured value to a median
ratio of 0.886 (p10–p90 0.76–1.07) using only the (D, X, Y) spectra. Structureless arms sit at 0.379 of
the classical `0.6745/√(n−R)` scale. Estimator sweep: ≤ 2.4% over fold count, ≤ 3% over α ∈ {0.01, 1},
30–53% reduction at α = 100 — the predeclared 25% bar fails.

### In plain terms

Removing the same background variables from two unrelated measurements accidentally makes the leftovers
look related. The formula for this has been known since 1907 and our pipeline reproduces it to fifteen
decimal places. The question was whether the size we see — about 0.08 — is just the familiar
statistical cost of using up degrees of freedom, which would make it old news, or something with more
substance.

It has more substance. Take the exact same nuisance model and shuffle its rows so it keeps its size but
loses any connection to the patients, and the effect falls from 0.075 to 0.0037 — twenty times smaller.
Make the nuisance model fifteen times bigger and the effect barely moves. Collect thirteen times more
patients and it barely moves. What does move it is how much of *both* the image side and the molecular
side the nuisance model explains — and for us that is almost entirely cancer type, not the 478 hospital
columns.

We also got things wrong and are saying so. We predicted the effect would grow with the size of the
nuisance model and shrink with the number of patients. It does neither. And we predicted that turning
the knobs of our own correction would barely change the number; turning the shrinkage knob all the way
up halves it — for the same reason everything else here works, namely that the effect tracks how much
the correction actually removes.

### Meaning for the claim

* **The magnitude claim survives the prior-art check and is now defensible**, provided it is worded as
  *quantifying, in correlation units, the effect whose inferential consequence Winkler et al. (2020)
  characterised* — never as discovering the phenomenon — and cited alongside Yule/FWL and
  Muirhead/Anderson.
* **The structural-vs-degrees-of-freedom separation is the content.** Without P5 passing, the classical
  N−R result would cover everything we have. With it, we can say the effect is 20× the residual
  classical term and does not vanish with n.
* **The floor becomes portable.** P3 lets a new cohort compute its expected induced correlation from its
  own design and data spectra before the instrument is run — the certificate field P4 needs and the
  cohort-sizing rule Track 3 needs.
* **The ledger row "confound adjustment does not destroy signal (attenuation 0.94–1.23)" is hardened**:
  its stated falsifier was "attenuation far from 1 under a differently constructed confound design at
  comparable rank", and eleven designs spanning 15× in rank all sit at the same induced baseline with
  every cell passing the null-like gate.
* **P3's D3 is helped:** since adding columns changes the induced baseline by ~1%, any channel change D3
  sees when purity enters the design is attributable to purity rather than to rank.
* **A correction we must publish about ourselves:** the plan's `k/n` prediction and the module's first
  two derivations were all wrong. Design rank is the wrong axis.
