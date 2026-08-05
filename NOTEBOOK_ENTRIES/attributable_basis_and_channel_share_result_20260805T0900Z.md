## 2026-08-05 09:00 UTC — The 29 certified axes carry **less** channel than a random 29, and the honest alternative bases certify **fewer** than PCA — but a basis fitted to cross-line replication certifies **101 of 128**, holds 100 on held-out atoms, survives an atom-fold swap, and is not reproduced by a Haar-random rotation

**Logged:** 2026-08-05 09:00 UTC. **Predeclared in**
`NOTEBOOK_ENTRIES/PREDECLARED_attributable_basis_and_channel_share_20260805T0710Z.md`, committed
(`cef4f39`) **before any rotation was fitted or any channel number read**, and amended (`cbeaa41`)
while the self-tests were being made to pass — before the perturbation data was opened — for a reason
recorded in the amendment itself.

**How obtained.** Lambda box `150.136.45.194`, fresh workspace `~/ws_p3rot` built from
`git archive HEAD` and verified against a hash manifest generated from the canonical checkout. The
six modules this work touches were re-verified against the **committed HEAD blobs after the runs** —
`_require_workspace_matches` → `{"checked": 6, "differ": 0, "missing": 0}` — so every number below was
produced by code that is in the repository. **CPU only**, `OMP_NUM_THREADS` 4–6, the GPU was at
31.8/81.9 GB for another agent's training run throughout and was not touched. Nothing installed into
`~/venv`. `claim_guards.py`, `claim_evidence.json`, other agents' `PREDECLARED_*` files and
`paper/P2_RANK_DRAFT.md` were not edited.

**No statistic is defined by this work.** The certificate is
`causal_attribution.{CERTIFICATE, gene_fold_ridge_r2, atom_cosines, attributed_set_coherence,
certifiable_attribution}` unchanged, at the identical thresholds, with the identical Haar-random
direction null and gene-label shuffle drawn from the same `default_rng(0)` in the same order — the
random-direction null's maximum is **1.352856369195e-03 in every one of the ten arms, spread exactly
0.000e+00**, which is the predeclared check that the nulls did not move. The channel is
`calibra.spectral.heldout_top_cca` on `calibra.residualise.cross_fitted_residuals` under
`confound_design(cancer, pooled_tissue_source_site)` — `run_calibra._channel_measurement`'s
`channel_statistic`, the same object `p1_evidence/baseline_paired_bootstrap.py` scores every target
block with. The one extension to a canonical function is an **optional, default-off** return of the
out-of-fold prediction from `gene_fold_ridge_r2`; a test asserts the returned dictionary is identical
with the flag off.

---

# Bad news first

## 1. The 29 certified PCA axes carry LESS of the channel than a random 29 — and the 29 *most legible* axes carry essentially all of it

This is the number that was never measured, and it is worse than the prediction. Size-matched
comparison, 200 random draws of 29 of the 128 frozen PCA score columns, `n_components = 16`, test
partition, 2,766 patients:

| cell | certified 29 | random-29 median | **certified's percentile in the random null** | most-legible 29 | all 128 |
|---|---:|---:|---:|---:|---:|
| `d2_h wsi_biology` | 0.3341 | 0.4317 | **0.150** | 0.5599 | 0.5520 |
| `d2_h full_biology` | 0.7083 | 0.7820 | **0.200** | 0.8754 | 0.8776 |
| `d2_i wsi_biology` | 0.3325 | 0.4044 | **0.300** | 0.5002 | 0.4905 |
| `d2_i full_biology` | 0.7935 | 0.8130 | **0.415** | 0.8775 | 0.8714 |

The certified set is below the random median in **all four cells**, at `k = 4, 8, 16, 29` alike (at
`k = 4` it falls to the 8.5th–41.5th percentile). The random null is informative — its interquartile
range is 0.03–0.20, far above the predeclared 0.01 kill-switch.

**The most legible 29 carry the whole block.** `channel(most legible 29) / channel(all 128)` is
**0.9888, 1.0142, 0.9974, 1.0198, 1.0070, …** — between 0.94 and 1.07 in *every* cell at *every* `k`,
sixteen values. Twenty-nine of 128 axes chosen by legibility **are** the channel. Twenty-nine chosen
by causal certification are 0.61, 0.81, 0.68 and 0.91 of it. The two selections overlap in 5 axes.

**Why, and how much of it survives the obvious confound.** The certified axes have median axis index
**78.0** against 63.5 for the block — they sit in the low-variance tail, because cross-line
replication is uncorrelated with explained variance (+0.000) while the shuffle condition is pushed the
other way by it. The channel statistic favours high-variance directions. Re-running the null
**stratified on explained-variance decile**, so the random draw matches the certified set's variance
profile, moves the certified percentile up to **0.235, 0.430, 0.385, 0.555** — still at or below the
median, but no longer extreme. So: *most* of the deficit is the variance profile the certificate
happens to select, not the certificate actively selecting against the channel. What is left after
stratification is "indistinguishable from a variance-matched random 29", which is still not "these are
the axes P4 needs".

**Predeclared predictions, graded.** (1) most-legible > certified in all four cells — **held**, by a
wide margin. (2) certified between the 20th and 70th percentile in ≥3 of 4 cells — **held at
`k = 16` (0.150, 0.200, 0.300, 0.415 — three of four at or above 0.20) and FAILED at `k = 4`**
(0.085, 0.165, 0.370, 0.415 — two of four); the point guess "slightly below the random median" was
right in every cell. (3) `certified/all_128` at `k = 16` in [0.55, 0.90] — **held in 3 of 4**
(0.605, 0.807, 0.678, and 0.911 just outside). (4) The favourable-direction distrust check
(percentile > 0.90) never fired, so the stratified null it was declared for is reported as an
explanatory measurement rather than a rescue.

This corroborates, from a completely different statistic, the independent P4 finding logged the same
morning in `composed_readout_and_causal_name_bridge_result_20260805T0745Z.md`: *the 29 certified
causal names do not predict what their axes read.* Two orthogonal measurements now say the certified
set is not the useful set.

## 2. Every basis chosen without looking at the certificate certifies FEWER axes than PCA

Identical certificate, identical thresholds, identical nulls, identical subspace — `mean_squared_cosine`
between each rotated span and span(PCA) is **1.0000** and `max |RᵀR − I|` is below 3.2e-15 for every
arm, so these are basis choices and nothing else.

| arm | what it optimises | **certified / 128** | fold-A | fold-B | fails: recon / shuffle / cross-line / coherence |
|---|---|---:|---:|---:|---|
| `none` (PCA) | variance | **29** | 25 | 31 | 0 / 34 / 78 / 5 |
| `varimax` | sparse gene loadings | **21** | 23 | 24 | 1 / 60 / 97 / 7 |
| `ica` | non-Gaussian independence | **21** | 21 | 21 | 3 / 59 / 97 / 6 |
| `r2opt` | mean gene-fold R²_cv | **27** | 23 | 33 | 0 / 3 / 100 / 6 |

`none` reproduces **29 exactly** — the predeclared kill-switch for the whole of Q1 did not fire, and
`matches_frozen_pca_block = True` at 5.96e-08 against `pca_targets.npz`.

**The objective the question names is provably almost inert, and the data confirms it.** Predeclared
before the run: the ridge is linear in the target block, so `R²_k = 1 − (RᵀAR)_kk/(RᵀBR)_kk`, and when
the target columns are equal-norm (`B = σ²I`) the sum `Σ_k (RᵀAR)_kk = tr(A)` is rotation-invariant —
the mean R² can only be redistributed, never raised. Measured: `B`'s diagonal has a relative spread of
**0.0127**, and the best rotation from three starts moved the mean R²_cv from **0.11508 to 0.11860**,
a change of **+0.0035**, inside the predeclared "< 0.01" band. Since the R² condition already passed
**128 of 128**, `r2opt` could not raise the count through what it optimises — and it did not (27,
against 29). What it *did* do is instructive: it drove shuffle failures from 34 to **3** (the rotated
axes are almost perfectly shuffle-sensitive) while pushing cross-line failures from 78 to **100**.

**`varimax` and `ica` are the honest arms and both land at 21.** Sparsity and independence are
plausible priors for "where a real causal programme points" and both make things *worse*, by the same
mechanism: median cross-line Spearman falls from 0.262 to 0.215 / 0.207 and the shuffle failures
roughly double. Predeclared interval [15, 45] with a point guess slightly below 29 — **held**, though
21 is 8 below 29 and therefore outside the ±6 band I named as "PCA is the ceiling". That miss is on
the side that makes PCA look *better*, so it does not change the reading, and is recorded as a miss.

**Prediction 5 was wrong.** I predicted certified rotated axes would lie mostly *inside* the span of
the 29 certified PCA axes (median total squared cosine > 0.5). Measured: **0.237 (varimax), 0.252
(ica), 0.226 (r2opt), 0.239 (xline)** against a chance value of exactly 29/128 = **0.2266**. The
certified rotated axes are at chance with respect to the certified-PCA span — neither "the same signal
re-expressed" nor "new directions PCA missed", but diffuse mixtures whose overlap with the certified
span is what a random direction would give. Median `max_j |cos(rotated axis, PCA_j)|` over certified
axes is 0.248–0.274, i.e. no rotated certified axis is close to any single PCA axis.

**The degenerate-concentration distrust check did not fire.** Median explained-variance ratio of an
arm's certified axes against its own block median: 0.0050/0.0045 (varimax), 0.0041/0.0044 (ica),
0.0056/0.0053 (r2opt), 0.0049/0.0051 (xline) — ratios of 0.9–1.1, nowhere near the 5× that would have
meant an arm bought its count by concentrating variance.

---

# The one result that goes the other way, and everything done to break it

## 3. A basis fitted to maximise cross-line replication certifies 101 of 128 — and it survives four attempts to explain it away

`xline` maximises a smooth count of axes clearing the certificate's own 0.30 cross-line bar, fitted
**only on atom fold A** (882 of the 1,764 shared K562/RPE1 atoms).

| arm | certified | **fold-B (held-out atoms)** | cross-line Spearman: p10 / median / p90 / sd | fails: recon / shuffle / x-line / coherence |
|---|---:|---:|---|---|
| `none` (PCA) | 29 | 31 | +0.122 / +0.262 / +0.408 / 0.109 | 0 / 34 / 78 / 5 |
| `random` seed 0/1/2 | 31 / 31 / 30 | 33 / 34 / 31 | +0.245 … +0.265 median | 0 / 36–40 / 84 / 0–2 |
| `xline_mean` | 96 | 97 | +0.359 / +0.453 / +0.508 / 0.082 | 0 / 27 / 7 / 0 |
| **`xline`** | **101** | **100** | +0.391 / **+0.431** / +0.453 / 0.070 | 0 / 27 / **4** / 0 |
| `xline` fitted on fold **B** | 102 | 101 (fold A) | +0.424 median / 0.066 sd | 0 / 24 / 4 / 2 |

**The plain count (101) is circular and has no evidential value** — it was declared so in advance. The
numbers that matter:

* **Held-out atoms: 100 of 128.** The rotation never saw fold B's 882 atoms. Predeclared rule: "if
  `xline` beats `none` on fold B by more than the optimiser's across-start spread, cross-line
  replication really was being wasted by the PCA basis; if fold B collapses to at or below `none`, the
  certificate was merely optimisable." Fold B is **100 against 31**, and the across-start spread of the
  objective is 0.9324/0.9362/0.9366. It does not collapse.
* **Atom-fold swap.** Which half is the fitting half is a free choice, so it was swept. Fitted on
  fold B, the certificate reads **101 on fold A** and 102 overall. Symmetric; not fold luck.
* **A Haar-random rotation of the same subspace, fitted to nothing, does NOT do this.** This was the
  control most likely to kill the result, and it was added *after* seeing it, which is stated wherever
  it is quoted. The mechanism it tests: a rotated axis is a dense mixture of PCA directions, and the
  cross-line correlation of a *mixture* of two noisy profiles is higher than its components' by
  Spearman–Brown alone, because independent noise averages down while shared signal adds. Three seeds
  give **31, 31, 30** certified (fold-B 33, 34, 31) and median cross-line Spearman **0.245, 0.265,
  0.249** — indistinguishable from PCA's 29 and 0.262. Density is not the mechanism.
* **The other three conditions did not weaken to pay for it.** Reconstruction still passes 128/128;
  coherence failures fall from 5 to **0** (median attributed-set coherence 0.206 against a random
  atom-set null median of ~0.023); shuffle failures fall from 34 to **27**. The certified axes' names
  still read as machines: `PCA_002` → NACA, TAF12, YEATS4, TAF5, KAT8, MED7, TAF2 (TFIID/mediator);
  `PCA_003` → MFAP1, ECD, TUT1, PRPF4B, SMN2, NAA38, GEMIN5 (spliceosome/SMN); `PCA_000` → POLR1D,
  DHX33, UTP18, DMAP1 (Pol I / nucleolar).

**A predeclared prediction that was wrong, and the algebra behind it.** The amendment predicted
`xline_mean` would be **inert**, because the mean per-axis cross-line correlation looked conserved
(`Σ_k r_kᵀ M r_k = tr(M)`). It moved the mean cross-line Spearman from **0.272 to 0.484** and
certified 96. The premise fails because the per-axis normaliser `sqrt(p_k q_k)` is *not* constant
across rotated directions — the conserved quantity is the unnormalised `tr(M)`, not the mean
correlation. Recorded as wrong, and the module docstring now carries the refutation rather than the
claim.

## 4. …but the extra 72 names do not buy channel *selection*, only channel *coverage*

The obvious next question, and the one that decides whether §3 is useful to P4: do the 101 certified
axes of the `xline` basis carry the channel? Same statistic, random null size-matched to each basis's
**own** certified count (`k = 16`):

| basis | n certified | certified channel | random-n median | percentile | most legible n | all 128 | certified / all |
|---|---:|---:|---:|---:|---:|---:|---:|
| PCA | 29 | 0.3341 | 0.4317 | 0.150 | 0.5599 | 0.5520 | 0.605 |
| `random0` | 31 | 0.4410 | 0.4825 | 0.130 | 0.5239 | 0.5520 | 0.799 |
| `varimax` | 21 | 0.5089 | 0.4571 | **0.870** | 0.5504 | 0.5520 | 0.922 |
| `xline_mean` | 96 | 0.5055 | 0.5375 | 0.065 | 0.5538 | 0.5520 | 0.916 |
| `xline` | 101 | 0.4956 | 0.5422 | 0.045 | 0.5514 | 0.5520 | 0.898 |

(`d2_h wsi_biology`; the other three cells agree in direction. `heldout_top_cca` is invariant to an
invertible reparametrisation of the *whole* block, so `all 128` is 0.5520 in every basis by
construction — only the subsets move.)

Two readings, both true and pulling opposite ways:

* **Selection is still anti-correlated with the channel, in four of five bases.** In `xline` the
  certified 101 sit at the **4.5th percentile** of a random 101 — the certificate is picking, within
  every basis, the axes that carry slightly *less*. The one exception is `varimax` (87th percentile,
  and 74th/84th/98th in the other cells), which is the arm with the *fewest* certified axes.
* **Coverage is genuinely better.** `certified/all_128` goes from **0.605** in the PCA basis to
  **0.898** in the `xline` basis, because naming 101 of 128 axes leaves less unnamed. For a P4
  interface, "101 of 128 axes carry a certified causal name and those axes carry 90% of the channel"
  is a materially better offer than "29 of 128, carrying 61%" — even though within each basis the
  certification is still not selecting *for* the channel.

---

# What this means

**Q1 — is 29/128 a ceiling on the biology or an artefact of the basis?** It is substantially an
artefact of the basis, but not in the direction the question anticipated. No basis chosen *without
reference to cross-line replication* beats PCA: variance (29), Haar-random (30±1), sparsity (21),
independence (21), and aggregate R² (27) all land in a narrow band, and the R² objective the question
names is provably near-inert because the condition it targets already passes 128/128. The binding
constraint is cross-line replication, and **that** is highly basis-dependent: rotating to spread
cross-line agreement across axes rather than concentrate it takes the count from 29 to 101, holds 100
on atoms the rotation never saw, is symmetric under swapping the atom halves, and is not reproduced by
a random rotation of the identical subspace. **29/128 is a property of the PCA coordinate system, not
a coverage limit on the interventional resource.**

The honest caveat, which must travel with that sentence: the `xline` basis was fitted to the
certificate's own condition 3, and the evidence that it is more than optimisation is a held-out-*atom*
split, not a held-out cell line, cohort or gene set. A rotation fitted to make two cell lines agree,
validated by showing two cell lines agree on different perturbations, is a weaker holdout than
this project usually accepts. **The next test, named here and not run: hold out genes, not atoms** —
refit the rotation on half the 6,207 shared genes and score the certificate on the other half. If the
count survives that, the claim is solid; if it collapses, the 101 is a statement about the atom
sample.

**Q2 — do the 29 certified axes carry the channel?** No. They carry less of it than a size-matched
random draw in all four cells, and the 29 most legible axes carry essentially the entire 128-block
channel. Most of the gap is explained by where certification lands in the variance ordering (the
stratified null moves the certified percentile from 0.15–0.42 to 0.24–0.56), but even
variance-matched they are at or below a random draw. **"29 of 128 axes are named" and "the named axes
are what P4 needs" are different claims and only the first is true.** The honest P4 exposure is a
coverage number *and* a channel-share number, not a coverage number alone.

---

# Meaning for the claims, and the exact prose to change

* **P3 — `NOTEBOOK_ENTRIES/post_pbs_constructions_result_20260804T2300Z.md` §1, "Verdict on
  construction 1"** currently reads *"a usable P4 primitive with an honest coverage number (23%)"*.
  That sentence needs two additions and no deletion: (i) 23% is the coverage **in the PCA basis**, and
  the same certificate in a cross-line-fitted rotation of the identical subspace reaches **79%**
  (101/128, 100/128 on held-out atoms); (ii) the 29 certified axes carry **less** channel than a
  random 29. Both belong in the same paragraph as the 29/128, because quoting 29/128 without them
  overstates the ceiling and overstates the usefulness at the same time.
* **Same entry, the table of per-condition failures.** The line "cross-line replication … is an
  independent axis of evidence, which is why it does most of the rejecting" is confirmed and should
  now say *and it is the one condition a change of basis can move* — 78 failures under PCA, 84 under a
  random rotation, 97 under varimax/ICA, 100 under an R²-optimising rotation, and **4** under a
  cross-line-fitted one.
* **Same entry, "Meaning for the claims → P4".** *"The coverage number (29/128) is the honest thing to
  expose in the interface"* should become: the honest thing to expose is the coverage number, **the
  basis it was measured in**, and the channel share of the certified set. A promptable interface that
  offers 29 named axes is offering the 61% of the channel those axes carry, not the 100% the caller
  will assume.
* **`PROJECT_GUIDE.md` §3, P3 "Open"** — both bullets are now closed and are updated in place by this
  entry's commit.
* **P1/P4 claim registry: nothing to change.** No claim record is affected; `claim_evidence.json` was
  not touched. The `no_external_cohort` blocker on per-axis claims is, if anything, reinforced — a
  basis rotation that changes the certified count from 29 to 101 is exactly the kind of analyst
  degree of freedom that a per-axis external replication would have to survive.

---

### In plain terms

We asked two questions about the "29 of 128 axes have a certified causal name" result.

The first was whether 29 is a fact about the biology or a fact about the particular set of directions
we happened to score. The answer is: mostly the latter, but only along one axis. Four different ways
of choosing directions that had nothing to do with the certificate — including a purely random
re-orientation — all gave between 21 and 31, so there is nothing magic about PCA. But choosing
directions specifically so that the two cell lines agree about them takes the count to 101 of 128, and
that holds up when we score it on perturbations the fitting never saw and when we swap which half of
the perturbations we fit on. So the "only a quarter of axes can be named" figure was a limitation of
our coordinate system, not of the perturbation data — with the honest caveat that we validated it by
holding out perturbations rather than a whole new cell line.

The second question was whether the 29 named axes are the ones that actually matter for reading
molecules off an image. They are not. They carry less of that signal than 29 axes picked at random,
and the 29 axes the images read *best* carry essentially all of it. Most of that gap is because
certification tends to land on directions carrying little of the cohort's variation in the first
place. Either way, "we can name a quarter of the axes" was never the same statement as "we can name
the quarter that matters", and we should stop letting the first imply the second.

---

### Files / commits

Modules (all six re-verified against committed HEAD blobs after the runs, `differ: 0, missing: 0`):
`v2/attributable_basis.py` (new), `v2/channel_share.py` (new),
`v2/tests/test_attributable_basis.py` (new, 17 tests), plus `development_pca`,
`cross_line_alignment`, `validated_rotation`, the `--rotation` flag and the two per-atom-fold
cross-line columns in `v2/causal_attribution.py`, `atom_folds` in
`v2/perturbation_basis_common.py`, and one allowlist entry in
`v2/tests/test_effective_rank_canonical.py` (`_polar` computes an SVD and discards the singular
values; it is the projection onto O(n), not a spectrum statistic — reasoning recorded inline).
Commits `cef4f39` (predeclaration), `a3937e7`, `cbeaa41` (amendment + einsum fix), `3716a05`
(bounded ascent), `5a50bab` (random control, fold swap, rotated-block channel share).

Results persisted to NFS at
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p3_attributable_basis/` —
`rotations/{none,random0,random1,random2,varimax,ica,r2opt,xline_mean,xline,xline_foldb}.{npz,json}`,
`results/attribution_*/{axis_attribution.csv,attribution_summary.json}`,
`results/channel_share*/{channel_share.csv,channel_share_summary.json}`, `results/arm_table.csv`,
run scripts, `collate.py`, `arms.py`, `checks.py` and all logs.

**Two bugs found and recorded, both in this work's own new code, both caught by a planted test rather
than by inspection.**
1. `np.einsum("ij,jk,ik->k", R, M, R)` for `diag(RᵀMR)` contracts to `Σ_j (RᵀR)_{jk} M_{jk} =
   diag(M)` — **constant for every rotation**. Both iterative arms would have reported "the objective
   cannot be improved from any start", which is indistinguishable from the conservation argument I had
   predeclared being confirmed. Caught only because a planted test demanded an *increase* and got
   equality to fifteen decimal places from three independent random starts.
2. The ascent's step doubled on every accepted step, so the polar retraction was eventually handed a
   matrix of order 1e300 and `np.linalg.svd` died. Fixed by normalising the tangent and capping the
   step. A separate, genuine LAPACK `gesdd` non-convergence on a well-conditioned 128×128 matrix on
   the loaded box now falls back to the `gesvd` driver.

**Suite** (fresh workspace `~/ws_p3rot2`, built from `git archive HEAD` at `b8ea611` and verified
**855/855 tracked files, differ: 0, missing: 0**, thread-capped, `morpheus/v2/tests morpheus/tests`):
**678 passed, 2 warnings, 27 errors in 76.80s** — **0 failures**, and all 27 errors are
`test_p2_figures.py` setup errors from the known matplotlib import, verified by count
(`grep -c "^ERROR"` = 27; the same grep filtered to non-`test_p2_figures` lines = **0**).
`test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree` passes with the one
new allowlist entry and its inline reasoning. Nothing was installed into `~/venv`.

Related: [[PREDECLARED_attributable_basis_and_channel_share_20260805T0710Z]],
[[post_pbs_constructions_result_20260804T2300Z]],
[[composed_readout_and_causal_name_bridge_result_20260805T0745Z]],
[[PREDECLARED_past_pbs_constructions_20260804T2240Z]]
