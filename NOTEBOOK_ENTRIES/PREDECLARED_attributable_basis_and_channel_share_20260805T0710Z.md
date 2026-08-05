## 2026-08-05 07:10 UTC — PREDECLARED: is 29/128 a property of the PCA basis or of the biology, and do the 29 certified axes carry the channel?

**Logged:** 2026-08-05 07:10 UTC, **before any rotation is fitted, any certificate is re-run, or any
channel number is read.** Written from
`NOTEBOOK_ENTRIES/post_pbs_constructions_result_20260804T2300Z.md` §1 (the 29/128 result being
interrogated) and from the source it names: `v2/causal_attribution.py`,
`v2/perturbation_basis_common.py`, `v2/calibra/spectral.py`, `v2/calibra/residualise.py`, and
`p1_evidence/baseline_paired_bootstrap.py`. Nothing has been computed at the time of writing beyond
reading the committed `axis_attribution.csv` / `attribution_summary.json` from
`/lambda/nfs/.../p3_post_pbs/results/attribution/` to confirm the column names and the certified list.

### The two questions

**Q1 (basis).** The 29/128 certificate count was measured on the 128 development-fit PCA axes. PCA
axes are variance-maximising and mutually orthogonal *by construction*, and nothing requires a real
causal programme to align with one. Does a basis chosen to be **attributable**, spanning the
**identical** 128-dimensional subspace, certify more than 29?

**Q2 (channel share).** Never measured: of the morphology→molecular channel the block actually
carries, how much lives in the 29 certified axes? Against a size-matched random draw of 29 of the
128, and against the 29 most legible axes (which §1 already shows are mostly *not* the certified
ones — the overlap is 5 of 29 on `d2_h wsi_biology`).

---

### Q1 — design, fixed now

Every arm is an **orthogonal rotation `R` (128×128, `RᵀR = I`) of the frozen PCA loadings**:
`loadings' = loadings @ R`, `scores' = scores @ R`. The span of the 128 gene-space directions is
therefore **exactly** the span the 29/128 was measured on — mean squared canonical cosine against
span(PCA) is 1.0 by construction and will be asserted as such, not assumed. Information content is
held fixed; only *which directions inside it get scored* changes. That is what makes this a fair test
of "does basis choice matter" rather than a different experiment.

The certificate is `causal_attribution.CERTIFICATE` **unchanged**: same four conditions, same
thresholds (`shuffle_rank_spearman <= 0.20`, `cross_line_rank_spearman >= 0.30`,
`attributed_set_coherence_percentile > 0.95`, `r2_cv >` the random-direction null's max), same
`gene_fold_ridge_r2`, same `atom_cosines`, same `attributed_set_coherence`, same
`certifiable_attribution`, same seed 0, same 8,403 K562 atoms / 7,072 genes, same 1,764 shared atoms /
6,207 shared genes for RPE1. The Haar-random direction null and the gene-label shuffle permutation are
drawn from the same `default_rng(seed)` in the same order, so **the nulls are literally identical
across arms**. No statistic is reimplemented; the rotation is applied to the loadings and the scores
immediately after the PCA fit and after the `matches_frozen_pca_block` verification, and everything
downstream is the existing code path.

**Arms.**

| arm | rotation | sees a certificate quantity? |
|---|---|---|
| `none` | identity — reproduces the published 29/128 | — |
| `varimax` | classical orthogonal varimax on the 128 gene loadings | **no** |
| `r2opt` | orthogonal rotation maximising Σ_k R²_cv(k) on the gene-fold out-of-fold residuals | yes (condition 1) |
| `ica` | FastICA on the 128 loadings with symmetric decorrelation, whitening off (orthogonal unmixing) | **no** |
| `xline` | orthogonal rotation maximising aggregate K562↔RPE1 atom-cosine agreement | **yes (condition 3 — the binding one)** |

`varimax` and `ica` are the honest arms: neither ever sees `r2_cv`, the shuffle, the cross-line
Spearman or the coherence percentile. `r2opt` is the arm the task literally names. `xline` targets the
binding constraint and is **circular by construction** — see the distrust section, where its
non-circular form is specified in advance.

**A structural fact stated before the numbers, because it constrains what `r2opt` can possibly do.**
The gene-fold ridge is linear in the target block: with a fixed alpha and fixed folds the out-of-fold
prediction operator `P` is fixed, so the residual of a rotated target is `E@R` where `E = T − P(T)`.
Per-axis `R²_k = 1 − (RᵀAR)_kk / (RᵀBR)_kk` with `A = EᵀE`, `B = T_cᵀT_c`. PCA loadings are unit-norm
and mutually orthogonal, so restricted to the 7,072 aligned genes `B ≈ σ²I`; when `B` is exactly
`σ²I`, `Σ_k (RᵀAR)_kk = tr(A)` is **rotation-invariant**, i.e. *the mean R²_cv cannot be raised by any
rotation at all — it can only be redistributed across axes.* Combined with the fact that condition 1
already passes **128 of 128**, this means **`r2opt` cannot raise the certified count through the
condition it optimises.** It is run anyway, and reported whatever it does, because it is the arm the
question names and because "the obvious objective is provably the wrong one" is itself the answer to
part of the question.

**Predictions (Q1).**

1. `none` reproduces **29** certified axes exactly, and `matches_frozen_pca_block = True` with max
   absolute difference < 1e-3 against `pca_targets.npz`. If it does not, nothing else in this entry
   is readable and I stop.
2. `r2opt` changes mean R²_cv by **less than 0.01 in absolute value** (the invariance argument above),
   and certifies **between 20 and 40** axes. I predict it does **not** beat 29 by a margin I would
   call real.
3. `varimax` certifies **between 15 and 45**, i.e. I predict *no decisive win*, with my point guess
   slightly **below** 29. Reason: cross-line replication is the binding constraint and is
   **uncorrelated with explained variance (+0.000) and near-uncorrelated with R²_cv (+0.06)** in the
   source result, so it is not obviously a quantity a variance-blind sparsity criterion can buy; and
   varimax concentrates loadings on few genes, which *reduces* the number of genes over which the
   cross-line cosine ranking is computed and should if anything make that ranking noisier.
4. `ica` certifies **between 15 and 45**, same reasoning.
5. The **union** over all five arms of "certified directions that are not close to any certified PCA
   axis" is what would show PCA missed something. I predeclare the measurement: for each certified
   rotated axis, `max_j |cos(axis, PCA_j certified)|` and the total squared cosine onto the span of
   the 29 certified PCA axes. I predict the certified rotated axes lie **mostly inside** that
   29-dimensional span — median total squared cosine onto span(certified PCA) **above 0.5** — i.e.
   the same underlying signal re-expressed, not new directions.
6. **What would mean "PCA is the ceiling, basis choice doesn't matter":** every honest arm
   (`varimax`, `ica`) lands within ±6 of 29, *and* certified rotated axes sit mostly inside the span of
   the certified PCA axes (prediction 5). **What would mean "PCA was undercounting":** an honest arm
   certifies **≥ 45** of 128 (a >50% relative increase, well outside the ±6 I would call noise) *and*
   a material share of those certified directions lie outside the certified-PCA span (median total
   squared cosine onto it **below 0.5**), *and* the certificate's own internal structure still looks
   right — condition 1 still passes ~128/128 and cross-line is still the modal failure. A count that
   rises while the failure profile inverts is a bug signature, not a discovery.

**What would make a favourable Q1 result untrustworthy — stated before the run.**

* **Circularity on the binding condition.** `xline` optimises the exact quantity condition 3 scores.
  Its plain certificate count is therefore an **upper bound with no evidential value**, and will be
  labelled as such wherever it appears. Its non-circular form, fixed now: the 1,764 shared atoms are
  split 50/50 by a seeded permutation; the rotation is fitted **only** on fold A, and a
  **held-out-atom cross-line Spearman** is computed on fold B. That held-out column is computed for
  **every** arm including `none`, so the comparison is like-for-like. Only the fold-B number may be
  read as evidence for `xline`. If `xline` beats `none` on the plain count but not on the fold-B
  count, that is a demonstration that the certificate is optimisable — worth reporting, but not a
  finding about biology.
* **A rotation that is not a rotation.** If `RᵀR = I` fails to 1e-10, or if the mean squared canonical
  cosine between span(loadings@R) and span(loadings) is not 1.0 to 1e-10, the arm is void: it changed
  the subspace and is no longer a basis-choice test. Asserted for every arm before its count is read.
* **A count that rises because the nulls moved.** The certified count can be inflated by a
  random-direction null max that drifts down, or by a coherence null that drifts. The
  random-direction null and its max, the shuffle null median, and the coherence null median are
  reported per arm; if the random-direction null max differs across arms by more than 1e-9 something
  is wrong, because the null does not depend on the rotation.
* **Degenerate concentration.** A rotation is free to produce a few very high-`R²` axes and 120 near-zero
  ones. If the certified axes of an arm have a median explained-variance share more than 5× the
  certified PCA axes', I will report that the arm bought its count by concentrating variance, and
  quote the per-arm distribution of `r2_cv` (median, p10, p90) alongside the count.
* **Optimiser noise masquerading as signal.** `r2opt` and `xline` are fitted by iterative ascent on the
  orthogonal group. Both are run from **three** starts (identity plus two seeded random rotations) and
  the spread of certified counts across starts is reported. A win smaller than that spread is not a win.

---

### Q2 — design, fixed now

**Statistic: the project's canonical channel measurement, imported, not restated.**
`calibra.spectral.heldout_top_cca(cross_fitted_residuals(x, design, seed), cross_fitted_residuals(y,
design, seed), n_components=k, seed=seed)` — which is exactly
`run_calibra._channel_measurement`'s `channel_statistic`, and exactly what
`p1_evidence/baseline_paired_bootstrap.py` scores every target block with. Confound design:
`confound_design(cancer, pooled_tissue_source_site)`, the 108-column cancer + pooled-TSS design;
partition `test`, 2,766 patients; artifacts `d2_h_seed42.npz` and `d2_i_seed42.npz`; states
`wsi_biology` and `full_biology`; seed 42; `n_components=16` as the headline, the value the block
comparison uses. `cross_fitted_residuals` is column-wise, so residualising the 128-column block once
and then subsetting columns is *identical* to residualising each subset — asserted in a test rather
than assumed.

**Arms (subsets of the 128 frozen PCA score columns, all size 29):**

* `certified` — the 29 axes with `causal_name_certified == True`.
* `random` — **200** size-29 draws without replacement, seeded.
* `most_legible` — the 29 axes with the largest `legibility__<artifact>__<state>`, chosen **per cell**
  from that cell's own legibility column.
* `all_128` — the full block, as the denominator.

Reported per cell: the channel value for each arm, the certified arm's **percentile within the
200-draw random null**, and the ratio `channel(subset) / channel(all_128)`. Because the statistic
whitens `y` to `n_components` directions, a 29-column block at `k=16` is read through only its own
top-16 principal subspace; so the same table is reported at `k ∈ {4, 8, 16, 29}`, with `k=29` the
value that uses all 29 columns. Any reading that holds only at one `k` is reported as holding only at
that `k`.

**Prediction (Q2).** I expect the **awkward** direction. The certified axes are skewed to high indices
(their median axis index is 80.5 against 63.5 for the block), i.e. to **low-variance** axes, because
condition 3 is uncorrelated with explained variance while conditions 1 and 2 are pushed the other way
by it. The channel statistic is dominated by high-variance directions. So I predeclare:

1. `most_legible` **> ** `certified` in all four cells, by a wide margin.
2. `certified` lands **between the 20th and 70th percentile** of the 200-draw random null — i.e.
   statistically indistinguishable from a random 29 — in at least 3 of the 4 cells. My point guess is
   *slightly below* the random median.
3. `certified / all_128` at `k=16` lands in **[0.55, 0.90]** in each cell.

**If prediction 2 is wrong in the favourable direction** (certified above the 90th percentile of the
random null), the first thing I distrust is the axis-index skew: I would check whether the certified
set is simply enriched for high-variance axes after all, by re-running the null **stratified on
explained-variance decile** so the random draw matches the certified set's variance profile. That
stratified null is predeclared now so it cannot be introduced after seeing a number. If prediction 2
is wrong in the *unfavourable* direction (certified below the 10th percentile) I report it as the
headline: it would mean the certificate actively selects against the channel.

**What Q2 decides.** If `certified` is at or below a random 29, then "29 of 128 axes are named" is
**not** the same thing as "the named axes carry what P4 needs", and the honest P4 exposure is a
coverage number *plus* a channel-share number, not a coverage number alone. If `certified` is well
above random, 29 named axes may already be sufficient for P4 and the coverage number understates it.

---

### Kill switches

* Q1 is void if `none` does not reproduce 29 certified axes.
* Q1 is void per-arm if that arm's rotation fails the orthogonality or the span-identity assertion.
* Q2 is void if the residualise-then-subset vs subset-then-residualise equivalence test fails.
* Q2's random null is reported as uninformative if its interquartile range at `k=16` is below 0.01,
  because then no subset could be distinguished from any other.

### Rules I am holding myself to

No statistic is reimplemented: `gene_fold_ridge_r2`, `atom_cosines`, `attributed_set_coherence`,
`certifiable_attribution`, `subspace_alignment`, `heldout_top_cca`, `confound_design`,
`cross_fitted_residuals`, `pooled_tissue_source_site` are all imported from where they already live.
The only extension to an existing canonical function is an **optional, default-off** return of the
out-of-fold prediction from `gene_fold_ridge_r2`, needed to form `A = EᵀE` for `r2opt` without a
second copy of the ridge; with the flag off its return value is unchanged, which a test asserts.
Any new module that computes an SVD will be allowlisted in
`v2/tests/test_effective_rank_canonical.py` **with its reasoning written inline**, as the previous
post-PBS modules were.

Related: [[post_pbs_constructions_result_20260804T2300Z]],
[[PREDECLARED_past_pbs_constructions_20260804T2240Z]]
