## 2026-08-04 22:40 UTC — PREDECLARED: four post-PBS constructions, their predictions, and what would make a favourable result untrustworthy

**Logged:** 2026-08-04 22:40 UTC, **before any of the four is built or scored.**
**How obtained:** written from `NOTEBOOK_ENTRIES/decision_iterate_past_pbs_and_p4_negatives_20260804T2230Z.md`
(the four candidate constructions, in the order given there) and from
`NOTEBOOK_ENTRIES/t11_t12_must_beat_baselines_20260803T0440Z.md` (the loss being iterated past).
Nothing has been run at the time of writing beyond loading the three Perturb-seq h5ads to read their
shapes.

### What is already measured, and is not being re-litigated

PBS (`pbs_targets_k128_v2.npz`, 128 columns, K562 gwps, 8,403 atoms, 7,072 genes shared with the
prepared PanCan RNA table, 6,427 patients) **loses to ordinary PCA of the same expression matrix** in
3 of 4 artifact×state cells with a paired-bootstrap CI excluding zero and ties the fourth; it beats a
size- and spectrum-matched random dictionary and a gene-label-shuffled copy of itself. That test was
fair, leak-safe and capacity-matched, and is not being repeated on the same construction.

### The harness every construction below is graded through — fixed in advance, identical to the loss

* Target blocks written through `v2/baseline_target_common.write_target_block`, cohort/split/gene
  order/expression transform **inherited** from `pbs_targets_k128_v2.npz` rather than re-derived.
* Any basis is fit on **development rows only** (`split != "test"`) via
  `build_pbs_targets.fit_development_expression_transform`, exactly as the PCA baseline was.
* **128 columns** for every construction: capacity-matched to PBS and to PCA.
* Scored with `p1_evidence/baseline_paired_bootstrap.py`, `--partition test --n-components 16
  --n-boot 400 --seed 42`, artifacts `d2_h_seed42.npz` and `d2_i_seed42.npz`, states `wsi_biology`
  and `full_biology`, 108-column cancer + pooled-TSS confound design, 2,766 held-out patients — the
  same four cells and the same numbers the PBS loss was read on.
* Reused canonical statistics only: `calibra.spectral.heldout_top_cca`,
  `heldout_single_direction_correlation`, `calibra.residualise.confound_design` /
  `cross_fitted_residuals`. No inline restatement of any of them.

### A structural fact that constrains what can possibly change the number

`heldout_top_cca` is **invariant to any invertible linear reparametrisation of the target block**.
Therefore *no rotation, rescaling or reordering of the PBS basis can move the result* — only a change
of the 128-dimensional **subspace of gene space** can. This is written down now because it rules out
a whole family of cosmetic "new constructions" in advance, and because two of the four below (2 and 4)
are only meaningful if they genuinely move the subspace. Each construction therefore also reports the
**principal angles / mean squared canonical cosine between its span and both span(PBS) and span(PCA)**;
a construction whose span is ≥0.99 aligned with span(PCA) is reported as *not a distinct construction*,
whatever its score.

---

### Construction 1 — Attribution, not competition (`v2/causal_attribution.py`)

**What.** Invert the direction of supervision. Take the axes that are *already* legible — the 128
development-fit PCA gene loadings that beat PBS — and ask, per axis, which single-gene CRISPRi
signatures reproduce that gene-space direction.

**Statistics, fixed now.**
1. *Subspace attribution*: `R2_cv(k)` — the gene-fold cross-validated ridge reconstruction of PCA
   loading `q_k` from the 8,403 atom signatures (design = atoms as columns, **genes as rows/samples**,
   5 folds over genes, ridge path fixed in advance at `alpha ∈ {1e-1,1,1e1,1e2,1e3,1e4}`, alpha
   chosen by inner CV). Cross-validation is over **genes**, because with 8,403 atoms and 7,072 genes
   an in-sample fit is guaranteed to be perfect and meaningless.
2. *Atom attribution*: per axis, the ranked cosine between `q_k` and each atom's response signature;
   the named top atoms are the "causal label".
3. *Legibility*: per-axis `heldout_single_direction_correlation(wsi_residual, pca_score_k)` on the
   test partition through the same confound design — a per-axis quantity on the same footing as the
   block-level comparison.
4. *Cross-cell-line stability*: attribution recomputed with the RPE1 atoms (2,326 atoms) and compared
   to the K562 attribution by Spearman of per-atom cosine over the shared atoms.

**Prediction (what I expect).** A minority of PCA axes — I predeclare **fewer than half**, and I
expect the strongest attributions on the proliferation/cell-cycle-flavoured axes — will have
`R2_cv` clearly above the random-direction null; most will not. I expect **no strong positive
correlation across the 128 axes between legibility and attribution strength** (Spearman between
`R2_cv(k)` and legibility(k) between −0.3 and +0.3), because the reason PBS lost is precisely that
the cohort-variance axes and the interventional directions answer different questions. A strong
positive correlation here would be a *better* result than I expect and is exactly the result I would
distrust first (see below).

**What would make a favourable result untrustworthy.**
* `R2_cv` high for **Haar-random gene-space directions** too — i.e. 8,403 atoms spanning enough of
  the 7,072-dimensional gene space that anything is reconstructable. This is the dominant risk and is
  why the random-direction null is mandatory, not optional. **If the random-direction null median
  `R2_cv` exceeds 0.5, the whole subspace-attribution statistic is reported as uninformative and no
  axis is claimed to be attributed.**
* Attribution that **survives permuting the gene labels of the perturbation matrix**. The gene-label
  shuffle must collapse the *atom-level* attribution (Spearman between true and shuffled per-atom
  cosine ranking ≤ 0.05 with a CI covering zero). If it does not, the causal *names* are void even if
  the subspace statistic is fine — the same two-part separation `t15_gene_label_shuffle` applied to PBS.
* Attribution that does not replicate across cell lines (K562 vs RPE1 Spearman near zero) — then the
  "causal name" is a K562 artefact and must be reported as one.
* Legibility and attribution both being driven by the same handful of proliferation genes, so the
  "causal name" of every legible axis is the same name. Reported explicitly as the degenerate case.

**Note this construction does not need to win any contest.** It is graded on whether the
interventional resource supplies a *stable, non-trivial, shuffle-sensitive* causal label for axes PCA
already found legible. A clean negative here ("no legible axis gets a stable causal name") is a
reportable result and will be reported as one.

---

### Construction 2 — Joint basis by CCA between the two covariances (`--construction joint_cca`)

**What.** Stop treating span(PBS) and span(PCA) as an OR. Let `P` be the 128 orthonormal PBS gene
directions and `Q` the 128 orthonormal development-fit PCA gene loadings. SVD `PᵀQ = U S Vᵀ` — this
*is* the CCA between the two covariance row spaces, with `S` the canonical cosines. The joint basis
is `z_k = normalise(P u_k + Q v_k)`, k = 0…127: each column is an equal-weight blend of a matched
canonical pair, so every direction has a declared component inside the interventional span and a
declared component inside the cohort-variance span. Patient codes are `E_dev-scaled @ Z`.

**Prediction.** The canonical cosine spectrum decides this before the CCA number does. I predict the
top few cosines are high (>0.8) and the bulk are low, so `Z` sits nearer PCA than PBS. I predict the
joint block **closes part but not all of the gap to PCA and does not beat it** — my prior is a
difference against PCA in `[-0.03, 0.00]` on `d2_h wsi_biology`, i.e. a tie-to-marginal-loss, and I
would be mildly surprised by a win.

**What would make a favourable result untrustworthy.** (a) mean squared canonical cosine between
span(Z) and span(Q) ≥ 0.99 — then it beat/tied PCA by *being* PCA, and is reported as such rather
than as a construction; (b) a win driven entirely by the handful of high-cosine columns, checkable by
rebuilding with only the low-cosine half and seeing the win vanish; (c) any leak — `Q` is fit on
development rows only and the manifest digest of the fit population is recorded.

---

### Construction 3 — Consensus / denoised interventional basis (`--construction consensus`)

**What.** Separate "the causal framing is wrong" from "this perturbation resource is too noisy". Keep
only atoms whose transcriptional response **reproduces across K562 and RPE1** before the SVD. Atoms
are matched by their Perturb-seq row identifier (verified to share the identifier format across the
two files: 2,326 RPE1 atoms, 8,403 K562 atoms). Retain atom `a` if
`corr(δ_K562,a , δ_RPE1,a)` over the shared genes exceeds the **95th percentile of the mismatched-pair
null** built from 20,000 random `a ≠ b` cross-cell-line pairs — a null that controls for the shared
global structure of the two response matrices, which a bare threshold on `r` would not. The basis is
the top 128 right singular vectors of the **mean** of the two cell lines' responses over the retained
atoms.

**Prediction.** I expect the retained set to be a few hundred atoms — the reproducible core — and I
expect the consensus basis to be **closer to PCA than PBS was but still to lose**, difference against
PCA on `d2_h wsi_biology` in `[-0.05, 0.00]`. Reason for the prior: cross-line-reproducible
perturbations are enriched for core-essential/proliferation genes, which is the one part of the
interventional resource that *does* overlap cohort variance, so this should help — but the retained
atoms are also the ones whose responses are most collinear, so the 128-dimensional basis will be
supported by fewer effectively-independent directions.

**What would make a favourable result untrustworthy.** (a) fewer than ~200 retained atoms, so the
128-component basis is nearly saturated by the retained set and the "consensus" step has become a
rank restriction rather than a denoising step — reported with the retained count and the effective
rank of the retained response matrix, always; (b) the retention filter selecting on a quantity that
correlates with expression magnitude in TCGA (checked by correlating per-atom retention with the
atom's response norm); (c) a win that disappears when the filter threshold is moved one percentile —
so the threshold sweep is run and reported regardless of outcome.

---

### Construction 4 — Domain-adapted projection (`--construction domain_adapted`)

**What.** `build_pbs_targets` divides the CRISPRi delta by the **TCGA development gene SD**
(`reference_delta / development_scale`) and standardises patient expression by the same SD — so the
patient side ends up unit-variance per gene and the perturbation side does **not**. That is the
declared-but-unproven shared-scale assumption named in `fit_development_expression_transform`. The
adapted construction puts both sides on the same footing: the perturbation matrix is z-scored per
gene by **its own across-atom SD**, so both matrices are unit-variance per gene before the SVD. The
128 right singular vectors of that matrix are the basis; the patient encoding is unchanged.

**Prediction.** This changes which genes dominate the SVD (genes that are variable in TCGA but quiet
across perturbations gain weight, and vice versa). I predict it **moves the subspace materially**
(mean squared canonical cosine with span(PBS) below 0.7) and that it **does not beat PCA**;
difference against PCA in `[-0.06, 0.00]`. I give this the lowest prior of the four because
re-weighting genes cannot manufacture alignment with cohort structure that the perturbation responses
do not contain.

**What would make a favourable result untrustworthy.** (a) a per-gene SD floor being hit for many
genes, so the rescaling is dominated by division by near-zero — the count of genes at the floor is
reported; (b) mean squared canonical cosine with span(PCA) ≥ 0.99, i.e. it "won" by becoming PCA;
(c) any improvement that is not stable to using the K562-essential file instead of gwps.

---

### Order, and the stopping rule

1 → 2 → 3 → 4, stopping to report fully **only** if one produces a controlled win (a paired-bootstrap
CI on the difference against PCA strictly above zero in at least the `d2_h wsi_biology` cell, with no
cell showing a CI strictly below zero). Otherwise every construction is run and every result — win,
tie or loss — is written up. **Construction 1 is reported in full regardless**, because it is not
competing for the same slot.

### Compute

All four are CPU. No retraining, no GPU: the largest object is a 8,403 × 7,072 float64 response
matrix (~475 MB) and a 7,072 × 7,072 gene Gram matrix; the D2 artifacts are read, never refit.

Related: [[decision_iterate_past_pbs_and_p4_negatives_20260804T2230Z]],
[[t11_t12_must_beat_baselines_20260803T0440Z]]
