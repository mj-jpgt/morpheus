## 2026-08-04 23:00 UTC — Iterating past the PBS loss: attribution works and is nameable for 29 of 128 axes; a joint CCA basis closes 54-86% of the gap to PCA but still never beats it; a cross-line consensus basis does not help at all

**Logged:** 2026-08-04 23:00 UTC. **How obtained:** Lambda box `150.136.45.194`, fresh workspace
`~/ws_pbs2` built from a `git archive HEAD` tarball and **verified against a hash manifest generated
from the canonical checkout** (`_require_workspace_matches` → `{"checked": 681, "differ": 0,
"missing": 0}`), not from a git-diff-based sync. All CPU; no GPU, no retraining. Predeclared in
`NOTEBOOK_ENTRIES/PREDECLARED_past_pbs_constructions_20260804T2240Z.md` **before** any construction
was built. Outputs under `/home/ubuntu/pbs2_out/`, persisted to
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p3_post_pbs/`.

Everything below is scored through the harness that produced the PBS loss, unchanged:
`p1_evidence/baseline_paired_bootstrap.py --partition test --n-components 16 --n-boot 400 --seed 42`,
artifacts `d2_h_seed42.npz` / `d2_i_seed42.npz`, states `wsi_biology` / `full_biology`, 108-column
cancer + pooled-TSS design, 2,766 held-out patients, `calibra.spectral.heldout_top_cca` and
`calibra.residualise.{confound_design,cross_fitted_residuals}` reused rather than restated.

---

### 0. The control that makes the rest readable

`--construction pbs_rebuild` rebuilds the frozen PBS block through the new code path. Against
`pbs_targets_k128_v2.npz`: **max absolute score difference 1.17e-05 (1.5e-07 relative), minimum
per-column |correlation| 0.9999999999996, subspace overlap 1.0000**. The new path is the old path.
A negative from the other three constructions is therefore a fact about the constructions.

The PCA block was re-derived identically too: `matches_frozen_pca_block = True`, max absolute
difference **1.86e-09** against `p1_evidence/inputs/pca_targets.npz`. The axes attributed in §1 are
literally the block that beat PBS.

---

### 1. Construction 1 — attribution, not competition

Take the 128 development-fit PCA gene loadings — the axes that beat PBS — and ask which single-gene
CRISPRi signatures reproduce each one. 8,403 K562 atoms, 7,072 genes, gene-fold cross-validated dual
ridge (folds over **genes**, because with more atoms than genes the in-sample fit is exactly 1).

**The subspace-level result is clean and decisive.**

| arm | selected alpha | median R²_cv | p10 | p90 | max |
|---|---:|---:|---:|---:|---:|
| PCA axes | 1e3 | **+0.1011** | +0.0454 | +0.1974 | +0.4050 |
| Haar-random gene-space directions | 1e4 | −0.0016 | −0.0030 | −0.0002 | +0.0014 |
| gene-label-shuffled | 1e4 | −0.0017 | −0.0030 | −0.0001 | +0.0013 |

**All 128 of 128 axes exceed the random-direction null's *maximum*** (worst axis +0.0244 against a
null max of +0.0014). The predeclared kill-switch — "if the random-direction null median R²_cv
exceeds 0.5 the statistic is uninformative" — did not fire: the null sits at zero. So the atom set is
**not** a generic spanning set of gene space; it specifically reconstructs the cohort-variance axes,
and permuting which gene each loading refers to destroys that completely. **The interventional
resource does supply something about a PCA axis that PCA does not contain.**

**The predeclared distrust check fired, and it changes the reading.** I predicted no strong
association between an axis's legibility and its attribution strength (|ρ| < 0.3) and said that a
strong positive association was the result I would distrust first. The raw association *is* strong
and positive — Spearman **+0.446** (CI95 [+0.292, +0.575]) on `d2_h wsi_biology`, +0.463, +0.376,
+0.401 on the other three cells, all CIs excluding zero. It does not survive:

* R²_cv vs the axis's own explained-variance ratio: Spearman **+0.925**.
* legibility vs explained-variance ratio: +0.48 to +0.57.
* legibility vs R²_cv **partialling out explained variance: −0.180, −0.213, −0.201, −0.236** —
  the association reverses sign.

**So attribution strength measures how much variance an axis carries, not how morphology-legible it
is.** The naive positive correlation is a signal-to-noise artefact and must not be reported as
"legible axes are the causally attributable ones."

**The atom-level result — the causal *name* — is much weaker, and the binding constraint is
cross-cell-line replication.** Per-axis conjunctive certificate
(`causal_attribution.CERTIFICATE`, four conditions each killing a different failure mode):

| condition | axes failing (of 128) |
|---|---:|
| R²_cv above the random-direction null's maximum | **0** |
| atom ranking destroyed by a gene-label shuffle (Spearman ≤ 0.20) | 34 |
| atom ranking replicated from RPE1 (Spearman ≥ 0.30) | **78** |
| named atom set more coherent than size-matched random atom sets (>95th pct) | 5 |
| **all four — certified** | **29 pass** |

Supporting distributions over the 128 axes: top-atom cosine median +0.097 against a gene-label-shuffle
null of +0.050 and a random-direction null of +0.047 (only 24 axes clear the shuffle null's *max*);
shuffle rank-Spearman median +0.058 but p90 +0.336 and max +0.661; cross-line rank-Spearman median
+0.262, p90 +0.408. Attributed-set coherence is the one control that passes overwhelmingly: median
mean-|pairwise-r| **0.229** among an axis's top-10 atoms against **0.023** for random size-10 atom
sets, with **96.1%** of axes above the 95th percentile of that null, and 0.020 for random directions.
The named sets are real machines, not ten unrelated genes.

**Which legible axes get a causal name? Essentially at random.** Median legibility of certified axes
0.0897 vs uncertified 0.0874 (`d2_h wsi_biology`; the other three cells agree to within 0.006), and
**only 2 of the 15 most legible axes are certified.** The mechanism is visible in the components:
shuffle-survival correlates **+0.498** with explained variance, so the high-variance — and therefore
most legible — axes are exactly the ones whose atom ranking survives permuting the gene labels, and
they fail on that condition. Cross-line replication, by contrast, is **uncorrelated with explained
variance (+0.000)** and with R²_cv (+0.06): it is an independent axis of evidence, which is why it
does most of the rejecting.

The 29 certified names are chemically coherent. A sample, with held-out legibility on
`d2_h wsi_biology`, R²_cv, cross-line Spearman and top attributed perturbations:

| axis | legibility | R²_cv | x-line | top perturbed genes |
|---|---:|---:|---:|---|
| PCA_004 | +0.167 | 0.405 | +0.31 | DHX36, RSL1D1, LAS1L, TEX10, ZNHIT1, RPL24, NOL8, RPL14, GNL2 — large-subunit ribosome biogenesis / nucleolar |
| PCA_013 | +0.206 | 0.185 | +0.36 | SSBP1, PNPT1, LRPPRC, TFAM, LONP1, TAMM41 — mitochondrial nucleoid / mtRNA turnover |
| PCA_047 | +0.134 | 0.160 | +0.31 | RAD21, SMC3, SRP68, OXA1L, HSPA9 — cohesin + import machinery |
| PCA_072 | +0.146 | 0.091 | +0.40 | SRP72, SRP68, DDX21, DKC1, MRPL4, RPP14, RCL1 — SRP / nucleolar RNP |
| PCA_007 | +0.291 | 0.146 | +0.45 | PTPN1, DYNLL2, NELFB, SUPT6H, TSC1, LDB1, PSMD13 |

Uncertified axes include the most legible ones: PCA_014 (legibility +0.439), PCA_002 (+0.427, atoms
SRP72/HSPA9/SEC61A1/TOMM22 — a coherent translocon set that nonetheless fails cross-line at +0.11),
PCA_003 (+0.354, ribosomal, fails the shuffle at +0.498).

**Verdict on construction 1.** The interventional resource supplies a real, null-controlled,
shuffle-sensitive causal label for **29 of 128** axes that PCA cannot supply for any of them — and
the labels land on coherent protein complexes. It does **not** supply that label preferentially for
the axes morphology can see, and the raw appearance that it does is an explained-variance artefact.
That is a usable P4 primitive with an honest coverage number (23%), not a legibility win.

---

### 2–4. Constructions 2, 3 and 4 — the same harness, the same four cells

Reference block = the construction; baseline = the block it is being compared with; the difference is
**construction − baseline** on the *same* bootstrap resample.

**Against ordinary PCA of the same expression matrix (the opponent PBS lost to). No construction
wins in any cell.**

| construction | d2_h wsi_biology | d2_h full_biology | d2_i wsi_biology | d2_i full_biology |
|---|---|---|---|---|
| PBS (rebuilt here) | **−0.0488** [−0.0734, −0.0183] LOSS | **−0.0359** [−0.0483, −0.0236] LOSS | −0.0300 [−0.0429, +0.0053] tie | **−0.0080** [−0.0233, −0.0001] LOSS |
| **joint_cca** | −0.0121 [−0.0195, +0.0021] tie | **−0.0051** [−0.0081, −0.0011] LOSS | −0.0138 [−0.0180, +0.0071] tie | −0.0018 [−0.0078, +0.0049] tie |
| consensus | **−0.0559** [−0.0842, −0.0137] LOSS | **−0.0394** [−0.0723, −0.0318] LOSS | −0.0005 [−0.0286, +0.0327] tie | +0.0032 [−0.0214, +0.0078] tie |
| domain_adapted | **−0.0372** [−0.0554, −0.0044] LOSS | **−0.0245** [−0.0291, −0.0115] LOSS | −0.0272 [−0.0400, +0.0093] tie | −0.0095 [−0.0181, +0.0020] tie |

**Against PBS itself (does the construction improve on the thing it replaces?).**

| construction | d2_h wsi_biology | d2_h full_biology | d2_i wsi_biology | d2_i full_biology |
|---|---|---|---|---|
| **joint_cca** | **+0.0367** [+0.0145, +0.0621] WIN | **+0.0308** [+0.0206, +0.0429] WIN | +0.0162 [−0.0077, +0.0344] tie | +0.0062 [−0.0014, +0.0214] tie |
| consensus | −0.0071 [−0.0318, +0.0215] tie | −0.0035 [−0.0371, +0.0019] tie | +0.0295 [−0.0079, +0.0493] tie | +0.0112 [−0.0040, +0.0168] tie |
| domain_adapted | **+0.0117** [+0.0002, +0.0335] WIN | **+0.0115** [+0.0085, +0.0222] WIN | +0.0028 [−0.0125, +0.0184] tie | −0.0016 [−0.0033, +0.0089] tie |

**The `pbs_rebuild` row reproduces the 2026-08-03 table digit for digit** — −0.0488 [−0.0734, −0.0183],
−0.0359 [−0.0483, −0.0236], −0.0300 [−0.0429, +0.0053], −0.0080 [−0.0233, −0.0001] — from a
different workspace, and its difference against the frozen PBS block is 0.0000 in all four cells.
The comparison really is apples-to-apples with what was already measured.

**Distinctness, checked before any score was read (predeclared, because held-out top-CCA is invariant
to reparametrisation and a rotation of PBS could not move the number):**

Mean squared canonical cosine (1.0 = same span, 0.0 = orthogonal), reported in **both** the gene
space the basis is built in and the **patient-code space the readout actually sees** on the 2,766
test rows:

| construction | gene space vs PBS | gene space vs PCA | **code space vs PBS** | **code space vs PCA** |
|---|---:|---:|---:|---:|
| PBS | 1.000 | 0.116 | 1.000 | **0.558** |
| joint_cca | 0.629 | 0.629 | 0.665 | **0.953** |
| consensus | 0.240 | — | 0.525 | **0.534** |
| domain_adapted | 0.753 | — | 0.833 | **0.552** |

None reached the predeclared 0.99 "it won by being PCA" bar, so all three are distinct constructions
by the rule set in advance — **but `joint_cca` at 0.953 in code space is the caveat this table exists
to surface, and it should be quoted with the result.** In the space the statistic reads, the joint
basis is 95% PCA; its 75–86% recovery of PBS's deficit is largely the price of getting there. Two
further things the table says: PBS's gene-space overlap with PCA is only 0.116 yet its code-space
overlap is 0.558, because projecting *any* gene directions onto real patients pulls them toward the
dominant cohort directions; and `domain_adapted` moved the code space away from PBS (0.833) without
moving it toward PCA (0.552), which is consistent with its small, real, but partial gain.

`domain_adapted`'s gene-space overlap with PBS came in at 0.753, above the <0.70 I predicted;
recorded as a miss on my own prediction, not adjusted after the fact.

**Construction 2 — joint CCA basis.** Canonical cosines between the two spans: top 0.900, then 0.856,
0.804, 0.794, 0.751, median 0.198, min 7.4e-05 — i.e. a handful of shared directions and a long tail
of directions each basis has to itself, which is why the equal-weight blend is a genuinely new
subspace. **This is the one basis-side construction that moved anything.** It closes 75%, 86%, 54%
and 78% of PBS's deficit against PCA in the four cells; it beats PBS outright in 2 of 4 with CIs
excluding zero and ties the other 2; against PCA it goes from PBS's 3 losses + 1 tie to **1 marginal
loss (−0.0051) + 3 ties**. It never beats PCA. My prediction was "a difference against PCA in
[−0.03, 0.00] on `d2_h wsi_biology`, tie-to-marginal-loss"; measured −0.0121, tie. Prediction held.

**Construction 3 — cross-line consensus basis.** 1,764 atoms are shared between K562 gwps and RPE1
over 6,207 shared genes; matched cross-line correlation median **0.152** against a mismatched-pair
null median of **0.042**, threshold at the null's 95th percentile **0.255**, **514 atoms retained**.
The predeclared rank-starvation check passes — retained-response effective rank **404.4**, far above
the 128 components taken. The predeclared selection check is a partial concern and is reported as
such: retention correlates **+0.272** with the atom's own response norm, so the filter is mildly
enriched for loud atoms. **The construction does not help.** Against PCA it loses in both `d2_h`
cells — and by *more* than PBS does (−0.0559 vs PBS's −0.0488; −0.0394 vs −0.0359) — and ties both
`d2_i` cells. Against PBS it is a tie in all four cells: filtering to the cross-line-reproducible
core neither helped nor hurt. My prediction was "closer to PCA than PBS was but still losing,
difference in [−0.05, 0.00]"; measured −0.0559 on `d2_h wsi_biology`, i.e. **worse than PBS and
outside the interval I predicted.** Recorded as a wrong prediction. The reading I take from it: the
gap is not explained by the perturbation resource being cross-cell-line noisy — denoising it that way
buys nothing.

**Construction 4 — domain-adapted projection.** The perturbation matrix z-scored per gene by its own
across-atom SD instead of by the TCGA development SD, so both sides carry unit per-gene variance
before the SVD. **0 genes hit the scale floor**, so the predeclared division-by-noise failure did not
occur. It beats PBS in both `d2_h` cells with CIs excluding zero (+0.0117, +0.0115) and ties `d2_i`,
so the shared-scale assumption in `fit_development_expression_transform` **is** costing something
measurable — but only about a quarter of what the joint basis recovers, and it still loses to PCA in
both `d2_h` cells (−0.0372, −0.0245). My prediction was "difference against PCA in [−0.06, 0.00],
lowest prior of the four"; measured −0.0372, inside the interval.

---

### In plain terms

We already knew our interventional gene dictionary reads tumour images *worse* than plain PCA of the
same expression data. We tried four different ways forward and reported all four.

The one that worked is not a competition at all. Instead of asking the perturbation data to be a
better set of axes, we asked it to *name* the axes PCA already found: for each PCA direction, which
single-gene CRISPRi experiments produce that same transcriptional pattern? Every one of the 128 axes
turns out to be reconstructable from perturbation signatures far beyond what a random direction or a
gene-scrambled version achieves, and the genes named for an axis are members of the same cellular
machine — ribosome assembly, the mitochondrial genome, cohesin — rather than an arbitrary list. But
when we insist the name also be reproducible in a second cell line and be destroyed by scrambling
gene labels, only 29 of 128 survive, and they are not the axes the images read best. So: we can put
a causal name on about a quarter of these axes, and we should say "about a quarter", not "we can name
them".

The three attempts to build a *better basis* did not beat PCA. Blending the interventional and
variance bases got most of the way there and clearly beat the original dictionary, but did not
overtake plain PCA in any of the four comparisons — and by the time it had got that far it had become
95% the same thing as PCA anyway, measured in the space the test looks at. Filtering the perturbations
down to the ones that reproduce across two cell lines did not help at all, which at least tells us the
problem is not that the perturbation data is noisy. Putting both datasets on a common per-gene scale
helped a little — a real gain over the original dictionary — but nowhere near enough to overtake PCA.

---

### Meaning for the claims

* **P3.** The withdrawal stands: no construction tested here beats ordinary PCA on legibility. What
  is *new* and defensible is the attribution result — an interventional causal label for a legible
  axis, with three nulls and a stated 29/128 coverage. That is a narrower claim than the original
  P3 and a different one; it is closer to "the perturbation atlas annotates a representation" than to
  "interventional coordinates are more legible".
* **The joint basis is the only basis-side direction worth continuing, and it must be quoted with its
  code-space caveat.** It beat PBS in 2 of 4 cells with CIs excluding zero and tied the other 2, and
  it reduced the deficit against PCA from −0.0488 to −0.0121 on `d2_h wsi_biology` — but it is 0.953
  aligned with PCA in patient-code space, so "a causally-grounded basis that is nearly as legible as
  PCA" is not yet distinguishable from "PCA with a causal component bolted on". The obvious next move
  — sweeping the blend weight between the two spans — is deliberately **not** run here, because
  selecting the weight on the test-partition outcome would be circular; if it is run it must be
  predeclared with a held-out selection fold, and it must report code-space overlap at every weight.
* **P4** gains a concrete certificate field. `causal_attribution.CERTIFICATE` is four conditions with
  numbers, a pass count and a per-condition failure count, which is what an axis needs before a
  promptable interface may answer "axis 4 matches perturbing the large ribosomal subunit". The
  coverage number (29/128) is the honest thing to expose in the interface, and the 78 axes that fail
  cross-line replication are the named failures.
* **Three predeclared predictions were wrong and are recorded as wrong.** (i) The
  legibility–attribution association was strongly positive rather than flat — though it did not
  survive the confound check I predeclared for exactly that reason. (ii) `consensus` came in *worse*
  than PBS against PCA (−0.0559), outside the [−0.05, 0.00] I predicted. (iii) `domain_adapted`
  moved the subspace less than I said it would (mean squared cosine 0.753 against span(PBS), vs the
  <0.70 I predicted). Two predictions held: `joint_cca`'s difference against PCA landed inside
  [−0.03, 0.00], and `domain_adapted`'s inside [−0.06, 0.00].

### Files / commits

Modules: `v2/perturbation_basis_common.py`, `v2/build_causal_basis_targets.py`,
`v2/causal_attribution.py`, `v2/tests/test_post_pbs_constructions.py`, plus
`development_expression_moments` in `v2/baseline_target_common.py` and two allowlist entries in
`v2/tests/test_effective_rank_canonical.py`. All six were re-verified against the committed HEAD
blobs after the runs — `_require_workspace_matches` → `{"checked": 6, "differ": 0, "missing": 0}` —
so the numbers above were produced by the code that is in the repository, not by a drifted copy.
Results, persisted to NFS at
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p3_post_pbs/` —
`results/construction_table.csv`, `results/bootstrap_{joint_cca,consensus,domain_adapted,pbs_rebuild}.{json,csv}`,
`results/attribution/{axis_attribution.csv,attribution_summary.json}`, target blocks in `inputs/`,
run scripts and logs alongside.

**Suite** (box `~/ws_pbs2`, thread-capped, `morpheus/v2/tests morpheus/tests`):
**568 passed, 27 errors in 61.78s** — all 27 errors are `v2/tests/test_p2_figures.py` setup errors
from the known matplotlib import, verified by count (`grep -c` = 27, no non-`test_p2_figures` errors).
Nothing was installed into `~/venv`. The first run of the suite was **1 failed, 567 passed, 27
errors**: `test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree` correctly
flagged the two new modules for containing `linalg.svd`. Both were inspected — principal-angle
cosines and the CCA that defines the construction, neither a spectrum statistic, and the one
effective rank the module needs is imported from `calibra.spectral` — and allowlisted with that
reasoning recorded inline. **The AST scan firing on brand-new code is itself worth recording** — it
does not distinguish a substitution from a legitimate SVD, which is the correct design: it forces the
justification to be written down rather than assumed.

Related: [[PREDECLARED_past_pbs_constructions_20260804T2240Z]],
[[decision_iterate_past_pbs_and_p4_negatives_20260804T2230Z]],
[[t11_t12_must_beat_baselines_20260803T0440Z]]
