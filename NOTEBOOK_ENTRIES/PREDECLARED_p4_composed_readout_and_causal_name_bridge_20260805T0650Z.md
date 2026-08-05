## PREDECLARED 2026-08-05 06:50 UTC — Two tests that separate "no signal" from "no readout", and ask whether a certified causal name generalises

**Status: predeclared. Committed before any of the measuring code exists.** Nothing below is a
result. Every threshold, every null, every fairness control and every condition under which I would
distrust a favourable answer is fixed here, in this commit, and graded verbatim in the result entry.

**Author's standing claim to be tested against:** I expect **(A)**, and I expect the causal-name
bridge to **fail**. Both tests are built so that the outcome I expect is the easy one to report and
the outcome I do not expect has to clear controls I would not be able to add afterwards.

---

## 0. The measured state of play this is written against

From `NOTEBOOK_ENTRIES/p4_inductive_adjustment_measured_20260804T2300Z.md` (the real numbers, not
the 20:00 counterfactual): on the **inductive** arm, `d2_h_seed42::wsi_biology`, test partition,
discovery fraction 0.5, seed 42, exposure fold **n = 1,382** — Policy C answers **0 of 90**;
**with condition 3 relaxed, 28 of 90**, of which **18 of 50** are `hallmark_in_training` and
**2 of 24** are `heldout_pathway`.

Two diagnoses have never been separated:

* **(A) No signal.** The representation carries nothing about untrained pathways.
* **(B) No readout.** The signal is there but distributed over existing axes, and the query rule
  fails because it demands *one dedicated axis per target* instead of composing directions.

Condition 3 is **relaxed throughout both tests below**, exactly as in the entry above, because at
full policy every cell is 0 and no comparison is possible. Nothing here is evidence about
condition 3, and no count below may be quoted as a certified answer count.

---

# TEST 1 — Composed readout vs the single-axis rule, each against its own floor and null

## 1.1 Setting, fixed

Artifact `runs/d2_final/artifacts/d2_h_seed42.npz`, state `wsi_biology`, partition `test`,
`--adjustment inductive`, `discovery_fraction = 0.5`, `seed = 42`, `min_site_count = 10`. This is
the identical prepared state as the measured entry; the harness calls `p4_certify.prepare_state`
unchanged so that it *is* the same rows and the same operator (reference digest
`2060a635fa83756a1c3b7aa8506b7b19fcc4431f5d1a303da39b3cb2bf9d62ce`).

**Targets: all 180 columns of `frozen_rna_targets.npz`**, not 90. The 90 `random_control` columns —
size- and expression-variance-matched random gene sets, one per real signature, named
`RANDOM_CONTROL__<parent>__0` — are the capacity control this test is built around, and they were
excluded from every previous P4 run. **24 of them are the paired controls of the 24
`heldout_pathway` targets**, so the primary endpoint has a per-target matched negative control.

**Primary target set: the 24 `heldout_pathway` targets.** Secondary: `hallmark_in_training` (50),
`immune_tme` (8), `tumour_state` (8). Control: `random_control` (90).

## 1.2 The two readouts. No statistic is defined by this work.

* **Incumbent (single-axis).** `spectral.heldout_single_direction_correlation` on the one supporting
  axis, the axis chosen as `argmax` over 256 of `|·|` — the published rule, byte-for-byte the
  `cmd_competitor` path.
* **Composed.** `spectral.heldout_top_cca(X_all_256_axes, y_target[:, None], n_components = k,
  seed = 42, train_fraction = 0.5)` — canonical held-out CCA: whitening and canonical directions fit
  on the train half, correlation read on the held-out half. This is a **multivariate readout over
  the representation's existing directions**: the target's gene-set score is projected onto the span
  of the axes rather than onto one of them.
* **Capacity ladder, predeclared:** `k ∈ {1, 2, 4, 8, 16, 32, 64}`. **k = 32 is the primary**
  (the canonical default of `heldout_top_cca`). `k = 1` is the same readout family, the same split
  and the same evaluation rows at capacity one — it isolates capacity from readout family.

No new statistic is written. `heldout_top_cca`, `heldout_top_cca_indexed`,
`heldout_cca_projection`, `paired_absolute_correlation`, `heldout_single_direction_correlation`,
`calibration.spike_targets`, `calibration.floors_from_recovery`, `calibration.spike_recovery_curve`,
`residualise.cross_fitted_residuals`, `confound_certificate.within_stratum_permutations` and
`inductive_adjustment.ConfoundAdjustmentOperator` are all used unchanged. The AST scan
(`test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree`) is expected to
pass without a new allowlist entry; if it fires, the offending code is deleted, not allowlisted.

## 1.3 THE FAIRNESS CONTROL — this is the crux of the test

A composed readout has more capacity than a single axis. Grading it against the single axis's null
and the single axis's floor would rig the comparison in its own favour and the result would be
worthless. **Both are matched, in three separate ways, and all three are reported.**

**(i) Matched permutation null.** For each target and each readout, the null is produced by the
**identical readout** applied to targets permuted **within cancer type** by
`confound_certificate.within_stratum_permutations`, 200 draws, seed 42 — the same permuter and the
same draw count as the published run. For the composed arm the whitening and the canonical
directions are **refit inside every permutation**, so the entire capacity of the readout is present
in its own null. A target is "beats its null" iff its observed statistic exceeds the **95th
percentile** of its own 200-draw null, which is the published rule.

**(ii) Matched detection floor — planting the identical signal and changing only the readout.**
`calibration.spike_targets(X_full_256, y_target, level, image_direction = e_a)` plants a spike of
known strength `level` into the **target**, oriented along the supporting axis `a` of the incumbent
rule (`e_a` = the unit vector on that axis in the 256-column space). Every level and draw therefore
puts **the same signal, in the same place, into the same data** for both arms. The spiked target is
residualised by `cross_fitted_residuals` on the scored rows' own design, once, and then read by
**each arm's own readout**. The recovery matrix is turned into a floor by
`calibration.floors_from_recovery` — the same function, the same 80%-of-draws rule, the same
level grid `(0, .01, .02, .05, .10, .20, .40)` and the same 25 draws as `spike_recovery_curve`.
Because the level-0 row is scored by the arm's own readout, a readout with more capacity gets a
higher level-0 upper tail and therefore a **higher floor**, automatically and without my choosing
it. That is the whole mechanism, and it is the reason the floor cannot be borrowed between arms.
The shipped `spike_recovery_curve` floor is reported alongside, unchanged, for continuity with the
published numbers — it is **not** used to grade the composed arm.

**(iii) Matched biological null — the paired random control.** Every count is reported for the
90 `random_control` targets under the identical pipeline. For the 24 primary targets the control is
**paired**: `RANDOM_CONTROL__<that pathway>__0`, same size, variance-matched. If the composed
readout answers random gene sets at a rate near its rate on real held-out pathways, controls (i) and
(ii) have failed and no composed count may be quoted.

**A known asymmetry, declared now rather than discovered later.** `heldout_top_cca` evaluates on a
held-out half (n ≈ 691) while `heldout_single_direction_correlation` is 5-fold and evaluates on all
1,382 rows. The composed arm therefore has the **noisier** evaluation, which widens its null and
makes its bar **harder**, not easier. This biases against the composed arm, i.e. against my
non-expected outcome, and is left in place for that reason. The `k = 1` rung of the ladder is the
comparator that removes it entirely.

**A second known asymmetry, also declared now, and it favours the incumbent.** The incumbent's
supporting axis is chosen by `argmax` over 256 axes *using the rows it is then scored on*, and the
published null holds that axis fixed. That is a selection the published null does not charge for.
I therefore also compute a **selection-aware null** for the incumbent — the full 256-axis grid
recomputed inside every permutation, taking the same `argmax` — and report the incumbent's count
under both. **The primary comparison uses the incumbent's published, generous grading.** If the
composed arm wins, it will not be because I handicapped the incumbent.

## 1.4 Primary endpoint

For each arm: **the number of the 24 `heldout_pathway` targets that clear their own matched
detection floor AND exceed their own matched permutation null p95** (condition 3 relaxed).
Incumbent, published grading: expected to reproduce **2 of 24**; a reproduction failure voids the
run.

## 1.5 What counts as (A) and what counts as (B) — fixed now

| composed (k = 32) answers, of 24 `heldout_pathway` | verdict |
|---|---|
| **≤ 4** | **(A) NO SIGNAL.** Capacity, honestly graded, buys essentially nothing. The representation does not carry a readable channel for untrained pathways, and a better query layer will not create one. Broader supervision is the only route. |
| **5 – 7** | **INTERMEDIATE.** Reported as a partial readout gain and explicitly as neither A nor B. No headline. |
| **≥ 8** | **(B) NO READOUT** — *conditional on the two clauses below.* The signal is present and distributed; the single-axis rule is the wrong shape and is discarding answerable queries. |

**(B) requires all three, not just the count:**

1. composed answers ≥ 8 of 24 `heldout_pathway`, **and**
2. composed answers **≤ 9 of 90** `random_control` (≤ 10%), **and**
3. the composed advantage survives on the **paired** control: strictly more of the 24 real held-out
   pathways answered than of their own 24 matched random controls, with the difference exceeding
   the 95th percentile of a 10,000-draw label-flip null over the 24 pairs.

If (1) holds but (2) or (3) fails, the reported verdict is **"the fairness control failed"** and no
gain is claimed. That is a legitimate outcome of this test and it will be reported as the headline
if it occurs.

## 1.6 What would make me distrust a favourable (B) result

* The `random_control` answer rate rises with `k` at all — even below 10%. A capacity-matched floor
  and null should hold it flat. A rising control rate means my floor is under-penalising capacity.
* The gain is concentrated in targets whose incumbent `|r|` was already within 0.02 of its
  threshold, i.e. the composed readout only tipped near-misses over. Reported per target.
* The gain does not survive `k = 16` and `k = 64` (a knife-edge at one capacity is a fit to the
  level grid, not a channel).
* The composed count on `hallmark_in_training` rises by a *larger proportion* than on
  `heldout_pathway`. That would say composition recovers supervised structure better, which is the
  finding already on record and not a new one.
* Fewer than 3 distinct targets drive the gain.

## 1.7 What would make me distrust an (A) result

* The composed statistic does not exceed the incumbent's on **any** target, including the
  `hallmark_in_training` ones the representation demonstrably does carry. That would mean the
  composed readout is broken, not that there is no signal. **Sanity gate: the composed arm must
  answer at least 25 of the 50 `hallmark_in_training` targets** — the channel that is known to exist
  must be visible through this readout, or the readout is not being tested. If this gate fails the
  run is void and no (A) is reported.
* The composed arm's level-0 floors come back NaN for most targets, i.e. no floor resolves and
  everything fails for want of a floor rather than for want of signal. Reported as a count; if
  > 12 of the 24 primary targets have a NaN composed floor, the (A) reading is downgraded to
  "not measurable at this floor resolution".

---

# TEST 2 — Do P3's certified causal names predict which untrained targets an axis can read?

## 2.1 The hypothesis, stated so it can fail

P3 certified **29 of 128** PCA axes with cross-cell-line-replicated causal perturbation names
(`causal_attribution.CERTIFICATE`, four conjunctive conditions). If a name means anything beyond
being a label, **an axis named for a biological process should read queries about that process and
about processes adjacent to it — including processes never trained on.** If names carry no
generalisation, name-adjacency will not predict readout strength once axis and target main effects
are removed.

## 2.2 Statistic — axis j's image channel, read against target t

`p1_evidence/inputs/pca_targets.npz` (128 PCA axis scores, same 6,427 patients),
`d2_h_seed42::wsi_biology`, partition `test` (n = 2,766), transductively adjusted — the footing
`causal_attribution._axis_legibility` used, so the legibility column in `axis_attribution.csv` is
the matching quantity and no second legibility is computed.

For axis j: `spectral.heldout_cca_projection(X_image, pca_scores[:, j][:, None], train, test)` gives
`px`, the image projected onto the canonical direction **fitted for axis j on the train half only**.
Then

> **r(j, t) = `spectral.paired_absolute_correlation(px, target_t[test_rows])`**

i.e. *how well the image channel that reads axis j also reads target t*, on rows neither the
direction nor the target saw. Canonical machinery only.

Targets: the **40 untrained targets** — `heldout_pathway` (24) + `immune_tme` (8) +
`tumour_state` (8) — this project's standard untrained-40 set.

## 2.3 Adjacency between a name and a target — two definitions, both fixed now

**A1 (primary, objective, no human judgement).** Axis j's name is its top-10 attributed
perturbations (`top_perturbed_genes` in `axis_attribution.csv`). Take those atoms' response vectors
from the K562 genome-wide perturbation matrix via `perturbation_basis_common.load_aligned_response`
and score them against target t's gene set with `causal_attribution.atom_cosines`, the target's
direction being its mean-centred gene-set indicator over the shared gene universe.
**adjacency(j, t) = mean over the axis's 10 atoms of |cosine|.** This asks the biologically correct
question — *does perturbing the genes this axis is named for move the genes of that pathway* —
rather than the sparse question of whether the ten names literally sit inside a 5-to-21-gene KEGG
set. Gene sets: `heldout_pathway` from
`data/processed/genesets/msigdb_discovery_2024.1.Hs.gmt` (all 24 present, verified); `immune_tme`
and `tumour_state` from the frozen signature manifest recorded in the targets' `metadata_json`.

**A1-lit (reported, expected to be underpowered).** Hypergeometric overlap of the same 10 genes with
the target's gene set. With 10 genes against sets of 5–21 over ~7,000, this will almost never be
significant. **Predeclared: if fewer than 10 of the 29 × 40 pairs reach p < 0.05, A1-lit is declared
underpowered and is reported as such, not as a negative result.** It cannot rescue or overturn A1.

**A2 (secondary, judgement-based, written into THIS file before any r(j, t) exists).** A
process-family label for each of the 29 certified axes, assigned from its top-10 gene list alone,
and for each of the 40 untrained targets, assigned from its name alone. adjacency = same family or
an explicitly listed adjacent pair. The map is below and is frozen by this commit.

**Axis families (from `top_perturbed_genes`, this commit):**
`RIBO` ribosome/translation — PCA_004, PCA_024, PCA_075, PCA_078;
`MITO` mitochondrial/OXPHOS — PCA_013, PCA_094, PCA_097, PCA_098, PCA_102;
`REPL` replication/cell-cycle — PCA_031, PCA_032, PCA_085, PCA_108, PCA_111;
`TXN` transcription/chromatin — PCA_038, PCA_068, PCA_104, PCA_109, PCA_112, PCA_122, PCA_126;
`SECR` SRP/secretory/protein targeting — PCA_047, PCA_051, PCA_063, PCA_072;
`SPLC` splicing/RNA processing — PCA_080, PCA_114;
`OTHER` — PCA_007, PCA_027.

**Target families (from names, this commit):**
`MITO` — `..._ELECTRON_TRANSFER_IN_COMPLEX_II`, `..._ELECTRON_TRANSFER_IN_COMPLEX_IV`,
`state_glycolysis`, `state_hypoxia`;
`REPL` — `..._NNK_NNN_TO_CHRNA7_E2F_SIGNALING_PATHWAY`, `..._NNK_TO_DNA_ADDUCTS`,
`..._DCE_TO_DNA_ADDUCTS`, `state_proliferation`, `state_dna_repair`;
`TXN` — `..._EBV_EBNA1_TO_P53_MEDIATED_TRANSCRIPTION`,
`..._EBV_EBNA2_TO_RBP_JK_MEDIATED_TRANSCRIPTION`, `..._TCDD_TO_AHR_SIGNALING_PATHWAY`,
`..._E2_TO_NUCLEAR_INITIATED_ESTROGEN_SIGNALING_PATHWAY`;
`APOP` — `..._EBV_BARF1_TO_INTRINSIC_APOPTOTIC_PATHWAY`,
`..._PARAQUAT_TO_FAS_JNK_SIGNALING_PATHWAY`, `state_apoptosis_senescence`;
`XENO` — `..._BENZO_A_PYRENRE_TO_CYP_MEDIATED_METABOLISM`,
`..._METALS_TO_KEAP1_NRF2_SIGNALIG_PATHWAY`;
`SIGNAL` — the seven RAS/ERK, JNK, NFKB, PI3K and JAK/STAT pathways;
`TRANSPORT` — the two `ANTEROGRADE_AXONAL_TRANSPORT` targets;
`IMMUNE` — all 8 `immune_tme`;
`OTHER_STATE` — `state_angiogenesis`, `state_emt`, `state_mechanotransduction`.

**A2 adjacency = (same family)** OR one of these three explicitly adjacent pairs, listed now:
`RIBO`↔`REPL` (both proliferative), `SECR`↔`MITO` (co-translational targeting to mitochondria),
`REPL`↔`APOP` (damage–death axis). Every other cross-family pair is non-adjacent. `IMMUNE`,
`SIGNAL`, `TRANSPORT`, `XENO` and `OTHER_STATE` have **no** adjacent certified axis family and are
therefore non-adjacent to every certified axis by construction — that is a fact about what P3
certified, and it is part of the result.

## 2.4 Removing the two nuisances

Legible axes read everything; some targets are readable by everything. Both main effects are removed
by **double-centring** the 29 × 40 matrix of `r(j, t)` (subtract row means, then column means, then
add back the grand mean) before any association is taken. The association statistic is the Spearman
correlation between the double-centred `r` and the double-centred adjacency over all pairs.

## 2.5 THE MATCHED-LEGIBILITY CONTROL — predeclared in full

Certified axes must be compared with uncertified axes **of similar legibility**, or the test merely
rediscovers that legible axes read things well.

**Matching rule, fixed now.** Sort the 29 certified axes by
`legibility__d2_h_seed42__wsi_biology` descending. For each in turn, take the nearest unmatched
uncertified axis by absolute legibility difference, **without replacement**, **caliper 0.02**. A
certified axis with no uncertified partner inside the caliper is **dropped from both arms** and the
number dropped is reported. Balance is reported as mean and max |Δ legibility| across matched pairs,
and the matched control arm must satisfy **mean |Δ| ≤ 0.005 and max |Δ| ≤ 0.02** or the control is
declared failed and no comparative claim is made. The same A1/A2 adjacency, the same double-centring
and the same statistic are then computed on the matched uncertified arm.

## 2.6 What counts as the bridge WORKING

All three:

1. **ρ_certified ≥ +0.15**, on A1, over the double-centred pairs; **and**
2. **permutation p < 0.05** against 10,000 draws that shuffle the adjacency vector **across the 40
   targets within each axis** (this preserves both main effects exactly and destroys only the
   name–target correspondence); **and**
3. **ρ_certified − ρ_matched_uncertified > 0** with a 95% bootstrap CI over resampled *axes*
   excluding zero.

## 2.7 What counts as the bridge FAILING

Any of: ρ_certified < +0.15; permutation p ≥ 0.05; or the certified arm not exceeding the matched
uncertified arm with a CI excluding zero. **A failure is a real and publishable result** — it says
the causal names P3 certified are labels on axes, not transferable readout priors, and that P4 may
not use "this axis is named for ribosome biogenesis" as grounds to answer a ribosome query. That is
what I expect and it will be reported as the headline if it occurs.

## 2.8 What would make me distrust a favourable Test 2 result

* **Leave-one-axis-out** changes the verdict for any single axis. Reported for all 29.
* The association is present on the **target-centred-only** matrix but disappears under full
  double-centring — that would mean it is a target main effect (some targets are simply readable and
  happen to be adjacent to something), not a name effect.
* The matched control fails its balance bars (§2.5), or more than 6 of the 29 certified axes are
  dropped for want of a partner inside the caliper.
* A2 and A1 disagree in sign. Both are reported whatever happens; agreement is required before the
  bridge may be called working.
* The certified arm's ρ is driven by `MITO`↔`MITO` pairs alone (the one family with an obvious and
  strong target counterpart). Reported as a per-family breakdown.

---

## 3. Point predictions, to be graded verbatim

| # | prediction |
|---|---|
| 1 | Test 1 incumbent reproduces **2 of 24** `heldout_pathway` |
| 2 | Test 1 composed k = 32 answers **2–4 of 24** `heldout_pathway` → verdict **(A)** |
| 3 | Test 1 composed k = 32 answers **≤ 6 of 90** `random_control` |
| 4 | Test 1 composed k = 32 clears the `hallmark_in_training` sanity gate (≥ 25 of 50) |
| 5 | Test 1 composed matched floors are **higher** than the incumbent matched floors on ≥ 80% of targets |
| 6 | Test 1 incumbent under the **selection-aware** null answers **0–2 of 24** (i.e. the published null is generous but not by much) |
| 7 | Test 2 ρ_certified on A1 lands in **[−0.05, +0.10]**, permutation p ≥ 0.05 → **bridge fails** |
| 8 | Test 2 ρ_matched_uncertified lands within 0.05 of ρ_certified |
| 9 | Test 2 A1-lit is declared **underpowered** (< 10 significant pairs of 1,160) |
| 10 | Test 2 fewer than 3 of the 29 certified axes are dropped by the caliper |

---

## 4. Rules of engagement

* CPU only. The GPU is saturated by another agent's training run and is not touched. Nothing is
  installed into `~/venv`.
* `claim_guards.py`, `claim_evidence.json`, other agents' `PREDECLARED_*` files and
  `paper/P2_RANK_DRAFT.md` are not edited by this work.
* No number here alters a published claim. The result entry names the exact P4/P3 prose locations
  that should change and leaves the edits to their owners.
* The full suite is run before the final commit and its pass/fail counts quoted verbatim.
  `test_p2_figures.py` matplotlib errors are the known condition of `~/venv`.

Related: [[p4_certification_end_to_end_20260804T2000Z]],
[[p4_inductive_adjustment_measured_20260804T2300Z]],
[[post_pbs_constructions_result_20260804T2300Z]]
