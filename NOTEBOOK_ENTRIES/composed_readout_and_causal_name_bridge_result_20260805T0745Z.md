## 2026-08-05 07:45 UTC — It is (A), and the control that decides it was never in the query set: the incumbent rule answers 2 of 24 held-out pathways and **6 of their own size-and-variance-matched random gene sets**. A composed readout raises both together. The 29 certified causal names do not predict what their axes read

**Logged:** 2026-08-05 07:45 UTC. **Predeclared in**
`NOTEBOOK_ENTRIES/PREDECLARED_p4_composed_readout_and_causal_name_bridge_20260805T0650Z.md`,
committed (`b83ac7c`) **before either measuring driver existed** — the drivers landed at `78a4278`,
`9a1f9a7` and `c20ab1a`.

**How obtained.** Lambda box `150.136.45.194`, fresh workspaces `~/ws_p4c` and `~/ws_p4c2` built from
`git -c core.autocrlf=false archive HEAD` and verified **827/827 tracked files by git blob SHA-1**
(`differ: 0, missing: 0`) at the commit each run was launched from. **CPU only**,
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`; the GPU was at
**99%** for another agent's training run throughout and was not touched. Nothing was installed into
`~/venv`. `claim_guards.py`, `claim_evidence.json`, other agents' `PREDECLARED_*` files and
`paper/P2_RANK_DRAFT.md` were not edited.

**No statistic is defined by this work.** Both drivers import
`spectral.{heldout_top_cca, heldout_top_cca_indexed, heldout_cca_projection,
paired_absolute_correlation, heldout_single_direction_correlation}`,
`calibration.{spike_targets, floors_from_recovery, spike_recovery_curve}`,
`residualise.{confound_design, cross_fitted_residuals, pooled_tissue_source_site}`,
`confound_certificate.within_stratum_permutations`, `causal_attribution.atom_cosines`,
`perturbation_basis_common.load_aligned_response`, and `p4_certify.{prepare_state, channel_grid,
load_state, load_targets}` unchanged. **`test_effective_rank_canonical.py::
test_no_second_definition_exists_in_the_tree` passes with no new allowlist entry** — that is the
check that would have caught a parallel implementation, and it did not fire.

---

# Bad news first

**1. The answer is (A) — no signal — and the evidence for it is stronger and more embarrassing than
the A/B question asked for.** The query set every P4 number has ever been computed on
(`frozen_rna_targets.npz`, 90 non-control targets) **excludes the file's own 90 `random_control`
columns** — size- and expression-variance-matched random gene sets, one per real signature. Put them
through the identical, published, condition-3-relaxed grading and:

| | targets | incumbent single-axis rule, **published grading** | composed readout, k = 32 |
|---|---:|---:|---:|
| `hallmark_in_training` | 50 | **18 (36.0%)** | 28 (56.0%) |
| **`heldout_pathway`** | **24** | **2 (8.3%)** | **7 (29.2%)** |
| `immune_tme` | 8 | 5 (62.5%) | 5 (62.5%) |
| `tumour_state` | 8 | 3 (37.5%) | 4 (50.0%) |
| **`random_control`** | **90** | **21 (23.3%)** | **33 (36.7%)** |

**The 24 genuinely untrained pathways are answered at 8.3%. Random gene sets drawn from the same
expression matrix are answered at 23.3%, by the same rule, on the same patients.** The one target
class P4 most needs to answer is the class it answers *least* often — less often than noise.

**2. The pairing makes it exact, and it is not a rate artefact.** Each held-out pathway has its own
matched control, `RANDOM_CONTROL__<that pathway>__0`. Over those 24 pairs, at every rung of the
predeclared capacity ladder:

| grading rule | real pathways answered | **their own matched random sets answered** | difference | null p95 | p |
|---|---:|---:|---:|---:|---:|
| **incumbent, published (the 2 of 24 on record)** | **2** | **6** | **−4** | 4.0 | 0.984 |
| incumbent, matched-readout floor | 5 | 6 | −1 | 5.0 | 0.750 |
| composed k = 1 | 3 | 5 | −2 | 4.0 | 0.894 |
| composed k = 2 | 4 | 8 | −4 | 6.0 | 0.942 |
| composed k = 4 | 7 | 9 | −2 | 6.0 | 0.831 |
| composed k = 8 | 9 | 8 | +1 | 5.0 | 0.499 |
| composed k = 16 | 7 | 7 | 0 | 6.0 | 0.614 |
| **composed k = 32 (primary)** | **7** | **8** | **−1** | 5.0 | 0.710 |
| composed k = 64 | 9 | 7 | +2 | 6.0 | 0.384 |

**Nine grading rules, nine failures to beat a matched random gene set.** The best point estimate in
the table is +2 of 24 at p = 0.38. On magnitudes rather than counts the same holds: a real held-out
pathway's `|r|` exceeds its own matched control's in **10 of 24** pairs under the single-axis
readout (median 0.0935 vs 0.1127, Wilcoxon p = 0.264) and **13 of 24** under the composed readout
(median 0.1959 vs 0.2036, p = 0.790).

**3. So the composed readout's 2 → 7 rise is capacity, not biology, and the run says so through the
control I predeclared for exactly this.** The `random_control` answer rate rises monotonically with
capacity — 17.8% at k = 1 to 36.7% at k = 32 — which is the **first** listed distrust condition of
§1.6 ("the `random_control` answer rate rises with `k` at all"). It fired, and it fires against the
(B) reading, not for it.

**4. The measured verdict is not the one I predicted, and it is worse for my instrument than the one
I predicted.** I predicted composed k = 32 would answer 2–4 of 24 → a clean (A) by count. It answered
**7**, which is the predeclared **INTERMEDIATE** band (5–7), so **by the letter of §1.5 the count
rule returns INTERMEDIATE, not (A)** — and that is recorded as a missed prediction, not smoothed
over. But (B) is excluded outright: it required ≥ 8 *and* ≤ 9 of 90 random controls *and* a paired
win, and it got 7, 33 and −1. §1.5's own escape clause therefore applies —
*"the reported verdict is 'the fairness control failed' and no gain is claimed"* — and the paired
control decides the substance: **there is no readout gain that survives its own control, so the
representation is not withholding a distributed signal about untrained pathways. It is (A).**

**5. The floor half of my fairness control did much less work than I claimed it would, and the null
half did almost none.** I predicted the capacity-matched floor would come out higher than the
single-axis matched floor on ≥ 80% of targets. Measured: **43.3%**, with both medians at 0.200 and
almost identical NaN counts. And the composed readout's own permutation null is *lower* than the
single-axis one (median p95 0.0756 vs 0.0824) while its observed statistic is **2.1× larger** on
`heldout_pathway`. So neither matched instrument absorbs the capacity inflation. **The thing that
does absorb it is the matched random gene set** — because a random gene-set score is not a random
vector: it is a real molecular measurement with real covariance with the image, so it inflates
exactly as much as a real pathway does. That is the methodological finding of Test 1 and it
generalises past this run: *a permutation null answers "is there any patient-level association",
and it is the wrong question. The question a promptable interface has to answer is "is the
association specific to this pathway", and only a matched-gene-set control asks it.*

**6. The causal-name bridge fails, on both adjacency definitions, and the legibility-matched
uncertified axes do it slightly better.** P3's 29 certified names do not predict which untrained
targets their axes can read.

---

# TEST 1 — the composed readout

## 1.1 Setting

`runs/d2_final/artifacts/d2_h_seed42.npz`, state `wsi_biology`, partition `test`,
`--adjustment inductive`, discovery fraction 0.5, seed 42, `min_site_count = 10`. Exposure fold
**n = 1,382**, 57 design columns, operator reference digest `2060a635fa83756a…` — **the identical
prepared state as `p4_inductive_adjustment_measured_20260804T2300Z.md`**, produced by calling
`p4_certify.prepare_state` unchanged. Condition 3 relaxed throughout, as in that entry, because at
full policy every cell is 0 and no comparison exists.

**Reproduction check, and it passes.** Under the published grading the incumbent answers
**28 of the 90 non-control targets** — 18 `hallmark_in_training` / **2** `heldout_pathway` /
5 `immune_tme` / 3 `tumour_state`. That is the measured entry's table to the digit. The run is
graded against a reproduced incumbent, not a re-derived one.

## 1.2 The three fairness controls, and what each turned out to do

**(i) Matched permutation null.** Each arm's null is that arm's own readout on within-cancer
permuted targets (`within_stratum_permutations`, 200 draws, seed 42); for the composed arm the
whitening and the canonical directions are refit inside every draw. **Verdict: it does not control
capacity.** Composed null p95 medians are *below* single-axis ones on every group
(`heldout_pathway` 0.0755 vs 0.0838; `random_control` 0.0760 vs 0.0833) while composed observed
statistics are 1.3–2.1× larger. Destroying the patient pairing destroys the composed readout's
train→test direction transfer as well, so the null collapses rather than widening.

**(ii) Matched detection floor, identical planted signal.** `spike_targets` plants the spike into
the target along `image_direction = e_a` — the incumbent's own supporting axis, as a unit vector in
the full 256-column space — so **both arms receive byte-identical spiked data** and differ only in
the readout. The recovery matrix is turned into a floor by `floors_from_recovery`, the same 80%-of-
draws rule, the same level grid `(0, .01, .02, .05, .10, .20, .40)`, the same 25 draws.
**Verdict: it barely discriminates.** Composed floor > single matched floor on 43.3% of gradeable
targets; medians identical at 0.200; 150 vs 151 of 180 targets resolve a finite floor. The floors
are also *coarse* — a seven-level grid on which most targets land at 0.20 or 0.40 — which is a
known property of the shipped instrument and is why the floor cannot carry this comparison alone.

**(iii) Matched biological null — the paired random gene set.** Verdict: **this is the control that
works**, and it is the one no previous P4 run included. See "Bad news first" §2.

**A fourth control, run because the incumbent's published null does not charge for axis selection.**
The supporting axis is an `argmax` over 256 axes taken on the rows it is then scored on, and the
published null holds it fixed. The selection-aware null recomputes the whole 256-axis grid inside
every one of the 200 permutations and takes the same `argmax`. It raises the bar (median null p95 on
`heldout_pathway` 0.0926 vs 0.0838) and changes almost nothing: `heldout_pathway` **2 → 2**,
`random_control` 21 → 21, `immune_tme` 5 → 5, `tumour_state` 3 → 2. **The published null was
generous, and by an amount that does not matter.** (Prediction 6 held: 0–2 predicted, 2 measured.
`hallmark_in_training` was excluded from this arm on cost and is reported as ungraded, not as
failed.)

## 1.3 The capacity ladder in full

| k | `heldout_pathway` /24 | `random_control` /90 | `hallmark_in_training` /50 | `immune_tme` /8 | `tumour_state` /8 |
|---:|---:|---:|---:|---:|---:|
| incumbent (single axis) | 2 (8.3%) | 21 (23.3%) | 18 (36.0%) | 5 | 3 |
| 1 | 3 (12.5%) | 16 (17.8%) | 10 (20.0%) | 4 | 1 |
| 2 | 4 (16.7%) | 26 (28.9%) | 18 (36.0%) | 5 | 4 |
| 4 | 7 (29.2%) | 30 (33.3%) | 24 (48.0%) | 5 | 3 |
| 8 | 9 (37.5%) | 29 (32.2%) | 25 (50.0%) | 5 | 4 |
| 16 | 7 (29.2%) | 29 (32.2%) | 24 (48.0%) | 5 | 4 |
| **32 (primary)** | **7 (29.2%)** | **33 (36.7%)** | **28 (56.0%)** | **5** | **4** |
| 64 | 9 (37.5%) | 33 (36.7%) | 25 (50.0%) | 5 | 5 |

**The `heldout_pathway` column never separates from the `random_control` column at any capacity.**
The `hallmark_in_training` column does — 56.0% against 36.7% at k = 32 — which is the supervised
channel doing what it was trained to do, and it is also the **sanity gate** of §1.7 (≥ 25 of 50
required for an (A) to be reportable at all; 28 measured, gate passes, prediction 4 held). The
`immune_tme` column is flat at 5 of 8 across every rule tested, incumbent included: that channel is
real, it is not composed out of anything, and one axis already carries it.

## 1.4 What the composed readout actually bought, target by target

Five of the 24 held-out pathways are newly answered at k = 32. Their incumbent margins to their own
bar were −0.005, −0.024, −0.047, −0.120, −0.137 — so two were near-misses and three were not, which
means §1.6's "the gain is only tipped near-misses" distrust condition does **not** fire. The gain is
real *as a gain in the statistic*: median `|r|` on `heldout_pathway` goes 0.0935 → 0.1959 (2.10×).
It simply is not specific: the same readout takes `random_control` from 0.1550 → 0.2436 (1.57×), and
`random_control`'s single-axis median (0.1550) was **already higher** than `heldout_pathway`'s
(0.0935) before any composition.

The strongest composed reads among the 24 are `..._NNK_NNN_TO_CHRNA7_E2F_SIGNALING_PATHWAY` (0.322,
floor unresolvable), `..._BENZO_A_PYRENRE_TO_CYP_MEDIATED_METABOLISM` (0.274, floor unresolvable) and
`..._EBV_BARF1_TO_INTRINSIC_APOPTOTIC_PATHWAY` (0.272 against a floor of 0.40). **The three largest
associations in the untrained set are all refused for want of a floor, not for want of a
correlation** — the same pathology the 20:00 entry recorded, unmoved by composition.

## 1.5 Grading against §1.5 and §1.7 verbatim

* Composed k = 32 = **7 of 24** → predeclared band **INTERMEDIATE (5–7)**. Recorded as such.
* (B) clause 1 (≥ 8): **FAIL** (7). Clause 2 (≤ 9 of 90 random controls): **FAIL** (33). Clause 3
  (paired win above the label-flip null p95): **FAIL** (−1, p = 0.710). **(B) is excluded on all
  three.**
* §1.7's (A) sanity gate (≥ 25 of 50 `hallmark_in_training` under the composed readout): **PASS**
  (28), so the readout is demonstrably able to see a channel that exists.
* §1.7's floor-resolution downgrade (> 12 of 24 with a NaN composed floor): **does not fire** —
  5 of 24. The (A) reading is not "not measurable"; it is measured.
* §1.6's first distrust condition (`random_control` rate rises with k): **FIRES**, against (B).

**The substantive verdict is (A).** The count rule lands one target inside INTERMEDIATE and that is
reported as a missed prediction; every control that distinguishes signal from capacity says the
same thing, which is that composing existing directions does not recover an untrained-pathway
channel because there is no untrained-pathway channel to recover. **Scaling the query layer will
not fix a 2-in-24. Neither will scaling it 32-fold, which is what this measured.**

---

# TEST 2 — the causal names do not generalise

## 2.1 Setting and statistic

`p1_evidence/inputs/pca_targets.npz` (128 PCA axis scores), `d2_h_seed42::wsi_biology`, partition
`test`, **n = 2,766**, transductively adjusted — the footing `causal_attribution._axis_legibility`
uses, so the `legibility__d2_h_seed42__wsi_biology` column of `axis_attribution.csv` is the matching
variable and no second legibility was computed. Split 1,383 train / 1,383 test, seed 42.

For axis j, `heldout_cca_projection` gives `px`, the image projected onto the canonical direction
**fitted for axis j on the train half only**; `r(j, t) = paired_absolute_correlation(px,
target_t[test])` over the **40 untrained targets** (24 `heldout_pathway` + 8 `immune_tme` +
8 `tumour_state`).

**Sanity check on the readout:** `|corr(px, py)|` — the axis reading *itself* — has Spearman
**+0.717** against the published per-axis legibility over all 128 axes. The channel this test scores
is the channel P3 measured.

## 2.2 The matched-legibility control passed its balance bars comfortably

1:1 nearest neighbour on legibility, without replacement, caliper 0.02, greedy in descending
certified legibility: **28 pairs, 1 certified axis dropped** (`PCA_098`, the only certified axis with
negative legibility, −0.0038, with no uncertified partner inside the caliper),
**mean |Δ legibility| = 0.00114, max = 0.00357** against the predeclared bars of ≤ 0.005 and ≤ 0.02.
Prediction 10 held (< 3 dropped).

## 2.3 Result — the bridge fails on both adjacency definitions

Spearman on the double-centred 28 × N matrices (row means, then column means, then grand mean added
back), permutation null shuffling adjacency across targets **within each axis**, 10,000 draws,
which preserves both main effects exactly:

| adjacency | arm | targets | ρ | null p95 | permutation p | leave-one-axis-out range |
|---|---|---:|---:|---:|---:|---|
| **A1** (perturbation-response cosine) | **certified** | 36 | **+0.0129** | 0.0522 | **0.334** | −0.005 … +0.035 |
| A1 | matched uncertified | 36 | **+0.0693** | 0.0545 | **0.018** | +0.039 … +0.083 |
| **A2** (predeclared process-family map) | **certified** | 40 | **−0.0210** | 0.0408 | **0.844** | −0.035 … −0.008 |
| A2 | matched uncertified (name-transplant) | 40 | −0.0114 | 0.0283 | 0.543 | −0.026 … +0.002 |

Difference certified − matched uncertified, 2,000-draw bootstrap over resampled axes:
**A1 −0.056 [−0.167, +0.058]**, CI includes zero; **A2 −0.010 [−0.091, +0.080]**, CI includes zero.

**Every one of §2.6's three conditions fails on A1 and on A2**: ρ is not ≥ +0.15, p is not < 0.05,
and the certified arm does not exceed the matched uncertified arm. **§2.7's verdict applies: the
bridge fails.** Predictions 7 (ρ in [−0.05, +0.10], p ≥ 0.05) and 8 (arms within 0.05 — measured
0.056 apart on A1, 0.010 on A2) both held. Prediction 9 held: A1-lit — literal hypergeometric
overlap of the ten named genes with the target gene set — reaches p < 0.05 in **5 of 1,044** certified
pairs and is declared **underpowered**, as the predeclaration said it would be, and is therefore not
read as evidence either way.

**The uncertified arm reads nominally higher on A1 and its own permutation p is 0.018.** Its CI
against the certified arm includes zero, so no claim is made that uncertified names transfer
*better*. What cannot be said is that certification selects for names that transfer: it does not.

## 2.4 The single most legible fact in Test 2

The strongest reads of untrained targets by certified axes are, in order:

| certified axis | its causal name | target it reads best | r | adjacent under A2? |
|---|---|---|---:|:---:|
| PCA_013 | SSBP1, PNPT1, LRPPRC, TFAM, LONP1 — **mitochondrial nucleoid** | `immune_cytolytic_activity` | **0.336** | **no** |
| PCA_013 | (same) | `immune_t_cell_inflammation` | 0.328 | no |
| PCA_013 | (same) | `immune_ifng` | 0.320 | no |
| PCA_032 | NEDD8, CSE1L, TPX2, TACC3, RFC3 — mitosis / neddylation | `immune_ifng` | 0.277 | no |
| PCA_051 | SRP68, SRP72, PHB, TIMM23B — SRP / import | `immune_t_cell_inflammation` | 0.230 | no |

**The ten strongest certified-axis reads in the whole 28 × 40 grid are all immune targets, and not
one certified axis is named for anything immune.** Across the certified arm the mean `r` on
A2-adjacent pairs is **0.0460** and on non-adjacent pairs **0.0542** — adjacency is, before any
centring, mildly *anti*-predictive. Per-family means of (adjacent − non-adjacent) are
SECR +0.019, RIBO +0.003, REPL −0.004, MITO −0.008, TXN −0.022: no family carries the effect, so
§2.8's "driven by MITO↔MITO alone" distrust check has nothing to fire on.

**What an axis is named for and what an axis reads are, on this representation, unrelated.** The
mitochondrial-nucleoid axis is the best *immune* reader in the certified set.

## 2.5 One departure from the predeclaration, recorded as a departure

§2.3's A2 process-family map assigns a family only to the 29 axes P3 certified, so the
legibility-matched **uncertified** arm had no A2 adjacency at all and came back degenerate on the
first run. The addition, made before the full run and labelled in the output JSON
(`a2_control_note`): the control axis is scored against **its certified partner's** adjacency
vector — a name-transplant control, which asks whether a name predicts readout for the axis that
earned it but not for a transplanted one of equal legibility. A1 needed none of this (every axis has
its own top atoms and therefore its own A1 adjacency, which is §2.5 verbatim), **and A1 is the
predeclared primary**, so the verdict does not rest on the addition.

Four of the 40 targets have no usable gene set in the aligned universe and are excluded from A1 (36
remain): `KEGG_MEDICUS_ENV_FACTOR_NNK_TO_DNA_ADDUCTS`, `immune_cytolytic_activity`,
`immune_t_cell_inflammation`, `stroma_caf`. A2 uses all 40.

---

## 3. How the ten predictions did

| # | predicted | measured | |
|---|---|---|:---:|
| 1 | incumbent reproduces 2 of 24 `heldout_pathway` | **2** (and 28 of 90 overall) | ✓ |
| 2 | composed k = 32 answers 2–4 of 24 → (A) | **7** → INTERMEDIATE band | ✗ |
| 3 | composed k = 32 answers ≤ 6 of 90 `random_control` | **33** | ✗ **badly** |
| 4 | composed clears the `hallmark_in_training` sanity gate (≥ 25/50) | **28** | ✓ |
| 5 | composed matched floors higher than single on ≥ 80% of targets | **43.3%** | ✗ |
| 6 | incumbent under the selection-aware null answers 0–2 of 24 | **2** | ✓ |
| 7 | Test 2 ρ_certified on A1 in [−0.05, +0.10], p ≥ 0.05 → bridge fails | **+0.0129, p = 0.334** | ✓ |
| 8 | ρ_matched_uncertified within 0.05 of ρ_certified | **0.056 apart** (A1); 0.010 (A2) | ~ |
| 9 | A1-lit declared underpowered (< 10 significant of 1,160) | **5 of 1,044** | ✓ |
| 10 | fewer than 3 certified axes dropped by the caliper | **1** (`PCA_098`) | ✓ |

**Three clear misses, and they cluster.** I over-trusted my own capacity controls: I predicted the
matched floor would penalise the composed readout (it barely does, #5), and I therefore predicted
the composed count would stay low (it did not, #2) and that random controls would stay quiet (they
did not, #3). The prediction the conclusion rests on — that there is no untrained-pathway channel to
find — was not the one I graded it by, and it survived a control I would not have thought to build
if #3 had come out as predicted.

---

## 4. What this changes, and the exact prose that should change

No edit is made here to `claim_evidence.json`, `claim_guards.py` or any paper file. Named for their
owners:

1. **`paper/P3_P4_PLAN.md` §9.2 (the "28 of 90" table and the paragraph beneath it, ~L436–450).**
   The sentence *"Of the 24 genuinely untrained `heldout_pathway` targets, one survives"* is correct
   and is now measured at **2** on the real inductive operator (already recorded in
   `p4_inductive_adjustment_measured_20260804T2300Z.md`). What must be **added** is the control:
   *the same rule answers 21 of 90 (23.3%) size- and expression-variance-matched random gene sets,
   and 6 of the 24 random sets matched to those very pathways — so the untrained-pathway rate is
   below the noise rate and the two are not separable (paired label-flip p = 0.98).* Until that
   sentence is in §9.2, the "28 of 90" table overstates what the interface can do, because it is
   quoted without the only control that makes an answer count interpretable.
2. **`paper/P3_P4_PLAN.md` §9.2, the Figure 1 description.** P4's figure should carry a
   **random-control band** — the same panel, same sorting, with the matched random gene sets' floor
   and correlation distribution drawn behind the real targets. The refusal-reason colouring already
   planned is necessary and no longer sufficient.
3. **`paper/P3_P4_PLAN.md` §10, gate item 4 (the abstention curve).** The curve must be measured
   **against the matched random-control set as well as against the real targets**, because a
   refusal rule that answers noise and real biology at the same rate has an abstention curve that
   looks identical to a working one. This is a change to the gate's definition, not to its
   threshold.
4. **`paper/P3_P4_PLAN.md` §11, blocker 5 (`composition_attribution` undischarged).** The blocker
   is understated. It currently says a certified axis may not be *named*. Test 2 shows something
   stronger and independent: **for the 29 axes that do carry a certified causal name, the name does
   not predict what the axis reads** (ρ = +0.013, p = 0.334, and the legibility-matched uncertified
   axes read nominally higher). So "certify the name, then use the name to route queries" is not
   available as a design, and the blocker should say that a name is a label on an axis, not a prior
   over what that axis can answer.
5. **`paper/P3_P4_PLAN.md` §2–§3 (what P3 becomes).** The attribution result
   (`post_pbs_constructions_result_20260804T2300Z.md`, 29 of 128) is not yet in the plan. When it is
   added, it must be added with this entry's negative attached: the names are chemically coherent
   and null-controlled, **and they do not transfer to queries.** "The perturbation atlas annotates a
   representation" survives; "the annotation tells you what the axis can be asked" does not.
6. **A methodological note that belongs wherever CALIBRA's floors and nulls are described
   (`paper/P1_CALIBRA_DRAFT.md`, the detection-floor section).** On this cohort neither the
   detection floor nor the within-stratum permutation null controls readout capacity: a 32-component
   held-out CCA reads 2.1× the single-axis statistic while its own permutation null reads *lower*
   and its own direction-matched floor is higher on only 43% of targets. **A matched-gene-set
   control does control it.** This is a limitation of the two shipped instruments, measured, and it
   is separate from the three limits already on record (first-moment-only adjustment, classifier
   family, fitting rows).

---

## 5. Suite

Run on `~/ws_p4c2` at the final commit, thread-capped:
`pytest morpheus/v2/tests morpheus/tests --ignore=morpheus/v2/tests/test_p2_figures.py -q` →
**652 passed in 71.90s** (0 failed).
`pytest morpheus/v2/tests/test_p2_figures.py -q` → **1 passed, 27 errors in 2.61s**, every error
`ModuleNotFoundError: No module named 'matplotlib'` — the known condition of `~/venv`. **Nothing was
installed into that environment.**
`pytest morpheus/v2/tests/test_effective_rank_canonical.py -q` → **13 passed in 4.45s**, including
`test_no_second_definition_exists_in_the_tree`, **with no allowlist entry added for either new
module**.

## 6. Files / provenance

Drivers `v2/research/rebase/nature/p4_certification/composed_readout.py` and
`causal_name_bridge.py` (commits `78a4278`, `9a1f9a7`, `c20ab1a`). Neither defines a statistic.
Outputs `~/p4c_out/{T1_composed_inductive_f50.json, T2_causal_name_bridge.json}`, vendored into
`v2/research/rebase/nature/p4_certification/out/` and persisted to NFS at
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p4_composed_readout/results/`
together with the run script and both logs.

Test 1: n = 1,382, 180 targets, 256 axes, 200 within-cancer permutations, 25 spike draws per
(target, arm), 7 capacity rungs, ~42 min wall on 20 CPU workers.
Test 2: n = 2,766, 128 axes × 40 untrained targets, 10,000 permutations, 2,000 bootstrap draws,
57 s wall. Perturbation matrix `K562_gwps_normalized_bulk_01.h5ad`, `tcga_sd` scaling, aligned
through `load_aligned_response` exactly as `causal_attribution.attribution_report` aligns it.

Related: [[PREDECLARED_p4_composed_readout_and_causal_name_bridge_20260805T0650Z]],
[[p4_certification_end_to_end_20260804T2000Z]],
[[p4_inductive_adjustment_measured_20260804T2300Z]],
[[post_pbs_constructions_result_20260804T2300Z]]
