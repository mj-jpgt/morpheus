## 2026-08-04 06:40 UTC — PREDECLARATION: is D2's headline an artifact of the evaluation's coordinate system?

**Written before any of the three tests below was run.** Nothing in this file is a result.

**Workspace.** `~/ws_d2sym/morpheus` on the A100 (`150.136.45.194`), built from
`git -c core.autocrlf=false archive HEAD` at commit **80a14d6** and verified **452/452 files by git
blob SHA-1**, zero mismatched, zero missing, zero extra on disk. HEAD moved to `064baea` while this
workspace was being built (a concurrent agent added P2 figure scripts);
`git diff --name-only 80a14d6 064baea -- v2/research/rebase/d2_compare.py v2/calibra/
v2/paired_bootstrap.py v2/survival_evaluation.py v2/tests/` is **empty**, so every module these tests
import is identical at both commits. CPU only, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`; a GPU chain is running and is not contended for.

**No statistic is computed inline.** Every rank/channel number comes from
`morpheus.v2.calibra.spectral`, `morpheus.v2.calibra.residualise`, `morpheus.v2.paired_bootstrap` or
`morpheus.v2.survival_evaluation`, imported. This is the discipline that caught three statistic
substitutions on this project.

### The confound

D2's published headline is that supervising with a 128-D Perturb-seq dictionary (arm I, "PBS") gives a
worse morphology→molecular channel than supervising with 50 curated Hallmark gene sets (arm H):
Δ = −0.1325 / −0.1089 / −0.1226 on seeds 42/43/44, on the "untrained 40"
(`heldout_pathway` + `immune_tme` + `tumour_state`), both bootstrap CIs excluding zero 3/3.

What was verified: those 40 targets do not **overlap** arm H's 50 training targets. What was never
verified: that they are the same **kind** of object. All 180 rows of `frozen_rna_targets.npz` are gene-set
scores. The untrained 40 are gene-set scores. Arm H is trained to produce gene-set scores. Arm I is
trained to produce Perturb-seq dictionary codes — a different space. So the exam may be written in
Hallmark's coordinate system, and "PBS is a worse supervision target" may reduce to "PBS was marked in
someone else's units".

Confirmed before predeclaring: arm I's supervision file is `pbs_targets_k128_v2.npz`,
sha256 `4f7d6f409988a8191bd41a84cea0a2e12096ac50f41bfe287e72ef6f7e40fd40`, recorded as
`pbs_target_sha256` in `D2_PAIR_MANIFEST.json` for all three seeds and byte-identical to
`~/e0_run/data/pbs_targets_k128_v2.npz`. The symmetric exam below therefore uses **exactly** the codes
arm I was trained on.

---

### Test 1 — the symmetric evaluation

Both arms scored against the **128 PBS dictionary codes** as targets, through
`morpheus.v2.research.rebase.d2_compare` **unmodified**: cancer + pooled-TSS residualisation
(`min_site_count=10`), `top_canonical_correlation` at 16 components, paired patient and cancer-cluster
bootstrap at 2,000 repeats, `--seed 42`, held-out `test` partition, the same three seed-matched
artifact pairs from `~/e0_run/d2_v3/d2_v3_s4{2,3,4}/artifacts/`. Only `--targets` changes.

An **anchor run** goes first: the published gene-set comparison re-run on this workspace. If it does
not reproduce −0.1325 / −0.1089 / −0.1226 the workspace is wrong and nothing below is admissible.

**The asymmetry is the point and is stated up front.** Arm I trained *on* these 128 codes. This exam is
generous to arm I in exactly the way the gene-set exam is generous to arm H. Neither exam is neutral;
the question is whether the published claim survives when the generosity is pointed the other way.

**Readings, fixed now:**

- **Arm I wins on PBS codes while arm H wins on gene sets** → the result is coordinate-system
  dependent. The headline collapses to "each arm is better in its own space", which is close to
  vacuous. **Report as a refutation of the published claim** and D2_RESULT.md §6 must be withdrawn as
  stated.
- **Arm H wins on both** → the claim survives its strongest available objection, on an exam tilted
  towards its opponent, and is **considerably stronger than it currently reads**.
- **Neither arm separates on PBS codes** (both CIs cover zero, or the two arms are within the
  −0.0099…−0.0280 band the `random_control` negative control already occupies) → **uninformative**.
  Say so. Do not read a null as support for either arm.

"Wins" means: point Δ in the stated direction **and** the patient CI₉₅ excluding zero, reported per
seed with the cancer CI beside it. A 2/3-seed split is reported as a 2/3 split, not rounded to a verdict.

### Test 2 — a space neither arm trained on

Predeclared preference order was: a clinical endpoint (survival, stage, or the pre-registered ER status
of P1's T1.7(b)); then drug response; then `random_control` gene sets as a floor. **What is actually on
disk decides this, and an absence is a finding**, not a prompt to substitute.

Established before predeclaring (inventory only, no statistic computed):

- **ER / PR status is unusable.** `p1_evidence/inputs/tcga_clinical_covariates.parquet` labels 690 ER
  and 688 PR patients; **all are BRCA, and BRCA is a development cancer in the maximal split — 585
  train / 105 val / 0 test.** The D2 comparison lives entirely on the test partition. The
  pre-registered clinical control has **zero** coverage there and cannot be run on this contrast at
  all. `PREREG_known_covariate.json` already declared the in-distribution weakness; the split makes it
  not merely weak but empty.
- **Overall survival is usable.** `runs/v22_a10_11v21_20260725/discovery_inputs/tcga_cdr_outcomes.parquet`
  carries `os_time_days` / `os_event` for **2,765 of the 2,766 test patients** and `pfi_*` for 2,159.
  It is a clinical endpoint, in neither coordinate system, on the exact held-out partition.
- **No drug-response data** ties to these patients.

So Test 2 is **overall survival**, primary; PFI reported beside it. Instrument:
`morpheus.v2.survival_evaluation.paired_cindex_bootstrap` — Harrell C-index, paired
challenger-minus-teacher, patient and cancer bootstrap, 2,000 repeats — on out-of-fold held-out risks
from the same frozen `wsi_biology` blocks, with the same cancer + pooled-TSS residualisation.
Top-CCA is **not** used: it is a multivariate maximum against a 1-column target and is inflated by
capacity alone, which `spectral.heldout_single_direction_correlation`'s docstring already records as a
way to manufacture a channel out of noise.

**Readings, fixed now:**

- **Arm H beats arm I on survival with a CI excluding zero** → the D2 direction is not a coordinate
  artifact; it reproduces on an endpoint in neither training space. Strongest possible support.
- **Arm I beats arm H on survival with a CI excluding zero** → the D2 direction reverses off the
  gene-set exam. Refutation.
- **Neither separates** → the arms are indistinguishable on the one neutral endpoint that exists, and
  the D2 gap is only observable inside one of the two training coordinate systems. This is
  **uninformative about which supervision is better** and must be reported as a limit on what D2 can
  claim, not as support for either arm. Note in advance that C-index differences on a frozen
  representation are low-powered and a null here is weak evidence.

### Test 3 — is the PCA comparison circular?

The second headline (`t11_t12_must_beat_baselines_20260803T0440Z.md`) is that ordinary PCA of the
expression matrix, at matched capacity, beats PBS in 3 of 4 cells and never loses.

**A correction to the brief, recorded before running.** The brief supposes PCA is a *predictor* handed
the targets' basis. It is not: in T1.1 the PCA block and the PBS block are both **targets**, scored
against a fixed frozen representation (`d2_h` / `d2_i`, `wsi_biology` and `full_biology`). So "PCA is
handed the gene-set targets' basis" is not the available failure mode. The failure mode that *is*
available, and that this test measures, is:

`spectral.cca_spectrum` PCA-whitens **both** sides to `n_components` before taking singular values, so
top-CCA at 16 components depends on a 128-column target block **only through that block's top-16
principal subspace**. The PCA block's top-16 PCs are, by construction, the 16 leading variance
directions of the expression matrix. If that is what drives the win, then (a) "capacity-matched at 128"
is not what was compared — 112 of the 128 columns never enter the statistic — and (b) the result says
where expression variance sits, not that an interventional basis is inferior.

Run, on the same `runs/d2_final/artifacts/d2_{h,i}_seed42.npz` the T1.1 table used, with the same
residualisation and a paired bootstrap over target blocks:

1. **Reproduce** the four T1.1 PBS-vs-PCA cells at `n_components = 16`. If they do not reproduce, stop.
2. **Sweep the readout budget** `n_components ∈ {8, 16, 32, 64, 128}` for PBS, PCA and the
   size- and spectrum-matched random dictionary.
3. **Report the share of each block's variance carried by its top-16 PCs**, which is the quantity the
   16-component readout is actually selecting on.

**Readings, fixed now:**

- **PCA's advantage holds across the budget sweep, including at 128** → not a readout artifact. The
  T1.1 verdict stands and the deflationary conclusion about PBS stands with it.
- **PCA's advantage is present at 16 and gone or reversed by 64/128** → the comparison was decided by
  the readout budget, not by capacity. "Capacity-matched at 128" is then **false as written** and the
  T1.1 claim must be withdrawn or restated as budget-conditional.
- **Every block converges as the budget grows** (all differences shrink toward zero together) → the
  statistic saturates and the comparison is uninformative at high budget; the 16-component result is
  then the only readable one and must be quoted with the subspace caveat attached.

### What is not being claimed

None of these tests can decide whether the Perturb-seq dictionary is a good scientific object. They
decide a narrower thing: whether D2's published number measures a property of the supervision or a
property of the exam. A refutation here does not vindicate PBS — it removes a claim.

### Files

Tests write to `~/ws_d2sym/out/` and are copied back under
`v2/research/rebase/nature/d2_coordinate_system/`. `NOTEBOOK.md`, `paper/P2_RANK_DRAFT.md` and
`v2/research/rebase/nature/D2_RESULT.md` are not edited by this work.
