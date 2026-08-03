# PRE-REGISTRATION — D3 (purity sensitivity) and D2.3 (per-axis proliferation / essentiality)

**Logged:** 2026-08-03 13:00 UTC. **This file is written and committed BEFORE either analysis is
run.** Nothing below is fitted to a result. Where a number appears it is either (a) a property of an
input table, computed without touching any outcome, and labelled as such, or (b) a threshold.

**How obtained:** planning only. Inputs inspected: `~/e0_run/d2_v3/d2_v3_s4{2,3,4}/artifacts/`,
`~/e0_run/data/pbs_targets_k128_v2.npz{,.axis_annotations.csv}`,
`~/e0_run/data/tcga_pancan_rna.sample_map.parquet`, and the ABSOLUTE purity table downloaded today
(provenance below). No representation state has been correlated with any target yet.

---

## 0. The purity table — blocker resolved with a real one, not a proxy

The ledger records "no TCGA consensus purity table on disk", and a previous agent substituted
dx_year/age + TSS pooling. **That substitution is not needed.** The box has outbound network access
and the PanCanAtlas ABSOLUTE calls are openly downloadable:

- **File:** `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` (PanCanAtlas ABSOLUTE purity/ploidy
  master calls; Taylor et al. 2018, *Cancer Cell* 33:676, aneuploidy companion release).
- **URL:** `https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5`
- **sha256:** `f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd`
- **Rows:** 10,786 samples; 10,642 with finite purity; range 0.08–1.00, median 0.65.
- **`purity_source` to be declared:** `absolute`. **Not** `expression_derived`. ABSOLUTE is inferred
  from SNP-array copy number and allelic fractions, i.e. **from DNA, not from the RNA that the
  targets are built from.** The circularity caveat that would attach to an ESTIMATE-style
  expression-derived score therefore **does not apply here**, and the fallback in
  `HANDOFF_PHASE_D.md` §D3 is not exercised.

**Sample selection rule, fixed now.** Purity is a per-*sample* quantity; our cohort is keyed by
12-character patient. Rather than guess, each cohort patient is matched to the **exact RNA source
sample** it was built from via `tcga_pancan_rna.sample_map.parquet`, joined to ABSOLUTE on the
15-character `TCGA-XX-XXXX-NN` key. 59 patients (0.6%) still carry more than one purity-bearing
sample; for these the rule is *lowest sample-type code first* (01 primary before 03/06/02/05), then
lexicographically smallest barcode. The within-patient purity spread for those 59 has median 0.16,
which is recorded as a known limitation rather than smoothed away.

Coverage of the D2 cohort under this rule (input property, no outcome touched):
**test 95.81%**, train 96.18%, val 97.05% — above `run_calibra`'s 80% `purity_coverage` gate.
The analysis is **complete-case**: `run_calibra` intersects the mask, and the before/after
comparison is run on the *identical* patient set, so cohort composition cannot masquerade as a
purity effect.

---

## 1. D3 — what would count as pass / fail

**Run.** `morpheus.v2.calibra.run_calibra --artifacts d2_h_seed42.npz d2_i_seed42.npz --targets
frozen_rna_targets.npz --partition test --levels 0.0,0.01,0.02,0.05,0.10,0.20,0.30,0.40,0.50
--n-draws 40 --n-components 16 --n-permutations 2000 --seed 42 --purity-table <above>
--purity-source absolute --require-rna-positive-control` — deliberately the **same** readout
settings as the 2026-08-03 02:30 UTC T1.2 ledger, so the "before" arm is comparable to a number
already on the ledger.

**Primary statistic.** `excess_over_null_median` for **`wsi_biology`** — the morphology→molecular
channel. It is `observed − null_median`, and it, not `adjusted_top_cca`, is the paired quantity,
because *the null is not zero*: the T1.2 ledger measured the 16-component within-cancer permutation
null at **0.1463–0.1483**, and adding a design column moves that null. Absolute levels are not
quoted across designs; only the paired within-run before/after difference is.

Secondary, reported before/after regardless: `adjusted_top_cca`, `heldout_top_cca`,
`permutation_null_median`, `permutation_null_p95`, `permutation_p`, `attenuation_slope`,
`detection_floor`, `effective_rank`.

**Declared outcomes for `wsi_biology`, decided on the after-purity arm:**

| verdict | condition |
|---|---|
| **FALSIFIER FIRES — channel dies** | after-purity `adjusted_top_cca` ≤ its own `permutation_null_p95`, **or** `excess_over_null_median` retains **< 50%** of its before value, **or** `permutation_p` leaves the 1/2001 resolution floor |
| **PARTIAL — channel is partly tumour content** | `excess_over_null_median` retains 50–80% of before |
| **SURVIVES** | retains ≥ 80% of before, observed still > `null_p95`, `permutation_p` still at floor |

This is escalation item 5 in `HANDOFF_PHASE_D.md` §5. If the falsifier fires it is reported, not
buried, and nothing is built on top of it.

**Mandatory companion — rank-matched placebo.** Adding `purity` adds exactly one column to the
108-column design (complete-case, so `confound_design` emits no missingness indicator). A drop in
the channel could therefore be *design rank increasing*, not *purity being removed* — this is
exactly the confusion T2.2 in `NOTEBOOK.md` exists to prevent. So a second, otherwise identical run
is made with a **placebo table: the same patients, the same 108+1 design rank, the purity values
randomly permuted across patients** (seed 20260803). Interpretation, fixed now:

- placebo delta ≈ 0 and real delta large ⇒ the change is **purity**.
- placebo delta ≈ real delta ⇒ the change is **design rank**, and D3 has measured nothing about
  purity. This reading is binding even if it is the inconvenient one.

`purity_source` for the placebo will be declared as a distinct value so it can never be read as a
purity adjustment; this requires a one-line widening of the `--purity-source` enum plus a test.

**Positive control.** `rna_biology` (RNA→RNA, circular by construction) must clear its own pairing
null in **both** the before and after arms, in the same run — `--require-rna-positive-control`
enforces this. If it fails, the run measured nothing and no D3 verdict is issued.

---

## 2. D2.3 — what would count as pass / fail

**The question.** Are the legible axes just proliferation?
**Ledger falsifier:** *every legible axis coming back proliferation-loaded.*

### 2a. A defect in the pre-built annotation, recorded in advance

The 128 axes are annotated, but `proliferation_loading` as built by `v2/build_pbs_targets.py:107` is
the **|loading|-weighted mean over all 7,072 basis genes** of a binary proliferation flag. The SVD
basis is dense, so those weights are near-uniform and **every axis is compressed towards the
background rate**. Measured on the annotation alone, no outcome involved:

- background proliferation fraction over the 7,072 basis genes: **0.0841** (576 genes).
- per-axis `proliferation_loading`: min 0.0738, median 0.0820, max 0.1453 — i.e. the *median axis
  sits at 0.98× background* and the whole spread is 0.85× background.
- the same axes scored as *proliferation fraction among their top-100 |loading| genes*: min 0.021,
  median 0.101, max 0.378 — spread **4.24× background**.
- the two statistics agree only at **Spearman 0.577**.

So the pre-built column is a **weak instrument**: it is diluted almost to background and it does not
rank axes the same way a concentrated statistic does. This is stated now, before any legibility
number exists, so that it cannot later be produced as an excuse for whichever result arrives.

**Consequence, fixed now.** Both statistics are carried through the whole analysis:
- `prol_wmean` — the pre-built `proliferation_loading`, **primary**, because it is the column the
  ledger names.
- `prol_top100` — proliferation fraction among the axis's top-100 |loading| genes, **co-primary**,
  because the primary is demonstrably diluted.
A verdict is only called "clean" if the two agree. Disagreement is itself reported.
`essentiality_loading` (and its top-100 analogue) is carried alongside throughout.

### 2b. Legibility

Per-axis legibility `L_a` = **held-out single-direction correlation** between `wsi_biology` and PBS
axis *a*, cross-fitted 5-fold ridge (`score_target_block_per_column`, the same code path CALIBRA
already uses), on the **test split**, both sides residualised on the same cancer+TSS design by
`cross_fitted_residuals`. Primary artifact is **`d2_h_seed42`** (Hallmark-supervised), because it
was **not** trained on the PBS axes and so its legibility of them is not a restatement of its own
supervision. `d2_i_seed42` (PBS-supervised) is reported as a secondary, explicitly labelled circular.
Seeds 43 and 44 are run for stability; only paired within-seed quantities are quoted.

**Null.** Per-axis legibility has a non-zero null for the same capacity reason the 16-component CCA
does. Null = the identical pipeline with patient labels permuted **within cancer stratum**, 200
permutations, giving a per-axis null distribution. "**Legible**" = `L_a` above the axis's own
within-cancer permutation null p95. The row-shuffle null is *not* substituted for it.

**Proliferation-loaded** = axis in the **top quartile** of the proliferation statistic (n=32 of 128).
Chosen as a distributional threshold precisely so that it needs no knowledge of the outcome.

**Declared outcomes:**

| verdict | condition |
|---|---|
| **FALSIFIER FIRES — legibility is proliferation** | ≥ 90% of legible axes are proliferation-loaded (top quartile), **and** Spearman(`L`, proliferation) is positive with a 95% CI excluding 0 |
| **PARTIAL** | 50–90% of legible axes proliferation-loaded, or Spearman CI excludes 0 with \|ρ\| ≥ 0.3 |
| **DISCHARGED — legibility is not proliferation** | < 50% of legible axes proliferation-loaded (chance is 25%), Spearman 95% CI includes 0 or \|ρ\| < 0.3, **and** legible non-proliferation axes exist in number |

**Distribution, not a summary.** Required in the report regardless of verdict: the full 128-axis
legibility distribution; a bimodality check (dip statistic + the gap between the top and bottom
mode); the share of total legibility carried by the top 5 and top 10 axes; and an explicit statement
if legibility is dominated by a few axes. A single headline correlation is not an acceptable answer.

**Guard-test consequence, acknowledged in advance.** `v2/calibra/claim_guards.py` pins
`proliferation_deflation` as an undischarged blocker on E0's `transfer` claim, and
`tests/test_claim_guards.py::test_current_e0_result_is_not_yet_an_admissible_transfer_claim`
asserts it. **If and only if the DISCHARGED verdict is reached**, that test is updated deliberately,
with the reasoning in the commit message — not worked around, and not touched at all under the
FALSIFIER or PARTIAL verdicts. Note also that discharging `proliferation_deflation` alone does not
make E0 admissible: `single_platform` remains, and is not in scope here.

---

## 3. Things that would invalidate either analysis, listed now

1. `rna_positive_control_passed` false, or any `channel_gate_failures` in the D3 run → no verdict.
2. `purity_coverage` < 0.80 → no verdict (currently 0.958 on test).
3. A degenerate permutation null (`permutation_null_non_degenerate` false) in either arm → no verdict.
4. D2.3: if the within-cancer permutation null for legibility is degenerate, or if fewer than 10
   axes are legible, the falsifier test has no power and that is reported as "no verdict", not as
   a discharge.

### Files / commits

- Purity table: `~/e0_run/d3/purity/{abs.txt,tcga_absolute_purity_by_patient.csv,placebo_shuffled_purity_by_patient.csv}`
- Annotation sharpness check: `~/e0_run/d3/axis_annotation_sharpness.csv`
- This pre-registration, committed before any run.
