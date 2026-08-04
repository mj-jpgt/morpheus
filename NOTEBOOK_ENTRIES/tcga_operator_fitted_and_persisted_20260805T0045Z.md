# The TCGA confound-adjustment operator, fitted and persisted — in two arms, the second forced by ALCHEMIST

**Code:** `v2/research/rebase/nature/fit_tcga_operator.py`.
**Artifacts:** `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/runs_misc/tcga_operators/`
— 14 `.npz` operators (11 MB), 14 sibling `*.provenance.json`, and
`fit_tcga_operator_report.json`. Persistent storage only.
**Workspace:** `git -c core.autocrlf=false archive HEAD`, 662 files verified by git blob SHA-1, 0
mismatches. Threads capped to 1.

`v2/calibra/inductive_adjustment.ConfoundAdjustmentOperator` existed and was bit-identical to the
transductive path on its reference cohort, but its own notebook entry recorded that **no TCGA
operator had ever been fitted or persisted**, so the thing it was built for — adjusting an external
cohort with *TCGA's* operator rather than its own, which is the only way to put two cohorts in one
coordinate system — was still not possible. It is now.

---

## 1. What was fitted

Reference cohort: the **2,530** `test`-partition patients of
`runs_misc/calibra_run/artifacts/diagnostic_full_seed42.npz` that also carry RNA targets — the
`run_calibra.py` mask reproduced exactly, so the operator is fitted on the same patients every
published CALIBRA number was measured on. 21 cancer levels, 342 raw TSS codes, 75 sites at or above
`min_site_count = 10`. Residualiser `n_splits = 5, alpha = 1.0, seed = 42`.

Seven representation states × two arms = **14 operators**:

| arm | columns | design columns | for |
|---|---|---|---|
| `full` | `["cancer", "tss"]`, site pooled at count 10 | **99** | TCGA-internal work; the design every published CALIBRA number used |
| `cancer_only` | `["cancer"]` | **22** | any comparison with ALCHEMIST |

## 2. Why the second arm exists, and what it costs

**ALCHEMIST has no `tissue_source_site` field at all.** The GDC facet returns `_missing` for all
1,176 cases; `v2/research/external/build_alchemist_manifest.py` writes `tss = "NA"` for every row with
the comment that ALCHEMIST publishes no site and the column is *"explicitly NA rather than a
fabricated code"*; and the paired manifest on the box carries NaN in all 1,106 rows. Worse for the
inductive path specifically: ALCHEMIST barcodes are `ALCH-B0CW-…`, whose second field is a *per-case*
code, so `pooled_tissue_source_site` would collapse every ALCHEMIST patient to a single `OTHER` level
— a constant column, which carries no information and cannot be adjusted for.

So the site half of the design cannot be applied to ALCHEMIST, and it must therefore be dropped from
**both** arms or the two cohorts are not adjusted by the same map. Stated plainly:

> **The external comparison is cancer-adjusted only, and is therefore weaker than TCGA's internal
> protocol.** Site is the confound this project's own leave-sites-out work treats as the one that
> matters most, and in the ALCHEMIST comparison it goes unadjusted on both sides. A `cancer_only`
> result is not evidence that the channel survives site adjustment; it is evidence about a
> cancer-adjusted channel and nothing more. The 99-column `full` operator is the one to quote for
> anything TCGA-internal, and the two must never be reported as if they were the same measurement.

## 3. Verification — all 14 operators, all three checks

| check | result |
|---|---|
| **V1 — identity on the fitting cohort.** `adjust_reference` vs `cross_fitted_residuals` on the same frozen design | **14/14 bit-for-bit** (`np.array_equal`), max abs difference **0.0** in every case |
| **V2 — unseen site.** Rows whose TSS code the reference never saw | **7/7 `full` operators** pool all such rows to `OTHER`, report `sites_pooled_to_other`, and carry the policy string *"unseen site -> OTHER (the reference pooling rule at a count of zero)"*. Adjusted values are finite, and land on **exactly the same coordinates** as an explicitly-`OTHER` row — so the policy is applied, not merely reported. Not applicable to the `cancer_only` arm, which has no site column |
| **V3 — persistence round-trip.** `load(save(op))` | **14/14** reproduce the reference path *and* the new-row path under `np.array_equal`, and carry the provenance (`reference_digest`, `cohort_name`) unchanged |

The unseen-site policy is not new and is not invented for the inductive path: a site the reference
never saw has reference count 0, which is below every `min_site_count`, so
`pooled_tissue_source_site`'s own rule already says `OTHER`. The alternative — refusing such rows —
would make the operator unusable on exactly the cohorts it exists for.

## 4. Provenance persisted

Each `.npz` carries a JSON `meta` block and each has a sibling `*.provenance.json` so the digest is
greppable without loading the array file. Fields: `format_version` (1), `cohort_name`
(`TCGA::diagnostic_full_seed42::<state>::test::<arm>`), `reference_digest` (sha256 over the design,
the patient barcodes and the matrix shape — e.g. `a15853670043c9b6…`), `reference_design_digest`,
`reference_patient_digest`, `n_reference_rows` (2530), `n_representation_columns` (256),
`source_columns`, `design_columns`, `n_design_columns`, `site_pooling` (column, `min_site_count`, the
75 frequent sites by name), `residualiser`, `numpy_version`, `sklearn_version`, `created_utc`. Every
applied adjustment reports `applied_by.reference_digest`, so an adjusted slide can always be traced
to what fitted it.

## 5. A second blocker for the external comparison, found and reported rather than defaulted

**ALCHEMIST's cancer labels are disjoint from TCGA's.** ALCHEMIST carries
`diagnoses.primary_diagnosis` strings (`Adenocarcinoma_NOS` 808, `Squamous_cell_carcinoma_NOS` 209,
`Adenosquamous_carcinoma` 21, …, 16 levels); TCGA carries `LUAD`/`LUSC`. A TCGA operator therefore
**refuses every ALCHEMIST row** with `UnseenLevelError`, in both arms — verified, with the error text
recorded in the report.

That refusal is correct. A zero one-hot row is not "no cancer type", it is the reference's implicit
baseline level, so `on_unseen_level="zero"` would adjust every ALCHEMIST slide as though it belonged
to one arbitrary TCGA cancer type. The `zero` path was exercised and its report recorded, and it is
**not** adopted here. The external comparison needs an explicit, declared label mapping —
`Adenocarcinoma_NOS → LUAD`, `Squamous_cell_carcinoma_NOS → LUSC`, and a stated policy for the 14
remaining levels (≈ 65 cases) — and that mapping is a scientific decision, not a default. It is
flagged for whoever runs the comparison.

## 6. What is still missing before ALCHEMIST can actually be adjusted

* **The ALCHEMIST representation does not exist yet.** Extraction was at 562 of 1,106 slides when this
  was written; `v2/research/external/alchemist_channel.py` would build it but has not been run to
  completion. Nothing here is blocked on that, but the operator cannot be *applied* until it exists.
* **The operator is per representation matrix.** These 14 are fitted on the CALIBRA artifact's seven
  states. The ALCHEMIST comparison uses a different matrix (raw H-optimus, pooled and PCA-256), so
  when that matrix exists the *same* `fit_tcga_operator.py` call must be re-run against it. The
  design, the pooling rule and the verification are unchanged; only the matrix differs.
* **`alchemist_channel.py` still fits PCA within each cohort separately**, which is the transductive
  error this operator exists to remove. Wiring it to the persisted operator is the next step and is
  not done here.
