# An inductive confound-adjustment operator: a new slide can now be adjusted

**Code:** `v2/calibra/inductive_adjustment.py`, with a small refactor of
`v2/calibra/residualise.py`.
**Tests:** `v2/tests/test_inductive_adjustment.py`, 17 tests, all passing.
**Suite:** 534 passed on a verified workspace (517 collected at the last pre-existing commit
`9b9d070`, so the delta is exactly these 17 and no existing test changed).

---

## 1. What was blocking

`cross_fitted_residuals` is **transductive**, and in three separate places, not one:

1. the K nuisance models are fit inside the call and discarded;
2. which model produced a given row's residual depends on
   `KFold(shuffle=True, random_state=seed).split(matrix)`, i.e. on `len(matrix)` — undefined for a
   row that was not in the cohort;
3. the return is `residual - residual.mean(axis=0)`, a **cohort mean**, so adding one row moves every
   other row's coordinates.

Consequences: no adjusted coordinates can be produced for a new slide, so no claim about a deployable
system survives; and two cohorts adjusted separately end up in two different coordinate systems, so
ALCHEMIST and TCGA could never be compared on the same footing.

Point 3 is asserted directly rather than merely described —
`test_transductive_residualiser_does_NOT_have_that_property` runs the transductive function on a
block and on a 40-row subset of that block and shows the answers differ.

## 2. The operator

`ConfoundAdjustmentOperator.fit(matrix, frame, columns, patient_ids=..., site_column=...)` persists:

* a **`DesignSpec`** — the frozen `confound_design`. Two things in that function are cohort
  statistics and had to be frozen or new rows land on a different scale: a categorical's one-hot
  levels (`pd.get_dummies` sorts whatever it is given, so a new cohort produces different columns in
  a different order — handled by `reindex` onto the reference's columns), and a numeric column's
  **mean and s.d.**, which `confound_design` takes from the frame it is handed. Standardising a new
  slide by *its own cohort's* mean is precisely the cross-cohort error the operator exists to
  prevent; a test asserts the reference's purity mean is used and the new cohort's is not. The
  missingness-indicator column is likewise emitted iff the *reference* column had missing values.
* a **`SitePooling`** — the frozen `frequent` site set and `min_site_count`.
* an **`InductiveResidualiser`** — the K ridge models (coefficients + intercepts), the per-row fold
  map, and the reference residual centre.
* **provenance** — reference cohort digest (SHA-256 over the design, the patient ids and the matrix
  shape), separate design and patient digests, the design column names, the pooling threshold and
  frequent set, the residualiser's `n_splits`/`alpha`/`seed`, numpy and sklearn versions, and a
  creation timestamp. `adjust()` returns a report whose `applied_by` block carries the reference
  digest, so an applied adjustment can always name what fitted it.

Persistence is a single `.npz` (arrays plus a JSON `meta` string), `save`/`load`, no pickle.

## 3. Identity on the fitting cohort

`adjust_reference` is **bit-identical** to `cross_fitted_residuals`, asserted with `np.array_equal`
(not `allclose`) at `(n_splits, alpha, seed)` of `(3, 1.0, 7)`, `(5, 0.5, 42)` and `(10, 2.0, 1)`,
plus the zero-design degenerate branch. Bit-identity is the requirement because every published
CALIBRA number came out of the transductive function and a merely-close operator would silently make
all of them non-comparable.

Two implementation details are load-bearing and are commented as such in the source:

* the stored `centre` is subtracted rather than recomputed (it is the cohort mean of point 3 above);
* predictions go back through `Ridge.predict` on rebuilt estimators rather than a hand-written
  `design @ coef.T + intercept`, so the arithmetic path is the same one that produced the original
  numbers.

**A deliberate discontinuity, asserted rather than hidden.** A reference row keeps the one model that
did *not* see it (cross-fitting). A new row was seen by none of them and gets the **mean of the K
fold models' predictions**. The two therefore do not agree exactly, so a slide's coordinates depend
on whether it was in the reference cohort;
`test_reference_and_new_row_paths_differ_and_that_is_documented` asserts both the inequality and that
the two stay in the same coordinate system (correlation > 0.95). The alternatives were rejected
explicitly: pushing reference rows through the new-row path breaks identity and invalidates the
published numbers; refitting a single model on all rows introduces a model no reference row was ever
scored by and which no identity test can validate. The K-model ensemble is the only option that both
preserves identity and uses exclusively models that are already fitted and persisted.

## 4. The unseen-site policy, and why it is not a new rule

**An unseen site is mapped to `OTHER`.**

The justification is that this is not a new policy at all. `pooled_tissue_source_site` already sends
every site with fewer than `min_site_count` patients to `OTHER`. A site the reference cohort never
saw has reference count **zero**, and zero is below every threshold, so the reference's own rule
already answers the question. The implementation makes that literal: `pooled_tissue_source_site` is
refactored into `tissue_source_site` (barcode parsing) plus
`apply_pooled_tissue_source_site(patient_ids, frequent)` (the rule), and the transductive path now
calls the same rule the inductive path does, so the two cannot drift apart in how they parse a
barcode. `test_unseen_site_policy_equals_the_existing_pooling_rule_at_count_zero` asserts that
applying the frozen rule to the reference's own barcodes reproduces `pooled_tissue_source_site`
exactly.

The two alternatives were considered and rejected:

* **Refuse.** This would make the operator unusable on exactly the cohorts it exists for — ALCHEMIST
  shares almost no TSS codes with TCGA, so every row would be refused.
* **Nearest site.** This requires a similarity between site codes. TSS codes are registry
  identifiers, not coordinates; there is no metric on them that is not invented.

**An unseen level in any *other* categorical is REFUSED**, raising `UnseenLevelError`. Different
question, different answer: a zero one-hot row is not "no confound", it is the reference's implicit
baseline level, so silently zero-encoding an unknown cancer type would adjust that slide as though it
belonged to whichever level the design dropped. A caller who genuinely wants that must pass
`on_unseen_level="zero"`, and the choice is recorded in the report. Site codes never reach this branch
because pooling has already mapped them.

## 5. Tests

The four required cases, plus the ones that make them meaningful:

| test | what it pins |
|---|---|
| `test_operator_is_numerically_identical_to_cross_fitted_residuals` | bit-identity, `np.array_equal` |
| `test_identity_holds_for_every_residualiser_setting` | identity at 3 `(n_splits, alpha, seed)` settings |
| `test_operator_design_matches_confound_design_on_the_reference` | the frozen design equals `confound_design` |
| `test_no_design_columns_degenerates_to_centring_exactly` | the zero-design branch |
| `test_unseen_site_is_pooled_to_other_and_reported` | unseen site → `OTHER`, named in the report |
| `test_unseen_site_policy_equals_the_existing_pooling_rule_at_count_zero` | the policy *is* the existing rule |
| `test_unseen_cancer_level_is_refused_not_silently_zeroed` | refusal, and the opt-in escape |
| `test_single_new_row_is_adjustable_at_all` | one slide, no cohort |
| `test_one_row_at_a_time_equals_the_whole_block` | **no cohort statistic leaks into the new-row path** |
| `test_transductive_residualiser_does_NOT_have_that_property` | the motivating defect, asserted |
| `test_reference_and_new_row_paths_differ_and_that_is_documented` | the discontinuity, stated not hidden |
| `test_round_trip_through_persistence_is_exact` | save → load → identical on reference and new rows |
| `test_persisted_operator_carries_its_provenance` | digest, columns, threshold, versions survive |
| `test_reference_digest_separates_cohorts` | the digest identifies the cohort, not the run |
| `test_applied_report_names_the_fitting_cohort` | traceability of an applied adjustment |
| `test_a_second_cohort_can_be_adjusted_with_the_first_cohorts_operator` | the ALCHEMIST case end to end |

`test_one_row_at_a_time_equals_the_whole_block` is the property that distinguishes an inductive
operator from a transductive one and is checked at `atol=0, rtol=0`.

## 6. Not done

* **Not yet wired into any pipeline.** `run_calibra.py`, `p4_certify.py` and the ALCHEMIST comparison
  still call `cross_fitted_residuals` directly. Because `adjust_reference` is bit-identical, swapping
  them over cannot move a number — but the swap is not made here and no existing call site changed.
* **No TCGA operator has been fitted and persisted.** The obvious next step is to fit one on the
  2,530-patient test cohort with the 99-column cancer + pooled-TSS design and store it beside the
  artifacts, so ALCHEMIST can be adjusted with TCGA's operator rather than its own.
* **The fold-ensemble new-row rule is a choice, not a derivation.** It is unbiased in the sense that
  every constituent model is out-of-sample for a new row, but its variance differs from the
  cross-fitted reference path. If that matters for a downstream statistic, it should be measured
  rather than assumed.
