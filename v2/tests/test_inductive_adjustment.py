"""Tests for the inductive confound-adjustment operator.

The first test is the load-bearing one: if the operator is not numerically
identical to ``cross_fitted_residuals`` on the fitting cohort, every published
CALIBRA number becomes non-comparable and nothing else here matters.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.calibra.inductive_adjustment import (COHORT_STATISTIC_LEAK_SMALLEST_MEASURED,
                                                      NEW_ROW_BLOCK_INVARIANCE_ATOL,
                                                      NEW_ROW_BLOCK_INVARIANCE_MEASURED_ATOL,
                                                      NEW_ROW_BLOCK_INVARIANCE_RTOL,
                                                      ConfoundAdjustmentOperator, DesignSpec,
                                                      SitePooling, UnseenLevelError,
                                                      block_invariance_deviation)
from morpheus.v2.calibra.residualise import (apply_pooled_tissue_source_site, confound_design,
                                             cross_fitted_residuals, pooled_tissue_source_site)

SITES = ["01", "02", "05", "06", "09", "13", "17", "22", "31", "44"]
CANCERS = ["BRCA", "LUAD", "COAD", "GBM"]


def _cohort(n=240, p=12, seed=0, sites=SITES, cancers=CANCERS, rare_site=True):
    """Barcodes, a confound frame and a representation with real confound structure."""
    rng = np.random.default_rng(seed)
    site = rng.choice(sites, size=n)
    if rare_site:                      # guarantee at least one site below the pooling bar
        site[:3] = "99"
    cancer = rng.choice(cancers, size=n)
    patient_ids = np.asarray([f"TCGA-{s}-{i:04d}" for i, s in enumerate(site)])
    site_effect = {s: rng.normal(size=p) for s in np.unique(site)}
    cancer_effect = {c: rng.normal(size=p) for c in np.unique(cancer)}
    matrix = (rng.normal(size=(n, p))
              + np.stack([site_effect[s] for s in site])
              + np.stack([cancer_effect[c] for c in cancer]))
    purity = rng.uniform(0.2, 0.9, size=n)
    purity[rng.random(n) < 0.1] = np.nan          # exercise the missingness indicator
    frame = pd.DataFrame({"cancer": cancer, "tss": "", "purity": purity})
    return patient_ids, frame, matrix


def _fit(cohort_seed=0, **kwargs):
    """``cohort_seed`` draws the data; ``kwargs`` (n_splits/alpha/seed) tune the
    residualiser. Keeping the two names distinct matters -- collapsing them silently
    fits the operator at its default seed while the expectation uses another."""
    patient_ids, frame, matrix = _cohort(seed=cohort_seed)
    operator = ConfoundAdjustmentOperator.fit(
        matrix, frame, ["cancer", "tss", "purity"], patient_ids=patient_ids,
        site_column="tss", min_site_count=10, cohort_name="reference", **kwargs)
    return operator, patient_ids, frame, matrix


# --- 1. identity on the fitting cohort ------------------------------------------------

def test_operator_is_numerically_identical_to_cross_fitted_residuals():
    """THE test. Bit-for-bit, not merely close: every published number depends on it."""
    operator, patient_ids, frame, matrix = _fit()

    reference_frame = frame.copy()
    reference_frame["tss"], _ = pooled_tissue_source_site(patient_ids, min_site_count=10)
    design = confound_design(reference_frame, ["cancer", "tss", "purity"])
    expected = cross_fitted_residuals(design=design, matrix=matrix, n_splits=5, alpha=1.0, seed=42)

    actual = operator.adjust_reference(matrix, frame, patient_ids=patient_ids)
    assert actual.shape == expected.shape
    assert np.array_equal(actual, expected), "inductive reference path is not bit-identical"


def test_operator_design_matches_confound_design_on_the_reference():
    operator, patient_ids, frame, matrix = _fit()
    reference_frame = frame.copy()
    reference_frame["tss"], _ = pooled_tissue_source_site(patient_ids, min_site_count=10)
    expected = confound_design(reference_frame, ["cancer", "tss", "purity"])
    actual, _ = operator._frame_and_design(frame, patient_ids, "refuse")
    assert np.array_equal(actual, expected)
    assert len(operator.design_spec.design_columns) == expected.shape[1]


def test_identity_holds_for_every_residualiser_setting():
    for n_splits, alpha, seed in ((3, 1.0, 7), (5, 0.5, 42), (10, 2.0, 1)):
        operator, patient_ids, frame, matrix = _fit(n_splits=n_splits, alpha=alpha, seed=seed)
        reference_frame = frame.copy()
        reference_frame["tss"], _ = pooled_tissue_source_site(patient_ids, min_site_count=10)
        design = confound_design(reference_frame, ["cancer", "tss", "purity"])
        expected = cross_fitted_residuals(matrix, design, n_splits=n_splits, alpha=alpha, seed=seed)
        actual = operator.adjust_reference(matrix, frame, patient_ids=patient_ids)
        assert np.array_equal(actual, expected), (n_splits, alpha, seed)


def test_no_design_columns_degenerates_to_centring_exactly():
    rng = np.random.default_rng(3)
    matrix = rng.normal(size=(50, 6))
    frame = pd.DataFrame(index=range(50))
    operator = ConfoundAdjustmentOperator.fit(matrix, frame, [])
    expected = cross_fitted_residuals(matrix, np.zeros((50, 0)))
    assert np.array_equal(operator.adjust_reference(matrix, frame), expected)


# --- 2. an unseen site ----------------------------------------------------------------

def test_unseen_site_is_pooled_to_other_and_reported():
    """The predeclared policy: an unseen site is a site with reference count zero, and
    the reference's own rule already sends every below-threshold site to OTHER."""
    operator, patient_ids, frame, matrix = _fit()
    frequent = set(operator.site_pooling.frequent)
    assert "99" not in frequent                       # rare in the reference, pooled there too

    rng = np.random.default_rng(11)
    new_ids = np.asarray(["TCGA-ZZ-9001", "TCGA-YY-9002", "TCGA-01-9003"])
    new_frame = pd.DataFrame({"cancer": ["BRCA", "LUAD", "COAD"], "tss": "",
                              "purity": [0.5, 0.6, np.nan]})
    new_matrix = rng.normal(size=(3, matrix.shape[1]))

    adjusted, report = operator.adjust(new_matrix, new_frame, patient_ids=new_ids)
    assert adjusted.shape == (3, matrix.shape[1])
    assert np.isfinite(adjusted).all()
    pooling = report["site_pooling"]
    assert set(pooling["sites_pooled_to_other"]) == {"ZZ", "YY"}
    assert pooling["n_pooled_to_other"] == 2         # ZZ and YY; 01 is frequent in the reference
    assert "count of zero" in pooling["policy"]


def test_unseen_site_policy_equals_the_existing_pooling_rule_at_count_zero():
    """The policy is not a new invention: applying the frozen rule to reference
    barcodes reproduces ``pooled_tissue_source_site`` exactly."""
    patient_ids, _, _ = _cohort()
    expected, frequent = pooled_tissue_source_site(patient_ids, min_site_count=10)
    pooling = SitePooling.fit(patient_ids, min_site_count=10)
    assert set(pooling.frequent) == frequent
    actual, _ = pooling.apply(patient_ids)
    assert np.array_equal(actual, expected)
    assert np.array_equal(apply_pooled_tissue_source_site(patient_ids, frequent), expected)


def test_unseen_cancer_level_is_refused_not_silently_zeroed():
    """A zero one-hot row is the reference's baseline level, not "no confound"."""
    operator, patient_ids, frame, matrix = _fit()
    new_ids = np.asarray(["TCGA-01-9001"])
    new_frame = pd.DataFrame({"cancer": ["PAAD"], "tss": "", "purity": [0.5]})
    new_matrix = np.zeros((1, matrix.shape[1]))

    with pytest.raises(UnseenLevelError, match="PAAD"):
        operator.adjust(new_matrix, new_frame, patient_ids=new_ids)

    adjusted, report = operator.adjust(new_matrix, new_frame, patient_ids=new_ids,
                                       on_unseen_level="zero")
    assert np.isfinite(adjusted).all()
    assert report["unseen_levels"]["cancer"] == ["PAAD"]


# --- 3. a single new row --------------------------------------------------------------

def test_single_new_row_is_adjustable_at_all():
    """The whole point: one slide, no cohort, adjusted coordinates."""
    operator, patient_ids, frame, matrix = _fit()
    new_ids = np.asarray(["TCGA-05-9999"])
    new_frame = pd.DataFrame({"cancer": ["BRCA"], "tss": "", "purity": [0.42]})
    new_matrix = np.random.default_rng(5).normal(size=(1, matrix.shape[1]))

    adjusted, report = operator.adjust(new_matrix, new_frame, patient_ids=new_ids)
    assert adjusted.shape == (1, matrix.shape[1])
    assert np.isfinite(adjusted).all()
    assert report["n_rows_adjusted"] == 1
    assert report["path"] == "inductive_fold_ensemble"


def _new_rows(matrix, n=12, seed=23, with_purity=True):
    rng = np.random.default_rng(seed)
    ids = np.asarray([f"TCGA-0{i % 2 + 1}-8{i:03d}" for i in range(n)])
    columns = {"cancer": rng.choice(CANCERS, n), "tss": ""}
    if with_purity:
        columns["purity"] = rng.uniform(0.2, 0.9, n)
    return ids, pd.DataFrame(columns), rng.normal(size=(n, matrix.shape[1]))


def test_one_row_at_a_time_equals_the_whole_block():
    """No cohort statistic leaks into the new-row path: adjusting rows one by one must
    equal adjusting them together. This is the property the transductive function
    lacks, and the reason it could not be deployed.

    Asserted in two pieces, because they have different strengths and conflating
    them is what made the original single ``atol=0, rtol=0`` assertion untenable.

    *Everything this operator controls is bit-exact* -- the design encoding is
    row-wise, and a row's answer does not move when the OTHER rows of a block of the
    same size are replaced with different data. Those are the assertions a leaked
    cohort statistic breaks, and they use ``array_equal``.

    *What BLAS controls is not.* ``numpy`` sends ``(1, d) @ (d, k)`` to a GEMV and
    ``(m, d) @ (d, k)`` to a blocked GEMM, which sum the ``d`` products in different
    orders, so changing the block's SIZE moves the last bit. This is not a property
    of this module: ``A[i:i+1] @ B != (A @ B)[i:i+1]`` for dense operands under
    numpy 2.2.6 and 2.4.3 alike (both scipy-openblas). It is bounded here rather
    than waved through -- see
    ``test_the_block_size_tolerance_cannot_absorb_a_cohort_statistic``.
    """
    operator, patient_ids, frame, matrix = _fit()
    new_ids, new_frame, new_matrix = _new_rows(matrix)
    deviation = block_invariance_deviation(operator, new_matrix, new_frame,
                                           patient_ids=new_ids)

    # (a) exact, in everything the operator owns.
    assert deviation["design_bit_exact"], "the design encoding is not row-wise"
    assert deviation["composition_bit_exact"], \
        "a row's answer moved when the rest of the block changed: a cohort statistic leaked"

    # (b) at the last bit, in what BLAS owns. Both a scale-free bound (the deviation
    #     must stay within a few ulps of the values themselves) and the flat tolerance
    #     the published contract quotes.
    assert deviation["max_abs"] <= 64 * deviation["eps"] * max(1.0, deviation["max_abs_value"])
    assert deviation["max_abs"] <= NEW_ROW_BLOCK_INVARIANCE_MEASURED_ATOL

    block, _ = operator.adjust(new_matrix, new_frame, patient_ids=new_ids)
    for i in range(12):
        single, _ = operator.adjust(new_matrix[i:i + 1], new_frame.iloc[i:i + 1],
                                    patient_ids=new_ids[i:i + 1])
        assert np.allclose(single, block[i:i + 1], atol=NEW_ROW_BLOCK_INVARIANCE_ATOL,
                           rtol=NEW_ROW_BLOCK_INVARIANCE_RTOL)


def test_a_row_is_unmoved_by_whatever_else_is_in_its_block():
    """The bit-exact half, stated on its own because it is the real guarantee.

    Same block size, row 0 held fixed, every other row replaced by a different
    cancer, a different site and values 100x larger. Nothing the operator does may
    notice: the levels, the numeric mean/s.d. and the centre are all frozen against
    the reference. A permutation of the block must likewise return every row
    unchanged.
    """
    operator, patient_ids, frame, matrix = _fit()
    new_ids, new_frame, new_matrix = _new_rows(matrix)
    block, _ = operator.adjust(new_matrix, new_frame, patient_ids=new_ids)

    rng = np.random.default_rng(999)
    other_ids, other_frame, other_matrix = new_ids.copy(), new_frame.copy(), new_matrix.copy()
    other_ids[1:] = np.asarray([f"TCGA-{'ZZ' if i % 2 else '05'}-7{i:03d}" for i in range(1, 12)])
    other_frame.loc[1:, "cancer"] = rng.choice(CANCERS, 11)
    other_frame.loc[1:, "purity"] = rng.uniform(0.2, 0.9, 11)
    other_matrix[1:] = rng.normal(size=(11, matrix.shape[1])) * 100.0
    alternative, _ = operator.adjust(other_matrix, other_frame, patient_ids=other_ids)
    assert np.array_equal(alternative[0:1], block[0:1])

    order = np.random.default_rng(4).permutation(12)
    permuted, _ = operator.adjust(new_matrix[order], new_frame.iloc[order].reset_index(drop=True),
                                  patient_ids=new_ids[order])
    assert np.array_equal(permuted[np.argsort(order)], block)


def test_a_categorical_only_design_is_bit_exact_block_versus_single():
    """The deployed case, and the reason no published number rests on the tolerance.

    All 14 persisted operators at ``runs_misc/tcga_operators/`` have purely
    categorical designs (``cancer`` and ``tss``), so every design row is one-hot:
    the products are exactly representable and the remaining terms are zeros, which
    re-associate exactly. Measured over those 14 operators the deviation is
    ``0.000e+00``; reproduced here without the box.
    """
    patient_ids, frame, matrix = _cohort()
    operator = ConfoundAdjustmentOperator.fit(
        matrix, frame[["cancer", "tss"]], ["cancer", "tss"], patient_ids=patient_ids,
        site_column="tss", min_site_count=10, cohort_name="categorical_only")
    new_ids, new_frame, new_matrix = _new_rows(matrix, with_purity=False)
    deviation = block_invariance_deviation(operator, new_matrix, new_frame[["cancer", "tss"]],
                                           patient_ids=new_ids)
    assert deviation["design_bit_exact"] and deviation["composition_bit_exact"]
    assert deviation["max_abs"] == 0.0, "a categorical-only design must be bit-exact"
    assert deviation["n_differing"] == 0


def test_the_block_size_tolerance_cannot_absorb_a_cohort_statistic():
    """A tolerance is only defensible against the thing it must still catch.

    Two real leaks, constructed on the same block the tolerance is measured on:
    re-centring on the new cohort (what the transductive function does) and
    standardising a numeric covariate by the new cohort's own mean and s.d. The
    smaller of the two is ~1.3e11 times the tolerance, so nothing that would change
    a published number can hide underneath it.
    """
    operator, patient_ids, frame, matrix = _fit()
    new_ids, new_frame, new_matrix = _new_rows(matrix)
    block, _ = operator.adjust(new_matrix, new_frame, patient_ids=new_ids)

    recentred = float(np.max(np.abs(block - (block - block.mean(axis=0)))))

    reference_mean, reference_std, _ = operator.design_spec.numeric_stats["purity"]
    purity = new_frame["purity"].to_numpy(dtype=float)
    on_reference_scale = (purity - reference_mean) / reference_std
    on_new_cohort_scale = (purity - purity.mean()) / purity.std()
    column = operator.design_spec.design_columns.index("purity")
    ensemble_coefficient = np.mean([c[:, column] for c in operator.residualiser.coefficients],
                                   axis=0)
    restandardised = float(np.max(np.abs(np.outer(on_new_cohort_scale - on_reference_scale,
                                                  ensemble_coefficient))))

    smallest_leak = min(recentred, restandardised)
    assert smallest_leak == pytest.approx(COHORT_STATISTIC_LEAK_SMALLEST_MEASURED, rel=1e-3)
    assert smallest_leak > 1e10 * NEW_ROW_BLOCK_INVARIANCE_ATOL
    assert NEW_ROW_BLOCK_INVARIANCE_ATOL > 1e3 * NEW_ROW_BLOCK_INVARIANCE_MEASURED_ATOL


def test_transductive_residualiser_does_NOT_have_that_property():
    """The motivating defect, asserted rather than asserted-about: run the transductive
    function on a block and on one row of it and the answers differ."""
    patient_ids, frame, matrix = _cohort()
    frame = frame.copy()
    frame["tss"], _ = pooled_tissue_source_site(patient_ids, min_site_count=10)
    design = confound_design(frame, ["cancer", "tss", "purity"])
    block = cross_fitted_residuals(matrix, design)
    subset = cross_fitted_residuals(matrix[:40], design[:40])
    assert not np.allclose(subset, block[:40])


def test_reference_and_new_row_paths_differ_and_that_is_documented():
    """The deliberate discontinuity. Reference rows keep their own held-out fold model
    (which is what makes identity exact); a new row gets the K-model ensemble."""
    operator, patient_ids, frame, matrix = _fit()
    reference = operator.adjust_reference(matrix, frame, patient_ids=patient_ids)
    through_new_path, _ = operator.adjust(matrix, frame, patient_ids=patient_ids)
    assert not np.array_equal(reference, through_new_path)
    # ...but the two must stay in the same coordinate system, not drift apart.
    assert np.corrcoef(reference.ravel(), through_new_path.ravel())[0, 1] > 0.95


# --- 4. round trip through persistence ------------------------------------------------

def test_round_trip_through_persistence_is_exact(tmp_path):
    operator, patient_ids, frame, matrix = _fit()
    path = operator.save(tmp_path / "operator.npz")
    restored = ConfoundAdjustmentOperator.load(path)

    assert np.array_equal(restored.adjust_reference(matrix, frame, patient_ids=patient_ids),
                          operator.adjust_reference(matrix, frame, patient_ids=patient_ids))

    rng = np.random.default_rng(31)
    new_ids = np.asarray(["TCGA-02-7001", "TCGA-QQ-7002"])
    new_frame = pd.DataFrame({"cancer": ["GBM", "BRCA"], "tss": "", "purity": [0.3, np.nan]})
    new_matrix = rng.normal(size=(2, matrix.shape[1]))
    before, before_report = operator.adjust(new_matrix, new_frame, patient_ids=new_ids)
    after, after_report = restored.adjust(new_matrix, new_frame, patient_ids=new_ids)
    assert np.array_equal(before, after)
    assert before_report["site_pooling"]["sites_pooled_to_other"] == \
        after_report["site_pooling"]["sites_pooled_to_other"]


def test_persisted_operator_carries_its_provenance(tmp_path):
    operator, patient_ids, frame, matrix = _fit()
    restored = ConfoundAdjustmentOperator.load(operator.save(tmp_path / "operator.npz"))
    provenance = restored.provenance
    assert provenance["reference_digest"] == operator.provenance["reference_digest"]
    assert len(provenance["reference_digest"]) == 64
    assert provenance["design_columns"] == operator.provenance["design_columns"]
    assert provenance["site_pooling"]["min_site_count"] == 10
    assert provenance["site_pooling"]["n_frequent_sites"] == len(operator.site_pooling.frequent)
    assert provenance["residualiser"] == {"n_splits": 5, "alpha": 1.0, "seed": 42}
    assert provenance["n_reference_rows"] == len(matrix)
    for key in ("numpy_version", "sklearn_version", "created_utc", "cohort_name"):
        assert key in provenance


def test_reference_digest_separates_cohorts():
    """An applied adjustment must be traceable to the cohort that fitted it."""
    first, _, _, _ = _fit(cohort_seed=0)
    second, patient_ids, frame, matrix = _fit(cohort_seed=1)
    assert first.provenance["reference_digest"] != second.provenance["reference_digest"]
    # Same cohort, same digest -- so the digest identifies the cohort, not the run.
    again, _, _, _ = _fit(cohort_seed=0)
    assert first.provenance["reference_digest"] == again.provenance["reference_digest"]


def test_applied_report_names_the_fitting_cohort():
    operator, patient_ids, frame, matrix = _fit()
    _, report = operator.adjust(np.zeros((1, matrix.shape[1])),
                                pd.DataFrame({"cancer": ["BRCA"], "tss": "", "purity": [0.5]}),
                                patient_ids=np.asarray(["TCGA-06-9000"]))
    assert report["applied_by"]["reference_digest"] == operator.provenance["reference_digest"]
    assert report["applied_by"]["cohort_name"] == "reference"


# --- cross-cohort use, which is the point ---------------------------------------------

def test_a_second_cohort_can_be_adjusted_with_the_first_cohorts_operator():
    """ALCHEMIST adjusted with TCGA's operator: a disjoint cohort, disjoint sites,
    adjusted into the reference's coordinate system without refitting anything."""
    operator, patient_ids, frame, matrix = _fit()
    other_ids, other_frame, other_matrix = _cohort(n=90, seed=99,
                                                   sites=["AA", "AB", "AC"], rare_site=False)
    adjusted, report = operator.adjust(other_matrix, other_frame, patient_ids=other_ids)
    assert adjusted.shape == other_matrix.shape
    assert np.isfinite(adjusted).all()
    assert report["site_pooling"]["n_pooled_to_other"] == 90        # every site is unseen
    assert set(report["site_pooling"]["sites_pooled_to_other"]) == {"AA", "AB", "AC"}
    # The operator, not the new cohort, supplied the numeric standardisation -- otherwise
    # the two cohorts would be adjusted in two different coordinate systems, which is
    # exactly the comparison error this operator exists to prevent.
    reference_mean = float(np.nanmean(frame["purity"].to_numpy(dtype=float)))
    other_mean = float(np.nanmean(other_frame["purity"].to_numpy(dtype=float)))
    assert operator.design_spec.numeric_stats["purity"][0] == pytest.approx(reference_mean)
    assert operator.design_spec.numeric_stats["purity"][0] != pytest.approx(other_mean)


def test_design_spec_refuses_an_unknown_unseen_level_policy():
    spec = DesignSpec.fit(pd.DataFrame({"cancer": ["BRCA", "LUAD"]}), ["cancer"])
    with pytest.raises(ValueError, match="refuse"):
        spec.transform(pd.DataFrame({"cancer": ["BRCA"]}), on_unseen_level="nearest")
