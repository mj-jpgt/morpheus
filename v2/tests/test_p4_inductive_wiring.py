"""Tests for P4's discovery/exposure split and its use of the inductive operator.

The point of the wiring is that the state P4 certifies must be producible for a patient
who was not in the cohort the nuisance model was fitted on. Two tests carry that:

* ``test_an_exposure_row_is_unmoved_by_the_other_exposure_rows`` -- the property the
  whole exercise exists for, asserted with ``array_equal``; and
* ``test_the_matched_transductive_control_does_NOT_have_that_property`` -- the same
  assertion on the control, which must fail, so the first test cannot be passing for a
  trivial reason.

``test_transductive_mode_still_equals_cross_fitted_residuals`` guards the published
2026-08-04 20:00 numbers: the default path must not have moved.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.calibra.confound_certificate import certify_axes
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.research.rebase.nature.p4_certification.p4_certify import (ADJUSTMENTS,
                                                                            certificate_arms,
                                                                            design_for,
                                                                            exposure_split,
                                                                            prepare_state)

SITES = ["01", "02", "05", "06", "09", "13", "17", "22", "31", "44", "55", "61"]
CANCERS = ["BRCA", "LUAD", "COAD", "GBM", "STAD"]


def _cohort(n: int = 300, p: int = 10, seed: int = 0):
    """Barcodes, cancers and a representation carrying real site and cancer structure."""
    rng = np.random.default_rng(seed)
    site = rng.choice(SITES, size=n)
    site[:4] = "99"                       # a site below any sensible pooling bar
    cancers = rng.choice(CANCERS, size=n)
    patient_ids = np.asarray([f"TCGA-{s}-{i:04d}" for i, s in enumerate(site)])
    site_effect = {s: rng.normal(size=p) for s in np.unique(site)}
    cancer_effect = {c: rng.normal(size=p) for c in np.unique(cancers)}
    features = (rng.normal(size=(n, p))
                + np.stack([site_effect[s] for s in site])
                + np.stack([cancer_effect[c] for c in cancers]))
    targets = rng.normal(size=(n, 4)) + features[:, :4] * 0.5
    return features, targets, patient_ids, cancers


def _prepare(features, targets, patient_ids, cancers, adjustment, *, fraction=0.5,
             min_site_count=5, seed=42):
    return prepare_state(features, targets, patient_ids, cancers, adjustment=adjustment,
                         discovery_fraction=fraction, seed=seed, min_site_count=min_site_count)


# ------------------------------------------------------------------ the split


def test_the_split_shares_no_patient_and_keeps_every_cancer_on_both_sides():
    _, _, _, cancers = _cohort()
    discovery = exposure_split(cancers, discovery_fraction=0.5, seed=42)
    assert discovery.any() and (~discovery).any()
    assert not (set(np.flatnonzero(discovery)) & set(np.flatnonzero(~discovery)))
    # Every cancer on both sides: a cancer missing from the discovery fold is an unseen
    # one-hot level, which DesignSpec.transform refuses outright.
    assert set(cancers[discovery]) == set(cancers) == set(cancers[~discovery])


def test_the_split_is_deterministic_and_moves_with_the_seed():
    _, _, _, cancers = _cohort()
    a = exposure_split(cancers, discovery_fraction=0.5, seed=42)
    assert np.array_equal(a, exposure_split(cancers, discovery_fraction=0.5, seed=42))
    assert not np.array_equal(a, exposure_split(cancers, discovery_fraction=0.5, seed=43))


@pytest.mark.parametrize("fraction", [0.0, 1.0, -0.2, 1.5])
def test_a_degenerate_discovery_fraction_is_refused(fraction):
    _, _, _, cancers = _cohort()
    with pytest.raises(ValueError):
        exposure_split(cancers, discovery_fraction=fraction, seed=42)


# ------------------------------------------------------- the published path is unmoved


def test_transductive_mode_still_equals_cross_fitted_residuals():
    """The default path must reproduce the 2026-08-04 20:00 run bit for bit."""
    features, targets, patient_ids, cancers = _cohort()
    state = _prepare(features, targets, patient_ids, cancers, "transductive")
    design = design_for(patient_ids, cancers, 5)
    assert np.array_equal(state["adjusted_features"],
                          cross_fitted_residuals(features, design, seed=42))
    assert np.array_equal(state["adjusted_targets"],
                          cross_fitted_residuals(targets, design, seed=42))
    assert state["rows"].tolist() == list(range(len(features)))
    assert state["exposable"] is False


def test_only_the_inductive_arm_calls_its_state_exposable():
    features, targets, patient_ids, cancers = _cohort()
    exposable = {mode: _prepare(features, targets, patient_ids, cancers, mode)["exposable"]
                 for mode in ADJUSTMENTS}
    assert exposable == {"transductive": False, "transductive_exposure": False,
                         "inductive": True}


# --------------------------------------------------- the property that makes it exposable


def _perturb_other_exposure_rows(features, patient_ids, cancers, *, fraction=0.5, seed=42):
    """Replace every exposure row but the first with different data, leaving the split fixed.

    The split is derived from ``cancers`` alone, so this cannot move which rows are
    scored -- only what is in them. A row whose adjusted coordinates depend on its
    companions is a row no deployed system can produce.
    """
    discovery = exposure_split(cancers, discovery_fraction=fraction, seed=seed)
    rows = np.flatnonzero(~discovery)
    perturbed = np.array(features, dtype=np.float64, copy=True)
    perturbed[rows[1:]] = perturbed[rows[1:]] * 100.0 + 7.0
    return perturbed, rows


def test_an_exposure_row_is_unmoved_by_the_other_exposure_rows():
    features, targets, patient_ids, cancers = _cohort()
    baseline = _prepare(features, targets, patient_ids, cancers, "inductive")
    perturbed_features, rows = _perturb_other_exposure_rows(features, patient_ids, cancers)
    perturbed = _prepare(perturbed_features, targets, patient_ids, cancers, "inductive")
    assert np.array_equal(baseline["rows"], perturbed["rows"])
    assert np.array_equal(baseline["adjusted_features"][0:1],
                          perturbed["adjusted_features"][0:1])


def test_the_matched_transductive_control_does_NOT_have_that_property():
    """The motivating defect, asserted rather than described."""
    features, targets, patient_ids, cancers = _cohort()
    baseline = _prepare(features, targets, patient_ids, cancers, "transductive_exposure")
    perturbed_features, _ = _perturb_other_exposure_rows(features, patient_ids, cancers)
    perturbed = _prepare(perturbed_features, targets, patient_ids, cancers,
                         "transductive_exposure")
    assert not np.allclose(baseline["adjusted_features"][0:1],
                           perturbed["adjusted_features"][0:1])


def test_the_inductive_operator_never_saw_the_rows_it_scores():
    features, targets, patient_ids, cancers = _cohort()
    state = _prepare(features, targets, patient_ids, cancers, "inductive")
    discovery = exposure_split(cancers, discovery_fraction=0.5, seed=42)
    fit_ids = set(patient_ids[discovery].tolist())
    scored_ids = set(patient_ids[state["rows"]].tolist())
    assert not (fit_ids & scored_ids)
    assert state["report"]["n_discovery"] == int(discovery.sum())
    assert state["report"]["n_scored"] == len(scored_ids)


def test_the_inductive_report_names_the_site_coverage_it_actually_achieved():
    """§4 of the predeclaration: the operator can only site-adjust rows from sites it saw."""
    features, targets, patient_ids, cancers = _cohort()
    report = _prepare(features, targets, patient_ids, cancers, "inductive")["report"]
    covered = report["n_exposure_rows_with_a_site_adjustment"]
    assert covered + report["n_exposure_rows_pooled_to_OTHER"] == report["n_scored"]
    assert 0.0 <= report["fraction_exposure_rows_with_a_site_adjustment"] <= 1.0
    assert report["n_frequent_sites_in_discovery_fold"] >= 1
    assert report["operator_reference_digest"]


def test_the_operators_design_is_the_published_confound_design_on_the_discovery_fold():
    features, targets, patient_ids, cancers = _cohort()
    state = _prepare(features, targets, patient_ids, cancers, "inductive")
    discovery = exposure_split(cancers, discovery_fraction=0.5, seed=42)
    site, _ = pooled_tissue_source_site(patient_ids[discovery], min_site_count=5)
    expected = confound_design(pd.DataFrame({"cancer": cancers[discovery], "tss": site}),
                               ["cancer", "tss"])
    assert state["report"]["n_operator_design_columns"] == expected.shape[1]


# ------------------------------------------------------------------ the audit and the certificate


def test_the_adjustment_audit_would_expose_a_null_adjustment():
    """§6.1: a certificate PASS on a state the adjustment did not move is worth nothing."""
    features, targets, patient_ids, cancers = _cohort()
    real = _prepare(features, targets, patient_ids, cancers, "inductive")["audit"]
    # A representation with no confound structure at all: the adjustment can only be a
    # near-identity, and the audit must say so.
    rng = np.random.default_rng(1)
    null_features = rng.normal(size=features.shape)
    null_audit = _prepare(null_features, targets, patient_ids, cancers, "inductive")["audit"]
    assert null_audit["per_axis_corr_raw_adjusted_median"] > \
        real["per_axis_corr_raw_adjusted_median"]
    assert real["residual_variance_ratio_median"] < 1.0


def test_certificate_arms_agree_with_the_internal_residualise_branch():
    """The transductive arms lean on ``certify_axes(residualise=True)``; the inductive arm
    hands over a pre-adjusted block. This pins that the two routes are the same route."""
    features, _, patient_ids, cancers = _cohort(n=160, p=6)

    class _Args:
        min_site_count, n_permutations, seed = 5, 12, 42
        n_boot, n_boot_axes, n_jobs = 20, 2, 1

    state = _prepare(features, None, patient_ids, cancers, "transductive", min_site_count=5)
    internal = certificate_arms(state, _Args())["adjusted"]
    external = certify_axes(state["adjusted_features"], patient_ids, cancers,
                            min_site_count=5, n_permutations=12, seed=42, n_boot=20,
                            n_boot_axes=2, residualise=False, n_jobs=1)
    assert internal["joint_lda_balanced_accuracy"] == \
        pytest.approx(external["joint_lda_balanced_accuracy"], abs=1e-12)
    assert np.allclose(internal["per_axis_balanced_accuracy"],
                       external["per_axis_balanced_accuracy"], atol=1e-12)
