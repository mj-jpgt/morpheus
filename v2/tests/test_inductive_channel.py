"""Tests for the inductive channel / ceiling driver.

Three things carry this run, and each is a test rather than an argument:

* ``test_the_inductive_adjuster_reproduces_p4s_prepare_state_bit_for_bit`` -- the driver's
  adjuster and P4's ``prepare_state`` must produce the *same* adjusted block, or the
  channel measured here is not on the state P4 measured joint site LDA 0.2643 on;
* ``test_a_scored_row_is_unmoved_by_the_other_scored_rows`` and its negative twin
  ``test_the_transductive_adjuster_does_NOT_have_that_property`` -- the inductive adjuster
  must be a fixed map, so that permuting ``y`` underneath it in the pairing null is a
  well-defined operation; the twin exists so the first cannot pass trivially;
* ``test_adjust_y_defaults_to_adjust_byte_for_byte`` -- the new ``adjust_y`` argument on
  ``channel_under_adjustment`` must leave every existing caller unmoved.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.calibra.nonlinear_adjustment import (channel_under_adjustment,
                                                      labels_only_ceiling, make_adjuster)
from morpheus.v2.calibra.residualise import confound_design, pooled_tissue_source_site
from morpheus.v2.research.rebase.nature.p1_evidence.inductive_channel import (
    fit_inductive_adjuster)
from morpheus.v2.research.rebase.nature.p4_certification.p4_certify import (exposure_split,
                                                                            prepare_state)

SITES = ["01", "02", "05", "06", "09", "13", "17", "22", "31", "44"]
CANCERS = ["BRCA", "LUAD", "COAD", "GBM"]


def _cohort(n: int = 240, p: int = 8, k: int = 5, seed: int = 0):
    rng = np.random.default_rng(seed)
    site = rng.choice(SITES, size=n)
    cancers = rng.choice(CANCERS, size=n)
    patient_ids = np.asarray([f"TCGA-{s}-{i:04d}" for i, s in enumerate(site)])
    site_effect = {s: rng.normal(size=p) for s in np.unique(site)}
    cancer_effect = {c: rng.normal(size=p) for c in np.unique(cancers)}
    x = (rng.normal(size=(n, p)) + np.stack([site_effect[s] for s in site])
         + np.stack([cancer_effect[c] for c in cancers]))
    y = rng.normal(size=(n, k)) + x[:, :k] * 0.6
    return x, y, patient_ids, cancers


def _folds(cancers, fraction=0.5, seed=42):
    discovery = exposure_split(cancers, discovery_fraction=fraction, seed=seed)
    return np.flatnonzero(discovery), np.flatnonzero(~discovery)


# --- the contract with P4's state ---------------------------------------------------------

def test_the_inductive_adjuster_reproduces_p4s_prepare_state_bit_for_bit():
    x, y, ids, cancers = _cohort()
    fit_rows, score_rows = _folds(cancers)
    adjust_x, _ = fit_inductive_adjuster(x, fit_rows, score_rows, cancers, ids, min_site_count=5)
    adjust_y, _ = fit_inductive_adjuster(y, fit_rows, score_rows, cancers, ids, min_site_count=5)
    state = prepare_state(x, y, ids, cancers, adjustment="inductive", discovery_fraction=0.5,
                          seed=42, min_site_count=5)
    assert np.array_equal(adjust_x(x[score_rows]), state["adjusted_features"])
    assert np.array_equal(adjust_y(y[score_rows]), state["adjusted_targets"])
    assert np.array_equal(ids[score_rows], state["patient_ids"])


def test_the_operator_never_saw_the_rows_it_scores():
    _, _, ids, cancers = _cohort()
    fit_rows, score_rows = _folds(cancers)
    assert not (set(ids[fit_rows].tolist()) & set(ids[score_rows].tolist()))


# --- the fixed-map property, and its negative twin ----------------------------------------

def test_a_scored_row_is_unmoved_by_the_other_scored_rows():
    """The property that makes permuting ``y`` under the adjuster well-defined."""
    x, _, ids, cancers = _cohort()
    fit_rows, score_rows = _folds(cancers)
    adjust, _ = fit_inductive_adjuster(x, fit_rows, score_rows, cancers, ids, min_site_count=5)
    block = x[score_rows]
    reference = adjust(block)
    disturbed = block.copy()
    rng = np.random.default_rng(7)
    disturbed[1:] = rng.normal(size=disturbed[1:].shape) * 100.0
    assert np.array_equal(adjust(disturbed)[0], reference[0])


def test_the_transductive_adjuster_does_NOT_have_that_property():
    """So the test above cannot be passing for a trivial reason."""
    x, _, ids, cancers = _cohort()
    _, score_rows = _folds(cancers)
    site, _ = pooled_tissue_source_site(ids[score_rows], min_site_count=5)
    design = confound_design(pd.DataFrame({"cancer": cancers[score_rows], "tss": site}),
                             ["cancer", "tss"])
    adjust = make_adjuster("ridge", design=design, n_splits=5, seed=42)
    block = x[score_rows]
    reference = adjust(block)
    disturbed = block.copy()
    rng = np.random.default_rng(7)
    disturbed[1:] = rng.normal(size=disturbed[1:].shape) * 100.0
    assert not np.array_equal(adjust(disturbed)[0], reference[0])


def test_the_adjuster_refuses_a_block_of_the_wrong_length():
    x, _, ids, cancers = _cohort()
    fit_rows, score_rows = _folds(cancers)
    adjust, _ = fit_inductive_adjuster(x, fit_rows, score_rows, cancers, ids, min_site_count=5)
    with pytest.raises(ValueError):
        adjust(x[score_rows][:-1])


def test_the_operator_carries_a_reference_digest_of_its_fitting_rows_only():
    x, _, ids, cancers = _cohort()
    fit_rows, score_rows = _folds(cancers)
    _, operator = fit_inductive_adjuster(x, fit_rows, score_rows, cancers, ids, min_site_count=5)
    assert operator.provenance["n_reference_rows"] == len(fit_rows)
    perturbed = x.copy()
    perturbed[score_rows] *= 3.0            # only the SCORED rows move
    _, other = fit_inductive_adjuster(perturbed, fit_rows, score_rows, cancers, ids,
                                      min_site_count=5)
    assert operator.provenance["reference_digest"] == other.provenance["reference_digest"]


def test_the_frozen_design_spec_encodes_every_row_in_the_fitting_rows_columns():
    """The third labels encoding: the operator's OWN design, evaluated on the whole partition.

    It has to be defined on the discovery rows (which the inductive operator is fitted on)
    and on the exposure rows (which it scores), in the same columns, or a labels block built
    from it is not one block. Row-wise encoding with frozen levels and frozen site pooling is
    what makes that true, and it is asserted rather than assumed.
    """
    x, _, ids, cancers = _cohort()
    fit_rows, score_rows = _folds(cancers)
    _, operator = fit_inductive_adjuster(x, fit_rows, score_rows, cancers, ids, min_site_count=5)
    whole, _ = operator._frame_and_design(pd.DataFrame({"cancer": cancers}), ids, "refuse")
    on_fit, _ = operator._frame_and_design(pd.DataFrame({"cancer": cancers[fit_rows]}),
                                           ids[fit_rows], "refuse")
    assert whole.shape[1] == operator.provenance["n_design_columns"]
    assert np.array_equal(whole[fit_rows], on_fit)


# --- the new adjust_y argument ------------------------------------------------------------

def test_adjust_y_defaults_to_adjust_byte_for_byte():
    x, y, ids, cancers = _cohort()
    site, _ = pooled_tissue_source_site(ids, min_site_count=5)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": site}), ["cancer", "tss"])
    adjust = make_adjuster("ridge", design=design, n_splits=5, seed=42)
    default = channel_under_adjustment(x, y, adjust, strata=cancers, n_permutations=20,
                                       n_components=3, seed=42)
    explicit = channel_under_adjustment(x, y, adjust, strata=cancers, n_permutations=20,
                                        n_components=3, seed=42, adjust_y=adjust)
    assert default == explicit


def test_adjust_y_actually_routes_the_target_block_through_the_second_adjuster():
    x, y, ids, cancers = _cohort()
    fit_rows, score_rows = _folds(cancers)
    adjust_x, _ = fit_inductive_adjuster(x, fit_rows, score_rows, cancers, ids, min_site_count=5)
    adjust_y, _ = fit_inductive_adjuster(y, fit_rows, score_rows, cancers, ids, min_site_count=5)
    x_e, y_e = x[score_rows], y[score_rows]
    # The x operator has 8 output columns and y has 5, so a run that ignored ``adjust_y``
    # would raise rather than silently mis-adjust.
    record = channel_under_adjustment(x_e, y_e, adjust_x, strata=cancers[score_rows],
                                      n_permutations=10, n_components=3, seed=42,
                                      adjust_y=adjust_y)
    assert np.isfinite(record["observed_top_cca"])
    with pytest.raises(ValueError):
        channel_under_adjustment(x_e, y_e, adjust_x, strata=cancers[score_rows],
                                 n_permutations=0, n_components=3, seed=42)


def test_labels_only_ceiling_threads_adjust_y_through():
    x, y, ids, cancers = _cohort()
    fit_rows, score_rows = _folds(cancers)
    site, _ = pooled_tissue_source_site(ids, min_site_count=5)
    labels = confound_design(pd.DataFrame({"cancer": cancers, "tss": site}), ["cancer", "tss"])
    adjust_labels, _ = fit_inductive_adjuster(labels, fit_rows, score_rows, cancers, ids,
                                              min_site_count=5)
    adjust_y, _ = fit_inductive_adjuster(y, fit_rows, score_rows, cancers, ids, min_site_count=5)
    record = labels_only_ceiling(labels[score_rows], y[score_rows], n_components=3, seed=42,
                                 adjust=adjust_labels, strata=cancers[score_rows],
                                 n_permutations=10, adjust_y=adjust_y)
    assert np.isfinite(record["raw_top_cca"])
    assert np.isfinite(record["adjusted"]["observed_top_cca"])
    del x
