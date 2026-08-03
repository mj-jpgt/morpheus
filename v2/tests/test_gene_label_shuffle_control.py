"""Self-tests for the T1.5 gene-label-shuffle control.

A must-FAIL control is only evidence if it can be shown to react. Each test below
plants the situation the control exists to catch and asserts it reacts -- in both
directions, because a statistic that only ever says "collapsed" is not measuring
anything.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.build_shuffled_pbs_targets import reconstruct_scores
from morpheus.v2.calibra.gene_label_shuffle_control import (attribution_collapse, rank_correlation,
                                                            subspace_persistence)


def test_rank_correlation_matches_known_values():
    values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert rank_correlation(values, values) == pytest.approx(1.0)
    assert rank_correlation(values, -values) == pytest.approx(-1.0)
    # monotone but non-linear: Spearman is 1 where Pearson would not be
    assert rank_correlation(values, np.exp(values)) == pytest.approx(1.0)


def test_attribution_collapse_reports_no_collapse_when_the_shuffle_did_not_take_effect():
    """The direction that actually protects us: an ineffective shuffle must be visible."""
    rng = np.random.default_rng(0)
    basis = rng.normal(size=(400, 8))
    verdict = attribution_collapse(basis, basis, n_boot=200, seed=0)
    assert verdict["median_abs_spearman_by_axis_index"] == pytest.approx(1.0)
    assert verdict["attribution_collapsed"] is False
    assert verdict["median_best_match_abs_spearman"] == pytest.approx(1.0)


def test_attribution_collapse_detects_a_real_row_permutation():
    rng = np.random.default_rng(1)
    basis = rng.normal(size=(2000, 16))
    permuted = basis[rng.permutation(len(basis))]
    verdict = attribution_collapse(basis, permuted, n_boot=200, seed=1)
    assert verdict["attribution_collapsed"] is True
    assert verdict["median_abs_spearman_by_axis_index"] < 0.05
    # the strictly harder statistic must also stay low, otherwise "collapsed" only
    # means "axis k moved", not "the attribution is gone"
    assert verdict["median_best_match_abs_spearman"] < 0.10, verdict


def test_attribution_collapse_survives_a_pure_axis_reordering():
    """Permuting which COLUMN is axis k is not a gene-label shuffle. The by-index
    statistic must drop, and the best-match statistic must not -- that is the whole
    reason both are reported."""
    rng = np.random.default_rng(2)
    basis = rng.normal(size=(500, 12))
    reordered = basis[:, rng.permutation(basis.shape[1])]
    verdict = attribution_collapse(basis, reordered, n_boot=200, seed=2)
    assert verdict["median_abs_spearman_by_axis_index"] < 0.2
    assert verdict["median_best_match_abs_spearman"] == pytest.approx(1.0)


def test_attribution_collapse_refuses_a_shape_mismatch():
    assert attribution_collapse(np.zeros((10, 3)), np.zeros((10, 4)))["status"] == "unavailable_shape_mismatch"


def _paired_blocks(n=400, seed=0):
    rng = np.random.default_rng(seed)
    latent = rng.normal(size=(n, 4))
    x = latent @ rng.normal(size=(4, 12)) + rng.normal(scale=0.5, size=(n, 12))
    y = latent @ rng.normal(size=(4, 6)) + rng.normal(scale=0.5, size=(n, 6))
    return rng, x, y


def test_subspace_persists_under_a_rotation_that_preserves_the_information():
    """A rotated target block carries the same subspace, so the control must say so."""
    rng, x, y = _paired_blocks(seed=3)
    rotation = np.linalg.qr(rng.normal(size=(y.shape[1], y.shape[1])))[0]
    verdict = subspace_persistence(x, y, y @ rotation, np.zeros((len(x), 0)),
                                   n_components=4, n_boot=40, seed=3)
    assert verdict["subspace_persists"] is True
    assert verdict["difference_ci_excludes_zero"] is False


def test_subspace_does_not_persist_when_the_shuffled_block_is_noise():
    """The other direction. If the control cannot say 'destroyed', it says nothing."""
    rng, x, y = _paired_blocks(seed=4)
    noise = rng.normal(size=y.shape)
    verdict = subspace_persistence(x, y, noise, np.zeros((len(x), 0)),
                                   n_components=4, n_boot=40, seed=4)
    assert verdict["heldout_top_cca_true"] > verdict["heldout_top_cca_shuffled"] + 0.2
    assert verdict["subspace_persists"] is False
    assert verdict["difference_ci_excludes_zero"] is True


def test_subspace_persistence_runs_with_a_real_confound_design():
    """The design is not decorative: the statistic must survive residualisation."""
    rng, x, y = _paired_blocks(n=300, seed=5)
    design = np.column_stack([np.repeat([1.0, 0.0], 150), np.repeat([0.0, 1.0], 150)])
    verdict = subspace_persistence(x, y, y, design, n_components=3, n_boot=25, seed=5)
    assert verdict["status"] == "scored"
    # identical blocks: the difference is exactly zero on every bootstrap draw
    assert verdict["paired_difference_true_minus_shuffled"] == pytest.approx(0.0, abs=1e-12)
    assert verdict["subspace_persists"] is True


def test_reconstruct_scores_matches_the_frozen_form_up_to_a_constant():
    """``build_shuffled_pbs_targets`` drops ``gene_mean``; the claim is that only a
    per-column additive constant is lost. That claim is asserted, not assumed."""
    rng = np.random.default_rng(6)
    expression = rng.normal(size=(120, 30))
    basis = rng.normal(size=(30, 5))
    gene_mean = rng.normal(size=30)
    frozen = (expression - gene_mean) @ basis
    rebuilt = reconstruct_scores(expression, basis)
    difference = rebuilt - frozen
    assert np.allclose(difference, difference[0], atol=1e-10)
    centred = lambda a: a - a.mean(axis=0, keepdims=True)
    np.testing.assert_allclose(centred(rebuilt), centred(frozen), atol=1e-10)
