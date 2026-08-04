"""Self-tests for the new-modality lineage positive control.

A control is evidence only if it can be shown to react in both directions. Each
test below plants either a block that genuinely carries lineage or a block that
provably cannot, and asserts the verdict moves accordingly. A control that
always passes is a decoration; one that always fails is a broken gate.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.modality_block_control import (balanced_accuracy, evaluate_modality_block,
                                                        out_of_fold_lineage_prediction)


def _lineage_block(n_per_class: int = 60, n_classes: int = 5, n_features: int = 20,
                   separation: float = 3.0, seed: int = 0):
    rng = np.random.default_rng(seed)
    cancers = np.repeat([f"C{i}" for i in range(n_classes)], n_per_class)
    centres = rng.normal(size=(n_classes, n_features)) * separation
    features = np.repeat(centres, n_per_class, axis=0) + rng.normal(size=(len(cancers), n_features))
    return features, cancers


def test_balanced_accuracy_is_unweighted_so_a_rare_class_cannot_be_ignored():
    labels = np.array(["A"] * 90 + ["B"] * 10)
    all_a = np.array(["A"] * 100)
    value, per_class = balanced_accuracy(labels, all_a)
    # a pooled accuracy would be 0.90 here; the statistic we grade must be 0.50
    assert value == pytest.approx(0.5)
    assert per_class == {"A": pytest.approx(1.0), "B": pytest.approx(0.0)}


def test_out_of_fold_prediction_never_scores_a_patient_with_its_own_row():
    """Planted memorisable noise: with no shared structure, out-of-fold accuracy
    must collapse. An in-fold leak would make this near-perfect."""
    rng = np.random.default_rng(1)
    cancers = np.asarray([f"C{i % 4}" for i in range(160)])
    features = rng.normal(size=(160, 200))
    predicted = out_of_fold_lineage_prediction(features, cancers, seed=0)
    value, _ = balanced_accuracy(cancers, predicted)
    assert value < 0.45


def test_a_block_that_carries_lineage_passes():
    features, cancers = _lineage_block()
    verdict = evaluate_modality_block(features, cancers, n_permutations=30, seed=0)
    assert verdict["balanced_accuracy"] > 0.9
    assert verdict["passed"] is True
    assert verdict["above_measured_null"] and verdict["above_absolute_bar"]
    assert verdict["permutation_p"] <= verdict["permutation_resolution"] * 2


def test_a_block_that_carries_nothing_fails_and_the_null_is_measured_not_assumed():
    rng = np.random.default_rng(2)
    cancers = np.asarray([f"C{i % 5}" for i in range(300)])
    verdict = evaluate_modality_block(rng.normal(size=(300, 20)), cancers, n_permutations=30, seed=0)
    assert verdict["passed"] is False
    assert np.isfinite(verdict["null_balanced_accuracy_p95"])
    # the measured chance level is reported, and it is in the neighbourhood of 1/5
    assert 0.10 < verdict["null_balanced_accuracy_median"] < 0.35


def test_a_real_block_that_is_too_weak_fails_the_absolute_bar_even_if_it_beats_chance():
    """The two halves of the criterion are separable, and both are required."""
    features, cancers = _lineage_block(separation=0.12, n_per_class=90, seed=3)
    verdict = evaluate_modality_block(features, cancers, n_permutations=30, seed=0,
                                      minimum_balanced_accuracy=0.95)
    assert verdict["above_absolute_bar"] is False
    assert verdict["passed"] is False


def test_a_constant_block_is_refused_rather_than_scored_as_chance():
    cancers = np.asarray([f"C{i % 3}" for i in range(90)])
    verdict = evaluate_modality_block(np.zeros((90, 4)), cancers, n_permutations=5, seed=0)
    # a constant block predicts one class for everyone: balanced accuracy 1/3, not a pass
    assert verdict["passed"] is False
    assert verdict["balanced_accuracy"] == pytest.approx(1.0 / 3.0)


def test_non_finite_and_misaligned_blocks_are_refused():
    features, cancers = _lineage_block(n_per_class=10, n_classes=2, n_features=3)
    with pytest.raises(ValueError, match="non-finite"):
        broken = features.copy()
        broken[0, 0] = np.nan
        evaluate_modality_block(broken, cancers, n_permutations=2)
    with pytest.raises(ValueError, match="aligned"):
        evaluate_modality_block(features, cancers[:-1], n_permutations=2)


def test_the_verdict_is_reproducible_at_a_fixed_seed():
    features, cancers = _lineage_block(seed=4)
    first = evaluate_modality_block(features, cancers, n_permutations=10, seed=7)
    second = evaluate_modality_block(features, cancers, n_permutations=10, seed=7)
    assert first["balanced_accuracy"] == second["balanced_accuracy"]
    assert first["null_balanced_accuracy_p95"] == second["null_balanced_accuracy_p95"]
