"""Unit tests for the confound-aware molecular-prompting metrics (F-R4)."""
from __future__ import annotations

import math

import numpy as np

from morpheus.v2.honest_metrics import (control_adjusted_specificity,
                                        macro_group_pearson)


def test_macro_group_is_mean_of_within_group_pearsons() -> None:
    # Two cancers, each with a known within-group correlation.
    pred = np.array([1.0, 2.0, 3.0, 4.0, 10.0, 20.0, 30.0, 40.0])
    target = np.array([1.0, 2.0, 3.0, 4.0, 40.0, 30.0, 20.0, 10.0])  # group A: r=+1, group B: r=-1
    groups = np.array(["A", "A", "A", "A", "B", "B", "B", "B"])
    macro = macro_group_pearson(pred, target, groups)
    assert math.isclose(macro, 0.0, abs_tol=1e-9)  # mean of (+1, -1)


def test_macro_strips_cross_group_shortcut() -> None:
    # A pure cross-group signal (constant offset per group, flat within group) has
    # high pooled Pearson but ~zero within-group Pearson.
    rng = np.random.default_rng(0)
    groups = np.array(["A"] * 30 + ["B"] * 30)
    base = np.where(groups == "A", 0.0, 10.0)
    pred = base + rng.normal(scale=0.01, size=60)
    target = base + rng.normal(scale=0.01, size=60)
    pooled = np.corrcoef(pred, target)[0, 1]
    macro = macro_group_pearson(pred, target, groups)
    assert pooled > 0.99          # confounded pooled metric looks great
    assert abs(macro) < 0.5       # within-cancer reveals there is little real signal


def test_small_groups_are_dropped_and_all_small_returns_nan() -> None:
    pred = np.array([1.0, 2.0, 3.0, 4.0])
    target = np.array([1.0, 2.0, 3.0, 4.0])
    groups = np.array(["A", "A", "B", "B"])  # every group has <3 pairs
    assert math.isnan(macro_group_pearson(pred, target, groups, min_group=3))
    # With one qualifying group, only that group counts.
    pred2 = np.array([1.0, 2.0, 3.0, 9.0])
    target2 = np.array([1.0, 2.0, 3.0, 1.0])
    groups2 = np.array(["A", "A", "A", "B"])
    assert not math.isnan(macro_group_pearson(pred2, target2, groups2, min_group=3))


def test_control_adjusted_specificity_subtracts_random_baseline() -> None:
    rng = np.random.default_rng(1)
    groups = np.array(["A"] * 20 + ["B"] * 20)
    target = rng.normal(size=40)
    real = target + rng.normal(scale=0.3, size=40)      # correlated with target
    control = rng.normal(size=40)                         # random gene set
    spec = control_adjusted_specificity(real, target, control, groups)
    assert spec > 0.2                                     # real beats random within-cancer
    # A "real" signature that is actually random has ~zero specificity.
    null_spec = control_adjusted_specificity(rng.normal(size=40), target, control, groups)
    assert abs(null_spec) < 0.4
