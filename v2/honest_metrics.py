"""Confound-aware molecular-prompting metrics (F-R4).

The pooled (cross-cancer) Pearson correlation reported by most WSI->molecular
benchmarks is confounded: on our data ~46-49% of it is cross-cancer cohort
structure, and a size-matched random gene set already reaches ~0.30-0.32. The
honest, primary number is the **control-adjusted within-cancer specificity**:
the per-cancer (within-group) Pearson of a real signature, minus the same for a
matched random signature.

These are pure NumPy utilities so they can be unit-tested on CPU without the full
evaluation pipeline; the evaluator (`v21_evaluation.py`, `select_v21_profile.py`)
is expected to report `macro_group_pearson` alongside pooled Pearson and to gate
on the control-adjusted quantity. See v2/research/B2_implementation_and_audit.md
(R4) and the paper's metrics section.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError("pearson expects two 1-D arrays of equal length")
    if len(a) < 2:
        return float("nan")
    a = a - a.mean()
    b = b - b.mean()
    denom = np.sqrt((a * a).sum() * (b * b).sum())
    if denom <= 0:
        return float("nan")
    return float((a * b).sum() / denom)


def macro_group_pearson(prediction: Sequence[float], target: Sequence[float],
                        groups: Sequence, min_group: int = 3) -> float:
    """Mean of within-group (per-cancer) Pearson correlations.

    Groups with fewer than ``min_group`` paired observations are dropped (a
    two-patient cancer cannot yield a meaningful within-cancer correlation).
    Returns NaN when no group qualifies. This is the within-cancer number that
    strips the cross-cancer cohort shortcut out of pooled Pearson.
    """
    prediction = np.asarray(prediction, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    groups = np.asarray(groups)
    if not (len(prediction) == len(target) == len(groups)):
        raise ValueError("prediction, target and groups must be the same length")
    per_group: list[float] = []
    for value in np.unique(groups):
        mask = groups == value
        if int(mask.sum()) < min_group:
            continue
        r = _pearson(prediction[mask], target[mask])
        if not np.isnan(r):
            per_group.append(r)
    if not per_group:
        return float("nan")
    return float(np.mean(per_group))


def control_adjusted_specificity(prediction: Sequence[float], target: Sequence[float],
                                 control_prediction: Sequence[float], groups: Sequence,
                                 min_group: int = 3) -> float:
    """Within-cancer specificity of a real signature over a matched random control.

    ``macro_group_pearson`` of the real prediction minus that of a size/expression
    matched random-gene-set prediction against the same target and grouping. This
    is the paper's primary reported quantity (~+0.07 and method-invariant on our
    data).
    """
    real = macro_group_pearson(prediction, target, groups, min_group)
    control = macro_group_pearson(control_prediction, target, groups, min_group)
    if np.isnan(real) or np.isnan(control):
        return float("nan")
    return float(real - control)
