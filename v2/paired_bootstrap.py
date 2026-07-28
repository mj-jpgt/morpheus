"""Paired patient and cancer-cluster bootstrap comparisons for frozen probes."""
from __future__ import annotations

from typing import Callable
import numpy as np


def pearson_metric(actual: np.ndarray, predicted: np.ndarray) -> float:
    actual, predicted = np.asarray(actual, float), np.asarray(predicted, float)
    valid = np.isfinite(actual) & np.isfinite(predicted)
    if valid.sum() < 3 or np.std(actual[valid]) < 1e-12 or np.std(predicted[valid]) < 1e-12:
        return float("nan")
    return float(np.corrcoef(actual[valid], predicted[valid])[0, 1])


def _interval(values: np.ndarray) -> dict[str, float]:
    values = values[np.isfinite(values)]
    return {"n_valid": int(len(values)), "mean_delta": float(np.mean(values)) if len(values) else float("nan"),
            "ci95_low": float(np.quantile(values, .025)) if len(values) else float("nan"),
            "ci95_high": float(np.quantile(values, .975)) if len(values) else float("nan"),
            "p_improve": float(np.mean(values > 0)) if len(values) else float("nan")}


def paired_bootstrap_difference(
    metric: Callable[[np.ndarray, np.ndarray], float], actual: np.ndarray, teacher: np.ndarray,
    challenger: np.ndarray, *, cancers: np.ndarray | None = None, repeats: int = 2000,
    seed: int = 42, mode: str = "patient",
) -> dict[str, float | str]:
    """Paired resampling of a challenger-minus-teacher metric difference.

    ``patient`` resamples patients IID. ``cancer`` samples cancer clusters and
    then patients within each selected cluster, preserving the paired rows and
    reflecting held-out-cancer uncertainty.
    """
    y, ref, candidate = map(lambda value: np.asarray(value), (actual, teacher, challenger))
    if not (len(y) == len(ref) == len(candidate)) or y.ndim != 1 or ref.ndim != 1 or candidate.ndim != 1:
        raise ValueError("actual, teacher, and challenger must be equal-length vectors")
    if repeats < 100:
        raise ValueError("repeats must be at least 100")
    if mode not in {"patient", "cancer"}:
        raise ValueError("mode must be 'patient' or 'cancer'")
    cancer = None if cancers is None else np.asarray(cancers).astype(str)
    if mode == "cancer" and (cancer is None or len(cancer) != len(y)):
        raise ValueError("cancer bootstrap requires one cancer label per patient")
    rng, deltas = np.random.default_rng(seed), []
    labels = None if cancer is None else np.unique(cancer)
    for _ in range(repeats):
        if mode == "patient":
            indices = rng.integers(0, len(y), len(y))
        else:
            assert labels is not None and cancer is not None
            sampled = rng.choice(labels, len(labels), replace=True)
            chunks = [rng.choice(np.flatnonzero(cancer == label), int((cancer == label).sum()), replace=True) for label in sampled]
            indices = np.concatenate(chunks)
        deltas.append(metric(y[indices], candidate[indices]) - metric(y[indices], ref[indices]))
    result = _interval(np.asarray(deltas, float))
    result.update({"mode": mode, "repeats": repeats, "point_delta": metric(y, candidate) - metric(y, ref)})
    return result


def paired_patient_and_cancer_bootstrap(metric, actual, teacher, challenger, cancers, *, repeats: int = 2000, seed: int = 42) -> dict[str, object]:
    """Return both required uncertainty views under deterministic distinct seeds."""
    return {"patient": paired_bootstrap_difference(metric, actual, teacher, challenger, repeats=repeats, seed=seed, mode="patient"),
            "cancer": paired_bootstrap_difference(metric, actual, teacher, challenger, cancers=cancers, repeats=repeats, seed=seed + 1, mode="cancer")}
