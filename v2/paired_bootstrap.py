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


def paired_multivariate_patient_and_cancer_bootstrap(metric, actual, teacher, challenger, cancers, *, repeats: int = 2000, seed: int = 42) -> dict[str, object]:
    """Paired row bootstrap for a multivariate channel statistic.

    All inputs are [patient, feature] matrices and are resampled along axis 0
    with the same indices.  This is intentionally separate from the scalar
    helper above: accepting matrices there would silently make one-dimensional
    assumptions inside unrelated probe comparisons.
    """
    y, ref, candidate = (np.asarray(value, float) for value in (actual, teacher, challenger))
    if any(value.ndim != 2 for value in (y, ref, candidate)) or not (len(y) == len(ref) == len(candidate)):
        raise ValueError("multivariate bootstrap needs aligned finite [patient,feature] matrices")
    if not all(np.isfinite(value).all() for value in (y, ref, candidate)):
        raise ValueError("multivariate bootstrap inputs must be finite")
    labels = np.asarray(cancers).astype(str)
    if len(labels) != len(y) or len(np.unique(labels)) < 2 or repeats < 100:
        raise ValueError("multivariate cancer bootstrap needs >=2 cancer labels and >=100 repeats")

    def one(mode: str, local_seed: int) -> dict[str, float]:
        rng, values = np.random.default_rng(local_seed), []
        groups = np.unique(labels)
        for _ in range(repeats):
            if mode == "patient":
                index = rng.integers(0, len(y), len(y))
            else:
                sampled = rng.choice(groups, len(groups), replace=True)
                index = np.concatenate([rng.choice(np.flatnonzero(labels == group), int((labels == group).sum()), replace=True)
                                        for group in sampled])
            values.append(float(metric(y[index], candidate[index]) - metric(y[index], ref[index])))
        result = _interval(np.asarray(values, float))
        result.update({"mode": mode, "repeats": repeats, "point_delta": float(metric(y, candidate) - metric(y, ref))})
        return result
    return {"patient": one("patient", seed), "cancer": one("cancer", seed + 1)}
