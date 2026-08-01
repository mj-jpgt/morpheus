"""Leakage-safe survival probes for frozen MORPHEUS representations.

The evaluator never treats clinical fields as input features.  It fits all
scaling and Coxnet alpha selection within development cancers, exports one
held-out risk per participant, and represents every unsupported endpoint with
an explicit status instead of silently omitting survival.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from .tasks import harrell_cindex


DEFAULT_ALPHA_GRID = tuple(float(x) for x in np.logspace(-4, 1, 20))


@dataclass(frozen=True)
class SurvivalEvaluation:
    rows: list[dict[str, Any]]
    risks: pd.DataFrame
    selected_alpha: float | None
    status: str


def _structured_target(event: np.ndarray, time: np.ndarray) -> np.ndarray:
    return np.asarray([(bool(e), float(t)) for e, t in zip(event, time)], dtype=[("event", "?"), ("time", "<f8")])


def _bootstrap_cindex(time: np.ndarray, event: np.ndarray, risk: np.ndarray, cancers: np.ndarray, *, repeats: int, seed: int, mode: str) -> dict[str, float | int | str]:
    if repeats < 100:
        raise ValueError("bootstrap repeats must be at least 100")
    rng = np.random.default_rng(seed)
    values: list[float] = []
    labels = np.unique(cancers)
    for _ in range(repeats):
        if mode == "patient":
            index = rng.integers(0, len(time), len(time))
        elif mode == "cancer":
            sampled = rng.choice(labels, len(labels), replace=True)
            pieces = [rng.choice(np.flatnonzero(cancers == label), int((cancers == label).sum()), replace=True) for label in sampled]
            index = np.concatenate(pieces)
        else:
            raise ValueError("mode must be patient or cancer")
        values.append(harrell_cindex(time[index], event[index], risk[index]))
    values = np.asarray(values, float); values = values[np.isfinite(values)]
    return {"mode": mode, "repeats": repeats, "n_valid": int(len(values)),
            "point": harrell_cindex(time, event, risk),
            "ci95_low": float(np.quantile(values, .025)) if len(values) else float("nan"),
            "ci95_high": float(np.quantile(values, .975)) if len(values) else float("nan")}


def paired_cindex_bootstrap(
    time: np.ndarray,
    event: np.ndarray,
    teacher_risk: np.ndarray,
    challenger_risk: np.ndarray,
    cancers: np.ndarray,
    *,
    repeats: int = 2000,
    seed: int = 42,
) -> dict[str, dict[str, float | int | str]]:
    """Paired challenger-minus-teacher C-index intervals.

    Both risk vectors are evaluated on precisely the same held-out patients.
    We report an IID patient bootstrap and a cancer-cluster bootstrap rather
    than comparing two independently bootstrapped C-indices, which would
    exaggerate uncertainty and discard the paired design.
    """
    time = np.asarray(time, float)
    event = np.asarray(event, bool)
    teacher_risk = np.asarray(teacher_risk, float)
    challenger_risk = np.asarray(challenger_risk, float)
    cancers = np.asarray(cancers).astype(str)
    if not (len(time) == len(event) == len(teacher_risk) == len(challenger_risk) == len(cancers)):
        raise ValueError("paired survival bootstrap inputs must have equal length")
    if repeats < 100:
        raise ValueError("repeats must be at least 100")

    def summarize(values: list[float], point: float, mode: str) -> dict[str, float | int | str]:
        finite = np.asarray(values, float)
        finite = finite[np.isfinite(finite)]
        return {
            "mode": mode,
            "repeats": repeats,
            "n_valid": int(len(finite)),
            "point_delta": point,
            "mean_delta": float(np.mean(finite)) if len(finite) else float("nan"),
            "ci95_low": float(np.quantile(finite, .025)) if len(finite) else float("nan"),
            "ci95_high": float(np.quantile(finite, .975)) if len(finite) else float("nan"),
            "p_improve": float(np.mean(finite > 0)) if len(finite) else float("nan"),
        }

    point = harrell_cindex(time, event, challenger_risk) - harrell_cindex(time, event, teacher_risk)
    labels = np.unique(cancers)
    payload: dict[str, dict[str, float | int | str]] = {}
    for offset, mode in enumerate(("patient", "cancer")):
        rng = np.random.default_rng(seed + offset)
        values: list[float] = []
        for _ in range(repeats):
            if mode == "patient":
                index = rng.integers(0, len(time), len(time))
            else:
                selected = rng.choice(labels, len(labels), replace=True)
                index = np.concatenate([
                    rng.choice(np.flatnonzero(cancers == label), int((cancers == label).sum()), replace=True)
                    for label in selected
                ])
            values.append(
                harrell_cindex(time[index], event[index], challenger_risk[index])
                - harrell_cindex(time[index], event[index], teacher_risk[index])
            )
        payload[mode] = summarize(values, point, mode)
    return payload


def _macro_cindex(time: np.ndarray, event: np.ndarray, risk: np.ndarray, cancers: np.ndarray, *, min_samples: int, min_events: int) -> float:
    values = []
    for cancer in np.unique(cancers):
        keep = cancers == cancer
        if keep.sum() >= min_samples and event[keep].sum() >= min_events:
            value = harrell_cindex(time[keep], event[keep], risk[keep])
            if np.isfinite(value):
                values.append(value)
    return float(np.mean(values)) if values else float("nan")


def _select_alpha(features: np.ndarray, time: np.ndarray, event: np.ndarray, development: np.ndarray, cancers: np.ndarray,
                  *, alpha_grid: Sequence[float], folds: int, l1_ratio: float) -> float | None:
    try:
        from sksurv.linear_model import CoxnetSurvivalAnalysis
    except ImportError:
        return None
    keep = development & np.isfinite(time) & np.isfinite(event) & (time > 0)
    groups = cancers[keep]
    if keep.sum() < 8 or len(np.unique(groups)) < 2:
        return None
    scores = {float(alpha): [] for alpha in alpha_grid}
    splitter = GroupKFold(n_splits=min(folds, len(np.unique(groups))))
    x, t, e = features[keep], time[keep], event[keep]
    for train_rows, valid_rows in splitter.split(x, groups=groups):
        if e[train_rows].sum() < 2 or e[valid_rows].sum() < 1:
            continue
        scaler = StandardScaler().fit(x[train_rows])
        try:
            model = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=np.asarray(alpha_grid)).fit(
                scaler.transform(x[train_rows]), _structured_target(e[train_rows], t[train_rows]))
            xv = scaler.transform(x[valid_rows])
            for alpha in alpha_grid:
                prediction = model.predict(xv, alpha=float(alpha))
                scores[float(alpha)].append(harrell_cindex(t[valid_rows], e[valid_rows].astype(bool), prediction))
        except Exception:
            continue
    viable = [(alpha, np.nanmean(values)) for alpha, values in scores.items() if values and np.isfinite(np.nanmean(values))]
    return max(viable, key=lambda item: (item[1], -item[0]))[0] if viable else None


def evaluate_coxnet_endpoint(
    features: np.ndarray,
    metadata: pd.DataFrame,
    outcomes: pd.DataFrame,
    *,
    endpoint: str,
    method: str,
    representation_state: str,
    alpha_grid: Sequence[float] = DEFAULT_ALPHA_GRID,
    l1_ratio: float = 0.5,
    inner_folds: int = 3,
    bootstrap_repeats: int = 2000,
    seed: int = 42,
    min_development_patients: int = 30,
    min_development_events: int = 10,
    min_test_patients: int = 20,
    min_test_events: int = 5,
) -> SurvivalEvaluation:
    """Evaluate a single CDR endpoint with cancer-grouped inner CV.

    ``metadata`` must contain `patient_id`, `cancer`, `split` in exactly the
    same order as `features`. `outcomes` must contain `patient_id` and
    `{endpoint}_time_days` / `{endpoint}_event` from the canonical CDR builder.
    """
    features = np.asarray(features, dtype=float)
    required_meta = {"patient_id", "cancer", "split"}
    required_outcomes = {"patient_id", f"{endpoint}_time_days", f"{endpoint}_event"}
    if not required_meta.issubset(metadata) or not required_outcomes.issubset(outcomes):
        raise ValueError("metadata/outcomes lack required endpoint columns")
    if len(features) != len(metadata) or metadata.patient_id.astype(str).duplicated().any() or outcomes.patient_id.astype(str).duplicated().any():
        raise ValueError("features and metadata must be aligned with unique patient IDs")
    labels = metadata[["patient_id", "cancer", "split"]].copy()
    labels["patient_id"] = labels.patient_id.astype(str)
    merged = labels.merge(outcomes[["patient_id", f"{endpoint}_time_days", f"{endpoint}_event"]].assign(patient_id=outcomes.patient_id.astype(str)),
                          on="patient_id", how="left", validate="one_to_one")
    time = pd.to_numeric(merged[f"{endpoint}_time_days"], errors="coerce").to_numpy(float)
    event = pd.to_numeric(merged[f"{endpoint}_event"], errors="coerce").to_numpy(float)
    development = merged.split.astype(str).isin(("train", "val")).to_numpy()
    test = merged.split.astype(str).eq("test").to_numpy()
    valid = np.isfinite(time) & np.isfinite(event) & (time > 0)
    base = {"method": method, "representation_state": representation_state, "task": f"survival_{endpoint}"}
    counts = {"development_patients": int((development & valid).sum()), "development_events": int((development & valid & event.astype(bool)).sum()),
              "test_patients": int((test & valid).sum()), "test_events": int((test & valid & event.astype(bool)).sum())}
    if counts["development_patients"] < min_development_patients or counts["development_events"] < min_development_events or counts["test_patients"] < min_test_patients or counts["test_events"] < min_test_events:
        return SurvivalEvaluation([{**base, "metric": "status", "value": np.nan, "note": "unavailable_insufficient_endpoint_coverage_or_events", **counts}], pd.DataFrame(), None, "unavailable_insufficient_endpoint_coverage_or_events")
    alpha = _select_alpha(features, time, event, development, merged.cancer.astype(str).to_numpy(), alpha_grid=alpha_grid, folds=inner_folds, l1_ratio=l1_ratio)
    if alpha is None:
        return SurvivalEvaluation([{**base, "metric": "status", "value": np.nan, "note": "unavailable_sksurv_or_inner_cv_fit_failure", **counts}], pd.DataFrame(), None, "unavailable_sksurv_or_inner_cv_fit_failure")
    try:
        from sksurv.linear_model import CoxnetSurvivalAnalysis
        scaler = StandardScaler().fit(features[development & valid])
        model = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=np.asarray(alpha_grid)).fit(
            scaler.transform(features[development & valid]), _structured_target(event[development & valid], time[development & valid]))
        test_valid = test & valid
        risks = model.predict(scaler.transform(features[test_valid]), alpha=alpha)
    except Exception as exc:
        return SurvivalEvaluation([{**base, "metric": "status", "value": np.nan, "note": f"unavailable_coxnet_final_fit_failed:{type(exc).__name__}", **counts}], pd.DataFrame(), alpha, "unavailable_coxnet_final_fit_failed")
    test_time, test_event, test_cancer = time[test_valid], event[test_valid].astype(bool), merged.cancer.astype(str).to_numpy()[test_valid]
    overall = harrell_cindex(test_time, test_event, risks)
    macro = _macro_cindex(test_time, test_event, risks, test_cancer, min_samples=5, min_events=1)
    patient_boot = _bootstrap_cindex(test_time, test_event, risks, test_cancer, repeats=bootstrap_repeats, seed=seed, mode="patient")
    cancer_boot = _bootstrap_cindex(test_time, test_event, risks, test_cancer, repeats=bootstrap_repeats, seed=seed + 1, mode="cancer")
    rows = [{**base, "metric": "harrell_cindex", "value": overall, "note": "development_only_grouped_inner_cv", "selected_alpha": alpha, **counts},
            {**base, "metric": "macro_cancer_harrell_cindex", "value": macro, "note": "test_cancers_with_at_least_5_patients_and_1_event", "selected_alpha": alpha, **counts}]
    for prefix, payload in (("patient_bootstrap", patient_boot), ("cancer_bootstrap", cancer_boot)):
        for key in ("ci95_low", "ci95_high", "point", "n_valid"):
            rows.append({**base, "metric": f"{prefix}_{key}", "value": payload[key], "note": "heldout_test_only", "selected_alpha": alpha, **counts})
    risks_table = pd.DataFrame({"patient_id": merged.loc[test_valid, "patient_id"].to_numpy(), "cancer": test_cancer,
                                "endpoint": endpoint, "method": method, "representation_state": representation_state,
                                "time_days": test_time, "event": test_event.astype(int), "risk": risks})
    return SurvivalEvaluation(rows, risks_table, alpha, "ok")
