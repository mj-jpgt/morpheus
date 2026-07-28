from __future__ import annotations

import numpy as np
import pandas as pd

from morpheus.v2.survival_evaluation import _bootstrap_cindex, evaluate_coxnet_endpoint, paired_cindex_bootstrap


def _inputs() -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    n = 40
    ids = [f"TCGA-AA-{i:04d}" for i in range(n)]
    metadata = pd.DataFrame({"patient_id": ids, "cancer": np.repeat(["A", "B", "C", "D"], 10),
                             "split": ["train"] * 10 + ["val"] * 10 + ["test"] * 20})
    outcomes = pd.DataFrame({"patient_id": ids, "os_time_days": np.arange(1, n + 1, dtype=float),
                             "os_event": np.array([1, 0] * 20, dtype=float), "pfi_time_days": np.nan,
                             "pfi_event": np.nan})
    return np.random.default_rng(2).normal(size=(n, 6)), metadata, outcomes


def test_survival_availability_gate_is_explicit() -> None:
    features, metadata, outcomes = _inputs()
    outcomes.loc[outcomes.index[-20:], "os_event"] = 0.0
    result = evaluate_coxnet_endpoint(features, metadata, outcomes, endpoint="os", method="m", representation_state="wsi", min_test_events=5)
    assert result.status == "unavailable_insufficient_endpoint_coverage_or_events"
    assert result.rows[0]["note"] == result.status


def test_survival_bootstrap_returns_both_cluster_aware_intervals() -> None:
    time = np.arange(1, 41, dtype=float); event = np.array([1, 0] * 20, dtype=bool); risk = -time
    patient = _bootstrap_cindex(time, event, risk, np.repeat(["A", "B", "C", "D"], 10), repeats=100, seed=4, mode="patient")
    cancer = _bootstrap_cindex(time, event, risk, np.repeat(["A", "B", "C", "D"], 10), repeats=100, seed=5, mode="cancer")
    assert patient["point"] > .95 and cancer["n_valid"] > 0


def test_paired_cindex_bootstrap_preserves_method_pairing() -> None:
    payload = paired_cindex_bootstrap(
        np.asarray([1, 2, 3, 4, 5, 6], dtype=float),
        np.asarray([1, 1, 1, 0, 1, 0], dtype=bool),
        np.asarray([.10, .20, .30, .40, .50, .60]),
        np.asarray([.20, .10, .40, .30, .60, .50]),
        np.asarray(["A", "A", "A", "B", "B", "B"]),
        repeats=100,
    )
    assert set(payload) == {"patient", "cancer"}
    assert payload["patient"]["n_valid"] > 0
    assert "point_delta" in payload["cancer"]
