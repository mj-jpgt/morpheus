"""Focused tests for leakage-safe RNA target construction."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.discovery_targets import (
    TargetBuildError,
    build_rna_target_bundle,
    fit_train_only_standardizer,
    read_gmt,
)


def test_bundle_keeps_canonical_order_and_per_target_masks() -> None:
    rna = pd.DataFrame(
        {
            "patient_id": ["p1", "p2", "p3"],
            "G1": [1.0, 2.0, 3.0],
            "G2": [3.0, np.nan, 5.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "patient_id": ["p1", "p2", "p3", "p4"],
            "split": ["train", "train", "test", "test"],
            "cancer": ["A", "B", "C", "C"],
        }
    )
    bundle = build_rna_target_bundle(
        rna,
        {"present": ["G1", "G2"], "absent": ["NOT_IN_MATRIX"]},
        metadata,
        canonical_patient_ids=["p4", "p2", "p1", "p3"],
    )

    assert bundle.raw_scores["patient_id"].tolist() == ["p4", "p2", "p1", "p3"]
    assert bundle.availability["present"].tolist() == [False, True, True, True]
    assert bundle.availability["absent"].tolist() == [False, False, False, False]
    assert bundle.raw_scores["absent"].isna().all()
    assert bundle.manifest["available_patients_by_target"] == {"present": 3, "absent": 0}


def test_standardization_does_not_fit_or_center_an_unseen_test_cancer() -> None:
    scores = pd.DataFrame(
        {"patient_id": ["p1", "p2", "p3", "p4"], "programme": [1.0, 3.0, 100.0, 2.0]}
    )
    metadata = pd.DataFrame(
        {
            "patient_id": ["p1", "p2", "p3", "p4"],
            "split": ["train", "train", "test", "val"],
            "cancer": ["A", "B", "C", "A"],
        }
    )
    fitted = fit_train_only_standardizer(scores, metadata)
    transformed = fitted.transform(scores, metadata).set_index("patient_id")["programme"]

    # Training cancer offsets are fit from p1/p2 only. Test cancer C is unseen,
    # so its extreme test value cannot be centered by test-derived information.
    assert fitted.global_mean["programme"] == pytest.approx(2.0)
    assert fitted.global_std["programme"] == pytest.approx(1.0)
    assert fitted.cancer_offset["programme"] == {"A": -1.0, "B": 1.0}
    assert transformed["p1"] == pytest.approx(0.0)
    assert transformed["p2"] == pytest.approx(0.0)
    assert transformed["p3"] == pytest.approx(98.0)
    assert transformed["p4"] == pytest.approx(1.0)


def test_gmt_parser_and_invalid_insufficient_training_target() -> None:
    assert read_gmt({"sig": ["A", "A", "B"]}) == {"sig": ("A", "B")}
    scores = pd.DataFrame({"patient_id": ["p1", "p2"], "x": [1.0, np.nan]})
    metadata = pd.DataFrame(
        {"patient_id": ["p1", "p2"], "split": ["train", "test"], "cancer": ["A", "B"]}
    )
    with pytest.raises(TargetBuildError, match="fewer than two"):
        fit_train_only_standardizer(scores, metadata)
