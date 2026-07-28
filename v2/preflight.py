"""Strict primary-protocol checks shared by local and Lambda V2 launches."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


def sha256_json(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def validate_runtime_split(split_path: str | Path, patient_ids, cancers, labels,
                           expected_development_cancers: int = 11,
                           expected_heldout_cancers: int = 22) -> dict[str, object]:
    """Reject malformed or non-held-out-cancer primary splits.

    The caller passes the *actual paired cohort* after data loading, preventing
    a split manifest from claiming patients absent from the training tensors.
    """
    payload = json.loads(Path(split_path).read_text(encoding="utf-8"))
    declared = payload.get("patient_ids")
    if not isinstance(declared, dict) or not {"train", "val", "test"}.issubset(declared):
        raise ValueError("primary V2 split needs patient_ids for train, val, and test")
    sets = {name: set(map(str, declared[name])) for name in ("train", "val", "test")}
    if any(not values for values in sets.values()):
        raise ValueError("primary V2 split cannot contain an empty partition")
    if sets["train"] & sets["val"] or sets["train"] & sets["test"] or sets["val"] & sets["test"]:
        raise ValueError("patient appears in more than one primary split")
    actual_ids = list(map(str, patient_ids))
    if len(set(actual_ids)) != len(actual_ids):
        raise ValueError("loaded paired cohort contains duplicate patient identifiers")
    actual = {str(pid): str(label) for pid, label in zip(patient_ids, labels)}
    if any(actual.get(pid) != partition for partition, rows in sets.items() for pid in rows if pid in actual):
        raise ValueError("loaded cohort labels disagree with split manifest")
    declared_universe = set().union(*sets.values())
    if not set(actual).issubset(declared_universe):
        raise ValueError("loaded paired cohort contains unassigned patients")
    if declared_universe != set(actual):
        missing = len(declared_universe - set(actual)); extra = len(set(actual) - declared_universe)
        raise ValueError(f"primary split/cohort mismatch (missing={missing}, extra={extra}); rebuild an explicit paired split")
    cancer_by_patient = {str(pid): str(cancer) for pid, cancer in zip(patient_ids, cancers)}
    train_cancers = {cancer_by_patient[pid] for pid in sets["train"] if pid in cancer_by_patient}
    val_cancers = {cancer_by_patient[pid] for pid in sets["val"] if pid in cancer_by_patient}
    test_cancers = {cancer_by_patient[pid] for pid in sets["test"] if pid in cancer_by_patient}
    if not test_cancers or test_cancers & train_cancers or test_cancers & val_cancers:
        raise ValueError("outer-test cancer types must be disjoint from development cancers")
    if not val_cancers.issubset(train_cancers):
        raise ValueError("validation cancers must be a subset of development cancers")
    expected = payload.get("actual_train_test_cancer_counts")
    if expected and list(expected) != [len(train_cancers), len(test_cancers)]:
        raise ValueError("split cancer counts disagree with manifest")
    if len(train_cancers) != expected_development_cancers or len(test_cancers) != expected_heldout_cancers:
        raise ValueError(f"strict primary protocol requires {expected_development_cancers} development and {expected_heldout_cancers} held-out cancers, got {len(train_cancers)}/{len(test_cancers)}")
    return {
        "split_path": str(Path(split_path)), "split_digest": sha256_json({k: sorted(v) for k, v in sets.items()}),
        "cohort_digest": sha256_json({"patient_ids": list(map(str, patient_ids)), "cancers": list(map(str, cancers)), "labels": list(map(str, labels))}),
        "counts": {name: len(rows) for name, rows in sets.items()},
        "cancer_counts": {"development": len(train_cancers), "heldout": len(test_cancers)},
        "protocol": f"heldout_cancer_{expected_development_cancers}v{expected_heldout_cancers}",
    }
