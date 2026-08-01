"""Strict primary-protocol checks shared by local and Lambda V2 launches."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


def sha256_json(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def split_universe(split_path: str | Path) -> set[str]:
    payload = json.loads(Path(split_path).read_text(encoding="utf-8"))
    declared = payload.get("patient_ids")
    if not isinstance(declared, dict) or not {"train", "val", "test"}.issubset(declared):
        raise ValueError("primary V2 split needs patient_ids for train, val, and test")
    return {str(pid) for name in ("train", "val", "test") for pid in declared[name]}


def split_file_digest(split_path: str | Path) -> str:
    """The same digest `validate_runtime_split` records, computable from the file alone."""
    payload = json.loads(Path(split_path).read_text(encoding="utf-8"))
    declared = payload["patient_ids"]
    return sha256_json({name: sorted(map(str, declared[name])) for name in ("train", "val", "test")})


def restrict_cohort_to_split(data, split_path: str | Path,
                             max_dropped_fraction: float = 0.05) -> tuple[object, dict[str, object]]:
    """Make the split file the authoritative cohort lookup.

    ``validate_runtime_split`` requires the loaded cohort and the declared split
    to be the *same* set, which is what stops a stale manifest from claiming
    patients that no longer have features.  But the loader's cohort also grows
    on its own -- a newly extracted slide adds a patient nobody assigned -- and
    a downstream artifact (PBS codes) may legitimately cover fewer patients than
    the loader returns.  This helper resolves that in the only safe direction:

    * every patient the split **declares** must exist in the loaded cohort, or
      this raises.  That is the protection ``validate_runtime_split`` provides
      and it is preserved exactly.
    * loaded patients the split **omits** are dropped, by name, with the count
      recorded.  They were already labelled ``excluded`` by ``_split_labels``,
      so they never entered the train-fold standardisation; dropping them here
      is therefore equivalent to never having loaded them, not a second
      normalisation population.

    The drop is also BOUNDED.  Without a ceiling this flag would silently convert
    "you are using a stale split" -- previously a hard failure -- into a quiet
    training run on whatever subset the stale file happened to name.  Beyond
    `max_dropped_fraction` the split is not narrower, it is wrong.

    Callers must record the returned manifest.  A silent cohort change is the
    failure this whole module exists to prevent.
    """
    import dataclasses

    universe = split_universe(split_path)
    actual = [str(patient) for patient in data.patient_ids]
    if len(set(actual)) != len(actual):
        raise ValueError("loaded paired cohort contains duplicate patient identifiers")
    absent = sorted(universe - set(actual))
    if absent:
        raise ValueError(
            f"split declares {len(absent)} patients the loaded cohort does not contain, so the "
            f"split is stale rather than merely narrower; examples={absent[:5]}")
    keep = np.asarray([patient in universe for patient in actual], dtype=bool)
    dropped = [patient for patient, flag in zip(actual, keep) if not flag]
    fraction = len(dropped) / max(len(actual), 1)
    if fraction > max_dropped_fraction:
        raise ValueError(
            f"split omits {len(dropped)} of {len(actual)} loaded patients ({fraction:.3f} > "
            f"{max_dropped_fraction:.3f}); that is a different cohort, not a narrower one. "
            f"Rebuild the split with morpheus.v2.build_paired_split.")
    # `dataclasses.replace` reconstructs the object, so any attribute attached
    # AFTER loading (the _v2_* target matrices, teacher anchors) would be lost
    # without warning. Restriction must therefore happen before attachment.
    attached = sorted(name for name in vars(data) if name.startswith("_"))
    if attached:
        raise ValueError(f"restrict the cohort before attaching targets; already attached: {attached}")
    manifest = {
        "enabled": True,
        "loaded_cohort_count": len(actual),
        "retained_count": int(keep.sum()),
        "dropped_count": len(dropped),
        "dropped_patients": dropped,
        "dropped_fraction": fraction,
        "max_dropped_fraction": float(max_dropped_fraction),
        "rationale": "split file is the authoritative cohort; dropped rows were already "
                     "'excluded' and never entered train-fold standardisation",
    }
    if keep.all():
        return data, manifest
    rows = np.where(keep)[0]
    changed: dict[str, object] = {}
    # Column metadata is never row-aligned; filtering it by a length coincidence
    # would silently truncate the target names.
    not_row_aligned = {"hallmark_names", "clinical_names", "hoptimus_store"}
    for field in dataclasses.fields(data):
        if field.name in not_row_aligned:
            continue
        value = getattr(data, field.name)
        if isinstance(value, np.ndarray) and value.shape[:1] == (len(actual),):
            changed[field.name] = value[keep]
        elif isinstance(value, list) and len(value) == len(actual):
            changed[field.name] = [value[index] for index in rows]
    if "patient_ids" not in changed or len(changed["patient_ids"]) != int(keep.sum()):
        raise AssertionError("cohort restriction failed to filter the patient axis")
    return dataclasses.replace(data, **changed), manifest


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
