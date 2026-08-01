"""The cohort is the experiment. These tests exist because a silent cohort
change is unfalsifiable after the fact: every number still computes, the gates
still pass, and nothing in the artifact says the patients moved.

Covered:
  * the maximal split inherits the held-out protocol instead of re-drawing it;
  * no source assignment is ever moved;
  * newly available patients land on the correct side of the boundary;
  * `restrict_cohort_to_split` filters only the patient axis, and fail-closes
    when the split is stale rather than merely narrower.
"""
from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from morpheus.v2.build_paired_split import build_split
from morpheus.v2.preflight import restrict_cohort_to_split, validate_runtime_split


DEV = [f"DEV{index:03d}" for index in range(60)]
HELD = [f"HLD{index:03d}" for index in range(40)]


def _cohort(extra_dev: int = 0, extra_held: int = 0):
    """11 development and 21 held-out cancers, matching the real V2 protocol."""
    patients, cancers = [], []
    for index, patient in enumerate(DEV):
        patients.append(patient); cancers.append(f"DEVCANCER{index % 11}")
    for index, patient in enumerate(HELD):
        patients.append(patient); cancers.append(f"HELDCANCER{index % 21}")
    for index in range(extra_dev):
        patients.append(f"NEWDEV{index:03d}"); cancers.append(f"DEVCANCER{index % 11}")
    for index in range(extra_held):
        patients.append(f"NEWHLD{index:03d}"); cancers.append(f"HELDCANCER{index % 21}")
    return patients, cancers


def _source_split(tmp_path: Path) -> Path:
    payload = {"patient_ids": {"train": DEV[:50], "val": DEV[50:], "test": HELD}}
    path = tmp_path / "source_split.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_maximal_split_absorbs_new_patients_without_moving_old_ones(tmp_path):
    source = _source_split(tmp_path)
    patients, cancers = _cohort(extra_dev=22, extra_held=33)
    result = build_split(patient_ids=patients, cancers=cancers, source_split=source)

    assert result["paired_cohort_count"] == len(patients)
    assert result["added_patient_count"] == 55
    sets = result["patient_ids"]
    # Not one source patient may change partition.
    assert set(DEV[:50]).issubset(sets["train"])
    assert set(DEV[50:]).issubset(sets["val"])
    assert set(HELD).issubset(sets["test"])
    # Every new held-out-cancer patient goes to test; no new patient is dropped.
    assert {f"NEWHLD{index:03d}" for index in range(33)}.issubset(sets["test"])
    assert set(sets["train"]) | set(sets["val"]) | set(sets["test"]) == set(patients)
    assert not (set(sets["train"]) & set(sets["val"]))


def test_new_development_patients_never_land_in_test(tmp_path):
    source = _source_split(tmp_path)
    patients, cancers = _cohort(extra_dev=40)
    sets = build_split(patient_ids=patients, cancers=cancers, source_split=source)["patient_ids"]
    assert not any(patient.startswith("NEWDEV") for patient in sets["test"])
    assert all(patient.startswith(("DEV", "NEWDEV")) for patient in sets["train"] + sets["val"])


def test_assignment_is_deterministic_and_seed_free(tmp_path):
    source = _source_split(tmp_path)
    patients, cancers = _cohort(extra_dev=40, extra_held=10)
    first = build_split(patient_ids=patients, cancers=cancers, source_split=source)
    order = np.random.default_rng(0).permutation(len(patients))
    shuffled = build_split(patient_ids=[patients[i] for i in order],
                           cancers=[cancers[i] for i in order], source_split=source)
    assert first["patient_ids"] == shuffled["patient_ids"]


def test_eligibility_exclusion_is_named_not_counted(tmp_path):
    source = _source_split(tmp_path)
    patients, cancers = _cohort(extra_dev=10)
    dropped = {"NEWDEV000", "NEWDEV001", DEV[0]}
    result = build_split(patient_ids=patients, cancers=cancers, source_split=source,
                         eligible=set(patients) - dropped)
    assert result["eligibility_excluded_count"] == 3
    assert set(result["eligibility_excluded_patients"]) == dropped
    # A SOURCE patient lost to eligibility is a protocol change, so it is reported
    # separately from a merely-new patient that never had an assignment.
    assert result["source_patients_excluded_by_eligibility"] == [DEV[0]]
    assert DEV[0] not in result["patient_ids"]["train"]


def test_patient_from_an_unroled_cancer_is_refused(tmp_path):
    source = _source_split(tmp_path)
    patients, cancers = _cohort()
    patients.append("ORPHAN001"); cancers.append("CANCER_NEVER_SEEN")
    with pytest.raises(ValueError, match="never gave a role"):
        build_split(patient_ids=patients, cancers=cancers, source_split=source)


def test_maximal_split_satisfies_the_runtime_validator(tmp_path):
    source = _source_split(tmp_path)
    patients, cancers = _cohort(extra_dev=22, extra_held=33)
    result = build_split(patient_ids=patients, cancers=cancers, source_split=source)
    destination = tmp_path / "maximal.json"
    destination.write_text(json.dumps(result), encoding="utf-8")
    labels = {patient: name for name, members in result["patient_ids"].items() for patient in members}
    manifest = validate_runtime_split(destination, patients, cancers,
                                      [labels[patient] for patient in patients],
                                      expected_development_cancers=11, expected_heldout_cancers=21)
    assert manifest["counts"]["test"] == len(HELD) + 33
    assert manifest["cancer_counts"] == {"development": 11, "heldout": 21}


# --------------------------------------------------------------------------
# restrict_cohort_to_split
# --------------------------------------------------------------------------

@dataclasses.dataclass
class _Cohort:
    patient_ids: list
    cancers: list
    split: np.ndarray
    rna: np.ndarray
    hallmark_names: list
    hoptimus_store: object = None


def _tiny(tmp_path, declared):
    path = tmp_path / "split.json"
    path.write_text(json.dumps({"patient_ids": declared}), encoding="utf-8")
    return path


def test_restriction_drops_only_unassigned_patients(tmp_path):
    split = _tiny(tmp_path, {"train": ["a", "b"], "val": ["c"], "test": ["d"]})
    # Three column names, deliberately equal in length to nothing else, plus a
    # five-row cohort of which one row is unassigned.
    data = _Cohort(patient_ids=["a", "b", "c", "d", "e"], cancers=["X", "X", "X", "Y", "X"],
                   split=np.asarray(["train", "train", "val", "test", "excluded"]),
                   rna=np.arange(15, dtype=np.float32).reshape(5, 3),
                   hallmark_names=["H1", "H2", "H3"])
    restricted, manifest = restrict_cohort_to_split(data, split, max_dropped_fraction=0.25)
    assert restricted.patient_ids == ["a", "b", "c", "d"]
    assert restricted.rna.shape == (4, 3)
    assert manifest["dropped_patients"] == ["e"]
    assert manifest["retained_count"] == 4
    # Column metadata must survive untouched even though it is a same-typed list.
    assert restricted.hallmark_names == ["H1", "H2", "H3"]


def test_restriction_preserves_row_order(tmp_path):
    split = _tiny(tmp_path, {"train": ["c", "a"], "val": ["e"], "test": ["d"]})
    data = _Cohort(patient_ids=["a", "b", "c", "d", "e"], cancers=list("XXXYX"),
                   split=np.asarray(["train", "excluded", "train", "test", "val"]),
                   rna=np.arange(5, dtype=np.float32).reshape(5, 1), hallmark_names=["H"])
    restricted, _ = restrict_cohort_to_split(data, split, max_dropped_fraction=0.25)
    assert restricted.patient_ids == ["a", "c", "d", "e"]
    assert restricted.rna.ravel().tolist() == [0.0, 2.0, 3.0, 4.0]


def test_restriction_refuses_a_stale_split(tmp_path):
    """The protective direction is not relaxed: a split may be narrower than the
    cohort, never wider. A wider split means the features are gone."""
    split = _tiny(tmp_path, {"train": ["a", "GHOST"], "val": ["c"], "test": ["d"]})
    data = _Cohort(patient_ids=["a", "c", "d"], cancers=list("XXY"),
                   split=np.asarray(["train", "val", "test"]),
                   rna=np.zeros((3, 2), dtype=np.float32), hallmark_names=["H", "I"])
    with pytest.raises(ValueError, match="stale rather than merely narrower"):
        restrict_cohort_to_split(data, split)


def test_restriction_is_a_no_op_when_the_split_already_covers_the_cohort(tmp_path):
    split = _tiny(tmp_path, {"train": ["a"], "val": ["c"], "test": ["d"]})
    data = _Cohort(patient_ids=["a", "c", "d"], cancers=list("XXY"),
                   split=np.asarray(["train", "val", "test"]),
                   rna=np.zeros((3, 2), dtype=np.float32), hallmark_names=["H", "I"])
    restricted, manifest = restrict_cohort_to_split(data, split)
    assert restricted is data
    assert manifest["dropped_count"] == 0


def test_restriction_refuses_a_drop_larger_than_the_ceiling(tmp_path):
    """A stale split used to fail loudly. The flag must not turn that into a quiet
    run on whatever subset the stale file happened to name."""
    split = _tiny(tmp_path, {"train": ["a"], "val": ["c"], "test": ["d"]})
    data = _Cohort(patient_ids=["a", "b", "c", "d", "e"], cancers=list("XXXYX"),
                   split=np.asarray(["train", "excluded", "val", "test", "excluded"]),
                   rna=np.zeros((5, 2), dtype=np.float32), hallmark_names=["H", "I"])
    with pytest.raises(ValueError, match="different cohort, not a narrower one"):
        restrict_cohort_to_split(data, split)


def test_restriction_refuses_to_run_after_targets_are_attached(tmp_path):
    """`dataclasses.replace` would drop the attached matrices without a word."""
    split = _tiny(tmp_path, {"train": ["a"], "val": ["c"], "test": ["d"]})
    data = _Cohort(patient_ids=["a", "b", "c", "d"], cancers=list("XXXY"),
                   split=np.asarray(["train", "excluded", "val", "test"]),
                   rna=np.zeros((4, 2), dtype=np.float32), hallmark_names=["H", "I"])
    data._v2_programmes = np.zeros((4, 3), dtype=np.float32)
    with pytest.raises(ValueError, match="before attaching targets"):
        restrict_cohort_to_split(data, split, max_dropped_fraction=0.5)
