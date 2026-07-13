import numpy as np
import pandas as pd
import pytest

from morpheus.src.eval.baseline_registry import baseline_registry_frame
from morpheus.src.eval.tcga_ut_patch_benchmarks import (
    evaluate_protocol_benchmark,
    load_aligned_protocol_archive,
    protocol_splits_from_npz,
)


def _archive(path, **extra):
    values = {
        "patient_ids": np.asarray(["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003"]),
        "wsi": np.eye(3, dtype=np.float32),
        "rna": np.eye(3, dtype=np.float32),
    }
    values.update(extra)
    np.savez(path, **values)


def test_protocol_archive_uses_explicit_internal_membership(tmp_path):
    path = tmp_path / "embeddings.npz"
    _archive(path, internal_split=np.asarray(["train", "val", "test"]), external_split=np.asarray(["test", "train", "val"]))
    loaded = load_aligned_protocol_archive(path, "internal")
    assert loaded["split"].tolist() == ["train", "val", "test"]


def test_protocol_archive_rejects_ambiguous_legacy_split(tmp_path):
    path = tmp_path / "embeddings.npz"
    _archive(path, split=np.asarray(["train", "val", "test"]))
    with np.load(path) as data:
        with pytest.raises(KeyError, match="ambiguous"):
            protocol_splits_from_npz(data, "internal")


def test_protocol_archive_rejects_wrong_declared_protocol(tmp_path):
    path = tmp_path / "embeddings.npz"
    _archive(path, split=np.asarray(["train", "val", "test"]), split_protocol=np.asarray("external"))
    with np.load(path) as data:
        with pytest.raises(ValueError, match="not 'internal'"):
            protocol_splits_from_npz(data, "internal")


def test_protocol_archive_rejects_duplicate_patients(tmp_path):
    path = tmp_path / "embeddings.npz"
    _archive(path, patient_ids=np.asarray(["TCGA-AA-0001", "TCGA-AA-0001", "TCGA-AA-0003"]), internal_split=np.asarray(["train", "val", "test"]))
    with pytest.raises(ValueError, match="one row per patient"):
        load_aligned_protocol_archive(path, "internal")


def test_baseline_registry_contains_required_same_data_controls():
    registry = baseline_registry_frame()
    assert isinstance(registry, pd.DataFrame)
    assert {"H-Optimus-0 patch mean", "ABMIL on H-Optimus-0", "CLAM on H-Optimus-0", "Token-aware BioQueryFormer"}.issubset(set(registry["name"]))


def test_protocol_evaluation_writes_train_only_prompting_result(tmp_path):
    ids = np.asarray([f"TCGA-AA-{i:04d}" for i in range(9)])
    archive = tmp_path / "aligned.npz"
    _archive(
        archive,
        patient_ids=ids,
        wsi=np.eye(9, dtype=np.float32),
        rna=np.eye(9, dtype=np.float32),
        internal_split=np.asarray(["train"] * 5 + ["val"] * 2 + ["test"] * 2),
    )
    targets = tmp_path / "targets.csv"
    pd.DataFrame({"patient_id": ids, "hallmark": np.arange(9, dtype=float)}).to_csv(targets, index=False)
    result = evaluate_protocol_benchmark(archive, targets, "internal", tmp_path / "result")
    payload = __import__("json").loads(result.read_text(encoding="utf-8"))
    assert payload["splits"]["test"]["molecular_prompting"]["n_train_reference"] == 5
    assert payload["splits"]["test"]["molecular_prompting"]["n_eval"] == 2
