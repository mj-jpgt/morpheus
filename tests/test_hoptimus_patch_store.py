import h5py
import numpy as np
import pandas as pd
import pytest

from morpheus.src.data.hoptimus_patch_store import (
    EXPECTED_HOPTIMUS_DIM,
    HoptimusPatchStore,
    build_dual_split_catalog,
    write_canonical_patch_store,
)


def _metadata():
    return pd.DataFrame([
        {"patient_id": "TCGA-AA-0001", "slide_id": "TCGA-AA-0001-01Z-DX1", "patch_id": "a", "source_id": "s:a", "source_sha256": "a" * 64, "internal_split": "train", "external_split": "valid", "cancer_type": "BRCA"},
        {"patient_id": "TCGA-AA-0001", "slide_id": "TCGA-AA-0001-01Z-DX1", "patch_id": "b", "source_id": "s:b", "source_sha256": "b" * 64, "internal_split": "train", "external_split": "valid", "cancer_type": "BRCA"},
        {"patient_id": "TCGA-AA-0001", "slide_id": "TCGA-AA-0001-01Z-DX2", "patch_id": "c", "source_id": "s:c", "source_sha256": "c" * 64, "internal_split": "train", "external_split": "valid", "cancer_type": "BRCA"},
        {"patient_id": "TCGA-BB-0002", "slide_id": "TCGA-BB-0002-01Z-DX1", "patch_id": "d", "source_id": "s:d", "source_sha256": "d" * 64, "internal_split": "test", "external_split": "test", "cancer_type": "LUAD"},
    ])


def test_store_validates_and_loads_slide_balanced_tokens(tmp_path):
    values = np.arange(4 * EXPECTED_HOPTIMUS_DIM, dtype=np.float32).reshape(4, EXPECTED_HOPTIMUS_DIM)
    store = write_canonical_patch_store(tmp_path, values, _metadata(), {"dataset_revision": "fixture"})
    assert store.validate()["n_patches"] == 4
    tokens, rows = store.load_patient_tokens("TCGA-AA-0001", max_tokens=2, seed=7, slide_balanced=True, required_membership="external:valid")
    assert tokens.shape == (2, EXPECTED_HOPTIMUS_DIM)
    assert set(rows.slide_id) == {"TCGA-AA-0001-01Z-DX1", "TCGA-AA-0001-01Z-DX2"}
    assert np.array_equal(tokens[:, 0], rows.row_idx.to_numpy() * EXPECTED_HOPTIMUS_DIM)
    missing, missing_rows = store.load_patient_tokens("TCGA-BB-0002", required_membership="internal:train")
    assert missing.shape == (0, EXPECTED_HOPTIMUS_DIM)
    assert missing_rows.empty
    assert (tmp_path / "hoptimus_patient_patch_index.parquet").exists()


def test_store_rejects_unversioned_or_wrong_dimension_h5(tmp_path):
    store = write_canonical_patch_store(tmp_path, np.zeros((4, EXPECTED_HOPTIMUS_DIM), np.float32), _metadata(), {})
    with h5py.File(store.embeddings_path, "r+") as handle:
        handle.attrs["store_version"] = "not-supported"
    with pytest.raises(ValueError, match="Unsupported patch store version"):
        store.validate()


def _patch(root, slide, patch, image=b"same-image", label='{"cancer_type":"BRCA"}'):
    directory = root / slide
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"{patch}.jpg").write_bytes(image)
    (directory / f"{patch}.json").write_text(label, encoding="utf-8")


def test_catalog_merges_duplicate_payload_split_memberships(tmp_path):
    internal, external = tmp_path / "internal", tmp_path / "external"
    _patch(internal, "TCGA-AA-0001-01Z-DX1", "0_0_255")
    _patch(external, "TCGA-AA-0001-01Z-DX1", "0_0_255")
    frame = build_dual_split_catalog([(internal, ("internal", "train")), (external, ("external", "valid"))])
    assert len(frame) == 1
    assert frame.loc[0, "patient_id"] == "TCGA-AA-0001"
    assert frame.loc[0, "internal_split"] == "train"
    assert frame.loc[0, "external_split"] == "valid"
