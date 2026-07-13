import numpy as np

from morpheus.src.data.wsi_patch_bags import inspect_patch_npz, load_patch_bag


def test_patch_bag_inspection_and_sampling(tmp_path):
    path = tmp_path / "slide.npz"
    np.savez(
        path,
        feats=np.arange(40, dtype=np.float32).reshape(10, 4),
        coords=np.arange(20, dtype=np.float32).reshape(10, 2),
        patient_id=np.asarray("TCGA-XX-0001"),
        slide_id=np.asarray("slide-1"),
        cancer_type=np.asarray("BRCA"),
        edge_index_8=np.zeros((2, 4), dtype=np.int64),
    )
    row = inspect_patch_npz(path)
    assert row["patient_id"] == "TCGA-XX-0001"
    assert row["n_tokens"] == 10
    assert row["feature_dim"] == 4
    assert row["has_edge_index_8"]

    bag = load_patch_bag(path, max_tokens=5, seed=1)
    assert bag.feats.shape == (5, 4)
    assert bag.coords.shape == (5, 2)
    assert bag.edge_index is None
