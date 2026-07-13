import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from morpheus.src.data.tumor_state_registry import build_tumor_state_registry, create_tumor_state_splits


def test_tumor_state_registry_uses_disease_labels(tmp_path: Path):
    root = tmp_path
    processed = root / "processed"
    processed.mkdir()
    meta_dir = root / "meta-intersurv" / "data" / "embeddings" / "patient_master"
    meta_dir.mkdir(parents=True)
    pd.DataFrame({"patient_id": ["TCGA-AA-0001", "TCGA-AA-0002"], "cancer_type": ["BRCA", "LUAD"], "time": [10.0, 20.0], "event": [1, 0]}).to_parquet(meta_dir / "survival_labels.parquet", index=False)
    pd.DataFrame({"patient_id": ["TCGA-AA-0001", "TCGA-AA-0002"], "age": [1.0, 2.0]}).to_parquet(meta_dir / "tabular_features.parquet", index=False)
    rna_dir = root / "rna"
    rna_dir.mkdir()
    pd.DataFrame({"patient_id": ["TCGA-AA-0001", "TCGA-AA-0002"], **{f"bulkformer_{i:03d}": [0.0, 1.0] for i in range(512)}}).to_parquet(rna_dir / "bulk.parquet", index=False)
    pd.DataFrame({"patient_id": ["TCGA-AA-0001", "TCGA-AA-0002"], "sample_id": ["TCGA-AA-0001-01", "TCGA-AA-0002-01"], "cancer_type": ["BRCA", "LUAD"]}).to_parquet(rna_dir / "tcga_bulkformer_embedding_metadata.parquet", index=False)
    hallmark = root / "hallmark.parquet"
    pd.DataFrame({"patient_id": ["TCGA-AA-0001"], "HALLMARK_X": [0.5]}).to_parquet(hallmark, index=False)
    wsi_dir = root / "wsi"
    wsi_dir.mkdir()
    with h5py.File(wsi_dir / "tcga_ut_hoptimus0_patient_embeddings.h5", "w") as handle:
        handle.create_dataset("patient_ids", data=np.asarray(["TCGA-AA-0001", "TCGA-AA-0002"], dtype="S"))
        handle.create_dataset("embeddings", data=np.zeros((2, 1536), dtype=np.float32))
    cfg = {
        "project_root": str(root),
        "seed": 1,
        "paths": {
            "processed_dir": "processed",
            "outputs_dir": "outputs",
            "rna_bulkformer_embeddings": "rna/bulk.parquet",
            "wsi_standard_dir": "wsi",
            "hallmark_scores": "hallmark.parquet",
            "snv_processed": "missing_snv.parquet",
            "cnv_processed": "missing_cnv.parquet",
        },
        "splits": {"val_fraction": 0.2, "test_fraction": 0.2},
    }
    cfg_path = root / "config.json"
    cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
    outputs = build_tumor_state_registry(cfg_path)
    splits = create_tumor_state_splits(cfg_path)
    registry = pd.read_parquet(outputs["registry"])
    assert set(registry["cancer_type"]) == {"BRCA", "LUAD"}
    assert splits["stratified"].exists()

