"""Local BulkFormer artifact and embedding-store utilities."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


EXPECTED_BULKFORMER_DIM = 512


def selected_artifact_paths(manifest_path: str | Path) -> dict[str, Path]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    return {category: Path(record["local_path"]) for category, record in manifest.get("selected", {}).items()}


def require_bulkformer_artifacts(manifest_path: str | Path) -> dict[str, Path]:
    paths = selected_artifact_paths(manifest_path)
    required = ("checkpoint", "gene_info", "gene_length", "tcga_h5ad", "graph", "graph_weight", "gene_embedding", "interested_gene_list")
    missing = [name for name in required if name not in paths or not paths[name].exists()]
    if missing:
        raise FileNotFoundError(f"BulkFormer artifacts are incomplete: {', '.join(missing)}")
    return paths


def write_embedding_store(patient_ids: Iterable[str], embeddings: np.ndarray, output_path: str | Path) -> Path:
    matrix = np.asarray(embeddings, dtype=np.float32)
    if matrix.ndim != 2:
        raise ValueError(f"BulkFormer embeddings must be 2D, got {matrix.shape}")
    if matrix.shape[1] != EXPECTED_BULKFORMER_DIM:
        raise ValueError(f"BulkFormer embeddings must be {EXPECTED_BULKFORMER_DIM}D, got {matrix.shape[1]}")
    ids = list(patient_ids)
    if len(ids) != matrix.shape[0]:
        raise ValueError(f"patient_id count {len(ids)} does not match embedding rows {matrix.shape[0]}")
    cols = [f"bulkformer_{i:03d}" for i in range(matrix.shape[1])]
    frame = pd.DataFrame(matrix, columns=cols)
    frame.insert(0, "patient_id", ids)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(output, index=False)
    return output


def load_embedding_store(path: str | Path) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    if "patient_id" not in frame:
        raise ValueError(f"BulkFormer embedding store missing patient_id: {path}")
    numeric = frame.select_dtypes(include=["number"]).columns
    if len(numeric) != EXPECTED_BULKFORMER_DIM:
        raise ValueError(f"BulkFormer embedding store must have {EXPECTED_BULKFORMER_DIM} numeric columns, got {len(numeric)}")
    return frame
