"""Standardize existing H-Optimus-0 embeddings into root data stores."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from morpheus.src.utils.ids import normalize_patient_id, parse_tcga_barcode
from morpheus.src.utils.provenance import write_json


def _write_h5(path: Path, embeddings: np.ndarray, ids: list[str], id_name: str = "ids") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as handle:
        handle.create_dataset("embeddings", data=np.asarray(embeddings, dtype=np.float32), compression="gzip")
        handle.create_dataset(id_name, data=np.asarray(ids, dtype=object), dtype=string_dtype)


def load_hoptimus_patient_splits(source_dir: str | Path) -> tuple[np.ndarray, pd.DataFrame]:
    source = Path(source_dir)
    arrays = []
    metas = []
    for split in ("train", "valid", "test"):
        emb_path = source / f"{split}.emb.npy"
        meta_path = source / f"{split}.meta.parquet"
        if not emb_path.exists() or not meta_path.exists():
            continue
        arr = np.load(emb_path).astype(np.float32)
        meta = pd.read_parquet(meta_path).copy()
        if len(meta) != arr.shape[0]:
            raise ValueError(f"{split} metadata rows {len(meta)} != embeddings rows {arr.shape[0]}")
        meta["source_split"] = split
        meta["row_idx"] = np.arange(len(meta)) + sum(len(m) for m in metas)
        arrays.append(arr)
        metas.append(meta)
    if not arrays:
        raise FileNotFoundError(f"No H-Optimus-0 split embeddings found in {source}")
    matrix = np.concatenate(arrays, axis=0)
    meta = pd.concat(metas, ignore_index=True)
    patient_col = next((c for c in meta.columns if "patient" in str(c).lower()), None)
    if not patient_col:
        raise ValueError("Could not find patient column in H-Optimus metadata")
    meta["patient_id"] = meta[patient_col].map(normalize_patient_id)
    meta["embedding_dim"] = int(matrix.shape[1])
    meta["encoder_name"] = "H-Optimus-0"
    meta["embedding_level"] = "patient"
    return matrix, meta


def build_standard_hoptimus_store(
    source_dir: str | Path = "meta-intersurv/data/tcga_ut/embeddings_hoptimus0",
    patch_dir: str | Path = "meta-intersurv/data/wsi/wsi_embeddings",
    output_dir: str | Path = "data/processed/wsi",
) -> dict[str, str | int]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    matrix, meta = load_hoptimus_patient_splits(source_dir)
    patch_rows = []
    patch_root = Path(patch_dir)
    if patch_root.exists():
        for fp in sorted(patch_root.glob("*.npz")):
            barcode = parse_tcga_barcode(fp.name)
            patch_rows.append({"patient_id": barcode.patient_id, "slide_id": fp.stem, "source_path": str(fp)})
    patch_meta = pd.DataFrame(patch_rows)
    patch_meta_path = output / "tcga_ut_hoptimus0_patch_metadata.parquet"
    patch_meta.to_parquet(patch_meta_path, index=False)

    patient_rows = []
    patient_vectors = []
    for patient_id, group in meta.dropna(subset=["patient_id"]).groupby("patient_id", sort=True):
        idx = group.index.to_numpy()
        patient_vectors.append(matrix[idx].mean(axis=0))
        patient_rows.append(patient_id)
    patient_matrix = np.vstack(patient_vectors).astype(np.float32)
    patient_path = output / "tcga_ut_hoptimus0_patient_embeddings.h5"
    _write_h5(patient_path, patient_matrix, patient_rows, "patient_ids")

    slide_path = output / "tcga_ut_hoptimus0_slide_embeddings.h5"
    _write_h5(slide_path, matrix, meta["patient_id"].fillna("").astype(str).tolist(), "patient_ids")
    slide_meta_path = output / "tcga_ut_hoptimus0_patch_embeddings_metadata.parquet"
    meta.to_parquet(slide_meta_path, index=False)

    manifest = {
        "source_dir": str(source_dir),
        "patch_dir": str(patch_dir),
        "patient_embeddings": str(patient_path),
        "slide_embeddings": str(slide_path),
        "patch_metadata": str(patch_meta_path),
        "split_metadata": str(slide_meta_path),
        "n_rows": int(matrix.shape[0]),
        "n_patients": int(len(patient_rows)),
        "embedding_dim": int(matrix.shape[1]),
        "note": "Patch-level tensor h5 is intentionally not materialized by default; source npz paths are preserved in metadata.",
    }
    manifest_path = output / "tcga_ut_hoptimus0_manifest.json"
    write_json(manifest_path, manifest)
    return {**manifest, "manifest": str(manifest_path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", default="meta-intersurv/data/tcga_ut/embeddings_hoptimus0")
    parser.add_argument("--patch-dir", default="meta-intersurv/data/wsi/wsi_embeddings")
    parser.add_argument("--output-dir", default="data/processed/wsi")
    args = parser.parse_args()
    print(build_standard_hoptimus_store(args.source_dir, args.patch_dir, args.output_dir)["manifest"])


if __name__ == "__main__":
    main()
