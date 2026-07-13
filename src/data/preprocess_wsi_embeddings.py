"""Normalize existing WSI embeddings into a v1 feature store."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.src.utils.config import load_config
from morpheus.src.utils.ids import normalize_patient_id, parse_tcga_barcode
from morpheus.src.utils.provenance import base_manifest, write_json


def _load_hoptimus(path: Path) -> tuple[np.ndarray, pd.DataFrame] | None:
    train = path / "train.emb.npy"
    if not train.exists():
        return None
    arrays = []
    metas = []
    for split in ("train", "valid", "test"):
        emb = path / f"{split}.emb.npy"
        meta = path / f"{split}.meta.parquet"
        if not emb.exists() or not meta.exists():
            continue
        arr = np.load(emb).astype(np.float32)
        df = pd.read_parquet(meta).copy()
        if len(df) != arr.shape[0]:
            raise ValueError(f"Metadata/embedding row mismatch for {split}: {len(df)} vs {arr.shape[0]}")
        df["source_split"] = split
        arrays.append(arr)
        metas.append(df)
    if not arrays:
        return None
    matrix = np.concatenate(arrays, axis=0)
    meta = pd.concat(metas, ignore_index=True)
    patient_col = next((c for c in meta.columns if "patient" in str(c).lower()), None)
    if patient_col is None:
        raise ValueError(f"Could not find patient column in {path}")
    meta["patient_id"] = meta[patient_col].map(normalize_patient_id)
    meta["encoder_name"] = "H-Optimus0"
    meta["embedding_level"] = "patient"
    meta["coordinate_available"] = False
    return matrix, meta


def _load_patch_npz(path: Path, max_patches_per_slide: int = 1500) -> tuple[np.ndarray, pd.DataFrame]:
    rows = []
    vectors = []
    for fp in sorted(path.glob("*.npz")):
        data = np.load(fp, allow_pickle=False)
        feats = np.asarray(data["feats"], dtype=np.float32)
        if feats.ndim != 2:
            continue
        if len(feats) > max_patches_per_slide:
            feats = feats[:max_patches_per_slide]
        vector = feats.mean(axis=0)
        barcode = parse_tcga_barcode(fp.name)
        vectors.append(vector)
        rows.append(
            {
                "patient_id": barcode.patient_id,
                "slide_id": fp.stem,
                "source_path": str(fp),
                "encoder_name": "precomputed_wsi_npz",
                "embedding_level": "slide_mean_from_patches",
                "coordinate_available": "coords" in data.files,
                "n_patches": int(feats.shape[0]),
                "feature_dim": int(feats.shape[1]),
            }
        )
    if not vectors:
        return np.empty((0, 0), dtype=np.float32), pd.DataFrame()
    return np.vstack(vectors).astype(np.float32), pd.DataFrame(rows)


def build_wsi_feature_store(config_path: str = "morpheus/configs/v1.json") -> Path:
    cfg = load_config(config_path)
    out_dir = cfg.path("processed_dir") / "features"
    out_dir.mkdir(parents=True, exist_ok=True)
    hoptimus = _load_hoptimus(cfg.path("wsi_hoptimus_dir"))
    if hoptimus is not None:
        matrix, meta = hoptimus
        source = "hoptimus_patient_embeddings"
    else:
        matrix, meta = _load_patch_npz(cfg.path("wsi_patch_dir"))
        source = "patch_npz_mean_pool"
    if matrix.size == 0:
        raise FileNotFoundError("No WSI embeddings found in configured H-Optimus0 or patch directories")
    feature_path = out_dir / "wsi_features.npz"
    meta_path = out_dir / "wsi_metadata.parquet"
    np.savez_compressed(feature_path, embeddings=matrix)
    meta.to_parquet(meta_path, index=False)
    manifest = base_manifest(cfg.project_root, cfg.config_path, int(cfg.raw.get("seed", 42)))
    manifest.update(
        {
            "source": source,
            "feature_path": str(feature_path),
            "metadata_path": str(meta_path),
            "shape": list(matrix.shape),
            "patients": int(meta["patient_id"].nunique()),
        }
    )
    write_json(out_dir / "wsi_features.manifest.json", manifest)
    return feature_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    args = parser.parse_args()
    print(build_wsi_feature_store(args.config))


if __name__ == "__main__":
    main()
