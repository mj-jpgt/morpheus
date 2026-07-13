"""Registry and loading utilities for WSI patch-bag feature stores."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from numpy.lib import format as npy_format

from morpheus.src.utils.provenance import write_json


REQUIRED_PATCH_KEYS = ("feats", "coords", "patient_id", "slide_id")


@dataclass(frozen=True)
class PatchBag:
    patient_id: str
    slide_id: str
    feats: np.ndarray
    coords: np.ndarray | None
    edge_index: np.ndarray | None


def _scalar_str(value) -> str:
    arr = np.asarray(value)
    return str(arr.item() if arr.shape == () else arr)


def _npz_member_shape(path: Path, key: str) -> tuple[int, ...]:
    with zipfile.ZipFile(path) as zf:
        with zf.open(f"{key}.npy") as handle:
            version = npy_format.read_magic(handle)
            if version == (1, 0):
                shape, _, _ = npy_format.read_array_header_1_0(handle)
            elif version == (2, 0):
                shape, _, _ = npy_format.read_array_header_2_0(handle)
            else:
                shape, _, _ = npy_format.read_array_header_2_0(handle)
    return tuple(int(x) for x in shape)


def inspect_patch_npz(path: str | Path) -> dict:
    fp = Path(path)
    with zipfile.ZipFile(fp) as zf:
        keys = {Path(name).stem for name in zf.namelist() if name.endswith(".npy")}
    missing = [key for key in ("feats", "coords") if key not in keys]
    if missing:
        raise ValueError(f"{fp} is missing required patch-bag keys: {missing}")
    feat_shape = _npz_member_shape(fp, "feats")
    coord_shape = _npz_member_shape(fp, "coords") if "coords" in keys else None
    if len(feat_shape) != 2:
        raise ValueError(f"{fp} feats must be 2D, got {feat_shape}")
    if coord_shape is not None and (len(coord_shape) != 2 or coord_shape[0] != feat_shape[0]):
        raise ValueError(f"{fp} coords shape {coord_shape} does not match feats {feat_shape}")
    slide_id = fp.stem
    patient_id = slide_id[:12] if slide_id.startswith("TCGA-") and len(slide_id) >= 12 else None
    if patient_id is None:
        data = np.load(fp, allow_pickle=True)
        patient_id = _scalar_str(data["patient_id"]) if "patient_id" in data.files else slide_id
        slide_id = _scalar_str(data["slide_id"]) if "slide_id" in data.files else slide_id
    return {
        "patient_id": patient_id,
        "slide_id": slide_id,
        "source_path": str(fp),
        "feature_source": "existing_wsi_patchbag_2048",
        "encoder_name": "unknown_precomputed_patch_encoder",
        "embedding_level": "patch",
        "n_tokens": int(feat_shape[0]),
        "feature_dim": int(feat_shape[1]),
        "has_coords": bool("coords" in keys),
        "has_edge_index_8": bool("edge_index_8" in keys),
        "cancer_type": None,
    }


def inspect_patch_filename(path: str | Path) -> dict:
    fp = Path(path)
    slide_id = fp.stem
    patient_id = slide_id[:12] if slide_id.startswith("TCGA-") and len(slide_id) >= 12 else slide_id
    return {
        "patient_id": patient_id,
        "slide_id": slide_id,
        "source_path": str(fp),
        "feature_source": "existing_wsi_patchbag_2048",
        "encoder_name": "unknown_precomputed_patch_encoder",
        "embedding_level": "patch",
        "n_tokens": None,
        "feature_dim": 2048,
        "has_coords": None,
        "has_edge_index_8": None,
        "cancer_type": None,
        "registry_mode": "filename_fast_unverified",
    }


def build_patch_bag_registry(
    patch_dir: str | Path = "meta-intersurv/data/wsi/wsi_embeddings",
    output_path: str | Path = "data/processed/wsi/tcga_ut_patch_bag_registry.parquet",
    limit: int | None = None,
    offset: int = 0,
    progress_every: int = 250,
    fast_filename_registry: bool = False,
) -> Path:
    patch_root = Path(patch_dir)
    paths = sorted(patch_root.glob("*.npz"))
    if offset > 0:
        paths = paths[offset:]
    if limit is not None and limit > 0:
        paths = paths[:limit]
    rows = []
    for i, path in enumerate(paths, start=1):
        rows.append(inspect_patch_filename(path) if fast_filename_registry else inspect_patch_npz(path))
        if progress_every > 0 and i % progress_every == 0:
            print(f"inspected {i}/{len(paths)} patch bags", flush=True)
    if not rows:
        raise FileNotFoundError(f"No patch-bag npz files found in {patch_root}")
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_parquet(out, index=False)
    manifest = {
        "patch_dir": str(patch_root),
        "registry": str(out),
        "n_slides": int(len(frame)),
        "n_patients": int(frame["patient_id"].nunique()),
        "limit": int(limit) if limit is not None else None,
        "offset": int(offset),
        "feature_dims": sorted(int(x) for x in frame["feature_dim"].unique()),
        "feature_source": "existing_wsi_patchbag_2048",
        "fast_filename_registry": bool(fast_filename_registry),
        "note": "These are local precomputed patch bags. They are intentionally not labeled H-Optimus-0 unless provenance is confirmed. Fast filename registries do not verify per-file shapes.",
    }
    write_json(out.with_suffix(".manifest.json"), manifest)
    return out


def load_patch_bag(path: str | Path, max_tokens: int = 512, seed: int = 42) -> PatchBag:
    fp = Path(path)
    data = np.load(fp, allow_pickle=True)
    feats = np.asarray(data["feats"], dtype=np.float32)
    coords = np.asarray(data["coords"], dtype=np.float32) if "coords" in data.files else None
    edge_index = np.asarray(data["edge_index_8"], dtype=np.int64) if "edge_index_8" in data.files else None
    if max_tokens > 0 and feats.shape[0] > max_tokens:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(feats.shape[0], size=max_tokens, replace=False))
        feats = feats[idx]
        coords = coords[idx] if coords is not None else None
        edge_index = None
    return PatchBag(
        patient_id=_scalar_str(data["patient_id"]),
        slide_id=_scalar_str(data["slide_id"]),
        feats=feats,
        coords=coords,
        edge_index=edge_index,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--patch-dir", default="meta-intersurv/data/wsi/wsi_embeddings")
    parser.add_argument("--output-path", default="data/processed/wsi/tcga_ut_patch_bag_registry.parquet")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=250)
    parser.add_argument("--fast-filename-registry", action="store_true")
    args = parser.parse_args()
    print(build_patch_bag_registry(args.patch_dir, args.output_path, args.limit, args.offset, args.progress_every, args.fast_filename_registry))


if __name__ == "__main__":
    main()
