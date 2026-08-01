"""Open PLIP semantic-teacher preparation and strict cache validation.

PLIP is optional.  This module is intentionally isolated from core V2 so a
strict WSI--RNA run never needs an external vision-language dependency.
"""
from __future__ import annotations

import argparse
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


PLIP_MODEL_ID = "vinid/plip"


def _normalise(value: np.ndarray) -> np.ndarray:
    return value / np.maximum(np.linalg.norm(value, axis=-1, keepdims=True), 1e-12)


def canonical_patient_source_digests(metadata: pd.DataFrame) -> dict[str, str]:
    required = {"patient_id", "source_sha256"}
    if not required.issubset(metadata.columns):
        raise ValueError("canonical metadata lacks source checksum provenance")
    return {str(patient): "|".join(sorted(group.source_sha256.astype(str).unique()))
            for patient, group in metadata.groupby("patient_id", sort=False)}


def load_plip_teacher_cache(path: str | Path, patient_ids: list[str], source_digests: dict[str, str] | None = None) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Align a cache by patient ID and expose absent patients explicitly."""
    with np.load(path, allow_pickle=False) as raw:
        required = {"patient_ids", "embeddings", "model_id", "source_sha256"}
        if not required.issubset(raw.files):
            raise ValueError(f"PLIP cache missing {sorted(required - set(raw.files))}")
        ids = raw["patient_ids"].astype(str)
        if len(np.unique(ids)) != len(ids):
            raise ValueError("PLIP cache contains duplicate patient IDs")
        vectors = raw["embeddings"].astype(np.float32)
        if vectors.ndim != 2 or len(vectors) != len(ids) or not np.isfinite(vectors).all():
            raise ValueError("PLIP cache embeddings are invalid")
        model_id = str(raw["model_id"].item())
        if model_id != PLIP_MODEL_ID:
            raise ValueError(f"expected {PLIP_MODEL_ID!r}, found {model_id!r}")
        hashes = raw["source_sha256"].astype(str)
    lookup = {patient: row for row, patient in enumerate(ids)}
    rows = np.asarray([lookup.get(str(patient), -1) for patient in patient_ids])
    present = rows >= 0
    if source_digests is not None:
        mismatched = [str(patient) for patient, row in zip(patient_ids, rows) if row >= 0 and hashes[row] != source_digests.get(str(patient), "")]
        if mismatched:
            raise ValueError(f"PLIP cache checksum mismatch for {len(mismatched)} active patients")
    aligned = np.zeros((len(patient_ids), vectors.shape[1]), dtype=np.float32)
    aligned[present] = _normalise(vectors[rows[present]])
    return aligned, present, {"model_id": model_id, "n_cached_patients": int(len(ids)), "n_aligned_patients": int(present.sum()), "source_checksum_count": int(len(np.unique(hashes)),), "checksum_verified": source_digests is not None}


def _open_image(source: str, member: str):
    from PIL import Image
    source_path = Path(source)
    if source_path.is_dir():
        return Image.open(source_path / member).convert("RGB")
    with tarfile.open(source_path, "r:*") as archive:
        handle = archive.extractfile(member)
        if handle is None:
            raise FileNotFoundError(f"missing {source_path}:{member}")
        image = Image.open(handle).convert("RGB")
        return image.copy()


def build_patient_teacher(metadata_path: str | Path, output: str | Path, batch_size: int = 128, device: str = "cuda") -> Path:
    """Embed exact source JPEGs and aggregate normalized PLIP tokens by patient.

    Metadata must originate from the canonical patch catalog and therefore
    include source path/member/checksum provenance rather than inferred files.
    """
    required = {"patient_id", "source_path", "image_member", "source_sha256"}
    table = pd.read_parquet(metadata_path)
    if not required.issubset(table.columns):
        raise ValueError(f"metadata missing {sorted(required - set(table.columns))}")
    if table.source_sha256.duplicated().any():
        table = table.drop_duplicates("source_sha256").copy()
    try:
        import torch
        from transformers import CLIPModel, CLIPProcessor
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("PLIP teacher extraction requires transformers, torch, and pillow") from exc
    model = CLIPModel.from_pretrained(PLIP_MODEL_ID).to(device).eval()
    processor = CLIPProcessor.from_pretrained(PLIP_MODEL_ID)
    vectors: list[np.ndarray] = []
    rows: list[pd.DataFrame] = []
    for start in range(0, len(table), batch_size):
        subset = table.iloc[start : start + batch_size]
        images = [_open_image(str(row.source_path), str(row.image_member)) for row in subset.itertuples()]
        encoded = processor(images=images, return_tensors="pt")
        encoded = {key: value.to(device) for key, value in encoded.items()}
        with torch.inference_mode():
            value = model.get_image_features(**encoded).float().cpu().numpy()
        vectors.append(_normalise(value.astype(np.float32))); rows.append(subset)
    values = np.vstack(vectors); records = pd.concat(rows, ignore_index=True)
    patient_ids, patient_vectors, patient_hashes = [], [], []
    for patient, group in records.assign(_row=np.arange(len(records))).groupby("patient_id", sort=True):
        patient_ids.append(str(patient)); patient_vectors.append(_normalise(values[group._row.to_numpy()].mean(0, keepdims=True))[0])
        patient_hashes.append("|".join(sorted(group.source_sha256.astype(str))))
    target = Path(output); target.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(target, patient_ids=np.asarray(patient_ids), embeddings=np.asarray(patient_vectors, dtype=np.float32),
                        source_sha256=np.asarray(patient_hashes), model_id=np.asarray(PLIP_MODEL_ID), aggregation=np.asarray("mean_normalized_patch_embeddings"))
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description="Build checksum-audited PLIP patient semantic targets")
    parser.add_argument("--metadata", required=True); parser.add_argument("--output", required=True)
    parser.add_argument("--batch-size", type=int, default=128); parser.add_argument("--device", default="cuda")
    print(build_patient_teacher(parser.parse_args().metadata, parser.parse_args().output, parser.parse_args().batch_size, parser.parse_args().device))


if __name__ == "__main__":
    main()
