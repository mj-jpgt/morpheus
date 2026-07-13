"""Strict, versioned H-Optimus-0 patch-store primitives for TCGA-UT."""
from __future__ import annotations

import hashlib
import json
import os
import tarfile
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, Iterator, Sequence

import h5py
import numpy as np
import pandas as pd

EXPECTED_HOPTIMUS_DIM = 1536
STORE_FORMAT = "morpheus.hoptimus_patch_store"
STORE_VERSION = "1.0"
EMBEDDINGS_DATASET = "embeddings"
VALID_PARTITIONS = frozenset(("train", "valid", "test"))
REQUIRED_METADATA_COLUMNS = frozenset((
    "row_idx", "patient_id", "slide_id", "patch_id", "source_id",
    "source_sha256", "internal_split", "external_split", "cancer_type",
))


def normalize_tcga_patient_id(value: str) -> str:
    """Return the canonical 12-character TCGA patient barcode or raise."""
    parts = str(value).strip().upper().split("-")
    if len(parts) < 3 or parts[0] != "TCGA" or len(parts[1]) != 2 or len(parts[2]) != 4:
        raise ValueError(f"Not a TCGA patient/slide identifier: {value!r}")
    return "-".join(parts[:3])


def split_membership_column(scheme: str) -> str:
    scheme = str(scheme).strip().lower()
    if scheme not in {"internal", "external"}:
        raise ValueError("Split scheme must be 'internal' or 'external'")
    return f"{scheme}_split"


def parse_membership(value: str) -> tuple[str, str]:
    try:
        scheme, partition = value.split(":", 1)
    except ValueError as exc:
        raise ValueError("Membership must be internal:train or external:test") from exc
    partition = partition.strip().lower()
    if partition not in VALID_PARTITIONS:
        raise ValueError(f"Unknown split partition {partition!r}")
    return scheme.strip().lower(), partition


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _cancer_type(payload: bytes) -> str | None:
    try:
        value = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("Patch sidecar JSON is invalid") from exc
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        for key in ("cancer_type", "label", "class", "tumor_type"):
            if value.get(key) is not None:
                return str(value[key])
    return None


def _identity(member: str) -> tuple[str, str, str]:
    path = Path(member.replace("\\", "/"))
    if path.suffix.lower() not in {".jpg", ".jpeg"} or len(path.parts) < 2:
        raise ValueError(f"Expected <slide>/<patch>.jpg, got {member!r}")
    slide_id = path.parts[-2]
    return normalize_tcga_patient_id(slide_id), slide_id, path.stem


@dataclass(frozen=True)
class PatchSource:
    source_path: Path
    image_member: str
    json_member: str
    patient_id: str
    slide_id: str
    patch_id: str
    cancer_type: str | None
    source_sha256: str
    internal_split: str | None = None
    external_split: str | None = None

    @property
    def source_id(self) -> str:
        return f"{self.source_path.name}:{self.image_member}"


def _source_kwargs(membership: tuple[str, str] | None) -> dict[str, str | None]:
    values: dict[str, str | None] = {"internal_split": None, "external_split": None}
    if membership:
        values[split_membership_column(membership[0])] = membership[1]
    return values


def _iter_directory(root: Path, membership: tuple[str, str] | None) -> Iterator[PatchSource]:
    images = sorted(list(root.rglob("*.jpg")) + list(root.rglob("*.jpeg")))
    for image in images:
        sidecar = image.with_suffix(".json")
        if not sidecar.exists():
            raise ValueError(f"Missing JSON sidecar for {image}")
        relative = image.relative_to(root).as_posix()
        patient, slide, patch = _identity(relative)
        yield PatchSource(root, relative, sidecar.relative_to(root).as_posix(), patient, slide, patch,
                          _cancer_type(sidecar.read_bytes()), _sha256(image.read_bytes()), **_source_kwargs(membership))


def _iter_tar(shard: Path, membership: tuple[str, str] | None) -> Iterator[PatchSource]:
    with tarfile.open(shard, "r:*") as archive:
        members = {m.name: m for m in archive.getmembers() if m.isfile()}
        for name in sorted(k for k in members if Path(k).suffix.lower() in {".jpg", ".jpeg"}):
            json_name = str(Path(name).with_suffix(".json")).replace("\\", "/")
            if json_name not in members:
                raise ValueError(f"Missing JSON sidecar for {shard}:{name}")
            image_f, json_f = archive.extractfile(members[name]), archive.extractfile(members[json_name])
            if image_f is None or json_f is None:
                raise ValueError(f"Unreadable archive member {shard}:{name}")
            image_payload, json_payload = image_f.read(), json_f.read()
            patient, slide, patch = _identity(name)
            yield PatchSource(shard, name, json_name, patient, slide, patch, _cancer_type(json_payload),
                              _sha256(image_payload), **_source_kwargs(membership))


def iter_patch_sources(path: str | Path, membership: tuple[str, str] | None = None) -> Iterator[PatchSource]:
    """Enumerate validated paired JPEG/JSON records from a tar or extracted shard."""
    source = Path(path)
    if source.is_dir():
        yield from _iter_directory(source, membership)
    elif source.is_file():
        # TCGA-UT shard names can carry a numeric suffix (for example
        # ``.tar110``); tarfile inspects content rather than relying on suffix.
        yield from _iter_tar(source, membership)
    else:
        raise FileNotFoundError(f"TCGA-UT source must be a directory or tar: {source}")


def build_dual_split_catalog(sources: Sequence[tuple[str | Path, tuple[str, str]]], output_path: str | Path | None = None) -> pd.DataFrame:
    """Deduplicate byte-identical tiles and merge their internal/external memberships."""
    records: dict[str, dict[str, Any]] = {}
    for path, membership in sources:
        for source in iter_patch_sources(path, membership):
            record = records.get(source.source_sha256)
            if record is None:
                record = {"patient_id": source.patient_id, "slide_id": source.slide_id, "patch_id": source.patch_id,
                          "source_id": source.source_id, "source_path": str(source.source_path),
                          "image_member": source.image_member, "json_member": source.json_member,
                          "source_sha256": source.source_sha256, "cancer_type": source.cancer_type,
                          "internal_split": source.internal_split, "external_split": source.external_split}
                records[source.source_sha256] = record
            else:
                if any(record[k] != getattr(source, k) for k in ("patient_id", "slide_id", "patch_id")):
                    raise ValueError(f"Checksum identity conflict for {source.source_id}")
                for col in ("internal_split", "external_split"):
                    value = getattr(source, col)
                    if value and record[col] and value != record[col]:
                        raise ValueError(f"Patch occurs in incompatible {col} partitions: {source.source_id}")
                    record[col] = record[col] or value
    if not records:
        raise ValueError("No paired TCGA-UT images found")
    frame = pd.DataFrame(records.values()).sort_values(["patient_id", "slide_id", "patch_id", "source_sha256"]).reset_index(drop=True)
    frame.insert(0, "row_idx", np.arange(len(frame), dtype=np.int64))
    if output_path:
        target = Path(output_path); target.parent.mkdir(parents=True, exist_ok=True); frame.to_parquet(target, index=False)
    return frame


def _sample_slide_balanced(rows: pd.DataFrame, budget: int, seed: int) -> pd.DataFrame:
    if budget <= 0 or len(rows) <= budget:
        return rows
    rng, selected = np.random.default_rng(seed), []
    groups = [g.index.to_numpy() for _, g in rows.groupby("slide_id", sort=True)]
    # First draw ensures multi-slide patients are represented, then random fill.
    order = rng.permutation(len(groups))
    for index in order[:budget]:
        selected.append(rng.choice(groups[index]))
    if len(selected) < budget:
        remain = np.setdiff1d(rows.index.to_numpy(), np.asarray(selected), assume_unique=False)
        selected.extend(rng.choice(remain, budget - len(selected), replace=False).tolist())
    return rows.loc[np.asarray(selected)].sort_values("row_idx")


@dataclass(frozen=True)
class HoptimusPatchStore:
    """Validated reader for frozen H-Optimus tokens; never falls back to NPZ."""
    embeddings_path: Path
    metadata_path: Path
    index_path: Path | None = None
    expected_version: str = STORE_VERSION

    def __post_init__(self) -> None:
        object.__setattr__(self, "embeddings_path", Path(self.embeddings_path))
        object.__setattr__(self, "metadata_path", Path(self.metadata_path))
        if self.index_path is not None:
            object.__setattr__(self, "index_path", Path(self.index_path))

    @property
    def feature_h5(self) -> Path:
        """Backward-compatible name for the embedding file."""
        return self.embeddings_path

    @property
    def metadata_parquet(self) -> Path:
        return self.metadata_path

    @cached_property
    def _metadata_frame(self) -> pd.DataFrame:
        return pd.read_parquet(self.metadata_path)

    def _metadata(self) -> pd.DataFrame:
        """Cached immutable metadata; never reread the parquet per patient."""
        return self._metadata_frame

    @cached_property
    def _patient_rows(self) -> dict[str, pd.DataFrame]:
        return {str(patient): group.reset_index(drop=True) for patient, group in self._metadata_frame.groupby("patient_id", sort=False)}

    @staticmethod
    def _membership_filter(value: str) -> tuple[str, str]:
        scheme, partition = parse_membership(value)
        return split_membership_column(scheme), partition

    def validate(self, required_membership: str | None = None) -> dict[str, Any]:
        if not self.embeddings_path.exists() or not self.metadata_path.exists():
            raise FileNotFoundError("H-Optimus store requires an HDF5 embedding file and Parquet metadata file")
        metadata = self._metadata()
        missing = sorted(REQUIRED_METADATA_COLUMNS - set(metadata.columns))
        if missing:
            raise ValueError(f"Patch metadata is missing required columns: {missing}")
        expected_idx = np.arange(len(metadata), dtype=np.int64)
        if not np.array_equal(metadata.row_idx.to_numpy(dtype=np.int64), expected_idx):
            raise ValueError("metadata.row_idx must be ordered and contiguous from 0 to N-1")
        for col in ("internal_split", "external_split"):
            invalid = set(metadata[col].dropna().astype(str).str.lower()) - VALID_PARTITIONS
            if invalid:
                raise ValueError(f"Invalid {col} values: {sorted(invalid)}")
        if required_membership:
            col, partition = self._membership_filter(required_membership)
            if not (metadata[col].astype("string").str.lower() == partition).any():
                raise ValueError(f"No tokens exist for required membership {required_membership!r}")
        with h5py.File(self.embeddings_path, "r") as handle:
            if handle.attrs.get("store_format", "") != STORE_FORMAT:
                raise ValueError("Not a canonical MORPHEUS H-Optimus patch store")
            if str(handle.attrs.get("store_version", "")) != self.expected_version:
                raise ValueError(f"Unsupported patch store version {handle.attrs.get('store_version')!r}")
            if EMBEDDINGS_DATASET not in handle:
                raise KeyError("Patch HDF5 missing 'embeddings' dataset")
            dataset = handle[EMBEDDINGS_DATASET]
            if dataset.ndim != 2 or dataset.shape[1] != EXPECTED_HOPTIMUS_DIM:
                raise ValueError(f"Expected [N,{EXPECTED_HOPTIMUS_DIM}] embeddings, got {dataset.shape}")
            if int(handle.attrs.get("embedding_dim", -1)) != EXPECTED_HOPTIMUS_DIM:
                raise ValueError("HDF5 embedding_dim attribute does not identify H-Optimus-0 tokens")
            count = int(dataset.shape[0])
        if count != len(metadata):
            raise ValueError(f"Metadata rows ({len(metadata)}) != embedding rows ({count})")
        if self.index_path is not None:
            if not self.index_path.exists():
                raise FileNotFoundError(f"Missing patient token index {self.index_path}")
            index = pd.read_parquet(self.index_path)
            required_index = {"patient_id", "start_row", "end_row", "n_tokens"}
            if not required_index.issubset(index.columns):
                raise ValueError("Patient token index is missing required span columns")
        return {"n_patches": count, "n_slides": int(metadata.slide_id.nunique()),
                "n_patients": int(metadata.patient_id.nunique()), "embedding_dim": EXPECTED_HOPTIMUS_DIM,
                "store_version": self.expected_version}

    def load_patient_tokens(self, patient_id: str, max_tokens: int | None = 512, seed: int = 42,
                            slide_balanced: bool = True, required_membership: str | None = None) -> tuple[np.ndarray, pd.DataFrame]:
        """Return finite float32 ``[T,1536]`` tokens and aligned patch metadata.

        ``max_tokens`` of ``None`` or <=0 returns all available tokens.  The
        same seed gives the same selection; slide-balanced sampling represents
        every available slide before filling remaining budget.
        """
        return self.load_many_patient_tokens([patient_id], max_tokens=max_tokens, seed=seed,
                                             slide_balanced=slide_balanced, required_membership=required_membership)[0]

    def load_many_patient_tokens(self, patient_ids: Sequence[str], max_tokens: int | None = 512, seed: int = 42,
                                 slide_balanced: bool = True, required_membership: str | None = None) -> list[tuple[np.ndarray, pd.DataFrame]]:
        """Read a dynamic patient batch with one metadata cache and one HDF5 handle.

        This is the uncapped-bag fast path used by V2.  It preserves the exact
        semantics of :meth:`load_patient_tokens` while avoiding one Parquet and
        one HDF5 open per patient.
        """
        selected: list[pd.DataFrame] = []
        for offset, patient_id in enumerate(patient_ids):
            patient = normalize_tcga_patient_id(patient_id)
            rows = self._patient_rows.get(patient, self._metadata_frame.iloc[0:0]).copy()
            if required_membership:
                col, partition = self._membership_filter(required_membership)
                rows = rows.loc[rows[col].astype("string").str.lower() == partition].copy()
            budget = len(rows) if max_tokens is None else int(max_tokens)
            if budget > 0 and len(rows) > budget:
                rows = _sample_slide_balanced(rows, budget, seed + offset) if slide_balanced else rows.sample(n=budget, random_state=seed + offset).sort_values("row_idx")
            selected.append(rows.reset_index(drop=True))
        result: list[tuple[np.ndarray, pd.DataFrame]] = []
        with h5py.File(self.embeddings_path, "r") as handle:
            dataset = handle[EMBEDDINGS_DATASET]
            for rows in selected:
                if rows.empty:
                    result.append((np.zeros((0, EXPECTED_HOPTIMUS_DIM), dtype=np.float32), rows))
                    continue
                indices = rows.row_idx.to_numpy(dtype=np.int64)
                tokens = np.asarray(dataset[indices], dtype=np.float32)
                if tokens.shape != (len(rows), EXPECTED_HOPTIMUS_DIM) or not np.isfinite(tokens).all():
                    raise ValueError("Invalid token read: non-finite or metadata/embedding misalignment")
                result.append((tokens, rows))
        return result

    def load_patient(self, patient_id: str, max_patches: int = 512, seed: int = 42) -> tuple[np.ndarray, pd.DataFrame]:
        """Compatibility alias for older code."""
        return self.load_patient_tokens(patient_id, max_tokens=max_patches, seed=seed)


def build_patient_patch_index(metadata: pd.DataFrame) -> pd.DataFrame:
    """Build contiguous patient spans after canonical patient/slide/patch ordering."""
    if not np.array_equal(metadata.row_idx.to_numpy(dtype=np.int64), np.arange(len(metadata))):
        raise ValueError("Cannot index non-contiguous metadata")
    rows = []
    for patient, group in metadata.groupby("patient_id", sort=True):
        values = group.row_idx.to_numpy(dtype=np.int64)
        if not np.array_equal(values, np.arange(values[0], values[-1] + 1)):
            raise ValueError(f"Patient {patient} does not occupy one contiguous row span")
        rows.append({"patient_id": patient, "start_row": int(values[0]), "end_row": int(values[-1] + 1),
                     "n_tokens": int(len(values)), "n_slides": int(group.slide_id.nunique())})
    return pd.DataFrame(rows)


def write_canonical_patch_store(output_dir: str | Path, embeddings: np.ndarray, metadata: pd.DataFrame,
                                provenance: dict[str, Any], overwrite: bool = False) -> HoptimusPatchStore:
    """Atomically write a canonical store and its patient-span index."""
    target = Path(output_dir); target.mkdir(parents=True, exist_ok=True)
    h5_path = target / "hoptimus_patch_embeddings.h5"
    meta_path = target / "hoptimus_patch_metadata.parquet"
    index_path = target / "hoptimus_patient_patch_index.parquet"
    manifest_path = target / "hoptimus_patch_manifest.json"
    if not overwrite and any(p.exists() for p in (h5_path, meta_path, index_path, manifest_path)):
        raise FileExistsError(f"Refusing to overwrite patch store in {target}")
    values = np.asarray(embeddings, dtype=np.float32)
    if values.ndim != 2 or values.shape[1] != EXPECTED_HOPTIMUS_DIM or not np.isfinite(values).all():
        raise ValueError(f"Embeddings must be finite [N,{EXPECTED_HOPTIMUS_DIM}] float32")
    frame = metadata.copy().reset_index(drop=True)
    frame["_input_row"] = np.arange(len(frame), dtype=np.int64)
    frame["row_idx"] = np.arange(len(frame), dtype=np.int64)
    missing = REQUIRED_METADATA_COLUMNS - set(frame.columns)
    if missing:
        raise ValueError(f"Cannot write store; metadata missing {sorted(missing)}")
    if len(frame) != len(values):
        raise ValueError("Embedding and metadata row counts differ")
    # Contiguous grouping is required for inexpensive index reads.  Reorder the
    # embedding matrix by the same input-row permutation, never metadata alone.
    frame = frame.sort_values(["patient_id", "slide_id", "patch_id", "source_sha256"]).reset_index(drop=True)
    values = values[frame["_input_row"].to_numpy(dtype=np.int64)]
    frame = frame.drop(columns="_input_row")
    frame["row_idx"] = np.arange(len(frame), dtype=np.int64)
    tmp_h5, tmp_meta, tmp_index = h5_path.with_suffix(".tmp"), meta_path.with_suffix(".tmp"), index_path.with_suffix(".tmp")
    with h5py.File(tmp_h5, "w") as handle:
        handle.attrs.update({"store_format": STORE_FORMAT, "store_version": STORE_VERSION,
                             "embedding_dim": EXPECTED_HOPTIMUS_DIM, "encoder_name": "H-Optimus-0"})
        handle.create_dataset(EMBEDDINGS_DATASET, data=values, dtype="float32",
                              chunks=(min(max(len(values), 1), 1024), EXPECTED_HOPTIMUS_DIM), compression="gzip", shuffle=True)
    frame.to_parquet(tmp_meta, index=False)
    build_patient_patch_index(frame).to_parquet(tmp_index, index=False)
    os.replace(tmp_h5, h5_path); os.replace(tmp_meta, meta_path); os.replace(tmp_index, index_path)
    manifest = {"store_format": STORE_FORMAT, "store_version": STORE_VERSION, "encoder_name": "H-Optimus-0",
                "embedding_dim": EXPECTED_HOPTIMUS_DIM, "n_patches": int(len(frame)),
                "files": {"embeddings": h5_path.name, "metadata": meta_path.name, "patient_index": index_path.name},
                "provenance": provenance}
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    store = HoptimusPatchStore(h5_path, meta_path, index_path); store.validate(); return store


def write_patch_extraction_manifest(*args: Any, **kwargs: Any) -> Path:
    """Deprecated planning helper retained for the old manifest script."""
    output = Path(kwargs.get("output_path", args[0] if args else "data/processed/wsi/tcga_ut_hoptimus0_patch_extraction_manifest.json"))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"store_format": STORE_FORMAT, "store_version": STORE_VERSION,
                                  "status": "use scripts/extract_hoptimus_tcga_ut.py"}, indent=2), encoding="utf-8")
    return output
