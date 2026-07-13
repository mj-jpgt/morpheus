"""Hugging Face dataset inventory without bulk downloads."""

from __future__ import annotations

import json
import os
from pathlib import Path
from collections import deque
from typing import Any
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas as pd

from morpheus.src.utils.provenance import write_json


DEFAULT_DATASET = "mj0jpgg/tcga_cptac_cgga"
HF_ROOT = "https://huggingface.co"


def classify_modality(path: str) -> str:
    lower = path.lower()
    suffix = Path(lower).suffix
    if suffix in {".svs", ".tif", ".tiff", ".jpg", ".jpeg", ".png"} or "wsi" in lower or "slide" in lower:
        return "wsi"
    if "phospho" in lower:
        return "phosphoproteomics"
    if "acetyl" in lower:
        return "acetylproteomics"
    if "proteom" in lower or "protein" in lower:
        return "proteomics"
    if "mutation" in lower or "maf" in lower or "snv" in lower or "variant" in lower:
        return "snv"
    if "cnv" in lower or "copy" in lower:
        return "cnv"
    if "mirna" in lower:
        return "mirna"
    if "circrna" in lower or "circ_rna" in lower or "circular" in lower:
        return "circrna"
    if "rna" in lower or "transcript" in lower or "expression" in lower:
        return "rna"
    if "hla" in lower:
        return "hla"
    if "cibersort" in lower or "xcell" in lower or "immune" in lower:
        return "immune"
    if "purity" in lower:
        return "purity"
    if "clinical" in lower or "follow" in lower or "history" in lower or "survival" in lower:
        return "clinical"
    return "unknown"


def is_table_file(path: str) -> bool:
    return Path(path.lower()).suffix in {".csv", ".tsv", ".txt", ".json", ".parquet", ".xlsx", ".h5ad", ".h5"}


def _read_json_url(url: str, token: str | None = None) -> Any:
    headers = {"Accept": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urlopen(Request(url, headers=headers), timeout=90) as response:
        return json.loads(response.read().decode("utf-8"))


def _tree_url(repo_id: str, revision: str, path: str | None = None, expand: bool = True) -> str:
    encoded_repo = quote(repo_id, safe="/")
    encoded_revision = quote(revision, safe="")
    suffix = "" if not path else "/" + quote(path.strip("/"), safe="/")
    expand_flag = "?expand=1" if expand else ""
    return f"{HF_ROOT}/api/datasets/{encoded_repo}/tree/{encoded_revision}{suffix}{expand_flag}"


def _is_wsi_heavy_path(path: str) -> bool:
    lower = path.lower()
    return any(token in lower for token in ("wsi", "slide", "svs", "image", "tiles", "patch"))


def _item_rows(repo_id: str, revision: str, item: dict[str, Any]) -> dict[str, Any] | None:
    if item.get("type") not in {None, "file"}:
        return None
    path = str(item.get("path", ""))
    if not path:
        return None
    lfs = item.get("lfs") or {}
    size = item.get("size", lfs.get("size"))
    return {
        "repo_id": repo_id,
        "revision": revision,
        "path": path,
        "filename": Path(path).name,
        "extension": Path(path).suffix.lower(),
        "size_bytes": int(size) if size is not None else None,
        "oid": item.get("oid"),
        "lfs_oid": lfs.get("oid"),
        "modality": classify_modality(path),
        "is_table_file": is_table_file(path),
        "is_wsi_like": classify_modality(path) == "wsi",
        "download_url": f"{HF_ROOT}/datasets/{repo_id}/resolve/{revision}/{path}",
    }


def inventory_hf_dataset(
    repo_id: str = DEFAULT_DATASET,
    revision: str = "main",
    token: str | None = None,
    limit: int | None = None,
    max_files: int = 5000,
    include_wsi_dirs: bool = False,
) -> pd.DataFrame:
    queue: deque[str | None] = deque([None])
    rows: list[dict[str, Any]] = []
    visited_dirs: set[str | None] = set()
    while queue and len(rows) < max_files:
        path_prefix = queue.popleft()
        if path_prefix in visited_dirs:
            continue
        visited_dirs.add(path_prefix)
        payload = _read_json_url(_tree_url(repo_id, revision, path_prefix), token)
        if not isinstance(payload, list):
            raise RuntimeError(f"Unexpected Hugging Face tree response for {repo_id}: {type(payload).__name__}")
        for item in payload:
            item_type = item.get("type")
            item_path = str(item.get("path", ""))
            if item_type == "directory":
                if include_wsi_dirs or not _is_wsi_heavy_path(item_path):
                    queue.append(item_path)
                else:
                    rows.append(
                        {
                            "repo_id": repo_id,
                            "revision": revision,
                            "path": item_path,
                            "filename": Path(item_path).name,
                            "extension": "",
                            "size_bytes": None,
                            "oid": item.get("oid"),
                            "lfs_oid": None,
                            "modality": "wsi",
                            "is_table_file": False,
                            "is_wsi_like": True,
                            "download_url": "",
                            "inventory_note": "directory_skipped_by_default",
                        }
                    )
                continue
            row = _item_rows(repo_id, revision, item)
            if row is None:
                continue
            row["inventory_note"] = "file"
            rows.append(row)
            if limit is not None and len(rows) >= limit:
                break
            if len(rows) >= max_files:
                break
        if limit is not None and len(rows) >= limit:
            continue
    return pd.DataFrame(rows)


def write_inventory_outputs(
    frame: pd.DataFrame,
    inventory_path: str | Path,
    summary_path: str | Path | None = None,
) -> None:
    inventory = Path(inventory_path)
    inventory.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(inventory, index=False)
    if summary_path:
        by_modality = (
            frame.groupby("modality", dropna=False)
            .agg(files=("path", "count"), bytes=("size_bytes", "sum"))
            .reset_index()
            .sort_values("modality")
        )
        write_json(
            Path(summary_path),
            {
                "inventory_path": str(inventory),
                "n_files": int(len(frame)),
                "total_size_bytes": int(frame["size_bytes"].fillna(0).sum()) if "size_bytes" in frame else 0,
                "modalities": by_modality.to_dict(orient="records"),
            },
        )


def hf_token_from_env() -> str | None:
    return os.environ.get("HF_TOKEN") or os.environ.get("HUGGINGFACE_TOKEN") or None
