"""Zenodo acquisition helpers for BulkFormer artifacts."""

from __future__ import annotations

import hashlib
import json
import os
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen

from morpheus.src.utils.provenance import write_json


ZENODO_API = "https://zenodo.org/api/records/{record_id}"
DEFAULT_PRIMARY_RECORD = "15744294"
DEFAULT_COMPAT_RECORD = "15559368"
MANIFEST_NAME = "bulkformer_zenodo_manifest.json"

REQUIRED_CATEGORIES = (
    "tcga_h5ad",
    "gene_info",
    "gene_length",
    "checkpoint",
    "graph",
    "graph_weight",
    "gene_embedding",
    "interested_gene_list",
)


@dataclass(frozen=True)
class ZenodoFile:
    record_id: str
    record_title: str
    record_doi: str | None
    version: str | None
    filename: str
    size: int | None
    checksum: str | None
    download_url: str
    category: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "record_id": self.record_id,
            "record_title": self.record_title,
            "record_doi": self.record_doi,
            "version": self.version,
            "filename": self.filename,
            "size": self.size,
            "checksum": self.checksum,
            "download_url": self.download_url,
            "category": self.category,
        }


def _read_json_url(url: str, token: str | None = None) -> dict[str, Any]:
    headers = {"Accept": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urlopen(Request(url, headers=headers), timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def categorize_bulkformer_file(filename: str) -> str:
    name = filename.lower()
    if name.endswith(".pt") and ("bulkformer" in name or "ckpt" in name or "epoch" in name):
        return "checkpoint"
    if name == "g_tcga.pt":
        return "graph"
    if name == "g_tcga_weight.pt":
        return "graph_weight"
    if name == "esm2_feature_concat.pt":
        return "gene_embedding"
    if name == "interested_gene_list.pt" or name == "high_var_gene_list.pt":
        return "interested_gene_list"
    if name == "gene_length_df.csv":
        return "gene_length"
    if name.endswith(".csv") and ("gene" in name or "vocab" in name):
        return "gene_info"
    if name.endswith(".h5ad") and "tcga" in name:
        return "tcga_h5ad"
    if name.endswith(".zip") and ("code" in name or "bulkformer" in name):
        return "code_archive"
    if name.endswith(".h5ad"):
        return "h5ad_optional"
    if name.endswith(".zip"):
        return "zip_optional"
    return "optional"


def inventory_zenodo_record(record_id: str, token: str | None = None) -> list[ZenodoFile]:
    payload = _read_json_url(ZENODO_API.format(record_id=record_id), token)
    title = str(payload.get("metadata", {}).get("title", ""))
    doi = payload.get("doi")
    version = payload.get("metadata", {}).get("version")
    files: list[ZenodoFile] = []
    for item in payload.get("files", []):
        filename = str(item.get("key") or item.get("filename") or "")
        links = item.get("links", {})
        download_url = links.get("self") or links.get("download")
        if not filename or not download_url:
            continue
        files.append(
            ZenodoFile(
                record_id=str(record_id),
                record_title=title,
                record_doi=doi,
                version=version,
                filename=filename,
                size=item.get("size"),
                checksum=item.get("checksum"),
                download_url=str(download_url),
                category=categorize_bulkformer_file(filename),
            )
        )
    return files


def inventory_bulkformer_records(
    primary_record_id: str = DEFAULT_PRIMARY_RECORD,
    compat_record_id: str | None = DEFAULT_COMPAT_RECORD,
    token: str | None = None,
) -> list[ZenodoFile]:
    record_ids = [primary_record_id]
    if compat_record_id and compat_record_id != primary_record_id:
        record_ids.append(compat_record_id)
    files: list[ZenodoFile] = []
    for record_id in record_ids:
        files.extend(inventory_zenodo_record(record_id, token))
    return files


def select_bulkformer_files(files: list[ZenodoFile], primary_record_id: str = DEFAULT_PRIMARY_RECORD) -> dict[str, ZenodoFile]:
    grouped: dict[str, list[ZenodoFile]] = {}
    for file in files:
        grouped.setdefault(file.category, []).append(file)

    selected: dict[str, ZenodoFile] = {}
    for category in REQUIRED_CATEGORIES:
        candidates = grouped.get(category, [])
        if not candidates:
            continue
        candidates = sorted(
            candidates,
            key=lambda f: (
                f.record_id != primary_record_id,
                "v2" not in f.filename.lower(),
                f.filename.lower(),
            ),
        )
        selected[category] = candidates[0]
    return selected


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_file(file: ZenodoFile, destination: Path, token: str | None = None) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and file.size is not None and destination.stat().st_size == int(file.size):
        return "already_present"
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urlopen(Request(file.download_url, headers=headers), timeout=120) as response, destination.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    return "downloaded"


def download_bulkformer_artifacts(
    cache_dir: str | Path,
    primary_record_id: str = DEFAULT_PRIMARY_RECORD,
    compat_record_id: str | None = DEFAULT_COMPAT_RECORD,
    token: str | None = None,
    dry_run: bool = False,
) -> Path:
    cache = Path(cache_dir)
    files = inventory_bulkformer_records(primary_record_id, compat_record_id, token)
    selected = select_bulkformer_files(files, primary_record_id)
    downloads: dict[str, dict[str, Any]] = {}
    for category, file in selected.items():
        destination = cache / file.record_id / file.filename
        status = "dry_run" if dry_run else _download_file(file, destination, token)
        downloads[category] = {**file.to_dict(), "local_path": str(destination), "status": status}
    missing = [category for category in REQUIRED_CATEGORIES if category not in selected]
    manifest = {
        "primary_record_id": primary_record_id,
        "compat_record_id": compat_record_id,
        "records": sorted({file.record_id for file in files}),
        "selected": downloads,
        "missing_required": missing,
        "mixed_release": sorted({item["record_id"] for item in downloads.values()}),
        "code_source": {
            "kind": "github",
            "repo": "KangBoming/BulkFormer",
            "url": "https://github.com/KangBoming/BulkFormer",
            "note": "Zenodo v2/v1 API inventory did not expose a BulkFormer code archive; use the upstream GitHub repository for adapter inspection.",
        },
        "all_files": [file.to_dict() for file in files],
    }
    manifest_path = cache / MANIFEST_NAME
    write_json(manifest_path, manifest)
    return manifest_path


def verify_download_manifest(manifest_path: str | Path) -> dict[str, Any]:
    path = Path(manifest_path)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    selected = manifest.get("selected", {})
    missing_required = [category for category in REQUIRED_CATEGORIES if category not in selected]
    missing_files: list[str] = []
    checksum_failures: list[str] = []
    size_failures: list[str] = []
    for category, record in selected.items():
        local_path = Path(record.get("local_path", ""))
        if not local_path.exists():
            missing_files.append(category)
            continue
        expected_size = record.get("size")
        if expected_size is not None and local_path.stat().st_size != int(expected_size):
            size_failures.append(category)
        checksum = record.get("checksum")
        if isinstance(checksum, str) and checksum.startswith("md5:"):
            expected = checksum.split(":", 1)[1].lower()
            if _md5(local_path).lower() != expected:
                checksum_failures.append(category)
    return {
        "manifest_path": str(path),
        "missing_required": sorted(set(missing_required + list(manifest.get("missing_required", [])))),
        "missing_files": missing_files,
        "size_failures": size_failures,
        "checksum_failures": checksum_failures,
        "ready": not missing_required and not missing_files and not size_failures and not checksum_failures,
        "mixed_release": manifest.get("mixed_release", []),
    }


def extract_code_archives(manifest_path: str | Path, output_dir: str | Path) -> list[Path]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    selected = manifest.get("selected", {})
    extracted: list[Path] = []
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for category, record in selected.items():
        if category != "code_archive":
            continue
        archive = Path(record["local_path"])
        if not archive.exists():
            raise FileNotFoundError(archive)
        with zipfile.ZipFile(archive) as zf:
            for member in zf.infolist():
                target = (output / member.filename).resolve()
                if not str(target).startswith(str(output.resolve())):
                    raise ValueError(f"Unsafe archive member: {member.filename}")
            zf.extractall(output)
            extracted.append(output / archive.stem)
    return extracted


def env_token() -> str | None:
    return os.environ.get("ZENODO_ACCESS_TOKEN") or None
