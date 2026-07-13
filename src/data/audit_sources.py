"""Audit existing v1 data sources without mutating source data."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import base_manifest, write_json


def _file_record(path: Path) -> dict[str, Any]:
    exists = path.exists()
    record: dict[str, Any] = {
        "path": str(path),
        "exists": exists,
        "is_file": path.is_file() if exists else False,
        "is_dir": path.is_dir() if exists else False,
    }
    if exists and path.is_file():
        record["bytes"] = path.stat().st_size
    return record


def _count_files(path: Path, patterns: tuple[str, ...]) -> dict[str, Any]:
    if not path.exists():
        return {"path": str(path), "exists": False, "counts": {}, "bytes": 0}
    counts: dict[str, int] = {}
    total = 0
    for pattern in patterns:
        files = [p for p in path.rglob(pattern) if p.is_file()]
        counts[pattern] = len(files)
        total += sum(p.stat().st_size for p in files)
    return {"path": str(path), "exists": True, "counts": counts, "bytes": total}


def build_inventory(config_path: str = "morpheus/configs/v1.json") -> dict[str, Any]:
    cfg = load_config(config_path)
    seed = int(cfg.raw.get("seed", 42))
    inventory = base_manifest(cfg.project_root, cfg.config_path, seed)
    inventory["configured_sources"] = {
        key: _file_record(cfg.path(key))
        for key in cfg.raw.get("paths", {})
        if key not in {"processed_dir", "outputs_dir", "morpheus_root"}
    }
    inventory["modalities"] = {
        "wsi_hoptimus": _count_files(cfg.path("wsi_hoptimus_dir"), ("*.npy", "*.parquet")),
        "wsi_patch_npz": _count_files(cfg.path("wsi_patch_dir"), ("*.npz",)),
        "processed_omics": _count_files(cfg.path("meta_intersurv_data") / "omics", ("*.parquet", "*.tsv", "*.gz", "*.h5ad")),
        "gene_sets": _count_files(cfg.path("meta_intersurv_data") / "msigdb", ("*.gmt", "*.csv")),
    }
    rna_model = dict(cfg.section("foundation_models").get("rna_primary", {}))
    inventory["bulkformer_status"] = {
        "repo_id_configured": bool(rna_model.get("repo_id")),
        "repo_id": rna_model.get("repo_id", ""),
        "official_v1_blocked_until_configured": not bool(rna_model.get("repo_id")),
    }
    return inventory


def write_readiness_report(inventory: dict[str, Any], output_dir: Path) -> None:
    lines = [
        "# V1 Data Readiness Report",
        "",
        f"Git commit: `{inventory.get('git_commit')}`",
        "",
        "## Source Status",
    ]
    for name, record in inventory["configured_sources"].items():
        status = "present" if record["exists"] else "missing"
        lines.append(f"- `{name}`: {status} - `{record['path']}`")
    lines.extend(["", "## Modality Counts"])
    for name, record in inventory["modalities"].items():
        lines.append(f"- `{name}`: {record.get('counts', {})}, bytes={record.get('bytes', 0)}")
    if inventory["bulkformer_status"]["official_v1_blocked_until_configured"]:
        lines.extend(
            [
                "",
                "## Blocking Item",
                "BulkFormer is required for the official v1 proof, but `foundation_models.rna_primary.repo_id` is blank.",
                "Existing Geneformer/processed RNA may be used only for fallback sanity checks.",
            ]
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "v1_readiness_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    args = parser.parse_args()
    cfg = load_config(args.config)
    out = cfg.path("outputs_dir") / "v1_data_audit"
    inventory = build_inventory(args.config)
    write_json(out / "data_inventory.json", inventory)
    write_readiness_report(inventory, out)
    print(f"Wrote audit to {out}")


if __name__ == "__main__":
    main()
