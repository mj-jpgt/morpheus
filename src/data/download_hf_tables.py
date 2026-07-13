"""Download selected non-WSI tables from the HF multimodal mirror."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from morpheus.src.utils.provenance import write_json


DEFAULT_MODALITIES = [
    "clinical",
    "purity",
    "immune",
    "hla",
    "proteomics",
    "phosphoproteomics",
    "acetylproteomics",
    "mirna",
    "circrna",
    "snv",
    "cnv",
    "rna",
]


def _download(url: str, out_path: Path, token: str | None = None) -> str:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        return "already_present"
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urlopen(Request(url, headers=headers), timeout=120) as response, out_path.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    return "downloaded"


def select_table_files(
    inventory_path: str | Path,
    modalities: list[str],
    max_total_bytes: int | None = None,
) -> pd.DataFrame:
    inv = pd.read_parquet(inventory_path)
    table = inv[
        inv["modality"].isin(modalities)
        & inv["is_table_file"].fillna(False)
        & ~inv["is_wsi_like"].fillna(False)
        & inv["download_url"].astype(str).str.len().gt(0)
    ].copy()
    table = table.sort_values(["modality", "size_bytes", "path"], na_position="last")
    if max_total_bytes is not None:
        keep = []
        total = 0
        for idx, row in table.iterrows():
            size = int(row["size_bytes"]) if pd.notna(row["size_bytes"]) else 0
            if total + size > max_total_bytes:
                continue
            keep.append(idx)
            total += size
        table = table.loc[keep]
    return table


def download_selected_tables(
    inventory_path: str | Path,
    output_dir: str | Path = "data/raw/hf_tcga_cptac_cgga",
    modalities: list[str] | None = None,
    max_total_bytes: int | None = None,
    dry_run: bool = False,
    token: str | None = None,
) -> Path:
    modalities = modalities or DEFAULT_MODALITIES
    selected = select_table_files(inventory_path, modalities, max_total_bytes)
    root = Path(output_dir)
    records = []
    for _, row in selected.iterrows():
        rel = Path(str(row["path"]))
        out = root / rel
        status = "dry_run" if dry_run else _download(str(row["download_url"]), out, token)
        records.append(
            {
                "path": str(row["path"]),
                "modality": str(row["modality"]),
                "size_bytes": int(row["size_bytes"]) if pd.notna(row["size_bytes"]) else None,
                "download_url": str(row["download_url"]),
                "local_path": str(out),
                "status": status,
            }
        )
    manifest = root / "download_manifest.json"
    write_json(
        manifest,
        {
            "inventory_path": str(inventory_path),
            "output_dir": str(root),
            "modalities": modalities,
            "dry_run": dry_run,
            "n_files": len(records),
            "total_size_bytes": int(sum(r["size_bytes"] or 0 for r in records)),
            "files": records,
        },
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", default="data/manifests/hf_tcga_cptac_cgga_inventory.parquet")
    parser.add_argument("--output-dir", default="data/raw/hf_tcga_cptac_cgga")
    parser.add_argument("--modalities", default="clinical,purity,immune,hla")
    parser.add_argument("--max-total-gb", type=float, default=1.0)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    modalities = [x.strip() for x in args.modalities.split(",") if x.strip()]
    max_bytes = int(args.max_total_gb * 1024**3) if args.max_total_gb is not None else None
    token = os.environ.get("HF_TOKEN") or os.environ.get("HUGGINGFACE_TOKEN")
    print(download_selected_tables(args.inventory, args.output_dir, modalities, max_bytes, args.dry_run, token))


if __name__ == "__main__":
    main()
