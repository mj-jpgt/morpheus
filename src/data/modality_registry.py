"""Build modality registry and initial patient table outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from morpheus.src.utils.ids import normalize_patient_id
from morpheus.src.utils.provenance import write_json


MODALITY_ORDER = [
    "rna",
    "snv",
    "cnv",
    "clinical",
    "proteomics",
    "phosphoproteomics",
    "acetylproteomics",
    "mirna",
    "circrna",
    "hla",
    "purity",
    "immune",
    "wsi",
    "unknown",
]


def build_modality_registry(inventory: pd.DataFrame) -> dict[str, Any]:
    registry: dict[str, Any] = {"modalities": {}}
    for modality in MODALITY_ORDER:
        subset = inventory[inventory["modality"] == modality] if "modality" in inventory else pd.DataFrame()
        if subset.empty:
            registry["modalities"][modality] = {
                "available_in_inventory": False,
                "n_files": 0,
                "size_bytes": 0,
                "sample_paths": [],
                "table_files_first": [],
            }
            continue
        table_subset = subset[subset["is_table_file"].fillna(False)] if "is_table_file" in subset else subset.iloc[0:0]
        registry["modalities"][modality] = {
            "available_in_inventory": True,
            "n_files": int(len(subset)),
            "size_bytes": int(subset["size_bytes"].fillna(0).sum()) if "size_bytes" in subset else 0,
            "sample_paths": subset["path"].head(10).tolist(),
            "table_files_first": table_subset["path"].head(25).tolist(),
            "download_policy": "tables_first_no_wsi" if modality != "wsi" else "metadata_only_until_mapping",
        }
    return registry


def write_modality_registry(inventory_path: str | Path, output_path: str | Path) -> dict[str, Any]:
    inventory = pd.read_parquet(inventory_path)
    registry = build_modality_registry(inventory)
    write_json(Path(output_path), registry)
    return registry


def _existing_master(existing_master_path: Path | None) -> pd.DataFrame:
    if existing_master_path and existing_master_path.exists():
        master = pd.read_parquet(existing_master_path).copy()
        if "patient_id" in master:
            master["patient_id"] = master["patient_id"].map(normalize_patient_id)
        return master.dropna(subset=["patient_id"]).drop_duplicates("patient_id")
    return pd.DataFrame(columns=["patient_id"])


def build_master_patient_table_v0(
    output_path: str | Path,
    existing_master_path: str | Path | None = "morpheus/data/processed/master_patient_table.parquet",
    hoptimus_patient_index: str | Path | None = "meta-intersurv/data/tcga_ut/embeddings_hoptimus0/patient_index.parquet",
) -> Path:
    master = _existing_master(Path(existing_master_path) if existing_master_path else None)
    if "patient_id" not in master:
        master["patient_id"] = []
    wsi_patients: set[str] = set()
    if hoptimus_patient_index and Path(hoptimus_patient_index).exists():
        wsi = pd.read_parquet(hoptimus_patient_index)
        col = next((c for c in wsi.columns if "patient" in str(c).lower()), None)
        if col:
            wsi_patients = {pid for pid in wsi[col].map(normalize_patient_id) if pid}
    if wsi_patients:
        master = pd.concat([master, pd.DataFrame({"patient_id": sorted(wsi_patients)})], ignore_index=True)
    master = master.dropna(subset=["patient_id"]).drop_duplicates("patient_id").sort_values("patient_id")
    master["has_wsi_hoptimus0"] = master["patient_id"].isin(wsi_patients)
    for col in ["has_rna", "has_snv", "has_cnv", "has_clinical", "has_wsi"]:
        if col not in master:
            master[col] = False
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    master.to_parquet(output, index=False)
    return output
