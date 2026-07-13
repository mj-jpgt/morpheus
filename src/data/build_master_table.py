"""Build a patient-level master table from existing TCGA artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd

from morpheus.src.utils.config import load_config
from morpheus.src.utils.ids import normalize_patient_id, parse_tcga_barcode
from morpheus.src.utils.provenance import base_manifest, write_json


def _read_parquet_if_exists(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_parquet(path)


def _patient_column(df: pd.DataFrame) -> str | None:
    candidates = ["Patient ID", "patient_id", "Patient_ID", "case_id", "submitter_id"]
    for col in candidates:
        if col in df.columns:
            return col
    for col in df.columns:
        lower = str(col).lower()
        if "patient" in lower or "barcode" in lower:
            return col
    return None


def _patients_from_table(path: Path) -> set[str]:
    df = _read_parquet_if_exists(path)
    if df.empty:
        return set()
    col = _patient_column(df)
    if col is None:
        return set()
    return {pid for pid in (normalize_patient_id(v) for v in df[col]) if pid}


def _wsi_hoptimus_patients(path: Path) -> set[str]:
    patients: set[str] = set()
    for meta_name in ("train.meta.parquet", "valid.meta.parquet", "test.meta.parquet", "patient_index.parquet"):
        meta = path / meta_name
        if not meta.exists():
            continue
        df = pd.read_parquet(meta)
        col = _patient_column(df)
        if col is None:
            continue
        patients.update(pid for pid in (normalize_patient_id(v) for v in df[col]) if pid)
    return patients


def _wsi_patch_rows(path: Path) -> pd.DataFrame:
    rows = []
    for fp in sorted(path.glob("*.npz")):
        bc = parse_tcga_barcode(fp.name)
        if bc.patient_id:
            rows.append(
                {
                    "patient_id": bc.patient_id,
                    "slide_id": fp.stem,
                    "wsi_patch_path": str(fp),
                }
            )
    return pd.DataFrame(rows)


def _first_existing_column(df: pd.DataFrame, names: Iterable[str]) -> str | None:
    for name in names:
        if name in df.columns:
            return name
    return None


def build_master_table(config_path: str = "morpheus/configs/v1.json") -> dict[str, Path]:
    cfg = load_config(config_path)
    processed = cfg.path("processed_dir")
    processed.mkdir(parents=True, exist_ok=True)

    rna_patients = _patients_from_table(cfg.path("rna_processed"))
    snv_patients = _patients_from_table(cfg.path("snv_processed"))
    cnv_patients = _patients_from_table(cfg.path("cnv_processed"))
    clinical_patients = _patients_from_table(cfg.path("clinical_processed"))
    geneformer_patients = _patients_from_table(cfg.path("rna_geneformer_embeddings"))
    hoptimus_patients = _wsi_hoptimus_patients(cfg.path("wsi_hoptimus_dir"))
    wsi_rows = _wsi_patch_rows(cfg.path("wsi_patch_dir"))
    patch_patients = set(wsi_rows["patient_id"]) if not wsi_rows.empty else set()

    all_patients = sorted(
        rna_patients
        | snv_patients
        | cnv_patients
        | clinical_patients
        | geneformer_patients
        | hoptimus_patients
        | patch_patients
    )
    master = pd.DataFrame({"patient_id": all_patients})
    master["has_rna"] = master["patient_id"].isin(rna_patients)
    master["has_rna_geneformer"] = master["patient_id"].isin(geneformer_patients)
    master["has_snv"] = master["patient_id"].isin(snv_patients)
    master["has_cnv"] = master["patient_id"].isin(cnv_patients)
    master["has_clinical"] = master["patient_id"].isin(clinical_patients)
    master["has_wsi_hoptimus"] = master["patient_id"].isin(hoptimus_patients)
    master["has_wsi_patch"] = master["patient_id"].isin(patch_patients)
    master["has_wsi"] = master["has_wsi_hoptimus"] | master["has_wsi_patch"]

    clinical = _read_parquet_if_exists(cfg.path("clinical_processed"))
    if not clinical.empty:
        col = _patient_column(clinical)
        if col:
            clinical = clinical.copy()
            clinical["patient_id"] = clinical[col].map(normalize_patient_id)
            cancer_col = _first_existing_column(
                clinical,
                ["Cancer Type Acronym", "cancer_type", "project_id", "type"],
            )
            time_col = next((c for c in clinical.columns if "Overall Survival" in str(c) and "Month" in str(c)), None)
            event_col = next((c for c in clinical.columns if "Overall Survival Status" in str(c)), None)
            keep = ["patient_id"]
            for maybe in (cancer_col, time_col, event_col):
                if maybe:
                    keep.append(maybe)
            clinical_small = clinical[keep].drop_duplicates("patient_id")
            rename = {}
            if cancer_col:
                rename[cancer_col] = "cancer_type"
            if time_col:
                rename[time_col] = "survival_time"
            if event_col:
                rename[event_col] = "survival_event_raw"
            master = master.merge(clinical_small.rename(columns=rename), on="patient_id", how="left")
    if "cancer_type" not in master.columns:
        master["cancer_type"] = master["patient_id"].str[5:7]
    if "survival_event_raw" in master.columns:
        raw = master["survival_event_raw"].astype(str).str.lower()
        master["survival_event"] = raw.str.contains("deceased|dead|1:true|true|event").fillna(False).astype(int)

    master["modality_count"] = master[["has_rna", "has_snv", "has_cnv", "has_clinical", "has_wsi"]].sum(axis=1)
    master_path = processed / "master_patient_table.parquet"
    master.to_parquet(master_path, index=False)

    missingness = (
        master.groupby(["has_wsi", "has_rna", "has_snv", "has_cnv", "has_clinical"], dropna=False)
        .size()
        .reset_index(name="n_patients")
    )
    missingness_path = processed / "modality_missingness.parquet"
    missingness.to_parquet(missingness_path, index=False)

    slides_path = processed / "patient_to_slides.parquet"
    if wsi_rows.empty:
        pd.DataFrame(columns=["patient_id", "slide_id", "wsi_patch_path"]).to_parquet(slides_path, index=False)
    else:
        wsi_rows.to_parquet(slides_path, index=False)

    manifest = base_manifest(cfg.project_root, cfg.config_path, int(cfg.raw.get("seed", 42)))
    manifest["outputs"] = {k: str(v) for k, v in {"master": master_path, "missingness": missingness_path, "slides": slides_path}.items()}
    manifest["counts"] = {
        "patients": int(len(master)),
        "wsi_and_rna": int((master["has_wsi"] & master["has_rna"]).sum()),
        "wsi_rna_clinical": int((master["has_wsi"] & master["has_rna"] & master["has_clinical"]).sum()),
        "full_wsi_rna_snv_cnv_clinical": int(
            (master["has_wsi"] & master["has_rna"] & master["has_snv"] & master["has_cnv"] & master["has_clinical"]).sum()
        ),
    }
    write_json(processed / "master_patient_table.manifest.json", manifest)
    return {"master": master_path, "missingness": missingness_path, "slides": slides_path}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    args = parser.parse_args()
    outputs = build_master_table(args.config)
    print(f"Wrote master table: {outputs['master']}")


if __name__ == "__main__":
    main()
