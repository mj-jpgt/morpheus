"""Build patient registry and locked splits for tumor-state modeling."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd

from morpheus.src.utils.config import load_config
from morpheus.src.utils.ids import normalize_patient_id
from morpheus.src.utils.provenance import base_manifest, write_json


def _read_parquet(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.exists() else pd.DataFrame()


def _patient_ids_from_index_or_column(frame: pd.DataFrame) -> set[str]:
    if frame.empty:
        return set()
    if "patient_id" in frame.columns:
        values = frame["patient_id"]
    else:
        values = frame.index.to_series()
    return {pid for pid in values.map(normalize_patient_id).dropna().astype(str)}


def _wsi_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with h5py.File(path, "r") as handle:
        return {x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in handle["patient_ids"][:]}


def _bulkformer_metadata(rna_path: Path) -> pd.DataFrame:
    meta_path = rna_path.with_name("tcga_bulkformer_embedding_metadata.parquet")
    if not meta_path.exists():
        return pd.DataFrame(columns=["patient_id", "cancer_type"])
    meta = pd.read_parquet(meta_path).copy()
    meta["patient_id"] = meta["patient_id"].map(normalize_patient_id)
    meta["sample_type"] = meta["sample_id"].astype(str).str.split("-").str[3].str[:2]
    preferred = meta[meta["sample_type"].isin(["01", "03"])].copy()
    if preferred.empty:
        preferred = meta
    rows = []
    for patient_id, group in preferred.dropna(subset=["patient_id"]).groupby("patient_id"):
        mode = group["cancer_type"].astype(str).mode()
        rows.append({"patient_id": patient_id, "rna_cancer_type": mode.iloc[0] if not mode.empty else str(group["cancer_type"].iloc[0])})
    return pd.DataFrame(rows)


def _safe_numeric_ids(path: Path) -> set[str]:
    return _patient_ids_from_index_or_column(_read_parquet(path))


def build_tumor_state_registry(config_path: str = "morpheus/configs/v1.json") -> dict[str, Path]:
    cfg = load_config(config_path)
    processed = cfg.path("processed_dir")
    processed.mkdir(parents=True, exist_ok=True)

    survival_path = cfg.project_root / "meta-intersurv" / "data" / "embeddings" / "patient_master" / "survival_labels.parquet"
    survival = _read_parquet(survival_path).rename(columns={"time": "survival_time", "event": "survival_event"})
    if not survival.empty:
        survival["patient_id"] = survival["patient_id"].map(normalize_patient_id)
        survival = survival.dropna(subset=["patient_id"]).drop_duplicates("patient_id")

    rna_path = cfg.path("rna_bulkformer_embeddings")
    rna_ids = _safe_numeric_ids(rna_path)
    rna_meta = _bulkformer_metadata(rna_path)
    wsi_path = cfg.path("wsi_standard_dir") / "tcga_ut_hoptimus0_patient_embeddings.h5"
    wsi_ids = _wsi_ids(wsi_path)
    hallmark_ids = _safe_numeric_ids(cfg.path("hallmark_scores"))
    clinical_path = cfg.project_root / "meta-intersurv" / "data" / "embeddings" / "patient_master" / "tabular_features.parquet"
    clinical_ids = _safe_numeric_ids(clinical_path)
    snv_ids = _safe_numeric_ids(cfg.path("snv_processed"))
    cnv_ids = _safe_numeric_ids(cfg.path("cnv_processed"))

    all_ids = sorted(set(survival.get("patient_id", pd.Series(dtype=str)).astype(str)) | rna_ids | wsi_ids | hallmark_ids | clinical_ids | snv_ids | cnv_ids)
    registry = pd.DataFrame({"patient_id": all_ids})
    registry = registry.merge(survival[["patient_id", "cancer_type", "survival_time", "survival_event"]] if not survival.empty else pd.DataFrame(columns=["patient_id", "cancer_type", "survival_time", "survival_event"]), on="patient_id", how="left")
    registry = registry.merge(rna_meta, on="patient_id", how="left")
    registry["cancer_type"] = registry["cancer_type"].fillna(registry["rna_cancer_type"])
    registry = registry.drop(columns=["rna_cancer_type"])
    registry["has_wsi_hoptimus0"] = registry["patient_id"].isin(wsi_ids)
    registry["has_rna_bulkformer"] = registry["patient_id"].isin(rna_ids)
    registry["has_hallmark_scores"] = registry["patient_id"].isin(hallmark_ids)
    registry["has_clinical_tabular"] = registry["patient_id"].isin(clinical_ids)
    registry["has_snv"] = registry["patient_id"].isin(snv_ids)
    registry["has_cnv"] = registry["patient_id"].isin(cnv_ids)
    registry["has_survival"] = registry["survival_time"].notna() & registry["survival_event"].notna() & (registry["survival_time"] > 0)
    registry["modality_count"] = registry[["has_wsi_hoptimus0", "has_rna_bulkformer", "has_hallmark_scores", "has_clinical_tabular", "has_snv", "has_cnv"]].sum(axis=1)

    out_path = processed / "patient_feature_registry.parquet"
    registry.to_parquet(out_path, index=False)
    counts = {
        "patients": int(len(registry)),
        "disease_labeled": int(registry["cancer_type"].notna().sum()),
        "disease_types": int(registry["cancer_type"].dropna().nunique()),
        "wsi_rna": int((registry["has_wsi_hoptimus0"] & registry["has_rna_bulkformer"]).sum()),
        "wsi_rna_hallmark": int((registry["has_wsi_hoptimus0"] & registry["has_rna_bulkformer"] & registry["has_hallmark_scores"]).sum()),
        "wsi_rna_clinical": int((registry["has_wsi_hoptimus0"] & registry["has_rna_bulkformer"] & registry["has_clinical_tabular"]).sum()),
    }
    manifest = base_manifest(cfg.project_root, cfg.config_path, int(cfg.raw.get("seed", 42)))
    manifest.update({"output": str(out_path), "counts": counts, "sources": {"survival": str(survival_path), "rna": str(rna_path), "wsi": str(wsi_path), "hallmark": str(cfg.path("hallmark_scores")), "clinical": str(clinical_path)}})
    write_json(processed / "patient_feature_registry.manifest.json", manifest)
    return {"registry": out_path, "manifest": processed / "patient_feature_registry.manifest.json"}


def _split_ids(ids: list[str], rng: np.random.Generator, val_fraction: float, test_fraction: float) -> dict[str, list[str]]:
    ids = ids[:]
    rng.shuffle(ids)
    n = len(ids)
    n_test = int(round(n * test_fraction))
    n_val = int(round(n * val_fraction))
    return {
        "train": sorted(ids[: max(0, n - n_val - n_test)]),
        "val": sorted(ids[max(0, n - n_val - n_test) : max(0, n - n_test)]),
        "test": sorted(ids[max(0, n - n_test) :]),
    }


def _payload(cfg, name: str, patient_ids: dict[str, list[str]], extra: dict[str, Any] | None = None) -> dict[str, Any]:
    payload = base_manifest(cfg.project_root, cfg.config_path, int(cfg.raw.get("seed", 42)))
    payload.update({"split_name": name, "patient_ids": patient_ids})
    if extra:
        payload.update(extra)
    return payload


def create_tumor_state_splits(config_path: str = "morpheus/configs/v1.json") -> dict[str, Path]:
    cfg = load_config(config_path)
    registry_path = cfg.path("processed_dir") / "patient_feature_registry.parquet"
    if not registry_path.exists():
        build_tumor_state_registry(config_path)
    registry = pd.read_parquet(registry_path)
    eligible = registry[registry["cancer_type"].notna() & (registry["has_wsi_hoptimus0"] | registry["has_rna_bulkformer"] | registry["has_clinical_tabular"])].copy()
    seed = int(cfg.raw.get("seed", 42))
    rng = np.random.default_rng(seed)
    out_dir = cfg.path("processed_dir") / "splits"
    out_dir.mkdir(parents=True, exist_ok=True)

    stratified = {"train": [], "val": [], "test": []}
    split_cfg = cfg.section("splits")
    for _, group in eligible.groupby("cancer_type", dropna=True):
        part = _split_ids(sorted(group["patient_id"].astype(str).unique()), rng, float(split_cfg.get("val_fraction", 0.15)), float(split_cfg.get("test_fraction", 0.15)))
        for key in stratified:
            stratified[key].extend(part[key])
    stratified = {k: sorted(v) for k, v in stratified.items()}

    counts = eligible["cancer_type"].value_counts()
    train_cancers = counts.head(11).index.astype(str).tolist()
    test_cancers = counts.loc[~counts.index.isin(train_cancers)].index.astype(str).tolist()
    heldout = {
        "train": sorted(eligible.loc[eligible["cancer_type"].isin(train_cancers), "patient_id"].astype(str).tolist()),
        "val": [],
        "test": sorted(eligible.loc[eligible["cancer_type"].isin(test_cancers), "patient_id"].astype(str).tolist()),
    }
    heldout_part = _split_ids(heldout["train"], rng, float(split_cfg.get("val_fraction", 0.15)), 0.0)
    heldout["train"], heldout["val"] = heldout_part["train"], heldout_part["val"]

    masks = {
        "wsi_only": ["wsi"],
        "rna_only": ["rna"],
        "clinical_only": ["clinical"],
        "wsi_rna": ["wsi", "rna"],
        "wsi_clinical": ["wsi", "clinical"],
        "rna_clinical": ["rna", "clinical"],
        "wsi_rna_clinical": ["wsi", "rna", "clinical"],
        "full_available": ["wsi", "rna", "clinical", "snv", "cnv", "hallmark"],
    }

    outputs = {
        "stratified": out_dir / "tumor_state_stratified.json",
        "heldout_cancer": out_dir / "tumor_state_heldout_cancer.json",
        "missingness": out_dir / "tumor_state_missingness.json",
    }
    write_json(outputs["stratified"], _payload(cfg, "tumor_state_stratified", stratified, {"eligible_n": int(len(eligible))}))
    write_json(outputs["heldout_cancer"], _payload(cfg, "tumor_state_heldout_cancer", heldout, {"train_cancers": train_cancers, "test_cancers": test_cancers, "actual_train_test_cancer_counts": [len(train_cancers), len(test_cancers)]}))
    write_json(outputs["missingness"], {**base_manifest(cfg.project_root, cfg.config_path, seed), "split_name": "tumor_state_missingness", "masks": masks})
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    args = parser.parse_args()
    outputs = build_tumor_state_registry(args.config)
    outputs.update(create_tumor_state_splits(args.config))
    print({k: str(v) for k, v in outputs.items()})


if __name__ == "__main__":
    main()
