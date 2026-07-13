"""Create locked patient-level split files for v1."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import base_manifest, write_json


def _split_ids(ids: list[str], rng: np.random.Generator, val_fraction: float, test_fraction: float) -> dict[str, list[str]]:
    shuffled = ids[:]
    rng.shuffle(shuffled)
    n = len(shuffled)
    n_test = int(round(n * test_fraction))
    n_val = int(round(n * val_fraction))
    return {
        "train": sorted(shuffled[: max(0, n - n_val - n_test)]),
        "val": sorted(shuffled[max(0, n - n_val - n_test) : max(0, n - n_test)]),
        "test": sorted(shuffled[max(0, n - n_test) :]),
    }


def create_splits(config_path: str = "morpheus/configs/v1.json") -> dict[str, Path]:
    cfg = load_config(config_path)
    split_cfg = cfg.section("splits")
    seed = int(cfg.raw.get("seed", 42))
    rng = np.random.default_rng(seed)
    master_path = cfg.path("processed_dir") / "master_patient_table.parquet"
    if not master_path.exists():
        raise FileNotFoundError(f"Build master table first: {master_path}")
    master = pd.read_parquet(master_path)
    eligible = master[master["has_clinical"] | master["has_rna"] | master["has_wsi"]].copy()
    ids = sorted(eligible["patient_id"].astype(str).unique())
    out_dir = cfg.path("processed_dir") / "splits"
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs: dict[str, Path] = {}
    random_split = {
        **base_manifest(cfg.project_root, cfg.config_path, seed),
        "split_name": "random_patient_split",
        "patient_ids": _split_ids(ids, rng, float(split_cfg.get("val_fraction", 0.15)), float(split_cfg.get("test_fraction", 0.15))),
    }
    outputs["random"] = out_dir / "random_patient_split.json"
    write_json(outputs["random"], random_split)

    stratified = {"train": [], "val": [], "test": []}
    for _, group in eligible.groupby("cancer_type", dropna=False):
        group_ids = sorted(group["patient_id"].astype(str).unique())
        part = _split_ids(group_ids, rng, float(split_cfg.get("val_fraction", 0.15)), float(split_cfg.get("test_fraction", 0.15)))
        for key in stratified:
            stratified[key].extend(part[key])
    strat_payload = {
        **base_manifest(cfg.project_root, cfg.config_path, seed),
        "split_name": "stratified_pan_cancer_split",
        "patient_ids": {k: sorted(v) for k, v in stratified.items()},
    }
    outputs["stratified"] = out_dir / "stratified_pan_cancer_split.json"
    write_json(outputs["stratified"], strat_payload)

    heldout: dict[str, dict[str, list[str]]] = {}
    for cancer in split_cfg.get("heldout_cancers", []):
        test_ids = sorted(eligible.loc[eligible["cancer_type"].astype(str) == cancer, "patient_id"].astype(str).unique())
        train_val = sorted(set(ids) - set(test_ids))
        part = _split_ids(train_val, rng, float(split_cfg.get("val_fraction", 0.15)), 0.0)
        heldout[cancer] = {"train": part["train"], "val": part["val"], "test": test_ids}
    outputs["heldout"] = out_dir / "heldout_cancer_splits.json"
    write_json(outputs["heldout"], {**base_manifest(cfg.project_root, cfg.config_path, seed), "split_name": "heldout_cancer_splits", "splits": heldout})

    adaptation = {}
    fractions = [float(x) for x in split_cfg.get("adaptation_fractions", [0.05, 0.1, 0.25, 0.5, 1.0])]
    for cancer in split_cfg.get("domain_adaptation_cancers", []):
        target = sorted(eligible.loc[eligible["cancer_type"].astype(str) == cancer, "patient_id"].astype(str).unique())
        rng.shuffle(target)
        adaptation[cancer] = {str(frac): sorted(target[: max(1, int(round(len(target) * frac)))]) for frac in fractions if target}
    outputs["domain_adaptation"] = out_dir / "domain_adaptation_brca_luad.json"
    write_json(outputs["domain_adaptation"], {**base_manifest(cfg.project_root, cfg.config_path, seed), "split_name": "domain_adaptation_brca_luad", "adaptation_subsets": adaptation})

    masks = {
        "clinical_only": ["clinical"],
        "wsi_only": ["wsi"],
        "rna_only": ["rna"],
        "clinical_wsi": ["clinical", "wsi"],
        "clinical_rna": ["clinical", "rna"],
        "wsi_rna": ["wsi", "rna"],
        "full": ["clinical", "wsi", "rna", "snv", "cnv"],
    }
    outputs["missingness"] = out_dir / "missingness_eval_masks.json"
    write_json(outputs["missingness"], {**base_manifest(cfg.project_root, cfg.config_path, seed), "split_name": "missingness_eval_masks", "masks": masks})
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    args = parser.parse_args()
    print(json.dumps({k: str(v) for k, v in create_splits(args.config).items()}, indent=2))


if __name__ == "__main__":
    main()
