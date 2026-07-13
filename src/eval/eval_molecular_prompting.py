"""Soft-kNN molecular prompting evaluation from aligned WSI-RNA embeddings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from morpheus.src.eval.retrieval_metrics import l2_normalize
from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import base_manifest, write_json


def soft_knn_predict(query: np.ndarray, reference: np.ndarray, targets: np.ndarray, k: int = 5, tau: float = 0.1) -> np.ndarray:
    q = l2_normalize(query)
    r = l2_normalize(reference)
    scores = q @ r.T
    top = np.argsort(-scores, axis=1)[:, :k]
    preds = []
    for i, idx in enumerate(top):
        weights = np.exp(scores[i, idx] / tau)
        weights = weights / max(float(weights.sum()), 1e-12)
        preds.append(weights @ targets[idx])
    return np.vstack(preds).astype(np.float32)


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray, names: list[str]) -> tuple[pd.DataFrame, dict]:
    rows = []
    for j, name in enumerate(names):
        truth = y_true[:, j]
        pred = y_pred[:, j]
        ss_res = float(np.sum((truth - pred) ** 2))
        ss_tot = float(np.sum((truth - truth.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        pearson = float(np.corrcoef(truth, pred)[0, 1]) if len(truth) > 1 and np.std(truth) > 0 and np.std(pred) > 0 else 0.0
        spearman = float(spearmanr(truth, pred).correlation) if len(truth) > 1 and np.std(truth) > 0 and np.std(pred) > 0 else 0.0
        rows.append({"gene_set": name, "r2": r2, "pearson": pearson, "spearman": 0.0 if np.isnan(spearman) else spearman})
    df = pd.DataFrame(rows)
    return df, {"mean_r2": float(df["r2"].mean()), "mean_pearson": float(df["pearson"].mean()), "mean_spearman": float(df["spearman"].mean()), "n_gene_sets": int(len(df))}


def _load_targets(cfg, patient_ids: np.ndarray, rna: np.ndarray) -> tuple[np.ndarray, list[str], str]:
    score_path = cfg.path("hallmark_scores")
    if score_path.exists():
        scores = pd.read_parquet(score_path)
        if "patient_id" not in scores:
            raise ValueError(f"Gene-set score table must include patient_id: {score_path}")
        scores = scores.set_index("patient_id")
        names = [c for c in scores.columns if pd.api.types.is_numeric_dtype(scores[c])]
        aligned = scores.reindex(patient_ids)
        mask = ~aligned[names].isna().any(axis=1).to_numpy()
        targets = aligned[names].to_numpy(dtype=np.float32)
        return targets, names, "hallmark_gene_set_scores"
    target_names = [f"rna_dim_{i}" for i in range(min(64, rna.shape[1]))]
    return rna[:, : len(target_names)], target_names, "rna_embedding_dimensions_fallback"


def _empty_prompt_metrics(n_eval: int, n_train: int, target_source: str) -> dict:
    return {
        "n_eval": int(n_eval),
        "n_train_reference": int(n_train),
        "mean_r2": None,
        "mean_pearson": None,
        "mean_spearman": None,
        "n_gene_sets": 0,
        "target_source": target_source,
    }


def _evaluate_prompt_split(
    split_name: str,
    query_mask: np.ndarray,
    train_mask: np.ndarray,
    patient_ids: np.ndarray,
    wsi: np.ndarray,
    rna: np.ndarray,
    targets: np.ndarray,
    target_names: list[str],
    target_source: str,
    k: int,
    tau: float,
    out_dir: Path,
) -> dict:
    if not train_mask.any() or not query_mask.any():
        return _empty_prompt_metrics(int(query_mask.sum()), int(train_mask.sum()), target_source)
    pred = soft_knn_predict(wsi[query_mask], rna[train_mask], targets[train_mask], k, tau)
    per, global_metrics = regression_metrics(targets[query_mask], pred, target_names)
    per.to_csv(out_dir / f"{split_name}_per_geneset_metrics.csv", index=False)
    pd.DataFrame(pred, columns=target_names).assign(patient_id=patient_ids[query_mask]).to_parquet(out_dir / f"{split_name}_predicted_geneset_scores.parquet", index=False)
    return {
        **global_metrics,
        "n_eval": int(query_mask.sum()),
        "n_train_reference": int(train_mask.sum()),
        "target_source": target_source,
    }


def evaluate_molecular_prompting(config_path: str, aligned_npz: str | Path, output_dir: str | Path | None = None) -> dict[str, Path]:
    cfg = load_config(config_path)
    data = np.load(aligned_npz, allow_pickle=True)
    patient_ids = data["patient_ids"].astype(str)
    wsi = data["wsi"].astype(np.float32)
    rna = data["rna"].astype(np.float32)
    aligned_split = data["split"].astype(str) if "split" in data.files else None
    split_path = cfg.path("processed_dir") / "splits" / "stratified_pan_cancer_split.json"
    if not split_path.exists():
        raise FileNotFoundError("Create splits before molecular prompting evaluation")
    split = json.load(open(split_path, encoding="utf-8"))["patient_ids"]
    train_mask_raw = aligned_split == "train" if aligned_split is not None else np.isin(patient_ids, np.asarray(split["train"]))
    val_mask_raw = aligned_split == "val" if aligned_split is not None else np.isin(patient_ids, np.asarray(split.get("val", [])))
    test_mask_raw = aligned_split == "test" if aligned_split is not None else np.isin(patient_ids, np.asarray(split["test"]))
    if not train_mask_raw.any() or not (val_mask_raw.any() or test_mask_raw.any()):
        raise ValueError("Need train and val/test aligned patients for molecular prompting")
    targets, target_names, target_source = _load_targets(cfg, patient_ids, rna)
    valid_target_mask = ~np.isnan(targets).any(axis=1)
    train_mask = train_mask_raw & valid_target_mask
    val_mask = val_mask_raw & valid_target_mask
    test_mask = test_mask_raw & valid_target_mask
    mp = cfg.section("molecular_prompting")
    out_dir = Path(output_dir) if output_dir else cfg.path("outputs_dir") / "v1_molecular_prompting"
    out_dir.mkdir(parents=True, exist_ok=True)
    global_path = out_dir / "global_metrics.json"
    k = int(mp.get("k", 5))
    tau = float(mp.get("tau", 0.1))
    val_metrics = _evaluate_prompt_split("val", val_mask, train_mask, patient_ids, wsi, rna, targets, target_names, target_source, k, tau, out_dir)
    test_metrics = _evaluate_prompt_split("test", test_mask, train_mask, patient_ids, wsi, rna, targets, target_names, target_source, k, tau, out_dir)
    payload = base_manifest(cfg.project_root, cfg.config_path, int(cfg.raw.get("seed", 42)))
    payload.update(test_metrics)
    payload["val_metrics"] = val_metrics
    payload["test_metrics"] = test_metrics
    payload["aligned_npz"] = str(aligned_npz)
    payload["target_source"] = target_source
    payload["claim_level"] = "retrieval_based_molecular_prompting" if target_source == "hallmark_gene_set_scores" else "retrieval_based_embedding_score_recovery"
    payload["note"] = "Soft-kNN over aligned RNA training references; target table is never fit on test predictions."
    write_json(global_path, payload)
    return {"global": global_path}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--aligned-npz", required=True)
    parser.add_argument("--output-dir")
    args = parser.parse_args()
    print(evaluate_molecular_prompting(args.config, args.aligned_npz, args.output_dir))


if __name__ == "__main__":
    main()
