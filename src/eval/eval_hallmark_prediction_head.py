"""Train/evaluate a common Hallmark prediction head on aligned embeddings."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.src.eval.eval_molecular_prompting import _load_targets, regression_metrics
from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import base_manifest, write_json


def _empty_head_metrics(n_eval: int, n_train: int, target_source: str) -> dict:
    return {
        "n_eval": int(n_eval),
        "n_train": int(n_train),
        "mean_r2": None,
        "mean_pearson": None,
        "mean_spearman": None,
        "n_gene_sets": 0,
        "target_source": target_source,
    }


def _fit_predict_ridge(train_x: np.ndarray, train_y: np.ndarray, eval_x: np.ndarray) -> np.ndarray:
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    model = make_pipeline(
        StandardScaler(),
        RidgeCV(alphas=np.asarray([0.01, 0.1, 1.0, 10.0, 100.0], dtype=np.float32)),
    )
    model.fit(train_x, train_y)
    return model.predict(eval_x).astype(np.float32)


def _evaluate_split(
    split_name: str,
    x: np.ndarray,
    targets: np.ndarray,
    target_names: list[str],
    target_source: str,
    train_mask: np.ndarray,
    eval_mask: np.ndarray,
    out_dir: Path,
) -> dict:
    if not train_mask.any() or not eval_mask.any():
        return _empty_head_metrics(int(eval_mask.sum()), int(train_mask.sum()), target_source)
    pred = _fit_predict_ridge(x[train_mask], targets[train_mask], x[eval_mask])
    per, metrics = regression_metrics(targets[eval_mask], pred, target_names)
    per.to_csv(out_dir / f"{split_name}_head_per_geneset_metrics.csv", index=False)
    pd.DataFrame(pred, columns=target_names).to_parquet(out_dir / f"{split_name}_head_predicted_geneset_scores.parquet", index=False)
    return {
        **metrics,
        "n_eval": int(eval_mask.sum()),
        "n_train": int(train_mask.sum()),
        "target_source": target_source,
    }


def evaluate_hallmark_prediction_head(
    config_path: str,
    aligned_npz: str | Path,
    output_dir: str | Path,
    feature_key: str = "wsi",
) -> Path:
    cfg = load_config(config_path)
    data = np.load(aligned_npz, allow_pickle=True)
    if feature_key not in data.files:
        raise KeyError(f"{feature_key!r} not found in {aligned_npz}; available keys: {data.files}")
    patient_ids = data["patient_ids"].astype(str)
    split = data["split"].astype(str)
    x = data[feature_key].astype(np.float32)
    targets, target_names, target_source = _load_targets(cfg, patient_ids, x)
    valid_target_mask = ~np.isnan(targets).any(axis=1)
    train_mask = (split == "train") & valid_target_mask
    val_mask = (split == "val") & valid_target_mask
    test_mask = (split == "test") & valid_target_mask
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    val_metrics = _evaluate_split("val", x, targets, target_names, target_source, train_mask, val_mask, out_dir)
    test_metrics = _evaluate_split("test", x, targets, target_names, target_source, train_mask, test_mask, out_dir)
    payload = base_manifest(cfg.project_root, cfg.config_path, int(cfg.raw.get("seed", 42)))
    payload.update(test_metrics)
    payload.update(
        {
            "aligned_npz": str(aligned_npz),
            "feature_key": feature_key,
            "head": "StandardScaler + RidgeCV",
            "val_metrics": val_metrics,
            "test_metrics": test_metrics,
            "claim_level": "common_supervised_hallmark_prediction_head",
            "note": "The same train-only ridge head is fit on each method's selected embedding feature key.",
        }
    )
    out_path = out_dir / "head_global_metrics.json"
    write_json(out_path, payload)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--aligned-npz", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--feature-key", default="wsi")
    args = parser.parse_args()
    print(evaluate_hallmark_prediction_head(args.config, args.aligned_npz, args.output_dir, args.feature_key))


if __name__ == "__main__":
    main()
