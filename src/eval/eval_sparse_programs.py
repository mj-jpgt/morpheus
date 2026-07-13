"""Sparse biological program discovery baselines for V2."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.src.eval.eval_molecular_prompting import _load_targets
from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import base_manifest, write_json


def _load_embedding_matrix(path: str | Path, feature_key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.load(path, allow_pickle=True)
    if feature_key not in data.files:
        raise KeyError(f"{feature_key!r} not found in {path}; available keys: {data.files}")
    return data["patient_ids"].astype(str), data["split"].astype(str), data[feature_key].astype(np.float32)


def _standardize_train(x: np.ndarray, train: np.ndarray) -> np.ndarray:
    mean = x[train].mean(axis=0, keepdims=True)
    std = x[train].std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return ((x - mean) / std).astype(np.float32)


def _fit_programs(method: str, x: np.ndarray, train: np.ndarray, n_programs: int, seed: int) -> np.ndarray:
    if method == "pca":
        from sklearn.decomposition import PCA

        model = PCA(n_components=n_programs, random_state=seed)
        model.fit(x[train])
        return model.transform(x).astype(np.float32)
    if method == "ica":
        from sklearn.decomposition import FastICA

        model = FastICA(n_components=n_programs, random_state=seed, max_iter=1000, whiten="unit-variance")
        model.fit(x[train])
        return model.transform(x).astype(np.float32)
    if method == "nmf":
        from sklearn.decomposition import NMF
        from sklearn.preprocessing import MinMaxScaler

        scaler = MinMaxScaler()
        x_pos = scaler.fit_transform(x)
        model = NMF(n_components=n_programs, random_state=seed, max_iter=1000, init="nndsvda")
        model.fit(x_pos[train])
        return model.transform(x_pos).astype(np.float32)
    if method == "kmeans":
        from sklearn.cluster import KMeans

        model = KMeans(n_clusters=n_programs, random_state=seed, n_init=10)
        model.fit(x[train])
        labels = model.predict(x)
        return np.eye(n_programs, dtype=np.float32)[labels]
    if method == "random_sparse":
        rng = np.random.default_rng(seed)
        weights = rng.normal(size=(x.shape[1], n_programs)).astype(np.float32)
        mask = rng.random(weights.shape) < 0.1
        weights = weights * mask
        return (x @ weights).astype(np.float32)
    raise ValueError(f"Unknown sparse program method: {method}")


def _mean_abs_corr(a: np.ndarray, b: np.ndarray) -> float:
    vals = []
    for i in range(a.shape[1]):
        for j in range(b.shape[1]):
            if np.std(a[:, i]) > 0 and np.std(b[:, j]) > 0:
                vals.append(abs(float(np.corrcoef(a[:, i], b[:, j])[0, 1])))
    return float(np.mean(vals)) if vals else 0.0


def _hallmark_association(programs: np.ndarray, targets: np.ndarray, valid: np.ndarray) -> dict:
    if not valid.any():
        return {"mean_abs_hallmark_corr": None, "max_abs_hallmark_corr": None}
    vals = []
    for i in range(programs.shape[1]):
        for j in range(targets.shape[1]):
            p = programs[valid, i]
            t = targets[valid, j]
            if np.std(p) > 0 and np.std(t) > 0:
                vals.append(abs(float(np.corrcoef(p, t)[0, 1])))
    return {
        "mean_abs_hallmark_corr": float(np.mean(vals)) if vals else None,
        "max_abs_hallmark_corr": float(np.max(vals)) if vals else None,
    }


def _cancer_leakage(programs: np.ndarray, cancers: np.ndarray, split: np.ndarray) -> dict:
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.preprocessing import LabelEncoder, StandardScaler
    from sklearn.pipeline import make_pipeline

    train = split == "train"
    test = split == "test"
    if train.sum() < 10 or test.sum() < 10 or len(np.unique(cancers[train])) < 2:
        return {"cancer_balanced_accuracy": None}
    enc = LabelEncoder()
    y_train = enc.fit_transform(cancers[train])
    known = np.isin(cancers[test], enc.classes_)
    if known.sum() < 10:
        return {"cancer_balanced_accuracy": None}
    y_test = enc.transform(cancers[test][known])
    model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, class_weight="balanced"))
    model.fit(programs[train], y_train)
    pred = model.predict(programs[test][known])
    return {"cancer_balanced_accuracy": float(balanced_accuracy_score(y_test, pred))}


def evaluate_sparse_programs(
    config_path: str,
    embedding_npz: str | Path,
    output_dir: str | Path,
    feature_key: str = "wsi_biology",
    methods: tuple[str, ...] = ("pca", "nmf", "ica", "kmeans", "random_sparse"),
    n_programs: int = 16,
    seeds: tuple[int, ...] = (1, 2, 3),
) -> Path:
    cfg = load_config(config_path)
    patient_ids, split, x = _load_embedding_matrix(embedding_npz, feature_key)
    data = np.load(embedding_npz, allow_pickle=True)
    cancers = data["cancers"].astype(str) if "cancers" in data.files else np.asarray(["unknown"] * len(patient_ids))
    train = split == "train"
    x = _standardize_train(x, train)
    targets, target_names, target_source = _load_targets(cfg, patient_ids, x)
    valid_targets = ~np.isnan(targets).any(axis=1)
    rows = []
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for method in methods:
        seed_programs = []
        for seed in seeds:
            programs = _fit_programs(method, x, train, n_programs, seed)
            seed_programs.append(programs)
        programs = seed_programs[0]
        pd.DataFrame(programs, columns=[f"program_{i:02d}" for i in range(programs.shape[1])]).assign(patient_id=patient_ids, split=split, cancer=cancers).to_parquet(out / f"{method}_program_scores.parquet", index=False)
        stability = float(np.mean([_mean_abs_corr(seed_programs[0], other) for other in seed_programs[1:]])) if len(seed_programs) > 1 else 1.0
        sparsity = float(np.mean(np.abs(programs) < np.percentile(np.abs(programs), 25)))
        hallmark = _hallmark_association(programs, targets, valid_targets)
        leakage = _cancer_leakage(programs, cancers, split)
        rows.append(
            {
                "method": method,
                "feature_key": feature_key,
                "n_programs": int(programs.shape[1]),
                "stability_mean_abs_corr": stability,
                "activation_sparsity_q25": sparsity,
                **hallmark,
                **leakage,
                "target_source": target_source,
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "sparse_program_summary.csv", index=False)
    payload = base_manifest(cfg.project_root, cfg.config_path, int(cfg.raw.get("seed", 42)))
    payload.update(
        {
            "embedding_npz": str(embedding_npz),
            "feature_key": feature_key,
            "n_patients": int(len(patient_ids)),
            "methods": list(methods),
            "n_programs": int(n_programs),
            "summary_csv": str(out / "sparse_program_summary.csv"),
            "claim_level": "computational_program_characterization",
            "note": "Programs are not discoveries unless stable and independently grounded.",
        }
    )
    manifest_path = out / "sparse_program_manifest.json"
    write_json(manifest_path, payload)
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--embedding-npz", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--feature-key", default="wsi_biology")
    parser.add_argument("--n-programs", type=int, default=16)
    args = parser.parse_args()
    print(evaluate_sparse_programs(args.config, args.embedding_npz, args.output_dir, args.feature_key, n_programs=args.n_programs))


if __name__ == "__main__":
    main()
