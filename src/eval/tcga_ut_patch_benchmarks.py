"""Leakage-safe internal/external TCGA-UT benchmark helpers.

The patch store has two valid, non-interchangeable split schemes.  This module
keeps their choice explicit and evaluates an already-produced aligned embedding
archive.  In particular, molecular prompting *always* uses RNA references from
the selected protocol's training patients only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from morpheus.src.eval.eval_molecular_prompting import regression_metrics, soft_knn_predict
from morpheus.src.eval.retrieval_metrics import paired_retrieval_metrics
from morpheus.src.utils.ids import normalize_patient_id
from morpheus.src.utils.provenance import write_json


VALID_PROTOCOLS = ("internal", "external")
VALID_SPLITS = ("train", "val", "test")


def _scalar_string(value: Any) -> str:
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError("split_protocol must be a scalar value")
    return str(array.reshape(-1)[0]).strip().lower()


def protocol_splits_from_npz(data: np.lib.npyio.NpzFile, protocol: str) -> np.ndarray:
    """Return the requested split membership without silently reusing another protocol."""
    protocol = protocol.lower()
    if protocol not in VALID_PROTOCOLS:
        raise ValueError(f"protocol must be one of {VALID_PROTOCOLS}, got {protocol!r}")
    explicit_key = f"{protocol}_split"
    if explicit_key in data.files:
        split = data[explicit_key]
    elif "split" in data.files and "split_protocol" in data.files:
        declared = _scalar_string(data["split_protocol"])
        if declared != protocol:
            raise ValueError(f"Archive declares split_protocol={declared!r}, not {protocol!r}")
        split = data["split"]
    else:
        raise KeyError(
            f"Archive must contain {explicit_key!r}, or both 'split' and 'split_protocol'; "
            "refusing an ambiguous split assignment."
        )
    split = np.asarray(split).astype(str)
    unknown = sorted(set(split) - set(VALID_SPLITS))
    if unknown:
        raise ValueError(f"Unexpected split labels: {unknown}; expected {VALID_SPLITS}")
    return split


def _canonical_patient_ids(values: np.ndarray) -> np.ndarray:
    result = np.asarray([normalize_patient_id(v) for v in values], dtype=object)
    if any(v is None for v in result):
        raise ValueError("Every benchmark patient_id must be non-empty")
    result = result.astype(str)
    duplicates = pd.Series(result).duplicated(keep=False)
    if duplicates.any():
        repeated = sorted(pd.Series(result)[duplicates].unique())[:5]
        raise ValueError(f"Benchmark archive must have one row per patient; duplicates include {repeated}")
    return result


def load_aligned_protocol_archive(path: str | Path, protocol: str) -> dict[str, np.ndarray]:
    """Load a paired WSI/RNA archive and verify its protocol-specific membership."""
    with np.load(path, allow_pickle=False) as data:
        required = {"patient_ids", "wsi", "rna"}
        missing = sorted(required - set(data.files))
        if missing:
            raise KeyError(f"Aligned archive missing keys: {missing}")
        patient_ids = _canonical_patient_ids(data["patient_ids"].astype(str))
        wsi = np.asarray(data["wsi"], dtype=np.float32)
        rna = np.asarray(data["rna"], dtype=np.float32)
        split = protocol_splits_from_npz(data, protocol)
    if wsi.ndim != 2 or rna.ndim != 2 or len(wsi) != len(patient_ids) or len(rna) != len(patient_ids):
        raise ValueError("patient_ids, wsi, and rna must have the same first dimension; embeddings must be matrices")
    if not np.isfinite(wsi).all() or not np.isfinite(rna).all():
        raise ValueError("Aligned archive contains NaN or Inf embeddings")
    return {"patient_ids": patient_ids, "wsi": wsi, "rna": rna, "split": split}


def load_targets(path: str | Path, patient_ids: np.ndarray) -> tuple[np.ndarray, list[str]]:
    """Join a patient-level numeric target table with strict one-row cardinality."""
    source = Path(path)
    targets = pd.read_parquet(source) if source.suffix.lower() in {".parquet", ".pq"} else pd.read_csv(source)
    if "patient_id" not in targets:
        raise ValueError(f"Target table must include patient_id: {source}")
    targets = targets.copy()
    targets["patient_id"] = _canonical_patient_ids(targets["patient_id"].to_numpy())
    targets = targets.set_index("patient_id")
    names = [column for column in targets if pd.api.types.is_numeric_dtype(targets[column])]
    if not names:
        raise ValueError("Target table has no numeric target columns")
    joined = targets.reindex(patient_ids)[names]
    return joined.to_numpy(dtype=np.float32), names


def _fit_predict_ridge(train_x: np.ndarray, train_y: np.ndarray, query_x: np.ndarray) -> np.ndarray:
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    model = make_pipeline(StandardScaler(), RidgeCV(alphas=np.asarray([0.01, 0.1, 1.0, 10.0, 100.0])))
    model.fit(train_x, train_y)
    prediction = model.predict(query_x).astype(np.float32)
    return prediction[:, None] if prediction.ndim == 1 else prediction


def _evaluate_target_split(
    name: str,
    query_mask: np.ndarray,
    train_mask: np.ndarray,
    patient_ids: np.ndarray,
    wsi: np.ndarray,
    rna: np.ndarray,
    targets: np.ndarray,
    target_names: list[str],
    output_dir: Path,
    k: int,
    tau: float,
) -> dict[str, Any]:
    if np.any(query_mask & train_mask):
        raise AssertionError("Prompting query patients overlap train RNA references")
    valid = np.isfinite(targets).all(axis=1)
    query_mask = query_mask & valid
    train_mask = train_mask & valid
    if not query_mask.any() or not train_mask.any():
        return {"n_eval": int(query_mask.sum()), "n_train_reference": int(train_mask.sum()), "soft_knn": None, "ridge": None}
    soft_pred = soft_knn_predict(wsi[query_mask], rna[train_mask], targets[train_mask], k=k, tau=tau)
    ridge_pred = _fit_predict_ridge(wsi[train_mask], targets[train_mask], wsi[query_mask])
    soft_per, soft_metrics = regression_metrics(targets[query_mask], soft_pred, target_names)
    ridge_per, ridge_metrics = regression_metrics(targets[query_mask], ridge_pred, target_names)
    soft_per.to_csv(output_dir / f"{name}_soft_knn_per_target.csv", index=False)
    ridge_per.to_csv(output_dir / f"{name}_ridge_per_target.csv", index=False)
    prediction_ids = patient_ids[query_mask]
    pd.DataFrame(soft_pred, columns=target_names).assign(patient_id=prediction_ids).to_parquet(output_dir / f"{name}_soft_knn_predictions.parquet", index=False)
    pd.DataFrame(ridge_pred, columns=target_names).assign(patient_id=prediction_ids).to_parquet(output_dir / f"{name}_ridge_predictions.parquet", index=False)
    return {
        "n_eval": int(query_mask.sum()),
        "n_train_reference": int(train_mask.sum()),
        "soft_knn": soft_metrics,
        "ridge": ridge_metrics,
    }


def evaluate_protocol_benchmark(
    aligned_npz: str | Path,
    targets_path: str | Path,
    protocol: str,
    output_dir: str | Path,
    k: int = 5,
    tau: float = 0.1,
) -> Path:
    """Evaluate paired retrieval plus train-only soft-kNN/ridge molecular prompting."""
    if k < 1 or tau <= 0:
        raise ValueError("k must be >= 1 and tau must be positive")
    archive = load_aligned_protocol_archive(aligned_npz, protocol)
    patient_ids, wsi, rna, split = (archive[key] for key in ("patient_ids", "wsi", "rna", "split"))
    targets, target_names = load_targets(targets_path, patient_ids)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    train_mask = split == "train"
    payload: dict[str, Any] = {
        "protocol": protocol,
        "aligned_npz": str(aligned_npz),
        "targets_path": str(targets_path),
        "n_patients": int(len(patient_ids)),
        "embedding_dimensions": {"wsi": int(wsi.shape[1]), "rna": int(rna.shape[1])},
        "n_targets": int(len(target_names)),
        "prompting_reference_rule": "RNA references and supervised heads are fit on selected protocol train patients only.",
        "splits": {},
    }
    for split_name in ("val", "test"):
        query_mask = split == split_name
        retrieval = paired_retrieval_metrics(wsi[query_mask], rna[query_mask]) if query_mask.any() else None
        molecular = _evaluate_target_split(split_name, query_mask, train_mask, patient_ids, wsi, rna, targets, target_names, output, k, tau)
        payload["splits"][split_name] = {"identity_retrieval": retrieval, "molecular_prompting": molecular}
    path = output / f"{protocol}_patch_benchmark_metrics.json"
    write_json(path, payload)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate an internal or external TCGA-UT patch-token embedding archive")
    parser.add_argument("--aligned-npz", required=True)
    parser.add_argument("--targets", required=True, help="Parquet/CSV with patient_id and numeric molecular targets")
    parser.add_argument("--protocol", choices=VALID_PROTOCOLS, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--tau", type=float, default=0.1)
    args = parser.parse_args()
    print(evaluate_protocol_benchmark(args.aligned_npz, args.targets, args.protocol, args.output_dir, args.k, args.tau))


if __name__ == "__main__":
    main()
