"""Evaluate V2 BioQueryFormer typed surfaces with the fair molecular protocol."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.src.eval.eval_hallmark_prediction_head import evaluate_hallmark_prediction_head
from morpheus.src.eval.eval_molecular_prompting import evaluate_molecular_prompting


DEFAULT_JOBS = (
    ("patient_stratified", Path("morpheus/outputs/v2_bio_query_former_patient_stratified_full/bio_query_former_embeddings.npz")),
    ("patient_heldout", Path("morpheus/outputs/v2_bio_query_former_patient_heldout_full/bio_query_former_embeddings.npz")),
    ("patch_stratified", Path("morpheus/outputs/v2_bio_query_former_patch_stratified_full/bio_query_former_embeddings.npz")),
    ("patient_stratified_stronger_all", Path("morpheus/outputs/v2_bio_query_former_patient_stratified_stronger_full/bio_query_former_embeddings.npz")),
    ("patient_stratified_hybrid_fn", Path("morpheus/outputs/v2_bio_query_former_patient_stratified_hybrid_fn_full/bio_query_former_embeddings.npz")),
    ("patient_heldout_hybrid_fn", Path("morpheus/outputs/v2_bio_query_former_patient_heldout_hybrid_fn_full/bio_query_former_embeddings.npz")),
)


FEATURE_PAIRS = {
    "identity": ("wsi_identity", "rna_identity"),
    "biology": ("wsi_biology", "rna_biology"),
    "full_biology_to_rna_biology": ("full_biology", "rna_biology"),
}


def _read_json(path: Path) -> dict:
    return json.load(open(path, encoding="utf-8"))


def _write_aligned_npz(source_npz: Path, output_path: Path, wsi_key: str, rna_key: str) -> Path:
    data = np.load(source_npz, allow_pickle=True)
    missing = [key for key in ("patient_ids", "split", wsi_key, rna_key) if key not in data.files]
    if missing:
        raise KeyError(f"{source_npz} is missing required keys for V2 molecular evaluation: {missing}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        patient_ids=data["patient_ids"],
        split=data["split"],
        cancers=data["cancers"] if "cancers" in data.files else np.asarray(["unknown"] * len(data["patient_ids"])),
        wsi=data[wsi_key].astype(np.float32),
        rna=data[rna_key].astype(np.float32),
    )
    return output_path


def _fmt(value) -> str:
    return "NA" if value is None else f"{float(value):.4f}"


def _run_one(config_path: str, output_root: Path, run_name: str, source_npz: Path, feature_name: str, wsi_key: str, rna_key: str) -> dict:
    aligned_npz = _write_aligned_npz(source_npz, output_root / "aligned" / run_name / f"{feature_name}_aligned_embeddings.npz", wsi_key, rna_key)
    out_base = output_root / f"{run_name}_{feature_name}"
    soft_path = evaluate_molecular_prompting(config_path, aligned_npz, out_base / "soft_knn")["global"]
    head_path = evaluate_hallmark_prediction_head(config_path, aligned_npz, out_base / "common_head", "wsi")
    soft = _read_json(Path(soft_path))
    head = _read_json(Path(head_path))
    return {
        "run": run_name,
        "feature": feature_name,
        "source_npz": str(source_npz),
        "aligned_npz": str(aligned_npz),
        "wsi_key": wsi_key,
        "rna_key": rna_key,
        "val_soft_knn_pearson": soft.get("val_metrics", {}).get("mean_pearson"),
        "test_soft_knn_pearson": soft.get("test_metrics", {}).get("mean_pearson"),
        "val_common_head_pearson": head.get("val_metrics", {}).get("mean_pearson"),
        "test_common_head_pearson": head.get("test_metrics", {}).get("mean_pearson"),
        "val_n": head.get("val_metrics", {}).get("n_eval"),
        "test_n": head.get("test_metrics", {}).get("n_eval"),
    }


def _write_report(rows: list[dict], output_root: Path) -> Path:
    frame = pd.DataFrame(rows)
    output_root.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_root / "v2_molecular_benchmark_summary.csv", index=False)
    lines = [
        "# V2 BioQueryFormer Molecular Benchmarks",
        "",
        "Rows use the same soft-kNN and train-only StandardScaler + RidgeCV protocol as the fair baseline report.",
        "",
        "| run | feature | val soft-kNN | test soft-kNN | val common head | test common head | test n |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["run"]),
                    str(row["feature"]),
                    _fmt(row["val_soft_knn_pearson"]),
                    _fmt(row["test_soft_knn_pearson"]),
                    _fmt(row["val_common_head_pearson"]),
                    _fmt(row["test_common_head_pearson"]),
                    str(row["test_n"]),
                ]
            )
            + " |"
        )
    path = output_root / "v2_molecular_benchmark_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_v2_molecular_benchmarks(config_path: str, output_root: str | Path, jobs: tuple[tuple[str, Path], ...] = DEFAULT_JOBS) -> Path:
    out = Path(output_root)
    rows = []
    for run_name, source_npz in jobs:
        if not source_npz.exists():
            continue
        for feature_name, (wsi_key, rna_key) in FEATURE_PAIRS.items():
            rows.append(_run_one(config_path, out, run_name, source_npz, feature_name, wsi_key, rna_key))
    if not rows:
        raise FileNotFoundError("No V2 embedding files were found for molecular benchmarking")
    return _write_report(rows, out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--output-root", default="morpheus/outputs/v2_molecular_benchmarks")
    args = parser.parse_args()
    print(run_v2_molecular_benchmarks(args.config, args.output_root))


if __name__ == "__main__":
    main()
