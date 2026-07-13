"""Evaluate Hallmark prompting on corrected fair retrieval benchmark embeddings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from morpheus.src.eval.eval_hallmark_prediction_head import evaluate_hallmark_prediction_head
from morpheus.src.eval.eval_molecular_prompting import evaluate_molecular_prompting


def _read_json(path: Path) -> dict:
    return json.load(open(path, encoding="utf-8"))


def _embedding_jobs(retrieval_root: Path) -> list[tuple[str, str, Path]]:
    jobs = []
    qf = {
        "stratified": Path("morpheus/outputs/v1_query_former_no_clinical_stratified/query_former_aligned_embeddings.npz"),
        "heldout": Path("morpheus/outputs/v1_query_former_no_clinical_heldout/query_former_aligned_embeddings.npz"),
    }
    for split, path in qf.items():
        if path.exists():
            jobs.append((split, "query_former_no_clinical", path))
    for npz in sorted(retrieval_root.glob("*_*/*_aligned_embeddings.npz")):
        stem = npz.stem.removesuffix("_aligned_embeddings")
        split = npz.parent.name.split("_", 1)[0]
        method = stem
        jobs.append((split, method, npz))
    return jobs


def _fmt(value) -> str:
    return "NA" if value is None else f"{float(value):.4f}"


def _run_job(config_path: str, output_root: Path, split: str, method: str, aligned_npz: Path) -> dict:
    out_base = output_root / f"{split}_{method}"
    soft_path = evaluate_molecular_prompting(config_path, aligned_npz, out_base / "soft_knn")["global"]
    head_path = evaluate_hallmark_prediction_head(config_path, aligned_npz, out_base / "common_head", "wsi")
    soft = _read_json(Path(soft_path))
    head = _read_json(Path(head_path))
    return {
        "split": split,
        "method": method,
        "aligned_npz": str(aligned_npz),
        "val_soft_knn_pearson": soft.get("val_metrics", {}).get("mean_pearson"),
        "test_soft_knn_pearson": soft.get("test_metrics", {}).get("mean_pearson"),
        "val_common_head_pearson": head.get("val_metrics", {}).get("mean_pearson"),
        "test_common_head_pearson": head.get("test_metrics", {}).get("mean_pearson"),
        "val_n": head.get("val_metrics", {}).get("n_eval"),
        "test_n": head.get("test_metrics", {}).get("n_eval"),
    }


def _write_report(rows: list[dict], output_root: Path) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_root / "fair_molecular_benchmark_summary.csv", index=False)
    lines = [
        "# Corrected Hallmark Molecular Benchmarks",
        "",
        "All rows use corrected no-shared-clinical WSI/RNA retrieval embeddings. Hallmark heads are train-only StandardScaler + RidgeCV.",
        "",
        "| split | method | val soft-kNN | test soft-kNN | val common head | test common head | test n |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["split"]),
                    str(row["method"]),
                    _fmt(row["val_soft_knn_pearson"]),
                    _fmt(row["test_soft_knn_pearson"]),
                    _fmt(row["val_common_head_pearson"]),
                    _fmt(row["test_common_head_pearson"]),
                    str(row["test_n"]),
                ]
            )
            + " |"
        )
    path = output_root / "fair_molecular_benchmark_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--retrieval-root", default="morpheus/outputs/fair_retrieval_benchmarks")
    parser.add_argument("--output-root", default="morpheus/outputs/fair_molecular_benchmarks")
    args = parser.parse_args()
    retrieval_root = Path(args.retrieval_root)
    output_root = Path(args.output_root)
    rows = [_run_job(args.config, output_root, split, method, npz) for split, method, npz in _embedding_jobs(retrieval_root)]
    print(_write_report(rows, output_root))


if __name__ == "__main__":
    main()
