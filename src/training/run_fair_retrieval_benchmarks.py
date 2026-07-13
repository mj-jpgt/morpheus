"""Run corrected WSI-RNA retrieval benchmarks without shared clinical side channels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from morpheus.src.training.train_wsi_rna_alignment import run_alignment
from morpheus.src.utils.config import load_config


SPLITS = {
    "stratified": "tumor_state_stratified.json",
    "heldout": "tumor_state_heldout_cancer.json",
}

DEFAULT_METHODS = (
    "ridge_grid",
    "cca_grid",
    "pls_grid",
    "procrustes_grid",
    "mlp_siglip",
    "mlp_debiased",
)


def _read_json(path: Path) -> dict:
    return json.load(open(path, encoding="utf-8"))


def _fmt(value) -> str:
    if value is None:
        return "NA"
    return f"{float(value):.4f}"


def _run_method(config_path: str, split_name: str, method: str, output_root: Path) -> dict:
    cfg = load_config(config_path)
    split_file = cfg.path("processed_dir") / "splits" / SPLITS[split_name]
    out_dir = output_root / f"{split_name}_{method}"
    metrics_path = run_alignment(config_path, method, split_file, out_dir)
    payload = _read_json(Path(metrics_path))
    return {
        "split": split_name,
        "method": method,
        "output_dir": str(out_dir),
        "details": payload.get("details", {}),
        "train_r10": payload.get("train_metrics", {}).get("retrieval_r10"),
        "train_mrr": payload.get("train_metrics", {}).get("retrieval_mrr"),
        "val_r10": payload.get("val_metrics", {}).get("retrieval_r10"),
        "val_mrr": payload.get("val_metrics", {}).get("retrieval_mrr"),
        "test_r10": payload.get("test_metrics", {}).get("retrieval_r10"),
        "test_mrr": payload.get("test_metrics", {}).get("retrieval_mrr"),
        "test_r1": payload.get("test_metrics", {}).get("retrieval_r1"),
        "test_median_rank": payload.get("test_metrics", {}).get("median_rank"),
        "test_same_cancer_at_10": payload.get("test_metrics", {}).get("same_cancer_at_10"),
        "n_train": payload.get("n_train"),
        "n_val": payload.get("n_val"),
        "n_test": payload.get("n_test"),
    }


def _load_queryformer_rows(output_root: Path) -> list[dict]:
    rows = []
    candidates = [
        ("stratified", "query_former_no_clinical", Path("morpheus/outputs/v1_query_former_no_clinical_stratified/test_metrics.json")),
        ("heldout", "query_former_no_clinical", Path("morpheus/outputs/v1_query_former_no_clinical_heldout/test_metrics.json")),
    ]
    for split_name, method, path in candidates:
        if not path.exists():
            continue
        payload = _read_json(path)
        rows.append(
            {
                "split": split_name,
                "method": method,
                "output_dir": str(path.parent),
                "details": {"modality_specific_clinical": payload.get("modality_specific_clinical")},
                "train_r10": None,
                "train_mrr": None,
                "val_r10": payload.get("val_metrics", {}).get("retrieval_r10"),
                "val_mrr": payload.get("val_metrics", {}).get("retrieval_mrr"),
                "test_r10": payload.get("test_metrics", {}).get("retrieval_r10"),
                "test_mrr": payload.get("test_metrics", {}).get("retrieval_mrr"),
                "test_r1": payload.get("test_metrics", {}).get("retrieval_r1"),
                "test_median_rank": None,
                "test_same_cancer_at_10": payload.get("test_metrics", {}).get("same_cancer_at_10"),
                "n_train": payload.get("n_train"),
                "n_val": payload.get("n_val"),
                "n_test": payload.get("n_test"),
            }
        )
    return rows


def _write_report(rows: list[dict], output_root: Path) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    csv_path = output_root / "fair_retrieval_benchmark_summary.csv"
    frame.to_csv(csv_path, index=False)
    lines = [
        "# Corrected WSI-RNA Retrieval Benchmarks",
        "",
        "Old high QueryFormer retrieval numbers are quarantined because the modality-specific views shared clinical tokens. This report uses no shared clinical side channel for WSI-vs-RNA retrieval.",
        "",
        "| split | method | val R@10 | val MRR | test R@1 | test R@10 | test MRR | same cancer@10 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["split"]),
                    str(row["method"]),
                    _fmt(row["val_r10"]),
                    _fmt(row["val_mrr"]),
                    _fmt(row["test_r1"]),
                    _fmt(row["test_r10"]),
                    _fmt(row["test_mrr"]),
                    _fmt(row["test_same_cancer_at_10"]),
                ]
            )
            + " |"
        )
    report_path = output_root / "fair_retrieval_benchmark_report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--output-root", default="morpheus/outputs/fair_retrieval_benchmarks")
    parser.add_argument("--splits", nargs="+", default=["stratified", "heldout"], choices=sorted(SPLITS))
    parser.add_argument("--methods", nargs="+", default=list(DEFAULT_METHODS))
    parser.add_argument("--include-queryformer", action="store_true")
    args = parser.parse_args()

    output_root = Path(args.output_root)
    rows = []
    if args.include_queryformer:
        rows.extend(_load_queryformer_rows(output_root))
    for split_name in args.splits:
        for method in args.methods:
            rows.append(_run_method(args.config, split_name, method, output_root))
    report = _write_report(rows, output_root)
    print(report)


if __name__ == "__main__":
    main()
