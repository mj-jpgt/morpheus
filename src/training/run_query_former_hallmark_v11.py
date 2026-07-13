"""Run fixed-split QueryFormer Hallmark V1.1 experiments and fair evaluations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from morpheus.src.eval.eval_hallmark_prediction_head import evaluate_hallmark_prediction_head
from morpheus.src.eval.eval_molecular_prompting import evaluate_molecular_prompting
from morpheus.src.training.train_query_former import run_query_former_training
from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import write_json


SPLITS = {
    "stratified": "tumor_state_stratified.json",
    "heldout": "tumor_state_heldout_cancer.json",
}

VARIANTS = (
    "wsi_hallmark_loss",
    "wsi_rna_hallmark_loss",
    "neighborhood_distill",
    "teacher_distill",
    "combined_best",
    "retentive_hallmark",
    "strong_teacher_hallmark",
    "teacher_hallmark_x1",
    "teacher_hallmark_x5",
    "wsi_hallmark_strong",
)


def _default_teacher_npz(output_root: Path, split_name: str) -> Path:
    baseline_name = "baselines_heldout" if split_name == "heldout" else "baselines_stratified"
    return output_root / "v1_same_split_comparison" / baseline_name / "cca_aligned_embeddings.npz"


def _read_json(path: Path) -> dict:
    return json.load(open(path, encoding="utf-8"))


def _mean_pearson(path: Path, split: str = "test") -> float | None:
    payload = _read_json(path)
    metrics = payload.get(f"{split}_metrics", payload if split == "test" else {})
    return metrics.get("mean_pearson")


def _run_one(
    config_path: str,
    split_name: str,
    variant: str,
    epochs: int,
    batch_size: int,
    device: str,
    smoke: bool,
    output_root: Path,
) -> dict:
    cfg = load_config(config_path)
    split_file = cfg.path("processed_dir") / "splits" / SPLITS[split_name]
    out_dir = output_root / f"query_former_{variant}_{split_name}"
    teacher_npz = None
    if variant in {"teacher_distill", "combined_best", "strong_teacher_hallmark", "teacher_hallmark_x1", "teacher_hallmark_x5"}:
        teacher_path = _default_teacher_npz(cfg.path("outputs_dir"), split_name)
        if teacher_path.exists():
            teacher_npz = str(teacher_path)
    metrics_path = run_query_former_training(
        config_path=config_path,
        split_file=str(split_file),
        output_dir=out_dir,
        epochs=epochs,
        batch_size=batch_size,
        smoke=smoke,
        device_name=device,
        variant=variant,
        teacher_npz=teacher_npz,
        teacher_key="wsi",
    )
    aligned_npz = out_dir / "query_former_aligned_embeddings.npz"
    molecular_dir = output_root / f"molecular_{variant}_{split_name}"
    head_dir = output_root / f"hallmark_head_{variant}_{split_name}"
    molecular_global = evaluate_molecular_prompting(config_path, aligned_npz, molecular_dir)["global"]
    head_global = evaluate_hallmark_prediction_head(config_path, aligned_npz, head_dir, "wsi")
    train_payload = _read_json(Path(metrics_path))
    molecular_payload = _read_json(Path(molecular_global))
    head_payload = _read_json(Path(head_global))
    summary = {
        "split": split_name,
        "variant": variant,
        "output_dir": str(out_dir),
        "checkpoint": train_payload.get("checkpoint"),
        "teacher_npz": teacher_npz,
        "epochs_run": train_payload.get("epochs_run"),
        "val_retrieval_r10": train_payload.get("val_metrics", {}).get("retrieval_r10"),
        "test_retrieval_r10": train_payload.get("test_metrics", {}).get("retrieval_r10"),
        "val_retrieval_mrr": train_payload.get("val_metrics", {}).get("retrieval_mrr"),
        "test_retrieval_mrr": train_payload.get("test_metrics", {}).get("retrieval_mrr"),
        "val_model_wsi_hallmark_pearson": train_payload.get("val_metrics", {}).get("hallmark_wsi_pearson"),
        "test_model_wsi_hallmark_pearson": train_payload.get("test_metrics", {}).get("hallmark_wsi_pearson"),
        "val_soft_knn_pearson": molecular_payload.get("val_metrics", {}).get("mean_pearson"),
        "test_soft_knn_pearson": molecular_payload.get("test_metrics", {}).get("mean_pearson"),
        "val_common_head_pearson": head_payload.get("val_metrics", {}).get("mean_pearson"),
        "test_common_head_pearson": head_payload.get("test_metrics", {}).get("mean_pearson"),
        "n_test_hallmark": head_payload.get("test_metrics", {}).get("n_eval"),
    }
    write_json(out_dir / "v11_fair_eval_summary.json", summary)
    return summary


def _write_report(rows: list[dict], out_path: Path) -> Path:
    frame = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path.with_suffix(".csv"), index=False)
    lines = [
        "# QueryFormer Hallmark V1.1 Results",
        "",
        "Fixed splits were used throughout. The common Hallmark head is train-only StandardScaler + RidgeCV; soft-kNN uses train references only.",
        "",
        "| split | variant | val soft-kNN | test soft-kNN | val common head | test common head | test R@10 | test MRR |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        def fmt(value):
            return "NA" if value is None else f"{float(value):.4f}"

        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["split"]),
                    str(row["variant"]),
                    fmt(row["val_soft_knn_pearson"]),
                    fmt(row["test_soft_knn_pearson"]),
                    fmt(row["val_common_head_pearson"]),
                    fmt(row["test_common_head_pearson"]),
                    fmt(row["test_retrieval_r10"]),
                    fmt(row["test_retrieval_mrr"]),
                ]
            )
            + " |"
        )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--output-root", default="morpheus/outputs/v1_query_former_hallmark_v11")
    parser.add_argument("--splits", nargs="+", default=["stratified", "heldout"], choices=sorted(SPLITS))
    parser.add_argument("--variants", nargs="+", default=list(VARIANTS), choices=list(VARIANTS))
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    output_root = Path(args.output_root)
    rows = []
    for split_name in args.splits:
        for variant in args.variants:
            rows.append(_run_one(args.config, split_name, variant, args.epochs, args.batch_size, args.device, args.smoke, output_root))
    report = _write_report(rows, output_root / "hallmark_v11_report.md")
    print(report)


if __name__ == "__main__":
    main()
