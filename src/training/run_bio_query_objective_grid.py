"""Run MORPHEUS BioQueryFormer objective ablations on fixed splits."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd


GRID = {
    "baseline_infonce": [],
    "hybrid_fn": ["--alignment-loss", "hybrid", "--false-negative-aware"],
    "siglip_fn": ["--alignment-loss", "siglip", "--false-negative-aware"],
    "program_head": ["--program-head-weight", "0.25"],
    "hybrid_fn_program": ["--alignment-loss", "hybrid", "--false-negative-aware", "--program-head-weight", "0.25"],
    "cancer_adv_light": ["--cancer-adv-weight", "0.02"],
    "pls_light": ["--pls-distill-weight", "0.05", "--pls-components", "16"],
}


def _command(args, name: str, extra: list[str]) -> list[str]:
    return [
        args.python,
        "-u",
        "-m",
        "morpheus.src.training.train_bio_query_former",
        "--config",
        args.config,
        "--split-file",
        args.split_file,
        "--wsi-mode",
        args.wsi_mode,
        "--epochs",
        str(args.epochs),
        "--batch-size",
        str(args.batch_size),
        "--device",
        args.device,
        "--hidden-dim",
        str(args.hidden_dim),
        "--num-layers",
        str(args.num_layers),
        "--num-heads",
        str(args.num_heads),
        "--eval-batch-size",
        str(args.eval_batch_size),
        "--output-dir",
        str(Path(args.output_root) / name),
        *extra,
    ]


def _read_metrics(path: Path) -> dict:
    if not path.exists():
        return {}
    payload = json.load(open(path, encoding="utf-8"))
    test = payload.get("test_metrics", {})
    val = payload.get("val_metrics", {})
    return {
        "best_epoch": payload.get("best_epoch"),
        "val_r10": val.get("retrieval_r10"),
        "val_mrr": val.get("retrieval_mrr"),
        "val_hallmark": val.get("hallmark_wsi_biology_pearson"),
        "test_r10": test.get("retrieval_r10"),
        "test_mrr": test.get("retrieval_mrr"),
        "test_hallmark": test.get("hallmark_wsi_biology_pearson"),
    }


def run_grid(args) -> Path:
    out = Path(args.output_root)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    selected = args.jobs or list(GRID)
    for name in selected:
        if name not in GRID:
            raise KeyError(f"Unknown grid job {name!r}; choices: {sorted(GRID)}")
        cmd = _command(args, name, GRID[name])
        job_dir = out / name
        row = {"job": name, "command": " ".join(cmd), "output_dir": str(job_dir)}
        if not args.dry_run:
            subprocess.run(cmd, check=True)
            row.update(_read_metrics(job_dir / "test_metrics.json"))
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary_path = out / "objective_grid_summary.csv"
    summary.to_csv(summary_path, index=False)
    (out / "objective_grid_manifest.json").write_text(json.dumps({"jobs": rows}, indent=2), encoding="utf-8")
    return summary_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=".\\PRISM\\.venv\\Scripts\\python.exe")
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--split-file", default="morpheus/data/processed/splits/tumor_state_stratified.json")
    parser.add_argument("--output-root", default="morpheus/outputs/v2_objective_grid")
    parser.add_argument("--wsi-mode", choices=["patient", "patch"], default="patient")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", default="cuda", choices=["auto", "cuda", "cpu"])
    parser.add_argument("--hidden-dim", type=int, default=256)
    parser.add_argument("--num-layers", type=int, default=1)
    parser.add_argument("--num-heads", type=int, default=8)
    parser.add_argument("--eval-batch-size", type=int, default=32)
    parser.add_argument("--jobs", nargs="*")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    print(run_grid(args))


if __name__ == "__main__":
    main()
