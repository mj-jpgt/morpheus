"""Train/export strictly matched E1 decorrelation arms.

The legacy diagnostics are useful for locating the observation, but only this
driver can produce the matched three-seed artifacts required for an E1 claim.
It deliberately uses a fixed, predeclared epoch count and a development-only
refit: the outer held-out cancers are never consulted during selection.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def _command(*, module: str, values: list[str]) -> list[str]:
    return [sys.executable, "-m", module, *values]


def _read_liveness(training_dir: Path) -> dict:
    path = training_dir / "train_metrics.jsonl"
    if not path.is_file():
        raise RuntimeError(f"training did not write metrics: {path}")
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    losses = [float(row["train_loss"]) for row in rows if "train_loss" in row]
    finite = bool(losses) and all(value == value and abs(value) != float("inf") for value in losses)
    reduced = bool(finite and len(losses) >= 2 and losses[-1] <= 0.8 * losses[0])
    return {"n_epochs_logged": len(rows), "loss_finite": finite, "loss_initial": losses[0] if losses else None,
            "loss_final": losses[-1] if losses else None, "loss_reduced_20_percent": reduced}


def run_arm(args: argparse.Namespace, *, seed: int, arm: str, decorrelation_weight: float) -> dict:
    root = Path(args.output).resolve(); training = root / "training" / f"seed{seed}" / arm
    artifact = root / "artifacts" / f"e1_seed{seed}_{arm}.npz"
    training.mkdir(parents=True, exist_ok=True); artifact.parent.mkdir(parents=True, exist_ok=True)
    runner = _command(module="morpheus.v2.runner", values=[
        "--data-config", args.data_config, "--split-file", args.split_file, "--output-dir", str(training),
        "--epochs", str(args.epochs), "--token-budget", str(args.token_budget), "--seed", str(seed),
        "--objective-profile", "programme_only", "--decorrelation-weight", str(decorrelation_weight),
        "--gradient-diagnostics-every", str(args.gradient_diagnostics_every), "--fit-development",
        "--expected-development-cancers", str(args.expected_development_cancers),
        "--expected-heldout-cancers", str(args.expected_heldout_cancers), "--device", args.device,
    ])
    if not (training / "last.pt").is_file():
        with (training / "runner.log").open("w", encoding="utf-8") as handle:
            result = subprocess.run(runner, stdout=handle, stderr=subprocess.STDOUT, text=True)
        if result.returncode:
            raise RuntimeError(f"E1 {arm} seed {seed} training failed; see {training / 'runner.log'}")
    liveness = _read_liveness(training)
    # Do not quietly accept an inert arm.  A failed liveness check is a failed
    # experiment, not evidence that the regularizer had no effect.
    if not liveness["loss_finite"] or not liveness["loss_reduced_20_percent"]:
        raise RuntimeError(f"E1 {arm} seed {seed} fails G2 loss liveness: {liveness}")
    if not artifact.is_file():
        exporter = _command(module="morpheus.v2.export", values=[
            "--data-config", args.data_config, "--split-file", args.split_file,
            "--checkpoint", str(training / "last.pt"), "--output", str(artifact),
            "--token-budget", str(args.token_budget), "--device", args.device,
        ])
        with (training / "export.log").open("w", encoding="utf-8") as handle:
            result = subprocess.run(exporter, stdout=handle, stderr=subprocess.STDOUT, text=True)
        if result.returncode:
            raise RuntimeError(f"E1 {arm} seed {seed} export failed; see {training / 'export.log'}")
    return {"seed": seed, "arm": arm, "decorrelation_weight": decorrelation_weight,
            "training_dir": str(training), "artifact": str(artifact), "liveness": liveness}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-config", required=True)
    parser.add_argument("--split-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--epochs", type=int, default=12)
    parser.add_argument("--token-budget", type=int, default=65536)
    parser.add_argument("--before-weight", type=float, default=0.0)
    parser.add_argument("--after-weight", type=float, default=0.04)
    parser.add_argument("--gradient-diagnostics-every", type=int, default=25)
    parser.add_argument("--expected-development-cancers", type=int, default=11)
    parser.add_argument("--expected-heldout-cancers", type=int, default=21)
    parser.add_argument("--device", default="auto")
    args = parser.parse_args()
    if args.epochs < 2 or args.token_budget <= 0:
        raise ValueError("E1 needs at least two epochs and a positive token budget")
    if args.before_weight != 0.0 or args.after_weight <= 0.0:
        raise ValueError("E1 intervention must compare 0.0 decorrelation against a positive predeclared weight")
    seeds = tuple(int(value) for value in args.seeds.split(",") if value.strip())
    if len(seeds) != 3 or len(set(seeds)) != 3:
        raise ValueError("E1 requires exactly three independent seeds")
    records = []
    for seed in seeds:
        records.append(run_arm(args, seed=seed, arm="before", decorrelation_weight=args.before_weight))
        records.append(run_arm(args, seed=seed, arm="after", decorrelation_weight=args.after_weight))
    output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    (output / "training_manifest.json").write_text(json.dumps({"experiment": "E1_matched_training", "seeds": seeds,
        "epochs": args.epochs, "token_budget": args.token_budget, "objective_profile": "programme_only",
        "intervention": {"before": args.before_weight, "after": args.after_weight}, "records": records}, indent=2))
    print(output / "training_manifest.json")


if __name__ == "__main__":
    main()
