"""Run the post-training E1 and E2 evidence stages transactionally.

This is intentionally a Python controller rather than a shell tail-chain.  It
does not decide a scientific outcome: each analysis owns its health gates and
returns non-zero only for an invalid measurement.  The controller merely makes
the dependency explicit: a complete, matched E1 training bundle is required
before any E1/E2 conclusion is written.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path
from typing import Any


def validate_training_manifest(path: Path) -> dict[str, Any]:
    """Fail closed unless the six planned E1 artifacts and liveness records exist."""
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("experiment") != "E1_matched_training":
        raise ValueError("training manifest is not an E1 matched-training bundle")
    seeds = tuple(int(seed) for seed in value.get("seeds", ()))
    if len(seeds) != 3 or len(set(seeds)) != 3:
        raise ValueError("E1 analysis requires exactly three independent training seeds")
    expected = {(seed, arm) for seed in seeds for arm in ("before", "after")}
    observed: set[tuple[int, str]] = set()
    for record in value.get("records", []):
        key = (int(record.get("seed", -1)), str(record.get("arm", "")))
        observed.add(key)
        artifact = Path(str(record.get("artifact", "")))
        liveness = record.get("liveness", {})
        if not artifact.is_file():
            raise FileNotFoundError(f"missing declared E1 artifact: {artifact}")
        if not (liveness.get("loss_finite") and liveness.get("loss_reduced_20_percent")):
            raise ValueError(f"E1 training liveness failed for {key}: {liveness}")
    if observed != expected:
        raise ValueError(f"incomplete E1 arm set: expected={sorted(expected)} observed={sorted(observed)}")
    return value


def _run(command: list[str], *, log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}); see {log}")


def _atomic_directory(final: Path, action) -> None:
    final.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{final.name}.", dir=final.parent))
    try:
        action(staging)
        if final.exists():
            shutil.rmtree(final)
        os.replace(staging, final)
    except Exception:
        failure = staging / "FAILED.json"
        failure.write_text(json.dumps({"error": traceback.format_exc(), "failed_at": time.time()}, indent=2), encoding="utf-8")
        if final.exists():
            shutil.rmtree(final)
        os.replace(staging, final)
        raise


def run_pipeline(args: argparse.Namespace) -> dict[str, Any]:
    training_root = Path(args.training_root).resolve()
    manifest = validate_training_manifest(training_root / "training_manifest.json")
    targets = Path(args.targets).resolve()
    if not targets.is_file():
        raise FileNotFoundError(targets)
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=True)
    records = {(int(row["seed"]), str(row["arm"])): row for row in manifest["records"]}
    seeds = tuple(int(seed) for seed in manifest["seeds"])

    def one_seed(seed: int) -> tuple[int, str]:
        final = out / "e1" / f"seed{seed}"
        def action(staging: Path) -> None:
            command = [sys.executable, "-m", "morpheus.v2.calibra.e1_rank_information",
                       "--before-artifact", str(records[(seed, "before")]["artifact"]),
                       "--after-artifact", str(records[(seed, "after")]["artifact"]),
                       "--targets", str(targets), "--output", str(staging), "--state", args.state,
                       "--partition", "test", "--n-components", str(args.n_components),
                       "--n-draws", str(args.n_draws), "--n-permutations", str(args.n_permutations),
                       "--n-bootstrap", str(args.n_bootstrap), "--seed", str(seed)]
            if args.official_gate_log:
                command += ["--official-gate-log", args.official_gate_log]
            _run(command, log=staging / "stage.log")
            gate = json.loads((staging / "gate_summary.json").read_text(encoding="utf-8"))
            if not gate.get("gates_pass", False):
                raise RuntimeError(f"E1 seed {seed} did not pass health gates")
        _atomic_directory(final, action)
        return seed, str(final)

    workers = max(1, min(int(args.e1_workers), len(seeds)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        finished = list(pool.map(one_seed, seeds))
    e1_paths = [Path(path) for _, path in sorted(finished)]

    def aggregate_action(staging: Path) -> None:
        command = [sys.executable, "-m", "morpheus.v2.calibra.aggregate_e1", "--runs",
                   *(str(path) for path in e1_paths), "--output", str(staging)]
        _run(command, log=staging / "stage.log")
        gate = json.loads((staging / "gate_summary.json").read_text(encoding="utf-8"))
        if not gate.get("gates_pass", False):
            raise RuntimeError("three-seed E1 aggregation missing valid gates")
    _atomic_directory(out / "e1_summary", aggregate_action)

    def e2_action(staging: Path) -> None:
        command = [sys.executable, "-m", "morpheus.v2.calibra.e2_expressible_intersection",
                   "--features", str(records[(seeds[0], "after")]["artifact"]), "--feature-key", args.state,
                   "--targets", str(targets), "--target-key", args.target_key, "--output", str(staging),
                   "--development-splits", "train,val", "--ks", args.e2_ks, "--seeds", ",".join(map(str, seeds)),
                   "--weight-decays", args.e2_weight_decays, "--latent-dim", str(args.e2_latent_dim),
                   "--epochs", str(args.e2_epochs), "--device", args.e2_device]
        if args.official_gate_log:
            command += ["--official-gate-log", args.official_gate_log]
        _run(command, log=staging / "stage.log")
        gate = json.loads((staging / "gate_summary.json").read_text(encoding="utf-8"))
        if not gate.get("gates_pass", False):
            raise RuntimeError("E2 did not pass health gates")
    _atomic_directory(out / "e2", e2_action)

    result = {"experiment": "E1_E2_evidence_pipeline", "training_manifest": str(training_root / "training_manifest.json"),
              "targets": str(targets), "seeds": seeds, "e1_workers": workers,
              "e1_summary": str(out / "e1_summary"), "e2": str(out / "e2"), "completed_at": time.time()}
    (out / "SUCCESS.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-root", required=True); parser.add_argument("--targets", required=True); parser.add_argument("--output", required=True)
    parser.add_argument("--state", default="wsi_biology"); parser.add_argument("--target-key", default="scores")
    parser.add_argument("--e1-workers", type=int, default=2); parser.add_argument("--n-components", type=int, default=32)
    parser.add_argument("--n-draws", type=int, default=25); parser.add_argument("--n-permutations", type=int, default=1000); parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--e2-ks", default="5,10,20,50,100"); parser.add_argument("--e2-weight-decays", default="0.0001,0")
    parser.add_argument("--e2-latent-dim", type=int, default=128); parser.add_argument("--e2-epochs", type=int, default=300); parser.add_argument("--e2-device", default="cuda")
    parser.add_argument("--official-gate-log", default="")
    args = parser.parse_args()
    if args.e1_workers < 1 or args.n_permutations < 100 or args.n_bootstrap < 100:
        raise ValueError("E1 needs >=1 worker and headline null/bootstrap resolution of at least 100")
    print(json.dumps(run_pipeline(args), indent=2))


if __name__ == "__main__":
    main()
