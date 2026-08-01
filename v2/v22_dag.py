"""Manifested, resumable execution controller for the MORPHEUS V2.2 run.

The controller is deliberately configuration-driven: expensive scientific
stages remain normal Python CLIs, while this module owns the reproducibility
contract, atomic promotion, and resume semantics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1


def _digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _digest_json(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)


@dataclass(frozen=True)
class Stage:
    name: str
    command: tuple[str, ...]
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    required: bool = True

    @classmethod
    def parse(cls, payload: dict[str, Any]) -> "Stage":
        required = {"name", "command", "inputs", "outputs"}
        missing = required - set(payload)
        if missing:
            raise ValueError(f"stage is missing {sorted(missing)}")
        command = payload["command"]
        if not isinstance(command, list) or not all(isinstance(item, str) for item in command):
            raise ValueError("stage command must be a list of strings")
        return cls(str(payload["name"]), tuple(command), tuple(map(str, payload["inputs"])),
                   tuple(map(str, payload["outputs"])), bool(payload.get("required", True)))


def _fingerprint(stage: Stage, root: Path, code_digest: str, config_digest: str) -> dict[str, Any]:
    inputs: dict[str, str] = {}
    for value in stage.inputs:
        candidate = Path(value)
        path = candidate if candidate.is_absolute() else root / candidate
        if not path.is_file():
            raise FileNotFoundError(f"stage {stage.name} input is missing: {path}")
        inputs[value] = _digest_file(path)
    return {"schema_version": SCHEMA_VERSION, "stage": stage.name, "command": list(stage.command),
            "input_digests": inputs, "code_digest": code_digest, "config_digest": config_digest}


def _complete(stage_dir: Path, fingerprint: dict[str, Any]) -> bool:
    marker = stage_dir / "SUCCESS.json"
    if not marker.is_file():
        return False
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    if payload.get("fingerprint") != fingerprint:
        return False
    return all((stage_dir / output).is_file() for output in payload.get("outputs", []))


def _run_stage(stage: Stage, root: Path, state_root: Path, fingerprint: dict[str, Any], resume: bool) -> dict[str, Any]:
    final = state_root / stage.name
    if resume and _complete(final, fingerprint):
        return {"stage": stage.name, "status": "resumed", "stage_dir": str(final)}
    if final.exists():
        shutil.rmtree(final)
    staging = Path(tempfile.mkdtemp(prefix=f".{stage.name}.", dir=state_root))
    try:
        command = [item.replace("{stage_dir}", str(staging)).replace("{root}", str(root)) for item in stage.command]
        log = staging / "stage.log"
        with log.open("w", encoding="utf-8") as handle:
            result = subprocess.run(command, cwd=root, stdout=handle, stderr=subprocess.STDOUT, text=True)
        if result.returncode:
            raise RuntimeError(f"exit={result.returncode}; see {log}")
        missing = [output for output in stage.outputs if not (staging / output).is_file()]
        if missing:
            raise RuntimeError(f"stage did not produce required outputs: {missing}")
        _atomic_json(staging / "SUCCESS.json", {"fingerprint": fingerprint, "outputs": list(stage.outputs), "completed_at": time.time()})
        os.replace(staging, final)
        return {"stage": stage.name, "status": "complete", "stage_dir": str(final)}
    except Exception as exc:
        _atomic_json(staging / "FAILED.json", {"fingerprint": fingerprint, "error": str(exc), "traceback": traceback.format_exc(), "failed_at": time.time()})
        os.replace(staging, final)
        if stage.required:
            raise
        return {"stage": stage.name, "status": "unavailable", "stage_dir": str(final), "error": str(exc)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path, help="immutable V2.2 DAG JSON")
    parser.add_argument("--run-root", required=True, type=Path)
    # Hash the V2 source package, not the whole repository: recursively
    # walking a repository-level data directory makes a lightweight preflight
    # unexpectedly scan terabytes of WSI artifacts.
    parser.add_argument("--code-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    stages = [Stage.parse(item) for item in config.get("stages", [])]
    if not stages or len({stage.name for stage in stages}) != len(stages):
        raise ValueError("DAG config needs uniquely named stages")
    root = args.run_root.resolve(); root.mkdir(parents=True, exist_ok=True)
    state_root = root / "stages"; state_root.mkdir(exist_ok=True)
    config_digest = _digest_file(args.config)
    # V2 is intentionally a flat package.  A non-recursive scan avoids slow
    # or hanging traversal through cloud-synchronised cache/data directories.
    code_files = sorted(args.code_root.glob("*.py"))
    code_digest = _digest_json({str(path.relative_to(args.code_root)): _digest_file(path) for path in code_files})
    _atomic_json(root / "run_manifest.json", {"schema_version": SCHEMA_VERSION, "config": str(args.config),
                                                "config_digest": config_digest, "code_digest": code_digest,
                                                "python": sys.version, "stages": [stage.name for stage in stages]})
    results = []
    for stage in stages:
        fingerprint = _fingerprint(stage, root, code_digest, config_digest)
        results.append(_run_stage(stage, root, state_root, fingerprint, args.resume))
        _atomic_json(root / "dag_state.json", {"schema_version": SCHEMA_VERSION, "results": results,
                                                 "current_stage": stage.name, "complete": False})
    _atomic_json(root / "dag_state.json", {"schema_version": SCHEMA_VERSION, "results": results,
                                             "current_stage": None, "complete": True})


if __name__ == "__main__":
    main()
