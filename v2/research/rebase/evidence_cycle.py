"""Atomic, resumable controller for the CLD evidence cycle.

This deliberately owns orchestration only.  Scientific stages remain normal
CLIs, which makes their inputs and outputs inspectable without asking a shell
pipeline to silently decide what a failed experiment means.
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_digest(value: Any) -> str:
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
    gate_file: str = ""

    @classmethod
    def parse(cls, value: dict[str, Any]) -> "Stage":
        missing = {"name", "command", "inputs", "outputs"} - set(value)
        if missing:
            raise ValueError(f"stage missing keys: {sorted(missing)}")
        command = value["command"]
        if not isinstance(command, list) or not all(isinstance(part, str) for part in command):
            raise ValueError("stage command must be a string list")
        return cls(
            name=str(value["name"]), command=tuple(command),
            inputs=tuple(map(str, value["inputs"])), outputs=tuple(map(str, value["outputs"])),
            required=bool(value.get("required", True)), gate_file=str(value.get("gate_file", "")),
        )


def _resolve(repo: Path, item: str) -> Path:
    candidate = Path(item)
    return candidate if candidate.is_absolute() else repo / candidate


def _code_digest(repo: Path) -> str:
    """Hash executable source only; never traverse data, outputs, or cloud mounts."""
    roots = (repo / "v2", repo / "src")
    files = sorted(path for root in roots if root.is_dir() for path in root.rglob("*.py"))
    return _json_digest({str(path.relative_to(repo)): _sha256(path) for path in files})


def _fingerprint(stage: Stage, repo: Path, config_digest: str, code_digest: str) -> dict[str, Any]:
    inputs: dict[str, str] = {}
    for item in stage.inputs:
        path = _resolve(repo, item)
        if not path.is_file():
            raise FileNotFoundError(f"{stage.name}: missing required input {path}")
        inputs[item] = _sha256(path)
    return {"schema_version": SCHEMA_VERSION, "stage": stage.name, "command": list(stage.command),
            "input_digests": inputs, "config_digest": config_digest, "code_digest": code_digest}


def _completed(directory: Path, fingerprint: dict[str, Any]) -> bool:
    marker = directory / "SUCCESS.json"
    if not marker.is_file():
        return False
    try:
        value = json.loads(marker.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    return value.get("fingerprint") == fingerprint and all((directory / output).is_file() for output in value.get("outputs", []))


def _stage_gate_passed(directory: Path, gate_file: str) -> bool:
    if not gate_file:
        return True
    path = directory / gate_file
    if not path.is_file():
        return False
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    return bool(value.get("gates_pass", False))


def _run(stage: Stage, repo: Path, state_root: Path, fingerprint: dict[str, Any], resume: bool) -> dict[str, Any]:
    final = state_root / stage.name
    if resume and _completed(final, fingerprint):
        return {"stage": stage.name, "status": "resumed", "path": str(final)}
    if final.exists():
        shutil.rmtree(final)
    staging = Path(tempfile.mkdtemp(prefix=f".{stage.name}.", dir=state_root))
    try:
        replacements = {"{stage_dir}": str(staging), "{repo}": str(repo), "{python}": sys.executable}
        command = []
        for part in stage.command:
            for source, target in replacements.items():
                part = part.replace(source, target)
            command.append(part)
        with (staging / "stage.log").open("w", encoding="utf-8") as handle:
            result = subprocess.run(command, cwd=repo, stdout=handle, stderr=subprocess.STDOUT, text=True)
        if result.returncode:
            raise RuntimeError(f"stage command failed with exit={result.returncode}")
        missing = [output for output in stage.outputs if not (staging / output).is_file()]
        if missing:
            raise RuntimeError(f"stage did not emit required outputs: {missing}")
        if not _stage_gate_passed(staging, stage.gate_file):
            raise RuntimeError(f"stage gates failed or missing: {stage.gate_file}")
        _atomic_json(staging / "SUCCESS.json", {"fingerprint": fingerprint, "outputs": list(stage.outputs),
                                                  "completed_at": time.time()})
        os.replace(staging, final)
        return {"stage": stage.name, "status": "complete", "path": str(final)}
    except Exception as exc:
        _atomic_json(staging / "FAILED.json", {"fingerprint": fingerprint, "error": str(exc),
                                                 "traceback": traceback.format_exc(), "failed_at": time.time()})
        os.replace(staging, final)
        if stage.required:
            raise
        return {"stage": stage.name, "status": "unavailable", "path": str(final), "error": str(exc)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    config = json.loads(args.config.read_text(encoding="utf-8"))
    stages = [Stage.parse(value) for value in config.get("stages", [])]
    if not stages or len({stage.name for stage in stages}) != len(stages):
        raise ValueError("config needs a non-empty list of uniquely named stages")
    run_root = args.run_root.resolve(); run_root.mkdir(parents=True, exist_ok=True)
    state_root = run_root / "stages"; state_root.mkdir(exist_ok=True)
    config_digest, code_digest = _sha256(args.config), _code_digest(repo)
    _atomic_json(run_root / "run_manifest.json", {"schema_version": SCHEMA_VERSION, "config": str(args.config),
                                                    "config_digest": config_digest, "code_digest": code_digest,
                                                    "python": sys.version, "stages": [s.name for s in stages]})
    completed: list[dict[str, Any]] = []
    for stage in stages:
        fingerprint = _fingerprint(stage, repo, config_digest, code_digest)
        completed.append(_run(stage, repo, state_root, fingerprint, args.resume))
        _atomic_json(run_root / "dag_state.json", {"complete": False, "current_stage": stage.name,
                                                     "results": completed})
    _atomic_json(run_root / "dag_state.json", {"complete": True, "current_stage": None, "results": completed})


if __name__ == "__main__":
    main()
