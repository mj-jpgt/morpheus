"""Configuration loading for MORPHEUS v1 proof scripts."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


@dataclass(frozen=True)
class MorpheusConfig:
    """Resolved project configuration."""

    raw: dict[str, Any]
    config_path: Path
    project_root: Path

    @classmethod
    def load(cls, path: str | Path = "morpheus/configs/v1.json") -> "MorpheusConfig":
        config_path = Path(path)
        with config_path.open("r", encoding="utf-8") as handle:
            raw = json.load(handle)
        project_root = Path(raw.get("project_root", ".")).resolve()
        return cls(raw=raw, config_path=config_path.resolve(), project_root=project_root)

    def path(self, key: str) -> Path:
        value = self.raw.get("paths", {}).get(key)
        if value is None:
            raise KeyError(f"Missing configured path: paths.{key}")
        path = Path(value)
        if not path.is_absolute():
            path = self.project_root / path
        return path

    def section(self, key: str) -> Mapping[str, Any]:
        value = self.raw.get(key, {})
        if not isinstance(value, Mapping):
            raise TypeError(f"Config section {key!r} must be an object")
        return value

    def ensure_dirs(self) -> None:
        for key in ("processed_dir", "outputs_dir"):
            self.path(key).mkdir(parents=True, exist_ok=True)


def load_config(path: str | Path = "morpheus/configs/v1.json") -> MorpheusConfig:
    return MorpheusConfig.load(path)
