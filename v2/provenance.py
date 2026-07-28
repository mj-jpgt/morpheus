"""Immutable source and run provenance for MORPHEUS V2 artifacts.

The project has historically been run from more than one checkout.  This
module gives every controller, checkpoint, and representation export one
stable description of the source tree that produced it.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
from typing import Iterable


PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def _git(args: list[str]) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=PACKAGE_ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip() or None
    except (OSError, subprocess.CalledProcessError):
        return None


def source_tree_digest(root: str | Path = PACKAGE_ROOT) -> str:
    """Hash tracked V2 Python/configuration files independent of file mtimes."""
    base = Path(root)
    files: Iterable[Path] = sorted(
        path for path in base.rglob("*")
        if path.is_file() and path.suffix in {".py", ".yaml", ".yml", ".json", ".sh"}
        and "__pycache__" not in path.parts
    )
    digest = sha256()
    for path in files:
        relative = path.relative_to(base).as_posix().encode("utf-8")
        digest.update(relative + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()


def source_manifest(*, configuration: dict[str, object] | None = None) -> dict[str, object]:
    """Return serialisable provenance to embed in every run artifact."""
    return {
        "package": "morpheus.v2",
        "package_root": str(PACKAGE_ROOT),
        "git_commit": _git(["rev-parse", "HEAD"]),
        "git_dirty": bool(_git(["status", "--porcelain"])),
        "source_tree_sha256": source_tree_digest(),
        "configuration_sha256": sha256(
            json.dumps(configuration or {}, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
    }


def write_source_manifest(destination: str | Path, *, configuration: dict[str, object] | None = None) -> dict[str, object]:
    """Write and return a frozen source manifest before a controller starts work."""
    manifest = source_manifest(configuration=configuration)
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
