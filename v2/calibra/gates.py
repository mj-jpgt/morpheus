"""Small, append-only health-gate ledger used by CLD experiments."""
from __future__ import annotations

import csv
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class GateLedger:
    columns = ("timestamp_utc", "experiment", "gate", "value", "threshold", "status", "note")

    def __init__(self, output: str | Path, experiment: str):
        self.output = Path(output); self.output.mkdir(parents=True, exist_ok=True)
        self.path = self.output / "gate_rows.csv"; self.experiment = experiment
        self.rows: list[dict[str, Any]] = []

    def add(self, gate: str, value: Any, threshold: Any, passed: bool, note: str = "") -> None:
        self.rows.append({"timestamp_utc": datetime.now(timezone.utc).isoformat(), "experiment": self.experiment,
                          "gate": gate, "value": value, "threshold": threshold,
                          "status": "PASS" if passed else "FAIL", "note": note})

    def artifact(self, gate: str, path: str | Path) -> None:
        p = Path(path); stat = p.stat()
        digest = hashlib.sha256()
        with p.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        digest = digest.hexdigest()[:16]
        self.add(gate, f"size={stat.st_size};mtime={stat.st_mtime_ns};sha={digest}", "readable", p.is_file(), str(p))

    def write(self) -> bool:
        with self.path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.columns); writer.writeheader(); writer.writerows(self.rows)
        return all(row["status"] == "PASS" for row in self.rows)
