"""Small, append-only health-gate ledger used by CLD experiments."""
from __future__ import annotations

import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class GateLedger:
    columns = ("timestamp_utc", "experiment", "gate", "value", "threshold", "status", "note")

    def __init__(self, output: str | Path, experiment: str, *, official_log: str | Path | None = None):
        self.output = Path(output); self.output.mkdir(parents=True, exist_ok=True)
        self.path = self.output / "gate_rows.csv"; self.experiment = experiment
        self.official_log = Path(official_log) if official_log else None
        self.rows: list[dict[str, Any]] = []

    def add(self, gate: str, value: Any, threshold: Any, passed: bool, note: str = "") -> None:
        self.rows.append({"timestamp_utc": datetime.now(timezone.utc).isoformat(), "experiment": self.experiment,
                          "gate": gate, "value": value, "threshold": threshold,
                          "status": "PASS" if passed else "FAIL", "note": note})

    def observe(self, gate: str, value: Any, expectation: Any, note: str = "") -> None:
        """Record a SCIENTIFIC OUTCOME, which must never decide whether a run is valid.

        A gate answers "is this run trustworthy"; an observation answers "what did
        the experiment find". Registering a finding as a gate makes a *true negative*
        indistinguishable from a broken pipeline -- the run gets marked FAIL and its
        output quarantined precisely when the science came back "no", which is the
        one answer we most need to be able to read. Observations are written to the
        same ledger for provenance but are excluded from ``write()``'s verdict.
        """
        self.rows.append({"timestamp_utc": datetime.now(timezone.utc).isoformat(), "experiment": self.experiment,
                          "gate": gate, "value": value, "threshold": expectation,
                          "status": "OBSERVED", "note": note})

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
        # The run-local CSV is convenient for machine parsing.  The experiment
        # protocol also requires an append-only human-readable project ledger;
        # write the exact same rows there rather than making a second, lossy
        # summary by hand.
        if self.official_log:
            self.official_log.parent.mkdir(parents=True, exist_ok=True)
            new_file = not self.official_log.exists()
            with self.official_log.open("a", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=self.columns)
                if new_file:
                    writer.writeheader()
                writer.writerows(self.rows)
        # FAIL is disqualifying; OBSERVED rows are findings and carry no verdict.
        return not any(row["status"] == "FAIL" for row in self.rows)
