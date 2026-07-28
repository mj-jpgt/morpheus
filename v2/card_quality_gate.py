"""Predeclared quantitative gate for scientific hypothesis-card rendering."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


@dataclass(frozen=True)
class CardQualityGate:
    """Minimum evidence required before an LLM may verbalize a card.

    Values are deliberately conservative defaults.  They are recorded in the
    gate artifact and may be made stricter, but never relaxed after outer-test
    results are seen.
    """
    programme_ndcg_min: float = 0.20
    direction_accuracy_min: float = 0.60
    confidence_brier_max: float = 0.25
    selective_risk_max: float = 0.35

    def assess(self, metrics: Mapping[str, Any]) -> dict[str, Any]:
        values = {
            "programme_ndcg_at_k": float(metrics.get("programme_ndcg_at_k", float("nan"))),
            "direction_accuracy": float(metrics.get("direction_accuracy", float("nan"))),
            "confidence_brier": float(metrics.get("confidence_brier", float("nan"))),
            "selective_risk": float(metrics.get("selective_risk", float("nan"))),
        }
        checks = {
            "programme_ranking": values["programme_ndcg_at_k"] >= self.programme_ndcg_min,
            "direction": values["direction_accuracy"] >= self.direction_accuracy_min,
            "calibration": values["confidence_brier"] <= self.confidence_brier_max,
            "selective_risk": values["selective_risk"] <= self.selective_risk_max,
        }
        return {"schema_version": 1, "passed": bool(all(checks.values())), "checks": checks,
                "metrics": values, "thresholds": self.__dict__}


def load_quality_gate(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("passed"), bool):
        raise ValueError("quality gate must be a JSON object with boolean 'passed'")
    return payload
