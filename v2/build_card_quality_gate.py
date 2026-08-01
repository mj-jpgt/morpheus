"""Materialize the predeclared quantitative gate for hypothesis cards."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .card_quality_gate import CardQualityGate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True, help="JSON object from structured-card evaluation")
    parser.add_argument("--output", required=True)
    parser.add_argument("--programme-ndcg-min", type=float, default=0.20)
    parser.add_argument("--direction-accuracy-min", type=float, default=0.60)
    parser.add_argument("--confidence-brier-max", type=float, default=0.25)
    parser.add_argument("--selective-risk-max", type=float, default=0.35)
    args = parser.parse_args()
    metrics = json.loads(Path(args.metrics).read_text(encoding="utf-8"))
    gate = CardQualityGate(args.programme_ndcg_min, args.direction_accuracy_min,
                           args.confidence_brier_max, args.selective_risk_max).assess(metrics)
    Path(args.output).write_text(json.dumps(gate, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
