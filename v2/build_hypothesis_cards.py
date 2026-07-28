"""Build deterministic WSI-only hypothesis cards from frozen probe predictions."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .card_quality_gate import load_quality_gate
from .hypothesis_cards import PROGRAMMES, abstain_hypothesis_cards, generate_hypothesis_cards


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, help="v21 molecular_predictions.parquet")
    parser.add_argument("--method", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition-patients", default="", help="optional newline-delimited patient allowlist")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--quality-gate", default="", help="required quantitative card-gate JSON for non-abstaining cards")
    args = parser.parse_args()
    frame = pd.read_parquet(args.predictions)
    frame = frame.loc[(frame.method.astype(str) == args.method) & frame.target.astype(str).isin(PROGRAMMES)].copy()
    if args.partition_patients:
        allowed = {line.strip() for line in Path(args.partition_patients).read_text(encoding="utf-8").splitlines() if line.strip()}
        frame = frame.loc[frame.patient_id.astype(str).isin(allowed)]
    if frame.empty:
        raise ValueError("no mechanism-programme predictions were found for the requested method")
    scores = frame.pivot(index="patient_id", columns="target", values="prediction")
    # A row may contain only the programme targets measurable for that patient;
    # the deterministic card builder preserves that missingness rather than
    # inventing an evidence value.
    ids = scores.index.astype(str).tolist()
    gate = load_quality_gate(args.quality_gate) if args.quality_gate else {"passed": False, "reason": "missing_quality_gate"}
    cards = (generate_hypothesis_cards(ids, scores.to_numpy(float), scores.columns.astype(str).tolist(), threshold=args.threshold)
             if gate["passed"] else abstain_hypothesis_cards(ids, reason=str(gate.get("reason", "quantitative_card_gate_failed"))))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(json.dumps(card, sort_keys=True) + "\n" for card in cards), encoding="utf-8")
    output.with_suffix(output.suffix + ".manifest.json").write_text(json.dumps({"method": args.method, "n_cards": len(cards), "programmes": scores.columns.astype(str).tolist(), "threshold": args.threshold, "quality_gate": gate, "quality_gate_path": args.quality_gate or None}, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
