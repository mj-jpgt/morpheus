"""Aggregate three independently gated E1 runs without re-reading representations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED = ("delta_effective_rank", "delta_direction_count_above_floor", "delta_information_density")


def _metric(rows: pd.DataFrame, name: str) -> float:
    found = rows.loc[(rows.arm == "paired") & (rows.metric == name), "value"]
    if len(found) != 1 or not np.isfinite(float(found.iloc[0])):
        raise ValueError(f"missing or non-finite paired metric {name}")
    return float(found.iloc[0])


def aggregate(paths: list[Path], *, equivalence_margin: float = 0.10) -> tuple[pd.DataFrame, dict]:
    if len(paths) != 3:
        raise ValueError("E1 headline aggregation requires exactly three seed runs")
    records = []
    for path in paths:
        gate = json.loads((path / "gate_summary.json").read_text(encoding="utf-8"))
        if not gate.get("gates_pass", False):
            raise ValueError(f"E1 run is not gate-valid: {path}")
        rows = pd.read_csv(path / "task_rows.csv")
        before = rows.loc[(rows.arm == "before") & (rows.metric == "direction_count_above_floor"), "value"]
        if len(before) != 1 or not np.isfinite(float(before.iloc[0])):
            raise ValueError(f"missing before direction count: {path}")
        record = {"run": str(path), "before_direction_count": float(before.iloc[0])}
        record.update({metric: _metric(rows, metric) for metric in REQUIRED})
        records.append(record)
    frame = pd.DataFrame(records)
    rank_positive = bool((frame.delta_effective_rank > 0).all())
    # A zero pre-intervention count is still meaningful: the post-arm must
    # remain zero to meet the stated equivalence condition.
    allowed = equivalence_margin * np.maximum(frame.before_direction_count.to_numpy(float), 1.0)
    count_equivalent = bool((np.abs(frame.delta_direction_count_above_floor.to_numpy(float)) <= allowed).all())
    summary = {
        "experiment": "E1_rank_vs_information_three_seed", "n_seeds": 3,
        "equivalence_margin_relative_to_before_count": equivalence_margin,
        "all_seeds_rank_increase": rank_positive,
        "all_seeds_direction_count_equivalent": count_equivalent,
        "supports_rank_without_detectable_information": bool(rank_positive and count_equivalent),
        "mean": {key: float(frame[key].mean()) for key in frame.columns if key not in {"run"}},
        "seed_rows": records,
    }
    return frame, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs=3, required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--equivalence-margin", type=float, default=0.10)
    args = parser.parse_args()
    if not 0 <= args.equivalence_margin <= 1:
        raise ValueError("equivalence margin must be in [0,1]")
    frame, summary = aggregate(args.runs, equivalence_margin=args.equivalence_margin)
    args.output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output / "e1_seed_summary.csv", index=False)
    (args.output / "e1_claim_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (args.output / "gate_summary.json").write_text(json.dumps({"experiment": summary["experiment"],
        "gates_pass": True, "claim_supported": summary["supports_rank_without_detectable_information"]}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
