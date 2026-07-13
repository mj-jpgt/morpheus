"""Aggregate shared discovery-support task results across seeds/method variants."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


SEED_RE = re.compile(r"_seed(\d+)$")


def split_method_seed(method: str) -> tuple[str, int | None]:
    match = SEED_RE.search(method)
    if match:
        return method[: match.start()], int(match.group(1))
    if method.startswith("v2_seed_"):
        return "morpheus_v2", int(method.rsplit("_", 1)[-1])
    return method, None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    frames = [pd.read_csv(path) for path in args.inputs]
    raw = pd.concat(frames, ignore_index=True)
    parsed = raw["method"].astype(str).map(split_method_seed)
    raw["method_family"] = [x[0] for x in parsed]
    raw["seed"] = [x[1] for x in parsed]

    group_cols = ["method_family", "task", "metric"]
    summary = (
        raw.groupby(group_cols, dropna=False)["value"]
        .agg(["count", "mean", "std", "min", "max"])
        .reset_index()
        .sort_values(["task", "metric", "mean"], ascending=[True, True, False])
    )

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    raw.to_csv(out / "task_suite_with_families.csv", index=False)
    summary.to_csv(out / "task_suite_summary.csv", index=False)
    (out / "task_suite_summary.json").write_text(summary.to_json(orient="records", indent=2))
    print(out / "task_suite_summary.csv")


if __name__ == "__main__":
    main()
