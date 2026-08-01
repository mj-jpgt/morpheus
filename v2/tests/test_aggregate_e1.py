import json
from pathlib import Path

import pandas as pd
import pytest

from morpheus.v2.calibra.aggregate_e1 import aggregate


def _run(path: Path, rank: float, count: float = 0.0) -> None:
    path.mkdir()
    pd.DataFrame([
        {"arm": "before", "metric": "direction_count_above_floor", "value": 3.0},
        {"arm": "paired", "metric": "delta_effective_rank", "value": rank},
        {"arm": "paired", "metric": "delta_direction_count_above_floor", "value": count},
        {"arm": "paired", "metric": "delta_information_density", "value": -0.1},
    ]).to_csv(path / "task_rows.csv", index=False)
    (path / "gate_summary.json").write_text(json.dumps({"gates_pass": True}))


def test_aggregate_requires_three_valid_consistent_runs(tmp_path: Path) -> None:
    runs = [tmp_path / f"seed{seed}" for seed in (42, 43, 44)]
    for run in runs: _run(run, 10.0)
    _, summary = aggregate(runs)
    assert summary["supports_rank_without_detectable_information"]
    _run(tmp_path / "fourth", 2.0)
    with pytest.raises(ValueError, match="exactly three"):
        aggregate(runs + [tmp_path / "fourth"])


def test_aggregate_rejects_count_gain_beyond_equivalence_margin(tmp_path: Path) -> None:
    runs = [tmp_path / f"seed{seed}" for seed in (42, 43, 44)]
    for run in runs: _run(run, 10.0, count=1.0)
    _, summary = aggregate(runs, equivalence_margin=.10)
    assert not summary["supports_rank_without_detectable_information"]
