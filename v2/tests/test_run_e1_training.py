from argparse import Namespace
from pathlib import Path
import json

import pytest

from morpheus.v2.research.rebase.run_e1_training import _read_liveness


def test_liveness_requires_real_loss_reduction(tmp_path: Path) -> None:
    path = tmp_path / "train_metrics.jsonl"
    path.write_text("\n".join(json.dumps({"train_loss": value}) for value in (2.0, 1.5)))
    result = _read_liveness(tmp_path)
    assert result["loss_finite"] and result["loss_reduced_20_percent"]


def test_liveness_rejects_nonfinite_or_short_records(tmp_path: Path) -> None:
    (tmp_path / "train_metrics.jsonl").write_text(json.dumps({"train_loss": float("nan")}))
    result = _read_liveness(tmp_path)
    assert not result["loss_finite"] and not result["loss_reduced_20_percent"]
