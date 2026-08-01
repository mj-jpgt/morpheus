import json
from pathlib import Path

import pytest

from morpheus.v2.research.rebase.run_e1_e2_analysis import validate_training_manifest


def _record(root: Path, seed: int, arm: str) -> dict:
    artifact = root / f"e1_seed{seed}_{arm}.npz"
    artifact.write_bytes(b"artifact")
    return {"seed": seed, "arm": arm, "artifact": str(artifact),
            "liveness": {"loss_finite": True, "loss_reduced_20_percent": True}}


def test_training_manifest_requires_all_matched_arms(tmp_path: Path):
    rows = [_record(tmp_path, seed, arm) for seed in (42, 43, 44) for arm in ("before", "after")]
    manifest = tmp_path / "training_manifest.json"
    manifest.write_text(json.dumps({"experiment": "E1_matched_training", "seeds": [42, 43, 44], "records": rows}))
    assert validate_training_manifest(manifest)["seeds"] == [42, 43, 44]
    manifest.write_text(json.dumps({"experiment": "E1_matched_training", "seeds": [42, 43, 44], "records": rows[:-1]}))
    with pytest.raises(ValueError, match="incomplete E1 arm set"):
        validate_training_manifest(manifest)
