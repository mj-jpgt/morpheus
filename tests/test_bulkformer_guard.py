import json
from pathlib import Path

import pytest

from morpheus.src.encoders.encode_bulkformer import ensure_bulkformer_configured


def test_bulkformer_guard_requires_zenodo_record_id(tmp_path: Path):
    cfg = {
        "project_root": ".",
        "paths": {},
        "foundation_models": {"rna_primary": {"zenodo_record_id": ""}},
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    with pytest.raises(RuntimeError) as exc:
        ensure_bulkformer_configured(path)
    assert "zenodo_record_id is blank" in str(exc.value)


def test_bulkformer_guard_requires_manifest(tmp_path: Path):
    cfg = {
        "project_root": str(tmp_path),
        "paths": {"bulkformer_cache_dir": "cache"},
        "foundation_models": {"rna_primary": {"zenodo_record_id": "15744294"}},
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    with pytest.raises(RuntimeError) as exc:
        ensure_bulkformer_configured(path)
    assert "download_bulkformer_zenodo.py" in str(exc.value)
