import json
from pathlib import Path

from morpheus.v2.research.rebase.evidence_cycle import Stage, _completed, _fingerprint, _json_digest
from morpheus.v2.research.rebase.evidence_preflight import _manifest_summary


def test_stage_requires_complete_contract(tmp_path: Path) -> None:
    stage = Stage.parse({"name": "one", "command": ["python", "-V"], "inputs": ["input.txt"], "outputs": ["out.json"]})
    (tmp_path / "input.txt").write_text("input")
    fingerprint = _fingerprint(stage, tmp_path, "config", "code")
    directory = tmp_path / "stage"; directory.mkdir()
    (directory / "out.json").write_text("{}")
    (directory / "SUCCESS.json").write_text(json.dumps({"fingerprint": fingerprint, "outputs": ["out.json"]}))
    assert _completed(directory, fingerprint)
    (directory / "out.json").unlink()
    assert not _completed(directory, fingerprint)


def test_stage_parse_rejects_missing_contract() -> None:
    try:
        Stage.parse({"name": "broken"})
    except ValueError as error:
        assert "missing" in str(error)
    else:
        raise AssertionError("missing stage contract was accepted")


def test_json_digest_is_key_order_invariant() -> None:
    assert _json_digest({"a": 1, "b": 2}) == _json_digest({"b": 2, "a": 1})


def test_preflight_manifest_summary_excludes_patient_rows() -> None:
    summary = _manifest_summary({"seed": 42, "epochs": 10, "patient_ids": ["P1"], "cancers": ["X"]})
    assert summary["seed"] == 42
    assert "patient_ids" not in summary and "cancers" not in summary
