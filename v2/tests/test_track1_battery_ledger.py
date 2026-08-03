"""Self-tests for the T1.8 battery ledger.

The property under test is not "does it write rows" -- it is the separation the
battery exists to demonstrate: a **losing baseline must not be able to mark the
run FAILED**, because that would make a true negative indistinguishable from a
broken pipeline; while a must-FAIL control that fails to fail must mark it FAILED
loudly. Both directions are asserted.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from morpheus.v2.calibra.gates import GateLedger
from morpheus.v2.calibra.track1_battery_ledger import assemble


def _certificate_rows(method: str, state: str, joint: float, null_p95: float, breaching: float) -> list[dict]:
    common = dict(method=method, representation_state=state, task="confound_certificate",
                  target="tissue_source_site", note="scored")
    return [dict(common, metric=m, value=v) for m, v in
            (("joint_lda_balanced_accuracy", joint), ("joint_null_p95", null_p95),
             ("n_breaching_axes", breaching), ("chance_rate", 0.0118))]


def _write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def _root(tmp_path: Path, *, joint: float, null_p95: float, breaching: float) -> Path:
    _write_rows(tmp_path / "track1" / "certificate_raw" / "task_rows.csv",
                _certificate_rows("art", "wsi_biology", joint, null_p95, breaching))
    return tmp_path


def _run(root: Path, tmp_path: Path):
    ledger = GateLedger(tmp_path / "out", "test_battery", official_log=tmp_path / "GATE_LOG.md")
    assemble(root, ledger)
    return ledger, ledger.write()


def test_a_site_leak_fails_the_run(tmp_path):
    """The must-FAIL control that did not fail must be disqualifying."""
    ledger, verdict = _run(_root(tmp_path, joint=0.36, null_p95=0.15, breaching=17.0), tmp_path)
    assert verdict is False
    failed = [r for r in ledger.rows if r["status"] == "FAIL"]
    assert any("T1.3_site_certificate_raw" in r["gate"] for r in failed)


def test_a_clean_certificate_passes(tmp_path):
    ledger, verdict = _run(_root(tmp_path, joint=0.011, null_p95=0.052, breaching=0.0), tmp_path)
    certificate = [r for r in ledger.rows if "T1.3_site_certificate_raw::" in r["gate"]]
    assert certificate and all(r["status"] == "PASS" for r in certificate)


def _baseline_block(root: Path, block: str, value: float) -> None:
    _write_rows(root / "track1" / f"calibra_{block}" / "task_rows.csv",
                [dict(method="art", representation_state="wsi_biology", task="calibra",
                      target="", metric=m, value=v, note="")
                 for m, v in (("adjusted_top_cca", value), ("heldout_top_cca", value),
                              ("baseline_recovered_median", 0.1), ("detection_floor", 0.3))])


def test_a_losing_baseline_is_an_observation_and_cannot_move_the_verdict(tmp_path):
    """THE structural test. PBS losing badly to a random dictionary is a finding, not a defect.

    The comparison is made against the SAME fixture without the baseline blocks, so
    the assertion is about the baseline's contribution and not about whichever other
    controls happen to be present.
    """
    bare = _root(tmp_path / "bare", joint=0.011, null_p95=0.052, breaching=0.0)
    _, bare_verdict = _run(bare, tmp_path / "bare")

    withbase = _root(tmp_path / "withbase", joint=0.011, null_p95=0.052, breaching=0.0)
    _baseline_block(withbase, "pbs", 0.40)          # ours
    _baseline_block(withbase, "randdict", 0.95)     # the baseline, winning by a mile
    ledger, verdict = _run(withbase, tmp_path / "withbase")

    baseline_rows = [r for r in ledger.rows if "T1.2_baseline_block::" in r["gate"]]
    assert baseline_rows, "the baseline table was not written at all"
    assert all(r["status"] == "OBSERVED" for r in baseline_rows), \
        "a baseline comparison was registered as a pass/fail gate"
    assert verdict == bare_verdict, \
        "adding a baseline we lose to changed the run verdict; a true negative is now " \
        "indistinguishable from a broken pipeline"
    assert not any(r["status"] == "FAIL" for r in baseline_rows)


def test_baseline_observations_are_excluded_from_the_verdict_even_when_all_gates_pass(tmp_path):
    """The same property in the all-clear case, where nothing else can mask it."""
    root = _root(tmp_path, joint=0.011, null_p95=0.052, breaching=0.0)
    _baseline_block(root, "randdict", 0.99)
    ledger = GateLedger(tmp_path / "out", "test_battery", official_log=tmp_path / "GATE_LOG.md")
    # only the rows this test is about, so "everything else is missing" cannot decide it
    for row in [r for r in _run(root, tmp_path)[0].rows
                if "T1.3_site_certificate_raw::" in r["gate"] or "T1.2_baseline_block::" in r["gate"]]:
        ledger.rows.append(row)
    assert ledger.write() is True
    assert any(r["status"] == "OBSERVED" for r in ledger.rows)


def test_a_missing_control_is_recorded_as_a_failure_not_omitted(tmp_path):
    """A missing row counts as FAIL; silence is not a pass."""
    ledger, verdict = _run(tmp_path, tmp_path)          # nothing on disk at all
    gates = {r["gate"]: r for r in ledger.rows if r["status"] != "OBSERVED"}
    for expected in ("T1.3_site_certificate_raw", "T1.4_random_gene_sets",
                     "T1.6_modality_shuffled_pairing", "T1.7a_rna_positive_control",
                     "T1.5_gene_label_shuffle"):
        assert expected in gates, f"{expected} vanished instead of failing"
        assert gates[expected]["status"] == "FAIL"
        assert str(gates[expected]["note"]).startswith("inadmissible_")
    assert verdict is False


def test_the_official_log_is_append_only(tmp_path):
    root = _root(tmp_path, joint=0.011, null_p95=0.052, breaching=0.0)
    _run(root, tmp_path)
    first = (tmp_path / "GATE_LOG.md").read_text(encoding="utf-8").splitlines()
    _run(root, tmp_path)
    second = (tmp_path / "GATE_LOG.md").read_text(encoding="utf-8").splitlines()
    assert len(second) > len(first)
    assert second[:len(first)] == first, "the ledger rewrote history instead of appending"
    assert sum(1 for line in second if line.startswith("timestamp_utc")) == 1


def test_written_rows_carry_the_full_schema(tmp_path):
    root = _root(tmp_path, joint=0.36, null_p95=0.15, breaching=17.0)
    ledger, _ = _run(root, tmp_path)
    with (Path(ledger.path)).open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows and set(rows[0]) == set(GateLedger.columns)
    assert all(r["experiment"] == "test_battery" for r in rows)
