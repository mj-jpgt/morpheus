import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.card_quality_gate import CardQualityGate
from morpheus.v2.canonical_registry import RegistryError, audit_survival_alignment
from morpheus.v2.hypothesis_cards import abstain_hypothesis_cards


def test_survival_alignment_reports_split_coverage_and_rejects_duplicates():
    outcomes = pd.DataFrame({
        "patient_id": ["TCGA-AA-0001", "TCGA-BB-0002"],
        "os_time_days": [10.0, 20.0], "os_event": [1.0, 0.0],
        "pfi_time_days": [8.0, 18.0], "pfi_event": [1.0, 0.0],
    })
    audit = audit_survival_alignment(outcomes, ["TCGA-AA-0001", "TCGA-BB-0002", "TCGA-CC-0003"], ["train", "test", "test"])
    assert audit["matched_outcome_count"] == 2
    assert audit["split_alignment"]["test"]["unmatched_patients"] == ["TCGA-CC-0003"]
    with pytest.raises(RegistryError):
        audit_survival_alignment(pd.concat([outcomes, outcomes.iloc[[0]]]), ["TCGA-AA-0001"], ["test"])


def test_card_evaluator_emits_selective_risk_so_the_gate_can_pass():
    """Regression: the quality gate requires a selective_risk metric, but
    evaluate_hypothesis_cards never produced one -> nan -> the gate could NEVER
    pass and every discovery card was silently forced to abstention. The evaluator
    must now emit selective_risk, and a confident+correct card set must clear the gate."""
    import numpy as np
    from morpheus.v2.hypothesis_cards import (PROGRAMMES, evaluate_hypothesis_cards,
                                              generate_hypothesis_cards)
    names = sorted(PROGRAMMES)[:6]
    pids = [f"TCGA-AA-{i:04d}" for i in range(10)]
    scores = np.random.default_rng(0).normal(size=(len(pids), len(names))).astype(float)
    cards = generate_hypothesis_cards(pids, scores, names)
    # Evaluate against the SAME score matrix as targets => predictions are confident
    # and direction-correct => selective risk ~0.
    metrics = evaluate_hypothesis_cards(cards, pids, scores, names)
    assert "selective_risk" in metrics
    assert metrics["selective_risk"] == metrics["selective_risk"]  # finite, not nan
    gate = CardQualityGate().assess(metrics)
    assert gate["passed"], (gate["checks"], metrics)
    assert gate["checks"]["selective_risk"]


def test_failed_card_gate_emits_only_valid_abstentions():
    gate = CardQualityGate().assess({"programme_ndcg_at_k": 0.1, "direction_accuracy": 0.5,
                                     "confidence_brier": 0.5, "selective_risk": 0.8})
    assert not gate["passed"]
    card = abstain_hypothesis_cards(["TCGA-AA-0001"], reason="gate_failed")[0]
    assert card["uncertainty"]["abstain"]
    assert not card["predicted_programmes"]


def test_primary_selector_never_promotes_identity_only(tmp_path):
    rows = pd.DataFrame([
        {"method": "teacher", "task": "retrieval", "metric": "recall_at_10", "value": .50, "target": ""},
        {"method": "full", "task": "retrieval", "metric": "recall_at_10", "value": .50, "target": ""},
        {"method": "identity", "task": "retrieval", "metric": "recall_at_10", "value": .50, "target": ""},
        {"method": "full", "task": "molecular_prompting", "metric": "pearson", "value": .31, "target": "pathway"},
        {"method": "full", "task": "molecular_prompting", "metric": "pearson", "value": .30, "target": "RANDOM_CONTROL__pathway__0"},
        {"method": "identity", "task": "molecular_prompting", "metric": "pearson", "value": .70, "target": "pathway"},
        {"method": "identity", "task": "molecular_prompting", "metric": "pearson", "value": .30, "target": "RANDOM_CONTROL__pathway__0"},
    ])
    csv = tmp_path / "rows.csv"; rows.to_csv(csv, index=False)
    panel = tmp_path / "panel.json"; panel.write_text(json.dumps({"eligible_targets": ["pathway"]}))
    candidates = tmp_path / "candidates.json"; candidates.write_text(json.dumps({"full": "full", "identity_only": "identity"}))
    output = tmp_path / "selection.json"
    subprocess.run([sys.executable, "-m", "morpheus.v2.select_v21_profile", "--rows", str(csv), "--panel", str(panel),
                    "--candidate-map", str(candidates), "--teacher-method", "teacher", "--output", str(output)], check=True)
    payload = json.loads(output.read_text())
    assert payload["selected_profile"] == "full"
    assert payload["control_profiles_passing_diagnostic_gate"] == ["identity_only"]


def test_v22_dag_promotes_atomic_stage_and_resumes(tmp_path):
    source = tmp_path / "source.txt"; source.write_text("input")
    config = tmp_path / "dag.json"
    command = [sys.executable, "-c", "from pathlib import Path; import sys; Path(sys.argv[1]).mkdir(parents=True, exist_ok=True); Path(sys.argv[1], 'done.txt').write_text('ok')", "{stage_dir}"]
    config.write_text(json.dumps({"stages": [{"name": "smoke", "command": command, "inputs": [str(source)], "outputs": ["done.txt"]}]}))
    run_root = tmp_path / "run"
    invocation = [sys.executable, "-m", "morpheus.v2.v22_dag", "--config", str(config), "--run-root", str(run_root)]
    subprocess.run(invocation, check=True)
    assert (run_root / "stages" / "smoke" / "SUCCESS.json").is_file()
    subprocess.run([*invocation, "--resume"], check=True)
    state = json.loads((run_root / "dag_state.json").read_text())
    assert state["complete"] and state["results"][0]["status"] == "resumed"
