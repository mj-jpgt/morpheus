"""Fail-closed preflight for E1/E2 and the PBS feasibility branch."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.gates import GateLedger


def _manifest(raw: np.lib.npyio.NpzFile) -> dict:
    if "manifest_json" not in raw.files:
        return {}
    value = raw["manifest_json"]
    if np.asarray(value).shape == ():
        value = np.asarray(value).item()
    try:
        return json.loads(str(value))
    except (json.JSONDecodeError, TypeError):
        return {"_invalid": True}


def _manifest_summary(value: dict) -> dict:
    """Keep provenance useful without duplicating patient-level arrays into a log."""
    if not isinstance(value, dict):
        return {"_invalid": True}
    allowed = {
        "artifact_version", "seed", "epochs", "learning_rate", "lr", "token_budget",
        "objective_profile", "architecture", "model_config", "split_digest", "cohort_digest",
        "anchor_teacher_source", "anchor_teacher_sha256", "trained_states",
    }
    result = {key: value[key] for key in sorted(allowed & set(value))}
    result["manifest_keys"] = sorted(value)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--raw-rna", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--official-gate-log", type=Path, default=None)
    args = parser.parse_args()
    output = args.output; output.mkdir(parents=True, exist_ok=True)
    ledger = GateLedger(output, "E1_E2_PREFLIGHT", official_log=args.official_gate_log)
    for name, path in [("G0.3_targets", args.targets), ("G1_raw_rna", args.raw_rna)]:
        ledger.artifact(name, path)
    target = np.load(args.targets, allow_pickle=True)
    required_target = {"patient_ids", "scores", "target_names"}
    ledger.add("G1.1_target_schema", sorted(required_target - set(target.files)), "all required fields", required_target <= set(target.files))
    target_ids = np.asarray([str(x) for x in target["patient_ids"]]) if required_target <= set(target.files) else np.asarray([])
    ledger.add("G1.2_target_unique_patients", int(len(target_ids) - len(np.unique(target_ids))), "0 duplicates", len(target_ids) == len(np.unique(target_ids)))
    reference: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
    artifact_manifest: dict[str, dict] = {}
    for path in args.artifacts:
        ledger.artifact("G0.3_artifact", path)
        raw = np.load(path, allow_pickle=True)
        required = {"patient_ids", "cancers", "split", "trained_states"}
        ledger.add("G0.1_artifact_schema", f"{path.name}:{sorted(required - set(raw.files))}", "all required fields", required <= set(raw.files))
        if required <= set(raw.files):
            values = tuple(np.asarray(raw[name]).astype(str) for name in ("patient_ids", "cancers", "split"))
            artifact_manifest[path.name] = _manifest_summary(_manifest(raw))
            if reference is None:
                reference = values
            else:
                ledger.add("G1.3_artifact_alignment", path.name, "identical patient/cancer/split rows",
                           all(np.array_equal(left, right) for left, right in zip(reference, values)))
            patients = values[0]
            ledger.add("G1.4_artifact_unique_patients", f"{path.name}:{len(patients) - len(np.unique(patients))}", "0 duplicates", len(patients) == len(np.unique(patients)))
            declared = {str(x) for x in raw["trained_states"]}
            present = {name for name in declared if name in raw.files}
            ledger.add("G2.0_declared_states_present", f"{path.name}:{sorted(declared - present)}", "all declared states present", declared == present)
        else:
            artifact_manifest[path.name] = {"_missing_schema": sorted(required - set(raw.files))}
    if reference is not None:
        overlap = sum(patient in set(target_ids) for patient in reference[0])
        ledger.add("G1.5_target_coverage", f"{overlap}/{len(reference[0])}", ">=50 paired", overlap >= 50)
    gates_pass = ledger.write()
    (output / "preflight.json").write_text(json.dumps({"gates_pass": gates_pass, "artifacts": artifact_manifest,
                                                          "n_target_patients": int(len(target_ids))}, indent=2))
    pd.DataFrame(ledger.rows).to_csv(output / "task_rows.csv", index=False)
    if not gates_pass:
        raise SystemExit("preflight gate failure")


if __name__ == "__main__":
    main()
