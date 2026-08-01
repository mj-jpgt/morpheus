"""Predeclared inner-validation promotion gate for V2.1 objective profiles.

``identity_only`` is an anchored preservation control, not a competing V2
method.  In particular, it is expected to preserve a teacher-trained
retrieval geometry and therefore must not be allowed to win a biology profile
selection simply by retaining that geometry.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", required=True, help="v21_evaluation task_rows.csv from development diagnostics")
    parser.add_argument("--panel", required=True, help="curated_panel_manifest.json")
    parser.add_argument("--candidate-map", required=True, help="JSON mapping profile -> artifact stem")
    parser.add_argument("--teacher-method", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--retrieval-tolerance", type=float, default=0.02)
    parser.add_argument(
        "--allow-control-selection",
        action="store_true",
        help=("Diagnostic-only escape hatch. Never use for a primary V2 claim; "
              "the output remains explicitly labelled as a control selection."),
    )
    args = parser.parse_args()
    rows = pd.read_csv(args.rows)
    panel = json.loads(Path(args.panel).read_text(encoding="utf-8"))
    candidates = json.loads(Path(args.candidate_map).read_text(encoding="utf-8"))
    if "full" not in candidates:
        raise ValueError("candidate-map must contain the declared primary 'full' V2 profile")
    eligible_targets = set(panel["eligible_targets"])
    control_targets = {value for value in rows.get("target", pd.Series(dtype=str)).dropna().astype(str) if value.startswith("RANDOM_CONTROL__")}
    teacher = rows.loc[(rows.method == args.teacher_method) & (rows.task == "retrieval") & (rows.metric == "recall_at_10"), "value"]
    if len(teacher) != 1:
        raise ValueError("teacher must have exactly one inner-validation recall_at_10")
    teacher_r10 = float(teacher.iloc[0])
    summary: dict[str, dict[str, object]] = {}
    for profile, method in candidates.items():
        pearson = rows.loc[(rows.method == method) & (rows.task == "molecular_prompting") & (rows.metric == "pearson")]
        biological = pearson.loc[pearson.target.astype(str).isin(eligible_targets), "value"].to_numpy(float)
        controls = pearson.loc[pearson.target.astype(str).isin(control_targets), "value"].to_numpy(float)
        r10 = rows.loc[(rows.method == method) & (rows.task == "retrieval") & (rows.metric == "recall_at_10"), "value"].to_numpy(float)
        bio_mean = float(np.nanmean(biological)) if len(biological) else float("nan")
        control_mean = float(np.nanmean(controls)) if len(controls) else float("nan")
        retrieval = float(r10[0]) if len(r10) == 1 else float("nan")
        delta = bio_mean - control_mean
        eligible = bool(np.isfinite(delta) and np.isfinite(retrieval) and delta > 0 and retrieval >= teacher_r10 - args.retrieval_tolerance)
        role = ("primary_method" if profile == "full" else
                "anchored_preservation_control" if profile == "identity_only" else
                "ablation_control")
        summary[profile] = {"method": method, "role": role,
                            "eligible_for_primary_claim": profile == "full",
                            "biological_mean_pearson": bio_mean, "control_mean_pearson": control_mean,
                            "biological_minus_control": delta, "retrieval_r10": retrieval, "teacher_r10": teacher_r10,
                            "eligible": eligible}
    # The primary method may only select an epoch/profile within ``full``.
    # Controls are always reported but cannot substitute for a failed full
    # method.  This prevents an anchor-equivalent state from being mislabeled
    # as a full V2 success.
    full_eligible = bool(summary["full"]["eligible"])
    control_promoted = [name for name, values in summary.items()
                        if name != "full" and bool(values["eligible"])]
    selected = "full" if full_eligible else None
    if args.allow_control_selection and not full_eligible and control_promoted:
        selected = max(control_promoted, key=lambda name: (
            float(summary[name]["biological_minus_control"]),
            float(summary[name]["biological_mean_pearson"]),
        ))
    payload = {
        "schema_version": 2,
        "promoted": selected is not None,
        "selected_profile": selected,
        "primary_profile": "full",
        "primary_claim_eligible": bool(selected == "full" and full_eligible),
        "control_selection_enabled": bool(args.allow_control_selection),
        "control_profiles_passing_diagnostic_gate": control_promoted,
        "selection": summary,
        "rule": "full V2 only: positive curated-biological-minus-matched-control mean Pearson and R@10 within tolerance of seed-matched teacher",
        "claim_safeguard": (
            "identity_only and programme_only are controls/ablations. They cannot replace full V2 for a primary claim; "
            "a failed full profile is reported as a negative result."
        ),
    }
    Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    if selected is None:
        raise SystemExit("The declared full V2 profile did not meet the predeclared promotion gate")


if __name__ == "__main__":
    main()
