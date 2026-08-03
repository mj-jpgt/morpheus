"""T1.8 -- assemble the negative-control battery into one append-only ledger.

THE STRUCTURAL POINT THIS FILE EXISTS TO MAKE, and it is a P1 exhibit in its own
right: **must-FAIL and must-PASS controls are validity conditions and go through
``GateLedger.add``; must-beat baseline comparisons are scientific outcomes and go
through ``GateLedger.observe``.** Wiring "the random dictionary beat PBS" into a
pass/fail gate would mark the run FAILED and quarantine its output at exactly the
moment the science came back "no" -- making a true negative indistinguishable
from a broken pipeline. ``GateLedger.write`` deliberately excludes OBSERVED rows
from its verdict for this reason (``v2/calibra/gates.py``).

A missing row counts as FAIL. Controls whose inputs are absent are written with
``value = NaN`` and ``note = "inadmissible_<code>"`` rather than omitted, so an
unrun control is visible instead of silent.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .gates import GateLedger

__all__ = ["assemble"]

_NAN = float("nan")


def _load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def _rows_frame(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path) if path.is_file() else None


def _value(frame: pd.DataFrame, method: str, state: str, metric: str) -> float:
    hit = frame[(frame.method == method) & (frame.representation_state == state) & (frame.metric == metric)]
    return float(hit.value.iloc[0]) if len(hit) else _NAN


def assemble(root: Path, ledger: GateLedger) -> dict:
    """Write every battery row. Returns a small summary for the caller to print."""
    summary: dict[str, object] = {}

    # --- T1.3 confound certificate, both arms -------------------------------
    for arm, folder in (("raw", "certificate_raw"), ("adjusted", "certificate_adjusted")):
        frame = _rows_frame(root / "track1" / folder / "task_rows.csv")
        if frame is None:
            ledger.add(f"T1.3_site_certificate_{arm}", _NAN, "joint<=null_p95", False,
                       f"inadmissible_missing_{folder}")
            continue
        for (method, state), _ in frame.groupby(["method", "representation_state"]):
            joint = _value(frame, method, state, "joint_lda_balanced_accuracy")
            null_p95 = _value(frame, method, state, "joint_null_p95")
            breaching = _value(frame, method, state, "n_breaching_axes")
            chance = _value(frame, method, state, "chance_rate")
            passed = bool(np.isfinite(joint) and np.isfinite(null_p95) and joint <= null_p95 and breaching == 0)
            ledger.add(f"T1.3_site_certificate_{arm}::{method}::{state}", joint,
                       f"<= joint_null_p95 ({null_p95:.4f}) and 0 breaching axes", passed,
                       f"chance={chance:.4f}; breaching_axes={breaching:.0f}; "
                       f"must-FAIL control: the axes must be UNABLE to predict TSS")
        summary[f"certificate_{arm}"] = "written"

    # --- T1.4 / T1.6 / T1.7(a,c) from the main CALIBRA run -------------------
    gates = _load_json(root / "track1" / "calibra_frozen_rna" / "calibra_gates.json")
    frame = _rows_frame(root / "track1" / "calibra_frozen_rna" / "task_rows.csv")
    if gates is None or frame is None:
        ledger.add("T1.4_random_gene_sets", _NAN, "median<floor and <=5% exceed", False,
                   "inadmissible_missing_calibra_frozen_rna")
        ledger.add("T1.6_modality_shuffled_pairing", _NAN, "perm p at stated resolution", False,
                   "inadmissible_missing_calibra_frozen_rna")
        ledger.add("T1.7a_rna_positive_control", _NAN, "passes in the same run", False,
                   "inadmissible_missing_calibra_frozen_rna")
    else:
        for key, verdict in gates.get("random_control_verdicts", {}).items():
            ledger.add(f"T1.4_random_gene_sets::{key}", verdict["median"],
                       f"median < detection_floor ({verdict['detection_floor']}) and "
                       f"exceedance_fraction <= {verdict['exceedance_ceiling']}",
                       bool(verdict["passed"]),
                       f"statistic={verdict['statistic']}; n_controls={verdict['n_controls']}; "
                       f"exceedances={verdict['n_exceedances']}; must-FAIL control")
            fitted = verdict.get("fitted_direction", {})
            if fitted:
                # SCIENTIFIC OUTCOME, not a validity condition: the fitted-direction
                # statistic is not on the floor's scale, so grading it against the
                # floor would be a category error. It is the number a reader cares
                # about, so it is recorded as an observation.
                ratio = (fitted["median"] / fitted["real_target_median"]
                         if fitted.get("real_target_median") else _NAN)
                ledger.observe(f"T1.4_random_vs_real_fitted_direction::{key}", ratio,
                               "well below 1 if the channel is pathway-specific",
                               f"random_control_median={fitted['median']:.4f}; "
                               f"real_target_median={fitted['real_target_median']:.4f}; "
                               f"NOT on the detection-floor scale")

        # T1.6: the permutation null IS the modality-shuffled pairing distribution.
        for (method, state), _ in frame[frame.task == "calibra"].groupby(["method", "representation_state"]):
            perm_p = _value(frame, method, state, "permutation_p")
            observed = _value(frame, method, state, "adjusted_top_cca")
            null_median = _value(frame, method, state, "permutation_null_median")
            null_p95 = _value(frame, method, state, "permutation_null_p95")
            passed = bool(np.isfinite(observed) and np.isfinite(null_p95) and observed > null_p95)
            ledger.add(f"T1.6_modality_shuffled_pairing::{method}::{state}", null_median,
                       f"shuffled pairing collapses to the null; observed ({observed:.4f}) > null_p95 ({null_p95:.4f})",
                       passed,
                       f"permutation_p={perm_p:.6f} at resolution 1/2001; null_median is the "
                       f"16-component capacity floor, NOT zero")
            if state.startswith("rna_"):
                ledger.add(f"T1.7a_rna_positive_control::{method}::{state}", observed,
                           "must clear its own pairing null in the same run", passed,
                           "circular by construction; HANDOFF_GATES governing rule")
            transmission = _value(frame, method, state, "transmission_floor")
            detection = _value(frame, method, state, "detection_floor")
            ledger.add(f"T1.7c_spike_recovery::{method}::{state}", transmission,
                       "transmission_floor finite (a planted spike is recovered)",
                       bool(np.isfinite(transmission)),
                       f"detection_floor={detection}; floor_scale=targeted_single_direction")
        summary["calibra_frozen_rna"] = "written"

    # --- T1.5 gene-label shuffle --------------------------------------------
    shuffle = _load_json(root / "track1" / "gene_label_shuffle" / "gene_label_shuffle_control.json")
    if shuffle is None:
        ledger.add("T1.5_gene_label_shuffle", _NAN, "subspace persists and attribution collapses", False,
                   "inadmissible_missing_gene_label_shuffle")
    else:
        for result in shuffle["results"]:
            name = Path(result.get("shuffled_targets", "?")).name
            sub, attr = result.get("subspace", {}), result.get("attribution", {})
            ledger.add(f"T1.5_attribution_collapses::{name}", attr.get("median_abs_spearman_by_axis_index", _NAN),
                       "<= 0.05", bool(attr.get("attribution_collapsed")),
                       "TRUE BY CONSTRUCTION (rows permuted after the fit); build-integrity check only")
            ledger.add(f"T1.5_subspace_persists::{name}", sub.get("heldout_top_cca_shuffled", _NAN),
                       f"inside CI95 of the unshuffled value {sub.get('ci95_true')}",
                       bool(sub.get("subspace_persists")),
                       f"true={sub.get('heldout_top_cca_true'):.4f}; paired difference CI="
                       f"{sub.get('ci95_paired_difference')}")
            # The scientific reading of a PASS, which the gate cannot express.
            ledger.observe(f"T1.5_shuffled_minus_true::{name}",
                           -sub.get("paired_difference_true_minus_shuffled", _NAN),
                           "negative if the fitted dictionary beats a spectrum-matched random projection",
                           "a pass here means any spectrum-matched projection is equally legible")
        summary["gene_label_shuffle"] = "written"

    # --- T1.7(b) known covariate --------------------------------------------
    for name in ("covariate_er", "covariate_pr"):
        frame = _rows_frame(root / "track1" / name / "task_rows.csv")
        payload = _load_json(root / "track1" / name / "known_covariate_control.json")
        if frame is None and payload is None:
            ledger.add(f"T1.7b_known_covariate::{name}", _NAN, "CI overlaps the pre-registered interval",
                       False, f"inadmissible_missing_{name}")
            continue
        # The covariate module keys its payload by "method::state::adjustment".
        records = [(key, value) for key, value in (payload or {}).items()
                   if isinstance(value, dict) and value.get("status") == "scored"]
        for key, record in records:
            ledger.add(f"T1.7b_known_covariate::{name}::{key}",
                       record.get("within_cancer_auroc", _NAN),
                       f"CI95 overlaps pre-registered {record.get('expected_interval')}",
                       bool(record.get("passed")),
                       f"ci={record.get('within_cancer_auroc_ci95')}; null_p95={record.get('null_p95')}; "
                       f"covariate={record.get('covariate')}; pre-registered before the run")
        summary[name] = f"{len(records)} rows"

    # --- T1.1/T1.2 must-beat baseline table: OBSERVE, never gate -------------
    reference = None
    for block in ("pbs", "randdict", "pca", "pbsshuf1", "frozen_rna"):
        frame = _rows_frame(root / "track1" / f"calibra_{block}" / "task_rows.csv")
        if frame is None:
            ledger.observe(f"T1.2_baseline_block::{block}", _NAN, "scored on the same instrument",
                           f"unavailable_missing_calibra_{block}")
            continue
        for (method, state), _ in frame[frame.task == "calibra"].groupby(["method", "representation_state"]):
            value = _value(frame, method, state, "adjusted_top_cca")
            ledger.observe(f"T1.2_baseline_block::{block}::{method}::{state}", value,
                           "compared against the PBS block on the identical instrument",
                           f"heldout={_value(frame, method, state, 'heldout_top_cca'):.4f}; "
                           f"induced_baseline={_value(frame, method, state, 'baseline_recovered_median'):.4f}; "
                           f"detection_floor={_value(frame, method, state, 'detection_floor')}; "
                           f"SCIENTIFIC OUTCOME -- never a pass/fail gate")
            if block == "pbs" and state == "wsi_biology":
                reference = value
    summary["baseline_reference_pbs_wsi_biology"] = reference
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--official-log", required=True)
    parser.add_argument("--experiment", default="P1_track1_negative_control_battery")
    args = parser.parse_args()

    ledger = GateLedger(args.output, args.experiment, official_log=args.official_log)
    summary = assemble(Path(args.evidence_root), ledger)
    passed = ledger.write()
    gates = [r for r in ledger.rows if r["status"] != "OBSERVED"]
    failed = [r for r in gates if r["status"] == "FAIL"]
    print(json.dumps({"summary": summary, "n_rows": len(ledger.rows), "n_gates": len(gates),
                      "n_observations": len(ledger.rows) - len(gates),
                      "n_failed_gates": len(failed),
                      "failed_gates": [r["gate"] for r in failed],
                      "verdict_all_gates_pass": passed}, indent=2))


if __name__ == "__main__":
    main()
