"""Track 2 analysis: turn the induced-correlation sweep tables into the reported law.

Reads the CSVs written by ``induced_correlation_sweep`` and emits (a) the design
x n table, (b) the scaling-law fit with residuals for every declared predictor,
(c) the matched-rank falsifier contrast, (d) the estimator-robustness table, and
(e) both floors as functions of (design rank, n). Findings are written to the
GateLedger with ``observe`` -- they are scientific outcomes and must never decide
whether a run is valid.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .gates import GateLedger
from .induced_correlation_sweep import fit_scaling_law

_PREDICTORS = ("predicted_induced_correlation_p1", "predicted_induced_correlation_p2",
               "predicted_induced_correlation_p3")


def _agg(frame: pd.DataFrame, keys: list[str], value: str = "induced_correlation_median") -> pd.DataFrame:
    """Median over seeds, with the seed spread kept -- a law quoted without its
    seed-to-seed spread cannot be checked by anyone else."""
    grouped = frame.groupby(keys, dropna=False)
    spec = {"n_seeds": (value, "size"), "induced": (value, "median"),
            "induced_min": (value, "min"), "induced_max": (value, "max"),
            "k_columns": ("n_confound_columns", "median"), "k_eff": ("k_eff", "median"),
            "k_eff_shared": ("k_eff_shared", "median"),
            "r2_x": ("design_r2_x_median", "median"), "r2_y": ("design_r2_y_median", "median"),
            "null_scale": ("baseline_null_scale", "median"),
            "gate_threshold": ("baseline_gate_threshold", "median"),
            "closed_form_err": ("closed_form_max_abs_error", "max")}
    out = grouped.agg(**{k: v for k, v in spec.items() if v[0] in frame.columns}).reset_index()
    for key in _PREDICTORS:
        if key in frame.columns:
            out[key.replace("predicted_induced_correlation_", "pred_")] = grouped[key].median().to_numpy()
    return out


def analyse(main: pd.DataFrame, knobs: pd.DataFrame | None, floors: pd.DataFrame | None) -> dict:
    real = main[main.design_mode == "real"].copy()
    synthetic = main[main.design_mode != "real"].copy()
    report: dict = {}

    report["design_by_n"] = _agg(main, ["design", "design_mode", "n_requested"]).to_dict("records")
    report["design_marginal"] = _agg(main, ["design", "design_mode"]).to_dict("records")
    report["n_marginal_anchor_design"] = _agg(
        main[main.design == "cancer_tss_pool10"], ["n_requested"]).to_dict("records")

    # --- T2.4 ---------------------------------------------------------------
    rows_real = real.to_dict("records")
    report["scaling_law_real_designs"] = fit_scaling_law(rows_real)
    report["scaling_law_all_designs"] = fit_scaling_law(main.to_dict("records"))

    # n-dependence at FIXED design: the single number that decides P0 vs the
    # derivation. P0 (|r| ~ k/n) needs a slope near -1; the derivation needs ~0.
    per_design = []
    for name, block in real.groupby("design"):
        block = block[block.induced_correlation_median > 0]
        if block.n_requested.nunique() < 3:
            continue
        slope = float(np.polyfit(np.log(block.n_requested.to_numpy(dtype=float)),
                                 np.log(block.induced_correlation_median.to_numpy()), 1)[0])
        per_design.append({"design": name, "n_exponent": slope,
                           "k_eff": float(block.k_eff.median()),
                           "induced_at_min_n": float(block.loc[block.n_requested.idxmin(), "induced_correlation_median"]),
                           "induced_at_max_n": float(block.loc[block.n_requested.idxmax(), "induced_correlation_median"])})
    report["n_exponent_by_design"] = per_design
    if per_design:
        report["n_exponent_median_over_designs"] = float(np.median([d["n_exponent"] for d in per_design]))

    # rank-dependence at FIXED n, real designs only
    per_n = []
    for n_value, block in real.groupby("n_requested"):
        block = block[(block.induced_correlation_median > 0) & (block.k_eff > 0)]
        if block.k_eff.nunique() < 3:
            continue
        slope = float(np.polyfit(np.log(block.k_eff.to_numpy()),
                                 np.log(block.induced_correlation_median.to_numpy()), 1)[0])
        per_n.append({"n": int(n_value), "k_eff_exponent": slope,
                      "k_eff_min": float(block.k_eff.min()), "k_eff_max": float(block.k_eff.max())})
    report["k_eff_exponent_by_n"] = per_n
    if per_n:
        report["k_eff_exponent_median_over_n"] = float(np.median([d["k_eff_exponent"] for d in per_n]))

    # --- the falsifier ------------------------------------------------------
    contrast = []
    for n_value in sorted(main.n_requested.unique()):
        block = main[main.n_requested == n_value]
        real_value = block[block.design == "cancer_tss_pool10"].induced_correlation_median.median()
        for name in ("permuted_cancer_tss_pool10", "permuted_cancer_tss_pool1", "gaussian_k99", "gaussian_k600"):
            match = block[block.design == name]
            if match.empty:
                continue
            value = float(match.induced_correlation_median.median())
            contrast.append({"n": int(n_value), "arm": name, "induced": value,
                             "real_anchor": float(real_value),
                             "ratio_to_real": float(value / real_value) if real_value else float("nan"),
                             "sampling_null_3_over_sqrt_n": float(match.baseline_null_scale.median())})
    report["matched_rank_falsifier"] = contrast
    if contrast:
        report["matched_rank_falsifier_max_ratio"] = float(max(c["ratio_to_real"] for c in contrast))

    # --- T2.5 ---------------------------------------------------------------
    if knobs is not None and len(knobs):
        table = []
        for (design, splits, alpha), block in knobs.groupby(["design", "residualiser_n_splits", "residualiser_alpha"]):
            table.append({"design": design, "n_splits": int(splits), "alpha": float(alpha),
                          "induced": float(block.induced_correlation_median.median()),
                          "r2_x": float(block.design_r2_x_median.median()),
                          "r2_y": float(block.design_r2_y_median.median()),
                          "closed_form_max_abs_error": float(block.closed_form_max_abs_error.max())})
        report["estimator_robustness"] = table
        summary = {}
        for design, block in knobs.groupby("design"):
            default = block[(block.residualiser_alpha == 1.0) & (block.residualiser_n_splits == 5)]
            if default.empty:
                continue
            reference = float(default.induced_correlation_median.median())
            folds = block[block.residualiser_alpha == 1.0]
            mild = block[(block.residualiser_alpha >= 0.01) & (block.residualiser_alpha <= 10.0)]
            summary[design] = {
                "reference_default": reference,
                "fold_count_span_ratio": float(folds.induced_correlation_median.max()
                                               / max(folds.induced_correlation_median.min(), 1e-12)),
                "mild_alpha_span_ratio": float(mild.induced_correlation_median.max()
                                               / max(mild.induced_correlation_median.min(), 1e-12)),
                "max_closed_form_error": float(block.closed_form_max_abs_error.max()),
            }
        report["estimator_robustness_summary"] = summary

    # --- T2.6 ---------------------------------------------------------------
    if floors is not None and len(floors):
        table = []
        for (design, n_value), block in floors.groupby(["design", "n_requested"]):
            table.append({
                "design": design, "n": int(n_value),
                "k_columns": float(block.n_confound_columns.median()),
                "k_eff": float(block.k_eff.median()),
                "induced": float(block.induced_correlation_median.median()),
                "transmission_floor": float(block.transmission_floor.median()),
                "detection_floor": float(block.detection_floor.median()),
                "detection_floor_finite_fraction": float(np.isfinite(block.detection_floor).mean()),
                "attenuation_slope": float(block.attenuation_slope.median()),
            })
        report["floors"] = table
    return report


def main_cli() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main", required=True)
    parser.add_argument("--knobs", default="")
    parser.add_argument("--floors", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--official-log", default="")
    args = parser.parse_args()

    read = lambda path: pd.read_csv(path) if path and Path(path).exists() else None
    report = analyse(pd.read_csv(args.main), read(args.knobs), read(args.floors))
    output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    (output / "track2_report.json").write_text(json.dumps(report, indent=2, default=float), encoding="utf-8")

    ledger = GateLedger(output, "track2_induced_correlation",
                        official_log=args.official_log or None)
    anchor = [r for r in report["design_by_n"] if r["design"] == "cancer_tss_pool10"]
    for row in anchor:
        ledger.observe("induced_correlation", f"{row['induced']:.4f}",
                       f"k_eff={row['k_eff']:.1f};n={row['n_requested']}",
                       f"cancer+tss pooled@10; seeds={row['n_seeds']}")
    ledger.observe("induced_correlation_n_exponent",
                   f"{report.get('n_exponent_median_over_designs', float('nan')):.3f}",
                   "plan guess -1.0 (k/n); derivation 0.0")
    ledger.observe("induced_correlation_k_eff_exponent",
                   f"{report.get('k_eff_exponent_median_over_n', float('nan')):.3f}",
                   "plan guess +1.0; derivation -0.5")
    ledger.observe("matched_rank_falsifier_max_ratio",
                   f"{report.get('matched_rank_falsifier_max_ratio', float('nan')):.3f}",
                   "<0.25 means a structureless design of matched rank induces nothing")
    ledger.write()
    print(json.dumps({k: report[k] for k in report if not isinstance(report[k], list)}, indent=2, default=float))
    print(f"[done] -> {output}")


if __name__ == "__main__":
    main_cli()
