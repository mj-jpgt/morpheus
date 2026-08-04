"""Does Track 1's random-control verdict survive the sign fix?

``run_calibra.random_direction_column_correlation`` used to return
``np.median(scores, axis=0)`` over per-draw **signed** correlations, and
``grade_random_controls`` then compared that against a strictly positive
``detection_floor``. With ``u`` drawn at random the sign of each draw is random, so
the signed median collapses towards zero for **every** column whatever it carries.
A control that always reads ~0 always passes: the T1.4 negative control was passing
partly by construction rather than on merit, and that had never been investigated.

The fix (schema 2) takes the median of ``|correlation|``. That moves the controls
**up**, i.e. it makes our own negative control harder to pass. This script measures
by how much, on the shipped states, and re-grades the T1.4 verdict both ways. It is
run in full knowledge that the honest direction of the fix is against us, and both
verdicts are reported whatever they say.

The real (non-control) target block is scored through the identical statistic, because
"the controls are below the floor" is only informative next to what the real targets
read on the same scale.

Every statistic is imported from ``v2/calibra``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,  # noqa: E402
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.run_calibra import (grade_random_controls,  # noqa: E402
                                             random_direction_column_correlation)


def signed_random_direction_column_correlation(x_residual, y_residual, *, n_draws, seed):
    """The SCHEMA-1 statistic, reproduced here and NOWHERE else in the library.

    Kept in this script rather than left in ``run_calibra`` so that the defective
    convention cannot be called by accident, while the comparison that shows why it
    is defective is still reproducible. It is the same loop as the fixed function
    with the ``np.abs`` removed, which is exactly the one-token difference under test.
    """
    x_residual = np.asarray(x_residual, dtype=np.float64)
    y_residual = np.asarray(y_residual, dtype=np.float64)
    rng = np.random.default_rng(seed)
    scores = np.empty((n_draws, y_residual.shape[1]))
    y_std = (y_residual - y_residual.mean(axis=0)) / np.maximum(y_residual.std(axis=0), 1e-12)
    for draw in range(n_draws):
        u = rng.normal(size=x_residual.shape[1])
        u /= np.linalg.norm(u)
        projected = x_residual @ u
        projected = (projected - projected.mean()) / max(projected.std(), 1e-12)
        scores[draw] = (projected[:, None] * y_std).mean(axis=0)
    return np.median(scores, axis=0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--task-rows", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n-draws", type=int, default=40)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--exceedance-ceiling", type=float, default=0.05)
    args = parser.parse_args()

    rows = pd.read_csv(args.task_rows)
    shipped = rows[rows.task == "calibra"].pivot_table(
        index=["method", "representation_state"], columns="metric", values="value")

    targets = np.load(args.targets, allow_pickle=True)
    names = np.asarray(targets["target_names"]).astype(str)
    is_control = np.char.startswith(names, "RANDOM_CONTROL__")
    all_scores = np.asarray(targets["scores"], dtype=np.float64)
    target_index = {str(pid): i for i, pid in enumerate(targets["patient_ids"])}

    results = []
    for path in args.artifacts:
        path = Path(path)
        artifact = np.load(path, allow_pickle=True)
        patient_ids = artifact["patient_ids"]
        cancers = artifact["cancers"]
        split = artifact["split"]
        aligned = np.asarray([target_index.get(str(pid), -1) for pid in patient_ids])
        mask = (split == "test") & (aligned >= 0)
        tss, _ = pooled_tissue_source_site(patient_ids[mask], min_site_count=args.min_site_count)
        design = confound_design(pd.DataFrame({"cancer": cancers[mask], "tss": tss}),
                                 ["cancer", "tss"])
        y = all_scores[aligned[mask]][:, ~is_control]
        control = all_scores[aligned[mask]][:, is_control]
        y_res = cross_fitted_residuals(y, design, seed=args.seed)
        control_res = cross_fitted_residuals(control, design, seed=args.seed)

        for state in sorted({str(s) for s in artifact["trained_states"]}):
            key = (path.stem, state)
            if key not in shipped.index:
                continue
            floor = float(shipped.loc[key]["detection_floor"])
            x = np.asarray(artifact[state], dtype=np.float64)[mask]
            x_res = cross_fitted_residuals(x, design, seed=args.seed)

            record = {"method": path.stem, "state": state, "detection_floor": floor,
                      "n_controls": int(is_control.sum()), "n_real_targets": int((~is_control).sum())}
            for label, fn in (("schema1_signed", signed_random_direction_column_correlation),
                              ("schema2_abs", random_direction_column_correlation)):
                controls = fn(x_res, control_res, n_draws=args.n_draws, seed=args.seed)
                real = fn(x_res, y_res, n_draws=args.n_draws, seed=args.seed)
                verdict = grade_random_controls(controls, names[is_control], detection_floor=floor,
                                                exceedance_ceiling=args.exceedance_ceiling)
                record[label] = {
                    "control_median": verdict["median"], "control_p95": verdict["p95"],
                    "control_max": verdict["max"],
                    "n_exceedances": verdict["n_exceedances"],
                    "exceedance_fraction": verdict["exceedance_fraction"],
                    "passed": bool(verdict["passed"]),
                    "real_target_median": float(np.nanmedian(real)),
                    "real_target_max": float(np.nanmax(real)),
                    # The number that says whether the control is informative: if the
                    # real targets read the same as the controls, "controls below the
                    # floor" is a statement about a weak readout, not about controls.
                    "real_over_control_median": float(np.nanmedian(real) / verdict["median"])
                    if verdict["median"] > 1e-12 else float("nan"),
                }
            record["verdict_changed_by_sign_fix"] = bool(
                record["schema1_signed"]["passed"] != record["schema2_abs"]["passed"])
            results.append(record)
            print(f"[{path.stem}::{state}] floor={floor} "
                  f"signed: med={record['schema1_signed']['control_median']:.4f} "
                  f"exc={record['schema1_signed']['exceedance_fraction']:.3f} "
                  f"pass={record['schema1_signed']['passed']} | "
                  f"abs: med={record['schema2_abs']['control_median']:.4f} "
                  f"exc={record['schema2_abs']['exceedance_fraction']:.3f} "
                  f"pass={record['schema2_abs']['passed']} | "
                  f"changed={record['verdict_changed_by_sign_fix']}", flush=True)

    Path(args.out).write_text(json.dumps(results, indent=2))
    changed = sum(r["verdict_changed_by_sign_fix"] for r in results)
    print(f"[written] {args.out} ({len(results)} states, {changed} verdicts changed)", flush=True)


if __name__ == "__main__":
    main()
