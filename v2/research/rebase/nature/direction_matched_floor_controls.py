"""H6 (the blocking refactor control) and the mechanism behind the matched-floor failure.

H6, declared in ``NOTEBOOK_ENTRIES/PREDECLARED_direction_matched_floor_20260804T2230Z.md``
as blocking: the extension adds ``image_direction`` to ``spike_targets``,
``direction_pairs`` to ``spike_recovery_curve``, and lifts the floor rule into
``floors_from_recovery``. Re-running the shipped path -- ``spike_recovery_curve`` with
**no direction arguments at all**, at the shipped protocol -- must reproduce the
shipped ``detection_floor`` and ``transmission_floor`` for all 13 states exactly. If
any shipped floor moves, the refactor changed the instrument and every direction-matched
number is void.

Note what this is NOT. The ``F_rand_real`` cell of ``direction_matched_floor.py``
supplies random pairs drawn by the *script*, so it is a **re-draw** of the
random-direction floor, not a reproduction of the shipped one. It answers a different
and separately interesting question -- how stable is the floor to the direction draw --
and it is not the H6 control. This script is.

The mechanism block measures why a fitted direction pair breaks the spike construction.
``spike_targets`` builds ``a_perp = standardise(a - rho*s)`` with ``rho = corr(s, a)``
taken in **raw** space. For a random pair ``rho`` is near zero and the subtraction is a
no-op; for a **fitted** pair ``rho`` is large by construction, so the subtraction removes
a large raw-space component that residualisation would not have removed, and the level-0
readout lands far from zero. That is a property of the instrument, not of the channel,
and it is what makes the level-0 baselines of the matched cells 2-6x the random ones.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from morpheus.v2.calibra.calibration import spike_recovery_curve, spike_targets  # noqa: E402
from morpheus.v2.calibra.residualise import cross_fitted_residuals  # noqa: E402
from morpheus.v2.calibra.spectral import paired_absolute_correlation  # noqa: E402

from direction_matched_floor import (build_state, fitted_pairs, random_pairs,  # noqa: E402
                                     stratified_order)


def raw_and_residual_alignment(x, y, design, pairs, *, seed):
    """|corr| along each pair, before and after residualisation, plus the level-0 readout.

    ``rho_raw`` is exactly the quantity ``spike_targets`` subtracts. ``rho_residual`` is
    what the readout would see. When the two differ, orthogonalising in raw space
    leaves a residue in residual space -- and ``level0_readout`` measures that residue
    directly, by planting r_true = 0 and reading the planted axis.
    """
    x_res = cross_fitted_residuals(x, design, seed=seed)
    raw, residual, level0 = [], [], []
    for u, v in pairs:
        raw.append(paired_absolute_correlation(x @ u, y @ v))
        y_res = cross_fitted_residuals(y, design, seed=seed)
        residual.append(paired_absolute_correlation(x_res @ u, y_res @ v))
        spiked = spike_targets(x, y, 0.0, rng=np.random.default_rng(0),
                               image_direction=u, molecular_direction=v)
        spiked_res = cross_fitted_residuals(spiked, design, seed=seed)
        level0.append(paired_absolute_correlation(x_res @ u, spiked_res @ v))
    return {"rho_raw_median": float(np.nanmedian(raw)),
            "rho_residual_median": float(np.nanmedian(residual)),
            "level0_readout_abs_median": float(np.nanmedian(level0))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--task-rows", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--levels", default="0.0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.3,0.4,0.5,0.6")
    parser.add_argument("--n-draws", type=int, default=40)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-mechanism-pairs", type=int, default=8)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--skip-h6", action="store_true",
                        help="mechanism block only; H6 itself is unaffected by it")
    args = parser.parse_args()

    levels = tuple(float(v) for v in args.levels.split(","))
    rows = pd.read_csv(args.task_rows)
    shipped = rows[rows.task == "calibra"].pivot_table(
        index=["method", "representation_state"], columns="metric", values="value")

    targets = np.load(args.targets, allow_pickle=True)
    names = np.asarray(targets["target_names"]).astype(str)
    keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
    scores = np.asarray(targets["scores"], dtype=np.float64)[:, keep]
    target_index = {str(pid): i for i, pid in enumerate(targets["patient_ids"])}

    results = []
    for path in args.artifacts:
        path = Path(path)
        artifact = np.load(path, allow_pickle=True)
        for state in sorted({str(s) for s in artifact["trained_states"]}):
            key = (path.stem, state)
            if key not in shipped.index:
                continue
            row = shipped.loc[key]
            x, y, design, strata, _ = build_state(artifact, state, target_index, scores,
                                                  args.min_site_count)

            # H6: the shipped path, no direction arguments, shipped protocol.
            if args.skip_h6:
                summary = {"detection_floor": float("nan"), "transmission_floor": float("nan"),
                           "observed_matched_direction": float("nan"),
                           "baseline_recovered_median": float("nan")}
                detection_ok = transmission_ok = matched_ok = False
            else:
                result = spike_recovery_curve(x, y, design, levels=levels, n_draws=args.n_draws,
                                              n_components=args.n_components, seed=args.seed,
                                              n_jobs=args.n_jobs)
                summary = result.summary()
                detection_ok = bool(np.isclose(summary["detection_floor"],
                                               float(row["detection_floor"]), equal_nan=True))
                transmission_ok = bool(np.isclose(summary["transmission_floor"],
                                                  float(row["transmission_floor"]), equal_nan=True))
                matched_ok = bool(np.isclose(summary["observed_matched_direction"],
                                             float(row["observed_matched_direction"]), atol=1e-9))

            x_res = cross_fitted_residuals(x, design, seed=args.seed)
            y_res = cross_fitted_residuals(y, design, seed=args.seed)
            n = args.n_mechanism_pairs
            mechanism = {
                "fitted": raw_and_residual_alignment(
                    x, y, design,
                    fitted_pairs(x_res, y_res, n_draws=n, n_components=args.n_components,
                                 seed=args.seed), seed=args.seed),
                "random": raw_and_residual_alignment(
                    x, y, design,
                    random_pairs(x.shape[1], y.shape[1], n_draws=n, seed=args.seed),
                    seed=args.seed),
            }
            # The SAME fitted pairs on the pairing-destroyed cohort. This is the cell
            # the primary floors were measured in, so if its level-0 baseline is also
            # inflated the cause must still be a large raw-space rho -- cancer-mediated
            # rather than channel-mediated, since a within-cancer permutation preserves
            # cancer structure on both sides. Measured, not asserted.
            y_perm = y[stratified_order(strata, np.random.default_rng(args.seed))]
            mechanism["fitted_on_permuted_cohort"] = raw_and_residual_alignment(
                x, y_perm, design,
                fitted_pairs(x_res, y_res, n_draws=n, n_components=args.n_components,
                             seed=args.seed), seed=args.seed)
            mechanism["random_on_permuted_cohort"] = raw_and_residual_alignment(
                x, y_perm, design,
                random_pairs(x.shape[1], y.shape[1], n_draws=n, seed=args.seed), seed=args.seed)

            record = {
                "method": path.stem, "state": state,
                "H6_shipped_detection_floor": float(row["detection_floor"]),
                "H6_reproduced_detection_floor": summary["detection_floor"],
                "H6_detection_floor_reproduced": detection_ok,
                "H6_shipped_transmission_floor": float(row["transmission_floor"]),
                "H6_reproduced_transmission_floor": summary["transmission_floor"],
                "H6_transmission_floor_reproduced": transmission_ok,
                "H6_shipped_observed_matched_direction": float(row["observed_matched_direction"]),
                "H6_reproduced_observed_matched_direction": summary["observed_matched_direction"],
                "H6_observed_matched_direction_reproduced": matched_ok,
                "H6_shipped_baseline_recovered_median": float(row["baseline_recovered_median"]),
                "H6_reproduced_baseline_recovered_median": summary["baseline_recovered_median"],
                "mechanism": mechanism,
            }
            results.append(record)
            print(f"[{path.stem}::{state}] H6 floor {row['detection_floor']} -> "
                  f"{summary['detection_floor']} ok={detection_ok} | trans "
                  f"{row['transmission_floor']} -> {summary['transmission_floor']} "
                  f"ok={transmission_ok} | matched ok={matched_ok} | "
                  f"rho_raw fitted={mechanism['fitted']['rho_raw_median']:.3f} "
                  f"random={mechanism['random']['rho_raw_median']:.3f} | "
                  f"level0 fitted={mechanism['fitted']['level0_readout_abs_median']:.3f} "
                  f"random={mechanism['random']['level0_readout_abs_median']:.3f} | PERM "
                  f"rho_raw fitted={mechanism['fitted_on_permuted_cohort']['rho_raw_median']:.3f} "
                  f"random={mechanism['random_on_permuted_cohort']['rho_raw_median']:.3f} "
                  f"level0 fitted={mechanism['fitted_on_permuted_cohort']['level0_readout_abs_median']:.3f} "
                  f"random={mechanism['random_on_permuted_cohort']['level0_readout_abs_median']:.3f}", flush=True)
            Path(args.out).write_text(json.dumps(results, indent=2))

    failed = [r for r in results if not (r["H6_detection_floor_reproduced"]
                                         and r["H6_transmission_floor_reproduced"])]
    Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"[written] {args.out} ({len(results)} states, H6 FAILURES: {len(failed)})", flush=True)


if __name__ == "__main__":
    main()
