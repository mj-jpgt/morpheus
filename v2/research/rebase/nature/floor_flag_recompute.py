"""Recompute the floor comparison correctly on the shipped CALIBRA artifacts.

The shipped ``observed_above_floor`` is broken in two independent ways, both
demonstrated on constructed cases with known answers by
``floor_flag_diagnostic.py``:

  1. it compares a **signed** statistic against a positive threshold, so a
     channel of the right magnitude and the wrong sign reads FAIL; and
  2. for a multi-column representation the statistic is a correlation along a
     **random** direction pair, which carries essentially none of the channel and
     is ~0 whatever the channel's strength.

This script rebuilds x, y and the design for every state of a CALIBRA run
exactly as ``run_calibra.py`` does, reproduces the shipped ``heldout_top_cca``
as a check on the reconstruction, and grades each state against its own
``detection_floor`` on the corrected comparator, with the guards predeclared in
``NOTEBOOK_ENTRIES/PREDECLARED_observed_above_floor_20260804T1843Z.md``:

  G1  the corrected comparator under DESTROYED cross-modal pairing (y rows
      permuted within cancer strata) must not clear the floor;
  G3  it must exceed the confound-induced level-0 baseline;
  G4  the verdict must hold across seeds.

Every statistic is imported. Nothing is computed inline.
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
from morpheus.v2.calibra.spectral import heldout_top_cca, paired_absolute_correlation  # noqa: E402

try:
    from joblib import Parallel, delayed
except ImportError:                                   # pragma: no cover
    Parallel = None


def build_state(artifact, state: str, target_index: dict, scores: np.ndarray,
                min_site_count: int):
    """x, y, design, strata for one state -- the run_calibra.py construction, verbatim."""
    patient_ids = artifact["patient_ids"]
    cancers = artifact["cancers"]
    split = artifact["split"]
    aligned = np.asarray([target_index.get(str(pid), -1) for pid in patient_ids])
    mask = (split == "test") & (aligned >= 0)

    x = np.asarray(artifact[state], dtype=np.float64)[mask]
    y = scores[aligned[mask]]
    tss, frequent = pooled_tissue_source_site(patient_ids[mask], min_site_count=min_site_count)
    frame = pd.DataFrame({"cancer": cancers[mask], "tss": tss})
    design = confound_design(frame, ["cancer", "tss"])
    return x, y, design, np.asarray(cancers[mask]), len(frequent)


def random_direction_magnitude(x_res, y_res, *, n_draws: int, seed: int) -> float:
    """Median |corr| along RANDOM direction pairs -- the shipped comparator with its
    sign defect repaired and nothing else. Reported to show that repairing the sign
    alone does NOT rescue the flag: a random direction through a 256-column
    representation sees the induced baseline, not the channel."""
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_draws):
        u = rng.normal(size=x_res.shape[1]); u /= np.linalg.norm(u)
        v = rng.normal(size=y_res.shape[1]); v /= np.linalg.norm(v)
        values.append(paired_absolute_correlation(x_res @ u, y_res @ v))
    return float(np.nanmedian(values))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--task-rows", required=True, help="task_rows.csv of the run being graded")
    parser.add_argument("--out", required=True)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--seeds", type=int, default=3, help="G4: how many seeds")
    parser.add_argument("--n-permutations", type=int, default=100, help="G1")
    parser.add_argument("--n-random-draws", type=int, default=40)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    rows = pd.read_csv(args.task_rows)
    shipped = rows[rows.task == "calibra"].pivot_table(
        index=["method", "representation_state"], columns="metric", values="value")

    targets = np.load(args.targets, allow_pickle=True)
    target_names = np.asarray(targets["target_names"]).astype(str)
    keep = ~np.char.startswith(target_names, "RANDOM_CONTROL__")
    scores = np.asarray(targets["scores"], dtype=np.float64)[:, keep]
    target_index = {str(pid): i for i, pid in enumerate(targets["patient_ids"])}

    results = []
    for path in args.artifacts:
        path = Path(path)
        artifact = np.load(path, allow_pickle=True)
        for state in sorted({str(s) for s in artifact["trained_states"]}):
            key = (path.stem, state)
            if key not in shipped.index:
                print(f"[skip] {key} not in task_rows", flush=True)
                continue
            row = shipped.loc[key]
            floor = float(row["detection_floor"])
            x, y, design, strata, n_sites = build_state(artifact, state, target_index, scores,
                                                        args.min_site_count)

            x_res = cross_fitted_residuals(x, design, seed=args.seed)
            y_res = cross_fitted_residuals(y, design, seed=args.seed)

            # Reconstruction check: must reproduce the shipped heldout_top_cca.
            reproduced = heldout_top_cca(x_res, y_res, n_components=args.n_components,
                                         seed=args.seed)

            # G4 -- seeds. The seed moves BOTH the residualiser folds and the CCA split.
            per_seed = []
            for seed in range(args.seed, args.seed + args.seeds):
                xr = cross_fitted_residuals(x, design, seed=seed)
                yr = cross_fitted_residuals(y, design, seed=seed)
                per_seed.append(heldout_top_cca(xr, yr, n_components=args.n_components, seed=seed))

            # G1 -- destroyed cross-modal pairing, permuted WITHIN cancer strata.
            permute_rng = np.random.default_rng(args.seed)
            orders = []
            for _ in range(args.n_permutations):
                order = np.arange(len(y))
                for level in np.unique(strata):
                    idx = np.flatnonzero(strata == level)
                    order[idx] = permute_rng.permutation(idx)
                orders.append(order)

            def _one(order):
                permuted = cross_fitted_residuals(y[order], design, seed=args.seed)
                return heldout_top_cca(x_res, permuted, n_components=args.n_components,
                                       seed=args.seed)

            if args.n_jobs != 1 and Parallel is not None:
                null = np.asarray(Parallel(n_jobs=args.n_jobs, prefer="processes")(
                    delayed(_one)(o) for o in orders), dtype=np.float64)
            else:
                null = np.asarray([_one(o) for o in orders], dtype=np.float64)

            magnitude_random = random_direction_magnitude(
                x_res, y_res, n_draws=args.n_random_draws, seed=args.seed)

            record = {
                "method": path.stem, "state": state, "n_patients": int(len(x)),
                "n_design_columns": int(design.shape[1]), "n_sites_kept": int(n_sites),
                "detection_floor": floor,
                "transmission_floor": float(row["transmission_floor"]),
                "confound_induced_baseline": float(row["baseline_recovered_median"]),
                "shipped_observed_matched_direction": float(row["observed_matched_direction"]),
                "shipped_observed_above_floor": bool(row["observed_above_floor"]),
                "shipped_heldout_top_cca": float(row["heldout_top_cca"]),
                "reproduced_heldout_top_cca": float(reproduced),
                "reconstruction_delta": float(abs(reproduced - float(row["heldout_top_cca"]))),
                "sign_fix_only_random_direction_magnitude": magnitude_random,
                "sign_fix_only_clears_floor": bool(magnitude_random > floor),
                "corrected_flag": bool(np.isfinite(floor) and reproduced > floor),
                "corrected_ratio_over_floor": float(reproduced / floor) if floor > 0 else float("nan"),
                "G4_per_seed_heldout_top_cca": [float(v) for v in per_seed],
                "G4_min_over_seeds": float(np.min(per_seed)),
                "G4_all_seeds_clear": bool(np.all(np.asarray(per_seed) > floor)),
                "G1_null_median": float(np.median(null)),
                "G1_null_p95": float(np.percentile(null, 95)),
                "G1_null_max": float(np.max(null)),
                "G1_null_max_clears_floor": bool(np.max(null) > floor),
                "G1_permutation_p": float((int(np.sum(null >= reproduced)) + 1)
                                          / (args.n_permutations + 1)),
                "G3_clears_induced_baseline": bool(reproduced > float(row["baseline_recovered_median"])),
                "n_permutations": int(args.n_permutations),
            }
            results.append(record)
            print(f"[{path.stem}::{state}] floor={floor} heldout={reproduced:.4f} "
                  f"(shipped {float(row['heldout_top_cca']):.4f}, d={record['reconstruction_delta']:.2e}) "
                  f"CORRECTED={record['corrected_flag']} ratio={record['corrected_ratio_over_floor']:.2f} "
                  f"G1max={record['G1_null_max']:.4f} G1clears={record['G1_null_max_clears_floor']} "
                  f"G4min={record['G4_min_over_seeds']:.4f} "
                  f"signfixonly={magnitude_random:.4f}", flush=True)

    Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"[written] {args.out} ({len(results)} states)", flush=True)


if __name__ == "__main__":
    main()
