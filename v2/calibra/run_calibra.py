"""CALIBRA Phase 1 runner: recovery curve + calibrated channel measurement.

Runs the instrument over every declared representation state in one or more
frozen artifacts, against the RNA-derived molecular targets, under a fixed
confound adjustment (cancer type + tissue source site, plus any extra numeric
covariates supplied).

Emits the repo-standard flat task-row schema, including the mandatory
``metric="status" / value=NaN / note="unavailable_..."`` convention so that a
state we could not measure is visible rather than silently absent.

Example
-------
    PYTHONPATH=$WS python -m morpheus.v2.calibra.run_calibra \
        --artifacts .../diagnostic_full_seed42.npz \
        --targets   .../frozen_rna_targets.npz \
        --output    runs/calibra_phase1 --n-draws 10
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .calibration import permutation_null, spike_recovery_curve
from .residualise import confound_design, cross_fitted_residuals
from .spectral import cca_spectrum, effective_rank, heldout_top_cca


def _row(**kwargs) -> dict:
    base = {"method": "", "representation_state": "", "task": "", "target": "",
            "metric": "", "value": np.nan, "note": ""}
    base.update(kwargs)
    return base


def _unavailable(method: str, state: str, task: str, reason: str) -> dict:
    return _row(method=method, representation_state=state, task=task,
                metric="status", value=np.nan, note=f"unavailable_{reason}")


def load_artifact(path: Path):
    raw = np.load(path, allow_pickle=True)
    declared = set()
    if "trained_states" in raw.files:
        declared = {str(s) for s in raw["trained_states"]}
    return raw, declared


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True, help="frozen_rna_targets.npz")
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition", default="test", choices=("test", "val", "all"))
    parser.add_argument("--levels", default="0.0,0.01,0.02,0.05,0.10,0.20,0.40")
    parser.add_argument("--n-draws", type=int, default=10)
    parser.add_argument("--n-components", type=int, default=32)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-permutations", type=int, default=50)
    parser.add_argument("--min-site-count", type=int, default=10,
                        help="pool TSS sites with fewer patients into OTHER")
    parser.add_argument("--target-group", default="", help="restrict to one target_group")
    parser.add_argument("--n-jobs", type=int, default=1,
                        help="parallel workers; affects wall-clock only, never the numbers")
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    levels = tuple(float(v) for v in args.levels.split(","))

    targets_raw = np.load(args.targets, allow_pickle=True)
    target_ids = np.asarray([str(p) for p in targets_raw["patient_ids"]])
    target_names = np.asarray([str(t) for t in targets_raw["target_names"]])
    target_groups = np.asarray([str(g) for g in targets_raw["target_groups"]])
    scores = np.asarray(targets_raw["scores"], dtype=np.float64)
    keep_targets = np.ones(len(target_names), dtype=bool)
    if args.target_group:
        keep_targets = target_groups == args.target_group
    # Exclude the random-control columns from the molecular block itself; they are
    # a null for a different question and would dilute the channel estimate.
    keep_targets &= ~np.char.startswith(target_names, "RANDOM_CONTROL__")
    scores = scores[:, keep_targets]
    target_index = {pid: i for i, pid in enumerate(target_ids)}

    rows: list[dict] = []
    summaries: dict[str, dict] = {}

    for artifact_path in args.artifacts:
        path = Path(artifact_path)
        method = path.stem
        raw, declared = load_artifact(path)
        if not declared:
            rows.append(_unavailable(method, "", "calibra", "no_declared_states"))
            continue
        patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
        cancers = np.asarray([str(c) for c in raw["cancers"]])
        split = np.asarray([str(s) for s in raw["split"]])
        mask = np.ones(len(patient_ids), dtype=bool) if args.partition == "all" else (split == args.partition)
        # Align to the molecular targets.
        aligned = np.asarray([target_index.get(pid, -1) for pid in patient_ids])
        mask &= aligned >= 0
        if mask.sum() < 50:
            rows.append(_unavailable(method, "", "calibra", f"insufficient_paired_patients_{int(mask.sum())}"))
            continue

        y = scores[aligned[mask]]
        tss_raw = np.asarray([pid.split("-")[1] if len(pid.split("-")) > 1 else "NA"
                              for pid in patient_ids[mask]])
        # Pool rare sites. ~600 TSS codes with 145 singletons means a rare-site dummy
        # is effectively a per-patient indicator, and residualising it would delete
        # that patient's signal entirely rather than adjust for a site effect.
        unique, counts = np.unique(tss_raw, return_counts=True)
        frequent = {u for u, c in zip(unique, counts) if c >= args.min_site_count}
        tss = np.asarray([s if s in frequent else "OTHER" for s in tss_raw])
        frame = pd.DataFrame({"cancer": cancers[mask], "tss": tss})
        design = confound_design(frame, ["cancer", "tss"])
        rows.append(_row(method=method, task="calibra", metric="n_patients", value=float(mask.sum()),
                         note=f"partition={args.partition}"))
        rows.append(_row(method=method, task="calibra", metric="n_confound_columns",
                         value=float(design.shape[1]), note="cancer+tss_pooled"))
        rows.append(_row(method=method, task="calibra", metric="n_distinct_sites_kept",
                         value=float(len(frequent)), note=f"min_site_count={args.min_site_count}"))

        for state in sorted(declared):
            if state not in raw.files:
                rows.append(_unavailable(method, state, "calibra", "state_declared_but_absent"))
                continue
            x = np.asarray(raw[state], dtype=np.float64)[mask]
            if x.ndim != 2 or not np.isfinite(x).all():
                rows.append(_unavailable(method, state, "calibra", "state_not_finite_2d"))
                continue

            unadjusted = cca_spectrum(x, y, n_components=args.n_components)
            result = spike_recovery_curve(x, y, design, levels=levels, n_draws=args.n_draws,
                                          n_components=args.n_components, seed=args.seed,
                                          n_jobs=args.n_jobs)
            summary = result.summary()
            # The chance level. Without this the adjusted top-CCA is uninterpretable:
            # it is a multivariate maximum and is inflated by capacity alone.
            null = permutation_null(x, y, design, strata=frame["cancer"].to_numpy(),
                                    n_permutations=args.n_permutations,
                                    n_components=args.n_components, seed=args.seed,
                                    n_jobs=args.n_jobs)
            summary.update(null)
            summaries[f"{method}::{state}"] = summary

            emit = {
                "effective_rank": effective_rank(x),
                "unadjusted_top_cca": float(unadjusted[0]) if unadjusted.size else np.nan,
                "adjusted_top_cca": summary["observed"],
                "detection_floor": summary["detection_floor"],
                "attenuation_slope": summary["attenuation_slope"],
                "null_reference_p90": summary["null_reference_p90"],
                "observed_above_floor": float(summary["observed_above_floor"]),
                # Targeted-readout diagnostics. baseline_* is the REAL-DATA GATE:
                # with a direction-matched readout the level-0 value must sit near
                # the sampling null, not near 1. A baseline near 1 means the readout
                # is reading ambient structure and every floor below it is void.
                "baseline_recovered_median": summary["baseline_recovered_median"],
                "baseline_is_null_like": float(summary["baseline_is_null_like"]),
                "baseline_null_scale": summary["baseline_null_scale"],
                "observed_matched_direction": summary["observed_matched_direction"],
                "permutation_null_median": null["null_median"],
                "permutation_null_p95": null["null_p95"],
                "excess_over_null_median": null["excess_over_null_median"],
                "permutation_p": null["permutation_p"],
                "heldout_top_cca": heldout_top_cca(
                    cross_fitted_residuals(x, design, seed=args.seed),
                    cross_fitted_residuals(y, design, seed=args.seed),
                    n_components=args.n_components, seed=args.seed),
            }
            for metric, value in emit.items():
                rows.append(_row(method=method, representation_state=state, task="calibra",
                                 metric=metric, value=float(value),
                                 note=f"n_draws={args.n_draws};k={args.n_components}"))
            for level, med in zip(summary["levels"], summary["recovered_median"]):
                rows.append(_row(method=method, representation_state=state, task="calibra_recovery",
                                 target=f"r_true={level:g}", metric="recovered_median",
                                 value=float(med), note="spike_recovery_curve"))
            gate = "OK" if summary["baseline_is_null_like"] else "GATE-FAIL(baseline_not_null_like)"
            print(f"[{method}::{state}] adj_cca={summary['observed']:.4f} "
                  f"floor={summary['detection_floor']} "
                  f"atten={summary['attenuation_slope']:.3f} "
                  f"base={summary['baseline_recovered_median']:.4f} {gate}", flush=True)

    frame = pd.DataFrame(rows)
    frame.to_csv(output / "task_rows.csv", index=False)
    (output / "calibra_summary.json").write_text(json.dumps(summaries, indent=2))
    (output / "calibra_protocol.json").write_text(json.dumps({
        "levels": list(levels), "n_draws": args.n_draws, "n_components": args.n_components,
        "seed": args.seed, "partition": args.partition, "confounds": ["cancer", "tss"],
        "n_permutations": args.n_permutations, "permutation_strata": "cancer",
        "n_targets": int(scores.shape[1]), "target_group": args.target_group or "all_non_control",
        "recovery_fraction": 0.8, "n_jobs": args.n_jobs,
        "readout": "targeted_single_direction",
        "readout_note": ("recovery scored on the planted (u,v) axis, not a top-CCA maximum; "
                         "floor units are single-direction correlation and are NOT comparable "
                         "to the multivariate adjusted/held-out top-CCA"),
    }, indent=2))
    print(f"[done] {len(frame)} rows -> {output}", flush=True)


if __name__ == "__main__":
    main()
