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
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .calibration import permutation_null, spike_recovery_curve
from .residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from .spectral import cca_spectrum, effective_rank, heldout_top_cca


def _canonical_tcga_patient(value: object) -> str:
    """Accept a TCGA aliquot/sample ID only by reducing a valid TCGA prefix."""
    text = str(value).strip().upper()
    pieces = text.split("-")
    if len(pieces) < 3 or pieces[0] != "TCGA" or any(not piece for piece in pieces[:3]):
        raise ValueError(f"purity table contains a non-TCGA identifier: {value!r}")
    return "-".join(pieces[:3])


def load_purity_table(path: str, patient_column: str = "", purity_column: str = "", *,
                      units: str = "fraction", source: str = "", reference: str = "") -> tuple[dict[str, float], dict[str, object]]:
    """Read one published/expression-derived purity value per canonical patient."""
    if source not in {"published_consensus", "absolute", "expression_derived"} or not reference.strip():
        raise ValueError("purity source and a concrete method/citation reference are mandatory")
    table = pd.read_parquet(path) if str(path).lower().endswith(".parquet") else pd.read_csv(path, sep=None, engine="python")
    patient_candidates = ([patient_column] if patient_column else []) + ["patient_id", "Patient ID", "bcr_patient_barcode", "sample_id", "Sample"]
    purity_candidates = ([purity_column] if purity_column else []) + ["purity", "tumor_purity", "Tumor_Purity", "ABSOLUTE_PURITY"]
    patient = next((name for name in patient_candidates if name in table.columns), None)
    purity = next((name for name in purity_candidates if name in table.columns), None)
    if patient is None or purity is None:
        raise ValueError("purity table needs an explicit patient/sample ID column and a numeric purity column")
    values = pd.to_numeric(table[purity], errors="coerce")
    if units not in {"fraction", "percent"}:
        raise ValueError("purity units must be fraction or percent")
    upper = 1.0 if units == "fraction" else 100.0
    records: dict[str, float] = {}
    for raw_id, value in zip(table[patient], values):
        if not np.isfinite(value):
            continue
        if value < 0.0 or value > upper:
            raise ValueError(f"purity value {value} is outside declared {units} range")
        canonical = _canonical_tcga_patient(raw_id)
        if canonical in records:
            raise ValueError(f"purity table has duplicate canonical patient {canonical}; choose one consensus value explicitly")
        records[canonical] = float(value / upper)
    if not records:
        raise ValueError("purity table has no finite canonical TCGA values")
    return records, {"path": str(Path(path).resolve()), "sha256": sha256(Path(path).read_bytes()).hexdigest(),
                     "patient_column": patient, "purity_column": purity, "units_input": units,
                     "units_internal": "fraction", "source": source, "reference": reference,
                     "circularity": "expression_derived_from_same_RNA_targets" if source == "expression_derived" else "none",
                     "n_canonical_patients": len(records)}


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


def _channel_measurement(x: np.ndarray, y: np.ndarray, design: np.ndarray, strata: np.ndarray, *,
                         levels: tuple[float, ...], n_draws: int, n_components: int,
                         n_permutations: int, seed: int, n_jobs: int) -> tuple[dict, object]:
    """Run every CALIBRA channel gate through one shared code path.

    D3 compares two confound designs.  A before/after contrast is scientific
    evidence only if both designs receive the same target integrity, recovery,
    pairing-null and held-out checks.  Keeping this operation indivisible makes
    it impossible for a future caller to accidentally report a weaker `before`.
    """
    result = spike_recovery_curve(x, y, design, levels=levels, n_draws=n_draws,
                                  n_components=n_components, seed=seed, n_jobs=n_jobs)
    summary = result.summary()
    null = permutation_null(x, y, design, strata=strata, n_permutations=n_permutations,
                            n_components=n_components, seed=seed, n_jobs=n_jobs)
    summary.update(null)
    summary["permutation_null_non_degenerate"] = bool(
        all(np.isfinite(float(null[key])) for key in ("null_median", "null_p95", "null_max", "permutation_p"))
        and float(null["null_max"]) - float(null["null_median"]) > 1e-12)
    unadjusted = cca_spectrum(x, y, n_components=n_components)
    summary.update({
        "effective_rank": effective_rank(x),
        "unadjusted_top_cca": float(unadjusted[0]) if unadjusted.size else float("nan"),
        "heldout_top_cca": heldout_top_cca(
            cross_fitted_residuals(x, design, seed=seed),
            cross_fitted_residuals(y, design, seed=seed), n_components=n_components, seed=seed),
    })
    return summary, result


_CHANNEL_METRICS = (
    "adjusted_top_cca", "heldout_top_cca", "permutation_null_median", "permutation_null_p95",
    "permutation_p", "excess_over_null_median", "baseline_is_null_like", "detection_floor",
    "transmission_floor", "attenuation_slope", "observed_matched_direction", "permutation_null_non_degenerate",
)


def _channel_value(summary: dict, metric: str) -> float:
    aliases = {"adjusted_top_cca": "observed", "permutation_null_median": "null_median",
               "permutation_null_p95": "null_p95"}
    return float(summary.get(aliases.get(metric, metric), float("nan")))


def validate_evaluated_targets(y: np.ndarray, names: np.ndarray) -> None:
    """G1 target-matrix guard for the exact cohort reaching CALIBRA."""
    value, labels = np.asarray(y, dtype=np.float64), np.asarray(names).astype(str)
    if value.ndim != 2 or value.shape[1] != len(labels) or not np.isfinite(value).all():
        raise ValueError("evaluated RNA targets must be a finite [patient,target] matrix")
    std = np.std(value, axis=0)
    dead = labels[np.where(~np.isfinite(std) | (std < 1e-8))[0][:5]].tolist()
    if dead:
        raise ValueError(f"evaluated RNA targets contain constant columns: {dead}")


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
    parser.add_argument("--purity-table", default="", help="published TCGA consensus/ABSOLUTE or declared expression-derived purity table")
    parser.add_argument("--purity-patient-column", default="")
    parser.add_argument("--purity-column", default="")
    parser.add_argument("--purity-source", default="", choices=("", "published_consensus", "absolute", "expression_derived"))
    parser.add_argument("--purity-units", default="fraction", choices=("fraction", "percent"))
    parser.add_argument("--purity-reference", default="", help="dataset DOI/URL or the declared expression-derived method/version")
    parser.add_argument("--require-rna-positive-control", action="store_true",
                        help="fail the stage unless an rna_* state clears its same-run pairing null")
    parser.add_argument("--require-channel-gates", action="store_true",
                        help="fail the stage if an evaluated state fails CALIBRA's targeted-readout gate")
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    levels = tuple(float(v) for v in args.levels.split(","))
    purity_by_patient, purity_manifest = ({}, {"enabled": False})
    if args.purity_table:
        purity_by_patient, purity_manifest = load_purity_table(
            args.purity_table, args.purity_patient_column, args.purity_column,
            units=args.purity_units, source=args.purity_source, reference=args.purity_reference,
        )
        purity_manifest.update({"enabled": True})

    targets_raw = np.load(args.targets, allow_pickle=True)
    target_ids = np.asarray([str(p) for p in targets_raw["patient_ids"]])
    target_names = np.asarray([str(t) for t in targets_raw["target_names"]])
    target_groups = np.asarray([str(g) for g in targets_raw["target_groups"]])
    scores = np.asarray(targets_raw["scores"], dtype=np.float64)
    if (scores.ndim != 2 or scores.shape != (len(target_ids), len(target_names))
            or target_groups.shape != target_names.shape or len(set(target_ids)) != len(target_ids)):
        raise ValueError("target artifact violates G1: unique patient IDs and aligned [patient,target] scores are required")
    keep_targets = np.ones(len(target_names), dtype=bool)
    if args.target_group:
        keep_targets = target_groups == args.target_group
    # Exclude the random-control columns from the molecular block itself; they are
    # a null for a different question and would dilute the channel estimate.
    keep_targets &= ~np.char.startswith(target_names, "RANDOM_CONTROL__")
    scores = scores[:, keep_targets]
    target_names = target_names[keep_targets]
    if scores.shape[1] == 0 or not np.isfinite(scores).all() or len(set(target_names)) != len(target_names):
        raise ValueError("selected RNA targets are empty, non-finite, or have duplicate target names")
    target_index = {pid: i for i, pid in enumerate(target_ids)}

    rows: list[dict] = []
    summaries: dict[str, dict] = {}
    channel_gate_failures: list[str] = []
    rna_positive_controls: list[tuple[str, str, str, bool]] = []

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
        purity_values = None
        if purity_by_patient:
            purity_values = np.asarray([purity_by_patient.get(patient, np.nan) for patient in patient_ids])
            purity_coverage = float(np.isfinite(purity_values[mask]).mean()) if mask.any() else 0.0
            rows.append(_row(method=method, task="calibra", metric="purity_coverage", value=purity_coverage,
                             note=f"source={args.purity_source}"))
            if purity_coverage < 0.80:
                rows.append(_unavailable(method, "", "calibra", "insufficient_purity_coverage"))
                continue
            # This is an explicit complete-case confound-adjusted analysis.
            # Missingness is not imputed because doing so would turn a partly
            # unadjusted comparison into an apparently purity-adjusted result.
            mask &= np.isfinite(purity_values)
            rows.append(_row(method=method, task="calibra", metric="n_purity_complete", value=float(mask.sum()),
                             note="complete_case_purity_adjustment"))
            if mask.sum() < 50:
                rows.append(_unavailable(method, "", "calibra", "insufficient_complete_purity_cases"))
                continue
            observed_purity = purity_values[mask & np.isfinite(purity_values)]
            if observed_purity.size < 2 or np.std(observed_purity) < 1e-8:
                rows.append(_unavailable(method, "", "calibra", "constant_purity"))
                continue
            rows.append(_row(method=method, task="calibra", metric="purity_min", value=float(observed_purity.min()),
                             note=f"source={args.purity_source}"))
            rows.append(_row(method=method, task="calibra", metric="purity_max", value=float(observed_purity.max()),
                             note=f"source={args.purity_source}"))
            if args.purity_source == "expression_derived":
                rows.append(_row(method=method, task="calibra", metric="purity_circularity", value=1.0,
                                 note="expression_derived_from_same_RNA_targets; interpret as partial-circularity sensitivity only"))
        if mask.sum() < 50:
            rows.append(_unavailable(method, "", "calibra", f"insufficient_paired_patients_{int(mask.sum())}"))
            continue

        y = scores[aligned[mask]]
        # G1.1/G1.2 apply to the actual complete-case cohort, not merely the
        # source table.  A constant target invalidates CCA and can turn a broken
        # pipeline into a superficially plausible NaN-free result.
        try:
            validate_evaluated_targets(y, target_names)
        except ValueError as exc:
            rows.append(_unavailable(method, "", "calibra", f"nonfinite_or_constant_evaluated_targets_{exc}"))
            continue
        # Pool rare sites. ~600 TSS codes with 145 singletons means a rare-site dummy
        # is effectively a per-patient indicator, and residualising it would delete
        # that patient's signal entirely rather than adjust for a site effect.
        tss, frequent = pooled_tissue_source_site(patient_ids[mask], min_site_count=args.min_site_count)
        frame = pd.DataFrame({"cancer": cancers[mask], "tss": tss})
        base_confounds = ["cancer", "tss"]
        base_design = confound_design(frame, base_confounds)
        confounds = base_confounds.copy()
        if purity_values is not None:
            frame["purity"] = purity_values[mask]
            confounds.append("purity")
        design = confound_design(frame, confounds)
        rows.append(_row(method=method, task="calibra", metric="n_patients", value=float(mask.sum()),
                         note=f"partition={args.partition}"))
        rows.append(_row(method=method, task="calibra", metric="n_confound_columns",
                         value=float(design.shape[1]), note="+".join(confounds) + "_pooled"))
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

            summary, result = _channel_measurement(
                x, y, design, frame["cancer"].to_numpy(), levels=levels, n_draws=args.n_draws,
                n_components=args.n_components, n_permutations=args.n_permutations,
                seed=args.seed, n_jobs=args.n_jobs,
            )
            if not summary["baseline_is_null_like"]:
                channel_gate_failures.append(f"{method}::{state}")
            if not summary["permutation_null_non_degenerate"]:
                channel_gate_failures.append(f"{method}::{state}::degenerate_permutation_null")
            if purity_values is not None:
                # The before/after purity comparison must use the exact same
                # complete-case patients; otherwise a cohort-composition shift
                # can masquerade as a confound effect.
                base_summary, base_result = _channel_measurement(
                    x, y, base_design, frame["cancer"].to_numpy(), levels=levels, n_draws=args.n_draws,
                    n_components=args.n_components, n_permutations=args.n_permutations,
                    seed=args.seed, n_jobs=args.n_jobs,
                )
                if not base_summary["baseline_is_null_like"]:
                    channel_gate_failures.append(f"{method}::{state}::before_purity")
                if not base_summary["permutation_null_non_degenerate"]:
                    channel_gate_failures.append(f"{method}::{state}::before_purity::degenerate_permutation_null")
                summaries[f"{method}::{state}::before_purity"] = base_summary
                if state.startswith("rna_"):
                    rna_positive_controls.append((method, state, "before", bool(
                        np.isfinite(base_summary["heldout_top_cca"])
                        and base_summary["observed"] > base_summary["null_p95"]
                        and base_summary["heldout_top_cca"] > 0.10)))
                for metric in _CHANNEL_METRICS:
                    before = _channel_value(base_summary, metric)
                    after = _channel_value(summary, metric)
                    rows.extend([
                        _row(method=method, representation_state=state, task="calibra_purity_paired",
                             target="before", metric=metric, value=before,
                             note="same_complete_case_patients;cancer+tss;full_channel_measurement"),
                        _row(method=method, representation_state=state, task="calibra_purity_paired",
                             target="after", metric=metric, value=after,
                             note="same_complete_case_patients;cancer+tss+purity;full_channel_measurement"),
                        _row(method=method, representation_state=state, task="calibra_purity_paired",
                             target="after_minus_before", metric=f"{metric}_delta", value=float(after - before),
                             note="paired_complete_case_comparison"),
                    ])
            summaries[f"{method}::{state}"] = summary

            if state.startswith("rna_"):
                rna_positive_controls.append((method, state, "after", bool(
                    np.isfinite(summary["heldout_top_cca"])
                    and summary["observed"] > summary["null_p95"]
                    and summary["heldout_top_cca"] > 0.10)))

            emit = {
                "effective_rank": summary["effective_rank"],
                "unadjusted_top_cca": summary["unadjusted_top_cca"],
                "adjusted_top_cca": summary["observed"],
                "detection_floor": summary["detection_floor"],
                "transmission_floor": summary["transmission_floor"],
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
                "permutation_null_median": summary["null_median"],
                "permutation_null_p95": summary["null_p95"],
                "excess_over_null_median": summary["excess_over_null_median"],
                "permutation_p": summary["permutation_p"],
                "heldout_top_cca": summary["heldout_top_cca"],
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
                  f"floor={summary['detection_floor']} trans={summary['transmission_floor']} "
                  f"atten={summary['attenuation_slope']:.3f} "
                  f"base={summary['baseline_recovered_median']:.4f} {gate}", flush=True)

    frame = pd.DataFrame(rows)
    frame.to_csv(output / "task_rows.csv", index=False)
    (output / "calibra_summary.json").write_text(json.dumps(summaries, indent=2))
    (output / "calibra_protocol.json").write_text(json.dumps({
        "levels": list(levels), "n_draws": args.n_draws, "n_components": args.n_components,
        "seed": args.seed, "partition": args.partition,
        "confounds": ["cancer", "tss", *( ["purity"] if purity_by_patient else [])], "purity": purity_manifest,
        "n_permutations": args.n_permutations, "permutation_strata": "cancer",
        "n_targets": int(scores.shape[1]), "target_group": args.target_group or "all_non_control",
        "recovery_fraction": 0.8, "n_jobs": args.n_jobs,
        "purity_complete_case_paired": bool(purity_by_patient),
        "readout": "targeted_single_direction",
        "readout_note": ("recovery scored on the planted (u,v) axis, not a top-CCA maximum; "
                         "floor units are single-direction correlation and are NOT comparable "
                         "to the multivariate adjusted/held-out top-CCA"),
    }, indent=2))
    positive_by_condition: dict[str, dict[tuple[str, str], bool]] = {}
    for method, state, condition, passed in rna_positive_controls:
        positive_by_condition.setdefault(condition, {})[(method, state)] = bool(passed)
    if purity_by_patient:
        positive_passed = any(before and positive_by_condition.get("after", {}).get(key, False)
                              for key, before in positive_by_condition.get("before", {}).items())
    else:
        positive_passed = any(positive_by_condition.get("after", {}).values())
    gate_status = {"channel_gate_failures": channel_gate_failures,
                   "rna_positive_controls": [{"method": method, "state": state, "condition": condition, "passed": passed}
                                             for method, state, condition, passed in rna_positive_controls],
                   "rna_positive_control_passed": positive_passed,
                   "gates_pass": not channel_gate_failures and (not args.require_rna_positive_control or positive_passed)}
    (output / "calibra_gates.json").write_text(json.dumps(gate_status, indent=2), encoding="utf-8")
    if args.require_channel_gates and channel_gate_failures:
        raise RuntimeError(f"CALIBRA targeted-readout gate failed for {channel_gate_failures[:3]}")
    if args.require_rna_positive_control and not positive_passed:
        raise RuntimeError("G4.1 failed: no rna_* state passed the same-run positive control")
    print(f"[done] {len(frame)} rows -> {output}", flush=True)


if __name__ == "__main__":
    main()
