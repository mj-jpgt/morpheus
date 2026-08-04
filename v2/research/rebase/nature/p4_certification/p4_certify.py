"""P4 system tests: the abstention rule, the five-point rule end to end, and the competitor gap.

Predeclared in ``NOTEBOOK_ENTRIES/PREDECLARED_p4_certification_system_tests_20260804T1750Z.md``,
committed before this file was written.

P4's claim is *you cannot prompt what you cannot certify*. The interface does not exist, but the
set of axes an interface would be **allowed** to answer from is computable today, from artifacts
already on disk. That set is what this script computes.

**No statistic is defined here.** Every number comes from ``v2/calibra``:

* ``confound_certificate.certify_axes``            -- condition 3 of MULTIMODAL_EXPANSION §1
* ``confound_certificate.within_stratum_permutations`` -- the within-cancer permuter
* ``calibration.spike_recovery_curve``             -- condition 2, the CALIBRA detection floor
* ``spectral.heldout_single_direction_correlation``-- the per-axis, per-target channel readout
* ``spectral.paired_absolute_correlation``         -- the within-cohort reference magnitude
* ``leave_sites_out.site_folds`` / ``evaluate_fold``-- condition 4a, whole sites held out
* ``residualise.confound_design`` / ``cross_fitted_residuals`` / ``pooled_tissue_source_site``
* ``inductive_adjustment.ConfoundAdjustmentOperator`` -- the exposable adjusted state (below)

The three subcommands correspond one-to-one to Tests A, B and C of the predeclaration.

The adjustment mode, added 2026-08-04 22:45
-------------------------------------------
The 20:00 end-to-end run found that condition 3 fails on the state a user would be shown, and named
the reason: ``cross_fitted_residuals`` fits its nuisance model on folds of the very rows it scores,
so **the state that certifies cannot be exposed** to a patient outside the fitting cohort. That entry
could only report a counterfactual (28 of 90) obtained by relaxing condition 3.

``--adjustment`` replaces the counterfactual with a measurement. It selects which rows are scored and
how their adjusted coordinates are produced:

``transductive`` (default)
    Every row of the partition, adjusted in sample by ``cross_fitted_residuals``. Unchanged; this is
    the published 20:00 path and it reproduces those numbers exactly. The adjusted state is **not**
    exposable, so the exposed state stays ``raw``.
``transductive_exposure``
    The exposure fold only, adjusted in sample. The matched-n control: same rows, same n and same
    design as the inductive arm, differing only in that the nuisance model saw the rows it scored.
    Without it, a drop in the inductive arm confounds *inductive vs transductive* with *n halved*.
``inductive``
    A ``ConfoundAdjustmentOperator`` is fitted on the **discovery** fold and applied unchanged to the
    exposure fold. Those coordinates are producible for a patient the operator never saw, so the
    adjusted state becomes the **exposed** state and condition 3 is decided on it.

Why the discovery fold is carved out of the certification cohort rather than taken from the
representation's own train/val split: the D2 split holds out whole *cancers*, and TCGA assigns a TSS
code per site per disease, so ``train+val`` and ``test`` share **0 of 32 cancer types and 0 of 610
TSS codes**. An operator fitted on the representation's discovery fold would see every test cancer as
an unseen level and every test site as ``OTHER``, and its adjustment would degenerate to subtracting
a constant. Conditions 1 and 3 therefore cannot be discharged by the same split; see
``NOTEBOOK_ENTRIES/PREDECLARED_p4_inductive_adjustment_20260804T2245Z.md`` §2.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.calibration import spike_recovery_curve
from morpheus.v2.calibra.confound_certificate import certify_axes, within_stratum_permutations
from morpheus.v2.calibra.inductive_adjustment import ConfoundAdjustmentOperator
from morpheus.v2.calibra.leave_sites_out import evaluate_fold, site_folds
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import (heldout_single_direction_correlation,
                                          paired_absolute_correlation)

PLANTED_NAMES = ("PLANT_site", "PLANT_cancer", "PLANT_noise")


# --------------------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------------------

def load_state(artifact: Path, state: str, partition: str = "test") -> dict:
    raw = np.load(artifact, allow_pickle=True)
    split = np.asarray([str(s) for s in raw["split"]])
    mask = np.ones(len(split), dtype=bool) if partition == "all" else (split == partition)
    return {
        "features": np.asarray(raw[state], dtype=np.float64)[mask],
        "patient_ids": np.asarray([str(p) for p in raw["patient_ids"]])[mask],
        "cancers": np.asarray([str(c) for c in raw["cancers"]])[mask],
        "split": split,
        "mask": mask,
        "n_train": int((split == "train").sum()),
        "n_val": int((split == "val").sum()),
        "n_test": int((split == "test").sum()),
    }


def load_targets(path: Path, patient_ids: np.ndarray, groups: tuple[str, ...]) -> dict:
    raw = np.load(path, allow_pickle=True)
    ids = np.asarray([str(p) for p in raw["patient_ids"]])
    names = np.asarray([str(n) for n in raw["target_names"]])
    target_groups = np.asarray([str(g) for g in raw["target_groups"]])
    order = {pid: i for i, pid in enumerate(ids)}
    rows = np.asarray([order[p] for p in patient_ids], dtype=np.int64)
    keep = np.flatnonzero(np.isin(target_groups, list(groups)))
    return {"scores": np.asarray(raw["scores"], dtype=np.float64)[np.ix_(rows, keep)],
            "names": names[keep], "groups": target_groups[keep]}


def design_for(patient_ids: np.ndarray, cancers: np.ndarray, min_site_count: int = 10) -> np.ndarray:
    site, _ = pooled_tissue_source_site(patient_ids, min_site_count=min_site_count)
    return confound_design(pd.DataFrame({"cancer": cancers, "tss": site}), ["cancer", "tss"])


# --------------------------------------------------------------------------------------
# the discovery / exposure split, and the state that is actually exposed
# --------------------------------------------------------------------------------------

ADJUSTMENTS = ("transductive", "transductive_exposure", "inductive")


def exposure_split(cancers: np.ndarray, *, discovery_fraction: float, seed: int) -> np.ndarray:
    """Boolean mask, ``True`` = discovery fold. Permuted WITHIN cancer type.

    Stratifying by cancer is not cosmetic: the operator's design is one-hot over cancer,
    so a cancer present in the exposure fold and absent from the discovery fold is an
    unseen level that ``DesignSpec.transform`` refuses outright. Every cancer with at
    least two patients therefore contributes to both folds, and the refusal is left
    switched on (``on_unseen_level="refuse"``) so that a split which broke this rule
    would stop the run rather than silently adjust those rows as the baseline level.
    """
    cancers = np.asarray(cancers).astype(str)
    if not 0.0 < discovery_fraction < 1.0:
        raise ValueError("discovery_fraction must lie strictly between 0 and 1")
    rng = np.random.default_rng(seed)
    discovery = np.zeros(len(cancers), dtype=bool)
    for cancer in sorted(set(cancers.tolist())):
        rows = rng.permutation(np.flatnonzero(cancers == cancer))
        k = int(round(len(rows) * discovery_fraction))
        if len(rows) >= 2:
            k = min(max(k, 1), len(rows) - 1)
        discovery[rows[:k]] = True
    return discovery


def _adjustment_audit(raw: np.ndarray, adjusted: np.ndarray) -> dict:
    """How much the adjustment actually moved the state.

    Required by §6.1 of the predeclaration: a certificate PASS on a state the
    adjustment barely touched is a statement about a null adjustment, so the
    per-axis raw-vs-adjusted correlation and the residual variance ratio are
    reported next to every verdict rather than left to be asked for.
    """
    raw = np.asarray(raw, dtype=np.float64)
    adjusted = np.asarray(adjusted, dtype=np.float64)
    a = raw - raw.mean(axis=0)
    b = adjusted - adjusted.mean(axis=0)
    denominator = np.sqrt((a * a).sum(axis=0) * (b * b).sum(axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = np.where(denominator > 1e-30, (a * b).sum(axis=0) / denominator, np.nan)
        variance_ratio = np.where(raw.var(axis=0) > 1e-30,
                                  adjusted.var(axis=0) / raw.var(axis=0), np.nan)
    return {
        "per_axis_corr_raw_adjusted_median": float(np.nanmedian(correlation)),
        "per_axis_corr_raw_adjusted_min": float(np.nanmin(correlation)),
        "per_axis_corr_raw_adjusted_max": float(np.nanmax(correlation)),
        "n_axes_corr_above_0_99": int(np.sum(np.nan_to_num(correlation, nan=0.0) > 0.99)),
        "residual_variance_ratio_median": float(np.nanmedian(variance_ratio)),
        "residual_variance_ratio_min": float(np.nanmin(variance_ratio)),
        "residual_variance_ratio_max": float(np.nanmax(variance_ratio)),
    }


def prepare_state(features: np.ndarray, target_scores, patient_ids: np.ndarray,
                  cancers: np.ndarray, *, adjustment: str, discovery_fraction: float,
                  seed: int, min_site_count: int) -> dict:
    """Select the scored rows and produce their adjusted coordinates.

    Returns the raw block, the adjusted block, the confound design on the scored rows,
    and a provenance report. ``target_scores`` may be ``None`` (Test A' plants axes but
    reads no target).

    The detection floor is computed downstream from ``design`` — the scored rows' own
    design — in **every** mode, so the two exposure arms are handed identical floors by
    construction and no difference between them can be a floor difference.
    """
    if adjustment not in ADJUSTMENTS:
        raise ValueError(f"adjustment must be one of {ADJUSTMENTS}")
    features = np.asarray(features, dtype=np.float64)
    patient_ids = np.asarray(patient_ids).astype(str)
    cancers = np.asarray(cancers).astype(str)
    targets = None if target_scores is None else np.asarray(target_scores, dtype=np.float64)

    if adjustment == "transductive":
        rows = np.arange(len(features))
        design = design_for(patient_ids, cancers, min_site_count)
        adjusted = cross_fitted_residuals(features, design, seed=seed)
        adjusted_targets = (None if targets is None
                            else cross_fitted_residuals(targets, design, seed=seed))
        report = {"adjustment": adjustment, "n_scored": int(len(rows)),
                  "n_design_columns": int(design.shape[1]),
                  "exposable": False,
                  "note": "in-sample cross-fitted residuals on every row of the partition; the "
                          "published 20:00 path. Not producible for a new patient, so the exposed "
                          "state stays raw."}
        return {"rows": rows, "features": features, "targets": targets, "design": design,
                "adjusted_features": adjusted, "adjusted_targets": adjusted_targets,
                "patient_ids": patient_ids, "cancers": cancers,
                "certify_internally": True, "exposable": False,
                "audit": _adjustment_audit(features, adjusted), "report": report}

    discovery = exposure_split(cancers, discovery_fraction=discovery_fraction, seed=seed)
    rows = np.flatnonzero(~discovery)
    fit_rows = np.flatnonzero(discovery)
    overlap = sorted(set(patient_ids[fit_rows].tolist()) & set(patient_ids[rows].tolist()))
    if overlap:                       # §6.3 of the predeclaration: the run is void if this fires
        raise AssertionError(f"{len(overlap)} patient(s) in both discovery and exposure folds")
    design = design_for(patient_ids[rows], cancers[rows], min_site_count)

    if adjustment == "transductive_exposure":
        adjusted = cross_fitted_residuals(features[rows], design, seed=seed)
        adjusted_targets = (None if targets is None
                            else cross_fitted_residuals(targets[rows], design, seed=seed))
        report = {"adjustment": adjustment, "n_discovery": int(len(fit_rows)),
                  "n_scored": int(len(rows)), "n_design_columns": int(design.shape[1]),
                  "discovery_fraction": float(discovery_fraction), "exposable": False,
                  "note": "matched-n control: the exposure fold adjusted in sample. Same rows and "
                          "same design as the inductive arm; the nuisance model saw them."}
        return {"rows": rows, "features": features[rows],
                "targets": None if targets is None else targets[rows], "design": design,
                "adjusted_features": adjusted, "adjusted_targets": adjusted_targets,
                "patient_ids": patient_ids[rows], "cancers": cancers[rows],
                "certify_internally": True, "exposable": False,
                "audit": _adjustment_audit(features[rows], adjusted), "report": report}

    # --- inductive: fit on the discovery fold, apply unchanged to the exposure fold ---
    fit_frame = pd.DataFrame({"cancer": cancers[fit_rows]})
    operator = ConfoundAdjustmentOperator.fit(
        features[fit_rows], fit_frame, ["cancer", "tss"], patient_ids=patient_ids[fit_rows],
        site_column="tss", min_site_count=min_site_count, n_splits=5, alpha=1.0, seed=seed,
        cohort_name="d2_test_discovery_fold")
    apply_frame = pd.DataFrame({"cancer": cancers[rows]})
    adjusted, apply_report = operator.adjust(features[rows], apply_frame,
                                             patient_ids=patient_ids[rows],
                                             on_unseen_level="refuse")
    adjusted_targets, target_report = None, None
    if targets is not None:
        target_operator = ConfoundAdjustmentOperator.fit(
            targets[fit_rows], pd.DataFrame({"cancer": cancers[fit_rows]}), ["cancer", "tss"],
            patient_ids=patient_ids[fit_rows], site_column="tss", min_site_count=min_site_count,
            n_splits=5, alpha=1.0, seed=seed, cohort_name="d2_test_discovery_fold_targets")
        adjusted_targets, target_report = target_operator.adjust(
            targets[rows], pd.DataFrame({"cancer": cancers[rows]}),
            patient_ids=patient_ids[rows], on_unseen_level="refuse")

    pooling = operator.site_pooling
    covered = int(np.sum(np.isin(
        np.asarray([p.split("-")[1] if len(p.split("-")) > 1 else "NA"
                    for p in patient_ids[rows]]), np.asarray(pooling.frequent))))
    report = {
        "adjustment": adjustment, "n_discovery": int(len(fit_rows)), "n_scored": int(len(rows)),
        "discovery_fraction": float(discovery_fraction), "exposable": True,
        "n_operator_design_columns": int(operator.provenance["n_design_columns"]),
        "n_exposure_design_columns": int(design.shape[1]),
        "n_frequent_sites_in_discovery_fold": len(pooling.frequent),
        "n_exposure_rows_with_a_site_adjustment": covered,
        "fraction_exposure_rows_with_a_site_adjustment": float(covered / max(len(rows), 1)),
        "n_exposure_rows_pooled_to_OTHER": int(len(rows) - covered),
        "operator_reference_digest": operator.provenance["reference_digest"],
        "operator_provenance": operator.provenance,
        "apply_report": apply_report,
        "target_apply_report": target_report,
        "note": "nuisance model fitted on the discovery fold ONLY and applied unchanged; these "
                "coordinates are producible for a patient the operator never saw, so this is the "
                "state condition 3 is decided on.",
    }
    return {"rows": rows, "features": features[rows],
            "targets": None if targets is None else targets[rows], "design": design,
            "adjusted_features": adjusted, "adjusted_targets": adjusted_targets,
            "patient_ids": patient_ids[rows], "cancers": cancers[rows],
            "certify_internally": False, "exposable": True,
            "audit": _adjustment_audit(features[rows], adjusted), "report": report}


def certificate_arms(state: dict, args, *, n_boot_axes: int | None = None) -> dict:
    """``certify_axes`` on the raw block and on the adjusted block of a prepared state.

    For the transductive arms the adjusted arm is produced by ``certify_axes``'s own
    ``residualise=True`` branch, which is the identical design and seed — keeping the
    published path byte-for-byte. For the inductive arm the block is *already* adjusted
    by an out-of-sample operator, so it is handed over with ``residualise=False``.
    """
    arms = {}
    for arm in ("raw", "adjusted"):
        internal = state["certify_internally"] and arm == "adjusted"
        block = state["features"] if (arm == "raw" or internal) else state["adjusted_features"]
        arms[arm] = certify_axes(block, state["patient_ids"], state["cancers"],
                                 min_site_count=args.min_site_count,
                                 n_permutations=args.n_permutations, seed=args.seed,
                                 n_boot=args.n_boot,
                                 n_boot_axes=(args.n_boot_axes if n_boot_axes is None
                                              else n_boot_axes),
                                 residualise=internal, n_jobs=args.n_jobs)
    return arms


# --------------------------------------------------------------------------------------
# Test A -- the abstention rule, and whether the certificate discriminates at all
# --------------------------------------------------------------------------------------

def planted_axes(patient_ids: np.ndarray, cancers: np.ndarray, *, seed: int = 42,
                 min_site_count: int = 10, snr: float = 1.0) -> np.ndarray:
    """Three axes of KNOWN character, appended to a state before it is certified.

    These are data, not statistics: a site code, a cancer code and pure noise. The
    certificate must refuse the first and accept the other two; any other pattern
    means it does not discriminate. Constructed from labels already on the artifact.

    ``snr`` is the ratio of signal standard deviation to noise standard deviation.
    The predeclaration fixed ``snr = 1.0``; the ladder run over several values is a
    declared addition, because a single-axis nearest-class-mean rule at 85 classes
    is a weak detector by construction and a binary verdict at one SNR would not say
    *how strong* a site code must be before it is refused. ``snr = inf`` is the
    noiseless site code -- the strongest single-axis leak that can exist.
    """
    rng = np.random.default_rng(seed)
    site, _ = pooled_tissue_source_site(patient_ids, min_site_count=min_site_count)
    _, site_index = np.unique(site, return_inverse=True)
    _, cancer_index = np.unique(np.asarray(cancers).astype(str), return_inverse=True)
    columns = []
    for index in (site_index, cancer_index):
        signal = rng.normal(size=int(index.max()) + 1)[index]
        signal = signal - signal.mean()
        scale = float(signal.std())
        noise = 0.0 if not np.isfinite(snr) else (scale if scale > 1e-12 else 1.0) / max(snr, 1e-12)
        columns.append(signal + rng.normal(size=len(signal)) * noise)
    columns.append(rng.normal(size=len(site_index)))
    return np.column_stack(columns)


def answerable(result: dict) -> dict:
    """Fraction of axes a query layer may answer from, under condition 3 alone.

    Strict = per-axis bar AND per-axis CI AND the state's joint test (A1 and A2 and
    A3 of the predeclaration). Per-axis-only drops A3. They differ, and the
    difference is the point: a state can pass every axis one at a time and still
    carry the confound in a direction no per-axis screen looks along.
    """
    if result.get("status") != "scored":
        return {"status": result.get("status", "unavailable")}
    n_axes = int(result["n_axes"])
    breaching = {int(a["axis"]) for a in result["breaching_axes"]}
    # breaching_axes is truncated to 25 for reporting; n_breaching_axes is the count.
    n_breaching = int(result["n_breaching_axes"])
    joint_ok = bool(result["joint_certified"])
    per_axis_only = (n_axes - n_breaching) / n_axes
    return {
        "status": "scored",
        "n_axes": n_axes,
        "n_breaching": n_breaching,
        "joint_certified": joint_ok,
        "joint_lda_balanced_accuracy": float(result["joint_lda_balanced_accuracy"]),
        "joint_null_p95": float(result["joint_null_p95"]),
        "chance_rate": float(result["chance_rate"]),
        "answerable_fraction_per_axis_only": float(per_axis_only),
        "answerable_fraction_strict": float(per_axis_only if joint_ok else 0.0),
        "named_breaching_axes_top": sorted(breaching),
    }


def cmd_abstention(args) -> None:
    out: dict = {"artifacts": {}, "planted": {}}
    targets_design = None
    for artifact_path in args.artifacts:
        path = Path(artifact_path)
        for state in args.states:
            block = load_state(path, state, args.partition)
            key = f"{path.stem}::{state}"
            for adjusted in (False, True):
                result = certify_axes(block["features"], block["patient_ids"], block["cancers"],
                                      min_site_count=args.min_site_count, n_splits=args.n_splits,
                                      n_permutations=args.n_permutations, seed=args.seed,
                                      n_boot=args.n_boot, n_boot_axes=args.n_boot_axes,
                                      residualise=adjusted, n_jobs=args.n_jobs)
                arm = "adjusted" if adjusted else "raw"
                out["artifacts"].setdefault(key, {})[arm] = answerable(result)
                print(f"[{key}/{arm}] {out['artifacts'][key][arm]}", flush=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[done] -> {args.output}", flush=True)


def cmd_plant(args) -> None:
    """Test A': does the certificate discriminate, or does it pass everything?"""
    path = Path(args.artifact)
    block = load_state(path, args.state, args.partition)
    base = block["features"].shape[1]
    out: dict = {"state": f"{path.stem}::{args.state}", "n_patients": int(len(block["patient_ids"])),
                 "predeclared_snr": 1.0, "ladder": {}}
    for snr in args.snr:
        planted = planted_axes(block["patient_ids"], block["cancers"], seed=args.seed,
                               min_site_count=args.min_site_count, snr=float(snr))
        augmented = np.column_stack([block["features"], planted])
        # The planted columns take the SAME adjustment as the real ones. For the inductive
        # arm that means they are part of the matrix the operator is fitted on over the
        # discovery fold and are adjusted out of sample on the exposure fold; otherwise the
        # ladder would grade a state nobody is exposed to.
        planted_state = prepare_state(augmented, None, block["patient_ids"], block["cancers"],
                                      adjustment=args.adjustment,
                                      discovery_fraction=args.discovery_fraction, seed=args.seed,
                                      min_site_count=args.min_site_count)
        out["state_report"] = {k: v for k, v in planted_state["report"].items()
                               if k not in ("operator_provenance", "apply_report",
                                            "target_apply_report")}
        # every planted axis must get a bootstrap CI
        for arm, result in certificate_arms(
                planted_state, args, n_boot_axes=args.n_boot_axes + planted.shape[1]).items():
            breaching = {int(a["axis"]) for a in result["breaching_axes"]}
            accuracy = np.asarray(result["per_axis_balanced_accuracy"])
            null_p95 = np.asarray(result["per_axis_null_p95"])
            axis_p = np.asarray(result["per_axis_permutation_p"])
            record = {
                "n_axes_total": int(result["n_axes"]),
                "chance_rate": float(result["chance_rate"]),
                "joint_certified": bool(result["joint_certified"]),
                "joint_lda_balanced_accuracy": float(result["joint_lda_balanced_accuracy"]),
                "joint_null_p95": float(result["joint_null_p95"]),
                "n_breaching_axes": int(result["n_breaching_axes"]),
                "planted": {
                    name: {"axis": base + i,
                           "balanced_accuracy": float(accuracy[base + i]),
                           "null_p95": float(null_p95[base + i]),
                           "permutation_p": float(axis_p[base + i]),
                           "exceeds_null_p95": bool(accuracy[base + i] > null_p95[base + i]),
                           "breaches": bool((base + i) in breaching)}
                    for i, name in enumerate(PLANTED_NAMES)},
            }
            out["ladder"].setdefault(str(snr), {})[arm] = record
            print(f"[snr={snr}/{arm}] {json.dumps(record['planted'])}", flush=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[done] -> {args.output}", flush=True)


# --------------------------------------------------------------------------------------
# the per-axis x per-target channel grid, shared by Tests B and C
# --------------------------------------------------------------------------------------

def channel_grid(features: np.ndarray, targets: np.ndarray, *, seed: int, n_jobs: int) -> np.ndarray:
    """[axis, target] out-of-fold single-direction correlation, SIGNED.

    ``heldout_single_direction_correlation`` is the only statistic used. It is the
    one expressed in the same units as the CALIBRA detection floor, which is why a
    per-target readout may be graded against that floor at all.
    """
    from joblib import Parallel, delayed

    def _one(axis: int):
        column = features[:, axis][:, None]
        return [heldout_single_direction_correlation(column, targets[:, t], seed=seed)
                for t in range(targets.shape[1])]

    rows = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(_one)(a) for a in range(features.shape[1]))
    return np.asarray(rows, dtype=np.float64)


# --------------------------------------------------------------------------------------
# Test B -- the five conditions, end to end, on one axis
# --------------------------------------------------------------------------------------

def cmd_endtoend(args) -> None:
    path = Path(args.artifact)
    block = load_state(path, args.state, "test")
    targets = load_targets(Path(args.targets), block["patient_ids"], tuple(args.target_groups))
    state = prepare_state(block["features"], targets["scores"], block["patient_ids"],
                          block["cancers"], adjustment=args.adjustment,
                          discovery_fraction=args.discovery_fraction, seed=args.seed,
                          min_site_count=args.min_site_count)
    design = state["design"]
    adjusted_features, adjusted_targets = state["adjusted_features"], state["adjusted_targets"]
    print(f"[state] {json.dumps({k: v for k, v in state['report'].items() if k not in ('operator_provenance', 'apply_report', 'target_apply_report')})}", flush=True)
    print(f"[audit] {json.dumps(state['audit'])}", flush=True)

    grid = channel_grid(adjusted_features, adjusted_targets, seed=args.seed, n_jobs=args.n_jobs)
    flat = np.abs(grid)
    flat[~np.isfinite(flat)] = -np.inf
    axis, target = np.unravel_index(int(np.argmax(flat)), flat.shape)
    selection = {"axis": int(axis), "target": str(targets["names"][target]),
                 "target_group": str(targets["groups"][target]),
                 "heldout_single_direction_correlation": float(grid[axis, target]),
                 "rule": "argmax |heldout_single_direction_correlation| over 256 axes x "
                         f"{targets['scores'].shape[1]} non-control targets, adjusted block"}
    print(f"[selected] {json.dumps(selection)}", flush=True)

    conditions: dict[str, dict] = {}

    # --- condition 1: operator estimated on a discovery fold only -------------------
    split = block["split"]
    raw = np.load(path, allow_pickle=True)
    all_ids = np.asarray([str(p) for p in raw["patient_ids"]])
    test_ids = set(all_ids[split == "test"].tolist())
    dev_ids = set(all_ids[split != "test"].tolist())
    leak = sorted(test_ids & dev_ids)
    conditions["1_discovery_fold_only"] = {
        "verdict": "PASS" if not leak else "FAIL",
        "n_train": block["n_train"], "n_val": block["n_val"], "n_test": block["n_test"],
        "n_patients_in_both_test_and_development": len(leak),
        "readout_fit_out_of_fold": True,
        "adjustment": args.adjustment,
        "nuisance_operator_discovery_fold": state["report"].get("n_discovery"),
        "nuisance_operator_scored_rows": state["report"]["n_scored"],
        "adjusted_state_is_exposable": bool(state["exposable"]),
        "note": "representation trained on train+val, axis read on test; the readout direction is "
                "fit out of fold inside heldout_single_direction_correlation, by construction. "
                "The NUISANCE operator has its own discovery fold, carved out of the test "
                "partition, because train+val shares 0 cancer types and 0 TSS codes with test and "
                "an operator fitted there is inapplicable by construction.",
    }

    # --- condition 2: clears the CALIBRA detection floor ----------------------------
    curve = spike_recovery_curve(state["features"][:, axis][:, None],
                                 state["targets"][:, target][:, None], design,
                                 n_draws=args.n_draws, seed=args.seed, n_jobs=args.n_jobs)
    # SCHEMA 2 (2026-08-04): the amendment below is now the library's own rule, so the
    # channel statistic is handed to the curve and ``observed_above_floor`` is a real
    # verdict instead of the constant it was. ``magnitude`` is computed first for that
    # reason; the criterion itself is unchanged.
    magnitude = float(abs(grid[axis, target]))
    curve.channel_statistic = magnitude
    curve.channel_statistic_name = "abs(heldout_single_direction_correlation)"
    summary = curve.summary()
    floor = summary["detection_floor"]
    # AMENDED CRITERION, see the amendment appended to the predeclaration: the shipped
    # observed_above_floor flag compares a SIGNED correlation along a random direction
    # pair against the floor, and for a one-column x that sign is a coin flip. The
    # magnitude of the single-direction statistic is what the floor is expressed in.
    clears = bool(np.isfinite(floor) and magnitude > floor)
    conditions["2_clears_detection_floor"] = {
        "verdict": "PASS" if clears else "FAIL",
        "criterion": "amended: abs(heldout_single_direction_correlation) > detection_floor",
        "abs_single_direction_correlation": magnitude,
        "detection_floor": floor,
        "transmission_floor": summary.get("transmission_floor"),
        "observed_matched_direction": summary.get("observed_matched_direction"),
        # At schema 1 this was a constant and disagreed with the amended criterion; at
        # schema 2 the library applies the amended criterion itself, so the two agree
        # by construction and the key is kept only so old and new runs stay diffable.
        "shipped_flag_observed_above_floor": bool(summary["observed_above_floor"]),
        "observed_above_floor_status": summary["observed_above_floor_status"],
        "summary_schema_version": int(summary["summary_schema_version"]),
        "floor_is_nan": bool(not np.isfinite(floor)),
        "levels": summary["levels"], "recovered_median": summary["recovered_median"],
        "delta_median": summary["delta_median"],
        "note": "spike_recovery_curve on the single axis against the single target, cancer+pooled "
                "TSS design. A NaN floor is scored FAIL: an unresolvable floor is not a pass.",
    }

    # --- condition 3: the confound certificate --------------------------------------
    certificate: dict[str, dict] = {}
    for arm, result in certificate_arms(state, args).items():
        breaching = {int(a["axis"]) for a in result["breaching_axes"]}
        accuracy = np.asarray(result["per_axis_balanced_accuracy"])
        null_p95 = np.asarray(result["per_axis_null_p95"])
        certificate[arm] = {
            "axis_balanced_accuracy": float(accuracy[axis]),
            "axis_null_p95": float(null_p95[axis]),
            "axis_exceeds_null_p95": bool(accuracy[axis] > null_p95[axis]),
            "axis_breaches": bool(axis in breaching),
            "joint_certified": bool(result["joint_certified"]),
            "joint_lda_balanced_accuracy": float(result["joint_lda_balanced_accuracy"]),
            "joint_null_p95": float(result["joint_null_p95"]),
            "chance_rate": float(result["chance_rate"]),
            "state_n_breaching_axes": int(result["n_breaching_axes"]),
        }
    exposed_name = "adjusted" if state["exposable"] else "raw"
    exposed = certificate[exposed_name]
    conditions["3_confound_certificate"] = {
        "verdict": "PASS" if (not exposed["axis_breaches"] and exposed["joint_certified"]) else "FAIL",
        "exposed_state": exposed_name,
        "raw": certificate["raw"], "adjusted": certificate["adjusted"],
        "adjustment_audit": state["audit"],
        "note": "the exposed state decides the verdict. With a transductive adjustment the adjusted "
                "coordinates cannot be produced for a new patient, so a user is shown a RAW axis; "
                "with the inductive operator the adjusted coordinates ARE producible out of sample, "
                "so the adjusted state is the exposed one.",
    }

    # --- condition 4a: replication in untouched patients (whole sites held out) ------
    sites_raw, _ = pooled_tissue_source_site(state["patient_ids"], min_site_count=1)
    pooled, _ = pooled_tissue_source_site(state["patient_ids"], min_site_count=args.min_site_count)
    fold = site_folds(sites_raw, state["cancers"], n_folds=args.n_folds)
    x_axis = state["features"][:, axis][:, None]
    y_target = state["targets"][:, target][:, None]
    within = paired_absolute_correlation(adjusted_features[:, axis], adjusted_targets[:, target])
    folds = []
    for f in range(args.n_folds):
        held = fold == f
        result = evaluate_fold(x_axis, y_target, state["cancers"], sites_raw, pooled, held,
                               n_components=1, adjust=True, n_permutations=args.n_permutations,
                               n_boot=args.n_boot_sites, seed=args.seed + f)
        result["fold"] = int(f)
        if result.get("status") == "scored":
            rows = np.flatnonzero(held)
            result["signed_within_heldout_sites"] = heldout_single_direction_correlation(
                adjusted_features[rows, axis][:, None], adjusted_targets[rows, target],
                seed=args.seed)
        folds.append(result)
    scored = [f for f in folds if f.get("status") == "scored"]
    observed = np.asarray([f["observed"] for f in scored], dtype=np.float64)
    median_out_of_site = float(np.nanmedian(observed)) if observed.size else float("nan")
    retained = float(median_out_of_site / within) if np.isfinite(within) and within > 1e-12 else float("nan")
    signs = [np.sign(f.get("signed_within_heldout_sites", np.nan)) for f in scored]
    reference_sign = float(np.sign(grid[axis, target]))
    same_sign = bool(scored) and all(s == reference_sign for s in signs if np.isfinite(s))
    conditions["4a_untouched_patients"] = {
        "verdict": "PASS" if (np.isfinite(retained) and retained >= args.retention_bar and same_sign)
                   else "FAIL",
        "within_cohort_paired_absolute_correlation": float(within),
        "median_out_of_site_observed": median_out_of_site,
        "retained_fraction": retained,
        "predeclared_retention_bar": args.retention_bar,
        "same_sign_in_every_scored_fold": same_sign,
        "n_folds_scored": len(scored),
        "n_folds_clearing_null": int(sum(1 for f in scored if f["clears_null"])),
        "n_folds_ci_clearing_null": int(sum(1 for f in scored if f["ci_clears_null"])),
        "folds": folds,
        "note": "whole sites held out via leave_sites_out.site_folds; the fold-level verdicts "
                "clears_null / ci_clears_null are the helper's own predeclared bars and are "
                "reported alongside the retention bar declared for THIS test.",
    }

    # --- condition 4b: external cohort ----------------------------------------------
    conditions["4b_external_cohort"] = {
        "verdict": "UNEVALUABLE",
        "external_material_on_box": args.external_note,
        "note": "claim_guards.no_external_cohort is undischarged for every morphology result on "
                "the project. Scored UNEVALUABLE, which counts as NOT PASSED for the gate.",
    }

    # --- condition 5: failures exposed alongside successes ---------------------------
    from morpheus.v2.calibra.claim_guards import validate_claim
    verdict = validate_claim({"kind": "legible_axis"})
    conditions["5_failures_exposed"] = {
        "verdict": "PASS" if (not verdict.admissible and verdict.as_rows()) else "FAIL",
        "kind": "procedural",
        "certificate_names_breaching_axes": True,
        "inadmissible_claim_emits_status_rows": len(verdict.as_rows()),
        "blockers_named": [b.code for b in verdict.blockers],
        "note": "procedural, scored by inspection: certify_axes names breaching_axes rather than "
                "dropping them, and ClaimVerdict.as_rows emits a visible NaN status row per blocker.",
    }

    verdicts = {name: c["verdict"] for name, c in conditions.items()}
    gate = {"conditions_passed": sum(1 for v in verdicts.values() if v == "PASS"),
            "conditions_failed": sum(1 for v in verdicts.values() if v == "FAIL"),
            "conditions_unevaluable": sum(1 for v in verdicts.values() if v == "UNEVALUABLE"),
            "certifiable": all(v == "PASS" for v in verdicts.values())}
    out = {"selection": selection, "state": state["report"], "adjustment_audit": state["audit"],
           "verdicts": verdicts, "gate": gate, "conditions": conditions}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"[gate] {json.dumps(verdicts)} -> {json.dumps(gate)}", flush=True)
    print(f"[done] -> {args.output}", flush=True)


# --------------------------------------------------------------------------------------
# Test C -- the competitor gap
# --------------------------------------------------------------------------------------

def cmd_competitor(args) -> None:
    path = Path(args.artifact)
    block = load_state(path, args.state, "test")
    targets = load_targets(Path(args.targets), block["patient_ids"], tuple(args.target_groups))
    state = prepare_state(block["features"], targets["scores"], block["patient_ids"],
                          block["cancers"], adjustment=args.adjustment,
                          discovery_fraction=args.discovery_fraction, seed=args.seed,
                          min_site_count=args.min_site_count)
    design = state["design"]
    adjusted_features, adjusted_targets = state["adjusted_features"], state["adjusted_targets"]
    print(f"[state] {json.dumps({k: v for k, v in state['report'].items() if k not in ('operator_provenance', 'apply_report', 'target_apply_report')})}", flush=True)
    print(f"[audit] {json.dumps(state['audit'])}", flush=True)

    grid = channel_grid(adjusted_features, adjusted_targets, seed=args.seed, n_jobs=args.n_jobs)
    magnitude = np.abs(grid)
    magnitude[~np.isfinite(magnitude)] = -np.inf
    best_axis = np.argmax(magnitude, axis=0)                         # one supporting axis per target

    # The certificate is read on the EXPOSED state: raw when the adjustment is transductive
    # (its adjusted coordinates do not exist for a new patient), the inductively adjusted
    # block when it is not.
    arms = certificate_arms(state, args)
    exposed_name = "adjusted" if state["exposable"] else "raw"
    certificate = arms[exposed_name]
    breaching = {int(a["axis"]) for a in certificate["breaching_axes"]}
    accuracy = np.asarray(certificate["per_axis_balanced_accuracy"])
    null_p95 = np.asarray(certificate["per_axis_null_p95"])
    joint_ok = bool(certificate["joint_certified"])

    from joblib import Parallel, delayed

    permutations = within_stratum_permutations(np.arange(len(state["cancers"])), state["cancers"],
                                               n_permutations=args.n_null, seed=args.seed)

    def _one_target(t: int):
        axis = int(best_axis[t])
        column = adjusted_features[:, axis][:, None]
        null = np.asarray([heldout_single_direction_correlation(column,
                                                                adjusted_targets[order, t],
                                                                seed=args.seed)
                           for order in permutations], dtype=np.float64)
        curve_result = spike_recovery_curve(state["features"][:, axis][:, None],
                                            state["targets"][:, t][:, None], design,
                                            n_draws=args.n_draws, seed=args.seed, n_jobs=1)
        curve_result.channel_statistic = float(abs(grid[axis, t]))
        curve_result.channel_statistic_name = "abs(heldout_single_direction_correlation)"
        curve = curve_result.summary()
        return {"target": str(targets["names"][t]), "target_group": str(targets["groups"][t]),
                "axis": axis,
                "correlation": float(grid[axis, t]),
                "abs_correlation": float(abs(grid[axis, t])),
                "null_p95": float(np.nanpercentile(np.abs(null), 95)),
                "null_median": float(np.nanmedian(np.abs(null))),
                "detection_floor": float(curve["detection_floor"]),
                "observed_matched_direction": float(curve.get("observed_matched_direction", np.nan)),
                "shipped_flag_observed_above_floor": bool(curve["observed_above_floor"]),
                "summary_schema_version": int(curve["summary_schema_version"]),
                # AMENDED: magnitude against the floor, per the appended amendment.
                "clears_detection_floor": bool(np.isfinite(curve["detection_floor"])
                                               and abs(float(grid[axis, t])) > curve["detection_floor"]),
                "axis_breaches_certificate": bool(axis in breaching),
                "axis_balanced_accuracy": float(accuracy[axis]),
                "axis_null_p95": float(null_p95[axis])}

    rows = Parallel(n_jobs=args.n_jobs, prefer="processes")(
        delayed(_one_target)(t) for t in range(adjusted_targets.shape[1]))

    for row in rows:
        row["refused_site"] = bool(row["axis_breaches_certificate"] or not joint_ok)
        row["refused_floor"] = bool(not row["clears_detection_floor"])
        row["refused_null"] = bool(not (row["abs_correlation"] > row["null_p95"]))
        row["answered_by_policy_C"] = bool(not (row["refused_site"] or row["refused_floor"]
                                                or row["refused_null"]))
        # The 20:00 entry's counterfactual, computed here on every arm so the projection
        # and the measurement are the same quantity: condition 3 dropped, floor and null kept.
        row["answered_with_condition_3_relaxed"] = bool(not (row["refused_floor"]
                                                             or row["refused_null"]))

    groups = sorted({r["target_group"] for r in rows})

    def _by_group(key: str) -> dict:
        return {g: int(sum(1 for r in rows if r["target_group"] == g and r[key])) for g in groups}

    n_n = len(rows)
    n_c = int(sum(r["answered_by_policy_C"] for r in rows))
    n_relaxed = int(sum(r["answered_with_condition_3_relaxed"] for r in rows))
    out = {
        "state": f"{path.stem}::{args.state}",
        "adjustment": args.adjustment,
        "adjustment_state": state["report"],
        "adjustment_audit": state["audit"],
        "exposed_state": exposed_name,
        "n_patients": int(len(state["patient_ids"])),
        "joint_certified": joint_ok,
        "joint_lda_balanced_accuracy": float(certificate["joint_lda_balanced_accuracy"]),
        "joint_null_p95": float(certificate["joint_null_p95"]),
        "chance_rate": float(certificate["chance_rate"]),
        "n_queries": n_n,
        "n_answered_policy_N_competitor_style": n_n,
        "n_answered_policy_C_certified": n_c,
        "n_answered_with_condition_3_relaxed": n_relaxed,
        "gap": n_n - n_c,
        "n_queries_by_group": {g: int(sum(1 for r in rows if r["target_group"] == g))
                               for g in groups},
        "answered_policy_C_by_group": _by_group("answered_by_policy_C"),
        "answered_with_condition_3_relaxed_by_group": _by_group("answered_with_condition_3_relaxed"),
        "n_targets_with_unresolvable_floor": int(sum(1 for r in rows
                                                     if not np.isfinite(r["detection_floor"]))),
        "refusals_by_reason": {
            "site_certificate": int(sum(r["refused_site"] for r in rows)),
            "below_detection_floor": int(sum(r["refused_floor"] for r in rows)),
            "inside_permutation_null": int(sum(r["refused_null"] for r in rows)),
            "union": int(sum(not r["answered_by_policy_C"] for r in rows)),
        },
        "supporting_correlation_quantiles": {
            q: float(np.nanpercentile([r["abs_correlation"] for r in rows], v))
            for q, v in (("p05", 5), ("p25", 25), ("p50", 50), ("p75", 75), ("p95", 95))},
        "n_null": args.n_null, "n_draws": args.n_draws,
        "rows": rows,
    }
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"[competitor/{args.adjustment}] N={n_n} C={n_c} relaxed={n_relaxed} gap={out['gap']} "
          f"by_group={json.dumps(out['answered_with_condition_3_relaxed_by_group'])} "
          f"reasons={json.dumps(out['refusals_by_reason'])}", flush=True)
    print(f"[done] -> {args.output}", flush=True)


# --------------------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--output", required=True)
    common.add_argument("--min-site-count", type=int, default=10)
    common.add_argument("--n-splits", type=int, default=5)
    common.add_argument("--n-permutations", type=int, default=1000)
    common.add_argument("--n-boot", type=int, default=200)
    common.add_argument("--n-boot-axes", type=int, default=8)
    common.add_argument("--seed", type=int, default=42)
    common.add_argument("--n-jobs", type=int, default=1)
    common.add_argument("--targets", default="")
    common.add_argument("--adjustment", choices=ADJUSTMENTS, default="transductive",
                        help="how the adjusted state is produced; see the module docstring. "
                             "'transductive' is the published 20:00 path and reproduces it.")
    common.add_argument("--discovery-fraction", type=float, default=0.5,
                        help="share of the partition held back to FIT the nuisance operator; "
                             "ignored by --adjustment transductive")
    common.add_argument("--target-groups", nargs="+",
                        default=["hallmark_in_training", "heldout_pathway", "immune_tme",
                                 "tumour_state"])

    a = sub.add_parser("abstention", parents=[common])
    a.add_argument("--artifacts", nargs="+", required=True)
    a.add_argument("--states", nargs="+", default=["wsi_biology", "rna_biology", "full_biology"])
    a.add_argument("--partition", default="test")
    a.set_defaults(func=cmd_abstention)

    p = sub.add_parser("plant", parents=[common])
    p.add_argument("--artifact", required=True)
    p.add_argument("--state", default="wsi_biology")
    p.add_argument("--partition", default="test")
    p.add_argument("--snr", nargs="+", type=float, default=[0.5, 1.0, 2.0, float("inf")])
    p.set_defaults(func=cmd_plant)

    b = sub.add_parser("endtoend", parents=[common])
    b.add_argument("--artifact", required=True)
    b.add_argument("--state", default="wsi_biology")
    b.add_argument("--n-draws", type=int, default=25)
    b.add_argument("--n-folds", type=int, default=5)
    b.add_argument("--n-boot-sites", type=int, default=1000)
    b.add_argument("--retention-bar", type=float, default=0.5)
    b.add_argument("--external-note", default="external/cptac_gdc_rna is RNA-only STAR counts; "
                                              "HEST-1k spatial exists but no D2 axis has been read "
                                              "on it. No paired external H&E+RNA cohort on the box.")
    b.set_defaults(func=cmd_endtoend)

    c = sub.add_parser("competitor", parents=[common])
    c.add_argument("--artifact", required=True)
    c.add_argument("--state", default="wsi_biology")
    c.add_argument("--n-draws", type=int, default=25)
    c.add_argument("--n-null", type=int, default=200)
    c.set_defaults(func=cmd_competitor)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
