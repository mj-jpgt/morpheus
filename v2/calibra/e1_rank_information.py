"""E1: test whether decorrelation-created rank carries molecular information.

E1 compares *matched* frozen representation artifacts before and after the
feature-decorrelation term.  Rank alone is intentionally not an endpoint:
each arm is measured on the same held-out patients, against the same RNA
targets and confound design, with a paired spike-calibrated detection floor.

The module is deliberately an analysis-only CLI.  It never trains a model and
will refuse a pair whose patient registry, split, cancer labels, or manifest
differs beyond the explicitly allowed decorrelation setting.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .calibration import permutation_null, spike_recovery_curve
from .residualise import confound_design, cross_fitted_residuals
from .spectral import effective_rank


def _row(**kwargs: object) -> dict[str, object]:
    row: dict[str, object] = {"method": "e1_rank_information", "representation_state": "",
                              "arm": "", "task": "e1", "target": "", "metric": "",
                              "value": np.nan, "note": ""}
    row.update(kwargs)
    return row


def _digest(value: object) -> str:
    return sha256(json.dumps(value, sort_keys=True, default=str, separators=(",", ":")).encode()).hexdigest()


def stable_rank(x: np.ndarray) -> float:
    """Stable rank from centred singular values; zero for a constant batch."""
    value = np.asarray(x, dtype=np.float64)
    singular = np.linalg.svd(value - value.mean(axis=0, keepdims=True), compute_uv=False)
    if singular.size == 0 or singular[0] <= 1e-12:
        return 0.0
    return float(np.square(singular).sum() / np.square(singular[0]))


def _set_path(value: dict, dotted: str) -> None:
    current: object = value
    parts = dotted.split(".")
    for part in parts[:-1]:
        if not isinstance(current, dict) or part not in current:
            return
        current = current[part]
    if isinstance(current, dict):
        current.pop(parts[-1], None)


def manifest_without(manifest: dict, allowed_paths: Iterable[str]) -> dict:
    """Copy a manifest and remove the only fields permitted to differ."""
    copied = json.loads(json.dumps(manifest))
    for path in allowed_paths:
        _set_path(copied, path)
    return copied


def parse_manifest(raw: np.lib.npyio.NpzFile) -> dict:
    if "manifest_json" not in raw.files:
        raise ValueError("artifact is missing manifest_json; matched-arm provenance cannot be verified")
    try:
        value = json.loads(str(raw["manifest_json"].item()))
    except (ValueError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("artifact manifest_json is not valid JSON") from exc
    if not isinstance(value, dict):
        raise ValueError("artifact manifest_json must encode an object")
    return value


def validate_matched_artifacts(before_path: str | Path, after_path: str | Path, *,
                               allowed_manifest_paths: Iterable[str] =
                               ("decorrelation_weight", "config.decorrelation_weight")) -> dict:
    """Hard-stop unless the two E1 arms are the same experiment except the term.

    Exact registry equality is stronger than an inner join: silently dropping a
    patient would make a rank difference non-comparable.  Manifest comparison
    is likewise exact after removing only the preregistered intervention field.
    """
    with np.load(before_path, allow_pickle=True) as before, np.load(after_path, allow_pickle=True) as after:
        required = ("patient_ids", "cancers", "split")
        for key in required:
            if key not in before.files or key not in after.files:
                raise ValueError(f"matched artifacts require {key}")
            left = np.asarray(before[key]).astype(str)
            right = np.asarray(after[key]).astype(str)
            if not np.array_equal(left, right):
                raise ValueError(f"matched artifacts differ in {key}; E1 comparison is invalid")
        ids = np.asarray(before["patient_ids"]).astype(str)
        if len(ids) == 0 or len(np.unique(ids)) != len(ids):
            raise ValueError("artifact patient_ids must be non-empty and unique")
        left_manifest, right_manifest = parse_manifest(before), parse_manifest(after)
        allowed = tuple(allowed_manifest_paths)
        if manifest_without(left_manifest, allowed) != manifest_without(right_manifest, allowed):
            raise ValueError("matched artifacts have non-intervention manifest differences")
        return {
            "patient_ids": ids,
            "cancers": np.asarray(before["cancers"]).astype(str),
            "split": np.asarray(before["split"]).astype(str),
            "before_manifest": left_manifest,
            "after_manifest": right_manifest,
            "registry_digest": _digest({key: np.asarray(before[key]).astype(str).tolist() for key in required}),
            "matched_manifest_digest": _digest(manifest_without(left_manifest, allowed)),
            "allowed_manifest_paths": list(allowed),
        }


def _whiten_map(a: np.ndarray, n_components: int) -> tuple[np.ndarray, np.ndarray]:
    centre = a.mean(axis=0, keepdims=True)
    _, singular, vectors = np.linalg.svd(a - centre, full_matrices=False)
    keep = singular > (singular.max() * 1e-10) if singular.size and singular.max() > 0 else np.zeros_like(singular, bool)
    k = min(int(n_components), int(keep.sum()))
    if k == 0:
        return centre, np.zeros((a.shape[1], 0), dtype=np.float64)
    return centre, vectors[:k].T / singular[:k]


def heldout_directional_cca(x: np.ndarray, y: np.ndarray, *, n_components: int = 32,
                            seed: int = 42, train_fraction: float = 0.5) -> np.ndarray:
    """Canonical correlations scored one direction at a time on held-out rows.

    Directions are learned exclusively on one random half.  Unlike an
    in-sample CCA spectrum this cannot manufacture a high-dimensional signal by
    choosing directions after seeing the test rows.
    """
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if len(x) != len(y):
        raise ValueError("x and y must have equal rows")
    order = np.random.default_rng(seed).permutation(len(x))
    cut = int(len(x) * train_fraction)
    train, test = order[:cut], order[cut:]
    if len(train) < 10 or len(test) < 10:
        return np.empty(0, dtype=np.float64)
    cx, wx = _whiten_map(x[train], n_components)
    cy, wy = _whiten_map(y[train], n_components)
    if wx.shape[1] == 0 or wy.shape[1] == 0:
        return np.empty(0, dtype=np.float64)
    ux, uy = (x[train] - cx) @ wx, (y[train] - cy) @ wy
    left, _, right_t = np.linalg.svd(ux.T @ uy, full_matrices=False)
    px = ((x[test] - cx) @ wx) @ left
    py = ((y[test] - cy) @ wy) @ right_t.T
    denom = np.sqrt((np.square(px - px.mean(axis=0)).sum(axis=0)) *
                    (np.square(py - py.mean(axis=0)).sum(axis=0)))
    corr = np.divide(((px - px.mean(axis=0)) * (py - py.mean(axis=0))).sum(axis=0), denom,
                     out=np.full(px.shape[1], np.nan), where=denom > 1e-12)
    return np.abs(corr)


def _pool_tss(patient_ids: np.ndarray, min_site_count: int) -> tuple[np.ndarray, int]:
    sites = np.asarray([pid.split("-")[1] if len(pid.split("-")) > 1 else "NA" for pid in patient_ids])
    unique, counts = np.unique(sites, return_counts=True)
    frequent = {site for site, count in zip(unique, counts) if count >= min_site_count}
    return np.asarray([site if site in frequent else "OTHER" for site in sites]), len(frequent)


def _direction_null(x: np.ndarray, y: np.ndarray, design: np.ndarray, cancers: np.ndarray, *,
                    n_components: int, n_permutations: int, seed: int) -> np.ndarray:
    """Within-cancer pairing null for per-direction CCA, preserving lineage mix."""
    rng = np.random.default_rng(seed)
    x_res = cross_fitted_residuals(x, design, seed=seed)
    values: list[np.ndarray] = []
    for draw in range(n_permutations):
        order = np.arange(len(y))
        for cancer in np.unique(cancers):
            idx = np.flatnonzero(cancers == cancer)
            order[idx] = rng.permutation(idx)
        y_res = cross_fitted_residuals(y[order], design, seed=seed + draw + 1)
        values.append(heldout_directional_cca(x_res, y_res, n_components=n_components, seed=seed + draw + 1))
    width = min((len(value) for value in values), default=0)
    return np.stack([value[:width] for value in values]) if width else np.empty((0, 0))


def _load_targets(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    with np.load(path, allow_pickle=True) as raw:
        for key in ("patient_ids", "scores", "target_names"):
            if key not in raw.files:
                raise ValueError(f"targets is missing {key}")
        ids = np.asarray(raw["patient_ids"]).astype(str)
        score = np.asarray(raw["scores"], dtype=np.float64)
        names = np.asarray(raw["target_names"]).astype(str)
    if score.ndim != 2 or score.shape != (len(ids), len(names)) or not np.isfinite(score).all():
        raise ValueError("targets scores must be finite and aligned [patient, target]")
    keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
    if not keep.any():
        raise ValueError("targets contains no non-control molecular columns")
    return ids, score[:, keep], names[keep], {pid: i for i, pid in enumerate(ids)}


def analyse_arm(x: np.ndarray, y: np.ndarray, design: np.ndarray, cancers: np.ndarray, *,
                n_components: int, levels: tuple[float, ...], n_draws: int,
                n_permutations: int, seed: int) -> dict[str, object]:
    """Compute all E1 quantities for one arm through a shared analysis path."""
    x_res, y_res = cross_fitted_residuals(x, design, seed=seed), cross_fitted_residuals(y, design, seed=seed)
    directional = heldout_directional_cca(x_res, y_res, n_components=n_components, seed=seed)
    # n_components=1 makes the calibration statistic a one-direction CCA score;
    # the same finite floor is then applied to each independently held-out
    # canonical direction. The recovery draws are paired by construction.
    recovery = spike_recovery_curve(x, y, design, levels=levels, n_draws=n_draws,
                                    n_components=1, seed=seed)
    floor = float(recovery.detection_floor)
    above = int(np.sum(directional > floor)) if np.isfinite(floor) else 0
    null_directional = _direction_null(x, y, design, cancers, n_components=n_components,
                                       n_permutations=n_permutations, seed=seed)
    null_counts = (null_directional > floor).sum(axis=1) if np.isfinite(floor) and null_directional.size else np.empty(0)
    pair_null = permutation_null(x, y, design, strata=cancers, n_permutations=n_permutations,
                                 n_components=n_components, seed=seed)
    return {
        "effective_rank": effective_rank(x), "stable_rank": stable_rank(x),
        "directional_cca": directional, "detection_floor": floor,
        "direction_count_above_floor": above,
        "information_density": float(above / effective_rank(x)) if effective_rank(x) > 0 else np.nan,
        "recovery": recovery.summary(), "directional_null": null_directional,
        "direction_count_null_median": float(np.median(null_counts)) if null_counts.size else np.nan,
        "direction_count_null_p95": float(np.percentile(null_counts, 95)) if null_counts.size else np.nan,
        "direction_count_permutation_p": float((np.sum(null_counts >= above) + 1) / (len(null_counts) + 1)) if null_counts.size else np.nan,
        "pairing_null": pair_null,
    }


def run_e1(before_path: str | Path, after_path: str | Path, targets_path: str | Path, output: str | Path, *,
           state: str = "wsi_biology", partition: str = "test", n_components: int = 32,
           levels: tuple[float, ...] = (0.0, 0.02, 0.05, 0.10, 0.20), n_draws: int = 10,
           n_permutations: int = 100, min_site_count: int = 10, seed: int = 42,
           allowed_manifest_paths: Iterable[str] = ("decorrelation_weight", "config.decorrelation_weight")) -> pd.DataFrame:
    """Run E1 and atomically emit rows plus protocol/provenance manifests."""
    matched = validate_matched_artifacts(before_path, after_path, allowed_manifest_paths=allowed_manifest_paths)
    target_ids, scores, names, index = _load_targets(targets_path)
    split = matched["split"]
    mask = split == partition
    if mask.sum() < 50:
        raise ValueError(f"partition {partition!r} has fewer than 50 artifact patients")
    mapped = np.asarray([index.get(pid, -1) for pid in matched["patient_ids"]])
    if np.any(mapped[mask] < 0):
        raise ValueError("all selected artifact patients must be present in frozen targets")
    ids, cancers = matched["patient_ids"][mask], matched["cancers"][mask]
    y = scores[mapped[mask]]
    tss, n_frequent = _pool_tss(ids, min_site_count)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    rows: list[dict[str, object]] = []
    payload: dict[str, object] = {}
    with np.load(before_path, allow_pickle=True) as before, np.load(after_path, allow_pickle=True) as after:
        for arm, raw in (("before", before), ("after", after)):
            if state not in raw.files:
                raise ValueError(f"{arm} artifact lacks E1 state {state!r}")
            x = np.asarray(raw[state], dtype=np.float64)[mask]
            if x.ndim != 2 or not np.isfinite(x).all():
                raise ValueError(f"{arm} {state} is not a finite 2D state")
            result = analyse_arm(x, y, design, cancers, n_components=n_components, levels=levels,
                                 n_draws=n_draws, n_permutations=n_permutations, seed=seed)
            payload[arm] = result
            common = {"arm": arm, "representation_state": state}
            for metric in ("effective_rank", "stable_rank", "detection_floor", "direction_count_above_floor",
                           "information_density", "direction_count_null_median", "direction_count_null_p95",
                           "direction_count_permutation_p"):
                rows.append(_row(**common, task="e1_rank_information", metric=metric,
                                 value=float(result[metric]), note="heldout_test;cancer+tss"))
            for direction, value in enumerate(np.asarray(result["directional_cca"])):
                rows.append(_row(**common, task="e1_directional_cca", target=f"canonical_{direction + 1}",
                                 metric="heldout_cca", value=float(value), note="fit_half_score_half"))
            for metric, value in result["pairing_null"].items():
                if isinstance(value, (int, float, np.floating)):
                    rows.append(_row(**common, task="e1_pairing_null", metric=metric, value=float(value),
                                     note="within_cancer_permutation"))
    before_result, after_result = payload["before"], payload["after"]
    for metric, key in (("delta_effective_rank", "effective_rank"), ("delta_stable_rank", "stable_rank"),
                        ("delta_direction_count_above_floor", "direction_count_above_floor"),
                        ("delta_information_density", "information_density")):
        rows.append(_row(arm="paired", representation_state=state, task="e1_paired", metric=metric,
                         value=float(after_result[key] - before_result[key]), note="after_minus_before;matched_registry"))
    frame = pd.DataFrame(rows)
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "task_rows.csv", index=False)
    protocol = {"experiment": "E1_rank_vs_information", "state": state, "partition": partition,
                "n_components": n_components, "levels": list(levels), "n_draws": n_draws,
                "n_permutations": n_permutations, "seed": seed, "confounds": ["cancer", "tss"],
                "min_site_count": min_site_count, "direction_metric": "heldout_per_canonical_direction_cca",
                "floor_metric": "paired_spike_recovery_one_direction_cca",
                "target_names_digest": _digest(names.tolist()), "n_targets": int(len(names)),
                "n_patients": int(mask.sum()), "n_frequent_tss": int(n_frequent),
                "allowed_manifest_paths": list(allowed_manifest_paths)}
    (out / "e1_protocol.json").write_text(json.dumps(protocol, indent=2, sort_keys=True))
    manifest = {**{key: value for key, value in matched.items() if key not in {"patient_ids", "cancers", "split", "before_manifest", "after_manifest"}},
                "before_path": str(Path(before_path).resolve()), "after_path": str(Path(after_path).resolve()),
                "targets_path": str(Path(targets_path).resolve()), "targets_patient_digest": _digest(target_ids.tolist()),
                "before_manifest": matched["before_manifest"], "after_manifest": matched["after_manifest"]}
    (out / "e1_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    np.savez_compressed(out / "e1_arrays.npz", before_directional_cca=before_result["directional_cca"],
                        after_directional_cca=after_result["directional_cca"],
                        before_direction_null=before_result["directional_null"],
                        after_direction_null=after_result["directional_null"])
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-artifact", required=True)
    parser.add_argument("--after-artifact", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--state", default="wsi_biology")
    parser.add_argument("--partition", choices=("train", "val", "test"), default="test")
    parser.add_argument("--n-components", type=int, default=32)
    parser.add_argument("--levels", default="0.0,0.02,0.05,0.10,0.20")
    parser.add_argument("--n-draws", type=int, default=10)
    parser.add_argument("--n-permutations", type=int, default=100)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--allowed-manifest-path", action="append",
                        default=["decorrelation_weight", "config.decorrelation_weight"])
    args = parser.parse_args()
    rows = run_e1(args.before_artifact, args.after_artifact, args.targets, args.output, state=args.state,
                  partition=args.partition, n_components=args.n_components,
                  levels=tuple(float(value) for value in args.levels.split(",")), n_draws=args.n_draws,
                  n_permutations=args.n_permutations, min_site_count=args.min_site_count, seed=args.seed,
                  allowed_manifest_paths=args.allowed_manifest_path)
    print(f"[done] {len(rows)} E1 rows -> {args.output}", flush=True)


if __name__ == "__main__":
    main()
