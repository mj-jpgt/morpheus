"""TEST 1 -- does a COMPOSED readout over existing axes answer untrained pathways
that the single-axis rule refuses? Separates "no signal" (A) from "no readout" (B).

Predeclared in
``NOTEBOOK_ENTRIES/PREDECLARED_p4_composed_readout_and_causal_name_bridge_20260805T0650Z.md``,
committed (``b83ac7c``) before this file was written.

**No statistic is defined here.** Every number comes from ``v2/calibra`` and from
``p4_certify``, both used unchanged:

* ``p4_certify.prepare_state``                       -- the identical inductive exposure fold
* ``p4_certify.channel_grid``                        -- the incumbent [axis, target] grid
* ``spectral.heldout_single_direction_correlation``  -- the incumbent single-axis readout
* ``spectral.heldout_top_cca``                       -- the COMPOSED readout, capacity ``k``
* ``calibration.spike_targets``                      -- plants the spike (identical for both arms)
* ``calibration.floors_from_recovery``               -- turns a recovery matrix into a floor,
  with the same 80%-of-draws rule ``spike_recovery_curve`` uses; it was lifted out of that
  function verbatim for exactly this call site (a caller scoring the spike with a different
  readout)
* ``calibration.spike_recovery_curve``               -- the SHIPPED floor, reported for continuity
* ``residualise.cross_fitted_residuals``             -- the residualiser, in the floor pipeline
* ``confound_certificate.within_stratum_permutations`` -- the within-cancer permuter

The fairness control (§1.3 of the predeclaration)
-------------------------------------------------
A composed readout has more capacity than a single axis, so grading it against the single
axis's floor and null would rig the comparison. Three matches, all reported:

1. **Matched null.** Every arm's null is that arm's OWN readout applied to within-cancer
   permuted targets, 200 draws. For the composed arm the whitening and the canonical
   directions are refit inside every permutation, so the readout's entire capacity is in
   its own null.
2. **Matched floor.** ``spike_targets`` plants the spike into the TARGET along
   ``image_direction = e_a``, the unit vector on the incumbent's supporting axis, in the
   full 256-column space. Both arms therefore receive **the identical planted signal in the
   identical data** and differ only in the readout that reads it back. Because the level-0
   row is scored by each arm's own readout, a higher-capacity readout gets a higher level-0
   upper tail and hence a higher floor -- automatically, not by choice.
3. **Matched biological null.** The 90 ``random_control`` targets (size- and
   expression-variance-matched random gene sets, one per real signature) go through the
   identical pipeline. 24 of them are the PAIRED controls of the 24 ``heldout_pathway``
   targets.

The incumbent additionally gets a **selection-aware null** -- the full 256-axis grid recomputed
inside every permutation, taking the same argmax -- because its published null holds the
argmax-selected axis fixed and does not charge for the selection. The primary comparison still
uses the incumbent's published, generous grading.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from morpheus.v2.calibra.calibration import (floors_from_recovery, spike_recovery_curve,
                                             spike_targets)
from morpheus.v2.calibra.confound_certificate import within_stratum_permutations
from morpheus.v2.calibra.residualise import cross_fitted_residuals
from morpheus.v2.calibra.spectral import (heldout_single_direction_correlation, heldout_top_cca)
from morpheus.v2.research.rebase.nature.p4_certification.p4_certify import (channel_grid,
                                                                            load_state,
                                                                            load_targets,
                                                                            prepare_state)

#: The capacity ladder, fixed in the predeclaration. 32 is primary (``heldout_top_cca``'s own
#: default); 1 is the same readout family, the same split and the same evaluation rows at
#: capacity one, which isolates capacity from readout family.
K_LADDER = (1, 2, 4, 8, 16, 32, 64)
K_PRIMARY = 32

#: The level grid and draw count of ``spike_recovery_curve``, unchanged, so that the
#: matched-readout floor and the shipped floor are read off the same ladder.
LEVELS = (0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40)

TARGET_GROUPS = ("hallmark_in_training", "heldout_pathway", "immune_tme", "tumour_state",
                 "random_control")

#: Groups the (expensive) selection-aware incumbent null is computed for. The 50
#: ``hallmark_in_training`` targets are excluded on cost; they are not the primary endpoint
#: and the published null is quoted for them.
SELECTION_AWARE_GROUPS = ("heldout_pathway", "immune_tme", "tumour_state", "random_control")


def _composed(x: np.ndarray, y_column: np.ndarray, k: int, seed: int) -> float:
    """The composed readout at capacity ``k``. One call, no statistic defined here."""
    return float(heldout_top_cca(x, np.asarray(y_column, dtype=np.float64)[:, None],
                                 n_components=k, seed=seed, train_fraction=0.5))


def _observed(adjusted_features, adjusted_targets, best_axis, seed, n_jobs):
    """Both arms' observed statistics on every target."""
    from joblib import Parallel, delayed

    def _one(t: int):
        return {k: _composed(adjusted_features, adjusted_targets[:, t], k, seed) for k in K_LADDER}

    composed = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(_one)(t) for t in range(adjusted_targets.shape[1]))
    return composed


def _nulls(adjusted_features, adjusted_targets, best_axis, permutations, seed, n_jobs,
           selection_aware_for):
    """Per-target nulls: composed (per k), incumbent fixed-axis, incumbent selection-aware."""
    from joblib import Parallel, delayed
    n_axes = adjusted_features.shape[1]

    def _one(t: int):
        axis = int(best_axis[t])
        column = adjusted_features[:, axis][:, None]
        fixed, composed = [], {k: [] for k in K_LADDER}
        selection = [] if t in selection_aware_for else None
        for order in permutations:
            permuted = adjusted_targets[order, t]
            fixed.append(heldout_single_direction_correlation(column, permuted, seed=seed))
            for k in K_LADDER:
                composed[k].append(_composed(adjusted_features, permuted, k, seed))
            if selection is not None:
                selection.append(max(
                    abs(heldout_single_direction_correlation(adjusted_features[:, a][:, None],
                                                             permuted, seed=seed))
                    for a in range(n_axes)))
        out = {"single_fixed_axis_p95": float(np.nanpercentile(np.abs(fixed), 95)),
               "single_fixed_axis_median": float(np.nanmedian(np.abs(fixed))),
               "composed_p95": {k: float(np.nanpercentile(np.abs(composed[k]), 95))
                                for k in K_LADDER},
               "composed_median": {k: float(np.nanmedian(np.abs(composed[k]))) for k in K_LADDER}}
        out["single_selection_aware_p95"] = (float(np.nanpercentile(selection, 95))
                                             if selection is not None else float("nan"))
        return out

    return Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(_one)(t) for t in range(adjusted_targets.shape[1]))


def _matched_floors(raw_features, raw_targets, design, best_axis, *, n_draws, seed, n_jobs):
    """The capacity-matched floors: ONE planted signal, read back by each arm's own readout.

    The spike is planted into the target along ``e_a`` -- the incumbent's own supporting axis --
    in the full 256-column space, so the two arms are handed byte-identical spiked data. The
    residualisation, the level grid, the draw count and the 80%-of-draws rule
    (``floors_from_recovery``) are ``spike_recovery_curve``'s, unchanged.
    """
    from joblib import Parallel, delayed
    levels = np.asarray(LEVELS, dtype=np.float64)
    x_residual = cross_fitted_residuals(raw_features, design, seed=seed)

    def _one(t: int):
        axis = int(best_axis[t])
        direction = np.zeros(raw_features.shape[1], dtype=np.float64)
        direction[axis] = 1.0
        y = raw_targets[:, t][:, None]
        rng = np.random.default_rng(seed)
        draw_seeds = [int(rng.integers(1 << 31)) for _ in range(n_draws)]
        single = np.full((len(levels), n_draws), np.nan)
        composed = {k: np.full((len(levels), n_draws), np.nan) for k in K_LADDER}
        for d, draw_seed in enumerate(draw_seeds):
            for i, level in enumerate(levels):
                spiked = spike_targets(raw_features, y, float(level),
                                       rng=np.random.default_rng(draw_seed),
                                       image_direction=direction)
                spiked_residual = cross_fitted_residuals(spiked, design, seed=seed)
                single[i, d] = heldout_single_direction_correlation(
                    x_residual[:, axis][:, None], spiked_residual[:, 0], seed=seed)
                for k in K_LADDER:
                    composed[k][i, d] = _composed(x_residual, spiked_residual[:, 0], k, seed)
        return {"single": floors_from_recovery(levels, np.abs(single)),
                "composed": {k: floors_from_recovery(levels, np.abs(composed[k]))
                             for k in K_LADDER}}

    return Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(_one)(t) for t in range(raw_targets.shape[1]))


def _shipped_floor(raw_features, raw_targets, design, best_axis, grid, *, n_draws, seed, n_jobs):
    """The SHIPPED floor, exactly as ``cmd_competitor`` computes it. Continuity only."""
    from joblib import Parallel, delayed

    def _one(t: int):
        axis = int(best_axis[t])
        result = spike_recovery_curve(raw_features[:, axis][:, None], raw_targets[:, t][:, None],
                                      design, n_draws=n_draws, seed=seed, n_jobs=1)
        result.channel_statistic = float(abs(grid[axis, t]))
        result.channel_statistic_name = "abs(heldout_single_direction_correlation)"
        return float(result.summary()["detection_floor"])

    return Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(_one)(t) for t in range(raw_targets.shape[1]))


def _answered(statistic: float, floor: float, null_p95: float) -> bool:
    """Condition 3 relaxed: clears its OWN floor and exceeds its OWN null p95."""
    return bool(np.isfinite(statistic) and np.isfinite(floor)
                and abs(statistic) > floor and abs(statistic) > null_p95)


def _paired_control_test(rows, seed: int, n_draws: int = 10000) -> dict:
    """Clause 3 of (B): more real held-out pathways answered than their OWN matched controls.

    The pairing is exact -- ``RANDOM_CONTROL__<pathway>__0`` is that pathway's size- and
    variance-matched random gene set. The null flips, independently per pair, which member of
    the pair is called "real", which is the only assignment the biology could have changed.
    """
    pairs = []
    by_name = {r["target"]: r for r in rows}
    for r in rows:
        if r["target_group"] != "heldout_pathway":
            continue
        control = by_name.get(f"RANDOM_CONTROL__{r['target']}__0")
        if control is None:
            continue
        pairs.append((bool(r[f"answered_composed_k{K_PRIMARY}"]),
                      bool(control[f"answered_composed_k{K_PRIMARY}"])))
    if not pairs:
        return {"status": "no_pairs_found"}
    real = np.asarray([a for a, _ in pairs], dtype=float)
    fake = np.asarray([b for _, b in pairs], dtype=float)
    observed = float(real.sum() - fake.sum())
    rng = np.random.default_rng(seed)
    flips = rng.integers(0, 2, size=(n_draws, len(pairs))).astype(bool)
    a = np.where(flips, fake[None, :], real[None, :]).sum(axis=1)
    b = np.where(flips, real[None, :], fake[None, :]).sum(axis=1)
    null = a - b
    return {"status": "scored", "n_pairs": len(pairs),
            "n_real_answered": int(real.sum()), "n_control_answered": int(fake.sum()),
            "observed_difference": observed,
            "null_p95": float(np.percentile(null, 95)),
            "permutation_p": float((int(np.sum(null >= observed)) + 1) / (n_draws + 1)),
            "clears_null_p95": bool(observed > float(np.percentile(null, 95)))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--state", default="wsi_biology")
    parser.add_argument("--targets", required=True)
    parser.add_argument("--adjustment", default="inductive")
    parser.add_argument("--discovery-fraction", type=float, default=0.5)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-null", type=int, default=200)
    parser.add_argument("--n-draws", type=int, default=25)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--skip-selection-aware", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    path = Path(args.artifact)
    block = load_state(path, args.state, "test")
    targets = load_targets(Path(args.targets), block["patient_ids"], TARGET_GROUPS)
    state = prepare_state(block["features"], targets["scores"], block["patient_ids"],
                          block["cancers"], adjustment=args.adjustment,
                          discovery_fraction=args.discovery_fraction, seed=args.seed,
                          min_site_count=args.min_site_count)
    adjusted_features, adjusted_targets = state["adjusted_features"], state["adjusted_targets"]
    names, groups = targets["names"], targets["groups"]
    print(f"[state] n={len(state['patient_ids'])} axes={adjusted_features.shape[1]} "
          f"targets={adjusted_targets.shape[1]} exposable={state['exposable']}", flush=True)

    grid = channel_grid(adjusted_features, adjusted_targets, seed=args.seed, n_jobs=args.n_jobs)
    magnitude = np.abs(grid)
    magnitude[~np.isfinite(magnitude)] = -np.inf
    best_axis = np.argmax(magnitude, axis=0)
    print("[grid] done", flush=True)

    composed_observed = _observed(adjusted_features, adjusted_targets, best_axis,
                                  args.seed, args.n_jobs)
    print("[observed] done", flush=True)

    permutations = within_stratum_permutations(np.arange(len(state["cancers"])), state["cancers"],
                                               n_permutations=args.n_null, seed=args.seed)
    selection_aware_for = set() if args.skip_selection_aware else {
        t for t in range(len(names)) if str(groups[t]) in SELECTION_AWARE_GROUPS}
    nulls = _nulls(adjusted_features, adjusted_targets, best_axis, permutations, args.seed,
                   args.n_jobs, selection_aware_for)
    print("[nulls] done", flush=True)

    floors = _matched_floors(state["features"], state["targets"], state["design"], best_axis,
                             n_draws=args.n_draws, seed=args.seed, n_jobs=args.n_jobs)
    print("[matched floors] done", flush=True)
    shipped = _shipped_floor(state["features"], state["targets"], state["design"], best_axis,
                             grid, n_draws=args.n_draws, seed=args.seed, n_jobs=args.n_jobs)
    print("[shipped floors] done", flush=True)

    rows = []
    for t in range(len(names)):
        axis = int(best_axis[t])
        single = float(grid[axis, t])
        row = {"target": str(names[t]), "target_group": str(groups[t]), "axis": axis,
               "single_correlation": single, "abs_single_correlation": float(abs(single)),
               "single_null_p95_published": nulls[t]["single_fixed_axis_p95"],
               "single_null_p95_selection_aware": nulls[t]["single_selection_aware_p95"],
               "single_floor_shipped": float(shipped[t]),
               "single_floor_matched": float(floors[t]["single"]["detection_floor"])}
        row["answered_single_published"] = _answered(
            single, row["single_floor_shipped"], row["single_null_p95_published"])
        row["answered_single_matched_floor"] = _answered(
            single, row["single_floor_matched"], row["single_null_p95_published"])
        row["answered_single_selection_aware"] = _answered(
            single, row["single_floor_shipped"], row["single_null_p95_selection_aware"])
        for k in K_LADDER:
            statistic = float(composed_observed[t][k])
            floor = float(floors[t]["composed"][k]["detection_floor"])
            p95 = float(nulls[t]["composed_p95"][k])
            row[f"composed_k{k}"] = statistic
            row[f"composed_floor_k{k}"] = floor
            row[f"composed_null_p95_k{k}"] = p95
            row[f"answered_composed_k{k}"] = _answered(statistic, floor, p95)
        rows.append(row)

    group_names = sorted({r["target_group"] for r in rows})

    def _counts(key: str) -> dict:
        return {g: int(sum(1 for r in rows if r["target_group"] == g and r[key]))
                for g in group_names}

    summary = {
        "n_by_group": {g: int(sum(1 for r in rows if r["target_group"] == g))
                       for g in group_names},
        "answered_single_published": _counts("answered_single_published"),
        "answered_single_matched_floor": _counts("answered_single_matched_floor"),
        "answered_single_selection_aware": _counts("answered_single_selection_aware"),
        "answered_composed": {f"k{k}": _counts(f"answered_composed_k{k}") for k in K_LADDER},
        "nan_floor_single_matched": {
            g: int(sum(1 for r in rows if r["target_group"] == g
                       and not np.isfinite(r["single_floor_matched"]))) for g in group_names},
        "nan_floor_composed": {
            f"k{k}": {g: int(sum(1 for r in rows if r["target_group"] == g
                                 and not np.isfinite(r[f"composed_floor_k{k}"])))
                      for g in group_names} for k in K_LADDER},
        "floor_composed_above_single_fraction": {
            f"k{k}": float(np.mean([r[f"composed_floor_k{k}"] > r["single_floor_matched"]
                                    for r in rows
                                    if np.isfinite(r[f"composed_floor_k{k}"])
                                    and np.isfinite(r["single_floor_matched"])]))
            for k in K_LADDER},
        "paired_control_test": _paired_control_test(rows, args.seed),
    }
    out = {"state": f"{path.stem}::{args.state}", "adjustment": args.adjustment,
           "adjustment_state": state["report"], "adjustment_audit": state["audit"],
           "n_patients": int(len(state["patient_ids"])), "n_null": args.n_null,
           "n_draws": args.n_draws, "k_ladder": list(K_LADDER), "k_primary": K_PRIMARY,
           "levels": list(LEVELS), "seed": args.seed,
           "summary": summary, "rows": rows}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"[summary] {json.dumps(summary, default=float)}", flush=True)
    print(f"[done] -> {args.output}", flush=True)


if __name__ == "__main__":
    main()
