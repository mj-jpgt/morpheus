"""Measure a DIRECTION-MATCHED detection floor, and the like-for-like floor.

Predeclaration: ``NOTEBOOK_ENTRIES/PREDECLARED_direction_matched_floor_20260804T2230Z.md``.

The open hole. Every shipped CALIBRA ``detection_floor`` is measured by planting
a spike along a direction pair that is **random** on at least the image side --
``spike_targets`` drew ``u`` from its rng and, until now, took no image-direction
argument at all. The channel that floor grades, ``spectral.heldout_top_cca``, uses
a pair **fitted** on both sides. Two independent mismatches hide inside "random vs
fitted":

  M1  direction geometry -- the confound-induced level-0 correlation depends on how
      much of (u, v) lies in the 99-column design span, and a fitted direction need
      not have the same overlap as a random one. Sign unknown a priori.
  M2  oracle vs estimated readout -- the floor's readout KNOWS (u, v); the channel's
      must estimate them. An oracle detects a smaller planted signal, so an oracle
      floor is a LOWER BOUND on the floor that applies to the channel's statistic.
      M2 runs in our favour and direction-matching alone does not touch it.

So this script measures a 2x2 (plus the two references), each on all 13 states of a
CALIBRA run, each with the shipped protocol (13-level grid, n_draws=40, k=16,
recovery_fraction 0.8, 2,530 test patients, 99-column cancer + pooled-TSS design):

  F_rand_real   random pair,  real cohort,              oracle readout   <- the shipped floor
  F_match_real  fitted pair,  real cohort,              oracle readout
  F_rand_perm   random pair,  pairing-destroyed cohort, oracle readout
  F_match_perm  fitted pair,  pairing-destroyed cohort, oracle readout
  F_lfl_rand    random pair,  pairing-destroyed cohort, heldout_top_cca readout
  F_lfl_match   fitted pair,  pairing-destroyed cohort, heldout_top_cca readout  <- PRIMARY

Why the pairing-destroyed cohort for the primary cells. Two reasons, both declared in
advance. (a) ``spike_targets`` orthogonalises the planted component against ``x @ u``
in RAW space; on a fitted pair the raw and residualised correlations differ, so a
residue of the pre-existing channel can survive along the planted axis at r_true = 0
and inflate the level-0 baseline -- an instrument artifact, not a floor. Permuting
``y`` within cancer strata removes the channel while preserving the design geometry
and both marginal covariances. (b) A direction fitted on the same rows the floor is
then measured on is circular; the permuted cohort breaks that.

Draw-to-draw variability, which the UNPAIRED floor rule needs, comes from the split
seed: draw k fits its direction pair on a random half of patients under seed
``seed + k``, and the like-for-like readout uses split seed ``seed + k`` too. A
pinned direction with no other randomness gives a degenerate, noiseless floor --
itself a declared failure mode.

Every statistic is imported. The floor rule itself is ``calibration.floors_from_recovery``,
the SAME rule the shipped instrument uses, so the like-for-like cell is not a
re-implementation.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))

from morpheus.v2.calibra.calibration import (floors_from_recovery, spike_recovery_curve,  # noqa: E402
                                             spike_targets)
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,  # noqa: E402
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import heldout_cca_directions, heldout_top_cca  # noqa: E402

try:
    from joblib import Parallel, delayed
except ImportError:                                   # pragma: no cover
    Parallel = None


def build_state(artifact, state: str, target_index: dict, scores: np.ndarray,
                min_site_count: int):
    """x, y, design, strata for one state -- the run_calibra.py construction, verbatim.

    Identical to ``floor_flag_recompute.build_state``; kept in step with it by the
    reconstruction check (the rebuilt pipeline must reproduce the shipped
    ``heldout_top_cca``), which is asserted in the output of every run.
    """
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


def stratified_order(strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """A within-strata permutation of row indices -- ``permutation_null``'s scheme.

    Cancer-level structure is preserved and only the patient-level cross-modal
    pairing is destroyed; permuting globally would conflate the pairing with the
    cohort effect residualisation already removes.
    """
    order = np.arange(len(strata))
    for level in np.unique(strata):
        idx = np.flatnonzero(strata == level)
        order[idx] = rng.permutation(idx)
    return order


def fitted_pairs(x_residual: np.ndarray, y_residual: np.ndarray, *, n_draws: int,
                 n_components: int, seed: int) -> list:
    """One fitted top-canonical direction pair per draw, each from its own half-split.

    The pair is what ``heldout_top_cca`` uses, fit on train rows only. Refitting per
    draw is what gives the pinned-direction floor its draw-to-draw variability; a
    single pair reused for every draw would produce a noiseless, degenerate floor.
    """
    pairs = []
    for draw in range(n_draws):
        order = np.random.default_rng(seed + draw).permutation(len(x_residual))
        train = order[:len(order) // 2]
        u, v = heldout_cca_directions(x_residual, y_residual, train, n_components=n_components)
        if u.size == 0 or v.size == 0:
            raise RuntimeError(f"degenerate CCA fit at draw {draw}")
        pairs.append((u, v))
    return pairs


def random_pairs(p: int, q: int, *, n_draws: int, seed: int) -> list:
    """Random unit pairs, drawn the way ``spike_targets`` draws them (u first, then v).

    Supplied explicitly rather than left to ``spike_recovery_curve``'s internal draw
    so that the random and fitted cells differ in exactly one thing: the pair.
    """
    pairs = []
    for draw in range(n_draws):
        rng = np.random.default_rng(seed + draw)
        u = rng.normal(size=p); u /= np.linalg.norm(u)
        v = rng.normal(size=q); v /= np.linalg.norm(v)
        pairs.append((u, v))
    return pairs


def oracle_floor(x, y, design, pairs, *, levels, n_draws, n_components, seed, n_jobs):
    """The shipped instrument, with the direction pair supplied instead of drawn."""
    result = spike_recovery_curve(x, y, design, levels=levels, n_draws=n_draws,
                                  n_components=n_components, seed=seed,
                                  direction_pairs=pairs, n_jobs=n_jobs)
    summary = result.summary()
    return {"detection_floor": summary["detection_floor"],
            "transmission_floor": summary["transmission_floor"],
            "null_reference_p90": summary["null_reference_p90"],
            "level0_baseline_abs": summary["confound_induced_baseline"],
            "recovered_median": summary["recovered_median"],
            "unpaired_hit_rate": summary["unpaired_hit_rate"],
            "direction_source": summary["direction_source"],
            "readout": "oracle_pinned_pair"}


def _lfl_column(x, y, design, u, v, levels, n_components, seed, resid_seed):
    """One draw of the LIKE-FOR-LIKE cell: plant on (u, v), read with the CHANNEL's
    own statistic, whose directions are refit on train rows of the spiked data."""
    x_residual = cross_fitted_residuals(x, design, seed=resid_seed)
    column = np.full(len(levels), np.nan)
    for i, level in enumerate(levels):
        spiked = spike_targets(x, y, float(level), rng=np.random.default_rng(0),
                               image_direction=u, molecular_direction=v)
        spiked_residual = cross_fitted_residuals(spiked, design, seed=resid_seed)
        column[i] = heldout_top_cca(x_residual, spiked_residual,
                                    n_components=n_components, seed=seed)
    return column


def likeforlike_floor(x, y, design, pairs, *, levels, n_draws, n_components, seed,
                      recovery_fraction, n_jobs):
    """The nearest well-posed comparison: same statistic on both sides of the grade.

    The spike is planted along the supplied pair and scored with
    ``spectral.heldout_top_cca`` at the same component budget the channel uses, so
    the floor is expressed in the units of the number it grades -- directions fitted
    on train rows, correlation scored on rows those directions never saw. The floor
    rule is ``calibration.floors_from_recovery``, unchanged from the instrument.
    """
    items = [(x, y, design, pairs[draw % len(pairs)][0], pairs[draw % len(pairs)][1],
              levels, n_components, seed + draw, seed) for draw in range(n_draws)]
    if n_jobs != 1 and Parallel is not None:
        columns = Parallel(n_jobs=n_jobs, prefer="processes")(
            delayed(_lfl_column)(*item) for item in items)
    else:
        columns = [_lfl_column(*item) for item in items]
    recovered = np.column_stack(columns)
    floors = floors_from_recovery(levels, recovered, recovery_fraction=recovery_fraction)
    return {"detection_floor": floors["detection_floor"],
            "transmission_floor": floors["transmission_floor"],
            "null_reference_p90": floors["null_reference_p90"],
            "level0_baseline_abs": float(np.nanmedian(np.abs(recovered[0]))),
            "recovered_median": np.nanmedian(recovered, axis=1).tolist(),
            "unpaired_hit_rate": floors["unpaired_hit_rate"],
            "direction_source": "supplied_pairs",
            "readout": f"heldout_top_cca(k={n_components})"}


def pair_design_overlap(pairs, design: np.ndarray, x: np.ndarray, y: np.ndarray) -> dict:
    """How much of each side's planted axis the confound design can explain.

    H1's declared discriminator. If a direction-matched floor comes out LOWER than
    the random-direction one, that must be because the fitted axis lies further out
    of the design span -- and this is the number that says whether it does. Reported
    as the R^2 of the design against the projected score, via the residual the
    pipeline itself produces, so it is the same quantity residualisation removes.
    """
    out = {}
    for name, matrix, index in (("image", x, 0), ("molecular", y, 1)):
        explained = []
        for pair in pairs:
            score = matrix @ pair[index]
            residual = cross_fitted_residuals(score[:, None], design, seed=42).ravel()
            total = float(np.var(score))
            explained.append(1.0 - float(np.var(residual)) / total if total > 1e-12 else np.nan)
        out[f"{name}_design_r2_median"] = float(np.nanmedian(explained))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--task-rows", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--levels", default="0.0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.3,0.4,0.5,0.6")
    parser.add_argument("--n-draws", type=int, default=40)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--recovery-fraction", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--permutation-reps", type=int, default=3,
                        help="independent within-cancer permutations for the two matched "
                             "permuted cells; 1 disables the sensitivity check")
    parser.add_argument("--only-state", default="")
    parser.add_argument("--only-method", default="")
    args = parser.parse_args()

    levels = tuple(float(v) for v in args.levels.split(","))
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
        if args.only_method and path.stem != args.only_method:
            continue
        artifact = np.load(path, allow_pickle=True)
        for state in sorted({str(s) for s in artifact["trained_states"]}):
            if args.only_state and state != args.only_state:
                continue
            key = (path.stem, state)
            if key not in shipped.index:
                print(f"[skip] {key} not in task_rows", flush=True)
                continue
            row = shipped.loc[key]
            x, y, design, strata, n_sites = build_state(artifact, state, target_index, scores,
                                                        args.min_site_count)
            x_res = cross_fitted_residuals(x, design, seed=args.seed)
            y_res = cross_fitted_residuals(y, design, seed=args.seed)

            # H6 / reconstruction: the channel we grade must be the shipped one.
            channel = heldout_top_cca(x_res, y_res, n_components=args.n_components,
                                      seed=args.seed)
            shipped_channel = float(row["heldout_top_cca"])

            pairs_fitted = fitted_pairs(x_res, y_res, n_draws=args.n_draws,
                                        n_components=args.n_components, seed=args.seed)
            pairs_random = random_pairs(x.shape[1], y.shape[1], n_draws=args.n_draws,
                                        seed=args.seed)

            # Pairing-destroyed cohort. ONE within-cancer permutation, shared by every
            # permuted cell, so the six cells differ only in the two things under study
            # (which pair, which readout). ``spike_recovery_curve`` takes a single y, so
            # a per-draw permutation is not expressible in the oracle cells and using
            # one there and many elsewhere would make the cells non-comparable. The
            # obvious objection -- that a single permutation could be lucky -- is
            # answered by ``--permutation-reps``, which re-measures the two matched
            # permuted cells at independent permutations and reports the spread.
            y_perm = y[stratified_order(strata, np.random.default_rng(args.seed))]

            cells = {}
            common = dict(levels=levels, n_draws=args.n_draws,
                          n_components=args.n_components, seed=args.seed, n_jobs=args.n_jobs)
            cells["F_rand_real"] = oracle_floor(x, y, design, pairs_random, **common)
            cells["F_match_real"] = oracle_floor(x, y, design, pairs_fitted, **common)
            cells["F_rand_perm"] = oracle_floor(x, y_perm, design, pairs_random, **common)
            cells["F_match_perm"] = oracle_floor(x, y_perm, design, pairs_fitted, **common)
            lfl = dict(levels=levels, n_draws=args.n_draws, n_components=args.n_components,
                       seed=args.seed, recovery_fraction=args.recovery_fraction,
                       n_jobs=args.n_jobs)
            cells["F_lfl_rand"] = likeforlike_floor(x, y_perm, design, pairs_random, **lfl)
            cells["F_lfl_match"] = likeforlike_floor(x, y_perm, design, pairs_fitted, **lfl)

            # Permutation sensitivity for the two matched permuted cells (the ones the
            # verdict rests on): independent within-cancer permutations, same everything
            # else. A floor that moves with the permutation is a floor with a hidden
            # degree of freedom and must be quoted as a range.
            permutation_reps = {"F_match_perm": [], "F_lfl_match": []}
            for rep in range(1, max(args.permutation_reps, 1)):
                y_rep = y[stratified_order(strata, np.random.default_rng(args.seed + 1000 * rep))]
                permutation_reps["F_match_perm"].append(
                    oracle_floor(x, y_rep, design, pairs_fitted, **common)["detection_floor"])
                permutation_reps["F_lfl_match"].append(
                    likeforlike_floor(x, y_rep, design, pairs_fitted, **lfl)["detection_floor"])

            for name, cell in cells.items():
                floor = cell["detection_floor"]
                cell["channel"] = float(channel)
                cell["clears"] = bool(np.isfinite(floor) and abs(channel) > floor)
                cell["ratio"] = float(abs(channel) / floor) if np.isfinite(floor) and floor > 0 \
                    else float("nan")

            record = {
                "method": path.stem, "state": state, "n_patients": int(len(x)),
                "n_design_columns": int(design.shape[1]), "n_sites_kept": int(n_sites),
                "n_image_columns": int(x.shape[1]), "n_target_columns": int(y.shape[1]),
                "shipped_detection_floor": float(row["detection_floor"]),
                "shipped_transmission_floor": float(row["transmission_floor"]),
                "shipped_heldout_top_cca": shipped_channel,
                "reproduced_heldout_top_cca": float(channel),
                "reconstruction_delta": float(abs(channel - shipped_channel)),
                "overlap_fitted": pair_design_overlap(pairs_fitted, design, x, y),
                "overlap_random": pair_design_overlap(pairs_random, design, x, y),
                "permutation_reps_detection_floor": permutation_reps,
                "cells": cells,
            }
            results.append(record)
            summary_line = " ".join(
                f"{name}={cells[name]['detection_floor']}"
                f"({cells[name]['ratio']:.2f}x{'PASS' if cells[name]['clears'] else 'FAIL'})"
                for name in ("F_rand_real", "F_match_real", "F_rand_perm", "F_match_perm",
                             "F_lfl_rand", "F_lfl_match"))
            print(f"[{path.stem}::{state}] shipped_floor={record['shipped_detection_floor']} "
                  f"channel={channel:.4f} d={record['reconstruction_delta']:.2e} {summary_line}",
                  flush=True)
            Path(args.out).write_text(json.dumps(results, indent=2))

    Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"[written] {args.out} ({len(results)} states)", flush=True)


if __name__ == "__main__":
    main()
