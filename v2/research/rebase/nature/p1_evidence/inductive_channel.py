"""Is P1's morphology->molecular channel a transductive artefact? The channel and the
labels-only ceiling, re-measured under an INDUCTIVE confound adjustment.

Predeclared in ``NOTEBOOK_ENTRIES/PREDECLARED_inductive_channel_and_ceiling_20260804T2315Z.md``,
committed before this file was written.

**No statistic is defined here.** Every number comes from code that already existed:

* ``nonlinear_adjustment.channel_under_adjustment`` -- P1 4.4's S1 (``top_canonical_correlation``
  at 16 components), S2 (``heldout_top_cca``), ``effective_rank`` and the within-cancer pairing
  permutation null; pinned by ``test_nonlinear_adjustment.py`` to reproduce
  ``calibration.permutation_null`` under the incumbent adjuster
* ``nonlinear_adjustment.labels_only_ceiling`` -- the bound from the other side
* ``nonlinear_adjustment.retention_of_excess`` / ``adjuster_agreement`` / ``cross_fitted_r2``
* ``nonlinear_adjustment._load_block``  -- ``run_calibra``'s own cohort assembly, imported rather
  than re-derived so the cohort is the one P1 4.4 was measured on
* ``inductive_adjustment.ConfoundAdjustmentOperator`` -- the operator, used unchanged
* ``p4_certify.exposure_split`` / ``prepare_state`` / ``_adjustment_audit`` -- the discovery /
  exposure split the P4 inductive entry measured joint site LDA 0.2643 on, imported so that this
  run is on *that* state and not on a re-derived look-alike

Why the split has to be carved out of the certification cohort
--------------------------------------------------------------
D2's own train+val/test split holds out whole cancers, and TCGA assigns a TSS code per site per
disease, so the two sides share **0 cancer types and 0 TSS codes**. An operator fitted on the
representation's own discovery fold would see every test cancer as an unseen level and every test
site as ``OTHER``. See ``PREDECLARED_p4_inductive_adjustment_20260804T2245Z.md`` 2(b).

Why ``adjust_y`` exists
-----------------------
A transductive adjuster is a closure over a design and is column-count agnostic, so one callable
adjusts both blocks. An inductive operator is a *fitted object* whose ridge coefficients have one
row per output column, so the 256-column image block and the 90-column target block need two
operators. ``channel_under_adjustment(..., adjust_y=...)`` takes the second one; its default keeps
every existing caller byte-identical.

The one asymmetry this run must not hide
----------------------------------------
A transductive ``adjust`` **refits** its nuisance model on ``y[order]`` inside every permutation, so
correlation induced by shared residualisation is regenerated inside the null (P1 4.6). An inductive
one applies a **fixed** map keyed to the scored rows' own design, so a permuted patient's un-removed
confound is no longer aligned with that design and is not regenerated. The two arms' nulls are
therefore of different construction. Both nulls are reported for every arm, and the predeclaration
routes the verdict through the ceiling -- which is read against nulls built the same way as the
channel's -- rather than through retention alone.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.inductive_adjustment import ConfoundAdjustmentOperator
from morpheus.v2.calibra.nonlinear_adjustment import (_load_block, adjuster_agreement, cell_design,
                                                      channel_under_adjustment, cross_fitted_r2,
                                                      labels_only_ceiling, make_adjuster,
                                                      retention_of_excess)
from morpheus.v2.research.rebase.nature.p4_certification.p4_certify import (_adjustment_audit,
                                                                           design_for,
                                                                           exposure_split,
                                                                           load_state,
                                                                           prepare_state)

#: The image operator's digest recorded by ``p4_inductive_adjustment_measured_20260804T2300Z.md``
#: for the f = 0.5 discovery-fold operator as fitted for its Tests B and C. If this run's operator
#: does not reproduce it, this run is not on P4's state and the comparison to its 0.2643 is void.
P4_REFERENCE_OPERATOR_DIGEST = "2060a635fa83756a1c3b7aa8506b7b19fcc4431f5d1a303da39b3cb2bf9d62ce"


def fit_inductive_adjuster(full_matrix, fit_rows, score_rows, cancers, patient_ids, *,
                           min_site_count: int = 10, n_splits: int = 5, alpha: float = 1.0,
                           seed: int = 42, cohort_name: str = ""):
    """A ``ConfoundAdjustmentOperator`` fitted on ``fit_rows``, wrapped as an adjuster callable.

    The callable is bound to ``score_rows``' frame and patient ids, because that is what an
    inductive adjustment *is*: a fixed map applied to a fixed set of rows, whose values may then
    be permuted underneath it. ``on_unseen_level="refuse"`` is left on, so a split that put a
    cancer in the scored fold and not in the fitting fold would stop the run rather than silently
    encode those rows as the design's baseline level.
    """
    full_matrix = np.asarray(full_matrix, dtype=np.float64)
    cancers = np.asarray(cancers).astype(str)
    patient_ids = np.asarray(patient_ids).astype(str)
    operator = ConfoundAdjustmentOperator.fit(
        full_matrix[fit_rows], pd.DataFrame({"cancer": cancers[fit_rows]}), ["cancer", "tss"],
        patient_ids=patient_ids[fit_rows], site_column="tss", min_site_count=min_site_count,
        n_splits=n_splits, alpha=alpha, seed=seed, cohort_name=cohort_name)
    frame = pd.DataFrame({"cancer": cancers[score_rows]})
    ids = patient_ids[score_rows]

    def adjust(matrix):
        matrix = np.asarray(matrix, dtype=np.float64)
        if len(matrix) != len(ids):
            raise ValueError(f"inductive adjuster is bound to {len(ids)} scored rows, "
                             f"got {len(matrix)}")
        return operator.adjust(matrix, frame, patient_ids=ids, on_unseen_level="refuse")[0]

    return adjust, operator


def resolve_split_seed(seed: int, split_seed: int) -> int:
    """Which seed carves the discovery/exposure partition.

    ``-1`` means "the published command", i.e. reuse ``--seed`` for the split as
    ``p4_certify.prepare_state`` does. Any other value isolates the partition: it moves which
    patients land in which fold and moves **nothing else**, because the ridge ``KFold``, the
    operator fit, the S2 held-out split and the permutation RNG all stay on ``--seed``.
    """
    return int(seed) if int(split_seed) < 0 else int(split_seed)


def _channel(x, y, adjust, cancers, args, *, adjust_y=None, with_global: bool = True) -> dict:
    """One arm: the within-cancer pairing null, and the global one beside it."""
    started = time.time()
    record = channel_under_adjustment(x, y, adjust, strata=cancers,
                                      n_permutations=args.n_permutations,
                                      n_components=args.n_components, seed=args.seed,
                                      n_jobs=args.n_jobs, adjust_y=adjust_y)
    record["null_strata"] = "cancer"
    if with_global:
        record["global_pairing_null"] = channel_under_adjustment(
            x, y, adjust, strata=None, n_permutations=args.n_permutations,
            n_components=args.n_components, seed=args.seed, n_jobs=args.n_jobs,
            with_heldout=False, adjust_y=adjust_y)
    record["wall_seconds"] = round(time.time() - started, 1)
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--state", default="wsi_biology")
    parser.add_argument("--partition", default="test")
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--n-permutations", type=int, default=2000)
    parser.add_argument("--discovery-fraction", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--split-seed", type=int, default=-1,
                        help="seed for the discovery/exposure PARTITION only; -1 (default) reuses "
                             "--seed, which is the published command. Set it to isolate a change "
                             "of partition from a change of the ridge KFold, the S2 held-out split "
                             "and the permutation RNG, all of which stay on --seed.")
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--skip-full", action="store_true",
                        help="skip the n = 2,766 reproduction gate (for smoke runs only)")
    parser.add_argument("--no-global-null", action="store_true")
    parser.add_argument("--labels-blocks", nargs="+",
                        default=["additive_design", "saturated_cell_design"],
                        choices=["additive_design", "saturated_cell_design",
                                 "frozen_discovery_design"])
    args = parser.parse_args()

    with_global = not args.no_global_null
    path = Path(args.artifact)
    block = _load_block(path, Path(args.targets), args.state, args.partition, args.min_site_count)
    key = f"{path.stem}::{args.state}"
    x_all, y_all = block["x"], block["y"]
    cancers, patient_ids = block["cancers"], block["patient_ids"]

    split_seed = resolve_split_seed(args.seed, args.split_seed)
    discovery = exposure_split(cancers, discovery_fraction=args.discovery_fraction, seed=split_seed)
    fit_rows = np.flatnonzero(discovery)
    score_rows = np.flatnonzero(~discovery)
    overlap = sorted(set(patient_ids[fit_rows].tolist()) & set(patient_ids[score_rows].tolist()))
    if overlap:
        raise AssertionError(f"{len(overlap)} patient(s) in both discovery and exposure folds")

    results: dict = {
        "config": {k: v for k, v in vars(args).items()},
        "state": key,
        "cohort": {
            "n_partition": block["n"], "n_targets": block["n_targets"],
            "n_axes": int(x_all.shape[1]),
            "n_design_columns": int(block["design"].shape[1]),
            "design_rank": int(np.linalg.matrix_rank(block["design"])),
            "n_cells": int(len(block["cell_levels"])), "n_sites_kept": block["n_sites_kept"],
            "n_discovery": int(len(fit_rows)), "n_exposure": int(len(score_rows)),
            "split_seed": int(split_seed),
            "n_cancers_discovery": int(len(set(cancers[fit_rows].tolist()))),
            "n_cancers_exposure": int(len(set(cancers[score_rows].tolist()))),
            "exposure_patient_digest": hashlib.sha256(
                "\x1f".join(sorted(patient_ids[score_rows].tolist())).encode()).hexdigest(),
            "n_patients_in_both_folds": 0,
            "duplicate_patient_ids": int(len(patient_ids) - len(set(patient_ids.tolist()))),
        },
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    def _flush():
        output.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")

    _flush()
    print(f"[cohort] {json.dumps(results['cohort'])}", flush=True)

    # ---- the exposure fold, and the three adjusters that act on it -----------------------
    x_e, y_e = x_all[score_rows], y_all[score_rows]
    cancers_e, ids_e = cancers[score_rows], patient_ids[score_rows]
    design_e = design_for(ids_e, cancers_e, args.min_site_count)

    none_adjust = make_adjuster("none")
    transductive_e = make_adjuster("ridge", design=design_e, n_splits=args.n_splits, seed=args.seed)
    inductive_x, operator_x = fit_inductive_adjuster(
        x_all, fit_rows, score_rows, cancers, patient_ids, min_site_count=args.min_site_count,
        n_splits=args.n_splits, seed=args.seed, cohort_name="d2_test_discovery_fold")
    inductive_y, operator_y = fit_inductive_adjuster(
        y_all, fit_rows, score_rows, cancers, patient_ids, min_site_count=args.min_site_count,
        n_splits=args.n_splits, seed=args.seed, cohort_name="d2_test_discovery_fold_targets")

    # ---- provenance: is this P4's state, bit for bit? ------------------------------------
    # The comparison is only defined when the partition is the one ``prepare_state`` builds --
    # ``prepare_state`` derives its own split from its single ``seed``, so on a run whose
    # ``--split-seed`` differs it would be comparing two different sets of patients. Skipping is
    # recorded as a skip rather than reported as a False.
    provenance = {
        "operator_x_reference_digest": operator_x.provenance["reference_digest"],
        "p4_recorded_digest": P4_REFERENCE_OPERATOR_DIGEST,
        "operator_x_digest_matches_p4": bool(
            operator_x.provenance["reference_digest"] == P4_REFERENCE_OPERATOR_DIGEST),
        "operator_y_reference_digest": operator_y.provenance["reference_digest"],
        "operator_x_n_design_columns": int(operator_x.provenance["n_design_columns"]),
        "operator_y_n_design_columns": int(operator_y.provenance["n_design_columns"]),
        "n_frequent_sites_in_discovery_fold": len(operator_x.site_pooling.frequent),
        "n_exposure_design_columns": int(design_e.shape[1]),
    }
    # Site coverage of the exposure fold by the discovery fold's frozen pooling rule: the
    # 43.8% of `p4_inductive_adjustment_measured_20260804T2300Z.md` §4. Read from
    # SitePooling.apply's own report -- nothing is counted here that the operator does not
    # already count when it adjusts a row.
    _, pooling_report = operator_x.site_pooling.apply(ids_e)
    n_covered = int(len(ids_e) - pooling_report["n_pooled_to_other"])
    provenance["site_coverage"] = {
        "n_exposure_rows": int(len(ids_e)),
        "n_exposure_rows_site_adjustable": n_covered,
        "fraction_exposure_rows_site_adjustable": float(n_covered / max(len(ids_e), 1)),
        "n_frequent_sites_in_discovery_fold": len(operator_x.site_pooling.frequent),
        "n_exposure_sites_unknown_to_discovery_fold":
            int(pooling_report["n_sites_unknown_to_reference"]),
    }
    if split_seed == args.seed:
        p4_block = load_state(path, args.state, args.partition)
        p4_state = prepare_state(p4_block["features"], y_all, p4_block["patient_ids"],
                                 p4_block["cancers"], adjustment="inductive",
                                 discovery_fraction=args.discovery_fraction, seed=args.seed,
                                 min_site_count=args.min_site_count)
        provenance.update({
            "adjusted_features_equal_p4_prepare_state": bool(
                np.array_equal(inductive_x(x_e), p4_state["adjusted_features"])),
            "adjusted_targets_equal_p4_prepare_state": bool(
                np.array_equal(inductive_y(y_e), p4_state["adjusted_targets"])),
            "exposure_rows_equal_p4_prepare_state": bool(
                np.array_equal(ids_e, p4_state["patient_ids"])),
        })
    else:
        provenance["p4_prepare_state_comparison"] = (
            f"skipped: --split-seed {split_seed} != --seed {args.seed}, so prepare_state's own "
            f"partition is a different set of patients and the comparison is undefined")
    results["provenance"] = provenance
    print(f"[provenance] {json.dumps(provenance)}", flush=True)
    _flush()

    # ---- arms ----------------------------------------------------------------------------
    arms: list[tuple[str, dict]] = []

    if not args.skip_full:
        full_adjust = make_adjuster("ridge", design=block["design"], n_splits=args.n_splits,
                                    seed=args.seed)
        arms.append(("transductive_full", {"x": x_all, "y": y_all, "adjust": full_adjust,
                                           "adjust_y": None, "cancers": cancers,
                                           "raw": x_all}))
    arms.extend([
        ("none_exposure", {"x": x_e, "y": y_e, "adjust": none_adjust, "adjust_y": None,
                           "cancers": cancers_e, "raw": x_e}),
        ("transductive_exposure", {"x": x_e, "y": y_e, "adjust": transductive_e, "adjust_y": None,
                                   "cancers": cancers_e, "raw": x_e}),
        ("inductive_exposure", {"x": x_e, "y": y_e, "adjust": inductive_x,
                                "adjust_y": inductive_y, "cancers": cancers_e, "raw": x_e}),
    ])

    for name, spec in arms:
        record = _channel(spec["x"], spec["y"], spec["adjust"], spec["cancers"], args,
                          adjust_y=spec["adjust_y"], with_global=with_global)
        adjusted_x = spec["adjust"](spec["x"])
        record["adjustment_audit"] = _adjustment_audit(spec["raw"], adjusted_x)
        record["r2_of_image_removed"] = cross_fitted_r2(spec["x"], spec["adjust"])
        record["r2_of_targets_removed"] = cross_fitted_r2(
            spec["y"], spec["adjust"] if spec["adjust_y"] is None else spec["adjust_y"])
        record["arm"] = name
        results[f"channel::{name}"] = record
        print(f"[channel {name}] S1={record['observed_top_cca']:.4f} "
              f"null_median={record.get('null_median', float('nan')):.4f} "
              f"excess={record.get('excess_over_null_median', float('nan')):.4f} "
              f"S2={record['heldout_top_cca']:.4f} "
              f"rank={record['effective_rank_x_adjusted']:.2f} p={record.get('permutation_p')} "
              f"({record['wall_seconds']:.0f} s)", flush=True)
        _flush()

    # ---- the labels-only ceiling, under both adjustments ---------------------------------
    # The "representation" is the confound labels themselves, encoded on the WHOLE partition so
    # the discovery and exposure rows land in identical columns and so the widths match the
    # published 108 / 105.
    #
    # Three encodings, because the published 6.0-11.2% figure was measured with the labels
    # block and the adjustment design being the SAME 108 columns, and that coincidence does
    # not survive the split: the exposure fold's own pooling keeps 57 columns and the
    # discovery fold's keeps 55, so ~50 columns of a partition-wide labels block lie outside
    # BOTH arms' adjustment spans and pass through untouched. ``frozen_discovery_design``
    # restores the published coincidence -- it is the operator's OWN frozen design spec
    # evaluated on every row of the partition, so it is defined on the discovery rows (which
    # the inductive operator must be fitted on) and spans what both arms adjust against. The
    # first two encodings are kept because the comparison this run exists to make is
    # inductive vs the matched transductive control, and that comparison must not depend on
    # which of the three encodings is chosen.
    available = {"additive_design": lambda: block["design"],
                 "saturated_cell_design": lambda: cell_design(block["codes"]),
                 "frozen_discovery_design": lambda: operator_x._frame_and_design(
                     pd.DataFrame({"cancer": cancers}), patient_ids, "refuse")[0]}
    labels_blocks = {name: available[name]() for name in args.labels_blocks}
    ceiling_arms: list[tuple[str, str]] = []
    if not args.skip_full:
        ceiling_arms.append(("transductive_full", "full"))
    ceiling_arms.extend([("transductive_exposure", "exposure"),
                         ("inductive_exposure", "exposure")])

    for label, labels_all in labels_blocks.items():
        for arm, scope in ceiling_arms:
            started = time.time()
            if scope == "full":
                adjust, adjust_y = make_adjuster("ridge", design=block["design"],
                                                 n_splits=args.n_splits, seed=args.seed), None
                matrix, targets, strata = labels_all, y_all, cancers
            elif arm == "transductive_exposure":
                adjust, adjust_y = transductive_e, None
                matrix, targets, strata = labels_all[score_rows], y_e, cancers_e
            else:
                adjust, _ = fit_inductive_adjuster(
                    labels_all, fit_rows, score_rows, cancers, patient_ids,
                    min_site_count=args.min_site_count, n_splits=args.n_splits, seed=args.seed,
                    cohort_name=f"d2_test_discovery_fold_labels_{label}")
                adjust_y = inductive_y
                matrix, targets, strata = labels_all[score_rows], y_e, cancers_e
            record = labels_only_ceiling(matrix, targets, n_components=args.n_components,
                                         seed=args.seed, adjust=adjust, strata=strata,
                                         n_permutations=args.n_permutations, n_jobs=args.n_jobs,
                                         adjust_y=adjust_y)
            record["wall_seconds"] = round(time.time() - started, 1)
            record["arm"] = arm
            record["labels_block"] = label
            results[f"ceiling::{label}::{arm}"] = record
            adjusted = record["adjusted"]
            print(f"[ceiling {label} {arm}] raw_S1={record['raw_top_cca']:.4f} "
                  f"adjusted_S1={adjusted['observed_top_cca']:.4f} "
                  f"null_median={adjusted.get('null_median', float('nan')):.4f} "
                  f"excess={adjusted.get('excess_over_null_median', float('nan')):.4f} "
                  f"p={adjusted.get('permutation_p')} ({record['wall_seconds']:.0f} s)", flush=True)
            _flush()

    # ---- the two derived comparisons, from the module's own definitions -------------------
    derived: dict = {}
    control = results.get("channel::transductive_exposure")
    inductive = results.get("channel::inductive_exposure")
    if control and inductive:
        derived["retention_inductive_vs_matched_transductive"] = retention_of_excess(inductive,
                                                                                     control)
        derived["adjuster_agreement_inductive_vs_matched_transductive"] = adjuster_agreement(
            inductive_x(x_e), transductive_e(x_e))
    for label in labels_blocks:
        for arm in ("transductive_exposure", "inductive_exposure"):
            ceiling = results.get(f"ceiling::{label}::{arm}")
            channel = results.get(f"channel::{arm}")
            if not (ceiling and channel):
                continue
            derived[f"ceiling_share_of_channel_excess::{label}::{arm}"] = retention_of_excess(
                ceiling["adjusted"], channel)
            derived[f"ceiling_S1_minus_channel_null_median::{label}::{arm}"] = float(
                ceiling["adjusted"]["observed_top_cca"] - channel["null_median"])
    results["derived"] = derived
    print(f"[derived] {json.dumps(derived, default=float)}", flush=True)
    _flush()
    print(f"[done] -> {output}", flush=True)


if __name__ == "__main__":
    main()
