"""Paired D2 H-vs-PBS comparison on frozen held-out patient rows.

This is deliberately a post-export analysis: it reads the exact representation
artifacts that CALIBRA received and recomputes the adjusted molecular channel
under paired patient and cancer-cluster resampling.  It never selects a
checkpoint or fits a representation using test rows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from morpheus.v2.calibra.spectral import top_canonical_correlation
from morpheus.v2.paired_bootstrap import paired_multivariate_patient_and_cancer_bootstrap


def _load(path: str):
    raw = np.load(path, allow_pickle=False)
    required = {"patient_ids", "split", "cancers", "trained_states", "wsi_biology"}
    if not required.issubset(raw.files) or "wsi_biology" not in set(raw["trained_states"].astype(str)):
        raise ValueError(f"D2 comparison needs a declared wsi_biology artifact: {path}")
    return raw


def _targets(path: str, groups: list[str] | None = None) -> tuple[dict[str, int], np.ndarray]:
    raw = np.load(path, allow_pickle=False)
    ids = np.asarray(raw["patient_ids"]).astype(str); values = np.asarray(raw["scores"], dtype=np.float64)
    names = np.asarray(raw["target_names"]).astype(str)
    if groups:
        # `hallmark_in_training` is 50 of the 90 non-control targets and is the
        # Hallmark arm's OWN supervision. Scoring both arms on it hands one arm a
        # structural advantage, so the fair contrast restricts to groups neither
        # arm was trained on. Kept opt-in so the unrestricted number stays quotable.
        #
        # Selection is BY GROUP ONLY here, deliberately. The `RANDOM_CONTROL__`
        # name-prefix drop below is the default "all non-control" behaviour;
        # applying it on top of an explicit group request made
        # `--target-groups random_control` select nothing and raise, so the
        # negative control -- the one readout that proves the instrument is not
        # manufacturing a channel out of noise -- could not be run at all.
        # The prefix and the `random_control` group label are in exact
        # correspondence in the frozen artifact (all 90 prefixed names are the
        # 90 rows of that group and no others), so for every non-control group
        # this is identical to the previous mask.
        available = np.asarray(raw["target_groups"]).astype(str)
        unknown = sorted(set(groups) - set(available.tolist()))
        if unknown:
            raise ValueError(f"unknown target groups {unknown}; available: {sorted(set(available.tolist()))}")
        keep = np.isin(available, list(groups))
    else:
        keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
    values, names = values[:, keep], names[keep]
    if not values.shape[1]:
        raise ValueError(f"target-group selection {groups} left no targets")
    if values.ndim != 2 or values.shape != (len(ids), len(names)) or len(set(ids)) != len(ids) or len(set(names)) != len(names):
        raise ValueError("frozen RNA target table fails D2 ID/target uniqueness requirements")
    return {identifier: index for index, identifier in enumerate(ids)}, values


def load_deflation_block(npz_path: str, axes_file: str) -> tuple[dict[str, int], np.ndarray, list[str]]:
    """Patient-level scores for the axes whose subspace is to be regressed out.

    D2.3's remedy for ``claim_guards.proliferation_deflation`` is "re-run with the
    proliferation programme regressed out". The programme is defined by a list of
    PBS axes, and it is removed by appending their patient scores to the confound
    design -- so ``cross_fitted_residuals`` strips that subspace from the targets
    and from BOTH arms identically. Deleting the axes from the readout block
    instead would answer a different question (can the arms predict
    non-proliferation *targets*) and would silently change the pre-registered
    target set.

    The axis list comes from a file rather than the command line so the exact cut
    is a durable artifact next to the results.
    """
    raw = np.load(npz_path, allow_pickle=False)
    ids = np.asarray(raw["patient_ids"]).astype(str)
    names = np.asarray(raw["target_names"]).astype(str)
    values = np.asarray(raw["scores"], dtype=np.float64)
    if values.shape != (len(ids), len(names)) or len(set(ids)) != len(ids) or len(set(names)) != len(names):
        raise ValueError(f"deflation artifact {npz_path} fails ID/axis uniqueness requirements")
    wanted = [line.strip() for line in Path(axes_file).read_text().splitlines() if line.strip()]
    if not wanted:
        raise ValueError(f"deflation axis file {axes_file} is empty; omit --deflate-npz instead")
    if len(set(wanted)) != len(wanted):
        raise ValueError(f"deflation axis file {axes_file} lists an axis twice")
    position = {name: index for index, name in enumerate(names)}
    missing = sorted(set(wanted) - set(position))
    if missing:
        raise ValueError(f"deflation axes absent from {npz_path}: {missing[:5]}")
    columns = values[:, [position[name] for name in wanted]]
    if not np.isfinite(columns).all() or np.any(np.std(columns, axis=0) < 1e-8):
        raise ValueError("deflation axes are non-finite or constant; residualising them is a no-op")
    return {identifier: index for index, identifier in enumerate(ids)}, columns, wanted


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hallmark-artifacts", nargs="+", required=True)
    parser.add_argument("--pbs-artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True); parser.add_argument("--output", required=True)
    parser.add_argument("--repeats", type=int, default=2000); parser.add_argument("--seed", type=int, default=42)
    # D1 reuses this exact comparison for programme_only vs programme_free. The
    # arithmetic is identical; only the arm names change, and a result labelled
    # "pbs_minus_h" for a D1 run would be a lie in the output file.
    parser.add_argument("--label-a", default="hallmark", help="name of the reference arm in the output")
    parser.add_argument("--label-b", default="pbs", help="name of the contrast arm in the output")
    parser.add_argument("--experiment", default="D2")
    parser.add_argument("--target-groups", nargs="*", default=None,
                        help="restrict the readout to these target_groups (e.g. heldout_pathway "
                             "immune_tme tumour_state). Omit to score every non-control target.")
    parser.add_argument("--deflate-npz", default="",
                        help="axis-score artifact (e.g. pbs_targets_k128_v2.npz) whose named axes are "
                             "regressed out of the targets and BOTH arms before the comparison")
    parser.add_argument("--deflate-axes-file", default="",
                        help="text file, one axis name per line, naming the subspace to remove")
    parser.add_argument("--deflate-label", default="",
                        help="name of this cut in the output, e.g. prol_top100 or placebo_random")
    args = parser.parse_args()
    if bool(args.deflate_npz) != bool(args.deflate_axes_file):
        raise ValueError("--deflate-npz and --deflate-axes-file must be given together")
    # The output keys are built from these labels, so a collision would overwrite
    # an arm's path or point estimate with the other's and still emit valid JSON.
    reserved = {"n_test", "confounds", "pairs", "experiment", "arm_a", "arm_b"}
    if args.label_a == args.label_b or {args.label_a, args.label_b} & reserved:
        raise ValueError(f"--label-a/--label-b must differ and avoid the reserved keys {sorted(reserved)}")
    if len(args.hallmark_artifacts) != len(args.pbs_artifacts):
        raise ValueError(f"{args.experiment} requires one seed-matched artifact pair per seed")
    target_index, target_scores = _targets(args.targets, args.target_groups)
    deflate_index, deflate_scores, deflate_names = ({}, None, [])
    if args.deflate_npz:
        deflate_index, deflate_scores, deflate_names = load_deflation_block(
            args.deflate_npz, args.deflate_axes_file)
    result: dict[str, object] = {"experiment": args.experiment, "arm_a": args.label_a,
                                 "arm_b": args.label_b,
                                 "target_groups": args.target_groups or "all_non_control",
                                 "n_targets": int(target_scores.shape[1]),
                                 "deflation": {"enabled": bool(args.deflate_npz),
                                               "label": args.deflate_label,
                                               "npz": args.deflate_npz,
                                               "axes_file": args.deflate_axes_file,
                                               "n_axes": len(deflate_names),
                                               "axes": deflate_names,
                                               "applied_to": "confound design, so removed from the "
                                                             "targets and BOTH arms identically"},
                                 "pairs": []}
    for pair_index, (h_path, i_path) in enumerate(zip(args.hallmark_artifacts, args.pbs_artifacts)):
        h, i = _load(h_path), _load(i_path)
        for field in ("patient_ids", "split", "cancers"):
            if not np.array_equal(h[field].astype(str), i[field].astype(str)):
                raise ValueError(f"paired frozen artifacts disagree on {field}")
        test = h["split"].astype(str) == "test"; ids = h["patient_ids"].astype(str)[test]
        rows = np.asarray([target_index.get(identifier, -1) for identifier in ids])
        if (rows < 0).any():
            raise ValueError("frozen RNA target artifact does not cover a held-out test patient")
        y, hx, ix = target_scores[rows], h["wsi_biology"].astype(np.float64)[test], i["wsi_biology"].astype(np.float64)[test]
        if not np.isfinite(y).all() or np.any(np.std(y, axis=0) < 1e-8):
            raise ValueError("evaluated RNA targets are non-finite or constant")
        cancers = h["cancers"].astype(str)[test]
        tss, frequent_sites = pooled_tissue_source_site(ids, min_site_count=10)
        def paired_metric(actual: np.ndarray, representation: np.ndarray) -> float:
            # labels are selected by the bootstrap function only implicitly, so
            # this function is reserved for the IID patient bootstrap below.
            return top_canonical_correlation(representation, actual, n_components=16)
        # Bootstrap the exact residualised matrices.  Confound residualisation
        # is fitted once on held-out rows for both arms identically; resampling
        # then quantifies the paired channel difference rather than re-selecting.
        frame = pd.DataFrame({"cancer": cancers, "tss": tss})
        columns = ["cancer", "tss"]
        if deflate_scores is not None:
            deflate_rows = np.asarray([deflate_index.get(identifier, -1) for identifier in ids])
            if (deflate_rows < 0).any():
                raise ValueError("deflation artifact does not cover a held-out test patient")
            for offset, name in enumerate(deflate_names):
                frame[name] = deflate_scores[deflate_rows, offset]
                columns.append(name)
        design = confound_design(frame, columns)
        residual_y = cross_fitted_residuals(y, design, seed=args.seed)
        residual_h = cross_fitted_residuals(hx, design, seed=args.seed)
        residual_i = cross_fitted_residuals(ix, design, seed=args.seed)
        comparison = paired_multivariate_patient_and_cancer_bootstrap(
            paired_metric, residual_y, residual_h, residual_i, cancers, repeats=args.repeats, seed=args.seed + pair_index)
        result["pairs"].append({args.label_a: str(Path(h_path).resolve()), args.label_b: str(Path(i_path).resolve()),
                                "n_test": int(test.sum()), f"point_{args.label_a}": paired_metric(residual_y, residual_h),
                                f"point_{args.label_b}": paired_metric(residual_y, residual_i),
                                f"{args.label_b}_minus_{args.label_a}": comparison,
                                "confounds": {"columns": ["cancer", "tss"], "min_site_count": 10,
                                              "n_distinct_sites_kept": len(frequent_sites),
                                              "n_deflated_axes": len(deflate_names),
                                              "n_design_columns": int(design.shape[1])}})
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(result, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
