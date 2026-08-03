"""T1.5 -- must-FAIL control 3: shuffled gene labels on the PBS dictionary.

Two claims are being separated, and the shuffle is what separates them.

  (i)  SUBSPACE. "There is a molecular subspace this representation can read."
       Permuting which gene each dictionary loading names must not destroy it,
       because the encoder still consumes the same expression matrix with the
       same spectrum -- only the labels moved. Pass = the shuffled top-CCA lies
       INSIDE the bootstrap CI of the unshuffled value.

  (ii) ATTRIBUTION. "Gene g drives axis k." This must collapse. Pass = the
       Spearman between true and shuffled per-axis gene loadings is <= 0.05 with
       a CI covering zero.

READ THIS BEFORE QUOTING (ii). ``build_pbs_targets`` applies the shuffle AFTER
the dictionary is fit, by permuting the rows of ``gene_mean``/``gene_basis``.
So (ii) is true *by construction* -- a permuted loading vector is uncorrelated
with itself -- and is a build-integrity check, not an empirical finding. It is
still run and still reported, because a NON-zero result would mean the shuffle
never took effect, but it may not be presented as evidence that our attribution
is shuffle-sensitive in any deeper sense.

And (i) has a sharper consequence than it first appears. After the permutation
the target block is a spectrum-matched RANDOM projection of the same expression
matrix -- i.e. it is the T1.1 random-dictionary baseline arriving by another
road. If (i) "passes", the honest reading is NOT "the subspace claim survived":
it is "any spectrum-matched projection of this expression matrix is equally
legible to this representation, so the readout is reading expression variance,
not the dictionary". That reading is reported alongside the verdict rather than
left for a reviewer to find.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from .spectral import heldout_top_cca

__all__ = ["rank_correlation", "attribution_collapse", "subspace_persistence"]


def _rank(values: np.ndarray) -> np.ndarray:
    return np.argsort(np.argsort(values, kind="stable"), kind="stable").astype(np.float64)


def rank_correlation(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman without scipy: Pearson on ranks."""
    ra, rb = _rank(np.asarray(a, dtype=np.float64)), _rank(np.asarray(b, dtype=np.float64))
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denominator = float(np.linalg.norm(ra) * np.linalg.norm(rb))
    return float(ra @ rb / denominator) if denominator > 1e-30 else float("nan")


def attribution_collapse(true_basis: np.ndarray, shuffled_basis: np.ndarray, *,
                         n_boot: int = 2000, seed: int = 42) -> dict:
    """Per-axis gene-ranking agreement, by axis index AND by best match.

    The by-index statistic is the one T1.5 names. The best-match statistic is
    strictly harder: it asks whether ANY shuffled axis reproduces the gene
    ranking of a true axis, which is the form a reviewer would reach for after
    noticing that axis k of two different fits need not be the same axis.
    """
    true_basis = np.asarray(true_basis, dtype=np.float64)
    shuffled_basis = np.asarray(shuffled_basis, dtype=np.float64)
    if true_basis.shape != shuffled_basis.shape:
        return {"status": "unavailable_shape_mismatch",
                "true_shape": list(true_basis.shape), "shuffled_shape": list(shuffled_basis.shape)}
    n_axes = true_basis.shape[1]
    by_index = np.asarray([rank_correlation(true_basis[:, k], shuffled_basis[:, k]) for k in range(n_axes)])
    ranks_true = np.column_stack([_rank(true_basis[:, k]) for k in range(n_axes)])
    ranks_shuffled = np.column_stack([_rank(shuffled_basis[:, k]) for k in range(n_axes)])
    centred_t = ranks_true - ranks_true.mean(axis=0, keepdims=True)
    centred_s = ranks_shuffled - ranks_shuffled.mean(axis=0, keepdims=True)
    centred_t = centred_t / np.maximum(np.linalg.norm(centred_t, axis=0, keepdims=True), 1e-30)
    centred_s = centred_s / np.maximum(np.linalg.norm(centred_s, axis=0, keepdims=True), 1e-30)
    cross = centred_t.T @ centred_s                       # [true axis, shuffled axis]
    best_match = np.max(np.abs(cross), axis=1)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, n_axes, size=(int(n_boot), n_axes))
    boot_index = np.median(np.abs(by_index)[draws], axis=1)
    boot_best = np.median(best_match[draws], axis=1)
    return {
        "status": "scored", "n_axes": int(n_axes), "n_genes": int(true_basis.shape[0]),
        "median_abs_spearman_by_axis_index": float(np.median(np.abs(by_index))),
        "median_spearman_by_axis_index": float(np.median(by_index)),
        "max_abs_spearman_by_axis_index": float(np.max(np.abs(by_index))),
        "ci95_median_abs_by_index": [float(np.percentile(boot_index, 2.5)),
                                     float(np.percentile(boot_index, 97.5))],
        "median_best_match_abs_spearman": float(np.median(best_match)),
        "max_best_match_abs_spearman": float(np.max(best_match)),
        "ci95_median_best_match": [float(np.percentile(boot_best, 2.5)),
                                   float(np.percentile(boot_best, 97.5))],
        "collapse_threshold": 0.05,
        "attribution_collapsed": bool(float(np.median(np.abs(by_index))) <= 0.05),
        "by_construction_caveat": "the shuffle permutes gene_basis rows after the fit, so a null here is "
                                  "guaranteed; a NON-null would mean the shuffle failed to take effect",
    }


def subspace_persistence(x: np.ndarray, y_true: np.ndarray, y_shuffled: np.ndarray, design: np.ndarray,
                         *, n_components: int = 32, n_boot: int = 500, seed: int = 42) -> dict:
    """Held-out top-CCA for true vs shuffled targets, with a patient bootstrap CI.

    The CI is on the TRUE value, because the T1.5 bar is "the shuffled value lies
    inside the unshuffled CI". A paired difference CI is reported too: it is the
    statistic that actually decides whether the two differ, and a containment
    test on its own can be passed by a wide CI rather than by agreement.
    """
    x_res = cross_fitted_residuals(np.asarray(x, dtype=np.float64), design, seed=seed)
    true_res = cross_fitted_residuals(np.asarray(y_true, dtype=np.float64), design, seed=seed)
    shuffled_res = cross_fitted_residuals(np.asarray(y_shuffled, dtype=np.float64), design, seed=seed)
    observed_true = heldout_top_cca(x_res, true_res, n_components=n_components, seed=seed)
    observed_shuffled = heldout_top_cca(x_res, shuffled_res, n_components=n_components, seed=seed)
    rng = np.random.default_rng(seed)
    boot_true, boot_shuffled = [], []
    n = len(x_res)
    for _ in range(int(n_boot)):
        rows = rng.integers(0, n, size=n)
        boot_true.append(heldout_top_cca(x_res[rows], true_res[rows], n_components=n_components, seed=seed))
        boot_shuffled.append(heldout_top_cca(x_res[rows], shuffled_res[rows], n_components=n_components, seed=seed))
    boot_true = np.asarray(boot_true, dtype=np.float64)
    boot_shuffled = np.asarray(boot_shuffled, dtype=np.float64)
    difference = boot_true - boot_shuffled
    ci_true = [float(np.nanpercentile(boot_true, 2.5)), float(np.nanpercentile(boot_true, 97.5))]
    return {
        "status": "scored", "n_patients": int(n), "n_components": int(n_components),
        "n_bootstrap": int(n_boot),
        "heldout_top_cca_true": float(observed_true),
        "heldout_top_cca_shuffled": float(observed_shuffled),
        "ci95_true": ci_true,
        "ci95_shuffled": [float(np.nanpercentile(boot_shuffled, 2.5)),
                          float(np.nanpercentile(boot_shuffled, 97.5))],
        "paired_difference_true_minus_shuffled": float(observed_true - observed_shuffled),
        "ci95_paired_difference": [float(np.nanpercentile(difference, 2.5)),
                                   float(np.nanpercentile(difference, 97.5))],
        "difference_ci_excludes_zero": bool(np.nanpercentile(difference, 2.5) > 0
                                            or np.nanpercentile(difference, 97.5) < 0),
        "subspace_persists": bool(ci_true[0] <= float(observed_shuffled) <= ci_true[1]),
        "reading_if_it_persists": "a spectrum-matched random projection of the same expression matrix is as "
                                  "legible as the fitted dictionary; the readout is reading expression "
                                  "variance, not the dictionary",
    }


def _load_targets(path: Path):
    raw = np.load(path, allow_pickle=True)
    ids = np.asarray([str(p) for p in raw["patient_ids"]])
    names = np.asarray([str(t) for t in raw["target_names"]])
    return ids, np.asarray(raw["scores"], dtype=np.float64), names, np.asarray(raw["gene_basis"], dtype=np.float64)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--state", default="wsi_biology")
    parser.add_argument("--true-targets", required=True)
    parser.add_argument("--shuffled-targets", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition", default="test", choices=("test", "val", "all"))
    parser.add_argument("--n-components", type=int, default=32)
    parser.add_argument("--n-boot", type=int, default=500)
    parser.add_argument("--n-boot-axes", type=int, default=2000)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    import pandas as pd

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    raw = np.load(args.artifact, allow_pickle=True)
    patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
    cancers = np.asarray([str(c) for c in raw["cancers"]])
    split = np.asarray([str(s) for s in raw["split"]])
    true_ids, true_scores, _, true_basis = _load_targets(Path(args.true_targets))
    index = {pid: i for i, pid in enumerate(true_ids)}
    aligned = np.asarray([index.get(pid, -1) for pid in patient_ids])
    mask = aligned >= 0
    if args.partition != "all":
        mask &= split == args.partition
    x = np.asarray(raw[args.state], dtype=np.float64)[mask]
    tss, _ = pooled_tissue_source_site(patient_ids[mask], min_site_count=args.min_site_count)
    design = confound_design(pd.DataFrame({"cancer": cancers[mask], "tss": tss}), ["cancer", "tss"])
    y_true = true_scores[aligned[mask]]

    results = []
    for path in args.shuffled_targets:
        shuffled_ids, shuffled_scores, _, shuffled_basis = _load_targets(Path(path))
        shuffled_index = {pid: i for i, pid in enumerate(shuffled_ids)}
        shuffled_aligned = np.asarray([shuffled_index.get(pid, -1) for pid in patient_ids[mask]])
        if (shuffled_aligned < 0).any():
            results.append({"shuffled_targets": str(path), "status": "unavailable_cohort_mismatch"})
            continue
        row = {
            "shuffled_targets": str(Path(path).resolve()), "state": args.state,
            "partition": args.partition, "n_confound_columns": int(design.shape[1]),
            "subspace": subspace_persistence(x, y_true, shuffled_scores[shuffled_aligned], design,
                                             n_components=args.n_components, n_boot=args.n_boot,
                                             seed=args.seed),
            "attribution": attribution_collapse(true_basis, shuffled_basis,
                                                n_boot=args.n_boot_axes, seed=args.seed),
        }
        results.append(row)
        print(json.dumps(row, indent=2), flush=True)

    (output / "gene_label_shuffle_control.json").write_text(
        json.dumps({"artifact": str(Path(args.artifact).resolve()),
                    "true_targets": str(Path(args.true_targets).resolve()),
                    "results": results}, indent=2), encoding="utf-8")
    print(f"[done] -> {output}", flush=True)


if __name__ == "__main__":
    main()
