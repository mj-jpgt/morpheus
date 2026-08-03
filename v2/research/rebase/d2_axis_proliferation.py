"""D2.3 — is per-axis legibility just proliferation?

Answers `claim_guards.proliferation_deflation`, the blocker that makes E0's
`transfer` claim inadmissible. The falsifier stated in the ledger is **every
legible axis coming back proliferation-loaded**.

Two things this module is careful about, both of which could turn a null into a
confident wrong answer:

1. **The null for per-axis legibility is not zero.** A cross-validated ridge
   direction through a 256-column representation finds structure by capacity
   alone. The null here is built with the *same* within-cancer row permutation
   that `calibra.calibration.permutation_null` uses (rows of the target block
   permuted inside each cancer stratum, `x` residuals held fixed), so that
   cancer-level structure survives and only patient-level pairing is destroyed.
   A row-shuffle null is a different number and is not substituted for it.

2. **The pre-built `proliferation_loading` is a weak instrument.**
   `build_pbs_targets` defines it as the |loading|-weighted mean over *all*
   basis genes of a binary proliferation flag. The SVD basis is dense, so the
   weights are near-uniform and every axis is compressed towards the background
   rate. This module therefore carries a second, concentrated statistic — the
   proliferation fraction among each axis's top-k |loading| genes — as
   co-primary, and reports the two side by side rather than picking whichever
   agrees with the conclusion.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.run_calibra import score_target_block_per_column

PROLIFERATION_STATISTICS = ("prol_wmean", "prol_top100")


def axis_gene_statistics(gene_basis: np.ndarray, genes: np.ndarray,
                         annotations: pd.DataFrame, *, top_k: int = 100) -> pd.DataFrame:
    """Per-axis proliferation / essentiality, both diluted and concentrated.

    ``*_wmean`` reproduces ``build_pbs_targets``'s |loading|-weighted mean over
    every basis gene. ``*_top{k}`` restricts to the axis's k largest-|loading|
    genes, which is what "this axis is a proliferation axis" actually means.
    """
    basis = np.asarray(gene_basis, dtype=np.float64)
    symbols = np.asarray([str(g).upper() for g in genes])
    table = annotations.copy()
    gene_column = next(c for c in ("gene", "gene_symbol", "symbol") if c in table.columns)
    table[gene_column] = table[gene_column].astype(str).str.upper()
    aligned = table.drop_duplicates(gene_column).set_index(gene_column).reindex(symbols)
    records = []
    for axis in range(basis.shape[1]):
        weight = np.abs(basis[:, axis])
        order = np.argsort(-weight)[:top_k]
        record = {"axis": f"PBS_{axis:03d}", "axis_index": axis,
                  "top_k_loading_mass": float(weight[order].sum() / max(weight.sum(), 1e-12))}
        for name, column in (("prol", "proliferation_loading"), ("ess", "essentiality_loading")):
            values = aligned[column].to_numpy(dtype=np.float64)
            finite = np.isfinite(values)
            record[f"{name}_wmean"] = (float(np.average(values[finite], weights=weight[finite]))
                                       if finite.any() else np.nan)
            record[f"{name}_top{top_k}"] = float(np.nanmean(values[order]))
        records.append(record)
    return pd.DataFrame(records)


def per_axis_legibility(x: np.ndarray, y: np.ndarray, design: np.ndarray, strata: np.ndarray, *,
                        n_permutations: int = 200, seed: int = 42, n_jobs: int = 1) -> dict:
    """Held-out per-axis correlation, with its own within-cancer permutation null.

    The observed statistic and every permutation replicate go through the
    identical code path (`score_target_block_per_column`), so the null cannot be
    made easier to beat than the observation by an accident of implementation.
    """
    from joblib import Parallel, delayed

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    strata = np.asarray(strata)
    x_residual = cross_fitted_residuals(x, design, seed=seed)
    observed = score_target_block_per_column(x_residual, cross_fitted_residuals(y, design, seed=seed),
                                             seed=seed)
    rng = np.random.default_rng(seed)
    orders = []
    for _ in range(n_permutations):
        order = np.arange(len(y))
        for level in np.unique(strata):
            idx = np.flatnonzero(strata == level)
            order[idx] = rng.permutation(idx)
        orders.append(order)

    def _one(order):
        return score_target_block_per_column(x_residual, cross_fitted_residuals(y[order], design, seed=seed),
                                             seed=seed)

    if n_jobs == 1:
        null = np.asarray([_one(o) for o in orders])
    else:
        null = np.asarray(Parallel(n_jobs=n_jobs)(delayed(_one)(o) for o in orders))
    return {"observed": observed,
            "null_median": np.median(null, axis=0),
            "null_p95": np.percentile(null, 95, axis=0),
            "null_matrix": null,
            "permutation_p": (np.sum(null >= observed[None, :], axis=0) + 1) / (n_permutations + 1),
            "n_permutations": int(n_permutations)}


def bimodality_report(values: np.ndarray) -> dict:
    """Is the legibility distribution one mode or two?

    Reported because "the mean axis is not proliferative" is compatible with a
    small proliferative cluster carrying everything. Two independent readings:
    a 2- vs 1-component Gaussian mixture by BIC, and the largest gap in the
    sorted values relative to the spread.
    """
    from sklearn.mixture import GaussianMixture

    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if len(finite) < 8:
        return {"status": "too_few_values"}
    column = finite[:, None]
    bic = {}
    for k in (1, 2, 3):
        model = GaussianMixture(n_components=k, random_state=42, n_init=5).fit(column)
        bic[k] = float(model.bic(column))
    ordered = np.sort(finite)
    gaps = np.diff(ordered)
    spread = float(ordered[-1] - ordered[0]) or 1.0
    return {"status": "scored", "bic": bic, "best_n_components": int(min(bic, key=bic.get)),
            "bic_1_minus_2": float(bic[1] - bic[2]),
            "largest_gap": float(gaps.max()), "largest_gap_over_spread": float(gaps.max() / spread),
            "largest_gap_at_quantile": float((np.argmax(gaps) + 1) / len(ordered)),
            "bimodal_by_bic": bool(min(bic, key=bic.get) != 1)}


def spearman_with_bootstrap_ci(a: np.ndarray, b: np.ndarray, *, n_boot: int = 2000,
                               seed: int = 42) -> dict:
    """Spearman over AXES, with the CI resampled over axes (the unit of analysis)."""
    from scipy.stats import spearmanr

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    rho = float(spearmanr(a, b).statistic)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(a), len(a))
        if np.std(a[idx]) < 1e-12 or np.std(b[idx]) < 1e-12:
            continue
        draws.append(spearmanr(a[idx], b[idx]).statistic)
    draws = np.asarray(draws, dtype=np.float64)
    return {"rho": rho, "ci95_low": float(np.percentile(draws, 2.5)),
            "ci95_high": float(np.percentile(draws, 97.5)), "n_axes": int(len(a)),
            "n_boot": int(len(draws)),
            "ci_excludes_zero": bool(np.percentile(draws, 2.5) > 0 or np.percentile(draws, 97.5) < 0)}


def verdict(frame: pd.DataFrame, statistic: str, *, quartile: float = 0.75) -> dict:
    """The pre-declared D2.3 decision, applied to one proliferation statistic.

    Thresholds are those fixed in NOTEBOOK_ENTRIES/d3_d2p3_preregistration_
    20260803T1300Z.md before any legibility number existed. `no_verdict` is a
    real outcome: with fewer than 10 legible axes the falsifier has no power and
    must not be reported as a discharge.
    """
    legible = frame["is_legible"].to_numpy(dtype=bool)
    cut = float(frame[statistic].quantile(quartile))
    loaded = (frame[statistic].to_numpy(dtype=np.float64) >= cut)
    n_legible = int(legible.sum())
    share = float(loaded[legible].mean()) if n_legible else float("nan")
    correlation = spearman_with_bootstrap_ci(frame["legibility"].to_numpy(), frame[statistic].to_numpy())
    if n_legible < 10:
        call = "no_verdict_too_few_legible_axes"
    elif share >= 0.90 and correlation["ci_excludes_zero"] and correlation["rho"] > 0:
        call = "FALSIFIER_FIRES_legibility_is_proliferation"
    elif share >= 0.50 or (correlation["ci_excludes_zero"] and abs(correlation["rho"]) >= 0.30):
        call = "PARTIAL"
    else:
        call = "DISCHARGED_legibility_is_not_proliferation"
    return {"statistic": statistic, "verdict": call, "n_legible": n_legible,
            "top_quartile_cut": cut, "share_of_legible_that_are_loaded": share,
            "chance_share": float(1.0 - quartile), "spearman": correlation,
            "n_legible_and_loaded": int((legible & loaded).sum()),
            "n_legible_not_loaded": int((legible & ~loaded).sum())}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--pbs-targets", required=True)
    parser.add_argument("--annotations", required=True, help="per-GENE annotation parquet/csv")
    parser.add_argument("--output", required=True)
    parser.add_argument("--state", default="wsi_biology")
    parser.add_argument("--partition", default="test")
    parser.add_argument("--n-permutations", type=int, default=200)
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    targets = np.load(args.pbs_targets, allow_pickle=True)
    target_ids = np.asarray([str(p) for p in targets["patient_ids"]])
    axis_names = np.asarray([str(t) for t in targets["target_names"]])
    scores = np.asarray(targets["scores"], dtype=np.float64)
    index = {pid: i for i, pid in enumerate(target_ids)}
    annotations = (pd.read_parquet(args.annotations) if args.annotations.endswith(".parquet")
                   else pd.read_csv(args.annotations))
    gene_statistics = axis_gene_statistics(targets["gene_basis"], targets["genes"], annotations,
                                           top_k=args.top_k)
    prol_column = f"prol_top{args.top_k}"

    rows, reports = [], {}
    for artifact_path in args.artifacts:
        path = Path(artifact_path)
        method = path.stem
        raw = np.load(path, allow_pickle=True)
        patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
        cancers = np.asarray([str(c) for c in raw["cancers"]])
        split = np.asarray([str(s) for s in raw["split"]])
        aligned = np.asarray([index.get(pid, -1) for pid in patient_ids])
        mask = (split == args.partition) if args.partition != "all" else np.ones(len(patient_ids), bool)
        mask &= aligned >= 0
        x = np.asarray(raw[args.state], dtype=np.float64)[mask]
        y = scores[aligned[mask]]
        tss, _ = pooled_tissue_source_site(patient_ids[mask], min_site_count=args.min_site_count)
        frame = pd.DataFrame({"cancer": cancers[mask], "tss": tss})
        design = confound_design(frame, ["cancer", "tss"])
        result = per_axis_legibility(x, y, design, frame["cancer"].to_numpy(),
                                     n_permutations=args.n_permutations, seed=args.seed,
                                     n_jobs=args.n_jobs)
        table = gene_statistics.copy()
        table["method"] = method
        table["axis_name"] = axis_names
        table["legibility"] = result["observed"]
        table["null_median"] = result["null_median"]
        table["null_p95"] = result["null_p95"]
        table["excess_over_null_median"] = result["observed"] - result["null_median"]
        table["permutation_p"] = result["permutation_p"]
        table["is_legible"] = result["observed"] > result["null_p95"]
        rows.append(table)

        legibility = table["legibility"].to_numpy()
        ordered = np.sort(legibility)[::-1]
        positive = np.clip(legibility, 0, None)
        reports[method] = {
            "n_patients": int(mask.sum()), "state": args.state, "partition": args.partition,
            "n_axes": int(len(axis_names)), "n_confound_columns": int(design.shape[1]),
            "n_permutations": int(result["n_permutations"]),
            "null_median_over_axes": float(np.median(result["null_median"])),
            "null_p95_over_axes": float(np.median(result["null_p95"])),
            "legibility_quantiles": {q: float(np.percentile(legibility, q)) for q in (0, 5, 25, 50, 75, 95, 100)},
            "n_legible": int(table["is_legible"].sum()),
            "top5_share_of_total": float(ordered[:5].sum() / max(positive.sum(), 1e-12)),
            "top10_share_of_total": float(ordered[:10].sum() / max(positive.sum(), 1e-12)),
            "dominated_by_a_few_axes": bool(ordered[:10].sum() / max(positive.sum(), 1e-12) > 0.5),
            "bimodality": bimodality_report(legibility),
            "verdicts": {statistic: verdict(table, statistic)
                         for statistic in ("prol_wmean", prol_column)},
            "essentiality_check": {statistic: spearman_with_bootstrap_ci(
                table["legibility"].to_numpy(), table[statistic].to_numpy())
                for statistic in ("ess_wmean", f"ess_top{args.top_k}")},
            "annotation_agreement_spearman": spearman_with_bootstrap_ci(
                table["prol_wmean"].to_numpy(), table[prol_column].to_numpy()),
        }
        print(f"[{method}] n_legible={reports[method]['n_legible']}/128 "
              f"top10_share={reports[method]['top10_share_of_total']:.3f} "
              f"verdict[{prol_column}]={reports[method]['verdicts'][prol_column]['verdict']}", flush=True)

    pd.concat(rows, ignore_index=True).to_csv(output / "axis_table.csv", index=False)
    (output / "d2_3_report.json").write_text(json.dumps(reports, indent=2))
    print(f"[done] -> {output}", flush=True)


if __name__ == "__main__":
    main()
