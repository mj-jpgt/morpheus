"""How much of the morphology→molecular channel do the *certified* axes carry?

The attribution certificate names 29 of 128 PCA axes. "29 of 128 axes have a
causal name" and "the named axes are the ones a promptable interface needs" are
different claims, and only the first has ever been measured. This module measures
the second: the canonical channel statistic, restricted to the 29 certified
columns, against a size-matched random draw of 29 of the 128 and against the 29
most legible columns.

The statistic is not restated here. It is
``calibra.spectral.heldout_top_cca`` on ``calibra.residualise.cross_fitted_residuals``
of both sides under ``confound_design(cancer, pooled_tissue_source_site)`` —
which is exactly ``run_calibra._channel_measurement``'s ``channel_statistic`` and
exactly what ``p1_evidence/baseline_paired_bootstrap.py`` scores every target
block with. Only the set of target columns changes.

Two things about the statistic that shape how the table must be read, stated here
rather than discovered later:

* it whitens the target side to ``n_components`` directions, so a 29-column block
  read at ``k=16`` is seen through its own top-16 principal subspace only. The
  table is therefore reported across ``k``, with ``k=29`` the value that uses every
  selected column;
* it is a correlation, not a variance decomposition, so ``channel(subset) /
  channel(all)`` is a ratio of top canonical correlations and is **not** a share
  of anything additive. It is reported as a ratio and named as one.

The random null is drawn two ways: uniformly over the 128, and **stratified on the
explained-variance decile** so the draw matches the certified set's variance
profile. The certified axes are skewed towards low-variance directions and the
channel statistic favours high-variance ones, so the unstratified null and the
stratified null answer different questions — "is this set special at all" and "is
this set special for a reason other than where it sits in the variance ordering".
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .calibra.residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from .calibra.spectral import heldout_top_cca

__all__ = ["SUBSET_SIZE", "load_target_block", "subset_channel", "channel_share_report"]

#: The certified count. The random draw is matched to it exactly; a comparison at
#: any other size is a different question.
SUBSET_SIZE = 29


def load_target_block(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(patient_ids, scores, target_names)`` with the random-control columns dropped.

    The drop rule is `baseline_paired_bootstrap._targets`'s: a block may carry
    ``RANDOM_CONTROL__`` columns which are not part of the block being graded, and
    a channel measured with them in is not the channel that block was scored on.
    """
    raw = np.load(Path(path), allow_pickle=True)
    ids = np.asarray([str(p) for p in raw["patient_ids"]])
    names = np.asarray([str(t) for t in raw["target_names"]])
    keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
    return ids, np.asarray(raw["scores"], dtype=np.float64)[:, keep], names[keep]


def subset_channel(x: np.ndarray, y: np.ndarray, columns, *, n_components: int, seed: int) -> float:
    """The canonical channel statistic on a column subset of an already-residualised block.

    ``cross_fitted_residuals`` fits one ridge per target column against a shared
    design, so residualising the whole block and then selecting columns is exactly
    residualising the selected columns — asserted in
    ``v2/tests/test_channel_share.py`` rather than assumed, because if it were
    false every number here would be measured on a different residual from the
    one the block comparison uses.
    """
    return float(heldout_top_cca(x, y[:, np.asarray(columns, dtype=np.int64)],
                                 n_components=n_components, seed=seed))


def _decile_strata(values: np.ndarray, n_strata: int = 10) -> np.ndarray:
    order = np.argsort(-np.asarray(values, dtype=float))
    strata = np.empty(len(values), dtype=np.int64)
    strata[order] = np.minimum((np.arange(len(values)) * n_strata) // len(values), n_strata - 1)
    return strata


def _stratified_draw(rng: np.random.Generator, strata: np.ndarray, target: np.ndarray) -> np.ndarray:
    """A draw with the same per-stratum counts as ``target``, disjointness respected."""
    picked = []
    for stratum in np.unique(strata):
        pool = np.where(strata == stratum)[0]
        want = int((strata[target] == stratum).sum())
        picked.append(rng.choice(pool, size=min(want, len(pool)), replace=False))
    return np.concatenate(picked) if picked else np.zeros(0, dtype=np.int64)


def channel_share_report(*, pca_targets: str, attribution_table: str, artifacts: tuple[str, ...],
                         states: tuple[str, ...], output: str, partition: str = "test",
                         components: tuple[int, ...] = (4, 8, 16, 29), n_random: int = 200,
                         seed: int = 42) -> dict:
    table = pd.read_csv(attribution_table)
    certified = np.where(table["causal_name_certified"].to_numpy(bool))[0]
    variance = table["explained_variance_ratio"].to_numpy(float)
    strata = _decile_strata(variance)
    patient_ids, scores, names = load_target_block(pca_targets)
    if scores.shape[1] != len(table):
        raise ValueError(f"block has {scores.shape[1]} columns, attribution table has {len(table)}")

    rows: list[dict] = []
    summary: dict = {"n_certified": int(len(certified)), "subset_size": SUBSET_SIZE,
                     "certified_axes": table.loc[table["causal_name_certified"], "axis"].tolist(),
                     "certified_axis_indices": certified.tolist(),
                     "median_certified_axis_index": float(np.median(certified)),
                     "median_all_axis_index": float(np.median(np.arange(len(table)))),
                     "pca_targets": str(Path(pca_targets).resolve()),
                     "attribution_table": str(Path(attribution_table).resolve()),
                     "n_random": int(n_random), "seed": int(seed), "partition": partition,
                     "statistic": f"heldout_top_cca(k=*,seed={seed}) on cancer+pooled-TSS residuals",
                     "equivalence_check": {}}
    if len(certified) != SUBSET_SIZE:
        summary["warning_subset_size"] = f"certified count is {len(certified)}, not {SUBSET_SIZE}"

    index = {pid: i for i, pid in enumerate(patient_ids)}
    for artifact in artifacts:
        raw = np.load(artifact, allow_pickle=True)
        artifact_ids = np.asarray([str(p) for p in raw["patient_ids"]])
        cancers = np.asarray([str(c) for c in raw["cancers"]])
        split = np.asarray([str(s) for s in raw["split"]])
        aligned = np.asarray([index.get(pid, -1) for pid in artifact_ids])
        mask = aligned >= 0
        if partition != "all":
            mask &= split == partition
        tss, _ = pooled_tissue_source_site(artifact_ids[mask])
        design = confound_design(pd.DataFrame({"cancer": cancers[mask], "tss": tss}), ["cancer", "tss"])
        block = scores[aligned[mask]]
        y = cross_fitted_residuals(block, design, seed=seed)
        # The claim that lets every subset below be read off one residualisation.
        direct = cross_fitted_residuals(block[:, certified], design, seed=seed)
        summary["equivalence_check"][Path(artifact).stem] = float(np.abs(direct - y[:, certified]).max())

        for state in states:
            x = cross_fitted_residuals(np.asarray(raw[state], dtype=np.float64)[mask], design, seed=seed)
            legibility_column = f"legibility__{Path(artifact).stem}__{state}"
            legible = (np.argsort(-table[legibility_column].to_numpy(float))[:SUBSET_SIZE]
                       if legibility_column in table.columns else np.zeros(0, dtype=np.int64))
            rng = np.random.default_rng(seed)
            uniform = [rng.choice(len(table), size=SUBSET_SIZE, replace=False) for _ in range(n_random)]
            stratified = [_stratified_draw(rng, strata, certified) for _ in range(n_random)]
            for k in components:
                null = np.asarray([subset_channel(x, y, draw, n_components=k, seed=seed)
                                   for draw in uniform])
                null_stratified = np.asarray([subset_channel(x, y, draw, n_components=k, seed=seed)
                                              for draw in stratified])
                value = subset_channel(x, y, certified, n_components=k, seed=seed)
                whole = subset_channel(x, y, np.arange(len(table)), n_components=k, seed=seed)
                legible_value = (subset_channel(x, y, legible, n_components=k, seed=seed)
                                 if len(legible) else float("nan"))
                row = {"artifact": Path(artifact).stem, "state": state, "n_components": int(k),
                       "channel_all_128": whole, "channel_certified": value,
                       "channel_most_legible": legible_value,
                       "ratio_certified_over_all": value / whole if whole else float("nan"),
                       "random_median": float(np.median(null)),
                       "random_p05": float(np.percentile(null, 5)),
                       "random_p95": float(np.percentile(null, 95)),
                       "random_iqr": float(np.percentile(null, 75) - np.percentile(null, 25)),
                       "certified_percentile_in_random_null": float((null < value).mean()),
                       "stratified_random_median": float(np.median(null_stratified)),
                       "certified_percentile_in_stratified_null": float((null_stratified < value).mean()),
                       "ratio_most_legible_over_all": (legible_value / whole
                                                       if len(legible) and whole else float("nan")),
                       "n_patients": int(mask.sum()), "n_confound_columns": int(design.shape[1])}
                rows.append(row)
                print(f"[{row['artifact']}::{state} k={k}] certified {value:.4f} "
                      f"random med {row['random_median']:.4f} (pct {row['certified_percentile_in_random_null']:.2f}) "
                      f"legible {row['channel_most_legible']:.4f} all128 {whole:.4f}", flush=True)

    frame = pd.DataFrame(rows)
    directory = Path(output)
    directory.mkdir(parents=True, exist_ok=True)
    frame.to_csv(directory / "channel_share.csv", index=False)
    summary["rows"] = rows
    (directory / "channel_share_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=float) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pca-targets", required=True)
    parser.add_argument("--attribution-table", required=True)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--states", nargs="*", default=["wsi_biology"])
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition", default="test")
    parser.add_argument("--components", nargs="*", type=int, default=[4, 8, 16, 29])
    parser.add_argument("--n-random", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    summary = channel_share_report(
        pca_targets=args.pca_targets, attribution_table=args.attribution_table,
        artifacts=tuple(args.artifacts), states=tuple(args.states), output=args.output,
        partition=args.partition, components=tuple(args.components), n_random=args.n_random,
        seed=args.seed)
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=2, sort_keys=True,
                     default=float))


if __name__ == "__main__":
    main()
