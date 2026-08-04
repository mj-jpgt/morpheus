"""Construction 1: attribute an already-legible axis back to perturbation atoms.

The interventional dictionary lost a legibility contest to ordinary PCA of the
same expression matrix. This module inverts the direction of supervision instead
of re-entering that contest: take the axes PCA *already* found legible and ask,
per axis, **which single-gene CRISPRi signatures reproduce this gene-space
direction**. A PCA component cannot supply a causal name; the perturbation
resource can, and that is a different contribution from "our basis is more
legible", which is falsified.

Four quantities per axis, each with its own null, because "axis 46 looks like
perturbing gene G" is worth exactly nothing unless all four hold:

``r2_cv``
    Gene-fold cross-validated ridge reconstruction of the axis's gene loading
    from the 8,403 atom signatures. Cross-validated **over genes**: with more
    atoms than genes the in-sample fit is exactly 1.0 and says nothing. Graded
    against Haar-random gene-space directions put through the identical
    procedure — if a random direction is reconstructable too, the atom set simply
    spans gene space and no axis is attributed.
``top_atom`` / ``atom_cosine``
    The ranked cosine between the axis loading and each atom's response. This is
    the causal *name*. Graded against a **gene-label shuffle** of the response
    matrix: the name must not survive permuting which gene each loading refers
    to, or it is not a name.
``legibility``
    :func:`calibra.spectral.heldout_single_direction_correlation` between the
    residualised image states and that single axis on the held-out partition —
    the per-axis counterpart of the block-level comparison, on the same confound
    design.
``cross_line_spearman``
    The same atom ranking recomputed from RPE1 rather than K562. An attribution
    that does not replicate across cell lines is a K562 artefact.

The headline this module exists to produce is the cross-tab: *do the axes that
are morphology-legible get a strong, stable causal attribution, and which do not?*
It does not need to win a legibility contest to answer that.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .baseline_target_common import development_expression_moments, load_reference_targets
from .perturbation_basis_common import load_aligned_response

__all__ = ["RIDGE_ALPHAS", "gene_fold_ridge_r2", "atom_cosines", "attribution_report"]

#: Predeclared ridge path for the gene-fold reconstruction. Fixed before the run.
RIDGE_ALPHAS = (1e-1, 1.0, 1e1, 1e2, 1e3, 1e4)


def _digest(values: np.ndarray) -> str:
    return sha256(np.ascontiguousarray(np.asarray(values, dtype=np.float32)).view(np.uint8)).hexdigest()


def gene_fold_ridge_r2(gram: np.ndarray, targets: np.ndarray, *, alphas=RIDGE_ALPHAS,
                       n_folds: int = 5, seed: int = 0) -> dict:
    """Out-of-fold R² of reconstructing each column of ``targets`` from the atom span.

    ``gram`` is ``X Xᵀ`` for the ``[gene, atom]`` design ``X``; passing the Gram
    rather than ``X`` is what makes the three arms (true axes, random directions,
    gene-label shuffle) share one 8-hundred-gigaflop matrix product instead of
    three, and — more importantly — makes them share it *exactly*, so an arm
    cannot differ from another through a different design matrix.

    Folds are over **genes**. The dual (kernel) ridge is used because there are
    more atoms than genes; it is the same estimator as the primal form, not an
    approximation of it.

    One alpha is selected, shared across all target columns, as the one maximising
    the mean out-of-fold R². Every arm selects its own alpha under the identical
    rule and the full alpha curve is returned, so the selection cannot advantage
    one arm over another.
    """
    gram = np.asarray(gram, dtype=np.float64)
    targets = np.asarray(targets, dtype=np.float64)
    n_genes = gram.shape[0]
    if gram.shape != (n_genes, n_genes) or targets.shape[0] != n_genes:
        raise ValueError("gene-fold ridge needs a [gene, gene] Gram and [gene, target] targets")
    folds = np.array_split(np.random.default_rng(seed).permutation(n_genes), n_folds)
    curve = np.zeros((len(alphas), targets.shape[1]), dtype=np.float64)
    for index, alpha in enumerate(alphas):
        prediction = np.zeros_like(targets)
        for fold in folds:
            test = np.asarray(fold, dtype=np.int64)
            train = np.setdiff1d(np.arange(n_genes), test, assume_unique=False)
            k_train = gram[np.ix_(train, train)] + float(alpha) * np.eye(len(train))
            dual = np.linalg.solve(k_train, targets[train] - targets[train].mean(axis=0, keepdims=True))
            prediction[test] = gram[np.ix_(test, train)] @ dual + targets[train].mean(axis=0, keepdims=True)
        residual = ((targets - prediction) ** 2).sum(axis=0)
        total = ((targets - targets.mean(axis=0, keepdims=True)) ** 2).sum(axis=0)
        curve[index] = 1.0 - residual / np.maximum(total, 1e-300)
    best = int(np.argmax(curve.mean(axis=1)))
    return {"alphas": [float(a) for a in alphas], "r2_curve": curve, "selected_alpha": float(alphas[best]),
            "selected_index": best, "r2": curve[best].copy(), "n_folds": int(n_folds), "seed": int(seed)}


def atom_cosines(response: np.ndarray, directions: np.ndarray) -> np.ndarray:
    """``[atom, direction]`` cosine between each atom's response and each direction."""
    response = np.asarray(response, dtype=np.float64)
    directions = np.asarray(directions, dtype=np.float64)
    atom_norm = np.maximum(np.linalg.norm(response, axis=1, keepdims=True), 1e-300)
    direction_norm = np.maximum(np.linalg.norm(directions, axis=0, keepdims=True), 1e-300)
    return (response / atom_norm) @ (directions / direction_norm)


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import spearmanr
    if len(a) < 3:
        return float("nan")
    return float(spearmanr(a, b).statistic)


def _spearman_with_interval(a: np.ndarray, b: np.ndarray, *, n_boot: int = 2000,
                            seed: int = 0) -> dict:
    """Spearman plus a bootstrap interval over the resampled *axes*.

    The 128 axes are the sample here, not the patients: the question is whether
    the across-axis association between legibility and attribution is bigger than
    the noise in a set of 128 axes.
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    finite = np.isfinite(a) & np.isfinite(b)
    a, b = a[finite], b[finite]
    rng = np.random.default_rng(seed)
    draws = [_spearman(a[d], b[d]) for d in (rng.integers(0, len(a), size=len(a)) for _ in range(n_boot))]
    draws = np.asarray(draws, dtype=float)
    return {"spearman": _spearman(a, b), "n": int(len(a)),
            "ci95_low": float(np.nanpercentile(draws, 2.5)),
            "ci95_high": float(np.nanpercentile(draws, 97.5))}


def _partial_spearman(a: np.ndarray, b: np.ndarray, control: np.ndarray) -> float:
    """Spearman between ``a`` and ``b`` with ``control`` linearly removed, in rank space.

    The confound this exists for, named in the predeclaration as the first thing
    to distrust: a high-variance PCA axis carries more signal, so it is *both*
    easier to see from morphology *and* easier to reconstruct from anything. If
    the legibility-attribution association is only the axis's own signal-to-noise
    ratio talking, it disappears once the explained-variance rank is partialled out.
    """
    from scipy.stats import rankdata
    finite = np.isfinite(a) & np.isfinite(b) & np.isfinite(control)
    if finite.sum() < 4:
        return float("nan")
    ranks = [rankdata(np.asarray(v, dtype=float)[finite]) for v in (a, b, control)]
    design = np.column_stack([np.ones(int(finite.sum())), ranks[2]])
    residual = [v - design @ np.linalg.lstsq(design, v, rcond=None)[0] for v in ranks[:2]]
    if residual[0].std() < 1e-12 or residual[1].std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(residual[0], residual[1])[0, 1])


def attributed_set_coherence(response: np.ndarray, cosine: np.ndarray, *, n_top: int = 10,
                             n_random: int = 200, seed: int = 0) -> dict:
    """Is an axis's top-attributed atom set more mutually coherent than a random set?

    A per-axis attribution that names ten unrelated genes is a ranking artefact; one
    that names ten members of the same machine is a claim about a mechanism. The
    statistic is the mean absolute pairwise correlation among the top ``n_top``
    atoms' response vectors, expressed as a percentile of the same quantity for
    ``n_random`` size-matched random atom sets. It needs no external complex
    annotation, which is the point — none is on disk.
    """
    response = np.asarray(response, dtype=np.float64)
    centred = response - response.mean(axis=1, keepdims=True)
    unit = centred / np.maximum(np.linalg.norm(centred, axis=1, keepdims=True), 1e-300)
    rng = np.random.default_rng(seed)

    def coherence(rows: np.ndarray) -> float:
        block = unit[rows]
        gram = np.abs(block @ block.T)
        upper = gram[np.triu_indices(len(rows), k=1)]
        return float(upper.mean()) if upper.size else float("nan")

    null = np.asarray([coherence(rng.choice(len(unit), size=n_top, replace=False))
                       for _ in range(n_random)])
    observed, percentile = [], []
    for k in range(cosine.shape[1]):
        rows = np.argsort(-np.abs(cosine[:, k]))[:n_top]
        value = coherence(rows)
        observed.append(value)
        percentile.append(float((null < value).mean()))
    return {"coherence": np.asarray(observed), "percentile": np.asarray(percentile),
            "null_median": float(np.median(null)), "n_top": int(n_top), "n_random": int(n_random)}


def _perturbed_gene(atom_id: str) -> str:
    """``10005_ZBTB4_P1_ENSG00000174282`` -> ``ZBTB4``. The causal name itself."""
    parts = str(atom_id).split("_")
    return parts[1] if len(parts) > 2 else str(atom_id)


def _axis_legibility(artifacts, reference, scores, states, *, partition: str, seed: int) -> dict:
    """Per-axis held-out single-direction correlation, on the block harness's footing."""
    from .calibra.residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
    from .calibra.spectral import heldout_single_direction_correlation

    out: dict[str, dict[str, list[float]]] = {}
    index = {pid: i for i, pid in enumerate(reference["patient_ids"])}
    for artifact in artifacts:
        raw = np.load(artifact, allow_pickle=True)
        patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
        cancers = np.asarray([str(c) for c in raw["cancers"]])
        split = np.asarray([str(s) for s in raw["split"]])
        aligned = np.asarray([index.get(pid, -1) for pid in patient_ids])
        mask = aligned >= 0
        if partition != "all":
            mask &= split == partition
        tss, _ = pooled_tissue_source_site(patient_ids[mask])
        design = confound_design(pd.DataFrame({"cancer": cancers[mask], "tss": tss}), ["cancer", "tss"])
        target = cross_fitted_residuals(scores[aligned[mask]], design, seed=seed)
        for state in states:
            x = cross_fitted_residuals(np.asarray(raw[state], dtype=np.float64)[mask], design, seed=seed)
            out.setdefault(Path(artifact).stem, {})[state] = [
                float(heldout_single_direction_correlation(x, target[:, k], seed=seed))
                for k in range(target.shape[1])]
    return out


def attribution_report(*, pbs_targets: str, rna_table: str, perturbation: str, output: str,
                       secondary_perturbation: str = "", pca_targets: str = "",
                       artifacts: tuple[str, ...] = (), states: tuple[str, ...] = ("wsi_biology",),
                       n_components: int = 0, partition: str = "test", n_top_atoms: int = 10,
                       seed: int = 0) -> dict:
    reference = load_reference_targets(pbs_targets)
    expression, transform, _mean, scale = development_expression_moments(rna_table, reference)
    genes = [str(g) for g in reference["genes"]]
    width = int(n_components) if n_components else int(reference["gene_basis"].shape[1])
    development = reference["split"] != "test"

    from sklearn.decomposition import PCA
    pca = PCA(n_components=width, random_state=0).fit(expression[development])
    loadings = np.asarray(pca.components_, dtype=np.float64).T     # [gene, axis]
    axis_scores = np.asarray(pca.transform(expression), dtype=np.float32)
    verification = {"pca_targets": str(Path(pca_targets).resolve()) if pca_targets else "",
                    "scores_digest": _digest(axis_scores), "matches_frozen_pca_block": None}
    if pca_targets:
        frozen = np.load(pca_targets, allow_pickle=True)
        frozen_scores = np.asarray(frozen["scores"], dtype=np.float32)
        verification["matches_frozen_pca_block"] = bool(
            frozen_scores.shape == axis_scores.shape
            and float(np.abs(frozen_scores - axis_scores).max()) < 1e-3)
        verification["max_abs_difference_vs_frozen"] = float(np.abs(frozen_scores - axis_scores).max())

    primary = load_aligned_response(perturbation, genes, scale, scaling="tcga_sd")
    gene_index, response, atom_ids = primary["gene_index"], primary["response"], primary["atom_ids"]
    design = response.T                                             # [gene, atom]
    gram = design @ design.T
    targets = loadings[gene_index]

    rng = np.random.default_rng(seed)
    random_directions = np.linalg.qr(rng.normal(size=(len(gene_index), width)))[0]
    shuffle_order = rng.permutation(len(gene_index))

    arms = {"axes": gene_fold_ridge_r2(gram, targets, seed=seed),
            "random_direction_null": gene_fold_ridge_r2(gram, random_directions, seed=seed),
            "gene_label_shuffle_null": gene_fold_ridge_r2(gram, targets[shuffle_order], seed=seed)}

    cosine = atom_cosines(response, targets)
    cosine_shuffled = atom_cosines(response, targets[shuffle_order])
    cosine_random = atom_cosines(response, random_directions)

    cross_line: dict[str, object] = {"status": "not_requested"}
    cross_line_spearman = np.full(width, np.nan)
    if secondary_perturbation:
        secondary = load_aligned_response(secondary_perturbation, genes, scale, scaling="tcga_sd")
        shared_genes = np.intersect1d(gene_index, secondary["gene_index"])
        shared_atoms = np.intersect1d(atom_ids, secondary["atom_ids"])
        primary_order, secondary_order = np.argsort(atom_ids), np.argsort(secondary["atom_ids"])
        rows_a = primary_order[np.searchsorted(atom_ids, shared_atoms, sorter=primary_order)]
        rows_b = secondary_order[np.searchsorted(secondary["atom_ids"], shared_atoms, sorter=secondary_order)]
        cols_a = np.searchsorted(gene_index, shared_genes)
        cols_b = np.searchsorted(secondary["gene_index"], shared_genes)
        directions = loadings[shared_genes]
        cosine_a = atom_cosines(response[np.ix_(rows_a, cols_a)], directions)
        cosine_b = atom_cosines(secondary["response"][np.ix_(rows_b, cols_b)], directions)
        cross_line_spearman = np.asarray([_spearman(cosine_a[:, k], cosine_b[:, k]) for k in range(width)])
        cross_line = {"status": "scored", "n_shared_atoms": int(len(shared_atoms)),
                      "n_shared_genes": int(len(shared_genes)),
                      "secondary": secondary["provenance"]}

    shuffle_spearman = np.asarray([_spearman(cosine[:, k], cosine_shuffled[:, k]) for k in range(width)])
    coherence = attributed_set_coherence(response, cosine, n_top=n_top_atoms, seed=seed)
    coherence_null = attributed_set_coherence(response, cosine_random, n_top=n_top_atoms, seed=seed)
    legibility = _axis_legibility(artifacts, reference, axis_scores, states,
                                  partition=partition, seed=seed) if artifacts else {}

    records = []
    for k in range(width):
        order = np.argsort(-np.abs(cosine[:, k]))[:n_top_atoms]
        record = {"axis": f"PCA_{k:03d}", "axis_index": k,
                  "explained_variance_ratio": float(pca.explained_variance_ratio_[k]),
                  "r2_cv": float(arms["axes"]["r2"][k]),
                  "r2_cv_random_direction_null": float(arms["random_direction_null"]["r2"][k]),
                  "r2_cv_gene_label_shuffle_null": float(arms["gene_label_shuffle_null"]["r2"][k]),
                  "top_atom_cosine": float(np.abs(cosine[:, k]).max()),
                  "top_atom_cosine_shuffle_null": float(np.abs(cosine_shuffled[:, k]).max()),
                  "top_atom_cosine_random_direction_null": float(np.abs(cosine_random[:, k]).max()),
                  "shuffle_rank_spearman": float(shuffle_spearman[k]),
                  "cross_line_rank_spearman": float(cross_line_spearman[k]),
                  "attributed_set_coherence": float(coherence["coherence"][k]),
                  "attributed_set_coherence_percentile": float(coherence["percentile"][k]),
                  "attributed_set_coherence_random_direction_null": float(coherence_null["coherence"][k]),
                  "top_atoms": ";".join(str(atom_ids[i]) for i in order),
                  "top_perturbed_genes": ";".join(_perturbed_gene(atom_ids[i]) for i in order)}
        for artifact_name, per_state in legibility.items():
            for state, values in per_state.items():
                record[f"legibility__{artifact_name}__{state}"] = float(values[k])
        records.append(record)
    table = pd.DataFrame(records)

    summary = {
        "n_axes": int(width), "n_atoms": int(len(atom_ids)), "n_genes": int(len(gene_index)),
        "pbs_targets": str(Path(pbs_targets).resolve()),
        "perturbation": primary["provenance"], "expression_transform": transform,
        "pca_verification": verification, "cross_line": cross_line,
        "ridge": {arm: {"selected_alpha": value["selected_alpha"], "alphas": value["alphas"],
                        "mean_r2_by_alpha": value["r2_curve"].mean(axis=1).round(6).tolist(),
                        "median_r2": float(np.median(value["r2"]))}
                  for arm, value in arms.items()},
        "random_direction_null_median_r2": float(np.median(arms["random_direction_null"]["r2"])),
        "attribution_statistic_usable": bool(np.median(arms["random_direction_null"]["r2"]) <= 0.5),
        "shuffle_rank_spearman_median": float(np.nanmedian(shuffle_spearman)),
        "cross_line_rank_spearman_median": float(np.nanmedian(cross_line_spearman)),
        "states": list(states), "partition": partition, "seed": int(seed),
    }
    summary["attributed_set_coherence"] = {
        "median": float(np.median(coherence["coherence"])),
        "random_atom_set_null_median": coherence["null_median"],
        "random_direction_axes_median": float(np.median(coherence_null["coherence"])),
        "fraction_of_axes_above_95th_percentile_of_random_atom_sets": float(
            (coherence["percentile"] > 0.95).mean()),
        "n_top": coherence["n_top"], "n_random_sets": coherence["n_random"]}
    variance = table["explained_variance_ratio"].to_numpy(dtype=float)
    r2 = table["r2_cv"].to_numpy(dtype=float)
    cosine_column = table["top_atom_cosine"].to_numpy(dtype=float)
    summary["explained_variance_confound"] = {
        "spearman_r2_vs_explained_variance": _spearman(r2, variance),
        "spearman_top_cosine_vs_explained_variance": _spearman(cosine_column, variance)}
    for artifact_name, per_state in legibility.items():
        for state in per_state:
            column = f"legibility__{artifact_name}__{state}"
            finite = table[column].to_numpy(dtype=float)
            summary[f"legibility_vs_r2_spearman__{artifact_name}__{state}"] = _spearman_with_interval(
                finite, r2, seed=seed)
            summary[f"legibility_vs_top_cosine_spearman__{artifact_name}__{state}"] = _spearman_with_interval(
                finite, cosine_column, seed=seed)
            # The predeclared distrust check: is the association only the axis's own
            # signal-to-noise ratio, which raises legibility and reconstructability alike?
            summary[f"legibility_vs_r2_partial_spearman__{artifact_name}__{state}"] = _partial_spearman(
                finite, r2, variance)
            summary[f"legibility_vs_explained_variance_spearman__{artifact_name}__{state}"] = _spearman(
                finite, variance)

    directory = Path(output)
    directory.mkdir(parents=True, exist_ok=True)
    table.to_csv(directory / "axis_attribution.csv", index=False)
    (directory / "attribution_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n",
                                                        encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pbs-targets", required=True)
    parser.add_argument("--rna-table", required=True)
    parser.add_argument("--perturbation", required=True)
    parser.add_argument("--secondary-perturbation", default="")
    parser.add_argument("--pca-targets", default="", help="frozen PCA block; verifies the axes attributed are the ones that beat PBS")
    parser.add_argument("--artifacts", nargs="*", default=[])
    parser.add_argument("--states", nargs="*", default=["wsi_biology"])
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-components", type=int, default=0)
    parser.add_argument("--partition", default="test")
    parser.add_argument("--n-top-atoms", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(attribution_report(
        pbs_targets=args.pbs_targets, rna_table=args.rna_table, perturbation=args.perturbation,
        output=args.output, secondary_perturbation=args.secondary_perturbation,
        pca_targets=args.pca_targets, artifacts=tuple(args.artifacts), states=tuple(args.states),
        n_components=args.n_components, partition=args.partition, n_top_atoms=args.n_top_atoms,
        seed=args.seed), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
