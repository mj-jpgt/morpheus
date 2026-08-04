"""Post-PBS interventional constructions, on the PBS block's own footing.

`t11_t12_must_beat_baselines_20260803T0440Z` measured the interventional
dictionary losing to ordinary PCA of the same expression matrix in 3 of 4 cells.
`decision_iterate_past_pbs_and_p4_negatives_20260804T2230Z` reads that
structurally rather than as a verdict on the framing, and names three further
*constructions* to test with the identical harness. This module builds them, all
128 columns, all bound by ``baseline_target_common`` to the same cohort, split,
gene order and development-only expression transform as the block they lost with.

``--construction pbs_rebuild``
    The control, and the reason to trust the other three: rebuild the frozen PBS
    block through this module's code path. It must reproduce
    ``pbs_targets_k128_v2.npz`` to floating-point noise. If it does not, every
    other construction here is measuring a re-derivation bug rather than a
    construction.
``--construction joint_cca``
    A *joint* basis instead of an OR. With ``P`` the PBS gene directions and ``Q``
    the development-fit PCA gene loadings, ``PᵀQ = U S Vᵀ`` is the CCA between the
    two covariance row spaces and ``S`` are the canonical cosines; the joint basis
    is the equal-weight blend ``normalise(P u_k + Q v_k)``. Every column then has
    a declared component in the interventional span and a declared component in
    the cohort-variance span.
``--construction consensus``
    A denoised interventional basis: only atoms whose response reproduces across
    K562 and RPE1 enter the SVD, retained against a **mismatched-pair null**
    rather than a bare correlation threshold, because the two response matrices
    share global structure that would make any fixed threshold look impressive.
    Separates "the causal framing is wrong" from "this resource is noisy".
``--construction domain_adapted``
    ``build_pbs_targets`` divides the CRISPRi delta by the TCGA development gene
    SD, which leaves the patient side unit-variance per gene and the perturbation
    side not. Here the perturbation matrix is z-scored by its own across-atom SD
    so both sides share a scale. The current choice is declared in
    ``fit_development_expression_transform``; this is the test of it.

The SVD itself is always ``pbs.ReferenceDictionary.fit`` — the same primitive PBS
uses — so what differs between these blocks is the atom set, the gene scaling or
the basis blend, never a second implementation of the decomposition.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .baseline_target_common import (development_expression_moments, load_reference_targets,
                                     write_target_block)
from .pbs import ReferenceDictionary
from .perturbation_basis_common import load_aligned_response, subspace_alignment

__all__ = ["CONSTRUCTIONS", "build_causal_basis_targets", "cross_line_consensus_atoms"]

CONSTRUCTIONS = ("pbs_rebuild", "joint_cca", "consensus", "domain_adapted")


def _digest(values: np.ndarray) -> str:
    return sha256(np.ascontiguousarray(np.asarray(values, dtype=np.float32)).view(np.uint8)).hexdigest()


def _row_correlation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Per-row Pearson correlation between two aligned ``[row, column]`` matrices."""
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    denominator = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
    return np.where(denominator > 0, np.einsum("ij,ij->i", a, b) / np.maximum(denominator, 1e-300), 0.0)


def cross_line_consensus_atoms(primary: np.ndarray, secondary: np.ndarray, *,
                               percentile: float = 95.0, n_null: int = 20000,
                               seed: int = 0) -> dict:
    """Retain atoms whose response reproduces across two cell lines.

    ``primary`` and ``secondary`` are ``[atom, gene]`` on the SAME atom order and
    the SAME gene order. The bar is the ``percentile`` of a **mismatched-pair**
    null — the correlation between atom ``a`` in one line and a different atom
    ``b`` in the other. A bare threshold on ``r`` would be graded against zero,
    which is the wrong reference: both response matrices carry the same global
    transcriptional structure, so unrelated atoms already correlate.
    """
    if primary.shape != secondary.shape or primary.ndim != 2 or len(primary) < 2:
        raise ValueError("consensus needs two aligned [atom, gene] matrices with >=2 atoms")
    matched = _row_correlation(primary, secondary)
    rng = np.random.default_rng(seed)
    left = rng.integers(0, len(primary), size=n_null)
    right = rng.integers(0, len(primary), size=n_null)
    distinct = left != right
    null = _row_correlation(primary[left[distinct]], secondary[right[distinct]])
    threshold = float(np.percentile(null, percentile))
    retained = matched > threshold
    return {"matched": matched, "null": null, "threshold": threshold, "retained": retained,
            "n_retained": int(retained.sum()), "n_atoms": int(len(matched)),
            "percentile": float(percentile), "n_null_pairs": int(distinct.sum()),
            "null_median": float(np.median(null)), "matched_median": float(np.median(matched))}


def build_causal_basis_targets(*, pbs_targets: str, rna_table: str, perturbation: str,
                               output: str, construction: str, n_components: int = 0,
                               secondary_perturbation: str = "", consensus_percentile: float = 95.0,
                               seed: int = 0) -> dict[str, object]:
    if construction not in CONSTRUCTIONS:
        raise ValueError(f"unknown construction {construction!r}; expected one of {CONSTRUCTIONS}")
    reference = load_reference_targets(pbs_targets)
    expression, transform, _mean, scale = development_expression_moments(rna_table, reference)
    genes = [str(g) for g in reference["genes"]]
    width = int(n_components) if n_components else int(reference["gene_basis"].shape[1])
    development = reference["split"] != "test"
    detail: dict[str, object] = {}

    if construction == "joint_cca":
        from sklearn.decomposition import PCA
        pca = PCA(n_components=width, random_state=seed).fit(expression[development])
        q = np.asarray(pca.components_, dtype=np.float64).T
        p = np.asarray(reference["gene_basis"], dtype=np.float64)
        left, cosines, right_t = np.linalg.svd(p.T @ q, full_matrices=False)
        blended = p @ left + q @ right_t.T
        basis = blended / np.maximum(np.linalg.norm(blended, axis=0, keepdims=True), 1e-12)
        scores = expression @ basis
        detail = {"canonical_cosines": np.asarray(cosines, dtype=float).round(6).tolist(),
                  "canonical_cosine_mean": float(np.mean(cosines)),
                  "canonical_cosine_top": float(cosines[0]),
                  "alignment_vs_pbs_span": subspace_alignment(basis, p),
                  "alignment_vs_pca_span": subspace_alignment(basis, q),
                  "pbs_vs_pca_span": subspace_alignment(p, q),
                  "explained_variance_ratio_sum": float(pca.explained_variance_ratio_.sum()),
                  "basis_digest": _digest(basis), "n_genes_used": int(len(genes))}
        prefix, group = "JOINTCCA", "JOINT_CCA_BASIS"
    else:
        scaling = "own_sd" if construction == "domain_adapted" else "tcga_sd"
        primary = load_aligned_response(perturbation, genes, scale, scaling=scaling)
        gene_index = primary["gene_index"]
        response, atom_ids = primary["response"], primary["atom_ids"]
        detail["primary_perturbation"] = primary["provenance"]
        if construction == "consensus":
            if not secondary_perturbation:
                raise ValueError("consensus needs --secondary-perturbation (the second cell line)")
            secondary = load_aligned_response(secondary_perturbation, genes, scale, scaling=scaling)
            shared_genes = np.intersect1d(gene_index, secondary["gene_index"])
            primary_columns = np.searchsorted(gene_index, shared_genes)
            secondary_columns = np.searchsorted(secondary["gene_index"], shared_genes)
            shared_atoms = np.intersect1d(atom_ids, secondary["atom_ids"])
            if len(shared_atoms) < width * 2:
                raise ValueError(f"only {len(shared_atoms)} atoms are shared across the two cell lines; "
                                 f"a {width}-component consensus basis would be rank-starved")
            primary_order = np.argsort(atom_ids)
            primary_rows = primary_order[np.searchsorted(atom_ids, shared_atoms, sorter=primary_order)]
            secondary_order = np.argsort(secondary["atom_ids"])
            secondary_rows = secondary_order[np.searchsorted(secondary["atom_ids"], shared_atoms,
                                                             sorter=secondary_order)]
            a = response[np.ix_(primary_rows, primary_columns)]
            b = secondary["response"][np.ix_(secondary_rows, secondary_columns)]
            consensus = cross_line_consensus_atoms(a, b, percentile=consensus_percentile, seed=seed)
            retained = consensus["retained"]
            if int(retained.sum()) < width * 2:
                raise ValueError(f"consensus retained only {int(retained.sum())} atoms at percentile "
                                 f"{consensus_percentile}; a {width}-component basis would be rank-starved")
            response = 0.5 * (a[retained] + b[retained])
            atom_ids = shared_atoms[retained]
            gene_index = shared_genes
            norms = np.linalg.norm(a, axis=1)
            detail.update({
                "secondary_perturbation": secondary["provenance"],
                "n_shared_atoms": int(len(shared_atoms)), "n_shared_genes": int(len(shared_genes)),
                "n_retained_atoms": int(retained.sum()),
                "retention_threshold": consensus["threshold"],
                "matched_correlation_median": consensus["matched_median"],
                "mismatched_null_median": consensus["null_median"],
                "retention_percentile": consensus["percentile"],
                # Untrustworthiness check (b) of the predeclaration: is retention just
                # selecting the loudest atoms?
                "retention_vs_response_norm_pointbiserial": float(
                    np.corrcoef(retained.astype(float), norms)[0, 1]),
                "retained_response_effective_rank": None})
        gene_names = [genes[i] for i in gene_index]
        dictionary = ReferenceDictionary.fit(response, gene_names, atom_ids.tolist(), n_components=width)
        scores = dictionary.encode_expression(expression[:, gene_index], gene_names)
        basis = np.asarray(dictionary.gene_basis, dtype=np.float64)
        if construction == "consensus":
            from .calibra.spectral import effective_rank
            detail["retained_response_effective_rank"] = float(effective_rank(response))
        embedded = np.zeros((len(genes), basis.shape[1]))
        embedded[gene_index] = basis
        detail.update({"n_atoms_used": int(response.shape[0]), "n_genes_used": int(len(gene_index)),
                       "basis_digest": _digest(basis),
                       "singular_values": np.asarray(dictionary.singular_values, dtype=float).round(6).tolist(),
                       "alignment_vs_pbs_span": subspace_alignment(
                           embedded, np.asarray(reference["gene_basis"], dtype=np.float64))})
        prefix = {"pbs_rebuild": "PBSREBUILD", "consensus": "CONSENSUS",
                  "domain_adapted": "DOMADAPT"}[construction]
        group = {"pbs_rebuild": "PBS_REBUILD", "consensus": "CONSENSUS_DICTIONARY",
                 "domain_adapted": "DOMAIN_ADAPTED_DICTIONARY"}[construction]

    names = np.asarray([f"{prefix}_{i:03d}" for i in range(scores.shape[1])])
    manifest = {"target_kind": f"post_pbs_{construction}_coordinates", "construction": construction,
                "n_components": int(scores.shape[1]), "seed": int(seed),
                "fit_population": "development_rows_only_split_not_test",
                "n_development_rows": int(development.sum()),
                "pbs_targets": str(Path(pbs_targets).resolve()),
                "pbs_manifest_digest": reference["manifest_digest"],
                "gene_count": int(len(genes)), "gene_digest": reference["gene_digest"],
                "expression_transform": transform,
                "scores_digest": _digest(scores), **detail}
    return write_target_block(output, reference, scores, names, group, manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pbs-targets", required=True)
    parser.add_argument("--rna-table", required=True)
    parser.add_argument("--perturbation", default="", help="E0-validated Perturb-seq h5ad (not needed for joint_cca)")
    parser.add_argument("--secondary-perturbation", default="", help="second cell line, for --construction consensus")
    parser.add_argument("--construction", required=True, choices=CONSTRUCTIONS)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-components", type=int, default=0, help="default: match the reference dictionary")
    parser.add_argument("--consensus-percentile", type=float, default=95.0)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(build_causal_basis_targets(
        pbs_targets=args.pbs_targets, rna_table=args.rna_table, perturbation=args.perturbation,
        output=args.output, construction=args.construction, n_components=args.n_components,
        secondary_perturbation=args.secondary_perturbation,
        consensus_percentile=args.consensus_percentile, seed=args.seed), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
