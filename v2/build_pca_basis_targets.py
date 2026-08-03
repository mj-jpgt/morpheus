"""Must-beat baseline: an ordinary PCA basis of the expression matrix.

This is the most deflationary alternative in §9 and the one a reviewer reaches
for first: if the top 128 principal components of the same expression matrix are
as morphology-legible as the interventional dictionary, then the dictionary has
contributed nothing that ordinary variance decomposition does not already
supply, and the interventional framing is decoration.

The leak discipline is copied from PBS exactly rather than approximated: the
components are fit on development rows only (``split != "test"``), through
``build_pbs_targets.fit_development_expression_transform``, so no held-out
cancer influences the basis. Fitting PCA on all rows would flatter this baseline
and, worse, would make a *loss* to it unreadable — nobody could tell whether PCA
won because it is genuinely competitive or because it saw the test cancers.

NMF is offered on the same footing (``--basis nmf``). It is a strictly
non-negative parts decomposition, which is the form usually argued to be more
"programme-like" than PCA, so including it removes the objection that the PCA
comparison was chosen because PCA is the weaker opponent.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .baseline_target_common import (development_expression, load_reference_targets,
                                     write_target_block)

__all__ = ["build_pca_basis_targets"]


def build_pca_basis_targets(*, pbs_targets: str, rna_table: str, output: str,
                            n_components: int = 0, basis: str = "pca",
                            seed: int = 0) -> dict[str, object]:
    if basis not in {"pca", "nmf"}:
        raise ValueError("basis must be pca or nmf")
    reference = load_reference_targets(pbs_targets)
    expression, transform = development_expression(rna_table, reference)
    width = int(n_components) if n_components else int(reference["gene_basis"].shape[1])
    development = reference["split"] != "test"
    if basis == "pca":
        from sklearn.decomposition import PCA
        model = PCA(n_components=width, random_state=seed).fit(expression[development])
        scores = model.transform(expression)
        detail = {"explained_variance_ratio_sum": float(model.explained_variance_ratio_.sum()),
                  "explained_variance_ratio_first": float(model.explained_variance_ratio_[0])}
    else:
        from sklearn.decomposition import NMF
        # NMF needs non-negative input; the development-standardised matrix is
        # signed, so it is shifted by the development minimum. The shift is fit
        # on development rows only, exactly like every other transform here.
        shift = float(expression[development].min())
        model = NMF(n_components=width, init="nndsvd", random_state=seed, max_iter=400).fit(
            np.maximum(expression[development] - shift, 0.0))
        scores = model.transform(np.maximum(expression - shift, 0.0))
        detail = {"reconstruction_err": float(model.reconstruction_err_), "development_shift": shift,
                  "n_iter": int(model.n_iter_)}
    prefix = "PCA" if basis == "pca" else "NMF"
    names = np.asarray([f"{prefix}_{i:03d}" for i in range(scores.shape[1])])
    manifest = {"target_kind": f"{basis}_expression_basis_coordinates", "basis": basis,
                "n_components": int(scores.shape[1]), "seed": int(seed),
                "fit_population": "development_rows_only_split_not_test",
                "n_development_rows": int(development.sum()),
                "pbs_targets": str(Path(pbs_targets).resolve()),
                "pbs_manifest_digest": reference["manifest_digest"],
                "gene_count": int(len(reference["genes"])), "gene_digest": reference["gene_digest"],
                "expression_transform": transform,
                "scores_digest": sha256(np.ascontiguousarray(
                    np.asarray(scores, dtype=np.float32)).view(np.uint8)).hexdigest(), **detail}
    return write_target_block(output, reference, scores, names, f"{prefix}_BASIS", manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pbs-targets", required=True)
    parser.add_argument("--rna-table", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-components", type=int, default=0, help="default: match the reference dictionary")
    parser.add_argument("--basis", default="pca", choices=("pca", "nmf"))
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(build_pca_basis_targets(
        pbs_targets=args.pbs_targets, rna_table=args.rna_table, output=args.output,
        n_components=args.n_components, basis=args.basis, seed=args.seed), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
