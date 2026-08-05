"""Rotations of the frozen 128-axis PCA subspace, chosen to be *attributable*.

The causal-attribution certificate passes for 29 of 128 PCA axes. That count was
measured in a basis picked for a reason that has nothing to do with attribution:
PCA axes are variance-maximising and mutually orthogonal because of how PCA is
defined, and no causal programme is obliged to line up with one. So 29/128 could
be a ceiling on the biology, or it could be an artefact of *where the axes point*.

This module answers that by holding the subspace fixed and moving only the basis
inside it. Every rotation here is an orthogonal ``[axis, axis]`` matrix ``R``;
the rotated loadings are ``loadings @ R`` and span **exactly** the same
128 dimensions (``subspace_alignment`` = 1.0, asserted by
:func:`morpheus.v2.causal_attribution.attribution_report` before it scores
anything). Information content is therefore constant across arms and a difference
in the certified count is a fact about basis choice alone.

Four rotations, two of which never see a certificate quantity:

``varimax``
    Classical orthogonal varimax on the gene loadings. Seeks *simple structure* —
    each axis loading on few genes. It is the honest arm: sparsity is a plausible
    prior for "a real causal programme" and is measured without reference to the
    perturbation data at all.
``ica``
    FastICA with symmetric decorrelation on the same loadings. A different
    inductive bias again — non-Gaussian, statistically independent directions
    rather than sparse ones — and equally blind to the certificate.
``r2opt``
    Maximises the mean gene-fold cross-validated R² across axes, subject to
    orthogonality. **A structural warning belongs with this one, not after it.**
    The ridge is linear in the target block, so the rotated out-of-fold residual is
    ``E @ R`` and ``R²_k = 1 - (RᵀAR)_kk / (RᵀBR)_kk`` with ``A = EᵀE``,
    ``B = T_cᵀT_c``. PCA loadings are unit-norm, so ``B ≈ σ²I``, and when ``B`` is
    exactly ``σ²I`` the sum ``Σ_k (RᵀAR)_kk = tr(A)`` is **rotation-invariant**:
    the mean R² cannot be raised by any rotation, only redistributed. Since the R²
    condition already passes 128 of 128, this arm cannot raise the certified count
    through the quantity it optimises. It is here because it is the objective the
    question names, and because a provably-inert objective is part of the answer.
``xline_mean`` / ``xline``
    Maximise K562-vs-RPE1 agreement of the atom-cosine profile — the condition that
    does the actual rejecting (78 of 99 non-certified axes). Both are **circular
    against the plain certificate** and must never be quoted from it: they are
    fitted on one half of the shared atoms (:func:`perturbation_basis_common.atom_folds`)
    and the only number that may be read as evidence for them is the certificate
    scored on the other half, which the report computes for every arm including
    the identity.

    They differ in *what* they maximise, and the difference matters more than it
    looks. ``xline_mean`` maximises the mean per-axis correlation between the two
    lines — and that quantity is very nearly conserved under rotation for the same
    reason the mean R² is: with ``p_k = q_k = 1`` the per-axis correlation is
    ``r_kᵀ M r_k`` and its sum over an orthonormal basis is ``tr(M)``. A rotation
    can move cross-line agreement *between* axes; it cannot manufacture more of
    it. ``xline`` therefore maximises a smooth surrogate for the quantity the
    certificate actually counts — **the number of axes whose agreement clears
    0.30** — which under a conserved total is maximised by *spreading* the
    agreement evenly rather than concentrating it. If a certified count can be
    raised that way and the raise survives on held-out atoms, the 29/128 is a
    statement about the PCA basis; if it collapses on held-out atoms, the
    certificate was merely optimisable.

Both iterative arms ascend on the orthogonal group with a polar retraction, from
several starts, and report the spread across starts: a win smaller than the
optimiser's own scatter is not a win.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .causal_attribution import atom_cosines, cross_line_alignment, development_pca, gene_fold_ridge_r2
from .perturbation_basis_common import atom_folds, load_aligned_response, subspace_alignment

__all__ = ["ROTATIONS", "varimax_rotation", "ica_rotation", "r2_optimising_rotation",
           "cross_line_rotation", "fit_rotation"]

#: The named rotations. ``none`` is the identity and exists so the baseline arm is
#: run through the identical driver rather than trusted from a previous entry.
ROTATIONS = ("none", "varimax", "ica", "r2opt", "xline_mean", "xline")

#: The certificate's cross-line bar. `xline` maximises a smooth count of axes above
#: it, so it is read from the certificate rather than chosen here.
CROSS_LINE_FLOOR = 0.30


def _polar(matrix: np.ndarray) -> np.ndarray:
    """The nearest orthogonal matrix to ``matrix`` — the retraction onto O(n).

    ``U Vᵀ`` from ``M = U S Vᵀ``. This is an SVD used as a *projection*, not as a
    spectrum statistic: nothing is computed from ``S``, which is discarded.
    """
    u, _, vt = np.linalg.svd(np.asarray(matrix, dtype=np.float64), full_matrices=False)
    return u @ vt


def _random_rotation(width: int, seed: int) -> np.ndarray:
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(width, width)))
    return q * np.sign(np.diag(r))[None, :]


def varimax_rotation(loadings: np.ndarray, *, gamma: float = 1.0, max_iter: int = 1000,
                     tol: float = 1e-10) -> dict:
    """Classical orthogonal varimax on ``[gene, axis]`` loadings.

    Maximises the sum over axes of the variance of that axis's squared loadings —
    "few large loadings and many near-zero ones" — over the orthogonal group. The
    update is the standard one (Kaiser 1958): the polar factor of
    ``Lᵀ (B³ - (γ/p) B diag(diag(BᵀB)))`` with ``B = L R``.

    Sees no perturbation data and no certificate quantity whatsoever.
    """
    loadings = np.asarray(loadings, dtype=np.float64)
    n_genes, width = loadings.shape
    rotation = np.eye(width)
    previous = 0.0
    criterion = []
    for iteration in range(max_iter):
        rotated = loadings @ rotation
        gradient = loadings.T @ (rotated ** 3
                                 - (gamma / n_genes) * rotated @ np.diag(np.diag(rotated.T @ rotated)))
        rotation = _polar(gradient)
        value = float(((loadings @ rotation) ** 4).sum()
                      - (gamma / n_genes) * (((loadings @ rotation) ** 2).sum(axis=0) ** 2).sum())
        criterion.append(value)
        if previous != 0.0 and abs(value - previous) < tol * abs(previous):
            break
        previous = value
    return {"rotation": rotation, "n_iterations": iteration + 1, "criterion": criterion[-1],
            "criterion_at_identity": criterion[0] if criterion else float("nan")}


def ica_rotation(loadings: np.ndarray, *, seed: int = 0, max_iter: int = 2000,
                 tol: float = 1e-6) -> dict:
    """FastICA with symmetric decorrelation on the ``[gene, axis]`` loadings.

    The loadings are already an orthonormal basis, so they are white up to the
    ``sqrt(n_genes)`` scale; FastICA is therefore run with whitening **off**, which
    is what makes its unmixing matrix orthogonal and hence a genuine rotation of
    the same subspace rather than a new subspace. The residual departure from
    orthogonality (finite iteration count) is projected away with :func:`_polar`
    and reported, so "it is a rotation" is a measured quantity.
    """
    from sklearn.decomposition import FastICA

    loadings = np.asarray(loadings, dtype=np.float64)
    n_genes, width = loadings.shape
    whitened = loadings * np.sqrt(n_genes)
    model = FastICA(n_components=width, whiten=False, algorithm="parallel", fun="logcosh",
                    max_iter=max_iter, tol=tol, random_state=seed)
    model.fit(whitened)
    candidate = np.asarray(model.components_, dtype=np.float64).T
    error = float(np.abs(candidate.T @ candidate - np.eye(width)).max())
    rotation = _polar(candidate)
    return {"rotation": rotation, "n_iterations": int(getattr(model, "n_iter_", -1)),
            "orthogonality_error_before_projection": error}


def _orthogonal_ascent(objective, gradient, start: np.ndarray, *, max_iter: int = 500,
                       step: float = 1.0, tol: float = 1e-12) -> dict:
    """Projected-gradient ascent on the orthogonal group with a polar retraction.

    The Riemannian gradient at ``R`` is ``R skew(Rᵀ G)``; the retraction is the
    nearest orthogonal matrix to ``R + η T``. Backtracking on ``η`` makes every
    accepted step a genuine increase, so the reported trajectory is monotone and a
    non-improving arm cannot be mistaken for a failed optimiser.
    """
    rotation = np.asarray(start, dtype=np.float64)
    value = float(objective(rotation))
    initial, iteration = value, 0
    for iteration in range(1, max_iter + 1):
        raw = gradient(rotation)
        skew = rotation.T @ raw
        tangent = rotation @ ((skew - skew.T) / 2.0)
        if np.abs(tangent).max() < 1e-14:
            break
        improved = False
        for _ in range(40):
            candidate = _polar(rotation + step * tangent)
            trial = float(objective(candidate))
            if trial > value + tol:
                rotation, value, improved = candidate, trial, True
                step *= 2.0
                break
            step *= 0.5
        if not improved:
            break
    return {"rotation": rotation, "objective": value, "objective_at_start": initial,
            "n_iterations": iteration}


def r2_optimising_rotation(residual: np.ndarray, targets: np.ndarray, *, seeds=(0, 1, 2),
                           max_iter: int = 500) -> dict:
    """Orthogonal rotation maximising the mean gene-fold R² across axes.

    ``residual`` is the out-of-fold residual ``targets - prediction`` at the
    selected alpha, taken from :func:`causal_attribution.gene_fold_ridge_r2`
    itself — the ridge is not re-fitted here, which is the whole reason that
    function grew a ``return_prediction`` flag.
    """
    residual = np.asarray(residual, dtype=np.float64)
    centred = np.asarray(targets, dtype=np.float64)
    centred = centred - centred.mean(axis=0, keepdims=True)
    a_matrix, b_matrix = residual.T @ residual, centred.T @ centred
    width = a_matrix.shape[0]

    def value(rotation):
        a = np.einsum("ik,ij,jk->k", rotation, a_matrix, rotation)
        b = np.einsum("ik,ij,jk->k", rotation, b_matrix, rotation)
        return float(np.mean(1.0 - a / np.maximum(b, 1e-300)))

    def grad(rotation):
        a = np.einsum("ik,ij,jk->k", rotation, a_matrix, rotation)
        b = np.maximum(np.einsum("ik,ij,jk->k", rotation, b_matrix, rotation), 1e-300)
        return (-2.0 / width) * (a_matrix @ rotation / b[None, :]
                                 - b_matrix @ rotation * (a / b ** 2)[None, :])

    runs = [_orthogonal_ascent(value, grad, start, max_iter=max_iter)
            for start in [np.eye(width)] + [_random_rotation(width, s) for s in seeds[1:]]]
    best = max(runs, key=lambda run: run["objective"])
    return {"rotation": best["rotation"], "objective": best["objective"],
            "objective_at_identity": value(np.eye(width)),
            "objective_by_start": [run["objective"] for run in runs],
            "n_iterations_by_start": [run["n_iterations"] for run in runs],
            "trace_invariance_check": {"trace_a": float(np.trace(a_matrix)),
                                       "trace_b": float(np.trace(b_matrix)),
                                       "b_diagonal_relative_spread": float(
                                           np.std(np.diag(b_matrix)) / np.mean(np.diag(b_matrix)))}}


def cross_line_rotation(cosine_a: np.ndarray, cosine_b: np.ndarray, rows: np.ndarray, *,
                        objective: str = "count", floor: float = CROSS_LINE_FLOOR,
                        tau: float = 0.05, seeds=(0, 1, 2), max_iter: int = 500) -> dict:
    """Orthogonal rotation maximising K562-vs-RPE1 atom-profile agreement.

    ``cosine_a`` / ``cosine_b`` are the ``[atom, axis]`` unnormalised cosine
    profiles of the two cell lines, restricted to ``rows`` — the fitting half of
    the shared atoms.

    ``objective='mean'`` maximises the mean per-axis **Pearson** correlation.
    ``objective='count'`` maximises ``mean_k sigmoid((rho_k - floor)/tau)``, a
    smooth count of axes clearing the certificate's bar. The certificate scores
    **Spearman**, so both are surrogates for the graded quantity even before the
    circularity is considered; the fold-B certificate is the only number either
    may be judged on.
    """
    if objective not in ("mean", "count"):
        raise ValueError(f"unknown cross-line objective {objective!r}")
    a = np.asarray(cosine_a, dtype=np.float64)[rows]
    b = np.asarray(cosine_b, dtype=np.float64)[rows]
    a = a - a.mean(axis=0, keepdims=True)
    b = b - b.mean(axis=0, keepdims=True)
    cross = (a.T @ b + b.T @ a) / 2.0
    saa, sbb = a.T @ a, b.T @ b
    width = cross.shape[0]

    def correlation(rotation):
        """Per-axis cross-line Pearson correlation, and its per-column gradient."""
        m = np.einsum("ik,ij,jk->k", rotation, cross, rotation)
        p = np.maximum(np.einsum("ik,ij,jk->k", rotation, saa, rotation), 1e-300)
        q = np.maximum(np.einsum("ik,ij,jk->k", rotation, sbb, rotation), 1e-300)
        root = np.sqrt(p * q)
        direction = (2.0 * (cross @ rotation) / root[None, :]
                     - (saa @ rotation * (m * q / root ** 3)[None, :]
                        + sbb @ rotation * (m * p / root ** 3)[None, :]))
        return m / root, direction

    def value(rotation):
        rho, _ = correlation(rotation)
        if objective == "mean":
            return float(np.mean(rho))
        return float(np.mean(1.0 / (1.0 + np.exp(-(rho - floor) / tau))))

    def grad(rotation):
        rho, direction = correlation(rotation)
        if objective == "mean":
            return direction / width
        sigmoid = 1.0 / (1.0 + np.exp(-(rho - floor) / tau))
        return direction * (sigmoid * (1.0 - sigmoid) / tau)[None, :] / width

    runs = [_orthogonal_ascent(value, grad, start, max_iter=max_iter)
            for start in [np.eye(width)] + [_random_rotation(width, s) for s in seeds[1:]]]
    best = max(runs, key=lambda run: run["objective"])
    rho_identity, _ = correlation(np.eye(width))
    rho_best, _ = correlation(best["rotation"])
    return {"rotation": best["rotation"], "objective": best["objective"],
            "objective_at_identity": value(np.eye(width)),
            "objective_by_start": [run["objective"] for run in runs],
            "n_iterations_by_start": [run["n_iterations"] for run in runs],
            "objective_kind": objective, "floor": float(floor), "tau": float(tau),
            # The conservation the docstring claims, measured on the real matrices
            # rather than asserted: if these two are equal the mean agreement is a
            # budget the rotation can only redistribute.
            "mean_correlation_at_identity": float(np.mean(rho_identity)),
            "mean_correlation_at_best": float(np.mean(rho_best)),
            "n_above_floor_at_identity": int((rho_identity >= floor).sum()),
            "n_above_floor_at_best": int((rho_best >= floor).sum()),
            "n_fitting_atoms": int(len(rows))}


def fit_rotation(*, method: str, pbs_targets: str, rna_table: str, perturbation: str, output: str,
                 secondary_perturbation: str = "", pca_targets: str = "", n_components: int = 0,
                 seed: int = 0, max_iter: int = 500) -> dict:
    """Fit one named rotation and write it, with its diagnostics, to ``output``."""
    if method not in ROTATIONS:
        raise ValueError(f"unknown rotation {method!r}; expected one of {ROTATIONS}")
    block = development_pca(pbs_targets, rna_table, n_components=n_components, pca_targets=pca_targets)
    loadings, width = block["loadings"], block["width"]
    diagnostics: dict = {"method": method, "width": int(width), "seed": int(seed),
                         "pca_verification": block["verification"]}

    if method == "none":
        fitted = {"rotation": np.eye(width)}
    elif method == "varimax":
        fitted = varimax_rotation(loadings)
    elif method == "ica":
        fitted = ica_rotation(loadings, seed=seed)
    elif method == "r2opt":
        primary = load_aligned_response(perturbation, block["genes"], block["scale"], scaling="tcga_sd")
        design = primary["response"].T
        targets = loadings[primary["gene_index"]]
        ridge = gene_fold_ridge_r2(design @ design.T, targets, seed=seed, return_prediction=True)
        diagnostics["ridge_selected_alpha"] = ridge["selected_alpha"]
        diagnostics["mean_r2_at_identity"] = float(ridge["r2"].mean())
        fitted = r2_optimising_rotation(targets - ridge["prediction"], targets, max_iter=max_iter)
    else:
        if not secondary_perturbation:
            raise ValueError("the cross-line rotation needs --secondary-perturbation")
        primary = load_aligned_response(perturbation, block["genes"], block["scale"], scaling="tcga_sd")
        secondary = load_aligned_response(secondary_perturbation, block["genes"], block["scale"],
                                          scaling="tcga_sd")
        aligned = cross_line_alignment(primary, secondary, loadings)
        norms = np.linalg.norm(aligned["directions"], axis=0)
        cosine_a = atom_cosines(aligned["response_a"], aligned["directions"]) * norms[None, :]
        cosine_b = atom_cosines(aligned["response_b"], aligned["directions"]) * norms[None, :]
        fold_a, fold_b = atom_folds(len(aligned["shared_atoms"]), seed=seed)
        diagnostics["n_shared_atoms"] = int(len(aligned["shared_atoms"]))
        diagnostics["n_atoms_fold_a"] = int(len(fold_a))
        diagnostics["n_atoms_fold_b"] = int(len(fold_b))
        fitted = cross_line_rotation(cosine_a, cosine_b, fold_a, max_iter=max_iter,
                                     objective="mean" if method == "xline_mean" else "count")

    rotation = np.asarray(fitted.pop("rotation"), dtype=np.float64)
    diagnostics.update({key: value for key, value in fitted.items()})
    diagnostics["max_abs_orthogonality_error"] = float(np.abs(rotation.T @ rotation - np.eye(width)).max())
    diagnostics["mean_squared_cosine_vs_pca_span"] = subspace_alignment(
        loadings, loadings @ rotation)["mean_squared_cosine"]
    diagnostics["mean_abs_diagonal"] = float(np.abs(np.diag(rotation)).mean())
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, rotation=rotation, method=np.asarray(method))
    path.with_suffix(".json").write_text(json.dumps(diagnostics, indent=2, sort_keys=True, default=float)
                                         + "\n", encoding="utf-8")
    return diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", required=True, choices=list(ROTATIONS))
    parser.add_argument("--pbs-targets", required=True)
    parser.add_argument("--rna-table", required=True)
    parser.add_argument("--perturbation", required=True)
    parser.add_argument("--secondary-perturbation", default="")
    parser.add_argument("--pca-targets", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-components", type=int, default=0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-iter", type=int, default=500)
    args = parser.parse_args()
    print(json.dumps(fit_rotation(
        method=args.method, pbs_targets=args.pbs_targets, rna_table=args.rna_table,
        perturbation=args.perturbation, output=args.output,
        secondary_perturbation=args.secondary_perturbation, pca_targets=args.pca_targets,
        n_components=args.n_components, seed=args.seed, max_iter=args.max_iter),
        indent=2, sort_keys=True, default=float))


if __name__ == "__main__":
    main()
