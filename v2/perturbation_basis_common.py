"""The reference CRISPRi response matrix, aligned to a frozen PBS target block.

Every construction that iterates past PBS — a joint basis, a cross-cell-line
consensus basis, a domain-adapted rescaling, or an attribution of an already
legible axis back onto perturbation atoms — needs the *same* object:
``[atom, gene]`` responses in the frozen gene order of ``pbs_targets_*.npz``, on
the same development-fitted scale ``build_pbs_targets`` used.

Re-deriving that alignment per construction is how a comparison silently stops
being apples-to-apples, so it is derived once here. The perturbation file itself
is read through the one E0-validated loader
(:func:`morpheus.v2.calibra.e0_basis_transfer._load_perturbation`); this module
does not open an h5ad.

Two scalings are offered and both are *named*, because the choice is a declared
assumption rather than a fact:

``tcga_sd``
    ``delta / development_gene_sd`` — exactly what ``build_pbs_targets`` does.
    The patient side is unit-variance per gene after standardisation; the
    perturbation side is not.
``own_sd``
    ``delta / across-atom sd of that gene's delta`` — both sides then carry unit
    per-gene variance. This is the domain-adapted alternative; it exists to be
    tested against ``tcga_sd``, not to replace it by assertion.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import numpy as np

__all__ = ["SCALINGS", "load_aligned_response", "subspace_alignment", "atom_folds"]


def atom_folds(n_atoms: int, *, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """A seeded 50/50 split of atom ROW indices into (fold_a, fold_b).

    Lives here, next to the loader, because two callers must agree on it exactly
    and neither may re-roll it: a rotation fitted to maximise cross-cell-line
    agreement is fitted on fold A, and the only cross-line number that may be read
    as evidence for that rotation is computed on fold B. If the fitter and the
    scorer derived the split separately they could drift, and the held-out claim
    would be false without anything failing.
    """
    order = np.random.default_rng(seed).permutation(int(n_atoms))
    cut = int(n_atoms) // 2
    return np.sort(order[:cut]), np.sort(order[cut:])

#: The declared perturbation-side scalings. Named so a deviation cannot be silent.
SCALINGS = ("tcga_sd", "own_sd", "raw")


def _digest(values: np.ndarray) -> str:
    return sha256(np.ascontiguousarray(np.asarray(values, dtype=np.float32)).view(np.uint8)).hexdigest()


def load_aligned_response(perturbation: str | Path, genes, development_scale=None, *,
                          scaling: str = "tcga_sd", min_gene_overlap: int = 1000) -> dict:
    """Return ``[atom, gene]`` responses aligned to ``genes``, plus provenance.

    ``genes`` is the frozen gene order of the PBS block. A perturbation file that
    does not carry every one of them (RPE1 does not: its gene set differs from
    K562's) is aligned on the intersection and the retained gene *positions* are
    returned, so the caller can restrict the patient matrix to the identical
    columns rather than pretending the missing genes were zero.
    """
    if scaling not in SCALINGS:
        raise ValueError(f"unknown perturbation scaling {scaling!r}; expected one of {SCALINGS}")
    from .calibra.e0_basis_transfer import _load_perturbation

    bundle = _load_perturbation(Path(perturbation))
    frozen = [str(g) for g in genes]
    position = {str(gene): index for index, gene in enumerate(bundle.genes)}
    kept_gene_index = np.asarray([i for i, gene in enumerate(frozen) if gene in position], dtype=np.int64)
    if len(kept_gene_index) < min_gene_overlap:
        raise ValueError(f"{Path(perturbation).name} shares only {len(kept_gene_index)} genes with the "
                         f"frozen PBS gene order; that cannot support a stable basis")
    columns = np.asarray([position[frozen[i]] for i in kept_gene_index], dtype=np.int64)
    response = np.asarray(bundle.x, dtype=np.float64)[:, columns]
    if scaling == "tcga_sd":
        if development_scale is None:
            raise ValueError("scaling='tcga_sd' needs the development gene sd the PBS block was built with")
        scale = np.asarray(development_scale, dtype=np.float64)[kept_gene_index]
        n_at_floor = 0
    elif scaling == "own_sd":
        raw_scale = response.std(axis=0)
        floor = float(np.median(raw_scale)) * 1e-3
        n_at_floor = int((raw_scale < floor).sum())
        scale = np.maximum(raw_scale, floor)
    else:
        scale = np.ones(response.shape[1], dtype=np.float64)
        n_at_floor = 0
    response = response / scale[None, :]
    provenance = {
        "perturbation": str(Path(perturbation).resolve()),
        "n_atoms": int(response.shape[0]), "n_genes_aligned": int(response.shape[1]),
        "n_genes_frozen": int(len(frozen)),
        "n_genes_absent_from_perturbation_file": int(len(frozen) - len(kept_gene_index)),
        "scaling": scaling, "n_genes_at_scale_floor": n_at_floor,
        "response_digest": _digest(response), "scale_digest": _digest(scale),
    }
    return {"response": response, "atom_ids": np.asarray([str(a) for a in bundle.row_ids]),
            "gene_index": kept_gene_index, "scale": scale, "provenance": provenance}


def subspace_alignment(a: np.ndarray, b: np.ndarray) -> dict:
    """Principal-angle summary between two column spaces of the same ambient dimension.

    Held-out top-CCA is invariant to any invertible linear reparametrisation of a
    target block, so two blocks spanning the same subspace *must* score identically
    and a "new construction" that only rotates an old one is not one. This reports
    the canonical cosines between the two spans so that claim is checkable rather
    than assumed: ``mean_squared_cosine`` is the normalised subspace overlap
    (1.0 = identical span, 0.0 = orthogonal).
    """
    qa = np.linalg.qr(np.asarray(a, dtype=np.float64))[0]
    qb = np.linalg.qr(np.asarray(b, dtype=np.float64))[0]
    cosines = np.clip(np.linalg.svd(qa.T @ qb, compute_uv=False), 0.0, 1.0)
    return {"n_angles": int(cosines.size),
            "mean_squared_cosine": float(np.mean(cosines ** 2)),
            "max_cosine": float(cosines.max()) if cosines.size else float("nan"),
            "min_cosine": float(cosines.min()) if cosines.size else float("nan"),
            "cosines": cosines.tolist()}
