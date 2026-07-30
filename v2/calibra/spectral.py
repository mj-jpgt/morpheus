"""Spectral utilities for CALIBRA: effective rank and cross-fitted CCA.

Single source of truth for `effective_rank` (previously duplicated in
`v2/tests/test_stress_collapse.py` and `v2/run_rank_ablation.py`).
"""
from __future__ import annotations

import numpy as np

__all__ = ["effective_rank", "cca_spectrum", "top_canonical_correlation"]


def effective_rank(x) -> float:
    """Roy-Vetterli effective rank: exp(entropy of L1-normalised singular values).

    Accepts a numpy array or anything exposing ``.detach().cpu().numpy()``
    (e.g. a torch tensor). Column-centres first. Returns 0.0 for a constant batch.
    """
    if hasattr(x, "detach"):
        x = x.detach().cpu().numpy()
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean(axis=0, keepdims=True)
    singular = np.linalg.svd(x, compute_uv=False)
    singular = singular[singular > 1e-12]
    if singular.size == 0:
        return 0.0
    p = singular / singular.sum()
    return float(np.exp(-(p * np.log(p)).sum()))


def _whiten(a: np.ndarray, n_components: int, eps: float = 1e-8):
    """Centre + PCA-whiten to at most ``n_components`` directions."""
    a = np.asarray(a, dtype=np.float64)
    a = a - a.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(a, full_matrices=False)
    keep = s > (s.max() * 1e-10) if s.size and s.max() > 0 else np.zeros_like(s, dtype=bool)
    k = int(min(n_components, int(keep.sum())))
    if k == 0:
        return np.zeros((a.shape[0], 0)), vt[:0], s[:0]
    return u[:, :k], vt[:k], s[:k]


def cca_spectrum(x: np.ndarray, y: np.ndarray, *, n_components: int = 32) -> np.ndarray:
    """Canonical correlations between x and y after PCA-whitening both sides.

    Whitening to a fixed component budget is what keeps the statistic comparable
    across encoders of very different width (a 1536-d encoder is otherwise
    trivially advantaged over a 256-d one).
    """
    ux, _, _ = _whiten(x, n_components)
    uy, _, _ = _whiten(y, n_components)
    if ux.shape[1] == 0 or uy.shape[1] == 0:
        return np.zeros(0)
    # Columns of ux/uy are orthonormal, so the cross-product singular values ARE
    # the canonical correlations.
    singular = np.linalg.svd(ux.T @ uy, compute_uv=False)
    return np.clip(singular, 0.0, 1.0)


def top_canonical_correlation(x: np.ndarray, y: np.ndarray, *, n_components: int = 32) -> float:
    spectrum = cca_spectrum(x, y, n_components=n_components)
    return float(spectrum[0]) if spectrum.size else float("nan")
