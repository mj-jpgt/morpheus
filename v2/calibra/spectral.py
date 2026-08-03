"""Spectral utilities for CALIBRA: effective rank and cross-fitted CCA.

Single source of truth for `effective_rank` (previously duplicated in
`v2/tests/test_stress_collapse.py` and `v2/run_rank_ablation.py`).
"""
from __future__ import annotations

import numpy as np

__all__ = ["effective_rank", "cca_spectrum", "top_canonical_correlation", "heldout_top_cca",
           "heldout_single_direction_correlation"]


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


def _whiten_map(a: np.ndarray, n_components: int):
    """Return (centre, W) so that ``(a - centre) @ W`` has orthonormal columns."""
    centre = a.mean(axis=0, keepdims=True)
    centred = a - centre
    _, s, vt = np.linalg.svd(centred, full_matrices=False)
    keep = s > (s.max() * 1e-10) if s.size and s.max() > 0 else np.zeros_like(s, dtype=bool)
    k = int(min(n_components, int(keep.sum())))
    if k == 0:
        return centre, np.zeros((a.shape[1], 0))
    return centre, (vt[:k].T / s[:k])


def heldout_top_cca(x: np.ndarray, y: np.ndarray, *, n_components: int = 32,
                    seed: int = 42, train_fraction: float = 0.5) -> float:
    """Top canonical correlation with directions FIT ON TRAIN, SCORED ON HELD-OUT rows.

    In-sample CCA is a multivariate maximum and is badly upward-biased at finite n
    — quoting it is like reporting a score on questions you already saw. Here the
    canonical directions are estimated on one split and the correlation is measured
    on the other, giving an unbiased (and much smaller) absolute number.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = len(x)
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    cut = int(n * train_fraction)
    train, test = order[:cut], order[cut:]
    if len(train) < 10 or len(test) < 10:
        return float("nan")

    cx, wx = _whiten_map(x[train], n_components)
    cy, wy = _whiten_map(y[train], n_components)
    if wx.shape[1] == 0 or wy.shape[1] == 0:
        return float("nan")
    ux, uy = (x[train] - cx) @ wx, (y[train] - cy) @ wy
    a, _, bt = np.linalg.svd(ux.T @ uy, full_matrices=False)

    px = ((x[test] - cx) @ wx) @ a[:, 0]
    py = ((y[test] - cy) @ wy) @ bt[0]
    if px.std() < 1e-12 or py.std() < 1e-12:
        return float("nan")
    return float(abs(np.corrcoef(px, py)[0, 1]))


def heldout_single_direction_correlation(x: np.ndarray, y: np.ndarray, *, n_splits: int = 5,
                                         alpha: float = 1.0, seed: int = 42) -> float:
    """Out-of-fold correlation between ONE fitted image direction and ONE target column.

    Why this statistic and not a top-CCA. The CALIBRA ``detection_floor`` is
    expressed in *single-direction* correlation units — it is the smallest
    ``corr(X_res u, Y_res v)`` the pipeline reliably resolves for a fixed
    direction pair. Grading a per-target negative control against that floor
    therefore requires a per-target statistic on the same scale. A canonical
    correlation between a 256-column X and a 1-column y is a maximum over 256
    in-sample directions and is inflated by capacity alone: on pure noise at
    n=2,530 it is nowhere near zero, which would make *every* random control
    "clear the floor" and turn a working negative control into a fake failure.

    Here the direction is a ridge predictor fit on the training folds only and
    scored on the rows it never saw, so it is a single direction chosen without
    reference to the values it is graded on — the same contract ``heldout_top_cca``
    applies to the multivariate case, restricted to one target column.

    SIGNED, for the reason given in ``calibration._correlation``: an
    out-of-fold predictor that anti-correlates with its target is not evidence
    that the target is legible, and taking ``|r|`` would score it as if it were.
    """
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).ravel()
    if x.ndim != 2 or len(x) != len(y) or len(x) < 2 * n_splits:
        return float("nan")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or np.std(y) < 1e-12:
        return float("nan")
    prediction = np.empty(len(y), dtype=np.float64)
    splitter = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for train_idx, test_idx in splitter.split(x):
        model = Ridge(alpha=alpha, fit_intercept=True).fit(x[train_idx], y[train_idx])
        prediction[test_idx] = model.predict(x[test_idx])
    if np.std(prediction) < 1e-12:
        return float("nan")
    return float(np.corrcoef(prediction, y)[0, 1])
