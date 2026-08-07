"""Alternative multivariate association metrics -- WS-A3, complaint #5.

Why this module exists
-----------------------
Before this module, the entire CALIBRA instrument reads association off two blocks with
exactly one statistic: Pearson correlation on top-CCA projections
(``spectral.top_canonical_correlation`` / ``heldout_top_cca``). A repo-wide grep for HSIC,
distance correlation, RV coefficient and kernel CCA returned zero hits before this file was
written. Complaint #5 ("only top-CCA and single-direction r; no RV/HSIC/kCCA/dCor") is fully
valid, and the predeclared question this module was built to answer
(``NOTEBOOK_ENTRIES/PREDECLARED_ws_a3_association_metrics_20260807T2211Z.md``) is sharper than
"do other metrics also see structure": **the induced correlation described in
``induced_correlation_sweep.py`` and P1 §4 is a projection artifact of Frisch-Waugh-Lovell, a
statement about the RESIDUAL SUBSPACES, not about the statistic used to read association off
them.** If it is real, it must appear under every metric that can detect nonlinear or
multivariate dependence, not only under Pearson-on-a-direction. If HSIC or distance
correlation shows no such floor, the FWL account in the paper is incomplete.

Four metrics, one interface
----------------------------
``rv_coefficient``, ``distance_correlation``, ``hsic``, ``kernel_cca`` all take
``(X, Y) -> float`` on two [n, p] / [n, q] matrices (1-D arrays are accepted and treated as a
single column). None of them duplicates a statistic already defined in ``spectral.py``:
top-CCA is a linear, whitened, capacity-matched maximum over projections; these four are
different mathematical objects (a normalised Frobenius inner product of Gram matrices, a
distance-covariance statistic, a kernel independence statistic, and a *regularised, kernel*
canonical correlation). Where this module needs residualisation or a permutation null it
imports ``residualise.cross_fitted_residuals`` and reuses the pairing-null pattern of
``calibration.permutation_null`` and the floor rule of ``calibration.floors_from_recovery`
rather than re-deriving either.

Kernel and bandwidth rule -- declared in the predeclaration BEFORE any number was measured
--------------------------------------------------------------------------------------------
RBF kernel on both sides for HSIC and kernel CCA, **median heuristic** bandwidth:
``gamma = 1 / (2 * median(pairwise squared Euclidean distance))``, computed independently for
each side from its own subsample. HSIC uses Gretton et al. (2005)'s **biased** V-statistic
estimator, ``HSIC = (1/n^2) * sum(Kc * Lc)`` for doubly-centred kernel matrices ``Kc``, ``Lc``
-- algebraically identical to ``(1/n^2) trace(K H L H)`` for the centring matrix
``H = I - (1/n) 11^T``, verified in ``test_association.py``. Kernel CCA (Bach & Jordan 2002)
solves the regularised problem via a **generalised eigenvalue formulation with
``scipy.linalg.solve``/``eigh``, never ``numpy.linalg.svd``/``svdvals``** -- a legitimate
alternate solution route for kernel CCA (whitening by ``(K + kappa I)`` and taking the leading
eigenvalue of the resulting product is mathematically the same top canonical correlation an SVD
of the whitened cross-kernel would give), chosen specifically so this module needs no entry in
``test_effective_rank_canonical.py``'s ``SVD_ALLOWLIST``. Ridge ``kappa = ridge_scale * n``,
``ridge_scale = 1e-3`` by default (Bach & Jordan's own scaling convention).

Subsampling rule -- declared before any n=6,427-scale run, never applied silently
------------------------------------------------------------------------------------
HSIC, distance correlation and kernel CCA are O(n^2) memory / O(n^2)-O(n^3) compute (RV is not
-- it is computed from p x p / q x q cross-product matrices and needs no subsampling).
:func:`subsample_rows` implements the declared rule: if ``n > 2000``, subsample to 2000 rows via
``np.random.default_rng(seed).choice(n, 2000, replace=False)``, seed fixed and reported. Every
metric that needs it takes ``seed``/``max_n`` and calls this function; nothing truncates rows any
other way.
"""
from __future__ import annotations

import numpy as np

from .calibration import _map, spike_targets
from .residualise import cross_fitted_residuals

__all__ = [
    "subsample_rows", "median_heuristic_gamma",
    "rv_coefficient", "distance_correlation", "hsic", "kernel_cca",
    "ASSOCIATION_METRICS", "compute_all_metrics",
    "metric_permutation_null", "metric_recovery_curve",
    "SUBSAMPLE_THRESHOLD", "SUBSAMPLE_TARGET", "KERNEL_CCA_RIDGE_SCALE",
]

#: Declared subsampling rule (predeclaration §2). Applied by dCor/HSIC/kernel-CCA only.
SUBSAMPLE_THRESHOLD = 2000
SUBSAMPLE_TARGET = 2000
#: Bach & Jordan's own regularisation scaling: ridge = KERNEL_CCA_RIDGE_SCALE * n.
KERNEL_CCA_RIDGE_SCALE = 1e-3


def subsample_rows(n: int, *, seed: int, threshold: int = SUBSAMPLE_THRESHOLD,
                   target: int = SUBSAMPLE_TARGET) -> np.ndarray:
    """Declared, seeded subsampling rule. Never a silent truncation.

    Returns ``arange(n)`` unchanged when ``n <= threshold`` -- the common case for every real
    P1 cohort at the ``test`` partition (n approx 2,530), so no real run this session actually
    exercises the subsample branch; it exists and is unit-tested (``test_association.py``) for
    the ``n=6,427`` (``all`` partition) case a future run may hit.
    """
    if n <= threshold:
        return np.arange(n)
    return np.sort(np.random.default_rng(int(seed)).choice(n, size=int(target), replace=False))


def _as_2d(a) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    return a[:, None] if a.ndim == 1 else a


def _center_cols(a: np.ndarray) -> np.ndarray:
    return a - a.mean(axis=0, keepdims=True)


def _pairwise_sq_dist(a: np.ndarray) -> np.ndarray:
    square = np.sum(a * a, axis=1)
    d2 = square[:, None] + square[None, :] - 2.0 * (a @ a.T)
    np.maximum(d2, 0.0, out=d2)
    return d2


def _double_center(m: np.ndarray) -> np.ndarray:
    """``H M H`` for the centring matrix ``H = I - (1/n) 11^T``, computed without forming H."""
    return m - m.mean(axis=1, keepdims=True) - m.mean(axis=0, keepdims=True) + m.mean()


def median_heuristic_gamma(a: np.ndarray) -> float:
    """RBF bandwidth ``gamma = 1 / (2 * median(pairwise squared distance))``.

    Returns ``1.0`` for a degenerate (all-identical-row) input rather than raising or dividing
    by zero, which only ever arises on a constant synthetic column in a must-fail control.
    """
    a = _as_2d(a)
    n = len(a)
    if n < 2:
        return 1.0
    d2 = _pairwise_sq_dist(a)
    upper = d2[np.triu_indices(n, k=1)]
    median = float(np.median(upper)) if upper.size else 0.0
    return float(1.0 / (2.0 * median)) if median > 1e-300 else 1.0


def _rbf_kernel(a: np.ndarray, gamma: float) -> np.ndarray:
    return np.exp(-float(gamma) * _pairwise_sq_dist(a))


# --------------------------------------------------------------------------------------- #
# the four metrics                                                                          #
# --------------------------------------------------------------------------------------- #

def rv_coefficient(X: np.ndarray, Y: np.ndarray) -> float:
    """RV coefficient (Escoufier 1973, via the Josse & Holmes 2016 review).

    ``RV(X,Y) = trace(Sx Sy) / sqrt(trace(Sx^2) trace(Sy^2))`` for ``Sx = Xc Xc^T``,
    ``Sy = Yc Yc^T``. Computed from the ``p x p`` / ``q x q`` cross-product matrices via
    ``trace(Sx Sy) = ||Xc^T Yc||_F^2`` (a standard trace-cyclic identity), so the cost is
    ``O(n (p^2 + q^2 + pq))`` rather than the naive ``O(n^2 (p+q))`` -- no ``n x n`` matrix is
    ever formed and no subsampling is needed at any cohort size in this project.
    """
    Xc, Yc = _center_cols(_as_2d(X)), _center_cols(_as_2d(Y))
    cross = Xc.T @ Yc
    cxx = Xc.T @ Xc
    cyy = Yc.T @ Yc
    numerator = float(np.sum(cross * cross))
    denominator = float(np.sqrt(np.sum(cxx * cxx) * np.sum(cyy * cyy)))
    return float(numerator / denominator) if denominator > 1e-300 else float("nan")


def distance_correlation(X: np.ndarray, Y: np.ndarray, *, seed: int = 42,
                         max_n: int = SUBSAMPLE_THRESHOLD) -> float:
    """Distance correlation (Szekely, Rizzo & Bakirov 2007), biased V-statistic form.

    ``dCor = sqrt(dCov^2 / sqrt(dVar_X^2 dVar_Y^2))`` from doubly-centred Euclidean distance
    matrices. ``O(n^2)`` memory: subject to :func:`subsample_rows`.
    """
    X, Y = _as_2d(X), _as_2d(Y)
    idx = subsample_rows(len(X), seed=seed, threshold=max_n, target=max_n)
    a = _double_center(np.sqrt(_pairwise_sq_dist(X[idx])))
    b = _double_center(np.sqrt(_pairwise_sq_dist(Y[idx])))
    n2 = float(len(idx)) ** 2
    dcov2 = float(np.sum(a * b)) / n2
    dvarx2 = float(np.sum(a * a)) / n2
    dvary2 = float(np.sum(b * b)) / n2
    denominator = np.sqrt(max(dvarx2, 0.0) * max(dvary2, 0.0))
    if denominator < 1e-300:
        return 0.0
    return float(np.sqrt(max(dcov2, 0.0) / denominator))


def hsic(X: np.ndarray, Y: np.ndarray, *, seed: int = 42, max_n: int = SUBSAMPLE_THRESHOLD,
        gamma_x: float | None = None, gamma_y: float | None = None) -> float:
    """Hilbert-Schmidt Independence Criterion (Gretton et al. 2005), biased estimator.

    ``HSIC = (1/n^2) sum(Kc * Lc)`` for doubly-centred RBF kernel matrices -- algebraically
    ``(1/n^2) trace(K H L H)``, pinned by ``test_association.py``. Median-heuristic bandwidth
    per side unless overridden. ``O(n^2)`` memory: subject to :func:`subsample_rows`.
    """
    X, Y = _as_2d(X), _as_2d(Y)
    idx = subsample_rows(len(X), seed=seed, threshold=max_n, target=max_n)
    Xs, Ys = X[idx], Y[idx]
    gx = float(gamma_x) if gamma_x is not None else median_heuristic_gamma(Xs)
    gy = float(gamma_y) if gamma_y is not None else median_heuristic_gamma(Ys)
    Kc = _double_center(_rbf_kernel(Xs, gx))
    Lc = _double_center(_rbf_kernel(Ys, gy))
    return float(np.sum(Kc * Lc)) / (float(len(idx)) ** 2)


def kernel_cca(X: np.ndarray, Y: np.ndarray, *, seed: int = 42, max_n: int = SUBSAMPLE_THRESHOLD,
              gamma_x: float | None = None, gamma_y: float | None = None,
              ridge_scale: float = KERNEL_CCA_RIDGE_SCALE) -> float:
    """Regularised kernel CCA (Bach & Jordan 2002): top kernel canonical correlation.

    Solves the regularised generalised eigenproblem
    ``[[0, Kx Ky], [Ky Kx, 0]] w = rho [[(Kx+kI)^2, 0], [0, (Ky+kI)^2]] w``,
    ``k = ridge_scale * n``, by whitening with ``Mx = Kx + kI``, ``My = Ky + kI`` (both
    symmetric positive-definite for ``k > 0``) via ``scipy.linalg.solve`` and taking the
    leading eigenvalue of ``T T^T``, ``T = Mx^{-1} Kx Ky My^{-1}``, via ``scipy.linalg.eigh`` --
    the top singular value of ``T`` without ever calling ``numpy.linalg.svd``/``svdvals``
    (see the module docstring for why that is deliberate, not incidental). ``O(n^2)`` memory
    and effectively ``O(n^3)`` compute (two ``n x n`` solves and one ``n x n`` eigh): subject to
    :func:`subsample_rows`, and materially the most expensive of the four metrics at a given n.
    """
    from scipy.linalg import eigh, solve

    X, Y = _as_2d(X), _as_2d(Y)
    idx = subsample_rows(len(X), seed=seed, threshold=max_n, target=max_n)
    Xs, Ys = X[idx], Y[idx]
    n2 = len(idx)
    if n2 < 4:
        return float("nan")
    gx = float(gamma_x) if gamma_x is not None else median_heuristic_gamma(Xs)
    gy = float(gamma_y) if gamma_y is not None else median_heuristic_gamma(Ys)
    Kx = _double_center(_rbf_kernel(Xs, gx))
    Ky = _double_center(_rbf_kernel(Ys, gy))
    kappa = float(ridge_scale) * n2
    Mx = Kx + kappa * np.eye(n2)
    My = Ky + kappa * np.eye(n2)
    A = solve(Mx, Kx, assume_a="pos")           # Mx^{-1} Kx
    B = solve(My, Ky, assume_a="pos")            # My^{-1} Ky
    t = A @ B.T                                  # Mx^{-1} Kx Ky My^{-1}  (since Ky, My symmetric)
    eigenvalues = eigh(t @ t.T, eigvals_only=True)
    top = float(np.sqrt(max(float(eigenvalues[-1]), 0.0)))
    return float(np.clip(top, 0.0, 1.0))


ASSOCIATION_METRICS = {
    "rv": rv_coefficient,
    "dcor": distance_correlation,
    "hsic": hsic,
    "kernel_cca": kernel_cca,
}


def compute_all_metrics(X: np.ndarray, Y: np.ndarray, *, seed: int = 42,
                        max_n: int = SUBSAMPLE_THRESHOLD) -> dict[str, float]:
    """Every metric on the same pair, same seed -- the row a results table actually wants."""
    out = {"rv": rv_coefficient(X, Y)}
    for name in ("dcor", "hsic", "kernel_cca"):
        out[name] = ASSOCIATION_METRICS[name](X, Y, seed=seed, max_n=max_n)
    return out


# --------------------------------------------------------------------------------------- #
# reused instrument: permutation null and the floor, generalised over the metric            #
# --------------------------------------------------------------------------------------- #

def metric_permutation_null(x: np.ndarray, y: np.ndarray, design: np.ndarray, metric, *,
                            strata=None, n_permutations: int = 100, seed: int = 42,
                            n_jobs: int = 1, metric_kwargs: dict | None = None) -> dict:
    """``calibration.permutation_null``, generalised over the association statistic.

    Identical structure to the canonical pairing null: ``x`` is residualised once, ``y`` is
    re-residualised inside every permutation (so any correlation induced by shared
    residualisation is regenerated in the null exactly as it is for the top-CCA channel), rows
    of ``y`` are permuted within ``strata``. Only the final scoring step -- ``metric`` instead
    of ``top_canonical_correlation`` -- differs, and ``cross_fitted_residuals`` is imported, not
    reimplemented.
    """
    metric_kwargs = metric_kwargs or {}
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    rng = np.random.default_rng(seed)
    x_residual = cross_fitted_residuals(x, design, seed=seed)
    strata = np.zeros(len(x), dtype=int) if strata is None else np.asarray(strata)

    observed = float(metric(x_residual, cross_fitted_residuals(y, design, seed=seed), **metric_kwargs))
    orders = []
    for _ in range(n_permutations):
        order = np.arange(len(y))
        for level in np.unique(strata):
            idx = np.flatnonzero(strata == level)
            order[idx] = rng.permutation(idx)
        orders.append((order,))

    def _one(order):
        return float(metric(x_residual, cross_fitted_residuals(y[order], design, seed=seed),
                            **metric_kwargs))

    null = np.asarray(_map(_one, orders, n_jobs), dtype=np.float64) if n_permutations > 0 else np.zeros(0)
    record = {"observed": observed, "n_permutations": int(n_permutations)}
    if n_permutations > 0:
        exceed = int(np.sum(null >= observed))
        record.update({
            "null_median": float(np.median(null)), "null_p95": float(np.percentile(null, 95)),
            "null_max": float(np.max(null)), "excess_over_null_median": float(observed - np.median(null)),
            "permutation_p": float((exceed + 1) / (n_permutations + 1)),
        })
    return record


def metric_recovery_curve(x: np.ndarray, y: np.ndarray, design: np.ndarray, metric, *,
                          levels=(0.0, 0.05, 0.10, 0.20, 0.40), n_draws: int = 10,
                          n_splits: int = 5, alpha: float = 1.0, seed: int = 42,
                          recovery_fraction: float = 0.8, n_jobs: int = 1,
                          metric_kwargs: dict | None = None) -> dict:
    """The calibrated floor (``calibration.spike_recovery_curve``), generalised over the metric.

    A signal of known strength ``r_true`` is planted on one direction pair per draw
    (:func:`calibration.spike_targets`, imported not copied -- identical spike construction the
    shipped floor uses), the *whole spiked block* is pushed through the same
    ``cross_fitted_residuals`` residualisation, and the recovery is read with ``metric`` on the
    two WHOLE residual blocks rather than a single projected direction -- the natural readout
    for a genuinely multivariate statistic, and a materially different (and in principle harder)
    test than the shipped single-direction floor: a whole-block statistic must detect one
    directional signal against the ambient structure of every other column, not read a single
    known axis. The two floors -- ``detection_floor`` (unpaired, conservative) and
    ``transmission_floor`` (paired) -- are computed by ``calibration.floors_from_recovery``,
    reused unchanged so the floor RULE is identical to the shipped instrument's.
    """
    metric_kwargs = metric_kwargs or {}
    from .calibration import floors_from_recovery

    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    levels = np.asarray(levels, dtype=np.float64)
    rng = np.random.default_rng(seed)
    x_residual = cross_fitted_residuals(x, design, n_splits=n_splits, alpha=alpha, seed=seed)

    draw_seeds = [int(rng.integers(1 << 31)) for _ in range(n_draws)]

    def _one_draw(draw_seed):
        column = np.full(len(levels), np.nan)
        for i, level in enumerate(levels):
            spiked = spike_targets(x, y, float(level), rng=np.random.default_rng(draw_seed))
            spiked_residual = cross_fitted_residuals(spiked, design, n_splits=n_splits, alpha=alpha,
                                                     seed=seed)
            column[i] = float(metric(x_residual, spiked_residual, **metric_kwargs))
        return column

    columns = _map(_one_draw, [(s,) for s in draw_seeds], n_jobs)
    recovered = np.column_stack(columns)
    floors = floors_from_recovery(levels, recovered, recovery_fraction=recovery_fraction)
    return {
        "levels": levels.tolist(), "n_draws": int(n_draws),
        "recovered_median": np.nanmedian(recovered, axis=1).tolist(),
        "recovered_p10": np.nanpercentile(recovered, 10, axis=1).tolist(),
        "recovered_p90": np.nanpercentile(recovered, 90, axis=1).tolist(),
        "baseline_median": float(np.nanmedian(recovered[0])),
        **floors,
    }
