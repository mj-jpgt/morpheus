"""Spectral utilities for CALIBRA: effective rank and cross-fitted CCA.

**Single source of truth for `effective_rank`.** This module holds the ONLY
implementation in the repository. Three mutually incompatible statistics were
previously implemented under that one name (`spectral.py`, `d1_audit.py`,
`d1_geometry_probe.py`), plus byte-identical torch duplicates in
`run_rank_ablation.py`, `test_stress_collapse.py` and `e0_basis_transfer.py` and
a fourth variant inline in `training.py`. They are now all this function, and
`v2/tests/test_effective_rank_canonical.py` fails if a second definition
reappears anywhere in the tree.
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

__all__ = ["effective_rank", "RankVariant", "CANONICAL", "RANK_VARIANTS", "cca_spectrum",
           "top_canonical_correlation", "heldout_top_cca", "heldout_cca_projection",
           "heldout_top_cca_indexed", "paired_absolute_correlation",
           "heldout_single_direction_correlation", "heldout_cca_directions"]


class RankVariant(NamedTuple):
    """A named preprocessing + statistic choice, so a deviation cannot be silent."""
    centre: bool
    normalise_rows: bool
    order: int

    @property
    def label(self) -> str:
        return (f"{'centred' if self.centre else 'uncentred'}"
                f"{'+rownorm' if self.normalise_rows else ''}"
                f"|order{self.order}")


#: The canonical definition. Every call site uses this unless it names a deviation.
#:
#: Both source papers were re-verified at full text on 2026-08-05 (Roy & Vetterli
#: from the Zenodo PDF of DOI 10.5281/zenodo.40328, pp. 606-610; RankMe from
#: arXiv:2210.02885v3), and the verification changes what can be claimed:
#:
#: - ``order=1`` -> Roy & Vetterli Definition 1 exactly, ``erank(A) = exp{H(p)}``,
#:   ``H = -sum p_k log p_k``, ``p_k = sigma_k / ||sigma||_1``, natural log, and
#:   Property 1 ``1 <= erank(A) <= rank(A) <= Q``. This is the published statistic
#:   and the one RankMe adopts; it is the only variant comparable to any number
#:   outside this repository. The other two statistics in this repository are
#:   order-2 Hill numbers of the same spectrum and are NOT effective rank.
#: - ``centre=True`` -> the SVD is taken of the COLUMN-CENTRED matrix. **This is a
#:   deliberate deviation and must be stated wherever a value is quoted.** Roy &
#:   Vetterli take the SVD of the raw matrix A ("a complex-valued non-all-zero
#:   matrix A of size M x N whose singular value decomposition is given by
#:   A = UDV*") with no preprocessing anywhere in the paper; RankMe is silent on
#:   centring and applies none. We centre because for a representation the column
#:   mean is a rank-1 component carrying no between-sample information, and the
#:   collapse family documented on this project is ``z_i = m + a_i * u`` (one shared
#:   offset plus one direction), whose uncentred effective rank is ~2 and whose
#:   centred effective rank is ~1. Uncentred rank charges the mean vector as a
#:   dimension. It is also what the majority of this project's quoted values used,
#:   so centring is the choice that keeps the historical record recomputable.
#: - ``normalise_rows=False`` -> rows are used at their own norms. Row L2
#:   normalisation appears in neither source paper, discards norm variation that is
#:   part of the representation, and changes the statistic's response to collapse.
#:
#: Not reconciled, and no number here is comparable to a published RankMe value:
#: RankMe uses ``p_k = sigma_k/||sigma||_1 + eps`` with eps OUTSIDE the division
#: (verified at glyph coordinates in the v3 PDF), so its p_k sum to
#: ``1 + min(N,K)*eps`` rather than to 1 and its statistic is not exp of a Shannon
#: entropy; the paper never states the eps used in that equation. Roy & Vetterli
#: use no eps at all, only the convention ``0 log 0 = 0``.
CANONICAL = RankVariant(centre=True, normalise_rows=False, order=1)

#: The three historical statistics of `paper/P2_RANK_DRAFT.md` §3.1, expressed as
#: variants of the one function, so that they can be measured side by side.
RANK_VARIANTS = {
    "R1": CANONICAL,                                                   # spectral.py
    "R2": RankVariant(centre=True, normalise_rows=False, order=2),     # d1_audit.py
    "R3": RankVariant(centre=True, normalise_rows=True, order=2),      # d1_geometry_probe.py, training.py
    "R1_uncentred": RankVariant(centre=False, normalise_rows=False, order=1),
    "R1_rownorm": RankVariant(centre=True, normalise_rows=True, order=1),
    "R2_uncentred": RankVariant(centre=False, normalise_rows=False, order=2),
}

#: Singular values are retained when ``sigma > sigma_max * max(n, p) * eps`` -- the
#: standard LAPACK/numpy relative rank cut at float64, matching the argument already
#: made at `e0_basis_transfer.py:_full_dictionary_rank`. Roy & Vetterli adopt the
#: ``0 log 0 = 0`` convention and set no tolerance; RankMe adds ``eps`` OUTSIDE the
#: normalisation, ``p_k = sigma_k/||sigma||_1 + eps``. This implementation previously
#: used an ABSOLUTE cut of 1e-12, which breaks the scale invariance the statistic
#: otherwise has -- multiplying the matrix by 1e-9 changed the answer. No value in
#: this repository is comparable to a published RankMe number under any of the three.
EPS = float(np.finfo(np.float64).eps)

#: A matrix whose post-preprocessing spectrum is this small relative to its own
#: pre-preprocessing scale is degenerate and scores 0.0. This replaces the old
#: absolute 1e-12 (which was implicitly relative to O(1) data) with the same number
#: made scale-invariant, so that a representation collapsed to within float noise
#: still reports 0 rather than the dimensionality of the noise.
DEGENERACY_TOLERANCE = 1e-12


def _finalise(singular, order: int, exp, log, total, square_sum):
    """Shared tail: L1-normalise the retained spectrum and take the Hill number."""
    if order == 1:
        p = singular / total
        return float(exp(-(p * log(p)).sum()))
    if order == 2:
        # (sum sigma)^2 / sum sigma^2 == 1 / sum p^2, the order-2 Hill number of the
        # same distribution. Because Hill numbers are non-increasing in the order,
        # this is <= the order-1 value for every matrix, with equality iff the
        # spectrum is flat. That is why R2/R3 read systematically lower than R1 and
        # why the two must never be compared.
        return float((total ** 2) / square_sum)
    raise ValueError(f"order must be 1 or 2, got {order!r}")


def effective_rank(x, *, centre: bool = CANONICAL.centre,
                   normalise_rows: bool = CANONICAL.normalise_rows,
                   order: int = CANONICAL.order,
                   tolerance: float | None = None,
                   variant: RankVariant | None = None) -> float:
    """Effective rank of a sample x feature matrix. Canonical = Roy & Vetterli, centred.

    Roy & Vetterli, "The effective rank: a measure of effective dimensionality",
    EUSIPCO 2007, Definition 1::

        erank(A) = exp(H(p_1, ..., p_Q)),  H = -sum p_k log p_k,  p_k = sigma_k / ||sigma||_1

    Defaults are ``CANONICAL`` and every call site should use them. The keyword
    ``tolerance`` overrides the relative singular-value cut, expressed as a fraction
    of the largest singular value. It exists for exactly one call site --
    `e0_basis_transfer._full_dictionary_rank`, which reports the effective rank
    beside a numerical rank taken at the FLOAT32 storage precision of its data and
    must use the same cut for both or the two disagree. Leave it None everywhere else.

    The other keyword
    arguments exist so that the three historical statistics of this repository can
    be measured against each other in one place, and so that a call site that must
    deviate (there is exactly one: the in-run rank tripwire in `training.py`, whose
    abort threshold was calibrated on ``R3``) has to say so in the call.

    Accepts a numpy array or a torch tensor. Torch inputs are computed on-device in
    float64 without a host transfer, because the training tripwire evaluates this
    every optimisation step; `test_effective_rank_canonical.py` asserts the two
    backends agree. Returns 0.0 when nothing survives the tolerance cut (e.g. a
    constant batch after centring).
    """
    if variant is not None:
        centre, normalise_rows, order = variant.centre, variant.normalise_rows, variant.order

    if hasattr(x, "detach") and hasattr(x, "device"):  # torch tensor: stay on device
        import torch
        value = x.detach().to(torch.float64)
        if value.ndim != 2:
            raise ValueError(f"effective_rank expects a 2-D matrix, got shape {tuple(value.shape)}")
        if normalise_rows:
            value = value / torch.linalg.vector_norm(value, dim=-1, keepdim=True).clamp_min(1e-300)
        scale = float(torch.linalg.matrix_norm(value))
        if centre:
            value = value - value.mean(dim=0, keepdim=True)
        singular = torch.linalg.svdvals(value)
        if singular.numel() == 0 or float(singular[0]) <= scale * DEGENERACY_TOLERANCE:
            return 0.0
        cut = max(value.shape) * EPS if tolerance is None else float(tolerance)
        singular = singular[singular > singular[0] * cut]
        if singular.numel() == 0:
            return 0.0
        return _finalise(singular, order, torch.exp, torch.log,
                         singular.sum(), singular.square().sum())

    value = np.asarray(x, dtype=np.float64)
    if value.ndim != 2:
        raise ValueError(f"effective_rank expects a 2-D matrix, got shape {value.shape}")
    if normalise_rows:
        norms = np.linalg.norm(value, axis=-1, keepdims=True)
        value = value / np.maximum(norms, 1e-300)
    scale = float(np.linalg.norm(value))
    if centre:
        value = value - value.mean(axis=0, keepdims=True)
    singular = np.linalg.svd(value, compute_uv=False)
    if singular.size == 0 or singular[0] <= scale * DEGENERACY_TOLERANCE:
        return 0.0
    cut = max(value.shape) * EPS if tolerance is None else float(tolerance)
    singular = singular[singular > singular[0] * cut]
    if singular.size == 0:
        return 0.0
    return _finalise(singular, order, np.exp, np.log,
                     singular.sum(), np.square(singular).sum())


def effective_rank_all_variants(x) -> dict[str, float]:
    """Every named variant of `RANK_VARIANTS` on one matrix, for side-by-side reporting."""
    return {name: effective_rank(x, variant=variant) for name, variant in RANK_VARIANTS.items()}


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


def _heldout_cca_fit(x: np.ndarray, y: np.ndarray, train: np.ndarray, n_components: int):
    """Shared fit for the held-out CCA family: whitening maps + the top singular pair.

    Returns ``(cx, wx, cy, wy, a, bt)`` or ``None`` when either side whitens to
    zero components. Split out so that :func:`heldout_cca_projection` and
    :func:`heldout_cca_directions` cannot drift: the projection is
    ``(x[test] - cx) @ wx @ a[:, 0]`` and the direction is ``wx @ a[:, 0]``, the
    same two factors associated differently. The projection keeps its original
    association order so that no published held-out number moves by a float ulp.
    """
    cx, wx = _whiten_map(x[train], n_components)
    cy, wy = _whiten_map(y[train], n_components)
    if wx.shape[1] == 0 or wy.shape[1] == 0:
        return None
    ux, uy = (x[train] - cx) @ wx, (y[train] - cy) @ wy
    a, _, bt = np.linalg.svd(ux.T @ uy, full_matrices=False)
    return cx, wx, cy, wy, a, bt


def heldout_cca_directions(x: np.ndarray, y: np.ndarray, train: np.ndarray, *,
                           n_components: int = 32) -> tuple[np.ndarray, np.ndarray]:
    """The top canonical direction PAIR fit on ``train`` rows, in x/y column space.

    Returns ``(u, v)`` with ``u`` of length ``x.shape[1]`` and ``v`` of length
    ``y.shape[1]``, each L2-normalised, such that ``x @ u`` and ``y @ v`` are the
    canonical scores up to an additive constant (correlation is shift-invariant,
    so dropping the whitening centres changes nothing downstream).

    Why this exists. The CALIBRA ``detection_floor`` is measured by planting a
    spike along a direction pair and reading it back on that same pair, and until
    now the pair was always **random** on at least the image side, while the
    channel it grades (:func:`heldout_top_cca`) uses a pair that was **fitted**.
    Handing this pair to ``calibration.spike_targets`` is what makes a
    direction-matched floor measurable at all. The pair is fitted on ``train``
    rows only, so it does not see the rows a held-out statistic is scored on.

    Sign is arbitrary (an SVD singular pair is defined up to a joint sign flip);
    the pair is jointly consistent, which is all a planted-correlation
    construction needs.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    train = np.asarray(train, dtype=np.int64)
    if len(train) < 10:
        return np.zeros(0), np.zeros(0)
    fit = _heldout_cca_fit(x, y, train, n_components)
    if fit is None:
        return np.zeros(0), np.zeros(0)
    _, wx, _, wy, a, bt = fit
    u, v = wx @ a[:, 0], wy @ bt[0]
    norm_u, norm_v = np.linalg.norm(u), np.linalg.norm(v)
    if norm_u < 1e-12 or norm_v < 1e-12:
        return np.zeros(0), np.zeros(0)
    return u / norm_u, v / norm_v


def heldout_cca_projection(x: np.ndarray, y: np.ndarray, train: np.ndarray, test: np.ndarray, *,
                           n_components: int = 32) -> tuple[np.ndarray, np.ndarray]:
    """Project ``test`` rows onto the top canonical pair FIT ON ``train`` rows.

    Returns ``(px, py)``, the two held-out score vectors whose correlation *is*
    :func:`heldout_top_cca_indexed`. Exposing the projection rather than only the
    scalar is what makes a permutation null and a bootstrap over the held-out rows
    affordable: the whitening and the SVD are the expensive part and they depend
    only on ``train``, so a null that permutes or resamples ``test`` reuses them
    unchanged. That is not merely an optimisation — recomputing the fit per
    permutation would refit the directions to each permuted sample and silently
    test a different hypothesis.

    Returns two empty arrays when either side whitens to zero components.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    train = np.asarray(train, dtype=np.int64)
    test = np.asarray(test, dtype=np.int64)
    if len(train) < 10 or len(test) < 10:
        return np.zeros(0), np.zeros(0)

    fit = _heldout_cca_fit(x, y, train, n_components)
    if fit is None:
        return np.zeros(0), np.zeros(0)
    cx, wx, cy, wy, a, bt = fit

    px = ((x[test] - cx) @ wx) @ a[:, 0]
    py = ((y[test] - cy) @ wy) @ bt[0]
    return px, py


def paired_absolute_correlation(px: np.ndarray, py: np.ndarray) -> float:
    """``|corr(px, py)|``, NaN when either side is constant. The channel's readout."""
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    if px.size < 2 or py.size != px.size or px.std() < 1e-12 or py.std() < 1e-12:
        return float("nan")
    return float(abs(np.corrcoef(px, py)[0, 1]))


def heldout_top_cca_indexed(x: np.ndarray, y: np.ndarray, train: np.ndarray, test: np.ndarray, *,
                            n_components: int = 32) -> float:
    """Top canonical correlation, directions fit on ``train`` rows, scored on ``test`` rows.

    The row-index form of :func:`heldout_top_cca`. It exists because the split that
    matters is not always random: a leave-sites-out design has to hand this function
    the rows belonging to held-out tissue source sites, and a caller that re-rolled
    the split itself would be a second implementation of the statistic.
    """
    return paired_absolute_correlation(*heldout_cca_projection(x, y, train, test,
                                                               n_components=n_components))


def heldout_top_cca(x: np.ndarray, y: np.ndarray, *, n_components: int = 32,
                    seed: int = 42, train_fraction: float = 0.5) -> float:
    """Top canonical correlation with directions FIT ON TRAIN, SCORED ON HELD-OUT rows.

    In-sample CCA is a multivariate maximum and is badly upward-biased at finite n
    — quoting it is like reporting a score on questions you already saw. Here the
    canonical directions are estimated on one split and the correlation is measured
    on the other, giving an unbiased (and much smaller) absolute number.

    The split is a uniformly random one. Everything else is
    :func:`heldout_top_cca_indexed`, which is the single implementation.
    """
    n = len(np.asarray(x))
    order = np.random.default_rng(seed).permutation(n)
    cut = int(n * train_fraction)
    return heldout_top_cca_indexed(x, y, order[:cut], order[cut:], n_components=n_components)


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
