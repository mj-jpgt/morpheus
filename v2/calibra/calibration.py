"""Spike-recovery calibration — the CALIBRA methods contribution.

The problem this solves. Every confound-adjusted cross-modal result in this
field is reported without a sensitivity statement: when a paper residualises
site/cohort/purity out of an H&E->molecular analysis and finds a small effect,
nobody can say whether the adjustment *also* destroyed real signal, so a null is
uninterpretable. Three independent adversarial reviews of this project killed
three separate theses on exactly that objection.

The fix, borrowed in spirit from ERCC spike-ins in genomics and simulation-based
calibration in Bayesian methodology, and (to our knowledge) never applied to
certify a cross-modal morphology-molecular claim: inject a synthetic signal of
*known* strength into the molecular side, push it through the **identical**
pipeline including residualisation, and measure what comes out. That yields

  * a recovery curve  r_hat(r_true),
  * an empirical detection floor (smallest r_true recovered in >=80% of draws),
  * an attenuation slope,

against which the observed real-data effect can finally be read. It converts
"we adjusted and found nothing" into "our floor is 0.031 and we measured 0.068".

Spike construction. For a random unit direction ``v`` in molecular space and
``u`` in image space, let ``s = standardise(X u)`` and ``a = standardise(Y v)``.
We *replace* the v-component of Y with a signal having exactly the requested
correlation with s::

    a_perp = standardise(a - corr(s, a) * s)
    a_new  = r_true * s + sqrt(1 - r_true^2) * a_perp

so ``corr(s, a_new) == r_true`` by construction, before any residualisation.
Structured spikes (``direction="structured"``) use a real programme loading
vector for ``v`` instead of a random one, because a random direction is the
favourable case and would flatter the instrument.

TARGETED READOUT (2026-07-30 fix -- read this before changing the readout back).
The first implementation scored recovery with ``top_canonical_correlation``, a
*maximum* over ``n_components`` directions per side. The spike, however, lives on
exactly one known direction pair ``(u, v)``. On real data the ambient top-CCA sits
near 0.97, so a spike of r_true<=0.2 was invisible against it: every detection
floor came back NaN, level-0 recovery read 0.97 instead of ~0, and at r_true=0.2
the measured value *fell* (the replacement destroys pre-existing structure along
``v`` faster than the weak spike restores it). The instrument's headline output was
therefore unusable on the only data that matters.

The readout is now **direction-matched**: after pushing the spiked targets through
the identical residualisation, we score ``|corr(X_res u, Y_spiked_res v)|`` -- the
planted axis itself. Level 0 then reads ~0 (a genuine null, since ``a_new`` is
constructed orthogonal to ``s``), the curve is monotone, and the floor is a real
detection threshold.

Scale warning. The floor is in *single-direction* correlation units. The headline
real-data number (adjusted / held-out top-CCA) is a *multivariate maximum* over 16
components and is inflated by capacity. **The two are not on the same scale and must
never be compared directly.** ``observed_matched_direction`` is provided as the
same-units comparator; ``observed`` remains the multivariate figure the real
analysis reports.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .residualise import cross_fitted_residuals
from .spectral import top_canonical_correlation

__all__ = ["SpikeRecoveryResult", "spike_targets", "spike_recovery_curve", "permutation_null"]

try:                                          # optional; the serial path is identical
    from joblib import Parallel, delayed
except ImportError:                           # pragma: no cover
    Parallel = None


def _map(fn, items, n_jobs: int):
    """Run ``fn(*item)`` over ``items``. Results are order-preserving and seed-fixed,
    so n_jobs changes wall-clock only -- never the numbers."""
    if n_jobs == 1 or Parallel is None or len(items) < 2:
        return [fn(*item) for item in items]
    return Parallel(n_jobs=n_jobs, prefer="processes")(delayed(fn)(*item) for item in items)


def permutation_null(x: np.ndarray, y: np.ndarray, design: np.ndarray, *, strata=None,
                     n_permutations: int = 100, n_components: int = 32, seed: int = 42,
                     n_jobs: int = 1) -> dict:
    """Top-CCA under destroyed cross-modal pairing — the chance level.

    Essential companion to the recovery curve. A top canonical correlation is a
    multivariate *maximum* over ``n_components`` directions per side, so at finite
    n it is substantially inflated by capacity alone: an observed 0.45 is
    meaningless until you know that chance gives, say, 0.42. Rows of ``y`` are
    permuted **within strata** (normally cancer type) so that cancer-level
    structure is preserved and only the patient-level pairing is destroyed —
    permuting globally would conflate the pairing with the cohort effect the
    residualisation is already meant to remove.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    rng = np.random.default_rng(seed)
    x_residual = cross_fitted_residuals(x, design, seed=seed)
    strata = np.zeros(len(x), dtype=int) if strata is None else np.asarray(strata)

    observed = top_canonical_correlation(x_residual, cross_fitted_residuals(y, design, seed=seed),
                                         n_components=n_components)
    # Draw every permutation order up front: identical numbers at any n_jobs.
    orders = []
    for _ in range(n_permutations):
        order = np.arange(len(y))
        for level in np.unique(strata):
            idx = np.flatnonzero(strata == level)
            order[idx] = rng.permutation(idx)
        orders.append((order,))

    def _one_permutation(order):
        permuted = cross_fitted_residuals(y[order], design, seed=seed)
        return top_canonical_correlation(x_residual, permuted, n_components=n_components)

    null = np.asarray(_map(_one_permutation, orders, n_jobs), dtype=np.float64)
    exceed = int(np.sum(null >= observed))
    return {
        "observed_top_cca": float(observed),
        "null_median": float(np.median(null)),
        "null_p95": float(np.percentile(null, 95)),
        "null_max": float(np.max(null)),
        "excess_over_null_median": float(observed - np.median(null)),
        "permutation_p": float((exceed + 1) / (n_permutations + 1)),
        "n_permutations": int(n_permutations),
    }


def _standardise(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=np.float64)
    centred = vector - vector.mean()
    scale = centred.std()
    return centred / scale if scale > 1e-12 else centred


def _correlation(a: np.ndarray, b: np.ndarray) -> float:
    """|Pearson r| between two 1-D scores, via standardisation (no covariance matrix)."""
    a = _standardise(a)
    b = _standardise(b)
    if a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float(abs(np.dot(a, b) / len(a)))


def spike_targets(x: np.ndarray, y: np.ndarray, r_true: float, *, rng: np.random.Generator,
                  molecular_direction: np.ndarray | None = None,
                  return_directions: bool = False):
    """Return a copy of ``y`` whose ``v``-component correlates with ``Xu`` at exactly ``r_true``.

    With ``return_directions=True`` returns ``(y_spiked, u, v)``. The caller needs
    ``(u, v)`` to score the spike on its own axis -- see the TARGETED READOUT note
    in the module docstring. The internal draw order (``u`` first, then ``v`` when
    no ``molecular_direction`` is supplied) is part of the contract: tests
    reconstruct the directions by re-seeding an identical generator.
    """
    if not 0.0 <= r_true < 1.0:
        raise ValueError("r_true must lie in [0, 1)")
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    u = rng.normal(size=x.shape[1])
    u /= np.linalg.norm(u)
    if molecular_direction is None:
        v = rng.normal(size=y.shape[1])
    else:
        v = np.asarray(molecular_direction, dtype=np.float64).copy()
    norm = np.linalg.norm(v)
    if norm < 1e-12:
        raise ValueError("molecular_direction must be non-zero")
    v /= norm

    s = _standardise(x @ u)
    a = _standardise(y @ v)
    rho = float(np.clip(np.dot(s, a) / len(s), -1.0, 1.0))
    a_perp = _standardise(a - rho * s)
    a_new = r_true * s + np.sqrt(max(0.0, 1.0 - r_true ** 2)) * a_perp
    # Replace, don't add: swap the v-component of Y for the constructed signal.
    # The update MUST be rescaled by the raw component's own s.d. ``a`` is
    # standardised while ``y @ v`` is not, so an unscaled ``a_new - a`` leaves a
    # residual (sigma - 1) * a in the target. That leftover carries the ambient
    # correlation rho, which then leaks into the readout: it put the level-0
    # baseline at 0.099 instead of ~0 and attenuated an r_true=0.6 spike to 0.27.
    # With the rescale, ``standardise(spiked @ v) == a_new`` exactly.
    sigma = float(np.std(y @ v))
    spiked = y + np.outer((a_new - a) * (sigma if sigma > 1e-12 else 1.0), v)
    return (spiked, u, v) if return_directions else spiked


@dataclass
class SpikeRecoveryResult:
    levels: np.ndarray
    recovered: np.ndarray            # (n_levels, n_draws)
    detection_floor: float
    attenuation_slope: float
    observed: float = float("nan")
    n_components: int = 32
    recovery_fraction: float = 0.8
    meta: dict = field(default_factory=dict)

    delta: np.ndarray | None = None   # per-draw increment over that draw's level-0 baseline

    def summary(self) -> dict:
        median = np.nanmedian(self.recovered, axis=1)
        delta = self.delta if self.delta is not None else self.recovered - self.recovered[0][None, :]
        matched = self.meta.get("observed_matched_direction", float("nan"))
        return {
            "levels": self.levels.tolist(),
            "recovered_median": median.tolist(),
            "delta_median": np.nanmedian(delta, axis=1).tolist(),
            "delta_p10": np.nanpercentile(delta, 10, axis=1).tolist(),
            "recovered_p10": np.nanpercentile(self.recovered, 10, axis=1).tolist(),
            "recovered_p90": np.nanpercentile(self.recovered, 90, axis=1).tolist(),
            "detection_floor": self.detection_floor,
            "attenuation_slope": self.attenuation_slope,
            "observed": self.observed,
            # SCALE-SAFE: compares like with like. ``observed`` is a multivariate
            # maximum; the floor is a single-direction correlation. Comparing those
            # two was the original error, so the flag uses the matched statistic.
            "observed_above_floor": bool(np.isfinite(matched) and np.isfinite(self.detection_floor)
                                         and matched > self.detection_floor),
            "floor_scale": "targeted_single_direction",
            "n_components": self.n_components,
            "recovery_fraction": self.recovery_fraction,
            **self.meta,
        }


def spike_recovery_curve(x: np.ndarray, y: np.ndarray, design: np.ndarray, *,
                         levels=(0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40),
                         n_draws: int = 25, n_components: int = 32, seed: int = 42,
                         recovery_fraction: float = 0.8,
                         molecular_directions: np.ndarray | None = None,
                         n_jobs: int = 1) -> SpikeRecoveryResult:
    """Push known-strength spikes through the full pipeline and measure recovery.

    ``design`` is the confound design matrix; residualisation is applied to the
    spiked data exactly as it is to the real data — that identity is the whole
    point, so do not "optimise" it away.

    The detection floor is the smallest non-zero level whose recovered statistic
    exceeds the null (level 0) 90th percentile in at least ``recovery_fraction``
    of draws. Using the level-0 upper tail as the reference is what makes the
    floor a *detection* threshold rather than a bias estimate.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    levels = np.asarray(levels, dtype=np.float64)
    rng = np.random.default_rng(seed)

    x_residual = cross_fitted_residuals(x, design, seed=seed)
    y_residual = cross_fitted_residuals(y, design, seed=seed)

    # Draw every (direction, seed) up front so the result is identical regardless of
    # how the draws are scheduled across workers.
    plan = []
    for _ in range(n_draws):
        direction = None
        if molecular_directions is not None and len(molecular_directions):
            direction = molecular_directions[rng.integers(len(molecular_directions))]
        plan.append((int(rng.integers(1 << 31)), direction))

    def _one_draw(draw_seed: int, direction):
        """PAIRED: the same (u, v) is reused across every level within a draw."""
        column = np.full(len(levels), np.nan)
        matched = float("nan")
        for i, level in enumerate(levels):
            spiked, u, v = spike_targets(x, y, float(level),
                                         rng=np.random.default_rng(draw_seed),
                                         molecular_direction=direction,
                                         return_directions=True)
            spiked_residual = cross_fitted_residuals(spiked, design, seed=seed)
            # TARGETED readout: score the planted axis, not a maximum over the whole
            # subspace. See the module docstring -- a max readout is swamped by
            # ambient structure and returns NaN floors on real data.
            column[i] = _correlation(x_residual @ u, spiked_residual @ v)
            if i == 0:
                # Same directions, same units, but on the UNSPIKED real targets.
                matched = _correlation(x_residual @ u, y_residual @ v)
        return column, matched

    results = _map(_one_draw, plan, n_jobs)
    recovered = np.column_stack([column for column, _ in results])
    observed_matched = float(np.nanmedian([m for _, m in results]))

    # Level 0 is now a genuine null (a_new is constructed orthogonal to s), so the
    # floor is read directly off the recovered values against the level-0 upper tail.
    null_reference = float(np.nanpercentile(recovered[0], 90)) if np.isfinite(recovered[0]).any() else 0.0
    detection_floor = float("nan")
    for i, level in enumerate(levels):
        if level <= 0:
            continue
        hits = np.nanmean(recovered[i] > max(null_reference, 1e-9))
        if np.isfinite(hits) and hits >= recovery_fraction:
            detection_floor = float(level)
            break

    median_recovered = np.nanmedian(recovered, axis=1)
    finite = np.isfinite(levels) & np.isfinite(median_recovered)
    slope = float("nan")
    if finite.sum() >= 2:
        # d(recovered)/d(r_true): 1.0 means the adjustment cost nothing, 0 means it
        # destroyed the signal entirely. This is the attenuation the reviewers ask for.
        slope = float(np.polyfit(levels[finite], median_recovered[finite], 1)[0])

    delta = recovered - recovered[0][None, :]
    observed = top_canonical_correlation(x_residual, y_residual, n_components=n_components)

    # REAL-DATA GATE. With a targeted readout, level 0 must sit near the sampling
    # null (~1/sqrt(n)), not near 1. A baseline near 1 means the readout is picking
    # up ambient structure instead of the spike -- exactly the defect this fix
    # repairs -- and every floor below it would be meaningless.
    baseline = float(np.nanmedian(recovered[0]))
    null_scale = 3.0 / np.sqrt(max(len(x), 2))
    return SpikeRecoveryResult(levels=levels, recovered=recovered, detection_floor=detection_floor,
                               attenuation_slope=slope, observed=float(observed),
                               n_components=n_components, recovery_fraction=recovery_fraction,
                               delta=delta,
                               meta={"null_reference_p90": null_reference, "n_draws": n_draws,
                                     "n_patients": int(len(x)),
                                     "baseline_recovered_median": baseline,
                                     "baseline_null_scale": float(null_scale),
                                     "baseline_is_null_like": bool(baseline < max(null_scale, 0.05)),
                                     "observed_matched_direction": observed_matched,
                                     "observed_multivariate_top_cca": float(observed),
                                     "baseline_top_cca": float(top_canonical_correlation(
                                         x_residual, y_residual, n_components=n_components))})
