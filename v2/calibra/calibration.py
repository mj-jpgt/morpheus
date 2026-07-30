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
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .residualise import cross_fitted_residuals
from .spectral import top_canonical_correlation

__all__ = ["SpikeRecoveryResult", "spike_targets", "spike_recovery_curve", "permutation_null"]


def permutation_null(x: np.ndarray, y: np.ndarray, design: np.ndarray, *, strata=None,
                     n_permutations: int = 100, n_components: int = 32, seed: int = 42) -> dict:
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
    null = np.empty(n_permutations)
    for i in range(n_permutations):
        order = np.arange(len(y))
        for level in np.unique(strata):
            idx = np.flatnonzero(strata == level)
            order[idx] = rng.permutation(idx)
        permuted = cross_fitted_residuals(y[order], design, seed=seed)
        null[i] = top_canonical_correlation(x_residual, permuted, n_components=n_components)
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


def spike_targets(x: np.ndarray, y: np.ndarray, r_true: float, *, rng: np.random.Generator,
                  molecular_direction: np.ndarray | None = None) -> np.ndarray:
    """Return a copy of ``y`` whose ``v``-component correlates with ``Xu`` at exactly ``r_true``."""
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
    return y + np.outer(a_new - a, v)


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
            "observed_above_floor": bool(np.isfinite(self.observed) and np.isfinite(self.detection_floor)
                                         and self.observed > self.detection_floor),
            "n_components": self.n_components,
            "recovery_fraction": self.recovery_fraction,
            **self.meta,
        }


def spike_recovery_curve(x: np.ndarray, y: np.ndarray, design: np.ndarray, *,
                         levels=(0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40),
                         n_draws: int = 25, n_components: int = 32, seed: int = 42,
                         recovery_fraction: float = 0.8,
                         molecular_directions: np.ndarray | None = None) -> SpikeRecoveryResult:
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
    recovered = np.full((len(levels), n_draws), np.nan)
    # PAIRED design: the same (u, v) is reused across every level within a draw, so
    # the level-0 value is that draw's own baseline. An unpaired design cannot detect
    # a small spike, because the absolute top-CCA is dominated by pre-existing
    # structure (and by CCA's own capacity) and swamps the injected effect.
    for draw in range(n_draws):
        direction = None
        if molecular_directions is not None and len(molecular_directions):
            direction = molecular_directions[rng.integers(len(molecular_directions))]
        draw_seed = int(rng.integers(1 << 31))
        for i, level in enumerate(levels):
            spiked = spike_targets(x, y, float(level), rng=np.random.default_rng(draw_seed),
                                   molecular_direction=direction)
            spiked_residual = cross_fitted_residuals(spiked, design, seed=seed)
            recovered[i, draw] = top_canonical_correlation(x_residual, spiked_residual,
                                                           n_components=n_components)

    # Work in per-draw INCREMENTS over that draw's own level-0 baseline.
    delta = recovered - recovered[0][None, :]
    null_reference = float(np.nanpercentile(delta[0], 90)) if np.isfinite(delta[0]).any() else 0.0
    detection_floor = float("nan")
    for i, level in enumerate(levels):
        if level <= 0:
            continue
        hits = np.nanmean(delta[i] > max(null_reference, 1e-9))
        if np.isfinite(hits) and hits >= recovery_fraction:
            detection_floor = float(level)
            break

    median_delta = np.nanmedian(delta, axis=1)
    finite = np.isfinite(levels) & np.isfinite(median_delta)
    slope = float("nan")
    if finite.sum() >= 2:
        slope = float(np.polyfit(levels[finite], median_delta[finite], 1)[0])

    observed = top_canonical_correlation(x_residual, cross_fitted_residuals(y, design, seed=seed),
                                         n_components=n_components)
    return SpikeRecoveryResult(levels=levels, recovered=recovered, detection_floor=detection_floor,
                               attenuation_slope=slope, observed=float(observed),
                               n_components=n_components, recovery_fraction=recovery_fraction,
                               delta=delta,
                               meta={"null_reference_p90": null_reference, "n_draws": n_draws,
                                     "n_patients": int(len(x)),
                                     "baseline_top_cca": float(np.nanmedian(recovered[0]))})
