"""Induced-correlation sweep — Track 2 / phase-gate item C.

WHAT THIS MEASURES
------------------
Residualising two *orthogonal* signals through a shared confound design induces a
correlation between them. CALIBRA's spike is built orthogonal to the image score
by construction (``a_perp`` in :mod:`calibration`), so a level-0 recovery curve
must read 0 -- and on real data it reads 0.067-0.140 for the 99-column
cancer+TSS design at n=2,530. That single observation is what this module turns
into a reported *function* of (design, design rank, n, residualiser knobs, seed).

The measured quantity is already computed and stored by
``spike_recovery_curve``: ``meta["confound_induced_baseline"]`` (magnitude),
``baseline_recovered_median_signed`` (signed), ``baseline_null_scale`` (the
sampling null) and ``baseline_gate_threshold``. This module supplies the driver
that was missing, plus the closed form below.

THE MECHANISM, DERIVED (this is the prediction T2.4 asks to be written down)
---------------------------------------------------------------------------
Let ``s = standardise(X u)`` and ``a`` the level-0 spike target, built so that
``corr(s, a) = 0`` exactly. Residualisation replaces each by its cross-fitted
residual, ``s_res = s - s_hat``, ``a_res = a - a_hat``, where the hats are
predictions from the shared design D. Because the SAME D produces both hats and
the fit is cross-fitted (so the hat carries no fold-local noise),

    cov(s_res, a_res) = cov(s, a) - cov(s, a_hat) - cov(s_hat, a) + cov(s_hat, a_hat)
                      = 0 - c - c + c = -c ,      c = cov(s_hat, a_hat)

i.e. **the induced covariance is exactly minus the covariance of the two
design-explained parts.** Standardising,

    r_induced = - R_s R_a cos(theta) / sqrt((1 - R_s^2)(1 - R_a^2))          (1)

with ``R_s^2, R_a^2`` the design's cross-fitted R^2 for the two scores and
``cos(theta) = corr(s_hat, a_hat)``. Define the *residualisation gain*

    kappa = R_s R_a / sqrt((1 - R_s^2)(1 - R_a^2)) .                          (2)

Equation (1) is an **identity**, not a model: ``closed_form_max_abs_error``
checks it against the pipeline draw for draw and it holds to ~1e-16. Everything
below is the separate, genuinely uncertain question of how its two factors scale.

``u`` and ``v`` are drawn at random, so ``s_hat`` and ``a_hat`` are vectors
inside the design span and ``cos(theta)`` has mean 0 and variance ``1 / k_eff``,
where ``k_eff`` is the *participation ratio* of the ridge-shrunk design spectrum
(``design_participation_ratio``), not the raw column count. That gives

    P1:   median |r_induced|  =  0.6745 * kappa / sqrt(k_eff) .               (3a)

P1 is the naive reading and it OVER-predicts, for a reason worth stating in
advance. The spike target ``a`` is built orthogonal to ``s`` *globally*, so its
design part is ``a_hat ~ (a0_hat - rho * s_hat)``: the construction removes from
``a_hat`` exactly the component along ``s_hat`` that (1) is about to charge for.
Carrying that through cancels the ``(1 - R_s^2)`` deflation in kappa and leaves

    P2:   median |r_induced|  =  0.6745 * R_s * R_a / sqrt(k_eff) .           (3b)

Both are declared BEFORE the sweep and both are reported with residuals,
alongside the plan's own guess

    P0:   |r_induced| ~ k / n .

PROVENANCE OF P3 -- STATED, NOT HIDDEN. P1 and P2 both fail on the real rank
ladder: the measured value is essentially FLAT from k=33 to k=501 while both
predict a 1/sqrt(k) decay. P3 below was written after seeing that failure on one
n (2,530, seed 42) and is therefore a *post hoc* form on the rank axis; it is
genuinely out of sample on every other n and seed in the grid, and it is
reported as such. The reason P1/P2 fail is that they assume the two
design-explained parts spread evenly over the design span. They do not: on a
real cohort almost everything the design explains about BOTH modalities lives on
cancer type, so 400 extra tissue-source-site columns raise the rank without
raising the dimension the two explained parts actually share. Replacing the
denominator with that shared dimension (``shared_design_geometry``) gives

    P3:   median |r_induced|  =  0.6745 * R_x * R_y / sqrt(k_eff_shared) ,     (3c)

where R_x, R_y and k_eff_shared come from the (D, X, Y) spectra alone -- no
draws, no spike, no recovery curve.

Two consequences that make (3b)/(3c) falsifiable, and that P0 gets wrong:

  * ``|r_induced|`` **decreases** with design rank at fixed R (as 1/sqrt(k_eff)),
    it does not grow with it. Widening a design makes each individual induced
    correlation smaller even as the design explains more;
  * ``|r_induced|`` has **no explicit n**. R_s, R_a and k_eff converge as n grows,
    so the induced correlation converges to a non-zero constant. It is a *bias*,
    not a sampling fluctuation, and more patients do not remove it. P0 says it
    vanishes as 1/n; that is the disagreement the n-sweep decides.

WHY THIS IS NOT AN ARTEFACT OF OUR DESIGN MATRIX
------------------------------------------------
Equation (1) says the effect is driven by kappa -- how much of *both* scores the
design explains -- and not by rank per se. The falsifier is therefore built in:
a design of matched rank and matched n that explains neither score (a Gaussian
design, or the real design with its rows permuted) must induce ~0. Both are run
as first-class cells of the sweep, with ``design_mode`` recording which.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .calibration import _correlation, _standardise, spike_recovery_curve
from .residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site

__all__ = [
    "DesignSpec", "stratified_subsample", "build_design", "design_participation_ratio",
    "shared_design_geometry",
    "planted_directions", "closed_form_induced_correlation", "predict_induced_correlation",
    "sweep_cell", "fit_scaling_law",
]

_MEDIAN_ABS_NORMAL = 0.6744897501960817          # median |N(0,1)|


# --------------------------------------------------------------------------- #
# cohort                                                                       #
# --------------------------------------------------------------------------- #

def stratified_subsample(strata, n_target: int, *, seed: int = 0) -> np.ndarray:
    """Indices of a size-``n_target`` subsample that preserves stratum shares.

    Composition is held fixed on purpose. An unstratified subsample changes the
    cancer mix as n moves, which would confound the n-sweep with a cohort-effect
    sweep -- and cohort effects are exactly what the design residualises, so the
    two would be inseparable. Allocation is proportional with largest-remainder
    rounding, so the realised n is exactly ``n_target`` whenever that is feasible.
    """
    strata = np.asarray(strata)
    total = len(strata)
    if n_target >= total:
        return np.arange(total)
    if n_target < 1:
        raise ValueError("n_target must be >= 1")
    rng = np.random.default_rng(seed)
    levels, counts = np.unique(strata, return_counts=True)
    exact = counts * (n_target / total)
    take = np.floor(exact).astype(int)
    remainder = n_target - int(take.sum())
    if remainder:
        order = np.argsort(-(exact - take), kind="stable")
        take[order[:remainder]] += 1
    take = np.minimum(take, counts)
    chosen: list[np.ndarray] = []
    for level, want in zip(levels, take):
        idx = np.flatnonzero(strata == level)
        if want:
            chosen.append(rng.choice(idx, size=int(want), replace=False))
    picked = np.sort(np.concatenate(chosen)) if chosen else np.empty(0, dtype=int)
    # Largest-remainder can still fall short once per-stratum caps bite.
    if len(picked) < n_target:
        rest = np.setdiff1d(np.arange(total), picked, assume_unique=False)
        picked = np.sort(np.concatenate([picked, rng.choice(rest, size=n_target - len(picked), replace=False)]))
    return picked


# --------------------------------------------------------------------------- #
# designs                                                                      #
# --------------------------------------------------------------------------- #

class DesignSpec:
    """One point on the design-rank ladder.

    ``mode`` selects the falsifier arm:
      ``real``     -- the design as built from the patients' own covariates;
      ``permuted`` -- the same matrix with its ROWS permuted: identical rank,
                      identical column marginals, zero relationship to any
                      patient. Equation (1) predicts ~0 here;
      ``gaussian`` -- an i.i.d. normal matrix of the same width. Matched rank,
                      no structure at all.
    The two synthetic arms are the ledger's stated falsifier ("induced
    correlation ~= 0 at matched design rank and n on a second design") run
    deliberately rather than hoped for.
    """

    def __init__(self, name: str, columns: tuple[str, ...] = (), *, min_site_count: int = 10,
                 mode: str = "real", width: int | None = None, needs: tuple[str, ...] = ()):
        if mode not in {"real", "permuted", "gaussian"}:
            raise ValueError(f"unknown design mode {mode!r}")
        if mode == "gaussian" and not width:
            raise ValueError("a gaussian design needs an explicit width")
        self.name, self.columns, self.min_site_count, self.mode, self.width = (
            name, tuple(columns), int(min_site_count), mode, width)
        self.needs = tuple(needs)          # extra covariates that must be supplied

    def __repr__(self) -> str:            # pragma: no cover - debugging aid
        return f"DesignSpec({self.name!r}, {self.columns}, min_site_count={self.min_site_count}, mode={self.mode!r})"


def build_design(spec: DesignSpec, patient_ids: np.ndarray, cancers: np.ndarray, *,
                 seed: int = 0, extra: dict | None = None) -> tuple[np.ndarray, dict]:
    """Materialise ``spec`` for one cohort. Never assume the column count -- read it."""
    import pandas as pd

    n = len(patient_ids)
    if spec.mode == "gaussian":
        design = np.random.default_rng(seed + 977).normal(size=(n, int(spec.width)))
        return design, {"n_confound_columns": int(spec.width), "n_distinct_sites_kept": 0}

    tss, frequent = pooled_tissue_source_site(patient_ids, min_site_count=spec.min_site_count)
    frame = pd.DataFrame({"cancer": np.asarray(cancers).astype(str), "tss": tss})
    for key, values in (extra or {}).items():
        frame[key] = values
    design = confound_design(frame, list(spec.columns))
    if spec.mode == "permuted":
        # Rows shuffled: rank and marginals are preserved exactly, the link to the
        # patient is destroyed. This is the matched-rank falsifier.
        design = design[np.random.default_rng(seed + 4231).permutation(n)]
    return design, {"n_confound_columns": int(design.shape[1]), "n_distinct_sites_kept": int(len(frequent))}


def shared_design_geometry(design: np.ndarray, x: np.ndarray, y: np.ndarray, *,
                           alpha: float = 1.0) -> dict:
    """``k_eff_shared``, ``R_x`` and ``R_y`` from the (D, X, Y) spectra alone.

    P1 and P2 use the design's own participation ratio as the denominator, which
    silently assumes the two design-explained parts spread evenly over the whole
    design span. On real cohorts they do not: cancer type carries most of what
    the design explains about *both* modalities, and 400 extra tissue-source-site
    columns add rank that neither score loads on. The denominator that matters is
    therefore the effective dimension of the **jointly explained** subspace.

    In the design's eigenbasis ``U`` (from the SVD of the centred design) with
    ridge shrinkage ``f_i = lam_i / (lam_i + alpha)``, write

        A = diag(f) U' X_c   (k x p) ,      B = diag(f) U' Y_c   (k x q)

    so the design-explained parts of the two scores are ``A u`` and ``B v``. For
    ``u``, ``v`` uniform on their spheres this gives, exactly,

        E[cos^2(theta)] = ||A' B||_F^2 / (||A||_F^2 ||B||_F^2)
        k_eff_shared    = 1 / E[cos^2(theta)]                                (4)
        R_x^2 = ||A||_F^2 / ||X_c||_F^2 ,   R_y^2 = ||B||_F^2 / ||Y_c||_F^2

    The cross term ``A'B`` is what makes (4) sensitive to *alignment* and not
    merely to how much variance each modality puts on the design. A design that
    explains both modalities strongly but on DIFFERENT directions has ``A'B ~ 0``,
    hence a huge ``k_eff_shared`` and a predicted induced correlation near zero --
    which is the correct answer and which a variance-only denominator gets wrong.

    Nothing here touches the spike, the recovery curve, or a single random draw:
    it is computable before the instrument is run, which is what makes P3 a
    prediction rather than a summary.
    """
    design = np.asarray(design, dtype=np.float64)
    out = {"k_eff_shared": 0.0, "predicted_r2_x": 0.0, "predicted_r2_y": 0.0}
    if design.size == 0 or design.shape[1] == 0:
        return out
    centred = design - design.mean(axis=0, keepdims=True)
    u_basis, singular, _ = np.linalg.svd(centred, full_matrices=False)
    keep = singular > 1e-8 * max(singular.max(), 1e-30)
    if not keep.any():
        return out
    u_basis, lam = u_basis[:, keep], singular[keep] ** 2
    f = (lam / (lam + alpha))[:, None]
    xc = np.asarray(x, dtype=np.float64) - np.mean(x, axis=0, keepdims=True)
    yc = np.asarray(y, dtype=np.float64) - np.mean(y, axis=0, keepdims=True)
    a_block, b_block = f * (u_basis.T @ xc), f * (u_basis.T @ yc)
    norm_a, norm_b = float(np.sum(a_block ** 2)), float(np.sum(b_block ** 2))
    total_x, total_y = float(np.sum(xc ** 2)), float(np.sum(yc ** 2))
    if norm_a > 0 and norm_b > 0:
        cross = float(np.sum((a_block.T @ b_block) ** 2))
        out["k_eff_shared"] = float(norm_a * norm_b / cross) if cross > 0 else float("inf")
    out["predicted_r2_x"] = float(norm_a / total_x) if total_x > 0 else 0.0
    out["predicted_r2_y"] = float(norm_b / total_y) if total_y > 0 else 0.0
    return out


def design_participation_ratio(design: np.ndarray, *, alpha: float = 1.0) -> dict:
    """``k_eff`` from the design ALONE -- no reference to X or Y.

    The relevant dimension is not the column count. One-hot blocks are rank
    deficient (a dropped-level collinearity), rare levels contribute almost no
    variance, and ridge shrinks small directions away. The projection weight in
    eigen-direction i is ``f_i = lam_i / (lam_i + alpha)``; the effective number
    of directions a random vector can spread over is the participation ratio
    ``(sum f)^2 / sum f^2``. This is the ``k_eff`` of equation (3) and it is a
    prediction, not a fitted quantity.
    """
    design = np.asarray(design, dtype=np.float64)
    if design.size == 0 or design.shape[1] == 0:
        return {"design_rank": 0, "k_eff": 0.0, "trace_hat": 0.0}
    centred = design - design.mean(axis=0, keepdims=True)
    lam = np.linalg.svd(centred, compute_uv=False) ** 2
    lam = lam[lam > 1e-9 * max(lam.max(), 1e-30)]
    f = lam / (lam + alpha)
    return {"design_rank": int(np.linalg.matrix_rank(centred)),
            "k_eff": float(f.sum() ** 2 / np.sum(f ** 2)) if f.size else 0.0,
            "trace_hat": float(f.sum())}


# --------------------------------------------------------------------------- #
# closed form                                                                  #
# --------------------------------------------------------------------------- #

def planted_directions(n_features_x: int, n_features_y: int, *, n_draws: int, seed: int):
    """Reproduce the exact (u, v) pairs ``spike_recovery_curve`` will plant.

    The draw order is part of ``spike_targets``' documented contract (u first,
    then v). Reconstructing it here means the closed form is compared to the
    pipeline **draw for draw**, not merely distribution to distribution -- which
    is the difference between "consistent with" and "is".
    """
    rng = np.random.default_rng(seed)
    draw_seeds = [int(rng.integers(1 << 31)) for _ in range(n_draws)]
    pairs = []
    for draw_seed in draw_seeds:
        gen = np.random.default_rng(draw_seed)
        u = gen.normal(size=n_features_x); u /= np.linalg.norm(u)
        v = gen.normal(size=n_features_y); v /= np.linalg.norm(v)
        pairs.append((u, v))
    return pairs


def closed_form_induced_correlation(x: np.ndarray, y: np.ndarray, design: np.ndarray, *,
                                    n_draws: int = 40, seed: int = 42, n_splits: int = 5,
                                    alpha: float = 1.0) -> dict:
    """Equation (1) evaluated on the same (u, v) the pipeline plants.

    Deliberately NOT a reimplementation of the pipeline: it never builds a spiked
    target matrix, never residualises the 256-dim or 90-target blocks, and works
    only on the two 1-D scores. It reuses ``cross_fitted_residuals`` unchanged, so
    if it agrees with the pipeline the agreement is about the mathematics rather
    than about shared code.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    rows = []
    for u, v in planted_directions(x.shape[1], y.shape[1], n_draws=n_draws, seed=seed):
        s = _standardise(x @ u)
        a_raw = _standardise(y @ v)
        rho = float(np.clip(np.dot(s, a_raw) / len(s), -1.0, 1.0))
        a = _standardise(a_raw - rho * s)                 # the level-0 spike target
        pair = np.column_stack([s, a])
        res = cross_fitted_residuals(pair, design, n_splits=n_splits, alpha=alpha, seed=seed)
        s_res, a_res = res[:, 0], res[:, 1]
        s_hat, a_hat = s - s_res - (s - s_res).mean(), a - a_res - (a - a_res).mean()
        # Cross-fitted R^2 is genuinely negative when the design predicts worse
        # than the mean -- that is the honest reading for a Gaussian design, and
        # it must map to a PREDICTED induced correlation of zero, not to NaN.
        r2_s = float(np.clip(1.0 - np.var(s_res) / max(np.var(s), 1e-30), 0.0, 0.999999))
        r2_a = float(np.clip(1.0 - np.var(a_res) / max(np.var(a), 1e-30), 0.0, 0.999999))
        gain = float(np.sqrt(r2_s * r2_a) / np.sqrt((1 - r2_s) * (1 - r2_a)))
        rows.append((_correlation(s_res, a_res), r2_s, r2_a, gain, _correlation(s_hat, a_hat)))
    r, r2_s, r2_a, gain, cos = (np.asarray(c, dtype=np.float64) for c in zip(*rows))
    finite_cos = cos[np.isfinite(cos)]
    return {
        "closed_form_induced_median": float(np.nanmedian(np.abs(r))),
        "closed_form_induced_median_signed": float(np.nanmedian(r)),
        "closed_form_induced_p10": float(np.nanpercentile(np.abs(r), 10)),
        "closed_form_induced_p90": float(np.nanpercentile(np.abs(r), 90)),
        "design_r2_x_median": float(np.nanmedian(r2_s)),
        "design_r2_y_median": float(np.nanmedian(r2_a)),
        "residualisation_gain_median": float(np.nanmedian(gain)),
        "fitted_part_cos_rms": float(np.sqrt(np.mean(finite_cos ** 2))) if finite_cos.size else float("nan"),
        "k_eff_empirical": float(1.0 / np.mean(finite_cos ** 2)) if finite_cos.size and np.mean(finite_cos ** 2) > 0 else float("nan"),
        "per_draw_closed_form": r.tolist(),
    }


def predict_induced_correlation(closed: dict, geometry: dict) -> dict:
    """Equations (3a) and (3b): numbers computed before the pipeline is ever run.

    Uses only the design's own participation ratio and the design's R^2 for the
    two scores. No free parameter is fitted anywhere in this function, and both
    candidate forms are returned so that neither can be selected after the fact.
    """
    k_eff = float(geometry.get("k_eff", 0.0))
    gain = float(closed.get("residualisation_gain_median", float("nan")))
    r2_x = max(float(closed.get("design_r2_x_median", float("nan"))), 0.0)
    r2_y = max(float(closed.get("design_r2_y_median", float("nan"))), 0.0)
    k_shared = float(geometry.get("k_eff_shared", 0.0))
    out = {"predicted_induced_correlation_p1": 0.0, "predicted_induced_correlation_p2": 0.0,
           "predicted_induced_correlation_p3": 0.0}
    if k_eff > 0:
        scale = _MEDIAN_ABS_NORMAL / np.sqrt(k_eff)
        out["predicted_induced_correlation_p1"] = float(scale * gain) if np.isfinite(gain) else float("nan")
        out["predicted_induced_correlation_p2"] = float(scale * np.sqrt(r2_x * r2_y))
    if k_shared > 0:
        # P3 uses the SHARED effective dimension and R^2 values derived from the
        # (D, X, Y) spectra -- no draws, no spike, no recovery curve.
        out["predicted_induced_correlation_p3"] = float(
            _MEDIAN_ABS_NORMAL * np.sqrt(max(geometry.get("predicted_r2_x", 0.0), 0.0)
                                         * max(geometry.get("predicted_r2_y", 0.0), 0.0)) / np.sqrt(k_shared))
    return out


# --------------------------------------------------------------------------- #
# one cell                                                                     #
# --------------------------------------------------------------------------- #

def sweep_cell(x: np.ndarray, y: np.ndarray, design: np.ndarray, *, levels=(0.0,),
               n_draws: int = 40, n_components: int = 16, seed: int = 42,
               n_splits: int = 5, alpha: float = 1.0, n_jobs: int = 1,
               with_closed_form: bool = True) -> dict:
    """Measure the induced correlation through the real instrument, then predict it.

    ``levels=(0.0,)`` is the cheap mode: the induced baseline is a level-0
    quantity, so the rest of the recovery curve is not needed for T2.1-T2.5.
    Passing the full level grid additionally yields both floors (T2.6) from the
    identical code path -- the floors must never come from a second
    implementation.
    """
    result = spike_recovery_curve(x, y, design, levels=tuple(levels), n_draws=n_draws,
                                  n_components=n_components, seed=seed, n_splits=n_splits,
                                  alpha=alpha, n_jobs=n_jobs)
    summary = result.summary()
    row = {
        "n_patients": int(len(x)),
        "n_draws": int(n_draws),
        "seed": int(seed),
        "residualiser_n_splits": int(n_splits),
        "residualiser_alpha": float(alpha),
        "induced_correlation_median": float(summary["confound_induced_baseline"]),
        "induced_correlation_median_signed": float(summary["baseline_recovered_median_signed"]),
        "induced_correlation_p10": float(np.nanpercentile(np.abs(result.recovered[0]), 10)),
        "induced_correlation_p90": float(np.nanpercentile(np.abs(result.recovered[0]), 90)),
        "baseline_null_scale": float(summary["baseline_null_scale"]),
        "baseline_gate_threshold": float(summary["baseline_gate_threshold"]),
        "baseline_is_null_like": bool(summary["baseline_is_null_like"]),
        "observed_multivariate_top_cca": float(summary["observed_multivariate_top_cca"]),
        "observed_matched_direction": float(summary["observed_matched_direction"]),
        "attenuation_slope": float(summary["attenuation_slope"]),
        "transmission_floor": float(summary["transmission_floor"]),
        "detection_floor": float(summary["detection_floor"]),
        "levels": list(summary["levels"]),
        "recovered_median": list(summary["recovered_median"]),
    }
    row.update(design_participation_ratio(design, alpha=alpha))
    row.update(shared_design_geometry(design, x, y, alpha=alpha))
    if with_closed_form:
        closed = closed_form_induced_correlation(x, y, design, n_draws=n_draws, seed=seed,
                                                 n_splits=n_splits, alpha=alpha)
        per_draw = np.asarray(closed.pop("per_draw_closed_form"), dtype=np.float64)
        row.update(closed)
        row.update(predict_induced_correlation(closed, row))
        measured = np.asarray(result.recovered[0], dtype=np.float64)
        both = np.isfinite(measured) & np.isfinite(per_draw)
        # PAIRED agreement: same draw, same (u, v). This is the mechanism check.
        row["closed_form_max_abs_error"] = float(np.max(np.abs(measured[both] - per_draw[both]))) if both.any() else float("nan")
        row["closed_form_pearson"] = (float(np.corrcoef(measured[both], per_draw[both])[0, 1])
                                      if both.sum() >= 3 and np.std(measured[both]) > 1e-12 else float("nan"))
    return row


# --------------------------------------------------------------------------- #
# scaling law                                                                  #
# --------------------------------------------------------------------------- #

def fit_scaling_law(rows, *, y_key: str = "induced_correlation_median") -> dict:
    """Fit ``log|r| = c + b1 log k_eff + b2 log n`` and report residuals.

    The plan's stated guess is ``|r| ~ k / n`` -- i.e. ``(b1, b2) = (+1, -1)``.
    Equation (3) predicts ``b1 = -1/2`` at fixed kappa and ``b2 = 0``. kappa is
    NOT held fixed across real designs (a wider design explains more), so the
    fitted ``b1`` is the sum of the -1/2 geometric term and however kappa moves;
    the parameter-free comparison in ``predicted_vs_measured`` is the decisive
    test and this fit is descriptive.
    """
    rows = [r for r in rows
            if np.isfinite(r.get(y_key, np.nan)) and r.get(y_key, 0) > 0
            and r.get("k_eff", 0) > 0 and r.get("n_patients", 0) > 0]
    if len(rows) < 4:
        return {"n_cells": len(rows), "status": "insufficient_cells"}
    log_k = np.log(np.asarray([r["k_eff"] for r in rows]))
    log_n = np.log(np.asarray([r["n_patients"] for r in rows], dtype=np.float64))
    log_r = np.log(np.asarray([r[y_key] for r in rows]))
    a = np.column_stack([np.ones_like(log_k), log_k, log_n])
    coef, *_ = np.linalg.lstsq(a, log_r, rcond=None)
    residual = log_r - a @ coef
    ss_tot = float(np.sum((log_r - log_r.mean()) ** 2))
    log_ks = np.log(np.maximum([r.get("k_eff_shared", 0.0) for r in rows], 1e-9))
    a_shared = np.column_stack([np.ones_like(log_ks), log_ks, log_n])
    coef_shared, *_ = np.linalg.lstsq(a_shared, log_r, rcond=None)
    residual_shared = log_r - a_shared @ coef_shared
    out = {
        "n_cells": len(rows),
        "intercept": float(coef[0]), "exponent_k_eff": float(coef[1]), "exponent_n": float(coef[2]),
        "exponent_k_eff_shared": float(coef_shared[1]), "exponent_n_given_k_shared": float(coef_shared[2]),
        "r2_with_k_eff_shared": float(1.0 - np.sum(residual_shared ** 2)
                                      / max(float(np.sum((log_r - log_r.mean()) ** 2)), 1e-30)),
        "r2": float(1.0 - np.sum(residual ** 2) / ss_tot) if ss_tot > 0 else float("nan"),
        "residual_rms_log": float(np.sqrt(np.mean(residual ** 2))),
        "residual_max_abs_log": float(np.max(np.abs(residual))),
        "plan_prediction_k_over_n": {"exponent_k_eff": 1.0, "exponent_n": -1.0},
        "derived_prediction_eq3": {"exponent_k_eff": -0.5, "exponent_n": 0.0},
    }
    measured = np.asarray([r[y_key] for r in rows])
    for key in ("predicted_induced_correlation_p1", "predicted_induced_correlation_p2",
                "predicted_induced_correlation_p3"):
        predicted = np.asarray([r.get(key, np.nan) for r in rows], dtype=np.float64)
        ok = np.isfinite(predicted) & (predicted > 0)
        if ok.sum() < 3:
            continue
        ratio = measured[ok] / predicted[ok]
        out[key.replace("predicted_induced_correlation_", "predicted_vs_measured_")] = {
            "n_cells": int(ok.sum()),
            "median_ratio": float(np.median(ratio)),
            "p10_ratio": float(np.percentile(ratio, 10)),
            "p90_ratio": float(np.percentile(ratio, 90)),
            # A predictor is only useful if its ERROR is small AND stable; the
            # spread matters more than the centre, because a constant offset can
            # be absorbed into the constant while a drifting one cannot.
            "log10_rms_error": float(np.sqrt(np.mean(np.log10(ratio) ** 2))),
            "log10_sd_error": float(np.std(np.log10(ratio))),
        }
    return out


# --------------------------------------------------------------------------- #
# driver                                                                       #
# --------------------------------------------------------------------------- #

def _load_cohort(artifact: Path, targets: Path, state: str, partition: str):
    raw = np.load(artifact, allow_pickle=True)
    tgt = np.load(targets, allow_pickle=True)
    target_ids = np.asarray([str(p) for p in tgt["patient_ids"]])
    names = np.asarray([str(t) for t in tgt["target_names"]])
    keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
    scores = np.asarray(tgt["scores"], dtype=np.float64)[:, keep]
    index = {pid: i for i, pid in enumerate(target_ids)}

    patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
    cancers = np.asarray([str(c) for c in raw["cancers"]])
    split = np.asarray([str(s) for s in raw["split"]])
    aligned = np.asarray([index.get(pid, -1) for pid in patient_ids])
    mask = (aligned >= 0)
    if partition != "all":
        mask &= (split == partition)
    x = np.asarray(raw[state], dtype=np.float64)[mask]
    return x, scores[aligned[mask]], patient_ids[mask], cancers[mask]


def _load_clinical_covariates(path: str, patient_ids: np.ndarray) -> dict[str, np.ndarray]:
    """Align ``dx_year``/``age`` to the cohort. Returns {} when the table is absent.

    An empty return is a real answer: the designs that need these covariates are
    then skipped and their absence appears in the sweep table as missing rows,
    which is the recorded-substitution behaviour the plan asks for.
    """
    if not path:
        return {}
    import pandas as pd

    table = pd.read_parquet(path)
    table["patient_id"] = table["patient_id"].astype(str).str.upper()
    lookup = table.drop_duplicates("patient_id").set_index("patient_id")
    ids = np.asarray([str(p).upper() for p in patient_ids])
    out: dict[str, np.ndarray] = {}
    for column in ("dx_year", "age"):
        if column in lookup.columns:
            out[column] = lookup[column].reindex(ids).to_numpy(dtype=np.float64)
    if "dx_year" in out:
        # Categorical form: same variable, ~35 levels instead of 1 numeric column.
        # Holding the VARIABLE fixed while changing only its encoding isolates
        # design rank from confound content -- the cleanest single rank knob here.
        year = out["dx_year"]
        finite = np.isfinite(year)
        labels = np.full(len(year), "NA", dtype=object)
        labels[finite] = np.round(year[finite]).astype(np.int64).astype(str)
        out["dx_year_cat"] = labels.astype(str)
    return out


def _default_designs() -> list[DesignSpec]:
    """The rank ladder.

    RECORDED SUBSTITUTION. The plan named ``["cancer","tss","purity"]`` and
    ``["cancer","tss","dx_year"]`` as the third and fourth rank points. **No TCGA
    purity table exists on either machine** (``~/e0_run/data`` holds none, and
    ``run_calibra --purity-table`` therefore has nothing to read) -- which is the
    same reason P3's D3 is NOT STARTED. Purity is dropped, not faked; an
    expression-derived surrogate was rejected because it is computed from the
    very RNA targets that form Y and would inflate ``R_a`` by construction, i.e.
    it would manufacture the effect under study.

    ``dx_year`` IS available (``year_of_initial_pathologic_diagnosis`` in the
    TCGA PanCan clinical mirror, 98.8% coverage) and is used, in both its numeric
    form (rank +2 including the missingness indicator) and its categorical form
    (rank +~35), together with ``age``. The TSS pooling threshold supplies the
    rest of the ladder: it moves design rank over more than an order of magnitude
    using only covariates already on disk, and -- unlike adding a variable -- it
    changes rank while holding the *kind* of confound fixed, which is what the
    scaling question actually needs.
    """
    return [
        DesignSpec("none", ()),
        DesignSpec("cancer", ("cancer",)),
        DesignSpec("tss_pool50", ("tss",), min_site_count=50),
        DesignSpec("tss_pool10", ("tss",), min_site_count=10),
        DesignSpec("cancer_tss_pool50", ("cancer", "tss"), min_site_count=50),
        DesignSpec("cancer_tss_pool10", ("cancer", "tss"), min_site_count=10),   # the anchor
        DesignSpec("cancer_tss_pool3", ("cancer", "tss"), min_site_count=3),
        DesignSpec("cancer_tss_pool1", ("cancer", "tss"), min_site_count=1),
        DesignSpec("cancer_tss_dxyear", ("cancer", "tss", "dx_year"), needs=("dx_year",)),
        DesignSpec("cancer_tss_dxyear_age", ("cancer", "tss", "dx_year", "age"), needs=("dx_year", "age")),
        DesignSpec("cancer_tss_dxyearcat", ("cancer", "tss", "dx_year_cat"), needs=("dx_year_cat",)),
        DesignSpec("permuted_cancer_tss_pool10", ("cancer", "tss"), min_site_count=10, mode="permuted"),
        DesignSpec("permuted_cancer_tss_pool1", ("cancer", "tss"), min_site_count=1, mode="permuted"),
        DesignSpec("gaussian_k99", mode="gaussian", width=99),
        DesignSpec("gaussian_k600", mode="gaussian", width=600),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--state", default="wsi_biology")
    parser.add_argument("--partition", default="all", choices=("test", "val", "all"))
    parser.add_argument("--n-grid", default="500,1000,2000,2530,4000,6427")
    parser.add_argument("--designs", default="", help="comma-separated subset of the default ladder")
    parser.add_argument("--levels", default="0.0", help="'0.0' for the cheap induced-only sweep")
    parser.add_argument("--n-draws", type=int, default=40)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--n-splits-grid", default="5")
    parser.add_argument("--alpha-grid", default="1.0")
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--tag", default="sweep")
    parser.add_argument("--clinical-covariates", default="",
                        help="parquet with patient_id + dx_year/age; enables the dx_year rank points")
    args = parser.parse_args()

    output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    x_all, y_all, ids_all, cancers_all = _load_cohort(Path(args.artifact), Path(args.targets),
                                                      args.state, args.partition)
    covariates_all = _load_clinical_covariates(args.clinical_covariates, ids_all)
    designs = _default_designs()
    if args.designs:
        wanted = [d.strip() for d in args.designs.split(",") if d.strip()]
        by_name = {d.name: d for d in designs}
        missing = [w for w in wanted if w not in by_name]
        if missing:
            raise SystemExit(f"unknown design(s): {missing}")
        designs = [by_name[w] for w in wanted]

    levels = tuple(float(v) for v in args.levels.split(","))
    n_grid = sorted({min(int(v), len(x_all)) for v in args.n_grid.split(",")})
    seeds = [int(s) for s in args.seeds.split(",")]
    splits_grid = [int(s) for s in args.n_splits_grid.split(",")]
    alpha_grid = [float(a) for a in args.alpha_grid.split(",")]

    rows: list[dict] = []
    for n_target in n_grid:
        for seed in seeds:
            idx = stratified_subsample(cancers_all, n_target, seed=seed)
            x, y = x_all[idx], y_all[idx]
            ids, cancers = ids_all[idx], cancers_all[idx]
            extra = {key: values[idx] for key, values in covariates_all.items()}
            for spec in designs:
                if any(key not in extra for key in spec.needs):
                    continue                    # covariate not on disk; absence is reported, not patched
                design, meta = build_design(spec, ids, cancers, seed=seed, extra=extra)
                for n_splits in splits_grid:
                    for alpha in alpha_grid:
                        row = {"design": spec.name, "design_mode": spec.mode,
                               "design_columns": "+".join(spec.columns) or "none",
                               "min_site_count": spec.min_site_count,
                               "state": args.state, "partition": args.partition,
                               "n_requested": int(n_target), **meta}
                        row.update(sweep_cell(x, y, design, levels=levels, n_draws=args.n_draws,
                                              n_components=args.n_components, seed=seed,
                                              n_splits=n_splits, alpha=alpha, n_jobs=args.n_jobs))
                        rows.append(row)
                        print(f"[{args.tag}] {spec.name:28s} mode={spec.mode:8s} n={row['n_patients']:5d} "
                              f"k={row['n_confound_columns']:4d} k_eff={row['k_eff']:7.2f} "
                              f"seed={seed} splits={n_splits} alpha={alpha:g} "
                              f"induced={row['induced_correlation_median']:.4f} "
                              f"p1={row.get('predicted_induced_correlation_p1', float('nan')):.4f} "
                              f"p2={row.get('predicted_induced_correlation_p2', float('nan')):.4f} "
                              f"p3={row.get('predicted_induced_correlation_p3', float('nan')):.4f} "
                              f"kshared={row.get('k_eff_shared', float('nan')):.2f}", flush=True)

    import pandas as pd
    frame = pd.DataFrame(rows)
    frame.to_csv(output / f"{args.tag}_rows.csv", index=False)
    (output / f"{args.tag}_law.json").write_text(json.dumps({
        "all_cells": fit_scaling_law(rows),
        "real_designs_only": fit_scaling_law([r for r in rows if r["design_mode"] == "real"]),
        "protocol": {"artifact": str(args.artifact), "state": args.state, "partition": args.partition,
                     "levels": list(levels), "n_draws": args.n_draws, "seeds": seeds,
                     "n_grid": n_grid, "n_splits_grid": splits_grid, "alpha_grid": alpha_grid,
                     "designs": [d.name for d in designs]},
    }, indent=2), encoding="utf-8")
    print(f"[done] {len(rows)} cells -> {output}", flush=True)


if __name__ == "__main__":
    main()
