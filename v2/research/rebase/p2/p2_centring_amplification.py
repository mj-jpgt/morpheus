"""Does removing a dominant shared direction amplify a spectral statistic's retraining floor?

PREDECLARED at `NOTEBOOK_ENTRIES/PREDECLARED_centring_amplification_law_20260804T1750Z.md`,
commit `d4e344c`, written and pushed before this module existed. The derivation,
the assumptions and the nine falsifiers (S1-S4 synthetic, P1-P5 real) are there
and are NOT restated here in a form that could drift from them.

WHAT IS BEING TESTED. Draft §4.1a records two facts that line up: RankMe's floor
on the exported `wsi_biology` raw block is 1.811x against our column-centred R1's
3.111x on the same five runs, and the same five runs' centred floor is 3.295x on
`wsi_biology` but 1.019x on `rna_biology`. The candidate law is that both are the
same phenomenon -- removing a dominant, re-run-stable shared component transfers
the statistic's reproducibility from the original spectrum to the residual one
and amplifies its estimation variance.

THE DERIVED FORM, in one line, so this module's variable names are readable::

    ln D_1(uncentred spectrum) = h(t) + (1 - t) * ln D_1(centred spectrum)

with ``h`` the binary entropy and ``t = s_dot / (s_dot + ||Z_c||_*)`` the L1
(nuclear-norm) share of the mean term, ``s_dot = sqrt(n) ||z_bar||``. Under the
assumption that ``t`` itself does not move across re-runs, the amplification of
the log-scale dispersion is ``A_1 = 1/(1-t)``, and for the general order-a Hill
number ``A_a = 1/(1-w_a)`` with ``w_a`` as in ``predicted_amplification`` below.

AN ALGEBRAIC POINT WORTH STATING, because it removes a free parameter. Writing
``f = n||z_bar||^2 / ||Z||_F^2`` for the energy share, ``||Z||_F^2 = s_dot^2 +
||Z_c||_F^2`` exactly, and ``R2_centred = (sum sigma)^2 / sum sigma^2 =
||Z_c||_*^2 / ||Z_c||_F^2``. Hence

    t = sqrt(f) / ( sqrt(f) + sqrt(R2_centred) * sqrt(1-f) )

is an IDENTITY, not an approximation: the ``f``-form and the ``t``-form are the
same number. The whole derivation therefore rests on exactly one approximation --
that the uncentred spectrum is ``{s_dot}`` union the centred spectrum -- and P1
measures it directly instead of assuming it.

NOTHING IS COMPUTED INLINE. Every statistic comes from `v2/calibra/spectral.py`
or from `p2_competing_metrics.py`. Four inline-formula substitutions have been
caught in this paper and the tell was identical each time: an arithmetic
expression where an import belonged. The quantities this module does compute
directly -- ``||z_bar||``, ``||Z_c||_F``, ``||Z_c||_*`` -- are vector and matrix
norms taken with `numpy.linalg.norm`, not rank statistics, and ``R2_centred`` is
read from `RANK_VARIANTS["R2"]` rather than from those norms so that the identity
above is a CHECK (`t` against `t_from_f`) and not a tautology in code.

THE UNCENTRED FORM OF A STATISTIC THAT HAS NO FLAG. Only `rankme` is uncentred;
`participation_ratio`, `stable_rank` and `alpha_req` centre internally and expose
no keyword, and `p2_competing_metrics.py` is vendored byte-for-byte and must not
be edited. So instead of reimplementing them, ``doubled(X) = [X; -X]`` has exact
zero column mean (centring is a no-op on it) and Gram ``2 X^T X``, so its spectrum
is ``sqrt(2)`` times ``X``'s uncentred spectrum. Every statistic here is
scale-invariant, so ``stat_centred(doubled(X)) == stat_uncentred(X)`` exactly.
``check_doubling_is_exact`` asserts this against the two statistics that DO have
an uncentred form of their own (`R1_uncentred`, `rankme`); the predeclaration
makes a failure a stopping condition.

usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 python3 p2_centring_amplification.py \
        --reps '/home/ubuntu/e0_run/d1_envelope/rep*.npz' \
        --targets /home/ubuntu/e0_run/data/frozen_rna_targets.npz \
        --output ~/ws_amp/out/P2_CENTRING_AMPLIFICATION.json

    ... --synthetic-only --output ~/ws_amp/out/P2_CENTRING_SYNTHETIC.json
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2.p2_competing_metrics import (alpha_req, participation_ratio,
                                                                 rankme, stable_rank)
from morpheus.v2.research.rebase.p2.p2_envelope_floors import ALPHA_INDEX_RANGE, UNTRAINED_40, VIEWS


# --------------------------------------------------------------------------
# the exact reduction from "uncentred statistic" to "centred statistic"
# --------------------------------------------------------------------------
def doubled(x: np.ndarray) -> np.ndarray:
    """``[X; -X]`` -- zero column mean, spectrum ``sqrt(2) * spec(X)``.

    Lets a centre-only statistic be evaluated uncentred without reimplementing it.
    """
    x = np.asarray(x, dtype=np.float64)
    return np.concatenate([x, -x], axis=0)


#: Statistic label -> (callable on one matrix, spectrum it lives on, Hill order).
#:
#: ``spectrum`` is "singular" for statistics taken on the L1-normalised SINGULAR
#: values (Roy-Vetterli's family) and "eigen" for those taken on the eigenvalues
#: of the covariance. The distinction decides which share the dominant component
#: contributes: ``t`` (nuclear-norm share) on the singular spectrum and ``f``
#: (energy share) on the eigen spectrum, because ``sum lambda = ||sigma||_2^2``
#: and the spike's own eigenvalue is ``s_dot^2``. Getting this wrong is the
#: difference between a prediction of 1.6x and one of 5x.
#:
#: ``order`` is ``None`` for a statistic that is not a Hill number of either
#: spectrum; it is reported but no amplification is predicted for it.
STATISTICS = {
    "R1": (lambda x: effective_rank(x, variant=RANK_VARIANTS["R1"]), "singular", 1),
    "R2": (lambda x: effective_rank(x, variant=RANK_VARIANTS["R2"]), "singular", 2),
    "R3": (lambda x: effective_rank(x, variant=RANK_VARIANTS["R3"]), "singular", None),
    "PR": (participation_ratio, "eigen", 2),
    "stable_rank": (stable_rank, "eigen", math.inf),
    "hard_rank": (lambda x: float(np.linalg.matrix_rank(np.asarray(x, dtype=np.float64))),
                  "singular", 0),
    "alpha_req_abs_dev": (lambda x: abs(alpha_req(x, *ALPHA_INDEX_RANGE)["alpha"] - 1.0),
                          "eigen", None),
}


def binary_entropy(t: float) -> float:
    """``h(t) = -t ln t - (1-t) ln(1-t)``, with the ``0 log 0 = 0`` convention."""
    if t <= 0.0 or t >= 1.0:
        return 0.0
    return float(-t * math.log(t) - (1.0 - t) * math.log1p(-t))


def predicted_amplification(order, share: float, centred_value: float) -> float | None:
    """``A_a = 1/(1-w_a)`` -- the derived amplification, Theorem 2 of the predeclaration.

    ``share`` is the dominant component's share of the relevant spectrum's L1 mass
    (``t`` for the singular spectrum, ``f`` for the eigen spectrum) and
    ``centred_value`` is ``D_a`` of the residual spectrum.

    Order 1 gives ``w_1 = share`` exactly. Order ``inf`` is the max-share case and
    is BINARY, not continuous: once the spike's own share exceeds the residual's
    largest, the uncentred statistic is ``1/share`` and stops responding to the
    residual at all, so the amplification is infinite (returned as ``None`` with
    the branch recorded by the caller, because an infinity in a JSON is a trap).
    """
    if order is None or not np.isfinite(centred_value) or centred_value <= 0:
        return None
    if not (0.0 < share < 1.0):
        return None
    if order == math.inf:
        # D_inf = 1/max_k p_k. The mixture's largest atom is max(share, (1-share)/D_inf_centred).
        return None if share * centred_value > (1.0 - share) else 1.0
    if order == 1:
        w = share
    else:
        spike = share ** order
        bulk = ((1.0 - share) ** order) * (centred_value ** (1.0 - order))
        if spike + bulk <= 0:
            return None
        w = spike / (spike + bulk)
    return float("inf") if w >= 1.0 else float(1.0 / (1.0 - w))


# --------------------------------------------------------------------------
# per-matrix readout
# --------------------------------------------------------------------------
def score_matrix(x: np.ndarray) -> dict:
    """Every statistic centred and uncentred on one block, plus f, t and the P1 residual."""
    x = np.asarray(x, dtype=np.float64)
    n = x.shape[0]
    mean = x.mean(axis=0)
    centred = x - mean

    s_dot = float(math.sqrt(n) * np.linalg.norm(mean))
    frob_total = float(np.linalg.norm(x))
    frob_centred = float(np.linalg.norm(centred))
    nuclear_centred = float(np.linalg.norm(centred, ord="nuc"))

    f = float(s_dot ** 2 / frob_total ** 2) if frob_total > 0 else float("nan")
    t = float(s_dot / (s_dot + nuclear_centred)) if (s_dot + nuclear_centred) > 0 else float("nan")

    # the mean patient-to-patient cosine, the one-line quantity the observation quotes
    rows = x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-300)
    gram = rows @ rows.T
    cos_bar = float((gram.sum() - np.trace(gram)) / (n * (n - 1)))

    entry = {"n": int(n), "d": int(x.shape[1]), "f": f, "t": t,
             "mean_offdiag_cosine": cos_bar, "s_dot": s_dot,
             "frobenius": frob_total, "frobenius_centred": frob_centred,
             "nuclear_centred": nuclear_centred, "statistics": {}}

    x_doubled = doubled(x)
    for label, (fn, spectrum, order) in STATISTICS.items():
        centred_value = float(fn(x))
        uncentred_value = float(fn(x_doubled))
        share = t if spectrum == "singular" else f
        pred = predicted_amplification(order, share, centred_value)
        entry["statistics"][label] = {
            "centred": centred_value, "uncentred": uncentred_value,
            "spectrum": spectrum, "order": None if order is None else (
                "inf" if order == math.inf else order),
            "share_used": share,
            "predicted_amplification": pred,
            "predicted_amplification_is_infinite": (order == math.inf and pred is None),
        }

    # RankMe as published, on the raw matrix -- the uncentred order-1 comparator §4.1a quotes
    entry["statistics"]["RankMe"] = {"centred": None, "uncentred": float(rankme(x)),
                                     "spectrum": "singular", "order": 1, "share_used": t,
                                     "predicted_amplification": None,
                                     "predicted_amplification_is_infinite": False}

    # --- P1: is the uncentred spectrum {s_dot} union the centred spectrum? ---
    d1_centred = entry["statistics"]["R1"]["centred"]
    d1_uncentred = entry["statistics"]["R1"]["uncentred"]
    predicted_log = binary_entropy(t) + (1.0 - t) * math.log(d1_centred) if d1_centred > 0 else None
    entry["identity"] = {
        "ln_D1_uncentred_observed": math.log(d1_uncentred) if d1_uncentred > 0 else None,
        "ln_D1_uncentred_predicted": predicted_log,
        "relative_error": (abs(math.log(d1_uncentred) - predicted_log) / abs(math.log(d1_uncentred))
                           if (predicted_log is not None and d1_uncentred > 1.0) else None),
        "_": "Theorem 1 of the predeclaration. Exact iff the mean direction is orthogonal to the "
             "centred row space; the ONLY approximation in the derivation.",
    }
    # --- the f-form against the t-form: an identity, so a mismatch is a code bug ---
    r2_centred = entry["statistics"]["R2"]["centred"]
    if 0.0 < f < 1.0 and r2_centred > 0:
        t_from_f = math.sqrt(f) / (math.sqrt(f) + math.sqrt(r2_centred) * math.sqrt(1.0 - f))
        entry["t_from_f"] = float(t_from_f)
        entry["t_identity_absolute_error"] = float(abs(t_from_f - t))
    else:
        entry["t_from_f"] = None
        entry["t_identity_absolute_error"] = None
    return entry


def check_doubling_is_exact(x: np.ndarray, tolerance: float = 1e-9) -> dict:
    """The stopping condition of §2 of the predeclaration, on the two statistics that can check it."""
    x = np.asarray(x, dtype=np.float64)
    xd = doubled(x)
    own = effective_rank(x, variant=RANK_VARIANTS["R1_uncentred"])
    via = effective_rank(xd, variant=RANK_VARIANTS["R1"])
    own2 = effective_rank(x, variant=RANK_VARIANTS["R2_uncentred"])
    via2 = effective_rank(xd, variant=RANK_VARIANTS["R2"])
    own_rankme, via_rankme = rankme(x), rankme(xd)
    out = {
        "R1_uncentred_direct": float(own), "R1_uncentred_via_doubling": float(via),
        "R1_relative_error": float(abs(via - own) / abs(own)),
        "R2_uncentred_direct": float(own2), "R2_uncentred_via_doubling": float(via2),
        "R2_relative_error": float(abs(via2 - own2) / abs(own2)),
        "rankme_direct": float(own_rankme), "rankme_via_doubling": float(via_rankme),
        "rankme_relative_error": float(abs(via_rankme - own_rankme) / abs(own_rankme)),
        "tolerance": tolerance,
    }
    out["passed"] = bool(max(out["R1_relative_error"], out["R2_relative_error"],
                             out["rankme_relative_error"]) <= tolerance)
    return out


# --------------------------------------------------------------------------
# dispersion across repeats
# --------------------------------------------------------------------------
def dispersion(values: list[float]) -> dict:
    """Both dispersions the predeclaration names: sd(ln .) and ln(max/min)."""
    v = np.asarray([float(x) for x in values], dtype=np.float64)
    if not np.isfinite(v).all() or (v <= 0).any():
        return {"n": int(v.size), "sd_log": None, "log_fold": None, "fold": None,
                "why": "a non-positive or non-finite value -- a log dispersion is not defined"}
    logs = np.log(v)
    return {"n": int(v.size), "sd_log": float(logs.std(ddof=1)),
            "log_fold": float(logs.max() - logs.min()),
            "fold": float(v.max() / v.min()),
            "min": float(v.min()), "max": float(v.max())}


def variance_decomposition(t_values, d1_centred, d1_uncentred) -> dict:
    """P3b: how much of the UNCENTRED statistic's re-run variance is the spike's own wobble?

    Theorem 2 says ``delta ln D_1(s) = (1-t) delta ln D_1(sigma) + c delta t`` with
    ``c = ln((1-t)/t) - ln D_1(sigma)``. A2 -- the "stable" in "dominant stable
    component" -- is the claim that the second term is negligible. If it is not,
    the law's premise fails on this data and its prediction is not being tested.
    """
    t = np.asarray(t_values, dtype=np.float64)
    lc = np.log(np.asarray(d1_centred, dtype=np.float64))
    lu = np.log(np.asarray(d1_uncentred, dtype=np.float64))
    if t.size < 3 or not (np.isfinite(t).all() and np.isfinite(lc).all() and np.isfinite(lu).all()):
        return {"defined": False}
    t_bar, lc_bar = float(t.mean()), float(lc.mean())
    if not (0.0 < t_bar < 1.0):
        return {"defined": False}
    c = math.log((1.0 - t_bar) / t_bar) - lc_bar
    var_total = float(lu.var(ddof=1))
    var_bulk = float(((1.0 - t_bar) ** 2) * lc.var(ddof=1))
    var_spike = float((c ** 2) * t.var(ddof=1))
    return {"defined": True, "t_mean": t_bar, "c": float(c),
            "sd_t": float(t.std(ddof=1)), "sd_ln_D1_centred": float(lc.std(ddof=1)),
            "sd_ln_D1_uncentred": float(lu.std(ddof=1)),
            "var_uncentred_observed": var_total,
            "var_from_bulk_term": var_bulk, "var_from_spike_term": var_spike,
            "spike_share_of_uncentred_variance": (float(var_spike / var_total)
                                                  if var_total > 0 else None),
            "var_predicted_sum_ignoring_covariance": var_bulk + var_spike}


def amplification(centred: list[float], uncentred: list[float]) -> dict:
    """Observed A = dispersion(centred)/dispersion(uncentred), on both dispersion measures."""
    a, b = dispersion(centred), dispersion(uncentred)
    out = {"centred": a, "uncentred": b}
    for key in ("sd_log", "log_fold"):
        num, den = a.get(key), b.get(key)
        out[f"A_from_{key}"] = (float(num / den) if (num is not None and den not in (None, 0.0))
                                else None)
    return out


# --------------------------------------------------------------------------
# synthetic sweep -- S1..S4
# --------------------------------------------------------------------------
#: The three perturbation conditions. Only the first instantiates the
#: predeclaration's A2; the other two are what "the mean is stable" actually buys.
CONDITIONS = ("t_pinned", "stable_mean", "unstable_mean")


def synthetic_sweep(*, n: int = 2766, d: int = 256, runs: int = 25,
                    f_grid=None, noise: float = 0.05, gain_noise: float = 0.15,
                    seed: int = 42, condition: str = "stable_mean") -> list[dict]:
    """``Z_r = a 1 mu^T + W_r`` at controlled energy share ``f``.

    ``W_r`` is an anisotropic bulk whose *directional gains* are re-drawn per run,
    ``decay_k * exp(gain_noise * g_rk)``, on top of a fixed set of latent scores,
    plus an additive perturbation at relative size ``noise``.

    WHY THE GAINS MOVE AND NOT JUST AN ADDITIVE NOISE. The first version of this
    sweep perturbed a *fixed* bulk additively. At n = 2766 that self-averages: the
    centred statistic's dispersion came out at ``sd(ln D_1) ~ 1e-4``, five orders
    below the 1.02x-3.3x folds the real repeats show, so every amplification ratio
    was a ratio of two numbers that were not moving. A re-run is not a fixed
    encoder plus noise; it is a different encoder. ``gain_noise`` is what makes the
    residual spectrum's *shape* move, and it is set so that the centred statistic's
    fold lands in the range the five real repeats occupy. The dispersion of the
    centred statistic is reported per row so this calibration is checkable rather
    than asserted.

    THREE CONDITIONS, AND WHY THERE ARE THREE RATHER THAN THE TWO PREDECLARED.
    The predeclaration called its stable-mean condition "A2 holds by
    construction". **That is wrong, and the smoke test found it before the sweep
    was run.** A2 is a statement about ``t = s_dot/(s_dot + ||Z_c||_*)``, and
    holding the mean term fixed pins only the numerator: the residual's own
    nuclear norm still moves from run to run, so ``t`` still moves. Holding the
    shared direction stable is **not sufficient** for A2, which is a defect in the
    candidate law's own statement of itself.

      * ``t_pinned``      -- each run's centred part is rescaled so that
                             ``||Z_c||_*`` matches run 1. Every statistic here is
                             scale-invariant on the centred matrix, so this changes
                             no centred value and pins ``t`` EXACTLY. This is the
                             only condition in which A2 holds, and it is a
                             POST-HOC ADDITION, added after the predeclared
                             stable-mean condition failed, and reported as such.
      * ``stable_mean``   -- the predeclared condition: ``a`` identical in every
                             run, and ``E_r``'s own column mean removed so that
                             the column mean is literally constant.
      * ``unstable_mean`` -- the predeclared A2-broken condition: ``E_r``
                             uncentred AND ``a_r = a (1 + eta_r)``.
    """
    if condition not in CONDITIONS:
        raise ValueError(f"condition must be one of {CONDITIONS}, got {condition!r}")
    unstable_mean = condition == "unstable_mean"
    if f_grid is None:
        f_grid = [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.65, 0.70,
                  0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.94, 0.95, 0.97, 0.99]
    rng = np.random.default_rng(seed)
    mu = rng.normal(size=d)
    mu /= np.linalg.norm(mu)
    decay = np.exp(-np.arange(d) / 40.0)                 # anisotropic bulk
    basis = np.linalg.qr(rng.normal(size=(d, d)))[0]
    scores = rng.normal(size=(n, d))                      # fixed latent scores

    out = []
    for f in f_grid:
        # energy share f  <=>  ||a mu term||_F^2 / total = f, with the bulk at unit Frobenius
        amplitude = math.sqrt(f / (1.0 - f)) / math.sqrt(n) if f < 1.0 else float("nan")
        per_run = {"centred": {k: [] for k in STATISTICS},
                   "uncentred": {k: [] for k in STATISTICS}}
        shares = {"f": [], "t": []}
        pinned_nuclear = None
        for r in range(runs):
            bulk = (scores * (decay * np.exp(gain_noise * rng.normal(size=d)))) @ basis.T
            bulk /= np.linalg.norm(bulk)                          # unit Frobenius
            noise_r = rng.normal(size=(n, d))
            if not unstable_mean:
                bulk -= bulk.mean(axis=0, keepdims=True)
                noise_r -= noise_r.mean(axis=0, keepdims=True)   # column mean literally constant
            noise_r *= noise / np.linalg.norm(noise_r)
            scale = amplitude * (1.0 + noise * rng.normal()) if unstable_mean else amplitude
            z = scale * np.ones((n, 1)) @ mu[None, :] + bulk + noise_r
            if condition == "t_pinned":
                centred_part = z - z.mean(axis=0, keepdims=True)
                nuc_r = float(np.linalg.norm(centred_part, ord="nuc"))
                if pinned_nuclear is None:
                    pinned_nuclear = nuc_r
                # scale-invariance of every centred statistic means this changes no
                # centred value; it pins t exactly, which is what A2 asserts
                z = z.mean(axis=0, keepdims=True) + centred_part * (pinned_nuclear / nuc_r)
            zd = doubled(z)
            for label, (fn, _, _) in STATISTICS.items():
                per_run["centred"][label].append(float(fn(z)))
                per_run["uncentred"][label].append(float(fn(zd)))
            mean = z.mean(axis=0)
            s_dot = math.sqrt(n) * float(np.linalg.norm(mean))
            nuc = float(np.linalg.norm(z - mean, ord="nuc"))
            shares["f"].append(float(s_dot ** 2 / float(np.linalg.norm(z)) ** 2))
            shares["t"].append(float(s_dot / (s_dot + nuc)))
        row = {"f_target": f, "f_measured": float(np.mean(shares["f"])),
               "t_measured": float(np.mean(shares["t"])),
               "t_sd": float(np.std(shares["t"], ddof=1)), "runs": runs,
               "condition": condition, "unstable_mean": unstable_mean, "statistics": {},
               "A2_check": variance_decomposition(shares["t"], per_run["centred"]["R1"],
                                                  per_run["uncentred"]["R1"])}
        for label, (_, spectrum, order) in STATISTICS.items():
            obs = amplification(per_run["centred"][label], per_run["uncentred"][label])
            share = row["t_measured"] if spectrum == "singular" else row["f_measured"]
            centred_mean = float(np.mean(per_run["centred"][label]))
            obs["predicted_amplification"] = predicted_amplification(order, share, centred_mean)
            obs["predicted_amplification_is_infinite"] = (
                order == math.inf and obs["predicted_amplification"] is None
                and 0.0 < share < 1.0)
            obs["naive_1_over_1_minus_f"] = (float(1.0 / (1.0 - row["f_measured"]))
                                             if row["f_measured"] < 1.0 else None)
            obs["centred_mean"] = centred_mean
            row["statistics"][label] = obs
        out.append(row)
        print(f"  f={f:.2f} t={row['t_measured']:.4f} "
              f"A1_obs={row['statistics']['R1']['A_from_sd_log']} "
              f"A1_pred={row['statistics']['R1']['predicted_amplification']}", flush=True)
    return out


# --------------------------------------------------------------------------
# real data
# --------------------------------------------------------------------------
def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def score_repeat(path: str, seed: int = 42) -> dict:
    """Every view, raw and residualised, for one exported repeat."""
    raw = np.load(path, allow_pickle=False)
    test = np.asarray(raw["split"]).astype(str) == "test"
    ids = np.asarray(raw["patient_ids"]).astype(str)[test]
    cancers = np.asarray(raw["cancers"]).astype(str)[test]
    tss, _ = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])

    record = {"path": str(Path(path).resolve()), "sha256": _sha256(path),
              "n_test": int(test.sum()), "views": {}, "guards": {}}
    for view in VIEWS:
        x = np.asarray(raw[view], dtype=np.float64)[test]
        blocks = {"raw": x, "residualised": cross_fitted_residuals(x, design, seed=seed)}
        record["views"][view] = {b: score_matrix(m) for b, m in blocks.items()}
        record["guards"][view] = check_doubling_is_exact(blocks["raw"])
    return record


def collect(reps: dict[str, dict]) -> dict:
    """Across-repeat dispersion, observed against predicted amplification."""
    names = sorted(reps)
    out: dict = {}
    for view in VIEWS:
        out[view] = {}
        for block in ("raw", "residualised"):
            cell = {"shares": {}, "identity": {}, "statistics": {}}
            for key in ("f", "t", "mean_offdiag_cosine", "t_from_f", "t_identity_absolute_error"):
                cell["shares"][key] = {r: reps[r]["views"][view][block][key] for r in names}
            cell["identity"] = {r: reps[r]["views"][view][block]["identity"] for r in names}
            for label in list(STATISTICS) + ["RankMe"]:
                per = [reps[r]["views"][view][block]["statistics"][label] for r in names]
                centred = [p["centred"] for p in per]
                uncentred = [p["uncentred"] for p in per]
                entry = (amplification(centred, uncentred) if centred[0] is not None
                         else {"centred": None, "uncentred": dispersion(uncentred),
                               "A_from_sd_log": None, "A_from_log_fold": None})
                preds = [p["predicted_amplification"] for p in per]
                entry["predicted_amplification_per_rep"] = preds
                finite = [p for p in preds if p is not None and np.isfinite(p)]
                entry["predicted_amplification"] = float(np.mean(finite)) if finite else None
                entry["predicted_amplification_is_infinite"] = bool(
                    any(p["predicted_amplification_is_infinite"] for p in per))
                cell["statistics"][label] = entry
            cell["A2_check"] = variance_decomposition(
                [reps[r]["views"][view][block]["t"] for r in names],
                [reps[r]["views"][view][block]["statistics"]["R1"]["centred"] for r in names],
                [reps[r]["views"][view][block]["statistics"]["R1"]["uncentred"] for r in names])
            # the same, with the divergent repeat removed -- n=5 with one outlier is
            # otherwise a statement about one run
            rest = [r for r in names if r != "rep2"]
            cell["excluding_rep2"] = {}
            for label in STATISTICS:
                per = [reps[r]["views"][view][block]["statistics"][label] for r in rest]
                cell["excluding_rep2"][label] = amplification([p["centred"] for p in per],
                                                              [p["uncentred"] for p in per])
            out[view][block] = cell
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", default=None)
    ap.add_argument("--output", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--synthetic-only", action="store_true")
    ap.add_argument("--skip-synthetic", action="store_true")
    ap.add_argument("--synthetic-runs", type=int, default=25)
    ap.add_argument("--conditions", nargs="+", default=list(CONDITIONS),
                    help="run a subset of the synthetic conditions, so the three can go "
                         "in parallel processes at one thread each")
    ap.add_argument("--synthetic-n", type=int, default=2766)
    ap.add_argument("--synthetic-d", type=int, default=256)
    a = ap.parse_args(argv)

    out: dict = {"_": "Predeclared at NOTEBOOK_ENTRIES/PREDECLARED_centring_amplification_law_"
                      "20260804T1750Z.md, commit d4e344c. n = 5 same-seed repeats, no interval.",
                 "_config": {"seed": a.seed, "untrained_40_groups": list(UNTRAINED_40),
                             "alpha_index_range": list(ALPHA_INDEX_RANGE)}}

    if not a.skip_synthetic:
        out["synthetic"] = {"config": {"n": a.synthetic_n, "d": a.synthetic_d,
                                       "runs": a.synthetic_runs, "seed": a.seed,
                                       "conditions": list(a.conditions)}}
        for condition in a.conditions:
            print(f"synthetic, {condition}:", flush=True)
            # the seed offset is keyed to the condition, NOT to its position in the
            # requested subset, so splitting the three across processes reproduces
            # exactly what running them together would have produced
            out["synthetic"][condition] = synthetic_sweep(
                n=a.synthetic_n, d=a.synthetic_d, runs=a.synthetic_runs,
                seed=a.seed + CONDITIONS.index(condition), condition=condition)

    if not a.synthetic_only:
        paths = sorted(glob.glob(a.reps))
        if not paths:
            raise SystemExit(f"no repeats matched {a.reps!r}")
        reps = {}
        for p in paths:
            name = Path(p).stem
            print(f"scoring {name} ...", flush=True)
            reps[name] = score_repeat(p, seed=a.seed)
        out["repeats"] = reps
        out["collected"] = collect(reps)

    Path(a.output).expanduser().parent.mkdir(parents=True, exist_ok=True)
    Path(a.output).expanduser().write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("wrote", a.output)
    return out


if __name__ == "__main__":
    main()
