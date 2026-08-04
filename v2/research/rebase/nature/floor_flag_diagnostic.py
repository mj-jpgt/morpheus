"""Diagnostic: does ``observed_above_floor`` report the truth on cases with known answers?

Constructed cases only -- no real data. The planted association is built with
``calibra.calibration.spike_targets`` itself, i.e. with the *identical*
construction the detection floor is calibrated against, so the true
single-direction correlation carried by the "real" data is known exactly and is
on the floor's own scale by construction. That removes the units argument from
the question entirely: if a dataset carrying a planted single-direction
association of 0.60 is graded against a floor of 0.05 and the flag reads FAIL,
the flag is broken irrespective of what units anything is in.

Every statistic is imported from ``v2/calibra`` or ``v2/calibra/spectral``.
Nothing is computed inline except the construction of the synthetic matrices.

Predeclaration: NOTEBOOK_ENTRIES/PREDECLARED_observed_above_floor_20260804T1843Z.md
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))

from morpheus.v2.calibra.calibration import spike_recovery_curve, spike_targets  # noqa: E402
from morpheus.v2.calibra.residualise import cross_fitted_residuals  # noqa: E402
from morpheus.v2.calibra.spectral import (heldout_single_direction_correlation,  # noqa: E402
                                          heldout_top_cca, paired_absolute_correlation)

LEVELS = (0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60)


def make_cohort(n: int, p: int, q: int, n_sites: int, rng: np.random.Generator):
    """X, Y, design, strata with a confound that acts on BOTH modalities.

    The confound matters: the whole reason the floor is not zero is that
    residualising through a shared design manufactures correlation between
    orthogonal signals, so a diagnostic on an unconfounded cohort would not
    exercise the code path under test.
    """
    site = rng.integers(0, n_sites, size=n)
    cancer = rng.integers(0, 5, size=n)
    design = np.hstack([(site[:, None] == np.arange(n_sites)[None, :]).astype(np.float64),
                        (cancer[:, None] == np.arange(5)[None, :]).astype(np.float64)])
    site_x = rng.normal(size=(n_sites, p))
    site_y = rng.normal(size=(n_sites, q))
    cancer_x = rng.normal(size=(5, p))
    cancer_y = rng.normal(size=(5, q))
    x = rng.normal(size=(n, p)) + 1.5 * site_x[site] + 1.5 * cancer_x[cancer]
    y = rng.normal(size=(n, q)) + 1.5 * site_y[site] + 1.5 * cancer_y[cancer]
    return x, y, design, cancer


def plant(x, y, rho: float, seed: int, negative: bool):
    """Return (y_planted, u0, v0): y carrying corr(x@u0, y@v0) == +/- rho exactly."""
    planted, u0, v0 = spike_targets(x, y, abs(rho), rng=np.random.default_rng(seed),
                                    return_directions=True)
    if negative:
        # Exact reflection in the v0 coordinate: flips the sign of the association
        # along (u0, v0) and leaves its magnitude, and every other direction, alone.
        planted = planted - 2.0 * np.outer(planted @ v0, v0)
    return planted, u0, v0


def case(name: str, *, n: int, p: int, q: int, rho: float, negative: bool, seed: int,
         n_draws: int, n_components: int, n_jobs: int, n_sites: int = 20) -> dict:
    rng = np.random.default_rng(seed)
    x, y, design, strata = make_cohort(n, p, q, n_sites, rng)
    y_real, u0, v0 = plant(x, y, rho, seed, negative)

    curve = spike_recovery_curve(x, y_real, design, levels=LEVELS, n_draws=n_draws,
                                 n_components=min(n_components, p, q), seed=seed, n_jobs=n_jobs)
    summary = curve.summary()
    floor = summary["detection_floor"]

    x_res = cross_fitted_residuals(x, design, seed=seed)
    y_res = cross_fitted_residuals(y_real, design, seed=seed)

    # TRUTH, on the floor's own scale: the correlation along the planted direction
    # pair after the identical residualisation. Not a fitted direction -- the pair
    # is the one the signal was planted on, known before any data was looked at.
    truth = paired_absolute_correlation(x_res @ u0, y_res @ v0)

    # CANDIDATE CORRECTED COMPARATORS, both out-of-fold / held-out, both from the library.
    heldout = heldout_top_cca(x_res, y_res, n_components=min(n_components, p, q), seed=seed)
    single = heldout_single_direction_correlation(x_res, y_res @ v0, seed=seed)

    # G1, destroyed pairing: permute y rows within cancer strata and repeat.
    order = np.arange(n)
    permute_rng = np.random.default_rng(seed + 9_000)
    for level in np.unique(strata):
        idx = np.flatnonzero(strata == level)
        order[idx] = permute_rng.permutation(idx)
    y_perm_res = cross_fitted_residuals(y_real[order], design, seed=seed)
    heldout_permuted = heldout_top_cca(x_res, y_perm_res, n_components=min(n_components, p, q),
                                       seed=seed)
    truth_permuted = paired_absolute_correlation(x_res @ u0, y_perm_res @ v0)

    matched = summary["observed_matched_direction"]
    return {
        "case": name, "n": n, "p": p, "q": q, "seed": seed,
        "planted_rho": rho, "negative": negative,
        "detection_floor": floor,
        "transmission_floor": summary["transmission_floor"],
        "confound_induced_baseline": summary["confound_induced_baseline"],
        "truth_planted_direction_abs": truth,
        "truth_over_floor": float(truth / floor) if np.isfinite(floor) and floor > 0 else float("nan"),
        "shipped_observed_matched_direction": matched,
        "shipped_observed_matched_direction_abs": float(abs(matched)),
        "shipped_observed_above_floor": bool(summary["observed_above_floor"]),
        "corrected_flag_truth_direction": bool(np.isfinite(floor) and truth > floor),
        "heldout_top_cca": heldout,
        "corrected_flag_heldout": bool(np.isfinite(floor) and np.isfinite(heldout) and heldout > floor),
        "heldout_single_direction_on_v0": single,
        "G1_heldout_top_cca_pairing_destroyed": heldout_permuted,
        "G1_truth_direction_pairing_destroyed": truth_permuted,
        "G1_clears": bool(np.isfinite(floor) and np.isfinite(heldout_permuted)
                          and heldout_permuted > floor),
        "observed_multivariate_top_cca": summary["observed_multivariate_top_cca"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n", type=int, default=800)
    parser.add_argument("--n-draws", type=int, default=40)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    results: list[dict] = []

    # --- A-E: multivariate channel, both signs, above / at / below the floor ---------
    grid = [("A_no_association", 0.0, False),
            ("B_far_above_floor_positive", 0.60, False),
            ("C_far_above_floor_negative", 0.60, True),
            ("D_mid", 0.20, False),
            ("E_far_below_floor", 0.01, False)]
    for name, rho, negative in grid:
        results.append(case(name, n=args.n, p=16, q=16, rho=rho, negative=negative,
                            seed=args.seed, n_draws=args.n_draws,
                            n_components=args.n_components, n_jobs=args.n_jobs))
        print(f"[{name}] " + json.dumps({k: results[-1][k] for k in (
            "detection_floor", "truth_planted_direction_abs", "truth_over_floor",
            "shipped_observed_matched_direction", "shipped_observed_above_floor",
            "corrected_flag_truth_direction", "heldout_top_cca",
            "G1_heldout_top_cca_pairing_destroyed")}), flush=True)

    # --- B3: seed lottery on the ONE-COLUMN case P4 certification actually uses ------
    for seed in range(args.seed, args.seed + 8):
        results.append(case(f"F_single_column_seed{seed}", n=args.n, p=1, q=1, rho=0.60,
                            negative=False, seed=seed, n_draws=args.n_draws,
                            n_components=1, n_jobs=args.n_jobs))
        print(f"[F seed={seed}] " + json.dumps({k: results[-1][k] for k in (
            "detection_floor", "truth_planted_direction_abs",
            "shipped_observed_matched_direction", "shipped_observed_above_floor")}), flush=True)

    # --- B4: does the shipped comparator move with the real channel's strength? ------
    for rho in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
        results.append(case(f"G_sweep_rho{rho:g}", n=args.n, p=16, q=16, rho=rho,
                            negative=False, seed=args.seed, n_draws=args.n_draws,
                            n_components=args.n_components, n_jobs=args.n_jobs))
        print(f"[G rho={rho:g}] " + json.dumps({k: results[-1][k] for k in (
            "truth_planted_direction_abs", "shipped_observed_matched_direction",
            "shipped_observed_matched_direction_abs", "heldout_top_cca",
            "shipped_observed_above_floor")}), flush=True)

    # --- G4: seed stability of the corrected verdict on the far-above-floor case -----
    for seed in range(args.seed, args.seed + 4):
        results.append(case(f"H_stability_seed{seed}", n=args.n, p=16, q=16, rho=0.60,
                            negative=False, seed=seed, n_draws=args.n_draws,
                            n_components=args.n_components, n_jobs=args.n_jobs))
        print(f"[H seed={seed}] " + json.dumps({k: results[-1][k] for k in (
            "detection_floor", "truth_planted_direction_abs",
            "shipped_observed_above_floor", "corrected_flag_truth_direction",
            "corrected_flag_heldout")}), flush=True)

    Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"[written] {args.out} ({len(results)} cases)", flush=True)


if __name__ == "__main__":
    main()
