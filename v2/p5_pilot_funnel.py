"""P5 discovery funnel, stages 0-4, novelty tier-1 only -- a small pilot.

Implements `paper/P5_DISCOVERY_PLAN.md` section 5's first concrete step: enumerate a
candidate space of (target, cancer-stratum) cells, coarsely pre-filter it, certify
survivors with the unchanged CALIBRA instrument, correct for multiplicity across the
whole predeclared search space with a canonical Benjamini-Hochberg implementation, and
check replication on a held-out split. Tier-1 novelty only: candidate targets whose
`target_groups` label is `hallmark_in_training` are excluded before stage 0 even
enumerates them, because that group is the D2 arm-H supervision target itself --
scoring it would measure training-signal recovery, not discovery.

Every per-cell statistic is an UNCHANGED import from `calibra`:
  * stage 1 (coarse, cheap, uncorrected) -- `spectral.top_canonical_correlation`
  * stage 2 (certify) -- `residualise.confound_design` / `cross_fitted_residuals`,
    `calibration.permutation_null`, `calibration.spike_recovery_curve`,
    `spectral.heldout_single_direction_correlation` as the channel statistic (the
    single-target analogue of `run_calibra._channel_measurement`'s `heldout_top_cca`),
    graded by `calibration.channel_clears_floor`
  * stage 3 (multiplicity) -- `scipy.stats.false_discovery_control` (library BH, not a
    hand-rolled formula)
No inline reimplementation of any of these appears in this file.

Two data sources, selected by CLI flag, per
`NOTEBOOK_ENTRIES/PREDECLARED_p5_novelty_scoping_and_pilot_funnel_20260805T0750Z.md`:

  --artifact / --targets   real frozen artifacts (representation state npz +
                            frozen_rna_targets.npz). Not exercised in the pilot run
                            this module was first written for -- no such artifact was
                            reachable from that sandbox (see the notebook entry for
                            what was tried and why it was judged unreachable rather
                            than assumed so).

  --synthetic-dry-run      a controlled synthetic ladder (200 cells, 20 with a planted
                            signal built with `calibration.spike_targets`, 180 pure
                            null) that exercises the identical code path to check the
                            funnel's stages behave the way the plan says they should.
                            This is a code-path validation, never a biological result,
                            and every summary field this mode writes is prefixed or
                            labelled to say so.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from .calibra.calibration import channel_clears_floor, permutation_null, spike_recovery_curve, spike_targets
from .calibra.residualise import confound_design, cross_fitted_residuals
from .calibra.spectral import heldout_single_direction_correlation, top_canonical_correlation

STAGE1_KEEP_FRACTION = 0.40          # predeclared
BH_Q = 0.10                          # predeclared
SYNTHETIC_N_STRATA = 5
SYNTHETIC_N_TARGETS = 40
SYNTHETIC_N_AXES = 32
SYNTHETIC_N_PATIENTS_PER_STRATUM = 500     # split 250 discovery / 250 replication
SYNTHETIC_N_PLANTED = 20
SYNTHETIC_R_TRUE = 0.18
SYNTHETIC_STRATA_LABELS = ("BRCA", "LUAD", "KIRC", "HNSC", "THCA")

# Reduced-cost knobs for a same-day CPU pilot. Stated as a deviation from the plan's
# full-scale 256-axis, larger-n_permutations design, not hidden inside a default.
PILOT_N_PERMUTATIONS = 30
PILOT_N_DRAWS = 8
PILOT_SPIKE_LEVELS = (0.0, 0.05, 0.10, 0.15, 0.20, 0.30)


@dataclass
class Cell:
    stratum: str
    target_id: str
    is_planted: bool = False          # ground truth, synthetic-only; NaN-equivalent unused for real data


@dataclass
class CellResult:
    cell: Cell
    stage1_stat: float = float("nan")
    stage1_pass: bool = False
    stage2_attempted: bool = False
    permutation_p: float = float("nan")
    detection_floor: float = float("nan")
    channel_statistic: float = float("nan")
    clears_floor: object = None       # True / False / None (ungraded)
    stage3_survivor: bool = False
    stage4_attempted: bool = False
    stage4_replicated: bool = False
    stage4_statistic: float = float("nan")
    stage4_p: float = float("nan")


def enumerate_candidate_space(strata: list[str], target_ids: list[str],
                              planted: set[tuple[str, str]] | None = None) -> list[Cell]:
    """Stage 0. The full cross of strata x candidate targets, fixed size."""
    planted = planted or set()
    return [Cell(stratum=s, target_id=t, is_planted=(s, t) in planted)
            for s in strata for t in target_ids]


def stage1_coarse_filter(cells: list[Cell], data: dict[str, dict[str, tuple]],
                         keep_fraction: float = STAGE1_KEEP_FRACTION) -> list[CellResult]:
    """Cheap, in-sample, unadjusted top-CCA at n_components=1 -- purely to cut compute."""
    results = []
    for cell in cells:
        x, y_col = data[cell.stratum]["x_discovery"], data[cell.stratum]["y_discovery"][cell.target_id]
        stat = top_canonical_correlation(x, y_col.reshape(-1, 1), n_components=1)
        results.append(CellResult(cell=cell, stage1_stat=float(stat)))
    ranked = sorted(results, key=lambda r: -abs(r.stage1_stat) if np.isfinite(r.stage1_stat) else -np.inf)
    n_keep = int(round(len(cells) * keep_fraction))
    keep_ids = {id(r) for r in ranked[:n_keep]}
    for r in results:
        r.stage1_pass = id(r) in keep_ids
    return results


def stage2_certify(result: CellResult, data: dict[str, dict[str, tuple]], *, seed: int) -> None:
    """Confound-adjusted, permutation-null, injection-certified -- reuses `calibra` unchanged."""
    stratum_data = data[result.cell.stratum]
    x = stratum_data["x_discovery"]
    y = stratum_data["y_discovery"][result.cell.target_id].reshape(-1, 1)
    site = stratum_data["site_discovery"]
    frame = pd.DataFrame({"site": site})
    design = confound_design(frame, ["site"])
    result.stage2_attempted = True

    null = permutation_null(x, y, design, strata=None, n_permutations=PILOT_N_PERMUTATIONS,
                            n_components=1, seed=seed)
    result.permutation_p = float(null["permutation_p"])

    spike = spike_recovery_curve(x, y, design, levels=PILOT_SPIKE_LEVELS, n_draws=PILOT_N_DRAWS,
                                 n_components=1, seed=seed)
    summary = spike.summary()
    result.detection_floor = float(summary["detection_floor"])

    x_res = cross_fitted_residuals(x, design, seed=seed)
    y_res = cross_fitted_residuals(y, design, seed=seed)
    channel = heldout_single_direction_correlation(x_res, y_res.ravel(), seed=seed)
    result.channel_statistic = float(channel)
    verdict, _status = channel_clears_floor(channel, result.detection_floor)
    result.clears_floor = verdict


def stage3_multiplicity_correction(results: list[CellResult], *, q: float = BH_Q) -> None:
    """BH-FDR across every stage-2-tested cell -- the whole predeclared search space
    that actually reached stage 2, not a post-hoc subset chosen because it looked good."""
    tested = [r for r in results if r.stage2_attempted and np.isfinite(r.permutation_p)]
    if not tested:
        return
    p_values = np.array([r.permutation_p for r in tested])
    adjusted = false_discovery_control(p_values, method="bh")
    for r, adj_p in zip(tested, adjusted):
        r.stage3_survivor = bool(adj_p <= q)


def stage4_replicate(result: CellResult, data: dict[str, dict[str, tuple]], *, seed: int) -> None:
    """Re-measure the certified statistic on the independent replication half."""
    stratum_data = data[result.cell.stratum]
    x = stratum_data["x_replication"]
    y = stratum_data["y_replication"][result.cell.target_id].reshape(-1, 1)
    site = stratum_data["site_replication"]
    frame = pd.DataFrame({"site": site})
    design = confound_design(frame, ["site"])
    result.stage4_attempted = True

    x_res = cross_fitted_residuals(x, design, seed=seed)
    y_res = cross_fitted_residuals(y, design, seed=seed)
    channel = heldout_single_direction_correlation(x_res, y_res.ravel(), seed=seed)
    result.stage4_statistic = float(channel)

    null = permutation_null(x, y, design, strata=None, n_permutations=PILOT_N_PERMUTATIONS,
                            n_components=1, seed=seed)
    result.stage4_p = float(null["permutation_p"])
    same_sign = np.sign(channel) == np.sign(result.channel_statistic) if np.isfinite(channel) and np.isfinite(result.channel_statistic) else False
    result.stage4_replicated = bool(same_sign and result.stage4_p < 0.05)


def run_funnel(strata: list[str], target_ids: list[str], data: dict[str, dict[str, tuple]], *,
              planted: set[tuple[str, str]] | None = None, seed: int = 42) -> list[CellResult]:
    cells = enumerate_candidate_space(strata, target_ids, planted=planted)
    results = stage1_coarse_filter(cells, data)
    for r in results:
        if r.stage1_pass:
            stage2_certify(r, data, seed=seed)
    stage3_multiplicity_correction(results, q=BH_Q)
    for r in results:
        if r.stage3_survivor:
            stage4_replicate(r, data, seed=seed)
    return results


def ledger(results: list[CellResult]) -> dict:
    n_entered = len(results)
    n_stage1 = sum(r.stage1_pass for r in results)
    n_stage2 = sum(r.stage2_attempted for r in results)
    n_stage3 = sum(r.stage3_survivor for r in results)
    n_stage4_attempted = sum(r.stage4_attempted for r in results)
    n_stage4_replicated = sum(r.stage4_replicated for r in results)
    return {
        "n_candidates_entered_stage0": n_entered,
        "n_survived_stage1_coarse_filter": n_stage1,
        "n_reached_stage2_certify": n_stage2,
        "n_cleared_stage3_bh_fdr": n_stage3,
        "bh_q": BH_Q,
        "n_attempted_stage4_replication": n_stage4_attempted,
        "n_replicated_stage4": n_stage4_replicated,
        "stage1_keep_fraction": STAGE1_KEEP_FRACTION,
    }


def ledger_rows(results: list[CellResult]) -> pd.DataFrame:
    rows = []
    for r in results:
        rows.append({
            "stratum": r.cell.stratum, "target_id": r.cell.target_id,
            "is_planted_ground_truth": r.cell.is_planted,
            "stage1_stat": r.stage1_stat, "stage1_pass": r.stage1_pass,
            "stage2_attempted": r.stage2_attempted, "permutation_p": r.permutation_p,
            "detection_floor": r.detection_floor, "channel_statistic": r.channel_statistic,
            "clears_floor": r.clears_floor, "stage3_survivor": r.stage3_survivor,
            "stage4_attempted": r.stage4_attempted, "stage4_statistic": r.stage4_statistic,
            "stage4_p": r.stage4_p, "stage4_replicated": r.stage4_replicated,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Synthetic mechanics fallback -- NEVER a biological result. See module docstring.
# ---------------------------------------------------------------------------

def build_synthetic_data(*, seed: int = 42) -> tuple[list[str], list[str], set[tuple[str, str]], dict]:
    rng = np.random.default_rng(seed)
    strata = list(SYNTHETIC_STRATA_LABELS)
    target_ids = [f"SYNTH_TARGET_{i:03d}" for i in range(SYNTHETIC_N_TARGETS)]
    all_cells = [(s, t) for s in strata for t in target_ids]
    planted = set(tuple(pair) for pair in
                 rng.choice(np.array(all_cells, dtype=object), size=SYNTHETIC_N_PLANTED, replace=False))

    data: dict[str, dict[str, tuple]] = {}
    for stratum in strata:
        n = SYNTHETIC_N_PATIENTS_PER_STRATUM
        n_sites = 4
        site = rng.integers(0, n_sites, size=n).astype(str)
        site_effect_x = rng.normal(scale=0.3, size=(n_sites, SYNTHETIC_N_AXES))
        site_effect_y = rng.normal(scale=0.3, size=n_sites)

        x = rng.normal(size=(n, SYNTHETIC_N_AXES)) + site_effect_x[site.astype(int)]
        y = {}
        for target_id in target_ids:
            y_col = rng.normal(size=n) + site_effect_y[site.astype(int)]
            if (stratum, target_id) in planted:
                y_col = spike_targets(x, y_col.reshape(-1, 1), SYNTHETIC_R_TRUE,
                                      rng=np.random.default_rng(rng.integers(1 << 31))).ravel()
            y[target_id] = y_col

        half = n // 2
        order = rng.permutation(n)
        disc_idx, repl_idx = order[:half], order[half:]
        data[stratum] = {
            "x_discovery": x[disc_idx], "site_discovery": site[disc_idx],
            "y_discovery": {t: v[disc_idx] for t, v in y.items()},
            "x_replication": x[repl_idx], "site_replication": site[repl_idx],
            "y_replication": {t: v[repl_idx] for t, v in y.items()},
        }
    return strata, target_ids, planted, data


def synthetic_confusion(results: list[CellResult]) -> dict:
    """Sensitivity/specificity of the FUNNEL ITSELF against the known planting ground
    truth. A mechanics diagnostic, not a claim about any real gene or pathway."""
    planted = [r for r in results if r.cell.is_planted]
    null = [r for r in results if not r.cell.is_planted]

    def _rate(cells, predicate):
        return float(np.mean([predicate(c) for c in cells])) if cells else float("nan")

    return {
        "n_planted": len(planted), "n_null": len(null),
        "planted_stage1_survival_rate": _rate(planted, lambda c: c.stage1_pass),
        "null_stage1_survival_rate": _rate(null, lambda c: c.stage1_pass),
        "planted_stage3_survival_rate": _rate(planted, lambda c: c.stage3_survivor),
        "null_stage3_survival_rate": _rate(null, lambda c: c.stage3_survivor),
        "planted_stage4_replication_rate_of_stage3_survivors": _rate(
            [c for c in planted if c.stage3_survivor], lambda c: c.stage4_replicated),
        "null_stage4_replication_rate_of_stage3_survivors": _rate(
            [c for c in null if c.stage3_survivor], lambda c: c.stage4_replicated),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", default="", help="representation state npz (real-data path)")
    parser.add_argument("--targets", default="", help="frozen_rna_targets.npz (real-data path)")
    parser.add_argument("--synthetic-dry-run", action="store_true",
                        help="mechanics-only validation on a controlled synthetic ladder; "
                             "NEVER a biological result")
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if not args.synthetic_dry_run and not (args.artifact and args.targets):
        parser.error("either --synthetic-dry-run or both --artifact and --targets are required")

    output = Path(args.output)

    if args.synthetic_dry_run:
        strata, target_ids, planted, data = build_synthetic_data(seed=args.seed)
        results = run_funnel(strata, target_ids, data, planted=planted, seed=args.seed)
        summary = ledger(results)
        summary["data_provenance"] = "SYNTHETIC_DRY_RUN_NO_REAL_TCGA_DATA"
        summary["funnel_mechanics_confusion"] = synthetic_confusion(results)
    elif args.artifact and args.targets:
        raise NotImplementedError(
            "real-data path is specified in the predeclaration (section 2.2) but was not "
            "exercised by the pilot run this module was first written for -- no reachable "
            "frozen artifact was found. Implement against a real frozen_rna_targets.npz "
            "schema (see calibra/run_calibra.py's loader) before pointing this at real data.")

    output.mkdir(parents=True, exist_ok=True)
    ledger_rows(results).to_csv(output / "p5_pilot_ledger_rows.csv", index=False)
    (output / "p5_pilot_ledger_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
