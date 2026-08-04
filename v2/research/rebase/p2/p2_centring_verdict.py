"""Grade the centring-amplification law against the thresholds fixed BEFORE it was run.

Reads the outputs of `p2_centring_amplification.py` and returns a verdict per
falsifier. **Every threshold here is a transcription of
`NOTEBOOK_ENTRIES/PREDECLARED_centring_amplification_law_20260804T1750Z.md` §4
and §5, commit `d4e344c`, and `v2/tests/test_p2_centring_amplification.py`
asserts that the constants below match the ones in that file.** The point of the
module is that the verdicts cannot be re-read after the numbers arrive: a
threshold edited here breaks the build.

A verdict is one of:
  * ``"falsified"``   -- the predeclared condition for falsification was met;
  * ``"not falsified"`` -- it was not;
  * ``"not applicable"`` -- the quantity the falsifier needs does not exist
    (e.g. a dispersion ratio when neither dispersion is non-zero). This is
    neither a pass nor a failure and is counted separately, exactly as
    `p2_floor_audit.py` counts an unjudgeable row.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

#: Every threshold, transcribed from the predeclaration. The test asserts the
#: numbers appear verbatim in that file.
THRESHOLDS = {
    "P1_relative_error_limit": 0.05,
    "P1_max_cells_allowed_over": 2,
    "P2_relative_error_limit": 0.50,
    "P3_spearman_floor": 0.0,
    "P3b_spike_share_limit": 0.50,
    "P5_uncentred_share_limit": 0.70,
    "S1_median_relative_error_limit": 0.10,
    "S1_single_point_limit": 0.20,
    "S1_f_ceiling": 0.90,
    "S3_max_points_allowed_over": 2,
    "S4_max_points_allowed_over": 2,
}

VIEWS = ("wsi_biology", "rna_biology", "full_biology")


def _spearman(a, b) -> float:
    """Rank correlation, computed on ranks with `numpy.argsort` -- no scipy on the box."""
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    if a.size < 3:
        return float("nan")
    ra = np.empty(a.size); ra[np.argsort(a)] = np.arange(a.size)
    rb = np.empty(b.size); rb[np.argsort(b)] = np.arange(b.size)
    return float(np.corrcoef(ra, rb)[0, 1])


def _jackknife_amplification(centred, uncentred) -> tuple[float, float] | None:
    """Leave-one-repeat-out range of ``sd(ln centred)/sd(ln uncentred)``.

    n = 5 with one divergent run is not enough for a bootstrap to mean anything;
    the leave-one-out range at least says whether a verdict survives dropping any
    single run. It is a range, not an interval, and is never called one.
    """
    c, u = np.asarray(centred, dtype=np.float64), np.asarray(uncentred, dtype=np.float64)
    if c.size < 4 or (c <= 0).any() or (u <= 0).any():
        return None
    out = []
    for i in range(c.size):
        cc = np.log(np.delete(c, i)); uu = np.log(np.delete(u, i))
        du = uu.std(ddof=1)
        if du <= 0:
            return None
        out.append(cc.std(ddof=1) / du)
    return float(min(out)), float(max(out))


def evaluate_real(real: dict) -> dict:
    """P1, P2, P3, P3b, P4, P5 against the predeclared thresholds."""
    reps, col = real["repeats"], real["collected"]
    names = sorted(reps)
    out: dict = {}

    # --- P1: does the transfer identity hold? (tests assumption A1) -----------
    cells = []
    for view in VIEWS:
        for block in ("raw", "residualised"):
            for r in names:
                err = reps[r]["views"][view][block]["identity"]["relative_error"]
                cells.append({"view": view, "block": block, "rep": r,
                              "relative_error": err})
    over = [c for c in cells if c["relative_error"] is not None
            and c["relative_error"] > THRESHOLDS["P1_relative_error_limit"]]
    finite = [c["relative_error"] for c in cells if c["relative_error"] is not None]
    out["P1"] = {"cells": len(cells), "cells_over_limit": len(over),
                 "worst_relative_error": max(finite) if finite else None,
                 "worst_on_raw_block": max(c["relative_error"] for c in cells
                                           if c["block"] == "raw"),
                 "verdict": ("falsified" if len(over) > THRESHOLDS["P1_max_cells_allowed_over"]
                             else "not falsified")}

    # --- P2: the law itself, order 1, raw block ------------------------------
    rows = {}
    for view in VIEWS:
        e = col[view]["raw"]["statistics"]["R1"]
        obs, pred = e["A_from_sd_log"], e["predicted_amplification"]
        rel = (abs(obs - pred) / pred) if (obs is not None and pred) else None
        rows[view] = {"A_observed_sd_log": obs, "A_observed_log_fold": e["A_from_log_fold"],
                      "A_predicted": pred, "relative_error": rel,
                      "t_mean": float(np.mean(list(col[view]["raw"]["shares"]["t"].values()))),
                      "f_mean": float(np.mean(list(col[view]["raw"]["shares"]["f"].values())))}
        rows[view]["naive_1_over_1_minus_f"] = (
            1.0 / (1.0 - rows[view]["f_mean"]) if rows[view]["f_mean"] < 1 else None)
    by_obs = sorted(VIEWS, key=lambda v: -rows[v]["A_observed_sd_log"])
    by_pred = sorted(VIEWS, key=lambda v: -rows[v]["A_predicted"])
    any_outside = any(r["relative_error"] is not None
                      and r["relative_error"] > THRESHOLDS["P2_relative_error_limit"]
                      for r in rows.values())
    out["P2"] = {"per_view": rows, "order_by_observed": by_obs, "order_by_predicted": by_pred,
                 "order_agrees": by_obs == by_pred,
                 "verdict": "falsified" if (any_outside or by_obs != by_pred) else "not falsified"}

    # --- P3: does the order-a prediction track across the statistic family? ---
    pts = []
    for view in VIEWS:
        for label, entry in col[view]["raw"]["statistics"].items():
            obs, pred = entry["A_from_sd_log"], entry["predicted_amplification"]
            if obs is None or pred is None or not np.isfinite(pred):
                continue
            pts.append({"view": view, "statistic": label, "observed": obs, "predicted": pred})
    rho = _spearman([p["predicted"] for p in pts], [p["observed"] for p in pts])
    within = {v: _spearman([p["predicted"] for p in pts if p["view"] == v],
                           [p["observed"] for p in pts if p["view"] == v]) for v in VIEWS}
    out["P3"] = {"n_points": len(pts), "spearman": rho, "within_view_spearman": within,
                 "points": pts,
                 "verdict": ("falsified" if (not np.isfinite(rho)
                                             or rho <= THRESHOLDS["P3_spearman_floor"])
                             else "not falsified")}

    # --- P3b: does A2 -- "the dominant component is stable" -- hold here? -----
    a2 = {}
    for view in VIEWS:
        c = col[view]["raw"]["A2_check"]
        share = c.get("spike_share_of_uncentred_variance") if c.get("defined") else None
        a2[view] = {"spike_share_of_uncentred_variance": share,
                    "sd_t": c.get("sd_t"), "t_mean": c.get("t_mean"),
                    "relative_sd_of_t": (c["sd_t"] / c["t_mean"]
                                         if c.get("defined") and c.get("t_mean") else None),
                    "coherent": (share is not None and 0.0 <= share <= 1.0)}
    broken = [v for v, c in a2.items()
              if c["spike_share_of_uncentred_variance"] is not None
              and c["spike_share_of_uncentred_variance"] > THRESHOLDS["P3b_spike_share_limit"]]
    out["P3b"] = {"per_view": a2, "views_where_A2_is_broken": broken,
                  "verdict": "falsified" if broken else "not falsified",
                  "_": "'falsified' here means A2 is FALSE on this data, i.e. the law's own "
                       "premise fails -- not that the law's arithmetic is wrong."}

    # --- P4: the sign. Does centring ever IMPROVE reproducibility? -----------
    below = []
    for view in VIEWS:
        for label in col[view]["raw"]["statistics"]:
            entry = col[view]["raw"]["statistics"][label]
            obs = entry["A_from_sd_log"]
            if obs is None or obs >= 1.0:
                continue
            per = [reps[r]["views"][view]["raw"]["statistics"][label] for r in names]
            jack = _jackknife_amplification([p["centred"] for p in per],
                                            [p["uncentred"] for p in per])
            below.append({"view": view, "statistic": label, "A_observed": obs,
                          "leave_one_out_range": jack,
                          "survives_dropping_any_one_repeat": bool(jack and jack[1] < 1.0)})
    survivors = [b for b in below if b["survives_dropping_any_one_repeat"]]
    out["P4"] = {"cells_below_one": below, "cells_below_one_robust": survivors,
                 "verdict": "falsified" if survivors else "not falsified",
                 "_": "NOT a blind test -- the predeclaration records that R1 = 1.0228x against "
                      "RankMe = 1.0299x on rna_biology raw had already been read."}

    # --- P5: is the VIEW difference in floors a centring effect at all? ------
    logs = {}
    for view in VIEWS:
        e = col[view]["raw"]["statistics"]["R1"]
        logs[view] = {"uncentred": math.log(e["uncentred"]["fold"]),
                      "centred": math.log(e["centred"]["fold"]),
                      "residualised": math.log(
                          col[view]["residualised"]["statistics"]["R1"]["centred"]["fold"])}
    spread_c = max(l["centred"] for l in logs.values()) - min(l["centred"] for l in logs.values())
    spread_u = (max(l["uncentred"] for l in logs.values())
                - min(l["uncentred"] for l in logs.values()))
    share = spread_u / spread_c if spread_c > 0 else None
    out["P5"] = {"log_folds": logs, "cross_view_log_spread_centred": spread_c,
                 "cross_view_log_spread_uncentred": spread_u,
                 "fraction_of_the_view_effect_present_before_centring": share,
                 "verdict": ("falsified" if (share is not None
                                             and share > THRESHOLDS["P5_uncentred_share_limit"])
                             else "not falsified")}
    return out


def evaluate_synthetic(rows: list[dict]) -> dict:
    """S1, S2, S3, S4 for one condition."""
    ceiling = THRESHOLDS["S1_f_ceiling"]
    errs, naive_errs, s3_violations, s4_violations = [], [], 0, 0
    detail = []
    for row in rows:
        r1 = row["statistics"]["R1"]
        obs, pred = r1["A_from_sd_log"], r1["predicted_amplification"]
        naive = r1["naive_1_over_1_minus_f"]
        rel = abs(obs - pred) / pred if (obs is not None and pred) else None
        rel_naive = abs(obs - naive) / naive if (obs is not None and naive) else None
        if row["f_measured"] <= ceiling and rel is not None:
            errs.append(rel)
            if rel_naive is not None:
                naive_errs.append(rel_naive)
        if obs is not None and pred is not None and obs >= pred:
            s3_violations += 1
        a = {k: row["statistics"][k]["A_from_sd_log"] for k in
             ("R2", "R1", "stable_rank", "hard_rank")}
        ordered = [a["R2"], a["R1"], a["stable_rank"], a["hard_rank"]]
        if any(x is None for x in ordered) or any(
                ordered[i] < ordered[i + 1] - 1e-9 for i in range(3)):
            s4_violations += 1
        detail.append({"f": row["f_measured"], "t": row["t_measured"], "t_sd": row.get("t_sd"),
                       "A_observed": obs, "A_predicted": pred, "naive": naive,
                       "relative_error": rel, "naive_relative_error": rel_naive,
                       "centred_fold": r1["centred"]["fold"],
                       "uncentred_fold": r1["uncentred"]["fold"],
                       "spike_share": (row["A2_check"] or {}).get(
                           "spike_share_of_uncentred_variance")})
    median = float(np.median(errs)) if errs else None
    worst = max(errs) if errs else None
    median_naive = float(np.median(naive_errs)) if naive_errs else None
    return {
        "condition": rows[0]["condition"] if rows else None,
        "detail": detail,
        "S1": {"median_relative_error": median, "worst_relative_error": worst,
               "n_points_at_or_below_f_ceiling": len(errs),
               "verdict": ("falsified" if (median is None
                                           or median > THRESHOLDS["S1_median_relative_error_limit"]
                                           or worst > THRESHOLDS["S1_single_point_limit"])
                           else "not falsified")},
        "S2": {"median_relative_error_derived": median, "median_relative_error_naive": median_naive,
               "verdict": ("falsified" if (median_naive is None or median is None
                                           or median_naive <= median)
                           else "not falsified")},
        "S3": {"points_where_observed_at_or_above_predicted": s3_violations,
               "n_points": len(rows),
               "verdict": ("falsified"
                           if s3_violations > THRESHOLDS["S3_max_points_allowed_over"]
                           else "not falsified"),
               "_": "S3 is only meaningful on the unstable-mean condition; it is reported for "
                    "all three so the contrast is visible."},
        "S4": {"points_violating_A2_ge_A1_ge_Ainf_ge_A0": s4_violations, "n_points": len(rows),
               "verdict": ("falsified"
                           if s4_violations > THRESHOLDS["S4_max_points_allowed_over"]
                           else "not falsified")},
    }


def evaluate(real_path, synthetic_paths) -> dict:
    real = json.loads(Path(real_path).expanduser().read_text(encoding="utf-8"))
    out = {"thresholds": THRESHOLDS, "real": evaluate_real(real), "synthetic": {}}
    for p in synthetic_paths:
        blob = json.loads(Path(p).expanduser().read_text(encoding="utf-8"))
        for condition, rows in blob["synthetic"].items():
            if condition == "config":
                continue
            out["synthetic"][condition] = evaluate_synthetic(rows)
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--synthetic", nargs="*", default=[])
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    out = evaluate(a.real, a.synthetic)
    for name, block in out["real"].items():
        print(f"{name:5s} {block['verdict']}")
    for condition, block in out["synthetic"].items():
        print(f"{condition:16s} " + "  ".join(f"{k}={block[k]['verdict']}"
                                              for k in ("S1", "S2", "S3", "S4")))
    Path(a.output).expanduser().parent.mkdir(parents=True, exist_ok=True)
    Path(a.output).expanduser().write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("wrote", a.output)
    return out


if __name__ == "__main__":
    main()
