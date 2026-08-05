"""Aggregate the multi-partition inductive-channel runs and report the SPREAD.

Predeclared in NOTEBOOK_ENTRIES/PREDECLARED_inductive_channel_split_stability_20260805T0015Z.md.

Defines no statistic. Every retention is
``nonlinear_adjustment.retention_of_excess(inductive_arm, matched_transductive_control)``, computed
inside each run and stored under ``derived``; every S1 / null median / excess / permutation p is
``channel_under_adjustment``'s. This file reads JSON and takes ranges.
"""
import glob
import json
from pathlib import Path

import numpy as np
from scipy import stats

P = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/split_stability"

DEGENERACY_BARS = {              # §5.2 of the predeclaration, checked per split
    "n_frequent_sites_min": 20,
    "op_design_cols_min": 45,
    "audit_corr_band": (0.65, 0.85),
    "audit_var_band": (0.50, 0.70),
    "audit_above99_max": 0,
}


def load():
    rows = []
    for f in sorted(glob.glob(f"{P}/*.json")):
        tag = Path(f).stem
        if "overlap" in tag or tag == "AGGREGATE":
            continue
        d = json.load(open(f))
        art = "_".join(tag.split("_")[:2])
        name = "_".join(tag.split("_")[2:])
        ch_i = d.get("channel::inductive_exposure")
        ch_t = d.get("channel::transductive_exposure")
        ch_n = d.get("channel::none_exposure")
        if not (ch_i and ch_t and "derived" in d):
            rows.append({"artifact": art, "tag": name, "status": "INCOMPLETE"})
            continue
        cov = d["provenance"]["site_coverage"]
        audit = ch_i["adjustment_audit"]
        r = {
            "artifact": art, "tag": name,
            "split_seed": d["cohort"]["split_seed"],
            "f": d["config"]["discovery_fraction"],
            "n_discovery": d["cohort"]["n_discovery"], "n_exposure": d["cohort"]["n_exposure"],
            "exposure_digest": d["cohort"]["exposure_patient_digest"][:10],
            "n_frequent_sites": cov["n_frequent_sites_in_discovery_fold"],
            "coverage": cov["fraction_exposure_rows_site_adjustable"],
            "n_covered_rows": cov["n_exposure_rows_site_adjustable"],
            "op_design_cols": d["provenance"]["operator_x_n_design_columns"],
            "exp_design_cols": d["provenance"]["n_exposure_design_columns"],
            "retention": d["derived"]["retention_inductive_vs_matched_transductive"],
            "ind_S1": ch_i["observed_top_cca"], "ind_null": ch_i["null_median"],
            "ind_excess": ch_i["excess_over_null_median"], "ind_p": ch_i["permutation_p"],
            "ind_S2": ch_i["heldout_top_cca"], "ind_rank": ch_i["effective_rank_x_adjusted"],
            "tr_S1": ch_t["observed_top_cca"], "tr_null": ch_t["null_median"],
            "tr_excess": ch_t["excess_over_null_median"], "tr_p": ch_t["permutation_p"],
            "tr_S2": ch_t["heldout_top_cca"], "tr_rank": ch_t["effective_rank_x_adjusted"],
            "audit_corr": audit["per_axis_corr_raw_adjusted_median"],
            "audit_var": audit["residual_variance_ratio_median"],
            "audit_above99": audit["n_axes_corr_above_0_99"],
            "null_ratio_ind_over_tr": ch_i["null_median"] / ch_t["null_median"],
            "agreement": d["derived"]["adjuster_agreement_inductive_vs_matched_transductive"],
        }
        if ch_n:
            r["none_retention"] = (ch_n["excess_over_null_median"]
                                   / ch_t["excess_over_null_median"])
            r["none_S1"] = ch_n["observed_top_cca"]
        g_i, g_t = ch_i.get("global_pairing_null"), ch_t.get("global_pairing_null")
        if g_i and g_t and abs(g_t["excess_over_null_median"]) > 1e-12:
            r["retention_global"] = (g_i["excess_over_null_median"]
                                     / g_t["excess_over_null_median"])
            r["ind_global_null"] = g_i["null_median"]
            r["tr_global_null"] = g_t["null_median"]
        for enc in ("additive_design", "saturated_cell_design", "frozen_discovery_design"):
            for arm, short in (("inductive_exposure", "ind"), ("transductive_exposure", "tr")):
                k = f"ceiling_share_of_channel_excess::{enc}::{arm}"
                if k in d["derived"]:
                    r[f"ceil::{enc}::{short}"] = d["derived"][k]
                    r[f"ceilp::{enc}::{short}"] = \
                        d[f"ceiling::{enc}::{arm}"]["adjusted"]["permutation_p"]
        # §5.2 degeneracy guards, per split
        flags = []
        if r["n_frequent_sites"] < DEGENERACY_BARS["n_frequent_sites_min"]:
            flags.append("few_frequent_sites")
        if r["op_design_cols"] < DEGENERACY_BARS["op_design_cols_min"]:
            flags.append("narrow_operator_design")
        lo, hi = DEGENERACY_BARS["audit_corr_band"]
        if not lo <= r["audit_corr"] <= hi:
            flags.append("audit_corr_out_of_band")
        lo, hi = DEGENERACY_BARS["audit_var_band"]
        if not lo <= r["audit_var"] <= hi:
            flags.append("audit_var_out_of_band")
        if r["audit_above99"] > DEGENERACY_BARS["audit_above99_max"]:
            flags.append("axes_untouched")
        if r["null_ratio_ind_over_tr"] < 0.70:
            flags.append("inductive_null_collapsed")     # §5.3
        r["degeneracy_flags"] = flags
        r["status"] = "OK"
        rows.append(r)
    return rows


def main():
    rows = load()
    Path(f"{P}/AGGREGATE.json").write_text(json.dumps(rows, indent=1, default=float))
    incomplete = [r for r in rows if r["status"] != "OK"]
    if incomplete:
        print("INCOMPLETE:", [(r["artifact"], r["tag"]) for r in incomplete])

    for art in ("d2_h", "d2_i"):
        primary = [r for r in rows if r.get("status") == "OK" and r["artifact"] == art
                   and r["f"] == 0.5]
        if not primary:
            continue
        primary.sort(key=lambda r: r["retention"])
        ret = np.array([r["retention"] for r in primary])
        cov = np.array([r["coverage"] for r in primary])
        print(f"\n################ {art}  n_splits={len(ret)}  (f = 0.5) ################")
        print(f"retention: median {np.median(ret):.4f}  min {ret.min():.4f}  max {ret.max():.4f}  "
              f"range {ret.max()-ret.min():.4f}  sd {ret.std(ddof=1):.4f}")
        seed42 = [r for r in primary if r["split_seed"] == 42][0]
        rank42 = 1 + sorted(ret).index(seed42["retention"])
        print(f"seed 42 = {seed42['retention']:.4f}, rank {rank42} of {len(ret)} "
              f"(1 = lowest); max - seed42 = {ret.max()-seed42['retention']:+.4f}")
        if len(ret) > 2:
            rho, p = stats.spearmanr(cov, ret)
            print(f"spearman(retention, discovery-fold site coverage) = {rho:+.3f} (p = {p:.3f}, "
                  f"n = {len(ret)} -- essentially no power, reported for completeness)")
        print(f"coverage: min {cov.min():.4f}  max {cov.max():.4f}")
        hdr = ("seed", "ret", "retGlob", "cov", "sites", "opCol", "indS1", "indNull", "indExc",
               "trS1", "trNull", "trExc", "indS2", "trS2", "nullRatio", "p", "flags")
        print("  " + "".join(f"{h:>9}" for h in hdr[:-1]) + "   flags")
        for r in primary:
            print("  " + "".join(f"{v:>9}" for v in (
                r["split_seed"], f"{r['retention']:.4f}",
                f"{r.get('retention_global', float('nan')):.4f}", f"{r['coverage']:.4f}",
                r["n_frequent_sites"], r["op_design_cols"], f"{r['ind_S1']:.4f}",
                f"{r['ind_null']:.4f}", f"{r['ind_excess']:.4f}", f"{r['tr_S1']:.4f}",
                f"{r['tr_null']:.4f}", f"{r['tr_excess']:.4f}", f"{r['ind_S2']:.4f}",
                f"{r['tr_S2']:.4f}", f"{r['null_ratio_ind_over_tr']:.4f}",
                f"{r['ind_p']:.5f}")) + "   " + (",".join(r["degeneracy_flags"]) or "-"))

        # ceilings
        print(f"  --- labels-only ceiling, share of the channel's excess ({art}) ---")
        for enc in ("additive_design", "saturated_cell_design", "frozen_discovery_design"):
            for short in ("tr", "ind"):
                vals = [r[f"ceil::{enc}::{short}"] for r in primary
                        if f"ceil::{enc}::{short}" in r]
                ps = [r[f"ceilp::{enc}::{short}"] for r in primary
                      if f"ceilp::{enc}::{short}" in r]
                if not vals:
                    continue
                v = np.array(vals)
                print(f"    {enc:<24} {short:<4} median {np.median(v):+.4f}  "
                      f"[{v.min():+.4f}, {v.max():+.4f}]  n_sig(p<0.05) "
                      f"{sum(1 for x in ps if x < 0.05)}/{len(ps)}")

        grad = [r for r in rows if r.get("status") == "OK" and r["artifact"] == art
                and r["f"] != 0.5]
        if grad:
            print(f"  --- deliberate coverage gradient ({art}) ---")
            for r in sorted(grad + [seed42], key=lambda r: r["f"]):
                print(f"    f={r['f']}  n_disc={r['n_discovery']:>5}  cov={r['coverage']:.4f}  "
                      f"sites={r['n_frequent_sites']:>3}  opCol={r['op_design_cols']:>3}  "
                      f"ret={r['retention']:.4f}  indS1={r['ind_S1']:.4f}  "
                      f"indNull={r['ind_null']:.4f}  trS1={r['tr_S1']:.4f}  "
                      f"trNull={r['tr_null']:.4f}  p={r['ind_p']:.5f}")


if __name__ == "__main__":
    main()
