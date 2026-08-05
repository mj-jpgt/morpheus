"""Tests for the P2 figure suite, `v2/research/rebase/p2/figures/`.

Three things are checked.

**Every script imports and runs on SYNTHETIC input.** `_synthetic_corpus` builds a
complete, self-consistent stand-in for everything the nine display items read -
the vendored box artifacts under `figures/data/` and the repository markdown they
parse - with numbers drawn from a fixed RNG and every derived table (the variance
decomposition, the selection-rule marks, the subsample gap/sd, the exact binomial
p) computed from those numbers rather than copied. Each script is then run against
it. Because the scripts assert their sources against one another before drawing,
a script that carried a hard-coded figure value would fail on this corpus rather
than quietly draw the real number over synthetic data.

**Every script runs on the vendored data and writes vector + raster output.**
That is the regression test for the artwork actually shipped.

**A pending measurement is loud.** `--strict` must turn every pending placeholder
into a raised `PendingMeasurement`, so a release build cannot emit a figure with a
hole in it by accident.

Matplotlib only: the suite also asserts that no figure script imports seaborn.
"""
from __future__ import annotations

import importlib
import json
import math
import random
import sys
from pathlib import Path

import pytest

FIGURES = Path(__file__).resolve().parents[1] / "research" / "rebase" / "p2" / "figures"
MODULES = ["fig_f1_envelope", "fig_f2_variance", "fig_f3_floor", "fig_f4_defeaters",
           "fig_f5_verdict", "fig_f6_necessity", "fig_f7_dilution", "fig_f8_collapse",
           "fig_f9_decorrelation", "tab_t1_selection_rule"]

D2_ARMS = ["H", "I"]
D1_ARMS = ["P", "F"]
SEEDS = (42, 43, 44)
DILUTION_LEVELS = ["000", "010", "020", "030", "040", "060", "080"]
METRIC_KEYS = ["effective_rank_raw", "effective_rank_residualised", "rankme_raw",
               "rankme_residualised", "participation_ratio_raw",
               "participation_ratio_residualised", "stable_rank_raw",
               "stable_rank_residualised"]
SELECTION_ROWS = [
    ("effective_rank_raw", "effective_rank_raw"),
    ("effective_rank_residualised", "effective_rank_residualised"),
    ("rankme_raw (as published)", "rankme_raw"),
    ("rankme_residualised", "rankme_residualised"),
    ("participation_ratio_raw", "participation_ratio_raw"),
    ("participation_ratio_residualised", "participation_ratio_residualised"),
    ("stable_rank_raw", "stable_rank_raw"),
    ("stable_rank_residualised", "stable_rank_residualised"),
    ("alpha_ReQ_raw |alpha-1|", None),
    ("alpha_ReQ_resid |alpha-1|", None),
    ("LiDAR_raw", None),
    ("LiDAR_residualised", None),
]


@pytest.fixture(scope="module")
def figmod():
    if str(FIGURES) not in sys.path:
        sys.path.insert(0, str(FIGURES))
    import p2fig  # noqa: PLC0415
    return p2fig


# --------------------------------------------------------------------------
# Synthetic corpus
# --------------------------------------------------------------------------
def _decompose(values, *, log):
    v = {a: [math.log(x) if log else x for x in xs] for a, xs in values.items()}
    flat = [x for xs in v.values() for x in xs]
    grand = sum(flat) / len(flat)
    ss_arm = sum(len(xs) * (sum(xs) / len(xs) - grand) ** 2 for xs in v.values())
    ss_seed = sum((x - sum(xs) / len(xs)) ** 2 for xs in v.values() for x in xs)
    return ss_arm, ss_seed, 100.0 * ss_arm / (ss_arm + ss_seed), (ss_arm / 3) / (ss_seed / 8)


def _binom(k, n=6):
    pmf = [math.comb(n, i) / 2 ** n for i in range(n + 1)]
    return sum(p for p in pmf if p <= pmf[k] + 1e-12)


def _artifact(rng, label):
    """One synthetic artifact: metrics, points, subsample, all mutually consistent."""
    er_raw = rng.uniform(6.0, 34.0)
    er_res = er_raw * rng.uniform(1.02, 1.20)
    cca = rng.uniform(0.42, 0.65)
    metrics = {
        "effective_rank_raw": er_raw,
        "effective_rank_residualised": er_res,
        "rankme_raw": er_raw * rng.uniform(0.1, 0.5),
        "rankme_residualised": er_res * rng.uniform(0.99, 1.01),
        "participation_ratio_raw": er_raw * rng.uniform(0.1, 0.2),
        "participation_ratio_residualised": er_res * rng.uniform(0.1, 0.2),
        "stable_rank_raw": rng.uniform(1.8, 3.5),
        "stable_rank_residualised": rng.uniform(1.8, 3.5),
    }
    points = {name: {"n_targets": n, "top_cca": cca * rng.uniform(0.95, 1.05),
                     "null_mean": rng.uniform(0.138, 0.142)}
              for name, n in (("unrestricted", 90), ("untrained40", 40),
                              ("random_control", 90), ("hallmark_in_training", 50))}
    points["untrained40"]["top_cca"] = cca
    sub = {"effective_rank_residualised": {"mean": er_res * rng.uniform(0.995, 1.005),
                                           "sd": er_res * 0.004, "n": 40, "frac": 0.8}}
    return {"path": f"/synthetic/{label}.npz", "n_test": 2766,
            "metrics": metrics, "points": points, "subsample": sub}


def _synthetic_corpus(root: Path) -> dict:
    rng = random.Random(20260804)
    root.mkdir(parents=True, exist_ok=True)

    d2 = {f"{a}{s}": _artifact(rng, f"{a}{s}") for a in D2_ARMS for s in SEEDS}
    d1 = {f"{a}{s}": _artifact(rng, f"{a}{s}") for a in D1_ARMS for s in SEEDS}
    # The information direction is a fact of the design, not of the draw: the
    # first-named arm carries the larger channel on every seed. F1(c) and F6
    # assert that direction, so the corpus has to honour it while leaving every
    # RANK value random - which is exactly the asymmetry the paper is about.
    for src_dict, hi, lo in ((d2, "H", "I"), (d1, "P", "F")):
        for s in SEEDS:
            a = src_dict[f"{hi}{s}"]["points"]["untrained40"]["top_cca"]
            b = src_dict[f"{lo}{s}"]["points"]["untrained40"]["top_cca"]
            if a <= b:
                (src_dict[f"{hi}{s}"]["points"]["untrained40"]["top_cca"],
                 src_dict[f"{lo}{s}"]["points"]["untrained40"]["top_cca"]) = (b, a)
    # D2 seed 44 is the pair T1 marks unresolvable and F4(a) annotates: its two
    # arms must sit inside the estimator's own sampling noise, or the figures'
    # own guards are testing nothing.
    for key in ("effective_rank_raw", "effective_rank_residualised"):
        d2["I44"]["metrics"][key] = d2["H44"]["metrics"][key] * 1.004
    d2["I44"]["subsample"]["effective_rank_residualised"]["mean"] = (
        d2["H44"]["subsample"]["effective_rank_residualised"]["mean"] * 1.0004)
    d2["_config"] = d1["_config"] = {"seed": 42, "subsamples": 40, "subsample_fraction": 0.8}

    def write(rel, payload):
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(payload if isinstance(payload, str)
                     else json.dumps(payload, indent=1), encoding="utf-8")

    write("ws_p2/out/P2_METRICS_D2.json", d2)
    write("ws_p2/out/P2_METRICS_D1.json", d1)

    # ---- rank variants + the selection marks they must agree with ---------
    variants, verdicts = {}, {}
    labels = [f"{a}{s}" for a in D2_ARMS for s in SEEDS] + [f"{a}{s}" for a in D1_ARMS for s in SEEDS]
    src = {**{k: v for k, v in d2.items() if k != "_config"},
           **{k: v for k, v in d1.items() if k != "_config"}}
    for label in labels:
        er = src[label]["metrics"]["effective_rank_residualised"]
        variants[label] = {"R1": er, "R2": er * 0.63, "R3": er * 0.66,
                           "PR": er * 0.18, "PR_rownorm": er * 0.21,
                           "cca40": src[label]["points"]["untrained40"]["top_cca"]}
    pairs = [("D2 s42", "H42", "I42"), ("D2 s43", "H43", "I43"), ("D2 s44", "H44", "I44"),
             ("D1 s42", "P42", "F42"), ("D1 s43", "P43", "F43"), ("D1 s44", "P44", "F44")]
    for stat in ("R1", "R2", "R3", "PR", "PR_rownorm"):
        marks = ["OK" if ((variants[a][stat] > variants[b][stat])
                          == (variants[a]["cca40"] > variants[b]["cca40"])) else "MISS"
                 for _lab, a, b in pairs]
        verdicts[stat] = {"marks": marks,
                          "D2": sum(m == "OK" for m in marks[:3]),
                          "D1": sum(m == "OK" for m in marks[3:]),
                          "ALL": sum(m == "OK" for m in marks)}
    write("ws_p2/out/P2_RANK_VARIANTS.json", {"artifacts": variants, "verdicts": verdicts})

    write("ws_p2/out/P2_ROBUSTNESS.json", {
        label: {view: {"eff_rank": variants[label]["R1"] * f,
                       "cca_insample": variants[label]["cca40"] * f,
                       "cca_heldout": variants[label]["cca40"] * f * 0.96}
                for view, f in (("wsi_biology", 1.0), ("rna_biology", 1.3), ("full_biology", 1.4))}
        for label in labels})

    # ---- the run log, every table derived from the JSON above -------------
    dec = []
    for name, scale, getter in (
            ("effective rank (residualised)", "log",
             lambda x: x["metrics"]["effective_rank_residualised"]),
            ("RankMe (raw, as published)", "log", lambda x: x["metrics"]["rankme_raw"]),
            ("GROUND TRUTH: untrained40 top-CCA", "raw",
             lambda x: x["points"]["untrained40"]["top_cca"])):
        by_arm = {a: [getter(src[f"{a}{s}"]) for s in SEEDS] for a in ("H", "I", "P", "F")}
        dec.append((name, scale) + _decompose(by_arm, log=(scale == "log")))

    lines = ["===== p2_selection_rule =====", "SELECTION-RULE SCORE"]
    lines.append("metric" + " " * 34 + "".join(f"{p[0]:>10}" for p in pairs)
                 + f"{'D2':>6}{'D1':>6}{'ALL':>6}{'p':>8}")
    for row_label, key in SELECTION_ROWS:
        if key is None:
            marks = verdicts["R1"]["marks"]
        else:
            marks = ["OK" if ((src[a]["metrics"][key] > src[b]["metrics"][key])
                              == (variants[a]["cca40"] > variants[b]["cca40"])) else "MISS"
                     for _lab, a, b in pairs]
        allk = sum(m == "OK" for m in marks)
        lines.append(f"{row_label:<40}" + "".join(f"{m:>10}" for m in marks)
                     + f"{sum(m == 'OK' for m in marks[:3])}/3".rjust(6)
                     + f"{sum(m == 'OK' for m in marks[3:])}/3".rjust(6)
                     + f"{allk}/6".rjust(6) + f"{_binom(allk):>8.3f}")
    lines.append("WITHIN-ARM SEED VARIABILITY")
    lines.append("  --- effective_rank_residualised ---")
    for lab, a, b in pairs:
        sa = src[a]["subsample"]["effective_rank_residualised"]
        sb = src[b]["subsample"]["effective_rank_residualised"]
        gap = abs(sa["mean"] - sb["mean"])
        lines.append(f"    {lab}: {a}={sa['mean']:>8.3f}+-{sa['sd']:>6.3f}  "
                     f"{b}={sb['mean']:>8.3f}+-{sb['sd']:>6.3f}  gap={gap:>7.3f}  "
                     f"gap/sd={gap / math.hypot(sa['sd'], sb['sd']):>6.2f}")
    lines += ["===== p2_necessity_and_variance =====", "(A) VARIANCE DECOMPOSITION",
              "quantity                                   scale      SS_arm     SS_seed"
              "   arm share    F(3,8)",
              "-" * 104]
    for name, scale, ss_a, ss_s, share, f in dec:
        lines.append(f"{name:<40}{scale:>10}{ss_a:>12.4f}{ss_s:>12.4f}{share:>11.1f}%{f:>10.2f}")
    lines += ["", "  2 violating pair(s) of 66 examined:", ""]
    for lo, hi, scope in (("P44", "I43", "cross-experiment"),
                          ("H44", "I43", "within-experiment, across arms")):
        lines.append(f"{lo:>24}{variants[lo]['R1']:>9.3f}{variants[lo]['cca40']:>9.4f}   |"
                     f"{hi:>22}{variants[hi]['R1']:>9.3f}{variants[hi]['cca40']:>9.4f}"
                     f"   {abs(variants[hi]['R1'] / variants[lo]['R1']):>6.2f}x"
                     f"{abs(variants[lo]['cca40'] - variants[hi]['cca40']):>+9.4f}  {scope}")
    lines += ["===== p2_rank_variants ====="]
    for stat, entry in verdicts.items():
        lines.append(f"{stat:>12}" + "".join(f"{m:>9}" for m in entry["marks"]))
    write("ws_p2/out/p2_run.log", "\n".join(lines) + "\n")

    # ---- canonical rank recomputation -------------------------------------
    def variant_block(er):
        r2 = er * 0.63
        return {"R1": er, "R2": r2, "R3": r2, "R1_uncentred": er * 0.99,
                "R1_rownorm": er, "R2_uncentred": r2 * 0.9}
    inst6 = {}
    for label in labels[:6]:
        raw = src[label]["metrics"]["effective_rank_raw"]
        res = src[label]["metrics"]["effective_rank_residualised"]
        inst6[label] = {"n_rows": 2766, "n_features": 256,
                        "raw": variant_block(raw),
                        "residualised": {**variant_block(res), "R3": res * 0.66}}
    inst6["H42_recovered_reexport"] = {
        "raw": variant_block(7.9), "residualised": {**variant_block(8.68), "R3": 6.0}}
    inst5 = {}
    for a, key in (("p", "P"), ("f", "F")):
        for s in SEEDS:
            label = f"{key}{s}"
            inst5[f"d1_{a}_seed{s}"] = {
                "raw": variant_block(src[label]["metrics"]["effective_rank_raw"]),
                "residualised": {
                    **variant_block(src[label]["metrics"]["effective_rank_residualised"]),
                    "R3": src[label]["metrics"]["effective_rank_residualised"] * 0.66}}
    dil_raw = [196.2 - 5.0 * i for i in range(7)]
    dil_res = [210.2 - 1.1 * i for i in range(7)]
    inst4 = {f"wsi_dilution_{lvl}": {"raw": variant_block(dil_raw[i]),
                                     "residualised": variant_block(dil_res[i])}
             for i, lvl in enumerate(DILUTION_LEVELS)}
    tripwire = {}
    for a, arm in (("p", 110.0), ("f", 8.0)):
        for i, s in enumerate(SEEDS):
            tripwire[f"d1_{a}_seed{s}"] = {"epoch": 11, "R3_at_step_200": arm + 3.0 * i,
                                           "source": f"/synthetic/d1_{a}_seed{s}.jsonl"}
    write("ws_rank/RANK_RECOMPUTE.json", {
        "canonical": "centred|order1",
        "variant_order": ["R1", "R2", "R3", "R1_uncentred", "R1_rownorm", "R2_uncentred"],
        "instances": {"instance6_D2": inst6, "instance4_dilution": inst4,
                      "instance5_D1B": inst5,
                      "instability_tripwire_step200_R3": tripwire}})
    write("ws_rank/RANK_RECOMPUTE_P1B.json", {
        "v21_release_20260720_retry3_resume_safe": {
            "diagnostic_full_seed42": {"wsi_biology": {"raw": variant_block(38.5),
                                                       "residualised": variant_block(43.4)}},
            "diagnostic_programme_only_seed42": {"wsi_biology": {"raw": variant_block(32.1),
                                                                 "residualised": variant_block(35.4)}}}})

    # ---- box readouts and logs --------------------------------------------
    write("e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json",
          {k: {"path": v["path"], "n_test": v["n_test"],
               "effective_rank_raw": v["metrics"]["effective_rank_raw"],
               "effective_rank_residualised": v["metrics"]["effective_rank_residualised"],
               "points": v["points"]}
           for k, v in d2.items() if k != "_config"})
    write("e0_run/d2_v3/RECOVERED_SEED42_READOUT.json", {
        "H42_recovered": {"path": "/synthetic/recovered.npz", "n_test": 2766,
                          "effective_rank_raw": 7.9, "effective_rank_residualised": 8.6809,
                          "points": {"unrestricted": {"n_targets": 90, "top_cca": 0.5861222}}}})

    def bootstrap(n_targets, groups, scale):
        return {"experiment": "D1", "arm_a": "programme_only", "arm_b": "programme_free",
                "target_groups": groups, "n_targets": n_targets,
                "pairs": [{"programme_only": f"/synthetic/d1_p_seed{s}.npz",
                           "programme_free": f"/synthetic/d1_f_seed{s}.npz",
                           "n_test": 2766,
                           "point_programme_only": d1[f"P{s}"]["points"]["untrained40"]["top_cca"],
                           "point_programme_free": d1[f"F{s}"]["points"]["untrained40"]["top_cca"],
                           "programme_free_minus_programme_only": {
                               mode: {"n_valid": 2000, "repeats": 2000,
                                      "point_delta": scale * (0.05 + 0.01 * i),
                                      "mean_delta": scale * (0.05 + 0.01 * i),
                                      "ci95_low": scale * (0.09 + 0.01 * i),
                                      "ci95_high": scale * (0.02 + 0.01 * i), "mode": mode,
                                      "p_improve": 0.0}
                               for mode in ("patient", "cancer")}}
                          for i, s in enumerate(SEEDS)]}
    write("e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json",
          bootstrap(40, ["heldout_pathway", "immune_tme", "tumour_state"], -1.0))
    write("e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json",
          bootstrap(90, ["random_control"], -0.3))
    write("e0_run/d1_v2/D1_PAIR_MANIFEST.json", {
        "objective_only_difference": True,
        "preregistered_prediction": "synthetic preregistered prediction -- escalate"})

    write("e0_run/collapse_diag.log",
          "--- A. synthetic arm\n"
          "  in-batch InfoNCE   3.4681 -> 2.7734   (chance 2.7726, perfect ~0)\n"
          "  retrieval acc@1    0.062 -> 0.000          (chance 0.062)\n"
          "  cross pos cos      0.0538 -> 0.9959\n"
          "  cross neg cos      0.0816 -> 0.9960\n"
          "  WSI within-modality offdiag cos 0.7089 -> 0.9999\n"
          "  z_biology matrix rank 16 -> 16  (max 16)\n\n--- B. other\n")
    write("e0_run/d1_diag/diag_d.log",
          "cosine margin required for in-batch InfoNCE <= 0.10 with 15 negatives: 0.3472\n"
          + "".join(f"  {st:>5} loss 2.4 acc 0.125 pos {p} worst-neg {p} min-margin {m} "
                    f"wsi-wsi 0.61 eff-rank {r} std 0.03\n"
                    for st, p, m, r in ((0, "0.1235", "-0.2190", "12.88"),
                                        (25, "0.9344", "-0.1103", "1.21"),
                                        (50, "0.9993", "-0.0001", "1.00"),
                                        (100, "0.9995", "-0.0001", "1.00"),
                                        (200, "0.9996", "-0.0001", "1.00"),
                                        (400, "0.9999", "-0.0000", "1.00"))))
    for m, vals in (("0.999", (7.15, 6.92, 7.25)), ("0", (1.80, 1.46, 1.98))):
        for rep, v in zip((1, 2, 3), vals):
            write(f"e0_run/d1_diag/probevar_m{m}_{rep}.log",
                  f"momentum={m}  decorrelation=0.04  capacity=4096  steps=200\n"
                  f"DONE momentum={m} decorr=0.04 final_eff_rank={v:.2f} rna_rna=0.89\n")

    # F9's decorrelation ablation, in the two-rank-statistic log format
    # `d1_momentum_probe.py` writes. The three arms MUST share their step-0 row,
    # because fig_f9 refuses to draw a sweep that varies more than the weight -
    # and the mseed log must share every header field except the step budget,
    # because that is what makes it a same-seed repeat rather than another arm.
    def probe_log(momentum, decorr, steps, series, rna):
        head = (f"momentum={momentum}  decorrelation={decorr}  capacity=4096  lr=0.0002  "
                f"seed=42  steps={steps}  probe=256\n"
                "tau=1000.0 steps   T=capacity/214=19.1 steps   tau/T=52.25\n"
                "  step   R3-rank  CANONICAL  feat-std  rna-rna  contrastive\n")
        rows = ""
        for i, step in enumerate(range(0, steps + 1, 50)):
            r3 = 67.55 if step == 0 else series * (1.0 + 0.03 * (i % 3))
            rows += (f"{step:>6}{r3:>10.2f}{r3 * 1.52:>11.2f}{0.03:>10.4f}"
                     f"{(0.3650 if step == 0 else rna):>9.4f}"
                     f"{'nan' if step == 0 else '4.5000':>13}\n")
            rows += f"        STALENESS age:cos  0:0.1{i}0\n"
        return head + rows + (f"\nDONE momentum={momentum} decorr={decorr} "
                              f"final_eff_rank={series:.2f} rna_rna={rna:.4f}\n")

    for decorr, r3, rna in (("0.0", 4.32, 0.4774), ("0.01", 6.22, 0.7657), ("0.04", 8.01, 0.8696)):
        write(f"e0_run/d1_diag/ablate_decorr{decorr}.log", probe_log("0.999", decorr, 400, r3, rna))
    write("e0_run/d1_diag/mseed_m0.999_s42.log", probe_log("0.999", "0.04", 500, 7.46, 0.8300))

    # ---- the probe block's own retraining floor ----------------------------
    # F9 draws THIS floor and not draft 4.1's exported-block 3.295x, because the
    # runs it plots are on the probe. The floor is per statistic and per reading
    # step by construction, so the synthetic blob is shaped that way and the two
    # statistics are given DIFFERENT folds - a script that read one and drew it on
    # both panels would pass a corpus where they agreed.
    def _floor_cell(lo, fold):
        return {"floor": fold, "floor_arm": "m0", "floor_min": lo, "floor_max": lo * fold,
                "why": None,
                "arm_floor": {arm: {"fold": fold if arm == "m0" else 1.1,
                                    "min": lo, "max": lo * fold,
                                    "shape": {"bimodal": False, "defined": True}}
                              for arm in ("m0999", "m0")}}

    _probe_stats = {"R1": (2.0, 1.5702), "R3": (1.5, 1.4489), "R2": (1.5, 1.44),
                    "PR": (1.2, 1.30), "PR_rownorm": (1.2, 1.30), "RankMe": (2.0, 1.30),
                    "stable_rank": (1.4, 1.34), "alpha_req": (2.0, 1.20),
                    "alpha_req_abs_dev": (1.0, 1.49), "hard_rank": (256.0, 1.0)}
    write("e0_run/d1_probefloor/out/P2_PROBE_FLOORS.json",
          {"_": "synthetic", "absent": [], "arms": {},
           "config": {"arms": ["m0999", "m0"], "repeats_per_arm": 5, "seed": 42,
                      "lr": 0.0002, "capacity": 4096, "decorrelation": 0.04,
                      "objective_profile": "programme_free",
                      "block": "fixed held-out probe, raw"},
           "floor_rule": "max/min over five same-seed repeats, pooled across the two arms",
           "shape_rule": "synthetic", "reps": {},
           "floors": {str(step): {view: {stat: _floor_cell(lo, fold)
                                         for stat, (lo, fold) in _probe_stats.items()}
                                  for view in ("wsi_biology", "rna_biology", "wsi_rna_paired")}
                      for step in (0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)}})

    # ---- extracted ---------------------------------------------------------
    runs = {}
    for a, arm, base in (("p", "programme_only", 110.0), ("f", "programme_free", 8.0)):
        for i, s in enumerate((42, 43, 44, 45, 46)):
            runs[f"d1_{a}_seed{s}"] = {
                "arm": arm, "seed": s, "value": base + 3.0 * i, "epoch": 11,
                "source": f"/synthetic/d1_{a}_seed{s}.jsonl",
                "source_sha256": "0" * 64, "source_bytes": 1}
    write("extracted/F3_TRIPWIRE_STEP200_R3_n5.json",
          {"statistic": "R3", "key": "train_rank_tripwire_observed", "epoch": 11,
           "global_step": 200, "runs": runs})
    # The controlled retraining repeat, as it looks once it has reported. The
    # shape is the real one: four repeats close together and one at a third of
    # them, because F1(a) is drawn around that bimodality and a smooth spread
    # would let a panel that averaged it away still pass.
    env_rank_res = [28.320, 8.834, 28.348, 29.106, 28.959]
    env_rank_raw = [24.481, 8.033, 24.504, 24.990, 24.912]
    env_chan = [0.6182, 0.5859, 0.6123, 0.6110, 0.6098]
    write("extracted/F1_RETRAINING_REPEAT.json",
          {"root": "/synthetic/d1_envelope", "expected": "rep{1..5}/", "complete": True,
           "predeclaration": "NOTEBOOK_ENTRIES/SYNTHETIC.md",
           "readout_log": "/synthetic/d1_envelope_readout.log",
           "readout_module": "v2/research/rebase/d1_envelope_readout.py",
           "statistic": "R1",
           "printed_spread": {
               key: {"min": min(vals), "max": max(vals),
                     "spread": round(max(vals) / min(vals), 3)}
               for key, vals in (("rank_residualised", env_rank_res),
                                 ("rank_raw", env_rank_raw), ("channel", env_chan))},
           "reps": {f"rep{i}": {"epochs_logged": 40, "last_epoch": 39,
                                "artifact_npz": f"/synthetic/d1_envelope/rep{i}.npz",
                                "rank_raw": env_rank_raw[i - 1],
                                "rank_residualised": env_rank_res[i - 1],
                                "channel": env_chan[i - 1]}
                    for i in range(1, 6)}})
    write("MANIFEST.json", {"fetched_utc": "synthetic", "box": "synthetic", "files": {}})
    return {"d2": d2, "d1": d1, "variants": variants, "dil_raw": dil_raw}


def _synthetic_repo(corpus) -> dict[str, str]:
    """Synthetic stand-ins for the repository markdown the figures parse."""
    d2, variants, dil_raw = corpus["d2"], corpus["variants"], corpus["dil_raw"]
    strat = "\n".join(
        f"| {s} | {d2[f'H{s}']['points']['untrained40']['top_cca']:.4f} | "
        f"{d2[f'I{s}']['points']['untrained40']['top_cca']:.4f} | "
        f"{d2[f'I{s}']['points']['untrained40']['top_cca'] - d2[f'H{s}']['points']['untrained40']['top_cca']:+.4f} "
        f"| [-0.16, -0.09] | 0.0000 | [-0.18, -0.06] | 0.0010 |" for s in SEEDS)
    d2_md = (
        "# D2\n\n## 2. THE STRATIFIED READOUT\n\n"
        "| seed | Hallmark | PBS | d | patient CI | p | cancer CI | p |\n"
        "|---|---|---|---|---|---|---|---|\n" + strat + "\n\n"
        "## 4. Determinism\n\n"
        "| | unrestricted top-CCA | effective rank |\n|---|---|---|\n"
        "| recorded in the lost run, H seed 42 | 0.5861 | - |\n"
        "| re-export | 0.58612 | 8.68 |\n")

    dilution_md = (
        "# Dilution\n\n## 2. The sweep\n\n"
        "| requested d | achieved d | adjusted | held-out | ratio | null-corrected | floor | "
        "attenuation | effective rank | perm p |\n"
        "|---|---|---|---|---|---|---|---|---|---|\n"
        + "\n".join(
            f"| {r:.2f} | {a:.3f} | 0.55 | {0.49 - 0.04 * i:.4f} | {1.0 - 0.07 * i:.3f} | "
            f"{1.0 - 0.10 * i:.3f} | 0.40 | {1.0 + 0.01 * i:.3f} | {dil_raw[i]:.1f} | 0.0033 |"
            for i, (r, a) in enumerate(zip((0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8),
                                           (0.0, 0.091, 0.211, 0.302, 0.4, 0.6, 0.8))))
        + "\n\nThe permutation null median is 0.145-0.147 at every level.\n")

    rank_entry = (
        "# rank\n\nthe maximum relative difference between the old absolute cut and\n"
        "the new relative one, over all 68 recomputed\nartifact x block combinations, "
        "is **0.000e+00**.\n\n"
        "| ratio | min | median | max |\n|---|---|---|---|\n"
        "| R2 / R1 | 0.338 | 0.629 | 0.813 |\n"
        "| R3 / R1 | 0.351 | 0.655 | 0.826 |\n"
        "| R1 uncentred / R1 | 0.116 | 0.995 | 1.000 |\n")

    criterion = ("> a pair (lo, hi) is a violation iff `eff_rank(hi)/eff_rank(lo) >= 2.0` **and**\n"
                 "> `CCA(lo) - CCA(hi) >= 0.0705`.\n")
    script = "RANK_FOLD = 2.0\nCCA_DELTA = 0.0705\n"

    d1b = ("verified: effective rank 67.55 at step 0\n\n"
           "| `decorrelation` | `consistency` | step 50 | 100 | 150 | 200 | 250 |\n"
           "|---|---|---|---|---|---|---|\n"
           "| 0.04 | 1.0 | 4.08 | 1.95 | 2.16 | 1.68 | 1.59 |\n"
           "| 0.0 | 1.0 | 2.62 | 2.16 | 2.47 | 1.94 | 2.17 |\n"
           "| 0.04 | 0.1 | 2.99 | 3.43 | - | - | - |\n"
           "| 0.0 | 0.1 | 2.97 | 2.00 | 2.50 | - | - |\n"
           "| 0.0 | 0.0 | 2.98 | 1.98 | 1.86 | - | - |\n\n"
           "at step 200 the RNA-view mutual cosine is 0.9813 with versus 0.4160 without.\n")
    d1_mid = ("| arm | epoch | rank | hard | cos | std |\n|---|---|---|---|---|---|\n"
              "| `programme_free` seed 42 | 23 | 1.76 | 9 | 0.977 | 0.0137 |\n")
    d1a_end = ("| arm | epoch | rank | hard | cos | std |\n|---|---|---|---|---|---|\n"
               "| `programme_free` seed 42 | 39 | 1.71 | 11 | 0.986 | 0.0156 |\n")

    # T1's block-dependence band. The gene-set block MUST reproduce the metrics
    # JSON's own untrained-40 contrast (T1 asserts exactly that), and one block
    # flips seed 43, which is the shape of the real finding.
    exam_panel = []
    for key in ("geneset_untrained40", "pbs_codes128_ARM_I_OWN", "pca_basis128",
                "pbs_shuffled_s1", "geneset_random_control", "random_dictionary128"):
        for seed in SEEDS:
            h = d2[f"H{seed}"]["points"]["untrained40"]["top_cca"]
            i = d2[f"I{seed}"]["points"]["untrained40"]["top_cca"]
            delta = i - h
            if seed == 43 and key in ("pbs_codes128_ARM_I_OWN", "pca_basis128"):
                delta = -delta
            exam_panel.append({"exam": key, "seed": seed, "block": "residualised",
                               "n_targets": 40, "point_hallmark": h, "point_pbs": i,
                               "pbs_minus_hallmark": delta})

    return {
        "v2/research/rebase/nature/d2_coordinate_system/out/EXAM_PANEL.json":
            json.dumps(exam_panel),
        "v2/research/rebase/nature/D2_RESULT.md": d2_md,
        "v2/research/rebase/nature/DILUTION_LOWER_BOUND.md": dilution_md,
        "NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed"
        "_20260804T0005Z.md": rank_entry,
        "NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md": criterion,
        "v2/research/rebase/p2/p2_necessity_and_variance.py": script,
        "NOTEBOOK_ENTRIES/d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md": d1b,
        "NOTEBOOK_ENTRIES/d1_programme_free_collapsing_in_training_20260803T1930Z.md": d1_mid,
        "NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner"
        "_20260804T0100Z.md": d1a_end,
    }


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", MODULES)
def test_module_imports(figmod, name):
    module = importlib.import_module(name)
    assert callable(module.main), name
    assert module.__doc__ and len(module.__doc__) > 200, f"{name} has no provenance docstring"


def test_no_seaborn():
    """Matplotlib only. The docstrings may SAY seaborn; nothing may import it."""
    for path in sorted(FIGURES.glob("*.py")):
        for line in path.read_text(encoding="utf-8").splitlines():
            code = line.split("#")[0].strip()
            assert not code.startswith(("import seaborn", "from seaborn")), (
                f"{path.name}: {line.strip()!r}; the brief is matplotlib only")


def test_axis_label_requires_a_statistic_and_a_block(figmod):
    assert "R1" in figmod.axis_label("R1_short", "resid")
    assert "residualised" in figmod.axis_label("R1_short", "resid")
    with pytest.raises(KeyError):
        figmod.axis_label("effective_rank", "resid")
    with pytest.raises(KeyError):
        figmod.axis_label("R1_short", "whatever")


def test_minus_sign_is_normalised(figmod):
    """A real minus sign must parse as negative; this guard protects F1(c)'s sign."""
    assert figmod.nums("−0.1325") == [-0.1325]
    assert figmod.ci("[−0.1605, −0.0993]") == (-0.1605, -0.0993)


def test_md_section_scopes_the_parse(figmod):
    text = "# A\n\n| x | 1 |\n\n## B\n\n| x | 2 |\n\n## C\n\n| x | 3 |\n"
    assert figmod.md_rows(figmod.md_section(text, "B"), n_cols=2) == [["x", "2"]]


@pytest.mark.parametrize("name", MODULES)
def test_runs_on_synthetic_input(figmod, tmp_path, monkeypatch, name):
    """Every script must render from a synthetic corpus and nothing else."""
    corpus = _synthetic_corpus(tmp_path / "data")
    repo = _synthetic_repo(corpus)
    monkeypatch.setattr(figmod, "DATA", tmp_path / "data")
    monkeypatch.setattr(figmod, "OUT", tmp_path / "out")
    monkeypatch.setattr(figmod, "repo_text",
                        lambda rel: repo.get(rel) or pytest.fail(f"unstubbed repo read: {rel}"))
    monkeypatch.setattr(sys, "argv", [name])

    module = importlib.import_module(name)
    importlib.reload(module)
    monkeypatch.setattr(module.P, "DATA", tmp_path / "data")
    monkeypatch.setattr(module.P, "OUT", tmp_path / "out")
    monkeypatch.setattr(module.P, "repo_text",
                        lambda rel: repo.get(rel) or pytest.fail(f"unstubbed repo read: {rel}"))
    assert module.main() == 0
    written = sorted(p.suffix for p in (tmp_path / "out").glob("*"))
    assert written == [".pdf", ".png", ".svg"], written


def test_make_all_writes_every_item(figmod, tmp_path, monkeypatch):
    """The shipped artwork: nine items, three formats each, from the vendored data."""
    monkeypatch.setattr(figmod, "OUT", tmp_path / "out")
    monkeypatch.setattr(sys, "argv", ["make_all.py"])
    import make_all  # noqa: PLC0415
    importlib.reload(make_all)
    for name in MODULES:
        module = importlib.import_module(name)
        importlib.reload(module)
        monkeypatch.setattr(module.P, "OUT", tmp_path / "out")
    assert make_all.main() == 0
    for ext in (".pdf", ".svg", ".png"):
        assert len(list((tmp_path / "out").glob(f"*{ext}"))) == len(MODULES), ext
    for path in (tmp_path / "out").glob("*.pdf"):
        assert path.stat().st_size > 5_000, path


def test_pending_measurement_fails_loudly(figmod):
    """--strict must raise rather than draw a placeholder.

    Every display item is currently PLOTTABLE -- F1(d)'s retraining repeat landed on
    2026-08-04 and was the last pending panel -- so this asserts the mechanism itself
    rather than a particular figure. The next pending measurement to be added must
    still fail a release build, and that is what is checked here.
    """
    from matplotlib import pyplot as plt  # noqa: PLC0415

    fig, ax = plt.subplots()
    try:
        with pytest.raises(figmod.PendingMeasurement):
            figmod.pending_panel(ax, title="(x) synthetic", marker="[PENDING]",
                                 blocked_on="a measurement", arrives_at="/nowhere",
                                 strict=True)
    finally:
        plt.close(fig)


def test_no_display_item_is_pending(figmod, tmp_path, monkeypatch):
    """A release build must emit nine complete items with no hole in any of them."""
    monkeypatch.setattr(figmod, "OUT", tmp_path / "out")
    monkeypatch.setattr(sys, "argv", ["make_all.py", "--strict"])
    import make_all  # noqa: PLC0415
    importlib.reload(make_all)
    for name in MODULES:
        module = importlib.import_module(name)
        importlib.reload(module)
        monkeypatch.setattr(module.P, "OUT", tmp_path / "out")
    assert make_all.main() == 0, figmod.PENDING_LOG


def test_vendored_data_matches_its_manifest(figmod):
    """Every vendored byte is the byte the manifest recorded from the box."""
    manifest = figmod.manifest()
    assert manifest["files"], "the data manifest is empty"
    for rel, entry in manifest["files"].items():
        assert figmod.sha256(rel) == entry["sha256"], (
            f"{rel} does not match data/MANIFEST.json; re-run extract_from_box.py")
