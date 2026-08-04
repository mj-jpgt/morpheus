"""Score every label-free metric AS A SELECTION RULE on the matched arm pairs.

VENDORED 2026-08-05 from `~/e0_run/p2_selection_rule.py`, unchanged apart from this paragraph.
Produces `paper/P2_RANK_DRAFT.md` §4.6's table and §4.4(1)'s subsampling table. Reads only the
JSON written by `p2_competing_metrics.py`, so it recomputes nothing and cannot drift from it.

The question the paper turns on is not "does the metric correlate with information" but
"given two matched arms, does it pick the one that actually carries more molecular
information". Correlation coefficients over pooled artifacts confound the within-arm and
between-arm variation; the paired comparison does not.

Ground truth = top canonical correlation with the 40 targets neither arm trained on
(heldout_pathway + immune_tme + tumour_state), residualised, from the same code path as
d2_readout.py.

D2 pairs: H (Hallmark) vs I (PBS). D1 pairs: P (programme_only) vs F (programme_free).
The "correct" arm is whichever has the larger untrained40 top-CCA -- established
independently of every metric under test.
"""
from __future__ import annotations

import json
import sys
from math import comb

import numpy as np

# metric -> (accessor, direction) ; direction=+1 means "larger is better quality"
#   alpha-ReQ is scored as PROXIMITY TO 1 (|alpha - 1|, smaller is better). NOTE: the paper does
#   NOT state a |alpha - 1| rule; it states a "Goldilocks zone" -- "the best representations are
#   those where alpha is in a range that is close to 1" and "either too high or too low an alpha
#   value implies poor generalization". |alpha - 1| is OUR operationalisation of that, and the
#   paper's own Algorithm 1 uses an adaptive interval instead. Flagged in the notebook entry.
METRICS = [
    ("effective_rank_raw", lambda m: m["effective_rank_raw"], +1),
    ("effective_rank_residualised", lambda m: m["effective_rank_residualised"], +1),
    ("rankme_raw (as published)", lambda m: m["rankme_raw"], +1),
    ("rankme_residualised", lambda m: m["rankme_residualised"], +1),
    ("participation_ratio_raw", lambda m: m["participation_ratio_raw"], +1),
    ("participation_ratio_residualised", lambda m: m["participation_ratio_residualised"], +1),
    ("stable_rank_raw", lambda m: m["stable_rank_raw"], +1),
    ("stable_rank_residualised", lambda m: m["stable_rank_residualised"], +1),
    ("alpha_ReQ_raw |alpha-1|", lambda m: abs(m["alpha_req_raw"]["alpha"] - 1.0), -1),
    ("alpha_ReQ_resid |alpha-1|", lambda m: abs(m["alpha_req_residualised"]["alpha"] - 1.0), -1),
    ("LiDAR_raw", lambda m: m["lidar_raw"]["lidar"], +1),
    ("LiDAR_residualised", lambda m: m["lidar_residualised"]["lidar"], +1),
]

PAIRS = [
    ("D2", "s42", "H42", "I42"), ("D2", "s43", "H43", "I43"), ("D2", "s44", "H44", "I44"),
    ("D1", "s42", "P42", "F42"), ("D1", "s43", "P43", "F43"), ("D1", "s44", "P44", "F44"),
]


def two_sided_binomial(k, n, p=0.5):
    """Exact two-sided binomial p-value against a fair coin (sum of outcomes at most as likely)."""
    probs = [comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(n + 1)]
    obs = probs[k]
    return float(sum(q for q in probs if q <= obs + 1e-12))


def main(paths):
    data = {}
    for p in paths:
        d = json.load(open(p))
        data.update({k: v for k, v in d.items() if not k.startswith("_")})

    truth = {k: v["points"]["untrained40"]["top_cca"] for k, v in data.items()}

    print("=" * 108)
    print("GROUND TRUTH -- top-CCA with the 40 untrained targets (residualised)")
    print("=" * 108)
    for exp, seed, a, b in PAIRS:
        d = truth[a] - truth[b]
        print(f"  {exp} {seed}: {a}={truth[a]:.4f}  {b}={truth[b]:.4f}   "
              f"delta({a}-{b})={d:+.4f}   winner={a if d > 0 else b}")

    print()
    print("=" * 108)
    print("SELECTION-RULE SCORE -- does the metric pick the arm that carries more molecular information?")
    print("=" * 108)
    hdr = f"{'metric':36s}" + "".join(f"{e+' '+s:>10s}" for e, s, _, _ in PAIRS) + f"{'D2':>6s}{'D1':>6s}{'ALL':>6s}{'p':>8s}"
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for name, fn, direction in METRICS:
        marks, correct = [], {"D2": 0, "D1": 0}
        for exp, seed, a, b in PAIRS:
            va, vb = fn(data[a]["metrics"]), fn(data[b]["metrics"])
            picked_a = (va > vb) if direction > 0 else (va < vb)
            true_a = truth[a] > truth[b]
            ok = picked_a == true_a
            correct[exp] += ok
            marks.append("OK" if ok else "MISS")
        tot = correct["D2"] + correct["D1"]
        rows.append((name, correct["D2"], correct["D1"], tot))
        p = two_sided_binomial(tot, len(PAIRS))
        print(f"{name:36s}" + "".join(f"{m:>10s}" for m in marks) +
              f"{correct['D2']:>4d}/3{correct['D1']:>4d}/3{tot:>4d}/6{p:>8.3f}")

    print()
    print("=" * 108)
    print("WITHIN-ARM SEED VARIABILITY -- the dynamic range of each metric under a nuisance change")
    print("(same arm, same objective, same data, same schedule; only the training seed differs.")
    print(" RankMe explicitly reserves this comparison for itself: 'should only be used to compare")
    print(" different runs of a given method'.)")
    print("=" * 108)
    arms = {"D2 H (Hallmark)": ["H42", "H43", "H44"], "D2 I (PBS)": ["I42", "I43", "I44"],
            "D1 P (programme_only)": ["P42", "P43", "P44"], "D1 F (programme_free)": ["F42", "F43", "F44"]}
    print(f"{'arm':24s}{'eff_rank (3 seeds)':>34s}{'fold':>7s}   {'untrained40 CCA (3 seeds)':>32s}{'rel.range':>10s}")
    for arm, ls in arms.items():
        er = [data[l]["metrics"]["effective_rank_residualised"] for l in ls]
        cc = [truth[l] for l in ls]
        print(f"{arm:24s}" + f"{'  '.join('%7.3f' % v for v in er):>34s}"
              f"{max(er)/min(er):>6.2f}x   " + f"{'  '.join('%7.4f' % v for v in cc):>32s}"
              f"{(max(cc)-min(cc))/np.mean(cc)*100:>9.1f}%")

    print()
    print("The between-arm information gap each metric is being asked to resolve:")
    for exp, seed, a, b in PAIRS:
        print(f"   {exp} {seed}: |delta CCA| = {abs(truth[a]-truth[b]):.4f} "
              f"({abs(truth[a]-truth[b])/np.mean([truth[a],truth[b]])*100:.1f}% relative)")

    print()
    print("=" * 108)
    print("POOLED ASSOCIATION across all 12 artifacts (reported for completeness; the paired")
    print("comparison above is the operative one, since pooling mixes within- and between-arm variance)")
    print("=" * 108)
    from scipy.stats import spearmanr, pearsonr
    labels = [l for _, _, a, b in PAIRS for l in (a, b)]
    y = np.array([truth[l] for l in labels])
    for name, fn, direction in METRICS:
        x = np.array([fn(data[l]["metrics"]) for l in labels])
        rho, prho = spearmanr(x, y)
        r, pr = pearsonr(x, y)
        print(f"{name:36s} spearman={rho:+.3f} (p={prho:.3f})   pearson={r:+.3f} (p={pr:.3f})")

    # split by experiment, because D1 and D2 have opposite structure
    print()
    for exp in ("D2", "D1"):
        ls = [l for e, _, a, b in PAIRS if e == exp for l in (a, b)]
        yy = np.array([truth[l] for l in ls])
        print(f"--- {exp} only (n=6) ---")
        for name, fn, direction in METRICS:
            xx = np.array([fn(data[l]["metrics"]) for l in ls])
            rho, prho = spearmanr(xx, yy)
            print(f"  {name:36s} spearman={rho:+.3f} (p={prho:.3f})")

    # noise floor, if the subsampled run supplied it
    if any("subsample" in v for v in data.values()):
        print()
        print("=" * 108)
        print("METRIC NOISE FLOOR vs BETWEEN-ARM GAP")
        print("(subsampling patients at 80%, 40 draws. If the between-arm gap in a metric is smaller")
        print(" than the metric's own sampling spread, that metric cannot resolve the pair at all.)")
        print("=" * 108)
        keys = ["effective_rank_residualised", "rankme_raw",
                "participation_ratio_residualised", "lidar_residualised"]
        for k in keys:
            print(f"\n  --- {k} ---")
            for exp, seed, a, b in PAIRS:
                if "subsample" not in data[a]:
                    continue
                sa, sb = data[a]["subsample"][k], data[b]["subsample"][k]
                gap = abs(sa["mean"] - sb["mean"])
                pooled = np.hypot(sa["sd"], sb["sd"])
                print(f"    {exp} {seed}: {a}={sa['mean']:8.3f}+-{sa['sd']:6.3f}  "
                      f"{b}={sb['mean']:8.3f}+-{sb['sd']:6.3f}  gap={gap:7.3f}  "
                      f"gap/sd={gap/pooled if pooled > 0 else float('nan'):6.2f}")


if __name__ == "__main__":
    main(sys.argv[1:])
