"""Two analyses that the selection-rule count (n=6, underpowered) cannot deliver.

VENDORED 2026-08-05 from `~/e0_run/p2_necessity_and_variance.py`, unchanged apart from this
paragraph. Produces `paper/P2_RANK_DRAFT.md` §4.2 (the variance decomposition — the paper's
most important display item) and §4.7 (the necessity test). Reads only the JSON written by
`p2_competing_metrics.py`. The two pre-declared thresholds below were fixed before the pair
list was inspected; see `NOTEBOOK_ENTRIES/PREDECLARED_D1_necessity_test_20260803T2300Z.md`.

(A) VARIANCE DECOMPOSITION. Across the 12 frozen artifacts (2 experiments x 2 arms x 3 seeds),
    how much of each metric's variance is carried by the ARM -- the factor that actually changes
    the molecular information -- and how much by the SEED, a nuisance factor that by construction
    changes nothing about the objective, data, architecture or schedule?

    A proxy is usable for arm selection only if the arm term dominates the seed term. This is a
    strictly stronger and better-powered question than "did it pick the right arm 5 times out of
    6", because it uses the magnitudes rather than only the signs, and because the seed term is
    estimated from 6 within-arm contrasts rather than 6 paired comparisons.

    Reported as omega-like ratio  SS_arm / (SS_arm + SS_seed_within_arm)  on the log scale for
    the rank-type metrics (they are multiplicative, spanning 6.4 to 34.1) and on the raw scale
    for the CCA ground truth.

(B) NECESSITY SCAN. RankMe's defence is that high rank is "a necessary (but not sufficient)
    condition for good downstream performance". Under that hedge, high-rank/low-information is
    PREDICTED and is not a counterexample. The only configuration that breaks it is
    LOW RANK carrying HIGH information. This scans every ordered pair of artifacts for that
    configuration and reports the strongest instances, separated by whether the pair is
    within-arm (the comparison RankMe explicitly reserves for itself), within-experiment, or
    cross-experiment.
"""
from __future__ import annotations

import json
import sys
from itertools import combinations

import numpy as np

ARM_OF = {"H": ("D2", "H"), "I": ("D2", "I"), "P": ("D1", "P"), "F": ("D1", "F")}
LABELS = ["H42", "H43", "H44", "I42", "I43", "I44", "P42", "P43", "P44", "F42", "F43", "F44"]

# Pre-declared thresholds for what counts as a necessity violation. These are fixed here
# BEFORE inspecting the pair list, and are chosen from quantities already established
# independently of this analysis:
#   RANK_FOLD  -- the smallest fold-change in effective rank this project has ever treated as
#                 meaningful. 2.0x is below the 2.10-3.75x that the SEED alone produces within a
#                 single arm, so it is a conservative floor: anything at or above it is at least
#                 as large as the metric's own nuisance range.
#   CCA_DELTA  -- 0.0705, the SMALLEST between-arm channel gap measured on matched arms in this
#                 project (D1 seed 42). A rank deficit paired with an information surplus at
#                 least this large is a gap the project has already accepted as real elsewhere.
RANK_FOLD = 2.0
CCA_DELTA = 0.0705


def main(paths):
    data = {}
    for p in paths:
        d = json.load(open(p))
        data.update({k: v for k, v in d.items() if not k.startswith("_")})

    er = {l: data[l]["metrics"]["effective_rank_residualised"] for l in LABELS}
    rm = {l: data[l]["metrics"]["rankme_raw"] for l in LABELS}
    cca = {l: data[l]["points"]["untrained40"]["top_cca"] for l in LABELS}

    # ---------------- (A) variance decomposition ----------------
    print("=" * 104)
    print("(A) VARIANCE DECOMPOSITION -- arm (signal) vs training seed (nuisance)")
    print("=" * 104)
    print("Design: 4 arms (D2-H, D2-I, D1-P, D1-F) x 3 seeds (42/43/44), one frozen artifact each.")
    print("The seed changes nothing about objective, architecture, data, split or schedule.")
    print()

    def decompose(values, log=False):
        v = {l: (np.log(values[l]) if log else values[l]) for l in LABELS}
        groups = {}
        for l in LABELS:
            groups.setdefault(l[0], []).append(v[l])
        grand = np.mean([v[l] for l in LABELS])
        ss_arm = sum(len(g) * (np.mean(g) - grand) ** 2 for g in groups.values())
        ss_seed = sum(sum((x - np.mean(g)) ** 2 for x in g) for g in groups.values())
        return ss_arm, ss_seed

    print(f"{'quantity':40s}{'scale':>8s}{'SS_arm':>12s}{'SS_seed':>12s}{'arm share':>12s}{'F(3,8)':>10s}")
    print("-" * 104)
    for name, vals, log in [
        ("effective rank (residualised)", er, True),
        ("RankMe (raw, as published)", rm, True),
        ("GROUND TRUTH: untrained40 top-CCA", cca, False),
    ]:
        a, s = decompose(vals, log=log)
        share = a / (a + s)
        f = (a / 3) / (s / 8) if s > 0 else float("inf")
        print(f"{name:40s}{'log' if log else 'raw':>8s}{a:>12.4f}{s:>12.4f}{share:>11.1%}{f:>10.2f}")

    print()
    print("Reading: the ground-truth information channel is almost entirely an ARM effect --")
    print("the seed barely moves it. If effective rank were a proxy for that channel it would have")
    print("to show the same structure. The arm share above is what says whether it does.")

    # per-arm coefficient of variation on the log scale = multiplicative seed noise
    print()
    print(f"{'arm':16s}{'eff.rank seeds':>34s}{'fold':>8s}{'CCA seeds':>30s}{'fold':>8s}")
    print("-" * 104)
    for a in ("H", "I", "P", "F"):
        ls = [l for l in LABELS if l[0] == a]
        e = [er[l] for l in ls]
        c = [cca[l] for l in ls]
        print(f"{a:16s}" + f"{'  '.join('%7.3f' % x for x in e):>34s}{max(e)/min(e):>7.2f}x"
              f"{'  '.join('%7.4f' % x for x in c):>30s}{max(c)/min(c):>7.3f}x")

    # ---------------- (B) necessity scan ----------------
    print()
    print("=" * 104)
    print("(B) NECESSITY SCAN -- looking for LOW RANK carrying HIGH INFORMATION")
    print("=" * 104)
    print("RankMe's hedge: high rank is 'a necessary (but not sufficient) condition for good")
    print("downstream performance'. High-rank/low-information is therefore PREDICTED by the hedge")
    print("and is NOT a counterexample. Only the reverse breaks it.")
    print()
    print(f"PRE-DECLARED violation criterion (fixed before inspecting the pairs):")
    print(f"  a pair (lo, hi) is a violation iff  eff_rank(hi)/eff_rank(lo) >= {RANK_FOLD}")
    print(f"  AND  CCA(lo) - CCA(hi) >= {CCA_DELTA}")
    print(f"  i.e. the LOWER-rank artifact carries at least as much more information as the")
    print(f"  smallest between-arm gap this project has accepted as real (D1 seed 42, +0.0705),")
    print(f"  while sitting at least {RANK_FOLD}x lower in rank.")
    print()

    def scope(a, b):
        if a[0] == b[0]:
            return "WITHIN-ARM (RankMe's own reserved scope)"
        if ARM_OF[a[0]][0] == ARM_OF[b[0]][0]:
            return "within-experiment, across arms"
        return "cross-experiment"

    hits = []
    for a, b in combinations(LABELS, 2):
        lo, hi = (a, b) if er[a] < er[b] else (b, a)
        fold = er[hi] / er[lo]
        dcca = cca[lo] - cca[hi]
        if fold >= RANK_FOLD and dcca >= CCA_DELTA:
            hits.append((dcca, fold, lo, hi))
    hits.sort(reverse=True)

    if not hits:
        print("  NO VIOLATIONS FOUND under the pre-declared criterion.")
    else:
        print(f"  {len(hits)} violating pair(s) of {len(list(combinations(LABELS, 2)))} examined:")
        print()
        print(f"  {'lower-rank artifact':>22s}{'rank':>9s}{'CCA':>9s}   |"
              f"{'higher-rank artifact':>22s}{'rank':>9s}{'CCA':>9s}   {'fold':>7s}{'dCCA':>9s}  scope")
        for dcca, fold, lo, hi in hits:
            print(f"  {lo:>22s}{er[lo]:>9.3f}{cca[lo]:>9.4f}   |{hi:>22s}{er[hi]:>9.3f}{cca[hi]:>9.4f}"
                  f"   {fold:>6.2f}x{dcca:>+9.4f}  {scope(lo, hi)}")

    # the weaker but in-scope statement: large rank change, NO information change
    print()
    print("-" * 104)
    print("SECONDARY (does not break necessity, but breaks USEFULNESS): large rank change with a")
    print("negligible information change, WITHIN a single arm -- i.e. two runs of the same method.")
    print("-" * 104)
    print(f"  {'pair':>20s}{'rank lo':>10s}{'rank hi':>10s}{'fold':>8s}{'CCA lo':>10s}{'CCA hi':>10s}{'dCCA':>10s}")
    for a in ("H", "I", "P", "F"):
        ls = [l for l in LABELS if l[0] == a]
        for x, y in combinations(ls, 2):
            lo, hi = (x, y) if er[x] < er[y] else (y, x)
            print(f"  {lo+' vs '+hi:>20s}{er[lo]:>10.3f}{er[hi]:>10.3f}{er[hi]/er[lo]:>7.2f}x"
                  f"{cca[lo]:>10.4f}{cca[hi]:>10.4f}{cca[hi]-cca[lo]:>+10.4f}")


if __name__ == "__main__":
    main(sys.argv[1:])
