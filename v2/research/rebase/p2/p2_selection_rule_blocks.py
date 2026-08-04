"""Re-score the selection rule with EACH of D2's six target blocks as ground truth.

`paper/P2_RANK_DRAFT.md` §4.6 scores every label-free metric as a selection rule against
one ground truth: the held-out top-CCA on the 40 gene-set targets neither arm trained on.
`NOTEBOOK_ENTRIES/d2_coordinate_system_result_20260804T0800Z.md` §1a then established that
the D2 arm contrast is **block-dependent** -- the −0.12 exists on the gene-set block and on
no other target space on disk, and rotating the exam into the dictionary's own coordinates
moves the contrast by the entire size of the published effect.

If the arm contrast is block-dependent then so is the ground truth, and so is every OK/MISS
mark in §4.6's table. That consequence had not been chased. This script chases it: it holds
the metrics fixed and swaps the truth, once per block, and reports whether the verdict is
stable.

**Nothing here computes a rank, a channel or an interval.** The metrics come from
`P2_METRICS_D{1,2}.json` (written by `p2_competing_metrics.py`); the per-block arm contrasts
come from `d2_coordinate_system/out/EXAM_PANEL.json` (written on the box by `exam_panel.py`
calling `calibra`). The rule itself is imported from `p2_selection_rule` so there is exactly
one definition of "does this metric pick the better arm".

**Only the D2 half can be re-scored.** The exam panel re-scored the two D2 arms on six target
blocks; the D1 arms were never re-scored on anything but the gene sets. The D1 column is
therefore held at the gene-set truth in every row and is labelled as not re-scored, rather
than being quietly carried as though it had been.

    OMP_NUM_THREADS=1 python -m morpheus.v2.research.rebase.p2.p2_selection_rule_blocks \\
        P2_METRICS_D2.json P2_METRICS_D1.json EXAM_PANEL.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict

from .p2_selection_rule import METRICS, two_sided_binomial

#: Exam blocks in `EXAM_PANEL.json`, in the order §1a of the coordinate-system entry
#: reports them: the published one first, then the four 128-column expression-derived
#: blocks, then the published negative control.
BLOCKS = [
    ("geneset_untrained40", "gene sets 40", "gene sets, untrained 40 (PUBLISHED TRUTH)"),
    ("pbs_codes128_ARM_I_OWN", "PBS codes", "PBS codes 128 (arm I's own supervision)"),
    ("pca_basis128", "PCA basis", "PCA basis 128"),
    ("pbs_shuffled_s1", "shuf. PBS", "gene-label-shuffled PBS 128"),
    ("geneset_random_control", "rand ctrl", "random_control gene sets 90 (NEGATIVE CONTROL)"),
    ("random_dictionary128", "rand dict", "size/spectrum-matched random dictionary 128"),
]

D2_PAIRS = [("D2 s42", 42, "H42", "I42"), ("D2 s43", 43, "H43", "I43"),
            ("D2 s44", 44, "H44", "I44")]
D1_PAIRS = [("D1 s42", "P42", "F42"), ("D1 s43", "P43", "F43"), ("D1 s44", "P44", "F44")]


def load_truth(exam_panel_path: str) -> dict[str, dict[int, float]]:
    """block -> seed -> Δ (PBS − Hallmark). Negative means arm H carries more."""
    rows = json.load(open(exam_panel_path))
    truth: dict[str, dict[int, float]] = defaultdict(dict)
    for r in rows:
        assert r["block"] == "residualised", r
        truth[r["exam"]][int(r["seed"])] = float(r["pbs_minus_hallmark"])
    for key, _short, _label in BLOCKS:
        assert sorted(truth[key]) == [42, 43, 44], (key, sorted(truth.get(key, {})))
    return truth


def marks_for(data: dict, truth: dict[str, dict[int, float]], block: str,
              geneset_d1: dict) -> dict[str, dict]:
    """OK/MISS per metric, with the D2 truth taken from `block`."""
    out = {}
    for name, fn, direction in METRICS:
        marks = []
        for _label, seed, a, b in D2_PAIRS:
            va, vb = fn(data[a]["metrics"]), fn(data[b]["metrics"])
            picked_a = (va > vb) if direction > 0 else (va < vb)
            # Δ is PBS − Hallmark, so Δ < 0 means the FIRST-named arm (H) is the
            # information winner. The sign convention is asserted, not assumed:
            # a source edit that reversed it would flip every mark silently.
            true_a = truth[block][seed] < 0
            marks.append("OK" if picked_a == true_a else "MISS")
        d1_marks = geneset_d1[name]
        out[name] = {"d2_marks": marks, "d1_marks": d1_marks,
                     "d2": sum(m == "OK" for m in marks),
                     "d1": sum(m == "OK" for m in d1_marks),
                     "all": sum(m == "OK" for m in marks + d1_marks)}
    return out


def d1_marks(data: dict) -> dict[str, list[str]]:
    """The D1 half, at the gene-set truth, which is the only truth it has."""
    out = {}
    for name, fn, direction in METRICS:
        marks = []
        for _label, a, b in D1_PAIRS:
            va, vb = fn(data[a]["metrics"]), fn(data[b]["metrics"])
            picked_a = (va > vb) if direction > 0 else (va < vb)
            true_a = (data[a]["points"]["untrained40"]["top_cca"]
                      > data[b]["points"]["untrained40"]["top_cca"])
            marks.append("OK" if picked_a == true_a else "MISS")
        out[name] = marks
    return out


def main(argv: list[str]) -> int:
    metrics_paths, panel_path = argv[:-1], argv[-1]
    data: dict = {}
    for p in metrics_paths:
        data.update({k: v for k, v in json.load(open(p)).items() if not k.startswith("_")})
    truth = load_truth(panel_path)
    fixed_d1 = d1_marks(data)

    per_block = {k: marks_for(data, truth, k, fixed_d1) for k, _s, _l in BLOCKS}

    def pattern_of(key: str) -> str:
        return "".join("H" if truth[key][s] < 0 else "I" for _l, s, _a, _b in D2_PAIRS)

    print("=" * 112)
    print("GROUND TRUTH IS BLOCK-DEPENDENT -- the D2 arm contrast on each of six target blocks")
    print("=" * 112)
    print(f"{'exam block used as ground truth':46s}"
          + "".join(f"{'s' + str(s):>12s}" for _l, s, _a, _b in D2_PAIRS)
          + "   winner")
    for key, _short, label in BLOCKS:
        print(f"{label:46s}"
              + "".join(f"{truth[key][s]:>+12.4f}" for _l, s, _a, _b in D2_PAIRS)
              + f"   {pattern_of(key)}")

    patterns = sorted({pattern_of(k) for k, _s, _l in BLOCKS})
    print(f"\n{len(patterns)} distinct winner pattern(s) across the six blocks: "
          f"{', '.join(patterns)}")
    flips = [s for _l, s, _a, _b in D2_PAIRS
             if len({truth[k][s] < 0 for k, _s, _l in BLOCKS}) > 1]
    print(f"seed(s) whose winner changes with the block: "
          f"{', '.join(str(s) for s in flips) if flips else 'none'}")

    print()
    print("=" * 112)
    print("SELECTION-RULE SCORE UNDER EACH GROUND TRUTH")
    print("D2 columns are re-scored per block. The D1 columns CANNOT be: the D1 arms were never")
    print("re-scored on any block but the gene sets, so D1 is held at the gene-set truth in every")
    print("row and is NOT evidence that the D1 half is block-stable.")
    print("=" * 112)
    hdr = f"{'metric':36s}" + "".join(f"{short:>12s}" for _k, short, _l in BLOCKS)

    for scope, n, note in (("d2", 3, "D2 pairs only, the half that can be re-scored"),
                           ("all", 6, "D2 re-scored + D1 held at the gene-set truth")):
        print(f"\n--- {scope.upper()} count out of {n}  ({note}) ---")
        print(hdr + "   stable?")
        print("-" * (len(hdr) + 12))
        for name, _fn, _d in METRICS:
            vals = [per_block[k][name][scope] for k, _s, _l in BLOCKS]
            stable = "yes" if len(set(vals)) == 1 else f"NO ({min(vals)}-{max(vals)})"
            print(f"{name:36s}" + "".join(f"{v}/{n}".rjust(12) for v in vals) + f"   {stable}")

    print()
    print("=" * 112)
    print("WHAT MOVES, AND WHAT IT COSTS THE TABLE")
    print("=" * 112)
    unstable = [name for name, _f, _d in METRICS
                if len({per_block[k][name]["d2"] for k, _s, _l in BLOCKS}) > 1]
    print(f"{len(unstable)} of {len(METRICS)} metric rows change their D2 count when the "
          f"ground-truth block changes.")

    print("\nExact two-sided binomial on the ALL column, best against worst block "
          "(6/6 is p = 0.031 -- 'significant' by this design's own bar):")
    for name, _f, _d in METRICS:
        row = {short: per_block[k][name]["all"] for k, short, _l in BLOCKS}
        best, worst = max(row.values()), min(row.values())
        if best == worst:
            continue
        print(f"   {name:36s} best {best}/6 p={two_sided_binomial(best, 6):.3f} "
              f"({[b for b, v in row.items() if v == best][0]})   "
              f"worst {worst}/6 p={two_sided_binomial(worst, 6):.3f} "
              f"({[b for b, v in row.items() if v == worst][0]})")

    print("\nOrdering between the two rows the draft's 4.6 names explicitly:")
    for k, _short, label in BLOCKS:
        er = per_block[k]["effective_rank_residualised"]["d2"]
        rm = per_block[k]["rankme_raw (as published)"]["d2"]
        verdict = ("RankMe ahead" if rm > er else
                   "canonical effective rank ahead" if er > rm else "tied")
        print(f"   {label:46s} canonical R1 {er}/3 vs RankMe {rm}/3   -> {verdict}")

    print("\nThe published table is the `geneset_untrained40` column, reproduced here unchanged:")
    for name, _f, _d in METRICS:
        r = per_block["geneset_untrained40"][name]
        print(f"   {name:36s} {' '.join(r['d2_marks']):>15s}   D2 {r['d2']}/3   "
              f"D1 {r['d1']}/3   ALL {r['all']}/6")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
