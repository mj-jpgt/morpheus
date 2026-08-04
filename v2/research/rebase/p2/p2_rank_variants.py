"""The three statistics this repository calls `effective_rank`, scored as selection rules.

Produces `paper/P2_RANK_DRAFT.md` §4.5(a) — the table showing that the functions this repository
calls `effective_rank` return different verdicts on the same matched pairs. The draft says 3 of
6; measured under the canonical statistics it is **2 of 6** (D2 s43 and D2 s44), for the reason
in note 2 below.

VENDORED 2026-08-04 from `~/e0_run/p2_rank_variants.py`. Two changes, both deliberate and both
recorded in `NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260804T0255Z.md`:

1. **The three rank statistics are no longer reimplemented inline.** The original defined its
   own `R1`, `R2`, `R3` and prefixed the file with `sys.path.insert(0, "/home/ubuntu/ws")` —
   a workspace whose `spectral.py` predated the rank canonicalisation entirely
   (`NOTEBOOK_ENTRIES/WORKSPACE_DRIFT_AUDIT_ALL_20260803T2359Z.md`). This version calls
   `effective_rank(..., variant=RANK_VARIANTS[...])`, the one implementation, so the numbers
   are produced by the same function object the rest of the paper uses and
   `v2/tests/test_effective_rank_canonical.py` polices.

2. **The original's `R2` and `R3` were not R2 and R3.** They were
   `(sum s^2)^2 / sum s^4` — the order-2 Hill number of the EIGENVALUE distribution — where
   `d1_audit.py`'s R2, and therefore `RANK_VARIANTS["R2"]`, is `(sum s)^2 / sum s^2`, the
   order-2 Hill number of the SINGULAR-VALUE distribution. They are different statistics.
   Rather than silently substitute one for the other, this script reports **both**: the
   canonical R1/R2/R3 columns and, beside them, `PR` and `PR_rownorm`, which are exactly what
   the pre-vendoring script computed for its "R2" and "R3" columns. Whether the published
   §4.5(a) verdicts move is then a measurement rather than an assumption.

   **They move.** On the twelve frozen artifacts, `PR` and `PR_rownorm` reproduce the draft's
   published "R2" and "R3" rows cell for cell, and the canonical R2 and R3 do not:
   canonical R2 scores D1 3/3 where the draft reports 2/3, and canonical R3 swaps its two D2
   verdicts. The canonical R2 and R3 *levels* in draft §4.5(b) are unaffected — they were
   already computed with the canonical function and reproduce exactly.

usage:
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python3 p2_rank_variants.py --targets <frozen_rna_targets.npz> LABEL=path.npz ...
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank, top_canonical_correlation
from morpheus.v2.research.rebase.p2.p2_competing_metrics import (participation_ratio,
                                                                 participation_ratio_rownorm)

#: The 40 targets neither arm was supervised on, per
#: `NOTEBOOK_ENTRIES/d1_readout_preregistration_20260803T1700Z.md`.
UNTRAINED_40 = ["heldout_pathway", "immune_tme", "tumour_state"]

#: Column label -> callable. R1/R2/R3 are the canonical variants; PR/PR_rownorm are the two
#: statistics the pre-vendoring script actually computed under the labels "R2" and "R3".
STATISTICS = {
    "R1": lambda x: effective_rank(x, variant=RANK_VARIANTS["R1"]),
    "R2": lambda x: effective_rank(x, variant=RANK_VARIANTS["R2"]),
    "R3": lambda x: effective_rank(x, variant=RANK_VARIANTS["R3"]),
    "PR": participation_ratio,
    "PR_rownorm": participation_ratio_rownorm,
}

DEFAULT_PAIRS = [("D2 s42", "H42", "I42"), ("D2 s43", "H43", "I43"), ("D2 s44", "H44", "I44"),
                 ("D1 s42", "P42", "F42"), ("D1 s43", "P43", "F43"), ("D1 s44", "P44", "F44")]


def load_targets(path, groups=UNTRAINED_40):
    raw = np.load(path, allow_pickle=False)
    index = {p: i for i, p in enumerate(np.asarray(raw["patient_ids"]).astype(str))}
    keep = np.isin(np.asarray(raw["target_groups"]).astype(str), list(groups))
    return index, np.asarray(raw["scores"], dtype=np.float64)[:, keep]


def score_artifact(path, index, values, seed=42, block="wsi_biology"):
    """Every statistic plus the ground-truth channel, on the residualised held-out block."""
    raw = np.load(path, allow_pickle=False)
    test = np.asarray(raw["split"]).astype(str) == "test"
    ids = np.asarray(raw["patient_ids"]).astype(str)[test]
    cancers = np.asarray(raw["cancers"]).astype(str)[test]
    tss, _ = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    rx = cross_fitted_residuals(np.asarray(raw[block], dtype=np.float64)[test], design, seed=seed)
    ry = cross_fitted_residuals(values[[index[i] for i in ids]], design, seed=seed)
    record = {name: float(fn(rx)) for name, fn in STATISTICS.items()}
    record["cca40"] = top_canonical_correlation(rx, ry, n_components=16)
    return record


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--output", default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("artifacts", nargs="+", help="LABEL=path.npz")
    a = ap.parse_args(argv)

    index, values = load_targets(a.targets)
    result = {}
    for item in a.artifacts:
        label, path = item.split("=", 1)
        result[label] = score_artifact(path, index, values, seed=a.seed)

    names = list(STATISTICS)
    print(f"{'artifact':>9}" + "".join(f"{n:>12}" for n in names) + f"{'CCA40':>10}")
    for label, record in result.items():
        print(f"{label:>9}" + "".join(f"{record[n]:>12.3f}" for n in names)
              + f"{record['cca40']:>10.4f}", flush=True)

    pairs = [p for p in DEFAULT_PAIRS if p[1] in result and p[2] in result]
    print("\nSELECTION RULE by rank statistic (does higher rank = more information?)")
    print(f"{'stat':>12}" + "".join(f"{n:>9}" for n, _, _ in pairs) + f"{'D2':>6}{'D1':>6}{'ALL':>6}")
    verdicts = {}
    for name in names:
        marks, d2, d1 = [], 0, 0
        for pair, a_label, b_label in pairs:
            ok = ((result[a_label][name] > result[b_label][name])
                  == (result[a_label]["cca40"] > result[b_label]["cca40"]))
            marks.append("OK" if ok else "MISS")
            if pair.startswith("D2"):
                d2 += ok
            else:
                d1 += ok
        verdicts[name] = {"marks": marks, "D2": d2, "D1": d1, "ALL": d2 + d1}
        print(f"{name:>12}" + "".join(f"{m:>9}" for m in marks) + f"{d2:>4}/3{d1:>4}/3{d2 + d1:>4}/6")

    if a.output:
        Path(a.output).parent.mkdir(parents=True, exist_ok=True)
        Path(a.output).write_text(
            json.dumps({"artifacts": result, "verdicts": verdicts}, indent=2), encoding="utf-8")
        print("wrote", a.output)
    return {"artifacts": result, "verdicts": verdicts}


if __name__ == "__main__":
    main(sys.argv[1:])
