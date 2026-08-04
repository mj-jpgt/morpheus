"""Two robustness checks a referee will demand of the ground truth and of the rank statistic.

VENDORED 2026-08-05 from `~/e0_run/p2_robustness.py`; produces `paper/P2_RANK_DRAFT.md` §4.5(c).
One change from the original: the output path was hardcoded to
`/home/ubuntu/e0_run/P2_ROBUSTNESS.json` and is now the first positional argument, so the
script can be run from a checkout on any machine. Nothing else was touched.

usage: OMP_NUM_THREADS=1 ... python3 p2_robustness.py <out.json> <targets.npz> <seed> LABEL=path...

(1) HELD-OUT ground truth. The headline channel statistic is an IN-SAMPLE top canonical
    correlation at a 16-component budget, which is a maximum over 16 in-sample directions and
    is upward-biased at finite n. If the arm ordering only survives in sample, the whole
    comparison is worthless. `heldout_top_cca` fits the canonical directions on one half of the
    patients and scores them on the other, which is unbiased. The arm ordering must survive.

(2) REPRESENTATION CHOICE. Both the rank and the channel are read off `wsi_biology`. If the
    conclusion flips on `rna_biology` or `full_biology` that is a finding, not a detail.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import effective_rank, top_canonical_correlation, heldout_top_cca

GROUP = ["heldout_pathway", "immune_tme", "tumour_state"]


def load_targets(path):
    raw = np.load(path, allow_pickle=False)
    ids = np.asarray(raw["patient_ids"]).astype(str)
    values = np.asarray(raw["scores"], dtype=np.float64)
    available = np.asarray(raw["target_groups"]).astype(str)
    keep = np.isin(available, GROUP)
    return {p: i for i, p in enumerate(ids)}, values[:, keep]


def main(output_path, targets_path, seed, items):
    index, tvalues = load_targets(targets_path)
    rows_out = {}
    for item in items:
        label, path = item.split("=", 1)
        raw = np.load(path, allow_pickle=False)
        test = raw["split"].astype(str) == "test"
        ids = raw["patient_ids"].astype(str)[test]
        cancers = raw["cancers"].astype(str)[test]
        tss, _ = pooled_tissue_source_site(ids, min_site_count=10)
        design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
        rows = np.asarray([index[p] for p in ids])
        ry = cross_fitted_residuals(tvalues[rows], design, seed=seed)
        rec = {}
        for view in ("wsi_biology", "rna_biology", "full_biology"):
            x = raw[view].astype(np.float64)[test]
            rx = cross_fitted_residuals(x, design, seed=seed)
            rec[view] = {
                "eff_rank": effective_rank(rx),
                "cca_insample": top_canonical_correlation(rx, ry, n_components=16),
                "cca_heldout": heldout_top_cca(rx, ry, n_components=16, seed=seed),
            }
        rows_out[label] = rec
        print(label, " ".join(
            f"[{v[:3]}] er={r['eff_rank']:7.3f} cca_in={r['cca_insample']:.4f} cca_out={r['cca_heldout']:.4f}"
            for v, r in rec.items()), flush=True)

    print()
    print("=" * 112)
    print("ARM ORDERING under each representation and each CCA estimator")
    print("(ground truth must survive the held-out estimator or it is not a ground truth)")
    print("=" * 112)
    pairs = [("D2 s42", "H42", "I42"), ("D2 s43", "H43", "I43"), ("D2 s44", "H44", "I44"),
             ("D1 s42", "P42", "F42"), ("D1 s43", "P43", "F43"), ("D1 s44", "P44", "F44")]
    for view in ("wsi_biology", "rna_biology", "full_biology"):
        print(f"\n--- {view} ---")
        print(f"{'pair':10s}{'cca_in A':>10s}{'cca_in B':>10s}{'d_in':>9s}"
              f"{'cca_out A':>11s}{'cca_out B':>11s}{'d_out':>9s}   {'rank A':>8s}{'rank B':>8s}"
              f"   {'info winner':>12s}{'rank winner':>12s}{'agree?':>8s}")
        for name, a, b in pairs:
            if a not in rows_out or b not in rows_out:
                continue
            ra, rb = rows_out[a][view], rows_out[b][view]
            din = ra["cca_insample"] - rb["cca_insample"]
            dout = ra["cca_heldout"] - rb["cca_heldout"]
            iw = a if dout > 0 else b
            rw = a if ra["eff_rank"] > rb["eff_rank"] else b
            print(f"{name:10s}{ra['cca_insample']:>10.4f}{rb['cca_insample']:>10.4f}{din:>+9.4f}"
                  f"{ra['cca_heldout']:>11.4f}{rb['cca_heldout']:>11.4f}{dout:>+9.4f}"
                  f"   {ra['eff_rank']:>8.3f}{rb['eff_rank']:>8.3f}   {iw:>12s}{rw:>12s}"
                  f"{('YES' if iw == rw else 'NO'):>8s}")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(json.dumps(rows_out, indent=2), encoding="utf-8")
    return rows_out


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4:])
