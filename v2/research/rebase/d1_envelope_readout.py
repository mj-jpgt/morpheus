"""Retraining envelope: rank AND channel, per repeat, on identical configurations.

The asymmetry is the point. If rank spreads widely across re-runs of the SAME
configuration while the channel does not, that is the paper's thesis measured on
one pair of arms rather than argued across tables.

Statistic: canonical Roy & Vetterli order 1, centred, rows at own norms, imported
from v2.calibra.spectral. Nothing computed inline -- a fourth statistic was found
hiding under the R2/R3 names, and an inline formula is the tell.
Block: exported artifact `wsi_biology`, held-out test partition, cancer +
pooled-TSS residualised, top-CCA at 16 components -- identical to d2_compare.
"""
import glob, sys
import numpy as np
import pandas as pd
from morpheus.v2.calibra.spectral import effective_rank, top_canonical_correlation
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)

TARGETS = "/home/ubuntu/e0_run/data/frozen_rna_targets.npz"
STRATIFIED = ("heldout_pathway", "immune_tme", "tumour_state")
pattern = sys.argv[1] if len(sys.argv) > 1 else "/home/ubuntu/e0_run/d1_envelope/rep*.npz"

raw = np.load(TARGETS, allow_pickle=False)
tids = np.asarray(raw["patient_ids"]).astype(str)
groups = np.asarray(raw["target_groups"]).astype(str)
keep = np.isin(groups, list(STRATIFIED))
scores = np.asarray(raw["scores"], dtype=np.float64)[:, keep]
index = {p: i for i, p in enumerate(tids)}
print(f"targets: {scores.shape[1]} ({', '.join(STRATIFIED)})\n")

rows = []
for path in sorted(glob.glob(pattern)):
    art = np.load(path, allow_pickle=False)
    test = art["split"].astype(str) == "test"
    ids = art["patient_ids"].astype(str)[test]
    states = np.asarray(art["wsi_biology"], dtype=np.float64)[test]
    y = scores[[index[i] for i in ids]]
    cancers = art["cancers"].astype(str)[test]
    tss, _ = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    ry = cross_fitted_residuals(y, design, seed=42)
    rx = cross_fitted_residuals(states, design, seed=42)
    rank_raw = effective_rank(states)
    rank_res = effective_rank(rx)
    channel = top_canonical_correlation(rx, ry, n_components=16)
    rows.append((path.split("/")[-1], rank_raw, rank_res, channel))
    print(f"{rows[-1][0]:>12}  rank_raw={rank_raw:7.3f}  rank_resid={rank_res:7.3f}  channel={channel:.4f}")

if len(rows) >= 2:
    for label, col in (("rank (raw)", 1), ("rank (residualised)", 2), ("channel", 3)):
        values = [r[col] for r in rows]
        lo, hi = min(values), max(values)
        print(f"\n{label:22s} min={lo:.4f}  max={hi:.4f}  spread={hi/max(lo,1e-12):.3f}x")
    print("\nThis envelope is a FLOOR: programme_only is the stable arm, and same-seed")
    print("repeats exclude seed variation entirely. Report it as a floor, not as THE envelope.")
