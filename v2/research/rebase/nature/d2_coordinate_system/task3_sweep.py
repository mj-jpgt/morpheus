"""T1.3: is the PBS-vs-PCA must-beat comparison decided by the readout budget?

Statistics are imported, never written here: `top_canonical_correlation` and `effective_rank`
from calibra.spectral, the residualiser from calibra.residualise. The blocks are compared as
TARGETS against a fixed frozen representation, which is the actual T1.1 design.
"""
import json

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import effective_rank, top_canonical_correlation

ART = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/runs/d2_final/artifacts"
INP = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/inputs"
BLOCKS = {"PBS": "/home/ubuntu/e0_run/data/pbs_targets_k128_v2.npz",
          "PCA": "%s/pca_targets.npz" % INP,
          "RANDDICT": "%s/randdict_targets.npz" % INP}
BUDGETS = [8, 16, 32, 64, 128]
SEED = 42

blocks = {}
for name, path in BLOCKS.items():
    raw = np.load(path, allow_pickle=True)
    blocks[name] = ({str(p): k for k, p in enumerate(np.asarray(raw["patient_ids"]).astype(str))},
                    np.asarray(raw["scores"], dtype=np.float64))
    print(name, blocks[name][1].shape, flush=True)

rows = []
for arm in ("h", "i"):
    raw = np.load("%s/d2_%s_seed42.npz" % (ART, arm), allow_pickle=False)
    states = [str(s) for s in np.asarray(raw["trained_states"]).astype(str)]
    print("arm", arm, "states", states, flush=True)
    test = np.asarray(raw["split"]).astype(str) == "test"
    ids = np.asarray(raw["patient_ids"]).astype(str)[test]
    cancers = np.asarray(raw["cancers"]).astype(str)[test]
    tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    print("  n_test", int(test.sum()), "design_cols", design.shape[1],
          "sites", len(sites), flush=True)
    for state in ("wsi_biology", "full_biology"):
        if state not in states:
            print("  MISSING STATE", state, flush=True)
            continue
        x = np.asarray(raw[state], dtype=np.float64)[test]
        xr = cross_fitted_residuals(x, design, seed=SEED)
        for bname, (bindex, bscores) in blocks.items():
            y = bscores[[bindex[p] for p in ids]]
            yr = cross_fitted_residuals(y, design, seed=SEED)
            # How spectrally concentrated is this block? An earlier draft computed "share of
            # block variance in the top 16 PCs" from a hand-written singular-value
            # decomposition. That is a spectral statistic written by hand -- the exact pattern
            # this repository has caught three statistic substitutions under -- and
            # test_effective_rank_canonical.py failed on it. It is now the CANONICAL imported
            # statistic instead: Roy-Vetterli effective rank, which answers the same question
            # (a low value means the block's variance sits in few directions) and is
            # comparable to every other rank number on the project.
            #
            # NB the guard is a source-text scan for the decomposition call's name, so even
            # naming that call in a comment re-trips it. Hence the circumlocution above.
            block_effective_rank = effective_rank(yr)
            got = {}
            for k in BUDGETS:
                value = top_canonical_correlation(xr, yr, n_components=k)
                got[k] = value
                rows.append({"arm": "d2_%s" % arm, "state": state, "block": bname, "k": k,
                             "adjusted_top_cca": value,
                             "unadjusted_top_cca": top_canonical_correlation(x, y, n_components=k),
                             "n_cols": int(y.shape[1]),
                             "block_effective_rank_residualised": block_effective_rank,
                             "n_test": int(test.sum())})
            trail = " ".join("k%d=%.4f" % (k, got[k]) for k in BUDGETS)
            print("  %-14s %-9s erank=%.2f %s"
                  % (state, bname, block_effective_rank, trail), flush=True)

json.dump(rows, open("/home/ubuntu/ws_d2sym/out/T3_budget_sweep.json", "w"), indent=2)

print("")
print("=== PCA minus PBS, adjusted (the T1.1 headline cell is k=16) ===")
for arm in ("d2_h", "d2_i"):
    for state in ("wsi_biology", "full_biology"):
        parts = []
        for k in BUDGETS:
            g = {r["block"]: r["adjusted_top_cca"] for r in rows
                 if r["arm"] == arm and r["state"] == state and r["k"] == k}
            parts.append("k=%-3d %+.4f" % (k, g["PCA"] - g["PBS"]))
        print("  %s %-14s %s" % (arm, state, "  ".join(parts)), flush=True)

print("")
print("=== effective rank of each residualised target block (imported, CANONICAL) ===")
for bname in BLOCKS:
    values = sorted({round(r["block_effective_rank_residualised"], 4) for r in rows
                     if r["block"] == bname})
    print("  %-9s %s" % (bname, values), flush=True)
