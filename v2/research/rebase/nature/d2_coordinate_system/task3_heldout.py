"""T1.3(b): the T1.1 main table's statistic is heldout_top_cca, not adjusted_top_cca.

The published table is headed "Held-out top-CCA ... (n_components = 16)"; its numbers
(d2_h wsi_biology: PBS 0.5032, PCA 0.5520, randdict 0.4551) are heldout_top_cca, while
the block-level supporting table in the same entry quotes adjusted_top_cca
(PBS 0.5504, randdict 0.5226) -- which task3_sweep.py already reproduced exactly.

Both statistics are imported from calibra.spectral. Nothing is written here.
"""
import json
import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import heldout_top_cca, top_canonical_correlation

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
    ids = np.asarray(raw["patient_ids"]).astype(str)
    blocks[name] = ({str(p): k for k, p in enumerate(ids)},
                    np.asarray(raw["scores"], dtype=np.float64))

rows = []
for arm in ("h", "i"):
    raw = np.load("%s/d2_%s_seed42.npz" % (ART, arm), allow_pickle=False)
    test = np.asarray(raw["split"]).astype(str) == "test"
    ids = np.asarray(raw["patient_ids"]).astype(str)[test]
    cancers = np.asarray(raw["cancers"]).astype(str)[test]
    tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    for state in ("wsi_biology", "full_biology"):
        x = np.asarray(raw[state], dtype=np.float64)[test]
        xr = cross_fitted_residuals(x, design, seed=SEED)
        for bname, (bindex, bscores) in blocks.items():
            y = bscores[[bindex[p] for p in ids]]
            yr = cross_fitted_residuals(y, design, seed=SEED)
            got = {}
            for k in BUDGETS:
                value = heldout_top_cca(xr, yr, n_components=k, seed=SEED)
                got[k] = value
                rows.append({"arm": "d2_%s" % arm, "state": state, "block": bname, "k": k,
                             "heldout_top_cca": value,
                             "adjusted_top_cca": top_canonical_correlation(xr, yr, n_components=k)})
            trail = " ".join("k%d=%.4f" % (k, got[k]) for k in BUDGETS)
            print("d2_%s %-14s %-9s heldout: %s" % (arm, state, bname, trail), flush=True)

json.dump(rows, open("/home/ubuntu/ws_d2sym/out/T3_heldout_sweep.json", "w"), indent=2)

print("")
print("=== PCA minus PBS ===  (positive = PCA beats the interventional dictionary)")
for stat in ("heldout_top_cca", "adjusted_top_cca"):
    print("-- %s --" % stat)
    for arm in ("d2_h", "d2_i"):
        for state in ("wsi_biology", "full_biology"):
            parts = []
            for k in BUDGETS:
                g = {r["block"]: r[stat] for r in rows
                     if r["arm"] == arm and r["state"] == state and r["k"] == k}
                parts.append("k=%-3d %+.4f" % (k, g["PCA"] - g["PBS"]))
            print("   %s %-14s %s" % (arm, state, "  ".join(parts)))

print("")
print("=== PBS minus RANDDICT ===  (positive = the dictionary beats random directions)")
for stat in ("heldout_top_cca", "adjusted_top_cca"):
    print("-- %s --" % stat)
    for arm in ("d2_h", "d2_i"):
        for state in ("wsi_biology", "full_biology"):
            parts = []
            for k in BUDGETS:
                g = {r["block"]: r[stat] for r in rows
                     if r["arm"] == arm and r["state"] == state and r["k"] == k}
                parts.append("k=%-3d %+.4f" % (k, g["PBS"] - g["RANDDICT"]))
            print("   %s %-14s %s" % (arm, state, "  ".join(parts)))
