"""T1.3: is the PBS-vs-PCA must-beat comparison decided by the readout budget?

Statistics are imported, never written here: `top_canonical_correlation` and
`cca_spectrum` from calibra.spectral, the residualiser from calibra.residualise.
The blocks are compared as TARGETS against a fixed frozen representation, which is
the actual T1.1 design.
"""
import json
import numpy as np, pandas as pd
from morpheus.v2.calibra.residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from morpheus.v2.calibra.spectral import top_canonical_correlation, cca_spectrum

ART = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/runs/d2_final/artifacts"
INP = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/inputs"
BLOCKS = {"PBS": "/home/ubuntu/e0_run/data/pbs_targets_k128_v2.npz",
          "PCA": f"{INP}/pca_targets.npz",
          "RANDDICT": f"{INP}/randdict_targets.npz"}
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
    raw = np.load(f"{ART}/d2_{arm}_seed42.npz", allow_pickle=False)
    states = [s for s in np.asarray(raw["trained_states"]).astype(str)]
    print("arm", arm, "states", states, flush=True)
    test = np.asarray(raw["split"]).astype(str) == "test"
    ids = np.asarray(raw["patient_ids"]).astype(str)[test]
    cancers = np.asarray(raw["cancers"]).astype(str)[test]
    tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    print("  n_test", int(test.sum()), "design_cols", design.shape[1], "sites", len(sites), flush=True)
    for state in ("wsi_biology", "full_biology"):
        if state not in states:
            print("  MISSING STATE", state, flush=True); continue
        x = np.asarray(raw[state], dtype=np.float64)[test]
        xr = cross_fitted_residuals(x, design, seed=SEED)
        for bname, (bindex, bscores) in blocks.items():
            y = bscores[[bindex[p] for p in ids]]
            yr = cross_fitted_residuals(y, design, seed=SEED)
            # share of the block variance the 16-component readout can even see
            s = np.linalg.svd(yr - yr.mean(0, keepdims=True), compute_uv=False)
            var = s**2
            share16 = float(var[:16].sum() / var.sum())
            got = {}
            for k in BUDGETS:
                value = top_canonical_correlation(xr, yr, n_components=k)
                got[k] = value
                rows.append({"arm": f"d2_{arm}", "state": state, "block": bname, "k": k,
                             "adjusted_top_cca": value,
                             "unadjusted_top_cca": top_canonical_correlation(x, y, n_components=k),
                             "n_cols": int(y.shape[1]), "top16_variance_share": share16,
                             "n_test": int(test.sum())})
            trail = " ".join("k%d=%.4f" % (k, got[k]) for k in BUDGETS)
            print("  %-14s %-9s share16=%.3f %s" % (state, bname, share16, trail), flush=True)
json.dump(rows, open("/home/ubuntu/ws_d2sym/out/T3_budget_sweep.json", "w"), indent=2)

print("\n=== PCA minus PBS, adjusted (the T1.1 headline cell is k=16) ===", flush=True)
for arm in ("d2_h", "d2_i"):
    for state in ("wsi_biology", "full_biology"):
        line = [f"{arm} {state:14s}"]
        for k in BUDGETS:
            g = {r["block"]: r["adjusted_top_cca"] for r in rows if r["arm"]==arm and r["state"]==state and r["k"]==k}
            if "PCA" in g and "PBS" in g:
                line.append(f"k={k:3d}: {g[PCA]-g[PBS]:+.4f}")
        print("  ".join(line), flush=True)
