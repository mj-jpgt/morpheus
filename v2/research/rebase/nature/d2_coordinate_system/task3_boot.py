"""Paired bootstrap on the PCA-vs-PBS target-block difference, at two readout budgets.

Design note. `paired_multivariate_patient_and_cancer_bootstrap` was written to hold the
TARGET fixed and compare two REPRESENTATIONS. Here the representation is fixed and the
two TARGET BLOCKS are compared, so the representation is passed as `actual` and the
blocks as `teacher`/`challenger`. That role swap is exact and not an approximation:
`cca_spectrum` whitens both sides symmetrically and takes the singular values of
`ux.T @ uy`, whose spectrum equals that of its transpose, so
`top_canonical_correlation(a, b) == top_canonical_correlation(b, a)`. The script asserts
this on the real matrices before using it.

Repeats = 400 to match the published T1.1 protocol (`baseline_paired_bootstrap.py --n-boot 400`).
Every statistic is imported from calibra.spectral / calibra.residualise.
"""
import json
import sys

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import top_canonical_correlation
from morpheus.v2.paired_bootstrap import paired_multivariate_patient_and_cancer_bootstrap

ART = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/runs/d2_final/artifacts"
INP = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/inputs"
BLOCKS = {"PBS": "/home/ubuntu/e0_run/data/pbs_targets_k128_v2.npz",
          "PCA": "%s/pca_targets.npz" % INP,
          "RANDDICT": "%s/randdict_targets.npz" % INP}
SEED = 42
REPEATS = 400

arm, state, budget = sys.argv[1], sys.argv[2], int(sys.argv[3])

blocks = {}
for name, path in BLOCKS.items():
    raw = np.load(path, allow_pickle=True)
    ids = np.asarray(raw["patient_ids"]).astype(str)
    blocks[name] = ({str(p): k for k, p in enumerate(ids)},
                    np.asarray(raw["scores"], dtype=np.float64))

raw = np.load("%s/d2_%s_seed42.npz" % (ART, arm), allow_pickle=False)
test = np.asarray(raw["split"]).astype(str) == "test"
ids = np.asarray(raw["patient_ids"]).astype(str)[test]
cancers = np.asarray(raw["cancers"]).astype(str)[test]
tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
x = np.asarray(raw[state], dtype=np.float64)[test]
xr = cross_fitted_residuals(x, design, seed=SEED)

resid = {}
for name, (bindex, bscores) in blocks.items():
    y = bscores[[bindex[p] for p in ids]]
    resid[name] = cross_fitted_residuals(y, design, seed=SEED)

# The role swap this script depends on, checked on the real matrices.
forward = top_canonical_correlation(xr, resid["PBS"], n_components=budget)
reverse = top_canonical_correlation(resid["PBS"], xr, n_components=budget)
assert abs(forward - reverse) < 1e-10, (forward, reverse)


def metric(actual, representation):
    return top_canonical_correlation(actual, representation, n_components=budget)


out = {"arm": "d2_%s" % arm, "state": state, "k": budget, "repeats": REPEATS,
       "n_test": int(test.sum()), "n_design_columns": int(design.shape[1]),
       "n_sites_kept": len(sites), "seed": SEED,
       "symmetry_check": {"forward": forward, "reverse": reverse},
       "points": {n: metric(xr, resid[n]) for n in BLOCKS},
       "contrasts": {}}
for teacher, challenger in (("PBS", "PCA"), ("RANDDICT", "PBS")):
    out["contrasts"]["%s_minus_%s" % (challenger, teacher)] = \
        paired_multivariate_patient_and_cancer_bootstrap(
            metric, xr, resid[teacher], resid[challenger], cancers,
            repeats=REPEATS, seed=SEED)

path = "/home/ubuntu/ws_d2sym/out/T3_boot_d2%s_%s_k%d.json" % (arm, state, budget)
json.dump(out, open(path, "w"), indent=2)
for name, payload in out["contrasts"].items():
    for mode in ("patient", "cancer"):
        p = payload[mode]
        print("d2_%s %s k=%d %-16s %-8s delta=%+.4f CI[%+.4f,%+.4f] p_improve=%.4f"
              % (arm, state, budget, name, mode, p["point_delta"], p["ci95_low"],
                 p["ci95_high"], p["p_improve"]), flush=True)
