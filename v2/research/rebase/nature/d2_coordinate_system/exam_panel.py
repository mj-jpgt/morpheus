"""Is arm I's parity on the PBS codes SPECIFIC to its own supervision, or generic?

Test 1 found the arms tie on the 128 PBS codes while arm H wins by ~0.12 on gene sets. Two
readings are still open:
  (i)  arm I's training bought it exactly enough to draw on its own coordinates -- specific; or
  (ii) both arms tie on ANY expression-derived 128-column block, and the gene-set exam is the
       outlier rather than the PBS exam.
(ii) is testable with blocks already on disk that neither arm trained on: the PCA basis, the
size- and spectrum-matched random dictionary, and gene-label-shuffled PBS. If the arms tie on
those too, the PBS tie says nothing about arm I's supervision.

Point estimates only (the 2,000-repeat bootstraps are budgeted for the two predeclared exams).
d2_compare's own loaders, calibra's residualiser, calibra's statistic. Nothing written here.
"""
import json

import numpy as np
import pandas as pd

from morpheus.v2.research.rebase.d2_compare import _load, _targets
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import top_canonical_correlation

DATA = "/home/ubuntu/e0_run/data"
INP = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/inputs"
A = "/home/ubuntu/e0_run/d2_v3"
SEED = 42

EXAMS = [
    ("geneset_untrained40", "%s/frozen_rna_targets.npz" % DATA,
     ["heldout_pathway", "immune_tme", "tumour_state"]),
    ("pbs_codes128_ARM_I_OWN", "%s/pbs_targets_k128_v2.npz" % DATA, None),
    ("pbs_shuffled_s1", "%s/pbs_shuffled_seed1.npz" % INP, None),
    ("pca_basis128", "%s/pca_targets.npz" % INP, None),
    ("random_dictionary128", "%s/randdict_targets.npz" % INP, None),
    ("geneset_random_control", "%s/frozen_rna_targets.npz" % DATA, ["random_control"]),
]

rows = []
for name, path, groups in EXAMS:
    index, scores = _targets(path, groups)
    for seed in (42, 43, 44):
        h = _load("%s/d2_v3_s%d/artifacts/d2_h_seed%d.npz" % (A, seed, seed))
        i = _load("%s/d2_v3_s%d/artifacts/d2_i_seed%d.npz" % (A, seed, seed))
        test = h["split"].astype(str) == "test"
        ids = h["patient_ids"].astype(str)[test]
        cancers = h["cancers"].astype(str)[test]
        y = scores[[index[p] for p in ids]]
        hx = h["wsi_biology"].astype(np.float64)[test]
        ix = i["wsi_biology"].astype(np.float64)[test]
        tss, _ = pooled_tissue_source_site(ids, min_site_count=10)
        design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
        yy = cross_fitted_residuals(y, design, seed=SEED)
        hh = cross_fitted_residuals(hx, design, seed=SEED)
        ii = cross_fitted_residuals(ix, design, seed=SEED)
        ph = top_canonical_correlation(yy, hh, n_components=16)
        pi = top_canonical_correlation(yy, ii, n_components=16)
        rows.append({"exam": name, "seed": seed, "block": "residualised",
                     "n_targets": int(y.shape[1]), "point_hallmark": ph, "point_pbs": pi,
                     "pbs_minus_hallmark": pi - ph})

json.dump(rows, open("/home/ubuntu/ws_d2sym/out/EXAM_PANEL.json", "w"), indent=2)
print("residualised block, top-CCA k=16, wsi_biology, n_test=2766")
print("%-24s %8s %8s %8s   %s" % ("exam", "s42", "s43", "s44", "mean"))
for name, _, _ in EXAMS:
    d = [r["pbs_minus_hallmark"] for s in (42, 43, 44)
         for r in rows if r["exam"] == name and r["seed"] == s]
    print("%-24s %+8.4f %+8.4f %+8.4f   %+.4f" % (name, d[0], d[1], d[2], float(np.mean(d))))
