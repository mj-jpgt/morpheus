"""How different ARE the two exams? Geometry of the gene-set block vs the PBS-code block.

If the two target spaces were near-identical, a +0.12 shift in the arm contrast from swapping
them would be inexplicable. If they are largely distinct spaces, the shift is what a
coordinate-system confound looks like. This quantifies which.

`cca_spectrum` and `effective_rank` are imported from calibra.spectral; the residualiser from
calibra.residualise. Nothing is computed by hand.
"""
import json

import numpy as np
import pandas as pd

from morpheus.v2.research.rebase.d2_compare import _load, _targets
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import cca_spectrum, effective_rank

DATA = "/home/ubuntu/e0_run/data"
A = "/home/ubuntu/e0_run/d2_v3"
SEED = 42

gs_index, gs_scores = _targets("%s/frozen_rna_targets.npz" % DATA,
                               ["heldout_pathway", "immune_tme", "tumour_state"])
pbs_index, pbs_scores = _targets("%s/pbs_targets_k128_v2.npz" % DATA, None)

raw = _load("%s/d2_v3_s42/artifacts/d2_h_seed42.npz" % A)
test = raw["split"].astype(str) == "test"
ids = raw["patient_ids"].astype(str)[test]
cancers = raw["cancers"].astype(str)[test]
tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])

gs = gs_scores[[gs_index[p] for p in ids]]
pbs = pbs_scores[[pbs_index[p] for p in ids]]
out = {"n_test": int(test.sum()), "n_design_columns": int(design.shape[1]),
       "n_sites_kept": len(sites), "blocks": {}}
for block, (g, p) in {"raw": (gs, pbs),
                      "residualised": (cross_fitted_residuals(gs, design, seed=SEED),
                                       cross_fitted_residuals(pbs, design, seed=SEED))}.items():
    spectrum = cca_spectrum(g, p, n_components=16)
    out["blocks"][block] = {
        "geneset40_effective_rank": effective_rank(g),
        "pbs128_effective_rank": effective_rank(p),
        "cca_spectrum_k16": [float(v) for v in spectrum],
        "top_cca": float(spectrum[0]),
        "mean_cca": float(np.mean(spectrum)),
        "n_canonical_above_0p9": int((spectrum > 0.9).sum()),
        "n_canonical_above_0p5": int((spectrum > 0.5).sum()),
    }
    b = out["blocks"][block]
    print("%-13s erank(geneset40)=%.2f erank(pbs128)=%.2f  top-CCA=%.4f mean=%.4f  "
          ">0.9: %d/16  >0.5: %d/16" % (block, b["geneset40_effective_rank"],
          b["pbs128_effective_rank"], b["top_cca"], b["mean_cca"],
          b["n_canonical_above_0p9"], b["n_canonical_above_0p5"]), flush=True)
    print("   spectrum: " + " ".join("%.3f" % v for v in spectrum), flush=True)

json.dump(out, open("/home/ubuntu/ws_d2sym/out/EXAM_GEOMETRY.json", "w"), indent=2)
