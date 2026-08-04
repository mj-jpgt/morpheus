"""Verify, rather than assert, what the 16-component readout actually sees of a 128-column block.

Claim under test: `top_canonical_correlation(x, y, n_components=k)` depends on the target block y
ONLY through y's top-k principal subspace. If true, then at k=16 the 112 remaining columns of a
128-column block never enter the statistic, and "capacity-matched at 128" is not what the T1.1
table compared.

Three probes, each using the imported `top_canonical_correlation`:
  (a) replace y by its own rank-k PCA reconstruction    -> must be identical
  (b) reparametrise that rank-k block by a random invertible k x k map -> must be identical
  (c) delete everything outside the top-k PCs and rescale each retained PC -> must be identical
and one negative control:
  (d) replace y by its rank-(k+16) reconstruction       -> must DIFFER, or the probe is inert.
"""
import json

import numpy as np
import pandas as pd

from morpheus.v2.research.rebase.d2_compare import _load, _targets
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import top_canonical_correlation

DATA = "/home/ubuntu/e0_run/data"
A = "/home/ubuntu/e0_run/d2_v3"
K = 16
SEED = 42

pbs_index, pbs_scores = _targets("%s/pbs_targets_k128_v2.npz" % DATA, None)
raw = _load("%s/d2_v3_s42/artifacts/d2_h_seed42.npz" % A)
test = raw["split"].astype(str) == "test"
ids = raw["patient_ids"].astype(str)[test]
cancers = raw["cancers"].astype(str)[test]
tss, _ = pooled_tissue_source_site(ids, min_site_count=10)
design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
x = cross_fitted_residuals(raw["wsi_biology"].astype(np.float64)[test], design, seed=SEED)
y = cross_fitted_residuals(pbs_scores[[pbs_index[p] for p in ids]], design, seed=SEED)

centre = y.mean(axis=0, keepdims=True)
u, s, vt = np.linalg.svd(y - centre, full_matrices=False)
rng = np.random.default_rng(0)
m = rng.normal(size=(K, K))

variants = {
    "original_128_columns": y,
    "a_rank16_reconstruction": centre + (u[:, :K] * s[:K]) @ vt[:K],
    "b_rank16_random_invertible_reparam": (u[:, :K] * s[:K]) @ m,
    "c_top16_PCs_rescaled": u[:, :K] * (1.0 + rng.random(K)),
    # A rank-32 reconstruction PRESERVES the top-16 PC subspace, so it is invariant too --
    # it is not a control at all. The controls that actually perturb the selected subspace are
    # these two: swap in the NEXT 16 PCs, and rotate the block so different directions become
    # its leading ones.
    "d_rank32_reconstruction_NOT_A_CONTROL": centre + (u[:, :2 * K] * s[:2 * K]) @ vt[:2 * K],
    "e_PCs_17_to_32_only_CONTROL": u[:, K:2 * K] * s[K:2 * K],
    "f_spectrum_reversed_CONTROL": u * s[::-1],
}
out = {"k": K, "n_test": int(test.sum()), "cond_of_random_map": float(np.linalg.cond(m)),
       "values": {}}
base = None
for name, block in variants.items():
    value = top_canonical_correlation(x, block, n_components=K)
    out["values"][name] = value
    if base is None:
        base = value
    print("%-42s top-CCA(k=16) = %.12f   delta_vs_original = %+.2e"
          % (name, value, value - base), flush=True)

invariant = [n for n in ("a_rank16_reconstruction", "b_rank16_random_invertible_reparam",
                         "c_top16_PCs_rescaled")
             if abs(out["values"][n] - base) < 1e-10]
out["invariant_probes_passing"] = invariant
out["controls"] = {n: bool(abs(out["values"][n] - base) > 1e-6)
                   for n in ("e_PCs_17_to_32_only_CONTROL", "f_spectrum_reversed_CONTROL")}
out["negative_control_differs"] = all(out["controls"].values())
print("")
print("invariant under %d/3 subspace-preserving transforms; both subspace-CHANGING controls "
      "differ: %s  %s" % (len(invariant), out["negative_control_differs"], out["controls"]))
json.dump(out, open("/home/ubuntu/ws_d2sym/out/SUBSPACE_CHECK.json", "w"), indent=2)
