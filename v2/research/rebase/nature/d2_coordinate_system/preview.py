"""Point estimates only, through d2_compare own loaders and calibra own statistic.

No statistic is written here: `_targets`/`_load` are d2_compare, the residualiser is
calibra.residualise, the channel number is spectral.top_canonical_correlation. This is
d2_compare minus the bootstrap, so the point estimates must agree with the full run.
"""
import json, sys
import numpy as np, pandas as pd
from morpheus.v2.research.rebase.d2_compare import _load, _targets
from morpheus.v2.calibra.residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from morpheus.v2.calibra.spectral import top_canonical_correlation

A = "/home/ubuntu/e0_run/d2_v3"
DATA = "/home/ubuntu/e0_run/data"
EXAMS = {
    "geneset_untrained40": (f"{DATA}/frozen_rna_targets.npz", ["heldout_pathway","immune_tme","tumour_state"]),
    "geneset_all90":       (f"{DATA}/frozen_rna_targets.npz", None),
    "geneset_random_control": (f"{DATA}/frozen_rna_targets.npz", ["random_control"]),
    "pbs_codes128":        (f"{DATA}/pbs_targets_k128_v2.npz", None),
}
out = []
for exam,(tpath,groups) in EXAMS.items():
    tindex, tscores = _targets(tpath, groups)
    for seed in (42,43,44):
        h = _load(f"{A}/d2_v3_s{seed}/artifacts/d2_h_seed{seed}.npz")
        i = _load(f"{A}/d2_v3_s{seed}/artifacts/d2_i_seed{seed}.npz")
        test = h["split"].astype(str)=="test"; ids = h["patient_ids"].astype(str)[test]
        rows = np.asarray([tindex[x] for x in ids])
        y = tscores[rows]
        hx = h["wsi_biology"].astype(np.float64)[test]; ix = i["wsi_biology"].astype(np.float64)[test]
        cancers = h["cancers"].astype(str)[test]
        tss,_ = pooled_tissue_source_site(ids, min_site_count=10)
        design = confound_design(pd.DataFrame({"cancer":cancers,"tss":tss}), ["cancer","tss"])
        # RAW block: no residualisation. RESIDUALISED block: exactly d2_compare design.
        ry = cross_fitted_residuals(y, design, seed=42)
        rh = cross_fitted_residuals(hx, design, seed=42)
        ri = cross_fitted_residuals(ix, design, seed=42)
        for block,(yy,hh,ii) in {"raw":(y,hx,ix), "residualised":(ry,rh,ri)}.items():
            ph = top_canonical_correlation(yy,hh,n_components=16)
            pi = top_canonical_correlation(yy,ii,n_components=16)
            out.append({"exam":exam,"seed":seed,"block":block,"n_targets":int(y.shape[1]),
                        "n_test":int(test.sum()),"point_hallmark":ph,"point_pbs":pi,
                        "pbs_minus_hallmark":pi-ph})
            print(f"{exam:24s} s{seed} {block:13s} H={ph:.4f} I={pi:.4f} d={pi-ph:+.4f}", flush=True)
json.dump(out, open("/home/ubuntu/ws_d2sym/out/PREVIEW_points.json","w"), indent=2)
