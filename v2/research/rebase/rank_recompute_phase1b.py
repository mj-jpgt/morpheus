"""Instance 2 (Phase 1b) -- artifacts located, recomputed under the canonical definition."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "/home/ubuntu/ws_rank")
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank

V = list(RANK_VARIANTS)
ROOTS = {"v21_release_20260720_retry3_resume_safe":
         Path("/lambda/nfs/geeg/biorag3_persistent_20260711/runs/v21_release_20260720_retry3_resume_safe/artifacts"),
         "v22_a10_11v21_20260725":
         Path("/lambda/nfs/geeg/biorag3_persistent_20260711/runs/v22_a10_11v21_20260725/artifacts")}
RECORDED = {("diagnostic_full_seed42", "wsi_biology"): 38.483394849435676,
            ("diagnostic_programme_only_seed42", "wsi_biology"): 32.05939849533399,
            ("diagnostic_full_seed42", "wsi_identity"): 191.0723820209697,
            ("diagnostic_identity_only_seed42", "wsi_identity"): 191.07122558003888,
            ("diagnostic_full_seed42", "full_biology"): 47.25778175996386,
            ("diagnostic_programme_only_seed42", "full_biology"): 38.710319195785296,
            ("diagnostic_full_seed42", "rna_biology"): 32.583767009997516,
            ("diagnostic_programme_only_seed42", "rna_biology"): 30.3290322739212}
out = {}
for root_name, root in ROOTS.items():
    for artifact in sorted(root.glob("diagnostic_*.npz")):
        raw = np.load(artifact, allow_pickle=True)
        test = np.asarray(raw["split"]).astype(str) == "test"
        ids = np.asarray(raw["patient_ids"]).astype(str)[test]
        cancers = np.asarray(raw["cancers"]).astype(str)[test]
        tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
        design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
        for key in [k for k in raw.files if k.endswith(("_biology", "_identity", "_patient"))]:
            arr = np.asarray(raw[key], dtype=np.float64)
            if arr.ndim != 2 or len(arr) != len(test):
                continue
            x = arr[test]
            entry = {"n_rows": int(x.shape[0]), "n_features": int(x.shape[1]),
                     "raw": {n: effective_rank(x, variant=RANK_VARIANTS[n]) for n in V}}
            rx = cross_fitted_residuals(x, design, seed=42)
            entry["residualised"] = {n: effective_rank(rx, variant=RANK_VARIANTS[n]) for n in V}
            recorded = RECORDED.get((artifact.stem, key))
            if recorded is not None:
                entry["recorded_2026_07"] = recorded
                entry["reproduces_recorded"] = bool(abs(entry["raw"]["R1"] - recorded) < 1e-6)
                print(f"{root_name:<40}{artifact.stem:<36}{key:<16}"
                      f"R1={entry['raw']['R1']:10.4f}  recorded={recorded:10.4f}  "
                      f"match={entry['reproduces_recorded']}", flush=True)
            out.setdefault(root_name, {}).setdefault(artifact.stem, {})[key] = entry
Path("/home/ubuntu/ws_rank/RANK_RECOMPUTE_P1B.json").write_text(json.dumps(out, indent=2))
print("written", flush=True)
