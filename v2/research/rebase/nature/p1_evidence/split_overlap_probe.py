"""Are the eight partitions eight different partitions? (5.1 of the predeclaration.)

Defines no statistic. Imports ``exposure_split`` and ``_load_block`` unchanged and reports the
pairwise Jaccard overlap of the eight exposure patient sets, plus each split fold size and each
fold cancer count. A pair above 0.9 voids the multi-split run.
"""
import itertools, json, sys
from pathlib import Path
import numpy as np
from morpheus.v2.calibra.nonlinear_adjustment import _load_block
from morpheus.v2.research.rebase.nature.p4_certification.p4_certify import exposure_split

P = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d"
SEEDS = [42, 7, 11, 23, 101, 555, 2718, 31337]
block = _load_block(Path(f"{P}/runs/d2_final/artifacts/d2_h_seed42.npz"),
                    Path(f"{P}/data/frozen_rna_targets.npz"), "wsi_biology", "test", 10)
cancers, ids = block["cancers"], block["patient_ids"]
sets, rows = {}, []
for s in SEEDS:
    d = exposure_split(cancers, discovery_fraction=0.5, seed=s)
    sets[s] = set(ids[~d].tolist())
    rows.append({"split_seed": s, "n_discovery": int(d.sum()), "n_exposure": int((~d).sum()),
                 "n_cancers_discovery": len(set(cancers[d].tolist())),
                 "n_cancers_exposure": len(set(cancers[~d].tolist()))})
pairs = []
for a, b in itertools.combinations(SEEDS, 2):
    j = len(sets[a] & sets[b]) / len(sets[a] | sets[b])
    pairs.append({"a": a, "b": b, "jaccard": round(j, 4)})
js = [p["jaccard"] for p in pairs]
out = {"folds": rows, "pairs": pairs,
       "jaccard_min": min(js), "jaccard_max": max(js), "jaccard_median": float(np.median(js)),
       "n_distinct_partitions": len({frozenset(v) for v in sets.values()})}
print(json.dumps(out, indent=1))
Path(f"{P}/p1_evidence/split_stability/split_overlap_probe.json").write_text(json.dumps(out, indent=1))
