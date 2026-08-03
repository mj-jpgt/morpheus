"""Recompute every surviving effective-rank instance under the ONE canonical definition.

CPU only, thread-capped. Reads frozen .npz artifacts; touches no GPU and no
training state. Reproduces each instance's ORIGINAL computation first (so the
recomputation can be checked against the published number) and then reports every
variant of `RANK_VARIANTS` side by side.

usage: OMP_NUM_THREADS=1 ... python3 recompute_rank.py <out.json>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/ubuntu/ws_rank")
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import CANONICAL, RANK_VARIANTS, effective_rank

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "RANK_RECOMPUTE.json")
E0 = Path("/home/ubuntu/e0_run")
VARIANTS = list(RANK_VARIANTS)


def all_variants(x):
    return {name: effective_rank(x, variant=RANK_VARIANTS[name]) for name in VARIANTS}


def legacy_r1_absolute_cut(x):
    """The pre-2026-08-05 `spectral.py`: centred, order 1, ABSOLUTE 1e-12 cut.

    Kept here so the recomputation can prove that replacing the absolute cut with a
    relative one did not move any historical value.
    """
    v = np.asarray(x, dtype=np.float64)
    v = v - v.mean(axis=0, keepdims=True)
    s = np.linalg.svd(v, compute_uv=False)
    s = s[s > 1e-12]
    if s.size == 0:
        return 0.0
    p = s / s.sum()
    return float(np.exp(-(p * np.log(p)).sum()))


def residualise(x, ids, cancers, seed=42):
    """d2_readout.py's residualisation, verbatim."""
    tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    return cross_fitted_residuals(x, design, seed=seed), len(sites)


def load_block(path, key="wsi_biology", split="test"):
    raw = np.load(path, allow_pickle=True)
    keep = np.asarray(raw["split"]).astype(str) == split
    return (np.asarray(raw[key], dtype=np.float64)[keep],
            np.asarray(raw["patient_ids"]).astype(str)[keep],
            np.asarray(raw["cancers"]).astype(str)[keep])


results = {"canonical": CANONICAL.label, "variant_order": VARIANTS, "instances": {}}


def record(group, label, x, *, residualised_too=False, ids=None, cancers=None, extra=None):
    entry = {"n_rows": int(x.shape[0]), "n_features": int(x.shape[1]),
             "raw": all_variants(x), "raw_legacy_R1_absolute_cut": legacy_r1_absolute_cut(x)}
    if residualised_too:
        rx, n_sites = residualise(x, ids, cancers)
        entry["residualised"] = all_variants(rx)
        entry["residualised_legacy_R1_absolute_cut"] = legacy_r1_absolute_cut(rx)
        entry["n_sites_kept"] = n_sites
    if extra:
        entry.update(extra)
    results["instances"].setdefault(group, {})[label] = entry
    print(f"  {group:<22} {label:<26} R1={entry['raw']['R1']:9.4f}"
          + (f" R1res={entry['residualised']['R1']:9.4f}" if residualised_too else ""), flush=True)


# ---------------------------------------------------------------- instance 6 (D2)
print("instance 6 -- D2 supervision-target ablation, 3 seeds x 2 arms", flush=True)
D2 = {
    "H42": E0 / "d2_v3/d2_v3_s42/artifacts/d2_h_seed42.npz",
    "I42": E0 / "d2_v3/d2_v3_s42/artifacts/d2_i_seed42.npz",
    "H43": E0 / "d2_v3/d2_v3_s43/artifacts/d2_h_seed43.npz",
    "I43": E0 / "d2_v3/d2_v3_s43/artifacts/d2_i_seed43.npz",
    "H44": E0 / "d2_v3/d2_v3_s44/artifacts/d2_h_seed44.npz",
    "I44": E0 / "d2_v3/d2_v3_s44/artifacts/d2_i_seed44.npz",
    "H42_recovered_reexport": E0 / "d2_v3/recovered_artifacts/d2_h_seed42.npz",
    "I42_recovered_epoch13": E0 / "d2_v3/recovered_artifacts/d2_i_seed42_EPOCH13.npz",
}
for label, path in D2.items():
    if not path.exists():
        results["instances"].setdefault("instance6_D2", {})[label] = {"status": "ARTIFACT MISSING"}
        continue
    x, ids, cancers = load_block(path)
    record("instance6_D2", label, x, residualised_too=True, ids=ids, cancers=cancers,
           extra={"path": str(path)})

# ---------------------------------------------------------------- instance 4 (dilution)
print("instance 4 -- foreign-tumour dilution, 7 levels, zero fitted parameters", flush=True)
DIL = Path("/home/ubuntu/p1_out/dilution/dilution_foreign_tumour_pca256.npz")
if DIL.exists():
    raw = np.load(DIL, allow_pickle=True)
    test = np.asarray(raw["split"]).astype(str) == "test"
    ids = np.asarray(raw["patient_ids"]).astype(str)[test]
    cancers = np.asarray(raw["cancers"]).astype(str)[test]
    for key in [k for k in raw.files if k.startswith("wsi_dilution_")]:
        x = np.asarray(raw[key], dtype=np.float64)[test]
        record("instance4_dilution", key, x, residualised_too=True, ids=ids, cancers=cancers,
               extra={"path": str(DIL)})
else:
    results["instances"]["instance4_dilution"] = {"status": "ARTIFACT MISSING"}

# ---------------------------------------------------------------- instance 5 (D1-B)
print("instance 5 -- D1-B programme_only vs programme_free, 3 seeds", flush=True)
for label in ("d1_p_seed42", "d1_p_seed43", "d1_p_seed44",
              "d1_f_seed42", "d1_f_seed43", "d1_f_seed44"):
    path = E0 / f"d1_v2/artifacts/{label}.npz"
    if not path.exists():
        results["instances"].setdefault("instance5_D1B", {})[label] = {"status": "ARTIFACT MISSING"}
        continue
    x, ids, cancers = load_block(path)
    record("instance5_D1B", label, x, residualised_too=True, ids=ids, cancers=cancers,
           extra={"path": str(path)})

# --------------------------------------------- instability: in-training R3 tripwire
print("instability -- in-run rank tripwire at global step 200 (statistic R3)", flush=True)
trip = {}
for label in ("d1_p_seed42", "d1_p_seed43", "d1_p_seed44",
              "d1_f_seed42", "d1_f_seed43", "d1_f_seed44"):
    metrics = E0 / f"d1_v2/{label}/train_metrics.jsonl"
    if not metrics.exists():
        continue
    for line in metrics.read_text().splitlines():
        row = json.loads(line)
        if "train_rank_tripwire_observed" in row:
            trip[label] = {"epoch": row["epoch"], "R3_at_step_200": row["train_rank_tripwire_observed"],
                           "source": str(metrics)}
results["instances"]["instability_tripwire_step200_R3"] = trip
for k, v in sorted(trip.items()):
    print(f"  {k:<14} R3@200 = {v['R3_at_step_200']:.3f}", flush=True)

OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
print(f"\nwritten: {OUT}", flush=True)
