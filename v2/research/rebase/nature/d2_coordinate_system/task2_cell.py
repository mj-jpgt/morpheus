"""Test 2, one cell: (endpoint, seed, block) -> both arms' risks + the PAIRED C-index contrast.

Why survival and not the pre-registered ER status. `tcga_clinical_covariates.parquet` labels
690 ER / 688 PR patients and every one is BRCA, a DEVELOPMENT cancer in the maximal split:
585 train / 105 val / **0 test**. The D2 contrast lives entirely on the test partition, so the
pre-registered clinical control has zero coverage on it. Overall survival from the canonical
CDR builder covers 2,765 of the 2,766 test patients and is in neither arm's coordinate system.

Instruments, all imported, nothing statistical written here:
  `survival_evaluation.evaluate_coxnet_endpoint`  -- development-only Coxnet fit, cancer-grouped
      inner CV for alpha, risks exported on held-out test rows only.
  `survival_evaluation.paired_cindex_bootstrap`   -- paired Harrell C-index difference,
      patient and cancer-cluster bootstrap.
  `calibra.residualise`                           -- the same cancer + pooled-TSS design D2 uses.

Run one process per (seed, block); the two arms are handled inside so the paired contrast sees
identical patients.
"""
import json
import sys

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.survival_evaluation import evaluate_coxnet_endpoint, paired_cindex_bootstrap

A = "/home/ubuntu/e0_run/d2_v3"
OUTCOMES = ("/lambda/nfs/geeg/biorag3_persistent_20260711/runs/"
            "v22_a10_11v21_20260725/discovery_inputs/tcga_cdr_outcomes.parquet")
SEED = 42
BOOT = 400

endpoint, seed, block = sys.argv[1], int(sys.argv[2]), sys.argv[3]

outcomes = pd.read_parquet(OUTCOMES)
outcomes["patient_id"] = outcomes["patient_id"].astype(str)
outcomes = outcomes.drop_duplicates("patient_id")

cell = {"endpoint": endpoint, "seed": seed, "block": block, "arms": {}}
tables = {}
for arm in ("h", "i"):
    raw = np.load("%s/d2_v3_s%d/artifacts/d2_%s_seed%d.npz" % (A, seed, arm, seed),
                  allow_pickle=False)
    ids = np.asarray(raw["patient_ids"]).astype(str)
    cancers = np.asarray(raw["cancers"]).astype(str)
    split = np.asarray(raw["split"]).astype(str)
    features = np.asarray(raw["wsi_biology"], dtype=np.float64)
    if block == "residualised":
        tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
        design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}),
                                 ["cancer", "tss"])
        # The design uses no outcome information, so residualising every row leaks
        # nothing about survival; it is the same subspace removal D2 applies.
        features = cross_fitted_residuals(features, design, seed=SEED)
        cell["n_design_columns"] = int(design.shape[1])
        cell["n_sites_kept"] = len(sites)
    metadata = pd.DataFrame({"patient_id": ids, "cancer": cancers, "split": split})
    result = evaluate_coxnet_endpoint(
        features, metadata, outcomes, endpoint=endpoint,
        method="d2_%s_seed%d" % (arm, seed), representation_state="wsi_biology",
        bootstrap_repeats=BOOT, seed=SEED)
    row = {r["metric"]: r["value"] for r in result.rows}
    cell["arms"][arm] = {"status": result.status, "selected_alpha": result.selected_alpha,
                         **{k: result.rows[0].get(k) for k in
                            ("development_patients", "development_events",
                             "test_patients", "test_events")},
                         **row}
    tables[arm] = result.risks
    print("%s s%d %-13s arm_%s status=%s cindex=%s"
          % (endpoint, seed, block, arm, result.status,
             ("%.4f" % row["harrell_cindex"]) if "harrell_cindex" in row else "NA"), flush=True)

if all(cell["arms"][a]["status"] == "ok" for a in ("h", "i")):
    left, right = tables["h"], tables["i"]
    merged = left.merge(right, on="patient_id", suffixes=("_h", "_i"), validate="one_to_one")
    assert (merged.time_days_h.to_numpy() == merged.time_days_i.to_numpy()).all()
    assert (merged.event_h.to_numpy() == merged.event_i.to_numpy()).all()
    cell["n_paired"] = int(len(merged))
    cell["paired_i_minus_h"] = paired_cindex_bootstrap(
        merged.time_days_h.to_numpy(float), merged.event_h.to_numpy(bool),
        teacher_risk=merged.risk_h.to_numpy(float),
        challenger_risk=merged.risk_i.to_numpy(float),
        cancers=merged.cancer_h.to_numpy().astype(str), repeats=BOOT, seed=SEED)
    for mode in ("patient", "cancer"):
        p = cell["paired_i_minus_h"][mode]
        print("%s s%d %-13s PAIRED arm_I minus arm_H %-8s delta=%+.4f CI[%+.4f,%+.4f] p_improve=%.4f"
              % (endpoint, seed, block, mode, p["point_delta"], p["ci95_low"],
                 p["ci95_high"], p["p_improve"]), flush=True)

path = "/home/ubuntu/ws_d2sym/out/T2_%s_s%d_%s.json" % (endpoint, seed, block)
json.dump(cell, open(path, "w"), indent=2, default=str)
