"""Where does the residual site signal in the inductively adjusted state live?

The operator can only site-adjust an exposure row whose site was frequent in the
discovery fold. This splits the exposure fold on exactly that line and runs the SAME
certificate on each side, for the inductive and the matched transductive arm. It defines
no statistic: certify_axes and prepare_state are imported unchanged.
"""
import json, sys
import numpy as np
from morpheus.v2.calibra.confound_certificate import certify_axes
from morpheus.v2.calibra.residualise import tissue_source_site
from morpheus.v2.research.rebase.nature.p4_certification.p4_certify import (exposure_split,
                                                                            load_state,
                                                                            prepare_state)

B = "/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d"
ART = B + "/runs/d2_final/artifacts/d2_h_seed42.npz"
block = load_state(ART, "wsi_biology", "test")

class A:
    min_site_count, n_permutations, seed, n_boot, n_boot_axes, n_jobs = 10, 1000, 42, 200, 8, 6

out = {}
for mode in ("inductive", "transductive_exposure"):
    state = prepare_state(block["features"], None, block["patient_ids"], block["cancers"],
                          adjustment=mode, discovery_fraction=0.5, seed=42, min_site_count=10)
    discovery = exposure_split(block["cancers"], discovery_fraction=0.5, seed=42)
    frequent_in_discovery = set(np.unique(tissue_source_site(block["patient_ids"][discovery]),
                                          return_counts=True)[0][
        np.unique(tissue_source_site(block["patient_ids"][discovery]), return_counts=True)[1] >= 10])
    site_of_row = tissue_source_site(state["patient_ids"])
    covered = np.asarray([s in frequent_in_discovery for s in site_of_row])
    for stratum, mask in (("all", np.ones(len(covered), bool)),
                          ("site_adjustable", covered), ("pooled_to_OTHER", ~covered)):
        for arm, matrix in (("raw", state["features"]),
                            ("adjusted", state["adjusted_features"])):
            result = certify_axes(matrix[mask], state["patient_ids"][mask],
                                  state["cancers"][mask], min_site_count=A.min_site_count,
                                  n_permutations=A.n_permutations, seed=A.seed, n_boot=A.n_boot,
                                  n_boot_axes=A.n_boot_axes, residualise=False, n_jobs=A.n_jobs)
            key = f"{mode}/{stratum}/{arm}"
            out[key] = {"n_rows": int(mask.sum()), "n_site_classes": int(result["n_classes"]),
                        "chance_rate": float(result["chance_rate"]),
                        "joint_lda": float(result["joint_lda_balanced_accuracy"]),
                        "joint_null_p95": float(result["joint_null_p95"]),
                        "joint_certified": bool(result["joint_certified"]),
                        "n_breaching_axes": int(result["n_breaching_axes"])}
            print(key, json.dumps(out[key]), flush=True)
json.dump(out, open("out/site_coverage_probe.json", "w"), indent=2)
print("PROBE_DONE")
