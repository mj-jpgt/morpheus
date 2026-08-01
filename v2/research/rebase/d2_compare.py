"""Paired D2 H-vs-PBS comparison on frozen held-out patient rows.

This is deliberately a post-export analysis: it reads the exact representation
artifacts that CALIBRA received and recomputes the adjusted molecular channel
under paired patient and cancer-cluster resampling.  It never selects a
checkpoint or fits a representation using test rows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from morpheus.v2.calibra.spectral import top_canonical_correlation
from morpheus.v2.paired_bootstrap import paired_multivariate_patient_and_cancer_bootstrap


def _load(path: str):
    raw = np.load(path, allow_pickle=False)
    required = {"patient_ids", "split", "cancers", "trained_states", "wsi_biology"}
    if not required.issubset(raw.files) or "wsi_biology" not in set(raw["trained_states"].astype(str)):
        raise ValueError(f"D2 comparison needs a declared wsi_biology artifact: {path}")
    return raw


def _targets(path: str) -> tuple[dict[str, int], np.ndarray]:
    raw = np.load(path, allow_pickle=False)
    ids = np.asarray(raw["patient_ids"]).astype(str); values = np.asarray(raw["scores"], dtype=np.float64)
    names = np.asarray(raw["target_names"]).astype(str)
    keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
    values, names = values[:, keep], names[keep]
    if values.ndim != 2 or values.shape != (len(ids), len(names)) or len(set(ids)) != len(ids) or len(set(names)) != len(names):
        raise ValueError("frozen RNA target table fails D2 ID/target uniqueness requirements")
    return {identifier: index for index, identifier in enumerate(ids)}, values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hallmark-artifacts", nargs="+", required=True)
    parser.add_argument("--pbs-artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True); parser.add_argument("--output", required=True)
    parser.add_argument("--repeats", type=int, default=2000); parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if len(args.hallmark_artifacts) != len(args.pbs_artifacts):
        raise ValueError("D2 requires one seed-matched Hallmark/PBS artifact pair")
    target_index, target_scores = _targets(args.targets)
    result: dict[str, object] = {"pairs": []}
    for pair_index, (h_path, i_path) in enumerate(zip(args.hallmark_artifacts, args.pbs_artifacts)):
        h, i = _load(h_path), _load(i_path)
        for field in ("patient_ids", "split", "cancers"):
            if not np.array_equal(h[field].astype(str), i[field].astype(str)):
                raise ValueError(f"H/PBS frozen artifacts disagree on {field}")
        test = h["split"].astype(str) == "test"; ids = h["patient_ids"].astype(str)[test]
        rows = np.asarray([target_index.get(identifier, -1) for identifier in ids])
        if (rows < 0).any():
            raise ValueError("frozen RNA target artifact does not cover a D2 test patient")
        y, hx, ix = target_scores[rows], h["wsi_biology"].astype(np.float64)[test], i["wsi_biology"].astype(np.float64)[test]
        if not np.isfinite(y).all() or np.any(np.std(y, axis=0) < 1e-8):
            raise ValueError("D2 evaluated RNA targets are non-finite or constant")
        cancers = h["cancers"].astype(str)[test]
        tss, frequent_sites = pooled_tissue_source_site(ids, min_site_count=10)
        def paired_metric(actual: np.ndarray, representation: np.ndarray) -> float:
            # labels are selected by the bootstrap function only implicitly, so
            # this function is reserved for the IID patient bootstrap below.
            return top_canonical_correlation(representation, actual, n_components=16)
        # Bootstrap the exact residualised matrices.  Confound residualisation
        # is fitted once on held-out rows for both arms identically; resampling
        # then quantifies the paired channel difference rather than re-selecting.
        design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
        residual_y = cross_fitted_residuals(y, design, seed=args.seed)
        residual_h = cross_fitted_residuals(hx, design, seed=args.seed)
        residual_i = cross_fitted_residuals(ix, design, seed=args.seed)
        comparison = paired_multivariate_patient_and_cancer_bootstrap(
            paired_metric, residual_y, residual_h, residual_i, cancers, repeats=args.repeats, seed=args.seed + pair_index)
        result["pairs"].append({"hallmark": str(Path(h_path).resolve()), "pbs": str(Path(i_path).resolve()),
                                "n_test": int(test.sum()), "point_h": paired_metric(residual_y, residual_h),
                                "point_i": paired_metric(residual_y, residual_i), "pbs_minus_h": comparison,
                                "confounds": {"columns": ["cancer", "tss"], "min_site_count": 10,
                                              "n_distinct_sites_kept": len(frequent_sites)}})
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(result, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
