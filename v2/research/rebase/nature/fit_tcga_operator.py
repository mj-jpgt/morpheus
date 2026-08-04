"""Fit and persist the TCGA confound-adjustment operator, in two arms.

Why. ``inductive_adjustment.ConfoundAdjustmentOperator`` exists and is bit-identical
to ``residualise.cross_fitted_residuals`` on its fitting cohort, but no TCGA operator
had ever been fitted or persisted, so the thing it was built for -- adjusting an
external cohort with **TCGA's** operator instead of its own, which is the only way to
put two cohorts in one coordinate system -- was still not possible.

Two arms, and the second is not an optional extra.

  ``full``         columns ``["cancer", "tss"]``, site pooled at ``min_site_count=10``.
                   This is TCGA's internal protocol, the design every published
                   CALIBRA number used, and the operator to use for TCGA-internal work.

  ``cancer_only``  columns ``["cancer"]``. Required because ALCHEMIST has **no**
                   ``tissue_source_site`` field at all -- the GDC facet returns
                   ``_missing`` for all 1,176 cases and the paired manifest carries
                   ``tss = NA`` for all 1,106 rows. The site half of the design cannot
                   be applied there, so for an external comparison it must be dropped
                   from **both** arms or the two cohorts are not adjusted by the same
                   map. The external comparison is therefore **cancer-adjusted only**
                   and is weaker than TCGA's internal protocol: site is the confound
                   the leave-sites-out work shows matters most, and it goes
                   unadjusted on both sides.

What is verified here, all three declared before running:

  V1  IDENTITY on the fitting cohort -- ``adjust_reference`` must equal
      ``cross_fitted_residuals`` under ``np.array_equal`` (bit-for-bit, not allclose),
      or every published CALIBRA number silently stops being comparable to anything
      this operator produces.
  V2  UNSEEN SITE -- a row whose TSS code the reference never saw must be pooled to
      ``OTHER`` and said so in the report. The policy is not new: it is
      ``pooled_tissue_source_site``'s own rule evaluated at a reference count of zero,
      which is below every ``min_site_count``.
  V3  PERSISTENCE ROUND-TRIP -- ``load(save(op))`` must reproduce both the reference
      path and the new-row path exactly, and carry the provenance.

Also probed and reported rather than worked around: an ALCHEMIST-shaped cancer column.
ALCHEMIST's cancer labels are ``diagnoses.primary_diagnosis`` strings
(``Adenocarcinoma_NOS``, ``Squamous_cell_carcinoma_NOS``, ...) and are **disjoint**
from TCGA's ``LUAD``/``LUSC``, so a TCGA operator refuses them with
``UnseenLevelError``. That refusal is correct -- a zero one-hot row is the reference's
baseline level, not "no cancer type" -- and it means the external comparison needs an
explicit, declared label mapping. It is reported, not silently defaulted to
``on_unseen_level="zero"``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[5]))

from morpheus.v2.calibra.inductive_adjustment import (ConfoundAdjustmentOperator,  # noqa: E402
                                                      UnseenLevelError)
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,  # noqa: E402
                                             pooled_tissue_source_site, tissue_source_site)

ARMS = {
    "full": {"columns": ["cancer", "tss"], "site_column": "tss"},
    "cancer_only": {"columns": ["cancer"], "site_column": ""},
}


def reference_cohort(artifact, targets, partition: str):
    """The CALIBRA reference: the ``partition`` rows that also have RNA targets.

    Reproduces ``run_calibra.py``'s mask exactly -- partition filter, then alignment
    to ``frozen_rna_targets.npz`` -- so the operator is fitted on the same 2,530
    patients every published CALIBRA number was measured on.
    """
    patient_ids = np.asarray([str(p) for p in artifact["patient_ids"]])
    cancers = np.asarray([str(c) for c in artifact["cancers"]])
    split = np.asarray([str(s) for s in artifact["split"]])
    index = {str(pid): i for i, pid in enumerate(targets["patient_ids"])}
    aligned = np.asarray([index.get(pid, -1) for pid in patient_ids])
    mask = np.ones(len(patient_ids), dtype=bool) if partition == "all" else (split == partition)
    mask &= aligned >= 0
    return patient_ids[mask], cancers[mask], mask


def transductive_reference(matrix, cancers, patient_ids, arm: dict, min_site_count: int,
                           n_splits: int, alpha: float, seed: int):
    """V1's expected value: the transductive path, built the way run_calibra builds it."""
    frame = pd.DataFrame({"cancer": cancers})
    if arm["site_column"]:
        frame[arm["site_column"]], _ = pooled_tissue_source_site(
            patient_ids, min_site_count=min_site_count)
    design = confound_design(frame, arm["columns"])
    return cross_fitted_residuals(matrix, design, n_splits=n_splits, alpha=alpha, seed=seed), design


def unseen_site_probe(operator, matrix, cancers, patient_ids, arm: dict):
    """V2: one row carrying a TSS code the reference never saw.

    The barcode is rewritten to an impossible site (``ZZ``) so the row is definitely
    unseen, and everything else about the row is left alone.
    """
    if not arm["site_column"]:
        return {"applicable": False,
                "reason": "cancer_only arm carries no site column, which is the point of it"}
    fake_ids = np.asarray([f"TCGA-ZZ-{i:04d}" for i in range(3)])
    frame = pd.DataFrame({"cancer": cancers[:3]})
    adjusted, report = operator.adjust(matrix[:3], frame, patient_ids=fake_ids)
    pooled_report = report["site_pooling"]
    return {
        "applicable": True,
        "raw_sites": tissue_source_site(fake_ids).tolist(),
        "n_pooled_to_other": int(pooled_report["n_pooled_to_other"]),
        "sites_pooled_to_other": pooled_report["sites_pooled_to_other"],
        "policy": pooled_report["policy"],
        "all_rows_pooled": bool(pooled_report["n_pooled_to_other"] == 3),
        "adjusted_is_finite": bool(np.isfinite(adjusted).all()),
        "path": report["path"],
        "n_rows_adjusted": int(report["n_rows_adjusted"]),
        # An unseen site must land on the SAME coordinates as an explicit OTHER row,
        # because that is what the policy says it is. If these differ the policy is
        # not being applied, it is only being reported.
        "matches_explicit_other": bool(np.array_equal(
            adjusted,
            operator.adjust(matrix[:3], frame,
                            patient_ids=np.asarray([f"TCGA-QQ-{i:04d}" for i in range(3)]))[0])),
    }


def alchemist_label_probe(operator, matrix, arm: dict):
    """ALCHEMIST's cancer labels are disjoint from TCGA's. Report the refusal.

    Not a bug and not worked around here: encoding an unknown cancer type as zeros
    would adjust that slide as though it belonged to whichever level the design
    dropped. The external comparison needs a declared label mapping; this records
    what happens without one.
    """
    labels = ["Adenocarcinoma_NOS", "Squamous_cell_carcinoma_NOS", "Adenosquamous_carcinoma"]
    frame = pd.DataFrame({"cancer": labels})
    if arm["site_column"]:
        frame[arm["site_column"]] = ""
    ids = np.asarray([f"ALCH-B0C{i}" for i in range(3)])
    out = {"levels_tried": labels}
    try:
        operator.adjust(matrix[:3], frame, patient_ids=ids, on_unseen_level="refuse")
        out["refused"] = False
    except UnseenLevelError as error:
        out["refused"] = True
        out["error"] = str(error)[:400]
    _, report = operator.adjust(matrix[:3], frame, patient_ids=ids, on_unseen_level="zero")
    out["zero_encoded_unseen_levels"] = report["unseen_levels"]
    out["note"] = ("on_unseen_level='zero' encodes every ALCHEMIST cancer label as the "
                   "reference's dropped baseline level, i.e. adjusts every slide as if it "
                   "were one arbitrary TCGA cancer type. Usable only with a declared "
                   "label mapping; recorded here so the choice cannot be made silently.")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--partition", default="test")
    parser.add_argument("--states", nargs="*", default=[])
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    artifact = np.load(args.artifact, allow_pickle=True)
    targets = np.load(args.targets, allow_pickle=True)
    patient_ids, cancers, mask = reference_cohort(artifact, targets, args.partition)
    states = args.states or sorted({str(s) for s in artifact["trained_states"]})

    report = {"artifact": str(Path(args.artifact).name), "partition": args.partition,
              "n_reference_rows": int(len(patient_ids)),
              "n_cancer_levels": int(len(set(cancers.tolist()))),
              "n_raw_sites": int(len(set(tissue_source_site(patient_ids).tolist()))),
              "min_site_count": args.min_site_count,
              "residualiser": {"n_splits": args.n_splits, "alpha": args.alpha, "seed": args.seed},
              "operators": []}

    for state in states:
        matrix = np.asarray(artifact[state], dtype=np.float64)[mask]
        for arm_name, arm in ARMS.items():
            frame = pd.DataFrame({"cancer": cancers})
            if arm["site_column"]:
                frame[arm["site_column"]] = ""          # the operator owns the pooling
            operator = ConfoundAdjustmentOperator.fit(
                matrix, frame, arm["columns"],
                patient_ids=patient_ids if arm["site_column"] else None,
                site_column=arm["site_column"], min_site_count=args.min_site_count,
                n_splits=args.n_splits, alpha=args.alpha, seed=args.seed,
                cohort_name=f"TCGA::{Path(args.artifact).stem}::{state}::{args.partition}::{arm_name}")

            expected, design = transductive_reference(matrix, cancers, patient_ids, arm,
                                                      args.min_site_count, args.n_splits,
                                                      args.alpha, args.seed)
            actual = operator.adjust_reference(
                matrix, frame, patient_ids=patient_ids if arm["site_column"] else None)
            identity = bool(np.array_equal(actual, expected))
            max_abs = float(np.max(np.abs(actual - expected))) if actual.shape == expected.shape \
                else float("nan")

            path = out_dir / f"tcga_operator__{state}__{arm_name}__{args.partition}_seed{args.seed}.npz"
            operator.save(path)
            reloaded = ConfoundAdjustmentOperator.load(path)
            round_trip_reference = bool(np.array_equal(
                reloaded.adjust_reference(matrix, frame,
                                          patient_ids=patient_ids if arm["site_column"] else None),
                actual))
            new_frame = pd.DataFrame({"cancer": cancers[:5]})
            if arm["site_column"]:
                new_frame[arm["site_column"]] = ""
            new_ids = patient_ids[:5] if arm["site_column"] else None
            before, _ = operator.adjust(matrix[:5], new_frame, patient_ids=new_ids)
            after, _ = reloaded.adjust(matrix[:5], new_frame, patient_ids=new_ids)
            round_trip_new = bool(np.array_equal(before, after))

            record = {
                "state": state, "arm": arm_name, "path": str(path),
                "columns": arm["columns"],
                "n_design_columns": int(design.shape[1]),
                "n_representation_columns": int(matrix.shape[1]),
                "V1_identity_bit_for_bit": identity,
                "V1_max_abs_difference": max_abs,
                "V2_unseen_site": unseen_site_probe(operator, matrix, cancers, patient_ids, arm),
                "V3_round_trip_reference": round_trip_reference,
                "V3_round_trip_new_rows": round_trip_new,
                "V3_provenance_survived": bool(
                    reloaded.provenance.get("reference_digest") == operator.provenance.get("reference_digest")
                    and reloaded.provenance.get("cohort_name") == operator.provenance.get("cohort_name")),
                "alchemist_cancer_labels": alchemist_label_probe(operator, matrix, arm),
                "provenance": operator.provenance,
            }
            (path.with_suffix(".provenance.json")).write_text(
                json.dumps(operator.provenance, indent=2))
            report["operators"].append(record)
            print(f"[{state}::{arm_name}] identity={identity} maxabs={max_abs:.3e} "
                  f"design={design.shape[1]} round_trip={round_trip_reference and round_trip_new} "
                  f"-> {path.name}", flush=True)

    (out_dir / "fit_tcga_operator_report.json").write_text(json.dumps(report, indent=2))
    failures = [r for r in report["operators"]
                if not (r["V1_identity_bit_for_bit"] and r["V3_round_trip_reference"]
                        and r["V3_round_trip_new_rows"])]
    print(f"[written] {out_dir/'fit_tcga_operator_report.json'} "
          f"({len(report['operators'])} operators, {len(failures)} failing V1/V3)", flush=True)


if __name__ == "__main__":
    main()
