"""Rebuild an explicit paired split over the *maximal* available paired cohort.

Why this was rewritten.  The previous version built the paired split by
*intersecting* a source manifest with the loaded cohort.  That is safe only
while the cohort is shrinking.  It is not: the H-Optimus patch store finished
extracting after the split was frozen, so the loader now returns 6,443 paired
patients against 6,192 assigned ones, and `preflight.validate_runtime_split`
fail-closes with ``loaded paired cohort contains unassigned patients``.  That
check is right -- training on unlabelled patients is how a held-out cancer leaks
-- so the split has to be regenerated, and regenerating it should *keep* the
patients rather than throw them away.

What this module guarantees, in order of importance:

1.  **The held-out-cancer protocol is inherited, never re-drawn.**  Each cancer's
    role (development vs held-out) is read off the source split, and a patient
    is assigned by *its cancer's* role.  No patient can cross the boundary that
    defines the protocol.  Re-drawing the split from scratch would silently
    redefine what "held-out" means and make every prior number incomparable.
2.  **Every existing assignment is preserved exactly** -- train stays train, val
    stays val, test stays test.  The output is a strict superset of the input.
    Asserted, not assumed.
3.  **The cohort is maximised.**  Every newly available patient is assigned:
    held-out-cancer patients to ``test``, development-cancer patients to
    ``train``/``val`` by a deterministic seed-free hash order that tops each
    cancer's validation fold up towards the source split's validation fraction.
4.  **Eligibility is explicit.**  ``--require-patient-table`` intersects the
    cohort with a downstream artifact's patient index (in practice the prepared
    PanCan RNA table that PBS target construction needs).  Excluded patients are
    listed by name in the manifest instead of vanishing somewhere downstream
    where the count would never be reported.

The written split is then the authoritative cohort lookup: pass it together with
``--restrict-to-split`` so the runner trains on precisely these patients.
"""
from __future__ import annotations

import argparse
import collections
from hashlib import sha256
import json
import os
from pathlib import Path

import numpy as np

from morpheus.src.training.train_bio_query_former import load_bio_query_data
from .preflight import restrict_cohort_to_split, sha256_json, validate_runtime_split


def _order_key(patient: str) -> str:
    """Stable, seed-free ordering; identical on every machine and every rerun."""
    return sha256(str(patient).encode("utf-8")).hexdigest()


def _patient_index(path: str | Path, patient_column: str = "patient_id") -> set[str]:
    """Patient identifiers an eligibility artifact actually carries."""
    source = Path(path)
    if source.suffix == ".npz":
        with np.load(source, allow_pickle=False) as raw:
            if "patient_ids" not in raw.files:
                raise ValueError(f"{source} has no patient_ids array")
            return set(np.asarray(raw["patient_ids"]).astype(str).tolist())
    import pandas as pd

    frame = pd.read_parquet(source)
    if patient_column in frame.columns:
        return set(frame[patient_column].astype(str))
    if frame.index.name == patient_column:
        return set(frame.index.astype(str))
    raise ValueError(f"{source} exposes no '{patient_column}' column or index")


def _cancer_roles(sets: dict[str, set[str]], cancer_of: dict[str, str]) -> dict[str, str]:
    """Map each cancer to `development` or `heldout` exactly as the SOURCE split defined it."""
    roles: dict[str, str] = {}
    conflicts: set[str] = set()
    for partition, members in sets.items():
        role = "heldout" if partition == "test" else "development"
        for patient in members:
            cancer = cancer_of.get(patient)
            if cancer is None:
                continue
            if roles.setdefault(cancer, role) != role:
                conflicts.add(cancer)
    if conflicts:
        raise ValueError(f"source split assigns these cancers to both roles: {sorted(conflicts)}")
    return roles


def build_split(*, patient_ids, cancers, source_split: str | Path,
                eligible: set[str] | None = None) -> dict:
    ids = [str(patient) for patient in patient_ids]
    if len(set(ids)) != len(ids):
        raise ValueError("cohort contains duplicate patient identifiers")
    cancer_of = {patient: str(cancer) for patient, cancer in zip(ids, cancers)}

    payload = json.loads(Path(source_split).read_text(encoding="utf-8"))
    declared = payload.get("patient_ids")
    if not isinstance(declared, dict) or not {"train", "val", "test"}.issubset(declared):
        raise ValueError("source split needs patient_ids for train, val, and test")
    source = {name: {str(p) for p in declared[name]} for name in ("train", "val", "test")}
    if source["train"] & source["val"] or source["train"] & source["test"] or source["val"] & source["test"]:
        raise ValueError("source split assigns a patient to more than one partition")

    roles = _cancer_roles(source, cancer_of)
    cohort = set(ids)
    keep = cohort if eligible is None else cohort & set(eligible)
    source_universe = set().union(*source.values())

    unknown = sorted({patient for patient in keep if cancer_of[patient] not in roles})
    if unknown:
        raise ValueError(
            f"{len(unknown)} eligible patients belong to cancers the source split never gave a role, "
            f"so the held-out protocol cannot be inherited; examples={unknown[:5]}")

    assignment: dict[str, str] = {}
    new_development: dict[str, list[str]] = collections.defaultdict(list)
    for patient in sorted(keep):
        if patient in source["train"]:
            assignment[patient] = "train"
        elif patient in source["val"]:
            assignment[patient] = "val"
        elif patient in source["test"]:
            assignment[patient] = "test"
        elif roles[cancer_of[patient]] == "heldout":
            assignment[patient] = "test"
        else:
            new_development[cancer_of[patient]].append(patient)

    development_source = len(source["train"] | source["val"])
    val_fraction = len(source["val"]) / development_source if development_source else 0.0
    # Top each cancer's validation fold up towards the source fraction using only
    # NEW patients.  An existing assignment is never moved to hit a target.
    for cancer, newcomers in new_development.items():
        # Count only patients that survive eligibility: a source patient dropped by
        # `eligible` is not in the output, so including it here would size the
        # validation fold against a cohort that does not exist.
        existing_train = sum(1 for p in source["train"] & keep if cancer_of[p] == cancer)
        existing_val = sum(1 for p in source["val"] & keep if cancer_of[p] == cancer)
        target_val = int(round(val_fraction * (existing_train + existing_val + len(newcomers))))
        wanted = max(0, min(len(newcomers), target_val - existing_val))
        for index, patient in enumerate(sorted(newcomers, key=_order_key)):
            assignment[patient] = "val" if index < wanted else "train"

    sets = {name: sorted(p for p, part in assignment.items() if part == name)
            for name in ("train", "val", "test")}

    # --- invariants: assert, never assume ------------------------------------
    if set(assignment) != keep:
        raise AssertionError("every eligible patient must receive exactly one partition")
    for name in ("train", "val", "test"):
        moved = sorted(p for p in source[name] & keep if assignment[p] != name)
        if moved:
            raise AssertionError(f"source {name} assignments must be preserved; moved={moved[:5]}")
    if any(not members for members in sets.values()):
        raise AssertionError("no partition may be empty")
    train_cancers = {cancer_of[p] for p in sets["train"]}
    val_cancers = {cancer_of[p] for p in sets["val"]}
    test_cancers = {cancer_of[p] for p in sets["test"]}
    if test_cancers & train_cancers or test_cancers & val_cancers:
        raise AssertionError("held-out cancers leaked into development")
    if not val_cancers.issubset(train_cancers):
        raise AssertionError("validation cancers must be a subset of training cancers")

    added = sorted(keep - source_universe)
    return {
        "patient_ids": sets,
        "split_name": "tumor_state_heldout_cancer_maximal",
        "paired_from_source_split": str(Path(source_split)),
        "paired_source_digest": sha256_json(payload.get("patient_ids", {})),
        "derivation": "inherit each cancer's development/held-out role from the source split, preserve "
                      "every source assignment, then assign all newly available eligible patients",
        "actual_train_test_cancer_counts": [len(train_cancers), len(test_cancers)],
        "paired_train_cancers": sorted(train_cancers),
        "paired_test_cancers": sorted(test_cancers),
        "paired_cohort_count": len(keep),
        "source_cohort_count": len(source_universe),
        "loaded_cohort_count": len(cohort),
        "counts": {name: len(members) for name, members in sets.items()},
        "added_patient_count": len(added),
        "added_by_partition": dict(collections.Counter(assignment[p] for p in added)),
        "added_by_cancer": dict(sorted(collections.Counter(cancer_of[p] for p in added).items(),
                                       key=lambda item: -item[1])),
        "eligibility_excluded_count": len(cohort - keep),
        "eligibility_excluded_patients": sorted(cohort - keep),
        "source_patients_excluded_by_eligibility": sorted((source_universe & cohort) - keep),
        "source_patients_missing_from_cohort": sorted(source_universe - cohort),
        "source_val_fraction_of_development": val_fraction,
        "assignment_order": "sha256(patient_id) ascending; seed-free and machine-independent",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Materialize the maximal paired cohort for strict V2",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-config", required=True)
    parser.add_argument("--source-split", required=True,
                        help="frozen split whose cancer roles and assignments are inherited")
    parser.add_argument("--output", required=True)
    parser.add_argument("--require-patient-table", nargs="*", default=[],
                        help="one or more parquet/npz artifacts whose patient index bounds eligibility. Every "
                             "downstream artifact the run needs should be listed: a patient present in the "
                             "split but absent from the frozen RNA targets or the PBS codes aborts the run "
                             "hours later, and the exclusion is far cheaper to make explicit here.")
    parser.add_argument("--patient-column", default="patient_id")
    parser.add_argument("--expected-development-cancers", type=int, default=11)
    parser.add_argument("--expected-heldout-cancers", type=int, default=21)
    args = parser.parse_args()

    data = load_bio_query_data(args.data_config, args.source_split, wsi_mode="hoptimus_patch")
    eligible: set[str] | None = None
    per_table: dict[str, int] = {}
    for table in args.require_patient_table:
        available = _patient_index(table, args.patient_column)
        per_table[str(table)] = len(available)
        eligible = available if eligible is None else (eligible & available)

    payload = build_split(patient_ids=data.patient_ids, cancers=data.cancers,
                          source_split=args.source_split, eligible=eligible)
    payload["require_patient_table"] = list(map(str, args.require_patient_table))
    payload["require_patient_table_sizes"] = per_table
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Validate a candidate and only then publish it. A split that fails the protocol
    # check must not be left at --output where the next command silently picks it up.
    candidate = destination.with_suffix(destination.suffix + ".candidate")
    candidate.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    try:
        # End-to-end proof that the artifact is one the runner will actually
        # accept, obtained through the runner's own code path.
        restricted, restriction = restrict_cohort_to_split(
            load_bio_query_data(args.data_config, str(candidate), wsi_mode="hoptimus_patch"), candidate)
        manifest = validate_runtime_split(candidate, restricted.patient_ids, restricted.cancers,
                                          restricted.split, args.expected_development_cancers,
                                          args.expected_heldout_cancers)
    except BaseException:
        candidate.unlink(missing_ok=True)
        raise
    os.replace(candidate, destination)
    manifest["split_path"] = str(destination)
    print(json.dumps({"written": str(destination),
                      "summary": {key: value for key, value in payload.items()
                                  if key not in {"patient_ids", "eligibility_excluded_patients"}},
                      "restriction": {k: v for k, v in restriction.items() if k != "dropped_patients"},
                      "runtime_split_manifest": manifest}, indent=2))


if __name__ == "__main__":
    main()
