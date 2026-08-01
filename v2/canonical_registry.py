"""Canonical TCGA participant registry and explicit CDR survival endpoints."""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd


class RegistryError(ValueError):
    pass


def canonical_tcga_patient_id(value: object) -> str:
    """Return the 12-character TCGA participant barcode or fail closed."""
    text = str(value).strip().upper()
    pieces = text.split("-")
    if len(pieces) < 3 or pieces[0] != "TCGA" or any(not part for part in pieces[:3]):
        raise RegistryError(f"not a TCGA participant/sample barcode: {value!r}")
    return "-".join(pieces[:3])


def _id_column(frame: pd.DataFrame) -> str:
    candidates = ("patient_id", "bcr_patient_barcode", "Patient ID", "patient", "submitter_id")
    found = [column for column in candidates if column in frame.columns]
    if not found:
        raise RegistryError(f"expected exactly one patient identifier column, found {found}")
    if len(found) == 1:
        return found[0]
    # TCGA CDR/PanCan tables legitimately carry both a canonical patient_id
    # and the original bcr_patient_barcode.  Accept that redundant schema only
    # when the overlapping non-null rows resolve to the same participant; a
    # disagreement is an alignment error, never a reason to choose silently.
    preferred = "bcr_patient_barcode" if "bcr_patient_barcode" in found else "patient_id"
    reference = frame[preferred]
    for other in found:
        if other == preferred:
            continue
        overlap = reference.notna() & frame[other].notna()
        if overlap.any():
            def equivalent(left: object, right: object) -> bool:
                canonical = canonical_tcga_patient_id(left)
                text = str(right).strip().upper()
                # PanCan CDR's auxiliary `patient_id` is frequently only the
                # four-character participant suffix (e.g. A5J1), while the
                # BCR barcode is the canonical identifier.  It is redundant,
                # not contradictory, precisely when that suffix agrees.
                try:
                    return canonical_tcga_patient_id(text) == canonical
                except RegistryError:
                    return text == canonical.rsplit("-", 1)[-1]
            if not all(equivalent(left, right) for left, right in zip(reference.loc[overlap], frame.loc[overlap, other])):
                raise RegistryError(f"conflicting participant identifier columns: {preferred!r} and {other!r}")
    return preferred


def build_canonical_registry(sources: Mapping[str, pd.DataFrame], *, cancer_column: str = "cancer") -> pd.DataFrame:
    """Merge modality source tables into one participant-level availability registry.

    Source rows are retained as pipe-separated provenance IDs.  Cancer labels
    may be absent, but contradictory non-null labels for one participant are a
    hard error rather than an arbitrary first-row choice.
    """
    records: list[dict[str, object]] = []
    for modality, table in sources.items():
        if not modality or not isinstance(table, pd.DataFrame):
            raise RegistryError("sources must map non-empty modality names to DataFrames")
        frame = table.copy()
        column = _id_column(frame)
        for _, row in frame.iterrows():
            raw = row[column]
            if pd.isna(raw):
                continue
            patient = canonical_tcga_patient_id(raw)
            cancer = row.get(cancer_column, row.get("cancer_type", np.nan))
            records.append({"patient_id": patient, "modality": modality, "source_id": str(raw),
                            "cancer": None if pd.isna(cancer) else str(cancer)})
    if not records:
        raise RegistryError("registry sources produced no TCGA patient records")
    detail = pd.DataFrame(records)
    contradictions = detail.dropna(subset=["cancer"]).groupby("patient_id")["cancer"].nunique()
    bad = contradictions[contradictions > 1]
    if len(bad):
        raise RegistryError(f"contradictory cancer labels for patients: {bad.index[:5].tolist()}")
    patients = sorted(detail.patient_id.unique())
    registry = pd.DataFrame({"patient_id": patients})
    cancer = detail.dropna(subset=["cancer"]).drop_duplicates("patient_id").set_index("patient_id")["cancer"]
    registry["cancer"] = registry.patient_id.map(cancer)
    for modality in sorted(sources):
        grouped = detail.loc[detail.modality.eq(modality)].groupby("patient_id")["source_id"].agg(lambda values: "|".join(sorted(set(values))))
        registry[f"has_{modality}"] = registry.patient_id.isin(grouped.index)
        registry[f"source_ids_{modality}"] = registry.patient_id.map(grouped)
    registry["registry_digest"] = sha256(registry.to_csv(index=False, na_rep="<NA>").encode("utf-8")).hexdigest()
    return registry


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame[column], errors="coerce") if column in frame else pd.Series(np.nan, index=frame.index)


def build_tcga_cdr_outcomes(clinical: pd.DataFrame) -> pd.DataFrame:
    """Build participant-level OS and PFI labels from TCGA clinical follow-up.

    OS is death time for deceased participants, otherwise last follow-up.
    PFI uses documented progression/new tumour event time when present; rows
    without a reliable PFI time remain unavailable rather than being guessed.
    """
    frame = clinical.copy()
    column = _id_column(frame)
    frame = frame.loc[frame[column].notna()].copy()
    frame["patient_id"] = frame[column].map(canonical_tcga_patient_id)
    vital = frame.get("vital_status", pd.Series("", index=frame.index)).astype(str).str.strip().str.upper()
    death = _numeric(frame, "days_to_death")
    followup = _numeric(frame, "days_to_last_followup")
    # PanCan follow-up uses both "Dead" and "Deceased" across releases.
    os_event = vital.isin({"DECEASED", "DEAD"}) & death.notna() & (death >= 0)
    os_time = np.where(os_event, death, followup)
    progression_time = _numeric(frame, "days_to_patient_progression_free")
    if progression_time.isna().all():
        progression_time = _numeric(frame, "days_to_tumor_progression")
    new_event_time = _numeric(frame, "days_to_new_tumor_event_after_initial_treatment")
    progression_time = progression_time.where(progression_time.notna(), new_event_time)
    progression = frame.get("patient_progression_status", pd.Series("", index=frame.index)).astype(str).str.upper()
    event_text = frame.get("new_tumor_event_after_initial_treatment", pd.Series("", index=frame.index)).astype(str).str.upper()
    pfi_event = (progression.str.contains("PROGRESS|WITH TUMOR", regex=True, na=False) |
                 event_text.isin({"YES", "TRUE", "1"})) & progression_time.notna() & (progression_time >= 0)
    pfi_time = np.where(pfi_event, progression_time, followup)
    out = pd.DataFrame({"patient_id": frame.patient_id, "os_time_days": os_time,
                        "os_event": os_event.astype(float), "pfi_time_days": pfi_time,
                        "pfi_event": pfi_event.astype(float), "source_patient_id": frame[column].astype(str)})
    # Multiple records can arise from follow-up tables.  Death/progression wins;
    # otherwise retain the longest observed follow-up, all without test labels.
    rows = []
    for patient, group in out.groupby("patient_id", sort=True):
        def pick(prefix: str) -> tuple[float, float]:
            time, event = group[f"{prefix}_time_days"].to_numpy(float), group[f"{prefix}_event"].to_numpy(float)
            event_rows = np.isfinite(time) & event.astype(bool)
            eligible = event_rows if event_rows.any() else np.isfinite(time)
            if not eligible.any():
                return np.nan, np.nan
            index = np.where(eligible)[0][np.argmax(time[eligible])]
            return float(time[index]), float(event[index])
        os_time_value, os_event_value = pick("os")
        pfi_time_value, pfi_event_value = pick("pfi")
        rows.append({"patient_id": patient, "os_time_days": os_time_value, "os_event": os_event_value,
                     "pfi_time_days": pfi_time_value, "pfi_event": pfi_event_value,
                     "source_patient_ids": "|".join(sorted(group.source_patient_id.unique()))})
    return pd.DataFrame(rows)


def survival_coverage(outcomes: pd.DataFrame, canonical_patient_ids: list[str], split: list[str]) -> dict[str, object]:
    """Report, never hide, split-specific survival availability and event counts."""
    if len(canonical_patient_ids) != len(split):
        raise RegistryError("canonical_patient_ids and split must have equal length")
    indexed = outcomes.set_index("patient_id")
    report: dict[str, object] = {}
    for endpoint in ("os", "pfi"):
        time = indexed.reindex(canonical_patient_ids)[f"{endpoint}_time_days"].to_numpy(float)
        event = indexed.reindex(canonical_patient_ids)[f"{endpoint}_event"].to_numpy(float)
        valid = np.isfinite(time) & np.isfinite(event) & (time > 0)
        report[endpoint] = {part: {"n": int(((np.asarray(split) == part) & valid).sum()),
                                  "events": int(((np.asarray(split) == part) & valid & event.astype(bool)).sum())}
                            for part in ("train", "val", "test")}
    report["outcome_digest"] = sha256(outcomes.to_csv(index=False, na_rep="<NA>").encode("utf-8")).hexdigest()
    return report


def audit_survival_alignment(
    outcomes: pd.DataFrame, canonical_patient_ids: list[str], split: list[str]
) -> dict[str, object]:
    """Fail-closed audit for survival labels before any Cox fit.

    The coverage summary alone is not enough: a malformed CDR table can have
    apparently adequate event counts while silently failing to align patients.
    This report makes every unmatched or duplicate identifier visible to the
    controller and is serialisable as a run manifest.
    """
    required = {"patient_id", "os_time_days", "os_event", "pfi_time_days", "pfi_event"}
    missing = required - set(outcomes.columns)
    if missing:
        raise RegistryError(f"survival outcome table is missing {sorted(missing)}")
    if len(canonical_patient_ids) != len(split):
        raise RegistryError("canonical_patient_ids and split must have equal length")
    ids = [canonical_tcga_patient_id(value) for value in canonical_patient_ids]
    if len(ids) != len(set(ids)):
        raise RegistryError("canonical cohort contains duplicate patient IDs")
    frame = outcomes.copy()
    frame["patient_id"] = frame["patient_id"].map(canonical_tcga_patient_id)
    duplicate_ids = sorted(frame.loc[frame.patient_id.duplicated(keep=False), "patient_id"].unique())
    if duplicate_ids:
        raise RegistryError(f"survival outcome table has duplicate canonical IDs: {duplicate_ids[:5]}")
    known, cohort = set(frame.patient_id), set(ids)
    matched = sorted(known & cohort)
    by_split = {}
    split_values = np.asarray(split, dtype=str)
    for name in ("train", "val", "test"):
        members = {patient for patient, part in zip(ids, split_values) if part == name}
        by_split[name] = {
            "cohort_patients": len(members),
            "matched_outcomes": len(members & known),
            "unmatched_patients": sorted(members - known),
        }
    report = {
        "schema_version": 1,
        "canonical_cohort_count": len(ids),
        "outcome_row_count": len(frame),
        "matched_outcome_count": len(matched),
        "cohort_patients_without_outcome": sorted(cohort - known),
        "outcome_patients_outside_cohort": sorted(known - cohort),
        "split_alignment": by_split,
        "coverage": survival_coverage(frame, ids, list(split)),
    }
    report["audit_digest"] = sha256(json.dumps(report, sort_keys=True).encode("utf-8")).hexdigest()
    return report
