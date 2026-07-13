"""TCGA identifier parsing helpers."""

from __future__ import annotations

import re
from dataclasses import dataclass

TCGA_PATIENT_RE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", re.IGNORECASE)


@dataclass(frozen=True)
class TcgaBarcode:
    """Parsed TCGA barcode components used for leakage-safe joins."""

    raw: str
    patient_id: str | None
    sample_id: str | None
    slide_id: str | None


def normalize_patient_id(value: object) -> str | None:
    """Return a canonical TCGA patient ID when one can be parsed."""
    if value is None:
        return None
    text = str(value).strip()
    match = TCGA_PATIENT_RE.search(text)
    if match:
        return match.group(1).upper()
    if text.upper().startswith("TCGA-") and len(text) >= 12:
        return text[:12].upper()
    return text.upper() if text else None


def parse_tcga_barcode(value: object) -> TcgaBarcode:
    text = "" if value is None else str(value).strip()
    patient_id = normalize_patient_id(text)
    sample_id = text[:15].upper() if text.upper().startswith("TCGA-") and len(text) >= 15 else None
    slide_id = text.upper() if text.upper().startswith("TCGA-") else None
    return TcgaBarcode(raw=text, patient_id=patient_id, sample_id=sample_id, slide_id=slide_id)
