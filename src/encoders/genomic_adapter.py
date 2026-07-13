"""Generic tabular genomic adapter for SNV/CNV feature matrices."""

from __future__ import annotations

from morpheus.src.encoders.clinical_adapter import ClinicalTokenAdapter


class GenomicTokenAdapter(ClinicalTokenAdapter):
    """Shape-safe adapter for available SNV/CNV matrices."""

    pass

