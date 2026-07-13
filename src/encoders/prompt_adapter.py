"""Adapter for molecular prompting outputs or pathway score vectors."""

from __future__ import annotations

from morpheus.src.encoders.clinical_adapter import ClinicalTokenAdapter


class PromptTokenAdapter(ClinicalTokenAdapter):
    """Pathway/prompt score adapter."""

    pass
