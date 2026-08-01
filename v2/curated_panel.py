"""Frozen cancer-relevant Reactome/KEGG molecular evaluation panel.

This module intentionally does *not* search a GMT alphabetically.  The panel
below is a reviewed list of biological processes that are plausible pathology
targets.  A concrete GMT release is resolved once, written to a manifest, and
must be reused verbatim for every compared method.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from typing import Mapping, Sequence


class PanelResolutionError(ValueError):
    """Raised when a frozen target cannot be resolved unambiguously."""


@dataclass(frozen=True)
class PanelEntry:
    name: str
    family: str
    candidates: tuple[str, ...]


# Names are deliberately explicit and reviewable.  Candidate aliases cover
# common MSigDB Reactome/KEGG naming differences; resolution records the exact
# source name so no later release can silently change the endpoint.
CURATED_CANCER_PANEL: tuple[PanelEntry, ...] = (
    PanelEntry("antigen_processing_presentation", "immune", ("KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",)),
    PanelEntry("interferon_signaling", "immune", ("REACTOME_INTERFERON_SIGNALING",)),
    PanelEntry("interferon_gamma_signaling", "immune", ("REACTOME_INTERFERON_GAMMA_SIGNALING",)),
    PanelEntry("adaptive_immune_system", "immune", ("REACTOME_ADAPTIVE_IMMUNE_SYSTEM",)),
    PanelEntry("tcr_signaling", "immune", ("REACTOME_TCR_SIGNALING",)),
    PanelEntry("chemokine_signaling", "immune", ("KEGG_CHEMOKINE_SIGNALING_PATHWAY",)),
    PanelEntry("complement_coagulation", "immune", ("KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",)),
    PanelEntry("toll_like_receptor", "immune", ("KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",)),
    PanelEntry("extracellular_matrix_organization", "stroma_ecm", ("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",)),
    PanelEntry("collagen_formation", "stroma_ecm", ("REACTOME_COLLAGEN_FORMATION",)),
    PanelEntry("integrin_cell_surface_interactions", "stroma_ecm", ("REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",)),
    PanelEntry("focal_adhesion", "stroma_ecm", ("KEGG_FOCAL_ADHESION",)),
    PanelEntry("tgfb_signaling", "stroma_ecm", ("KEGG_TGF_BETA_SIGNALING_PATHWAY", "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX")),
    PanelEntry("hypoxia_response", "metabolism_stress", ("REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",)),
    PanelEntry("glycolysis_gluconeogenesis", "metabolism_stress", ("KEGG_GLYCOLYSIS_GLUCONEOGENESIS",)),
    PanelEntry("oxidative_phosphorylation", "metabolism_stress", ("KEGG_OXIDATIVE_PHOSPHORYLATION",)),
    PanelEntry("hif1_signaling", "metabolism_stress", ("KEGG_HIF_1_SIGNALING_PATHWAY",)),
    PanelEntry("angiogenesis", "metabolism_stress", ("REACTOME_SIGNALING_BY_VEGF", "KEGG_VEGF_SIGNALING_PATHWAY")),
    PanelEntry("cell_cycle", "proliferation", ("KEGG_CELL_CYCLE",)),
    PanelEntry("mitotic_cell_cycle", "proliferation", ("REACTOME_CELL_CYCLE_MITOTIC",)),
    PanelEntry("dna_replication", "proliferation", ("REACTOME_DNA_REPLICATION", "KEGG_DNA_REPLICATION")),
    PanelEntry("p53_signaling", "damage_repair", ("KEGG_P53_SIGNALING_PATHWAY",)),
    PanelEntry("dna_repair", "damage_repair", ("REACTOME_DNA_REPAIR",)),
    PanelEntry("homologous_recombination", "damage_repair", ("KEGG_HOMOLOGOUS_RECOMBINATION",)),
    PanelEntry("mismatch_repair", "damage_repair", ("KEGG_MISMATCH_REPAIR",)),
    PanelEntry("base_excision_repair", "damage_repair", ("KEGG_BASE_EXCISION_REPAIR",)),
    PanelEntry("nucleotide_excision_repair", "damage_repair", ("KEGG_NUCLEOTIDE_EXCISION_REPAIR",)),
    PanelEntry("apoptosis", "damage_repair", ("KEGG_APOPTOSIS", "REACTOME_APOPTOTIC_EXECUTION_PHASE")),
    PanelEntry("mapk_signaling", "oncogenic_signaling", ("KEGG_MAPK_SIGNALING_PATHWAY",)),
    PanelEntry("pi3k_akt_signaling", "oncogenic_signaling", ("KEGG_PI3K_AKT_SIGNALING_PATHWAY",)),
    PanelEntry("wnt_signaling", "oncogenic_signaling", ("KEGG_WNT_SIGNALING_PATHWAY",)),
    PanelEntry("notch_signaling", "oncogenic_signaling", ("KEGG_NOTCH_SIGNALING_PATHWAY",)),
)


# These names are the contract between molecular scoring and the constrained
# card evaluator. They are a fixed research-hypothesis vocabulary, not model
# inputs, and are scored from hidden RNA only in WSI-only evaluations.
FROZEN_MECHANISM_PROGRAMMES: dict[str, tuple[str, ...]] = {
    "cytolytic_activity": ("NKG7", "PRF1", "GZMA", "GZMB", "GNLY", "CTSW", "CD8A", "CD8B"),
    "interferon_gamma": ("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "IDO1", "GBP1", "HLA-DRA"),
    "antigen_presentation": ("B2M", "TAP1", "TAP2", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "PSMB8", "PSMB9"),
    "myeloid_macrophage": ("LYZ", "LST1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CSF1R", "CD68"),
    "stroma_caf": ("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "FAP", "PDGFRB"),
    "tgf_beta_emt": ("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "VIM", "SNAI1", "ZEB1", "FN1"),
    "immune_exclusion": ("TGFB1", "VEGFA", "CXCL12", "FAP", "COL1A1", "COL3A1", "PDGFRB", "IL6", "CCL2"),
    "hypoxia": ("HIF1A", "CA9", "EGLN3", "BNIP3", "SLC2A1", "LDHA", "VEGFA", "ADM", "ENO1"),
    "glycolysis": ("HK2", "GPI", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "SLC2A1"),
    "angiogenesis": ("VEGFA", "KDR", "FLT1", "ANGPT2", "VWF", "ESM1", "PECAM1", "EMCN", "RAMP2"),
    "proliferation": ("MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "CDK1", "CCNB1", "AURKA"),
    "dna_repair": ("BRCA1", "BRCA2", "RAD51", "ATM", "ATR", "CHEK1", "CHEK2", "PARP1", "XRCC1", "FANCD2"),
    "apoptosis_senescence": ("BAX", "BBC3", "PMAIP1", "CASP3", "CASP8", "CDKN1A", "CDKN2A", "GADD45A", "SERPINE1"),
    "emt": ("VIM", "FN1", "CDH2", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "ITGA5"),
    "mechanotransduction": ("YAP1", "WWTR1", "TEAD1", "TEAD4", "ITGB1", "PXN", "PTK2", "TLN1", "VCL", "ACTA2"),
}


def frozen_mechanism_programme_manifest() -> dict[str, object]:
    payload = {name: list(genes) for name, genes in FROZEN_MECHANISM_PROGRAMMES.items()}
    return {"format_version": 1, "programmes": payload,
            "digest": sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()}


def freeze_curated_panel(
    signatures: Mapping[str, Sequence[str]],
    measured_genes: Sequence[str],
    *,
    training_signatures: Mapping[str, Sequence[str]] | None = None,
    min_mapped_genes: int = 10,
    max_training_jaccard: float = 0.35,
) -> dict[str, object]:
    """Resolve the reviewed panel against one immutable GMT release.

    Candidate order is a reviewable, fixed preference order (Reactome where
    listed first). Missing signatures, too-small mappings, and excessive direct
    overlap with V2 programme-training sets are recorded as unavailable, never
    substituted by a similarly named pathway.
    """
    genes = {str(g) for g in measured_genes}
    training = {name: {str(g) for g in values} for name, values in (training_signatures or {}).items()}
    targets: list[dict[str, object]] = []
    for entry in CURATED_CANCER_PANEL:
        available = [name for name in entry.candidates if name in signatures]
        source = available[0] if available else None
        mapped = [] if source is None else sorted({str(g) for g in signatures[source] if str(g) in genes})
        overlap = 0.0
        overlap_name = None
        if mapped and training:
            member = set(mapped)
            overlap_name, overlap = max(
                ((name, len(member & other) / len(member | other)) for name, other in training.items() if member | other),
                key=lambda item: item[1], default=(None, 0.0),
            )
        status = "eligible"
        if source is None:
            status = "unavailable_signature_missing"
        elif len(mapped) < min_mapped_genes:
            status = "unavailable_too_few_mapped_genes"
        elif overlap > max_training_jaccard:
            status = "transfer_diagnostic_overlap_with_training"
        targets.append({"name": entry.name, "family": entry.family, "source_signature": source,
                        "available_approved_candidates": available, "source_selection_rule": "fixed_priority_order",
                        "mapped_genes": mapped, "n_mapped_genes": len(mapped), "max_training_jaccard": overlap,
                        "max_overlap_training_signature": overlap_name, "status": status})
    digest = sha256(json.dumps(targets, sort_keys=True).encode("utf-8")).hexdigest()
    return {"format_version": 1, "panel": [asdict(x) for x in CURATED_CANCER_PANEL], "targets": targets,
            "eligible_targets": [x["name"] for x in targets if x["status"] == "eligible"], "digest": digest,
            "min_mapped_genes": min_mapped_genes, "max_training_jaccard": max_training_jaccard}
