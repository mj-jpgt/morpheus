"""Deterministic, evidence-grounded tumour hypothesis cards.

This module deliberately does *not* call an LLM.  It converts WSI-derived
programme scores into a small, auditable card whose vocabulary is closed.  A
separate system may render a validated card into prose, but this module is the
scientific record and is safe to use in automated evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, log2
from typing import Any, Mapping, Sequence

import numpy as np


CARD_MODE = "research_hypothesis_not_clinical_decision"
CARD_VERSION = "morpheus_hypothesis_card_v1"
MAX_EVIDENCE_IDS = 8

# This is intentionally a closed vocabulary.  New biological assertions must
# be added by editing and reviewing this registry, never by a text generator.
PROGRAMMES = frozenset(
    {
        "cytolytic_activity", "interferon_gamma", "antigen_presentation",
        "myeloid_macrophage", "stroma_caf", "tgf_beta_emt",
        "immune_exclusion", "hypoxia", "glycolysis", "angiogenesis",
        "proliferation", "dna_repair", "apoptosis_senescence", "emt",
        "mechanotransduction",
    }
)


@dataclass(frozen=True)
class MechanismTemplate:
    name: str
    required: tuple[tuple[str, str], ...]
    assay: str
    predicted_observation: str


MECHANISMS: tuple[MechanismTemplate, ...] = (
    MechanismTemplate(
        "stromal_immune_exclusion",
        (("stroma_caf", "high"), ("tgf_beta_emt", "high"),
         ("cytolytic_activity", "low")),
        "spatial_transcriptomics_or_multiplex_if",
        "CAF/TGF-beta-rich regions spatially separate from cytolytic T-cell regions.",
    ),
    MechanismTemplate(
        "immune_active_interferon",
        (("cytolytic_activity", "high"), ("interferon_gamma", "high"),
         ("antigen_presentation", "high")),
        "multiplex_if_or_rna_panel",
        "Cytolytic T-cell and interferon-response markers are co-localized with tumour regions.",
    ),
    MechanismTemplate(
        "hypoxic_glycolytic_adaptation",
        (("hypoxia", "high"), ("glycolysis", "high")),
        "spatial_transcriptomics_or_hypoxia_ihc",
        "Hypoxia and glycolysis markers concentrate in poorly perfused tumour regions.",
    ),
    MechanismTemplate(
        "proliferative_repair_stress",
        (("proliferation", "high"), ("dna_repair", "high")),
        "rna_panel_or_ihc",
        "Cell-cycle and DNA-repair markers are jointly elevated in tumour cells.",
    ),
    MechanismTemplate(
        "emt_invasion_programme",
        (("emt", "high"), ("tgf_beta_emt", "high")),
        "multiplex_if_or_spatial_transcriptomics",
        "EMT and TGF-beta-associated cells enrich at invasive tumour-stroma boundaries.",
    ),
    MechanismTemplate(
        "angiogenic_hypoxic_niche",
        (("angiogenesis", "high"), ("hypoxia", "high")),
        "spatial_transcriptomics_or_vascular_ihc",
        "Angiogenic and hypoxia signals co-occur in vascularly heterogeneous tumour regions.",
    ),
    MechanismTemplate(
        "myeloid_suppressed_inflammation",
        (("myeloid_macrophage", "high"), ("cytolytic_activity", "low")),
        "multiplex_if_or_spatial_transcriptomics",
        "Myeloid/macrophage markers are enriched while cytolytic T-cell markers remain sparse.",
    ),
)
MECHANISM_BY_NAME = {item.name: item for item in MECHANISMS}
ASSAYS = frozenset(item.assay for item in MECHANISMS) | {"no_assay_recommended"}


def _as_score_rows(
    patient_ids: Sequence[str], programme_scores: Any, programme_names: Sequence[str]
) -> dict[str, dict[str, float]]:
    """Normalise mapping or dense score matrix inputs without changing values."""
    ids = [str(value) for value in patient_ids]
    if len(ids) != len(set(ids)):
        raise ValueError("patient_ids must be unique")
    names = [str(value) for value in programme_names]
    if len(names) != len(set(names)):
        raise ValueError("programme_names must be unique")
    if isinstance(programme_scores, Mapping):
        if set(ids) - set(map(str, programme_scores)):
            raise ValueError("programme score mapping is missing a patient")
        rows: dict[str, dict[str, float]] = {}
        for patient_id in ids:
            row = programme_scores[patient_id]
            if not isinstance(row, Mapping):
                raise TypeError("mapping scores must map patient IDs to mappings")
            rows[patient_id] = {name: float(row[name]) for name in names if name in row}
        return rows
    values = np.asarray(programme_scores, dtype=float)
    if values.shape != (len(ids), len(names)):
        raise ValueError("score matrix shape must be (n_patients, n_programmes)")
    return {
        patient_id: {name: float(value) for name, value in zip(names, values[index])}
        for index, patient_id in enumerate(ids)
    }


def _direction(score: float, threshold: float) -> str:
    if score >= threshold:
        return "high"
    if score <= -threshold:
        return "low"
    return "neutral"


def _confidence(score: float) -> float:
    """A deterministic monotonic score, not a calibrated probability claim."""
    return float(min(0.99, max(0.01, 1.0 / (1.0 + exp(-abs(score))))))


def _mechanism_for(directions: Mapping[str, str], scores: Mapping[str, float]) -> MechanismTemplate | None:
    candidates = [
        template for template in MECHANISMS
        if all(directions.get(programme) == direction for programme, direction in template.required)
    ]
    if not candidates:
        return None
    # Tie break by registry order to make cards reproducible.
    return max(candidates, key=lambda item: sum(abs(scores.get(p, 0.0)) for p, _ in item.required))


def generate_hypothesis_cards(
    patient_ids: Sequence[str],
    programme_scores: Any,
    programme_names: Sequence[str],
    *,
    evidence_ids_by_patient: Mapping[str, Sequence[str]] | None = None,
    threshold: float = 0.5,
    max_programmes: int = 3,
) -> list[dict[str, Any]]:
    """Create deterministic cards from WSI-predicted programme score matrices.

    Unknown programmes are ignored, preventing unreviewed label names from
    entering generated output.  ``evidence_ids_by_patient`` is copied only;
    evidence is never invented by the generator.
    """
    if threshold < 0:
        raise ValueError("threshold must be non-negative")
    if max_programmes < 1:
        raise ValueError("max_programmes must be positive")
    rows = _as_score_rows(patient_ids, programme_scores, programme_names)
    cards: list[dict[str, Any]] = []
    for patient_id in map(str, patient_ids):
        scores = {name: score for name, score in rows[patient_id].items() if name in PROGRAMMES and np.isfinite(score)}
        directions = {name: _direction(score, threshold) for name, score in scores.items()}
        ranked = sorted(scores, key=lambda name: (-abs(scores[name]), name))[:max_programmes]
        predicted = [
            {"programme": name, "direction": directions[name], "score": float(scores[name]),
             "confidence": _confidence(scores[name])}
            for name in ranked
        ]
        mechanism = _mechanism_for(directions, scores)
        evidence = [] if evidence_ids_by_patient is None else list(evidence_ids_by_patient.get(patient_id, []))[:MAX_EVIDENCE_IDS]
        if not all(isinstance(item, str) and item for item in evidence):
            raise ValueError("evidence IDs must be non-empty strings")
        if mechanism is None:
            composite = {"template": "indeterminate", "supporting_programmes": []}
            falsification = {"recommended_assay": "no_assay_recommended", "predicted_observation": "No composite mechanism met the predeclared evidence rule."}
        else:
            composite = {"template": mechanism.name, "supporting_programmes": [p for p, _ in mechanism.required]}
            falsification = {"recommended_assay": mechanism.assay, "predicted_observation": mechanism.predicted_observation}
        confidence = float(np.mean([item["confidence"] for item in predicted])) if predicted else 0.01
        card = {
            "version": CARD_VERSION,
            "mode": CARD_MODE,
            "patient_id": patient_id,
            "input": "wsi_only",
            "predicted_programmes": predicted,
            "composite_hypothesis": composite,
            "falsification": falsification,
            "uncertainty": {"confidence": confidence, "abstain": mechanism is None,
                            "reason": "No composite template met the predeclared rule." if mechanism is None else "Predeclared programme evidence is concordant."},
            "evidence": {"top_patch_or_slide_evidence_ids": evidence,
                         "programme_scores": {name: float(scores[name]) for name in ranked}},
        }
        validate_hypothesis_card(card)
        cards.append(card)
    return cards


def abstain_hypothesis_cards(patient_ids: Sequence[str], *, reason: str) -> list[dict[str, Any]]:
    """Emit schema-valid abstentions when the quantitative card gate fails."""
    cards: list[dict[str, Any]] = []
    for patient_id in map(str, patient_ids):
        card = {
            "version": CARD_VERSION,
            "mode": CARD_MODE,
            "patient_id": patient_id,
            "input": "wsi_only",
            "predicted_programmes": [],
            "composite_hypothesis": {"template": "indeterminate", "supporting_programmes": []},
            "falsification": {"recommended_assay": "no_assay_recommended", "predicted_observation": "No hypothesis emitted because the quantitative card gate did not pass."},
            "uncertainty": {"confidence": 0.01, "abstain": True, "reason": str(reason)},
            "evidence": {"top_patch_or_slide_evidence_ids": [], "programme_scores": {}},
        }
        validate_hypothesis_card(card)
        cards.append(card)
    return cards


def validate_hypothesis_card(card: Mapping[str, Any]) -> None:
    """Raise ``ValueError`` unless a card conforms to the closed safety schema."""
    required = {"version", "mode", "patient_id", "input", "predicted_programmes", "composite_hypothesis", "falsification", "uncertainty", "evidence"}
    if set(card) != required:
        raise ValueError("card has missing or unexpected fields")
    if card["version"] != CARD_VERSION or card["mode"] != CARD_MODE or card["input"] != "wsi_only":
        raise ValueError("invalid card metadata")
    if not isinstance(card["patient_id"], str) or not card["patient_id"]:
        raise ValueError("patient_id must be a non-empty string")
    programmes = card["predicted_programmes"]
    if not isinstance(programmes, list) or len(programmes) > 3:
        raise ValueError("invalid predicted_programmes")
    seen: set[str] = set()
    for item in programmes:
        if set(item) != {"programme", "direction", "score", "confidence"} or item["programme"] not in PROGRAMMES:
            raise ValueError("unsupported programme claim")
        if item["programme"] in seen or item["direction"] not in {"high", "low", "neutral"}:
            raise ValueError("invalid programme claim")
        if not np.isfinite(item["score"]) or not 0.0 <= item["confidence"] <= 1.0:
            raise ValueError("invalid programme score/confidence")
        seen.add(item["programme"])
    composite = card["composite_hypothesis"]
    if set(composite) != {"template", "supporting_programmes"}:
        raise ValueError("invalid composite hypothesis")
    template = composite["template"]
    if template != "indeterminate" and template not in MECHANISM_BY_NAME:
        raise ValueError("unsupported mechanism")
    if not isinstance(composite["supporting_programmes"], list) or not set(composite["supporting_programmes"]).issubset(PROGRAMMES):
        raise ValueError("invalid supporting programmes")
    falsification = card["falsification"]
    if set(falsification) != {"recommended_assay", "predicted_observation"} or falsification["recommended_assay"] not in ASSAYS:
        raise ValueError("invalid falsification")
    uncertainty = card["uncertainty"]
    if set(uncertainty) != {"confidence", "abstain", "reason"} or not 0.0 <= uncertainty["confidence"] <= 1.0 or not isinstance(uncertainty["abstain"], bool):
        raise ValueError("invalid uncertainty")
    evidence = card["evidence"]
    if set(evidence) != {"top_patch_or_slide_evidence_ids", "programme_scores"}:
        raise ValueError("invalid evidence")
    if not isinstance(evidence["top_patch_or_slide_evidence_ids"], list) or len(evidence["top_patch_or_slide_evidence_ids"]) > MAX_EVIDENCE_IDS:
        raise ValueError("invalid evidence IDs")
    if not all(isinstance(value, str) and value for value in evidence["top_patch_or_slide_evidence_ids"]):
        raise ValueError("invalid evidence IDs")
    if set(evidence["programme_scores"]) - seen:
        raise ValueError("evidence scores must correspond to predicted programmes")


def verbalize_hypothesis_card(card: Mapping[str, Any]) -> str:
    """Safely render a validated card without adding external facts or entities."""
    validate_hypothesis_card(card)
    claims = card["predicted_programmes"]
    if claims:
        evidence = "; ".join(f"{item['programme']} is {item['direction']}" for item in claims)
    else:
        evidence = "No supported programme claim was emitted"
    hypothesis = card["composite_hypothesis"]["template"].replace("_", " ")
    assay = card["falsification"]["recommended_assay"].replace("_", " ")
    return (
        f"Research hypothesis only. WSI-derived programme evidence: {evidence}. "
        f"Composite hypothesis: {hypothesis}. Falsification assay: {assay}. "
        f"Expected observation: {card['falsification']['predicted_observation']} "
        f"Uncertainty: {card['uncertainty']['reason']}"
    )


def evaluate_hypothesis_cards(
    cards: Sequence[Mapping[str, Any]],
    patient_ids: Sequence[str],
    target_scores: Any,
    programme_names: Sequence[str],
    *,
    k: int = 3,
    threshold: float = 0.5,
) -> dict[str, float]:
    """Evaluate generated WSI-only claims against hidden programme targets.

    Relevance is the top-``k`` absolute target programmes for each patient;
    direction accuracy is computed only for emitted claims with an available
    target.  The function deliberately returns coverage so sparse panels cannot
    masquerade as zero performance.
    """
    if k < 1:
        raise ValueError("k must be positive")
    targets = _as_score_rows(patient_ids, target_scores, programme_names)
    by_id = {str(card["patient_id"]): card for card in cards}
    if set(by_id) != set(map(str, patient_ids)):
        raise ValueError("cards must contain exactly one card per patient")
    precision: list[float] = []
    recall: list[float] = []
    ndcg: list[float] = []
    direction_correct: list[float] = []
    confidence_outcomes: list[tuple[float, float]] = []
    composite_correct: list[float] = []
    for patient_id in map(str, patient_ids):
        card = by_id[patient_id]
        validate_hypothesis_card(card)
        available = {name: value for name, value in targets[patient_id].items() if name in PROGRAMMES and np.isfinite(value)}
        if not available:
            continue
        truth_ranked = sorted(available, key=lambda name: (-abs(available[name]), name))[:k]
        truth = set(truth_ranked)
        predicted = [item["programme"] for item in card["predicted_programmes"][:k] if item["programme"] in available]
        hits = [name in truth for name in predicted]
        precision.append(float(sum(hits) / max(1, len(predicted))))
        recall.append(float(sum(hits) / max(1, len(truth))))
        dcg = sum((1.0 if hit else 0.0) / log2(index + 2) for index, hit in enumerate(hits))
        ideal = sum(1.0 / log2(index + 2) for index in range(min(len(truth), len(predicted))))
        ndcg.append(float(dcg / ideal) if ideal else 0.0)
        for claim in card["predicted_programmes"]:
            programme = claim["programme"]
            if programme in available:
                correct = float(claim["direction"] == _direction(available[programme], threshold))
                direction_correct.append(correct)
                confidence_outcomes.append((float(claim["confidence"]), correct))
        truth_directions = {name: _direction(value, threshold) for name, value in available.items()}
        expected = _mechanism_for(truth_directions, available)
        composite_correct.append(float(card["composite_hypothesis"]["template"] == (expected.name if expected else "indeterminate")))
    def mean(values: Sequence[float]) -> float:
        return float(np.mean(values)) if values else float("nan")
    brier = mean([(confidence - outcome) ** 2 for confidence, outcome in confidence_outcomes])
    return {
        "n_cards": float(len(cards)), "n_target_covered": float(len(precision)),
        "programme_precision_at_k": mean(precision), "programme_recall_at_k": mean(recall),
        "programme_ndcg_at_k": mean(ndcg), "direction_accuracy": mean(direction_correct),
        "confidence_brier": brier, "composite_accuracy": mean(composite_correct),
    }
