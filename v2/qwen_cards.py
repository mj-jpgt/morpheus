"""Closed-RAG Qwen rendering for MORPHEUS research hypothesis cards.

Qwen is optional and never has access to paired RNA, outcome labels, external
tools, or arbitrary prompt text.  The deterministic card is the only patient
input.  The model may rank pre-approved claims and cite retrieved local
evidence, but a strict validator rejects any new programme, mechanism, assay,
or evidence reference.  A deterministic abstention is emitted on all failures.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any, Mapping, Protocol, Sequence

from .evidence_registry import EvidenceChunk, EvidenceRegistry
from .hypothesis_cards import _as_score_rows, _direction, validate_hypothesis_card, verbalize_hypothesis_card


QWEN_CARD_VERSION = "morpheus_qwen_card_v1"
QWEN_MODE = "research_hypothesis_only"
_PROHIBITED_TERMS = frozenset((
    "treat", "treatment", "therapy", "therapeutic", "prescribe", "diagnos", "clinical decision",
    "causes", "causal", "prognosis", "survival", "response to", "drug recommendation",
))
_NOVEL_RELATIONS = frozenset(("co_occurrence", "co_localisation", "spatial_separation", "inverse_association"))
_NOVEL_SAFE_TERMS = frozenset((
    "the", "and", "with", "is", "are", "may", "be", "associated", "observed", "in", "at", "within",
    "tumour", "tumor", "regions", "state", "states", "programme", "programmes", "program", "programs",
    "high", "low", "evidence", "suggests", "consistent", "for", "a", "an", "of", "to", "by", "between",
    "co", "occurrence", "localisation", "localization", "spatial", "separation", "inverse", "association",
))


class QwenBackend(Protocol):
    """Minimal backend contract; enables deterministic tests without a model."""

    def generate(self, prompt: str, *, max_new_tokens: int, temperature: float) -> str: ...


@dataclass(frozen=True)
class QwenGenerationConfig:
    model_id: str = "Qwen/Qwen2.5-7B-Instruct"
    max_new_tokens: int = 384
    temperature: float = 0.0
    evidence_limit: int = 4
    retry_once: bool = True
    local_files_only: bool = True


class TransformersQwenBackend:
    """Lazy local Qwen backend; no import, download, or GPU allocation at import time."""

    def __init__(self, config: QwenGenerationConfig = QwenGenerationConfig()):
        self.config = config
        try:
            import torch
            from transformers import AutoModelForCausalLM, AutoTokenizer
        except ImportError as exc:  # pragma: no cover - exercised on Lambda only
            raise RuntimeError("Qwen runtime requires torch and transformers") from exc
        self._torch = torch
        self._tokenizer = AutoTokenizer.from_pretrained(config.model_id, local_files_only=config.local_files_only)
        self._model = AutoModelForCausalLM.from_pretrained(
            config.model_id, torch_dtype=torch.bfloat16, device_map="auto", local_files_only=config.local_files_only,
        ).eval()

    def generate(self, prompt: str, *, max_new_tokens: int, temperature: float) -> str:  # pragma: no cover - GPU runtime
        messages = [{"role": "system", "content": "Return only the requested JSON object."}, {"role": "user", "content": prompt}]
        rendered = self._tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=True)
        inputs = self._tokenizer(rendered, return_tensors="pt").to(self._model.device)
        with self._torch.inference_mode():
            output = self._model.generate(
                **inputs, max_new_tokens=max_new_tokens, do_sample=temperature > 0,
                temperature=temperature if temperature > 0 else None, pad_token_id=self._tokenizer.eos_token_id,
            )
        return self._tokenizer.decode(output[0, inputs["input_ids"].shape[1]:], skip_special_tokens=True).strip()


def _extract_json(text: str) -> dict[str, Any]:
    """Accept one JSON object, optionally fenced, while rejecting surrounding prose."""
    stripped = text.strip()
    if stripped.startswith("```"):
        stripped = re.sub(r"^```(?:json)?\s*|\s*```$", "", stripped, flags=re.IGNORECASE)
    try:
        value = json.loads(stripped)
    except json.JSONDecodeError as exc:
        raise ValueError("Qwen response is not one valid JSON object") from exc
    if not isinstance(value, dict):
        raise ValueError("Qwen response must be a JSON object")
    return value


def _allowed_from_card(card: Mapping[str, Any]) -> tuple[dict[str, Mapping[str, Any]], str, str, tuple[str, ...]]:
    validate_hypothesis_card(card)
    allowed = {str(item["programme"]): item for item in card["predicted_programmes"]}
    mechanism = str(card["composite_hypothesis"]["template"])
    assay = str(card["falsification"]["recommended_assay"])
    patient_evidence = tuple(map(str, card["evidence"]["top_patch_or_slide_evidence_ids"]))
    return allowed, mechanism, assay, patient_evidence


def build_qwen_request(card: Mapping[str, Any], registry: EvidenceRegistry, *, evidence_limit: int = 4) -> dict[str, Any]:
    """Create a public, label-free request. Hidden RNA/outcome data are not accepted."""
    allowed, mechanism, assay, patch_ids = _allowed_from_card(card)
    chunks = registry.retrieve(list(allowed), assay=assay, limit=evidence_limit)
    return {
        "schema_version": QWEN_CARD_VERSION,
        "mode": QWEN_MODE,
        "patient_id": str(card["patient_id"]),
        "allowed_programme_claims": [
            {"programme": name, "direction": item["direction"], "confidence": float(item["confidence"])}
            for name, item in allowed.items()
        ],
        "allowed_mechanism": mechanism,
        "allowed_assay": assay,
        "patient_evidence_ids": list(patch_ids),
        "registry_version": registry.version,
        "registry_digest": registry.digest,
        "retrieved_evidence": [
            {"id": chunk.id, "title": chunk.title, "text": chunk.text,
             "supported_programmes": list(chunk.supported_programmes), "permitted_assays": list(chunk.permitted_assays)}
            for chunk in chunks
        ],
    }


def build_qwen_prompt(request: Mapping[str, Any], *, repair_error: str | None = None) -> str:
    """Construct the closed-world Qwen prompt; no free clinical wording is requested."""
    if repair_error:
        repair = f"Your previous JSON was rejected: {repair_error}. Correct it without adding information.\n"
    else:
        repair = ""
    return (
        "You are a constrained research-hypothesis formatter. This is NOT clinical decision support. "
        "Use only the allowed claims, mechanism, assay, and retrieved evidence IDs below. Do not introduce "
        "entities, genes, diseases, treatments, causal statements, prognosis, or patient facts. Return only one JSON object "
        "with exactly these keys: schema_version, mode, patient_id, abstained, programme_claims, mechanism, falsification, "
        "evidence_references, novel_hypothesis, uncertainty, abstention_reason. programme_claims must be an ordered subset of allowed programme names. "
        "mechanism must equal allowed_mechanism or indeterminate. falsification must equal allowed_assay or no_assay_recommended. "
        "evidence_references must be IDs from retrieved_evidence only. uncertainty must be one of low, moderate, high.\n"
        "novel_hypothesis must be null or an exploratory combination of two or more allowed programme names, a relation, "
        "retrieved evidence IDs, and a statement using only terms already present in the allowed programme names or evidence.\n"
        f"{repair}INPUT_JSON:\n{json.dumps(request, sort_keys=True)}"
    )


def _abstention(card: Mapping[str, Any], registry: EvidenceRegistry, *, reason: str) -> dict[str, Any]:
    _, _, _, patch_ids = _allowed_from_card(card)
    result = {
        "schema_version": QWEN_CARD_VERSION, "mode": QWEN_MODE, "patient_id": str(card["patient_id"]),
        "abstained": True, "programme_claims": [], "mechanism": "indeterminate",
        "falsification": "no_assay_recommended", "evidence_references": [],
        "novel_hypothesis": None,
        "uncertainty": "high", "abstention_reason": reason,
        "provenance": {"registry_version": registry.version, "registry_digest": registry.digest,
                       "patient_evidence_ids": list(patch_ids), "renderer": "validated_abstention"},
    }
    validate_qwen_card(result, card, registry)
    return result


def validate_qwen_card(result: Mapping[str, Any], source_card: Mapping[str, Any], registry: EvidenceRegistry) -> None:
    """Validate a Qwen card against its deterministic source and closed registry."""
    allowed, source_mechanism, source_assay, patch_ids = _allowed_from_card(source_card)
    expected = {
        "schema_version", "mode", "patient_id", "abstained", "programme_claims", "mechanism", "falsification",
        "evidence_references", "novel_hypothesis", "uncertainty", "abstention_reason", "provenance",
    }
    if set(result) != expected:
        raise ValueError("Qwen card has missing or unexpected fields")
    if result["schema_version"] != QWEN_CARD_VERSION or result["mode"] != QWEN_MODE:
        raise ValueError("invalid Qwen card metadata")
    if str(result["patient_id"]) != str(source_card["patient_id"]):
        raise ValueError("Qwen card patient mismatch")
    if not isinstance(result["abstained"], bool) or result["uncertainty"] not in {"low", "moderate", "high"}:
        raise ValueError("invalid abstention or uncertainty")
    if not isinstance(result["programme_claims"], list) or len(result["programme_claims"]) > len(allowed):
        raise ValueError("invalid Qwen programme claims")
    claims = result["programme_claims"]
    seen: set[str] = set()
    for claim in claims:
        if set(claim) != {"programme", "direction", "confidence"}:
            raise ValueError("invalid Qwen programme claim shape")
        programme = str(claim["programme"])
        if programme not in allowed or programme in seen:
            raise ValueError("Qwen invented an unsupported programme claim")
        if claim["direction"] != allowed[programme]["direction"] or float(claim["confidence"]) != float(allowed[programme]["confidence"]):
            raise ValueError("Qwen changed deterministic programme direction/confidence")
        seen.add(programme)
    permitted_mechanisms = {source_mechanism, "indeterminate"}
    if result["mechanism"] not in permitted_mechanisms:
        raise ValueError("Qwen invented an unsupported mechanism")
    permitted_assays = {source_assay, "no_assay_recommended"}
    if result["falsification"] not in permitted_assays:
        raise ValueError("Qwen invented an unsupported assay")
    refs = result["evidence_references"]
    if not isinstance(refs, list) or len(refs) > 4 or len(set(refs)) != len(refs):
        raise ValueError("invalid Qwen evidence references")
    for evidence_id in refs:
        chunk = registry.get(str(evidence_id))
        if not set(claim["programme"] for claim in claims) & set(chunk.supported_programmes) and result["falsification"] not in chunk.permitted_assays:
            raise ValueError("evidence does not support the Qwen card")
    novel = result["novel_hypothesis"]
    if novel is not None:
        if result["abstained"]:
            raise ValueError("abstained Qwen card cannot contain a novel hypothesis")
        if set(novel) != {"programme_refs", "relation", "evidence_refs", "statement"}:
            raise ValueError("invalid novel hypothesis shape")
        programmes = novel["programme_refs"]
        evidence_refs = novel["evidence_refs"]
        if not isinstance(programmes, list) or len(programmes) < 2 or len(programmes) != len(set(programmes)):
            raise ValueError("novel hypothesis needs at least two unique programme references")
        if not set(programmes).issubset(allowed):
            raise ValueError("novel hypothesis invented a programme")
        if novel["relation"] not in _NOVEL_RELATIONS:
            raise ValueError("novel hypothesis has unsupported relation")
        if not isinstance(evidence_refs, list) or not evidence_refs or not set(evidence_refs).issubset(refs):
            raise ValueError("novel hypothesis must cite selected evidence")
        supporting = set()
        for evidence_id in evidence_refs:
            supporting.update(registry.get(evidence_id).supported_programmes)
        if not set(programmes).issubset(supporting):
            raise ValueError("novel hypothesis evidence does not support all programmes")
        statement = novel["statement"]
        if not isinstance(statement, str) or not statement.strip() or any(term in statement.lower() for term in _PROHIBITED_TERMS):
            raise ValueError("unsafe or empty novel hypothesis statement")
        whitelist = set(_NOVEL_SAFE_TERMS)
        for programme in programmes:
            whitelist.update(re.findall(r"[a-z]+", programme.lower()))
        for evidence_id in evidence_refs:
            chunk = registry.get(evidence_id)
            whitelist.update(re.findall(r"[a-z]+", chunk.title.lower()))
            whitelist.update(re.findall(r"[a-z]+", chunk.text.lower()))
        statement_terms = set(re.findall(r"[a-z]+", statement.lower()))
        if statement_terms - whitelist:
            raise ValueError("novel hypothesis statement contains non-whitelisted entity or term")
    if result["abstained"] and (claims or result["mechanism"] != "indeterminate" or result["falsification"] != "no_assay_recommended"):
        raise ValueError("abstained Qwen card must make no claim")
    if not isinstance(result["abstention_reason"], str):
        raise ValueError("invalid Qwen abstention reason")
    provenance = result["provenance"]
    if set(provenance) != {"registry_version", "registry_digest", "patient_evidence_ids", "renderer"}:
        raise ValueError("invalid Qwen provenance")
    if provenance["registry_version"] != registry.version or provenance["registry_digest"] != registry.digest:
        raise ValueError("Qwen registry provenance mismatch")
    if list(provenance["patient_evidence_ids"]) != list(patch_ids) or provenance["renderer"] not in {"qwen_validated", "validated_abstention"}:
        raise ValueError("invalid Qwen source provenance")


def _normalise_qwen_response(raw: Mapping[str, Any], request: Mapping[str, Any], source_card: Mapping[str, Any], registry: EvidenceRegistry) -> dict[str, Any]:
    required_raw = {
        "schema_version", "mode", "patient_id", "abstained", "programme_claims", "mechanism", "falsification",
        "evidence_references", "novel_hypothesis", "uncertainty", "abstention_reason",
    }
    if set(raw) != required_raw:
        raise ValueError("Qwen response has missing or unexpected fields")
    allowed = {item["programme"]: item for item in request["allowed_programme_claims"]}
    values = raw.get("programme_claims", [])
    if not isinstance(values, list):
        raise ValueError("programme_claims must be a list")
    claims = []
    for programme in values:
        if not isinstance(programme, str) or programme not in allowed:
            raise ValueError("Qwen selected a non-allowed programme")
        item = allowed[programme]
        claims.append({"programme": programme, "direction": item["direction"], "confidence": item["confidence"]})
    abstention_reason = raw["abstention_reason"]
    if not isinstance(abstention_reason, str) or any(term in abstention_reason.lower() for term in _PROHIBITED_TERMS):
        raise ValueError("Qwen abstention reason contains prohibited clinical or causal language")
    result = {
        "schema_version": raw.get("schema_version"), "mode": raw.get("mode"), "patient_id": raw.get("patient_id"),
        "abstained": raw.get("abstained"), "programme_claims": claims, "mechanism": raw.get("mechanism"),
        "falsification": raw.get("falsification"), "evidence_references": raw.get("evidence_references", []),
        "novel_hypothesis": raw.get("novel_hypothesis"),
        "uncertainty": raw.get("uncertainty"), "abstention_reason": abstention_reason,
        "provenance": {"registry_version": registry.version, "registry_digest": registry.digest,
                       "patient_evidence_ids": request["patient_evidence_ids"], "renderer": "qwen_validated"},
    }
    validate_qwen_card(result, source_card, registry)
    return result


def render_qwen_card(card: Mapping[str, Any], registry: EvidenceRegistry, backend: QwenBackend | None, *, config: QwenGenerationConfig = QwenGenerationConfig()) -> dict[str, Any]:
    """Render one card with Qwen, retry once, otherwise return a safe abstention."""
    request = build_qwen_request(card, registry, evidence_limit=config.evidence_limit)
    if backend is None:
        return _abstention(card, registry, reason="Qwen backend unavailable; deterministic abstention.")
    error = ""
    attempts = 2 if config.retry_once else 1
    for attempt in range(attempts):
        try:
            raw = _extract_json(backend.generate(build_qwen_prompt(request, repair_error=error or None), max_new_tokens=config.max_new_tokens, temperature=config.temperature))
            return _normalise_qwen_response(raw, request, card, registry)
        except Exception as exc:  # validated safe fallback is intentional
            error = str(exc)[:300]
            if attempt + 1 == attempts:
                return _abstention(card, registry, reason=f"Qwen output rejected: {error}")
    raise AssertionError("unreachable")


def render_qwen_cards(cards: Sequence[Mapping[str, Any]], registry: EvidenceRegistry, backend: QwenBackend | None, *, config: QwenGenerationConfig = QwenGenerationConfig()) -> list[dict[str, Any]]:
    """Render a sequence independently so a malformed patient response cannot stop a run."""
    return [render_qwen_card(card, registry, backend, config=config) for card in cards]


def deterministic_qwen_summary(result: Mapping[str, Any], source_card: Mapping[str, Any], registry: EvidenceRegistry) -> str:
    """Render validated JSON without letting an LLM add untraceable prose."""
    validate_qwen_card(result, source_card, registry)
    if result["abstained"]:
        return "Insufficient WSI-derived evidence for a ranked research hypothesis; validate with additional molecular profiling."
    return verbalize_hypothesis_card(source_card)


def evaluate_qwen_cards(
    results: Sequence[Mapping[str, Any]],
    source_cards: Sequence[Mapping[str, Any]],
    patient_ids: Sequence[str],
    hidden_target_scores: Any,
    programme_names: Sequence[str],
    *,
    k: int = 3,
    threshold: float = 0.5,
    registry: EvidenceRegistry,
) -> dict[str, float]:
    """Score validated Qwen selections against hidden labels after generation.

    The target matrix intentionally enters only this evaluator.  It is not an
    argument to request construction, prompting, rendering, or validation.
    """
    if k < 1:
        raise ValueError("k must be positive")
    ids = [str(value) for value in patient_ids]
    source_by_id = {str(card["patient_id"]): card for card in source_cards}
    result_by_id = {str(card["patient_id"]): card for card in results}
    if set(source_by_id) != set(ids) or set(result_by_id) != set(ids):
        raise ValueError("source and Qwen cards must have exact patient coverage")
    targets = _as_score_rows(ids, hidden_target_scores, programme_names)
    precision: list[float] = []
    recall: list[float] = []
    ndcg: list[float] = []
    correct_confidence: list[tuple[float, float]] = []
    non_abstained = 0
    evidence_valid = 0
    for patient_id in ids:
        source = source_by_id[patient_id]
        result = result_by_id[patient_id]
        validate_qwen_card(result, source, registry)
        evidence_valid += 1
        if result["abstained"]:
            continue
        non_abstained += 1
        available = {name: value for name, value in targets[patient_id].items() if name in {x["programme"] for x in source["predicted_programmes"]}}
        if not available:
            continue
        truth_ranked = sorted(available, key=lambda name: (-abs(available[name]), name))[:k]
        truth = set(truth_ranked)
        predicted = result["programme_claims"][:k]
        hits = [item["programme"] in truth for item in predicted]
        precision.append(float(sum(hits) / max(1, len(predicted))))
        recall.append(float(sum(hits) / max(1, len(truth))))
        ideal = sum(1.0 / (index + 1) for index in range(min(len(truth), len(predicted))))
        ndcg.append(float(sum((1.0 if hit else 0.0) / (index + 1) for index, hit in enumerate(hits)) / ideal) if ideal else 0.0)
        for item in predicted:
            target = available.get(item["programme"])
            if target is not None:
                correct = float(item["direction"] == _direction(float(target), threshold))
                correct_confidence.append((float(item["confidence"]), correct))
    def _mean(values: Sequence[float]) -> float:
        return float(sum(values) / len(values)) if values else float("nan")
    return {
        "n_cards": float(len(ids)), "abstention_rate": float(1.0 - non_abstained / max(1, len(ids))),
        "non_abstained_cards": float(non_abstained), "evidence_validity": float(evidence_valid / max(1, len(ids))),
        "programme_precision_at_k": _mean(precision), "programme_recall_at_k": _mean(recall), "programme_ndcg_at_k": _mean(ndcg),
        "direction_accuracy": _mean([outcome for _, outcome in correct_confidence]),
        "confidence_brier": _mean([(confidence - outcome) ** 2 for confidence, outcome in correct_confidence]),
    }


def write_qwen_jsonl(records: Sequence[Mapping[str, Any]], path: str | Path) -> None:
    """Persist structured outputs only; prose is derived on demand after validation."""
    Path(path).write_text("".join(json.dumps(dict(record), sort_keys=True) + "\n" for record in records), encoding="utf-8")
