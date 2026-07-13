"""Leakage-safe prompted label and frozen-prototype cancer classification."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import torch
from torch import nn


DEFAULT_TEXT_MODEL = "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext"


def l2_normalize(values: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    return values / np.maximum(np.linalg.norm(values, axis=-1, keepdims=True), eps)


@dataclass(frozen=True)
class CancerDescription:
    cancer: str
    text: str
    source: str
    public_knowledge_only: bool = True


def validate_descriptions(descriptions: Iterable[CancerDescription], expected_cancers: Iterable[str]) -> dict[str, CancerDescription]:
    indexed = {item.cancer: item for item in descriptions}
    missing = sorted(set(expected_cancers) - set(indexed))
    if missing:
        raise ValueError(f"missing descriptions for: {missing}")
    if any(not item.public_knowledge_only or not item.source for item in indexed.values()):
        raise ValueError("every description must declare a public source and exclude patient-level data")
    return indexed


def embed_descriptions(descriptions: list[CancerDescription], output_path: str | Path, model_id: str = DEFAULT_TEXT_MODEL, revision: str | None = None, max_length: int = 256) -> Path:
    """Embed once and cache. Transformers is intentionally imported only here."""
    try:
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:  # pragma: no cover - dependency is optional in unit tests
        raise RuntimeError("install transformers and safetensors before generating description embeddings") from exc
    tokenizer = AutoTokenizer.from_pretrained(model_id, revision=revision)
    model = AutoModel.from_pretrained(model_id, revision=revision).eval()
    texts = [item.text for item in descriptions]
    encoded = tokenizer(texts, padding=True, truncation=True, max_length=max_length, return_tensors="pt")
    with torch.no_grad():
        hidden = model(**encoded).last_hidden_state
    weights = encoded["attention_mask"].unsqueeze(-1)
    vectors = (hidden * weights).sum(1) / weights.sum(1).clamp_min(1)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, cancers=np.asarray([item.cancer for item in descriptions]), embeddings=l2_normalize(vectors.cpu().numpy().astype(np.float32)), texts=np.asarray(texts), sources=np.asarray([item.source for item in descriptions]), model_id=model_id, revision="" if revision is None else revision)
    return output


class TextPrototypeMapper(nn.Module):
    """Small seen-cancer-only mapper from patient state to frozen text space."""

    def __init__(self, input_dim: int, text_dim: int, temperature: float = 0.07) -> None:
        super().__init__()
        self.mapper = nn.Sequential(nn.Linear(input_dim, text_dim), nn.LayerNorm(text_dim), nn.GELU(), nn.Linear(text_dim, text_dim))
        self.log_temperature = nn.Parameter(torch.tensor(float(np.log(temperature))))

    def logits(self, patient: torch.Tensor, text_prototypes: torch.Tensor) -> torch.Tensor:
        patient = nn.functional.normalize(self.mapper(patient), dim=-1)
        text_prototypes = nn.functional.normalize(text_prototypes, dim=-1)
        return patient @ text_prototypes.T / self.log_temperature.exp().clamp(0.01, 1.0)


def zero_shot_predict(patient_embeddings: np.ndarray, text_prototypes: np.ndarray, alpha: float = 1.0, support_prototypes: np.ndarray | None = None) -> np.ndarray:
    """Predict class indices from text-only or fixed text/support blended prototypes."""
    prototypes = text_prototypes if support_prototypes is None else alpha * text_prototypes + (1.0 - alpha) * support_prototypes
    return l2_normalize(patient_embeddings) @ l2_normalize(prototypes).T


def make_support_prototypes(embeddings: np.ndarray, labels: np.ndarray, classes: np.ndarray, support_indices: np.ndarray) -> np.ndarray:
    result = []
    for label in classes:
        selected = support_indices[labels[support_indices] == label]
        if len(selected) == 0:
            raise ValueError(f"support set has no examples for class {label}")
        result.append(embeddings[selected].mean(axis=0))
    return l2_normalize(np.vstack(result).astype(np.float32))


def few_shot_episodes(embeddings: np.ndarray, labels: np.ndarray, classes: np.ndarray, k_values: tuple[int, ...] = (1, 2, 5, 10, 20), episodes: int = 100, seed: int = 42) -> list[dict[str, object]]:
    """Generate reproducible, disjoint support/query episodes for frozen encoders."""
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for k in k_values:
        for episode in range(episodes):
            support: list[int] = []
            feasible = True
            for label in classes:
                candidates = np.where(labels == label)[0]
                if len(candidates) <= k:
                    feasible = False
                    break
                support.extend(rng.choice(candidates, size=k, replace=False).tolist())
            if not feasible:
                continue
            support_array = np.asarray(sorted(support), dtype=np.int64)
            query = np.setdiff1d(np.arange(len(labels)), support_array, assume_unique=True)
            prototypes = make_support_prototypes(embeddings, labels, classes, support_array)
            scores = zero_shot_predict(embeddings[query], prototypes)
            predicted = classes[np.argmax(scores, axis=1)]
            rows.append({"k": k, "episode": episode, "support_indices": support_array, "query_indices": query, "accuracy": float((predicted == labels[query]).mean())})
    return rows
