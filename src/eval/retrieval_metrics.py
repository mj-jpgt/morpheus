"""Retrieval metrics for paired WSI-RNA alignment."""

from __future__ import annotations

import numpy as np


def l2_normalize(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    denom = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.maximum(denom, eps)


def paired_retrieval_metrics(
    query: np.ndarray,
    reference: np.ndarray,
    k_values: tuple[int, ...] = (1, 5, 10),
    query_cancers: list[str] | None = None,
    reference_cancers: list[str] | None = None,
) -> dict[str, float]:
    """Compute paired retrieval assuming row i in query matches row i in reference."""
    if len(query) != len(reference):
        raise ValueError("Paired retrieval requires equal row counts")
    q = l2_normalize(query)
    r = l2_normalize(reference)
    scores = q @ r.T
    order = np.argsort(-scores, axis=1)
    ranks = np.empty(len(q), dtype=np.int64)
    for i in range(len(q)):
        ranks[i] = int(np.where(order[i] == i)[0][0]) + 1
    out: dict[str, float] = {
        "n": float(len(q)),
        "mrr": float(np.mean(1.0 / ranks)) if len(ranks) else 0.0,
        "median_rank": float(np.median(ranks)) if len(ranks) else 0.0,
        "matched_cosine_mean": float(np.mean(np.diag(scores))) if len(q) else 0.0,
    }
    mask = ~np.eye(len(q), dtype=bool)
    out["unmatched_cosine_mean"] = float(np.mean(scores[mask])) if mask.any() else 0.0
    for k in k_values:
        out[f"recall_at_{k}"] = float(np.mean(ranks <= k)) if len(ranks) else 0.0
    if query_cancers is not None and reference_cancers is not None:
        same = []
        for i, cancer in enumerate(query_cancers):
            top = order[i, : min(10, order.shape[1])]
            same.append(any(reference_cancers[j] == cancer for j in top if j != i))
        out["same_cancer_in_top10"] = float(np.mean(same)) if same else 0.0
    return out


def unpaired_retrieval_metrics(
    query: np.ndarray,
    reference: np.ndarray,
    true_reference_index: np.ndarray,
    k_values: tuple[int, ...] = (1, 5, 10),
) -> dict[str, float]:
    """Compute retrieval where each query has a known row index in reference."""
    q = l2_normalize(query)
    r = l2_normalize(reference)
    true_reference_index = np.asarray(true_reference_index, dtype=np.int64)
    if len(q) != len(true_reference_index):
        raise ValueError("true_reference_index must have one entry per query")
    scores = q @ r.T
    order = np.argsort(-scores, axis=1)
    ranks = np.empty(len(q), dtype=np.int64)
    for i, target in enumerate(true_reference_index):
        ranks[i] = int(np.where(order[i] == target)[0][0]) + 1
    out: dict[str, float] = {
        "n": float(len(q)),
        "mrr": float(np.mean(1.0 / ranks)) if len(ranks) else 0.0,
        "median_rank": float(np.median(ranks)) if len(ranks) else 0.0,
        "matched_cosine_mean": float(np.mean(scores[np.arange(len(q)), true_reference_index])) if len(q) else 0.0,
    }
    for k in k_values:
        out[f"recall_at_{k}"] = float(np.mean(ranks <= k)) if len(ranks) else 0.0
    return out
