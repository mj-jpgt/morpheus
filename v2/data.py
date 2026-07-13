"""Leakage-safe patient records and dynamic uncapped patch batching for V2."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Sequence

import numpy as np


@dataclass(frozen=True)
class PatientRecord:
    patient_id: str
    cancer: str
    split: str
    patch_count: int
    has_wsi: bool
    has_rna: bool
    has_clinical: bool = False
    has_snv: bool = False
    has_cnv: bool = False


class DynamicTokenBatchSampler:
    """Batch complete patient bags without a patch cap or dropped examples."""

    def __init__(self, token_counts: Sequence[int], token_budget: int, seed: int = 42) -> None:
        self.token_counts = np.asarray(token_counts, dtype=np.int64)
        if token_budget < 1 or np.any(self.token_counts < 1):
            raise ValueError("token budget and all patient token counts must be positive")
        self.token_budget = int(token_budget)
        self.seed = int(seed)

    def batches(self, epoch: int, shuffle: bool = True) -> Iterator[np.ndarray]:
        order = np.argsort(self.token_counts, kind="stable")
        if shuffle:
            rng = np.random.default_rng(self.seed + int(epoch))
            chunks = []
            for start in range(0, len(order), 64):
                chunk = order[start : start + 64].copy()
                rng.shuffle(chunk)
                chunks.append(chunk)
            order = np.concatenate(chunks) if chunks else order
        batch: list[int] = []
        used = 0
        for index in order:
            count = int(self.token_counts[index])
            if batch and used + count > self.token_budget:
                yield np.asarray(batch, dtype=np.int64)
                batch, used = [], 0
            batch.append(int(index))
            used += count
        if batch:
            yield np.asarray(batch, dtype=np.int64)

    def coverage(self, epoch: int) -> np.ndarray:
        seen = np.concatenate(list(self.batches(epoch)))
        return np.bincount(seen, minlength=len(self.token_counts))


@dataclass
class TrainOnlyStandardizer:
    """Transform with a recorded fit population for leakage audits."""

    mean: np.ndarray | None = None
    scale: np.ndarray | None = None
    fit_indices: tuple[int, ...] = ()

    def fit(self, values: np.ndarray, indices: np.ndarray) -> "TrainOnlyStandardizer":
        if len(indices) == 0:
            raise ValueError("cannot fit on an empty population")
        values = np.asarray(values, dtype=np.float32)
        self.mean = values[indices].mean(axis=0, keepdims=True)
        self.scale = values[indices].std(axis=0, keepdims=True)
        self.scale[self.scale < 1e-6] = 1.0
        self.fit_indices = tuple(map(int, indices))
        return self

    def transform(self, values: np.ndarray) -> np.ndarray:
        if self.mean is None or self.scale is None:
            raise RuntimeError("fit must be called before transform")
        return ((np.asarray(values, dtype=np.float32) - self.mean) / self.scale).astype(np.float32)
