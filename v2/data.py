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

    def __init__(self, indices_or_token_counts: Sequence[int], patch_counts_or_token_budget: Sequence[int] | int | None = None,
                 token_budget: int | None = None, seed: int = 42) -> None:
        """Create an uncapped dynamic sampler.

        Preferred form is ``DynamicTokenBatchSampler(token_counts, budget,
        seed)``.  The historical public form ``DynamicTokenBatchSampler(
        indices, patch_counts, token_budget=budget, seed=seed)`` remains
        supported for callers that own non-contiguous patient indices.
        """
        second_is_counts = patch_counts_or_token_budget is not None and not np.isscalar(patch_counts_or_token_budget)
        if second_is_counts:
            if token_budget is None:
                raise TypeError("legacy indices/patch_counts form requires token_budget")
            self.indices = np.asarray(indices_or_token_counts, dtype=np.int64)
            self.token_counts = np.asarray(patch_counts_or_token_budget, dtype=np.int64)
            budget = int(token_budget)
        else:
            if patch_counts_or_token_budget is None:
                raise TypeError("token budget is required")
            self.indices = np.arange(len(indices_or_token_counts), dtype=np.int64)
            self.token_counts = np.asarray(indices_or_token_counts, dtype=np.int64)
            budget = int(patch_counts_or_token_budget)
            # Backward-compatible third positional argument is the seed in
            # the preferred form, e.g. (counts, budget, seed).
            if token_budget is not None:
                seed = int(token_budget)
        if len(self.indices) != len(self.token_counts):
            raise ValueError("indices and token counts must have equal length")
        if budget < 1 or np.any(self.token_counts < 1):
            raise ValueError("token budget and all patient token counts must be positive")
        self.token_budget = budget
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
                yield self.indices[np.asarray(batch, dtype=np.int64)]
                batch, used = [], 0
            batch.append(int(index))
            used += count
        if batch:
            yield self.indices[np.asarray(batch, dtype=np.int64)]

    def __iter__(self) -> Iterator[np.ndarray]:
        return self.batches(self.seed, shuffle=True)

    def coverage(self, epoch: int | None = None) -> np.ndarray | list[int]:
        """Return per-input coverage; no-argument legacy calls return a list."""
        seen = np.concatenate(list(self.batches(self.seed if epoch is None else epoch)))
        lookup = {int(value): row for row, value in enumerate(self.indices)}
        local = np.asarray([lookup[int(value)] for value in seen], dtype=np.int64)
        result = np.bincount(local, minlength=len(self.token_counts))
        return result.tolist() if epoch is None else result


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
