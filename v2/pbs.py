"""Reference-context Perturbation-Basis Supervision primitives.

These primitives deliberately operate on dictionary *quotient coordinates*,
not individual CRISPRi atoms.  Highly correlated reference perturbations are
therefore represented together until a stability analysis supports a sharper
attribution.  This module supplies the one-seed feasibility pilot; it does not
make a patient-specific intervention prediction.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold


def _matrix(name: str, value: np.ndarray) -> np.ndarray:
    output = np.asarray(value, dtype=np.float64)
    if output.ndim != 2 or min(output.shape) == 0 or not np.isfinite(output).all():
        raise ValueError(f"{name} must be a finite nonempty [row, column] matrix")
    return output


@dataclass(frozen=True)
class ReferenceDictionary:
    """Frozen low-rank quotient of reference perturbation response signatures."""
    genes: tuple[str, ...]
    atom_ids: tuple[str, ...]
    gene_mean: np.ndarray
    gene_basis: np.ndarray  # [gene, quotient_component]
    singular_values: np.ndarray
    atom_coordinates: np.ndarray  # [atom, quotient_component]

    @classmethod
    def fit(cls, responses: np.ndarray, genes: Sequence[str], atom_ids: Sequence[str], *, n_components: int) -> "ReferenceDictionary":
        response = _matrix("responses", responses)
        genes, atom_ids = tuple(map(str, genes)), tuple(map(str, atom_ids))
        if len(genes) != response.shape[1] or len(atom_ids) != response.shape[0]:
            raise ValueError("dictionary identifiers do not match response matrix")
        if len(set(genes)) != len(genes) or len(set(atom_ids)) != len(atom_ids):
            raise ValueError("dictionary genes and atoms must be unique")
        mean = response.mean(axis=0, keepdims=True)
        _, singular, right_t = np.linalg.svd(response - mean, full_matrices=False)
        rank = int(np.sum(singular > singular[0] * 1e-8)) if singular.size and singular[0] > 0 else 0
        width = min(int(n_components), rank)
        if width < 1:
            raise ValueError("reference dictionary has zero numerical rank")
        basis = right_t[:width].T
        coordinates = (response - mean) @ basis
        return cls(genes, atom_ids, mean.ravel(), basis, singular[:width], coordinates)

    @property
    def components(self) -> int:
        return int(self.gene_basis.shape[1])

    def encode_expression(self, expression: np.ndarray, genes: Sequence[str]) -> np.ndarray:
        """Encode RNA without viewing WSI; exact gene identity is mandatory."""
        value = _matrix("expression", expression)
        observed = tuple(map(str, genes))
        if observed != self.genes:
            raise ValueError("RNA genes must exactly match frozen dictionary gene order")
        return (value - self.gene_mean) @ self.gene_basis


@dataclass(frozen=True)
class LegibilityOperator:
    """Frozen diagonal operator on quotient coordinates, estimated on a train fold."""
    weights: np.ndarray
    alpha: float
    fold_score: np.ndarray

    @classmethod
    def fit(cls, wsi: np.ndarray, codes: np.ndarray, cancers: Sequence[str], *,
            alphas=(1., 10., 100., 1e3, 1e4, 1e5, 1e6)) -> "LegibilityOperator":
        """Estimate per-axis legibility by grouped-CV held-out R2.

        The alpha grid must reach far beyond unity. `wsi` is a 1536-dimensional
        standardised patient vector against ~3.7k development rows, so at
        alpha<=10 the ridge is effectively unregularised: measured on the real
        cohort, EVERY Hallmark and PBS axis returned a NEGATIVE held-out R2
        (best -1.28), was clipped to zero, and the run aborted with "no nonzero
        development-fitted axes". At alpha=1e4 -- an interior optimum, selected
        by CV rather than assumed -- 45 of 50 Hallmark axes score positive.
        """
        x, y = _matrix("wsi", wsi), _matrix("codes", codes)
        groups = np.asarray(cancers).astype(str)
        if len(x) != len(y) or len(x) != len(groups) or len(np.unique(groups)) < 2:
            raise ValueError("legibility operator requires aligned rows from at least two cancer groups")
        folds = GroupKFold(n_splits=min(3, len(np.unique(groups))))
        errors: dict[float, list[float]] = {float(alpha): [] for alpha in alphas}
        per_component: dict[float, list[np.ndarray]] = {float(alpha): [] for alpha in alphas}
        for fit, valid in folds.split(x, y, groups):
            for alpha in errors:
                model = Ridge(alpha=alpha).fit(x[fit], y[fit])
                prediction = model.predict(x[valid])
                errors[alpha].append(float(np.mean((prediction - y[valid]) ** 2)))
                denom = np.var(y[valid], axis=0)
                score = 1.0 - np.mean((prediction - y[valid]) ** 2, axis=0) / np.maximum(denom, 1e-8)
                per_component[alpha].append(score)
        alpha = min(errors, key=lambda value: (float(np.mean(errors[value])), value))
        score = np.mean(per_component[alpha], axis=0)
        # Negative held-out R2 is not a visible coordinate.  The square root
        # makes W^(1/2) the multiplier consumed by the supervised loss.
        weights = np.clip(score, 0.0, 1.0)
        return cls(weights=weights, alpha=float(alpha), fold_score=score)

    def transform(self, code: np.ndarray) -> np.ndarray:
        return _matrix("code", code) * np.sqrt(self.weights)[None, :]


def weighted_code_loss(prediction, target, operator: LegibilityOperator):
    """MSE in the frozen, train-estimated legibility metric."""
    import torch
    if prediction.shape != target.shape or prediction.ndim != 2:
        raise ValueError("PBS prediction and target must be aligned [patient, component]")
    weights = torch.as_tensor(operator.weights, dtype=prediction.dtype, device=prediction.device).sqrt()
    return torch.nn.functional.mse_loss(prediction * weights, target * weights)
