"""V2 objectives. Biology-state losses intentionally exclude paired CLIP."""

from __future__ import annotations

import torch
from torch import nn


def _norm(x: torch.Tensor) -> torch.Tensor:
    return nn.functional.normalize(x, dim=-1)


def symmetric_infonce(left: torch.Tensor, right: torch.Tensor, temperature: float = 0.07) -> torch.Tensor:
    if len(left) < 2:
        return left.new_zeros(())
    logits = _norm(left) @ _norm(right).T / temperature
    labels = torch.arange(len(left), device=left.device)
    return 0.5 * (nn.functional.cross_entropy(logits, labels) + nn.functional.cross_entropy(logits.T, labels))


def centered_cross_covariance(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    if len(left) < 2:
        return left.new_zeros(())
    left = left - left.mean(dim=0, keepdim=True)
    right = right - right.mean(dim=0, keepdim=True)
    return ((left.T @ right) / max(len(left) - 1, 1)).square().mean()


def variance_floor(state: torch.Tensor, target_std: float = 1.0, eps: float = 1e-4) -> torch.Tensor:
    """Keep a specialised state from satisfying separation by collapsing."""
    if len(state) < 2:
        return state.new_zeros(())
    std = torch.sqrt(state.var(dim=0, unbiased=False) + eps)
    return torch.relu(float(target_std) - std).mean()


def whitened_cross_covariance(left: torch.Tensor, right: torch.Tensor, eps: float = 1e-4) -> torch.Tensor:
    """Cross-covariance after per-dimension whitening, robust to scale tricks."""
    if len(left) < 2:
        return left.new_zeros(())
    left = (left - left.mean(dim=0, keepdim=True)) / torch.sqrt(left.var(dim=0, unbiased=False, keepdim=True) + eps)
    right = (right - right.mean(dim=0, keepdim=True)) / torch.sqrt(right.var(dim=0, unbiased=False, keepdim=True) + eps)
    return ((left.T @ right) / max(len(left) - 1, 1)).square().mean()


def feature_decorrelation(state: torch.Tensor, min_batch: int = 8, eps: float = 1e-4) -> torch.Tensor:
    """Off-diagonal feature-*correlation* penalty (Barlow-Twins style).

    Unlike :func:`variance_floor`, which acts per dimension and cannot raise the
    rank of a representation, this penalises *cross-dimension* correlation and so
    directly discourages rank collapse of the biology head onto the low-rank
    programme-target manifold.

    The features are standardised to unit variance BEFORE the off-diagonal
    penalty (as in :func:`whitened_cross_covariance`).  This is deliberate and
    load-bearing: ``z_biology`` is L2-normalised, so its raw per-feature variance
    is O(1/dim) and a plain covariance penalty is numerically negligible (it
    silently no-ops at any sane weight).  Standardising makes the term a
    scale-invariant correlation penalty of O(1).  The batch correlation is a
    rank-``(B-1)`` estimator, so it is skipped below ``min_batch`` where the
    estimate is unreliable on small, ragged patient batches
    (F-R2; see v2/research/B2_implementation_and_audit.md).
    """
    if len(state) < min_batch:
        return state.new_zeros(())
    standardized = (state - state.mean(dim=0, keepdim=True)) / torch.sqrt(
        state.var(dim=0, unbiased=False, keepdim=True) + eps)
    correlation = (standardized.T @ standardized) / (len(state) - 1)
    off_diagonal = correlation - torch.diag(torch.diag(correlation))
    return off_diagonal.square().sum() / correlation.shape[0]


def programme_neighbourhood_loss(state: torch.Tensor, targets: torch.Tensor, temperature: float = 0.20) -> torch.Tensor:
    """Match state similarity to train-fold programme similarity."""
    if len(state) < 3:
        return state.new_zeros(())
    mask = ~torch.eye(len(state), dtype=torch.bool, device=state.device)
    target_logits = (_norm(targets) @ _norm(targets).T).masked_fill(~mask, -1e4)
    state_logits = (_norm(state) @ _norm(state).T / temperature).masked_fill(~mask, -1e4)
    return nn.functional.kl_div(torch.log_softmax(state_logits, dim=-1), torch.softmax(target_logits, dim=-1), reduction="batchmean")


def supervised_programme_contrastive(state: torch.Tensor, positive_mask: torch.Tensor, temperature: float = 0.20) -> torch.Tensor:
    if len(state) < 3 or not positive_mask.any():
        return state.new_zeros(())
    valid = positive_mask.bool() & ~torch.eye(len(state), dtype=torch.bool, device=state.device)
    # Self-pairs are invalid; off-diagonal programme-similar patients are the
    # only legal positives.  Masking the complement would turn every positive
    # into -1e4 and dominate the entire objective.
    logits = (_norm(state) @ _norm(state).T / temperature).masked_fill(torch.eye(len(state), dtype=torch.bool, device=state.device), -1e4)
    log_prob = torch.log_softmax(logits, dim=-1)
    losses = [-log_prob[row, valid[row]].mean() for row in range(len(state)) if valid[row].any()]
    return torch.stack(losses).mean() if losses else state.new_zeros(())


def gaussian_nll(mean: torch.Tensor, log_variance: torch.Tensor, target: torch.Tensor,
                 target_mask: torch.Tensor | None = None) -> torch.Tensor:
    """Heteroscedastic Gaussian NLL with optional per-programme availability.

    A patient can legitimately lack one RNA-derived programme while retaining
    others.  Reducing over a row-level mask would either discard usable labels
    or silently treat unavailable targets as zeros.
    """
    value = 0.5 * (log_variance + (target - mean).square() * torch.exp(-log_variance))
    if target_mask is None:
        return value.mean()
    if target_mask.shape != value.shape:
        raise ValueError("target_mask must have the same [batch, programme] shape as Gaussian NLL inputs")
    mask = target_mask.bool()
    if not mask.any():
        return value.new_zeros(())
    return value.masked_select(mask).mean()
