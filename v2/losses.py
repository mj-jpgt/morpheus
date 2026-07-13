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


def gaussian_nll(mean: torch.Tensor, log_variance: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    return 0.5 * (log_variance + (target - mean).square() * torch.exp(-log_variance)).mean()
