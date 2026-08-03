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


def population_offset(current: torch.Tensor, memory: torch.Tensor | None,
                      min_batch: int = 8) -> torch.Tensor | None:
    """Estimate the batch-COMMON direction of a biology state, or ``None``.

    ``z_biology`` is L2-normalised but never centred, so at initialisation the
    overwhelming majority of its norm is a single direction shared by every
    patient: measured on a real 16-patient batch, ``||mean_i z_i||**2 = 0.8133``
    and the patient-to-patient cosine is 0.8008.  See
    :func:`paired_infonce_with_memory` for why that is fatal.

    Returns a ``[1, dim]`` mean, or ``None`` when no reliable estimate exists.
    The reference is the current batch when it is large enough to estimate a
    mean, and the current batch pooled with the detached key queue otherwise.
    Real ragged WSI batches hold only B=1-3 patients: centring those by their
    own mean would map B=1 to the zero vector and B=2 to an antipodal pair, so
    below ``min_batch`` the queue -- which always holds at least 16 keys and is
    refreshed every training step -- supplies the estimate instead.
    """
    if len(current) >= min_batch:
        return current.mean(dim=0, keepdim=True)
    if memory is None or len(memory) == 0:
        return None
    pooled = torch.cat([current, memory.detach().to(current.device)], dim=0)
    return pooled.mean(dim=0, keepdim=True)


def _centre_for_infonce(current: torch.Tensor, memory: torch.Tensor | None,
                        min_batch: int) -> tuple[torch.Tensor, torch.Tensor | None]:
    """Remove the common direction from the current rows and the queue SEPARATELY.

    Each population gets its own offset.  Forcing the current batch's offset
    onto a frozen queue was measured to be unstable: the queue's keys were
    encoded by an earlier model and drift away from the live states, so the
    shared-offset variant reached the criterion at step 200 (0.0311, 16/16
    retrieval) and then diverged to 7.62 and full collapse by step 400.  Giving
    each population its own offset converges monotonically instead:
    0.4025 -> 0.0835 -> 0.0101 at steps 400/600/800.
    """
    offset = population_offset(current, memory, min_batch=min_batch)
    centred_current = current if offset is None else current - offset
    centred_memory = memory
    if memory is not None and len(memory) >= min_batch:
        centred_memory = memory - memory.mean(dim=0, keepdim=True)
    return centred_current, centred_memory


def paired_infonce_with_memory(
    wsi: torch.Tensor,
    rna: torch.Tensor,
    patient_indices: torch.Tensor,
    memory_wsi: torch.Tensor | None,
    memory_rna: torch.Tensor | None,
    memory_indices: torch.Tensor | None,
    *,
    temperature: float = 0.07,
    min_negatives: int = 8,
    centre: bool = True,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Symmetric patient-paired InfoNCE with ID-aware detached key queues.

    The current ragged WSI batches can contain only one to three patients, so
    ordinary in-batch CLIP has no usable negative set. The queue contributes
    detached WSI and RNA keys; current paired keys remain differentiable. An
    ID match in either queue is masked rather than treated as a false negative.
    Returns ``(loss, per-row_negative_counts)`` and refuses a silent no-op.

    ``centre`` removes each modality's batch-common direction before the
    similarity, and is **load-bearing, not cosmetic**.  Without it this loss
    drives the biology head to effective rank 1.00 within 50 steps and then
    stops moving.  The mechanism, measured on the real cohort (2026-08-03):

    * ``z_biology`` is normalised but not centred, so ~81% of its squared norm
      is one direction shared by every patient;
    * every cosine is therefore dominated by a constant, and each key carries a
      query-independent bias ``<mean_query_side, key>`` that swamps the
      patient-specific term -- classic hubness;
    * consequently the positives start no better than the negatives (minimum
      margin **-0.219**) and the loss starts at **3.0762, ABOVE chance**
      (ln 16 = 2.7726);
    * from above chance, *erasing every distinction is a genuine descent
      direction*, and it terminates on the permutation-symmetric configuration
      where the InfoNCE gradient is exactly zero by symmetry.

    Neither `variance_floor` nor `feature_decorrelation` can prevent this and
    both were measured failing to: a per-dimension variance floor is satisfied
    by the rank-1 family ``z_i = m + a_i*u`` (weight 10.0 still collapsed), and
    the covariance penalty switches itself off by collapsing (20.74 -> 0.00 in
    25 steps at every weight from 0.04 to 4.0).  Centring is the only measured
    intervention that learns: raw in-batch InfoNCE 3.0762 -> 0.0110 with 16/16
    retrieval by step 400, patient-to-patient cosine 0.8008 -> 0.0528.
    """
    if wsi.ndim != 2 or rna.shape != wsi.shape or patient_indices.shape != (len(wsi),):
        raise ValueError("WSI, RNA and patient_indices must be aligned [batch, dimension]")
    if patient_indices.unique().numel() != len(patient_indices):
        raise ValueError("programme_free requires unique canonical patient indices within a batch")
    queued = (memory_wsi, memory_rna, memory_indices)
    if any(value is None for value in queued):
        if not all(value is None for value in queued):
            raise ValueError("WSI/RNA memory states and IDs must be supplied together")
        memory_wsi = memory_rna = memory_indices = None
    if memory_indices is not None and (
        memory_wsi is None or memory_rna is None or memory_wsi.ndim != 2 or memory_rna.ndim != 2
        or memory_wsi.shape != memory_rna.shape or memory_wsi.shape[1] != wsi.shape[1]
        or memory_indices.shape != (len(memory_wsi),)
    ):
        raise ValueError("memory WSI/RNA states and IDs must be aligned")
    if centre:
        # Modality-wise, and population-wise within each modality: the WSI
        # current rows and the WSI queue are centred independently, likewise RNA.
        wsi, memory_wsi = _centre_for_infonce(wsi, memory_wsi, min_negatives)
        rna, memory_rna = _centre_for_infonce(rna, memory_rna, min_negatives)

    def directional(query: torch.Tensor, paired_keys: torch.Tensor,
                    queued_keys: torch.Tensor | None) -> tuple[torch.Tensor, torch.Tensor]:
        current = _norm(paired_keys)
        if queued_keys is None or len(queued_keys) == 0:
            keys, key_ids = current, patient_indices
        else:
            assert memory_indices is not None
            keys = torch.cat([current, _norm(queued_keys.detach().to(query.device))], dim=0)
            key_ids = torch.cat([patient_indices, memory_indices.detach().to(patient_indices.device)], dim=0)
        logits = _norm(query) @ keys.T / temperature
        valid = torch.ones_like(logits, dtype=torch.bool)
        if len(keys) > len(query):
            valid[:, len(query):] = patient_indices[:, None] != key_ids[None, len(query):]
        labels = torch.arange(len(query), device=query.device)
        valid[labels, labels] = True
        negatives = valid.sum(dim=1) - 1
        return nn.functional.cross_entropy(logits.masked_fill(~valid, -torch.inf), labels), negatives

    wsi_to_rna, wsi_negatives = directional(wsi, rna, memory_rna)
    rna_to_wsi, rna_negatives = directional(rna, wsi, memory_wsi)
    negatives = torch.minimum(wsi_negatives, rna_negatives)
    if (negatives < min_negatives).any():
        raise RuntimeError(
            f"programme_free InfoNCE has fewer than {min_negatives} effective negatives "
            f"at real batch size; counts={negatives.detach().cpu().tolist()}"
        )
    return 0.5 * (wsi_to_rna + rna_to_wsi), negatives


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

    .. warning::
       **This term alone does not prevent collapse -- total collapse is its
       global minimum.**  Measured on a 16x512 batch: a healthy representation
       scores 38.97, an identical-across-patients one scores 1.2e-17.  When every
       row is the same vector, ``state - mean`` is zero, so standardising yields
       ``0 / sqrt(0 + eps) = 0``, the correlation matrix is all zeros and the
       penalty vanishes.  This is the VICReg failure mode: a covariance term
       requires a **variance** term beside it or the optimiser will simply
       collapse the representation to switch it off.  It must always be paired
       with :func:`variance_floor` on the same state.  ``programme_free`` shipped
       without that pairing and collapsed to cosine 0.9999 (2026-08-02).
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
                 target_mask: torch.Tensor | None = None,
                 axis_weights: torch.Tensor | None = None) -> torch.Tensor:
    """Heteroscedastic Gaussian NLL with optional per-programme availability.

    A patient can legitimately lack one RNA-derived programme while retaining
    others.  Reducing over a row-level mask would either discard usable labels
    or silently treat unavailable targets as zeros.
    """
    if target_mask is None:
        target_mask = torch.ones_like(target, dtype=torch.bool)
    if target_mask.shape != target.shape or target.shape != mean.shape or log_variance.shape != mean.shape:
        raise ValueError("target_mask must have the same [batch, programme] shape as Gaussian NLL inputs")
    mask = target_mask.bool()
    # A padded/absent programme is not the numerical value zero.  Substitute a
    # detached prediction *before* arithmetic so 0 * NaN cannot contaminate
    # the weighted reduction used by fixed-width D2 heads.
    safe_target = torch.where(mask, target, mean.detach())
    value = 0.5 * (log_variance + (safe_target - mean).square() * torch.exp(-log_variance))
    if not mask.any():
        return value.new_zeros(())
    if axis_weights is None:
        return value.masked_select(mask).mean()
    if axis_weights.ndim != 1 or axis_weights.shape[0] != value.shape[1] or not torch.isfinite(axis_weights).all():
        raise ValueError("axis_weights must be a finite [programme] vector")
    weights = axis_weights.to(value.device, dtype=value.dtype).clamp_min(0.0)[None, :].expand_as(value)
    effective = weights * mask
    if effective.sum() <= 0:
        raise ValueError("all programme axis weights are zero; supervision would be a silent no-op")
    return (value * effective).sum() / effective.sum()
