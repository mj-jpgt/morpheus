"""WSI token adapter for frozen H-Optimus-0 patient embeddings."""

from __future__ import annotations

import torch
from torch import nn

from morpheus.src.encoders.adapter_common import ModalityTokens


class WSITokenAdapter(nn.Module):
    def __init__(self, input_dim: int = 1536, hidden_dim: int = 512):
        super().__init__()
        self.proj = nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(), nn.Linear(hidden_dim, hidden_dim))
        self.missing = nn.Parameter(torch.zeros(1, 1, hidden_dim))

    def forward(self, x: torch.Tensor | None, present: torch.Tensor | None = None) -> ModalityTokens:
        if x is None:
            if present is None:
                raise ValueError("present mask is required when WSI features are missing")
            tokens = self.missing.expand(present.shape[0], -1, -1)
            mask = torch.zeros(present.shape[0], 1, dtype=torch.bool, device=present.device)
            return ModalityTokens(tokens=tokens, mask=mask, modality_name="wsi")
        tokens = self.proj(x).unsqueeze(1)
        if present is None:
            present = torch.ones(x.shape[0], dtype=torch.bool, device=x.device)
        tokens = torch.where(present[:, None, None], tokens, self.missing.expand(x.shape[0], -1, -1))
        return ModalityTokens(tokens=tokens, mask=present[:, None], modality_name="wsi")


class WSIPatchTokenAdapter(nn.Module):
    """Adapter for variable-length precomputed WSI patch bags."""

    def __init__(self, input_dim: int = 2048, hidden_dim: int = 512, use_coords: bool = True, num_slots: int = 0):
        super().__init__()
        self.use_coords = use_coords
        self.patch_proj = nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(), nn.Linear(hidden_dim, hidden_dim))
        self.coord_proj = nn.Sequential(nn.Linear(2, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(), nn.Linear(hidden_dim, hidden_dim))
        self.slot_aggregator = WSISlotAggregator(hidden_dim, num_slots) if num_slots > 0 else None
        self.missing = nn.Parameter(torch.zeros(1, 1, hidden_dim))

    def forward(self, feats: torch.Tensor | None, mask: torch.Tensor | None = None, coords: torch.Tensor | None = None) -> ModalityTokens:
        if feats is None:
            if mask is None:
                raise ValueError("mask is required when WSI patch features are missing")
            tokens = self.missing.expand(mask.shape[0], 1, -1)
            return ModalityTokens(tokens=tokens, mask=torch.zeros(mask.shape[0], 1, dtype=torch.bool, device=mask.device), modality_name="wsi_patch")
        if feats.ndim != 3:
            raise ValueError(f"Patch features must be [batch, tokens, dim], got {tuple(feats.shape)}")
        if feats.shape[-1] != self.patch_proj[0].in_features:
            raise ValueError(
                f"Expected {self.patch_proj[0].in_features}-D WSI patch features, got {feats.shape[-1]}"
            )
        if mask is None:
            mask = torch.ones(feats.shape[:2], dtype=torch.bool, device=feats.device)
        tokens = self.patch_proj(feats)
        if self.use_coords and coords is not None:
            coord_scale = coords.abs().amax(dim=1, keepdim=True).clamp_min(1.0)
            tokens = tokens + self.coord_proj(coords / coord_scale)
        tokens = torch.where(mask[:, :, None], tokens, self.missing.expand(feats.shape[0], feats.shape[1], -1))
        if self.slot_aggregator is not None:
            tokens, slot_mask = self.slot_aggregator(tokens, mask)
            return ModalityTokens(tokens=tokens, mask=slot_mask, modality_name="wsi_patch")
        return ModalityTokens(tokens=tokens, mask=mask, modality_name="wsi_patch")


class WSISlotAggregator(nn.Module):
    """Normalized softmax pooling of variable-length patch bags into WSI slots.

    Slots are deliberately query-only (rather than learned pooling logits) so the
    returned attention is interpretable as cosine similarity to a normalized tissue
    token. ``last_attention`` is detached and intended for batch-level QC logging.
    """

    def __init__(self, hidden_dim: int, num_slots: int = 32):
        super().__init__()
        if num_slots < 1:
            raise ValueError("num_slots must be positive")
        self.queries = nn.Parameter(torch.empty(num_slots, hidden_dim))
        nn.init.orthogonal_(self.queries)
        self.token_projection = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.slot_projection = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.value_projection = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.log_temperature = nn.Parameter(torch.tensor(0.0))
        self.last_attention: torch.Tensor | None = None
        self.last_effective_token_count: torch.Tensor | None = None

    def forward(self, tokens: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if tokens.ndim != 3 or mask.shape != tokens.shape[:2]:
            raise ValueError("tokens must be [batch, tokens, hidden] and mask [batch, tokens]")
        token_keys = torch.nn.functional.normalize(self.token_projection(tokens), dim=-1)
        slot_queries = torch.nn.functional.normalize(self.slot_projection(self.queries), dim=-1)
        logits = torch.einsum("bth,sh->bst", token_keys, slot_queries)
        logits = logits * self.log_temperature.exp().clamp(0.01, 100.0)
        has_tokens = mask.any(dim=1)
        # Softmax cannot consume an all-masked row.  Give such rows a temporary
        # zero-valued token then retain a false output mask for downstream models.
        safe_mask = mask.clone()
        safe_mask[~has_tokens, 0] = True
        logits = logits.masked_fill(~safe_mask[:, None, :], torch.finfo(logits.dtype).min)
        attention = torch.softmax(logits, dim=-1)
        values = torch.nn.functional.normalize(self.value_projection(tokens), dim=-1)
        pooled = torch.einsum("bst,bth->bsh", attention, values)
        pooled = torch.where(has_tokens[:, None, None], pooled, torch.zeros_like(pooled))
        self.last_attention = attention.detach()
        self.last_effective_token_count = (1.0 / attention.square().sum(dim=-1).clamp_min(1e-12)).detach()
        return pooled, has_tokens[:, None].expand(-1, self.queries.shape[0])
