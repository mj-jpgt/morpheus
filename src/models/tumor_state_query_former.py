"""Tumor-State Query Former over frozen modality tokens."""

from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import nn

from morpheus.src.encoders.adapter_common import ModalityTokens


@dataclass(frozen=True)
class QueryFormerConfig:
    hidden_dim: int = 512
    num_layers: int = 2
    num_heads: int = 8
    dropout: float = 0.1
    shared_slots: int = 16
    wsi_residual_slots: int = 4
    rna_residual_slots: int = 4
    clinical_residual_slots: int = 4
    genomic_residual_slots: int = 4
    uncertainty_slots: int = 4
    task_slots: int = 4

    @property
    def total_slots(self) -> int:
        return (
            self.shared_slots
            + self.wsi_residual_slots
            + self.rna_residual_slots
            + self.clinical_residual_slots
            + self.genomic_residual_slots
            + self.uncertainty_slots
            + self.task_slots
        )


class QueryBlock(nn.Module):
    def __init__(self, hidden_dim: int, num_heads: int, dropout: float):
        super().__init__()
        self.cross_attn = nn.MultiheadAttention(hidden_dim, num_heads, dropout=dropout, batch_first=True)
        self.self_attn = nn.MultiheadAttention(hidden_dim, num_heads, dropout=dropout, batch_first=True)
        self.cross_norm = nn.LayerNorm(hidden_dim)
        self.self_norm = nn.LayerNorm(hidden_dim)
        self.ffn_norm = nn.LayerNorm(hidden_dim)
        self.ffn = nn.Sequential(nn.Linear(hidden_dim, hidden_dim * 4), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden_dim * 4, hidden_dim))
        self.dropout = nn.Dropout(dropout)

    def forward(self, queries: torch.Tensor, tokens: torch.Tensor, key_padding_mask: torch.Tensor | None):
        cross, cross_weights = self.cross_attn(queries, tokens, tokens, key_padding_mask=key_padding_mask, need_weights=False)
        queries = self.cross_norm(queries + self.dropout(cross))
        self_out, self_weights = self.self_attn(queries, queries, queries, need_weights=False)
        queries = self.self_norm(queries + self.dropout(self_out))
        queries = self.ffn_norm(queries + self.dropout(self.ffn(queries)))
        return queries, {"cross": cross_weights, "self": self_weights}


class TumorStateQueryFormer(nn.Module):
    def __init__(self, config: QueryFormerConfig | None = None):
        super().__init__()
        self.config = config or QueryFormerConfig()
        self.query = nn.Parameter(torch.randn(1, self.config.total_slots, self.config.hidden_dim) * 0.02)
        self.blocks = nn.ModuleList([QueryBlock(self.config.hidden_dim, self.config.num_heads, self.config.dropout) for _ in range(self.config.num_layers)])
        self.output_norm = nn.LayerNorm(self.config.hidden_dim)

    def forward(self, modalities: dict[str, ModalityTokens]) -> dict[str, torch.Tensor | list[dict[str, torch.Tensor]] | dict[str, torch.Tensor]]:
        if not modalities:
            raise ValueError("At least one modality token group is required")
        tokens = torch.cat([value.tokens for value in modalities.values()], dim=1)
        masks = torch.cat([value.mask for value in modalities.values()], dim=1).bool()
        key_padding_mask = ~masks
        batch = tokens.shape[0]
        queries = self.query.expand(batch, -1, -1)
        attention_maps = []
        for block in self.blocks:
            queries, attn = block(queries, tokens, key_padding_mask)
            attention_maps.append(attn)
        queries = self.output_norm(queries)
        slots = self._split_slots(queries)
        slots["attention_maps"] = attention_maps
        slots["modality_masks"] = {name: value.mask for name, value in modalities.items()}
        slots["z_patient"] = slots["z_shared"].mean(dim=1)
        return slots

    def _split_slots(self, queries: torch.Tensor) -> dict[str, torch.Tensor]:
        cfg = self.config
        cursor = 0

        def take(n: int) -> torch.Tensor:
            nonlocal cursor
            value = queries[:, cursor : cursor + n]
            cursor += n
            return value

        return {
            "z_shared": take(cfg.shared_slots),
            "z_wsi_resid": take(cfg.wsi_residual_slots),
            "z_rna_resid": take(cfg.rna_residual_slots),
            "z_clinical_resid": take(cfg.clinical_residual_slots),
            "z_genomic_resid": take(cfg.genomic_residual_slots),
            "z_uncertainty": take(cfg.uncertainty_slots),
            "z_task": take(cfg.task_slots),
        }
