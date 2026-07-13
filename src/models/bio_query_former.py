"""Typed biological Query Former for MORPHEUS V2."""

from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import nn

from morpheus.src.encoders.adapter_common import ModalityTokens
from morpheus.src.models.tumor_state_query_former import QueryBlock


@dataclass(frozen=True)
class BioQueryFormerConfig:
    hidden_dim: int = 512
    num_layers: int = 2
    num_heads: int = 8
    dropout: float = 0.1
    identity_slots: int = 8
    biology_slots: int = 8
    program_slots: int = 16
    wsi_residual_slots: int = 4
    rna_residual_slots: int = 4
    clinical_residual_slots: int = 4
    uncertainty_slots: int = 4
    hypothesis_slots: int = 4

    @property
    def total_slots(self) -> int:
        return (
            self.identity_slots
            + self.biology_slots
            + self.program_slots
            + self.wsi_residual_slots
            + self.rna_residual_slots
            + self.clinical_residual_slots
            + self.uncertainty_slots
            + self.hypothesis_slots
        )


class BioQueryFormer(nn.Module):
    """Query bottleneck with typed output surfaces for V2 tasks."""

    def __init__(self, config: BioQueryFormerConfig | None = None):
        super().__init__()
        self.config = config or BioQueryFormerConfig()
        # Typed queries must start distinct: otherwise the self-attention stack can
        # make identity and biology surfaces exchangeable from the first update.
        self.query = nn.Parameter(torch.empty(1, self.config.total_slots, self.config.hidden_dim))
        nn.init.orthogonal_(self.query[0])
        self.blocks = nn.ModuleList([QueryBlock(self.config.hidden_dim, self.config.num_heads, self.config.dropout) for _ in range(self.config.num_layers)])
        self.output_norm = nn.LayerNorm(self.config.hidden_dim)

    def forward(self, modalities: dict[str, ModalityTokens]) -> dict[str, torch.Tensor | dict[str, torch.Tensor] | list[dict[str, torch.Tensor]]]:
        if not modalities:
            raise ValueError("At least one modality token group is required")
        tokens = torch.cat([value.tokens for value in modalities.values()], dim=1)
        masks = torch.cat([value.mask for value in modalities.values()], dim=1).bool()
        key_padding_mask = ~masks
        queries = self.query.expand(tokens.shape[0], -1, -1)
        attention_maps = []
        for block in self.blocks:
            queries, attn = block(queries, tokens, key_padding_mask)
            attention_maps.append(attn)
        slots = self._split_slots(self.output_norm(queries))
        slots["z_identity"] = slots["identity_slots"].mean(dim=1)
        slots["z_biology"] = slots["biology_slots"].mean(dim=1)
        slots["z_programs"] = slots["program_slots"]
        slots["z_wsi_residual"] = slots["wsi_residual_slots"].mean(dim=1)
        slots["z_rna_residual"] = slots["rna_residual_slots"].mean(dim=1)
        slots["z_clinical_residual"] = slots["clinical_residual_slots"].mean(dim=1)
        slots["z_uncertainty"] = slots["uncertainty_slots"].mean(dim=1)
        slots["z_hypothesis"] = slots["hypothesis_slots"].mean(dim=1)
        slots["attention_maps"] = attention_maps
        slots["modality_masks"] = {name: value.mask for name, value in modalities.items()}
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
            "identity_slots": take(cfg.identity_slots),
            "biology_slots": take(cfg.biology_slots),
            "program_slots": take(cfg.program_slots),
            "wsi_residual_slots": take(cfg.wsi_residual_slots),
            "rna_residual_slots": take(cfg.rna_residual_slots),
            "clinical_residual_slots": take(cfg.clinical_residual_slots),
            "uncertainty_slots": take(cfg.uncertainty_slots),
            "hypothesis_slots": take(cfg.hypothesis_slots),
        }
