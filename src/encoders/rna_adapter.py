"""RNA token adapter for BulkFormer and explicit gene-set features."""

from __future__ import annotations

import torch
from torch import nn

from morpheus.src.encoders.adapter_common import ModalityTokens


class RNATokenAdapter(nn.Module):
    def __init__(self, bulkformer_dim: int = 512, gene_set_dim: int = 50, hidden_dim: int = 512):
        super().__init__()
        self.bulk_proj = nn.Sequential(nn.Linear(bulkformer_dim, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(), nn.Linear(hidden_dim, hidden_dim))
        self.gene_set_proj = nn.Sequential(nn.Linear(gene_set_dim, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(), nn.Linear(hidden_dim, hidden_dim))
        self.bulk_missing = nn.Parameter(torch.zeros(1, 1, hidden_dim))
        self.gene_set_missing = nn.Parameter(torch.zeros(1, 1, hidden_dim))

    def forward(
        self,
        bulkformer: torch.Tensor | None,
        gene_sets: torch.Tensor | None = None,
        bulk_present: torch.Tensor | None = None,
        gene_set_present: torch.Tensor | None = None,
    ) -> ModalityTokens:
        if bulkformer is None and gene_sets is None:
            if bulk_present is None and gene_set_present is None:
                raise ValueError("at least one presence mask is required when RNA features are missing")
            present = bulk_present if bulk_present is not None else gene_set_present
            batch = int(present.shape[0])
            device = present.device
            tokens = torch.cat([self.bulk_missing.expand(batch, -1, -1), self.gene_set_missing.expand(batch, -1, -1)], dim=1).to(device)
            mask = torch.zeros(batch, 2, dtype=torch.bool, device=device)
            return ModalityTokens(tokens=tokens, mask=mask, modality_name="rna")
        parts, masks = [], []
        if bulkformer is not None:
            if bulk_present is None:
                bulk_present = torch.ones(bulkformer.shape[0], dtype=torch.bool, device=bulkformer.device)
            tok = self.bulk_proj(bulkformer).unsqueeze(1)
            tok = torch.where(bulk_present[:, None, None], tok, self.bulk_missing.expand(bulkformer.shape[0], -1, -1))
            parts.append(tok)
            masks.append(bulk_present[:, None])
        if gene_sets is not None:
            if gene_set_present is None:
                gene_set_present = torch.ones(gene_sets.shape[0], dtype=torch.bool, device=gene_sets.device)
            tok = self.gene_set_proj(gene_sets).unsqueeze(1)
            tok = torch.where(gene_set_present[:, None, None], tok, self.gene_set_missing.expand(gene_sets.shape[0], -1, -1))
            parts.append(tok)
            masks.append(gene_set_present[:, None])
        return ModalityTokens(tokens=torch.cat(parts, dim=1), mask=torch.cat(masks, dim=1), modality_name="rna")

