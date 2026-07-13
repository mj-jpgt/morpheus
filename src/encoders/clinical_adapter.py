"""Clinical tabular token adapter."""

from __future__ import annotations

import torch
from torch import nn

from morpheus.src.encoders.adapter_common import ModalityTokens


class ClinicalTokenAdapter(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int = 512):
        super().__init__()
        self.proj = nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.LayerNorm(hidden_dim), nn.GELU(), nn.Linear(hidden_dim, hidden_dim))
        self.missing = nn.Parameter(torch.zeros(1, 1, hidden_dim))

    def forward(self, x: torch.Tensor | None, present: torch.Tensor | None = None) -> ModalityTokens:
        if x is None:
            if present is None:
                raise ValueError("present mask is required when clinical features are missing")
            return ModalityTokens(tokens=self.missing.expand(present.shape[0], -1, -1), mask=torch.zeros(present.shape[0], 1, dtype=torch.bool, device=present.device), modality_name="clinical")
        if present is None:
            present = torch.ones(x.shape[0], dtype=torch.bool, device=x.device)
        tokens = self.proj(x).unsqueeze(1)
        tokens = torch.where(present[:, None, None], tokens, self.missing.expand(x.shape[0], -1, -1))
        return ModalityTokens(tokens=tokens, mask=present[:, None], modality_name="clinical")

