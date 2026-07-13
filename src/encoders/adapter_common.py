"""Common adapter output types for tumor-state modeling."""

from __future__ import annotations

from dataclasses import dataclass

import torch


@dataclass
class ModalityTokens:
    """Tokenized modality evidence passed into the Query Former."""

    tokens: torch.Tensor
    mask: torch.Tensor
    modality_name: str

