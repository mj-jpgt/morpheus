"""Strict train-cancer-only self-supervision for the V2.1 slide aggregator.

The module consumes already-frozen H-Optimus patch embeddings; it never
re-extracts tiles and deliberately has no knowledge of RNA, cancer labels, or
outer-test bags.  The controller decides which patients are supplied.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import nn

from .model import HierarchicalWSIEncoder


@dataclass(frozen=True)
class SlidePretrainingConfig:
    mask_fraction: float = 0.30
    view_keep_fraction: float = 0.70
    target_dim: int = 128
    masked_bag_weight: float = 1.0
    consistency_weight: float = 1.0
    seed: int = 917


class SlidePretrainingObjective(nn.Module):
    """Masked-bag feature reconstruction plus two-view slide consistency.

    A randomly projected mean of held-out patch features is reconstructed from
    a visible-patch slide encoding.  This is intentionally called *masked bag*
    reconstruction, not token prediction: the current aggregator has no
    patch-to-patch decoder and should not claim a denser objective than it
    implements.
    """

    def __init__(self, patch_dim: int, hidden_dim: int, config: SlidePretrainingConfig | None = None) -> None:
        super().__init__()
        self.config = config or SlidePretrainingConfig()
        if not 0.0 < self.config.mask_fraction < 1.0 or not 0.0 < self.config.view_keep_fraction <= 1.0:
            raise ValueError("mask_fraction and view_keep_fraction must be in (0, 1]")
        generator = torch.Generator(device="cpu").manual_seed(self.config.seed)
        projection = torch.randn((patch_dim, self.config.target_dim), generator=generator)
        projection = nn.functional.normalize(projection, dim=0)
        self.register_buffer("target_projection", projection, persistent=True)
        self.masked_bag_decoder = nn.Sequential(
            nn.LayerNorm(hidden_dim), nn.Linear(hidden_dim, hidden_dim), nn.GELU(),
            nn.Linear(hidden_dim, self.config.target_dim),
        )

    @staticmethod
    def _masked_mean(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        weight = mask.to(values.dtype).unsqueeze(-1)
        return (values * weight).sum(1) / weight.sum(1).clamp_min(1.0)

    def loss(self, encoder: HierarchicalWSIEncoder, batch: dict[str, torch.Tensor]) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        patches, valid = batch["patches"], batch["patch_mask"].bool()
        generator = torch.Generator(device=patches.device).manual_seed(self.config.seed + int(valid.sum().item()))
        random = torch.rand(valid.shape, device=patches.device, generator=generator)
        masked = valid & (random < self.config.mask_fraction)
        # Preserve one visible tile and one target tile for short slides.
        visible = valid & ~masked
        for row in range(len(valid)):
            active = torch.where(valid[row])[0]
            if len(active) < 2:
                masked[row].zero_(); visible[row] = valid[row]
            elif not visible[row].any():
                visible[row, active[0]] = True; masked[row, active[0]] = False
            elif not masked[row].any():
                masked[row, active[-1]] = True; visible[row, active[-1]] = False
        coordinate_present = batch.get("coordinate_present")
        first, _ = encoder(patches, visible, batch["slide_ids"], batch.get("coordinates"), coordinate_present)
        raw_target = self._masked_mean(patches, masked) @ self.target_projection.to(patches.dtype)
        prediction = self.masked_bag_decoder(first.tokens.mean(1))
        masked_bag = nn.functional.mse_loss(prediction, raw_target.detach())

        view_a = valid & (torch.rand(valid.shape, device=patches.device, generator=generator) < self.config.view_keep_fraction)
        view_b = valid & (torch.rand(valid.shape, device=patches.device, generator=generator) < self.config.view_keep_fraction)
        for row in range(len(valid)):
            for view in (view_a, view_b):
                if not view[row].any():
                    view[row, torch.where(valid[row])[0][0]] = True
        out_a, _ = encoder(patches, view_a, batch["slide_ids"], batch.get("coordinates"), coordinate_present)
        out_b, _ = encoder(patches, view_b, batch["slide_ids"], batch.get("coordinates"), coordinate_present)
        consistency = 1.0 - nn.functional.cosine_similarity(out_a.tokens.mean(1), out_b.tokens.mean(1), dim=-1).mean()
        total = self.config.masked_bag_weight * masked_bag + self.config.consistency_weight * consistency
        return total, {
            "pretrain_masked_bag": masked_bag.detach(),
            "pretrain_view_consistency": consistency.detach(),
            "pretrain_masked_fraction": masked.float().sum().detach() / valid.float().sum().clamp_min(1.0),
        }


class SlidePretrainer:
    """Small training wrapper; callers must pass development-only loaders."""

    def __init__(self, encoder: HierarchicalWSIEncoder, objective: SlidePretrainingObjective,
                 optimizer: torch.optim.Optimizer, device: str = "cuda") -> None:
        self.encoder, self.objective, self.optimizer, self.device = encoder, objective, optimizer, device

    def train_epoch(self, loader) -> dict[str, float]:
        self.encoder.train(); self.objective.train()
        rows: list[dict[str, float]] = []
        for batch in loader:
            batch = {key: value.to(self.device, non_blocking=True) if isinstance(value, torch.Tensor) else value
                     for key, value in batch.items()}
            self.optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=self.device.startswith("cuda")):
                loss, metrics = self.objective.loss(self.encoder, batch)
            loss.backward(); nn.utils.clip_grad_norm_(self.encoder.parameters(), 1.0); self.optimizer.step()
            rows.append({key: float(value.cpu()) for key, value in metrics.items()} | {"pretrain_loss": float(loss.detach().cpu())})
        return {key: sum(row[key] for row in rows if key in row) / sum(key in row for row in rows)
                for key in {name for row in rows for name in row}}
