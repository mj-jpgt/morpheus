"""Reusable V2 training engine with staged objectives and resumable state."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn

from .losses import centered_cross_covariance, gaussian_nll, programme_neighbourhood_loss, supervised_programme_contrastive, symmetric_infonce
from .model import TumorStateV2


@dataclass(frozen=True)
class V2LossSchedule:
    warmup_epochs: int = 4
    identity_warmup: float = 1.0
    identity_after_warmup: float = 0.25
    programme_warmup: float = 0.50
    programme_after_warmup: float = 1.0
    neighbourhood_after_warmup: float = 0.20
    supcon_after_warmup: float = 0.20
    separation_after_warmup: float = 0.01
    rna_reconstruction_after_warmup: float = 0.03

    def weights(self, epoch: int) -> dict[str, float]:
        warmup = epoch < self.warmup_epochs
        return {
            "identity": self.identity_warmup if warmup else self.identity_after_warmup,
            "programme": self.programme_warmup if warmup else self.programme_after_warmup,
            "neighbourhood": 0.0 if warmup else self.neighbourhood_after_warmup,
            "supcon": 0.0 if warmup else self.supcon_after_warmup,
            "separation": 0.0 if warmup else self.separation_after_warmup,
            "rna_reconstruction": 0.0 if warmup else self.rna_reconstruction_after_warmup,
        }


@dataclass
class V2Trainer:
    model: TumorStateV2
    optimizer: torch.optim.Optimizer
    schedule: V2LossSchedule
    device: str
    amp_dtype: torch.dtype = torch.bfloat16

    def step(self, batch: dict[str, torch.Tensor], epoch: int) -> tuple[torch.Tensor, dict[str, float], dict[str, torch.Tensor]]:
        out_wsi = self.model(batch, view="wsi")
        out_rna = self.model(batch, view="rna")
        output = self.model(batch, view="full")
        weights = self.schedule.weights(epoch)
        loss = output["z_identity"].new_zeros(())
        metrics: dict[str, float] = {}
        if "rna" in batch:
            # Identity alignment is exclusively the WSI-only/RNA-only pair.
            identity = symmetric_infonce(out_wsi["z_identity"], out_rna["z_identity"])
            loss = loss + weights["identity"] * identity
            metrics["identity"] = float(identity.detach())
        if "programme_target" in batch:
            present = batch.get("programme_present", torch.ones(len(output["z_biology"]), device=self.device, dtype=torch.bool)).bool()
            if present.any():
                mean = output["programme_mean"][present]
                logvar = output["programme_log_variance"][present]
                target = batch["programme_target"][present]
                programme = gaussian_nll(mean, logvar, target)
                loss = loss + weights["programme"] * programme
                metrics["programme"] = float(programme.detach())
                if weights["neighbourhood"]:
                    neighbourhood = programme_neighbourhood_loss(output["z_biology"][present], target)
                    loss = loss + weights["neighbourhood"] * neighbourhood
                    metrics["neighbourhood"] = float(neighbourhood.detach())
                if weights["supcon"] and "programme_positive_mask" in batch:
                    contrastive = supervised_programme_contrastive(output["z_biology"][present], batch["programme_positive_mask"][present][:, present])
                    loss = loss + weights["supcon"] * contrastive
                    metrics["programme_supcon"] = float(contrastive.detach())
        if weights["separation"]:
            separation = centered_cross_covariance(output["z_identity"], output["z_biology"])
            loss = loss + weights["separation"] * separation
            metrics["separation"] = float(separation.detach())
        if weights["rna_reconstruction"] and "rna" in batch:
            present = batch.get("rna_present", torch.ones(len(output["z_identity"]), device=self.device, dtype=torch.bool)).bool()
            if present.any():
                reconstruction = nn.functional.mse_loss(output["rna_reconstruction"][present], batch["rna"][present])
                loss = loss + weights["rna_reconstruction"] * reconstruction
                metrics["rna_reconstruction"] = float(reconstruction.detach())
        metrics["loss"] = float(loss.detach())
        return loss, metrics, output

    def train_epoch(self, loader: Any, epoch: int) -> dict[str, float]:
        self.model.train()
        rows = []
        for batch in loader:
            batch = {key: value.to(self.device, non_blocking=True) if isinstance(value, torch.Tensor) else value for key, value in batch.items()}
            self.optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=self.amp_dtype, enabled=self.device.startswith("cuda")):
                loss, metrics, _ = self.step(batch, epoch)
            loss.backward()
            nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
            self.optimizer.step()
            rows.append(metrics)
        return {key: float(np.mean([row[key] for row in rows if key in row])) for key in {name for row in rows for name in row}}

    def save_checkpoint(self, path: str | Path, epoch: int, sampler_state: dict[str, Any], manifest: dict[str, Any], scheduler: Any | None = None) -> None:
        state = {
            "version": 2,
            "epoch": int(epoch),
            "model": self.model.state_dict(),
            "optimizer": self.optimizer.state_dict(),
            "scheduler": None if scheduler is None else scheduler.state_dict(),
            "schedule": asdict(self.schedule),
            "sampler_state": sampler_state,
            "manifest": manifest,
            "torch_rng": torch.get_rng_state(),
            "numpy_rng": np.random.get_state(),
        }
        torch.save(state, Path(path))

    def load_checkpoint(self, path: str | Path, scheduler: Any | None = None) -> dict[str, Any]:
        state = torch.load(path, map_location=self.device, weights_only=False)
        self.model.load_state_dict(state["model"])
        self.optimizer.load_state_dict(state["optimizer"])
        if scheduler is not None and state.get("scheduler") is not None:
            scheduler.load_state_dict(state["scheduler"])
        torch.set_rng_state(state["torch_rng"])
        np.random.set_state(state["numpy_rng"])
        return state
