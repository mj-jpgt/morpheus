"""Reusable V2 training engine with staged objectives and resumable state."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn

from .losses import (feature_decorrelation, gaussian_nll,
                     paired_infonce_with_memory,
                     programme_neighbourhood_loss,
                     supervised_programme_contrastive, symmetric_infonce,
                     variance_floor, whitened_cross_covariance)
from .model import TumorStateV2


@dataclass(frozen=True)
class V2LossSchedule:
    objective_profile: str = "full"
    warmup_epochs: int = 4
    identity_warmup: float = 1.0
    # Identity stays active; retrieval preservation is a validation gate, not
    # a fixed-epoch annealing assumption.
    identity_after_warmup: float = 1.0
    programme_warmup: float = 0.50
    programme_after_warmup: float = 1.0
    neighbourhood_after_warmup: float = 0.20
    supcon_after_warmup: float = 0.20
    separation_after_warmup: float = 0.01
    variance_after_warmup: float = 0.01
    # Off-diagonal covariance decorrelation on the biology head. The per-dim
    # variance floor cannot raise rank; this term directly counters the collapse
    # onto the low-rank programme manifold (F-R2).
    decorrelation_after_warmup: float = 0.04
    rna_reconstruction_after_warmup: float = 0.03
    fusion_identity_after_warmup: float = 0.25
    patient_consistency_after_warmup: float = 0.25
    semantic_after_warmup: float = 0.20
    fusion_identity_warmup: float = 0.10
    patient_consistency_warmup: float = 0.10
    semantic_warmup: float = 0.10

    def __post_init__(self) -> None:
        if self.objective_profile not in {"full", "identity_only", "programme_only", "programme_free"}:
            raise ValueError("objective_profile must be full, identity_only, programme_only, or programme_free")

    def weights(self, epoch: int) -> dict[str, float]:
        warmup = epoch < self.warmup_epochs
        weights = {
            "identity": self.identity_warmup if warmup else self.identity_after_warmup,
            "programme": self.programme_warmup if warmup else self.programme_after_warmup,
            "neighbourhood": 0.0 if warmup else self.neighbourhood_after_warmup,
            "supcon": 0.0 if warmup else self.supcon_after_warmup,
            "separation": 0.0 if warmup else self.separation_after_warmup,
            "variance": 0.0 if warmup else self.variance_after_warmup,
            "decorrelation": 0.0 if warmup else self.decorrelation_after_warmup,
            "rna_reconstruction": 0.0 if warmup else self.rna_reconstruction_after_warmup,
            "fusion_identity": self.fusion_identity_warmup if warmup else self.fusion_identity_after_warmup,
            "patient_consistency": self.patient_consistency_warmup if warmup else self.patient_consistency_after_warmup,
            "semantic": self.semantic_warmup if warmup else self.semantic_after_warmup,
        }
        if self.objective_profile == "identity_only":
            # Fusion consistency gives the explicit fused identity view a
            # gradient without introducing patient/biology objectives.
            return {key: (value if key in {"identity", "fusion_identity"} else 0.0) for key, value in weights.items()}
        if self.objective_profile == "programme_only":
            # decorrelation trains biology geometry, so it must stay active in
            # the profile that trains the biology head in isolation.
            return {key: (value if key in {"programme", "neighbourhood", "supcon", "decorrelation"} else 0.0)
                    for key, value in weights.items()}
        if self.objective_profile == "programme_free":
            # This is the D1 contrastive replacement arm.  It is deliberately
            # identical to programme_only with respect to decorrelation, but
            # receives no curated-programme supervision whatsoever.
            free = {key: 0.0 for key in weights}
            free["decorrelation"] = weights["decorrelation"]
            free["biology_contrastive"] = 1.0
            # `full_biology` is exported by D1, so it needs a declared
            # gradient path as well.  This is intentionally a WSI-targeted
            # consistency loss, not a hidden RNA or Hallmark target.
            free["biology_full_consistency"] = 1.0
            return free
        return weights


@dataclass
class ProgrammeMemoryBank:
    """Train-only queue that makes global programme neighbours usable with ragged batches."""

    capacity: int = 4096
    states: torch.Tensor | None = None
    indices: torch.Tensor | None = None
    cursor: int = 0
    size: int = 0

    def _ensure(self, width: int, device: torch.device) -> None:
        if self.states is None:
            self.states = torch.zeros((self.capacity, width), device=device, dtype=torch.float32)
            self.indices = torch.full((self.capacity,), -1, device=device, dtype=torch.long)

    def contrastive(self, states: torch.Tensor, neighbour_indices: torch.Tensor,
                    temperature: float = 0.10) -> tuple[torch.Tensor, int]:
        self._ensure(states.shape[1], states.device)
        if self.size == 0:
            return states.new_zeros(()), 0
        assert self.states is not None and self.indices is not None
        bank_states, bank_indices = self.states[:self.size], self.indices[:self.size]
        positive = (neighbour_indices[..., None] == bank_indices[None, None, :]).any(dim=1)
        active = positive.any(dim=1)
        if not active.any():
            return states.new_zeros(()), 0
        logits = nn.functional.normalize(states[active].float(), dim=-1) @ nn.functional.normalize(bank_states, dim=-1).T
        logits = logits / temperature
        positive = positive[active]
        return -(torch.logsumexp(logits.masked_fill(~positive, -torch.inf), dim=1) - torch.logsumexp(logits, dim=1)).mean(), int(active.sum())

    @torch.no_grad()
    def update(self, states: torch.Tensor, indices: torch.Tensor) -> None:
        self._ensure(states.shape[1], states.device)
        assert self.states is not None and self.indices is not None
        values, ids = nn.functional.normalize(states.detach().float(), dim=-1), indices.detach().long()
        for value, patient_index in zip(values, ids):
            self.states[self.cursor].copy_(value); self.indices[self.cursor] = patient_index
            self.cursor = (self.cursor + 1) % self.capacity; self.size = min(self.size + 1, self.capacity)

    def state_dict(self) -> dict[str, Any]:
        return {"capacity": self.capacity, "states": None if self.states is None else self.states.detach().cpu(),
                "indices": None if self.indices is None else self.indices.detach().cpu(), "cursor": self.cursor, "size": self.size}

    def load_state_dict(self, state: dict[str, Any], device: str) -> None:
        self.capacity, self.cursor, self.size = int(state["capacity"]), int(state["cursor"]), int(state["size"])
        self.states = None if state["states"] is None else state["states"].to(device=device, dtype=torch.float32)
        self.indices = None if state["indices"] is None else state["indices"].to(device=device, dtype=torch.long)


@dataclass
class PairedBiologyMemoryBank:
    """Detached WSI/RNA keys for D1's real-small-batch symmetric InfoNCE."""

    capacity: int = 4096
    wsi_states: torch.Tensor | None = None
    rna_states: torch.Tensor | None = None
    indices: torch.Tensor | None = None
    cursor: int = 0
    size: int = 0

    def _ensure(self, width: int, device: torch.device) -> None:
        if self.wsi_states is None:
            self.wsi_states = torch.zeros((self.capacity, width), device=device, dtype=torch.float32)
            self.rna_states = torch.zeros((self.capacity, width), device=device, dtype=torch.float32)
            self.indices = torch.full((self.capacity,), -1, device=device, dtype=torch.long)

    @torch.no_grad()
    def update(self, wsi_states: torch.Tensor, rna_states: torch.Tensor, indices: torch.Tensor) -> None:
        if wsi_states.shape != rna_states.shape or wsi_states.ndim != 2:
            raise ValueError("paired biology memory needs aligned WSI/RNA states")
        self._ensure(wsi_states.shape[1], wsi_states.device)
        assert self.wsi_states is not None and self.rna_states is not None and self.indices is not None
        wsi_values = nn.functional.normalize(wsi_states.detach().float(), dim=-1)
        rna_values = nn.functional.normalize(rna_states.detach().float(), dim=-1)
        for wsi_state, rna_state, patient_index in zip(wsi_values, rna_values, indices.detach().long()):
            self.wsi_states[self.cursor].copy_(wsi_state); self.rna_states[self.cursor].copy_(rna_state); self.indices[self.cursor] = patient_index
            self.cursor = (self.cursor + 1) % self.capacity; self.size = min(self.size + 1, self.capacity)

    def view(self) -> tuple[torch.Tensor | None, torch.Tensor | None, torch.Tensor | None]:
        if not self.size or self.wsi_states is None or self.rna_states is None or self.indices is None:
            return None, None, None
        return self.wsi_states[:self.size], self.rna_states[:self.size], self.indices[:self.size]

    @property
    def unique_count(self) -> int:
        return 0 if self.indices is None else int(torch.unique(self.indices[:self.size]).numel())

    def state_dict(self) -> dict[str, Any]:
        return {"capacity": self.capacity, "wsi_states": None if self.wsi_states is None else self.wsi_states.detach().cpu(),
                "rna_states": None if self.rna_states is None else self.rna_states.detach().cpu(),
                "indices": None if self.indices is None else self.indices.detach().cpu(), "cursor": self.cursor, "size": self.size}

    def load_state_dict(self, state: dict[str, Any], device: str) -> None:
        self.capacity, self.cursor, self.size = int(state["capacity"]), int(state["cursor"]), int(state["size"])
        self.wsi_states = None if state["wsi_states"] is None else state["wsi_states"].to(device=device, dtype=torch.float32)
        self.rna_states = None if state["rna_states"] is None else state["rna_states"].to(device=device, dtype=torch.float32)
        self.indices = None if state["indices"] is None else state["indices"].to(device=device, dtype=torch.long)


@dataclass
class V2Trainer:
    model: TumorStateV2
    optimizer: torch.optim.Optimizer
    schedule: V2LossSchedule
    device: str
    amp_dtype: torch.dtype = torch.bfloat16
    programme_memory: ProgrammeMemoryBank | None = None
    biology_memory: PairedBiologyMemoryBank | None = None
    gradient_diagnostics_every: int = 25
    decorrelation_bank_capacity: int = 512

    def __post_init__(self) -> None:
        if self.programme_memory is None:
            self.programme_memory = ProgrammeMemoryBank()
        if self.biology_memory is None:
            self.biology_memory = PairedBiologyMemoryBank()

    def liveness_parameter_groups(self) -> dict[str, list[nn.Parameter]]:
        """Named trainable groups used by the real-batch G2 telemetry.

        These are architectural groups, not merely Adam's single default
        parameter group: the latter would miss a detached biology head while
        still reporting a non-zero global norm.
        """
        groups = {
            "wsi": list(self.model.wsi.parameters()),
            "rna": list(self.model.rna.parameters()),
            "shared": [self.model.queries, *self.model.blocks.parameters(), *self.model.norm.parameters()],
            "biology_programme": [*self.model.biology.parameters(), *self.model.programme_mean.parameters(),
                                    *self.model.programme_log_variance.parameters()],
            "identity": list(self.model.identity.parameters()),
            "patient": list(self.model.patient.parameters()),
        }
        if self.schedule.objective_profile in {"programme_only", "programme_free"}:
            return {key: groups[key] for key in ("wsi", "rna", "shared", "biology_programme")}
        if self.schedule.objective_profile == "identity_only":
            return {key: groups[key] for key in ("wsi", "rna", "shared", "identity")}
        return groups

    def _gradient_group_norms(self) -> dict[str, float]:
        values: dict[str, float] = {}
        for name, parameters in self.liveness_parameter_groups().items():
            squared = sum((parameter.grad.detach().float().square().sum()
                           for parameter in parameters if parameter.grad is not None),
                          torch.zeros((), device=self.device))
            values[name] = float(torch.sqrt(squared).detach().cpu())
        return values

    @staticmethod
    def _cosine_consistency(prediction: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        return (1.0 - nn.functional.cosine_similarity(prediction, target.detach(), dim=-1)).mean()

    def _programme_loss(self, state: dict[str, torch.Tensor], batch: dict[str, torch.Tensor],
                        weights: dict[str, float], include_structure: bool) -> tuple[torch.Tensor, dict[str, float]]:
        zero = state["z_biology"].new_zeros(())
        if "programme_target" not in batch:
            return zero, {}, {}
        present = batch.get("programme_present", torch.ones(len(state["z_biology"]), device=self.device, dtype=torch.bool)).bool()
        if not present.any():
            return zero, {}, {}
        target_mask = batch.get("programme_target_mask")
        if target_mask is not None:
            target_mask = target_mask[present]
        programme = gaussian_nll(
            state["programme_mean"][present], state["programme_log_variance"][present],
            batch["programme_target"][present], target_mask,
        )
        nll_term = weights["programme"] * programme
        value, metrics = nll_term, {"programme": float(programme.detach())}
        # Keep the weighted biology sub-losses as graph-carrying tensors so the
        # gradient-conflict diagnostic can see NLL vs neighbour-KL vs supcon
        # separately — the exact pairs implicated in biology-head collapse.
        components: dict[str, torch.Tensor] = {"programme_nll": nll_term}
        structure_present = present
        if target_mask is not None:
            structure_present = present.clone()
            structure_present[present] = target_mask.all(dim=1)
        if include_structure and weights["neighbourhood"] and structure_present.any():
            neighbourhood = programme_neighbourhood_loss(state["z_biology"][structure_present], batch["programme_target"][structure_present])
            components["programme_neighbour"] = weights["neighbourhood"] * neighbourhood
            value = value + components["programme_neighbour"]
            metrics["neighbourhood"] = float(neighbourhood.detach())
        if include_structure and weights["supcon"] and "programme_positive_mask" in batch:
            positive_mask = batch["programme_positive_mask"][structure_present][:, structure_present]
            active = (positive_mask.bool() & ~torch.eye(len(positive_mask), dtype=torch.bool, device=positive_mask.device)).any(dim=1)
            terms: list[torch.Tensor] = []
            if active.any():
                terms.append(supervised_programme_contrastive(state["z_biology"][structure_present], positive_mask))
                metrics["programme_intrabatch_positive_anchors"] = float(active.sum().detach())
            if self.programme_memory is not None and "programme_neighbor_indices" in batch:
                cross_batch, n_cross_batch = self.programme_memory.contrastive(
                    state["z_biology"][structure_present], batch["programme_neighbor_indices"][structure_present]
                )
                if n_cross_batch:
                    terms.append(cross_batch)
                    metrics["programme_memory_positive_anchors"] = float(n_cross_batch)
            if terms:
                contrastive = torch.stack(terms).mean()
                components["programme_supcon"] = weights["supcon"] * contrastive
                value = value + components["programme_supcon"]
                metrics["programme_supcon"] = float(contrastive.detach())
        return value, metrics, components

    @torch.no_grad()
    def prime_biology_memory(self, loader: Any, *, minimum_unique_keys: int = 16) -> int:
        """Seed D1's detached RNA queue before the first optimisation step.

        This is train-only and is needed because the first ragged WSI batch
        may contain one patient.  Rather than quietly returning a zero InfoNCE
        loss, the subsequent loss asserts that at least eight negatives exist.
        """
        if self.biology_memory is None:
            raise RuntimeError("biology memory bank is not initialised")
        was_training = self.model.training
        self.model.eval()
        for raw_batch in loader:
            batch = {key: value.to(self.device, non_blocking=True) if isinstance(value, torch.Tensor) else value
                     for key, value in raw_batch.items()}
            out_wsi = self.model(batch, view="wsi")
            out_rna = self.model(batch, view="rna")
            self.biology_memory.update(out_wsi["z_biology"], out_rna["z_biology"], batch["indices"])
            if self.biology_memory.unique_count >= minimum_unique_keys:
                break
        self.model.train(was_training)
        if self.biology_memory.unique_count < minimum_unique_keys:
            raise RuntimeError(f"programme_free could not prime {minimum_unique_keys} distinct RNA memory keys from the real training loader")
        return self.biology_memory.size

    def step(self, batch: dict[str, torch.Tensor], epoch: int) -> tuple[torch.Tensor, dict[str, float], dict[str, torch.Tensor]]:
        out_wsi = self.model(batch, view="wsi")
        out_rna = self.model(batch, view="rna")
        output = self.model(batch, view="full")
        weights = self.schedule.weights(epoch)
        loss = output["z_identity"].new_zeros(())
        metrics: dict[str, float] = {}
        for key in ("anchor_residual_scale", "anchor_gate_mean", "anchor_correction_norm",
                    "local_slot_entropy", "patient_slot_effective_tokens", "coordinate_present_fraction"):
            if key in out_wsi:
                metrics[key] = float(out_wsi[key].detach().float().mean())
        if "rna" in batch:
            # Identity alignment is exclusively the WSI-only/RNA-only pair.
            identity = symmetric_infonce(out_wsi["z_identity"], out_rna["z_identity"])
            loss = loss + weights["identity"] * identity
            metrics["identity"] = float(identity.detach())
            if weights["fusion_identity"]:
                fusion_identity = 0.5 * (
                    self._cosine_consistency(output["z_identity"], out_wsi["z_identity"])
                    + self._cosine_consistency(output["z_identity"], out_rna["z_identity"])
                )
                loss = loss + weights["fusion_identity"] * fusion_identity
                metrics["fusion_identity"] = float(fusion_identity.detach())
            if weights["patient_consistency"]:
                target = nn.functional.normalize(out_wsi["z_identity"].detach() + out_rna["z_identity"].detach(), dim=-1)
                patient = self._cosine_consistency(output["z_patient"], target)
                loss = loss + weights["patient_consistency"] * patient
                metrics["patient_consistency"] = float(patient.detach())
        biology_contrastive_term = output["z_identity"].new_zeros(())
        biology_full_consistency_term = output["z_identity"].new_zeros(())
        if weights.get("biology_contrastive", 0.0):
            if "indices" not in batch:
                raise RuntimeError("programme_free requires canonical patient indices for its memory-bank negatives")
            assert self.biology_memory is not None
            memory_wsi, memory_rna, memory_indices = self.biology_memory.view()
            contrastive, negatives = paired_infonce_with_memory(
                out_wsi["z_biology"], out_rna["z_biology"], batch["indices"], memory_wsi, memory_rna, memory_indices,
            )
            biology_contrastive_term = weights["biology_contrastive"] * contrastive
            loss = loss + biology_contrastive_term
            metrics["biology_contrastive"] = float(contrastive.detach())
            metrics["biology_contrastive_effective_negatives_min"] = float(negatives.min().detach())
            metrics["biology_contrastive_effective_negatives_mean"] = float(negatives.float().mean().detach())
        if weights.get("biology_full_consistency", 0.0):
            # `full_biology` is an exported state under this profile.  Match it
            # to the WSI-only biology state with a stopped-gradient target: it
            # receives a direct gradient without turning the fused RNA input
            # into a trivial full→RNA objective or reintroducing Hallmarks.
            biology_full_consistency_term = self._cosine_consistency(output["z_biology"], out_wsi["z_biology"])
            biology_full_consistency_term = weights["biology_full_consistency"] * biology_full_consistency_term
            loss = loss + biology_full_consistency_term
            metrics["biology_full_consistency"] = float(biology_full_consistency_term.detach())
        # Every exported biology state must receive direct supervision.  WSI
        # receives the structural programme losses because it is the primary
        # molecular-prompting view.
        programme_total = output["z_biology"].new_zeros(())
        programme_component_totals: dict[str, torch.Tensor] = {}
        structure_for_all_biology_views = self.schedule.objective_profile == "programme_only"
        views = () if self.schedule.objective_profile == "programme_free" else (
            ("wsi", out_wsi, True, 1.0),
            ("rna", out_rna, structure_for_all_biology_views, 0.5),
            ("full", output, structure_for_all_biology_views, 1.0),
        )
        for name, state, structure, scale in views:
            programme_loss, programme_metrics, programme_components = self._programme_loss(state, batch, weights, structure)
            loss = loss + scale * programme_loss
            programme_total = programme_total + scale * programme_loss
            for component_name, component_value in programme_components.items():
                programme_component_totals[component_name] = (
                    programme_component_totals.get(component_name, 0.0) + scale * component_value
                )
            metrics.update({f"{name}_{key}": value for key, value in programme_metrics.items()})
        # Stable aggregate retained for legacy monitors; the namespaced terms
        # above identify which view contributed to it.
        metrics["programme"] = float(programme_total.detach())
        separation_term = output["z_identity"].new_zeros(())
        if weights["separation"]:
            separation = whitened_cross_covariance(output["z_identity"], output["z_biology"])
            separation_term = weights["separation"] * separation
            loss = loss + separation_term
            metrics["separation"] = float(separation.detach())
        variance_term = output["z_identity"].new_zeros(())
        if weights["variance"]:
            variance = variance_floor(output["z_identity"]) + variance_floor(output["z_biology"])
            variance_term = weights["variance"] * variance
            loss = loss + variance_term
            metrics["variance_floor"] = float(variance.detach())
        decorrelation_term = output["z_identity"].new_zeros(())
        if weights["decorrelation"]:
            # Decorrelate the WSI biology view (the view the collapse pressure and
            # the reported rank fingerprint sit on). Real uncapped-patch batches
            # hold only B~1-3 patients, far below a usable 256-D correlation
            # estimate, so we pool the current (gradient-carrying) rows with a
            # detached ring buffer of recent biology features. Gradient flows only
            # through the current rows; the bank supplies sample count. The bank is
            # updated in training only, so evaluation never contaminates it.
            current = out_wsi["z_biology"]
            bank = getattr(self, "_decorr_bank", None)
            pool = current if bank is None else torch.cat([current, bank.to(current.device)], dim=0)
            decorrelation = feature_decorrelation(pool)
            decorrelation_term = weights["decorrelation"] * decorrelation
            loss = loss + decorrelation_term
            metrics["decorrelation"] = float(decorrelation.detach())
            metrics["decorrelation_pool"] = float(pool.shape[0])
            if torch.is_grad_enabled():
                with torch.no_grad():
                    fresh = current.detach().float()
                    combined = fresh if bank is None else torch.cat([bank, fresh], dim=0)
                    self._decorr_bank = combined[-self.decorrelation_bank_capacity:]
        for name, value in (("identity", output["z_identity"]), ("biology", output["z_biology"])):
            metrics[f"{name}_feature_std"] = float(value.detach().float().std(dim=0).mean())
        reconstruction_term = output["z_identity"].new_zeros(())
        if weights["rna_reconstruction"] and "rna" in batch:
            present = batch.get("rna_present", torch.ones(len(output["z_identity"]), device=self.device, dtype=torch.bool)).bool()
            if present.any():
                reconstruction = nn.functional.mse_loss(out_wsi["rna_reconstruction"][present], batch["rna"][present])
                reconstruction_term = weights["rna_reconstruction"] * reconstruction
                loss = loss + reconstruction_term
                metrics["rna_reconstruction"] = float(reconstruction.detach())
        if weights["semantic"] and "semantic_target" in batch and "z_semantic" in out_wsi:
            present = batch.get("semantic_present", torch.ones(len(out_wsi["z_semantic"]), device=self.device, dtype=torch.bool)).bool()
            if present.any():
                semantic = self._cosine_consistency(out_wsi["z_semantic"][present], batch["semantic_target"][present])
                loss = loss + weights["semantic"] * semantic
                metrics["semantic"] = float(semantic.detach())
        metrics["loss"] = float(loss.detach())
        # Expose the biology sub-losses (NLL / neighbour-KL / supcon) as separate
        # components so the gradient-conflict diagnostic can detect intra-biology
        # conflict, not just biology-vs-identity. Falls back to the combined
        # programme term if a profile produced no split components.
        self._last_loss_components = {
            "identity": weights["identity"] * identity if "rna" in batch else loss.new_zeros(()),
            "reconstruction": reconstruction_term,
            "separation": separation_term,
            "decorrelation": decorrelation_term,
            "biology_contrastive": biology_contrastive_term,
            "biology_full_consistency": biology_full_consistency_term,
            "variance": variance_term,
            **(programme_component_totals or {"programme": programme_total}),
        }
        if weights.get("biology_contrastive", 0.0) and torch.is_grad_enabled() and "indices" in batch:
            assert self.biology_memory is not None
            self.biology_memory.update(out_wsi["z_biology"], out_rna["z_biology"], batch["indices"])
        return loss, metrics, output

    def _gradient_conflict_metrics(self) -> dict[str, float]:
        """Log objective-gradient cosines; do not silently alter optimisation.

        This is sampled periodically because it takes extra backwards through
        the graph.  A sustained negative cosine is evidence required before a
        future PCGrad/GradNorm experiment, not a reason to enable either by
        default.
        """
        components = getattr(self, "_last_loss_components", {})
        # This diagnostic retains the current forward graph while it obtains a
        # gradient for every declared objective.  On a large ragged WSI batch
        # that can require several hundred additional MiB.  Telemetry must
        # never turn an otherwise valid experiment into an OOM failure.  We
        # therefore record an explicit (rather than silent) skipped sample if
        # there is insufficient headroom; a dedicated lower-token calibration
        # pass supplies the actual gradient-cosine measurements.
        if self.device.startswith("cuda"):
            free_bytes, total_bytes = torch.cuda.mem_get_info(self.device)
            minimum_headroom = 768 * 1024 * 1024
            if free_bytes < minimum_headroom:
                return {
                    "gradient_diagnostics_skipped_low_memory": 1.0,
                    "gradient_diagnostics_free_gib": float(free_bytes / 1024**3),
                    "gradient_diagnostics_total_gib": float(total_bytes / 1024**3),
                }
        # Measure conflict at the SHARED TRUNK where the 32 query slots are jointly
        # produced and then split into the identity/biology heads — i.e. the query
        # bank + Query-Former blocks + final norm. The previous set used pre-trunk,
        # per-view input encoders (wsi.patch / rna.projection), which sit *before*
        # the head split and so confound the cosine and miss the trunk entirely.
        shared = [self.model.queries]
        shared += [p for p in self.model.blocks.parameters() if p.requires_grad]
        shared += [p for p in self.model.norm.parameters() if p.requires_grad]
        gradients: dict[str, list[torch.Tensor]] = {}
        try:
            for name, component in components.items():
                if not component.requires_grad or not torch.isfinite(component).all() or component.detach().abs().item() == 0:
                    continue
                values = torch.autograd.grad(component, shared, retain_graph=True, allow_unused=True)
                # Preserve parameter alignment across objectives; dropping an
                # unused gradient would make later zip() comparisons invalid.
                gradients[name] = [
                    torch.zeros_like(parameter, dtype=torch.float32) if value is None else value.detach().float()
                    for parameter, value in zip(shared, values)
                ]
        except torch.OutOfMemoryError:
            # Fragmentation can still defeat the preflight estimate.  Release
            # cached temporary blocks and make the skipped measurement visible
            # in the epoch log while allowing optimisation to continue.
            if self.device.startswith("cuda"):
                torch.cuda.empty_cache()
            return {"gradient_diagnostics_skipped_oom": 1.0}
        result: dict[str, float] = {}
        cosines: list[float] = []
        names = sorted(gradients)
        for index, left_name in enumerate(names):
            for right_name in names[index + 1:]:
                left, right = gradients[left_name], gradients[right_name]
                if not left or not right:
                    continue
                dot = sum((a * b).sum() for a, b in zip(left, right))
                norm_left = torch.sqrt(sum(a.square().sum() for a in left)).clamp_min(1e-12)
                norm_right = torch.sqrt(sum(b.square().sum() for b in right)).clamp_min(1e-12)
                cosine = float((dot / (norm_left * norm_right)).cpu())
                result[f"gradient_cosine_{left_name}_{right_name}"] = cosine
                cosines.append(cosine)
        if cosines:
            result["gradient_cosine_median"] = float(np.median(cosines))
        return result

    def train_epoch(self, loader: Any, epoch: int) -> dict[str, float]:
        self.model.train()
        rows = []
        first_gradients: dict[str, float] | None = None
        last_gradients: dict[str, float] | None = None
        for step, batch in enumerate(loader):
            batch = {key: value.to(self.device, non_blocking=True) if isinstance(value, torch.Tensor) else value for key, value in batch.items()}
            self.optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=self.amp_dtype, enabled=self.device.startswith("cuda")):
                loss, metrics, output = self.step(batch, epoch)
            if self.gradient_diagnostics_every > 0 and step % self.gradient_diagnostics_every == 0:
                metrics.update(self._gradient_conflict_metrics())
            loss.backward()
            group_norms = self._gradient_group_norms()
            if first_gradients is None:
                first_gradients = group_norms
            last_gradients = group_norms
            nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
            self.optimizer.step()
            if self.programme_memory is not None and "indices" in batch:
                self.programme_memory.update(output["z_biology"], batch["indices"])
            rows.append(metrics)
        result = {key: float(np.mean([row[key] for row in rows if key in row])) for key in {name for row in rows for name in row}}
        for name, value in (first_gradients or {}).items():
            result[f"gradient_norm_{name}_first"] = value
        for name, value in (last_gradients or {}).items():
            result[f"gradient_norm_{name}_last"] = value
        return result

    @torch.no_grad()
    def evaluate_epoch(self, loader: Any, epoch: int) -> dict[str, float]:
        """Evaluate the same declared objectives without dropout or updates."""
        self.model.eval()
        rows = []
        for batch in loader:
            batch = {key: value.to(self.device, non_blocking=True) if isinstance(value, torch.Tensor) else value for key, value in batch.items()}
            with torch.autocast(device_type="cuda", dtype=self.amp_dtype, enabled=self.device.startswith("cuda")):
                _, metrics, _ = self.step(batch, epoch)
            rows.append(metrics)
        return {key: float(np.mean([row[key] for row in rows if key in row])) for key in {name for row in rows for name in row}}

    def save_checkpoint(self, path: str | Path, epoch: int, sampler_state: dict[str, Any], manifest: dict[str, Any], scheduler: Any | None = None) -> None:
        state = {
            "version": 3,
            "epoch": int(epoch),
            "model": self.model.state_dict(),
            "optimizer": self.optimizer.state_dict(),
            "scheduler": None if scheduler is None else scheduler.state_dict(),
            "schedule": asdict(self.schedule),
            "sampler_state": sampler_state,
            "manifest": manifest,
            "torch_rng": torch.get_rng_state(),
            "numpy_rng": np.random.get_state(),
            "programme_memory": None if self.programme_memory is None else self.programme_memory.state_dict(),
            "biology_memory": None if self.biology_memory is None else self.biology_memory.state_dict(),
        }
        torch.save(state, Path(path))

    def load_checkpoint(self, path: str | Path, scheduler: Any | None = None) -> dict[str, Any]:
        state = torch.load(path, map_location=self.device, weights_only=False)
        self.model.load_state_dict(state["model"])
        self.optimizer.load_state_dict(state["optimizer"])
        if scheduler is not None and state.get("scheduler") is not None:
            scheduler.load_state_dict(state["scheduler"])
        # `map_location=self.device` can restore this CPU generator state as
        # a CUDA tensor (or an array in older checkpoints).  PyTorch's global
        # CPU generator accepts only a CPU uint8 tensor, so normalise it
        # explicitly before resuming an OOM-retried job.
        torch_rng = torch.as_tensor(state["torch_rng"], dtype=torch.uint8, device="cpu").contiguous()
        torch.set_rng_state(torch_rng)
        np.random.set_state(state["numpy_rng"])
        if self.programme_memory is not None and state.get("programme_memory") is not None:
            self.programme_memory.load_state_dict(state["programme_memory"], self.device)
        if self.biology_memory is not None and state.get("biology_memory") is not None:
            self.biology_memory.load_state_dict(state["biology_memory"], self.device)
        return state
