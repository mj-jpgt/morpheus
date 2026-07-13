"""Hierarchical typed Tumor-State Query Former for local V2 training."""

from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import nn

from morpheus.src.encoders.adapter_common import ModalityTokens


@dataclass(frozen=True)
class V2ModelConfig:
    patch_dim: int = 1536
    rna_dim: int = 512
    hidden_dim: int = 512
    heads: int = 8
    layers: int = 4
    local_slots: int = 8
    slide_slots: int = 8
    patient_slots: int = 16
    identity_slots: int = 4
    biology_slots: int = 4
    context_slots: int = 2
    residual_slots: int = 4
    uncertainty_slots: int = 2
    dropout: float = 0.10

    @property
    def query_slots(self) -> int:
        return self.identity_slots + self.biology_slots + self.context_slots + 5 * self.residual_slots + self.uncertainty_slots


class SoftSlotPool(nn.Module):
    """Orthogonal normalized slot pooling with softmax over valid tokens."""

    def __init__(self, hidden_dim: int, slots: int) -> None:
        super().__init__()
        self.queries = nn.Parameter(torch.empty(slots, hidden_dim))
        nn.init.orthogonal_(self.queries)
        self.token_key = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.slot_key = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.value = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.log_temperature = nn.Parameter(torch.tensor(0.0))

    def forward(self, tokens: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        if tokens.ndim != 3 or mask.shape != tokens.shape[:2]:
            raise ValueError("tokens must be [batch, token, hidden] and mask [batch, token]")
        available = mask.any(dim=1)
        safe_mask = mask.clone()
        safe_mask[~available, 0] = True
        keys = nn.functional.normalize(self.token_key(tokens), dim=-1)
        queries = nn.functional.normalize(self.slot_key(self.queries), dim=-1)
        logits = torch.einsum("bth,sh->bst", keys, queries) * self.log_temperature.exp().clamp(0.01, 100.0)
        logits = logits.masked_fill(~safe_mask[:, None, :], -1e4)
        attention = torch.softmax(logits, dim=-1)
        pooled = torch.einsum("bst,bth->bsh", attention, nn.functional.normalize(self.value(tokens), dim=-1))
        pooled = torch.where(available[:, None, None], pooled, torch.zeros_like(pooled))
        return pooled, available[:, None].expand(-1, self.queries.shape[0]), attention


class QueryBlock(nn.Module):
    def __init__(self, hidden: int, heads: int, dropout: float) -> None:
        super().__init__()
        self.cross = nn.MultiheadAttention(hidden, heads, dropout=dropout, batch_first=True)
        self.self_attention = nn.MultiheadAttention(hidden, heads, dropout=dropout, batch_first=True)
        self.norm1, self.norm2, self.norm3 = nn.LayerNorm(hidden), nn.LayerNorm(hidden), nn.LayerNorm(hidden)
        self.ffn = nn.Sequential(nn.Linear(hidden, hidden * 4), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden * 4, hidden))
        self.dropout = nn.Dropout(dropout)

    def forward(self, query: torch.Tensor, evidence: torch.Tensor, padding: torch.Tensor) -> torch.Tensor:
        value, _ = self.cross(query, evidence, evidence, key_padding_mask=padding, need_weights=False)
        query = self.norm1(query + self.dropout(value))
        value, _ = self.self_attention(query, query, query, need_weights=False)
        query = self.norm2(query + self.dropout(value))
        return self.norm3(query + self.dropout(self.ffn(query)))


class DenseAdapter(nn.Module):
    def __init__(self, input_dim: int, hidden: int, name: str) -> None:
        super().__init__()
        self.name = name
        self.projection = nn.Sequential(nn.Linear(input_dim, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Linear(hidden, hidden))
        self.missing = nn.Parameter(torch.zeros(1, 1, hidden))

    def forward(self, values: torch.Tensor | None, present: torch.Tensor) -> ModalityTokens:
        if values is None:
            return ModalityTokens(self.missing.expand(len(present), -1, -1), torch.zeros(len(present), 1, device=present.device, dtype=torch.bool), self.name)
        token = self.projection(values).unsqueeze(1)
        token = torch.where(present[:, None, None], token, self.missing.expand(len(present), -1, -1))
        return ModalityTokens(token, present[:, None].bool(), self.name)


class HierarchicalWSIEncoder(nn.Module):
    """Patch -> slide -> patient hierarchy; the sampler owns patch coverage."""

    def __init__(self, config: V2ModelConfig) -> None:
        super().__init__()
        h = config.hidden_dim
        self.patch = nn.Sequential(nn.Linear(config.patch_dim, h), nn.LayerNorm(h), nn.GELU(), nn.Linear(h, h))
        self.coordinates = nn.Sequential(nn.Linear(2, h), nn.LayerNorm(h), nn.GELU(), nn.Linear(h, h))
        self.local = SoftSlotPool(h, config.local_slots)
        self.slide = SoftSlotPool(h, config.slide_slots)
        self.patient = SoftSlotPool(h, config.patient_slots)

    def forward(self, patches: torch.Tensor, patch_mask: torch.Tensor, slide_ids: torch.Tensor, coordinates: torch.Tensor | None = None) -> tuple[ModalityTokens, dict[str, torch.Tensor]]:
        tokens = self.patch(patches)
        if coordinates is not None:
            scale = coordinates.abs().amax(1, keepdim=True).clamp_min(1.0)
            tokens = tokens + self.coordinates(coordinates / scale)
        rows, masks, entropies = [], [], []
        for batch_index in range(tokens.shape[0]):
            valid = patch_mask[batch_index]
            if not valid.any():
                rows.append(tokens.new_zeros((1, tokens.shape[-1])))
                masks.append(torch.zeros(1, device=tokens.device, dtype=torch.bool))
                continue
            current = tokens[batch_index, valid]
            ids = slide_ids[batch_index, valid]
            slides = []
            for slide_id in torch.unique(ids, sorted=True):
                group = current[ids == slide_id].unsqueeze(0)
                group_mask = torch.ones((1, group.shape[1]), device=tokens.device, dtype=torch.bool)
                pooled, _, attention = self.local(group, group_mask)
                slides.append(pooled.squeeze(0))
                entropies.append((-(attention * attention.clamp_min(1e-12).log()).sum(-1)).mean())
            slide_tokens = torch.cat(slides, dim=0).unsqueeze(0)
            slide_mask = torch.ones((1, slide_tokens.shape[1]), device=tokens.device, dtype=torch.bool)
            pooled, _, _ = self.slide(slide_tokens, slide_mask)
            rows.append(pooled.squeeze(0))
            masks.append(torch.ones(pooled.shape[1], device=tokens.device, dtype=torch.bool))
        width = max(row.shape[0] for row in rows)
        packed = tokens.new_zeros((len(rows), width, tokens.shape[-1]))
        packed_mask = torch.zeros((len(rows), width), device=tokens.device, dtype=torch.bool)
        for index, row in enumerate(rows):
            packed[index, : row.shape[0]] = row
            packed_mask[index, : row.shape[0]] = masks[index]
        pooled, pooled_mask, attention = self.patient(packed, packed_mask)
        return ModalityTokens(pooled, pooled_mask, "wsi"), {
            "local_slot_entropy": torch.stack(entropies).mean() if entropies else tokens.new_zeros(()),
            "patient_slot_effective_tokens": 1.0 / attention.square().sum(-1).clamp_min(1e-12),
        }


class TumorStateV2(nn.Module):
    """Typed fusion with dedicated identity, biology, context and residual states."""

    def __init__(self, config: V2ModelConfig | None = None, clinical_dim: int | None = None, snv_dim: int | None = None, cnv_dim: int | None = None, programme_dim: int = 50) -> None:
        super().__init__()
        self.config = config or V2ModelConfig()
        c = self.config
        self.wsi = HierarchicalWSIEncoder(c)
        self.rna = DenseAdapter(c.rna_dim, c.hidden_dim, "rna")
        self.clinical = DenseAdapter(clinical_dim, c.hidden_dim, "clinical") if clinical_dim else None
        self.snv = DenseAdapter(snv_dim, c.hidden_dim, "snv") if snv_dim else None
        self.cnv = DenseAdapter(cnv_dim, c.hidden_dim, "cnv") if cnv_dim else None
        self.queries = nn.Parameter(torch.empty(1, c.query_slots, c.hidden_dim))
        nn.init.orthogonal_(self.queries.squeeze(0))
        self.blocks = nn.ModuleList([QueryBlock(c.hidden_dim, c.heads, c.dropout) for _ in range(c.layers)])
        self.norm = nn.LayerNorm(c.hidden_dim)
        self.identity = nn.Linear(c.hidden_dim, 256)
        self.biology = nn.Linear(c.hidden_dim, 256)
        self.context = nn.Linear(c.hidden_dim, 128)
        self.patient = nn.Linear(c.hidden_dim, 256)
        self.rna_reconstruction = nn.Linear(256, c.rna_dim)
        self.programme_mean = nn.Linear(c.hidden_dim, programme_dim)
        self.programme_log_variance = nn.Linear(c.hidden_dim, programme_dim)

    def forward(self, batch: dict[str, torch.Tensor], view: str = "full") -> dict[str, torch.Tensor]:
        if view not in {"full", "wsi", "rna"}:
            raise ValueError("view must be one of 'full', 'wsi', or 'rna'")
        c = self.config
        modalities: dict[str, ModalityTokens] = {}
        diagnostics: dict[str, torch.Tensor] = {}
        if view != "rna":
            wsi, diagnostics = self.wsi(batch["patches"], batch["patch_mask"], batch["slide_ids"], batch.get("coordinates"))
            modalities["wsi"] = wsi
        else:
            first = batch["rna"]
            wsi = ModalityTokens(first.new_zeros((len(first), 0, c.hidden_dim)), torch.zeros((len(first), 0), device=first.device, dtype=torch.bool), "wsi")
        if view != "wsi":
            present = batch.get("rna_present", torch.ones(len(wsi.mask), device=wsi.mask.device, dtype=torch.bool))
            modalities["rna"] = self.rna(batch.get("rna"), present)
            for name, adapter in (("clinical", self.clinical), ("snv", self.snv), ("cnv", self.cnv)):
                if adapter is not None:
                    present = batch.get(f"{name}_present", torch.zeros(len(wsi.mask), device=wsi.mask.device, dtype=torch.bool))
                    modalities[name] = adapter(batch.get(name), present)
        evidence = torch.cat([item.tokens for item in modalities.values()], dim=1)
        mask = torch.cat([item.mask for item in modalities.values()], dim=1).bool()
        query = self.queries.expand(len(evidence), -1, -1)
        for block in self.blocks:
            query = block(query, evidence, ~mask)
        query = self.norm(query)
        cursor = 0
        def take(count: int) -> torch.Tensor:
            nonlocal cursor
            value = query[:, cursor : cursor + count].mean(1)
            cursor += count
            return value
        identity, biology, context = take(c.identity_slots), take(c.biology_slots), take(c.context_slots)
        residuals = {name: take(c.residual_slots) for name in ("wsi", "rna", "clinical", "snv", "cnv")}
        uncertainty = take(c.uncertainty_slots)
        z_identity = self.identity(identity)
        return {
            "z_identity": z_identity, "z_biology": self.biology(biology), "z_context": self.context(context), "z_patient": self.patient(query.mean(1)),
            "rna_reconstruction": self.rna_reconstruction(z_identity),
            "programme_mean": self.programme_mean(biology), "programme_log_variance": self.programme_log_variance(biology).clamp(-8, 8),
            "z_uncertainty": uncertainty, **{f"z_{name}_residual": value for name, value in residuals.items()}, **diagnostics,
        }
