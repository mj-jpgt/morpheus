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
    # A direct global path protects the strong frozen H-Optimus signal while
    # slots add regional and cross-slide evidence.  The old 8->8->16 hierarchy
    # compressed thousands of tiles before patient fusion.
    local_slots: int = 64
    slide_slots: int = 32
    patient_slots: int = 128
    identity_slots: int = 4
    biology_slots: int = 4
    context_slots: int = 2
    residual_slots: int = 4
    uncertainty_slots: int = 2
    semantic_dim: int = 0
    # When supplied, frozen all-patch MLP-CLIP embeddings are a retrieval
    # skip connection.  A zero-initialised residual coefficient makes the
    # anchored model exactly equal to the baseline at step zero.
    use_mlp_clip_anchor: bool = False
    # The refiner is deliberately bounded.  At initialization its global
    # coefficient is zero, so an anchored identity state is exactly the
    # supplied frozen MLP-CLIP state.
    anchor_residual_bound: float = 0.25
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
        self.global_summary = nn.Sequential(
            nn.Linear(config.patch_dim * 2, h), nn.LayerNorm(h), nn.GELU(), nn.Linear(h, h)
        )
        self.coordinates = nn.Sequential(nn.Linear(2, h), nn.LayerNorm(h), nn.GELU(), nn.Linear(h, h))
        self.local = SoftSlotPool(h, config.local_slots)
        self.slide = SoftSlotPool(h, config.slide_slots)
        self.patient = SoftSlotPool(h, config.patient_slots)

    def forward(self, patches: torch.Tensor, patch_mask: torch.Tensor,
                slide_ids: torch.Tensor, coordinates: torch.Tensor | None = None,
                coordinate_present: torch.Tensor | None = None) -> tuple[ModalityTokens, dict[str, torch.Tensor]]:
        """Encode an uncapped ragged patch bag.

        ``coordinate_present`` is patient-level on purpose.  Callers may pad
        coordinate tensors for a mixed batch, but an absent/constant metadata
        row must never be interpreted as a real all-zero spatial layout.
        """
        tokens = self.patch(patches)
        if coordinates is not None:
            if coordinate_present is None:
                coordinate_present = torch.ones(len(tokens), device=tokens.device, dtype=torch.bool)
            if coordinate_present.shape != (len(tokens),):
                raise ValueError("coordinate_present must be a patient-level [batch] mask")
            scale = coordinates.abs().amax(1, keepdim=True).clamp_min(1.0)
            encoded = self.coordinates(coordinates / scale)
            tokens = tokens + encoded * coordinate_present[:, None, None].to(encoded.dtype)
        rows, masks, globals_, entropies = [], [], [], []
        for batch_index in range(tokens.shape[0]):
            valid = patch_mask[batch_index]
            if not valid.any():
                rows.append(tokens.new_zeros((1, tokens.shape[-1])))
                masks.append(torch.zeros(1, device=tokens.device, dtype=torch.bool))
                globals_.append(tokens.new_zeros((1, tokens.shape[-1])))
                continue
            current = tokens[batch_index, valid]
            raw = patches[batch_index, valid]
            # Mean and dispersion are a lossless-in-coverage global skip path;
            # learned slots are a residual refinement, not the only evidence.
            summary = torch.cat((raw.mean(0), raw.std(0, unbiased=False)), dim=-1)
            globals_.append(self.global_summary(summary).unsqueeze(0))
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
        global_tokens = torch.stack(globals_, dim=0)
        tokens_out = torch.cat((pooled, global_tokens), dim=1)
        mask_out = torch.cat((pooled_mask, torch.ones((len(rows), 1), dtype=torch.bool, device=pooled.device)), dim=1)
        return ModalityTokens(tokens_out, mask_out, "wsi"), {
            "local_slot_entropy": torch.stack(entropies).mean() if entropies else tokens.new_zeros(()),
            "patient_slot_effective_tokens": 1.0 / attention.square().sum(-1).clamp_min(1e-12),
            "coordinate_present_fraction": (
                coordinate_present.float().mean() if coordinate_present is not None else tokens.new_zeros(())
            ),
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
        self.patient = nn.Sequential(nn.LayerNorm(c.hidden_dim), nn.Linear(c.hidden_dim, 256))
        self.rna_reconstruction = nn.Linear(256, c.rna_dim)
        # Programme regression is deliberately downstream of the exported
        # biology state.  The old heads consumed the pre-projection hidden
        # vector, leaving `z_biology` without a direct regression gradient
        # whenever neighbourhood positives were sparse.
        self.programme_mean = nn.Linear(256, programme_dim)
        self.programme_log_variance = nn.Linear(256, programme_dim)
        self.semantic = nn.Linear(c.hidden_dim, c.semantic_dim) if c.semantic_dim else None
        self.anchor_residual_scale = nn.Parameter(torch.zeros(())) if c.use_mlp_clip_anchor else None
        # This gate determines where a contextual correction is useful; the
        # separate zero-initialised global coefficient preserves the teacher
        # exactly at step zero.  It is intentionally not a second teacher loss.
        self.anchor_gate = nn.Sequential(nn.LayerNorm(c.hidden_dim), nn.Linear(c.hidden_dim, 1)) if c.use_mlp_clip_anchor else None

    def forward(self, batch: dict[str, torch.Tensor], view: str = "full") -> dict[str, torch.Tensor]:
        if view not in {"full", "wsi", "rna"}:
            raise ValueError("view must be one of 'full', 'wsi', or 'rna'")
        c = self.config
        modalities: dict[str, ModalityTokens] = {}
        diagnostics: dict[str, torch.Tensor] = {}
        if view != "rna":
            wsi, diagnostics = self.wsi(
                batch["patches"], batch["patch_mask"], batch["slide_ids"],
                batch.get("coordinates"), batch.get("coordinate_present"),
            )
            modalities["wsi"] = wsi
        else:
            first = batch["rna"]
            wsi = ModalityTokens(first.new_zeros((len(first), 0, c.hidden_dim)), torch.zeros((len(first), 0), device=first.device, dtype=torch.bool), "wsi")
        if view != "wsi":
            present = batch.get("rna_present", torch.ones(len(wsi.mask), device=wsi.mask.device, dtype=torch.bool))
            modalities["rna"] = self.rna(batch.get("rna"), present)
            # The RNA retrieval view must stay RNA-only.  Extra modalities are
            # admitted only to the explicitly requested fused patient view.
            if view == "full":
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
        identity_residual = nn.functional.normalize(self.identity(identity), dim=-1)
        anchor = None
        if self.anchor_residual_scale is not None:
            if view == "wsi":
                anchor = batch.get("mlp_clip_anchor_wsi")
            elif view == "rna":
                anchor = batch.get("mlp_clip_anchor_rna")
            else:
                wsi_anchor, rna_anchor = batch.get("mlp_clip_anchor_wsi"), batch.get("mlp_clip_anchor_rna")
                if wsi_anchor is not None and rna_anchor is not None:
                    anchor = nn.functional.normalize(wsi_anchor + rna_anchor, dim=-1)
            if anchor is None:
                raise ValueError("anchored V2 requires MLP-CLIP anchors for every requested view")
            # Keep the learned correction bounded so anchored retrieval cannot
            # drift arbitrarily far from the fair all-patch baseline.
            anchor = nn.functional.normalize(anchor, dim=-1)
            residual_scale = c.anchor_residual_bound * torch.tanh(self.anchor_residual_scale)
            assert self.anchor_gate is not None
            gate = torch.sigmoid(self.anchor_gate(identity)).squeeze(-1)
            correction = residual_scale * gate[:, None] * identity_residual
            z_identity = nn.functional.normalize(anchor + correction, dim=-1)
        else:
            z_identity = identity_residual
        # Normalise the biology state once so the programme regression heads
        # (Gaussian NLL) and the rank-spreading objectives (z_biology / KL /
        # supcon) operate on the SAME geometry.  Previously the NLL saw the raw
        # projection while z_biology/KL/supcon saw the unit-norm state, so the
        # dominant low-rank pull acted on a different manifold than the weak
        # spreaders (F-R1; see v2/research/B2_implementation_and_audit.md).
        biology_state = nn.functional.normalize(self.biology(biology), dim=-1)
        result = {
            "z_identity": z_identity,
            "z_biology": biology_state,
            "z_context": self.context(context),
            "z_patient": nn.functional.normalize(self.patient(query.mean(1)), dim=-1),
            "rna_reconstruction": self.rna_reconstruction(z_identity),
            "programme_mean": self.programme_mean(biology_state), "programme_log_variance": self.programme_log_variance(biology_state).clamp(-8, 8),
            "z_uncertainty": uncertainty, **{f"z_{name}_residual": value for name, value in residuals.items()}, **diagnostics,
        }
        if self.semantic is not None:
            # Only the WSI path is eligible for pathology-text semantic use.
            result["z_semantic"] = self.semantic(identity)
        if self.anchor_residual_scale is not None:
            result["anchor_residual_scale"] = (c.anchor_residual_bound * torch.tanh(self.anchor_residual_scale)).detach()
            result["anchor_gate_mean"] = gate.detach().mean()
            result["anchor_correction_norm"] = correction.detach().norm(dim=-1).mean()
        return result
