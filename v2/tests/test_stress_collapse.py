"""Fast CPU stress tests that fail before an expensive full run.

These guard the failure modes the existing suite does NOT catch (see
``v2/research/B2_implementation_and_audit.md`` §3): representational/rank
collapse, NaN in the KL / Gaussian-NLL chain, heteroscedastic mode collapse,
no-op losses, and anchor-gate saturation.

The load-bearing utility is :func:`effective_rank`; it is the instrument used to
attribute the biology-head collapse fix (F-R1 normalize, F-R2 covariance
decorrelation). Tiny synthetic fixtures cannot show an absolute rank-256 floor,
so the collapse checks assert *relative / mechanistic* properties on a locally
built, higher-rank synthetic batch, plus a deterministic detector self-test.
"""
from __future__ import annotations

import torch

from morpheus.v2 import losses as L
from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.training import V2LossSchedule, V2Trainer


def effective_rank(x: torch.Tensor) -> float:
    """Roy-Vetterli effective rank: exp(entropy of L1-normalised singular values).

    Computed on the column-centred batch (feature covariance spectrum). Returns
    0.0 for an all-constant batch.
    """
    x = x - x.mean(dim=0, keepdim=True)
    singular = torch.linalg.svdvals(x.double())
    singular = singular[singular > 1e-12]
    if singular.numel() == 0:
        return 0.0
    p = singular / singular.sum()
    return float(torch.exp(-(p * p.log()).sum()))


def _config(hidden: int = 32, anchor: bool = False) -> V2ModelConfig:
    return V2ModelConfig(
        patch_dim=8, rna_dim=4, hidden_dim=hidden, heads=4, layers=1,
        local_slots=4, slide_slots=2, patient_slots=2, use_mlp_clip_anchor=anchor,
    )


def _big_batch(n: int = 16, k: int = 8, seed: int = 0, anchor: bool = False,
               low_rank_targets: bool = False) -> dict[str, torch.Tensor]:
    """A B>=8 synthetic batch with genuine cross-patient variance.

    Mirrors the key/shape contract of ``tests.test_v21_model._batch`` but scaled
    up so the effective-rank instrument has something to measure.
    """
    g = torch.Generator().manual_seed(seed)
    patches = torch.randn(n, 6, 8, generator=g)
    patch_mask = torch.ones(n, 6, dtype=torch.bool)
    slide_ids = torch.zeros(n, 6, dtype=torch.long)
    slide_ids[:, 3:] = 1
    if low_rank_targets:
        basis = torch.randn(2, k, generator=g)
        targets = torch.randn(n, 2, generator=g) @ basis
    else:
        targets = torch.randn(n, k, generator=g)
    neighbours = torch.stack([torch.tensor([(i + 1) % n, (i - 1) % n]) for i in range(n)])
    batch = {
        "patches": patches, "patch_mask": patch_mask, "slide_ids": slide_ids,
        "rna": torch.randn(n, 4, generator=g), "rna_present": torch.ones(n, dtype=torch.bool),
        "programme_target": targets, "programme_present": torch.ones(n, dtype=torch.bool),
        "programme_target_mask": torch.ones(n, k, dtype=torch.bool),
        "programme_positive_mask": (torch.rand(n, n, generator=g) > 0.6),
        "programme_neighbor_indices": neighbours,
        "indices": torch.arange(n),
    }
    if anchor:
        batch["mlp_clip_anchor_wsi"] = torch.nn.functional.normalize(torch.randn(n, 256, generator=g), dim=-1)
        batch["mlp_clip_anchor_rna"] = torch.nn.functional.normalize(torch.randn(n, 256, generator=g), dim=-1)
    return batch


def _trainer(model: TumorStateV2, profile: str = "full") -> tuple[V2Trainer, torch.optim.Optimizer]:
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    trainer = V2Trainer(model, optimizer, V2LossSchedule(objective_profile=profile, warmup_epochs=0),
                        "cpu", gradient_diagnostics_every=0)
    return trainer, optimizer


# --- (R0 gate) the instrument itself: must FAIL on rank-deficient, PASS on spread ---

def test_effective_rank_detects_collapse() -> None:
    torch.manual_seed(0)
    collapsed = torch.randn(64, 2) @ torch.randn(2, 256)  # true rank 2
    spread = torch.randn(64, 256)
    assert effective_rank(collapsed) < 3.0
    assert effective_rank(spread) > 20.0
    assert effective_rank(torch.ones(16, 256)) == 0.0  # all-constant -> 0 after centring


# --- numerical health: NaN / Inf guards around the KL and Gaussian-NLL chains ---

def test_kl_is_finite_under_degenerate_rows() -> None:
    state = torch.randn(6, 16); state[0] = 0.0            # zero row -> _norm divides by ~0
    targets = torch.randn(6, 8); targets[1] = 0.0
    assert torch.isfinite(L.programme_neighbourhood_loss(state, targets))
    assert torch.isfinite(L.supervised_programme_contrastive(state, torch.rand(6, 6) > 0.5))


def test_gaussian_nll_is_finite_under_extreme_variance_and_masked_nan() -> None:
    mean = torch.zeros(4, 3)
    target = torch.randn(4, 3); target[0, 0] = float("nan")
    mask = torch.ones(4, 3, dtype=torch.bool); mask[0, 0] = False  # the NaN cell is unavailable
    log_variance = torch.full((4, 3), 40.0)
    assert torch.isfinite(L.gaussian_nll(mean, log_variance, target, mask))


def test_full_training_step_metrics_are_finite() -> None:
    torch.manual_seed(0)
    model = TumorStateV2(_config(), programme_dim=8)
    trainer, _ = _trainer(model)
    loss, _, metrics = trainer.step(_big_batch(), epoch=1)
    assert torch.isfinite(loss)
    for name, value in metrics.items():
        tensor = value if torch.is_tensor(value) else torch.tensor(float(value))
        assert torch.isfinite(tensor).all(), name


def test_programme_log_variance_not_pinned_at_clamp() -> None:
    """A biology head explaining every programme as pure noise pins log-variance at +8."""
    torch.manual_seed(0)
    model = TumorStateV2(_config(), programme_dim=8)
    out = model(_big_batch(), "full")
    assert out["programme_log_variance"].max().item() < 7.9


# --- collapse instrument on the real model (lenient; fixes should raise this) ---

def test_biology_head_effective_rank_is_not_degenerate() -> None:
    torch.manual_seed(0)
    model = TumorStateV2(_config(hidden=32), programme_dim=8)
    trainer, optimizer = _trainer(model, profile="programme_only")
    batch = _big_batch(n=16, k=8)
    for _ in range(15):
        optimizer.zero_grad()
        loss, _, _ = trainer.step(batch, epoch=1)
        loss.backward()
        optimizer.step()
    out = model(batch, "full")
    assert effective_rank(out["z_biology"]) > 1.5  # documents the instrument; F-R1/R2 raise it


# --- no-op / dead-loss guard: every declared-active head must get gradient ---

def test_every_active_head_receives_gradient_in_full_profile() -> None:
    torch.manual_seed(0)
    model = TumorStateV2(_config(), programme_dim=8)
    trainer, _ = _trainer(model)
    loss, _, _ = trainer.step(_big_batch(), epoch=1)
    loss.backward()
    for head in ("identity", "biology", "patient"):
        module = getattr(model, head)
        assert any(p.grad is not None and p.grad.abs().sum() > 0 for p in module.parameters()), head


# --- anchor gate must not be saturated; the residual must be able to move ---

def test_anchor_gate_not_saturated_and_residual_can_move() -> None:
    torch.manual_seed(0)
    model = TumorStateV2(_config(anchor=True), programme_dim=8)
    trainer, _ = _trainer(model)
    loss, _, _ = trainer.step(_big_batch(anchor=True), epoch=1)
    loss.backward()
    assert model.anchor_residual_scale.grad is not None  # a real gradient path -> can learn to move
    out = model(_big_batch(anchor=True), "wsi")
    gate = out["anchor_gate_mean"].item()
    assert 0.0 < gate < 1.0


# --- slot accounting: take() must consume every query slot exactly once ---

def test_query_slots_are_fully_accounted() -> None:
    config = V2ModelConfig()
    assert config.query_slots == (
        config.identity_slots + config.biology_slots + config.context_slots
        + 5 * config.residual_slots + config.uncertainty_slots
    )
    model = TumorStateV2(_config(), programme_dim=8)
    # self.queries is [1, query_slots, hidden], expanded over the batch in forward().
    assert model.queries.shape[-2] == _config().query_slots
