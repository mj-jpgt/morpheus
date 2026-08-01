from __future__ import annotations

import torch

from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.runner import _selected_checkpoint
from morpheus.v2.runner import _trained_states_for_profile
from morpheus.v2.slide_pretraining import SlidePretrainingConfig, SlidePretrainingObjective
from morpheus.v2.training import V2LossSchedule, V2Trainer


def _config(anchor: bool = False) -> V2ModelConfig:
    return V2ModelConfig(
        patch_dim=8, rna_dim=4, hidden_dim=16, heads=4, layers=1,
        local_slots=4, slide_slots=2, patient_slots=2, use_mlp_clip_anchor=anchor,
    )


def _batch(anchor: bool = False) -> dict[str, torch.Tensor]:
    value = {
        "patches": torch.randn(3, 7, 8),
        "patch_mask": torch.tensor([[1] * 7, [1, 1, 1, 0, 0, 0, 0], [1, 1, 1, 1, 1, 0, 0]], dtype=torch.bool),
        "slide_ids": torch.tensor([[0, 0, 1, 1, 1, 2, 2], [0, 0, 0, 0, 0, 0, 0], [0, 1, 1, 1, 2, 0, 0]]),
        "rna": torch.randn(3, 4), "rna_present": torch.ones(3, dtype=torch.bool),
        "programme_target": torch.randn(3, 3), "programme_present": torch.ones(3, dtype=torch.bool),
        "programme_target_mask": torch.tensor([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=torch.bool),
        "programme_positive_mask": torch.tensor([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=torch.bool),
        "programme_neighbor_indices": torch.tensor([[1, -1], [0, 2], [1, -1]]),
        "indices": torch.arange(3),
    }
    if anchor:
        value["mlp_clip_anchor_wsi"] = torch.nn.functional.normalize(torch.randn(3, 256), dim=-1)
        value["mlp_clip_anchor_rna"] = torch.nn.functional.normalize(torch.randn(3, 256), dim=-1)
    return value


def test_anchor_is_exact_baseline_at_initialization() -> None:
    model = TumorStateV2(_config(anchor=True), programme_dim=3).eval()
    batch = _batch(anchor=True)
    with torch.no_grad():
        output = model(batch, "wsi")
    assert torch.allclose(output["z_identity"], batch["mlp_clip_anchor_wsi"], atol=1e-6)
    assert output["anchor_correction_norm"].item() == 0.0


def test_missing_coordinates_are_not_zero_coordinate_layout() -> None:
    model = TumorStateV2(_config(), programme_dim=3).eval()
    batch = _batch()
    coordinates = torch.randn(3, 7, 2)
    batch["coordinates"] = coordinates
    batch["coordinate_present"] = torch.tensor([True, False, True])
    no_coordinates = {key: value for key, value in batch.items() if key not in {"coordinates", "coordinate_present"}}
    with torch.no_grad():
        mixed = model(batch, "wsi")["z_identity"]
        absent = model(no_coordinates, "wsi")["z_identity"]
    assert torch.allclose(mixed[1], absent[1], atol=1e-6)


def test_programme_and_patient_heads_receive_training_gradient() -> None:
    model = TumorStateV2(_config(), programme_dim=3)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    trainer = V2Trainer(model, optimizer, V2LossSchedule(warmup_epochs=0), "cpu", gradient_diagnostics_every=0)
    loss, _, _ = trainer.step(_batch(), epoch=1)
    loss.backward()
    assert model.programme_mean.weight.grad is not None and model.programme_mean.weight.grad.abs().sum() > 0
    assert model.patient[-1].weight.grad is not None and model.patient[-1].weight.grad.abs().sum() > 0


def test_slide_pretraining_uses_only_patch_batch_fields() -> None:
    model = TumorStateV2(_config(), programme_dim=3)
    objective = SlidePretrainingObjective(8, 16, SlidePretrainingConfig(target_dim=6))
    batch = _batch()
    loss, metrics = objective.loss(model.wsi, {key: batch[key] for key in ("patches", "patch_mask", "slide_ids")})
    assert torch.isfinite(loss)
    assert set(metrics) == {"pretrain_masked_bag", "pretrain_view_consistency", "pretrain_masked_fraction"}


def test_objective_profiles_disable_only_declared_loss_families() -> None:
    identity = V2LossSchedule(objective_profile="identity_only", warmup_epochs=0).weights(1)
    assert identity["identity"] > 0
    assert identity["fusion_identity"] > 0
    assert all(identity[key] == 0 for key in identity if key not in {"identity", "fusion_identity"})
    programme = V2LossSchedule(objective_profile="programme_only", warmup_epochs=0).weights(1)
    assert all(programme[key] > 0 for key in ("programme", "neighbourhood", "supcon"))
    assert all(programme[key] == 0 for key in ("identity", "fusion_identity", "patient_consistency", "rna_reconstruction", "separation", "semantic"))


def test_profile_declarations_match_loss_gradient_families() -> None:
    assert _trained_states_for_profile("identity_only") == ["wsi_identity", "rna_identity", "full_identity"]
    assert _trained_states_for_profile("programme_only") == ["wsi_biology", "rna_biology", "full_biology"]
    assert "full_patient" in _trained_states_for_profile("full")
    for profile, active, inactive in (
        ("identity_only", "identity", "biology"),
        ("programme_only", "biology", "identity"),
        ("full", "patient", None),
    ):
        model = TumorStateV2(_config(), programme_dim=3)
        trainer = V2Trainer(model, torch.optim.AdamW(model.parameters()), V2LossSchedule(objective_profile=profile, warmup_epochs=0), "cpu", gradient_diagnostics_every=0)
        loss, _, _ = trainer.step(_batch(), epoch=1)
        loss.backward()
        module = getattr(model, active)
        assert any(parameter.grad is not None and parameter.grad.abs().sum() > 0 for parameter in module.parameters())
        if inactive is not None:
            module = getattr(model, inactive)
            assert all(parameter.grad is None or parameter.grad.abs().sum() == 0 for parameter in module.parameters())


def test_each_programme_only_biology_view_has_direct_projection_gradient() -> None:
    """The exported biology projection cannot depend on sparse neighbour pairs."""
    for view in ("wsi", "rna", "full"):
        model = TumorStateV2(_config(), programme_dim=3)
        output = model(_batch(), view)
        loss = torch.nn.functional.mse_loss(output["programme_mean"], _batch()["programme_target"])
        loss.backward()
        assert any(parameter.grad is not None and parameter.grad.abs().sum() > 0 for parameter in model.biology.parameters()), view


def test_pretrain_only_stage_selects_its_encoder_checkpoint(tmp_path) -> None:
    checkpoint = tmp_path / "slide_pretrain_best.pt"
    checkpoint.write_bytes(b"checkpoint")
    selected = _selected_checkpoint(tmp_path, {"trained": True, "checkpoint": str(checkpoint)}, fine_tune_epochs=0)
    assert selected == checkpoint
