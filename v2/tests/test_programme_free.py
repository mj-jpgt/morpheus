"""D1 contract tests: no curated programme loss, no hidden small-batch no-op."""
from __future__ import annotations

import pytest
import torch

from morpheus.v2.losses import paired_infonce_with_memory
from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.runner import (_overfit_programme_only_actual,
                                _overfit_programme_free_contrastive,
                                _require_programme_free_overfit,
                                _trained_states_for_profile)
from morpheus.v2.training import V2LossSchedule, V2Trainer


def _config(hidden: int = 32) -> V2ModelConfig:
    return V2ModelConfig(patch_dim=8, rna_dim=4, hidden_dim=hidden, heads=4, layers=1,
                         local_slots=4, slide_slots=2, patient_slots=2)


def _big_batch(n: int, k: int, seed: int) -> dict[str, torch.Tensor]:
    generator = torch.Generator().manual_seed(seed)
    return {
        "patches": torch.randn(n, 6, 8, generator=generator),
        "patch_mask": torch.ones(n, 6, dtype=torch.bool),
        "slide_ids": torch.zeros(n, 6, dtype=torch.long),
        "rna": torch.randn(n, 4, generator=generator), "rna_present": torch.ones(n, dtype=torch.bool),
        "programme_target": torch.randn(n, k, generator=generator), "programme_present": torch.ones(n, dtype=torch.bool),
        "programme_target_mask": torch.ones(n, k, dtype=torch.bool),
        "programme_positive_mask": torch.ones(n, n, dtype=torch.bool),
        "programme_neighbor_indices": torch.zeros(n, 1, dtype=torch.long), "indices": torch.arange(n) + seed * 10,
    }


def test_programme_free_weight_contract_and_artifact_contract() -> None:
    free = V2LossSchedule(objective_profile="programme_free").weights(epoch=99)
    programme = V2LossSchedule(objective_profile="programme_only").weights(epoch=99)
    assert free["programme"] == free["neighbourhood"] == free["supcon"] == 0.0
    assert free["identity"] == free["fusion_identity"] == free["patient_consistency"] == free["semantic"] == 0.0
    assert free["decorrelation"] == programme["decorrelation"]
    assert free["biology_contrastive"] > 0.0
    assert free["biology_full_consistency"] > 0.0
    assert _trained_states_for_profile("programme_free") == ["wsi_biology", "rna_biology", "full_biology"]


def test_memory_paired_infonce_refuses_a_real_small_batch_without_negatives() -> None:
    wsi, rna, ids = torch.randn(1, 8), torch.randn(1, 8), torch.tensor([7])
    with pytest.raises(RuntimeError, match="effective negatives"):
        paired_infonce_with_memory(wsi, rna, ids, None, None, None)


def test_memory_paired_infonce_is_live_with_id_aware_queue() -> None:
    torch.manual_seed(3)
    wsi, rna = torch.randn(2, 8, requires_grad=True), torch.randn(2, 8, requires_grad=True)
    ids = torch.tensor([7, 8])
    queue_wsi, queue_rna, queue_ids = torch.randn(10, 8), torch.randn(10, 8), torch.arange(7, 17)
    loss, negatives = paired_infonce_with_memory(wsi, rna, ids, queue_wsi, queue_rna, queue_ids)
    assert torch.isfinite(loss) and int(negatives.min()) >= 8
    loss.backward()
    assert wsi.grad is not None and wsi.grad.abs().sum() > 0
    assert rna.grad is not None and rna.grad.abs().sum() > 0


def test_programme_free_step_has_live_contrastive_gradient_at_small_batch() -> None:
    torch.manual_seed(4)
    model = TumorStateV2(_config(hidden=32), programme_dim=8)
    trainer = V2Trainer(model, torch.optim.AdamW(model.parameters(), lr=1e-3),
                        V2LossSchedule(objective_profile="programme_free", warmup_epochs=0), "cpu")
    batches = [_big_batch(n=2, k=8, seed=seed) for seed in range(6)]
    assert trainer.prime_biology_memory(batches, minimum_unique_keys=9) >= 9
    loss, metrics, _ = trainer.step(_big_batch(n=2, k=8, seed=99), epoch=1)
    assert metrics["biology_contrastive"] > 1e-4
    assert metrics["biology_full_consistency"] > 1e-4
    assert metrics["biology_contrastive_effective_negatives_min"] >= 8
    loss.backward()
    assert all(value > 0.0 for value in trainer._gradient_group_norms().values())


def test_programme_free_overfit_uses_the_actual_model_path() -> None:
    torch.manual_seed(5)
    model = TumorStateV2(_config(hidden=32), programme_dim=8)
    schedule = V2LossSchedule(objective_profile="programme_free", warmup_epochs=0, decorrelation_after_warmup=0.0)
    batches = [_big_batch(n=2, k=8, seed=seed) for seed in range(12)]
    result = _overfit_programme_free_contrastive(model, schedule, batches, "cpu", steps=3, minimum_memory_keys=16)
    assert result["objective_scope"] == "actual_v2_encoder_and_biology_path_without_decorrelation_floor"
    assert result["memory_unique_keys"] >= 16
    assert all(value > 0.0 for value in result["gradient_norms_first"].values())


def test_programme_only_overfit_uses_the_actual_model_path() -> None:
    torch.manual_seed(6)
    model = TumorStateV2(_config(hidden=32), programme_dim=8)
    schedule = V2LossSchedule(objective_profile="programme_only", warmup_epochs=0,
                               neighbourhood_after_warmup=0.0, supcon_after_warmup=0.0,
                               decorrelation_after_warmup=0.0)
    result = _overfit_programme_only_actual(model, schedule, [_big_batch(n=2, k=8, seed=6)], "cpu", steps=3)
    assert result["objective_scope"] == "actual_v2_encoder_and_programme_path_without_decorrelation_floor"
    assert all(value > 0.0 for value in result["gradient_norms_first"].values())


def test_programme_free_overfit_gate_fails_closed() -> None:
    with pytest.raises(RuntimeError, match="G2.6"):
        _require_programme_free_overfit({
            "initial_loss": 3.0, "final_loss": 2.5, "relative_reduction": 0.1,
            "final_biology_contrastive": 2.4, "final_full_consistency": 0.1,
        })


def test_programme_free_rejects_a_detached_biology_state() -> None:
    class DetachedBiology(torch.nn.Module):
        def __init__(self, inner: torch.nn.Module) -> None:
            super().__init__(); self.inner = inner

        def forward(self, value: torch.Tensor) -> torch.Tensor:
            return self.inner(value).detach()

    model = TumorStateV2(_config(hidden=32), programme_dim=8)
    model.biology = DetachedBiology(model.biology)
    trainer = V2Trainer(model, torch.optim.AdamW(model.parameters(), lr=1e-3),
                        V2LossSchedule(objective_profile="programme_free", warmup_epochs=0,
                                       decorrelation_after_warmup=0.0), "cpu")
    batches = [_big_batch(n=2, k=8, seed=seed) for seed in range(6)]
    trainer.prime_biology_memory(batches, minimum_unique_keys=9)
    loss, _, _ = trainer.step(_big_batch(n=2, k=8, seed=99), epoch=1)
    # The complete model still computes unused auxiliary tensors, so checking
    # `loss.requires_grad` is insufficient.  The G2 assertion is that the
    # actual D1 biology group receives no gradient and would therefore fail
    # the real overfit/liveness check.
    loss.backward()
    assert trainer._gradient_group_norms()["biology_programme"] == 0.0
