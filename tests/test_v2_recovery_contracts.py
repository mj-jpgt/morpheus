"""Regression tests for the V2 recovery contract.

These tests are deliberately small enough to run locally without a patch store
or any downloaded foundation-model weights.
"""
from __future__ import annotations

import numpy as np
import pytest
import torch
import json

from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.plip import load_plip_teacher_cache
from morpheus.v2.training import ProgrammeMemoryBank, V2LossSchedule, V2Trainer
from morpheus.v2.preflight import validate_runtime_split
from morpheus.v2.contracts import require_state, validate_artifact


def _model_and_batch():
    config = V2ModelConfig(patch_dim=4, rna_dim=3, hidden_dim=8, heads=2, layers=1,
                           local_slots=2, slide_slots=2, patient_slots=2,
                           identity_slots=2, biology_slots=2, context_slots=1,
                           residual_slots=1, uncertainty_slots=1, use_mlp_clip_anchor=True)
    model = TumorStateV2(config, clinical_dim=2, programme_dim=2)
    anchors = torch.randn(3, 256)
    batch = {
        "patches": torch.randn(3, 4, 4), "patch_mask": torch.ones(3, 4, dtype=torch.bool),
        "slide_ids": torch.zeros(3, 4, dtype=torch.long), "rna": torch.randn(3, 3),
        "rna_present": torch.ones(3, dtype=torch.bool), "clinical": torch.randn(3, 2),
        "clinical_present": torch.ones(3, dtype=torch.bool), "mlp_clip_anchor_wsi": anchors,
        "mlp_clip_anchor_rna": torch.randn(3, 256), "programme_target": torch.randn(3, 2),
        "programme_present": torch.ones(3, dtype=torch.bool),
        "programme_positive_mask": torch.tensor([[False, True, False], [True, False, False], [False, False, False]]),
    }
    return model, batch


def test_anchored_wsi_identity_starts_as_frozen_mlp_clip_output():
    model, batch = _model_and_batch()
    with torch.no_grad():
        output = model(batch, "wsi")["z_identity"]
    expected = torch.nn.functional.normalize(batch["mlp_clip_anchor_wsi"], dim=-1)
    assert torch.allclose(output, expected, atol=1e-6)


def test_patient_and_biology_heads_receive_gradients():
    model, batch = _model_and_batch()
    trainer = V2Trainer(model, torch.optim.AdamW(model.parameters(), lr=1e-3),
                        V2LossSchedule(warmup_epochs=0), "cpu")
    loss, _, _ = trainer.step(batch, epoch=1)
    loss.backward()
    assert model.patient[1].weight.grad is not None and model.patient[1].weight.grad.abs().sum() > 0
    assert model.biology.weight.grad is not None and model.biology.weight.grad.abs().sum() > 0


def test_programme_memory_supplies_cross_batch_positives_and_round_trips():
    bank = ProgrammeMemoryBank(capacity=8)
    bank.update(torch.tensor([[1.0, 0.0], [0.0, 1.0]]), torch.tensor([10, 11]))
    loss, active = bank.contrastive(torch.tensor([[1.0, 0.1]]), torch.tensor([[10, -1]]))
    assert active == 1 and torch.isfinite(loss) and loss >= 0
    restored = ProgrammeMemoryBank(); restored.load_state_dict(bank.state_dict(), "cpu")
    recovered, recovered_active = restored.contrastive(torch.tensor([[1.0, 0.1]]), torch.tensor([[10, -1]]))
    assert recovered_active == 1 and torch.allclose(loss, recovered)


def test_strict_split_requires_exact_loaded_paired_cohort(tmp_path):
    train = [f"D{i}" for i in range(11)]
    val = [f"V{i}" for i in range(11)]
    test = [f"T{i}" for i in range(22)]
    path = tmp_path / "split.json"
    path.write_text(json.dumps({"patient_ids":{"train":train,"val":val,"test":test},
                                "actual_train_test_cancer_counts":[11,22]}))
    ids = train + val + test
    cancers = [f"C{i}" for i in range(11)] + [f"C{i}" for i in range(11)] + [f"H{i}" for i in range(22)]
    labels = ["train"] * 11 + ["val"] * 11 + ["test"] * 22
    assert validate_runtime_split(path, ids, cancers, labels)["cancer_counts"] == {"development":11,"heldout":22}
    with pytest.raises(ValueError, match="mismatch"):
        validate_runtime_split(path, ids[:-1], cancers[:-1], labels[:-1])


def test_legacy_baseline_identity_aliases_are_explicit_but_patient_is_not(tmp_path):
    path = tmp_path / "baseline.npz"
    np.savez_compressed(path, patient_ids=np.asarray(["P1", "P2"]), cancers=np.asarray(["A", "B"]),
                        split=np.asarray(["train", "test"]), wsi=np.ones((2, 3), np.float32), rna=np.ones((2, 3), np.float32))
    validate_artifact(path)
    with np.load(path, allow_pickle=False) as raw:
        assert require_state(raw, "wsi_identity") is not None
        assert require_state(raw, "full_patient") is None


def test_plip_cache_aligns_by_patient_id_and_rejects_wrong_model(tmp_path):
    cache = tmp_path / "teacher.npz"
    np.savez_compressed(cache, patient_ids=np.asarray(["P2", "P1"]), embeddings=np.eye(2, dtype=np.float32),
                        source_sha256=np.asarray(["b", "a"]), model_id=np.asarray("vinid/plip"))
    values, present, manifest = load_plip_teacher_cache(cache, ["P1", "P3", "P2"])
    assert present.tolist() == [True, False, True]
    assert values.shape == (3, 2) and manifest["n_aligned_patients"] == 2
    broken = tmp_path / "broken.npz"
    np.savez_compressed(broken, patient_ids=np.asarray(["P1"]), embeddings=np.ones((1, 2), np.float32),
                        source_sha256=np.asarray(["a"]), model_id=np.asarray("other/model"))
    with pytest.raises(ValueError, match="vinid/plip"):
        load_plip_teacher_cache(broken, ["P1"])
