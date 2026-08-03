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
from morpheus.v2.training import PairedBiologyMemoryBank, V2LossSchedule, V2Trainer


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


# --- REGRESSION: the padded-head structure gate ------------------------------

def test_structure_losses_survive_a_padded_programme_head():
    """--programme-head-dim (default 256) pads the 50-wide Hallmark target with NaN.
    A plain target_mask.all(dim=1) is then False for EVERY row, which silently
    disabled neighbour-KL and supcon in every programme_only run -- and those two
    losses ARE the collapse mechanism D1 exists to test. It then aborted at epoch 40
    on the G2.2 liveness check, i.e. after the GPU time was already spent."""
    import numpy as np, torch
    for head_dim, expected in ((50, 6), (128, 6), (256, 6)):
        targets = np.random.default_rng(0).normal(size=(6, 50)).astype(np.float32)
        padded = np.full((6, head_dim), np.nan, dtype=np.float32)
        padded[:, :50] = targets
        mask = torch.from_numpy(np.isfinite(padded))
        real_axis = mask.any(dim=0)
        gate = mask[:, real_axis].all(dim=1) if bool(real_axis.any()) else torch.zeros(6, dtype=torch.bool)
        assert int(gate.sum()) == expected, f"head_dim={head_dim} gated {int(gate.sum())}/6 rows"
        if head_dim > 50:
            assert int(mask.all(dim=1).sum()) == 0, "fixture must reproduce the original defect"


def test_a_row_with_a_genuinely_missing_axis_is_still_excluded():
    """The fix must not become permissive: a patient missing a REAL Hallmark axis
    must still be excluded from the structure losses."""
    import numpy as np, torch
    padded = np.full((4, 256), np.nan, dtype=np.float32)
    padded[:, :50] = np.random.default_rng(1).normal(size=(4, 50))
    padded[2, 7] = np.nan                      # one real axis genuinely missing
    mask = torch.from_numpy(np.isfinite(padded))
    real_axis = mask.any(dim=0)
    gate = mask[:, real_axis].all(dim=1)
    assert gate.tolist() == [True, True, False, True]


def _replay_queue_alignment(*, freeze: bool, steps: int = 12) -> float:
    """Replay ONE fixed batch; report how close the queue keys sit to that batch.

    The defect this measures: with a live queue, replaying one batch overwrites
    every slot with re-encoded copies of that batch's own states, so the InfoNCE
    "negatives" become the queries. Returns max cosine similarity between any
    queue key and any current batch state.
    """
    torch.manual_seed(11)
    model = TumorStateV2(_config(hidden=32), programme_dim=8)
    optimiser = torch.optim.AdamW(model.parameters(), lr=1e-3)
    trainer = V2Trainer(model, optimiser,
                        V2LossSchedule(objective_profile="programme_free", warmup_epochs=0), "cpu")
    trainer.biology_memory = PairedBiologyMemoryBank(capacity=16)
    trainer.prime_biology_memory([_big_batch(n=2, k=8, seed=seed) for seed in range(8)],
                                 minimum_unique_keys=9)
    trainer.freeze_biology_memory = freeze
    fixed = _big_batch(n=4, k=8, seed=99)
    for _ in range(steps):
        optimiser.zero_grad(set_to_none=True)
        loss, _, _ = trainer.step(fixed, epoch=1)
        loss.backward()
        optimiser.step()
    with torch.no_grad():
        current = torch.nn.functional.normalize(
            trainer.model(fixed, view="wsi")["z_biology"].float(), dim=-1)
        keys, _, _ = trainer.biology_memory.view()
        return float((keys @ current.T).max())


def test_biology_memory_refreshes_during_normal_training() -> None:
    """The freeze flag must not change the training path: default is live."""
    torch.manual_seed(12)
    model = TumorStateV2(_config(hidden=32), programme_dim=8)
    trainer = V2Trainer(model, torch.optim.AdamW(model.parameters(), lr=1e-3),
                        V2LossSchedule(objective_profile="programme_free", warmup_epochs=0), "cpu")
    assert trainer.freeze_biology_memory is False
    trainer.prime_biology_memory([_big_batch(n=2, k=8, seed=seed) for seed in range(6)],
                                 minimum_unique_keys=9)
    before = trainer.biology_memory.size
    trainer.step(_big_batch(n=4, k=8, seed=77), epoch=1)
    assert trainer.biology_memory.size == before + 4


def test_frozen_queue_does_not_absorb_the_replayed_batch() -> None:
    torch.manual_seed(13)
    model = TumorStateV2(_config(hidden=32), programme_dim=8)
    trainer = V2Trainer(model, torch.optim.AdamW(model.parameters(), lr=1e-3),
                        V2LossSchedule(objective_profile="programme_free", warmup_epochs=0), "cpu")
    trainer.prime_biology_memory([_big_batch(n=2, k=8, seed=seed) for seed in range(6)],
                                 minimum_unique_keys=9)
    trainer.freeze_biology_memory = True
    keys_before = trainer.biology_memory.indices[:trainer.biology_memory.size].clone()
    states_before = trainer.biology_memory.wsi_states[:trainer.biology_memory.size].clone()
    for _ in range(5):
        trainer.step(_big_batch(n=4, k=8, seed=77), epoch=1)
    size = trainer.biology_memory.size
    assert torch.equal(trainer.biology_memory.indices[:size], keys_before)
    assert torch.equal(trainer.biology_memory.wsi_states[:size], states_before)


def test_a_live_queue_turns_the_replayed_batch_into_its_own_negatives() -> None:
    """The mechanism, stated as a measurement rather than as an outcome.

    Whether this is *sufficient* to pin the contrastive term at chance depends on
    scale and on how collinear the WSI states are; it is not reproduced by a toy
    model with random features, so this test deliberately asserts the mechanism
    and NOT a loss improvement. The loss claim belongs to a real G2.6 run.
    """
    live = _replay_queue_alignment(freeze=False)
    frozen = _replay_queue_alignment(freeze=True)
    # The absolute number is not the claim and is trajectory-dependent: the key
    # is written during step N and probed after the optimiser has moved again,
    # so it drifts by however far one step travels.  It was 0.9924 against the
    # uncentred objective and is 0.9855 against the centred one (2026-08-03) for
    # exactly that reason.  The SEPARATION between live and frozen is the
    # mechanism, so assert that too rather than lean on a tight absolute bound.
    assert live > 0.95, f"a live queue should hold the batch's own states, got {live:.4f}"
    assert frozen < live, f"frozen keys must stay distinct from the queries: {frozen:.4f} vs {live:.4f}"
    assert live - frozen > 0.10, f"freezing must visibly separate keys from queries: {live:.4f} vs {frozen:.4f}"


def test_decorrelation_is_minimised_by_the_collapse_it_claims_to_prevent() -> None:
    """Documents the hazard that broke programme_free (2026-08-02).

    `feature_decorrelation` standardises before penalising off-diagonal
    correlation, so an all-rows-identical batch standardises to zero and the
    penalty vanishes. Total collapse is its global MINIMUM. This is not a bug to
    fix in isolation -- it is the VICReg contract -- but it means the term must
    never ship without a variance floor beside it.
    """
    from morpheus.v2.losses import feature_decorrelation, variance_floor
    torch.manual_seed(0)
    healthy = torch.randn(16, 512)
    collapsed = torch.randn(1, 512).repeat(16, 1)
    assert float(feature_decorrelation(collapsed)) < 1e-9 < float(feature_decorrelation(healthy))
    # The variance floor is the counter-force: maximal exactly where
    # decorrelation is minimal.
    target = 512 ** -0.5
    assert (float(variance_floor(collapsed, target_std=target))
            > float(variance_floor(healthy, target_std=target)))


def test_both_d1_arms_carry_identical_regularisation() -> None:
    """The D1 contrast must measure programme supervision and nothing else.

    SYMMETRY is the invariant, not any particular weight. An arm carrying
    different regularisation from the other makes the comparison measure more
    than the objective under test, whatever the weights happen to be.

    This test previously also asserted `decorrelation > 0` in both arms, on the
    theory that decorrelation was the anti-collapse force. Measurement falsified
    that on 2026-08-03, in the opposite direction: `feature_decorrelation` has
    total collapse as its global minimum (pinned separately below), it collapses
    the representation at every weight from 0.001 to 4.0 while switching itself
    off, and a per-dimension variance floor provably cannot stop it because the
    rank-1 family `z_i = m + a_i*u` satisfies such a floor. Requiring the term to
    be present is therefore not a safety property, and the assertion is removed.

    What survives is the standing hazard, kept as a CONDITIONAL: decorrelation
    must never appear without a variance floor beside it.
    """
    free = V2LossSchedule(objective_profile="programme_free", warmup_epochs=0).weights(1)
    only = V2LossSchedule(objective_profile="programme_only", warmup_epochs=0).weights(1)
    regularisers = ("decorrelation", "variance", "separation")
    for name in regularisers:
        assert free[name] == only[name], (
            f"D1 arms disagree on {name}: programme_free={free[name]} programme_only={only[name]}; "
            "the contrast would measure regularisation as well as supervision")
    for name, weights in (("programme_free", free), ("programme_only", only)):
        if weights["decorrelation"] > 0:
            assert weights["variance"] > 0, (
                f"{name} carries decorrelation without a variance floor; the covariance penalty is "
                "minimised by the collapse it claims to prevent")


def test_paired_infonce_removes_the_batch_common_direction_by_default() -> None:
    """Centring is load-bearing: without it the D1 objective starts ABOVE chance.

    `z_biology` is L2-normalised but never centred.  Measured on the real cohort
    at initialisation, 81% of its squared norm is one direction shared by every
    patient, the positives hold no advantage over the negatives (minimum margin
    -0.219) and the loss starts at 3.0762 against chance ln(16)=2.7726.  From
    above chance, erasing every distinction is a DESCENT direction, and the
    biology head reaches effective rank 1.00 within 50 steps.

    This reproduces the pathology in miniature: a batch whose patient-specific
    signal is buried under a large shared direction is at chance uncentred and
    solvable centred.
    """
    torch.manual_seed(11)
    dim, n = 16, 12
    common = torch.nn.functional.normalize(torch.randn(1, dim), dim=-1)
    signal = torch.nn.functional.normalize(torch.randn(n, dim), dim=-1)
    # 97% shared direction, 3% patient identity -- the geometry measured on the
    # real 16-patient batch, where the shared component is ~0.81 of the norm.
    wsi = 0.97 * common + 0.03 * signal
    rna = 0.97 * common + 0.03 * signal
    ids = torch.arange(n)
    queue_wsi = 0.97 * common + 0.03 * torch.nn.functional.normalize(torch.randn(40, dim), dim=-1)
    queue_rna = 0.97 * common + 0.03 * torch.nn.functional.normalize(torch.randn(40, dim), dim=-1)
    queue_ids = torch.arange(100, 140)
    raw, _ = paired_infonce_with_memory(wsi, rna, ids, queue_wsi, queue_rna, queue_ids, centre=False)
    centred, _ = paired_infonce_with_memory(wsi, rna, ids, queue_wsi, queue_rna, queue_ids)
    # Uncentred the shared direction swamps the identity signal; centred, the
    # identical WSI/RNA pairing is trivially recoverable.
    assert float(centred) < float(raw)
    assert float(centred) < 0.10 < float(raw)


def test_centring_never_zeroes_a_real_ragged_batch() -> None:
    """Real D1 batches hold B=1-3 patients; centring those by their own mean
    would map B=1 to the zero vector and B=2 to an antipodal pair.  Below
    `min_negatives` the estimate must come from the queue instead."""
    from morpheus.v2.losses import population_offset
    torch.manual_seed(12)
    for batch_size in (1, 2, 3):
        current = torch.randn(batch_size, 8)
        memory = torch.randn(64, 8)
        offset = population_offset(current, memory, min_batch=8)
        assert offset is not None and offset.shape == (1, 8)
        assert not torch.allclose(current - offset, torch.zeros_like(current))
        # With no queue to borrow from, refuse to centre rather than annihilate.
        assert population_offset(current, None, min_batch=8) is None
    assert population_offset(torch.randn(16, 8), None, min_batch=8) is not None


def test_g26_grades_the_uncentred_contrastive_number() -> None:
    """Centring changes the value of the quantity G2.6 thresholds at <= 0.10.

    The graded metric therefore keeps its historical definition -- uncentred,
    queue included -- so that a fix to the optimisation cannot be mistaken for
    a relaxed criterion.  The centred value is reported separately.
    """
    torch.manual_seed(13)
    # Dropout off so the metric can be reproduced exactly by a second forward.
    config = V2ModelConfig(patch_dim=8, rna_dim=4, hidden_dim=32, heads=4, layers=1,
                           local_slots=4, slide_slots=2, patient_slots=2, dropout=0.0)
    model = TumorStateV2(config, programme_dim=8)
    trainer = V2Trainer(model, torch.optim.AdamW(model.parameters(), lr=1e-3),
                        V2LossSchedule(objective_profile="programme_free", warmup_epochs=0), "cpu")
    batches = [_big_batch(n=2, k=8, seed=seed) for seed in range(6)]
    trainer.prime_biology_memory(batches, minimum_unique_keys=9)
    # Frozen so the queue the metric saw is the queue this test recomputes with;
    # a live queue absorbs the batch at the end of the step.
    trainer.freeze_biology_memory = True
    batch = _big_batch(n=2, k=8, seed=99)
    _, metrics, _ = trainer.step(batch, epoch=1)
    assert "biology_contrastive" in metrics and "biology_contrastive_centred" in metrics
    memory_wsi, memory_rna, memory_indices = trainer.biology_memory.view()
    with torch.no_grad():
        out_wsi, out_rna = model(batch, view="wsi"), model(batch, view="rna")
        expected, _ = paired_infonce_with_memory(
            out_wsi["z_biology"], out_rna["z_biology"], batch["indices"],
            memory_wsi, memory_rna, memory_indices, centre=False)
        centred, _ = paired_infonce_with_memory(
            out_wsi["z_biology"], out_rna["z_biology"], batch["indices"],
            memory_wsi, memory_rna, memory_indices)
    assert metrics["biology_contrastive"] == pytest.approx(float(expected), abs=1e-5)
    assert metrics["biology_contrastive_centred"] == pytest.approx(float(centred), abs=1e-5)
    assert metrics["biology_contrastive"] != pytest.approx(metrics["biology_contrastive_centred"], abs=1e-6)
