# Lane 6 — Fast stress tests that fail before an expensive full run

## Queries/searches run

- Web: "effective rank singular value entropy representation collapse detection contrastive learning stopgrad"
- Web: "dimensional collapse self-supervised learning singular value spectrum diagnostic unit test"
- Code: read `morpheus/v2/losses.py`, `model.py`, `training.py`, `export.py`, `runner.py`, `contracts.py`, and the two fixture files `tests/test_v21_model.py`, `morpheus/tests/test_v2_recovery_contracts.py`.
- Grep: `assert_module_used`, `effective_rank`, `_trained_states_for_profile` (confirmed no effective-rank helper exists anywhere in the live tree).

## Sources

### Web (methods)
- Roy & Vetterli, "The effective rank: a measure of effective dimensionality" — effective rank = exp(Shannon entropy of L1-normalized singular values). Basis-/label-free, sensitive to *partial* collapse. (surfaced via the two searches below)
- Jing et al., "Understanding Dimensional Collapse in Contrastive Self-supervised Learning", arXiv:2110.09348 — collapse diagnosed by the singular-value spectrum of the embedding covariance; negative samples alone do NOT prevent dimensional collapse (directly relevant: our biology head has no paired contrastive by design).
- "A Cookbook of Self-Supervised Learning", arXiv:2304.12210 — singular-value spectrum / effective rank as the standard collapse diagnostic (a steep log-singular-value drop = collapse).
- WERank, arXiv:2402.09586; Orthogonality-reg, arXiv:2411.00392 — per-dimension variance floors do NOT prevent rank collapse; off-diagonal covariance structure is what collapses. Confirms our `variance_floor` (losses.py:29-34) is the wrong instrument.

### Code (file:line — the actual mechanisms these tests must trip on)
- `model.py:286` z_biology = normalize(biology_state); `model.py:283` biology_state = self.biology(biology) is the UNNORMALIZED head feeding programme NLL. Effective rank collapse observed here (~5-6/256).
- `model.py:260,280` z_identity path; `model.py:276-279` anchor gate = sigmoid(anchor_gate(identity)), correction = 0.25*tanh(scale)*gate*residual. Observed residual ~0 (gate/scale saturation).
- `model.py:290` programme_log_variance clamped to (-8, 8) — the ONLY existing guard; nothing tests it is load-bearing, and NLL mean can still be Inf if mean/target are non-finite.
- `losses.py:51-52` `.masked_fill(~mask, -1e4)` inside KL neighbourhood loss; `losses.py:63` masked_fill in supcon. `-1e4` + softmax + kl_div is the NaN-prone site.
- `losses.py:77` gaussian_nll: `0.5*(log_variance + (target-mean)^2 * exp(-log_variance))`. exp(-log_variance) overflows if log_variance not clamped; unbounded if a caller bypasses the model clamp.
- `losses.py:29-34` variance_floor is per-dimension std — cannot see rank.
- `training.py:238-239` already logs `identity_feature_std`, `biology_feature_std` (per-dim mean std) — a hook exists but it is a scalar, blind to rank.
- `training.py:194-196` identity InfoNCE (only z_identity trained contrastively); `training.py:212-222` per-view programme loss.
- `export.py:54-55` widths dict: only wsi/rna/full identity+biology+patient (+semantic). z_context, z_uncertainty, 5x residual NEVER exported → dead heads.
- `contracts.py:76-85` `assert_module_used(module, loss)` helper EXISTS but is called by no shipped test in the live tree (only a snapshot copy). Free to reuse.
- Fixtures: `_batch()` at `tests/test_v21_model.py:19` (+`_config()` at :12, hidden_dim=16, programme_dim passed as 3); `_model_and_batch()` at `morpheus/tests/test_v2_recovery_contracts.py:20` (anchor=True, hidden_dim=8, programme_dim=2).

## Findings

The 31 existing CPU tests check gradient *existence* per head and `isfinite` on the final artifact, but nothing checks the *mechanism* failures we actually hit: rank collapse, NaN inside KL/NLL, log_variance mode collapse, dead losses, gate saturation, or slot accounting. All six can be caught in seconds on the tiny fixtures by running a handful of `_batch()` steps and inspecting the resulting tensors — no GPU, no real data.

Key design point grounded in the literature: **effective rank (Roy-Vetterli) on the embedding matrix is the correct, label-free, partial-collapse-sensitive diagnostic**, and per-dimension `variance_floor` provably cannot catch it (WERank / Jing 2110.09348). So test (a) must compute the singular-value-entropy effective rank, not std.

One realism caveat: the tiny fixtures are batch=3, hidden=16, so absolute effective-rank *floors of 256* cannot be asserted on them. The tests instead assert **relative / mechanistic** properties that hold at any scale (rank grows toward min(batch, dim) on random-separable synthetic data after N steps; biology rank not collapsing to 1; gate variance > 0). This is the smallest change that would still have tripped: our failure was rank ~5-6 while identity ~84 — the *ratio* and the *trend-toward-collapse* are visible even at small scale on a synthetic separable batch.

Recommend one new file `morpheus/v2/tests/test_stress_collapse.py` with six tests + two small shared helpers, all reusing `_batch()`/`_model_and_batch()` and `V2Trainer.step`. Runtime target < 5 s total.

## Recommended change (file:line, exact)

New file: `C:/Users/mobar/OneDrive/biorag/morpheus/v2/tests/test_stress_collapse.py`.

Two helpers first:

```python
import torch, torch.nn.functional as F
from morpheus.v2.model import TumorStateV2
from morpheus.v2.training import V2LossSchedule, V2Trainer
from morpheus.v2.tests.test_v21_model import _config, _batch  # reuse fixtures

def effective_rank(x: torch.Tensor, center: bool = True) -> float:
    # Roy & Vetterli: exp(entropy of L1-normalized singular values).
    x = x - x.mean(0, keepdim=True) if center else x
    s = torch.linalg.svdvals(x.float())
    p = (s / s.sum().clamp_min(1e-12)).clamp_min(1e-12)
    return float(torch.exp(-(p * p.log()).sum()))

def _train_n(model, batch, steps=20, profile="full"):
    opt = torch.optim.AdamW(model.parameters(), lr=1e-2)
    tr = V2Trainer(model, opt, V2LossSchedule(objective_profile=profile, warmup_epochs=0),
                   "cpu", gradient_diagnostics_every=0)
    for _ in range(steps):
        opt.zero_grad(set_to_none=True)
        loss, m, out = tr.step(batch, epoch=1)
        loss.backward(); opt.step()
    return m, out
```

### (a) Effective-rank / singular-value floor — fixture `_batch()`
`test_biology_effective_rank_does_not_collapse_below_identity`
- Build a **synthetic separable** batch: replace `_batch()` patches/rna with `batch_size >= 8` rows drawn as `k` well-separated clusters (so a healthy head *can* reach rank ~min(batch,256)). Train 20 steps, profile `"full"`.
- Assert `effective_rank(out["z_biology"]) >= 2.0` (not mode-collapsed to a line) AND `effective_rank(out["z_biology"]) >= 0.4 * effective_rank(out["z_identity"])`.
- Rationale: our real failure was biology rank / identity rank ≈ 5.5/84 ≈ 0.065 « 0.4. This ratio gate trips on exactly that pathology and is scale-free. Uses `svdvals` (Roy-Vetterli), the literature-standard partial-collapse detector.

### (b) NaN/Inf guards around masked_fill(-1e4) softmax (KL) and log-variance (NLL) — fixture `_batch()`
`test_kl_and_nll_are_finite_under_degenerate_inputs`
- Directly unit-test `losses.programme_neighbourhood_loss` and `losses.gaussian_nll` (not the model):
  - Neighbourhood: pass `state` and `targets` where one row is all-zeros (so `_norm` divides by ~0) and batch=3 (min size). Assert `torch.isfinite(programme_neighbourhood_loss(state, targets))`.
  - NLL: call `gaussian_nll(mean, log_variance, target)` with `log_variance = torch.full(..., 40.0)` (i.e. an UNclamped value that bypasses model.py:290) and with a NaN target row masked out. Assert result finite; assert that WITHOUT the `clamp(-8,8)` the raw `exp(-log_variance)` path would overflow (documents why the clamp is load-bearing).
- Also an integration check: `test_full_step_metrics_all_finite` — run 1 `V2Trainer.step` on `_batch()` and assert every value in `metrics` is finite (guards the `-1e4`→softmax→kl_div chain end to end).

### (c) Mode collapse: log_variance → +inf / mean → constant — fixture `_batch()`
`test_programme_log_variance_not_saturated_and_mean_not_constant`
- Train 20 steps profile `"programme_only"`. Then:
  - Assert `out["programme_log_variance"].max() < 7.9` (not pinned at the +8 clamp ceiling → the model is not "explaining away" everything as noise).
  - Assert `out["programme_mean"].std(0).mean() > 1e-3` (mean head not collapsed to a constant across the batch).
- Rationale: Gaussian-NLL's failure mode is driving log_variance up to make the squared-error term vanish; the clamp at +8 makes saturation *visible* as a pile-up at the ceiling, which this asserts against.

### (d) No-op / dead-loss + profile-toggles-right-heads — fixtures `_batch()` and `_model_and_batch()`
`test_every_declared_active_head_receives_gradient` (reuses `contracts.assert_module_used`, currently untested)
- For each profile in `("identity_only","programme_only","full")`: run `V2Trainer.step`, `loss.backward()`, then for each head the profile *declares* trained (map via `_trained_states_for_profile`) call `assert_module_used(head_module, loss)`; for each head it declares *inactive*, assert grad is None/zero (mirrors the existing but stricter than `test_profile_declarations_match_loss_gradient_families`).
`test_dead_heads_are_not_silently_exported`
- Assert the model computes `z_context`, `z_uncertainty`, `z_wsi_residual` (etc.) but that `export.widths` (import the dict literal or the state list from `_trained_states_for_profile("full")`) contains NONE of them — i.e. the DEAD wiring is *intentional and documented*, so a future accidental export of an untrained head fails the test. This freezes the contract rather than the bug.

### (e) Anchor-gate saturation — fixture `_model_and_batch()` (anchor=True) or `_batch(anchor=True)`
`test_anchor_gate_not_saturated_and_residual_can_move`
- Model with `use_mlp_clip_anchor=True`. Train 20 steps profile `"identity_only"` on a separable batch. Then:
  - Assert `out["anchor_gate_mean"]` is in `(0.02, 0.98)` (gate not stuck at ~0 or ~1) — computed at model.py:298.
  - Assert `out["anchor_correction_norm"] > 1e-4` after training (residual is non-degenerate; at init it is exactly 0 by design, so this specifically detects "residual never learns to move", our observed ~0 residual).
  - Assert `model.anchor_residual_scale` moved off 0 OR gate variance `> 1e-4` across the batch (the correction is input-dependent, not a dead constant).

### (f) Slot-consumption correctness — fixture `_batch()`
`test_take_consumes_every_query_slot_exactly_once`
- White-box the `take()` cursor in `model.py:252-259`: after a forward, assert `identity_slots + biology_slots + context_slots + 5*residual_slots + uncertainty_slots == config.query_slots` (matches the `query_slots` property, model.py:42-44) AND monkeypatch/instrument `take` to record the final cursor == `query.shape[1]` (every one of the 32 slots consumed, none double-counted, none dropped). Also assert the residual dict has exactly 5 keys `{wsi,rna,clinical,snv,cnv}` (model.py:258) so a slot-count edit that desyncs the header can't pass silently.

## Risks & scaling

- **Small-fixture rank ceilings**: batch=3, hidden=16 fixtures cannot exhibit rank-256 collapse. Mitigation baked into (a)/(c)/(e): use a locally-built `batch_size>=8`, k-cluster **separable** synthetic batch and assert *relative* rank (biology/identity ratio) and *trend* (rank > 2), which are the scale-free signatures of the real failure. Do NOT assert absolute floors like `rank >= 100`.
- **Step count vs runtime**: 20 AdamW steps × 3 forward views × 3 profiles is still < 5 s on CPU at hidden=16. Keep `gradient_diagnostics_every=0` (the extra autograd pass in `_gradient_conflict_metrics` is the slow part) — already handled in the helper.
- **Flakiness**: seed `torch.manual_seed` at the top of each test; the ratio/relative thresholds (0.4, 2.0) are chosen with margin against the observed 0.065 so a healthy small model passes and the pathological one fails. Tune once against a known-good and known-bad checkpoint if available.
- **False confidence**: these are *mechanism* tests, not performance tests. They catch collapse/NaN/dead-loss/gate-saturation but not "molecular prompting is only +0.07" — that needs the held-out protocol and is out of lane. Document that (a) is a *necessary not sufficient* gate.
- **Guard-coupling**: (b) deliberately tests `gaussian_nll` with an unclamped log_variance to prove the model.py:290 `clamp(-8,8)` is load-bearing; if a future refactor moves the clamp, this test tells you the loss is now unguarded. Keep the loss-level and model-level guards both, or the test will (correctly) fail.
