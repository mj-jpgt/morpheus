# Lane 1 — Anti-collapse remedy for the L2-normalized 256-D biology head

## Queries/searches run

Web (WebSearch/WebFetch):
- "VICReg variance covariance regularization prevent dimensional collapse L2 normalized embeddings"
- "effective rank dimensional collapse contrastive learning covariance regularizer weight"
- "VICReg covariance term hurts on L2 normalized features interaction normalization dimensional collapse contrastive"
- "Barlow Twins vs VICReg covariance off-diagonal loss weight sensitivity small batch 256 dimensions"
- WebFetch of the VICReg paper PDF (arxiv 2105.04906v2) for the exact covariance term formula c(Z) and default coefficients.

Code (Read/Grep):
- Read `morpheus/v2/losses.py` (all loss primitives).
- Read `morpheus/v2/training.py` (full training engine, `step`, loss wiring, weights).
- Read `morpheus/v2/model.py:1-60` (config) and `:240-300` (Query Former take() + head outputs + normalization).
- Grep `variance_floor|whitened_cross_covariance|centered_cross_covariance` across v2 (call sites).
- Grep `effective_rank|matrix_rank|svdvals` across v2 (none in code — only in research/*.md).
- Grep tests for `variance_floor|separation|z_biology|feature_std`; Glob `tests/*.py`; inspected `_batch()` fixture in `test_v21_model.py:19`.

## Sources

Web:
- VICReg (Bardes, Ponce, LeCun 2021), arxiv 2105.04906v2 — covariance term c(Z) = (1/d)·Σ_{i≠j} C(Z)²_ij; variance term v(Z)=max(0, γ−√(var+eps)), γ=1; default coeffs invariance 25 / variance 25 / covariance 1. https://arxiv.org/pdf/2105.04906v2
- "Understanding Dimensional Collapse in Contrastive SSL" / DirectCLR (Jing, Vincent, LeCun, Tian) — cosine/InfoNCE on normalized features drives collapse to a low-rank subspace; effective rank as the collapse metric. https://openreview.net/pdf?id=YevsQ05DEN7 , https://ai.meta.com/blog/understanding-dimensional-collapse/
- Barlow Twins (Zbontar et al. 2021) — off-diagonal cross-correlation penalty; robust at batch 256. https://proceedings.mlr.press/v139/zbontar21a/zbontar21a.pdf
- Ablation consensus (VICReg paper + eval studies): removing the **variance** term collapses to a point (catastrophic); removing the **covariance** (decorrelation) term costs ~1-2% and is precisely the term that governs *rank* / dimensional usage. Effective rank = exp(entropy of normalized covariance eigenvalues).

Code (file:line):
- `model.py:286` — `z_biology = normalize(biology_state, dim=-1)` (the head is L2-normalized).
- `model.py:283,290` — `biology_state = self.biology(biology)`; programme NLL heads consume the UNNORMALIZED `biology_state`.
- `losses.py:29-34` — `variance_floor` = `relu(target_std − per-dim std).mean()`: PER-DIMENSION only, no cross-dimension term.
- `losses.py:21-26` — `centered_cross_covariance`: mean of squared entries of a CROSS covariance (identity×biology). This is exactly the VICReg off-diagonal machinery, but cross-view, not within-biology.
- `losses.py:46-53` — `programme_neighbourhood_loss` (KL, temp 0.20) operates on `_norm(state)` cosine logits.
- `losses.py:56-66` — `supervised_programme_contrastive` (temp 0.20) operates on `_norm(state)` cosine logits.
- `training.py:232-237` — variance term = `variance_floor(z_identity)+variance_floor(z_biology)`, weight 0.01 (`variance_after_warmup`, training.py:32).
- `training.py:226-231` — separation = `whitened_cross_covariance(z_identity, z_biology)`, weight 0.01.
- `training.py:238-239` — only diagnostic logged for biology is `biology_feature_std` (per-dim std mean) — no rank metric.

## Findings

1. **Mechanism of the collapse.** The biology head is L2-normalized (model.py:286) and every biology objective that shapes its *geometry* is a cosine-similarity objective on that normalized vector: neighbourhood-KL (losses.py:51-52) and supcon (losses.py:63) both build logits from `_norm(z_biology) @ _norm(z_biology).T`. This is exactly the regime the dimensional-collapse literature identifies (DirectCLR / "Understanding Dimensional Collapse"): similarity-only objectives on unit-norm features are satisfiable inside a low-dimensional subspace, so the covariance matrix collapses to a few large eigenvalues → observed effective rank ~5-6 of 256. The Gaussian-NLL (weight 1.0) trains only `programme_mean/log_variance` linear read-outs off the *unnormalized* `biology_state` (model.py:290) and does not constrain the normalized 256-D geometry, so it does not restore rank.

2. **Why the existing variance floor does not fix it (confirmed, not just asserted).** `variance_floor` (losses.py:29-34) is `relu(1 − per-dimension std)`. It is diagonal: it can be fully satisfied while all mass lives in a rank-5 subspace as long as each *coordinate* has ≥1.0 std. It has no off-diagonal / decorrelation term, so it cannot raise rank. This matches the VICReg ablation finding that the **variance** term prevents *point* collapse but the **covariance** term is the one that prevents *dimensional* collapse. So the missing ingredient is precisely a covariance-decorrelation term on the biology features. Its weight is also tiny (0.01, training.py:32) and off during warmup (training.py:53).

3. **Method comparison for THIS head (L2-normalized, 256-D, small batch, partly KL/contrastive).**
   - **VICReg covariance term** — c(Z)=(1/d)Σ_{i≠j}C(Z)²_ij. Directly penalizes off-diagonal covariance → pushes the covariance toward diagonal → raises effective rank. O(B·d²) = 3·256² per view, negligible. **Best fit.** The codebase already contains the identical primitive (`centered_cross_covariance`), so the change is a 6-line variant + one call.
   - **VICReg variance term** — already present as `variance_floor`; keep it, it guards against point collapse, but it is not the rank fix.
   - **Barlow Twins** — needs a *paired second view* to form a cross-correlation matrix. Biology has "no paired contrastive by design" (there is no biology positive-pair view), so Barlow's cross-correlation is not naturally defined here. Would require inventing an augmentation pair. Rejected: larger change, worse fit.
   - **Whitening / W-MSE / DirectCLR** — DirectCLR removes the projector and back-props only a sub-vector; it is an *architecture* change to the whole head. Whitening needs a stable batch covariance inverse (Cholesky) — brittle at batch 3 in tests / small real batches. Rejected: heavier, less robust.
   - **Spectral / effective-rank penalty (log-det, nuclear)** — optimizes the target metric directly but needs an eigendecomposition every step (unstable gradients near degenerate spectra, costly). Overkill vs. the covariance term, which is a cheap surrogate that provably diagonalizes the covariance. Rejected as the primary lever; useful only as the *verification metric* (see below).
   - **EMA / stop-grad (BYOL/SimSiam style)** — prevents collapse via asymmetry, not a regularizer you add to an existing symmetric head; would restructure training. Rejected.

4. **Interaction with L2-normalization (the specific concern in the remit).** VICReg is designed on *unnormalized* features and its variance term substitutes for normalization. Here the head is already L2-normalized, so: (a) do NOT apply VICReg's variance term to the normalized vector as-is — on a unit vector, per-coordinate variance is bounded and the γ=1 target is unreachable (that is why the code sensibly uses `target_std=1.0` as a *floor* with relu, and applies it broadly). (b) The **covariance** term is well-behaved on normalized features: it only asks the batch covariance to be diagonal, which is achievable on the unit sphere (a well-spread set of unit vectors has near-diagonal covariance). Apply the covariance term to the **normalized** `z_biology` so it acts on the exact geometry the KL/supcon objectives shape. This is the clean, mechanism-matched choice.

5. **Weight.** VICReg's paper ratio is cov:inv = 1:25 (≈0.04) at d=... with an unnormalized invariance MSE, not transferable verbatim. For a normalized head where the competing "invariance"-like terms (neighbourhood 0.20, supcon 0.20) are already small, and to stay conservative alongside grad-clipping (training.py:318), start the covariance weight at **0.04** (4× the current variance floor weight of 0.01, same order as separation 0.01, and it should measurably dominate the tiny decorrelation pressure currently absent). It is a floor-style penalty (0 once decorrelated) so it will not fight the primary NLL/identity objectives once rank recovers. Keep the existing variance floor at 0.01 to retain point-collapse protection. Turn the new term on after warmup like the other structural terms (training.py:53).

## Recommended change (file:line, exact)

Single mechanism, smallest surgical change: add a **within-feature covariance decorrelation term on the normalized biology state** and give it a weight.

1. `losses.py` — add a within-feature covariance primitive (mirrors `centered_cross_covariance` at losses.py:21-26, off-diagonal only, VICReg c(Z)):

```python
def feature_decorrelation(state: torch.Tensor, eps: float = 1e-4) -> torch.Tensor:
    """VICReg covariance term: push the batch covariance toward diagonal.

    Per-dimension variance_floor cannot raise rank; this off-diagonal penalty
    is the term that prevents dimensional (rank) collapse of z_biology.
    """
    if len(state) < 2:
        return state.new_zeros(())
    centered = state - state.mean(dim=0, keepdim=True)
    cov = (centered.T @ centered) / max(len(state) - 1, 1)
    off_diagonal = cov - torch.diag(torch.diag(cov))
    return off_diagonal.square().sum() / state.shape[1]
```

2. `training.py:15` — import it alongside `variance_floor`.

3. `training.py:32` (V2LossSchedule) — add field `decorrelation_after_warmup: float = 0.04` and include `"decorrelation": 0.0 if warmup else self.decorrelation_after_warmup` in `weights()` (training.py:45-58); include it in the `full` and `programme_only` gate paths (training.py:63-65) so it is active whenever biology geometry is trained.

4. `training.py:232-237` — right after the variance-floor block, add (on the **normalized** state, model.py:286 already normalizes `output["z_biology"]`):

```python
if weights["decorrelation"]:
    decorrelation = feature_decorrelation(output["z_biology"])
    loss = loss + weights["decorrelation"] * decorrelation
    metrics["biology_decorrelation"] = float(decorrelation.detach())
```

Do NOT apply it to `z_identity` (rank ~84, healthy; adding pressure risks fighting the anchored InfoNCE).

## Verify rank recovery

Add an effective-rank diagnostic + one collapse-guard test (no rank metric exists anywhere in code today — Grep for `effective_rank|matrix_rank|svdvals` hit only research/*.md).

- **Metric to log** (in `training.py` diagnostics, evaluate-time on `z_biology` over an epoch's accumulated states, or a batch-level proxy): effective rank = exp(entropy of normalized covariance eigenvalues):
  ```python
  s = torch.linalg.svdvals(z_biology - z_biology.mean(0))
  p = (s / s.sum()).clamp_min(1e-12)
  eff_rank = float(torch.exp(-(p * p.log()).sum()))
  ```
  Log `biology_effective_rank` next to `biology_feature_std` (training.py:238-239).
- **Success criterion:** biology effective rank climbs from ~5-6 toward the identity head's regime (target ≥ ~30-50 of 256; identity sits at ~84). Watch that neighbourhood-KL and supcon (losses.py) do not degrade — decorrelation should be complementary, not antagonistic (log both).
- **New CPU test** (in `tests/test_v21_model.py`, reuse `_batch()` at :19): construct a deliberately rank-deficient `z_biology` (e.g. a batch spanning 2 directions) and assert `feature_decorrelation` is > a small threshold and its gradient is nonzero; assert a well-spread random unit batch gives a near-zero value. This is the missing "catches rank collapse" test the remit flags as absent.

## Risks & scaling

- **Batch-size dependence of covariance.** The covariance estimate is noisy at small batch (tests run batch 3; the term correctly returns 0 for len<2). On the A10 real batches this is fine; keep the weight modest (0.04) so noisy off-diagonal estimates do not inject high-variance gradients. Grad-clip at training.py:318 already caps blow-ups.
- **AMP / bf16.** `training.py:313` autocasts to bf16. `svdvals` and the covariance matmul are numerically touchier in bf16; compute the diagnostic in fp32 (`.float()`) — the loss term itself (matmul + square) is stable in bf16 but cast to fp32 for the eff-rank log.
- **Interaction with separation term** (whitened_cross_covariance, training.py:228, weight 0.01): both act on covariance structure but on different pairs (identity×biology vs biology×biology). No conflict; if anything decorrelating biology makes the cross-covariance easier to zero.
- **Over-decorrelation.** Too large a weight forces an isotropic covariance that can wash out genuine low-rank programme structure (the 50-D Hallmark target is itself ~50-D). 0.04 is a floor-style penalty (→0 once diagonal) so it will not over-flatten; if eff-rank overshoots far past ~50 with degraded neighbourhood-KL, lower to 0.02. Do not exceed ~0.1.
- **Scope discipline.** This is the minimal mechanism fix (missing decorrelation term). It does NOT touch the dead residual/context wiring, the detached memory bank, or head widths — out of lane and unnecessary for rank recovery.
