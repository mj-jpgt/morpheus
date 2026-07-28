# Lane 3: SigLIP sigmoid loss swap for the identity contrastive

Remit: replace the symmetric InfoNCE identity alignment at `training.py:194-196`
with the SigLIP (Zhai et al., ICCV 2023) sigmoid pairwise loss using a learnable
temperature `t` and bias `b`, grounded in the actual MORPHEUS V2 code.

## Queries/searches run

- WebSearch: "SigLIP Zhai 2023 sigmoid loss language image pretraining logit bias temperature initialization"
- WebSearch: "SigLIP sigmoid loss implementation logit_scale log 10 bias -10 label -1 pytorch small batch 4096 512"
- WebFetch arxiv 2303.15343 (PDF, binary — unreadable, fell back to secondary sources)
- WebFetch emergentmind SigLIP topic page (loss equation, decoupling claims)
- WebFetch github big_vision issue #130 (canonical JAX pseudocode + init values)
- WebFetch arxiv 2410.17427 SigCLR (small-batch sigmoid-vs-softmax evidence)
- Read: `morpheus/v2/losses.py` (full), `morpheus/v2/training.py:1-80, 120-160, 160-230`
- Read: `morpheus/v2/model.py:195-234`
- Grep: identity heads / nn.Parameter / temperature in `model.py`; optimizer/parameters/clip in `training.py`; test fixtures in `tests/`

## Sources

Web:
- Zhai, Mustafa, Kolesnikov, Beyer. "Sigmoid Loss for Language Image Pre-Training." ICCV 2023. https://arxiv.org/abs/2303.15343 , CVF PDF https://openaccess.thecvf.com/content/ICCV2023/papers/Zhai_Sigmoid_Loss_for_Language_Image_Pre-Training_ICCV_2023_paper.pdf
- big_vision reference implementation / chunked issue (canonical JAX pseudocode + init): https://github.com/google-research/big_vision/issues/130
- SigLIP topic summary (loss eqn, batch-size decoupling): https://www.emergentmind.com/topics/sigmoid-loss-for-language-image-pre-training-siglip
- SigCLR, "Sigmoid Contrastive Learning of Visual Representations," small-batch evidence: https://arxiv.org/pdf/2410.17427

Code (file:line):
- `morpheus/v2/losses.py:13-18` — current `symmetric_infonce` (temp 0.07 hardcoded).
- `morpheus/v2/training.py:194-196` — the call site to replace (identity alignment on `z_identity`).
- `morpheus/v2/model.py:260,280` — `z_identity` is L2-normalized before export (SigLIP assumes normalized `x,y`).
- `morpheus/v2/model.py:202` — `self.identity = nn.Linear(hidden,256)`; head lives inside the `TumorStateV2` nn.Module.
- `morpheus/v2/model.py:57,67,198,214` — precedent for `nn.Parameter` scalars in the model (`log_temperature`, `anchor_residual_scale`).
- `morpheus/v2/training.py:318` — `clip_grad_norm_(self.model.parameters(), 1.0)`; `:319` optimizer.step over model params.
- `morpheus/v2/training.py:341` — checkpoint saves `self.model.state_dict()`.
- `morpheus/v2/tests/test_v21_model.py:19,62,98` — `_batch()` fixture, `V2Trainer` construction, profile test.

## Findings

### 1. Exact SigLIP loss and init (Zhai 2023; big_vision reference)
For a batch of B paired embeddings, with L2-normalized `x_i` (WSI) and `y_j` (RNA):

    logits[i,j] = t * (x_i . y_j) + b
    labels[i,j] = +1 if i==j else -1            # 2*I - 1
    L = -(1/B) * sum_i sum_j  log_sigmoid( labels[i,j] * logits[i,j] )

`log_sigmoid(z) = -softplus(-z)`; use `binary_cross_entropy_with_logits`
with target `(labels+1)/2` (i.e. 1 on diagonal, 0 off) which is numerically identical.

Init (canonical, big_vision):
- Temperature parameterized in log-space to keep it positive: store
  `t_prime = log_temperature = log(10) ≈ 2.3026`, use `t = exp(t_prime)` (so t starts at 10).
- Bias `b = -10.0`, learnable.
The bias init of -10 counteracts the heavy positive/negative imbalance at t=0
(B-1 negatives vs 1 positive per row): it starts every pair near "negative" so
the first optimization steps are not dominated by the flood of easy negatives.

### 2. Why it decouples from batch/candidate size
InfoNCE (`losses.py:16-18`) applies a softmax over each row, so every logit is
normalized against all other in-batch candidates — the loss for a pair depends on
the composition and *count* of the other candidates. The sigmoid loss applies an
independent Bernoulli term per pair; there is no cross-pair normalization. The
learnable `b` absorbs the global positive/negative ratio instead of the softmax
denominator doing it implicitly. Consequence for MORPHEUS: our identity candidate
sets are modest, and InfoNCE's effective difficulty and gradient scale drift with
batch size (`symmetric_infonce` even returns 0 for B<2, `losses.py:14`). Sigmoid
removes that coupling.

### 3. Small-batch evidence (directly relevant — our sets are modest)
- SigLIP paper: sigmoid "performs better at smaller batch sizes"; with <32k batch
  it outperforms the softmax WebLI baseline; performance saturates by 32k. The win
  is specifically in the small/medium regime, not just at extreme scale.
- SigCLR (2410.17427): sigmoid contrastive is competitive-to-better than
  softmax/InfoNCE at small batches (128-512 studied), where softmax needs many
  negatives to be stable; the learnable bias `b` is what lets it adapt the decision
  boundary without large negative pools.
This is the regime MORPHEUS is in, so the mechanism (not just scale) motivates the swap.

### 4. Where the learnable params MUST live (code constraint)
`symmetric_infonce` is a free function. Grad clipping (`training.py:318`),
optimizer (`:319`), and checkpointing (`:341`) all operate on
`self.model.parameters()` / `self.model.state_dict()`. Therefore `t` and `b` must
be `nn.Parameter`s on the model (`TumorStateV2`), not tensors created inside the
loss. There is exact precedent: `log_temperature` (`model.py:57`),
`anchor_residual_scale` (`model.py:214`). Putting them in the loss fn would give
them no gradient, no optimizer state, and no checkpoint entry — a silent no-op bug
of exactly the kind the test suite does not catch. `z_identity` is already
L2-normalized (`model.py:260,280`), so the loss needs no internal normalize.

## Recommended change (file:line, exact)

Smallest change that fixes the mechanism, in 3 edits:

**(a) `model.py`** — add two learnable scalars to `TumorStateV2.__init__`
(next to `self.identity`, `model.py:202`; mirror `log_temperature` at `:57`):

    # SigLIP identity-alignment params (Zhai 2023): t = exp(log_temp), b learnable.
    self.identity_logit_scale = nn.Parameter(torch.tensor(2.3026))   # log(10) -> t=10
    self.identity_logit_bias  = nn.Parameter(torch.tensor(-10.0))

These are picked up automatically by the optimizer, grad-clip, and state_dict.
(Optional: clamp `exp(scale)` to <=100 in the loss to match `model.py:67`.)

**(b) `losses.py`** — add alongside `symmetric_infonce` (after `losses.py:18`):

    def sigmoid_infonce(left, right, logit_scale, logit_bias):
        """SigLIP pairwise sigmoid loss (Zhai 2023). left/right already L2-normed
        upstream, but normalize defensively. logit_scale is log-temperature."""
        if len(left) < 2:
            return left.new_zeros(())
        t = logit_scale.exp().clamp(max=100.0)
        logits = (_norm(left) @ _norm(right).T) * t + logit_bias
        labels = 2.0 * torch.eye(len(left), device=left.device) - 1.0   # +1 diag, -1 off
        return -nn.functional.logsigmoid(labels * logits).mean()

(`nn.functional.binary_cross_entropy_with_logits(logits, (labels+1)/2)` is an
equivalent one-liner; the logsigmoid form matches the paper's algorithm exactly.)

**(c) `training.py:194`** — swap the call (add import `sigmoid_infonce` at `:13-15`):

    # before (training.py:194):
    identity = symmetric_infonce(out_wsi["z_identity"], out_rna["z_identity"])
    # after:
    identity = sigmoid_infonce(out_wsi["z_identity"], out_rna["z_identity"],
                               self.model.identity_logit_scale, self.model.identity_logit_bias)

Also log `metrics["identity_logit_scale"] = float(self.model.identity_logit_scale.exp().detach())`
and `metrics["identity_logit_bias"] = float(self.model.identity_logit_bias.detach())`
next to `training.py:196` so the learned t/b are observable (they diagnose regime health).

Note: the paper sums over j then means over i (`nll = -sum(loglik, axis=-1); l = mean(nll)`),
giving a per-row scale of ~B. The `.mean()` form above is the per-element mean
(scale-invariant to B) — preferable here because the identity weight (`identity_after_warmup=1.0`,
`training.py:26`) is tuned against InfoNCE's mean-scaled magnitude. Keep `.mean()`
to avoid silently multiplying the identity weight by batch size.

## A/B it against InfoNCE fairly

- Single toggle: add `identity_loss: str = "infonce"` to `V2LossSchedule`
  (`training.py:19-40`) and branch at `training.py:194`. Do NOT fork the trainer.
- Hold everything else fixed: same `identity` weight schedule (1.0/1.0), same
  optimizer, same batches/seed, same cancer-held-out folds (11 dev / 21 held-out),
  same warmup. The only diff is the loss at `:194` (+ the two params, which are
  inert under the InfoNCE arm).
- Fair magnitude: InfoNCE uses fixed temp 0.07 (`losses.py:13`); SigLIP *learns* t.
  Do not also re-tune InfoNCE's temp — that biases the comparison. Report the
  learned SigLIP t at convergence for interpretability.
- Metric: the existing within-cancer identity retrieval gate (SigLIP already beats
  MLP-CLIP +0.005 within-cancer per the lane context) plus identity effective rank
  (currently ~84/256) — confirm the swap does not degrade rank. Run each arm >=3
  seeds; the effect size is small so report mean +/- sd, not single runs.
- Add a regression test (new, in `tests/test_v21_model.py` using `_batch()`/`V2Trainer`):
  assert `identity_logit_scale.grad` and `identity_logit_bias.grad` are non-None and
  finite after one `step()` — this catches the "params created in loss fn get no
  gradient" no-op the current suite would miss.

## Risks & scaling

- **No-op risk (highest):** if t/b are instantiated inside the loss instead of on
  the model, they never train and never checkpoint, and no existing test catches it
  (lane context: suite misses no-op losses). Mitigated by placing them as
  `nn.Parameter` on `TumorStateV2` + the grad-finite test above.
- **Loss-magnitude mismatch:** paper's sum-over-j convention scales ~B; using it
  would inflate the identity term ~batch-fold vs InfoNCE and break the fixed 1.0
  weight balance against the programme losses (`training.py:186-225`). Use `.mean()`.
- **Bias saturation:** if positives are very rare relative to B, b can drift very
  negative and stall. Our candidate sets are modest (few negatives), so this is
  low-risk; the -10 init is already conservative. Monitor the logged `identity_logit_bias`.
- **Small-batch degenerate guard:** the `len<2 -> zeros` early return (kept from
  `symmetric_infonce`) preserves current ragged-batch behavior; SigLIP's per-pair
  independence otherwise handles small B gracefully (its main advantage here).
- **Scaling up:** if batch/candidate size later grows, SigLIP needs no chunking at
  our scale (chunked all-gather is a >many-GPU concern); on a single A10 the BxB
  logits are trivial. No change needed as batch grows within one device.
- **Out of scope (do not touch):** biology gets no paired contrastive by design;
  this swap is `z_identity`-only. Do not extend sigmoid to the programme/supcon
  paths. Dead residual/context wiring is a separate lane.
