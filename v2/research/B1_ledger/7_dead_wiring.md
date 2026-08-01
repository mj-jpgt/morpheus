# Lane 7: Dead-wiring audit and remediation priority

## Queries/searches run
- WebSearch: "effective rank collapse contrastive representation normalized vs unnormalized embedding regression head dimensional collapse"
- WebSearch: "attention pooling vs mean pooling query slots set transformer effective rank"
- Read: morpheus/v2/model.py, losses.py, training.py, export.py, contracts.py (full)
- Read: morpheus/v2/runner.py (coordinate loading 43-149; profiles/targets 168-318)
- Grep: `coordinate_present|coordinates`, `z_context|z_uncertainty|residual|z_semantic|semantic` across morpheus/v2

## Sources
Web:
- Set Transformer (PMA) https://arxiv.org/pdf/1810.00825 — learned attention pooling (PMA) beats mean/max for permutation-invariant set aggregation.
- RankMe https://arxiv.org/pdf/2210.02885 — effective rank predicts downstream quality; low rank == under-used representation.
- Understanding Dimensional Collapse in Contrastive Learning https://openreview.net/pdf?id=YevsQ05DEN7 — collapse arises when the supervised subspace is low-dimensional.
- ContraNorm https://arxiv.org/pdf/2303.06562 — normalization/whitening choices directly change effective rank.

Code (file:line):
- model.py:254 `take()` = `query[:, cursor:cursor+count].mean(1)` — MEAN pool over slots.
- model.py:257-259 identity(4).mean, biology(4).mean, context(2), 5x residual(4)=20, uncertainty(2).
- model.py:260,286 z_identity, z_biology L2-normalized; model.py:287 z_context NOT normalized.
- model.py:283 `biology_state = self.biology(biology)` (UNNORMALIZED).
- model.py:290 `programme_mean/log_variance = Linear(biology_state)` — heads read the UNNORMALIZED state.
- model.py:291 `z_uncertainty`, `z_{name}_residual` emitted but never consumed.
- model.py:293-295 semantic head only if `semantic_dim>0` (default 0 -> off).
- model.py:117,132-139 coordinate path; runner.py:43-50 `_metadata_coordinates` returns None unless real varying x/y exist -> coordinate_present_fraction=0 in practice.
- losses.py:46-53 KL neighbourhood on `_norm(z_biology)`; losses.py:56-66 supcon on `_norm(z_biology)`; losses.py:37-43 separation whitened; losses.py:29-34 variance floor per-dim.
- training.py:104 memory bank stores `states.detach()` (no backprop); training.py:95 read side normalizes bank.
- training.py:148-151 Gaussian-NLL on `programme_mean/log_variance` (from unnormalized state), weight `programme`=1.0.
- export.py:54-55 widths dict has only identity/biology/patient/semantic — no context/uncertainty/residual.

## Findings

Ranked by plausible contribution to the observed biology collapse (effective rank ~5-6 of 256).

### FINDING 1 (PRIMARY) — Normalize inconsistency drives the rank collapse. BUG-ish, high impact.
The programme Gaussian-NLL (weight 1.0, the dominant biology term) is computed from the
**unnormalized** `biology_state` (model.py:283 -> 290 -> training.py:148). The target is a
50-D cancer-residualized Hallmark vector. NLL at weight 1.0 with a free `log_variance` head
directly regresses a 256-D linear projection onto a <=50-effective-dim target: the optimal
solution puts signal in ~<=50 directions and lets `log_variance` absorb the rest. This is a
textbook low-dimensional-supervision collapse (openreview YevsQ05DEN7; RankMe). Meanwhile
the rank-spreading structural terms (KL, supcon, separation, variance) all run on the
*post-normalization* `z_biology` (losses.py:46-66). Normalization is a nonlinear rescale;
gradients on the normalized state cannot restore rank that the strong unnormalized-NLL term
removed on the raw axis, and their weights are tiny (neighbourhood 0.20, supcon 0.20,
separation 0.01, variance 0.01) vs programme 1.0. Net: one strong low-rank objective on the
raw state, four weak spreaders on a different (normalized) view of it. The variance floor
(losses.py:29-34) is explicitly per-dimension and, as the codebase facts note, does NOT
prevent rank collapse — it only stops any single dim from having zero variance, which a
rank-6 solution easily satisfies. Identity avoids this because its only strong loss (InfoNCE)
is on the normalized state itself, matching where separation/variance act — consistent with
its healthy rank ~84.

### FINDING 2 — 20 residual slots + z_context + z_uncertainty are dead. INTENTIONAL but wasteful.
32 query slots: 24 (75%) feed dead heads. residual/context/uncertainty are computed
(model.py:258-259, 287, 291) but never in any loss (training.py has no reference) and never
exported (export.py:54-55 widths). They still consume cross-attention capacity in every
QueryBlock and dilute the shared query bank the gradient-conflict monitor tracks. Does NOT
directly cause collapse, but reclaiming slots for identity/biology is the cheapest capacity
win. Contracts.py deliberately refuses to export untrained states, so this is intentional
hygiene, not an accidental leak — but the slots are pure waste as wired.

### FINDING 3 — mean-pool over 4 slots. INTENTIONAL, minor. Not the collapse cause.
`take()` mean-pools 4 identity / 4 biology / 2 context slots (model.py:254). Literature
(Set Transformer PMA, arXiv:1810.00825) shows learned attention pooling > mean for set
aggregation, but mean over 4 already-attended slots is not the bottleneck: the slots pass
through 4 QueryBlocks of cross+self attention first, and identity (same pooling) has rank 84.
Mean-pooling does not cause the biology collapse; Finding 1 does. Low priority.

### FINDING 4 — detached memory bank. CORRECT as designed. Not a bug.
training.py:104 stores `states.detach()`. This is standard MoCo-style queue behaviour: the
bank is a set of negatives/positives from *past* steps whose graph is gone; back-propping into
stale activations is impossible and undesirable. The current-batch anchors (states[active],
training.py:95) DO carry gradient. No change needed. Note it also does not spread rank because
supcon is a weak 0.20 term on the normalized view (subsumed by Finding 1).

### FINDING 5 — coordinate / semantic paths gated off. INTENTIONAL, out of lane scope to force on.
coordinate_present_fraction=0 (runner.py:43-50 guards against manufactured zeros — correct
leakage hygiene); semantic off unless PLIP cache supplied. Neither is wired-wrong; both are
capability-gated. Not collapse-relevant. Leave as-is.

## Recommended change (file:line, exact)

Minimal, single-mechanism fix targeting Finding 1 — make the biology regression supervise the
**same normalized state** the structural/collapse-preventing losses act on, so all biology
gradients share one geometry and the rank-spreaders are not fighting a different view.

model.py:283 — change the programme heads to read the normalized biology state.
  Current:
    biology_state = self.biology(biology)
    ... "z_biology": nn.functional.normalize(biology_state, dim=-1),
    ... "programme_mean": self.programme_mean(biology_state),
        "programme_log_variance": self.programme_log_variance(biology_state).clamp(-8, 8),
  Change to (compute the normalized state once, feed it everywhere):
    biology_state = nn.functional.normalize(self.biology(biology), dim=-1)
    ... "z_biology": biology_state,
    ... "programme_mean": self.programme_mean(biology_state),
        "programme_log_variance": self.programme_log_variance(biology_state).clamp(-8, 8),

This is a ~2-line diff. It removes the raw-axis low-rank pull: the NLL head now regresses from
the unit sphere, and separation/variance/KL/supcon all operate on that identical vector, so the
weak rank-spreaders finally act on the same geometry the strong NLL shapes. Because the
programme target is residualized/scaled on the train fold (runner.py:168-172), the head can
still fit it from a normalized input (the two Linear heads absorb the scale).

Optional, only if Finding-1 fix alone under-delivers rank (do NOT bundle): raise the variance
floor from per-dim to a rank surrogate by adding an off-diagonal covariance penalty on
z_biology (a VICReg-style covariance term) at losses.py, weight ~0.04. Keep it a separate
follow-up so the primary fix is attributable.

Explicitly NOT recommended now: reclaiming the 20 residual slots (Finding 2), switching mean->
attention pool (Finding 3), un-detaching the bank (Finding 4) — none address the collapse
mechanism and each adds parameters/risk. Slot reclamation is a good but separate capacity PR.

## Risks & scaling
- Behaviour change on resume: the biology head now sees unit-norm input, so a checkpoint trained
  pre-fix is not weight-compatible in effect (numerically loads, but the head must re-fit the
  scale). Retrain from scratch or from the warmup; do not hot-swap into a late checkpoint.
- Gaussian-NLL on a normalized input: because ||biology_state||=1, per-programme target magnitude
  is now bounded by the head weights, not the state. If targets have large residualized scale the
  `programme_mean` weights grow — acceptable, but watch `programme_log_variance` doesn't peg at
  the +8 clamp (a mode-collapse signature). Add a test asserting mean log_variance stays inside
  (-8, 8) and that biology effective rank rises above ~20 on a small fixture.
- No test currently catches rank collapse or log_variance->clamp (contracts.py only checks
  isfinite; assert_module_used only checks a head got *any* gradient). Add a rank/log-var guard
  test using the existing `_batch()` fixture (tests/test_v21_model.py) before landing.
- Scaling: change is O(1), no memory/throughput impact on A10+; normalize is already computed for
  z_biology so there is zero added cost — it is literally reusing one existing tensor.
- Leakage protocol unaffected: no new data path, transforms still train-fold-fit (runner.py:168).
