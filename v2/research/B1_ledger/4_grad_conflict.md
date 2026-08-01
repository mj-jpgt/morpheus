# Lane 4 — Gradient-conflict instrumentation + mitigation

Remit: cheapest correct way to measure per-loss gradient conflict (pairwise cosine of
per-task gradients at the shared trunk) for one epoch, then the minimal mitigation.
Ground truth: partial impl already exists at `training.py:264-305`.

## Queries/searches run

- WebSearch: "PCGrad CAGrad GradNorm GradVac multi-task gradient conflict cosine similarity comparison cost"
- WebSearch: "PCGrad number of backward passes per step cost O(T) tasks GradVac EMA cosine target implementation"
- WebFetch: PCGrad NeurIPS 2020 PDF (Yu et al.) — binary/unreadable, fell back to known algorithm + secondary sources
- Read: `morpheus/v2/training.py:180-262` (step / loss assembly), `:264-305` (existing conflict metric), `:307-323` (train_epoch call site)
- Read: `morpheus/v2/losses.py:1-86` (all objective terms)
- Read: `morpheus/v2/model.py:52-72, 193-289` (shared trunk, slot split, heads)
- Grep: `gradient_diagnostics_every`, `_gradient_conflict_metrics`, trunk params — call sites and defaults

## Sources

Web:
- Yu et al., "Gradient Surgery for Multi-Task Learning" (PCGrad), NeurIPS 2020 — https://proceedings.neurips.cc/paper/2020/file/3fe78a8acf5fda99de95303940a2420c-Paper.pdf . PCGrad: for each ordered task pair with `cos(g_i,g_j)<0`, project `g_i -= (g_i·g_j / ||g_j||^2) g_j`. Cost = one backward per task (O(T) backward passes/step).
- GradVac (Wang et al., ICLR 2021) via LibMTL docs — https://www.aidoczh.com/libmtl/_modules/LibMTL/weighting/GradVac.html . Generalizes PCGrad: aligns whenever `cos < phi_ij` where `phi_ij` is an EMA target per pair; PCGrad is the special case `phi=0`. Same O(T)-backward cost + per-pair EMA state.
- Reference PyTorch PCGrad/GradVac impl — https://github.com/anzeyimana/Pytorch-PCGrad-GradVac-AMP-GradAccum (confirms per-task backward + AMP handling).
- Cost/comparison summary — https://arxiv.org/html/2509.16959v1 , https://web.stanford.edu/class/archive/cs/cs224n/cs224n.1244/final-projects/DavidSaykinKfirShmuelDolev.pdf : gradient-surgery (PCGrad/CAGrad) is O(T) backward; GradNorm is O(1) backward but rebalances *magnitudes only* (does not fix direction conflict).

Code (file:line):
- `training.py:255-261` — `_last_loss_components`: only 5 buckets: `identity`, `programme`, `reconstruction`, `separation`, `variance`. All three biology views + all three programme sub-losses are collapsed into one `programme` scalar (`programme_total`, `:257`).
- `training.py:212-221` — `programme_total` sums wsi/rna/full programme losses; each `_programme_loss` internally sums Gaussian-NLL + neighbour-KL + supcon.
- `training.py:264-305` — existing conflict metric. `:273-277` shared-param set = `[model.queries]` + `wsi.patch[-1].weight` + `rna.projection[-1].weight`.
- `training.py:282` — `torch.autograd.grad(component, shared, retain_graph=True)` per component: correct, cheap (see cost below).
- `training.py:315-316`, `:127` — sampled every `gradient_diagnostics_every` steps (default 25); `runner.py:322-323` guards `>=0`.
- `model.py:52,66, 198-199, 247-259, 250` — shared trunk = `self.queries` (orthogonal-init slots) → `self.blocks` cross-attn → `self.norm`, then `take()` splits the 32 slots into identity/biology/context/residual/uncertainty at `:257-259`. `self.identity`/`self.biology` heads at `:260,283` are the split point.
- `losses.py:46-53` (neighbour-KL, temp 0.20), `:56-66` (supcon), `:69-85` (Gaussian-NLL), `:13-18` (identity InfoNCE), `:37-43` (separation), `:29-34` (variance floor).

## Findings

1. **Measurement is 80% built but mis-scoped in two ways.** The mechanism (per-component
   `autograd.grad` against shared params, then pairwise cosine, `training.py:282-301`) is
   correct and cheap. Two defects make it unable to see the failure this lane targets:

   (a) **Wrong trunk params.** The shared trunk where the 32 slots are jointly produced and
   then split is `self.queries` + `self.blocks` + `self.norm` (`model.py:247-250`). The
   existing set (`training.py:273-277`) uses `model.queries` (correct) but adds
   `wsi.patch[-1].weight` and `rna.projection[-1].weight` — those are *pre-trunk, per-view
   input encoders*, not the shared trunk, and they dilute/confound the cosine. The
   biology-vs-identity split happens **after** the trunk, so conflict must be read at the
   trunk output, not at the view encoders. Fix: shared set = `self.model.queries` plus the
   parameters of `self.model.blocks` and `self.model.norm`.

   (b) **Granularity too coarse to see the real conflict.** The observed failure is biology
   effective rank ~5-6 (rank collapse) while identity is ~84. The plausible cause is
   conflict *among the three biology programme losses* (Gaussian-NLL vs neighbour-KL vs
   supcon) and between biology and identity. But `_last_loss_components` (`training.py:255-261`)
   collapses all biology programme terms into a single `programme` scalar (`programme_total`,
   `:257`). So the current metric can report identity↔programme cosine but is structurally
   blind to NLL↔KL, NLL↔supcon, KL↔supcon — exactly the pairs most likely to be
   collapsing the biology head. This is the load-bearing gap.

2. **Cost is negligible for measurement — do NOT gate it heavily.** Each extra component
   gradient is `autograd.grad` w.r.t. only the shared-trunk params (`self.queries` is
   32×512 ≈ 16K params; blocks/norm add a few M), with `retain_graph=True` reusing the
   single forward+loss graph already built by `step()`. This is **not** a second backward
   through the encoders — it is a partial backward to the trunk only, one per active
   component. With ~6 useful components that is ~6 cheap partial backwards vs the one full
   `loss.backward()`. On uncapped H-Optimus patches the forward (thousands of patch tokens
   through cross-attention) dominates; the trunk-only grads are a small fraction of one
   step. Conclusion: for a **one-epoch diagnostic** the honest thing is to run it **every
   step** (`gradient_diagnostics_every=1`) for that one epoch — the default 25 undersamples
   and median-over-sparse-samples hides variance. Sampling every step for one epoch is
   affordable; leave the default 25 for normal training.

3. **Mitigation ranking, fitted to this codebase (smallest fix first):**
   - **PCGrad** is the right first mitigation *if* measurement shows sustained negative
     cosine. It is drop-in at the existing loss-component boundary (`_last_loss_components`),
     needs no new hyperparameters, and is a strict superset of "do nothing" (it only edits
     gradients that actually conflict). Cost: O(T) partial backwards to the trunk per step —
     but full backward for the combined grad is still one pass, and PCGrad only needs
     per-component grads w.r.t. the shared trunk (same cheap grads already computed for the
     metric), then applies surgery before `optimizer.step()`. On uncapped patches the marginal
     cost is the same small fraction as the diagnostic.
   - **GradVac** = PCGrad + per-pair EMA target cosine. Strictly more machinery (EMA state
     per pair, a target φ) for a marginal gain; not justified until PCGrad is shown
     insufficient. Skip for the minimal change.
   - **CAGrad** solves a small QP per step for a worst-case-improvement direction. More
     compute + a `c` hyperparameter; over-engineered for a first cut.
   - **GradNorm** only rebalances loss *magnitudes* (O(1) backward) — it does **not** fix
     directional conflict, which is the measured mechanism. Wrong tool here.
   - **Plain loss reweighting** is the true floor: if the epoch of measurement shows the
     conflict is dominated by one over-weighted term (e.g. supcon or NLL swamping the trunk),
     a single weight change in `V2LossSchedule` beats any surgery. **Measure first; reweight
     if the conflict is magnitude-driven; PCGrad only if it is genuinely directional.**

4. **Recommendation: measurement first, mitigation second, in one small PR each.**
   The mechanism the mitigation must fix is not yet observed at the right granularity, so
   committing PCGrad now would be premature. The minimal correct move is to fix the metric
   (trunk params + biology-loss granularity), run one epoch at every-step sampling on
   uncapped patches, then decide reweight-vs-PCGrad from the actual cosines.

## Recommended change (file:line, exact)

**Change 1 (do now) — fix the measurement, `training.py:255-261` and `:273-277`.**

a) Expose biology programme sub-losses as separate components. In `_programme_loss`
   (return NLL, neighbour-KL, supcon separately) and in `step()` (`training.py:212-225`),
   populate `_last_loss_components` with `programme_nll`, `programme_neighbour`,
   `programme_supcon` (summed across the 3 views, keeping view-scale factors) instead of the
   single `programme` bucket at `:257`. Keep `identity`, `separation`, `variance`,
   `reconstruction` as-is. This turns the 5-bucket set into ~7 and makes the collapse-relevant
   pairs visible.

b) Fix the shared trunk set at `training.py:273-277` to the actual post-attention trunk:
```python
shared = [self.model.queries]
shared += [p for p in self.model.blocks.parameters() if p.requires_grad]
shared += [p for p in self.model.norm.parameters() if p.requires_grad]
```
   Remove `wsi.patch[-1].weight` and `rna.projection[-1].weight` (pre-trunk view encoders,
   `model.py:113,96` — they precede the split and confound the cosine).

c) For the one diagnostic epoch, run with `--gradient-diagnostics-every 1` (arg already
   exists, `runner.py:379`; default stays 25 for normal runs). No code change needed.

**Change 2 (only if Change-1 epoch shows sustained negative directional cosine, not
magnitude imbalance) — add PCGrad at the existing component boundary.** Reuse the exact
per-component trunk gradients already computed in `_gradient_conflict_metrics`; before
`optimizer.step()` (`training.py:319`), for each conflicting ordered pair
`cos(g_i,g_j)<0` apply `g_i -= (g_i·g_j)/(||g_j||^2 + eps) * g_j`, sum the surgically-
adjusted per-component grads, and write them into `.grad` on the shared-trunk params only
(leave view-encoder grads from the normal `loss.backward()` untouched). Gate behind a
`--pcgrad` flag defaulting off. If instead the epoch shows one term dominating by magnitude,
prefer a single weight edit in `V2LossSchedule` over any surgery.

## Risks & scaling

- **Autocast/AMP interaction:** the existing metric casts component grads to `float32`
  (`training.py:286`) — correct. PCGrad surgery must likewise run in fp32 to avoid the
  dot-product/`||g||^2` underflowing under bf16 autocast (the reference impl handles this).
- **`retain_graph=True` memory:** holding the graph for ~7 partial backwards raises peak
  activation memory for the sampled step. On uncapped H-Optimus patches (thousands of
  tokens) the forward activations already dominate; sampling every step for a single
  diagnostic epoch is fine, but leaving every-step *PCGrad* on for full training would keep
  the graph every step — validate peak memory on the A10 before enabling continuously.
  If tight, restrict PCGrad grads to `self.queries` + `self.norm` only (drop `blocks`) to
  shrink retained graph.
- **`allow_unused` / zero grads:** components that don't touch the trunk (e.g. reconstruction
  is off `z_identity` which flows through the trunk — it does) return `None`; existing code
  zero-fills (`training.py:286`). A zero-grad component contributes a degenerate cosine —
  the guard at `:280` (skip zero/non-finite components) already handles this; keep it.
- **Cosine-of-conflict is a signal, not a verdict:** rank collapse (biology rank ~5-6) can
  arise from the *supcon+KL objectives themselves* being low-rank targets (~50-D Hallmark,
  losses.py:46-53) rather than from gradient conflict. Measurement may show *little*
  trunk-level conflict — in which case the correct conclusion is "not a gradient-conflict
  problem, hand off to the rank-collapse / objective-design lane," and no mitigation should
  be added here. State this explicitly so a null result is not overridden by shipping PCGrad
  anyway.
- **Scaling to more losses:** cost is O(active components) partial backwards; adding the
  currently-dead objectives (context, residuals) later would grow this linearly. Fine at
  ~7 components; revisit if it exceeds ~15.
