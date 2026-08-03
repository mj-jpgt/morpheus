# A queued contrastive objective needs an independent reference frame, not fresher keys

*Draft section, 2026-08-03. Companion to `LIVENESS_GATE_DESIGN.md`. Every number traces to a named
entry in `NOTEBOOK_ENTRIES/` or a log under `~/e0_run/d1_diag/`. Claims not measured are marked as
such. This section reports a fix that works and an explanation that does **not**, and is written that
way deliberately.*

---

## The failing configuration

A patient-paired cross-modal InfoNCE with a queue of detached negative keys. The queue is refreshed
every step from the **query encoder itself** — the standard "end-to-end with a memory bank"
arrangement:

```python
self.biology_memory.update(out_wsi["z_biology"], out_rna["z_biology"], batch["indices"])
```

At cohort scale (3,118 training patients, 4,096-key queue, 214-patient batches) the representation
collapses from effective rank 67.55 at initialisation to **~2 within 150 steps** and does not recover.
The collapse is invisible to a memorisation gate: the same objective, checked on 16 fixed patients
against a *frozen* 64-key queue, reaches effective rank 5.81 with 16/16 retrieval.

The collapse is **not** a regularisation failure. Holding the queue fixed and sweeping the two
candidate terms — a covariance penalty and a fused-state consistency term — over
`decorrelation ∈ {0, 0.04}` × `consistency ∈ {0, 0.1, 1.0}`, every one of five configurations collapses
to effective rank 1.6–2.4 by step 150, including both terms at zero simultaneously.

## The fix

Give the queue its own **momentum key encoder**: an EMA copy of the query encoder,
`θ_k ← m·θ_k + (1−m)·θ_q`, kept out of the optimiser, used only to encode keys before enqueueing.

Measured with **capacity held at 4,096 in every arm**, so the number of negatives is identical and only
the key-writing encoder differs. All arms from the same initialisation (verified: 67.55 at step 0):

| key encoder | eff-rank step 50 | step 100 |
|---|---:|---:|
| none — query encoder writes keys | 4.07 | **2.58** |
| EMA, m = 0.99 | 8.63 | **6.65** |
| EMA, m = 0.999 | 9.33 | **6.89** |

**2.6× at fixed negative count.** *Longer runs to 1,500 steps — 2.5× the duration of the 40-epoch
training this objective is used for — are in progress and are a precondition for any claim that the
recovery is durable; see Limits.*

## The explanation we expected, and why it is wrong

The natural account is MoCo's (He et al., 2020): a long queue holds keys written by an encoder that
has since moved, the keys are therefore mutually inconsistent, and the cheapest escape from an
inconsistent key set is to stop distinguishing anything. We predicted, in advance, that a momentum
encoder would flatten a directly-measured key-staleness curve.

**It did not, and three independent measurements rule the account out.**

**1. The queue is never stale.** With 214-patient batches into 4,096 slots, the queue **turns over
completely every 19 steps**. No key is ever more than nineteen steps old. This is the argument that
does the most work and it is available to anyone by arithmetic, before any experiment: *a memory bank
is only "stale" relative to its turnover rate, and turnover is batch size over capacity.*

**2. Key-to-encoder agreement does not predict rank.** Measuring cosine(stored key, fresh re-encoding
of the same patient by the current query encoder):

| key encoder | agreement, step 100 | eff-rank, step 100 |
|---|---:|---:|
| none | 0.427 | 2.58 |
| m = 0.99 | **0.908** | 6.65 |
| m = 0.999 | 0.441 | **6.89** |

The best-agreeing arm does not have the best rank; the best-ranked arm agrees no better than the
failing one.

**3. The healthiest queue holds the freshest keys.** Sweeping capacity at fixed key encoder, the
strongest arm is capacity 64 — which, at 214 patients per batch, is *entirely overwritten every single
step*. If staleness were the mechanism this arm should be the worst; it is the best of the sweep
(effective rank 6.17 versus 2.16 at capacity 4,096).

## What we think is happening instead

When the query encoder also writes the keys, a transformation applied to **all patients at once** —
contracting toward a shared direction, say — moves the queries and the keys *identically*. The
similarity structure between them is preserved, so the loss does not penalise it. **Collapse is free.**

A decoupled key encoder holds a slowly-moving reference frame. A global transformation of the queries
now changes their similarity to that frame, and therefore costs loss. On this account the queue's
defect is that it supplies **no independent frame**, not that its contents are old — an *anchoring*
story rather than a *staleness* story.

This account is consistent with all three measurements above, and with the memorisation gate passing
on a frozen queue (a frozen queue is a perfect anchor). It also explains why no loss-weight setting
helped: the degenerate direction is invisible to the loss by construction, so no reweighting of terms
that are computed *within* that loss can see it either.

**It is not yet confirmed.** Its distinguishing prediction is about *how much* decoupling is needed. If
anchoring strength drives the effect, larger `m` should help monotonically. Measured, m = 0.999 gives
6.89 against m = 0.99's 6.65 — barely separated, suggesting **decoupling at all matters far more than
how much**. That further distinguishes the account from staleness, where the amount of lag would be
the entire story. An `m = 0.9` arm is running as the sharper test: if it also works, decoupling-at-all
is confirmed.

## Why this may matter beyond one codebase

The failing arrangement — a queue refreshed from the query encoder each step — is a common
simplification of MoCo, adopted precisely because it avoids maintaining a second encoder. If the
reason MoCo's second encoder is necessary is *anchoring* rather than *staleness*, then the
simplification is unsafe **even when the queue is short enough that staleness is impossible**, which is
exactly the regime where practitioners reason it should be harmless. The turnover-rate calculation in
§1 is the cheap diagnostic: it tells you staleness cannot be your problem, and on the anchoring account
that is no reassurance at all.

## Limits

- One objective, one architecture, one cohort. No claim about generality is made.
- **The durability of the fix is not yet established.** Recovery is measured to step 100; runs to 1,500
  steps are in progress. This project has twice had a fix look correct inside a short window and fail
  outside it, so this section may not be cited for the fix until those return.
- The anchoring account is **not** confirmed against a control. The decisive experiment — decoupling
  the key encoder while holding it at a *fixed* lag, versus a genuinely stale queue at matched lag —
  has not been run.
- The capacity effect (§3) remains confounded: capacity changes both anchoring quality and negative
  count, and we have not separated them.
- The staleness metric reported measures key-to-*query* lag, not key-to-*key* spread. With a 19-step
  turnover the latter is small in every arm, which is itself part of the argument, but it is not the
  same quantity MoCo discusses.
