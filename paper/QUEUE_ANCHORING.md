# A queued contrastive objective needs an independent reference frame, not fresher keys

> **SUPERSEDED AS A SUBMISSION UNIT, 2026-08-04.** This material is now **§5.2 of
> `paper/P2_RANK_DRAFT.md`**, by the author's decision. **Do not submit this file separately.**
>
> Three things the integration added and which this file should be read with. (i) **MoCo now carries an
> identifier** — He, Fan, Wu, Xie & Girshick, "Momentum Contrast for Unsupervised Visual Representation
> Learning", CVPR 2020, arXiv:1911.05722, verified at full text — where this repository previously had
> only a bare "He et al., 2020". (ii) **MoCo advances the staleness account twice as a hypothesis**
> (*"Our hypothesis is that…"*, *"We hypothesize that…"*), never as an established mechanism, and the
> falsification below is written to say so. (iii) MoCo ties the argument specifically to **queue** use
> (*"a slowly evolving key encoder is a core to making use of a queue"*), so the falsification is in
> scope only because ours is a queue setting — which is now stated rather than left to be inferred.
>
> **And the thing the integration exposed, updated 2026-08-04 now that the seed replication has
> reported.** The `m = 0` versus `m = 0.999` decision at step 600 rests on **2.81 against 7.42, a 2.64×
> ratio** — and this sweep is **one seed per momentum value**. The single-seed limitation was a
> **defect, not a design choice: the momentum harness had its seed hardcoded**, so the sweep could not
> have varied seeds had it been asked to; the "Limits" bullet below reading *"The sweep is one seed per
> momentum value"* is a property of the apparatus, not a judgement that one seed sufficed.
>
> **The seed replication has now run** (`~/e0_run/d1_diag/mseed_*`, three seeds per momentum, 500
> steps): canonical R1 **11.26 / 10.45 / 10.55** at m = 0.999 against **3.18 / 1.13 / 2.36** at m = 0;
> R3 **7.40 / 6.85 / 7.15** against **2.81 / 1.05 / 2.06**, all on the fixed held-out probe. **Every
> m = 0.999 seed exceeds every m = 0 seed on both statistics**, so the separation is not an artefact of
> the one seed the harness allowed.
>
> **But it does not clear P2 §4.1's bar either, and the file must not be read as though it did.** The
> worst-case separation is **10.45 / 3.18 = 3.29×** against §4.1's measured retraining floor of
> **3.295×** (canonical R1, residualised exported block; 3.111× raw). The earlier
> figure of 2.69× quoted in this header was the superseded n = 1 estimate and **must not be quoted as
> the floor**. Different arm, duration and block make the comparison indicative rather than a
> like-for-like disqualification, and all three mismatches point toward a *larger* floor in this regime
> rather than a smaller one. **The accurate verdict, and the word P2 §5.4 now uses throughout, is
> `unjudgeable` rather than "not resolvable"**: these runs are read on the fixed held-out probe, which
> has **no measured retraining floor** and cannot get one from the five exported repeats, so the
> criterion has not been applied to the comparison in either direction. §5.4's conclusion is unaffected
> — it rests on the binary outcome below.
>
> **And the mechanism is not momentum.** A predeclared learning-rate test (P2 §5.2a, `f68a7ac`) shows
> this collapse is primarily a **learning-rate** phenomenon: at `lr = 1e-3` the representation sits at
> 1.05–1.06 whatever the momentum, and at `lr = 4e-5` it reaches **12.30 with no momentum encoder at
> all**. Momentum is **neither necessary nor sufficient**; it does rescue the objective at the rate
> this project trains at, `2e-4`, and **lowering the learning rate would have solved the original
> problem more simply, and we did not try it.** Every threshold in `m` in this file is scoped to
> `lr = 2e-4`.
>
> **What the fix actually rests on is a binary outcome, not a ratio** — see **P2 §5.4**, which is now a
> subsection of its own. `programme_free` completed 40 epochs uncollapsed on **0 of 3** seeds before the
> fix (the one arm that reached epoch 40 sat at R3 **1.71** with RNA-view mutual cosine **0.986**, and
> no exports, CALIBRA or bootstrap were produced at all) and on **3 of 3** after it, with canonical R1
> 13.418 / 7.600 / 6.394, a held-out channel of 0.5412 / 0.5336 / 0.5126 above its own `random_control`,
> and paired bootstrap intervals per seed. Rank is how the collapse was **found**; it is not what
> establishes the repair, and no sentence in this file may be quoted as though it were. The choice of
> `m = 0.999` over `m = 0.99` (1.26×) is supported by nothing that clears the floor; only **momentum
> against none** is.

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
the key-writing encoder differs. All four arms start from the same initialisation (verified: centred
effective rank 67.55 at step 0). The full sweep, not just the chosen value — 40 epochs of the training
this objective is used for amounts to **583 steps**, so the table spans the real duration:

| step | m = 0 | m = 0.9 | m = 0.99 | m = 0.999 |
|---:|---:|---:|---:|---:|
| 0 | 67.55 | 67.55 | 67.55 | 67.55 |
| 50 | 4.10 | 3.88 | 8.65 | 9.35 |
| 100 | 1.62 | 3.51 | 6.49 | 7.03 |
| 150 | 1.62 | 2.15 | 4.56 | 6.99 |
| 200 | 2.26 | 1.65 | 5.70 | 7.60 |
| 300 | 3.32 | 2.70 | 6.01 | 7.33 |
| 400 | 2.18 | 2.31 | 5.50 | 7.84 |
| 500 | 2.43 | 2.34 | 5.50 | 7.61 |
| **600** | **2.81** | **2.23** | **5.88** | **7.42** |

Three things this table shows that a single number would not. The effect is **monotone in m** and
large — the per-step fold `m = 0.999 / m = 0` is **3.363× (200), 2.208× (300), 3.596× (400), 3.132×
(500), 2.641× (600)**, a range of **2.208×–3.596×** past step 150, and **4.340×** at step 100.
*An earlier version of this sentence said "2.6–3.3× at every step past 150"; **both ends of that range
are wrong against the table above** and it is corrected here (draft §4.1a, row 47).* It is **durable**:
both working arms are flat from step 200 to
600, spanning the full 583-step training duration, which matters because two earlier "fixes" on this
objective looked correct inside a short window and failed outside it. And `m = 0.9` **fails**, tracking
the no-momentum arm rather than the working ones, so this is not a matter of perturbing the key
encoder slightly — there is a threshold in m **at this learning rate**, and it lies between m = 0.9 and
m = 0.99. A predeclared learning-rate test has since shown that the threshold moves with the learning
rate and is not a property of m at all: see draft §5.2a.

**`m = 0.999` is used because it measured best in this sweep. No mechanism is claimed here.** That is a
weaker justification than a hyperparameter usually receives and it is stated deliberately: **three**
mechanisms were proposed for this effect and all three were falsified by measurement (below), and the
fourth — the momentum threshold this table shows — was falsified in turn by the learning-rate test in
P2 §5.2a. A reader can
see exactly how much to trust the choice, which is more than a value defended by an unverified story
would offer.

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

**It is not confirmed, and one sharpening of it has been tested and refuted.** We proposed, and
predeclared, that the required momentum should track queue turnover — that the EMA time constant
`τ = 1/(1−m)` must exceed `T = capacity/batch`, making the criterion dimensionless in `τ/T`. Five
predictions were committed before the runs. **Four were wrong, and the discriminating one produced no
effect at all:** the criterion required the arm with the *lowest* `τ/T` (0.52, capacity 8192, `m =
0.95`) to be the worst of its group and fixed *"fails ≤ 3.5"* in advance; it read **3.67** and so did
not fail. *An earlier version of this document called that an inversion, reading an ordering off 3.67
against 3.53 — a 1.04× difference, smaller than any floor this project has measured on any statistic
or any view. The accurate statement is that a predicted effect was **absent**; the same two numbers do
carry the equality claim "nearly flat in capacity" below, which an ordering claim cannot (draft §4.1a,
rows 51–52).*

What the data supports is narrower: rank is **monotone in `m` and nearly flat in capacity**. At fixed
`m = 0.95`, tripling the queue from 2,048 to 8,192 moves effective rank 3.53 → 3.67. At fixed capacity
4,096, raising `m` from 0.95 to 0.999 moves it 2.91 → 7.82. A threshold sits between `τ = 20` (fails)
and `τ = 100` (works) and appears at the same place for every capacity tested — an *absolute* time
constant, not a ratio **at the one learning rate this sweep was run at**. We recorded that as an
observation and explicitly did **not** advance it as a mechanism.

**And it has since been falsified, by the experiment the Limits section below named as the next
testable shape.** A predeclared learning-rate test (draft §5.2a, `f68a7ac`) ran six arms at 400 steps
from one initialisation. At `lr = 1e-3` the representation sits at **1.06 / 1.05 / 1.05** across
`m = 0 / 0.9 / 0.999` — the momentum threshold buys nothing, a spread of **1.01×**; at `lr = 4e-5` it
reaches **12.30** with `m = 0`, **27.88** at `m = 0.99` and **35.24** at `m = 0.999`. **Learning rate
predicts the outcome perfectly and momentum predicts none of it at the high rate.** So the `τ`
threshold above is a property of the learning rate it was measured at, not of `m`: **momentum is
neither necessary nor sufficient** for this rescue. At the rate this project actually trains at,
`2e-4`, momentum does rescue the objective and the seed-varied replication is unambiguous — but
**lowering the learning rate would have solved the original problem more simply, and we did not try
it.** Everything in this document about a threshold in `m` must be read with "at `lr = 2e-4`"
attached.

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
- Durability is established only to **step 600**, against a training duration of 583 steps. Nothing is
  known about behaviour at 10× that horizon.
- The sweep is **one seed per momentum value**. The arms are separated by **2.208×–3.596×** past step
  150 and are flat over 400
  steps, so the ordering is not in doubt, but no variability estimate is offered and none should be
  read in. *A seed replication has since run — three seeds per momentum at 500 steps, every m = 0.999
  seed above every m = 0 seed — and its worst-case separation is 3.29× against the nearest measured
  floor of 3.295×, which is on a different block. The block these runs are read on has no measured
  floor, so that comparison is **unjudgeable** by the paper's own criterion, not failing. Draft §5.4.*
- The anchoring account is **not** confirmed against a control, and as stated it makes no prediction
  about how much decoupling is required — which is simultaneously why it survived the turnover
  experiment and why it is currently unfalsifiable. It should not be cited as established.
- The turnover/`τ-over-T` sharpening was predeclared, tested and **refuted**; it is withdrawn.
- ~~The next testable shape, not yet run: vary the learning rate.~~ **RUN, predeclared, and it
  falsified the `τ` threshold** (draft §5.2a). The design was right — the learning rate moves per-step
  drift without touching `τ` or `T` — and neither of the two accounts it was built to separate
  survived: the outcome is a **pure learning-rate effect**, which the original design could not see
  because every high-rate arm it ran had sub-threshold momentum and every low-rate arm had
  supra-threshold momentum. The two missing cells were predeclared and run. This is the **fourth**
  account proposed for this collapse and the **first to survive a predeclared test**.
- The capacity effect (§3) remains confounded: capacity changes both anchoring quality and negative
  count, and we have not separated them.
- The staleness metric reported measures key-to-*query* lag, not key-to-*key* spread. With a 19-step
  turnover the latter is small in every arm, which is itself part of the argument, but it is not the
  same quantity MoCo discusses.
