# What a liveness gate certifies — and four ways we learned it certifies less than we assumed

*Draft section, 2026-08-03. Companion to `P1_CALIBRA_DRAFT.md`, which asks what an **analysis** would
have missed; this asks what a **training gate** would have missed. Every number traces to a named
entry in `NOTEBOOK_ENTRIES/` or a log under `~/e0_run/`. Numbers not measured are marked as such.*

---

## Claim

A liveness gate that certifies a model **can fit** a small, fixed, favourable problem does **not**
certify that the objective **will learn** at the scale and duration it is about to run at. The two
differ whenever the training regime contains a dynamic the gate's regime removes — and gates are
usually designed to remove exactly such dynamics, because that is what makes them fast and
deterministic.

We report four independent instances from a single objective over two days. Each was believed to be
a completed fix at the time; each was falsified by extending the gate's regime toward the run's.

## Instance 1 — the window was shorter than the failure

`programme_free`'s biology head was observed collapsing under a covariance penalty. A per-dimension
variance floor was added and the collapse disappeared: at 800 steps the contrastive term descended to
2.0875 with patient-to-patient cosine 0.7261. Extending the same run to 5,000 steps reversed the
verdict — 2.4626 at step 1,000, then pinned at ln(16) = 2.7726 from step 2,500 onward with
patient cosine **1.0000**.

The fix delayed the failure past the observation window. Nothing about the 800-step measurement was
wrong; it was answering "has it collapsed *yet*".

*Provenance: `g26_variance_floor_fix`, `g26_stepbudget_sweep`.*

## Instance 2 — the batch was smaller than the problem

The repaired gate memorises one fixed 16-patient batch. With the objective's regularisers excluded it
reaches contrastive 0.012–0.057, retrieval 16/16, patient cosine 0.0597 and effective rank 5.81 — a
clean pass on a criterion of ≤ 0.10, on three seeds.

The same objective at cohort scale, 3,118 streaming patients, collapses to effective rank ~1.8 by
epoch 21 (against the supervised arm's 7.38 / 7.35 on two seeds, measured on 282 held-out patients).

Memorising sixteen items is a capacity question. Representing three thousand is a learning question.
The gate answered the first and was read as answering the second.

*Provenance: `g26_passes`, `d1_programme_free_collapsing_in_training`.*

## Instance 3 — the gate removed the dynamic that causes the failure

The subtlest of the four, because the removal was deliberate, documented, and correct on its own terms.

The objective uses a queue of detached negative keys. Replaying one batch against a *live* queue makes
the queue fill with re-encodings of that same batch, so the negatives become the queries; the gate
therefore **freezes** the queue at 64 keys. That fix was correct and necessary — measured reduction
0.054 → 0.394, unique keys 16 → 64.

But a frozen queue is perfectly self-consistent, and training's queue is not: it holds 4,096 keys
written by an encoder that has moved substantially since the oldest were enqueued. Under a live queue
the same objective collapses to effective rank ~2 within 150 steps, from an initialisation of 67.55,
under **every** setting of both regularisers suspected of causing it — including both at zero:

| `decorrelation` | `biology_full_consistency` | step 50 | 100 | 150 | 200 | 250 |
|---|---|---:|---:|---:|---:|---:|
| 0.04 | 1.0 | 4.08 | 1.95 | 2.16 | 1.68 | 1.59 |
| 0.0 | 1.0 | 2.62 | 2.16 | 2.47 | 1.94 | 2.17 |
| 0.04 | 0.1 | 2.99 | 3.43 | — | — | — |
| 0.0 | 0.1 | 2.97 | 2.00 | 2.50 | — | — |
| 0.0 | 0.0 | 2.98 | 1.98 | 1.86 | — | — |

The pathology is therefore not in the objective's weighting. It is in the key set: giving the queue its
own slowly-moving EMA encoder, **at fixed capacity 4096 so the negative count is unchanged**, raises
effective rank from 2.58 to 6.89 by step 100.

The mechanism is *not* the one we predicted, and the check that established this is the part worth
copying. We expected MoCo's account (He et al., 2020) — a long queue holding keys written by a
substantially different encoder — and predicted in advance that a momentum encoder would flatten the
measured key-staleness curve. It did not. Measured directly, the queue turns over completely every
**19 steps** (214-patient batches into 4,096 slots), so no key is ever meaningfully old; and
key-to-encoder agreement fails to predict rank, the best-agreeing arm (0.908) having *lower* rank than
the worst-agreeing one (0.441). The capacity sweep's healthiest arm holds the *freshest* keys of all,
being fully overwritten every step.

What the evidence supports instead is that a queue written by the query encoder lets the loss be
reduced by moving queries and keys *together*, and that a decoupled key encoder removes that escape.
We report the fix as effective and the mechanism as open, rather than attaching a familiar explanation
that our own measurement contradicts.

The gate froze the queue to remove a *known* pathology, and in doing so removed the *unknown* one.

*Provenance: `d1b_premise_fails_all_five_arms_collapse`, `queue_size_implicates_the_key_set`,
`momentum_rescues_rank_but_staleness_is_not_the_mechanism`.*

## Instance 4 — the gate was read from a re-implementation of itself

The sharpest of the four, and the cheapest to have avoided.

To decide whether the objective was ready to launch, we ran the gate function
(`_overfit_programme_free_contrastive`) from a standalone harness that reconstructed its inputs:
cohort, split, schedule, model. Three seeds passed with margin, so the run was launched. The gate then
failed *inside the runner*.

For identical seeds and the identical 2,400-step budget, the same gate function returned:

| seed | standalone harness | inside the runner |
|---|---:|---:|
| 42 | 0.01871 | passed |
| 43 | 0.01206 | **0.50883** ✗ |
| 44 | 0.05666 | **2.14122** ✗ |

**Three of three pass in the harness; one of three in the runner.** Seed 44's in-runner value is close
to the chance value ln(16) = 2.7726 — not a marginal miss. The harness did not merely give optimistic
numbers, it *inverted the verdict* on two of three seeds. Nothing about either measurement is
incorrect. The runner reaches the gate having constructed the model from the
full experiment configuration and having consumed a different quantity of RNG, so
`copy.deepcopy(model)` starts the memorisation loop from a different initialisation and different
dropout draws. The harness reproduced the gate's *code* and not the gate's *caller*.

This cost two launches: one aborted at the first run, and one that lost **two of three** contrastive
arms mid-experiment after three arms had already trained to completion.

The gate was right and the harness was wrong: the arms it refused belong to an objective independently
measured, at full training duration, to collapse to effective rank 1.71 against the supervised arm's
9.81. A gate that rejects work the harness would have wasted is functioning correctly; the failure was
in reading it from the wrong place.

**Design rule: a liveness gate must be read from the process that will perform the training, never
from a harness that reconstructs the setup.** A harness is useful for developing a gate and is
evidence about the harness. If a gate is expensive enough that one is tempted to rehearse it
elsewhere, the correct response is to make the real gate cheaper — not to move it.

*Provenance: `d1_relaunch`, `queue_size_implicates_the_key_set`.*

## What this implies for gate design

0. **Read the gate from the process that will train.** Instance 4 is the cheapest of the four
   to avoid and cost the most: two launches.
1. **State what the gate certifies, in the gate.** "This model can memorise 16 patients against a
   frozen key set" is a different sentence from "this objective will learn", and only the first is
   supported. Our gate's own failure text said "expected a practically memorised actual-model
   objective" — accurate, and read as more than it said.
2. **A gate's regime must match the run's in every dimension the run's failure can live in.** Ours
   differed in three simultaneously: duration (800 vs 40 epochs), problem size (16 vs 3,118 patients),
   and key-set dynamics (frozen vs live). Each hid a distinct failure — and a fourth hid in the gap
   between the gate and a faithful re-implementation of it.
3. **Simplifications made to defeat one pathology are the first place to look for the next.** Each
   of instances 1–3 came from a change that was correct on its own terms.
4. **A gate that cannot fail on the training pathology should not be quoted as evidence about it.**
   The strongest version: report alongside each gate the dynamics it removes, so a reader can see
   what it cannot see.

## Honest limits of this section

- Instance 3's *cause* is established only as far as "the key set, not the objective's weighting".
  The specific mechanism is open: MoCo's staleness account was predicted, tested and **falsified**
  here, and no replacement has been confirmed. The proposed alternative — that a query-written queue
  permits query/key co-movement — is consistent with the measurements but has not itself been tested
  against a control.
- Whether the momentum fix holds beyond 100 steps, or at 40 epochs, is not measured.
- All four instances come from one objective on one architecture. We do not claim a general rate at
  which liveness gates mislead; we claim the failure mode is real, recurred three times in two days
  once we looked for it, and is cheap to test for.
- The effective-rank measurements are participation ratios of centred singular values on held-out
  patients; they are descriptive, not a certified statistic, and are reported here without inference.

## Not measured

- Whether a momentum key encoder resolves the collapse (implementation not begun).
- Whether the collapse reproduces on a second architecture or a second cohort.
- Whether any of the three instances would have been caught by an existing published gate design; we
  did not survey gate designs.
