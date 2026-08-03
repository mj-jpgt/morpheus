# What a liveness gate certifies — and three ways we learned it certifies less than we assumed

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

This is the sharpest of the three, because the removal was deliberate and documented.

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

The pathology is therefore not in the objective's weighting. The leading hypothesis — under test as of
this writing — is the mechanism MoCo (He et al., 2020) was introduced to address: a rapidly-changing
key encoder makes queued keys mutually inconsistent, and the cheapest escape from an inconsistent key
set is to stop distinguishing anything. `PairedBiologyMemoryBank.update()` enqueues keys from the
current encoder every step with no momentum, which is the configuration MoCo contrasts against.

The gate froze the queue to remove a *known* pathology, and in doing so removed the *unknown* one.

*Provenance: `d1b_premise_fails_all_five_arms_collapse`; queue-capacity sweep in progress.*

## Instance 4 — the gate was read from a re-implementation of itself

The sharpest instance, and the cheapest to avoid.

To decide whether the objective was ready to launch, we ran the gate function
(`_overfit_programme_free_contrastive`) from a standalone harness that reconstructed its inputs:
cohort, split, schedule, model. Three seeds passed with margin, so the run was launched. The gate then
failed *inside the runner*.

For the identical seed and the identical 2,400-step budget, the same gate function returned:

| | `programme_free` seed 43 |
|---|---|
| standalone harness | **0.0047** |
| inside the runner | **0.50883** |

Two orders of magnitude, against a threshold of 0.10 — a pass and a decisive fail. Nothing about
either measurement is incorrect. The runner reaches the gate having constructed the model from the
full experiment configuration and having consumed a different quantity of RNG, so
`copy.deepcopy(model)` starts the memorisation loop from a different initialisation and different
dropout draws. The harness reproduced the gate's *code* and not the gate's *caller*.

This cost two launches: one aborted at the first run, and one that lost an arm mid-experiment after
three arms had already trained to completion.

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
   and key-set dynamics (frozen vs live). Each hid a distinct failure.
3. **Simplifications made to defeat one pathology are the first place to look for the next.** Every
   one of the three above came from a change that was correct on its own terms.
4. **A gate that cannot fail on the training pathology should not be quoted as evidence about it.**
   The strongest version: report alongside each gate the dynamics it removes, so a reader can see
   what it cannot see.

## Honest limits of this section

- Instance 3's *cause* is a live hypothesis, not a result. The queue-capacity sweep that would
  implicate or exonerate the key set is running; the momentum-encoder fix is not implemented and not
  tested. If the sweep clears the queue, the objective itself becomes the candidate and this section's
  third instance still stands as stated — the gate did not predict training health — but its
  mechanism paragraph must be rewritten.
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
