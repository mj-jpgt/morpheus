## 2026-08-03 20:30 UTC — D1-B's premise fails: `programme_free` collapses at training scale under **every** setting of both suspected drivers

**Logged:** 2026-08-03 20:30 UTC. **How obtained:** `~/ws_d1/decorr_causal.py` on the A100
(`150.136.45.194`). Five throwaway `programme_free` models, all from the **same initialisation**
(verified: effective rank 67.55 at step 0), real streaming train batches, live 4096-key queue,
bf16 autocast and grad-clip 1.0 exactly as `V2Trainer.train_epoch`, lr 2e-4 as D1 runs. Centred
effective rank measured on a fixed 256-patient held-out probe. Logs `~/e0_run/d1_diag/causal_*`.

### Technical

| `decorrelation` | `biology_full_consistency` | step 50 | 100 | 150 | 200 | 250 |
|---|---|---:|---:|---:|---:|---:|
| 0.04 | 1.0 *(as D1-A runs)* | 4.08 | 1.95 | 2.16 | 1.68 | **1.59** |
| 0.0 | 1.0 *(as D1-B was specified)* | 2.62 | 2.16 | 2.47 | 1.94 | **2.17** |
| 0.04 | 0.1 | 2.99 | 3.43 | — | — | — |
| 0.0 | 0.1 | 2.97 | 2.00 | 2.50 | — | — |
| 0.0 | 0.0 | 2.98 | 1.98 | **1.86** | — | — |

**Every configuration collapses from 67.55 to effective rank ~2 within 150 steps.** Neither
`decorrelation` nor `biology_full_consistency` is the cause, at any setting including zero for both.

**Decorrelation aggravates but does not cause.** At step 250 it is 1.59 with the term versus 2.17
without, and at step 200 the RNA-view mutual cosine is 0.9813 with versus 0.4160 without. So its
removal remains warranted — it makes a bad situation worse — but it is not the fix, and D1-B as
specified would have produced a collapsed `programme_free` arm across all ten runs.

**Consistency at 1.0 is not the second driver either.** Setting it to 0.1, or removing it entirely,
changes nothing material: 1.86 at step 150 with the term fully off. The G2.6 isolation that put it on
par with decorrelation (1.847 vs 0.0034) was measured on the memorisation batch, and does not
transfer.

**The real finding is that G2.6, even repaired, does not predict training health.** The two regimes
differ more than the gate's design assumed:

| | G2.6 memorisation check | real D1 training |
|---|---|---|
| patients | 16, fixed, replayed | 3,118, streaming |
| queue | 64 keys, **frozen** | 4,096 keys, **live** |
| outcome | eff-rank 5.81, retrieval 16/16, contrastive 0.012–0.057 | eff-rank ~2, contrastive drifting toward chance |

Memorising sixteen fixed patients against a frozen queue is achievable; learning a general
representation over three thousand streaming patients against a queue that tracks the model is not.
The gate is a genuine liveness test and its repairs were genuine — but passing it is *necessary and
not sufficient*, and today it was treated as sufficient. That is the same error in a new place: the
gate's regime is not the run's regime, exactly as the memorisation batch's scale was not training's
scale.

**Consequence.** Queue item 2 (D1-B, `--decorrelation-weight 0`) should **not** launch as specified.
The prior superseding note already halted it pending these arms; this entry closes that question in
the negative. The open question is no longer "which regulariser causes this" but "**can this
contrastive objective learn a non-degenerate biology representation at cohort scale at all**", and
nothing measured today answers yes.

Candidate directions, none yet tested at training scale, listed so the next attempt is designed
rather than guessed: the live queue is the largest untested difference from the passing regime
(a frozen or slowly-updated queue is the obvious first probe); the learning rate 2e-4 was inherited
from D2 and never swept for this objective; and the centring fix that rescued the gate is applied
inside the contrastive loss only, not to the states the other training terms see.

### Honest note on process

Two errors of mine in this thread, both corrected, recorded because the project's norms ask for it.
I described the D1-B direction as a coordinator chat message when it is in fact a decision recorded in
`NOTEBOOK.md` Stage 3b by the main session; then, over-correcting, I declared I had fabricated it and
reverted the test rewrite — which that same decision explicitly authorises, by name, and which was
already green. The revert has been reverted and the tree is back to the intended state
(`5d7f8d2`). Net effect on the repository is nil; net effect on this log is two extra commits and
this paragraph. The lesson is narrow and worth keeping: **check `NOTEBOOK.md` for the standing
decision before asserting where an instruction came from — and before undoing work on the strength of
that assertion.**

### In plain terms

The plan was to re-run the experiment with one penalty term switched off. I tested that plan first
instead of trusting it, and it does not work: with the term off, and with the other suspect off too,
and with both off at once, the representation still collapses to about two dimensions within a
hundred and fifty steps. Ten runs would have produced ten collapsed models.

The deeper problem is that the health check we spent the day repairing asks the model to memorise
sixteen patients, and it can. Real training asks it to represent three thousand, and it cannot. Those
turn out to be different questions, so passing the check does not license the run.

### Meaning for the claim

The biology-head-without-programme-supervision arm still does not exist in usable form, and the route
planned for obtaining it is now measured not to work. D1-A remains valuable as the control that
documents the collapse at epoch 40 with A5 quantifying it. Nothing about P2's primary evidence should
be considered unblocked.

### Files / commits

- `~/ws_d1/decorr_causal.py`, `~/e0_run/d1_diag/decorr_causal_*.log`, `causal_d*_c*.log`
- `v2/research/rebase/d1_geometry_probe.py`
- Supersedes: `NOTEBOOK.md` Stage 3b queue item 2 and its "Decision, 2026-08-03: remove decorrelation"
- Prior: `d1_programme_free_collapsing_in_training_20260803T1930Z.md`
