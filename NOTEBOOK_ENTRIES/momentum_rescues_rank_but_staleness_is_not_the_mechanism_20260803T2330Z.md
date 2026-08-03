## 2026-08-03 23:30 UTC — A momentum key encoder rescues rank at fixed capacity — but the staleness explanation is falsified by our own measurement

**Logged:** 2026-08-03 23:30 UTC. **How obtained:** `~/ws_d1/momentum_test.py` on the A100
(`150.136.45.194`). Capacity held at **4096 throughout** — same negative count in every arm — varying
only which encoder writes the keys. All arms from the same initialisation (verified: effective rank
67.55 at step 0), real streaming batches, bf16 + grad-clip 1.0 as `train_epoch`, lr 2e-4, decorrelation
0.04 as D1 runs. Logs `~/e0_run/d1_diag/mom_*`.

### Technical

**The fix works.** Capacity fixed at 4096, so negative count is identical across arms:

| key encoder | eff-rank step 50 | step 100 |
|---|---:|---:|
| none — query encoder writes keys *(current behaviour)* | 4.07 | **2.58** |
| EMA, m = 0.99 | 8.63 | **6.65** |
| EMA, m = 0.999 | 9.33 | **6.89** |

A momentum key encoder raises effective rank **2.6×** at step 100 with the negative count unchanged.
That is the discriminating comparison the capacity sweep could not make, and it comes out in favour of
the key encoder rather than the number of negatives.

**The proposed mechanism does not survive contact with the measurement.** Two independent
observations kill it.

*First, the queue is not stale.* Batch size is 214 patients (measured: 4309.4 effective negatives
against a 4096 queue). A 4096-slot queue therefore **turns over completely every 19 steps**. The
oldest key in the queue is at most nineteen steps old. Classic MoCo staleness — keys written hundreds
or thousands of steps ago by a substantially different encoder — is not the regime we are in.

*Second, key-to-encoder agreement does not predict rank.* Measured directly, as
cosine(stored key, fresh re-encoding of the same patient by the current query encoder):

| key encoder | agreement s50 | s100 | eff-rank s100 |
|---|---:|---:|---:|
| none | 0.409 | 0.427 | 2.58 |
| m = 0.99 | 0.620 | **0.908** | 6.65 |
| m = 0.999 | 0.181 | 0.441 | **6.89** |

The arm with the **best** agreement (0.908) does **not** have the best rank, and the arm with the best
rank (m = 0.999) has agreement 0.441 — statistically indistinguishable from the *failing* arm's 0.427.
The predicted signature was that momentum would flatten this curve; it did not, and rank recovered
anyway. That was the pre-declared falsification condition and it has fired.

*Third, and pointing the same way:* the capacity sweep's best arm was capacity **64**, where the queue
is entirely overwritten **every single step** (64 slots, 214 patients per batch) — the freshest
possible keys — and it produced the highest rank of the sweep (6.17). If staleness were the mechanism,
that arm should have been the healthiest by construction and it was; but it is also the arm with the
fewest negatives, so it cannot separate the two either.

**What the evidence now supports.** Two distinct effects, neither of which is key age:

1. **Decoupling the key encoder from the query encoder** (momentum, at fixed capacity). Without it the
   keys are produced by essentially the same weights as the queries one step earlier, so the loss can
   be reduced by moving queries and keys *together* — the degenerate direction. A slowly-moving key
   encoder is a target the query encoder cannot drag along, which removes that escape. This is
   *related* to MoCo's motivation but is not its stated mechanism: MoCo argues from key-to-key
   inconsistency across a long queue, and our queue is 19 steps deep.
2. **Fewer negatives** (capacity sweep, at fixed key encoder). Still confounded, and still unexplained.

**Limitation of the staleness metric, stated so it is not over-read.** What was measured is key-to-*query*
lag, not key-to-*key* inconsistency. For m = 0.999 a low cosine is expected and reflects uniform lag of
the whole key set, which is harmless — MoCo's concern is *spread* across the key set, not mean offset.
With a 19-step turnover, key-to-key spread is small in every arm, which is itself the argument that
inconsistency is not what is being fixed. A cleaner metric would be the variance of key-to-fresh
cosine across the age range within a single queue; that is not measured here.

### In plain terms

Giving the negative examples their own slow-moving encoder more than doubles the health of the
representation, without changing how many negatives there are. So the fix is real and it is about
where the negatives come from, not how many there are.

The reason we thought it would work is wrong. The standard story is that a long queue holds
out-of-date examples. Ours does not: it completely refills every nineteen steps, and the setting that
worked best holds examples that are *more* out of date than the setting that failed. The likelier
explanation is simpler — when the same network produces both the questions and the answers, it can
cheat by moving both together, and a separate slow-moving network for the answers takes that away.

### Meaning for the claim

A fix that raises effective rank 2.6× at cohort scale exists and is contained. But it should not be
landed with the MoCo justification attached, because that justification is measurably not what is
happening here, and a wrong mechanism in the record is worse than an unexplained fix. The paper
section's instance 3 must have its mechanism paragraph rewritten accordingly — the *observation* (the
gate froze the queue and thereby removed the dynamic that fails) stands unchanged; the *explanation*
does not.

Not yet established: whether the momentum arms hold rank beyond 100 steps, whether the effect survives
at 40 epochs, and which of the two effects dominates. None of this licenses a D1-B launch yet.

### Files / commits

- `~/ws_d1/momentum_test.py`, `~/e0_run/d1_diag/mom_{0,0.99,0.999}_d0.04.log`
- `~/e0_run/d1_diag/qsweep_*.log` — the capacity series
- `v2/training.py` `PairedBiologyMemoryBank.update` — the no-momentum enqueue, unmodified
- Requires revision: `paper/LIVENESS_GATE_DESIGN.md` instance 3, mechanism paragraph
