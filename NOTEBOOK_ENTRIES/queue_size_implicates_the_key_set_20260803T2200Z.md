## 2026-08-03 22:00 UTC — Collapse severity tracks queue size: the key set is implicated, the objective is provisionally exonerated

**Logged:** 2026-08-03 22:00 UTC. **How obtained:** `~/ws_d1/decorr_causal.py` on the A100
(`150.136.45.194`), varying `PairedBiologyMemoryBank(capacity=...)` with everything else held fixed.
All runs from the **same initialisation** (verified: effective rank 67.55 at step 0), real streaming
train batches, live queue, bf16 + grad-clip 1.0 as `train_epoch`, lr 2e-4 as D1 runs. Centred
effective rank on a fixed 256-patient held-out probe. Logs `~/e0_run/d1_diag/qsweep_*`.

### Technical

Two independent series, each varying **only** queue capacity:

| capacity | `decorrelation 0.04` step 50 / 100 | `decorrelation 0.0` step 50 / 100 |
|---|---|---|
| 4096 *(as D1 runs)* | 4.08 / **1.95** | 2.62 / **2.16** |
| 512 | 5.05 / **2.55** | 3.62 / **3.64** |
| 64 | 9.49 / **6.12** | 6.50 / **5.14** |

**Monotonic in queue size, in both series, at both step counts — four independent confirmations.**
Shrinking the queue from 4096 to 64 raises effective rank 3.1× (decorrelation 0.04) and 2.4×
(decorrelation 0.0) at step 100.

This is the predicted signature. The pathology tracks the **key set**, not the objective's weighting —
consistent with every earlier observation: G2.6 passes with a *frozen* 64-key queue and reaches
effective rank 5.81; no setting of `decorrelation` or `biology_full_consistency`, including both at
zero, changed the outcome at capacity 4096.

The hypothesis this supports is the one MoCo (He et al., 2020) was introduced to address:
`PairedBiologyMemoryBank.update()` enqueues keys from the *current* encoder every step with no
momentum (`v2/training.py`), so a 4096-key queue spans a long stretch of encoder drift and its keys
are mutually inconsistent. The cheapest escape from an inconsistent key set is to stop distinguishing
anything.

**What this does NOT yet establish, and it matters.** Capacity changes *two* things at once: how
**inconsistent** the keys are, and how **many** negatives there are (~4300 versus ~64 effective
negatives). A smaller queue is both more self-consistent *and* an easier discrimination problem. This
sweep therefore implicates the key set but **cannot separate inconsistency from negative count**.

The disambiguation is the momentum-encoder test, precisely because it holds capacity at 4096 — same
number of negatives — while making the keys consistent. If collapse resolves there, inconsistency is
the cause. If it does not, negative count is, and the remedy is different (a smaller queue, or a
harder-negative-aware loss). That test is **not yet implemented**; this entry exists so the
implicating evidence and its limit are on record before any implementation is committed to.

### Concurrent: D1-A has partially failed, and for the familiar reason

`programme_free` seed 43 failed its **in-runner** G2.6 at **contrastive 0.50883** (reduction 0.890,
consistency 0.00488) against the unchanged ≤ 0.10 threshold. State of the six runs:

| run | outcome |
|---|---|
| `d1_p_seed42` | 40 epochs ✓ |
| `d1_f_seed42` | 40 epochs ✓ |
| `d1_p_seed43` | 40 epochs ✓ |
| `d1_f_seed43` | **gate failed, 0 epochs** |
| `d1_p_seed44`, `d1_f_seed44` | in flight |

`run_d1` will raise on that non-zero return code after the remaining runs finish, so it will **not**
produce exports, CALIBRA or the bootstrap. Seed 42 nonetheless retains a **complete paired arm at
epoch 40**, which is what the control needs: A5's effective-rank comparison can still be made on it.

The discrepancy is worth recording precisely. For seed 43's `programme_free` arm, the same gate at the
same 2400-step budget read **0.0047 in my standalone harness** and **0.50883 in the runner** — two
orders of magnitude apart, because the runner reaches the gate with a different model configuration
and RNG history. I flagged this as a methodological point when the first launch aborted; it has now
cost a second run. The standalone harness is not a proxy for the gate, and should not be used to
decide whether to launch.

### In plain terms

Shrinking the pool of "other patients" the model compares against — and changing nothing else — makes
the collapse two to three times less severe, every time, in two separate configurations. That points
at the comparison pool rather than at the loss function, which is where the last two days of fixes
were aimed.

The honest gap: making the pool smaller also makes the task easier, and this test cannot tell those
apart. The next test can, because it keeps the pool exactly as large while making its contents
mutually consistent.

Separately, the control run lost one of its six arms to the same health check that has been the
subject of all of today's repairs — the check failed the arm it had passed in rehearsal.

### Meaning for the claim

The objective is provisionally exonerated and the queue implicated, but the mechanism is not settled
and no fix is yet warranted on this evidence alone. D1-B remains correctly held. D1-A survives as a
partial control: one complete seed-42 pair at epoch 40.

### Files / commits

- `~/e0_run/d1_diag/qsweep_d{0.04,0.0}_cap{512,64}.log`, `decorr_causal_*.log`, `causal_d*_c*.log`
- `v2/research/rebase/d1_collapse_causal_test.py`, `d1_geometry_probe.py`
- `v2/training.py` `PairedBiologyMemoryBank.update` — the no-momentum enqueue
- `~/e0_run/d1_v1_launch.log`
