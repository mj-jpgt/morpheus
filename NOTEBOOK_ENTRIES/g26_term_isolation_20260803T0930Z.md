## 2026-08-03 09:30 UTC — Every term the G2.6 memorisation clone adds beyond the contrastive loss costs it an order of magnitude

**Logged:** 2026-08-03 09:30 UTC. **How obtained:** `~/ws_d1/diag_i.py` and `~/ws_d1/gate_check.py`
on the A100 (`150.136.45.194`). Every row runs the **actual** gate function
(`_overfit_programme_free_contrastive`), not a reimplementation, on the real cohort, seed 42,
800 steps, with `V2LossSchedule.weights` monkey-patched per arm. Logs in `~/e0_run/d1_diag/`.

### Technical

After the centring fix the gate solved its own objective but the graded metric stayed near chance.
Isolating each remaining term against a contrastive-only baseline. All values are the **raw,
uncentred, queue-included `biology_contrastive`** — the number G2.6 thresholds at ≤ 0.10. Chance is
ln 16 = 2.7726.

| `biology_full_consistency` | `decorrelation` | `variance` | raw graded contrastive | |
|---|---|---|---|---|
| 0 | 0 | 0 | **0.00340** | PASS |
| 0.01 | 0 | 0 | **0.04613** | PASS |
| 0.1 | 0 | 0 | **0.08343** | PASS |
| 1.0 | 0 | 0 | 1.84745 | FAIL |
| 0 | 0 | **0.01** | 0.53165 | FAIL |
| 0 | 0.004 | 0.01 | 0.13706 | FAIL |
| 0 | 0.001 | 0.01 | 2.77288 | FAIL |
| 0 | 0.04 | 0.01 | 2.60579 | FAIL |
| 1.0 | 0.04 | 0.01 (**as shipped**) | 2.63086 | FAIL |

Three separate findings.

**1. The variance floor alone costs a factor of 150.** At its shipped weight of 0.01, with nothing
else added, it moves the graded term from 0.00340 to 0.53165 — past a threshold set at 0.10, by a
term contributing at most 6.25e-4 to the loss value. It buys nothing in exchange: a per-dimension
floor cannot raise rank (its own docstring says so, and weight 10.0 was measured still collapsing).

**2. `decorrelation` blocks the gate at every weight tried, and not monotonically.** 0.04 → 2.61,
0.004 → 0.14, 0.001 → 2.78. The non-monotonicity is itself informative: these are unstable
configurations, not a dose-response curve, so no small-but-nonzero weight is safe.

**3. `biology_full_consistency` at weight 1.0 costs three orders of magnitude** (0.0034 → 1.85), and
is fine at ≤ 0.1. Nothing in the code justifies parity with the primary objective; its stated purpose
is only to give the exported `full_biology` state a declared gradient path.

**What was changed, and what deliberately was not.** Both G2.6 memorisation clones — `programme_only`
and `programme_free` — now exclude `decorrelation` **and** the variance floor. This is not new
policy: `_overfit_programme_only_actual` has always excluded decorrelation, on the stated grounds
that it "has a batch-statistics floor unrelated to memorisation". The same reasoning applies to the
variance floor and now has the same measurement behind it. The `programme_free` clone previously kept
decorrelation on, justified by a comment asserting decorrelation was "the only term opposing that
collapse" — **that assertion is now falsified in the opposite direction**, and the code's own
requirement that "the two arms' liveness checks are only comparable evidence if they are run
identically" was being violated by it.

Removing these terms cannot let a collapsed model through the gate, which is the thing to check
before touching a gate: **collapse scores ln 16 = 2.7726 on the graded contrastive term**, so the
contrastive criterion *is* the anti-collapse test. The regularisers were never what caught collapse.

`biology_full_consistency` was kept at 1.0. Instead the *comparison* it makes is now computed after
removing the batch-common direction, exactly as the contrastive term is — raw, it is trivially
satisfied by the shared direction that is 81% of `z_biology`'s norm, so it was rewarding the very
component that drives collapse. That is a stronger requirement than the raw cosine, not a weaker one,
and the gate still independently requires `full_consistency <= 0.02`.

**Training is untouched by all of this.** Both arms keep `decorrelation` 0.04 and `variance` 0.01 at
identical weights in the real objective, so the D1 contrast still measures programme supervision and
nothing else, and `test_both_d1_arms_pair_decorrelation_with_a_variance_floor` is unaffected.

**Current gate status**, real gate, three seeds, 800 steps:

| seed | shipped (before) | after centring + clone repairs |
|---|---|---|
| 42 | 2.63086 FAIL | **0.02345 PASS** |
| 43 | 2.43554 FAIL | 0.17656 FAIL |
| 44 | 2.66053 FAIL | **0.02939 PASS** |

Seed 43 fell from 2.44 to 0.18 but has not cleared 0.10. The centred objective was measured
descending monotonically well past step 800 in the standalone harness — 0.4025 at 400, 0.0835 at 600,
0.0101 at 800, 0.0030 at 3000 — so 800 steps lands mid-descent. Whether seed 43 is still descending
or has become unstable is being measured at 1600 and 2400 steps before any budget is changed; the
distinction matters, because "still descending" justifies a larger budget and "bounced" does not.

### In plain terms

The check that decides whether the model is healthy enough to train was being run with three extra
penalty terms switched on. Each of them, on its own, is enough to make a healthy model look dead —
the smallest of them, contributing well under a thousandth of the loss, moves the score by a factor
of 150 across a line drawn at 0.10. The sister arm of this same experiment had already switched one
of them off years-equivalent ago for exactly this reason; this arm had not, on the strength of a
comment that turns out to state the opposite of what happens.

Turning them off in the check cannot let a broken model through, because a collapsed model scores
2.77 on the number being graded and the line is at 0.10. The model that gets trained is completely
unaffected — the penalties are still there during real training, at identical strength in both arms.

Two of three seeds now pass with large margin. The third has come down from 2.44 to 0.18 and the
question is whether it simply needs more than 800 steps.

### Meaning for the claim

The ≤ 0.10 threshold has not moved and must not. Two independent configurations now clear it by more
than 3×, and the contrastive-only baseline clears it by 30×, so attainability is settled by
construction rather than by argument. What has changed is the *check*, in the direction of measuring
what it claims to measure: memorisation capacity, not the batch-statistics floors of regularisers
that are irrelevant to memorisation and, on measurement, actively opposed to it.

### Files / commits

- `v2/runner.py` — both G2.6 clones now `replace(schedule, decorrelation_after_warmup=0.0, variance_after_warmup=0.0)`
- `v2/training.py` — `biology_full_consistency` compared after removing the common direction
- `v2/research/rebase/phase_d.py` — `--max-parallel` for the six D1 runs
- `~/e0_run/d1_diag/diag_i_*.log`, `gate2_*.log`, `gate3_*.log`, `gateP_*.log`
