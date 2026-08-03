## 2026-08-04 05:00 UTC — D1-B blocked: G2.6 rejects the fixed objective, because the gate cannot see the fix

**Logged:** 2026-08-04 05:00 UTC. **How obtained:** `~/e0_run/d1_v2_launch.log` and code inspection of
`v2/runner.py` / `v2/training.py` on the A100 (`150.136.45.194`). **No gate has been modified. This
entry requests a decision rather than recording one.**

### Technical

D1-B launched with `--biology-key-momentum 0.999` verified present in all three runners' argv and in
the pair manifest. `programme_free` seed 42 then **failed G2.6 at contrastive 0.89237** (reduction
0.813, consistency 0.00336) against the unchanged ≤ 0.10 threshold.

**The gate does not exercise the momentum key encoder at all.** Two independent reasons, both
verified in source:

1. The gate builds its own trainer as
   `V2Trainer(clone, optimiser, clone_schedule, device, gradient_diagnostics_every=0)` —
   `biology_key_momentum` is not passed, so it takes the default of 0.0.
2. Even if it were passed, `update_biology_keys` returns immediately when `freeze_biology_memory` is
   set, and the gate sets it. The gate replays one fixed batch against a **frozen** queue, so there is
   no enqueueing for a key encoder to do.

So G2.6 is measuring the configuration we have just spent three rounds establishing does not work,
and rejecting arms that would now train under the configuration that does. **The gate has become a
false-negative filter on the fixed objective.**

This is the paper's instance 3 recurring one level up. The gate froze the queue to defeat a real
pathology; that freezing is what makes it blind to a fix whose entire mechanism is *who writes the
live queue*. A gate cannot certify a repair to a dynamic it removes.

**A second fact, and it is not about the fix.** Seed 42's `programme_free` arm **passed** this same
gate in D1-A and **fails** it now at 0.892. The gate is momentum-independent by construction, as
above, so the fix cannot be the cause. Nothing in the changed code consumes RNG before the gate runs.
The remaining explanation is that the gate is **not reproducible run to run at fixed seed** — which is
consistent with the standing note that training on this stack is not seed-reproducible, and with the
gate's own history on marginal arms (harness 3/3 pass, runner 1/3, values ranging 0.012 → 2.14 on
nominally identical configurations).

If that is right, then for arms near the threshold the gate's verdict is close to a coin flip, and
D1-B will keep losing contrastive arms for reasons unrelated to whether the objective learns.

### Current state

| run | status |
|---|---|
| `d1_f_seed42` | **G2.6 failed, 0.89237** |
| `d1_f_seed43`, `d1_p_seed42`, `d1_p_seed43` | running |
| seeds 44 | not yet started |

`run_d1` raises on the first non-zero return code, so as it stands D1-B **cannot complete** — no
exports, no CALIBRA, no bootstrap — regardless of what the remaining arms do. It is being left running
for now only because the remaining gate outcomes are cheap information about the gate's variance.

### The decision needed

Three options, none of which I am taking unilaterally, since all touch a gate:

1. **Give the gate the fix it is meant to certify.** Requires a memorisation check whose queue is
   live and momentum-written, which is a substantial redesign — and the queue was frozen for a real
   reason (a live queue on a replayed batch turns the negatives into the queries; that was defect 1).
2. **Accept the gate as a capacity check only and stop gating D1 on it**, recording that it certifies
   "can memorise 16 patients" and no longer matches the training configuration. This weakens a gate,
   which is precisely what I have been instructed not to do without a decision.
3. **Establish the non-reproducibility first** — run the same seed's gate N times and measure the
   spread. If a marginal arm's verdict is a coin flip, that is a property of the gate worth knowing
   before any redesign, and it is cheap: one gate run is ~40 minutes and they parallelise.

My recommendation is **(3) first**, because it is cheap, it is diagnostic rather than corrective, and
both other options are better informed by knowing whether the gate's verdict is stable. It also
directly tests the claim I have just made about run-to-run variance rather than leaving it as an
inference.

### In plain terms

The repair works on the part of training the health check deliberately switches off. So the check
cannot see the repair, is still measuring the broken arrangement, and is refusing to let the repaired
version train.

Separately, the same arm that the check passed a few hours ago now fails it, with nothing in between
that the check can see. That points at the check being unreliable near its threshold rather than at
anything about the repair — which would be worth knowing regardless of what we do next.

### Meaning for the claim

D1-B is blocked and P2 stays blocked behind it. Nothing measured today is invalidated: the momentum
fix's evidence comes from direct training runs, not from the gate, and stands on its own.

### Files / commits

- `~/e0_run/d1_v2_launch.log`, `~/e0_run/d1_v2/D1_PAIR_MANIFEST.json`
- `v2/runner.py` `_overfit_programme_free_contrastive` — builds its trainer without momentum
- `v2/training.py` `update_biology_keys` — early-returns under `freeze_biology_memory`
- `paper/LIVENESS_GATE_DESIGN.md` instance 3
