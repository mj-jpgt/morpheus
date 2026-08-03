## 2026-08-04 07:00 UTC — G2.6 is not reproducible: eight identical runs span 650×, and 2 of 8 fail

**Logged:** 2026-08-04 07:00 UTC. **How obtained:** `~/ws_d1/gate_check.py` on the A100
(`150.136.45.194`), run **eight times with identical inputs** — same seed (42), same cohort, same
split, same schedule, same 2,400-step budget, same code. Logs `~/e0_run/d1_diag/gatevar_*.log`.

### Technical

`final_biology_contrastive`, sorted, against the unchanged ≤ 0.10 threshold:

```
0.00859  0.01076  0.01770  0.02019  0.02407  0.03266  0.38009  5.58511
   PASS     PASS     PASS     PASS     PASS     PASS     FAIL     FAIL
```

| | |
|---|---|
| runs | 8, identical inputs |
| pass rate | **6/8 = 75%** |
| min / median / max | 0.00859 / 0.02213 / **5.58511** |
| range | **650×** |
| threshold | 0.10 (unchanged) |

**The gate is not reproducible at fixed seed.** The only source of variation is non-determinism in the
GPU kernels, compounded over 2,400 optimisation steps.

**The distribution is bimodal, not noisy.** Six runs cluster tightly at 0.009–0.033 — a spread of
under 4× and all comfortably passing. The two failures sit at 0.380 and **5.585**, the latter twice the
chance value of ln 16 = 2.7726. This is not measurement scatter around a mean that happens to straddle
a threshold; it is an optimisation that usually converges and occasionally diverges outright. When it
diverges, the gate reports a collapsed model that is not collapsed.

**Consequence for any experiment gated on it.** Each contrastive arm rolls this independently:

| arms gated | P(all pass) |
|---|---|
| 1 | 0.75 |
| 2 | 0.56 |
| **3** | **0.42** |

A three-seed D1 has a **58% chance of losing at least one contrastive arm** to gate non-determinism
alone, irrespective of whether the objective works. That matches the observed history exactly:

| run | F42 | F43 | F44 |
|---|---|---|---|
| D1-A | pass | fail (0.509) | fail (2.141) |
| D1-B | fail (0.892) | **pass** | not reached |

Two of the three comparable arms flipped verdict between runs with nothing changed that the gate can
see — the momentum encoder is inert inside the gate, verified in source. Previously attributed to a
harness/runner difference; that difference is real and separately evidenced, but **this** result shows
the gate is also unreliable against itself.

**What this does not say.** It does not say the gate is wrong about the objective. Every D1-A arm the
gate refused belonged to a configuration independently measured to collapse to effective rank 1.71 at
full training duration. It says the gate cannot be *relied upon* to reach the same verdict twice, so a
single reading of it — pass or fail — is weak evidence either way.

### In plain terms

The check that decides whether a run may start was given the same model, the same data and the same
settings eight times. Six times it said "healthy, comfortably". Twice it said "dead" — once
emphatically, reporting a score twice as bad as random guessing. Nothing differed between those runs
except the order in which the graphics card happened to add up its numbers.

So when this check fails an arm, that is now roughly a one-in-four event that carries little
information. And a three-repeat experiment has a better-than-even chance of losing an arm to it for no
reason at all.

### Meaning for the claim

The two D1-A contrastive arms rejected at 0.509 and 2.141, and the D1-B arm rejected at 0.892, cannot
be interpreted as evidence about the objective. Nor can the arms that passed be interpreted as
evidence that it works — which is consistent with what we independently established, that a passing
gate did not predict training health.

This is a fifth instance for `paper/LIVENESS_GATE_DESIGN.md`, and a different kind from the first
four: those were about a gate measuring the wrong thing, this is about a gate failing to measure the
same thing twice. Both matter, and the second is cheaper to test for than the first — eight repeated
runs is a few GPU-hours and would have saved considerably more here.

**No gate has been changed.** The replacement rank probe is being measured for the same defect before
anything is proposed; a replacement with this property would be no improvement.

### Files / commits

- `~/e0_run/d1_diag/gatevar_1..8.log`
- `v2/runner.py` `_overfit_programme_free_contrastive`
- Prior: `d1b_blocked_gate_does_not_exercise_the_fix`, `d1a_control_complete_and_gate_fails_2of3_in_runner`
