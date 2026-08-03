## 2026-08-03 17:00 UTC — Pre-registering the D1 readout: the default target set puts `programme_only`'s own supervision on the exam

**Logged:** 2026-08-03 17:00 UTC. **How obtained:** `~/e0_run/d1_v1/D1_LAUNCH_PLAN.json` and
`frozen_rna_targets.npz` on the A100 (`150.136.45.194`). **Written before any D1 result exists** —
training is at epoch 2 of 40 — precisely so that the choice of primary readout cannot be, or look
like, a reaction to the numbers.

### Technical

`frozen_rna_targets.npz` groups its targets as:

| group | targets |
|---|---|
| `hallmark_in_training` | **50** |
| `heldout_pathway` | 24 |
| `immune_tme` | 8 |
| `tumour_state` | 8 |
| `random_control` | 90 |

`d2_compare` with `--target-groups` omitted scores *every non-control target* — all 90 of
`hallmark_in_training + heldout_pathway + immune_tme + tumour_state`.

D1's `programme_only` arm is launched with **no** `--programme-targets`, so it trains on
`data.hallmark`: the 50 curated Hallmark programmes. Those are exactly the 50 targets labelled
`hallmark_in_training`.

**So the unrestricted D1 readout would score one arm on its own supervision for 50 of its 90
targets** — 56% of the exam. This is the identical defect that made D2's headline uninterpretable
(Stage-1 audit, item A2: *"50 of the 90 scored targets were `hallmark_in_training` — one arm's own
supervision"*), and for D1 it bites harder, because D1's entire question is whether programme
supervision helps or hurts the molecular channel. Scoring that question partly on the supervision
itself does not measure a channel; it measures memorisation of the training target, and it does so in
the direction that favours `programme_only`.

**Pre-registered readout plan, in order of authority:**

1. **Primary — stratified.** `--target-groups heldout_pathway immune_tme tumour_state`, the 40
   targets *neither* arm trained on. This is the D1 result. The pre-registered prediction recorded in
   `d1_pair_manifest` is `programme_free >= programme_only` on the held-out molecular channel; if
   `programme_only` wins here, the manifest's own instruction is to escalate, not to reframe.
2. **Negative control.** `--target-groups random_control`, 90 targets. Both arms must sit at chance.
   Per the audit, a channel on random controls voids every number on the project, not just D1's.
3. **Secondary, reported but not headlined — unrestricted.** All 90 non-control targets, quoted only
   alongside the explicit statement that 50 of them are `programme_only`'s training targets. A gap
   that exists here and vanishes in (1) is evidence about contamination, not about biology.
4. Seed agreement across 42/43/44, effective rank per seed per arm (**reported, not interpreted** —
   blocker 5), and the seed-42 re-export reproducibility check.

CALIBRA's G4 controls (`--require-rna-positive-control --require-channel-gates`) run inside
`run_d1` and must pass before any of the above is quoted.

**Implemented, also before the numbers exist**, as `v2/research/rebase/d1_audit.py` (added
2026-08-03 18:45 while training was at epoch 17/40). It runs the two readouts this run's pipeline
will not produce, then executes A1–A6 mechanically and writes `D1_AUDIT.json` and
`D1_READOUT_INDEX.json` into the run root. Two deliberate properties:

* **It refuses to compute any readout if A1 fails.** An incomplete set of runs exits with
  *"do not compare incomplete arms"* rather than producing a partial comparison.
* **A3 is not a test against zero.** The null for a 16-component canonical correlation is *not* zero
  — CCA is biased upward — so "both arms at chance" cannot be tested that way. The testable statement
  is that random targets must score **below** real ones for both arms and every seed, and even that
  is marked as necessary-but-not-sufficient: the margin needs a human read, because a control sitting
  just below the real score is an alarm, not a pass.

A5 is likewise computed and reported but explicitly **not** interpreted, per blocker 5.

### In plain terms

Half of the standard exam for this experiment is made up of the exact questions one of the two
students was taught the answers to. That student will score better on those questions whether or not
it understands anything, so the honest exam is the other forty questions, which neither student saw.
That is being fixed as the headline result now, before any marks are in, so that it cannot later look
like the marks were chosen to suit an argument.

### Meaning for the claim

Nothing yet — D1 has produced no numbers. This fixes which number will count as D1's answer. It also
means the D1 headline will rest on 40 targets rather than 90, which is a smaller and noisier basis,
and the confidence intervals should be expected to be correspondingly wider. That is the honest cost
of removing the contamination, and it is better paid than hidden.

### Files / commits

- `~/e0_run/d1_v1/D1_LAUNCH_PLAN.json`, `D1_PAIR_MANIFEST.json`
- `v2/research/rebase/d2_compare.py` (`--target-groups`), `frozen_rna_targets.npz`
- NOTEBOOK.md "Post-D2 sequence" Stage 1, items A2 and A3
