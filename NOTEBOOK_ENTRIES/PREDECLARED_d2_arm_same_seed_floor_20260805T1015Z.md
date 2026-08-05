# PREDECLARATION — the D2 arms' (Hallmark H, PBS I) own same-seed retraining floor

**Logged: 2026-08-05 10:15 UTC.** Committed **before** any D2 same-seed repeat has been launched and
before any number below exists. Nothing here is measured; every value quoted as already existing is
re-read from a named file that predates this entry. Execution is currently **blocked** (see the
companion entry logged immediately after this one) — this predeclaration is written anyway, per this
project's rule that predeclaration precedes measurement, so that whoever next has GPU access to the
training box can run exactly this protocol without re-deriving it.

---

## 0. Why this exists

`NOTEBOOK_ENTRIES/composing_the_momentum_and_unstable_arm_corrections_20260805T0900Z.md` §4 names the
gap directly: **"No D2 arm (Hallmark or PBS) has ever had its own retraining floor measured, on either
training arm."** Every floor currently used to judge a D2 pair (`rna_biology` 1.2305×, `full_biology`
1.4220×, `wsi_biology` 3.2947×, RankMe D2 3.5484×) is **transferred** from D1's arms
(`programme_only`/`programme_free`), which share the runner's architecture and objective machinery but
were never themselves trained as arm H or arm I — D2's own same-seed variance has never been measured
at all. The one D2 pair that currently clears anything, `rna_biology` seed 44 at fold **1.2380×**
against the transferred floor of **1.2305×**, clears by **0.61%** — the same shape as the
`m=0.999`-over-`m=0.99` momentum row, which cleared at n=5 by 5.6% and broke at n=10
(`NOTEBOOK_ENTRIES/limit2_breaks_at_ten_repeats_20260805T0330Z.md`).

## 1. Repeat count: n=10, not n=5 — the standard this project uses today

The momentum row above is the direct precedent this predeclaration follows: a floor measured at n=5
passed a fragile-margin row; the **same floor, measured at n=10, reversed the verdict**, and the
break was caused by a single one-in-ten repeat, not visible in the first five. Given that this D2
question is the *same shape* of fragile margin (0.61% here vs 5.6% there), n=5 is not an adequate
floor to report as final. **Primary protocol: n=10 same-seed repeats per D2 arm** (20 runs total: 10
of arm H, 10 of arm I, both at one fixed seed). Following the momentum entry's own discipline, the
first 5 of each arm are scored and reported as an interim checkpoint, then the second 5 are added and
the combined 10 is the number that governs the verdict — explicitly checking, as that entry did,
whether the first five was a favourable draw.

## 2. Exactly what will be run, and what will not change

**Protocol: the same envelope-floor protocol used for D1 (`~/chain_retrain_envelope.sh` /
`~/chain_unstable_envelope.sh`), reused unchanged except for which artifact is retrained.** D2 has no
`--objective-profile` knob of its own — `phase_d.run_d2` hard-codes `objective_profile=programme_only`
for both arms (`v2/research/rebase/phase_d.py::d2_pair_manifest`, `_runner_command`) — so **this
predeclaration's primary measurement is D2's floor on the only training arm D2 has ever used,
`programme_only`.** A `programme_free` extension of D2 (running the D2 runner invocation with
`--objective-profile programme_free` substituted, mirroring how D1's unstable-arm floor was obtained)
is named as a distinct, larger follow-up and is **not** part of this predeclaration's primary
measurement — it was never attempted for D2 before today and carries its own risk (D2's
`--fit-programme-legibility` / neighbourhood / supcon losses have never been run under
`programme_free`, and whether that combination is even a supported, gated configuration is unverified).

**Commands, read from the existing plan file as authority, not retyped** — the same discipline
`NOTEBOOK_ENTRIES/d1_recovery_procedure_20260803T1805Z.md` fixed for D1 recovery: read
`~/e0_run/d2_v3/d2_v3_s42/D2_LAUNCH_PLAN.json` (or the equivalent for whichever seed is chosen) and
reuse its exact recorded `argv` for arm H and arm I, changing only `--output-dir` per repeat. If that
file is unavailable, the equivalent command is reconstructed from `phase_d.py`'s own `_runner_command`
and its `d2` argparse defaults (`v2/research/rebase/phase_d.py` lines ~702–725), which for arm H is:

```
python -m morpheus.v2.runner \
  --data-config <the D2 v3 data config> --split-file <the D2 v3 split file> \
  --output-dir <run-root>/d2_envelope/h_seed42_rep$r \
  --objective-profile programme_only --epochs 40 --token-budget 32768 \
  --hidden-dim 512 --layers 4 --heads 8 --learning-rate 2e-4 --weight-decay 1e-2 \
  --decorrelation-weight 0.04 --loss-warmup-epochs 4 --seed 42 --device cuda \
  --fit-development --fixed-final-epoch --fit-programme-legibility \
  --programme-warmup-weight 0.50 --programme-weight 1.0 \
  --programme-neighbourhood-weight 0.20 --programme-supcon-weight 0.20 \
  --separation-weight 0.01 --variance-weight 0.01 --programme-head-dim 256 \
  --gradient-diagnostics-every 25 --pretrain-learning-rate 2e-4 --pretrain-mask-fraction 0.30 \
  --pretrain-view-keep-fraction 0.70 --pretrain-target-dim 128 \
  --expected-development-cancers 11 --expected-heldout-cancers 21 \
  --d2-pair-manifest <the D2 v3 pair manifest> --d2-arm H --d2-analysis-role primary \
  --d2-pbs-components 128 --restrict-to-split
```

and for arm I, identical plus `--d2-arm I` and `--programme-targets <the frozen PBS target artifact>`.

Held fixed and stated so a later reader can check them: **10 repeats per arm** (5 scored as an
interim checkpoint, 5 more added to reach the primary n=10); **one fixed seed** in every repeat of a
given arm (seed 42, matching D1's floor seed and the existing single H42/I42 pair, so the first repeat
of each arm can be cross-checked against the already-exported `d2_h_seed42`/`d2_i_seed42` artifacts as
an anchor, exactly as `d1_recovery_procedure` and `d2_coordinate_system_result` both anchor against a
published reproduction before trusting anything new); 40 epochs; the identical data config, split
file, architecture, optimiser, PBS target artifact (for arm I), pair manifest and export invocation
already on record for D2 v3; **as many repeats concurrent per arm as the card can take without
starving another agent's job** (§3 below governs launch).

**Readout, also unchanged.** `v2/research/rebase/p2/p2_envelope_floors.py --reps
'<run-root>/d2_envelope/h_seed42_rep*.npz' --targets <the frozen RNA target artifact> --output
D2_ENVELOPE_FLOORS_H.json`, and the same invocation with `i_seed42_rep*.npz` for arm I. This module is
generic over which exported artifacts it is pointed at — it does not assume D1's arms, only that the
inputs are `morpheus.v2.export` npz files with the three co-trained views — so **no code change is
needed**, only the `--reps` glob. Every statistic is imported (`calibra.spectral` for R1/R2/R3 and the
channel; `p2_competing_metrics` for RankMe/PR/stable rank/α-ReQ/LiDAR; `numpy.linalg.matrix_rank` for
the hard rank). **No statistic is written inline and no new module defines one.**

**Combining the two D2 arms, and combining with the existing D1-transferred floor.** The convention
already fixed in `p2_probe_floors.combine()` is reused, not re-invented: **`max` of the folds, with
the carrying arm/source named, never a pool and never an average.** Three quantities will exist per
view/statistic after this measurement — D2 arm H's own floor, D2 arm I's own floor, and the
already-recorded D1-transferred floor — and the number that governs any D2 pair's verdict is the `max`
of all three.

## 3. GPU occupancy — checked before anything is launched, and this predeclaration does not launch

Per `PROJECT_GUIDE.md` §2 rule 15, occupancy must be checked by `nvidia-smi` compute-app count and
memory, not a process-name guess, before any of the above is launched, and queued rather than forced
into a saturated card. **As of this entry, occupancy could not be checked at all**, for a reason
prior to occupancy: this checkout could not authenticate to the training box (see the companion entry
logged immediately after this one). No launch has occurred and none will until that is resolved.

## 4. Prediction, written before any run exists

The decision-relevant number is narrow: the fragile pass is `rna_biology` D2 seed 44, fold **1.2380×**,
against the current transferred floor of **1.2305×**. Whatever D2's own floor turns out to be, the
**combined** floor (max of D2's own and the transferred 1.2305×) only matters for this row if it
exceeds 1.2380× — anywhere at or below that, seed 44 still clears (by a shrinking margin as the floor
rises toward 1.2380×, none at all below the current 1.2305×).

| quantity | prediction | interval I would be surprised to fall outside | reasoning |
|---|---|---|---|
| D2 arm H, `rna_biology`, R1 residualised, own floor (n=10) | ≈ 1.15× | 1.02× – 1.6× | arm H trains with the same programme-anchoring character as D1's `programme_only` stable arm (1.019× there), but D2 additionally carries `--fit-programme-legibility` plus neighbourhood/supcon losses D1 never runs, which is new noise with no precedent to size |
| D2 arm I, `rna_biology`, own floor (n=10) | ≈ 1.20× | 1.02× – 1.8× | arm I trains against a 128-column PBS dictionary rather than a fixed Hallmark gene-set target — a higher-dimensional, less-anchored supervision signal that plausibly reproduces less tightly than arm H |
| **combined D2 floor, `rna_biology`** (max of H, I, and 1.2305×) | **≈ 1.20×, straddling the 1.2380× decision boundary** | 1.05× – 1.9× | this is the number this predeclaration is actually about |
| shape | **not confidently predicted either way** | — | the momentum precedent (favourable at n=5, broken at n=10 by a single repeat) argues for genuine caution rather than a confident "survives" call |

**I do not predict the fragile pass survives.** Given today's established base rate on this project —
the momentum row broke under exactly this kind of scrutiny, and the RankMe D2 selection separately
flipped clear→fail once its own floor was measured on both arms — I predict this is closer to a
coin flip than the interval above's point estimate suggests, and I flag in advance that I would not
be surprised by either outcome. The interval is wide because there is **no existing measurement of
D2's own same-seed spread on any view, on any arm** to anchor a tighter prediction; this predeclaration
exists specifically because that anchor does not exist yet.

## 5. What each outcome MEANS — fixed now so it cannot be chosen afterwards

**(A) Combined `rna_biology` floor (D2-own max D1-transferred) ≤ 1.2305×.** D2's own arms are no
noisier than the value already in use. Nothing about the current audit changes; the transferred floor
was already conservative for this row. Seed 44 continues to clear by 0.61%, and the margin is now
licensed by a same-seed measurement made on D2's own arms, not merely inherited. **The fragile pass
survives, and is no longer only a transfer.**

**(B) Combined floor in (1.2305×, 1.2380×).** Seed 44 still clears, but by less than 0.61% — report
the exact new margin, flag it more sharply than the antecedent already does, and do not call it a
clean pass at any margin under roughly 0.3%.

**(C) Combined floor ≥ 1.2380×.** Seed 44 **flips clear → fail**. `rna_biology`'s D2 count falls from
1/6 to 0/6 (all six `rna_biology` D2 pairs then fail; the view's total count of resolvable pairs falls
from 4/6 to 3/6, and the paper's overall 7/12 becomes 6/12 — every one of the six D1 pairs, none of
the six D2 pairs). **This is the outcome the momentum precedent and the RankMe-D2 precedent both make
the single most likely reading, and it will be reported as such plainly, not softened.** `floor_audit.json`
and `paper/P2_RANK_DRAFT.md` would need updating via `p2_floor_audit.py` per this project's
regeneration rule (never hand-edited), and PROJECT_GUIDE.md §3's P2 section would need its "7 of 12"
figure corrected to 6 of 12, with the D2 side of the view-count claim recorded as **zero surviving D2
pairs on any view**.

**(D) `full_biology` and `wsi_biology`, for completeness.** Every D2 pair on these two views already
fails against the transferred floor (0/6 each). A D2-own floor can only raise the bar further on these
views — it cannot flip a fail to a clear — so no outcome on these two views changes any verdict; they
are measured anyway because rule 16 (a per-view-only measurement flatters) applies here exactly as it
did to the original exported-floor predeclaration, and because the *size* of the D2-own floor on these
views is itself informative about how D2's arms compare to D1's, independent of any single row's
verdict.

**(E) One or more of the twenty repeats does not complete 40 epochs, trips the rank tripwire, or fails
a gate.** Per the same rule the `programme_free` predeclaration fixed: the ten (or five) repeats of
that arm are not ten (or five) repeats of one configuration, and **no floor may be computed from the
survivors alone** — that is the selection effect this project's own rule exists to prevent. The outcome
reported in that case is "D2 arm <H|I> does not reliably reach 40 epochs at n=<count>", the completion
count is the result, and the transferred floor stays in force with that caveat attached.

## 6. What this measurement cannot do

- It does not vary the seed; a same-seed floor is deliberately blind to seed variance, which this
  project has repeatedly found larger than the same-seed spread.
- It does not touch `programme_free` for D2 — named in §2 as a distinct, larger follow-up, not
  attempted here.
- It says nothing about the probe block, the in-run training batch, or the 16-patient gate batch for
  D2 — those blocks have never been probed for D2 and remain absent regardless of this measurement's
  outcome.
- n = 10 per arm is n = 10. No interval will be quoted, and "the floor is X×" will not be written as
  "D2's rank varies X×".

## 7. Files

- To be produced (once GPU access exists): `<run-root>/d2_envelope/{h,i}_seed42_rep{1..10}/`,
  `{h,i}_seed42_rep{1..10}.npz`, `D2_ENVELOPE_FLOORS_{H,I}.json`
- Reused unchanged: `v2/research/rebase/p2/p2_envelope_floors.py`, the D2 v3 pair manifest and PBS
  target artifact already on the training box, `phase_d.py`'s own D2 argument defaults as the
  reconstruction fallback if the recorded launch plan cannot be read
- Sources: `NOTEBOOK_ENTRIES/composing_the_momentum_and_unstable_arm_corrections_20260805T0900Z.md`
  §4; `NOTEBOOK_ENTRIES/unstable_arm_exported_floor_measured_20260805T0755Z.md`;
  `NOTEBOOK_ENTRIES/limit2_breaks_at_ten_repeats_20260805T0330Z.md`;
  `NOTEBOOK_ENTRIES/PREDECLARED_unstable_arm_exported_floor_20260805T0045Z.md` (the D1 protocol this
  reuses); `v2/research/rebase/phase_d.py` (`d2_pair_manifest`, `_runner_command`, `run_d2`)
- **Not touched:** `claim_guards.py`, `claim_evidence.json`, any other agent's `PREDECLARED_*`,
  `paper/P2_RANK_DRAFT.md`, `v2/research/rebase/p2/floor_audit.json`.
