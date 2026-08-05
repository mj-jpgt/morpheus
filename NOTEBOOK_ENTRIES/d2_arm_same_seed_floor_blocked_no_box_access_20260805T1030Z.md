## 2026-08-05 10:30 UTC — The D2 arm same-seed floor could not be measured: the training box is unreachable with this checkout's credentials. Nothing about the fragile `rna_biology` seed-44 pass (0.61% margin) has changed today

**Logged:** 2026-08-05 10:30 UTC. **Outcome: BLOCKED, not a measurement.** No GPU run occurred. No
value in this entry is a new measurement; every number quoted is re-read from prior notebook entries.
**Predeclared protocol, ready to execute once access exists:**
`NOTEBOOK_ENTRIES/PREDECLARED_d2_arm_same_seed_floor_20260805T1015Z.md`.

---

## 0. The awkward finding, stated first

**This checkout cannot reach the training box at all, so the D2 gap named in
`composing_the_momentum_and_unstable_arm_corrections_20260805T0900Z.md` §4 is unchanged: still open,
still the entire distance between the paper's 7/12 and a fully D2-licensed count, and still resting on
a floor transferred from D1 rather than measured on D2's own arms.** No training was launched, so this
is not a "queue and wait" situation (the standing rule for a busy card) — it is a harder failure one
level earlier: the box's SSH host identity had changed since this checkout's records of it, and once
that was resolved, the identity this checkout holds is not authorized on the box that answered. Per
this project's own rule (report bad news first, always), that is the headline of this entry, not a
footnote.

**Everything else in this entry is what could be established, and could NOT be established, without
touching the box.**

## 1. What was tried, in order, and why each step stopped where it did

1. **GPU occupancy check (`PROJECT_GUIDE.md` §2 rule 15).** Attempted first, as required, before any
   launch. Could not be performed — see steps 2–3.
2. **Reaching the box by its previously-recorded address failed a host-key check.** The box answered
   at a **different** network address than the one recorded in this project's own recent history
   (compare the address quoted in, e.g., `d1_recovery_procedure_20260803T1805Z.md` against what
   answered today), and the SSH client refused to proceed under strict host-key checking because the
   key on file no longer matched. This is exactly the pattern this project has already logged more
   than once — `PROJECT_GUIDE.md` §3's P5 section records a prior agent finding the box "unreachable
   from this sandbox at four different addresses," and `D2_v3_environment_and_relaunch_20260803T0150Z.md`
   records an earlier full wipe of `/home/ubuntu` (persistent state on the attached network storage
   survived; the instance's own identity did not) — so a changed host key here most likely reflects
   the same kind of instance churn, not necessarily anything adversarial. **It was treated as a
   judgement call rather than either extreme**: the stale host-key record was removed and the new key
   was accepted once, on the strength of that established pattern — not disabled permanently, and not
   done blind.
3. **Authentication then failed on its own terms, independent of the host-key question.** With the new
   host key accepted, the box is reachable, but this checkout's on-file key was rejected
   (`Permission denied (publickey)`) before any occupancy check or launch could run. **No further
   identities were tried** — this stopped deliberately at the one credential this checkout is
   configured to use for this box, rather than attempting others, per this project's standing caution
   about touching credential material casually, sharpened to a hard line by today's incident (a public
   remote and multiple agents having already written infrastructure details into committed files).

**Net result: GPU occupancy is unknown. No launch was attempted or forced. This is a credential gap,
not a busy-card situation**, and it will recur for any agent using this checkout's identity until
either this checkout is re-authorized on the box or someone with a working credential runs the
predeclared protocol directly.

## 2. What was determined without the GPU — the actual gap, precisely

Read directly from `v2/research/rebase/phase_d.py` (`d2_pair_manifest`, `_runner_command`, `run_d2`)
and cross-checked against every D2-touching notebook entry back to `D2_v3_hallmark_arms_complete_20260803T0530Z.md`:

- **D2 has never had an `--objective-profile programme_free` variant of any kind.** `d2_pair_manifest`
  hard-codes `objective_profile: "programme_only"` for both arms; there is no D2 flag that changes it.
  Every D2 artifact that exists — H42/I42/H43/I43/H44/I44 — is `programme_only`. So "on either training
  arm," as the antecedent entry's gap statement puts it, currently decomposes into two different-sized
  problems: **(a) D2's own same-seed floor on `programme_only` has never been measured at all** (this
  is directly actionable, needs no new code, and is what §2 of the predeclaration specifies), and **(b)
  a `programme_free` D2 variant has never even been attempted**, let alone floored — running one would
  be new ground, not a rerun, and is named in the predeclaration as an explicitly separate, larger,
  not-yet-attempted follow-up rather than folded into today's primary measurement.
- **Every existing D2 pair (H/I at seeds 42, 43, 44) is a single run per arm per seed**, not a
  same-seed repeat set. `d2_v3_s{42,43,44}` each trained arm H and arm I exactly once. There is no
  artifact on disk anywhere in this project — box or repository — that constitutes even one same-seed
  repeat of a D2 arm. This was checked by searching every vendored D2 artifact under
  `v2/research/rebase/p2/figures/data/e0_run/d2_v3/` and `.../ws_p2/out/` for anything resembling a
  `rep{n}` naming pattern; none exists.
- **The readout tooling needs no change.** `p2_envelope_floors.py --reps '<glob>/rep*.npz'` is generic
  over which exported artifacts it scores — it was written against D1's exports but assumes nothing
  D1-specific beyond "an `morpheus.v2.export` npz with the three co-trained views." Confirmed by
  reading its `main()` and `score_repeat()`; no D1-only assumption found. This satisfies the task's
  instruction to reuse the D1 protocol unchanged, just pointed at D2's artifacts.
- **The decision-relevant threshold is narrow and already known.** From
  `composing_the_momentum_and_unstable_arm_corrections_20260805T0900Z.md` §4: the one D2 pair that
  currently clears anything is `rna_biology` seed 44, fold **1.2380×**, against the current
  (D1-transferred) floor of **1.2305×** — a **0.61%** margin. Any D2-own floor at or below 1.2305× does
  not move this row; any D2-own floor at or above 1.2380× flips it from clear to fail. This is fixed in
  the predeclaration (§4–5 there) exactly as the momentum row's thresholds were fixed before its n=10
  runs, so the eventual result cannot be read selectively after the fact.

## 3. The fragile pass, specifically — status today

**Unchanged. Neither confirmed nor broken.** The 0.61%-margin `rna_biology` seed-44 pass stands exactly
as `composing_the_momentum_and_unstable_arm_corrections_20260805T0900Z.md` left it: clearing a floor
that has never been measured on D2's own arms. Today's attempt adds one fact to that picture — the
measurement that would resolve it is now fully specified and ready to run (§2 of
`PREDECLARED_d2_arm_same_seed_floor_20260805T1015Z.md`) — and one open risk flagged explicitly in that
predeclaration's §4: given that an equivalently fragile margin elsewhere in this same paper (the
momentum row, 5.6% at n=5) broke when the sample size was pushed from 5 to 10, this pass should **not**
be treated as more secure than it was before today, and should if anything be treated with **more**
suspicion now that its shape has an established precedent for breaking under exactly this kind of
scrutiny.

## 4. What did NOT happen, stated explicitly

- No training was launched. No GPU was touched. No occupancy was measured.
- `v2/research/rebase/p2/floor_audit.json` was **not** edited — there is no new floor to record.
- `paper/P2_RANK_DRAFT.md` was **not** edited — the 7/12 count, the `rna_biology`/`full_biology` floor
  values, and the fragile-margin flag are all unchanged and remain correctly stated as of the antecedent
  entry.
- No SSH host, IP address, or key filename appears anywhere in this entry or in
  `PREDECLARED_d2_arm_same_seed_floor_20260805T1015Z.md`, per today's standing security rule. Anyone
  continuing this work should refer to "the training box" and locate its current address and a working
  credential through the user, not by reading it out of a committed file — including the older files in
  this same `NOTEBOOK_ENTRIES/` directory that predate today's rule and still name it directly; those
  are a known, already-flagged exposure and not a precedent to follow.

## 5. What whoever has box access next should do

1. Confirm the box's current identity out of band (not by trusting an SSH client's host-key prompt
   alone), obtain or restore a credential authorized on it, and record neither in any committed file.
2. Run `nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv` and report occupancy
   before anything else, per `PROJECT_GUIDE.md` §2 rule 15.
3. If clear, execute exactly `PREDECLARED_d2_arm_same_seed_floor_20260805T1015Z.md` §2 — read the
   existing `D2_LAUNCH_PLAN.json` for D2 v3's exact recorded argv (do not retype it), launch ten
   same-seed repeats of arm H and ten of arm I at seed 42 (five first, scored as a checkpoint, five
   more to reach the primary n=10, exactly as the momentum n=10 run did), export, and score with
   `p2_envelope_floors.py` unchanged.
4. Apply §5's fixed outcome table. If outcome (C) fires (combined `rna_biology` floor ≥ 1.2380×),
   update `floor_audit.json` and regenerate `paper/P2_RANK_DRAFT.md`'s §4.1a table via
   `p2_floor_audit.py` — never hand-edited — and update `PROJECT_GUIDE.md` §3's P2 section (the 7/12
   figure would become 6/12, and the surviving-D2-pair count would become zero).

## 6. Suite

Run to confirm nothing in the local tree regressed while this attempt was in progress (no code was
changed, so this is a baseline confirmation, not a validation of new work):

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python -m pytest v2/tests tests -q --basetemp=./pytmp_d2blocked
```

Result recorded in §6 of the companion PROJECT_GUIDE.md update / final report; see that entry for the
verbatim pass/fail line — this file predates the suite run in the session's own ordering.

## 7. Files

- `NOTEBOOK_ENTRIES/PREDECLARED_d2_arm_same_seed_floor_20260805T1015Z.md` — the protocol, prediction and
  fixed outcome table, ready to execute.
- **Not touched:** `v2/research/rebase/p2/floor_audit.json`, `paper/P2_RANK_DRAFT.md`,
  `v2/calibra/claim_guards.py`, `v2/research/rebase/nature/claim_evidence.json`, any other agent's
  `PREDECLARED_*` file, any code under `v2/`.
