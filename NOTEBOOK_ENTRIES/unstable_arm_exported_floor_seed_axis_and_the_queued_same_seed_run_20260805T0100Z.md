## 2026-08-05 01:00 UTC — The arm effect on rank reproducibility **inverts with the view**: on the seed axis `programme_free` is 1.9× noisier than `programme_only` on exactly the two views §4.1b calls usable, and *less* noisy on the one it calls unusable. The predeclared same-seed floor on the unstable arm is queued behind another agent's GPU job

**Logged:** 2026-08-05 01:00 UTC. **Outcome: RESULT (supplementary, CPU) + RUN IN FLIGHT.**
**Predeclared in** `NOTEBOOK_ENTRIES/PREDECLARED_unstable_arm_exported_floor_20260805T0045Z.md`,
committed at `a392c0a` before anything below was computed and before the training was armed.

**Scope warning, first, because this entry contains two different quantities and only one of them is
the predeclared one.** §1 is the **seed axis** (five *different* seeds), which is **not** a retraining
floor and is not the quantity §4.1's 3.295× belongs to. The predeclared quantity — five **same-seed**
`programme_free` retrains at 40 epochs on the exported block — is **still running** and is §3. Nothing
in §1 may be substituted for it, and §2 gives the measured reason why it cannot be.

---

## 0. Why this exists

The completeness audit of 2026-08-04 20:30 filed as its highest-priority open item that **every
exported-block floor in this paper is `programme_only`** — the stable arm — while the probe block's
both-arm measurement showed the collapsed arm carries the floor by **~2×**
(`the_probe_block_has_a_floor_at_last_20260804T1620Z.md` §3). §4.1's **3.295×** is therefore measured
in exactly the configuration the paper's own argument says flatters a floor.

The predeclaration works out where the exposure actually is, and it is **not** where "the headline
floor may double" sounds like it is:

* On `wsi_biology` the six between-arm pairs span 1.004×–**3.246×** and the floor is 3.295×, so
  **0 of 6** are resolvable. **Any larger floor leaves that at 0 of 6.** The `wsi_biology` half of the
  claim cannot be hurt by this measurement.
* The **whole** exposure is the other half: **6 of 6** on `rna_biology` and **6 of 6** on
  `full_biology`, and those two counts rest entirely on floors of **1.019×** and **1.020×** —
  `programme_only`, same-seed, n = 5.

## 1. The awkward finding available today: the arm effect on the floor **inverts with the view**

**What was measured.** `p2_envelope_floors.py`, unchanged, on the ten exported 40-epoch D1 artifacts
that already existed — `programme_free` and `programme_only` at **seeds 42, 43, 44, 45, 46**
(`~/e0_run/d1_v2/artifacts/` and `~/e0_run/d1_seeds4546/artifacts/`). Same module, same statistic
table, same α-ReQ index range, same LiDAR ridge, same targets (`frozen_rna_targets.npz`, the 40
untrained), same `test` partition, same cancer + pooled-TSS residualisation. **CPU only, `nice -n 15`,
thread-capped; the GPU was not touched.** Outputs `~/e0_run/pf_seedaxis/out/SEEDAXIS_{f,p}.json`.

Canonical R1, five seeds, `max/min`:

| view / block | `programme_only` (stable) | `programme_free` (unstable) | unstable ÷ stable |
|---|---:|---:|---:|
| `wsi_biology`, residualised | **2.798×** | **2.151×** | **0.77** |
| `wsi_biology`, raw | 2.657× | 2.108× | 0.79 |
| **`rna_biology`, residualised** | **1.061×** | **1.978×** | **1.86** |
| `rna_biology`, raw | 1.070× | 1.967× | 1.84 |
| **`full_biology`, residualised** | **1.207×** | **2.158×** | **1.79** |
| `full_biology`, raw | 1.186× | 2.111× | 1.78 |
| channel (top-CCA 16) | 1.038× | 1.153× | 1.11 |

**The stable/unstable ordering the probe block established does not hold on the exported block; it
holds on two of the three views and reverses on the third.** On `wsi_biology` — the view §4.1's floor
is measured on and the view every §4 between-arm comparison is read from — the **stable** arm is the
noisier one, by 1.3×. On `rna_biology` and `full_biology` the unstable arm is noisier by **1.86×** and
**1.79×**, and those are precisely the two views whose 1.019× / 1.020× floors carry §4.1b's "12 of 12
resolvable" and §1.3's *"on two other views ... the floor is fifty times smaller and every between-arm
difference clears it."*

**The mechanism is visible in the levels and is not mysterious.** `programme_free`'s three views all
sit at R1 ≈ 6–14 and move together (its per-seed `wsi`/`rna`/`full` folds are 2.15/1.98/2.16 — one
number, three views). `programme_only`'s RNA view is pinned near 27.6–29.3 across every seed while its
WSI view swings 11.1–31.1. **`programme_only`'s RNA view is stable because programme supervision holds
it there; there is no general fact about the RNA view being reproducible.** The paper currently reads
1.019× as a property of the *view*. On this evidence it is at least as much a property of the *arm*,
and the paper has one arm.

**Nothing here is the predeclared floor, and it does not license a verdict.** Five seeds is a strictly
larger perturbation than five same-seed retrains and it is §4.2's axis, not §4.1's.

## 2. And the seed axis does **not** bound the same-seed axis — measured, not assumed

The tempting shortcut is "the seed spread is bigger, so it caps the retraining floor". **On the one
arm where both are measured, it is false in the direction that matters.** `programme_only`,
`wsi_biology`, residualised R1:

| axis | n | floor | shape |
|---|---:|---:|---|
| **same seed**, five retrains (§4.1) | 5 | **3.295×** | bimodal (rep2 at a third; rest 1.028×) |
| **five seeds** (§4.2's axis, measured here at n = 5) | 5 | **2.798×** | not bimodal (rest 2.208×) |

The same-seed floor is **larger** than the five-seed floor on this view, because one same-seed repeat
collapsed and no seed did. On `rna_biology` the order is the other way (1.019× same-seed against
1.061× five-seed). **So §1 cannot be read as an upper bound on §3, in either direction, and this run
cannot be skipped.** That is stated here, before §3 reports, so that a convenient §3 cannot be
rationalised afterwards by pointing at §1.

## 3. The predeclared run — armed, code-identity verified, queued behind another agent

**Status: RUNNING (waiting on the card). No number exists yet.**

* **Protocol.** `~/chain_unstable_envelope.sh` is `~/chain_retrain_envelope.sh` — the script that
  produced 3.295× — with **one flag changed**, `--objective-profile programme_only` →
  `programme_free`. Five repeats, **seed 42 in all five**, 40 epochs, identical data config, split
  file, architecture, optimiser, schedule, `--biology-key-momentum 0.999`, tripwire, and the identical
  `morpheus.v2.export` invocation. Five runs concurrent on one A100, as the stable arm's five were.
  Outputs `~/e0_run/d1_envelope_pf/rep{1..5}.npz`; readout by `d1_envelope_readout.py` **and**
  `p2_envelope_floors.py`, both unchanged, plus a re-scoring of the stable arm's five through the same
  invocation in the same session so the two arms' folds are comparable by construction.
* **Workspace.** `~/ws_uf/morpheus`, built with `git -c core.autocrlf=false -c core.eol=lf archive`
  from `a392c0a` and verified file by file by git blob SHA-1 against `git ls-tree -r`: **750/750
  files, 0 missing, 0 extra, 0 differing.** It was then `git init`-ed so the artifacts record a clean
  worktree; its tree object is **`b92991ad1e110c98b2dcee4e7195226317c07c96`, byte-identical to
  `a392c0a^{tree}`**.
* **The two arms are code-matched in the training path, verified rather than assumed.** The stable
  arm's five repeats ran at commit **`9cf6c84`** (`~/morpheus-rebase-d1`, per
  `rep1.npz.diagnostics.json`). Against this workspace, after CR normalisation, **`v2/runner.py`,
  `v2/model.py`, `v2/data.py`, `v2/training.py`, `v2/losses.py`, `v2/contracts.py`, `v2/preflight.py`,
  `v2/provenance.py`, `v2/plip.py`, `v2/pbs.py`, `v2/slide_pretraining.py`, `v2/export.py` and all of
  `src/training/` are identical.** The one module the runner imports that *did* change is
  `v2/calibra/spectral.py`, and the two symbols it imports — `effective_rank` and `RANK_VARIANTS` —
  are identical by AST (the changes are in the held-out CCA family, which the runner does not call).
  So the arm swap is the only difference in the training path, which is the whole point of the
  comparison.
* **Queue.** At 00:37 UTC the A100 was at 100% utilisation, 61.7 of 81.9 GB, running **ten**
  concurrent `d1_momentum_probe.py` 600-step jobs from another agent's `~/ws_j2`
  (`launch_pf10.sh`, wave 1 of 2, ~3 h per wave). **Nothing was launched into that.** The chain polls
  `nvidia-smi --query-compute-apps` — not a process-name match, so it cannot match itself — and starts
  only after **three consecutive clear checks 120 s apart**, which the seconds-long gap between that
  job's two waves cannot satisfy. Expected start ≈ 05:45 UTC, expected finish ≈ 08:00 UTC (the stable
  arm's five took 1 h 47 m, exports 17 m).
* **How to finish this if the session that armed it does not.** `grep UNSTABLE-ENVELOPE
  ~/e0_run/chain.log` for progress; the chain writes `~/e0_run/d1_envelope_pf/out/`
  `P2_ENVELOPE_FLOORS_PF.json`, `P2_ENVELOPE_FLOORS_PO_RECHECK.json` and
  `d1_envelope_pf_readout.log`. **Read them against §3 of the predeclaration, which fixes the four
  outcomes in advance.** Do not compute a pooled fold across the ten runs of the two arms: the
  convention already fixed in `p2_probe_floors.combine()` is `max` of the two arms' folds with the
  carrying arm named, because a pooled fold would score the arms' genuine difference as noise.

## 4. What the answer will and will not change — the arithmetic, fixed in advance

Every resolvable count below is `# of the six pair ratios > floor`, and the six ratios are already on
disk (`~/e0_run/P2_ROBUSTNESS.json`; the twelve D2/D1 artifacts).

| view | the six pair ratios, sorted | count as a function of the floor `f` |
|---|---|---|
| `wsi_biology` | 1.004, 1.186, 1.573, 1.738, 2.190, 3.246 | 0/6 at 3.295× and **0/6 at anything larger** |
| `rna_biology` | 1.116, 1.205, 1.238, 1.766, 2.852, 3.014 | ≤1.116 → 6/6; ≤1.205 → 5/6; ≤1.238 → 4/6; ≤1.766 → 3/6; ≤2.852 → 2/6; ≤3.014 → 1/6; else 0/6 |
| `full_biology` | 1.042, 1.140, 1.234, 2.248, 3.606, 5.250 | ≤1.042 → 6/6; ≤1.140 → 5/6; ≤1.234 → 4/6; ≤2.248 → 3/6; ≤3.606 → 2/6; ≤5.250 → 1/6; else 0/6 |

For calibration of how much is at stake: **if the unstable arm's same-seed `rna_biology` floor came
back anywhere near its five-seed value of 1.978×, the count would be 2 of 6, not 6 of 6** — and §1.3's
*"the floor is fifty times smaller and every between-arm difference clears it"*, §1.4's contribution 3,
§4.1b's box and §4.5(c) rows 27–28 would all need rewriting rather than renumbering. That is outcome
**(B)** or **(C)** of the predeclaration, and which one it is depends on a number that does not exist
yet.

## 5. A scope point the paper does not currently make, and should

§4.1a rows 26–28 apply this floor to **six** pairs, of which **three are D2** (arm H, Hallmark, against
arm I, PBS). Neither D2 arm is `programme_only` or `programme_free`. **No retraining floor of any kind
has ever been measured on a D2 arm**, so half of every "0 of 6 / 6 of 6" count in §4.1b is judged
against a floor measured on a different experiment's arms. The both-arm measurement in flight closes
this for the three D1 pairs and **cannot** close it for the three D2 pairs. That is a separate GPU cost
(five same-seed retrains per D2 arm) and it belongs on §6.2's list, which does not currently carry it.

## 6. Suite

Baseline before any of this, at `288a124`, thread-capped, `--basetemp=./pytmp`, the repository reached
as `morpheus/` through a junction outside the tree so nothing in the repo was modified to run it:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python -m pytest v2/tests tests -q --basetemp=./pytmp
→ 672 passed, 2 skipped, 445 warnings in 154.67s
```

No test, module or statistic was written for this entry: `p2_envelope_floors.py` was run unchanged,
and every rank statistic in it is imported from `v2/calibra/spectral.py` and
`v2/research/rebase/p2/p2_competing_metrics.py`.

## 7. Files

- Measured (this entry): `~/e0_run/pf_seedaxis/out/SEEDAXIS_f.json`, `SEEDAXIS_p.json`,
  `seedaxis_{f,p}.log`
- In flight: `~/chain_unstable_envelope.sh`, `~/ws_uf/morpheus` (750/750, tree `b92991ad`),
  `~/e0_run/d1_envelope_pf/`
- Read, unchanged: `v2/research/rebase/p2/p2_envelope_floors.py`,
  `v2/research/rebase/d1_envelope_readout.py`, `~/chain_retrain_envelope.sh`
- **Not touched:** `paper/P2_RANK_DRAFT.md`, `v2/research/rebase/p2/floor_audit.json`,
  `v2/calibra/claim_guards.py`, `claim_evidence.json`, any other agent's `PREDECLARED_*`, and
  `v2/research/rebase/p2/figures/extract_from_box.py` (another agent has it open with uncommitted
  work, so the vendoring of this entry's outputs is deliberately deferred rather than merged into
  their edit).
