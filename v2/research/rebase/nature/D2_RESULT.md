# D2 — Hallmark vs Perturbation-Basis Supervision (PBS), 3 seeds

**Logged:** 2026-08-03. **Status: COMPLETE AND QUOTABLE.** The stratified readout that was missing
has been run, and it decides the claim in the direction the unrestricted number already pointed.

**Provenance note (retained).** An earlier D2 sweep completed on 2026-08-01 but its outputs were
written to ephemeral storage and were destroyed when the Lambda instance stopped. The numbers from
that run were recorded from the live session transcript and are preserved below for comparison. The
result reported here is a **complete re-run** (`d2_v3`), not a recovery of those numbers: all six
arms were retrained and all outputs are on persistent NFS under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/e0_run/d2_v3/`.

**How obtained:** Lambda A100 80 GB (`150.136.45.194`).
`morpheus.v2.research.rebase.phase_d d2 --seeds {42|43|44} --restrict-to-split
--pbs-components 128 --analysis-role primary --token-budget 8192`, three concurrent run roots,
40 epochs/arm, 6,427-patient maximal paired split (3,118 train / 543 val / 2,766 test),
11 development cancers / 21 held out. Six arms, **6/6 `TRAIN_SUCCESS`**, all six representation
artifacts exported by `morpheus.v2.export`. CALIBRA G4 controls run per seed:
`gates_pass: true`, `rna_positive_control_passed: true`, `channel_gate_failures: []` for all three.
Readout `morpheus.v2.research.rebase.d2_compare`, cancer+TSS residualisation (84 sites kept,
min 10/site), top canonical correlation at 16 components, paired patient and cancer-cluster
bootstrap at 2,000 repeats. `n_test = 2,766`. The `p` columns are `p_improve`, the bootstrap
probability that PBS *beats* Hallmark.

## 1. Primary readout, unrestricted (all 90 non-control targets)

| seed | Hallmark | PBS | Δ (PBS−H) | patient CI₉₅ | p | cancer CI₉₅ | p |
|---|---:|---:|---:|:---:|---:|:---:|---:|
| 42 | 0.6214 | 0.4855 | **−0.1359** | [−0.1547, −0.0943] | 0.0000 | [−0.1589, −0.0517] | 0.0010 |
| 43 | 0.6080 | 0.5003 | **−0.1077** | [−0.1366, −0.0699] | 0.0000 | [−0.1373, −0.0014] | 0.0240 |
| 44 | 0.6075 | 0.4883 | **−0.1192** | [−0.1472, −0.0836] | 0.0000 | [−0.1537, −0.0303] | 0.0030 |

Reproduces the lost run in direction, magnitude and seed ordering (previously −0.1100 / −0.0896 /
−0.1117). Both CIs now exclude zero in **3/3** seeds, where the earlier run had the cancer CI graze
zero on seed 43.

## 2. THE STRATIFIED READOUT — 40 targets neither arm trained on

`--target-groups heldout_pathway immune_tme tumour_state` (24 + 8 + 8 = 40 targets).

| seed | Hallmark | PBS | Δ (PBS−H) | patient CI₉₅ | p | cancer CI₉₅ | p |
|---|---:|---:|---:|:---:|---:|:---:|---:|
| 42 | 0.6126 | 0.4800 | **−0.1325** | [−0.1605, −0.0993] | 0.0000 | [−0.1792, −0.0632] | 0.0010 |
| 43 | 0.5970 | 0.4882 | **−0.1089** | [−0.1460, −0.0749] | 0.0000 | [−0.1623, −0.0118] | 0.0105 |
| 44 | 0.5983 | 0.4757 | **−0.1226** | [−0.1502, −0.0866] | 0.0000 | [−0.1653, −0.0411] | 0.0000 |

**The gap survives, at full size, on targets neither arm was trained on.** Stratified minus
unrestricted is +0.0034 / −0.0012 / −0.0034 — i.e. **no material change in any seed**. Patient and
cancer CIs both exclude zero in 3/3.

### The contamination worry was real but ran the other way

The concern was that 50 of the 90 unrestricted targets are `hallmark_in_training`, the Hallmark
arm's own supervision, inflating its score. Scoring that group alone:

| seed | Hallmark | PBS | Δ on H's own targets | Δ on the 40 untrained |
|---|---:|---:|---:|---:|
| 42 | 0.6203 | 0.5112 | −0.1091 | **−0.1325** |
| 43 | 0.6009 | 0.5222 | −0.0787 | **−0.1089** |
| 44 | 0.6076 | 0.4927 | −0.1149 | **−0.1226** |

The gap is **smaller** on Hallmark's own supervision targets than on the untrained ones, in all three
seeds. Restricting to H's training targets *understates* PBS's deficit rather than manufacturing it.

## 3. Negative control — the 90 `random_control` targets

| seed | Hallmark | PBS | Δ (PBS−H) | patient CI₉₅ | p | cancer CI₉₅ | p |
|---|---:|---:|---:|:---:|---:|:---:|---:|
| 42 | 0.4671 | 0.4572 | −0.0099 | [−0.0591, +0.0123] | 0.1600 | [−0.1055, +0.0099] | 0.0645 |
| 43 | 0.4681 | 0.4400 | −0.0280 | [−0.0719, −0.0048] | 0.0110 | [−0.0905, +0.0232] | 0.0870 |
| 44 | 0.4637 | 0.4369 | −0.0268 | [−0.0697, +0.0003] | 0.0285 | [−0.0969, +0.0285] | 0.1960 |

**The instrument is not manufacturing arm differences.** The arm gap on random controls is
−0.0099 to −0.0280, i.e. **4–13× smaller** than the −0.109 to −0.133 on real targets. The cancer CI
includes zero in 3/3 and the patient CI includes zero in 2/3.

Two honest caveats, both recorded rather than buried:

- The residual −0.01 to −0.03 is small but not perfectly zero, and it points the same way. Some part
  of the headline gap — of order 10–20% of it — may be a generic representation-quality difference
  rather than supervision content. It does not come close to accounting for the effect.
- The **absolute** level on random controls is high (~0.44–0.47), not ~0. That is not a defect. A
  permutation null (shuffling patient rows of the residualised target matrix, preserving target count
  and covariance, 200 draws) puts chance at **0.140** for every group. Random *gene sets* are not
  random *numbers*: they still track global expression covariance, so a real channel onto them is
  expected. What matters for D2 is the paired arm difference, and that is what collapses.

## 4. Counter-claims

**Effective rank does not explain the gap — it contradicts it.** Roy–Vetterli effective rank of the
residualised held-out `wsi_biology` block (256 nominal dimensions):

| seed | Hallmark | PBS | higher rank | Δ (PBS−H) on untrained 40 |
|---|---:|---:|:---:|---:|
| 42 | 23.39 | 14.87 | H | −0.1325 |
| 43 | 28.77 | **34.12** | **PBS** | −0.1089 |
| 44 | 9.14 | 9.11 | ~equal | −0.1226 |

In seed 43 the PBS arm has **higher** effective rank than the Hallmark arm and still loses by
−0.1089; in seed 44 the two are equal to two decimals and PBS still loses by −0.1226. A capacity or
collapse explanation would have to produce a rank ordering that tracks the score ordering, and it
does not. Flagged separately: effective rank is **unstable across seeds** (9.1 to 34.1 for the same
configuration), so it is a poor summary of these representations, but the arm gap is stable across
exactly the seeds where rank is not.

**Permutation nulls.** Every arm × group is far above its null (0.140) at `permutation_p = 0.005`,
the floor for 200 draws. Both arms are reading a real molecular channel; they differ in how much.

**Seed-42 re-export reproduces the lost point estimate exactly — and retraining does not.**
This is the one number that contradicts a naive expectation and it needs stating plainly:

| | unrestricted top-CCA | effective rank |
|---|---:|---:|
| recorded in the lost run, H seed 42 | 0.5861 | — |
| **re-export of the surviving H seed-42 checkpoint** | **0.58612** | 8.68 |
| **retrained H seed 42 (`d2_v3`)** | **0.6214** | 23.39 |

Re-exporting the surviving checkpoint reproduces the lost number **to five significant figures**, so
the export/readout path is deterministic and the recorded numbers were not mis-transcribed. But
**retraining with the same seed and the same configuration does not reproduce the same model**:
0.6214 vs 0.5861, and effective rank 23.39 vs 8.68. Training is not seed-reproducible on this stack
(GPU non-determinism). The D2 contrast is unaffected because it is paired within a run, but **no
individual D2 point estimate should be quoted as reproducible from the seed alone.**

## 5. What could not be done, and why

**The seed-42 PBS checkpoint could not be re-exported.** The handoff recorded both surviving seed-42
checkpoints as complete with `manifest.json` and `TRAIN_SUCCESS.json`. `d2_i_seed42/` in fact had
neither, only `last.pt` and a 14-line `train_metrics.jsonl`, and the checkpoint's own internal
`epoch` field reads **13**, not 39 — it is a mid-training snapshot copied to persistent storage while
the run was still going. Exported for the record as
`recovered_artifacts/d2_i_seed42_EPOCH13.npz`, it scores 0.5170 unrestricted, nowhere near the
recorded 0.4762, exactly as an undertrained model should. **The plan to re-export both seed-42 arms
was therefore only half possible**, and seed 42 was retrained along with 43 and 44.

**`--target-groups random_control` could not select anything.** `d2_compare._targets` applied the
`RANDOM_CONTROL__` name-prefix drop unconditionally and *then* intersected it with the requested
groups, so the negative control selected the empty set and raised. Fixed to select by group when
`--target-groups` is given; provably identical for every non-control group, since prefix and group
label are in exact 90/90 correspondence in the frozen artifact. **No previous D2 result had ever had
a negative control applied to it.**

**`SUCCESS.json` is absent from the three run roots.** `phase_d`'s final convenience bootstrap ran
with unbounded BLAS threading, which on this box is ~23× slower than single-threaded (SVD of
2766×256: 0.205 s at 1 thread, 4.80 s at 4), and would not have completed for >16 h. Those three
subprocesses were killed and the readout recomputed correctly with `OMP_NUM_THREADS=1` and
parallelism across processes (13 jobs, ~13 minutes total). Everything upstream was verified intact
after the kill: 6/6 `TRAIN_SUCCESS`, all six artifacts, all three CALIBRA gate files.

**Wall clock.** Launch 01:13 UTC → 6/6 trained 07:55 → exports 08:04 → all readouts 12:06.
**6 h 42 min of training**, three seeds concurrent. Concurrency bought less than hoped: the box was
shared with another agent running up to 38 joblib workers, and per-epoch time tracked their load
between 2.2 and 5.6 min/epoch. GPU was never the constraint (~20 GB of 80 GB).

## 6. Verdict

**PBS underperforms Hallmark supervision on the held-out molecular channel, by ~0.11–0.13 top-CCA,
in 3/3 seeds, on 40 targets neither arm was trained on, with both patient and cancer-cluster CIs
excluding zero in 3/3, against a negative control that is 4–13× smaller and mostly indistinguishable
from zero.** P3's headline hypothesis is refuted by its own predeclared test. The refutation is not
an artifact of scoring the Hallmark arm on its own supervision — that restriction makes PBS look
*better*, not worse — and it is not explained by representation rank, which orders the wrong way in
seed 43.

---

# G2.6 — the `programme_free` queue defect, measured on real data

**Logged:** 2026-08-02. **How obtained:** `_overfit_programme_free_contrastive` invoked directly on
the real cohort (3,118 train patients, H-Optimus patch store), hidden 512 / 4 layers / 8 heads,
programme head 256, seed 42, 800 steps, lr 1e-3, 16-patient fixed batch, queue capacity 64. The two
arms differ **only** in the new `freeze_memory` flag.

| | contrastive start → end | reduction | unique queue keys |
|---|---|---:|---:|
| live queue (`freeze_memory=False`) | 4.5755 → 4.3306 | 0.0535 | **16** |
| frozen queue (`freeze_memory=True`) | 4.5755 → **2.7726** | **0.3940** | **64** |

`full_consistency` reached 1.4e-04 (live) and 5.7e-04 (frozen); both arms fit that term fine.

### Technical
The `unique queue keys` column is direct confirmation of the mechanism: replaying one 16-patient batch
against a live queue overwrites all 64 slots with re-encoded copies of those same 16 patients within
4 steps, so the InfoNCE negatives become the queries and the term cannot descend. Freezing the queue
after priming leaves 64 distinct real train patients as keys and yields **7.4× more descent**.

**But the gate still fails.** ln(16) = 2.772589; the frozen arm ends at 2.772559 — in-batch chance to
five decimals. The model defeated the static queue negatives entirely and then could not separate the
16 in-batch patients at all. G2.6 requires `biology_contrastive ≤ 0.10`.

### In plain terms
The model was being asked to tell patients apart while the thing it was being compared against kept
turning into a copy of the patients themselves. Fixing that helped a lot. What's left is a second,
separate problem: the image-side representations of different patients are nearly identical, so
there is nothing to tell apart.

### Meaning for the claim
The queue diagnosis was correct and the fix is necessary — but **not sufficient**, and D1 remains
blocked. The remaining blocker is representational collapse on the WSI biology head (measured
elsewhere at 0.736 mutual collinearity at initialisation), not the memory bank. D1 must not be
launched until the in-batch term can descend below chance.

---

## Paste-ready notebook block

**Logged:** 2026-08-03 12:10 UTC.

**How obtained:** Lambda A100 `150.136.45.194`. Full D2 re-run `d2_v3`: three concurrent
`phase_d d2` streams (`--seeds 42|43|44 --restrict-to-split --pbs-components 128 --token-budget 8192`),
40 epochs/arm, 6/6 `TRAIN_SUCCESS`, CALIBRA G4 `gates_pass: true` on all three seeds. Readouts via
`d2_compare` at 2,000 repeats, per seed-pair, `OMP_NUM_THREADS=1`. Outputs under
`~/e0_run/d2_v3/bootstrap/` and `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json` (persistent NFS).

### Technical
On the 40 targets neither arm trained on (`heldout_pathway` + `immune_tme` + `tumour_state`), PBS
loses to Hallmark by **−0.1325 / −0.1089 / −0.1226** (seeds 42/43/44), patient CI₉₅
[−0.1605,−0.0993] / [−0.1460,−0.0749] / [−0.1502,−0.0866] and cancer CI₉₅ [−0.1792,−0.0632] /
[−0.1623,−0.0118] / [−0.1653,−0.0411] — both exclude zero in 3/3. The unrestricted 90-target readout
gives −0.1359 / −0.1077 / −0.1192, so stratifying changes the answer by at most 0.0034. On the
`hallmark_in_training` group alone the gap is *smaller* (−0.1091 / −0.0787 / −0.1149). The
`random_control` negative control gives −0.0099 / −0.0280 / −0.0268 with the cancer CI covering zero
3/3. Permutation null for every group is 0.140. Effective rank of residualised `wsi_biology` is
H 23.39/28.77/9.14 vs PBS 14.87/34.12/9.11 — PBS is *higher* in seed 43 and still loses.
Re-exporting the surviving seed-42 Hallmark checkpoint reproduces the lost 0.5861 as **0.58612**;
retraining the same seed gives **0.6214**, so training is not seed-reproducible even though the
readout path is deterministic. The surviving seed-42 PBS checkpoint was a mid-training epoch-13
snapshot and could not be re-exported. `--target-groups random_control` had never been runnable
before this run because of a mask bug in `_targets`.

### In plain terms
The worry was that the earlier result only looked bad for the perturbation-based method because more
than half the exam was made of the other method's own homework. We rebuilt the experiment from
scratch and marked it again on forty questions neither method had ever seen. The gap is the same
size, in all three repeats. We also marked both methods on deliberately meaningless questions, and
there the two are nearly tied — so the marking scheme is not simply preferring one method. And when
we checked whether the winning method just had more room to express itself, in one of the three
repeats the losing method had *more* room and still lost.

### Meaning for the claim
D2 is settled and quotable: perturbation-basis supervision is genuinely worse than Hallmark
supervision for held-out molecular structure, and P3's headline hypothesis is refuted by its own
predeclared test. Two things must travel with the number: the negative control is small but not
exactly zero, so ~10–20% of the gap may be generic representation quality; and individual point
estimates are not reproducible from the seed alone, so quote the paired difference, never the level.
