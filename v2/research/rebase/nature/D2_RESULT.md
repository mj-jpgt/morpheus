# D2 — Hallmark vs Perturbation-Basis Supervision (PBS), 3 seeds

**Logged:** 2026-08-02, recovered from the live session transcript after the Lambda instance was
stopped. **The underlying artifacts were on ephemeral storage (`~/e0_run/d2_v2/`) and are presumed
lost**; the numbers below were read out of the completed run before shutdown and are recorded here
because re-deriving them costs ~6 GPU hours.

**How obtained:** `phase_d d2 --seeds 42,43,44 --restrict-to-split`, 40 epochs/arm, 6,427-patient
maximal paired split, 11 development cancers / 21 held-out. Six runs, all reaching `TRAIN_SUCCESS`.
Readout `morpheus.v2.research.rebase.d2_compare`, 3 seed-matched artifact pairs, cancer+TSS
residualisation (84 sites kept, min 10/site), top canonical correlation at 16 components,
paired patient and cancer bootstrap at 2,000 repeats. `n_test = 2,766`.

## Primary readout (as run: all 90 non-control targets)

| seed | Hallmark | PBS | Δ (PBS−H) | patient CI₉₅ | p | cancer CI₉₅ | p |
|---|---:|---:|---:|:---:|---:|:---:|---:|
| 42 | 0.5861 | 0.4762 | −0.1100 | [−0.1374, −0.0673] | 0.0000 | [−0.1444, −0.0158] | 0.0065 |
| 43 | 0.5879 | 0.4983 | −0.0896 | [−0.1230, −0.0527] | 0.0000 | [−0.1243, +0.0016] | 0.0280 |
| 44 | 0.6051 | 0.4935 | −0.1117 | [−0.1366, −0.0758] | 0.0000 | [−0.1426, −0.0144] | 0.0095 |

PBS loses to Hallmark by ~0.10 held-out molecular CCA in all three seeds. The patient-level CI
excludes zero in 3/3; the cancer-level CI excludes zero in 2/3 (seed 43 grazes zero at +0.0016).

## ⚠ This number is NOT yet interpretable

`frozen_rna_targets.npz` carries 180 targets in five groups. `d2_compare._targets` drops the 90
`RANDOM_CONTROL__` targets, leaving 90 — of which **50 are `hallmark_in_training`, i.e. the Hallmark
arm's own supervision targets.** 56% of the readout is one arm's training signal, which hands that arm
a structural advantage.

`--target-groups` was added to `d2_compare.py` for exactly this, and the fair contrast —
`heldout_pathway` (24) + `immune_tme` (8) + `tumour_state` (8) = **40 targets neither arm trained on** —
**was launched but did not complete before the instance was stopped.** It must be re-run before any
version of this result is quoted. It is CPU-only and takes minutes; it does not need the checkpoints
to be retrained *if* the six `.npz` artifacts survived to persistent storage. If they did not, this
costs a full D2 re-run.

**Verdict pending.** If the gap survives on the 40 untrained targets, PBS genuinely underperforms
Hallmark supervision and P3's headline hypothesis is refuted by its own predeclared test. If it
collapses, the −0.10 was largely measuring which arm's training targets were on the exam.

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
