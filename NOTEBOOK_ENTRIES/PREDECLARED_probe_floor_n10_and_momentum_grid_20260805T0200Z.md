## 2026-08-05 02:00 UTC — PREDECLARED: pushing §5.4 limit 2 (m = 0.999 over m = 0.99, step 600) until it breaks or holds

**Logged:** 2026-08-05 02:00 UTC, **before any run was launched and before any statistic beyond
R1/R3 was read.** **Machine:** A100 80 GB, `ubuntu@150.136.45.194`. GPU checked idle at 02:00 UTC
(`nvidia-smi`: 0 MiB / 81920 MiB used, 0 % utilisation; no `python` training process owned by any
agent). Nothing is queued behind another agent's job.

### 0. What is being pushed on, and why it is worth pushing on

`NOTEBOOK_ENTRIES/three_floors_close_the_last_three_unjudgeable_rows_20260805T0000Z.md` §4 closed
P2 §4.1a's last unjudgeable row, §5.4 limit 2 — **the comparison at the exact momentum this project
ships**. It closed it in a state that is recorded as fragile in the draft's §"Status" item 12:

| test | statistic | value | floor | verdict at n = 5 |
|---|---|---:|---:|---|
| **A** — the audit rule: published single-draw ratio > floor | R3 | 1.262× | 1.195× | **clears by 5.6 %** |
| **B** — worst case over repeats: min(m0.999) / max(m0.99) > floor | R3 | **1.138×** | 1.195× | **FAILS** |
| **B** | canonical R1 | **1.453×** | 1.155× | **clears** |

Test A is what `p2_floor_audit.check` computes and is the paper's rule. Test B is the stronger
statement the same ten runs support, and R1 and R3 disagree about it. This predeclaration fixes, in
advance, what would make that disagreement resolve.

### 1. The arithmetic constraint that makes half of this uninformative, stated first

**Both halves of Test B are monotone in the number of repeats and both move against the pass.** The
floor is `max/min` over an arm's repeats, so it is non-decreasing in n. The worst-case separation is
`min(high arm)/max(low arm)`, so it is non-increasing in n. **Adding repeats therefore cannot make
Test B easier for any statistic, ever.**

Consequences, predeclared so they cannot be spun afterwards:

1. **"R3's Test B got worse at n = 10" is NOT a finding and will not be reported as one.** It is
   arithmetically guaranteed not to improve.
2. **"R1 still passes Test B at n = 10" IS a finding**, because it is a test that had every
   opportunity to break and did not.
3. **Test A can flip in either direction** — the ratio 1.262× is fixed (it is a property of two
   specific published runs) while the floor grows — so a floor that grows past 1.262× is a **hard
   break of the row by the paper's own checker**, and that is the single sharpest falsifier available
   here.

### 2. What will be run

**Item 1 — more repeats at the shipped setting. 10 runs.** Five further same-seed repeats of each of
the two arms of §5.4 limit 2, taking each arm from n = 5 to **n = 10**:

```
d1_momentum_probe.py 0.999 0.04 600 4096 2e-4 42 ~/e0_run/d1_probefloor600/m0.999_rep{6..10}
d1_momentum_probe.py 0.99  0.04 600 4096 2e-4 42 ~/e0_run/d1_probefloor600/m0.99_rep{6..10}
```

Identical to repeats 1–5 in every argument, in the workspace (`~/ws_j2`), in the seed (42), in the
step budget (600), in the objective profile (`programme_free`) and in the launch concurrency (10
simultaneous runs on one A100 — concurrency is part of the conditions the original floor was
measured under and is held fixed deliberately). The only difference between any two repeats is GPU
non-determinism, which is the whole quantity being measured.

**Item 2 — a third rank statistic. 0 runs.** Runs first, because it needs no GPU. `RANK_VARIANTS`
carries six variants (`R1`=`CANONICAL`, `R2`, `R3`, `R1_uncentred`, `R1_rownorm`, `R2_uncentred`) and
`p2_envelope_floors.STATISTICS` carries five more published alternatives. **No new statistic will be
invented and none will be computed inline**; every number comes from
`effective_rank(..., variant=RANK_VARIANTS[...])` or from `p2_competing_metrics`, through
`p2_probe_floors.score_state`, off the already-saved `probe_step600.npz` states.

**Item 3 — adjacent momentum values. 10 runs.** Five same-seed repeats at each of **m = 0.995** and
**m = 0.98**, everything else identical, so the shipped comparison sits inside a four-point grid
{0, 0.98, 0.99, 0.995, 0.999} rather than being a two-point contrast:

```
d1_momentum_probe.py 0.995 0.04 600 4096 2e-4 42 ~/e0_run/d1_probefloor600/m0.995_rep{1..5}
d1_momentum_probe.py 0.98  0.04 600 4096 2e-4 42 ~/e0_run/d1_probefloor600/m0.98_rep{1..5}
```

**Total: 20 GPU runs**, two waves of 10 at ~3 h per wave (the original 10-concurrent step-600 wave ran
17:43 → 20:44 UTC on 2026-08-04). Item 1's wave goes first; if the box is claimed by another agent
before wave 2, item 3 is abandoned and reported as not run rather than run at a different
concurrency.

### 3. Workspace provenance, verified before launch rather than asserted

`~/ws_j2` — the workspace the original fifteen step-600 repeats were run from — was checked
file-by-file against a manifest generated from the local canonical checkout at `HEAD` (`6f28814`).
It differs from `HEAD` in 15 files and is missing 64, **none of which is on the training path or the
scoring path**. Every file that is was confirmed byte-identical after LF normalisation:

```
SAME  v2/research/rebase/d1_momentum_probe.py   (verified byte-identical in the 2026-08-05 00:00 entry)
SAME  v2/training.py   v2/runner.py   v2/model.py   v2/preflight.py
SAME  src/training/train_bio_query_former.py
SAME  v2/research/rebase/p2/{p2_probe_floors,p2_envelope_floors,p2_competing_metrics}.py
```

`v2/calibra/spectral.py` **is** in the differing list, and it is the one that matters, so it was
diffed rather than waved through: the entire block from `class RankVariant` to `def cca_spectrum` —
`RankVariant`, `CANONICAL`, `RANK_VARIANTS`, `EPS`, `DEGENERACY_TOLERANCE`, `_finalise`,
`effective_rank` — is byte-identical (10,050 characters, both sides). The only changes are in the
held-out-CCA family, which no number in this entry touches. **Repeats 6–10 are therefore the same
experiment as repeats 1–5, not a re-run of a moved target.**

### 4. What resolves the fragility, decided in advance

Every statistic in `p2_envelope_floors.STATISTICS` plus every variant in `RANK_VARIANTS` is scored on
`wsi_biology` at step 600, and each gets Test A and Test B at n = 10.

**BREAKS — the fragility is real and worse than recorded.** Any one of:

- **B1.** The n = 10 **R3 floor ≥ 1.262×**. The row then stops clearing under Test A, the paper's own
  rule, and §5.4 limit 2 goes back to failing rather than passing. This is the primary falsifier.
- **B2.** **R1 fails Test B at n = 10** (`min(m0.999) / max(m0.99) ≤ R1 floor`). The R1-passes side of
  the disagreement was then a five-repeat artifact, both statistics fail the strong test, and the
  pass rests on nothing but the single draw §5.2 happens to record.
- **B3.** A majority of the statistics scored (≥ 6 of the 11 in `STATISTICS`) fail Test A at n = 10.

**HOLDS — the fragility resolves in favour of the pass.** All of:

- **H1.** R1 still passes Test B at n = 10, having had ten repeats per arm to break it.
- **H2.** R3 still passes Test A at n = 10 (floor < 1.262×).
- **H3.** At least one statistic other than R1 and its algebraic duplicates passes Test B at n = 10 —
  i.e. the pass is not carried by one statistic alone.

**AMBIGUOUS — and this is a reportable outcome, not a failure to conclude.** R1 passes Test B, R3
fails it, Test A survives, and no third *independent* statistic breaks the tie. If this is what
happens it will be reported as "we doubled the repeats and it is still a coin flip between two
statistics", in those terms.

### 5. What would make me distrust a favourable result

1. **A batch effect between repeats 1–5 and 6–10.** They were launched a day apart. If, for either
   arm, the range of repeats 6–10 does not overlap the range of repeats 1–5 under R1 or R3, then the
   widened floor is a between-batch difference and not the within-batch GPU non-determinism the floor
   claims to be — and I will say so and will **not** treat the n = 10 floor as the same quantity as
   the n = 5 one. Reported either way.
2. **A run that did not train.** Any repeat whose step-600 `biology_contrastive` is at chance
   (ln 80 = 4.382) or whose R1 is at the collapsed arm's level (< 4) is excluded and named, with the
   n it leaves behind stated. It is not silently replaced.
3. **Concurrency drift.** If wave 1 does not run 10-at-a-time on an otherwise idle card, that is
   recorded beside the floor.
4. **A statistic that is not independent counted as if it were.** §6 of the 00:00 entry established
   that R2 and R3 coincide on this block because `z_biology` is L2-normalised at output, so R3's row
   normalisation is a no-op. Any variant that is algebraically the same statistic on this block will
   be **checked numerically and reported as a duplicate**, not counted as corroboration. I expect this
   to catch `R1_rownorm` (= R1) and `R2` (= R3) and I will verify rather than assume.
5. **A favourable result at n = 10 that a sixth-to-tenth-repeat-only floor contradicts.** As a
   secondary read, the floor computed from repeats 6–10 alone (an independent n = 5) is reported
   beside the n = 5 and n = 10 floors. If the independent five disagree with the original five about
   whether Test B passes, that is itself the answer and it is "unstable", not "passes".

### 6. What item 3 would show, decided in advance

With five same-seed repeats at each of m ∈ {0, 0.98, 0.99, 0.995, 0.999} at step 600:

- **STEP** — R1 (and R3) are flat within the floor across 0.98 / 0.99 / 0.995 and rise only at 0.999.
  The two-point comparison is then well-posed: `m = 0.999` is a genuine discontinuity and the closeness
  of 0.99 is not evidence against it.
- **SMOOTH** — R1 rises monotonically across the grid with each adjacent gap of the same order as the
  floor. The two-point comparison is then a coarse read of a smooth trend, the "fragility" of
  0.999-vs-0.99 is just adjacent-point closeness, and the honest claim is about the *trend*, not about
  `0.999` beating `0.99`.
- **NEITHER** — non-monotone in m. Reported as such; it would undercut the momentum story more
  broadly than §5.4 limit 2 does.

### 7. Rules this entry binds itself to

- No threshold in this file may be moved after a number is seen. The verdicts are computed by the
  same `p2_probe_floors.py` / `p2_floor_audit.check` path, at the arms and step written into the block
  string, and nothing is borrowed across arms, steps, learning rates or capacities.
- **`paper/P2_RANK_DRAFT.md` is not edited by this work.** A completeness-audit agent is editing it
  concurrently. Findings go to `NOTEBOOK_ENTRIES/` with the exact prose location flagged.
- `claim_guards.py`, `claim_evidence.json` and other agents' `PREDECLARED_*` files are not touched.
- The negative or awkward finding is reported first.
