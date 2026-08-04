# PREDECLARED — is the inductive channel result a property of *one* discovery/exposure split?

**Written and committed 2026-08-05 00:15 UTC, before any code for it was written and before any
number below exists.** Nothing in this file may be edited after the commit; the result entry grades
itself against it verbatim.

---

## 0. The narrowing this run exists to attack

`NOTEBOOK_ENTRIES/inductive_channel_and_ceiling_result_20260804T2345Z.md` reports that
`d2_h::wsi_biology`'s morphology→molecular channel retains **0.9966** of its matched transductive
control's excess when the nuisance model is fitted on a separate discovery fold and applied out of
sample (`d2_i`: **0.9710**). That entry's own §12 names the constraint first: **"One split.
Discovery fraction 0.5 at seed 42 only."**

The published sentence is therefore currently licensed to say *"the channel survives **this**
inductive adjustment"*. The project's standing instruction is that a narrowed claim is not accepted
and that the thing should be pushed until it fails. The specific failure mode being hunted: the
seed-42 partition is one draw of ~1,382 patients out of a 2,766-patient certification cohort, and if
retention is a partition-dependent quantity then 0.9966 is a draw from a distribution nobody has
looked at, quoted as if it were a constant. That is structurally the same defect the ceiling
calibration turned out to have (6.0% was a property of one labels encoding, not of the cohort), and
it is the defect this run is written to find.

---

## 1. What will be run

**The existing machinery, unchanged.** Operator `v2/calibra/inductive_adjustment.ConfoundAdjustmentOperator`;
split `p4_certify.exposure_split`; state `p4_certify.prepare_state`; channel
`nonlinear_adjustment.channel_under_adjustment` / `labels_only_ceiling` / `retention_of_excess` /
`adjuster_agreement` / `cross_fitted_r2`; cohort assembly `nonlinear_adjustment._load_block`
(`run_calibra`'s own path). Driver `v2/research/rebase/nature/p1_evidence/inductive_channel.py`.
**No second operator, no second split function, and no statistic defined inline.**

**The one code change permitted to this run**, and it is a driver change, not a library change:

1. `--split-seed` on `inductive_channel.py`, defaulting to `--seed`, threaded **only** into
   `exposure_split`. This is what isolates the question. `--seed` continues to drive the ridge
   `KFold`, the operator fit, the held-out S2 split and the permutation RNG, so that a difference
   between two runs is a difference of **partition** and of nothing else. When `--split-seed`
   differs from `--seed` the run's P4 `prepare_state` bit-identity comparison is *skipped and
   recorded as skipped* rather than reported as a False — `prepare_state` couples the two seeds and
   there is no P4 state to be identical to on a partition P4 never measured.
2. Site-coverage reporting in the provenance block, read from
   `SitePooling.apply`'s **existing** report: how many exposure rows have a site the discovery fold
   saw at least `min_site_count` times. This is the 43.8% figure of
   `p4_inductive_adjustment_measured_20260804T2300Z.md` §4, and it is the covariate §4 below
   regresses retention against. It is a count, not a statistic.

Both are default-preserving: `--split-seed` unset reproduces the published command exactly, and the
10 tests in `v2/tests/test_inductive_channel.py` import the module's functions, not `main`.

### The arms

**Primary — 8 partitions at the published discovery fraction.** `--discovery-fraction 0.5`,
`--seed 42` fixed, `--split-seed ∈ {42, 7, 11, 23, 101, 555, 2718, 31337}`, on **`d2_h_seed42`** and
**`d2_i_seed42`**, state `wsi_biology`, partition `test`, targets `frozen_rna_targets.npz`, 90
non-control targets, 16 components, 2,000-permutation within-cancer pairing null, `min_site_count`
10. Every split carries its own **matched transductive control at the same n on the same rows**, and
retention is always inductive-vs-that-control, never against a published n = 2,766 number.
`--split-seed 42` is the **reproduction gate**: it must return 0.9966 / 0.9710 to the digit, and the
n = 2,766 `transductive_full` arm is run in that one run only, where it must return 0.6052 / 0.1483 /
0.4569.

**Secondary — a deliberate coverage gradient.** `--split-seed 42`, `--discovery-fraction ∈ {0.3,
0.7}`, `d2_h` only. Eight random partitions at a fixed fraction will vary discovery-fold site
coverage only slightly (the coverage is mostly a property of the site-size distribution, not of the
draw), so a seed sweep alone cannot answer "does worse coverage cost retention". This arm moves
coverage on purpose. P4 measured 43.8% coverage at f = 0.5 and 54.8% at f = 0.7 for the site
certificate; f = 0.3 should fall well below 43.8%.

**Ceilings.** All three labels encodings (`additive_design`, `saturated_cell_design`,
`frozen_discovery_design`) on every split, both exposure arms, because
`inductive_channel_and_ceiling_result_20260804T2345Z.md` puts the verdict's weight on the ceiling and
not on retention alone. The ceiling's own split-stability is a second question this run answers for
free.

---

## 2. Predictions — point values, committed now

| quantity | prediction | p |
|---|---|---|
| `d2_h` retention, **median over the 8 partitions** | 0.975 | — |
| `d2_h` retention, **min over the 8 partitions** | 0.93 | — |
| `d2_h` retention, **max − min over the 8 partitions** | 0.07 | — |
| all 8 `d2_h` partitions land in 0.90–1.05 | yes | 0.60 |
| all 8 `d2_h` partitions land in 0.97–1.00 (the range the single result implies) | **no** | 0.75 |
| `d2_i` retention, median / min over the 8 partitions | 0.960 / 0.90 | — |
| every partition's `permutation_p` at the 0.0005 floor, every arm | yes | 0.90 |
| discovery-fold site coverage at f = 0.5, spread over the 8 partitions | within ±3 percentage points of 43.8% | 0.80 |
| retention at f = 0.3 **below** retention at f = 0.7 by more than 0.03 | **yes** — coverage buys retention | 0.45 |
| the `frozen_discovery_design` ceiling stays ≤ 0 or insignificant in the inductive arm on **all** 8 partitions | yes | 0.55 |

**The prediction the conclusion rests on is the fourth row: that not all eight partitions land in
0.97–1.00.** I expect the single-split number to be a favourable-to-typical draw from a distribution
about 0.03–0.07 wide, i.e. the *claim* survives and the *precision of the quoted figure* does not. If
all eight land inside 0.97–1.00 I was wrong about the width, and the single-split result was
representative in the strongest available sense.

---

## 3. What would make me say the single-split result WAS representative

All of the following, jointly:

* **R1.** The reproduction gate at `--split-seed 42` returns retention 0.9966 (`d2_h`) and 0.9710
  (`d2_i`) and S1 0.6052 / null 0.1483 / excess 0.4569 on the n = 2,766 arm, to the digit.
* **R2.** Every one of the 8 partitions returns `d2_h` retention **≥ 0.90**, and `d2_i` **≥ 0.85**.
* **R3.** `max − min` over the 8 partitions is **≤ 0.10** on `d2_h`. The claim then reads "retention
  is 0.9x with a spread of ≤ 0.10 over partitions", which is a claim about the cohort.
* **R4.** No partition's inductive arm loses the permutation floor while its matched control keeps
  it.
* **R5.** The seed-42 value is not an outlier of the 8: it sits within the observed range and is not
  the maximum by more than 0.02.

## 4. What would make me say it was NOT representative — the narrowing fires

Any one of these, and the result entry must say so in its first paragraph:

* **N1.** Any partition returns `d2_h` retention **< 0.90**, or `d2_i` **< 0.85**. The published
  0.9966 is then one tail of a distribution and the paper may not quote it as *the* retention.
* **N2.** `max − min` over the 8 partitions **> 0.10** on `d2_h`. The quantity is partition-dependent
  and must be quoted as a range with n_splits stated, never as a point.
* **N3.** The seed-42 partition is the **maximum** of the 8 by more than 0.02. That is the signature
  of a lucky draw, and it fires even if every value clears R2.
* **N4.** Any partition's inductive arm falls inside its own permutation null (p > 0.05) while its
  matched transductive control does not.
* **N5.** The `frozen_discovery_design` ceiling turns **significantly positive** (p < 0.05 with excess
  > 20% of the channel's) in the inductive arm on any partition. §6.2 of the source entry called that
  encoding's zero result "the strongest form of P1 §5's argument"; if it is partition-dependent, that
  sentence needs the same treatment as the 6.0%.

**If N1–N3 fire, the mechanism must be reported, not just the number.** Specifically: retention is
regressed against the split's own discovery-fold site coverage and design width — both recorded per
split — and against the f = 0.3 / 0.7 gradient arm. A spread that tracks coverage is a mechanistic
finding about what an inductive adjustment needs in its discovery fold; a spread that does not track
coverage is sampling noise in a 1,382-row measurement and must be reported as such, with no
mechanism invented for it.

## 5. What would make me distrust a FAVOURABLE result (all eight tight)

1. **The split not actually changing.** If two `--split-seed` values produce the same exposure row
   set, `--split-seed` is not wired to `exposure_split`. Guard: the exposure patient-id set of every
   pair of partitions is compared, and the pairwise Jaccard overlaps are reported. Expected ≈ 0.33
   for two independent halves; anything above 0.9 voids the run.
2. **The adjustment degenerating to a constant subtraction on a new partition.** A discovery fold
   whose pooling keeps no frequent sites would leave a cancer-only design, and the "inductive" arm
   would then be measuring something else entirely. Guard, checked **before** the retention number of
   that split is trusted, exactly as the source entry did: `n_frequent_sites_in_discovery_fold ≥ 20`,
   operator design width ≥ 45 columns, `_adjustment_audit` median raw-vs-adjusted correlation in
   0.65–0.85, residual variance ratio in 0.50–0.70, and **0** of 256 axes above 0.99. A split that
   fails any of these is reported with its numbers **and** its degeneracy flag, and is excluded from
   the spread only if it is excluded loudly.
3. **Retention flattered by a collapsing null.** Same check as the source entry §5.1: if a split's
   inductive null median falls more than 30% below its matched control's, its retention is
   uninterpretable and is reported as such. Both the within-cancer and the global pairing null are
   computed for every arm of every split, and both are tabulated.
4. **Fold leakage.** 0 patients in both folds is asserted in the driver for every split; a run that
   raises is a void run, not a repaired one.

## 6. What would make me distrust an UNFAVOURABLE result

1. **A partition that refuses.** `on_unseen_level="refuse"` stays on. A split that puts a cancer
   wholly in the exposure fold stops the run; that is a **defect of the split**, and it is reported
   and the seed replaced by the next one in a predeclared list ({4242, 9, 77} in that order), never
   silently retried until a favourable one appears.
2. **The transductive control moving too.** If a partition's retention is low because its *control*
   moved rather than its inductive arm, that is a statement about the partition's difficulty, not
   about inductive adjustment. Both arms' raw S1, null median and excess are tabulated for every
   split so this is readable, and the inductive arm's excess is additionally reported against the
   *pooled* control excess.
3. **n.** Every comparison is at matched n within a split. No cross-split comparison of a raw S1 is
   made.

## 7. Reporting rules

* Report the spread first — the range, the median, the min, the max and every individual value — and
  the individual numbers second. A table of eight retentions with no dispersion statement would
  repeat the defect this run exists to find.
* Report the awkward direction first whichever way it points, per the project rule.
* Grade §2's ten predictions verbatim, including the ones that were wrong.
* `claim_guards.py`, `claim_evidence.json`, other agents' `PREDECLARED_*` files and
  `paper/P1_CALIBRA_DRAFT.md` are **not** edited. The exact prose location that this result bears on
  — the residual-bound paragraph in **§4.2** of `paper/P1_CALIBRA_DRAFT.md`, which was updated once
  today with the single-split figure, and **§4.4**'s retention sentence — is named in the result
  entry for the main session to update, and left alone here.
* The full test suite is run before the final commit and its pass/fail counts quoted verbatim.

## 8. Cost and honesty about it

Per split, skipping the n = 2,766 gate: ~90 s of channel arms plus ~50 s of ceilings on 8 CPU
workers. 18 runs ≈ 45 min. CPU only; the GPU is not required and will not be touched. If the box is
contended the `d2_i` partitions are the arm that gets dropped, and the entry will say which of the
predeclared runs did not happen rather than quietly reporting fewer.
