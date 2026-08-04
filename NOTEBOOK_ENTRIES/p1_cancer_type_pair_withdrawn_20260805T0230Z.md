## 2026-08-05 02:30 UTC — the cohort behind `0.463 → 0.035` is found and hash-pinned; the numbers are not reproduced, and are withdrawn

**Logged:** 2026-08-05 02:30 UTC. **Predeclaration:**
`NOTEBOOK_ENTRIES/PREDECLARED_p1_cancer_type_balanced_accuracy_20260805T0130Z.md`, committed at
`6a2392f` before anything was computed. Outcome **C** of the three it names.
**Question:** P1 publishes cancer-type balanced accuracy `0.463 → 0.035`, chance `0.048`, n = 2,530 as
its evidence that confound adjustment works. A provenance audit
(`p1_42_artifact_identified_by_hash_20260804T2210Z.md` §4) concluded it could not be traced to any
artifact. Regenerate it, or withdraw it.
**How obtained:** box `150.136.45.194`, workspace `~/ws_p1prov` shipped as
`git -c core.autocrlf=false archive HEAD` — every tracked file, not a diff. `~/venv`, threads capped
to 1. Every statistic imported from `v2/calibra/`; nothing computed inline.

---

### The awkward part first

**The published pair is not a pair any single estimator produces on this cohort.** Across all 56
readings taken — 7 representation states × 8 estimator variants — not one has both endpoints; the
closest is 0.16 away in L1. Its "before", 0.463, lies in the *nonlinear* band; its "after", 0.035,
lies in the *linear* band. Under a consistent estimator the drop is either far larger than the
advertised 13.2× or far smaller:

| estimator (`wsi_biology`, n = 2,530, 21 classes, chance 0.047619) | raw | adjusted | drop | adjusted ÷ chance |
|---|---:|---:|---:|---:|
| **joint LDA, standardised** — the canonical certificate | **0.733915** | **0.030813** | **23.8×** | **0.65×** |
| joint LDA, unstandardised | 0.728295 | 0.033269 | 21.9× | 0.70× |
| k-NN, k = 15 | 0.444680 | 0.176587 | 2.5× | 3.71× |
| k-NN, k = 15, prior-corrected | 0.508317 | 0.228350 | 2.2× | 4.80× |
| random forest, 300 trees | 0.524007 | 0.268088 | 2.0× | 5.63× |
| RBF-SVM | 0.615676 | 0.295786 | 2.1× | 6.21× |
| per-axis nearest-class-mean, max over 256 axes | 0.164967 | 0.051644 | 3.2× | 1.08× |
| **withdrawn — published 2026-07-30** | *0.463* | *0.035* | *13.2×* | *0.74×* |

So the correction cuts both ways and neither way is flattering to the published figure. On the
**mean-based** reading the adjustment is *better* than advertised — it removes far more cancer
information than 0.463 → 0.035 implies, because the raw representation carries far more than 0.463 of
it (0.734). On the **nonlinear** reading it is *much worse* — the adjusted state still reads 3.7–6.2×
chance, where the published pair implies a residue below chance. The published number's error is not a
decimal: it makes a first-moment certificate look like a modest, honest drop when in fact the
first-moment drop is steep and the thing the certificate does not cover is large.

The last column is a **raw ratio to chance**, not netted against a null that regenerates the
adjustment inside each permutation. It is therefore *not* the same quantity as the 3.45× recorded for
the n = 2,766 arm in `tcga_nonlinear_confound_probe_result_20260804T2100Z.md`, and the two must not be
quoted interchangeably. The direction agrees, which is the point worth carrying.

---

### 1. The cohort was found by structure, not by matching the number

This is the distinction the previous agent's §4 refused to collapse, and it is worth keeping sharp:
identifying a **cohort** from a configuration is not the same move as attributing a **number** to
whichever file reproduces it.

`v2/calibra/run_calibra.py:main` selects its cohort as `split == args.partition` in the input
artifact, intersected with the patients present in the target table. Applying that rule to
`runs_misc/calibra_run/artifacts/diagnostic_full_seed42.npz`:

| quantity | measured | `PHASE1_RESULT.md` header |
|---|---:|---:|
| patients, `split == "test"`, aligned to targets | **2,530** | 2,530 |
| distinct cancer types | **21** | 21 |
| `1 / n_classes` | **0.047619** | chance 0.048 |
| confound columns, `cancer` + pooled `tss`, `min_site_count=10` | **99** | 99 |
| TSS sites kept | **75** | 75 |
| dev cancers (train ∪ val), disjoint from test | **11** | "14" ← **wrong, corrected** |

Six of six on the numbers the header asserts, and a seventh discrepancy found on the way: the header
said the 21 test cancers were disjoint from **14** dev cancers. There are **11**, which is the
project's documented 11-train / 21-test partition. Corrected in place, with the old value recorded.

The split totals are train 3,124 / val 538 / test 2,530 = **6,192** — the **pre-rebuild** cohort. That
is the entire explanation of the n = 2,530 vs n = 2,766 gap that made the pair look uninheritable from
the site arm: the site arm runs on the rebuilt 6,427-patient maximal split (entry 2026-08-01 09:05),
this does not. Nothing exotic; two different splits, five days apart.

Corroborating, independent of any of the above: the committed run `runs/calibra_v2_local/`
(commit `4c7166b`, the commit that created `PHASE1_RESULT.md`) records `n_patients 2530`,
`n_confound_columns 99`, `n_distinct_sites_kept 75` in `task_rows.csv`, and its `heldout_top_cca`
values reproduce the published table digit for digit (0.476787 → 0.477, 0.539273 → 0.539,
0.898314 → 0.898). And `runs_misc/calibra_run/run.log` on the box prints
`[diagnostic_full_seed42::wsi_biology] adj_cca=0.5200 …` against the file's `0.520`, for all seven
states. The artifact is not in doubt.

**Artifacts, by content hash** (under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/`):

| path | SHA-256 |
|---|---|
| `runs_misc/calibra_run/artifacts/diagnostic_full_seed42.npz` | `72dcefcf05482288e4a353f7697678b9f82f7648078e223345eb3f6444b82c71` |
| `runs_misc/calibra_run/artifacts/frozen_rna_targets.npz` | `d526a7adc7456ac4f0e5e3ff71c0ef2bac96dc8488435ea714ba9840d8b51fb2` |

Three files on the box are called `diagnostic_full_seed42.npz`. Two are byte-identical
(`runs_misc/calibra_run/` and `runs/v21_release_20260720_retry3_resume_safe/`, both `72dcefcf…`); the
third, `runs/v22_a10_11v21_20260725/` (`7674145216572fad…`), is a different file and is pinned in the
digest registry precisely so that citing it is caught. This is the same three-copies-one-filename
hazard the 22:10 entry found for `d2_h_seed42.npz`, in a second place.

### 2. The numbers were not reproduced

Primary reading, as predeclared — `wsi_biology`, `confound_certificate.lda_oof_balanced_accuracy`,
standardised as `certify_axes` standardises, 5 class-stratified folds, seed 42, shrinkage 0.1;
adjustment `residualise.cross_fitted_residuals` against the 99-column design:

* **raw 0.733915** against a published **0.463** — off by 0.271, a factor of 1.59;
* **adjusted 0.030813** against a published **0.035** — off by 0.004;
* **chance 0.047619**, published 0.048 — this, and only this, reproduces.

All seven declared states, joint LDA raw → adjusted:

| state | raw | adjusted |
|---|---:|---:|
| **wsi_biology** | **0.733915** | **0.030813** |
| wsi_identity | 0.818565 | 0.013768 |
| full_biology | 0.870859 | 0.035488 |
| full_identity | 0.909204 | 0.015169 |
| full_patient | 0.861133 | 0.036934 |
| rna_biology | 0.859967 | 0.081892 |
| rna_identity | 0.891230 | 0.080151 |

**`full_biology`'s adjusted value is 0.035488, which rounds to the published 0.035.** It is recorded
because it was measured, and it is **not** an identification. The predeclaration named `wsi_biology`
as the primary state before any of this was computed, for the reason that F1 is about `wsi_biology`
and the check exists to validate that state's adjustment. Picking `full_biology` after the fact
because its second number lands well — while its *first* number, 0.870859, misses 0.463 by more than
any other state in the table — would be the artifact-shopping move in a new costume. The pair still
does not reproduce under it.

### 3. Why it could never have been traced

* The pair was published on **2026-07-30** in the message of commit `4c7166b`.
* `v2/calibra/confound_certificate.py`, which defines `balanced_accuracy`,
  `lda_oof_balanced_accuracy` and `nearest_class_mean_oof` — the *only* balanced-accuracy functions in
  the repository — was created on **2026-08-02** (`942d3c2`), three days later.
* `run_calibra.py` has never emitted a balanced accuracy of any kind; the metric appears in none of the
  122 data rows of `runs/calibra_v2_local/task_rows.csv` nor anywhere in its `calibra_summary.json`.
* No file in git history, live or deleted, computes a cancer-type balanced accuracy on this cohort.
  Checked `--diff-filter=D` across 2026-07-25 → 2026-08-06.
* A `find` over `/home/ubuntu` and `/lambda/nfs/geeg/biorag3_persistent_20260711` returns no run
  output recording 0.463 or 0.035 as a balanced accuracy.

It was a session probe that never persisted, written against code that no longer exists. That is
Note 5 to future agents ("*Where a number came from a session probe that did not persist to disk, say
so explicitly*") being violated at the moment the number was created — not a later transcription
error.

### 4. What changed

**Withdrawn and replaced at all five sites**, with the withdrawn figure left visible rather than
deleted:

* `v2/research/rebase/nature/PHASE1_RESULT.md` — the origin. The `PROVENANCE UNRESOLVED` block becomes
  a withdrawal block that states the cohort, both hashes, the measured replacement and why the pair
  cannot be traced; the validity-check bullet now reads 0.734 → 0.031 with the seven-estimator table
  beside it; the header gains the artifact + digest and the 14 → 11 dev-cancer correction.
* `paper/P1_CALIBRA_DRAFT.md` — abstract (§0), §1.4 contributions, §4.2 and the conclusion. §4.2 gains
  the artifact hash block, the withdrawal blockquote and the estimator table.
* `paper/P1_FIGURES.md` — the artifact hash table gains the Phase 1 cohort row; panel (d) now plots
  the measured pair with the nonlinear overlay, and its caption must record the withdrawal.

**New, so the defect cannot recur in the same place:**

* `v2/research/rebase/nature/p1_cancer_type_certificate.py` — the check as a module, importing every
  statistic, writing its answer to a file. The original defect was not a wrong number; it was a number
  with nowhere to live.
* `v2/research/rebase/nature/p1_cancer_type/out/P1_CANCER_TYPE_CERTIFICATE.json` — the artifact, with
  both input digests, the cohort description, the estimator settings and all seven states.

**`v2/tests/test_paper_artifact_digests.py`:** three artifacts added to `ARTIFACTS` (the Phase 1
cohort, its target table, and the third same-named decoy); `0.7339` and `0.0308` pinned to the cohort
artifact; `diagnostic_full_seed42.npz` added to `HASH_REQUIRED_BASENAMES`; `UNIDENTIFIED_MARKERS`
extended with *withdrawn* / *could not be reproduced* / *unreproducible*, because a withdrawn number
that stays quoted needs the same guard an unattributed one does. `0.463` **keeps** `artifact: None`
and its `unidentified_because` is rewritten to record the failed regeneration, so
`test_the_registry_is_not_vacuous` is untouched and no assertion had to be relaxed.
`v2/tests/test_paper_paths_resolve.py` gains `runs_misc` as a declared box tree and
`diagnostic_full_seed42.npz` as a hash-pinned basename.

### In plain terms

The paper claimed that after we statistically remove cancer type from the image features, a
classifier can no longer tell which cancer a patient has: it drops from getting 46.3% right to 3.5%,
where random guessing gets 4.8%. Nobody could find the file that number came from. We have now found
the *group of patients* it was measured on — 2,530 of them, and everything about that group matches
the paper exactly, right down to the number of hospitals — but re-running the check on them gives
73.4% before and 3.1% after, not 46.3% and 3.5%. Six different classifiers were tried and not one
produces both of the published figures. So the old number is withdrawn and the new one, which is
written to a file with a fingerprint attached, replaces it.

The correction is not cosmetic and it does not simply flatter us. The image features carry far more
cancer information to begin with than we said (73% not 46%), so the removal step is doing much more
work than advertised. But that only holds for the *kind* of classifier the check uses — one that looks
at averages. Ask a classifier that can see shapes rather than averages, and it still identifies the
cancer four to six times better than chance after the removal. Both of those are worth saying out
loud, and the old number said neither.

### Meaning for the claim

P1's first headline — "the adjustment is verified rather than assumed" — **stands, and its evidence is
now on disk with a digest** rather than in prose. But the number that carried it was wrong in both
endpoints, and the corrected picture is sharper in two directions at once: the first-moment removal is
steeper than published (23.8×, to below chance), and what the certificate fails to cover is larger
than published (3.7–6.2× chance under three nonlinear families). The already-recorded sentence that
"cancer is gone" is refuted survives; what changes is that the *evidence for how much is left* is now
measured on this cohort rather than carried over from another.

**Falsifier:** the same cohort, artifact and estimator returning materially different values on a
re-run — the certificate is deterministic given the seed, so any movement is a defect. Or a persisted
2026-07-30 artifact turning up that records 0.463 as a cancer-type balanced accuracy, which would
supersede this reconstruction outright.

**Not claimed:** that the July probe used any particular estimator. Fifty-six readings were taken and
none reproduces the pair; the nearest single value to 0.463 is k-NN's raw 0.4447, and under k-NN the adjusted value is
0.1766, five times the published 0.035. The estimator remains unknown and is recorded as unknown.

### Suite

Run on the box, `~/venv`, threads capped to 1, `pytest morpheus/v2/tests morpheus/tests -q`:

| | result |
|---|---|
| baseline, commit `6a2392f` (predeclaration only, workspace `~/ws_p1base`) | `551 passed, 27 errors in 62.70s` |
| this work, commit `c3a25d9` (workspace `~/ws_p1prov`) | `551 passed, 27 errors in 62.24s` |

**Delta is zero**, and that is the expected number: no test function was added. The changes to
`test_paper_artifact_digests.py` and `test_paper_paths_resolve.py` are registry and allowlist data,
which run inside the existing parametrised tests. The 27 errors are `test_p2_figures` on missing
matplotlib and are the known, expected condition in a checkout.

One real failure was caught by the suite and fixed rather than allowlisted:
`test_every_repository_path_cited_in_a_draft_exists` rejected a bare `frozen_rna_targets.npz` in
§4.2, because the rule is that a bare filename must be unique inside the repository and this one
exists only on the box. Written with its `runs_misc/` directory instead (`c3a25d9`). This is the
guard behaving exactly as its docstring says it should — *"If this test fails, fix the citation. Do
not add the path to the allowlist"* — on the very commit that extends it.

### Files / commits

- `v2/research/rebase/nature/p1_cancer_type_certificate.py`
- `v2/research/rebase/nature/p1_cancer_type/out/P1_CANCER_TYPE_CERTIFICATE.json`
- `v2/research/rebase/nature/PHASE1_RESULT.md`, `paper/P1_CALIBRA_DRAFT.md`, `paper/P1_FIGURES.md`
- `v2/tests/test_paper_artifact_digests.py`, `v2/tests/test_paper_paths_resolve.py`
- Predeclaration `NOTEBOOK_ENTRIES/PREDECLARED_p1_cancer_type_balanced_accuracy_20260805T0130Z.md` (`6a2392f`)
- Prior: `NOTEBOOK_ENTRIES/p1_42_artifact_identified_by_hash_20260804T2210Z.md`
