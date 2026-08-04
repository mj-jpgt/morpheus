## PREDECLARATION — regenerate `0.463 → 0.035` (cancer-type balanced accuracy, chance 0.048, n = 2,530) or withdraw it

**Logged:** 2026-08-05 01:30 UTC. **Written before any balanced accuracy was computed.**
**Status at time of writing:** the pair is marked `artifact not identified` / `PROVENANCE UNRESOLVED`
at all four P1 sites and pinned with `"artifact": None` in `v2/tests/test_paper_artifact_digests.py`
(entry `p1_42_artifact_identified_by_hash_20260804T2210Z.md`, §4).

### What is already established, before computing anything

The cohort **is** identifiable, and it was found by structure rather than by hunting for a file that
reproduces the number. `v2/calibra/run_calibra.py` selects its cohort as
`split == "test"` in the input artifact, intersected with the patients present in the target table.
On `runs_misc/calibra_run/artifacts/diagnostic_full_seed42.npz`
(SHA-256 `72dcefcf05482288e4a353f7697678b9f82f7648078e223345eb3f6444b82c71`) that mask has exactly

| quantity | this artifact | `PHASE1_RESULT.md` states |
|---|---:|---:|
| patients, `split == "test"`, aligned to targets | **2,530** | 2,530 |
| distinct cancer types on that mask | **21** | 21 |
| `1 / n_classes` | **0.047619** | chance 0.048 |
| confound design columns (`cancer` + pooled `tss`, `min_site_count=10`) | **99** | 99 |
| TSS sites kept | **75** | 75 |
| split totals | train 3,124 / val 538 / test 2,530 = **6,192** | pre-rebuild cohort |

Six of six. The cohort differs from the site arm's 2,766 because that arm runs on the **rebuilt**
maximal split (6,427 patients, entry 2026-08-01 09:05); this is the **pre-rebuild** 6,192-patient
split. That is the whole of the n = 2,530 vs n = 2,766 discrepancy.

The committed run that produced the rest of `PHASE1_RESULT.md` is `runs/calibra_v2_local/`
(commit `4c7166b`, 2026-07-30), whose `task_rows.csv` records `n_patients 2530`,
`n_confound_columns 99`, `n_distinct_sites_kept 75`, and whose `heldout_top_cca` values reproduce the
published table (0.476787 → 0.477, 0.539273 → 0.539, …). **It contains no cancer-type balanced
accuracy of any kind** — the metric is absent from `calibra_summary.json` and from all 123 rows of
`task_rows.csv`. `run_calibra.py` has never emitted such a metric.

**The estimator that exists today did not exist when the number was published.**
`v2/calibra/confound_certificate.py`, which defines `lda_oof_balanced_accuracy` and
`nearest_class_mean_oof`, was created **2026-08-02** (`942d3c2`). The number was committed
**2026-07-30** (`4c7166b`, in the commit message). No script in the repository's history, deleted or
live, computes a cancer-type balanced accuracy on this cohort. So the published pair came from a
session probe that never persisted, written against code that no longer exists — and an exact
reproduction is *a priori* unlikely rather than expected. I am writing that down now so that a
non-reproduction cannot later be presented as a surprise, and a reproduction cannot be presented as
more than a coincidence of two independent implementations agreeing.

### What I will compute

On the mask above, for **every one of the seven declared representation states** in that artifact:

* **raw** — `confound_certificate.lda_oof_balanced_accuracy(scaled, cancer_index, 21)`, joint
  shrunk-covariance LDA, `n_splits=5`, `seed=42`, `shrinkage=0.1`, folds from
  `_stratified_folds`, per-axis standardisation exactly as `certify_axes` does it;
* **adjusted** — the same statistic on `residualise.cross_fitted_residuals(x, design, seed=42)`
  with the 99-column `cancer` + pooled-`tss` design, i.e. the identical adjustment CALIBRA applies
  before it measures the channel;
* **chance** — `1 / n_classes`.

Everything is **imported** from `v2/calibra/`. Nothing is reimplemented inline. (The project has had
five statistic substitutions from exactly that failure mode, and an AST scan test now guards it.)

**Primary state, named in advance: `wsi_biology`.** It is the state F1 is about and the state whose
adjustment the "confound removal verified" check exists to validate. All seven are reported, but
`wsi_biology` is the one whose agreement or disagreement decides the outcome. Choosing the state
*after* seeing which one lands nearest 0.463 would be the same epistemic move the previous agent
declined to make with the artifact, and I will not make it with the state.

Secondary readings, reported but **explicitly not admissible as a reproduction**: the unstandardised
variant, the per-axis `nearest_class_mean_oof` maximum and median, and the k-NN / forest / RBF-SVM
probes from `nonlinear_confound_probe.py`. If one of those matches while the primary does not, it is
a *candidate mechanism* for what the July probe might have run — it is not provenance, and the
withdrawal still happens.

### What I will conclude, per outcome — decided now

**Outcome A — `wsi_biology` reproduces 0.463 and 0.035 (and chance 0.048) to the three decimals
published.** I report a reproduction. I pin path + SHA-256
(`runs_misc/calibra_run/artifacts/diagnostic_full_seed42.npz`, `72dcefcf…`) in the style
`test_paper_artifact_digests.py` enforces, add the artifact to `ARTIFACTS`, flip the `0.463` entry's
`"artifact"` from `None` to that key, delete the `unidentified_because` reason, remove the
`PROVENANCE UNRESOLVED` blocks, and update
`test_the_registry_is_not_vacuous`'s `any(artifact is None)` assertion **as a deliberate decision in
the same commit**, per Note 4 to future agents. I will state that the reproduction is by a *different
estimator implementation* than the one that produced the original, because that is true.

**Outcome B — the numbers come back close but not equal (e.g. 0.461 vs 0.463, or 0.041 vs 0.035).**
This is **not** a reproduction and I will not report it as one. I report both figures side by side, in
full precision, state plainly that they differ, and take Outcome C's action. I will not round toward
the published value and I will not describe the difference as "consistent with".

**Outcome C — the numbers do not come back.** I withdraw. `0.463 → 0.035` is replaced at all four P1
sites plus `PHASE1_RESULT.md` with what I actually measured, attributed to artifact + SHA-256 + the
named estimator, and the prose says the previously published figure could not be reproduced. The
history is not deleted: the original pair stays visible as a withdrawn figure with the date and commit
it came from.

**In all three outcomes**, the direction of the finding is reported before its convenience. If the
adjustment turns out to work *less* well than `0.463 → 0.035` advertises — for example if the adjusted
balanced accuracy sits at or above chance rather than a seventh of it, or if the raw value is far
below 0.463 so that the advertised drop was never that large — that is the result and it goes into the
paper, not into a footnote. §4.2's surrounding prose already concedes that the stronger sentence
("cancer is gone") is refuted and that a nonlinear probe recovers cancer at 3.45× chance; a weaker
first-moment drop is consistent with that concession and is not a new embarrassment, but it does
change the number that gets quoted.

### What would make this entry wrong

* A persisted July-30 artifact recording a cancer-type balanced accuracy on this cohort turning up
  somewhere I did not search — that would supersede the reconstruction outright. Searched: the whole
  of `/home/ubuntu` and `/lambda/nfs/geeg/biorag3_persistent_20260711` on `150.136.45.194`, every
  committed file at HEAD, and every deleted file in git history since 2026-07-25.
* The `split == "test"` mask on a *different* artifact also giving exactly 2,530 / 21 / 99 / 75. That
  would make the cohort identification ambiguous, and I would report both readings rather than pick
  one. (`diagnostic_full_seed42.npz` exists in three places on the box; two share the hash above, the
  third — `runs/v22_a10_11v21_20260725/` — does not, and is checked.)
