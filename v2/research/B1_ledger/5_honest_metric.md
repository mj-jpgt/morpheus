# Lane 5 — Confound-aware / honest evaluation metric harness

## Queries/searches run
- WebSearch: "HESCAPE benchmark histopathology gene expression prediction batch effect pitfalls spatial transcriptomics"
- WebSearch: "TCGA batch effects tissue source site confounder machine learning site-preserved cross-validation Howard 2021"
- WebFetch: arxiv 2508.01490 (HESCAPE, PDF — methodology body not in excerpt)
- WebFetch: PMC11814321 (Nat Commun 2025 spatial-expression benchmark — HVG/SVG + sparsity confounder recommendation)
- Grep (codebase): `macro_cancer_pearson|within_cancer|per_cancer`, `RANDOM_CONTROL__|build_matched_random_controls`, `ComBat|combat|limma|batch.*effect|site.*preserv`, `site|tss|batch|plate|center|submitter`, `pearson_metric|paired_bootstrap`
- Read: paired_bootstrap.py, v21_evaluation.py, discovery_targets.py, select_v21_profile.py, build_paired_split.py, tests/test_v21_evaluation.py

## Sources
Web:
- Howard et al., "The impact of site-specific digital histology signatures on deep learning model accuracy and bias", Nat Commun 12:4423 (2021) — https://www.nature.com/articles/s41467-021-24698-1 — TCGA submitting-site batch effect; recommends never train+validate on same site; report per-site variation.
- HESCAPE, "A Large-Scale Benchmark of Cross-Modal Learning for Histology and Gene Expression in Spatial Transcriptomics", arXiv:2508.01490 — https://arxiv.org/pdf/2508.01490 — contrastive pretraining *degrades* direct expression prediction vs baseline encoders; **batch effects identified as the key factor interfering with cross-modal alignment**; calls for batch-robust evaluation.
- Nat Commun 2025 spatial-expression benchmark — https://pmc.ncbi.nlm.nih.gov/articles/PMC11814321/ — pooled correlation across all genes is masked by low/zero-expression genes; recommends restricting to variable genes and flags sparsity as the dominant confounder ("average masking the identification of biologically meaningful concordance").
- El Nahhas/Kather et al. "Multimodal deep learning: improvement or reflection of batch effect?", Cancer Cell (2022) — https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00522-0 — multimodal prognostic gains can be batch-effect artifacts.
- "Auditing Data Leakage in Whole-Slide Image Multimodal Benchmarks", arXiv:2607.12278 — https://arxiv.org/html/2607.12278v1.

Code (file:line):
- `paired_bootstrap.py:8` `pearson_metric` (POOLED); `:24-58` `paired_bootstrap_difference` with `mode="patient"|"cancer"` (cancer = cluster resample); `:61-64` `paired_patient_and_cancer_bootstrap`.
- `v21_evaluation.py:188-191` `_macro_cancer_pearson` (per-cancer mean of within-cancer pearson) — computed but only emitted as one report row (`:225`).
- `v21_evaluation.py:210,215-227` `_molecular_rows` fits Ridge on development, reports pooled `pearson` first; macro is a sibling row.
- `v21_evaluation.py:439,447-453` `_control_comparison_rows` computes biology-minus-matched-control gap on POOLED `pearson` only.
- `v21_evaluation.py:553` paired-vs-teacher bootstrap is called with `pearson_metric` (POOLED).
- `select_v21_profile.py:46-54` promotion gate: biological-minus-control delta computed on `metric == "pearson"` (POOLED); this is THE primary-claim gate.
- `discovery_targets.py:115-175` `score_gene_signatures` default `competitive_rank` = within-sample percentile (library-shift robust); `:178-271` `build_matched_random_controls` = size-matched (n = member count) controls matched on train-only mean/log-var/PC1, biology genes excluded, deterministic seed; `:330-411` `fit_train_only_standardizer` residualizes each target by **train-only per-cancer offset** (cancer-safe; unseen cancer -> offset 0).

## Findings
The harness already implements 3 of the 4 requested controls; the mechanism failure is that the PRIMARY reported/bootstrapped/gated quantity is **pooled Pearson**, which is exactly the confounded statistic the literature warns against.

1. **Within-group correlation exists but is not primary.** `_macro_cancer_pearson` (v21_evaluation.py:188) is the correct honest metric (average of within-cancer Pearson), but it is a single sibling row. The promotion gate (select_v21_profile.py:47), the paired bootstrap (v21_evaluation.py:553), and the control-gap (v21_evaluation.py:439) all key off pooled `pearson`. Pooled Pearson over 21 held-out cancers is inflated by cross-cancer tissue-type separation (a WSI feature that trivially tracks cancer type will correlate with any cancer-varying programme). This is the documented ~50% cross-cancer / +0.07 confound in the lane brief and the HESCAPE batch-effect mechanism.

2. **Size-matched random-gene null already correct** (build_matched_random_controls, discovery_targets.py:178) — n-gene-matched, covariate-matched, biology-excluded, train-only. No change needed to construction. It just needs to be differenced against the *macro* metric, not pooled, so the "specific over random" claim is itself within-cancer.

3. **Paired bootstrap already correct** (paired_bootstrap.py) — patient IID and cancer-cluster modes both present. Only the metric *callable* passed to it is wrong (pooled). A cancer-aware metric callable that ignores per-resample cancer labels cannot be threaded through the existing `metric(y[idx], pred[idx])` signature, so the metric must carry cancers via closure over the resampled indices — handled below.

4. **Batch/site (TSS) correction is genuinely absent** and, per Howard 2021, ideally wanted. BUT: (a) no TSS/site column exists in the split, metadata, or artifacts anywhere in v2 (grep found none); (b) the held-out-cancer protocol already prevents same-cancer train/test leakage; (c) RNA targets use within-sample competitive_rank + train-only cancer residualization, which removes the library-size and per-cancer-mean batch axes that ComBat/limma would target. Full ComBat/limma on H-Optimus patch features is a large, scaling-sensitive change (fit-on-train-fold ComBat over 21 cancers, must be re-fit inside every bootstrap to stay honest) that is out of proportion to the mechanism fix. **Recommendation: make within-cancer the primary metric now (cheap, correct, closes the confound the data actually shows); defer feature-space ComBat to a separate lane, and add TSS extraction from the TCGA barcode only if a site confound is later demonstrated.** Within-cancer macro already absorbs the dominant confound axis (cancer identity); site-within-cancer is second-order and unmeasurable without a TSS column.

## Recommended change (file:line, exact)
Smallest change that fixes the mechanism: make the within-cancer (macro) Pearson the primary, bootstrapped, and gated statistic. Four edits, one new helper, one test.

1. **New cancer-aware metric factory** — add to `paired_bootstrap.py` (after `pearson_metric`, ~line 14):
```python
def macro_group_metric(base: Callable[[np.ndarray, np.ndarray], float], groups: np.ndarray):
    """Return a metric(actual, predicted) that averages `base` within each group.

    `groups` is indexed by the SAME rows as the vectors the returned metric will
    receive, so it composes with the resampling in paired_bootstrap_difference:
    resample the row indices, then this closure re-partitions the resampled rows
    by their carried group label. Groups with <3 finite pairs are dropped.
    """
    groups = np.asarray(groups).astype(str)
    def _metric(actual, predicted):
        vals = [base(actual[groups == g], predicted[groups == g]) for g in np.unique(groups)]
        vals = [v for v in vals if np.isfinite(v)]
        return float(np.mean(vals)) if vals else float("nan")
    return _metric
```
Note: because `paired_bootstrap_difference` resamples an `indices` array and then calls `metric(y[indices], candidate[indices])` (paired_bootstrap.py:55), the group vector must be resampled in lockstep. Thread it by passing `cancers[indices]` into a per-iteration metric OR (cleaner, no signature change) build the closure over the ALREADY-resampled `cancers[indices]` inside the loop. Concretely, add a `group_metric: bool = False` branch in `paired_bootstrap_difference` (paired_bootstrap.py:55) that, when a `cancers` array is present, computes `m = macro_group_metric(metric, cancer[indices])` and uses `m(y[indices], ...)`. This keeps patient-mode IID resampling but scores it within-cancer — the honest paired delta.

2. **`v21_evaluation.py:553`** — pass the macro metric to the paired bootstrap so the primary teacher-vs-challenger delta is within-cancer, not pooled. Keep pooled as a secondary reported row for transparency.

3. **`v21_evaluation.py:439`** — in `_control_comparison_rows`, difference biology-vs-matched-control on `metric == "macro_cancer_pearson"` (add it alongside the existing pooled row so the "specific over random gene sets" claim is within-cancer). This directly tests the lane brief's +0.07 number under the honest metric.

4. **`select_v21_profile.py:46-53`** — change the promotion-gate metric filter from `metric == "pearson"` to `metric == "macro_cancer_pearson"` for both `biological` and `controls`. This makes the primary V2 claim gate = "within-cancer curated-biological Pearson exceeds size-matched random-gene controls, and R@10 within teacher tolerance." One-line filter change; the delta/eligibility logic is unchanged.

5. **Test** — extend `tests/test_v21_evaluation.py` (after line 54): assert a `macro_cancer_pearson` row exists for `molecular_prompting`, assert the paired-vs-teacher rows include a macro variant, and add a unit test in `tests/test_data_evaluation.py` that `macro_group_metric` over 2 groups equals the mean of the two within-group pearsons and returns NaN when every group has <3 pairs. Also add a regression asserting `_control_comparison_rows` emits the macro-based `mean_matched_biological_minus_control_pearson`.

Report BOTH pooled and within-cancer everywhere (never drop pooled), but define the PRIMARY/headline and the promotion gate as within-cancer. This is the exact HESCAPE + Howard recommendation instantiated with the metric the codebase already computes.

## Risks & scaling
- **Metric can go NaN / thin per cancer.** Within-cancer Pearson needs >=3 finite pairs per cancer (already enforced by `_correlation`). With 21 held-out cancers and small per-cancer test n, some cancers drop each bootstrap resample; the macro mean is over surviving cancers, so variance is higher than pooled. Mitigate: keep repeats=2000 (already default), and report `n_valid` cancers per bootstrap (the `_interval` helper already reports n_valid). Do not fall back to pooled silently.
- **Cancer-cluster bootstrap + macro interaction.** In `mode="cancer"` (paired_bootstrap.py:52-54) a cancer may be sampled 0 or multiple times; the macro-over-resampled-labels closure handles duplicates correctly (they weight that cancer more, which is the intended cluster-uncertainty behaviour). Verify the closure uses `cancer[indices]` (resampled), not the original vector.
- **Gate becomes harder to pass.** Switching the promotion gate to within-cancer will likely make `full` fail more often (the confound was inflating it). That is the point — a failed gate is the honest negative result the safeguard at select_v21_profile.py:86 already anticipates. No code risk, but stakeholders must expect lower headline numbers.
- **Deferred ComBat/limma.** If a site confound is later suspected, TSS is recoverable from the TCGA barcode (chars 6-7 of `TCGA-XX-....`), addable at canonical_registry.py where `bcr_patient_barcode` is already handled (canonical_registry.py:28). Feature-space ComBat must be re-fit on train fold and, for honest CIs, re-fit inside each bootstrap — O(repeats x cancers x dims); on A10 this is minutes-to-hours and belongs in its own lane, not this metric fix. Within-cancer macro is the 90% fix at ~0 extra compute.
- **No new heavy deps.** All edits are numpy/pandas already imported. CPU test suite impact negligible (adds ~1 unit test).
