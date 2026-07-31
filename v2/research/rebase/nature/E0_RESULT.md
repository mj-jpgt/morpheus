# E0 / E0b — basis transfer. **The engine survives, at ~10% of the achievable ceiling.**

**Run:** `runs/e0_20260731`, A100-40GB, 901 s wall, peak GPU 4.5 GB.
**Code:** commit `24d1bff`, clean worktree, provenance symlink verified.
**Scale:** 1,000 Haar draws, 1,000 bootstrap resamples, k ∈ {10, 25, 50, 100}, 2 RNA transforms.
**Gates:** 248 PASS, 89 OBSERVED, **0 FAIL.**

---

## 0. The result

| context | k | responsive | control | gap | paired 95% CI | Haar p95 | % of ceiling |
|---|---:|---:|---:|---:|---|---:|---:|
| K562 signed | 10 | 0.0822 | 0.0095 | +0.0727 | [+0.048, +0.068] | 0.0018 | 10.5% |
| K562 signed | 25 | 0.0743 | 0.0211 | +0.0532 | [+0.041, +0.053] | 0.0039 | 10.6% |
| K562 signed | 50 | 0.0714 | 0.0308 | +0.0405 | [+0.033, +0.041] | 0.0074 | 10.7% |
| K562 signed | 100 | 0.0725 | 0.0400 | +0.0325 | [+0.026, +0.031] | 0.0144 | 11.1% |
| K562 clip | 10–100 | 0.0833–0.0732 | 0.0095–0.0401 | +0.074 → +0.033 | all > 0 | — | 10.3–11.1% |
| **RPE1 (both)** | all | — | — | — | **UNDECIDABLE** | — | — |

**Verdict: `supported`** — 2 of 4 contexts decided, both K562. Responsive perturbations align with
tumour expression structure more than non-responsive ones at every k, in both transforms, with the
paired bootstrap CI entirely above zero.

**H3 is not refuted. The causal-dictionary engine may proceed.**

## 1. The control arm was not a formality — it absorbed most of the signal

**The non-responsive control scores 2–5× above the random floor** (0.0095 vs 0.0018 at k=10;
0.0400 vs 0.0144 at k=100), and its share of the raw alignment *grows* with k: at k=100 the control
accounts for **55%** of the responsive arm's total overlap.

That is the audit's warning made quantitative. The original Haar-only rule would have credited the
entire 0.0725 to perturbation biology. **Roughly half of it is generic expression structure that
perturbations with no detectable transcriptional effect reproduce just as well.**

The `haar_only_would_say` field is `True` everywhere — the old design would have said "proceed" here
too, but it would also have said "proceed" on data with no shared biology at all. The two rules agree
on this dataset; only one of them was capable of disagreeing.

## 2. Do not overclaim the size

- **Normalised alignment is 10.3–11.1% of the ceiling** (the TCGA self-split maximum, same k, same q,
  same estimator). The transferable structure is real and replicated, and it is **small**.
- Absolute overlaps are 0.07–0.08. This is not a "cell lines recapitulate tumours" result.
- The correct claim is: *a minority component of tumour expression geometry is shared with genuine
  perturbation responses, distinguishable from both random directions and from perturbations that
  did nothing.*

## 3. Not a library-size artifact

`pc1_share` is **0.089–0.214** — stripping PC1 removes only 9–21% of the alignment. Had the result
been mean-expression/library-size, this would have been near 1. The offset sweep {0,1,2,5} is in the
JSON. The headline statistic is computed with PC1 already removed.

## 4. **The main limitation: no cross-lineage replication**

RPE1 has only **50** non-responsive perturbations (against the estimator's 151-row minimum), so both
RPE1 contexts are `unavailable_insufficient_rows_50_lt_151` and are reported as
`contexts_undecidable`. They do not veto the K562 result, but they do not support it either.

**Agreement between two unrelated lineages was the evidence that the alignment is biological rather
than K562-specific. That evidence does not exist.** Every claim must currently be scoped to K562.
Options: lower `nonresponsive_p` to 0.2 for RPE1 (168 rows, but a contaminated control raises the
floor — errs safe), or add a third cell context.

## 5. E0b — the dictionary is ~132 effective dimensions, not 11,000

| | perturbations used | numerical rank | **effective rank** | stable rank | coherence |
|---|---:|---:|---:|---:|---:|
| K562 GWPS | 8,403 (of 11,258 raw) | 8,141 / 8,246 | **132.1** | 17.4 | 0.846 |
| RPE1 | 2,326 | 2,325 / 2,326 | **113.9** | 7.8 | 0.900 |

Computed in float64 with a `s[0]·max(n,p)·eps` cutoff — the previous float32 `1e-10` cutoff never
fired and reported 799 for a matrix of true rank 50.

**This settles the banned claim.** ~11,000 measured perturbations occupy an effective rank of ~132,
with a stable rank of 17. Convergent perturbations are biological equivalence under the measurement,
not a defect — but the catalogue's resolution ceiling is set by this quotient, and it is two orders of
magnitude below the perturbation count. Max off-diagonal atom correlation is 0.85–0.90.

*Caveat:* `n_equivalence_classes` returned 8,403 and 2,326 (i.e. n), which cannot be right alongside a
coherence of 0.85 and an effective rank of 132. The clustering threshold is almost certainly
mis-specified. **This number is not usable and must be recomputed before any catalogue claim.**

## 6. Provenance
- 2,341 of 11,258 K562 rows (21%) dropped as non-finite. Loss is unbiased w.r.t. the arm split
  (responsive 21.2%, non-responsive 20.7%), so it does not confound the control — but it is currently
  ungated and should be.
- 7,094 shared genes (K562↔TCGA), 7,621 (RPE1↔TCGA); zero constant columns; all five housekeeping
  genes survive; registry coverage 0.927; 29 cancer strata, smallest 28.

## 7. Known limitations of the statistic
- **Sign-blind.** `svdvals(Vaᵀ Vb)` is a subspace statistic: an *anti*-aligned perturbation response
  scores identically to an aligned one. E0 tests shared **directions**, never shared **sign**.
- `bootstrap_ci95` (marginal, per-arm) is **biased low by 0.04–0.09** and is flagged
  `bootstrap_ci95_is_biased_low`. The decision uses the *paired* difference, where the bias cancels.
  Never quote the marginal interval as a 95% CI.
- The responsive subsample is a single seeded draw of 956 from 3,132; the two transforms use different
  subsamples, so cross-transform agreement doubles as a 2-point stability check.
