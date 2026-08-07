## PREDECLARED — WS-A6: nested exchangeability blocks + Freedman-Lane (complaint #11, method half)

**Logged:** 2026-08-07 22:21 UTC, by agent W1-A. **Status: predeclared, not yet measured.** Committed
before the extended `calibration.permutation_null` / new `calibration.freedman_lane_null` are run
against any real artifact.

---

## 1. What is being built and why

`calibration.permutation_null` currently blocks within cancer type only. Two additive extensions:

1. **`permutation_null` gains an optional `nested_strata` keyword.** When supplied (pooled
   tissue-source-site, i.e. site-within-cancer), the permutation happens within the intersection of
   `strata` (cancer) and `nested_strata` (site) — a genuinely finer exchangeability block than cancer
   alone, matching Winkler et al. 2015's "sites nested within a scanning/acquisition factor" structure.
   A combined cell with fewer than 2 members cannot be permuted internally (nothing to swap); those
   rows are left in place and the count of such rows is reported explicitly rather than silently
   folded into the null as if they had been permuted. **`nested_strata=None` reproduces the exact
   existing behaviour** — this is additive, not a rewrite, per my file-ownership scope.
2. **New function `freedman_lane_null`.** The existing scheme permutes the RAW `y` rows (within
   strata) and re-fits `cross_fitted_residuals(y[order], design)` from scratch every permutation — this
   is closer to a "permute-then-refit" scheme. Freedman & Lane (1983) instead permutes the RESIDUAL
   `e = y − ŷ(design)` directly (computed once, real, unpermuted) and scores
   `top_canonical_correlation(x_residual, e[order])` — no per-permutation refit. For linear
   (ridge/OLS) residualisation this is algebraically the standard equivalent of "permute the reduced-
   model residual, add back the reduced fit, refit the full model" (the refit changes nothing beyond
   what the residual permutation already captures, because the design's column space is what is being
   projected out both times) — cited in Winkler et al. 2014 Table 1 as one of several exchangeable
   schemes for a GLM nuisance design. It also accepts `nested_strata` for the same nested-block
   comparison as (1).

Both import `cross_fitted_residuals` (`residualise.py`) and `top_canonical_correlation` (`spectral.py`)
unchanged; neither reimplements a canonical statistic.

## 2. What is measured, and against what real artifact

Same cohort/design/state as WS-A2 (`d2_h_seed42.npz`, `wsi_biology`, `test` partition, cancer+TSS
`min_site_count=10` design, `n_components=32` to match the paper's already-published headline cell).
Four nulls compared on the identical `(x_residual, y_residual)`:

1. existing `permutation_null`, `strata=cancer` only (today's shipped behaviour, unchanged);
2. `permutation_null`, `strata=cancer`, `nested_strata=pooled_site`;
3. `freedman_lane_null`, `strata=cancer` only;
4. `freedman_lane_null`, `strata=cancer`, `nested_strata=pooled_site`.

`n_permutations=200` for all four (matched, so a null_p95 comparison is apples-to-apples).

## 3. The falsifier — stated before any null is run

**Question:** does the null move enough to change any verdict currently reported in the paper?
Specifically: does `null_p95` under (2), (3), or (4) exceed the paper's published `observed_top_cca`
for the same cell (§4 of `P1_CALIBRA_DRAFT.md`, `wsi_biology` d2_h in-sample top-CCA), which would flip
a currently-significant channel to non-significant under permutation?

- **Bar for "the null moved, and it matters":** `null_p95` under any of (2)/(3)/(4) differs from (1) by
  more than 20% relative, OR the permutation `p` value crosses the conventional 0.05 boundary relative
  to (1) for the same observed statistic. Either condition triggers a full paragraph in the writeup,
  not a single sentence.
- **If none of (2)/(3)/(4) differ from (1) by more than 20% relative and no `p` crosses 0.05:** report
  this as the robustness result the spec anticipates ("if neither moves the null, say so — that is a
  robustness result and it is worth a sentence, not a section") — plainly, with the actual numbers
  shown in a small table, not asserted without the comparison visible.

## 4. What would make me distrust a FAVOURABLE (robustness) result

- If the nested-block singleton-cell fallback rate is high (many site-cells too small to permute), a
  "the null didn't move" finding could just mean nested blocking degenerated back toward cancer-only
  blocking for most rows, not that nesting is genuinely inert. The fallback rate is reported alongside
  the null comparison, not omitted.
- If Freedman-Lane's null is suspiciously IDENTICAL (not just similar) to the existing scheme, that is
  a red flag for a bug (e.g. accidentally permuting the same object twice) rather than evidence of
  robustness, and must be checked by hand on a small synthetic fixture before being reported as
  agreement.

## 5. What would make me distrust an UNFAVOURABLE (null moves) result

- A shift driven entirely by very few extreme permutation draws (`n_permutations=200` is modest); if
  `null_max` is a large outlier relative to `null_p95`/`null_median` under the moved scheme, the p95
  estimate itself is noisy and the comparison should be re-run at a larger `n_permutations` before
  being reported as a settled finding.

## 6. Scope note

Sibling predeclarations: `PREDECLARED_ws_a1_induced_correlation_theory_20260807T2221Z.md`,
`PREDECLARED_ws_a2_capacity_law_20260807T2221Z.md`. This file covers WS-A6 only. All edits to
`calibration.py` are additive (new function + new optional kwarg with a no-op default); no existing
call site's behaviour changes, verified by running the existing `v2/tests/test_calibra.py` suite
unchanged before and after.
