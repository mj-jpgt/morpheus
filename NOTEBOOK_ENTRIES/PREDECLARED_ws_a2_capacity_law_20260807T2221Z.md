## PREDECLARED — WS-A2: the capacity law, floor(k) vs channel(k) (complaints #2 and #6)

**Logged:** 2026-08-07 22:21 UTC, by agent W1-A. **Status: predeclared, not yet measured.** Committed
before `v2/calibra/capacity_sweep.py` is run against any real artifact.

**This is the highest-value item in the whole revision, and the honesty constraint on it is absolute:
the full sweep ships whatever it says. If the channel does not clear the floor at matched capacity,
that is the result, it leads the report, and `k` is never tuned after the fact to find a favourable
value.**

---

## 1. What is being measured

`v2/calibra/capacity_sweep.py` sweeps `n_components ∈ {1, 2, 4, 8, 16, 32, 64, 128}` on **identical
data and folds** (same cohort, same confound design, same seed, same train/test split machinery) and
reports, per `k`:

- the calibrated **detection floor** from `calibration.spike_recovery_curve` / `floors_from_recovery`
  (imported, not re-derived);
- the within-strata **permutation-null p95** of the channel at that `k`, from
  `calibration.permutation_null` (imported);
- the analytic **chance level** from the existing `hest_claims.capacity_floor_prediction(n, k,
  design_rank)` — `2*sqrt(k/n_eff)`, already implemented, not touched;
- the **observed channel** at that `k`: `spectral.heldout_top_cca(x_residual, y_residual,
  n_components=k)` (imported), the same fitted-direction, out-of-fold statistic
  `run_calibra._channel_measurement` already uses to grade `channel_clears_floor`.

## 2. Fixed protocol, so `k` cannot be picked after seeing results

- **Artifact:** `d2_h_seed42.npz` (`runs/d2_final/artifacts/` on the training box), the same artifact
  the paper's headline Track 1 numbers are drawn from.
- **Targets:** `frozen_rna_targets.npz`, RANDOM_CONTROL columns excluded (as `run_calibra.py` does).
- **State:** `wsi_biology`.
- **Partition:** `test` (n≈2,530), matching the paper's published headline cell (detection floor 0.30
  at `k=32`, §4.2 of `P1_CALIBRA_DRAFT.md`) — chosen specifically so `k=32` in this sweep is a direct
  consistency check against an already-published number, not a free choice.
- **Design:** cancer + pooled tissue-source-site, `min_site_count=10` (the paper's "anchor" design).
- **Seed:** 42 primary; repeated at 43 and 44 if the primary run's verdict is favourable to the channel,
  per PROJECT_GUIDE §2 rule 3 (push a favourable result until it breaks or repeats confirm it, not bank
  n=1).
- **Levels grid:** `0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50`; `n_draws=40`
  (matches Track 2's convention); `n_permutations=200`.

## 3. The falsifier — stated before any `k` is run

**Primary question:** at every `k` in the grid, does `heldout_top_cca(k) > detection_floor(k)`
(`calibration.channel_clears_floor`, magnitude comparison, imported)?

- **A pass at every `k`** is the favourable result, and per the binding honesty rule it is reported
  with the FULL curve shown, not just the passing tail. A pass at large `k` alongside a fail at small
  `k` (or vice versa) is not "the channel clears its floor" — it is "the channel clears its floor at
  capacities ≥ X", a qualified claim, and the qualification is not optional.
- **A fail at any `k`** is reported as a fail at that `k`, full stop, in the same table as every pass.
- The chance level `2*sqrt(k/n_eff)` is plotted alongside but is **not** the bar the channel is graded
  against — the detection floor is. The chance level exists so a reader can see how much of the floor
  is "you'd expect this from `k` alone" versus confound-specific.

## 4. What would make me distrust a FAVOURABLE result

- If the channel clears the floor at every `k` **only because the floor was measured with too few
  draws to resolve small effects** (a floor that reads `NaN`/censored at small `k` is not "cleared",
  it is unmeasured — must be reported as censored, not silently treated as 0).
- If `heldout_top_cca` at large `k` (64, 128) is inflated by the same capacity effect
  `hest_claims.capacity_floor_prediction` predicts for chance alone — i.e., if the *margin* between
  channel and floor does not grow with `k` in a way distinguishable from the chance-level curve also
  growing, the "clears the floor" verdict at large `k` is riding on capacity, not signal, and that must
  be stated even if the raw magnitude comparison passes.
- If the detection floor **rises with `k`** as fast as or faster than the channel does (i.e. `floor(k)`
  tracks `chance(k)` and the channel does too), that is evidence the whole comparison is capacity vs
  capacity, not signal vs noise, and is the headline finding regardless of which side numerically wins.

## 5. What would make me distrust an UNFAVOURABLE result (does not clear at some k)

- A single non-clearing `k` on `seed=42` alone is not reported as "the channel fails at capacity `k`"
  without at least the two repeat seeds (43, 44) named in §2, run precisely because a single-seed
  failure could be as unrepresentative as a single-seed pass. Both directions get the same scrutiny.

## 6. Scope note

Sibling predeclarations: `PREDECLARED_ws_a1_induced_correlation_theory_20260807T2221Z.md`,
`PREDECLARED_ws_a6_exchangeability_20260807T2221Z.md`. This file covers WS-A2 only.
