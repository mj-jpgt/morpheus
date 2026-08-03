# PRE-REGISTRATION — proliferation-deflated D2 readout (the named remedy for `proliferation_deflation`), and D3 seeds 43/44

**Logged:** 2026-08-03 14:00 UTC. **Written and committed BEFORE the deflated comparison runs.**
Every number below is either a property of an input (annotation-side only, no outcome touched) or a
threshold. Follows `NOTEBOOK_ENTRIES/d3_d2p3_preregistration_20260803T1300Z.md` (`cd9b056`).

**How obtained:** planning only. Inputs inspected: `~/e0_run/d2_v3/bootstrap/D2_*_seed4{2,3,4}.json`,
`~/e0_run/d3/d2_3/axis_table.csv`, `v2/research/rebase/d2_compare.py`, `v2/calibra/claim_guards.py`,
`tests/test_claim_guards.py`.

---

## 1. Which side gets cut, and why

`d2_compare` measures `top_canonical_correlation(wsi_biology_residual, target_residual)` at 16
components for each arm, then bootstraps the **paired** difference `pbs − hallmark`. There are two
things one could cut, and they answer different questions:

- Cutting the **readout target block** (dropping proliferation-ish columns from
  `frozen_rna_targets`) asks "can the arms still predict non-proliferation *targets*". That is **not**
  the blocker. It also destroys the pre-registered headline, whose target set is fixed at the 40
  `heldout_pathway + immune_tme + tumour_state` columns.
- Cutting the **PBS axes** asks "is the arm difference carried by the proliferation content of the
  perturbation basis". That **is** the blocker, whose text is about the responsive arm being selected
  on having an effect and therefore enriching for cell-cycle genes.

**So the cut is on PBS axes, and it is applied by residualisation, not by column deletion.** The
patient-level scores of the proliferation-loaded PBS axes are appended to the confound design, so
`cross_fitted_residuals` removes that subspace from **`y`, `hallmark_wsi` and `pbs_wsi` identically**
before any canonical correlation is taken. The readout target block is left exactly as
pre-registered. This is `claim_guards`' own remedy verbatim — *"re-run with proliferation/cell-cycle
programme regressed out"* — and the symmetry across arms is what makes the difference interpretable.

Implementation: two new `d2_compare` flags, `--deflate-npz` and `--deflate-axes-file`, feeding the
existing `confound_design` numeric path. The axis list is read from a file so that the exact cut is a
durable artifact rather than a shell argument.

## 2. Which proliferation definition governs — declared now

I previously showed these two disagree (Spearman 0.577), and that disagreement is already logged as a
FAIL row in `GATE_LOG.md`. So the governing definition is fixed here, in advance:

- **GOVERNING: `prol_top100`** — proliferation fraction among each axis's top-100 |loading| genes.
  It governs because `prol_wmean` is demonstrably diluted: on a dense basis it cannot distinguish an
  axis that leads entirely on proliferation genes from one that leads entirely on non-proliferation
  genes (pinned by `test_top_k_statistic_separates_a_proliferation_axis_that_the_weighted_mean_dilutes`).
- **CO-PRIMARY: `prol_wmean`** — the column the ledger names. Reported in full either way.
- **CONSERVATIVE: `union`** — the union of both top quartiles. A clean discharge requires this too,
  so the more favourable definition cannot become the headline by default.

Cut sets, computed from the annotation alone (no outcome involved), written to
`~/e0_run/d3/d2_deflate/axes/`:

| set | n axes | contains PBS_001? |
|---|---:|---|
| `prol_top100` (governing) | 33 | **yes** |
| `prol_wmean` (co-primary) | 32 | **yes** |
| `union` (conservative) | 46 | **yes** |
| intersection | 19 | yes |
| `placebo_random` (rank-matched) | 33 | no — drawn from the 82 axes in neither quartile |
| `placebo_random_union_sized` | 46 | no |

**PBS_001 — the 4.5×-background axis that is the most legible axis in 4/6 runs and 2nd in the other
2 — is in every real cut.** It is the test case, and it is removed in all three.

## 3. Rank-matched placebo, again

Deflating 33 axes costs the design 33 columns, and a shrinking gap could be that rather than
proliferation. So each real cut is paired with a **placebo cut of the same size drawn at random from
the 82 axes in neither top quartile** (seed 20260803). Fixed reading:

- real shrinkage ≫ placebo shrinkage ⇒ the change is **proliferation**.
- real shrinkage ≈ placebo shrinkage ⇒ the change is **design rank**, and the deflation has measured
  nothing about proliferation. Binding even if inconvenient.

## 4. What counts as pass / fail

**Baseline** (already on the ledger, `D2_RESULT.md` §2, stratified 40 targets, `pbs − hallmark`):
seed 42 **−0.1325**, seed 43 **−0.1089**, seed 44 **−0.1226**; patient CI₉₅ excludes zero in 3/3.

Primary readout: the same 40 untrained targets. Negative control: the 90 `random_control` targets,
run under the identical deflation. Everything quoted as **paired within-run differences**, never
levels, and baseline vs deflated compared **within the same seed**.

| verdict | condition, on the GOVERNING cut |
|---|---|
| **GAP SURVIVES** | \|Δ_deflated\| ≥ **70%** of \|Δ_baseline\| in **3/3** seeds, patient CI₉₅ excludes zero in **3/3**, **and** the `union` cut also clears both bars |
| **PARTIAL** | retains 40–70%, or CI₉₅ excludes zero in only 2/3, or governing clears but `union` does not |
| **GAP IS PROLIFERATION** | retains **< 40%**, or the patient CI₉₅ includes zero in ≥ 2 seeds |

Negative control requirement: the deflated `random_control` gap must stay far smaller than the
deflated untrained40 gap (undeflated it is −0.0099…−0.0280 against −0.109…−0.133). If deflation
inflates the random-control gap towards the real one, the instrument has been damaged by the
deflation and **no verdict is issued**.

## 5. The scope question, decided in advance so it cannot be fitted afterwards

Last round I declined to discharge `proliferation_deflation` partly on scope. That reasoning has to
be applied consistently whichever way the numbers land, so it is written down now.

**Structural fact, verified today:** nothing in production builds an E0 claim dict. The project's
record of E0's admissibility exists **only** as a hardcoded fixture at
`tests/test_claim_guards.py:135` (`{"kind": "transfer", "proliferation_controlled": False, ...}`).
`validate_claim` reads real evidence from nowhere. So "discharging the blocker" is, mechanically,
editing that pinned dict — which is exactly why it must be deliberate.

**Pre-declared readings, both directions:**

- If the gap **SURVIVES**: this discharges the proliferation confound **for the D2 comparison** — the
  finding that perturbation-basis supervision underperforms Hallmark is not a proliferation artefact.
- Whether that also discharges `proliferation_deflation` for **E0's `transfer` claim** turns on a
  direction argument that is fixed here in advance: **D2's result is that PBS supervision is *worse*
  than Hallmark.** Showing a *negative* result about the perturbation basis is not proliferation-driven
  does not establish that E0's *positive* alignment is more than proliferation. Those are different
  claims with opposite signs.
- Therefore, **stated before the numbers: a surviving D2 gap alone will NOT flip
  `proliferation_controlled` for E0.** The most it can do — combined with D2.3's finding that 85–95
  of the 95 non-proliferation axes stay legible at ~90% of median — is discharge the blocker for
  *the D2 gap and the legibility of the perturbation basis*, which is what I will record.
- If the gap **DOES NOT survive**, that is a substantive negative about D2's headline and is reported
  as such, not softened.

The one thing that would change this: if the deflated gap survives **and** the non-proliferation PBS
axes are shown to carry the E0 alignment itself. That is an E0 re-run, not this analysis, and it is
not in scope.

## 6. D3 seeds 43 and 44

Identical command to the seed-42 run (`d3_purity_result_20260803T1330Z.md`), only the artifacts
change: `--n-permutations 2000 --n-components 16 --n-draws 40 --purity-source absolute
--require-rna-positive-control`, complete-case ABSOLUTE purity, before/after on the same patients.

Pre-declared bar, same as seed 42: **`wsi_biology` retains ≥ 80% of `excess_over_null_median`,
observed still above `null_p95`, `permutation_p` still at the 1/2001 floor** — in each seed
separately. Seed 42 gave 94.2% (d2_h) and 97.9% (d2_i). Three seeds are reported individually; **no
pooling and no averaging across seeds**, because training is not seed-reproducible on this stack and
only within-run paired differences are quotable. A seed that misses the bar is reported as a miss.

### Files / commits

- Axis cut sets: `~/e0_run/d3/d2_deflate/axes/{prol_top100,prol_wmean,union,intersection,placebo_random,placebo_random_union_sized}.txt`, `axis_sets.json`
- This pre-registration, committed before the deflated comparison runs.
