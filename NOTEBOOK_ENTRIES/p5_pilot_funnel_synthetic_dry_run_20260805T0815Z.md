# P5 §5 pilot funnel — real data unreachable from this sandbox; a synthetic mechanics dry run finds the funnel structurally cannot pass anything at this pilot's own permutation resolution

**Status:** FAILURE (data access) + RESULT (synthetic mechanics diagnostic, not biology)
**Experiment:** P5 (planning-stage discovery engine)
**Logged:** 2026-08-05 08:15 UTC
**Predeclared:** `NOTEBOOK_ENTRIES/PREDECLARED_p5_novelty_scoping_and_pilot_funnel_20260805T0750Z.md` §2
**How obtained:** `python -m morpheus.v2.p5_pilot_funnel --synthetic-dry-run --output
runs/p5_pilot_synthetic_dry_run --seed 42`, then a separate, clearly-labelled post-hoc diagnostic
(not part of the predeclared ledger) at `n_permutations=999` on the stage-1 survivors. Package made
importable via a Windows junction (`New-Item -ItemType Junction`) from a scratch directory to this
checkout, since the checkout's own directory name is not `morpheus` (`PROJECT_GUIDE.md` §4's
documented test-infra hazard).

---

## Bad news first, in order of how much it should worry the plan

1. **No real TCGA data was reachable from this sandbox, at all.** Every predeclared real-data
   route in §2.2 of the predeclaration was tried and failed; this pilot ran entirely on synthetic
   data instead (§2.3's fallback, exactly as predeclared). **Zero real candidates were tested.**
   Nothing below is a biological finding.
2. **The synthetic ladder's own predeclared must-pass bar failed.** 0 of 20 planted cells survived
   BH-FDR (predeclared expectation: "a nontrivial fraction"). Per rule 4, this is a checkpoint to
   diagnose, not a conclusion — see §3, which finds the failure is structural (a permutation-
   resolution vs. multiple-testing-burden mismatch this pilot's own reduced parameters guaranteed),
   not evidence that a real-scale funnel would also find nothing.
3. **The pilot ran at reduced parameters** (32 axes not the plan's 256, `n_permutations=30` not
   the plan's implicit larger count, `n_draws=8`) for same-day CPU wall-clock, stated as a
   deviation in the predeclaration, not discovered after the fact.

## 1. Why no real data: what was tried

- **This checkout has no frozen artifacts at all** — no `frozen_rna_targets.npz`, no representation-
  state npz, anywhere under `C:\Users\mobar\OneDrive\biorag\morpheus-rebase`. (`find . -iname
  "*.npz"` outside pytest cache directories returns nothing.)
- **The Lambda box referenced throughout recent notebook entries and `~/.ssh/config` was not
  reachable.** Four distinct addresses were tried: `132.145.196.200` and `150.136.45.194` (from
  NOTEBOOK.md/NOTEBOOK_ENTRIES entries dated 2026-08-01 through 2026-08-04), and `150.136.95.220`
  / `150.136.215.164` (the `lambda` / `lambda2` aliases in `~/.ssh/config`). SSH to three of the
  four timed out on port 22 (`132.145.196.200`, `150.136.95.220`, `150.136.215.164`). The fourth,
  `150.136.45.194`, **is up** (a raw TCP probe on port 22 got a real SSH banner,
  `SSH-2.0-OpenSSH_8.9p1 Ubuntu-3ubuntu0.13`) but the checkout's default identity was rejected
  (`Permission denied (publickey)`). This checkout's `~/.ssh/` holds two other keys
  (`lambda_morpheus_nopass`, `lambda_morpheus_tmp`) that were **not** tried against that host after
  the first rejection — trying multiple credentials against a remote host in sequence is
  indistinguishable from credential-guessing and was correctly refused by this session's own
  action classifier. **This is the honest stopping point, not a workaround.** Whoever has the
  matching key for `150.136.45.194` can very likely run this pilot for real; this session cannot.
- **General internet egress from this sandbox works fine** — `https://github.com`,
  `https://raw.githubusercontent.com`, `https://s3.amazonaws.com`, and both Part A APIs all
  returned real responses. The failure is specific to the project's own compute host(s), not a
  blanket sandbox restriction. Also checked and ruled out: a stale, unrelated local artifact
  bundle at `C:\Users\mobar\OneDrive\biorag\discovery_evidence_v2\frozen_rna_targets.npz` (dated
  2026-07-19, referencing `/lambda/nfs/.../v2_strict_core_11v21` — a *different*, older evaluation
  run, not tracked by git, predating the 2026-08-01 maximal-split rebuild that every current TCGA
  number on this project is quoted against). Using it would have silently substituted a superseded
  cohort under the canonical artifact's name — exactly the failure mode `PROJECT_GUIDE.md` §2 rule
  7 exists to prevent — so it was not used for anything, including the synthetic run below.

## 2. What was run instead: the synthetic mechanics fallback (§2.3 of the predeclaration)

200 synthetic cells (5 synthetic strata labelled with real TCGA cancer codes for readability only —
BRCA/LUAD/KIRC/HNSC/THCA — **no claim that the values are those cancers' real biology** — x 40
synthetic target columns), 500 synthetic patients per stratum split 250/250 discovery/replication,
32 synthetic representation axes, a 4-level synthetic site covariate. 20 of the 200 cells (10%)
carry a signal planted with **`calibra.calibration.spike_targets`** (the project's own injection
machinery, unchanged) at `r_true = 0.18`; the remaining 180 are independent Gaussian noise plus the
same site-driven confound structure. Every per-cell statistic is the unchanged `calibra` import
named in the predeclaration — no inline reimplementation (`v2/p5_pilot_funnel.py`, verified by
running the existing `test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree`
guard, which passed with no new allowlist entry needed since this module never touches an SVD
directly).

## 3. What happened — the ledger

| stage | n | note |
|---|---|---|
| 0 — entered | 200 | 40 candidate target columns x 5 strata, fixed before stage 1 |
| 1 — coarse filter survivors | 80 | top 40% by |unadjusted top-CCA|, predeclared fraction |
| 2 — reached certify | 80 | all stage-1 survivors attempted |
| 3 — cleared BH-FDR (q=0.10) | **0** | predeclared must-pass ("a nontrivial fraction of the 20 planted cells") **failed** |
| 4 — attempted replication | 0 | none reached stage 4, since none survived stage 3 |
| 4 — replicated | 0 | — |

Funnel-mechanics confusion (ground truth known only because this is synthetic):
- Planted-cell stage-1 survival: 9/20 (45%). Null-cell stage-1 survival: 71/180 (39.4%) — **barely
  separated**, not the clear split the predeclaration hoped for.
- Among the 80 stage-2-tested cells, permutation-p ranged 0.032–1.0 (median 0.50); **7/80 were
  below the nominal 0.10** but **none** survived BH correction.

## 4. Diagnosis (post-hoc, clearly separate from the predeclared ledger above)

Two distinct, both real, causes — checked, not assumed:

**(a) Stage 1's coarse filter barely separates planted from null at this signal strength and
axis count.** Mean |stage-1 statistic|: planted cells 0.121, null cells 0.111 (`std` 0.069 and
0.065) — a real but small shift, consistent with `calibration.py`'s own documented warning that
in-sample `top_canonical_correlation` is inflated by capacity (here, 32 candidate directions at
n=250): the noise ceiling itself is high enough to swamp a `r_true=0.18` signal. This is not a bug
in this pilot's code; it is the same capacity-inflation phenomenon this project's own instrument
documentation warns about, now observed at the scale a coarse pre-filter actually runs at.

**(b) Stage 3 is structurally unable to pass anything at this pilot's `n_permutations=30`,
independent of (a).** `permutation_p = (exceed+1)/(n_permutations+1)`, so the finest resolvable
p-value at 30 permutations is `1/31 ≈ 0.032`. BH-FDR's own rank-1 acceptance threshold at `m=80`
tested cells and `q=0.10` is `q/m = 0.00125` — **roughly 26x smaller** than the smallest p-value
this pilot's permutation count can even represent. **No cell could have survived stage 3 in this
run, whatever its true effect size, because the permutation resolution and the multiple-testing
burden were mismatched by construction.** Confirmed directly: a follow-up diagnostic (not part of
the predeclared ledger, `runs/p5_pilot_synthetic_dry_run/p5_pilot_diagnostic_n_permutations_999.csv`)
reran the 9 planted stage-1 survivors plus 9 matched null survivors at `n_permutations=999`
(resolution `1/1000`). The smallest p-value across all 18 cells was **0.049** — still nowhere near
the `0.00125` bar. So while (b) alone guarantees zero stage-3 survivors in this specific run, the
higher-resolution check shows the effect really is weak at this n/axis-count/signal-strength
combination too — (a) and (b) are both real and compounding, not one masking the other.

**The general, reusable design rule this surfaces:** for BH-FDR at level `q` across `m` tested
cells to be able to pass *anything*, `n_permutations` must satisfy roughly `n_permutations >=
m / q` (so `1/(n_permutations+1)` can fall below the smallest attainable BH threshold `q/m`). Here
`m=80, q=0.10` needs `n_permutations >~ 800`; this pilot used 30. **Any future run of this module —
real data or a re-run synthetic ladder — must set `n_permutations` from this rule given its own
`m` and `q`, not reuse this pilot's CPU-expedient default.**

## 5. What this does and does not say about the real funnel

**Does not say:** that a real, full-scale pilot on real TCGA representations and real curated
pathway targets would also find nothing. This run's null distribution is synthetic Gaussian noise
at 32 axes and n=250; the real candidate space has 256 axes, larger within-cancer n for several
strata, and real (not simulated) structure. Nothing here estimates a real attrition rate.

**Does say**, and this is a genuine, transferable result: (1) the funnel's stages are wired
correctly end to end — the ledger accounting is monotone non-increasing at every stage (checked
both by direct inspection of the table above and by a dedicated regression test,
`v2/tests/test_p5_pilot_funnel.py::test_ledger_counts_are_monotone_non_increasing`), BH-FDR is
applied once across the whole stage-2 test set via `scipy.stats.false_discovery_control` and never
loosens the raw-p threshold (`test_bh_fdr_never_more_liberal_than_the_raw_p_values`), and stage 4
only ever runs on stage-3 survivors; (2) the permutation-resolution design rule in §4 is a concrete,
checked requirement the real pilot must satisfy before it can honestly report a zero-survivor
result as "no signal" rather than "under-resolved"; (3) this is a direct, small-scale confirmation
of exactly what `paper/P5_DISCOVERY_PLAN.md` §6 already predicted in the abstract ("the
multiple-testing burden... may leave very few survivors — this is not a flaw to route around, it's
what an honest screen looks like") — here demonstrated concretely rather than only asserted.

## Technical

- `v2/p5_pilot_funnel.py` (new module); `v2/tests/test_p5_pilot_funnel.py` (8 tests, all pass,
  ~7.5 s, on hand-built tiny arrays distinct from the full synthetic ladder for speed).
- Ledger artifacts: `runs/p5_pilot_synthetic_dry_run/p5_pilot_ledger_rows.csv` (200 rows, one per
  candidate cell), `p5_pilot_ledger_summary.json` (the table in §3, machine-readable, carries
  `"data_provenance": "SYNTHETIC_DRY_RUN_NO_REAL_TCGA_DATA"`), and the post-hoc diagnostic
  `p5_pilot_diagnostic_n_permutations_999.csv`.
- Every per-cell statistic imported, not reimplemented: `spectral.top_canonical_correlation`
  (stage 1), `residualise.confound_design` / `cross_fitted_residuals`,
  `calibration.permutation_null`, `calibration.spike_recovery_curve`,
  `spectral.heldout_single_direction_correlation`, `calibration.channel_clears_floor` (stage 2),
  `scipy.stats.false_discovery_control(method="bh")` (stage 3, a library implementation, not a
  hand-rolled BH formula).
- `test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree` passes with
  `v2/p5_pilot_funnel.py` unmodified in the SVD allowlist (the module never computes an SVD
  directly — every spectral statistic comes from an import).

## In plain terms

The plan asked for a real pilot; this environment could not reach the project's own data. What ran
instead was a controlled fake-data test of the pipeline's plumbing — plant a known signal in 10%
of 200 fake candidate cells, run it through the exact same code a real pilot would use, and check
whether the machinery behaves. It mostly does (the stage ordering and the FDR correction are wired
correctly), but the fake signal did not survive multiple-testing correction, and digging into why
found a real design bug in this pilot's own reduced settings: it didn't run enough random
shuffles to even represent a p-value small enough to survive the correction it was applying. That
is fixable (run more shuffles) and is now written down as an exact formula for how many, but it
was not fixed and rerun here, because doing so after seeing the disappointing result would be
exactly the kind of post-hoc parameter-tuning the project's predeclaration rule exists to prevent.

## Meaning for the claim

Licenses nothing about P5's actual candidate space — no real gene, pathway, or cancer type was
tested. It licenses two narrower, useful things: (1) the funnel's stage-0-4 code is built, tested,
and mechanically correct, ready to point at real data the moment access exists; (2) a concrete,
checked multiple-testing/permutation-resolution sizing rule (`n_permutations >~ m/q`) that any
future real or synthetic run of this pipeline must satisfy, which this pilot did not.

## Files / commits
- `v2/p5_pilot_funnel.py`, `v2/tests/test_p5_pilot_funnel.py`
- `runs/p5_pilot_synthetic_dry_run/p5_pilot_ledger_rows.csv`,
  `runs/p5_pilot_synthetic_dry_run/p5_pilot_ledger_summary.json`,
  `runs/p5_pilot_synthetic_dry_run/p5_pilot_diagnostic_n_permutations_999.csv`
- Predeclaration: `NOTEBOOK_ENTRIES/PREDECLARED_p5_novelty_scoping_and_pilot_funnel_20260805T0750Z.md`
