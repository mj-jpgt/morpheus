# Leave-sites-out: the channel survives unseen sites adjusted, and does not clear the bar unadjusted

**UTC** 2026-08-04T18:30Z
**Predeclaration** `NOTEBOOK_ENTRIES/PREDECLARED_leave_sites_out_20260804T1745Z.md` (commit `502f469`)
**Code** `v2/calibra/leave_sites_out.py`, `v2/calibra/spectral.py` (`heldout_top_cca_indexed`)
**Run** `/lambda/nfs/geeg/biorag3_persistent_20260711/runs/leave_sites_out_20260804/`
**Artifact** `.../morpheus_phase_d/runs/d2_final/artifacts/d2_h_seed42.npz`, `test` partition, n=2,766
**Targets** `frozen_rna_targets.npz`, 90 non-control targets
**Workspace** `~/ws_lso/morpheus`, 611/611 files verified against `git ls-tree -r HEAD`, 0 mismatches

The reading was fixed before the run and is applied as written. Nothing below reinterprets it.

---

## Bad news first

**The strong version of the claim is not established.** The predeclared bar for the
**unadjusted** arm returns **collapses** — 2 of 5 folds, and the rule needed 3.

**The `no_external_cohort` blocker is untouched and `claim_guards.py` is not edited.**
This is one institution network, one set of protocols, one era. It is evidence about
acquisition-condition generalisation, not an external cohort.

## What survives

**Adjusted arm, `wsi_biology` — SURVIVES, 5 of 5 folds, at full strength.**

The number that matters is not the channel itself but its ratio to the matched random
split: **1.010**. Holding out whole tissue source sites costs **nothing** measurable
relative to a random split of identical size and identical per-cancer composition.

## The two arms

`wsi_biology`, n_components=16, 1,000 permutations (resolution 1/1001), 1,000 bootstrap
resamples. `obs` = top canonical correlation with directions fit on the training sites and
scored on the held-out sites. `CI95` is the **site-cluster** bootstrap. `random` is the
matched random split.

### Adjusted (`cancer + pooled TSS`, cross-fitted, within each block)

| fold | held-out n | sites | obs | null p95 | p | site CI95 | matched random |
|---:|---:|---:|---:|---:|---:|---|---:|
| 0 | 864 | 32 | 0.4319 | 0.0927 | 0.0010 | [0.317, 0.635] | 0.5928 |
| 1 | 501 | 65 | 0.5797 | 0.0934 | 0.0010 | [0.506, 0.650] | 0.5860 |
| 2 | 482 | 82 | 0.6093 | 0.0994 | 0.0010 | [0.547, 0.664] | 0.5639 |
| 3 | 466 | 84 | 0.5423 | 0.0977 | 0.0010 | [0.453, 0.623] | 0.5428 |
| 4 | 453 | 89 | 0.5796 | 0.1012 | 0.0010 | [0.485, 0.657] | 0.5737 |

median **0.5796** vs matched random **0.5737**, ratio **1.010**, **5/5 survive**.
Every fold clears its null at the resolution floor and every site-cluster interval sits
entirely above the null p95.

### Unadjusted

| fold | obs | null p95 | p | site CI95 | clears bar | matched random |
|---:|---:|---:|---:|---|:--:|---:|
| 0 | 0.7676 | 0.7058 | 0.0010 | [0.608, 0.857] | no | 0.7996 |
| 1 | 0.7745 | 0.6623 | 0.0010 | [0.651, 0.848] | no | 0.7929 |
| 2 | 0.8108 | 0.7072 | 0.0010 | [0.707, 0.865] | no | 0.8135 |
| 3 | 0.8035 | 0.6861 | 0.0010 | [0.687, 0.868] | yes | 0.8002 |
| 4 | 0.7302 | 0.6067 | 0.0010 | [0.618, 0.808] | yes | 0.7595 |

median **0.7745** vs matched random **0.7996**, ratio **0.969**, **2/5 survive** -> collapses.

**Why it fails, stated precisely, because the point estimate is not what failed.** The
observed value clears the permutation null in **5 of 5** folds, at p = 0.0010 every time.
What fails is the second half of the predeclared bar: the site-cluster interval's lower
bound must also exceed the null p95, and in folds 0-2 it does not.

The unadjusted null p95 is **0.61-0.71** — enormous next to the adjusted arm's 0.09-0.10.
That is correct behaviour, not a defect. The null permutes the pairing *within cancer*, so
every patient keeps their cancer label; an unadjusted representation encodes cancer type
strongly, both projections therefore stay correlated under permutation, and the null
absorbs it. So most of the unadjusted 0.77 **is cancer identity**. The excess over the null
is small, and with only 32-89 independent sites in a fold the site-level uncertainty on
that excess is too wide to exclude a site-level fluctuation.

The honest statement: unadjusted, we cannot separate the residual excess from site-level
noise. Not "the channel is absent unadjusted" — the p-values say otherwise — but the
predeclared bar is the bar, and it was not cleared.

## Fold 0 is where site shift actually bites

Fold 0 is the only fold with a real gap: **0.4319 held-out sites vs 0.5928 matched random**,
a 27% shortfall. Folds 1-4 show none (ratios 0.99, 1.08, 1.00, 1.01).

This is a property of the fold builder, and it is informative rather than a flaw. The greedy
within-cancer assignment hands the largest remaining site to the emptiest fold, so fold 0
collects **the single biggest site of every cancer**: 32 sites carrying 864 patients, against
89 sites carrying 453 in fold 4. Holding out the largest sites is the hardest version of the
test, and it is the one version where the penalty is visible. Any future site-robustness
claim should be graded on a fold built this way, not on an average over easy folds.

## Fold composition

TSS is heavily imbalanced and the folds inherit it. All 352 raw sites, 2,766 patients:

| fold | held-out patients | held-out sites | training patients | training sites | held-out cancers | shared with training | largest held-out site |
|---:|---:|---:|---:|---:|---:|---:|---|
| 0 | 864 | 32 | 1,902 | 320 | 21 | 21 | 117 (13.5%) |
| 1 | 501 | 65 | 2,265 | 287 | 21 | 21 | 55 (11.0%) |
| 2 | 482 | 82 | 2,284 | 270 | 21 | 21 | 42 (8.7%) |
| 3 | 466 | 84 | 2,300 | 268 | 20 | 20 | 25 (5.4%) |
| 4 | 453 | 89 | 2,313 | 263 | 19 | 19 | 21 (4.6%) |

85 pooled site classes at `min_site_count=10` (84 frequent + OTHER) — the same 85 that give
the confound certificate its 1/85 = 0.0118 chance rate, so the two instruments agree about
what a site is. Folds 3 and 4 carry 20 and 19 cancers rather than 21 because ACC has only 3
sites and PCPG 5; a cancer cannot populate more folds than it has sites.

## Controls

- **Positive control** `rna_biology` (circular by construction): survives 5/5, median 0.8811
  adjusted. The harness reacts when the pairing is real.
- `full_biology`: survives 5/5, 0.8827 adjusted, ratio 0.996.
- Both are near-circular and are reported only as evidence the instrument works.

## What this does and does not buy

**Does:** the morphology->molecular channel is not an artifact of the specific sites it was
measured on. Fit the read-out on one set of sites, score it on sites never seen, and after
CALIBRA's own adjustment it returns essentially the same number as a random split of the
same size — 5/5 folds, ratio 1.010. That is the main thing an external cohort is for, and it
now has a measurement rather than a promise.

**Does not:**
- **Discharge `no_external_cohort`.** Same institution network, same protocols, same era,
  same extraction pipeline. The blocker stands and `claim_guards.py` is unedited.
- Establish the unadjusted claim. See above.
- Say anything about a *different* cohort. This project's own measurement is that a
  classifier separates TCGA from HEST at AUC 0.99999
  (`spatial_baselines_20260803T0620Z.md`). TCGA site-to-site shift is a far smaller thing
  than TCGA-to-elsewhere shift, and surviving the former does not predict surviving the latter.
- Test the encoder. The encoder is frozen; this is a **read-out** transfer test. It is kept
  honest by the fact that the encoder saw none of these patients, sites or cancers (the
  `test` partition holds out cancers, and site is nested in cancer), but a retrained-per-fold
  version would be strictly stronger and was not run.

## Naming

Two different generalisation tests, never merged:

- **held-out-cancer** — the existing `tumor_state_heldout_cancer` split.
- **held-out-site** — this entry.

## A limitation of my own criterion

The site-cluster interval is systematically wider on the site-held-out folds than on the
matched random folds, because the matched folds draw patients from many partial sites and so
have many more independent clusters. The CI half of the predeclared bar is therefore
stricter on the site arm than on the comparator. That asymmetry is *appropriate* — site-level
uncertainty is real and the comparator genuinely does not carry it — but it means the
"site arm fails CI, matched arm passes CI" contrast in the unadjusted table is partly
mechanical and should not be read as a second, independent piece of evidence.

## Reproduce

```
PYTHONPATH=$WS python -m morpheus.v2.calibra.leave_sites_out \
  --artifacts .../d2_final/artifacts/d2_h_seed42.npz \
  --targets   .../frozen_rna_targets.npz \
  --output    .../runs/leave_sites_out_20260804 \
  --partition test --n-folds 5 --n-components 16 \
  --n-permutations 1000 --n-boot 1000 --seed 42
```
Folds are seed-free by construction; `--seed` moves only the permutation, bootstrap and
matched-comparator draws.
