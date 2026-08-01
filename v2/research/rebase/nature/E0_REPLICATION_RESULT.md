# E0 cross-lineage replication — **it replicates, and the n-matching was load-bearing**

**Runs:** `runs/e0_replication_20260731/` — 4 runs, ~400 s each on A100, commit `94f0e53`, clean worktree.
**Settings:** k ∈ {10,25}, q=60, 1,000 Haar draws, 1,000 bootstrap resamples, both arms capped at n=168.
**Preregistered:** `nonresponsive_p = 0.2` — the threshold `E0_RESULT.md` §4 already named, before any of
these numbers existed. The other three runs are a **declared sensitivity analysis**; the verdict comes
from the preregistered run alone.

---

## 1. Result

All four runs: **`verdict = supported`, `gates_pass = True`, 4 of 4 contexts decided, none undecidable.**

Preregistered run (`nonresponsive_p = 0.2`, both arms n = 168), `signed_log1p`:

| lineage | k | responsive | control | gap | paired 95% CI |
|---|---:|---:|---:|---:|---|
| K562 | 10 | 0.0569 | 0.0073 | **+0.0496** | [+0.0302, +0.0463] |
| K562 | 25 | 0.0528 | 0.0141 | **+0.0387** | [+0.0245, +0.0346] |
| RPE1 | 10 | 0.0592 | 0.0205 | **+0.0387** | [+0.0220, +0.0611] |
| RPE1 | 25 | 0.0592 | 0.0198 | **+0.0394** | [+0.0298, +0.0500] |

**At k=25 the two lineages agree to within 2%** (+0.0387 vs +0.0394) — a near-triploid CML line and a
near-diploid retinal epithelial line, n-matched, producing the same effect.

Stability across the threshold sweep (k=25 gap, signed):

| `nonresponsive_p` | K562 pool | RPE1 pool | K562 gap | RPE1 gap |
|---|---:|---:|---:|---:|
| 0.05 | 4,251 | 285 | +0.0353 | +0.0361 |
| 0.10 | 3,500 | 227 | +0.0369 | +0.0397 |
| **0.20** (preregistered) | 2,501 | 168 | **+0.0387** | **+0.0394** |

K562's gap declines **monotonically** as the control is loosened (0.0387 → 0.0369 → 0.0353), exactly the
conservative direction predicted: admitting weak responders into the control raises its alignment and
shrinks the gap. The conclusion does not depend on the threshold.

## 2. The n-matching was not bookkeeping — it was the result

The uncapped run is the control on our own method. Same threshold, same everything, arms at their
natural sizes:

| | K562 arms n | K562 gap (k=10) | RPE1 arms n | RPE1 gap (k=10) |
|---|---:|---:|---:|---:|
| **n-matched** | 168 | +0.0496 | 168 | +0.0387 |
| **uncapped** | 2,501 | **+0.0671** | 168 | +0.0387 |

**Capping changes K562's gap by 35%** (0.0671 → 0.0496) with RPE1 untouched. Uncapped, we would have
reported K562 at +0.0671 against RPE1 at +0.0387 and concluded the effect is far stronger in K562 —
**a conclusion produced entirely by sample size.** The pre-run simulation predicted 35–40% attenuation;
the measured value is 35%.

This is why the fix mattered: without it, the run would have returned a confident, wrong answer about
lineage specificity.

## 3. What this licenses

**Licensed:**

> In two Replogle Perturb-seq lineages — K562 (near-triploid CML) and RPE1 (near-diploid retinal
> epithelium) — perturbations with a detectable transcriptional effect align with TCGA tumour expression
> geometry more than **n-matched** perturbations with no detectable effect, at k=10 and k=25, in both RNA
> transforms and at all three control thresholds tested. **The E0 alignment is not specific to K562.**

**NOT licensed:**

> The alignment reflects perturbation biology rather than a property of the Perturb-seq assay, or of
> perturbation strength itself.

**There is no control in this design that separates biology from a shared-platform or
effect-strength artifact, and cross-lineage replication cannot supply one.** Both lineages share one
platform, one normalisation pipeline, one pseudobulk procedure, and one `energy_test_p_value` statistic
that is *monotone in effect size*. Effect size tracks gene essentiality; essential-gene knockdown
produces a stereotyped proliferation/stress programme; tumours are proliferative. **That confound
replicates across lineages precisely because it is biology of essential genes rather than of lineage** —
so this result is exactly what it would produce. `pc1_share` does not rescue it: proliferation is not a
single principal component after per-column z-scoring of two different matrices.

See `claim_guards.proliferation_deflation` and `claim_guards.single_platform`. **E0 remains an
inadmissible transfer claim**, pinned by `test_current_e0_result_is_not_yet_an_admissible_transfer_claim`.

## 4. What must NOT be compared

- **Gap magnitudes against the q=150 primary run** — different n and different estimator.
- **`normalised_alignment` (% of ceiling) across runs** — numerator and denominator both move with n.
- **k=50 and k=100** — not evaluated here, because RPE1's control arm cannot support them. Note that in
  the primary run the control absorbed **55%** of the alignment at k=100, so **the dropped k values are
  the least favourable ones.** Stating that explicitly; the omission would otherwise read as selection.

## 5. Highest-value next control — a wrong-target arm

Nothing in this design tests whether the alignment is to **tumour** expression specifically, as opposed
to any large human bulk-RNA panel. The cheapest decisive addition is **GTEx normal tissue** as the target,
run through the identical path: if the alignment is tumour-specific, gap(TCGA) > gap(GTEx). If the gaps
match, "perturbations align with tumour biology" collapses to "perturbations align with human
transcriptome structure," and the claim is far weaker than it reads.

Then, in cost order: an **essentiality-stratified responsive arm** (DepMap common-essential vs
non-essential, matched on effect magnitude) for the proliferation confound; and a **second platform**
(LINCS L1000, or non-Replogle Perturb-seq) for the shared-assay confound.

*Requires data we do not hold: GTEx bulk RNA, DepMap essentiality. Both are open-access.*
