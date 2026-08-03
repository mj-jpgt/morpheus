# D2.3 — the falsifier does not fire: 117–127 of 128 PBS axes are legible from morphology and only ~25% of them are proliferation-loaded, exactly chance. But the single most proliferative axis is the single most legible one in 6/6 runs, the two pre-declared statistics disagree, and I am NOT discharging `proliferation_deflation`

**Logged:** 2026-08-03 13:45 UTC. Pre-registered in
`NOTEBOOK_ENTRIES/d3_d2p3_preregistration_20260803T1300Z.md`, committed `cd9b056`, **before** any
legibility number existed.

**How obtained:** Lambda box `150.136.45.194`, `~/ws_d3`, CPU only, thread-capped.
`python -m morpheus.v2.research.rebase.d2_axis_proliferation --artifacts d2_{h,i}_seed{42,43,44}.npz
--pbs-targets pbs_targets_k128_v2.npz --annotations gene_annotations.parquet --state wsi_biology
--partition test --n-permutations 200 --top-k 100 --seed 42 --n-jobs 6`.
Outputs `~/e0_run/d3/d2_3/{axis_table.csv,d2_3_report.json}`. n = 2,766 test patients, 105 confound
columns (cancer + pooled TSS), 6 artifacts × 128 axes.

---

### Technical

#### Lead finding: the pre-built annotation is a weak instrument, and this was recorded before the run

`build_pbs_targets.py:107` defines `proliferation_loading` as the **|loading|-weighted mean over all
7,072 basis genes** of a binary MSigDB proliferation flag. The SVD gene basis is dense — an axis's
top-100 genes carry a median of just **7.1%** of its |loading| mass — so the weights are near-uniform
and every axis is squeezed onto the background rate.

| statistic | min | median | max | spread ÷ background |
|---|---:|---:|---:|---:|
| background rate over the 7,072 genes | — | **0.0841** | — | — |
| `proliferation_loading` as built (`prol_wmean`) | 0.0738 | 0.0820 | 0.1453 | **0.85×** |
| proliferation fraction among top-100 |loading| genes | 0.021 | 0.101 | 0.378 | **4.24×** |

The median axis's shipped `proliferation_loading` sits at **0.98× background** — the column barely
varies. And the two statistics rank the axes differently: Spearman **0.577** [0.432, 0.697]. A unit
test (`test_top_k_statistic_separates_a_proliferation_axis_that_the_weighted_mean_dilutes`) pins the
mechanism: on a dense basis, an axis that leads entirely on proliferation genes and one that leads
entirely on non-proliferation genes get the *same* weighted mean.

So the ledger's "the analysis is free, the annotations already exist" is **half true**. The
annotations exist; one of them is close to uninformative. Both statistics are carried as co-primary
throughout, per the pre-registration.

#### Legibility, with a null that is not zero

Per-axis legibility = held-out cross-fitted ridge single-direction correlation from `wsi_biology` to
each PBS axis, both sides residualised on cancer+TSS. The null uses the **same within-cancer row
permutation** as `calibra.calibration.permutation_null` (200 permutations); a row-shuffle null was
not substituted. Measured null median over axes −0.022 to −0.007, null p95 **0.025–0.037**.

| artifact | legibility p25 / p50 / p75 / p100 | n legible | top-10 share | bimodal? |
|---|---|---:|---:|---|
| d2_h_seed42 | 0.098 / 0.136 / 0.194 / 0.438 | 120/128 | 0.165 | no (BIC k=1) |
| d2_h_seed43 | 0.107 / 0.155 / 0.199 / 0.445 | 127/128 | 0.153 | no |
| d2_h_seed44 | 0.088 / 0.123 / 0.171 / 0.403 | 118/128 | 0.176 | no |
| d2_i_seed42 | 0.077 / 0.124 / 0.184 / 0.400 | 117/128 | 0.178 | no |
| d2_i_seed43 | 0.083 / 0.132 / 0.181 / 0.417 | 120/128 | 0.169 | no |
| d2_i_seed44 | 0.072 / 0.120 / 0.171 / 0.392 | 117/128 | 0.182 | no |

**The distribution, as required, not a summary.** It is **unimodal** — a 3-way BIC comparison prefers
one Gaussian component in 6/6 artifacts — and it is **not dominated by a few axes**: the top 10 of
128 carry only 15–18% of total legibility (uniform would be 7.8%), the top 5 carry 8.5–10.2%. The
largest gap in the sorted values is 12–21% of the spread, with no clean separation. It is a smooth
gradient from ~0.02 to ~0.44, not two populations.

#### The falsifier

Ledger falsifier: *every legible axis coming back proliferation-loaded.* Proliferation-loaded =
top quartile, so chance is 0.25.

| artifact | statistic | legible & loaded | share | chance | Spearman(legibility, proliferation) 95% CI | verdict |
|---|---|---:|---:|---:|---|---|
| d2_h_s42 | `prol_wmean` | 32/120 | 0.267 | 0.25 | 0.307 [0.127, 0.471] | PARTIAL |
| d2_h_s42 | `prol_top100` | 30/120 | 0.250 | 0.25 | 0.237 [0.051, 0.413] | DISCHARGED |
| d2_h_s43 | `prol_wmean` | 32/127 | 0.252 | 0.25 | 0.338 [0.155, 0.500] | PARTIAL |
| d2_h_s43 | `prol_top100` | 32/127 | 0.252 | 0.25 | 0.261 [0.073, 0.430] | DISCHARGED |
| d2_h_s44 | `prol_wmean` | 31/118 | 0.263 | 0.25 | 0.290 [0.103, 0.453] | DISCHARGED |
| d2_h_s44 | `prol_top100` | 29/118 | 0.246 | 0.25 | 0.198 [0.005, 0.367] | DISCHARGED |
| d2_i_s42 | `prol_wmean` | 31/117 | 0.265 | 0.25 | 0.282 [0.090, 0.453] | DISCHARGED |
| d2_i_s42 | `prol_top100` | 32/117 | 0.274 | 0.25 | 0.228 [0.036, 0.400] | DISCHARGED |
| d2_i_s43 | `prol_wmean` | 31/120 | 0.258 | 0.25 | 0.303 [0.112, 0.470] | PARTIAL |
| d2_i_s43 | `prol_top100` | 31/120 | 0.258 | 0.25 | 0.257 [0.059, 0.428] | DISCHARGED |
| d2_i_s44 | `prol_wmean` | 31/117 | 0.265 | 0.25 | 0.291 [0.105, 0.456] | DISCHARGED |
| d2_i_s44 | `prol_top100` | 32/117 | 0.274 | 0.25 | 0.243 [0.050, 0.412] | DISCHARGED |

**The falsifier as stated does not fire, and it is not close.** The share of legible axes that are
proliferation-loaded is **0.246–0.274 against a chance rate of 0.25** — indistinguishable from
random — where the falsifier needed ≥0.90.

**But three things cut the other way and are reported here, not buried:**

1. **There is a real, reproducible proliferation gradient.** Spearman is positive and its 95% CI
   excludes zero in **12/12** tests, ρ = 0.198–0.338. More proliferative axes *are* somewhat more
   legible. Essentiality shows the same weaker pattern (ρ = 0.128–0.244, CI excludes zero in 8/12).
2. **The single most proliferation-loaded axis is the single most legible axis.** PBS_001 has the
   highest top-100 proliferation fraction of all 128 (0.378, i.e. 4.5× background) and ranks **1st in
   four artifacts and 2nd in the other two** for legibility. That is exactly the pattern the blocker
   describes, perfectly reproducible, and no amount of aggregate statistics makes it go away.
3. **The two pre-declared co-primary statistics disagree** — `prol_wmean` returns PARTIAL on 3 of 6
   artifacts, `prol_top100` returns DISCHARGED on 6 of 6. The pre-registration says a verdict is
   only clean if the two agree. They do not.

#### Proliferation deflation, done directly

Dropping the top-quartile proliferation axes and re-reading the rest:

| artifact | median legibility, all 128 | median, 96 non-proliferation axes | retained | n legible |
|---|---:|---:|---:|---|
| d2_h_seed42 | 0.1365 | 0.1253 | 91.8% | 90/95 |
| d2_h_seed43 | 0.1552 | 0.1418 | 91.4% | 95/95 |
| d2_h_seed44 | 0.1229 | 0.1175 | 95.6% | 89/95 |
| d2_i_seed42 | 0.1243 | 0.1111 | 89.4% | 85/95 |
| d2_i_seed43 | 0.1323 | 0.1213 | 91.7% | 89/95 |
| d2_i_seed44 | 0.1199 | 0.1098 | 91.6% | 85/95 |

**Delete every proliferation axis and 85–95 of the remaining 95 are still legible, at ~90% of the
median.** The morphology→axis channel is not carried by proliferation.

#### Post-hoc control, labelled as post-hoc

Not pre-declared. The proliferation gradient is confounded with **axis strength**: the PBS singular
value correlates with legibility (ρ 0.398–0.442) *and* with proliferation loading (ρ 0.453) — early
high-variance axes are both more prominent and more proliferative. Partial Spearman controlling for
the singular value collapses the gradient to **0.014–0.087** in all 6 artifacts on both statistics,
i.e. essentially nothing. So the gradient looks like axis prominence rather than proliferation per se.

**This is deliberately not used to upgrade the verdict.** It was found after seeing the result, and
converting a PARTIAL into a DISCHARGE with a post-hoc control is precisely the fitted-then-explained
move this project's process exists to prevent. It is recorded so a later pre-registered analysis can
test it properly.

### In plain terms

The fear was that when a model reads a slide and predicts molecular programmes, the only programme it
is really reading is "how fast are these cells dividing" — the most generic signal in cancer, which
would look identical to a real finding and mean far less.

That is not what we see. Almost every one of the 128 molecular axes is readable from the image, and
the readable ones are no more likely to be proliferation axes than chance. Throw away every
proliferation axis and nearly all the rest are still readable at nearly full strength.

Two honest dents. First, the axis that *is* most about proliferation is the one the image reads best
— every single time, in all six runs. Second, the proliferation annotation the ledger told us to use
turned out to be nearly flat: it averages over all 7,000 genes with almost equal weight, so it cannot
really tell a proliferation axis from any other. We had to build a sharper one to ask the question at
all, and the two versions do not fully agree.

### Meaning for the claim

**I am NOT discharging `claim_guards.proliferation_deflation`, and I did not touch
`tests/test_claim_guards.py`.** The suite is green as it stands. Three reasons, in increasing order
of importance:

1. My own pre-declared rule required the two co-primary statistics to agree. They do not (3/6).
2. The control that would rescue a clean discharge (partial correlation on axis strength) is post-hoc.
3. **Decisive: D2.3 does not perform the remedy the blocker names.** `claim_guards` scopes
   `proliferation_deflation` to the **`transfer`** claim, and its stated remedy is *"re-run with
   proliferation/cell-cycle programme regressed out, or with the responsive arm stratified by
   proliferation loading, and show **the gap** survives"* — the gap being the E0/D2 alignment gap
   between cell-line perturbation structure and tumour expression. What D2.3 measures is whether
   morphology-legible PBS axes are proliferation axes. That is real, relevant evidence about the
   mechanism, and it comes out favourably — but it is **not** the same quantity, and discharging a
   blocker with an adjacent measurement is how a caveat gets quietly lost.

**What would actually discharge it:** re-run `d2_compare`'s Hallmark-vs-PBS readout with the
proliferation-loaded component of the targets removed (or restricted to the non-proliferation axes)
across all three seeds, and show the −0.109 to −0.133 gap survives. That is CPU-only and cheap; it
was not in this task's scope.

What D2.3 *does* license: the statement that the morphology→PBS-axis channel is broad rather than
proliferation-confined, with the PBS_001 caveat attached. It also supplies the first per-axis
legibility distribution on this project, which **P4**'s per-axis certification will need.

### Files / commits

- `~/e0_run/d3/d2_3/{axis_table.csv,d2_3_report.json}`, `~/e0_run/d3/axis_annotation_sharpness.csv`
- `~/e0_run/d3/run_d23.sh`, `~/e0_run/d3/logs/d23.log`
- Code: `9bc7085` — `v2/research/rebase/d2_axis_proliferation.py` + `v2/tests/test_d2_axis_proliferation.py`
