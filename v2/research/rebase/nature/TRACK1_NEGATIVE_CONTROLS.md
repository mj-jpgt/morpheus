# Track 1 — the negative-control battery (T1.1–T1.8)

Phase-gate item **B** for P1 (CALIBRA, the instrument/methods paper).

Protocol common to everything below unless stated: artifacts
`~/e0_run/../runs/d2_final/artifacts/{d2_h_seed42,d2_i_seed42}.npz`; molecular block
`~/e0_run/data/frozen_rna_targets.npz`; partition `test` = **2,766 held-out patients**; confound
design **108 columns** (cancer + 84 pooled tissue-source sites, `min_site_count=10`); seed 42;
cross-fitted residualisation, 5 folds, Ridge α = 1.0. Machine: Lambda A100 box, CPU only, `~/ws_p1`.

**Headline: two controls came back the wrong way, and both are reported as defects.**

| # | control | direction required | verdict |
|---|---|---|---|
| T1.3 | site/scanner prediction, **raw** states | must FAIL | ❌ **DID NOT FAIL — defect** |
| T1.3 | site/scanner prediction, **adjusted** states | must FAIL | ✅ fails as required |
| T1.4 | random gene sets, floor-scale statistic | must FAIL | ✅ fails as required |
| T1.4 | random gene sets, fitted-direction statistic | *(observation)* | ⚠️ **reach 76–82% of real gene sets** |
| T1.5(ii) | shuffled gene labels, attribution | must collapse | ✅ collapses (by construction) |
| T1.5(i) | shuffled gene labels, subspace | must persist | ✅ persists — **and the pass is damaging** |
| T1.6 | modality-shuffled pairing | must FAIL | ✅ fails as required, at 1/2001 resolution |
| T1.7(a) | RNA→RNA circular control | must PASS | ✅ passes |
| T1.7(b) | known covariate at published strength | must PASS | ✅ passes (ER 0.878 vs pre-registered [0.78, 0.92]) |
| T1.7(c) | synthetic spike above the floor | must PASS | ✅ passes (transmission floor 0.01) |

---

## T1.3 — must-FAIL 1: site / scanner / batch prediction

Pass criterion, fixed in advance: per-axis out-of-fold balanced accuracy for pooled TSS must not
exceed the 95th percentile of a ≥1,000-draw **within-cancer** label-permutation null, and the axis
bootstrap CI must include the chance rate. Chance = 1/85 = **0.0118**. Permutation resolution
1/1001 = 0.000999.

| artifact | state | joint LDA **raw** | joint LDA **adjusted** | joint null p95 (adj) | per-axis max raw → adj | breaching axes raw → adj |
|---|---|---:|---:|---:|---|---|
| d2_h | **wsi_biology** | **0.3633** | 0.0118 | 0.0528 | 0.0532 → 0.0123 | 17 → 0 |
| d2_h | full_biology | 0.2630 | 0.0101 | 0.0668 | 0.0548 → 0.0107 | 60 → 0 |
| d2_h | rna_biology | 0.2563 | 0.0074 | 0.0654 | 0.0506 → 0.0139 | 58 → 0 |
| d2_i | **wsi_biology** | **0.2348** | 0.0052 | 0.0418 | 0.0511 → 0.0104 | 43 → 0 |
| d2_i | full_biology | 0.2689 | 0.0085 | 0.0758 | 0.0551 → 0.0102 | 61 → 0 |
| d2_i | rna_biology | 0.2744 | 0.0079 | 0.0732 | 0.0495 → 0.0106 | 48 → 0 |

Every raw joint permutation p is at the 1/1001 floor: **not one of 1,000 within-cancer label
permutations reached the observed joint accuracy in any state.**

Two things must be said about the raw arm and neither is optional.

1. **The leak is smeared, not concentrated.** The best single axis anywhere reaches 0.055 — 4.6×
   chance — and the median axis sits *below* its own permutation null p95. The joint discriminant over
   all 256 axes reaches 20–31× chance. **A per-axis-only certificate, which is literally what the T1.3
   spec asks for, would have passed this.** The joint test must be a required field of the certificate
   schema, not an optional extra. That is a finding about the certification rule itself.
2. **The adjustment fully discharges it.** After the same cancer + pooled-TSS cross-fitted
   residualisation CALIBRA applies before every measurement, joint accuracy falls 21–45× to at or below
   the chance rate, with zero breaching axes in all six state/artifact combinations. So the defect is a
   property of the **raw representation**, and no adjusted number on this project is reading site.

Consequences: P4's condition 4 blocks exposure of any raw axis; the certificate schema needs a
`certified_on = {raw | adjusted}` field, because certifying on the adjusted state and then showing the
raw axis is exactly the laundering P4 forbids. `e0_basis_transfer.py:923`'s
`G3.5 = unavailable_no_site_labels` can now be closed with two rows — raw FAIL, adjusted PASS.

## T1.4 — must-FAIL 2: random gene sets

90 `RANDOM_CONTROL__` columns matched on **training-only** per-gene mean, variance and PC1 loading —
stricter than same-size random sets. Graded against this run's own unpaired `detection_floor`
(0.30 for d2_h, 0.40 for d2_i).

**Primary, spec-literal statistic** — random image direction, i.e. the units the floor is actually
measured in:

| state | control median | control p95 | control max | exceedances / 90 |
|---|---:|---:|---:|---:|
| d2_h full_biology | −0.0695 | 0.0389 | 0.0719 | **0** |
| d2_h rna_biology | −0.0728 | — | — | **0** |
| d2_h wsi_biology | −0.0341 | — | — | **0** |
| d2_i full_biology | −0.0238 | — | — | **0** |
| d2_i rna_biology | −0.0233 | — | — | **0** |
| d2_i wsi_biology | −0.0132 | — | — | **0** |

0/90 in every state against a 5% ceiling. **The control fails as required.**

**Second statistic, an observation and not a gate** — a cross-validated *fitted* image direction per
target column. This is not on the floor's scale and grading it against the floor would be a category
error; it is reported because it is the number a reader cares about:

| state | random-control median | real-target median | **control / real** |
|---|---:|---:|---:|
| d2_h full_biology | 0.4821 | 0.6323 | **0.762** |
| d2_h rna_biology | 0.4788 | 0.6288 | **0.761** |
| d2_h **wsi_biology** | 0.2158 | 0.2804 | **0.770** |
| d2_i full_biology | 0.4740 | 0.5790 | **0.819** |
| d2_i rna_biology | 0.4762 | 0.5863 | **0.812** |
| d2_i **wsi_biology** | 0.1642 | 0.2164 | **0.759** |

**Covariate-matched random gene sets are read at 76–82% of the level real curated gene sets are read
at, in every state on both artifacts.** The gap is consistent and real, but three quarters of the
apparent per-target molecular legibility is reproduced by a gene set matched only on marginal mean,
variance and PC1 loading. Any per-gene-set legibility claim must henceforth be stated as a *difference
against the matched random control, with a CI*, never as an absolute correlation.

## T1.5 — must-FAIL 3: shuffled gene labels

Build note: `build_pbs_targets --shuffle-gene-labels` refits the dictionary and **can no longer run on
this machine** — the data config whose digest the frozen target file records is gone, and every
surviving config declares a cohort missing 249 split patients. `v2/build_shuffled_pbs_targets.py`
therefore **rebinds** the frozen `gene_basis` (permuting its rows, exactly what the flag does after the
fit) and refuses to write unless the unshuffled reconstruction matches the frozen scores at r ≥ 0.9999
per column. Measured: **0.99999999999999**.

**(ii) Attribution must collapse — PASSES.** Median |Spearman| between true and shuffled per-axis gene
rankings: 0.0069 / 0.0073 / 0.0077 over three shuffle draws (bar ≤ 0.05); the strictly harder
best-match statistic, max over all 128 shuffled axes, is 0.033. **This is true by construction** — the
shuffle permutes basis rows after the fit — and is a build-integrity check, not a finding. A non-null
here would have meant the shuffle never took effect.

**(i) Subspace must persist — PASSES, and the pass is the problem.**

| shuffle seed | held-out top-CCA, true | shuffled | paired difference | CI95 of difference | inside CI95 of true? |
|---|---:|---:|---:|---|---|
| 1 | 0.5411 | **0.5600** | −0.0189 | [−0.0489, +0.0384] | yes |
| 2 | 0.5411 | **0.5360** | +0.0051 | [−0.0564, +0.0418] | yes |
| 3 | 0.5411 | **0.4771** | +0.0640 | [−0.0260, +0.0988] | no |

CI95 of the true value [0.4874, 0.5962], 500 patient bootstrap draws. By the containment criterion the
control passes on 2/3 draws; by the paired difference — the statistic that actually decides whether
they differ — **all three CIs cover zero**, and on one draw the shuffled dictionary scores higher.

After the row permutation the target block **is** a spectrum-matched random projection of the same
expression matrix. So the honest reading of "the subspace persists" is: *any spectrum-matched
128-dimensional projection of this expression matrix is as legible as the fitted interventional
dictionary.* This is the T1.1 random-dictionary baseline arriving by a second road, and it agrees with
T1.4's 76–82%.

Method note: the containment test in the T1.5 spec is weak — a wide CI passes it by accident. The
paired difference CI is what should be quoted; both are reported here so the difference is visible
rather than chosen.

## T1.6 — must-FAIL 4: modality-shuffled pairing

`permutation_null` permutes `y` **within cancer strata**, so cancer-level structure survives and only
patient-level pairing is destroyed. Raised from the 50-permutation default to **2,000**, because
p is floored at 1/(n+1) and 1/51 = 0.0196 cannot support a headline.

| artifact | state | adjusted top-CCA | null median | null p95 | excess over null median | permutation p |
|---|---|---:|---:|---:|---:|---:|
| d2_h | full_biology | 0.8890 | 0.1465 | 0.1645 | 0.7425 | 0.0005 |
| d2_h | rna_biology | 0.8874 | 0.1463 | 0.1642 | 0.7411 | 0.0005 |
| d2_h | wsi_biology | 0.6052 | 0.1483 | 0.1685 | 0.4569 | 0.0005 |
| d2_i | full_biology | 0.8479 | 0.1466 | 0.1638 | 0.7013 | 0.0005 |
| d2_i | rna_biology | 0.8533 | 0.1468 | 0.1643 | 0.7065 | 0.0005 |
| d2_i | wsi_biology | 0.4703 | 0.1472 | 0.1659 | 0.3231 | 0.0005 |

`permutation_p = 0.0005 = 1/2001` throughout — **no permutation of two thousand reached the observed
value in any state.** The resolution is quoted with the p; it is never written "p < 0.05".

**The null median is 0.146–0.148, not zero, and this must accompany every channel number on the
project.** Destroying the pairing does not take a 16-component top-CCA to zero; it takes it to the
capacity floor that 16 components fitted on 2,766 patients produce by construction. Quoting
`wsi_biology = 0.4703` without saying that chance is 0.147 overstates the effect by a factor a
reviewer will find immediately.

## T1.7 — must-PASS controls

**(a) RNA→RNA, circular by construction.** `--require-rna-positive-control` passes on both artifacts;
`channel_gate_failures` is empty and `gates_pass` is true. Adjusted top-CCA 0.8874 / 0.8533.

**(b) Known-legible covariate, excluded from the adjustment set.** Substitution recorded, not silent:
MSI, TP53 and consensus subtype are all unusable from the data on disk — the TCGA PanCan clinical
mirror carries `microsatellite_instability` only as an assay-performed flag (7 `YES`, 74 `NO`, 6,171
`NONE`) with no MSI-H/MSI-L/MSS calls, and there is no mutation or subtype table on either machine.
Substituted **TCGA-BRCA ER status by IHC** (690 labelled, 528+/162−), with PR as a second anchor.

Pre-registered at 01:45 UTC, before the 01:47 run (`p1_evidence/inputs/PREREG_known_covariate.json`):
ER band **[0.78, 0.92]**, point estimate 0.86 — Naik et al. *Nat Commun* 2020;11:5727 (0.92 internal,
**0.86 TCGA external**); Rawat et al. *Sci Rep* 2020;10:7275; Shamai et al. *JAMA Netw Open*
2019;2:e197700; Couture et al. *npj Breast Cancer* 2018;4:30. PR band [0.70, 0.85].

| artifact | state | adjustment | within-cancer AUROC | CI95 | null p95 | verdict |
|---|---|---|---:|---|---:|---|
| d2_h | **wsi_biology** | raw | **0.8781** | [0.8457, 0.9115] | 0.546 | **PASS** |
| d2_h | wsi_biology | cancer+TSS | 0.8714 | [0.8379, 0.9055] | 0.546 | PASS |
| d2_i | **wsi_biology** | raw | **0.8667** | [0.8360, 0.8971] | 0.542 | **PASS** |
| d2_i | wsi_biology | cancer+TSS | 0.8644 | [0.8340, 0.8946] | 0.544 | PASS |
| d2_h | full_biology | raw | 0.9195 | [0.8901, 0.9455] | 0.544 | PASS |
| d2_i | full_biology | raw | 0.9401 | [0.9127, 0.9640] | 0.545 | PASS (marginal) |

The image-only numbers sit essentially on the literature point estimate. The measured within-cancer
chance level is **0.542–0.546, not 0.5**, so grading against an assumed 0.5 would have been wrong by
four points. Declared weaknesses, written into the pre-registration before the run: BRCA is a
*development* cancer in the maximal split, so this ran on `--partition all` and is in-distribution;
only one cancer carries the label, so within-cancer and pooled coincide and the lineage-guessing
protection is not exercised. The RNA-containing states reach 0.92–0.94, which is expected by
construction (ER status is close to a monotone function of *ESR1*) and carries no morphological claim.

**(c) Synthetic spike recovery.** `transmission_floor = 0.01` — the finest level on the grid — for
every state: the pipeline transmits a paired signal of r = 0.01 without destroying it.
`detection_floor` 0.30 (d2_h) / 0.40 (d2_i) is the conservative unpaired number and the only one
quotable as a detection limit. **Attenuation 0.974–1.039 across all six states**: the confound
adjustment does not destroy signal.

`observed_above_floor = 0` for every state and that is the correct answer, not a failure. The floor is
in single-random-direction units and `observed_matched_direction` (the real channel through a *random*
direction pair) is −0.028 to +0.036. The channel is concentrated in particular directions, so a random
pair sees nothing. `floor_scale = targeted_single_direction` is emitted for exactly this reason.

## T1.1 / T1.2 — the must-beat baseline table

Built and bound to the same 6,427-patient maximal cohort, gene order, log transform and
development-only fit discipline as the PBS block, via `v2/baseline_target_common.py`:

* **curated pathway** — `hallmark_scores_pancan.parquet` (pre-existing) and the 82 non-control columns
  of `frozen_rna_targets.npz`, scored above.
* **size- and spectrum-matched random dictionary** — `v2/build_random_dictionary_targets.py`, Haar
  orthonormal basis rescaled by the frozen dictionary's own singular values.
* **PCA expression basis** — `v2/build_pca_basis_targets.py`, fit on development rows only.
* **gene-label-shuffled PBS** — `v2/build_shuffled_pbs_targets.py`.

**Not built, and this is a gap:** the **text-prior (GenePT-style)** and **capacity-matched
cell-composition** blocks. Both need an external resource (a gene text-embedding table; a deconvolution
signature matrix) that is on neither machine. The plan flagged these two as the most likely to slip and
they slipped. `composition_attribution` in `claim_guards.py` therefore remains undischarged.

**Degenerate by construction:** the zero-parameter cancer-type-mean baseline. The maximal split holds
out **whole cancers**, so every test patient receives the global training mean and the state is
constant on the evaluated partition. That is a property of the split, not a bug; the comparison is only
meaningful on `--partition all` and is recorded as such.

All four blocks were scored through the identical instrument (`--partition test --levels
0.0,0.05,0.10,0.20,0.40 --n-draws 16 --n-components 16 --n-permutations 500 --seed 42`; these are
`GateLedger.observe` comparisons, not validity gates, so the 2,000-permutation requirement binding the
must-FAIL controls does not bind here), and **every difference carries a 400-draw paired bootstrap CI**
— PBS and the baseline scored on the *same* resample, so the CI is on the difference.

**Held-out top-CCA, PBS minus baseline:**

| artifact | state | baseline | PBS | baseline | difference | CI95 | verdict |
|---|---|---|---:|---:|---:|---|---|
| d2_h | **wsi_biology** | random dictionary | 0.5032 | 0.4551 | +0.0481 | [−0.0177, +0.0693] | TIE |
| d2_h | **wsi_biology** | **PCA basis** | 0.5032 | **0.5520** | **−0.0488** | **[−0.0734, −0.0183]** | **BASELINE WINS** |
| d2_h | **wsi_biology** | shuffled labels (s1) | 0.5032 | 0.5146 | −0.0114 | [−0.0601, +0.0445] | TIE |
| d2_h | **wsi_biology** | shuffled labels (s2) | 0.5032 | 0.5187 | −0.0155 | [−0.0561, +0.0318] | TIE |
| d2_h | full_biology | random dictionary | 0.8417 | 0.8102 | +0.0315 | [+0.0241, +0.0653] | PBS wins |
| d2_h | full_biology | **PCA basis** | 0.8417 | **0.8776** | **−0.0359** | **[−0.0483, −0.0236]** | **BASELINE WINS** |
| d2_h | full_biology | shuffled labels (s1) | 0.8417 | 0.8140 | +0.0277 | [+0.0197, +0.0660] | PBS wins |
| d2_h | full_biology | shuffled labels (s2) | 0.8417 | 0.8085 | +0.0332 | [+0.0116, +0.0632] | PBS wins |
| d2_i | **wsi_biology** | random dictionary | 0.4605 | 0.4108 | +0.0497 | [+0.0251, +0.1372] | PBS wins |
| d2_i | **wsi_biology** | **PCA basis** | 0.4605 | 0.4905 | −0.0300 | [−0.0429, +0.0053] | TIE |
| d2_i | **wsi_biology** | shuffled labels (s1) | 0.4605 | 0.4245 | +0.0360 | [+0.0111, +0.0926] | PBS wins |
| d2_i | **wsi_biology** | shuffled labels (s2) | 0.4605 | 0.4317 | +0.0288 | [−0.0030, +0.0887] | TIE |
| d2_i | full_biology | random dictionary | 0.8634 | 0.8487 | +0.0147 | [−0.0004, +0.0283] | TIE |
| d2_i | full_biology | **PCA basis** | 0.8634 | **0.8714** | **−0.0080** | **[−0.0233, −0.0001]** | **BASELINE WINS (marginal)** |
| d2_i | full_biology | shuffled labels (s1) | 0.8634 | 0.8408 | +0.0227 | [−0.0029, +0.0288] | TIE |
| d2_i | full_biology | shuffled labels (s2) | 0.8634 | 0.8374 | +0.0260 | [+0.0057, +0.0383] | PBS wins |

Reading down the three opponents:

* **vs. size- and spectrum-matched random dictionary** — PBS wins 2/4, ties 2/4, **never loses**.
* **vs. gene-label-shuffled PBS** — PBS wins 3/8, ties 5/8, **never loses**.
* **vs. ordinary PCA of the same expression matrix** — **PBS loses 3/4 with a CI excluding zero and ties
  the fourth. It never wins.**

**This is the deflationary result T1.1 was designed to force.** The PCA basis is fit on development rows
only, through the same transform PBS uses, and is capacity-matched at 128 columns, so neither leak nor
capacity explains it. The claim "interventional coordinates are legible in a way generic decompositions
are not" is contradicted by the strongest available comparison, with a CI. What survives is the much
weaker "interventional coordinates beat *random* projections of the same expression matrix" — which PCA
also satisfies.

**A readout dependence that qualifies §T1.5.** At `n_components = 32` with a 500-draw bootstrap the
gene-label-shuffled block was indistinguishable from PBS on all three shuffle draws. At
`n_components = 16` here, PBS beats it on 3 of 8 draws with a CI excluding zero. The defensible
statement is therefore *the shuffle costs the dictionary little and sometimes nothing, never more than a
few hundredths of a CCA* — still enough to void gene-level attribution, but not the flat tie the
n=32 readout alone suggested. Both readouts are reported; neither is selected after the fact.

Block-level context (adjusted top-CCA, and the same value minus the ~0.147 permutation null median):

| block | d2_h wsi_biology | above null | d2_i wsi_biology | above null | induced baseline (d2_h) |
|---|---:|---:|---:|---:|---:|
| curated pathway, 82 cols | 0.6052 | 0.4569 | 0.4703 | 0.3231 | 0.1062 |
| PBS, 128 | 0.5504 | 0.4050 | 0.5217 | 0.3730 | 0.0705 |
| random dictionary, 128 | 0.5226 | 0.3761 | 0.4415 | 0.2919 | 0.0806 |

The induced baseline differs by block, which is itself a Track 2 confirmation: blocks whose coordinates
the confound design explains less produce less induced correlation.

## T1.8 — the ledger

`v2/calibra/track1_battery_ledger.py` assembles every row into
`v2/research/rebase/nature/GATE_LOG.md` through `GateLedger`, and the add/observe separation is
**enforced and tested**, not merely described:

* must-FAIL and must-PASS controls → `GateLedger.add` (validity conditions, decide the verdict);
* must-beat baseline comparisons → `GateLedger.observe` (scientific outcomes, excluded from the
  verdict).

`v2/tests/test_track1_battery_ledger.py::test_a_losing_baseline_is_an_observation_and_cannot_move_the_verdict`
plants a baseline that beats us 0.95 to 0.40 and asserts the verdict is unchanged against the identical
fixture without it. Registering that as a gate would mark the run FAILED exactly when the science came
back "no", making a true negative indistinguishable from a broken pipeline — which is the separation
P1 claims to enforce, so the battery has to demonstrate it. A missing control is written with
`value = NaN`, `note = "inadmissible_<code>"` and status FAIL rather than omitted; silence is never a
pass.

Current ledger: **81 rows — 62 gates, 19 observations, 7 failed gates** (the six raw site
certificates, and gene-label-shuffle seed 3's containment).

## Tests

`v2/tests/test_track1_controls.py` (inherited, audited), `v2/tests/test_gene_label_shuffle_control.py`
(9 new), `v2/tests/test_track1_battery_ledger.py` (7 new). The inherited agent's reported "single-column
bug" was located and its fix **is present**: `run_calibra.score_target_block_per_column` reshapes
sklearn's single-output prediction (`.reshape(len(test_idx), -1)`), covered by
`test_block_scoring_matches_the_single_column_reference`.

---

**Logged:** 2026-08-03, 01:20–03:30 UTC.

**How obtained:** Lambda A100 box `ubuntu@150.136.45.194`, workspace `~/ws_p1` (`morpheus` →
`~/morpheus-rebase-p1`), `~/venv/bin/python`, CPU only. Commands:
`confound_certificate --partition test --n-permutations 1000 --n-boot 200 --n-boot-axes 8` (with and
without `--residualise`); `run_calibra --partition test --levels 0.0,…,0.50 --n-draws 40
--n-components 16 --n-permutations 2000 --seed 42 --score-random-controls
--require-rna-positive-control`; `gene_label_shuffle_control --n-boot 500`;
`known_covariate_control --partition all --n-boot 1000 --n-permutations 1000`;
`track1_battery_ledger`. All outputs under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/` (persistent storage).

### Technical

Six must-FAIL/must-PASS controls at ≥1,000 permutations on 2,766 held-out patients, 108 confound
columns. Raw states predict pooled tissue source site at joint out-of-fold balanced accuracy
0.235–0.363 against a chance rate of 0.0118 and a permutation-null p95 of 0.122–0.186, p = 1/1001 in
every state; the identical certificate on cancer+TSS cross-fitted residuals returns 0.0052–0.0118 with
zero breaching axes. Modality-shuffled pairing gives p = 1/2001 in every state against a null median of
0.1465–0.1483. Matched random gene-set controls give 0/90 exceedances of the detection floor in the
floor's own random-direction units, but reach 0.759–0.819 of the real targets' value under a fitted
direction. A gene-label-shuffled dictionary is statistically indistinguishable from the fitted one
(paired bootstrap difference CI covers zero on 3/3 draws). Pre-registered ER-status recovery: 0.8781
[0.8457, 0.9115] against [0.78, 0.92]. Transmission floor 0.01, attenuation 0.974–1.039.

### In plain terms

The instrument does the things a working instrument must do: scrambling which patient's slide goes with
which patient's RNA destroys the signal completely, planting a fake signal of known size gets it back
intact, randomly chosen gene sets do not clear the detection threshold, and a covariate whose
difficulty four published papers had already measured came back at almost exactly the published
difficulty.

Two checks came back the wrong way. Our raw axes *can* tell which hospital a slide came from — not
strongly on any one axis, but a simple classifier using all of them picks the right hospital about a
third of the time out of eighty-five. The correction we already apply before every measurement removes
this completely, so our published numbers are safe, but the raw axes must never be shown to anyone.
And when we scramble which gene each of our dictionary coordinates refers to, nothing gets worse —
which means the specific biology in the dictionary is not what the images are reading. Randomly chosen
gene sets score about three quarters of what real ones do, by an entirely separate measurement, and
agree.

### Meaning for the claim

* **P1 is strengthened.** A battery that only ever passes proves nothing. This one caught a real
  confound leak, located it in the raw representation rather than in the instrument, and showed the
  adjustment discharging it — and it did so through the joint test that the spec as written would have
  omitted. That last point is a contribution to the certification rule itself.
* **P3 and P4 are damaged in the same place, twice, by independent routes.** T1.5 (a shuffled
  dictionary is as legible as the real one) and T1.4 (matched random gene sets reach 76–82% of real
  ones) both say the channel is largely non-specific to the named biology. **P3 may not claim
  "interventional coordinates beat curated pathway scores" on this evidence**, and P4's `inspect_gene`
  cannot ship, until a comparison exists in which the fitted basis wins with a CI excluding zero.
* **P4 gains two required certificate fields:** the joint site test, and `certified_on = {raw |
  adjusted}`.
* **P2 is unblocked**: the positive controls passed in the same run through the same code path, which
  is what makes its negative result admissible; and it now has a measured chance level (null median
  0.147, p95 0.164–0.169) rather than an assumed zero.
* **Outstanding:** text-prior and cell-composition baseline blocks (external resources absent);
  the full baseline table (running); a second seed for every control.
