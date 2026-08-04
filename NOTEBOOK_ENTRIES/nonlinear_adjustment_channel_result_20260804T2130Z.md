## 2026-08-04 21:30 UTC — The channel does not move. A **saturated** cell-mean adjustment — the provable upper limit of every nonlinear adjustment on this design — takes P1 §4.4's `d2_h::wsi_biology` channel from **0.6052 to 0.6051**; a representation that is *nothing but* the confound labels, pushed through the identical pipeline, scores **0.0903** against the channel's own null median of **0.1483**. And the residual the probe found is real but is **2.4–3.3×** its own *measured* chance, not the 4.3–4.9× reported, because a cross-fitted adjustment manufactures ~1.5× of that multiple and the probe's null cannot see it

**Logged:** 2026-08-04 21:30 UTC. **Predeclared:**
`NOTEBOOK_ENTRIES/PREDECLARED_nonlinear_adjustment_channel_20260804T2015Z.md`, committed `010d9c4`,
**before** the instrument existed and before any number below existed. **How obtained:**
`v2/calibra/nonlinear_adjustment.py` (new, 25 tests) on the A100 box (150.136.45.194), workspaces
`~/ws_nla{,2,3,4}` and `~/ws_nlabase`, each deployed by `git -c core.autocrlf=false archive HEAD` and
verified file-by-file against `git ls-tree -r HEAD` blob SHA-1 — **669 files, 0 mismatches, every
time**. Thread caps `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1`, process parallelism, CPU only; the box
carried a co-tenant load of 20–50 on 30 cores throughout. Outputs
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/nonlinear_adjustment/`.

---

### 0. The answer, in the order the brief asked for it

**Does the channel survive nonlinear adjustment, and by how much? It survives entirely — retention of
excess-over-null 0.987–1.007 under every arm — but by the letter of the predeclaration this is
Reading 4 (NO TEST), not Reading 1, and that distinction is kept.** No arm cleared the predeclared
step-3 validity bar (k-NN ≤ 2.0× chance on both targets), and §6 shows *why no arm could*: the
measured floor for **any** cross-fitted adjustment on this cohort is 1.45–1.75× design chance,
because the adjustment itself puts structure there. The bar was unattainable when it was written.
That is stated as a defect of the bar, not as a pass.

The question the brief actually asks — *is the residual large enough to produce the channel?* — is
answered anyway, from three directions that do not depend on step 3:

1. **The saturated limit (§4).** `confound_design(["cancer","tss"])` is purely one-hot and purely
   additive, and TCGA site nests inside cancer. A linear model on the **saturated (cancer × site)
   cell design** therefore spans every function of the confound labels that exists, and upper-bounds
   what any kernel-ridge, forest or boosting estimate of `E[X | C]` can remove. Measured: the channel
   goes **0.6052 → 0.6051**, retention **0.998**. This is not "we tried a nonlinear model and it did
   not help"; it is "no conditional-mean adjustment on this design *can* help, and here is the limit".
2. **The ceiling from the other side (§5).** A representation that is *nothing but* the confound
   labels, pushed through the identical pipeline, scores **0.0903** (saturated cell design) or
   **0.1237** (additive design) — against the real channel's own null median of **0.1483**. In excess
   over each block's own null, that is **6.0%** and **11.2%** of the channel's excess.
3. **The residual, re-measured against a null that can see the adjustment (§6).** Still real, still at
   the permutation floor, but **2.4–3.3× chance for site** and **3.5–3.7× for cancer**, not 4.3–4.9×.

**Verdict against the predeclared bands: Reading 4 on the letter**, with the question answered by §4
and §5 instead. Reading 4 was predeclared as *not* a favourable result and is not written as one.

---

### 1. The reproduction gate, discharged before anything else was believed

The predeclaration made everything conditional on arm A1 — `residualise.cross_fitted_residuals`, the
incumbent — reproducing P1 §4.4's published row on the same artifact. It does, on **all three**
quantities and **both** artifacts, to four decimal places:

| artifact | state | statistic | published §4.4 | this run | match |
|---|---|---|---:|---:|---|
| d2_h | wsi_biology | adjusted top-CCA (S1) | 0.6052 | **0.6052** | ✔ |
| d2_h | wsi_biology | within-cancer pairing null median | 0.1483 | **0.1483** | ✔ |
| d2_h | wsi_biology | excess over null median | 0.4569 | **0.4569** | ✔ |
| d2_i | wsi_biology | adjusted top-CCA (S1) | 0.4703 | **0.4703** | ✔ |
| d2_i | wsi_biology | within-cancer pairing null median | 0.1472 | **0.1472** | ✔ |
| d2_i | wsi_biology | excess over null median | 0.3231 | **0.3231** | ✔ |

*Artifacts `runs/d2_final/artifacts/d2_{h,i}_seed42.npz` — the copies the probe entry §6 showed
reproduce §4.2, in preference to `~/e0_run/d2_v3/*`. Targets
`morpheus_phase_d/data/frozen_rna_targets.npz`, 180 columns of which the 90 `RANDOM_CONTROL__` ones
are excluded exactly as `run_calibra.py:364-371` does, leaving 90. `--partition test`, n = 2,766
after alignment, 0 duplicate patient IDs, design 108 columns of rank 105, 84 pooled sites kept at
`min_site_count=10`, 105 observed (cancer × site) cells, `n_components=16`, `seed=42`,
2,000 permutations, p floor 1/2001 = 0.0005.*

The **step-3 probe** reproduces the previous entry as exactly, on a fresh workspace: `d2_h`,
`adjusted` arm, site **0.0511 = 4.34× chance at k3**, global p 0.0050, global p95 **0.0154**,
within-cancer p95 **0.0399**; cancer **0.2126 = 4.46× at prior-corrected k5**, p 0.0050, p95
**0.0582**. Those are the previous entry's `d2_h | adjusted` rows digit for digit.

**Both statistics are reported throughout**, because the brief's "0.58–0.62 held-out against ~0.147"
merges two of them and P2_RANK_DRAFT §… already warns that 0.140 and 0.147 must not be interchanged:

* **S1** = `spectral.top_canonical_correlation(adjust(x), adjust(y), n_components=16)` — the §4.4
  headline; an in-sample maximum over 16 whitened directions per side, computed on the held-out
  patient partition. Its null is the within-cancer **pairing** permutation, median 0.147.
* **S2** = `spectral.heldout_top_cca(..., seed=42)` — directions fitted on one half of the held-out
  patients, scored on the other. **d2_h reads 0.5841**, which is also inside "0.58–0.62"; d2_i 0.4206.

---

### 2. Which null applies, on both sides, and why they differ

The brief flags that `within_stratum_permutations` is degenerate when the confound nests inside the
stratum. That is true of the **probe's label null** and **not** of the **channel's pairing null**, and
the two must not be conflated:

* **Probe (§6, §7).** Permutes a *label* inside cancer. TCGA site nests totally inside cancer
  (`nesting_diagnostic` on these exact rows: 21 strata, 1 with a single class, 0.0101 of rows
  unpermutable, **1 of 85 classes spanning more than one stratum**, `within_stratum_chance` =
  **0.1739** against a design chance of 0.011765). The **global** null is the applicable bar and it
  *measures* chance; the within-cancer figure is carried beside the site rows for decomposition only.
* **Channel (§3, §4).** Permutes the *patient pairing* between two blocks and re-adjusts `y` on every
  draw. Nothing nests, nothing is degenerate. The project's convention is within-cancer and it is
  quoted; the **global pairing null was also computed at 2,000 permutations for every arm** and the
  two agree to ~0.003 (ridge: within-cancer 0.1483, global 0.1457; excess 0.4569 vs 0.4596), so no
  conclusion here turns on the choice.

---

### 3. Step 1 — the channel re-measured under every arm

`wsi_biology`, test partition, n = 2,766, S1 at 16 components, within-cancer pairing null at 2,000
permutations, p floor 0.0005. **Retention** is `nonlinear_adjustment.retention_of_excess` — a ratio of
**excess over each arm's own null**, never of raw S1, for the reason P1 §4.4 gives.

**`d2_h_seed42`**

| arm | what it removes | S1 | own null median | excess | **retention** | S2 held-out | eff. rank | vs incumbent (per-axis r) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| none | column means only | 0.8155 | 0.7102 | 0.1053 | — | 0.8062 | 19.65 | 0.744 |
| **ridge (incumbent)** | additive one-hot, α=1 | **0.6052** | **0.1483** | **0.4569** | **1.000** | 0.5841 | 22.50 | 1.000 |
| **saturated cell mean** | **every function of the labels** | 0.6051 | 0.1491 | 0.4560 | **0.998** | 0.5835 | 22.50 | 0.9992 |
| kernel ridge α=1 γ=0.25 | RBF on the one-hot design | 0.6052 | 0.1473 | 0.4579 | **1.002** | 0.5837 | 22.50 | ≈1 |
| kernel ridge α=1 γ=0.5 | RBF on the one-hot design | 0.6052 | 0.1477 | 0.4575 | **1.001** | 0.5837 | 22.50 | 0.9995 |
| location–scale (2nd moment) | per-cell mean **and** SD | 0.5978 | 0.1468 | 0.4510 | **0.987** | 0.5719 | 42.41 | 0.9889 |
| in-sample saturated | cell means, exactly zero | 0.6054 | 0.1453 | 0.4601 | **1.007** | 0.5825 | 22.55 | 0.9943 |

**`d2_i_seed42`**

| arm | S1 | own null median | excess | retention | S2 held-out | eff. rank |
|---|---:|---:|---:|---:|---:|---:|
| none | 0.7700 | 0.6672 | 0.1028 | — | 0.7706 | 10.78 |
| **ridge (incumbent)** | **0.4703** | **0.1472** | **0.3231** | **1.000** | 0.4206 | 11.73 |

*(full d2_i arm set in the JSON; see §11.)*

Every `permutation_p` in both tables is **0.0005 = 1/2001**, the floor: no permutation of two thousand
reached the observed value under any adjustment.

**Predeclared distrust items 3, 5 and 6, discharged on this table.** (3) The nulls barely move across
arms (0.1453–0.1491), so grading on excess and grading on raw S1 give the same answer here — reported
anyway. (5) Effective rank *rises* under the location–scale arm (22.50 → 42.41) rather than falling
toward the 16-component budget, so no arm's S1 is propped up by a rank collapse. (6) S1 and S2 agree
in direction and in magnitude of change for every arm — S2 moves 0.5841 → 0.5719 over the same span
that S1 moves 0.6052 → 0.5978.

---

### 4. Why "adjust nonlinearly" has a provable answer on this design, not merely an empirical one

`confound_design(frame, ["cancer", "tss"])` produces 21 cancer indicators + 84 pooled-site indicators
+ `OTHER` + two `dummy_na` columns — **purely one-hot, purely additive**. For a categorical predictor
set a linear model on a one-hot design already spans the conditional-mean model, so the only thing a
nonlinear learner can add is the **cancer × site interaction**. And that interaction is nearly
degenerate here, because site nests inside cancer (§2): the site indicator already determines the cell
for all but the 1.0% of rows in pooled `OTHER`.

`saturated_cell_residuals` removes the exact conditional mean of every one of the 105 observed cells —
the column space of that design is *every* function of `(cancer, site)`. **It is an upper bound on
what any kernel ridge, random forest or gradient-boosting estimate of `E[X | C]` can remove, and it
changes the channel by 0.0001.** The kernel arm, run over the predeclared six-cell grid, agrees to
four decimals and its residuals correlate with the incumbent's at a per-axis median of 0.9992–0.9995
(predeclared distrust item 2: **these arms are relabelled incumbents, and are reported as such**).

This was written down as a prediction in the predeclaration §2 before any of it was run, precisely so
that it could not be offered afterwards as an excuse. It is now measured.

**Consequence for the brief's step 1.** "Build a cross-fitted nonlinear adjustment against the same
`cancer + pooled TSS` design and re-measure the channel" is a well-posed instruction whose answer on
this design is fixed by arithmetic: the channel cannot move, because there is no additional
conditional mean to remove. The informative version of the question is the one §6 asks instead.

---

### 5. Step 2 — the labels-only ceiling, the bound from the other side

How much channel is obtainable from `cancer + pooled TSS` **alone**, with no image features. The
confound design stands in for the image representation and the **same 16-component statistic** is
read, so capacity is matched rather than assumed comparable. This is a property of the cohort, not of
an artifact, and duly reads identically for `d2_h` and `d2_i`.

| "representation" | raw S1 | raw S2 | **after the incumbent adjustment: S1** | its own null median | its p | **excess** |
|---|---:|---:|---:|---:|---:|---:|
| additive design (108 cols, rank 105) | **0.9273** | 0.9234 | **0.1237** | 0.0723 | 0.0010 | 0.0514 |
| saturated cell design (105 cols, rank 105) | **0.7722** | 0.7628 | **0.0903** | 0.0631 | 0.0080 | 0.0272 |
| **`wsi_biology` (the real channel)** | 0.8155 | 0.8062 | **0.6052** | 0.1483 | 0.0005 | **0.4569** |

Cross-fitted R² of each block on the labels: **targets 0.4416**, image **0.4139** (`d2_h`) / **0.3623**
(`d2_i`).

Three things this says.

1. **Raw, the confound labels carry more channel than the image representation does** — 0.9273 against
   0.8155. Confounding on this cohort is not a subtle worry; before adjustment the labels beat the
   image. That is the strongest possible motivation for the adjustment CALIBRA applies.
2. **After adjustment they carry almost nothing.** 0.1237 and 0.0903 sit *below* the real channel's own
   null median of 0.1483. As a share of the channel's excess over its own null, the ceiling is
   **11.2%** (additive) and **6.0%** (saturated). The additive figure is the larger only because
   residualising a design on *itself* leaves ridge-shrinkage remainder at α = 1; the saturated figure
   is the one that carries genuine un-removed interaction, and it is the smaller.
3. Both remain **significant against their own nulls** (p = 0.0010 and 0.0080), so this is a bound, not
   a zero. Reported as a bound.

**The honest limit of this bound, stated rather than buried.** It bounds what a representation that is
a *function of the labels* can contribute. The residual §6 measures is not a function of the labels in
the mean — it is higher-moment. Two things narrow that gap: a top canonical correlation is a
**second-moment** statistic between two blocks, and conditional-variance heterogeneity that is
independent across the two blocks contributes nothing to their cross-covariance; and whatever it does
contribute is regenerated inside the channel's own null, because `permutation_null` re-adjusts `y` on
every permutation (`calibration.py:158`). Neither is a proof, and this is written as the weakest link
in §0's argument.

---

### 6. Step 3 — the validity check, and the reason it cannot be passed

**The check, as predeclared.** `nonlinear_confound_probe.probe_state` re-run on each arm-adjusted
`wsi_biology`: k ∈ {1,3,5,10,15,25,50}, plain **and** inverse-training-frequency vote, reading = max
over k and vote rule, 200 global permutations, p floor 0.0050. Bar: **verified** at ≤ 2.0× design
chance on both targets.

| arm | site k-NN max | × design chance | cancer k-NN max | × design chance | verdict |
|---|---:|---:|---:|---:|---|
| ridge (incumbent) | 0.0511 | 4.34 | 0.2126 | 4.46 | reference |
| saturated cell mean | 0.0497 | 4.23 | 0.2065 | 4.34 | **not verified** |
| kernel ridge α=1 γ=0.25 | 0.0546 | 4.64 | 0.2067 | 4.34 | **not verified** |
| kernel ridge α=1 γ=0.5 | 0.0518 | 4.40 | 0.2096 | 4.40 | **not verified** |
| location–scale | 0.0494 | 4.20 | 0.2102 | 4.41 | **not verified** |
| in-sample saturated | 0.0603 | 5.13 | 0.2125 | 4.46 | **not verified** |

Every global p is 0.0050, the floor. **No arm reduces the recovery; the strongest first-moment
adjustment that exists reduces it by 2.5%, and removing the second moment as well reduces it by 3.2%.**

**Why. The chance rate for an adjusted block is not `1/n_classes`, and the probe's null cannot
measure it.** The asymmetry is in the code:

* `calibration.permutation_null` — the **channel's** null — permutes the pairing and then
  **re-residualises `y`** on every permutation (`calibration.py:158`). Any correlation the shared
  residualisation *induces* is therefore regenerated inside the null. This is P1 §4.6's
  induced-correlation floor, and it is properly handled.
* `nonlinear_confound_probe.probe_state` — the **probe's** null — permutes the **labels** of an
  already-adjusted feature matrix. Structure that the adjustment tied to the *true* cells is broken in
  the null and intact in the observed, so it is scored as surviving confound.

The previous entry's §7 item 10 argues the opposite — *"the global permutation null is computed on the
same already-residualised features, so any structure the residualiser itself introduced is inside the
null as well as inside the observed"*. **That is not so, and the correction is one of this run's two
main results.** A cross-fitted residual leaves every (cell × fold) group displaced by that fold's
estimation error `μ_c − μ̂_c^(−f)`, an offset of order `σ/√n_c` **shared by every patient in the
group**; permuting labels destroys the correspondence between group and label in the null while
leaving it intact in the observed.

**The mechanism, measured** (`cross_fitting_offset_energy`, `d2_h::wsi_biology`; the expected share for
G groups with no structure is G/n, emitted beside every value):

| arm | cell-mean energy | expected | ratio | (cell × fold)-mean energy | expected | ratio |
|---|---:|---:|---:|---:|---:|---:|
| none (raw) | 0.46552 | 0.03796 | 12.3× | 0.55045 | 0.18981 | 2.90× |
| ridge | 0.00245 | 0.03796 | **0.065×** | 0.22788 | 0.18981 | **1.20×** |
| saturated | 0.00129 | 0.03796 | 0.034× | 0.23298 | 0.18981 | 1.23× |
| in-sample saturated | **0.00000** | 0.03796 | **0.000×** | 0.15891 | 0.18981 | 0.84× |

The signature is exact: the adjustment drives the cell means to **6.5%** of what random grouping would
give (the anti-correlation cross-fitting is known to induce), and simultaneously drives the
cell × fold means to **120%** of theirs. In-sample removal zeroes the cell means exactly — the normal
equations — and leaves cell × fold *below* baseline.

**The control.** Row-shuffle the block **before** the adjustment: same rows, same covariance, same
effective rank, no confound association left. Anything the probe still reads was put there by the
adjustment. `d2_h::wsi_biology`, incumbent arm, same probe, same nulls:

| block | site k-NN max | × chance | global p | cancer k-NN max | × chance | global p |
|---|---:|---:|---:|---:|---:|---:|
| **observed** (real block, adjusted) | 0.0511 | 4.34 | 0.0050 | 0.2126 | 4.46 | 0.0050 |
| **row-shuffled, then adjusted** | 0.0186 | **1.58** | **0.0249** | 0.0713 | **1.50** | **0.0149** |
| row-shuffled, **not** adjusted | 0.0158 | 1.34 | 0.1045 | 0.0609 | 1.28 | 0.0647 |

**The adjustment turns a non-significant reading into a significant one on a block with no confound
left in it at all.** (The 1.28–1.34× of the unadjusted shuffled row is the upward bias of a
max-over-fourteen-probes statistic, and is not significant against its own null.)

Independently, in simulation at the same shape — n = 2,766, 85 classes with a 30% `OTHER` class, 256
axes, **no confound whatsoever** — a cross-fitted residualisation against the cell design yields a
k-NN reading of **2.15× chance at p = 0.024** (its 1/41 floor), while the same noise block raw reads
**0.97× at p = 0.54** and in-sample-residualised reads **0.77× at p = 0.95**.

**The corrected null.** `regenerated_adjustment_null` permutes the **rows of the block before the
adjustment**, so every draw carries an adjustment artefact of its own size, and takes the max over the
same k grid and both vote rules in the null as in the observed. `corrected_multiple` is the observed
divided by that null's median — the **measured** chance rate under the adjustment rather than the
assumed one. 200 permutations, p floor 0.0050.

**`d2_h_seed42`, `wsi_biology`, test partition**

| arm | target | observed | × design chance | regenerated null median | × design chance | **corrected multiple** | p |
|---|---|---:|---:|---:|---:|---:|---:|
| ridge (incumbent) | site | 0.0511 | 4.34 | 0.0187 | 1.59 | **2.73** | 0.0050 |
| ridge (incumbent) | cancer | 0.2126 | 4.46 | 0.0580 | 1.22 | **3.67** | 0.0050 |
| saturated | site | 0.0497 | 4.23 | 0.0200 | 1.70 | **2.49** | 0.0050 |
| saturated | cancer | 0.2065 | 4.34 | 0.0582 | 1.22 | **3.55** | 0.0050 |
| kernel ridge α=1 γ=0.25 | site | 0.0546 | 4.64 | 0.0176 | 1.49 | **3.10** | 0.0050 |
| kernel ridge α=1 γ=0.25 | cancer | 0.2067 | 4.34 | 0.0575 | 1.21 | **3.59** | 0.0050 |
| kernel ridge α=1 γ=0.5 | site | 0.0518 | 4.40 | 0.0180 | 1.53 | **2.88** | 0.0050 |
| kernel ridge α=1 γ=0.5 | cancer | 0.2096 | 4.40 | 0.0579 | 1.22 | **3.62** | 0.0050 |
| location–scale | site | 0.0494 | 4.20 | 0.0205 | 1.75 | **2.41** | 0.0050 |
| location–scale | cancer | 0.2102 | 4.41 | 0.0594 | 1.25 | **3.54** | 0.0050 |
| in-sample saturated | site | 0.0603 | 5.13 | 0.0181 | 1.54 | **3.33** | 0.0050 |
| in-sample saturated | cancer | 0.2125 | 4.46 | 0.0571 | 1.20 | **3.72** | 0.0050 |

**`d2_i_seed42`, incumbent:** site observed 0.0507 (4.31× design chance) against a regenerated null
median of 0.0170 (1.45×), **corrected 2.98×**, p = 0.0050.

**What this changes and what it does not.**

* **It does not overturn the previous entry's direction.** The residual is real. Every corrected
  multiple is ≥ 2.4 and every p is at the floor, against a null that contains the adjustment.
* **It does change the magnitudes.** The publishable multiple for site on `d2_h` is **2.73×**, not
  4.34×; for cancer **3.67×**, not 4.46×. Between **27% and 40%** of the reported multiple is the
  adjustment reading itself back.
* **It explains why nothing passes step 3.** The floor is 1.45–1.75× design chance for *any*
  cross-fitted adjustment, so a bar of "≤ 2.0× chance" leaves a window of 0.25–0.55× that no
  adjustment can reach. The predeclared bar was not attainable. That is a defect of the bar and it is
  reported as one rather than converted into a pass.

---

### 7. Where the surviving residual is *not*, which narrows what it can be

Three adjustments bracket the residual from three sides, and none of them touches it:

* **First moment, maximal flexibility** — the saturated cell-mean arm removes every function of the
  labels. Recovery 4.34× → 4.23×.
* **First moment, exactly zero** — the in-sample saturated arm makes every cell mean exactly zero by
  the normal equations. Recovery 4.34× → **5.13×**, i.e. it goes *up*. (It also removes the
  cross-fitting artefact, and its corrected multiple 3.33× is correspondingly higher than the
  incumbent's 2.73×.)
* **Second moment, per axis** — the location–scale arm removes the per-cell SD as well. Recovery
  4.34× → 4.20×, corrected 2.73× → 2.41×.

So the surviving structure is **not** a conditional mean of any shape, and **not** per-axis conditional
variance. What remains available to it: conditional **covariance** (a rotation of the within-cell
cloud), higher moments, or conditional support/manifold structure. A per-cell whitening would test
covariance directly and **was not run**: 105 cells hold ~26 patients each against 256 axes, so the
per-cell covariance is rank-deficient by an order of magnitude and any shrinkage strong enough to make
it estimable would make the result uninterpretable. That is the concrete next experiment and it needs
a bigger cohort or a lower-dimensional block, not more compute.

---

### 8. The over-removal guard was not triggered, and why that is reportable

Predeclaration §6: if any arm reached Reading 2 or the lower half of Reading 3, its attenuation would
be measured through the identical adjustment with `calibration.spike_recovery_curve` before a collapse
could be reported as a confound finding. **No arm collapsed** — the lowest retention is 0.987 — so the
guard did not trigger and the spike curves were not run. Stated so that their absence is a
consequence of the result rather than an omission.

---

### 9. Every predeclared distrust item, discharged

The predeclaration listed eight ways a **favourable** result would be untrustworthy and three ways an
**unfavourable** one would be.

1. **Step 3 fails → Reading 4, not Reading 1.** Honoured. The headline sentence of §0 says the channel
   survives *and* says this is Reading 4.
2. **The adjuster is a relabelled incumbent.** `adjuster_agreement` is reported for every arm. The
   kernel arms sit at per-axis r-median 0.9992–0.9995 and **are** relabelled incumbents; the flag says
   so in the JSON. The location–scale arm is not (relative Frobenius difference 0.9919).
3. **The null moved.** It did not, materially: 0.1453–0.1491 across arms. Grading is on excess anyway.
4. **The ceiling is large.** It is not: 6.0–11.2% of the channel's excess, §5, and it is in the
   headline table rather than an appendix.
5. **Capacity.** Effective rank of every adjusted block reported; it rises rather than falls.
6. **S1 vs S2.** Both reported for every arm; they agree.
7. **Artifact disagreement.** `d2_h` and `d2_i` agree in direction on every arm run on both.
8. **Folds.** The adjuster's `KFold(5, seed=42)` and the probe's `_stratified_folds(seed=42)` are
   different partitions, so a probe test row can be an adjuster train row. Recorded, and — unlike the
   previous entry — **not** defended with the claim that the null contains it: §6 shows it does not,
   and `regenerated_adjustment_null` is the fix.

Unfavourable-side items: (1) over-removal — not triggered, §8; (2) A5 is not the same adjustment — the
location–scale arm is labelled a *different and stronger* operation everywhere it appears; (3) loss of
rank — measured, and it goes the other way.

---

### 10. Prose that is now wrong or incomplete, flagged and **not** edited

Per the rules for this run, `NOTEBOOK.md`, the paper drafts and `claim_guards.py` were not touched.

1. **`NOTEBOOK_ENTRIES/tcga_nonlinear_confound_probe_result_20260804T2100Z.md` §7 item 10:** *"the
   global permutation null is computed on the same already-residualised features, so any structure the
   residualiser itself introduced is inside the null as well as inside the observed."* **Not
   supported.** The null permutes labels of fixed features, so cell-tied structure is broken in the
   null and intact in the observed. Measured: a row-shuffled block with no confound association reads
   1.58× (site) / 1.50× (cancer) at p = 0.0249 / 0.0149 after the same adjustment.
2. **Same entry, every "× chance" for an *adjusted* arm** (§0, §3, §3b, §8 — the 4.3–4.9× site and
   3.4–4.9× cancer figures). The multiplier is against `1/n_classes`, which is the chance rate for a
   *raw* block. Against the adjustment's own measured rate the same readings are **2.4–3.3× (site)**
   and **3.5–3.7× (cancer)**. The *significance* is unchanged; the *magnitude* is overstated by
   27–40%. The raw-arm multiples are unaffected.
3. **Same entry §2, "Measured chance matches design chance, so the probe is not capacity-bound"** —
   correct as written (the global label-permutation null does sit at the design rate) but incomplete:
   that null measures the chance rate of the *estimator*, not of the *estimator applied to an adjusted
   block*. The sentence should name which.
4. **This run's own predeclaration §5**, the "verified ≤ 2.0× chance" bar, was not attainable by any
   cross-fitted adjustment. Recorded here rather than quietly restated.

Nothing in `paper/P1_CALIBRA_DRAFT.md` §4.4 or §4.6 is contradicted by this run; §4.4's row reproduces
exactly and §4.6's induced-correlation floor is precisely the phenomenon whose *probe-side* analogue
was missing. One **addition** worth making: `nonlinear_confound_probe`'s module docstring argues
correctly for measuring chance rather than assuming it, and should say that a label-permutation null on
an adjusted block measures the estimator's chance rate but **not** the adjustment's contribution to it.

---

### 11. Suite status, files, outputs

**Suite.** On the box in `~/venv`, `morpheus/v2/tests`: **467 passed, 1 failed, 27 errors** at this
run's HEAD, against **443 passed, 1 failed, 27 errors** at the parent commit `ff9a6f9` on the same
workspace machinery. The delta is exactly the **24** new tests (25 at the final commit). The 27 errors
are `test_p2_figures.py`, which needs matplotlib, absent from that venv by policy. **The 1 failure is
pre-existing and unrelated**: `test_inductive_adjustment.py::test_one_row_at_a_time_equals_the_whole_block`
fails identically at `ff9a6f9` on this box — an `np.allclose(..., atol=0, rtol=0)` exact-equality
assertion that this BLAS build misses by a float ULP. Reproduced at the parent commit in `~/ws_nlabase`
and flagged, not fixed, since it is not this run's file. Locally (numpy 2.4.3 / sklearn 1.7.0) it
passes at both commits, so it is an environment-sensitive exact-equality assertion.

**Files.** `v2/calibra/nonlinear_adjustment.py`, `v2/tests/test_nonlinear_adjustment.py` — commits
`52154ef`, `1ab5586`, `efee0f8`, `73caac7`, `b7c189c` on `research/rebase-vision`. Predeclaration
`010d9c4`. **No statistic is computed inline anywhere**: `top_canonical_correlation`,
`heldout_top_cca`, `effective_rank`, `cross_fitted_residuals`, `confound_design`,
`pooled_tissue_source_site`, `probe_state`, `knn_balanced_accuracy_oof`, `global_permutations`,
`null_summary`, `nesting_diagnostic`, `balanced_accuracy` and `calibration._map` are all imported, and
the two derived quantities of this entry are `retention_of_excess` and `corrected_multiple` in the
module rather than arithmetic in this file.

**Outputs**, under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/nonlinear_adjustment/`:
`gate_reproduction.json`, `channel_main.json` (§1, §3, §5), `artefact_control.json` (§6 shuffle
control), `regenerated_null.json` (§6 corrected null), `kernel_grid.json` (§3, §4, §6),
`forest_arm.json`, `smoke.json`; logs in `logs/`.

---

### 12. Honest constraints on every number above

* **One cohort, one partition, one pair of artifacts.** `wsi_biology` only; `full_biology` and
  `rna_biology` are RNA-derived and their channel to RNA targets is near-circular at ~0.89, so they
  are not a morphology→molecular measurement and were not run.
* **The labels-only ceiling bounds functions of the labels, not higher-moment residual structure.**
  §5 states the gap and the two arguments that narrow it; neither is a proof.
* **The per-cell whitening that would test conditional covariance was not run**, for the estimability
  reason in §7. "The residual is not a conditional mean and not a per-axis conditional variance" is
  therefore where this run stops.
* **`min_site_count` was left at the project default of 10** and no sensitivity sweep was run.
* **Gradient boosting was declared out of scope in advance** (single-output, so every permutation
  would refit per column); kernel ridge and the random forest are the two families run.
* **The corrected null runs at 200 permutations**, p floor 0.0050, against the channel's 2,000. A
  larger count would sharpen the null median, not the direction.
