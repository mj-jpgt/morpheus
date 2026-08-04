## 2026-08-04 23:45 UTC — P1's channel is **not** a transductive artefact: fitted on a separate discovery fold and applied out of sample, `d2_h::wsi_biology` reads 0.6145 against the matched control's 0.6173 on the identical 1,382 patients, **retention of excess 0.997**, and the labels-only ceiling does not rise on any of three label encodings. What does **not** survive is the ceiling's *calibration*: the published **6.0%** is a property of the labels block and the adjustment design being the same columns, and it becomes **19.2%** — in the **transductive** arm — the moment they are not

**Logged:** 2026-08-04 23:45 UTC. **Predeclared in**
`NOTEBOOK_ENTRIES/PREDECLARED_inductive_channel_and_ceiling_20260804T2315Z.md`, committed (`385b049`)
before the wiring was written and before any number below existed.

**How obtained.** Workspace `~/ws_p1ind/morpheus` on the A100 (`150.136.45.194`), built from
`git -c core.autocrlf=false archive HEAD` at commit `ce8f582` and verified **708/708 tracked files**
against an LF-normalised md5 manifest generated from that commit: **0 code files differ, 0 missing,
0 extra** (the 20 files reported differing are `.pdf`/`.png` figures, whose manifest digest is an
artefact of LF-normalising binary content; the `phase_d --workspace-manifest` rule checks
`v2/ tests/ src/ configs/` × `.py|.json|.yaml`, of which **0** differ). CPU only,
`OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1`, `--n-jobs 8` on a 30-core box carrying a co-tenant load of
6–12 throughout; the GPU was not touched and was not needed. Python 3.10.12, numpy 2.2.6,
scikit-learn 1.7.2, scipy 1.15.3.

Driver `v2/research/rebase/nature/p1_evidence/inductive_channel.py`, which **defines no statistic**.
Operator `v2/calibra/inductive_adjustment.ConfoundAdjustmentOperator`, used unchanged — no second
operator was written. Split and state from `p4_certify.exposure_split` / `prepare_state`, imported
unchanged. `claim_guards.py`, `claim_evidence.json` and `paper/P1_CALIBRA_DRAFT.md` are **not
edited**.

Artifacts `runs/d2_final/artifacts/d2_{h,i}_seed42.npz`, state `wsi_biology`, `test` partition;
targets `data/frozen_rna_targets.npz`, the same **90** non-control targets.

---

## Bad news first

**1. The predeclared ceiling bar was not attainable by the transductive control either, and the
reason is a coincidence in the published construction that nobody had named.** §3 of the
predeclaration required the labels-only ceiling under inductive adjustment to sit below the channel's
own null median and below 20% of the channel's excess, calibrated on the published 6.0% / 11.2%. On
the exposure fold the **matched transductive control** reads **19.2%** (saturated cell design) and
**35.8%** (additive design) — worse than the inductive arm on both. The bar therefore does not
discriminate what it was written to discriminate, and it is reported as a defect of the bar rather
than converted into a pass, exactly as `nonlinear_adjustment_channel_result_20260804T2130Z.md` §6 did
with its own step-3 bar.

**The mechanism, measured rather than argued.** The published ceiling used the **108-column
`confound_design`** as the labels "representation" *and* as the design the adjustment residualises
against — the same columns on both sides. That coincidence cannot survive a split: `min_site_count =
10` keeps 84 sites on 2,766 rows, **57 design columns** on the 1,382 exposure rows and **55** on the
1,384 discovery rows. A partition-wide labels block therefore carries ~50 columns that lie outside
**both** arms' adjustment spans and pass through untouched. Restore the coincidence — labels block =
the operator's **own frozen 55-column design**, evaluated on every row — and the ceiling collapses in
every arm:

| labels encoding (`d2_h`) | n = 2,766 transductive | n = 1,382 transductive | n = 1,382 **inductive** |
|---|---:|---:|---:|
| additive `confound_design`, 108 cols | **11.2%** | 35.8% | **27.1%** |
| saturated cell design, 105 cols | **6.0%** | 19.2% | **18.9%** |
| **operator's own frozen design, 55 cols** | **−0.3%** (p = 0.605) | **6.4%** (p = 0.086) | **−19.3%** (p = 1.000) |

*(share of the channel's excess over its own within-cancer pairing null, from
`nonlinear_adjustment.retention_of_excess`.)*

**So 6.0% is not a property of "confound labels versus this channel". It is a property of a labels
block the adjustment design happens to span.** On the encoding that the adjustment actually spans, a
confound-only representation carries **no measurable channel at all** after adjustment — its excess is
negative or insignificant in all three arms. On an encoding it does not span, the ceiling is 3–6×
larger. Both statements are true of the transductive arm and of the inductive arm alike.

**2. The within-cancer and the global pairing null, which agree to 0.003 in every published arm,
disagree by 46% in the inductive arm.** `d2_h`, n = 1,382: within-cancer null median 0.2067 (excess
0.4078) against **global** 0.3028 (excess 0.3117). Retention against the matched control is **0.997**
on the project's within-cancer convention and **0.757** on the global one — the latter below the
0.8 bar §4 set. §5 explains why the global null is the inflated one here and not the informative one,
and the number is reported either way rather than left out.

**3. Neither of those changes the answer to the question this run was commissioned to ask, and the
answer is that the channel survives.** Every comparison at matched n and matched encoding goes the
same way, on both artifacts.

---

## 1. Yes, the flagship result was computed transductively — verified against the code path

`nonlinear_adjustment.channel_under_adjustment` computes `x_adjusted = adjust(x)`,
`y_adjusted = adjust(y)` and `adjust(y[order])` per permutation, and every adjuster
`make_adjuster` builds — `ridge` (which **is** `residualise.cross_fitted_residuals`), `saturated`,
`kernel_ridge`, `forest`, `location_scale`, both `in_sample_*` — fits its nuisance model inside the
call on the rows handed to it. `labels_only_ceiling` routes through the same function. There was no
inductive path in that module. **So all eleven adjusted arms of the 21:30 entry, both legs of its
argument, and P1 §4.4's row are transductive numbers.** Confirmed by reading the code, not inferred
from the general finding.

## 2. This is P4's state, bit for bit

| check | result |
|---|---|
| image operator `reference_digest` | `2060a635fa83756a1c3b7aa8506b7b19fcc4431f5d1a303da39b3cb2bf9d62ce` — **equals** the digest recorded in `p4_inductive_adjustment_measured_20260804T2300Z.md` §8 |
| adjusted image block vs `p4_certify.prepare_state`'s | **bit-identical** (`np.array_equal`) |
| adjusted target block vs `prepare_state`'s | **bit-identical** |
| exposure patient ids vs `prepare_state`'s | **identical** |
| patients in both discovery and exposure folds | **0** (asserted in code) |
| fold sizes / design widths | 1,384 discovery / 1,382 exposure; operator design 55 cols, 31 frequent sites; exposure design 57 cols |

**So the channel numbers below are measured on the very state whose joint site LDA reads 0.2643
against 0.0109 in sample — a factor of 24.** The two entries are about the same 1,382 patients and the
same fitted operator.

## 3. Both reproduction gates, exact

Nothing downstream would be comparable to §4.4 without these, and they were run first.

| quantity | published | this run | |
|---|---:|---:|:---:|
| `d2_h::wsi_biology` adjusted top-CCA (S1) | 0.6052 | **0.6052** | ✔ |
| its within-cancer pairing null median | 0.1483 | **0.1483** | ✔ |
| excess over null median | 0.4569 | **0.4569** | ✔ |
| S2 held-out | 0.5841 | **0.5841** | ✔ |
| `d2_i::wsi_biology` S1 / null / excess | 0.4703 / 0.1472 / 0.3231 | **0.4703 / 0.1472 / 0.3231** | ✔ |
| labels-only ceiling, additive design: raw S1 → adjusted S1, its null, excess | 0.9273 → 0.1237, 0.0723, 0.0514 | **0.9273 → 0.1237, 0.0723, 0.0514** | ✔ |
| labels-only ceiling, saturated cell design | 0.7722 → 0.0903, 0.0631, 0.0272 | **0.7722 → 0.0903, 0.0631, 0.0273** | ✔ |

---

## 4. The channel, inductively adjusted

`wsi_biology`, `test` partition, S1 at 16 components, within-cancer **pairing** null at 2,000
permutations (p floor 1/2001 = 0.0005). Retention is `nonlinear_adjustment.retention_of_excess`, a
ratio of excess over each arm's own null, against the **matched transductive control** — never
against a published n = 2,766 number.

**`d2_h_seed42`**

| arm | n | exposed state | S1 | own null median | excess | **retention** | S2 held-out | eff. rank | p |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| `transductive_full` (gate) | 2,766 | in-sample residual | 0.6052 | 0.1483 | 0.4569 | — | 0.5841 | 22.50 | 0.0005 |
| `none_exposure` | 1,382 | column centring | 0.8234 | 0.7262 | 0.0972 | 0.238 | 0.7958 | 18.89 | 0.0005 |
| **`transductive_exposure`** (control) | 1,382 | in-sample residual | **0.6173** | **0.2080** | **0.4092** | **1.000** | 0.5805 | 21.76 | 0.0005 |
| **`inductive_exposure`** | 1,382 | **out-of-sample operator** | **0.6145** | **0.2067** | **0.4078** | **0.9966** | 0.5701 | 22.04 | 0.0005 |

**`d2_i_seed42`**

| arm | n | S1 | own null median | excess | **retention** | S2 | eff. rank | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `transductive_full` (gate) | 2,766 | 0.4703 | 0.1472 | 0.3231 | — | 0.4206 | 11.73 | 0.0005 |
| `none_exposure` | 1,382 | 0.7877 | 0.6829 | 0.1048 | 0.354 | 0.7682 | 10.63 | 0.0005 |
| **`transductive_exposure`** | 1,382 | 0.5020 | 0.2056 | 0.2963 | 1.000 | 0.4215 | 11.56 | 0.0005 |
| **`inductive_exposure`** | 1,382 | **0.4898** | **0.2021** | **0.2877** | **0.9710** | 0.4170 | 11.56 | 0.0005 |

**Retention 0.997 (`d2_h`) and 0.971 (`d2_i`).** Every `permutation_p` is at the 0.0005 floor: no
permutation of two thousand reached the observed value under any adjustment, inductive or otherwise.
S2 — directions fitted on one half of the exposure fold and scored on the other, so immune to the
in-sample maximisation S1 performs — moves the same way and by the same little: 0.5805 → 0.5701 on
`d2_h`, 0.4215 → 0.4170 on `d2_i`.

**The adjustment is not a no-op, and it is not the transductive one relabelled.** Per-axis
raw-vs-adjusted correlation median **0.7536** (inductive) against 0.7524 (transductive), residual
variance ratio **0.5896** against 0.5931, **0 of 256** axes above 0.99 in either arm — i.e. both
remove ~41% of the median axis's variance. Between the two adjusted blocks, per-axis correlation
median **0.9598** (min 0.8978) and relative Frobenius difference **0.285**;
`adjuster_agreement.is_relabelled_incumbent` is **False**. This reproduces the P4 entry's §6.1
finding from the channel side: the inductive adjustment removes *as much* as the transductive one and
removes *a different part* — and the channel does not notice.

## 5. Which null, and why the global one is the inflated one here

The 21:30 entry computed both pairing nulls for every arm and found them agreeing to ~0.003, so no
conclusion turned on the choice. That agreement is a property of a **transductive** adjuster, and it
breaks here:

| arm (`d2_h`, n = 1,382) | within-cancer null median | global null median | within-cancer excess | global excess |
|---|---:|---:|---:|---:|
| `none_exposure` | 0.7262 | 0.1980 | 0.0972 | 0.6254 |
| `transductive_exposure` | 0.2080 | 0.2056 | 0.4092 | 0.4117 |
| **`inductive_exposure`** | **0.2067** | **0.3028** | **0.4078** | **0.3117** |

Retention against the matched control is **0.997** within cancer and **0.757** globally.

**The mechanism, and it is the same one P1 §4.6 already owns.** An inductive adjustment is a *fixed*
map: `y_adj[i] = y[i] − f_y(d_i) − c`, with `d_i` the exposure row's own design. Permuting `y`
underneath it leaves the `− f_y(d_i)` term keyed to position `i`, and `x_adj` carries `− f_x(d_i)`,
also keyed to `i` — a deterministic function of the same design, shared by the two blocks. That is
induced correlation, and it lands in the null. Under a **within-cancer** permutation `y[order[i]]`
comes from the same cancer, so its own conditional mean largely cancels the cancer part of
`f_y(d_i)` and what is regenerated is the site-scale residual — which is what should be regenerated.
Under a **global** permutation `y[order[i]]` carries an unrelated cancer's mean, so `− f_y(d_i)`
survives at full cancer scale and correlates strongly with `− f_x(d_i)`; the null is measuring the
map's induced correlation at a magnitude the observed data does not have. The `none_exposure` row
shows the two nulls routinely differ by far more than this in the other direction (0.7262 vs 0.1980),
so "they agree" was never the general case.

**The within-cancer null is the project's convention, is the one P1 §4.4 quotes, and is the one
defended here. The global figure is reported beside it, and a reader who prefers it should read
retention 0.757 rather than 0.997 — still not a collapse, and still above the `none` arm's 0.238.**

## 6. The labels-only ceiling — the leg the verdict actually rests on

The question is not whether confound information decreased. It is whether the channel stays above
what confound labels alone can produce **through the same pipeline**. Under the inductive adjustment
the un-removed confound in each block is `m(d) − f̂_D(d)` — the true conditional mean minus the
discovery fold's estimate of it — which is **itself a function of the labels**. So the labels-only
ceiling is a *tighter* bound out of sample than in sample, where the published entry had to concede
(§5, "the honest limit of this bound") that the surviving residual was higher-moment and not a
function of the labels in the mean.

`d2_h`, n = 1,382 exposure fold, same 2,000-permutation within-cancer pairing null, same 16
components, the same inductively adjusted targets on both sides.

| labels encoding | arm | raw S1 | **adjusted S1** | its own null median | excess | **% of the channel's excess** | adjusted S2 | p |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| additive design, 108 cols (rank 105) | transductive, n = 2,766 | 0.9273 | 0.1237 | 0.0723 | 0.0514 | **11.2%** | 0.0050 | 0.0010 |
| additive design | transductive, n = 1,382 | 0.9309 | 0.3540 | 0.2073 | 0.1467 | **35.8%** | 0.1240 | 0.0005 |
| additive design | **inductive, n = 1,382** | 0.9309 | **0.3096** | 0.1991 | 0.1105 | **27.1%** | 0.1051 | 0.0005 |
| saturated cell design, 105 cols | transductive, n = 2,766 | 0.7722 | 0.0903 | 0.0631 | 0.0273 | **6.0%** | 0.0811 | 0.0080 |
| saturated cell design | transductive, n = 1,382 | 0.7316 | 0.2828 | 0.2042 | 0.0786 | **19.2%** | 0.1167 | 0.0025 |
| saturated cell design | **inductive, n = 1,382** | 0.7316 | **0.2780** | 0.2009 | 0.0771 | **18.9%** | 0.1338 | 0.0025 |
| **operator's own frozen design, 55 cols (rank 52)** | transductive, n = 2,766 | 0.9277 | 0.0599 | 0.0615 | −0.0016 | **−0.3%** | 0.0968 | **0.6047** |
| **frozen design** | transductive, n = 1,382 | 0.9304 | 0.1970 | 0.1706 | 0.0264 | **6.4%** | 0.0187 | **0.0860** |
| **frozen design** | **inductive, n = 1,382** | 0.9304 | **0.2599** | 0.3386 | **−0.0787** | **−19.3%** | 0.1346 | **1.0000** |

**`d2_i`.** The ceiling arms do not depend on the artifact — the labels blocks and the targets are
properties of the cohort — and duly reproduce **to every digit**. Against `d2_i`'s smaller channel
excess the shares are 49.5% / 38.4% (additive, transductive / inductive), 26.5% / 26.8% (saturated),
8.9% / −27.4% (frozen). Same ordering, same conclusion.

**Three things follow.**

1. **The ceiling does not rise under inductive adjustment on any encoding.** 27.1% against 35.8%,
   18.9% against 19.2%, −19.3% against 6.4%. The mechanism §4 of the predeclaration named — site
   information leaking back at 24× and dragging the confound-only representation up with it — **did
   not happen**. It was the outcome I put 0.7 on, and it is wrong.
2. **On the encoding the adjustment can span, a confound-only representation carries nothing.**
   Frozen design: p = 0.605 at n = 2,766, p = 0.086 at n = 1,382 transductive, p = 1.000 inductive
   (observed *below* its own null). That is the strongest form of P1 §5's argument, and it is the
   first time it has been shown out of sample.
3. **The absolute 6.0% / 11.2% figures are not stable, and the instability is not inductive.** It is
   the labels-block-vs-adjustment-design span mismatch of §"Bad news first" item 1, and it moves the
   transductive control by more than it moves the inductive arm.

**Raw, the labels still beat the image**, as published: 0.9309 (additive design) against the image
block's 0.8234 on the exposure fold. Confounding on this cohort is not a subtle worry before
adjustment, which is the whole motivation for applying one.

---

## 7. Every predeclared check, discharged

### §3 — what would make me say the claim SURVIVES

| | predeclared | measured | |
|---|---|---|:---:|
| S1 | `permutation_p` at the 0.0005 floor | **0.0005** on both artifacts | ✓ |
| S1 | retention ≥ 0.8 vs matched control | **0.9966** (`d2_h`), **0.9710** (`d2_i`) within cancer; **0.757** (`d2_h`) global | ✓ / see §5 |
| S2 | ceiling below the channel's own null median | 0.3096 / 0.2780 / 0.2599 against 0.2067 — **above** on all three | ✗ |
| S2 | ceiling excess ≤ 20% of the channel's excess | 27.1% ✗, 18.9% ✓, −19.3% ✓ | mixed |

### §4 — what would make me say it needs the same NARROWING

* **N1** fires by the letter on the additive encoding (27.1% > 20%) and on the null-median leg for
  all three. **It fires on the matched transductive control by more** (35.8%, and two of three above
  the null median). A bar that the in-sample control also fails is not evidence about out-of-sample
  adjustment; it is evidence that the bar was calibrated on a construction that the split destroys.
  **Reported as a defect of the bar. N1 is NOT treated as fired against the inductive arm.**
* **N2** does not fire on the project's null convention (0.997, 0.971) and does fire on the global one
  (0.757), for the reason and with the mechanism in §5.

### §5 — what would make me distrust this FAVOURABLE result

1. **The inductive null failing to regenerate what the transductive one does, inflating retention.**
   Predeclared: if the inductive null median falls >30% below the control's, retention is
   uninterpretable. **It did not fall — it fell 0.6% (0.2067 vs 0.2080) within cancer and *rose* 47%
   globally.** My predicted direction was wrong; the check that was written to catch a flattering
   result would have caught it, and there was nothing to catch. Retention is interpretable.
2. **The adjustment doing nothing.** 0 of 256 axes above 0.99, median raw-vs-adjusted correlation
   0.7536, variance ratio 0.5896, `is_relabelled_incumbent` False. Does not fire.
3. **The reproduction gates.** Four of four exact, §3. Does not fire.
4. **Provenance.** Digest match, bit-for-bit agreement with `prepare_state`, 0 fold overlap. Does not
   fire.
5. **The two artifacts disagreeing.** They do not: retention 0.997 and 0.971, ceiling ordering
   identical.
6. **Capacity / rank collapse.** Adjusted `effective_rank` 22.04 (inductive) vs 21.76 (transductive)
   vs 22.50 (published full), against a 16-component budget. It *rises* under the inductive
   adjustment rather than collapsing. Does not fire.

### §6 — what would make me distrust an UNFAVOURABLE result

1. Over-removal guard: **not triggered** — no arm collapsed (lowest retention 0.971), so the
   `spike_recovery_curve` attenuation measurement was not run. Stated so its absence is a consequence
   of the result and not an omission.
2. The ceiling is graded on excess over its own null **and** on raw S1 against the channel's null
   median; both are in §6's table.
3. Every exposure arm sits at n = 1,382, and **no inductive number in this entry is compared to a
   published n = 2,766 number without the matched control between them.** The matched control is what
   caught the ceiling's instability.

---

## 8. What this changes

**For P1, and this is the load-bearing sentence.** The 23:00 P4 entry narrowed what may be said about
an adjusted state's *confound content* — "the confound is removed from the first moment of the rows
the nuisance model was fitted on". That narrowing is correct and stands. **It does not propagate to
the channel claim.** Measured on the identical patients and the identical operator: the
morphology→molecular channel is 0.997 of its matched in-sample value when the nuisance model is fitted
on a separate discovery fold and applied out of sample, at the permutation floor, with the
labels-only ceiling no higher than in sample on any of three encodings and *zero* on the encoding the
adjustment spans. **P1 §4.4 does not need the narrowing the site-recoverability claim just got, and
this run is a strengthening of it, not a qualification.**

**One thing P1 §5 does need.** The 6.0% / 11.2% ceiling figures are properties of a labels block that
the adjustment design spans, not of the cohort. Quoted without that condition they are not
reproducible on any sub-cohort whose pooling keeps fewer sites. §5's existing sentence — *"the
additive figure is the larger only because residualising a design on itself leaves ridge-shrinkage
remainder at α = 1"* — already contains the observation in embryo; what it does not say is that the
number moves 3–6× when the labels block leaves the design's span.

### Prose flagged, and deliberately NOT edited (multiple agents are correcting these files)

1. **`paper/P1_CALIBRA_DRAFT.md` §4.4 / §5 and
   `NOTEBOOK_ENTRIES/nonlinear_adjustment_channel_result_20260804T2130Z.md` §0 item 2, §5.** Wherever
   the labels-only ceiling is quoted as **6.0%** / **11.2%**, it should say *measured with the labels
   block encoded in the same columns the adjustment residualises against*. Measured here: the same
   quantity at n = 1,382 with a partition-wide labels block reads **19.2%** / **35.8%** under the
   **identical transductive** adjustment, and **−0.3%** at n = 2,766 / **6.4%** at n = 1,382 with the
   frozen 55-column encoding.
2. **Same entry, §0 and §3, "the channel survives entirely".** True, and now demonstrably *not*
   restricted to a transductive adjustment: it survives an out-of-sample one at retention 0.997
   (`d2_h`) / 0.971 (`d2_i`). This is an **addition** available to §4.4, not a correction.
3. **`nonlinear_adjustment_channel_result_20260804T2130Z.md` §2**, *"the global pairing null was also
   computed at 2,000 permutations for every arm and the two agree to ~0.003 … so no conclusion here
   turns on the choice"*. Correct for every arm in that entry and **not general**: under an inductive
   adjuster the two disagree by 46% (0.2067 vs 0.3028), and under no adjustment at all by 267%
   (0.7262 vs 0.1980). The sentence should be scoped to transductively adjusted arms.
4. **`nonlinear_adjustment.py` module docstring** does not say that every adjuster it builds is
   transductive. It now has an inductive caller and should say which is which; the `adjust_y`
   docstring added in this run states the consequence for the null.

**For P4.** Nothing here rescues P4's condition 3 — the exposed adjusted state still does not certify
(joint site LDA 0.2643 against a null p95 of 0.1495). This entry says only that the *channel* on that
state is intact. **A state can carry recoverable site information and still carry a morphology→
molecular channel that confound labels cannot reproduce; those are different claims about the same
matrix and this run separates them for the first time.**

---

## 9. How the predictions did

| predicted (§4 of the predeclaration) | measured | |
|---|---|:---:|
| N1 fires — the ceiling rises under inductive adjustment (p = 0.7) | **it falls** on all three encodings | ✗ |
| N2 fires — retention < 0.8 (p = 0.35) | 0.997 / 0.971 within cancer; 0.757 global | ✗ / partly |
| the inductive null median sits *lower* than the control's | 0.6% lower within cancer; **47% higher** globally | ✗ |
| the adjustment audit shows a real adjustment | 0.7536 corr, 0.5896 variance ratio, 0 axes > 0.99 | ✓ |
| the reproduction gates hold | 4/4 exact | ✓ |

**The one prediction the conclusion rests on — that the ceiling would rise out of sample — was wrong,
and it was the prediction I was most confident in.** The mechanism I proposed for it (56.2% of
exposure rows get no site adjustment, so the labels pass through) is real and is exactly what makes
the *site certificate* fail at 24×; it simply does not translate into cross-block channel, because a
labels-only representation pushed through the same pipeline does not gain either.

---

## 10. Suite

Run on this workspace at the commit the runs were launched from (`ce8f582`):
`pytest morpheus/v2/tests morpheus/tests --ignore=morpheus/v2/tests/test_p2_figures.py -q` →
**610 passed, 0 failed in 65.25 s**. `test_p2_figures.py` run separately reads **1 passed, 27 errors
in 2.53 s**, every error `ModuleNotFoundError: No module named 'matplotlib'` — the known condition of
`~/venv`. **Nothing was installed into that environment.**

The 10 new tests are this run's whole test delta: `test_inductive_channel.py` alone reads **10
passed**. `test_nonlinear_adjustment.py`, whose module this run modified, still reads **25 passed**
unchanged, and `test_p4_inductive_wiring.py` **15 passed**.

## 11. Files / provenance

Driver `v2/research/rebase/nature/p1_evidence/inductive_channel.py` (commits `b7aadae`, `ce8f582`).
Tests `v2/tests/test_inductive_channel.py`, 10 tests. Library change:
`nonlinear_adjustment.channel_under_adjustment` / `labels_only_ceiling` gained an optional
`adjust_y`, default-preserving, pinned by `test_adjust_y_defaults_to_adjust_byte_for_byte`.

Outputs `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/inductive_channel/{d2_h_final.json,
d2_i_final.json}` plus `logs/`, vendored into `v2/research/rebase/nature/p1_evidence/out/`.
Operator digests: image `2060a635…d9c2ce` (equals P4's), targets `d7a37321…d847a7`.

**No statistic is computed inline anywhere in the driver.** `top_canonical_correlation`,
`heldout_top_cca`, `effective_rank`, `cross_fitted_residuals`, `confound_design`,
`pooled_tissue_source_site`, `cell_design`, `channel_under_adjustment`, `labels_only_ceiling`,
`retention_of_excess`, `adjuster_agreement`, `cross_fitted_r2`, `_load_block`, `exposure_split`,
`prepare_state` and `_adjustment_audit` are all imported unchanged.

## 12. Honest constraints on every number above

* **One split.** Discovery fraction 0.5 at seed 42 only — the one P4 measured on, chosen for
  comparability. The 0.7 sensitivity arm P4 ran was **not** repeated here.
* **One cohort, one partition, `wsi_biology` only.** `full_biology` and `rna_biology` are RNA-derived
  and near-circular at ~0.89, so they are not a morphology→molecular measurement.
* **`min_site_count` left at the project default of 10**, no sensitivity sweep — and §"Bad news
  first" shows the ceiling *is* sensitive to what that parameter does to the design width, so a sweep
  is the obvious next experiment.
* **The ceiling bounds what a representation that is a function of the labels can contribute.** Out
  of sample that is a tighter bound than in sample (§6), but it is still a bound on first-moment
  label structure and not on conditional covariance; the per-cell whitening the 21:30 entry declared
  unestimable is still unrun and still unestimable at 105 cells × 256 axes.
* **The confound certificate and the detection floors were not re-run**; P4's 23:00 entry carries
  them on this exact state and this run does not restate them.
* **`n = 1,382` throughout the exposure arms.** Every inductive claim is against the matched
  transductive control at the same n, never against a published n = 2,766 number.
