# 2026-08-04 18:30 UTC — The centring-amplification law is exact in its own premise and the premise is false on every view we have. It does not explain the `wsi`/`rna` floor gap, and on two views centring *improves* reproducibility

**Logged:** 2026-08-04 18:30 UTC. **Predeclared** at
`NOTEBOOK_ENTRIES/PREDECLARED_centring_amplification_law_20260804T1750Z.md`, commit **`d4e344c`**,
written and pushed **before the harness existed and before any number was taken**. **How obtained:**
CPU only, thread-capped (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1`), four background processes on `150.136.45.194` from `~/ws_amp`, a workspace
built with `git -c core.autocrlf=false archive HEAD` and verified **618/618 files by git blob SHA-1**
before anything ran. The A100 was at 99% on another chain throughout and was never touched. Scored by
`v2/research/rebase/p2/p2_centring_amplification.py`, graded by
`v2/research/rebase/p2/p2_centring_verdict.py`, outputs vendored at
`v2/research/rebase/p2/figures/data/ws_amp/out/`.

---

## 0. Headline, bad news first

**The candidate law is not a coincidence and it is also not usable here.** Its functional form is
derivable exactly, it is confirmed to a median relative error of **0.05%** across `f ∈ [0, 0.99]` in
the one synthetic condition where its premise holds — and **its premise is false on all three of our
views, by factors of 4×, 8× and 17×.** Consequently:

* **It does not explain the observation that motivated it.** `1/(1−f)` predicts **132×** for the
  exported `wsi_biology` block; the correctly derived form predicts **4.92×**; the measured
  amplification is **1.93×**. On `full_biology` the derived form predicts 1.16× and the measured value
  is **3.16×** — wrong by 173% and in the wrong direction relative to the other views. **The three
  views' predicted ordering is `wsi > rna > full`; the observed ordering is `full > wsi > rna`.**
  Predeclared falsifier **P2: FALSIFIED**, twice over.
* **The two observations that "line up suspiciously well" are two different phenomena.** Observation 1
  (RankMe 1.811× vs centred R1 3.111×) is a centring effect and the derivation gets its *sign and rough
  size* right. Observation 2 (3.295× on `wsi_biology` vs 1.019× on `rna_biology`) is **52.6% already
  present before anything is centred at all**, and it is measured on the *residualised* block where
  **`f = 0` to 20 decimal places** and the law therefore predicts an amplification of exactly 1. The
  law says nothing about observation 2.
* **Centring measurably IMPROVES reproducibility in nine of our own view × statistic cells**, robustly
  to dropping any single repeat — R1 on `rna_biology` at **0.785** (leave-one-out range 0.750–0.833),
  down to 0.328 for R2 on `full_biology`. Predeclared falsifier **P4: FALSIFIED**. The derivation
  predicted exactly this regime (§1.9 of the predeclaration) and it is where we live.
* **The one-line quantity the framing quotes is wrong for this block.** `‖mean_i z_i‖² = 0.8133` is
  from the 16-patient gate batch at initialisation. On the block the floors are actually measured on,
  **`f = 0.9917–0.9950`**, and — for L2-normalised rows, which these are — **`f` and the
  patient-to-patient cosine are the same number** (0.9917 against 0.9917 on rep1). Observation 3 is not
  independent evidence for observation 2; it is the same measurement.

**What survives, and it is worth keeping.** The transfer identity is exact, the naive `1/(1−f)` is
badly wrong in a way that is now derived rather than fitted, and the failure is *localised to a named,
measurable assumption* with a usable sensitivity criterion attached (§5). That is a smaller result than
a law, and it is a real one.

---

## 1. The derivation, and where `1/(1−f)` goes wrong

Full statement and assumptions in the predeclaration §1. In brief, with `Z = 1 z̄ᵀ + Z_c`,
`s_• = √n‖z̄‖`, and `t = s_• / (s_• + ‖Z_c‖_*)` the **L1 / nuclear-norm share** of the mean term:

> **Theorem 1.** If `z̄ ⊥ row(Z_c)` then `spec(Z) = {s_•} ∪ spec(Z_c)` exactly, and by the Shannon
> chain rule
> **`ln D₁(uncentred) = h(t) + (1 − t) · ln D₁(centred)`**, `h` the binary entropy.
>
> **Theorem 2.** `δ ln D_a(uncentred) = (1 − w_a) δ ln D_a(centred) + c_a δt`, with
> `w_a = tᵃ / (tᵃ + (1−t)ᵃ D_a(centred)^{1−a})`, `w₁ = t`. **Under `Var(δt) ≈ 0` (assumption A2) the
> amplification is `A_a = 1/(1 − w_a)`, and `A₁ = 1/(1 − t)`.**
>
> **Corollary.** `t = √f / (√f + √R2_centred · √(1−f))` is an **identity**, not an approximation
> (`‖Z‖_F² = s_•² + ‖Z_c‖_F²` and `R2_centred = ‖Z_c‖_*²/‖Z_c‖_F²`), so
> **`A₁ = 1 + √(f/(1−f)) / √R2_centred`.**

**Two structural disagreements with `1/(1−f)`, both now measured.**

1. **The divergence is a square root, not a pole.** `A₁ ~ (1−f)^{−1/2}`. At `f = 0.990` in the
   synthetic sweep the derived form gives **2.106** and the measured value is **2.107**; `1/(1−f)`
   gives **101.2**, wrong by **48×**.
2. **The constant is the residual spectrum's own order-2 statistic**, `√R2_centred`. A dominant
   component sitting on a rich residual costs far less than the same component on a thin one at
   identical `f`. `1/(1−f)` has no room for that term and it is worth a factor of 3–5 on our data.

**A defect in the derivation's own gloss, found by running it.** Predeclaration §1.8 asserted an
ordering `A₂ ≥ A₁ ≥ A_∞ ≥ A₀`. That is **wrong for `a → ∞`**: `D_∞ = 1/max_k p_k`, so once the spike's
share exceeds the residual's largest normalised eigenvalue the uncentred statistic is exactly `1/f` and
**stops responding to the residual at all** — the amplification is unbounded, not small. The general
formula `A_a = 1/(1−w_a)` already contains this (the branch returns "infinite"); only the informal
sentence was wrong. **S4 is therefore FALSIFIED in all three synthetic conditions, including the one in
which the derivation is exact to 0.05%.** It is recorded as falsified rather than repaired.

---

## 2. Synthetic sweep — the form is exact where its premise holds, and inverts where it does not

`n = 2766`, `d = 256`, 25 runs per point, 20 values of `f` from 0 to 0.99, three conditions. The bulk's
directional gains are re-drawn per run so that the centred statistic's fold lands at **1.03–1.06×**,
inside the range our four concordant repeats occupy; the first version of the sweep perturbed a fixed
bulk additively, which self-averages at `n = 2766` to `sd(ln D₁) ~ 1e-4` and made every ratio a ratio
of two numbers that were not moving. That is recorded in the module, not worked around.

| condition | what it does to `t` | S1 (does `A_obs = 1/(1−t)`?) | S2 (does it beat `1/(1−f)`?) |
|---|---|---|---|
| **`t_pinned`** *(post-hoc, see below)* | pinned exactly | **not falsified** — median rel. err **0.0005**, worst **0.0020** | **not falsified** — `1/(1−f)` median rel. err **0.4927** |
| **`stable_mean`** *(the predeclared A2 condition)* | `sd(t) ≤ 0.0028` | **FALSIFIED** — median **0.2195** | not falsified |
| **`unstable_mean`** | `sd(t) ≤ 0.0134` | **FALSIFIED** — median **0.5189** | not falsified |

**`t_pinned` is a post-hoc addition and is labelled as one.** The predeclaration described its
stable-mean condition as *"A2 holds by construction"*. **It does not, and that is a defect in the
candidate law's statement of itself, not a bookkeeping slip.** A2 is a claim about
`t = s_•/(s_• + ‖Z_c‖_*)`; holding the shared direction fixed pins only the numerator, and the
residual's own nuclear norm still moves. *Holding the dominant component stable is not sufficient for
the law's premise.* `t_pinned` rescales each run's centred part to a common nuclear norm — which
changes no centred value, every statistic here being scale-invariant on the centred matrix — and is the
only condition in which A2 actually holds.

**The observed amplification, `t_pinned`, `sd(ln·)` ratio, 25 runs per row (every third row):**

| `f` | `t` | **A₁ observed** | **A₁ derived** | `1/(1−f)` | A₂ obs | A₂ derived | A₀ obs | A₀ derived |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.000 | 0.0000 | 1.0000 | 1.0000 | 1.00 | 1.0000 | 1.000 | 1.0000 | 1.004 |
| 0.200 | 0.0524 | 1.0540 | 1.0553 | 1.25 | 1.2515 | 1.249 | 1.0000 | 1.004 |
| 0.497 | 0.0990 | 1.1094 | 1.1099 | 1.99 | 1.9875 | 1.987 | 1.0000 | 1.004 |
| 0.697 | 0.1441 | 1.1690 | 1.1683 | 3.30 | 3.2951 | 3.299 | 1.0000 | 1.004 |
| 0.850 | 0.2094 | 1.2655 | 1.2648 | 6.66 | 6.6422 | 6.660 | 1.0000 | 1.004 |
| 0.921 | 0.2742 | 1.3795 | 1.3778 | 12.62 | 12.6160 | 12.62 | 1.0000 | 1.004 |
| 0.969 | 0.3826 | 1.6208 | 1.6197 | 32.24 | 32.5943 | 32.25 | 1.0000 | 1.004 |
| **0.990** | **0.5251** | **2.1071** | **2.1057** | **101.16** | 96.7805 | 101.2 | 1.0000 | 1.004 |

Median `|relative error|` over all 20 rows: **R1 0.0006, R2 0.0020, hard rank 0.0039.** Three orders of
the same family, over two orders of magnitude in the amplification, with no fitted parameter.

**Where the derivation breaks even under `t_pinned`, reported because it is a boundary and not a
detail.** The eigenvalue-spectrum order-2 statistic (`PR`) is predicted 3.3 → 3.8×10⁵ as `f` goes
0.2 → 0.99 and observed 5.4 → 154 — median relative error **0.97**. The reason is visible in the
formula: for `PR` the mixture weight uses `f` rather than `t`, `w₂ → 1` almost immediately, and the
first-order expansion behind Theorem 2 is invalid once `1 − w_a` is small. **The derived form is
trustworthy on the singular-value family (orders 0, 1, 2) and not on the eigenvalue family above
`f ≈ 0.5`.** Stable rank's *qualitative* prediction — unbounded, the largest of any statistic — is
confirmed (observed 12 → 904) but its magnitude is not predicted at all.

**S3 — the regime where centring helps.** Under `unstable_mean`, `A₁` falls to **0.155** at `f = 0.99`:
**centring makes the statistic 6.4× MORE reproducible.** The mechanism is exactly the one the
derivation names — the uncentred statistic absorbs the spike's own fluctuation through the `c_a δt`
term — and it is not a curiosity, because §3 shows it is the regime our own artifacts are in.

---

## 3. Real data — five same-seed repeats, three views, both blocks

`~/e0_run/d1_envelope/rep{1..5}.npz`, test partition, `n = 2766`, `d = 256`, cancer + pooled-TSS
cross-fitted residualisation at seed 42 — the same block, split and residualisation
`p2_envelope_floors.py` reads. **The `[X; −X]` reduction used to evaluate the three centre-only
statistics uncentred was checked against `R1_uncentred`, `R2_uncentred` and `rankme` on every rep and
view: 15/15 passed at relative error ≤ 1e-9** (observed maximum ~1e-16). The predeclaration made a
failure a stopping condition.

### 3.1 P1 — assumption A1 holds. This is the part that works

The transfer identity's worst relative error over all 30 rep × view × block cells is **0.0077**, and on
the residualised block it is **0.0000** (`f = 0` exactly). **P1: not falsified.** So the mean direction
*does* detach as a clean spike, exactly as the finite-rank separation results describe, and the
derivation's one approximation is not what fails.

### 3.2 The shared direction, measured — and it is not 0.8133

| view | block | `f` (energy share of the mean) | mean patient-to-patient cosine | `t` (nuclear-norm share) |
|---|---|---|---|---|
| `wsi_biology` | raw | **0.9917 / 0.9950 / 0.9919 / 0.9917 / 0.9920** | 0.9917 … 0.9950 | 0.7686 … **0.8618** |
| `rna_biology` | raw | 0.2606 – 0.2775 | 0.2603 – 0.2772 | 0.1371 – 0.1438 |
| `full_biology` | raw | 0.2615 – 0.2828 | 0.2612 – 0.2826 | 0.1346 – 0.1410 |
| all three | residualised | **< 1e-20** | 0.0009 – 0.0025 | **0.0000** |

*Per rep, never a mean; rep2 is the second value in each list.* Two things follow immediately. **`f`
and the mean cosine are the same number to four decimals** — for L2-normalised rows
`f = ‖z̄‖²` and the mean off-diagonal cosine is `(n‖z̄‖² − 1)/(n − 1)`, so the framing's observation 3
is not a second measurement. And **residualisation removes the shared direction completely**, so on the
residualised block — which is where 3.295× and 1.019× are measured — there is nothing left for
centring to amplify and the law's prediction is exactly 1× for every view.

### 3.3 P2 — the law itself. FALSIFIED

Raw block, canonical R1, amplification as the ratio of `sd(ln ·)` across the five repeats:

| view | `f` | `t` | **A observed** | **A derived** `1/(1−t)` | rel. err | `1/(1−f)` |
|---|---:|---:|---:|---:|---:|---:|
| `wsi_biology` | 0.9925 | 0.7882 | **1.928** | 4.923 | −61% | 132.5 |
| `rna_biology` | 0.2688 | 0.1401 | **0.785** | 1.163 | −33% | 1.37 |
| `full_biology` | 0.2703 | 0.1373 | **3.163** | 1.159 | **+173%** | 1.37 |

*The `ln(max/min)` ratio agrees with the `sd(ln·)` ratio to within 7% on every row (1.911 / 0.764 /
2.950), so the verdict is not an artefact of which dispersion is used.*

**Falsified on the magnitude** (`full_biology` at +173%, over the predeclared ±50%) **and falsified on
the ordering** — predicted `wsi > rna > full`, observed `full > wsi > rna`. Either alone was
predeclared as sufficient.

### 3.4 P3b — why. Assumption A2 is false on all three views

Theorem 2's second term is `c₁ δt` with `c₁ = ln((1−t)/t) − ln D₁(centred)`. A2 says it is negligible.

| view | `t` | `sd(t)` | `sd(t)/t` | `\|c₁\|` | **spike term `\|c₁\|·sd(t)`** | **bulk term `(1−t)·sd(ln D₁)`** | ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| `wsi_biology` | 0.7882 | 0.0412 | 5.2% | 4.297 | **0.1769** | 0.1065 | **1.66** |
| `rna_biology` | 0.1401 | 0.0026 | 1.9% | 1.376 | 0.0036 | 0.0087 | **0.42** |
| `full_biology` | 0.1373 | 0.0026 | 1.9% | 1.438 | 0.0038 | 0.0049 | **0.78** |

For the law to hold to 10% the ratio must be ≤ 0.1. It is **17×, 4× and 8× too large.** By the
predeclared threshold (spike term > 50% of the uncentred variance) **P3b: FALSIFIED** on
`full_biology`, where the decomposition is incoherent (predicted variance 4.6× the observed, i.e. the
two terms are cancelling); `wsi_biology` sits at 46%, just under the line, and `rna_biology` at 8%.

**The plain-language version: "a dominant *stable* component" is not a thing our retraining produces.**
Rep2's WSI state is more collapsed onto the shared direction than its siblings — `f = 0.99498` against
0.9917–0.9920 — which moves `t` from 0.768 to 0.862. A 5% wobble in the shared direction's mass share
is enough to invert the sign of the predicted effect, and 5% is small compared with everything else
these five runs do.

### 3.5 P4 — centring improves reproducibility in nine cells. FALSIFIED

Raw block, `A < 1` and still `< 1` after dropping any single repeat:

| view | statistic | **A observed** | leave-one-out range |
|---|---|---:|---|
| `rna_biology` | **R1** | **0.785** | 0.750 – 0.833 |
| `rna_biology` | R2 / R3 | 0.772 | 0.730 – 0.850 |
| `rna_biology` | PR | 0.694 | 0.641 – 0.774 |
| `rna_biology` | stable rank | 0.490 | 0.322 – 0.772 |
| `full_biology` | R2 / R3 | **0.328** | 0.140 – 0.479 |
| `full_biology` | PR | 0.420 | 0.352 – 0.640 |
| `full_biology` | stable rank | 0.404 | 0.314 – 0.671 |

**This is not blind** — the predeclaration records that R1 = 1.0228× against RankMe = 1.0299× on
`rna_biology` raw had already been read, and says so in the falsifier itself. It is recomputed here
against a true `R1_uncentred` rather than against RankMe, and it survives.

### 3.6 P3 — the order-a prediction passes on the letter and its content is negative

Spearman `ρ = 0.367 > 0` over the nine raw-block statistic × view cells, so **P3: not falsified.** But
the within-view correlations are **+1.0 (`wsi`), −1.0 (`rna`), −0.5 (`full`)**: the positive overall
value is entirely between-view, and inside two of three views the prediction is *anti*-correlated with
what happens. The threshold was too weak and I am saying so rather than quoting the `ρ`.

### 3.7 P5 — half the view effect predates centring, and the half that does not is on the wrong block

The three nested preprocessing levels, canonical-R1-family fold across the five repeats:

| view | **nothing removed** (uncentred) | **mean removed** (centred) | **mean + cancer + TSS removed** (residualised) |
|---|---:|---:|---:|
| `wsi_biology` | **1.8112×** | 3.1110× | 3.2947× |
| `rna_biology` | 1.0299× | 1.0228× | 1.0193× |
| `full_biology` | 1.0047× | 1.0138× | 1.0200× |

Cross-view spread on the log scale: **1.1212 centred, 0.5893 uncentred → 52.6% of the view effect is
present before anything is removed.** By the predeclared threshold (>70%) **P5: not falsified**, and
its honest content is: *centring roughly doubles a floor gap that is already 20× before centring.* The
20× is the finding; the doubling is the law's contribution to it.

### 3.8 Four attempts to break it further, as asked

* **Is the hard rank a counterexample?** No — it is a **limit case the derivation gets right.** Order 0
  predicts `A₀ = 256/255 = 1.0039`, i.e. no amplification, and it predicts no variance to amplify when
  the residual is full rank in every run. Measured: **exactly 1.000× in every view and both blocks**,
  and the synthetic sweep returns `A₀ = 1.0000` at every `f`. It is pinned by a test.
* **Are the low-`f` views lower-dimensional, or on fewer effective samples?** **No, and the confound
  runs the other way.** All three views are `n = 2766 × d = 256` on the identical patient set. The
  low-`f` views have *higher* centred effective rank — `rna` 23.98–24.53 and `full` 26.33–26.69 against
  `wsi`'s 8.03–24.99. The tighter floors sit on the *richer* representations.
* **Is `f` at least a usable screen, even if the form is wrong?** On three views, with two distinct
  values of `f`, and measured on the same states as the outcome — so it cannot be separated from "rep2's
  WSI encoder partially collapsed and everything about that view is worse". **Three points is not a
  test and we do not claim one.**
* **Is the eps in RankMe doing any of this?** No. `R1_uncentred` and `rankme` agree to within
  0.02% on every rep and view (e.g. 3.5725 against 3.5732), so nothing here turns on it.

---

## 4. What the paper may and may not take from this

**May.** (i) `f` — one line, `n‖z̄‖²/‖Z‖_F²`, identical to the mean patient cosine for normalised rows
— is **0.992 on the block §4.1a's `wsi_biology` numbers are read from**, which is a sharper statement of
§3.1's "uncentred rank charges the mean vector as a dimension" than the draft currently makes: RankMe's
own 1.99–3.60 readings on that block are almost entirely the mean. (ii) The reason RankMe's floor is
narrower than ours there is now **derived**: `ln D₁(uncentred) = h(t) + (1−t) ln D₁(centred)`, so the
uncentred statistic is a *shrunk* image of the centred one, and at `t = 0.79` the shrinkage factor is
0.21. (iii) The residualised block has `f = 0`, so **§4.1a's remark that R1 and RankMe "coincide" on
the residualised block is not a coincidence but an identity**, and can be stated as one.

**May not.** Nothing here licenses predicting a floor from `f`. The derived amplification is wrong on
our data by −61%, −33% and +173%, and its ordering across views is wrong. **The `wsi`/`rna` floor gap
is not a centring effect**, and §4.1a's reading — "the catastrophic one-in-five is a property of that
run's WSI encoder" — is supported here and not displaced: rep2's WSI state is measurably the most
collapsed of the five *before any preprocessing* (`f = 0.99498`, uncentred R1 = 1.988 against
3.53–3.60).

---

## 5. The one thing worth carrying beyond this project

Not a law, a **precondition**, and it is checkable in one line from quantities anyone computing a
spectral summary already has:

> Before assuming that centring / residualising / batch-correcting will make a spectral quality metric
> *less* reproducible, check
>
> ```
> |ln((1−t)/t) − ln D_a(centred)| · sd(t)     ≪     (1 − t) · sd(ln D_a(centred))
> ```
>
> across your re-runs, with `t = s_• / (s_• + ‖Z_c‖_*)`. If it fails, the amplification factor
> `1/(1−t)` does not apply and the effect can have **either sign** — on our data it has the opposite
> sign on two views of three. On five same-seed retrains the left side is **1.66 / 0.42 / 0.78** times
> the right, i.e. 4×–17× above the ≤ 0.1 the approximation needs.

And, negatively but usefully: **`1/(1−f)` is wrong by 48× at `f = 0.99` even in the ideal case.** If
anyone is using an energy-share-based rule of thumb for this, the correct ideal-case form is
`1 + √(f/(1−f))/√R2_centred`, which is a square-root divergence damped by the residual spectrum's own
order-2 statistic.

---

## 6. Provenance, and what this rests on

| item | where |
|---|---|
| predeclaration (derivation, assumptions, all nine falsifiers) | `NOTEBOOK_ENTRIES/PREDECLARED_centring_amplification_law_20260804T1750Z.md` @ `d4e344c` |
| harness | `v2/research/rebase/p2/p2_centring_amplification.py` |
| grader (every threshold transcribed; a test asserts each appears in the predeclaration) | `v2/research/rebase/p2/p2_centring_verdict.py` |
| tests (27, incl. the identity asserted exact at rel 1e-10 under its own hypothesis) | `v2/tests/test_p2_centring_amplification.py` |
| real-data readout, 5 reps × 3 views × 2 blocks | `…/figures/data/ws_amp/out/P2_CENTRING_AMPLIFICATION_REAL.json` |
| synthetic sweeps, three conditions | `…/figures/data/ws_amp/out/SYNTH_{t_pinned,stable_mean,unstable_mean}.json` |
| graded verdicts | `…/figures/data/ws_amp/out/P2_CENTRING_VERDICT.json` |
| run logs | `…/figures/data/ws_amp/out/{real,synth_*}.log` |

**Nothing is computed inline.** R1/R2/R3 come from `v2/calibra/spectral.py`; PR, stable rank, α-ReQ and
RankMe from `p2_competing_metrics.py`; the hard rank is `numpy.linalg.matrix_rank`. Neither new module
contains `svdvals`, `linalg.svd`, `eigvalsh` or `linalg.eigh` **at all**, and a test asserts it. The
uncentred form of the three centre-only statistics is obtained by the exact `[X; −X]` reduction, not by
reimplementation.

**What this rests on, plainly.** Five repeats, one arm, one configuration, one seed, one stack; three
views; six statistic families; no interval. The synthetic sweep is synthetic and its perturbation model
is a choice. **Three views is not enough to establish a functional form and it is enough to refute
one**, which is the only thing done here.

### References retrieved for this work (URL for every one)

| reference | retrieved from | on |
|---|---|---|
| Benaych-Georges & Nadakuditi, *The singular values and vectors of low rank perturbations of large rectangular random matrices*, arXiv:1103.2221, doi:10.48550/arXiv.1103.2221 | <https://arxiv.org/abs/1103.2221> | 2026-08-04 |
| Baik, Ben Arous & Péché, *Phase transition of the largest eigenvalue for non-null complex sample covariance matrices*, arXiv:math/0403022, doi:10.48550/arXiv.math/0403022 | <https://arxiv.org/abs/math/0403022> | 2026-08-04 |
| Hill, *Diversity and Evenness: A Unifying Notation and Its Consequences*, **Ecology 54**:427–432 (1973), doi:10.2307/1934352 | <https://api.crossref.org/works/10.2307/1934352> | 2026-08-04 |
| Roy & Vetterli, EUSIPCO 2007, pp. 606–610, doi:10.5281/zenodo.40328 | already VERIFIED at full text, draft §2.6 | prior |
| RankMe — Garrido, Balestriero, Najman & LeCun, ICML 2023, arXiv:2210.02885v3 | already VERIFIED at full text, draft §2.6 | prior |

**The first two are ABSTRACT-LEVEL retrievals.** They are cited for the *existence and shape* of the
finite-rank separation phenomenon, which is stated verbatim in both abstracts. **No constant,
threshold value or theorem number is taken from either**, and none is needed: the separation assumption
(A1) is measured directly by P1 and holds to 0.77%.
