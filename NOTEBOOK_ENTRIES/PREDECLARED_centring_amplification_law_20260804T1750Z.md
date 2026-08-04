# 2026-08-04 17:50 UTC — PREDECLARED: a candidate law for how centring transfers a spectral statistic's reproducibility, derived before any of it is tested

**Logged:** 2026-08-04 17:50 UTC. **Status at logging:** the derivation below is written; **no synthetic
sweep has been run and no new measurement has been taken off the five repeats.** The falsifiers in §5
and §6 are fixed here, before either.

**What I had already seen, disclosed because it makes one of the falsifiers non-blind.** Before writing
this I read `paper/P2_RANK_DRAFT.md` §3.1/§4.1/§4.1a and the whole of
`v2/research/rebase/p2/figures/data/ws_floor/out/P2_ENVELOPE_FLOORS.json` — i.e. **every published
fold in the floors table, for all ten statistics, three views and two blocks.** So:

* **P4 below (the sign of the effect) is NOT blind.** The table already contains a centred/uncentred
  pair per view on the raw block (R1 against RankMe) and I have read it. P4 is therefore recorded as a
  *prediction the derivation makes*, tested on numbers I have seen, and it is labelled that way in the
  result rather than presented as a blind test.
* **P1, P2, P3, P5 and the whole synthetic programme S1–S3 are blind**: they need quantities
  (`f`, `t`, uncentred order-2/order-∞/α-ReQ variants, per-rep dispersions) that exist nowhere on disk
  and that I have not computed.

---

## 0. The claim being tested, stated as someone else's

> Removing a dominant stable component from a representation transfers the reproducibility of a
> spectral summary from the original spectrum to the residual one. If a fraction *f* of the variance
> sits in a common direction that is stable across re-runs, then centring (or residualising, or
> batch-correcting) amplifies the estimation variance of any subsequent spectral statistic by a factor
> increasing in *f* — plausibly ~1/(1−*f*).

**I am not testing `1/(1−f)`.** §1 derives the functional form from the statistic's definition; the
derived form is **not** `1/(1−f)` and disagrees with it structurally, not just in a constant. The
falsifiers test the *derived* form, and S2 additionally tests whether the derived form beats
`1/(1−f)` — because a derivation that fits no better than the naive guess it replaces has earned
nothing.

---

## 1. The derivation

### 1.1 Setup and notation

`Z ∈ R^{n×d}` a representation, rows `z_i`, column mean `z̄ = (1/n) Σ_i z_i`. Write

```
Z = 1 z̄ᵀ  +  Z_c ,        Z_c = Z − 1 z̄ᵀ   (the column-centred matrix)
```

* `σ = (σ_1 ≥ σ_2 ≥ …)` — singular values of `Z_c` (the **residual** spectrum),
* `s = (s_1 ≥ s_2 ≥ …)` — singular values of `Z` (the **original** spectrum),
* `s_• = √n ‖z̄‖₂` — the singular value of the rank-one mean term `1 z̄ᵀ`.

For a non-negative spectrum `v`, L1-normalise `p_k = v_k / ‖v‖₁` and take the **order-a Hill number**
(Hill, *Ecology* 54:427–432, 1973, doi:10.2307/1934352):

```
D_a(v) = ( Σ_k p_k^a )^{1/(1−a)} ,      D_1(v) = exp( − Σ_k p_k log p_k ).
```

This is the family every statistic in this paper's floors table belongs to:
`R1 = D_1(σ)` (Roy & Vetterli, order 1, centred — `v2/calibra/spectral.py CANONICAL`);
`RankMe ≈ D_1(s)` (order 1, **uncentred**);
`R2 = D_2(σ)`; the participation ratio `PR = D_2(σ²)` (order 2 of the *squared* spectrum);
stable rank `= ‖σ‖₂²/σ_1²`, the `a → ∞` member;
hard numerical rank `= D_0(σ)`, the count.

### 1.2 Lemma 1 — the mean term contributes one extra singular value and leaves the rest alone

`ZᵀZ = n z̄ z̄ᵀ + Z_cᵀ Z_c` (the cross terms vanish because `Z_cᵀ1 = 0`). **If `z̄` is orthogonal to the
row space of `Z_c`**, the two summands have orthogonal ranges, so

```
spec(Z)  =  { s_• }  ∪  spec(Z_c)      exactly.
```

### 1.3 Assumption A1 — spike separation (the regime, stated as a boundary not a hope)

Lemma 1's hypothesis never holds exactly. It holds to `O(‖P_{row(Z_c)} z̄‖ / ‖z̄‖)` when `s_•` is large
relative to the bulk edge of `spec(Z_c)`. This is precisely the regime described by

* **Benaych-Georges & Nadakuditi, "The singular values and vectors of low rank perturbations of large
  rectangular random matrices", arXiv:1103.2221** — retrieved 2026-08-04 from
  <https://arxiv.org/abs/1103.2221>. Finite-rank additive perturbation of a large rectangular random
  matrix: the bulk singular value distribution is unchanged, and the perturbation's own singular value
  detaches from the bulk edge **iff** it exceeds a critical threshold, with a corresponding transition
  in the singular vectors.
* the rectangular analogue of the transition first identified for spiked covariance by
  **Baik, Ben Arous & Péché, "Phase transition of the largest eigenvalue for non-null complex sample
  covariance matrices", arXiv:math/0403022** — retrieved 2026-08-04 from
  <https://arxiv.org/abs/math/0403022>.

**A1 therefore has a measurable failure mode with a name**: below the BBP threshold the mean direction
does *not* separate, Lemma 1 breaks, and everything below is void. P1 measures whether A1 holds on our
data rather than assuming it.

### 1.4 Theorem 1 — the exact transfer identity at order 1

Define the **L1 (nuclear-norm) share** of the dominant component

```
t  =  s_• / ( s_• + ‖σ‖₁ ) .
```

Under A1 the L1-normalised spectrum of `Z` is the mixture `(t, (1−t)·q)` where `q` is the
L1-normalised residual spectrum. The Shannon entropy chain rule gives, **exactly**,

```
        ln D_1(s)  =  h(t)  +  (1 − t) · ln D_1(σ) ,        h(t) = −t ln t − (1−t) ln(1−t).
```

i.e. `D_1(uncentred) = e^{h(t)} · D_1(centred)^{1−t}`.

### 1.5 Theorem 1a — general order a ≠ 1

```
        D_a(s)^{1−a}  =  t^a  +  (1 − t)^a · D_a(σ)^{1−a} .
```

### 1.6 Theorem 2 — variance transfer, and the law

Perturb: a re-run moves both `t` and the residual summary. Differentiating §1.4,

```
δ ln D_1(s)  =  (1 − t) · δ ln D_1(σ)  +  [ ln((1−t)/t) − ln D_1(σ) ] · δt
```

and differentiating §1.5, with

```
w_a  =  t^a / ( t^a + (1−t)^a · D_a(σ)^{1−a} )        (so that w_1 = t),
```

```
δ ln D_a(s)  =  (1 − w_a) · δ ln D_a(σ)  +  c_a · δt .
```

**Assumption A2 — the dominant component is stable.** `Var(δt)` is negligible against
`[(1−w_a)/c_a]² Var(δ ln D_a(σ))`. This is what the word *stable* in the candidate law has to mean, it
is the load-bearing assumption, and it is directly measurable (P3b).

**Under A1 and A2 the amplification of the log-scale estimation variance caused by removing the
dominant component is**

```
        A_a  =  sd( ln D_a(centred) ) / sd( ln D_a(uncentred) )  =  1 / (1 − w_a),
        A_1  =  1 / (1 − t).
```

### 1.7 Corollary — the form in terms of *f*, and why 1/(1−f) is wrong twice

Let `f = s_•² / ‖Z‖_F²` — the **squared-norm / energy** fraction in the shared direction, the quantity
the observation quotes as 0.8133. Let `κ = ‖σ‖₁/‖σ‖₂`, and note `κ² = (Σσ)²/Σσ² = D_2(σ) = R2` — the
*order-2 statistic of the residual spectrum itself*. Then `s_• = √f ‖Z‖_F` and `‖σ‖₁ = κ√(1−f) ‖Z‖_F`,
so

```
        t  =  √f / ( √f + √(R2_centred) · √(1−f) )
```

and therefore

```
   ┌──────────────────────────────────────────────────────────────┐
   │   A₁  =  1 / (1 − t)  =  1  +  √( f / (1−f) ) / √(R2_centred) │
   └──────────────────────────────────────────────────────────────┘
```

**Two structural disagreements with the naive `1/(1−f)`, and both are testable:**

1. **The exponent is halved.** As `f → 1`, `A₁ ~ (1−f)^{−1/2}`, not `(1−f)^{−1}`. The divergence is a
   square root.
2. **The constant is not free — it is `√R2` of the residual spectrum.** A dominant component sitting on
   a *rich* residual is far less damaging than the same component on a *thin* one, at identical `f`.
   `1/(1−f)` has no room for this term, and it is the term that will decide whether the form is a
   derivation or a fit.

For orientation only (this is arithmetic on the derivation, not a result): at `f = 0.81` and
`R2_centred ≈ 11`, `A₁ ≈ 1.63`; `1/(1−f) ≈ 5.3`.

### 1.8 What the derivation predicts at the limits, stated before they are checked

* **Order 0 — hard numerical rank.** `D_0` counts. `w_0 = 1/(1+r)`, so `A_0 = (1+r)/r ≈ 1.004` at
  `r = 255`: **no amplification**. Further, if the residual is full-rank in every run then
  `Var(ln D_0(σ)) = 0` and there is no variance to amplify. **The law predicts hard rank has a floor of
  exactly 1.000× and is unamplified by centring. It is a LIMIT CASE the law gets right, not a
  counterexample** — and this prediction is checkable against a number already published (1.000×).
* **`a → ∞` — stable rank.** `w_∞ → 0` once the spike is removed from the numerator's competition;
  amplification `→ 1 + O(σ_1²/‖σ‖₂²)`, i.e. small. Stable rank should be nearly unamplified.
* **Monotonicity in `a`.** `w_a` is increasing in `a` for `t > 1/(1+D)`, so **A_a should be ordered:
  order 2 amplified at least as much as order 1** whenever the spike dominates. This is a sign
  prediction with no free parameters.

### 1.9 When the law predicts centring *improves* reproducibility

If A2 fails — `δt` is large and correlated with `δ ln D_a(σ)` — the uncentred statistic absorbs the
spike's own fluctuation and `A_a` can fall **below 1**. The derivation therefore does not merely permit
the "centring helps" regime, it says exactly where to find it: **wherever the shared direction's mass
share is itself unstable across re-runs.** S3 constructs that regime deliberately.

---

## 2. How the uncentred variant of every statistic is obtained without writing a formula

Only `rankme` in `p2_competing_metrics.py` is uncentred; `participation_ratio`, `stable_rank` and
`alpha_req` all centre internally and expose no flag, and that file is vendored byte-for-byte and must
not be edited. **Four inline-formula substitutions have already been caught in this paper.** So no
statistic is reimplemented. Instead, for any `X`, form

```
        X′  =  [ X ; −X ]        (2n × d)
```

whose column mean is exactly zero, so centring is a no-op on it, and whose Gram is `X′ᵀX′ = 2 XᵀX`, so
`spec(X′) = √2 · spec(X)`. Every statistic here is scale-invariant, therefore

```
        stat_centred( [X;−X] )  ==  stat_uncentred( X )        exactly.
```

**Guard, declared as a stopping condition:** `effective_rank([X;−X], variant=R1)` must equal
`effective_rank(X, variant=R1_uncentred)` to 1e-9 relative, and `rankme([X;−X])` must equal
`rankme(X)` to 1e-9 relative, on every rep × view × block. If either fails, the trick is not exact on
this data and **nothing in the order-2/order-∞/α-ReQ columns may be used.**

---

## 3. What will be measured on the five repeats

`~/e0_run/d1_envelope/rep{1..5}.npz`, the same five exports `p2_envelope_floors.py` reads, same test
partition, same cancer + pooled-TSS cross-fitted residualisation at seed 42. CPU only, thread-capped.
Per rep × view (`wsi_biology`, `rna_biology`, `full_biology`) × block (`raw`, `residualised`):

| symbol | definition | how |
|---|---|---|
| `f` | `n‖z̄‖² / ‖Z‖_F²` | the energy share of the shared direction — **the one-line quantity the law is supposed to be predicted from** |
| `cos̄` | mean off-diagonal cosine of row-normalised `Z` | the patient-to-patient collinearity the observation quotes (0.8008 / 0.274) |
| `t` | `s_• / (s_• + ‖σ‖₁)` | measured from the actual centred spectrum |
| `t_model` | `√f/(√f + √R2·√(1−f))` | §1.7 — must agree with `t`; disagreement means A1 is failing |
| `D_a` centred / uncentred | R1, R2, PR, stable rank, α-ReQ, hard rank, RankMe | imported from `v2/calibra/spectral.py` and `p2_competing_metrics.py`; uncentred via §2 |

Dispersion across the five reps is reported **both** as `sd(ln ·)` (what Theorem 2 is about) and as
`ln(max/min)` (what the paper calls a floor), and **both on all five reps and on the four concordant
reps with rep2 removed**, because a ratio of dispersions on n = 5 with one outlier is otherwise a
statement about one run.

---

## 4. The synthetic sweep, and its falsifiers — fixed now

`Z_r = a·1μᵀ + W_r` for runs `r = 1…50`, `n = 2766`, `d = 256` (our shapes), `W_r` a fixed bulk
`W_0` plus an independent perturbation `E_r` of controlled relative size, `μ` a fixed unit vector
orthogonal to nothing in particular (drawn once). `a` swept so that `f ∈ {0, 0.05, …, 0.95}`.
Two conditions:

* **stable-mean** (A2 holds by construction): `a` identical in every run.
* **unstable-mean** (A2 deliberately broken): `a_r = a·(1+η_r)` with `η_r` at the same relative
  dispersion as the bulk perturbation.

| # | prediction | **falsified if** |
|---|---|---|
| **S1** | stable-mean: `A₁^obs` matches `1/(1−t̂)` | median relative error > 10% over `f ≤ 0.90`, or any single `f ≤ 0.90` off by > 20% |
| **S2** | the derived form beats the naive one | median \|rel. error\| of `1/(1−f)` **≤** that of `1/(1−t̂)`. *If S2 fails the derivation has earned nothing and I say so.* |
| **S3** | unstable-mean: `A₁^obs < 1/(1−t̂)`, and dips below 1 at small `f` | `A₁^obs ≥ 1/(1−t̂)` in more than 2 of the 20 `f` values |
| **S4** | order ordering (§1.8): `A_2 ≥ A_1 ≥ A_∞ ≥ A_0` under stable-mean | violated at more than 2 of the 20 `f` values |

---

## 5. The real-data falsifiers — fixed now

| # | prediction | **falsified if** |
|---|---|---|
| **P1** *(blind)* | **A1 holds**: `\|ln D_1(s) − h(t) − (1−t)ln D_1(σ)\| / ln D_1(s) ≤ 0.02` per rep × view × block | > 5% in more than 2 of the 30 cells. *A1 failing is a legitimate outcome and kills the law on this data rather than the law in general.* |
| **P2** *(blind)* | **the law itself**: on the raw block, `A₁^obs` within ±25% of `1/(1−t̄)` in all three views | any view outside ±50%, **or** the three views' rank order by `A^obs` disagrees with their rank order by `A^pred` |
| **P3** *(blind)* | **order dependence**: across the statistic family, `A_a^obs` tracks `A_a^pred` | Spearman ρ ≤ 0 over the raw-block statistic × view cells |
| **P3b** *(blind)* | **A2 holds on this data**: `sd(δt)` contributes < 25% of `Var(δ ln D_1(s))` | it contributes > 50% — in which case A2 is false here and P2 is not a test of the law |
| **P4** *(NOT blind — see header)* | **sign**: `A_a ≥ 1` wherever A2 holds | any raw-block cell with `A^obs < 1` by more than its own bootstrap error. *I have already read R1 = 1.0228× against RankMe = 1.0299× on `rna_biology` raw, which is `A^obs ≈ 0.76`. If that survives being recomputed against a true `R1_uncentred`, **P4 is already falsified** and the law does not hold unconditionally.* |
| **P5** *(blind)* | **the motivating observation**: the *view* difference in centred floors is what centring did — `ln(floor_centred) ≈ A₁(view) · ln(floor_uncentred)` across the three views | the uncentred floors already differ across views by more than 70% of the log-scale spread of the centred floors. *In that case the `wsi` vs `rna` contrast is not a centring effect at all and observation 2 has nothing to do with observation 1.* |

**P5 is the decisive one.** Observations 1 and 2 in the framing are two different comparisons —
1 is centred-vs-uncentred within one view, 2 is view-vs-view with centring held fixed — and the law is
only a single explanation of both if P5 passes.

---

## 6. What each outcome means, written before the outcome

* **P1 and P2 pass, P5 fails** → the transfer identity is real and useful, but it explains the
  RankMe/R1 gap **only**, and the `wsi`/`rna` floor difference is a separate phenomenon (most likely:
  one run's WSI encoder actually diverged). Headline: *the law is real and the observation that
  motivated it was two things, not one.*
* **P1 fails** → A1 does not hold on our states; the derivation may be correct and untestable here.
  Headline: *not tested, and here is why.*
* **P2 fails with S1 passing** → the form is right in the controlled setting and wrong on real
  re-runs, which localises the failure to A2. Headline: *the "stable" in "dominant stable component"
  is doing all the work and is not satisfied by retraining.*
* **S2 fails** → the derivation fits no better than the naive guess. **Report as a coincidence.**
* **Everything passes** → report it, and report that three views and ~6 statistics is still a small
  number of points, and that the sweep is synthetic.

**Nothing here may be re-read after the numbers arrive.** No third reading is invented after the fact;
if an outcome falls between the rows, the magnitudes are reported and no verdict is given.

---

## 7. References retrieved for this entry (URL given for every one)

| reference | retrieved from | on |
|---|---|---|
| Benaych-Georges & Nadakuditi, *The singular values and vectors of low rank perturbations of large rectangular random matrices*, arXiv:1103.2221, doi:10.48550/arXiv.1103.2221 | <https://arxiv.org/abs/1103.2221> | 2026-08-04 |
| Baik, Ben Arous & Péché, *Phase transition of the largest eigenvalue for non-null complex sample covariance matrices*, arXiv:math/0403022, doi:10.48550/arXiv.math/0403022 | <https://arxiv.org/abs/math/0403022> | 2026-08-04 |
| Hill, *Diversity and Evenness: A Unifying Notation and Its Consequences*, **Ecology 54**:427–432 (1973), doi:10.2307/1934352 | <https://api.crossref.org/works/10.2307/1934352> | 2026-08-04 |
| Roy & Vetterli, EUSIPCO 2007, pp. 606–610, doi:10.5281/zenodo.40328 | already **VERIFIED** at full text for this paper, §2.6 | 2026-08-05 (prior) |
| RankMe — Garrido, Balestriero, Najman & LeCun, ICML 2023, arXiv:2210.02885v3 | already **VERIFIED** at full text for this paper, §2.6 | prior |

**Abstract-level only** for the first two: I retrieved the arXiv abstract pages, not the full texts.
They are cited for the *existence and shape* of the finite-rank separation phenomenon, which is in both
abstracts verbatim; **no numerical constant, threshold value or theorem number is attributed to
either**, and none is used in the derivation. Nothing here needs a constant from them: A1 is checked
empirically by P1.
