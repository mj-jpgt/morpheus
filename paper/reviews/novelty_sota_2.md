# Adversarial Review — SOTA Claim / Diagnostic Framing (novelty_sota_2)

**Scope.** Does the paper claim or imply SOTA where it should not? Is the diagnostic/negative-result
framing defensible without a SOTA number? What is the minimum real-data T4 result needed for the paper
to stand? Files read: `paper/main.md`, `paper/05_experiments_results.md`,
`paper/07_implementation_status.md`.

---

## 1. Does the paper claim or imply SOTA anywhere it should not?

**Verdict on this axis: No SOTA is claimed — the paper is, if anything, exemplary in its restraint. The
real exposure is the opposite: one efficacy overclaim ("provable" fix) and a self-defeating benchmark posture.**

What I checked and where an implied SOTA could hide:

- **The headline model loses on its own benchmark, and the paper says so.** `morpheus_v2_anchored` is the
  worst method on global Pearson (0.327, the lowest row of T2), loses head-to-head to the MLP-CLIP baseline
  (−0.021 global / −0.022 within-cancer, winning only 32%/19% of targets in T3), and trails on retrieval
  (0.060 vs 0.066). A paper that wanted to sneak in a SOTA claim would not publish these tables. This is
  disarming honesty and directly neutralizes the "where's your win?" reflex. Good.
- **R6 explicitly disavows SOTA**: "which is exactly why we do **not** claim a state-of-the-art method."
  The entire C2 thesis (method-invariant +0.07 ceiling) is logically incompatible with a SOTA claim, so the
  paper cannot accidentally imply one without contradicting itself. This internal consistency is a strength.
- **No leaderboard language, no "outperforms," no "best" anywhere in the abstract/intro/discussion.** Correct
  for the genre.

**The one place the paper overclaims is not SOTA — it is the word "provable."**
Abstract C3 and §3.3 call the fix "prescriptive and **provable**." The Proposition (§3.2) only proves the
*negative* control: a per-dimension variance floor **cannot** restore rank. It does **not** prove the
*positive*: that the off-diagonal decorrelation term **will** restore rank on the real biology head. The
positive direction is asserted heuristically ("driving `L_cov` toward zero pushes `C` toward diagonal, which
raises erank") and supported only by a **synthetic** PoC (rank 10.3→21.2, §07). Selling a fix as "provable"
when only its negative control is proven is an efficacy overclaim that a hostile reviewer will catch
immediately, and it is more dangerous than a SOTA claim because it invites the referee to demand the missing
proof/experiment. **Fix: restrict "provable" to the variance-floor impossibility result; describe the
decorrelation fix as "prescribed and synthetically validated, real-data validation queued (T4)."**

Secondary, smaller: the synthetic PoC number (21.2 on a `hidden=64`, `B=24` batch) sits at roughly the
batch-size erank ceiling (erank ≤ min(n−1, d) ≈ 23). It is a legitimate "term is not a no-op" demonstration,
but it **cannot** distinguish "restored toward healthy" from "saturated at the batch cap," and should not be
read as evidence of the *magnitude* of real-data recovery. State this explicitly so it is not later attacked
as an inflated proxy.

---

## 2. Is the diagnostic / negative-result framing defensible without a SOTA number?

**Yes — the genre is well established and the paper fits it — but defensibility rests on C1+C2 being
*empirically solid*, and C1's geometry is currently single-seed, which is the real threat, not the absence
of SOTA.**

Precedent (negative-result / diagnostic / audit papers that landed at top venues with **zero** SOTA):

- **Venet et al. 2011, PLOS Comp Biol** — "Most random gene signatures are significantly associated with
  breast cancer outcome." Pure negative result, no new method; high-impact and field-defining. This is
  *directly* the lineage of the paper's random-gene control and its "the floor is high" argument.
- **Geirhos et al. 2020, Nature Machine Intelligence** — "Shortcut learning in deep neural networks." A
  Perspective/diagnostic with recommendations, no new model. Establishes that naming+characterizing a failure
  mode is a publishable contribution on its own.
- **"Are Neurons Actually Collapsed?" ICML 2023** — diagnostic analysis of (neural) collapse geometry, no
  SOTA. Same genre as C1.
- **Howard et al. 2021, Nature Communications** — site-signature audit; conceded prior art here, and proof
  that the confound-audit genre publishes without a competitive method.
- **"Auditing Data Leakage in WSI Multimodal Benchmarks" (2607.12278)** and **NRC1 (Andriopoulos, NeurIPS
  2024)** — recent, same-domain audit/prevalence papers, no new SOTA method.

So the framing is defensible **in principle**. The exposure is in the *strength of the positive diagnostic
evidence*, and there are two gaps a referee will press:

1. **C1 geometry is single-seed (conceded, §6.4).** A "fingerprint" is a *measurement claim*; a single-seed
   measurement is thin for a headline. Every precedent above earns its negative result through repetition
   (Venet: thousands of random draws; Howard: many sites/folds). A 14× gap (84.3 vs 6.0) is unlikely to be
   noise, but "unlikely" is not "shown." **This — not the missing SOTA — is the single biggest risk to the
   diagnostic framing.** The multi-seed geometry run is cheaper than T4 and shores up the actual headline.

2. **The "confounder-agnostic detector" claim is asserted, not validated as a detector.** The paper claims
   biology-head erank "flags cohort-confounded alignment without knowing the confounder a priori" (C1, §3.3,
   R1). But there is no demonstration that erank *discriminates* confounded from non-confounded setups: no
   case where a head trained on a genuinely non-confounded / higher-rank target *stays* high-rank, no
   detector ROC/threshold. The only "negative" is the identity head, which differs by *supervision type*
   (contrastive vs regression), not by *confound status*. As stated, the mechanism reduces to "regression
   onto a low-rank target collapses to that rank; contrastive supervision does not" — which is NRC1 (cited)
   plus Wang/Isola (cited). The *diagnostic-for-confounding* leap needs at least one positive control (a
   low-rank target that is **not** cohort-confounded but still collapses, or a confounded setup where erank
   stays high) to separate "low-rank-target collapse" from "confound fingerprint." Without it, R1/R4 ("this
   is DECAT / NRC1 restated") are not fully closed. This is the novelty-of-the-diagnostic question, and it is
   answerable with analysis the team likely already has, but it is not in the draft.

Net: the framing survives the *SOTA* objection cleanly. It is exposed on *reproducibility of C1* and
*validation that erank is a confound detector rather than a low-rank-target detector*.

---

## 3. Minimum quantitative real-data T4 needed for the paper to stand

Distinguish two thresholds, because they gate different claims.

**Threshold A — for the paper to publish as a *diagnostic* (C1+C2), T4 is NOT required.**
The load-bearing missing result at this threshold is not T4 at all; it is **multi-seed C1 geometry**:
biology-head vs identity-head erank across seeds {42,43,44} (± the other configs already trained: no-anchor,
v1), showing the collapse gap is seed-stable. That is the cheaper run and it directly removes the §6.4
single-seed caveat. With it, C1+C2 stand on their own (all C2 numbers are already multi-seed), C3 is demoted
to "prescription + synthetic PoC + queued validation," and "provable" is softened per §1. This is a
publishable diagnostic paper without any A100 T4 result.

**Threshold B — for the *fix* (C3) to stand as a demonstrated result rather than a promise, T4 must show a
joint, not a single number.** The minimum real-data T4 is the pre-registered triplet, and both legs are
mandatory:

1. **Rank recovery, separated from the negative control.** Biology-head erank (`full_biology` /
   `wsi_biology` / `rna_biology`) with **F-R2 on vs off vs variance-floor arm**, on the real held-out cohort.
   Minimum bar: (i) F-R2-on erank rises *clearly and seed-stably* well above the ~5–6 target-intrinsic-rank
   floor — directional-with-separation is the floor; a rise toward a substantial fraction of the identity
   regime is what makes it convincing; and (ii) the **variance-floor arm stays collapsed (~5–6)**, which is
   the empirical confirmation of the Proposition and is *non-optional* — it is what converts "provable
   negative control" from theory to demonstrated fact.
2. **Specificity decoupling (the falsifiable core).** Control-adjusted within-cancer specificity
   (`macro_cancer_pearson`, real − random-gene null) must stay statistically **unchanged at ~+0.07** when
   rank is restored. This is the more dangerous leg: if fixing rank *moves* Δ_spec, the entire "rank and task
   score are decoupled" thesis (C2 + the fingerprint's independence-from-score argument) breaks — a rise
   would falsify the "method-invariant ceiling," a fall would make the fix harmful. **The paper's headline
   survives only if rank↑ AND Δ_spec flat both hold.** This must be reported with a seed-level dispersion so
   "unchanged" is a tested null, not an eyeball.

Minimum viable if the A100 run is partial at submission: even a **single-seed** real-data
on/off/variance-floor triplet exhibiting rank↑(bio) with Δ_spec flat converts C3 from "synthetic-only +
argued" to "demonstrated on real data," which is the threshold for the fix claim to stand. Multi-seed
strengthens it and simultaneously discharges the C1 single-seed caveat — so the highest-value single run is
the multi-seed T4, which pays down Threshold A and B at once.

**What would sink the paper:** submitting with C3 as-is ("provable fix"), no T4, *and* single-seed C1 — that
is a diagnostic claim resting on one seed plus a fix claim resting on a synthetic toy, with an efficacy word
("provable") the evidence does not support. Any one of {soften "provable" + multi-seed C1} rescues the
diagnostic; the fix claim needs the real T4 triplet.

---

## Minor / adjacent

- **Panel-size framing wobble.** §5 evaluates over "180 MSigDB Hallmark targets"; abstract/method describe a
  "50-D Hallmark" target and "intrinsic rank ~5–6." The 180 is the full target panel (real + heldout +
  random-control groups); the 50-D/rank-5–6 is the supervision manifold. Both are internally correct but a
  fast reader conflates them — state the relationship once. (Consistency lane, not SOTA, flagged only because
  it touches the "what is the benchmark" question a SOTA-hunting reviewer asks.)
- **Elevate the no-SOTA argument out of the risks table.** The "we deliberately have no competitive method,
  and that is the finding" argument currently lives only in R6. For an ML audience with a SOTA reflex, put a
  one-sentence version in the intro contribution list so it frames the reader before the tables do.

---

## SUMMARY

- **No SOTA is claimed or implied** — the paper's own model deliberately underperforms every baseline in its
  own tables (T2/T3/retrieval) and R6 explicitly disavows SOTA; the framing is a well-precedented
  negative-result/diagnostic genre (Venet 2011, Geirhos 2020, Howard 2021, NRC1, WSI-leakage audits) that
  routinely publishes without a competitive number. The SOTA axis is a strength, not a liability.
- **The genuine overclaim is "provable fix," not SOTA**: the Proposition proves only the negative control
  (variance floor cannot fix rank); the positive efficacy of the decorrelation term rests on a synthetic PoC
  (rank 10.3→21.2, near the batch cap) — soften "provable" to cover only the impossibility result.
- **Minimum to stand**: the diagnostic (C1+C2) can publish *without* T4 if C1 geometry is made multi-seed
  (cheap, removes the sole single-seed caveat) and the erank-as-confound-detector claim gets one positive
  control; the *fix* (C3) requires the real-data T4 triplet showing biology erank↑ (F-R2 on), erank flat
  (variance-floor arm), AND Δ_spec unchanged at ~+0.07 — the rank↑/score-flat joint is the falsifiable core.

**VERDICT: paper_ready=conditional** — SOTA posture is sound (no unearned claim); acceptance is gated on
(1) softening "provable" to the negative-control result, (2) multi-seed C1 geometry to lift the diagnostic
off a single seed, and (3) the real-data T4 rank↑/Δ_spec-flat triplet before C3 can be stated as a
demonstrated fix rather than a queued promise. Items (1)+(2) rescue the diagnostic on their own; (3) is
required only for the fix claim.
