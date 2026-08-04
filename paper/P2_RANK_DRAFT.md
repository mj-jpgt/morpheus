# Effective rank cannot resolve its own re-measurement: a reproducibility floor inside the regime RankMe reserves for itself

**A negative result about a representation-geometry proxy, reported with the experiment that went
against it first, and with its own prior art ahead of its own evidence.**

*Working draft, 2026-08-04. Companion to `paper/P1_CALIBRA_DRAFT.md`, which asks what an **analysis**
would have missed; this asks what a **representation-geometry summary** fails to tell you. Every
number in this document traces to a named artifact, notebook entry or evidence file in this
repository; each table carries a `provenance` line. Numbers that do not exist are marked as not
measured rather than estimated. Every citation carries an explicit verification status (§2.6). Three
fabricated citations have previously contaminated this project; §2.6 is not a formality.*

> ### Status — read before anything else
>
> **1. This draft's previous claim has been falsified by our own experiment and has been removed.**
> Earlier versions were organised around *"effective rank does not track information content"* and
> around a preregistered necessity test that was expected to break RankMe's *necessary-but-not-
> sufficient* hedge. The test ran. It went the other way: `programme_only` (high rank) beats
> `programme_free` (low rank) on the measured molecular channel in **all three seeds** — interval-backed
> 3/3 on the patient bootstrap and **2/3 on the conservative cancer-cluster bootstrap** — under the
> canonical statistic — and under every canonical statistic we compute; an earlier draft's "2/3
> under R2" qualification rested on a mislabelled statistic and is withdrawn (§4.5a). Per the
> preregistered
> instruction in
> `NOTEBOOK_ENTRIES/PREDECLARED_D1_necessity_test_20260803T2300Z.md` (outcome **O2**), this is
> reported **at abstract prominence as a negative for the paper's generality** and is §4.7, ahead of
> every result that favours us.
>
> **2. The claim that survived is different, narrower, and better powered.** It is stated in §1.3 and
> carried by §4.1–§4.5. It does not depend on any sign count.
>
> **3. Pending slots — all three closed as of 2026-08-04, and two of the three cost us something.**
> (a) The D1 intervals were never missing: the stratified paired bootstrap existed all along and was
> hidden by a stale absolute path in the audit chain. §4.7.2 reports **both** estimators and the
> conservative one is less favourable — 3/3 on the patient bootstrap, **2/3 on the cancer-cluster
> bootstrap**, seed 43's interval touching zero at +0.0006.
> (b) **§4.1's retraining envelope is no longer n = 1.** Five identical same-seed retrains report a
> **3.295×** rank floor against a **1.055×** channel spread, bimodally distributed. It moves our own
> headline count from six of seven to seven of seven — *in our favour*, reported with the scepticism
> that requires — and it makes **§4.7's necessity result unresolvable rather than refuted**, which is
> the reading that was predeclared before the measurement existed.
> (c) **§5.2's momentum sweep now has its seed replication, and §5 has been rewritten around it**
> (`mseed_*`, three seeds per momentum, 500 steps). Every m = 0.999 seed exceeds every m = 0 seed, so
> §5.3's predeclared disjunction resolves in favour of **separation**, and the single-seed defect — a
> **hardcoded harness parameter**, not a design choice — is closed. **But the worst-case separation is
> 3.29× against §4.1's 3.295× floor**, so by §4.1's own criterion that fix's *rank* difference is **not
> resolvable either**. **That is now §5.4, a subsection of its own rather than a caveat**, and §5.2 no
> longer rests the fix on a rank ratio: it rests on a **binary training outcome** — `programme_free`
> produced no completed, uncollapsed, exportable run before the fix and three of three after it. The
> numbers are in `NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §3 and
> `NOTEBOOK_ENTRIES/p2_section5_rewritten_around_the_momentum_replication_20260804T2000Z.md`.
>
> **5. §4.6's ground truth is a coordinate choice, and §4.6a is new.** The D2 arm contrast the
> selection-rule table scores against exists on gene-set targets and on **none** of the five other
> molecular target blocks on disk. Re-scored against each block in turn, every metric row's count
> moves, the ordering between the two metrics the table is quoted for reverses, and our own statistic
> can be handed a nominally significant 6/6 by the choice of block alone. `D2_RESULT.md` §6 has been
> narrowed accordingly.
>
> **4. Every effective-rank number in this draft was recomputed under one canonical definition on
> 2026-08-04** and every surviving instance reproduces exactly. Two instances are marked
> `[NOT RECOMPUTABLE — artifact never existed]` and one `[NOT RECOMPUTED — needs a GPU forward pass]`;
> they are not carried forward as though they had been measured under the canonical definition. See
> §3.1 and
> `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`.
>
> **6. §4.1a is new: the paper audited exhaustively against its own criterion, and it costs us —
> including a scope error in the first version of the audit itself.** §4.1's floor had been applied
> to this paper's own numbers five times, each instance found separately and late. §4.1a enumerates
> **all 56** rank comparisons the paper makes or relies on — Results, the worked example, the
> appendices, the figure plan, `QUEUE_ANCHORING.md` and `LIVENESS_GATE_DESIGN.md` — each judged
> against the floor measured on **its own statistic and its own block**. The first version judged
> them all against the **one** floor that existed, canonical R1 on the exported `wsi_biology` block,
> and printed "fails" for twenty-four rows the criterion could not reach. **Every floor recoverable
> from the five same-seed repeats has now been measured** — ten statistics, three views, two blocks,
> from exports that already existed — and the counts are restated on that basis: **of the 25
> selections between candidate configurations, 13 fail a floor their own statistic and block license,
> 11 cannot be judged at all, and exactly 1 clears — RankMe as published, whose floor is 1.811×
> against our own statistic's 3.111× on the same five runs.** The floor turns out to be a property of
> the statistic (1.000× to 3.295× on one block) and of the **view** (3.295× on `wsi_biology`, 1.019×
> on `rna_biology`): the one divergent repeat lost its WSI-view rank and kept its RNA-view rank to
> within 2%. Four blocks still have no floor and none can be recovered without a GPU — **every rank
> number in §5 is on one of them**. It found a sixth instance the draft had not stated (§4.7.4's second violation, P44
> against I43, 3.07×), two §5.2 claims that read an ordering off a 1.04× difference and must be
> restated as the equalities they are, and one like-for-like same-seed pair in §5.4's own regime that
> §5.4 says does not exist — n = 2, concordant, and explicitly **not** quotable as a floor. It also
> resolves §5.1's instance 2, which had been flagged and left. The list is machine-checkable at
> `v2/research/rebase/p2/floor_audit.json` and a test fails if a ratio in the table disagrees with
> its source.
>
> **7. §4.9a is new, and so is F9.** *"`feature_decorrelation` is defective"* was **conditional on a
> query-written queue**. With a momentum key encoder the same term raises effective rank
> monotonically across three levels while the RNA-view mutual cosine — a direct measure of the
> collapse rank exists to detect — rises monotonically **with** it, co-measured on the identical
> runs. The rank change (1.85× R3 / 1.94× R1) is inside §4.1's floor, so the monotonicity and the
> cosine carry it and the magnitude does not.

---

## Abstract

Effective rank — the exponential of the Shannon entropy of a matrix's L1-normalised singular values
(Roy & Vetterli, EUSIPCO 2007, pp. 606–610) — is proposed as a label-free criterion for assessing
self-supervised representations and for selecting hyperparameters without labels (RankMe; Garrido,
Balestriero, Najman & LeCun, ICML 2023, arXiv:2210.02885). **We claim no discovery of the negative,
and we report our own failed prediction first.** We ran a preregistered test of RankMe's
*necessary-but-not-sufficient* hedge — two arms of one method differing only in objective, three
seeds, paired bootstraps — and the hedge held: the higher-rank arm carried the larger held-out
molecular channel in all three seeds, in the same direction as its rank advantage, interval-backed
3/3 under a patient bootstrap and **2/3 under the conservative cancer-cluster bootstrap** (both are
quoted; the conservative one is the one to weight). That is a negative for this paper's
generality and we state it here, at the prominence a confirmation would have received — including
after a subsequent measurement made it **unresolvable rather than refuted**: the rank differences that
test reads (1.74×–3.25×) all sit inside the retraining floor of (i) below, so the experiment decides
nothing about rank in **either** direction. That reading, and the refusal to treat it as a rescue, were
fixed in a predeclaration written before the floor was measured. It joins
published prior art we cannot claim novelty over: Aldeneh, Thilak, Higuchi, Theobald & Likhomanenko
(ICASSP 2025, DOI 10.1109/ICASSP49660.2025.10889651) already report that *"rank does not reliably
predict the best-performing layer for specific downstream tasks, as lower-ranked layers can
outperform higher-ranked ones"* — a within-method selection-rule failure with a published
low-rank/high-information instance.

What survives is a different and better-powered claim, about **usefulness rather than truth**:

> **Effective rank is unusable as a selection signal, because its between-arm differences are smaller
> than its own within-arm reproducibility floor — inside the regime its proponents explicitly reserve
> for it** (*"RankMe should however only be used to compare different runs of a given method"*).

Four measurements carry it, on 12 frozen artifacts from 4 arms × 3 seeds of two matched-arm
experiments in cross-modal morphology → bulk-transcriptome learning on 2,766 held-out TCGA patients.
**(i)** **The single cleanest observation is one configuration retrained five times at the same seed,
with GPU non-determinism the only source of variation: effective rank spreads 3.295× while the
molecular channel those same five runs carry spreads 1.055×.** The spread is **bimodal** rather than
smooth — four repeats agree to within 2% and one lands at a third of them, with its rank down 3.3× and
its channel down 5% — so the usable statement is that effective rank is **reproducible about 80% of the
time and catastrophically not about 20% of the time**, on identical inputs. **All seven** between-arm
rank differences this project has ever measured (**1.004× to 3.246×**) fall inside that floor, each
judged against the floor measured on its own preprocessing block. The number is a floor twice over —
it is measured on this project's *stable* arm, and same-seed repeats exclude seed variation entirely —
and it is one measurement, of five repeats, on one arm, at one seed, which is why (ii) carries the
argument. Its effect on our own headline count (six of seven to seven of seven) runs in our favour and
is reported with the scepticism a result running the other way would receive. **(ii)** A variance decomposition over the 12 artifacts puts effective rank at **34.5% arm /
65.5% training-seed nuisance**, with the arm term not significant (F(3,8) = 1.41); the molecular
channel those same artifacts carry is **98.0% arm**, F = 128.2. **(iii)** At a matched training step,
one arm spans **6.05×** across five seeds while its sibling spans **1.18×** — the floor is a
property of the *arm*, and is worst exactly where the arm is interesting. **(iv)** The instability is
in **training, not estimation**: patient-resampling gives a measurement SD of ≈0.1 on a rank of 25,
re-exporting a surviving checkpoint reproduces to five significant figures, and a controlled
short-horizon repeat at fixed seed spans 4.7%.

Two further findings constrain any use of the number at all. The between-arm *verdict* flips with the
choice of rank statistic (2 of 6 pairs, corrected from a published 3 of 6 — §4.5a), with whether the
block is confound-residualised (2 of 3 D2
seeds), and with which co-trained view of the same model is measured (2 of 6 pairs); the information
verdict never flips under any of them. And on a seven-level dose–response with a **zero-parameter**
representation, rank under-reports the information loss by a factor between **1.95× and 21.5×**
depending on those same choices — **21.5×** on the block the channel is actually read from, a
correction that moved the result *in our favour* relative to the 3.7× previously published and is
reported for that reason.

We report what cuts against us at equal prominence. Centred effective rank does fall to ~1 under
total collapse, so the collapse-diagnostic use is supported; §5 is a worked example of exactly that
use — not a separate contribution — in which a training gate and a queue fix are graded in rank. That
worked example returns the paper's standard against the paper: the queue fix's own seed-replicated rank
separation is **3.29×** against the **3.295×** floor of (i), so **it is inside our own floor too**, and
§5.4 rests that fix on a binary training outcome — an objective that never once completed training
uncollapsed, and does so on three of three seeds after the fix — rather than on the rank ratio that
first drew attention to the problem. At effective rank 9.1 of a nominal
256, two representations still read held-out channels of 0.5983 and 0.4757 against a permutation null
of 0.140, so "low rank means little information" is false outside total collapse. We evaluated three
published alternatives on the same artifacts and **none is better**: LiDAR picks the
information-poorer arm 0/3 on the in-scope contrast at every δ over eight orders of magnitude.
Selection-rule counts on six pairs are reported and are explicitly underpowered — a *flawless* 6/6
would give p = 0.031 — and they are also **not stable under the choice of ground truth**: the arm
contrast they are scored against exists on gene-set targets and on none of the five other target blocks
on disk, and re-scoring against each block in turn moves every metric row's count, reverses the ordering
between the two metrics the table is quoted for, and can hand our own statistic a nominally significant
6/6 by nothing but the coordinate system the exam is written in. We report that we can manufacture our
own significance this way because it is the strongest available reason not to read the counts at all;
the argument rests on the variance decomposition instead. And we
apply the paper's own standard to the paper: our smallest arm difference (0.0705) is about **half**
the real-versus-random-control margin of 0.139, so the effect we measure is itself small relative to
its instrument's floor. That is §4.1's argument arriving from the channel side, it points at us, and
we state it before it is put to us.

### Short abstract (~200 words)

Effective rank (Roy & Vetterli 2007) is proposed as a label-free criterion for representation quality
and hyperparameter selection (RankMe, ICML 2023). We first report a failed prediction: a preregistered
within-method test of RankMe's necessity hedge confirmed it in all three seeds — interval-backed 3/3
under a patient bootstrap and 2/3 under the conservative cluster bootstrap — though that test's own rank
differences later proved to sit inside the floor below, leaving it unresolvable rather than refuted.
What survives is a claim about
usefulness, not truth: **effective rank's between-arm differences are smaller than its own within-arm
reproducibility floor, inside the same-method regime RankMe reserves for itself.** Retraining one
configuration five times at a fixed seed spreads effective rank **3.295×** while spreading the molecular
channel those same runs carry **1.055×**; the spread is bimodal — four repeats within 2%, one at a third
— so rank is reproducible ~80% of the time and catastrophically not ~20% of the time. All seven
between-arm rank differences ever measured on this project (1.004×–3.246×) lie inside that floor. Over
12 artifacts (4 arms × 3 seeds), rank is
34.5% arm and 65.5% seed nuisance (arm term not significant, F(3,8) = 1.41) while the molecular
channel is 98.0% arm (F = 128.2). At a matched step one arm spans 6.05× across seeds and its sibling
1.18×, so the floor belongs to the arm. Measurement noise is negligible (SD ≈ 0.1 on rank 25): the
instability is in training. The between-arm verdict also flips with the rank statistic, the
preprocessing block and the measured view; the information verdict never does. Three published
alternatives were computed on the same artifacts; none is better.

---

## 1. Introduction

### 1.1 A scalar that is cheap, label-free, and therefore load-bearing

Evaluating a learned representation properly requires labels, a downstream task and a held-out split.
All three are expensive, and in cross-modal biological settings the downstream label is often the
scarce quantity the representation exists to predict. This creates strong demand for a *label-free*
scalar computed on the representation matrix alone.

The spectrum supplies the obvious candidate, and it is genuinely attractive. Given a patient × feature
matrix, take its singular values, normalise them to a probability vector, and report the exponential
of their Shannon entropy — a smooth interpolation between 1 (all mass on one direction) and the
matrix's rank (a flat spectrum). It is cheap, it has a clean spectral definition, it needs no labels,
and it returns a number in every circumstance. This is Roy & Vetterli's effective rank, and RankMe
adopts it verbatim as a label-free criterion for self-supervised representation quality and
hyperparameter selection (§2.1). A paper that opens by sneering at a metric loses the readers who
use it; the appeal is real and the rest of this paper is about what that appeal does not buy.

The intuition licenses three distinct practices, and they are not equally defensible:

1. **As a collapse diagnostic.** "Rank has fallen to ~1, so the model has collapsed." Supported by the
   evidence here (§4.10), the only use we end up recommending, and the use demonstrated at length in
   §5 — with the caveat that a cheaper statistic does the job better.
2. **As a training objective or a target of regularisation.** Anti-collapse regularisers penalise
   off-diagonal feature covariance in order to raise the occupied dimensionality (VICReg, ICLR 2022;
   Barlow Twins). §4.9 reports what happened when we did this.
3. **As a comparative selection signal** — "configuration A has higher effective rank than B,
   therefore A is the one to keep". This is the use RankMe formalises and restricts, and it is the
   only one this paper argues against.

Use (3) is the one cheap enough to be tempting and the one that, if unreliable, is *silently*
unreliable: a rank number never fails, never returns `NaN`, and never announces that the difference
you are reading is smaller than the difference you would get by training the same configuration
twice.

### 1.2 What is already known, and what we are not claiming

**We claim no discovery of the negative.** The prior art is explicit, some of it is peer-reviewed,
some of it is in-scope, and we state all of it before our own results rather than after (§2.2).

The single most constraining item is **Aldeneh, Thilak, Higuchi, Theobald & Likhomanenko, "Towards
Automatic Assessment of Self-Supervised Speech Models using Rank", ICASSP 2025**, which reports
verbatim that *"rank does not reliably predict the best-performing layer for specific downstream
tasks, **as lower-ranked layers can outperform higher-ranked ones**"*. That is a **within-method,
matched-arm, selection-rule failure** and a **published low-rank/high-information instance**, by an
author group that includes LiDAR's first author. It pre-empts both the "in-scope failure" framing and
the necessity-violation framing that earlier drafts of this paper were built on. Earlier drafts
asserted that *"none [of the prior negatives] tests the within-method, matched-arm, fixed-architecture
regime that RankMe reserves for itself"*. **That sentence was false and has been deleted.**

Also already published: Thilak et al. (LiDAR, ICLR 2024) report RankMe Spearman 0.3174 against
linear-probe accuracy on VICReg at 100 epochs and state that *"a high rank does not guarantee superior
performance"*; a RankMe co-author has co-signed the statement that *"current methods like RankMe fail
to adequately evaluate representation quality, making cross-validation without labels infeasible"*
(Otero, Mateus & Balestriero, arXiv:2410.04289) — which is itself a named *selection-rule* failure,
not merely a correlational one.

**We also do not claim more than RankMe claims.** RankMe hedges in four ways, and the hedge that
matters most to this paper is the one we tested and failed to break:

- Rank is **necessary but not sufficient**: *"a necessary (but not sufficient) condition"*. A
  high-rank representation with poor performance therefore does not contradict RankMe — RankMe
  predicts exactly that. §4.7 reports our attempt to produce the configuration that *would* contradict
  it, and its failure.
- **Same-method comparisons only**: *"RankMe should however only be used to compare different runs of
  a given method, since the embeddings' rank is not the only factor that affects performance."* This
  restriction is the one this paper takes seriously, because the reproducibility floor we measure is
  measured **inside** it: three seeds of one arm *are* different runs of a given method.
- **No monotone transfer to other datasets**; and **a named failure region**, *"Except for some
  degenerate solutions at full-rank…"*.

And **the quality-proxy claim must not be attributed to Roy & Vetterli**, who propose effective rank
as a real-valued relaxation of integer rank for signal-processing optimisation and make no claim about
representation quality anywhere. Nor to Jing et al. (ICLR 2022), whose paper is diagnostic and
mechanistic about dimensional collapse and contains no sentence claiming the singular-value spectrum
predicts downstream performance. Both attributions appear in earlier drafts on this project and in
P1 §2.6, and are corrected here and there (§2.4, §2.6).

### 1.3 The claim

Following Leavitt & Morcos (ICLR 2021), we state the claim as a property rather than as a performance
result, and we state the strongest form with its hedge attached:

> **Effective rank is unusable as a selection signal on this stack, because its between-arm
> differences are smaller than its own within-arm reproducibility floor — and this is true inside the
> same-method regime its proponents explicitly reserve for it. Its dynamic range is dominated by the
> one factor that carries no information: the training seed.**

Note what this claim is *not*. It is not "rank is uncorrelated with quality" — we did not measure that
and §4.7 points the other way. It is not "rank is wrong" — the metric is a precise measurement of the
thing it measures (§4.4). It is a claim about the **ratio of signal to nuisance in the quantity a
practitioner actually compares**, and it is the form of claim that survives whether or not rank is
right on average, because a criterion that cannot resolve its own re-measurement cannot select between
configurations regardless of what it would say if it could.

The practitioner-facing form, which is the sentence we would want quoted:

> *Before using a rank difference to select between configurations, retrain one configuration with the
> same seed several times and measure the rank spread. If the between-configuration difference does not
> exceed it, the comparison is uninformative.* On this stack that check would have disqualified **all
> seven** rank comparisons this project ever made. It says *several* repeats because the spread we
> measured is bimodal: any pair drawn from our four concordant repeats spans at most 1.028×, so a
> single retraining pair can return a floor that licenses everything.

### 1.4 Contributions

In descending order of how well evidenced they are.

1. **A measured reproducibility floor on effective rank, and the demonstration that every between-arm
   difference we have measured is inside it** (§4.1). Seven of seven differences, 1.004×–3.246×,
   against a **3.295×** floor measured from five identical same-seed retrains — a floor that is
   **bimodal**, four repeats within 2% and one at a third of them, on runs whose molecular channel
   spreads only 1.055×. We have not found either the floor or the bimodality reported.
2. **A variance decomposition separating the arm term from the seed term** (§4.2). Rank: 34.5% arm,
   65.5% seed, arm term not significant. Channel: 98.0% arm, F = 128.2. This is the contribution that
   does not depend on a sign count and it is what the paper rests on.
3. **The floor is a property of the arm, not of the statistic** (§4.3). At one training step, 6.05×
   across five seeds for one arm against 1.18× for its sibling — so a reproducibility envelope measured on
   a well-behaved configuration does not transfer to the one being diagnosed, which is the
   configuration a practitioner is looking at when they reach for rank.
4. **A defeater check** (§4.4), which the exemplar literature will demand and which earlier drafts
   lacked: the instability is in training and not in estimation, established four independent ways.
5. **The verdict is under-determined by choices nobody states** (§4.5): statistic (2/6 pairs), block
   (2/3 D2 seeds), measured view (2/6 pairs). The information verdict flips under none of them.
6. **A dose–response magnitude miscalibration, corrected against ourselves** (§4.8): −3.10% rank
   against −66.7% channel on matched preprocessing, i.e. 21.5×, replacing a previously published 3.74×
   whose preprocessing was mismatched and unstated.
7. **A negative for our own generality, preregistered and reported first** (§4.7).
8. **Three mutually incompatible statistics named `effective_rank` across ten call sites in one
   repository** (§3.1) — two of them live abort thresholds that kill training runs — reported **and**
   harmonised, with every recomputable historical instance recomputed under one definition.
9. **A fourth one, found in our own analysis code, in the section arguing the name is unreliable**
   (§4.5a). It changed a published count from 3 of 6 to 2 of 6 and withdrew a qualification from §4.7.
   It was invisible to review, to the test suite and to the authors, and surfaced **only** when the
   traceability rule was enforced by vendoring the producing code into the repository. The
   recommendation that follows is the one this paper is willing to make generally: **against this class
   of error the defence is mechanical provenance, not diligence.**

Plus a domain in which the within-method matched-arm evaluation of rank as a selection rule has not
been done: across **453 de-duplicated works** citing RankMe, Roy & Vetterli, LiDAR or α-ReQ, none
evaluates rank as a selection rule on matched arms in computational pathology or transcriptomics
(§2.2). That is where the novelty claim is placed, and it is narrow.

### 1.5 What this paper is not, and its relationship to P1

This is a methods and measurement paper. **No biological claim is made anywhere in it.** The
representation states audited here exist for other work on the same project; here they are specimens.

Its relationship to `paper/P1_CALIBRA_DRAFT.md` must be stated precisely, because the two share
evidence and would be **parallel archival submissions**. TMLR's editorial policy forbids reuse of
*results*, not merely of claims, between a submission and any paper "submitted in parallel at another
archival, peer-reviewed venue". P1 §4.11 previously carried a four-row table of rank dissociations and
P1 F11 plotted them. **Those have been removed from P1 and moved here**, per the four edits listed in
`paper/P2_FIGURES.md` §"Cross-paper deconfliction"; the edits were executed in this pass and are
recorded in §4 of `NOTEBOOK_ENTRIES/p2_rewritten_around_the_surviving_claim_20260804T1200Z.md`. P1
§4.11 is now a two-paragraph pointer with no table and no rank numbers. P1 §4.12(iv)
retains one sentence citing D2's rank values, because there it discharges an objection to P1's own
ablation rather than making a claim about rank; that sentence and §4.1 of this paper must not both be
described as a finding about effective rank.

---

## 2. Related work

### 2.1 Effective rank, and its proposal as a label-free quality criterion

**Definition — VERIFIED at full text.** Olivier Roy & Martin Vetterli, "The effective rank: a measure
of effective dimensionality", *15th European Signal Processing Conference (EUSIPCO 2007)*, Poznań,
Poland, 3–7 September 2007, IEEE, **pp. 606–610**; IEEE Xplore document **7098875**; DOI
**10.5281/zenodo.40328** (a Zenodo/EURASIP archive DOI — **no IEEE `10.1109/…` DOI exists on any
record checked, and none may be invented**). Definition 1, verbatim: *"The effective rank of the
matrix A, denoted erank(A), is defined as erank(A) = exp{H(p1, p2,…,pQ)}, where H(p1,p2,…,pQ) is the
(Shannon) entropy given by H(p1,p2,…,pQ) = − Σ p_k log p_k"*, with *"the singular value distribution
p_k = σ_k / ‖σ‖₁"* and *"all logarithms are to the base e"*. They prove *"1 ≤ erank(A) ≤ rank(A) ≤
Q"* (Property 1) and `erank(cA) = erank(A)` (Property 2).

Two consequences for this paper. **(i)** Roy & Vetterli take the SVD of the matrix as handed to them —
*"a complex-valued non-all-zero matrix A of size M × N whose singular value decomposition (SVD) is
given by A = UDV\*"* — and **no preprocessing of any kind appears anywhere in their paper**. Our
implementation **column-centres first**. That is a deliberate deviation with a stated reason (§3.1),
and **earlier drafts of this section wrongly said our implementation matched theirs exactly**. Their
treatment of near-zero singular values also differs from both ours and RankMe's: *"we adopt the
convention that 0 log 0 = 0"*, with **no ε and no tolerance**. **(ii)** Their abstract is about making
rank minimisation tractable — *"Since rank minimization is generally not practicable owing to its
integer nature, we propose a real-valued extension"* — and contains **no claim about representation
quality or downstream performance**. The closest it comes is "assess the loss incurred by
dimensionality reduction methods, such as PCA". Any text attributing a quality-proxy claim to them is
wrong.

**The quality-proxy proposal — VERIFIED at full text.** RankMe: Garrido, Balestriero, Najman & LeCun,
"RankMe: Assessing the Downstream Performance of Pretrained Self-Supervised Representations by Their
Rank", *Proceedings of the 40th International Conference on Machine Learning* (ICML 2023), PMLR 202,
pp. 10929–10974; arXiv:2210.02885 (v3, 2023-06-26 is the camera-ready). It adopts Roy & Vetterli's
definition explicitly, with an ε **outside** the division: `p_k = σ_k(Z)/‖σ(Z)‖₁ + ε`.

What it claims, verbatim — and the claim strengthens between abstract and body:

- Abstract, hedged: *"we develop a simple unsupervised criterion that is indicative of the quality of
  the learned JE-SSL representations: their effective rank"*; *"RankMe can be used for hyperparameter
  selection with nearly no reduction in final performance compared to the current selection method
  that involve a dataset's labels"*.
- Body, stronger: section headings *"RankMe Consistently Predicts Downstream performances From
  Representations"* and *"RankMe Predicts Linear Probing performance Even on Unseen Datasets"*;
  *"RankMe accurately predicts a model's performance"*.

What it **restricts**, verbatim — and this is the part a negative result must respect, and the part
this paper's claim is built inside:

- *"a necessary (but not sufficient) condition"*; *"having a high rank is a necessary condition for
  good downstream performance"*; *"This further highlights how maximal rank is only a necessary
  condition for good performance."*
- ***"RankMe should however only be used to compare different runs of a given method, since the
  embeddings' rank is not the only factor that affects performance."*** Three seeds of one arm, one
  objective, one architecture, one schedule, one cohort **are** different runs of a given method. §4.2
  is measured entirely inside this sentence.
- *"there is no inherent reason for the rank of embeddings to transfer in a monotonic way to them."*
- *"Except for some degenerate solutions at full-rank, RankMe values correlate well with
  in-distribution … and out-of-distribution … classification performance."* (Figure 1 caption.)
- A self-reported miss: *"We see than RankMe can improve OOD performance for VICReg, but leads to a
  small drop for SimCLR"* (sic).

**The ε discrepancy, unresolved and stated.** RankMe's `+ ε` sits **outside** the division (verified at
glyph coordinates in the v3 PDF: the `+ ε` glyphs sit to the right of the fraction rule, on the
fraction axis). Its `p_k` therefore sum to **`1 + min(N,K)·ε`**, not to 1, so **RankMe's statistic is
not the exponential of a Shannon entropy**. RankMe never states the ε used in that equation — the only
`10⁻⁷` in the paper belongs to the *contrasting* threshold-rank definition it argues against.
Consequently **no number in this paper is comparable to a published RankMe value, and no such
comparison is made.** The practical size of the gap is nonetheless small on our data: a faithful
RankMe implementation (ε = 1e-7, uncentred) reads 23.391 where the canonical centred statistic reads
23.387 on the same artifact, and equally close on all 12. **Everything that separates RankMe's score
from ours is column centring, not ε** — which matters, because "you evaluated a centred variant" is a
referee's best line of attack and it is answered in §4.4.

**Earlier proposal in the same family — VERIFIED.** α-ReQ: Agrawal, Mondal, Ghosh & Richards, *Advances
in Neural Information Processing Systems 35* (**NeurIPS 2022**), **pp. 17626–17638**, DOI
**10.52202/068431-1281**. *The `[COULD-NOT-VERIFY: NeurIPS 2022 venue]` marker carried by earlier
drafts is now cleared*: four independent records agree (DBLP `conf/nips/AgrawalMGR22`; the NeurIPS
proceedings page and its official BibTeX; OpenAlex W7133259885; the camera-ready PDF footer). Two
caveats travel with the citation: the PDF's own title page omits "in Self-Supervised Learning", present
only in proceedings metadata; and the NeurIPS BibTeX renders the first author "Agrawal, Kumar K" where
DBLP has "Kumar Krishna Agrawal". α-ReQ states verbatim that *"a task-agnostic measure like α is a
**necessary but not sufficient condition**"* — the same hedge as RankMe — and describes a *"Goldilocks
zone"* rather than a monotone rule. RankMe itself reports α-ReQ failing on embeddings.

### 2.2 Prior negative results — what we must not claim novelty over

**This section must be read before §4. A version of this paper that presents §4 as a discovery is not
submittable.** The census behind it is a Semantic Scholar citation-graph sweep of RankMe (159 citers),
Roy & Vetterli (581, of which 325 survived a deep-learning-context filter), LiDAR (33) and α-ReQ (15)
— **453 unique works de-duplicated by normalised title**, ~70 abstracts read in full
(`NOTEBOOK_ENTRIES/p2_prior_art_citation_graph_sweep_20260803T2326Z.md`).

#### The item that constrains us most, and it pre-empts two of our former framings

**Aldeneh, Thilak, Higuchi, Theobald & Likhomanenko, "Towards Automatic Assessment of Self-Supervised
Speech Models using Rank", ICASSP 2025.** DOI **10.1109/ICASSP49660.2025.10889651**,
arXiv:**2409.10787** (v2). **VERIFIED** independently from Crossref (title, five authors,
container-title *ICASSP 2025 - 2025 IEEE International Conference on Acoustics, Speech and Signal
Processing*, published 2025-04-06) and the arXiv Atom API.

> *"The findings indicate rank correlates with downstream performance **within encoder layers** across
> various downstream tasks and for in- and out-of-domain scenarios. **However, rank does not reliably
> predict the best-performing layer for specific downstream tasks, as lower-ranked layers can
> outperform higher-ranked ones.**"*

This hits two things earlier drafts claimed. It is a **within-method matched-arm selection-rule
failure** — one encoder, layers as the arms, rank picks the wrong one — so the claim that no prior
negative tests a within-method regime is false. And *"lower-ranked layers can outperform
higher-ranked ones"* is a **published low-rank/high-information instance**, which is precisely the
configuration our necessity test set out to find. Vimal Thilak is LiDAR's first author, so this is a
negative about rank published by the group proposing the leading replacement.

**What is still ours after this.** Aldeneh et al. vary **layer depth within one trained encoder**.
Layers of a single network are not independent training runs, and RankMe's reserved scope is
*"different runs of a given method"*. Our arms are **separately trained runs** differing in one
objective term, three seeds each, with a reproducibility floor **measured on the statistic itself**.
They report rank *correlating* with performance within layers and failing only at the argmax; we report
a variance decomposition in which the arm term is not significant at all. **The distinction is real and
it is narrow, and the paper states it rather than implying novelty.**

#### The rest of the prior negatives

| work | status | what it establishes | the regime it does *not* test |
|---|---|---|---|
| **LiDAR** — Thilak, Huang, Saremi, Dinh, Goh, Nakkiran, Susskind & Littwin, **ICLR 2024**, arXiv:2312.04000 | **VERIFIED** (content, full text); venue verified via DBLP `conf/iclr/Thilak0SDGNSL24`, OpenReview `f3g5XpL9Kb` | *"RankMe correlates poorly with downstream performance for most models"*; Spearman **0.3174**, Kendall τ **0.2056** vs LiDAR's 0.9161 / 0.8105 on VICReg at 100 epochs; *"a high rank does not guarantee superior performance"*; rank inflates early then diminishes, with *"peak downstream performance occurring far from the point of maximal rank"* | seed-level reproducibility of the statistic; it also *claims* within-method hyperparameter selection for LiDAR, so it already strains any "nobody tests in-scope" statement |
| **Otero, Mateus & Balestriero**, arXiv:2410.04289v1 | **PARTIALLY VERIFIED** (abstract; authors via arXiv Atom) | *"current methods like RankMe fail to adequately evaluate representation quality, making cross-validation without labels infeasible"* — a **named selection-rule failure**, over 250 experiments with backbone and SSL framework held fixed and class imbalance varied. Balestriero is a RankMe co-author | separately-trained matched arms; a reproducibility floor |
| **Zhang, Jiang, Gao, Willett & Maire**, "Residual Connections Harm Generative Representation Learning", arXiv:2404.10947v5 | **`[ABSTRACT-LEVEL ONLY — the sign must be confirmed in the body before this is cited as a necessity violation]`** | *"we find **correlation between low effective feature rank and downstream task performance**"*, matched-arm by construction (one MAE / ViT-B/16 scalar varied), better arm reported as the lower-rank one (kNN 27.4% → 63.9%) | — |
| **Kulkarni, Springer, Subramonian & Swayamdipta**, arXiv:2602.20433v2 | **PARTIALLY VERIFIED** (abstract) | *"While the best-performing models often exhibit a high effective rank, this trend is not universal… none can reliably predict downstream performance."* | LLM unembedding geometry, not joint-embedding selection |
| **Cheng**, arXiv:2607.13432v1 | **PARTIALLY VERIFIED** (abstract) | existing measures including effective rank *"lack theoretical grounding and correlate poorly with performance on new tasks"* | continual-learning plasticity |
| **Li, Agrawal, Ghosh, Teru, Santoro, Lajoie & Richards**, arXiv:2509.23024 | `[S2/abstract record only]` | measures **both** RankMe and α-ReQ across pretraining and post-training; *"consistent **non-monotonic** sequence of three geometric phases"*, with a compression phase that **reduces** dimensionality while *"marked with significant improvement in downstream task performance"* | **pre-empts the non-monotonicity limb** within a single run |
| **Jamali, Cheng & Vargas-Hernández**, arXiv:2510.14217 | `[S2/abstract record only]` | *"richer spectral features do not consistently yield better generalization performance, contradicting common representation heuristics used in self-supervised learning"*, with *"consistently negative correlations"* for local 3D representations | **the closest published relative in a molecular-science domain**; cites RankMe and α-ReQ |
| **Kokilepersaud, Prabhushankar & Al-Regib**, "AdaDim", arXiv:2505.12576 | `[S2/abstract record only]` | *"the best performing SSL models do not have the highest H(R) nor the lowest I(R;Z)"* | — |
| **Yusupov et al.**, arXiv:2509.25359 | `[S2/abstract record only]` | geometric metrics' *"apparent discriminative power **collapses once length is controlled**"* — our confound/matched-arm argument, executed for LLMs | — |
| **Adilova, Petzka, Fischer & Geiger**, arXiv:2606.21593 | `[S2/abstract record only]` | *"low MI does not reliably correspond to geometric compression … a negative and nonlinear relationship **that can reverse when varying training setup**"* | a sign reversal under a controlled training change — the same shape of claim as ours, one level up |

**Cite defensively:** Zaiem, Kemiche, Parcollet, Essid & Ravanelli (Interspeech 2023, arXiv:2306.00452,
DBLP `conf/interspeech/ZaiemKPER23`; and *Computer Speech and Language*, arXiv:2308.14456) report that
*"altering the downstream architecture structure leads to significant fluctuations in the performance
ranking of the evaluated models"*. That attacks the **ground truth** any rank-vs-downstream comparison
is measured against, **including ours**, and must be conceded in §6 rather than omitted.

#### Where this leaves the novelty claim

Re-scoped, and narrower than earlier drafts:

> **A within-method, matched-arm evaluation of rank as a selection rule in computational pathology /
> transcriptomics, with a measured seed-level reproducibility floor on the statistic and a variance
> decomposition separating the arm term from the seed term.**

Across all 453 de-duplicated citers, **none** evaluates rank as a selection rule on matched arms in
pathology or transcriptomics, and none isolates the seed term from the arm term as the explanatory
claim. Computational pathology and omics representation learning are otherwise live in this corpus
(Jaume et al., CVPR 2024 DOI 10.1109/CVPR52733.2024.00920 and ECCV 2024 arXiv:2408.02859; HEST-1k,
NeurIPS 2024 arXiv:2406.16192; DenAdel et al., **Nature Methods 2026**, DOI 10.1038/s41592-026-03120-y;
kaiko.ai arXiv:2404.15217; Gustafsson et al. arXiv:2410.00945 on WSI→gene expression, closest to our
setup).

`[SEARCH RISK — DOWNGRADED, NOT CLOSED.]` The citation-graph leg is now closed for RankMe, Roy &
Vetterli, LiDAR and α-ReQ. It remains **abstract-level only**: no full texts were read in triage, and
**12 of 159** RankMe citers had no abstract in the API response and were triaged on title alone. LiDAR
was nearly missed in exactly that way by the previous OpenAlex-only sweep, which saw roughly **8%** of
RankMe's citation graph. One item must be read at full text before submission: **Luisto et al.,
arXiv:2604.14815** (Finnish histopathological reports; *"predict the benefit of domain-specific
pre-training … from observing the **geometry of embedding changes**"*) — pathology plus a
rank-geometry selection question, small and unrefereed, and the one item in the sweep that could scoop
the framing.

### 2.3 The strongest defences of rank, and what we say to each

Earlier drafts cited **none** of these. A negative paper that engages only the weakest form of the
opposing case is not defensible, so the three strongest are stated in their own words first.

**C1 — the design we claim breaks, reported working, with a power law.** Deng, Sun, Dou & Xu, "Unify
Variables in Neural Scaling Laws for General Audio Representations via Embedding Effective Rank",
arXiv:**2510.10948**v1, published 2025-10-13. **VERIFIED** from the arXiv Atom API (title, four
authors, date). `[UNVERIFIED: peer-reviewed venue — none found; cite as an arXiv preprint.]`

> *"we present a systematic study of scaling laws for general audio representations by utilizing
> **embedding effective rank (RankMe) as a unifying metric** … allowing us to examine scaling
> behaviors across a wide hyper-parameter space, including **model size, training data volume,
> computational budget, architectural configurations**, etc. Our empirical findings reveal a
> **consistent power-law relationship between RankMe and representation quality**, suggesting that
> embedding effective rank serves as a **reliable proxy** for assessing and predicting model
> performance."*

**Our answer, in order of strength.** (i) The sweep is over **capacity-like** variables — model size,
data volume, compute — which move rank and quality together for reasons that have nothing to do with
rank being a quality proxy. Our arms are matched **on capacity** and differ only in an objective term
or a supervision target; that is the case that separates the two accounts, and it is the case a
practitioner faces when choosing between equally-sized configurations. (ii) They report **no
seed-level reproducibility floor** for RankMe, and our §4.2 measurement is that the seed term dominates
the arm term; a power law fitted across a capacity range far larger than the seed spread is not
evidence that rank resolves differences *inside* that spread. (iii) Audio, not pathology. We do not
claim their result is wrong, and we do not claim rank fails at capacity scale — **we have not measured
that**, and §6 records it as unmeasured.

**C2 — a peer-reviewed RankMe-works result inside our own domain.** Awasthi, Mend Mend Arachchige &
Zhu, "Unsupervised evaluation of pre-trained DNA language model embeddings", ***BMC Genomics*** **26**,
issue 1, article **710** (2025); DOI **10.1186/s12864-025-11913-2**, published 2025-08-01. **VERIFIED**
from Crossref (title, three authors, journal, volume, issue, article number, date). PMID 40751178 and
PMC12315385 are reported by the sweep and were **not** independently verified.

> *"We propose a framework to evaluate DLM embeddings using unsupervised numerical linear
> algebra-based metrics **RankMe, NESum, and StableRank** … we observed a **positive correlation
> between unsupervised metrics and supervised performance, supporting the utility of unsupervised
> metrics as effective proxies for model quality assessment**."*

**This is the single most awkward citation for a paper claiming rank fails on molecular
representations, and it cannot be left out.** Our answer is narrow and stated as such: it compares
**six different pre-trained models** — a **cross-method** comparison, which is exactly the regime
RankMe itself disclaims (*"should however only be used to compare different runs of a given method"*).
Our comparison is within-method, and our finding is that the same-method regime is where the seed
nuisance eats the signal. **The two results are compatible**, and a reader should take from the pair
that rank may separate genuinely different model families while failing to separate runs of one. We
have not reimplemented their pipeline and make no claim about their numbers.

**C3 — a matched-arm design in which rank tracks transfer, with a proof.** Ruan, Zhang, Wang & Zhang,
"Muon Learns More Robust and Transferable Features than Adam", arXiv:**2606.09658** (S2
CorpusId 289099544). `[NOT INDEPENDENTLY VERIFIED — Semantic Scholar record only; verify before
submission.]` Reports *"we **prove** that Muon attains larger margins and **higher effective rank**
than Adam and SGD"*, with the same architecture and the optimiser varied. Our answer: an optimiser
change is arguably a method change; and, more usefully, **their design is the one our §4.1 rule asks
for a control on** — a proof about expectations says nothing about whether a single pair of runs is
resolvable.

**Further defences, cited and not rebutted at length** — all `[S2 record only, not independently
verified]` unless marked: **Zhang, Deidda, Higham & Tudisco**, arXiv:2502.04591 (*"rank-based metrics
consistently capture oversmoothing, whereas energy-based metrics often fail"*); **Zhuo, Wang, Ma &
Wang**, ICLR 2023, arXiv:2303.02387 — the strongest *theoretical* case that rank is causally tied to
representation quality (*"This rank difference will **provably lead to an improvement of effective
dimensionality**"*); **Kim, Kokilepersaud, Prabhushankar & AlRegib**, **WACV 2026**, DOI
10.1109/WACV61042.2026.00461, arXiv:2511.06450 — peer-reviewed, rank used prescriptively in a
**multi-modal fusion** setting, which is also ours; **Sun, Lin, Zhang, Duan & Liu**, **IEEE TIP 2026**,
DOI 10.1109/TIP.2026.3682105 — effective rank as diagnostic *and* training objective, peer-reviewed
journal; **Billa**, arXiv:2602.15997 (*"only rankme reliably precedes capability acquisition for hard
tasks"*); **Gupta**, arXiv:2606.24903 — effective rank as a stopping rule that works (ρ_pool = 0.6366,
p = 2.9×10⁻⁵⁷; AUC 0.787), `[unrefereed preprint]`. And **Dai, Xu, Wen, Liu & Huang**,
arXiv:2510.17299, is the closest structural analogue to our contribution outside speech: a within-
method **checkpoint-selection** rule built from *"a class-relevance measure **and** an effective
dimensionality measure"*, +3.0 mIoU — a design that **concedes dimensionality alone is insufficient**.

**The honest summary of §2.3.** Rank works in several published settings and fails in several others,
and the published cases where it works are predominantly cross-method, cross-capacity or
cross-checkpoint. This paper adds one measurement that discriminates between those regimes rather than
a verdict on the metric.

### 2.4 Dimensional collapse and anti-collapse regularisation

**VERIFIED at full text.** Jing, Vincent, LeCun & Tian, "Understanding Dimensional Collapse in
Contrastive Self-Supervised Learning", ICLR 2022, arXiv:2110.09348. Verbatim: *"dimensional collapse,
whereby the embedding vectors end up spanning a lower-dimensional subspace instead of the entire
available embedding space"*, diagnosed from the spectrum — *"A number of singular values drop to zero,
indicating collapsed dimensions."*

**Correction to earlier drafts on this project and to P1.** A full-text search of Jing et al. for
sentences relating the spectrum to downstream performance returns none. The paper is diagnostic and
mechanistic — its own conclusion identifies "two mechanisms causing dimensional collapse: strong
augmentation and implicit regularization" — and **contains no claim that rank predicts quality.** P1
§2.6 grouped it among proposals of "geometric proxies for representation quality"; that grouping was
inaccurate and **has been corrected in P1 in this pass**.

**VERIFIED.** VICReg: Bardes, Ponce & LeCun, ICLR 2022, arXiv:2105.04906. Camera-ready abstract,
verbatim: *"a method that explicitly avoids the collapse problem with two regularizations terms applied
to both embeddings separately: (1) a term that maintains the variance of each embedding dimension
above a threshold, (2) a term that decorrelates each pair of variables"* (sic). **Note:** the arXiv
metadata abstract and the ICLR camera-ready abstract differ; quote whichever is cited and do not mix
them. VICReg's geometric terms are claimed to *prevent collapse*, not to *indicate downstream quality*.

**VERIFIED.** Wang & Isola, ICML 2020, PMLR 119, arXiv:2005.10242: *"Extensive experiments on standard
vision and language datasets confirm the strong agreement between both metrics and downstream task
performance."* Correlational and empirical; not predictive and not causal. Note the adjacent
precedent that this paper must differentiate itself from: Fang, Li, Sun & Wang, "Rethinking the
Uniformity Metric in Self-Supervised Learning", **ICLR 2024**, arXiv:2403.00642 (verified via DBLP;
**full text NOT RETRIEVED**) makes the same move one metric over, on *axiomatic* grounds. Ours is
empirical and about **reproducibility**, which is a different and independent objection. `[Read Fang
et al. at full text before submission — it is the paper a referee will say this one duplicates.]`

`[UNVERIFIED]` Barlow Twins (Zbontar et al. 2021) — PMLR URL recorded in this repository, no quote
retrieved. `[UNVERIFIED]` LDReg (Huang et al., arXiv:2401.10474) — *local* dimensionality, not the same
construct. In pathology specifically, the Robustness Index (de Jong et al., 2025) and
representational-similarity analysis (Mishra & Lotter, 2025) are proposed as confound-aware
representation diagnostics; both were spot-check verified in a 2026-07-29 literature audit and not
re-verified in this pass.

§4.9 reports what happened when we followed practice (2): adding a covariance-decorrelation term
raised effective rank by 107% while moving the then-benchmark statistic by 0.0001, and a later,
better-sourced measurement found that **our implementation** of the term has a **collapsed global
minimum** and self-extinguishes at every weight tested — **in the absence of a momentum key encoder**,
which is a condition we did not know we were assuming until §4.9a measured the same term with one and
found it raising rank monotonically instead. We state all of it as an observation about our
implementation *in one queue configuration*, not as a refutation of the published methods, which we
have not reimplemented faithfully or benchmarked.

### 2.5 The published alternatives, computed on our own artifacts

A referee is entitled to ask why a paper reporting rank's failure did not evaluate the published
alternatives. **The `[NOT MEASURED]` marker earlier drafts carried here is now closed.** LiDAR, α-ReQ,
RankMe-as-published, participation ratio and stable rank were all computed on the same 12 frozen
artifacts against the same ground truth (§4.6, and
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md`).

**On our data there is no better metric.** LiDAR — the strongest published alternative and the one
that most constrains our claim — scores **0/3 on the in-scope D2 contrast**, choosing the
information-poorer arm every time, at **every δ across eight orders of magnitude**. That is not a
gotcha: LiDAR's own limitations section reports *"instances where **the LiDAR metric exhibits a
negative correlation with probe accuracy**, particularly pronounced in scenarios like VICReg when
dealing with higher dimensional embeddings"*, so our observation is a documented failure mode of the
metric rather than an artifact of our adaptation. LiDAR also supplies a verbatim admission we use in
§4.8: *"This illustrates that high rank is a necessary but not a sufficient condition for high
performance."*

**Fidelity, and the one genuine gap.** Our LiDAR implementation matches arXiv:2312.04000v1 §3 Eqs.
(1)–(4) exactly. The modality-as-view adaptation is **explicitly licensed** by the paper's footnote 4
(*"Data augmentations, or otherwise data points which are treated as positive samples"*), which our
patient-paired objective satisfies; their `n > p` recommendation is met (n = 2,766 > p = 256). **The
gap is q**: they use q = 50 or q = 10, sweep q only over 10–100, and never discuss q = 2 or state a
minimum; we have **q = 2**, the bare minimum at which Σ_w is estimable. That is disclosed as the one
respect in which our LiDAR is outside the authors' tested regime. α-ReQ was re-implemented from the
authors' released `fastssl` code (`stringer_get_powerlaw`, weighted least squares over eigenvalue
ranks 4–100, weights 1/rank) rather than from the paper text, which is **not sufficient to reproduce a
number**; the index range must **not** be attributed to the paper. All 12 of our artifacts have
α between 2.6 and 4.8, far outside α-ReQ's "Goldilocks zone", so on its own terms it declines to
distinguish them and the `|α − 1|` rule we scored is **ours**, not theirs.

`[STILL NOT MEASURED]` — a labelled linear probe on every artifact, which is the ground truth LiDAR and
RankMe were validated against. Ours is a held-out canonical correlation against unsupervised molecular
targets (§3.2), which is a different reference standard and is the one Zaiem et al. would attack.

### 2.6 Reference verification status

Three fabricated citations have previously contaminated this project (`HANDOFF_BUILD_AGENT.md:156`).
The verification protocol is recorded in `NOTEBOOK_ENTRIES/winkler_prior_art_20260803T0120Z.md`:
retrieve the full text (not the abstract), make targeted passes over named sections plus every figure
caption and table header, tabulate what the paper *actually* reports, quote verbatim the closest it
comes to the claim being attributed to it, and state where it falls short.

| reference | status | retrieval | note |
|---|---|---|---|
| Roy & Vetterli, EUSIPCO 2007, pp. 606–610, IEEE Xplore 7098875, DOI 10.5281/zenodo.40328 | **VERIFIED** | full-text PDF + DBLP + OpenAlex | makes **no** quality claim; **no IEEE DOI exists** — do not invent one; `spectral.py` **deviates** by centring |
| RankMe — Garrido, Balestriero, Najman & LeCun, ICML 2023, PMLR 202:10929–10974, arXiv:2210.02885v3 | **VERIFIED** | full-text PDF + DBLP | claims and restrictions quoted in §2.1; ε **outside** the division |
| **Aldeneh, Thilak, Higuchi, Theobald & Likhomanenko, ICASSP 2025, DOI 10.1109/ICASSP49660.2025.10889651, arXiv:2409.10787v2** | **VERIFIED** | Crossref + arXiv Atom (abstract) | **leads §2.2**; body not read — `[full text not retrieved]` |
| LiDAR — Thilak et al., **ICLR 2024**, arXiv:2312.04000v1 | **VERIFIED** (content + venue) | full-text PDF; venue via DBLP `conf/iclr/Thilak0SDGNSL24`, OpenReview `f3g5XpL9Kb` | **earlier drafts' `[UNVERIFIED: venue] / cite as preprint` is now wrong and has been corrected.** ICLR camera-ready **not** retrieved (OpenReview bot challenge); all quotes are from arXiv v1. The paper is internally inconsistent about VICReg's n (10k main text vs 5,000 in Appendix 11.1) — do not cite a single n |
| α-ReQ — Agrawal, Mondal, Ghosh & Richards, **NeurIPS 2022**, pp. 17626–17638, DOI 10.52202/068431-1281 | **VERIFIED** | DBLP + NeurIPS proceedings BibTeX + OpenAlex + camera-ready PDF | **`[COULD-NOT-VERIFY: venue]` is now cleared.** Title and first-author-name discrepancies noted in §2.1 |
| **Deng, Sun, Dou & Xu, arXiv:2510.10948v1** | **VERIFIED** (metadata + abstract) | arXiv Atom | `[UNVERIFIED: peer-reviewed venue — none found]` |
| **Awasthi, Mend Mend Arachchige & Zhu, *BMC Genomics* 26(1):710 (2025), DOI 10.1186/s12864-025-11913-2** | **VERIFIED** | Crossref | PMID/PMCID reported by the sweep, **not** independently verified |
| Otero, Mateus & Balestriero, arXiv:2410.04289v1 | **PARTIALLY VERIFIED** | arXiv Atom (abstract) | body not read |
| Zhang, Jiang, Gao, Willett & Maire, arXiv:2404.10947v5 | **PARTIALLY VERIFIED** | arXiv Atom (abstract) | **`[ABSTRACT-LEVEL ONLY — sign must be confirmed in the body]`** |
| Ruan, Zhang, Wang & Zhang, arXiv:2606.09658 | `[NOT INDEPENDENTLY VERIFIED]` | S2 record only | verify before submission |
| Kulkarni et al. 2602.20433v2; Cheng 2607.13432v1 | **PARTIALLY VERIFIED** | abstracts | LLM / plasticity domains |
| Zhuo et al. ICLR 2023 2303.02387; Kim et al. WACV 2026 10.1109/WACV61042.2026.00461; Sun et al. IEEE TIP 2026 10.1109/TIP.2026.3682105; Zhang/Deidda/Higham/Tudisco 2502.04591; Billa 2602.15997; Gupta 2606.24903; Dai et al. 2510.17299; Li et al. 2509.23024; Jamali et al. 2510.14217; Kokilepersaud et al. 2505.12576; Yusupov et al. 2509.25359; Adilova et al. 2606.21593 | `[S2 RECORD ONLY]` | Semantic Scholar | **must be verified against Crossref/arXiv before submission**; none is load-bearing for any number |
| Zaiem et al., Interspeech 2023 arXiv:2306.00452 (DBLP `conf/interspeech/ZaiemKPER23`) and arXiv:2308.14456 | **PARTIALLY VERIFIED** | DBLP + abstract | cited **against** us in §6 |
| Jing, Vincent, LeCun & Tian, ICLR 2022, arXiv:2110.09348 | **VERIFIED** | full-text PDF | contains **no** rank→performance claim |
| MoCo — He, Fan, Wu, Xie & Girshick, CVPR 2020, arXiv:1911.05722 | **VERIFIED** | full-text PDF | see §5.2 for the corrections it forced |
| VICReg — Bardes, Ponce & LeCun, ICLR 2022, arXiv:2105.04906 | **VERIFIED** | full-text PDF + arXiv API | arXiv abstract ≠ camera-ready abstract |
| Wang & Isola, ICML 2020, PMLR 119, arXiv:2005.10242 | **VERIFIED** | full-text PDF + arXiv API | "strong agreement", correlational |
| Fang, Li, Sun & Wang, ICLR 2024, arXiv:2403.00642 | **VERIFIED** (bibliographic) | DBLP | **full text NOT RETRIEVED**; nearest neighbour to this paper's genre |
| Leavitt & Morcos, ICLR 2021, arXiv:2003.01262 | **VERIFIED** | OpenReview API + arXiv API | structural exemplar (§1.3, §4.4); the arXiv record carries no venue field |
| Barlow Twins — Zbontar et al. 2021 | `[UNVERIFIED]` | — | PMLR URL only |
| LDReg — Huang et al., arXiv:2401.10474 | `[UNVERIFIED]` | — | *local* dimensionality; different construct |
| de Jong et al. 2025; Mishra & Lotter 2025 | spot-check verified 2026-07-29 | — | not re-verified in this pass |
| Prior-art census for §2.2 | **DOWNGRADED FROM INCOMPLETE** | S2 citation graph (453 works) + OpenAlex + Crossref + arXiv + DBLP | abstract-level only; 12/159 RankMe citers had no abstract; Luisto et al. arXiv:2604.14815 must be read at full text |

**Could not retrieve, recorded rather than worked around:** Semantic Scholar `/paper/search` (HTTP 429
throughout); OpenReview forum pages for LiDAR and α-ReQ (bot challenge); numeric δ and ε for LiDAR
(absent from the paper); any IEEE `10.1109/…` DOI for Roy & Vetterli (none exists on the records
checked). **OpenAlex is unusable as a citation source for RankMe**: its two RankMe records return 12
unique citers combined against Semantic Scholar's 159, most likely because PMLR does not mint Crossref
DOIs.

---

## 3. Methods

### 3.1 The canonical definition, three statistics under one name, and four unstated degrees of freedom

#### The canonical definition, stated explicitly

> **Canonical effective rank = Roy & Vetterli 2007, Definition 1, order 1, computed on the
> column-centred matrix, with rows left at their own norms, and singular values cut at the standard
> LAPACK relative tolerance `σ > σ_max · max(n, p) · eps`.**
> `v2/calibra/spectral.py`, `CANONICAL.label == "centred|order1"`.

Every rank number in §4 is this statistic unless the row says otherwise, and every row also names the
**block** it was computed on (raw or confound-residualised), because that choice is worth more than the
choice of statistic.

Four choices, each stated because each changes the number:

1. **Order 1, not order 2.** R1 is the published statistic and the only one comparable to anything
   outside this repository.
2. **Centred — a deliberate deviation from both source papers, stated at every table.** Neither Roy &
   Vetterli nor RankMe centres (§2.1). We do, because *uncentred* effective rank is not a property of
   the representation's spread at all: it is a function of the column mean's magnitude relative to that
   spread, and it errs in **both** directions. A large shared offset drives the uncentred value of an
   isotropic, near-full-rank representation to ~1 — reading as total collapse when nothing has
   collapsed; and on the collapse family `zᵢ = m + aᵢ·u` documented on this project
   (`NOTEBOOK_ENTRIES/g26_centring_fix_20260803T0730Z.md`), with `m` comparable to the spread, it reads
   ~2 where there is exactly one direction of variation. The centred value is exactly invariant to a
   shift. All four statements are pinned by
   `test_effective_rank_canonical.py::test_uncentred_rank_is_a_function_of_the_mean_and_centred_rank_is_not`.
   §4.4 shows the choice does not carry the result.
3. **No row normalisation.** It appears in neither source paper, it discards norm variation that is
   part of the representation, and it is one of the choices that flips a verdict (§4.5).
4. **A relative, not absolute, singular-value cut.** The implementation formerly filtered at an
   **absolute** `> 1e-12`, which breaks the scale invariance Roy & Vetterli prove (Property 2): scaling
   a matrix by 1e-9 emptied the spectrum and returned 0. A separate *scale-relative* degeneracy floor at
   1e-12 of the input's own norm keeps a representation collapsed to within float noise scoring 0.
   **This change moved no number in this paper**: the maximum relative difference between the two cuts
   over all 68 recomputed artifact × block combinations is `0.000e+00`.

#### And how it differs from RankMe's

RankMe's `p_k = σ_k(Z)/‖σ(Z)‖₁ + ε` places the ε **outside** the division, so its probabilities sum to
**`1 + min(N,K)·ε`** and **its statistic is not the exponential of a Shannon entropy** (§2.1). Roy &
Vetterli use no ε and the `0 log 0 = 0` convention. **No number in this paper is comparable to a
published RankMe value.**

#### Three statistics under one name, across ten call sites

A full AST + SVD scan of the tree at the start of this work found **ten call sites carrying three
mutually incompatible statistics, all named `effective_rank`**. Two of the ten are **live abort
thresholds that kill training runs**.

| # | definition | implementation sites (pre-`85c0fa8`) | used for |
|---|---|---|---|
| **R1** | Roy & Vetterli, order 1, column-centred: `exp(−Σ pᵢ ln pᵢ)`, `p = σ/Σσ` | `v2/calibra/spectral.py:14` (declared "single source of truth"); torch duplicates at `v2/run_rank_ablation.py:35`, `v2/tests/test_stress_collapse.py:23`, `v2/calibra/e0_basis_transfer.py:432`; a fifth inline at `e0_basis_transfer.py:480` with a *different* (float32-eps) tolerance | every CALIBRA readout: §4.1, §4.8, §4.9 |
| **R2** | order-2 participation ratio of the centred singular values, `(Σσ)²/Σσ²` | `v2/research/rebase/d1_audit.py:149` | D1 audit check A5 — the statistic an earlier draft nominated for the D1 table |
| **R3** | R2 **after L2-normalising each patient row** | `v2/research/rebase/d1_geometry_probe.py:50`; **`v2/training.py:569` (the in-run rank tripwire, `--rank-tripwire-minimum 4.0`)**; **`v2/runner.py:942` (the gate probe, same bar)**; `v2/research/rebase/d1_collapse_causal_test.py:75` | all live-checkpoint geometry probes; the momentum sweep of §5.2; both admission gates |

Additionally, the historical instance discussed in §4.9 is reported in its source not as an effective
rank at all but as "`z_biology` matrix rank" — a **hard numerical rank**, maximal at the batch size of
16, and a fourth thing again.

**They are not close, and the difference has a fixed sign.** `(Σσ)²/Σσ² = 1/Σpᵢ²` is exactly the
order-2 Hill number of the same singular-value distribution, of which effective rank is the order-1
Hill number. Hill numbers are non-increasing in the order, so **R2 ≤ R1 and R3 ≤ R1 for every matrix**,
with equality only on a flat spectrum (asserted in
`test_effective_rank_canonical.py::test_order_two_is_never_above_order_one`). Measured over all 68
recomputed artifact × block combinations:

| ratio | min | median | max |
|---|---:|---:|---:|
| R2 / R1 | 0.338 | 0.629 | 0.813 |
| R3 / R1 | 0.351 | 0.655 | 0.826 |
| R1 uncentred / R1 | 0.116 | **0.995** | 1.000 |

*Provenance: `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`
§4; outputs `~/ws_rank/RANK_RECOMPUTE.json`, `~/ws_rank/RANK_RECOMPUTE_P1B.json`; scripts vendored into
the repository at commit `8609081`.*

**A rank quoted without naming its statistic is uncertain by up to a factor of three.** Note the third
row: centring, the choice a referee is most likely to attack, is the one that matters *least* at the
median — a fact §4.4 uses.

#### The fourth degree of freedom: which matrix

Beyond centring, row normalisation and the order there is **which matrix the SVD is taken of**. Every
CALIBRA readout computes `effective_rank(x)` on the **raw** representation block (`run_calibra.py:140`)
while the channel it is compared against is a `heldout_top_cca` on the **confound-residualised** block.
`d2_readout.py` reports both. Earlier drafts quoted the residualised value for one instance and the raw
value for two others **and said neither**. The choice is not cosmetic: it flips the R3 arm ordering in
two of three D2 seeds (§4.5) and changes the dilution headline magnitude by a factor of **5.8** (§4.8).
Note also that on the **raw** exported artifacts `R2 ≡ R3` exactly, because the model already
L2-normalises `z_biology`; the R2/R3 distinction exists only after residualisation, which is why it went
unnoticed.

#### One implementation, and a test that keeps it

There is now **one implementation**, imported by all ten sites.
`v2/tests/test_effective_rank_canonical.py` (13 tests) pins it against hand-computed values on a matrix
of known spectrum (`σ ∝ (2,1,1)` ⟹ `erank = 2√2 = 2.8284271247461903`, order-2 `= 8/3`), checked both
uncentred on `diag(2,1,1)` and through a non-zero column mean so that centring is simultaneously proven
applied; against Roy & Vetterli's Property 1 including both equality cases; against scale invariance;
against torch ≡ numpy on all six variants; asserts **object identity** across every importable call
site; and **fails on an AST + SVD scan of the source tree if a second definition, or any unallowlisted
SVD-based rank, reappears.** The two live abort thresholds retain **R3**, named explicitly at the call
site because their 4.0 bar was calibrated against R3 readings (6.92–7.25 healthy, 1.46–1.98 collapsed),
and both now log the canonical value beside it (`biology_effective_rank_canonical`) so the bar can be
recalibrated on evidence rather than by assumption. Suite: **317 passed** with thread caps (304 before;
+13).

We report the original state as a finding rather than an embarrassment because it is evidence for the
paper's own thesis — a scalar whose name is stable while its definition is not gets quoted across
contexts it does not survive — and because it took reading three implementations side by side, while
writing this paper, to notice.

**We are not exempt from the practice this paper criticises.**
`v2/calibra/e1_rank_information.py` is a preregistered, gate-enforced, three-seed experiment in this
repository whose aggregated endpoints include `delta_effective_rank`, and whose aggregator asserts
`rank_positive = (frame.delta_effective_rank > 0).all()` (`v2/calibra/aggregate_e1.py:38`). It further
computes an `information_density` defined as `direction_count_above_floor / effective_rank`
(`e1_rank_information.py:298`) — a quantity with effective rank in its denominator. Its own docstring
already hedges: *"Rank alone is intentionally not an endpoint."* **E1 has never been run**; no `E1_*`
rows appear in `v2/research/rebase/nature/GATE_LOG.md` and no E1 outputs exist under `runs/`.

### 3.2 The information measures, and their measured nulls

No result here compares rank to a bare correlation. Every information measure has an explicit chance
level that was *measured*, not assumed.

| measure | definition | chance level | where |
|---|---|---|---|
| **Held-out top canonical correlation ("the channel")** | canonical directions fit on one half of the held-out patients and scored on the other; 16 components per side; both blocks cross-fitted-residualised against a cancer + pooled-tissue-source-site design first | permutation null median **0.145–0.147** (dilution sweep, 2,766 patients); **0.140** (D2, 200 draws, `permutation_p` floor 0.005) | `v2/calibra/spectral.py:78-108`; `v2/calibra/run_calibra.py` |
| **Null-corrected channel ratio** | `(observed − null median) / (observed₀ − null median₀)` | 0 by construction | `DILUTION_LOWER_BOUND.md` §2 |
| **Within-cancer specificity** | the benchmark statistic of an earlier codebase generation; **definition not recoverable from any file now in this repository** | not recorded | §4.9, §6.4 |
| **Retrieval accuracy @1** | cross-modal nearest-neighbour retrieval within a fixed 16-patient batch | **0.062** = 1/16 | `NOTEBOOK.md` 2026-08-02 01:20 UTC |
| **In-batch InfoNCE** | symmetric cross-modal InfoNCE, temperature 0.07 | **ln 16 = 2.7726** (16-patient batch, frozen queue excluded from the candidate count); **ln 80 = 4.38** (earlier configuration, live 64-key queue contributing candidates); **ln 2576 = 7.854** and **ln 4310 = 8.369** at training scale | `v2/losses.py:13`; `NOTEBOOK.md:1554` |

**These chance levels must not be mixed.** `ln 80 = 4.38` belongs to the pre-fix gate configuration and
`ln 16 = 2.7726` to the post-fix one. Every InfoNCE number in this paper is quoted with its own chance
level attached.

**A null that is not zero.** The permutation null for a 16-component canonical correlation on ~2,700
patients sits between **0.140 and 0.147 depending on the cohort and the permutation procedure** (see
the footnote below — the two ends are different procedures and must not be substituted), because a
multivariate maximum over 16 fitted directions is upward-biased
at finite *n*. A raw channel ratio therefore flatters the surviving signal and the null-corrected column
is the one quoted throughout. This is also why "both arms at chance" is not a testable statement for
these readouts, and why §4.7's negative control is stated as "controls must score below real targets" —
necessary and not sufficient
(`NOTEBOOK_ENTRIES/d1_readout_preregistration_20260803T1700Z.md:59-63`).

> **Footnote — this project quotes at least three permutation nulls and they are not substitutable.**
> **0.140** is D2's, from a **200-draw row-shuffle of the residualised target matrix**
> (`NOTEBOOK_ENTRIES/D2_stratified_result_20260803T1210Z.md`); it is the comparator for every D1 and D2
> number in this paper, including §4.7.3's random controls. **0.145–0.147** is the dilution sweep's,
> from a **300-draw within-cancer permutation** at a different *n* (`DILUTION_LOWER_BOUND.md` §2); it is
> the comparator for §4.8 and nothing else. **0.151–0.158** appears elsewhere on the project at other
> cohort sizes and component counts. They differ in *n*, in component count and — for 0.140 — in the
> **permutation procedure itself**, and P1's audit records explicitly that the median is cohort- and
> capacity-specific and **must not be carried across**
> (`NOTEBOOK_ENTRIES/p1_submission_draft_20260803T1230Z.md` §5).
>
> An earlier version of §4.7.3 compared D2's random controls against **0.147**. That was wrong; it is
> corrected to **0.140**, and the conclusion is unchanged (0.44/0.147 = 3.0× against 0.44/0.140 = 3.2×).
> **We record the error and not only the fix, because it is the same failure mode as §3.1 and §4.5(a) —
> numbers that look interchangeable, are not, and carry no label saying which is which.** A related
> label conflict survives in our own notes:
> `NOTEBOOK_ENTRIES/d1_a3_verdict_and_effect_vs_floor_20260804T0300Z.md` calls 0.140 a "within-cancer"
> null where `p1_submission_draft` identifies it as a row-shuffle. **The value is agreed; the procedure
> label is not.** We follow the latter and flag the former rather than quietly picking one.

**In-sample bias, checked.** The headline channel is fitted at a 16-component budget and is therefore
upward-biased. Re-run with `heldout_top_cca` (directions fit on half the patients, scored on the other),
**every arm ordering survives on all three co-trained views, with equal or larger deltas** (e.g. D2 seed
42 `wsi_biology`: +0.1325 in-sample → +0.1541 held-out). The ground truth is not an in-sample artifact.

### 3.3 Cohorts, representations and artifacts

All measurements are TCGA. Two cohort configurations appear and are labelled at every table.

- **Maximal paired split — 6,427 patients** (3,118 train / 543 validation / 2,766 test), holding out
  whole cancers: 11 development cancers, 21 held out. Used by §4.1–§4.8.
- **Earlier configuration — 2,530 held-out patients over 21 test cancers.** Used by §4.9's instance 2.

Representations are one of:

- **Zero-parameter patch statistics.** `concat(mean, std)` over frozen H-Optimus-0 patch tokens
  (1,536-d), PCA-reduced to 256 dimensions refit on train rows only per level, retained variance
  0.879–0.923. Used by §4.8 — it has *no fitted parameters*, so neither the rank change nor the channel
  change can be attributed to different training runs.
- **Trained `wsi_biology` states,** 256 output features, exported to `.npz` diagnostic artifacts with a
  `split` column; test rows only. Each artifact additionally stores two further co-trained views,
  `rna_biology` and `full_biology`, which §4.5 exploits.

**The 12 artifacts that carry §4.1–§4.7.** Four arms × three seeds, all frozen:

| experiment | arms | differ in | artifacts |
|---|---|---|---|
| **D2** | **H** = Hallmark pathway supervision, **I** = 128 perturbation-basis coordinates | the supervision **target table** only; both arms `--objective-profile programme_only` | `~/e0_run/d2_v3/d2_v3_s{42,43,44}/artifacts/d2_{h,i}_seed{42,43,44}.npz` |
| **D1-B** | **P** = `programme_only` (biology head on 50 Hallmark scores), **F** = `programme_free` (no programme regression; patient-paired cross-modal InfoNCE on the biology view) | `--objective-profile` only | `~/e0_run/d1_v2/artifacts/d1_{p,f}_seed{42,43,44}.npz` |

All six D1-B runs are complete and live: `~/e0_run/d1_audit.log` records `[PASS]
A1_all_six_runs_complete` with `overfit_present: True` and `all_grads_positive: True` for all six.

### 3.4 What "matched" means here, and the three places it fails

- **D2 is matched by construction.** A single `D2_PAIR_MANIFEST.json` enumerates the 40 common
  arguments; both arms record the same `pair_manifest_sha256`
  (`ce1352e0ac7a98334e4fada8178986e8413fac1046ebb67a96f5c3cbc7c2fb0b`) and the same
  `common_config_sha256` (`b7b2441fd9d03a3a00152027efe8c7ada3bedc48e7939f1dfc0b320b02adf1fb`). Both
  arms use `--objective-profile programme_only`; **they differ only in the supervision target table.**
  This is a *within-method* comparison in RankMe's sense.
- **D1-B is matched by construction but differs in the objective.**
  `~/e0_run/d1_v2/D1_PAIR_MANIFEST.json` records `"objective_only_difference": true`. Because the arms
  differ in *objective*, D1 sits closer to a between-method comparison than D2 does and lands further
  outside RankMe's stated scope. §4.7 accounts for this in both directions, since the result there is a
  confirmation of RankMe rather than a refutation.
- **The three-seed replicates of any one arm are unambiguously in scope.** They differ in nothing but
  the seed and are therefore *"different runs of a given method"* on RankMe's own words. §4.2 and §4.3
  are measured there.
- **§4.9's instance 2 (Phase 1b) is not matched.** Its source states it: *"`full` vs `programme_only`
  manifests were not verified as matched on epochs/LR/budget in this run (G0.4). Until they are, the
  rank comparison in §5 is suggestive, not causal"* (`PHASE1B_TARGETED_READOUT.md:147-148`). The three
  diagnostic artifacts record only `configuration_sha256`, `git_commit`, `git_dirty`, `package`,
  `package_root`, `source_tree_sha256` — **no epochs, no learning rate, no token budget, no seed** — and
  `git_dirty` is `True` for all three.
- **§4.9's instance 1 cannot be assessed at all.** See §6.4.

### 3.5 Seed reproducibility on this stack, and the rule it forces

Training is **not seed-reproducible** on this hardware. Re-exporting a surviving checkpoint reproduces
its recorded readout to five significant figures (0.58612 against 0.5861 recorded), so the
export/readout path is deterministic. Retraining the same seed with the same configuration gives a
different model: held-out top-CCA **0.6214 versus 0.5861**, and canonical R1 effective rank **23.39
versus 8.68** (`D2_RESULT.md` §4).

**The rule this forces, and which this paper obeys: quote paired within-run differences for the channel,
never levels.** It does *not* rescue effective rank, because rank is quoted as a **level** in every
practice this paper is about — including RankMe's own, which selects between runs by comparing their
rank values. There is no paired-difference form of "this run has higher rank". §4.1 develops this
asymmetry into the paper's central measurement, and §4.4 establishes that it is a property of training
rather than of the estimator.

### 3.6 Reporting rules adopted for this paper

1. **The result that went against us is §4.7, and it is summarised in the abstract.** Its reading was
   fixed in advance (`NOTEBOOK_ENTRIES/PREDECLARED_D1_necessity_test_20260803T2300Z.md`), before any
   channel number was seen.
2. Every rank number names its **statistic** and its **block**. Numbers that cannot be recomputed under
   the canonical definition are marked `[NOT RECOMPUTABLE]` or `[NOT RECOMPUTED]` and are excluded from
   every count and summary statement rather than carried forward.
3. **Selection-rule counts carry no rhetorical weight** (§4.6). At n = 6 a flawless 6/6 gives an exact
   two-sided binomial p = 0.031; 5/6 gives 0.219; 4/6 gives 0.688. The design has just enough power to
   detect a *perfect* rule and none to detect a good one. No conclusion of the form "metric X beats
   metric Y" is drawn from six pairs, in either direction.
4. No number is compared across the statistics of §3.1, and no number is compared to a published RankMe
   value (§2.1, §2.6).
5. Where an instance's own source file corrects, caveats or withdraws it, that text is quoted.
6. Evidence that cuts against the claim is reported in §4.7, §4.10 and §5.4, at the same prominence as
   the evidence that supports it.

---

## 4. Results

### 4.1 All seven between-arm rank differences ever measured are inside the statistic's own retraining floor — and the floor is bimodal

This is the paper's headline. It is a statement about **usefulness**, and it is measured inside the
same-method regime RankMe reserves for itself.

**The floor, and the measurement it now rests on.** Until 2026-08-04 this section rested on **one**
retraining pair (D2 arm H, seed 42: re-export 8.6809 against retrain 23.3868, a 2.69× move) and said so
in place of a limitation. That n = 1 defect was named in §6.2 as the paper's most valuable missing
measurement, the repeat was predeclared in
`NOTEBOOK_ENTRIES/PREDECLARED_retraining_envelope_20260804T0330Z.md` — including how each of four
outcomes would be read — and it has now reported. **Five identical `programme_only` runs, same seed 42,
same configuration, same data, same schedule. The only source of variation is GPU non-determinism.**

| repeat | rank, raw block | **rank, residualised block** | **channel (top-CCA, 40 untrained targets)** |
|---|---:|---:|---:|
| 1 | 24.481 | 28.320 | 0.6182 |
| **2** | **8.033** | **8.834** | **0.5859** |
| 3 | 24.504 | 28.348 | 0.6123 |
| 4 | 24.990 | 29.106 | 0.6110 |
| 5 | 24.912 | 28.959 | 0.6098 |
| **spread (max / min)** | **3.111×** | **3.295×** | **1.055×** |

*Provenance: `~/e0_run/d1_envelope/rep{1..5}.npz`, readout `~/e0_run/d1_envelope_readout.log`, produced
by `v2/research/rebase/d1_envelope_readout.py` — canonical R1 and top-CCA both imported from
`v2/calibra`, nothing computed inline. Block: exported `wsi_biology`, held-out test partition,
cancer + pooled-TSS residualised, top-CCA at 16 components, the 40 targets neither D1 arm trained on —
identical to `d2_compare`. Reported in
`NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §1.*

**The distribution is bimodal, not a smooth spread, and that is the more useful statement.** Four
repeats agree to within 2% (28.32–29.11) and one lands at a third of them (8.83). This is the same
signature as the liveness gate's 6-of-8 pass rate over eight identical runs (§5.1). A range invites the
reader to imagine a distribution the data does not have; the shape says what actually happens.
**Effective rank on this stack is reproducible about 80% of the time and catastrophically not about 20%
of the time, on identical inputs.**

**Repeat 2 is the cleanest single observation this project has produced.** Same seed, same
configuration: **effective rank falls 3.3×; the channel falls 5%.** One quantity moves by a factor, the
other barely moves at all, and there is no difference between the two runs to attribute either move to.
That single row is this paper's argument without any arm contrast in it, and it is why it appears in the
abstract.

**The number is a FLOOR, twice over**, as predeclared before it existed. First, `programme_only` is this
project's *stable* arm — its channel varies 1.018–1.026× across seeds where `programme_free` varies
1.056×, and its step-200 tripwire rank spans 1.18× across five seeds against `programme_free`'s 6.05×
(§4.3). Measuring
the retraining spread on the stable arm understates it. Second, same-seed repeats exclude seed variation
**entirely**; §4.2 measures that other axis separately and finds it larger still. So 3.295× is a lower
bound on what a practitioner faces, and every quotation of it in this paper says "floor", not
"envelope".

**The comparisons a rank-based selection rule would act on.** Every between-arm rank comparison this
project has ever measured, all canonical R1, each with the block named in §4.9 or below. **Each ratio is
judged against the floor measured on its own block** — 3.295× residualised, 3.111× raw. Comparing a
residualised ratio against the raw floor would put D1-B seed 43 outside it, which is precisely the
raw/residualised confusion §4.5 is about:

| comparison | arms | rank ratio | block | inside that block's retraining floor? |
|---|---|---:|---|---|
| D2 seed 44 | H 9.1426 / I 9.1052 | **1.004×** | resid. | **yes** |
| D2 seed 43 | H 28.7715 / I 34.1168 | **1.186×** | resid. | **yes** |
| Phase 1b, single seed | `full` 38.4834 / `programme_only` 32.0594 | **1.200×** | raw | **yes** (3.111×) |
| D2 seed 42 | H 23.3868 / I 14.8675 | **1.573×** | resid. | **yes** |
| D1-B seed 44 | P 11.1148 / F 6.3937 | **1.738×** | resid. | **yes** |
| D1-B seed 42 | P 29.3813 / F 13.4184 | **2.190×** | resid. | **yes** |
| D1-B seed 43 | P 24.6730 / F 7.6003 | **3.246×** | resid. | **yes — by 1.5%** |

*Provenance: `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`
§5, §6; `v2/research/rebase/nature/D2_RESULT.md` §4; `PHASE1B_TARGETED_READOUT.md` §3, §5, §7. Artifacts
`~/e0_run/d2_v3/d2_v3_s{42,43,44}/artifacts/`, `~/e0_run/d2_v3/recovered_artifacts/`,
`~/e0_run/d1_v2/artifacts/`, and the Phase-1b artifacts at
`/lambda/nfs/geeg/biorag3_persistent_20260711/runs/v21_release_20260720_retry3_resume_safe/artifacts/`.
Every value recomputed under the canonical implementation on 2026-08-04 and reproducing the published
figures exactly; outputs `~/ws_rank/RANK_RECOMPUTE.json`, `~/ws_rank/RANK_RECOMPUTE_P1B.json`. The
raw-block D1-B ratios, quoted in the predeclaration and in §4.7, are 2.02× / 3.09× / 1.68×.*

**The count moves from six of seven to seven of seven — in our favour, and it is reported with the
scepticism a result going the other way would get.** The predeclaration fixed this in advance: an
envelope wider than §4.1's own threshold *"changes our own headline count, in our favour"* and *"is to
be reported with the same scepticism we would apply to a result going the other way, including that it
rests on this one measurement"*. So: the seventh point clears by **1.5%** (3.246× against 3.295×), which
is not a margin anyone should lean on; the floor is one measurement, of five repeats, on one arm, at one
seed, on one stack; and it cannot separate rank-specific variance from stack non-determinism,
architecture or schedule. **The correct quotation is "seven of seven fall inside a floor of 3.295×
measured on five same-seed repeats", never a bare "seven of seven".**

**Against the same measurement, the channel is resolvable.** On those five identical runs the channel
spreads **1.055×** while rank spreads **3.295×**. Against the D2 runs the paired channel difference is
**−0.1325 / −0.1089 / −0.1226** — the same sign 3/3, spread 0.024 on a mean of −0.121, with both
patient- and cancer-clustered bootstrap CIs excluding zero in all three seeds (§4.9, and see §4.6a for
what that contrast does and does not generalise over). The channel readout is **not** exempt from §3.5 —
0.5861 against 0.6214 is a 6% move at one seed — but the channel is quoted as a *paired within-run
difference* and rank is not, and the paired difference is stable where both levels are not. **That
asymmetry is a statement about how the two quantities are used, and it is the practical heart of this
paper.** It now rests on one measurement in which both quantities were computed on the *same five runs*,
rather than on a comparison across tables.

**The consequence for RankMe's own recommended use.** RankMe selects between runs by comparing their
rank values and restricts that comparison to runs of a given method — exactly the setting above. A
criterion whose value moves 3.3× when the same configuration is retrained cannot resolve a
between-configuration difference smaller than that. In particular, the largest D2 difference (1.573×,
seed 42) is comfortably inside the floor, which is an independent reason why that seed's "correct"
ordering carries no weight.

**Stated as a rule a practitioner can apply:** *before using a rank difference to select between
configurations, retrain one configuration with the same seed several times and measure the rank spread;
if the between-configuration difference does not exceed it, the comparison is uninformative.* We have
not seen this check proposed anywhere, and on this stack it would have disqualified **every** rank
comparison this project made. The bimodality is why the rule says *several* repeats and not one: a
single pair draws two runs, and any pair drawn from our four concordant repeats spans at most **1.028×**
— a measured floor of 1.03× would have licensed every comparison in the table above. One repeat is not
a cheaper version of this check; it is a different check that can return the wrong answer.

**What this rests on, stated plainly.** Five repeats, one arm, one configuration, one seed, one stack,
no interval. That is much better than the n = 1 estimate it replaces and it is still thin, and §4.2
exists because it is thin: the variance decomposition estimates the nuisance term from **8 within-arm
degrees of freedom** across four arms rather than from repeats of one, and reaches the same conclusion
without using this floor at all. The two should be read together, and if only one survives review it
should be §4.2.

### 4.1a The paper audited against its own criterion — every rank comparison it makes or relies on

§4.1's rule is the paper's central criterion, and the paper applies it to RankMe's recommended use.
It has also applied it to **itself**, five times, each instance discovered separately and late: the
momentum fix (3.29×, §5.4), the choice of `m = 0.999` over `m = 0.99` (1.26×, §5.4), §4.7.4's
surviving necessity violation (3.73×), the decorrelation ablation (1.85×, §4.9a) and §5.1's
instance 2 (≈3.2×, flagged in a notebook entry and left unrewritten until now). Finding them one at
a time is how a referee finds the sixth. **This subsection enumerates all of them at once**, and it
is a cheap thing for us to write and an expensive thing for a referee to construct independently.

**And this subsection was itself published with a scope error, which is repaired here — badly enough
that it is the first thing the subsection should say.** The first version of the audit judged fifty
comparisons against **one** floor, canonical R1 on the exported `wsi_biology` block, because that was
the only floor this project had ever measured. Thirteen rows sat on a statistic or a block with no
floor at all and were nonetheless printed with a verdict of "**no**"; a further eleven sat on the
fixed held-out probe and were judged against a floor measured on the exported artifact. **A
comparison with no measured floor has not failed the criterion. The criterion has not been applied to
it**, and printing those two outcomes in the same column flattered us: it made an unmeasured ruler
look like a failed test. The repair has two halves. First, **every floor recoverable from the five
repeats has now been measured** — a re-derivation from exports that already existed, not a new
training run. Second, a row whose statistic-and-block has no floor now records `clears: null`, prints
**unjudgeable**, and is counted separately; the coupling between "no floor" and "no verdict" is
enforced by `p2_floor_audit.check()` rather than left to convention.

#### The floors, measured — every statistic T1 scores, every view §4.5(c) reads

§4.1's five identical `programme_only` retrains (seed 42, GPU non-determinism the only source of
variation) export all three co-trained views, so a floor for every other statistic and every other
view was already sitting on disk. Measured by `v2/research/rebase/p2/p2_envelope_floors.py`,
CPU only, thread-capped, on a box workspace built with `git -c core.autocrlf=false archive` and
verified **543/543 files by git blob SHA-1** before it was run:

| statistic | block | min | max | **floor** | other four agree to | divergent run | bimodal? |
|---|---|---:|---:|---:|---:|---|:---:|
| **the exported `wsi_biology` block — the one §4.1 measures, every statistic T1 scores** | | | | | | | |
| R1 | residualised | 8.8340 | 29.1057 | **3.295×** | 1.028× | rep2 (low) | **yes** |
| RankMe (residualised) | residualised | 8.8359 | 29.1105 | **3.295×** | 1.028× | rep2 (low) | **yes** |
| R1 | raw | 8.0326 | 24.9895 | **3.111×** | 1.021× | rep2 (low) | **yes** |
| α-ReQ \|α−1\| | residualised | 1.1419 | 2.6247 | **2.299×** | 1.121× | rep2 (high) | no |
| R3 | residualised | 6.3663 | 14.5795 | **2.290×** | 1.028× | rep2 (low) | **yes** |
| R2 | residualised | 5.5972 | 12.4506 | **2.224×** | 1.025× | rep2 (low) | **yes** |
| R2 | raw | 5.1004 | 10.9157 | **2.140×** | 1.014× | rep2 (low) | **yes** |
| R3 | raw | 5.1004 | 10.9157 | **2.140×** | 1.014× | rep2 (low) | **yes** |
| α-ReQ \|α−1\| | raw | 1.2291 | 2.6304 | **2.140×** | 1.137× | rep2 (high) | no |
| RankMe (raw, uncentred, ε = 1e-7) | raw | 1.9883 | 3.6011 | **1.811×** | 1.020× | rep2 (low) | no |
| PR_rownorm | residualised | 3.7152 | 5.4474 | **1.466×** | 1.024× | rep2 (low) | no |
| PR | raw | 2.5913 | 3.7480 | **1.446×** | 1.047× | rep2 (low) | no |
| PR_rownorm | raw | 2.5913 | 3.7480 | **1.446×** | 1.047× | rep2 (low) | no |
| PR | residualised | 3.1643 | 4.4913 | **1.419×** | 1.051× | rep2 (low) | no |
| stable rank | raw | 1.7729 | 2.2082 | **1.246×** | 1.035× | rep2 (low) | no |
| stable rank | residualised | 2.1977 | 2.6898 | **1.224×** | 1.054× | rep2 (low) | no |
| top-CCA (16 components, 40 untrained targets) | residualised | 0.5859 | 0.6182 | **1.055×** | 1.014× | rep2 (low) | no |
| hard numerical matrix rank | raw | 256.0000 | 256.0000 | **1.000×** | 1.000× | rep1 (high) | no |
| hard numerical matrix rank | residualised | 256.0000 | 256.0000 | **1.000×** | 1.000× | rep1 (high) | no |
| **the `wsi_biology` + `rna_biology` positive pair — LiDAR's block** | | | | | | | |
| LiDAR | residualised | 38.7877 | 41.1002 | **1.060×** | 1.024× | rep2 (low) | no |
| LiDAR | raw | 39.6212 | 40.9816 | **1.034×** | 1.026× | rep2 (low) | no |
| **the `rna_biology` view, canonical R1 — same five runs, same statistic, other view** | | | | | | | |
| R1 | raw — `rna_biology` view | 23.9832 | 24.5300 | **1.023×** | 1.017× | rep2 (low) | no |
| R1 | residualised — `rna_biology` view | 27.2245 | 27.7497 | **1.019×** | 1.016× | rep2 (low) | no |
| **the `full_biology` view, canonical R1 — same five runs, same statistic, other view** | | | | | | | |
| R1 | residualised — `full_biology` view | 29.8042 | 30.3991 | **1.020×** | 1.015× | rep2 (high) | no |
| R1 | raw — `full_biology` view | 26.3285 | 26.6921 | **1.014×** | 1.009× | rep2 (high) | no |

*Provenance: `~/e0_run/d1_envelope/rep{1..5}.npz` (each file's SHA-256 recorded in the output),
targets `~/e0_run/data/frozen_rna_targets.npz`, cancer + pooled-TSS cross-fitted residualisation at seed 42,
α-ReQ at the authors' estimator over eigenvalue ranks 11–50 and LiDAR at q = 2, δ = 1e-4 — the
settings T1 is scored with. R1/R2/R3 and the channel are imported from `v2/calibra/spectral.py`, the
four published alternatives from `p2_competing_metrics.py`, the hard rank is
`numpy.linalg.matrix_rank`; **nothing is computed inline**. Output vendored at
`v2/research/rebase/p2/figures/data/ws_floor/out/P2_ENVELOPE_FLOORS.json`, and every value in the
table above is re-read out of it by `v2/tests/test_p2_floor_audit.py`. The two floors §4.1 already
published are now carried with a **second, independent source**: 3.295× and 3.111× are parsed from
`d1_envelope_readout.py`'s printed log, and this recomputation from the artifacts agrees with them to
four decimal places on both extremes.*

**"Other four agree to" and "divergent run" are the SHAPE, and they are reported per statistic
because whether the shape survives a change of statistic is on this paper's thesis.** The rule is
fixed in the module, not chosen per row: the divergent run is the one whose removal minimises the
remaining four's fold; a floor is called bimodal when those four agree to within 5% *and* the full
fold is at least twice theirs. Under that rule §4.1's own floor reads exactly as §4.1 describes it —
repeat 2, four others within 2.8%, separation 3.21×.

**Five things this table says. The first two cost us.**

1. **The floor is a property of the statistic, and it varies by as much as the effect the paper is
   about.** On one block — exported `wsi_biology`, residualised, the same five runs — it runs from
   **1.000×** (hard numerical rank, which does not move at all) through 1.224× (stable rank), 1.419×
   (participation ratio) and 2.224×/2.290× (R2/R3) to **3.295×** (R1). §4.3's heading says the floor
   is a property of the arm *and not of the statistic*; on this measurement it is a property of both,
   and §4.3 is corrected. Concretely: judging an R3 comparison against R1's floor, which the first
   version of this table did on fourteen rows, is **1.4× too strict**.
2. **RankMe as published is markedly more reproducible on our own artifacts than our centred
   statistic is.** On the raw block its floor is **1.811×** against canonical R1's **3.111×**, on the
   same five runs. The mechanism is the one §4.6 already names for RankMe's D2 advantage: the
   uncentred normalisation retains the mean-offset direction, and every row of every exported view
   has L2 norm exactly 1.000, so that direction is both large and stable. On the *residualised* block
   the mean is gone and the two statistics coincide — floors of 3.295× and 3.295×. **This is a result
   against our own instrument**, and it is why the RankMe row is the one selection in the whole audit
   that clears a floor its own statistic and block license.
3. **The floor is also a property of the VIEW, and the divergent run is divergent in only one of
   them.** The same five runs spread **3.295×** on `wsi_biology`, **1.019×** on `rna_biology` and
   **1.020×** on `full_biology`. Repeat 2 — whose WSI-view rank is a third of its siblings' — sits
   within 1.6% of them on the RNA view and is the *high* member on the full view. The catastrophic
   one-in-five is not a property of the run; it is a property of that run's **WSI encoder**. Every
   one of §4.5(c)'s twelve `rna_biology`/`full_biology` arm comparisons is resolvable against its own
   view's floor, and none of the six on `wsi_biology` is — same artifacts, same statistic, same runs.
4. **The bimodal shape is statistic-dependent; the divergence is not hidden from any statistic.**
   Every statistic that moves at all identifies repeat 2 as the outlier, and in the degradation
   direction (lower rank, higher α-ReQ α). What changes is the magnitude, and therefore the shape:
   the four-concordant-plus-a-factor signature §4.1 calls bimodal survives under R1, R2, R3 and
   residualised RankMe and under nothing else. The same run is 1.22× away under stable rank, 1.06×
   under LiDAR and 1.00× under the hard rank. **The divergence is a redistribution of spectral mass
   in the tail — which the entropy-based ranks are sensitive to by construction and the
   dominant-subspace statistics are not.** That is a narrower and better-supported reading of §4.1's
   asymmetry than "rank is unreliable", and it is the reading the paper should carry.
5. **A floor of 1.000× licenses everything, and one statistic has one.** The hard numerical rank is
   pinned at 256 in every repeat, every view and both blocks. §4.1 warns that a floor drawn from a
   single pair of concordant repeats would have been 1.028× and would have licensed every comparison
   in this paper; a statistic whose measured floor is exactly 1 is the same failure with five repeats
   instead of two. Nothing in the audit is judged against it.

**Four blocks still have no floor, and the absence is now a stated result rather than a silence.**
The five repeats were **exported, not probed**: `~/e0_run/d1_envelope/` holds `rep{n}.npz` and a
per-run `~/e0_run/d1_envelope/rep{n}/train_metrics.jsonl`, and neither carries a probe forward pass, a training batch's
activations, or a gate batch. So no floor exists, and none can be recovered without a GPU, for **the
fixed held-out probe**, **the in-run training batch**, **the 16-patient gate batch**, or **the
282-patient live checkpoint**. What each would cost is recorded per block in
`P2_ENVELOPE_FLOORS.json`'s `absent_blocks`, which a test asserts the audit's own list against.
**Every rank number in §5 is on the first of those four blocks**, as are §4.9a's decorrelation
ablation and §4.4(3)'s probe repeat. What §5 would need in order to become judgeable at all is one
specific run: five same-seed repeats of the `programme_free` / 500-step configuration with
`d1_momentum_probe.py` attached, read at a fixed step. Until that exists, §5's rank comparisons are
neither confirmed nor refuted by this paper's criterion — they are outside its reach, and §5.4 and
§6.2 should say so in those words.

**The caveat that travels with every row of that table.** Each floor is **n = 5**, one arm, one seed,
one stack, no interval — and it is a **floor twice over**, exactly as §4.1 says of the two it
extends: `programme_only` is this project's *stable* arm, and same-seed repeats exclude seed
variation entirely, which §4.2 measures as the larger term. A table of twenty-six numbers is not a
table of estimated distributions and must not be read as one; the right sentence is "a floor of
2.290× measured on five same-seed repeats", never "R3 varies 2.29×".

#### The audit, re-run against the floors that now exist

Every rank comparison this paper makes or relies on is below — Results, the worked example, the
appendices, the figure plan, `QUEUE_ANCHORING.md` and `LIVENESS_GATE_DESIGN.md`. For each: the two
values, the fold, the statistic, the **block**, which floor applies, whether it clears, and — where
it does not — what the claim rests on instead. Where a comparison legitimately falls outside the
criterion the exemption is **stated**, never taken silently.

**Read the statistic and block columns first.** A ratio may be judged only against a floor measured
on its own statistic *and* its own block, and both halves of that rule have already produced a wrong
verdict on this project. On the block: D1-B seed 43 is **3.246×** residualised and **3.091×** raw, and
its residualised figure against the *raw* floor of 3.111× reads as **outside** when on its own block
it is inside (rows 6 and 9). On the statistic: fourteen rows were judged against R1's 3.295× while
measuring R3, whose own floor is 2.290×. A ‡ marks a row that still carries a mismatch note — of
statistic, of block, or of kind — and the note travels with the row in
`v2/research/rebase/p2/floor_audit.json`.

56 rank comparisons are enumerated; 25 of them are selections between candidate configurations. **13 of those 25 do not clear the floor their own statistic and block license, 1 clears it, and 11 cannot be judged at all** because no floor has ever been measured for the block they sit on. 5 rows are exempt with the reason stated, 25 of the 56 are unjudgeable for want of a floor, and 35 carry an explicit statistic-, block- or kind-mismatch note.

| # | § | comparison | values | ratio | statistic | block | floor | clears? | what the claim rests on |
|---|---|---|---|---:|---|---|---|:---:|---|
| 1 | §4.1 | D2 seed 44 — arm H against arm I | 9.143 / 9.105 | **1.004×** | R1 | residualised | 3.295× | **no** | nothing — the section states it: the difference is not resolvable. |
| 2 | §4.1 | D2 seed 43 — arm H against arm I | 28.771 / 34.117 | **1.186×** | R1 | residualised | 3.295× | **no** | nothing — the section states it: the difference is not resolvable. |
| 3 | §4.1 | D2 seed 42 — arm H against arm I | 23.387 / 14.868 | **1.573×** | R1 | residualised | 3.295× | **no** | nothing — the section states it: the difference is not resolvable. |
| 4 | §4.1, §4.7.2 | D1-B seed 44 — `programme_only` against `programme_free` | 11.115 / 6.394 | **1.738×** | R1 | residualised | 3.295× | **no** | nothing — the section states it: the difference is not resolvable. |
| 5 | §4.1, §4.7.2 | D1-B seed 42 — `programme_only` against `programme_free` | 29.381 / 13.418 | **2.190×** | R1 | residualised | 3.295× | **no** | nothing — the section states it: the difference is not resolvable. |
| 6 | §4.1, §4.7.2 | D1-B seed 43 — `programme_only` against `programme_free` | 24.673 / 7.600 | **3.246×** | R1 | residualised | 3.295× | **no** | nothing — and it is the closest of the seven, clearing **by 1.5%** in the wrong direction (3.246× against 3.295×). |
| 7 | §4.1, §4.9 | Phase 1b, single seed — `full` against `programme_only` | 38.483 / 32.059 | **1.200×** | R1 | raw | 3.111× | **no** | nothing. Its arms were never verified matched (§3.4), its channel moved 0.002, and §4.9 already reads it as *neither quantity resolvable at one seed*. |
| 8 | §4.1 (provenance note), §4.7 | D1-B seed 42 — the same pair on the **raw** block | 25.307 / 12.532 | **2.019×** | R1 | raw | 3.111× | **no** | nothing. |
| 9 | §4.1 (provenance note), §4.7 | D1-B seed 43 — the same pair on the **raw** block | 21.480 / 6.950 | **3.091×** | R1 | raw | 3.111× | **no** | nothing — **and this row is why block-matching is load-bearing**: the same pair reads 3.246× residualised, and 3.246× judged against the *raw* floor of 3.111× would read as OUTSIDE when on its own block it is inside. |
| 10 | §4.1 (provenance note), §4.7 | D1-B seed 44 — the same pair on the **raw** block | 10.089 / 6.008 | **1.679×** | R1 | raw | 3.111× | **no** | nothing. |
| 11 | §4.1, §3.5 | the superseded n = 1 envelope — D2 H42 re-export against retrain | 23.387 / 8.681 | **2.694×** | R1 | residualised | 3.295× ‡ | **no** | nothing — superseded by the n = 5 measurement. It is listed because §4.1's own argument is that a single pair can return a floor that licenses everything, and 2.694× is the pair that did. |
| 12 | §4.2 | within-arm seed fold, D2 H (Hallmark) — three seeds, nothing else varied | 28.771 / 9.143 | **3.147×** | R1 | residualised | 3.295× ‡ | **no** | §4.2's cleanest single object: one arm, no arm contrast, 3.15× in rank against 0.0012 in channel; the channel over the same three seeds moves 1.8–5.6%. |
| 13 | §4.2 | within-arm seed fold, D2 I (perturbation basis) — three seeds, nothing else varied | 34.117 / 9.105 | **3.747×** | R1 | residualised | 3.295× ‡ | yes | the seed alone moves rank **further than five same-seed retrains do** (3.75× against 3.295×); the channel over the same three seeds moves 1.8–5.6%. |
| 14 | §4.2 | within-arm seed fold, D1 P (`programme_only`) — three seeds, nothing else varied | 29.381 / 11.115 | **2.643×** | R1 | residualised | 3.295× ‡ | **no** | the stable arm; the channel over the same three seeds moves 1.8–5.6%. |
| 15 | §4.2 | within-arm seed fold, D1 F (`programme_free`) — three seeds, nothing else varied | 13.418 / 6.394 | **2.099×** | R1 | residualised | 3.295× ‡ | **no** | the unstable arm; the channel over the same three seeds moves 1.8–5.6%. |
| 16 | §4.3 | `programme_free` step-200 spread across five seeds | 45.646 / 7.545 | **6.050×** | R3 | training batch, in-run | **none measured** | **unjudgeable** | five seeds, one step, no interval; the claim is that the spread is large. |
| 17 | §4.3 | `programme_only` step-200 spread across five seeds | 112.078 / 94.952 | **1.180×** | R3 | training batch, in-run | **none measured** | **unjudgeable** | five seeds, one step, no interval; the claim is that the spread is small. |
| 18 | §4.3 | **the two spreads against each other** — 6.05× against 1.18× | 6.050 / 1.180 | **5.125×** | R3 | training batch, in-run | **none measured** | **unjudgeable** | five seeds per arm and no interval on either spread. The claim survives because the two arms' step-200 *levels* are an order of magnitude apart and the supervised arm's five seeds agree to 18%, not because the ratio of the two spreads clears anything. |
| 19 | §4.4(3), §5.4 | controlled 200-step probe, **seed held fixed** — min(m = 0.999) against max(m = 0) | 6.920 / 1.980 | **3.495×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | it clears by **6%** — no more comfortable than the 1.5% §4.1 declines to lean on — and it holds the **seed** fixed, which §4.2 measures as the dominant term. The seed-varied version of this same check is row 4.9a/5.4 below and it does not clear. |
| 20 | §4.5(a) | the six arm pairs under **R1** — canonical Roy & Vetterli order 1 | 1.004× – 3.246× over 6 | **3.246×** | R1 | residualised | 3.295× | **no** | the **instability of the ordering** across statistics, not any magnitude. Against R1's own floor of 3.295×, **0 of the 6** pairs are resolvable. **The count of resolvable pairs is itself statistic-dependent — 0 under R1, 1 under R2 and R3, 2 under the two participation ratios — which is §4.5's own thesis applied to the criterion rather than to the verdict.** |
| 21 | §4.5(a) | the six arm pairs under **R2** — the order-2 Hill number of the singular values (`d1_audit.py`) | 1.014× – 2.381× over 6 | **2.381×** | R2 | residualised | 2.224× | yes | the **instability of the ordering** across statistics, not any magnitude. Against R2's own floor of 2.224×, **1 of the 6** pairs are resolvable (P43/F43). **The count of resolvable pairs is itself statistic-dependent — 0 under R1, 1 under R2 and R3, 2 under the two participation ratios — which is §4.5's own thesis applied to the criterion rather than to the verdict.** |
| 22 | §4.5(a) | the six arm pairs under **R3** — R2 on L2-normalised rows (the geometry probes and the in-run tripwire) | 1.042× – 2.476× over 6 | **2.476×** | R3 | residualised | 2.29× | yes | the **instability of the ordering** across statistics, not any magnitude. Against R3's own floor of 2.290×, **1 of the 6** pairs are resolvable (P43/F43). **The count of resolvable pairs is itself statistic-dependent — 0 under R1, 1 under R2 and R3, 2 under the two participation ratios — which is §4.5's own thesis applied to the criterion rather than to the verdict.** |
| 23 | §4.5(a) | the six arm pairs under **PR** — the eigenvalue participation ratio | 1.012× – 1.506× over 6 | **1.506×** | PR | residualised | 1.419× | yes | the **instability of the ordering** across statistics, not any magnitude. Against PR's own floor of 1.419×, **2 of the 6** pairs are resolvable (P43/F43, P44/F44). **The count of resolvable pairs is itself statistic-dependent — 0 under R1, 1 under R2 and R3, 2 under the two participation ratios — which is §4.5's own thesis applied to the criterion rather than to the verdict.** |
| 24 | §4.5(a) | the six arm pairs under **PR_rownorm** — the eigenvalue participation ratio on normalised rows | 1.020× – 1.652× over 6 | **1.652×** | PR_rownorm | residualised | 1.466× | yes | the **instability of the ordering** across statistics, not any magnitude. Against PR_rownorm's own floor of 1.466×, **2 of the 6** pairs are resolvable (P44/F44, P43/F43). **The count of resolvable pairs is itself statistic-dependent — 0 under R1, 1 under R2 and R3, 2 under the two participation ratios — which is §4.5's own thesis applied to the criterion rather than to the verdict.** |
| 25 | §4.5(b) | the four **R3** rows of §4.5(b)'s table — D2 seeds 43 and 44, raw and residualised | 1.012× – 1.079× over 4 | **1.079×** | R3 | raw and residualised | 2.14× ‡ | **no** | the **sign flipping with the block**, on differences of 1–8%. The claim is that the block choice is worth more than the difference it adjudicates — which is a statement that the difference is unresolvable, made in the other direction. |
| 26 | §4.5(c) | the six arm pairs on the co-trained **`wsi_biology`** view | 1.004× – 3.246× over 6 | **3.246×** | R1 | residualised | 3.295× | **no** | the information ordering being identical under all three views for all six pairs while the rank ordering is not — a sign claim. Against this view's own floor of **3.295×**, **0 of the 6** pairs are resolvable. This is the view §4.1's floor is measured on and the only one of the three on which no pair is resolvable. `rna_biology` is also partly circular, so the 2/9 count is already not quoted as a rate. |
| 27 | §4.5(c) | the six arm pairs on the co-trained **`rna_biology`** view | 1.116× – 3.014× over 6 | **3.014×** | R1 | residualised, `rna_biology` view | 1.019× | yes | the information ordering being identical under all three views for all six pairs while the rank ordering is not — a sign claim. Against this view's own floor of **1.019×**, **6 of the 6** pairs are resolvable. **The floor is a property of the VIEW.** The same five retrains that spread 3.295× on `wsi_biology` spread 1.019× on `rna_biology` and 1.020× on `full_biology`: the one divergent repeat lost its WSI-view rank and kept its RNA-view rank to within 2%. Every between-arm difference on this view is therefore resolvable, and none on `wsi_biology` is — on the same artifacts, the same statistic and the same runs. |
| 28 | §4.5(c) | the six arm pairs on the co-trained **`full_biology`** view | 1.042× – 5.250× over 6 | **5.250×** | R1 | residualised, `full_biology` view | 1.02× | yes | the information ordering being identical under all three views for all six pairs while the rank ordering is not — a sign claim. Against this view's own floor of **1.020×**, **6 of the 6** pairs are resolvable. **The floor is a property of the VIEW.** The same five retrains that spread 3.295× on `wsi_biology` spread 1.019× on `rna_biology` and 1.020× on `full_biology`: the one divergent repeat lost its WSI-view rank and kept its RNA-view rank to within 2%. Every between-arm difference on this view is therefore resolvable, and none on `wsi_biology` is — on the same artifacts, the same statistic and the same runs. |
| 29 | §4.6, §4.6a, T1 | **the six pairs every one of T1's twelve metric rows is scored on** | 1.004× – 3.246× over 6 | **3.246×** | R1 | residualised | 3.295× | **no** | nothing. §4.6 gives three reasons not to read its counts — n = 6, D2 s44's 1.4 sampling SDs, and §4.6a's coordinate choice. **This is a fourth, and the paper does not currently state it**: not one of the six pairs is resolvable under the one metric that has a measured floor. |
| 30 | §4.6, T1 | RankMe as published, the three D2 pairs — the row §4.6 quotes against ours | 1.248× – 3.382× over 3 | **3.382×** | RankMe (raw, uncentred, ε = 1e-7) | raw | 1.811× | yes | **the one selection in this paper that clears a floor its own statistic and block licenses — and it is not ours.** RankMe as published has a retraining floor of **1.811×** on the raw exported block, against canonical R1's 3.111× on the same five runs, because its uncentred normalisation is dominated by the mean-offset direction, which is the stable part. Only **1 of the 3** D2 pairs clears it (s43, 3.382×; s42 is 1.677× and s44 is 1.248×), so the 3/3 count still rests on two unresolvable orderings — but the criterion no longer refuses the row, and the paper must say that the published metric is the more reproducible one on our own artifacts. |
| 31 | §4.7.1 | **the predeclared violation threshold** (2.0× fold) against the floor | 3.295 / 2.000 | **1.647×** | R1 | residualised | 3.295× ‡ | **no** | the fact itself: the fold a pair had to exceed to count as a necessity violation was **1.65× smaller than the floor** and below §4.2's 2.10–3.75× seed band, so the 66-pair scan could only ever return violations that are unresolvable. Both violations it returned (rows below) are. |
| 32 | §4.7.4 | the surviving necessity violation — H44 against I43 | 34.117 / 9.143 | **3.732×** | R1 | residualised | 2.1–3.75× ‡ | **no** | +0.1101 more channel on the lower-rank artifact, a gap comparable to the headline D2 arm effect. It is presented as supporting, not load-bearing, and is partially pre-empted by Aldeneh et al. (ICASSP 2025). It clears the band by **0.5%**. |
| 33 | §4.7.4, F6(d) | the second violation of the 66-pair scan — P44 against I43 | 34.117 / 11.115 | **3.070×** | R1 | residualised | 3.295× | **no** | +0.1206 more channel. **The draft names this pair only in `P2_FIGURES.md` F6(d) and attaches no floor verdict to it**; it is inside the 3.295× floor and inside §4.2's band, and it is cross-experiment as well as cross-seed. |
| 34 | §4.8 | dilution d = 0 → 0.80, raw block | 196.187 / 161.226 | **1.217×** | R1 | raw | 3.111× ‡ | **exempt** | monotonicity over seven nested levels and the *ratio* of the two fractional changes (−3.10% rank against −66.7% channel), both read from the same representation through the same instrument. Single seed, single donor draw, no interval on any level-to-level difference. |
| 35 | §4.8 | dilution d = 0 → 0.80, residualised block | 210.179 / 203.667 | **1.032×** | R1 | residualised | 3.295× ‡ | **exempt** | monotonicity over seven nested levels and the *ratio* of the two fractional changes (−3.10% rank against −66.7% channel), both read from the same representation through the same instrument. Single seed, single donor draw, no interval on any level-to-level difference. |
| 36 | §4.9, §6.4 | the historical decorrelation instance — “+107% rank at flat benchmark” | 103.3 / 49.9 | **2.070×** | unknown — predates the consolidation | unknown | **none measured** | **exempt** | nothing. It is excluded from every count on **provenance**, not on the floor — and it happens also to be inside the floor, which is a second reason and not the operative one. |
| 37 | §4.9 | the “16/16” instance — pinned at the batch size in every arm | 16 / 16 | **1.000×** | hard numerical matrix rank | 16-patient train batch | **none measured** | **exempt** | the co-measured collapse evidence — patient cosine 0.9999, retrieval 0.000 *below* its 0.062 chance level — and the fact that the R3 rank of the same objective falls 12.88 → 1.00. The rank column carries nothing. |
| 38 | §4.9, §5.1, §5.4 | D1-A epoch 39 — `programme_only` 9.81 against `programme_free` 1.71 | 9.810 / 1.710 | **5.737×** | R3 | 282 held-out patients, live checkpoint | **none measured** | **unjudgeable** | a collapse verdict with co-measured evidence: RNA-view mutual cosine 0.986, hard rank 11, and no exported readout at all. Its own source entry forbids concluding anything about supervision from it. |
| 39 | §4.10 | clean in-batch InfoNCE, step 0 → 50 | 12.880 / 1.000 | **12.880×** | R3 | diagnostic script, train batch | **none measured** | **unjudgeable** | positive and worst-negative cosine both 0.9993 and a minimum margin of −0.0001, co-measured. §4.10 is where the paper says rank works. |
| 40 | §4.10 | the boundary — D2 s44's rank against its own nominal dimensionality (256) | 256.000 / 9.143 | **28.001×** | R1 | residualised | **none measured** | **unjudgeable** | held-out channels of 0.5983 and 0.4757 against a **measured** permutation null of 0.140 — a co-measured quantity, not a ratio. A representation at 3.6% of nominal rank was carrying a large, permutation-significant channel. |
| 41 | §5.1, instance 2 | **the residual**: gate rank 5.81 (16 memorised patients, frozen 64-key queue) against ~1.8 at cohort scale (3,118 streaming patients, live 4,096-key queue) | 5.810 / 1.800 | **3.228×** | R3 | TWO different blocks — see the note | **none measured** | **exempt** | two things, neither of them a rank difference. (i) The gate's own **binary** pass criterion — contrastive 0.012–0.057 against a ≤ 0.10 bar with retrieval 16/16 on three seeds. (ii) The cohort-scale arm being independently **collapsed**: RNA-view mutual cosine 0.977 / 0.986, hard rank 9 / 11, against the supervised sibling at 7.38 / 7.35. The claim — *a gate that certifies memorisation of 16 does not certify learning at 3,118* — is a statement about what the gate's regime removes, not a selection between two configurations. |
| 42 | §5.1, instance 3; S1 | the five-arm regulariser sweep — widest between-arm fold at a shared step (50) | 4.080 / 2.620 | **1.557×** | R3 | fixed 256-patient held-out probe | **none measured** | **unjudgeable** | an absolute-level claim, not a between-arm one: **all five** arms fall from a verified common initialisation of 67.55 to 1.59–3.43, i.e. ≥ 19.7×, including both regularisers at zero. The between-arm folds carry nothing and are not quoted. |
| 43 | §5.4 row 1, §5.2 | m = 0.999 against m = 0, **one seed**, step 600 | 7.420 / 2.810 | **2.641×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | nothing — §5.4 says so. |
| 44 | §5.4 row 2 | m = 0.999 against m = 0, **worst case over three seeds**, 500 steps | 10.450 / 3.180 | **3.286×** | R1 | fixed held-out probe | **none measured** | **unjudgeable** | a **binary training outcome**: `programme_free` completed 40 epochs uncollapsed 0 of 3 seeds before the fix and 3 of 3 after, with a channel and paired bootstrap intervals where no export existed at all. It fails the floor by **0.3%**. |
| 45 | §5.4 row 3 | the same replication under the tripwire statistic | 6.850 / 2.810 | **2.438×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | as row above. |
| 46 | §5.4 limit 2, Appendix C | **the value the project actually runs** — m = 0.999 over m = 0.99, step 600 | 7.420 / 5.880 | **1.262×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | nothing. The binary outcome supports **momentum against none**, not this value over its neighbour; `m = 0.999` is retained because it measured best in a sweep whose resolution the paper does not trust. |
| 47 | §5.2 prose | m = 0.999 against m = 0 at step 400 — the widest fold past step 150 | 7.840 / 2.180 | **3.596×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | nothing quotable. **§5.2's prose says the effect is “2.6–3.3× at every step past 150”, which disagrees with its own table**: the per-step folds are 3.363× (200), 2.208× (300), 3.596× (400), 3.132× (500), 2.641× (600), and 4.343× at step 100. Flagged, not substituted (§5 is another agent's). This single-seed fold does exceed the floor at one step, but it varies no seed, and the seed-varied worst case above does not. |
| 48 | §5.2 measurement 2 | m = 0.999 (6.89) against m = 0.99 (6.65) at step 100 | 6.890 / 6.650 | **1.036×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | **not the ordering.** What falsifies MoCo's staleness account here is that key-to-encoder agreement varies 2.06× (0.908 against 0.441) while rank is *indistinguishable* — an equality claim, which a 1.04× difference supports. The sentence “the best-agreeing arm does not have the best rank” reads an ordering off that 1.04× and should be restated as the equality it is. |
| 49 | §5.2 measurement 2 | m = 0.999 (6.89) against no momentum (2.58) at step 100 | 6.890 / 2.580 | **2.671×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | the same binary outcome as §5.4; one seed. |
| 50 | §5.2 measurement 3 | capacity 64 (6.17) against capacity 4,096 (2.16), fixed key encoder | 6.170 / 2.160 | **2.857×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | one seed, and the section's own admission that capacity confounds anchoring quality with negative count. It is one of three measurements said to rule the staleness account out, and it is inside the floor. |
| 51 | §5.2 turnover falsification | the discriminating τ/T prediction — P4 (3.67) against P2 (3.53) | 3.670 / 3.530 | **1.040×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | **the predeclared reading, not the ordering.** The criterion required P4 to be the worst of its group and fixed “fails ≤ 3.5” in advance; P4 read 3.67 and so did not fail. That is a predicted effect being ABSENT. The draft's “the discriminating one inverted” reads an ordering off a 1.04× difference and overstates what the data support. |
| 52 | §5.2 turnover falsification | “nearly flat in capacity” — 2,048 (3.53) against 8,192 (3.67) at m = 0.95 | 3.670 / 3.530 | **1.040×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | an **equality** claim, and the same 1.04× that cannot carry the ordering above does carry this. The two rows are the same two numbers used two ways, and are listed together deliberately. |
| 53 | §5.2 turnover falsification | m = 0.95 → 0.999 at fixed capacity 4,096 — 2.91 against 7.82 | 7.820 / 2.910 | **2.687×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | one seed per cell; the surviving statement is monotonicity in m across five values. |
| 54 | §4.9a | decorrelation 0.0 → 0.04 at m = 0.999, step 400 (R3 — the column quoted) | 8.010 / 4.320 | **1.854×** | R3 | fixed held-out probe | **none measured** | **unjudgeable** | **monotonicity across three levels and a co-measured contradicting quantity**: the RNA-view mutual cosine rises 0.4774 → 0.7657 → 0.8696 on the *same runs* as rank rises. One seed per level. The magnitude carries nothing; the direction and the cosine do. Figure F9. |
| 55 | §4.9a | the same three runs under the canonical statistic | 12.200 / 6.290 | **1.940×** | R1 | fixed held-out probe | **none measured** | **unjudgeable** | as above. Recorded because the 1.85× the notebook entry quotes is the **R3** column; statistic-matched to the floor it is 1.940×, and both are inside. |
| 56 | §4.1a, §5.4, §6.2 | **the like-for-like pair §5.4 says does not exist** — `ablate_decorr0.04` against `mseed_m0.999_s42`, identical configuration and seed, step 400 | 12.200 / 11.440 | **1.066×** | R1 | fixed held-out probe | **none measured** | **unjudgeable** | nothing may be read off it as a floor. It is recorded because §5.4 and §6.2 both say **no like-for-like measurement exists in this regime**, and one does: the two runs share momentum, decorrelation, capacity, learning rate, seed and step-0 state (67.55 / 101.38 / 0.0342 / 0.3650), and `d1_momentum_probe.py` has no schedule that depends on the step budget. Over the eight shared logged steps the widest canonical fold is 1.128× (step 150). |

‡ = the row carries a `floor_note`: its statistic, its block or its kind does not match the floor it
is judged against, and the note says which. *Provenance and machine-checkable form:*
`v2/research/rebase/p2/floor_audit.json`, checked by `v2/research/rebase/p2/p2_floor_audit.py`
and `v2/tests/test_p2_floor_audit.py`, which **re-reads every value out of the file it came from**
and fails if a ratio in this table disagrees with its source, if a floor disagrees with the
measurement file it was read from, or if a row records a verdict without a floor to record it
against. Both tables above are *generated* from that list; neither is maintained by hand.

**What the audit found that the paper did not already say.** In descending order of how much it costs
us.

1. **Nothing in this paper clears a rank floor except the metric we are arguing against.** Of the 25
   selections between candidate configurations, **13 fail a floor their own statistic and block
   license, 11 cannot be judged at all, and exactly 1 clears** — row 30, **RankMe as published**,
   whose floor (1.811× raw) is nearly half of ours (3.111× raw) on the same five runs. The two rows
   the previous version of this table reported as clearing (§4.4(3)'s fixed-seed probe repeat at
   3.495×, and §5.2's step-400 fold at 3.596×) are **both on the fixed held-out probe**, which has no
   measured floor: they cleared a floor that does not license them, and they are now unjudgeable. The
   count of failures fell from 23 to 13 — but not one of the ten that left the failing column moved
   into the clearing one.
2. **Every rank number in §5 is unjudgeable, and that is a sharper limitation than "it fails".**
   Eleven of the 25 selections, and 25 of the 56 rows, sit on the fixed held-out probe, an in-run
   training batch, the gate batch or a live checkpoint. The criterion cannot reach them in either
   direction. §6.2 now names the one run that would change this.
3. **§4.6's counts have a fourth defect, and it is the paper's own.** §4.6 refuses its counts on
   three grounds — n = 6, D2 s44's 1.4 sampling SDs, §4.6a's coordinate choice. To those add: **not
   one of the six pairs those twelve metric rows are scored on is resolvable** under canonical R1
   (row 29). The table is a record of orderings read off differences that instrument cannot see —
   though see finding 4, because that verdict is not statistic-invariant either.
4. **The resolvability verdict is itself under-determined, which is §4.5's thesis applied to the
   criterion rather than to the ordering.** Scored against each statistic's own floor, the number of
   §4.5(a)'s six arm pairs that are resolvable is **0 under R1, 1 under R2, 1 under R3, 2 under PR
   and 2 under PR_rownorm** (rows 20–24); scored against each view's own floor it is **0 on
   `wsi_biology` and 6 on each of `rna_biology` and `full_biology`** (rows 26–28). §4.5 says the
   *verdict* flips with the statistic, the block and the view. So does the question of whether there
   is a verdict to be had.
5. **§4.3's headline is the one claim the criterion cannot be applied to at all**, and it is now
   unjudgeable for two independent reasons. The compared quantity is a *spread* (6.05× against
   1.18×), and this project has measured floors on rank **levels** — now for ten statistics, three
   views and two blocks — and never on a spread; and the block is the in-run training batch, which
   has no floor either (rows 16–18). §4.3 survives on the size of the gap between the two arms'
   step-200 levels and on five seeds per arm, not on a ratio that clears anything.
6. **The predeclared violation criterion of §4.7.1 could only ever return unresolvable violations.**
   Its fold threshold was **2.0×** — 1.65× *smaller* than the floor and below §4.2's 2.10–3.75× seed
   band (row 31). Both violations the 66-pair scan returned are inside a floor: H44/I43 at 3.73×
   inside the seed band by 0.5% (already stated in §4.7.4), and **P44/I43 at 3.07×, inside the
   3.295× floor, which appears only in `P2_FIGURES.md` F6(d) and carries no floor verdict there**
   (row 33). That is the sixth instance, and it is the one this audit was built to find.
7. **Two claims in §5.2 read an ordering off a 1.04× rank difference** (rows 48, 51). The
   staleness falsification's second measurement — *"the best-agreeing arm does not have the best
   rank"* — compares **6.89 against 6.65**; and the turnover criterion's discriminating prediction —
   *"the discriminating one inverted"* — compares **3.67 against 3.53**. Both are on the probe block
   and are therefore unjudgeable against any measured floor, which does not rescue them: an ordering
   read off a 4% difference needs a floor smaller than 4% to be worth reading, and no floor this
   project has measured, on any statistic or any view, is that small. **Both claims survive when
   restated as what they actually are**: key-to-encoder agreement varies 2.06× while rank is
   *indistinguishable*, and the predeclared bar (*"fails ≤ 3.5"*) was simply **not met** — a
   predicted effect absent, not an inversion. Row 52 is the same two numbers doing exactly that.
8. **A like-for-like floor for §5.4's regime does exist, at n = 2, and it is concordant.** §5.4 and
   §6.2 both record that no floor has been measured for the `programme_free` / held-out-probe
   regime, and this audit confirms that none can be recovered from the exports. Two runs in
   `~/e0_run/d1_diag/` are the same configuration at the same seed — `ablate_decorr0.04` and
   `mseed_m0.999_s42` share momentum, decorrelation, capacity, learning rate, seed and step-0 state,
   and `d1_momentum_probe.py` runs a constant learning rate with no schedule keyed to the step budget
   — and they span **1.066×** at step 400, at most **1.128×** over the eight shared logged steps
   (row 56). **This may not be quoted as a floor and is not**: it is a *pair*, and §4.1's own
   argument is that any pair drawn from four concordant repeats spans at most 1.028×, so a floor
   measured that way would license everything in this paper.

**Exemptions, stated.** Five rows are exempt and each says why. §4.8's dose–response (rows 34–35):
the representation has no stochastically trained parameters — the only fitted step is a
deterministic per-level PCA — so there is no retraining and §3.5 does not apply; the claim rests on
monotonicity over seven nested levels and on the *ratio* of two fractional changes read from one
representation through one instrument. §5.1's instance 2 (row 41): 5.81 on a 16-patient *train*
batch against a frozen 64-key queue and ~1.8 on 282 *held-out* patients at cohort scale are **not
two configurations and not the same block**, so §3.1's own rule forbids forming the ratio; the claim
rests on the gate's binary pass criterion and on the cohort-scale arm being independently collapsed
(see §5.1). §4.9's decorrelation instance (row 36) and its "16/16" instance (row 37) are excluded on
provenance and on the statistic being a hard numerical rank, not on the floor — and the hard rank now
has a measured floor on the exported block, **1.000×**, which is a second reason nothing may be read
off it. **An exemption stated is fine; an omission is not**, which is why the nuisance measurements
(rows 11–18) and the collapse-regime readings (rows 38–40) are in the table too rather than left out
as obviously inapplicable.

**What this costs and what it buys.** It costs the paper any remaining licence to read a rank
ordering anywhere in it, and it costs more than the first version of this table admitted: **13 of the
25 selections between candidate configurations are inside a floor their own statistic and block
license, 11 cannot be judged at all, and the single one that clears is RankMe as published, not
ours.** What it buys is three things. §1.3's rule is now applied to this paper exhaustively and
mechanically rather than opportunistically; the next comparison added to the draft is checked by a
test rather than by a reviewer; and the criterion itself is no longer being applied outside the scope
of the one measurement that licenses it — which is the objection this paper makes to RankMe, made to
us, and it had been true of us for as long as this subsection existed.

### 4.2 Where rank's variance lives: 34.5% arm, 65.5% training seed — and the channel is 98.0% arm

Twelve artifacts, four arms (D2-H, D2-I, D1-P, D1-F) × three seeds. The seed changes **nothing** about
objective, architecture, data, split or schedule. Rank-type metrics are decomposed on the log scale
(they are multiplicative and span 6.4–34.1); the canonical-correlation ground truth on the raw scale.

| quantity | SS_arm | SS_seed | **arm share** | F(3,8) |
|---|---:|---:|---:|---:|
| **canonical effective rank (residualised)** | 1.3047 | 2.4762 | **34.5%** | **1.41** (n.s.) |
| RankMe as published (raw, uncentred, ε = 1e-7) | 0.9772 | 2.3817 | 29.1% | 1.09 (n.s.) |
| **ground truth: held-out top-CCA, 40 untrained targets** | 0.0353 | 0.0007 | **98.0%** | **128.20** |

*Provenance: `NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §4.1; scripts
`~/e0_run/p2_competing_metrics.py`, `p2_necessity_and_variance.py`; outputs `~/e0_run/P2_METRICS_D2.json`,
`P2_METRICS_D1.json`. Ground truth reproduces `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json` to 4 d.p.
and the recovered Δ values (−0.1325 / −0.1089 / −0.1226) exactly. All 12 artifacts frozen; CPU only,
thread-capped.*

*Independent reproduction, and a near miss worth recording.* A workspace-drift audit
(`NOTEBOOK_ENTRIES/WORKSPACE_DRIFT_AUDIT_ALL_20260803T2359Z.md`) established that the workspace these
scripts ran in carried a `spectral.py` **predating the rank canonicalisation entirely** — no
`CANONICAL`, no `RANK_VARIANTS` — so every number in this section was computed by a different function
object from the one now in the repository. Recomputed against a workspace verified current (0 files
differing, 0 missing), **all five values above reproduce exactly**: 34.5%, F = 1.41, 29.1%, 98.0%,
F = 128.20. They reproduce because the pre-consolidation implementation is numerically the same
function as the new `CANONICAL` — the consolidation added named variants without moving the default —
and `residualise.py` is byte-identical across every workspace, so the channel side was never exposed.
**That is luck rather than design, and it is stated as such**: had the consolidation chosen a different
default, every number here would have shifted silently and the reproduction check is the only thing
that would have caught it. It is also the reason `v2/tests/test_effective_rank_canonical.py` now
asserts object identity across call sites and fails the build if a second definition reappears.

**Two-thirds of the variation in effective rank across these artifacts is training-seed nuisance. Two
percent of the variation in the information channel is.** The arm effect on rank is not significant at
all; the arm effect on the information those same artifacts carry is overwhelming. Note that this holds
for **RankMe as published** as well as for the centred canonical statistic, and slightly more strongly —
so it is not an artifact of our centring deviation.

Per arm, holding everything but the seed fixed:

| arm | canonical effective rank across 3 seeds | fold | held-out top-CCA across 3 seeds | fold |
|---|---|---:|---|---:|
| D2 H (Hallmark) | 23.387 / 28.771 / 9.143 | **3.15×** | 0.6126 / 0.5970 / 0.5983 | 1.026× |
| D2 I (perturbation basis) | 14.868 / 34.117 / 9.105 | **3.75×** | 0.4800 / 0.4882 / 0.4757 | 1.026× |
| D1 P (`programme_only`) | 29.381 / 24.673 / 11.115 | **2.64×** | 0.6117 / 0.6198 / 0.6087 | 1.018× |
| D1 F (`programme_free`) | 13.418 / 7.600 / 6.394 | **2.10×** | 0.5412 / 0.5336 / 0.5126 | 1.056× |

**This is inside RankMe's own reserved scope.** Three seeds of one arm **are** *"different runs of a
given method"*. Across them the proxy moves **2.10–3.75×** while the quantity it is a proxy for moves
**1.8–5.6%**. The between-arm gaps the proxy is being asked to resolve are 12–24% relative. **The
nuisance-induced range of the proxy exceeds the signal it must resolve by roughly an order of magnitude,
in the one regime its authors reserve for it.**

The sharpest single instance is strictly within one arm — no arm contrast at all, only a seed:

> **D2 arm H, seed 44 against seed 43: 3.15× apart in effective rank (9.143 against 28.771), and
> 0.0012 apart in the molecular channel (0.5983 against 0.5970) — with the *lower*-rank run marginally
> ahead.** Same objective, same target table, same data, same split, same schedule, same architecture.

That single line is the paper's cleanest object. It contains no arm effect to argue about, and it is
inside RankMe's scope by RankMe's own sentence.

### 4.3 The floor is a property of the arm — and, §4.1a adds, of the statistic and of the view

The reproducibility envelope is not a constant of the metric that could be measured once and reused. At
the **same global training step**, on the **same configuration**, with only the seed differing:

| arm, D1-B, global step 200 (epoch 11) | seed 42 | seed 43 | seed 44 | seed 45 | seed 46 | spread |
|---|---:|---:|---:|---:|---:|---:|
| `programme_free` | 7.545 | **45.646** | 12.194 | 18.881 | 22.518 | **6.05×** |
| `programme_only` | 110.765 | 110.879 | 111.078 | 94.952 | 112.078 | **1.18×** |

*Provenance: `~/e0_run/d1_v2/d1_{p,f}_seed{42,43,44}/train_metrics.jsonl` and
`~/e0_run/d1_seeds4546/d1_{p,f}_seed{45,46}/train_metrics.jsonl`, key
`train_rank_tripwire_observed`, epoch 11; distilled with per-file SHA-256 into
`v2/research/rebase/p2/figures/data/extracted/F3_TRIPWIRE_STEP200_R3_n5.json`. **Statistic
R3** — the tripwire statistic. These are in-training measurements whose states were never saved, so they
are `[NOT RECOMPUTABLE]` under R1 and are not compared with any R1 number in this paper.
`NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md` §6.
**Seeds 45 and 46 landed after that entry was written**; on the original three seeds alone the
`programme_only` spread is 1.003× and the `programme_free` spread is unchanged at 6.05×. Earlier
versions of this table showed only the three, which understated the supervised arm's spread by a factor
of 60 and is corrected here.*

The supervised arm reproduces to **18%** across five seeds and to **0.3%** across the original three;
the contrastive arm spans **6.05×** either way. **Rank's reproducibility is a property of the arm being
measured, and it is worst exactly where the arm is interesting.** A reproducibility envelope measured on
a well-behaved configuration does not transfer to the one being diagnosed, and the one being diagnosed
is the configuration a practitioner is looking at when they reach for rank.

**An earlier version of this section added "and not of the statistic". That is now measurably wrong
and is withdrawn.** §4.1a measures the retraining floor for ten statistics on one block and finds it
running from 1.000× (hard numerical rank) to 3.295× (canonical R1) — a 3.3× spread between statistics
on the *same five runs* — and for the same statistic on three co-trained views, where it runs from
1.019× to 3.295×. The reproducibility of a rank number depends on the arm, on the statistic and on the
view, and the paper's argument does not need the stronger claim: **three axes of dependence make
"calibrate the envelope once, then use it" less workable, not more.**

This also disposes of the most natural repair: "calibrate the envelope once, then use it". There is no
one envelope, and §4.1a shows there is not even one per arm. On this stack the same statistic, at the
same step, is reproducible to eighteen percent on one arm and unusable on its sibling; on the same five
retrains of one arm it is reproducible to 2% under stable rank and unusable under canonical R1; and on
those same five artifacts it is reproducible to 2% on the `rna_biology` view and unusable on
`wsi_biology`. A practitioner who wants to use a rank difference has to calibrate the arm, the statistic
and the view they are actually using.

### 4.4 Defeater check — is the metric fine and the estimator bad?

The most obvious way this paper's result could be an artifact is that **our measurement of rank is
poor**, not that rank is unstable: a badly conditioned SVD, an unlucky definition, too few patients, a
centring choice the source papers do not make. Given that we report rank spanning 9.1 to 34.1 across
seeds, this is the first objection a referee will raise, and the exemplar for this genre (Leavitt &
Morcos, ICLR 2021, §4.2) devotes a section to the analogous check. Four independent measurements rule
it out.

**(1) Measurement noise on a fixed trained model is negligible.** Patient-subsampling at 80%, 40 draws,
per artifact — the sampling variability of the statistic **given the model**:

| pair | eff. rank, arm A | eff. rank, arm B | between-arm gap | gap / sd |
|---|---:|---:|---:|---:|
| D2 s42 | 23.325 ± 0.072 | 14.834 ± 0.037 | 8.491 | 105.3 |
| D2 s43 | 28.657 ± 0.115 | 33.912 ± 0.111 | 5.255 | 33.0 |
| **D2 s44** | **9.132 ± 0.021** | **9.093 ± 0.019** | **0.039** | **1.39** |
| D1 s42 | 29.244 ± 0.102 | 13.393 ± 0.027 | 15.850 | 150.3 |
| D1 s43 | 24.584 ± 0.069 | 7.586 ± 0.016 | 16.997 | 238.2 |
| D1 s44 | 11.098 ± 0.031 | 6.387 ± 0.011 | 4.710 | 145.1 |

*Provenance: `NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §3; output
`~/e0_run/P2_METRICS_ALL_SUBSAMPLED.json`.*

The standard deviation is **≈0.1 on a rank of 25**. Effective rank is a *precise* measurement of the
matrix it is handed. **The instability documented in §4.1–§4.3 is entirely in training, not in
estimation** — and that distinction must be made explicitly, because §4.1 can otherwise be misread as an
estimator problem. It also has a consequence for §4.6: at D2 seed 44 the between-arm gap is **1.4
sampling standard deviations**, an unresolvable tie, against a ground-truth gap of +0.1226 whose paired
bootstrap CI excludes zero. That "hit" is not a hit.

**(2) The readout path is deterministic.** Re-exporting a surviving checkpoint reproduces its recorded
channel to five significant figures (0.58612 against 0.5861) and gives a stable rank (8.6809); the
3.295× spread of §4.1 appears only when the model is retrained (§3.5).

**(3) At a fixed seed and a fixed short horizon, rank is tight.** Three repeats of a controlled 200-step
probe with **identical inputs**, real streaming batches, live queue: m = 0.999 gives **7.15 / 6.92 /
7.25**, a relative spread of **4.7%**; m = 0 gives 1.80 / 1.46 / 1.98 (30%), with an empty band from 1.98
to 6.92. *Provenance: `NOTEBOOK_ENTRIES/rank_probe_repeat_variance_20260804T0900Z.md`; logs
`~/e0_run/d1_diag/probevar_*.log`.*

**These repeats hold the seed fixed and vary only stack non-determinism; the 6.05× row of §4.3 is at the
same 200 steps on the same objective and varies the seed.** The two are not in conflict and the
difference is the point: rank at step 200 is reproducible to 4.7% when the seed is held, spans 6.05×
when it is not, and spans 1.003× on the sibling arm across those same three seeds (1.18× once seeds 45
and 46 are included, §4.3). **The variance is
neither intrinsic to the measurement nor merely accumulated over training: it is seed sensitivity that
is specific to the arm.**

*The honest limit on (3):* three repeats constrain the *typical* spread and essentially not the *tail*.
On the same stack, the liveness gate diverges at 2/8 = 25%; `P(0 in 3 | p = 0.25) = 0.42` and the exact
upper 95% bound from 0/3 is `p ≤ 0.63`. The planned design was ten repetitions and was cut to six, then
three per condition, because the ten-way launch exhausted GPU memory. That is a real weakening and is
recorded rather than glossed.

**(4) The definitional choices do not carry the result.** The canonicalisation changed **no** published
value: every surviving instance reproduces to the digits published (§4.9). The absolute-to-relative
tolerance change moved no number (max relative difference `0.000e+00` over 68 artifact × block
combinations). **Uncentred R1 is within 0.5% of centred R1 at the median** over those same 68
combinations, and a faithful RankMe (uncentred, ε = 1e-7) reads **23.391** where the canonical statistic
reads **23.387** — so the "you evaluated a centred variant, the published metric does better" objection
is answered numerically: on these artifacts the two are the same number to four significant figures.

**What (4) does *not* answer.** The choice of *statistic* and of *block* does change the between-arm
**verdict**, on a third of the pairs. That is §4.5 — where the table has itself been corrected once,
for exactly the reason this paper is about. It is a separate objection which we do not dismiss; we
adopt it.

**One thing this section cannot establish.** We have two configurations and one stack. We do not know
whether this variance is a property of effective rank, of this hardware's non-determinism, of this
architecture, or of the 40-epoch schedule. The controlled repeat design that would begin to separate
them — N retrainings of one configuration with rank and channel measured on each — **has now been run
at N = 5** and is §4.1's floor; it establishes that the variance is in retraining rather than in
estimation, and it still cannot attribute retraining variance to the metric rather than to the stack,
because there is only one stack. A labelled downstream probe remains the more valuable missing
measurement (§6.2).

### 4.5 The verdict is under-determined: it flips with the statistic, the block and the measured view

Even setting reproducibility aside, "which arm has the higher effective rank" does not have one answer
on our data, while "which arm carries more molecular information" always does.

**(a) It flips with the statistic — and this table has itself been corrected once, which is the
section's own point made at our expense.** All statistics computed on the same 12 artifacts, same
residualisation, same ground truth; "OK" means the higher-rank arm is the one that carries the larger
channel:

| statistic | D2 s42 | D2 s43 | D2 s44 | D1 s42 | D1 s43 | D1 s44 | D2 | D1 | ALL |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **R1** canonical — Roy & Vetterli, order 1 on σ | OK | **MISS** | **OK** | OK | OK | OK | 2/3 | 3/3 | 5/6 |
| **R2** canonical — `d1_audit.py`, order 2 on σ, `(Σσ)²/Σσ²` | OK | **OK** | **MISS** | OK | OK | OK | 2/3 | 3/3 | 5/6 |
| **R3** canonical — R2 after L2-normalising rows | OK | **MISS** | **OK** | OK | OK | OK | 2/3 | 3/3 | 5/6 |
| *PR* — order 2 on the **eigenvalue** distribution, `(Σσ²)²/Σσ⁴` *(published in error as "R2")* | OK | OK | MISS | **MISS** | OK | OK | 2/3 | 2/3 | 4/6 |
| *PR_rownorm* — PR after L2-normalising rows *(published in error as "R3")* | OK | OK | MISS | OK | OK | OK | 2/3 | 3/3 | 5/6 |

*Provenance: `v2/research/rebase/p2/p2_rank_variants.py`, vendored at commit `7b37dce` and re-run on a
workspace verified byte-equal to HEAD before execution (402/402 files by git blob SHA-1);
`NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260804T0255Z.md` §0. Originally published from
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §4.3. The first three rows
are the canonical `RANK_VARIANTS`; the last two are the statistics the pre-vendoring script actually
computed under the labels "R2" and "R3", retained under their true names because the numbers are correct
values of a different function. Note that the `PR` row is **identical, cell for cell, to §4.6's
"participation ratio" row** — that identity is the arithmetic signature of the substitution and is why
the two must not be counted as independent evidence.*

#### The correction to this table is itself the paper's result, and we present it that way

Earlier versions of this table labelled rows 2 and 3 **R2** and **R3** — the statistics of
`d1_audit.py` and `d1_geometry_probe.py`. **They were not.** The script that produced them carried its
own inline definitions and computed `(Σσ²)²/Σσ⁴`, the order-2 Hill number of the **eigenvalue**
distribution, where `d1_audit.py`'s R2 — and therefore `RANK_VARIANTS["R2"]` — is `(Σσ)²/Σσ²`, the
order-2 Hill number of the **singular-value** distribution. They are different statistics. The rows are
relabelled **PR** and **PR_rownorm**, which is exactly what they are; the numbers themselves are correct
values of that function and reproduce cell for cell.

**This is a fourth statistic living under a name §3.1 spends its length disambiguating — and it was
inside our own analysis code, in the section of this paper that argues the name is unreliable.** We
report it as evidence rather than as an erratum, because almost nothing else we could measure would
support the claim as directly.

**How it survived, which is the part that generalises.** The substitution was invisible to review — the
label said R2 and the number was plausible for R2. It was invisible to the test suite, because no test
compared the analysis script's rank to the repository's. It was invisible to the authors through two
drafts and a full recomputation pass. **It was found only when the traceability rule was enforced —
when the script that produced the numbers was vendored into the repository and rewritten to call the
single canonical implementation instead of its own.** Nothing about care, review or expertise caught
it; a mechanical provenance requirement did. **That is the recommendation this paper can actually
make: against this class of error the defence is mechanical provenance, not diligence.**

**Two details that make it usable rather than confessional.** First, the arithmetic fingerprint: the
`PR` row is **identical, cell for cell, to §4.6's "participation ratio" row**. That identity is what a
reader can look for to detect the same substitution in their own tables — two differently-named rows
that agree exactly are one statistic reported twice, not two statistics agreeing. Second, the direction
of the damage: the correction moved one count **against** us (3 of 6 → 2 of 6) and withdrew one
qualification **for** us (D1's "only 2/3 under R2" becomes 3/3, removing a hedge from a result that
already went against this paper). **An error scatters; a bias does not.** We note that not as a defence
but because it is the check a referee should run on any self-reported correction.

**What moves as a result, and it moves against us.** Re-measured with the canonical `RANK_VARIANTS`:

- **The headline count falls from 3 of 6 pairs to 2 of 6** (D2 s43 and D2 s44; **D1 s42 no longer
  disagrees**). The claim that the statistic changes the verdict survives; its size does not.
- **Canonical R2 scores D1 3/3**, where the mislabelled row reported 2/3. The consequence for §4.7.3 is
  recorded there: the qualification *"D1 is only 2/3 under the statistic an earlier draft nominated"*
  **does not survive**, and D1's confirmation of RankMe is correspondingly less qualified than the
  previous draft allowed.
- **Canonical R3 swaps its two D2 verdicts.**
- **§4.5(b)'s R2 and R3 *levels* are unaffected** — they were computed with the canonical function all
  along and reproduce exactly.

**A second, smaller error in the same section, also corrected.** §4.5's provenance note previously
explained the discrepancy between the two source entries' R3 rows as *"one is computed on the raw block
and the other on the residualised block"*. **That was not the reason.** Both are on the residualised
block; the difference was the statistic. The R3 *levels* quoted in (b) below (H43 14.746 / I43 15.915;
H44 6.564 / I44 6.302) are the canonical R3 and reproduce exactly — it was only (a)'s *verdicts* that
came from `PR_rownorm`.

All of PR, PR_rownorm, R1, R2 and R3 have been called `effective_rank` somewhere in this repository or
in the published literature, and the pairwise disagreements above are between functions that share
that name.

**(b) It flips with the block — 2 of 3 D2 seeds, for R3.** Raw versus confound-residualised, the choice
no version of this project stated:

| seed | statistic | H (raw) | I (raw) | higher | H (resid.) | I (resid.) | higher |
|---|---|---:|---:|:--:|---:|---:|:--:|
| 43 | R1 | 24.674 | 28.800 | I | 28.772 | 34.117 | I |
| 43 | R2 | 11.720 | 11.111 | **H** | 13.227 | 12.972 | **H** |
| 43 | **R3** | 11.720 | 11.111 | **H** | 14.746 | 15.915 | **I** |
| 44 | R1 | 8.447 | 8.313 | H | 9.143 | 9.105 | H |
| 44 | R2 | 5.385 | 5.449 | I | 5.733 | 5.815 | I |
| 44 | **R3** | 5.385 | 5.449 | **I** | 6.564 | 6.302 | **H** |

(Arm H carries the larger channel in all three seeds, so "higher = H" is the ordering a rank rule needs.)

**(c) It flips with which co-trained view of the same model you measure — 2 of 6 pairs.** Each artifact
stores three co-trained views: `wsi_biology` (the canonical readout view), `rna_biology`,
`full_biology`.

| pair | information winner (wsi / rna / full) | rank winner (wsi / rna / full) | rank verdict stable? |
|---|---|---|:---:|
| D2 s42 | H / H / H | **H / I / I** | **NO** |
| D2 s43 | H / H / H | I / I / I | yes (and wrong) |
| D2 s44 | H / H / H | **H / I / I** | **NO** |
| D1 s42 | P / P / P | P / P / P | yes |
| D1 s43 | P / P / P | P / P / P | yes |
| D1 s44 | P / P / P | P / P / P | yes |

**The information ordering is identical under all three views for all six pairs. The rank ordering is
not.** Aggregated over all 18 (pair × view) comparisons the canonical statistic is right 11/18
(p ≈ 0.48); restricted to D2, the matched-target-table contrast, it is right **2/9**.

*Caveat, stated plainly:* `rna_biology` → RNA-derived pathway targets is partly circular, and its
absolute CCA (0.79–0.85) must not be read as a clean image→molecular channel. That is why the canonical
readout is `wsi_biology`. The **rank** measurements on the other views are unaffected by that
circularity, and the instability of the rank *verdict* across views stands on its own — but the 2/9
count inherits the caveat and is not quoted as a rate.

*Provenance for (a)–(c): `NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md`
§4.2, §4.3 and `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`
§4. The two entries' R3 rows differ because one is computed on the raw block and the other on the
residualised block; the reconciliation is stated in the second entry's §7 and the block is named in every
row above. Outputs `~/e0_run/P2_ROBUSTNESS.json`, `~/ws_rank/RANK_RECOMPUTE.json`.*

**What (a)–(c) establish, and what they do not.** They do not show that rank is wrong; they show that
**a rank verdict is not a well-defined object until three implementation choices are stated, and that
those choices are worth more than the between-arm difference they are used to adjudicate.** The
information verdict, computed through a single documented pipeline, is invariant to all three. This is
the paper's own thesis applied to the paper's own instances, and it is reported here rather than in the
limitations for that reason.

### 4.6 Rank as a selection rule, against the published alternatives — reported, and underpowered

Six matched pairs, one rule: *pick the arm with the higher metric; is it the arm that carries more
molecular information?* Ground truth is the held-out top canonical correlation against the **40 targets
neither arm was supervised on** (`heldout_pathway` + `immune_tme` + `tumour_state`), reproducing
`d2_readout.py`'s residualisation, seed and component budget exactly.

| pair | arm A | arm B | A | B | Δ(A−B) | winner |
|---|---|---|---:|---:|---:|:---:|
| D2 s42 | H42 | I42 | 0.6126 | 0.4800 | +0.1325 | H |
| D2 s43 | H43 | I43 | 0.5970 | 0.4882 | +0.1089 | H |
| D2 s44 | H44 | I44 | 0.5983 | 0.4757 | +0.1226 | H |
| D1 s42 | P42 | F42 | 0.6117 | 0.5412 | +0.0705 | P |
| D1 s43 | P43 | F43 | 0.6198 | 0.5336 | +0.0863 | P |
| D1 s44 | P44 | F44 | 0.6087 | 0.5126 | +0.0961 | P |

**Sign convention, stated because it differs from §4.7.2's by design.** Δ here is `A − B` with the
better arm written first, so every value is positive and the "winner" column carries the meaning.
§4.7.2 uses the **predeclared** direction `channel(programme_free) − channel(programme_only)`, under
which the same three D1 gaps read **−0.0705 / −0.0863 / −0.0961** and their bootstrap intervals are
negative. **Same numbers, opposite sign, two tables** — flagged here rather than left for a reader to
trip over.

| metric | D2 s42 | D2 s43 | D2 s44 | D1 s42 | D1 s43 | D1 s44 | D2 | D1 | ALL | exact 2-sided binomial *p* |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|---:|
| canonical effective rank (raw) | OK | MISS | OK | OK | OK | OK | 2/3 | 3/3 | **5/6** | 0.219 |
| canonical effective rank (resid.) | OK | MISS | OK | OK | OK | OK | 2/3 | 3/3 | **5/6** | 0.219 |
| **RankMe (raw, as published)** | OK | OK | OK | MISS | OK | MISS | **3/3** | 1/3 | 4/6 | 0.688 |
| RankMe (residualised) | OK | MISS | OK | OK | OK | OK | 2/3 | 3/3 | 5/6 | 0.219 |
| participation ratio (raw / resid.) | OK | OK | MISS | MISS | OK | OK | 2/3 | 2/3 | 4/6 | 0.688 |
| stable rank (raw) | OK | OK | MISS | MISS | OK | MISS | 2/3 | 1/3 | 3/6 | 1.000 |
| stable rank (resid.) | MISS | OK | MISS | MISS | OK | OK | 1/3 | 2/3 | 3/6 | 1.000 |
| α-ReQ \|α−1\| (raw / resid.) | OK | MISS | OK | OK | OK | MISS | 2/3 | 2/3 | 4/6 | 0.688 |
| **LiDAR (raw)** | MISS | MISS | MISS | OK | OK | OK | **0/3** | 3/3 | 3/6 | 1.000 |
| LiDAR (residualised) | MISS | OK | MISS | OK | OK | OK | 1/3 | 3/3 | 4/6 | 0.688 |

*Provenance: `NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §3, §6;
scripts `~/e0_run/p2_competing_metrics.py`, `p2_selection_rule.py`; outputs `~/e0_run/P2_METRICS_D2.json`,
`P2_METRICS_D1.json`. LiDAR implemented from arXiv:2312.04000v1 §3 Eqs. (1)–(4) with q = 2 and δ swept
over 1e-8…1e0; α-ReQ from the authors' released `fastssl` estimator. Fidelity and the departures from
each source are stated in §2.5.*

**Four things must be said about this table, in this order.** The fourth is §4.6a and is the largest:
the ground truth this table scores against is one coordinate system out of six, and every count below
moves when it changes. Nothing in this section may be read without it.

**First, it is underpowered and must not carry rhetorical weight.** At n = 6, the exact two-sided
binomial gives **p = 0.031 for a flawless 6/6**, 0.219 for 5/6 and 0.688 for 4/6. The design therefore
has just enough power to detect a *perfect* rule and none at all to detect a merely good one. **No
conclusion of the form "metric X is better than metric Y" can be drawn from these six pairs, in either
direction**, and 5/6 against 4/6 is not evidence of anything. This is the main reason the paper rests on
§4.2's variance decomposition, which uses magnitudes and estimates the nuisance term from 8 within-arm
degrees of freedom, rather than on any count here.

**Second, one apparent hit is not one.** D2 s44's between-arm rank gap is 1.4 sampling standard
deviations (§4.4), i.e. unresolvable. Effective rank's honest D2 record is **1 clear hit, 1 clear miss,
1 pair it cannot resolve at all** — not 2/3.

**Third, the answer to "why did you not just use the better metric" is that on our data there is not
one.** LiDAR, the strongest published alternative, scores **0/3 on D2** — choosing the
information-poorer arm every time — and does so at **every δ across eight orders of magnitude**
(1e-8 → 1e0; the absolute LiDAR value moves from 176.1 to 7.0 over that range, but the *ordering*, which
is all a selection rule uses, is invariant). RankMe as published wins D2 3/3 and loses D1 1/3; the
mechanism of its D2 advantage is the **mean-offset direction** its uncentred normalisation retains, and
that mechanism does not survive to D1. α-ReQ's 4/6 is squeezed out of α differences its own authors
would call meaningless: all 12 artifacts have α between 2.6 and 4.8, far outside its "Goldilocks zone",
and `|α − 1|` is our operationalisation, not theirs. **Both halves of the RankMe result must be
reported**, because "you evaluated a centred variant and the published metric does better" is true on
D2 and false on D1.

**Fourth, and largest: all of the above is scored against one target block.** See §4.6a. The 3/3-versus-
2/3 comparison just made between RankMe and canonical effective rank **reverses** if the ground truth is
taken on the perturbation dictionary's own codes or on a plain PCA basis, and canonical effective rank
then reaches a nominally significant 6/6. Neither reading is quotable. The sentences above stand as a
description of what happens on the gene-set block and on no wider domain than that.

### 4.6a The ground truth §4.6 scores against is itself a coordinate choice, and every count in that table moves when it changes

§4.6's marks are scored against **one** ground truth: the held-out channel onto the 40 **gene-set**
targets neither arm trained on. On 2026-08-04 that ground truth stopped being a single number.
`NOTEBOOK_ENTRIES/d2_coordinate_system_result_20260804T0800Z.md` re-scored both D2 arms on every other
molecular target block that exists on this project's disk, and the **−0.12 arm contrast appears on
exactly one of six** — the gene sets. On the dictionary's own 128 supervision codes all six 2,000-repeat
intervals cover zero; every other block sits at or inside the published `random_control` negative
control. Rotating the exam's basis moves the arm contrast by **+0.1227 / +0.1177 / +0.1200**, which is
the entire published effect.

**If the arm contrast is block-dependent, so is §4.6's ground truth, and so is every OK/MISS in it.** We
had not chased that consequence and it is the counterfactual the table needs. Holding all twelve metrics
fixed and swapping only the truth, once per block
(`v2/research/rebase/p2/p2_selection_rule_blocks.py`, reading the same metrics JSON §4.6 reads and the
per-block contrasts from `v2/research/rebase/nature/d2_coordinate_system/out/EXAM_PANEL.json`):

| exam block used as ground truth | Δ s42 | Δ s43 | Δ s44 | arm ordering |
|---|---:|---:|---:|:---:|
| **gene sets, untrained 40 — the published truth** | −0.1325 | −0.1089 | −0.1226 | H H H |
| PBS codes 128 (arm I's own supervision) | −0.0098 | **+0.0088** | −0.0026 | H **I** H |
| PCA basis 128 | −0.0201 | **+0.0049** | −0.0284 | H **I** H |
| gene-label-shuffled PBS 128 | −0.0359 | −0.0057 | −0.0175 | H H H |
| `random_control` gene sets 90 (negative control) | −0.0099 | −0.0280 | −0.0268 | H H H |
| size/spectrum-matched random dictionary 128 | −0.0597 | −0.0132 | −0.0454 | H H H |

**Everything turns on seed 43, and seed 43 is the seed the coordinate system flips.** The D2 selection
count under each truth, for the two rows §4.6 quotes against one another:

| metric | gene sets 40 | PBS codes | PCA basis | shuffled PBS | rand control | rand dict |
|---|:--:|:--:|:--:|:--:|:--:|:--:|
| **canonical effective rank (resid.)**, D2 | 2/3 | **3/3** | **3/3** | 2/3 | 2/3 | 2/3 |
| **RankMe (raw, as published)**, D2 | **3/3** | 2/3 | 2/3 | **3/3** | **3/3** | **3/3** |
| — the same, as ALL / 6 (canonical R1) | 5/6 | **6/6**, *p* = 0.031 | **6/6**, *p* = 0.031 | 5/6 | 5/6 | 5/6 |
| — the same, as ALL / 6 (RankMe) | 4/6 | 3/6 | 3/6 | 4/6 | 4/6 | 4/6 |

*Provenance: `v2/research/rebase/p2/p2_selection_rule_blocks.py` over
`~/ws_p2/out/P2_METRICS_D{1,2}.json` and
`v2/research/rebase/nature/d2_coordinate_system/out/EXAM_PANEL.json`. The `geneset_untrained40` column
reproduces §4.6's published table exactly, in all twelve rows, which is what licenses the other five.
Tested in `v2/tests/test_p2_analysis_scripts.py`.*

Three things follow, and the first two are against us.

**First, every one of the twelve metric rows changes its D2 count.** Not one is stable across the six
blocks. §4.6 already refused to let a 5/6 against a 4/6 carry weight on power grounds; this is a second,
independent reason, and it is worse than the power problem because it is not fixed by more pairs.

**Second, the ordering between the two rows the section names explicitly reverses.** On the published
gene sets, RankMe-as-published beats canonical effective rank on D2, 3/3 against 2/3 — a fact §4.6
reports because it cuts against us. On the dictionary's own codes and on a plain PCA basis the two
**swap**, and canonical effective rank reaches **6/6 overall, exact two-sided p = 0.031** — "significant"
by this very design's own stated bar, produced by nothing whatsoever except the choice of which
coordinate system the exam is written in. **We report that we can manufacture our own significance by
choosing a target block, because that is the strongest available demonstration that no count in §4.6 may
be quoted without its block**, and it disqualifies the favourable reading as firmly as the unfavourable
one.

**Third, the D1 half cannot be re-scored at all, and is not treated as though it had been.** The
coordinate-system work re-scored the two **D2** arms; the D1 arms were never scored against any block
but the gene sets. So the D1 column is held fixed in every row above, and the ALL counts inherit that.
**This is not evidence that the D1 half is block-stable — it is the absence of the measurement**, and it
is named in §6.3 alongside the ground truth's other weaknesses.

**What survives §4.6, and it is less than it was.** The direction — arm H is never behind on any block
tested — is stable across all six. What is not stable is any *count*, any *ranking between metrics*, and
any *p*. §4.6 is therefore reported as a descriptive record of one coordinate system, and the paper's
argument continues to rest on §4.2's variance decomposition, which uses magnitudes rather than a
verdict and does not consult a ground truth at all.

### 4.7 The necessity test, which went against us — and which §4.1's floor has made unresolvable rather than refuted

**This is the result that falsified this paper's previous framing, and it is reported before the
instances that favour us. It stays at that prominence.**

> **Read this before §4.7.2.** The measured retraining floor of §4.1 (3.295× residualised, 3.111× raw)
> is **larger than all three of D1's rank ratios** — 2.190× / 3.246× / 1.738× residualised, 2.02× /
> 3.09× / 1.68× raw. The rank difference this section reads is therefore **inside the noise floor of
> the statistic that produced it**, and *D1 does not resolve whether effective rank tracks information,
> in either direction.*
>
> **That is not a rescue and is not written as one.** The reading was fixed in
> `PREDECLARED_retraining_envelope_20260804T0330Z.md` before the floor existed, in the band the floor
> landed in: *"No D1 rank difference is resolvable. D1 is uninformative about rank; the channel remains
> resolvable; the asymmetry is the finding. Explicitly **do not** claim the necessity result is
> refuted."* Applying it without amendment:
>
> - **The necessity result is NOT refuted.** `programme_only` — the higher-rank arm — carries the larger
>   channel in **3/3** seeds, with patient CIs excluding zero 3/3 and cancer-cluster CIs 2/3. Those
>   intervals are unaffected by anything in §4.1. A comparison inside the noise floor is not evidence
>   **for** rank's reliability; it is equally not evidence **against** it.
> - **It is therefore reported at full strength below and flagged as *not resolvable by this
>   comparison*, not deleted, not softened, and not moved out of the abstract.**
> - **The asymmetry is the finding, and it holds only because both halves are reported.** On the same
>   comparison, the channel difference clears its own instrument and the rank difference does not. A
>   version of this section that quoted the rank half alone would convert a demonstration into a
>   convenience.

#### 4.7.1 What was predeclared

RankMe's defence is that high rank is *"a necessary (but not sufficient) condition for good downstream
performance"*. Under that hedge, **high rank + low information is predicted and is not a
counterexample**; only **low rank + high information** breaks it. D1-B was set up to produce that
configuration: `programme_free` has the lower rank; does it carry a comparable or better molecular
channel?

The criterion was fixed in `p2_necessity_and_variance.py` before the pair list was inspected:

> a pair (lo, hi) is a violation iff `eff_rank(hi)/eff_rank(lo) ≥ 2.0` **and**
> `CCA(lo) − CCA(hi) ≥ 0.0705`.

Both thresholds come from quantities established independently of this analysis: 2.0× is *below* the
2.10–3.75× that the seed alone produces within a single arm (§4.2), so it is a conservative floor; and
0.0705 is the **smallest** between-arm channel gap this project has accepted as real. The reading of
each of four outcomes was committed in
`NOTEBOOK_ENTRIES/PREDECLARED_D1_necessity_test_20260803T2300Z.md`, written before `d1_audit.log` was
opened and before any channel value was seen.

#### 4.7.2 Result: necessity is not violated. It is confirmed — 3/3 on the patient bootstrap, 2/3 on the conservative one

Δ is the **predeclared** direction, `channel(programme_free) − channel(programme_only)`, so a negative Δ
means the lower-rank arm loses.

| seed | rank, `programme_only` | rank, `programme_free` | ratio | channel, `programme_only` | channel, `programme_free` | **Δ** | patient CI₉₅ | **cancer-cluster CI₉₅** | rank ordering correct? |
|---|---:|---:|---:|---:|---:|---:|:---:|:---:|:---:|
| 42 | **29.381** | **13.418** | 2.19× | 0.6117 | 0.5412 | **−0.0705** | [−0.0938, −0.0444] | [−0.0957, −0.0180] | **yes** |
| 43 | **24.673** | **7.600** | 3.25× | 0.6198 | 0.5336 | **−0.0863** | [−0.1186, −0.0522] | **[−0.1386, +0.0006]** | **yes** |
| 44 | **11.115** | **6.394** | 1.74× | 0.6087 | 0.5126 | **−0.0961** | [−0.1314, −0.0618] | [−0.1535, −0.0016] | **yes** |

*Provenance: rank — `~/e0_run/d1_v2/artifacts/d1_{p,f}_seed{42,43,44}.npz`, **canonical R1 on the
residualised held-out `wsi_biology` block**, recomputed 2026-08-04
(`NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md` §5).
Channel and intervals — `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json`, on the 40 targets neither
arm trained on (`heldout_pathway` + `immune_tme` + `tumour_state`) per the preregistration in
`NOTEBOOK_ENTRIES/d1_readout_preregistration_20260803T1700Z.md`; point estimates independently
reproduced by `NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §5.2.
**The bootstrap existed all along and was hidden by the audit chain's stale absolute path**
(`/lambda/nfs/.../d1_v2/artifacts/` against the real `/home/ubuntu/e0_run/...`) — a one-line bug in the
chain, not a fault in the data, and the reason earlier versions of this table carried
`[D1 PAIRED BOOTSTRAP PENDING]`. `D1_PAIRED_BOOTSTRAP.json` (unstratified) is **not** used here and must
not be: it scores all 90 non-control targets, of which 50 are `programme_only`'s own supervision.*

`programme_free` has both the lower rank **and** the lower channel, **3/3 seeds**, in the same direction
as its rank advantage. This is outcome **O2** of the predeclaration: *"Rank vindicated on this pair.
Goes in the paper as a limitation, with the same prominence a confirmation would have had. The broad
claim must then be narrowed or defended on other evidence."* It is, and it has been: the broad claim is
gone and §1.3 is what replaced it.

**Both bootstraps are quoted, and the conservative one is the one to weight.** On the patient bootstrap
the result is decisive in **3/3** seeds. On the **cancer-cluster** bootstrap — which resamples whole
cancer types rather than patients, and is therefore the estimator that survives the objection that our
2,766 test patients are not 2,766 independent observations — it is decisive in **2/3**, with seed 43's
interval touching zero at **+0.0006**. A referee will weight the cluster bootstrap, and so do we.
Reporting only the patient result would be exactly the kind of selective quotation this section exists
to avoid: **§4.6 refuses to let a 5/6 sign count carry weight, and it would be incoherent to then quote
the favourable one of two estimators on the paper's most load-bearing negative result.**

This does not rescue the paper's old claim. A cluster interval that includes zero by 0.0006 in one seed
of three is not evidence that `programme_free` matches `programme_only`; it is evidence that **one seed
is not decisive at the cancer-cluster level**, against a point estimate of −0.0863 pointing the same way
as the other two. The honest statement is: *the higher-rank arm wins in all three seeds; the win is
interval-backed in three of three at the patient level and two of three at the cancer level.* Necessity
is confirmed, slightly less strongly than the patient bootstrap alone would suggest.

**It also trips a preregistered escalation elsewhere.** `~/e0_run/d1_v2/D1_PAIR_MANIFEST.json` records
`"preregistered_prediction": "programme_free >= programme_only on the held-out molecular channel; if
programme_only wins, the collapse story is wrong -- escalate, do not proceed to D2"`. `programme_only`
wins 3/3. That escalation is flagged here and is **not** resolved by this paper.

#### 4.7.3 The deflations we tested, and which of them survive

We looked for reasons this result might not mean what it says, and report the outcome of each including
the ones that failed.

- **"`programme_free` simply collapsed again, so this is the collapse-diagnostic use, not a test of
  necessity."** *Rejected.* `programme_free` is **not** collapsed this time — canonical ranks 6.4–13.4
  against the ~1.7 of the defective earlier run — and its channel (0.513–0.541) sits clearly above its
  own `random_control` (0.443–0.474). This is the "both arms train" branch. The queue fix (§5.2) worked —
  and **this row is the evidence that fix actually rests on**, not its rank ratio, which §5.4 reports as
  inside §4.1's floor like every other rank comparison in this paper.
- **"The negative control passes, so the readout is clean."** *Only in the narrow sense, and the verdict
  must be recorded with its qualification.* Audit check A3 passes **on arm difference**: the
  `random_control` arm gap is **−0.0224 / −0.0072 / −0.0322** with all three CIs spanning zero, so
  **the instrument does not manufacture an arm difference where none should exist**. But the
  **absolute** level on random controls is **0.4425–0.4810**, against real targets' 0.51–0.62 and the
  D2 row-shuffle permutation null of **0.140** — so random gene sets sit at **≈3.2× the floor**
  (3.1–3.4× across the six control readings), and 0.44 is emphatically **not** the null. Two
  independent estimates agree on the size of that: these controls are **77–85%** of their own arm's
  real-target level, and T1.4 found covariate-matched random gene sets reproducing **76–82%** of the
  per-target channel. **The comparator is 0.140 and not the dilution sweep's 0.147** — an earlier
  version of this bullet used 0.147, which is a different permutation procedure at a different *n*; see
  §3.2's footnote. **The correct
  verdict is therefore "the instrument does not manufacture an arm difference; the absolute level is
  high and separately explained", not an unqualified pass** — and any future quotation of A3 must carry
  that sentence. It does not weaken the D1 result, which is a *paired arm difference* and is exactly the
  quantity the control clears; it does mean no absolute channel level in this paper may be read against
  an assumed null of zero, which §3.2 already requires and §6.3 now restates.
- **"`programme_only` wins because the readout is scored on its own supervision."** *Rejected by
  construction and by measurement.* The readout was pre-restricted, per
  `NOTEBOOK_ENTRIES/d1_readout_preregistration_20260803T1700Z.md`, to the **40 targets neither arm was
  supervised on**; the secondary 90-target readout, half of which *is* `programme_only`'s supervision, is
  explicitly excluded from the headline. And re-estimating with `heldout_top_cca` — directions fit on
  half the patients, scored on the other half — **preserves every arm ordering on all three views with
  equal or larger deltas** (§3.2).
- **"The gap is inside the noise."** *Resolved, and it is the qualification that costs the most.* The
  stratified paired bootstrap exists (§4.7.2). At the patient level the gap is decisive 3/3. At the
  **cancer-cluster** level it is decisive **2/3**, seed 43's interval reaching **+0.0006**. So the
  objection does not succeed — the point estimates agree 3/3 and two of three clear the conservative
  estimator — but it is not dismissed either, and the paper states the 2/3 wherever it states the 3/3.
- **"It depends on the rank statistic."** ***Withdrawn — this qualification does not survive, and its
  withdrawal costs us.*** Earlier versions of this section reported that D1 scores 3/3 under R1 and R3
  but only **2/3 under R2**, the statistic an earlier draft nominated for exactly this table, and
  offered *"2–3 of 3 depending on which function you call effective rank"* as the honest summary. **The
  row that produced the 2/3 was not R2.** It was the order-2 Hill number of the *eigenvalue*
  distribution, mislabelled (§4.5a). Re-measured with the canonical `RANK_VARIANTS`, **D1 scores 3/3
  under canonical R2 as well.** So the necessity result is *less* qualified than the previous draft
  claimed, not more: the higher-rank arm wins in all three seeds under every canonical statistic we
  compute. The only surviving qualification on D1 is the interval one — 3/3 on the patient bootstrap,
  **2/3 on the cancer-cluster bootstrap** (§4.7.2).
- **"The rank gap is ~9×, so this is a large-gap regime."** *Rejected; the 9× figure was ours and it was
  wrong.* The "12 versus 111" numbers are the **in-run tripwire** (R3, training batches, step 200), not
  the quantity the channel is computed on. Measured where the channel is measured — held-out patients,
  epoch 40, exported artifact — the gap is **1.74–3.25×**. Stating "9× lower rank" in the paper would
  have been a false statement, and it was caught by the predeclaration rather than by review.

#### 4.7.4 What this costs the paper, and what it does not

It costs the paper the claim it was previously organised around. *"Effective rank does not track
information content"* is not supported by our own best-matched three-seed experiment, and every version
of this draft asserting it has been withdrawn.

It does not touch §4.1–§4.5, and the reason is worth stating precisely: **the rank ratios in the table
above (1.74×, 2.19×, 3.25×) are not larger than the within-arm seed range (2.10–3.75×) of §4.2, and
since 2026-08-04 all three are also inside the *measured* retraining floor of §4.1 (3.295×).** So even
in the instance where rank gets the answer right, the size of the gap it is reading sits inside two
independently measured noise floors — one estimated across seeds, one across identical repeats at a
fixed seed. **A metric can be right on average and still be unusable for a single comparison**; that is
the whole content of this paper's surviving claim, and D1 is an example of it rather than a
counterexample to it.

**The strongest form of the objection, stated against ourselves.** A referee may say: *"you have made
your own inconvenient result disappear by measuring a noise floor wide enough to swallow it."* Three
things answer that, and none of them is that the result went away.

1. **The reading was predeclared before the floor was measured**, in a document that named this exact
   outcome as *"suspiciously convenient for us"* and foreclosed the rescue in advance. The floor landed
   in the `> 3.09×` band whose reading had already been written.
2. **The floor was measured on the arm that makes it smallest.** `programme_only` is the stable arm
   (§4.1, §4.3). A floor measured on `programme_free` would be larger, not smaller, so the objection
   would have to argue we chose the *less* favourable arm on purpose.
3. **The result is not deleted and its intervals are untouched.** `programme_free` still loses the
   channel 3/3 with patient CIs excluding zero. What §4.1 removes is the licence to read the *rank*
   half of this comparison, in either direction — which costs us the ability to point at D1 as evidence
   that our own thesis fails there, and equally denies anyone the ability to point at it as evidence
   that rank works.

**And the same argument applies to us, which we state before it is put to us.** The sharpest available
objection to this paper is that it holds rank to a strict standard — *your difference must exceed your
own noise floor* — and holds its own effect to a loose one. So: **the smallest D1 arm difference we
report (0.0705) is roughly half the real-versus-random-control margin of 0.139** (real targets
0.51–0.62 against random controls 0.4425–0.4810, §4.7.3; the margin is 0.6117 − 0.4728 = **0.1389** on
seed 42, the seed that supplies the 0.0705). That is §4.1's envelope argument arriving from
the channel side rather than the rank side, and it points at us. The effect we measure is *small
relative to its own instrument's floor*, and a reader is entitled to hold that against our channel
readings exactly as we hold the 3.295× floor against rank.

Three things distinguish the two cases, and none of them makes our margin large:

1. **The channel difference is a paired within-run difference and the rank difference is a level**
   (§3.5). The floor that matters for a paired difference is the spread of that difference, which is
   0.024 on a mean of −0.121 across the D2 seeds and same-signed 3/3; the floor that matters for a level
   is 3.295×, which **all seven** arm differences fail to clear (§4.1). That asymmetry is the practical
   heart of the paper and it is not an accident of which quantity we like better; on the five identical
   retrains that produced the floor, the two quantities were measured on the *same runs* and moved
   3.295× against 1.055×.
2. **The channel difference is interval-backed and the rank difference is not** — 3/3 at the patient
   level and 2/3 at the cancer level here, 3/3 on both in D2. No rank difference in this paper carries
   an interval at all, because no one computes one.
3. **The variance decomposition does not depend on any margin.** §4.2's 34.5% / 98.0% split is a
   statement about where each quantity's variance lives across 12 matched artifacts, and it would read
   the same if every absolute level were halved.

**What none of that buys is a claim that our effect is comfortably above its floor. It is not**, and
§6.3 records this as the paper's second-largest exposure after the single-stack limitation.

**One violation of necessity does exist in our data, and we deliberately do not lead with it.** Scanning
all 66 ordered pairs of the 12 artifacts against the predeclared criterion, two violate. The usable one
is **H44 against I43**: 3.73× lower effective rank (9.143 against 34.117) carrying **+0.110** more
molecular channel — a gap comparable to the headline D2 arm effect, same architecture, cohort, schedule
and modality pair, differing in arm and seed. It is a genuine low-rank/high-information instance and it
is the only configuration RankMe's hedge cannot absorb. **And §4.1's rule applies to it too, which we
state first rather than leaving it to the scope objection.** Its 3.73× clears the 3.295× retraining
floor by 13%, but the comparison varies the **seed**, and the nuisance range that governs a seed-varied
comparison is §4.2's within-arm **2.10–3.75×** — which 3.73× sits inside. So this instance is not a
resolvable rank difference either, by the standard this paper applies to everyone else's, and that is a
second reason it is supporting rather than load-bearing. **But it is a cross-arm, cross-seed comparison,
and RankMe reserves itself to "different runs of a given method"; a referee will argue H44 and I43 are
not that.** It is presented here as a *supporting* observation with its scope stated, never as the
load-bearing one — and, per §2.2, it is in any case **partially pre-empted**: Aldeneh et al. (ICASSP
2025) have already published *"lower-ranked layers can outperform higher-ranked ones"*. **A version of
this paragraph that presents low-rank/high-information as novel is not submittable.**

### 4.8 A dose–response magnitude miscalibration, corrected against ourselves

*(This instance does **not** contradict RankMe as stated. It contradicts the informal practice of
reading rank as a representation-health indicator.)*

Patch bags were contaminated with same-cancer, different-patient tumour patches at seven nested levels.
The representation is `concat(mean, std)` over frozen H-Optimus-0 tokens with **no fitted parameters**,
so nothing here is a different training run and §3.5 does not apply. Both quantities are read from the
same representation at each level, through the same instrument.

| requested d | achieved d | adjusted top-CCA | held-out top-CCA | raw ratio | **null-corrected ratio** | detection floor | attenuation | **R1, raw block** | **R1, residualised block** | perm *p* |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.000 | 0.5573 | 0.4932 | 1.000 | **1.000** | 0.20 | 1.130 | **196.2** | **210.2** | 0.0033 |
| 0.10 | 0.091 | 0.5571 | 0.5017 | 0.9996 | **0.999** | ≥ 0.40 | 0.985 | **194.1** | **210.2** | 0.0033 |
| 0.20 | 0.211 | 0.5447 | 0.5129 | 0.977 | **0.968** | ≥ 0.40 | 1.003 | **190.5** | **209.5** | 0.0033 |
| 0.30 | 0.302 | 0.5190 | 0.4986 | 0.931 | **0.905** | ≥ 0.40 | 1.057 | **187.5** | **208.4** | 0.0033 |
| 0.40 | 0.400 | 0.4774 | 0.4619 | 0.857 | **0.804** | ≥ 0.40 | 1.014 | **184.7** | **207.0** | 0.0033 |
| 0.60 | 0.600 | 0.3971 | 0.3680 | 0.713 | **0.607** | ≥ 0.40 | 0.855 | **176.5** | **205.1** | 0.0033 |
| 0.80 | 0.800 | 0.2844 | 0.1922 | 0.510 | **0.333** | ≥ 0.40 | 0.863 | **161.2** | **203.7** | 0.0033 |

*Provenance: `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2, §6;
`NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`. Cohort 6,427 patients, 238,610 tumour
patches, 7,644 slides; 2,766 evaluated on `test`. Instrument: 108-column cancer + pooled-TSS design, seed
42, 16 components, 20 draws, 300 permutations (resolution 1/301 = 0.0033). Permutation null median
0.145–0.147 at every level. Outputs under `p1_evidence/dilution/` **on the box** —
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/dilution/`, **not a path
in this repository**; the artifact this paper reads is
`~/p1_out/dilution/dilution_foreign_tumour_pca256.npz` and the build scripts are
`v2/research/dilution/`. **Statistic R1 (canonical), given for both
blocks.** Every value in both rank columns recomputed under the canonical implementation on 2026-08-04;
the raw column reproduces the originally published figures exactly (196.187, 194.102, 190.532, 187.494,
184.680, 176.523, 161.226).*

**The correction this paper must carry.** On the **raw** block — where the originally published numbers
came from — effective rank falls **196.2 → 161.2, i.e. −17.8%**, while the null-corrected channel falls
**1.000 → 0.333, i.e. −66.7%**. But **the channel is measured on the confound-residualised block, and the
two are different representations** (§3.1). Measured on the same block the channel is read from, the same
canonical statistic falls only **210.2 → 203.7, i.e. −3.10%**, and the miscalibration is **21.5×**, not
the 3.7× previously published.

| block | statistic | 0.00 → 0.80 | rank lost | miscalibration against the channel's −66.7% |
|---|---|---|---:|---:|
| raw | **R1 (canonical, as published)** | 196.187 → 161.226 | 17.82% | **3.74×** |
| raw | R2 | 147.039 → 96.854 | 34.13% | 1.95× |
| raw | R3 | 150.493 → 105.215 | 30.09% | 2.22× |
| raw | R1 uncentred | 193.752 → 155.116 | 19.94% | 3.34× |
| **residualised** | **R1 (canonical, matched to the channel)** | **210.179 → 203.667** | **3.10%** | **21.53×** |
| residualised | R2 | 170.674 → 160.946 | 5.70% | 11.70× |
| residualised | R3 | 173.348 → 165.545 | 4.50% | 14.82× |

*Provenance: `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`
§5; `~/ws_rank/RANK_RECOMPUTE.json`.*

**The direction survives every implementation choice; the magnitude does not.** Rank under-reports the
information loss in **every** cell, by between **1.95× and 21.5×**. The previously quoted constant of
"3.7×" is one cell of that table, and its provenance — raw-block rank against residualised-block channel
— was never stated. **This paper therefore quotes the range, with 21.5× as the matched-preprocessing
figure and 3.74× as the raw-block sensitivity.**

**This correction moved the result in our own favour, and is reported for that reason.** Matched
preprocessing makes the under-reporting nearly six times larger than the number we had been quoting. A
correction that flatters the correcting party is exactly the kind that gets left undone, so it is stated
here rather than absorbed silently — and the more flattering single number is not used anywhere.

**What it adds and what it does not contradict.** It does **not** contradict RankMe: high rank with
degraded information is precisely the necessary-not-sufficient case RankMe reserves, and LiDAR's own text
says so — *"This illustrates that high rank is a necessary but not a sufficient condition for high
performance."* What is added is a **magnitude**, which a two-point comparison cannot show. Level by level
on the raw block, the ratio (fraction of rank retained) / (fraction of channel retained) reads **1.000,
0.990, 1.003, 1.056, 1.171, 1.482, 2.467** — essentially 1 while almost nothing is happening, rising
steeply exactly when the damage becomes real. (The dip to 0.990 at d = 0.091 is a single-seed wobble on a
level where the channel changed by 0.001; nothing is claimed about strict monotonicity.) On the
residualised block the divergence is present from the first level and larger throughout. Either way it
defeats the natural retreat position — "rank is at least a rough guide" — which neither the
necessary-condition framing nor a correlation coefficient addresses.

**What travels with it,** quoted from the source: the detection floor is **censored** — the grid tops out
at 0.40 and the floor reads 0.40 from d = 0.09 onward, so it is "≥ 0.40"; the transmission floor reads
0.05 everywhere, the finest level, so it is censored from below too. *"The whole curve is one
representation"* — a trained attention aggregator could plausibly down-weight foreign patches, so the
number is a property of unweighted mean pooling, not of the modality. And it is a **single seed (42) and
a single draw of donor assignments**, which "gives no error bar on the level-to-level differences";
monotonicity over seven levels is what carries the result in the absence of intervals. Finally, the
source file is named `DILUTION_LOWER_BOUND.md` and **its own §4 withdraws the phrase "lower bound"**: the
measured quantity is *"the cost of preparation-matched, information-free contamination"*. Nothing here
turns on which bound it is, since both quantities are affected identically by the contaminant's nature.

### 4.9 The historical instances, recomputed — what survives, what is qualified, what is withdrawn

Four instances predate this paper and were previously presented as a set of "dissociations". They are
retained because they are the provenance of the hypothesis and because withdrawing them silently would
misrepresent how the claim was arrived at. **They are not a count, and none of them carries §1.3's
claim.**

| instance | manipulation | statistic & block | rank | information | verdict now |
|---|---|---|---|---|---|
| **D2** (§4.6) | supervision target table: Hallmark vs perturbation basis; one method, one architecture, 3 seeds, matched by construction | R1, residualised | 23.39 / 28.77 / 9.14 vs 14.87 / **34.12** / 9.11 | Δ channel −0.1325 / −0.1089 / −0.1226; both CIs exclude zero 3/3 | **Survives, but restated** — see below |
| **dilution** (§4.8) | patch-bag contamination, d = 0 → 0.80, 7 levels, zero fitted parameters | R1, **both blocks** | 196.2 → 161.2 (−17.8%) raw / 210.2 → 203.7 (**−3.10%**) resid. | null-corrected channel 1.000 → 0.333 (−66.7%) | **Survives with a corrected magnitude**; does not contradict RankMe |
| **Phase 1b** | objective profile `full` → `programme_only`, single seed, arms **not verified matched** | R1, raw | 38.4834 → 32.0594 (−16.7%) | held-out top-CCA 0.4768 → 0.4748 (−0.002) | **Neither quantity is resolvable at one seed** — see below |
| **"16/16"** | full schedule vs contrastive-only, 16-patient train batch | **hard `matrix_rank`** — none of R1/R2/R3 | pinned at **16/16** in every arm | patient cosine 0.7089 → 0.9999; retrieval 0.062 → **0.000**, below chance | **WITHDRAWN as a rank instance** |
| **decorrelation** | covariance-decorrelation term added | **unknown** | 49.9 → 103.3 (+107%) | "within-cancer specificity" 0.1366 → 0.1367 | **`[NOT RECOMPUTABLE — artifact never existed]`**; history only |

*Provenance: D2 — `D2_RESULT.md` §2, §4, outputs `~/e0_run/d2_v3/bootstrap/` and
`D2_PER_ARTIFACT_READOUT.json`. Dilution — as §4.8. Phase 1b — `PHASE1B_TARGETED_READOUT.md` §3, §5, §7;
readout `runs/calibra_v3_targeted`. "16/16" — `NOTEBOOK.md` 2026-08-02 01:20 UTC and
`NOTEBOOK_ENTRIES/g26_variance_floor_fix_20260803T0210Z.md`. Decorrelation —
`v2/research/rebase/ENGINE_CLD.md` §1 and `HANDOFF_BUILD_AGENT.md` §1–2 only. Recomputation:
`NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md` §5.*

**D2 seed 43, restated as implementation-dependent.** Earlier drafts called *"in seed 43 the higher-rank
arm loses (34.12 against 28.77)"* the paper's single strongest sentence. Under the canonical R1 it is
true and the arithmetic reproduces exactly. **It is false under R2**, which orders H above I (13.23
against 12.97), and under R3 it depends on the block (§4.5b). What survives all four combinations is
narrower and is what we now state:

> **The rank ordering is wrong in 1 of 3 seeds under every combination of statistic and block we
> computed — and it is never the same seed.**

What does *not* survive is any sentence of the form "in seed 43 the higher-rank arm loses" without
naming the statistic and the block. D2's other useful property is unchanged and does not depend on the
inversion: seed 44's arms are **equal in rank to two decimals** (9.1426 against 9.1052, a 0.4% gap that
§4.4 shows is 1.4 sampling SDs, i.e. no gap at all) while differing by **−0.1226** in channel with both
CIs excluding zero. Its rank column was recorded as audit check A5 of a predeclared post-training
checklist whose instruction reads *"reported, **not** interpreted"*, added for an unrelated worry, so it
was not collected to support any conclusion.

**Phase 1b, re-read honestly.** The originally reported reading was "the representation loses 17% of its
effective rank and its molecular channel is unchanged at the second decimal". Its `wsi_identity` rows are
a genuine internal control — that head is a frozen MLP-CLIP teacher passed through, `max|diff| = 2.6e-04`
between arms, with rank and channel identical to four significant figures, so the instrument does return
"nothing changed" when nothing changed. But given §3.5 — the same seed retrained moves a channel estimate
by 0.035 and a rank estimate by up to 3.3× (§4.1) — **a 0.002 channel difference and a 6.42-point rank difference are
both inside this stack's run-to-run noise. Instance 2 shows that neither quantity is resolvable at one
seed, not that one moved and the other did not.** Its own source adds that the arms were never verified
matched (G0.4) and that *"no CI on any between-run difference; a paired bootstrap on the biology gap is
still required"*, and that a neighbouring claim from the same run has already been withdrawn.

**The "16/16" instance is withdrawn, and it is evidence *for* the collapse-diagnostic use.** Arm A is
total representational collapse: every patient's `z_biology` converges to the same vector (off-diagonal
cosine 0.9999), cross-modal positive and negative pairs become indistinguishable (0.9959 against 0.9960 —
the negatives marginally *higher*), and retrieval falls to 0.000, **below** its chance level of 0.062;
extended to 5,000 steps a sibling arm reaches patient cosine 1.0000 with InfoNCE pinned at exactly
ln 16 = 2.7726 from step 2,500 onward. The rank column reads 16/16 in every arm including that one — and
**four reasons it must not be taken at face value**: (i) its source labels it "`z_biology` matrix rank",
a **hard numerical rank**, none of R1/R2/R3, and a 16 × 256 float matrix has full row rank under
essentially any perturbation; (ii) its maximum is 16 **because the batch is 16**, a structural ceiling set
by the experiment's design; (iii) **the centred effective rank of the same objective does fall — to
1.00**, recorded as `eff-rank 12.88 → 1.00 by step 50` on the same gate, cohort and seed with positive and
worst-negative cosines converging to 0.9993; (iv) it is a train batch of 16, not held-out. This project
previously listed it among its two strongest instances and **that description is withdrawn**, here and in
P1. What it does establish is that **hard/numerical rank is worthless as a collapse diagnostic** — and it
is nonetheless what a `matrix_rank` call returns, and `d1_geometry_probe.py:53` computes and prints one
beside the effective rank, so this is a live confusion in our own tooling.

**The decorrelation instance is history, not evidence.** `HANDOFF_BUILD_AGENT.md:98` cites
`paper/.../RESULTS.md`; **no such file exists in this repository**. "Within-cancer specificity" is defined
in no file now present and is none of the statistics in §3.2. The rank statistic predates `spectral.py`'s
consolidation and cannot be assigned to R1, R2 or R3. It is `[NOT RECOMPUTABLE]`, it is excluded from
every count and summary statement in this paper, and it belongs in a history paragraph rather than a
results table. One thing it does contribute, from a later and better-sourced measurement: the same
decorrelation term was independently characterised on a live objective and **its global minimum is a
collapsed state** (38.97 on a healthy batch against 1.19e-17 on an all-identical one); it self-extinguishes
at every weight from 0.001 to 4.0 within 25 steps (20.74 → 0.00); and a per-dimension variance floor
provably cannot stop it, because the rank-1 family `zᵢ = m + aᵢ·u` satisfies one. **A regulariser
introduced to raise rank whose own minimum is collapse is a sharper cautionary tale than the +107% number
ever was** — and we state it about **our implementation**, not about VICReg or Barlow Twins, which we have
not reimplemented faithfully or benchmarked.

**And every one of those sentences carries a condition we did not know we were assuming: *in the
absence of a momentum key encoder*.** Every measurement above — the self-extinction, the five-arm
sweep of §5.1's instance 3, the "decorrelation aggravates the collapse" reading — was taken on a queue
written by the **query encoder**. Re-run at `m = 0.999`, the same term does the opposite (§4.9a). No
claim about `feature_decorrelation` in this paper, in `paper/QUEUE_ANCHORING.md` or in
`paper/LIVENESS_GATE_DESIGN.md` may be quoted without that qualification attached.

**One further instance is `[NOT RECOMPUTED]` rather than recomputable.** The D1-A geometry probe's
`programme_only` 9.81 / 10.47 against `programme_free` 1.71 (epoch 39, 282 held-out patients) is **R3 on
live checkpoints**; recomputing it under R1 needs a GPU forward pass from the surviving checkpoints
(`~/e0_run/d1_v1/d1_p_seed{42,43,44}/last.pt` and `~/e0_run/d1_v1/d1_f_seed42/last.pt`; the
`d1_f_seed43` and `d1_f_seed44` directories exist but hold no checkpoint — the gate refused those arms) and the GPU was occupied by a training run that was not contended for. It is quoted nowhere in this
paper as a rank instance, and its source entry forbids the reading in any case: *"Nothing about programme
supervision may be concluded from it — the contrastive arm never trained, so the comparison measures a
defect, not an ablation."*

### 4.9a The decorrelation term's defect was conditional on the queue — and the same three runs dissociate rank from a co-measured collapse measure

Two findings, and the second is the stronger dissociation in this paper.

**(i) "`feature_decorrelation` is defective" was never unconditional.** It was conditional on a queue
written by the **query encoder** (§5.2). Re-run with a momentum key encoder at `m = 0.999`, on a fixed
256-patient held-out probe, 400 steps, one verified common initialisation (R3 67.55, canonical R1
101.38, RNA-view mutual cosine 0.3650 at step 0 in all three arms):

| `decorrelation` | **R3** | **canonical R1** | **RNA-view mutual cosine** |
|---|---:|---:|---:|
| 0.0 | 4.32 | 6.29 | 0.4774 |
| 0.01 | 6.22 | 9.32 | 0.7657 |
| **0.04** | **8.01** | **12.20** | **0.8696** |

*Provenance: `~/e0_run/d1_diag/ablate_decorr{0.0,0.01,0.04}.log`, produced by
`v2/research/rebase/d1_momentum_probe.py`, which imports both rank statistics from `v2/calibra` and
computes neither inline; vendored and hashed at
`v2/research/rebase/p2/figures/data/e0_run/d1_diag/`
(`v2/research/rebase/p2/figures/data/MANIFEST.json`). Values at step 400.
Reported in `NOTEBOOK_ENTRIES/lr_test_and_decorrelation_reversal_20260804T1130Z.md` §2. **One seed per
level.** Note that the "eff-rank" a reader would take from these logs' own `final_eff_rank=` line is
the **R3** column; the canonical column is given beside it because the floor of §4.1 is R1.*

Without momentum the same term *aggravated* the collapse — 1.59 against 2.17 at step 250 with
`m = 0` (§5.1, instance 3). With momentum it raises rank monotonically. **Every claim this project
has made about `feature_decorrelation` therefore needs *"in the absence of a momentum key encoder"*
attached**, and §2.4, §4.9, §5.1 and Appendix C now carry it.

**(ii) Rank rises while a direct collapse measure rises with it, monotonically, on the same runs.** As
decorrelation increases, effective rank goes **up** (4.32 → 8.01 under R3, 6.29 → 12.20 under
canonical R1) and the RNA-view patient-to-patient mutual cosine goes **up too** (0.4774 → 0.8696).
Rank says the representation is occupying more directions; a direct measurement of the condition rank
exists to detect says the patients' states are converging on one vector. At 0.8696 the RNA-view states
of different patients are nearly the same vector — the thing §4.10 defends rank for detecting — and
rank moves the wrong way as it worsens.

This is a stronger dissociation than anything in §4.9, for two reasons that are about the *shape* of
the evidence rather than its size. It is **monotone across three levels** rather than a single
contrast. And the contradicting quantity is **co-measured on the identical runs**, in the same log
lines, rather than inferred from a downstream readout in another table — which is the objection
§4.9's instances cannot answer.

**And the paper's own criterion applies to it, which is why the magnitude is not what carries it.**
The rank change is **1.854×** under R3 and **1.940×** under canonical R1, both **inside** §4.1's
3.295× floor (and on a third block again — the fixed held-out probe, which has no floor of its own;
§4.1a rows 48–49). **The monotonicity and the co-measured cosine carry this observation; the size of
the rank change does not, and no sentence here may be quoted as though it did.** One seed per level,
400 steps, one objective. **Figure F9.**

### 4.10 The use that survives, and where its boundary actually is

**Rank does detect total collapse.** Every measurement here of a representation collapsed to a single
direction shows a centred effective rank at or near 1–2:

| regime | statistic | value | co-measured evidence of collapse |
|---|---|---:|---|
| clean in-batch InfoNCE, step 50 | centred eff-rank (diagnostic script) | 12.88 → **1.00** | pos cos 0.9993, worst-neg cos 0.9993, min margin −0.0001 |
| `programme_free` at training scale, step 150 | R3 | 67.55 → **~2** | RNA-view mutual cosine 0.9813 |
| `programme_free`, epoch 21 / epoch 39, 282 held-out patients | R3 | **1.76** / **1.71** | RNA–RNA mutual cosine 0.977 / 0.986; hard rank 9 / 11 |

*Provenance: `NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md`;
`d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`;
`d1_programme_free_collapsing_in_training_20260803T1930Z.md`;
`d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`. **Statistic named per row**; the
first is a diagnostic-script centred effective rank and the rest are R3, so no two rows are compared.*

**But the boundary is far lower than "low rank".** At canonical R1 effective rank **9.11 and 9.14 of a
nominal 256** — 3.6% of ambient, a number any rank-monitoring practice would read as severe collapse —
the two D2 seed-44 arms still read held-out channels of **0.5983 and 0.4757** against a permutation null
of **0.140**. A representation at 3.6% of its nominal rank was carrying a large, permutation-significant
molecular channel.

So the defensible statement is not "low rank means little information". It is:

> Effective rank near its floor (≈ 1–2, with patient-to-patient mutual cosine ≈ 1) is reliable evidence
> of total collapse. Anywhere above that — including at 3.6% of nominal dimensionality — it is
> uninformative about the channel, in both directions, and non-monotone with respect to it.

**And even for the surviving use, rank is not the best instrument.** Patient-to-patient mutual cosine
reaches 0.977–0.999 in every collapsed arm above, is a single matrix product, needs no SVD, and has a
natural scale with a meaningful maximum. In every instance here where rank correctly signalled collapse,
mutual cosine signalled it too and more legibly. **As a partial remedy: if you report effective rank, at
least also report the patient-to-patient mutual cosine and the seed spread of both.**

**A second scalar fails too, and in the opposite direction**, so this is not a recommendation of
per-feature spread. `programme_free` at epoch 21 has **higher** mean per-feature standard deviation than
`programme_only` (0.0137 against 0.0044) and **lower** effective rank (1.76 against 7.38), because the
collapse is to the family `zᵢ = m + aᵢ·u` rather than to a point. At epoch 39: 0.0156 against 0.0056, and
1.71 against 9.81 / 10.47. Isotropic per-feature std for d = 256 is 0.0625.

**What we have not shown.** We have not found a case of total collapse that effective rank missed. The
one-sided claim is therefore asymmetric on the evidence as well as in its statement, and this section
should be read as *"we did not falsify the collapse-diagnostic use"*, not as *"we verified it"*.

---

## 5. A worked example: the metric used inside the one regime this paper says it works

The two sections that follow were drafted separately as `paper/LIVENESS_GATE_DESIGN.md` and
`paper/QUEUE_ANCHORING.md` and are integrated here by the author's decision. **This is a worked example
of the metric used in the regime the paper says it works, not a separate contribution**, and it is not
to be split out: it belongs here for a reason that is not merely convenience — **it is the same metric,
on the same stack, used throughout as a collapse diagnostic, the one use §4.10 defends — and it is
therefore the paper's demonstration of what the surviving recommendation looks like in practice.**

§5.3 states why the demonstration is not a contradiction. **§5.4 states the one place it fails the
paper's own test, and does not defend it away**: the seed-replicated momentum separation is 3.29×
against §4.1's 3.295× floor, so **the rank difference behind our own fix is not resolvable either** —
and §5.4 sets out what the fix does rest on, which is not a rank ratio.

### 5.1 What a liveness gate certifies, and four ways we learned it certifies less than we assumed

**Claim.** A liveness gate that certifies a model **can fit** a small, fixed, favourable problem does
**not** certify that the objective **will learn** at the scale and duration it is about to run at. The
two differ whenever the training regime contains a dynamic the gate's regime removes — and gates are
usually designed to remove exactly such dynamics, because that is what makes them fast and
deterministic.

Four independent instances from one objective over two days. Each was believed to be a completed fix at
the time; each was falsified by extending the gate's regime toward the run's.

**Instance 1 — the window was shorter than the failure.** `programme_free`'s biology head was observed
collapsing under a covariance penalty. A per-dimension variance floor was added and the collapse
disappeared: at 800 steps the contrastive term descended to 2.0875 with patient-to-patient cosine
0.7261. Extending the same run to 5,000 steps reversed the verdict — 2.4626 at step 1,000, then pinned
at ln 16 = 2.7726 from step 2,500 onward with patient cosine **1.0000**. The fix delayed the failure past
the observation window; nothing about the 800-step measurement was wrong, it was answering "has it
collapsed *yet*". *Provenance: `NOTEBOOK_ENTRIES/g26_variance_floor_fix_20260803T0210Z.md`,
`g26_stepbudget_sweep_20260803T0340Z.md`.*

**Instance 2 — the batch was smaller than the problem.** The repaired gate memorises one fixed
16-patient batch. With the objective's regularisers excluded it reaches contrastive 0.012–0.057,
retrieval 16/16, patient cosine 0.0597 and **effective rank 5.81** — a clean pass on a criterion of
≤ 0.10, on three seeds. The same objective at cohort scale, 3,118 streaming patients, collapses to
**effective rank ~1.8 by epoch 21**, against the supervised arm's 7.38 / 7.35 on two seeds, measured on
282 held-out patients. Memorising sixteen items is a capacity question; representing three thousand is a
learning question. *Provenance: `NOTEBOOK_ENTRIES/g26_passes_20260803T1100Z.md`,
`d1_programme_free_collapsing_in_training_20260803T1930Z.md`. Statistic **R3** throughout this section
unless stated.*

**Those two rank numbers are not a ratio, and the instance does not rest on one — stated here because
§4.1a's audit required it to be resolved rather than flagged again.** 5.81 is read on a **16-patient
train batch against a frozen 64-key queue**; ~1.8 is read on **282 held-out patients** with a live
4,096-key queue at cohort scale. Those are two different blocks and two different cohorts, and §3.1's
own rule forbids comparing rank across blocks — so the ≈3.2× a reader might form from them is not a
quantity this paper may quote, and §4.1a records it as exempt on that ground rather than as a
comparison that fails the floor. **What the instance rests on instead is two things, neither of them a
rank difference.** First, the gate's own **binary** verdict: contrastive 0.012–0.057 against a
predeclared ≤ 0.10 bar with retrieval 16/16, on three seeds — a pass, not a margin. Second, the
cohort-scale arm being **independently established as collapsed**, by quantities measured alongside
its rank: RNA-view patient-to-patient mutual cosine **0.977 / 0.986**, hard rank **9 / 11**, against
the supervised sibling at 7.38 / 7.35 on the same probe (§4.10). The claim — *a gate that certifies
memorisation of sixteen patients does not certify learning at three thousand* — is a statement about
which dynamics the gate's regime removes, and it would stand unchanged if neither rank number had been
recorded.

**Instance 3 — the gate removed the dynamic that causes the failure.** The subtlest of the four, because
the removal was deliberate, documented, and correct on its own terms. The objective uses a queue of
detached negative keys. Replaying one batch against a *live* queue makes the queue fill with
re-encodings of that same batch, so the negatives become the queries; the gate therefore **freezes** the
queue at 64 keys. That fix was correct and necessary — measured reduction 0.054 → 0.394, unique keys
16 → 64. But a frozen queue supplies a fixed reference frame and training's queue does not: its 4,096
keys are written by the query encoder itself, one step behind. Under a live queue the same objective
collapses to effective rank ~2 within 150 steps from an initialisation of 67.55, under **every** setting
of both regularisers suspected of causing it — including both at zero:

| `decorrelation` | `biology_full_consistency` | step 50 | 100 | 150 | 200 | 250 |
|---|---|---:|---:|---:|---:|---:|
| 0.04 | 1.0 | 4.08 | 1.95 | 2.16 | 1.68 | 1.59 |
| 0.0 | 1.0 | 2.62 | 2.16 | 2.47 | 1.94 | 2.17 |
| 0.04 | 0.1 | 2.99 | 3.43 | — | — | — |
| 0.0 | 0.1 | 2.97 | 2.00 | 2.50 | — | — |
| 0.0 | 0.0 | 2.98 | 1.98 | 1.86 | — | — |

*Provenance: `NOTEBOOK_ENTRIES/d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`; logs
`~/e0_run/d1_diag/`. **Centred R3 on a fixed 256-patient held-out probe**, one verified common
initialisation (67.55 at step 0).*

The pathology is therefore not in the objective's weighting; it is in the key set. **The gate froze the
queue to remove a *known* pathology, and in doing so removed the *unknown* one.**

**Instance 4 — the gate was read from a re-implementation of itself.** The sharpest of the four and the
cheapest to have avoided. To decide whether the objective was ready to launch, the gate function
(`_overfit_programme_free_contrastive`) was run from a standalone harness that reconstructed its inputs.
Three seeds passed with margin; the run was launched; the gate then failed *inside the runner*. For
identical seeds and an identical 2,400-step budget:

| seed | standalone harness | inside the runner |
|---|---:|---|
| 42 | 0.01871 | passed |
| 43 | 0.01206 | **0.50883** ✗ |
| 44 | 0.05666 | **2.14122** ✗ |

Three of three pass in the harness; one of three in the runner. Seed 44's in-runner value is close to
chance (ln 16 = 2.7726) — not a marginal miss. **The harness did not merely give optimistic numbers, it
inverted the verdict on two of three seeds.** Neither measurement is incorrect: the runner reaches the
gate having constructed the model from the full experiment configuration and having consumed a different
quantity of RNG, so `copy.deepcopy(model)` starts the memorisation loop from a different initialisation
and different dropout draws. The harness reproduced the gate's *code* and not the gate's *caller*. This
cost two launches: one aborted immediately, and one that lost two of three contrastive arms mid-experiment
after three arms had already trained to completion. *Provenance:
`NOTEBOOK_ENTRIES/d1_relaunch_20260803T1530Z.md`, `queue_size_implicates_the_key_set_20260803T2200Z.md`,
`d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`.*

**Design rule, in two halves — and the second half is the one usually missing.** *A liveness gate must be
read from the process that will perform the training, never from a harness that reconstructs the setup.*
**And: a gate that rejects runs a harness would have approved is not too strict — it is the component
working.** Every arm this gate refused belonged to an objective later measured, independently and at full
training duration, sitting at effective rank 1.71 against the supervised arm's 9.81; the gate's
strictness saved roughly twelve GPU-hours of training an objective that could not learn. The tempting
response to a gate that fails two of three seeds is to suspect the gate; here that would have been
exactly wrong. **The failure was where we read the gate, not how strict it was.**

**And the gate is itself not reproducible, which §4.1's argument applies to as well.** Eight runs with
identical inputs, identical seed and an identical 2,400-step budget gave `final_biology_contrastive`
values spanning **650×** (0.00859, 0.01076, 0.01770, 0.02019, 0.02407, 0.03266, 0.38009, 5.58511), a
**6/8 pass rate** against an unchanged 0.10 threshold, and a bimodal shape — six clustered within 4× and
two divergent (`NOTEBOOK_ENTRIES/g26_is_not_reproducible_20260804T0700Z.md`). At a 75% per-arm pass rate
the probability of all three contrastive arms clearing the gate is **0.42**. **Any table reporting arms
admitted by this gate must state which arms were admitted, because admission is a stochastic filter and
arms that fail it are not a random sample.** The replacement rank probe is far better behaved (§4.4(3):
4.7% spread, unimodal, 3.5× separation with an empty band from 1.98 to 6.92) but on *divergence rate*
specifically it currently rests on a thinner basis, not a better one — three repeats cannot rule out a
25% tail — and should be read as best-of-N rather than trusted for having looked tidy three times.

**What this implies for gate design.** (0) Read the gate from the process that will train, and do not
soften a gate that then fails. (1) State what the gate certifies, *in* the gate: "this model can memorise
16 patients against a frozen key set" is a different sentence from "this objective will learn", and only
the first is supported. (2) A gate's regime must match the run's in every dimension the run's failure can
live in; ours differed in three simultaneously — duration (800 steps vs 40 epochs), problem size (16 vs
3,118 patients) and key-set dynamics (frozen vs live) — and a fourth hid in the gap between the gate and
a faithful re-implementation of it. (3) Simplifications made to defeat one pathology are the first place
to look for the next. (4) A gate that cannot fail on the training pathology should not be quoted as
evidence about it.

### 5.2 A queued contrastive objective needs an independent reference frame, not fresher keys — found in rank, established on a binary outcome

**The failing configuration.** A patient-paired cross-modal InfoNCE with a queue of detached negative
keys, refreshed every step from the **query encoder itself** — the standard "end-to-end with a memory
bank" arrangement. At cohort scale (3,118 training patients, 4,096-key queue, 214-patient batches) the
representation collapses from effective rank 67.55 at initialisation to **~2 within 150 steps** and does
not recover. The collapse is invisible to the memorisation gate (§5.1, instance 3) and is **not** a
regularisation failure (the five-arm sweep above).

**The fix.** Give the queue its own **momentum key encoder**: an EMA copy of the query encoder,
`θ_k ← m·θ_k + (1−m)·θ_q`, kept out of the optimiser and used only to encode keys before enqueueing.
Measured with **capacity held at 4,096 in every arm**, so the number of negatives is identical and only
the key-writing encoder differs; all four arms start from the same verified initialisation. 40 epochs of
this objective is **583 steps**, so the table spans the real duration:

| step | m = 0 | m = 0.9 | m = 0.99 | m = 0.999 |
|---:|---:|---:|---:|---:|
| 0 | 67.55 | 67.55 | 67.55 | 67.55 |
| 50 | 4.10 | 3.88 | 8.65 | 9.35 |
| 100 | 1.62 | 3.51 | 6.49 | 7.03 |
| 150 | 1.62 | 2.15 | 4.56 | 6.99 |
| 200 | 2.26 | 1.65 | 5.70 | 7.60 |
| 300 | 3.32 | 2.70 | 6.01 | 7.33 |
| 400 | 2.18 | 2.31 | 5.50 | 7.84 |
| 500 | 2.43 | 2.34 | 5.50 | 7.61 |
| **600** | **2.81** | **2.23** | **5.88** | **7.42** |

*Provenance: `NOTEBOOK_ENTRIES/queue_size_implicates_the_key_set_20260803T2200Z.md`,
`momentum_rescues_rank_but_staleness_is_not_the_mechanism_20260803T2330Z.md`; logs `~/e0_run/d1_diag/`.
**Centred R3 on a fixed held-out probe.** One seed per momentum value.*

Three things this table shows that a single number would not. The effect is **monotone in m** and large —
2.6–3.3× at every step past 150. It is **durable**: both working arms are flat from step 200 to 600,
spanning the full training duration, which matters because two earlier "fixes" on this objective looked
correct inside a short window and failed outside it. And `m = 0.9` **fails**, tracking the no-momentum arm
rather than the working ones, so there is a threshold and it lies between 0.9 and 0.99. **`m = 0.999` is
used because it measured best in this sweep. No mechanism is claimed.**

**The sweep is one seed per momentum value, and that was a defect rather than a design choice.** The
momentum harness had its **seed hardcoded**, so the sweep could not have varied seeds had we asked it
to — the arms above are one seed because the apparatus admitted only one, not because a one-seed sweep
was judged sufficient. That is recorded for the same reason §5.1's instance 4 is: **a measurement
apparatus that silently constrains what can be measured produces something that reads as a decision and
is not one.** A seed-replicated sweep was armed against a disjunction predeclared before it ran (§5.4),
and it has since reported:

| statistic | m | seed 42 | seed 43 | seed 44 | within-arm spread |
|---|---|---:|---:|---:|---:|
| **canonical R1** | **0.999** | **11.26** | **10.45** | **10.55** | **1.08×** |
| **canonical R1** | **0** | **3.18** | **1.13** | **2.36** | **2.81×** |
| R3 | 0.999 | 7.40 | 6.85 | 7.15 | 1.08× |
| R3 | 0 | 2.81 | 1.05 | 2.06 | 2.68× |

*Provenance: `~/e0_run/d1_diag/mseed_m{0,0.999}_s{42,43,44}.log`, produced by
`v2/research/rebase/d1_momentum_probe.py`, 500 steps, three seeds per momentum. **Block: the fixed
held-out probe** — canonical R1 and R3 read there, **not** the residualised exported `wsi_biology` block
§4.1's floor is measured on, which is why §5.4 treats the two as non-substitutable. Reported in
`NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §3.*

**Every m = 0.999 seed exceeds every m = 0 seed, on both statistics.** The separation is therefore not an
artefact of the one seed the harness allowed, and the single-seed defect is closed. **What that does not
do is make the rank difference resolvable by this paper's own criterion** — the worst-case separation is
3.29× against a 3.295× floor — and §5.4 is about that, and about what the fix rests on instead.

**The explanation we expected, and why it is wrong.** The natural account is MoCo's — **He, Fan, Wu, Xie
& Girshick, "Momentum Contrast for Unsupervised Visual Representation Learning", CVPR 2020,
arXiv:1911.05722**, VERIFIED at full text. Three corrections that the standalone drafts required and that
are made here. **(i)** The repository previously contained *no* arXiv number, DOI or venue for MoCo
anywhere, only a bare "He et al., 2020"; the identifier is now stated. **(ii)** MoCo advances the account
**twice as a hypothesis** — *"Our hypothesis is that good features can be learned by a large dictionary
… while the encoder for the dictionary keys is kept as consistent as possible despite its evolution"*
and *"We hypothesize that such failure is caused by the rapidly changing encoder that reduces the key
representations' consistency"* — never as an established mechanism, and a falsification must say so.
**(iii)** MoCo ties the argument specifically to queue use — *"a slowly evolving key encoder is a core to
making use of a queue"* — so a falsification is only in scope if it is in a queue setting. **Ours is**,
and that is stated rather than left to be inferred.

**Three independent measurements rule the staleness account out.** (1) **The queue is never stale**: with
214-patient batches into 4,096 slots it **turns over completely every 19 steps**, so no key is ever more
than nineteen steps old — an argument available by arithmetic before any experiment. (2) **Key-to-encoder
agreement does not predict rank**: cosine(stored key, fresh re-encoding by the current query encoder) at
step 100 is 0.427 / **0.908** / 0.441 for no-momentum / m = 0.99 / m = 0.999 against effective ranks
2.58 / 6.65 / **6.89** — the best-agreeing arm does not have the best rank. (3) **The healthiest queue
holds the freshest keys**: sweeping capacity at fixed key encoder, the strongest arm is capacity 64,
entirely overwritten every single step (effective rank 6.17 against 2.16 at capacity 4,096).

**What we think is happening instead, stated as unconfirmed.** When the query encoder also writes the
keys, a transformation applied to **all patients at once** moves queries and keys *identically*; the
similarity structure between them is preserved, so the loss does not penalise it. **Collapse is free.** A
decoupled key encoder holds a slowly-moving reference frame, so a global transformation of the queries
now costs loss. On this account the queue's defect is that it supplies **no independent frame**, not that
its contents are old — an *anchoring* story rather than a *staleness* story. It is consistent with all
three measurements and with the memorisation gate passing on a frozen queue (a perfect anchor), and it
explains why no loss-weight setting helped. **It is not confirmed, and one sharpening of it was
predeclared, tested and refuted:** that the EMA time constant `τ = 1/(1−m)` must exceed queue turnover
`T = capacity/batch`. Five predictions were committed before the runs; **four were wrong and the
discriminating one inverted** — at m = 0.95 the arm with the *lowest* `τ/T` (0.52, capacity 8,192) was
the healthiest of its group at 3.67, where the criterion required it to be the worst. What the data
supports is narrower: rank is monotone in `m` and nearly flat in capacity, with a threshold between
`τ = 20` (fails) and `τ = 100` (works) that appears in the same place at every capacity tested. We record
that as an observation and explicitly do **not** advance it as a mechanism; it would be the fourth
proposed explanation, resting on the data that refuted the third.

**Why this may matter beyond one codebase.** The failing arrangement — a queue refreshed from the query
encoder each step — is a common simplification of MoCo, adopted precisely because it avoids maintaining a
second encoder. If the reason MoCo's second encoder is necessary is *anchoring* rather than *staleness*,
the simplification is unsafe **even when the queue is short enough that staleness is impossible**, which
is exactly the regime where practitioners reason it should be harmless. The turnover calculation is the
cheap diagnostic: it tells you staleness cannot be your problem, and on the anchoring account that is no
reassurance at all.

### 5.3 Why §5 is a demonstration and not a contradiction

Every rank number in §5.1 and §5.2 is used in the regime §4.10 defends: **as a collapse diagnostic, near
the floor, with independently co-measured evidence of collapse.** The readings that carry decisions are
`~1.6–2.9` (collapsed) against `5.9–9.4` (not collapsed), and in every case a second, cheaper measurement
agrees — patient-to-patient mutual cosine 0.977–1.0000, retrieval at or below chance, InfoNCE pinned at
`ln 16`. The gate's own bar (R3 ≥ 4.0) sits **inside an empty band** between the two distributions
(1.98 to 6.92), not near either. Nothing in §5 selects between two healthy configurations on a rank
difference, which is the practice §4.1 disqualifies.

Two features of §5.2 also distinguish it from the practice §4.1 disqualifies, and they are stated here
rather than in §5.4 because they are about the *shape* of the evidence and not about whether the ratio
clears a floor. **It is not a two-point comparison**: the effect is monotone across four values of m and
flat from step 200 to 600 in both working arms — nine time points × four arms, not one ratio. **And the
readings sit at the collapse floor, where §4.10 says rank is reliable**, with the failing arms
independently confirmed collapsed by mutual cosine; §4.1 constrains rank as a selection signal *between
healthy configurations*, and §5.2's failing arms are not healthy. **Neither of those observations makes
the ratio resolvable, and §5.4 does not let them be read as though they did.**

### 5.4 The paper's own standard, applied to the paper: our own fix's rank difference is inside our own floor

**The bad news first, and it is ours.** Every rank comparison behind the §5.2 fix is inside §4.1's
measured retraining floor. All three quantifications of the momentum effect, with the statistic and
block on which each was read:

| comparison | statistic, block | ratio | judged against | resolvable? |
|---|---|---:|---|:---:|
| m = 0.999 vs m = 0, **one seed**, step 600 (§5.2's original table) | R3, fixed held-out probe | **2.64×** (7.42 / 2.81) | 3.295× | **no** |
| m = 0.999 vs m = 0, **worst case over three seeds**, 500 steps | **canonical R1**, fixed held-out probe | **3.29×** (10.45 / 3.18) | **3.295×** | **no** |
| the same replication under the tripwire statistic | R3, fixed held-out probe | **2.44×** (6.85 / 2.81) | 3.295× | **no** |

*Provenance: `NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §3; logs
`~/e0_run/d1_diag/mseed_*`, `~/e0_run/d1_diag/`. The floor is §4.1's, **canonical R1 on the residualised
exported `wsi_biology` block**: 3.295× residualised, 3.111× raw.*

So the seed replication resolved the disjunction predeclared below in favour of **separation** — every m = 0.999 seed
above every m = 0 seed, on both statistics — and **did not** lift the comparison clear of the floor. By
§4.1's own criterion, **the momentum fix's rank difference is not resolvable either.** Stated the way we
would state it about someone else's paper: a 3.29× difference measured with an instrument whose
same-seed retraining spread is 3.295× is not a difference that instrument can see.

**Three qualifications are real, and none of them rescues it.** The floor was measured on
`programme_only`, at 40 epochs, on the exported residualised block; the momentum runs are
`programme_free`, at 500 steps, on a fixed held-out probe. **Different arm, different duration,
different block**, so the two quantities are not substitutable and this is *indicative* rather than a
like-for-like disqualification; **no like-for-like floor for this regime has been measured — though one same-seed PAIR now exists in it, and it is concordant.** §4.1a row 50 records it: `ablate_decorr0.04` and `mseed_m0.999_s42` are the same momentum, decorrelation, capacity, learning rate and seed, from the same step-0 state, and they span 1.066× at step 400 and at most 1.128× over the eight shared logged steps. **That is n = 2 and it may not be quoted as a floor**: §4.1's own argument is that any pair drawn from four concordant repeats spans at most 1.028×, so a floor measured that way would license everything in this paper. We record it because "no like-for-like measurement exists" is no longer the accurate sentence, and because the accurate one — *a pair exists, it is concordant, and a pair is not a floor* — is the one a referee will check. But a
qualification only rescues a number if it points the right way, and these do not. The arm mismatch runs
against us: §4.3 measures the floor as a property of the arm, and `programme_free` is the *unstable*
one — 6.05× across five seeds at step 200 against `programme_only`'s 1.18× — so a floor measured in this
regime would be expected to be **larger**, not smaller. The axis the replication varies is the **seed**,
which §4.2 measures as the dominant nuisance term for exactly this statistic on exactly this stack. And
the duration mismatch cuts the same way, since §4.1's floor accumulates over 40 epochs while these arms
have had 500 steps to diverge. **The honest position is that the comparison is not decisive and that
everything we know about its direction says a decisive version would be worse for us**, so we do not
treat the mismatch as a reprieve.

**The one rank measurement on this fix that does clear the floor holds fixed the variable that
matters, and we report it with that stated.** §4.4(3)'s controlled repeat — three repeats of a 200-step
probe with identical inputs, live queue — gives 7.15 / 6.92 / 7.25 at m = 0.999 against 1.80 / 1.46 /
1.98 at m = 0, an empty band from 1.98 to 6.92 and the **3.5×** separation §4.4(3) quotes (6.92 / 1.98;
R3, fixed held-out probe). That is outside the 3.295× floor, but by about **6%**, a margin no more comfortable than the
1.5% §4.1 declines to lean on for D1-B seed 43 — and, decisively, **those repeats hold the seed fixed**,
where §4.2 measures the seed as the dominant term and §4.3 measures it at 6.05× on this very arm. The
seed-varied replication is the version of that check that matters, and it reads 3.29× against 3.295×.
**A fixed-seed check that passes and a seed-varied check that fails are not two results to choose
between; the second is the one the rule in §1.3 asks for.**

**Block-matching is load-bearing here, which is why every rank in this paper is quoted with its
statistic and its block.** §4.1's own table makes the point: D1-B's arm ratios exist in **raw**
(2.02× / 3.09× / 1.68×) and **residualised** (2.190× / 3.246× / 1.738×) form, and seed 43's residualised
3.246× judged against the **raw** floor of 3.111× would read as *outside* the floor when on its own block
it is inside. A ratio and a floor from different blocks do not compare. The momentum numbers are on a
third block again — a fixed held-out probe with no floor of its own — and that is a limitation of the
comparison, not a defence of the fix.

#### The fix was never justified by a rank ratio, and this is what it does rest on

**Rank is how we noticed the problem, not how we established the repair.** What rank did was make an
invisible failure visible: the objective falls from 67.55 at initialisation to ~2 within 150 steps under
every setting of both regularisers suspected of causing it, while the memorisation gate — which freezes
the queue — passes it (§5.1, instance 3). That is the collapse-diagnostic use §4.10 defends, at readings
of ~2 where §4.10 says the diagnostic is reliable, and it is not a selection between two healthy
configurations.

**What established the repair is a binary outcome with a channel behind it.** Before the fix,
`programme_free` had never produced a completed, uncollapsed, exportable run in this project's history:

| | before the fix (D1-A, `~/e0_run/d1_v1/`) | after the fix (D1-B, `~/e0_run/d1_v2/`, `--biology-key-momentum 0.999`) |
|---|---|---|
| `programme_free` arms reaching epoch 40 | **1 of 3 seeds**; seeds 43 and 44 were refused by the liveness gate at contrastive **0.50883** and **2.14122** | **3 of 3** |
| state of the arm that did reach epoch 40 | **collapsed** — R3 **1.71** at epoch 39 on 282 held-out patients against `programme_only`'s **9.81 / 10.47**, RNA-view mutual cosine **0.986**, hard rank 11 | **not collapsed** — canonical R1 **13.418 / 7.600 / 6.394** on the residualised held-out block |
| exported artifacts, CALIBRA, paired bootstrap | **none.** `run_d1` raises on the first non-zero return code, so no exports, no CALIBRA and no bootstrap were produced at all | all six runs complete and live (`~/e0_run/d1_audit.log`, `[PASS] A1_all_six_runs_complete`), and **three paired bootstraps with CI₉₅ per seed** (§4.7.2) |
| held-out molecular channel | **not measurable** | **0.5412 / 0.5336 / 0.5126**, clearly above its own `random_control` at **0.4504 / 0.4738 / 0.4425** |

*Provenance: `NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`
(D1-A dispositions, epoch-40 geometry probe, gate outcomes);
`NOTEBOOK_ENTRIES/d1b_blocked_gate_does_not_exercise_the_fix_20260804T0500Z.md` (the momentum flag
verified in all three runners' argv and in the pair manifest); §3.3, §4.7.2 and §4.7.3 for the D1-B
values. **The pre-fix ranks are R3 and the post-fix ranks are canonical R1**, so the two columns'
rank rows are not a ratio and are not quoted as one — they are a collapsed/not-collapsed verdict, each
with independent co-measured evidence.*

**That is a change of kind, not a change of degree.** It is the difference between an experiment that
produces no readout and one that produces three, each with an interval; between a representation at
mutual cosine 0.986 across 282 patients and one carrying a channel above its own negative control.
**None of it is a ratio, and none of it is inside anyone's noise floor.** §4.7.3 states the same thing
from the other side, when it rejects the deflation that D1's necessity test merely re-measured the
collapse: *"`programme_free` is **not** collapsed this time… This is the 'both arms train' branch. The
queue fix (§5.2) worked."*

**Three limits on that argument, stated because they are the ones a referee will find.**

1. **It is a before/after across two launches, not a controlled contrast.**
   `~/e0_run/d1_v2/D1_PAIR_MANIFEST.json`'s
   `"objective_only_difference": true` matches the two *arms within* D1-B (§3.4); no manifest asserts
   that the momentum encoder is the only difference between D1-A and D1-B. The direction is unambiguous
   and the outcome is binary, but this is not the same object as a matched ablation.
2. **What the binary outcome supports is momentum versus none, not `m = 0.999` over `m = 0.99`.** At
   step 600 those two read 7.42 against 5.88 — a **1.26×** difference, deep inside every floor in this
   paper, with no seed replication of the intermediate values. **The specific value this project runs was
   selected by a rank comparison our own rule disqualifies**, and we say so rather than letting the
   binary outcome launder it. `m = 0.999` is retained because it measured best in a sweep whose
   resolution we do not trust and because the objective trains under it; that is the whole of the
   justification for the value.
3. **The seed replication is 500 steps, not 40 epochs.** Durability of the *replicated* arms past 500
   steps is unmeasured; the 40-epoch evidence is D1-B's completion, which is the binary outcome and not
   the rank curve.

#### What this costs the paper, and why it is stated at this length

**§4.1's rule has now disqualified every rank comparison this project has ever made** — the seven
between-arm differences of §4.1, and the hyperparameter choice of §5.2. We regard that as the rule
working rather than the rule failing, and the alternative — applying the floor to RankMe's use and
exempting our own — is precisely the double standard this paper accuses the practice of. **A rank ratio
inside our own measured floor cannot carry a claim, and we do not permit it to carry ours.**

What it does not cost is §5's purpose. §5 is a worked example of the metric used in the regime §4.10
defends, and the momentum episode is a *better* example for having failed the test: rank found a
collapse that a gate could not see, at readings where the collapse diagnostic is reliable, and then
proved unable to grade the repair it had prompted. **Diagnosis and selection are different uses, and
this section is the sharpest illustration in the paper that a metric can be good at the first and
useless at the second on the same objective, in the same week, on the same stack.**

**The disjunction we predeclared, and the way it under-specified its own outcome.** Before the
replication ran, this section committed to a reading of each branch: *if the m = 0 and m = 0.999
distributions separate across seeds, §5.2's fix clears §4.1's own bar and the tension disappears; if
they overlap, the momentum choice is a rank comparison this paper's rule disqualifies, and §5.2 must be
rewritten to rest on the objective's downstream behaviour rather than on its rank.* **The distributions
separated and the ratio still did not clear the bar.** The disjunction conflated two different
questions — whether the arms separate, and whether the separation exceeds the floor — and it was written
before the floor's value was known. We record it as under-specified rather than claim it resolved
cleanly, and we have taken the branch it assigned to overlap: **§5.2 now rests on the objective's
downstream behaviour, and the rank curve is reported as the diagnostic that found the problem.**

**Two further honest limits carried from the source drafts.** Instance 3's *cause* is established only as
far as "the key set, not the objective's weighting"; the mechanism is open, MoCo's staleness account was
predicted, tested and falsified, and no replacement has been confirmed. The capacity effect remains
confounded — capacity changes both anchoring quality and negative count, and we have not separated them.
And all of §5 is one objective, one architecture, one cohort; we claim the failure mode is real and cheap
to test for, not a general rate at which liveness gates mislead.

---

## 6. Limitations

Structured as named objections with answers, after Eklund et al. (2016), rather than as a list of
caveats.

### 6.1 "This is not novel"

**Largely correct, and the claim is stated at its narrowest.** The core negative — that effective rank is
an unreliable selection signal — **is already published**, including *within-method* (Aldeneh et al.,
ICASSP 2025), including as a named selection-rule failure (Otero, Mateus & Balestriero), including by the
authors of the leading replacement (LiDAR, ICLR 2024), and including a published low-rank/high-information
instance. **The necessity-violation framing this paper was previously built on is both pre-empted and, on
our own best experiment, unsupported (§4.7).**

What remains: (i) the **reproducibility floor** and the demonstration that essentially every between-arm
difference we measured is inside it (§4.1); (ii) the **arm/seed variance decomposition** (§4.2) and the
finding that the floor belongs to the arm (§4.3); (iii) the corrected **dose–response magnitude** (§4.8);
(iv) the domain. We have not found (i)–(iii) reported. The census that supports that statement is
abstract-level over 453 works and is a **lower bound on the prior art** (§2.2); if a paper reporting a
seed-level reproducibility floor on effective rank exists, contribution (i) collapses and this paper
should be withdrawn to a replication.

### 6.2 "You did not measure the thing that would settle it"

| would-be measurement | why it is absent |
|---|---|
| ~~A controlled repeat design for §4.1~~ | **CLOSED, at N = 5.** Five identical `programme_only` retrains at seed 42 give a **3.295×** rank floor against a **1.055×** channel spread, bimodally distributed. It replaced an n = 1 estimate of 2.69×. What it still cannot do is attribute retraining variance to the metric rather than to this stack, because there is only one stack; and it is measured on the *stable* arm at a fixed seed, so it is a floor rather than an envelope. §4.2 remains the contribution that does not depend on it. |
| ~~A seed replication of §5.2's momentum sweep~~ | **CLOSED, and it does not clear §4.1's bar.** Three seeds per momentum at 500 steps: every m = 0.999 seed exceeds every m = 0 seed (canonical R1 11.26 / 10.45 / 10.55 against 3.18 / 1.13 / 2.36), so the single-seed defect is closed and §5.3's disjunction resolves in favour of separation. But the worst-case ratio is **3.29× against §4.1's 3.295× floor**, so the fix's *rank* difference is not resolvable by this paper's own criterion. Different arm, duration and block, so indicative rather than a like-for-like disqualification — and all three mismatches point toward a larger floor in this regime, not a smaller one. **§5 has been rewritten around it: §5.4 states the failure against our own criterion in its own subsection and rests the fix on a binary training outcome (`programme_free` completing uncollapsed 0 of 3 seeds before, 3 of 3 after, with intervals) rather than on the rank ratio.** What is still missing is a like-for-like **floor** for the `programme_free` / 500-step / held-out-probe regime. One same-seed **pair** in that regime does exist and is concordant (1.066× at step 400, ≤ 1.128× over the shared steps; §4.1a row 50) — **n = 2 is not a floor** by §4.1's own argument, and it is recorded rather than used. |
| ~~A reproducibility floor for any statistic other than canonical R1, or on any view other than `wsi_biology`~~ | **CLOSED for the statistics and the views, at n = 5, from exports that already existed.** §4.1a now carries a floor for every statistic T1 scores (R1, R2, R3, PR, PR_rownorm, RankMe raw and residualised, stable rank, α-ReQ, LiDAR, hard numerical rank) on the exported block raw and residualised, and for canonical R1 on all three co-trained views — the same five same-seed repeats, re-read. It is a re-derivation, not a re-run, and it should have been done when the repeats landed. Three results came out of it and two of them cost us: the floor spans **1.000× to 3.295×** between statistics on one block, so judging one statistic against another's floor was never admissible; **RankMe as published is the more reproducible statistic on our own artifacts** (1.811× raw against our 3.111×); and the floor is a property of the **view** (1.019× on `rna_biology` against 3.295× on `wsi_biology`), so the divergent repeat's collapse is a property of its WSI encoder rather than of the run. What is *not* closed is the arm and the seed: every floor is measured on `programme_only` at seed 42, and is a floor twice over for that reason. |
| **A reproducibility floor on the fixed held-out probe, the in-run training batch, the 16-patient gate batch or a live checkpoint** | **Not measured, not recoverable from the exports, and this is now the paper's largest single gap.** The five repeats were *exported*, not probed: `~/e0_run/d1_envelope/` holds `rep{n}.npz` and a per-run `~/e0_run/d1_envelope/rep{n}/train_metrics.jsonl`, and neither carries a probe forward pass, a training batch's activations or a gate batch. **11 of §4.1a's 25 selections — and every rank number in §5 — sit on the first of those four blocks**, so this paper's central criterion cannot reach them in either direction: they are not failing it, they are outside it. **What would close it is one specific run: five same-seed repeats of the `programme_free` / 500-step configuration with `d1_momentum_probe.py` attached, read at a fixed step.** Until that exists, §5.4's statement that the momentum fix's rank difference "is not resolvable" is stronger than the evidence supports — the honest form is that it is *unjudgeable*, and §5.4 rests on a binary training outcome either way. |
| **A per-block ground truth for the D1 arms** | **Not run.** §4.6a re-scores the *D2* arms on six target blocks and finds the selection verdict unstable on all twelve metric rows. The D1 arms were never scored on any block but the gene sets, so §4.6a's D1 column is held fixed in every row and is **not** evidence that the D1 half is block-stable. |
| **A labelled linear probe on every artifact** | Not run. It is the reference standard RankMe and LiDAR were validated against; ours is a held-out canonical correlation against unsupervised molecular targets (§3.2), which is a different standard. |
| ~~The D1 paired bootstrap~~ | **CLOSED.** It existed all along and was hidden by the audit chain's stale absolute path. §4.7.2 now carries both estimators: decisive 3/3 on the patient bootstrap, **2/3 on the cancer-cluster bootstrap** with seed 43 at +0.0006. The stale path is still unfixed in the chain and should be. |
| An error bar on any dilution rank or channel value | Single seed, single donor draw; bootstrapping donor assignments is CPU-only and unrun. |
| An equivalence test on Phase 1b's channel difference | The paired bootstrap its own source says "is still required" was never run. §4.9 states that "unchanged" means the point estimates differ by 0.002 and nothing more. |
| Rank and channel under **one** statistic across **all** instances | **Done for D2, dilution, Phase 1b and D1-B** — one implementation, every surviving artifact recomputed, published values reproduced exactly. **Impossible** for the decorrelation instance (artifacts never existed; the cited source file is absent) and the "16/16" instance (a train-time batch of 16, never exported, and a hard `matrix_rank`). D1-A's 9.81 / 1.71 are `[NOT RECOMPUTED]` — a GPU forward pass from surviving checkpoints, and the GPU was in use. |
| Any instance on a **second architecture, cohort or modality pair** | All of it is one architecture family (transformer aggregator over frozen H-Optimus-0 patch tokens with a biology head), one cohort (TCGA), one modality pair (morphology → bulk expression). `claim_guards.no_external_cohort` is undischarged for every morphology result on this project. |
| **E1**, the preregistered rank-versus-information experiment in this repository | Built (`v2/calibra/e1_rank_information.py`, `aggregate_e1.py`, equivalence margin 0.10, three-seed requirement) and **never run**. It is the experiment this paper should have been built on. |
| Rank at **capacity scale**, where Deng et al. report a power law | Not measured. §2.3 argues their sweep is over capacity-like variables and ours is not; **we do not claim rank fails at capacity scale.** |
| A case where the collapse diagnostic **fails** | We have not found one (§4.10). |
| ~~The competing-metric and variance-decomposition scripts, vendored~~ | **CLOSED, and enforcing the rule found an error.** All five now live at `v2/research/rebase/p2/` with an end-to-end test (`v2/tests/test_p2_analysis_scripts.py`), vendored at commit `7b37dce`. Re-run from a byte-verified checkout (402/402 files by git blob SHA-1), **§4.2, §4.4(1), §4.5(b), §4.5(c), §4.6 and §4.7 reproduce to every published digit**. §4.5(a)'s second and third rows did **not**: they were a mislabelled statistic, are relabelled PR / PR_rownorm, and the count falls from 3 of 6 to 2 of 6 (§4.5a). The rule was worth enforcing precisely because it caught something. |
| Rank computed with RankMe's exact ε convention as a *published-comparable* number | Not attempted; §2.1's discrepancy means no number here is comparable to a published RankMe value. (A faithful RankMe *was* computed on our artifacts for §4.2 and §4.6 — that is an internal comparison, not a cross-paper one.) |

### 6.3 "Your ground truth is the weak link"

Conceded, and the citation that makes the objection sharpest is ours to supply: Zaiem et al.
(Interspeech 2023; *Computer Speech and Language*) report that *"altering the downstream architecture
structure leads to significant fluctuations in the performance ranking of the evaluated models"*. Our
reference standard is a 16-component held-out canonical correlation against confound-residualised
molecular targets, with a measured permutation null of 0.140 (D2 row-shuffle) or 0.145–0.147 (dilution
within-cancer) depending on the readout — §3.2's footnote, which also records that we got this wrong
once. Four things partially answer the
objection and none of them dismisses it: the ordering is **98.0% arm** with F = 128.2 (§4.2); it is
**identical across all three co-trained views for all six pairs** (§4.5c); it **survives held-out
re-estimation with equal or larger deltas** (§3.2); and the `random_control` arm difference spans zero
in all three seeds (§4.7.3).

**Two things are not answered, and the second is this paper's second-largest exposure after the
single-stack limitation.**

*First*, it is one readout, and a different downstream task could order the arms differently. **If the
ground truth moves, §4.6's counts move with it — but §4.1–§4.3 do not, because they are statements about
rank's own variance and require the ground truth only to establish that the information ordering is
stable, which it is under every perturbation we applied.**

*Second — and stated because §4.1 holds rank to this standard, so the paper must hold itself to it.*
**The margins we work with are not comfortably above the instrument's own floor.** The `random_control`
level is 0.4425–0.4810 against real targets' 0.51–0.62 and the D2 row-shuffle permutation null of
**0.140** (**not** the dilution sweep's 0.147 — §3.2's footnote), so random gene sets sit at ≈3.2× the
floor — consistent with T1.4's independent finding that covariate-matched random gene sets reproduce
76–82% of the per-target channel. The real-versus-random
margin is therefore about **0.139**, and the smallest arm difference we quote (0.0705) is **about half
of it** (§4.7.4). The negative control does what a negative control must — it does not manufacture an
arm *difference* — but its **absolute level is high and separately explained**, and no absolute channel
number anywhere in this paper may be read against an assumed null of zero. §4.7.4 sets out the three
respects in which our paired, interval-backed differences are nonetheless in a different position from
rank's uninterval-backed levels, and states plainly that none of the three makes our own margin large.

### 6.4 A defect in our own evidence base, reported not repaired

The decorrelation instance's cited source (`paper/.../RESULTS.md`) does not exist in this repository. We
have not attempted to reconstruct it. Its benchmark statistic, "within-cancer specificity", is defined in
no file now present, so we cannot state what 0.1366 measures, what its chance level was, or whether it
could resolve a change. The number is retained in §4.9 as history, is marked `[NOT RECOMPUTABLE]`, and is
excluded from every count and summary statement. **Any future citation of "+107% rank at flat benchmark"
from this project must carry this paragraph.**

### 6.5 "The claim is about a scalar, not about geometry"

Nothing here says representation geometry is uninformative. It says one scalar summary of the spectrum is
unusable *as a selection signal*, in one regime, on one stack. Alignment, uniformity, the full
singular-value profile, LiDAR's statistic, per-feature spread and mutual cosine may each behave
differently — and §4.10 recommends one of them over rank. Per-feature spread is *also* misleading in at
least one place, in the opposite direction (§4.10), which is reported as a second scalar failing rather
than as a recommendation.

### 6.6 Scope of the negative

This paper does not claim that anti-collapse regularisation is useless, that rank should never be
computed, that RankMe is wrong about what it claims, or that the published defences in §2.3 are wrong. On
our own best-matched three-seed experiment RankMe's necessity hedge **held** (§4.7). What is claimed is
narrower and is about usefulness: that on this stack, inside RankMe's own reserved same-method regime,
effective rank's between-arm differences are smaller than its own within-arm reproducibility floor; that
two-thirds of its variance across matched artifacts is training-seed nuisance while the information those
artifacts carry is 98% arm; that the floor is a property of the arm and so cannot be calibrated once and
reused; and that the between-arm verdict is not even well defined until three unstated implementation
choices are fixed.

Following the exemplar's phrasing: **one caveat to our results is that they are limited to one
cross-modal biomedical pipeline, one cohort and one architecture family. It is possible our findings are
due to idiosyncrasies of this stack's non-determinism and would not generalise.** The reason to report
them anyway is the same reason the metric is used everywhere: a scalar promoted as a general, label-free
criterion has to work in the places people will apply it, and showing that it does not resolve its own
re-measurement *anywhere* is informative about the practice of quoting it *everywhere*.

---

## 7. Conclusion

Effective rank is cheap, label-free and never fails to return a number. Those three properties are why it
is used to select between representations, and they are also why it is silently unreliable in that use.

**We begin with what went against us.** We ran a preregistered test designed to break RankMe's
necessary-but-not-sufficient hedge, on two arms of one method with three seeds, and the hedge held: the
higher-rank arm carried the larger held-out molecular channel in all three, interval-backed 3/3 under a
patient bootstrap and 2/3 under the conservative cancer-cluster one. Rank was right about which arm to
keep. That result is not a footnote to this paper; it removed the claim this paper used to make, and the
claim that replaced it had to be one that survives rank being right on average.

It is. **The metric's between-arm differences are smaller than its own within-arm reproducibility floor,
inside the regime its authors reserve for it.** All seven of the between-arm rank differences this project
ever measured — 1.004× to 3.246× — are inside the **3.295×** spread measured over five identical retrains
of one configuration at one seed: a floor whose distribution is **bimodal** (four repeats within 2%, one at
a third of them), and on whose same five runs the molecular channel spreads only **1.055×**. Across twelve
matched artifacts, 65.5% of the variance in effective rank is training-seed
nuisance and the arm term is not significant; across the same twelve, 98.0% of the variance in the
molecular channel is the arm, at F = 128.2. And the floor is not a constant that could be calibrated once
and reused: at a fixed training step, one arm spans 6.05× across seeds while its sibling spans
1.18×. The instability is in training, not in estimation — the statistic itself is precise to about 0.1
on a rank of 25 — which is what makes it fatal rather than fixable by measuring harder. So we propose the
check that would have caught it: **retrain one configuration, measure the rank spread, and require the
between-configuration difference to exceed it.**

We report at equal prominence the evidence that constrains this. Centred effective rank does fall to ~1
under total collapse, so the collapse-diagnostic use survives — and §5 is that use, at length, on a
training gate and a queue fix that this project actually made and that this paper grades in rank. **The
worked example ends by failing the paper's own test.** The queue fix's seed-replicated rank separation is
**3.29×** against the **3.295×** floor above, so our own hyperparameter decision is a rank comparison our
own rule disqualifies; §5.4 says so in its own subsection and rests that fix where it always actually
rested — on a binary training outcome, an objective that had never once completed training uncollapsed
and that does so on three of three seeds after the repair, with a measured channel and paired bootstrap
intervals where there had previously been no readout at all. Rank was how we noticed the failure; it was
never how we established the repair, and a ratio inside our own floor is not permitted to carry our claim
any more than it is permitted to carry anyone else's. A representation at
3.6% of its nominal dimensionality still carried a large, permutation-significant channel, so "low rank
means little information" is false outside total collapse. The instance this project once called its most
dramatic turns out to be a *hard* matrix rank at a structural ceiling, and that description is withdrawn
here and in P1. Three published alternatives were computed on the same artifacts and none is better — the
strongest of them picks the information-poorer arm 3/3 on our in-scope contrast at every value of its
unspecified constant — and every selection-rule count in this paper is reported with the statement that
six pairs cannot support such a count. Our own headline miscalibration constant moved from 3.74× to 21.5×
when we matched the preprocessing, in our favour, and the flattering number is not used.

**And the paper's clearest single piece of evidence is one we did not design and would rather not have
had.** While enforcing our own rule that every number trace to a file in the repository, we found a
*fourth* statistic living under the name `effective_rank` — in the analysis script behind the very
section that argues the name is unreliable. It moved a published count from 3 of 6 pairs to 2 of 6 and
removed a qualification from the result that had already gone against us. It was invisible to review,
invisible to the tests, and invisible to us across two drafts and a full recomputation; it surfaced
only when the code that produced the numbers was vendored and forced to call the one canonical
implementation. **We could not have argued the thesis better than that, and we did not argue it — we
tripped over it.** The practical form is the recommendation we would most like carried: the defence
against a scalar whose name is stable while its definition is not is *mechanical provenance*, not care.

Finally, we report that three mutually incompatible statistics were implemented under the name
`effective_rank`, across ten call sites, in the repository that produced these results — two of them
thresholds that abort training runs — and that our own historical instances were not all measured with the
same one. We found this while writing this paper. We then fixed it, and the fix is part of the result:
under one definition every surviving instance reproduces exactly, but the choice of statistic is worth up
to a factor of three, it reverses the arm ordering on a third of our pairs, a fourth choice nobody had been
stating — raw or confound-residualised block — moves our own headline magnitude by 5.8×, and simply
measuring a different co-trained view of the same model flips the verdict on two pairs in six. The
information verdict flips under none of them.

**It's critical to empirically verify that a cheap summary actually reflects a functionally relevant
property of the system being examined** — and the cheapest such verification, for any spectral summary, is
to run the thing twice.

---

## Appendix A — provenance index

| § | claim | source file(s) | box path |
|---|---|---|---|
| 2.1 | Roy & Vetterli definition; RankMe claims, restrictions and the ε placement | full-text PDFs, verified 2026-08-04 | — |
| 2.2 | Aldeneh et al. quote and identifiers | Crossref `10.1109/ICASSP49660.2025.10889651`; arXiv Atom `2409.10787` | — |
| 2.2 | LiDAR quotes, Spearman/Kendall values, ICLR 2024 venue | full-text PDF arXiv:2312.04000v1; DBLP `conf/iclr/Thilak0SDGNSL24` | — |
| 2.2 | the 453-work census and its limits | `NOTEBOOK_ENTRIES/p2_prior_art_citation_graph_sweep_20260803T2326Z.md` | — |
| 2.3 | Deng et al. and Awasthi et al. quotes and identifiers | arXiv Atom `2510.10948`; Crossref `10.1186/s12864-025-11913-2` | — |
| 2.5, 4.6 | competing metrics computed on our artifacts; LiDAR δ sweep; α-ReQ estimator | `NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §3, §6 | `~/e0_run/P2_METRICS_D2.json`, `P2_METRICS_D1.json`, `p2_competing_metrics.py`, `p2_selection_rule.py` |
| 3.1 | canonical definition; R1/R2/R3 as Hill numbers; the ten call sites; the raw/residualised choice | `v2/calibra/spectral.py` (`CANONICAL`, `RANK_VARIANTS`); `v2/tests/test_effective_rank_canonical.py`; `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md` §1–§4 | `~/ws_rank/RANK_RECOMPUTE.json` |
| 3.1 | E1 built, never run | `v2/calibra/e1_rank_information.py`, `aggregate_e1.py:38`; absence confirmed against `GATE_LOG.md` and `runs/` | — |
| 3.2 | channel statistic, nulls, held-out re-estimation | `v2/calibra/spectral.py:78-108`, `v2/calibra/run_calibra.py`; metrics entry §4.4 | — |
| 3.4 | D2 and D1 pair-manifest hashes and `objective_only_difference` | `D2_PAIR_MANIFEST.json`; `~/e0_run/d1_v2/D1_PAIR_MANIFEST.json` | `~/e0_run/d2_v3/`, `~/e0_run/d1_v2/` |
| 3.5, 4.1 | seed non-reproducibility; re-export vs retrain | `D2_RESULT.md` §4 | `~/e0_run/d2_v3/recovered_artifacts/` |
| 4.1 | the seven between-arm ratios | recomputation entry §6 | `~/ws_rank/RANK_RECOMPUTE.json` |
| 4.1 | the **3.295×** retraining floor, five same-seed repeats (**n = 5**) | `NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §1; `v2/research/rebase/d1_envelope_readout.py` | `~/e0_run/d1_envelope/rep{1..5}.npz`, `~/e0_run/d1_envelope_readout.log` |
| **4.1a** | **the floor audit — every rank comparison the paper makes or relies on, judged against the floor measured on its own block** | `v2/research/rebase/p2/floor_audit.json` (the list), `v2/research/rebase/p2/p2_floor_audit.py` (renders §4.1a's table and re-reads every value from its source), `v2/tests/test_p2_floor_audit.py` (fails if a ratio in the table disagrees with its source) | reads only `v2/research/rebase/p2/figures/data/` and repository markdown; no box access, nothing recomputed |
| **4.9a** | **the decorrelation reversal and the rank/cosine dissociation, three runs, one seed per level** | `NOTEBOOK_ENTRIES/lr_test_and_decorrelation_reversal_20260804T1130Z.md` §2; `v2/research/rebase/d1_momentum_probe.py`; figure `v2/research/rebase/p2/figures/fig_f9_decorrelation.py` | `~/e0_run/d1_diag/ablate_decorr{0.0,0.01,0.04}.log`, vendored and hashed at `v2/research/rebase/p2/figures/data/e0_run/d1_diag/` |
| 4.6a | the arm contrast on six target blocks, and the selection verdict under each as truth | `NOTEBOOK_ENTRIES/d2_coordinate_system_result_20260804T0800Z.md` §1, §1a; `v2/research/rebase/p2/p2_selection_rule_blocks.py` | `v2/research/rebase/nature/d2_coordinate_system/out/EXAM_PANEL.json` |
| **4.1–4.7** | **two independent recomputations agree.** (i) `~/ws_rank/` under the canonical `spectral.py`, 2026-08-04; (ii) `~/ws_p2/morpheus`, a fresh workspace verified byte-equal to HEAD *before* execution (402/402 files by git blob SHA-1) running the scripts vendored at `7b37dce`. Both reproduce §4.1, §4.2, §4.4(1), §4.5(b), §4.5(c), §4.6 and §4.7 to every published digit, and both independently identified the §4.5(a) statistic substitution **and the same corrected values** (3 of 6 → 2 of 6; D1 under canonical R2 is 3/3). | `NOTEBOOK_ENTRIES/effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`; `NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260804T0255Z.md` | `~/ws_rank/`, `~/ws_p2/out/` |
| 4.2 | variance decomposition and per-arm folds | metrics entry §4.1 | `~/e0_run/p2_necessity_and_variance.py` |
| 4.3 | step-200 in-flight R3, 6.05× vs 1.18× over five seeds | `~/e0_run/d1_v2/d1_{f,p}_seed{42,43,44}/` and `~/e0_run/d1_seeds4546/d1_{f,p}_seed{45,46}/train_metrics.jsonl`, `train_rank_tripwire_observed`; distilled to `v2/research/rebase/p2/figures/data/extracted/F3_TRIPWIRE_STEP200_R3_n5.json` | — |
| 4.4 | subsampling SDs; probe repeats; tolerance and centring insensitivity | metrics entry §3; `NOTEBOOK_ENTRIES/rank_probe_repeat_variance_20260804T0900Z.md`; recomputation entry §2, §4 | `~/e0_run/P2_METRICS_ALL_SUBSAMPLED.json`, `~/e0_run/d1_diag/probevar_*.log` |
| 4.5 | statistic / block / view flips | metrics entry §4.2, §4.3; recomputation entry §4, §7 | `~/e0_run/P2_ROBUSTNESS.json` |
| 4.7 | necessity test, predeclaration, escalation | `NOTEBOOK_ENTRIES/PREDECLARED_D1_necessity_test_20260803T2300Z.md`; metrics entry §5; `d1_readout_preregistration_20260803T1700Z.md` | `~/e0_run/d1_v2/`, `~/e0_run/d1_audit.log` |
| 4.7.2 | both paired bootstraps on the 40 untrained targets | `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json` | patient and cancer-cluster CI₉₅ per seed |
| 4.7.3, 6.3 | `random_control` arm gaps, absolute levels, and the qualified A3 verdict | D1 audit check A3; `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json`; T1.4's 76–82% random-gene-set finding | — |
| 4.8 | dilution table, both blocks, the 1.95–21.5× range | `DILUTION_LOWER_BOUND.md` §2, §6; `NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`; recomputation entry §5 | `p1_evidence/dilution/` **on the box**; `~/p1_out/dilution/dilution_foreign_tumour_pca256.npz`; build scripts `v2/research/dilution/` |
| 4.9 | D2, Phase 1b, "16/16", decorrelation | `D2_RESULT.md` §2, §4; `PHASE1B_TARGETED_READOUT.md` §3, §5, §7; `NOTEBOOK.md` 2026-08-02 01:20 UTC; `ENGINE_CLD.md` §1 + `HANDOFF_BUILD_AGENT.md` §1–2 | **decorrelation: no artifact; cited source does not exist** |
| 4.9, 4.10 | counter-measurement 12.88 → 1.00; collapse-regime values | `g26_rank_collapse_diagnosis_20260803T0500Z.md`; `d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`; `d1_programme_free_collapsing_in_training_20260803T1930Z.md`; `d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md` | `~/e0_run/d1_diag/`, `~/e0_run/d1_v1/` |
| 5.1 | four gate instances; gate non-reproducibility over 8 runs | `g26_variance_floor_fix`, `g26_stepbudget_sweep`, `g26_passes`, `d1b_premise_fails_all_five_arms_collapse`, `d1_relaunch`, `d1a_control_complete_and_gate_fails_2of3_in_runner`, `g26_is_not_reproducible_20260804T0700Z.md` | `~/e0_run/d1_diag/gatevar_{1..8}.log` |
| 5.2 | momentum sweep, staleness falsification, turnover | `queue_size_implicates_the_key_set_20260803T2200Z.md`; `momentum_rescues_rank_but_staleness_is_not_the_mechanism_20260803T2330Z.md`; `turnover_criterion_FALSIFIED_20260804T0330Z.md` | `~/e0_run/d1_diag/` |
| 5.2, 5.4 | **momentum seed replication**, three seeds per momentum, 500 steps, canonical R1 and R3 on the fixed held-out probe | `NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md` §3; `v2/research/rebase/d1_momentum_probe.py` | `~/e0_run/d1_diag/mseed_m{0,0.999}_s{42,43,44}.log` |
| 5.4 | the binary outcome the fix rests on: D1-A dispositions and epoch-40 collapse; D1-B completion, ranks, channel and intervals | `NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`; `d1b_blocked_gate_does_not_exercise_the_fix_20260804T0500Z.md`; §3.3, §4.7.2, §4.7.3 | `~/e0_run/d1_v1/`, `~/e0_run/d1_v2/`, `~/e0_run/d1_audit.log` |

## Appendix B — code index

| function | file:line | used by |
|---|---|---|
| `effective_rank` — **the only implementation**, canonical R1 by default, R2/R3 reachable as named variants | `v2/calibra/spectral.py` | every rank number in this paper |
| `RANK_VARIANTS` (R1/R2/R3 + uncentred / row-normalised combinations) | `v2/calibra/spectral.py` | §3.1, §4.5, §4.8 |
| single-definition + no-duplicate AST/SVD test | `v2/tests/test_effective_rank_canonical.py` | §3.1 |
| R3 retained under a calibrated abort threshold, canonical logged beside it | `v2/training.py:569` (tripwire), `v2/runner.py:942` (gate probe) | §3.1, §4.3, §5.1 |
| `heldout_top_cca` | `v2/calibra/spectral.py:78` | every channel number |
| `heldout_single_direction_correlation` | `v2/calibra/spectral.py:111` | detection-floor-scale controls |
| `cross_fitted_residuals`, `confound_design` | `v2/calibra/residualise.py` | all adjusted readouts |
| `permutation_null` | `v2/calibra/calibration.py` | all nulls in §3.2 |
| paired bootstrap | `v2/paired_bootstrap.py` | §4.6, §4.7 (pending), §4.9 |
| `symmetric_infonce` | `v2/losses.py:13` | §4.9, §5 |
| rank recomputation scripts, vendored | `v2/research/rebase/rank_recompute_all_instances.py`, `v2/research/rebase/rank_recompute_phase1b.py` (commit `8609081`) | §3.1, §4.1, §4.8, §4.9 |
| `stable_rank` | `v2/calibra/e1_rank_information.py` | §4.6 only |
| **the floor audit** — resolves every recorded comparison back to its source; **computes no rank statistic** | `v2/research/rebase/p2/p2_floor_audit.py`, `floor_audit.json` | §4.1a, and every section it audits |
| **F9** — the decorrelation dissociation | `v2/research/rebase/p2/figures/fig_f9_decorrelation.py` | §4.9a |

## Appendix C — the caveat that must travel with each number

Reproduced verbatim so that any future quotation can carry it.

- **The 3.295× retraining floor (§4.1).** **Five** identical retrains, one arm (`programme_only`, the
  *stable* one), one configuration, one seed, one stack, no interval. It cannot distinguish rank-specific
  variance from stack non-determinism, architecture or schedule. It is a **floor**, not an envelope, and
  its distribution is **bimodal** rather than a spread. Its effect on our own headline count (six of seven
  to seven of seven) runs in our favour and the seventh point clears it by 1.5%. The claim does not depend
  on it alone; §4.2 reaches the same conclusion from 8 within-arm degrees of freedom. The superseded n = 1
  estimate was 2.69× and no sentence in this paper should quote it as the floor.
- **Every count in §4.6.** Scored against **one** target block out of six. §4.6a re-scores against all six:
  every metric row's D2 count moves, the ordering between canonical effective rank and RankMe reverses on
  two blocks, and canonical effective rank reaches a nominally significant 6/6 on those two. **No
  selection-rule count in this paper may be quoted without its target block.** The D1 half was never
  re-scored, so its stability is unmeasured rather than established.
- **The variance decomposition (§4.2).** Four arms × three seeds = 12 artifacts, two experiments, one
  stack. Log-scale decomposition for rank-type metrics, raw scale for the channel.
- **The 6.05× / 1.18× pair (§4.3).** **Statistic R3**, in-training, states never saved, therefore
  `[NOT RECOMPUTABLE]` under R1 and never compared with an R1 number.
- **D2 (§4.6, §4.9).** *"REPORTED, NOT INTERPRETED (blocker 5). WSI states are ~0.80 collinear at init, so
  a narrower rank may reflect resistance to an already-collapsed view rather than dictionary content."*
  (`v2/research/rebase/d1_audit.py`, A5 note.) Two values for that collinearity exist in this evidence
  base and are **not the same measurement** — 0.7362 (std 0.0314) against RNA's 0.2740 on the full cohort
  (`NOTEBOOK.md:121`), and 0.80 against RNA's 0.62 on the fixed 16-patient gate batch. They must not be
  quoted interchangeably. The negative control is small but non-zero and same-signed, so of order 10–20%
  of the arm gap may be generic representation quality.
- **D1 (§4.7).** Quote **both** bootstraps or neither: decisive **3/3 on the patient** bootstrap and
  **2/3 on the cancer-cluster** bootstrap, seed 43's cluster interval reaching **+0.0006**. The cluster
  estimator is the conservative one and is the one to weight. `D1_PAIRED_BOOTSTRAP.json` (unstratified)
  must **not** be used — it scores all 90 non-control targets, 50 of which are `programme_only`'s own
  supervision; the stratified file on the 40 untrained targets is the only admissible one. Separately,
  the rank result is **3/3 under every canonical statistic** — R1, R2 and R3. The previously quoted
  "2/3 under R2" rested on a mislabelled statistic and is **withdrawn** (§4.5a).
- **The `random_control` verdict (§4.7.3, §6.3).** *"The instrument does not manufacture an arm
  difference; the absolute level is high and separately explained."* Never quote audit check A3 as an
  unqualified pass. Arm gap −0.0224 / −0.0072 / −0.0322, CIs spanning zero — but absolute level
  0.4425–0.4810 against real targets' 0.51–0.62 and the **D2 row-shuffle** null of **0.140**, i.e.
  ≈3.2× the floor, consistent with T1.4's 76–82%. **The comparator is 0.140, never 0.147** (§3.2's
  footnote). **No absolute channel level in this paper may be read against an assumed null of zero.**
- **Our own margin (§4.7.4, §6.3).** The smallest arm difference quoted (0.0705) is about **half** the
  real-versus-random-control margin (0.139). Any quotation of a D1 or D2 arm difference should carry
  this, for the same reason §4.1's envelope travels with every rank level.
- **Dilution (§4.8).** *"Single seed (42) and a single draw of donor assignments… gives no error bar on
  the level-to-level differences."* And: *"The whole curve is one representation."* Detection floor
  censored at ≥ 0.40 from d = 0.09; the source's own §4 **withdraws** the phrase "lower bound".
- **Phase 1b (§4.9).** *"`full` vs `programme_only` manifests were not verified as matched on
  epochs/LR/budget in this run (G0.4). Until they are, the rank comparison in §5 is suggestive, not
  causal."* Both differences are inside this stack's retraining noise.
- **"16/16" (§4.9).** The column is labelled *"`z_biology` matrix rank"* in its source. It is a **hard
  numerical rank** whose maximum is the batch size of 16, and the centred effective rank of the same
  objective falls 12.88 → 1.00. **The project's earlier description of this as a strong effective-rank
  instance is withdrawn.**
- **Decorrelation (§4.9, §4.9a).** Earlier codebase generation; benchmark statistic undefined in this
  repository; rank statistic unknown; **cited source file does not exist**. `[NOT RECOMPUTABLE]`.
  **And every claim this project makes about `feature_decorrelation` being defective — the collapsed
  global minimum, the self-extinction at every weight, the five-arm sweep of §5.1's instance 3 — was
  measured on a queue written by the query encoder and must be quoted with *"in the absence of a
  momentum key encoder"* attached.** At `m = 0.999` the same term raises rank monotonically
  (4.32 → 8.01 under R3, 6.29 → 12.20 under canonical R1) while the RNA-view mutual cosine also rises
  monotonically (0.4774 → 0.8696) on the same three runs — §4.9a and figure F9. That rank change is
  **1.85× / 1.94×, inside §4.1's 3.295× floor**: the monotonicity and the co-measured cosine carry
  the observation, not the size, and it is **one seed per level**.
- **The floor audit (§4.1a).** 56 rank comparisons, each judged against the floor measured on its
  own **statistic and** block; of the 25 selections between candidate configurations, **13 fail, 11
  cannot be judged at all, and 1 clears — and the one that clears is RankMe as published, not ours**.
  25 of the 56 rows sit on a block for which **no floor has ever been measured**, including every
  rank number in §5. **Unjudgeable is not a pass and not a failure**: the criterion has not been
  applied to those rows, and the audit records `clears: null` for them rather than a verdict.
  Machine-checkable at `v2/research/rebase/p2/floor_audit.json`, enforced by
  `v2/tests/test_p2_floor_audit.py`. Any quotation of a rank comparison from this paper should be
  checked against that table first — including the floor it is being judged against, which for eight
  of T1's twelve metric rows was measured only on 2026-08-04 and for four blocks does not exist.
- **D1-A's 9.81 / 1.71 (§4.9).** **Statistic R3**, `[NOT RECOMPUTED]` under R1. *"Nothing about programme
  supervision may be concluded from it — the contrastive arm never trained, so the comparison measures a
  defect, not an ablation."*
- **§5's momentum sweep (§5.2, §5.4).** The original sweep is **one seed per momentum value** — because
  the harness had its seed **hardcoded**, a defect rather than a judgement — centred R3, durability
  established only to step 600 against a 583-step training duration; its m = 0 against m = 0.999 ratio
  (2.64×) is **inside** §4.1's 3.295× disqualifying floor. **The seed replication has since run** (three
  seeds per momentum, 500 steps, canonical R1 and R3 on a **fixed held-out probe**, not §4.1's
  residualised exported block) and separates the two momenta 3/3 on both statistics, closing the
  single-seed defect — **but at a worst-case 3.29× against that same 3.295× floor**, so the *rank*
  difference is **not resolvable by this paper's own criterion**. Different arm, duration and block, so
  indicative rather than a like-for-like disqualification; all three mismatches point toward a larger
  floor in this regime, and no like-for-like floor has been measured. **The fix's justification does not
  rest on any rank ratio and must never be quoted as though it did**: it rests on the binary outcome of
  §5.4 — `programme_free` completing 40 epochs uncollapsed 0 of 3 seeds before the fix and 3 of 3 after,
  with a channel and paired bootstrap intervals where there had been no exported readout at all. **The
  choice of `m = 0.999` over `m = 0.99` (1.26×) is not supported by anything**; only momentum against
  none is. The anchoring account is **not confirmed** and should not be cited as established; the
  turnover sharpening was predeclared, tested and **refuted**.
- **All rank numbers.** RankMe's ε sits outside the division, so its statistic is not the exponential of a
  Shannon entropy and **no number in this paper is comparable to a published RankMe value**.
