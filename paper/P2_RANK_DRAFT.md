# Effective rank fails inside its own stated scope: six dissociations from a cross-modal biological pipeline, and a reproducibility floor that no rank comparison clears

**A negative result about a representation-geometry proxy, reported with its own weakest evidence and its own prior art first.**

*Working draft, 2026-08-04. Companion to `paper/P1_CALIBRA_DRAFT.md`, which asks what an **analysis**
would have missed; this asks what a **representation-geometry summary** fails to tell you. Every
number in this document traces to a named artifact, notebook entry or evidence file in this
repository; each table carries a `provenance` line. Numbers that do not exist are marked as not
measured rather than estimated. Every citation carries an explicit verification status (§2.6). Three
fabricated citations have previously contaminated this project; §2.6 is not a formality.*

> **Status.** Complete except for §4.8, which carries a marked `[D1 RESULTS PENDING]` slot. D1-B is in
> training at the time of writing (`~/e0_run/d1_v2/`: three arms at 40/40 epochs, one at 3/40, one at
> 4/40, one not started). Nothing else in this draft depends on it, and the conclusion does not change
> under any D1-B outcome — §4.8.4 states in advance what each outcome would and would not license.

---

## Abstract

Effective rank — the exponential of the Shannon entropy of a representation's L1-normalised singular
values (Roy & Vetterli, EUSIPCO 2007, pp. 606–610) — is proposed as a label-free criterion for
assessing self-supervised representations and selecting hyperparameters without labels (RankMe;
Garrido, Balestriero, Najman & LeCun, ICML 2023, arXiv:2210.02885). **We claim no discovery of the
negative.** Thilak et al. (arXiv:2312.04000) have already reported that RankMe correlates poorly with
downstream performance (Spearman 0.3174 on VICReg at 100 epochs) and that "a high rank does not
guarantee superior performance", and a RankMe co-author has since co-signed the statement that
"current methods like RankMe fail to adequately evaluate representation quality"
(arXiv:2410.04289). This paper confirms that finding in a domain it has not been tested in, and adds
three things the existing negatives do not supply.

First, a failure **inside RankMe's own stated scope.** RankMe restricts itself twice — high rank is
"necessary (but not sufficient)", and "RankMe should however only be used to compare different runs of
a given method". Most published negatives, and four of our own six dissociations, fall outside those
restrictions and therefore do not contradict it. One does not: in a three-seed supervision-target
ablation with arms matched by construction, one method, one architecture, in-distribution held-out
evaluation and 2,000-repeat paired patient- and cancer-clustered bootstraps on 2,766 patients, the
arm with the **higher** effective rank loses in one seed (34.12 vs 28.77, Δ = −0.1089) and the two
arms are **equal to two decimals** in another (9.11 vs 9.14) with the same arm losing by −0.1226. Both
confidence intervals exclude zero in 3/3 seeds. The rank ordering agrees with the performance ordering
in one seed of three.

Second, a **reproducibility floor** that we believe has not been reported and that constrains the
within-method use RankMe recommends. On this stack effective rank spans **9.11 to 34.12** across three
seeds of one configuration, and **8.68 versus 23.39** for one seed retrained with an identical
configuration — a 2.7× move — while the paired channel difference it would be used to explain is
stable at −0.1089 to −0.1325 with the same sign 3/3. A criterion used to select between configurations
must resolve an effect larger than its own run-to-run spread; here it does not.

Third, a **magnitude miscalibration** measurable only on a dose–response, not on a two-point
comparison. Under seven levels of information-free patch contamination of a representation with **zero
fitted parameters**, effective rank falls 196.2 → 161.2 (−18%) while the null-corrected cross-modal
channel falls 1.000 → 0.333 (−67%) — the two move together but are miscalibrated by a factor of 3.7.
This is consistent with, not contrary to, RankMe's necessary-not-sufficient framing, and we label it
as such.

We report the evidence that cuts against us at the same prominence. Centred effective rank **does**
fall to ~1 under total collapse (12.88 → 1.00 within 50 steps), so the collapse-diagnostic use is
supported. The instance previously described by this project as its most dramatic — rank pinned at
16/16 while patient-to-patient cosine reached 0.9999 and retrieval fell below chance to 0.000 — turns
out on inspection to be a *hard numerical* rank at a structural ceiling, not an effective rank, and we
withdraw the description. And at effective rank 9.1 of a nominal 256, two representations still read
held-out channels of 0.5983 and 0.4757 against a permutation null of 0.140, so "low rank means little
information" is also false outside total collapse. Separately we document that three mutually
incompatible statistics are implemented under the name `effective_rank` in the single repository that
produced these results, and that our own historical instances were not all measured with the same one.

### Short abstract (~200 words)

Effective rank (Roy & Vetterli 2007) is proposed as a label-free criterion for assessing
self-supervised representations and for label-free hyperparameter selection (RankMe, ICML 2023). That
it can fail is known: Thilak et al. (2023) report RankMe Spearman 0.32 against linear-probe accuracy
and that high rank does not guarantee good representations. We confirm this in cross-modal
morphology→transcriptome learning and add three things. (i) A failure **inside RankMe's own stated
scope** — one method, one architecture, arms matched by construction, three seeds, paired bootstraps:
the higher-rank arm loses (34.12 vs 28.77, Δ = −0.1089) and equal-rank arms differ by −0.1226, with
both CIs excluding zero 3/3 and rank ordering agreeing once in three. (ii) A **reproducibility floor**:
rank spans 9.11–34.12 across three seeds of one configuration and 8.68 vs 23.39 for one seed
retrained, while the paired difference is stable — so a within-method rank comparison cannot resolve
an effect smaller than 2.7×. (iii) A **magnitude miscalibration** visible only on a dose–response:
rank −18% while the channel falls −67%. Against this, centred rank does fall to ~1 under total
collapse, so the collapse-diagnostic use survives. We also document three incompatible statistics
named `effective_rank` in one codebase.

---

## 1. Introduction

### 1.1 A scalar that is cheap, label-free, and therefore load-bearing

Evaluating a learned representation properly requires labels, a downstream task and a held-out split.
All three are expensive, and in cross-modal biological settings the downstream label is often the
scarce quantity the representation exists to predict. This creates strong demand for a *label-free*
scalar computed on the representation matrix alone.

The spectrum supplies the obvious candidate. Given a patient × feature matrix, take its singular
values, normalise them to a probability vector, and report the exponential of their Shannon entropy —
a smooth interpolation between 1 (all mass on one direction) and the matrix's rank (a flat spectrum).
This is Roy & Vetterli's effective rank, and RankMe adopts it verbatim as a label-free criterion for
self-supervised representation quality and hyperparameter selection (§2.1).

The intuition licenses three distinct practices, and they are not equally defensible:

1. **As a collapse diagnostic.** "Rank has fallen to ~1, so the model has collapsed." Supported by
   the evidence here (§4.9), and the only use we end up recommending — with a caveat that a cheaper
   statistic does the job better.
2. **As a training objective or a target of regularisation.** Anti-collapse regularisers penalise
   off-diagonal feature covariance in order to raise the occupied dimensionality (VICReg, ICLR 2022;
   Barlow Twins). §4.5 reports what happened when we did this.
3. **As a comparative quality score** — "configuration A has higher effective rank than B, therefore
   A's representation is better". This is the use RankMe formalises, restricts, and this paper reports
   evidence against **within those restrictions**.

Use (3) is the one cheap enough to be tempting and the one that, if wrong, is silently wrong: a rank
number never fails, never returns `NaN`, and never announces that it was not measuring the thing you
cared about.

### 1.2 What is already known, and what we are not claiming

**We claim no discovery of the negative.** The prior art is explicit and we state it before our own
results rather than after.

Thilak et al. (arXiv:2312.04000) report that "RankMe correlates poorly with downstream performance for
most models", measure Spearman 0.3174 and Kendall τ 0.2056 for RankMe against linear-probe accuracy on
VICReg at 100 epochs, observe that "peak downstream performance occur[s] far from the point of maximal
rank", and state in their limitations that "a high rank does not guarantee superior performance". They
also give the constructive counterexample: appending independent noise dimensions raises rank while
adding no utility. Separately, a RankMe co-author has co-signed the statement that "current methods
like RankMe fail to adequately evaluate representation quality, making cross-validation without labels
infeasible" (arXiv:2410.04289). If our claim were "effective rank does not predict downstream
performance", it would not be novel.

**We also do not claim more than RankMe claims.** RankMe hedges in four ways that blunt most negative
results, including four of our own six (§4.1):

- Rank is **necessary but not sufficient**: "a necessary (but not sufficient) condition"; "having a
  high rank is a necessary condition for good downstream performance"; "maximal rank is only a
  necessary condition". *A high-rank representation with poor performance therefore does not
  contradict RankMe — RankMe predicts exactly that.*
- **Same-method comparisons only**: "RankMe should however only be used to compare different runs of a
  given method, since the embeddings' rank is not the only factor that affects performance."
- **No monotone transfer to other datasets**: "there is no inherent reason for the rank of embeddings
  to transfer in a monotonic way to them."
- **A named failure region**: "Except for some degenerate solutions at full-rank…"

And **the quality-proxy claim must not be attributed to Roy & Vetterli**, who propose effective rank
as a real-valued relaxation of rank for signal-processing optimisation and make no claim about
representation quality anywhere. Nor to Jing et al. (ICLR 2022), whose paper is diagnostic and
mechanistic about dimensional collapse and contains no sentence claiming the singular-value spectrum
predicts downstream performance. Both attributions appear in earlier drafts on this project and are
corrected here (§2.6).

### 1.3 What this paper adds

> **Effective rank fails inside RankMe's own stated scope, on a within-method comparison with matched
> arms, three seeds and intervals; and independently of whether it is right on average, it is less
> reproducible on this stack than the effect it would be used to select on.**

Three additions to the existing negative:

1. **An in-scope failure.** One method, one architecture, arms matched by construction, varying one
   supervision target, in-distribution held-out evaluation, non-degenerate ranks, three seeds with
   paired bootstraps (§4.3). This is the regime RankMe reserves for itself and it is the regime the
   existing negatives do not test: LiDAR compares across methods and across training checkpoints;
   Kulkarni et al. work on LLM unembeddings; Cheng works on continual-learning plasticity.
2. **A reproducibility floor on the statistic itself** (§4.7). We have not found this reported. It
   applies specifically to the within-method use RankMe recommends: if rank moves 2.7× when the same
   configuration is retrained with the same seed, no within-method comparison smaller than that can be
   read.
3. **A dose–response magnitude miscalibration** (§4.2), which a two-point comparison cannot show and
   which defeats the natural retreat position that rank is "at least a rough guide": rank and the
   channel move together and are miscalibrated by 3.7×.

Plus a new domain — cross-modal morphology → bulk transcriptome on human tumours — which we frame as
extension and replication, not discovery.

### 1.4 What this paper is not, and its relationship to P1

This is a methods and measurement paper. **No biological claim is made anywhere in it.** The
representation states audited here exist for other work on the same project; here they are specimens.

Its relationship to `paper/P1_CALIBRA_DRAFT.md` must be stated precisely, because the two share
evidence. P1 is an instrument paper about injection-certified transmission and detection floors for
confound-adjusted cross-modal analyses. P1 §4.11 currently carries a four-row table of rank
dissociations and closes it with: *"For this paper the point is narrow… A fuller treatment belongs to
a companion paper and is not claimed here."* This is that companion paper.

**Required resolution of the overlap, for the authors rather than the reader.** P1 §4.11 and P1 figure
F11 should be reduced to a two-sentence pointer — that a geometric quality metric is computed on the
representation rather than through the analysis pipeline whose null is in question, and therefore
cannot substitute for a sensitivity statement — with the table and the F11 panel moving here in full.
As currently drafted the two papers carry the same table with the same caveats, and P1's version now
contains two statements this paper withdraws (§2.6, §4.6). **§4 of this paper and §4.11 of P1 must not
both be submitted.** `paper/P2_FIGURES.md` §"Cross-paper deconfliction" states the four specific edits.

### 1.5 Contributions

In descending order of how well evidenced they are, and each labelled with whether it contradicts
RankMe *as stated* or only the looser practice.

1. **An in-scope within-method failure with intervals** (§4.3). Contradicts RankMe's body claims
   ("RankMe Consistently Predicts Downstream performances From Representations"; "a predictor of
   representations' performance"). Three seeds, matched by construction, both CIs excluding zero 3/3,
   rank ordering correct in 1/3. Collected as a routine audit check, not designed as a rank
   experiment.
2. **A reproducibility floor on effective rank** (§4.7). 3.7× across seeds of one configuration; 2.7×
   across retrainings at one seed; against a paired difference stable to ±0.012. Contradicts the
   practicability of RankMe's own recommended within-method use, without contradicting any claim it
   makes about expectation.
3. **A dose–response miscalibration of magnitude** (§4.2). Rank −18% against channel −67% over seven
   levels on a zero-parameter representation. **Consistent with** RankMe's necessary-not-sufficient
   framing; contradicts the informal practice of reading rank as a health indicator.
4. **The constraint on our own claim, reported at equal prominence** (§4.6, §4.9). Centred effective
   rank falls 12.88 → 1.00 and 67.55 → ~2 under total collapse; and at effective rank 9.1 of 256 the
   representation still reads a channel of 0.4757–0.5983 against a null of 0.140, so even "low rank
   means little information" fails outside total collapse.
5. **A withdrawal.** This project previously described "rank pinned at 16/16 while the representation
   collapsed to cosine 0.9999" as one of its two strongest instances. That column is a **hard
   numerical rank** at a structural ceiling equal to the batch size, and the centred effective rank of
   the same objective falls to 1.00. The description is withdrawn here (§4.6).
6. **Three mutually incompatible statistics named `effective_rank` in one repository** (§3.1), with
   our own historical instances not all measured with the same one. Reported rather than harmonised,
   because it is the most likely thing a referee would find and because it is itself evidence for the
   thesis: a scalar whose name is stable while its definition is not gets quoted across contexts it
   does not survive.

---

## 2. Related work

### 2.1 Effective rank, and its proposal as a label-free quality criterion

**Definition — VERIFIED at full text.** Roy & Vetterli, "The effective rank: a measure of effective
dimensionality", 15th European Signal Processing Conference (EUSIPCO 2007), Poznań, 3–7 September
2007, pp. 606–610; DOI 10.5281/zenodo.40328. Definition 1, verbatim: *"The effective rank of the
matrix A, denoted erank(A), is defined as erank(A) = exp{H(p1, p2,…,pQ)}, where H(p1,p2,…,pQ) is the
(Shannon) entropy given by H(p1,p2,…,pQ) = − Σ p_k log p_k"*, with *"the singular value distribution
p_k = σ_k / ‖σ‖₁"* and *"all logarithms are to the base e"*. They prove *"1 ≤ erank(A) ≤ rank(A) ≤
Q"*.

Two consequences for this paper. **(i)** Our `v2/calibra/spectral.py:14-29` implements exactly this —
raw singular values, L1-normalised, natural log — which retires the definitional defect recorded at
`v2/research/B1_ledger/1_collapse_remedy.md:26` (covariance eigenvalues) and at
`HANDOFF_BUILD_AGENT.md:123` ("that error reached a paper draft (mistake #1)"). **(ii)** Roy &
Vetterli's abstract is about making rank minimisation tractable — *"Since rank minimization is
generally not practicable owing to its integer nature, we propose a real-valued extension"* — and
contains **no claim about representation quality or downstream performance**. Any text on this
project attributing a quality-proxy claim to them is wrong and must be corrected.

**The quality-proxy proposal — VERIFIED at full text.** RankMe: Garrido, Balestriero, Najman & LeCun,
"RankMe: Assessing the Downstream Performance of Pretrained Self-Supervised Representations by Their
Rank", *Proceedings of the 40th International Conference on Machine Learning* (ICML 2023),
arXiv:2210.02885 (v3, 2023-06-26 is the camera-ready). It adopts Roy & Vetterli's definition
explicitly, with an ε inside the normalisation: `p_k = σ_k(Z)/‖σ(Z)‖₁ + ε`.

What it claims, verbatim, and the claim strengthens between abstract and body:

- Abstract, hedged: *"we develop a simple unsupervised criterion that is indicative of the quality of
  the learned JE-SSL representations: their effective rank"*; *"allows one to assess the performance of
  JE-SSL representations, even on different downstream datasets, without requiring any labels"*;
  *"RankMe can be used for hyperparameter selection with nearly no reduction in final performance
  compared to the current selection method that involve a dataset's labels"*.
- Body, stronger: section headings *"RankMe Consistently Predicts Downstream performances From
  Representations"* and *"RankMe Predicts Linear Probing performance Even on Unseen Datasets"*;
  *"RankMe accurately predicts a model's performance"*; *"RankMe computed on the embeddings of the
  source dataset is a predictor of representations' performance on target datasets"*.

What it **restricts**, verbatim — and this is the part a negative result must respect:

- *"a necessary (but not sufficient) condition"*; *"having a high rank is a necessary condition for
  good downstream performance"*; *"This further highlights how maximal rank is only a necessary
  condition for good performance."*
- *"RankMe should however only be used to compare different runs of a given method, since the
  embeddings' rank is not the only factor that affects performance."*
- *"there is no inherent reason for the rank of embeddings to transfer in a monotonic way to them."*
- *"Except for some degenerate solutions at full-rank, RankMe values correlate well with
  in-distribution … and out-of-distribution … classification performance."* (Figure 1 caption.)
- A self-reported miss: *"We see than RankMe can improve OOD performance for VICReg, but leads to a
  small drop for SimCLR"* (sic).

**Earlier proposal in the same family — PARTIALLY VERIFIED.** α-ReQ (Agrawal, Mondal, Ghosh &
Richards), assessing representation quality by eigenspectrum decay. The related preprint
"Investigating Power laws in Deep Representation Learning" (Ghosh, Mondal, Agrawal & Richards,
arXiv:2202.05808) states verbatim: *"proximity of α to 1, is strongly correlated to the downstream
generalization performance… α is computable from the representations without knowledge of any labels,
thereby offering a framework to evaluate the quality of representations in unlabelled datasets."*
`[COULD-NOT-VERIFY: the NeurIPS 2022 venue. It rests on LiDAR's bibliography entry (a secondary
source); the OpenReview record `ii9X4vtZGTZ` could not be retrieved. Do not state the venue as
verified.]` Note that RankMe itself reports α-ReQ failing: *"the power law prior of α-ReQ fails on the
embeddings and as such those results must be interpreted with care."*

### 2.2 Prior negative results — what we must not claim novelty over

**LiDAR — VERIFIED at full text, and the paper that most constrains our claim.** Thilak, Huang,
Saremi, Dinh, Goh, Nakkiran, Susskind & Littwin, "LiDAR: Sensing Linear Probing Performance in Joint
Embedding SSL Architectures", arXiv:2312.04000v1, 2023-12-07. `[UNVERIFIED: peer-reviewed venue. The
arXiv comments field says only "Technical report" and there is no journal-ref. Cite as an arXiv
preprint.]`

Verbatim:

- *"We observe that (1) RankMe correlates poorly with downstream performance for most models, with the
  exception of the worst performing model, and (2) LiDAR correlates highly with downstream performance
  for all models."* (Figure 1 caption.)
- *"…the rank of a representation inflates early on in training, only to rapidly diminish, with peak
  downstream performance occurring far from the point of maximal rank… This makes any standard data
  covariance rank based measures extremely limited in their ability to spot the point of optimal
  performance in a training run."*
- *"We highlight that a low regularization loss is achievable by random embeddings, which might be
  high rank, but devoid of any utility as for downstream tasks. A measure of representation quality
  which is based on covariance rank alone would therefore, theoretically, fail to capture this failure
  case."*
- *"These two factors are not necessarily causally linked; a high rank does not guarantee superior
  performance."* (Limitations.)
- Quantitatively: RankMe Spearman **0.3174**, Kendall τ **0.2056** for VICReg at 100 epochs, against
  LiDAR's 0.9161 / 0.8105.

**Further negatives — PARTIALLY VERIFIED (abstract only).**

- Otero, Mateus & Balestriero, "Self-Supervised Anomaly Detection in the Wild: Favor Joint Embeddings
  Methods", arXiv:2410.04289v1 (2024-10-05): *"current methods like RankMe fail to adequately evaluate
  representation quality, making cross-validation without labels infeasible."* Balestriero is a RankMe
  co-author, so this is a partial retraction of applicability by an original author.
- Kulkarni, Springer, Subramonian & Swayamdipta, "Disentangling Geometry, Performance, and Training in
  Language Models", arXiv:2602.20433v2: *"While the best-performing models often exhibit a high
  effective rank, this trend is not universal across tasks and training setups… none can reliably
  predict downstream performance."* Domain is LLM unembedding geometry.
- Cheng, "Local Redundancy: An Information-Theoretic Measure of Plasticity from Synthetic
  Memorization", arXiv:2607.13432v1: existing measures including effective rank *"lack theoretical
  grounding and correlate poorly with performance on new tasks."* Domain is continual-learning
  plasticity.

**Where this leaves us.** Every one of these tests rank *across* methods, *across* checkpoints, or in
domains other than joint-embedding representation selection. None tests the within-method,
matched-arm, fixed-architecture regime that RankMe reserves for itself, and none reports a
reproducibility floor on the statistic. Those two gaps are §4.3 and §4.7 and they are what this paper
is for. **This section must be read before §4; a version of this paper that presents §4 as a discovery
is not submittable.**

`[SEARCH INCOMPLETE — must be closed before submission.]` The prior-art sweep behind §2.2 ran through
the arXiv Atom API and OpenAlex only; the natural-language search path was unavailable. A
Semantic-Scholar citation-graph sweep of RankMe's citing papers has not been run. Treat §2.2 as a
lower bound on the prior art. If a paper reporting an *in-scope, within-method* RankMe failure exists,
contribution (1) of §1.5 collapses to a replication in a new domain and the paper must be reframed
around contribution (2).

### 2.3 Dimensional collapse and anti-collapse regularisation

**VERIFIED at full text.** Jing, Vincent, LeCun & Tian, "Understanding Dimensional Collapse in
Contrastive Self-Supervised Learning", ICLR 2022, arXiv:2110.09348. Definition, verbatim:
*"dimensional collapse, whereby the embedding vectors end up spanning a lower-dimensional subspace
instead of the entire available embedding space"*, diagnosed from the spectrum — *"A number of
singular values drop to zero, indicating collapsed dimensions."*

**Correction to earlier drafts on this project.** A full-text search of Jing et al. for sentences
relating the spectrum to downstream performance returns none. The paper is diagnostic and mechanistic
— its own conclusion identifies "two mechanisms causing dimensional collapse: strong augmentation and
implicit regularization" — and **contains no claim that rank predicts quality.** P1 §2.6 lists it
among proposals of "geometric proxies for representation quality"; that grouping is inaccurate for
this reference and should be amended.

**VERIFIED.** VICReg: Bardes, Ponce & LeCun, "VICReg: Variance-Invariance-Covariance Regularization
for Self-Supervised Learning", ICLR 2022, arXiv:2105.04906. Camera-ready abstract, verbatim: *"a
method that explicitly avoids the collapse problem with two regularizations terms applied to both
embeddings separately: (1) a term that maintains the variance of each embedding dimension above a
threshold, (2) a term that decorrelates each pair of variables"* (sic). **Note:** the arXiv metadata
abstract and the ICLR camera-ready abstract differ; quote whichever is cited and do not mix them.
VICReg's geometric terms are claimed to *prevent collapse*, not to *indicate downstream quality* —
which is exactly the gap LiDAR exploits with its high-rank-random-embedding argument.

**VERIFIED.** Wang & Isola, "Understanding Contrastive Representation Learning through Alignment and
Uniformity on the Hypersphere", ICML 2020, PMLR 119, arXiv:2005.10242. Verbatim: *"Extensive
experiments on standard vision and language datasets confirm the strong agreement between both metrics
and downstream task performance."* Correlational and empirical; not predictive and not causal.

`[UNVERIFIED]` Barlow Twins (Zbontar et al. 2021) — PMLR URL recorded in this repository, no quote
retrieved. `[UNVERIFIED]` LDReg (Huang et al., arXiv:2401.10474). In pathology specifically, the
Robustness Index (de Jong et al., 2025) and representational-similarity analysis (Mishra & Lotter,
2025) are proposed as confound-aware representation diagnostics; both were spot-check verified in a
2026-07-29 literature audit and not re-verified in this pass.

§4.5 reports what happened when we followed practice (2): adding a covariance-decorrelation term
raised effective rank by 107% while moving the benchmark statistic by 0.0001, and a later independent
measurement found that our implementation of the term has a **collapsed global minimum** and
self-extinguishes at every weight tested. We state that as an observation about our implementation,
not as a refutation of the published method, which we have not reimplemented faithfully or
benchmarked.

### 2.4 Momentum encoders and memory banks

A companion draft in this repository (`paper/QUEUE_ANCHORING.md`) reports a fix to a queued contrastive
objective and states that it falsifies a mechanism attributed to MoCo. The reference is now verified
and the companion draft must be corrected in three specific ways before it is submitted.

**VERIFIED at full text.** He, Fan, Wu, Xie & Girshick, "Momentum Contrast for Unsupervised Visual
Representation Learning", CVPR 2020, arXiv:1911.05722. The momentum-encoder argument, verbatim:

> *"Our hypothesis is that good features can be learned by a large dictionary that covers a rich set of
> negative samples, while the encoder for the dictionary keys is kept as consistent as possible despite
> its evolution."*
>
> *"removing the oldest mini-batch can be beneficial, because its encoded keys are the most outdated
> and thus the least consistent with the newest ones."*
>
> *"We hypothesize that such failure is caused by the rapidly changing encoder that reduces the key
> representations' consistency."*
>
> *"though the keys in the queue are encoded by different encoders (in different mini-batches), the
> difference among these encoders can be made small. In experiments, a relatively large momentum (e.g.,
> m = 0.999, our default) works much better than a smaller value (e.g., m = 0.9), suggesting that a
> slowly evolving key encoder is a core to making use of a queue."*

Three required corrections to `paper/QUEUE_ANCHORING.md` and `paper/LIVENESS_GATE_DESIGN.md`:
(i) add the identifier — the repository currently contains **no** arXiv number, DOI or venue for MoCo
anywhere, only a bare "He et al., 2020"; (ii) the account is advanced **twice as a hypothesis**
("Our hypothesis is that…", "We hypothesize that…"), never as an established mechanism, and a
falsification must say so; (iii) MoCo ties the argument specifically to queue use ("a core to making
use of a queue"), so a falsification is only in scope if it is in a queue setting — ours is, which
should be stated rather than left to be inferred. The substantive characterisation in
`NOTEBOOK_ENTRIES/momentum_rescues_rank_but_staleness_is_not_the_mechanism_20260803T2330Z.md` — that
MoCo argues from key-to-key inconsistency across a long queue — is **correct** as verified above.

This paper makes no claim about MoCo. The note is recorded here because §2.6 is a project-wide
requirement and this was the worst outstanding item on it.

### 2.5 What existing work supplies that would have answered our question better

LiDAR proposes a replacement criterion, not merely a critique, and reports Spearman 0.9161 where
RankMe reaches 0.3174. `[NOT MEASURED — we have not computed LiDAR, or any published alternative
criterion, on any representation in this paper.]`
A referee is entitled to ask why a paper reporting rank's failure did not evaluate the
published alternative; the honest answer is that the alternative was identified during this draft's
prior-art sweep and computing it requires a labelled probe on every artifact, which has not been run.
It is named in §5.2 as the most valuable missing measurement.

### 2.6 Reference verification status

Three fabricated citations have previously contaminated this project (`HANDOFF_BUILD_AGENT.md:156`).
The verification protocol is recorded in `NOTEBOOK_ENTRIES/winkler_prior_art_20260803T0120Z.md`:
retrieve the full text (not the abstract), make targeted passes over named sections plus every figure
caption and table header, tabulate what the paper *actually* reports, quote verbatim the closest it
comes to the claim being attributed to it, and state where it falls short.

| reference | status | retrieval | note |
|---|---|---|---|
| Roy & Vetterli, EUSIPCO 2007, pp. 606–610, DOI 10.5281/zenodo.40328 | **VERIFIED** | full-text PDF | definition matches `spectral.py` exactly; makes **no** quality claim |
| RankMe — Garrido, Balestriero, Najman & LeCun, ICML 2023, arXiv:2210.02885v3 | **VERIFIED** | full-text PDF | claims and restrictions both quoted in §2.1 |
| LiDAR — Thilak et al., arXiv:2312.04000v1 | **VERIFIED (content)**; `[UNVERIFIED: venue]` | full-text PDF | arXiv comments say "Technical report"; no journal-ref. Cite as preprint |
| Jing, Vincent, LeCun & Tian, ICLR 2022, arXiv:2110.09348 | **VERIFIED** | full-text PDF | contains **no** rank→performance claim; P1 §2.6's grouping is inaccurate |
| MoCo — He, Fan, Wu, Xie & Girshick, CVPR 2020, arXiv:1911.05722 | **VERIFIED** | full-text PDF | see §2.4 for the three corrections required elsewhere in the repo |
| VICReg — Bardes, Ponce & LeCun, ICLR 2022, arXiv:2105.04906 | **VERIFIED** | full-text PDF + arXiv API | arXiv abstract ≠ camera-ready abstract; do not mix |
| Wang & Isola, ICML 2020, PMLR 119, arXiv:2005.10242 | **VERIFIED** | full-text PDF + arXiv API | "strong agreement", correlational |
| Otero, Mateus & Balestriero, arXiv:2410.04289v1 | **PARTIALLY VERIFIED** | abstract only | body not read |
| Kulkarni et al., arXiv:2602.20433v2 | **PARTIALLY VERIFIED** | abstract only | LLM domain |
| Cheng, arXiv:2607.13432v1 | **PARTIALLY VERIFIED** | abstract only | plasticity domain; "ICML 2026 (Spotlight)" per comments field, not independently checked |
| α-ReQ — Agrawal, Mondal, Ghosh & Richards | **PARTIALLY VERIFIED**; `[COULD-NOT-VERIFY: NeurIPS 2022 venue]` | related preprint arXiv:2202.05808 full text; venue via LiDAR's bibliography | OpenReview `ii9X4vtZGTZ` blocked |
| Barlow Twins — Zbontar et al. 2021 | `[UNVERIFIED]` | — | PMLR URL only |
| LDReg — Huang et al., arXiv:2401.10474 | `[UNVERIFIED]` | — | *local* dimensionality; not the same construct |
| de Jong et al. 2025; Mishra & Lotter 2025 | spot-check verified 2026-07-29 | — | not re-verified in this pass |
| Prior-art census for §2.2 | **INCOMPLETE** | arXiv + OpenAlex only | natural-language search unavailable; citation-graph sweep not run |

**One implementation discrepancy to resolve.** RankMe uses `p_k = σ_k/‖σ‖₁ + ε` (ε outside the
division). Roy & Vetterli use no ε and adopt the `0 log 0 = 0` convention. `v2/calibra/spectral.py:25`
uses neither: it **filters** singular values at `> 1e-12` before normalising. On near-collapsed spectra
these three differ measurably. If any number in this paper is to be compared to a RankMe value, the
implementations must be reconciled first; no such comparison is made here.

---

## 3. Methods

### 3.1 Three statistics named `effective_rank`, and why it matters here

This repository implements three mutually incompatible functions under the name `effective_rank`.
They are not variants of one statistic; they have different ranges, different maxima, and different
sensitivity to the collapse mode this paper is about.

| # | definition | implementation | used for |
|---|---|---|---|
| **R1** | Roy & Vetterli exactly: `exp(−Σ pᵢ ln pᵢ)`, `p = σ/Σσ`, on the column-centred matrix, singular values filtered at `> 1e-12` | `v2/calibra/spectral.py:14-29` (declared "single source of truth"); byte-identical torch duplicates at `v2/run_rank_ablation.py:35-42` and `v2/tests/test_stress_collapse.py:23-35` | every CALIBRA readout: instances 2, 4 and 6 |
| **R2** | participation ratio of the centred singular values, `(Σσ)²/Σσ²` | `v2/research/rebase/d1_audit.py:149-153` | D1 audit check A5 — the statistic the pending §4.8 table will be produced with |
| **R3** | participation ratio of the centred singular values **after L2-normalising each patient row** | `v2/research/rebase/d1_geometry_probe.py:50-53` | all live-checkpoint geometry probes: 67.55, 9.81, 10.47, 1.71, 7.38, 1.76, and the momentum sweep in `paper/QUEUE_ANCHORING.md` |

Additionally, instance 3 (§4.6) is reported in its source not as an effective rank at all but as
"`z_biology` matrix rank" — a **hard numerical rank**, maximal at the batch size of 16.

**Consequences, stated plainly.** R1 is the published definition and is the only one comparable to any
external number. R2 and R3 are participation ratios, bounded by `min(n_patients, n_features)` rather
than by the number of non-negligible singular values, and R3's row normalisation removes the norm
variation R1 and R2 retain. Comparing an R3 value of 1.71 on 282 patients against an R1 value of
196.2 on 2,766 patients is meaningless. **Every quantitative statement in §4 names which statistic
produced it, and no statement is made that requires comparing across them.**

We report this as a finding rather than an embarrassment because it is evidence for the paper's own
thesis, and because it took reading three implementations side by side, while writing this paper, to
notice.

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
| **Within-cancer specificity** | the benchmark statistic of an earlier codebase generation; **definition not recoverable from any file now in this repository** | not recorded | §4.5, §5.4 |
| **Retrieval accuracy @1** | cross-modal nearest-neighbour retrieval within a fixed 16-patient batch | **0.062** = 1/16 | `NOTEBOOK.md` 2026-08-02 01:20 UTC |
| **In-batch InfoNCE** | symmetric cross-modal InfoNCE, temperature 0.07 | **ln 16 = 2.7726** (16-patient batch, frozen queue excluded from the candidate count); **ln 80 = 4.38** (earlier configuration, live 64-key queue contributing candidates); **ln 2576 = 7.854** and **ln 4310 = 8.369** at training scale | `v2/losses.py:13`; `NOTEBOOK.md:1554` |

**These chance levels must not be mixed.** `ln 80 = 4.38` belongs to the pre-fix gate configuration
and `ln 16 = 2.7726` to the post-fix one. Every InfoNCE number in §4 is quoted with its own chance
level attached.

**A null that is not zero.** The permutation null for a 16-component canonical correlation on ~2,700
patients sits at 0.140–0.147, because a multivariate maximum over 16 fitted directions is upward-biased
at finite *n*. A raw channel ratio therefore flatters the surviving signal and the null-corrected
column is the one quoted throughout. This is also why "both arms at chance" is not a testable
statement for these readouts, and why §4.3's negative control is stated as "controls must score below
real targets" — necessary and not sufficient
(`NOTEBOOK_ENTRIES/d1_readout_preregistration_20260803T1700Z.md:59-63`).

### 3.3 Cohorts, representations and artifacts

All measurements are TCGA. Two cohort configurations appear and are labelled at every table.

- **Maximal paired split — 6,427 patients** (3,118 train / 543 validation / 2,766 test), holding out
  whole cancers: 11 development cancers, 21 held out. Used by §4.2, §4.3, §4.8.
- **Earlier configuration — 2,530 held-out patients over 21 test cancers.** Used by §4.4.

Representations are one of:

- **Zero-parameter patch statistics.** `concat(mean, std)` over frozen H-Optimus-0 patch tokens
  (1,536-d), PCA-reduced to 256 dimensions refit on train rows only per level, retained variance
  0.879–0.923. Used by §4.2 — it has *no fitted parameters*, so the rank change and the channel change
  cannot be attributed to different training runs.
- **Trained `wsi_biology` states,** 256 output features, exported to `.npz` diagnostic artifacts with a
  `split` column; test rows only.

### 3.4 What "matched" means here, and the two places it fails

- **§4.3 (D2) is matched by construction.** A single `D2_PAIR_MANIFEST.json` enumerates the 40 common
  arguments; both arms record the same `pair_manifest_sha256`
  (`ce1352e0ac7a98334e4fada8178986e8413fac1046ebb67a96f5c3cbc7c2fb0b`) and the same
  `common_config_sha256` (`b7b2441fd9d03a3a00152027efe8c7ada3bedc48e7939f1dfc0b320b02adf1fb`). Both
  arms use `--objective-profile programme_only`; **they differ only in the supervision target table.**
  This is what makes §4.3 a *within-method* comparison in RankMe's sense.
- **§4.4 (Phase 1b) is not matched.** Its source states it: *"`full` vs `programme_only` manifests were
  not verified as matched on epochs/LR/budget in this run (G0.4). Until they are, the rank comparison
  in §5 is suggestive, not causal"* (`PHASE1B_TARGETED_READOUT.md:147-148`). The three diagnostic
  artifacts record only `configuration_sha256`, `git_commit`, `git_dirty`, `package`, `package_root`,
  `source_tree_sha256` — **no epochs, no learning rate, no token budget, no seed** — and `git_dirty` is
  `True` for all three (`HANDOFF_PHASE_D.md` §D1.0).
- **§4.5 cannot be assessed at all.** See §5.4.

### 3.5 Seed reproducibility on this stack, and the rule it forces

Training is **not seed-reproducible** on this hardware. Re-exporting a surviving checkpoint reproduces
its recorded readout to five significant figures (0.58612 against 0.5861 recorded), so the
export/readout path is deterministic. Retraining the same seed with the same configuration gives a
different model: held-out top-CCA 0.6214 versus 0.5861, and R1 effective rank 23.39 versus 8.68
(`D2_RESULT.md` §4).

**The rule this forces, and which this paper obeys: quote paired within-run differences, never
levels.** It applies to §4.3 and §4.8. It does *not* rescue effective rank, because rank is quoted as a
level in every practice this paper is about — including RankMe's own, which selects between runs by
comparing their rank values. There is no paired-difference form of "this run has higher rank". §4.7
develops this into the paper's second contribution.

### 3.6 Reporting rules adopted for this paper

1. §4.1's overview table carries an explicit strength ranking **and** a column stating whether the
   instance contradicts RankMe *as stated* or only the looser practice. Instances that do not
   contradict RankMe are labelled as such in their own section heading.
2. Every instance carries its caveats in the same table row, not in a later paragraph.
3. Where an instance's own source file corrects, caveats or withdraws it, that text is quoted.
4. Evidence that cuts against the claim is reported in §4.6 and §4.9 at the same prominence.
5. No number is compared across the three rank statistics of §3.1, and no number is compared to a
   published RankMe value (§2.6, implementation discrepancy).

---

## 4. Results

### 4.1 Overview: six dissociations, ranked, and what each does and does not contradict

Instance numbers 1–4 are inherited from `P1_CALIBRA_DRAFT.md` §4.11 so the two papers can be
cross-read; instance 5 is D1 (pending, §4.8) and instance 6 is the D2 audit (§4.3). The **strength**
column is our assessment of evidential weight and is *not* the instance number. The final column is
the one a referee will read first.

| # | manipulation | rank stat. | rank change | information change | strength | contradicts RankMe as stated? |
|---|---|---|---|---|:---:|---|
| **6** | supervision target: Hallmark vs perturbation-basis, one method, one architecture, 3 seeds, matched by construction | R1 | 23.39/28.77/9.14 vs 14.87/**34.12**/9.11 | Δ channel −0.1325 / −0.1089 / −0.1226; both CIs exclude zero 3/3 | **1st** | **Yes.** Within-method, in-distribution, non-degenerate ranks — RankMe's own reserved regime. Rank ordering agrees 1/3 |
| **4** | patch-bag contamination, d = 0 → 0.80, 7 levels, zero-parameter representation | R1 | 196.2 → 161.2, **−18%** | null-corrected channel 1.000 → 0.333, **−67%** | **2nd** | **No.** High rank with degraded information is exactly the necessary-not-sufficient case. Contradicts the informal health-monitor practice, and adds a magnitude miscalibration (3.7×) that a two-point comparison cannot show |
| **2** | objective profile: `full` → `programme_only` | R1 | 38.48 → 32.06, **−17%** | held-out top-CCA 0.4768 → 0.4748, **−0.002** | **3rd** | **No.** Different objectives are arguably different methods; and both differences are inside this stack's retraining noise (§3.5, §4.4) |
| **3** | full training schedule vs contrastive-only, 16-patient train batch | **hard matrix rank** | pinned at **16/16** in every arm | patient cosine 0.7089 → 0.9999; pos/neg 0.9959 vs 0.9960; retrieval 0.062 → **0.000** | **4th** | **No** — and it is not about effective rank at all. Hard rank at a structural ceiling of 16; the centred effective rank of the same objective falls **12.88 → 1.00** |
| **1** | covariance-decorrelation term added | unknown | 49.9 → 103.3, **+107%** | within-cancer specificity 0.1366 → 0.1367, **flat** | **5th** | **No.** Earlier codebase generation; benchmark statistic undefined; **cited source file does not exist in this repository** |
| **5** | D1: `programme_only` vs `programme_free` | R2 (pending) | `[D1 RESULTS PENDING]` | `[D1 RESULTS PENDING]` | — | see §4.8.4 |

*Provenance: instance 6 — `v2/research/rebase/nature/D2_RESULT.md` §2, §4; outputs
`~/e0_run/d2_v3/bootstrap/` and `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json` on persistent NFS.
Instance 4 — `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2, §6 and
`NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`; outputs `p1_evidence/dilution/`.
Instance 2 — `v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §3, §5, §7; run
`runs/calibra_v3_targeted`. Instance 3 — `NOTEBOOK.md` entry 2026-08-02 01:20 UTC and
`NOTEBOOK_ENTRIES/g26_variance_floor_fix_20260803T0210Z.md`; counter-measurement
`NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md`. Instance 1 —
`v2/research/rebase/ENGINE_CLD.md` §1 and `HANDOFF_BUILD_AGENT.md` §1–2 only; see §5.4.*

**Read this table with §2.2 in hand.** Four of six instances do not contradict RankMe as written.
That is not a weakness of the evidence; it is the reason §4.3 and §4.7 are the paper and the others
are context.

### 4.2 Instance 4 — a dose–response in which rank under-reports the loss by 3.7× *(does not contradict RankMe; contradicts the health-monitor practice)*

Patch bags were contaminated with same-cancer, different-patient tumour patches at seven nested
levels. The representation is `concat(mean, std)` over frozen H-Optimus-0 tokens with **no fitted
parameters**, so nothing here is a different training run. Both quantities are read from the same
representation at each level, through the same instrument.

| requested d | achieved d | adjusted top-CCA | held-out top-CCA | raw ratio | **null-corrected ratio** | detection floor | attenuation | **R1 effective rank** | perm *p* |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.000 | 0.5573 | 0.4932 | 1.000 | **1.000** | 0.20 | 1.130 | **196.2** | 0.0033 |
| 0.10 | 0.091 | 0.5571 | 0.5017 | 0.9996 | **0.999** | ≥ 0.40 | 0.985 | **194.1** | 0.0033 |
| 0.20 | 0.211 | 0.5447 | 0.5129 | 0.977 | **0.968** | ≥ 0.40 | 1.003 | **190.5** | 0.0033 |
| 0.30 | 0.302 | 0.5190 | 0.4986 | 0.931 | **0.905** | ≥ 0.40 | 1.057 | **187.5** | 0.0033 |
| 0.40 | 0.400 | 0.4774 | 0.4619 | 0.857 | **0.804** | ≥ 0.40 | 1.014 | **184.7** | 0.0033 |
| 0.60 | 0.600 | 0.3971 | 0.3680 | 0.713 | **0.607** | ≥ 0.40 | 0.855 | **176.5** | 0.0033 |
| 0.80 | 0.800 | 0.2844 | 0.1922 | 0.510 | **0.333** | ≥ 0.40 | 0.863 | **161.2** | 0.0033 |

*Provenance: `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2, §6;
`NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`. Cohort 6,427 patients, 238,610 tumour
patches, 7,644 slides; 2,766 evaluated on `test`. Instrument: 108-column cancer + pooled-TSS design,
seed 42, 16 components, 20 draws, 300 permutations (resolution 1/301 = 0.0033). Permutation null
median 0.145–0.147 at every level. Outputs under `p1_evidence/dilution/`. Rank statistic R1.*

**The result.** Over the full sweep effective rank falls **196.2 → 161.2, i.e. −18%**, while the
null-corrected channel falls **1.000 → 0.333, i.e. −67%**. Two thirds of the cross-modal information
is destroyed while the spectrum's occupied dimensionality is nearly preserved. Someone monitoring rank
would see an 18% drift and conclude the representation was substantially intact.

**What this does and does not contradict.** It does **not** contradict RankMe: high rank with degraded
information is precisely the necessary-not-sufficient case RankMe reserves, and LiDAR's noise-dimension
argument already establishes that high rank can be uninformative. What it adds is a **magnitude**: the
two quantities move *together* and are miscalibrated by a factor of **3.74** over the full sweep
(17.84% of the rank lost against 66.7% of the channel). Level by level, the ratio (fraction of rank
retained) / (fraction of channel retained) reads **1.000, 0.990, 1.003, 1.056, 1.171, 1.482, 2.467** —
essentially 1 for the first two levels and then rising steeply, i.e. rank tracks the channel while
almost nothing is happening and diverges from it exactly when the damage becomes real. (The dip to
0.990 at d = 0.091 is a single-seed wobble on a level where the channel changed by 0.001; nothing is
claimed about strict monotonicity of the ratio.) That defeats the natural retreat position — "rank is
at least a rough guide" — which neither the necessary-condition framing nor a correlation coefficient
addresses.

**Why it is nonetheless the second-strongest instance.** Three properties no other instance has
together. (i) It is a **dose–response over seven monotone levels**, not a two-point comparison.
(ii) The representation has **zero fitted parameters** — no training run, no seed, no optimiser
non-determinism, so §3.5 does not apply. (iii) Both quantities pass through the **same instrument at
the same level**, with attenuation 0.855–1.130 (≈ 1) throughout, so the instrument is not itself
degrading with dilution.

**What travels with it,** quoted from the source: the detection floor is **censored** — the grid tops
out at 0.40 and the floor reads 0.40 from d = 0.09 onward, so it is "≥ 0.40"; the transmission floor
reads 0.05 everywhere, the finest level, so it is censored from below too. *"The whole curve is one
representation"* — a trained attention aggregator could plausibly down-weight foreign patches, so the
number is a property of unweighted mean pooling, not of the modality. And it is a **single seed (42)
and a single draw of donor assignments**, which "gives no error bar on the level-to-level differences";
monotonicity over seven levels is what carries the result in the absence of intervals.

A labelling caveat the source imposes on itself and which travels with any quotation: the file is named
`DILUTION_LOWER_BOUND.md` and its own §4 **withdraws** the phrase "lower bound". The measured quantity
is *"the cost of preparation-matched, information-free contamination"*. Nothing in this paper turns on
which bound it is — both quantities are affected identically by the contaminant's nature.

### 4.3 Instance 6 — the rank ordering inverts inside RankMe's own reserved regime

Two arms were trained differing **only** in the molecular target table they are supervised on: curated
Hallmark pathway scores (arm H) versus 128 perturbation-basis coordinates (arm I). Both use
`--objective-profile programme_only`, the same architecture, the same 40 shared arguments and the same
manifest hash (§3.4). Three seeds each, 40 epochs, the 6,427-patient split. The readout is a
16-component held-out top-CCA with a 2,000-repeat paired patient- and cancer-clustered bootstrap on
2,766 held-out test patients, restricted to the **40 targets neither arm trained on**.

**This is a within-method, single-architecture, in-distribution comparison at non-degenerate ranks —
the regime RankMe reserves for itself** (*"RankMe should however only be used to compare different runs
of a given method"*).

| seed | channel, arm H | channel, arm I | Δ (I − H) | patient CI₉₅ | cancer CI₉₅ | **R1 rank, H** | **R1 rank, I** | higher rank | rank ordering correct? |
|---|---:|---:|---:|:---:|:---:|---:|---:|:---:|:---:|
| 42 | 0.6126 | 0.4800 | **−0.1325** | [−0.1605, −0.0993] | [−0.1792, −0.0632] | 23.39 | 14.87 | H | ✅ |
| 43 | 0.5970 | 0.4882 | **−0.1089** | [−0.1460, −0.0749] | [−0.1623, −0.0118] | 28.77 | **34.12** | **I** | ❌ inverted |
| 44 | 0.5983 | 0.4757 | **−0.1226** | [−0.1502, −0.0866] | [−0.1653, −0.0411] | 9.14 | 9.11 | ~equal | ❌ no signal |

*Provenance: `v2/research/rebase/nature/D2_RESULT.md` §2, §4. Run `d2_v3`; outputs under
`~/e0_run/d2_v3/bootstrap/` and `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json`. Rank statistic R1 on
the residualised held-out `wsi_biology` block, 256 nominal dimensions. Permutation null 0.140 for
every group at `permutation_p = 0.005`, the floor for 200 draws. Negative control on 90
`random_control` targets: Δ = −0.0099 / −0.0280 / −0.0268, 4–13× smaller, cancer CI covering zero 3/3
and patient CI 2/3.*

**The result.** A rank-based selection rule applied to these three seeds would pick the correct arm
once, pick the wrong arm once, and be unable to choose once. In seed 43 the losing arm has **higher**
effective rank (34.12 against 28.77) and still loses by −0.1089 with both CIs excluding zero. In seed
44 the two ranks are equal to two decimals (9.11 against 9.14) and the same arm still loses by
−0.1226. The performance ordering is stable in all three seeds; the rank ordering is not.

**Why this is the strongest instance.**

- **It is in scope.** One method, one architecture, matched by construction, one variable changed,
  in-distribution held-out evaluation, ranks well away from the full-rank degeneracy RankMe excludes.
  Neither LiDAR nor the other published negatives test this regime.
- **It has intervals on the information side.** Three seeds, paired patient- and cancer-clustered
  bootstraps, both CIs excluding zero 3/3. No other instance in this paper has that.
- **It was not designed to produce this conclusion.** The rank column exists because it is audit check
  A5 of a predeclared post-training checklist whose instruction reads *"reported, **not
  interpreted**"* (`NOTEBOOK.md:1114`), added for an unrelated worry — that WSI biology states are
  strongly mutually collinear at initialisation, so a narrower rank might reflect resistance to an
  already-collapsed view. The numbers were recorded before anyone asked whether they would order
  correctly. *(Two values for that collinearity appear in this evidence base and are not the same
  measurement: **0.7362 (std 0.0314)** against RNA's 0.2740 on the full cohort, measured 2026-08-01
  and recorded at `NOTEBOOK.md:121`; and **0.80** for WSI against 0.62 for RNA on the fixed
  16-patient gate batch, `NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md`. The A5
  note quotes the second. Nothing in this paper depends on which is used, but they should not be
  quoted interchangeably.)*
- **It contains its own within-seed control.** Seed 44's arms are equal in rank and differ by −0.1226
  in channel — a matched-rank pair with a large, interval-backed information difference. That single
  line is the cleanest refutation of "rank explains the gap" in this evidence base.

**What travels with it.** The negative control is small but not exactly zero and points the same way,
so of order 10–20% of the arm gap may be generic representation quality rather than supervision
content — which, if anything, tightens the argument, since that is the component a rank explanation
would be entitled to. Per §3.5, no individual point estimate here is reproducible from the seed alone;
only the paired within-run difference is quoted, and the rank column is quoted as a **level** precisely
because levels are what rank practice depends on. Three seeds is a small number; the finding is that
rank ordering is unreliable, not a quantified error rate.

**A second, independent observation from the same table** is developed in §4.7: the rank column spans
9.11 to 34.12 across three seeds of one configuration while the channel difference is stable at
−0.1089 to −0.1325.

### 4.4 Instance 2 — a 17% rank change with the channel unmoved *(and both differences inside the noise)*

Two exported artifacts of one architecture, differing in objective profile, scored through the CALIBRA
targeted readout on 2,530 held-out-cancer TCGA patients over 21 test cancers, with 13 spike levels, 40
draws, 2,000 permutations, 16 components per side and a 99-column cancer + pooled-TSS adjustment,
seed 42.

| run | state | held-out top-CCA | **R1 effective rank** | detection floor | attenuation | perm *p* |
|---|---|---:|---:|---:|---:|---:|
| `full` | `wsi_biology` | 0.4768 | **38.48** | 0.2 | 1.086 | 0.0005 |
| `programme_only` | `wsi_biology` | **0.4748** | **32.06** | 0.2 | 1.103 | 0.0005 |
| **change** | | **−0.002** | **−17%** | — | — | — |
| `full` | `wsi_identity` | 0.5393 | 191.07 | 0.3 | 1.228 | 0.0005 |
| `identity_only` | `wsi_identity` | 0.5393 | 191.07 | 0.3 | 1.228 | 0.0005 |

*Provenance: `v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §3, §5, §7; run
`runs/calibra_v3_targeted`. `permutation_p = 0.0005 = 1/2001` throughout. Rank statistic R1.*

**The result as originally read.** The representation changes materially between arms — direct array
comparison gives `max|diff| = 1.4e-01` on `wsi_biology` — loses 17% of its effective rank, and its
molecular channel is unchanged at the second decimal.

**The `wsi_identity` rows are the instance's own internal control.** That head is the frozen MLP-CLIP
teacher passed through: `max|diff| = 2.6e-04` between arms, i.e. numerically invariant, and both its
rank and its channel are identical to four significant figures. So the instrument does return "nothing
changed" when nothing changed, in the same run on the same artifacts, which rules out the reading that
the flat channel is an insensitivity of the readout.

**Why it ranks third, and the honest re-reading.** Three caveats, all from the source file.

1. **The arms were never verified matched** on epochs, learning rate or step budget (gate G0.4).
   Quoted: *"Until they are, the rank comparison in §5 is **suggestive, not causal**."*
2. **Single seed, and no interval on the difference.** Quoted: *"No CI on any between-run difference;
   a paired bootstrap on the biology gap is still required before the C2 numbers are quoted as a
   difference."* Given §3.5 — the same seed retrained moves a channel estimate by 0.035 and a rank
   estimate by 2.7× — a 0.002 channel difference and a 6.42-point rank difference are **both** inside
   this stack's run-to-run noise. **Instance 2 should therefore be read as showing that neither
   quantity is resolvable at one seed, not as showing that one moved and the other did not.** We state
   this rather than letting "unchanged" do work it has not earned, and we retain the instance because
   it is the origin of the project's "both directions" framing and withdrawing it silently would
   misrepresent how the claim was arrived at.
3. **A claim from the same run has already been withdrawn.** The neighbouring F2 claim — "the head
   trained for biology is worse at biology than the head trained for identity" — was withdrawn when E3
   showed the identity head to be a passed-through frozen teacher
   (`PHASE1B_TARGETED_READOUT.md` §4). That is why this paper treats the run's other between-arm
   statements conservatively.

### 4.5 Instance 1 — the weakest, and why we still report it

Adding a covariance-decorrelation term to the biology head raised effective rank **49.9 → 103.3**
(+107%, ~2.1×, reported over 3 seeds) while the benchmark statistic — within-cancer specificity —
moved **0.1366 → 0.1367**.

*Provenance: `v2/research/rebase/ENGINE_CLD.md` §1; `HANDOFF_BUILD_AGENT.md` §1 and §2.*

**This is the weakest item in the paper and must be labelled as such wherever it appears.** Four
reasons:

1. **The primary source does not exist.** `HANDOFF_BUILD_AGENT.md:98` cites `paper/.../RESULTS.md`.
   There is no `RESULTS.md` anywhere in this repository. The numbers survive only as prose in two
   summary documents written after the fact. No per-seed values, no run output, no artifact.
2. **The benchmark statistic is not recoverable.** "Within-cancer specificity" is defined in no file
   now present and is none of the statistics in §3.2. We cannot state what 0.1366 measures, what its
   chance level was, or whether it could resolve a change.
3. **The rank statistic is unknown.** It predates `spectral.py`'s consolidation and cannot be assigned
   to R1, R2 or R3.
4. **It comes from an earlier codebase generation** whose defects are documented — the same generation
   produced the instrument failure that is P1's motivating example.

**Why we report it anyway.** It is the observation that started this line of work — `ENGINE_CLD.md` §1
names it O1, "we manufactured capacity with no information in it" — and omitting it would present the
project's reasoning as cleaner than it was. It is **provenance for the hypothesis, not evidence for the
conclusion**, and it is excluded from every count and summary statement in this paper. If the paper is
submitted, instance 1 belongs in a history paragraph, not a results table.

**One thing it does contribute, from a later and better-sourced measurement.** The same decorrelation
term was independently characterised on a live objective: its global minimum **is** a collapsed state
(38.97 for a healthy batch against 1.19e-17 for an all-identical one); it collapses at every weight
from 0.001 to 4.0 within 25 steps while driving itself to zero (20.74 → 0.00); and a per-dimension
variance floor provably cannot stop it, because the rank-1 family `zᵢ = m + aᵢ·u` satisfies one
(`NOTEBOOK.md` decision of 2026-08-03; `NOTEBOOK_ENTRIES/g26_term_isolation_20260803T0930Z.md`;
`NOTEBOOK_ENTRIES/g26_centring_fix_20260803T0730Z.md`). **We state this about our implementation, not
about VICReg or Barlow Twins, which we have not reimplemented faithfully or benchmarked.** A
regulariser introduced to raise rank whose own minimum is collapse is a sharper cautionary tale than
the +107% number ever was.

### 4.6 Instance 3 — a withdrawal: what "pinned at 16/16" actually was

This is the most dramatic *collapse* in the evidence base and the **weakest of the four historical
*rank* instances**. This project previously listed it among its two strongest. We withdraw that
description here, before the numbers rather than after.

Three arms were trained on one fixed 16-patient batch from the real cohort, 800 steps, seed 42, hidden
512 / 4 layers / 8 heads, programme head 256, frozen 64-key memory queue, differing only in loss
weights.

| arm | in-batch InfoNCE (chance 2.7726) | retrieval acc@1 (chance 0.062) | cross-modal pos cos | cross-modal neg cos | WSI within-modality off-diag cos | `z_biology` **matrix rank** |
|---|---|---|---|---:|---:|:---:|
| **A** full `programme_free` schedule | 3.4681 → **2.7734** | 0.062 → **0.000** | 0.0538 → 0.9959 | 0.0816 → 0.9960 | 0.7089 → **0.9999** | **16/16** |
| **B** contrastive only | 3.4681 → **2.0789** | 0.062 → **0.188** | 0.0538 → 0.9988 | 0.0816 → **0.4922** | 0.7089 → **0.4946** | **16/16** |
| **C** contrastive only, lr 3e-3 | 3.4681 → 2.7724 | 0.062 → 0.125 | 0.0538 → 0.9998 | 0.0816 → 0.9997 | 0.7089 → 0.9998 | **16/16** |

*Provenance: `NOTEBOOK.md` entry 2026-08-02 01:20 UTC, source `scratchpad/collapse_diag.py` on the
A100; `NOTEBOOK_ENTRIES/g26_variance_floor_fix_20260803T0210Z.md`. Raw inputs are separable: raw
H-Optimus patient-mean off-diagonal cosine 0.3265, raw RNA 0.0587.*

**Arm A is total representational collapse.** Every patient's `z_biology` converges to the same vector
(off-diagonal cosine 0.9999); cross-modal positive and negative pairs become indistinguishable (0.9959
versus 0.9960 — the negatives marginally *higher*); retrieval falls to 0.000, **below** its chance
level of 0.062. Extended to 5,000 steps a sibling arm reaches patient cosine 1.0000 with InfoNCE pinned
at exactly ln 16 = 2.7726 from step 2,500 onward
(`NOTEBOOK_ENTRIES/g26_stepbudget_sweep_20260803T0340Z.md`). No measure applied to this representation
finds information in it.

**And the rank column reads 16/16 in every arm, including the collapsed one.** Taken at face value
this is the most striking single line the project has produced. **Four reasons it must not be taken at
face value.**

1. **That column is a hard numerical matrix rank, not an effective rank.** Its source labels it
   "`z_biology` matrix rank". It is none of R1, R2 or R3 (§3.1). A hard rank counts non-zero singular
   values, and a 16 × 256 float matrix has full row rank under essentially any perturbation, including
   the residual variation that survives a collapse to cosine 0.9999. "Hard rank is insensitive to
   near-collapse" is true, unsurprising, and not a finding about effective rank.
2. **Its maximum is 16 because the batch is 16** — a structural ceiling set by the experiment's
   design, not by the representation.
3. **The centred effective rank of the same objective does fall — to 1.00.** A later diagnostic on the
   same gate, same cohort, same seed, running clean in-batch InfoNCE, records `eff-rank 12.88 → 1.00
   by step 50`, with positive and worst-negative cosines converging to 0.9993 and minimum margin to
   −0.0001 (`NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md`). **Effective rank
   detected this collapse.** Reporting instance 3 as an effective-rank failure while that measurement
   sits in the same evidence base would be a selection error.
4. **It is a train batch, not held-out,** with n = 16.

**What instance 3 does establish, and it is worth keeping.** Two things, neither of them "effective
rank fails". First, **hard/numerical rank is worthless as a collapse diagnostic** — and it is
nonetheless what a `matrix_rank` call returns, and `d1_geometry_probe.py:53` computes and prints one
alongside the effective rank, so this is a live confusion in our own tooling. Second, **effective rank
works in this regime**, which is the boundary of our claim; instance 3 is the strongest *positive*
evidence in this paper for the collapse-diagnostic use defended in §4.9.

**Recorded as a correction.** P1 §4.11 lists instance 3 as one of the two "strongest" instances and
describes it as "rank pinned while information collapses". On the evidence assembled here that
description is not sustainable for effective rank, and the P1 text must be amended when P1 §4.11 is
reduced per §1.4.

### 4.7 Effective rank is less reproducible than the difference it would select on

Two measurements, both from `D2_RESULT.md` §4, on the R1 statistic. This is the contribution we
believe is not in the prior art, and it applies specifically to the within-method use RankMe
recommends.

**Across seeds of one configuration.** Arm H's effective rank reads 23.39 / 28.77 / 9.14 across seeds
42 / 43 / 44 — a **3.1× spread**. Arm I reads 14.87 / 34.12 / 9.11 — a **3.7× spread**. Over both arms
and all seeds the statistic spans **9.11 to 34.12**. Across exactly those runs the paired channel
difference is **−0.1325 / −0.1089 / −0.1226**: a spread of 0.024 on a mean of −0.121, with overlapping
CIs and the same sign 3/3.

**Across retrainings at a fixed seed.**

| | held-out top-CCA | **R1 effective rank** |
|---|---:|---:|
| recorded in the original run, arm H seed 42 | 0.5861 | — |
| re-export of the surviving arm-H seed-42 checkpoint | **0.58612** | **8.68** |
| **retrained** arm H seed 42, identical configuration | **0.6214** | **23.39** |

Re-export is deterministic to five significant figures; the export and readout path is not the source
of the variance. Retraining with the same seed and the same configuration is not deterministic: the
channel moves by 0.035 and the effective rank by **2.7×**.

**The consequence for RankMe's own recommended use.** RankMe selects between runs by comparing their
rank values, and restricts that comparison to runs of a given method — which is exactly the setting
above. A criterion whose value moves 2.7× when the same configuration is retrained cannot resolve a
between-configuration difference smaller than that. In §4.3 the largest rank difference on offer is
1.6× (23.39 against 14.87 in seed 42), comfortably inside the statistic's own reproducibility envelope
— which is a second, independent reason seed 42's "correct" ordering carries no weight.

**Stated as a rule a practitioner can apply:** *before using a rank difference to select between
configurations, retrain one configuration and measure the rank spread. If the between-configuration
difference does not exceed it, the comparison is uninformative.* We have not seen this check proposed
anywhere, and on this stack it would have disqualified every rank comparison this project made.

**The asymmetry is deliberate and we state it.** The channel readout is **not** exempt from §3.5
either — 0.5861 versus 0.6214 is a 6% move on the same seed. But the channel is quoted as a *paired
within-run difference* and rank is not, and the paired difference is stable where both levels are not.
That is a statement about how the two quantities are used, and it is the practical heart of this
paper.

**What is not established.** Two configurations and one stack. We do not know whether this variance is
a property of effective rank, of this hardware's non-determinism, of this architecture, or of the
40-epoch schedule. A controlled repeat design — N retrainings of one configuration, rank and channel
measured on each — has not been run and is named in §5.2. A short-horizon controlled probe on this
same stack finds rank *tight* (three repeats at 200 steps: 7.15 / 6.92 / 7.25, relative spread 4.7%;
`NOTEBOOK_ENTRIES/rank_probe_repeat_variance_20260804T0900Z.md`), which shows the variance is not
intrinsic to the measurement and is accumulated over training. That distinction matters and we do not
have it resolved.

### 4.8 `[D1 RESULTS PENDING]`

> **This section is a placeholder with a fixed shape. It will be filled from D1-B when the run
> completes. Nothing elsewhere in this paper depends on it, and §6 does not change under any outcome.
> §4.8.4 states in advance what each outcome would and would not license, so the section cannot be
> written to suit the result.**

#### 4.8.1 What D1 is

D1 trains two arms from a single command differing only in `--objective-profile`: `programme_only`
(biology head supervised on 50 Hallmark pathway scores) and `programme_free` (no programme regression;
a patient-paired cross-modal InfoNCE on the biology view instead), seeds 42/43/44, 40 epochs, on the
6,427-patient maximal split. It is the arm the withdrawn F2 claim never had
(`PHASE1B_TARGETED_READOUT.md` §4: *"there is no biology head trained without programme supervision
anywhere on disk"*), and it is the only place on this project where a rank gap and a channel gap will
be measured on arms **matched by construction, with three seeds and a paired bootstrap on both**.

Note for §4.1: because the two arms differ in *objective*, D1 sits closer to a between-method
comparison than §4.3 does, and therefore lands further outside RankMe's stated scope. §4.8.4 accounts
for this.

#### 4.8.2 What is already measured, from D1-A, and why it is not the result

D1-A was the first launch. Its `programme_free` arms were subject to a collapse defect since diagnosed
and fixed, so **D1-A is a control documenting the defect, not the ablation.**

Measured on the final epoch-39 checkpoints, 282 held-out test patients, **statistic R3**
(`v2/research/rebase/d1_geometry_probe.py`):

| arm | seed | epoch | **R3 effective rank** | hard rank @1e-3 | RNA–RNA mutual cos | feature std |
|---|---:|---:|---:|---:|---:|---:|
| `programme_only` | 42 | 39 | **9.81** | 136 | 0.258 | 0.0056 |
| `programme_only` | 43 | 39 | **10.47** | 158 | 0.265 | 0.0063 |
| `programme_free` | 42 | 39 | **1.71** | 11 | 0.986 | 0.0156 |

*Provenance: `NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`.
Statistic R3; maximum is `n_patients` = 282. At epoch 21 the same probe read 7.38 / 7.35 for
`programme_only` and 1.76 for `programme_free`.*

The gap is **5.7–6.1×**, with the two `programme_only` seeds agreeing (9.81 / 10.47).

**Four reasons this is not instance 5 and must not be quoted as it.**

1. **The `programme_free` arm never trained.** Seeds 43 and 44 were refused by the in-runner liveness
   gate (contrastive 0.50883 and 2.14122 against a threshold of 0.10; the latter close to chance
   ln 16 = 2.7726), and `programme_only` seed 44 stopped at 27 epochs. There is one `programme_free`
   seed, not three.
2. **The source entry forbids the reading.** Quoted: *"Nothing about programme supervision may be
   concluded from it — the contrastive arm never trained, so the comparison measures a defect, not an
   ablation."*
3. **No channel was measured.** `run_d1` raises on the first non-zero return code, so **no exports, no
   CALIBRA readout and no bootstrap were produced.** The information side of instance 5 does not exist
   in any form.
4. **A rank gap between a trained arm and a collapsed arm is not a dissociation.** It is the
   collapse-diagnostic use working correctly (§4.9). If it is all D1 ever delivers, it belongs in §4.6
   as further support for the one-sided claim, not in §4.1 as a sixth dissociation.

#### 4.8.3 The slot

D1-B is the rerun with the queue fix (a momentum key encoder). At the time of writing, from
`~/e0_run/d1_v2/`: `d1_p_seed42` 40/40, `d1_p_seed43` 40/40, `d1_f_seed42` 40/40, `d1_f_seed43` 3/40,
`d1_p_seed44` 4/40, `d1_f_seed44` not started.

**The table to be filled**, on the 40 targets neither arm trained on (`heldout_pathway` + `immune_tme`
+ `tumour_state`), with the 90 `random_control` targets as the negative control, per the
preregistration in `NOTEBOOK_ENTRIES/d1_readout_preregistration_20260803T1700Z.md`:

| seed | channel, `programme_only` | channel, `programme_free` | Δ | patient CI₉₅ | cancer CI₉₅ | rank, `programme_only` | rank, `programme_free` | rank ratio | ordering correct? |
|---|---:|---:|---:|:---:|:---:|---:|---:|---:|:---:|
| 42 | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` |
| 43 | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` |
| 44 | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` | `[PENDING]` |

*Provenance to be recorded on completion: `~/e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json`
(headline), `D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json` (negative control), `D1_PAIRED_BOOTSTRAP.json`
(secondary — **do not headline**, it scores all 90 non-control targets of which 50 are
`programme_only`'s own supervision), `D1_AUDIT.json` check A5 for the rank column, and
`D1_READOUT_INDEX.json`. The rank column will be **statistic R2** (`d1_audit.py:149-153`), which is
**not** the statistic that produced the 9.81 / 1.71 numbers in §4.8.2 (R3). The two must not be
compared; if both are reported the table must carry both labels.*

#### 4.8.4 Preregistered reading of each outcome

| outcome | what it licenses | what it does not |
|---|---|---|
| Both arms train; large rank gap, **large** channel gap in the same direction | Rank tracked information in this instance. **Report as a negative for this paper's generality**, in §4.1 and the abstract, at the same prominence as the six that dissociate. | It does not overturn §4.3 (in-scope, three seeds) or §4.7 (reproducibility floor), and one agreeing instance does not restore a proxy. |
| Both arms train; large rank gap, **small or absent** channel gap | A seventh dissociation, and the first with matched-by-construction arms, three seeds and a paired bootstrap on **both** quantities. Would be promoted to 1st in §4.1 — *but see the scope note*: the arms differ in objective, so it lands further outside RankMe's stated scope than §4.3, and §4.3 stays the in-scope result. | Still one architecture, one cohort, one modality pair. |
| `programme_free` collapses again (R2 rank ≈ 2, channel at or near the permutation null ~0.14) | Further support for the **collapse-diagnostic** use (§4.9) and nothing more. Belongs in §4.6. | **No** statement about programme supervision, per §4.8.2 item 2. |
| Fewer than three `programme_free` seeds complete | Report the completed seeds with the count in the table caption and mark instance 5 as under-powered. | Do not quote a two-seed gap as if it were three. |

**A standing constraint that applies to every outcome.** The liveness gate that admits arms to D1 has
been measured to be non-reproducible: eight runs with identical inputs, identical seed and an identical
2,400-step budget gave `final_biology_contrastive` values spanning **650×** (0.00859, 0.01076, 0.01770,
0.02019, 0.02407, 0.03266, 0.38009, 5.58511), a **6/8 pass rate** against an unchanged 0.10 threshold,
and a bimodal shape — six clustered within 4× and two divergent
(`NOTEBOOK_ENTRIES/g26_is_not_reproducible_20260804T0700Z.md`). At a 75% per-arm pass rate the
probability of all three contrastive arms clearing the gate is **0.42**. Any D1 table must state which
arms were admitted, because admission is itself a stochastic filter and arms that fail it are not a
random sample.

### 4.9 The use that survives, and where its boundary actually is

**Rank does detect total collapse.** Every measurement here of a representation collapsed to a single
direction shows a centred effective rank at or near 1–2:

| regime | rank statistic | value | co-measured evidence of collapse |
|---|---|---:|---|
| clean in-batch InfoNCE, step 50 | centred eff-rank (diagnostic script) | 12.88 → **1.00** | pos cos 0.9993, worst-neg cos 0.9993, min margin −0.0001 |
| `programme_free` at training scale, step 150 | R3 | 67.55 → **~2** | RNA-view mutual cosine 0.9813 |
| `programme_free`, epoch 21 / epoch 39, 282 held-out patients | R3 | **1.76** / **1.71** | RNA–RNA mutual cosine 0.977 / 0.986; hard rank 9 / 11 |

*Provenance: `NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md`;
`NOTEBOOK_ENTRIES/d1b_premise_fails_all_five_arms_collapse_20260803T2030Z.md`;
`NOTEBOOK_ENTRIES/d1_programme_free_collapsing_in_training_20260803T1930Z.md`;
`NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md`.*

**But the boundary is much lower than "low rank".** At R1 effective rank **9.11 and 9.14** of a nominal
256 dimensions — 3.6% of ambient, a number any rank-monitoring practice would read as severe collapse —
the two D2 seed-44 arms still read held-out channels of **0.5983 and 0.4757** against a permutation
null of **0.140** (§4.3). A representation at 3.6% of its nominal rank was carrying a large,
permutation-significant molecular channel.

So the defensible statement is not "low rank means little information". It is:

> Effective rank near its floor (≈ 1–2, with patient-to-patient mutual cosine ≈ 1) is reliable
> evidence of total collapse. Anywhere above that — including at 3.6% of nominal dimensionality — it
> is uninformative about the channel, in both directions, and non-monotone with respect to it.

**And even for the surviving use, rank is not the best instrument.** Patient-to-patient mutual cosine
reaches 0.977–0.999 in every collapsed arm above, is a single matrix product, needs no SVD, and has a
natural scale with a meaningful maximum. In every instance here where rank correctly signalled
collapse, mutual cosine signalled it too and more legibly. We therefore recommend mutual cosine as the
primary collapse diagnostic, with rank as secondary confirmation.

**What we have not shown.** We have not found a case of total collapse that effective rank missed. Our
one-sided claim is therefore asymmetric on the evidence as well as in its statement, and §4.9 should
be read as "we did not falsify the collapse-diagnostic use", not as "we verified it".

---

## 5. Limitations

### 5.1 Novelty, stated at its narrowest

The core negative — that effective rank does not predict downstream quality — **is already published**
(§2.2). Thilak et al. report RankMe Spearman 0.3174 and state that "a high rank does not guarantee
superior performance"; a RankMe co-author has co-signed the claim that RankMe "fail[s] to adequately
evaluate representation quality". This paper's claim to novelty is confined to (i) the in-scope
within-method failure of §4.3, (ii) the reproducibility floor of §4.7, (iii) the dose–response
magnitude miscalibration of §4.2, and (iv) the domain.

The prior-art census behind §2.2 is **incomplete**: the natural-language search path was unavailable
and no citation-graph sweep of RankMe's citing papers has been run. **If a paper reporting an
in-scope, within-method RankMe failure exists, contribution (i) collapses to a replication and this
paper must be reframed around (ii).** Closing that census is the first pre-submission task.

### 5.2 What is not measured

| would-be measurement | why it is absent |
|---|---|
| **LiDAR, or any published alternative criterion, computed on our artifacts** | Identified during this draft's prior-art sweep; requires a labelled probe on every artifact; not run. **The most valuable missing measurement in this paper** (§2.5). |
| A controlled repeat design for §4.7 (N retrainings of one configuration, rank and channel on each) | Not run. §4.7 rests on three seeds plus one accidental retrain, and cannot distinguish rank-specific variance from stack non-determinism, architecture or schedule. |
| An error bar on any dilution rank or channel value | Single seed, single donor draw; bootstrapping donor assignments is CPU-only and unrun. |
| An equivalence test on instance 2's channel difference | The paired bootstrap the source says "is still required" was never run. §4.4 states that "unchanged" means the point estimates differ by 0.002 and nothing more. |
| Rank and channel measured with **one** statistic across all instances | Three implementations exist (§3.1); the historical numbers cannot be recomputed without re-running the original code generations, and for instance 1 the artifacts do not exist. |
| Any instance on a **second architecture, cohort or modality pair** | All six are one architecture family (transformer aggregator over frozen H-Optimus-0 patch tokens with a biology head), one cohort (TCGA), one modality pair (morphology → bulk expression). `claim_guards.no_external_cohort` is undischarged for every morphology result on this project. |
| **E1**, the preregistered rank-versus-information experiment | Built (`v2/calibra/e1_rank_information.py`, `aggregate_e1.py`, equivalence margin 0.10, three-seed requirement) and **never run**. It is the experiment this paper should have been built on. |
| A case where the collapse diagnostic **fails** | We have not found one (§4.9). |
| Rank computed with RankMe's ε convention | Our implementation filters at `> 1e-12` instead (§2.6). No number here is comparable to a published RankMe value. |

### 5.3 Instance-specific limitations

Stated in each instance's own section and gathered here: §4.2 — censored floors, single seed and donor
draw, one pooling scheme, "lower bound" label withdrawn by the source. §4.3 — negative control small
but non-zero and same-signed, so 10–20% of the gap may be generic quality; point estimates not
reproducible from the seed; three seeds is a small number and the finding is unreliability, not a
quantified error rate. §4.4 — G0.4 not discharged, single seed, no interval, both differences inside
the stack's retraining noise. §4.5 — primary source file does not exist; benchmark statistic
undefined; rank statistic unknown. §4.6 — hard rank not effective rank, structural ceiling of 16,
train batch of 16, and a same-objective effective-rank measurement that cuts the other way. §4.7 — two
configurations, one stack, no controlled repeat design. §4.8 — pending, with a non-reproducible
admission gate.

### 5.4 A defect in our own evidence base, reported not repaired

Instance 1's cited source (`paper/.../RESULTS.md`) does not exist in this repository. We have not
attempted to reconstruct it. The number is retained in §4.5 as history and is excluded from every count
and summary statement. Any future citation of "+107% rank at flat benchmark" from this project must
carry this paragraph.

### 5.5 The claim is about a scalar, not about geometry

Nothing here says representation geometry is uninformative. It says one scalar summary of the spectrum
is, in one regime. Alignment, uniformity, the full singular-value profile, LiDAR's statistic,
per-feature spread and mutual cosine may each behave differently, and §4.9 recommends one of them over
rank. Per-feature spread is *also* misleading in at least one place: `programme_free` at epoch 21 has
**higher** per-feature std than `programme_only` (0.0137 against 0.0044) and **lower** rank, because
the collapse is to the family `zᵢ = m + aᵢ·u` rather than to a point
(`NOTEBOOK_ENTRIES/d1_programme_free_collapsing_in_training_20260803T1930Z.md`). We report that as a
second scalar failing, not as a recommendation.

### 5.6 Scope of the negative

This paper does not claim that anti-collapse regularisation is useless, that rank should never be
computed, that RankMe is wrong about what it claims, or that the published methods in §2.3 are wrong.
It claims that within RankMe's own reserved regime the criterion misordered two of three seeds on our
data; that its run-to-run variance on this stack exceeds the differences it would be used to select
on; and that the informal reading of rank as a representation-health indicator is contradicted by a
dose–response in which rank under-reports information loss by 3.7×.

---

## 6. Conclusion

Effective rank is cheap, label-free and never fails to return a number. Those three properties are why
it is used to select between representations, and they are also why it is silently unreliable in that
use. That it can be unreliable is known: Thilak et al. reported it in 2023, and a RankMe co-author has
since agreed.

What we add is a failure inside the boundary RankMe draws around itself. In a within-method,
single-architecture, in-distribution comparison with arms matched by construction and three seeds with
paired bootstraps, a rank-based selection rule picks the right arm once, the wrong arm once, and cannot
choose once — while the performance ordering is stable in all three. And independently of whether rank
is right on average, on this stack it moves 2.7× when one configuration is retrained with the same
seed and spans 3.7× across three seeds, against a paired difference stable to ±0.012. A criterion that
cannot resolve its own re-measurement cannot select between configurations, and we propose the check
that would have caught this: retrain one configuration, measure the rank spread, and require the
between-configuration difference to exceed it.

We report the evidence that constrains this at equal prominence. Centred effective rank does fall to
~1 under total collapse, so the collapse-diagnostic use survives — though a patient-to-patient cosine
is cheaper and more legible for the same job. The instance this project previously called its most
dramatic turns out to be a *hard* matrix rank at a structural ceiling, and we withdraw that
description. And a representation at 3.6% of its nominal dimensionality still carried a large,
permutation-significant channel, so "low rank means little information" is also false outside total
collapse.

Finally, we report that three mutually incompatible statistics are implemented under the name
`effective_rank` in the repository that produced these results, and that our own historical instances
were not all measured with the same one. We found this while writing this paper. It is the best single
illustration of the thesis we could have asked for, and we did not design it.

---

## Appendix A — provenance index

| § | claim | source file(s) | box path |
|---|---|---|---|
| 2.1 | Roy & Vetterli definition; RankMe claims and restrictions | full-text PDFs, verified 2026-08-04 | — |
| 2.2 | LiDAR quotes and Spearman/Kendall values | full-text PDF, arXiv:2312.04000v1 | — |
| 3.1 | R1 definition | `v2/calibra/spectral.py:14-29` | — |
| 3.1 | R2 definition | `v2/research/rebase/d1_audit.py:149-153` | — |
| 3.1 | R3 definition | `v2/research/rebase/d1_geometry_probe.py:50-53` | — |
| 3.1 | E1 built, never run | `v2/calibra/e1_rank_information.py`, `aggregate_e1.py:38`; absence confirmed against `GATE_LOG.md` and `runs/` | — |
| 3.2 | channel statistic and null | `v2/calibra/spectral.py:78-108`, `v2/calibra/run_calibra.py` | — |
| 3.4 | D2 pair-manifest hashes | `D2_PAIR_MANIFEST.json`; `NOTEBOOK.md` 2026-08-01 20:35 UTC | `~/e0_run/d2_v3/` |
| 3.5, 4.7 | seed non-reproducibility | `D2_RESULT.md` §4 | `~/e0_run/d2_v3/` |
| 4.2 | dilution table | `DILUTION_LOWER_BOUND.md` §2, §6; `NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md` | `p1_evidence/dilution/` |
| 4.3 | D2 table | `D2_RESULT.md` §2, §4 | `~/e0_run/d2_v3/bootstrap/`, `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json` |
| 4.4 | Phase 1b table | `PHASE1B_TARGETED_READOUT.md` §3, §5, §7 | `runs/calibra_v3_targeted` |
| 4.5 | +107% instance | `ENGINE_CLD.md` §1; `HANDOFF_BUILD_AGENT.md` §1–2 | **no artifact; cited source does not exist** |
| 4.5 | decorrelation minimum is collapse | `NOTEBOOK.md` 2026-08-03 decision; `NOTEBOOK_ENTRIES/g26_term_isolation_20260803T0930Z.md` | `~/e0_run/d1_diag/` |
| 4.6 | three-arm collapse table | `NOTEBOOK.md` 2026-08-02 01:20 UTC; `NOTEBOOK_ENTRIES/g26_variance_floor_fix_20260803T0210Z.md` | `scratchpad/collapse_diag.py` on the A100 |
| 4.6 | counter-measurement 12.88 → 1.00 | `NOTEBOOK_ENTRIES/g26_rank_collapse_diagnosis_20260803T0500Z.md` | `~/e0_run/d1_diag/` |
| 4.7 | short-horizon probe repeats | `NOTEBOOK_ENTRIES/rank_probe_repeat_variance_20260804T0900Z.md` | `~/ws_d1/probevar_*.log` |
| 4.8 | D1-A epoch-39 geometry | `NOTEBOOK_ENTRIES/d1a_control_complete_and_gate_fails_2of3_in_runner_20260804T0100Z.md` | `~/e0_run/d1_v1/` |
| 4.8 | gate non-reproducibility, 8 runs | `NOTEBOOK_ENTRIES/g26_is_not_reproducible_20260804T0700Z.md` | `~/ws_d1/gatevar_*.log` |
| 4.9 | collapse-regime rank values | `g26_rank_collapse_diagnosis`, `d1b_premise_fails_all_five_arms_collapse`, `d1_programme_free_collapsing_in_training`, `d1a_control_complete_and_gate_fails_2of3_in_runner` | `~/e0_run/d1_diag/`, `~/e0_run/d1_v1/` |

## Appendix B — code index

| function | file:line | used by |
|---|---|---|
| `effective_rank` (R1, = Roy & Vetterli) | `v2/calibra/spectral.py:14` | §4.2, §4.3, §4.4, §4.7 |
| `effective_rank` (R2) | `v2/research/rebase/d1_audit.py:149` | §4.8 pending table |
| effective rank (R3, inline) | `v2/research/rebase/d1_geometry_probe.py:50` | §4.8.2, §4.9 |
| `heldout_top_cca` | `v2/calibra/spectral.py:78` | every channel number |
| `heldout_single_direction_correlation` | `v2/calibra/spectral.py:111` | detection-floor-scale controls |
| `cross_fitted_residuals`, `confound_design` | `v2/calibra/residualise.py` | all adjusted readouts |
| `permutation_null` | `v2/calibra/calibration.py` | all nulls in §3.2 |
| paired bootstrap | `v2/paired_bootstrap.py` | §4.3, §4.8 |
| `symmetric_infonce` | `v2/losses.py:13` | §4.6 |
| `stable_rank` | `v2/calibra/e1_rank_information.py` | not used in this paper |

## Appendix C — the caveat that must travel with each instance

Reproduced verbatim so that any future quotation of a number can carry it.

- **Instance 6 (§4.3).** *"REPORTED, NOT INTERPRETED (blocker 5). WSI states are ~0.80 collinear at
  init, so a narrower rank may reflect resistance to an already-collapsed view rather than dictionary
  content."* (`v2/research/rebase/d1_audit.py`, A5 note.)
- **Instance 4 (§4.2).** *"Single seed (42) and a single draw of donor assignments… gives no error bar
  on the level-to-level differences."* And: *"The whole curve is one representation."*
- **Instance 2 (§4.4).** *"`full` vs `programme_only` manifests were not verified as matched on
  epochs/LR/budget in this run (G0.4). Until they are, the rank comparison in §5 is suggestive, not
  causal."*
- **Instance 3 (§4.6).** The column is labelled *"`z_biology` matrix rank"* in its source. It is a hard
  numerical rank whose maximum is the batch size of 16, and the centred effective rank of the same
  objective falls 12.88 → 1.00. **The project's earlier description of this as a strong effective-rank
  instance is withdrawn.**
- **Instance 1 (§4.5).** Earlier codebase generation; different benchmark statistic, undefined in this
  repository; rank statistic unknown; **cited source file does not exist**.
- **Instance 5 (§4.8).** Pending. D1-A's 9.81 / 1.71 numbers are a control documenting a defect:
  *"Nothing about programme supervision may be concluded from it — the contrastive arm never trained,
  so the comparison measures a defect, not an ablation."*
- **All instances.** RankMe restricts itself to same-method comparisons and to a necessary-not-
  sufficient reading. Four of the six instances do not contradict it as written, and each is labelled
  in §4.1.
