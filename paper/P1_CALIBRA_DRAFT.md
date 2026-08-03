# What would this analysis have missed? Injection-certified transmission and detection floors for confound-adjusted morphology–molecular analyses

**CALIBRA — a calibration instrument for cross-modal claims.**

*Submission draft, 2026-08-03. Every number in this document traces to a named artifact, notebook
entry or evidence file in this repository; each table carries a `provenance` line. Numbers that do
not exist are marked as not measured rather than estimated. Citations that have not been verified
against a live bibliographic source in this drafting pass are marked `[UNVERIFIED]`; missing
citations are marked `[CITATION NEEDED]`. Nothing in either category may enter a submitted
manuscript without resolution — see §2.7.*

---

## Abstract

Computational pathology routinely reports that tissue morphology predicts molecular state. Because
tissue source site, scanner and cancer type confound these analyses, the standard remedy is to
residualise both modalities against a shared covariate design before measuring association. Such
analyses report a point estimate and sometimes a permutation *p*, but almost never the complementary
quantity: the effect size the analysis would have failed to detect. Without it a small or absent
association cannot be distinguished from an analysis that was never sensitive enough to see one.

We describe CALIBRA, an instrument that injects a correlation of known strength onto a named
direction pair spanning the image and molecular blocks, pushes it through the *identical* pipeline
including the residualisation, and reports what came back. On 2,530–6,427 TCGA patients with frozen
H-Optimus-0 whole-slide features and bulk expression targets we report four things. (i) The
adjustment works and is verified rather than assumed, twice and on two different confounders: raw
cancer-type balanced accuracy falls 0.463 → 0.035 against a chance rate of 0.048, and joint
tissue-source-site accuracy falls 0.3633 → 0.0118 against a chance rate of 0.0118, a 21–45× drop with
zero breaching axes in six representation states. (ii) It does not destroy the signal it is meant to
leave alone: injected-signal attenuation is 0.974–1.039 across six states and 1.07–1.12 across twelve
design × sample-size cells. (iii) The floors themselves — a paired *transmission* floor censored at
≤ 0.01, and an unpaired *detection* floor of 0.25–0.40 that is the only quotable detection limit —
with a permutation *p* of 1/2001 for every state. (iv) The detection floor is set by a
residualisation-induced correlation floor, not by sample size: it is pinned at 0.25–0.30 from
n = 1,000 to n = 6,427 while a structureless design of comparable rank falls 0.050 → 0.010.

The induced correlation is not a new phenomenon. It is the Yule (1907) / Frisch–Waugh–Lovell partial
correlation identity, which our pipeline reproduces draw-for-draw to 8.6 × 10⁻¹⁶ over 270 cells, and
it has been published as a warning in fMRI, GWAS and genomics. What we add is a magnitude in
correlation units for a cross-modal, cross-fitted, correctly-specified adjustment — 0.0748 for the
anchor design at n = 2,530 against 0.0037 for its own row-permuted twin and 0.0035 for a Gaussian
design of matched width — together with the demonstration that it is structural rather than the
classical degrees-of-freedom term, that it is invariant to design rank over a 14.5× range, and that
it does not shrink with n.

We report the instrument's failures as fully as its successes: a must-fail site control that did not
fail on raw representations, a per-axis certificate that would have issued a clean bill of health on
a real leak, a predeclared estimator-robustness prediction that failed at 55.7% against a 25% bar,
and the finding that 76–82% of the per-target channel is reproduced by covariate-matched random gene
sets, invariant to contamination. No external cohort has been through the instrument; every number
here is TCGA.

### Short abstract (~200 words, for a length-capped methods venue)

Confound-adjusted cross-modal analyses in computational pathology report a point estimate and
sometimes a permutation *p*, but not the effect size they would have failed to detect. Without it a
small association cannot be distinguished from an insensitive analysis. We inject a correlation of
known strength onto a named direction pair spanning image and molecular blocks and push it through
the identical pipeline, including a cross-fitted residualisation against cancer type and tissue source
site. On 2,530–6,427 TCGA patients the adjustment removes what it claims to (joint site accuracy
0.3633 → 0.0118 against chance 0.0118, zero breaching axes in six states) and costs the signal
essentially nothing (attenuation 0.974–1.039), but cannot detect a single-direction correlation below
0.25–0.40. That floor does not improve from n = 1,000 to n = 6,427, because it is set by the
residualisation-induced correlation — the Yule/Frisch–Waugh–Lovell partial correlation, reproduced by
our pipeline to 8.6 × 10⁻¹⁶ — which is 20–35× larger for a real confound design than for a
matched-rank structureless one and does not shrink with n. We report a must-fail control that did not
fail, a per-axis certificate that would have passed a real site leak, a failed predeclared robustness
prediction, and the finding that covariate-matched random gene sets reach 76–82% of real ones. No
external cohort.

---

## 1. Introduction

### 1.1 A literature of point estimates without sensitivity statements

The claim that tumour morphology carries molecular information is now routine. Deep models predict
point mutations, microsatellite instability and molecular subtype directly from haematoxylin-and-eosin
slides (Kather et al., *Nature Cancer* 2020) `[UNVERIFIED]`; slide-level foundation models report
molecular-profile and outcome prediction across cohorts (CHIEF, Wang et al., *Nature* 2024;
Prov-GigaPath, Xu et al., *Nature* 2024) `[UNVERIFIED]`; and encoders are pretrained against paired
genomic and transcriptomic profiles so that molecular meaning is internalised into the image
representation (THREADS, Vaidya et al., arXiv:2501.16652, 2025). Benchmarks have followed, most
prominently HEST-1k (Jaume et al., NeurIPS 2024, arXiv:2406.16192).

The reporting convention across this literature is a point estimate of association, sometimes with a
bootstrap interval, occasionally with a permutation null. What is almost never reported is the answer
to the complementary question: *what effect size would this analysis have failed to detect?* Without
that number a null result is uninterpretable. A study that residualises confounds, measures a small
cross-modal correlation and concludes that morphology carries little information about a given
programme cannot distinguish "there is no signal here" from "this pipeline could not have recovered a
signal of this size even if it were present". The two conclusions have opposite scientific
consequences and the same published appearance.

This is not a hypothetical failure mode in our own work. Three separate theses on this project were
retired by adversarial review on exactly this objection — that a confound adjustment might have
destroyed the effect it was meant to isolate — and none could be defended, because no sensitivity
number existed to defend them with.

### 1.2 Confounding, and why the standard remedy has a null of its own

The confounding in this setting is documented and substantial. Howard et al. (*Nature Communications*
12:4423, 2021) showed that the TCGA tissue-submitting site is detectable from the slide image,
survives colour normalisation, and correlates with clinical and genomic labels. Subsequent work found
that the problem is not solved by scale: hospital signatures persist in foundation-model embeddings
(Kömen et al., arXiv:2411.05489, 2024) `[UNVERIFIED]`; across ten public pathology foundation models
medical-centre identity is encoded more strongly than tissue or cancer type (de Jong et al.,
arXiv:2501.18055, 2025); representational-similarity analysis finds pronounced slide-dependence and
weak disease-dependence (Mishra & Lotter, arXiv:2509.15482, 2025); scanner hardware is a separate
axis again (Carloni et al., arXiv:2507.22092, 2025) `[UNVERIFIED]`.

We reproduce this on our own representations rather than citing it. A joint linear discriminant over
all 256 axes of a raw representation state recovers the pooled tissue source site at balanced
accuracy 0.2348–0.3633 against a chance rate of 0.0118, with *p* = 1/1001 in every one of six
state × artifact combinations (§4.2).

The standard remedy is to remove the covariate — correct features post hoc (ComBat on deep features,
Murchan et al., *J. Pathol. Inform.* 15:100396, 2024), or residualise both modalities against a design
matrix before measuring association. This is the right instinct, and it has a consequence that is
well known in several other literatures and almost never carried into this one. Residualisation is a
projection. Projecting two signals onto the orthogonal complement of a shared design does not leave
their relationship untouched: for signals constructed exactly orthogonal in-sample, the residual
correlation is precisely the multivariate partial correlation

```
corr(Mu, Mv) = − R_u · R_v · cos θ / sqrt((1 − R_u²)(1 − R_v²)),    M = I − X(XᵀX)⁻¹Xᵀ
```

a direct consequence of the Yule–Frisch–Waugh–Lovell identity (Yule, *Proc. R. Soc. Lond. A* 79:182,
1907; Frisch & Waugh, *Econometrica* 1(4):387, 1933; Lovell, *JASA* 58(304):993, 1963). The practical
consequence is that the null hypothesis of "no cross-modal association after adjustment" is not
"correlation zero"; it is a non-zero, design-dependent, draw-dependent quantity that the pathology
literature does not compute. **We claim no discovery here.** The phenomenon is textbook and it has
been published as an applied warning at least five times (§2.4). What did not exist, and what we
supply, is its magnitude in correlation units for a cross-modal, cross-fitted, correctly-specified
adjustment at TCGA scale, together with evidence that it is structural rather than the classical
degrees-of-freedom term.

### 1.3 What this paper is, and what it is not

This is a methods and instrument paper. Its deliverable is an audit procedure and a set of measured
floors, not a finding about tumours. **No biological claim is made anywhere in this manuscript.** The
representation states audited here exist for other work on the same project; here they are specimens,
and the strongest results in the paper are about measurement, not about biology.

Two consequences follow, and both are load-bearing.

First, the paper reports every control that came back the wrong way with the same prominence as the
ones that passed, because a battery that only ever passes proves nothing. The two most valuable
outputs of the negative-control battery are a must-fail control that did not fail (§4.2) and a
must-fail control that fails on one statistic and reaches 76–82% of the real signal on another
(§4.8.2).

Second, the paper is careful about the difference between an instrument that is *demonstrated* and
one that is *validated*. No external cohort has been through CALIBRA. Every number here is TCGA, with
its documented site and scanner effects. The `no_external_cohort` blocker in
`v2/calibra/claim_guards.py` is undischarged for every morphology result on this project, and it is
enforced in code: `validate_claim` refuses to mark a `legible_axis` or `gene_attribution` claim
publishable, and no such claim appears here.

### 1.4 Contributions

In descending order of how well evidenced they are.

1. **A confound adjustment verified rather than asserted, on two different confounders, with the
   verification reported as a certificate.** Cancer-type balanced accuracy 0.463 → 0.035 (chance
   0.048); joint tissue-source-site balanced accuracy 0.3633 → 0.0118 (chance 0.0118), a 21–45× drop
   with zero breaching axes across six states.
2. **A demonstration that the same adjustment does not destroy signal.** Injected-signal attenuation
   0.974–1.039 across six states, 1.07–1.12 across twelve design × n cells, and 0.944–1.228 in the
   original 40-draw sweep — i.e. ≈ 1 in every measurement we have. This is the number that answers
   the objection that retired three earlier theses on this project.
3. **Injection-certified floors, and the requirement that they be reported as a pair.** A paired
   transmission floor (censored at ≤ 0.01) that answers "does the pipeline transmit a signal of this
   size?", and an unpaired detection floor (0.25–0.40 here) that is the only quantity quotable as a
   detection limit. Quoting the paired floor as a detection limit understates the threshold by more
   than an order of magnitude.
4. **The measured magnitude of the residualisation-induced correlation floor in a cross-modal
   setting, with the phenomenon explicitly attributed to Yule/FWL and to the five literatures that
   have already warned about it.** 0.0748 at the anchor design and n = 2,530, against 0.0037 for the
   row-permuted twin of the same design and 0.0035 for a Gaussian design of matched width — 20.4×,
   widening to 35.4× at n = 6,427. Magnitude only; nothing broader is claimed.
5. **The operational consequence: more patients do not buy sensitivity.** For every real design the
   detection floor is pinned at 0.25–0.30 from n = 1,000 to n = 6,427, while a structureless design of
   comparable rank falls 0.050 → 0.010. Cohort sizing for this class of analysis must be written
   against the induced-correlation predictor, not against a power calculation in n.
6. **A finding about certification rules, not about data: a per-axis certificate is not sufficient.**
   The site leak we found is smeared, not concentrated — best single axis 0.055, median axis *below*
   its own permutation null — while the joint discriminant reaches 0.235–0.363. A per-axis-only
   screen, which is literally what our own specification asked for, would have passed this
   representation. The joint test must be a required field of the certificate schema.
7. **An executable admissibility layer and a ledger that separates health gates from observations,**
   with the separation enforced by a test that plants a baseline beating us 0.95 to 0.40 and asserts
   the run verdict does not move.
8. **A recorded instrument failure as the paper's motivating example:** a spike readout that passed
   all eleven synthetic self-tests while returning undefined detection floors on every real dataset,
   caused by three nested defects.

---

## 2. Related work

### 2.1 Morphology-to-molecular prediction

Tile-level self-supervised encoders (UNI, Chen et al., *Nature Medicine* 2024; Virchow, Vorontsov et
al., arXiv:2309.07778, 2023; Phikon-v2, Filiot et al., arXiv:2409.09173, 2024; H-Optimus-0, Bioptimus
open weights release, 2024) and slide-level encoders (Prov-GigaPath, Xu et al., *Nature* 2024; TITAN,
Ding et al., arXiv:2411.19666, 2024) have made frozen pathology features a commodity input. The
morphology-to-molecular claim proper begins with Kather et al. (*Nature Cancer* 2020) and continues
through CHIEF (Wang et al., *Nature* 2024) and THREADS (Vaidya et al., arXiv:2501.16652, 2025).
HEST-1k (Jaume et al., NeurIPS 2024) is the reference benchmark for expression prediction from
morphology. Tizhoosh (arXiv:2510.23807, 2025) argues the failures in this area are structural rather
than a matter of scale.

None of this work reports what effect size its analysis would have missed. Our contribution is
orthogonal to all of it: we do not propose an encoder or a benchmark, but an audit that any of them
could be run through.

### 2.2 Confounding and site effects in TCGA-based studies

Howard et al. (*Nature Communications* 12:4423, 2021, DOI 10.1038/s41467-021-24698-1) established
that submitting-site signatures are learnable from H&E, survive stain normalisation, and inflate
apparent accuracy on survival, mutation and stage; they propose preserved-site cross-validation.
Notably, Howard et al. also run a synthetic manipulation — varying ER negativity of target slides
from 0 to 100% and applying an artificial staining artifact to 0–100% of them — but to *demonstrate*
confounding, not to certify a detection floor after adjustment. It is the nearest thing in this
domain to what we do. Schmitt et al. (*JMIR* 2021) `[UNVERIFIED]` broadened the confounder list to
scanner type, institution and preparation date. Dawood et al. ("Buyer Beware", bioRxiv
2024.06.23.600257) show that per-biomarker models predict a correlated bundle rather than an isolated
biomarker.

Mitigations divide into image-space normalisation, feature-space correction (Murchan et al.,
*J. Pathol. Inform.* 15:100396, 2024, DOI 10.1016/j.jpi.2024.100396), adversarial removal, and
covariate residualisation. Murchan et al. is the closest published instance of the exact analysis
this paper audits — ComBat-harmonising WSI deep features against tissue source site and *then*
predicting molecular features — with no consideration of an induced floor and no sensitivity
statement.

### 2.3 Injection and spike-in certification

Spike-in certification of detection limits is long established **within a single modality**. The ERCC
synthetic spike-ins yield sensitivity limits for RNA-seq (Jiang et al., *Genome Research* 21(9):1543,
2011, DOI 10.1101/gr.121095.111), and Munro et al. (*Nature Communications* 5:5125, 2014,
DOI 10.1038/ncomms6125) define a named **limit of detection of ratio (LODR)**, estimated across 12
laboratories and 3 measurement processes. Injecting known-strength signal into *real* data expressly
to benchmark confound-adjustment methods is established for differential expression: Gerard's
`seqgendiff` (*BMC Bioinformatics* 21:206, 2020, DOI 10.1186/s12859-020-3450-9) adds a known amount of
signal to real RNA-seq counts and lets the user control the correlation between observed covariates
and unobserved surrogate variables, specifically to benchmark SVA/RUV/CATE.

The structural precedent for our transmission floor is the **injection–recovery completeness curve**
used for instrument characterisation outside biology — LIGO hardware injections (Biwer et al.,
*Phys. Rev. D* 95:062002, 2017, DOI 10.1103/PhysRevD.95.062002), Kepler DR25 flux-level transit
injection tests, eROSITA source-detection simulations. Naming the connection costs nothing and
inoculates the paper against the most likely methods-venue objection.

**What we may therefore claim, and no more:** adapting that logic, we inject a correlation of known
strength onto a named direction pair spanning **two modalities** and push it through the identical
cross-fitted residualisation, reporting both the fraction transmitted and the smallest true
correlation reliably recovered. To our knowledge this paired transmission/detection reporting has not
previously been applied to a confound-adjusted cross-modal biological analysis. The load-bearing
qualifiers — *"within a single modality" → "two modalities"*, *"to our knowledge"*, *"has not
previously been applied to"* rather than *"is the first"* — are deliberate. Two independent
literature sweeps (a spike-in/calibration sweep and a neuroimaging sweep) separately failed to find
any paper that injects a signal of known strength and measures how much survives nuisance regression
specifically; a Europe PMC query for `("simulated signal" OR "synthetic signal" OR "known ground
truth") AND ("nuisance regression" OR "confound") AND "recover"` returned **0 hits**. That is absence
from two sweeps, not proof of absence, and it is stated that way.

*Provenance: `v2/research/rebase/nature/NOVELTY_SEARCH.md` §Q2. Verdict recorded there: PARTIALLY
ANTICIPATED.*

### 2.4 Residualisation-induced correlation — prior art, and the narrow claim that survives

**This phenomenon is not novel and we do not claim it.** An adversarial prior-art sweep returned the
verdict **ALREADY REPORTED**, in at least five literatures:

| literature | reference | what it already says |
|---|---|---|
| statistics / econometrics | Yule 1907; Frisch & Waugh 1933; Lovell 1963; history in Basu, arXiv:2307.00369 | residualising two variables on a shared regressor set yields exactly their partial correlation; a zero raw correlation maps to a non-zero partial correlation whenever both load on the design |
| fMRI | Murphy et al., *NeuroImage* 44(3):893, 2009, DOI 10.1016/j.neuroimage.2008.09.036; Saad et al., *Brain Connect.* 2(1):25, 2012; Fox et al., *J. Neurophysiol.* 101(6):3270, 2009; Murphy & Fox, *NeuroImage* 154:169, 2017 | global-signal regression mathematically *mandates* negative correlations; "calls into question the interpretation of negatively correlated regions" |
| cross-modal neuroimaging | **Winkler et al., *NeuroImage* 220:117065, 2020**, DOI 10.1016/j.neuroimage.2020.117065 | CCA on imaging × non-imaging residuals: "residualisation introduces dependencies among the observations that violate the exchangeability assumption", leading to inflated error rates — plus a fix |
| cross-modal neuroimaging | Alfaro-Almagro et al., *NeuroImage* 224:117002, 2021, DOI 10.1016/j.neuroimage.2020.117002 | "spurious associations can be induced between pairs of otherwise independent variables **if the unconfounding is not carried out correctly**" |
| genomics | Li et al., *Biostatistics* 24(3):635, 2023, DOI 10.1093/biostatistics/kxab039; Nygaard et al., *Biostatistics* 17(1):29, 2016, DOI 10.1093/biostatistics/kxv027; Zindler et al., *BMC Bioinformatics* 21:271, 2020 | two-step batch correction "introduces a correlation structure into the adjusted data", derived as `Σ = (I − X(XᵀX)⁻¹Xᵀ)σ²`; increasing the number of corrected factors gives an exponential increase in false positives |
| GWAS | Aschard et al., *AJHG* 96(2):329, 2015; Dahl et al., *Genetics* 211(4):1179, 2019, DOI 10.1534/genetics.118.301768 | "conditioning on genomic PCs can cause, rather than remove, bias"; "induces spurious correlation" |
| folklore | Gordon Smyth, Bioconductor support #133791, 2020, on `limma::removeBatchEffect` | "the batch correction has introduced correlations" |

*Provenance: `v2/research/rebase/nature/NOVELTY_SEARCH.md` §Q1 (Google Scholar / Europe PMC REST /
arXiv / PubMed-PMC / publisher pages / Bioconductor archive; 200-call web budget exhausted).*

Two distinctions survive that sweep, and every claim in this paper is phrased around them.

**(a) Correctness.** Alfaro-Almagro et al. condition their statement on *incorrect* unconfounding
("if the confounds were not demeaned first"). Our statement is that *correctly performed,
cross-fitted* residualisation of *exactly orthogonal* signals still yields a non-zero residual
correlation, as a geometric fact. This is the sharpest surviving distinction.

**(b) Units.** Winkler et al. 2020 is the most threatening prior work — it is cross-modal, it warns,
and it supplies a fix. We retrieved the full text from PMC7573815 and audited it in five targeted
passes over Theory §2.6, Simulations, Results §4.3/4.5, Discussion §5.3, every figure caption and
every table header. **It does not report an induced-correlation magnitude.** Everything it reports is
in error-rate, power or p-value units: Table 6 per-comparison error rate (%), Table 7 pcer (simple
residualisation 83.85% [82.17–85.40] vs ~5% for Huh–Jhun/Theil), Table 8 pcer, Table 9 FWER against
N = 100…1000, Table 10 power. The closest statement is a stochastic-ordering claim with no magnitude:
"the sample canonical correlations in the unpermuted case are *stochastically larger* than in the
permuted". The only place the design rank R enters a formula is the Bartlett/Wilks degrees-of-freedom
correction, used to *undo* the effect, never inverted into a magnitude.

*Provenance: `v2/research/rebase/nature/TRACK2_INDUCED_CORRELATION.md` §0 (full text from PMC7573815);
the check was recorded as settled in `P1_PREDECLARATION.md` §C before the Track 2 sweep was run. This
closes residual coverage gap #5 of `NOVELTY_SEARCH.md`, which had flagged that if Winkler's tables
contained a magnitude, our numeric residue would disappear entirely.*

We must therefore cite, so that no reviewer can present them as omissions: (i) Yule/FWL for the
identity; (ii) Winkler et al. 2020 for the inferential consequence and the remedy; (iii) **Muirhead
(1982), *Aspects of Multivariate Statistical Theory* / Anderson (2003), *An Introduction to
Multivariate Statistical Analysis***, for the classical result that rank-R residualisation leaves an
effective sample size of N − R, which Winkler's own correction uses `[UNVERIFIED — edition, chapter
and page not checked in this pass]`.

We must also concede, in the text and before a reviewer says it, that the identity is one line of
algebra (§1.2), and that empirical cross-modal noise floors are an established sub-industry with a
*different* mechanism — spatial autocorrelation rather than a confound design (Markello & Misic,
*NeuroImage* 236:118052, 2021; Burt et al., *NeuroImage* 220:117038, 2020; Fulcher et al.,
*Nat. Commun.* 12, 2021, DOI 10.1038/s41467-021-22862-1). Our floor is *additional and orthogonal* to
that family, not a predecessor of it. Finally we cite Marek et al. (*Nature* 603:654, 2022,
DOI 10.1038/s41586-022-04492-9) ourselves rather than letting a reviewer introduce it: a
confound-induced floor of 0.07–0.10 would *swamp* the median |r| = 0.02–0.03 typical of large
cross-modal association studies.

### 2.5 Negative-control culture in genomics

Genomics has a mature negative-control culture that computational pathology has partly not imported.
Venet, Dhanasekaran & Sotiriou (*PLoS Comput. Biol.* 2011) showed that most random gene-expression
signatures are significantly associated with breast-cancer outcome, establishing the random-signature
null as a required control `[UNVERIFIED]`. Random and size-matched gene-set nulls are codified in
standard tooling. On the perturbation side, simple baselines beat deep models under held-out-
perturbation splits (Wong, Hill & Moccia, *Bioinformatics* 41(6):btaf317, 2025). Covariate adjustment
itself has a long literature — surrogate variable analysis (Leek & Storey 2007), ComBat (Johnson, Li
& Rabinowitz 2007), removal of unwanted variation (Gagnon-Bartsch & Speed, *Biostatistics*
13(3):539, 2012) — evaluated on whether the unwanted factor is removed and the factor of interest
preserved `[UNVERIFIED for Leek & Storey and Johnson et al. — author initials and year not
re-checked]`. Note that RUV-type methods use controls known **not** to change, a *negative*-control
lineage; our positive-injection construction is a genuine and worth-stating distinction. Simulation-
based calibration (Talts et al., arXiv:1804.06788; Cook, Gelman & Rubin, *JCGS* 15(3):675, 2006)
measures uniformity of rank statistics with no injection into real data and no detection floor; the
distinction from CALIBRA is real.

### 2.6 Representation-geometry metrics

A parallel line proposes geometric proxies for representation quality: alignment and uniformity on
the hypersphere (Wang & Isola, ICML 2020), dimensional collapse in contrastive learning (Jing et al.,
ICLR 2022), the variance and covariance terms of VICReg (Bardes, Ponce & LeCun, arXiv:2105.04906),
and effective rank as a scalar summary of a spectrum (Roy & Vetterli, 2007). In pathology
specifically, the Robustness Index (de Jong et al., 2025) and representational-similarity analysis
(Mishra & Lotter, 2025) are proposed as confound-aware representation diagnostics.
`[CITATION NEEDED: RankMe, or an equivalent proposal of effective rank as a label-free representation
quality score. The claim "effective rank is used as a quality proxy" needs one canonical citation and
we do not currently hold a verified one.]`

These metrics describe geometry. Whether geometry tracks the molecular channel is a separate
empirical question, and §4.11 reports four independent instances in which it does not. For this paper
the relevance is narrower: a geometric quality metric cannot substitute for a sensitivity statement,
because it is computed on the representation rather than through the analysis pipeline whose null is
in question.

### 2.7 Reference verification status — a standing pre-submission requirement

Three fabricated citations have previously contaminated this project. Every reference in this section
must be verified against a live bibliographic API before the manuscript is submitted.

Current status:

* **Verified with verbatim quotes and DOIs** during the adversarial prior-art sweep: all references in
  §2.3 and §2.4 (`NOVELTY_SEARCH.md`, `TRACK2_INDUCED_CORRELATION.md` §0). Winkler et al. 2020 was
  verified at *full-text* level.
* **Spot-check verified** during literature quality audits (`v2/research/rebase/quality/*.md`,
  2026-07-29): HEST-1k, THREADS, TITAN, de Jong Robustness Index, Mishra & Lotter, Wong et al.,
  Tizhoosh, plus a further set judged real-but-paywalled (Howard, Dawood, CHIEF, Prov-GigaPath, UNI,
  Virchow, Phikon-v2, PLIP, CONCH, PathChat). Those audits record **0 fabrications detected** across
  the lanes checked, but they are spot-checks, not exhaustive, and venue/award labels were explicitly
  flagged as lower-confidence.
* **Marked `[UNVERIFIED]` above** and not yet checked in any pass: Kather et al. 2020 (exact
  journal/year), Kömen et al. 2024, Carloni et al. 2025, Schmitt et al. 2021, Venet et al. 2011, Leek
  & Storey 2007, Johnson et al. 2007, Muirhead 1982 / Anderson 2003 chapter and page.
* **Marked `[CITATION NEEDED]`**: the effective-rank-as-quality-proxy citation (§2.6); the four papers
  underlying the pre-registered ER-status band (§3.9) are recorded in
  `p1_evidence/inputs/PREREG_known_covariate.json` but have not been re-verified in this pass.

---

## 3. Methods

### 3.1 Cohort, representation and the analysed field of view

All measurements use TCGA. The maximal paired split contains **6,427 patients** (3,118 train / 543
validation / 2,766 test), holding out **whole cancers**: 11 development cancers, 21 held out. An
earlier configuration used 2,530 held-out patients over 21 test cancers; both cohort sizes appear
below and are labelled per table, because they are not interchangeable.

Whole-slide images are represented by frozen **H-Optimus-0** patch tokens (1,536-d). The dilution
experiment (§4.10) uses `concat(mean, std)` over those tokens with **no fitted parameters**, reduced
to 256 dimensions by PCA refit per level on training rows only. Other experiments use exported
256-dimensional representation states from trained models (`d2_h_seed42`, `d2_i_seed42`), which serve
here as specimens rather than as objects of study.

**A correction to the patch specification.** Our patch spec nominally states a 128 µm field. It is
wrong. H-Optimus-0's `timm` `pretrained_cfg` is 224 px at `crop_pct` 0.875, `crop_mode` center, so
`create_transform` centre-crops our 256 px patch to 224. **The encoder sees 112 µm, not 128 µm** —
the analysed field is overstated by 14% by the nominal figure. The transform is resolved from the
same config for all 271,710 TCGA patches, so internal comparability is untouched, but both numbers
now ride in every artifact manifest and the corrected value is what should be quoted.

*Provenance: `v2/calibra/hest.py:60–106`, `v2/tests/test_hest.py::test_effective_field_accounts_for_the_encoder_centre_crop`,
`NOTEBOOK_ENTRIES/spatial_baselines_20260803T0620Z.md`.*

### 3.2 Confound design and cross-fitted residualisation — `v2/calibra/residualise.py`

Categorical covariates are one-hot encoded and numeric covariates standardised with an explicit
missingness indicator, so that "missing" is itself adjustable rather than silently imputed. TCGA
tissue source site is derived from the patient barcode and rare sites are pooled once, centrally,
across all analyses, to prevent singleton-site dummies acting as per-patient indicators.

Residualisation is **cross-fitted**: the nuisance model is never fit on the rows whose residuals it
produces. This is not cosmetic. In-sample residualisation removes more than the confound, and §4.6.5
shows that cross-fitting suppresses the classical degrees-of-freedom inflation to about a third of
its in-sample scale — which is precisely why the structural effect we report is visible rather than
mixed with it.

Three design sizes appear in this paper and must not be conflated:

| label | columns | pooling | cohort | where used |
|---|---:|---|---:|---|
| 99-column | 99 | cancer + TSS, 75 sites kept, `min_site_count = 10` | n = 2,530 | §4.1, §4.3 (Phase 1/1b) |
| 108-column | 108 | cancer + 84 pooled sites, `min_site_count = 10` | n = 2,766 | §4.2, §4.4, §4.5, §4.8–4.10 (Track 1, dilution) |
| 109-column ("anchor") | 109 | cancer + `tss_pool10` | n = 500–6,427 | §4.6, §4.7 (Track 2 sweep) |

*Provenance: `PHASE1_RESULT.md` (99); `TRACK1_NEGATIVE_CONTROLS.md` header and `DILUTION_LOWER_BOUND.md`
§1 (108); `TRACK2_INDUCED_CORRELATION.md` §4 (109). The three designs differ by cohort and by the
number of sites surviving pooling at the respective n and partition; the pooling rule
(`min_site_count = 10`) is identical throughout. Two further bookkeeping differences are recorded
rather than hidden: the earlier n = 2,530 configuration holds out 21 test cancers against **14**
development cancers, while the maximal split holds out 21 against **11**, so the two cohorts are not
nested; and the site-dummy counts (75 / 84 / 76) differ because they are computed on different
partitions. No result in this paper compares a number from one design directly against a number from
another.*

### 3.3 Channel statistic and its null — `v2/calibra/run_calibra.py`

The channel is the **top canonical correlation** at 16 whitened components per side, with the
canonical directions fit on one half and scored on the held-out half, which removes the in-sample
maximisation bias. Both blocks are residualised against the confound design first.

The null destroys **only** patient pairing: `permutation_null` permutes the molecular rows *within
cancer strata*, so cohort-level structure is preserved. Permutation *p* is reported with its
resolution — `p = 1/(n_perm + 1)` — and is never written as "p < 0.05". Headline runs use 2,000
permutations (`p` floor 1/2001 = 0.0005); the dilution sweep uses 300 (`p` floor 1/301 = 0.0033).

**A 16-component top canonical correlation has a capacity floor.** At n = 2,766 the within-cancer
shuffled-pairing null median is **0.1465–0.1483**, not zero (§4.4). Every channel number in this paper
is quoted against that median, never against zero.

### 3.4 Spike construction — `v2/calibra/calibration.py::spike_targets`

For a random unit direction `u` in image space and a direction `v` in molecular space, the
`v`-component of the molecular matrix is **replaced** by a signal correlating with the standardised
image score at exactly `r_true`, before any adjustment. The replacement is rescaled by the raw
component's own standard deviation; without that rescale a residual term carries the ambient
correlation into the readout (see §4.1, defect 2). The spiked targets then pass through the
*identical* residualisation and the recovery is scored **on the planted axis**, `corr(X_res·u,
Y_spiked_res·v)`, not on a maximum over the recovered subspace. Within a draw the same `(u, v)` pair
is reused at every spike level, which is what makes the paired transmission floor computable. The
statistic is **signed**; magnitudes are taken only at reporting time (see §4.1, defect 3).

### 3.5 The two floors

| | pairing | question answered | how computed | reportable as a detection limit? |
|---|---|---|---|---|
| **transmission floor** | paired within draw | *Does the pipeline transmit a signal of this size, or destroy it?* | smallest level whose increment over the same draw's level-zero value is seen in at least the required fraction of draws | **No** |
| **detection floor** | unpaired across draws | *Smallest `r_true` reliably distinguishable given the draw-to-draw variability a single-shot analysis faces* | smallest level clearing the level-zero upper tail across draws | **Yes — this is the quotable one** |

Pairing cancels the draw-specific induced baseline, so the transmission floor is close to noiseless
by construction and typically resolves at the finest level on the grid. A real single-shot analysis
has no paired baseline to subtract. Because both quantities are naturally called "the floor", and
because the paired number is the flattering one, we treat the distinction as a reporting requirement
rather than an implementation detail.

Both floors are **censored** at the grid edges and are reported as inequalities where that applies.

### 3.6 The confound certificate — `v2/calibra/confound_certificate.py`

A representation state is certified against a confound (here pooled tissue source site) by two
tests, both required:

* **per-axis**: out-of-fold balanced accuracy for each axis must not exceed the 95th percentile of a
  ≥ 1,000-draw **within-cancer** label-permutation null, and the axis bootstrap CI must include the
  chance rate;
* **joint**: a linear discriminant over *all* axes must satisfy the same condition.

Section 4.2 shows why the joint test is not optional. The certificate records `certified_on =
{raw | adjusted}`, because certifying on the adjusted state and then exposing the raw axis is exactly
the laundering the certification rule forbids.

### 3.7 Ledger and claim admissibility — `v2/calibra/gates.py`, `v2/calibra/claim_guards.py`

`GateLedger.add` records a graded pass/fail on a **run-validity** condition; `GateLedger.observe`
records a **scientific outcome** with an expectation but no verdict. The separation exists because a
finding wired into a pass/fail gate makes a true negative indistinguishable from a broken pipeline —
the one answer we most need to be able to read. Only FAIL rows are disqualifying. A missing control
is written with `value = NaN`, `note = "inadmissible_<code>"` and status FAIL rather than omitted;
silence is never a pass.

The separation is enforced by a test, not by convention:
`v2/tests/test_track1_battery_ledger.py::test_a_losing_baseline_is_an_observation_and_cannot_move_the_verdict`
plants a baseline that beats us 0.95 to 0.40 and asserts the verdict is identical to the same fixture
without it.

`validate_claim` refuses to mark a claim publishable while any blocker required by its claim kind is
undischarged. Six blockers are encoded — `composition_attribution`, `purity_confound`, `sign_blind`,
`proliferation_deflation`, `single_platform`, `no_external_cohort` — each recording the mechanism by
which the claim goes wrong *while the numbers still look fine*, plus the evidence that discharges it.
An unknown claim kind is inadmissible by default.

### 3.8 Uncertainty on differences — `v2/paired_bootstrap.py`

Between-arm and between-configuration differences are reported with a paired bootstrap resampling the
same rows for both arms, in two modes: patients resampled independently, and cancer clusters
resampled before patients within cluster. **The marginal per-arm interval is not used for decisions;
the paired difference is.** Where a marginal interval is known to be biased (E0's `bootstrap_ci95` is
biased low by 0.04–0.09) it is flagged in the emitted record rather than quoted.

### 3.9 Predeclaration

Every Track 2 and dilution prediction in this paper was written into
`v2/research/rebase/nature/P1_PREDECLARATION.md` and committed (`1c4b4b5`) **before** the
corresponding sweep ran; grading was performed by `p1_evidence/grade_t2.py` against that file. The
predeclaration also records, in advance, the one provenance caveat it does not cover: the P3
predictor form was already on disk and its own docstring states it was written after P1 and P2 failed
on the rank ladder at n = 2,530 / seed 42. **P3 is therefore post hoc on the rank axis** and is
reported as such throughout; it is out of sample on every other n and every other seed, and that is
the only sense in which it is tested here.

The ER-status positive control was pre-registered at 01:45 UTC, before the 01:47 run, in
`p1_evidence/inputs/PREREG_known_covariate.json`, with the acceptance band **[0.78, 0.92]** and point
estimate 0.86 derived from four published papers. The declared weaknesses were also written into the
pre-registration before the run: BRCA is a *development* cancer in the maximal split, so this ran on
`--partition all` and is in-distribution; and only one cancer carries the label, so within-cancer and
pooled coincide.

### 3.10 Compute and reproducibility

All analyses in §4.2–§4.10 are CPU linear algebra on a Lambda A100 box (30 cores; the core count is
what mattered, not the GPU). Draw plans and permutation orders are fixed before dispatch so that
parallelism changes wall-clock only, never the numbers.

**Threading.** Unbounded BLAS threading is ~23× *slower* than single-threaded on this stack (SVD of
2,766 × 256: 0.205 s at 1 thread, 4.80 s at 4). All runs set `OMP_NUM_THREADS=1` and parallelise
across processes.

**Training is not seed-reproducible on this stack, and the analysis path is.** Re-exporting a
surviving checkpoint reproduces its recorded value to five significant figures (0.5861 recorded,
0.58612 re-exported), so the export and readout path is deterministic. Retraining the *same seed* with
the *same configuration* does not reproduce the same model: 0.6214 versus 0.5861, effective rank 23.39
versus 8.68. **No individual point estimate from a trained artifact is quoted as reproducible from
the seed alone; only paired differences within a run are quoted.**

*Provenance: `v2/research/rebase/nature/D2_RESULT.md` §4, §5.*

Test suite: **275 passed** in 44 s with thread caps
(`pytest v2/tests tests -q`, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1`), verified at the time of writing.

---

## 4. Results

### 4.1 The instrument's own failure, and why it is the motivating example

The first version of the readout produced `detection_floor = NaN` and `observed_above_floor = false`
for **every** state on real data — i.e. the instrument's headline capability produced nothing — while
passing **all eleven** of its synthetic self-tests.

Three nested defects, each found by fixing the one above it:

1. **The readout was a maximum, not a measurement.** Recovery was scored with
   `top_canonical_correlation`, a maximum over 16 components per side, while the spike lives on
   exactly one known direction pair. Ambient top-CCA in this data is ~0.97, so a spike of
   `r_true ≤ 0.2` was invisible: level-0 recovery read **0.97 instead of ~0**, and at `r_true = 0.2`
   the measured value *fell*, because replacing the `v`-component destroys existing structure faster
   than a weak spike restores it. *Fix:* score the planted axis.
2. **The spike only partially replaced its target.** `spike_targets` computed `y + outer(a_new − a, v)`
   with `a` standardised but `y·v` not, leaving a residual `(σ − 1)·a` carrying the ambient
   correlation. This put level 0 at 0.099 instead of ~0 and attenuated an `r_true = 0.6` spike to
   0.27. *Fix:* rescale the update by the raw component's own standard deviation.
3. **The readout took an absolute value before pairing.** Confound-induced correlation has a random
   sign, so under `|r|` a draw with a negative baseline gets *smaller* when a positive spike is added,
   destroying the paired comparison and reintroducing undefined floors. *Fix:* signed statistic,
   magnitudes only at reporting time.

*Provenance: `v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §0; `PHASE1_RESULT.md` header;
`HANDOFF_PHASE_D.md` §0; code in `v2/calibra/calibration.py`.*

The methodological point is not that we had a bug. It is that **a test suite a broken instrument
passes is not a control.** Eleven synthetic self-tests all passed because they were written against
synthetic data with no ambient cross-modal structure — exactly the condition under which a
maximum-over-components readout and a direction-matched readout agree. A spike pushed through the
real pipeline, against real ambient structure, is the control that discriminates. Every self-test in
this class should be understood as calibrating the *test*, not the instrument.

### 4.2 The adjustment is verified, not assumed — twice, on two different confounders

**Cancer type.** From the residualised representation, out-of-fold balanced accuracy for cancer type
falls from **0.463 raw to 0.035 adjusted**, against a chance rate of **0.048** (1/21 test cancers) at
n = 2,530. The adjusted value is *below* chance, which is the expected behaviour of cross-fitted
residualisation against the variable being removed.

*Provenance: `v2/research/rebase/nature/PHASE1_RESULT.md`, "Validity checks passed".*

**Tissue source site.** The identical certificate, run on raw and on cancer+TSS cross-fitted
residuals of the same states, n = 2,766, 85 pooled TSS classes (chance = 1/85 = 0.0118), 108-column
design, ≥ 1,000 within-cancer label permutations:

| artifact | state | joint LDA **raw** | joint LDA **adjusted** | joint null p95 (adj) | per-axis max raw → adj | breaching axes raw → adj |
|---|---|---:|---:|---:|---|---|
| d2_h | **wsi_biology** | **0.3633** | 0.0118 | 0.0528 | 0.0532 → 0.0123 | 17 → 0 |
| d2_h | full_biology | 0.2630 | 0.0101 | 0.0668 | 0.0548 → 0.0107 | 60 → 0 |
| d2_h | rna_biology | 0.2563 | 0.0074 | 0.0654 | 0.0506 → 0.0139 | 58 → 0 |
| d2_i | **wsi_biology** | **0.2348** | 0.0052 | 0.0418 | 0.0511 → 0.0104 | 43 → 0 |
| d2_i | full_biology | 0.2689 | 0.0085 | 0.0758 | 0.0551 → 0.0102 | 61 → 0 |
| d2_i | rna_biology | 0.2744 | 0.0079 | 0.0732 | 0.0495 → 0.0106 | 48 → 0 |

*Provenance: `v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md` §T1.3;
`NOTEBOOK_ENTRIES/t13_confound_certificate_20260803T0152Z.md`;
`v2/research/rebase/nature/GATE_LOG.md` rows `T1.3_site_certificate_{raw,adjusted}::*`; run outputs
`p1_evidence/track1/certificate_{raw,adjusted}/` on persistent NFS. Command:
`confound_certificate --partition test --n-permutations 1000 --n-boot 200 --n-boot-axes 8`, with and
without `--residualise`.*

Two things must be said and neither is optional.

**(1) The must-fail control did not fail on the raw representation. This is reported as a defect.**
Every raw joint permutation *p* is at the 1/1001 resolution floor: not one of 1,000 within-cancer
label permutations reached the observed joint accuracy in any state. The raw representations are
partly a site code. Consequently no raw axis of these artifacts may be exposed to a user, and
`e0_basis_transfer.py:923`'s `G3.5 = unavailable_no_site_labels` is closed with two rows — raw FAIL,
adjusted PASS.

**(2) A per-axis certificate is not sufficient, and this is a finding about certification rules.**
The leak is **smeared, not concentrated**. The best single axis anywhere reaches 0.055 — 4.6× chance —
and the *median* axis (0.030–0.033) sits *below* its own permutation-null p95. The joint discriminant
over all 256 axes reaches 0.235–0.363, i.e. **20–31× chance** and roughly 2× the joint null p95. Had
the T1.3 criterion been implemented literally as our own specification wrote it — per-axis only —
`wsi_biology` on `d2_h` would have shown 17 of 256 axes breaching and a per-axis maximum of 0.053,
arguably dismissible, while the real leak is 31× chance in the joint direction. **The joint test must
be a required field of the certificate schema, not an optional extra.**

The adjustment discharges the leak completely: joint accuracy falls **21–45×** to at or below the
chance rate, with **zero** breaching axes in all six state × artifact combinations. The defect is
therefore a property of the raw representation, and no adjusted number in this paper is reading site.

### 4.3 …and the adjustment does not destroy signal

Injected-signal attenuation — the slope of recovered against injected correlation on the planted axis
— is ≈ 1 in every measurement we have:

| run | design | cohort | draws | attenuation |
|---|---|---:|---:|---|
| Track 1 negative-control battery, six states | 108-column cancer + pooled TSS | n = 2,766 | 40 | **0.974 – 1.039** |
| Track 2 floor sweep, twelve real design × n cells | cancer, cancer+TSS at four pooling thresholds | n = 1,000 / 2,530 / 6,427 | 40 | **1.07 – 1.12** |
| Track 2, Gaussian structureless control | `gaussian_k99` | n = 1,000 / 2,530 / 6,427 | 40 | **1.000 – 1.001** |
| Phase 1b targeted readout, seven states | 99-column cancer + TSS | n = 2,530 | 40 | 0.944 – 1.228 |
| Dilution sweep, seven levels | 108-column cancer + pooled TSS | n = 2,766 | 20 | 0.855 – 1.130 |

*Provenance: `TRACK1_NEGATIVE_CONTROLS.md` §T1.7(c); `TRACK2_INDUCED_CORRELATION.md` §8;
`PHASE1B_TARGETED_READOUT.md` §3; `DILUTION_LOWER_BOUND.md` §2. The five rows are different runs on
different cohorts and design sizes and are not pooled.*

The Phase 1b range is the widest because it predates the Track 1 protocol; the Track 1 range
(0.974–1.039) is the tightest and is the one to quote. Slopes slightly above 1 mean the adjustment
removes variance that was diluting the planted axis. The dilution row falls below 1 only at the two
heaviest contamination levels (0.855 at d = 0.60, 0.863 at d = 0.80), which is expected.

**This single measurement retires the objection that killed three earlier theses on this project.**
The predeclared falsifier for it — "attenuation far from 1 under a differently constructed confound
design at comparable rank" — was tested by eleven designs spanning 15× in rank, and all eleven sit at
the same induced baseline with attenuation ≈ 1.

### 4.4 Chance is 0.147, not 0

Destroying patient pairing does not take a 16-component top canonical correlation to zero. It takes it
to the capacity floor that 16 components fitted on 2,766 patients produce by construction:

| artifact | state | adjusted top-CCA | null median | null p95 | excess over null median | permutation *p* |
|---|---|---:|---:|---:|---:|---:|
| d2_h | full_biology | 0.8890 | 0.1465 | 0.1645 | 0.7425 | 0.0005 |
| d2_h | rna_biology | 0.8874 | 0.1463 | 0.1642 | 0.7411 | 0.0005 |
| d2_h | wsi_biology | 0.6052 | 0.1483 | 0.1685 | 0.4569 | 0.0005 |
| d2_i | full_biology | 0.8479 | 0.1466 | 0.1638 | 0.7013 | 0.0005 |
| d2_i | rna_biology | 0.8533 | 0.1468 | 0.1643 | 0.7065 | 0.0005 |
| d2_i | wsi_biology | 0.4703 | 0.1472 | 0.1659 | 0.3231 | 0.0005 |

*Provenance: `TRACK1_NEGATIVE_CONTROLS.md` §T1.6; `GATE_LOG.md` rows
`T1.6_modality_shuffled_pairing::*`. 2,000 within-cancer permutations, `p` floor 1/2001 = 0.0005.*

`permutation_p = 0.0005 = 1/2001` throughout: **no permutation of two thousand reached the observed
value in any state.** The resolution is quoted with the *p*; it is never written "p < 0.05".

**Quoting `wsi_biology = 0.4703` without saying that chance is 0.147 overstates the effect by a factor
a reviewer will find immediately.** Every channel number in this paper is quoted against the
permutation null median. The median is cohort- and capacity-specific: 0.1465–0.1483 at n = 2,766 with
16 components, 0.145–0.147 in the dilution sweep at the same setting, 0.151–0.158 at n = 2,530 with 60
permutations, and 0.140 for the D2 readout, which shuffles patient rows of the residualised target
matrix rather than permuting within cancer strata. These are different procedures at different n; do
not carry one across.

### 4.5 Injection-certified floors — and the scale on which they are *not* comparable to the channel

| quantity | value | censoring |
|---|---|---|
| transmission floor (paired), all six Track 1 states | **0.01** — the finest level on the grid | censored from below: report as **≤ 0.01** |
| detection floor (unpaired), Track 1 main run | **0.30** (`d2_h`) / **0.40** (`d2_i`) against `frozen_rna_targets` | — |
| detection floor, across five target blocks on the same artifacts | **0.20 – 0.40** | — |
| detection floor, twelve real design × n cells | **0.25 – 0.30** | — |
| detection floor, structureless `gaussian_k99` at n = 1,000 / 2,530 / 6,427 | **0.050 / 0.015 / 0.010** | — |
| permutation *p*, all states | **0.0005 = 1/2001** | resolution floor |

*Provenance: `TRACK1_NEGATIVE_CONTROLS.md` §T1.7(c); `GATE_LOG.md` rows `T1.7c_spike_recovery::*` and
`T1.2_baseline_block::*` (the per-block `detection_floor=` field, which ranges 0.2–0.4 on the same
artifacts); `TRACK2_INDUCED_CORRELATION.md` §8.*

The detection floor is a property of the **(representation × target block × design × n)** cell, not of
the representation alone. On `d2_h::wsi_biology` alone it reads 0.4 for the PBS block, 0.2 for the
random-dictionary block, 0.4 for the PCA block and 0.3 for the curated-pathway block. Reporting a
single "the detection floor of this representation" would be wrong.

**The one thing the floors do not license, stated plainly.** `observed_above_floor = 0` for every
state, and that is the *correct* answer rather than a failure. The floor is in **single random
direction** correlation units — `floor_scale = targeted_single_direction` is emitted in every summary
for exactly this reason — while the headline channel is a **16-component multivariate maximum**.
Measured through a *random* direction pair, the real channel reads `observed_matched_direction` of
**−0.028 to +0.036** (Track 1; −0.09 to +0.02 in Phase 1b), i.e. a random pair sees nothing, because
the channel is concentrated in particular directions. **This repository therefore contains no
measurement of the real channel, in the floor's own units, that exceeds the floor.** Comparing the
0.47–0.61 channel against the 0.30–0.40 floor directly was the original defect of §4.1 and must not be
done. What the floors certify is the *sensitivity of the pipeline*, not the significance of the
observed channel; the latter is certified by the permutation null of §4.4.

### 4.6 The induced-correlation floor

#### 4.6.1 The mechanism is an identity, verified draw for draw

For `corr(u, v) = 0` the residualised correlation is the multivariate partial-correlation identity

```
r_induced = − R_s · R_a · cos θ / sqrt((1 − R_s²)(1 − R_a²))
```

with `R_s², R_a²` the design's cross-fitted R² for the two scores and `cos θ = corr(ŝ, â)`. Across
**all 270 sweep cells**, the maximum absolute disagreement between the pipeline's measured value and
this closed form evaluated on the *same* planted `(u, v)` pairs is **8.6 × 10⁻¹⁶** (median
3.1 × 10⁻¹⁶), with a per-draw Pearson correlation of **1.000000**. An independent verification against
a standalone simulation (300 draws, n = 2,530, 99-column dummy design) gives a maximum absolute error
of **7.4 × 10⁻¹⁶**.

There is no empirical question about the mechanism. **We concede this immediately and in the text:
the result is one line of algebra, derivable on a whiteboard by anyone with a first-year regression
course.** The open questions are how its factors scale, and whether the magnitude exceeds the
classical degrees-of-freedom term.

*Provenance: `TRACK2_INDUCED_CORRELATION.md` §1; `NOVELTY_SEARCH.md` §1.0. The standalone
verification scripts (`sim.py`, `sim3.py`) were scratchpad artifacts and **should be reproduced into
the repository before submission** — see §5.*

#### 4.6.2 It is structural, not degrees-of-freedom bookkeeping — the decisive predeclared falsifier passes

If a structureless design of matched rank at matched n reproduced the effect, the phenomenon would be
fully covered by the classical N − R result and we would have nothing. This was predeclared as P5,
the falsifier that decides whether the track has content.

| n | k_eff | real cancer+TSS | **same design, rows permuted** | **Gaussian, k = 99** | real / permuted |
|---:|---:|---:|---:|---:|---:|
| 500 | 32.7 | 0.0866 | 0.0087 | 0.0164 | 9.9× |
| 1,000 | 45.7 | 0.0804 | 0.0065 | 0.0101 | 12.4× |
| 2,000 | 87.4 | 0.0809 | 0.0040 | 0.0057 | 20.2× |
| **2,530** | **104.3** | **0.0748** | **0.0037** | **0.0035** | **20.4×** |
| 4,000 | 144.4 | 0.0789 | 0.0031 | 0.0026 | 25.7× |
| 6,427 | 215.1 | 0.0718 | 0.0020 | 0.0016 | 35.4× |

*Provenance: `TRACK2_INDUCED_CORRELATION.md` §2; predeclared in `P1_PREDECLARATION.md` §B (commit
`1c4b4b5`); graded by `p1_evidence/grade_t2.py`. Outputs under `p1_evidence/track2/`.*

The predeclared bar was: structureless arms ≤ 0.025 at n = 2,530 **and** ≥ 3× smaller than real.
Measured 0.0037 and 0.0035, ratio **20.4×**. The row-permuted design has **identical rank and
identical column marginals** and differs only in having no relationship to any patient. The
separation *widens* with n because the structureless arms decay like a sampling term while the real
design does not.

#### 4.6.3 It is a bias, not a sampling fluctuation — and our own scaling prediction was wrong

| n | 500 | 1,000 | 2,000 | 2,530 | 4,000 | 6,427 |
|---|---:|---:|---:|---:|---:|---:|
| induced \|r\| (median of 3 seeds) | 0.0866 | 0.0804 | 0.0809 | 0.0748 | 0.0789 | 0.0718 |

`|r|(6,427) / |r|(2,530) = 0.960`, against a predeclared bar of ≥ 0.60 and a pure-sampling expectation
of √(2530/6427) = 0.627. A 17% decline over 13× more patients, not the 72% decline a sampling term
would give. **More data does not remove it** (predeclared P4: PASS).

**P0 — the plan's own guess, `|r| ~ k/n` — is falsified.** Fitted log–log exponents on real designs
are `b_k_eff = +0.288` (predicted +1.0) and `b_n = −0.180` (predicted −1.0), R² = 0.398. The module's
own derived alternative (`b_k = −0.5`, `b_n = 0`) is also wrong on the rank axis. **We were wrong
twice and both are reported.**

#### 4.6.4 Design rank is the wrong axis; jointly-explained variance is the right one

Rank ladder at n = 2,530, cohort held fixed, rank moved over 15×:

| design | k | k_eff | k_eff_shared | R²_x | R²_y | **measured \|r\|** |
|---|---:|---:|---:|---:|---:|---:|
| none | 0 | 0 | 0 | 0 | 0 | 0.0000 |
| tss_pool50 | 2 | 0.0 | 0.0 | 0.00 | 0.00 | 0.0003 |
| **cancer** | 33 | 31.0 | 9.04 | 0.297 | 0.473 | **0.0844** |
| cancer + tss_pool50 | 35 | 31.0 | 9.05 | 0.297 | 0.474 | 0.0844 |
| **tss_pool10 (no cancer)** | 76 | 74.0 | 9.73 | 0.166 | 0.245 | **0.0379** |
| **cancer + tss_pool10 (anchor)** | 109 | 104.3 | 10.26 | 0.321 | 0.470 | **0.0748** |
| cancer + tss + dx_year | 111 | 106.3 | 10.27 | 0.321 | 0.470 | 0.0744 |
| cancer + tss + dx_year + age | 113 | 108.3 | 10.29 | 0.320 | 0.469 | 0.0732 |
| cancer + tss + dx_year_cat | 141 | 132.3 | 10.45 | 0.313 | 0.464 | 0.0743 |
| cancer + tss_pool3 | 283 | 269.5 | 11.28 | 0.309 | 0.450 | 0.0729 |
| cancer + tss_pool1 | 511 | 450.9 | 11.85 | 0.305 | 0.445 | 0.0727 |

*Provenance: `TRACK2_INDUCED_CORRELATION.md` §4. Single artifact `d2_h_seed42`, single state
`wsi_biology`, three seeds per cell.*

1. **Rank is nearly irrelevant.** `k_eff` 31 → 451 (14.5×, 478 extra columns) moves the effect −14%.
2. **Cancer type alone produces essentially all of it** — 0.0844 against an anchor of 0.0748; the
   anchor is in fact slightly *lower*. The 478 tissue-source-site columns add nothing.
3. **What tracks it is how much the design explains of BOTH sides.** `tss_pool10` has more than twice
   the rank of `cancer` but half the induced correlation, because its R² over the two modalities is
   0.166/0.245 instead of 0.297/0.473. `k_eff_shared` — the effective dimension of the *jointly*
   explained subspace — moves only 9.0 → 11.9 across the whole ladder, which is why the effect is flat.

The practical warning is therefore **not** "keep your nuisance model small". It is: *the induced
correlation is set by how much of both modalities the confound explains, not by how many columns you
spend on it, and it does not shrink with n.*

A predictor of the form `0.6745 · R_x · R_y / √k_eff_shared` reproduces the measured value to a median
ratio of **0.886** (p10–p90 0.76–1.07, log₁₀ RMS error 0.079) over 169 real-design cells, using only
the (design, X, Y) spectra — no draws, no spike, no recovery curve. Competing a-priori forms
under-predict by 2–3× and drift (median ratios 1.98 and 2.98, log₁₀ RMS 0.365 and 0.526). This makes
an induced-correlation floor computable for a new cohort *before* the instrument is run. **Its
provenance is stated and not hidden: it is post hoc on the rank axis** (§3.9); what is genuinely out
of sample is every other n (500–6,427) and every other seed (43, 44), and it holds there.

#### 4.6.5 Cross-fitting is why the structural term is visible

For structureless designs the classical N − R sampling scale is `0.6745/√(n − R)`. Measured over 24
structureless cells, the ratio measured/classical has median **0.379** (p10 0.238, p90 0.699).
Cross-fitted residualisation delivers only about **a third** of the in-sample degrees-of-freedom
inflation the Winkler/Muirhead literature describes. The structural effect we report is ~20× what
remains of the classical mechanism *because* we residualise correctly; an in-sample analysis would
see the two mixed together. The one regime where cross-fitting cannot help is a design that nearly
spans the sample: `gaussian_k600` at n = 500 (`k_eff` = 499) gives 0.055, the largest structureless
value in the grid, though still far below the classical prediction of 0.617 there.

#### 4.6.6 Estimator robustness — the predeclared prediction FAILED

| n | n_splits | α = 0.01 | α = 1.0 | α = 100 |
|---:|---:|---:|---:|---:|
| 2,530 | 2 | 0.0828 | 0.0807 | 0.0380 |
| 2,530 | 5 | 0.0825 | 0.0817 | 0.0531 |
| 2,530 | 10 | 0.0805 | 0.0802 | 0.0573 |
| 2,530 | 20 | 0.0815 | 0.0809 | 0.0592 |
| 6,427 | 2 | 0.0922 | 0.0910 | 0.0632 |
| 6,427 | 5 | 0.0901 | 0.0904 | 0.0725 |
| 6,427 | 10 | 0.0905 | 0.0907 | 0.0760 |
| 6,427 | 20 | 0.0893 | 0.0894 | 0.0772 |

*Provenance: `TRACK2_INDUCED_CORRELATION.md` §7; predeclared as P6 in `P1_PREDECLARATION.md`.*

Predeclared: **< 25% relative movement** across the full grid. **Measured 55.7% at n = 2,530. FAILED.**
The falsifier's own ">2× movement" threshold is also marginally breached (max/min = 2.18 at
n = 2,530), and that is reported rather than rounded down.

Decomposed: fold count is irrelevant (≤ 2.4% spread over `n_splits` ∈ {2, 5, 10, 20} at fixed α);
shrinkage over the usable range is irrelevant (α 0.01 vs 1.0 differ by ≤ 3%); **the entire failure is
α = 100**, a 30–53% reduction. The mechanism is not an artefact and **must not be sold as
robustness**: α = 100 on a one-hot design at n = 2,530 heavily under-fits the nuisance model, the
design then explains less of *both* modalities, and the identity requires the induced correlation to
fall with `R_s R_a`. **Under-adjusting reduces the induced correlation and leaves the confound in.
That is a trade-off, not robustness.**

The quotable statement is narrower than the prediction was: *the induced correlation is invariant to
the cross-fitting scheme and to shrinkage in the range anyone would use (0.07–0.09 for α ∈ [0.01, 1]
at any fold count), and it scales with how much the nuisance model actually removes.*

#### 4.6.7 Replication on a second artifact and across target blocks

The full Track 2 grid is one artifact (`d2_h_seed42`) and one state (`wsi_biology`); the `d2_i` sweep
was queued and did not complete. However, the induced baseline is emitted with every Track 1 block
measurement, giving an independent replication of the magnitude on a second artifact and across five
target blocks:

| block | induced baseline, d2_h `wsi_biology` | induced baseline, d2_i `wsi_biology` |
|---|---:|---:|
| curated pathway (82 cols) | 0.1062 | 0.0984 |
| PBS (128) | 0.0705 | 0.0768 |
| random dictionary (128) | 0.0806 | 0.0937 |
| PCA basis (128) | 0.0833 | 0.0604 |
| gene-label-shuffled (128) | 0.0834 | 0.0720 |

*Provenance: `GATE_LOG.md`, the `induced_baseline=` field of every `T1.2_baseline_block::*` row.
Across all 101 ledger rows the induced baseline spans 0.0557–0.1062.*

The baseline differs by block, and it does so in the direction the mechanism predicts: blocks whose
coordinates the confound design explains less produce less induced correlation. This is a
within-paper confirmation of §4.6.4 obtained from a different experiment.

**Reconciliation of two ranges that appear in our own records.** Phase 1b reported "0.067–0.140" for
the 99-column design at n = 2,530; Track 2 reports 0.0748 for the anchor at the same n. These are
different statistics on different runs: the first is the observed **spread across draws** within one
run, the second is the **median over 40 draws**. The Track 1 block table above shows the same spread
(0.056–0.106) across blocks. Quote the median with its draw spread; do not quote the range as if it
were an interval estimate.

### 4.7 More patients will not buy sensitivity

| design | n | k_eff | induced | **detection floor** | transmission floor | attenuation | ambient top-CCA |
|---|---:|---:|---:|---:|---:|---:|---:|
| cancer | 1,000 | 30.9 | 0.0847 | **0.25** | ≤ 0.01 | 1.124 | 0.684 |
| cancer | 2,530 | 31.0 | 0.0857 | **0.30** | ≤ 0.01 | 1.109 | 0.668 |
| cancer | 6,427 | 31.0 | 0.0844 | **0.30** | ≤ 0.01 | 1.099 | 0.663 |
| cancer + tss_pool50 | 1,000 | 30.9 | 0.0847 | 0.25 | ≤ 0.01 | 1.124 | 0.684 |
| cancer + tss_pool50 | 2,530 | 31.0 | 0.0857 | 0.30 | ≤ 0.01 | 1.109 | 0.668 |
| cancer + tss_pool50 | 6,427 | 50.0 | 0.0842 | 0.30 | ≤ 0.01 | 1.101 | 0.663 |
| cancer + tss_pool10 | 1,000 | 45.3 | 0.0799 | 0.25 | ≤ 0.01 | 1.112 | 0.683 |
| cancer + tss_pool10 | 2,530 | 102.3 | 0.0774 | 0.25 | ≤ 0.01 | 1.091 | 0.664 |
| cancer + tss_pool10 | 6,427 | 215.1 | 0.0811 | 0.25 | ≤ 0.01 | 1.086 | 0.660 |
| cancer + tss_pool1 | 1,000 | 316.6 | 0.0806 | 0.25 | ≤ 0.01 | 1.070 | 0.677 |
| cancer + tss_pool1 | 2,530 | 446.2 | 0.0790 | 0.30 | ≤ 0.01 | 1.086 | 0.661 |
| cancer + tss_pool1 | 6,427 | 579.9 | 0.0831 | 0.30 | ≤ 0.01 | 1.080 | 0.654 |
| **gaussian_k99** | 1,000 | 99.0 | 0.0097 | **0.050** | ≤ 0.01 | 1.000 | 0.877 |
| **gaussian_k99** | 2,530 | 99.0 | 0.0040 | **0.015** | ≤ 0.01 | 1.000 | 0.860 |
| **gaussian_k99** | 6,427 | 99.0 | 0.0014 | **0.010** | ≤ 0.01 | 1.001 | 0.863 |

*Provenance: `TRACK2_INDUCED_CORRELATION.md` §8. Full level grid 0.0…0.50, 40 draws, 2 seeds, same
driver as every other floor in this paper, so the floors cannot come from a second implementation.*

**The detection floor is set by the induced correlation, not by n.** For every real design it sits at
0.25–0.30 and does *not* improve as n goes from 1,000 to 6,427 — a 6.4× increase in sample size buys
nothing. For the structureless Gaussian design of comparable rank it falls 0.050 → 0.015 → 0.010,
roughly like the sampling scale, exactly as it should when there is nothing structural to hit.

This is the most operationally consequential result in the paper. **An external cohort cannot buy a
lower floor by recruiting more patients**, if its confound design explains both modalities the way
cancer type does here. What lowers the floor is reducing `R_x·R_y/√k_eff_shared` — a cohort whose
nuisance structure is less predictive of both modalities — not a bigger n. Cohort sizing for this
class of analysis should be written against the induced-correlation predictor, not against a power
calculation in n.

### 4.8 Negative controls

Ten controls were specified with their pass criteria fixed in advance. Summary:

| # | control | direction required | verdict |
|---|---|---|---|
| T1.3 | site/scanner prediction, **raw** states | must FAIL | ❌ **DID NOT FAIL — reported as a defect** (§4.2) |
| T1.3 | site/scanner prediction, **adjusted** states | must FAIL | ✅ fails as required |
| T1.4 | random gene sets, floor-scale statistic | must FAIL | ✅ fails as required, 0/90 |
| T1.4 | random gene sets, fitted-direction statistic | observation | ⚠️ **reach 76–82% of real gene sets** |
| T1.5(ii) | shuffled gene labels, attribution | must collapse | ✅ collapses (by construction) |
| T1.5(i) | shuffled gene labels, subspace | must persist | ✅ persists — **and the pass is damaging** |
| T1.6 | modality-shuffled pairing | must FAIL | ✅ fails as required, at 1/2001 resolution |
| T1.7(a) | RNA→RNA circular control | must PASS | ✅ passes |
| T1.7(b) | known covariate at published strength | must PASS | ✅ passes (§4.9) |
| T1.7(c) | synthetic spike above the floor | must PASS | ✅ passes |

*Provenance: `TRACK1_NEGATIVE_CONTROLS.md`, headline table; all rows in `GATE_LOG.md`.*

#### 4.8.1 Modality-shuffled pairing

Reported in §4.4. `p = 1/2001` in every state, against a null median of 0.1465–0.1483.

#### 4.8.2 Random gene sets — the control that answers differently on two statistics

Ninety `RANDOM_CONTROL__` columns matched on **training-only** per-gene mean, variance and PC1
loading — stricter than same-size random sets.

On the **primary, spec-literal statistic** (random image direction, i.e. the units the floor is
actually measured in), graded against this run's own unpaired detection floor:

| state | control median | control p95 | control max | exceedances / 90 |
|---|---:|---:|---:|---:|
| d2_h full_biology | −0.0695 | 0.0389 | 0.0719 | **0** |
| d2_h rna_biology | −0.0728 | — | — | **0** |
| d2_h wsi_biology | −0.0341 | — | — | **0** |
| d2_i full_biology | −0.0238 | — | — | **0** |
| d2_i rna_biology | −0.0233 | — | — | **0** |
| d2_i wsi_biology | −0.0132 | — | — | **0** |

0/90 in every state against a 5% ceiling: the control fails as required.

On a **second statistic** — a cross-validated *fitted* image direction per target column, recorded as
an observation rather than a gate because grading it against the floor would be a category error —
the answer is very different:

| state | random-control median | real-target median | **control / real** |
|---|---:|---:|---:|
| d2_h full_biology | 0.4821 | 0.6323 | **0.762** |
| d2_h rna_biology | 0.4788 | 0.6288 | **0.761** |
| d2_h **wsi_biology** | 0.2158 | 0.2804 | **0.770** |
| d2_i full_biology | 0.4740 | 0.5790 | **0.819** |
| d2_i rna_biology | 0.4762 | 0.5863 | **0.812** |
| d2_i **wsi_biology** | 0.1642 | 0.2164 | **0.759** |

*Provenance: `TRACK1_NEGATIVE_CONTROLS.md` §T1.4; `GATE_LOG.md` rows `T1.4_random_gene_sets::*` and
`T1.4_random_vs_real_fitted_direction::*`.*

**Covariate-matched random gene sets are read at 76–82% of the level real curated gene sets are read
at, in every state on both artifacts.** The gap is consistent and real, but three quarters of the
apparent per-target molecular legibility is reproduced by a gene set matched only on marginal mean,
variance and PC1 loading. Two consequences, both central and both against us:

* **Any per-gene-set legibility claim must be stated as a difference against the matched random
  control, with a CI — never as an absolute correlation.**
* The ratio is **invariant to contamination**: across dilution levels d = 0 → 0.80 it reads 0.815,
  0.810, 0.799, 0.792, 0.798, 0.777, 0.727 (§4.10). Contamination removes real and random-control
  signal at essentially the same rate, so the non-specificity is a **property of the readout, not of
  patch quality**, and cleaning up the patches would not fix it.

#### 4.8.3 Shuffled gene labels — the pass is the problem

Attribution collapses as required: median |Spearman| between true and shuffled per-axis gene rankings
is 0.0069 / 0.0073 / 0.0077 over three shuffle draws against a bar of ≤ 0.05, and the strictly harder
best-match statistic (max over all 128 shuffled axes) is 0.033. This is true **by construction** — the
shuffle permutes basis rows after the fit — and is a build-integrity check, not a finding.

The subspace persists, and that is damaging:

| shuffle seed | held-out top-CCA, true | shuffled | paired difference | CI95 of difference |
|---|---:|---:|---:|---|
| 1 | 0.5411 | 0.5600 | −0.0189 | [−0.0489, +0.0384] |
| 2 | 0.5411 | 0.5360 | +0.0051 | [−0.0564, +0.0418] |
| 3 | 0.5411 | 0.4771 | +0.0640 | [−0.0260, +0.0988] |

*Provenance: `TRACK1_NEGATIVE_CONTROLS.md` §T1.5; 500-patient bootstrap; the true value's own CI95 is
[0.4874, 0.5962]. Build integrity: the unshuffled reconstruction matches the frozen scores at
r = 0.99999999999999 per column, enforced by `v2/build_shuffled_pbs_targets.py`.*

**All three paired-difference CIs cover zero**, and on one draw the shuffled dictionary scores higher.
After the row permutation the target block **is** a spectrum-matched random projection of the same
expression matrix, so the honest reading is: *any spectrum-matched 128-dimensional projection of this
expression matrix is as legible as the fitted dictionary.* This agrees with §4.8.2's 76–82% by an
entirely independent route.

A method note that we report rather than choose between: the containment test in the original
specification is weak — a wide CI passes it by accident — and by containment the control passes on 2/3
draws while by paired difference all 3 cover zero. The paired difference is what should be quoted.
There is also a **readout dependence**: at `n_components = 32` with a 500-draw bootstrap the shuffled
block was indistinguishable from the fitted one on all three draws, whereas at `n_components = 16` the
fitted block wins on 3 of 8 comparisons with a CI excluding zero (§4.13). The defensible statement is
*the shuffle costs the dictionary little and sometimes nothing, never more than a few hundredths of a
canonical correlation* — enough to void gene-level attribution, but not the flat tie the n = 32
readout alone suggested. Both readouts are reported; neither was selected after the fact.

### 4.9 Positive controls

**(a) RNA→RNA, circular by construction.** `--require-rna-positive-control` passes on both artifacts;
`channel_gate_failures` is empty and `gates_pass` is true. Adjusted top-CCA 0.8874 / 0.8533.

**(b) A known-legible covariate at an independently published strength.** MSI, TP53 and consensus
subtype were all unusable from the data on disk — the TCGA PanCan clinical mirror carries
`microsatellite_instability` only as an assay-performed flag (7 `YES`, 74 `NO`, 6,171 `NONE`) with no
MSI-H/MSI-L/MSS calls, and no mutation or subtype table exists on either machine. The substitution is
recorded, not silent: **TCGA-BRCA ER status by IHC** (690 labelled, 528+/162−), with PR as a second
anchor.

| artifact | state | adjustment | within-cancer AUROC | CI95 | null p95 | verdict |
|---|---|---|---:|---|---:|---|
| d2_h | **wsi_biology** | raw | **0.8781** | [0.8457, 0.9115] | 0.546 | **PASS** |
| d2_h | wsi_biology | cancer+TSS | 0.8714 | [0.8379, 0.9055] | 0.546 | PASS |
| d2_i | **wsi_biology** | raw | **0.8667** | [0.8360, 0.8971] | 0.542 | **PASS** |
| d2_i | wsi_biology | cancer+TSS | 0.8644 | [0.8340, 0.8946] | 0.544 | PASS |
| d2_h | full_biology | raw | 0.9195 | [0.8901, 0.9455] | 0.544 | PASS |
| d2_i | full_biology | raw | 0.9401 | [0.9127, 0.9640] | 0.545 | PASS (marginal) |

*Provenance: `TRACK1_NEGATIVE_CONTROLS.md` §T1.7(b); `GATE_LOG.md` rows `T1.7b_known_covariate::*`;
pre-registration `p1_evidence/inputs/PREREG_known_covariate.json`, written 01:45 UTC before the 01:47
run. Band [0.78, 0.92], point estimate 0.86, from four published papers recorded in that file: Naik et
al. *Nat Commun* 2020;11:5727 (0.92 internal, 0.86 TCGA external); Rawat et al. *Sci Rep*
2020;10:7275; Shamai et al. *JAMA Netw Open* 2019;2:e197700; Couture et al. *npj Breast Cancer*
2018;4:30. `[CITATION NEEDED: these four have not been re-verified against a live bibliographic API in
this drafting pass.]`*

The image-only numbers sit essentially on the literature point estimate. Two observations matter more
than the pass itself:

* **The measured within-cancer chance level is 0.542–0.546, not 0.5.** Grading against an assumed 0.5
  would have been wrong by four points. Chance rates in this setting must be measured, not assumed.
* The declared weaknesses were written into the pre-registration **before** the run: BRCA is a
  *development* cancer in the maximal split, so this ran on `--partition all` and is
  **in-distribution**; and only one cancer carries the label, so within-cancer and pooled coincide and
  the lineage-guessing protection is not exercised. The RNA-containing states reach 0.92–0.94, which
  is expected by construction (ER status is close to a monotone function of *ESR1*) and carries no
  morphological claim.

**(c) Synthetic spike above the floor.** Transmission floor ≤ 0.01 for every state; the pipeline
transmits a paired signal of r = 0.01 without destroying it. Attenuation 0.974–1.039 (§4.3).

**The admissibility rule these controls exist for.** A negative result on this project is reportable
only if the positive control passed **in the same run, on the same data, through the same code path**.
This is stricter than "we also ran a control", because it forbids inheriting a control from a prior
run, a prior configuration or a prior commit, and it is enforced in code rather than in prose.

### 4.10 A dose–response demonstration: contamination with information-free patches

Each patient's patch bag is contaminated with **same-cancer, different-patient** tumour patches from
the same H-Optimus store, drawn donor-slide-first (a contiguous stretch of the wrong tissue, not an
unrealistic mosaic), with levels **nested** so the curve is not confounded by independent draw noise.
Representation: `concat(mean, std)` over 1,536-d tokens, **no fitted parameters**, PCA-reduced to 256
dimensions refit per level on train rows only (retained variance 0.879–0.923 at every level). Cohort:
6,427 patients, 238,610 tumour patches, 7,644 tumour slides; 2,766 evaluated on `test`.

| requested d | achieved d | adjusted top-CCA | held-out top-CCA | ratio to d = 0 | **null-corrected ratio** | detection floor | attenuation | effective rank | perm *p* |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.000 | **0.5573** | 0.4932 | 1.000 | **1.000** | 0.20 | 1.130 | 196.2 | 0.0033 |
| 0.10 | 0.091 | 0.5571 | 0.5017 | 0.9996 | **0.999** | ≥ 0.40 | 0.985 | 194.1 | 0.0033 |
| 0.20 | 0.211 | 0.5447 | 0.5129 | 0.977 | **0.968** | ≥ 0.40 | 1.003 | 190.5 | 0.0033 |
| 0.30 | 0.302 | 0.5190 | 0.4986 | 0.931 | **0.905** | ≥ 0.40 | 1.057 | 187.5 | 0.0033 |
| 0.40 | 0.400 | 0.4774 | 0.4619 | 0.857 | **0.804** | ≥ 0.40 | 1.014 | 184.7 | 0.0033 |
| 0.60 | 0.600 | 0.3971 | 0.3680 | 0.713 | **0.607** | ≥ 0.40 | 0.855 | 176.5 | 0.0033 |
| 0.80 | 0.800 | **0.2844** | 0.1922 | 0.510 | **0.333** | ≥ 0.40 | 0.863 | 161.2 | 0.0033 |

*Provenance: `v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2;
`NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md`. Null-corrected ratio =
`(observed − null median) / (observed₀ − null median₀)`; the permutation null median is 0.145–0.147 at
every level. 300 within-cancer permutations, `p` floor 1/301 = 0.0033. Outputs under
`p1_evidence/dilution/`.*

The **raw ratio flatters the surviving channel**; the null-corrected column is the one to quote. The
half-loss point on the null-corrected curve is **d ≈ 0.68**: roughly two thirds of the bag must be
replaced with same-cancer tumour from other patients before half the channel is lost. Every level
clears its own within-cancer shuffled-pairing null at the run's resolution floor.

**Two of four predeclared predictions were falsified, both in the same direction — the channel is far
more robust than predicted.** D2 (proportional decline, `channel(d)/channel(0) ≈ 1 − d` to ±0.15 for
d ≤ 0.40) is falsified for d ≥ 0.21: measured 1.000 / 0.977 / 0.931 / 0.857 against predicted 0.909 /
0.789 / 0.698 / 0.600. D4 (at d = 0.60, retains > 0 but < 55%) is falsified: 71.3% raw, 60.7%
null-corrected. D1 (monotonicity) passes. D3 (non-decreasing detection floor) passes but the
measurement is **censored** — the level grid tops out at 0.40 and the floor reads 0.40 from d = 0.09
onward, so it is "≥ 0.40", and D3's monotonicity is supported only over the one step it could resolve.

**On the word "lower bound".** The file recording this experiment is named `DILUTION_LOWER_BOUND.md`
and its own §4 withdraws that phrasing. What is **measured** is the cost of contaminating a bag with
same-cancer, different-patient tumour at matched stain, matched preparation (100% FFPE diagnostic on
both sides) and matched pipeline — there is no domain shift. What is **assumed and untested** is that
this is the most benign possible non-informative contaminant. A mechanism argues the other way:
normal tissue is off-manifold similarly across *all* patients, so it adds a near-constant offset to a
mean-pooled bag, damaging between-patient variation less than a patient-specific random tumour shift
does. If that mechanism dominates, this arm is an **upper** bound. **The correct phrasing is "the cost
of preparation-matched, information-free contamination", not "a lower bound",** and the number must
not be quoted with the words "lower bound" attached unless this paragraph travels with it. The
`pooled`, `matched` and `dx_normal` normal-tissue arms that would settle it need GPU re-embedding and
were not run.

### 4.11 Effective rank does not track information content — four independent instances

| # | manipulation | effective rank | information measure | source |
|---|---|---|---|---|
| 1 | covariance-decorrelation term added | 49.9 → 103.3, **+107%** | within-cancer specificity 0.1366 → 0.1367, **flat** | `v2/research/rebase/ENGINE_CLD.md`; `HANDOFF_BUILD_AGENT.md` §"Covariance-decorrelation fix" (3 seeds) |
| 2 | `full` → `programme_only` supervision | 38.48 → 32.06, **−17%** | held-out molecular top-CCA 0.4768 → 0.4748, **−0.002** | `PHASE1B_TARGETED_READOUT.md` §5 |
| 3 | full training schedule vs contrastive-only | `z_biology` matrix rank **pinned at 16/16 in every arm** | representation collapsed to a point (within-modality off-diagonal cosine 0.7089 → 0.9999; retrieval 0.062 → **0.000**, i.e. *below* chance 0.062) | `NOTEBOOK.md` 2026-08-02 01:20 UTC |
| 4 | patch-bag contamination d = 0 → 0.80 | 196.2 → 161.2, **−18%** | null-corrected channel 1.000 → 0.333, **−67%** | `DILUTION_LOWER_BOUND.md` §6 |

*Provenance: as in the final column. Instances 3 and 4 are the two limbs of the dissociation — rank
pinned while information collapses, and rank drifting down gently while information collapses — so
both directions now have evidence.*

Two honesty notes attach to this table and both must travel with it. Instance 2's two arms were
**not verified as matched on epochs, learning rate and step budget** (gate G0.4 in our own protocol),
so it is suggestive rather than causal. Instance 1 comes from an earlier generation of the codebase
and a different benchmark statistic. Instances 3 and 4 are the strongest. Separately, effective rank
is **unstable across seeds** for the same configuration (9.14 to 34.12 in three seeds of one arm),
which is itself an argument against using it as a summary.

**For this paper the point is narrow:** a geometric quality metric is computed on the representation,
not through the analysis pipeline whose null is in question, and therefore cannot substitute for a
sensitivity statement. A fuller treatment belongs to a companion paper and is not claimed here.

### 4.12 An application: a supervision-target ablation read through the instrument

To show the instrument doing work, we report an ablation it was used to decide. Two training arms
differ only in the molecular targets they are supervised on — curated Hallmark pathway scores versus
perturbation-basis coordinates — with three seeds each, 40 epochs, the 6,427-patient split, and the
readout run at 16 components with a 2,000-repeat paired patient and cancer-cluster bootstrap on
n = 2,766 held-out patients.

**On 40 targets neither arm trained on** (`heldout_pathway` + `immune_tme` + `tumour_state`):

| seed | Hallmark | PBS | Δ (PBS − H) | patient CI₉₅ | cancer CI₉₅ |
|---|---:|---:|---:|:---:|:---:|
| 42 | 0.6126 | 0.4800 | **−0.1325** | [−0.1605, −0.0993] | [−0.1792, −0.0632] |
| 43 | 0.5970 | 0.4882 | **−0.1089** | [−0.1460, −0.0749] | [−0.1623, −0.0118] |
| 44 | 0.5983 | 0.4757 | **−0.1226** | [−0.1502, −0.0866] | [−0.1653, −0.0411] |

**Negative control — the 90 `random_control` targets:**

| seed | Hallmark | PBS | Δ (PBS − H) | patient CI₉₅ | cancer CI₉₅ |
|---|---:|---:|---:|:---:|:---:|
| 42 | 0.4671 | 0.4572 | −0.0099 | [−0.0591, +0.0123] | [−0.1055, +0.0099] |
| 43 | 0.4681 | 0.4400 | −0.0280 | [−0.0719, −0.0048] | [−0.0905, +0.0232] |
| 44 | 0.4637 | 0.4369 | −0.0268 | [−0.0697, +0.0003] | [−0.0969, +0.0285] |

*Provenance: `v2/research/rebase/nature/D2_RESULT.md` §2, §3. Run `d2_v3`, outputs under
`~/e0_run/d2_v3/bootstrap/` and `~/e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json` on persistent NFS.
Unrestricted 90-target readout gives −0.1359 / −0.1077 / −0.1192, so stratifying to untrained targets
changes the answer by at most 0.0034. Permutation null for every group is 0.140 at
`permutation_p = 0.005`, the floor for 200 draws.*

Four properties make this a usable demonstration of the instrument rather than a biological claim.
(i) The gap survives at full size on targets neither arm trained on. (ii) It is *smaller* on the
Hallmark arm's own supervision targets (−0.1091 / −0.0787 / −0.1149) than on the untrained ones, so
scoring on one arm's own homework *understates* the gap rather than manufacturing it. (iii) The
negative control is 4–13× smaller, with the cancer CI covering zero 3/3 and the patient CI 2/3.
(iv) A capacity explanation is **contradicted**: effective rank is H 23.39 / 28.77 / 9.14 against PBS
14.87 / 34.12 / 9.11 — in seed 43 the losing arm has *higher* rank, and in seed 44 the two are equal
to two decimals and the losing arm still loses.

Two caveats travel with the number. The negative control is small but not exactly zero and points the
same way, so of order 10–20% of the gap may be generic representation quality rather than supervision
content. And per §3.10 no individual point estimate here is reproducible from the seed alone; only the
paired within-run difference is quoted.

**This is presented as an ablation supporting the instrument's utility. It is not a biological claim
and no biological conclusion is drawn from it in this paper.**

### 4.13 Target blocks scored through one instrument, with the losses reported

Four target blocks were built against the same 6,427-patient cohort, gene order, log transform and
development-only fit discipline (`v2/baseline_target_common.py`) and scored through the identical
instrument (`--partition test --levels 0.0,0.05,0.10,0.20,0.40 --n-draws 16 --n-components 16
--n-permutations 500 --seed 42`). Every difference carries a 400-draw paired bootstrap in which both
blocks are scored on the *same* resample, so the interval is on the difference. These are
`GateLedger.observe` rows — scientific outcomes, never validity gates.

Held-out top-CCA, PBS block minus baseline block:

| artifact | state | baseline block | PBS | baseline | difference | CI95 | verdict |
|---|---|---|---:|---:|---:|---|---|
| d2_h | wsi_biology | random dictionary | 0.5032 | 0.4551 | +0.0481 | [−0.0177, +0.0693] | TIE |
| d2_h | wsi_biology | **PCA basis** | 0.5032 | **0.5520** | **−0.0488** | **[−0.0734, −0.0183]** | **BASELINE WINS** |
| d2_h | wsi_biology | shuffled labels (s1) | 0.5032 | 0.5146 | −0.0114 | [−0.0601, +0.0445] | TIE |
| d2_h | wsi_biology | shuffled labels (s2) | 0.5032 | 0.5187 | −0.0155 | [−0.0561, +0.0318] | TIE |
| d2_h | full_biology | random dictionary | 0.8417 | 0.8102 | +0.0315 | [+0.0241, +0.0653] | PBS wins |
| d2_h | full_biology | **PCA basis** | 0.8417 | **0.8776** | **−0.0359** | **[−0.0483, −0.0236]** | **BASELINE WINS** |
| d2_h | full_biology | shuffled labels (s1) | 0.8417 | 0.8140 | +0.0277 | [+0.0197, +0.0660] | PBS wins |
| d2_h | full_biology | shuffled labels (s2) | 0.8417 | 0.8085 | +0.0332 | [+0.0116, +0.0632] | PBS wins |
| d2_i | wsi_biology | random dictionary | 0.4605 | 0.4108 | +0.0497 | [+0.0251, +0.1372] | PBS wins |
| d2_i | wsi_biology | **PCA basis** | 0.4605 | 0.4905 | −0.0300 | [−0.0429, +0.0053] | TIE |
| d2_i | wsi_biology | shuffled labels (s1) | 0.4605 | 0.4245 | +0.0360 | [+0.0111, +0.0926] | PBS wins |
| d2_i | wsi_biology | shuffled labels (s2) | 0.4605 | 0.4317 | +0.0288 | [−0.0030, +0.0887] | TIE |
| d2_i | full_biology | random dictionary | 0.8634 | 0.8487 | +0.0147 | [−0.0004, +0.0283] | TIE |
| d2_i | full_biology | **PCA basis** | 0.8634 | **0.8714** | **−0.0080** | **[−0.0233, −0.0001]** | **BASELINE WINS (marginal)** |
| d2_i | full_biology | shuffled labels (s1) | 0.8634 | 0.8408 | +0.0227 | [−0.0029, +0.0288] | TIE |
| d2_i | full_biology | shuffled labels (s2) | 0.8634 | 0.8374 | +0.0260 | [+0.0057, +0.0383] | PBS wins |

*Provenance: `TRACK1_NEGATIVE_CONTROLS.md` §T1.1/T1.2; `GATE_LOG.md` rows `T1.2_baseline_block::*`.
Note that the source table's header row is mis-labelled (it prints "baseline" twice); the column order
is (baseline name, PBS value, baseline value, difference), verified against the `heldout=` field of the
corresponding ledger rows. **This draft prints the corrected header.***

Reading down the three opponents: against a size- and spectrum-matched random dictionary, PBS wins
2/4 and ties 2/4, never losing; against a gene-label-shuffled version of itself, PBS wins 3/8 and ties
5/8, never losing; **against ordinary PCA of the same expression matrix, PBS loses 3/4 with a CI
excluding zero and ties the fourth — it never wins.** The PCA basis is fit on development rows only,
through the same transform, and is capacity-matched at 128 columns, so neither leak nor capacity
explains it.

**Not built, and this is a gap:** the **text-prior (GenePT-style)** and **capacity-matched
cell-composition** blocks, both of which need an external resource (a gene text-embedding table; a
deconvolution signature matrix) that is on neither machine. The plan flagged these two as the most
likely to slip and they slipped. `claim_guards.composition_attribution` therefore remains
undischarged. The zero-parameter cancer-type-mean baseline is **degenerate by construction** on this
split, because whole cancers are held out and every test patient receives the global training mean;
that is a property of the split, recorded as such, and the comparison is meaningful only on
`--partition all`.

### 4.14 The ledger

`v2/calibra/track1_battery_ledger.py` assembles every row above into
`v2/research/rebase/nature/GATE_LOG.md`. Current contents: **101 rows — 62 gates, 39 observations, 7
failed gates** (the six raw site certificates, and gene-label-shuffle seed 3's containment test). The
39 observations include the entire baseline table of §4.13, so the PCA loss is recorded in the same
append-only file as the gates **without ever having been able to quarantine the run**. That is the
whole point of the separation: registering a scientific outcome as a pass/fail gate would mark the run
FAILED exactly when the science came back "no", making a true negative indistinguishable from a broken
pipeline.

---

## 5. Limitations

Stated in descending order of how much they constrain the paper.

1. **No external cohort.** Every measurement here is TCGA, which carries documented site and scanner
   effects. This is a deliberate scope decision, not an oversight, but it means the instrument is
   **demonstrated, not validated**. `claim_guards.no_external_cohort` is undischarged for every
   morphology result on this project; `legible_axis` and `gene_attribution` claims are inadmissible
   and none is made. A related hard constraint we measured: a logistic classifier on frozen
   H-Optimus-0 embeddings separates TCGA from an external spatial cohort (HEST) at **AUC 0.99999**,
   while a within-TCGA split of the same size and classifier gives **0.5012**. Any cross-cohort
   transfer must declare that number.
   *Provenance: `NOTEBOOK_ENTRIES/spatial_baselines_20260803T0620Z.md`.*

2. **The induced-correlation phenomenon is not ours.** It is the Yule/FWL partial-correlation
   identity, already published as an applied warning in fMRI (Murphy 2009), GWAS (Aschard 2015; Dahl
   2019), genomics (Nygaard 2016; Li 2023), and folklore for `removeBatchEffect` (Smyth 2020). Winkler
   et al. 2020 warns about it in exactly the cross-modal setting and supplies a fix. **Only the
   magnitude characterisation is ours**, under correct cross-fitted residualisation of exactly
   orthogonal signals, and it must be claimed that narrowly. The identity is one line of algebra and
   we concede that in the text rather than letting a reviewer point it out.

3. **~76–82% of the per-target channel is reproduced by covariate-matched random gene sets**, and the
   ratio is invariant to contamination, so it is a property of the readout rather than of patch
   quality. This guts per-target pathway specificity. Any per-gene-set legibility claim must be a
   difference against the matched control with a CI. An independent route — gene-label shuffling —
   agrees: all three paired-difference CIs cover zero.

4. **The predeclared estimator-robustness prediction (P6) failed**, at 55.7% against a 25% bar, and
   the falsifier's own 2× threshold was marginally breached (2.18). The failure is entirely on the
   Ridge shrinkage axis (α = 100) and has a mechanistic explanation, but under-adjusting reducing the
   induced correlation is a **trade-off, not robustness**, and is reported as such.

5. **Two of our own a-priori scaling predictions were falsified.** P0 (`|r| ~ k/n`) is wrong on both
   axes: exponents +0.288 and −0.180 against +1 and −1. The module's own second derivation is also
   wrong on the rank axis. Design rank is the wrong axis. The predictor that does work is **post hoc
   on the rank axis** and is labelled as such wherever it appears.

6. **Training is not seed-reproducible on this stack.** Retraining seed 42 with an identical
   configuration gave 0.6214 against a recorded 0.5861, and effective rank 23.39 against 8.68, while
   re-exporting the surviving checkpoint reproduced the recorded value to five significant figures.
   The analysis path is deterministic; the training path is not. **Only paired within-run differences
   are quoted anywhere in this paper.**

7. **The floors are not comparable to the channel.** `observed_above_floor = 0` for every state, which
   is correct: the floors are in single-random-direction correlation units and the channel is a
   16-component multivariate maximum. There is no measurement in this repository of the real channel,
   in the floor's own units, that exceeds the floor. The paper therefore certifies pipeline
   sensitivity, not the significance of the channel; the latter rests on the permutation null.

8. **Censoring.** The transmission floor is censored from below at the finest grid level (≤ 0.01
   generally, ≤ 0.05 in the dilution grid) and the dilution detection floor is censored from above at
   0.40. Neither can be published as a function of its independent variable without a finer grid at
   both ends.

9. **Single-artifact, single-state coverage for the main sweep.** The 270-cell Track 2 grid is one
   artifact (`d2_h_seed42`) and one state (`wsi_biology`), with three seeds per cell; the second-artifact
   sweep was queued and did not complete. The magnitude replicates on the second artifact through the
   Track 1 block measurements (§4.6.7), but the rank ladder, the n ladder and the estimator sweep do
   not.

10. **Recorded substitutions.** No TCGA purity table exists on either machine, so the
    `["cancer","tss","purity"]` rank point could not be built; an expression-derived surrogate was
    **rejected** because it is computed from the very RNA targets that form Y and would inflate `R_a`
    by construction, manufacturing the effect under study. Replaced by `dx_year` and `age`, plus the
    TSS pooling threshold. MSI, TP53 and consensus subtype were unavailable for the known-covariate
    control and ER/PR status were substituted. The text-prior and cell-composition baseline blocks
    were not built.

11. **The dilution arm is one representation and one seed.** `raw_hoptimus_meanstd` has no fitted
    parameters, which is right for isolating the effect of the patches, but a trained attention
    aggregator could plausibly down-weight foreign patches and be more robust. The number is a
    property of unweighted mean pooling, not of the modality. Single seed (42), single draw of donor
    assignments; nesting makes the curve internally consistent but gives no error bar on
    level-to-level differences.

12. **The strongest positive control is in-distribution.** ER status is only labelled in BRCA, which
    is a *development* cancer in the maximal split, so the control ran on `--partition all`, and
    within-cancer and pooled evaluation coincide. Declared in the pre-registration before the run.

13. **The analysed field of view is 112 µm, not 128 µm** (§3.1). Any statement about the
    window-to-assay area ratio must use 5.28×, not 6.90×.

14. **The closed-form verification scripts are not in the repository.** The 7.4 × 10⁻¹⁶ standalone
    verification was run from scratchpad files (`sim.py`, `sim3.py`). The 8.6 × 10⁻¹⁶ in-pipeline
    verification over 270 cells *is* reproducible from repository code. The scratchpad scripts must be
    reproduced into the repo before the identity is quoted in a submitted manuscript.

15. **Reference verification is incomplete** (§2.7). Three fabricated citations have previously
    contaminated this project; every reference must be verified against a live bibliographic API
    before submission, and the `[UNVERIFIED]` / `[CITATION NEEDED]` markers resolved or the
    corresponding sentences removed.

---

## 6. Conclusion

A confound-adjusted cross-modal analysis has three numbers it usually does not report and cannot be
interpreted without: what its adjustment actually removed, what its adjustment cost the signal, and
what effect size it could not have seen. We built an instrument that measures all three by injecting a
correlation of known strength onto a named direction pair and pushing it through the identical
pipeline.

Applied to TCGA morphology and expression, the answers are: the adjustment removes what it claims to
(cancer-type balanced accuracy 0.463 → 0.035 against chance 0.048; joint site accuracy 0.3633 → 0.0118
against chance 0.0118, with zero breaching axes in six states); it costs the signal essentially
nothing (attenuation 0.974–1.039); and it could not have seen a single-direction correlation below
0.25–0.40. The last number is not improved by recruiting patients — it is pinned across a 6.4×
increase in n — because it is set by a residualisation-induced correlation floor that is structural
rather than a sampling term, 20–35× larger than what a matched-rank structureless design produces.

The negative-control battery caught a real confound leak, located it in the raw representation,
showed the adjustment discharging it, and — most usefully — showed that the per-axis certificate our
own specification called for would have passed that leak. It also showed that three quarters of the
apparent per-target molecular legibility is reproduced by covariate-matched random gene sets, and that
this is invariant to patch contamination.

What we do **not** claim: that residualisation-induced correlation is a new phenomenon; that any
representation state here is biologically interpretable; that any per-gene or per-pathway attribution
survives; or that the instrument has been validated on a cohort other than TCGA. The value on offer
is a protocol: report a paired transmission floor and an unpaired detection floor, quote every channel
number against its measured permutation null rather than against zero, require a joint confound test
alongside the per-axis one, and compute the induced-correlation baseline your own design will produce
before interpreting the correlation you measure.

---

## Appendix A — predeclared predictions and their grades

All written into `v2/research/rebase/nature/P1_PREDECLARATION.md`, committed `1c4b4b5`, before any of
the corresponding sweeps ran; graded by `p1_evidence/grade_t2.py`.

| # | prediction | bar | measured | verdict |
|---|---|---|---|---|
| P0 | `\|r_induced\| ~ k/n`; log-log exponents (+1, −1) | — | +0.288, −0.180 (R² 0.398) | **FALSIFIED** |
| P1 | `0.6745·κ/√k_eff` | ratio ≈ 1 | median 1.98 (p10–p90 1.15–3.98) | fails |
| P2 | `0.6745·R_s R_a/√k_eff` | ratio ≈ 1 | median 2.98 (1.65–6.41) | fails |
| P3 | `0.6745·R_x R_y/√k_eff_shared` **(post hoc on rank)** | ratio ≈ 1 | median 0.886 (0.76–1.07) | holds, with provenance caveat |
| P4 | effect is a bias: `\|r\|(6427) ≥ 0.6·\|r\|(2530)` | ≥ 0.60 | 0.960 (sampling would give 0.627) | **PASS** |
| P5 | structureless matched-rank design induces ≤ 0.025 and ≥ 3× less than real | ≤ 0.025, ≥ 3× | 0.0037 / 0.0035, 20.4× | **PASS (decisive)** |
| P6 | estimator knobs move `\|r\|` by < 25% | < 25% | 55.7% at n = 2,530; max/min 2.18 | **FAILED** |
| D1 | channel declines monotonically in d | no non-monotone step beyond CI | 0.5573 → 0.2844, flat then strictly monotone | **PASS** |
| D2 | `channel(d)/channel(0) ≈ (1−d)` ± 0.15 for d ≤ 0.40 | ±0.15 | 1.000 / 0.977 / 0.931 / 0.857 vs 0.909 / 0.789 / 0.698 / 0.600 | **FALSIFIED** (d ≥ 0.21) |
| D3 | detection floor non-decreasing in d | monotone | 0.20 → 0.40 then flat | PASS, **censored** |
| D4 | at d = 0.60 channel retains > 0 but < 55% | < 55% | 71.3% raw, 60.7% null-corrected | **FALSIFIED** |

Three of eleven predeclared predictions failed and one holds only with a stated provenance caveat.

## Appendix B — provenance index

| section | primary evidence file | run outputs |
|---|---|---|
| §4.1 instrument failure | `v2/research/rebase/nature/PHASE1B_TARGETED_READOUT.md` §0; `PHASE1_RESULT.md` | `runs/calibra_v3_targeted` |
| §4.2 confound certificate | `TRACK1_NEGATIVE_CONTROLS.md` §T1.3; `NOTEBOOK_ENTRIES/t13_confound_certificate_20260803T0152Z.md`; `PHASE1_RESULT.md` | `p1_evidence/track1/certificate_{raw,adjusted}/` |
| §4.3 attenuation | `TRACK1_NEGATIVE_CONTROLS.md` §T1.7(c); `TRACK2_INDUCED_CORRELATION.md` §8 | `p1_evidence/track1/`, `p1_evidence/track2/` |
| §4.4 permutation null | `TRACK1_NEGATIVE_CONTROLS.md` §T1.6; `GATE_LOG.md` | `p1_evidence/track1/` |
| §4.5 floors | `TRACK1_NEGATIVE_CONTROLS.md` §T1.7(c); `GATE_LOG.md` | `p1_evidence/track1/` |
| §4.6 induced correlation | `TRACK2_INDUCED_CORRELATION.md`; `NOVELTY_SEARCH.md`; `P1_PREDECLARATION.md` | `p1_evidence/track2/` |
| §4.7 floors vs n | `TRACK2_INDUCED_CORRELATION.md` §8 | `p1_evidence/track2/` |
| §4.8 negative controls | `TRACK1_NEGATIVE_CONTROLS.md` §T1.4, §T1.5, §T1.6 | `p1_evidence/track1/` |
| §4.9 positive controls | `TRACK1_NEGATIVE_CONTROLS.md` §T1.7 | `p1_evidence/track1/`, `p1_evidence/inputs/PREREG_known_covariate.json` |
| §4.10 dilution | `DILUTION_LOWER_BOUND.md`; `NOTEBOOK_ENTRIES/dilution_foreign_tumour_20260803T0355Z.md` | `p1_evidence/dilution/` |
| §4.11 effective rank | `ENGINE_CLD.md`; `PHASE1B_TARGETED_READOUT.md` §5; `NOTEBOOK.md` 2026-08-02 01:20; `DILUTION_LOWER_BOUND.md` §6 | — |
| §4.12 supervision ablation | `D2_RESULT.md` | `~/e0_run/d2_v3/` |
| §4.13 baseline blocks | `TRACK1_NEGATIVE_CONTROLS.md` §T1.1/T1.2; `GATE_LOG.md` | `p1_evidence/track1/` |
| §4.14 ledger | `v2/calibra/track1_battery_ledger.py`; `GATE_LOG.md` | — |

All `p1_evidence/` paths are relative to
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/` on the Lambda A100 box
(`ubuntu@150.136.45.194`), persistent NFS.

## Appendix C — code index

| component | module |
|---|---|
| confound design, cross-fitted residualisation | `v2/calibra/residualise.py` |
| spike construction, recovery curve, both floors, permutation null | `v2/calibra/calibration.py` |
| orchestration and manifest | `v2/calibra/run_calibra.py` |
| confound certificate (per-axis + joint) | `v2/calibra/confound_certificate.py` |
| induced-correlation sweep and predictors | `v2/calibra/induced_correlation_sweep.py`, `v2/calibra/analyse_induced_correlation.py` |
| health-gate / observation ledger | `v2/calibra/gates.py`, `v2/calibra/track1_battery_ledger.py` |
| claim admissibility | `v2/calibra/claim_guards.py` |
| paired bootstrap | `v2/paired_bootstrap.py` |
| gene-label shuffle control | `v2/calibra/gene_label_shuffle_control.py`, `v2/build_shuffled_pbs_targets.py` |
| known-covariate control | `v2/calibra/known_covariate_control.py` |
| baseline target blocks | `v2/baseline_target_common.py`, `v2/build_pca_basis_targets.py`, `v2/build_random_dictionary_targets.py` |
| dilution construction and reduction | `v2/research/dilution/build_dilution_artifact.py`, `v2/research/dilution/reduce_dilution.py` |
| field-of-view correction | `v2/calibra/hest.py` |
