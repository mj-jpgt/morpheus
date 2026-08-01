# P1 — CALIBRA: instrument / methods paper

*Front-matter draft. Every number below traces to `HANDOFF_PHASE_D.md` §0 or to the P1 evidence ledger
in `NOTEBOOK.md`. Nothing here is a submission; sections marked PENDING have no data behind them yet.*

---

## 1. Working title

**Calibrated auditing of morphology–molecular claims: transmission and detection floors for
confound-adjusted cross-modal analyses.**

Alternatives:

1. *What effect size would this analysis have missed? A spike-recovery instrument for confound-adjusted
   multimodal representation claims.*
2. *Shared residualisation induces correlation between orthogonal signals: a calibration instrument for
   morphology–molecular studies.*

---

## 2. Abstract (~200 words)

Computational pathology increasingly reports that tissue morphology predicts molecular state. Because
tissue source site, scanner and batch confound these analyses, the standard remedy is to residualise
both modalities against a covariate design. Neither the positive nor the negative results of such
analyses are usually accompanied by a sensitivity statement, so a small or absent effect cannot be
distinguished from an analysis that was not sensitive enough to see one. We describe CALIBRA, a
calibration instrument that injects a synthetic signal of known strength into the molecular side, pushes
it through the identical pipeline including residualisation, and reports the effect size the analysis
would have missed. The instrument separates two quantities that are routinely conflated: a paired
*transmission floor*, which asks whether the pipeline destroys signal, and an unpaired *detection floor*,
which is the quotable detection limit. Applying it produced a methodological observation we have not
seen reported: residualising two orthogonal signals through a shared confound design induces correlation
between them, measured at 0.067–0.140 for a 99-column cancer-plus-tissue-source-site design at n = 2,530.
Any study that residualises two modalities against shared covariates and then correlates the residuals
has a non-zero null it is not computing. We also report the instrument's own failure: an earlier readout
passed all eleven synthetic self-tests while returning undefined floors on real data.

---

## 3. Introduction (~1200 words)

### 3.1 A field of claims without sensitivity statements

The claim that tumour morphology carries molecular information is now routine. Deep models predict point
mutations, microsatellite instability and molecular subtype directly from haematoxylin-and-eosin slides
(Kather et al., *Nature Cancer* 2020); slide-level foundation models report molecular-profile and
outcome prediction across cohorts (CHIEF, Wang et al., *Nature* 2024; Prov-GigaPath, Xu et al., *Nature*
2024); and encoders are now pretrained against paired genomic and transcriptomic profiles specifically so
that molecular meaning is internalised into the image representation (THREADS, Vaidya et al.,
arXiv:2501.16652, 2025). Benchmarks have followed, most prominently HEST-1k (Jaume et al., NeurIPS 2024,
arXiv:2406.16192), which standardises "predict expression from morphology" into a comparable task.

The reporting convention across this literature is a point estimate of association, sometimes with a
bootstrap interval, occasionally with a permutation null. What is almost never reported is the answer to
the complementary question: *what effect size would this analysis have failed to detect?* Without that
number, a null result is uninterpretable. A study that residualises confounds, measures a small
cross-modal correlation, and concludes that morphology carries little information about a given programme
cannot distinguish "there is no signal here" from "this pipeline could not have recovered a signal of
this size even if it were present". The two conclusions have opposite scientific consequences and the
same published appearance.

This is not a hypothetical failure mode in our own work. Three separate theses on this project were
retired by adversarial review on exactly this objection — that a confound adjustment might have destroyed
the effect it was meant to isolate — and none of them could be defended, because no sensitivity number
existed to defend them with.

### 3.2 Confounding makes the problem worse, not better

The confounding in this setting is documented and substantial. Howard et al. (*Nature Communications*
2021) showed that the TCGA tissue-submitting site is detectable from the slide image, survives colour
normalisation, and correlates with clinical and genomic labels, so that reported accuracy for survival,
driver mutations and stage is in part a site detector. Subsequent work found that the problem is not
solved by scale: hospital signatures persist in foundation-model embeddings and dominate feature-space
distances (Kömen et al., arXiv:2411.05489, 2024); across ten public pathology foundation models, medical
centre identity is encoded more strongly than tissue or cancer type (de Jong et al., arXiv:2501.18055,
2025); representational-similarity analysis finds pronounced slide-dependence and weak disease-dependence
in the same models (Mishra & Lotter, arXiv:2509.15482, 2025); scanner hardware is a separate confounder
again (Carloni et al., arXiv:2507.22092, 2025); and models readily learn scanner type, institution and
slide-preparation date as latent variables (Schmitt et al., *JMIR* 2021). Dawood et al. ("Buyer Beware",
bioRxiv 2024.06.23.600257; *Nature Biomedical Engineering* 2026) add that per-biomarker models predict a
correlated bundle rather than an isolated biomarker.

The standard remedy is to remove the covariate: correct features post hoc (ComBat on deep features,
Murchan et al., *Journal of Pathology Informatics* 2024), or residualise both modalities against a design
matrix of cancer type, site and related covariates before measuring association. This is the right
instinct, and it has a side effect that we have not found quantified anywhere. Residualisation is a
projection. Projecting two signals onto the orthogonal complement of a shared design does not leave their
relationship untouched: it can create one. The larger the design, the more of each signal is removed, and
the more of what remains is determined by the same shared geometry.

The practical consequence is that the null hypothesis of "no cross-modal association after adjustment" is
not "correlation zero". It is some non-zero, design-dependent, draw-dependent quantity that nobody in
this literature computes. A study that residualises morphology and expression against a shared covariate
design and then correlates the residuals is testing against the wrong reference.

### 3.3 A calibration instrument

We take the approach used for spike-in controls in genomics and for simulation-based calibration in
Bayesian methodology — inject something of known magnitude and check what comes back — and apply it to
certify a cross-modal morphology–molecular analysis, which to our knowledge has not been done.
[CITATION NEEDED: the ERCC External RNA Controls Consortium spike-in reference.]
[CITATION NEEDED: the canonical simulation-based-calibration reference.]

The construction is deliberately unclever. For a random unit direction `u` in image space and a direction
`v` in molecular space, we replace the `v`-component of the molecular matrix with a signal correlating
with the image score at exactly the requested strength `r_true`, before any adjustment. The spiked
targets then pass through the *identical* pipeline — the same cross-fitted residualisation against the
same design — and we score the correlation on the planted axis. Sweeping `r_true` gives a recovery curve,
from which we read an attenuation slope (how much the pipeline costs a signal it does receive) and a
floor (the smallest signal the pipeline reliably surfaces). The observed real-data effect can then be read
against that floor rather than against zero.

Applied to a 99-column cancer-plus-tissue-source-site design at n = 2,530, the attenuation slope is
0.94–1.23, i.e. approximately one. The confound adjustment does not destroy the signal. That single
measurement retires the objection that had killed three earlier theses on this project — and it is the
kind of statement no study in this area currently makes about its own adjustment.

### 3.4 Two floors, routinely conflated

The instrument reports two distinct quantities and they are not interchangeable.

The **transmission floor** is paired: within a draw, the same `(u, v)` is used at every spike level, and
the test asks whether a draw's recovered value exceeds *its own* level-zero value. Pairing cancels the
draw-specific induced baseline, so this floor is close to noiseless and typically resolves at the finest
level on the grid. It answers the over-residualisation question — does the pipeline transmit a signal of
this size, or destroy it? It is *not* a statement about statistical power, because a real single-shot
analysis has no paired baseline to subtract.

The **detection floor** is unpaired: it asks whether the recovered value clears the level-zero upper tail
across draws. It carries the full draw-to-draw variability an actual analysis faces, and it is therefore
the conservative, quotable detection limit — approximately 0.2 in single-direction correlation units for
the whole-slide-image channel measured here.

Quoting the paired floor as a detection limit would understate the detection threshold by an order of
magnitude. Because both are naturally called "the floor", and because the paired number is the flattering
one, we treat the distinction as a reporting requirement rather than an implementation detail.

### 3.5 The central methodological finding

The spike is constructed orthogonal to the image score. Level zero should therefore read approximately
zero. On real data it does not. Residualising two orthogonal signals through a shared confound design
induces correlation between them: 0.067–0.140 for the 99-column cancer-plus-tissue-source-site design at
n = 2,530. The magnitude is draw-specific, because it depends on how much of the drawn `(u, v)` pair lies
in the span of the design, and its sign is random.

This is one measurement, on one design, at one n, on one cohort. We state it that way throughout. Its
generality is the paper's main open question, not one of its results. But the implication for practice
does not depend on the exact magnitude: any analysis that residualises two modalities against shared
covariates and then correlates the residuals has a non-zero null it is not computing, and that null
depends on its own design. Reporting a small post-adjustment correlation as evidence of weak biology, or
a moderate one as evidence of strong biology, requires knowing where the induced baseline sits.

The random sign is what forced the paired test. It also caused a defect: taking an absolute value before
pairing meant that a draw whose induced baseline was negative got *smaller* when a positive spike was
added, which destroyed the paired comparison and reintroduced undefined floors.

### 3.6 Admissibility

The final component is a reporting rule. A negative result is reportable only if the positive control
passed in the same run, on the same data, through the same code path. This is stricter than the usual
"we also ran a control" because it forbids inheriting a control from a prior run, a prior configuration
or a prior commit, and it is enforced in code rather than in prose. We also separate health gates, which
decide whether a run is trustworthy, from observations, which record what the run found — because
registering a scientific outcome as a pass/fail gate makes a true negative indistinguishable from a
broken pipeline, which is the one answer we most need to be able to read.

We include the instrument's own failure as evidence for the same argument. An earlier version of the
readout scored recovery as a maximum over sixteen canonical components while the spike lived on one known
direction. Ambient correlation in this data sits near 0.97, so every detection floor came back undefined
on real data — while all eleven synthetic self-tests passed. Three nested defects were involved. A test
suite that a broken instrument passes is not a control; a spike pushed through the real pipeline is.

---

## 4. Related work (~800 words)

### 4.1 Pathology foundation models and morphology-to-molecular prediction

Tile-level self-supervised encoders (UNI, Chen et al., *Nature Medicine* 2024; Virchow, Vorontsov et al.,
arXiv:2309.07778, 2023; Phikon-v2, Filiot et al., arXiv:2409.09173, 2024; H-Optimus, Bioptimus open
release, 2024–2025) and slide-level encoders (Prov-GigaPath, Xu et al., *Nature* 2024; TITAN, Ding et al.,
arXiv:2411.19666, 2024) have made frozen pathology features a commodity input. Vision–language models
extend this to text-conditioned inference (PLIP, Huang et al., *Nature Medicine* 2023; CONCH, Lu et al.,
*Nature Medicine* 2024; MUSK, Xiang et al., *Nature* 2025; PathChat, Lu et al., *Nature* 2024).

The morphology-to-molecular claim proper begins with Kather et al. (*Nature Cancer* 2020) and continues
through CHIEF (Wang et al., *Nature* 2024) and THREADS (Vaidya et al., arXiv:2501.16652, 2025), which
aligns slide embeddings to paired genomics and transcriptomics during pretraining. HEST-1k (Jaume et al.,
NeurIPS 2024) is the reference benchmark for expression prediction from morphology. Phikon-v2 reports
that margins between leading encoders on slide-level biomarker tasks are statistically non-significant,
and Tizhoosh (arXiv:2510.23807, 2025) argues the failures are structural rather than a matter of scale.
None of this work reports what effect size its analysis would have missed; the reporting convention is a
point estimate against a downstream task, not a sensitivity curve. Our contribution is orthogonal to all
of these: we do not propose an encoder or a benchmark, but an audit that any of them could be run through.

### 4.2 Confounding and site effects in TCGA-based studies

Howard et al. (*Nature Communications* 2021) established that submitting-site signatures are learnable
from H&E, survive stain normalisation, and inflate apparent accuracy on survival, mutation and stage.
Schmitt et al. (*JMIR* 2021) broadened the confounder list to scanner type, institution and preparation
date. More recent work shows the problem is not dissolved by foundation-model scale: batch signatures
persist in embeddings (Kömen et al., arXiv:2411.05489, 2024), centre identity is encoded more strongly
than biology (de Jong et al., arXiv:2501.18055, 2025), and representational geometry is organised by
slide rather than disease (Mishra & Lotter, arXiv:2509.15482, 2025). Scanner is a distinct axis (Carloni
et al., arXiv:2507.22092, 2025). Dawood et al. (bioRxiv 2024.06.23.600257; *Nature Biomedical
Engineering* 2026) show that per-biomarker predictions capture co-dependent bundles.

Mitigations divide into image-space normalisation, feature-space correction (ComBat on deep features,
Murchan et al., *Journal of Pathology Informatics* 2024), adversarial removal, and covariate
residualisation. This literature quantifies how much confound signal a remedy removes. It does not
quantify what the remedy does to the relationship between two adjusted modalities — which is the gap this
paper addresses. Wang et al. (*Nature Communications* 2025, PMID 39934114) is the closest in spirit,
judging expression-from-histology methods on cross-study generalisation and downstream survival rather
than within-slide accuracy, and noting explicitly that mean-expression, random-gene and permutation nulls
were not run.

### 4.3 Calibration, sensitivity analysis and negative controls in genomics

Genomics has a mature negative-control culture that computational pathology has partly not imported.
Venet, Dhanasekaran & Sotiriou (*PLoS Computational Biology* 2011) showed that most random gene-expression
signatures are significantly associated with breast-cancer outcome, establishing the random-signature null
as a required control; the Domany group (2018, PMC5839591) extended the effect across TCGA cancers and
attributed it to cohort sub-structure. Random and size-matched gene-set nulls are codified in standard
tooling and protocols (STAR Protocols 2022; NARGAB 2024; Zhang et al., *eLife* 2022, which quantifies
per-method false-positive rates on null signatures). On the perturbation side, simple baselines beat
deep models under held-out-perturbation splits (Wong, Hill & Moccia, *Bioinformatics* 41(6), 2025;
Ahlmann-Eltze, Huber & Anders, 2025, PMC12202205), with the field criticised for lacking appropriate
controls.

Covariate adjustment itself has a long methodological literature — surrogate variable analysis (Leek &
Storey, 2007), ComBat (Johnson, Li & Rabinovitch, 2007), and removal of unwanted variation
(Gagnon-Bartsch & Speed, 2012). These methods are evaluated on whether they remove the unwanted factor
and preserve the factor of interest. We are aware of no work in this line that measures the correlation a
*shared* adjustment induces between two independently adjusted modalities.
[CITATION NEEDED: any prior report that residualising two orthogonal signals against a shared design
induces correlation between the residuals — a null-result search that must be run before submission.]
[CITATION NEEDED: prior use of injected known-strength signal to certify the sensitivity of a
confound-adjusted cross-modal analysis, if any exists.]

### 4.4 Representation-quality metrics

A parallel line proposes geometric proxies for representation quality: alignment and uniformity on the
hypersphere (Wang & Isola, ICML 2020), dimensional collapse in contrastive learning (Jing et al., ICLR
2022), the variance and covariance terms of VICReg (Bardes, Ponce & LeCun, arXiv:2105.04906), neural
collapse in multivariate regression (arXiv:2409.04180, NeurIPS 2024), and effective rank as a scalar
summary of a spectrum (Roy & Vetterli, 2007). In pathology specifically, the Robustness Index (de Jong
et al., 2025) and representational-similarity analysis (Mishra & Lotter, 2025) are proposed as
confound-aware representation diagnostics.

These metrics describe geometry. Whether geometry tracks the molecular channel is a separate empirical
question, and the answer measured on this project is that it does not — a result that belongs to a
companion paper and is not claimed here. For the present paper the relevance is narrower: a geometric
quality metric cannot substitute for a sensitivity statement, because it is computed on the
representation rather than through the analysis pipeline whose null is in question.
[CITATION NEEDED: RankMe or an equivalent proposal of effective rank as a label-free representation
quality score.]

---

## 5. Contributions

1. A spike-recovery calibration instrument that pushes a synthetic signal of known strength through the
   identical analysis pipeline, including residualisation, and reports the effect size that analysis
   would have missed — converting an uninterpretable null into a null read against a measured floor.
2. The measurement that residualising two orthogonal signals through a shared confound design induces
   correlation between them (0.067–0.140; 99-column cancer-plus-tissue-source-site design; n = 2,530),
   which means every study correlating two shared-design residuals is testing against the wrong null.
3. The separation of a paired transmission floor from an unpaired detection floor, with the demonstration
   that the paired quantity is close to noiseless and must never be quoted as a detection limit.
4. A direction-matched readout for spike recovery, replacing a maximum-over-components readout that is
   swamped by ambient cross-modal structure, together with the recorded failure in which the
   maximum-based readout passed all eleven synthetic self-tests while returning undefined floors on real
   data.
5. An executable admissibility layer that encodes claim-level caveats as code rather than prose, refusing
   to mark a claim publishable while a named blocker is undischarged, and a ledger that separates health
   gates (which decide whether a run is trustworthy) from observations (which record what it found).
6. A worked case series of health-gate failures in which the gate was correct and the obvious remedy was
   wrong, supporting the design rule that a gate which cannot fail cleanly certifies nothing.

---

## 6. Methods outline

### 6.1 Confound design and cross-fitted residualisation — `v2/calibra/residualise.py`
Categorical covariates are one-hot encoded and numeric covariates standardised with an explicit
missingness indicator, so that "missing" is itself adjustable rather than silently imputed. TCGA tissue
source site is derived from the barcode and rare sites are pooled once, centrally, across all analyses.
Residualisation is cross-fitted: the nuisance model is never fit on the rows whose residuals it produces,
because in-sample residualisation removes more than the confound.

### 6.2 Spike construction — `v2/calibra/calibration.py` (`spike_targets`)
For a random unit direction `u` in image space and a direction `v` in molecular space, the `v`-component
of the molecular matrix is *replaced* by a signal correlating with the standardised image score at
exactly `r_true`. Replacement is rescaled by the raw component's own standard deviation; without that
rescale a residual term carries the ambient correlation into the readout. Structured spikes use a real
programme loading vector for `v`, since a random direction is the favourable case.

### 6.3 The recovery curve and the two floors — `v2/calibra/calibration.py` (`spike_recovery_curve`)
Spiked targets pass through the identical residualisation and are scored on the planted axis, not on a
maximum over the recovered subspace. Within a draw the same `(u, v)` is reused at every level, which is
what makes the paired transmission floor computable. The transmission floor is the smallest level whose
increment over the same draw's level-zero value is seen in at least the required fraction of draws; the
detection floor is the smallest level clearing the level-zero upper tail across draws.

### 6.4 Permutation null — `v2/calibra/calibration.py` (`permutation_null`)
Patient pairing is destroyed by permuting rows *within strata*, normally cancer type, so that cohort-level
structure is preserved and only the cross-modal pairing is broken. This is the companion to the recovery
curve: a top canonical correlation is a multivariate maximum and is inflated by capacity at finite n, so
an observed value is uninterpretable until chance level under the same design is known.

### 6.5 Health-gate ledger — `v2/calibra/gates.py`
`GateLedger.add` records a graded pass/fail on a run-validity condition; `GateLedger.observe` records a
scientific outcome with an expectation but no verdict. The separation exists because a finding wired into
a pass/fail gate makes a true negative indistinguishable from a broken pipeline. Only FAIL rows are
disqualifying; observations are written to the same append-only ledger for provenance.

### 6.6 Claim admissibility — `v2/calibra/claim_guards.py`
`validate_claim` refuses to mark a claim publishable while any blocker required by its claim kind is
undischarged. Each of the six blockers records the mechanism by which the claim goes wrong *while the
numbers still look fine*, plus the specific evidence that discharges it. An unknown claim kind is
inadmissible by default, and an inadmissible verdict is emitted as a visible status row rather than
silently dropped.

### 6.7 Uncertainty on differences — `v2/paired_bootstrap.py`
Between-arm and between-configuration differences are reported with a paired bootstrap resampling the
same rows for both arms, in two modes: patients resampled independently, and cancer clusters resampled
before patients within cluster. The marginal per-arm interval is not used for decisions; the paired
difference is.

### 6.8 Orchestration and provenance — `v2/calibra/run_calibra.py`
One entry point builds the confound design, runs the channel measurement, the recovery curve and the
permutation null on the same complete-case patient set, and writes a manifest recording the confound
columns, permutation strata, draw and permutation counts, and the input artefact identities.

---

## 7. Results outline

### Main figures and tables

| # | Item | Status | Source |
|---|---|---|---|
| F1 | The instrument's own failure: pre-fix recovery curve returning undefined floors on real data against ambient top-CCA ≈ 0.97, while 11/11 synthetic self-tests pass; the three nested defects and the post-fix curve | **DONE** | `HANDOFF_PHASE_D.md` §0; `v2/calibra/` |
| F2 | Post-fix recovery curve with attenuation slope 0.94–1.23 on the 99-column cancer + tissue-source-site design at n = 2,530 — the confound adjustment does not destroy signal | **DONE** | `HANDOFF_PHASE_D.md` §0 |
| F3 | Induced correlation at spike level zero under shared residualisation: 0.067–0.140, draw-specific in magnitude and random in sign. **One design, one n, one cohort** | **DONE**, single setting | `HANDOFF_PHASE_D.md` §0; P1 ledger |
| T1 | The two floors side by side: paired transmission floor (near-noiseless, resolves at the finest grid level, not a power statement) versus unpaired detection floor (≈ 0.2 in single-direction correlation units for the WSI channel) | **DONE** | `HANDOFF_PHASE_D.md` §0 |
| T2 | Admissibility layer: six blockers with mechanism and discharge condition, 15 tests; gate-versus-observation separation in the ledger API | **DONE** | `v2/calibra/claim_guards.py`, `v2/calibra/gates.py`, `tests/test_claim_guards.py` |
| T3 | Worked health-gate failure case series: a liveness gate straddling the objective warm-up boundary; an overfit gate reading a post-divergence loss; the same gate graded against a stale key queue; the same gate undecidable at 280 patients | **DONE** | `NOTEBOOK.md`, entries 15:31–17:45 |
| F4 | External / second-dataset demonstration of both floors | **PENDING** — not started; no external cohort has been through the instrument | P1 ledger |
| T4 | Negative-control battery: must-fail controls (site/scanner/batch prediction from certified axes; random gene sets clearing the floor; shuffled gene labels leaving attribution intact; modality-shuffled pairing preserving agreement) and must-pass positive controls, reported including losses | **PENDING** — specification only | `MULTIMODAL_EXPANSION.md` §9 |
| F5 | Induced correlation at a second design rank and a second n, to establish whether the effect is general or specific to this design | **PENDING** — not run | P1 ledger, known blockers |

### Supplementary

| # | Item | Status | Source |
|---|---|---|---|
| S1 | Worked legibility ceiling for a curated-pathway readout: at an interior cross-validated ridge penalty of 1e4, 45 of 50 Hallmark axes are predictable from the WSI vector above chance, with a maximum held-out R² of +0.161. **One measurement at one setting**, with a random-CV control ruling out cross-cancer extrapolation as the cause | **DONE** | `NOTEBOOK.md`, 2026-08-01 09:18 |
| S2 | Scale note and worked error: the detection floor is in single-direction correlation units while the headline channel figure is a multivariate maximum over 16 components; comparing them directly was the original defect | **DONE** | `v2/calibra/calibration.py` |
| S3 | Reproducibility: draw plans and permutation orders are fixed before dispatch so that parallelism changes wall-clock only, never the numbers | **DONE** | `v2/calibra/calibration.py` |

**Not in this paper.** No claim is made here about whether representation-geometry metrics track molecular
information, and no claim is made about perturbation-basis supervision. Both belong to companion papers
and neither is required for any statement above.

---

## 8. Open gaps

Everything that must close before submission.

1. **No external cohort has been through the instrument.** Every measurement reported here is TCGA, which
   carries documented site and scanner effects. `claim_guards.no_external_cohort` is undischarged for
   every morphology result on this project. Until a second cohort is run, the instrument is demonstrated,
   not validated, and the draft must say so in the abstract, the discussion and the limitations.
2. **The induced-correlation observation rests on one design at one n.** 0.067–0.140 was measured on a
   single 99-column cancer-plus-tissue-source-site design at n = 2,530. A methods paper asserting a
   general phenomenon needs the same measurement at a second design rank and a second sample size. Until
   then the finding is stated as one measurement at one setting, and the falsifier is explicit: induced
   correlation ≈ 0 at matched design rank and n on a second design would mean the effect is specific to
   this design rather than to shared residualisation.
3. **The negative-control battery is a specification, not run output.** `MULTIMODAL_EXPANSION.md` §9 lists
   the must-fail controls and must-pass positive controls; none has been executed. The paper's central
   admissibility argument — a negative result is reportable only if the positive control passed in the
   same run — is currently made without the paper itself having run that battery.
4. **The detection floor is quoted for one channel.** The ≈ 0.2 figure is the WSI channel under this
   design. Floors for other channels and other designs have not been measured, and the paper must not
   generalise from the one that has.
5. **No analytic account of the induced correlation.** We report it empirically. A derivation predicting
   its magnitude from design rank, n and the geometry of `(u, v)` relative to the design span would
   convert an observation into a result, and would make the second-design experiment a confirmation
   rather than the whole evidential basis.
6. **Two related-work searches have not been run to a null.** Whether anyone has previously reported
   shared-design-induced correlation between residuals, and whether anyone has previously certified a
   confound-adjusted cross-modal analysis with an injected known-strength signal. Both are marked
   `[CITATION NEEDED]` in §4.3 and must be resolved before any novelty language is used.
7. **Citation verification.** Three fabricated citations have previously contaminated this project.
   Every reference in §4 must be verified against a live API before it enters a manuscript, and the
   `[CITATION NEEDED]` markers in §3.3, §4.3 and §4.4 must be either filled from a verified source or the
   corresponding sentence removed.
8. **Venue framing.** The deliverable is an audit procedure, not a finding about tumours, and the paper
   should be written for a methods venue accordingly. No biological claim from this project is admissible
   as of this draft.
