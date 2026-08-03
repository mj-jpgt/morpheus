# Novelty / prior-art search — the two claims in `PHASE1B_TARGETED_READOUT.md`

**Date:** 2026-08-02 · **Scope:** adversarial prior-art sweep, brief was *find the prior art, not reassure*.
**Databases:** Google Scholar / general web (200-call budget exhausted), Europe PMC REST, arXiv, PubMed/PMC,
publisher pages, Bioconductor support archive. Statistics/econometrics, genomics/bioinformatics,
neuroimaging, computational pathology, and (for Q2) astronomy/gravitational-wave.

**Two verdicts up front:**

| | claim under test | verdict | single worst prior work |
|---|---|---|---|
| **Q1** | residualising two orthogonal signals against a shared confound design induces correlation between them; *"Nobody reports this"* | **ALREADY REPORTED** | Winkler et al. 2020, *NeuroImage* 220:117065 |
| **Q2** | injected known-strength spike-in certifying the sensitivity/detection floor of a confound-adjusted cross-modal analysis in biology | **PARTIALLY ANTICIPATED** | Gerard 2020, *BMC Bioinformatics* 21:206 (`seqgendiff`) |

**Required edits to `PHASE1B_TARGETED_READOUT.md` §1:**

- Delete *"Nobody reports this."* It is false in at least five literatures.
- Delete *"A real finding that fell out of the fix."* The phenomenon is a two-line algebraic identity
  (§1.0), verified here to machine precision.
- Keep the measurement. The **number** — 0.067–0.140 for this design at this n, in a cross-modal
  setting, under *correct* residualisation of *exactly orthogonal* signals — is the only part not
  already in print, and it is a calibration, not a discovery.
- Keep the operational conclusion (paired test required). That is a protocol claim and survives intact.

---

# Q1 — confound residualisation inducing cross-signal correlation

## Verdict: **ALREADY REPORTED**

Not "partially anticipated". The mechanism is textbook, it has a closed form, and it has been
independently rediscovered-and-warned-about in fMRI (2009), GWAS (2015), genomics (2016, 2019, 2023),
and named as folklore by the author of `limma` (2020).

**Critically, part (b) of the brief is also gone.** The brief hoped the contribution might be "warning
that it inflates *cross-modal* confound-adjusted measurements". That warning is in print twice, in
exactly the cross-modal setting: **Winkler et al. 2020** (CCA on imaging × non-imaging residuals —
*"residualisation introduces dependencies … leads to inflated error rates"*, plus a fix) and
**Alfaro-Almagro et al. 2021** (UK Biobank imaging × non-imaging — *"spurious associations can be induced
between pairs of otherwise independent variables"*). What remains unpublished is a **number**, not a
finding. See §1.1(2) and §1.1(3).

## 1.0 It is an exact algebraic identity, not an empirical finding

Let `M = I − X(XᵀX)⁻¹Xᵀ` be the residual maker for the design `X`, and `P = I − M`. For any two vectors
with `uᵀv = 0` exactly (the paper's construction),

```
uᵀMv = uᵀv − uᵀPv = −(Pu)ᵀ(Pv)
```

so the residual correlation is

```
corr(Mu, Mv) = − R_u · R_v · ρ / sqrt( (1 − R_u²)(1 − R_v²) )
```

where `R_u`, `R_v` are the multiple correlations of `u`, `v` with the design and `ρ = cos∠(Pu, Pv)`.
This is precisely the multivariate partial-correlation formula
`r_uv·X = (r_uv − R_u R_v ρ)/sqrt((1−R_u²)(1−R_v²))` evaluated at `r_uv = 0`.

**Verified numerically** (scratchpad `sim3.py`, 300 draws, n = 2,530, 99-column dummy design):
closed form vs. observed residual correlation, **max abs error 7.4 × 10⁻¹⁶** — machine precision.
Induced `|r|` median 0.056, range 0.000–0.243. The paper's 0.067–0.140 is a numerical instance of this
identity, and the observed draw-to-draw variation is exactly the variation in `ρ` and `R_u R_v`.

A second simulation confirms the effect is **structural, not finite-sample**: for two *independent random*
vectors residualised on the same 99-column design at n = 2,530, induced `|r|` is mean **0.0031**, max 0.0106
— two orders of magnitude below the reported floor. So the paper's floor is entirely driven by
`(u, v)` loading on the design span, i.e. by `R_u R_v ρ`, which is what §1 already says. This *supports*
the paper's mechanistic account and simultaneously confirms it is the classical partial correlation, not a
degrees-of-freedom artefact that might have been novel.

**Consequence:** any reviewer with a first-year regression course can derive the result on the board.
It cannot be a headline contribution, and the claim of discovery is a liability.

## 1.1 Strongest prior art

### (1) Yule (1907) → the Yule–Frisch–Waugh–Lovell theorem — *the textbook root*

- G. U. Yule, "On the theory of correlation for any number of variables, treated by a new system of
  notation", *Proc. R. Soc. Lond. A* **79** (1907): 182–193.
- Frisch & Waugh, *Econometrica* 1(4) (1933): 387–401. Lovell, *JASA* 58(304) (1963): 993–1010.
- History and priority: D. Basu, "The Yule-Frisch-Waugh-Lovell Theorem", arXiv:2307.00369.
  https://arxiv.org/abs/2307.00369 — *"the coefficients on any subset of covariates in a multiple
  regression is equal to the coefficients in a regression of the residualized outcome variable on the
  residualized subset of covariates"*; argues the theorem be renamed *"to recognize the pioneering
  contribution of the statistician G. Udny Yule"*.

**Covers:** that residualising two variables on a shared regressor set yields exactly their partial
correlation; hence that a zero raw correlation maps to a *non*-zero partial correlation whenever both
variables load on the design. Standard textbook treatment (e.g. Greene, *Econometrics I*, "Partial
Regression and Correlation"; Cohen & Cohen on suppression).
**Does not cover:** any applied warning, any magnitude, anything cross-modal.

### (2) Winkler, Renaud, Smith & Nichols (2020) — *the most threatening work overall: cross-modal, and it warns AND fixes*

"Permutation inference for canonical correlation analysis." *NeuroImage* **220**:117065.
DOI 10.1016/j.neuroimage.2020.117065 · arXiv:2002.10046 · PMC7573815
https://doi.org/10.1016/j.neuroimage.2020.117065

Abstract, verbatim: *"Canonical correlation analysis (CCA) has become a key tool for population
neuroimaging, allowing investigation of associations between many imaging and non-imaging measurements.
As other variables are often a source of variability not of direct interest, previous work has used CCA
on residuals from a model that removes these effects, then proceeded directly to permutation inference.
**We show that such a simple permutation test leads to inflated error rates. The reason is that
residualisation introduces dependencies among the observations that violate the exchangeability
assumption.**"*

**Covers:** the paper's claim (b) in full — residualising against a shared nuisance model, in an
explicitly **cross-modal** analysis (imaging block × non-imaging block), induces dependence that was not
there, and this **inflates the resulting cross-modal association statistics**. It is quantified as
Type-I error inflation, and they supply a remedy (project residuals onto a lower-dimensional exchangeable
basis). A reviewer holding this paper will say the warning was issued in 2020 and already solved.
**Does not cover:** the induced dependence is characterised **among observations** (exchangeability),
not as a correlation between two signal directions; no induced-|r| magnitude is reported.

### (3) Alfaro-Almagro, McCarthy, Afyouni, Andersson, Bastiani, Miller, Nichols & Smith (2021)

"Confound modelling in UK Biobank brain imaging." *NeuroImage* **224**:117002.
DOI 10.1016/j.neuroimage.2020.117002 · PMID 32502668 · PMC7610719

Introduction, verbatim: *"**spurious associations can be induced between pairs of otherwise independent
variables** if the unconfounding is not carried out correctly (e.g., if the confounds were not demeaned
first)."* Abstract: *"the resulting high statistical power also raises the sensitivity to confound
effects."* Cross-modal by construction (imaging-derived phenotypes × lifestyle/health/body variables);
they report that Bonferroni-significant cross-modal correlations move from 105,122 (simple confounds) to
53,995 (full confound model).

**Covers:** the sentence "adjustment induces association between otherwise-independent variables", in
print, in a cross-modal imaging-vs-non-imaging setting, in the field's reference deconfounding paper.
**Does not cover — and this is the one real escape hatch:** their statement is **conditioned on
incorrect unconfounding** ("if the confounds were not demeaned first"). The paper's claim is stronger —
that *correctly performed* residualisation of *exactly orthogonal* signals still yields a non-zero
residual correlation. That distinction must be made explicitly and loudly or this paper reads as a scoop.
No numeric floor (they report %VE and −log10 P).

### (4) Murphy, Birn, Handwerker, Jones & Bandettini (2009) — *the strongest quantified mechanism*

"The impact of global signal regression on resting state correlations: are anti-correlated networks
introduced?" *NeuroImage* **44**(3): 893–905. DOI 10.1016/j.neuroimage.2008.09.036 · PMID 18976716
https://doi.org/10.1016/j.neuroimage.2008.09.036

Abstract, verbatim: *"This global signal regression method has been shown to introduce negative
activation measures in standard fMRI analyses. … Here we show that, after global signal regression,
**correlation values to a seed voxel must sum to a negative value**. Simulations also show that small
phase differences between regions can lead to **spurious negative correlation values**. … These results
**call into question the interpretation** of negatively correlated regions in the brain when using global
signal regression as an initial processing step."*

**Covers, precisely, the structure of the paper's §1 claim:** regressing a *shared* nuisance out of two
signals mathematically forces a non-zero (here negative) correlation between the residuals; it is
quantified (a summation constraint, plus simulated magnitudes); and it is presented explicitly as a
**warning that the artefact contaminates a downstream correlation measurement in that literature**.
That is claim (b) in the brief — an applied, quantified inflation warning — already delivered in 2009.
Reinforced by: **Saad et al.**, *Brain Connect.* 2(1) (2012): 25–32, DOI 10.1089/brain.2012.0080 — GSR
*"can alter local and long-range correlations, potentially spreading underlying group differences to
regions that may never have had any"*, demonstrated on 3-ROI simulations with **no** true group
difference; **Fox et al.**, *J. Neurophysiol.* 101(6) (2009): 3270–3283, DOI 10.1152/jn.90777.2008 —
*"global signal regression mandates negative correlations at the single subject level"*; Anderson et al.
(2010); **Murphy & Fox**, *NeuroImage* 154 (2017): 169–173 — *"The mathematics of global signal
regression (GSR) mandate that functional connectivity analyses performed using this processing step show
both positive and negative values that average to zero."*

**Does not cover:** cross-modality — GSR is voxel-vs-voxel within one modality; and the nuisance is a
single regressor, not a 99-column categorical design.

### (5) Li, Zhang, Patil & Johnson (2023) — *the closed form, in genomics*

"Overcoming the impacts of two-step batch effect correction on gene expression estimation and inference."
*Biostatistics* **24**(3): 635–652. DOI 10.1093/biostatistics/kxab039 · PMC10449015
https://academic.oup.com/biostatistics/article/24/3/635/6459158

Verbatim: *"removing batch effects with two-step methods (such as ComBat) **introduces a correlation
structure into the adjusted data**"*; *"the adjusted data points within each batch are correlated with
each other, because they are functions of all the other data from the batch."*
They derive the induced correlation matrix as **`Σ = (I − X(XᵀX)⁻¹Xᵀ)σ²`** — the residual-maker matrix,
i.e. the exact object the paper's floor is measuring. Quantified downstream: FPR 18.3% vs 5% nominal.

**Covers:** publication of the identical projection operator as *the* source of adjustment-induced
correlation, with a false-positive-rate consequence.
**Does not cover:** the induced correlation here is **between samples**, not between two signals — the
dual of the paper's statement. Not cross-modal. No numeric correlation floor.

### (6) Dahl, Guillemot, Mefford, Aschard & Zaitlen (2019) — *owns the verb*

"Adjusting for principal components of molecular phenotypes induces replicating false positives."
*Genetics* **211**(4): 1179–1189. DOI 10.1534/genetics.118.301768 · PMID 30692194
https://academic.oup.com/genetics/article/211/4/1179/5931516

Verbatim: *"The key problem is that the effect of x on Y leads PCs to partially capture x, **analogous to
unshielded colliders** in a directed graphical model. That is, **conditioning on genomic PCs can cause,
rather than remove, bias**."* · *"conditioning on the collider u₁ **induces spurious correlation** with
all other y_p."*
Analytic bias approximation `Bias_p ≈ −a·V_p1·V_q1`; >10-fold FPR inflation; *"nearly all false positives
replicating."*

**Covers:** adjustment-induces-association, in biology, with an analytic magnitude, framed as a warning.
**Does not cover:** cross-modal; the confound is estimated PCs, not an observed design.

### (7) The "everyone already knows" tier — most damaging to a discovery framing

- **Nygaard, Rødland & Hovig (2016)**, *Biostatistics* **17**(1): 29–39, DOI 10.1093/biostatistics/kxv027,
  PMID 26272994. Verbatim: *"batch-adjusted values are no longer independent … it may also **induce new
  dependencies**"*; abstract: *"this approach may systematically induce incorrect group differences in
  downstream analyses … **The scientific community seems to be largely unaware of how this approach may
  lead to false discoveries.**"* (Note: the "community is unaware" line was already claimed in 2016.)
- **Gordon Smyth**, Bioconductor support #133791, 8 Sep 2020, on `limma::removeBatchEffect`:
  *"the batch correction has **introduced correlations**"*; *"There is no algorithm that can remake the
  data entirely as if the batch effects never existed."* https://support.bioconductor.org/p/133791/
  The `removeBatchEffect` man page carries the standing warning that it is *"not intended to be used prior
  to linear modelling."*
- **Aschard, Vilhjálmsson, Joshi, Price & Kraft (2015)**, "Adjusting for heritable covariates can bias
  effect estimates in genome-wide association studies", *AJHG* **96**(2): 329–339,
  DOI 10.1016/j.ajhg.2014.12.021, PMID 25640676: *"an unintended bias is introduced with respect to the
  primary outcome as a result of the adjustment, and this bias might lead to false positives."*
- **Zindler et al. (2020)**, *BMC Bioinformatics* **21**:271, DOI 10.1186/s12859-020-03559-6, PMC7328269:
  ComBat produces false positives on *random* data *"in unbalanced **as well as in balanced** sample
  distributions"* and *"**Increasing the number of corrected factors led to an exponential increase** in
  the number of false positive signals."* — i.e. the "99 columns" dependency is already published.
- **Micheletti, Schlauch, Quackenbush & Ben Guebila (2024)**, "Higher-order correction of persistent batch
  effects in correlation networks", *Bioinformatics* **40**(9): btae531, DOI 10.1093/bioinformatics/btae531,
  PMC11441315: *"artifactual DC can skew the correlation structure, leading to the identification of false,
  non-biological associations, **even when the input data are corrected using standard batch correction**."*
  This is the closest published framing of *a spurious-correlation floor surviving confound adjustment* —
  but gene–gene, within-modality, and no absolute magnitude.
- **Collider/selection-bias canon**, for the "textbook" citation the brief asked for bluntly:
  Hernán, Hernández-Díaz & Robins, "A structural approach to selection bias", *Epidemiology* 15(5) (2004):
  615–625; Cole et al., "Illustrating bias due to conditioning on a collider", *IJE* 39(2) (2010): 417–420,
  PMID 19926667; Elwert & Winship, "Endogenous selection bias: the problem of conditioning on a collider
  variable", *Annu. Rev. Sociol.* 40 (2014): 31–53, PMID 30111904.
- **Hamdan et al. (2023)**, "Confound-leakage: confound removal in machine learning leads to leakage",
  *GigaScience* 12:giad071, DOI 10.1093/gigascience/giad071: *"this common approach to confound removal
  biases ML models … this common deconfounding approach can leak information such that what are null or
  moderate effects become amplified to near-perfect prediction."*
- **Chyzhyk, Varoquaux, Milham & Thirion (2022)**, *GigaScience* 11:giac014,
  DOI 10.1093/gigascience/giac014: *"conditioning on [a confound] can reverse the correlation between 2
  variables … a phenomenon known as the Berkson or Simpson statistical paradox"*; empirical
  ρ = −0.67 → −0.07 (CamCan). **The explicit collider/Berkson framing the brief asked about, already
  applied to confound-adjusted brain analyses.**
- **Snoek, Miletić & Scholte (2019)**, "How to control for confounds in decoding analyses of neuroimaging
  data", *NeuroImage* 184: 741–760, DOI 10.1016/j.neuroimage.2018.09.074: most confound-control methods
  "yield biased results"; whole-data confound regression drives below-chance performance.
- **Pain, Dudbridge & Ronald (2018)**, "Are your covariates under control? How normalization can
  re-introduce covariate effects", *Eur. J. Hum. Genet.* 26(8): 1194–1201, DOI 10.1038/s41431-018-0159-6:
  *"applying rank-based INT to the dependent variable residuals after regressing out covariates
  re-introduces a linear correlation between the dependent variable and covariates, increasing type-I
  errors and reducing power."*
- **Yu, Zoh, Fluharty et al. (2024)**, "Misstatements, misperceptions, and mistakes in controlling for
  covariates in observational research", *eLife* 12:e82268, DOI 10.7554/eLife.82268 — "Misperception 8"
  is residualisation: *"this procedure can introduce bias in the magnitude of the coefficients (effect
  sizes) … and not just the test of statistical significance."* Clean, citable, general-statistics
  "residualisation biases effect sizes" reference. See also **Darlington & Smulders (2001)**, "Problems
  with residual analysis", *Anim. Behav.* 62: 599–602, DOI 10.1006/anbe.2001.1806.
- **Nalci, Rao & Liu (2019)**, *NeuroImage* 184: 1005–1031 and 202:116005 — *"nuisance regression does not
  necessarily eliminate the relationship between DFC estimates and nuisance norms."*
- **Spisak (2022)**, "Statistical quantification of confounding bias in machine learning models",
  *GigaScience* 11:giac082 — *"state-of-the-art confound mitigation approaches may fail preventing
  confounder bias in several cases"* (N = 1,865).

## 1.1b Pre-empting the two attacks that will actually be made

**Attack A — "your floor is just the partial-correlation null, i.e. Fisher's df correction."** A reviewer
will observe that residualising on a rank-q design leaves df = n − q − 2, so a partial correlation has null
sd ≈ 1/√(n−q−1), and will claim 0.067–0.140 is a restatement of that. **This attack fails and we can kill
it with a number.** Simulated: two *independent random* vectors residualised on the same 99-column design
at n = 2,530 give induced `|r|` of mean **0.0031**, p95 0.0069, max 0.0106 — 20× below the reported floor,
and consistent with 1/√(2530−99) ≈ 0.020 being an *upper* bound on the df contribution. The reported
0.067–0.140 is therefore **structural** — driven by `R_u R_v ρ`, the overlap of the signal pair with the
design span — exactly as §1 of the readout states. **Include this control in the manuscript.** It is the
single most useful defensive asset the search produced.

**Attack B — "this is one line of algebra."** It is (§1.0). Concede immediately and in the text; do not
let a reviewer be the one to point it out.

## 1.2 What is genuinely unoccupied

Narrow, and not headline material. **And note first what is *not* unoccupied:** the framing "we compute
an empirical noise floor for image-vs-molecular correlations" is an established sub-industry with a
different mechanism (spatial autocorrelation rather than a confound design). Do **not** present the floor
as the first cross-modal noise floor; present it as an *additional, orthogonal* floor source:

- **Markello & Misic (2021)**, "Comparing spatial null models for brain maps", *NeuroImage* 236:118052,
  DOI 10.1016/j.neuroimage.2021.118052 — *"naive null models that do not preserve spatial autocorrelation
  consistently yield elevated false positive rates and unrealistically liberal statistical estimates."*
- **Burt, Helmer, Shinn, Anticevic & Murray (2020)**, *NeuroImage* 220:117038; **Fulcher, Arnatkeviciute &
  Fornito (2021)**, *Nat. Commun.* 12, DOI 10.1038/s41467-021-22862-1 (>500-fold elevated FPRs in
  gene-category enrichment against brain maps); **Cao et al. (2024)**, *NeuroImage*,
  DOI 10.1016/j.neuroimage.2024.120622; **Bazinet, Liu & Misic (2025)**, *Imaging Neurosci.*,
  DOI 10.1162/imag.a.118.
- **Marek et al. (2022)**, "Reproducible brain-wide association studies require thousands of individuals",
  *Nature* 603: 654–660, DOI 10.1038/s41586-022-04492-9 — canonical cross-modal effect-size floor,
  median |r| = 0.02–0.03 at N ≈ 50,000. **Cite this ourselves, as ammunition:** a confound-induced floor
  of 0.067–0.140 would *swamp* typical real cross-modal effects. Letting a reviewer introduce it is worse.
- **Helmer et al. (2024)**, "On the stability of canonical correlation analysis and partial least squares
  with application to brain-behavior associations", *Commun. Biol.*, DOI 10.1038/s42003-024-05869-4 —
  direct competitor for "we quantify a floor for cross-modal analyses".

With that conceded, what remains:

1. No work states the induced-correlation floor for a **cross-modal** (image-embedding × molecular)
   confound-adjusted analysis. The closest is **Nersisyan, Loher & Rigoutsos (2025)**, CorrAdjust,
   *NAR* **53**(10): gkaf444, DOI 10.1093/nar/gkaf444 — genuinely cross-modal (miRNA–mRNA, 25,063
   TCGA/GTEx/Geuvadis datasets) and residualisation-based, but argues the *opposite* direction
   (residualisation *removes* spurious correlation: r 0.83 → 0.08) and never considers an induced floor.
2. No work reports a magnitude for a **cancer-type + tissue-source-site design at TCGA scale**.
2b. No work states that ***correctly performed*** residualisation (properly demeaned, properly fitted) of
   ***exactly orthogonal*** signals still yields non-zero residual correlation as a stated geometric fact.
   Alfaro-Almagro 2021 conditions on *incorrect* unconfounding; Winkler 2020 frames it at the level of
   observation exchangeability; the GSR papers derive it for a *single* shared regressor within one
   modality. **This is the sharpest surviving distinction and every claim should be phrased around it.**
3. In computational pathology specifically the gap is real but shallow — the field does site-adjusted
   cross-modal analysis without a floor. **Howard et al. (2021)**, *Nat. Commun.* **12**:4423,
   DOI 10.1038/s41467-021-24698-1, PMC8292530, establishes TCGA site signatures and proposes
   preserved-site cross-validation, but does not residualise-then-correlate and reports no floor.
   **Murchan et al. (2024)**, *J. Pathol. Inform.* **15**:100396, DOI 10.1016/j.jpi.2024.100396,
   PMC11470259, ComBat-harmonises WSI deep features against TSS and *then* predicts molecular features
   (MSI, BRAF, TP53) — i.e. exactly a confound-adjusted cross-modal analysis — with **no consideration of
   an induced floor**. Also **Dehkharghanian et al. (2021/2023)**, "Biased data, biased AI",
   DOI 10.21203/rs.3.rs-943804/v1.

## 1.3 Defensible wording

**Do not write:** "we report that residualisation induces spurious correlation", "nobody reports this",
"a real finding", or anything implying discovery.

**Write instead (safe):**

> Residualising two signals against a shared confound design does not leave them uncorrelated. For
> signals constructed orthogonal in-sample, the residual correlation is exactly the partial correlation
> `−R_u R_v ρ / sqrt((1−R_u²)(1−R_v²))` — a direct consequence of the Yule–Frisch–Waugh–Lovell identity
> [Yule 1907; Frisch & Waugh 1933; Lovell 1963], the same mechanism by which global-signal regression
> forces spurious anticorrelation in resting-state fMRI [Murphy et al. 2009; Saad et al. 2012], by which
> two-step batch correction imposes a residual-maker correlation structure on adjusted expression data
> [Li et al. 2023; Nygaard et al. 2016], and the reason simple permutation tests on residualised
> cross-modal CCA are anticonservative [Winkler et al. 2020]. **Unlike the induced associations described
> for imaging-vs-non-imaging analyses, which arise when unconfounding is performed incorrectly
> [Alfaro-Almagro et al. 2021], this floor is present under correct residualisation of exactly orthogonal
> signals, and its magnitude in a cross-modal design has not to our knowledge been reported.** For the
> 99-column cancer-type + tissue-source-site design at n = 2,530 used here, a signal pair drawn orthogonal
> to the image score reads 0.067–0.140 rather than 0, varying with how much of the pair lies in the design
> span; the same design applied to independent random signals yields only 0.003, confirming the floor is
> structural rather than a degrees-of-freedom artefact. This exceeds typical reported cross-modal effect
> sizes [Marek et al. 2022], so confound-adjusted image×molecular correlations in this range are not
> interpretable without a matched null. We therefore read our floor with a paired test.

**Also safe, and stronger as a contribution framing:** present it not as a phenomenon but as an
*operational requirement* — that cross-modal confound-adjusted analyses need a per-draw paired baseline,
which is a protocol claim, not a discovery claim. Cite Micheletti et al. 2024 as the nearest precedent
(a floor surviving correction, within-modality) and Murchan et al. 2024 as the field practice that omits it.
Position the floor as **orthogonal and additive to** the spatial-autocorrelation null family
[Markello & Misic 2021; Burt et al. 2020; Fulcher et al. 2021], not as a replacement for or predecessor of it.

**Budget:** at most one paragraph plus one supplementary figure. This is a methods caveat, not the
headline. If the manuscript's headline currently rests on §1 of `PHASE1B_TARGETED_READOUT.md`, it needs a
different headline.

---

# Q2 — spike-in certification of a detection floor

## Verdict: **PARTIALLY ANTICIPATED**

Every component has strong biological precedent, and two of the four co-occur in single papers. Only the
**cross-modal conjunction** is unoccupied. The claim survives, badly narrowed, and only if the cited works
below are cited pre-emptively.

Component audit of the claim:

| component | precedent | status |
|---|---|---|
| inject known-strength signal | ERCC, Munro 2014, Gerard 2020 | **occupied** |
| certify a **detection floor** by injection, in biology | Munro 2014 (LODR), Jiang 2011 | **occupied** |
| inject into **real** data to certify a **confound-adjusted** pipeline | Gerard 2020 (seqgendiff) | **occupied** |
| injected quantity is a **correlation on a direction pair**, **cross-modal** | — | **open** |
| report **transmission** *and* **detection** floors as a pair | astronomy injection-recovery completeness (concept), no biology instance | **weakly open** |

## 2.1 Strongest prior art

### (1) Munro et al. (2014) — *the most threatening; a named spike-in-certified detection floor in biology*

"Assessing technical performance in differential gene expression experiments with external spike-in RNA
control ratio mixtures." *Nature Communications* **5**:5125. DOI 10.1038/ncomms6125
https://www.nature.com/articles/ncomms6125

Verbatim: *"These control ratio mixtures with defined abundance ratios enable assessment of diagnostic
performance of differentially expressed transcript lists, **limit of detection of ratio (LODR)
estimates** and expression ratio variability and measurement bias. … An interlaboratory study using
identical samples shared among 12 laboratories with three different measurement processes demonstrates
generally consistent diagnostic power across 11 laboratories."*
Tooling: Lund, Pine, Salit & Munro, "The erccdashboard", *J. Biomol. Tech.* 25 (2014), PMC4162230.

**Covers:** known-strength signal injected into real biological samples; a named, published
**detection floor** (LODR = minimum reliably-detectable ratio); estimated *across 12 labs and 3
measurement processes*, i.e. under site/protocol nuisance variation.
**Does not cover:** no residualisation against a confound design matrix; single modality; the injected
quantity is an abundance ratio, not a correlation.
**Kills outright** the phrasing "nobody has used an injected spike-in to certify a detection floor in biology."

### (2) Gerard (2020), `seqgendiff` — *the injection-into-real-data-to-benchmark-confound-adjustment paper*

"Data-based RNA-seq simulations by binomial thinning." *BMC Bioinformatics* **21**(1):206.
DOI 10.1186/s12859-020-3450-9 · PMID 32448189 · docs
https://dcgerard.github.io/seqgendiff/reference/thin_diff.html

Verbatim (abstract): *"Rather than generate data from a theoretical model, in this paper we develop
methods to **add signal to real RNA-seq datasets**."*
Verbatim (`thin_diff`): *"Given a matrix of real RNA-seq counts, this function will add a **known amount
of signal** to the count matrix. … The user may also **control for the amount of correlation between the
observed covariates and any unobserved surrogate variables**."*

**Covers:** injection of known-strength signal into real biological data, *expressly designed* to
benchmark confound-adjustment methods (SVA/RUV/CATE) — components (i) + (iii) + (v) in one paper.
**Does not cover:** no named floor (power/FDR/AUC curves only); not cross-modal; the injected object is a
differential-expression effect, not a correlation on a paired direction.
**This is the paper that must be cited.** It owns "inject to certify a confound-adjusted pipeline."

### (3) Jiang et al. (2011) — the ERCC foundation

"Synthetic spike-in standards for RNA-seq experiments." *Genome Research* **21**(9): 1543–1551.
DOI 10.1101/gr.121095.111 · PMC3166838
Verbatim: *"We used a newly developed pool of 96 synthetic RNAs … as spike-in controls to **measure
sensitivity**, accuracy, and biases in RNA-seq experiments … By using the control RNAs, we **derive limits
for the discovery and detection of rare transcripts** in RNA-seq experiments."*
Also: External RNA Controls Consortium, *BMC Genomics* 6:150 (2005); SEQC/MAQC-III, *Nat. Biotechnol.* 32
(2014): 903–914.
**Covers:** the canonical "spike-ins → detection limit" idiom the paper's own code header already invokes
(`v2/calibra/calibration.py:10`, *"borrowed in spirit from ERCC spike-ins"*). Keep that acknowledgement.

### (4) Injection-recovery completeness — *the unnamed structural precedent*

The paper's **transmission floor** (fraction of injected signal surviving the pipeline, as a function of
injected strength) is formally an **injection–recovery completeness curve**, a standard instrument-
characterisation method outside biology:

- Biwer et al., "Validating gravitational-wave detections: the Advanced LIGO hardware injection system",
  *Phys. Rev. D* **95**:062002 (2017), DOI 10.1103/PhysRevD.95.062002, arXiv:1612.07864 — signals of known
  parameters physically injected; sensitivity characterised *"by looking for discrepancies between the
  injected and recovered signals."*
- Kepler DR25 flux-level transit injection tests (Christiansen et al., NASA NTRS 20170009551); eROSITA
  source-detection simulations, *A&A* 661 (2022), A3.

**Not biology, not confound-adjusted, not cross-modal** — but a reviewer at a methods venue will know
these. Failing to name the connection is the most likely objection; naming it costs nothing and
inoculates the paper.

### (5) Confirmed non-threats — the distinctions hold

- **Simulation-based calibration**: Talts, Betancourt, Simpson, Vehtari & Gelman, arXiv:1804.06788;
  Cook, Gelman & Rubin, *JCGS* 15(3) (2006): 675–692. The measured object is **uniformity of rank
  statistics / posterior calibration**, with no injection into real data and no detection floor. The
  paper's distinction from SBC is correct and safe to assert.
- **Negative-control methods** (RUV): Gagnon-Bartsch & Speed, *Biostatistics* 13(3) (2012): 539–552;
  Gerard & Stephens, *Statistica Sinica* (2021), DOI 10.5705/ss.202018.0345. These use controls known
  **not** to change — a *negative* control lineage. The positive-injection distinction is genuine and
  worth stating explicitly.
- **scRNA-seq integration benchmarks**: Luecken et al., *Nat. Methods* 19 (2022): 41–50,
  DOI 10.1038/s41592-021-01336-8; Tran et al., *Genome Biol.* 21:12 (2020). Score bio-conservation vs.
  batch removal on real labelled data — **no injected known-strength signal, no floor**.
- **Chyzhyk, Varoquaux, Milham & Thirion (2022)**, "How to remove or control confounds in predictive
  models, with applications to brain biomarkers", *GigaScience* 11:giac014, DOI 10.1093/gigascience/giac014.
  States the *concern* — *"using a more powerful predictive model … may remove signal of interest,
  unrelated to the confound"* — but no injection, no floor, not cross-modal. **Best citation for the
  motivation of the transmission floor**: frame the contribution as quantifying what this paper only warns about.
- **Computational pathology / radiogenomics**: nothing found that injects a synthetic signal to test
  detection limits in histopathology-embedding-vs-molecular analyses. HE2RNA, hist2RNA, THItoGene, GHIST
  are predictive-performance oriented. Howard et al. (2021) *does* run a synthetic manipulation —
  *"We varied the ER negativity of the 23 target slides from 0 to 100% … we applied an artificial staining
  artifact to 0–100% of the target slides"* — but to **demonstrate confounding**, not to certify a
  detection floor after adjustment. This is the paper's **strongest empty quadrant** and it is worth
  noting that Howard's manipulation is the nearest thing in the domain.
- **Two independent sweeps** (the spike-in/calibration silo and the neuroimaging silo) separately failed
  to find any paper that injects a signal of known strength and measures how much survives *nuisance /
  confound regression specifically*. A Europe PMC query for `("simulated signal" OR "synthetic signal" OR
  "known ground truth") AND ("nuisance regression" OR "confound") AND "recover"` returned **0 hits**.
  This is the strongest positive evidence available for Q2's residual novelty, and it is weak evidence —
  absence from two sweeps, not proof of absence.

## 2.2 Defensible wording

**Do not write:** "the first use of a spike-in to certify a detection floor", "no prior work injects a
known-strength signal to calibrate a confound-adjusted analysis", or any unqualified priority claim.

**Write instead (safe):**

> Spike-in certification of detection limits is long established within a single modality — ERCC controls
> yield sensitivity limits for RNA-seq [Jiang et al. 2011] and a named limit of detection of ratio
> [Munro et al. 2014] — and injecting known-strength signal into *real* data to benchmark
> confound-adjustment methods is established for differential expression [Gerard 2020]. Adapting that
> logic, we inject a correlation of known strength onto a named direction pair spanning **two modalities**
> and push it through the identical residualisation, reporting the fraction transmitted
> (an injection–recovery curve, in the sense used for pipeline characterisation in
> [Biwer et al. 2017]) and the smallest true correlation reliably recovered. To our knowledge this
> paired transmission/detection reporting has not previously been applied to a confound-adjusted
> cross-modal biological analysis, where the field's practice [e.g. Murchan et al. 2024] adjusts image
> features against tissue-source site and correlates with molecular endpoints without characterising
> either floor.

Note the load-bearing qualifiers: *within a single modality* → *two modalities*; *"to our knowledge"*;
*"has not previously been applied to"* rather than *"is the first"*. Keep the existing acknowledgement in
`v2/calibra/calibration.py:10` and promote it into the manuscript.

---

## Residual coverage gaps

Two probes went unrun (WebSearch budget exhausted at 200/200; sweep completed via Europe PMC REST and
publisher pages):

1. **fMRI papers that add synthetic activation to real resting-state data and measure loss specifically
   through nuisance regression.** A Europe PMC query for
   `("simulated signal" OR "synthetic signal" OR "known ground truth") AND ("nuisance regression" OR
   "confound") AND "recover"` returned **0 hits**, which is weak evidence of absence, but neuRosim
   (Welvaert & Rosseel, *J. Stat. Softw.* 44(10), 2011) and the Eklund et al. (2016) lineage should be
   checked by hand before relying on the negative.
2. **DNA-methylation array positive-control probes for limit of detection after cell-composition
   adjustment.**
3. **Bright, Tench & Murphy (2017)**, *NeuroImage* (PMID 28025128), on degrees of freedom in nuisance
   regression — flagged as a plausible source of a df-based quantified correlation floor;
   **unverified** (ScienceDirect 403). Worth a second pass given Attack A in §1.1b.
4. **Dinga et al. (2020)**, bioRxiv, "Controlling for effects of confounding variables on machine
   learning predictions" — **unverified** (bioRxiv full text 403).
5. **Winkler et al. 2020 FPR tables** — abstract verified verbatim; the numeric inflation magnitudes were
   not extracted (PDF extraction failed). If those tables contain an induced-correlation magnitude rather
   than only error rates, Q1's residual novelty (a number) disappears entirely. **Check this first.**

None of 1–4, if found, would change the Q1 verdict. Item 5 could remove even the numeric residue.
Items 1–2 could push Q2 from PARTIALLY ANTICIPATED toward ALREADY REPORTED.

## Artefacts

- Closed-form verification: scratchpad `sim3.py` (max abs error 7.4 × 10⁻¹⁶ over 300 draws).
- Noise-only control: scratchpad `sim.py` (independent vectors, same design → induced `|r|` mean 0.0031).
  These should be reproduced into the repo if the identity is quoted in the manuscript.
