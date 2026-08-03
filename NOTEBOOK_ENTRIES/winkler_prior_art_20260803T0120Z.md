## 2026-08-03 01:20 UTC — Winkler 2020 does NOT report an induced-correlation magnitude; Track 2's magnitude framing survives, but the classical N−R result must be cited

**Logged:** 2026-08-03 01:20 UTC. **How obtained:** full-text retrieval of PMC7573815 (open access) and PMC4010955, five targeted passes over Theory §2.6, Simulations, Results §4.3/4.5, Discussion §5.3, every figure caption and every table header; plus a secondary prior-art sweep.

### Technical

**Verdict: NO.** Winkler AM, Renaud O, Smith SM, Nichols TE, "Permutation inference for canonical
correlation analysis", *NeuroImage* 2020;220:117065 reports **inflated error rates, power and
p-value miscalibration — never a correlation-unit magnitude** for the dependence residualisation
induces between two otherwise-independent blocks.

What the paper actually reports, by exhibit:

| Exhibit | Units |
|---|---|
| Fig. 1 | conceptual diagram, no axes |
| Fig. 2 | x = canonical correlation, **y = p-value** (calibration plot; the scenario shown has *no* nuisance variables) |
| Tables 1–3 | taxonomy (partial / part / bipartial CCA), pseudocode, Theil vs Huh–Jhun — no data |
| Table 6 | per-comparison error rate, **%** |
| Table 7 | pcer **%** (simple residualisation 83.85% [82.17–85.40] vs ~5% for Huh–Jhun/Theil) |
| Table 8 | pcer **%** across distributions |
| Table 9 | **FWER %** vs N = 100…1000 (simple residuals 96.60% at N=100 → 7.70% at N=1000) |
| Table 10 | observed power, **%** |

No figure or table anywhere plots or tabulates a correlation value obtained from null/independent
data after residualisation. No numeric correlation coefficient appears in the narrative at all —
canonical correlations are referred to only symbolically (r_k, ρ_k).

The closest the paper comes, verbatim from §2.6 (Nuisance variables):

> "While **Y** occupies an N-dimensional space, **Ỹ** occupies a smaller one; its dimensions are, at
> most, of a size given by the rank of R_Z, which is N−R assuming **Y** and **Z** are of full rank."

> "With fewer effective observations determined by this lower space after residualisation, and the
> same number of variables, the sample canonical correlations in the unpermuted case are
> **stochastically larger** than in the permuted, which in turn leads to an excess of spuriously
> small p-values."

That is a stochastic-ordering statement with no magnitude, and it is framed *relative to the
permutation distribution* (an exchangeability failure), not as "residualising two orthogonal sets
induces correlation ρ".

Discussion §5.3, the paper's own summary of its headline finding, is entirely in error-rate terms:

> "Merely regressing out such nuisance variables from all other variables that are subjected to cca,
> then proceeding to a simple permutation test, leads to inflated error rates and an invalid test…
> This inflated error rate, even after multiple testing correction, is the probably the most striking
> finding of the current study… particularly if the number of nuisance variables is relatively large
> compared to the sample size, as shown in Section 4.5."

The R/N dependence is therefore asserted **qualitatively** and demonstrated only through error rates.
Table 9 / §4.5 is the nearest analogue to our n-sweep (N varied at fixed R) but its response variable
is FWER in percent.

The only place R enters a formula is the Bartlett/Wilks degrees-of-freedom correction —
`λ_k = −(N−C−P+Q+3/2) ln(∏(1−r_i²))` with "C=R for partial or part, and C=max(R,S) for bipartial
cca". That encodes the N−R effective-sample-size adjustment in order to *undo* the effect for a
parametric p-value; it is never inverted into a statement of induced magnitude.

**Winkler et al. 2014** (GLM permutation paper, PMC4010955) takes the same posture — "residuals
induce dependence and any ee or ise assumptions on ϵ will not be conveyed to ϵ̂_Z" — and reports Type
I error rates and power percentages only. Its cited antecedents (Anderson & Robinson 2001; Kennedy
1995; Huh & Jhun 2001) are all permutation-validity literature, none cited for a magnitude.

**Strongest genuine prior art is textbook, not a paper, and we should pre-empt it ourselves:** the
classical result that under H₀ the canonical correlations of rank-R-residualised data follow the null
CCA distribution with N replaced by an effective N−R (Muirhead 1982; Anderson, *Introduction to
Multivariate Statistical Analysis*, 3rd ed. 2003), which implies E[r²] ≈ 1/(N−R−1) in the univariate
case. Winkler 2020 *uses* exactly this fact (C=R in the Wilks statistic) without converting it into a
correlation magnitude.

Caveat on the sweep: the agent's WebSearch budget was exhausted, so the secondary prior-art sweep ran
through a weaker HTML-search path and is **low confidence**. The residualisation-methodology hits it
did surface (García García et al. 2020, *J. Applied Statistics*; arXiv:2410.17680) are about
multicollinearity and coefficient interpretation within a single regression, not induced correlation
between two independent blocks. A proper Scholar/PubMed sweep should be run before the novelty
sentence is finalised.

### In plain terms

We were worried that the one thing we might still own — an actual *number* for how much correlation
gets manufactured when you statistically remove the same nuisance variables from two unrelated
measurements — had already been published by the most relevant paper in the field. It has not.
Winkler and colleagues showed that doing this **breaks your statistics** (their false-positive rate hit
96.6% in one setting), and they showed how to fix it. They never said how big the manufactured
correlation itself is. That number is still ours to report.

But the underlying mathematics is classical and we must say so: statisticians have known since at
least Muirhead (1982) that residualising away R nuisance variables leaves you behaving as though you
had N−R subjects instead of N. Winkler's own formula uses that fact. Our contribution is the size of
the effect on real cross-modal data at realistic R and N — not the discovery that the effect exists.

### Meaning for the claim

* Track 2 may proceed with the magnitude framing. The claim must be worded as **"we quantify, in
  correlation units, the effect whose inferential consequence Winkler et al. (2020) characterised"** —
  never as the discovery of the phenomenon.
* The write-up must cite (a) Yule 1907 / Frisch–Waugh–Lovell for the identity, (b) Winkler et al.
  2020 for the inferential consequence and the fix, and (c) Muirhead 1982 / Anderson 2003 for the
  N−R effective-sample-size result, so that no reviewer can present any of the three as an omission.
* This makes the *structural vs degrees-of-freedom* distinction load-bearing rather than decorative.
  If our 0.067–0.140 were just the N−R sampling scale, the classical result already covers it
  completely and we would have nothing. Predeclared falsifier P5 (a matched-rank structureless design
  must induce ≈0) is therefore the test that decides whether Track 2 has any content at all.

### Files / commits

`v2/research/rebase/nature/P1_PREDECLARATION.md` (§C); to be folded into
`v2/research/rebase/nature/TRACK2_INDUCED_CORRELATION.md`.
Sources: [PMC7573815](https://pmc.ncbi.nlm.nih.gov/articles/PMC7573815/),
[PMC4010955](https://pmc.ncbi.nlm.nih.gov/articles/PMC4010955/).
