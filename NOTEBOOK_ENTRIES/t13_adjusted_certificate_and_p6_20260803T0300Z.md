## 2026-08-03 03:00 UTC — The adjustment fully discharges the site confound (joint accuracy 0.363 → 0.012, at a chance rate of 0.0118); and P6 estimator-robustness FAILS as predeclared, on the Ridge shrinkage axis only

**Logged:** 2026-08-03 03:00 UTC. **How obtained:** `python -m morpheus.v2.calibra.confound_certificate --artifacts d2_h_seed42.npz d2_i_seed42.npz --partition test --residualise --n-permutations 1000 --n-boot 200 --n-boot-axes 8 --n-jobs 6`; and `induced_correlation_sweep --n-splits-grid 2,5,10,20 --alpha-grid 0.01,1.0,100.0 --designs cancer_tss_pool10,cancer,gaussian_k99 --n-grid 2530,6427`, both on the Lambda box `~/ws_p1`.

### Technical

#### Part 1 — the adjusted confound certificate. The defect is confined to the raw representation.

Same certificate, same 1,000 within-cancer permutations, same 2,766 patients, run on
`cross_fitted_residuals(state, cancer + 84 pooled TSS)` — the exact adjustment CALIBRA applies before
it measures any channel. Chance rate 1/85 = 0.0118.

| artifact | state | joint LDA raw | **joint LDA adjusted** | joint null p95 (adj) | per-axis max raw | per-axis max adj | breaching axes raw → adj | verdict |
|---|---|---:|---:|---:|---:|---:|---|---|
| d2_h | **wsi_biology** | 0.3633 | **0.0118** | 0.0528 | 0.0532 | 0.0123 | 17 → **0** | **CERTIFIED** |
| d2_h | full_biology | 0.2630 | **0.0101** | 0.0668 | 0.0548 | 0.0107 | 60 → **0** | CERTIFIED |
| d2_h | rna_biology | 0.2563 | **0.0074** | 0.0654 | 0.0506 | 0.0139 | 58 → **0** | CERTIFIED |
| d2_i | **wsi_biology** | 0.2348 | **0.0052** | 0.0418 | 0.0511 | 0.0104 | 43 → **0** | **CERTIFIED** |
| d2_i | full_biology | 0.2689 | **0.0085** | 0.0758 | 0.0551 | 0.0102 | 61 → **0** | CERTIFIED |
| d2_i | rna_biology | 0.2744 | **0.0079** | 0.0732 | 0.0495 | 0.0106 | 48 → **0** | CERTIFIED |

Every adjusted joint accuracy is **at or below the chance rate of 0.0118**, and every one is well below
its own permutation null p95. Every per-axis maximum falls to ~0.010–0.014, i.e. to chance. Zero
breaching axes in every state.

The joint accuracy drops by a factor of **21–45×** (0.3633 → 0.0118 for the image-only channel). This
is not partial attenuation; the site signal is gone.

**This resolves the question the raw arm left open, in the favourable direction.** The defect recorded
at 01:52 UTC is a property of the **raw representation**, not of the instrument or of any adjusted
number. Concretely:

* Every CALIBRA channel number in this project is measured on adjusted states and is therefore **not**
  reading tissue source site.
* P4's condition-4 exposure rule still blocks the raw axes — an axis shown to a user is a raw axis.
  Certification must be issued against the adjusted state, and the certificate schema must record
  *which*.
* G3.5 in `e0_basis_transfer.py`, currently `unavailable_no_site_labels`, can be closed with two rows
  rather than one: raw FAIL, adjusted PASS.

#### Part 2 — P6, estimator robustness. Predeclared bar FAILED; the failure is informative.

Predeclared (`P1_PREDECLARATION.md`): "varying `n_splits` over {2, 5, 10, 20} and Ridge `alpha` over
{0.01, 1.0, 100.0} moves the anchor `|r_induced|` by less than 25% relative. Falsifier: >2× movement
means we are measuring our residualiser, not residualisation."

Anchor design `cancer+tss_pool10`, induced correlation:

| n | n_splits | alpha=0.01 | alpha=1.0 | alpha=100 |
|---:|---:|---:|---:|---:|
| 2,530 | 2 | 0.0828 | 0.0807 | 0.0380 |
| 2,530 | 5 | 0.0825 | 0.0817 | 0.0531 |
| 2,530 | 10 | 0.0805 | 0.0802 | 0.0573 |
| 2,530 | 20 | 0.0815 | 0.0809 | 0.0592 |
| 6,427 | 2 | 0.0922 | 0.0910 | 0.0632 |
| 6,427 | 5 | 0.0901 | 0.0904 | 0.0725 |
| 6,427 | 10 | 0.0905 | 0.0907 | 0.0760 |
| 6,427 | 20 | 0.0893 | 0.0894 | 0.0772 |

Decomposed:

* **Fold count: irrelevant.** At fixed alpha, the spread over `n_splits` ∈ {2, 5, 10, 20} is
  **≤ 2.4%** relative at both n. The cross-fitting scheme is not producing the effect.
* **Shrinkage over the sane range: irrelevant.** alpha 0.01 vs 1.0 differ by **≤ 3%**.
* **alpha = 100: 30–53% reduction.** Full-grid relative spread is **55.7%** at n = 2,530 and 32.3% at
  n = 6,427, against a predeclared bar of 25%. **The prediction as written is falsified.**

The falsifier's own threshold — ">2× movement" — is *not* breached (max/min = 2.18 at n = 2,530,
1.46 at n = 6,427; 2.18 is marginally over 2, so this is reported as a breach of the letter of the
falsifier too and not rounded down). But the mechanism is not an artefact and should not be presented
as one: alpha = 100 on a one-hot design at n = 2,530 heavily under-fits the nuisance model, so the
design explains less of *both* modalities (R_s and R_a fall), and Equation (1) then requires the
induced correlation to fall. The estimator knob is moving the very quantity the mechanism says drives
the effect.

**The honest statement, which is now what goes in the paper:** the induced correlation is invariant to
how the residualisation is cross-fitted and to shrinkage in the range anyone would actually use, but it
scales with *how much the nuisance model actually removes*. Under-adjusting reduces the induced
correlation and leaves the confound in — that is a genuine trade-off, not a robustness result, and
reporting it as robustness would have been wrong.

For the structureless `gaussian_k99` arm the relative spread looks enormous (1.70) but the absolute
values are 0.0010–0.0079, i.e. sampling noise on a quantity that is essentially zero. Relative spread
is not a meaningful statistic there and is reported only so its absence of meaning is explicit.

### In plain terms

**Good news.** The check that failed earlier — our axes being able to tell which hospital a slide came
from — is completely fixed by the correction we already apply before every measurement. The hospital
signal drops from "identifiable a third of the time out of 85" to "no better than guessing". So the
problem is in the raw representation, and none of our measured numbers are contaminated by it. What we
must not do is show a user a raw axis.

**Less good news.** We predicted in advance that fiddling with the knobs of our statistical correction
would change the induced-correlation number by less than a quarter. It changes by more than half if you
turn one knob — the shrinkage — all the way up. We are reporting that as a failed prediction, because
it was written down in advance and it failed. But it fails for a reason that supports the main claim
rather than undermining it: turning shrinkage up means the correction removes less, and the whole
theory says the effect is proportional to how much gets removed.

### Meaning for the claim

* **T1.3 becomes a two-row gate, and both rows are needed.** Raw: FAIL, 0.235–0.363 joint accuracy.
  Adjusted: PASS, ≤ 0.0118. A certificate that reports only one of these is uninformative in a
  different way in each direction.
* **Every channel number on this project survives.** The `wsi_biology` adjusted top-CCA of 0.605 (d2_h)
  is not a site code.
* **P4's certificate schema gains a required field:** `certified_on = {raw | adjusted}`, plus the
  joint-test row. Certifying on the adjusted state and then exposing the raw axis would be exactly the
  laundering P4 forbids.
* **Track 2 write-up must carry the alpha caveat**, phrased as a trade-off, not as robustness. The
  quoted magnitude 0.07–0.09 is for `alpha` ∈ [0.01, 1] and any fold count; at `alpha` = 100 it is
  0.038–0.077 and the adjustment is under-fitting.

### Files / commits

`v2/calibra/confound_certificate.py`, `v2/calibra/induced_correlation_sweep.py` (both inherited).
Results: `p1_evidence/track1/certificate_adjusted/{confound_certificate.json,task_rows.csv}`,
`p1_evidence/track2/{knobs_rows.csv,knobs_law.json}` under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/`.

---

## CORRECTION APPENDED 2026-08-04 21:40 UTC — "the site signal is gone" is refuted; the entry's numbers stand, its absolute sentences do not

This entry is append-only, so nothing above has been altered. Read this block as binding wherever it
contradicts the text above.

**What is withdrawn.** Four sentences in this entry are absolute and are refuted by measurement:

| line | withdrawn text |
|---:|---|
| 1 (title) | "The adjustment **fully discharges** the site confound" |
| 27 | "this is not partial attenuation; the site signal is gone" |
| 34 | "Every CALIBRA channel number in this project is measured on adjusted states and is therefore **not** reading tissue source site" |
| 89 | "…is **completely fixed** by the correction we already apply before every measurement" |
| 91 | "none of our measured numbers are contaminated by it" |
| 106–107 | "**Every channel number on this project survives.** The `wsi_biology` adjusted top-CCA of 0.605 (d2_h) is not a site code" |

**What still stands.** Every number above is reproduced and unchanged: joint accuracy 0.3633 → 0.0118
against a chance rate of 0.0118, the 21–45× drop, zero breaching axes in every state, per-axis maxima
falling to ~0.010–0.014. The entry's *finding about certification rules* — that a per-axis-only
certificate would have passed this leak, and that the joint test must be a required field of the
schema — is unaffected and is if anything reinforced below.

**What the measurement actually shows.** The adjustment removes the confound from the **first
moment**. LDA is a mean-based scorer, so a certificate built on it can certify exactly that and no
more. Probing the same adjusted block with classifiers that are not functions of class means recovers
both confounds: **3.15× chance for tissue source site and 3.45× for cancer type**, netted against a
null that regenerates the adjustment inside every permutation so each draw carries an adjustment
artefact of its own size, *p* at the 1/201 resolution floor, with k-NN, random forest and RBF-SVM all
agreeing. A synthetic block with exactly equal class means and unequal conditional variances shows the
blind spot directly: LDA reads 0.231 and the per-axis maximum 0.244 against chance 0.250, while k-NN
reads 0.554, the forest 0.625 and the SVM 0.658.

**The replacement claim is stronger than the one it replaces, and the entry's conclusion survives it.**
Three measurements bound the residual:

1. A **saturated cancer × site cell design** — one free parameter per occupied cell — spans every
   function of the confound labels and therefore **upper-bounds any conditional-mean adjustment
   whatsoever**, kernel, forest or boosted. It moves the `d2_h::wsi_biology` channel by 0.0001:
   **0.6052 → 0.6051**, retention 0.998.
2. The **confound labels alone** reach a channel of 0.1237 (additive design) and 0.0903 (saturated),
   both *below* the real channel's own null median of 0.1483 — **11.2%** and **6.0%** of its excess
   over that null.
3. The residual is therefore not a conditional mean of any shape.

So the 0.605 top-CCA is not a site code — but that is now a *measured and bounded* statement rather
than an inference from the certificate, and the correct sentence is "the confound is removed from the
first moment and cannot account for more than ~10% of the channel", not "the site signal is gone".

**A finding this entry could not have had: the certificate's raw row anti-predicts what survives.**
D1-B `programme_free` (`d1_f`) has the **lowest** raw joint site LDA of the twelve artifacts measured —
0.1071–0.1449, against `d1_p`'s 0.1778–0.3764 and `d2_i_seed43`'s 0.4735 — and the **highest**
post-adjustment nonlinear reading of all twelve, 6.21–7.38× chance for site and 4.98–6.69× for cancer.
A low raw joint LDA does not predict a low post-adjustment reading; on this cohort it anti-predicts it.
The raw joint row must therefore not be used as a proxy for "how confounded is this representation".

**Sources.** `NOTEBOOK_ENTRIES/tcga_nonlinear_confound_probe_result_20260804T2100Z.md` §7b (the
3.15/3.45 figures, stated there as excess over the regenerated-null median in × chance) and §8 (the
twelve-artifact breadth table); `NOTEBOOK_ENTRIES/nonlinear_adjustment_channel_result_20260804T2130Z.md`
§3–§5 (the saturated-cell arm, the labels-only ceiling). **Disagreement between sources, reported
rather than resolved:** the 21:30 entry expresses the same correction as a *ratio*
(`corrected_multiple` = observed ÷ regenerated-null median), giving **2.73× site / 3.66× cancer** on
this block and **2.4–3.3× / 3.1–3.7×** across six adjuster arms, where the 21:00 entry uses the
*difference* (observed × chance − null median × chance) giving 3.15× / 3.45×. Both are defensible and
they are not interchangeable; any quotation must name which.

**Not edited here.** `NOTEBOOK.md:1278` carries this entry's index row with the same absolute verb
("**adjustment fully discharges site**"). `NOTEBOOK.md` is out of scope for this correction pass and
is flagged rather than changed.

Corrected documents: `paper/P1_CALIBRA_DRAFT.md` (abstract, §1.4, §4.2, §6),
`paper/P1_FIGURES.md` (F2), `v2/research/rebase/nature/PHASE1_RESULT.md` (F1, validity checks),
`v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md` §T1.3.
