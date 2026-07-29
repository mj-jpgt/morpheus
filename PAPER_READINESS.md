# Diagnostic Paper — Readiness Verdict

**Paper:** *Rank Collapse as a Fingerprint of Cohort-Confounded WSI→Molecular Alignment.*
**Branch:** `paper/rank-collapse-diagnostic`. **Verdict: submittable as a diagnostic/mechanism
paper on its three core claims.** One enhancement experiment would upgrade it to the stronger
"confound-robustness" claim; it is not required for the diagnostic contribution.

## The three claims are all backed by real data
- **C1 (mechanism).** Matched dual-head rank asymmetry on real held-out data (seed-42 export):
  biology heads ~33–47 vs. identity heads ~142–191 effective rank (Roy–Vetterli, of 256). Fig. 1.
- **C2 (coupling).** ~46–49% of pooled molecular-prompting Pearson is cross-cancer cohort structure;
  control-adjusted within-cancer specificity is a **method-invariant ~+0.07** across all methods incl.
  baseline. Figs. 2–3. And the decisive coupling result: doubling biology rank leaves specificity
  **unchanged** (0.1366→0.1367) → rank is **decoupled** from the confounded score.
- **C3 (fix).** Multi-seed (42/43/44): a covariance-decorrelation term recovers biology rank
  **49.9±0.8 → 103.3±1.7 (+53.3±2.1, ~2.1×)** where the per-dimension variance floor cannot. The
  variance-floor arm is the confirmed negative control. Figs. 4–5.

## Review must-fixes — all cleared
- [x] Multi-seed **T4** rank recovery (was the #1 BLOCKER).
- [x] Real-data T4 **triplet** complete: rank↑ + variance-floor-flat + specificity-unchanged.
- [x] Metric corrected to **Roy–Vetterli** (the σ vs σ² bug) throughout + fig1.
- [x] Citation mis-attributions fixed (Wang et al. 2022; Dawood 2024).
- [x] "provable" softened to the proven impossibility result.
- [x] F-R2 wiring fix (feature bank) — without it the ON arm reads ~50, not ~103.
- [x] Internal consistency: 5/5 figures present+referenced; numbers consistent; no dangling
      placeholders except the author line.

## Two known gaps (enhancements, not blockers for the diagnostic core)
1. **Multi-seed *identity* geometry (T1) for seeds 43/44.** Only the seed-42 full-profile export has
   the identity head; the multi-seed run is `programme_only` (biology only). *To close:* one
   full-profile export per seed (needs the MLP-CLIP anchor), then recompute identity rank. Cheap-ish,
   strengthens C1's error bars.
2. **The site-stratified downstream causal probe (the stronger claim).** Freeze `z_biology`
   (collapsed), `z_biology` post-F-R2 (recovered), and `z_identity`; probe survival/subtype/mutation
   under **site-stratified vs. naive splits**. If rank recovery helps site-stratified transfer, the
   paper upgrades from "diagnostic" to "a label-free criterion that yields confound-robust
   representations." *This is the highest-value next experiment* but a larger build (needs site labels
   + saved model reps); it is what the rebase fleet independently flagged (PATH 1/BioELK-Bench is its
   generalized form).

## Recommendation
Submit the diagnostic/mechanism version now to **NeurIPS D&B or MICCAI** (its restraint — no SOTA
claim — is appropriate for the genre). Hold the "confound-robust representations" framing until the
site-stratified probe lands. The honest one-line status: **the core science is done, real, and
reproducible; the paper is submittable; the causal-robustness upgrade is optional and scoped.**
