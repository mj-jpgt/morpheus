# Review Synthesis — The Effective-Rank Fingerprint (MORPHEUS)

**Synthesized:** 2026-07-28 · **Inputs:** 7 adversarial reviews (3 novelty/SOTA, 1 SOTA-experiment, 2 code, 1 citations) in `paper/reviews/`
**Target venues:** NeurIPS D&B / MICCAI

---

## Executive Verdict (200 words)

**CONDITIONAL — not submittable today, but a real diagnostic paper is within reach.** Two of three novelty reviewers land at *conditional* (`novelty_sota_1`, `novelty_sota_2`); the harshest (`novelty_sota_3`) says *no*. The disagreement is not about quality — the paper is honest, well-written, and claims no SOTA (`novelty_sota_2` calls the restraint "exemplary"). It is about **the coupling claim being asserted, not shown**. All three agree every mechanism is conceded prior art (C1≈NRC1, C2's goal≈DECAT, C3≈VICReg's own motivation); the only novel object is "biology-head effective rank as a confound fingerprint," and it currently rests on **a single seed of one unreleased model with the load-bearing experiment (T4) unrun**. `novelty_sota_3`'s W1 is the sharpest threat: erank tracks *target rank*, not confounding, and this is never severed.

Two independent code reviews rate the fix lane **HIGH severity**: as wired, F-R2 decorrelation is silently inert on real A100 batches, so T4 could show "no recovery" for an engineering reason. Citations has two confirmed mis-attributions that must be corrected. The flip-to-accept path is concrete: soften "provable," add multi-seed geometry, add a positive/negative control, fix the code, and run the confound-controlled downstream test.

---

## 1. Overall Verdict — CONDITIONAL (leaning accept-as-diagnostic *iff* conditions met)

| Reviewer | Verdict | Core reason |
|---|---|---|
| `novelty_sota_1` | conditional | Thesis *not scooped*, but strongest-reject wins as-is: mechanisms conceded, one new claim is single-seed + T4 unrun. |
| `novelty_sota_2` | conditional | No SOTA overclaim (a strength); real exposure is the word "provable" and single-seed C1. |
| `novelty_sota_3` | **no** | Coupling (the only novelty) is a non-sequitur — erank tracks target rank, never severed from confound (W1–W4). |
| `sota_experiment` | (constructive) | "So what?" objection is the blocker; a confound-controlled payoff test kills it. |
| `citations` | fix-required | 2 confirmed mis-attributions propagate across abstract/related-work/risks. |
| `code_leakage_scaling` | HIGH | F-R2 inert/degenerate on real uncapped batches; T4 could fail for an implementation reason. |
| `code_correctness` | HIGH | F1 (min_batch no-op), F2 (wrong view decorrelated), F3 (no test catches either). |

**Why conditional, not accept:** The paper's sole novelty is a *coupling* (collapse ↔ confound), and every novelty reviewer independently finds the coupling is narrative, not demonstrated. `novelty_sota_3`-W1 (erank explained entirely by low-rank target = NRC1, with no reference to cohort structure) and W3 (the paper's own method-invariance shows rank is *decoupled* from the only honest metric) are the load-bearing objections — neither is cured by the queued T4.

**Why not reject:** `novelty_sota_1` confirms the exact thesis is **not scooped**, and `novelty_sota_2` confirms the negative-result/diagnostic genre publishes without SOTA (Venet 2011, Geirhos 2020, Howard 2021). The instrument — a confound read from one forward pass, no confounder label, no held-out cohort — is a real methodological reduction. It becomes publishable once the evidence is completed and the causal claim is either proven or honestly downgraded.

---

## 2. Ranked MUST-FIX List (before submission)

1. **[BLOCKER] Sever "erank fingerprints confounding" from "erank tracks target rank."** (`novelty_sota_3`-W1/W4, `novelty_sota_2`-§2) Add at least one positive/negative control: a low-rank target that is *not* cohort-confounded but still collapses, OR a de-confounded/site-held-out model whose biology head stays high-rank. Without this the central novel claim is, by three reviewers, applied NRC1.
2. **[BLOCKER] Multi-seed C1 geometry.** (`novelty_sota_1`, `novelty_sota_2`, `novelty_sota_3`-W13) The 84.3 vs 6.0 gap is single-seed on one proprietary model. Report seeds {42,43,44} with CIs. Cheapest, highest-value fix — a "fingerprint" is a measurement claim and cannot rest on n=1.
3. **[BLOCKER] Soften "provable."** (`novelty_sota_2`-§1, `novelty_sota_3`-W6) The Proposition proves only the *negative* control (variance floor cannot restore rank). Restrict "provable" to that impossibility result; describe decorrelation as "prescribed, synthetically validated, real-data validation queued (T4)."
4. **[BLOCKER] Fix the two citation mis-attributions** (see §3) — they touch abstract, related work, concession list, and risks table.
5. **[HIGH] Fix F-R2 wiring before running T4** (see §4) — otherwise T4's null result is uninterpretable.
6. **[HIGH] Run the real-data T4 triplet** (`novelty_sota_2`-§3): biology erank↑ (F-R2 on) AND variance-floor arm flat (~5–6) AND Δ_spec unchanged at ~+0.07. Required only for the *fix* (C3) claim; the diagnostic (C1+C2) can stand without it if 1–3 land.
7. **[MEDIUM] Add a linear cancer-type / tumour-purity probe of `z_biology`** (`novelty_sota_3`-W9) — substantiates the "absorbed cohort structure" mechanism directly.
8. **[MEDIUM] Scope "method-invariant."** (`novelty_sota_3`-W7) The compared set is a near-clonal CLIP family sharing one teacher. Either add a non-CLIP / dedicated regressor (HE2RNA-style) baseline or restate as "invariant across the tested embedding family."
9. **[MEDIUM] Coherence-matched null.** (`novelty_sota_3`-W8) The random-gene null matches set *size* but not co-expression coherence, so Δ_spec may be partly coherence, not morphology-resolved biology.
10. **[LOW] Reconcile 50-D vs 180-target and define the intrinsic-rank measure** (`novelty_sota_2`-minor, `novelty_sota_3`-W10/W11); add Roy–Vetterli 2007 and Tirosh 2016 to the bib.

---

## 3. Confirmed Citation Problems to Correct

Both confirmed by `citations.md` via WebSearch/WebFetch; no hallucinated refs, but two real entries are attached to claims made by *other* real papers.

- **M1 — `wang2020understanding` (Wang & Isola, ICML 2020) is the wrong paper.** The "minimal-sufficient representation discards task-relevant info" claim (§2.4, §5.4, §6.4, abstract) belongs to **Haoqing Wang et al., "Rethinking Minimal Sufficient Representation in Contrastive Learning," CVPR 2022 (arXiv:2203.07004)**. The in-text label "Wang & Isola, CVPR 2022" is a garbled hybrid (ICML-2020 authors + CVPR-2022 venue). Replace the entry with `wang2022rethinking`; relabel to "Wang et al., CVPR 2022" everywhere.
- **M2 — `tizhoosh2026rethinking` (NBME 2026) is the wrong paper for "Buyer Beware."** The confounder/TMB-inflation claim (§2.2, abstract, R2) belongs to **Dawood et al., "Buyer Beware…," bioRxiv 2024 (doi:10.1101/2024.06.23.600257)**. Add `dawood2024buyer`, retarget all "Buyer Beware" cites to it, and drop the misleading "NBME 2026" label and the `"buyer beware"` note field that drove the conflation.
- **Minor:** add prose-only Roy & Vetterli (EUSIPCO 2007) and Tirosh et al. (Science 2016 / Seurat `AddModuleScore`) to `references.bib`. The other 15 entries verified correct.

---

## 4. Code Findings to Act On (severity-ranked)

Two independent reviews converge on the **same HIGH-severity failure**: the F-R2 fix is inert on real data. Leakage controls are otherwise clean (both reviews confirm `feature_decorrelation` and `honest_metrics` are batch-local; within-group Pearson is even affine-robust to eval-time residualisation).

- **[HIGH] F1 — `min_batch=8` silently no-ops F-R2 on real uncapped batches.** (`code_correctness`-F1, `code_leakage_scaling`-Q2) `DynamicTokenBatchSampler` packs uncapped H-Optimus bags; at token_budget 16384–32768, B≈1–3 per batch, below the guard. The decorrelation term returns 0.0 on most/all batches, so the "fix arm" ≈ the `--decorrelation-weight 0.0` baseline. OOM handling *lowers* the budget, further disabling it. **The queued T4 could show "no rank recovery" for an implementation reason, not a scientific one.** Fix via feature bank / cross-microbatch covariance accumulation (mirror the existing `ProgrammeMemoryBank`).
- **[HIGH] F2 — decorrelation is applied to `full_biology`, but collapse pressure (and the likely-reported fingerprint) sits on `wsi_biology`.** (`code_correctness`-F2) Anti-collapse term and collapse-inducing term are on different views. Apply F-R2 to whichever biology view the paper reports rank on (most likely `out_wsi["z_biology"]`), or state explicitly the fingerprint is `full_biology`.
- **[HIGH-adjacent] Instrumentation gap.** (both) Skips are silent; logged `decorrelation` mean is dominated by 0.0s. Emit `decorrelation_active_fraction` and mean-B per epoch. If active-fraction ≈ 0 on the real run, T4 is invalid.
- **[MEDIUM] F3 — no test would fail if F-R2 were disabled.** (`code_correctness`-F3) Add an A/B test: same fixture, weight 0.0 vs 0.04, assert trained erank materially higher with it on (guards F1+F2).
- **[MEDIUM] Degenerate estimator when active.** (`code_leakage_scaling`-Consequence 4) 256-D correlation from ~7 effective samples is noise-dominated; also fires only on the small-WSI stratum (a cohort-correlated confound). Compute on a lower-dim projection or require N ≫ dim.
- **[LOW] `--expected-heldout-cancers` defaults to 22, protocol says 21.** Pass `21` explicitly before the run or preflight asserts the wrong contract.
- **[LOW] F4/F5** — `variance_floor(target_std=1.0)` is near-unsatisfiable after L2-norm (mildly fights F-R1); biased/unbiased N vs N-1 mismatch (~+7% at N=8). Cosmetic; note or retune.

**Gate before trusting T4:** land F1 (feature-bank accumulation) + instrumentation at minimum; align the decorrelated view with the measured view (F2).

---

## 5. Single Best Experiment to Run Next

**Proposal 1 (`sota_experiment.md`) — the confound-controlled downstream causal test.** Freeze three embeddings from the same model — `z_biology` (collapsed, erank~5), `z_biology_decorrelated` (post-T4, restored), `z_identity` (~84) — plus frozen UNI2/Virchow2/H-Optimus baselines. Probe three real clinical endpoints (TCGA-CDR survival Cox C-index; PAM50/TCGA molecular subtype; MC3 driver-mutation AUROC), each under **site-stratified (Howard 2021) *and* naive-random splits, reported side by side.**

**Falsifiable primary hypothesis:** rank restoration does *not* help on naive splits (confound-saturated, consistent with T4's +0.07 prediction) but *does* help on site-stratified splits — i.e. **fixing collapse buys confound-robustness specifically.**

**Why this over Proposals 2/3:** It directly kills the objection every reviewer raises — `novelty_sota_3`-W3 ("rank is decoupled from any honest metric; the fix by your own prediction changes nothing") and `sota_experiment`'s "so what?". It reuses on-disk data (no new dataset, no retraining for the frozen-probe version), rides the queued A100/T4 window, and is **publishable whatever the outcome**: a gain reframes the paper as "a label-free geometric criterion that produces confound-robust representations" (SOTA-relevant training principle); no gain sharpens C2/C3 into "collapse is genuinely orthogonal to utility." Proposal 2 (high-rank targets → HEST leaderboard) is the SOTA-upside follow-on; Proposal 3 (cross-model erank as a selection metric) is the integrative reframing that consumes Prop 1's site-stratified gap as its target variable.
