# MORPHEUS Rank-Collapse Study — Results Explained

**Status:** living results document. Geometry (T1), confound (T2/T3), and the seed-42 fix
ablation (T4) are complete; a multi-seed sweep with the molecular-prompting readout is running
on the A100 and this doc + the paper will be updated when it lands.

**Repos/branches (github.com/mj-jpgt/morpheus):** code fixes on `fix/biology-collapse`; paper +
figures + this doc on `paper/rank-collapse-diagnostic`.

---

## 1. One-paragraph summary

MORPHEUS V2 aligns whole-slide images (WSI) to a molecular-programme (MSigDB Hallmark) embedding
with a dual-head design: an **identity** head (paired image–RNA contrastive, anchored on a frozen
MLP-CLIP teacher) and a **biology** head (regressed/aligned to the Hallmark target). We find the
biology head uses **far less of its representational space** than its architecturally identical
sibling — a **~4–5× effective-rank gap** — and that this is a legible *fingerprint* of a benchmark
that is itself dominated by cross-cancer cohort structure. We prove a per-dimension variance floor
cannot repair this and show an off-diagonal **covariance-decorrelation** term can, confirmed on real
held-out data: it lifts biology-head effective rank **+49%** over the collapse baseline.

---

## 2. Setup

- **Model:** hierarchical patch-pooling Query-Former over uncapped H-Optimus WSI patches + BulkFormer
  RNA. Heads: `z_identity` (256-D, L2-normed, anchored bounded-residual on frozen MLP-CLIP),
  `z_biology` (256-D, L2-normed, supervised by Hallmark programme targets via Gaussian-NLL +
  programme-neighbour KL + supervised contrastive).
- **Protocol:** leakage-controlled cancer-held-out (11 development / 21 held-out cancers); all fits
  on the training fold only; held-out test split n = 2530.
- **Effective rank:** Roy–Vetterli — `exp(entropy(σ / Σσ))` on the singular values of the centred
  batch (a smooth count of "effectively used" dimensions; d = 256 for isotropic, → 1 for a single
  axis).

---

## 3. Results

### T1 — Dual-head rank geometry (the fingerprint)  ·  Fig. 1

On the held-out test split (seed 42), Roy–Vetterli effective rank of each head:

| head | effective rank / 256 |
|---|---:|
| wsi_identity | 191.1 |
| rna_identity | 141.6 |
| full_identity | 176.3 |
| **wsi_biology** | **38.5** |
| **rna_biology** | **32.6** |
| **full_biology** | **47.3** |
| full_patient | 60.5 |
| *(Hallmark target itself)* | *~92* |

The biology heads (~33–47) sit ~4–5× below the identity heads (~142–191) — and below the target's
own rank (~92), so the biology head is *contracting* under the regression objective, not merely
inheriting the target rank. The two heads share encoder, batch, and optimizer; they differ only in
what they are asked to match. The collapsed biology space also has a wider modality gap (0.475 vs.
0.296), consistent with two modalities pinned to a shared low-rank ridge.

> **Correction note.** An earlier internal analysis reported biology ~5–6 / identity ~84 — those used
> covariance *eigenvalues* (σ²), which over-concentrate the spectrum. The paper cites Roy–Vetterli
> (singular values); the numbers above are the correct ones. The qualitative story (biology ≪ identity;
> variance floor can't fix; covariance term can) is unchanged; the magnitudes are more modest.

### T2 — The benchmark is confounded  ·  Fig. 2

WSI→Hallmark molecular-prompting Pearson, pooled (cross-cancer) vs. within-cancer (macro-cancer),
mean over 180 Hallmark targets, seeds 42/43/44:

| method | pooled r | within-cancer r | lost to cohort |
|---|---:|---:|---:|
| MLP-CLIP | 0.348 | 0.188 | −46% |
| SigLIP | 0.349 | 0.193 | −45% |
| MORPHEUS-v2 (anchored) | 0.327 | 0.166 | −49% |
| MORPHEUS-v2 (no-anchor) | 0.338 | 0.185 | −45% |

**~46–49% of the pooled Pearson is cross-cancer cohort structure**, for every method including the
baseline. A size-matched random-gene-set null already reaches ~0.30–0.32 pooled / ~0.154 within-cancer.

### T3 — The genuine signal is small and method-invariant  ·  Fig. 3

Control-adjusted within-cancer specificity (real minus random-gene null, within-cancer):

| method | Δ specificity |
|---|---:|
| MLP-CLIP | +0.068 |
| SigLIP | +0.067 |
| MORPHEUS-v2 (anchored) | +0.068 |
| MORPHEUS-v2 (no-anchor) | +0.069 |
| MLP-CLIP hard-neg | +0.066 |
| MORPHEUS-v1 | +0.070 |

Genuine, morphology-resolved biology is **~+0.07 Pearson and identical across all methods** — the
richer models buy no real biological gain over the baseline. (SigLIP's +0.005 within-cancer edge is
inside the seed-noise band.) **Rank collapse and respectable scores coexist because the task rewards
cohort structure a low-rank embedding can still supply** — this is the coupling (C2).

### T4 — The fix recovers biology-head rank on real data  ·  Fig. 4

Real held-out cohort, `programme_only`, 25 epochs, A100. Both arms init at biology-head effective
rank 146.0:

| arm | final biology effective rank | vs. baseline |
|---|---:|---:|
| **baseline** (variance floor only, decorr = 0) | **68.0** | — |
| **F-R2** (covariance decorrelation, decorr = 0.04) | **101.1** | **+33.1 (+49%)** |

The baseline **collapses**; F-R2 **holds**. This confirms the primary prediction: the off-diagonal
covariance term recovers rank on real data where the per-dimension variance floor cannot.
The variance-floor-only baseline *is* the negative control (predicted flat/collapsed — confirmed).

**Multi-seed rank recovery (done).**  ·  Fig. 5. Seeds {42,43,44} × {baseline, F-R2},
`programme_only`, 15 epochs each on the A100. All runs init at rank 146.0:

| seed | baseline (decorr=0) | F-R2 (decorr=0.04) | Δ |
|---|---:|---:|---:|
| 42 | 49.0 | 102.8 | +53.7 |
| 43 | 51.0 | 101.5 | +50.6 |
| 44 | 49.8 | 105.5 | +55.7 |
| **mean ± s.d.** | **49.9 ± 0.8** | **103.3 ± 1.7** | **+53.3 ± 2.1** |

The recovery is **large (~2.1×) and remarkably tight across seeds** (baseline s.d. 0.8, F-R2 s.d. 1.7),
so the T4 effect is not a seed artifact. This closes the review's multi-seed BLOCKER for the
rank-recovery claim (C3-i).

**Specificity readout (ii) — done, and it confirms the prediction.** A ridge probe of
`wsi_biology → regenerated held-out Hallmark targets` (full test coverage, n=2530) gives within-cancer
Pearson **0.1366 (baseline, rank 49)** vs **0.1367 (F-R2, rank 103)**; pooled 0.217 vs 0.219. **Doubling
the biology-head rank does not move molecular-prompting specificity** — exactly the pre-registered
"rank rises, specificity unchanged." This is the pivotal coupling result (C2): effective rank is
**decoupled from the confounded benchmark score**, so the collapse is a pathology the metric is blind
to — which is *why* a rank fingerprint is diagnostic where the benchmark isn't. (Earlier this probe
returned nan because the loader's hallmark is a train-fold-only constant placeholder on the held-out
test split; fixed by scoring against regenerated `frozen_rna_targets`, commit `337b2b1`.)

Interpretation caveat (ties to the reviewer critique): "specificity unchanged" also means the fix does
**not** improve real molecular prediction on this benchmark — consistent with C2 (the benchmark can't
reward the fix). Whether recovering rank helps on a *de-confounded / site-stratified external* task is
the separate, stronger experiment still to run.

---

## 4. What the adversarial process caught (and we fixed)

The work was checked by an 8-agent adversarial workflow. It found real problems, all fixed:

- **F-R2 was a silent no-op on real batches.** Uncapped-patch batches hold only B≈1–3 patients —
  below a usable 256-D correlation estimate — so the decorrelation term returned 0 and the "ON arm"
  equalled the baseline. Fixed with a **biology feature bank** (pool current gradient-carrying rows
  with a detached ring buffer of recent features); a guard test now fails on the old no-op. *Without
  this fix the T4 ON arm would read 68, not 101.*
- **A metric bug** in the headline numbers (eigenvalue vs. singular-value effective rank) — corrected
  to Roy–Vetterli throughout (see T1 note).
- **Two mis-attributed citations** (Wang&Isola 2020 → Wang et al. CVPR 2022; "Buyer Beware" → Dawood
  et al. 2024) — corrected.
- **Overclaimed "provable."** Restricted "provable" to the impossibility result (variance floor cannot
  restore rank); the remedy is proven-necessary + synthetically validated + now real-data confirmed.

---

## 5. Paper-readiness

**Diagnostic-paper scope (C1 geometry + C2 confound + C3 fix):** supported, with the real-data T4
now in hand. Remaining before submission:

- **Multi-seed T4 rank recovery** — ✅ **done** (49.9 ± 0.8 → 103.3 ± 1.7, +53.3 ± 2.1 over 3 seeds; closes the multi-seed BLOCKER for C3-i).
- **T4 specificity readout (ii)** — pending one follow-up sweep (probe fixed in code); completes the T4 triplet.
- **Multi-seed *identity/biology geometry* (T1)** — only seed-42 exports exist; the multi-seed run above is `programme_only` (biology-only), so multi-seed *identity*-vs-biology geometry still needs full-profile exports for 43/44.
- **The causal downstream probe (the stronger novelty).** Freeze `z_biology` (collapsed),
  `z_biology` post-F-R2 (recovered), and `z_identity`, then probe real endpoints (survival C-index,
  molecular subtype, driver-mutation AUROC) under **site-stratified vs. naive splits**. If rank
  recovery helps on site-stratified splits, the contribution upgrades from "diagnostic" to
  "a label-free geometric criterion that yields confound-robust representations." This is the #1
  reviewer-identified experiment and the recommended next run (needs saved model representations).

**Honest caveats.** T4 is `programme_only` (the baseline plateaus ~68 vs. the full multi-loss export's
~47); the benchmark critique (C2) is the least novel part (conceded to Howard 2021 / Dawood 2024 /
DECAT); the novel object is the *coupled* fingerprint, which the causal probe would most strengthen.

---

## 6. Reproducibility

- Fix commits (`fix/biology-collapse`): stress suite + `effective_rank` (`6ac5e7e`), normalize
  (`ab10be0`), decorrelation (`13b6337`), feature-bank fix (`960e0ec`), ablation runner + probe.
- Ablation: `morpheus/v2/run_rank_ablation.py`; stress tests: `morpheus/v2/tests/test_stress_collapse.py`
  (10, all green on the A100 box); honest metrics: `morpheus/v2/honest_metrics.py`.
- Paper: `paper/main.md` + `paper/figures/` (generators included) + `paper/results/t4_{off,on}.json`.
- Full adversarial review: `paper/REVIEW_SYNTHESIS.md`.
