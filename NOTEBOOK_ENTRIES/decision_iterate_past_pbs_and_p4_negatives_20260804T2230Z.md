# Decision: negative results are checkpoints, not conclusions — reorienting P3 and P4

2026-08-04 · main session. Corrects how the last several entries were being read.

## The correction

`t11_t12_must_beat_baselines_20260803T0440Z.md` (PBS loses to PCA) and
`p4_certification_end_to_end_20260804T2000Z.md` (0/90 answerable) were reported and then treated as
near-final verdicts on their respective directions. That is the wrong read of what a fair-test loss
means on this project. **A negative result under a controlled comparison is the point at which the
construction gets redesigned, not the point at which the direction gets written up as dead.** The
goal was never "PBS in its current form" or "the current P4 pipeline" specifically — it was a working
method for morphology-legible causal structure, and a working promptable multimodal system. Neither
of those goals is falsified by one construction losing a fair fight.

## P3 — what the loss actually says, and the next construction

`v2/pbs.py::ReferenceDictionary.fit` already does something smarter than "one atom per perturbation":
it SVDs the K562/RPE1 CRISPRi response matrix into a low-rank quotient basis, then projects patient
RNA onto it. It lost to plain PCA of the patient expression matrix. Read structurally rather than as
"the idea is bad": **PCA of patient bulk RNA captures the dominant axes of variation across a tumour
*cohort*** (proliferation, immune infiltration, stromal content — exactly what shows up in an image).
**The interventional dictionary captures directions activated by perturbing single genes in cell
lines**, which need not align with cohort-level tumour heterogeneity at all. The two bases are
answering different questions; testing them as competing candidates for the same slot was the
mismatch, not the interventional framing itself.

**Next constructions to test, in order of promise, all reusing the existing capacity-matched
leak-safe harness (`baseline_paired_bootstrap.py`, `run_calibra`, the `test != development` split
discipline) so the comparison stays apples-to-apples with what's already measured:**

1. **Attribution, not competition.** Stop asking the interventional dictionary to *be* the legible
   basis. Instead, take the axes PCA (or the trained representation) already finds legible, and
   regress/project each onto the perturbation-atom space to ask *which single-gene perturbations
   reproduce this direction*. This reframes P3's contribution from "our basis is more legible"
   (falsified) to "we can name a legible axis in causal terms a PCA component cannot supply" — closer
   to novel-in-application than the original framing, and it is exactly what a promptable system (P4)
   needs: not just "axis 46 correlates with target X" but "axis 46 matches the transcriptional
   signature of perturbing gene G."
2. **A joint basis**, via CCA or PLS between the perturbation-response covariance and patient
   expression covariance, rather than an OR between two orthogonal candidates. Tests whether *some*
   causally-grounded subspace is legible, even if the raw SVD basis wasn't.
3. **A denoised/consensus interventional basis** — filter to perturbation atoms with reproducible
   effect across both K562 and RPE1 before the SVD, or shrink toward a background covariance model,
   to separate "the causal framing is wrong" from "this particular perturbation resource is noisy and
   cross-cell-line noisy."
4. **Domain-adapted projection** — CRISPRi delta (cell line, single-cell) and bulk tumour RNA
   (steady-state, mixed cell type) are different measurement regimes; align them before projecting
   rather than assuming a shared coordinate scale (current code already flags this as a declared
   choice in `fit_development_expression_transform`, not a proof it's the right one).

Option 5, explicitly flagged rather than pursued first: using the interventional dictionary as a
*training-time* regularizer on the biology head. Deferred because it risks reproducing the D1/D2
finding that molecular supervision degrades the molecular channel — worth trying, but only after 1–3
have a result, so a repeat of that failure mode doesn't get mistaken for a new one.

## P4 — the same correction, plus the modality gap the user actually asked for

The inductive-adjustment operator (agreed last turn, agent already tasked) is the mechanical fix
inside the *current* pipeline. It is necessary but it is not the ambition. The stated goal is a
promptable system across **WSI, RNA, CNV, SNV, and proteomics**, prompted in natural language — today
the pipeline is WSI+RNA only.

**What's already on disk, worth scoping rather than assuming:** `known_covariate_control.py` already
carries mutation-call targets (BRAF/KRAS/APC/PIK3CA/dMMR read as clinical-band checks in an earlier
entry), which means SNV is at least partially wired as a control, not a supervised modality — the
gap to a real SNV channel may be smaller than starting from zero. CNV and proteomics have no code
path found yet. TCGA carries RPPA proteomics (~200–250 antibodies, all cohorts with tissue) and
GISTIC/segment-level CNV without needing CPTAC or any new external access — both should be scoped
before assuming a new resource is required.

## What this does not change

Both prior entries stand as measured — the numbers are not being revised or softened, only what
happens next. `t11_t12_must_beat_baselines` and `p4_certification_end_to_end` are not edited.

Related: [[t11_t12_must_beat_baselines_20260803T0440Z]], [[p4_certification_end_to_end_20260804T2000Z]]
