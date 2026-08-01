# HANDOFF — PHASE D: build. **This supersedes the E-series as the active task list.**

**Read in this order:** this file → `HANDOFF_GATES.md` (mandatory; G2 *liveness* matters most here because
Phase D trains) → `v2/research/rebase/ENGINE_CLD.md` (why). Branch `research/rebase-vision`, currently
`e33a853`. `HANDOFF_EXPERIMENTS_NOW.md` (E0–E5) is **closed** — see §1.

---

## 0. CONTEXT UPDATE — what changed since you last had context

### The instrument was broken and is now fixed
CALIBRA's spike readout scored recovery with `top_canonical_correlation`, a **maximum over 16
components**, while the spike lives on **one** known direction. Ambient sits at ~0.97, so every
detection floor came back `NaN` on real data while all 11 synthetic self-tests passed. Three nested
defects (max-readout; partial replacement because `a` was standardised and `y·v` was not; absolute value
taken *before* pairing, which destroys the paired comparison since induced correlation has random sign).

**Consequences you must know:**
- **The confound adjustment does NOT destroy signal** — attenuation 0.94–1.23, i.e. ≈1. This retires the
  objection that killed three earlier theses.
- **Residualising two orthogonal signals through a shared confound design INDUCES correlation between
  them** (0.067–0.140 for the 99-column cancer+TSS design at n=2,530). Nobody reports this. It is why the
  floor is read with a paired test.
- Two floors now exist and are **not interchangeable**: `transmission_floor` (paired, near-noiseless,
  **never quote as a detection limit**) and `detection_floor` (unpaired, conservative, quotable ≈0.2 WSI).

### F2 IS WITHDRAWN
E3 ran. `wsi_identity` changes by **2.6e-04** between `full` and `identity_only`, against **1.4e-01** for
the biology head. The identity head is the **frozen MLP-CLIP teacher passed through**. So "the head
trained for biology is worse at biology" restated *"the frozen teacher beats our biology head"* — a
distillation observation, not a claim about objectives. **Do not cite F2.** Any draft text asserting
"molecular supervision degrades the molecular channel" must be removed until Task D1 provides the arm
that was never trained.

### What replaced it is stronger
Effective rank **−17%** (38.48 → 32.06) with the molecular channel **unchanged** (0.4768 → 0.4748).
Combined with the earlier +107% rank at flat specificity, **rank ≠ information is now demonstrated in
both directions** across independent experiments.

### E0 / E0b are done; the E-series is closed
- E0 as originally written **could not return a negative**: its "matched-spectrum null" ignored the
  spectrum and was a Haar-uniform random subspace (~k/n_genes) while any two real expression matrices
  score 0.4–0.8 from generic co-expression. A zero-shared-biology construction passed every gate. Fixed
  with a **biological control arm** (`energy_test_p_value`: responsive p<0.01 vs non-responsive p>0.5).
- **E0 result:** `supported`, K562, **~10% of the achievable ceiling**, `pc1_share` 0.09–0.21 (not library
  size). The control arm absorbs **55%** of the raw overlap at k=100 — half of what a random-floor rule
  would have credited to perturbation biology is generic expression structure.
- **Replication:** cross-lineage `supported`, 4/4 contexts. At k=25, K562 +0.0387 vs RPE1 +0.0394 — **within
  2%** when n-matched. The n-matching was load-bearing: uncapped, K562 reads +0.0671 vs +0.0496, a **35%
  inflation from sample size alone** that would have been reported as lineage specificity.
- **E0b:** 8,403 perturbations have **effective rank 132.1** (RPE1 113.9), stable rank 17.4, coherence
  0.85. *"11,000 independent causal directions" is dead.* **`n_equivalence_classes` returned n and is
  mis-specified — do not use it.**

### Decision: stop validating, start building
Further E0-phase work has diminishing returns because **E0 does not measure the operation we are
building.** It is a feasibility screen ("can this dictionary express tumour biology at all?" — yes,
partially, replicated). The causal claim is carried by **per-axis certification downstream**, not by E0.
GTEx is **dropped** — it interrogates E0's statistic, not the system.
**The proliferation confound is answered better by Task D2 than by any further E0 control** (see §3).

### Claims are now enforced in code
`v2/calibra/claim_guards.py` + `tests/test_claim_guards.py` (15 tests). Run `validate_claim()` before
any write-up. **E0 is currently an INADMISSIBLE transfer claim**, blocked on `proliferation_deflation`
and `single_platform`, pinned by `test_current_e0_result_is_not_yet_an_admissible_transfer_claim`. If
you discharge one, that test fails and you must update it **deliberately**.

---

## 1. Your tasks, in order

### D1 — The missing arm. **Does molecular supervision help or hurt?**

This question has never been answered because **the arm does not exist on disk.** `identity_only`
declares no biology states (`v2/runner.py:351`) — by design, since the biology head has no gradient path
under identity-only losses. So "biology head without programme supervision" was never trained.

**Build a new objective profile `programme_free`** (`v2/training.py:47`, `v2/runner.py:344`):
- biology head trained with **RNA-paired InfoNCE** — a rank-preserving contrastive signal against the
  paired RNA view;
- **no** Hallmark regression, **no** programme neighbour-KL, **no** supcon (these are the diagnosed
  collapse mechanism: they pin a 256-D head to a ~50-D manifold);
- declares `wsi_biology`, `rna_biology`, `full_biology`;
- **no MLP-CLIP anchor** (`runner.py:390` already skips the anchor for `programme_only`; do the same
  here, or D1 inherits the F2 anchoring artifact and answers nothing).

**The comparison.** `programme_only` (exists) vs `programme_free` (new). **Matched on epochs, LR, seed,
token budget, batch schedule — verify by diffing `manifest_json`, do not assume.** ≥3 seeds. If they are
not matched, the result is *suggestive, not causal*, and you say so.

**Measure:** CALIBRA held-out channel (`run_calibra.py` with `--targets-npz frozen_rna_targets.npz` —
`data.hallmark` is train-fold-only and **constant on the test split**, which is where the original `nan`
came from), effective rank, and `honest_metrics` within-cancer specificity. **Paired bootstrap on the
difference** — the F2 write-up had no CI on its headline gap and that is partly how it survived.

**Prediction to pre-register before running:** if `programme_free` ≥ `programme_only` on the molecular
channel, molecular supervision is not helping and may be harming → F2 is restored *as an objective
claim* and PBS is motivated. If `programme_only` wins, the collapse story needs rewriting and you
escalate.

### D2 — PBS head-to-head. **Interventional coordinates vs curated pathways.**

**This is the real test, and it subsumes the proliferation confound.** Hallmark already contains
proliferation programmes. If our dictionary's only real content is proliferation, it **cannot** beat
Hallmark — there is nothing left to win on. If it does beat Hallmark, it carries something curated
pathways lack, which is exactly the claim we want.

**Dictionary resolution — do not use 8,403 atoms.** E0b measured effective rank **132**, and the
equivalence-class clustering is broken. Use the **top-k right singular vectors of P** (k ≈ 100–132) as the
supervision basis. Per-gene attribution comes later and requires fixing the clustering first.

**Targets.** Project each patient's RNA onto that basis: `a_i = argmin ‖y_i − D a‖² + λΩ(a)`.
Reuse `v2/pbs.py` (`ReferenceDictionary`, `LegibilityOperator`, `weighted_code_loss`) — already committed,
do not reimplement.

**Arms**, identical architecture / budget / seeds (≥3):
- **H** — biology head regressed onto ~50 Hallmark scores (= current baseline).
- **I** — biology head regressed onto interventional coordinates `a_i`.

**Report:** held-out CALIBRA channel, effective rank, within-cancer specificity, and the standard
benchmark. **Plus the field that answers proliferation for free:** per-axis proliferation/essentiality
loading. You will have per-axis gene loadings anyway — tag every axis. *If every legible axis comes back
proliferation-loaded, that is the deflation, and you will see it without a separate experiment.*

**Escalate immediately if I ≈ H.** That means the interventional dictionary's content is already inside
curated pathways, and the central premise of the rebase is in trouble.

### D3 — Purity into the adjustment set. **Do this while D1/D2 train.**

Open since Phase 1, never closed, and it is a **complete alternative explanation** for the
morphology↔molecular channel: bulk RNA is a 30–90% tumour mixture, dictionary atoms are pure
populations, so coefficients absorb purity — and purity is one of the most visually obvious features on
a slide.

**No TCGA purity table is on disk.** (CPTAC has `tumor_purity__washu.parquet` per cohort — relevant later
for the external cohort, not for TCGA.) Two routes:
1. **Preferred:** obtain published TCGA consensus purity (ABSOLUTE / PanCanAtlas). Open access; ask.
2. **Fallback:** compute an ESTIMATE-style stromal+immune score from the expression we hold. **Flag it as
   expression-derived** and state the partial circularity — you are residualising expression on a
   quantity computed from that same expression.

Add to `confound_design` in `run_calibra.py`, re-run the Phase-1 channel measurement, and report the
effect **before and after**. If the channel dies when purity enters, that is a finding, not a failure —
report it.

---

## 2. Gates that bite hardest here

Phase D **trains**, so the G2 liveness family applies in full and has bitten us before:
- **G2.6 overfit-one-batch** — if the head cannot memorise one fixed batch, nothing downstream is
  interpretable. Cheapest decisive test there is; run it first on every new profile.
- **G2.1** ‖Δθ‖/‖θ‖ > 1e-2 — `feature_decorrelation` once contributed ~1e-6 with a 4e-5 param delta.
- **G2.2** log **every loss term's magnitude separately**; any term < 1e-4 × the largest is effectively off.
- **G2.7** guards must fire at the *real* batch size — an `min_batch=8` guard silently no-oped on every
  real step because uncapped H-Optimus batches hold B≈1–3.
- **G4.1 positive control** (`rna_*` → RNA targets) **must pass in the same run**, or you have measured
  nothing. **G1.1**: no constant columns *on the evaluated split*.
- **G3.1** effective rank via `calibra.spectral.effective_rank` — **singular values**, never covariance
  eigenvalues (that 6× error reached a paper draft).

---

## 3. Do not

- Do not cite F2, or the marginal `bootstrap_ci95` as a 95% CI (it is biased low by 0.04–0.09 and flagged
  `bootstrap_ci95_is_biased_low`; the decision uses the **paired** difference).
- Do not compare gap magnitudes across E0 runs with different n or q, or `normalised_alignment` across runs.
- Do not use `n_equivalence_classes`, or claim "11,000 independent causal directions".
- Do not claim E0 shows the alignment is *biological* — both lineages share one platform and one
  effect-size-monotone statistic; that confound replicates across lineages **because it is biology of
  essential genes, not of lineage**.
- Do not run further E0 controls (GTEx, second platform) — superseded by D2.

## 4. Escalate immediately
1. G4.1 positive control fails, or G2.6 overfit-one-batch fails.
2. **D2 returns I ≈ H** — the premise is in trouble.
3. **D1 returns `programme_only` > `programme_free`** — the collapse story is wrong.
4. D1's arms cannot be matched on epochs/LR/budget.
5. The channel vanishes when purity enters (D3) — a finding, but stop and report before building on it.

## 5. Self-audit
Max **3 agents including yourself**. After each result whose gates pass, spawn **one** adversarial
auditor instructed to **refute**, defaulting to "artifact" under uncertainty, returning
ARTIFACT / INCONCLUSIVE / SOUND. Third agent for tiebreak only. ARTIFACT verdicts are logged as
**defects, not findings**. Log gates to `GATE_LOG.md` and results — including negatives — to
`EXPERIMENT_LOG.md`. Commit and push per task once gates and audit both pass.

**GPU:** required for D1/D2 (this is the first thing on this project that genuinely needs it). D3 is CPU.
Ask for the SSH login when you need it. Runs must be launched from a **clean worktree** — E0's provenance
gate fails the run and quarantines output into `FAILED_*.json` otherwise.
