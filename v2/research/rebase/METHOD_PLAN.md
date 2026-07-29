# MORPHEUS — Formalized Method Plan

*Written after: a 499-paper rebase sweep (15 lanes), an 8-candidate adversarial novelty vote, and a
targeted 109-agent gap-closing study on the two un-voted directions + identifiability feasibility.
This document states the **contribution** and the **methods** precisely enough to build and to defend.*

---

## 0. Executive statement

**What we are building.** A tumor representation whose biological programmes are *individually
addressable* — you can ask it, in natural language, about one programme and get that programme back
rather than a confounded bundle — where the addressability is **earned from perturbation data**, its
limits are **measured rather than asserted**, and the system **abstains** when a question falls
outside what it can address.

**Why (the wound we are repairing).** Our own result this session: the biology head silently collapsed
to low effective rank, a covariance term recovered it **+53 (~2.1×) across 3 seeds**, and the
benchmark score **did not move at all** (within-cancer specificity 0.1366 → 0.1367). The field's
standard evaluation is structurally blind to representation quality. Everything below follows from
taking that seriously.

**What we may NOT claim** (settled by the sweep — these are retired as headlines):
"an NL-promptable multimodal cancer model" (Med-PaLM M, PathChat, ChatNT); "pathway-addressable
tokens" as architecture (SurvPath); "perturbation as a query on a latent" (CPA/GEARS/biolord/scGen);
"encode-vs-retrieve" as a mechanism (kNN-LM/RAG/RETRO); "our model has emergent biological knowledge"
(Schaeffer mirage; GenePT/Kedzierska nulls).

---

## 1. The hard constraint that reshapes the architecture

**Finding (gap study, GAP 3).** Identifiability theorems require *interventional* data. A **frozen,
observationally-trained WSI+RNA trunk cannot inherit an identifiability guarantee** — not even a block
one. The honest ceiling for such a trunk is *cross-seed stability* (MCC / Roeder linear
identifiability). Block/group identifiability is defensible **only for the sub-block actually trained
under real interventions** (cf. PerturbedVAE Thm 1: responsive latents recovered up to
permutation+scaling; invariant latents only up to an arbitrary linear block map).

**Our interventional budget (verified on disk, `PRISM/data/perturbseq/`):**

| Dataset | Perturbations | Genes |
|---|---:|---:|
| K562 genome-wide (Replogle GWPS) | **11,258** | 8,248 |
| K562 essential | 2,285 | 8,563 |
| RPE1 | 2,679 | 8,749 |

→ Satisfies the **Squires regime (~1 intervention/node)** across much of the transcriptome; does
**not** satisfy Varıcı (2/node — no double perturbations at scale). Two cell lines give a genuine
cross-context stability test.

**Architectural consequence.** The method must be **two-stage and honest about the seam**:

- **Stage A (interventional, cell line):** learn programme slots under real perturbations → *block
  identifiability is claimable here.*
- **Stage B (observational, patient tumor):** align/transfer those slots onto the frozen WSI+RNA
  tumor trunk → *identifiability is NOT claimable here; what is claimable is **measured stability and
  measured transfer**.*

This seam is not a weakness to hide — **it is a research question no one has answered: does
identifiability learned in cell lines survive transfer to patient tumors?** We measure it rather than
assume it.

---

## 2. Contributions (stated as defensible claims)

**C-I (primary, method).** *Interventionally-identified programme slots, transferred to an
observational tumor trunk, with the transfer's identifiability loss measured rather than assumed.*
The novelty is not slots and not identification; it is the **two-stage construction plus the honest
seam measurement** (cross-seed MCC, Roeder linear identifiability, and a per-programme statement of
what survived transfer).

**C-II (primary, evaluation).** *A confound-certified addressability test.* A prompt for programme *P*
must return *P* and not its correlated bundle. Scored with: (i) a Dawood-style **co-dependence
leakage** measure, (ii) a **validity certificate** — the recovered direction must *fail* to predict
site/scanner while succeeding on biology, (iii) a decisive margin over a **GenePT text-prior** and a
**PCA/MOFA+** null. *This is the instrument that would have caught our rank collapse; it is the most
defensible contribution because evaluation is hardest to scoop.*

**C-III (secondary, interface).** *Abstention over identified slots.* Auto-detect the latent
biological task from an NL query; return multiple interpretations when under-specified; **abstain when
no identified slot addresses the question**. Distinguished from input-quality triage: our latent is the
*biological task*, over *identified* slots, evaluated on **held-out task families**.

**C-IV (secondary, modality routing).** *A per-modality encode / retrieve / marginalize rule for
structured molecular modalities.* See §4 for the reframe that keeps this alive.

---

## 3. Methods

### 3.1 Stage A — interventionally-identified slots (the identification engine)
- **Bias must be named** (Locatello: unsupervised disentanglement is impossible). Ours: *perturbation
  interventions* + *pathway-membership grouping* + *metadata-conditioned prior*.
- **Mechanism:** conditional-prior identification (iVAE-style, auxiliary `u` = cell line / batch /
  perturbation target) + **sparse-mechanism-shift** (sVAE+ / SAMS-VAE) over Replogle Perturb-seq.
- **Grouping hedge** (Morioka & Hyvärinen): known pathway membership gives group identifiability
  *without* interventions — our fallback where perturbation coverage is thin.
- **Causally-dependent latents** (CauCA), not independent slots — pathways interact.
- **Claim ceiling:** block/programme-group only, reported against the Squires budget we actually meet.

### 3.2 Stage B — transfer to the frozen tumor trunk
- Lightweight adapter (Q-Former-style query tokens) mapping the frozen WSI+RNA trunk into the Stage-A
  slot basis; trunk stays frozen (cheap, and preserves the fair baseline).
- Keep the **anti-collapse covariance term** we already built and validated (`feature_decorrelation`,
  +53 rank, 3 seeds) — without it slots collapse and addressability is meaningless.
- **Measure, don't assume:** per-programme cross-seed MCC, Roeder linear identifiability across runs,
  and slot-wise transfer fidelity from cell line → tumor.

### 3.3 The addressability / confound test (C-II)
- Co-dependence leakage: prompt for a single programme, measure contamination by its correlated bundle.
- **Validity certificate:** train site/scanner predictors on each recovered slot direction — require
  *failure* on site, *success* on biology (report a de Jong Robustness-Index-style scalar).
- Nulls that must be beaten decisively: GenePT text-prior, PCA/MOFA+ factors, and a random-programme
  control (Hewitt–Liang selectivity).
- **Amnesic causal-use gate:** erase a programme direction (LEACE/INLP) and confirm the downstream
  prediction *changes* — i.e. the programme is **used**, not merely decodable. Wire as the accept/reject
  rule on generated hypothesis cards.

### 3.4 The NL interface (C-III)
- Task text → conditions the query tokens → posterior over identified slots → route / multi-interpret /
  **abstain** below a calibrated threshold. Wire PRISM's unused `tqi.py` (`ScopeDetector` + head
  registry) onto real slots.
- Evaluate on **held-out task families** (FLAN-style) — the credibility test against memorized heads —
  plus calibration and accuracy-coverage curves.

### 3.5 Modality routing (C-IV)
- Per modality {proteomics, phospho, CNV, SNV, bulk RNA}: decide **encode / retrieve / marginalize**.
- Mandatory baselines: **In-Context-RALM** (zero-training retrieval arm), **GenePT** (text-prior arm).
- **Survives-distillation gate:** a "must-retrieve" verdict only counts if the gap persists after
  distilling the retrieval advantage into a better-trained encoder (Xu/Alon/Neubig).
- Train with modality dropout; serve any observed subset (MoPoE-style marginalization).

---

## 4. Novelty verdicts from the gap study (and the required reframes)

| Direction | Verdict | Prior art that constrains it | Surviving delta (the reframe) |
|---|---|---|---|
| **Modality Encodability Score** (C-IV) | **NOVEL ONLY IF REFRAMED** | He et al., JMLR 25 (2024) — per-modality utility + greedy submodular selection, **verified**: include-vs-exclude only, no retrieval option. Singh et al. 2604.00715 — encode-vs-retrieve crossover (D/N≈4.14) + substitutability σ, but over undifferentiated text tokens/model scale. Adaptive-RAG/Self-RAG/SeaKR gate per **query**, not per **source**. | A **per-modality** (not per-query, not include/exclude) routing rule for **structured molecular modalities with no natural text or image form** — the surveyed open problem — pre-registered, with the two named baselines and the distillation gate. |
| **Identification-buys-intervention** (part of C-I) | **NOVEL ONLY IF REFRAMED** | Zhang/Squires/Uhler (NeurIPS 2023) already state the theory *and* instantiate on Perturb-seq. | Only the **controlled ablation** survives: identified slots vs a **matched non-identified trunk** vs an **additive/linear mean-shift null**. Bar is brutal — Ahlmann-Eltze et al. (Nat. Methods 2025) report deep models losing to a trivial additive baseline on held-out doubles. *Do not lead with this.* |
| **Frozen-trunk identifiability guarantee** | **NOT CLAIMABLE as stated** | Identifiability requires interventions; observational trunks inherit none. | Restated as C-I: identify in Stage A, **measure** what survives transfer in Stage B. |

**Unresolved (honest gap).** The five "reframe-or-die" near-collisions — Nguyen triage (C-III),
SurvPath (C-I), Winter 2606.29949 (Stage B), Probing/Fusion 2606.17115 + VCBench (C-II), Decode-gLM /
amnesic probing (C-II gate) — **were NOT resolved**: no differentiation claim survived verification in
the gap study. Positioning against these five is a **prerequisite for writing**, not a detail. Resolve
by reading the five papers directly (not by search), since search has now twice failed to settle them.

---

## 5. What to build, in order

1. **C-II first (the instrument).** Cheapest, most defensible, and it makes everything else falsifiable.
   It also already has its motivating result (the rank-decoupling finding).
2. **Stage A slots** on Replogle (data verified present) with the claim ceiling fixed at block/group.
3. **Stage B transfer + seam measurement** — the genuinely open scientific question.
4. **C-III abstention interface** on top of real slots.
5. **C-IV modality routing** as a study once slots exist (needs CPTAC wiring, currently inventory-only).

**Non-negotiable guardrails:** cancer-held-out validation; per-sample missingness; earn-your-place
modality gating; every capability claim paired with the null it beats (text-prior, PCA, additive,
random-programme) and the confound it survives (site/scanner).

---

## 6. Honest risk register

- **The additive-baseline null (Ahlmann-Eltze)** may beat us on perturbation prediction. Mitigation:
  do not make interventional prediction a headline claim; report it as an ablation.
- **Identifiability may not survive Stage-B transfer.** This is a real possible negative result — and
  it is publishable as such *only if* C-II exists to measure it credibly.
- **The five unresolved near-collisions** could each shrink a contribution to a footnote.
- **Verification caveat:** parts of the gap study ran while the safety classifier was unavailable. The
  He et al. citation was directly verified; the Ahlmann-Eltze null is corroborated by an independent
  earlier sweep but was not re-verified (web budget exhausted). Re-verify before submission.
