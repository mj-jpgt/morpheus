> # ⚠ ACTIVE TASK LIST IS NOW `HANDOFF_PHASE_D.md`
> The E-series (E0–E5) is **closed**. E0/E0b/E3 have run; **F2 is WITHDRAWN** (the identity head
> is the frozen MLP-CLIP teacher, 2.6e-04). Claims are enforced in `v2/calibra/claim_guards.py`.
> Read `HANDOFF_PHASE_D.md` §0 for the full context update before anything else.

# HANDOFF → Build Agent

**Read this fully before touching anything.** Your last context is the training runs (all seeds, all
baselines). Everything after that — an audit, a literature programme, and a full change of direction —
is new to you. This file is the shortest path to being useful without breaking something.

**Branch to work on: `research/rebase-vision`.** It contains all fixes plus the new work.
Other branches: `fix/biology-collapse` (collapse fixes), `paper/rank-collapse-diagnostic` (a finished
diagnostic paper). Do not start from `v2` or `main` — they predate everything below.

---

## 0. THE DESTINATION — read this before anything else

**We are not fixing an evaluation bug. We are building toward a Nature-tier paper.** Everything in
Milestones A–C is *instrumentation and evidence*; the paper is the thing at the end. If you only take
one thing from this file, take this section.

### What we are ultimately building
A **pathway-addressable tumour representation** whose individual biological programmes can be
*queried* — where the addressability is **earned from perturbation data**, its limits are **measured
rather than asserted**, and it **beats prior representations on a metric the standard benchmark cannot
see.**

### The architecture we are moving to (this is the real change)
The current model has a single dense 256-D "biology" vector supervised by regression onto ~50 Hallmark
scores. That design is what collapsed. The replacement is **two-stage**:

- **Stage A — interventionally-identified slots (cell line).** Replace the dense biology head with
  **`(batch, n_programme, D)` addressable slots**, identified under *real perturbations* using
  genome-scale Perturb-seq (Replogle: K562 GWPS 11,258 perturbations, K562 essential 2,285, RPE1 2,679
  — on disk at `PRISM/data/perturbseq/`). This coverage satisfies the **Squires ~1-intervention-per-node**
  regime, so **block/programme-group identifiability is defensible**. It does *not* satisfy Varıcı
  (2/node, no doubles), so **never claim per-dimension identifiability.**
- **Stage B — transfer to the frozen tumour trunk (patients).** Map the frozen WSI+RNA trunk into the
  Stage-A slot basis with a lightweight adapter. **Identifiability does NOT transfer** — a frozen
  *observational* trunk cannot inherit an interventional guarantee. What we claim here is **measured
  stability and measured transfer**, not a guarantee. The seam is itself an open question nobody has
  answered: *does identifiability learned in cell lines survive transfer to patient tumours?*

### The losses that change
- **Remove / down-weight:** the programme neighbour-KL + supervised-contrastive terms that pin a 256-D
  head to a ~50-D manifold. **This is the diagnosed collapse mechanism.**
- **Keep:** `losses.feature_decorrelation` — built and tested, +53 effective rank (2.1×, 3 seeds).
  Must standardise features first; on L2-normalised embeddings a raw covariance penalty silently no-ops.
- **Add:** conditional-prior / sparse-mechanism-shift identification (iVAE-style, metadata as the
  auxiliary variable) for Stage A.
- **Consider:** RNA-paired InfoNCE for the biology head — give it the rank-preserving contrastive
  signal that keeps the identity head healthy, instead of regression onto a low-rank target.

### The contributions we intend to claim
1. **C-I (method)** — interventionally-identified programme slots transferred to an observational
   tumour trunk, **with the transfer's identifiability loss measured rather than assumed**.
2. **C-II (evaluation)** — **CALIBRA**: a confound-certified addressability instrument that reports
   *what effect size an analysis would have missed*. Spike-recovery calibration has never been applied
   to certify a cross-modal morphology↔molecular claim. **This is the most defensible contribution and
   the one hardest to scoop.**
3. **C-III** — abstention over identified slots (answer / multi-interpret / refuse).
4. **C-IV** — a per-modality **encode-vs-retrieve** rule for structured molecular modalities.

### Where F2 fits — this is the point
**F2 is not the destination. F2 is the warrant for the architecture change.**
If molecular supervision demonstrably *degrades* the molecular channel, then the case for tearing out
the dense low-rank-regression head and replacing it with identified, addressable slots stops being a
design preference and becomes an **empirically forced move**. That is the difference between "we
refactored our model" and "we showed the standard way of supervising these models is harmful, and here
is what to do instead." F2 is the bridge from *audit* to *new architecture*.

### Honest ceiling (do not inflate this)
- **A–C alone → Nature Methods / Nature Biomedical Engineering.** The instrument + the defect finding.
- **+ D (the win) and the biology leg → Nature Cancer / Nature Medicine.**
- ***Nature* main is NOT reachable** and is not the target. That slot goes to scale-currency papers
  (Prov-GigaPath 1.3B tiles, CHIEF 60,530 WSIs); we cannot buy that on one A100 with open data.
  Aiming there is how this becomes a *Nature Communications* paper instead of a *Nature Methods* one.

---

## 1. What changed (the pivot, in one paragraph)

We stopped trying to make MORPHEUS V2 beat MLP-CLIP. Investigation showed **the benchmark cannot see
representation quality**: ~46–49% of WSI→molecular Pearson is cross-cancer cohort structure, and the
genuine within-cancer, random-control-adjusted signal is **~+0.07 for every method including the
baselines**. We also found the biology head had collapsed, built a fix that doubles its effective rank
— **and the benchmark score did not move at all** (0.1366 → 0.1367). So the project pivoted to building
**CALIBRA**, an instrument that reports *what effect size an analysis would have missed*, and then using
it to make a scientific claim. Do not re-litigate this; it is settled and documented.

## 2. Established results — treat as given, do not redo

| Finding | Value | Where |
|---|---|---|
| Biology head collapse | biology eff. rank ~38 vs identity ~191 (of 256) | `v2/research/rebase/` |
| Covariance-decorrelation fix works | rank 49.9→103.3 (+53, ~2.1×), 3 seeds | `paper/.../RESULTS.md` |
| …but changes nothing measurable | within-cancer specificity 0.1366 → 0.1367 | same |
| Benchmark is confounded | ~47% cross-cancer; genuine ≈ +0.07, method-invariant | same |
| **F1 (new)** channel survives full adjustment | **held-out CCA 0.477 vs 0.151 chance** | `v2/research/rebase/nature/PHASE1_RESULT.md` |
| **F2 (new)** supervision looks harmful | **`wsi_identity` 0.539 > `wsi_biology` 0.477** | same |

**F2 is the current scientific target**: the head trained *for* molecular programmes carries *less*
molecular signal than the head trained for retrieval.

## 3. What is already built and tested — reuse, do not rewrite

`v2/calibra/` (11/11 tests in `v2/tests/test_calibra.py`):
- `spectral.py` — `effective_rank` (single source of truth), `cca_spectrum`, **`heldout_top_cca`**
- `residualise.py` — `confound_design`, `cross_fitted_residuals` (fold-safe)
- `calibration.py` — `spike_recovery_curve` (recovery curve + detection floor), `permutation_null`
- `run_calibra.py` — Phase-1 CLI, emits the repo-standard task-row schema

Also reuse (do **not** reimplement): `v2/honest_metrics.py` (`macro_group_pearson`,
`control_adjusted_specificity`), `v2/paired_bootstrap.py`, `v2/discovery_targets.py`
(`build_matched_random_controls`, defaults to 20 draws), `v2/curated_panel.py`, `v2/contracts.py`
(`require_state` — never substitute an undeclared state), `v2/v21_evaluation.py`
(`_fit_predict_ridge`, `_select_alpha` — cancer-grouped CV).

## 4. THE MISTAKE LIST — every one of these has already bitten us

1. **Effective rank = SINGULAR values (Roy–Vetterli).** Using covariance eigenvalues (σ²) reports ~5–6
   instead of ~38 — a ~6× error that sat in a paper draft. Always call
   `morpheus.v2.calibra.spectral.effective_rank`. Never write a third copy.
2. **The workspace must be a real junction, not a copy.** A stale scratchpad copy silently validated
   old code for several runs. Before trusting any test result:
   `wc -l <workspace>/morpheus/v2/calibra/spectral.py` and confirm it matches the worktree.
3. **The loader's molecular targets are TRAIN-FOLD ONLY.** `data.hallmark` and `data._v2_programmes`
   are a **constant placeholder on the held-out test split** (std exactly 0) → every test correlation
   is NaN. Use `frozen_rna_targets.npz` (has real test variance). This cost hours.
4. **RNA-input states are circular.** `rna_biology`, `rna_identity`, `full_*` predicting RNA-derived
   targets is near-tautological. Reference only; **never** in a headline claim. They are also the
   correct *positive control* — if they don't show a strong channel, the pipeline is broken.
5. **The split is held-out-CANCER**: 14 dev / 21 test cancers, **disjoint**. BRCA, GBM, OV are
   *train-only*. TSS is perfectly nested in cancer type, so **cross-split site prediction is vacuous**
   (guaranteed failure, zero information). Site analyses must be *within-cancer*.
6. **Pool rare TSS sites.** ~600 codes, 145 singletons. A singleton dummy is effectively a per-patient
   indicator, so residualising it *deletes that patient* rather than adjusting. Use
   `--min-site-count 10` (600 → 75 sites).
7. **In-sample CCA is badly inflated** (pure noise reads >0.3). Use `heldout_top_cca` for any absolute
   number; the permutation null is what makes the *excess* meaningful.
8. **`frozen_rna_targets.npz` has only ONE matched control per target** (`RANDOM_CONTROL__T__0`, not
   `__00..__19`) and is **not reproducible by current code**. Regenerate with `controls_per_target=20`
   if you need a low-variance null, and say so if you don't.
9. **`v21_evaluation._control_comparison_rows` (~line 439) differences POOLED Pearson**, i.e. the
   confounded metric. The honest quantity is `macro_cancer_pearson`. Fixed on `feature/discovery-nl`;
   verify before relying on it.
10. **Corrupted data files — do not trust:** `morpheus/data/processed/master_patient_table.parquet`
    `cancer_type` contains **TSS codes**, not cancer acronyms; `data_cache/eval_labels/
    thorsson_panimmune.xlsx` is a byte-identical copy of the TCGA-CDR file (no Thorsson columns);
    `patient_feature_registry.parquet` has bare-integer `patient_id`s. Derive TSS from the barcode
    (`pid.split("-")[1]`) and take cancer type from `TCGA-CDR-SupplementalTableS1.xlsx`.
11. **Repo inconsistencies**: the split is written as 14/21, 11/22 *and* 11/21 in different files;
    cohort size appears as 6192 / 6427 / 6443. Pick one, document it, don't silently mix.
12. **Three citations in `v2/research/rebase/lit/*.md` are FABRICATED** ("Nguyen et al.
    Answer/Clarify/Abstain", "Decode-gLM", "VCBench (single-cell)"). Never cite from those ledgers
    without opening the actual paper. See `NEAR_COLLISIONS.md`.
13. **Never compare our numbers to the literature's.** Held-out-cancer macro Pearson ~0.13 is not
    comparable to a pooled-random-split paper's ~0.4. Different protocol.
14. **Emit unavailability, never silence.** Repo convention: a row with `metric="status"`,
    `value=NaN`, `note="unavailable_<reason>"`. A missing row reads as "not measured"; a silent
    substitution reads as a lie.

## 5. Your next tasks — **SUPERSEDED, see `HANDOFF_EXPERIMENTS_NOW.md`**

> **⚠ AUTHORITATIVE TASK LIST = `HANDOFF_EXPERIMENTS_NOW.md` (E0, E0b, E1–E5).**
> **⚠ MANDATORY HEALTH GATES = `HANDOFF_GATES.md`.** Governing rule: *a negative result is reportable
> only if the positive control passed in the same run.* Every failure mode listed there has already
> happened on this project once.
> The Milestone A1/A2/A3 framing below predates the CLD engine and **overlaps it**:
> **A1 ≡ E3** (same experiment, do it once), **A2 → E4**, **A3 → E5**.
> Read §5 below for *context on why* these matter, then execute the E-series. Do not run both.

## 5-legacy. Milestone A context (everything is already on disk, **no GPU needed**)

> ### KEY QUESTION to answer inside A1 — the anchoring caveat
> `z_identity` is **anchored on the frozen MLP-CLIP teacher with a learned residual of ≈0**
> (`anchor_residual_scale ≈ -0.001`), so identity is *essentially the teacher itself*. That means F2
> ("identity 0.539 > biology 0.477") may partly restate **"MLP-CLIP beats our biology head"** — which
> we already knew, and which would deflate F2 from a discovery to a restatement.
> **The test is free, because `diagnostic_programme_only_seed42.npz` has no anchor at all.**
> - If `programme_only`'s biology channel is still weak → the effect is about the **objective**. F2 holds.
> - If it is strong → F2 was an **anchoring artifact**; report that immediately and escalate.
>
> Report this comparison explicitly in the A1 output. Do not bury it.

**A1 — the objective ablation (highest value).** Three artifacts share architecture, data and split
and differ *only* in training objective:
`diagnostic_identity_only_seed42.npz`, `diagnostic_programme_only_seed42.npz`,
`diagnostic_full_seed42.npz` (in `discovery_evidence_v2/runs/v21_release_20260720_retry3_resume_safe/
artifacts/`). Run the CALIBRA channel measurement on all three.
**Pre-registered prediction: adding molecular supervision does NOT increase — and may reduce — the
measured molecular channel.**
*Before drawing any conclusion*, read each artifact's `manifest_json` and confirm epochs / LR / seed /
token budget match. If they differ materially, the ablation is **suggestive, not causal** — report that
plainly rather than overclaiming.

**A2 — encoder breadth.** Add the on-disk baselines (`raw_hoptimus_meanstd`, `ridge_alignment`,
`cca_alignment`, `mlp_clip_seed42`). Question: is the channel encoder-invariant? If a raw mean-pooled
H-Optimus baseline matches the trained model, that is the headline.

**A3 — reconcile 0.477 vs +0.07.** On the *same* residualised data, compute per-target within-cancer
specificity (`honest_metrics`) beside the multivariate channel. They measure different things; the gap
is itself the finding (per-target analysis dilutes a real multivariate channel). Needs the 20-draw
nulls from mistake #8.

### Run command that works today
```bash
# workspace must be a JUNCTION to the worktree (see mistake #2)
PYTHONPATH=$WS python -B -m morpheus.v2.calibra.run_calibra \
  --artifacts <artifact.npz> [<more.npz> ...] \
  --targets   discovery_evidence_v2/frozen_rna_targets.npz \
  --output    runs/<name> \
  --levels 0.0,0.02,0.05,0.10,0.20 --n-draws 5 --n-components 16 \
  --n-permutations 100 --min-site-count 10
```
Outputs `task_rows.csv` (per-metric rows) and `calibra_summary.json`.
**Note:** `heldout_top_cca` lands in `task_rows.csv`, *not* in the summary JSON — reading only the JSON
will make it look like NaN (this already fooled me once).

## 5b. The bigger picture — where F2 fits, and how we actually WIN

**Be clear-eyed: MORPHEUS V2 lost to MLP-CLIP.** Retrieval 0.060 vs 0.066; molecular prompting −0.021
global / −0.022 within-cancer; and on the honest metric every method ties at +0.07. We never beat the
baseline on anything. Nothing below is a claim that we secretly won.

**Why the confound still matters** (the objection "everyone uses TCGA, so it's a fair fight"):
a shared confound partly cancels for *ranking*, but (i) it does not boost all methods equally — it
rewards whoever best encodes cancer identity; (ii) it does not transfer to a new hospital, so it
mispredicts deployment; and (iii) **it is demonstrably blind** — we doubled biology effective rank
(49.9 → 103.3, 3 seeds) and the benchmark moved 0.1366 → 0.1367. A large representational change,
invisible. That is our own evidence, not borrowed.

**So the instrument is not a consolation prize — it is the measuring device that makes a win
provable.** We built CALIBRA because the benchmark could not tell us whether we had improved. The
endgame is to use it to demonstrate an improvement the old metric would have missed.

### The chain
```
F2 (does molecular supervision hurt?)   ← A1 answers this
      │  if YES
      ▼
D1  remove/replace the harmful supervision   →   D2 retrain   →   D3 measure with CALIBRA
      │                                                              │
      └── the fix is already built and tested ──────────────────────┘
```

### Milestone D — the win (needs a GPU; out of scope for you until A1 returns)
If A1 confirms the objective is the problem, the fix is already specified and partly built:
- **Drop or down-weight the low-rank programme supervision** — the neighbour-KL + supcon pin a 256-D
  head to a ~50-D Hallmark manifold (this is the diagnosed collapse mechanism).
- **Keep `losses.feature_decorrelation`** (built, tested; +53 rank / 2.1× over 3 seeds).
- **Consider RNA-paired InfoNCE for the biology head** — give biology the same rank-preserving
  contrastive signal that keeps identity healthy, instead of regression onto a low-rank target.

**The win condition, stated in advance (do not move these goalposts):**
1. beat MLP-CLIP on the **calibrated channel** (held-out CCA excess over permutation null), **and**
2. beat it on **within-cancer control-adjusted specificity** (the honest +0.07 metric), **and**
3. hold retrieval within **5%** of the anchor (the pre-existing success criterion).

Winning on the *confounded* pooled metric is explicitly **not** a goal. If we beat MLP-CLIP on the
calibrated channel while the old benchmark shows nothing, that is the strongest possible result —
it is simultaneously a method win and a demonstration that the field's metric missed it.

### The full milestone plan (A → B → C → D)
- **A — the F2 result** (on disk, no GPU): A1 objective ablation (+ the anchoring key question),
  A2 encoder breadth across the 4 baselines, A3 reconcile multivariate 0.477 vs univariate +0.07.
  *New:* `v2/calibra/ablation.py`, `v2/calibra/univariate.py`.
- **B — the audit frame** (makes A defensible): B1 N/I/C/B decomposition (Naive = zero-parameter
  cancer-mean for bulk / **per-slide mean** for spots; Identity = site oracle; Composition =
  capacity-matched; Biology = residual, all with bootstrap CIs). B2 multi-confound certificate
  {TSS within-cancer, dx-year, n_patches} with a random-direction floor, full-state ceiling, BH
  correction, mRNAsi as stratifier not gate; emit `robustness_index` + certified/failed/indeterminate.
  *Power limit: TSS evaluable in ~10 of 21 test cancers — report `n_cancers_evaluable`, never certify
  from <3.* *New:* `v2/calibra/decomposition.py`, `v2/calibra/certificate.py`.
- **C — spatial replication** (HEST-1k, ~1,200 paired H&E+ST, CC BY): re-run A1+A3 with
  spatially-localized targets to show the effect is a property of the channel, not an artifact of bulk
  RNA averaging. Ship the **per-slide-mean** zero-parameter baseline — absent from the HEST leaderboard
  and from the 46+ models the 2026 *Brief Bioinform* survey catalogues. *New:* `v2/calibra/hest.py`.
  **User action: download HEST-1k (+ `cellvit_seg` for B1).**
- **D — the win** (GPU): as above.

**Target paper:** *molecular supervision measurably degrades the molecular channel; the effect is
encoder-invariant, survives a calibrated confound audit, replicates at spot-level resolution, and
removing it yields a representation that wins on a metric the standard benchmark cannot see.*
Honest venue ceiling for A–C alone: Nature Methods / Nature BME. D is what lifts it.

**If A1 inverts the thesis** (supervision *helps*), the paper becomes "why the collapsed head still
works" — still publishable, but escalate before proceeding; the milestone chain changes.

## 6. Do NOT do these
- Do **not** retrain anything for Milestones A–B; they are frozen-artifact analyses.
- Do **not** chase a benchmark win over MLP-CLIP. That target is retired — the benchmark is confounded.
- Do **not** add regularisation broadly "to fix collapse". The specific, tested fix is
  `losses.feature_decorrelation` (must standardise features first — on L2-normalised embeddings a raw
  covariance penalty is a silent no-op).
- Do **not** report an absolute CCA without its permutation null and held-out counterpart.
- Do **not** invent or infer a citation. If you cannot open it, write COULD-NOT-VERIFY.

## 7. Escalate to the mastermind (do not decide alone)
- A1 shows supervision *helps* (inverts the paper's thesis — still publishable, but a redirect).
- The three ablation artifacts turn out not to be comparable.
- Any result that would change a claim in `PHASE1_RESULT.md`.
- Anything requiring a GPU, new data downloads, or an external cohort.

**Orientation docs, in reading order:** `v2/research/rebase/nature/PHASE1_RESULT.md` (current results) →
`v2/research/rebase/METHOD_PLAN.md` (the method + honest ceiling) →
`v2/research/rebase/nature/NATURE_THESIS.md` (why this direction) →
`v2/research/rebase/NEAR_COLLISIONS.md` (prior art, incl. the fabricated citations).

---

## CLAIM GUARDS — the caveats that produce good numbers while the claim is wrong

**Enforced in code: `v2/calibra/claim_guards.py`, tested in `tests/test_claim_guards.py` (15 tests).
Route every claim through `validate_claim()` before it is written up. An inadmissible claim emits
`metric="status", value=NaN, note="inadmissible_<code>"` — visible, never dropped.**

These are NOT completeness complaints. Partial coverage of biology is normal and is extended by
adding scales. Each item below is a case where a **causal attribution silently becomes wrong** while
the numbers look fine — which is the only failure mode a certified catalogue actually has.

| code | the silent failure | what discharges it |
|---|---|---|
| **`composition_attribution`** | The dictionary is built from **pure cell-line populations** — no immune cells, stroma, vasculature or architecture, which is much of what an H&E slide shows. In patients, tumour-intrinsic programmes **correlate** with composition (proliferative → immune-cold, mesenchymal → stroma-rich), so the fit succeeds numerically while the coefficients absorb the composition signal the dictionary cannot represent. The catalogue then reports "gene *g* is morphologically legible" when the truth is "*g* correlates with an infiltration pattern, and the pattern is what is visible." **Adding modalities multiplies this rather than fixing it.** | cell-of-origin label (tumour-autonomous / immune-mediated / stromal / composition-driven) from spatial or deconvolution evidence, **or** beating a capacity-matched composition baseline |
| **`purity_confound`** | TCGA bulk RNA is a **30–90% tumour mixture**; dictionary atoms are pure populations, so coefficients absorb purity — which is also one of the most visually obvious features on a slide. "Morphology predicts molecular state" has a trivial explanation that reproduces the result exactly. **Flagged open since Phase 1 and still not closed.** | purity / mRNAsi in the confound design, effect survives |
| **`sign_blind`** | `svdvals(Vaᵀ Vb)` is invariant to response sign: an **anti**-aligned effect scores identically to an aligned one. Separately, CRISPRi measures **loss** of function while tumours are largely **gain** of function. A discovery claim can come out inverted. | a signed statistic, plus explicit statement of the knockdown-vs-GoF asymmetry |
| **`proliferation_deflation`** | The responsive arm is selected on *having a detectable effect*, enriching for essential / ribosome / cell-cycle genes. If E0's alignment is **proliferation matching proliferation**, it is the most generic axis in cancer biology and deflates to near-nothing — while looking identical to a real finding. **Untested.** | re-run with proliferation regressed out or stratified; gap must survive |
| **`single_platform`** | K562 and RPE1 are different **lineages** but the **same Perturb-seq protocol**. A shared platform artifact replicates across them exactly as readily as shared biology. Cross-lineage agreement rules out "a K562 quirk", **not** "a Perturb-seq quirk". | a perturbation resource on a different platform, or scope the claim to "replicates across lineages within Perturb-seq" |
| **`no_external_cohort`** | Every morphology result so far is TCGA, with well-documented site/scanner effects. Confound removal was verified for cancer type, **not across cohorts**, so a site-specific artifact survives the current checks intact. | replicate in CPTAC / HEST-1k / equivalent |

**Pinned status (commit `24d1bff`): the E0 result is NOT yet an admissible transfer claim.**
It is blocked on `proliferation_deflation` and `single_platform`. `test_current_e0_result_is_not_yet_an_admissible_transfer_claim` fails if anyone discharges these, so the project state has to be updated deliberately rather than drifting.

**Ranked by how much each would cost us if it turns out to be the explanation:**
1. `proliferation_deflation` — cheapest to test, and the one most likely to be true. **Do this first.**
2. `purity_confound` — cheap, already overdue, and would explain the morphology result outright.
3. `composition_attribution` — needs the spatial stage; the reason cell-of-origin is non-negotiable there.
4. `sign_blind` — does not threaten E0's existence, but blocks every per-gene directional claim.
5. `single_platform` / `no_external_cohort` — scope limits, dischargeable by data we can obtain.
