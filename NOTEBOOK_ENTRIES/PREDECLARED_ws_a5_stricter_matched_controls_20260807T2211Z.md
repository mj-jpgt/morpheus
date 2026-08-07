# PREDECLARED — WS-A5, stricter matched gene-set controls

**Written before `discovery_targets.build_matched_random_controls` is extended and before any new
matching mode has been run.** Agent W1-B, executing `paper/P1_REVISION_SPEC.md` §5 WS-A5
(complaint #7). Read in full: `v2/discovery_targets.py`, `v2/build_discovery_inputs.py`,
`v2/curated_panel.py`.

## 0. Provenance of the number being challenged — verified against a real artifact and a real hash
(PROJECT_GUIDE §2.7), BEFORE any new measurement

The training box (`ssh p1box`) is reachable this session (reprovisioned; verified with
`ssh p1box hostname` and an idle-GPU check). The 76–82% figure in `paper/P1_CALIBRA_DRAFT.md` §4.8.2
traces to:

* `run_calibra_ledger.sh` on the training box, which calls
  `python -m morpheus.v2.calibra.run_calibra --artifacts d2_h_seed42.npz d2_i_seed42.npz
  --targets frozen_rna_targets.npz --score-random-controls ...` — the exact script that produced
  `p1_evidence/track1/calibra_frozen_rna/calibra_gates.json`.
* **Hashes** (SHA-256, computed on the training box this session):
  `d2_h_seed42.npz = 4a18b94f1017b85dd576f30ee8e3caf92d7897630a7054efb70166191cbe69e3`,
  `d2_i_seed42.npz = 028e8635465dd3c6d3dbead25a8c204ca1ae0cee4aabb20e5412847fb147b665`,
  `frozen_rna_targets.npz = d526a7adc7456ac4f0e5e3ff71c0ef2bac96dc8488435ea714ba9840d8b51fb2` — identical
  across the three on-disk copies found (`data/`, `e0_run/data/`, `runs_misc/calibra_run/artifacts/`),
  so there is no `d2_h_seed42`-style multi-file-same-name drift here (the incident PROJECT_GUIDE §2.7
  itself warns about).
* **Recomputed directly from `calibra_gates.json`'s `random_control_verdicts`, independent of the
  draft's own arithmetic:** d2_h::full_biology 0.7624, d2_h::rna_biology 0.7615, d2_h::wsi_biology
  0.7698, d2_i::full_biology 0.8187, d2_i::rna_biology 0.8122, d2_i::wsi_biology 0.7588 — matches the
  draft's quoted 76–82% range exactly. **The claim is verified as traced to a real artifact and stands
  as stated; this predeclaration is about what happens to it under stricter matching, not a correction
  to it.**
* `frozen_rna_targets.npz` holds 90 real (curated Hallmark/KEGG_MEDICUS/custom-mechanism) target
  columns and 90 `RANDOM_CONTROL__` columns (one per real target — `controls_per_target=1` in the
  shipped build), n=6,427 patients, matched via the **existing** mean/log-variance/PC1 rule.

## 1. What is being extended, and how

`build_matched_random_controls` gains an additive `matching_mode` keyword (default unchanged —
`"mean_logvar_pc1"`, byte-identical to current behaviour, verified by a regression test before any
new mode is trusted). Three new modes, declared before running:

* **`multi_pc`** — match on the first `m=5` PC loadings (declared in advance, not tuned) instead of
  PC1 alone, in the same median/MAD-standardised nearest-neighbour scheme.
* **`coexpression_degree`** — match on each candidate gene's degree in a thresholded gene–gene
  Pearson-correlation graph computed on the **training-only** expression matrix, threshold `|r| >= 0.5`
  (declared in advance). A control gene embedded in the same co-expression neighbourhood density as the
  real target's genes is harder to distinguish on a network-structure statistic than one matched only
  on marginal mean/variance/PC1.
* **`pathway_overlap`** — match on each candidate gene's membership count across the same background
  pathway collection (the discovery GMT, `msigdb_discovery_2024.1.Hs.gmt` on the training box) used to
  build the real targets, so a control drawn from a gene that co-occurs in many of the same pathways as
  the target is not silently favoured or disfavoured relative to genuinely pathway-independent genes.

All three add covariates to the existing standardised nearest-neighbour draw; none replaces the
existing mean/log-variance/PC1 covariates, so `multi_pc`/`coexpression_degree`/`pathway_overlap` are
**supersets**, strictly stricter than the shipped rule, matching the spec's framing of "stricter", not
an alternative, weaker rule.

## 2. Real-data re-measurement protocol

* **Real target definitions**: the exact 90 real target names read directly from the shipped
  `frozen_rna_targets.npz` (not re-derived from panel-eligibility logic, which removes one source of
  possible drift). Gene membership for `HALLMARK_*`/`KEGG_MEDICUS_*` names comes from the real
  `msigdb_discovery_2024.1.Hs.gmt` on the training box; the remaining custom `immune_*`/`state_*`/
  `stroma_*`/`tgfb_emt` names come from `v2/curated_panel.py`'s `FROZEN_MECHANISM_PROGRAMMES`
  (in-repo, not re-derived either).
* **Real RNA table**: `tcga_pancan_rna.parquet` on the training box (same file the original build
  used, confirmed present at `morpheus_phase_d/e0_run/data/` and `morpheus_phase_d/data/`).
  **Real metadata** (patient_id, cancer, split): taken directly from `frozen_rna_targets.npz`'s own
  `patient_ids`/`cancers`/`split` arrays, so the fitting population is identical to the one that
  produced the number under test.
* **Fresh baseline + new arms, same `controls_per_target`.** Because the shipped build used
  `controls_per_target=1` (one draw per target — a small sample for a *median* statistic), and because
  a fair A/B between the existing rule and three new stricter rules needs both computed the same way,
  this predeclaration regenerates a **fresh baseline** (`mean_logvar_pc1`, `controls_per_target=10`,
  seed 42) alongside the three stricter arms (`controls_per_target=10` each, same seed) rather than
  reusing the shipped one-draw-per-target control block. The freshly regenerated baseline is checked
  against the shipped 76–82% figure as its own counterfactual (§4) before any stricter-mode number is
  trusted.
* **Readout**: identical statistic and identical grading as the shipped figure —
  `run_calibra.py --score-random-controls`, `heldout_single_direction_correlation`-based
  `fitted_direction` verdict, `d2_h_seed42`/`d2_i_seed42`, states `wsi_biology`/`rna_biology`/
  `full_biology`, test partition — run unmodified, not reimplemented, against the newly built target
  bundles.

## 3. Predeclared expectation and falsifier

**Expected (per spec, and per intuition): stricter matching pushes the ratio UP**, not down — a control
that is harder to distinguish from the real target on more covariates should read *closer* to the real
target, not further, because non-specificity is about how much of the legibility any gene set with
similar first- and second-order statistics carries, and the new covariates make that similarity
stronger, not weaker.

**Falsifier / the more interesting result, stated in advance and NOT to be explained away if it fires:**
if any stricter mode moves the control/real ratio **DOWN** relative to the freshly regenerated
baseline, in a majority of the six (2 artifacts x 3 states) cells, that is reported as the headline
finding of this work item — it would mean the *stricter* covariates (co-expression degree, pathway
membership) are themselves picking out something informative that the original mean/logvar/PC1 rule
was diluting by drawing from a wider, noisier candidate pool, which is a different and more interesting
claim about *why* the non-specificity exists than "any 3-moment-matched gene set reproduces it."

## 4. Counterfactuals (rule §2.6)

* **Regression / must-match control**: the extended function called with the default
  `matching_mode="mean_logvar_pc1"` and no other arguments must reproduce byte-identical control gene
  sets to the current unextended function on the same inputs (same digest) — checked by unit test
  before any stricter mode is trusted.
* **Positive control**: the freshly regenerated baseline arm, run through the real
  `run_calibra --score-random-controls` pipeline, must land inside (or very close to) the shipped
  76–82% band. If it does not, the freshly built target bundle does not match the original construction
  closely enough to trust a stricter-mode comparison against it, and that mismatch is reported before
  any stricter-mode ratio is quoted as a finding.
* **Must-fail / sanity control**: a `matching_mode` that matches on a covariate provably unrelated to
  the real target's identity (a pure random re-shuffle of the standardised covariate rows before
  nearest-neighbour matching) must **not** systematically raise the ratio — if a shuffled "stricter"
  control raises the ratio as much as a genuinely stricter one, the effect is a property of the matching
  machinery (e.g. pool-size shrinkage) rather than of the added covariates, and this must be reported
  as a confound in the new instrument itself.

## 5. Reporting rules

* All six cells (2 artifacts x 3 states) reported for every arm, no cell dropped for being unfavourable.
* Bad news (falsifier firing, or the freshly regenerated baseline missing the shipped band) leads the
  result entry.
* `discovery_targets.build_matched_random_controls` is extended additively; existing callers
  (`build_discovery_inputs.py`, `v2/tests/test_discovery_targets.py`) must keep passing unchanged.
