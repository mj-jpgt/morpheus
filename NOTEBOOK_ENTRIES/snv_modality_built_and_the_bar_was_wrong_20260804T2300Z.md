# The SNV block is built and is not broken — and the positive control's absolute bar was mine, not measured

**Logged:** 2026-08-04 23:00 UTC. **Pre-registration:**
`NOTEBOOK_ENTRIES/PREDECLARED_snv_modality_channel_20260804T2234Z.md`, committed `8f7ed63`, written
before the artifact existed. **Code at the run:** `17c76b4`, clean worktree
(`git_dirty: false` is recorded inside the artifact manifest).

**Read the awkward part first.** The pre-registered positive control **did not pass as written**. It
required balanced accuracy `> null p95` **and** `>= 0.25`; it returned **0.2265**. That is a fail, it
is recorded as a fail, and the bar is **not** being lowered after the fact. Everything below is what
the fail turned out to mean.

---

## 1. What was scoped, and the negative findings, first

For each missing modality: is the data real, and does it cover the patients we already have?

| modality | source | patients | ∩ cohort (6,427) | ∩ test (2,766) | verdict |
|---|---|---:|---:|---:|---|
| **SNV** | MC3 public MAF, on disk since 2026-07-11 | 10,224 sequenced | **6,108 (95.0%)** | **2,686 (97.1%)** | ready; built here |
| **CNV** | GISTIC2 gene-level thresholded, on disk; SNP6 segments also on disk | 6,277 | **6,275 (97.6%)** | **2,655 (96.0%)** | ready; **not built** |
| **proteomics (RPPA)** | UCSC Xena PanCanAtlas; **was not on disk**, downloaded here | 7,701 | **4,655 (72.4%)** | **1,972 (71.3%)** | acquirable; **not built** |
| **mass-spec proteomics** | — | — | — | — | **not available to this project** |

Negatives, plainly:

- **Proteomics is the weak one and should be described as such.** RPPA is 258 antibodies, not a
  proteome; 14.1% of its matrix entries are missing and a sample carries 222 antibodies on average;
  and it covers **25 percentage points fewer** of the held-out partition than SNV or CNV. Any
  eventual system spanning it must say "RPPA proteomics", never "proteomics".
- **There was no proteomics anything.** A `find` over `/lambda/nfs/geeg` (maxdepth 5, patch stores
  pruned) for `*rppa*`, `*proteom*`, `*protein*` returned two unrelated markdown files. The claim in
  `decision_iterate_past_pbs_and_p4_negatives` that "TCGA carries RPPA proteomics … without needing
  CPTAC" is **correct as a statement about TCGA and was wrong as a statement about this disk**: it
  had to be fetched. It now exists at
  `raw_tcga_v1/TCGA-RPPA-pancan-clean.xena.gz`, sha256 `2dca46e0…1868d`, 5.6 MB.
- **The SNV wiring that "already existed" was a session script.** `known_covariate_control.py` grades
  mutation *labels*, but the labels came from an uncommitted script
  (`e1_clinical_endpoints_from_morphology_20260804T1930Z.md` §6: "No code was modified by this work").
  The one persisted per-patient matrix, `data/tcga_v1/tcga_mc3_snv_binary.parquet`, records its input
  sha256 and **no variant filter at all**. The gap to a first-class SNV modality was smaller than
  starting from zero, but it was not "already wired".
- **No GPU was needed and none was used.** Everything here is CPU plus one 5.6 MB download.

---

## 2. What was built

`v2/build_snv_targets.py` (`4269e79`) → `morpheus_phase_d/data/snv_targets_v1.npz`, sha256
**`f09729027c76f5abf21ca02b6d1ebc7137b93e62183c708c38f1dfa08a13ced4`**, plus
`.manifest.json`. Schema is exactly `frozen_rna_targets.npz`'s, so `run_calibra` and
`known_covariate_control` read it **unmodified** — the SNV block goes through the same instrument,
the same `cancer + pooled TSS` adjustment and the same floors as RNA. A modality with its own
statistics would not be comparable to the one already measured.

| quantity | value |
|---|---:|
| MAF rows read | 3,600,963 |
| …in the cohort | 1,823,179 |
| …passing FILTER ∈ {PASS, wga, native_wga_mix} | 1,597,224 |
| …of which strict `PASS` only | 1,568,806 (the rescue tags add 1.8%) |
| …protein-altering after that | 1,007,421 |
| patients emitted | **6,108** (train 2,917 / val 505 / **test 2,686**) |
| patients excluded, named in the manifest, never imputed to wild-type | 319 |
| patients with >1 sequenced aliquot (unioned) | 57 |
| candidate genes seen mutated | 19,461 |
| genes kept (≥2% and ≥100 **development** patients) | **350** |
| targets | 352 real (350 genes + 2 burden) + 352 matched permuted controls |

Sanity anchors that were not tuned to: TP53 prevalence **0.3716**, KRAS **0.0678**, mean
protein-altering burden 164.9 per patient — all in their published pan-cancer neighbourhoods. The
test-partition count **2,686** reproduces the independent count in
`PREDECLARED_E1_clinical_endpoints_20260804T1745Z.md` exactly. No target column is constant in any
partition (checked in all three), so nothing fails `run_calibra`'s G1 guard for a hidden reason.

---

## 3. M0 — the positive control, and why its failure was the bar's fault

`v2/calibra/modality_block_control.py` (`17c76b4`): ridge-onto-one-hot lineage, argmax readout, one
fixed alpha, 5 folds, `seed=42`; chance **measured** by 100 lineage-label permutations through the
identical code path, never assumed to be 1/32.

| arm | block | balanced accuracy | measured null p95 | permutation p |
|---|---|---:|---:|---:|
| **M0 (graded)** | `snv_targets_v1` real targets, n=6,108, d=352 | **0.2265** | 0.0338 | 0.0099 (resolution floor) |
| must-fail | its own matched permuted controls, d=352 | **0.0282** | 0.0434 | 0.9901 |
| diagnosis A | `frozen_rna_targets` real targets, n=6,427, d=90 | **0.6554** | 0.0334 | 0.0099 |
| diagnosis A2 | same, restricted to the 6,108 SNV patients | **0.6373** | 0.0329 | 0.0099 |
| diagnosis B | **v1's** `tcga_mc3_snv_binary.parquet`, same 350 genes, same 6,108 patients | **0.1994** | 0.0335 | 0.0099 |
| diagnosis C | this build, `snv_gene` only, d=350 | **0.2017** | 0.0338 | 0.0099 |
| diagnosis C | this build, `snv_burden` only, d=2 | **0.0931** | 0.0315 | 0.0099 |

**The verdict, in the order the pre-registration requires.**

1. **M0 failed as written.** 0.2265 < 0.25. Recorded as a fail. Not relaxed.
2. **The block is nevertheless not broken**, and four independent readings say so, each through the
   same code path: it clears its measured null by **6.7×** at the permutation resolution floor; its
   matched permuted controls sit *below* their own null and correctly fail; an **independent parse of
   the same MAF scores lower (0.1994) than this one**, so this build is not the weak one; and the
   probe is demonstrably capable of 0.6554 on a known-good block, so the bar was not above the
   probe's ceiling.
3. **The failure is therefore in the bar, and the bar was mine.** 0.25 came from a literature
   intuition about mutation-based tissue-of-origin classifiers, which use richer features (position,
   context, signatures) and nonlinear models. It was never a measurement, and I had no measurement to
   set it from. The right bar is relative — to the measured null, or to a known-good block through
   the same probe — and setting it now, after seeing 0.2265, would be exactly the move the notebook's
   rule 3 forbids. **It is left failed**, and a corrected criterion must be written in a *new*
   pre-registration, before its next use, and never applied backwards to this run.
4. **Consequence, honoured:** the pre-registration says that if M0 fails, M1 "is not interpretable and
   will not be reported as a result about biology". M1's numbers are in §4 as **measurements**, with
   that sentence attached. They are not a finding about SNV and morphology and must not be quoted as
   one.

**Independently useful, and it is a real result about the data rather than about the bar:** the
declared variant filter changes what a `1` means and changes it in the expected direction. This build
calls **4.82%** of gene×patient cells mutated against v1's **6.55%**; the two parses agree on
**98.27%** of cells; and the extra calls v1 carries — silent, UTR, and FILTER-failing rows — cost it
**0.02** balanced accuracy. That is the concrete cost of an unrecorded filter, measured rather than
argued.

---

## 4. M1 — cross-modal measurement (NOT a result about biology; see §3.4)

`run_calibra` unmodified, `--partition test` (**n = 2,686**), `cancer + pooled TSS` cross-fitted,
108 confound columns, 84 sites kept at `min_site_count=10`, `n_components=16`, `n_draws=10`,
`n_permutations=200` (resolution 1/201 = 0.00498), `seed=42`, `--score-random-controls`, 352 real
targets and 352 matched controls. Output `results/m1_calibra_snv/` (8,694 rows).

| artifact | state | adjusted top-CCA | permutation null median / p95 | p | held-out top-CCA | detection floor | random-control exceedance |
|---|---|---:|---:|---:|---:|---:|---:|
| d2_h_seed42 | **wsi_biology** | **0.2657** | 0.1431 / 0.1661 | 0.00995 | 0.2114 | 0.10 | 0.000 |
| d2_h_seed42 | rna_biology | 0.3445 | 0.1585 / 0.1781 | 0.00498 | 0.3231 | 0.10 | 0.000 |
| d2_h_seed42 | full_biology | 0.3418 | 0.1570 / 0.1754 | 0.00498 | 0.3214 | 0.10 | 0.000 |
| d2_i_seed42 | **wsi_biology** | **0.2094** | 0.1445 / 0.1638 | 0.00498 | 0.1575 | 0.05 | 0.000 |
| d2_i_seed42 | rna_biology | 0.3158 | 0.1563 / 0.1749 | 0.00498 | 0.2817 | 0.20 | 0.000 |
| d2_i_seed42 | full_biology | 0.3178 | 0.1551 / 0.1750 | 0.00498 | 0.2846 | 0.20 | 0.000 |

Adjustment attenuation 0.813–0.852, i.e. in the ≈1 band P1 already reports; transmission floor 0.01
everywhere; **not one** of the 352 matched permuted-control columns reached the observed channel in
any state (`random_control_exceedance_fraction` = 0.000 in all six).

**Which of the three pre-declared endings this is: (a) — the one I said I expected least.** The
morphology-only state's channel sits above its within-cancer pairing null and above its detection
floor in both arms. I predicted (b). It is recorded that way round.

**And now the four reasons it must not be read as "morphology reads mutation".**

1. **M0 failed as written, so by the pre-registration's own rule this is a measurement, not a result
   about biology.** That rule was written before the numbers and is honoured here.
2. **This cannot separate an SNV channel from residual lineage, and the project has already measured
   that the residual is real.** `tcga_nonlinear_confound_probe_result_20260804T2100Z.md` found that
   after this exact `cancer + pooled TSS` adjustment an out-of-fold k-NN still recovers cancer type at
   3.4–4.9× chance. §3 of this entry shows the SNV block is itself strongly lineage-patterned
   (balanced accuracy 0.2265 against a 0.0338 null). Two lineage-carrying blocks correlating after an
   adjustment that demonstrably does not remove lineage is the expected reading, not the exciting one.
   Nothing here rules it out, and nothing here was designed to.
3. **It is not comparable to P1's RNA number as a ratio.** P1 §4.4 quotes 0.6052 for the same
   artifact/state/partition against the RNA block, but that block has **90** targets to this one's
   **352** and n = 2,766 to this one's 2,686. A top-CCA maximum grows with the number of directions
   available to maximise over, so more targets flatters SNV, and `NOTEBOOK.md` already records that n
   sets the null. Quoting "SNV is 44% of RNA" would be quoting across two capacities.
4. **The channel is not decomposed.** Whether it lives in the 350 gene columns, in the two burden
   columns, or in the lineage residual, is **untested**. §3's diagnosis C shows burden alone carries
   3× chance lineage information on its own, so the TMB-shaped explanation is live and unexamined.

`observed_above_floor` reads 1.0 in all six states. Recorded with its standing caveat: that flag's
history is in `observed_above_floor_is_broken_and_every_channel_clears_20260804T2115Z.md`, and the
shipped `calibra_protocol.json` still carries the note that floor units are single-direction
correlation and "are NOT comparable to the multivariate adjusted/held-out top-CCA". The ratios that
entry reads instead — held-out top-CCA over detection floor — are **2.11×** (d2_h wsi_biology) and
**3.15×** (d2_i wsi_biology), inside the 1.80–2.74× band it reported for the 13 RNA states, with d2_i
above it against a floor that landed at 0.05 rather than 0.10.

**Nothing here is offered to `validate_claim()` and no blocker is proposed for discharge.**

---

## 5. Ordered build plan for the rest

Fixed in the pre-registration before the first number, repeated here so it is not re-chosen after it.

1. **SNV — done to artifact + M0.** 97.1% test coverage, one file, no network.
2. **CNV — next, and expected to be the stronger of the two.** 96.0% test coverage; the GISTIC2
   gene-level thresholded matrix is a `[patient, gene]` matrix in the same shape as the SNV one, so
   `build_snv_targets.py` generalises with a different reader. Two declared changes are needed and
   neither is optional: (a) the raw file is keyed by **sample** (`TCGA-A5-A0GI-01`, 10,845 columns),
   and v1's per-patient collapse produced fractional values like −0.714, i.e. it *averaged aliquots
   and thresholds together* — a per-patient rule must be declared instead of inherited; (b) a signed
   −2…+2 column read by a correlation is treated as monotone when the biology is not, so
   amplification (`≥ 1`) and deletion (`≤ −1`) become two separate binary target families. Rationale
   for expecting it to beat SNV morphologically: aneuploidy and focal amplification have visible
   correlates (nuclear atypia, grade); point mutations largely do not.
3. **RPPA — last, with the caveat in its name.** 72.4% cohort coverage and 14% missing entries mean
   it needs a missing-data policy declared in advance, which the other two do not.
4. **Not proposed:** CPTAC (excluded by decision), anything dbGaP-tier, and any GPU work.

---

## 6. Files

- Pre-registration `NOTEBOOK_ENTRIES/PREDECLARED_snv_modality_channel_20260804T2234Z.md` (`8f7ed63`).
- Code `v2/build_snv_targets.py`, `v2/calibra/modality_block_control.py`;
  tests `v2/tests/test_build_snv_targets.py` (10), `v2/tests/test_modality_block_control.py` (8).
- Artifact `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/data/snv_targets_v1.npz`
  (+ `.manifest.json`), sha256 `f09729027c76f5abf21ca02b6d1ebc7137b93e62183c708c38f1dfa08a13ced4`.
- New raw input `/lambda/nfs/geeg/biorag3_persistent_20260711/raw_tcga_v1/TCGA-RPPA-pancan-clean.xena.gz`,
  sha256 `2dca46e09280d2c26f9f83e82ee98c8c41d1b896f2966918529e913bfad1868d`.
- Results `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/snv_modality/results/`
  (`m0_lineage_all/`, `m0_diagnosis/`, `m1_calibra_snv/`); isolated workspace
  `.../snv_modality/ws/morpheus` at `17c76b4`, clean.
- `claim_guards.py` and `claim_evidence.json` untouched. No claim is made by this entry.

**Suite at the run commit**, `pytest v2/tests tests -q` on the local checkout:
**1 failed, 623 passed, 1 skipped** in 149.30s. The one failure is
`test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree`, which names
`v2/build_causal_basis_targets.py` and `v2/perturbation_basis_common.py` — both introduced by another
agent's commit `134365c` while this work was in flight. Neither of the two modules added here
contains any SVD or effective-rank code, so the failure is not attributable to this work; it is
reported rather than worked around, and it is that other agent's to resolve.
(`test_p2_figures` did not error in this run: matplotlib is present in the local environment, unlike
the box venv where those 27–28 errors are the documented norm.)
