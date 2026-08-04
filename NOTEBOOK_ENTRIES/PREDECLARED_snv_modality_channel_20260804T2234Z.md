# PREDECLARED — does an SNV block carry a channel at all, and at what floor?

**Written:** 2026-08-04 22:34 UTC (box clock, `date -u` on `150.136.45.194` = `2026-08-04T22:27:44Z`;
the local repo clock agrees to the minute). **Written before the SNV target artifact exists** and
therefore before any number in it exists. Scope: the *missing modalities* half of P4 — SNV, CNV and
proteomics — not the inductive-adjustment operator, which is a separate agent's task.

---

## 0. What was found on disk before this file was written (inventory, not measurement)

This section records file-existence and patient-ID set intersections only. No model, no probe, no
statistic that could come out for or against anything. It is stated here because the build below is
justified by it.

| modality | file on disk | patients | ∩ paired cohort (6,427) | ∩ test partition (2,766) |
|---|---|---:|---:|---:|
| SNV | `raw_tcga_v1/mc3.v0.2.8.PUBLIC (1).maf.gz` (753 MB) → `data/tcga_v1/tcga_mc3_snv_binary.parquet` (6,110 × 21,098, binary) | 6,110 | **6,108 (95.0%)** | **2,686 (97.1%)** |
| CNV | `raw_tcga_v1/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes (1).gz` (85 MB) → `data/tcga_v1/tcga_gistic_gene_cnv.parquet` (6,277 × 23,700) | 6,277 | **6,275 (97.6%)** | **2,655 (96.0%)** |
| CNV (segment) | `raw_tcga_v1/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg` (166 MB) | not parsed | — | — |
| proteomics (RPPA) | **absent from disk.** Downloaded 2026-08-04 22:32 UTC to `raw_tcga_v1/TCGA-RPPA-pancan-clean.xena.gz`, sha256 `2dca46e09280d2c26f9f83e82ee98c8c41d1b896f2966918529e913bfad1868d`, 5.6 MB, 258 antibodies × 7,754 samples | 7,701 | **4,655 (72.4%)** | **1,972 (71.3%)** |
| mass-spec proteomics | **none.** Not on disk, not acquirable for this project (CPTAC excluded by decision) | — | — | — |

`find` over `/lambda/nfs/geeg` (maxdepth 5, patch stores pruned) for `*rppa*`, `*proteom*`,
`*protein*` returned only two unrelated markdown research files. There was no proteomics code path
and no proteomics data; the RPPA file above is new to this project.

**Negative finding stated first:** proteomics has the worst coverage of the three by ~25 percentage
points, is the narrowest assay (258 antibodies, ~14% of matrix entries missing, mean 222 antibodies
non-missing per sample), and is the only one of the three that needed a new download. RPPA is *not*
mass-spec proteomics and a system built on it must never be described as spanning "proteomics"
without that qualifier.

---

## 1. What is being built

`v2/build_snv_targets.py` → `snv_targets_<version>.npz` + `.manifest.json`, in exactly the schema
`frozen_rna_targets.npz` uses (`patient_ids`, `cancers`, `split`, `target_names`, `target_groups`,
`scores` `float32 [patient, target]`, `metadata_json`), so it drops into `v2/calibra/run_calibra.py`
and `v2/calibra/known_covariate_control.py` without either being modified. Hashing follows
`v2/build_pbs_targets.py`: `_file_sha256` of every input, `_digest_array` of every emitted array,
`provenance.source_manifest()` for the code state.

**Declared construction choices, fixed now, before any output exists:**

1. **Source.** The MC3 public MAF is re-parsed from the `.gz` directly. The existing v1 parquet
   `tcga_mc3_snv_binary.parquet` is **not** inherited: its manifest records the input sha256 but
   records no variant-classification filter, so what "1" means in it is unknown. An unknown filter is
   not a modality.
2. **Variant filter.** Keep rows whose `Variant_Classification` is in
   {`Missense_Mutation`, `Nonsense_Mutation`, `Frame_Shift_Del`, `Frame_Shift_Ins`, `Splice_Site`,
   `In_Frame_Del`, `In_Frame_Ins`, `Nonstop_Mutation`, `Translation_Start_Site`} — i.e. protein-
   altering coding variants. Silent/UTR/flank/intron/IGR/RNA are dropped. Additionally require
   `FILTER` ∈ {`PASS`, `wga`, `native_wga_mix`, `PASS,wga`, `wga,native_wga_mix`} (MC3's own accepted
   set; the strict-PASS-only subset is recorded in the manifest as a count so the effect of this
   choice is visible).
3. **Patient identity.** `Tumor_Sample_Barcode` reduced to the 3-field canonical TCGA patient
   (`TCGA-XX-YYYY`), matching every other patient join on this project. A patient with more than one
   sequenced aliquot is the **union** of its variants; the duplicate count is recorded.
4. **Cohort and split.** `paired_split_maximal.json` (sha256
   `3c29cd98f534f1699a6c859aaabe671c7ce3e9cbf0572cdab2394fe2d884b29e`) is authoritative. Patients with
   no MC3 record are **excluded from the artifact and listed by name in the manifest**, never imputed
   to wild-type — absence of sequencing is not absence of mutation.
5. **Gene selection is fit on the development partition only.** Genes are ranked by mutation
   prevalence among `train`+`val` patients; a gene is kept if it is mutated in **≥ 2%** of
   development patients **and** at least 100 development patients. Test rows do not influence which
   genes exist. The kept-gene count is whatever that rule returns; it is not tuned afterwards.
6. **Targets.** Three groups:
   - `snv_gene` — binary protein-altering mutation status per kept gene.
   - `snv_burden` — two columns: `SNV_BURDEN_LOG1P` = log1p(protein-altering variant count) and
     `SNV_BURDEN_LOG1P_ALL` = log1p(all-classification variant count).
   - `random_control` — one `RANDOM_CONTROL__<target>__0` per real target, constructed by permuting
     that target's values across patients with `np.random.default_rng(42)`. Prevalence/marginal is
     preserved exactly; the patient pairing is destroyed. This is the same naming
     `run_calibra.py:364` already detects.
7. **No standardisation, no cohort-level centring** inside the builder. `run_calibra` residualises
   downstream; a second normalisation here would double-adjust.

---

## 2. What will be measured, and what each outcome means

Two measurements. **M0 is a precondition: if M0 fails, M1 is not interpretable and will not be
reported as a result about biology.**

### M0 — positive control: the SNV block must recover cancer type

Mutation profiles are strongly lineage-informative (BRAF/SKCM, APC+KRAS/COADREAD, VHL/KIRC). A
per-patient SNV matrix that cannot recover lineage is broken, not uninformative.

- **Statistic.** Out-of-fold multinomial-vs-cancer probe on the `snv_gene` + `snv_burden` columns of
  the artifact (controls excluded), 5-fold stratified on cancer, `seed=42`; read as **balanced
  accuracy** over the cancers present.
- **Null.** 100 cancer-label permutations, same folds, same code path; the 95th percentile is the
  chance level, measured not assumed.
- **Pass criterion, fixed now:** balanced accuracy **> null p95** *and* **≥ 0.25** in absolute terms.
- **Fail** means the builder is wrong. It gets diagnosed; it does not get relaxed.

### M1 — the real question: does the SNV block share a channel with the WSI representation, above the floor?

- **Machinery.** `v2/calibra/run_calibra.py`, unmodified, `--targets snv_targets_<version>.npz`,
  `--artifacts` the anchor block
  `morpheus_phase_d/runs/d2_final/artifacts/d2_h_seed42.npz` and `d2_i_seed42.npz`,
  `--partition test`, adjustment `cancer + pooled TSS (min_site_count=10)` cross-fitted, `seed=42`,
  `n_components=16`. No new statistic is invented for this modality.
- **Read out, all three, none of them alone:** the calibrated channel; its within-cancer pairing
  permutation null; and the `detection_floor` / `transmission_floor` pair, so "no channel" is
  distinguishable from "no power to see one".
- **The `RANDOM_CONTROL__` columns are scored in the same run** (`--score-random-controls`). A
  channel that the matched controls also clear is not a channel.

### The three endings, written before the run

- **(a) Channel above null and above the detection floor.** The SNV block carries morphology-readable
  structure at block level. Given §0 of `e1_clinical_endpoints_from_morphology_20260804T1930Z.md` —
  where BRAF, KRAS, APC, PIK3CA and dMMR all covered chance after cancer-type adjustment, and the
  four endpoints that did clear were flagged as possibly one atypia/TMB axis — **this is the outcome
  I expect least**, and if it lands it is most likely the burden columns, i.e. the same TMB-like axis
  under a new name, not per-gene mutation legibility. It would be reported that way.
- **(b) Channel at or below the null.** The honest and, on the prior evidence, most likely ending.
  Reported as: SNV is buildable and clean (M0 passed) but does not read out of this WSI
  representation at block level after cancer adjustment. This does **not** license "SNV is not a
  useful modality" — it licenses "SNV is not readable *from morphology* in this representation",
  which is a different claim and is the only one that would be written.
- **(c) Channel above null but below the detection floor.** Undecidable at this n. Reported as
  undecidable, with the floor quoted, not as a negative.

**What is explicitly NOT being claimed by this work regardless of outcome:** nothing about whether a
promptable multimodal system works, nothing about CNV or proteomics (not built here), and no paper
claim of any kind. `claim_guards.validate_claim()` has not been consulted because no claim is being
made; the outputs are an artifact and two measured numbers.

---

## 3. Ordered build plan for the remaining modalities (scoped, not built here)

Recorded now so the ordering is on record before the first result, not chosen after it.

1. **SNV — first.** Best-covered raw source of the three (97.1% of the test partition), already
   parsed once on this project, one file, no network. Built here.
2. **CNV — second.** 96.0% test coverage, already on disk in two forms (GISTIC2 gene-level thresholded
   and SNP6 segments), and the gene-level table needs no new method — it is a `[patient, gene]`
   matrix in the same shape as the SNV one, so `build_snv_targets.py` generalises to it with a
   different reader and a different discretisation rule (declared: amplification `≥ 1` and deletion
   `≤ −1` as two separate binary target families, because a signed single column would be read by a
   correlation as monotone when the biology is not). Expected to be the *stronger* of the two
   morphologically: aneuploidy and focal amplification have visible correlates (nuclear atypia,
   grade), where point mutations largely do not. Same M0/M1 protocol.
3. **Proteomics (RPPA) — third, and with a caveat attached to its name.** 72.4% cohort coverage,
   258 antibodies, 14% missing entries, and it is antibody-based, not mass-spec. It needs a
   missing-data policy (declared in advance, not chosen from the results) that the other two do not.
   Do it last, and never call the result "proteomics" unqualified.
4. **Not proposed:** CPTAC (excluded by decision), any dbGaP-tier resource, and any new GPU work. All
   three modalities above are CPU + one 5.6 MB download.

---

## 4. Files this predeclaration binds

- Builder: `v2/build_snv_targets.py` (does not exist at the time of writing).
- Artifact: `morpheus_phase_d/data/snv_targets_v1.npz` + `.manifest.json` (do not exist).
- M0: `v2/tests/` coverage for the builder, and the probe run recorded in the result entry.
- M1: `run_calibra` output directory, recorded in the result entry.
- Nothing in `claim_guards.py` or `claim_evidence.json` is touched by any of it.
