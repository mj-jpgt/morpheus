# Data/evaluation interfaces

These modules are deliberately model-independent and can be called by the
controller before a training/evaluation DAG starts.

* `canonical_registry.build_canonical_registry(sources)` returns one
  participant-level registry with modality availability and source IDs.
* `canonical_registry.build_tcga_cdr_outcomes(clinical)` creates OS/PFI labels;
  `survival_coverage(outcomes, patient_ids, split)` is the required coverage
  gate before a Cox probe.
* `curated_panel.freeze_curated_panel(gmt, measured_genes, training_signatures)`
  freezes the reviewed Reactome/KEGG panel and writes its JSON-compatible
  manifest. Only entries whose status is `eligible` are primary endpoints.
  `FROZEN_MECHANISM_PROGRAMMES` is the exact 15-name hidden-RNA vocabulary
  required by Qwen/card scoring; `build_discovery_inputs` writes it as
  `mechanism_programme_manifest.json` alongside target tables.
* `discovery_targets.build_rna_target_bundle(..., score_mode="competitive_rank")`
  creates per-target masked, development-train-standardised target tables.
* `baseline_exports.export_raw_hoptimus_baseline`,
  `export_ridge_alignment_baseline`, and `export_cca_alignment_baseline`
  write the shared frozen NPZ artifact schema.
  `python -m morpheus.v2.export_baselines` is the unattended CLI that derives
  exact uncapped per-patient H-Optimus mean/std vectors from the paired patch
  store and emits raw, Ridge, and CCA artifacts under one manifest.
* `paired_bootstrap.paired_patient_and_cancer_bootstrap` returns paired
  challenger-minus-teacher intervals for one held-out prediction vector.

All caller-provided arrays must be in the exact canonical patient order. Every
module fails on duplicates or alignment ambiguity rather than silently sorting
or dropping a patient.
