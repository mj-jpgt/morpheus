# V2.1 recovery runbook

`run_v21_recovery_on_lambda.sh` is the only supported unattended entry point
for this recovery run.  It runs a strict-core WSI+RNA experiment; clinical,
SNV, CNV, protein, and spatial inputs are intentionally not silently enabled.

## Required inputs

Set these paths on Lambda.  Paths are recorded in the run source manifest;
the script never reads or prints a Hugging Face token.

```bash
export ROOT=/lambda/nfs/geeg/biorag3_persistent_20260711
export CODE="$ROOT/code"
export ANCHOR_ARTIFACT=/absolute/path/mlp_clip_seed42.npz
export LEGACY_BASELINES_FILE=/absolute/path/legacy_baselines.txt
export TEACHER_BY_SEED_FILE=/absolute/path/teacher_by_seed.json
export RNA_TABLE=/absolute/path/patient_by_gene_rna.tsv
export GMT=/absolute/path/reviewed_msigdb_reactome_kegg.gmt
export TRAINING_GMT=/absolute/path/v2_training_programmes.gmt
export CLINICAL_CDR=/absolute/path/TCGA-CDR-SupplementalTableS1.xlsx
```

`legacy_baselines.txt` contains one split-matched artifact path per line.
`teacher_by_seed.json` is a JSON mapping such as
`{"42":"/…/mlp_clip_seed42.npz","43":"/…/mlp_clip_seed43.npz","44":"/…/mlp_clip_seed44.npz"}`.
All artifacts must have exactly the canonical patient/cancer/split universe;
the evaluator refuses partial or differently split artifacts.

Each seed-matched MLP-CLIP teacher must itself be a final refit on
development train+validation patients and declare
`fit_population: "development_train_val"` in `manifest_json`. This is a hard
gate because anchored V2 cannot be called a fair final comparison when its
teacher was trained on fewer development patients. Set
`REQUIRE_FINAL_TEACHERS=0` only for an explicitly exploratory dry run.

The Lambda environment must already contain `torch`, `numpy`, `pandas`,
`pyarrow`, `scikit-learn`, `scikit-survival`, and, when `ENABLE_QWEN=1`,
`transformers`. The controller verifies these before writing training outputs;
missing survival/Qwen dependencies stop immediately instead of producing a
partial report.

## Launch

```bash
cd "$CODE"
bash -n morpheus/v2/run_v21_recovery_on_lambda.sh
nohup bash morpheus/v2/run_v21_recovery_on_lambda.sh \
  > "$ROOT/runs/v21_recovery_11v21/launch.log" 2>&1 &
```

Monitor `runs/v21_recovery_11v21/state/status.json`,
`logs/controller.log`, and `logs/gpu_monitor.csv`.  Token budget begins at
65,536 and is halved only after a detected CUDA OOM, retaining every patch and
resuming from the checkpoint.  GPU utilisation is observed rather than
artificially inflated.

## Pipeline and outputs

The controller performs, in order:

1. Compile/split/runtime preflight and source manifests.
2. Canonical RNA targets, curated panel, random controls, and CDR outcomes.
3. Raw pooled H-Optimus, Ridge, and CCA baseline exports.
4. Development-train-only slide self-supervision.
5. Validation-only identity/programme/full profile selection.
6. Final seed-matched anchored and optional no-anchor V2 training.
7. Outer molecular prompting, biology-neighbour retrieval, development-
   thresholded molecular phenotypes, WSI/RNA/fused survival references,
   random-control, and paired teacher comparisons.
8. Closed-RAG Qwen card rendering and hidden-RNA scoring, labelled
   exploratory.

The default uses five train-only matched random controls per biological target
to keep the complete multi-method run tractable. Override
`CONTROLS_PER_TARGET` only if the measured controller throughput supports a
larger number; no target or patient is removed by this setting.

The mandatory deliverables are under `final/evaluation/`:

- `task_rows.csv` — metric/status rows, including explicit unavailable states;
- `molecular_predictions.parquet` and `survival_risks.parquet` — patient-level
  paired outputs;
- `evaluation_protocol.json` and `source_manifest.json` — cohort, fitting,
  bootstrap, and code provenance;
- `../qwen_source_cards.jsonl`, `../qwen_cards.jsonl`, and
  `../qwen_evaluation.json` — exploratory cards only.

Do not make a biology, survival, or generative claim from a run that fails the
profile promotion gate, the curated-vs-random-control comparison, or the
separate card quality gates.
