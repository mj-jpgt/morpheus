#!/usr/bin/env bash
# Autonomous strict-core V2.1 recovery DAG.  It is intentionally fail-closed:
# invalid targets, unavailable baselines, or a failed diagnostic promotion gate
# stop before the expensive three-seed stage.
set -Eeuo pipefail

ROOT="${ROOT:-/lambda/nfs/geeg/biorag3_persistent_20260711}"
CODE="${CODE:-$ROOT/code}"
RUN="${RUN:-$ROOT/runs/v21_recovery_11v21}"
PY="${PY:-$ROOT/.venv-morpheus/bin/python}"
CFG="${CFG:-$CODE/morpheus/configs/v1.json}"
SOURCE_SPLIT="${SOURCE_SPLIT:-$CODE/morpheus/data/processed/splits/tumor_state_heldout_cancer.json}"
SPLIT="${SPLIT:-$RUN/state/paired_split.json}"
EXPECTED_DEVELOPMENT_CANCERS="${EXPECTED_DEVELOPMENT_CANCERS:-11}"
EXPECTED_HELDOUT_CANCERS="${EXPECTED_HELDOUT_CANCERS:-21}"
LEGACY_BASELINES_FILE="${LEGACY_BASELINES_FILE:-}"
RAW_RNA_TABLE="${RAW_RNA_TABLE:-$ROOT/raw_rna/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp (1).tsv}"
RNA_TABLE="${RNA_TABLE:-$RUN/state/pancan_patient_by_gene.parquet}"
GMT="${GMT:?Frozen MSigDB GMT containing reviewed Reactome/KEGG signatures}"
TRAINING_GMT="${TRAINING_GMT:?GMT used for V2 programme training, for overlap audit}"
CLINICAL_CDR="${CLINICAL_CDR:?TCGA-CDR clinical table}"
SEEDS="${SEEDS:-42 43 44}"
DIAGNOSTIC_EPOCHS="${DIAGNOSTIC_EPOCHS:-12}"
MLP_DIAGNOSTIC_EPOCHS="${MLP_DIAGNOSTIC_EPOCHS:-40}"
MLP_FINAL_EPOCHS="${MLP_FINAL_EPOCHS:-40}"
BASELINE_GPU_CONCURRENCY="${BASELINE_GPU_CONCURRENCY:-2}"
PRETRAIN_EPOCHS="${PRETRAIN_EPOCHS:-10}"
PRETRAIN_CHECKPOINT_OVERRIDE="${PRETRAIN_CHECKPOINT_OVERRIDE:-}"
TOKEN_BUDGET="${TOKEN_BUDGET:-65536}"
MIN_TOKEN_BUDGET="${MIN_TOKEN_BUDGET:-16384}"
BOOTSTRAP_REPEATS="${BOOTSTRAP_REPEATS:-2000}"
CONTROLS_PER_TARGET="${CONTROLS_PER_TARGET:-5}"
MIN_ELIGIBLE_CURATED_TARGETS="${MIN_ELIGIBLE_CURATED_TARGETS:-12}"
MIN_MECHANISM_TEST_COVERAGE="${MIN_MECHANISM_TEST_COVERAGE:-0.80}"
RUN_NO_ANCHOR="${RUN_NO_ANCHOR:-1}"
# Qwen is not part of the scientific pipeline until the structured-card gate
# has passed.  It must be opted into with an explicit gate artifact.
ENABLE_QWEN="${ENABLE_QWEN:-0}"
CARD_QUALITY_GATE="${CARD_QUALITY_GATE:-}"
QWEN_MODEL_ID="${QWEN_MODEL_ID:-Qwen/Qwen2.5-7B-Instruct}"

mkdir -p "$RUN/logs" "$RUN/state" "$RUN/artifacts" "$RUN/diagnostic" "$RUN/final"
cd "$CODE"
export PYTHONPATH="$CODE:${PYTHONPATH:-}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export HF_HUB_CACHE="${HF_HUB_CACHE:-$ROOT/hf_cache/hub}"
export HF_DATASETS_CACHE="${HF_DATASETS_CACHE:-$ROOT/hf_cache/datasets}"
export TOKENIZERS_PARALLELISM=false
exec > >(tee -a "$RUN/logs/controller.log") 2>&1

status() {
  "$PY" - "$RUN/state/status.json" "$1" "${2:-}" <<'PY'
import json, sys, time
open(sys.argv[1], "w", encoding="utf-8").write(json.dumps({"stage": sys.argv[2], "detail": sys.argv[3], "time": time.strftime("%FT%T%z")}, indent=2))
PY
  echo "[$(date -Is)] $1 ${2:-}"
}

monitor() {
  echo "time,gpu_name,mem_used_mb,mem_total_mb,gpu_util_pct,power_w" > "$RUN/logs/gpu_monitor.csv"
  while true; do
    nvidia-smi --query-gpu=timestamp,name,memory.used,memory.total,utilization.gpu,power.draw --format=csv,noheader,nounits >> "$RUN/logs/gpu_monitor.csv" || true
    sleep 10
  done
}
monitor & MONITOR_PID=$!
trap 'kill "$MONITOR_PID" 2>/dev/null || true' EXIT

train_with_oom_retry() {
  local output="$1" profile="$2" seed="$3" epochs="$4" anchor="$5" pretrain="$6" fit_development="${7:-0}"
  local budget="$TOKEN_BUDGET" log="$RUN/logs/$(basename "$output").log"
  while true; do
    local command=("$PY" -m morpheus.v2.runner --data-config "$CFG" --split-file "$SPLIT" --output-dir "$output" --epochs "$epochs" --token-budget "$budget" --seed "$seed" --objective-profile "$profile" --pretrain-epochs 0 --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS")
    if [ -n "$anchor" ]; then command+=(--mlp-clip-anchor "$anchor"); fi
    if [ -n "$pretrain" ]; then command+=(--pretrain-checkpoint "$pretrain"); fi
    if [ "$fit_development" = "1" ]; then command+=(--fit-development); fi
    if [ -f "$output/last.pt" ]; then command+=(--resume "$output/last.pt"); fi
    set +e
    "${command[@]}" 2>&1 | tee -a "$log" >&2
    local result=${PIPESTATUS[0]}
    set -e
    if [ "$result" -eq 0 ]; then
      tail -n 1 "$log"
      return 0
    fi
    if grep -qi "cuda out of memory" "$log" && [ "$budget" -gt "$MIN_TOKEN_BUDGET" ]; then
      budget=$((budget / 2))
      status oom_retry "$(basename "$output") token_budget=$budget" >&2
      continue
    fi
    return "$result"
  done
}

export_checkpoint() {
  local checkpoint="$1" output="$2" anchor="$3"
  local command=("$PY" -m morpheus.v2.export --data-config "$CFG" --split-file "$SPLIT" --checkpoint "$checkpoint" --output "$output" --token-budget "$TOKEN_BUDGET" --diagnostics-output "$output.diagnostics.json")
  if [ -n "$anchor" ]; then command+=(--mlp-clip-anchor "$anchor"); fi
  "${command[@]}"
}

dedupe_paths() {
  declare -A seen=()
  local path
  for path in "$@"; do
    if [ -z "${seen[$path]:-}" ]; then
      seen[$path]=1
      printf '%s\n' "$path"
    fi
  done
}

status preflight "compile, source, split, and device checks"
test -x "$PY"; test -f "$RAW_RNA_TABLE"
PYTHONPYCACHEPREFIX=/tmp/morpheus_pyc "$PY" -m py_compile morpheus/v2/{runner,export,export_baselines,refit_mlp_clip,prepare_pancan_rna,build_discovery_inputs,v21_evaluation,select_v21_profile,build_hypothesis_cards,run_qwen_cards,evaluate_qwen_cards}.py
"$PY" - <<'PY'
import importlib
for package in ("numpy", "pandas", "pyarrow", "sklearn", "sksurv", "torch"):
    importlib.import_module(package)
PY
if [ "$ENABLE_QWEN" = "1" ]; then
  "$PY" - <<'PY'
import importlib
importlib.import_module("transformers")
PY
fi
"$PY" -m morpheus.v2.build_paired_split --data-config "$CFG" --source-split "$SOURCE_SPLIT" --output "$SPLIT" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS"
"$PY" -m morpheus.v2.runtime_preflight --data-config "$CFG" --split-file "$SPLIT" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS"
nvidia-smi

status rna_prepare "streaming raw PanCan matrix conversion in background"
if [ ! -s "$RNA_TABLE" ]; then
  "$PY" -m morpheus.v2.prepare_pancan_rna --source "$RAW_RNA_TABLE" --output "$RNA_TABLE" --work-dir "$RUN/state/rna_work" > "$RUN/logs/rna_prepare.log" 2>&1 &
  RNA_PREP_PID=$!
else
  RNA_PREP_PID=""
fi

status diagnostic_teacher "uncapped all-patch MLP-CLIP inner-fit teacher"
ANCHOR_ARTIFACT="$RUN/artifacts/teachers_inner/mlp_clip_seed42.npz"
mkdir -p "$RUN/artifacts/teachers_inner"
if [ ! -s "$ANCHOR_ARTIFACT" ]; then
  "$PY" -m morpheus.v2.refit_mlp_clip --data-config "$CFG" --split-file "$SPLIT" --fit-population train_only --seed 42 --epochs "$MLP_DIAGNOSTIC_EPOCHS" --batch-size 256 --output "$ANCHOR_ARTIFACT" --checkpoint "$RUN/artifacts/teachers_inner/mlp_clip_seed42.pt" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS"
fi

status metadata "canonical artifact order and frozen RNA/CDR targets"
if [ -n "$RNA_PREP_PID" ]; then wait "$RNA_PREP_PID"; fi
test -s "$RNA_TABLE"
if [ ! -s "$RUN/state/metadata.parquet" ]; then
  "$PY" -m morpheus.v2.build_evaluation_metadata --artifact "$ANCHOR_ARTIFACT" --output "$RUN/state/metadata.parquet"
fi
if [ ! -s "$RUN/discovery_inputs/targets/rna_targets.parquet" ] || [ ! -s "$RUN/discovery_inputs/survival_alignment_audit.json" ]; then
  "$PY" -m morpheus.v2.build_discovery_inputs --rna "$RNA_TABLE" --gmt "$GMT" --training-gmt "$TRAINING_GMT" --metadata "$RUN/state/metadata.parquet" --clinical "$CLINICAL_CDR" --controls-per-target "$CONTROLS_PER_TARGET" --output "$RUN/discovery_inputs"
fi
TARGETS="$RUN/discovery_inputs/targets/rna_targets.parquet"
OUTCOMES="$RUN/discovery_inputs/tcga_cdr_outcomes.parquet"
PANEL="$RUN/discovery_inputs/curated_panel_manifest.json"
"$PY" - "$PANEL" "$RUN/discovery_inputs/targets/rna_target_availability.parquet" "$RUN/state/metadata.parquet" "$MIN_ELIGIBLE_CURATED_TARGETS" "$MIN_MECHANISM_TEST_COVERAGE" <<'PY'
import json, sys
import pandas as pd
from morpheus.v2.curated_panel import FROZEN_MECHANISM_PROGRAMMES

panel_path, availability_path, metadata_path, minimum, coverage = sys.argv[1:]
panel = json.load(open(panel_path, encoding="utf-8"))
if len(panel.get("eligible_targets", [])) < int(minimum):
    raise SystemExit(f"curated panel has only {len(panel.get('eligible_targets', []))} eligible targets; need {minimum}")
availability = pd.read_parquet(availability_path)
metadata = pd.read_parquet(metadata_path)[["patient_id", "split"]]
joined = metadata.merge(availability, on="patient_id", how="left", validate="one_to_one")
test = joined.split.astype(str).eq("test")
missing = set(FROZEN_MECHANISM_PROGRAMMES) - set(joined.columns)
if missing:
    raise SystemExit(f"mechanism target table is missing {sorted(missing)}")
bad = {name: float(joined.loc[test, name].astype(bool).mean()) for name in FROZEN_MECHANISM_PROGRAMMES
       if float(joined.loc[test, name].astype(bool).mean()) < float(coverage)}
if bad:
    raise SystemExit(f"insufficient held-out mechanism target coverage: {bad}")
PY

status baselines "uncapped raw H-Optimus, Ridge, and CCA final development refits"
if [ ! -s "$RUN/artifacts/fair_baselines/baseline_export_manifest.json" ]; then
  "$PY" -m morpheus.v2.export_baselines --data-config "$CFG" --split-file "$SPLIT" --output-dir "$RUN/artifacts/fair_baselines" --fit-development --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS"
fi
if [ -n "$LEGACY_BASELINES_FILE" ]; then
  mapfile -t LEGACY_INPUTS < <(grep -Ev '^\s*(#|$)' "$LEGACY_BASELINES_FILE")
  for path in "${LEGACY_INPUTS[@]}"; do test -f "$path"; done
else
  LEGACY_INPUTS=()
fi

status diagnostic_pretrain "one development-only self-supervised slide encoder"
PRETRAIN_DIR="$RUN/diagnostic/pretrain_full"
PRETRAIN_LOG="$RUN/logs/pretrain_full.log"
if [ -n "$PRETRAIN_CHECKPOINT_OVERRIDE" ]; then
  PRETRAIN_CHECKPOINT="$PRETRAIN_CHECKPOINT_OVERRIDE"
  status diagnostic_pretrain_reuse "using verified prior development-only pretrain checkpoint"
elif [ -s "$PRETRAIN_DIR/slide_pretrain_best.pt" ]; then
  PRETRAIN_CHECKPOINT="$PRETRAIN_DIR/slide_pretrain_best.pt"
  status diagnostic_pretrain_reuse "resuming verified completed development-only pretrain checkpoint"
else
  "$PY" -m morpheus.v2.runner --data-config "$CFG" --split-file "$SPLIT" --output-dir "$PRETRAIN_DIR" --epochs 0 --token-budget "$TOKEN_BUDGET" --seed 42 --objective-profile full --mlp-clip-anchor "$ANCHOR_ARTIFACT" --pretrain-epochs "$PRETRAIN_EPOCHS" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS" 2>&1 | tee "$PRETRAIN_LOG"
  PRETRAIN_CHECKPOINT="$PRETRAIN_DIR/slide_pretrain_best.pt"
fi
test -f "$PRETRAIN_CHECKPOINT"

status diagnostic_train "one-seed objective-profile ablations"
DIAGNOSTIC_INPUTS=("$ANCHOR_ARTIFACT")
for profile in identity_only programme_only full; do
  out="$RUN/diagnostic/$profile"
  artifact="$RUN/artifacts/diagnostic_${profile}_seed42.npz"
  if [ ! -s "$artifact" ]; then
    checkpoint="$(train_with_oom_retry "$out" "$profile" 42 "$DIAGNOSTIC_EPOCHS" "$ANCHOR_ARTIFACT" "$PRETRAIN_CHECKPOINT" 0)"
    export_checkpoint "$checkpoint" "$artifact" "$ANCHOR_ARTIFACT"
  fi
  DIAGNOSTIC_INPUTS+=("$artifact")
done
status diagnostic_evaluate "validation-only promotion gate"
if [ ! -s "$RUN/diagnostic/evaluation/task_rows.csv" ]; then
  "$PY" -m morpheus.v2.v21_evaluation --artifacts "${DIAGNOSTIC_INPUTS[@]}" --targets "$TARGETS" --output "$RUN/diagnostic/evaluation" --evaluation-partition val --bootstrap-repeats "$BOOTSTRAP_REPEATS" --curated-panel-manifest "$PANEL" --matched-control-manifest "$RUN/discovery_inputs/matched_random_control_manifest.json"
fi
"$PY" - "$RUN/diagnostic/candidates.json" <<'PY'
import json, sys
json.dump({"identity_only":"diagnostic_identity_only_seed42", "programme_only":"diagnostic_programme_only_seed42", "full":"diagnostic_full_seed42"}, open(sys.argv[1], "w"), indent=2)
PY
ANCHOR_METHOD="$(basename "$ANCHOR_ARTIFACT" .npz)"
if [ ! -s "$RUN/diagnostic/profile_selection.json" ]; then
  "$PY" -m morpheus.v2.select_v21_profile --rows "$RUN/diagnostic/evaluation/task_rows.csv" --panel "$PANEL" --candidate-map "$RUN/diagnostic/candidates.json" --teacher-method "$ANCHOR_METHOD" --output "$RUN/diagnostic/profile_selection.json"
fi
SELECTED_PROFILE="$("$PY" - "$RUN/diagnostic/profile_selection.json" <<'PY'
import json, sys
result = json.load(open(sys.argv[1]))
if not result["promoted"]: raise SystemExit("diagnostic gate rejected every profile")
print(result["selected_profile"])
PY
)"
SELECTED_EPOCHS="$("$PY" - "$RUN/diagnostic/$SELECTED_PROFILE/selection.json" <<'PY'
import json, sys
selection = json.load(open(sys.argv[1]))
epoch = selection.get("best_pareto_epoch")
if epoch is None:
    epoch = selection.get("best_programme_epoch")
if epoch is None:
    raise SystemExit("selected diagnostic profile did not record a validation-selected epoch")
print(int(epoch) + 1)
PY
)"

status final_pretrain "self-supervised refit on all development patients"
FINAL_PRETRAIN_DIR="$RUN/final/pretrain_development"
if [ ! -s "$FINAL_PRETRAIN_DIR/slide_pretrain_best.pt" ]; then
  "$PY" -m morpheus.v2.runner --data-config "$CFG" --split-file "$SPLIT" --output-dir "$FINAL_PRETRAIN_DIR" --epochs 0 --fit-development --token-budget "$TOKEN_BUDGET" --seed 42 --objective-profile full --pretrain-epochs "$PRETRAIN_EPOCHS" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS" 2>&1 | tee "$RUN/logs/pretrain_development.log"
fi
FINAL_PRETRAIN_CHECKPOINT="$FINAL_PRETRAIN_DIR/slide_pretrain_best.pt"
test -f "$FINAL_PRETRAIN_CHECKPOINT"

status final_train "selected profile=$SELECTED_PROFILE selected_epochs=$SELECTED_EPOCHS development refits"
status final_teachers "seed-matched all-patch MLP-CLIP, SigLIP, and hard-negative CLIP development refits"
mkdir -p "$RUN/artifacts/teachers_final"
TEACHER_BY_SEED_FILE="$RUN/final/teacher_by_seed.json"
POOLED_CACHE="$RUN/artifacts/fair_baselines/uncapped_all_patch_meanstd_cache.npz"
test -f "$POOLED_CACHE"
run_final_teacher() {
  local seed="$1" family="$2" suffix="$3"
  local output="$RUN/artifacts/teachers_final/${suffix}_seed${seed}.npz"
  local checkpoint="$RUN/artifacts/teachers_final/${suffix}_seed${seed}.pt"
  local log="$RUN/logs/${suffix}_seed${seed}.log"
  if [ ! -s "$output" ]; then
    "$PY" -m morpheus.v2.refit_mlp_clip --data-config "$CFG" --split-file "$SPLIT" --fit-population development_train_val --seed "$seed" --epochs "$MLP_FINAL_EPOCHS" --batch-size 256 --loss-family "$family" --pooled-cache "$POOLED_CACHE" --output "$output" --checkpoint "$checkpoint" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS" > "$log" 2>&1
  fi
}
# These dense 256-batch projection baselines use only a small fraction of an
# A100 individually.  Two independent jobs share the GPU to improve measured
# throughput; the 65k-token V2/pretraining jobs remain exclusive because they
# are memory-bound.  Every child is waited on, so any failure stops the DAG.
declare -a BASELINE_PIDS=()
launch_final_teacher() {
  run_final_teacher "$1" "$2" "$3" &
  BASELINE_PIDS+=("$!")
  if [ "${#BASELINE_PIDS[@]}" -ge "$BASELINE_GPU_CONCURRENCY" ]; then
    wait "${BASELINE_PIDS[0]}"
    BASELINE_PIDS=("${BASELINE_PIDS[@]:1}")
  fi
}
for seed in $SEEDS; do
  launch_final_teacher "$seed" clip mlp_clip
  launch_final_teacher "$seed" siglip mlp_siglip
  launch_final_teacher "$seed" hard_negative_clip mlp_clip_hardneg
done
for pid in "${BASELINE_PIDS[@]}"; do wait "$pid"; done
"$PY" - "$TEACHER_BY_SEED_FILE" "$RUN/artifacts/teachers_final" "$SEEDS" <<'PY'
import json, sys
destination, root, seeds = sys.argv[1:]
json.dump({str(seed): f"{root}/mlp_clip_seed{seed}.npz" for seed in seeds.split()}, open(destination, "w"), indent=2, sort_keys=True)
PY
mapfile -t TEACHER_INPUTS < <("$PY" - "$TEACHER_BY_SEED_FILE" <<'PY'
import json, sys
for path in json.load(open(sys.argv[1])).values(): print(path)
PY
)
VARIANT_INPUTS=()
for seed in $SEEDS; do
  VARIANT_INPUTS+=("$RUN/artifacts/teachers_final/mlp_siglip_seed${seed}.npz" "$RUN/artifacts/teachers_final/mlp_clip_hardneg_seed${seed}.npz")
done
# Raw H-Optimus, Ridge, and CCA are deterministic given the fixed feature
# cache and grouped fit protocol; one artifact is the correct seed-invariant
# comparison, not three fabricated replicates.  Every stochastic neural
# baseline is independently refit for seeds 42/43/44 above.
mapfile -t BASELINE_INPUTS < <(dedupe_paths "$RUN/artifacts/fair_baselines/raw_hoptimus_meanstd.npz" "$RUN/artifacts/fair_baselines/ridge_alignment.npz" "$RUN/artifacts/fair_baselines/cca_alignment.npz" "${TEACHER_INPUTS[@]}" "${VARIANT_INPUTS[@]}" "${LEGACY_INPUTS[@]}")
FINAL_INPUTS=("${BASELINE_INPUTS[@]}")
for seed in $SEEDS; do
  teacher="$($PY - "$TEACHER_BY_SEED_FILE" "$seed" <<'PY'
import json, sys
mapping = json.load(open(sys.argv[1])); print(mapping[str(sys.argv[2])])
PY
)"
  test -f "$teacher"
  out="$RUN/final/anchored_seed$seed"
  checkpoint="$(train_with_oom_retry "$out" "$SELECTED_PROFILE" "$seed" "$SELECTED_EPOCHS" "$teacher" "$FINAL_PRETRAIN_CHECKPOINT" 1)"
  artifact="$RUN/artifacts/v21_anchor_${SELECTED_PROFILE}_seed${seed}.npz"
  export_checkpoint "$checkpoint" "$artifact" "$teacher"
  FINAL_INPUTS+=("$artifact")
  if [ "$RUN_NO_ANCHOR" = "1" ]; then
    no_anchor="$RUN/final/no_anchor_seed$seed"
    checkpoint="$(train_with_oom_retry "$no_anchor" "$SELECTED_PROFILE" "$seed" "$SELECTED_EPOCHS" "" "$FINAL_PRETRAIN_CHECKPOINT" 1)"
    artifact="$RUN/artifacts/v21_no_anchor_${SELECTED_PROFILE}_seed${seed}.npz"
    export_checkpoint "$checkpoint" "$artifact" ""
    FINAL_INPUTS+=("$artifact")
  fi
done

status final_evaluate "outer held-out molecular, retrieval, and co-primary survival"
"$PY" - "$RUN/final/teacher_map.json" "$TEACHER_BY_SEED_FILE" "$SELECTED_PROFILE" "$RUN_NO_ANCHOR" <<'PY'
import json, os, sys
destination, source, profile, no_anchor = sys.argv[1:]
teachers = json.load(open(source)); mapping = {}
for seed, path in teachers.items():
    teacher = os.path.splitext(os.path.basename(path))[0]
    mapping[f"v21_anchor_{profile}_seed{seed}"] = teacher
    if no_anchor == "1": mapping[f"v21_no_anchor_{profile}_seed{seed}"] = teacher
json.dump(mapping, open(destination, "w"), indent=2, sort_keys=True)
PY
"$PY" -m morpheus.v2.v21_evaluation --artifacts "${FINAL_INPUTS[@]}" --targets "$TARGETS" --outcomes "$OUTCOMES" --teacher-map "$RUN/final/teacher_map.json" --output "$RUN/final/evaluation" --bootstrap-repeats "$BOOTSTRAP_REPEATS" --curated-panel-manifest "$PANEL" --matched-control-manifest "$RUN/discovery_inputs/matched_random_control_manifest.json"

if [ "$ENABLE_QWEN" = "1" ]; then
  test -n "$CARD_QUALITY_GATE"; test -f "$CARD_QUALITY_GATE"
  status qwen "closed-RAG exploratory hypothesis cards; not a scientific claim endpoint"
  QWEN_METHOD="v21_anchor_${SELECTED_PROFILE}_seed42"
  "$PY" - "$RUN/final/test_patients.txt" "$RUN/state/metadata.parquet" <<'PY'
import pandas as pd, sys
pd.read_parquet(sys.argv[2]).query("split == 'test'").patient_id.astype(str).to_csv(sys.argv[1], index=False, header=False)
PY
  "$PY" -m morpheus.v2.build_hypothesis_cards --predictions "$RUN/final/evaluation/molecular_predictions.parquet" --method "$QWEN_METHOD" --partition-patients "$RUN/final/test_patients.txt" --quality-gate "$CARD_QUALITY_GATE" --output "$RUN/final/qwen_source_cards.jsonl"
  "$PY" -m morpheus.v2.run_qwen_cards --cards "$RUN/final/qwen_source_cards.jsonl" --output "$RUN/final/qwen_cards.jsonl" --quality-gate "$CARD_QUALITY_GATE" --enable-qwen --allow-download --model-id "$QWEN_MODEL_ID"
  "$PY" -m morpheus.v2.evaluate_qwen_cards --qwen-cards "$RUN/final/qwen_cards.jsonl" --source-cards "$RUN/final/qwen_source_cards.jsonl" --targets "$TARGETS" --metadata "$RUN/state/metadata.parquet" --output "$RUN/final/qwen_evaluation.json"
fi

status complete "strict V2.1 outputs ready; Qwen results are exploratory until separate card gates pass"
