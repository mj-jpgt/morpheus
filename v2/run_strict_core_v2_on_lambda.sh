#!/usr/bin/env bash
# Strict core: uncapped H-Optimus + BulkFormer only.  No remote source edits.
set -Eeuo pipefail

ROOT="${ROOT:-/lambda/nfs/geeg/biorag3_persistent_20260711}"
CODE="${CODE:-$ROOT/code}"
RUN="${RUN:-$ROOT/runs/v2_strict_core_11v21}"
CFG="${CFG:-$CODE/morpheus/configs/v1.json}"
SOURCE_SPLIT="${SOURCE_SPLIT:-$CODE/morpheus/data/processed/splits/tumor_state_heldout_cancer.json}"
SPLIT="${SPLIT:-$RUN/state/paired_tumor_state_heldout_cancer.json}"
ANCHOR_ARTIFACT="${ANCHOR_ARTIFACT:?Set ANCHOR_ARTIFACT to the split-matched all-patch MLP-CLIP artifact}"
BASELINE_INPUTS="${BASELINE_INPUTS:-}"
BASELINE_INPUT_FILE="${BASELINE_INPUT_FILE:-}"
SEEDS="${SEEDS:-42 43 44}"
EPOCHS="${EPOCHS:-40}"
TOKEN_BUDGET="${TOKEN_BUDGET:-32768}"
EXPECTED_DEVELOPMENT_CANCERS="${EXPECTED_DEVELOPMENT_CANCERS:-11}"
EXPECTED_HELDOUT_CANCERS="${EXPECTED_HELDOUT_CANCERS:-21}"
PY="${PY:-$ROOT/.venv-morpheus/bin/python}"

mkdir -p "$RUN/logs" "$RUN/state"
export PYTHONPATH="$CODE:${PYTHONPATH:-}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export HF_HUB_CACHE="${HF_HUB_CACHE:-$ROOT/hf_cache/hub}"
export HF_DATASETS_CACHE="${HF_DATASETS_CACHE:-$ROOT/hf_cache/datasets}"
cd "$CODE"
exec > >(tee -a "$RUN/logs/controller.log") 2>&1

status() {
  "$PY" - "$RUN/state/status.json" "$1" "${2:-}" <<'PY'
import json, sys, time
open(sys.argv[1], "w").write(json.dumps({"stage":sys.argv[2],"detail":sys.argv[3],"time":time.strftime("%FT%T%z")}, indent=2))
PY
  echo "[$(date -Is)] $1 ${2:-}"
}
monitor() {
  echo "time,gpu_name,mem_used_mb,mem_total_mb,gpu_util_pct,power_w" > "$RUN/logs/gpu_monitor.csv"
  while true; do nvidia-smi --query-gpu=timestamp,name,memory.used,memory.total,utilization.gpu,power.draw --format=csv,noheader,nounits >> "$RUN/logs/gpu_monitor.csv" || true; sleep 10; done
}
monitor & MONITOR_PID=$!; trap 'kill "$MONITOR_PID" 2>/dev/null || true' EXIT

status preflight "strict source, split, and device checks"
test -f "$ANCHOR_ARTIFACT"; test -f "$SOURCE_SPLIT"
PYTHONPYCACHEPREFIX=/tmp/morpheus_pyc "$PY" -m py_compile morpheus/v2/preflight.py morpheus/v2/build_paired_split.py morpheus/v2/runtime_preflight.py morpheus/v2/runner.py morpheus/v2/export.py morpheus/v2/comprehensive_evaluation.py
PYTHONPYCACHEPREFIX=/tmp/morpheus_pyc "$PY" -m morpheus.v2.build_paired_split --data-config "$CFG" --source-split "$SOURCE_SPLIT" --output "$SPLIT" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS"
PYTHONPYCACHEPREFIX=/tmp/morpheus_pyc "$PY" -m morpheus.v2.runtime_preflight --data-config "$CFG" --split-file "$SPLIT" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS"
nvidia-smi

INPUTS=()
if [ -z "$BASELINE_INPUT_FILE" ] && [ -z "$BASELINE_INPUTS" ]; then
  echo "Set BASELINE_INPUT_FILE or BASELINE_INPUTS to split-matched baseline/V1 artifacts" >&2
  exit 2
fi
if [ -n "$BASELINE_INPUT_FILE" ]; then
  test -f "$BASELINE_INPUT_FILE"
  mapfile -t BASELINES < "$BASELINE_INPUT_FILE"
else
  read -r -a BASELINES <<< "$BASELINE_INPUTS" || true
fi
for item in "${BASELINES[@]}"; do test -f "$item"; INPUTS+=("$item"); done

for variant in anchored no_anchor; do
  for seed in $SEEDS; do
    OUT="$RUN/${variant}_seed_${seed}"
    REP="$OUT/representations.npz"
    INPUTS+=("$REP")
    if [ -f "$REP" ]; then status skip "${variant} seed=$seed already exported"; continue; fi
    mkdir -p "$OUT"
    status train "${variant} seed=$seed"
    TRAIN=("$PY" -m morpheus.v2.runner --data-config "$CFG" --split-file "$SPLIT" --output-dir "$OUT" --epochs "$EPOCHS" --token-budget "$TOKEN_BUDGET" --seed "$seed" --expected-development-cancers "$EXPECTED_DEVELOPMENT_CANCERS" --expected-heldout-cancers "$EXPECTED_HELDOUT_CANCERS")
    if [ "$variant" = anchored ]; then TRAIN+=(--mlp-clip-anchor "$ANCHOR_ARTIFACT"); fi
    CHECKPOINT="$("${TRAIN[@]}" | tee "$RUN/logs/${variant}_seed_${seed}.log" | tail -n 1)"
    test -f "$CHECKPOINT"
    status export "${variant} seed=$seed checkpoint=$(basename "$CHECKPOINT")"
    EXPORT=("$PY" -m morpheus.v2.export --data-config "$CFG" --split-file "$SPLIT" --checkpoint "$CHECKPOINT" --output "$REP" --token-budget "$TOKEN_BUDGET")
    if [ "$variant" = anchored ]; then EXPORT+=(--mlp-clip-anchor "$ANCHOR_ARTIFACT"); fi
    "${EXPORT[@]}" | tee "$RUN/logs/export_${variant}_seed_${seed}.log"
  done
done

status evaluate "shared strict-core task suite"
EVAL=("$PY" -m morpheus.v2.comprehensive_evaluation --root "$ROOT" --inputs "${INPUTS[@]}" --output "$RUN/task_suite")
if [ -n "${PLIP_PROMPTS:-}" ]; then EVAL+=(--plip-prompts "$PLIP_PROMPTS"); fi
"${EVAL[@]}"
status complete "strict core completed"
