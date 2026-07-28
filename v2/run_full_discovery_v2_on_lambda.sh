#!/usr/bin/env bash
exec "$(dirname "$0")/run_strict_core_v2_on_lambda.sh" "$@"

# Legacy controller retained below only for historical provenance.  The exec
# above makes it unreachable so no future launch can invoke stale src modules.
set -Eeuo pipefail

ROOT="${ROOT:-/lambda/nfs/geeg/biorag3_persistent_20260711}"
RUN="${RUN:-$ROOT/runs/v2_full_discovery_11v22}"
CFG="${CFG:-$ROOT/morpheus/configs/v2_hoptimus.json}"
SPLIT="${SPLIT:-$ROOT/morpheus/data/processed/splits/tumor_state_heldout_cancer.json}"
EPOCHS="${EPOCHS:-30}"
BATCH_SIZE="${BATCH_SIZE:-256}"
SEEDS="${SEEDS:-42 43 44}"

mkdir -p "$RUN/logs" "$RUN/state"
LOG="$RUN/logs/controller.log"
MONITOR_LOG="$RUN/logs/gpu_monitor.csv"

cd "$ROOT"
PY="$ROOT/.venv-morpheus/bin/python"
export PYTHONPATH="$ROOT:${PYTHONPATH:-}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export HF_HUB_CACHE="${HF_HUB_CACHE:-$ROOT/hf_cache/hub}"
export HF_DATASETS_CACHE="${HF_DATASETS_CACHE:-$ROOT/hf_cache/datasets}"

exec > >(tee -a "$LOG") 2>&1

status() {
  local stage="$1"
  local detail="${2:-}"
  "$PY" - "$RUN/state/status.json" "$stage" "$detail" <<'PY'
import json, sys, time
path, stage, detail = sys.argv[1:4]
data = {"stage": stage, "detail": detail, "time": time.strftime("%Y-%m-%d %H:%M:%S %Z")}
open(path, "w").write(json.dumps(data, indent=2))
PY
  echo "[$(date -Is)] $stage $detail"
}

monitor_gpu() {
  echo "time,gpu_name,mem_used_mb,mem_total_mb,gpu_util_pct,power_w" > "$MONITOR_LOG"
  while true; do
    nvidia-smi --query-gpu=timestamp,name,memory.used,memory.total,utilization.gpu,power.draw --format=csv,noheader,nounits >> "$MONITOR_LOG" || true
    sleep 10
  done
}

monitor_gpu &
MONITOR_PID=$!
trap 'kill "$MONITOR_PID" 2>/dev/null || true' EXIT

status preflight "compile and device check"
"$PY" -m py_compile \
  morpheus/v2/export_frozen_representations.py \
  morpheus/src/eval/evaluate_tumor_state_suite.py \
  morpheus/src/training/train_v2_tumor_state.py \
  morpheus/src/models/v2_tumor_state.py \
  morpheus/v2/summarize_task_suite.py
nvidia-smi

status frozen_export "original MORPHEUS and fair all-patch baselines"
"$PY" morpheus/v2/export_frozen_representations.py

mapfile -t FROZEN_INPUTS < <(find "$ROOT/benchmarks/v2_frozen_exports" -maxdepth 1 -name '*.npz' | sort)
if [ "${#FROZEN_INPUTS[@]}" -lt 6 ]; then
  echo "Expected frozen MORPHEUS/baseline representation files, found ${#FROZEN_INPUTS[@]}" >&2
  exit 2
fi

status frozen_eval "shared discovery-support task suite"
"$PY" -m morpheus.src.eval.evaluate_tumor_state_suite \
  --root "$ROOT" \
  --inputs "${FROZEN_INPUTS[@]}" \
  --output "$RUN/frozen_discovery_eval"
"$PY" morpheus/v2/summarize_task_suite.py \
  --inputs "$RUN/frozen_discovery_eval/task_suite.csv" \
  --output "$RUN/frozen_discovery_eval"

V2_INPUTS=()
for seed in $SEEDS; do
  OUT="$RUN/v2_seed_$seed"
  V2_INPUTS+=("$OUT/representations.npz")
  if [ -f "$OUT/representations.npz" ]; then
    status v2_train_skip "seed=$seed already has representations"
  else
    status v2_train "seed=$seed epochs=$EPOCHS batch_size=$BATCH_SIZE"
    "$PY" -m morpheus.src.training.train_v2_tumor_state \
      --config "$CFG" \
      --split-file "$SPLIT" \
      --output-dir "$OUT" \
      --epochs "$EPOCHS" \
      --batch-size "$BATCH_SIZE" \
      --seed "$seed" 2>&1 | tee "$RUN/logs/v2_seed_${seed}.log"
  fi

  status v2_eval "seed=$seed discovery-support task suite"
  "$PY" -m morpheus.src.eval.evaluate_tumor_state_suite \
    --root "$ROOT" \
    --inputs "$OUT/representations.npz" \
    --output "$OUT/discovery_eval"
  "$PY" morpheus/v2/summarize_task_suite.py \
    --inputs "$OUT/discovery_eval/task_suite.csv" \
    --output "$OUT/discovery_eval"
done

status combined_eval "frozen baselines plus V2 seeds"
"$PY" -m morpheus.src.eval.evaluate_tumor_state_suite \
  --root "$ROOT" \
  --inputs "${FROZEN_INPUTS[@]}" "${V2_INPUTS[@]}" \
  --output "$RUN/combined_discovery_eval"
"$PY" morpheus/v2/summarize_task_suite.py \
  --inputs "$RUN/combined_discovery_eval/task_suite.csv" \
  --output "$RUN/combined_discovery_eval"

status complete "all frozen evaluations, V2 training, and combined evaluations finished"
