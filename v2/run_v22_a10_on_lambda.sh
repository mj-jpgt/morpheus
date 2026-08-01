#!/usr/bin/env bash
# Launch the strict V2.2 recovery protocol on a single 23 GB A10.
#
# This wrapper sets hardware-safe defaults only.  The underlying recovery
# script performs the split, target, baseline, training, and evaluation work;
# V2.2 source/manifests make each scientific artifact auditable.
set -Eeuo pipefail

ROOT="${ROOT:-/lambda/nfs/geeg/biorag3_persistent_20260711}"
CODE="${CODE:-$ROOT/code_snapshots/v22_implementation_20260724}"
RUN="${RUN:-$ROOT/runs/v22_a10_11v21}"
PY="${PY:-$ROOT/.venv-morpheus/bin/python}"

export ROOT CODE RUN PY
# 16k is the first measured A10 candidate.  The recovery runner halves this
# to 8k after OOM; it never discards patches.
export TOKEN_BUDGET="${TOKEN_BUDGET:-16384}"
export MIN_TOKEN_BUDGET="${MIN_TOKEN_BUDGET:-8192}"
export BASELINE_GPU_CONCURRENCY="${BASELINE_GPU_CONCURRENCY:-1}"
export ENABLE_QWEN="${ENABLE_QWEN:-0}"
export SEEDS="${SEEDS:-42 43 44}"

test -x "$PY"
test -f "$CODE/morpheus/v2/run_v21_recovery_on_lambda.sh"
exec bash "$CODE/morpheus/v2/run_v21_recovery_on_lambda.sh"
