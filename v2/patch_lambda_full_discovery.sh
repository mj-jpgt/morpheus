#!/usr/bin/env bash
set -Eeuo pipefail

cd /lambda/nfs/geeg/biorag3_persistent_20260711
chmod +x morpheus/v2/run_full_discovery_v2_on_lambda.sh

perl -0pi -e 's/batch_ids\(np\.arange\(len\(data\.patient_ids\)\),32\)/batch_ids(np.arange(len(data.patient_ids)),256)/g; s/batch_ids\(np\.arange\(len\(data\.patient_ids\)\),64\)/batch_ids(np.arange(len(data.patient_ids)),256)/g' morpheus/v2/export_frozen_representations.py
perl -0pi -e "s/method=Path\(npz_path\)\.stem/method=Path(npz_path).parent.name if Path(npz_path).stem == 'representations' else Path(npz_path).stem/" morpheus/src/eval/evaluate_tumor_state_suite.py

python -m py_compile \
  morpheus/v2/export_frozen_representations.py \
  morpheus/src/eval/evaluate_tumor_state_suite.py \
  morpheus/v2/summarize_task_suite.py \
  morpheus/src/training/train_v2_tumor_state.py \
  morpheus/src/models/v2_tumor_state.py
