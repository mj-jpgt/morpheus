"""Materialize a frozen representation artifact on an exact paired split."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter a frozen artifact without changing its learned vectors")
    parser.add_argument("--source", required=True)
    parser.add_argument("--paired-split", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    split = json.loads(Path(args.paired_split).read_text(encoding="utf-8"))
    ids = np.asarray(sum((list(map(str, split["patient_ids"][part])) for part in ("train", "val", "test")), []), dtype=str)
    labels = np.asarray(sum(([part] * len(split["patient_ids"][part]) for part in ("train", "val", "test")), []), dtype=str)
    with np.load(args.source, allow_pickle=False) as source:
        id_key = "patient_ids" if "patient_ids" in source.files else "patient_id"
        original_ids = source[id_key].astype(str)
        lookup = {patient: row for row, patient in enumerate(original_ids)}
        missing = [patient for patient in ids if patient not in lookup]
        if missing:
            raise ValueError(f"{args.source} lacks {len(missing)} paired patients; first={missing[0]}")
        order = np.asarray([lookup[patient] for patient in ids], dtype=np.int64)
        payload: dict[str, np.ndarray] = {}
        for name in source.files:
            value = source[name]
            if value.ndim and len(value) == len(original_ids):
                payload[name] = value[order]
            else:
                payload[name] = value
        payload["patient_ids"] = ids
        payload["split"] = labels
        payload["filtered_protocol"] = np.asarray("heldout_cancer_11v21")
        payload["filtered_source"] = np.asarray(str(Path(args.source)))
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **payload)
    print(args.output)


if __name__ == "__main__":
    main()
