"""Build a persistent, exact paired-cohort split from a source split manifest."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from morpheus.src.training.train_bio_query_former import load_bio_query_data
from .preflight import sha256_json, validate_runtime_split


def main() -> None:
    parser = argparse.ArgumentParser(description="Materialize the exact paired cohort for strict V2")
    parser.add_argument("--data-config", required=True)
    parser.add_argument("--source-split", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-development-cancers", type=int, default=11)
    parser.add_argument("--expected-heldout-cancers", type=int, default=22)
    args = parser.parse_args()
    source = json.loads(Path(args.source_split).read_text(encoding="utf-8"))
    data = load_bio_query_data(args.data_config, args.source_split, wsi_mode="hoptimus_patch")
    actual = set(map(str, data.patient_ids))
    patient_ids = {name: sorted(actual & set(map(str, source["patient_ids"][name])))
                   for name in ("train", "val", "test")}
    output = dict(source)
    output["patient_ids"] = patient_ids
    output["paired_from_source_split"] = str(Path(args.source_split))
    output["paired_source_digest"] = sha256_json(source.get("patient_ids", {}))
    output["paired_cohort_count"] = len(actual)
    train_cancers = {str(cancer) for cancer, label in zip(data.cancers, data.split) if str(label) == "train"}
    test_cancers = {str(cancer) for cancer, label in zip(data.cancers, data.split) if str(label) == "test"}
    output["actual_train_test_cancer_counts"] = [len(train_cancers), len(test_cancers)]
    output["paired_train_cancers"] = sorted(train_cancers)
    output["paired_test_cancers"] = sorted(test_cancers)
    destination = Path(args.output); destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(output, indent=2), encoding="utf-8")
    manifest = validate_runtime_split(destination, data.patient_ids, data.cancers, data.split,
                                      args.expected_development_cancers, args.expected_heldout_cancers)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
