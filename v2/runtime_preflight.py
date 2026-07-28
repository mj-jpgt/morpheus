"""Read-only runtime gate for the strict paired WSI--RNA protocol."""
from __future__ import annotations

import argparse
import json

from morpheus.src.training.train_bio_query_former import load_bio_query_data
from .preflight import validate_runtime_split


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the loaded strict V2 paired cohort")
    parser.add_argument("--data-config", required=True)
    parser.add_argument("--split-file", required=True)
    parser.add_argument("--expected-development-cancers", type=int, default=11)
    parser.add_argument("--expected-heldout-cancers", type=int, default=22)
    args = parser.parse_args()
    data = load_bio_query_data(args.data_config, args.split_file, wsi_mode="hoptimus_patch")
    print(json.dumps(validate_runtime_split(args.split_file, data.patient_ids, data.cancers, data.split,
                                             args.expected_development_cancers, args.expected_heldout_cancers), indent=2))


if __name__ == "__main__":
    main()
