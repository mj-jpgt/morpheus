"""Materialize canonical evaluation metadata from a validated representation artifact."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .contracts import validate_artifact
from .preflight import sha256_json
from .provenance import write_source_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True, help="split-matched anchor/baseline NPZ")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    validate_artifact(args.artifact)
    with np.load(args.artifact, allow_pickle=False) as raw:
        patient_ids = raw["patient_ids"].astype(str)
        cancers = raw["cancers"].astype(str)
        split = raw["split"].astype(str)
    if len(set(patient_ids)) != len(patient_ids):
        raise ValueError("artifact contains duplicate patient IDs")
    rows = [
        {"patient_id": patient, "cancer": cancer, "split": label}
        for patient, cancer, label in zip(patient_ids, cancers, split)
    ]
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.suffix.lower() in {".parquet", ".pq"}:
        import pandas as pd
        pd.DataFrame(rows).to_parquet(destination, index=False)
    else:
        destination.write_text("patient_id,cancer,split\n" + "".join(f"{x['patient_id']},{x['cancer']},{x['split']}\n" for x in rows), encoding="utf-8")
    manifest = {"artifact": str(Path(args.artifact)), "cohort_digest": sha256_json(rows), "n_patients": len(rows)}
    destination.with_suffix(destination.suffix + ".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    write_source_manifest(destination.with_suffix(destination.suffix + ".source.json"), configuration=manifest)


if __name__ == "__main__":
    main()
