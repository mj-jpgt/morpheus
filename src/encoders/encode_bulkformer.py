"""Guarded BulkFormer setup entrypoint.

BulkFormer is a Zenodo-first dependency for the official v1 RNA encoder. This
module keeps the old entrypoint name, but it now validates local Zenodo
artifacts instead of trying to download from Hugging Face.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from morpheus.src.utils.config import load_config
from morpheus.src.data.bulkformer_zenodo import verify_download_manifest


def ensure_bulkformer_configured(config_path: str = "morpheus/configs/v1.json") -> Path:
    cfg = load_config(config_path)
    rna = dict(cfg.section("foundation_models").get("rna_primary", {}))
    record_id = str(rna.get("zenodo_record_id", "")).strip()
    if not record_id:
        raise RuntimeError(
            "BulkFormer is required for official v1, but foundation_models.rna_primary.zenodo_record_id is blank. "
            "Set the BulkFormer Zenodo record ID before running this step."
        )
    manifest_path = cfg.path("bulkformer_cache_dir") / "bulkformer_zenodo_manifest.json"
    if not manifest_path.exists():
        raise RuntimeError(
            f"BulkFormer Zenodo manifest not found at {manifest_path}. "
            "Run scripts/download_bulkformer_zenodo.py first."
        )
    result = verify_download_manifest(manifest_path)
    if result["missing_required"]:
        missing = ", ".join(result["missing_required"])
        raise RuntimeError(f"BulkFormer artifacts are incomplete: missing {missing}")
    if result["checksum_failures"]:
        failed = ", ".join(result["checksum_failures"])
        raise RuntimeError(f"BulkFormer artifact checksum failures: {failed}")
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    args = parser.parse_args()
    manifest_path = ensure_bulkformer_configured(args.config)
    print(manifest_path)


if __name__ == "__main__":
    main()
