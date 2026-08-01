"""CLI that freezes curated RNA targets, matched controls, and CDR outcomes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from .canonical_registry import audit_survival_alignment, build_tcga_cdr_outcomes, survival_coverage
from .curated_panel import (FROZEN_MECHANISM_PROGRAMMES, freeze_curated_panel,
                            frozen_mechanism_programme_manifest)
from .discovery_targets import build_matched_random_controls, build_rna_target_bundle, read_gmt, read_rna_table


def _read_table(path: str) -> pd.DataFrame:
    source = Path(path)
    if source.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(source)
    kwargs = {"sep": "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","}
    try:
        return pd.read_csv(source, **kwargs)
    except UnicodeDecodeError:
        # TCGA PanCan clinical mirrors commonly use latin-1 / cp1252 text.
        return pd.read_csv(source, encoding="latin-1", **kwargs)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna", required=True); parser.add_argument("--gmt", required=True)
    parser.add_argument("--metadata", required=True, help="patient_id, cancer, split in canonical artifact order")
    parser.add_argument("--output", required=True); parser.add_argument("--training-gmt", default="")
    parser.add_argument("--clinical", default=""); parser.add_argument("--controls-per-target", type=int, default=20)
    args = parser.parse_args(); output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    metadata, rna, signatures = _read_table(args.metadata), read_rna_table(args.rna), read_gmt(args.gmt)
    training = read_gmt(args.training_gmt) if args.training_gmt else {}
    panel = freeze_curated_panel(signatures, [column for column in rna.columns if column != "patient_id"], training_signatures=training)
    (output / "curated_panel_manifest.json").write_text(json.dumps(panel, indent=2, sort_keys=True), encoding="utf-8")
    eligible = {item["name"]: signatures[item["source_signature"]] for item in panel["targets"] if item["status"] == "eligible"}
    mechanisms = dict(FROZEN_MECHANISM_PROGRAMMES)
    mechanism_manifest = frozen_mechanism_programme_manifest()
    (output / "mechanism_programme_manifest.json").write_text(json.dumps(mechanism_manifest, indent=2, sort_keys=True), encoding="utf-8")
    all_biological_targets = {**eligible, **mechanisms}
    controls, control_manifest = build_matched_random_controls(rna, all_biological_targets, metadata, controls_per_target=args.controls_per_target)
    (output / "matched_random_control_manifest.json").write_text(json.dumps(control_manifest, indent=2, sort_keys=True), encoding="utf-8")
    bundle = build_rna_target_bundle(rna, {**all_biological_targets, **controls}, metadata,
                                     canonical_patient_ids=metadata["patient_id"].astype(str).tolist(), score_mode="competitive_rank")
    bundle.write(output / "targets")
    if args.clinical:
        outcomes = build_tcga_cdr_outcomes(_read_table(args.clinical)); outcomes.to_parquet(output / "tcga_cdr_outcomes.parquet", index=False)
        report = survival_coverage(outcomes, metadata.patient_id.astype(str).tolist(), metadata.split.astype(str).tolist())
        (output / "survival_coverage.json").write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
        audit = audit_survival_alignment(outcomes, metadata.patient_id.astype(str).tolist(), metadata.split.astype(str).tolist())
        (output / "survival_alignment_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
