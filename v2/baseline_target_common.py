"""Shared cohort/transform loading for the T1.1 must-beat baseline target blocks.

Every baseline in the T1.1 table has to be bound to *the same* 6,427-patient
maximal split, the same gene set, the same expression transform and the same
train/development leak discipline as the PBS targets it is being compared with.
If each builder re-derived those independently, a difference in the baseline
table could come from a cohort difference rather than from the representation —
and that difference would be invisible in the numbers.

So the cohort is not re-derived at all: it is *inherited* from the frozen PBS
target file, whose manifest already records the split digest, the exact gene
order, the log transform and the fit population. A baseline built through this
module is bound by construction to the block it is competing with.
"""
from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

__all__ = ["load_reference_targets", "development_expression", "write_target_block"]


def _digest_strings(values) -> str:
    return sha256("\n".join(map(str, values)).encode("utf-8")).hexdigest()


def load_reference_targets(path: str | Path) -> dict:
    """Read cohort, split, gene order and dictionary geometry from a frozen PBS npz."""
    raw = np.load(path, allow_pickle=True)
    required = {"patient_ids", "split", "cancers", "genes", "gene_basis", "singular_values"}
    missing = required - set(raw.files)
    if missing:
        raise ValueError(f"reference target file is missing {sorted(missing)}; it cannot anchor a baseline")
    manifest = json.loads(str(raw["manifest_json"])) if "manifest_json" in raw.files else {}
    genes = np.asarray(raw["genes"]).astype(str)
    return {
        "patient_ids": np.asarray([str(p) for p in raw["patient_ids"]]),
        "split": np.asarray([str(s) for s in raw["split"]]),
        "cancers": np.asarray([str(c) for c in raw["cancers"]]),
        "genes": genes, "gene_digest": _digest_strings(genes),
        "gene_basis": np.asarray(raw["gene_basis"], dtype=np.float64),
        "singular_values": np.asarray(raw["singular_values"], dtype=np.float64),
        "rna_log_transform": str(manifest.get("rna_log_transform", "none")),
        "fit_population": str(manifest.get("cohort_fit_population", "development")),
        "manifest_digest": sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest(),
        "source": str(Path(path).resolve()),
    }


def development_expression(rna_table: str | Path, reference: dict) -> tuple[np.ndarray, dict]:
    """Return the patient-by-gene matrix on the SAME scale the PBS targets used.

    Development rows are ``split != "test"`` — the same rule
    ``build_pbs_targets`` applies for ``fit_population="development"``. The
    centring and scaling are estimated on those rows only, so no held-out cancer
    can influence the transform any baseline is built through.
    """
    from .build_pbs_targets import fit_development_expression_transform

    genes = reference["genes"]
    import pyarrow.parquet as pq
    schema = pq.ParquetFile(str(rna_table)).schema_arrow.names
    lookup = {str(c).upper(): str(c) for c in schema if str(c) != "patient_id"}
    missing = [g for g in genes if g not in lookup]
    if missing:
        raise ValueError(f"{len(missing)} reference genes are absent from the RNA table (e.g. {missing[:5]})")
    table = pd.read_parquet(rna_table, columns=["patient_id", *[lookup[g] for g in genes]])
    table["patient_id"] = table["patient_id"].astype(str)
    indexed = table.set_index("patient_id")
    absent = [p for p in reference["patient_ids"] if p not in indexed.index]
    if absent:
        raise ValueError(f"{len(absent)} cohort patients have no RNA row (e.g. {absent[:5]})")
    values = indexed.reindex(reference["patient_ids"])[[lookup[g] for g in genes]].to_numpy(dtype=np.float32)
    transform = reference["rna_log_transform"]
    if transform == "signed_log1p":
        values = np.sign(values) * np.log2(np.abs(values) + 1.0)
    elif transform == "clip_log1p":
        values = np.log2(np.maximum(values, 0.0) + 1.0)
    elif transform != "none":
        raise ValueError(f"unknown rna_log_transform {transform!r}")
    if not np.isfinite(values).all():
        raise ValueError("RNA matrix is non-finite after transform; the reference gene set should already exclude those")
    fit_rows = reference["split"] != "test"
    scaled, mean, scale = fit_development_expression_transform(values, fit_rows)
    return np.asarray(scaled, dtype=np.float64), {
        "rna_log_transform": transform, "fit_population": "development_rows_split_not_test",
        "n_fit_rows": int(fit_rows.sum()), "n_patients": int(len(values)),
        "mean_digest": sha256(np.ascontiguousarray(mean).view(np.uint8)).hexdigest(),
        "scale_digest": sha256(np.ascontiguousarray(scale).view(np.uint8)).hexdigest(),
        "rna_table": str(Path(rna_table).resolve()),
    }


def write_target_block(output: str | Path, reference: dict, scores: np.ndarray,
                       target_names: np.ndarray, group: str, manifest: dict) -> dict:
    """Write a CALIBRA-readable target npz bound to the inherited cohort."""
    scores = np.asarray(scores, dtype=np.float32)
    if scores.shape != (len(reference["patient_ids"]), len(target_names)):
        raise ValueError("target block must be [patient, target] on the inherited cohort")
    standard_deviation = scores.std(axis=0)
    if not np.isfinite(scores).all() or float(standard_deviation.min()) < 1e-8:
        raise ValueError("target block is non-finite or has a constant column; CALIBRA would refuse it")
    manifest = dict(manifest)
    manifest.update({"schema_version": "1.0", "target_group": group,
                     "cohort_source": reference["source"],
                     "n_patients": int(len(reference["patient_ids"])),
                     "patient_id_digest": _digest_strings(reference["patient_ids"]),
                     "split_digest": _digest_strings(reference["split"]),
                     "code_std_min": float(standard_deviation.min()),
                     "code_std_max": float(standard_deviation.max())})
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".npz", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        np.savez_compressed(temporary, patient_ids=reference["patient_ids"], split=reference["split"],
                            cancers=reference["cancers"], scores=scores,
                            target_names=np.asarray(target_names),
                            target_groups=np.asarray([group] * len(target_names)),
                            manifest_json=np.asarray(json.dumps(manifest, sort_keys=True)))
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    path.with_suffix(path.suffix + ".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
