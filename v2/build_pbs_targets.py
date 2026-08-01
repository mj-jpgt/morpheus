"""Build immutable Perturbation-Basis Supervision (PBS) targets.

The reference perturbation matrix is read through the one E0-validated loader;
this module does not invent a second interpretation of the Perturb-seq data.
It learns the dictionary only from the external perturbation resource, aligns
gene identity explicitly to the prepared TCGA RNA table, and writes codes for
the canonical paired cohort.  The cohort-dependent residualisation and
legibility operator are deliberately deferred to ``v2.runner`` where they can
be fit on the active development split only.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from morpheus.src.training.train_bio_query_former import load_bio_query_data
from .calibra.e0_basis_transfer import _load_perturbation
from .pbs import ReferenceDictionary


def _digest_strings(values: list[str] | tuple[str, ...] | np.ndarray) -> str:
    return sha256("\n".join(map(str, values)).encode("utf-8")).hexdigest()


def _digest_array(values: np.ndarray) -> str:
    return sha256(np.ascontiguousarray(values).view(np.uint8)).hexdigest()


def _file_sha256(path: str | Path, chunk_bytes: int = 4 * 1024 * 1024) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_bytes), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fit_development_expression_transform(expression: np.ndarray, train_rows: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Center/scale bulk expression only from development patients.

    Perturb-seq stores a control-centred delta, whereas PanCan RNA is a bulk
    abundance measurement.  Projection is meaningful only after both are in
    a declared common coordinate scale.  Test rows never affect the transform.
    """
    value, train = np.asarray(expression, dtype=np.float32), np.asarray(train_rows, dtype=bool)
    if value.ndim != 2 or train.shape != (len(value),) or train.sum() < 2 or not np.isfinite(value).all():
        raise ValueError("PBS expression transform requires finite [patient,gene] values and >=2 development rows")
    if np.abs(value).max() > 1e4:
        raise ValueError("PBS RNA scale is implausibly unlogged (max > 1e4)")
    mean = value[train].mean(axis=0)
    scale = value[train].std(axis=0)
    if np.any(scale < 1e-6):
        raise ValueError("PBS RNA has development-constant genes after exact gene alignment")
    return ((value - mean[None, :]) / scale[None, :]).astype(np.float32), mean.astype(np.float32), scale.astype(np.float32)


def annotate_dictionary_axes(dictionary: ReferenceDictionary, output: Path,
                             annotation_table: str = "") -> dict[str, object]:
    """Write a recoverable, per-axis loading report.

    PBS is only interpretable if a future analysis can inspect the *actual*
    frozen basis used to make the patient targets.  This report never invents
    an essentiality/proliferation score: missing annotations are explicitly
    unavailable rather than silently replaced by an unlabelled axis.
    """
    genes = np.asarray(dictionary.genes, dtype=str)
    basis = np.asarray(dictionary.gene_basis, dtype=np.float64)
    annotations: pd.DataFrame | None = None
    annotation_columns: list[str] = []
    status = "unavailable_missing_gene_annotations"
    if annotation_table:
        annotations = pd.read_parquet(annotation_table) if annotation_table.lower().endswith(".parquet") else pd.read_csv(annotation_table, sep=None, engine="python")
        gene_column = next((column for column in ("gene", "gene_symbol", "symbol") if column in annotations.columns), None)
        if gene_column is None or annotations[gene_column].astype(str).duplicated().any():
            raise ValueError("PBS gene annotation table needs unique gene/gene_symbol/symbol identifiers")
        annotations = annotations.copy()
        annotations["__gene"] = annotations[gene_column].astype(str).str.upper()
        annotation_columns = [column for column in ("proliferation_loading", "essentiality_loading")
                              if column in annotations.columns]
        if not annotation_columns:
            raise ValueError("PBS annotation table must contain proliferation_loading and/or essentiality_loading")
        for column in annotation_columns:
            annotations[column] = pd.to_numeric(annotations[column], errors="coerce")
        annotations = annotations.set_index("__gene")
        status = "annotated" if {"proliferation_loading", "essentiality_loading"}.issubset(annotation_columns) else "partial_annotation"
    records: list[dict[str, object]] = []
    for component in range(dictionary.components):
        loading = basis[:, component]
        order_positive = np.argsort(-loading)[:20]
        order_negative = np.argsort(loading)[:20]
        record: dict[str, object] = {
            "axis": f"PBS_{component:03d}",
            "axis_index": component,
            "top_positive_genes": ";".join(genes[order_positive]),
            "top_negative_genes": ";".join(genes[order_negative]),
            "loading_l1": float(np.abs(loading).sum()),
            "annotation_status": status,
        }
        for column in ("proliferation_loading", "essentiality_loading"):
            if annotations is None or column not in annotation_columns:
                record[column] = np.nan
                continue
            values = annotations.reindex(genes)[column].to_numpy(dtype=np.float64)
            finite = np.isfinite(values)
            record[column] = (float(np.average(values[finite], weights=np.abs(loading[finite])))
                              if finite.any() and np.abs(loading[finite]).sum() > 0 else np.nan)
        records.append(record)
    report = pd.DataFrame(records)
    report.to_csv(output, index=False)
    return {"path": str(output.resolve()), "sha256": _file_sha256(output), "status": status,
            "annotation_columns": annotation_columns, "n_axes": int(dictionary.components)}


def build_pbs_targets(*, data_config: str, split_file: str, rna_table: str,
                      perturbation: str, output: str, n_components: int = 128,
                      fit_population: str = "train", gene_annotations: str = "") -> dict[str, object]:
    """Create patient-ID keyed PBS codes without fitting on held-out patients."""
    if n_components not in {64, 128, 256}:
        raise ValueError("PBS component sensitivity is predeclared at 64, 128, or 256")
    if fit_population not in {"train", "development"}:
        raise ValueError("PBS fit_population must be train or development")
    data = load_bio_query_data(data_config, split_file, wsi_mode="hoptimus_patch")
    reference = _load_perturbation(Path(perturbation))
    import pyarrow.parquet as pq
    schema_names = pq.ParquetFile(rna_table).schema_arrow.names
    if "patient_id" not in schema_names:
        raise ValueError("prepared TCGA RNA parquet requires patient_id")
    gene_to_column = {str(column).upper(): str(column) for column in schema_names if str(column) != "patient_id"}
    if len(gene_to_column) != len(schema_names) - 1:
        raise ValueError("prepared TCGA RNA parquet has duplicate gene symbols after canonicalisation")
    reference_genes = np.asarray(reference.genes, dtype=str)
    keep = np.asarray([gene in gene_to_column for gene in reference_genes], dtype=bool)
    if int(keep.sum()) < max(1000, n_components * 4):
        raise ValueError(f"only {int(keep.sum())} reference/TCGA genes overlap; cannot build a stable PBS dictionary")
    genes = reference_genes[keep].tolist()
    columns = ["patient_id", *[gene_to_column[gene] for gene in genes]]
    table = pd.read_parquet(rna_table, columns=columns)
    table["patient_id"] = table["patient_id"].astype(str)
    if table["patient_id"].duplicated().any():
        raise ValueError("prepared TCGA RNA parquet contains duplicate canonical patient IDs")
    lookup = table.set_index("patient_id")
    missing = [str(patient) for patient in data.patient_ids if str(patient) not in lookup.index]
    if missing:
        raise ValueError(f"prepared RNA misses {len(missing)} canonical paired patients; examples={missing[:5]}")
    expression = lookup.reindex(data.patient_ids)[[gene_to_column[gene] for gene in genes]].to_numpy(dtype=np.float32)
    if not np.isfinite(expression).all():
        raise ValueError("PBS target RNA has non-finite values after canonical patient alignment")
    fit_rows = (np.asarray(data.split).astype(str) == "train" if fit_population == "train"
                else np.asarray(data.split).astype(str) != "test")
    transformed_expression, development_mean, development_scale = fit_development_expression_transform(expression, fit_rows)
    # Reference rows are already E0-validated control-centred deltas, so only
    # the development-fitted gene scale applies; centering them at TCGA bulk
    # abundance would erase the intervention interpretation.
    reference_response = np.asarray(reference.x[:, keep], dtype=np.float32) / development_scale[None, :]
    dictionary = ReferenceDictionary.fit(reference_response, genes, reference.row_ids, n_components=n_components)
    scores = dictionary.encode_expression(transformed_expression, genes).astype(np.float32)
    if not np.isfinite(scores).all() or np.all(np.std(scores, axis=0) < 1e-8):
        raise ValueError("PBS codes are degenerate; inspect RNA scale/gene alignment rather than training")
    output_path = Path(output); output_path.parent.mkdir(parents=True, exist_ok=True)
    target_names = np.asarray([f"PBS_{index:03d}" for index in range(dictionary.components)])
    annotation_manifest = annotate_dictionary_axes(dictionary, output_path.with_suffix(output_path.suffix + ".axis_annotations.csv"),
                                                    gene_annotations)
    manifest = {
        "schema_version": "1.0", "target_kind": "external_perturbation_dictionary_coordinates",
        "n_components": dictionary.components, "dictionary_fit_population": "external_reference_responses_scaled_by_development_TCGA_gene_sd",
        "cohort_fit_population": fit_population + ": bulk expression centering/scaling; runner residualisation",
        "perturbation_path": str(Path(perturbation).resolve()), "rna_table": str(Path(rna_table).resolve()),
        "perturbation_sha256": _file_sha256(perturbation), "rna_table_sha256": _file_sha256(rna_table),
        "data_config_sha256": _file_sha256(data_config), "split_file_sha256": _file_sha256(split_file),
        "reference_atoms": len(reference.row_ids), "reference_gene_count": len(reference.genes),
        "overlap_gene_count": len(genes), "overlap_gene_digest": _digest_strings(genes),
        "canonical_patient_count": len(data.patient_ids), "patient_id_digest": _digest_strings(data.patient_ids),
        "split_digest": _digest_strings(data.split.astype(str)), "target_names": target_names.tolist(),
        "fit_patient_id_digest": _digest_strings(np.asarray(data.patient_ids).astype(str)[fit_rows]),
        "code_std_min": float(scores.std(axis=0).min()), "code_std_max": float(scores.std(axis=0).max()),
        "expression_transform": "(bulk_expression-development_train_gene_mean)/development_train_gene_sd; reference_delta/development_train_gene_sd",
        "development_mean_digest": _digest_array(development_mean), "development_scale_digest": _digest_array(development_scale),
        "scores_sha256": _digest_array(scores), "singular_values_sha256": _digest_array(dictionary.singular_values.astype(np.float32)),
        "gene_basis_sha256": _digest_array(dictionary.gene_basis.astype(np.float32)),
        "atom_coordinates_sha256": _digest_array(dictionary.atom_coordinates.astype(np.float32)),
        "atom_id_digest": _digest_strings(dictionary.atom_ids),
        "axis_annotations": annotation_manifest,
    }
    with tempfile.NamedTemporaryFile(dir=output_path.parent, suffix=".npz", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        np.savez_compressed(temporary, patient_ids=np.asarray(data.patient_ids), split=np.asarray(data.split),
                            scores=scores, target_names=target_names, target_groups=np.asarray(["PBS"] * len(target_names)),
                            genes=np.asarray(genes), singular_values=dictionary.singular_values.astype(np.float32),
                            gene_basis=dictionary.gene_basis.astype(np.float32),
                            atom_coordinates=dictionary.atom_coordinates.astype(np.float32),
                            atom_ids=np.asarray(dictionary.atom_ids),
                            manifest_json=np.asarray(json.dumps(manifest, sort_keys=True)))
        os.replace(temporary, output_path)
    finally:
        temporary.unlink(missing_ok=True)
    output_path.with_suffix(output_path.suffix + ".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-config", default="morpheus/configs/v1.json")
    parser.add_argument("--split-file", required=True)
    parser.add_argument("--rna-table", required=True, help="prepared patient-by-gene RNA parquet")
    parser.add_argument("--perturbation", required=True, help="E0-validated K562 Perturb-seq h5ad")
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-components", type=int, default=128, choices=(64, 128, 256))
    parser.add_argument("--fit-population", default="train", choices=("train", "development"),
                        help="fit train-only for inner diagnostics; rebuild with development for final refit")
    parser.add_argument("--gene-annotations", default="", help="optional gene-level proliferation/essentiality annotation table; missing data is reported unavailable")
    print(json.dumps(build_pbs_targets(**vars(parser.parse_args())), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
