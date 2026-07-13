"""Prepare TCGA RNA inputs for BulkFormer compatibility checks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from morpheus.src.utils.provenance import write_json


def read_gene_info(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path)


def inspect_h5ad_genes(h5ad_path: str | Path) -> list[str]:
    try:
        import anndata as ad
    except ImportError as exc:
        raise RuntimeError("anndata is required to inspect h5ad gene ordering") from exc
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        return [str(x) for x in adata.var_names]
    finally:
        adata.file.close()


def write_gene_mapping_report(h5ad_path: str | Path, gene_info_path: str | Path, output_path: str | Path) -> Path:
    gene_info = read_gene_info(gene_info_path)
    h5ad_genes = inspect_h5ad_genes(h5ad_path)
    gene_col = next((c for c in gene_info.columns if "gene" in str(c).lower() or "symbol" in str(c).lower()), gene_info.columns[0])
    bulkformer_genes = [str(x) for x in gene_info[gene_col].dropna()]
    overlap = sorted(set(h5ad_genes) & set(bulkformer_genes))
    payload = {
        "h5ad_path": str(h5ad_path),
        "gene_info_path": str(gene_info_path),
        "h5ad_genes": len(h5ad_genes),
        "bulkformer_genes": len(bulkformer_genes),
        "overlap": len(overlap),
        "overlap_fraction_of_bulkformer": len(overlap) / max(1, len(set(bulkformer_genes))),
        "gene_column": gene_col,
    }
    write_json(Path(output_path), payload)
    return Path(output_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--gene-info", required=True)
    parser.add_argument("--output", default="data/processed/rna/bulkformer/gene_mapping_report.json")
    args = parser.parse_args()
    print(write_gene_mapping_report(args.h5ad, args.gene_info, args.output))


if __name__ == "__main__":
    main()
