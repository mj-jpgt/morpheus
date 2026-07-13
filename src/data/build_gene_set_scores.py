"""Build Hallmark/Reactome-style gene-set activity targets from RNA tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.src.utils.ids import normalize_patient_id
from morpheus.src.utils.provenance import write_json


def parse_gmt(path: str | Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                gene_sets[parts[0]] = sorted({g.upper() for g in parts[2:] if g})
    return gene_sets


def _gene_symbol(column: object) -> str:
    text = str(column)
    if "|" in text:
        text = text.split("|", 1)[0]
    return text.strip().upper()


def load_rna_matrix(path: str | Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    patient_col = next((c for c in df.columns if "patient" in str(c).lower()), df.columns[0])
    df = df.copy()
    df["patient_id"] = df[patient_col].map(normalize_patient_id)
    numeric = df.select_dtypes(include=["number"]).columns.tolist()
    rename = {col: _gene_symbol(col) for col in numeric}
    expr = df[["patient_id", *numeric]].rename(columns=rename).dropna(subset=["patient_id"])
    # Collapse duplicate symbols after stripping Entrez suffixes.
    expr = expr.groupby("patient_id", as_index=True).mean()
    expr = expr.T.groupby(level=0).mean().T
    expr = np.log2(expr.clip(lower=0.0) + 1.0)
    return expr.astype(np.float32)


def _mean_zscore_scores(expr: pd.DataFrame, gene_sets: dict[str, list[str]], min_genes: int) -> tuple[pd.DataFrame, dict]:
    z = (expr - expr.mean(axis=0)) / expr.std(axis=0).replace(0, np.nan)
    z = z.fillna(0.0)
    rows: dict[str, pd.Series] = {}
    coverage = {}
    genes = set(z.columns)
    for name, members in gene_sets.items():
        overlap = [g for g in members if g in genes]
        coverage[name] = {"n_genes": len(members), "n_overlap": len(overlap)}
        if len(overlap) >= min_genes:
            rows[name] = z[overlap].mean(axis=1)
    return pd.DataFrame(rows, index=expr.index), coverage


def _ssgsea_scores(expr: pd.DataFrame, gene_sets: dict[str, list[str]], min_genes: int, threads: int) -> tuple[pd.DataFrame, dict]:
    import gseapy as gp

    usable = {}
    genes = set(expr.columns)
    coverage = {}
    for name, members in gene_sets.items():
        overlap = [g for g in members if g in genes]
        coverage[name] = {"n_genes": len(members), "n_overlap": len(overlap)}
        if len(overlap) >= min_genes:
            usable[name] = overlap
    result = gp.ssgsea(
        data=expr.T,
        gene_sets=usable,
        outdir=None,
        sample_norm_method="rank",
        min_size=min_genes,
        max_size=500,
        permutation_num=0,
        threads=threads,
        no_plot=True,
        seed=42,
        verbose=False,
    )
    res = result.res2d.copy()
    if {"Name", "Term", "NES"}.issubset(res.columns):
        scores = res.pivot(index="Name", columns="Term", values="NES")
    else:
        scores = pd.DataFrame(result.resultsOnSamples).T
    scores.index = scores.index.map(normalize_patient_id)
    return scores.astype(np.float32), coverage


def build_gene_set_scores(
    rna_path: str | Path,
    gmt_path: str | Path,
    output_path: str | Path,
    manifest_path: str | Path,
    method: str = "ssgsea",
    min_genes: int = 10,
    threads: int = 4,
    max_samples: int | None = None,
) -> Path:
    expr = load_rna_matrix(rna_path)
    if max_samples is not None:
        expr = expr.iloc[:max_samples]
    gene_sets = parse_gmt(gmt_path)
    if method == "ssgsea":
        scores, coverage = _ssgsea_scores(expr, gene_sets, min_genes, threads)
    elif method == "mean_zscore":
        scores, coverage = _mean_zscore_scores(expr, gene_sets, min_genes)
    else:
        raise ValueError(f"Unknown gene-set scoring method: {method}")
    scores = scores.reset_index().rename(columns={"index": "patient_id", "Name": "patient_id"})
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    scores.to_parquet(output, index=False)
    write_json(
        Path(manifest_path),
        {
            "rna_path": str(rna_path),
            "gmt_path": str(gmt_path),
            "output_path": str(output),
            "method": method,
            "transform": "log2(x+1) before scoring",
            "n_patients": int(len(scores)),
            "n_gene_sets": int(scores.shape[1] - 1),
            "min_genes": min_genes,
            "coverage": coverage,
        },
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", default="meta-intersurv/data/omics/processed/rna_features.parquet")
    parser.add_argument("--gmt", default="meta-intersurv/data/msigdb/h.all.v2024.1.Hs.symbols.gmt")
    parser.add_argument("--output", default="data/processed/genesets/msigdb_hallmark_scores.parquet")
    parser.add_argument("--manifest", default="data/processed/genesets/msigdb_hallmark_scores.manifest.json")
    parser.add_argument("--method", choices=["ssgsea", "mean_zscore"], default="ssgsea")
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-samples", type=int, default=None)
    args = parser.parse_args()
    print(build_gene_set_scores(args.rna, args.gmt, args.output, args.manifest, args.method, args.min_genes, args.threads, args.max_samples))


if __name__ == "__main__":
    main()
