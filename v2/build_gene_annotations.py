"""Per-gene proliferation and essentiality annotations for the PBS axis report (D2.3).

Why this exists. `claim_guards.proliferation_deflation` is the caveat most likely to be
true and cheapest to test: the responsive arm is selected on *having a detectable effect*,
which enriches for essential / ribosome / cell-cycle genes, so an alignment that is really
"proliferation matching proliferation" would look identical to a real transfer result. The
D2 axis report answers it for free -- but only if every axis carries a real annotation.

Neither score is invented. Both come from public measurements, and a gene absent from a
source is recorded as missing rather than imputed:

  essentiality_loading -- mean DepMap CRISPR (Chronos) gene effect across cell lines,
      SIGN-FLIPPED so that larger means MORE essential. Chronos is negative for genes whose
      knockout impairs fitness, so -mean(effect) is the natural "essentiality" direction.

  proliferation_loading -- membership in the union of the MSigDB Hallmark proliferation
      programmes (E2F_TARGETS, G2M_CHECKPOINT, MYC_TARGETS_V1, MYC_TARGETS_V2,
      MITOTIC_SPINDLE). Binary per gene, so the axis-level aggregate is the readable
      quantity "what fraction of this axis's top-loading genes are proliferation genes".

Deliberately NOT done: no imputation, no rank-normalisation, no combining the two into a
single score. The report must let a reader see that an axis is proliferation-driven.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROLIFERATION_SETS = (
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_MITOTIC_SPINDLE",
)


def _symbol(value: str) -> str:
    """'A1BG (1)' -> 'A1BG'. DepMap columns carry the Entrez id in parentheses."""
    return str(value).split("(")[0].strip().upper()


def read_hallmark_proliferation(gmt_path: str | Path) -> tuple[set[str], dict[str, int]]:
    members: set[str] = set()
    per_set: dict[str, int] = {}
    with Path(gmt_path).open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            name = fields[0].strip().upper()
            if name in PROLIFERATION_SETS:
                genes = {g.strip().upper() for g in fields[2:] if g.strip()}
                per_set[name] = len(genes)
                members |= genes
    missing = [name for name in PROLIFERATION_SETS if name not in per_set]
    if missing:
        raise ValueError(f"GMT is missing declared proliferation sets: {missing}")
    return members, per_set


def read_depmap_essentiality(effect_path: str | Path) -> pd.Series:
    """Mean Chronos effect per gene, sign-flipped so larger == more essential."""
    frame = pd.read_csv(effect_path, index_col=0)
    values = frame.apply(pd.to_numeric, errors="coerce")
    mean_effect = values.mean(axis=0, skipna=True)
    mean_effect.index = [_symbol(column) for column in mean_effect.index]
    # A duplicated symbol (two Entrez ids) is averaged rather than dropped, and the
    # collision count is reported so the choice is visible.
    essentiality = (-mean_effect).groupby(level=0).mean()
    return essentiality[np.isfinite(essentiality.to_numpy())]


def build(*, effect: str, gmt: str, output: str, restrict_to: str = "") -> dict:
    proliferation_members, per_set = read_hallmark_proliferation(gmt)
    essentiality = read_depmap_essentiality(effect)

    genes = sorted(set(essentiality.index))
    if restrict_to:
        wanted = {str(g).strip().upper() for g in pd.read_csv(restrict_to)["gene"]}
        genes = sorted(set(genes) & wanted)
        if not genes:
            raise ValueError("no overlap between DepMap genes and the restriction list")

    table = pd.DataFrame({
        "gene": genes,
        "essentiality_loading": [float(essentiality[g]) for g in genes],
        "proliferation_loading": [1.0 if g in proliferation_members else 0.0 for g in genes],
    })
    if table["gene"].duplicated().any():
        raise ValueError("annotation table must have unique gene symbols")
    if not np.isfinite(table[["essentiality_loading", "proliferation_loading"]].to_numpy()).all():
        raise ValueError("annotation table contains non-finite values; refusing to emit")

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(output, index=False)
    manifest = {
        "n_genes": int(len(table)),
        "n_proliferation_genes": int(table["proliferation_loading"].sum()),
        "proliferation_sets": {name: int(size) for name, size in per_set.items()},
        "proliferation_union_size": int(len(proliferation_members)),
        "essentiality_source": str(effect),
        "essentiality_definition": "-mean(DepMap CRISPR Chronos gene effect) across cell lines; larger = more essential",
        "proliferation_definition": "binary membership in the union of the five Hallmark proliferation programmes",
        "gmt_source": str(gmt),
        "essentiality_min": float(table["essentiality_loading"].min()),
        "essentiality_max": float(table["essentiality_loading"].max()),
        "imputation": "none; genes absent from a source are omitted, never filled",
    }
    Path(str(output) + ".manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crispr-gene-effect", required=True)
    parser.add_argument("--hallmark-gmt", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--restrict-to", default="", help="optional CSV with a 'gene' column")
    args = parser.parse_args()
    print(json.dumps(build(effect=args.crispr_gene_effect, gmt=args.hallmark_gmt,
                           output=args.output, restrict_to=args.restrict_to), indent=2))


if __name__ == "__main__":
    main()
