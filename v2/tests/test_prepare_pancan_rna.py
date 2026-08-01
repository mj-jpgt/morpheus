"""Focused regression tests for the streaming PanCan RNA converter."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.prepare_pancan_rna import RNAPreparationError, prepare_pancan_rna


def test_streaming_converter_collapses_samples_and_duplicate_symbols(tmp_path) -> None:
    source = tmp_path / "pancan.tsv"
    source.write_text(
        '"gene"\t"TCGA-AA-0001-01A"\t"TCGA-AA-0001-11A"\t"TCGA-BB-0002-01A"\n'
        '"A1BG|1"\t1\t3\t5\n'
        '"ENSG00000121410|A1BG|1"\t3\t5\t7\n'
        '"B2M|567"\t2\t\t8\n'
        '"?|0"\t9\t9\t9\n',
        encoding="utf-8",
    )
    output = tmp_path / "patient_by_gene.parquet"
    manifest = prepare_pancan_rna(source, output, chunk_genes=2, patient_chunk=1)
    result = pd.read_parquet(output).set_index("patient_id")
    # A1BG has four values for AA across its two source rows: (1+3+3+5)/4.
    assert result.loc["TCGA-AA-0001", "A1BG"] == pytest.approx(3.0)
    assert result.loc["TCGA-BB-0002", "A1BG"] == pytest.approx(6.0)
    assert result.loc["TCGA-AA-0001", "B2M"] == pytest.approx(2.0)
    assert result.loc["TCGA-BB-0002", "B2M"] == pytest.approx(8.0)
    assert result.dtypes["A1BG"] == np.dtype("float32")
    assert manifest["canonical_patient_count"] == 2
    assert manifest["duplicate_gene_symbols_collapsed"] == 1
    assert manifest["dropped_unmapped_gene_rows"] == 1
    sidecar = json.loads((tmp_path / "patient_by_gene.manifest.json").read_text(encoding="utf-8"))
    assert sidecar["source_sha256"] == manifest["source_sha256"]
    assert len(pd.read_parquet(tmp_path / "patient_by_gene.sample_map.parquet")) == 3


def test_converter_rejects_ambiguous_gene_identity(tmp_path) -> None:
    source = tmp_path / "ambiguous.tsv"
    source.write_text("gene\tTCGA-AA-0001-01A\nA1BG|B2M|1\t1\n", encoding="utf-8")
    with pytest.raises(RNAPreparationError, match="ambiguous gene label"):
        prepare_pancan_rna(source, tmp_path / "out.parquet")
