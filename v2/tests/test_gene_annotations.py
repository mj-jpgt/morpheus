"""D2.3's annotation table decides whether the proliferation deflation is visible."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.build_gene_annotations import (PROLIFERATION_SETS, _symbol, build,
                                                read_depmap_essentiality,
                                                read_hallmark_proliferation)


def _gmt(tmp_path: Path, extra: str = "") -> Path:
    lines = [f"{name}\thttp://x\tGENE{i}\tSHARED{i}" for i, name in enumerate(PROLIFERATION_SETS)]
    lines.append("HALLMARK_APOPTOSIS\thttp://x\tNOTPROLIF1\tNOTPROLIF2")
    path = tmp_path / "h.gmt"; path.write_text("\n".join(lines) + extra + "\n", encoding="utf-8")
    return path


def _effect(tmp_path: Path) -> Path:
    frame = pd.DataFrame({
        "GENE0 (1)": [-2.0, -1.0], "SHARED0 (2)": [-0.5, -0.5],
        "NOTPROLIF1 (3)": [0.1, -0.1], "ONLYDEPMAP (4)": [-1.5, -1.5]},
        index=["ACH-1", "ACH-2"])
    path = tmp_path / "eff.csv"; frame.to_csv(path); return path


def test_depmap_column_symbols_are_parsed():
    assert _symbol("A1BG (1)") == "A1BG" and _symbol(" tp53 (7157) ") == "TP53"


def test_essentiality_is_sign_flipped_so_larger_is_more_essential(tmp_path):
    ess = read_depmap_essentiality(_effect(tmp_path))
    # GENE0 mean effect -1.5 (strongly essential) must outrank NOTPROLIF1 mean 0.0
    assert ess["GENE0"] == pytest.approx(1.5)
    assert ess["NOTPROLIF1"] == pytest.approx(0.0)
    assert ess["GENE0"] > ess["NOTPROLIF1"], "sign flip lost: essential genes must score HIGHER"


def test_missing_proliferation_set_is_rejected(tmp_path):
    bad = tmp_path / "bad.gmt"
    bad.write_text("HALLMARK_E2F_TARGETS\tx\tA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing declared proliferation sets"):
        read_hallmark_proliferation(bad)


def test_build_emits_finite_annotations_and_flags_proliferation(tmp_path):
    out = tmp_path / "ann.parquet"
    manifest = build(effect=str(_effect(tmp_path)), gmt=str(_gmt(tmp_path)), output=str(out))
    table = pd.read_parquet(out)
    assert set(table.columns) == {"gene", "essentiality_loading", "proliferation_loading"}
    assert np.isfinite(table[["essentiality_loading", "proliferation_loading"]].to_numpy()).all()
    assert not table["gene"].duplicated().any()
    flags = dict(zip(table["gene"], table["proliferation_loading"]))
    assert flags["GENE0"] == 1.0 and flags["SHARED0"] == 1.0
    assert flags["NOTPROLIF1"] == 0.0, "a non-proliferation gene must not be flagged"
    assert flags["ONLYDEPMAP"] == 0.0
    assert manifest["imputation"] == "none; genes absent from a source are omitted, never filled"
    assert json.loads((tmp_path / "ann.parquet.manifest.json").read_text())["n_genes"] == len(table)


def test_genes_absent_from_depmap_are_omitted_not_imputed(tmp_path):
    """A gene in the GMT but not in DepMap has no measured essentiality. Inventing one
    would silently fabricate the very quantity the proliferation caveat turns on."""
    table = pd.read_parquet(tmp_path / "x.parquet") if False else None
    out = tmp_path / "x.parquet"
    build(effect=str(_effect(tmp_path)), gmt=str(_gmt(tmp_path)), output=str(out))
    table = pd.read_parquet(out)
    assert "GENE1" not in set(table["gene"]), "GMT-only gene must be omitted, not imputed"
