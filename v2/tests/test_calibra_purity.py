from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.calibra.run_calibra import load_purity_table, validate_evaluated_targets


def test_purity_loader_canonicalises_sample_ids_and_rejects_duplicate_patients(tmp_path):
    path = tmp_path / "purity.tsv"
    pd.DataFrame({"sample_id": ["TCGA-AB-1234-01A", "TCGA-CD-5678"], "purity": [.4, .8]}).to_csv(path, sep="\t", index=False)
    values, manifest = load_purity_table(str(path), source="published_consensus", reference="doi:test")
    assert values == {"TCGA-AB-1234": .4, "TCGA-CD-5678": .8}
    assert manifest["n_canonical_patients"] == 2
    pd.DataFrame({"patient_id": ["TCGA-AB-1234", "TCGA-AB-1234-01A"], "purity": [.4, .5]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="duplicate canonical"):
        load_purity_table(str(path), source="published_consensus", reference="doi:test")


def test_purity_loader_requires_declared_units_source_and_reference(tmp_path):
    path = tmp_path / "purity.tsv"
    pd.DataFrame({"patient_id": ["TCGA-AB-1234"], "purity": [80.]}).to_csv(path, sep="\t", index=False)
    values, manifest = load_purity_table(str(path), units="percent", source="expression_derived", reference="ESTIMATE v1")
    assert values["TCGA-AB-1234"] == .8
    assert manifest["circularity"] == "expression_derived_from_same_RNA_targets"
    with pytest.raises(ValueError, match="mandatory"):
        load_purity_table(str(path))


def test_evaluated_targets_reject_constant_or_nonfinite_columns():
    with pytest.raises(ValueError, match="constant"):
        validate_evaluated_targets(np.asarray([[1., 3.], [2., 3.]]), np.asarray(["A", "B"]))
    with pytest.raises(ValueError, match="finite"):
        validate_evaluated_targets(np.asarray([[1., np.nan], [2., 3.]]), np.asarray(["A", "B"]))
