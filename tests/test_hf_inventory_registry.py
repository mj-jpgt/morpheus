import pandas as pd

from morpheus.src.data.hf_inventory import classify_modality
from morpheus.src.data.modality_registry import build_modality_registry


def test_modality_classification_core_modalities():
    assert classify_modality("clinical/follow_up.tsv") == "clinical"
    assert classify_modality("omics/cnv_table.parquet") == "cnv"
    assert classify_modality("somatic_mutation.maf.gz") == "snv"
    assert classify_modality("proteomics/phosphoproteome.csv") == "phosphoproteomics"
    assert classify_modality("slides/TCGA.svs") == "wsi"


def test_registry_summarizes_modalities():
    frame = pd.DataFrame(
        {
            "path": ["rna/a.csv", "wsi/a.svs"],
            "modality": ["rna", "wsi"],
            "is_table_file": [True, False],
            "size_bytes": [10, 100],
        }
    )
    registry = build_modality_registry(frame)
    assert registry["modalities"]["rna"]["available_in_inventory"]
    assert registry["modalities"]["rna"]["table_files_first"] == ["rna/a.csv"]
    assert registry["modalities"]["wsi"]["download_policy"] == "metadata_only_until_mapping"
