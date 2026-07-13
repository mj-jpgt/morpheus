from morpheus.src.data.bulkformer_zenodo import ZenodoFile, categorize_bulkformer_file, select_bulkformer_files


def _file(record: str, name: str) -> ZenodoFile:
    return ZenodoFile(
        record_id=record,
        record_title="Bulkformer",
        record_doi=None,
        version=None,
        filename=name,
        size=1,
        checksum=None,
        download_url=f"https://example.test/{name}",
        category=categorize_bulkformer_file(name),
    )


def test_bulkformer_file_categorization():
    assert categorize_bulkformer_file("TCGA_cancer_data.h5ad") == "tcga_h5ad"
    assert categorize_bulkformer_file("bulkformer_gene_info.csv") == "gene_info"
    assert categorize_bulkformer_file("Bulkformer_ckpt_epoch_29.pt") == "checkpoint"
    assert categorize_bulkformer_file("BulkFormerCode.zip") == "code_archive"
    assert categorize_bulkformer_file("G_tcga.pt") == "graph"
    assert categorize_bulkformer_file("G_tcga_weight.pt") == "graph_weight"
    assert categorize_bulkformer_file("esm2_feature_concat.pt") == "gene_embedding"
    assert categorize_bulkformer_file("interested_gene_list.pt") == "interested_gene_list"
    assert categorize_bulkformer_file("gene_length_df.csv") == "gene_length"


def test_selection_prefers_primary_when_category_exists():
    files = [
        _file("15744294", "TCGA_survival.h5ad"),
        _file("15559368", "TCGA_cancer_data.h5ad"),
        _file("15559368", "Bulkformer_ckpt_epoch_29.pt"),
    ]
    selected = select_bulkformer_files(files, "15744294")
    assert selected["tcga_h5ad"].record_id == "15744294"
    assert selected["checkpoint"].record_id == "15559368"
