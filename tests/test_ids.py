from morpheus.src.utils.ids import normalize_patient_id, parse_tcga_barcode


def test_normalize_patient_id_from_slide_filename():
    assert normalize_patient_id("TCGA-02-0001-01Z-00-DX1.npz") == "TCGA-02-0001"


def test_parse_tcga_barcode_components():
    parsed = parse_tcga_barcode("TCGA-AB-1234-01Z-00-DX1.svs")
    assert parsed.patient_id == "TCGA-AB-1234"
    assert parsed.sample_id == "TCGA-AB-1234-01"
    assert parsed.slide_id.startswith("TCGA-AB-1234")
