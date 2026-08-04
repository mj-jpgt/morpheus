"""Self-tests for the SNV target builder.

The failure modes this module exists to prevent are all silent ones: a mangled
patient barcode joins to nothing and looks like missing data; an unsequenced
patient imputed to wild-type looks like a real negative; a gene selected using
held-out patients looks like a legitimate target.  Each is planted below and
asserted to be caught, and where a rule has two directions both are tested --- a
check that only ever refuses is not checking anything.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pytest

from morpheus.v2.build_snv_targets import (ACCEPTED_FILTER_TAGS, PROTEIN_ALTERING, accepted_filter,
                                           build_snv_targets, canonical_patient,
                                           matched_random_controls, select_development_genes)

_HEADER = "Hugo_Symbol\tVariant_Classification\tTumor_Sample_Barcode\tFILTER\n"


def _write_maf(path: Path, rows: list[tuple[str, str, str, str]]) -> Path:
    with gzip.open(path, "wt") as handle:
        handle.write(_HEADER)
        for row in rows:
            handle.write("\t".join(row) + "\n")
    return path


def _write_reference(path: Path, patients, cancers, splits) -> Path:
    np.savez_compressed(path, patient_ids=np.asarray(patients, dtype=str),
                        cancers=np.asarray(cancers, dtype=str), split=np.asarray(splits, dtype=str),
                        target_names=np.asarray(["X"]), target_groups=np.asarray(["g"]),
                        scores=np.zeros((len(patients), 1), dtype=np.float32))
    return path


def test_canonical_patient_reduces_an_aliquot_and_refuses_a_non_tcga_identifier():
    assert canonical_patient("TCGA-02-0047-01A-01D-1490-08") == "TCGA-02-0047"
    assert canonical_patient("tcga-05-4397-01a") == "TCGA-05-4397"
    for bad in ("SAMPLE-1", "TCGA-02", "TCGA--0047", ""):
        with pytest.raises(ValueError):
            canonical_patient(bad)


def test_accepted_filter_takes_the_mc3_rescue_tags_and_rejects_everything_else():
    assert accepted_filter("PASS")
    assert accepted_filter("wga")
    assert accepted_filter("PASS,wga")
    assert accepted_filter("wga,native_wga_mix")
    # a single unaccepted tag disqualifies the row even when combined with PASS
    assert not accepted_filter("PASS,oxog")
    assert not accepted_filter("StrandBias")
    assert not accepted_filter("")
    assert not accepted_filter("nan")
    assert ACCEPTED_FILTER_TAGS == {"PASS", "wga", "native_wga_mix"}


def test_select_development_genes_ignores_test_rows_in_both_directions():
    """A gene recurrent only in test must not be selected; one recurrent only in
    development must be, even if it is absent from test entirely."""
    development = np.array([True] * 100 + [False] * 100)
    only_test = {"TESTONLY": set(range(100, 200))}
    with pytest.raises(ValueError):
        select_development_genes(only_test, development, min_prevalence=0.1, min_count=5)
    only_dev = {"DEVONLY": set(range(0, 40))}
    assert select_development_genes(only_dev, development, min_prevalence=0.1, min_count=5) == ["DEVONLY"]


def test_select_development_genes_orders_by_development_count_then_symbol():
    development = np.ones(100, dtype=bool)
    genes = {"LOW": set(range(10)), "HIGH": set(range(50)), "MID_B": set(range(20)),
             "MID_A": set(range(20))}
    assert select_development_genes(genes, development, min_prevalence=0.05,
                                    min_count=5) == ["HIGH", "MID_A", "MID_B", "LOW"]


def test_matched_random_controls_preserve_the_marginal_and_destroy_the_pairing():
    rng = np.random.default_rng(0)
    scores = np.column_stack([rng.integers(0, 2, 500).astype(np.float32), rng.normal(size=500).astype(np.float32)])
    control, names = matched_random_controls(scores, ["A", "B"], seed=42)
    assert names == ["RANDOM_CONTROL__A__0", "RANDOM_CONTROL__B__0"]
    for column in range(scores.shape[1]):
        assert np.allclose(np.sort(control[:, column]), np.sort(scores[:, column]))
    # the pairing is gone: correlation with the real column is not 1
    assert abs(np.corrcoef(control[:, 1], scores[:, 1])[0, 1]) < 0.5


def test_matched_random_controls_permute_each_column_independently():
    """A single shared permutation would preserve the cross-target covariance and
    make the control block look like the real one under a multivariate readout."""
    values = np.column_stack([np.arange(200.0), np.arange(200.0)]).astype(np.float32)
    control, _ = matched_random_controls(values, ["A", "B"], seed=1)
    assert not np.allclose(control[:, 0], control[:, 1])


def test_build_emits_the_frozen_target_schema_and_excludes_unsequenced_patients(tmp_path: Path):
    patients = [f"TCGA-AA-{i:04d}" for i in range(120)]
    splits = ["train"] * 80 + ["test"] * 40
    cancers = ["COAD"] * 60 + ["LUAD"] * 60
    reference = _write_reference(tmp_path / "reference.npz", patients, cancers, splits)

    rows = []
    for index, patient in enumerate(patients):
        if index >= 110:          # last ten patients are never sequenced
            continue
        rows.append(("TP53", "Missense_Mutation", f"{patient}-01A-11D-0001-01", "PASS"))
        if index % 2 == 0:
            rows.append(("KRAS", "Nonsense_Mutation", f"{patient}-01A-11D-0001-01", "wga"))
        rows.append(("QUIET", "Silent", f"{patient}-01A-11D-0001-01", "PASS"))
        rows.append(("DROPPED", "Missense_Mutation", f"{patient}-01A-11D-0001-01", "oxog"))
    maf = _write_maf(tmp_path / "test.maf.gz", rows)

    output = tmp_path / "snv_targets_test.npz"
    manifest = build_snv_targets(maf, reference, output, min_prevalence=0.1, min_development_count=5,
                                 seed=42, chunk_rows=17)

    raw = np.load(output, allow_pickle=True)
    assert set(raw.files) == {"patient_ids", "cancers", "split", "target_names", "target_groups",
                              "scores", "metadata_json"}
    assert len(raw["patient_ids"]) == 110
    assert manifest["cohort"]["excluded_no_mc3_record_count"] == 10
    assert manifest["cohort"]["excluded_no_mc3_record_patients"] == patients[110:]

    names = [str(n) for n in raw["target_names"]]
    groups = [str(g) for g in raw["target_groups"]]
    assert "SNV_TP53" in names and "SNV_KRAS" in names
    # a Silent-only gene and a gene whose only rows failed FILTER must not be targets
    assert "SNV_QUIET" not in names and "SNV_DROPPED" not in names
    assert names.count("SNV_BURDEN_LOG1P") == 1 and names.count("SNV_BURDEN_LOG1P_ALL") == 1
    assert groups.count("random_control") == len(names) - groups.count("random_control")
    assert raw["scores"].dtype == np.float32
    assert raw["scores"].shape == (110, len(names))
    assert np.isfinite(raw["scores"]).all()

    scores = raw["scores"]
    assert scores[:, names.index("SNV_TP53")].sum() == 110
    assert scores[:, names.index("SNV_KRAS")].sum() == 55
    # the FILTER-rejected and Silent rows are excluded from the protein-altering burden
    burden = scores[:, names.index("SNV_BURDEN_LOG1P")]
    assert burden[0] == pytest.approx(np.log1p(2.0))
    assert burden[1] == pytest.approx(np.log1p(1.0))
    # ...but the all-classification burden still counts the Silent row it kept
    assert scores[0, names.index("SNV_BURDEN_LOG1P_ALL")] == pytest.approx(np.log1p(3.0))

    embedded = json.loads(str(raw["metadata_json"]))
    assert embedded["gene_selection"]["partition"].startswith("train+val")
    assert embedded["variant_filter"]["protein_altering"] == list(PROTEIN_ALTERING)
    assert embedded["maf_sha256"] and embedded["reference_artifact_sha256"]
    assert Path(str(output) + ".manifest.json").is_file()


def test_build_is_byte_identical_across_chunk_sizes(tmp_path: Path):
    """Chunking is an implementation detail; if it moves a digest it is a bug."""
    patients = [f"TCGA-BB-{i:04d}" for i in range(60)]
    reference = _write_reference(tmp_path / "reference.npz", patients, ["BRCA"] * 60,
                                 ["train"] * 40 + ["test"] * 20)
    rows = [("PIK3CA", "Frame_Shift_Del", f"{p}-01A-11D-0001-01", "PASS") for p in patients]
    maf = _write_maf(tmp_path / "test.maf.gz", rows)
    digests = []
    for chunk in (7, 1000):
        manifest = build_snv_targets(maf, reference, tmp_path / f"out_{chunk}.npz",
                                     min_prevalence=0.1, min_development_count=5, chunk_rows=chunk)
        digests.append(manifest["digests"])
    assert digests[0] == digests[1]


def test_build_refuses_a_cohort_that_shares_no_patient_with_the_maf(tmp_path: Path):
    reference = _write_reference(tmp_path / "reference.npz", [f"TCGA-CC-{i:04d}" for i in range(30)],
                                 ["GBM"] * 30, ["train"] * 30)
    maf = _write_maf(tmp_path / "test.maf.gz",
                     [("TP53", "Missense_Mutation", "TCGA-ZZ-9999-01A-11D-0001-01", "PASS")])
    with pytest.raises(ValueError, match="no cohort patient"):
        build_snv_targets(maf, reference, tmp_path / "out.npz", min_prevalence=0.1,
                          min_development_count=1)


def test_build_refuses_when_no_gene_clears_the_development_rule(tmp_path: Path):
    """An empty block would otherwise be emitted and read downstream as a modality."""
    patients = [f"TCGA-DD-{i:04d}" for i in range(50)]
    reference = _write_reference(tmp_path / "reference.npz", patients, ["LUSC"] * 50,
                                 ["train"] * 30 + ["test"] * 20)
    rows = [("RARE", "Missense_Mutation", f"{patients[0]}-01A-11D-0001-01", "PASS")]
    rows += [("OTHER", "Silent", f"{p}-01A-11D-0001-01", "PASS") for p in patients]
    maf = _write_maf(tmp_path / "test.maf.gz", rows)
    with pytest.raises(ValueError, match="no gene met"):
        build_snv_targets(maf, reference, tmp_path / "out.npz", min_prevalence=0.5,
                          min_development_count=20)
