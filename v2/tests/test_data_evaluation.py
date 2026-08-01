from __future__ import annotations

import numpy as np
import pandas as pd

from morpheus.v2.baseline_exports import (export_cca_alignment_baseline,
                                           export_raw_hoptimus_baseline,
                                           export_ridge_alignment_baseline)
from morpheus.v2.canonical_registry import (build_canonical_registry,
                                             build_tcga_cdr_outcomes,
                                             survival_coverage)
from morpheus.v2.curated_panel import (CURATED_CANCER_PANEL,
                                       FROZEN_MECHANISM_PROGRAMMES,
                                       freeze_curated_panel,
                                       frozen_mechanism_programme_manifest)
from morpheus.v2.discovery_targets import build_matched_random_controls, score_gene_signatures
from morpheus.v2.paired_bootstrap import paired_patient_and_cancer_bootstrap, pearson_metric
from morpheus.v2.export_baselines import uncapped_patient_mean_std


def test_panel_only_resolves_its_explicit_reviewed_candidates() -> None:
    genes = [f"G{i}" for i in range(20)]
    first = CURATED_CANCER_PANEL[0]
    panel = freeze_curated_panel({first.candidates[0]: genes, "KEGG_MEDICUS_ARSENIC_EXPOSURE": genes}, genes)
    item = next(x for x in panel["targets"] if x["name"] == first.name)
    assert item["source_signature"] == first.candidates[0]
    assert item["status"] == "eligible"
    # An environmental/pathogen name cannot become an endpoint merely because
    # it appears first in a GMT file.
    assert all(x["source_signature"] != "KEGG_MEDICUS_ARSENIC_EXPOSURE" for x in panel["targets"])
    approved = {candidate for entry in CURATED_CANCER_PANEL for candidate in entry.candidates}
    assert all(x["source_signature"] in approved or x["source_signature"] is None for x in panel["targets"])
    assert set(FROZEN_MECHANISM_PROGRAMMES) == {"cytolytic_activity", "interferon_gamma", "antigen_presentation", "myeloid_macrophage", "stroma_caf", "tgf_beta_emt", "immune_exclusion", "hypoxia", "glycolysis", "angiogenesis", "proliferation", "dna_repair", "apoptosis_senescence", "emt", "mechanotransduction"}
    assert frozen_mechanism_programme_manifest()["digest"]


def test_competitive_scores_and_train_only_matched_controls() -> None:
    rna = pd.DataFrame({"patient_id": ["p1", "p2", "p3", "p4"],
                        "A": [10., 9., 1., 1.], "B": [9., 8., 1., 1.],
                        "C": [1., 2., 8., 9.], "D": [2., 1., 9., 8.]})
    scores, available, _ = score_gene_signatures(rna, {"signal": ["A", "B"]}, min_present_genes=2)
    assert available.signal.all() and scores.signal.notna().all()
    metadata = pd.DataFrame({"patient_id": rna.patient_id, "split": ["train", "train", "train", "test"]})
    controls, manifest = build_matched_random_controls(rna, {"signal": ["A", "B"]}, metadata, controls_per_target=2)
    assert len(controls) == 2
    assert all(not set(value) & {"A", "B"} for value in controls.values())
    assert manifest["fit_population"] == "development_train_only"


def test_matched_controls_ignore_rna_only_patients_outside_artifact_cohort() -> None:
    rna = pd.DataFrame({"patient_id": ["p1", "p2", "p3", "rna_only"],
                        "A": [10., 9., 1., 99.], "B": [9., 8., 1., 99.],
                        "C": [1., 2., 8., 0.], "D": [2., 1., 9., 0.]})
    metadata = pd.DataFrame({"patient_id": ["p1", "p2", "p3"],
                             "split": ["train", "train", "train"]})
    controls, manifest = build_matched_random_controls(rna, {"signal": ["A", "B"]}, metadata,
                                                        controls_per_target=1)
    assert len(controls) == 1
    assert manifest["n_source_rna_patients"] == 4
    assert manifest["n_canonical_rna_patients"] == 3
    assert manifest["n_noncanonical_rna_patients_excluded"] == 1


def test_registry_cdr_survival_and_baseline_exports(tmp_path) -> None:
    registry = build_canonical_registry({"wsi": pd.DataFrame({"patient_id": ["TCGA-AA-0001-01A", "TCGA-BB-0002"]}),
                                        "rna": pd.DataFrame({"patient_id": ["TCGA-AA-0001", "TCGA-CC-0003"], "cancer": ["A", "C"]})})
    assert registry.patient_id.tolist() == ["TCGA-AA-0001", "TCGA-BB-0002", "TCGA-CC-0003"]
    clinical = pd.DataFrame({"bcr_patient_barcode": ["TCGA-AA-0001-01A", "TCGA-BB-0002"],
                             "patient_id": ["0001", "0002"],
                             "vital_status": ["Dead", "LIVING"], "days_to_death": [100, np.nan],
                             "days_to_last_followup": [90, 200], "patient_progression_status": ["WITH TUMOR", "NO TUMOR"],
                             "days_to_patient_progression_free": [50, np.nan]})
    outcomes = build_tcga_cdr_outcomes(clinical)
    coverage = survival_coverage(outcomes, ["TCGA-AA-0001", "TCGA-BB-0002", "TCGA-CC-0003"], ["train", "test", "test"])
    assert coverage["os"]["test"]["n"] == 1
    ids = np.array([f"TCGA-AA-{i:04d}" for i in range(8)]); cancer = np.array(["A"] * 4 + ["B"] * 4); split = np.array(["train"] * 4 + ["test"] * 4)
    wsi = np.arange(64, dtype=np.float32).reshape(8, 8); rna = wsi[:, ::-1].copy()
    raw = export_raw_hoptimus_baseline(tmp_path / "raw.npz", patient_ids=ids, cancers=cancer, split=split, pooled_mean=wsi)
    ridge = export_ridge_alignment_baseline(tmp_path / "ridge.npz", patient_ids=ids, cancers=cancer, split=split, wsi=wsi, rna=rna, components=3)
    cca = export_cca_alignment_baseline(tmp_path / "cca.npz", patient_ids=ids, cancers=cancer, split=split, wsi=wsi, rna=rna, components=3)
    assert set(np.load(raw).files) >= {"wsi_identity", "trained_states"}
    assert set(np.load(ridge).files) >= {"wsi_identity", "rna_identity"}
    assert set(np.load(cca).files) >= {"wsi_identity", "rna_identity"}


def test_paired_patient_and_cancer_bootstrap_detects_improvement() -> None:
    rng = np.random.default_rng(0); y = rng.normal(size=80); teacher = rng.normal(size=80); challenger = y + rng.normal(scale=.05, size=80)
    result = paired_patient_and_cancer_bootstrap(pearson_metric, y, teacher, challenger, np.repeat(["A", "B", "C", "D"], 20), repeats=200, seed=3)
    assert result["patient"]["point_delta"] > 0.5
    assert result["cancer"]["p_improve"] > .9


def test_uncapped_patch_pooling_reads_every_token() -> None:
    class Store:
        def load_many_patient_tokens(self, patient_ids, max_tokens, seed, slide_balanced=True):
            assert max_tokens is None
            return [(np.arange((i + 2) * 3, dtype=np.float32).reshape(i + 2, 3), pd.DataFrame()) for i, _ in enumerate(patient_ids)]
    mean, std, count = uncapped_patient_mean_std(Store(), ["p1", "p2", "p3"], read_batch=2)
    assert count.tolist() == [2, 3, 2]  # each batch is independently handled
    assert mean.shape == std.shape == (3, 3)
