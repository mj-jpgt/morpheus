import torch

from morpheus.src.encoders.clinical_adapter import ClinicalTokenAdapter
from morpheus.src.encoders.rna_adapter import RNATokenAdapter
from morpheus.src.encoders.wsi_adapter import WSITokenAdapter
from morpheus.src.models.tumor_state_query_former import QueryFormerConfig, TumorStateQueryFormer


def test_adapters_and_query_former_handle_missing_modalities():
    batch = 3
    hidden = 64
    wsi = WSITokenAdapter(8, hidden)
    rna = RNATokenAdapter(6, 4, hidden)
    clinical = ClinicalTokenAdapter(5, hidden)
    tokens = {
        "wsi": wsi(torch.randn(batch, 8), torch.tensor([True, False, True])),
        "rna": rna(torch.randn(batch, 6), torch.randn(batch, 4), gene_set_present=torch.tensor([True, True, False])),
        "clinical": clinical(torch.randn(batch, 5), torch.tensor([False, True, True])),
    }
    model = TumorStateQueryFormer(QueryFormerConfig(hidden_dim=hidden, num_layers=1, num_heads=4, shared_slots=2, wsi_residual_slots=1, rna_residual_slots=1, clinical_residual_slots=1, genomic_residual_slots=1, uncertainty_slots=1, task_slots=1))
    out = model(tokens)
    for key in ["z_shared", "z_wsi_resid", "z_rna_resid", "z_clinical_resid", "z_genomic_resid", "z_uncertainty", "z_task", "z_patient", "modality_masks"]:
        assert key in out
    assert out["z_patient"].shape == (batch, hidden)
    assert out["z_shared"].shape == (batch, 2, hidden)
