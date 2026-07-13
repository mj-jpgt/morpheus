import torch

from morpheus.src.encoders.rna_adapter import RNATokenAdapter
from morpheus.src.encoders.wsi_adapter import WSIPatchTokenAdapter
from morpheus.src.models.bio_query_former import BioQueryFormer, BioQueryFormerConfig


def test_bio_query_former_exposes_typed_surfaces():
    batch = 2
    hidden = 64
    wsi = WSIPatchTokenAdapter(8, hidden)
    rna = RNATokenAdapter(6, 4, hidden)
    tokens = {
        "wsi": wsi(
            torch.randn(batch, 5, 8),
            torch.tensor([[True, True, True, False, False], [True, True, True, True, True]]),
            torch.randn(batch, 5, 2),
        ),
        "rna": rna(torch.randn(batch, 6), torch.randn(batch, 4), bulk_present=torch.tensor([True, False])),
    }
    model = BioQueryFormer(
        BioQueryFormerConfig(
            hidden_dim=hidden,
            num_layers=1,
            num_heads=4,
            identity_slots=2,
            biology_slots=2,
            program_slots=3,
            wsi_residual_slots=1,
            rna_residual_slots=1,
            clinical_residual_slots=1,
            uncertainty_slots=1,
            hypothesis_slots=1,
        )
    )
    out = model(tokens)
    for key in [
        "identity_slots",
        "biology_slots",
        "program_slots",
        "z_identity",
        "z_biology",
        "z_programs",
        "z_wsi_residual",
        "z_rna_residual",
        "z_uncertainty",
        "z_hypothesis",
        "modality_masks",
    ]:
        assert key in out
    assert out["z_identity"].shape == (batch, hidden)
    assert out["z_biology"].shape == (batch, hidden)
    assert out["z_programs"].shape == (batch, 3, hidden)
