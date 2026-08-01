from __future__ import annotations

import json

import numpy as np

from morpheus.v2.contracts import require_state, validate_artifact
from morpheus.v2.refit_mlp_clip import MLPClipRefitConfig, refit_and_export, train_mlp_clip


def _data(n: int = 18):
    rng = np.random.default_rng(17)
    ids = np.asarray([f"TCGA-AA-{index:04d}" for index in range(n)])
    cancers = np.asarray(["A"] * 6 + ["B"] * 6 + ["C"] * (n - 12))
    split = np.asarray(["train"] * 6 + ["val"] * 6 + ["test"] * (n - 12))
    latent = rng.normal(size=(n, 3)).astype(np.float32)
    mean = np.concatenate([latent @ rng.normal(size=(3, 5)), rng.normal(size=(n, 3))], axis=1).astype(np.float32)
    std = np.abs(rng.normal(size=mean.shape)).astype(np.float32)
    rna = np.concatenate([latent @ rng.normal(size=(3, 7)), rng.normal(size=(n, 2))], axis=1).astype(np.float32)
    return ids, cancers, split, mean, std, rna


def test_train_mlp_clip_is_normalized_and_uses_only_fit_rows() -> None:
    _, cancers, _, mean, std, rna = _data()
    config = MLPClipRefitConfig(embedding_dim=4, pca_dim=4, hidden_multiplier=1, epochs=2, batch_size=4, seed=3, use_bf16=False)
    z_wsi, z_rna, training, state = train_mlp_clip(wsi=np.c_[mean, std], rna=rna, cancers=cancers, fit_rows=np.arange(6), config=config, device="cpu")
    assert z_wsi.shape == z_rna.shape == (18, 4)
    assert np.allclose(np.linalg.norm(z_wsi, axis=1), 1.0, atol=1e-5)
    assert training["n_fit_patients"] == 6
    assert state["fit_rows"] == list(range(6))


def test_siglip_and_hard_negative_variants_export_valid_embeddings() -> None:
    _, cancers, _, mean, std, rna = _data()
    for family in ("siglip", "hard_negative_clip"):
        config = MLPClipRefitConfig(embedding_dim=4, pca_dim=4, hidden_multiplier=1, epochs=1,
                                    batch_size=4, seed=3, use_bf16=False, loss_family=family)
        z_wsi, z_rna, _, _ = train_mlp_clip(wsi=np.c_[mean, std], rna=rna, cancers=cancers,
                                            fit_rows=np.arange(6), config=config, device="cpu")
        assert np.isfinite(z_wsi).all() and np.isfinite(z_rna).all()
        assert np.allclose(np.linalg.norm(z_wsi, axis=1), 1.0, atol=1e-5)


def test_final_refit_exports_declared_teacher_contract(tmp_path) -> None:
    ids, cancers, split, mean, std, rna = _data()
    path = tmp_path / "teacher.npz"
    checkpoint = tmp_path / "teacher.pt"
    result = refit_and_export(patient_ids=ids, cancers=cancers, split=split, pooled_mean=mean, pooled_std=std, rna=rna,
                              output=path, fit_population="development_train_val",
                              config=MLPClipRefitConfig(embedding_dim=4, pca_dim=4, hidden_multiplier=1, epochs=2, batch_size=4, use_bf16=False),
                              checkpoint=checkpoint)
    assert result["fit_population"] == "development_train_val" and checkpoint.exists()
    schema = validate_artifact(path)
    assert schema["trained_states"] == ["rna_identity", "wsi_identity"]
    with np.load(path, allow_pickle=False) as raw:
        manifest = json.loads(raw["manifest_json"].item())
        assert manifest["fit_population"] == "development_train_val"
        assert manifest["config"]["pooling"] == "uncapped_patient_mean_plus_std"
        assert require_state(raw, "wsi_identity").shape == (18, 4)
        assert require_state(raw, "full_patient") is None
