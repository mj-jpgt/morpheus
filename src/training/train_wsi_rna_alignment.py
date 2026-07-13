"""Train/evaluate WSI-RNA alignment baselines for v1 proof."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np
import pandas as pd

from morpheus.src.eval.retrieval_metrics import paired_retrieval_metrics
from morpheus.src.encoders.bulkformer_encoder import load_embedding_store
from morpheus.src.utils.config import load_config
from morpheus.src.utils.ids import normalize_patient_id
from morpheus.src.utils.provenance import base_manifest, write_json


@dataclass
class PairedMatrix:
    patient_ids: list[str]
    cancers: list[str]
    wsi: np.ndarray
    rna: np.ndarray
    split: list[str]


def _load_wsi_features(cfg) -> pd.DataFrame:
    standard_path = cfg.path("wsi_standard_dir") / "tcga_ut_hoptimus0_patient_embeddings.h5"
    if standard_path.exists():
        import h5py

        with h5py.File(standard_path, "r") as handle:
            arr = handle["embeddings"][:].astype(np.float32)
            ids = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in handle["patient_ids"][:]]
        return pd.DataFrame({"patient_id": ids, "wsi_vector": [row for row in arr]})
    feature_path = cfg.path("processed_dir") / "features" / "wsi_features.npz"
    meta_path = cfg.path("processed_dir") / "features" / "wsi_metadata.parquet"
    if not feature_path.exists() or not meta_path.exists():
        raise FileNotFoundError("Run preprocess_wsi_embeddings first")
    arr = np.load(feature_path)["embeddings"]
    meta = pd.read_parquet(meta_path).copy()
    meta["row_idx"] = np.arange(len(meta))
    # Patient-level H-Optimus path may already be unique. Patch path can have multiple slides.
    grouped = meta.groupby("patient_id", dropna=True)["row_idx"].apply(list)
    rows = []
    for patient_id, idxs in grouped.items():
        vec = arr[np.asarray(idxs)].mean(axis=0)
        rows.append({"patient_id": patient_id, "wsi_vector": vec})
    return pd.DataFrame(rows)


def _load_rna_features(cfg) -> pd.DataFrame:
    bulkformer = cfg.path("rna_bulkformer_embeddings")
    if bulkformer.exists():
        df = load_embedding_store(bulkformer)
        numeric = df.select_dtypes(include=["number"]).columns.tolist()
        frame = pd.DataFrame(
            {
                "patient_id": df["patient_id"].map(normalize_patient_id),
                "rna_vector": [row.astype(np.float32) for row in df[numeric].to_numpy()],
            }
        ).dropna(subset=["patient_id"])
        return _collapse_vectors(frame, "rna_vector")
    fallback = cfg.path("rna_geneformer_embeddings")
    if fallback.exists():
        df = pd.read_parquet(fallback)
        col = next((c for c in df.columns if "patient" in str(c).lower() or str(c).lower() in {"pid", "patient_id"}), None)
        if col is None:
            col = df.columns[0]
        numeric = df.select_dtypes(include=["number"]).columns.tolist()
        if not numeric:
            raise ValueError(f"No numeric RNA embedding columns found in {fallback}")
        frame = pd.DataFrame(
            {
                "patient_id": df[col].map(normalize_patient_id),
                "rna_vector": [row.astype(np.float32) for row in df[numeric].to_numpy()],
            }
        ).dropna(subset=["patient_id"])
        return _collapse_vectors(frame, "rna_vector")
    processed = cfg.path("rna_processed")
    df = pd.read_parquet(processed)
    col = next((c for c in df.columns if "patient" in str(c).lower()), df.columns[0])
    numeric = df.select_dtypes(include=["number"]).columns.tolist()[:4096]
    frame = pd.DataFrame({"patient_id": df[col].map(normalize_patient_id), "rna_vector": [row.astype(np.float32) for row in df[numeric].to_numpy()]}).dropna(subset=["patient_id"])
    return _collapse_vectors(frame, "rna_vector")


def _collapse_vectors(frame: pd.DataFrame, vector_col: str) -> pd.DataFrame:
    rows = []
    for patient_id, group in frame.groupby("patient_id", dropna=True):
        matrix = np.vstack(group[vector_col].to_numpy()).astype(np.float32)
        rows.append({"patient_id": patient_id, vector_col: matrix.mean(axis=0)})
    return pd.DataFrame(rows)


def load_paired_matrix(config_path: str = "morpheus/configs/v1.json", split_file: str | Path | None = None) -> PairedMatrix:
    if split_file is not None:
        from morpheus.src.training.train_query_former import load_query_former_data

        data = load_query_former_data(config_path, split_file)
        return PairedMatrix(
            patient_ids=data.patient_ids,
            cancers=data.cancers,
            wsi=data.wsi,
            rna=data.rna,
            split=data.split.tolist(),
        )
    cfg = load_config(config_path)
    master = pd.read_parquet(cfg.path("processed_dir") / "master_patient_table.parquet")
    wsi = _load_wsi_features(cfg)
    rna = _load_rna_features(cfg)
    pairs = master[["patient_id", "cancer_type"]].merge(wsi, on="patient_id", how="inner").merge(rna, on="patient_id", how="inner")
    if pairs.empty:
        raise ValueError("No paired WSI-RNA patients found")
    split = _split_labels(cfg, pairs["patient_id"].astype(str).tolist())
    return PairedMatrix(
        patient_ids=pairs["patient_id"].astype(str).tolist(),
        cancers=pairs["cancer_type"].astype(str).tolist(),
        wsi=np.vstack(pairs["wsi_vector"].to_numpy()).astype(np.float32),
        rna=np.vstack(pairs["rna_vector"].to_numpy()).astype(np.float32),
        split=split,
    )


def _split_labels(cfg, patient_ids: list[str]) -> list[str]:
    split_path = cfg.path("processed_dir") / "splits" / "stratified_pan_cancer_split.json"
    labels = {pid: "train" for pid in patient_ids}
    if not split_path.exists():
        return [labels[pid] for pid in patient_ids]
    split = json.load(open(split_path, encoding="utf-8")).get("patient_ids", {})
    for name in ("train", "val", "test"):
        for pid in split.get(name, []):
            labels[str(pid)] = name
    return [labels.get(pid, "train") for pid in patient_ids]


def _fit_pca(train: np.ndarray, n_components: int):
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    n = min(n_components, train.shape[0] - 1, train.shape[1])
    scaler = StandardScaler()
    train_s = scaler.fit_transform(train)
    pca = PCA(n_components=max(1, n), random_state=42)
    pca.fit(train_s)
    return scaler, pca


def _transform_pca(scaler, pca, x: np.ndarray) -> np.ndarray:
    return pca.transform(scaler.transform(x)).astype(np.float32)


def _train_mask(data: PairedMatrix) -> np.ndarray:
    return np.asarray([s == "train" for s in data.split], dtype=bool)


def _split_mask(data: PairedMatrix, split_name: str) -> np.ndarray:
    return np.asarray([s == split_name for s in data.split], dtype=bool)


def align_none(data: PairedMatrix) -> tuple[np.ndarray, np.ndarray, dict]:
    dim = min(data.wsi.shape[1], data.rna.shape[1])
    return data.wsi[:, :dim], data.rna[:, :dim], {"method": "dimension_truncated_cosine", "dim": dim}


def align_ridge(data: PairedMatrix, dim: int) -> tuple[np.ndarray, np.ndarray, dict]:
    from sklearn.linear_model import Ridge

    train = _train_mask(data)
    sx, px = _fit_pca(data.wsi[train], min(dim, data.wsi.shape[1]))
    sy, py = _fit_pca(data.rna[train], min(dim, data.rna.shape[1]))
    x = _transform_pca(sx, px, data.wsi)
    y = _transform_pca(sy, py, data.rna)
    model = Ridge(alpha=1.0)
    model.fit(x[train], y[train])
    return model.predict(x).astype(np.float32), y.astype(np.float32), {"method": "ridge", "alpha": 1.0, "dim": x.shape[1]}


def align_ridge_grid(data: PairedMatrix, dim: int) -> tuple[np.ndarray, np.ndarray, dict]:
    from sklearn.linear_model import Ridge

    train = _train_mask(data)
    val = _split_mask(data, "val")
    best = None
    for pca_dim in (64, 128, 256, 512):
        sx, px = _fit_pca(data.wsi[train], min(pca_dim, data.wsi.shape[1]))
        sy, py = _fit_pca(data.rna[train], min(pca_dim, data.rna.shape[1]))
        x = _transform_pca(sx, px, data.wsi)
        y = _transform_pca(sy, py, data.rna)
        for alpha in (0.01, 0.1, 1.0, 10.0, 100.0):
            model = Ridge(alpha=alpha)
            model.fit(x[train], y[train])
            z_wsi = model.predict(x).astype(np.float32)
            z_rna = y.astype(np.float32)
            score = _selection_score(data, z_wsi, z_rna, val)
            if best is None or score > best[0]:
                best = (score, z_wsi, z_rna, {"method": "ridge_grid", "alpha": alpha, "pca_dim": int(x.shape[1]), "selected_by": "val_retrieval_r10_plus_mrr"})
    assert best is not None
    return best[1], best[2], best[3]


def align_cca(data: PairedMatrix, dim: int) -> tuple[np.ndarray, np.ndarray, dict]:
    from sklearn.cross_decomposition import CCA

    train = _train_mask(data)
    pca_dim = 128
    cca_dim = 50
    sx, px = _fit_pca(data.wsi[train], min(pca_dim, data.wsi.shape[1]))
    sy, py = _fit_pca(data.rna[train], min(pca_dim, data.rna.shape[1]))
    x = _transform_pca(sx, px, data.wsi)
    y = _transform_pca(sy, py, data.rna)
    n = min(cca_dim, x.shape[1], y.shape[1], max(1, int(train.sum()) - 2))
    cca = CCA(n_components=n, max_iter=1000)
    cca.fit(x[train], y[train])
    x_c, y_c = cca.transform(x, y)
    return x_c.astype(np.float32), y_c.astype(np.float32), {"method": "cca", "pca_dim": pca_dim, "cca_dim": n}


def align_cca_grid(data: PairedMatrix, dim: int) -> tuple[np.ndarray, np.ndarray, dict]:
    from sklearn.cross_decomposition import CCA

    train = _train_mask(data)
    val = _split_mask(data, "val")
    best = None
    for pca_dim in (64, 128, 256):
        sx, px = _fit_pca(data.wsi[train], min(pca_dim, data.wsi.shape[1]))
        sy, py = _fit_pca(data.rna[train], min(pca_dim, data.rna.shape[1]))
        x = _transform_pca(sx, px, data.wsi)
        y = _transform_pca(sy, py, data.rna)
        for cca_dim in (16, 32, 50, 64, 100):
            n = min(cca_dim, x.shape[1], y.shape[1], max(1, int(train.sum()) - 2))
            try:
                cca = CCA(n_components=n, max_iter=2000)
                cca.fit(x[train], y[train])
                x_c, y_c = cca.transform(x, y)
            except Exception:
                continue
            z_wsi = x_c.astype(np.float32)
            z_rna = y_c.astype(np.float32)
            score = _selection_score(data, z_wsi, z_rna, val)
            if best is None or score > best[0]:
                best = (score, z_wsi, z_rna, {"method": "cca_grid", "pca_dim": int(x.shape[1]), "cca_dim": int(n), "selected_by": "val_retrieval_r10_plus_mrr"})
    if best is None:
        return align_cca(data, dim)
    return best[1], best[2], best[3]


def align_pls_grid(data: PairedMatrix, dim: int) -> tuple[np.ndarray, np.ndarray, dict]:
    from sklearn.cross_decomposition import PLSRegression

    train = _train_mask(data)
    val = _split_mask(data, "val")
    best = None
    for pca_dim in (64, 128, 256):
        sx, px = _fit_pca(data.wsi[train], min(pca_dim, data.wsi.shape[1]))
        sy, py = _fit_pca(data.rna[train], min(pca_dim, data.rna.shape[1]))
        x = _transform_pca(sx, px, data.wsi)
        y = _transform_pca(sy, py, data.rna)
        for pls_dim in (16, 32, 50, 64):
            n = min(pls_dim, x.shape[1], y.shape[1], max(1, int(train.sum()) - 1))
            model = PLSRegression(n_components=n, scale=False, max_iter=1000)
            model.fit(x[train], y[train])
            z_wsi = model.transform(x).astype(np.float32)
            z_rna = model.transform(x, y)[1].astype(np.float32)
            score = _selection_score(data, z_wsi, z_rna, val)
            if best is None or score > best[0]:
                best = (score, z_wsi, z_rna, {"method": "pls_grid", "pca_dim": int(x.shape[1]), "pls_dim": int(n), "selected_by": "val_retrieval_r10_plus_mrr"})
    assert best is not None
    return best[1], best[2], best[3]


def align_procrustes_grid(data: PairedMatrix, dim: int) -> tuple[np.ndarray, np.ndarray, dict]:
    from scipy.linalg import orthogonal_procrustes

    train = _train_mask(data)
    val = _split_mask(data, "val")
    best = None
    for pca_dim in (64, 128, 256, 512):
        sx, px = _fit_pca(data.wsi[train], min(pca_dim, data.wsi.shape[1]))
        sy, py = _fit_pca(data.rna[train], min(pca_dim, data.rna.shape[1]))
        x = _transform_pca(sx, px, data.wsi)
        y = _transform_pca(sy, py, data.rna)
        n = min(x.shape[1], y.shape[1])
        x = x[:, :n]
        y = y[:, :n]
        rotation, scale = orthogonal_procrustes(x[train], y[train])
        z_wsi = (x @ rotation).astype(np.float32)
        z_rna = y.astype(np.float32)
        score = _selection_score(data, z_wsi, z_rna, val)
        if best is None or score > best[0]:
            best = (score, z_wsi, z_rna, {"method": "procrustes_grid", "pca_dim": int(n), "scale": float(scale), "selected_by": "val_retrieval_r10_plus_mrr"})
    assert best is not None
    return best[1], best[2], best[3]


class ResidualProjectionHead:
    def __init__(self, in_dim: int, out_dim: int):
        import torch
        from torch import nn

        self.torch = torch
        self.nn = nn
        self.model = nn.ModuleDict(
            {
                "main": nn.Sequential(
                    nn.Linear(in_dim, out_dim),
                    nn.LayerNorm(out_dim),
                    nn.GELU(),
                    nn.Linear(out_dim, out_dim),
                ),
                "skip": nn.Identity() if in_dim == out_dim else nn.Linear(in_dim, out_dim, bias=False),
                "norm": nn.LayerNorm(out_dim),
            }
        )

    def parameters(self):
        return self.model.parameters()

    def __call__(self, x):
        return self.nn.functional.normalize(self.model["norm"](self.model["main"](x) + self.model["skip"](x)), dim=1)


def align_mlp_regression(data: PairedMatrix, dim: int, epochs: int = 50) -> tuple[np.ndarray, np.ndarray, dict]:
    import torch

    train = _train_mask(data)
    sx, px = _fit_pca(data.wsi[train], min(dim, data.wsi.shape[1], 512))
    sy, py = _fit_pca(data.rna[train], min(dim, data.rna.shape[1], 512))
    x = _transform_pca(sx, px, data.wsi)
    y = _transform_pca(sy, py, data.rna)
    torch.manual_seed(42)
    model = ResidualProjectionHead(x.shape[1], y.shape[1]).model
    opt = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-2)
    xt = torch.tensor(x, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.float32)
    train_t = torch.tensor(train)
    for _ in range(epochs):
        opt.zero_grad()
        pred = model["norm"](model["main"](xt[train_t]) + model["skip"](xt[train_t]))
        loss = torch.nn.functional.mse_loss(pred, yt[train_t])
        loss.backward()
        opt.step()
    with torch.no_grad():
        pred = model["norm"](model["main"](xt) + model["skip"](xt)).numpy().astype(np.float32)
    return pred, y.astype(np.float32), {"method": "mlp_regression", "epochs": epochs, "loss": float(loss.item())}


def align_mlp_clip(data: PairedMatrix, dim: int, epochs: int = 50, hardneg: bool = False) -> tuple[np.ndarray, np.ndarray, dict]:
    import torch

    train = _train_mask(data)
    sx, px = _fit_pca(data.wsi[train], min(512, data.wsi.shape[1]))
    sy, py = _fit_pca(data.rna[train], min(512, data.rna.shape[1]))
    wx = _transform_pca(sx, px, data.wsi)
    rx = _transform_pca(sy, py, data.rna)
    torch.manual_seed(42)
    w_proj = ResidualProjectionHead(wx.shape[1], dim).model
    r_proj = ResidualProjectionHead(rx.shape[1], dim).model
    opt = torch.optim.AdamW(list(w_proj.parameters()) + list(r_proj.parameters()), lr=3e-4, weight_decay=1e-2)
    w = torch.tensor(wx, dtype=torch.float32)
    r = torch.tensor(rx, dtype=torch.float32)
    train_idx = np.where(train)[0]
    if hardneg:
        train_idx = np.asarray(sorted(train_idx, key=lambda i: data.cancers[i]))
    temperature = 0.07
    batch_size = min(256, len(train_idx))
    for _ in range(epochs):
        epoch_loss = 0.0
        for start in range(0, len(train_idx), batch_size):
            idx = train_idx[start : start + batch_size]
            if len(idx) < 2:
                continue
            opt.zero_grad()
            zw = _forward_head(w_proj, w[idx])
            zr = _forward_head(r_proj, r[idx])
            logits = zw @ zr.T / temperature
            labels = torch.arange(len(idx))
            clip_loss = (torch.nn.functional.cross_entropy(logits, labels) + torch.nn.functional.cross_entropy(logits.T, labels)) / 2
            cancer = np.asarray([data.cancers[i] for i in idx])
            pos = torch.tensor(cancer[:, None] == cancer[None, :])
            sup_loss = (_multi_positive_loss(logits, pos) + _multi_positive_loss(logits.T, pos.T)) / 2
            loss = clip_loss + 0.15 * sup_loss
            loss.backward()
            opt.step()
            epoch_loss = float(loss.item())
    with torch.no_grad():
        zw = _forward_head(w_proj, w).numpy().astype(np.float32)
        zr = _forward_head(r_proj, r).numpy().astype(np.float32)
    return zw, zr, {"method": "mlp_clip_hardneg" if hardneg else "mlp_clip", "epochs": epochs, "loss": epoch_loss, "temperature": temperature, "supervised_contrastive_weight": 0.15}


def align_mlp_siglip(data: PairedMatrix, dim: int, epochs: int = 75) -> tuple[np.ndarray, np.ndarray, dict]:
    import torch

    train = _train_mask(data)
    sx, px = _fit_pca(data.wsi[train], min(512, data.wsi.shape[1]))
    sy, py = _fit_pca(data.rna[train], min(512, data.rna.shape[1]))
    wx = _transform_pca(sx, px, data.wsi)
    rx = _transform_pca(sy, py, data.rna)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    torch.manual_seed(42)
    w_proj = ResidualProjectionHead(wx.shape[1], dim).model.to(device)
    r_proj = ResidualProjectionHead(rx.shape[1], dim).model.to(device)
    logit_scale = torch.nn.Parameter(torch.tensor(2.659, device=device))
    bias = torch.nn.Parameter(torch.tensor(-10.0, device=device))
    opt = torch.optim.AdamW(list(w_proj.parameters()) + list(r_proj.parameters()) + [logit_scale, bias], lr=3e-4, weight_decay=1e-2)
    w = torch.tensor(wx, dtype=torch.float32, device=device)
    r = torch.tensor(rx, dtype=torch.float32, device=device)
    train_idx = np.where(train)[0]
    rng = np.random.default_rng(42)
    batch_size = min(256, len(train_idx))
    epoch_loss = 0.0
    for _ in range(epochs):
        rng.shuffle(train_idx)
        for start in range(0, len(train_idx), batch_size):
            idx = train_idx[start : start + batch_size]
            if len(idx) < 2:
                continue
            opt.zero_grad()
            zw = _forward_head(w_proj, w[idx])
            zr = _forward_head(r_proj, r[idx])
            logits = logit_scale.exp().clamp(max=100.0) * (zw @ zr.T) + bias
            labels = torch.eye(len(idx), dtype=torch.float32, device=device)
            loss = torch.nn.functional.binary_cross_entropy_with_logits(logits, labels) + torch.nn.functional.binary_cross_entropy_with_logits(logits.T, labels)
            loss = loss / 2
            loss.backward()
            opt.step()
            epoch_loss = float(loss.item())
    with torch.no_grad():
        zw = _forward_head(w_proj, w).detach().cpu().numpy().astype(np.float32)
        zr = _forward_head(r_proj, r).detach().cpu().numpy().astype(np.float32)
    return zw, zr, {"method": "mlp_siglip", "epochs": epochs, "loss": epoch_loss, "logit_scale": float(logit_scale.detach().cpu()), "bias": float(bias.detach().cpu())}


def align_mlp_debiased(data: PairedMatrix, dim: int, epochs: int = 75) -> tuple[np.ndarray, np.ndarray, dict]:
    import torch

    train = _train_mask(data)
    sx, px = _fit_pca(data.wsi[train], min(512, data.wsi.shape[1]))
    sy, py = _fit_pca(data.rna[train], min(512, data.rna.shape[1]))
    wx = _transform_pca(sx, px, data.wsi)
    rx = _transform_pca(sy, py, data.rna)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    torch.manual_seed(42)
    w_proj = ResidualProjectionHead(wx.shape[1], dim).model.to(device)
    r_proj = ResidualProjectionHead(rx.shape[1], dim).model.to(device)
    opt = torch.optim.AdamW(list(w_proj.parameters()) + list(r_proj.parameters()), lr=3e-4, weight_decay=1e-2)
    w = torch.tensor(wx, dtype=torch.float32, device=device)
    r = torch.tensor(rx, dtype=torch.float32, device=device)
    train_idx = np.where(train)[0]
    cancers = np.asarray(data.cancers)
    rng = np.random.default_rng(42)
    batch_size = min(256, len(train_idx))
    temperature = 0.07
    epoch_loss = 0.0
    for _ in range(epochs):
        rng.shuffle(train_idx)
        for start in range(0, len(train_idx), batch_size):
            idx = train_idx[start : start + batch_size]
            if len(idx) < 2:
                continue
            opt.zero_grad()
            zw = _forward_head(w_proj, w[idx])
            zr = _forward_head(r_proj, r[idx])
            logits = zw @ zr.T / temperature
            same_cancer = torch.tensor(cancers[idx, None] == cancers[None, idx], dtype=torch.bool, device=device)
            identity = torch.eye(len(idx), dtype=torch.bool, device=device)
            positive = same_cancer | identity
            loss = (_multi_positive_loss(logits, positive) + _multi_positive_loss(logits.T, positive.T)) / 2
            loss.backward()
            opt.step()
            epoch_loss = float(loss.item())
    with torch.no_grad():
        zw = _forward_head(w_proj, w).detach().cpu().numpy().astype(np.float32)
        zr = _forward_head(r_proj, r).detach().cpu().numpy().astype(np.float32)
    return zw, zr, {"method": "mlp_debiased", "epochs": epochs, "loss": epoch_loss, "temperature": temperature, "same_cancer_positive": True}


def _forward_head(model, x):
    import torch

    return torch.nn.functional.normalize(model["norm"](model["main"](x) + model["skip"](x)), dim=1)


def _multi_positive_loss(logits, positive_mask):
    import torch

    log_prob = logits - torch.logsumexp(logits, dim=1, keepdim=True)
    pos = positive_mask.to(logits.device).float()
    return -((log_prob * pos).sum(dim=1) / pos.sum(dim=1).clamp_min(1.0)).mean()


def _retrieval_for_mask(data: PairedMatrix, z_wsi: np.ndarray, z_rna: np.ndarray, mask: np.ndarray, k_values: tuple[int, ...]) -> dict:
    if not mask.any():
        return {"n": 0}
    cancers = np.asarray(data.cancers)[mask].tolist()
    metrics = paired_retrieval_metrics(z_wsi[mask], z_rna[mask], k_values, cancers, cancers)
    return {
        "n": int(mask.sum()),
        "retrieval_r1": metrics.get("recall_at_1", 0.0),
        "retrieval_r5": metrics.get("recall_at_5", 0.0),
        "retrieval_r10": metrics.get("recall_at_10", 0.0),
        "retrieval_mrr": metrics.get("mrr", 0.0),
        "same_cancer_at_10": metrics.get("same_cancer_in_top10", 0.0),
        "median_rank": metrics.get("median_rank", 0.0),
        "matched_cosine_mean": metrics.get("matched_cosine_mean", 0.0),
        "unmatched_cosine_mean": metrics.get("unmatched_cosine_mean", 0.0),
    }


def _selection_score(data: PairedMatrix, z_wsi: np.ndarray, z_rna: np.ndarray, mask: np.ndarray) -> float:
    if not mask.any():
        mask = _train_mask(data)
    cancers = np.asarray(data.cancers)[mask].tolist()
    metrics = paired_retrieval_metrics(z_wsi[mask], z_rna[mask], (10,), cancers, cancers)
    return float(metrics.get("recall_at_10", 0.0)) + float(metrics.get("mrr", 0.0))


def run_alignment(config_path: str, method: str, split_file: str | Path | None = None, output_dir: str | Path | None = None) -> Path:
    cfg = load_config(config_path)
    data = load_paired_matrix(config_path, split_file)
    dim = int(cfg.section("alignment").get("embedding_dim", 256))
    epochs = int(cfg.section("alignment").get("max_epochs", 50))
    if method == "none":
        z_wsi, z_rna, details = align_none(data)
    elif method == "ridge":
        z_wsi, z_rna, details = align_ridge(data, dim)
    elif method == "ridge_grid":
        z_wsi, z_rna, details = align_ridge_grid(data, dim)
    elif method == "cca":
        z_wsi, z_rna, details = align_cca(data, dim)
    elif method == "cca_grid":
        z_wsi, z_rna, details = align_cca_grid(data, dim)
    elif method == "pls_grid":
        z_wsi, z_rna, details = align_pls_grid(data, dim)
    elif method == "procrustes_grid":
        z_wsi, z_rna, details = align_procrustes_grid(data, dim)
    elif method == "mlp_regression":
        z_wsi, z_rna, details = align_mlp_regression(data, dim, epochs)
    elif method == "mlp_clip":
        z_wsi, z_rna, details = align_mlp_clip(data, dim, epochs, False)
    elif method == "mlp_clip_hardneg":
        z_wsi, z_rna, details = align_mlp_clip(data, dim, epochs, True)
    elif method == "mlp_siglip":
        z_wsi, z_rna, details = align_mlp_siglip(data, dim, max(epochs, 75))
    elif method == "mlp_debiased":
        z_wsi, z_rna, details = align_mlp_debiased(data, dim, max(epochs, 75))
    else:
        raise ValueError(f"Unknown method: {method}")
    k_values = tuple(cfg.section("alignment").get("k_values", [1, 5, 10]))
    masks = {name: _split_mask(data, name) for name in ("train", "val", "test")}
    metrics_by_split = {name: _retrieval_for_mask(data, z_wsi, z_rna, mask, k_values) for name, mask in masks.items()}
    out_dir = Path(output_dir) if output_dir else cfg.path("outputs_dir") / "v1_alignment"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_dir / f"{method}_aligned_embeddings.npz", patient_ids=np.asarray(data.patient_ids), split=np.asarray(data.split), cancers=np.asarray(data.cancers), wsi=z_wsi, rna=z_rna)
    payload = base_manifest(cfg.project_root, cfg.config_path, int(cfg.raw.get("seed", 42)))
    rna_source = "bulkformer" if cfg.path("rna_bulkformer_embeddings").exists() else "fallback_geneformer_or_processed_rna"
    payload.update(
        {
            "method": method,
            "details": details,
            "split_file": str(split_file) if split_file is not None else str(cfg.path("processed_dir") / "splits" / "stratified_pan_cancer_split.json"),
            "metrics": metrics_by_split["test"],
            "train_metrics": metrics_by_split["train"],
            "val_metrics": metrics_by_split["val"],
            "test_metrics": metrics_by_split["test"],
            "n_pairs": len(data.patient_ids),
            "n_train": int(masks["train"].sum()),
            "n_val": int(masks["val"].sum()),
            "n_test": int(masks["test"].sum()),
            "rna_source": rna_source,
        }
    )
    out_path = out_dir / f"{method}_metrics.json"
    write_json(out_path, payload)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument(
        "--method",
        choices=["none", "ridge", "ridge_grid", "cca", "cca_grid", "pls_grid", "procrustes_grid", "mlp_regression", "mlp_clip", "mlp_clip_hardneg", "mlp_siglip", "mlp_debiased"],
        required=True,
    )
    parser.add_argument("--split-file")
    parser.add_argument("--output-dir")
    args = parser.parse_args()
    print(run_alignment(args.config, args.method, args.split_file, args.output_dir))


if __name__ == "__main__":
    main()
