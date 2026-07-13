"""Train the V1 Tumor-State Query Former on frozen foundation features."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from morpheus.src.data.tumor_state_registry import build_tumor_state_registry, create_tumor_state_splits
from morpheus.src.eval.retrieval_metrics import paired_retrieval_metrics
from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import base_manifest, write_json


@dataclass
class QueryFormerData:
    patient_ids: list[str]
    cancers: list[str]
    split: np.ndarray
    wsi: np.ndarray
    rna: np.ndarray
    hallmark: np.ndarray
    hallmark_present: np.ndarray
    clinical: np.ndarray
    clinical_present: np.ndarray
    hallmark_names: list[str]
    clinical_names: list[str]


def _variant_loss_weights(variant: str) -> dict[str, float]:
    variants = {
        "v1": {
            "clip": 1.0,
            "hallmark_fused": 0.5,
            "hallmark_wsi": 0.0,
            "hallmark_rna": 0.0,
            "rna_recon": 0.2,
            "wsi_recon": 0.0,
            "neighborhood": 0.0,
            "teacher": 0.0,
            "vicreg": 0.0,
        },
        "wsi_hallmark_loss": {
            "clip": 1.0,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 1.0,
            "hallmark_rna": 0.0,
            "rna_recon": 0.2,
            "wsi_recon": 0.0,
            "neighborhood": 0.0,
            "teacher": 0.0,
            "vicreg": 0.0,
        },
        "wsi_rna_hallmark_loss": {
            "clip": 1.0,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 1.0,
            "hallmark_rna": 0.5,
            "rna_recon": 0.2,
            "wsi_recon": 0.0,
            "neighborhood": 0.0,
            "teacher": 0.0,
            "vicreg": 0.0,
        },
        "neighborhood_distill": {
            "clip": 0.9,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 1.0,
            "hallmark_rna": 0.5,
            "rna_recon": 0.15,
            "wsi_recon": 0.0,
            "neighborhood": 0.35,
            "teacher": 0.0,
            "vicreg": 0.01,
        },
        "teacher_distill": {
            "clip": 0.9,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 1.0,
            "hallmark_rna": 0.5,
            "rna_recon": 0.15,
            "wsi_recon": 0.0,
            "neighborhood": 0.0,
            "teacher": 0.1,
            "vicreg": 0.01,
        },
        "combined_best": {
            "clip": 0.9,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 1.0,
            "hallmark_rna": 0.5,
            "rna_recon": 0.15,
            "wsi_recon": 0.0,
            "neighborhood": 0.35,
            "teacher": 0.08,
            "vicreg": 0.02,
        },
        "retentive_hallmark": {
            "clip": 1.0,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 1.0,
            "hallmark_rna": 0.25,
            "rna_recon": 0.15,
            "wsi_recon": 0.25,
            "neighborhood": 0.15,
            "teacher": 0.0,
            "vicreg": 0.01,
        },
        "strong_teacher_hallmark": {
            "clip": 1.0,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 0.75,
            "hallmark_rna": 0.25,
            "rna_recon": 0.1,
            "wsi_recon": 0.0,
            "neighborhood": 0.0,
            "teacher": 20.0,
            "vicreg": 0.01,
        },
        "teacher_hallmark_x1": {
            "clip": 1.0,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 0.75,
            "hallmark_rna": 0.25,
            "rna_recon": 0.1,
            "wsi_recon": 0.0,
            "neighborhood": 0.0,
            "teacher": 1.0,
            "vicreg": 0.01,
        },
        "teacher_hallmark_x5": {
            "clip": 1.0,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 0.75,
            "hallmark_rna": 0.25,
            "rna_recon": 0.1,
            "wsi_recon": 0.0,
            "neighborhood": 0.0,
            "teacher": 5.0,
            "vicreg": 0.01,
        },
        "wsi_hallmark_strong": {
            "clip": 0.75,
            "hallmark_fused": 0.25,
            "hallmark_wsi": 5.0,
            "hallmark_rna": 0.25,
            "rna_recon": 0.05,
            "wsi_recon": 0.0,
            "neighborhood": 0.15,
            "teacher": 0.0,
            "vicreg": 0.02,
        },
    }
    if variant not in variants:
        raise ValueError(f"Unknown QueryFormer variant: {variant}")
    return variants[variant]


def _load_wsi(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        arr = handle["embeddings"][:].astype(np.float32)
        ids = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in handle["patient_ids"][:]]
    return pd.DataFrame({"patient_id": ids, "wsi_vector": [row for row in arr]})


def _load_rna(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    numeric = [c for c in df.select_dtypes(include=["number"]).columns if c != "patient_id"]
    rows = []
    for patient_id, group in df.groupby("patient_id", dropna=True):
        rows.append({"patient_id": str(patient_id), "rna_vector": group[numeric].to_numpy(dtype=np.float32).mean(axis=0)})
    return pd.DataFrame(rows)


def _load_hallmark(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_parquet(path)
    names = [c for c in df.columns if c != "patient_id" and pd.api.types.is_numeric_dtype(df[c])]
    return df[["patient_id", *names]].copy(), names


def _load_clinical(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_parquet(path)
    if "patient_id" not in df.columns:
        df = df.reset_index().rename(columns={df.index.name or "index": "patient_id"})
    names = [c for c in df.columns if c != "patient_id" and pd.api.types.is_numeric_dtype(df[c]) or c != "patient_id" and pd.api.types.is_bool_dtype(df[c])]
    df[names] = df[names].astype(np.float32)
    return df[["patient_id", *names]].copy(), names


def _split_labels(split_path: Path, patient_ids: list[str]) -> np.ndarray:
    payload = json.load(open(split_path, encoding="utf-8"))
    labels = {pid: "excluded" for pid in patient_ids}
    for name in ("train", "val", "test"):
        for pid in payload["patient_ids"].get(name, []):
            labels[str(pid)] = name
    return np.asarray([labels[pid] for pid in patient_ids])


def _scale_trainfit(matrix: np.ndarray, train: np.ndarray) -> tuple[np.ndarray, dict]:
    mean = matrix[train].mean(axis=0, keepdims=True)
    std = matrix[train].std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return ((matrix - mean) / std).astype(np.float32), {"mean": mean.squeeze(0).tolist(), "std": std.squeeze(0).tolist()}


def load_query_former_data(config_path: str, split_file: str | Path | None = None) -> QueryFormerData:
    cfg = load_config(config_path)
    registry_path = cfg.path("processed_dir") / "patient_feature_registry.parquet"
    if not registry_path.exists():
        build_tumor_state_registry(config_path)
        create_tumor_state_splits(config_path)
    registry = pd.read_parquet(registry_path)
    wsi = _load_wsi(cfg.path("wsi_standard_dir") / "tcga_ut_hoptimus0_patient_embeddings.h5")
    rna = _load_rna(cfg.path("rna_bulkformer_embeddings"))
    hallmark, hallmark_names = _load_hallmark(cfg.path("hallmark_scores"))
    clinical_path = cfg.project_root / "meta-intersurv" / "data" / "embeddings" / "patient_master" / "tabular_features.parquet"
    clinical, clinical_names = _load_clinical(clinical_path)
    frame = registry[["patient_id", "cancer_type"]].dropna(subset=["cancer_type"]).merge(wsi, on="patient_id", how="inner").merge(rna, on="patient_id", how="inner")
    frame = frame.merge(hallmark, on="patient_id", how="left").merge(clinical, on="patient_id", how="left")
    frame = frame.sort_values("patient_id").reset_index(drop=True)
    split_path = Path(split_file) if split_file else cfg.path("processed_dir") / "splits" / "tumor_state_stratified.json"
    if not split_path.exists():
        create_tumor_state_splits(config_path)
    split = _split_labels(split_path, frame["patient_id"].astype(str).tolist())
    clinical_values = frame[clinical_names].to_numpy(dtype=np.float32)
    clinical_present = ~np.isnan(clinical_values).any(axis=1)
    clinical_values = np.nan_to_num(clinical_values, nan=0.0).astype(np.float32)
    hallmark_values = frame[hallmark_names].to_numpy(dtype=np.float32)
    hallmark_present = ~np.isnan(hallmark_values).any(axis=1)
    hallmark_values = np.nan_to_num(hallmark_values, nan=0.0).astype(np.float32)
    train = split == "train"
    wsi_arr, _ = _scale_trainfit(np.vstack(frame["wsi_vector"].to_numpy()).astype(np.float32), train)
    rna_arr, _ = _scale_trainfit(np.vstack(frame["rna_vector"].to_numpy()).astype(np.float32), train)
    clinical_values, _ = _scale_trainfit(clinical_values, train & clinical_present) if clinical_present.any() else (clinical_values, {})
    hallmark_values, _ = _scale_trainfit(hallmark_values, train & hallmark_present) if hallmark_present.any() else (hallmark_values, {})
    return QueryFormerData(
        patient_ids=frame["patient_id"].astype(str).tolist(),
        cancers=frame["cancer_type"].astype(str).tolist(),
        split=split,
        wsi=wsi_arr,
        rna=rna_arr,
        hallmark=hallmark_values,
        hallmark_present=hallmark_present,
        clinical=clinical_values,
        clinical_present=clinical_present,
        hallmark_names=hallmark_names,
        clinical_names=clinical_names,
    )


def _clip_loss(a, b, temperature: float = 0.07):
    import torch

    a = torch.nn.functional.normalize(a, dim=1)
    b = torch.nn.functional.normalize(b, dim=1)
    logits = a @ b.T / temperature
    labels = torch.arange(a.shape[0], device=a.device)
    return (torch.nn.functional.cross_entropy(logits, labels) + torch.nn.functional.cross_entropy(logits.T, labels)) / 2


def _hallmark_neighborhood_loss(z, hallmark, temperature: float = 0.15):
    import torch

    if z.shape[0] < 3:
        return z.new_tensor(0.0)
    z = torch.nn.functional.normalize(z, dim=1)
    hallmark = torch.nn.functional.normalize(hallmark, dim=1)
    logits = z @ z.T / temperature
    target_logits = hallmark @ hallmark.T / temperature
    eye = torch.eye(z.shape[0], dtype=torch.bool, device=z.device)
    logits = logits.masked_fill(eye, -1e4)
    target_logits = target_logits.masked_fill(eye, -1e4)
    target_prob = torch.softmax(target_logits, dim=1)
    log_prob = torch.log_softmax(logits, dim=1)
    return torch.nn.functional.kl_div(log_prob, target_prob, reduction="batchmean")


def _vicreg_loss(z, eps: float = 1e-4):
    import torch

    if z.shape[0] < 2:
        return z.new_tensor(0.0)
    z = z - z.mean(dim=0, keepdim=True)
    std = torch.sqrt(z.var(dim=0) + eps)
    variance = torch.relu(1.0 - std).mean()
    cov = (z.T @ z) / max(z.shape[0] - 1, 1)
    off_diag = cov.flatten()[:-1].view(z.shape[1] - 1, z.shape[1] + 1)[:, 1:].flatten()
    covariance = off_diag.pow(2).mean()
    return variance + covariance


def _mean_gene_pearson(truth: np.ndarray, pred: np.ndarray) -> tuple[float | None, int]:
    vals = []
    for j in range(truth.shape[1]):
        if np.std(truth[:, j]) > 0 and np.std(pred[:, j]) > 0:
            vals.append(float(np.corrcoef(truth[:, j], pred[:, j])[0, 1]))
    return (float(np.nanmean(vals)) if vals else None, len(vals))


def _load_teacher_vectors(teacher_npz: str | Path, data: QueryFormerData, key: str = "wsi") -> tuple[np.ndarray, np.ndarray]:
    teacher = np.load(teacher_npz, allow_pickle=True)
    if key not in teacher.files:
        raise KeyError(f"{key!r} not found in teacher embeddings {teacher_npz}; available keys: {teacher.files}")
    teacher_ids = teacher["patient_ids"].astype(str)
    values = teacher[key].astype(np.float32)
    by_id = {pid: i for i, pid in enumerate(teacher_ids)}
    dim = values.shape[1]
    aligned = np.zeros((len(data.patient_ids), dim), dtype=np.float32)
    present = np.zeros(len(data.patient_ids), dtype=bool)
    for i, pid in enumerate(data.patient_ids):
        row = by_id.get(pid)
        if row is not None:
            aligned[i] = values[row]
            present[i] = True
    train_present = (data.split == "train") & present
    if train_present.any():
        aligned, _ = _scale_trainfit(aligned, train_present)
    return aligned.astype(np.float32), present


def _batch_iter(indices: np.ndarray, batch_size: int, rng: np.random.Generator):
    order = indices.copy()
    rng.shuffle(order)
    for start in range(0, len(order), batch_size):
        yield order[start : start + batch_size]


def _model_outputs(model, adapters, batch, device, use_wsi: bool = True, use_rna: bool = True, use_clinical: bool = True):
    import torch

    wsi = torch.tensor(batch["wsi"], dtype=torch.float32, device=device)
    rna = torch.tensor(batch["rna"], dtype=torch.float32, device=device)
    clinical = torch.tensor(batch["clinical"], dtype=torch.float32, device=device)
    clinical_present = torch.tensor(batch["clinical_present"], dtype=torch.bool, device=device)
    wsi_present = torch.ones(wsi.shape[0], dtype=torch.bool, device=device) if use_wsi else torch.zeros(wsi.shape[0], dtype=torch.bool, device=device)
    rna_present = torch.ones(rna.shape[0], dtype=torch.bool, device=device) if use_rna else torch.zeros(rna.shape[0], dtype=torch.bool, device=device)
    clinical_present = clinical_present & use_clinical
    modalities = {
        "wsi": adapters["wsi"](wsi, wsi_present),
        "rna": adapters["rna"](rna, None, bulk_present=rna_present),
        "clinical": adapters["clinical"](clinical, clinical_present),
    }
    return model(modalities)


def _pool(x):
    return x.mean(dim=1)


def run_query_former_training(
    config_path: str = "morpheus/configs/v1.json",
    split_file: str | None = None,
    output_dir: str | Path = "morpheus/outputs/v1_query_former",
    epochs: int = 30,
    batch_size: int = 64,
    smoke: bool = False,
    device_name: str = "auto",
    variant: str = "v1",
    teacher_npz: str | None = None,
    teacher_key: str = "wsi",
) -> Path:
    import torch
    from torch import nn

    from morpheus.src.encoders.clinical_adapter import ClinicalTokenAdapter
    from morpheus.src.encoders.rna_adapter import RNATokenAdapter
    from morpheus.src.encoders.wsi_adapter import WSITokenAdapter
    from morpheus.src.models.tumor_state_query_former import QueryFormerConfig, TumorStateQueryFormer

    cfg = load_config(config_path)
    seed = int(cfg.raw.get("seed", 42))
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    data = load_query_former_data(config_path, split_file)
    teacher_vectors = None
    teacher_present = None
    if teacher_npz:
        teacher_vectors, teacher_present = _load_teacher_vectors(teacher_npz, data, teacher_key)
    if smoke:
        keep = np.r_[np.where(data.split == "train")[0][:128], np.where(data.split == "val")[0][:64], np.where(data.split == "test")[0][:64]]
        data = QueryFormerData(
            [data.patient_ids[i] for i in keep],
            [data.cancers[i] for i in keep],
            data.split[keep],
            data.wsi[keep],
            data.rna[keep],
            data.hallmark[keep],
            data.hallmark_present[keep],
            data.clinical[keep],
            data.clinical_present[keep],
            data.hallmark_names,
            data.clinical_names,
        )
        if teacher_vectors is not None and teacher_present is not None:
            teacher_vectors = teacher_vectors[keep]
            teacher_present = teacher_present[keep]
    device = "cuda" if device_name == "auto" and torch.cuda.is_available() else ("cpu" if device_name == "auto" else device_name)
    qcfg = QueryFormerConfig(hidden_dim=512, num_layers=2, num_heads=8, dropout=0.1)
    model = TumorStateQueryFormer(qcfg).to(device)
    adapters = nn.ModuleDict(
        {
            "wsi": WSITokenAdapter(data.wsi.shape[1], qcfg.hidden_dim),
            "rna": RNATokenAdapter(data.rna.shape[1], data.hallmark.shape[1], qcfg.hidden_dim),
            "clinical": ClinicalTokenAdapter(data.clinical.shape[1], qcfg.hidden_dim),
        }
    ).to(device)
    heads = nn.ModuleDict(
        {
            "hallmark": nn.Linear(qcfg.hidden_dim, data.hallmark.shape[1]),
            "rna_recon": nn.Linear(qcfg.hidden_dim, data.rna.shape[1]),
        }
    ).to(device)
    if teacher_vectors is not None:
        heads["teacher"] = nn.Linear(qcfg.hidden_dim, teacher_vectors.shape[1]).to(device)
    if _variant_loss_weights(variant)["wsi_recon"] > 0:
        heads["wsi_recon"] = nn.Linear(qcfg.hidden_dim, data.wsi.shape[1]).to(device)
    params = list(model.parameters()) + list(adapters.parameters()) + list(heads.parameters())
    opt = torch.optim.AdamW(params, lr=1e-4, weight_decay=1e-2)
    loss_weights = _variant_loss_weights(variant)
    train_idx = np.where(data.split == "train")[0]
    val_idx = np.where(data.split == "val")[0]
    test_idx = np.where(data.split == "test")[0]
    history = []
    best_val = -np.inf
    best_state = None
    patience = 10
    stale = 0

    for epoch in range(epochs):
        model.train()
        adapters.train()
        heads.train()
        losses = []
        for idx in _batch_iter(train_idx, batch_size, rng):
            if len(idx) < 2:
                continue
            batch = {k: getattr(data, k)[idx] for k in ("wsi", "rna", "hallmark", "hallmark_present", "clinical", "clinical_present")}
            if teacher_vectors is not None and teacher_present is not None:
                batch["teacher"] = teacher_vectors[idx]
                batch["teacher_present"] = teacher_present[idx]
            opt.zero_grad()
            out_full = _model_outputs(model, adapters, batch, device)
            out_wsi = _model_outputs(model, adapters, batch, device, use_rna=False, use_clinical=False)
            out_rna = _model_outputs(model, adapters, batch, device, use_wsi=False, use_clinical=False)
            z_shared = out_full["z_patient"]
            z_wsi = out_wsi["z_patient"]
            z_rna = out_rna["z_patient"]
            hallmark_target = torch.tensor(batch["hallmark"], dtype=torch.float32, device=device)
            hallmark_present = torch.tensor(batch["hallmark_present"], dtype=torch.bool, device=device)
            rna_target = torch.tensor(batch["rna"], dtype=torch.float32, device=device)
            wsi_target = torch.tensor(batch["wsi"], dtype=torch.float32, device=device)
            loss_terms = {}
            loss = loss_weights["clip"] * _clip_loss(z_wsi, z_rna)
            loss_terms["clip"] = float(loss.detach().cpu().item())
            if hallmark_present.any():
                hallmark_loss = torch.nn.functional.mse_loss(heads["hallmark"](z_shared)[hallmark_present], hallmark_target[hallmark_present])
                wsi_hallmark_loss = torch.nn.functional.mse_loss(heads["hallmark"](z_wsi)[hallmark_present], hallmark_target[hallmark_present])
                rna_hallmark_loss = torch.nn.functional.mse_loss(heads["hallmark"](z_rna)[hallmark_present], hallmark_target[hallmark_present])
                loss = loss + loss_weights["hallmark_fused"] * hallmark_loss
                loss = loss + loss_weights["hallmark_wsi"] * wsi_hallmark_loss
                loss = loss + loss_weights["hallmark_rna"] * rna_hallmark_loss
                if loss_weights["neighborhood"] > 0:
                    loss = loss + loss_weights["neighborhood"] * _hallmark_neighborhood_loss(z_wsi[hallmark_present], hallmark_target[hallmark_present])
            loss = loss + loss_weights["rna_recon"] * torch.nn.functional.mse_loss(heads["rna_recon"](out_wsi["z_patient"]), rna_target)
            if "wsi_recon" in heads and loss_weights["wsi_recon"] > 0:
                loss = loss + loss_weights["wsi_recon"] * torch.nn.functional.mse_loss(heads["wsi_recon"](out_wsi["z_patient"]), wsi_target)
            if teacher_vectors is not None and "teacher" in heads and loss_weights["teacher"] > 0:
                teacher_target = torch.tensor(batch["teacher"], dtype=torch.float32, device=device)
                teacher_mask = torch.tensor(batch["teacher_present"], dtype=torch.bool, device=device)
                if teacher_mask.any():
                    teacher_pred = heads["teacher"](z_wsi)[teacher_mask]
                    teacher_target = teacher_target[teacher_mask]
                    teacher_loss = torch.nn.functional.mse_loss(teacher_pred, teacher_target)
                    teacher_cosine = 1.0 - torch.nn.functional.cosine_similarity(teacher_pred, teacher_target, dim=1).mean()
                    loss = loss + loss_weights["teacher"] * (teacher_loss + teacher_cosine)
            if loss_weights["vicreg"] > 0:
                loss = loss + loss_weights["vicreg"] * (_vicreg_loss(z_wsi) + _vicreg_loss(z_rna))
            loss.backward()
            opt.step()
            losses.append(float(loss.item()))
        metrics = _evaluate_query_former(data, model, adapters, heads, val_idx, device)
        row = {"epoch": epoch, "train_loss": float(np.mean(losses)) if losses else float("nan"), **{f"val_{k}": v for k, v in metrics.items() if isinstance(v, (float, int))}}
        history.append(row)
        hallmark_score = metrics.get("hallmark_wsi_pearson") if variant != "v1" else metrics.get("hallmark_pearson")
        val_score = float(metrics.get("retrieval_r10", 0.0)) + (0.0 if hallmark_score is None else float(hallmark_score))
        if val_score > best_val:
            best_val = val_score
            best_state = {
                "model": model.state_dict(),
                "adapters": adapters.state_dict(),
                "heads": heads.state_dict(),
                "config": qcfg.__dict__,
                "variant": variant,
                "loss_weights": loss_weights,
                "teacher_dim": None if teacher_vectors is None else int(teacher_vectors.shape[1]),
                "wsi_recon_dim": int(data.wsi.shape[1]) if "wsi_recon" in heads else None,
                "teacher_key": teacher_key,
            }
            stale = 0
        else:
            stale += 1
        if stale >= patience:
            break
    out_dir = Path(output_dir)
    (out_dir / "checkpoints").mkdir(parents=True, exist_ok=True)
    if best_state:
        model.load_state_dict(best_state["model"])
        adapters.load_state_dict(best_state["adapters"])
        heads.load_state_dict(best_state["heads"])
        torch.save(best_state, out_dir / "checkpoints" / "best.pt")
    pd.DataFrame(history).to_csv(out_dir / "train_log.csv", index=False)
    val_metrics = _evaluate_query_former(data, model, adapters, heads, val_idx, device)
    test_metrics = _evaluate_query_former(data, model, adapters, heads, test_idx, device)
    _export_embeddings(data, model, adapters, out_dir, device)
    _export_alignment_views(data, model, adapters, out_dir, device, heads)
    payload = base_manifest(cfg.project_root, cfg.config_path, seed)
    payload.update({"model": "TumorStateQueryFormer", "variant": variant, "loss_weights": loss_weights, "teacher_npz": teacher_npz, "teacher_key": teacher_key, "modality_specific_clinical": False, "device": device, "epochs_run": len(history), "n_train": int(len(train_idx)), "n_val": int(len(val_idx)), "n_test": int(len(test_idx)), "val_metrics": val_metrics, "test_metrics": test_metrics, "checkpoint": str(out_dir / "checkpoints" / "best.pt")})
    write_json(out_dir / "val_metrics.json", {**payload, "metrics": val_metrics})
    write_json(out_dir / "test_metrics.json", {**payload, "metrics": test_metrics})
    return out_dir / "test_metrics.json"


def _build_query_former_modules(data: QueryFormerData, qcfg, teacher_dim: int | None = None, wsi_recon_dim: int | None = None):
    from torch import nn

    from morpheus.src.encoders.clinical_adapter import ClinicalTokenAdapter
    from morpheus.src.encoders.rna_adapter import RNATokenAdapter
    from morpheus.src.encoders.wsi_adapter import WSITokenAdapter
    from morpheus.src.models.tumor_state_query_former import TumorStateQueryFormer

    model = TumorStateQueryFormer(qcfg)
    adapters = nn.ModuleDict(
        {
            "wsi": WSITokenAdapter(data.wsi.shape[1], qcfg.hidden_dim),
            "rna": RNATokenAdapter(data.rna.shape[1], data.hallmark.shape[1], qcfg.hidden_dim),
            "clinical": ClinicalTokenAdapter(data.clinical.shape[1], qcfg.hidden_dim),
        }
    )
    heads = nn.ModuleDict(
        {
            "hallmark": nn.Linear(qcfg.hidden_dim, data.hallmark.shape[1]),
            "rna_recon": nn.Linear(qcfg.hidden_dim, data.rna.shape[1]),
        }
    )
    if teacher_dim is not None:
        heads["teacher"] = nn.Linear(qcfg.hidden_dim, teacher_dim)
    if wsi_recon_dim is not None:
        heads["wsi_recon"] = nn.Linear(qcfg.hidden_dim, wsi_recon_dim)
    return model, adapters, heads


def _checkpoint_teacher_dim(state: dict) -> int | None:
    teacher_dim = state.get("teacher_dim")
    if teacher_dim is not None:
        return int(teacher_dim)
    heads = state.get("heads", {})
    weight = heads.get("teacher.weight")
    return None if weight is None else int(weight.shape[0])


def _checkpoint_wsi_recon_dim(state: dict) -> int | None:
    heads = state.get("heads", {})
    weight = heads.get("wsi_recon.weight")
    return None if weight is None else int(weight.shape[0])


def evaluate_query_former_checkpoint(
    checkpoint: str | Path,
    config_path: str = "morpheus/configs/v1.json",
    split_file: str | None = None,
    output_dir: str | Path | None = None,
    device_name: str = "auto",
) -> Path:
    import torch

    from morpheus.src.models.tumor_state_query_former import QueryFormerConfig

    cfg = load_config(config_path)
    seed = int(cfg.raw.get("seed", 42))
    data = load_query_former_data(config_path, split_file)
    device = "cuda" if device_name == "auto" and torch.cuda.is_available() else ("cpu" if device_name == "auto" else device_name)
    state = torch.load(checkpoint, map_location=device)
    qcfg = QueryFormerConfig(**state.get("config", {}))
    model, adapters, heads = _build_query_former_modules(data, qcfg, _checkpoint_teacher_dim(state), _checkpoint_wsi_recon_dim(state))
    model.load_state_dict(state["model"])
    adapters.load_state_dict(state["adapters"])
    heads.load_state_dict(state["heads"])
    model.to(device)
    adapters.to(device)
    heads.to(device)
    metrics_by_split = {}
    for name in ("train", "val", "test"):
        idx = np.where(data.split == name)[0]
        metrics_by_split[name] = _evaluate_query_former(data, model, adapters, heads, idx, device)
    out_dir = Path(output_dir) if output_dir else Path(checkpoint).parents[1]
    _export_alignment_views(data, model, adapters, out_dir, device, heads)
    payload = base_manifest(cfg.project_root, cfg.config_path, seed)
    payload.update(
        {
            "model": "TumorStateQueryFormer",
            "variant": state.get("variant", "v1"),
            "modality_specific_clinical": False,
            "checkpoint": str(checkpoint),
            "split_file": str(split_file) if split_file else str(cfg.path("processed_dir") / "splits" / "tumor_state_stratified.json"),
            "device": device,
            "n_train": int(np.sum(data.split == "train")),
            "n_val": int(np.sum(data.split == "val")),
            "n_test": int(np.sum(data.split == "test")),
            "train_metrics": metrics_by_split["train"],
            "val_metrics": metrics_by_split["val"],
            "test_metrics": metrics_by_split["test"],
        }
    )
    out_path = out_dir / "eval_metrics.json"
    write_json(out_path, payload)
    return out_path


def _export_alignment_views(data: QueryFormerData, model, adapters, out_dir: Path, device: str, heads=None) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    wsi_arrays = _encode_all(data, model, adapters, device, use_rna=False, use_clinical=False)
    rna_arrays = _encode_all(data, model, adapters, device, use_wsi=False, use_clinical=False)
    full_arrays = _encode_all(data, model, adapters, device)
    out_path = out_dir / "query_former_aligned_embeddings.npz"
    payload = {
        "patient_ids": np.asarray(data.patient_ids),
        "split": np.asarray(data.split),
        "cancers": np.asarray(data.cancers),
        "wsi": wsi_arrays["z_patient"],
        "rna": rna_arrays["z_patient"],
        "full": full_arrays["z_patient"],
        "wsi_task": wsi_arrays["z_task"],
        "rna_task": rna_arrays["z_task"],
        "full_task": full_arrays["z_task"],
    }
    if heads is not None and "wsi_recon" in heads:
        payload["wsi_recon"] = _apply_head_all(heads["wsi_recon"], wsi_arrays["z_patient"], device)
    if heads is not None and "teacher" in heads:
        wsi_teacher = _apply_head_all(heads["teacher"], wsi_arrays["z_patient"], device)
        payload["wsi_teacher"] = wsi_teacher
        payload["wsi_bio"] = np.concatenate([wsi_arrays["z_patient"], wsi_teacher], axis=1).astype(np.float32)
    np.savez_compressed(out_path, **payload)
    return out_path


def _apply_head_all(head, x: np.ndarray, device: str, batch_size: int = 512) -> np.ndarray:
    import torch

    head.eval()
    out = []
    with torch.no_grad():
        for start in range(0, x.shape[0], batch_size):
            stop = min(start + batch_size, x.shape[0])
            batch = torch.tensor(x[start:stop], dtype=torch.float32, device=device)
            out.append(head(batch).detach().cpu().numpy().astype(np.float32))
    return np.vstack(out)


def _encode_all(data: QueryFormerData, model, adapters, device: str, batch_size: int = 256, use_wsi: bool = True, use_rna: bool = True, use_clinical: bool = True) -> dict[str, np.ndarray]:
    model.eval()
    adapters.eval()
    arrays: dict[str, list[np.ndarray]] = {"z_patient": [], "z_shared": [], "z_wsi_resid": [], "z_rna_resid": [], "z_clinical_resid": [], "z_uncertainty": [], "z_task": []}
    with __import__("torch").no_grad():
        for start in range(0, len(data.patient_ids), batch_size):
            stop = min(start + batch_size, len(data.patient_ids))
            idx = np.arange(start, stop)
            batch = {k: getattr(data, k)[idx] for k in ("wsi", "rna", "hallmark", "hallmark_present", "clinical", "clinical_present")}
            out = _model_outputs(model, adapters, batch, device, use_wsi=use_wsi, use_rna=use_rna, use_clinical=use_clinical)
            for key in arrays:
                val = out[key].mean(dim=1) if out[key].ndim == 3 else out[key]
                arrays[key].append(val.detach().cpu().numpy().astype(np.float32))
    return {k: np.vstack(v) for k, v in arrays.items()}


def _evaluate_query_former(data: QueryFormerData, model, adapters, heads, idx: np.ndarray, device: str) -> dict:
    if len(idx) < 2:
        return {"n": int(len(idx))}
    wsi_arrays = _encode_all(data, model, adapters, device, use_rna=False, use_clinical=False)
    rna_arrays = _encode_all(data, model, adapters, device, use_wsi=False, use_clinical=False)
    full_arrays = _encode_all(data, model, adapters, device)
    z_wsi = wsi_arrays["z_patient"][idx]
    z_rna = rna_arrays["z_patient"][idx]
    retrieval = paired_retrieval_metrics(z_wsi, z_rna, (1, 5, 10), [data.cancers[i] for i in idx], [data.cancers[i] for i in idx])
    present = data.hallmark_present[idx]
    pearson = None
    hallmark_gene_count = 0
    if present.any():
        import torch

        with torch.no_grad():
            full_pred = heads["hallmark"](torch.tensor(full_arrays["z_patient"][idx][present], dtype=torch.float32, device=device)).detach().cpu().numpy()
            wsi_pred = heads["hallmark"](torch.tensor(wsi_arrays["z_patient"][idx][present], dtype=torch.float32, device=device)).detach().cpu().numpy()
            rna_pred = heads["hallmark"](torch.tensor(rna_arrays["z_patient"][idx][present], dtype=torch.float32, device=device)).detach().cpu().numpy()
        truth = data.hallmark[idx][present]
        pearson, hallmark_gene_count = _mean_gene_pearson(truth, full_pred)
        wsi_pearson, wsi_gene_count = _mean_gene_pearson(truth, wsi_pred)
        rna_pearson, rna_gene_count = _mean_gene_pearson(truth, rna_pred)
    else:
        wsi_pearson = None
        rna_pearson = None
        wsi_gene_count = 0
        rna_gene_count = 0
    return {
        "n": int(len(idx)),
        "retrieval_r1": retrieval["recall_at_1"],
        "retrieval_r5": retrieval["recall_at_5"],
        "retrieval_r10": retrieval["recall_at_10"],
        "retrieval_mrr": retrieval["mrr"],
        "same_cancer_at_10": retrieval.get("same_cancer_in_top10", 0.0),
        "hallmark_pearson": pearson,
        "hallmark_wsi_pearson": wsi_pearson,
        "hallmark_rna_pearson": rna_pearson,
        "hallmark_n": int(present.sum()),
        "hallmark_gene_count": int(hallmark_gene_count),
        "hallmark_wsi_gene_count": int(wsi_gene_count),
        "hallmark_rna_gene_count": int(rna_gene_count),
    }


def _export_embeddings(data: QueryFormerData, model, adapters, out_dir: Path, device: str) -> None:
    arrays = _encode_all(data, model, adapters, device)
    out_dir.mkdir(parents=True, exist_ok=True)
    emb = pd.DataFrame(arrays["z_patient"], columns=[f"z_patient_{i:03d}" for i in range(arrays["z_patient"].shape[1])])
    frame = pd.concat([pd.DataFrame({"patient_id": data.patient_ids, "cancer_type": data.cancers, "split": data.split}), emb], axis=1)
    frame.to_parquet(out_dir / "test_patient_embeddings.parquet", index=False)
    with h5py.File(out_dir / "test_embeddings.h5", "w") as handle:
        handle.create_dataset("patient_ids", data=np.asarray(data.patient_ids, dtype="S"))
        handle.create_dataset("split", data=np.asarray(data.split, dtype="S"))
        for key, value in arrays.items():
            handle.create_dataset(key, data=value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--split-file")
    parser.add_argument("--output-dir", default="morpheus/outputs/v1_query_former")
    parser.add_argument("--checkpoint")
    parser.add_argument("--eval-only", action="store_true")
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    parser.add_argument("--variant", default="v1", choices=["v1", "wsi_hallmark_loss", "wsi_rna_hallmark_loss", "neighborhood_distill", "teacher_distill", "combined_best", "retentive_hallmark", "strong_teacher_hallmark", "teacher_hallmark_x1", "teacher_hallmark_x5", "wsi_hallmark_strong"])
    parser.add_argument("--teacher-npz")
    parser.add_argument("--teacher-key", default="wsi")
    args = parser.parse_args()
    if args.eval_only:
        if not args.checkpoint:
            raise SystemExit("--checkpoint is required with --eval-only")
        print(evaluate_query_former_checkpoint(args.checkpoint, args.config, args.split_file, args.output_dir, args.device))
    else:
        print(run_query_former_training(args.config, args.split_file, args.output_dir, args.epochs, args.batch_size, args.smoke, args.device, args.variant, args.teacher_npz, args.teacher_key))


if __name__ == "__main__":
    main()
