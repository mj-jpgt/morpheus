"""Train survival probes on frozen MORPHEUS foundation embeddings."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd

from morpheus.src.encoders.bulkformer_encoder import load_embedding_store
from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import base_manifest, write_json


@dataclass
class SurvivalDataset:
    patient_ids: list[str]
    cancers: list[str]
    time: np.ndarray
    event: np.ndarray
    features: dict[str, np.ndarray]


def _load_wsi(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        embeddings = handle["embeddings"][:].astype(np.float32)
        patient_ids = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in handle["patient_ids"][:]]
    return pd.DataFrame({"patient_id": patient_ids, "wsi_vector": [row for row in embeddings]})


def _load_rna(embedding_path: Path, metadata_path: Path, primary_tumor_only: bool = True) -> pd.DataFrame:
    embeddings = load_embedding_store(embedding_path)
    metadata = pd.read_parquet(metadata_path).copy()
    metadata["sample_type"] = metadata["sample_id"].astype(str).str.split("-").str[3].str[:2]
    joined = metadata.merge(embeddings, on="patient_id", how="inner")
    if primary_tumor_only:
        joined = joined[joined["sample_type"].isin(["01", "03"])].copy()
    numeric = [c for c in embeddings.columns if c != "patient_id"]
    rows = []
    for patient_id, group in joined.groupby("patient_id", dropna=True):
        matrix = group[numeric].to_numpy(dtype=np.float32)
        cancer = str(group["cancer_type"].mode().iloc[0] if not group["cancer_type"].mode().empty else group["cancer_type"].iloc[0])
        rows.append({"patient_id": patient_id, "cancer_type": cancer, "rna_vector": matrix.mean(axis=0)})
    return pd.DataFrame(rows)


def load_survival_dataset(config_path: str, modalities: Iterable[str]) -> SurvivalDataset:
    cfg = load_config(config_path)
    modalities = set(modalities)
    survival_path = cfg.project_root / "meta-intersurv" / "data" / "embeddings" / "patient_master" / "survival_labels.parquet"
    if survival_path.exists():
        master = pd.read_parquet(survival_path).rename(columns={"time": "survival_time", "event": "survival_event"})
        master = master[["patient_id", "cancer_type", "survival_time", "survival_event"]].dropna().copy()
    else:
        master = pd.read_parquet(cfg.path("processed_dir") / "master_patient_table.parquet")
        master = master[["patient_id", "cancer_type", "survival_time", "survival_event"]].dropna().copy()
    master = master[(master["survival_time"] > 0) & master["survival_event"].isin([0, 1])].copy()
    rna_path = cfg.path("rna_bulkformer_embeddings")
    rna_meta_path = rna_path.with_name("tcga_bulkformer_embedding_metadata.parquet")
    rna = _load_rna(rna_path, rna_meta_path)
    frame = master.merge(rna[["patient_id", "rna_vector"]], on="patient_id", how="inner")
    features: dict[str, np.ndarray] = {}
    if "wsi" in modalities:
        wsi = _load_wsi(cfg.path("wsi_standard_dir") / "tcga_ut_hoptimus0_patient_embeddings.h5")
        frame = frame.merge(wsi, on="patient_id", how="inner")
    frame = frame.sort_values("patient_id").reset_index(drop=True)
    if "rna" in modalities:
        features["rna"] = np.vstack(frame["rna_vector"].to_numpy()).astype(np.float32)
    if "wsi" in modalities:
        features["wsi"] = np.vstack(frame["wsi_vector"].to_numpy()).astype(np.float32)
    if not features:
        raise ValueError("At least one modality is required")
    return SurvivalDataset(
        patient_ids=frame["patient_id"].astype(str).tolist(),
        cancers=frame["cancer_type"].astype(str).tolist(),
        time=frame["survival_time"].to_numpy(dtype=np.float32),
        event=frame["survival_event"].to_numpy(dtype=bool),
        features=features,
    )


def _top_train_split(cancers: list[str], n_train_cancers: int) -> tuple[list[str], list[str]]:
    counts = pd.Series(cancers).value_counts()
    train = counts.head(n_train_cancers).index.astype(str).tolist()
    test = counts.loc[~counts.index.isin(train)].index.astype(str).tolist()
    return train, test


def _available_zero_shot_split(cancers: list[str]) -> tuple[list[str], list[str]]:
    counts = pd.Series(cancers).value_counts()
    n_train = max(1, len(counts) // 2 + len(counts) % 2)
    train = counts.head(n_train).index.astype(str).tolist()
    test = counts.iloc[n_train:].index.astype(str).tolist()
    return train, test


def _make_masks(cancers: list[str], train_cancers: list[str], test_cancers: list[str], seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cancer_arr = np.asarray(cancers)
    train_domain = np.isin(cancer_arr, train_cancers)
    test = np.isin(cancer_arr, test_cancers)
    rng = np.random.default_rng(seed)
    train = np.zeros(len(cancers), dtype=bool)
    val = np.zeros(len(cancers), dtype=bool)
    for cancer in train_cancers:
        idx = np.where(cancer_arr == cancer)[0]
        rng.shuffle(idx)
        n_val = max(1, int(round(len(idx) * 0.2))) if len(idx) >= 5 else 0
        val[idx[:n_val]] = True
        train[idx[n_val:]] = True
    train &= train_domain
    return train, val, test


def _prepare_features(data: SurvivalDataset, modalities: list[str], train_mask: np.ndarray, max_components: int) -> tuple[np.ndarray, dict]:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    matrix = np.concatenate([data.features[m] for m in modalities], axis=1).astype(np.float32)
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(matrix[train_mask])
    scaled = scaler.transform(matrix)
    n_components = min(max_components, train_scaled.shape[0] - 1, scaled.shape[1])
    pca = None
    if n_components >= 2 and scaled.shape[1] > n_components:
        pca = PCA(n_components=n_components, random_state=42)
        pca.fit(train_scaled)
        transformed = pca.transform(scaled).astype(np.float32)
    else:
        transformed = scaled.astype(np.float32)
    return transformed, {"input_dim": int(matrix.shape[1]), "feature_dim": int(transformed.shape[1]), "pca_components": int(transformed.shape[1]) if pca is not None else None}


def _cox_loss(risk, time, event):
    import torch

    order = torch.argsort(time, descending=True)
    risk = risk[order]
    event = event[order].float()
    log_cumsum = torch.logcumsumexp(risk, dim=0)
    observed = event.sum().clamp_min(1.0)
    return -((risk - log_cumsum) * event).sum() / observed


def _fit_cox_mlp(x: np.ndarray, time: np.ndarray, event: np.ndarray, train_mask: np.ndarray, val_mask: np.ndarray, epochs: int, lr: float, seed: int) -> tuple[np.ndarray, dict]:
    import torch
    from torch import nn

    torch.manual_seed(seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    hidden = min(128, max(16, x.shape[1]))
    model = nn.Sequential(nn.Linear(x.shape[1], hidden), nn.ReLU(), nn.Dropout(0.15), nn.Linear(hidden, 1)).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-2)
    xt = torch.tensor(x, dtype=torch.float32, device=device)
    tt = torch.tensor(time, dtype=torch.float32, device=device)
    et = torch.tensor(event, dtype=torch.bool, device=device)
    train_idx = torch.tensor(np.where(train_mask)[0], dtype=torch.long, device=device)
    best_state = None
    best_val = -np.inf
    patience = 25
    stale = 0
    history = []
    for epoch in range(epochs):
        model.train()
        opt.zero_grad()
        risk = model(xt[train_idx]).squeeze(1)
        loss = _cox_loss(risk, tt[train_idx], et[train_idx])
        loss.backward()
        opt.step()
        with torch.no_grad():
            model.eval()
            pred = model(xt).squeeze(1).detach().cpu().numpy()
        val_c = _harrell_cindex(event[val_mask], time[val_mask], pred[val_mask]) if val_mask.any() else _harrell_cindex(event[train_mask], time[train_mask], pred[train_mask])
        history.append({"epoch": epoch, "train_loss": float(loss.item()), "val_cindex": float(val_c)})
        if val_c > best_val:
            best_val = val_c
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
        if stale >= patience:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    with torch.no_grad():
        model.eval()
        risk = model(xt).squeeze(1).detach().cpu().numpy().astype(np.float32)
    return risk, {"epochs_run": len(history), "best_val_cindex": float(best_val), "device": device, "history": history}


def _surv_array(event: np.ndarray, time: np.ndarray):
    return np.asarray([(bool(e), float(t)) for e, t in zip(event, time)], dtype=[("event", "?"), ("time", "<f8")])


def _harrell_cindex(event: np.ndarray, time: np.ndarray, risk: np.ndarray) -> float:
    from sksurv.metrics import concordance_index_censored

    if len(risk) < 2 or np.unique(event).size < 2:
        return float("nan")
    try:
        return float(concordance_index_censored(event.astype(bool), time.astype(float), risk.astype(float))[0])
    except Exception:
        return float("nan")


def _ipcw_cindex(train_event: np.ndarray, train_time: np.ndarray, event: np.ndarray, time: np.ndarray, risk: np.ndarray) -> float | None:
    from sksurv.metrics import concordance_index_ipcw

    if len(risk) < 2 or np.unique(event).size < 2:
        return None
    train_surv = _surv_array(train_event, train_time)
    test_surv = _surv_array(event, time)
    tau = min(float(np.max(train_time)) - 1e-6, float(np.max(time)) - 1e-6)
    try:
        return float(concordance_index_ipcw(train_surv, test_surv, risk.astype(float), tau=tau)[0])
    except Exception:
        return None


def _bootstrap_cindex(event: np.ndarray, time: np.ndarray, risk: np.ndarray, n_bootstrap: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    vals = []
    idx = np.arange(len(risk))
    for _ in range(n_bootstrap):
        sample = rng.choice(idx, size=len(idx), replace=True)
        if np.unique(event[sample]).size < 2:
            continue
        val = _harrell_cindex(event[sample], time[sample], risk[sample])
        if not np.isnan(val):
            vals.append(val)
    if not vals:
        return {"harrell_cindex_ci95": None, "n_bootstrap_used": 0}
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return {"harrell_cindex_ci95": [float(lo), float(hi)], "n_bootstrap_used": int(len(vals))}


def _metrics(data: SurvivalDataset, risk: np.ndarray, train_mask: np.ndarray, eval_mask: np.ndarray, n_bootstrap: int, seed: int) -> dict:
    event = data.event[eval_mask]
    time = data.time[eval_mask]
    out = {
        "n": int(eval_mask.sum()),
        "events": int(event.sum()),
        "censored": int((~event).sum()),
        "harrell_cindex": _harrell_cindex(event, time, risk[eval_mask]),
        "ipcw_cindex": _ipcw_cindex(data.event[train_mask], data.time[train_mask], event, time, risk[eval_mask]),
    }
    out.update(_bootstrap_cindex(event, time, risk[eval_mask], n_bootstrap, seed))
    return out


def _per_cancer_metrics(data: SurvivalDataset, risk: np.ndarray, mask: np.ndarray, min_n: int) -> pd.DataFrame:
    rows = []
    cancer_arr = np.asarray(data.cancers)
    for cancer in sorted(set(cancer_arr[mask])):
        cmask = mask & (cancer_arr == cancer)
        if cmask.sum() < min_n:
            continue
        rows.append(
            {
                "cancer_type": cancer,
                "n": int(cmask.sum()),
                "events": int(data.event[cmask].sum()),
                "harrell_cindex": _harrell_cindex(data.event[cmask], data.time[cmask], risk[cmask]),
            }
        )
    return pd.DataFrame(rows)


def run_survival_probe(
    config_path: str,
    modalities: list[str],
    split_mode: str,
    n_train_cancers: int,
    output_dir: str | Path,
    epochs: int,
    max_components: int,
    n_bootstrap: int,
    allow_available_fallback: bool,
) -> Path:
    cfg = load_config(config_path)
    seed = int(cfg.raw.get("seed", 42))
    data = load_survival_dataset(config_path, modalities)
    cancers = sorted(set(data.cancers))
    requested_split = {"mode": split_mode, "n_train_cancers": n_train_cancers}
    infeasible_reason = None
    if split_mode == "train11_test22":
        train_cancers, test_cancers = _top_train_split(data.cancers, n_train_cancers)
        expected_test_cancers = 22
        if len(train_cancers) < n_train_cancers or len(test_cancers) != expected_test_cancers:
            infeasible_reason = f"Requested {n_train_cancers}/{expected_test_cancers} cancer split, but eligible modalities have {len(train_cancers)}/{len(test_cancers)} across {len(cancers)} cancer types: {cancers}"
            if not allow_available_fallback:
                raise ValueError(infeasible_reason)
            split_mode = f"available_top{len(train_cancers)}_test{len(test_cancers)}"
    elif split_mode == "available_zero_shot":
        train_cancers, test_cancers = _available_zero_shot_split(data.cancers)
    else:
        raise ValueError(f"Unknown split mode: {split_mode}")
    train_mask, val_mask, test_mask = _make_masks(data.cancers, train_cancers, test_cancers, seed)
    if not train_mask.any() or not test_mask.any():
        raise ValueError("Split produced empty train or test cohort")
    x, feature_info = _prepare_features(data, modalities, train_mask, max_components)
    risk, training = _fit_cox_mlp(x, data.time, data.event, train_mask, val_mask, epochs, 3e-4, seed)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    predictions = pd.DataFrame(
        {
            "patient_id": data.patient_ids,
            "cancer_type": data.cancers,
            "split": np.where(train_mask, "train", np.where(val_mask, "val", np.where(test_mask, "test", "excluded"))),
            "survival_time": data.time,
            "survival_event": data.event.astype(int),
            "risk": risk,
        }
    )
    predictions.to_parquet(out_dir / "predictions.parquet", index=False)
    per_cancer = _per_cancer_metrics(data, risk, test_mask, min_n=10)
    per_cancer.to_csv(out_dir / "per_cancer_metrics.csv", index=False)
    payload = base_manifest(cfg.project_root, cfg.config_path, seed)
    payload.update(
        {
            "model": "frozen_foundation_multimodal_cox_mlp",
            "modalities": modalities,
            "requested_split": requested_split,
            "actual_split_mode": split_mode,
            "infeasible_requested_split_reason": infeasible_reason,
            "train_cancers": train_cancers,
            "test_cancers": test_cancers,
            "eligible_cancers": cancers,
            "feature_info": feature_info,
            "training": {k: v for k, v in training.items() if k != "history"},
            "n_total": len(data.patient_ids),
            "n_train": int(train_mask.sum()),
            "n_val": int(val_mask.sum()),
            "n_test": int(test_mask.sum()),
            "train_metrics": _metrics(data, risk, train_mask, train_mask, n_bootstrap, seed),
            "val_metrics": _metrics(data, risk, train_mask, val_mask, n_bootstrap, seed + 1),
            "test_metrics": _metrics(data, risk, train_mask, test_mask, n_bootstrap, seed + 2),
            "notes": [
                "Foundation encoders are frozen; only a small Cox MLP head is trained.",
                "Cancer labels are taken from BulkFormer metadata, not the current master-table tissue-source code field.",
            ],
        }
    )
    write_json(out_dir / "metrics.json", payload)
    pd.DataFrame(training["history"]).to_csv(out_dir / "train_log.csv", index=False)
    return out_dir / "metrics.json"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--modalities", default="wsi,rna")
    parser.add_argument("--split-mode", choices=["train11_test22", "available_zero_shot"], default="train11_test22")
    parser.add_argument("--n-train-cancers", type=int, default=11)
    parser.add_argument("--output-dir", default="morpheus/outputs/v1_survival/multimodal_train11_test22")
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--max-components", type=int, default=128)
    parser.add_argument("--n-bootstrap", type=int, default=200)
    parser.add_argument("--allow-available-fallback", action="store_true")
    args = parser.parse_args()
    modalities = [m.strip() for m in args.modalities.split(",") if m.strip()]
    print(
        run_survival_probe(
            args.config,
            modalities,
            args.split_mode,
            args.n_train_cancers,
            args.output_dir,
            args.epochs,
            args.max_components,
            args.n_bootstrap,
            args.allow_available_fallback,
        )
    )


if __name__ == "__main__":
    main()
