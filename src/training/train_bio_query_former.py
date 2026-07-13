"""Train MORPHEUS V2 BioQueryFormer with typed identity and biology states."""

from __future__ import annotations

import argparse
import copy
import json
from dataclasses import dataclass
from typing import Any
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from morpheus.src.data.wsi_patch_bags import build_patch_bag_registry, load_patch_bag
from morpheus.src.eval.retrieval_metrics import paired_retrieval_metrics
from morpheus.src.training.train_query_former import _clip_loss, _scale_trainfit, _split_labels, load_query_former_data
from morpheus.src.utils.config import load_config
from morpheus.src.utils.provenance import base_manifest, write_json


@dataclass
class BioQueryData:
    patient_ids: list[str]
    cancers: list[str]
    split: np.ndarray
    wsi_patient: np.ndarray
    rna: np.ndarray
    hallmark: np.ndarray
    hallmark_present: np.ndarray
    clinical: np.ndarray
    clinical_present: np.ndarray
    patch_paths: list[list[str]]
    hallmark_names: list[str]
    clinical_names: list[str]
    # A canonical 1536-D H-Optimus store.  Kept separate from patch_paths so
    # hoptimus mode can never accidentally read legacy 2048-D NPZ features.
    hoptimus_store: Any | None = None


def _load_wsi_patient(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        arr = handle["embeddings"][:].astype(np.float32)
        ids = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in handle["patient_ids"][:]]
    return pd.DataFrame({"patient_id": ids, "wsi_patient_vector": [row for row in arr]})


def _load_patch_registry(path: Path) -> pd.DataFrame:
    if not path.exists():
        build_patch_bag_registry(output_path=path, fast_filename_registry=True)
    return pd.read_parquet(path)


def _open_hoptimus_store(cfg) -> Any:
    """Open the canonical token store without importing it for legacy runs."""
    try:
        from morpheus.src.data.hoptimus_patch_store import HoptimusPatchStore
    except ImportError as exc:  # pragma: no cover - depends on extraction install
        raise RuntimeError("hoptimus_patch mode requires morpheus.src.data.hoptimus_patch_store") from exc
    root = Path(cfg.raw.get("hoptimus_patch_store_dir", cfg.path("wsi_standard_dir")))
    candidates = {
        "embeddings_path": root / "hoptimus_patch_embeddings.h5",
        "metadata_path": root / "hoptimus_patch_metadata.parquet",
        "index_path": root / "hoptimus_patient_patch_index.parquet",
    }
    missing = [str(path) for path in candidates.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("H-Optimus patch store is incomplete: " + ", ".join(missing))
    return HoptimusPatchStore(**candidates)


def load_bio_query_data(config_path: str, split_file: str | Path, wsi_mode: str = "patient") -> BioQueryData:
    base = load_query_former_data(config_path, split_file)
    cfg = load_config(config_path)
    wsi_patient = _load_wsi_patient(cfg.path("wsi_standard_dir") / "tcga_ut_hoptimus0_patient_embeddings.h5")
    patch_paths: list[list[str]] = [[] for _ in base.patient_ids]
    hoptimus_store = None
    if wsi_mode == "patch":
        patch_registry = _load_patch_registry(cfg.path("wsi_standard_dir") / "tcga_ut_patch_bag_registry.parquet")
        patch_by_patient = patch_registry.groupby("patient_id")["source_path"].apply(list).to_dict()
        patch_paths = [patch_by_patient.get(pid, []) for pid in base.patient_ids]
    elif wsi_mode == "hoptimus_patch":
        hoptimus_store = _open_hoptimus_store(cfg)
        patch_paths = [[pid] if hoptimus_store.load_patient_tokens(pid, max_tokens=1)[0].shape[0] else [] for pid in base.patient_ids]
    elif wsi_mode != "patient":
        raise ValueError(f"Unknown wsi_mode: {wsi_mode}")
    keep = np.ones(len(base.patient_ids), dtype=bool)
    if wsi_mode in {"patch", "hoptimus_patch"}:
        keep = np.asarray([len(paths) > 0 for paths in patch_paths], dtype=bool)
    patient_vectors = wsi_patient.set_index("patient_id").reindex(base.patient_ids)["wsi_patient_vector"].to_numpy()
    wsi_matrix = np.vstack(patient_vectors).astype(np.float32)
    wsi_matrix, _ = _scale_trainfit(wsi_matrix, base.split == "train")
    if not keep.all():
        base = BioQueryData(
            [base.patient_ids[i] for i in np.where(keep)[0]],
            [base.cancers[i] for i in np.where(keep)[0]],
            base.split[keep],
            wsi_matrix[keep],
            base.rna[keep],
            base.hallmark[keep],
            base.hallmark_present[keep],
            base.clinical[keep],
            base.clinical_present[keep],
            [patch_paths[i] for i in np.where(keep)[0]],
            base.hallmark_names,
            base.clinical_names,
            hoptimus_store,
        )
        return base
    return BioQueryData(
        base.patient_ids,
        base.cancers,
        base.split,
        wsi_matrix,
        base.rna,
        base.hallmark,
        base.hallmark_present,
        base.clinical,
        base.clinical_present,
        patch_paths,
        base.hallmark_names,
        base.clinical_names,
        hoptimus_store,
    )


def _batch_iter(indices: np.ndarray, batch_size: int, rng: np.random.Generator):
    order = indices.copy()
    rng.shuffle(order)
    for start in range(0, len(order), batch_size):
        yield order[start : start + batch_size]


def _balanced_batch_iter(indices: np.ndarray, cancers: list[str], batch_size: int, rng: np.random.Generator):
    by_cancer: dict[str, list[int]] = {}
    for idx in indices:
        by_cancer.setdefault(cancers[int(idx)], []).append(int(idx))
    pools = {cancer: rng.permutation(values).tolist() for cancer, values in by_cancer.items()}
    cancer_order = list(pools)
    emitted = 0
    while emitted < len(indices):
        rng.shuffle(cancer_order)
        batch = []
        for cancer in cancer_order:
            if not pools[cancer]:
                pools[cancer] = rng.permutation(by_cancer[cancer]).tolist()
            batch.append(pools[cancer].pop())
            if len(batch) == batch_size:
                break
        emitted += len(batch)
        yield np.asarray(batch, dtype=int)


def _siglip_loss(a, b, temperature: float = 0.07, ignore_negative_mask=None):
    import torch

    a = torch.nn.functional.normalize(a, dim=1)
    b = torch.nn.functional.normalize(b, dim=1)
    logits = a @ b.T / temperature
    labels = torch.eye(a.shape[0], dtype=torch.float32, device=a.device)
    weights = torch.ones_like(labels)
    if ignore_negative_mask is not None:
        weights = weights.masked_fill(ignore_negative_mask & ~labels.bool(), 0.0)
    loss_ab = torch.nn.functional.binary_cross_entropy_with_logits(logits, labels, weight=weights, reduction="sum") / weights.sum().clamp_min(1.0)
    loss_ba = torch.nn.functional.binary_cross_entropy_with_logits(logits.T, labels, weight=weights.T, reduction="sum") / weights.sum().clamp_min(1.0)
    return (loss_ab + loss_ba) / 2


def _false_negative_mask(cancers: list[str], hallmark: np.ndarray, hallmark_present: np.ndarray, device, threshold: float = 0.85):
    import torch

    n = len(cancers)
    mask = np.equal.outer(np.asarray(cancers), np.asarray(cancers))
    valid = hallmark_present.astype(bool)
    if valid.sum() > 1:
        h = hallmark.astype(np.float32)
        denom = np.linalg.norm(h, axis=1, keepdims=True)
        h = h / np.maximum(denom, 1e-12)
        sim = h @ h.T
        mask |= (sim >= threshold) & np.outer(valid, valid)
    np.fill_diagonal(mask, False)
    return torch.tensor(mask, dtype=torch.bool, device=device)


def _alignment_loss(a, b, mode: str, cancers: list[str], hallmark: np.ndarray, hallmark_present: np.ndarray, device, false_negative_aware: bool):
    mask = _false_negative_mask(cancers, hallmark, hallmark_present, device) if false_negative_aware else None
    if mode == "infonce":
        return _clip_loss(a, b)
    if mode == "siglip":
        return _siglip_loss(a, b, ignore_negative_mask=mask)
    if mode == "hybrid":
        return 0.5 * _clip_loss(a, b) + 0.5 * _siglip_loss(a, b, ignore_negative_mask=mask)
    raise ValueError(f"Unknown alignment loss: {mode}")


class _GradientReverse:
    @staticmethod
    def apply(x, scale: float):
        import torch

        class _Fn(torch.autograd.Function):
            @staticmethod
            def forward(ctx, value):
                ctx.scale = scale
                return value.view_as(value)

            @staticmethod
            def backward(ctx, grad_output):
                return -ctx.scale * grad_output

        return _Fn.apply(x)


def _encode_cancers(cancers: list[str], train_idx: np.ndarray) -> tuple[np.ndarray, list[str]]:
    train_classes = sorted({cancers[int(i)] for i in train_idx})
    mapping = {name: i for i, name in enumerate(train_classes)}
    labels = np.asarray([mapping.get(name, -1) for name in cancers], dtype=np.int64)
    return labels, train_classes


def _fit_pls_distillation_targets(data: BioQueryData, train_idx: np.ndarray, n_components: int) -> tuple[np.ndarray, np.ndarray] | None:
    if n_components <= 0 or len(train_idx) <= 2:
        return None
    from sklearn.cross_decomposition import PLSRegression

    n_comp = min(int(n_components), data.wsi_patient.shape[1], data.rna.shape[1], max(len(train_idx) - 1, 1))
    if n_comp <= 0:
        return None
    model = PLSRegression(n_components=n_comp, scale=True)
    model.fit(data.wsi_patient[train_idx], data.rna[train_idx])
    wsi_scores, rna_scores = model.transform(data.wsi_patient, data.rna)
    wsi_scores, _ = _scale_trainfit(wsi_scores.astype(np.float32), data.split == "train")
    rna_scores, _ = _scale_trainfit(rna_scores.astype(np.float32), data.split == "train")
    return wsi_scores.astype(np.float32), rna_scores.astype(np.float32)


def fit_train_only_residual_hallmarks(data: BioQueryData, train_idx: np.ndarray, k: int = 5) -> tuple[np.ndarray, dict[int, set[int]]]:
    """Residualize Hallmarks with train-only cancer statistics and build positives.

    Missing Hallmark rows remain zero and are excluded by ``hallmark_present`` at
    loss time.  Statistics for a cancer unseen in training use the train-global
    fallback, which is deliberately an OOD diagnostic rather than a headline
    within-cancer result.
    """
    values = data.hallmark.astype(np.float32)
    valid_train = np.asarray(train_idx)[data.hallmark_present[train_idx]]
    if len(valid_train) == 0:
        return np.zeros_like(values), {}
    global_mean = values[valid_train].mean(axis=0)
    global_scale = values[valid_train].std(axis=0).clip(1e-6)
    stats: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for cancer in sorted({data.cancers[int(i)] for i in valid_train}):
        rows = valid_train[np.asarray([data.cancers[int(i)] == cancer for i in valid_train])]
        # Tiny cancer cohorts are too unstable for a separate residual centre.
        if len(rows) >= 2:
            stats[cancer] = (values[rows].mean(axis=0), values[rows].std(axis=0).clip(1e-6))
    residual = np.zeros_like(values)
    for i, cancer in enumerate(data.cancers):
        mean, scale = stats.get(cancer, (global_mean, global_scale))
        residual[i] = (values[i] - mean) / scale
    neighbors: dict[int, set[int]] = {}
    for cancer in sorted({data.cancers[int(i)] for i in valid_train}):
        rows = valid_train[np.asarray([data.cancers[int(i)] == cancer for i in valid_train])]
        if len(rows) < 2:
            continue
        z = residual[rows]
        z = z / np.maximum(np.linalg.norm(z, axis=1, keepdims=True), 1e-12)
        sim = z @ z.T
        for pos, row in enumerate(rows):
            order = np.argsort(-sim[pos])
            nearest = [int(rows[j]) for j in order if j != pos][: min(k, len(rows) - 1)]
            neighbors[int(row)] = set(nearest)
    return residual, neighbors


def _residual_neighborhood_loss(z: "torch.Tensor", residual: "torch.Tensor", temperature: float = 0.2) -> "torch.Tensor":
    import torch

    if z.shape[0] < 2:
        return z.new_zeros(())
    target = torch.nn.functional.normalize(residual, dim=1)
    target_logits = target @ target.T / temperature
    pred = torch.nn.functional.normalize(z, dim=1)
    pred_logits = pred @ pred.T / temperature
    diagonal = torch.eye(z.shape[0], dtype=torch.bool, device=z.device)
    target_logits = target_logits.masked_fill(diagonal, torch.finfo(z.dtype).min)
    pred_logits = pred_logits.masked_fill(diagonal, torch.finfo(z.dtype).min)
    return -(torch.softmax(target_logits, dim=1) * torch.log_softmax(pred_logits, dim=1)).sum(dim=1).mean()


def _residual_supcon_loss(z: "torch.Tensor", indices: np.ndarray, neighbors: dict[int, set[int]], temperature: float = 0.2) -> "torch.Tensor":
    import torch

    if z.shape[0] < 2:
        return z.new_zeros(())
    sim = torch.nn.functional.normalize(z, dim=1) @ torch.nn.functional.normalize(z, dim=1).T / temperature
    diagonal = torch.eye(len(indices), dtype=torch.bool, device=z.device)
    sim = sim.masked_fill(diagonal, torch.finfo(z.dtype).min)
    losses = []
    for i, global_idx in enumerate(indices):
        positive = [j for j, other_idx in enumerate(indices) if int(other_idx) in neighbors.get(int(global_idx), set())]
        if positive:
            losses.append(-torch.log_softmax(sim[i], dim=0)[positive].mean())
    return torch.stack(losses).mean() if losses else z.new_zeros(())


def _centered_crosscov_loss(a: "torch.Tensor", b: "torch.Tensor") -> "torch.Tensor":
    import torch

    if a.shape[0] < 2:
        return a.new_zeros(())
    a = torch.nn.functional.normalize(a - a.mean(dim=0, keepdim=True), dim=0)
    b = torch.nn.functional.normalize(b - b.mean(dim=0, keepdim=True), dim=0)
    return ((a.T @ b) / a.shape[0]).square().mean()


def _load_patch_batch(paths_by_patient: list[list[str]], max_tokens: int, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    feats_list, coords_list, masks = [], [], []
    for i, paths in enumerate(paths_by_patient):
        if not paths:
            feats = np.zeros((1, 2048), dtype=np.float32)
            coords = np.zeros((1, 2), dtype=np.float32)
            mask = np.zeros(1, dtype=bool)
        else:
            bag = load_patch_bag(paths[(seed + i) % len(paths)], max_tokens=max_tokens, seed=seed + i)
            feats = bag.feats.astype(np.float32)
            coords = np.zeros((feats.shape[0], 2), dtype=np.float32) if bag.coords is None else bag.coords.astype(np.float32)
            mask = np.ones(feats.shape[0], dtype=bool)
        feats_list.append(feats)
        coords_list.append(coords)
        masks.append(mask)
    width = max(x.shape[0] for x in feats_list)
    feat_batch = np.zeros((len(feats_list), width, 2048), dtype=np.float32)
    coord_batch = np.zeros((len(feats_list), width, 2), dtype=np.float32)
    mask_batch = np.zeros((len(feats_list), width), dtype=bool)
    for i, (feats, coords, mask) in enumerate(zip(feats_list, coords_list, masks)):
        n = feats.shape[0]
        feat_batch[i, :n] = feats
        coord_batch[i, :n] = coords
        mask_batch[i, :n] = mask
    return feat_batch, coord_batch, mask_batch


def _load_hoptimus_patch_batch(store: Any, patient_ids: list[str], max_tokens: int, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load explicit 1536-D H-Optimus bags; never fall back to NPZ bags."""
    if store is None:
        raise RuntimeError("hoptimus_patch mode requires an open H-Optimus token store")
    bags = []
    for offset, patient_id in enumerate(patient_ids):
        feats, _ = store.load_patient_tokens(patient_id, max_tokens=max_tokens, seed=seed + offset, slide_balanced=True)
        feats = np.asarray(feats, dtype=np.float32)
        if feats.ndim != 2 or feats.shape[1] != 1536:
            raise ValueError(f"H-Optimus store returned {feats.shape} for {patient_id}; expected [tokens, 1536]")
        bags.append(feats)
    width = max(1, max((x.shape[0] for x in bags), default=0))
    feat_batch = np.zeros((len(bags), width, 1536), dtype=np.float32)
    mask_batch = np.zeros((len(bags), width), dtype=bool)
    for row, feats in enumerate(bags):
        if len(feats):
            feat_batch[row, : len(feats)] = feats
            mask_batch[row, : len(feats)] = True
    return feat_batch, np.zeros((len(bags), width, 2), dtype=np.float32), mask_batch


def _model_outputs(model, adapters, batch: dict, device: str, wsi_mode: str, use_wsi: bool = True, use_rna: bool = True, use_clinical: bool = False):
    import torch

    modalities = {}
    if wsi_mode == "patch":
        feats = torch.tensor(batch["patch_feats"], dtype=torch.float32, device=device)
        coords = torch.tensor(batch["patch_coords"], dtype=torch.float32, device=device)
        mask = torch.tensor(batch["patch_mask"], dtype=torch.bool, device=device)
        if not use_wsi:
            mask = torch.zeros_like(mask)
        adapter_key = "wsi_hoptimus_patch" if wsi_mode == "hoptimus_patch" else "wsi_patch"
        modalities["wsi"] = adapters[adapter_key](feats, mask, coords)
    else:
        wsi = torch.tensor(batch["wsi_patient"], dtype=torch.float32, device=device)
        present = torch.ones(wsi.shape[0], dtype=torch.bool, device=device) if use_wsi else torch.zeros(wsi.shape[0], dtype=torch.bool, device=device)
        modalities["wsi"] = adapters["wsi_patient"](wsi, present)
    rna = torch.tensor(batch["rna"], dtype=torch.float32, device=device)
    rna_present = torch.ones(rna.shape[0], dtype=torch.bool, device=device) if use_rna else torch.zeros(rna.shape[0], dtype=torch.bool, device=device)
    modalities["rna"] = adapters["rna"](rna, None, bulk_present=rna_present)
    if use_clinical:
        clinical = torch.tensor(batch["clinical"], dtype=torch.float32, device=device)
        clinical_present = torch.tensor(batch["clinical_present"], dtype=torch.bool, device=device)
        modalities["clinical"] = adapters["clinical"](clinical, clinical_present)
    return model(modalities)


def _encode_all(
    data: BioQueryData,
    model,
    adapters,
    device: str,
    wsi_mode: str,
    batch_size: int = 64,
    max_patch_tokens: int = 256,
    use_wsi: bool = True,
    use_rna: bool = True,
    use_clinical: bool = False,
    indices: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    import torch

    model.eval()
    adapters.eval()
    arrays = {key: [] for key in ("z_identity", "z_biology", "z_programs", "z_wsi_residual", "z_rna_residual", "z_uncertainty", "z_hypothesis")}
    selected = np.arange(len(data.patient_ids)) if indices is None else np.asarray(indices, dtype=int)
    with torch.no_grad():
        for start in range(0, len(selected), batch_size):
            idx = selected[start : start + batch_size]
            batch = {k: getattr(data, k)[idx] for k in ("wsi_patient", "rna", "hallmark", "hallmark_present", "clinical", "clinical_present")}
            if wsi_mode == "patch":
                feats, coords, mask = _load_patch_batch([data.patch_paths[i] for i in idx], max_patch_tokens, seed=start)
                batch.update({"patch_feats": feats, "patch_coords": coords, "patch_mask": mask})
            elif wsi_mode == "hoptimus_patch":
                feats, coords, mask = _load_hoptimus_patch_batch(data.hoptimus_store, [data.patient_ids[i] for i in idx], max_patch_tokens, seed=start)
                batch.update({"patch_feats": feats, "patch_coords": coords, "patch_mask": mask})
            out = _model_outputs(model, adapters, batch, device, wsi_mode, use_wsi, use_rna, use_clinical)
            for key in arrays:
                arrays[key].append(out[key].detach().cpu().numpy().astype(np.float32))
    return {key: np.vstack(value) for key, value in arrays.items()}


def _representation_diagnostics(identity: np.ndarray, biology: np.ndarray) -> dict[str, float]:
    if len(identity) < 2:
        return {}
    ident = identity / np.maximum(np.linalg.norm(identity, axis=1, keepdims=True), 1e-12)
    bio = biology / np.maximum(np.linalg.norm(biology, axis=1, keepdims=True), 1e-12)
    centered_i, centered_b = ident - ident.mean(0, keepdims=True), bio - bio.mean(0, keepdims=True)
    hsic = np.linalg.norm(centered_i.T @ centered_b, ord="fro") ** 2
    cka = float(hsic / max(np.linalg.norm(centered_i.T @ centered_i, ord="fro") * np.linalg.norm(centered_b.T @ centered_b, ord="fro"), 1e-12))
    sim_i, sim_b = ident @ ident.T, bio @ bio.T
    np.fill_diagonal(sim_i, -np.inf)
    np.fill_diagonal(sim_b, -np.inf)
    k = min(10, len(identity) - 1)
    overlap = np.mean([len(set(np.argpartition(row_i, -k)[-k:]) & set(np.argpartition(row_b, -k)[-k:])) / k for row_i, row_b in zip(sim_i, sim_b)])
    return {
        "identity_biology_abs_cosine": float(np.mean(np.abs(np.sum(ident * bio, axis=1)))),
        "identity_biology_linear_cka": cka,
        "identity_biology_nn10_overlap": float(overlap),
    }


def _evaluate(data: BioQueryData, model, adapters, heads, idx: np.ndarray, device: str, wsi_mode: str, max_patch_tokens: int, batch_size: int = 64, residual_hallmark: np.ndarray | None = None, probe_train_idx: np.ndarray | None = None) -> dict:
    if len(idx) < 2:
        return {"n": int(len(idx))}
    wsi_arrays = _encode_all(data, model, adapters, device, wsi_mode, batch_size=batch_size, max_patch_tokens=max_patch_tokens, use_rna=False, use_clinical=False, indices=idx)
    rna_arrays = _encode_all(data, model, adapters, device, wsi_mode, batch_size=batch_size, max_patch_tokens=max_patch_tokens, use_wsi=False, use_clinical=False, indices=idx)
    retrieval = paired_retrieval_metrics(wsi_arrays["z_identity"], rna_arrays["z_identity"], (1, 5, 10), [data.cancers[i] for i in idx], [data.cancers[i] for i in idx])
    present = data.hallmark_present[idx]
    pearson = None
    if present.any():
        import torch

        with torch.no_grad():
            x = torch.tensor(wsi_arrays["z_biology"][present], dtype=torch.float32, device=device)
            pred = heads["hallmark"](x).detach().cpu().numpy()
        truth = (residual_hallmark if residual_hallmark is not None else data.hallmark)[idx][present]
        vals = [float(np.corrcoef(truth[:, j], pred[:, j])[0, 1]) for j in range(truth.shape[1]) if np.std(truth[:, j]) > 0 and np.std(pred[:, j]) > 0]
        pearson = float(np.mean(vals)) if vals else None
    result = {
        "n": int(len(idx)),
        "retrieval_r1": retrieval["recall_at_1"],
        "retrieval_r5": retrieval["recall_at_5"],
        "retrieval_r10": retrieval["recall_at_10"],
        "retrieval_mrr": retrieval["mrr"],
        "same_cancer_at_10": retrieval.get("same_cancer_in_top10", 0.0),
        "hallmark_wsi_biology_pearson": pearson,
        "hallmark_n": int(present.sum()),
    }
    result.update(_representation_diagnostics(wsi_arrays["z_identity"], wsi_arrays["z_biology"]))
    if probe_train_idx is not None and len(probe_train_idx) >= 4 and len(idx) >= 1:
        try:
            from sklearn.linear_model import LogisticRegression

            train_arrays = _encode_all(data, model, adapters, device, wsi_mode, batch_size=batch_size, max_patch_tokens=max_patch_tokens, use_rna=False, use_clinical=False, indices=probe_train_idx)
            train_labels = np.asarray([data.cancers[i] for i in probe_train_idx])
            test_labels = np.asarray([data.cancers[i] for i in idx])
            if len(set(train_labels)) > 1:
                for name, key in (("identity", "z_identity"), ("biology", "z_biology")):
                    probe = LogisticRegression(max_iter=500, class_weight="balanced")
                    probe.fit(train_arrays[key], train_labels)
                    result[f"cancer_probe_{name}_accuracy"] = float((probe.predict(wsi_arrays[key]) == test_labels).mean())
        except ImportError:
            pass
    return result


def run_bio_query_former_training(
    config_path: str = "morpheus/configs/v1.json",
    split_file: str = "morpheus/data/processed/splits/tumor_state_stratified.json",
    output_dir: str | Path = "morpheus/outputs/v2_bio_query_former",
    wsi_mode: str = "patient",
    epochs: int = 20,
    batch_size: int = 32,
    max_patch_tokens: int = 256,
    smoke: bool = False,
    device_name: str = "auto",
    hidden_dim: int = 512,
    num_layers: int = 2,
    num_heads: int = 8,
    eval_batch_size: int = 64,
    torch_threads: int | None = None,
    alignment_loss: str = "infonce",
    false_negative_aware: bool = False,
    balanced_batches: bool = False,
    cancer_adv_weight: float = 0.0,
    pls_distill_weight: float = 0.0,
    pls_components: int = 32,
    biology_identity_decorrelation: float = 0.01,
    program_head_weight: float = 0.0,
    identity_warmup_epochs: int = 3,
    identity_reduced_weight: float = 0.2,
    rna_recon_weight: float = 0.03,
    wsi_slots: int = 32,
) -> Path:
    import torch
    from torch import nn

    from morpheus.src.encoders.clinical_adapter import ClinicalTokenAdapter
    from morpheus.src.encoders.rna_adapter import RNATokenAdapter
    from morpheus.src.encoders.wsi_adapter import WSIPatchTokenAdapter, WSITokenAdapter
    from morpheus.src.models.bio_query_former import BioQueryFormer, BioQueryFormerConfig

    cfg = load_config(config_path)
    seed = int(cfg.raw.get("seed", 42))
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    if torch_threads is not None and torch_threads > 0:
        torch.set_num_threads(int(torch_threads))
    data = load_bio_query_data(config_path, split_file, wsi_mode)
    if smoke:
        keep = np.r_[np.where(data.split == "train")[0][:96], np.where(data.split == "val")[0][:32], np.where(data.split == "test")[0][:32]]
        data = BioQueryData(
            [data.patient_ids[i] for i in keep],
            [data.cancers[i] for i in keep],
            data.split[keep],
            data.wsi_patient[keep],
            data.rna[keep],
            data.hallmark[keep],
            data.hallmark_present[keep],
            data.clinical[keep],
            data.clinical_present[keep],
            [data.patch_paths[i] for i in keep],
            data.hallmark_names,
            data.clinical_names,
            data.hoptimus_store,
        )
    device = "cuda" if device_name == "auto" and torch.cuda.is_available() else ("cpu" if device_name == "auto" else device_name)
    bcfg = BioQueryFormerConfig(
        hidden_dim=int(hidden_dim),
        num_layers=int(num_layers),
        num_heads=int(num_heads),
        program_slots=min(50, max(8, data.hallmark.shape[1])),
    )
    model = BioQueryFormer(bcfg).to(device)
    adapters = nn.ModuleDict(
        {
            "wsi_patient": WSITokenAdapter(data.wsi_patient.shape[1], bcfg.hidden_dim),
            "wsi_patch": WSIPatchTokenAdapter(2048, bcfg.hidden_dim),
            "wsi_hoptimus_patch": WSIPatchTokenAdapter(1536, bcfg.hidden_dim, use_coords=False, num_slots=int(wsi_slots)),
            "rna": RNATokenAdapter(data.rna.shape[1], data.hallmark.shape[1], bcfg.hidden_dim),
            "clinical": ClinicalTokenAdapter(data.clinical.shape[1], bcfg.hidden_dim),
        }
    ).to(device)
    heads = nn.ModuleDict(
        {
            "hallmark": nn.Linear(bcfg.hidden_dim, data.hallmark.shape[1]),
            "program_hallmark": nn.Linear(bcfg.hidden_dim, 1),
            "rna_recon": nn.Linear(bcfg.hidden_dim, data.rna.shape[1]),
        }
    ).to(device)
    train_idx = np.where(data.split == "train")[0]
    val_idx = np.where(data.split == "val")[0]
    test_idx = np.where(data.split == "test")[0]
    cancer_labels, cancer_classes = _encode_cancers(data.cancers, train_idx=train_idx)
    if cancer_adv_weight > 0 and len(cancer_classes) > 1:
        heads["cancer_adv"] = nn.Linear(bcfg.hidden_dim, len(cancer_classes)).to(device)
    pls_targets = _fit_pls_distillation_targets(data, train_idx, pls_components) if pls_distill_weight > 0 else None
    if pls_targets is not None:
        heads["pls_distill"] = nn.Linear(bcfg.hidden_dim, pls_targets[0].shape[1]).to(device)
    opt = torch.optim.AdamW(list(model.parameters()) + list(adapters.parameters()) + list(heads.parameters()), lr=1e-4, weight_decay=1e-2)
    history = []
    best_score = -np.inf
    best_state = None
    residual_hallmark, residual_neighbors = fit_train_only_residual_hallmarks(data, train_idx)
    identity_weight = 1.0
    previous_val_r10: float | None = None
    stable_retrieval_epochs = 0
    for epoch in range(epochs):
        model.train()
        adapters.train()
        heads.train()
        losses = []
        batches = _balanced_batch_iter(train_idx, data.cancers, batch_size, rng) if balanced_batches else _batch_iter(train_idx, batch_size, rng)
        for idx in batches:
            if len(idx) < 2:
                continue
            batch = {k: getattr(data, k)[idx] for k in ("wsi_patient", "rna", "hallmark", "hallmark_present", "clinical", "clinical_present")}
            if wsi_mode == "patch":
                feats, coords, mask = _load_patch_batch([data.patch_paths[i] for i in idx], max_patch_tokens, seed + epoch)
                batch.update({"patch_feats": feats, "patch_coords": coords, "patch_mask": mask})
            elif wsi_mode == "hoptimus_patch":
                feats, coords, mask = _load_hoptimus_patch_batch(
                    data.hoptimus_store, [data.patient_ids[i] for i in idx], max_patch_tokens, seed + epoch
                )
                batch.update({"patch_feats": feats, "patch_coords": coords, "patch_mask": mask})
            opt.zero_grad()
            out_wsi = _model_outputs(model, adapters, batch, device, wsi_mode, use_rna=False, use_clinical=False)
            out_rna = _model_outputs(model, adapters, batch, device, wsi_mode, use_wsi=False, use_clinical=False)
            out_full = _model_outputs(model, adapters, batch, device, wsi_mode, use_wsi=True, use_rna=True, use_clinical=True)
            hallmark_target = torch.tensor(residual_hallmark[idx], dtype=torch.float32, device=device)
            hallmark_present = torch.tensor(batch["hallmark_present"], dtype=torch.bool, device=device)
            rna_target = torch.tensor(batch["rna"], dtype=torch.float32, device=device)
            loss = identity_weight * _clip_loss(out_wsi["z_identity"], out_rna["z_identity"])
            if hallmark_present.any():
                loss = loss + torch.nn.functional.mse_loss(heads["hallmark"](out_wsi["z_biology"])[hallmark_present], hallmark_target[hallmark_present])
                loss = loss + torch.nn.functional.mse_loss(heads["hallmark"](out_rna["z_biology"])[hallmark_present], hallmark_target[hallmark_present])
                loss = loss + 0.25 * torch.nn.functional.mse_loss(heads["hallmark"](out_full["z_biology"])[hallmark_present], hallmark_target[hallmark_present])
                loss = loss + 0.5 * _residual_neighborhood_loss(out_wsi["z_biology"][hallmark_present], hallmark_target[hallmark_present])
                loss = loss + 0.25 * _residual_supcon_loss(out_wsi["z_biology"][hallmark_present], idx[hallmark_present.detach().cpu().numpy()], residual_neighbors)
                if program_head_weight > 0 and out_wsi["z_programs"].shape[1] == hallmark_target.shape[1]:
                    program_pred = heads["program_hallmark"](out_wsi["z_programs"]).squeeze(-1)
                    loss = loss + program_head_weight * torch.nn.functional.mse_loss(program_pred[hallmark_present], hallmark_target[hallmark_present])
            loss = loss + rna_recon_weight * torch.nn.functional.mse_loss(heads["rna_recon"](out_wsi["z_identity"]), rna_target)
            if pls_targets is not None:
                pls_wsi = torch.tensor(pls_targets[0][idx], dtype=torch.float32, device=device)
                pls_rna = torch.tensor(pls_targets[1][idx], dtype=torch.float32, device=device)
                loss = loss + pls_distill_weight * (
                    torch.nn.functional.mse_loss(heads["pls_distill"](out_wsi["z_biology"]), pls_wsi)
                    + torch.nn.functional.mse_loss(heads["pls_distill"](out_rna["z_biology"]), pls_rna)
                )
            if cancer_adv_weight > 0 and "cancer_adv" in heads:
                labels_np = cancer_labels[idx]
                valid = labels_np >= 0
                if valid.any():
                    cancer_target = torch.tensor(labels_np[valid], dtype=torch.long, device=device)
                    reversed_biology = _GradientReverse.apply(out_wsi["z_biology"][valid], cancer_adv_weight)
                    loss = loss + torch.nn.functional.cross_entropy(heads["cancer_adv"](reversed_biology), cancer_target)
            if biology_identity_decorrelation > 0:
                loss = loss + biology_identity_decorrelation * _centered_crosscov_loss(out_wsi["z_identity"], out_wsi["z_biology"])
            loss.backward()
            opt.step()
            losses.append(float(loss.item()))
        metrics = _evaluate(data, model, adapters, heads, val_idx, device, wsi_mode, max_patch_tokens, eval_batch_size, residual_hallmark, train_idx)
        val_r10 = float(metrics.get("retrieval_r10", 0.0))
        stable_retrieval_epochs = stable_retrieval_epochs + 1 if previous_val_r10 is not None and abs(val_r10 - previous_val_r10) < 0.005 else 0
        previous_val_r10 = val_r10
        if epoch + 1 >= identity_warmup_epochs and (stable_retrieval_epochs >= 2 or epoch + 1 >= 6):
            identity_weight = identity_reduced_weight
        score = float(metrics.get("retrieval_r10", 0.0)) + float(metrics.get("hallmark_wsi_biology_pearson") or 0.0)
        history.append({"epoch": epoch, "train_loss": float(np.mean(losses)) if losses else None, "identity_weight": identity_weight, **{f"val_{k}": v for k, v in metrics.items()}})
        print(
            f"epoch={epoch + 1}/{epochs} loss={history[-1]['train_loss']} "
            f"val_r10={metrics.get('retrieval_r10')} val_mrr={metrics.get('retrieval_mrr')} "
            f"val_hallmark={metrics.get('hallmark_wsi_biology_pearson')}",
            flush=True,
        )
        if score > best_score:
            best_score = score
            best_state = {
                "model": copy.deepcopy(model.state_dict()),
                "adapters": copy.deepcopy(adapters.state_dict()),
                "heads": copy.deepcopy(heads.state_dict()),
                "config": bcfg.__dict__,
                "wsi_mode": wsi_mode,
                "epoch": int(epoch + 1),
                "score": float(score),
                "objective": {
                    "alignment_loss": alignment_loss,
                    "false_negative_aware": bool(false_negative_aware),
                    "balanced_batches": bool(balanced_batches),
                    "cancer_adv_weight": float(cancer_adv_weight),
                    "pls_distill_weight": float(pls_distill_weight),
                    "pls_components": int(pls_components),
                    "biology_identity_decorrelation": float(biology_identity_decorrelation),
                    "program_head_weight": float(program_head_weight),
                },
            }
    out_dir = Path(output_dir)
    (out_dir / "checkpoints").mkdir(parents=True, exist_ok=True)
    if best_state is not None:
        model.load_state_dict(best_state["model"])
        adapters.load_state_dict(best_state["adapters"])
        heads.load_state_dict(best_state["heads"])
        torch.save(best_state, out_dir / "checkpoints" / "best.pt")
    pd.DataFrame(history).to_csv(out_dir / "train_log.csv", index=False)
    val_metrics = _evaluate(data, model, adapters, heads, val_idx, device, wsi_mode, max_patch_tokens, eval_batch_size, residual_hallmark, train_idx)
    test_metrics = _evaluate(data, model, adapters, heads, test_idx, device, wsi_mode, max_patch_tokens, eval_batch_size, residual_hallmark, train_idx)
    _export_embeddings(data, model, adapters, out_dir, device, wsi_mode, max_patch_tokens, eval_batch_size, heads)
    payload = base_manifest(cfg.project_root, cfg.config_path, seed)
    payload.update(
        {
            "model": "BioQueryFormer",
            "wsi_mode": wsi_mode,
            "device": device,
            "max_patch_tokens": int(max_patch_tokens),
            "hidden_dim": int(hidden_dim),
            "num_layers": int(num_layers),
            "num_heads": int(num_heads),
            "eval_batch_size": int(eval_batch_size),
            "torch_threads": int(torch_threads) if torch_threads is not None else None,
            "alignment_loss": alignment_loss,
            "false_negative_aware": bool(false_negative_aware),
            "balanced_batches": bool(balanced_batches),
            "cancer_adv_weight": float(cancer_adv_weight),
            "pls_distill_weight": float(pls_distill_weight),
            "pls_components": int(pls_components),
            "biology_identity_decorrelation": float(biology_identity_decorrelation),
            "program_head_weight": float(program_head_weight),
            "n_train": int(len(train_idx)),
            "n_val": int(len(val_idx)),
            "n_test": int(len(test_idx)),
            "best_epoch": int(best_state["epoch"]) if best_state is not None else None,
            "best_score": float(best_state["score"]) if best_state is not None else None,
            "val_metrics": val_metrics,
            "test_metrics": test_metrics,
            "checkpoint": str(out_dir / "checkpoints" / "best.pt"),
        }
    )
    out_path = out_dir / "test_metrics.json"
    write_json(out_path, payload)
    return out_path


def _program_scores(arrays: dict[str, np.ndarray], heads, device: str) -> np.ndarray:
    import torch

    if heads is None or "program_hallmark" not in heads:
        return np.empty((arrays["z_identity"].shape[0], 0), dtype=np.float32)
    values = []
    heads.eval()
    with torch.no_grad():
        for start in range(0, arrays["z_programs"].shape[0], 128):
            x = torch.tensor(arrays["z_programs"][start : start + 128], dtype=torch.float32, device=device)
            values.append(heads["program_hallmark"](x).squeeze(-1).detach().cpu().numpy().astype(np.float32))
    return np.vstack(values)


def _export_embeddings(data: BioQueryData, model, adapters, out_dir: Path, device: str, wsi_mode: str, max_patch_tokens: int, batch_size: int = 64, heads=None) -> Path:
    wsi = _encode_all(data, model, adapters, device, wsi_mode, batch_size=batch_size, max_patch_tokens=max_patch_tokens, use_rna=False, use_clinical=False)
    rna = _encode_all(data, model, adapters, device, wsi_mode, batch_size=batch_size, max_patch_tokens=max_patch_tokens, use_wsi=False, use_clinical=False)
    full = _encode_all(data, model, adapters, device, wsi_mode, batch_size=batch_size, max_patch_tokens=max_patch_tokens, use_clinical=True)
    wsi_program_scores = _program_scores(wsi, heads, device)
    full_program_scores = _program_scores(full, heads, device)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "bio_query_former_embeddings.npz"
    np.savez_compressed(
        out_path,
        patient_ids=np.asarray(data.patient_ids),
        split=np.asarray(data.split),
        cancers=np.asarray(data.cancers),
        program_names=np.asarray(data.hallmark_names[: wsi_program_scores.shape[1]]),
        wsi_identity=wsi["z_identity"],
        rna_identity=rna["z_identity"],
        wsi_biology=wsi["z_biology"],
        rna_biology=rna["z_biology"],
        wsi_programs=wsi["z_programs"],
        rna_programs=rna["z_programs"],
        full_programs=full["z_programs"],
        wsi_program_scores=wsi_program_scores,
        full_program_scores=full_program_scores,
        full_biology=full["z_biology"],
        wsi_residual=wsi["z_wsi_residual"],
        rna_residual=rna["z_rna_residual"],
        uncertainty=full["z_uncertainty"],
        hypothesis=full["z_hypothesis"],
    )
    return out_path


def evaluate_bio_query_former_checkpoint(
    checkpoint_path: str | Path,
    config_path: str = "morpheus/configs/v1.json",
    split_file: str = "morpheus/data/processed/splits/tumor_state_stratified.json",
    output_dir: str | Path = "morpheus/outputs/v2_bio_query_former_eval",
    wsi_mode: str | None = None,
    max_patch_tokens: int = 256,
    device_name: str = "auto",
    eval_batch_size: int = 64,
) -> Path:
    import torch
    from torch import nn

    from morpheus.src.encoders.clinical_adapter import ClinicalTokenAdapter
    from morpheus.src.encoders.rna_adapter import RNATokenAdapter
    from morpheus.src.encoders.wsi_adapter import WSIPatchTokenAdapter, WSITokenAdapter
    from morpheus.src.models.bio_query_former import BioQueryFormer, BioQueryFormerConfig

    cfg = load_config(config_path)
    seed = int(cfg.raw.get("seed", 42))
    device = "cuda" if device_name == "auto" and torch.cuda.is_available() else ("cpu" if device_name == "auto" else device_name)
    state = torch.load(checkpoint_path, map_location=device, weights_only=False)
    mode = wsi_mode or state.get("wsi_mode", "patient")
    data = load_bio_query_data(config_path, split_file, mode)
    bcfg = BioQueryFormerConfig(**state["config"])
    model = BioQueryFormer(bcfg).to(device)
    adapters = nn.ModuleDict(
        {
            "wsi_patient": WSITokenAdapter(data.wsi_patient.shape[1], bcfg.hidden_dim),
            "wsi_patch": WSIPatchTokenAdapter(2048, bcfg.hidden_dim),
            "rna": RNATokenAdapter(data.rna.shape[1], data.hallmark.shape[1], bcfg.hidden_dim),
            "clinical": ClinicalTokenAdapter(data.clinical.shape[1], bcfg.hidden_dim),
        }
    ).to(device)
    heads = nn.ModuleDict(
        {
            "hallmark": nn.Linear(bcfg.hidden_dim, data.hallmark.shape[1]),
            "program_hallmark": nn.Linear(bcfg.hidden_dim, 1),
            "rna_recon": nn.Linear(bcfg.hidden_dim, data.rna.shape[1]),
        }
    ).to(device)
    head_state = state.get("heads", {})
    if "cancer_adv.weight" in head_state:
        heads["cancer_adv"] = nn.Linear(bcfg.hidden_dim, head_state["cancer_adv.weight"].shape[0]).to(device)
    if "pls_distill.weight" in head_state:
        heads["pls_distill"] = nn.Linear(bcfg.hidden_dim, head_state["pls_distill.weight"].shape[0]).to(device)
    model.load_state_dict(state["model"])
    adapters.load_state_dict(state["adapters"])
    heads.load_state_dict(state["heads"], strict=False)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    val_idx = np.where(data.split == "val")[0]
    test_idx = np.where(data.split == "test")[0]
    val_metrics = _evaluate(data, model, adapters, heads, val_idx, device, mode, max_patch_tokens, eval_batch_size)
    test_metrics = _evaluate(data, model, adapters, heads, test_idx, device, mode, max_patch_tokens, eval_batch_size)
    _export_embeddings(data, model, adapters, out_dir, device, mode, max_patch_tokens, eval_batch_size, heads)
    payload = base_manifest(cfg.project_root, cfg.config_path, seed)
    payload.update(
        {
            "model": "BioQueryFormer",
            "checkpoint": str(checkpoint_path),
            "wsi_mode": mode,
            "device": device,
            "max_patch_tokens": int(max_patch_tokens),
            "hidden_dim": int(bcfg.hidden_dim),
            "num_layers": int(bcfg.num_layers),
            "num_heads": int(bcfg.num_heads),
            "eval_batch_size": int(eval_batch_size),
            "n_train": int((data.split == "train").sum()),
            "n_val": int(len(val_idx)),
            "n_test": int(len(test_idx)),
            "best_epoch": int(state["epoch"]) if "epoch" in state else None,
            "best_score": float(state["score"]) if "score" in state else None,
            "val_metrics": val_metrics,
            "test_metrics": test_metrics,
        }
    )
    out_path = out_dir / "test_metrics.json"
    write_json(out_path, payload)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    parser.add_argument("--split-file", default="morpheus/data/processed/splits/tumor_state_stratified.json")
    parser.add_argument("--output-dir", default="morpheus/outputs/v2_bio_query_former")
    parser.add_argument("--wsi-mode", choices=["patient", "patch", "hoptimus_patch"], default="patient")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--max-patch-tokens", type=int, default=256)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    parser.add_argument("--hidden-dim", type=int, default=512)
    parser.add_argument("--num-layers", type=int, default=2)
    parser.add_argument("--num-heads", type=int, default=8)
    parser.add_argument("--eval-batch-size", type=int, default=64)
    parser.add_argument("--torch-threads", type=int, default=None)
    parser.add_argument("--eval-only-checkpoint", default=None)
    parser.add_argument("--alignment-loss", choices=["infonce", "siglip", "hybrid"], default="infonce")
    parser.add_argument("--false-negative-aware", action="store_true")
    parser.add_argument("--balanced-batches", action="store_true")
    parser.add_argument("--cancer-adv-weight", type=float, default=0.0)
    parser.add_argument("--pls-distill-weight", type=float, default=0.0)
    parser.add_argument("--pls-components", type=int, default=32)
    parser.add_argument("--biology-identity-decorrelation", type=float, default=0.01)
    parser.add_argument("--program-head-weight", type=float, default=0.0)
    parser.add_argument("--identity-warmup-epochs", type=int, default=3)
    parser.add_argument("--identity-reduced-weight", type=float, default=0.2)
    parser.add_argument("--rna-recon-weight", type=float, default=0.03)
    parser.add_argument("--wsi-slots", type=int, default=32)
    args = parser.parse_args()
    if args.eval_only_checkpoint:
        print(
            evaluate_bio_query_former_checkpoint(
                args.eval_only_checkpoint,
                args.config,
                args.split_file,
                args.output_dir,
                args.wsi_mode,
                args.max_patch_tokens,
                args.device,
                args.eval_batch_size,
            )
        )
        return
    print(
        run_bio_query_former_training(
            config_path=args.config, split_file=args.split_file, output_dir=args.output_dir,
            wsi_mode=args.wsi_mode, epochs=args.epochs, batch_size=args.batch_size,
            max_patch_tokens=args.max_patch_tokens, smoke=args.smoke, device_name=args.device,
            hidden_dim=args.hidden_dim, num_layers=args.num_layers, num_heads=args.num_heads,
            eval_batch_size=args.eval_batch_size, torch_threads=args.torch_threads,
            alignment_loss=args.alignment_loss, false_negative_aware=args.false_negative_aware,
            balanced_batches=args.balanced_batches, cancer_adv_weight=args.cancer_adv_weight,
            pls_distill_weight=args.pls_distill_weight, pls_components=args.pls_components,
            biology_identity_decorrelation=args.biology_identity_decorrelation,
            program_head_weight=args.program_head_weight,
            identity_warmup_epochs=args.identity_warmup_epochs,
            identity_reduced_weight=args.identity_reduced_weight,
            rna_recon_weight=args.rna_recon_weight, wsi_slots=args.wsi_slots,
        )
    )


if __name__ == "__main__":
    main()
