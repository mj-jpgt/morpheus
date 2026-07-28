"""Refit the all-patch pooled H-Optimus MLP-CLIP teacher without leakage.

The previous MLP-CLIP artifacts were useful exploratory baselines, but their
fit population was not always explicit.  This entry point is the canonical
teacher builder for anchored V2.1: it reads *every* H-Optimus token belonging
to a patient, pools mean and standard deviation, fits preprocessing and the
two projection heads only on the requested development population, and writes
the versioned V2 representation contract.

The exported NPZ is deliberately sufficient for anchoring and frozen task
evaluation.  It does not export an untrained ``patient`` or ``biology`` view.

Example (final teacher refit):

``python -m morpheus.v2.refit_mlp_clip --data-config ... --split-file ... \
    --fit-population development_train_val --seed 42 --output mlp_seed42.npz``
"""
from __future__ import annotations

import argparse
import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from morpheus.src.training.train_bio_query_former import load_bio_query_data

from .baseline_exports import _write
from .contracts import require_state, validate_artifact
from .export_baselines import uncapped_patient_mean_std
from .preflight import sha256_json, validate_runtime_split
from .provenance import source_manifest


class MLPClipRefitError(RuntimeError):
    """Raised before an invalid teacher can be exported."""


@dataclass(frozen=True)
class MLPClipRefitConfig:
    """Fully recorded hyperparameters for a deterministic teacher refit."""

    embedding_dim: int = 256
    pca_dim: int = 512
    hidden_multiplier: int = 2
    epochs: int = 40
    batch_size: int = 256
    learning_rate: float = 3e-4
    weight_decay: float = 1e-2
    temperature: float = 0.07
    cancer_supcon_weight: float = 0.15
    loss_family: str = "clip"
    hard_negative_weight: float = 0.20
    hard_negative_margin: float = 0.10
    seed: int = 42
    use_bf16: bool = True


def _normalise(value):
    import torch

    return torch.nn.functional.normalize(value, dim=-1, eps=1e-8)


def _cancer_multi_positive_loss(logits, same_cancer):
    """Stable multi-positive contrastive loss used by the original MLP-CLIP.

    The identity diagonal remains a positive.  Rows therefore always have at
    least one target even if a minibatch contains one patient from a cancer.
    """
    import torch

    mask = same_cancer.to(dtype=torch.bool)
    log_denominator = torch.logsumexp(logits, dim=1)
    masked = logits.masked_fill(~mask, float("-inf"))
    log_numerator = torch.logsumexp(masked, dim=1) - torch.log(mask.sum(dim=1).to(logits.dtype))
    return -(log_numerator - log_denominator).mean()


def _siglip_loss(logits):
    """Pairwise sigmoid alignment with diagonal positives only."""
    import torch
    labels = torch.eye(len(logits), device=logits.device, dtype=logits.dtype).mul(2).sub(1)
    return -torch.nn.functional.logsigmoid(labels * logits).mean()


def _hard_negative_margin_loss(logits, margin: float):
    """Penalise the most confusable non-matching item in each direction."""
    import torch
    positive = logits.diagonal()
    negative = logits.masked_fill(torch.eye(len(logits), device=logits.device, dtype=torch.bool), -torch.inf).max(dim=1).values
    return torch.relu(float(margin) - positive + negative).mean()


def _fit_pca(features: np.ndarray, fit_rows: np.ndarray, maximum_components: int, seed: int):
    """Fit standardisation/PCA on the fit population only and transform all rows."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    if len(fit_rows) < 3:
        raise MLPClipRefitError("MLP-CLIP needs at least three fit-population patients")
    n_components = min(int(maximum_components), len(fit_rows) - 1, features.shape[1])
    if n_components < 2:
        raise MLPClipRefitError("MLP-CLIP PCA needs at least two components")
    scaler = StandardScaler(copy=True).fit(features[fit_rows])
    # Explicit randomized SVD keeps the Lambda preflight bounded for the
    # 3,072-D mean/std H-Optimus view while remaining seed-deterministic.
    pca = PCA(n_components=n_components, random_state=int(seed), svd_solver="randomized", iterated_power=5).fit(scaler.transform(features[fit_rows]))
    metadata = {
        "input_dim": int(features.shape[1]), "pca_dim": int(n_components),
        "scaler_mean_sha256": sha256_json(np.asarray(scaler.mean_, dtype=np.float64).round(10).tolist()),
        "pca_components_sha256": sha256_json(np.asarray(pca.components_, dtype=np.float64).round(10).tolist()),
    }
    # This state is deliberately checkpoint-only, not copied into the frozen
    # artifact.  It makes an optional model checkpoint genuinely reusable for
    # future *new-patient* inference rather than a misleading head-only file.
    checkpoint_state = {
        "scaler_mean": np.asarray(scaler.mean_, dtype=np.float32),
        "scaler_scale": np.asarray(scaler.scale_, dtype=np.float32),
        "pca_mean": np.asarray(pca.mean_, dtype=np.float32),
        "pca_components": np.asarray(pca.components_, dtype=np.float32),
    }
    return pca.transform(scaler.transform(features)).astype(np.float32), metadata, checkpoint_state


def _set_seed(seed: int) -> None:
    import torch

    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        # These flags materially speed dense projection training on Ampere;
        # they do not change the frozen input data or evaluation protocol.
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True


def _projection_head(in_dim: int, config: MLPClipRefitConfig):
    import torch
    from torch import nn

    width = max(config.embedding_dim, int(config.embedding_dim * config.hidden_multiplier))
    return nn.Sequential(
        nn.LayerNorm(in_dim),
        nn.Linear(in_dim, width),
        nn.GELU(),
        nn.Dropout(0.0),  # explicit: teacher training is deterministic, no hidden stochastic regulariser
        nn.Linear(width, config.embedding_dim),
    )


def _encode(head, values):
    return _normalise(head(values))


def train_mlp_clip(*, wsi: np.ndarray, rna: np.ndarray, cancers: Sequence[str], fit_rows: np.ndarray,
                   config: MLPClipRefitConfig, device: str | None = None) -> tuple[np.ndarray, np.ndarray, dict[str, object], dict[str, object]]:
    """Fit a symmetric MLP-CLIP teacher and encode every cohort patient.

    ``wsi`` must already be the uncapped all-patch mean/std representation.
    Only ``fit_rows`` affect preprocessing, gradients, or mini-batch order.
    """
    import torch

    wsi = np.asarray(wsi, dtype=np.float32); rna = np.asarray(rna, dtype=np.float32)
    fit_rows = np.asarray(fit_rows, dtype=np.int64)
    labels = np.asarray(cancers).astype(str)
    if wsi.ndim != 2 or rna.ndim != 2 or len(wsi) != len(rna) or len(labels) != len(wsi):
        raise MLPClipRefitError("WSI, RNA, and cancer labels must be aligned 2-D patient matrices")
    if len(np.unique(fit_rows)) != len(fit_rows) or (fit_rows < 0).any() or (fit_rows >= len(wsi)).any():
        raise MLPClipRefitError("fit rows are not a unique subset of the cohort")
    if not np.isfinite(wsi).all() or not np.isfinite(rna).all():
        raise MLPClipRefitError("teacher inputs must be finite")
    _set_seed(config.seed)
    wsi_pca, wsi_transform, wsi_checkpoint_transform = _fit_pca(wsi, fit_rows, config.pca_dim, config.seed)
    rna_pca, rna_transform, rna_checkpoint_transform = _fit_pca(rna, fit_rows, config.pca_dim, config.seed)
    resolved_device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    if resolved_device.startswith("cuda") and not torch.cuda.is_available():
        raise MLPClipRefitError("CUDA was requested but is unavailable")
    torch_device = torch.device(resolved_device)
    w_head = _projection_head(wsi_pca.shape[1], config).to(torch_device)
    r_head = _projection_head(rna_pca.shape[1], config).to(torch_device)
    optimizer = torch.optim.AdamW(list(w_head.parameters()) + list(r_head.parameters()),
                                  lr=config.learning_rate, weight_decay=config.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=max(1, config.epochs))
    features_wsi = torch.from_numpy(wsi_pca)
    features_rna = torch.from_numpy(rna_pca)
    fit_tensor = torch.as_tensor(fit_rows, dtype=torch.long)
    generator = torch.Generator(device="cpu").manual_seed(config.seed)
    history: list[float] = []
    use_amp = bool(config.use_bf16 and torch_device.type == "cuda" and torch.cuda.is_bf16_supported())
    for _epoch in range(config.epochs):
        w_head.train(); r_head.train()
        order = fit_tensor[torch.randperm(len(fit_tensor), generator=generator)]
        epoch_losses: list[float] = []
        for start in range(0, len(order), max(2, config.batch_size)):
            rows = order[start:start + max(2, config.batch_size)]
            if len(rows) < 2:
                continue
            x = features_wsi[rows].to(torch_device, non_blocking=True)
            y = features_rna[rows].to(torch_device, non_blocking=True)
            cancer = labels[rows.numpy()]
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type=torch_device.type, dtype=torch.bfloat16, enabled=use_amp):
                zw, zr = _encode(w_head, x), _encode(r_head, y)
                logits = (zw @ zr.T) / config.temperature
                target = torch.arange(len(rows), device=torch_device)
                if config.loss_family == "siglip":
                    # SigLIP uses independent pair labels rather than a batch
                    # softmax.  It intentionally has no cancer-label term.
                    loss = _siglip_loss(logits)
                else:
                    identity = (torch.nn.functional.cross_entropy(logits, target) +
                                torch.nn.functional.cross_entropy(logits.T, target)) / 2
                    same_cancer = torch.as_tensor(cancer[:, None] == cancer[None, :], device=torch_device)
                    supcon = (_cancer_multi_positive_loss(logits, same_cancer) +
                              _cancer_multi_positive_loss(logits.T, same_cancer.T)) / 2
                    loss = identity + config.cancer_supcon_weight * supcon
                    if config.loss_family == "hard_negative_clip":
                        hard_negative = (_hard_negative_margin_loss(logits, config.hard_negative_margin) +
                                         _hard_negative_margin_loss(logits.T, config.hard_negative_margin)) / 2
                        loss = loss + config.hard_negative_weight * hard_negative
            loss.backward()
            torch.nn.utils.clip_grad_norm_(list(w_head.parameters()) + list(r_head.parameters()), 5.0)
            optimizer.step()
            epoch_losses.append(float(loss.detach().float().cpu()))
        scheduler.step()
        if not epoch_losses:
            raise MLPClipRefitError("no MLP-CLIP batches contained two fit-population patients")
        history.append(float(np.mean(epoch_losses)))
    w_head.eval(); r_head.eval()
    with torch.no_grad():
        all_wsi = _encode(w_head, features_wsi.to(torch_device)).float().cpu().numpy().astype(np.float32)
        all_rna = _encode(r_head, features_rna.to(torch_device)).float().cpu().numpy().astype(np.float32)
    state = {
        "wsi_head": {key: value.detach().cpu() for key, value in w_head.state_dict().items()},
        "rna_head": {key: value.detach().cpu() for key, value in r_head.state_dict().items()},
        "wsi_transform": wsi_checkpoint_transform, "rna_transform": rna_checkpoint_transform,
        "fit_rows": fit_rows.tolist(), "config": asdict(config),
    }
    training = {"device": str(torch_device), "autocast_bf16": use_amp, "loss_history": history,
                "final_loss": history[-1], "epochs_completed": len(history),
                "n_fit_patients": int(len(fit_rows)), "n_fit_cancers": int(len(np.unique(labels[fit_rows])))}
    return all_wsi, all_rna, training, state


def refit_and_export(*, patient_ids: Sequence[str], cancers: Sequence[str], split: Sequence[str],
                     pooled_mean: np.ndarray, pooled_std: np.ndarray, rna: np.ndarray,
                     output: str | Path, fit_population: str, config: MLPClipRefitConfig,
                     source_provenance: dict[str, object] | None = None, device: str | None = None,
                     checkpoint: str | Path | None = None) -> dict[str, object]:
    """Train and write one contract-compliant MLP-CLIP teacher artifact."""
    ids = np.asarray(patient_ids).astype(str); cancer = np.asarray(cancers).astype(str); parts = np.asarray(split).astype(str)
    if pooled_mean.shape != pooled_std.shape or len(ids) != len(pooled_mean):
        raise MLPClipRefitError("uncapped mean/std must have equal aligned patient dimensions")
    if len(set(ids)) != len(ids):
        raise MLPClipRefitError("teacher cohort has duplicate patient IDs")
    masks = {"train_only": parts == "train", "development_train_val": np.isin(parts, ("train", "val"))}
    if fit_population not in masks:
        raise MLPClipRefitError("fit_population must be train_only or development_train_val")
    fit_rows = np.where(masks[fit_population])[0]
    # Mean and standard deviation are both exact sufficient statistics of the
    # complete patch bag.  No top-k selection or random cap is applied.
    wsi = np.concatenate([np.asarray(pooled_mean, dtype=np.float32), np.asarray(pooled_std, dtype=np.float32)], axis=1)
    z_wsi, z_rna, training, state = train_mlp_clip(wsi=wsi, rna=rna, cancers=cancer, fit_rows=fit_rows, config=config, device=device)
    config_manifest = {
        "family": "all_patch_pooled_mlp_clip", "pooling": "uncapped_patient_mean_plus_std",
        "patch_cap": None, "fit_population": fit_population, "fit_patient_ids_sha256": sha256_json(ids[fit_rows].tolist()),
        "fit_cancers": sorted(set(cancer[fit_rows].tolist())), "preprocessing": "fit_population_standard_scaler_plus_pca",
        "training": {**asdict(config), **training},
    }
    source = dict(source_provenance or {})
    source.update({"refit_source": source_manifest(configuration=config_manifest),
                   "patch_count_not_available": True})
    path = _write(output, patient_ids=ids, cancers=cancer, split=parts,
                  states={"wsi_identity": z_wsi, "rna_identity": z_rna}, method="mlp_clip_all_patch_meanstd_refit",
                  config=config_manifest, source_provenance=source)
    if checkpoint is not None:
        import torch

        checkpoint_path = Path(checkpoint); checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        torch.save({"format_version": 1, "state": state, "manifest": config_manifest,
                    "patient_ids": ids[fit_rows].tolist()}, checkpoint_path)
    return {"artifact": str(path), "fit_population": fit_population, "training": training,
            "trained_states": ["wsi_identity", "rna_identity"], "checkpoint": str(checkpoint) if checkpoint else None}


def build_teacher_from_canonical_data(*, data_config: str, split_file: str, output: str, fit_population: str,
                                      config: MLPClipRefitConfig, patient_read_batch: int = 8,
                                      expected_development_cancers: int = 11, expected_heldout_cancers: int = 21,
                                      device: str | None = None, checkpoint: str | None = None,
                                      pooled_cache: str | None = None) -> dict[str, object]:
    """Load one canonical paired cohort, then build a clean final teacher."""
    data = load_bio_query_data(data_config, split_file, wsi_mode="hoptimus_patch")
    split_manifest = validate_runtime_split(split_file, data.patient_ids, data.cancers, data.split,
                                            expected_development_cancers, expected_heldout_cancers)
    if pooled_cache:
        with np.load(pooled_cache, allow_pickle=False) as cached:
            required = {"patient_ids", "cancers", "split", "pooled_mean", "pooled_std", "patch_counts"}
            if not required.issubset(cached.files):
                raise MLPClipRefitError(f"pooled cache is missing {sorted(required - set(cached.files))}")
            cached_ids, cached_cancers, cached_split = (cached[name].astype(str) for name in ("patient_ids", "cancers", "split"))
            # ``load_bio_query_data`` deliberately exposes registry fields as
            # Python lists in some configurations.  Normalize at this cache
            # boundary so validation is representation-independent rather
            # than assuming ndarray methods on canonical metadata.
            patient_ids = np.asarray(data.patient_ids, dtype=str)
            cancers = np.asarray(data.cancers, dtype=str)
            split = np.asarray(data.split, dtype=str)
            if not (np.array_equal(cached_ids, patient_ids) and np.array_equal(cached_cancers, cancers) and np.array_equal(cached_split, split)):
                raise MLPClipRefitError("pooled cache patient/cancer/split order differs from canonical cohort")
            mean, std, counts = (cached[name].astype(np.float32 if name != "patch_counts" else np.int64) for name in ("pooled_mean", "pooled_std", "patch_counts"))
            store = {"source": "validated_uncapped_feature_cache", "path": str(Path(pooled_cache).resolve())}
    else:
        if data.hoptimus_store is None:
            raise MLPClipRefitError("canonical patch store is unavailable")
        store = data.hoptimus_store.validate()
        mean, std, counts = uncapped_patient_mean_std(data.hoptimus_store, data.patient_ids, read_batch=patient_read_batch)
    if (counts < 1).any():
        raise MLPClipRefitError("canonical paired cohort contains an empty H-Optimus bag")
    provenance = {"data_config": str(Path(data_config)), "split_file": str(Path(split_file)),
                  "split_digest": split_manifest["split_digest"], "cohort_digest": split_manifest["cohort_digest"],
                  "canonical_hoptimus_store": store, "patch_count_digest": sha256_json(counts.tolist()),
                  "token_count": {"min": int(counts.min()), "max": int(counts.max()), "total": int(counts.sum())}}
    result = refit_and_export(patient_ids=data.patient_ids, cancers=data.cancers, split=data.split,
                              pooled_mean=mean, pooled_std=std, rna=np.asarray(data.rna, dtype=np.float32), output=output,
                              fit_population=fit_population, config=config, source_provenance=provenance,
                              device=device, checkpoint=checkpoint)
    validate_artifact(output)
    with np.load(output, allow_pickle=False) as artifact:
        if require_state(artifact, "wsi_identity") is None or require_state(artifact, "rna_identity") is None:
            raise MLPClipRefitError("MLP-CLIP artifact did not export both directly trained identities")
    return {**result, "split": split_manifest, "source_provenance": provenance}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-config", required=True); parser.add_argument("--split-file", required=True); parser.add_argument("--output", required=True)
    parser.add_argument("--fit-population", choices=("train_only", "development_train_val"), required=True)
    parser.add_argument("--seed", type=int, default=42); parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--batch-size", type=int, default=256); parser.add_argument("--embedding-dim", type=int, default=256)
    parser.add_argument("--pca-dim", type=int, default=512); parser.add_argument("--hidden-multiplier", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=3e-4); parser.add_argument("--weight-decay", type=float, default=1e-2)
    parser.add_argument("--temperature", type=float, default=.07); parser.add_argument("--cancer-supcon-weight", type=float, default=.15)
    parser.add_argument("--loss-family", choices=("clip", "siglip", "hard_negative_clip"), default="clip")
    parser.add_argument("--hard-negative-weight", type=float, default=.20); parser.add_argument("--hard-negative-margin", type=float, default=.10)
    parser.add_argument("--patient-read-batch", type=int, default=8); parser.add_argument("--device", default=None); parser.add_argument("--checkpoint", default=None); parser.add_argument("--pooled-cache", default="")
    parser.add_argument("--no-bf16", action="store_true")
    parser.add_argument("--expected-development-cancers", type=int, default=11); parser.add_argument("--expected-heldout-cancers", type=int, default=21)
    args = parser.parse_args()
    if args.epochs < 1 or args.batch_size < 2 or args.embedding_dim < 2 or args.pca_dim < 2 or args.hard_negative_weight < 0:
        parser.error("epochs >= 1, batch-size >= 2, embedding-dim >= 2, and pca-dim >= 2 are required")
    config = MLPClipRefitConfig(embedding_dim=args.embedding_dim, pca_dim=args.pca_dim, hidden_multiplier=args.hidden_multiplier,
                                epochs=args.epochs, batch_size=args.batch_size, learning_rate=args.learning_rate,
                                weight_decay=args.weight_decay, temperature=args.temperature,
                                cancer_supcon_weight=args.cancer_supcon_weight, loss_family=args.loss_family,
                                hard_negative_weight=args.hard_negative_weight, hard_negative_margin=args.hard_negative_margin,
                                seed=args.seed, use_bf16=not args.no_bf16)
    result = build_teacher_from_canonical_data(data_config=args.data_config, split_file=args.split_file, output=args.output,
                                               fit_population=args.fit_population, config=config,
                                               patient_read_batch=args.patient_read_batch,
                                               expected_development_cancers=args.expected_development_cancers,
                                               expected_heldout_cancers=args.expected_heldout_cancers,
                                               device=args.device, checkpoint=args.checkpoint,
                                               pooled_cache=args.pooled_cache or None)
    print(json.dumps(result, indent=2, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
