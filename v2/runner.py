"""Canonical uncapped TCGA V2 training path.

This module deliberately reads the canonical H-Optimus store directly.  It
does not call the legacy capped patch-batch helper.
"""

from __future__ import annotations

import argparse
import json
import warnings
from dataclasses import asdict
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import torch

from morpheus.src.training.train_bio_query_former import load_bio_query_data

from .data import DynamicTokenBatchSampler, TrainOnlyStandardizer
from .model import TumorStateV2, V2ModelConfig
from .training import V2LossSchedule, V2Trainer
from .contracts import trainable_state_names
from .plip import canonical_patient_source_digests, load_plip_teacher_cache
from .preflight import validate_runtime_split
from .provenance import source_manifest
from .slide_pretraining import SlidePretrainingConfig, SlidePretrainingObjective, SlidePretrainer
from morpheus.src.eval.retrieval_metrics import paired_retrieval_metrics


def residualise_programmes(values: np.ndarray, cancers: np.ndarray, train: np.ndarray) -> tuple[np.ndarray, dict[str, list[float]]]:
    """Cancer-residualise and scale programmes using development patients only."""
    global_mean, global_scale = values[train].mean(0), values[train].std(0).clip(1e-6)
    stats = {str(c): (values[train & (cancers == c)].mean(0), values[train & (cancers == c)].std(0).clip(1e-6))
             for c in np.unique(cancers[train])}
    fitted = [stats.get(str(c), (global_mean, global_scale)) for c in cancers]
    means = {c: mean.tolist() for c, (mean, _) in stats.items()}
    return np.vstack([(values[i] - mean) / scale for i, (mean, scale) in enumerate(fitted)]).astype(np.float32), means


def _metadata_coordinates(metadata: pd.DataFrame) -> np.ndarray | None:
    """Return validated tile coordinates or None; never manufacture zeros."""
    for x_name, y_name in (("x", "y"), ("x_coord", "y_coord"), ("coord_x", "coord_y")):
        if {x_name, y_name}.issubset(metadata.columns):
            values = metadata[[x_name, y_name]].to_numpy(np.float32)
            if np.isfinite(values).all() and np.ptp(values, axis=0).any():
                return values
    return None


def _attach_numeric_table(data, path: str, name: str) -> dict[str, object]:
    """Patient-ID align a tabular modality and fit all transforms on train only."""
    table = pd.read_parquet(path) if str(path).lower().endswith("parquet") else pd.read_csv(path, sep=None, engine="python")
    patient_column = next((column for column in ("patient_id", "Patient ID", "bcr_patient_barcode") if column in table.columns), None)
    if patient_column is None:
        raise ValueError(f"{name} table needs a patient_id-like column")
    table = table.copy(); table[patient_column] = table[patient_column].astype(str)
    if table[patient_column].duplicated().any():
        raise ValueError(f"{name} table has duplicate patient IDs")
    numeric = [column for column in table.columns if column != patient_column and pd.api.types.is_numeric_dtype(table[column])]
    if not numeric:
        raise ValueError(f"{name} table has no numeric feature columns")
    values = table.set_index(patient_column).reindex(data.patient_ids)[numeric].to_numpy(np.float32)
    present = np.isfinite(values).any(axis=1)
    train = np.asarray(data.split) == "train"
    fit = np.where(train & present)[0]
    if len(fit) < 2:
        raise ValueError(f"{name} has insufficient train coverage")
    medians = np.nanmedian(values[fit], axis=0)
    values = np.where(np.isfinite(values), values, medians[None, :]).astype(np.float32)
    standardizer = TrainOnlyStandardizer().fit(values, fit)
    setattr(data, f"_v2_{name}", standardizer.transform(values))
    setattr(data, f"_v2_{name}_present", present.astype(bool))
    return {"columns": numeric, "n_present": int(present.sum()), "fit_patients": [data.patient_ids[i] for i in fit]}


def _standardize_clinical(data, fit_mask: np.ndarray | None = None) -> dict[str, object]:
    train = np.asarray(data.split) == "train" if fit_mask is None else np.asarray(fit_mask, dtype=bool)
    present = np.asarray(data.clinical_present, dtype=bool)
    values = np.asarray(data.clinical, dtype=np.float32)
    fit = np.where(train & present)[0]
    if len(fit) < 2:
        return {"enabled": False, "reason": "insufficient_train_coverage"}
    medians = np.nanmedian(values[fit], axis=0)
    values = np.where(np.isfinite(values), values, medians[None, :]).astype(np.float32)
    data.clinical = TrainOnlyStandardizer().fit(values, fit).transform(values)
    return {"enabled": True, "n_present": int(present.sum()), "fit_patients": [data.patient_ids[i] for i in fit]}


class UncappedHoptimusBatches:
    """Load every canonical patch for a dynamic patient batch exactly once."""

    def __init__(self, data, indices: np.ndarray, token_budget: int, seed: int, shuffle: bool = True,
                 include_clinical: bool = False) -> None:
        if data.hoptimus_store is None:
            raise ValueError("V2 requires the canonical H-Optimus patch store")
        self.data, self.indices, self.seed, self.shuffle = data, np.asarray(indices, dtype=np.int64), int(seed), bool(shuffle)
        self.include_clinical = bool(include_clinical)
        index = data.hoptimus_store.index_path
        rows = __import__("pandas").read_parquet(index).set_index("patient_id")
        counts = [int(rows.loc[data.patient_ids[i], "n_tokens"]) for i in self.indices]
        self.sampler = DynamicTokenBatchSampler(counts, token_budget, seed)

    def __iter__(self) -> Iterator[dict[str, torch.Tensor]]:
        for local_group in self.sampler.batches(self.seed, shuffle=self.shuffle):
            group = self.indices[local_group]
            bags, slides, coords = [], [], []
            loaded = self.data.hoptimus_store.load_many_patient_tokens(
                [self.data.patient_ids[index] for index in group], max_tokens=None, seed=self.seed
            )
            for index, (features, metadata) in zip(group, loaded):
                if len(features) == 0:
                    raise RuntimeError(f"patient {self.data.patient_ids[index]} has no canonical tokens")
                codes, _ = __import__("pandas").factorize(metadata.slide_id.astype(str), sort=True)
                bags.append(features); slides.append(codes.astype(np.int64)); coords.append(_metadata_coordinates(metadata))
            width = max(len(bag) for bag in bags)
            patches = np.zeros((len(group), width, 1536), dtype=np.float32)
            mask = np.zeros((len(group), width), dtype=bool)
            slide_ids = np.zeros((len(group), width), dtype=np.int64)
            # Preserve spatial information for valid patients in a mixed
            # ragged batch.  A missing coordinate table is represented by an
            # explicit patient-level mask, never an artificial zero layout.
            coordinate_present = np.asarray([xy is not None for xy in coords], dtype=bool)
            coordinate_batch = np.zeros((len(group), width, 2), dtype=np.float32) if coordinate_present.any() else None
            for row, (bag, ids, xy) in enumerate(zip(bags, slides, coords)):
                patches[row, :len(bag)] = bag; mask[row, :len(bag)] = True
                slide_ids[row, :len(bag)] = ids
                if coordinate_batch is not None and xy is not None:
                    coordinate_batch[row, :len(bag)] = xy
            batch = {
                "indices": torch.from_numpy(group.copy()),
                "patches": torch.from_numpy(patches), "patch_mask": torch.from_numpy(mask),
                "slide_ids": torch.from_numpy(slide_ids),
                "rna": torch.from_numpy(self.data.rna[group].astype(np.float32)),
                "rna_present": torch.ones(len(group), dtype=torch.bool),
                "programme_target": torch.from_numpy(self.data._v2_programmes[group]),
                "programme_present": torch.from_numpy(self.data.hallmark_present[group].astype(bool)),
                "programme_target_mask": torch.from_numpy(np.isfinite(self.data._v2_programmes[group])),
                "programme_positive_mask": torch.from_numpy(self.data._v2_positive[np.ix_(group, group)]),
                "programme_neighbor_indices": torch.from_numpy(self.data._v2_neighbour_indices[group]),
            }
            if self.include_clinical:
                batch["clinical"] = torch.from_numpy(self.data.clinical[group].astype(np.float32))
                batch["clinical_present"] = torch.from_numpy(self.data.clinical_present[group].astype(bool))
            if coordinate_batch is not None:
                batch["coordinates"] = torch.from_numpy(coordinate_batch)
                batch["coordinate_present"] = torch.from_numpy(coordinate_present)
            for name in ("snv", "cnv"):
                if hasattr(self.data, f"_v2_{name}"):
                    batch[name] = torch.from_numpy(getattr(self.data, f"_v2_{name}")[group])
                    batch[f"{name}_present"] = torch.from_numpy(getattr(self.data, f"_v2_{name}_present")[group])
            if hasattr(self.data, "_mlp_clip_anchor_wsi"):
                batch["mlp_clip_anchor_wsi"] = torch.from_numpy(self.data._mlp_clip_anchor_wsi[group])
                batch["mlp_clip_anchor_rna"] = torch.from_numpy(self.data._mlp_clip_anchor_rna[group])
            if hasattr(self.data, "_plip_teacher"):
                batch["semantic_target"] = torch.from_numpy(self.data._plip_teacher[group])
                batch["semantic_present"] = torch.from_numpy(self.data._plip_present[group])
            yield batch

    def state_dict(self) -> dict[str, object]:
        return {"indices": self.indices.tolist(), "token_budget": self.sampler.token_budget, "seed": self.seed,
                "include_clinical": self.include_clinical, "coverage": self.sampler.coverage(self.seed).tolist()}


def attach_v2_targets(data, fit_mask: np.ndarray | None = None) -> dict[str, list[float]]:
    """Attach residual programme targets using the active fit population only."""
    train = np.asarray(data.split) == "train" if fit_mask is None else np.asarray(fit_mask, dtype=bool)
    if train.shape != (len(data.patient_ids),) or not train.any():
        raise ValueError("fit_mask must select at least one canonical patient")
    programmes, means = residualise_programmes(np.asarray(data.hallmark, dtype=np.float32), np.asarray(data.cancers), train)
    data._v2_programmes = programmes
    # Fixed top-k neighbours always supply a usable train-only positive graph;
    # unlike a hard similarity threshold it does not silently disappear for a
    # heterogeneous cancer cohort.  Positives may cross cancer type by design.
    normalized = programmes / np.maximum(np.linalg.norm(programmes, axis=1, keepdims=True), 1e-6)
    similarity = normalized @ normalized.T
    data._v2_positive = np.zeros_like(similarity, dtype=bool)
    data._v2_neighbour_indices = np.full((len(programmes), 8), -1, dtype=np.int64)
    train_rows = np.where(train & np.asarray(data.hallmark_present, dtype=bool))[0]
    for row in train_rows:
        candidates = train_rows[train_rows != row]
        if len(candidates):
            selected = candidates[np.argsort(-similarity[row, candidates])[: min(8, len(candidates))]]
            data._v2_positive[row, selected] = True
            data._v2_neighbour_indices[row, :len(selected)] = selected
    return means


def attach_mlp_clip_anchor(data, path: str) -> None:
    teacher = np.load(path, allow_pickle=False)
    required = {"patient_ids", "split"}
    wsi_key = "wsi" if "wsi" in teacher.files else "wsi_identity"
    rna_key = "rna" if "rna" in teacher.files else "rna_identity"
    if not required.issubset(teacher.files) or wsi_key not in teacher.files or rna_key not in teacher.files:
        raise ValueError("invalid MLP-CLIP teacher artifact")
    lookup = {str(pid): row for row, pid in enumerate(teacher["patient_ids"].astype(str))}
    rows = np.asarray([lookup.get(str(pid), -1) for pid in data.patient_ids])
    if (rows < 0).any() or teacher[wsi_key].shape[1] != 256 or teacher[rna_key].shape[1] != 256:
        raise ValueError("teacher does not exactly cover the paired cohort with 256-D embeddings")
    if not np.array_equal(teacher["split"][rows].astype(str), np.asarray(data.split).astype(str)):
        raise ValueError("teacher split labels differ from active split")
    data._mlp_clip_anchor_wsi = teacher[wsi_key][rows].astype(np.float32)
    data._mlp_clip_anchor_rna = teacher[rna_key][rows].astype(np.float32)
    data._mlp_clip_anchor_source = str(Path(path).resolve())


@torch.no_grad()
def _validation_selection(model: TumorStateV2, loader, device: str, cancers: np.ndarray) -> dict[str, float]:
    """Selection metrics evaluated only on development validation patients."""
    model.eval(); wsi, rna, target, prediction, present, labels = [], [], [], [], [], []
    for batch in loader:
        indices = batch.pop("indices").numpy()
        moved = {key: value.to(device, non_blocking=True) if isinstance(value, torch.Tensor) else value for key, value in batch.items()}
        out_wsi, out_rna = model(moved, "wsi"), model(moved, "rna")
        wsi.append(out_wsi["z_identity"].cpu().float().numpy()); rna.append(out_rna["z_identity"].cpu().float().numpy())
        target.append(moved["programme_target"].cpu().float().numpy()); prediction.append(out_wsi["programme_mean"].cpu().float().numpy())
        present.append(moved["programme_present"].cpu().numpy().astype(bool)); labels.extend(cancers[indices].astype(str).tolist())
    if not wsi:
        return {"retrieval_r10": float("nan"), "programme_mean_pearson": float("nan")}
    retrieval = paired_retrieval_metrics(np.vstack(wsi), np.vstack(rna), (10,), labels, labels)
    truth, pred, keep = np.vstack(target), np.vstack(prediction), np.concatenate(present)
    correlations = [np.corrcoef(truth[keep, j], pred[keep, j])[0, 1] for j in range(truth.shape[1])
                    if keep.sum() > 2 and np.std(truth[keep, j]) > 1e-8 and np.std(pred[keep, j]) > 1e-8]
    return {"retrieval_r10": float(retrieval.get("recall_at_10", float("nan"))),
            "programme_mean_pearson": float(np.nanmean(correlations)) if correlations else float("nan")}


def _load_slide_pretraining_checkpoint(model: TumorStateV2, path: str, split_digest: str, device: str) -> dict[str, object]:
    """Load an encoder-only checkpoint and reject foreign development splits."""
    state = torch.load(path, map_location=device, weights_only=False)
    manifest = state.get("manifest", {})
    if manifest.get("split_digest") and manifest["split_digest"] != split_digest:
        raise ValueError("slide pretraining checkpoint belongs to a different split")
    if "encoder" not in state:
        raise ValueError("slide pretraining checkpoint lacks encoder state")
    model.wsi.load_state_dict(state["encoder"], strict=True)
    return {"loaded": str(Path(path).resolve()), "epoch": int(state.get("epoch", -1)), "manifest": manifest}


def _run_slide_pretraining(model: TumorStateV2, data, train: np.ndarray, args: argparse.Namespace,
                           output: Path, split_manifest: dict[str, object], device: str) -> dict[str, object]:
    """Pretrain only the WSI aggregator on development-cancer bags.

    This stage intentionally creates no RNA/cancer-label loss and writes a
    standalone resume-friendly encoder checkpoint before fusion fine-tuning.
    """
    config = SlidePretrainingConfig(
        mask_fraction=args.pretrain_mask_fraction,
        view_keep_fraction=args.pretrain_view_keep_fraction,
        target_dim=args.pretrain_target_dim,
        seed=args.seed,
    )
    objective = SlidePretrainingObjective(model.config.patch_dim, model.config.hidden_dim, config).to(device)
    optimizer = torch.optim.AdamW(list(model.wsi.parameters()) + list(objective.parameters()),
                                  lr=args.pretrain_learning_rate, weight_decay=args.weight_decay)
    trainer = SlidePretrainer(model.wsi, objective, optimizer, device)
    manifest = {
        "stage": "train_cancer_only_slide_self_supervision",
        "split_digest": split_manifest["split_digest"],
        "fit_split": "train",
        "uses_rna": False, "uses_cancer_labels": False,
        "config": asdict(config),
        "source_manifest": source_manifest(configuration={"stage": "slide_pretraining", **asdict(config)}),
    }
    best = float("inf")
    metrics_path = output / "slide_pretrain_metrics.jsonl"
    for epoch in range(args.pretrain_epochs):
        loader = UncappedHoptimusBatches(data, train, args.token_budget, args.seed + epoch, include_clinical=False)
        metrics = trainer.train_epoch(loader)
        metrics_path.open("a", encoding="utf-8").write(json.dumps({"epoch": epoch, **metrics}) + "\n")
        checkpoint = {"version": 1, "epoch": epoch, "encoder": model.wsi.state_dict(),
                      "objective": objective.state_dict(), "optimizer": optimizer.state_dict(), "manifest": manifest}
        torch.save(checkpoint, output / "slide_pretrain_last.pt")
        if metrics.get("pretrain_loss", float("inf")) < best:
            best = metrics["pretrain_loss"]
            torch.save(checkpoint, output / "slide_pretrain_best.pt")
    return {"trained": True, "epochs": args.pretrain_epochs, "checkpoint": str((output / "slide_pretrain_best.pt").resolve()),
            "best_loss": best, "manifest": manifest}


def _selected_checkpoint(output: Path, slide_pretraining: dict[str, object], fine_tune_epochs: int,
                         *, allow_last: bool = False) -> Path:
    """Resolve the terminal artifact for both pretrain-only and full runs."""
    if fine_tune_epochs == 0 and slide_pretraining.get("trained"):
        checkpoint = Path(str(slide_pretraining["checkpoint"]))
        if checkpoint.exists():
            return checkpoint
        raise RuntimeError("slide pretraining reported a checkpoint that was not written")
    if (output / "best_pareto.pt").exists():
        return output / "best_pareto.pt"
    if (output / "best_retrieval.pt").exists():
        return output / "best_retrieval.pt"
    if allow_last and (output / "last.pt").exists():
        # A final development refit intentionally has no validation patients;
        # its epoch count was fixed by the preceding inner-validation stage.
        return output / "last.pt"
    raise RuntimeError("no validation checkpoint was selected; refuse to export an unselected last checkpoint")


def _trained_states_for_profile(objective_profile: str, semantic_dim: int = 0) -> list[str]:
    """Declare only representations with an explicit gradient path.

    The artifact contract is deliberately stricter than tensor availability:
    `TumorStateV2` can compute every view, but diagnostics must never export a
    state that its selected objective did not train.
    """
    if objective_profile == "identity_only":
        return ["wsi_identity", "rna_identity", "full_identity"]
    if objective_profile == "programme_only":
        return ["wsi_biology", "rna_biology", "full_biology"]
    if objective_profile != "full":
        raise ValueError(f"unknown objective profile {objective_profile!r}")
    states = ["wsi_identity", "rna_identity", "full_identity", "wsi_biology", "rna_biology", "full_biology", "full_patient"]
    if semantic_dim:
        states.append("semantic_wsi")
    return states


def run(args: argparse.Namespace) -> Path:
    if args.gradient_diagnostics_every < 0:
        raise ValueError("--gradient-diagnostics-every must be non-negative")
    if args.pretrain_epochs < 0:
        raise ValueError("--pretrain-epochs must be non-negative")
    if args.pretrain_epochs and args.pretrain_checkpoint:
        raise ValueError("choose either --pretrain-epochs or --pretrain-checkpoint, not both")
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    data = load_bio_query_data(args.data_config, args.split_file, wsi_mode="hoptimus_patch")
    split_manifest = validate_runtime_split(args.split_file, data.patient_ids, data.cancers, data.split,
                                            args.expected_development_cancers, args.expected_heldout_cancers)
    fit_mask = np.asarray(data.split).astype(str) != "test" if args.fit_development else np.asarray(data.split).astype(str) == "train"
    clinical_manifest = (_standardize_clinical(data, fit_mask) if args.include_clinical
                         else {"enabled": False, "reason": "strict_core_wsi_rna"})
    means = attach_v2_targets(data, fit_mask)
    optional_manifest = {"clinical": clinical_manifest}
    for name, path in (("snv", args.snv_features), ("cnv", args.cnv_features)):
        if path:
            optional_manifest[name] = _attach_numeric_table(data, path, name)
    anchor_path = args.mlp_clip_anchor or args.mlp_clip_teacher
    use_anchor = bool(anchor_path) and args.objective_profile != "programme_only"
    if anchor_path and not use_anchor:
        warnings.warn("programme_only ignores the MLP-CLIP anchor by design", RuntimeWarning)
    if use_anchor:
        attach_mlp_clip_anchor(data, anchor_path)
    semantic_dim = 0
    if args.plip_teacher:
        source_digests = canonical_patient_source_digests(data.hoptimus_store._metadata())
        data._plip_teacher, data._plip_present, plip_manifest = load_plip_teacher_cache(args.plip_teacher, data.patient_ids, source_digests)
        semantic_dim = int(data._plip_teacher.shape[1])
        optional_manifest["plip_semantic"] = plip_manifest
    train = np.where(fit_mask)[0]
    validation = np.asarray([], dtype=np.int64) if args.fit_development else np.where(np.asarray(data.split) == "val")[0]
    device = "cuda" if args.device == "auto" and torch.cuda.is_available() else args.device
    config = V2ModelConfig(rna_dim=int(data.rna.shape[1]), hidden_dim=args.hidden_dim, layers=args.layers, heads=args.heads,
                           semantic_dim=semantic_dim, use_mlp_clip_anchor=use_anchor)
    model = TumorStateV2(config, clinical_dim=int(data.clinical.shape[1]) if args.include_clinical else 0,
                         snv_dim=None if not hasattr(data, "_v2_snv") else int(data._v2_snv.shape[1]),
                         cnv_dim=None if not hasattr(data, "_v2_cnv") else int(data._v2_cnv.shape[1]),
                         programme_dim=int(data.hallmark.shape[1])).to(device)
    output = Path(args.output_dir); output.mkdir(parents=True, exist_ok=True)
    if args.pretrain_checkpoint:
        slide_pretraining = _load_slide_pretraining_checkpoint(
            model, args.pretrain_checkpoint, split_manifest["split_digest"], device
        )
    elif args.pretrain_epochs:
        slide_pretraining = _run_slide_pretraining(model, data, train, args, output, split_manifest, device)
    else:
        slide_pretraining = {"trained": False, "reason": "not_requested"}
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=max(args.epochs, 1))
    if args.teacher_warmup_epochs:
        warnings.warn(
            "--teacher-warmup-epochs is deprecated and ignored: V2.1 uses an explicit "
            "baseline-preserving anchor residual, not a separate teacher loss.", RuntimeWarning,
        )
    trainer = V2Trainer(
        model, optimizer, V2LossSchedule(objective_profile=args.objective_profile), device,
        gradient_diagnostics_every=args.gradient_diagnostics_every,
    )
    trained_states = _trained_states_for_profile(args.objective_profile, semantic_dim)
    run_configuration = {
        "seed": args.seed, "epochs": args.epochs, "token_budget": args.token_budget,
        "hidden_dim": args.hidden_dim, "layers": args.layers, "heads": args.heads,
        "learning_rate": args.learning_rate, "weight_decay": args.weight_decay,
        "gradient_diagnostics_every": args.gradient_diagnostics_every,
        "objective_profile": args.objective_profile,
        "pretrain_epochs": args.pretrain_epochs, "pretrain_checkpoint": args.pretrain_checkpoint,
        "anchor_path": anchor_path if use_anchor else "", "anchor_requested_but_disabled": bool(anchor_path) and not use_anchor,
        "fit_population": "development_train_val" if args.fit_development else "train_only_inner_fit",
        "include_clinical": args.include_clinical,
        "snv_features": args.snv_features, "cnv_features": args.cnv_features,
    }
    manifest = {"artifact_version": 4, "protocol": split_manifest["protocol"], "patch_cap": None, "patient_ids": list(data.patient_ids), "cancers": list(data.cancers), "split": list(data.split), "residual_means_fit_on": "train_only", "residual_means": means,
                "trained_states": list(trainable_state_names(*trained_states)), "selection_metric": "validation_pareto",
                "coordinates": "metadata_validated_or_disabled", "semantic_dim": semantic_dim,
                "modal_dims": {"clinical": int(data.clinical.shape[1]) if args.include_clinical else 0, "snv": 0 if not hasattr(data, "_v2_snv") else int(data._v2_snv.shape[1]), "cnv": 0 if not hasattr(data, "_v2_cnv") else int(data._v2_cnv.shape[1])},
                "optional_modalities": optional_manifest, "mlp_clip_anchor": use_anchor,
                "teacher_mode": "anchor_residual_refiner" if use_anchor else "none",
                "anchor_teacher_source": getattr(data, "_mlp_clip_anchor_source", None) if use_anchor else None,
                "model_config": asdict(config), "split_manifest": split_manifest,
                "strict_core": not args.include_clinical and not bool(args.snv_features) and not bool(args.cnv_features),
                "fit_population": "development_train_val" if args.fit_development else "train_only_inner_fit",
                "slide_pretraining": slide_pretraining,
                "state_activation_epoch": {name: 0 for name in trained_states},
                "objective_profile": args.objective_profile,
                "source_manifest": source_manifest(configuration=run_configuration),
                "run_configuration": run_configuration}
    start_epoch = 0
    if args.resume:
        restored = trainer.load_checkpoint(args.resume, scheduler)
        prior_digest = restored.get("manifest", {}).get("split_manifest", {}).get("split_digest")
        if prior_digest and prior_digest != split_manifest["split_digest"]:
            raise ValueError("resume checkpoint belongs to a different split")
        start_epoch = int(restored["epoch"]) + 1
    anchor_validation_r10 = float("nan")
    if use_anchor and len(validation):
        labels = np.asarray(data.cancers)[validation].tolist()
        anchor_validation_r10 = float(paired_retrieval_metrics(data._mlp_clip_anchor_wsi[validation], data._mlp_clip_anchor_rna[validation], (10,), labels, labels).get("recall_at_10", float("nan")))
    best_retrieval = best_programme = best_pareto = -float("inf")
    best_retrieval_epoch = best_programme_epoch = best_pareto_epoch = None
    for epoch in range(start_epoch, args.epochs):
        loader = UncappedHoptimusBatches(data, train, args.token_budget, args.seed + epoch, include_clinical=args.include_clinical)
        metrics = trainer.train_epoch(loader, epoch)
        validation_metrics = {}; selection_metrics = {}
        if len(validation):
            validation_loader = UncappedHoptimusBatches(data, validation, args.token_budget, args.seed, shuffle=False, include_clinical=args.include_clinical)
            validation_metrics = trainer.evaluate_epoch(validation_loader, epoch)
            selection_loader = UncappedHoptimusBatches(data, validation, args.token_budget, args.seed, shuffle=False, include_clinical=args.include_clinical)
            selection_metrics = _validation_selection(model, selection_loader, device, np.asarray(data.cancers))
        scheduler.step()
        row = {"epoch": epoch, "learning_rate": scheduler.get_last_lr()[0], "anchor_validation_r10": anchor_validation_r10, **{f"train_{k}": v for k, v in metrics.items()}, **{f"validation_{k}": v for k, v in validation_metrics.items()}, **{f"selection_{k}": v for k, v in selection_metrics.items()}}
        (output / "train_metrics.jsonl").open("a", encoding="utf-8").write(json.dumps(row) + "\n")
        trainer.save_checkpoint(output / "last.pt", epoch, loader.state_dict(), manifest, scheduler)
        r10, programme = float(selection_metrics.get("retrieval_r10", float("nan"))), float(selection_metrics.get("programme_mean_pearson", float("nan")))
        if np.isfinite(r10) and r10 > best_retrieval:
            best_retrieval = r10; best_retrieval_epoch = epoch; trainer.save_checkpoint(output / "best_retrieval.pt", epoch, loader.state_dict(), manifest, scheduler)
        if np.isfinite(programme) and programme > best_programme:
            best_programme = programme; best_programme_epoch = epoch; trainer.save_checkpoint(output / "best_programme.pt", epoch, loader.state_dict(), manifest, scheduler)
        eligible = np.isfinite(r10) and (not np.isfinite(anchor_validation_r10) or r10 >= 0.95 * anchor_validation_r10)
        if eligible and np.isfinite(programme) and programme > best_pareto:
            best_pareto = programme; best_pareto_epoch = epoch; trainer.save_checkpoint(output / "best_pareto.pt", epoch, loader.state_dict(), manifest, scheduler)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (output / "selection.json").write_text(json.dumps({"anchor_validation_r10": anchor_validation_r10, "best_retrieval_r10": best_retrieval, "best_programme_pearson": best_programme, "best_pareto_programme_pearson": best_pareto, "best_retrieval_epoch": best_retrieval_epoch, "best_programme_epoch": best_programme_epoch, "best_pareto_epoch": best_pareto_epoch, "fit_population": "development_train_val" if args.fit_development else "train_only_inner_fit"}, indent=2), encoding="utf-8")
    return _selected_checkpoint(output, slide_pretraining, args.epochs, allow_last=args.fit_development)


def main() -> None:
    parser = argparse.ArgumentParser(description="Train MORPHEUS V2 with uncapped canonical H-Optimus bags")
    parser.add_argument("--data-config", default="morpheus/configs/v1.json")
    parser.add_argument("--split-file", required=True)
    parser.add_argument("--output-dir", default="morpheus/outputs/v2")
    parser.add_argument("--epochs", type=int, default=40); parser.add_argument("--token-budget", type=int, default=32768)
    parser.add_argument("--hidden-dim", type=int, default=512); parser.add_argument("--layers", type=int, default=4); parser.add_argument("--heads", type=int, default=8)
    parser.add_argument("--learning-rate", type=float, default=2e-4); parser.add_argument("--weight-decay", type=float, default=1e-2)
    parser.add_argument("--seed", type=int, default=42); parser.add_argument("--device", default="auto")
    parser.add_argument("--mlp-clip-teacher", default="")
    parser.add_argument("--mlp-clip-anchor", default="", help="paired frozen all-patch MLP-CLIP artifact")
    parser.add_argument("--teacher-warmup-epochs", type=int, default=0,
                        help="deprecated compatibility flag; V2.1 uses an explicit anchor residual instead")
    parser.add_argument("--gradient-diagnostics-every", type=int, default=25,
                        help="log objective-gradient cosines every N train batches; 0 disables")
    parser.add_argument("--objective-profile", choices=("full", "identity_only", "programme_only"), default="full",
                        help="diagnostic objective ablation; preserves data/model while selectively zeroing losses")
    parser.add_argument("--pretrain-epochs", type=int, default=0,
                        help="development-train-only slide self-supervision epochs before V2 fine-tuning")
    parser.add_argument("--pretrain-checkpoint", default="",
                        help="previous development-split slide-pretraining checkpoint to load before fine-tuning")
    parser.add_argument("--pretrain-learning-rate", type=float, default=2e-4)
    parser.add_argument("--pretrain-mask-fraction", type=float, default=0.30)
    parser.add_argument("--pretrain-view-keep-fraction", type=float, default=0.70)
    parser.add_argument("--pretrain-target-dim", type=int, default=128)
    parser.add_argument("--snv-features", default="")
    parser.add_argument("--cnv-features", default="")
    parser.add_argument("--plip-teacher", default="", help="checksum-audited frozen PLIP patient target cache")
    parser.add_argument("--include-clinical", action="store_true", help="enable audited clinical tokens; strict core leaves them absent")
    parser.add_argument("--resume", default="", help="resume a V2 checkpoint with matching split provenance")
    parser.add_argument("--fit-development", action="store_true",
                        help="final refit on train+validation only; epoch count must already be selected")
    parser.add_argument("--expected-development-cancers", type=int, default=11)
    parser.add_argument("--expected-heldout-cancers", type=int, default=22)
    print(run(parser.parse_args()))


if __name__ == "__main__":
    main()
