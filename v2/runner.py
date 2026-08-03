"""Canonical uncapped TCGA V2 training path.

This module deliberately reads the canonical H-Optimus store directly.  It
does not call the legacy capped patch-batch helper.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import warnings
from dataclasses import asdict, replace
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import torch

from morpheus.src.training.train_bio_query_former import load_bio_query_data

from .data import DynamicTokenBatchSampler, TrainOnlyStandardizer
from .model import TumorStateV2, V2ModelConfig
from .training import PairedBiologyMemoryBank, V2LossSchedule, V2Trainer
from .contracts import trainable_state_names
from .plip import canonical_patient_source_digests, load_plip_teacher_cache
from .preflight import restrict_cohort_to_split, validate_runtime_split
from .provenance import source_manifest
from .slide_pretraining import SlidePretrainingConfig, SlidePretrainingObjective, SlidePretrainer
from .pbs import LegibilityOperator
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
            targets = np.asarray(self.data._v2_programmes[group], dtype=np.float32)
            head_dim = int(getattr(self.data, "_v2_programme_head_dim", targets.shape[1]))
            if head_dim < targets.shape[1]:
                raise ValueError("programme head dimension cannot truncate frozen supervision coordinates")
            if head_dim != targets.shape[1]:
                padded = np.full((len(group), head_dim), np.nan, dtype=np.float32)
                padded[:, :targets.shape[1]] = targets
                targets = padded
            batch = {
                "indices": torch.from_numpy(group.copy()),
                "patches": torch.from_numpy(patches), "patch_mask": torch.from_numpy(mask),
                "slide_ids": torch.from_numpy(slide_ids),
                "rna": torch.from_numpy(self.data.rna[group].astype(np.float32)),
                "rna_present": torch.ones(len(group), dtype=torch.bool),
                "programme_target": torch.from_numpy(targets),
                "programme_present": torch.from_numpy(self.data._v2_programme_present[group].astype(bool)),
                "programme_target_mask": torch.from_numpy(np.isfinite(targets)),
                "programme_positive_mask": torch.from_numpy(self.data._v2_positive[np.ix_(group, group)]),
                "programme_neighbor_indices": torch.from_numpy(self.data._v2_neighbour_indices[group]),
            }
            if hasattr(self.data, "_v2_programme_axis_weights"):
                weights = self.data._v2_programme_axis_weights.astype(np.float32)
                if head_dim != len(weights):
                    weights = np.pad(weights, (0, head_dim - len(weights)))
                batch["programme_axis_weights"] = torch.from_numpy(weights)
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


def _attach_programme_matrix(data, values: np.ndarray, present: np.ndarray,
                             fit_mask: np.ndarray | None = None) -> dict[str, list[float]]:
    """Attach a finite programme matrix and its train-only positive graph."""
    train = np.asarray(data.split) == "train" if fit_mask is None else np.asarray(fit_mask, dtype=bool)
    if train.shape != (len(data.patient_ids),) or not train.any():
        raise ValueError("fit_mask must select at least one canonical patient")
    values, present = np.asarray(values, dtype=np.float32), np.asarray(present, dtype=bool)
    if values.ndim != 2 or values.shape[0] != len(data.patient_ids) or present.shape != (len(data.patient_ids),):
        raise ValueError("programme targets must be [canonical_patient, target] with a row presence mask")
    if not np.isfinite(values[present]).all():
        raise ValueError("present programme rows must be finite")
    if not present[train].all():
        raise ValueError("programme supervision must cover every development fitting patient")
    # Residualisation must never silently turn missing rows into a numerical
    # target.  D2's canonical paired RNA matrix is required to cover all rows;
    # the explicit branch retains the invariant for future optional targets.
    filled = np.where(np.isfinite(values), values, 0.0)
    programmes, means = residualise_programmes(filled, np.asarray(data.cancers), train)
    programmes[~present] = np.nan
    development_std = np.nanstd(programmes[train & present], axis=0)
    dead = np.where(~np.isfinite(development_std) | (development_std < 1e-8))[0]
    if len(dead):
        raise ValueError(f"programme targets contain development-constant axes; examples={dead[:10].tolist()}")
    data._v2_programmes = programmes
    data._v2_programme_present = present
    # Fixed top-k neighbours always supply a usable train-only positive graph;
    # unlike a hard similarity threshold it does not silently disappear for a
    # heterogeneous cancer cohort.  Positives may cross cancer type by design.
    normalized = programmes / np.maximum(np.linalg.norm(programmes, axis=1, keepdims=True), 1e-6)
    similarity = normalized @ normalized.T
    data._v2_positive = np.zeros_like(similarity, dtype=bool)
    data._v2_neighbour_indices = np.full((len(programmes), 8), -1, dtype=np.int64)
    train_rows = np.where(train & present)[0]
    for row in train_rows:
        candidates = train_rows[train_rows != row]
        if len(candidates):
            selected = candidates[np.argsort(-similarity[row, candidates])[: min(8, len(candidates))]]
            data._v2_positive[row, selected] = True
            data._v2_neighbour_indices[row, :len(selected)] = selected
    return means


def fit_programme_legibility_operator(data, fit_mask: np.ndarray) -> dict[str, object]:
    """Fit the same grouped-CV quotient-axis operator for H and PBS targets.

    The operation is deliberately target-agnostic: D2 uses it in *both* arms,
    so the H-versus-I comparison changes only the target coordinate system.
    It is fit on development patients and frozen before model optimisation.
    """
    fit = np.asarray(fit_mask, dtype=bool)
    operator = LegibilityOperator.fit(data.wsi_patient[fit], data._v2_programmes[fit],
                                      np.asarray(data.cancers)[fit])
    weights = np.asarray(operator.weights, dtype=np.float32)
    if not np.isfinite(weights).all() or not (weights > 1e-8).any():
        raise ValueError("programme legibility operator has no nonzero development-fitted axes")
    data._v2_programme_axis_weights = weights
    identifiers = np.asarray(data.patient_ids).astype(str)[fit]
    return {"fit_population": "development_only", "fit_patient_id_digest": __import__("hashlib").sha256(
                "\n".join(identifiers).encode("utf-8")).hexdigest(),
            "split_digest": __import__("hashlib").sha256("\n".join(np.asarray(data.split).astype(str)[fit]).encode("utf-8")).hexdigest(),
            "operator": "diagonal_grouped_ridge", "alpha": operator.alpha,
            "nonzero_axes": int((weights > 1e-8).sum()), "weights": weights.tolist(),
            "fold_score": operator.fold_score.tolist()}


def attach_v2_targets(data, fit_mask: np.ndarray | None = None) -> dict[str, list[float]]:
    """Attach Hallmark targets using the active fit population only."""
    return _attach_programme_matrix(data, np.asarray(data.hallmark, dtype=np.float32),
                                    np.asarray(data.hallmark_present, dtype=bool), fit_mask)


def attach_external_programme_targets(data, path: str, fit_mask: np.ndarray | None = None) -> dict[str, object]:
    """Load immutable D2 targets by patient ID, never by positional row.

    The builder writes the raw intervention coordinates for every canonical
    paired patient.  This adapter performs the *only* cohort-dependent fit
    (cancer residualisation and neighbour graph) on development patients.
    """
    raw = np.load(path, allow_pickle=False)
    required = {"patient_ids", "split", "scores", "target_names", "target_groups", "genes", "singular_values",
                "gene_basis", "atom_coordinates", "atom_ids", "manifest_json"}
    if not required.issubset(raw.files):
        raise ValueError("programme target NPZ needs patient_ids, scores, and target_names")
    try:
        manifest = json.loads(str(np.asarray(raw["manifest_json"]).item()))
    except (ValueError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("programme target NPZ has an unreadable manifest_json") from exc
    if manifest.get("target_kind") != "external_perturbation_dictionary_coordinates":
        raise ValueError("programme target NPZ is not a PBS intervention-dictionary artifact")
    patient_ids = np.asarray(raw["patient_ids"]).astype(str)
    scores = np.asarray(raw["scores"], dtype=np.float32)
    names = np.asarray(raw["target_names"]).astype(str)
    groups = np.asarray(raw["target_groups"]).astype(str)
    genes = np.asarray(raw["genes"]).astype(str)
    singular_values = np.asarray(raw["singular_values"], dtype=np.float32)
    gene_basis = np.asarray(raw["gene_basis"], dtype=np.float32)
    atom_coordinates = np.asarray(raw["atom_coordinates"], dtype=np.float32)
    atom_ids = np.asarray(raw["atom_ids"]).astype(str)
    if scores.ndim != 2 or scores.shape != (len(patient_ids), len(names)) or len(set(patient_ids)) != len(patient_ids):
        raise ValueError("programme target NPZ has invalid dimensions or duplicate patient IDs")
    expected_names = np.asarray([f"PBS_{index:03d}" for index in range(scores.shape[1])])
    if (scores.shape[1] not in {64, 128, 256} or not np.array_equal(names, expected_names)
            or groups.shape != names.shape or not np.all(groups == "PBS")
            or singular_values.shape != names.shape or genes.ndim != 1 or len(set(genes)) != len(genes)
            or gene_basis.shape != (len(genes), len(names))
            or atom_coordinates.shape != (len(atom_ids), len(names))
            or len(atom_ids) == 0 or len(set(atom_ids)) != len(atom_ids)):
        raise ValueError("programme target NPZ violates the frozen PBS schema")
    digest_strings = lambda values: __import__("hashlib").sha256("\n".join(map(str, values)).encode("utf-8")).hexdigest()
    digest_array = lambda values: __import__("hashlib").sha256(np.ascontiguousarray(values).view(np.uint8)).hexdigest()
    if (manifest.get("patient_id_digest") != digest_strings(patient_ids)
            or manifest.get("split_digest") != digest_strings(np.asarray(raw["split"]).astype(str))
            or manifest.get("overlap_gene_digest") != digest_strings(genes)
            or manifest.get("scores_sha256") != digest_array(scores)
            or manifest.get("singular_values_sha256") != digest_array(singular_values)
            or manifest.get("gene_basis_sha256") != digest_array(gene_basis)
            or manifest.get("atom_coordinates_sha256") != digest_array(atom_coordinates)
            or manifest.get("atom_id_digest") != digest_strings(atom_ids)):
        raise ValueError("programme target NPZ provenance digest mismatch")
    fit = np.asarray(data.split).astype(str) != "test" if fit_mask is None else np.asarray(fit_mask, dtype=bool)
    expected_fit_digest = digest_strings(np.asarray(data.patient_ids).astype(str)[fit])
    if manifest.get("fit_patient_id_digest") != expected_fit_digest:
        raise ValueError("PBS target transform was fit on a different patient population than this runner")
    lookup = {patient: row for row, patient in enumerate(patient_ids)}
    rows = np.asarray([lookup.get(str(patient), -1) for patient in data.patient_ids], dtype=np.int64)
    if (rows < 0).any():
        missing = [str(data.patient_ids[index]) for index in np.where(rows < 0)[0][:5]]
        raise ValueError(f"programme target NPZ does not cover canonical cohort; examples={missing}")
    if "split" not in raw.files or not np.array_equal(np.asarray(raw["split"])[rows].astype(str), np.asarray(data.split).astype(str)):
        raise ValueError("programme target split labels differ from active split")
    values = scores[rows]
    present = np.isfinite(values).all(axis=1)
    means = _attach_programme_matrix(data, values, present, fit_mask)
    data._v2_programme_target_names = names.tolist()
    return {"source": str(Path(path).resolve()), "target_names": names.tolist(),
            "target_dimension": int(len(names)), "fit_residual_means": means,
            "coverage": int(present.sum()), "manifest": manifest,
            "axis_annotation_status": manifest.get("axis_annotations", {}).get("status", "unavailable_missing_axis_annotation_manifest")}


def _validate_d2_pair(args: argparse.Namespace) -> dict[str, object] | None:
    """Fail closed unless H and PBS runs were created from one paired manifest.

    D2 is causal only when *all* training choices other than target coordinates
    are identical.  The manifest is generated by ``phase_d.py d2`` and each
    runner process independently verifies its own invocation against it.
    """
    if not args.d2_pair_manifest:
        if args.d2_arm:
            raise ValueError("--d2-arm requires --d2-pair-manifest")
        return None
    if args.d2_arm not in {"H", "I"}:
        raise ValueError("--d2-pair-manifest requires --d2-arm H or I")
    path = Path(args.d2_pair_manifest)
    try:
        pair = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("D2 pair manifest is unreadable") from exc
    if pair.get("schema_version") != 1 or pair.get("experiment") != "D2_H_vs_I":
        raise ValueError("not a recognised D2 H-vs-I pair manifest")
    expected = pair.get("common_args")
    if not isinstance(expected, dict):
        raise ValueError("D2 pair manifest is missing common_args")
    current = _d2_common_args(args)
    if set(expected) != set(current):
        missing, extra = sorted(set(current) - set(expected)), sorted(set(expected) - set(current))
        raise ValueError(f"D2 pair manifest does not bind the exhaustive training configuration; missing={missing}, extra={extra}")
    if expected != current:
        changed = [name for name in expected if expected[name] != current[name]]
        raise ValueError(f"D2 paired-arm asymmetry in immutable configuration: {changed}")
    expected_digest = hashlib.sha256(json.dumps(expected, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    if pair.get("common_config_sha256") != expected_digest:
        raise ValueError("D2 pair manifest common configuration digest is invalid")
    if int(args.seed) not in {int(seed) for seed in pair.get("seeds", [])}:
        raise ValueError("D2 seed was not predeclared in the paired manifest")
    expected_target = str(pair.get("targets", {}).get(args.d2_arm, ""))
    actual_target = str(Path(args.programme_targets).resolve()) if args.programme_targets else ""
    if args.d2_arm == "H":
        if expected_target or actual_target:
            raise ValueError("D2 Hallmark arm must have no external programme target artifact")
    elif not expected_target or actual_target != str(Path(expected_target).resolve()):
        raise ValueError("D2 PBS arm target does not match the paired manifest")
    if args.d2_arm == "I" and hashlib.sha256(Path(actual_target).read_bytes()).hexdigest() != pair.get("pbs_target_sha256"):
        raise ValueError("D2 PBS target bytes differ from the immutable paired-manifest digest")
    if args.d2_arm == "I":
        with np.load(actual_target, allow_pickle=False) as raw:
            target_manifest = json.loads(str(np.asarray(raw["manifest_json"]).item()))
        if int(target_manifest.get("n_components", -1)) != int(pair.get("pbs_components", -1)):
            raise ValueError("D2 PBS component count differs from the paired manifest")
    if args.d2_analysis_role != pair.get("analysis_role") or int(args.d2_pbs_components) != int(pair.get("pbs_components", -1)):
        raise ValueError("D2 primary/sensitivity role differs from the immutable paired manifest")
    # Only the source target is allowed to differ between arms; whether to fit
    # the grouped-CV loss metric is itself a scientific choice and is frozen.
    if bool(args.fit_programme_legibility) != bool(pair.get("fit_programme_legibility")):
        raise ValueError("D2 arms disagree with the frozen legibility-operator policy")
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "arm": args.d2_arm, "fit_programme_legibility": bool(args.fit_programme_legibility)}


def _d2_common_args(args: argparse.Namespace) -> dict[str, object]:
    """Every model- or data-changing runner setting, except arm target/output/seed.

    D2's interpretation rests on one invariant: supervision coordinates are the
    only between-arm difference.  This deliberately enumerates even settings
    currently at their defaults so a future CLI option cannot escape the pair
    contract silently.
    """
    file_digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return {
        "data_config": str(Path(args.data_config).resolve()), "data_config_sha256": file_digest(args.data_config),
        "split_file": str(Path(args.split_file).resolve()), "split_file_sha256": file_digest(args.split_file),
        "epochs": int(args.epochs), "token_budget": int(args.token_budget), "hidden_dim": int(args.hidden_dim),
        "layers": int(args.layers), "heads": int(args.heads), "learning_rate": float(args.learning_rate),
        "weight_decay": float(args.weight_decay), "device": str(args.device),
        "mlp_clip_teacher": str(args.mlp_clip_teacher), "mlp_clip_anchor": str(args.mlp_clip_anchor),
        "teacher_warmup_epochs": int(args.teacher_warmup_epochs),
        "gradient_diagnostics_every": int(args.gradient_diagnostics_every),
        "objective_profile": str(args.objective_profile), "decorrelation_weight": float(args.decorrelation_weight),
        "loss_warmup_epochs": int(args.loss_warmup_epochs), "programme_warmup_weight": float(args.programme_warmup_weight),
        "programme_weight": float(args.programme_weight), "programme_neighbourhood_weight": float(args.programme_neighbourhood_weight),
        "programme_supcon_weight": float(args.programme_supcon_weight), "separation_weight": float(args.separation_weight),
        "variance_weight": float(args.variance_weight), "pretrain_epochs": int(args.pretrain_epochs),
        "programme_head_dim": int(args.programme_head_dim),
        "pretrain_checkpoint": str(args.pretrain_checkpoint), "pretrain_learning_rate": float(args.pretrain_learning_rate),
        "pretrain_mask_fraction": float(args.pretrain_mask_fraction), "pretrain_view_keep_fraction": float(args.pretrain_view_keep_fraction),
        "pretrain_target_dim": int(args.pretrain_target_dim), "snv_features": str(args.snv_features),
        "cnv_features": str(args.cnv_features), "plip_teacher": str(args.plip_teacher),
        "include_clinical": bool(args.include_clinical), "resume": str(args.resume),
        "fit_development": bool(args.fit_development), "fixed_final_epoch": bool(args.fixed_final_epoch),
        "expected_development_cancers": int(args.expected_development_cancers),
        "expected_heldout_cancers": int(args.expected_heldout_cancers),
        "fit_programme_legibility": bool(args.fit_programme_legibility),
        # The cohort is a data-changing setting: two arms trained on different
        # patient sets are not a target-only ablation, so it is bound here too.
        "restrict_to_split": bool(args.restrict_to_split),
        "d2_analysis_role": str(args.d2_analysis_role), "d2_pbs_components": int(args.d2_pbs_components),
    }


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


def _parameter_snapshot(trainer: V2Trainer) -> dict[str, list[torch.Tensor]]:
    return {name: [parameter.detach().cpu().float().clone() for parameter in parameters]
            for name, parameters in trainer.liveness_parameter_groups().items()}


def _parameter_deltas(trainer: V2Trainer, initial: dict[str, list[torch.Tensor]]) -> dict[str, float]:
    result: dict[str, float] = {}
    for name, parameters in trainer.liveness_parameter_groups().items():
        before = initial[name]
        numerator = sum(((parameter.detach().cpu().float() - reference).square().sum()
                         for parameter, reference in zip(parameters, before)), torch.zeros(()))
        denominator = sum((reference.square().sum() for reference in before), torch.zeros(())).clamp_min(1e-12)
        result[name] = float(torch.sqrt(numerator / denominator))
    return result


def _truncate_batch(batch: dict, limit: int) -> dict:
    """Shrink a training batch to `limit` patients for the overfit liveness test.

    "Overfit one batch" is only decisive on a SMALL batch. The token sampler
    hands out up to ~280 patients here (most patients carry ~30 tokens), and
    driving a heteroscedastic NLL over 280 x 50 targets to near zero in 300
    steps is not achievable even for a healthy model -- the real run measured a
    0.306 reduction against a 0.80 gate. Truncating makes the test strictly
    harder to pass spuriously: a model that cannot memorise 16 patients is
    broken beyond argument.
    """
    import torch as _torch

    size = len(batch["indices"])
    if limit >= size:
        return batch
    out = {}
    for key, value in batch.items():
        if isinstance(value, _torch.Tensor) and value.shape[:1] == (size,):
            value = value[:limit]
            if value.ndim >= 2 and value.shape[1] == size:  # [patient, patient] positive masks
                value = value[:, :limit]
        out[key] = value
    return out


def _overfit_programme_only_actual(model: TumorStateV2, schedule: V2LossSchedule,
                                   loader: UncappedHoptimusBatches, device: str,
                                   steps: int = 800, overfit_patients: int = 16) -> dict[str, object]:
    """G2.6 for the Hallmark D1 arm on the real model/trainer path.

    Training copied heads on frozen features is not a liveness test for the
    slide encoder.  This clone instead runs the exact three-view programme
    objective through `V2Trainer.step`, including the WSI/RNA/shared/biology
    groups that D1 compares.  The decorrelation term is omitted here only
    because it has a batch-statistics floor unrelated to memorisation.
    """
    # lr=3e-4. Swept against a stalling seed on the real runner at 800 steps:
    #   3e-3 FAIL (reduction 0.525)  1e-3 FAIL (flat, best at step 0)  3e-4 PASS
    # The direction is one-sided: lower is more stable. Same pathology as the
    # divergence recorded below, which was measured at a different seed.
    # Originally 1e-2. Measured on the real runner at 1e-2 the memorisation
    # DESCENDS then DIVERGES inside the window the gate reads:
    #   step 0 2.4154 -> 150 1.2822 -> 250 0.9432 (best) -> 299 2.1575
    # and by ~2000 steps it reaches nan outright. The heteroscedastic Gaussian
    # NLL's variance head collapses under AdamW at that rate, so the gate was
    # reading a post-divergence value and calling a healthy model dead. Both D1
    # arms use the same setting: their liveness checks are only comparable
    # evidence if they are run identically.
    clone = copy.deepcopy(model).to(device)
    clone_schedule = replace(schedule, decorrelation_after_warmup=0.0)
    optimiser = torch.optim.AdamW(clone.parameters(), lr=3e-4, weight_decay=0.0)
    trial = V2Trainer(clone, optimiser, clone_schedule, device, gradient_diagnostics_every=0)
    try:
        fixed_raw = next(iter(loader))
    except StopIteration as exc:
        raise RuntimeError("programme_only G2.6 received no real train batch") from exc
    fixed_batch = {key: value.to(device, non_blocking=True) if isinstance(value, torch.Tensor) else value
                   for key, value in fixed_raw.items()}
    fixed_batch = _truncate_batch(fixed_batch, overfit_patients)

    def objective() -> tuple[torch.Tensor, dict[str, float]]:
        optimiser.zero_grad(set_to_none=True)
        loss, metrics, _ = trial.step(fixed_batch, epoch=clone_schedule.warmup_epochs)
        if not loss.requires_grad:
            raise RuntimeError("programme_only G2.6 loss is detached from the actual model")
        return loss, metrics

    trajectory: list[dict[str, float]] = []
    initial_loss, initial_metrics = objective()
    initial = float(initial_loss.detach().cpu())
    # Holding `initial_loss` would pin one full activation graph for all 300
    # steps. At the real token budget that graph is gigabytes, and it was
    # enough to OOM a 40 GB A100 before the first optimiser step.
    del initial_loss
    first_gradients: dict[str, float] | None = None
    for step in range(steps):
        loss, _ = objective(); loss.backward()
        gradients = trial._gradient_group_norms()
        if step == 0:
            first_gradients = gradients
            missing = [name for name, value in gradients.items() if value <= 0.0]
            if missing:
                raise RuntimeError(f"programme_only G2.6 has detached trainable groups: {missing}")
        optimiser.step()
        # Record the descent in situ. A single terminal number cannot tell a
        # dead implementation ("flat from step 0") apart from an unconverged
        # one ("still falling at the last step") apart from an unstable one
        # ("fell, then reversed"), and those three demand different responses.
        if step % 25 == 0 or step == steps - 1:
            trajectory.append({"step": int(step), "loss": float(loss.detach().cpu())})
    final_loss, final_metrics = objective()
    final = float(final_loss.detach().cpu())
    tail = [row["loss"] for row in trajectory[-4:]]
    return {
        "batch_patients": int(len(fixed_batch["indices"])),
        "trajectory": trajectory,
        "still_descending": bool(len(tail) >= 2 and tail[-1] < tail[0]),
        "best_loss": float(min(row["loss"] for row in trajectory)) if trajectory else float("nan"),
        "best_step": int(min(trajectory, key=lambda row: row["loss"])["step"]) if trajectory else -1,
        "initial_loss": initial, "final_loss": final,
        "initial_programme": float(initial_metrics.get("programme", float("nan"))),
        "final_programme": float(final_metrics.get("programme", float("nan"))),
        "gradient_norms_first": first_gradients or {},
        "relative_reduction": float((initial - final) / max(abs(initial), 1e-12)),
        "objective_scope": "actual_v2_encoder_and_programme_path_without_decorrelation_floor",
    }


def _overfit_programme_free_contrastive(model: TumorStateV2, schedule: V2LossSchedule,
                                        loader: UncappedHoptimusBatches, device: str,
                                        steps: int = 800, minimum_memory_keys: int = 16,
                                        overfit_patients: int = 16,
                                        overfit_memory_capacity: int = 64,
                                        freeze_memory: bool = True) -> dict[str, object]:
    """G2.6 for D1 on a clone of the *actual* model and trainer path.

    This deliberately does not train post-hoc linear probes on frozen states.
    It exercises the WSI encoder, RNA encoder, shared slots, biology head and
    queue exactly as D1 does, while keeping the experiment model untouched.
    The rank-spreading regularizer is excluded only from this memorisation
    check because it has an irreducible batch-statistics floor; the two new
    D1 terms are the quantities that must be driven near zero.
    """
    # lr=1e-3, not 1e-2. Measured on the real runner at 1e-2 the memorisation
    # DESCENDS then DIVERGES inside the window the gate reads:
    #   step 0 2.4154 -> 150 1.2822 -> 250 0.9432 (best) -> 299 2.1575
    # and by ~2000 steps it reaches nan outright. The heteroscedastic Gaussian
    # NLL's variance head collapses under AdamW at that rate, so the gate was
    # reading a post-divergence value and calling a healthy model dead. Both D1
    # arms use the same setting: their liveness checks are only comparable
    # evidence if they are run identically.
    clone = copy.deepcopy(model).to(device)
    # decorrelation STAYS ON here, unlike the Hallmark arm. Measured at
    # initialisation on a real 16-patient batch: the WSI biology states are
    # 0.736 mutually collinear (unit-norm, std 0.031 across patients) against
    # 0.274 for RNA. InfoNCE cannot separate patients whose WSI vectors nearly
    # coincide, and decorrelation is the only term opposing that collapse --
    # zeroing it disabled the anti-collapse force and then graded the collapse,
    # pinning the contrastive term at chance (ln(80)=4.38; measured 4.27).
    # Its batch-statistics floor is handled by grading the two D1 terms
    # directly rather than the total, below.
    clone_schedule = schedule
    optimiser = torch.optim.AdamW(clone.parameters(), lr=1e-3, weight_decay=0.0)
    trial = V2Trainer(clone, optimiser, clone_schedule, device, gradient_diagnostics_every=0)
    # The queue is sized to the CHECK, not to training. G2.6 asks "can this model
    # memorise one small batch", and its criterion is contrastive <= 0.10. Against
    # the training queue's 4096 DETACHED keys -- encoded by the pre-optimisation
    # model and never refreshed -- that criterion is unreachable by construction:
    # measured, the term sat at 5.62 with a 0.070 reduction over 800 steps while
    # full_consistency reached 0.00023. The queue's gradient path is still
    # exercised (priming below, and the step-0 detached-group check), so nothing
    # that G2.6 exists to catch is lost by matching the queue to the batch scale.
    trial.biology_memory = PairedBiologyMemoryBank(capacity=overfit_memory_capacity)
    iterator = iter(loader)
    priming_batches: list[dict[str, torch.Tensor]] = []
    observed_ids: set[int] = set()
    while len(observed_ids) < minimum_memory_keys:
        try:
            raw = next(iterator)
        except StopIteration as exc:
            raise RuntimeError(
                f"programme_free G2.6 received fewer than {minimum_memory_keys} distinct real train patients"
            ) from exc
        priming_batches.append(raw)
        observed_ids.update(int(value) for value in raw["indices"].tolist())
    trial.prime_biology_memory(priming_batches, minimum_unique_keys=minimum_memory_keys)
    # FREEZE the queue for the memorisation loop. Training refreshes it every
    # step, which is correct there because every step sees new patients. Here one
    # fixed batch is replayed, so a live queue overwrites all `capacity` slots
    # with re-encoded copies of those same patients within `capacity/batch` steps:
    # the negatives then move with the queries and the contrastive term cannot
    # descend. Measured pinned at chance (4.27 vs ln(80)=4.38) across four
    # unrelated interventions. Frozen, the keys stay real held-out train patients
    # from the priming batches, and the queue's gradient path is still exercised.
    # `freeze_memory=False` exists only so the two arms can be measured against
    # each other on real data; the gate itself always runs frozen.
    trial.freeze_biology_memory = freeze_memory
    try:
        fixed_batch = next(iterator)
    except StopIteration:
        # A tiny cohort may fit entirely in the queue.  Reusing one real batch
        # is valid: its matching key is masked and all other queue keys remain
        # real train-patient negatives.
        fixed_batch = priming_batches[-1]
    fixed_batch = {key: value.to(device, non_blocking=True) if isinstance(value, torch.Tensor) else value
                   for key, value in fixed_batch.items()}
    # Same truncation as the Hallmark arm: the two D1 liveness checks must be
    # run at the same batch size or they are not comparable evidence. 16 still
    # clears the InfoNCE min_negatives=8 floor, and the memory queue is primed
    # from the full untruncated batches above.
    fixed_batch = _truncate_batch(fixed_batch, overfit_patients)

    def one_step(*, update: bool) -> tuple[torch.Tensor, dict[str, float]]:
        optimiser.zero_grad(set_to_none=True)
        loss, metrics, _ = trial.step(fixed_batch, epoch=clone_schedule.warmup_epochs)
        if not loss.requires_grad:
            raise RuntimeError("programme_free G2.6 loss is detached from the actual model")
        if update:
            loss.backward()
            gradients = trial._gradient_group_norms()
            missing = [name for name, value in gradients.items() if value <= 0.0]
            if missing:
                raise RuntimeError(f"programme_free G2.6 has detached trainable groups: {missing}")
            optimiser.step()
        return loss.detach(), metrics

    initial, initial_metrics = one_step(update=False)
    first_gradients: dict[str, float] | None = None
    for step in range(steps):
        optimiser.zero_grad(set_to_none=True)
        loss, metrics, _ = trial.step(fixed_batch, epoch=clone_schedule.warmup_epochs)
        loss.backward()
        gradients = trial._gradient_group_norms()
        if step == 0:
            first_gradients = gradients
            missing = [name for name, value in gradients.items() if value <= 0.0]
            if missing:
                raise RuntimeError(f"programme_free G2.6 has detached trainable groups: {missing}")
        optimiser.step()
    final, final_metrics = one_step(update=False)
    initial_value, final_value = float(initial.cpu()), float(final.cpu())
    return {
        "batch_patients": int(len(fixed_batch["indices"])),
        "memory_unique_keys": int(trial.biology_memory.unique_count if trial.biology_memory else 0),
        "memory_frozen": bool(trial.freeze_biology_memory),
        "initial_loss": initial_value, "final_loss": final_value,
        "initial_biology_contrastive": float(initial_metrics.get("biology_contrastive", float("nan"))),
        "final_biology_contrastive": float(final_metrics.get("biology_contrastive", float("nan"))),
        "biology_contrastive_reduction": float(
            (initial_metrics.get("biology_contrastive", float("nan")) - final_metrics.get("biology_contrastive", float("nan")))
            / max(abs(initial_metrics.get("biology_contrastive", 1.0)), 1e-12)),
        "initial_full_consistency": float(initial_metrics.get("biology_full_consistency", float("nan"))),
        "final_full_consistency": float(final_metrics.get("biology_full_consistency", float("nan"))),
        "gradient_norms_first": first_gradients or {},
        "relative_reduction": float((initial_value - final_value) / max(abs(initial_value), 1e-12)),
        "objective_scope": "actual_v2_encoder_and_biology_path_without_decorrelation_floor",
    }


def _require_programme_free_overfit(result: dict[str, object]) -> None:
    """Fail closed on D1 G2.6 before the expensive 40-epoch job starts."""
    initial, final = float(result["initial_loss"]), float(result["final_loss"])
    reduction = float(result["relative_reduction"])
    contrastive = float(result["final_biology_contrastive"])
    consistency = float(result["final_full_consistency"])
    if not all(np.isfinite(value) for value in (initial, final, reduction, contrastive, consistency)):
        raise RuntimeError(f"G2.6 programme_free overfit produced non-finite values: {result}")
    # The G2.6 target is not merely a decreasing curve: the two new D1 terms
    # must be practically memorised. This is intentionally strict enough that
    # a detached or tiny loss cannot be carried into a full GPU run.
    # `final` is the TOTAL, and with decorrelation active it carries an
    # irreducible batch-statistics floor unrelated to memorisation. Grade the
    # two D1 terms themselves, and require the contrastive term -- the one
    # under test -- to have actually descended.
    contrastive_reduction = float(result.get("biology_contrastive_reduction", float("nan")))
    if contrastive > 0.10 or consistency > 0.02 or not np.isfinite(contrastive_reduction) or contrastive_reduction < 0.80:
        raise RuntimeError(
            "G2.6 programme_free overfit failed before training: "
            f"loss={final:.5f}, contrastive={contrastive:.5f} (reduction {contrastive_reduction:.3f}), "
            f"full_consistency={consistency:.5f}, "
            f"reduction={reduction:.3f}; expected near-zero terms and >=80% reduction"
        )


def _require_programme_only_overfit(result: dict[str, object]) -> None:
    """Fail closed on the matched Hallmark-arm G2.6 check."""
    initial, final, reduction = (float(result[name]) for name in ("initial_loss", "final_loss", "relative_reduction"))
    programme = float(result["final_programme"])
    if not all(np.isfinite(value) for value in (initial, final, reduction, programme)):
        raise RuntimeError(f"G2.6 programme_only overfit produced non-finite values: {result}")
    # Gaussian NLL can legitimately become negative as fitted variance falls;
    # hence a small or negative terminal value plus a large reduction is the
    # appropriate analogue of the contrastive arm's near-zero check.
    if final > 0.10 or reduction < 0.80:
        raise RuntimeError(
            "G2.6 programme_only overfit failed before training: "
            f"loss={final:.5f}, programme={programme:.5f}, reduction={reduction:.3f}; "
            f"best={result.get('best_loss')} at step {result.get('best_step')}, "
            f"still_descending={result.get('still_descending')}, "
            f"trajectory={[(row['step'], round(row['loss'], 4)) for row in result.get('trajectory', [])][::3]}; "
            "expected a practically memorised actual-model objective"
        )


def _require_d1_liveness(liveness: dict[str, object], profile: str) -> None:
    """Enforce G2.1--G2.5 for D1 before returning an artifact checkpoint."""
    deltas = {str(name): float(value) for name, value in dict(liveness["parameter_relative_delta"]).items()}
    stalled = [name for name, value in deltas.items() if not np.isfinite(value) or value <= 1e-2]
    if stalled:
        raise RuntimeError(f"G2.1 failed for {profile}: parameter groups did not move >1e-2: {stalled}")
    for phase in ("gradient_norms_first", "gradient_norms_last"):
        gradients = {str(name): float(value) for name, value in dict(liveness[phase]).items()}
        missing = [name for name, value in gradients.items() if not np.isfinite(value) or value <= 0.0]
        if missing:
            raise RuntimeError(f"G2.3 failed for {profile} ({phase}): detached groups {missing}")
    initial, final, reduction = (float(liveness[name]) for name in ("loss_initial", "loss_final", "loss_relative_reduction"))
    if not all(np.isfinite(value) for value in (initial, final, reduction)) or reduction < 0.20:
        raise RuntimeError(f"G2.4 failed for {profile}: final loss did not decrease >=20% ({initial} -> {final})")
    active = {str(name): float(value) for name, value in dict(liveness["active_terms_final"]).items()}
    declared = [str(name) for name in liveness["declared_active_terms"]]
    absent = [name for name in declared if name not in active]
    largest = float(liveness["largest_active_term"])
    # DEAD vs SMALL. G2.2 exists to catch a term that never fired, not one that
    # succeeded: biology_full_consistency is a consistency loss driven toward 0 on
    # purpose, and the overfit gate demands it reach <=0.02, so at epoch 40 it can sit
    # below 1e-4 x largest precisely BECAUSE it worked. Aborting there destroys a
    # healthy arm at the finish line and, under phase_d, the whole sweep. Note also
    # that active_terms mixes weighted and unweighted values, so the ratio is not a
    # contribution comparison. Raise only on genuinely dead terms; record the rest.
    dead = [name for name, value in active.items() if not np.isfinite(value) or value == 0.0]
    small = [name for name, value in active.items() if name not in dead and abs(value) < 1e-4 * largest]
    liveness["g2_2_small_terms"] = small
    if absent or not np.isfinite(largest) or largest <= 0.0 or dead:
        raise RuntimeError(f"G2.2 failed for {profile}: absent={absent}, dead={dead}, largest={largest}")
    slope = float(liveness["tail_loss_slope"])
    # A terminal slope more negative than 1% of the initial loss over the
    # five-epoch tail means the run is visibly still improving; it must be
    # extended before a comparison is interpreted (G2.5).
    # RECORDED, not raised. The bar is scaled by |initial_loss|, which differs by
    # orders of magnitude between a Gaussian-NLL arm (unbounded below, and expected to
    # still be falling at epoch 40 as the fitted variance shrinks) and an InfoNCE arm
    # (~ln|queue|) -- so the "same" threshold is a different threshold per arm and would
    # preferentially destroy one side of the comparison. G2.5 is an interpretation
    # caveat: it must reach the log and the write-up, not kill the artifact.
    liveness["g2_5_still_falling"] = bool(np.isfinite(slope) and slope < -0.01 * max(abs(initial), 1e-6))
    liveness["g2_5_tail_slope"] = slope


def _trained_states_for_profile(objective_profile: str, semantic_dim: int = 0) -> list[str]:
    """Declare only representations with an explicit gradient path.

    The artifact contract is deliberately stricter than tensor availability:
    `TumorStateV2` can compute every view, but diagnostics must never export a
    state that its selected objective did not train.
    """
    if objective_profile == "identity_only":
        return ["wsi_identity", "rna_identity", "full_identity"]
    if objective_profile in {"programme_only", "programme_free"}:
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
    if args.loss_warmup_epochs < 0:
        raise ValueError("--loss-warmup-epochs must be non-negative")
    if any(value < 0.0 for value in (args.programme_warmup_weight, args.programme_weight,
                                     args.programme_neighbourhood_weight, args.programme_supcon_weight,
                                     args.separation_weight, args.variance_weight,
                                     args.decorrelation_weight)):
        raise ValueError("loss weights must be non-negative")
    if args.pretrain_epochs < 0:
        raise ValueError("--pretrain-epochs must be non-negative")
    if args.pretrain_epochs and args.pretrain_checkpoint:
        raise ValueError("choose either --pretrain-epochs or --pretrain-checkpoint, not both")
    if args.objective_profile in {"programme_only", "programme_free"} and not args.fixed_final_epoch:
        raise ValueError(
            "D1 diagnostic profiles require --fixed-final-epoch; validation Hallmark selection is prohibited"
        )
    d2_pair_manifest = _validate_d2_pair(args)
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    data = load_bio_query_data(args.data_config, args.split_file, wsi_mode="hoptimus_patch")
    restriction: dict[str, object] = {"enabled": False}
    if args.restrict_to_split:
        data, restriction = restrict_cohort_to_split(data, args.split_file)
    split_manifest = validate_runtime_split(args.split_file, data.patient_ids, data.cancers, data.split,
                                            args.expected_development_cancers, args.expected_heldout_cancers)
    split_manifest["cohort_restriction"] = restriction
    fit_mask = np.asarray(data.split).astype(str) != "test" if args.fit_development else np.asarray(data.split).astype(str) == "train"
    clinical_manifest = (_standardize_clinical(data, fit_mask) if args.include_clinical
                         else {"enabled": False, "reason": "strict_core_wsi_rna"})
    programme_target_manifest: dict[str, object]
    if args.programme_targets:
        programme_target_manifest = attach_external_programme_targets(data, args.programme_targets, fit_mask)
        means = dict(programme_target_manifest["fit_residual_means"])
    else:
        means = attach_v2_targets(data, fit_mask)
        programme_target_manifest = {"source": "hallmark", "target_names": list(data.hallmark_names),
                                     "target_dimension": int(data.hallmark.shape[1]),
                                     "coverage": int(np.asarray(data.hallmark_present, dtype=bool).sum())}
    if args.fit_programme_legibility:
        programme_target_manifest["legibility_operator"] = fit_programme_legibility_operator(data, fit_mask)
    if args.programme_head_dim < data._v2_programmes.shape[1]:
        raise ValueError("--programme-head-dim cannot be smaller than target width")
    data._v2_programme_head_dim = int(args.programme_head_dim)
    optional_manifest = {"clinical": clinical_manifest}
    for name, path in (("snv", args.snv_features), ("cnv", args.cnv_features)):
        if path:
            optional_manifest[name] = _attach_numeric_table(data, path, name)
    anchor_path = args.mlp_clip_anchor or args.mlp_clip_teacher
    use_anchor = bool(anchor_path) and args.objective_profile not in {"programme_only", "programme_free"}
    if anchor_path and not use_anchor:
        warnings.warn(f"{args.objective_profile} ignores the MLP-CLIP anchor by design", RuntimeWarning)
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
                         programme_dim=int(data._v2_programme_head_dim)).to(device)
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
        model, optimizer,
        V2LossSchedule(
            objective_profile=args.objective_profile,
            warmup_epochs=args.loss_warmup_epochs,
            programme_warmup=args.programme_warmup_weight,
            programme_after_warmup=args.programme_weight,
            neighbourhood_after_warmup=args.programme_neighbourhood_weight,
            supcon_after_warmup=args.programme_supcon_weight,
            separation_after_warmup=args.separation_weight,
            variance_after_warmup=args.variance_weight,
            decorrelation_after_warmup=args.decorrelation_weight,
        ),
        device,
        gradient_diagnostics_every=args.gradient_diagnostics_every,
    )
    biology_memory_keys = 0
    if args.objective_profile == "programme_free":
        # Prime on training patients only.  This makes the new contrastive
        # loss live even when the first real ragged WSI batch has B=1.
        biology_memory_keys = trainer.prime_biology_memory(
            UncappedHoptimusBatches(data, train, args.token_budget, args.seed, shuffle=False,
                                    include_clinical=args.include_clinical),
            minimum_unique_keys=16,
        )
    # Run the decisive G2.6 liveness gate before any long optimization.  It
    # trains only a deep copy, so no model/optimizer state from the experiment
    # itself is mutated.
    d1_overfit: dict[str, object] | None = None
    # The G2.6 batch is SHUFFLED (deterministically, by seed). Unshuffled, the
    # loader returns patients in identifier order, so the first 16 are one
    # cancer -- and after cancer-residualisation their targets are within-cancer
    # deviations while their slides are near-identical. Measured: that batch
    # plateaus at a 0.419 reduction by step 800 and then REVERSES to 0.282 by
    # 1200, while a shuffled batch of the same size reaches 1.217 by step 300.
    # "Can the encoder separate 16 near-identical same-cancer slides by their
    # residual noise" is not the question a liveness gate asks.
    if args.objective_profile == "programme_free":
        d1_overfit = _overfit_programme_free_contrastive(
            model, trainer.schedule,
            UncappedHoptimusBatches(data, train, args.token_budget, args.seed, shuffle=True,
                                    include_clinical=args.include_clinical),
            device,
        )
        _require_programme_free_overfit(d1_overfit)
    elif args.objective_profile == "programme_only":
        d1_overfit = _overfit_programme_only_actual(
            model, trainer.schedule,
            UncappedHoptimusBatches(data, train, args.token_budget, args.seed, shuffle=True,
                                    include_clinical=args.include_clinical),
            device,
        )
        _require_programme_only_overfit(d1_overfit)
    initial_group_parameters = _parameter_snapshot(trainer)
    trained_states = _trained_states_for_profile(args.objective_profile, semantic_dim)
    run_configuration = {
        "seed": args.seed, "epochs": args.epochs, "token_budget": args.token_budget,
        "hidden_dim": args.hidden_dim, "layers": args.layers, "heads": args.heads,
        "learning_rate": args.learning_rate, "weight_decay": args.weight_decay,
        "fixed_final_epoch": bool(args.fixed_final_epoch),
        # This is the sole intervention in the E1 matched-arm experiment.  It
        # must be part of the immutable artifact manifest; otherwise a pair of
        # exports can look provenance-identical even when their optimisation
        # differed, making the rank comparison unauditable.
        "decorrelation_weight": args.decorrelation_weight,
        "loss_warmup_epochs": args.loss_warmup_epochs,
        "programme_warmup_weight": args.programme_warmup_weight,
        "programme_weight": args.programme_weight,
        "programme_neighbourhood_weight": args.programme_neighbourhood_weight,
        "programme_supcon_weight": args.programme_supcon_weight,
        "separation_weight": args.separation_weight,
        "variance_weight": args.variance_weight,
        "gradient_diagnostics_every": args.gradient_diagnostics_every,
        "programme_head_dim": args.programme_head_dim,
        "objective_profile": args.objective_profile,
        "programme_targets": args.programme_targets,
        "biology_contrastive_memory_keys_before_train": biology_memory_keys,
        "pretrain_epochs": args.pretrain_epochs, "pretrain_checkpoint": args.pretrain_checkpoint,
        "anchor_path": anchor_path if use_anchor else "", "anchor_requested_but_disabled": bool(anchor_path) and not use_anchor,
        "fit_population": "development_train_val" if args.fit_development else "train_only_inner_fit",
        "include_clinical": args.include_clinical,
        "snv_features": args.snv_features, "cnv_features": args.cnv_features,
    }
    manifest = {"artifact_version": 4, "protocol": split_manifest["protocol"], "patch_cap": None, "patient_ids": list(data.patient_ids), "cancers": list(data.cancers), "split": list(data.split), "residual_means_fit_on": "train_only", "residual_means": means,
                "trained_states": list(trainable_state_names(*trained_states)),
                "selection_metric": "fixed_final_epoch" if args.fixed_final_epoch else "validation_pareto",
                "coordinates": "metadata_validated_or_disabled", "semantic_dim": semantic_dim,
                "modal_dims": {"clinical": int(data.clinical.shape[1]) if args.include_clinical else 0, "snv": 0 if not hasattr(data, "_v2_snv") else int(data._v2_snv.shape[1]), "cnv": 0 if not hasattr(data, "_v2_cnv") else int(data._v2_cnv.shape[1])},
                "optional_modalities": optional_manifest, "programme_targets": programme_target_manifest,
                "d2_paired_comparison": d2_pair_manifest,
                "mlp_clip_anchor": use_anchor,
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
    epoch_rows: list[dict[str, object]] = []
    for epoch in range(start_epoch, args.epochs):
        loader = UncappedHoptimusBatches(data, train, args.token_budget, args.seed + epoch, include_clinical=args.include_clinical)
        metrics = trainer.train_epoch(loader, epoch)
        validation_metrics = {}; selection_metrics = {}
        if len(validation):
            validation_loader = UncappedHoptimusBatches(data, validation, args.token_budget, args.seed, shuffle=False, include_clinical=args.include_clinical)
            validation_metrics = trainer.evaluate_epoch(validation_loader, epoch)
            if not args.fixed_final_epoch:
                selection_loader = UncappedHoptimusBatches(data, validation, args.token_budget, args.seed, shuffle=False, include_clinical=args.include_clinical)
                selection_metrics = _validation_selection(model, selection_loader, device, np.asarray(data.cancers))
        scheduler.step()
        row = {"epoch": epoch, "learning_rate": scheduler.get_last_lr()[0], "anchor_validation_r10": anchor_validation_r10, **{f"train_{k}": v for k, v in metrics.items()}, **{f"validation_{k}": v for k, v in validation_metrics.items()}, **{f"selection_{k}": v for k, v in selection_metrics.items()}}
        (output / "train_metrics.jsonl").open("a", encoding="utf-8").write(json.dumps(row) + "\n")
        epoch_rows.append(row)
        trainer.save_checkpoint(output / "last.pt", epoch, loader.state_dict(), manifest, scheduler)
        r10, programme = float(selection_metrics.get("retrieval_r10", float("nan"))), float(selection_metrics.get("programme_mean_pearson", float("nan")))
        if not args.fixed_final_epoch and np.isfinite(r10) and r10 > best_retrieval:
            best_retrieval = r10; best_retrieval_epoch = epoch; trainer.save_checkpoint(output / "best_retrieval.pt", epoch, loader.state_dict(), manifest, scheduler)
        if not args.fixed_final_epoch and np.isfinite(programme) and programme > best_programme:
            best_programme = programme; best_programme_epoch = epoch; trainer.save_checkpoint(output / "best_programme.pt", epoch, loader.state_dict(), manifest, scheduler)
        eligible = np.isfinite(r10) and (not np.isfinite(anchor_validation_r10) or r10 >= 0.95 * anchor_validation_r10)
        if not args.fixed_final_epoch and eligible and np.isfinite(programme) and programme > best_pareto:
            best_pareto = programme; best_pareto_epoch = epoch; trainer.save_checkpoint(output / "best_pareto.pt", epoch, loader.state_dict(), manifest, scheduler)
    parameter_deltas = _parameter_deltas(trainer, initial_group_parameters)
    liveness_loader = UncappedHoptimusBatches(data, train, args.token_budget, args.seed, shuffle=False,
                                                include_clinical=args.include_clinical)
    overfit = d1_overfit if args.objective_profile in {"programme_only", "programme_free"} else {"status": "not_applicable_identity_only"}
    last = epoch_rows[-1] if epoch_rows else {}
    first = epoch_rows[0] if epoch_rows else {}
    # G2.4 must compare LIKE WITH LIKE. `train_loss` at epoch 0 is a different
    # objective from `train_loss` at epoch 39: neighbourhood and supcon are off
    # during the warmup and switch on at --loss-warmup-epochs, so the total
    # necessarily JUMPS there (measured: 0.528 at epoch 3 -> 5.330 at epoch 4).
    # Comparing across that boundary made a run that fell 5.33 -> 3.19, a 40%
    # reduction in the objective it actually optimised, read as "did not
    # decrease >=20%" and destroyed it after all 40 epochs had been paid for.
    # Epoch 0 is still recorded; the GATE reads the first comparable epoch.
    comparable = [row for row in epoch_rows if int(row.get("epoch", 0)) >= int(args.loss_warmup_epochs)]
    first_comparable = comparable[0] if comparable else first
    loss_initial = float(first.get("train_loss", float("nan")))
    loss_initial_comparable = float(first_comparable.get("train_loss", float("nan")))
    loss_final = float(last.get("train_loss", float("nan")))
    declared_weights = trainer.schedule.weights(max(args.epochs - 1, 0))
    metric_for_loss = {
        "programme": "train_programme", "neighbourhood": "train_wsi_neighbourhood",
        "supcon": "train_wsi_programme_supcon", "separation": "train_separation",
        "variance": "train_variance_floor", "decorrelation": "train_decorrelation",
        "biology_contrastive": "train_biology_contrastive",
        "biology_full_consistency": "train_biology_full_consistency",
    }
    active_terms = {name: float(last[metric]) for name, metric in metric_for_loss.items()
                    if declared_weights.get(name, 0.0) > 0.0 and metric in last}
    largest_term = max((abs(value) for value in active_terms.values()), default=0.0)
    liveness = {
        "parameter_relative_delta": parameter_deltas,
        "gradient_norms_first": {key.removeprefix("train_gradient_norm_").removesuffix("_first"): value for key, value in first.items() if key.startswith("train_gradient_norm_") and key.endswith("_first")},
        "gradient_norms_last": {key.removeprefix("train_gradient_norm_").removesuffix("_last"): value for key, value in last.items() if key.startswith("train_gradient_norm_") and key.endswith("_last")},
        "loss_initial": loss_initial, "loss_final": loss_final,
        "loss_initial_comparable": loss_initial_comparable,
        "loss_comparable_from_epoch": int(first_comparable.get("epoch", 0)),
        "loss_warmup_epochs": int(args.loss_warmup_epochs),
        "loss_relative_reduction_from_epoch_zero": float(
            (loss_initial - loss_final) / max(abs(loss_initial), 1e-12)),
        "loss_relative_reduction": float(
            (loss_initial_comparable - loss_final) / max(abs(loss_initial_comparable), 1e-12)),
        "tail_loss_slope": float(np.polyfit(np.arange(min(5, len(epoch_rows))),
                                              [float(row.get("train_loss", float("nan"))) for row in epoch_rows[-min(5, len(epoch_rows)):]], 1)[0])
        if len(epoch_rows) >= 2 else float("nan"),
        "declared_active_terms": [name for name, value in declared_weights.items() if value > 0.0],
        "active_terms_final": active_terms, "largest_active_term": largest_term,
        "overfit_one_batch": overfit,
    }
    (output / "liveness.json").write_text(json.dumps(liveness, indent=2), encoding="utf-8")
    if args.objective_profile in {"programme_only", "programme_free"}:
        _require_d1_liveness(liveness, args.objective_profile)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (output / "selection.json").write_text(json.dumps({"anchor_validation_r10": anchor_validation_r10, "best_retrieval_r10": best_retrieval, "best_programme_pearson": best_programme, "best_pareto_programme_pearson": best_pareto, "best_retrieval_epoch": best_retrieval_epoch, "best_programme_epoch": best_programme_epoch, "best_pareto_epoch": best_pareto_epoch, "selection": "fixed_final_epoch" if args.fixed_final_epoch else "validation_selected", "fit_population": "development_train_val" if args.fit_development else "train_only_inner_fit"}, indent=2), encoding="utf-8")
    return _selected_checkpoint(output, slide_pretraining, args.epochs, allow_last=args.fit_development or args.fixed_final_epoch)


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
    parser.add_argument("--objective-profile", choices=("full", "identity_only", "programme_only", "programme_free"), default="full",
                        help="diagnostic objective ablation; preserves data/model while selectively zeroing losses")
    parser.add_argument("--decorrelation-weight", type=float, default=0.04,
                        help="F-R2 biology feature-decorrelation weight; set 0.0 for the collapse-baseline arm of the rank ablation")
    # Explicit knobs for controlled mechanism experiments.  The defaults
    # reproduce the canonical V2 schedule; E1 can instead make the compared
    # objective stationary so its G2 loss-liveness gate is meaningful.
    parser.add_argument("--loss-warmup-epochs", type=int, default=4)
    parser.add_argument("--programme-warmup-weight", type=float, default=0.50)
    parser.add_argument("--programme-weight", type=float, default=1.0)
    parser.add_argument("--programme-neighbourhood-weight", type=float, default=0.20)
    parser.add_argument("--programme-supcon-weight", type=float, default=0.20)
    parser.add_argument("--programme-head-dim", type=int, default=256,
                        help="fixed-width masked programme output adapter; D2 uses 256 for both H and PBS arms")
    parser.add_argument("--programme-targets", default="",
                        help="immutable patient-ID keyed target NPZ (e.g. D2 PBS codes); defaults to Hallmarks")
    parser.add_argument("--fit-programme-legibility", action="store_true",
                        help="fit the same grouped-CV target-axis legibility operator on development patients in either H or PBS arm")
    parser.add_argument("--d2-pair-manifest", default="",
                        help="immutable paired H-vs-PBS manifest; required with --d2-arm to enforce target-only difference")
    parser.add_argument("--d2-arm", default="", choices=("", "H", "I"),
                        help="paired D2 arm label; only valid with --d2-pair-manifest")
    parser.add_argument("--d2-analysis-role", default="", choices=("", "primary", "sensitivity"))
    parser.add_argument("--d2-pbs-components", type=int, default=0, choices=(0, 64, 128, 256))
    parser.add_argument("--separation-weight", type=float, default=0.01)
    parser.add_argument("--variance-weight", type=float, default=0.01)
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
    parser.add_argument("--restrict-to-split", action="store_true",
                        help="treat the split file as the authoritative cohort: drop loaded patients it "
                             "does not assign, and record every dropped identifier in the manifest")
    parser.add_argument("--fit-development", action="store_true",
                        help="final refit on train+validation only; epoch count must already be selected")
    parser.add_argument("--fixed-final-epoch", action="store_true",
                        help="emit last.pt after a predeclared epoch count; do not select a checkpoint with programme labels")
    parser.add_argument("--expected-development-cancers", type=int, default=11)
    parser.add_argument("--expected-heldout-cancers", type=int, default=22)
    print(run(parser.parse_args()))


if __name__ == "__main__":
    main()
