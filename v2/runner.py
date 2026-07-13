"""Canonical uncapped TCGA V2 training path.

This module deliberately reads the canonical H-Optimus store directly.  It
does not call the legacy capped patch-batch helper.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterator

import numpy as np
import torch

from morpheus.src.training.train_bio_query_former import load_bio_query_data

from .data import DynamicTokenBatchSampler
from .model import TumorStateV2, V2ModelConfig
from .training import V2LossSchedule, V2Trainer


def residualise_programmes(values: np.ndarray, cancers: np.ndarray, train: np.ndarray) -> tuple[np.ndarray, dict[str, list[float]]]:
    """Fit cancer intercepts on training patients only; unseen cancers use a train global mean."""
    global_mean = values[train].mean(0)
    means = {str(c): values[train & (cancers == c)].mean(0) for c in np.unique(cancers[train])}
    fitted = np.vstack([means.get(str(c), global_mean) for c in cancers])
    return (values - fitted).astype(np.float32), {c: v.tolist() for c, v in means.items()}


class UncappedHoptimusBatches:
    """Load every canonical patch for a dynamic patient batch exactly once."""

    def __init__(self, data, indices: np.ndarray, token_budget: int, seed: int) -> None:
        if data.hoptimus_store is None:
            raise ValueError("V2 requires the canonical H-Optimus patch store")
        self.data, self.indices, self.seed = data, np.asarray(indices, dtype=np.int64), int(seed)
        index = data.hoptimus_store.index_path
        rows = __import__("pandas").read_parquet(index).set_index("patient_id")
        counts = [int(rows.loc[data.patient_ids[i], "n_tokens"]) for i in self.indices]
        self.sampler = DynamicTokenBatchSampler(counts, token_budget, seed)

    def __iter__(self) -> Iterator[dict[str, torch.Tensor]]:
        for local_group in self.sampler.batches(self.seed, shuffle=True):
            group = self.indices[local_group]
            bags, slides, coords = [], [], []
            loaded = self.data.hoptimus_store.load_many_patient_tokens(
                [self.data.patient_ids[index] for index in group], max_tokens=None, seed=self.seed
            )
            for index, (features, metadata) in zip(group, loaded):
                if len(features) == 0:
                    raise RuntimeError(f"patient {self.data.patient_ids[index]} has no canonical tokens")
                codes, _ = __import__("pandas").factorize(metadata.slide_id.astype(str), sort=True)
                bags.append(features); slides.append(codes.astype(np.int64)); coords.append(np.zeros((len(features), 2), dtype=np.float32))
            width = max(len(bag) for bag in bags)
            patches = np.zeros((len(group), width, 1536), dtype=np.float32)
            mask = np.zeros((len(group), width), dtype=bool)
            slide_ids = np.zeros((len(group), width), dtype=np.int64)
            coordinate_batch = np.zeros((len(group), width, 2), dtype=np.float32)
            for row, (bag, ids, xy) in enumerate(zip(bags, slides, coords)):
                patches[row, :len(bag)] = bag; mask[row, :len(bag)] = True
                slide_ids[row, :len(bag)] = ids; coordinate_batch[row, :len(bag)] = xy
            batch = {
                "patches": torch.from_numpy(patches), "patch_mask": torch.from_numpy(mask),
                "slide_ids": torch.from_numpy(slide_ids), "coordinates": torch.from_numpy(coordinate_batch),
                "rna": torch.from_numpy(self.data.rna[group].astype(np.float32)),
                "rna_present": torch.ones(len(group), dtype=torch.bool),
                "programme_target": torch.from_numpy(self.data._v2_programmes[group]),
                "programme_present": torch.from_numpy(self.data.hallmark_present[group].astype(bool)),
                "programme_positive_mask": torch.from_numpy(self.data._v2_positive[np.ix_(group, group)]),
            }
            if hasattr(self.data, "_teacher_wsi"):
                batch["teacher_wsi"] = torch.from_numpy(self.data._teacher_wsi[group])
                batch["teacher_rna"] = torch.from_numpy(self.data._teacher_rna[group])
            yield batch

    def state_dict(self) -> dict[str, object]:
        return {"indices": self.indices.tolist(), "token_budget": self.sampler.token_budget, "seed": self.seed, "coverage": self.sampler.coverage(self.seed).tolist()}


def attach_v2_targets(data) -> dict[str, list[float]]:
    train = np.asarray(data.split) == "train"
    programmes, means = residualise_programmes(np.asarray(data.hallmark, dtype=np.float32), np.asarray(data.cancers), train)
    data._v2_programmes = programmes
    normalized = programmes / np.maximum(np.linalg.norm(programmes, axis=1, keepdims=True), 1e-6)
    similarity = normalized @ normalized.T
    data._v2_positive = (similarity >= 0.60) & train[:, None] & train[None, :]
    np.fill_diagonal(data._v2_positive, False)
    return means


def attach_mlp_clip_teacher(data, path: str) -> None:
    teacher = np.load(path, allow_pickle=False)
    required = {"patient_ids", "split", "wsi", "rna"}
    if not required.issubset(teacher.files):
        raise ValueError("invalid MLP-CLIP teacher artifact")
    lookup = {str(pid): row for row, pid in enumerate(teacher["patient_ids"].astype(str))}
    rows = np.asarray([lookup.get(str(pid), -1) for pid in data.patient_ids])
    if (rows < 0).any() or teacher["wsi"].shape[1] != 256 or teacher["rna"].shape[1] != 256:
        raise ValueError("teacher does not exactly cover the paired cohort with 256-D embeddings")
    if not np.array_equal(teacher["split"][rows].astype(str), np.asarray(data.split).astype(str)):
        raise ValueError("teacher split labels differ from active split")
    data._teacher_wsi = teacher["wsi"][rows].astype(np.float32)
    data._teacher_rna = teacher["rna"][rows].astype(np.float32)


def run(args: argparse.Namespace) -> Path:
    torch.manual_seed(args.seed); np.random.seed(args.seed)
    data = load_bio_query_data(args.data_config, args.split_file, wsi_mode="hoptimus_patch")
    means = attach_v2_targets(data)
    if args.mlp_clip_teacher:
        attach_mlp_clip_teacher(data, args.mlp_clip_teacher)
    train = np.where(np.asarray(data.split) == "train")[0]
    device = "cuda" if args.device == "auto" and torch.cuda.is_available() else args.device
    config = V2ModelConfig(rna_dim=int(data.rna.shape[1]), hidden_dim=args.hidden_dim, layers=args.layers, heads=args.heads)
    model = TumorStateV2(config, programme_dim=int(data.hallmark.shape[1])).to(device)
    trainer = V2Trainer(model, torch.optim.AdamW(model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay), V2LossSchedule(teacher_warmup_epochs=args.teacher_warmup_epochs), device)
    output = Path(args.output_dir); output.mkdir(parents=True, exist_ok=True)
    manifest = {"protocol": "heldout_cancer_11v22", "patch_cap": None, "patient_ids": list(data.patient_ids), "cancers": list(data.cancers), "split": list(data.split), "residual_means_fit_on": "train_only", "residual_means": means}
    for epoch in range(args.epochs):
        loader = UncappedHoptimusBatches(data, train, args.token_budget, args.seed + epoch)
        metrics = trainer.train_epoch(loader, epoch)
        (output / "train_metrics.jsonl").open("a", encoding="utf-8").write(json.dumps({"epoch": epoch, **metrics}) + "\n")
        trainer.save_checkpoint(output / "last.pt", epoch, loader.state_dict(), manifest)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return output / "last.pt"


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
    parser.add_argument("--teacher-warmup-epochs", type=int, default=4)
    print(run(parser.parse_args()))


if __name__ == "__main__":
    main()
