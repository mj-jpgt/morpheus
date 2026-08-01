"""Shared, leakage-checked evaluation for V1, V2, and frozen baselines."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import json
import numpy as np
import torch
from .tasks import mean_column_correlation, multiclass_metrics, retrieval_task, train_ridge_predict
from .text_prototypes import TextPrototypeMapper, few_shot_episodes, zero_shot_predict


@dataclass(frozen=True)
class RepresentationArtifact:
    patient_id: np.ndarray
    cancer: np.ndarray
    split: np.ndarray
    states: dict[str, np.ndarray]
    method: str

    def validate(self) -> None:
        n = len(self.patient_id)
        if len(set(self.patient_id.astype(str))) != n:
            raise ValueError("representation artifact contains duplicate patients")
        if any(len(v) != n or v.ndim != 2 for v in self.states.values()):
            raise ValueError("every representation state must be a [patient, feature] matrix")
        if not {"train", "test"}.issubset(set(self.split.astype(str))):
            raise ValueError("artifact must contain explicit train and test membership")


def save_artifact(path: str | Path, artifact: RepresentationArtifact) -> Path:
    artifact.validate(); output = Path(path); output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, patient_id=artifact.patient_id, cancer=artifact.cancer, split=artifact.split, method=artifact.method, **{"state__" + k: v for k, v in artifact.states.items()})
    return output


def load_artifact(path: str | Path) -> RepresentationArtifact:
    raw = np.load(path, allow_pickle=False)
    artifact = RepresentationArtifact(raw["patient_id"].astype(str), raw["cancer"].astype(str), raw["split"].astype(str), {k[7:]: raw[k].astype(np.float32) for k in raw.files if k.startswith("state__")}, str(raw["method"].item()))
    artifact.validate(); return artifact


def evaluate_programmes(state: np.ndarray, programmes: np.ndarray, split: np.ndarray) -> dict[str, float]:
    train, test = split == "train", split == "test"
    predicted = train_ridge_predict(state, programmes, train, test)
    return {"test_mean_pearson": mean_column_correlation(programmes[test], predicted), "n_train": int(train.sum()), "n_test": int(test.sum())}


def _fit_mapper(state: np.ndarray, cancer: np.ndarray, train: np.ndarray, classes: np.ndarray, text: np.ndarray, epochs: int = 300) -> TextPrototypeMapper:
    lookup = {label: i for i, label in enumerate(classes)}
    keep = train & np.isin(cancer, classes)
    target = np.asarray([lookup[x] for x in cancer[keep]], dtype=np.int64)
    mapper = TextPrototypeMapper(state.shape[1], text.shape[1])
    opt = torch.optim.AdamW(mapper.parameters(), lr=1e-3, weight_decay=1e-3)
    x, y, prototype = torch.tensor(state[keep]), torch.tensor(target), torch.tensor(text)
    mapper.train()
    for _ in range(epochs):
        opt.zero_grad(); loss = torch.nn.functional.cross_entropy(mapper.logits(x, prototype), y); loss.backward(); opt.step()
    return mapper.eval()


def evaluate_zero_shot(state: np.ndarray, cancer: np.ndarray, split: np.ndarray, seen: Iterable[str], seen_text_prototypes: np.ndarray, candidates: Iterable[str], candidate_text_prototypes: np.ndarray) -> dict[str, float]:
    candidates = np.asarray(sorted(candidates)); seen = np.asarray(sorted(seen))
    if len(candidates) != len(candidate_text_prototypes) or len(seen) != len(seen_text_prototypes):
        raise ValueError("each text-prototype matrix must follow its corresponding label order")
    mapper = _fit_mapper(state, cancer, split == "train", seen, seen_text_prototypes)
    # Mapping is trained only with seen label descriptions; all candidate descriptions
    # are frozen and may include unseen labels at test time.
    test = (split == "test") & np.isin(cancer, candidates)
    with torch.no_grad():
        projected = mapper.mapper(torch.tensor(state[test])).numpy()
    scores = zero_shot_predict(projected, candidate_text_prototypes)
    return multiclass_metrics(cancer[test], scores, candidates)


def evaluate_direct_zero_shot(image_space_state: np.ndarray, cancer: np.ndarray, split: np.ndarray,
                              candidates: Iterable[str], candidate_text_prototypes: np.ndarray) -> dict[str, float]:
    """Evaluate frozen image--text cosine directly, with no fitted label mapper."""
    candidates = np.asarray(sorted(candidates))
    if len(candidates) != len(candidate_text_prototypes):
        raise ValueError("candidate text prototypes must follow sorted candidate labels")
    test = (split == "test") & np.isin(cancer, candidates)
    if not test.any():
        return {"accuracy": float("nan")}
    return multiclass_metrics(cancer[test], zero_shot_predict(image_space_state[test], candidate_text_prototypes), candidates)


def evaluate_few_shot(state: np.ndarray, cancer: np.ndarray, split: np.ndarray, candidates: Iterable[str], k_values: tuple[int, ...] = (1, 2, 5, 10, 20), episodes: int = 100, seed: int = 42) -> list[dict[str, float]]:
    candidates = np.asarray(sorted(candidates)); test = (split == "test") & np.isin(cancer, candidates)
    rows = few_shot_episodes(state[test], cancer[test], candidates, k_values, episodes, seed)
    out = []
    for k in k_values:
        values = np.asarray([row["accuracy"] for row in rows if row["k"] == k], dtype=float)
        out.append({"k": k, "episodes": int(len(values)), "accuracy_mean": float(np.mean(values)) if len(values) else float("nan"), "accuracy_ci95_low": float(np.quantile(values, .025)) if len(values) else float("nan"), "accuracy_ci95_high": float(np.quantile(values, .975)) if len(values) else float("nan")})
    return out


def evaluate_retrieval_artifact(artifact: RepresentationArtifact, wsi_key: str = "wsi_identity", rna_key: str = "rna_identity") -> dict[str, float]:
    test = artifact.split == "test"
    result = retrieval_task(artifact.states[wsi_key][test], artifact.states[rna_key][test], artifact.cancer[test])
    return {**result.global_metrics, "within_cancer_r10": result.within_cancer_r10}


def write_task_result(path: str | Path, payload: dict) -> Path:
    output = Path(path); output.parent.mkdir(parents=True, exist_ok=True); output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"); return output
