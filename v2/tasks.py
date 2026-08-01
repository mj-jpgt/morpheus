"""Shared task metrics for V2 and frozen baseline representation artifacts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from sklearn.calibration import calibration_curve
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, balanced_accuracy_score, brier_score_loss, f1_score, roc_auc_score, top_k_accuracy_score

from morpheus.src.eval.retrieval_metrics import paired_retrieval_metrics


def train_ridge_predict(features: np.ndarray, targets: np.ndarray, train: np.ndarray, evaluate: np.ndarray) -> np.ndarray:
    mean = features[train].mean(0, keepdims=True)
    scale = np.maximum(features[train].std(0, keepdims=True), 1e-6)
    target_mean = targets[train].mean(0, keepdims=True)
    target_scale = np.maximum(targets[train].std(0, keepdims=True), 1e-6)
    model = RidgeCV(alphas=np.logspace(-3, 4, 16)).fit((features[train] - mean) / scale, (targets[train] - target_mean) / target_scale)
    return model.predict((features[evaluate] - mean) / scale) * target_scale + target_mean


def mean_column_correlation(actual: np.ndarray, predicted: np.ndarray) -> float:
    values = [np.corrcoef(actual[:, i], predicted[:, i])[0, 1] for i in range(actual.shape[1]) if np.std(actual[:, i]) > 1e-8 and np.std(predicted[:, i]) > 1e-8]
    return float(np.nanmean(values)) if values else float("nan")


def expected_calibration_error(probability: np.ndarray, target: np.ndarray, bins: int = 15) -> float:
    probability, target = np.asarray(probability), np.asarray(target)
    edges = np.linspace(0.0, 1.0, bins + 1)
    result = 0.0
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (probability >= left) & (probability < right if right < 1 else probability <= right)
        if mask.any():
            result += mask.mean() * abs(probability[mask].mean() - target[mask].mean())
    return float(result)


def binary_metrics(target: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    target = np.asarray(target, dtype=int)
    probability = np.asarray(probability, dtype=float)
    prediction = probability >= 0.5
    if len(np.unique(target)) < 2:
        return {key: float("nan") for key in ("auroc", "auprc", "balanced_accuracy", "macro_f1", "brier", "ece")}
    return {
        "auroc": float(roc_auc_score(target, probability)),
        "auprc": float(average_precision_score(target, probability)),
        "balanced_accuracy": float(balanced_accuracy_score(target, prediction)),
        "macro_f1": float(f1_score(target, prediction, average="macro")),
        "brier": float(brier_score_loss(target, probability)),
        "ece": expected_calibration_error(probability, target),
    }


def multiclass_metrics(target: np.ndarray, scores: np.ndarray, labels: np.ndarray, top_k: tuple[int, ...] = (1, 3, 5)) -> dict[str, float]:
    target = np.asarray(target)
    labels = np.asarray(labels)
    prediction = labels[np.argmax(scores, axis=1)]
    out = {"balanced_accuracy": float(balanced_accuracy_score(target, prediction)), "macro_f1": float(f1_score(target, prediction, average="macro")), "accuracy": float(accuracy_score(target, prediction))}
    encoded = np.searchsorted(labels, target)
    for k in top_k:
        if k <= len(labels):
            out[f"top_{k}_accuracy"] = float(top_k_accuracy_score(encoded, scores, k=k, labels=np.arange(len(labels))))
    return out


def bootstrap_patients(metric: Callable[[np.ndarray], float], count: int, repeats: int = 1000, seed: int = 42) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = [metric(rng.integers(0, count, count)) for _ in range(repeats)]
    return float(np.nanquantile(values, 0.025)), float(np.nanquantile(values, 0.975))


@dataclass(frozen=True)
class RetrievalResult:
    global_metrics: dict[str, float]
    within_cancer_r10: float


def retrieval_task(wsi: np.ndarray, rna: np.ndarray, cancers: np.ndarray) -> RetrievalResult:
    global_metrics = paired_retrieval_metrics(wsi, rna, (1, 5, 10), cancers.tolist(), cancers.tolist())
    values = []
    for cancer in np.unique(cancers):
        indices = np.where(cancers == cancer)[0]
        if len(indices) > 1:
            values.append(paired_retrieval_metrics(wsi[indices], rna[indices], (10,))["recall_at_10"])
    return RetrievalResult(global_metrics, float(np.mean(values)) if values else float("nan"))


def phenotype_prediction(features: np.ndarray, target: np.ndarray, train: np.ndarray, test: np.ndarray) -> dict[str, float]:
    """Train a train-only standardized logistic probe for every representation."""
    mean, scale = features[train].mean(0), np.maximum(features[train].std(0), 1e-6)
    classifier = LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced", random_state=42)
    classifier.fit((features[train] - mean) / scale, target[train])
    return binary_metrics(target[test], classifier.predict_proba((features[test] - mean) / scale)[:, 1])


def harrell_cindex(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    """Dependency-free Harrell concordance for a fixed outer-test prediction."""
    concordant = comparable = 0.0
    for i in range(len(time)):
        for j in range(i + 1, len(time)):
            if time[i] == time[j]:
                continue
            first, second = (i, j) if time[i] < time[j] else (j, i)
            if not event[first]:
                continue
            comparable += 1
            if risk[first] > risk[second]: concordant += 1
            elif risk[first] == risk[second]: concordant += .5
    return float(concordant / comparable) if comparable else float("nan")


def missingness_selective_risk(error: np.ndarray, uncertainty: np.ndarray, fractions: tuple[float, ...] = (.5, .7, .9, 1.0)) -> dict[str, float]:
    """Risk after retaining the least-uncertain patients; no retraining required."""
    order = np.argsort(np.asarray(uncertainty))
    result = {}
    for fraction in fractions:
        count = max(1, int(round(len(order) * fraction)))
        result[f"risk_at_{fraction:.2f}"] = float(np.mean(np.asarray(error)[order[:count]]))
    return result


def cancer_controlled_discordance(wsi_score: np.ndarray, rna_score: np.ndarray, cancers: np.ndarray) -> dict[str, np.ndarray | float]:
    """Return patient discordance after centering both immune scores within cancer."""
    wsi, rna, cancer = map(np.asarray, (wsi_score, rna_score, cancers))
    residual_wsi, residual_rna = np.zeros_like(wsi, dtype=float), np.zeros_like(rna, dtype=float)
    for label in np.unique(cancer):
        idx = cancer == label
        residual_wsi[idx] = wsi[idx] - wsi[idx].mean()
        residual_rna[idx] = rna[idx] - rna[idx].mean()
    discordance = residual_wsi - residual_rna
    return {"discordance": discordance, "mean_absolute_discordance": float(np.mean(np.abs(discordance))), "correlation": float(np.corrcoef(residual_wsi, residual_rna)[0, 1]) if len(wsi) > 1 else float("nan")}
