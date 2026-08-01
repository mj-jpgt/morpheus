"""Strict V2.1 molecular, survival, and paired-comparison evaluator.

This is intentionally separate from the legacy comprehensive evaluator.  It
has one canonical patient order, trains every probe on development cancers
only, and writes both scalar metrics and patient-level predictions so method
comparisons are genuinely paired.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             f1_score, r2_score, roc_auc_score)
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from .contracts import require_state, validate_artifact
from .paired_bootstrap import paired_patient_and_cancer_bootstrap, pearson_metric
from .preflight import sha256_json
from .provenance import write_source_manifest
from .tasks import retrieval_task


ALPHAS = tuple(float(value) for value in np.logspace(-3, 4, 16))
PHENOTYPE_TARGETS = {
    "immune_hot": "cytolytic_activity",
    "hypoxia_high": "hypoxia",
    "emt_high": "emt",
    "proliferation_high": "proliferation",
}


class EvaluationError(ValueError):
    """Raised when a strict evaluation precondition is not satisfied."""


def _file_digest(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_table(path: str | Path) -> pd.DataFrame:
    source = Path(path)
    if source.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(source)
    return pd.read_csv(source, sep="\t" if source.suffix.lower() in {".tsv", ".txt"} else ",")


def _validate_target_manifests(
    targets: pd.DataFrame,
    curated_panel_path: str = "",
    matched_control_path: str = "",
) -> dict[str, str]:
    """Validate frozen biological endpoints and matched controls before probes.

    This prevents an evaluator from silently scoring whichever pathway columns
    happened to be present in a target table.  Both manifests are optional
    only for legacy diagnostics; a strict V2.2 invocation supplies both.
    """
    digests: dict[str, str] = {}
    available = set(targets.columns) - {"patient_id"}
    if curated_panel_path:
        panel = json.loads(Path(curated_panel_path).read_text(encoding="utf-8"))
        eligible = set(map(str, panel.get("eligible_targets", [])))
        if not eligible:
            raise EvaluationError("curated panel has no eligible targets")
        missing = sorted(eligible - available)
        if missing:
            raise EvaluationError(f"target table lacks curated endpoints: {missing[:8]}")
        digests["curated_panel"] = _file_digest(curated_panel_path)
    if matched_control_path:
        controls = json.loads(Path(matched_control_path).read_text(encoding="utf-8"))
        # The control builder records names in either a flat mapping or a
        # list of records across historical run versions; validate both.
        encoded = json.dumps(controls, sort_keys=True)
        declared = set()
        if isinstance(controls, dict):
            for value in controls.values():
                if isinstance(value, list):
                    declared.update(map(str, value))
                elif isinstance(value, dict):
                    declared.update(map(str, value.get("controls", value.get("names", []))))
        declared.update(name for name in available if name.startswith("RANDOM_CONTROL__"))
        if not any(name.startswith("RANDOM_CONTROL__") for name in declared):
            raise EvaluationError("matched-control manifest/target table has no random control endpoints")
        # A biological endpoint is evaluable only when its matched controls
        # are actually present; this is checked per target in the row writer.
        digests["matched_control_manifest"] = sha256(encoded.encode("utf-8")).hexdigest()
    return digests


def _artifact(path: str | Path, patient_ids: np.ndarray) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    validate_artifact(path)
    with np.load(path, allow_pickle=False) as raw:
        source_ids = raw["patient_ids"].astype(str)
        lookup = {value: index for index, value in enumerate(source_ids)}
        missing = set(patient_ids) - set(lookup)
        if missing:
            raise EvaluationError(f"{path} is missing canonical patients: {sorted(missing)[:5]}")
        order = np.asarray([lookup[value] for value in patient_ids], dtype=np.int64)
        states: dict[str, np.ndarray] = {}
        for state in ("wsi_biology", "wsi_identity", "rna_biology", "rna_identity", "full_biology", "full_identity", "full_patient"):
            value = require_state(raw, state)
            if value is not None:
                states[state] = value[order]
        manifest = json.loads(str(raw["manifest_json"].item())) if "manifest_json" in raw.files else {}
    if not states:
        raise EvaluationError(f"{path} has no declared evaluation state")
    return manifest, states


def _assert_same_cohort(path: str | Path, patient_ids: np.ndarray, cancers: np.ndarray, split: np.ndarray) -> None:
    """Fail closed on any patient/cancer/split disagreement with the anchor."""
    with np.load(path, allow_pickle=False) as raw:
        source_ids = raw["patient_ids"].astype(str)
        lookup = {value: index for index, value in enumerate(source_ids)}
        if set(lookup) != set(patient_ids):
            raise EvaluationError(f"{path} patient universe differs from canonical anchor")
        order = np.asarray([lookup[value] for value in patient_ids], dtype=np.int64)
        if not np.array_equal(raw["cancers"].astype(str)[order], cancers) or not np.array_equal(raw["split"].astype(str)[order], split):
            raise EvaluationError(f"{path} cancer or split assignment differs from canonical anchor")


def _feature_for_wsi(states: dict[str, np.ndarray]) -> tuple[str, np.ndarray]:
    for state in ("wsi_biology", "wsi_identity"):
        if state in states:
            return state, states[state]
    raise EvaluationError("artifact lacks a WSI-trained state")


def _correlation(actual: np.ndarray, predicted: np.ndarray, *, rank: bool = False) -> float:
    valid = np.isfinite(actual) & np.isfinite(predicted)
    if valid.sum() < 3:
        return float("nan")
    x, y = actual[valid], predicted[valid]
    if rank:
        x = pd.Series(x).rank(method="average").to_numpy()
        y = pd.Series(y).rank(method="average").to_numpy()
    return pearson_metric(x, y)


def _calibration(actual: np.ndarray, predicted: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(actual) & np.isfinite(predicted)
    if valid.sum() < 3 or np.std(predicted[valid]) < 1e-12:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(predicted[valid], actual[valid], 1)
    return float(slope), float(intercept)


def _select_alpha(features: np.ndarray, target: np.ndarray, development: np.ndarray, cancers: np.ndarray) -> float:
    keep = development & np.isfinite(target)
    groups = cancers[keep]
    unique = np.unique(groups)
    if keep.sum() < 8 or len(unique) < 2:
        return 1.0
    folds = GroupKFold(n_splits=min(3, len(unique)))
    scores: dict[float, list[float]] = {alpha: [] for alpha in ALPHAS}
    for fit_rows, valid_rows in folds.split(features[keep], target[keep], groups):
        x_train, x_valid = features[keep][fit_rows], features[keep][valid_rows]
        y_train, y_valid = target[keep][fit_rows], target[keep][valid_rows]
        scaler = StandardScaler().fit(x_train)
        for alpha in ALPHAS:
            prediction = Ridge(alpha=alpha).fit(scaler.transform(x_train), y_train).predict(scaler.transform(x_valid))
            scores[alpha].append(_correlation(y_valid, prediction))
    return max(ALPHAS, key=lambda alpha: np.nanmean(scores[alpha]))


def _fit_predict_ridge(features: np.ndarray, target: np.ndarray, development: np.ndarray, test: np.ndarray, cancers: np.ndarray) -> tuple[np.ndarray, float] | None:
    keep = development & np.isfinite(target)
    if keep.sum() < 8 or test.sum() < 3:
        return None
    alpha = _select_alpha(features, target, development, cancers)
    scaler = StandardScaler().fit(features[keep])
    model = Ridge(alpha=alpha).fit(scaler.transform(features[keep]), target[keep])
    return model.predict(scaler.transform(features[test])), alpha


def _macro_cancer_pearson(actual: np.ndarray, predicted: np.ndarray, cancers: np.ndarray) -> float:
    values = [_correlation(actual[cancers == cancer], predicted[cancers == cancer]) for cancer in np.unique(cancers)]
    values = [value for value in values if np.isfinite(value)]
    return float(np.mean(values)) if values else float("nan")


def _finite_r2(actual: np.ndarray, predicted: np.ndarray) -> float:
    valid = np.isfinite(actual) & np.isfinite(predicted)
    if valid.sum() < 3:
        return float("nan")
    return float(r2_score(actual[valid], predicted[valid]))


def _molecular_rows(method: str, state_name: str, features: np.ndarray, targets: pd.DataFrame, metadata: pd.DataFrame, *, evaluation_partition: str) -> tuple[list[dict[str, Any]], pd.DataFrame]:
    development_labels = ("train",) if evaluation_partition == "val" else ("train", "val")
    development = metadata["split"].astype(str).isin(development_labels).to_numpy()
    test = metadata["split"].astype(str).eq(evaluation_partition).to_numpy()
    cancers = metadata["cancer"].astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    predictions: list[pd.DataFrame] = []
    for target_name in [column for column in targets.columns if column != "patient_id"]:
        target = pd.to_numeric(targets[target_name], errors="coerce").to_numpy(float)
        fitted = _fit_predict_ridge(features, target, development, test, cancers)
        base = {"method": method, "representation_state": state_name, "task": "molecular_prompting", "target": target_name}
        if fitted is None:
            rows.append({**base, "metric": "status", "value": np.nan, "note": "unavailable_insufficient_development_or_test_support"})
            continue
        predicted, alpha = fitted
        actual = target[test]
        test_cancers = cancers[test]
        slope, intercept = _calibration(actual, predicted)
        values = {
            "pearson": _correlation(actual, predicted),
            "spearman": _correlation(actual, predicted, rank=True),
            "r2": _finite_r2(actual, predicted),
            "calibration_slope": slope,
            "calibration_intercept": intercept,
            "macro_cancer_pearson": _macro_cancer_pearson(actual, predicted, test_cancers),
            "selected_ridge_alpha": alpha,
        }
        rows.extend({**base, "metric": key, "value": value} for key, value in values.items())
        predictions.append(pd.DataFrame({"patient_id": metadata.loc[test, "patient_id"].to_numpy(), "cancer": test_cancers,
                                         "method": method, "representation_state": state_name, "target": target_name,
                                         "actual": actual, "prediction": predicted}))
    return rows, pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame()


def _normalised_similarity(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    values = values / np.maximum(np.linalg.norm(values, axis=1, keepdims=True), 1e-12)
    return values @ values.T


def _neighbour_metrics(predicted_features: np.ndarray, programme_targets: np.ndarray,
                       cancers: np.ndarray, *, within_cancer: bool, k: int = 10) -> dict[str, float]:
    """Compare WSI neighbours with hidden RNA-programme neighbours on test only."""
    n = len(predicted_features)
    if n < 3:
        return {"recall_at_10": float("nan"), "ndcg_at_10": float("nan"), "n_queries": 0.0}
    pred_similarity = _normalised_similarity(predicted_features)
    truth_similarity = _normalised_similarity(programme_targets)
    recalls: list[float] = []
    ndcgs: list[float] = []
    for row in range(n):
        candidates = np.arange(n) if not within_cancer else np.flatnonzero(cancers == cancers[row])
        candidates = candidates[candidates != row]
        if not len(candidates):
            continue
        count = min(k, len(candidates))
        predicted_rank = candidates[np.argsort(-pred_similarity[row, candidates])[:count]]
        truth_rank = candidates[np.argsort(-truth_similarity[row, candidates])[:count]]
        truth_set = set(truth_rank.tolist())
        hits = np.asarray([int(item in truth_set) for item in predicted_rank], dtype=float)
        recalls.append(float(hits.mean()))
        ideal = sum(1.0 / np.log2(index + 2) for index in range(count))
        ndcgs.append(float(sum(hit / np.log2(index + 2) for index, hit in enumerate(hits)) / ideal))
    return {"recall_at_10": float(np.mean(recalls)) if recalls else float("nan"),
            "ndcg_at_10": float(np.mean(ndcgs)) if ndcgs else float("nan"),
            "n_queries": float(len(recalls))}


def _programme_neighbour_rows(method: str, state_name: str, features: np.ndarray,
                              targets: pd.DataFrame, metadata: pd.DataFrame,
                              *, evaluation_partition: str) -> list[dict[str, Any]]:
    """Evaluate biology-neighbour retrieval without letting RNA enter the model."""
    target_columns = [column for column in targets if column != "patient_id" and not str(column).startswith("RANDOM_CONTROL__")]
    test = metadata["split"].astype(str).eq(evaluation_partition).to_numpy()
    actual = targets.loc[:, target_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)[test]
    x = np.asarray(features, float)[test]
    cancers = metadata.loc[test, "cancer"].astype(str).to_numpy()
    # Do not impute a patient's hidden RNA programme profile.  A row must have
    # all predeclared biological endpoints available to participate.
    valid = np.isfinite(actual).all(axis=1)
    base = {"method": method, "representation_state": state_name, "task": "programme_neighbour_retrieval"}
    if valid.sum() < 3:
        return [{**base, "metric": "status", "value": np.nan,
                 "note": "unavailable_insufficient_complete_hidden_programme_profiles",
                 "n_targets": len(target_columns), "n_patients": int(valid.sum())}]
    rows: list[dict[str, Any]] = []
    for scope, within in (("global", False), ("within_cancer", True)):
        metrics = _neighbour_metrics(x[valid], actual[valid], cancers[valid], within_cancer=within)
        for metric, value in metrics.items():
            rows.append({**base, "metric": f"{scope}_{metric}", "value": value,
                         "n_targets": len(target_columns), "n_patients": int(valid.sum()),
                         "note": "hidden_rna_programme_similarity_test_only"})
    return rows


def _select_logistic_c(features: np.ndarray, labels: np.ndarray, development: np.ndarray,
                       cancers: np.ndarray) -> float | None:
    keep = development & np.isfinite(labels)
    groups = cancers[keep]
    if keep.sum() < 12 or len(np.unique(labels[keep])) < 2 or len(np.unique(groups)) < 2:
        return None
    candidates = tuple(float(value) for value in np.logspace(-3, 3, 9))
    scores: dict[float, list[float]] = {value: [] for value in candidates}
    for fit_rows, valid_rows in GroupKFold(n_splits=min(3, len(np.unique(groups)))).split(features[keep], labels[keep], groups):
        y_fit, y_valid = labels[keep][fit_rows], labels[keep][valid_rows]
        if len(np.unique(y_fit)) < 2 or len(np.unique(y_valid)) < 2:
            continue
        scaler = StandardScaler().fit(features[keep][fit_rows])
        for c in candidates:
            prediction = LogisticRegression(C=c, max_iter=2000, class_weight="balanced").fit(
                scaler.transform(features[keep][fit_rows]), y_fit
            ).predict(scaler.transform(features[keep][valid_rows]))
            scores[c].append(float(balanced_accuracy_score(y_valid, prediction)))
    viable = [(c, np.mean(values)) for c, values in scores.items() if values]
    return max(viable, key=lambda item: (item[1], -item[0]))[0] if viable else None


def _phenotype_rows(method: str, state_name: str, features: np.ndarray, targets: pd.DataFrame,
                    metadata: pd.DataFrame, *, evaluation_partition: str) -> list[dict[str, Any]]:
    """Development-thresholded WSI phenotype probes with explicit availability."""
    development_labels = ("train",) if evaluation_partition == "val" else ("train", "val")
    development = metadata["split"].astype(str).isin(development_labels).to_numpy()
    test = metadata["split"].astype(str).eq(evaluation_partition).to_numpy()
    cancers = metadata["cancer"].astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    for phenotype, target_name in PHENOTYPE_TARGETS.items():
        base = {"method": method, "representation_state": state_name, "task": "molecular_phenotype", "target": phenotype}
        if target_name not in targets:
            rows.append({**base, "metric": "status", "value": np.nan, "note": f"unavailable_target:{target_name}"})
            continue
        target = pd.to_numeric(targets[target_name], errors="coerce").to_numpy(float)
        train_valid = development & np.isfinite(target)
        test_valid = test & np.isfinite(target)
        if train_valid.sum() < 12 or test_valid.sum() < 3:
            rows.append({**base, "metric": "status", "value": np.nan, "note": "unavailable_insufficient_target_coverage"})
            continue
        # The high/low threshold is frozen from development patients only.
        threshold = float(np.nanmedian(target[train_valid]))
        labels = np.where(np.isfinite(target), (target >= threshold).astype(int), -1)
        c = _select_logistic_c(features, labels, development, cancers)
        if c is None or len(np.unique(labels[test_valid])) < 2:
            rows.append({**base, "metric": "status", "value": np.nan, "note": "unavailable_grouped_cv_or_test_class_support"})
            continue
        scaler = StandardScaler().fit(features[train_valid])
        classifier = LogisticRegression(C=c, max_iter=2000, class_weight="balanced").fit(
            scaler.transform(features[train_valid]), labels[train_valid]
        )
        probability = classifier.predict_proba(scaler.transform(features[test_valid]))[:, 1]
        prediction = (probability >= .5).astype(int)
        truth = labels[test_valid]
        values = {
            "auroc": float(roc_auc_score(truth, probability)),
            "auprc": float(average_precision_score(truth, probability)),
            "balanced_accuracy": float(balanced_accuracy_score(truth, prediction)),
            "macro_f1": float(f1_score(truth, prediction, average="macro")),
            "brier": float(np.mean((probability - truth) ** 2)),
            "threshold_train_median": threshold,
            "selected_logistic_c": c,
        }
        rows.extend({**base, "metric": metric, "value": value,
                     "note": "development_only_threshold_and_cancer_grouped_cv"}
                    for metric, value in values.items())
    return rows


def _survival_rows(method: str, state_name: str, features: np.ndarray, outcomes: pd.DataFrame | None, metadata: pd.DataFrame, *, repeats: int) -> tuple[list[dict[str, Any]], pd.DataFrame]:
    if outcomes is None:
        return ([{"method": method, "representation_state": state_name, "task": "survival", "metric": "status", "value": np.nan, "note": "unavailable_no_cdr_outcomes"}], pd.DataFrame())
    from .survival_evaluation import evaluate_coxnet_endpoint
    rows: list[dict[str, Any]] = []
    risks: list[pd.DataFrame] = []
    for endpoint in ("os", "pfi"):
        result = evaluate_coxnet_endpoint(features, metadata, outcomes, endpoint=endpoint, method=method,
                                          representation_state=state_name, bootstrap_repeats=repeats)
        rows.extend(result.rows)
        if not result.risks.empty:
            risks.append(result.risks)
    return rows, pd.concat(risks, ignore_index=True) if risks else pd.DataFrame()


def _paired_survival_rows(
    challenger: pd.DataFrame,
    teacher: pd.DataFrame,
    *,
    method: str,
    teacher_method: str,
    repeats: int,
) -> list[dict[str, Any]]:
    """Return paired C-index deltas for every shared survival endpoint."""
    if challenger.empty or teacher.empty:
        return []
    from .survival_evaluation import paired_cindex_bootstrap

    rows: list[dict[str, Any]] = []
    # Pair the WSI-only risks only.  A full/RNA risk has no same-input teacher
    # comparator and remains a separately reported reference, not an invalid
    # cross-input superiority comparison.
    challenger = challenger.loc[challenger.representation_state.astype(str).str.startswith("wsi")]
    teacher = teacher.loc[teacher.representation_state.astype(str).str.startswith("wsi")]
    for endpoint in sorted(set(challenger.endpoint.astype(str)) & set(teacher.endpoint.astype(str))):
        left = challenger.loc[challenger.endpoint.astype(str).eq(endpoint)]
        right = teacher.loc[teacher.endpoint.astype(str).eq(endpoint)]
        joined = left[["patient_id", "cancer", "time_days", "event", "risk"]].merge(
            right[["patient_id", "time_days", "event", "risk"]],
            on=["patient_id", "time_days", "event"],
            suffixes=("_challenger", "_teacher"),
            validate="one_to_one",
        )
        if len(joined) < 3:
            continue
        result = paired_cindex_bootstrap(
            joined.time_days.to_numpy(float),
            joined.event.to_numpy(bool),
            joined.risk_teacher.to_numpy(float),
            joined.risk_challenger.to_numpy(float),
            joined.cancer.to_numpy(str),
            repeats=repeats,
        )
        for mode, values in result.items():
            for key, value in values.items():
                if key not in {"mode", "repeats"}:
                    rows.append({
                        "method": method,
                        "representation_state": "wsi",
                        "task": "paired_survival_vs_teacher",
                        "target": endpoint,
                        "metric": f"{mode}_{key}",
                        "value": value,
                        "note": f"teacher={teacher_method}",
                    })
    return rows


def _control_comparison_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Summarise each method's biology-versus-matched-control gap explicitly."""
    frame = pd.DataFrame(rows)
    if frame.empty:
        return []
    molecular = frame.loc[(frame.task == "molecular_prompting") & (frame.metric == "pearson")].copy()
    payload: list[dict[str, Any]] = []
    for method, group in molecular.groupby("method", sort=True):
        values = {str(row.target): float(row.value) for _, row in group.dropna(subset=["target", "value"]).iterrows()}
        deltas = []
        for target, value in values.items():
            if target.startswith("RANDOM_CONTROL__"):
                continue
            controls = [control for name, control in values.items() if name.startswith(f"RANDOM_CONTROL__{target}__")]
            if controls:
                deltas.append(value - float(np.mean(controls)))
        if deltas:
            payload.extend([
                {"method": method, "representation_state": "wsi", "task": "random_signature_control", "metric": "mean_matched_biological_minus_control_pearson", "value": float(np.mean(deltas)), "n_targets": len(deltas), "note": "per-target_matched_controls"},
                {"method": method, "representation_state": "wsi", "task": "random_signature_control", "metric": "fraction_biological_beats_matched_controls", "value": float(np.mean(np.asarray(deltas) > 0)), "n_targets": len(deltas), "note": "per-target_matched_controls"},
            ])
        else:
            payload.append({"method": method, "representation_state": "wsi", "task": "random_signature_control", "metric": "status", "value": np.nan, "note": "unavailable_no_matched_control_pairs"})
    return payload


def evaluate(artifacts: list[str], targets_path: str, output: str, *, outcomes_path: str = "", teacher_map_path: str = "", repeats: int = 2000, evaluation_partition: str = "test", curated_panel_path: str = "", matched_control_path: str = "") -> None:
    if evaluation_partition not in {"val", "test"}:
        raise EvaluationError("evaluation_partition must be val or test")
    targets = _read_table(targets_path)
    if "patient_id" not in targets:
        raise EvaluationError("targets need patient_id")
    targets["patient_id"] = targets["patient_id"].astype(str)
    if targets.patient_id.duplicated().any():
        raise EvaluationError("targets contain duplicate patient IDs")
    target_manifest_digests = _validate_target_manifests(targets, curated_panel_path, matched_control_path)
    first = Path(artifacts[0])
    with np.load(first, allow_pickle=False) as raw:
        patient_ids = raw["patient_ids"].astype(str)
        cancers = raw["cancers"].astype(str)
        split = raw["split"].astype(str)
    metadata = pd.DataFrame({"patient_id": patient_ids, "cancer": cancers, "split": split}).merge(targets, on="patient_id", how="left", validate="one_to_one")
    # Target availability is inherently per programme: a canonical WSI
    # participant can lack RNA for every derived target while other patients
    # retain valid labels.  Downstream probes already mask each target and
    # write an explicit unavailable row when support is insufficient.  A
    # row-wise complete-case gate would contradict that contract and can make
    # a valid evaluation fail solely because an optional modality is missing.
    outcomes = _read_table(outcomes_path) if outcomes_path else None
    output_path = Path(output); output_path.mkdir(parents=True, exist_ok=True)
    artifact_digests = {str(Path(path)): _file_digest(path) for path in artifacts}
    input_digests = {"targets": _file_digest(targets_path), **artifact_digests}
    if outcomes_path:
        input_digests["outcomes"] = _file_digest(outcomes_path)
    write_source_manifest(output_path / "source_manifest.json", configuration={"artifacts": artifacts, "targets": targets_path, "outcomes": outcomes_path, "repeats": repeats, "input_digests": input_digests})
    all_rows: list[dict[str, Any]] = []
    all_predictions: list[pd.DataFrame] = []
    all_risks: list[pd.DataFrame] = []
    for artifact_path in artifacts:
        _assert_same_cohort(artifact_path, patient_ids, cancers, split)
        manifest, states = _artifact(artifact_path, patient_ids)
        method = Path(artifact_path).stem
        state_name, wsi = _feature_for_wsi(states)
        if "wsi_identity" in states and "rna_identity" in states:
            test_rows = metadata["split"].astype(str).eq(evaluation_partition).to_numpy()
            retrieval = retrieval_task(states["wsi_identity"][test_rows], states["rna_identity"][test_rows], metadata.loc[test_rows, "cancer"].astype(str).to_numpy())
            all_rows.extend(
                {"method": method, "representation_state": "wsi_identity", "task": "retrieval", "metric": key, "value": value}
                for key, value in {**retrieval.global_metrics, "within_cancer_r10": retrieval.within_cancer_r10}.items()
            )
        target_frame = metadata[["patient_id", *[x for x in targets.columns if x != "patient_id"]]]
        rows, prediction = _molecular_rows(method, state_name, wsi, target_frame, metadata[["patient_id", "cancer", "split"]], evaluation_partition=evaluation_partition)
        all_rows.extend(rows)
        all_rows.extend(_programme_neighbour_rows(method, state_name, wsi, target_frame,
                                                   metadata[["patient_id", "cancer", "split"]],
                                                   evaluation_partition=evaluation_partition))
        all_rows.extend(_phenotype_rows(method, state_name, wsi, target_frame,
                                        metadata[["patient_id", "cancer", "split"]],
                                        evaluation_partition=evaluation_partition))
        if not prediction.empty:
            all_predictions.append(prediction)
        if evaluation_partition == "test":
            survival_rows, risks = _survival_rows(method, state_name, wsi, outcomes, metadata[["patient_id", "cancer", "split"]], repeats=repeats)
            all_rows.extend(survival_rows)
            if not risks.empty:
                all_risks.append(risks)
            # RNA-only and fused states are references for prognosis only.
            # They never enter WSI-only molecular-prompting claims.
            for reference_state in ("rna_identity", "full_patient"):
                if reference_state in states:
                    reference_rows, reference_risks = _survival_rows(
                        method, reference_state, states[reference_state], outcomes,
                        metadata[["patient_id", "cancer", "split"]], repeats=repeats,
                    )
                    all_rows.extend(reference_rows)
                    if not reference_risks.empty:
                        all_risks.append(reference_risks)
        all_rows.append({"method": method, "representation_state": state_name, "task": "artifact", "metric": "declared_state_count", "value": len(states), "note": json.dumps(manifest, sort_keys=True)})
    prediction_table = pd.concat(all_predictions, ignore_index=True) if all_predictions else pd.DataFrame(columns=["patient_id", "cancer", "method", "representation_state", "target", "actual", "prediction"])
    prediction_table.to_parquet(output_path / "molecular_predictions.parquet", index=False)
    risk_table = pd.concat(all_risks, ignore_index=True) if all_risks else pd.DataFrame()
    risk_table.to_parquet(output_path / "survival_risks.parquet", index=False)
    teacher_map = json.loads(Path(teacher_map_path).read_text(encoding="utf-8")) if teacher_map_path else {}
    if teacher_map:
        for method, teacher in teacher_map.items():
            challenger = prediction_table.loc[prediction_table.method.eq(method)]
            reference = prediction_table.loc[prediction_table.method.eq(teacher)]
            for target in sorted(set(challenger.target) & set(reference.target)):
                # Random controls are aggregated in the explicit control
                # report below.  Thousands of redundant per-control
                # bootstraps add hours of CPU time without adding a new
                # inferential unit to the predeclared biological comparison.
                if str(target).startswith("RANDOM_CONTROL__"):
                    continue
                joined = challenger.loc[challenger.target.eq(target), ["patient_id", "cancer", "actual", "prediction"]].merge(
                    reference.loc[reference.target.eq(target), ["patient_id", "actual", "prediction"]], on=["patient_id", "actual"], suffixes=("_challenger", "_teacher"), validate="one_to_one"
                )
                if len(joined) < 3:
                    continue
                result = paired_patient_and_cancer_bootstrap(pearson_metric, joined.actual.to_numpy(), joined.prediction_teacher.to_numpy(), joined.prediction_challenger.to_numpy(), joined.cancer.to_numpy(), repeats=repeats)
                for mode, values in result.items():
                    for key, value in values.items():
                        if key not in {"mode", "repeats"}:
                            all_rows.append({"method": method, "representation_state": "wsi", "task": "paired_vs_teacher", "target": target, "metric": f"{mode}_{key}", "value": value, "note": f"teacher={teacher}"})
            all_rows.extend(_paired_survival_rows(
                risk_table.loc[risk_table.method.eq(method)] if not risk_table.empty else risk_table,
                risk_table.loc[risk_table.method.eq(teacher)] if not risk_table.empty else risk_table,
                method=method,
                teacher_method=teacher,
                repeats=repeats,
            ))
    all_rows.extend(_control_comparison_rows(all_rows))
    pd.DataFrame(all_rows).to_csv(output_path / "task_rows.csv", index=False)
    protocol = {"patient_count": len(metadata), "cohort_digest": sha256_json(metadata[["patient_id", "cancer", "split"]].to_dict("list")), "development_partitions": ["train"] if evaluation_partition == "val" else ["train", "val"], "evaluation_partition": evaluation_partition, "ridge_alphas": list(ALPHAS), "bootstrap_repeats": repeats, "teacher_map": teacher_map, "input_digests": input_digests, "target_manifest_digests": target_manifest_digests}
    (output_path / "evaluation_protocol.json").write_text(json.dumps(protocol, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--outcomes", default="")
    parser.add_argument("--teacher-map", default="")
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    parser.add_argument("--evaluation-partition", choices=("val", "test"), default="test")
    parser.add_argument("--curated-panel-manifest", default="")
    parser.add_argument("--matched-control-manifest", default="")
    args = parser.parse_args()
    evaluate(args.artifacts, args.targets, args.output, outcomes_path=args.outcomes, teacher_map_path=args.teacher_map, repeats=args.bootstrap_repeats, evaluation_partition=args.evaluation_partition, curated_panel_path=args.curated_panel_manifest, matched_control_path=args.matched_control_manifest)


if __name__ == "__main__":
    main()
