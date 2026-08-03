"""Leakage-safe frozen baseline exports in the V2 NPZ artifact contract."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

from .contracts import ARTIFACT_VERSION, trainable_state_names


class BaselineExportError(ValueError):
    pass


def _normalise(value: np.ndarray) -> np.ndarray:
    value = np.asarray(value, dtype=np.float32)
    return value / np.maximum(np.linalg.norm(value, axis=1, keepdims=True), 1e-8)


def _validate_inputs(patient_ids, cancers, split, wsi, rna: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    ids, cancer, parts = map(lambda value: np.asarray(value).astype(str), (patient_ids, cancers, split))
    x = np.asarray(wsi, dtype=np.float32)
    y = None if rna is None else np.asarray(rna, dtype=np.float32)
    if len(ids) == 0 or len(set(ids)) != len(ids) or not (len(ids) == len(cancer) == len(parts) == len(x)):
        raise BaselineExportError("patient IDs, cancer labels, splits, and WSI features must be unique and aligned")
    if x.ndim != 2 or not np.isfinite(x).all() or (y is not None and (y.ndim != 2 or len(y) != len(x) or not np.isfinite(y).all())):
        raise BaselineExportError("feature matrices must be finite aligned [patient, feature] arrays")
    if "train" not in set(parts) or "test" not in set(parts):
        raise BaselineExportError("baseline export needs explicit train and test patients")
    return ids, cancer, parts, x, y


def _write(output: str | Path, *, patient_ids, cancers, split, states: dict[str, np.ndarray], method: str,
           config: dict[str, object], source_provenance: dict[str, object] | None = None) -> Path:
    target = Path(output); target.parent.mkdir(parents=True, exist_ok=True)
    trained = np.asarray(trainable_state_names(*sorted(states)), dtype=str)
    manifest = {"artifact_version": ARTIFACT_VERSION, "method": method, "trained_states": trained.tolist(),
                "fit_population": config.get("fit_population", "none_no_fitted_representation"), "config": config,
                "source_provenance": source_provenance or {},
                "cohort_digest": sha256(json.dumps({"patient_ids": list(map(str, patient_ids)), "cancers": list(map(str, cancers)), "split": list(map(str, split))}, sort_keys=True).encode()).hexdigest()}
    np.savez_compressed(target, patient_ids=np.asarray(patient_ids, dtype=str), cancers=np.asarray(cancers, dtype=str),
                        split=np.asarray(split, dtype=str), trained_states=trained, artifact_version=np.asarray(ARTIFACT_VERSION),
                        manifest_json=np.asarray(json.dumps(manifest, sort_keys=True)), **states)
    return target


def export_raw_hoptimus_baseline(output: str | Path, *, patient_ids: Sequence[str], cancers: Sequence[str], split: Sequence[str],
                                 pooled_mean: np.ndarray, pooled_std: np.ndarray | None = None,
                                 source_provenance: dict[str, object] | None = None) -> Path:
    """Export raw frozen H-Optimus pooled representations with no fitted step."""
    ids, cancer, parts, mean, _ = _validate_inputs(patient_ids, cancers, split, pooled_mean)
    state = mean if pooled_std is None else np.concatenate([mean, np.asarray(pooled_std, dtype=np.float32)], axis=1)
    if pooled_std is not None and np.asarray(pooled_std).shape != mean.shape:
        raise BaselineExportError("pooled_std must have the same shape as pooled_mean")
    return _write(output, patient_ids=ids, cancers=cancer, split=parts, states={"wsi_identity": _normalise(state)},
                  method="raw_pooled_hoptimus", config={"pooling": "mean" if pooled_std is None else "mean_plus_std", "fit_population": "none_no_fitted_representation"},
                  source_provenance=source_provenance)


def export_zero_parameter_baseline(output: str | Path, *, patient_ids: Sequence[str], cancers: Sequence[str],
                                   split: Sequence[str], targets: np.ndarray,
                                   source_provenance: dict[str, object] | None = None) -> Path:
    """Export the cancer-type mean: the zero-fitted-parameter naive baseline.

    Every patient of cancer c receives the *training-fold* mean target vector for
    cancer c and nothing else. There is no fitted parameter beyond a group mean,
    and no patient-specific information whatsoever, so this is the floor any
    representation must clear before "it predicts molecular state" means anything.

    It is also the sharpest available test of the confound adjustment itself.
    Cancer type is in CALIBRA's design matrix, so after residualisation this
    baseline must collapse to nothing. If it does not, the adjustment is not
    removing the cancer effect it claims to remove, and every adjusted number in
    the run is reading lineage. A baseline that is *supposed* to score zero is
    worth more than one that is supposed to score high.

    Cancers absent from the training fold receive the global training mean rather
    than their own test-fold mean, which would be leakage dressed as a baseline.
    """
    ids, cancer, parts, y, _ = _validate_inputs(patient_ids, cancers, split, targets)
    train = parts == "train"
    if train.sum() < 2:
        raise BaselineExportError("the zero-parameter baseline needs at least two training patients")
    global_mean = y[train].mean(axis=0)
    state = np.repeat(global_mean[None, :], len(ids), axis=0).astype(np.float32)
    seen: list[str] = []
    for group in np.unique(cancer):
        rows = (cancer == group) & train
        if rows.sum() < 2:
            continue
        state[cancer == group] = y[rows].mean(axis=0)
        seen.append(str(group))
    # A representation that is constant everywhere is not measurable; that is a
    # property of the split, not a bug, and must surface as a refusal.
    if float(np.std(state, axis=0).max()) < 1e-8:
        raise BaselineExportError("cancer-type means are constant across the cohort; nothing to measure")
    return _write(output, patient_ids=ids, cancers=cancer, split=parts,
                  states={"wsi_identity": state},
                  method="zero_parameter_cancer_mean",
                  config={"rule": "training_fold_mean_target_vector_of_the_patient_cancer_type",
                          "n_cancers_with_training_mean": len(seen),
                          "fallback": "global_training_mean_for_cancers_absent_from_train",
                          "n_targets": int(y.shape[1]),
                          "fit_population": "train_fold_group_means_only_zero_fitted_parameters"},
                  source_provenance=source_provenance)


def _reduce_pair(wsi: np.ndarray, rna: np.ndarray, train: np.ndarray, components: int) -> tuple[np.ndarray, np.ndarray, int]:
    maximum = min(components, int(train.sum()) - 1, wsi.shape[1], rna.shape[1])
    if maximum < 2:
        raise BaselineExportError("at least three training patients and two shared components are required")
    # PCA is fit strictly on development patients; centering is retained in
    # sklearn's transform path for all val/test patients.
    px = PCA(n_components=maximum, whiten=True, random_state=0).fit(wsi[train])
    py = PCA(n_components=maximum, whiten=True, random_state=0).fit(rna[train])
    return px.transform(wsi), py.transform(rna), maximum


def _select_grouped_ridge_alpha(features: np.ndarray, target: np.ndarray, train: np.ndarray,
                                cancers: np.ndarray) -> float:
    """Select Ridge regularisation without mixing development cancer groups.

    ``RidgeCV`` defaults to ordinary sample folds.  That would be a different
    selection protocol from V2's cancer-grouped probes and can exaggerate a
    baseline's performance when lineage dominates.  The fallback is fixed and
    explicitly recorded for a synthetic/small cohort with fewer than two
    development cancer groups.
    """
    alphas = tuple(float(value) for value in np.logspace(-4, 4, 17))
    groups = cancers[train]
    if len(np.unique(groups)) < 2:
        return 1.0
    splitter = GroupKFold(n_splits=min(3, len(np.unique(groups))))
    scores: dict[float, list[float]] = {alpha: [] for alpha in alphas}
    for fit_rows, valid_rows in splitter.split(features[train], target[train], groups):
        x_fit, x_valid = features[train][fit_rows], features[train][valid_rows]
        y_fit, y_valid = target[train][fit_rows], target[train][valid_rows]
        for alpha in alphas:
            prediction = Ridge(alpha=alpha).fit(x_fit, y_fit).predict(x_valid)
            # CCA/Ridge alignment is multivariate; squared-error selection is
            # stable and does not inspect outer patients.
            scores[alpha].append(float(np.mean((prediction - y_valid) ** 2)))
    return min(alphas, key=lambda alpha: (float(np.mean(scores[alpha])), alpha))


def export_ridge_alignment_baseline(output: str | Path, *, patient_ids: Sequence[str], cancers: Sequence[str], split: Sequence[str],
                                    wsi: np.ndarray, rna: np.ndarray, components: int = 256,
                                    fit_development: bool = False,
                                    source_provenance: dict[str, object] | None = None) -> Path:
    """Fit development-only Ridge WSI→RNA alignment and export shared states."""
    ids, cancer, parts, x, y = _validate_inputs(patient_ids, cancers, split, wsi, rna)
    train = np.isin(parts, ("train", "val")) if fit_development else parts == "train"
    xr, yr, width = _reduce_pair(x, y, train, components)
    alpha = _select_grouped_ridge_alpha(xr, yr, train, cancer)
    ridge = Ridge(alpha=alpha).fit(xr[train], yr[train])
    states = {"wsi_identity": _normalise(ridge.predict(xr)), "rna_identity": _normalise(yr)}
    return _write(output, patient_ids=ids, cancers=cancer, split=parts, states=states, method="ridge_alignment",
                  config={"components": width, "alpha": float(alpha), "selection": "development_cancer_grouped_cv", "fit_population": "development_train_val" if fit_development else "train_only_inner_fit"}, source_provenance=source_provenance)


def export_cca_alignment_baseline(output: str | Path, *, patient_ids: Sequence[str], cancers: Sequence[str], split: Sequence[str],
                                  wsi: np.ndarray, rna: np.ndarray, components: int = 128, max_iter: int = 1000,
                                  fit_development: bool = False,
                                  source_provenance: dict[str, object] | None = None) -> Path:
    """Fit development-only PCA+CCA alignment and export canonical states."""
    ids, cancer, parts, x, y = _validate_inputs(patient_ids, cancers, split, wsi, rna)
    train = np.isin(parts, ("train", "val")) if fit_development else parts == "train"
    xr, yr, width = _reduce_pair(x, y, train, components)
    cca = CCA(n_components=width, max_iter=max_iter, scale=True).fit(xr[train], yr[train])
    wx, ry = cca.transform(xr, yr)
    states = {"wsi_identity": _normalise(wx), "rna_identity": _normalise(ry)}
    return _write(output, patient_ids=ids, cancers=cancer, split=parts, states=states, method="cca_alignment",
                  config={"components": width, "max_iter": max_iter, "fit_population": "development_train_val" if fit_development else "train_only_inner_fit"}, source_provenance=source_provenance)
