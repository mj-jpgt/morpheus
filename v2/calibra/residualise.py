"""Cross-fitted residualisation of confounds from both modalities.

Fold-safe by construction: the nuisance model is never fit on the rows whose
residuals it produces. This matters because in-sample residualisation removes
*more* than the confound (it eats real signal), which is precisely the
over-residualisation failure the CALIBRA instrument exists to quantify.
"""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

__all__ = ["confound_design", "cross_fitted_residuals"]


def confound_design(frame, columns) -> np.ndarray:
    """Build a numeric design matrix: categoricals one-hot, numerics standardised.

    NaNs in a numeric column become 0 after standardisation (mean imputation)
    plus an explicit missingness indicator column, so "missing" is itself
    adjustable rather than silently imputed.
    """
    import pandas as pd

    blocks: list[np.ndarray] = []
    for column in columns:
        values = frame[column]
        if values.dtype == object or str(values.dtype).startswith("category") or values.dtype.kind in "OUS":
            dummies = pd.get_dummies(values.astype(str), prefix=column, dummy_na=True)
            blocks.append(dummies.to_numpy(dtype=np.float64))
            continue
        numeric = values.to_numpy(dtype=np.float64)
        missing = ~np.isfinite(numeric)
        filled = np.where(missing, np.nan, numeric)
        mean = np.nanmean(filled) if np.isfinite(filled).any() else 0.0
        std = np.nanstd(filled) if np.isfinite(filled).any() else 0.0
        centred = (np.nan_to_num(filled, nan=mean) - mean) / (std if std > 1e-12 else 1.0)
        blocks.append(centred[:, None])
        if missing.any():
            blocks.append(missing.astype(np.float64)[:, None])
    if not blocks:
        return np.zeros((len(frame), 0))
    return np.hstack(blocks)


def cross_fitted_residuals(matrix: np.ndarray, design: np.ndarray, *, n_splits: int = 5,
                           alpha: float = 1.0, seed: int = 42) -> np.ndarray:
    """Return ``matrix`` with ``design`` regressed out, cross-fitted over folds."""
    matrix = np.asarray(matrix, dtype=np.float64)
    if design.shape[1] == 0:
        return matrix - matrix.mean(axis=0, keepdims=True)
    residual = np.empty_like(matrix)
    splitter = KFold(n_splits=min(n_splits, len(matrix)), shuffle=True, random_state=seed)
    for train_idx, test_idx in splitter.split(matrix):
        model = Ridge(alpha=alpha, fit_intercept=True)
        model.fit(design[train_idx], matrix[train_idx])
        residual[test_idx] = matrix[test_idx] - model.predict(design[test_idx])
    return residual - residual.mean(axis=0, keepdims=True)
