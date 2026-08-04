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

__all__ = ["confound_design", "cross_fitted_residuals", "pooled_tissue_source_site",
           "tissue_source_site", "apply_pooled_tissue_source_site"]


def tissue_source_site(patient_ids: np.ndarray) -> np.ndarray:
    """The raw TCGA TSS code (second barcode field), before any pooling.

    Split out so that the transductive and inductive pooling rules cannot drift
    apart in how they parse a barcode.
    """
    identifiers = np.asarray(patient_ids).astype(str)
    return np.asarray([identifier.split("-")[1] if len(identifier.split("-")) > 1 else "NA"
                       for identifier in identifiers])


def pooled_tissue_source_site(patient_ids: np.ndarray, *, min_site_count: int = 10) -> tuple[np.ndarray, set[str]]:
    """Derive TCGA TSS and pool rare sites exactly once across CALIBRA analyses."""
    raw = tissue_source_site(patient_ids)
    unique, counts = np.unique(raw, return_counts=True)
    frequent = {site for site, count in zip(unique, counts) if count >= min_site_count}
    return apply_pooled_tissue_source_site(patient_ids, frequent), frequent


def apply_pooled_tissue_source_site(patient_ids: np.ndarray, frequent) -> np.ndarray:
    """Apply an ALREADY-FITTED pooling rule to (possibly new) patients.

    A site absent from ``frequent`` becomes ``OTHER``. That covers a site that
    was rare in the reference cohort AND a site the reference cohort never saw
    at all -- the two are the same case, because an unseen site is a site whose
    reference count is zero, and zero is below every ``min_site_count``. So the
    unseen-site policy is not a new rule invented for the inductive path: it is
    :func:`pooled_tissue_source_site`'s own rule evaluated at a count of zero.

    Applying it to the reference cohort's own patients reproduces
    ``pooled_tissue_source_site(patient_ids, min_site_count=...)[0]`` exactly,
    which is what makes the inductive design identical on the fitting cohort.
    """
    frequent = set(frequent)
    return np.asarray([site if site in frequent else "OTHER"
                       for site in tissue_source_site(patient_ids)])


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
        # sklearn ravels a SINGLE-column target: Ridge.predict returns (n,) for a
        # (n, 1) matrix and (n, k) for k >= 2. Without this reshape the one-column
        # case broadcasts to (n, n) and raises -- which is why per-AXIS residualisation
        # (one column at a time, the unit P4 certification works in) was unreachable.
        # A strict no-op for every k >= 2 call site, so no existing number moves.
        prediction = np.asarray(model.predict(design[test_idx]))
        residual[test_idx] = matrix[test_idx] - prediction.reshape(matrix[test_idx].shape)
    return residual - residual.mean(axis=0, keepdims=True)
