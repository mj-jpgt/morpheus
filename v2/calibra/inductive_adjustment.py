"""An INDUCTIVE confound-adjustment operator: fit once, apply to unseen slides.

The problem this solves. :func:`~morpheus.v2.calibra.residualise.cross_fitted_residuals`
is **transductive**. Three separate things in it are properties of the cohort
you hand it, not of a fitted map:

  1. the K nuisance models are fit inside the call and thrown away;
  2. which model produced a given row's residual depends on a ``KFold`` shuffle
     over ``len(matrix)``, so it is not defined for a row that was not there;
  3. the final ``residual - residual.mean(axis=0)`` re-centres on the cohort,
     so adding one row changes every other row's coordinates.

The consequence is that there is no way to give a *new* slide its adjusted
coordinates, and therefore no way to state what a deployed system would do with
one. It also blocks the only fair way to compare two cohorts: adjusting
ALCHEMIST with **TCGA's** operator rather than with its own, so that the two sit
in the same coordinate system rather than in two separately-fitted ones.

What this module adds. :class:`ConfoundAdjustmentOperator` persists the whole
map -- the design encoding, the site-pooling rule, the K fitted ridge models,
the fold assignment, and the cohort centre -- with the provenance needed to
trace an applied adjustment back to what fitted it.

**Identity on the fitting cohort is the load-bearing property.** ``adjust_reference``
reproduces ``cross_fitted_residuals`` bit-for-bit (asserted in
``v2/tests/test_inductive_adjustment.py``), because every published CALIBRA
number was produced by the transductive function and a "nearly identical"
operator would silently make all of them non-comparable. That is why the
reference rows keep their own fold's model rather than being pushed through the
new-row path.

**Reference rows and new rows are adjusted by different rules, deliberately.**
A reference row uses the one model that did not see it (cross-fitting). A new
row was seen by none of them, so it uses the mean of all K predictions. The two
therefore do not agree exactly, and a slide's coordinates depend on whether it
was in the reference cohort. That discontinuity is real and is asserted by a
test rather than hidden: the alternatives are to break identity on the fitting
cohort (which invalidates the published numbers) or to refit on all rows (which
introduces a model no reference row was ever scored by, and which no identity
test can validate). The K-model ensemble is the only option that both preserves
identity and uses exclusively models that are already fitted and persisted.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

from .residualise import apply_pooled_tissue_source_site, pooled_tissue_source_site

__all__ = ["DesignSpec", "SitePooling", "InductiveResidualiser", "ConfoundAdjustmentOperator",
           "UnseenLevelError"]

#: Written into the persisted file so a future reader can tell which contract it holds.
OPERATOR_FORMAT_VERSION = 1


class UnseenLevelError(ValueError):
    """A new cohort carries a categorical level the reference operator never saw.

    Raised rather than silently zero-encoded. A zero row in a one-hot block is
    not "no confound", it is "the reference's baseline level", so quietly
    encoding an unknown cancer type as zeros would adjust that slide as though
    it belonged to whichever level the design dropped. Site codes never reach
    this, because :class:`SitePooling` maps them to ``OTHER`` first.
    """


# --------------------------------------------------------------------------- site pooling

@dataclass(frozen=True)
class SitePooling:
    """The reference cohort's frozen tissue-source-site pooling rule.

    ``pooled_tissue_source_site`` pools every site with fewer than
    ``min_site_count`` patients into ``OTHER``. Freezing ``frequent`` is what
    makes that rule inductive: a new patient's site is kept iff it was frequent
    **in the reference**, and is otherwise ``OTHER``.

    This is also the answer to the unseen-site question, and it is not a new
    policy -- it is the existing one evaluated at a count of zero. A site the
    reference never saw has reference count 0 < ``min_site_count``, so the
    reference's own rule already says ``OTHER``. Refusing such a row would make
    the operator unusable on exactly the cohorts it exists for (ALCHEMIST shares
    almost no TSS codes with TCGA), and a nearest-site rule would require a
    similarity between site codes that does not exist -- TSS codes are
    registry identifiers, not coordinates.
    """
    column: str
    min_site_count: int
    frequent: tuple[str, ...]

    @classmethod
    def fit(cls, patient_ids, *, column: str = "tss", min_site_count: int = 10) -> "SitePooling":
        _, frequent = pooled_tissue_source_site(patient_ids, min_site_count=min_site_count)
        return cls(column=column, min_site_count=int(min_site_count),
                   frequent=tuple(sorted(frequent)))

    def apply(self, patient_ids) -> tuple[np.ndarray, dict]:
        """Pooled site labels, plus a report naming every site that was pooled."""
        from .residualise import tissue_source_site

        raw = tissue_source_site(patient_ids)
        pooled = apply_pooled_tissue_source_site(patient_ids, set(self.frequent))
        unseen = sorted({site for site in raw.tolist() if site not in set(self.frequent)})
        report = {
            "n_rows": int(len(raw)),
            "n_pooled_to_other": int(np.sum(pooled == "OTHER")),
            "sites_pooled_to_other": unseen,
            "n_sites_unknown_to_reference": len(unseen),
            "policy": "unseen site -> OTHER (the reference pooling rule at a count of zero)",
        }
        return pooled, report


# --------------------------------------------------------------------------- design

@dataclass(frozen=True)
class DesignSpec:
    """``confound_design`` frozen against a reference cohort, so new rows land in the
    SAME columns in the SAME order.

    Two things in ``confound_design`` are cohort statistics and must be frozen
    or the new rows are encoded on a different scale:

    * a categorical column's one-hot levels (``pd.get_dummies`` sorts whatever
      levels are present, so a new cohort would otherwise produce different
      columns in a different order); and
    * a numeric column's mean and s.d., which ``confound_design`` takes from the
      frame it is given. Standardising a new slide by *its own cohort's* mean
      would put it in a different coordinate system from the reference -- the
      exact error this module exists to prevent.

    The missingness-indicator column is likewise a property of the reference: it
    is emitted iff the reference column had any missing value, whatever the new
    cohort looks like.
    """
    columns: tuple[str, ...]
    kinds: tuple[str, ...]                              # "categorical" | "numeric"
    categorical_columns: dict = field(default_factory=dict)   # source col -> dummy col names
    numeric_stats: dict = field(default_factory=dict)         # source col -> (mean, std, has_missing)
    design_columns: tuple[str, ...] = ()

    @classmethod
    def fit(cls, frame, columns) -> "DesignSpec":
        import pandas as pd

        kinds, categorical, numeric, design_columns = [], {}, {}, []
        for column in columns:
            values = frame[column]
            if (values.dtype == object or str(values.dtype).startswith("category")
                    or values.dtype.kind in "OUS"):
                dummies = pd.get_dummies(values.astype(str), prefix=column, dummy_na=True)
                kinds.append("categorical")
                categorical[column] = tuple(str(c) for c in dummies.columns)
                design_columns.extend(categorical[column])
                continue
            numeric_values = values.to_numpy(dtype=np.float64)
            missing = ~np.isfinite(numeric_values)
            filled = np.where(missing, np.nan, numeric_values)
            mean = float(np.nanmean(filled)) if np.isfinite(filled).any() else 0.0
            std = float(np.nanstd(filled)) if np.isfinite(filled).any() else 0.0
            kinds.append("numeric")
            numeric[column] = (mean, std, bool(missing.any()))
            design_columns.append(column)
            if missing.any():
                design_columns.append(f"{column}__missing")
        return cls(columns=tuple(str(c) for c in columns), kinds=tuple(kinds),
                   categorical_columns=categorical, numeric_stats=numeric,
                   design_columns=tuple(design_columns))

    def transform(self, frame, *, on_unseen_level: str = "refuse") -> tuple[np.ndarray, dict]:
        """Encode ``frame`` into the reference's design columns.

        ``on_unseen_level`` is ``"refuse"`` (default) or ``"zero"``. Refusing is
        the default because a zero one-hot row is not a neutral encoding, it is
        the reference's implicit baseline level; a caller that genuinely wants
        that has to ask for it, and the choice is recorded in the report.
        """
        import pandas as pd

        if on_unseen_level not in ("refuse", "zero"):
            raise ValueError("on_unseen_level must be 'refuse' or 'zero'")
        blocks, unseen_report = [], {}
        for column, kind in zip(self.columns, self.kinds):
            values = frame[column]
            if kind == "categorical":
                expected = list(self.categorical_columns[column])
                dummies = pd.get_dummies(values.astype(str), prefix=column, dummy_na=True)
                extra = [str(c) for c in dummies.columns if str(c) not in set(expected)]
                if extra:
                    unseen_report[column] = sorted(
                        name[len(column) + 1:] for name in extra)
                    if on_unseen_level == "refuse":
                        raise UnseenLevelError(
                            f"column {column!r} carries levels absent from the reference "
                            f"operator: {unseen_report[column]}. Pass on_unseen_level='zero' "
                            f"to encode them as the reference baseline, or pool them first.")
                aligned = dummies.reindex(columns=expected, fill_value=0)
                blocks.append(aligned.to_numpy(dtype=np.float64))
                continue
            mean, std, has_missing = self.numeric_stats[column]
            numeric_values = values.to_numpy(dtype=np.float64)
            missing = ~np.isfinite(numeric_values)
            filled = np.where(missing, np.nan, numeric_values)
            centred = (np.nan_to_num(filled, nan=mean) - mean) / (std if std > 1e-12 else 1.0)
            blocks.append(centred[:, None])
            if has_missing:
                blocks.append(missing.astype(np.float64)[:, None])
        design = np.hstack(blocks) if blocks else np.zeros((len(frame), 0))
        report = {"n_design_columns": int(design.shape[1]),
                  "unseen_levels": unseen_report,
                  "on_unseen_level": on_unseen_level}
        return design, report


# --------------------------------------------------------------------------- residualiser

@dataclass
class InductiveResidualiser:
    """The K cross-fitting ridge models, the fold map, and the cohort centre.

    ``fold_of_row`` and ``centre`` are what make the reference path exact:
    ``cross_fitted_residuals`` gives row *i* the model fit on the folds *not*
    containing *i*, and then subtracts the residual matrix's own column mean.
    Recomputing either from a new cohort would move every published number.
    """
    n_splits: int
    alpha: float
    seed: int
    n_reference_rows: int
    n_design_columns: int
    fold_of_row: np.ndarray
    centre: np.ndarray
    coefficients: list = field(default_factory=list)     # per fold, (k_out, d)
    intercepts: list = field(default_factory=list)       # per fold, (k_out,)

    @classmethod
    def fit(cls, matrix: np.ndarray, design: np.ndarray, *, n_splits: int = 5,
            alpha: float = 1.0, seed: int = 42) -> "InductiveResidualiser":
        matrix = np.asarray(matrix, dtype=np.float64)
        design = np.asarray(design, dtype=np.float64)
        n = len(matrix)
        if design.shape[1] == 0:
            # cross_fitted_residuals' degenerate branch: pure column centring.
            centre = matrix.mean(axis=0)
            return cls(n_splits=0, alpha=float(alpha), seed=int(seed), n_reference_rows=n,
                       n_design_columns=0, fold_of_row=np.full(n, -1, dtype=np.int64),
                       centre=centre)

        effective = min(n_splits, n)
        splitter = KFold(n_splits=effective, shuffle=True, random_state=seed)
        residual = np.empty_like(matrix)
        fold_of_row = np.full(n, -1, dtype=np.int64)
        coefficients, intercepts = [], []
        for fold, (train_idx, test_idx) in enumerate(splitter.split(matrix)):
            model = Ridge(alpha=alpha, fit_intercept=True)
            model.fit(design[train_idx], matrix[train_idx])
            prediction = np.asarray(model.predict(design[test_idx]))
            residual[test_idx] = matrix[test_idx] - prediction.reshape(matrix[test_idx].shape)
            fold_of_row[test_idx] = fold
            coefficients.append(np.atleast_2d(model.coef_).copy())
            intercepts.append(np.atleast_1d(model.intercept_).copy())
        return cls(n_splits=int(effective), alpha=float(alpha), seed=int(seed),
                   n_reference_rows=n, n_design_columns=int(design.shape[1]),
                   fold_of_row=fold_of_row, centre=residual.mean(axis=0),
                   coefficients=coefficients, intercepts=intercepts)

    def _model(self, fold: int) -> Ridge:
        """Rebuild fold ``fold``'s estimator.

        The prediction goes back through ``Ridge.predict`` rather than through a
        hand-written ``design @ coef.T + intercept`` so that the reference path
        is bit-identical to ``cross_fitted_residuals``, which is what the
        identity test asserts.
        """
        model = Ridge(alpha=self.alpha, fit_intercept=True)
        model.coef_ = self.coefficients[fold]
        model.intercept_ = self.intercepts[fold]
        model.n_features_in_ = int(self.n_design_columns)
        return model

    def _predict(self, design: np.ndarray, fold: int, n_outputs: int) -> np.ndarray:
        prediction = np.asarray(self._model(fold).predict(design))
        return prediction.reshape(len(design), n_outputs)

    def transform_reference(self, matrix: np.ndarray, design: np.ndarray) -> np.ndarray:
        """Reproduce ``cross_fitted_residuals(matrix, design, ...)`` exactly."""
        matrix = np.asarray(matrix, dtype=np.float64)
        design = np.asarray(design, dtype=np.float64)
        if len(matrix) != self.n_reference_rows:
            raise ValueError(f"reference path expects the {self.n_reference_rows} fitting rows "
                             f"in their fitting order, got {len(matrix)}")
        if self.n_design_columns == 0:
            return matrix - self.centre
        residual = np.empty_like(matrix)
        for fold in range(self.n_splits):
            test_idx = np.flatnonzero(self.fold_of_row == fold)
            residual[test_idx] = matrix[test_idx] - self._predict(design[test_idx], fold,
                                                                  matrix.shape[1])
        return residual - self.centre

    def transform(self, matrix: np.ndarray, design: np.ndarray) -> np.ndarray:
        """Adjust rows the operator never saw, using the mean of the K fold models."""
        matrix = np.asarray(matrix, dtype=np.float64)
        design = np.asarray(design, dtype=np.float64)
        if self.n_design_columns == 0:
            return matrix - self.centre
        if design.shape[1] != self.n_design_columns:
            raise ValueError(f"design has {design.shape[1]} columns, operator fitted on "
                             f"{self.n_design_columns}")
        prediction = np.mean([self._predict(design, fold, matrix.shape[1])
                              for fold in range(self.n_splits)], axis=0)
        return matrix - prediction - self.centre


# --------------------------------------------------------------------------- operator

def _digest(*parts) -> str:
    """A stable content digest over the reference cohort's identity."""
    hasher = hashlib.sha256()
    for part in parts:
        array = np.ascontiguousarray(part)
        hasher.update(str(array.dtype).encode())
        hasher.update(str(array.shape).encode())
        hasher.update(array.tobytes() if array.dtype.kind != "U"
                      else "\x1f".join(array.ravel().tolist()).encode())
    return hasher.hexdigest()


@dataclass
class ConfoundAdjustmentOperator:
    """A fitted, persistable confound->representation map.

    Fit it on a reference cohort, save it, and apply it to slides that cohort
    never contained. See the module docstring for the identity guarantee and for
    why reference and new rows take different paths.
    """
    design_spec: DesignSpec
    residualiser: InductiveResidualiser
    site_pooling: SitePooling | None = None
    provenance: dict = field(default_factory=dict)

    # -- fitting ------------------------------------------------------------------
    @classmethod
    def fit(cls, matrix: np.ndarray, frame, columns, *, patient_ids=None,
            site_column: str = "", min_site_count: int = 10, n_splits: int = 5,
            alpha: float = 1.0, seed: int = 42, cohort_name: str = "") -> "ConfoundAdjustmentOperator":
        """Fit on a reference cohort.

        ``frame`` is the confound table and ``columns`` the design columns, exactly
        as they would be passed to ``confound_design``. If ``site_column`` is given,
        ``patient_ids`` must be too, and that column of ``frame`` is (re)built from
        the barcodes through the frozen pooling rule -- so the operator owns the
        pooling rather than trusting the caller to reapply it consistently.
        """
        import pandas as pd

        frame = pd.DataFrame(frame).reset_index(drop=True)
        pooling = None
        if site_column:
            if patient_ids is None:
                raise ValueError("site_column requires patient_ids")
            pooling = SitePooling.fit(patient_ids, column=site_column,
                                      min_site_count=min_site_count)
            frame[site_column], _ = pooling.apply(patient_ids)

        spec = DesignSpec.fit(frame, columns)
        design, _ = spec.transform(frame, on_unseen_level="refuse")
        residualiser = InductiveResidualiser.fit(matrix, design, n_splits=n_splits,
                                                 alpha=alpha, seed=seed)
        identifiers = (np.asarray(patient_ids).astype(str) if patient_ids is not None
                       else np.asarray([], dtype="U1"))
        provenance = {
            "format_version": OPERATOR_FORMAT_VERSION,
            "cohort_name": cohort_name,
            "reference_digest": _digest(design, identifiers,
                                        np.asarray(matrix.shape, dtype=np.int64)),
            "reference_design_digest": _digest(design),
            "reference_patient_digest": _digest(identifiers),
            "n_reference_rows": int(len(matrix)),
            "n_representation_columns": int(np.asarray(matrix).shape[1]),
            "source_columns": list(spec.columns),
            "design_columns": list(spec.design_columns),
            "n_design_columns": int(design.shape[1]),
            "site_pooling": (None if pooling is None else
                             {"column": pooling.column, "min_site_count": pooling.min_site_count,
                              "n_frequent_sites": len(pooling.frequent),
                              "frequent_sites": list(pooling.frequent)}),
            "residualiser": {"n_splits": residualiser.n_splits, "alpha": residualiser.alpha,
                             "seed": residualiser.seed},
            "numpy_version": np.__version__,
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        try:                                             # provenance only, never numerics
            import sklearn
            provenance["sklearn_version"] = sklearn.__version__
        except Exception:                                # pragma: no cover
            provenance["sklearn_version"] = "unknown"
        return cls(design_spec=spec, residualiser=residualiser, site_pooling=pooling,
                   provenance=provenance)

    # -- applying -----------------------------------------------------------------
    def _frame_and_design(self, frame, patient_ids, on_unseen_level: str):
        import pandas as pd

        frame = pd.DataFrame(frame).reset_index(drop=True)
        report = {}
        if self.site_pooling is not None:
            if patient_ids is None:
                raise ValueError("this operator pools sites and therefore requires patient_ids")
            frame[self.site_pooling.column], report["site_pooling"] = \
                self.site_pooling.apply(patient_ids)
        design, design_report = self.design_spec.transform(frame, on_unseen_level=on_unseen_level)
        report.update(design_report)
        return design, report

    def adjust_reference(self, matrix: np.ndarray, frame, *, patient_ids=None) -> np.ndarray:
        """Adjust the fitting cohort. Bit-identical to ``cross_fitted_residuals``."""
        design, _ = self._frame_and_design(frame, patient_ids, "refuse")
        return self.residualiser.transform_reference(matrix, design)

    def adjust(self, matrix: np.ndarray, frame, *, patient_ids=None,
               on_unseen_level: str = "refuse") -> tuple[np.ndarray, dict]:
        """Adjust unseen rows. Returns ``(adjusted, report)``.

        The report names every site pooled to ``OTHER`` and every refused or
        zero-encoded level, so an applied adjustment can always say what it did
        to a row it had not seen before.
        """
        design, report = self._frame_and_design(frame, patient_ids, on_unseen_level)
        adjusted = self.residualiser.transform(matrix, design)
        report["applied_by"] = {"reference_digest": self.provenance.get("reference_digest"),
                                "cohort_name": self.provenance.get("cohort_name"),
                                "n_reference_rows": self.provenance.get("n_reference_rows")}
        report["n_rows_adjusted"] = int(len(matrix))
        report["path"] = "inductive_fold_ensemble"
        return adjusted, report

    # -- persistence --------------------------------------------------------------
    def save(self, path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "meta": json.dumps({
                "provenance": self.provenance,
                "design_spec": {
                    "columns": list(self.design_spec.columns),
                    "kinds": list(self.design_spec.kinds),
                    "categorical_columns": {k: list(v) for k, v
                                            in self.design_spec.categorical_columns.items()},
                    "numeric_stats": {k: list(v) for k, v
                                      in self.design_spec.numeric_stats.items()},
                    "design_columns": list(self.design_spec.design_columns),
                },
                "site_pooling": (None if self.site_pooling is None else {
                    "column": self.site_pooling.column,
                    "min_site_count": self.site_pooling.min_site_count,
                    "frequent": list(self.site_pooling.frequent)}),
                "residualiser": {
                    "n_splits": self.residualiser.n_splits, "alpha": self.residualiser.alpha,
                    "seed": self.residualiser.seed,
                    "n_reference_rows": self.residualiser.n_reference_rows,
                    "n_design_columns": self.residualiser.n_design_columns,
                    "n_folds_stored": len(self.residualiser.coefficients)},
            }),
            "fold_of_row": self.residualiser.fold_of_row,
            "centre": self.residualiser.centre,
        }
        for fold, (coef, intercept) in enumerate(zip(self.residualiser.coefficients,
                                                     self.residualiser.intercepts)):
            payload[f"coef_{fold}"] = coef
            payload[f"intercept_{fold}"] = intercept
        np.savez(path, **payload)
        return path

    @classmethod
    def load(cls, path) -> "ConfoundAdjustmentOperator":
        raw = np.load(Path(path), allow_pickle=False)
        meta = json.loads(str(raw["meta"]))
        spec_meta = meta["design_spec"]
        spec = DesignSpec(
            columns=tuple(spec_meta["columns"]), kinds=tuple(spec_meta["kinds"]),
            categorical_columns={k: tuple(v) for k, v in spec_meta["categorical_columns"].items()},
            numeric_stats={k: (float(v[0]), float(v[1]), bool(v[2]))
                           for k, v in spec_meta["numeric_stats"].items()},
            design_columns=tuple(spec_meta["design_columns"]))
        residual_meta = meta["residualiser"]
        residualiser = InductiveResidualiser(
            n_splits=int(residual_meta["n_splits"]), alpha=float(residual_meta["alpha"]),
            seed=int(residual_meta["seed"]),
            n_reference_rows=int(residual_meta["n_reference_rows"]),
            n_design_columns=int(residual_meta["n_design_columns"]),
            fold_of_row=np.asarray(raw["fold_of_row"]), centre=np.asarray(raw["centre"]),
            coefficients=[np.asarray(raw[f"coef_{i}"])
                          for i in range(int(residual_meta["n_folds_stored"]))],
            intercepts=[np.asarray(raw[f"intercept_{i}"])
                        for i in range(int(residual_meta["n_folds_stored"]))])
        pooling_meta = meta["site_pooling"]
        pooling = (None if pooling_meta is None else
                   SitePooling(column=pooling_meta["column"],
                               min_site_count=int(pooling_meta["min_site_count"]),
                               frequent=tuple(pooling_meta["frequent"])))
        return cls(design_spec=spec, residualiser=residualiser, site_pooling=pooling,
                   provenance=meta["provenance"])
