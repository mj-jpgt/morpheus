"""Leakage-safe construction of RNA-derived discovery targets.

This module deliberately has no dependency on the MORPHEUS model or evaluator.
It turns a patient-by-gene RNA table and a GMT signature collection into a
canonical patient-aligned target bundle.  The bundle retains raw scores,
per-target availability, and a train-only fitted standardization so downstream
evaluators cannot accidentally use an outer-test cancer to fit a transform.

The expected RNA layout is *one row per patient* and *one gene per column*.
Patient identifiers can be supplied in a column or in the dataframe index.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd


class TargetBuildError(ValueError):
    """Raised when a target table cannot be built without ambiguity."""


def _require_unique(frame: pd.DataFrame, column: str, label: str) -> None:
    if column not in frame.columns:
        raise TargetBuildError(f"{label} has no {column!r} column")
    duplicated = frame[column].duplicated(keep=False)
    if duplicated.any():
        examples = frame.loc[duplicated, column].astype(str).head(5).tolist()
        raise TargetBuildError(f"{label} contains duplicate patient IDs: {examples}")


def _frame_digest(frame: pd.DataFrame) -> str:
    """Return a stable digest suitable for a provenance manifest."""
    stable = frame.copy()
    stable.columns = stable.columns.astype(str)
    payload = stable.to_csv(index=False, na_rep="<NA>", lineterminator="\n")
    return sha256(payload.encode("utf-8")).hexdigest()


def read_rna_table(
    source: str | Path | pd.DataFrame,
    *,
    patient_id_col: str = "patient_id",
) -> pd.DataFrame:
    """Read a wide patient-by-gene RNA table from TSV, CSV, parquet, or memory.

    If ``patient_id_col`` is absent but the index has a name, the index is
    promoted to that column.  Gene columns are retained verbatim; callers may
    normalise gene symbols before this function if their source needs it.
    """
    if isinstance(source, pd.DataFrame):
        frame = source.copy()
    else:
        path = Path(source)
        suffix = path.suffix.lower()
        if suffix in {".parquet", ".pq"}:
            frame = pd.read_parquet(path)
        elif suffix in {".tsv", ".txt", ".tab"}:
            frame = pd.read_csv(path, sep="\t")
        elif suffix == ".csv":
            frame = pd.read_csv(path)
        else:
            raise TargetBuildError(
                f"Unsupported RNA input {path}; expected TSV, CSV, or parquet"
            )

    if patient_id_col not in frame.columns and frame.index.name is not None:
        frame = frame.reset_index().rename(columns={frame.index.name: patient_id_col})
    _require_unique(frame, patient_id_col, "RNA table")
    if frame.shape[1] < 2:
        raise TargetBuildError("RNA table must contain a patient ID and at least one gene")
    frame[patient_id_col] = frame[patient_id_col].astype(str)
    return frame


def read_gmt(source: str | Path | Mapping[str, Sequence[str]]) -> dict[str, tuple[str, ...]]:
    """Read a GMT file or normalise an in-memory signature mapping.

    GMT's second field is a free-text description and is intentionally not used
    as a gene. Duplicate genes are removed while preserving order.
    """
    if isinstance(source, Mapping):
        raw = source
    else:
        raw: dict[str, Sequence[str]] = {}
        with Path(source).open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 3 or not parts[0].strip():
                    raise TargetBuildError(f"Malformed GMT row {line_number}")
                if parts[0] in raw:
                    raise TargetBuildError(f"Duplicate GMT signature {parts[0]!r}")
                raw[parts[0]] = parts[2:]

    normalised: dict[str, tuple[str, ...]] = {}
    for name, genes in raw.items():
        clean_name = str(name).strip()
        clean_genes = tuple(dict.fromkeys(str(g).strip() for g in genes if str(g).strip()))
        if not clean_name or not clean_genes:
            raise TargetBuildError(f"Signature {name!r} has no usable genes")
        if clean_name in normalised:
            raise TargetBuildError(f"Duplicate signature {clean_name!r}")
        normalised[clean_name] = clean_genes
    if not normalised:
        raise TargetBuildError("No signatures were supplied")
    return normalised


def score_gene_signatures(
    rna: pd.DataFrame,
    signatures: Mapping[str, Sequence[str]],
    *,
    patient_id_col: str = "patient_id",
    min_present_genes: int = 1,
    score_mode: str = "competitive_rank",
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    """Compute patient-level signature scores without cross-patient fitting.

    ``competitive_rank`` is the default because it asks whether a set is high
    *relative to the rest of that patient's transcriptome*.  It is therefore
    less confounded by library-wide expression shifts than a raw mean.  The
    score is the mean within-sample percentile rank of member genes minus the
    mean percentile rank of all non-member genes.  It is deliberately a
    deterministic, dependency-light ssGSEA-family score rather than a claim
    that raw average expression is pathway activity.  ``mean_expression`` is
    retained only for reproducing legacy artifacts.

    A signature with fewer than ``min_present_genes`` in the RNA matrix is
    represented as all ``NaN`` and false in the availability table.  Individual
    patient scores are also unavailable if every contributing gene is missing.
    This per-target representation prevents one missing programme from deleting
    otherwise evaluable patients.
    """
    if min_present_genes < 1:
        raise TargetBuildError("min_present_genes must be at least one")
    if score_mode not in {"competitive_rank", "mean_expression"}:
        raise TargetBuildError("score_mode must be 'competitive_rank' or 'mean_expression'")
    rna = read_rna_table(rna, patient_id_col=patient_id_col)
    signatures = read_gmt(signatures)
    gene_columns = [column for column in rna.columns if column != patient_id_col]
    genes = set(gene_columns)
    expression = rna.loc[:, gene_columns].apply(pd.to_numeric, errors="coerce")
    # percentile=True makes each sample comparable despite differing numbers
    # of measured genes.  NaNs remain unavailable rather than becoming zeros.
    percentile_ranks = expression.rank(axis=1, method="average", pct=True)
    scores = pd.DataFrame({patient_id_col: rna[patient_id_col].astype(str)})
    availability = pd.DataFrame({patient_id_col: rna[patient_id_col].astype(str)})
    matched_genes: dict[str, list[str]] = {}

    for name, signature_genes in signatures.items():
        present = [gene for gene in signature_genes if gene in genes]
        matched_genes[name] = present
        if len(present) < min_present_genes:
            scores[name] = np.nan
            availability[name] = False
            continue
        if score_mode == "mean_expression":
            score = expression.loc[:, present].mean(axis=1, skipna=True)
        else:
            member = percentile_ranks.loc[:, present]
            # Use the whole measured transcriptome as the competitive
            # reference.  This stays defined when a small smoke matrix happens
            # to contain every member gene; in that case the score is exactly
            # zero rather than incorrectly becoming unavailable.
            background = percentile_ranks
            score = member.mean(axis=1, skipna=True) - background.mean(axis=1, skipna=True)
        scores[name] = score.astype(float)
        availability[name] = score.notna()
    return scores, availability, matched_genes


def build_matched_random_controls(
    rna_source: str | Path | pd.DataFrame,
    biological_signatures: Mapping[str, Sequence[str]],
    metadata: pd.DataFrame,
    *,
    patient_id_col: str = "patient_id",
    train_mask: Sequence[bool] | pd.Series | None = None,
    split_col: str = "split",
    train_label: str = "train",
    controls_per_target: int = 20,
    seed: int = 42,
) -> tuple[dict[str, tuple[str, ...]], dict[str, object]]:
    """Create deterministic controls matched on *training-only* gene covariates.

    Candidate genes exclude every biological-panel member.  Each sampled set
    matches its target's per-gene mean expression, variance, and first-PC
    loading using nearest-neighbour cost in robustly standardised covariate
    space.  This is deliberately stricter than random same-size controls and
    never looks at outer-test RNA values.
    """
    if controls_per_target < 1:
        raise TargetBuildError("controls_per_target must be positive")
    rna = read_rna_table(rna_source, patient_id_col=patient_id_col)
    n_source_rna_patients = int(len(rna))
    _require_unique(metadata, patient_id_col, "metadata")
    metadata_index = metadata.copy(); metadata_index[patient_id_col] = metadata_index[patient_id_col].astype(str)
    # The PanCan RNA source deliberately contains a broader universe than a
    # particular WSI representation artifact.  Matched-control covariates must
    # be estimated only on the artifact's canonical cohort; treating unrelated
    # RNA-only patients as an alignment error both rejects valid inputs and
    # obscures the actual fit population.  The intersection is deterministic
    # and is applied before the development-train mask is constructed.
    canonical_patients = set(metadata_index[patient_id_col])
    rna = rna.loc[rna[patient_id_col].astype(str).isin(canonical_patients)].copy()
    if rna.empty:
        raise TargetBuildError("no RNA patients overlap the canonical metadata cohort")
    aligned = metadata_index.set_index(patient_id_col).reindex(rna[patient_id_col].astype(str))
    if train_mask is None:
        if split_col not in aligned:
            raise TargetBuildError("Provide train_mask or metadata containing split")
        train = aligned[split_col].astype(str).str.lower().eq(train_label.lower()).to_numpy()
    else:
        train = np.asarray(train_mask, dtype=bool)
    if train.shape != (len(rna),) or train.sum() < 3:
        raise TargetBuildError("matched controls require at least three training RNA rows")
    expression = rna.drop(columns=[patient_id_col]).apply(pd.to_numeric, errors="coerce")
    genes = expression.columns.astype(str).tolist()
    values = expression.to_numpy(float)[train]
    mean = np.nanmean(values, axis=0); variance = np.nanvar(values, axis=0)
    centered = np.where(np.isfinite(values), values, mean) - mean
    scale = np.maximum(np.sqrt(variance), 1e-8)
    # The first right singular vector provides a deterministic co-expression
    # descriptor without fitting anything on held-out patients.
    _, _, vt = np.linalg.svd(centered / scale, full_matrices=False)
    pc1 = vt[0] if len(vt) else np.zeros(len(genes))
    covariates = np.column_stack([mean, np.log1p(variance), pc1])
    median = np.nanmedian(covariates, axis=0); mad = np.nanmedian(np.abs(covariates - median), axis=0)
    standard = (covariates - median) / np.maximum(mad, 1e-8)
    index = {gene: i for i, gene in enumerate(genes)}
    biological = {str(g) for values in biological_signatures.values() for g in values}
    candidates = np.asarray([i for i, gene in enumerate(genes) if gene not in biological], dtype=int)
    if len(candidates) == 0:
        raise TargetBuildError("no non-biological genes are available for matched controls")
    controls: dict[str, tuple[str, ...]] = {}
    details: dict[str, object] = {}
    for target, members in read_gmt(biological_signatures).items():
        present = [index[g] for g in members if g in index]
        if not present:
            details[target] = {"status": "unavailable_no_mapped_genes"}
            continue
        count = len(present)
        if count > len(candidates):
            details[target] = {"status": "unavailable_insufficient_candidates", "n_required": count}
            continue
        target_center = standard[present].mean(axis=0)
        distance = np.sum((standard[candidates] - target_center) ** 2, axis=1)
        # Draw from a close candidate pool, then rank deterministically.  This
        # avoids the same exact control repeated 20 times while preserving a
        # transparent, seed-controlled sampling rule.
        pool_size = min(len(candidates), max(count * 8, count))
        pool = candidates[np.argsort(distance, kind="stable")[:pool_size]]
        for control_index in range(controls_per_target):
            rng = np.random.default_rng(seed + int.from_bytes(sha256(f"{target}:{control_index}".encode()).digest()[:4], "little"))
            selected = rng.choice(pool, size=count, replace=False)
            name = f"RANDOM_CONTROL__{target}__{control_index:02d}"
            controls[name] = tuple(genes[i] for i in selected)
        details[target] = {"status": "eligible", "n_genes": count, "candidate_pool": int(pool_size),
                           "controls": [f"RANDOM_CONTROL__{target}__{i:02d}" for i in range(controls_per_target)]}
    digest = sha256(json.dumps(controls, sort_keys=True).encode("utf-8")).hexdigest()
    return controls, {"format_version": 1, "controls_per_target": controls_per_target, "seed": seed,
                      "n_source_rna_patients": n_source_rna_patients,
                      "n_canonical_rna_patients": int(len(rna)),
                      "n_noncanonical_rna_patients_excluded": int(n_source_rna_patients - len(rna)),
                      "fit_population": "development_train_only", "details": details, "digest": digest}


@dataclass(frozen=True)
class TrainOnlyStandardizer:
    """A serialisable target standardizer fit solely on development training rows.

    Cancer offsets are estimated only for cancers present in the fitting subset.
    An unseen outer-test cancer receives offset zero rather than an offset fit
    from its own test examples.  This is the important cancer-safe behaviour.
    """

    patient_id_col: str
    cancer_col: str | None
    targets: tuple[str, ...]
    global_mean: dict[str, float]
    global_std: dict[str, float]
    cancer_offset: dict[str, dict[str, float]]
    train_counts: dict[str, int]

    def transform(self, scores: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
        _require_unique(scores, self.patient_id_col, "score table")
        _require_unique(metadata, self.patient_id_col, "metadata")
        frame = scores[[self.patient_id_col, *self.targets]].copy()
        metadata_index = metadata.set_index(self.patient_id_col)
        missing = set(frame[self.patient_id_col]) - set(metadata_index.index.astype(str))
        if missing:
            raise TargetBuildError("Metadata is missing score-table patient IDs")
        if self.cancer_col is not None and self.cancer_col not in metadata.columns:
            raise TargetBuildError(f"Metadata has no {self.cancer_col!r} column")
        cancer_values = (
            metadata_index.loc[frame[self.patient_id_col].astype(str), self.cancer_col]
            .astype(str)
            .tolist()
            if self.cancer_col is not None
            else [None] * len(frame)
        )
        for target in self.targets:
            raw = pd.to_numeric(frame[target], errors="coerce")
            offsets = np.array(
                [self.cancer_offset[target].get(cancer, 0.0) for cancer in cancer_values],
                dtype=float,
            )
            transformed = (raw.to_numpy(dtype=float) - self.global_mean[target] - offsets)
            frame[target] = transformed / self.global_std[target]
        return frame

    def to_dict(self) -> dict[str, object]:
        return {
            "patient_id_col": self.patient_id_col,
            "cancer_col": self.cancer_col,
            "targets": list(self.targets),
            "global_mean": self.global_mean,
            "global_std": self.global_std,
            "cancer_offset": self.cancer_offset,
            "train_counts": self.train_counts,
        }


def fit_train_only_standardizer(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    patient_id_col: str = "patient_id",
    cancer_col: str | None = "cancer",
    train_mask: Sequence[bool] | pd.Series | None = None,
    split_col: str = "split",
    train_label: str = "train",
    allow_unavailable_targets: bool = False,
) -> TrainOnlyStandardizer:
    """Fit global scale and optional cancer offsets using only training patients."""
    _require_unique(scores, patient_id_col, "score table")
    _require_unique(metadata, patient_id_col, "metadata")
    targets = tuple(column for column in scores.columns if column != patient_id_col)
    if not targets:
        raise TargetBuildError("Score table contains no targets")
    indexed_metadata = metadata.copy()
    indexed_metadata[patient_id_col] = indexed_metadata[patient_id_col].astype(str)
    score_ids = scores[patient_id_col].astype(str)
    aligned_meta = indexed_metadata.set_index(patient_id_col).reindex(score_ids)
    if aligned_meta.index.isna().any() or aligned_meta.shape[0] != len(scores):
        raise TargetBuildError("Metadata is missing score-table patient IDs")
    if cancer_col is not None and cancer_col not in aligned_meta.columns:
        raise TargetBuildError(f"Metadata has no {cancer_col!r} column")

    if train_mask is None:
        if split_col not in aligned_meta.columns:
            raise TargetBuildError("Provide train_mask or metadata containing a split column")
        mask = aligned_meta[split_col].astype(str).str.lower().eq(train_label.lower()).to_numpy()
    else:
        mask = np.asarray(train_mask, dtype=bool)
        if mask.shape != (len(scores),):
            raise TargetBuildError("train_mask must have exactly one value per score-table row")
    if not mask.any():
        raise TargetBuildError("No training rows were supplied to fit the standardizer")

    global_mean: dict[str, float] = {}
    global_std: dict[str, float] = {}
    cancer_offset: dict[str, dict[str, float]] = {}
    train_counts: dict[str, int] = {}
    cancers = aligned_meta[cancer_col].astype(str).to_numpy() if cancer_col is not None else None
    for target in targets:
        values = pd.to_numeric(scores[target], errors="coerce").to_numpy(dtype=float)
        valid_train = mask & np.isfinite(values)
        count = int(valid_train.sum())
        if count < 2 and not allow_unavailable_targets:
            raise TargetBuildError(
                f"Target {target!r} has fewer than two available training values"
            )
        if count < 2:
            # A per-target unavailable signature must not delete every other
            # valid target.  Preserve it as unavailable and record neutral
            # transform parameters only in bundle-building mode.
            global_mean[target] = 0.0
            global_std[target] = 1.0
            train_counts[target] = count
            cancer_offset[target] = {}
            continue
        mean = float(np.mean(values[valid_train]))
        std = float(np.std(values[valid_train], ddof=0))
        # A constant target should remain finite and neutral after transform.
        if not np.isfinite(std) or std < 1e-12:
            std = 1.0
        global_mean[target] = mean
        global_std[target] = std
        train_counts[target] = count
        offsets: dict[str, float] = {}
        if cancers is not None:
            for cancer in np.unique(cancers[valid_train]):
                cancer_mask = valid_train & (cancers == cancer)
                offsets[str(cancer)] = float(np.mean(values[cancer_mask]) - mean)
        cancer_offset[target] = offsets
    return TrainOnlyStandardizer(
        patient_id_col=patient_id_col,
        cancer_col=cancer_col,
        targets=targets,
        global_mean=global_mean,
        global_std=global_std,
        cancer_offset=cancer_offset,
        train_counts=train_counts,
    )


@dataclass
class TargetBundle:
    """Canonical patient-aligned RNA targets and their provenance."""

    raw_scores: pd.DataFrame
    standardized_scores: pd.DataFrame
    availability: pd.DataFrame
    standardizer: TrainOnlyStandardizer
    manifest: dict[str, object]

    def write(self, output_dir: str | Path) -> None:
        """Write interoperable tabular outputs and a JSON manifest."""
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        self.raw_scores.to_parquet(output / "raw_rna_targets.parquet", index=False)
        self.standardized_scores.to_parquet(output / "rna_targets.parquet", index=False)
        self.availability.to_parquet(output / "rna_target_availability.parquet", index=False)
        with (output / "rna_target_manifest.json").open("w", encoding="utf-8") as handle:
            json.dump(self.manifest, handle, indent=2, sort_keys=True)


def build_rna_target_bundle(
    rna_source: str | Path | pd.DataFrame,
    gmt_source: str | Path | Mapping[str, Sequence[str]],
    metadata: pd.DataFrame,
    *,
    canonical_patient_ids: Sequence[str] | None = None,
    patient_id_col: str = "patient_id",
    cancer_col: str | None = "cancer",
    train_mask: Sequence[bool] | pd.Series | None = None,
    split_col: str = "split",
    train_label: str = "train",
    min_present_genes: int = 1,
    score_mode: str = "competitive_rank",
) -> TargetBundle:
    """Build canonical, masked, train-only-standardised RNA target tables."""
    rna = read_rna_table(rna_source, patient_id_col=patient_id_col)
    signatures = read_gmt(gmt_source)
    metadata = metadata.copy()
    _require_unique(metadata, patient_id_col, "metadata")
    metadata[patient_id_col] = metadata[patient_id_col].astype(str)

    raw_scores, availability, matched_genes = score_gene_signatures(
        rna, signatures, patient_id_col=patient_id_col, min_present_genes=min_present_genes,
        score_mode=score_mode,
    )
    if canonical_patient_ids is None:
        canonical_ids = metadata[patient_id_col].tolist()
    else:
        canonical_ids = [str(patient_id) for patient_id in canonical_patient_ids]
        if len(set(canonical_ids)) != len(canonical_ids):
            raise TargetBuildError("canonical_patient_ids contains duplicates")
    canonical = pd.DataFrame({patient_id_col: canonical_ids})
    raw_scores = canonical.merge(raw_scores, on=patient_id_col, how="left", validate="one_to_one")
    availability = canonical.merge(availability, on=patient_id_col, how="left", validate="one_to_one")
    target_columns = [name for name in signatures]
    # Assign column-wise to avoid pandas object downcasting and keep the
    # availability contract explicitly boolean for every target.
    for target in target_columns:
        availability[target] = availability[target].map(lambda value: bool(value) if pd.notna(value) else False).astype(bool)
    metadata_aligned = canonical.merge(
        metadata, on=patient_id_col, how="left", validate="one_to_one"
    )
    if metadata_aligned.drop(columns=[patient_id_col]).isna().all(axis=1).any():
        raise TargetBuildError("Metadata is missing one or more canonical patient IDs")

    standardizer = fit_train_only_standardizer(
        raw_scores,
        metadata_aligned,
        patient_id_col=patient_id_col,
        cancer_col=cancer_col,
        train_mask=train_mask,
        split_col=split_col,
        train_label=train_label,
        allow_unavailable_targets=True,
    )
    standardized = standardizer.transform(raw_scores, metadata_aligned)
    signature_digest = sha256(
        json.dumps(signatures, sort_keys=True).encode("utf-8")
    ).hexdigest()
    manifest: dict[str, object] = {
        "format_version": 1,
        "patient_id_col": patient_id_col,
        "target_columns": target_columns,
        "n_canonical_patients": len(canonical),
        "n_rna_patients": len(rna),
        "n_patients_with_rna": int(canonical[patient_id_col].isin(rna[patient_id_col]).sum()),
        "raw_scores_digest": _frame_digest(raw_scores),
        "standardized_scores_digest": _frame_digest(standardized),
        "availability_digest": _frame_digest(availability),
        "metadata_digest": _frame_digest(metadata_aligned),
        "signature_digest": signature_digest,
        "score_mode": score_mode,
        "matched_genes": matched_genes,
        "available_patients_by_target": {
            target: int(availability[target].sum()) for target in target_columns
        },
        "standardizer": standardizer.to_dict(),
    }
    return TargetBundle(
        raw_scores=raw_scores,
        standardized_scores=standardized,
        availability=availability,
        standardizer=standardizer,
        manifest=manifest,
    )
