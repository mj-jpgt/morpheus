"""E0/E0b: a fail-closed perturbation-basis transfer experiment.

E0 asks a deliberately narrow question: after matching genes and removing the
trivial leading axis, does a Replogle perturbation-response subspace align with
TCGA bulk RNA *because of perturbation biology*?

The matched-spectrum random-rotation floor alone cannot answer that.  A Haar
random subspace of gene space sits at ``k / n_genes`` (~0.005-0.05), while any
two real gene-expression matrices score 0.4-0.8 from generic co-expression
(ribosomal, cell-cycle, MHC modules) that has nothing to do with perturbation.
The principal-angle statistic is a function of the gene covariance alone, so no
covariance-preserving random null can separate the two.  The control therefore
has to be *biological*: Replogle's ``energy_test_p_value`` marks which
perturbations produced any detectable transcriptional effect at all, so

    responsive arm      = energy_test_p_value < 0.01
    NON-responsive arm  = energy_test_p_value > 0.5   (the negative control)

both run through the identical pipeline on an identical number of rows.  A
non-responsive perturbation still carries the generic covariance and none of
the perturbation biology, so whatever alignment *it* reaches is the real floor.

E0 supports perturbation transfer only when the responsive arm exceeds the
non-responsive arm with non-overlapping bootstrap uncertainty -- never merely
by beating the Haar floor.  This is *not* a causal transfer result.  It is the
gate for deciding whether the perturbation-basis branch is worth building.
"""
from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
import time
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from .gates import GateLedger

HOUSEKEEPING = {"ACTB", "GAPDH", "TUBB", "RPL13A", "B2M"}
# Replogle's documented `gene_transcript` index has stable first three fields
# (row, target, guide).  The fourth field is a transcript/source identifier
# whose spelling is not uniformly Ensembl-only, so validate its presence rather
# than incorrectly discarding otherwise well-formed perturbations.
TARGET_RE = re.compile(r"^(?P<row>\d+)_(?P<target>[^_]+)_(?P<guide>[^_]+)_(?P<transcript>.+)$")
REQUIRED_GATE_PREFIXES = {
    "G0.1", "G0.2", "G0.3", "G0.4", "G1.1", "G1.2", "G1.3", "G1.4", "G1.5", "G1.6",
    "G3.1", "G3.2", "G3.3", "G3.4", "G3.5", "G3.6", "G4.1", "G4.2", "G4.3", "G4.4",
    "G4.5", "G4.6", "E0.a", "E0.b", "E0.c", "E0.d", "E0.e", "E0.f", "E0.g", "E0.h", "E0.i",
}
# Distinct per-transform offsets so the Haar draws, the bootstrap resamples and
# the split-half assignment are independent across transforms.  With a
# transform-independent seed, "replication across transforms" would be
# correlated by construction.
TRANSFORM_SEED_OFFSET = {"signed_log1p": 0, "clip_log1p": 500_000}
# Measured float32 QR orthonormality defect at production shapes (8248x101) is
# ~5e-7; a rotation that is actually broken is O(0.1-1).  1e-4 leaves ~100x
# headroom against float32 noise while still being a real, failable check.
ROTATION_TOLERANCE = 1e-4


@dataclass(frozen=True)
class TransferConfig:
    """Every knob the observed statistic, its null and its bootstrap share.

    ``q`` is the randomized-SVD rank budget.  It must exceed
    ``max(offsets) + max(ks)`` by a real oversampling margin: randomized SVD
    attenuates the trailing retained directions, and that attenuation applies to
    the observed statistic and the ceiling but *not* to the exact Haar null, so
    a thin margin biases the comparison against the observed value.
    """

    ks: tuple[int, ...] = (10, 25, 50, 100)
    offsets: tuple[int, ...] = (0, 1, 2, 5)
    primary_offset: int = 1
    q: int = 150
    min_q: int = 101
    draws: int = 100
    bootstrap_draws: int = 200
    responsive_p: float = 0.01
    nonresponsive_p: float = 0.5

    @property
    def min_rows(self) -> int:
        """Rows an arm needs before ``_right_svd`` can honour the full budget."""
        return self.q + 1


def _symbol(value: object) -> str:
    """Canonical symbol from either SYMBOL or SYMBOL|ENTREZ input columns."""
    return str(value).strip().upper().split("|")[0]


def _git(root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()
    except Exception:
        return "unavailable"


def _meaningful_dirty(root: Path) -> list[str]:
    """Ignore only protocol-required append-only ledgers, never source changes."""
    allowed = {"v2/research/rebase/nature/GATE_LOG.md", "v2/research/rebase/nature/EXPERIMENT_LOG.md"}
    lines = [line for line in _git(root, "status", "--porcelain").splitlines() if line]
    dirty: list[str] = []
    for line in lines:
        path = line[3:].replace("\\", "/") if len(line) > 3 else line
        if path not in allowed: dirty.append(line)
    return dirty


def _right_svd(x: torch.Tensor, q: int, seed: int, *, min_q: int = 101) -> tuple[torch.Tensor, torch.Tensor]:
    """Top gene-space right singular vectors and singular values on GPU/CPU.

    ``min_q`` is the production floor (101 = offset 1 + k 100).  Unit tests lower
    it so the identical estimator can run on small synthetic matrices; nothing on
    the real-data path may pass it.
    """
    q = min(q, x.shape[0] - 1, x.shape[1] - 1)
    if q < min_q:
        raise ValueError(f"need >={min_q} rank candidates for E0; got q={q}")
    # pca_lowrank is deterministic after this explicit device generator seed.
    torch.manual_seed(seed)
    _, s, v = torch.pca_lowrank(x, q=q, center=False, niter=7)
    return v, s


def _overlap(a: torch.Tensor, b: torch.Tensor, k: int, offset: int = 0) -> float:
    singular = torch.linalg.svdvals(a[:, offset:offset + k].T @ b[:, offset:offset + k])
    return float(singular.square().mean().cpu())


def _matched_spectrum_rotated_null(
    source_singular: torch.Tensor, n_genes: int, target: torch.Tensor, *, k: int, draws: int, seed: int, offset: int = 1,
    spectrum_check_draws: int = 3,
) -> tuple[np.ndarray, dict[str, float]]:
    """Randomize P's retained right singular vectors while exactly keeping its spectrum.

    The rank-q matrix ``U diag(source_singular) Q.T`` has the same non-zero
    singular values as the rank-q approximation of P for every Haar-orthogonal
    ``Q``.  Principal-angle overlap depends only on Q, so constructing Q rather
    than materialising 100 dense 11k x 8k matrices is exactly equivalent and
    avoids an otherwise needless 40+ GB allocation.

    IMPORTANT INTERPRETATION LIMIT: this is a *Haar* floor.  It sits at
    ``k / n_genes`` and is therefore only a check that the arithmetic and the
    gene universe are sane.  It is not evidence of perturbation transfer, and
    the E0 verdict does not use it -- see ``_decision``.

    The returned invariant is measured, not restated: the singular values of the
    materialised ``Q diag(s)`` and the orthonormality defect of ``Q`` are both
    computed for the first ``spectrum_check_draws`` draws, so a non-orthonormal
    or mis-scaled rotation makes the reported error grow.
    """
    q = offset + k
    if source_singular.numel() < q:
        raise ValueError("source SVD does not contain enough components")
    generator = torch.Generator(device=target.device).manual_seed(seed)
    values: list[float] = []
    s = source_singular[:q].double()
    reference_energy = float(s.square().sum().cpu())
    orthonormality_error, spectrum_relative_error, checked = 0.0, 0.0, 0
    null_energy = float("nan")
    for draw in range(draws):
        q_basis, _ = torch.linalg.qr(torch.randn((n_genes, q), device=target.device, generator=generator), mode="reduced")
        if draw < spectrum_check_draws:
            wide = q_basis.double()
            gram = wide.T @ wide
            identity = torch.eye(q, dtype=gram.dtype, device=gram.device)
            orthonormality_error = max(orthonormality_error, float((gram - identity).abs().max().cpu()))
            # Singular values of the *materialised* null factor, measured.  The
            # previous implementation compared `s` against itself, so its
            # "spectrum error" was abs(x - x) == 0 for any Q whatsoever.
            null_singular = torch.linalg.svdvals(wide * s)
            spectrum_relative_error = max(spectrum_relative_error, float(((null_singular - s).abs().max() / s.max()).cpu()))
            null_energy = float(null_singular.square().sum().cpu())
            checked += 1
        values.append(_overlap(q_basis, target, k, offset))
    return np.asarray(values, dtype=np.float64), {
        "retained_rank": float(q),
        "source_retained_energy": reference_energy,
        "null_retained_energy": null_energy,
        "spectrum_max_relative_singular_error": spectrum_relative_error,
        "rotation_orthonormality_max_abs_error": orthonormality_error,
        "spectrum_checked_draws": float(checked),
    }


def _gene_label_shuffle_null(source: torch.Tensor, target: torch.Tensor, *, k: int, draws: int, seed: int, offset: int = 1) -> np.ndarray:
    """Preserve the perturbation basis but randomize its gene correspondence."""
    generator = torch.Generator(device=source.device).manual_seed(seed)
    return np.asarray([
        _overlap(source[torch.randperm(source.shape[0], device=source.device, generator=generator)], target, k, offset)
        for _ in range(draws)
    ], dtype=np.float64)


def _parse_targets(index: pd.Index, controls: np.ndarray) -> tuple[list[str | None] | None, dict[str, object]]:
    """Read the documented gene_transcript index schema, or mark E0b unavailable.

    These AnnData files do not carry a separate target column.  We only derive a
    target from the named ``gene_transcript`` index after *every* non-control
    row validates against its documented four-field schema; otherwise guide
    retrieval is explicitly unavailable rather than guessed from a split.
    """
    if str(index.name) != "gene_transcript":
        return None, {"guide_target_status": "unavailable_index_schema", "index_name": str(index.name)}
    parsed = [TARGET_RE.match(str(value)) for value in index[~controls]]
    valid_fraction = float(np.mean([item is not None for item in parsed]))
    if valid_fraction < .95:
        return None, {"guide_target_status": "unavailable_invalid_index_schema", "invalid_rows": int(sum(x is None for x in parsed)), "valid_fraction": valid_fraction}
    targets: list[str | None] = [item.group("target").upper() if item else None for item in parsed]
    counts = pd.Series([item for item in targets if item is not None]).value_counts()
    replicated = int((counts >= 2).sum())
    if replicated == 0:
        return None, {"guide_target_status": "unavailable_no_replicated_targets", "n_targets": int(len(counts))}
    return targets, {"guide_target_status": "validated_gene_transcript_index", "n_targets": int(len(counts)), "n_replicated_targets": replicated, "index_schema_valid_fraction": valid_fraction}


def _aggregate_duplicate_symbols(x: np.ndarray, genes: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """Mean technical duplicate columns once, preventing duplicated-gene weighting."""
    positions: dict[str, list[int]] = {}
    for index, gene in enumerate(genes.tolist()): positions.setdefault(str(gene), []).append(index)
    duplicates = int(sum(len(indices) - 1 for indices in positions.values()))
    if not duplicates: return x, genes, 0
    ordered = list(positions)
    merged = np.empty((x.shape[0], len(ordered)), dtype=np.float32)
    for index, gene in enumerate(ordered):
        merged[:, index] = x[:, positions[gene]].mean(axis=1)
    return merged, np.asarray(ordered), duplicates


@dataclass
class MatrixBundle:
    x: np.ndarray
    genes: list[str]
    row_ids: list[str]
    meta: dict[str, object]
    targets: list[str | None] | None = None
    groups: list[str] | None = None
    # Replogle's per-perturbation energy test.  Row-aligned with ``x``; None when
    # the source file does not carry it (then the control arm is unavailable and
    # E0 may not certify transfer).
    energy_p: np.ndarray | None = None


def _load_perturbation(path: Path) -> MatrixBundle:
    data = ad.read_h5ad(path, backed="r")
    raw = np.asarray(data.X[:], dtype=np.float32)
    controls = data.obs["core_control"].fillna(False).astype(bool).to_numpy()
    if not controls.any():
        raise ValueError(f"{path}: no core_control rows")
    genes = np.asarray([_symbol(v) for v in data.var["gene_name"].to_numpy()])
    finite_rows = np.isfinite(raw).all(axis=1)
    all_zero_rows = np.all(np.nan_to_num(raw, nan=0.0) == 0.0, axis=1)
    keep_rows = finite_rows & ~all_zero_rows
    # Evaluate column finiteness on the ROWS WE KEEP. Computing it on the unfiltered
    # matrix drops a gene for a NaN that occurs only in a row keep_rows discards
    # anyway: measured 73 genes in K562 and 2 in RPE1, ALL of them finite in the
    # retained rows. n_genes_dropped_nonfinite now means what it says.
    finite_cols = np.isfinite(raw[keep_rows]).all(axis=0)
    if not keep_rows[controls].any():
        raise ValueError(f"{path}: no usable control rows after integrity filtering")
    # The control rows test delta semantics before any optional centring.  The
    # documented control-expression fields must also agree with the measured
    # target-gene direction in X; a near-zero centroid by itself is not enough.
    raw = raw[:, finite_cols]; genes = genes[finite_cols]
    raw, genes, duplicates = _aggregate_duplicate_symbols(raw, genes)
    control_centroid = raw[controls & keep_rows].mean(axis=0)
    response = raw[(~controls) & keep_rows]
    response_norm = np.linalg.norm(response, axis=1)
    centroid_ratio = float(np.linalg.norm(control_centroid) / max(np.median(response_norm), 1e-12))
    fold = data.obs["fold_expr"].to_numpy(dtype=float)
    pct = data.obs["pct_expr"].to_numpy(dtype=float)
    control_expr = data.obs["control_expr"].to_numpy(dtype=float)
    valid_fold = np.isfinite(fold) & np.isfinite(pct)
    fold_relation_error = float(np.max(np.abs(pct[valid_fold] - (fold[valid_fold] - 1.0)))) if valid_fold.any() else float("inf")
    targets, target_meta = _parse_targets(data.obs_names, controls)
    target_vector_corr = float("nan"); target_vector_mae = float("inf"); target_vector_coverage = 0.0; target_gene_hits = 0
    if targets is not None:
        symbols = {str(gene): i for i, gene in enumerate(genes)}
        raw_noncontrol = np.flatnonzero(~controls)
        target_cols = np.asarray([symbols.get(target, -1) if target is not None else -1 for target in targets])
        target_gene_hits = int((target_cols >= 0).sum())
        valid_target = (target_cols >= 0) & np.isfinite(fold[raw_noncontrol]) & (fold[raw_noncontrol] >= 0)
        target_vector_coverage = float(valid_target.mean())
        if valid_target.sum() >= 30:
            rows = raw_noncontrol[valid_target]; columns = target_cols[valid_target]
            observed_target_effect = raw[rows, columns]
            # `control_expr` is the measured control baseline and `fold_expr`
            # is perturbed/control expression.  In this normalized-bulk file,
            # the matching X entry must therefore equal the log expression
            # change implied by those two stored quantities.
            expected_target_effect = np.log2(control_expr[rows] * fold[rows] + 1.0) - np.log2(control_expr[rows] + 1.0)
            finite_effect = np.isfinite(observed_target_effect) & np.isfinite(expected_target_effect)
            observed_target_effect, expected_target_effect = observed_target_effect[finite_effect], expected_target_effect[finite_effect]
            target_vector_coverage = float(len(observed_target_effect) / max(1, target_gene_hits))
            if len(observed_target_effect) >= 30:
                target_vector_corr = float(np.corrcoef(observed_target_effect, expected_target_effect)[0, 1])
                target_vector_mae = float(np.mean(np.abs(observed_target_effect - expected_target_effect)))
    control_expr_coverage = float(np.isfinite(control_expr[~controls]).mean())
    control_expr_nonnegative = bool(np.all(control_expr[~controls][np.isfinite(control_expr[~controls])] >= 0))
    min_target_gene_hits = max(30, int(np.ceil(.10 * max(1, (~controls).sum()))))
    matrix_is_delta = (centroid_ratio <= 0.10 and fold_relation_error <= 1e-5 and control_expr_coverage >= .80
                       and control_expr_nonnegative and target_gene_hits >= min_target_gene_hits and target_vector_coverage >= .70 and target_vector_corr >= .80 and target_vector_mae <= .30)
    if not matrix_is_delta:
        raise ValueError(f"{path}: cannot assert delta response (centroid_ratio={centroid_ratio:.3g}, control_expr_coverage={control_expr_coverage:.3g}, target_gene_hits={target_gene_hits}, target_coverage={target_vector_coverage:.3g}, target_control_ref_r={target_vector_corr:.3g}, target_control_ref_mae={target_vector_mae:.3g}, target_examples={(targets or [])[:3]}, gene_examples={genes[:3].tolist()})")
    usable_controls = controls & keep_rows
    noncontrol = (~controls) & keep_rows
    # targets are indexed over raw non-controls, so retain only usable
    # non-controls in the same order.  No invalid non-control rows are silently
    # retained in this dataset, but retain this explicit alignment safeguard.
    if targets is not None and np.any(~keep_rows[~controls]):
        # Keep target labels in exactly the same non-control row order as X.
        # This is an explicit, auditable index operation—not an assumption that
        # filtering never occurred.
        targets = [target for target, keep in zip(targets, keep_rows[~controls]) if keep]
    x = raw[noncontrol]
    # The energy test is the only per-perturbation statement in the file about
    # whether a perturbation did anything at all.  It is what makes a biological
    # negative control arm possible, so it is loaded row-aligned with X or
    # explicitly declared unavailable.
    energy_p: np.ndarray | None = None
    energy_meta: dict[str, object] = {"energy_test_status": "unavailable_missing_energy_test_p_value"}
    if "energy_test_p_value" in data.obs.columns:
        energy_p = data.obs["energy_test_p_value"].to_numpy(dtype=float)[noncontrol]
        finite_energy = np.isfinite(energy_p)
        energy_meta = {"energy_test_status": "available", "energy_p_finite_fraction": float(finite_energy.mean()),
                       "n_energy_responsive_p_lt_0p01": int((finite_energy & (energy_p < .01)).sum()),
                       "n_energy_nonresponsive_p_gt_0p5": int((finite_energy & (energy_p > .5)).sum())}
    return MatrixBundle(x=x, genes=genes.tolist(), row_ids=data.obs_names[noncontrol].astype(str).tolist(), targets=targets, energy_p=energy_p,
        meta={"source": str(path), "shape_raw": list(data.shape), "n_control_rows_raw": int(controls.sum()),
              "n_control_rows_used": int(usable_controls.sum()), "n_rows_dropped_nonfinite": int((~finite_rows).sum()),
              "n_rows_dropped_all_zero": int((all_zero_rows & finite_rows).sum()), "n_genes_dropped_nonfinite": int((~finite_cols).sum()),
              "n_duplicate_symbols_aggregated": duplicates, "n_perturbation_rows": int(noncontrol.sum()), "control_centroid_ratio": centroid_ratio,
              "fold_pct_relation_max_abs_error": fold_relation_error, "control_expr_coverage": control_expr_coverage,
              "control_expr_nonnegative": control_expr_nonnegative, "target_vector_coverage": target_vector_coverage,
              "target_gene_hits": target_gene_hits, "min_target_gene_hits": min_target_gene_hits,
              "target_vector_control_referenced_correlation": target_vector_corr, "target_vector_control_referenced_mae": target_vector_mae,
              "delta_status": "validated_control_centred" if matrix_is_delta else "FAILED", **target_meta, **energy_meta})


def _load_tcga(path: Path, transform: str) -> MatrixBundle:
    frame = pd.read_parquet(path)
    if "patient_id" not in frame:
        raise ValueError("TCGA parquet is missing patient_id")
    genes = np.asarray([_symbol(c) for c in frame.columns if c != "patient_id"])
    if len(genes) != len(set(genes)):
        raise ValueError("TCGA has duplicate canonical gene symbols")
    raw = frame.drop(columns="patient_id").to_numpy(dtype=np.float32, copy=False)
    finite_rows = np.isfinite(raw).all(axis=1)
    all_zero_rows = np.all(np.nan_to_num(raw, nan=0.0) == 0.0, axis=1)
    keep_rows = finite_rows & ~all_zero_rows
    # See the note in _load_perturbation: column finiteness must be judged on the
    # retained rows, not the raw matrix.
    finite_cols = np.isfinite(raw[keep_rows]).all(axis=0)
    raw = raw[keep_rows][:, finite_cols]
    raw_min, raw_max = float(raw.min()), float(raw.max())
    if transform == "signed_log1p":
        x = np.sign(raw) * np.log2(np.abs(raw) + 1.0)
    elif transform == "clip_log1p":
        x = np.log2(np.maximum(raw, 0.0) + 1.0)
    else:
        raise ValueError(transform)
    zero_fraction = float((x == 0).mean())
    flat = x.reshape(-1)
    mean, std = float(flat.mean()), float(flat.std())
    skew = float(np.mean(((flat - mean) / max(std, 1e-12)) ** 3))
    return MatrixBundle(x=x, genes=genes[finite_cols].tolist(), row_ids=frame.loc[keep_rows, "patient_id"].astype(str).tolist(), meta={
        "source": str(path), "n_patients_raw": int(len(frame)), "n_patients_used": int(keep_rows.sum()), "n_rows_dropped_nonfinite": int((~finite_rows).sum()),
        "n_rows_dropped_all_zero": int((all_zero_rows & finite_rows).sum()), "n_genes_dropped_nonfinite": int((~finite_cols).sum()),
        "raw_min": raw_min, "raw_max": raw_max, "transform": transform, "post_transform_min": float(x.min()), "post_transform_max": float(x.max()),
        "post_transform_zero_fraction": zero_fraction, "post_transform_skew": skew})


def _restrict_tcga_to_registry(bundle: MatrixBundle, registry_path: Path) -> MatrixBundle:
    """Use the canonical patient/cancer registry for a composition-checked ceiling."""
    registry = pd.read_parquet(registry_path)
    cancer_column = "cancer_type" if "cancer_type" in registry else "cancer"
    if "patient_id" not in registry or cancer_column not in registry:
        raise ValueError(f"{registry_path}: requires patient_id and cancer_type/cancer")
    registry = registry[["patient_id", cancer_column]].dropna().drop_duplicates("patient_id")
    lookup = dict(zip(registry.patient_id.astype(str), registry[cancer_column].astype(str)))
    groups = [lookup.get(patient) for patient in bundle.row_ids]
    keep = np.asarray([g is not None for g in groups])
    matrix_coverage = float(keep.mean())
    registry_coverage = float(len(set(np.asarray(bundle.row_ids)[keep])) / max(1, len(lookup)))
    # The registry is the declared evaluation universe; it is expected to be a
    # subset of PanCan RNA.  The join gate must therefore measure how much of
    # *that* universe was recovered, not penalise the intentionally excluded
    # PanCan-only patients in the denominator.
    if registry_coverage < .80:
        raise ValueError(f"canonical registry coverage {registry_coverage:.3f}<0.80; refuse unstratified ceiling")
    return MatrixBundle(x=bundle.x[keep], genes=bundle.genes, row_ids=np.asarray(bundle.row_ids)[keep].tolist(),
        groups=np.asarray(groups, dtype=object)[keep].astype(str).tolist(), meta={**bundle.meta, "registry_source": str(registry_path),
        "registry_patients": int(len(registry)), "registry_joined": int(keep.sum()), "registry_coverage": registry_coverage,
        "full_matrix_registry_overlap": matrix_coverage, "n_patients_used": int(keep.sum())})


def _align(p: MatrixBundle, t: MatrixBundle) -> tuple[np.ndarray, np.ndarray, list[str], dict[str, int]]:
    target_index = {gene: i for i, gene in enumerate(t.genes)}
    matches = [(i, target_index[g], g) for i, g in enumerate(p.genes) if g in target_index]
    if len(matches) < 100:
        raise ValueError(f"only {len(matches)} shared genes")
    pi, ti, genes = zip(*matches)
    px, tx = p.x[:, pi], t.x[:, ti]
    ps, ts = px.std(axis=0), tx.std(axis=0)
    constant = (ps < 1e-8) | (ts < 1e-8) | ~np.isfinite(ps) | ~np.isfinite(ts)
    if constant.any():
        px, tx = px[:, ~constant], tx[:, ~constant]
        genes = tuple(np.asarray(genes)[~constant])
    px = (px - px.mean(axis=0, keepdims=True)) / px.std(axis=0, keepdims=True)
    tx = (tx - tx.mean(axis=0, keepdims=True)) / tx.std(axis=0, keepdims=True)
    return px.astype(np.float32), tx.astype(np.float32), list(genes), {"n_left": len(p.genes), "n_right": len(t.genes), "n_matched_pre_constant": len(matches), "n_constant_columns_dropped": int(constant.sum()), "n_matched_evaluated": len(genes)}


def _effective_rank_torch(x: torch.Tensor) -> float:
    """Roy-Vetterli effective rank computed on the active device."""
    centred = x - x.mean(dim=0, keepdim=True)
    singular = torch.linalg.svdvals(centred)
    singular = singular[singular > 1e-12]
    if singular.numel() == 0: return 0.0
    weights = singular / singular.sum()
    return float(torch.exp(-(weights * torch.log(weights)).sum()).cpu())


def _health(x: torch.Tensor) -> dict[str, float]:
    var = x.var(dim=0, unbiased=False); mean_var = var.mean().clamp_min(1e-12)
    dead = float((var < 1e-6 * mean_var).float().mean().cpu())
    centred = x - x.mean(0, keepdim=True)
    normalized = centred / centred.std(0, keepdim=True).clamp_min(1e-12)
    corr = (normalized.T @ normalized) / max(1, x.shape[0])
    mask = ~torch.eye(corr.shape[0], dtype=torch.bool, device=x.device)
    duplicate = float((corr[mask].abs() > .99).float().mean().cpu())
    norms = torch.linalg.vector_norm(x, dim=1)
    return {"effective_rank": _effective_rank_torch(x), "dead_fraction": dead, "duplicate_pair_fraction": duplicate,
            "min_sample_std": float(x.std(0).min().cpu()), "mean_norm": float(norms.mean().cpu()), "median_norm": float(norms.median().cpu()),
            "n_nonfinite": int((~torch.isfinite(x)).sum().cpu())}


def _full_dictionary_rank(p: MatrixBundle, device: torch.device) -> dict[str, object]:
    """Exact full-spectrum E0b rank report for the filtered response dictionary.

    The rank cut is the standard LAPACK/numpy one, ``s[0] * max(n, p) * eps``,
    taken at the *storage* precision of the data (float32 here) while the
    spectrum itself is computed in float64.  The previous ``s[0] * 1e-10`` cut
    was below float32 eps (~1.19e-7), so it never fired and every dictionary was
    reported at algebraic rank ``min(n, p)`` -- exactly the overcomplete-rank
    overclaim E0b exists to correct.
    """
    x = p.x.astype(np.float32, copy=False)
    std = x.std(axis=0); keep = np.isfinite(std) & (std > 1e-8)
    if not keep.all(): x = x[:, keep]
    x = (x - x.mean(axis=0, keepdims=True)) / x.std(axis=0, keepdims=True)
    pt = torch.as_tensor(x, device=device).double()
    singular = torch.linalg.svdvals(pt)
    relative_tolerance = float(max(x.shape) * np.finfo(np.float32).eps)
    nonzero = singular[singular > singular[0] * relative_tolerance]
    if nonzero.numel() == 0: nonzero = singular[:1]
    weights = nonzero / nonzero.sum()
    return {"rank_mode": f"full_{device.type}_svd_float64", "n_rows": int(x.shape[0]), "n_genes": int(x.shape[1]),
            "constant_columns_dropped": int((~keep).sum()), "rank_relative_tolerance": relative_tolerance,
            "rank_tolerance_basis": "max(n,p)*finfo(float32).eps on a float64 spectrum",
            "numerical_rank": int(nonzero.numel()), "max_possible_rank": int(min(x.shape)),
            "effective_rank": float(torch.exp(-(weights * torch.log(weights)).sum()).cpu()),
            "stable_rank": float((singular.square().sum() / singular[0].square()).cpu())}


def _dictionary_metrics(x: torch.Tensor, targets: list[str | None] | None, threshold: float = .95) -> dict[str, object]:
    z = x / torch.linalg.vector_norm(x, dim=1, keepdim=True).clamp_min(1e-12)
    n, block = z.shape[0], 256; maximum = 0.0; edges = 0
    parent = list(range(n))
    def find(i: int) -> int:
        while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
        return i
    def union(i: int, j: int) -> None:
        i, j = find(i), find(j)
        if i != j: parent[j] = i
    hits = eligible = 0
    usable_targets = bool(targets) and all(target is not None for target in targets)
    counts = pd.Series(targets).value_counts().to_dict() if usable_targets else {}
    for start in range(0, n, block):
        sim = z[start:start + block] @ z.T
        for local, row in enumerate(sim):
            i = start + local
            absolute = row.abs()
            # Exclude self from every dictionary quantity.  Masking before
            # `abs()` would turn a negative sentinel back into a high score.
            coherence = absolute.clone(); coherence[i] = 0.0
            maximum = max(maximum, float(coherence.max().cpu()))
            for j in torch.where(coherence >= threshold)[0].cpu().tolist(): edges += 1; union(i, int(j))
            nearest = coherence.clone(); nearest[i] = -1.0
            if usable_targets and counts[targets[i]] > 1:
                eligible += 1; hits += int(targets[i] == targets[int(nearest.argmax().cpu())])
    return {"dictionary_coherence_abs": maximum, "equivalence_threshold_abs": threshold, "equivalence_edges_directed": edges,
            "n_equivalence_classes": len({find(i) for i in range(n)}), "guide_same_target_retrieval_at1": hits / eligible if eligible else None,
            "guide_retrieval_eligible": eligible, "guide_retrieval_status": "available" if usable_targets else "unavailable"}


def _pc1_gate(observed: float, null_p95: float) -> bool:
    return bool(np.isfinite(observed) and np.isfinite(null_p95) and observed > null_p95)


def _shuffle_gate(observed: float, shuffled_p95: float) -> bool:
    """True correspondence must beat the gene-label permutation floor."""
    return bool(np.isfinite(observed) and np.isfinite(shuffled_p95) and observed > shuffled_p95)


def _delta_gate(meta: dict[str, object]) -> bool:
    return meta["delta_status"] == "validated_control_centred"


def _split_half(x: torch.Tensor, groups: list[str], q: int, seed: int, *, min_q: int = 101) -> tuple[torch.Tensor, torch.Tensor, dict[str, object]]:
    """Seeded cancer-stratified split, with every sufficiently large cancer in both halves."""
    if len(groups) != x.shape[0]: raise ValueError("missing TCGA cancer groups for split-half ceiling")
    rng = np.random.default_rng(seed); groups_a = np.asarray(groups); left: list[int] = []; right: list[int] = []
    for group in sorted(set(groups)):
        idx = np.flatnonzero(groups_a == group); idx = rng.permutation(idx)
        cut = len(idx) // 2
        if cut == 0: raise ValueError(f"cancer {group} has fewer than two samples")
        left.extend(idx[:cut]); right.extend(idx[cut:])
    a, b = torch.as_tensor(left, device=x.device), torch.as_tensor(right, device=x.device)
    va, _ = _right_svd(x[a], q, seed + 1, min_q=min_q); vb, _ = _right_svd(x[b], q, seed + 2, min_q=min_q)
    all_counts = pd.Series(groups).value_counts(); a_counts = pd.Series(groups_a[np.asarray(left)]).value_counts(); b_counts = pd.Series(groups_a[np.asarray(right)]).value_counts()
    return va, vb, {"n_half_a": int(a.numel()), "n_half_b": int(b.numel()), "n_cancers": int(len(all_counts)),
        "max_group_count_imbalance": int(max(abs(a_counts.reindex(all_counts.index, fill_value=0) - b_counts.reindex(all_counts.index, fill_value=0))))}


ARM_SEED_OFFSET = {"all_perturbations": 0, "responsive_full": 13_000, "responsive_matched": 26_000, "nonresponsive": 39_000}


def _unavailable_arm(arm: str, reason: str, n_rows: int = 0) -> dict[str, object]:
    """Repo convention: an arm we could not measure is visible, never silent."""
    note = reason if reason.startswith("unavailable_") else f"unavailable_{reason}"
    return {"arm": arm, "available": False, "n_rows": int(n_rows), "metric": "status", "value": float("nan"),
            "status": note, "note": note}


def _energy_arms(energy_p: np.ndarray | None, *, cfg: TransferConfig, seed: int) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    """Split perturbations into responsive and NON-responsive (control) arms.

    The responsive arm is subsampled, deterministically, to exactly the
    non-responsive arm's row count.  Subspace estimation quality is a function
    of n, so an unmatched responsive-vs-control comparison measures sample size
    as much as biology and is not usable as a decision rule.
    """
    if energy_p is None:
        return {}, {"arm_status": "unavailable_no_energy_test_p_value"}
    energy_p = np.asarray(energy_p, dtype=float)
    finite = np.isfinite(energy_p)
    responsive_full = np.flatnonzero(finite & (energy_p < cfg.responsive_p))
    nonresponsive = np.flatnonzero(finite & (energy_p > cfg.nonresponsive_p))
    rng = np.random.default_rng(seed)
    if len(nonresponsive) and len(responsive_full) > len(nonresponsive):
        matched = np.sort(rng.choice(responsive_full, size=len(nonresponsive), replace=False))
    else:
        matched = responsive_full.copy()
    meta = {"arm_status": "available", "n_perturbation_rows": int(len(energy_p)),
            "n_energy_p_nonfinite": int((~finite).sum()), "responsive_p_threshold": cfg.responsive_p,
            "nonresponsive_p_threshold": cfg.nonresponsive_p, "arm_subsample_seed": int(seed),
            "n_responsive_full": int(len(responsive_full)), "n_responsive_matched": int(len(matched)),
            "n_nonresponsive": int(len(nonresponsive)),
            "arms_are_n_matched": bool(len(matched) == len(nonresponsive) and len(nonresponsive) > 0)}
    return {"responsive_full": responsive_full, "responsive_matched": matched, "nonresponsive": nonresponsive}, meta


def _rezscore(x: torch.Tensor) -> torch.Tensor:
    """Per-column standardisation used identically by every arm and bootstrap draw."""
    return (x - x.mean(0, keepdim=True)) / x.std(0, keepdim=True).clamp_min(1e-8)


def _bootstrap_overlap(x: torch.Tensor, vt: torch.Tensor, *, cfg: TransferConfig, seed: int) -> dict[int, np.ndarray]:
    """Percentile bootstrap over PERTURBATION ROWS of the observed statistic.

    This is the only genuine uncertainty statement E0 makes about the observed
    alignment.  The null-derived interval is a statement about the Haar floor,
    mirrored around the observation; using it as an interval for the
    observation is the same test as E0.b written twice.
    """
    generator = torch.Generator(device=x.device).manual_seed(seed)
    n = int(x.shape[0])
    values: dict[int, list[float]] = {k: [] for k in cfg.ks}
    for draw in range(cfg.bootstrap_draws):
        rows = torch.randint(0, n, (n,), device=x.device, generator=generator)
        vd, _ = _right_svd(_rezscore(x[rows]), cfg.q, seed + 1 + draw, min_q=cfg.min_q)
        for k in cfg.ks:
            values[k].append(_overlap(vd, vt, k, cfg.primary_offset))
    return {k: np.asarray(v, dtype=np.float64) for k, v in values.items()}


def _arm_result(x: torch.Tensor, vt: torch.Tensor, *, ceiling: dict[int, float], n_genes: int, cfg: TransferConfig,
                draws: int, seed: int, arm: str) -> tuple[dict[str, object], torch.Tensor | None, torch.Tensor | None]:
    """Run one arm through the pipeline every other arm runs through.

    Same alignment, same standardisation, same offset, same k grid, same q, same
    Haar null, same bootstrap.  The only difference between arms is which rows
    of P they contain -- which is the entire point.
    """
    n_rows = int(x.shape[0])
    if n_rows < cfg.min_rows:
        return _unavailable_arm(arm, f"unavailable_insufficient_rows_{n_rows}_lt_{cfg.min_rows}", n_rows), None, None
    vp, sp = _right_svd(x, cfg.q, seed, min_q=cfg.min_q)
    bootstrap = _bootstrap_overlap(x, vt, cfg=cfg, seed=seed + 5_000)
    result: dict[str, object] = {"arm": arm, "available": True, "status": "available", "note": "", "n_rows": n_rows,
                                 "n_zero_variance_columns": int((x.std(0) <= 1e-8).sum().cpu()),
                                 "retained_q": int(sp.numel()), "bootstrap_draws": int(cfg.bootstrap_draws)}
    for k in cfg.ks:
        overlaps = {int(offset): _overlap(vp, vt, k, offset) for offset in cfg.offsets if int(sp.numel()) >= offset + k}
        base, stripped = overlaps.get(0), overlaps.get(cfg.primary_offset)
        # G5/E0.b is named for exactly this question and the run has never
        # answered it: how much of the alignment is the trivial leading axis?
        pc1_share = float(1.0 - stripped / base) if (base is not None and base > 1e-12 and stripped is not None) else float("nan")
        observed = stripped if stripped is not None else float("nan")
        null, invariant = _matched_spectrum_rotated_null(sp, n_genes, vt, k=k, draws=draws, seed=seed + 100 + k, offset=cfg.primary_offset)
        null_median = float(np.median(null))
        ceiling_k = float(ceiling.get(k, float("nan")))
        gap = ceiling_k - null_median
        normalised = float((observed - null_median) / gap) if np.isfinite(gap) and gap > 1e-6 else float("nan")
        effect = observed - null
        draw_values = bootstrap[k]
        result[f"k{k}"] = {
            "pc1_removed_overlap": observed, "overlap_by_offset": {str(o): v for o, v in sorted(overlaps.items())},
            "pc1_share": pc1_share, "pc1_removed_ceiling": ceiling_k,
            "normalised_alignment": normalised,
            "normalised_alignment_status": "available" if np.isfinite(normalised) else "unavailable_degenerate_ceiling_null_gap",
            "ceiling_minus_null_median": float(gap),
            "null_median": null_median, "null_p95": float(np.quantile(null, .95)), "null_std": float(np.std(null)),
            "null_p": float((1 + np.sum(null >= observed)) / (draws + 1)),
            "effect_vs_null": float(np.mean(effect)),
            # Explicitly the dispersion of the NULL mirrored around the
            # observation.  Kept for continuity, never used as the uncertainty
            # of the observed statistic.
            "null_effect_ci95": [float(np.quantile(effect, .025)), float(np.quantile(effect, .975))],
            "bootstrap_median": float(np.median(draw_values)),
            # NOT a confidence interval for the observation, and must never be
            # described as one. Resampling rows with replacement leaves ~63% unique
            # rows, so every resampled subspace is estimated from an effectively
            # smaller sample and its overlap regresses downward; measured bias is
            # -0.04 to -0.09 and the observed value can lie ABOVE this interval
            # entirely. Worse for a kill decision, the bias GROWS with an arm's
            # overlap, so it penalises the responsive arm more than the control --
            # i.e. it is biased toward a false NEGATIVE. Use bootstrap_diff_ci95
            # on the paired difference (below), where the common-mode bias cancels.
            "bootstrap_ci95": [float(np.quantile(draw_values, .025)), float(np.quantile(draw_values, .975))],
            "bootstrap_ci95_is_biased_low": True,
            "bootstrap_samples": [float(v) for v in draw_values],
            "bootstrap_std": float(np.std(draw_values)), "bootstrap_draws": int(draw_values.size),
            "beats_haar_floor": bool(np.isfinite(observed) and observed > float(np.quantile(null, .95))),
            "theoretical_random_overlap": k / n_genes, **invariant}
    return result, vp, sp


def _decision(responsive: dict[str, object], control: dict[str, object], k: int) -> dict[str, object]:
    """THE E0 decision rule.

    Transfer is supported at k only when the n-matched responsive arm's
    bootstrap interval lies entirely above the non-responsive control arm's.
    Beating the Haar floor is explicitly *not* sufficient: two unrelated real
    expression matrices clear that floor from generic co-expression alone.
    """
    key = f"k{k}"
    if not responsive.get("available") or not control.get("available"):
        blocked = responsive if not responsive.get("available") else control
        reason = str(blocked.get("status", "unavailable_unknown"))
        return {"responsive_exceeds_nonresponsive": False, "decision_status": reason, "metric": "status",
                "value": float("nan"), "note": reason, "blocking_arm": blocked.get("arm"),
                "responsive_available": bool(responsive.get("available")), "control_available": bool(control.get("available"))}
    r, c = responsive[key], control[key]
    # n-matching is part of the rule, not a report: subspace estimation quality
    # is a function of n, so an unmatched win is not a win.
    n_matched = bool(responsive["n_rows"] == control["n_rows"])
    # PAIRED difference bootstrap. Draw i of each arm resamples the same number of
    # rows under the same generator schedule, so the downward resampling bias is
    # common-mode and largely cancels in resp_i - ctrl_i. Requiring two MARGINAL
    # intervals not to overlap is both badly biased (see bootstrap_ci95) and far
    # less powerful -- it is roughly an alpha of 0.005 rather than 0.05.
    rs, cs = np.asarray(r.get("bootstrap_samples", []), float), np.asarray(c.get("bootstrap_samples", []), float)
    paired = rs[:min(len(rs), len(cs))] - cs[:min(len(rs), len(cs))]
    if paired.size:
        diff_ci = [float(np.quantile(paired, .025)), float(np.quantile(paired, .975))]
        diff_median = float(np.median(paired))
    else:
        diff_ci, diff_median = [float("nan")] * 2, float("nan")
    exceeds = bool(n_matched and np.isfinite(diff_ci[0]) and diff_ci[0] > 0.0)
    # The stricter non-overlap criterion is retained for comparison, never as the rule.
    marginal_exceeds = bool(n_matched and np.isfinite(r["bootstrap_ci95"][0]) and np.isfinite(c["bootstrap_ci95"][1])
                            and r["bootstrap_ci95"][0] > c["bootstrap_ci95"][1])
    return {"responsive_exceeds_nonresponsive": exceeds, "decision_status": "available",
            "responsive_overlap": r["pc1_removed_overlap"], "nonresponsive_overlap": c["pc1_removed_overlap"],
            "responsive_ci95": r["bootstrap_ci95"], "nonresponsive_ci95": c["bootstrap_ci95"],
            "bootstrap_diff_ci95": diff_ci, "bootstrap_diff_median": diff_median,
            "n_paired_bootstrap_draws": int(paired.size),
            "responsive_exceeds_nonresponsive_marginal_ci": marginal_exceeds,
            "overlap_gap": float(r["pc1_removed_overlap"] - c["pc1_removed_overlap"]),
            "responsive_normalised_alignment": r["normalised_alignment"], "nonresponsive_normalised_alignment": c["normalised_alignment"],
            "n_rows": [int(responsive["n_rows"]), int(control["n_rows"])], "arms_are_n_matched": n_matched,
            # Recorded so a reader can see the difference between the old,
            # insufficient criterion and the one E0 now decides on.
            "responsive_beats_haar_floor_only": bool(r["beats_haar_floor"]),
            "control_beats_haar_floor": bool(c["beats_haar_floor"])}


def _context(name: str, p: MatrixBundle, t: MatrixBundle, full_dictionary_rank: dict[str, object], *, device: torch.device,
             draws: int, seed: int, cfg: TransferConfig | None = None) -> dict[str, object]:
    cfg = cfg or TransferConfig()
    px, tx, genes, join = _align(p, t)
    pt, tt = torch.as_tensor(px, device=device), torch.as_tensor(tx, device=device)
    vt, _ = _right_svd(tt, cfg.q, seed + 1, min_q=cfg.min_q)
    va, vb, split_meta = _split_half(tt, t.groups or [], cfg.q, seed + 2, min_q=cfg.min_q)
    ceiling = {k: _overlap(va, vb, k, cfg.primary_offset) for k in cfg.ks}
    ceiling_by_offset = {k: {str(o): _overlap(va, vb, k, o) for o in cfg.offsets if cfg.q >= o + k} for k in cfg.ks}
    all_arm, vp, sp = _arm_result(pt, vt, ceiling=ceiling, n_genes=len(genes), cfg=cfg, draws=draws,
                                  seed=seed + ARM_SEED_OFFSET["all_perturbations"], arm="all_perturbations")
    if vp is None or sp is None:
        raise ValueError(f"{name}: perturbation matrix has too few usable rows ({pt.shape[0]}) for E0")
    arm_rows, arm_meta = _energy_arms(p.energy_p, cfg=cfg, seed=seed + 11)
    arms: dict[str, object] = {"all_perturbations": all_arm}
    for arm_name in ("responsive_full", "responsive_matched", "nonresponsive"):
        rows = arm_rows.get(arm_name)
        if rows is None or len(rows) == 0:
            reason = str(arm_meta.get("arm_status", "unavailable_no_energy_test_p_value"))
            arms[arm_name] = _unavailable_arm(arm_name, reason if reason.startswith("unavailable_") else "unavailable_empty_arm",
                                              0 if rows is None else len(rows))
            continue
        arm_x = _rezscore(pt[torch.as_tensor(np.asarray(rows), device=device)])
        arms[arm_name], _, _ = _arm_result(arm_x, vt, ceiling=ceiling, n_genes=len(genes), cfg=cfg, draws=draws,
                                           seed=seed + ARM_SEED_OFFSET[arm_name], arm=arm_name)
    p_latent, t_latent = pt @ vp, tt @ vt
    p_health, t_health = _health(p_latent), _health(t_latent)
    # E0 uses randomized rank-q PCA for the principal-angle test.  A full
    # NumPy SVD here would silently move the run back to CPU for minutes/hours.
    # Report rank *within that explicitly retained q subspace* instead; this is
    # the relevant representation health quantity, not a claim about the
    # uncomputed full matrix rank.
    out: dict[str, object] = {"context": name, "n_shared_genes": len(genes), "join": join, "perturbation": p.meta, "tcga": t.meta,
        "config": {"ks": list(cfg.ks), "offsets": list(cfg.offsets), "primary_offset": cfg.primary_offset, "q": cfg.q,
                   "oversampling": cfg.q - (cfg.primary_offset + max(cfg.ks)), "draws": draws, "bootstrap_draws": cfg.bootstrap_draws},
        "housekeepers_present": sorted(HOUSEKEEPING.intersection(genes)), "perturbation_health": p_health, "tcga_health": t_health,
        "rank": {"retained_components": int(sp.numel()), "effective_rank_retained": p_health["effective_rank"],
                 "stable_rank_retained_estimate": float((sp.square().sum() / sp[0].square()).cpu()), "dictionary_full": full_dictionary_rank},
        "dictionary": _dictionary_metrics(pt, p.targets), "split": split_meta, "arm_split": arm_meta, "arms": arms,
        "shapes": {"p_rows": int(pt.shape[0]), "p_columns": int(pt.shape[1]), "t_rows": int(tt.shape[0]), "t_columns": int(tt.shape[1]),
                   "declared_n_perturbations": int(p.meta.get("n_perturbation_rows", -1)), "declared_n_patients": int(t.meta.get("n_patients_used", -1))}}
    for k in cfg.ks:
        base = dict(all_arm[f"k{k}"])
        observed = base["pc1_removed_overlap"]
        # The positive ceiling is evaluated against an independently held-out
        # TCGA half and its own matched-spectrum randomized-P floor.  Thus the
        # positive and negative controls use the exact same principal-angle
        # routine, PC stripping, and gene universe.
        ceiling_null, _ = _matched_spectrum_rotated_null(sp, len(genes), vb, k=k, draws=draws, seed=seed + 300 + k, offset=cfg.primary_offset)
        shuffled = _gene_label_shuffle_null(vp, vt, k=k, draws=draws, seed=seed + 700 + k, offset=cfg.primary_offset)
        heldout_a, heldout_b = _overlap(vp, va, k, cfg.primary_offset), _overlap(vp, vb, k, cfg.primary_offset)
        matched_decision = _decision(arms["responsive_matched"], arms["nonresponsive"], k)
        full_decision = _decision(arms["responsive_full"], arms["nonresponsive"], k)
        base.update({"ceiling_by_offset": ceiling_by_offset[k], "heldout_half_overlap_mean": (heldout_a + heldout_b) / 2,
            "gene_label_shuffle_median": float(np.median(shuffled)), "gene_label_shuffle_p95": float(np.quantile(shuffled, .95)),
            "gene_label_shuffle_std": float(np.std(shuffled)), "gene_label_shuffle_p": float((1 + np.sum(shuffled >= observed)) / (draws + 1)),
            "ceiling_null_median": float(np.median(ceiling_null)), "ceiling_null_p95": float(np.quantile(ceiling_null, .95)),
            "decision_n_matched": matched_decision, "decision_responsive_full": full_decision,
            "responsive_exceeds_nonresponsive": bool(matched_decision["responsive_exceeds_nonresponsive"])})
        out[f"k{k}"] = base
    return out


def _self_test() -> None:
    """Static/synthetic preflight: broken controls must fail before GPU spending."""
    device = torch.device("cpu"); generator = torch.Generator().manual_seed(7)
    g, q, k = 128, 12, 5
    v, _ = torch.linalg.qr(torch.randn(g, q, generator=generator), mode="reduced")
    s = torch.arange(q, 0, -1, dtype=torch.float32)
    null, invariant = _matched_spectrum_rotated_null(s, g, v, k=k, draws=100, seed=3, offset=1)
    # Measured, not restated: a non-orthonormal rotation would grow both errors.
    assert invariant["spectrum_max_relative_singular_error"] < ROTATION_TOLERANCE and invariant["rotation_orthonormality_max_abs_error"] < ROTATION_TOLERANCE
    assert invariant["spectrum_checked_draws"] > 0 and null.std() > 0.0
    assert abs(invariant["null_retained_energy"] - invariant["source_retained_energy"]) < ROTATION_TOLERANCE * invariant["source_retained_energy"]
    assert _gene_label_shuffle_null(v, v, k=k, draws=100, seed=4).std() > 0.0
    # The Haar floor must sit at k/n_genes: it is a sanity check on the gene
    # universe, never evidence of transfer.
    assert abs(float(np.median(null)) - k / g) < 0.5 * k / g
    try: _matched_spectrum_rotated_null(s[:3], g, v, k=k, draws=100, seed=3, offset=1)
    except ValueError: pass
    else: raise AssertionError("truncated spectrum did not fail")
    try: _parse_targets(pd.Index(["bad"], name="gene_transcript"), np.asarray([False]))
    except Exception: raise
    assert _parse_targets(pd.Index(["bad"], name="gene_transcript"), np.asarray([False]))[0] is None
    # Deliberately corrupted controls must fail the same predicates used by the
    # live ledger, not merely a nearby toy assertion.
    assert _pc1_gate(.2, .1) and not _pc1_gate(.1, .1)
    assert _shuffle_gate(.2, .1) and not _shuffle_gate(.05, .1)
    assert _delta_gate({"delta_status": "validated_control_centred"}) and not _delta_gate({"delta_status": "FAILED"})
    assert 1.0 <= _effective_rank_torch(torch.eye(8)) <= 8.0
    toy_rank = _full_dictionary_rank(MatrixBundle(np.asarray([[1., 0.], [0., 1.], [1., 1.]], dtype=np.float32), ["A", "B"], ["x", "y", "z"], {}), device)
    assert toy_rank["rank_mode"] == "full_cpu_svd_float64" and toy_rank["numerical_rank"] == 2
    # A deficient-rank dictionary must report its true rank, not min(n, p).  The
    # previous 1e-10 relative cut sat below float32 eps and never fired.
    deficient_generator = np.random.default_rng(11)
    deficient = (deficient_generator.normal(size=(120, 8)) @ deficient_generator.normal(size=(8, 60))).astype(np.float32)
    deficient_rank = _full_dictionary_rank(MatrixBundle(deficient, [f"G{i}" for i in range(60)], [f"r{i}" for i in range(120)], {}), device)
    assert deficient_rank["numerical_rank"] <= 9, deficient_rank["numerical_rank"]
    # Arms: an under-sized arm is explicitly unavailable and never certifies.
    tiny = _unavailable_arm("nonresponsive", "unavailable_insufficient_rows_80_lt_151", 80)
    assert tiny["metric"] == "status" and np.isnan(tiny["value"]) and tiny["note"].startswith("unavailable_")
    assert _decision({"available": True, "n_rows": 5, "k5": {}}, tiny, 5)["responsive_exceeds_nonresponsive"] is False
    arm_config = TransferConfig(ks=(5,), q=8, min_q=4, draws=10, bootstrap_draws=10)
    _, arm_meta = _energy_arms(np.asarray([*np.zeros(50), *np.ones(7)]), cfg=arm_config, seed=1)
    assert arm_meta["n_responsive_matched"] == arm_meta["n_nonresponsive"] == 7 and arm_meta["arms_are_n_matched"]
    # Same-target labels deliberately disagree with nearest non-self vectors;
    # this catches accidental self-neighbour retrieval immediately.
    self_leak = _dictionary_metrics(torch.tensor([[1., 0.], [0., 1.], [1., 0.], [0., 1.]]), ["A", "A", "B", "B"])
    assert self_leak["dictionary_coherence_abs"] == 1.0 and self_leak["guide_same_target_retrieval_at1"] == 0.0
    # Exercise the real H5AD loader with deliberately inconsistent X versus
    # (control_expr, fold_expr).  A cosmetic metadata check would pass this;
    # the E0.a data-semantic guard must not.
    with tempfile.TemporaryDirectory() as temporary:
        obs = pd.DataFrame({"core_control": [False, False, False, True], "control_expr": [2., 2., 2., np.nan],
                            "fold_expr": [.1, .1, .1, np.nan], "pct_expr": [-.9, -.9, -.9, np.nan]},
                           index=pd.Index(["0_A_P1_ENSG000001", "1_A_P2_ENSG000001", "2_A_P3_ENSG000001", "3_non-targeting_non-targeting_non-targeting"], name="gene_transcript"))
        toy = ad.AnnData(X=np.asarray([[0., 1.], [0., 1.], [0., 1.], [0., .01]], dtype=np.float32), obs=obs, var=pd.DataFrame({"gene_name": ["A", "B"]}))
        path = Path(temporary) / "inconsistent.h5ad"; toy.write_h5ad(path)
        try: _load_perturbation(path)
        except ValueError as error: assert "cannot assert delta response" in str(error)
        else: raise AssertionError("control_expr/X inconsistency did not fail E0.a")
        # Technical duplicate target columns are collapsed to their mean once,
        # and the target/control semantic check is still applied to that one
        # canonical gene rather than indexing the pre-aggregation layout.
        folds = np.linspace(.2, .8, 30); expected = np.log2(2.0 * folds + 1.0) - np.log2(3.0)
        valid_obs = pd.DataFrame({"core_control": [*np.zeros(len(folds), dtype=bool), True], "control_expr": [*np.full(len(folds), 2.), np.nan],
                                  "fold_expr": [*folds, np.nan], "pct_expr": [*(folds - 1.0), np.nan]},
                                 index=pd.Index([*(f"{i}_A_P{i}_ENSG000001" for i in range(len(folds))), f"{len(folds)}_non-targeting_non-targeting_non-targeting"], name="gene_transcript"))
        valid = ad.AnnData(X=np.vstack([np.column_stack([expected - .2, expected + .2, np.ones(len(folds))]), [0., 0., .01]]).astype(np.float32), obs=valid_obs, var=pd.DataFrame({"gene_name": ["A", "A", "B"]}))
        valid_path = Path(temporary) / "duplicate_valid.h5ad"; valid.write_h5ad(valid_path)
        duplicate_bundle = _load_perturbation(valid_path)
        assert duplicate_bundle.genes == ["A", "B"] and duplicate_bundle.meta["n_duplicate_symbols_aggregated"] == 1
        assert np.allclose(duplicate_bundle.x[:, 0], expected) and duplicate_bundle.meta["delta_status"] == "validated_control_centred"
    # A gate manifest with only one row must fail coverage: this protects
    # against accidentally declaring a run healthy because a new gate was
    # simply never wired into the ledger.
    assert _missing_gate_prefixes(["G0.1_workspace_real"]) != []
    print("E0 static/synthetic preflight PASS")


def _missing_gate_prefixes(gates: list[str]) -> list[str]:
    return sorted(prefix for prefix in REQUIRED_GATE_PREFIXES if not any(gate.startswith(prefix) for gate in gates))


def _add_gates(ledger: GateLedger, context: dict[str, object], label: str, draws: int, ks: tuple[int, ...] = (10, 25, 50, 100)) -> None:
    p, t, join = context["perturbation"], context["tcga"], context["join"]
    shapes, arms, arm_split = context["shapes"], context["arms"], context["arm_split"]
    # The stated threshold and the tested condition must be the same quantity.
    # This gate previously declared "0 dropped" and tested the *matched* count,
    # so 5000 dropped columns still passed.
    constant_fraction = join["n_constant_columns_dropped"] / max(1, join["n_matched_pre_constant"])
    ledger.add("G1.1_no_constant_columns", json.dumps({"dropped":join["n_constant_columns_dropped"],"pre_constant":join["n_matched_pre_constant"],"dropped_fraction":constant_fraction,"evaluated":join["n_matched_evaluated"]}),
               "n_constant_columns_dropped <= 0.01*n_matched_pre_constant AND n_matched_evaluated > 100",
               constant_fraction <= .01 and join["n_matched_evaluated"] > 100, label)
    ledger.add("G1.2_real_gene_join", f"{join['n_left']},{join['n_right']},{join['n_matched_pre_constant']}", ">=0.8*min(left,right)", join["n_matched_pre_constant"] >= .8 * min(join["n_left"], join["n_right"]), label)
    # Finding, not validity. "Observed does not beat the gene-shuffle null" is
    # ambiguous between a broken join and a true null, and the audit showed the
    # shuffle null is numerically indistinguishable from the Haar null anyway
    # (p95 0.00606 vs 0.00611). Join validity is carried by G1.2 (real gene join)
    # and G1.6 (housekeeping genes survive), which test the join itself.
    ledger.observe("G1.3_join_not_scrambled", "per-k observed vs gene-label permutation",
                   "recorded; join validity is enforced by G1.2/G1.6, not by this statistic", label)
    scale = {k: t[k] for k in ("raw_max", "post_transform_max", "post_transform_zero_fraction", "post_transform_skew")}
    ledger.add("G1.4_scale_sanity", json.dumps(scale), "all finite; post_transform_max>0; 0<=zero_fraction<1",
               all(np.isfinite(v) for v in scale.values()) and scale["post_transform_max"] > 0 and 0 <= scale["post_transform_zero_fraction"] < 1, label)
    drops = {"P_nonfinite":p["n_rows_dropped_nonfinite"],"T_nonfinite":t["n_rows_dropped_nonfinite"],"P_zero":p["n_rows_dropped_all_zero"],"T_zero":t["n_rows_dropped_all_zero"]}
    ledger.add("G1.5_explicit_row_filtering", json.dumps(drops), "every drop counted and P/T drops leave the declared row counts intact",
               all(v >= 0 for v in drops.values()) and shapes["p_rows"] == shapes["declared_n_perturbations"] and shapes["t_rows"] == shapes["declared_n_patients"], label)
    ledger.add("G1.6_housekeeping_mapping", ",".join(context["housekeepers_present"]), "all five present", set(context["housekeepers_present"]) == HOUSEKEEPING, label)
    ledger.add("E0.a_delta_not_absolute", json.dumps({"status":p["delta_status"],"ratio":p["control_centroid_ratio"],"fold_pct_error":p["fold_pct_relation_max_abs_error"],"control_expr_coverage":p["control_expr_coverage"],"target_control_ref_r":p["target_vector_control_referenced_correlation"],"target_control_ref_mae":p["target_vector_control_referenced_mae"]}), "X target effect agrees with control_expr x fold_expr reference", _delta_gate(p), label)
    # Assert the measured shapes of the matrices actually handed to the SVD, not
    # the literal True.  A transpose error reads as null, so this must be a test.
    orientation_ok = (shapes["p_rows"] == shapes["declared_n_perturbations"] and shapes["t_rows"] == shapes["declared_n_patients"]
                      and shapes["p_columns"] == shapes["t_columns"] == context["n_shared_genes"])
    ledger.add("E0.c_orientation", json.dumps(shapes), "P rows==n_perturbations, T rows==n_patients, both column counts==n_shared_genes", orientation_ok, label)
    ledger.add("G1.2b_split_composition", json.dumps(context["split"]), "cancer-stratified; imbalance<=1", context["split"]["max_group_count_imbalance"] <= 1, label)
    for who in ("perturbation_health", "tcga_health"):
        h = context[who]; prefix = "P" if who.startswith("perturbation") else "T"
        ledger.add(f"G3.1_effective_rank_{prefix}", h["effective_rank"], "finite", np.isfinite(h["effective_rank"]), label)
        ledger.add(f"G3.2_dead_dimensions_{prefix}", h["dead_fraction"], "report", h["dead_fraction"] <= .5, label)
        ledger.add(f"G3.3_duplicate_dimensions_{prefix}", h["duplicate_pair_fraction"], "report", np.isfinite(h["duplicate_pair_fraction"]), label)
        ledger.add(f"G3.4_sample_variation_{prefix}", h["min_sample_std"], ">0", h["min_sample_std"] > 0, label)
        # RECORDED OBSERVATION, NOT A TEST: neither input carries site labels, so
        # site degeneracy cannot be measured here.  Declared unavailable rather
        # than asserted.
        ledger.add(f"G3.5_site_degeneracy_{prefix}", "unavailable_no_site_labels", "RECORDED OBSERVATION, NOT A TEST: no site labels in either input", True,
                   f"{label}; recorded observation, not a test")
        ledger.add(f"G3.6_norm_sanity_{prefix}", json.dumps({"mean":h["mean_norm"],"median":h["median_norm"],"nonfinite":h["n_nonfinite"]}), "zero nonfinite", h["n_nonfinite"] == 0, label)
    d = context["dictionary"]
    full_rank = context["rank"]["dictionary_full"]
    ledger.add("E0b_full_dictionary_rank", json.dumps({k:full_rank[k] for k in ("rank_mode","numerical_rank","max_possible_rank","rank_relative_tolerance","effective_rank","stable_rank")}),
               "float64 spectrum with a max(n,p)*eps cut; finite ranks; rank <= min(n,p)",
               str(full_rank["rank_mode"]).endswith("_float64") and full_rank["numerical_rank"] <= full_rank["max_possible_rank"]
               and all(np.isfinite(full_rank[k]) for k in ("numerical_rank","effective_rank","stable_rank")), label)
    ledger.add("E0b_dictionary_coherence", d["dictionary_coherence_abs"], "finite", np.isfinite(d["dictionary_coherence_abs"]), label)
    ledger.add("E0b_equivalence_classes", d["n_equivalence_classes"], ">0", d["n_equivalence_classes"] > 0, label)
    ledger.add("E0b_guide_retrieval", d["guide_retrieval_status"], "available or explicit unavailable", d["guide_retrieval_status"] in {"available","unavailable"}, label)
    # NOT a gate. RPE1 has only 80 non-responsive perturbations, far below the rank
    # this estimator needs, so its control arm is legitimately unavailable. Failing
    # the whole run for that would quarantine a perfectly valid K562 answer. The
    # verdict already refuses "supported" for any context whose control arm is
    # missing, which is where that constraint belongs.
    ledger.observe("E0.h_control_arm_present", json.dumps(arm_split),
                   "responsive and non-responsive arms built from energy_test_p_value and n-matched", label)
    for k in ks:
        r = context[f"k{k}"]; tag = f"k{k}_{label}"; decision = r["decision_n_matched"]
        ledger.add(f"G4.1_positive_control_split_half_{tag}", r["pc1_removed_ceiling"], "> same-path heldout matched-null p95", r["pc1_removed_ceiling"] > r["ceiling_null_p95"], label)
        ledger.add(f"G4.2_negative_control_rotated_P_{tag}", r["null_median"], "near k/n_genes", abs(r["null_median"] - r["theoretical_random_overlap"]) < 0.5 * r["theoretical_random_overlap"], label)
        ledger.observe(f"G4.2b_negative_control_gene_shuffle_{tag}", json.dumps({"median":r["gene_label_shuffle_median"],"p95":r["gene_label_shuffle_p95"],"p":r["gene_label_shuffle_p"]}), "recorded; near-identical to the Haar null by construction, so not independent corroboration", label)
        ledger.add(f"G4.3_null_sanity_{tag}", json.dumps({"std":r["null_std"],"p95":r["null_p95"]}), "std>0; nondegenerate", r["null_std"] > 0, label)
        ledger.add(f"G4.4_heldout_vs_insample_{tag}", json.dumps({"in":r["pc1_removed_overlap"],"heldout":r["heldout_half_overlap_mean"]}), "reported; no fitted mapping", np.isfinite(r["heldout_half_overlap_mean"]), label)
        ledger.add(f"G4.5_null_resolution_{tag}", 1/(draws+1), "draws>=100", draws >= 100, label)
        # Uncertainty on the OBSERVED statistic is the row bootstrap, not the
        # mirrored null dispersion (which restates E0.b).
        ledger.observe(f"G4.6_effect_ci_{tag}", json.dumps({"bootstrap_ci95":r["bootstrap_ci95"],"null_p95":r["null_p95"]}),
                       "finding, not validity: an extreme negative here is a RESULT and must not fail the run", label)
        # SPLIT. "Is the trivial axis dominating?" is validity -> gate.
        # "Does the observed beat the Haar floor?" is a finding -> observation.
        # Fused, an extreme negative would fail the run and quarantine its output.
        ledger.add(f"E0.b_pc1_not_dominant_{tag}", r["pc1_share"],
                   "pc1_share < 0.9 so the result is not merely the mean-expression axis",
                   bool(np.isfinite(r["pc1_share"]) and r["pc1_share"] < 0.9), label)
        ledger.observe(f"E0.b_observed_vs_haar_floor_{tag}",
                       json.dumps({"observed": r["pc1_removed_overlap"], "null_p95": r["null_p95"]}),
                       "recorded for comparison; beating the Haar floor is explicitly NOT the decision rule", label)
        ledger.add(f"E0.d_ceiling_{tag}", r["pc1_removed_ceiling"], ">null p95", r["pc1_removed_ceiling"] > r["null_p95"], label)
        ledger.add(f"E0.e_matched_spectrum_floor_{tag}", json.dumps({"spectrum_rel_error":r["spectrum_max_relative_singular_error"],"orthonormality_error":r["rotation_orthonormality_max_abs_error"],"checked_draws":r["spectrum_checked_draws"]}),
                   f"measured rotation: relative singular error<{ROTATION_TOLERANCE} and ||QtQ-I||<{ROTATION_TOLERANCE} on >=1 draw",
                   r["spectrum_checked_draws"] >= 1 and r["spectrum_max_relative_singular_error"] < ROTATION_TOLERANCE
                   and r["rotation_orthonormality_max_abs_error"] < ROTATION_TOLERANCE, label)
        # THE DECISION -- and deliberately NOT a gate. "Responsive does not exceed
        # non-responsive" is E0 returning a true negative, which is a result the
        # handoff explicitly plans for. Registering it as a gate would mark the run
        # FAIL and quarantine its output into FAILED_*.json exactly when the science
        # came back "no" -- making a real null indistinguishable from a broken
        # pipeline, the one confusion this experiment exists to prevent.
        ledger.observe(f"E0.g_responsive_exceeds_nonresponsive_{tag}", json.dumps(decision),
                       "n-matched responsive bootstrap CI lower > non-responsive CI upper", label)
        # Validity: is the ceiling-null gap wide enough for the normalisation to mean
        # anything? That IS a gate. The normalised VALUE it produces is an observation.
        ledger.add(f"E0.i_normalised_alignment_interpretable_{tag}", r["ceiling_minus_null_median"],
                   "ceiling-null gap > 1e-6 so (observed-null)/(ceiling-null) is interpretable",
                   r["normalised_alignment_status"] == "available", label)
        ledger.observe(f"E0.i_normalised_alignment_value_{tag}", json.dumps({"normalised": r["normalised_alignment"], "status": r["normalised_alignment_status"]}),
                       "reported for scale; no pass/fail threshold is predeclared", label)


def _append_pending_experiment_log(root: Path, output: Path, result: dict[str, object]) -> None:
    """Record a completed gate-valid run without pre-empting the adversarial verdict."""
    log = root / "v2/research/rebase/nature/EXPERIMENT_LOG.md"; log.parent.mkdir(parents=True, exist_ok=True)
    line = (f"\n## {datetime.now(timezone.utc).isoformat()} — E0/E0b pending adversarial audit\n"
            f"- Code: `{result['code_sha']}`; device: `{result['device']}`; gates: `PASS`.\n"
            f"- Output: `{output}`. This entry is a provenance record, not a scientific conclusion until the required adversarial audit is logged.\n")
    with log.open("a", encoding="utf-8") as handle: handle.write(line)


def _verdict(results: dict[str, object], ks: tuple[int, ...]) -> dict[str, object]:
    """Separate 'the run is invalid' from 'the answer is no'.

    A negative E0 is a *result*, not a malfunction, and the handoff explicitly
    plans for it.  The verdict therefore records which contexts could be decided
    at all, and refuses "supported" wherever the biological control arm was not
    available -- that combination is the exact overclaim this experiment exists
    to prevent.
    """
    per_context: dict[str, object] = {}
    for key, context in results["contexts"].items():
        decisions = [context[f"k{k}"]["decision_n_matched"] for k in ks]
        available = [d for d in decisions if d["decision_status"] == "available"]
        if not available:
            per_context[key] = {"verdict": "indeterminate", "reason": decisions[0]["decision_status"],
                                "metric": "status", "value": float("nan"), "note": decisions[0]["decision_status"]}
            continue
        exceeds = [bool(d["responsive_exceeds_nonresponsive"]) for d in available]
        per_context[key] = {"verdict": "supported" if all(exceeds) else "not_supported",
                            "k_decided": len(available), "k_supported": int(sum(exceeds)),
                            "haar_floor_only_would_have_said": "supported" if all(bool(d.get("responsive_beats_haar_floor_only")) for d in available) else "not_supported"}
    # Aggregate over DECIDABLE contexts only. RPE1 has just 50 non-responsive
    # perturbations -- far below the rank this estimator needs -- so its control arm
    # is permanently unavailable. Letting that veto the aggregate made
    # verdict=="supported" literally unreachable: the strongest possible K562 result
    # would have printed as "indeterminate". The undecidable contexts stay visible
    # in contexts_undecidable rather than silently outvoting the evidence.
    decided = {v["verdict"] for v in per_context.values() if v["verdict"] != "indeterminate"}
    undecidable = sorted(k for k, v in per_context.items() if v["verdict"] == "indeterminate")
    overall = ("indeterminate" if not decided
               else "supported" if decided == {"supported"} else "not_supported")
    return {"transfer_supported": overall == "supported", "verdict": overall, "per_context": per_context, "contexts_undecidable": undecidable, "n_contexts_decided": len(per_context) - len(undecidable),
            "rule": "paired row-bootstrap of (responsive - nonresponsive) overlap has 95% CI lower bound > 0, on n-matched arms, at every k"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k562"); parser.add_argument("--rpe1"); parser.add_argument("--tcga"); parser.add_argument("--tcga-registry"); parser.add_argument("--output")
    parser.add_argument("--draws", type=int, default=100); parser.add_argument("--bootstrap-draws", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42); parser.add_argument("--device", default="cuda")
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2])); parser.add_argument("--workspace-link"); parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test: _self_test(); return
    if not all((args.k562,args.rpe1,args.tcga,args.tcga_registry,args.output,args.workspace_link)): parser.error("k562, rpe1, tcga, tcga-registry, output, and workspace-link are required")
    if args.draws < 100: parser.error("--draws must be >=100 for E0's predeclared floor control")
    if args.bootstrap_draws < 200: parser.error("--bootstrap-draws must be >=200 for a usable percentile interval")
    cfg = TransferConfig(draws=args.draws, bootstrap_draws=args.bootstrap_draws)
    started = time.monotonic(); root = Path(args.repo_root).resolve(); output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    device = torch.device(args.device if args.device == "cpu" or torch.cuda.is_available() else "cpu")
    ledger = GateLedger(output, "E0_E0b", official_log=root / "v2/research/rebase/nature/GATE_LOG.md")
    source = Path(__file__).resolve(); workspace_link = Path(args.workspace_link)
    link_is_real = workspace_link.is_symlink() and workspace_link.resolve() == root and root in source.parents
    ledger.add("G0.1_workspace_real", f"link={workspace_link};target={workspace_link.resolve() if workspace_link.exists() else 'missing'}", "declared runtime link resolves to exact git root", link_is_real, "Linux symlink is the remote equivalent of the required Windows Junction; prevents stale copied source")
    dirty = _meaningful_dirty(root); ledger.add("G0.2_code_identity", _git(root,"rev-parse","HEAD"), "clean source/config worktree (protocol ledgers allowed)", not dirty, f"dirty={dirty or 'none'}")
    for label, value in (("G0.3_k562_identity",args.k562),("G0.3_rpe1_identity",args.rpe1),("G0.3_tcga_identity",args.tcga),("G0.3_tcga_registry_identity",args.tcga_registry)): ledger.artifact(label, value)
    # RECORDED OBSERVATION, NOT A TEST: E0 compares no trained artifacts, so
    # there is no manifest to verify.  The input identities that *can* be tested
    # are the G0.3 artifact digests above.
    ledger.add("G0.4_manifest_read", "data-source experiment; no trained artifacts compared", "RECORDED OBSERVATION, NOT A TEST: not applicable, explicitly declared", True,
               "recorded observation, not a test; input identities emitted in input_manifest.json and gated by G0.3")
    p_k, p_r = _load_perturbation(Path(args.k562)), _load_perturbation(Path(args.rpe1))
    # E0b is intentionally expensive and exact: two full GPU spectra, once
    # per cell context, reused across the RNA transform sensitivity analysis.
    full_dictionary_ranks = {"K562": _full_dictionary_rank(p_k, device), "RPE1": _full_dictionary_rank(p_r, device)}
    results: dict[str, object] = {"schema_version":"3.0", "experiment":"E0_E0b", "device":str(device), "draws":args.draws,
        "bootstrap_draws":args.bootstrap_draws, "seed":args.seed, "config":{"ks":list(cfg.ks),"offsets":list(cfg.offsets),"q":cfg.q,"min_q":cfg.min_q,"primary_offset":cfg.primary_offset},
        "command":sys.argv, "code_sha":_git(root,"rev-parse","HEAD"), "code_dirty":dirty, "python":sys.version, "torch":torch.__version__, "cuda":torch.version.cuda, "platform":platform.platform(), "contexts":{}}
    for transform in ("signed_log1p", "clip_log1p"):
        tcga = _restrict_tcga_to_registry(_load_tcga(Path(args.tcga), transform), Path(args.tcga_registry))
        for name, perturb, base_seed in (("K562",p_k,args.seed),("RPE1",p_r,args.seed+1000)):
            # The transform enters the seed.  With a transform-independent seed
            # both transforms drew the SAME Haar matrices and the same bootstrap
            # resamples, so "replication across transforms" was correlated by
            # construction rather than evidence.
            seed = base_seed + TRANSFORM_SEED_OFFSET[transform]
            key = f"{name}_{transform}"; context = _context(name, perturb, tcga, full_dictionary_ranks[name], device=device, draws=args.draws, seed=seed, cfg=cfg)
            context["seed"] = seed
            results["contexts"][key] = context; _add_gates(ledger, context, key, args.draws, cfg.ks)
    # Transform robustness is predeclared, and it is now the biological decision
    # rule that has to replicate -- not the Haar-floor comparison.
    concordant = all(c[f"k{k}"]["responsive_exceeds_nonresponsive"] for c in results["contexts"].values() for k in cfg.ks)
    # Replication of a FINDING, not a validity check -- see E0.g.
    ledger.observe("E0.f_cross_context_transform_replication", json.dumps({"concordant": bool(concordant)}),
                   "responsive arm exceeds non-responsive arm in every decidable context/transform/k")
    missing = _missing_gate_prefixes([str(row["gate"]) for row in ledger.rows])
    ledger.add("G0.5_required_gate_coverage", ",".join(missing) if missing else "complete", "no mandatory gate family missing", not missing)
    results["verdict_detail"] = _verdict(results, cfg.ks)
    results["verdict"] = results["verdict_detail"]["verdict"]
    results["wall_seconds"] = time.monotonic() - started
    results["gates_pass"] = ledger.write()
    (output / "input_manifest.json").write_text(json.dumps({"inputs":{k:results[k] for k in ("command","code_sha","code_dirty","device","draws","bootstrap_draws","seed")},"k562":p_k.meta,"rpe1":p_r.meta},indent=2))
    # Fail closed: a verdict-bearing file that says PASS-shaped things while its
    # gates failed is the worst possible artifact to leave on disk.
    name = "e0_basis_transfer.json" if results["gates_pass"] else "FAILED_e0_basis_transfer.json"
    (output / name).write_text(json.dumps(results, indent=2))
    if results["gates_pass"]: _append_pending_experiment_log(root, output, results)
    print(json.dumps({"output":str(output/name),"gates_pass":results["gates_pass"],"verdict":results["verdict"],
                      "wall_seconds":results["wall_seconds"]},indent=2), flush=True)
    if not results["gates_pass"]: sys.exit(1)


if __name__ == "__main__": main()
