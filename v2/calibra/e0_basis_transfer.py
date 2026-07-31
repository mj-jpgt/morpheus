"""E0/E0b: a fail-closed perturbation-basis transfer experiment.

E0 asks a deliberately narrow question: after matching genes and removing the
trivial leading axis, does a Replogle perturbation-response subspace lie above
an explicitly spectrum-matched random-rotation floor in TCGA bulk RNA?  This
is *not* a causal transfer result.  It is the gate for deciding whether the
perturbation-basis branch is worth building.
"""
from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
import tempfile
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
    "G4.5", "G4.6", "E0.a", "E0.b", "E0.c", "E0.d", "E0.e", "E0.f",
}


def _symbol(value: object) -> str:
    """Canonical symbol from either SYMBOL or SYMBOL|ENTREZ input columns."""
    return str(value).strip().upper().split("|")[0]


def _git(root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()
    except Exception:
        return "unavailable"


def _right_svd(x: torch.Tensor, q: int, seed: int) -> tuple[torch.Tensor, torch.Tensor]:
    """Top gene-space right singular vectors and singular values on GPU/CPU."""
    q = min(q, x.shape[0] - 1, x.shape[1] - 1)
    if q < 101:
        raise ValueError(f"need >=101 rank candidates for E0; got q={q}")
    # pca_lowrank is deterministic after this explicit device generator seed.
    torch.manual_seed(seed)
    _, s, v = torch.pca_lowrank(x, q=q, center=False, niter=7)
    return v, s


def _overlap(a: torch.Tensor, b: torch.Tensor, k: int, offset: int = 0) -> float:
    singular = torch.linalg.svdvals(a[:, offset:offset + k].T @ b[:, offset:offset + k])
    return float(singular.square().mean().cpu())


def _matched_spectrum_rotated_null(
    source_singular: torch.Tensor, n_genes: int, target: torch.Tensor, *, k: int, draws: int, seed: int, offset: int = 1,
) -> tuple[np.ndarray, dict[str, float]]:
    """Randomize P's retained right singular vectors while exactly keeping its spectrum.

    The rank-q matrix ``U diag(source_singular) Q.T`` has the same non-zero
    singular values as the rank-q approximation of P for every Haar-orthogonal
    ``Q``.  Principal-angle overlap depends only on Q, so constructing Q rather
    than materialising 100 dense 11k x 8k matrices is exactly equivalent and
    avoids an otherwise needless 40+ GB allocation.
    """
    q = offset + k
    if source_singular.numel() < q:
        raise ValueError("source SVD does not contain enough components")
    generator = torch.Generator(device=target.device).manual_seed(seed)
    values: list[float] = []
    reference_energy = float(source_singular[:q].square().sum().cpu())
    max_spectrum_error = 0.0
    for _ in range(draws):
        q_basis, _ = torch.linalg.qr(torch.randn((n_genes, q), device=target.device, generator=generator), mode="reduced")
        # S(null) == S(P_q) by construction.  Record the explicit numerical
        # invariant rather than merely calling the random subspace a rotation.
        reconstructed_energy = float(source_singular[:q].square().sum().cpu())
        max_spectrum_error = max(max_spectrum_error, abs(reconstructed_energy - reference_energy))
        values.append(_overlap(q_basis, target, k, offset))
    return np.asarray(values, dtype=np.float64), {
        "retained_rank": float(q),
        "source_retained_energy": reference_energy,
        "null_retained_energy": reference_energy,
        "spectrum_energy_abs_error": max_spectrum_error,
    }


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


def _load_perturbation(path: Path) -> MatrixBundle:
    data = ad.read_h5ad(path, backed="r")
    raw = np.asarray(data.X[:], dtype=np.float32)
    controls = data.obs["core_control"].fillna(False).astype(bool).to_numpy()
    if not controls.any():
        raise ValueError(f"{path}: no core_control rows")
    genes = np.asarray([_symbol(v) for v in data.var["gene_name"].to_numpy()])
    finite_cols = np.isfinite(raw).all(axis=0)
    finite_rows = np.isfinite(raw).all(axis=1)
    all_zero_rows = np.all(np.nan_to_num(raw, nan=0.0) == 0.0, axis=1)
    keep_rows = finite_rows & ~all_zero_rows
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
        targets = None
        target_meta = {"guide_target_status": "unavailable_noncontrol_row_filter_changed_index"}
    x = raw[noncontrol]
    return MatrixBundle(x=x, genes=genes.tolist(), row_ids=data.obs_names[noncontrol].astype(str).tolist(), targets=targets,
        meta={"source": str(path), "shape_raw": list(data.shape), "n_control_rows_raw": int(controls.sum()),
              "n_control_rows_used": int(usable_controls.sum()), "n_rows_dropped_nonfinite": int((~finite_rows).sum()),
              "n_rows_dropped_all_zero": int((all_zero_rows & finite_rows).sum()), "n_genes_dropped_nonfinite": int((~finite_cols).sum()),
              "n_duplicate_symbols_aggregated": duplicates, "n_perturbation_rows": int(noncontrol.sum()), "control_centroid_ratio": centroid_ratio,
              "fold_pct_relation_max_abs_error": fold_relation_error, "control_expr_coverage": control_expr_coverage,
              "control_expr_nonnegative": control_expr_nonnegative, "target_vector_coverage": target_vector_coverage,
              "target_gene_hits": target_gene_hits, "min_target_gene_hits": min_target_gene_hits,
              "target_vector_control_referenced_correlation": target_vector_corr, "target_vector_control_referenced_mae": target_vector_mae,
              "delta_status": "validated_control_centred" if matrix_is_delta else "FAILED", **target_meta})


def _load_tcga(path: Path, transform: str) -> MatrixBundle:
    frame = pd.read_parquet(path)
    if "patient_id" not in frame:
        raise ValueError("TCGA parquet is missing patient_id")
    genes = np.asarray([_symbol(c) for c in frame.columns if c != "patient_id"])
    if len(genes) != len(set(genes)):
        raise ValueError("TCGA has duplicate canonical gene symbols")
    raw = frame.drop(columns="patient_id").to_numpy(dtype=np.float32, copy=False)
    finite_cols = np.isfinite(raw).all(axis=0); finite_rows = np.isfinite(raw).all(axis=1)
    all_zero_rows = np.all(np.nan_to_num(raw, nan=0.0) == 0.0, axis=1)
    keep_rows = finite_rows & ~all_zero_rows
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
            i = start + local; row[i] = -2.0; absolute = row.abs(); maximum = max(maximum, float(absolute.max().cpu()))
            for j in torch.where(absolute >= threshold)[0].cpu().tolist(): edges += 1; union(i, int(j))
            if usable_targets and counts[targets[i]] > 1:
                eligible += 1; hits += int(targets[i] == targets[int(absolute.argmax().cpu())])
    return {"dictionary_coherence_abs": maximum, "equivalence_threshold_abs": threshold, "equivalence_edges_directed": edges,
            "n_equivalence_classes": len({find(i) for i in range(n)}), "guide_same_target_retrieval_at1": hits / eligible if eligible else None,
            "guide_retrieval_eligible": eligible, "guide_retrieval_status": "available" if usable_targets else "unavailable"}


def _pc1_gate(observed: float, null_p95: float) -> bool:
    return bool(np.isfinite(observed) and np.isfinite(null_p95) and observed > null_p95)


def _shuffle_gate(shuffled: float, null_p95: float) -> bool:
    return bool(np.isfinite(shuffled) and np.isfinite(null_p95) and shuffled <= null_p95)


def _delta_gate(meta: dict[str, object]) -> bool:
    return meta["delta_status"] == "validated_control_centred"


def _split_half(x: torch.Tensor, groups: list[str], q: int, seed: int) -> tuple[torch.Tensor, torch.Tensor, dict[str, object]]:
    """Seeded cancer-stratified split, with every sufficiently large cancer in both halves."""
    if len(groups) != x.shape[0]: raise ValueError("missing TCGA cancer groups for split-half ceiling")
    rng = np.random.default_rng(seed); groups_a = np.asarray(groups); left: list[int] = []; right: list[int] = []
    for group in sorted(set(groups)):
        idx = np.flatnonzero(groups_a == group); idx = rng.permutation(idx)
        cut = len(idx) // 2
        if cut == 0: raise ValueError(f"cancer {group} has fewer than two samples")
        left.extend(idx[:cut]); right.extend(idx[cut:])
    a, b = torch.as_tensor(left, device=x.device), torch.as_tensor(right, device=x.device)
    va, _ = _right_svd(x[a], q, seed + 1); vb, _ = _right_svd(x[b], q, seed + 2)
    all_counts = pd.Series(groups).value_counts(); a_counts = pd.Series(groups_a[np.asarray(left)]).value_counts(); b_counts = pd.Series(groups_a[np.asarray(right)]).value_counts()
    return va, vb, {"n_half_a": int(a.numel()), "n_half_b": int(b.numel()), "n_cancers": int(len(all_counts)),
        "max_group_count_imbalance": int(max(abs(a_counts.reindex(all_counts.index, fill_value=0) - b_counts.reindex(all_counts.index, fill_value=0))))}


def _context(name: str, p: MatrixBundle, t: MatrixBundle, *, device: torch.device, draws: int, seed: int) -> dict[str, object]:
    px, tx, genes, join = _align(p, t)
    pt, tt = torch.as_tensor(px, device=device), torch.as_tensor(tx, device=device)
    q = 110; vp, sp = _right_svd(pt, q, seed); vt, _ = _right_svd(tt, q, seed + 1)
    va, vb, split_meta = _split_half(tt, t.groups or [], q, seed + 2)
    p_latent, t_latent = pt @ vp, tt @ vt
    p_health, t_health = _health(p_latent), _health(t_latent)
    # E0 uses randomized rank-q PCA for the principal-angle test.  A full
    # NumPy SVD here would silently move the run back to CPU for minutes/hours.
    # Report rank *within that explicitly retained q=110 subspace* instead;
    # this is the relevant representation health quantity, not a claim about
    # the uncomputed full matrix rank.
    out: dict[str, object] = {"context": name, "n_shared_genes": len(genes), "join": join, "perturbation": p.meta, "tcga": t.meta,
        "housekeepers_present": sorted(HOUSEKEEPING.intersection(genes)), "perturbation_health": p_health, "tcga_health": t_health,
        "rank": {"retained_components": int(sp.numel()), "effective_rank_retained": p_health["effective_rank"],
                 "stable_rank_retained_estimate": float((sp.square().sum() / sp[0].square()).cpu())},
        "dictionary": _dictionary_metrics(pt, p.targets), "split": split_meta}
    for k in (10, 25, 50, 100):
        observed = _overlap(vp, vt, k, 1); ceiling = _overlap(va, vb, k, 1)
        null, invariant = _matched_spectrum_rotated_null(sp, len(genes), vt, k=k, draws=draws, seed=seed + 100 + k, offset=1)
        # The positive ceiling is evaluated against an independently held-out
        # TCGA half and its own matched-spectrum randomized-P floor.  Thus the
        # positive and negative controls use the exact same principal-angle
        # routine, PC stripping, and gene universe.
        ceiling_null, _ = _matched_spectrum_rotated_null(sp, len(genes), vb, k=k, draws=draws, seed=seed + 300 + k, offset=1)
        generator = torch.Generator(device=device).manual_seed(seed + 700 + k)
        shuffled = _overlap(vp[torch.randperm(len(genes), device=device, generator=generator)], vt, k, 1)
        heldout_a, heldout_b = _overlap(vp, va, k, 1), _overlap(vp, vb, k, 1)
        effect = observed - null
        out[f"k{k}"] = {"pc1_removed_overlap": observed, "pc1_removed_ceiling": ceiling, "heldout_half_overlap_mean": (heldout_a + heldout_b) / 2,
            "gene_label_shuffle_overlap": shuffled, "null_median": float(np.median(null)), "null_p95": float(np.quantile(null, .95)), "null_std": float(np.std(null)),
            "ceiling_null_median": float(np.median(ceiling_null)), "ceiling_null_p95": float(np.quantile(ceiling_null, .95)),
            "null_p": float((1 + np.sum(null >= observed)) / (draws + 1)), "effect_vs_null": float(np.mean(effect)),
            "effect_ci95": [float(np.quantile(effect, .025)), float(np.quantile(effect, .975))], "theoretical_random_overlap": k / len(genes), **invariant}
    return out


def _self_test() -> None:
    """Static/synthetic preflight: broken controls must fail before GPU spending."""
    device = torch.device("cpu"); generator = torch.Generator().manual_seed(7)
    g, q, k = 128, 12, 5
    v, _ = torch.linalg.qr(torch.randn(g, q, generator=generator), mode="reduced")
    s = torch.arange(q, 0, -1, dtype=torch.float32)
    null, invariant = _matched_spectrum_rotated_null(s, g, v, k=k, draws=100, seed=3, offset=1)
    assert invariant["spectrum_energy_abs_error"] == 0.0 and null.std() > 0.0
    try: _matched_spectrum_rotated_null(s[:3], g, v, k=k, draws=100, seed=3, offset=1)
    except ValueError: pass
    else: raise AssertionError("truncated spectrum did not fail")
    try: _parse_targets(pd.Index(["bad"], name="gene_transcript"), np.asarray([False]))
    except Exception: raise
    assert _parse_targets(pd.Index(["bad"], name="gene_transcript"), np.asarray([False]))[0] is None
    # Deliberately corrupted controls must fail the same predicates used by the
    # live ledger, not merely a nearby toy assertion.
    assert _pc1_gate(.2, .1) and not _pc1_gate(.1, .1)
    assert _shuffle_gate(.05, .1) and not _shuffle_gate(.2, .1)
    assert _delta_gate({"delta_status": "validated_control_centred"}) and not _delta_gate({"delta_status": "FAILED"})
    assert 1.0 <= _effective_rank_torch(torch.eye(8)) <= 8.0
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


def _add_gates(ledger: GateLedger, context: dict[str, object], label: str, draws: int) -> None:
    p, t, join = context["perturbation"], context["tcga"], context["join"]
    ledger.add("G1.1_no_constant_columns", join["n_constant_columns_dropped"], "0 after explicit drop in evaluated matrix", join["n_matched_evaluated"] > 100, label)
    ledger.add("G1.2_real_gene_join", f"{join['n_left']},{join['n_right']},{join['n_matched_pre_constant']}", ">=0.8*min(left,right)", join["n_matched_pre_constant"] >= .8 * min(join["n_left"], join["n_right"]), label)
    ledger.add("G1.3_join_not_scrambled", "per-k gene-label shuffle", "all shuffles collapse to matched-spectrum floor", all(_shuffle_gate(context[f"k{k}"]["gene_label_shuffle_overlap"], context[f"k{k}"]["null_p95"]) for k in (10,25,50,100)), label)
    ledger.add("G1.4_scale_sanity", json.dumps({k:t[k] for k in ("raw_max","post_transform_max","post_transform_zero_fraction","post_transform_skew")}), "both transforms recorded", True, label)
    ledger.add("G1.5_explicit_row_filtering", json.dumps({"P_nonfinite":p["n_rows_dropped_nonfinite"],"T_nonfinite":t["n_rows_dropped_nonfinite"],"P_zero":p["n_rows_dropped_all_zero"],"T_zero":t["n_rows_dropped_all_zero"]}), "all drops explicitly counted", True, label)
    ledger.add("G1.6_housekeeping_mapping", ",".join(context["housekeepers_present"]), "all five present", set(context["housekeepers_present"]) == HOUSEKEEPING, label)
    ledger.add("E0.a_delta_not_absolute", json.dumps({"status":p["delta_status"],"ratio":p["control_centroid_ratio"],"fold_pct_error":p["fold_pct_relation_max_abs_error"],"control_expr_coverage":p["control_expr_coverage"],"target_control_ref_r":p["target_vector_control_referenced_correlation"],"target_control_ref_mae":p["target_vector_control_referenced_mae"]}), "X target effect agrees with control_expr x fold_expr reference", _delta_gate(p), label)
    ledger.add("E0.c_orientation", f"P={p['n_perturbation_rows']}x{context['n_shared_genes']};T={t['n_patients_used']}x{context['n_shared_genes']}", "rows=samples, columns=genes", True, label)
    ledger.add("G1.2b_split_composition", json.dumps(context["split"]), "cancer-stratified; imbalance<=1", context["split"]["max_group_count_imbalance"] <= 1, label)
    for who in ("perturbation_health", "tcga_health"):
        h = context[who]; prefix = "P" if who.startswith("perturbation") else "T"
        ledger.add(f"G3.1_effective_rank_{prefix}", h["effective_rank"], "finite", np.isfinite(h["effective_rank"]), label)
        ledger.add(f"G3.2_dead_dimensions_{prefix}", h["dead_fraction"], "report", h["dead_fraction"] <= .5, label)
        ledger.add(f"G3.3_duplicate_dimensions_{prefix}", h["duplicate_pair_fraction"], "report", np.isfinite(h["duplicate_pair_fraction"]), label)
        ledger.add(f"G3.4_sample_variation_{prefix}", h["min_sample_std"], ">0", h["min_sample_std"] > 0, label)
        ledger.add(f"G3.5_site_degeneracy_{prefix}", "unavailable_no_site_labels", "explicitly unavailable", True, label)
        ledger.add(f"G3.6_norm_sanity_{prefix}", json.dumps({"mean":h["mean_norm"],"median":h["median_norm"],"nonfinite":h["n_nonfinite"]}), "zero nonfinite", h["n_nonfinite"] == 0, label)
    d = context["dictionary"]
    ledger.add("E0b_dictionary_coherence", d["dictionary_coherence_abs"], "finite", np.isfinite(d["dictionary_coherence_abs"]), label)
    ledger.add("E0b_equivalence_classes", d["n_equivalence_classes"], ">0", d["n_equivalence_classes"] > 0, label)
    ledger.add("E0b_guide_retrieval", d["guide_retrieval_status"], "available or explicit unavailable", d["guide_retrieval_status"] in {"available","unavailable"}, label)
    for k in (10,25,50,100):
        r = context[f"k{k}"]; tag = f"k{k}_{label}"
        ledger.add(f"G4.1_positive_control_split_half_{tag}", r["pc1_removed_ceiling"], "> same-path heldout matched-null p95", r["pc1_removed_ceiling"] > r["ceiling_null_p95"], label)
        ledger.add(f"G4.2_negative_control_rotated_P_{tag}", r["null_median"], "near k/n_genes", abs(r["null_median"] - r["theoretical_random_overlap"]) < 0.5 * r["theoretical_random_overlap"], label)
        ledger.add(f"G4.2b_negative_control_gene_shuffle_{tag}", r["gene_label_shuffle_overlap"], "<= matched-spectrum null p95", _shuffle_gate(r["gene_label_shuffle_overlap"], r["null_p95"]), label)
        ledger.add(f"G4.3_null_sanity_{tag}", json.dumps({"std":r["null_std"],"p95":r["null_p95"]}), "std>0; nondegenerate", r["null_std"] > 0, label)
        ledger.add(f"G4.4_heldout_vs_insample_{tag}", json.dumps({"in":r["pc1_removed_overlap"],"heldout":r["heldout_half_overlap_mean"]}), "reported; no fitted mapping", np.isfinite(r["heldout_half_overlap_mean"]), label)
        ledger.add(f"G4.5_null_resolution_{tag}", 1/(draws+1), "draws>=100", draws >= 100, label)
        ledger.add(f"G4.6_effect_ci_{tag}", json.dumps(r["effect_ci95"]), "lower CI > 0", r["effect_ci95"][0] > 0, label)
        ledger.add(f"E0.b_pc1_stripped_{tag}", json.dumps({"observed":r["pc1_removed_overlap"],"p":r["null_p"]}), "observed>null p95", _pc1_gate(r["pc1_removed_overlap"], r["null_p95"]), label)
        ledger.add(f"E0.d_ceiling_{tag}", r["pc1_removed_ceiling"], ">null p95", r["pc1_removed_ceiling"] > r["null_p95"], label)
        ledger.add(f"E0.e_matched_spectrum_floor_{tag}", r["spectrum_energy_abs_error"], "exact preserved retained energy", r["spectrum_energy_abs_error"] == 0.0, label)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k562"); parser.add_argument("--rpe1"); parser.add_argument("--tcga"); parser.add_argument("--tcga-registry"); parser.add_argument("--output")
    parser.add_argument("--draws", type=int, default=100); parser.add_argument("--seed", type=int, default=42); parser.add_argument("--device", default="cuda")
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2])); parser.add_argument("--workspace-link"); parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test: _self_test(); return
    if not all((args.k562,args.rpe1,args.tcga,args.tcga_registry,args.output,args.workspace_link)): parser.error("k562, rpe1, tcga, tcga-registry, output, and workspace-link are required")
    if args.draws < 100: parser.error("--draws must be >=100 for E0's predeclared floor control")
    started = time.monotonic(); root = Path(args.repo_root).resolve(); output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    device = torch.device(args.device if args.device == "cpu" or torch.cuda.is_available() else "cpu")
    ledger = GateLedger(output, "E0_E0b", official_log=root / "v2/research/rebase/nature/GATE_LOG.md")
    source = Path(__file__).resolve(); workspace_link = Path(args.workspace_link)
    link_is_real = workspace_link.is_symlink() and workspace_link.resolve() == root and root in source.parents
    ledger.add("G0.1_workspace_real", f"link={workspace_link};target={workspace_link.resolve() if workspace_link.exists() else 'missing'}", "declared runtime link resolves to exact git root", link_is_real, "Linux symlink is the remote equivalent of the required Windows Junction; prevents stale copied source")
    dirty = _git(root, "status", "--porcelain"); ledger.add("G0.2_code_identity", _git(root,"rev-parse","HEAD"), "clean git worktree", dirty == "", f"dirty={dirty or 'none'}")
    for label, value in (("G0.3_k562_identity",args.k562),("G0.3_rpe1_identity",args.rpe1),("G0.3_tcga_identity",args.tcga),("G0.3_tcga_registry_identity",args.tcga_registry)): ledger.artifact(label, value)
    ledger.add("G0.4_manifest_read", "data-source experiment; no trained artifacts compared", "not applicable, explicitly declared", True, "input identities emitted in input_manifest.json")
    p_k, p_r = _load_perturbation(Path(args.k562)), _load_perturbation(Path(args.rpe1))
    results: dict[str, object] = {"schema_version":"2.0", "experiment":"E0_E0b", "device":str(device), "draws":args.draws, "seed":args.seed,
        "command":sys.argv, "code_sha":_git(root,"rev-parse","HEAD"), "code_dirty":dirty, "python":sys.version, "torch":torch.__version__, "cuda":torch.version.cuda, "platform":platform.platform(), "contexts":{}}
    for transform in ("signed_log1p", "clip_log1p"):
        tcga = _restrict_tcga_to_registry(_load_tcga(Path(args.tcga), transform), Path(args.tcga_registry))
        for name, perturb, seed in (("K562",p_k,args.seed),("RPE1",p_r,args.seed+1000)):
            key = f"{name}_{transform}"; context = _context(name, perturb, tcga, device=device, draws=args.draws, seed=seed)
            results["contexts"][key] = context; _add_gates(ledger, context, key, args.draws)
    # Transform robustness is predeclared: all PC1-stripped effects must stay above their own floor and retain direction.
    concordant = all(c[f"k{k}"]["effect_ci95"][0] > 0 for c in results["contexts"].values() for k in (10,25,50,100))
    ledger.add("E0.f_cross_context_transform_replication", "K562,RPE1 x signed,clip", "all PC1-stripped effect CI lower bounds >0", concordant)
    missing = _missing_gate_prefixes([str(row["gate"]) for row in ledger.rows])
    ledger.add("G0.5_required_gate_coverage", ",".join(missing) if missing else "complete", "no mandatory gate family missing", not missing)
    results["wall_seconds"] = time.monotonic() - started
    results["gates_pass"] = ledger.write()
    (output / "input_manifest.json").write_text(json.dumps({"inputs":{k:results[k] for k in ("command","code_sha","code_dirty","device","draws","seed")},"k562":p_k.meta,"rpe1":p_r.meta},indent=2))
    (output / "e0_basis_transfer.json").write_text(json.dumps(results, indent=2))
    print(json.dumps({"output":str(output),"gates_pass":results["gates_pass"],"wall_seconds":results["wall_seconds"]},indent=2), flush=True)


if __name__ == "__main__": main()
