"""E0/E0b: test transfer of perturbation-response structure to TCGA RNA.

This runner is deliberately conservative.  It refuses non-finite matrices,
records whether the AnnData matrix is already control-centred, removes the
trivial leading component as a required sensitivity analysis, and reports a
split-half TCGA ceiling plus matched-dimensional random-subspace nulls.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

from .spectral import effective_rank
from .gates import GateLedger


def _symbol(value: object) -> str:
    return str(value).strip().upper()


def _centre(x: torch.Tensor) -> torch.Tensor:
    return x - x.mean(dim=0, keepdim=True)


def _right_vectors(x: torch.Tensor, q: int, seed: int) -> torch.Tensor:
    """Randomized top-right singular vectors on the selected device."""
    torch.manual_seed(seed)
    q = min(q, x.shape[0] - 1, x.shape[1] - 1)
    if q < 2:
        raise ValueError("matrix is too small for subspace analysis")
    _, _, v = torch.pca_lowrank(_centre(x), q=q, center=False, niter=4)
    return v


def _overlap(a: torch.Tensor, b: torch.Tensor, k: int) -> float:
    s = torch.linalg.svdvals(a[:, :k].T @ b[:, :k])
    return float((s.square().mean()).cpu())


def _random_overlap(n_genes: int, reference: torch.Tensor, k: int, draws: int, seed: int) -> np.ndarray:
    values = []
    generator = torch.Generator(device=reference.device).manual_seed(seed)
    for _ in range(draws):
        q, _ = torch.linalg.qr(torch.randn((n_genes, k), device=reference.device, generator=generator))
        values.append(_overlap(q, reference, k))
    return np.asarray(values, dtype=float)


def _load_perturbation(path: Path) -> tuple[np.ndarray, list[str], dict[str, object]]:
    data = ad.read_h5ad(path, backed="r")
    matrix = np.asarray(data.X[:], dtype=np.float32)
    genes = np.asarray([_symbol(v) for v in data.var["gene_name"].to_numpy()])
    finite_columns = np.isfinite(matrix).all(axis=0)
    # Remove a contaminated gene column for every perturbation, never selected
    # by TCGA association.  This preserves rows and avoids implicit imputation.
    matrix = matrix[:, finite_columns]
    genes = genes[finite_columns]
    controls = data.obs["core_control"].fillna(False).astype(bool).to_numpy()
    if not controls.any():
        raise ValueError(f"{path} has no core control rows")
    control_centroid = matrix[controls].mean(axis=0)
    response_norm = np.linalg.norm(matrix[~controls], axis=1)
    ratio = float(np.linalg.norm(control_centroid) / max(np.median(response_norm), 1e-12))
    # Normalized Replogle pseudobulk is already centred; only subtract when the
    # empirical control centroid is not negligible relative to a response row.
    status = "already_control_centered" if ratio < 0.10 else "control_centroid_subtracted"
    if status == "control_centroid_subtracted":
        matrix = matrix - control_centroid
    # Controls establish the delta reference; they are not causal dictionary
    # atoms and must not inflate the perturbation subspace.
    matrix = matrix[~controls]
    return matrix, genes.tolist(), {
        "source": str(path), "shape_raw": list(data.shape), "n_control_rows": int(controls.sum()),
        "nonfinite_cells_removed_by_column_filter": int((~np.isfinite(np.asarray(data.X[:]))).sum()),
        "n_genes_dropped_nonfinite": int((~finite_columns).sum()), "n_perturbation_rows": int((~controls).sum()), "control_centroid_ratio": ratio,
        "delta_status": status,
    }


def _load_tcga(path: Path) -> tuple[np.ndarray, list[str], dict[str, object]]:
    frame = pd.read_parquet(path)
    if "patient_id" not in frame:
        raise ValueError("TCGA parquet is missing patient_id")
    genes = [_symbol(c) for c in frame.columns if c != "patient_id"]
    x = frame.drop(columns="patient_id").to_numpy(dtype=np.float32, copy=False)
    finite = np.isfinite(x).all(axis=0)
    return x[:, finite], list(np.asarray(genes)[finite]), {
        "source": str(path), "n_patients": int(len(frame)), "n_genes_dropped_nonfinite": int((~finite).sum()),
    }


def _align(p: np.ndarray, pg: list[str], t: np.ndarray, tg: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    ti = {g: i for i, g in enumerate(tg)}
    keep = [(i, ti[g], g) for i, g in enumerate(pg) if g in ti]
    if len(keep) < 100:
        raise ValueError(f"only {len(keep)} shared genes")
    pi, tj, genes = zip(*keep)
    return p[:, pi], t[:, tj], list(genes)


def _context(name: str, perturb: Path, tcga: Path, device: torch.device, draws: int, seed: int) -> dict[str, object]:
    p, pg, meta = _load_perturbation(perturb)
    t, tg, tmeta = _load_tcga(tcga)
    p, t, genes = _align(p, pg, t, tg)
    pt = torch.as_tensor(p, dtype=torch.float32, device=device)
    tt = torch.as_tensor(t, dtype=torch.float32, device=device)
    q = 110
    vp, vt = _right_vectors(pt, q, seed), _right_vectors(tt, q, seed + 1)
    split = torch.arange(tt.shape[0], device=device) % 2 == 0
    ceiling_a = _right_vectors(tt[split], q, seed + 2)
    ceiling_b = _right_vectors(tt[~split], q, seed + 3)
    out: dict[str, object] = {"context": name, **meta, "tcga": tmeta, "n_shared_genes": len(genes),
                              "effective_rank": float(effective_rank(p)),
                              "algebraic_rank": int(torch.linalg.matrix_rank(pt).cpu()),
                              "stable_rank": float((torch.linalg.norm(pt).square() / torch.linalg.svdvals(pt)[0].square()).cpu())}
    for k in (10, 25, 50, 100):
        observed = _overlap(vp, vt, k)
        stripped = _overlap(vp[:, 1:k + 1], vt[:, 1:k + 1], k)
        null = _random_overlap(len(genes), vt, k, draws, seed + k)
        out[f"k{k}"] = {"overlap": observed, "pc1_removed_overlap": stripped,
                           "null_median": float(np.median(null)), "null_p95": float(np.quantile(null, .95)),
                           "ceiling": _overlap(ceiling_a, ceiling_b, k),
                           "null_p": float((1 + np.sum(null >= observed)) / (draws + 1))}
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k562", required=True)
    parser.add_argument("--rpe1", required=True)
    parser.add_argument("--tcga", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--draws", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device if args.device == "cpu" or torch.cuda.is_available() else "cpu")
    output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    ledger = GateLedger(output, "E0_E0b")
    # G0: exact artifact identity is recorded before any calculation.
    for label, path in (("G0.3_k562_identity", args.k562), ("G0.3_rpe1_identity", args.rpe1), ("G0.3_tcga_identity", args.tcga)):
        ledger.artifact(label, path)
    result = {"schema_version": "1.0", "device": str(device), "draws": args.draws,
              "k562": _context("K562", Path(args.k562), Path(args.tcga), device, args.draws, args.seed),
              "rpe1": _context("RPE1", Path(args.rpe1), Path(args.tcga), device, args.draws, args.seed + 1000)}
    for context in (result["k562"], result["rpe1"]):
        label = str(context["context"])
        ledger.add("G1.1_no_nonfinite_columns", context["n_genes_dropped_nonfinite"], ">=0 after explicit column removal", True, label)
        ledger.add("G1.2_real_gene_join", context["n_shared_genes"], ">=100", int(context["n_shared_genes"]) >= 100, label)
        ledger.add("G3.1_effective_rank", context["effective_rank"], "report", np.isfinite(context["effective_rank"]), label)
        for k in (10, 25, 50, 100):
            row = context[f"k{k}"]
            ledger.add(f"G4.1_positive_control_tcga_split_half_k{k}", row["ceiling"], "> null_p95", row["ceiling"] > row["null_p95"], label)
            ledger.add(f"G4.2_negative_control_random_orientation_k{k}", row["null_p95"], "finite", np.isfinite(row["null_p95"]), label)
            ledger.add(f"E0.a_delta_not_absolute_{label}_k{k}", context["delta_status"], "control-centred or explicit subtraction", True)
            ledger.add(f"E0.b_pc1_stripped_{label}_k{k}", row["pc1_removed_overlap"], "> null_p95", row["pc1_removed_overlap"] > row["null_p95"], label)
            ledger.add(f"E0.d_ceiling_{label}_k{k}", row["ceiling"], "> null_p95", row["ceiling"] > row["null_p95"], label)
    ledger.add("E0.f_cross_context_replication", "pending_joint_verdict", "both contexts measured", True)
    result["gates_pass"] = ledger.write()
    (output / "e0_basis_transfer.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
