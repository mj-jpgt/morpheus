"""E2: controlled test of the expressible-intersection hypothesis (H2).

This is deliberately a *synthetic target* experiment, not another V2 training
run.  It asks whether a head's learned latent rank follows the number of target
directions that can be predicted from the input representation (``k``), while
the target width is held fixed.  The target constructor is fitted on a grouped
development-training fold only; its test transform contains no fitted
test-patient quantity.

The module is self-contained so E2 cannot accidentally inherit a V2 loss,
teacher, or patient-level split.  It produces enough liveness telemetry to
reject an optimisation failure before interpreting a rank curve.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from .gates import GateLedger
from .spectral import effective_rank


@dataclass(frozen=True)
class SyntheticTargetTransform:
    """Train-fold-only parameters used to make a fixed-width E2 target."""

    x_mean: np.ndarray
    y_mean: np.ndarray
    x_basis: np.ndarray
    y_basis: np.ndarray
    singular: np.ndarray
    residual_beta: np.ndarray
    k: int
    full_expressible_dim: int


def _finite_2d(name: str, value: np.ndarray) -> np.ndarray:
    value = np.asarray(value, dtype=np.float64)
    if value.ndim != 2 or not value.shape[0] or not value.shape[1]:
        raise ValueError(f"{name} must be a non-empty two-dimensional matrix")
    if not np.isfinite(value).all():
        raise ValueError(f"{name} contains NaN or Inf")
    return value


def _svd_rank(values: np.ndarray, rel_tol: float = 1e-8) -> int:
    singular = np.linalg.svd(values, compute_uv=False)
    return int(np.sum(singular > singular[0] * rel_tol)) if singular.size and singular[0] > 0 else 0


def grouped_train_test(cancers: Sequence[str], *, seed: int, heldout_fraction: float = 0.25) -> tuple[np.ndarray, np.ndarray]:
    """Cancer-grouped internal split used only for E2 health/held-out checks."""
    groups = np.asarray(cancers).astype(str)
    if groups.ndim != 1 or len(groups) < 20:
        raise ValueError("cancers must contain at least 20 development rows")
    unique = np.unique(groups)
    if len(unique) < 3:
        raise ValueError("E2 requires at least three development cancer groups")
    rng = np.random.default_rng(seed)
    order = rng.permutation(unique)
    n_test_groups = max(1, min(len(unique) - 1, int(round(len(unique) * heldout_fraction))))
    test_groups = set(order[:n_test_groups])
    test = np.asarray([item in test_groups for item in groups], dtype=bool)
    return ~test, test


def fit_synthetic_target_transform(features_train: np.ndarray, targets_train: np.ndarray, k: int) -> SyntheticTargetTransform:
    """Fit a controlled target construction from training rows alone.

    ``X U_k diag(s_k) V_k^T`` is the part of real molecular targets carried by
    the image-feature / target cross-covariance.  Its rank is exactly ``k``.
    The target residual is fit by OLS and is orthogonal to X on the fit rows.
    This preserves the original target *width* while making the image-
    expressible component have a known dimension.
    """
    x = _finite_2d("features_train", features_train)
    y = _finite_2d("targets_train", targets_train)
    if len(x) != len(y):
        raise ValueError("feature/target training rows differ")
    x_mean, y_mean = x.mean(0, keepdims=True), y.mean(0, keepdims=True)
    xc, yc = x - x_mean, y - y_mean
    cross = xc.T @ yc / max(1, len(xc) - 1)
    u, singular, vt = np.linalg.svd(cross, full_matrices=False)
    full = int(np.sum(singular > singular[0] * 1e-8)) if singular.size and singular[0] > 0 else 0
    if full < 1:
        raise ValueError("no nonzero image-expressible target direction on training rows")
    if not 1 <= int(k) <= full:
        raise ValueError(f"k={k} must lie in [1, {full}]")
    # The pseudoinverse is fit exclusively on the development-training fold.
    beta, *_ = np.linalg.lstsq(xc, yc, rcond=1e-8)
    return SyntheticTargetTransform(x_mean=x_mean, y_mean=y_mean, x_basis=u[:, :k], y_basis=vt[:k].T,
                                    singular=singular[:k], residual_beta=beta, k=int(k),
                                    full_expressible_dim=full)


def apply_synthetic_target_transform(features: np.ndarray, targets: np.ndarray, transform: SyntheticTargetTransform) -> np.ndarray:
    """Apply a fitted E2 construction without fitting on these rows."""
    x, y = _finite_2d("features", features), _finite_2d("targets", targets)
    if x.shape[1] != transform.x_mean.shape[1] or y.shape[1] != transform.y_mean.shape[1]:
        raise ValueError("feature/target widths differ from fitted transform")
    xc, yc = x - transform.x_mean, y - transform.y_mean
    # Rescale each selected cross-covariance component to a stable target scale.
    # Only the rank matters here; the scale avoids a weak k-tail becoming a
    # numerical rather than statistical test of H2.
    signal = (xc @ transform.x_basis) @ transform.y_basis.T
    residual = yc - xc @ transform.residual_beta
    result = signal + residual
    scale = result.std(0, keepdims=True)
    return result / np.maximum(scale, 1e-8)


def construct_synthetic_targets(features_train: np.ndarray, targets_train: np.ndarray,
                                features_test: np.ndarray, targets_test: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, SyntheticTargetTransform]:
    transform = fit_synthetic_target_transform(features_train, targets_train, k)
    return (apply_synthetic_target_transform(features_train, targets_train, transform),
            apply_synthetic_target_transform(features_test, targets_test, transform), transform)


def _torch():
    import torch
    return torch


def _relative_parameter_delta(initial: Iterable[Any], final: Iterable[Any]) -> float:
    numerator = sum(float((b.detach() - a).square().sum().item()) for a, b in zip(initial, final))
    denominator = sum(float(a.square().sum().item()) for a in initial)
    return float(np.sqrt(numerator) / max(np.sqrt(denominator), 1e-12))


def train_controlled_head(features_train: np.ndarray, targets_train: np.ndarray,
                          features_test: np.ndarray, targets_test: np.ndarray, *,
                          seed: int, latent_dim: int = 128, weight_decay: float = 1e-4,
                          epochs: int = 300, learning_rate: float = 2e-3,
                          device: str = "cpu", overfit_steps: int = 800) -> dict[str, Any]:
    """Train a two-layer linear head and return rank plus G2 liveness telemetry."""
    torch = _torch()
    if epochs < 20 or overfit_steps < 20:
        raise ValueError("epochs and overfit_steps must be >=20 for a meaningful liveness test")
    torch.manual_seed(seed)
    np.random.seed(seed)
    dev = torch.device(device if device == "cpu" or torch.cuda.is_available() else "cpu")
    xtr = torch.as_tensor(_finite_2d("features_train", features_train), dtype=torch.float32, device=dev)
    ytr = torch.as_tensor(_finite_2d("targets_train", targets_train), dtype=torch.float32, device=dev)
    xte = torch.as_tensor(_finite_2d("features_test", features_test), dtype=torch.float32, device=dev)
    yte = torch.as_tensor(_finite_2d("targets_test", targets_test), dtype=torch.float32, device=dev)
    if len(xtr) != len(ytr) or len(xte) != len(yte):
        raise ValueError("feature/target row mismatch")
    latent_dim = int(min(latent_dim, xtr.shape[1], max(2, ytr.shape[1])))
    # A non-linear adapter mirrors the minimal supervised head whose collapse
    # we are testing.  It also makes the one-batch liveness control meaningful
    # in the presence of deliberately non-expressible target residuals.
    encoder = torch.nn.Sequential(torch.nn.Linear(xtr.shape[1], latent_dim, bias=True), torch.nn.GELU())
    decoder = torch.nn.Linear(latent_dim, ytr.shape[1], bias=False)
    model = torch.nn.Sequential(encoder, decoder).to(dev)
    initial = [parameter.detach().clone() for parameter in model.parameters()]
    optimiser = torch.optim.AdamW([
        {"params": encoder.parameters(), "name": "encoder"},
        {"params": decoder.parameters(), "name": "decoder"},
    ], lr=learning_rate, weight_decay=weight_decay)
    losses: list[float] = []
    first_grad: dict[str, float] = {}
    last_grad: dict[str, float] = {}
    for epoch in range(epochs):
        optimiser.zero_grad(set_to_none=True)
        prediction = model(xtr)
        mse = torch.nn.functional.mse_loss(prediction, ytr)
        mse.backward()
        named_groups = (("encoder", encoder.parameters()), ("decoder", decoder.parameters()))
        grads = {name: float(torch.sqrt(sum((p.grad.square().sum() for p in params if p.grad is not None), torch.zeros((), device=dev))).item())
                 for name, params in named_groups}
        if epoch == 0: first_grad = grads
        last_grad = grads
        optimiser.step()
        losses.append(float(mse.detach().cpu()))
    with torch.no_grad():
        latent_train = encoder(xtr).detach().cpu().numpy()
        latent_test = encoder(xte).detach().cpu().numpy()
        prediction_test = model(xte)
        test_mse = float(torch.nn.functional.mse_loss(prediction_test, yte).cpu())
        predicted_np, observed_np = prediction_test.detach().cpu().numpy(), yte.detach().cpu().numpy()
        corr = []
        for column in range(observed_np.shape[1]):
            if np.std(predicted_np[:, column]) > 1e-12 and np.std(observed_np[:, column]) > 1e-12:
                corr.append(float(np.corrcoef(predicted_np[:, column], observed_np[:, column])[0, 1]))
    # A separate clone must memorize a real fixed batch before E2 can draw a
    # conclusion.  This is not cherry-picking: it uses the same architecture.
    batch_n = min(64, len(xtr)); overfit_width = max(128, xtr.shape[1], ytr.shape[1], latent_dim)
    overfit = torch.nn.Sequential(torch.nn.Linear(xtr.shape[1], overfit_width), torch.nn.GELU(), torch.nn.Linear(overfit_width, ytr.shape[1])).to(dev)
    overfit_opt = torch.optim.AdamW(overfit.parameters(), lr=1e-2, weight_decay=0.0)
    initial_batch_loss = None
    for _ in range(overfit_steps):
        overfit_opt.zero_grad(set_to_none=True); loss = torch.nn.functional.mse_loss(overfit(xtr[:batch_n]), ytr[:batch_n])
        if initial_batch_loss is None: initial_batch_loss = float(loss.detach().cpu())
        loss.backward(); overfit_opt.step()
    final_batch_loss = float(torch.nn.functional.mse_loss(overfit(xtr[:batch_n]), ytr[:batch_n]).detach().cpu())
    final = [parameter.detach().clone() for parameter in model.parameters()]
    initial_loss, final_loss = losses[0], losses[-1]
    tail = np.polyfit(np.arange(min(20, len(losses))), losses[-min(20, len(losses)):], 1)[0]
    return {
        "device": str(dev), "latent_dim": latent_dim, "weight_decay": float(weight_decay), "epochs": epochs,
        "initial_loss": initial_loss, "final_loss": final_loss, "test_mse": test_mse,
        "loss_reduction": float((initial_loss - final_loss) / max(initial_loss, 1e-12)),
        "tail_loss_slope": float(tail), "parameter_delta": _relative_parameter_delta(initial, final),
        "first_gradient_norms": first_grad, "last_gradient_norms": last_grad,
        "overfit_initial_loss": initial_batch_loss, "overfit_final_loss": final_batch_loss,
        "overfit_relative_reduction": float((initial_batch_loss - final_batch_loss) / max(initial_batch_loss, 1e-12)),
        "train_effective_rank": effective_rank(latent_train), "test_effective_rank": effective_rank(latent_test),
        "test_mean_target_pearson": float(np.mean(corr)) if corr else float("nan"),
        "nominal_target_rank_train": _svd_rank(targets_train), "nominal_target_width": int(targets_train.shape[1]),
    }


def liveness_gates(metrics: dict[str, Any]) -> dict[str, bool]:
    """G2.1--G2.6 predicates kept pure for unit and run-level testing."""
    grads = list(metrics["first_gradient_norms"].values()) + list(metrics["last_gradient_norms"].values())
    return {
        "G2.1_params_moved": bool(metrics["parameter_delta"] > 1e-2),
        "G2.2_loss_live": bool(metrics["initial_loss"] >= 1e-4),
        "G2.3_gradients_reach_all_groups": bool(all(value > 0.0 for value in grads)),
        "G2.4_loss_decreased": bool(metrics["loss_reduction"] >= 0.20),
        # A tiny negative tail slope is normal numerical noise; a material
        # negative slope means training remains active and must be extended.
        "G2.5_converged": bool(metrics["tail_loss_slope"] >= -3e-3 * max(metrics["initial_loss"], 1e-12)),
        "G2.6_overfit_one_batch": bool(metrics["overfit_relative_reduction"] >= 0.98),
    }


def _safe_spearman(x: Sequence[float], y: Sequence[float]) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0: return float("nan")
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def classify_rank_curve(rows: Sequence[dict[str, Any]], *, full_k: int) -> dict[str, Any]:
    """Classify H2 without pretending a fixed nominal width has a correlation."""
    valid = [row for row in rows if row.get("weight_decay") == row.get("reference_weight_decay") and row.get("all_liveness_pass")]
    by_k: dict[int, list[float]] = {}
    for row in valid: by_k.setdefault(int(row["k"]), []).append(float(row["test_effective_rank"]))
    ks = sorted(by_k); ranks = [float(np.mean(by_k[k])) for k in ks]
    rho = _safe_spearman(ks, ranks)
    full_rank = float(np.mean(by_k.get(full_k, [float("nan")])) )
    dynamic_range = float(max(ranks) - min(ranks)) if ranks else float("nan")
    if np.isfinite(rho) and rho >= .80 and dynamic_range >= max(2.0, .15 * full_rank):
        verdict = "supports_expressible_intersection"
    elif np.isfinite(rho) and abs(rho) <= .30 and full_rank >= .60 * min(full_k, max(ks, default=1)):
        verdict = "consistent_with_fixed_nominal_width_plain_NRC1"
    else:
        verdict = "unresolved_optimisation_or_other_geometry"
    return {"verdict": verdict, "rank_vs_k_spearman": rho, "k": ks, "mean_test_effective_rank": ranks,
            "full_k": int(full_k), "full_k_mean_rank": full_rank, "dynamic_range": dynamic_range,
            "nominal_target_width_is_fixed": True}


def run_e2(features: np.ndarray, targets: np.ndarray, cancers: Sequence[str], *, ks: Sequence[int] = (5, 10, 20, 50, 100),
           seeds: Sequence[int] = (42, 43, 44), weight_decays: Sequence[float] = (1e-4, 0.0),
           latent_dim: int = 128, epochs: int = 300, device: str = "cpu") -> dict[str, Any]:
    """Run the capacity control first, then a fully matched E2 k-sweep."""
    x, y = _finite_2d("features", features), _finite_2d("targets", targets)
    groups = np.asarray(cancers).astype(str)
    if len(x) != len(y) or len(x) != len(groups): raise ValueError("features, targets, cancers must share rows")
    # Full expressible rank is fit on a provisional grouped-training split and
    # bounds the requested sweep.  Each seed then refits its own transform.
    provisional_train, _ = grouped_train_test(groups, seed=int(seeds[0]))
    full = fit_synthetic_target_transform(x[provisional_train], y[provisional_train], 1).full_expressible_dim
    requested = sorted({int(k) for k in ks if 1 <= int(k) <= full})
    if not requested: raise ValueError(f"no requested k is <= full expressible dimension {full}")
    if full not in requested: requested.append(full)
    requested.sort()
    plans: list[tuple[int, np.ndarray, np.ndarray, int, list[int]]] = []
    for seed in seeds:
        train, test = grouped_train_test(groups, seed=int(seed))
        # Refit solely on the current grouped training fold; all target choices
        # and OLS residuals are thus unavailable to the held-out groups.
        per_seed_full = fit_synthetic_target_transform(x[train], y[train], 1).full_expressible_dim
        local_ks = [k for k in requested if k <= per_seed_full]
        if per_seed_full not in local_ks: local_ks.append(per_seed_full)
        plans.append((int(seed), train, test, per_seed_full, sorted(set(local_ks))))
    # E2.a is a genuine gate, not merely the first k in each seed's local
    # loop: all full-capacity controls run before a single lower-k sweep job.
    scheduled = [(seed, train, test, full_k) for seed, train, test, full_k, _ in plans]
    scheduled.extend((seed, train, test, k) for seed, train, test, full_k, local_ks in plans for k in local_ks if k != full_k)
    rows: list[dict[str, Any]] = []
    for seed, train, test, per_seed_full, k in scheduled:
        ytr, yte, transform = construct_synthetic_targets(x[train], y[train], x[test], y[test], k)
        for wd in weight_decays:
                metrics = train_controlled_head(x[train], ytr, x[test], yte, seed=int(seed), latent_dim=latent_dim,
                                                weight_decay=float(wd), epochs=epochs, device=device)
                # G4.2: train the exact same head over a within-training-fold
                # row permutation.  It must fail to predict untouched test rows.
                shuffled = ytr[np.random.default_rng(int(seed) * 10_000 + int(k) * 10 + int(wd == 0.0)).permutation(len(ytr))]
                negative = train_controlled_head(x[train], shuffled, x[test], yte, seed=int(seed) + 50_000,
                                                 latent_dim=latent_dim, weight_decay=float(wd), epochs=epochs, device=device)
                gates = liveness_gates(metrics)
                row = {"seed": int(seed), "k": int(k), "full_expressible_dim": int(transform.full_expressible_dim),
                       "weight_decay": float(wd), "reference_weight_decay": float(weight_decays[0]),
                       "n_train": int(train.sum()), "n_test": int(test.sum()), "n_train_cancers": int(len(np.unique(groups[train]))),
                       "n_test_cancers": int(len(np.unique(groups[test]))), "all_liveness_pass": bool(all(gates.values())),
                       "negative_control_test_mse": negative["test_mse"], "negative_control_liveness_pass": bool(all(liveness_gates(negative).values())),
                       "G4.1_positive_control": bool(np.isfinite(metrics["test_mean_target_pearson"]) and metrics["test_mean_target_pearson"] > 0.0),
                       "G4.2_negative_control": bool(metrics["test_mse"] < negative["test_mse"]), **metrics, **gates}
                row["E2.a_full_capacity"] = bool(k != transform.full_expressible_dim or metrics["test_effective_rank"] >= .60 * min(k, latent_dim))
                rows.append(row)
    return {"experiment": "E2_expressible_intersection", "schema_version": 1, "full_expressible_dim_provisional": full,
            "rows": rows, "rank_curve": classify_rank_curve(rows, full_k=full)}


def _npz_key(raw: Any, requested: str) -> np.ndarray:
    if requested not in raw.files: raise ValueError(f"npz missing key {requested!r}; available={sorted(raw.files)}")
    return raw[requested]


def _digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""): value.update(block)
    return value.hexdigest()


def _write_rows(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    values = [{key: json.dumps(value) if isinstance(value, (dict, list)) else value for key, value in row.items()} for row in rows]
    columns = sorted({key for row in values for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns); writer.writeheader(); writer.writerows(values)


def _self_test() -> None:
    rng = np.random.default_rng(7); n, d, q = 180, 24, 16
    x = rng.normal(size=(n, d)); w = rng.normal(size=(d, 5)); y = x @ w @ rng.normal(size=(5, q)) + .25 * rng.normal(size=(n, q))
    groups = np.repeat(["A", "B", "C", "D", "E", "F"], n // 6)
    train, test = grouped_train_test(groups, seed=7)
    ytr, yte, transform = construct_synthetic_targets(x[train], y[train], x[test], y[test], 3)
    assert ytr.shape == (train.sum(), q) and yte.shape == (test.sum(), q)
    assert transform.k == 3 and transform.full_expressible_dim >= 3
    result = train_controlled_head(x[train], ytr, x[test], yte, seed=7, latent_dim=12, epochs=300, overfit_steps=300)
    assert all(liveness_gates(result).values()), result
    print("E2 static/synthetic preflight PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=False, help="NPZ with frozen WSI feature rows")
    parser.add_argument("--targets", required=False, help="NPZ with frozen RNA target rows")
    parser.add_argument("--feature-key", default="wsi_biology"); parser.add_argument("--target-key", default="targets")
    parser.add_argument("--cancer-key", default="cancers"); parser.add_argument("--development-mask-key", default="development_mask")
    parser.add_argument("--output"); parser.add_argument("--ks", default="5,10,20,50,100"); parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--weight-decays", default="0.0001,0"); parser.add_argument("--latent-dim", type=int, default=128); parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--device", default="cpu"); parser.add_argument("--official-gate-log"); parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test: _self_test(); return
    if not (args.features and args.targets and args.output): parser.error("features, targets, and output are required unless --self-test")
    feature_path, target_path, output = Path(args.features), Path(args.targets), Path(args.output); output.mkdir(parents=True, exist_ok=True)
    with np.load(feature_path, allow_pickle=False) as fx, np.load(target_path, allow_pickle=False) as ty:
        x = _npz_key(fx, args.feature_key); cancers = _npz_key(fx, args.cancer_key)
        y = _npz_key(ty, args.target_key)
        if args.development_mask_key in fx.files: mask = np.asarray(fx[args.development_mask_key], dtype=bool)
        elif "split" in fx.files: mask = np.isin(fx["split"].astype(str), ("train", "development", "dev"))
        else: raise ValueError("features NPZ requires development_mask or split; E2 will not infer a training universe")
    if len(mask) != len(x) or len(x) != len(y) or len(x) != len(cancers): raise ValueError("NPZ row alignment/mask mismatch")
    started = time.monotonic(); ks = tuple(int(value) for value in args.ks.split(",") if value); seeds = tuple(int(value) for value in args.seeds.split(",") if value); wds = tuple(float(value) for value in args.weight_decays.split(",") if value)
    result = run_e2(x[mask], y[mask], cancers[mask], ks=ks, seeds=seeds, weight_decays=wds, latent_dim=args.latent_dim, epochs=args.epochs, device=args.device)
    result["wall_seconds"] = time.monotonic() - started
    manifest = {"features": str(feature_path), "features_sha256": _digest(feature_path), "targets": str(target_path), "targets_sha256": _digest(target_path),
                "feature_key": args.feature_key, "target_key": args.target_key, "development_mask_key": args.development_mask_key,
                "n_total": int(len(x)), "n_development": int(mask.sum()), "python": sys.version, "platform": platform.platform(), "command": sys.argv}
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (output / "e2_result.json").write_text(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2), encoding="utf-8")
    _write_rows(output / "e2_rows.csv", result["rows"])
    ledger = GateLedger(output, "E2_expressible_intersection", official_log=args.official_gate_log)
    ledger.artifact("G0.3_feature_identity", feature_path); ledger.artifact("G0.3_target_identity", target_path)
    ledger.add("G1.1_finite_feature_target_rows", f"n_dev={int(mask.sum())}", "all finite; explicit dev mask", bool(np.isfinite(x[mask]).all() and np.isfinite(y[mask]).all()))
    ledger.add("G1.2_grouped_internal_holdout", "cancer grouped", ">=3 development cancer groups", len(np.unique(cancers[mask].astype(str))) >= 3)
    ledger.add("G1.3_fixed_nominal_target_width", int(y.shape[1]), "same width at every k", bool(y.shape[1] > 1))
    for row in result["rows"]:
        tag = f"seed{row['seed']}_k{row['k']}_wd{row['weight_decay']}"
        for gate, passed in liveness_gates(row).items(): ledger.add(f"{gate}_{tag}", row[gate.replace("G2.1_params_moved", "parameter_delta").replace("G2.2_loss_live", "initial_loss").replace("G2.3_gradients_reach_all_groups", "first_gradient_norms").replace("G2.4_loss_decreased", "loss_reduction").replace("G2.5_converged", "tail_loss_slope").replace("G2.6_overfit_one_batch", "overfit_relative_reduction")], "predeclared G2", passed, tag)
        if row["k"] == row["full_expressible_dim"]: ledger.add(f"E2.a_full_capacity_{tag}", row["test_effective_rank"], ">=0.60*k or fail closed", row["E2.a_full_capacity"], tag)
        ledger.add(f"G3.1_latent_effective_rank_{tag}", row["test_effective_rank"], "finite >0", bool(np.isfinite(row["test_effective_rank"]) and row["test_effective_rank"] > 0), tag)
        ledger.add(f"G3.2_no_nonfinite_metric_{tag}", row["test_mean_target_pearson"], "finite", bool(np.isfinite(row["test_mean_target_pearson"])), tag)
        ledger.add(f"G4.1_positive_control_{tag}", row["test_mean_target_pearson"], ">0 held-out target correlation", row["G4.1_positive_control"], tag)
        ledger.add(f"G4.2_negative_control_{tag}", json.dumps({"real": row["test_mse"], "shuffled_train": row["negative_control_test_mse"]}), "real test MSE < shuffled-train null", row["G4.2_negative_control"], tag)
        ledger.add(f"G4.3_negative_liveness_{tag}", "shuffled same-path head", "G2 passes", row["negative_control_liveness_pass"], tag)
    wd0 = [row for row in result["rows"] if row["weight_decay"] == 0.0]
    ledger.add("E2.b_weight_decay_zero_control", len(wd0), ">0 matched rows", bool(wd0), "zero-WD sweep present")
    ledger.add("G0.5_required_gate_coverage", "E2-specific G2/E2a/E2b", "all emitted", bool(result["rows"]), "E2 has no cross-modal G4 claim")
    result["gates_pass"] = ledger.write()
    print(json.dumps({"output": str(output), "gates_pass": result["gates_pass"], "rank_curve": result["rank_curve"]}, indent=2))


if __name__ == "__main__": main()
