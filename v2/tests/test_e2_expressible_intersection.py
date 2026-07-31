"""Focused contract tests for E2's synthetic expressible-intersection sweep."""
from __future__ import annotations

import numpy as np

try:  # Lambda runs this checkout through a ``morpheus`` workspace symlink.
    from morpheus.v2.calibra.e2_expressible_intersection import (
        classify_rank_curve, construct_synthetic_targets, fit_synthetic_target_transform,
        grouped_train_test, liveness_gates, train_controlled_head,
    )
except ModuleNotFoundError:  # Local checkout direct test path.
    from v2.calibra.e2_expressible_intersection import (
        classify_rank_curve, construct_synthetic_targets, fit_synthetic_target_transform,
        grouped_train_test, liveness_gates, train_controlled_head,
    )


def _data(seed=3, n=180, features=20, targets=14):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, features))
    # Five deliberately image-expressible molecular directions plus residual
    # target variation makes the expected test construction non-degenerate.
    x_to_shared = rng.normal(size=(features, 5))
    shared_to_y = rng.normal(size=(5, targets))
    y = x @ x_to_shared @ shared_to_y + .2 * rng.normal(size=(n, targets))
    cancer = np.repeat(np.array(["A", "B", "C", "D", "E", "F"]), n // 6)
    return x, y, cancer


def test_grouped_split_never_leaks_cancer_identity():
    _, _, cancers = _data()
    train, test = grouped_train_test(cancers, seed=12)
    assert not set(cancers[train]).intersection(set(cancers[test]))
    assert train.any() and test.any()


def test_synthetic_targets_are_fixed_width_and_train_only_fit():
    x, y, cancers = _data()
    train, test = grouped_train_test(cancers, seed=7)
    train_target, test_target, transform = construct_synthetic_targets(x[train], y[train], x[test], y[test], 3)
    assert train_target.shape == (train.sum(), y.shape[1])
    assert test_target.shape == (test.sum(), y.shape[1])
    assert transform.k == 3
    # Altering held-out labels cannot change train-fold transform parameters.
    changed = y.copy(); changed[test] += 1_000
    refit = fit_synthetic_target_transform(x[train], changed[train], 3)
    assert np.allclose(transform.x_basis, refit.x_basis)
    assert np.allclose(transform.y_basis, refit.y_basis)
    assert np.allclose(transform.residual_beta, refit.residual_beta)


def test_invalid_k_fails_closed():
    x, y, cancers = _data(); train, _ = grouped_train_test(cancers, seed=1)
    full = fit_synthetic_target_transform(x[train], y[train], 1).full_expressible_dim
    try:
        fit_synthetic_target_transform(x[train], y[train], full + 1)
    except ValueError as error:
        assert "must lie" in str(error)
    else:
        raise AssertionError("out-of-range k must fail rather than silently clamp")


def test_head_liveness_and_one_batch_overfit_on_realistic_synthetic_data():
    x, y, cancers = _data(); train, test = grouped_train_test(cancers, seed=3)
    ytr, yte, _ = construct_synthetic_targets(x[train], y[train], x[test], y[test], 3)
    metrics = train_controlled_head(x[train], ytr, x[test], yte, seed=3, latent_dim=12, epochs=300, overfit_steps=400)
    assert all(liveness_gates(metrics).values()), metrics
    assert metrics["test_effective_rank"] > 1.0


def test_rank_curve_requires_a_k_association_for_h2():
    rows = []
    for k, rank in ((5, 4), (10, 8), (20, 15), (40, 28)):
        rows.append({"k": k, "test_effective_rank": rank, "weight_decay": 1e-4,
                     "reference_weight_decay": 1e-4, "all_liveness_pass": True})
    verdict = classify_rank_curve(rows, full_k=40)
    assert verdict["verdict"] == "supports_expressible_intersection"
    # A flat high-rank curve does not accidentally become an H2 result.
    for row in rows: row["test_effective_rank"] = 30
    assert classify_rank_curve(rows, full_k=40)["verdict"] != "supports_expressible_intersection"
