"""Self-tests for the attributable-basis rotations and the channel-share statistic.

The two results these support are both *comparisons*, so the thing that has to be
proved is not that a number comes out — it is that the arms differ only in the
one way they are supposed to. Concretely: that a rotation is a rotation (the span
does not move), that the ridge really is linear in the target block (which is what
makes an R²-optimising rotation computable without a second copy of the ridge and
what makes the invariance argument true rather than plausible), and that
residualising a block and then selecting columns is the same as residualising the
columns — because every channel-share number is read off one residualisation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.attributable_basis import (ROTATIONS, cross_line_rotation, ica_rotation,
                                            r2_optimising_rotation, varimax_rotation)
from morpheus.v2.calibra.residualise import confound_design, cross_fitted_residuals
from morpheus.v2.causal_attribution import (certifiable_attribution, gene_fold_ridge_r2,
                                            validated_rotation)
from morpheus.v2.channel_share import _decile_strata, _stratified_draw, subset_channel
from morpheus.v2.perturbation_basis_common import atom_folds, subspace_alignment


def _orthonormal(n_rows: int, width: int, seed: int) -> np.ndarray:
    return np.linalg.qr(np.random.default_rng(seed).normal(size=(n_rows, width)))[0]


# --- the atom split both the fitter and the scorer must agree on ---------


def test_atom_folds_are_disjoint_exhaustive_and_deterministic():
    a, b = atom_folds(101, seed=3)
    assert len(np.intersect1d(a, b)) == 0
    assert np.array_equal(np.union1d(a, b), np.arange(101))
    assert len(a) == 50 and len(b) == 51
    assert np.array_equal(a, atom_folds(101, seed=3)[0])
    assert not np.array_equal(a, atom_folds(101, seed=4)[0])


# --- the ridge is linear in the target block ----------------------------


def test_return_prediction_is_off_by_default_and_changes_nothing():
    """The flag must be inert: the published 29/128 was produced without it."""
    rng = np.random.default_rng(0)
    design = rng.normal(size=(60, 40))
    gram = design @ design.T
    targets = rng.normal(size=(60, 4))
    plain = gene_fold_ridge_r2(gram, targets, seed=1)
    withpred = gene_fold_ridge_r2(gram, targets, seed=1, return_prediction=True)
    assert "prediction" not in plain
    assert set(withpred) - set(plain) == {"prediction"}
    for key in plain:
        if isinstance(plain[key], np.ndarray):
            assert np.array_equal(plain[key], withpred[key]), key
        else:
            assert plain[key] == withpred[key], key


def test_rotated_r2_is_recoverable_from_the_out_of_fold_residual():
    """R²_k(T@R) = 1 - (RᵀAR)_kk/(RᵀBR)_kk with A = EᵀE, B = T_cᵀT_c.

    This identity is the entire basis of `r2_optimising_rotation`; if it were
    false the arm would be optimising a quantity unrelated to the one scored.
    """
    rng = np.random.default_rng(5)
    design = rng.normal(size=(80, 50))
    gram = design @ design.T
    targets = design @ rng.normal(size=(50, 6)) + 0.4 * rng.normal(size=(80, 6))
    fit = gene_fold_ridge_r2(gram, targets, alphas=(1e2,), seed=2, return_prediction=True)
    residual = targets - fit["prediction"]
    centred = targets - targets.mean(axis=0, keepdims=True)
    a_matrix, b_matrix = residual.T @ residual, centred.T @ centred

    rotation = np.linalg.qr(rng.normal(size=(6, 6)))[0]
    closed_form = 1.0 - (np.diag(rotation.T @ a_matrix @ rotation)
                         / np.diag(rotation.T @ b_matrix @ rotation))
    measured = gene_fold_ridge_r2(gram, targets @ rotation, alphas=(1e2,), seed=2)["r2"]
    assert np.allclose(measured, closed_form, atol=1e-10)


def test_mean_r2_is_rotation_invariant_when_the_target_columns_are_equal_norm():
    """The structural fact predeclared before the run, checked rather than argued.

    With `B = sigma^2 I` the sum of `diag(R^T A R)` is `trace(A)` for every
    orthogonal R, so the MEAN R2 cannot be raised by any rotation -- only
    redistributed across axes. Planted exactly, then perturbed to show the
    invariance degrades gracefully rather than being an artefact of the plant.
    """
    rng = np.random.default_rng(7)
    basis = _orthonormal(200, 8, 11)                    # exactly equal-norm columns
    residual = rng.normal(size=(200, 8)) * 0.3
    a_matrix = residual.T @ residual
    b_matrix = basis.T @ basis                          # = I

    def mean_r2(rotation):
        a = np.diag(rotation.T @ a_matrix @ rotation)
        b = np.diag(rotation.T @ b_matrix @ rotation)
        return float(np.mean(1.0 - a / b))

    identity = mean_r2(np.eye(8))
    for seed in range(5):
        rotation = np.linalg.qr(rng.normal(size=(8, 8)))[0]
        assert mean_r2(rotation) == pytest.approx(identity, abs=1e-10), seed
    # ... and the per-axis spread DOES move, which is the only thing a rotation
    # can do to this statistic.
    rotation = np.linalg.qr(rng.normal(size=(8, 8)))[0]
    spread = np.diag(rotation.T @ a_matrix @ rotation).std()
    assert spread != pytest.approx(np.diag(a_matrix).std(), rel=1e-3)


# --- the rotations themselves -------------------------------------------


@pytest.mark.parametrize("name", ["varimax", "ica"])
def test_honest_rotations_are_orthogonal_and_do_not_move_the_span(name):
    loadings = _orthonormal(300, 6, 21)
    fitted = varimax_rotation(loadings) if name == "varimax" else ica_rotation(loadings, seed=0)
    rotation = fitted["rotation"]
    assert np.abs(rotation.T @ rotation - np.eye(6)).max() < 1e-9
    overlap = subspace_alignment(loadings, loadings @ rotation)["mean_squared_cosine"]
    assert overlap == pytest.approx(1.0, abs=1e-9)


def test_varimax_finds_a_sparse_basis_that_is_hidden_from_the_identity():
    """Plant a sparse basis, rotate it away, and check varimax rotates it back."""
    rng = np.random.default_rng(31)
    width, n_genes = 5, 400
    sparse = np.zeros((n_genes, width))
    for k in range(width):
        rows = rng.choice(n_genes, size=12, replace=False)
        sparse[rows, k] = rng.normal(size=12)
    sparse = np.linalg.qr(sparse)[0]
    hidden = sparse @ np.linalg.qr(rng.normal(size=(width, width)))[0]

    def kurtosis(matrix):
        return float(((matrix ** 2).sum(axis=0) ** 2).sum() - (matrix ** 4).sum())

    recovered = hidden @ varimax_rotation(hidden)["rotation"]
    assert kurtosis(recovered) < kurtosis(hidden), "varimax did not increase simple structure"
    # the recovered axes match the planted ones up to sign and order
    overlap = np.abs(recovered.T @ sparse)
    assert overlap.max(axis=1).min() > 0.9, overlap.max(axis=1)


def test_r2_optimising_rotation_never_decreases_its_own_objective():
    rng = np.random.default_rng(41)
    targets = _orthonormal(150, 6, 43)
    residual = targets * rng.uniform(0.1, 0.9, size=(1, 6)) + 0.05 * rng.normal(size=(150, 6))
    fitted = r2_optimising_rotation(residual, targets, max_iter=60)
    assert fitted["objective"] >= fitted["objective_at_identity"] - 1e-12
    assert np.abs(fitted["rotation"].T @ fitted["rotation"] - np.eye(6)).max() < 1e-9
    assert len(fitted["objective_by_start"]) == 3


def _two_line_profiles(weights, n_atoms, seed):
    """Two ``[atom, axis]`` profiles whose per-axis cross-line correlation is ``w²``."""
    rng = np.random.default_rng(seed)
    weights = np.asarray(weights, dtype=float)
    shared = rng.normal(size=(n_atoms, len(weights)))
    residual = np.sqrt(np.maximum(1.0 - weights ** 2, 0.0))
    return (shared * weights + rng.normal(size=shared.shape) * residual,
            shared * weights + rng.normal(size=shared.shape) * residual)


def test_mean_cross_line_agreement_is_a_budget_a_rotation_can_only_redistribute():
    """The conservation the `xline` design turns on, measured rather than argued.

    With unit-variance profiles the per-axis cross-line correlation is
    ``r_kᵀ M r_k`` and its sum over an orthonormal basis is ``tr(M)``, so no
    rotation can raise the MEAN. This is why `xline_mean` is expected to be inert
    and `xline` maximises a count instead.
    """
    a, b = _two_line_profiles([1.0, 0.7, 0.5, 0.2, 0.0], 4000, 51)
    fitted = cross_line_rotation(a, b, np.arange(len(a)), objective="mean", max_iter=200)
    assert fitted["mean_correlation_at_best"] == pytest.approx(
        fitted["mean_correlation_at_identity"], abs=0.02)


def test_the_count_objective_spreads_a_concentrated_agreement_over_more_axes():
    """Budget ~1.5 concentrated on two axes; the bar is 0.30, so ~5 axes can clear it."""
    a, b = _two_line_profiles([1.0, np.sqrt(0.5), 0.0, 0.0, 0.0], 6000, 53)
    rows = np.arange(len(a))
    fitted = cross_line_rotation(a, b, rows, objective="count", floor=0.30, max_iter=400)
    assert fitted["n_above_floor_at_identity"] == 2, fitted["n_above_floor_at_identity"]
    assert fitted["n_above_floor_at_best"] > fitted["n_above_floor_at_identity"], fitted
    assert fitted["objective"] > fitted["objective_at_identity"]
    # ... and it did it by redistributing, not by creating agreement.
    assert fitted["mean_correlation_at_best"] <= fitted["mean_correlation_at_identity"] + 0.03


def test_the_registry_names_every_arm_that_was_predeclared():
    assert set(ROTATIONS) == {"none", "varimax", "ica", "r2opt", "xline_mean", "xline"}


# --- the guard that makes a basis-choice test a basis-choice test --------


def test_a_non_orthogonal_or_span_moving_matrix_is_refused(tmp_path):
    loadings = _orthonormal(50, 4, 61)
    good = np.linalg.qr(np.random.default_rng(0).normal(size=(4, 4)))[0]
    np.savez(tmp_path / "good.npz", rotation=good)
    matrix, summary = validated_rotation(str(tmp_path / "good.npz"), loadings)
    assert summary["mean_squared_cosine_vs_pca_span"] == pytest.approx(1.0, abs=1e-10)
    assert np.array_equal(matrix, good)

    np.savez(tmp_path / "scaled.npz", rotation=good @ np.diag([1.0, 1.0, 1.0, 2.0]))
    with pytest.raises(ValueError, match="not orthogonal"):
        validated_rotation(str(tmp_path / "scaled.npz"), loadings)

    np.savez(tmp_path / "wide.npz", rotation=np.eye(4)[:, :3])
    with pytest.raises(ValueError, match="rotation must be"):
        validated_rotation(str(tmp_path / "wide.npz"), loadings)


def test_a_rank_deficient_rotation_that_moves_the_span_is_refused(tmp_path):
    """Orthogonal to 1e-9 but not to 1e-8: the guard's bar, not a hand-wave."""
    loadings = _orthonormal(50, 4, 62)
    nearly = np.eye(4)
    nearly[0, 1] = 1e-6                                   # breaks RᵀR = I at 1e-6
    np.savez(tmp_path / "nearly.npz", rotation=nearly)
    with pytest.raises(ValueError, match="not orthogonal"):
        validated_rotation(str(tmp_path / "nearly.npz"), loadings)


# --- the certificate, re-scored on a substituted cross-line column -------


def test_certificate_thresholds_are_unchanged_by_the_column_substitution():
    """Swapping in the held-out-atom column must go through the same function.

    The held-out count is only comparable with the headline count if it applies
    the identical four bars, so it is computed by calling
    `certifiable_attribution` on a substituted frame rather than by restating the
    conditions.
    """
    frame = pd.DataFrame({
        "r2_cv": [0.5, 0.5, 0.5],
        "r2_cv_random_direction_null": [0.01, 0.01, 0.01],
        "shuffle_rank_spearman": [0.1, 0.1, 0.9],
        "cross_line_rank_spearman": [0.9, 0.1, 0.9],
        "cross_line_rank_spearman_atom_fold_b": [0.1, 0.9, 0.9],
        "attributed_set_coherence_percentile": [0.99, 0.99, 0.99]})
    assert certifiable_attribution(frame).tolist() == [True, False, False]
    swapped = frame.assign(cross_line_rank_spearman=frame["cross_line_rank_spearman_atom_fold_b"])
    assert certifiable_attribution(swapped).tolist() == [False, True, False]


# --- the channel-share statistic ----------------------------------------


def test_residualise_then_subset_equals_subset_then_residualise():
    """Every channel-share number is read off one residualisation of the block."""
    rng = np.random.default_rng(71)
    frame = pd.DataFrame({"cancer": rng.choice(list("abcd"), size=300),
                          "tss": rng.choice(list("xyz"), size=300)})
    design = confound_design(frame, ["cancer", "tss"])
    block = rng.normal(size=(300, 12)) + design[:, :1] * 2.0
    columns = np.asarray([0, 3, 5, 9])
    whole = cross_fitted_residuals(block, design, seed=42)
    part = cross_fitted_residuals(block[:, columns], design, seed=42)
    assert np.abs(whole[:, columns] - part).max() < 1e-9


def test_subset_channel_sees_the_columns_that_carry_the_signal():
    """A subset containing the planted direction must beat one that does not."""
    rng = np.random.default_rng(81)
    n = 600
    latent = rng.normal(size=n)
    x = np.column_stack([latent + 0.5 * rng.normal(size=n), rng.normal(size=(n, 5))])
    y = np.column_stack([rng.normal(size=(n, 4)), latent + 0.5 * rng.normal(size=n),
                         rng.normal(size=(n, 3))])
    with_signal = subset_channel(x, y, [3, 4, 5], n_components=3, seed=42)
    without = subset_channel(x, y, [0, 1, 2], n_components=3, seed=42)
    assert with_signal > 0.5 > without, (with_signal, without)


def test_stratified_draw_matches_the_target_variance_profile():
    variance = np.sort(np.random.default_rng(91).random(120))[::-1]
    strata = _decile_strata(variance)
    assert len(np.unique(strata)) == 10
    target = np.arange(90, 119)                     # all from the low-variance tail
    draw = _stratified_draw(np.random.default_rng(0), strata, target)
    assert len(draw) == len(target)
    for stratum in np.unique(strata):
        assert (strata[draw] == stratum).sum() == (strata[target] == stratum).sum()
