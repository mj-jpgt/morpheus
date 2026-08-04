"""Self-tests for the post-PBS constructions and for the attribution statistic.

Each test plants the situation the primitive exists to detect and asserts it
reacts. The point of the battery is that a *negative* result from these
constructions has to be readable as "the construction did not help" rather than
"the code did not work", and that is only true if each primitive is shown to move
when it should.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.build_causal_basis_targets import (CONSTRUCTIONS, _row_correlation,
                                                    cross_line_consensus_atoms)
from morpheus.v2.causal_attribution import atom_cosines, gene_fold_ridge_r2
from morpheus.v2.perturbation_basis_common import SCALINGS, load_aligned_response, subspace_alignment


# --- subspace alignment: the invariance guard ----------------------------


def test_alignment_is_one_for_a_rotation_and_zero_for_an_orthogonal_span():
    """Held-out top-CCA cannot see a rotation, so the guard must not either."""
    rng = np.random.default_rng(0)
    basis = np.linalg.qr(rng.normal(size=(40, 6)))[0]
    rotation = np.linalg.qr(rng.normal(size=(6, 6)))[0]
    assert subspace_alignment(basis, basis @ rotation)["mean_squared_cosine"] == pytest.approx(1.0, abs=1e-10)
    full = np.linalg.qr(rng.normal(size=(40, 12)))[0]
    assert subspace_alignment(full[:, :6], full[:, 6:])["mean_squared_cosine"] < 1e-16


def test_alignment_is_intermediate_for_a_partial_overlap():
    rng = np.random.default_rng(1)
    full = np.linalg.qr(rng.normal(size=(40, 8)))[0]
    overlap = subspace_alignment(full[:, :4], full[:, 2:6])["mean_squared_cosine"]
    assert 0.4 < overlap < 0.6, overlap


# --- cross-cell-line consensus -------------------------------------------


def test_row_correlation_agrees_with_numpy_row_by_row():
    rng = np.random.default_rng(2)
    a, b = rng.normal(size=(7, 30)), rng.normal(size=(7, 30))
    expected = np.asarray([np.corrcoef(a[i], b[i])[0, 1] for i in range(len(a))])
    assert np.allclose(_row_correlation(a, b), expected)


def test_consensus_retains_reproducible_atoms_and_drops_noise_atoms():
    """Atoms 0-19 share a signal across the two lines; atoms 20-99 do not."""
    rng = np.random.default_rng(3)
    shared = rng.normal(size=(20, 200))
    primary = np.vstack([shared + 0.3 * rng.normal(size=(20, 200)), rng.normal(size=(80, 200))])
    secondary = np.vstack([shared + 0.3 * rng.normal(size=(20, 200)), rng.normal(size=(80, 200))])
    result = cross_line_consensus_atoms(primary, secondary, percentile=95.0, n_null=5000, seed=3)
    assert result["retained"][:20].all(), "reproducible atoms were dropped"
    assert result["retained"][20:].mean() < 0.10, result["retained"][20:].mean()


def test_consensus_threshold_is_graded_against_a_mismatched_null_not_against_zero():
    """Two matrices sharing global structure must not all pass on a bare r > 0 bar."""
    rng = np.random.default_rng(4)
    background = rng.normal(size=(1, 200))
    primary = background + 0.6 * rng.normal(size=(120, 200))
    secondary = background + 0.6 * rng.normal(size=(120, 200))
    result = cross_line_consensus_atoms(primary, secondary, percentile=95.0, n_null=5000, seed=4)
    assert result["null_median"] > 0.05, "the shared-structure null collapsed; the test is not testing anything"
    assert result["n_retained"] < 0.25 * result["n_atoms"], result["n_retained"]


def test_consensus_refuses_misaligned_matrices():
    with pytest.raises(ValueError):
        cross_line_consensus_atoms(np.zeros((4, 5)), np.zeros((4, 6)))


# --- the attribution statistic -------------------------------------------


def _design(n_genes: int, n_atoms: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).normal(size=(n_genes, n_atoms))


def test_gene_fold_ridge_recovers_a_direction_inside_the_atom_span_and_not_one_outside():
    """The statistic's reason to exist: in-span is reconstructable, out-of-span is not."""
    rng = np.random.default_rng(5)
    design = _design(300, 40, 5)
    inside = design @ rng.normal(size=(40, 1))
    projector = design @ np.linalg.pinv(design)
    outside = rng.normal(size=(300, 1))
    outside = outside - projector @ outside
    result = gene_fold_ridge_r2(design @ design.T, np.hstack([inside, outside]), n_folds=5, seed=5)
    assert result["r2"][0] > 0.9, result["r2"][0]
    assert result["r2"][1] < 0.1, result["r2"][1]


def test_gene_fold_ridge_is_cross_validated_not_in_sample():
    """With more atoms than genes the in-sample fit is exactly 1; the CV fit is not."""
    rng = np.random.default_rng(6)
    design = _design(120, 400, 6)
    target = rng.normal(size=(120, 1))
    result = gene_fold_ridge_r2(design @ design.T, target, n_folds=5, seed=6)
    assert result["r2"][0] < 0.5, f"a pure-noise target scored {result['r2'][0]} out of fold"


def test_gene_fold_ridge_returns_the_full_alpha_curve_and_selects_from_it():
    design = _design(200, 30, 7)
    target = design @ np.random.default_rng(7).normal(size=(30, 2))
    result = gene_fold_ridge_r2(design @ design.T, target, seed=7)
    assert result["r2_curve"].shape == (len(result["alphas"]), 2)
    assert result["selected_alpha"] in result["alphas"]
    assert np.allclose(result["r2"], result["r2_curve"][result["selected_index"]])


def test_gene_fold_ridge_refuses_a_non_square_gram():
    with pytest.raises(ValueError):
        gene_fold_ridge_r2(np.zeros((4, 5)), np.zeros((4, 2)))


def test_atom_cosines_find_the_planted_atom():
    rng = np.random.default_rng(8)
    response = rng.normal(size=(50, 120))
    direction = (response[7] / np.linalg.norm(response[7]))[:, None]
    cosine = atom_cosines(response, direction)
    assert int(np.argmax(np.abs(cosine[:, 0]))) == 7
    assert cosine[7, 0] == pytest.approx(1.0, abs=1e-10)


def test_atom_cosines_are_invariant_to_atom_and_direction_scale():
    rng = np.random.default_rng(9)
    response, direction = rng.normal(size=(12, 40)), rng.normal(size=(40, 3))
    scaled = atom_cosines(response * 7.0, direction * -1.0)
    assert np.allclose(scaled, -atom_cosines(response, direction))


# --- declared names ------------------------------------------------------


def test_scaling_and_construction_names_are_closed_sets():
    assert set(SCALINGS) == {"tcga_sd", "own_sd", "raw"}
    assert set(CONSTRUCTIONS) == {"pbs_rebuild", "joint_cca", "consensus", "domain_adapted"}


def test_load_aligned_response_refuses_an_undeclared_scaling(tmp_path):
    with pytest.raises(ValueError, match="unknown perturbation scaling"):
        load_aligned_response(tmp_path / "absent.h5ad", ["A"], scaling="whatever")
