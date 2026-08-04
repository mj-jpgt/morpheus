"""Tests for the Track 1 negative-control battery.

These are the battery's own self-tests. A must-FAIL control is only evidence if
it *can* fail when it should and *does* pass when it should, so each test below
plants the situation the control exists to catch and asserts the control reacts.
A control that cannot be shown to react is decoration.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from morpheus.v2.calibra.confound_certificate import (balanced_accuracy, certify_axes,
                                                      lda_oof_balanced_accuracy,
                                                      nearest_class_mean_oof,
                                                      within_stratum_permutations)
from morpheus.v2.calibra.known_covariate_control import (_auroc, evaluate_known_covariate,
                                                         out_of_fold_scores, within_cancer_auroc)
from morpheus.v2.calibra.run_calibra import (grade_random_controls, random_direction_column_correlation,
                                             score_target_block_per_column)
from morpheus.v2.calibra.spectral import heldout_single_direction_correlation


def _patients(n: int, n_sites: int = 6, seed: int = 0):
    """Synthetic TCGA-shaped barcodes whose TSS field is a real, recoverable site."""
    rng = np.random.default_rng(seed)
    site = rng.integers(0, n_sites, size=n)
    ids = np.asarray([f"TCGA-{10 + s:02d}-{i:04d}" for i, s in enumerate(site)])
    cancers = np.asarray(["BRCA" if i % 2 else "LUAD" for i in range(n)])
    return ids, cancers, site


# --- T1.3 confound certificate -------------------------------------------


def test_balanced_accuracy_is_chance_for_a_constant_prediction():
    true = np.repeat(np.arange(4), 25)
    assert abs(balanced_accuracy(true, np.zeros(100, dtype=np.int64), 4) - 0.25) < 1e-12
    assert balanced_accuracy(true, true, 4) == 1.0


def test_nearest_class_mean_recovers_a_planted_site_axis_and_ignores_noise():
    """THE certificate's self-test: axis 0 IS the site code, axis 1 is noise."""
    rng = np.random.default_rng(1)
    n, n_sites = 600, 5
    site = rng.integers(0, n_sites, size=n)
    features = np.column_stack([site * 4.0 + rng.normal(scale=0.3, size=n), rng.normal(size=n)])
    accuracy = nearest_class_mean_oof(features, site, n_sites, seed=1)
    assert accuracy[0] > 0.9, f"planted site axis was not detected: {accuracy[0]}"
    assert accuracy[1] < 0.35, f"pure noise axis scored {accuracy[1]}, well above chance {1 / n_sites}"


def test_joint_lda_sees_site_spread_across_axes_that_no_single_axis_shows():
    """The reason the joint test is not optional: site can hide in a combination."""
    rng = np.random.default_rng(2)
    n, n_sites, d = 800, 4, 6
    site = rng.integers(0, n_sites, size=n)
    # A site code written into a random direction, so every SINGLE coordinate is
    # a weak, noisy view of it while the combination is not.
    direction = rng.normal(size=(n_sites, d))
    features = direction[site] * 1.2 + rng.normal(size=(n, d)) * 2.0
    per_axis = nearest_class_mean_oof(features, site, n_sites, seed=2)
    joint = lda_oof_balanced_accuracy(features, site, n_sites, seed=2)
    assert joint > per_axis.max() + 0.10, (joint, per_axis.max())


def test_within_stratum_permutation_preserves_the_cancer_site_association():
    ids, cancers, site = _patients(400, seed=3)
    permuted = within_stratum_permutations(site, cancers, n_permutations=3, seed=3)
    for draw in permuted:
        assert sorted(draw.tolist()) == sorted(site.tolist())
        for cancer in np.unique(cancers):
            rows = cancers == cancer
            # the multiset of labels inside each cancer is untouched
            assert sorted(draw[rows].tolist()) == sorted(site[rows].tolist())


def test_certificate_refuses_a_state_that_is_a_site_code_and_clears_pure_noise():
    """Both directions. A site code must NOT be certified; noise must be."""
    ids, cancers, site = _patients(500, n_sites=5, seed=4)
    site_from_id = np.asarray([int(str(i).split("-")[1]) for i in ids])
    rng = np.random.default_rng(4)
    leaking = np.column_stack([site_from_id * 5.0 + rng.normal(scale=0.2, size=len(ids)),
                               rng.normal(size=len(ids))])
    clean = rng.normal(size=(len(ids), 2))
    bad = certify_axes(leaking, ids, cancers, min_site_count=5, n_permutations=60, n_boot=25,
                       n_boot_axes=2, seed=4)
    good = certify_axes(clean, ids, cancers, min_site_count=5, n_permutations=60, n_boot=25,
                        n_boot_axes=2, seed=4)
    assert bad["certified"] is False and bad["n_breaching_axes"] >= 1
    assert bad["joint_certified"] is False
    assert good["certified"] is True, good["breaching_axes"]


def test_certificate_residualisation_removes_a_pure_site_code():
    """The adjustment must discharge the confound it claims to adjust for."""
    ids, cancers, _ = _patients(500, n_sites=5, seed=5)
    site_from_id = np.asarray([int(str(i).split("-")[1]) for i in ids])
    rng = np.random.default_rng(5)
    leaking = np.column_stack([site_from_id * 5.0 + rng.normal(scale=0.2, size=len(ids))] * 2)
    raw = certify_axes(leaking, ids, cancers, min_site_count=5, n_permutations=60, n_boot=20,
                       n_boot_axes=1, seed=5)
    adjusted = certify_axes(leaking, ids, cancers, min_site_count=5, n_permutations=60, n_boot=20,
                            n_boot_axes=1, seed=5, residualise=True)
    assert raw["joint_lda_balanced_accuracy"] > adjusted["joint_lda_balanced_accuracy"] + 0.2


# --- T1.4 random gene-set controls ---------------------------------------


def test_block_scoring_matches_the_single_column_reference():
    """The vectorised block path must be the per-column path, exactly."""
    rng = np.random.default_rng(6)
    x = rng.normal(size=(300, 8))
    y = np.column_stack([x @ rng.normal(size=8) + rng.normal(size=300), rng.normal(size=300)])
    block = score_target_block_per_column(x, y, seed=6)
    single = [heldout_single_direction_correlation(x, y[:, j], seed=6) for j in range(y.shape[1])]
    np.testing.assert_allclose(block, single, rtol=1e-10, atol=1e-12)


def test_fitted_direction_separates_signal_from_noise_targets():
    rng = np.random.default_rng(7)
    x = rng.normal(size=(400, 10))
    signal = x @ rng.normal(size=10) + rng.normal(scale=0.5, size=400)
    values = score_target_block_per_column(x, np.column_stack([signal, rng.normal(size=400)]), seed=7)
    assert values[0] > 0.6 and abs(values[1]) < 0.2, values


def test_random_direction_statistic_is_null_like_even_for_a_readable_target():
    """The scale gap, as a test. A target a FITTED direction reads at ~0.9 is read
    far lower by a RANDOM direction, so grading one against a floor measured in the
    other is a comparison that has to be argued for, not assumed."""
    rng = np.random.default_rng(8)
    x = rng.normal(size=(500, 20))
    y = (x @ rng.normal(size=20))[:, None]
    fitted = score_target_block_per_column(x, y, seed=8)[0]
    random_u = random_direction_column_correlation(x, y, n_draws=20, seed=8)[0]
    assert fitted > 0.9, fitted
    assert abs(random_u) < 0.25, random_u


def test_random_direction_statistic_is_a_magnitude_not_a_signed_median():
    """The sign fix (schema 2), as a regression test.

    ``random_direction_column_correlation`` is graded against a strictly positive
    ``detection_floor``. With ``u`` drawn at random the sign of each draw is random,
    so a SIGNED median collapses towards zero for every column whatever it carries
    -- which made the negative control pass structurally rather than on merit. The
    statistic must be non-negative, and on a column a random direction genuinely
    does see, it must NOT sit at zero."""
    rng = np.random.default_rng(4)
    x = rng.normal(size=(400, 3))            # low dimension: a random u overlaps a lot
    y = (x @ np.array([1.0, 0.0, 0.0]))[:, None]
    values = random_direction_column_correlation(x, y, n_draws=40, seed=4)
    assert np.all(values >= 0.0), values
    assert values[0] > 0.3, values          # a signed median would sit near 0 here
    noise = random_direction_column_correlation(x, rng.normal(size=(400, 5)),
                                                n_draws=40, seed=4)
    assert np.all(noise >= 0.0) and float(np.median(noise)) < 0.1, noise


def test_grade_random_controls_fails_when_too_many_controls_clear_the_floor():
    names = np.asarray([f"RANDOM_CONTROL__t__{i:02d}" for i in range(20)])
    below = grade_random_controls(np.full(20, 0.05), names, detection_floor=0.2)
    assert below["passed"] is True and below["n_exceedances"] == 0
    breached = np.full(20, 0.05)
    breached[:3] = 0.9                       # 15% exceedance, over the 5% ceiling
    verdict = grade_random_controls(breached, names, detection_floor=0.2)
    assert verdict["passed"] is False
    assert verdict["n_exceedances"] == 3 and verdict["exceedance_fraction"] > 0.05
    # the identity of every exceedance is reported, never just the count
    assert [name for name, _ in verdict["exceedances"]] == names[:3].tolist()


def test_grade_random_controls_fails_on_a_high_median_even_within_the_ceiling():
    names = np.asarray([f"RANDOM_CONTROL__t__{i:02d}" for i in range(20)])
    verdict = grade_random_controls(np.full(20, 0.19), names, detection_floor=0.2)
    assert verdict["median_below_floor"] is True and verdict["passed"] is True
    verdict = grade_random_controls(np.full(20, 0.201), names, detection_floor=0.2)
    assert verdict["median_below_floor"] is False and verdict["passed"] is False


# --- T1.7(b) known-covariate positive control -----------------------------


def test_auroc_matches_a_known_value_and_handles_ties():
    assert abs(_auroc(np.array([0, 0, 1, 1]), np.array([1.0, 2.0, 3.0, 4.0])) - 1.0) < 1e-12
    assert abs(_auroc(np.array([0, 0, 1, 1]), np.array([4.0, 3.0, 2.0, 1.0]))) < 1e-12
    assert abs(_auroc(np.array([0, 0, 1, 1]), np.ones(4)) - 0.5) < 1e-12


def test_within_cancer_auroc_ignores_lineage_that_pooled_auroc_would_reward():
    """The reason within-cancer is the primary statistic: a feature that only
    encodes cancer type scores high pooled and exactly chance within cancer."""
    n = 600
    cancers = np.asarray(["A"] * (n // 2) + ["B"] * (n // 2))
    rng = np.random.default_rng(9)
    # prevalence differs sharply by cancer; the feature IS the cancer indicator
    labels = np.concatenate([rng.random(n // 2) < 0.9, rng.random(n // 2) < 0.1])
    score = (cancers == "A").astype(float) + rng.normal(scale=1e-6, size=n)
    pooled = _auroc(labels, score)
    within, _ = within_cancer_auroc(labels, score, cancers)
    assert pooled > 0.85, pooled
    assert abs(within - 0.5) < 0.05, within


def test_known_covariate_control_fails_an_expectation_it_cannot_meet():
    """'We recovered something' is not a pass, in both directions."""
    rng = np.random.default_rng(10)
    n = 800
    cancers = np.asarray(["A", "B"] * (n // 2))
    ids = np.asarray([f"TCGA-{10 + (i % 4):02d}-{i:04d}" for i in range(n)])
    labels = rng.random(n) < 0.4
    features = np.column_stack([labels * 2.0 + rng.normal(size=n), rng.normal(size=(n, 4)).T]).T \
        if False else np.column_stack([labels * 2.0 + rng.normal(size=n), rng.normal(size=(n, 4))])
    strong = evaluate_known_covariate(features, labels, cancers, ids, expected_low=0.60,
                                      expected_high=0.80, n_boot=80, n_permutations=80, seed=10)
    assert strong["status"] == "scored"
    # A covariate recovered FAR above its published strength is a leak, not a win.
    impossible = evaluate_known_covariate(features, labels, cancers, ids, expected_low=0.10,
                                          expected_high=0.20, n_boot=80, n_permutations=80, seed=10)
    assert impossible["passed"] is False and impossible["above_expected"] is True
    # A covariate nothing can read must not pass an ordinary expectation.
    noise = evaluate_known_covariate(rng.normal(size=(n, 5)), labels, cancers, ids,
                                     expected_low=0.60, expected_high=0.80, n_boot=80,
                                     n_permutations=80, seed=10)
    assert noise["passed"] is False and noise["below_expected"] is True


def test_out_of_fold_scores_do_not_see_their_own_rows():
    """A perfectly memorisable label must NOT be reproduced out of fold."""
    rng = np.random.default_rng(11)
    n = 200
    labels = rng.random(n) < 0.5
    features = np.eye(n)                     # one free parameter per patient
    assert abs(_auroc(labels, out_of_fold_scores(features, labels, seed=11)) - 0.5) < 0.12


# --- T1.5 gene-label shuffle ----------------------------------------------


def test_gene_label_shuffle_preserves_geometry_and_destroys_attribution():
    """Both halves of the T1.5 criterion, on the operation itself."""
    from dataclasses import replace
    from morpheus.v2.pbs import ReferenceDictionary
    rng = np.random.default_rng(12)
    genes = [f"G{i}" for i in range(60)]
    responses = rng.normal(size=(40, 60))
    dictionary = ReferenceDictionary.fit(responses, genes, [f"A{i}" for i in range(40)], n_components=6)
    order = rng.permutation(len(genes))
    shuffled = replace(dictionary, gene_mean=dictionary.gene_mean[order],
                       gene_basis=dictionary.gene_basis[order])
    # (i) the directions, their norms and their mutual angles are untouched
    np.testing.assert_allclose(np.sort(np.abs(shuffled.gene_basis), axis=0),
                               np.sort(np.abs(dictionary.gene_basis), axis=0), atol=1e-12)
    np.testing.assert_allclose(shuffled.gene_basis.T @ shuffled.gene_basis,
                               dictionary.gene_basis.T @ dictionary.gene_basis, atol=1e-10)
    # (ii) per-axis gene attribution collapses
    from scipy.stats import spearmanr
    rho = [abs(spearmanr(dictionary.gene_basis[:, k], shuffled.gene_basis[:, k]).statistic)
           for k in range(dictionary.components)]
    assert max(rho) < 0.5 and float(np.median(rho)) < 0.2, rho


# --- T1.6 pairing null resolution -----------------------------------------


def test_permutation_p_cannot_beat_its_own_resolution():
    """G4.5, as an assertion. A headline 'must fail' cannot rest on 50 draws."""
    from morpheus.v2.calibra.calibration import permutation_null
    rng = np.random.default_rng(13)
    z = rng.normal(size=300)
    x = rng.normal(size=(300, 8)); x[:, 0] += 3 * z
    y = rng.normal(size=(300, 6)); y[:, 0] += 3 * z
    design = np.zeros((300, 0))
    coarse = permutation_null(x, y, design, n_permutations=20, n_components=4, seed=13)
    fine = permutation_null(x, y, design, n_permutations=200, n_components=4, seed=13)
    assert coarse["permutation_p"] >= 1.0 / 21.0
    assert fine["permutation_p"] >= 1.0 / 201.0
    assert fine["permutation_p"] < coarse["permutation_p"]
    # and the null must not be degenerate, or the p means nothing
    assert fine["null_max"] > fine["null_median"] > 0.0


def test_pairing_null_collapses_a_real_signal_to_the_capacity_floor():
    """Shuffled pairing must destroy the agreement, and the residue it leaves is
    the capacity floor -- not zero. Quoting zero is the error this guards."""
    from morpheus.v2.calibra.calibration import permutation_null
    rng = np.random.default_rng(14)
    n = 400
    z = rng.normal(size=n)
    x = rng.normal(size=(n, 20)); x[:, 0] += 4 * z
    y = rng.normal(size=(n, 20)); y[:, 0] += 4 * z
    result = permutation_null(x, y, np.zeros((n, 0)), n_permutations=200, n_components=8, seed=14)
    assert result["observed_top_cca"] > result["null_max"], result
    assert result["null_median"] > 0.05, "capacity floor is not zero and must be reported as such"


# --- zero-parameter baseline ----------------------------------------------


def test_zero_parameter_baseline_is_constant_within_cancer_and_uses_train_only():
    from morpheus.v2.baseline_exports import export_zero_parameter_baseline
    import tempfile
    rng = np.random.default_rng(15)
    n = 200
    cancers = np.asarray(["A", "B"] * (n // 2))
    split = np.asarray(["train"] * 120 + ["test"] * 80)
    ids = np.asarray([f"TCGA-{10 + (i % 3):02d}-{i:04d}" for i in range(n)])
    targets = rng.normal(size=(n, 5)) + (cancers == "A")[:, None] * 3.0
    with tempfile.TemporaryDirectory() as directory:
        path = export_zero_parameter_baseline(f"{directory}/zp.npz", patient_ids=ids, cancers=cancers,
                                              split=split, targets=targets)
        state = np.load(path, allow_pickle=True)["wsi_identity"]
    for cancer in ("A", "B"):
        rows = state[cancers == cancer]
        assert np.allclose(rows, rows[0]), "the naive baseline must be constant within cancer"
    # and it must equal the TRAIN-fold mean, never the test-fold mean
    train_mean = targets[(cancers == "A") & (split == "train")].mean(axis=0)
    np.testing.assert_allclose(state[cancers == "A"][0], train_mean, rtol=1e-5)
