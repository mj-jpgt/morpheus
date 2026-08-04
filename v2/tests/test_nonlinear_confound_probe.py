"""The probes must (a) agree with the estimator the spatial run used, (b) see a confound
that lives outside the first moment, and (c) report a null that is a chance rate only
when it is one."""
import numpy as np
import pytest

from morpheus.v2.calibra.confound_certificate import (lda_oof_balanced_accuracy,
                                                      nearest_class_mean_oof,
                                                      within_stratum_permutations)
from morpheus.v2.calibra.hest_claims import knn_balanced_accuracy_oof as hest_knn
from morpheus.v2.calibra.nonlinear_confound_probe import (forest_balanced_accuracy_oof,
                                                          global_permutations,
                                                          knn_balanced_accuracy_oof,
                                                          nesting_diagnostic, null_summary,
                                                          probe_state,
                                                          rbf_svm_balanced_accuracy_oof,
                                                          within_stratum_chance)


def _mean_shift_fixture(n_classes=4, per_class=40, dim=6, seed=0):
    rng = np.random.default_rng(seed)
    labels = np.repeat(np.arange(n_classes), per_class)
    centres = rng.normal(scale=3.0, size=(n_classes, dim))
    features = centres[labels] + rng.normal(scale=1.0, size=(len(labels), dim))
    return features, labels, n_classes


def _second_moment_fixture(n_classes=4, per_class=120, dim=4, seed=1):
    """Classes with IDENTICAL means, differing only in variance.

    This is the object the whole question is about: a confound that a mean-based
    certificate cannot see and a mean-removing adjustment cannot remove.
    """
    rng = np.random.default_rng(seed)
    labels = np.repeat(np.arange(n_classes), per_class)
    scales = np.linspace(0.25, 4.0, n_classes)
    features = rng.normal(size=(len(labels), dim)) * scales[labels][:, None]
    # Force the class means to be exactly equal so no residue of the first moment remains.
    for value in range(n_classes):
        rows = labels == value
        features[rows] -= features[rows].mean(axis=0)
    return features, labels, n_classes


# --- (a) the plain vote is the SAME estimator the spatial claim used --------------------

def test_plain_vote_reproduces_the_canonical_hest_implementation_exactly():
    features, labels, n_classes = _mean_shift_fixture()
    for k in (1, 5, 15):
        mine = knn_balanced_accuracy_oof(features, labels, n_classes, k=k, seed=42,
                                         prior_correction=False)
        theirs = hest_knn(features, labels, n_classes, k=k, seed=42)
        assert mine == theirs, f"k={k}: probe drifted from the spatial estimator"


def test_prior_correction_changes_the_answer_under_imbalance_and_helps_the_small_class():
    rng = np.random.default_rng(3)
    # 300 rows of class 0, 20 each of classes 1 and 2, all separable.
    labels = np.concatenate([np.zeros(300, int), np.ones(20, int), np.full(20, 2)])
    centres = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]])
    features = centres[labels] + rng.normal(scale=0.9, size=(len(labels), 2))
    plain = knn_balanced_accuracy_oof(features, labels, 3, k=15, prior_correction=False)
    weighted = knn_balanced_accuracy_oof(features, labels, 3, k=15, prior_correction=True)
    assert weighted > plain, "inverse-frequency voting must not lose to the majority vote here"


def test_knn_is_deterministic_under_a_fixed_seed():
    features, labels, n_classes = _mean_shift_fixture()
    a = knn_balanced_accuracy_oof(features, labels, n_classes, k=5, seed=7)
    b = knn_balanced_accuracy_oof(features, labels, n_classes, k=5, seed=7)
    assert a == b


# --- (b) the probes see what the certificate cannot -------------------------------------

def test_the_certificate_is_blind_to_a_variance_only_confound_and_the_probes_are_not():
    features, labels, n_classes = _second_moment_fixture()
    chance = 1.0 / n_classes
    lda = lda_oof_balanced_accuracy(features, labels, n_classes)
    per_axis = float(np.max(nearest_class_mean_oof(features, labels, n_classes)))
    # The certificate's two classifiers must be AT chance on this fixture, not merely below
    # the probes: if they were not, the fixture would still carry a first moment and the
    # test would prove nothing about estimator families.
    assert lda < 1.2 * chance, f"fixture is not mean-free: LDA reads {lda}"
    assert per_axis < 1.2 * chance, f"fixture is not mean-free: per-axis reads {per_axis}"
    for name, value in (("knn", knn_balanced_accuracy_oof(features, labels, n_classes, k=15,
                                                          prior_correction=True)),
                        ("forest", forest_balanced_accuracy_oof(features, labels, n_classes,
                                                                n_estimators=120)),
                        ("rbf_svm", rbf_svm_balanced_accuracy_oof(features, labels, n_classes))):
        assert value > 2.0 * chance, f"{name} missed a variance-only confound: {value}"


def test_every_probe_recovers_a_planted_mean_shift_and_returns_a_balanced_accuracy():
    features, labels, n_classes = _mean_shift_fixture()
    for value in (knn_balanced_accuracy_oof(features, labels, n_classes, k=5),
                  forest_balanced_accuracy_oof(features, labels, n_classes, n_estimators=100),
                  rbf_svm_balanced_accuracy_oof(features, labels, n_classes)):
        assert 0.0 <= value <= 1.0
        assert value > 0.7


def test_probes_sit_at_chance_when_the_label_carries_no_information():
    rng = np.random.default_rng(11)
    features = rng.normal(size=(240, 5))
    labels = np.repeat(np.arange(4), 60)
    chance = 0.25
    for value in (knn_balanced_accuracy_oof(features, labels, 4, k=15, prior_correction=True),
                  forest_balanced_accuracy_oof(features, labels, 4, n_estimators=120)):
        assert abs(value - chance) < 0.12, value


# --- (c) the nulls, and when one of them stops meaning "chance" -------------------------

def test_the_global_null_measures_chance_rather_than_assuming_it():
    features, labels, n_classes = _mean_shift_fixture()
    permutations = global_permutations(labels, n_permutations=40, seed=42)
    null = np.asarray([knn_balanced_accuracy_oof(features, p, n_classes, k=5)
                       for p in permutations])
    assert abs(float(np.median(null)) - 1.0 / n_classes) < 0.08
    assert len(permutations) == 40
    assert all(np.array_equal(np.sort(p), np.sort(labels)) for p in permutations)


def test_a_within_stratum_null_is_degenerate_when_a_stratum_holds_one_class():
    labels = np.array([0, 0, 1, 1, 2, 2, 3, 3])
    strata = np.array(["a", "a", "a", "a", "b", "b", "b", "b"])
    lone = np.array(["a", "a", "a", "a", "b", "b", "b", "b"])
    labels_lone = np.array([0, 0, 1, 1, 2, 2, 2, 2])   # stratum "b" holds only class 2
    healthy = nesting_diagnostic(labels, strata, 4)
    degenerate = nesting_diagnostic(labels_lone, lone, 3)
    assert healthy["fraction_of_rows_unpermutable"] == 0.0
    assert degenerate["fraction_of_rows_unpermutable"] == 0.5
    permuted = within_stratum_permutations(labels_lone, lone, n_permutations=5, seed=0)
    assert all(np.array_equal(p[4:], labels_lone[4:]) for p in permuted), \
        "a single-class stratum must be permuted by the identity - that is the trap"


def test_nesting_is_detected_and_raises_the_within_stratum_chance_above_the_design_chance():
    # Four classes, two strata, each class living in exactly one stratum: total nesting.
    labels = np.repeat(np.arange(4), 10)
    strata = np.where(labels < 2, "s0", "s1")
    diagnostic = nesting_diagnostic(labels, strata, 4)
    assert diagnostic["confound_nests_inside_stratum"] is True
    assert diagnostic["n_classes_spanning_multiple_strata"] == 0
    assert diagnostic["design_chance_rate"] == pytest.approx(0.25)
    # knowing the stratum halves the guessing problem: 1/2, not 1/4
    assert diagnostic["within_stratum_chance"] == pytest.approx(0.5)
    assert within_stratum_chance(labels, strata, 4) == pytest.approx(0.5)


def test_a_class_spread_over_strata_is_not_reported_as_nested():
    labels = np.repeat(np.arange(3), 8)
    strata = np.tile(np.array(["s0", "s1"]), 12)
    diagnostic = nesting_diagnostic(labels, strata, 3)
    assert diagnostic["confound_nests_inside_stratum"] is False
    assert diagnostic["n_classes_spanning_multiple_strata"] == 3
    assert diagnostic["within_stratum_chance"] == pytest.approx(1.0 / 3.0)


def test_null_summary_floors_the_p_value_and_reports_its_own_resolution():
    summary = null_summary(np.zeros(99), 1.0)
    assert summary["permutation_p"] == pytest.approx(1.0 / 100.0)
    assert summary["permutation_resolution"] == pytest.approx(1.0 / 100.0)
    assert summary["n_permutations"] == 99


# --- the assembled call ------------------------------------------------------------------

def test_probe_state_reports_both_nulls_the_linear_reference_and_every_k():
    features, labels, n_classes = _mean_shift_fixture(per_class=30)
    strata = np.where(labels < 2, "s0", "s1")
    result = probe_state(features, labels, n_classes, strata=strata, k_grid=(3, 9),
                         forest_trees=40, knn_permutations=6, model_permutations=4)
    assert set(result["probes"]) == {"knn_k3", "knn_prior_k3", "knn_k9", "knn_prior_k9", "forest"}
    assert result["chance_rate"] == pytest.approx(1.0 / n_classes)
    assert "joint_lda_balanced_accuracy" in result["linear_reference"]
    for name, record in result["probes"].items():
        assert {"global", "within_stratum"} <= set(record), name
        assert record["global"]["n_permutations"] == (6 if name.startswith("knn") else 4)
    assert result["probes"]["knn_k3"]["multiple_of_chance"] == pytest.approx(
        result["probes"]["knn_k3"]["balanced_accuracy"] * n_classes)


def test_probe_state_omits_the_within_stratum_null_when_no_stratum_is_given():
    features, labels, n_classes = _mean_shift_fixture(per_class=25)
    result = probe_state(features, labels, n_classes, strata=None, k_grid=(5,), forest_trees=0,
                         knn_permutations=4, model_permutations=0)
    assert set(result["probes"]) == {"knn_k5", "knn_prior_k5"}
    assert "within_stratum" not in result["probes"]["knn_k5"]
    assert "global" in result["probes"]["knn_k5"]
