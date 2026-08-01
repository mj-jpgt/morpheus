import numpy as np

from morpheus.v2.paired_bootstrap import paired_multivariate_patient_and_cancer_bootstrap


def test_paired_multivariate_bootstrap_resamples_patient_axis_for_both_modes():
    rng = np.random.default_rng(7)
    y = rng.normal(size=(24, 5))
    teacher = y + .6 * rng.normal(size=(24, 5))
    challenger = y + .2 * rng.normal(size=(24, 5))
    metric = lambda actual, representation: float(np.corrcoef(actual.ravel(), representation.ravel())[0, 1])
    result = paired_multivariate_patient_and_cancer_bootstrap(
        metric, y, teacher, challenger, np.repeat(["A", "B", "C"], 8), repeats=100, seed=3)
    assert result["patient"]["n_valid"] == 100
    assert result["cancer"]["n_valid"] == 100
    assert result["patient"]["point_delta"] > 0
