import numpy as np
import pytest

from morpheus.v2.pbs import LegibilityOperator, ReferenceDictionary


def test_reference_dictionary_requires_exact_gene_identity_for_image_free_codes():
    rng = np.random.default_rng(0)
    dictionary = ReferenceDictionary.fit(rng.normal(size=(20, 8)), [f"G{i}" for i in range(8)], [f"A{i}" for i in range(20)], n_components=4)
    codes = dictionary.encode_expression(rng.normal(size=(6, 8)), [f"G{i}" for i in range(8)])
    assert codes.shape == (6, 4)
    with pytest.raises(ValueError, match="exactly match"):
        dictionary.encode_expression(rng.normal(size=(6, 8)), [f"G{i}" for i in reversed(range(8))])


def test_legibility_operator_is_cross_cancer_fit_and_nonnegative():
    rng = np.random.default_rng(1)
    wsi = rng.normal(size=(60, 6)); codes = wsi[:, :3] @ rng.normal(size=(3, 4)) + .1 * rng.normal(size=(60, 4))
    operator = LegibilityOperator.fit(wsi, codes, np.repeat(["A", "B", "C"], 20))
    assert operator.weights.shape == (4,)
    assert np.all((0 <= operator.weights) & (operator.weights <= 1))
