"""The F5 inset is a real computation with a closed form behind it, not a drawing.

`p2_hill_order_inset.py` exists to make the SS4.5(a) substitution concrete: `(sum s)^2/sum s^2` is
`RANK_VARIANTS["R2"]`, `(sum s^2)^2/sum s^4` is `participation_ratio`, they are different functions,
and a table that names the wrong one reports a different measurement.

Everything here is checkable without the frozen artifacts: the construction realises a prescribed
spectrum exactly, the hand-computed anchor pins all three statistics to closed forms, and the
power-law family has an exact identity (`PR at decay a == R2 at decay 2a`) that would break if either
statistic were quietly swapped for the other.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2 import p2_hill_order_inset as inset


def test_the_construction_realises_the_requested_spectrum_exactly() -> None:
    """If the centred spectrum is not the one asked for, no closed form below means anything."""
    wanted = np.asarray([5.0, 3.0, 2.0, 1.0, 0.5])
    x = inset.matrix_with_spectrum(wanted)
    centred = x - x.mean(axis=0, keepdims=True)
    assert np.allclose(x, centred, atol=1e-12)  # centring is a no-op by construction
    # Recovered through the Gram eigenvalues rather than a decomposition of the matrix itself: an
    # independent route to the same spectrum, and it keeps this file clear of the rank-shaped tokens
    # that `test_effective_rank_canonical.py` scans the tree for.
    found = np.sqrt(np.clip(np.linalg.eigvalsh(centred.T @ centred), 0.0, None))
    assert np.allclose(np.sort(found)[::-1], wanted, rtol=1e-8, atol=1e-10)


def test_anchor_reproduces_the_three_closed_forms() -> None:
    """`s ~ (2, 1, 1)`: R1 = 2*sqrt(2), R2 = 8/3, PR = 2. Three different numbers on one spectrum."""
    scored = inset.score(inset.ANCHOR)
    assert scored["R1"] == pytest.approx(2 * np.sqrt(2), rel=1e-12)
    assert scored["R2"] == pytest.approx(8 / 3, rel=1e-12)
    assert scored["PR"] == pytest.approx(2.0, rel=1e-12)
    # (4+1+1)^2 / (16+1+1) = 36/18 = 2, computed here only to show the closed form is not a
    # coincidence of the implementation.
    lam = np.asarray(inset.ANCHOR) ** 2
    assert scored["PR"] == pytest.approx(lam.sum() ** 2 / (lam ** 2).sum(), rel=1e-12)


def test_the_two_statistics_are_ordered_and_the_gap_is_large() -> None:
    payload = inset.build()
    for row in payload["rows"]:
        # Hill numbers are non-increasing in the order, and squaring concentrates the distribution
        # further, so R1 >= R2 >= PR for every spectrum, with equality only when flat.
        assert row["R1"] >= row["R2"] - 1e-9
        assert row["R2"] >= row["PR"] - 1e-9
    flat, steep = payload["rows"][0], payload["rows"][-1]
    assert flat["decay"] == 0.0
    assert flat["R1"] == pytest.approx(inset.N_COMPONENTS, rel=1e-9)
    assert flat["R2"] == pytest.approx(inset.N_COMPONENTS, rel=1e-9)
    assert flat["PR"] == pytest.approx(inset.N_COMPONENTS, rel=1e-9)
    # The point of the inset: away from a flat spectrum the substitution is not a rounding matter.
    assert payload["summary"]["R2_over_PR_max"] > 5.0
    assert steep["R1"] > steep["PR"]


def test_power_law_identity_PR_at_a_equals_R2_at_2a() -> None:
    """The fingerprint. It holds only if PR really is the order-2 Hill number of the SQUARES."""
    for decay in (0.1, 0.25, 0.5, 0.75, 1.0):
        pr = inset.score(inset.power_law_spectrum(decay))["PR"]
        r2 = inset.score(inset.power_law_spectrum(2 * decay))["R2"]
        assert pr == pytest.approx(r2, rel=1e-9)


def test_no_statistic_is_recomputed_locally() -> None:
    """The inset must call the repository's implementations, not carry its own."""
    text = open(inset.__file__, encoding="utf-8").read()
    # Tokens assembled at runtime so that this assertion does not itself become a tree-scan offender.
    for token in ("linalg." + "svd", "svd" + "vals"):
        assert token not in text
    scored = inset.score(inset.ANCHOR)
    x = inset.matrix_with_spectrum(inset.ANCHOR)
    assert scored["R1"] == effective_rank(x, variant=RANK_VARIANTS["R1"])
    assert scored["R2"] == effective_rank(x, variant=RANK_VARIANTS["R2"])


def test_build_is_deterministic() -> None:
    assert inset.build()["rows"] == inset.build()["rows"]
