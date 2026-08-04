"""The centring measurement: what the two cosines are, and that the verdict is computed.

`v2/research/rebase/p2/p2_centred_cosine.py` exists to settle one question the
draft refuses to adjudicate — whether §5.2a's dissociation between a flat centred
rank and a moving uncentred cosine is a real insensitivity of rank or an artefact
of the fact that one statistic centres and the other does not.

Four things have to hold or the measurement cannot carry that:

* **the uncentred form is the harness's form.** `d1_momentum_probe.geometry()`
  prints the mean off-diagonal cosine of the row-normalised RNA states; if this
  module's `centre=False` branch is anything else, it is not recomputing the
  number the draft quotes;
* **centring is the only difference between the two forms**, and it is the
  difference the question is about: the centred form must be blind to a shared
  offset and the uncentred form must not be;
* **the collapse family behaves as §4.10 says it does.** `z_i = m + a_i·u` is the
  family this project documents, and the whole mean-offset account is the claim
  that a change confined to `m` moves the uncentred cosine and nothing else. That
  is pinned here on synthetic data, so the account is falsifiable independently of
  what the GPU returned;
* **the verdict is computed, not asserted.** The three branches of the
  predeclared rule are exercised on synthetic spreads.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2 import p2_centred_cosine as CC


def _harness_cosine(x: np.ndarray) -> float:
    """`d1_momentum_probe.geometry()`'s `rna-rna` column, written out longhand.

    Deliberately an independent restatement rather than an import: this is the
    one place a second implementation is wanted, because its whole job is to
    disagree if the module's has drifted.
    """
    unit = x / np.linalg.norm(x, axis=-1, keepdims=True)
    gram = unit @ unit.T
    eye = np.eye(len(x), dtype=bool)
    return float(gram[~eye].mean())


@pytest.fixture(scope="module")
def rng():
    return np.random.default_rng(20260804)


# --------------------------------------------------------------------------
def test_the_uncentred_form_is_the_column_the_harness_prints(rng):
    x = rng.normal(size=(64, 32))
    assert CC.mutual_cosine(x, centre=False) == pytest.approx(_harness_cosine(x), abs=1e-12)


def test_centring_is_the_only_difference_and_it_is_the_one_that_matters(rng):
    """A shared offset moves the uncentred cosine and cannot move the centred one.

    This is the mean-offset account stated as a property rather than as a story:
    if it did not hold, the measurement would not separate the two accounts at
    all.
    """
    x = rng.normal(size=(128, 32))
    offset = rng.normal(size=(1, 32)) * 5.0
    shifted = x + offset

    assert CC.mutual_cosine(shifted, centre=True) == pytest.approx(
        CC.mutual_cosine(x, centre=True), abs=1e-10)
    assert abs(CC.mutual_cosine(shifted, centre=False)
               - CC.mutual_cosine(x, centre=False)) > 0.5


def test_the_collapse_family_this_project_documents_behaves_as_4_10_says(rng):
    """`z_i = m + a_i·u`: uncentred cosine near 1, centred rank near 1.

    §4.10 states the family and the reason we centre: the column mean is a rank-1
    component carrying no between-sample information, and this family's uncentred
    effective rank is ~2 where its centred rank is ~1. Three offsets of
    decreasing size are scored, and the uncentred cosine falls sharply across
    them while the centred rank does not move — which is exactly the pattern
    §5.2a reports and account (B) attributes to the offset.
    """
    n, d = 128, 32
    u = rng.normal(size=d)
    u /= np.linalg.norm(u)
    a = rng.normal(size=(n, 1))
    mean_direction = rng.normal(size=d)
    mean_direction /= np.linalg.norm(mean_direction)

    uncentred, centred, ranks = [], [], []
    for scale in (20.0, 4.0, 0.5):
        z = scale * mean_direction[None, :] + a * u[None, :]
        uncentred.append(CC.mutual_cosine(z, centre=False))
        centred.append(CC.mutual_cosine(z, centre=True))
        ranks.append(effective_rank(z, variant=RANK_VARIANTS["R1"]))

    assert uncentred[0] > 0.99, uncentred
    assert uncentred[0] - uncentred[-1] > 0.4, uncentred
    # The centred rank is ~1 at every offset: centring removes the offset and
    # what is left is one direction.
    assert max(ranks) / min(ranks) < 1.05, ranks
    assert max(ranks) < 1.2, ranks
    # And the centred cosine does not track the offset either.
    assert max(centred) - min(centred) < 0.05, centred


def test_the_mean_offset_ratio_is_the_size_of_the_shared_offset(rng):
    n, d = 256, 16
    u = rng.normal(size=(n, d))
    u -= u.mean(axis=0, keepdims=True)
    u /= np.sqrt((u ** 2).sum(axis=1).mean())
    offset = rng.normal(size=d)
    offset /= np.linalg.norm(offset)
    for scale in (0.5, 2.0, 8.0):
        assert CC.mean_offset_ratio(u + scale * offset[None, :]) == pytest.approx(
            scale, rel=1e-6)


# --------------------------------------------------------------------------
# The verdict is computed from the numbers, not asserted beside them
# --------------------------------------------------------------------------
def _synthetic(across: dict[str, float], within: dict[str, float],
               uncentred_across: float = 0.47, uncentred_within: float = 0.01) -> dict:
    """A `within`/`across` pair in the shape `verdict()` reads, at one step."""
    step, view = "200", CC.COSINE_VIEW

    def entry(lo, hi):
        return {"min": lo, "max": hi, "fold": None, "values": {}, "shape": {}}

    w = {arm: {step: {view: {
        "mutual_cosine_centred": entry(0.0, spread),
        "mutual_cosine_uncentred": entry(0.0, uncentred_within)}}}
        for arm, spread in within.items()}
    a = {rep: {step: {view: {
        "mutual_cosine_centred": {"spread": spread, "min": 0.0, "max": spread},
        "mutual_cosine_uncentred": {"spread": uncentred_across, "min": 0.0,
                                    "max": uncentred_across}}}}
        for rep, spread in across.items()}
    return CC.verdict(w, a, 200)


def test_the_moves_branch_needs_the_smallest_repeat_draw_to_qualify():
    """(A) is taken only if every repeat draw shows the movement and clears the floor."""
    got = _synthetic({"rep1": 0.31, "rep2": 0.28, "rep3": 0.35}, {"m0": 0.02, "m0999": 0.03})
    assert got["reading"].startswith("A —"), got["reading"]

    # one flat draw is enough to refuse it
    got = _synthetic({"rep1": 0.31, "rep2": 0.05, "rep3": 0.35}, {"m0": 0.02, "m0999": 0.03})
    assert got["reading"].startswith("neither"), got["reading"]

    # movement that does not clear its own three-repeat floor is not a movement
    got = _synthetic({"rep1": 0.31, "rep2": 0.28, "rep3": 0.35}, {"m0": 0.02, "m0999": 0.40})
    assert got["reading"].startswith("neither"), got["reading"]


def test_the_dissolves_branch_needs_the_largest_repeat_draw_to_stay_flat():
    """(B) is taken only if even the widest draw is inside 0.10 and inside the floor."""
    got = _synthetic({"rep1": 0.03, "rep2": 0.02, "rep3": 0.04}, {"m0": 0.06, "m0999": 0.05})
    assert got["reading"].startswith("B —"), got["reading"]

    # inside 0.10 but OUTSIDE its own floor: not flat, not adjudicated
    got = _synthetic({"rep1": 0.03, "rep2": 0.02, "rep3": 0.09}, {"m0": 0.01, "m0999": 0.01})
    assert got["reading"].startswith("neither"), got["reading"]

    # and (B) requires the uncentred movement to be real in the first place
    got = _synthetic({"rep1": 0.03, "rep2": 0.02, "rep3": 0.04}, {"m0": 0.06, "m0999": 0.05},
                     uncentred_across=0.05)
    assert got["reading"].startswith("neither"), got["reading"]


def test_the_thresholds_are_the_predeclared_ones():
    got = _synthetic({"rep1": 0.01}, {"m0": 0.5})
    assert got["moves_threshold"] == 0.20 and got["flat_threshold"] == 0.10
    assert "PREDECLARED_centred_cosine_20260804T1700Z" in got["rule"]
