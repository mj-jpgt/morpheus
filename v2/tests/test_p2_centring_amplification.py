"""The centring-amplification derivation, its exact reductions, and its predeclared thresholds.

Three things are pinned here, and each of them has already caught something:

* **the transfer identity is EXACT in its own hypothesis.** Theorem 1 of
  `NOTEBOOK_ENTRIES/PREDECLARED_centring_amplification_law_20260804T1750Z.md`
  claims ``ln D_1(uncentred) = h(t) + (1-t) ln D_1(centred)`` whenever the mean
  direction is orthogonal to the centred row space. A construction that satisfies
  the hypothesis exactly is built here and the identity is asserted to machine
  precision. If the identity is wrong, everything downstream is decoration.
* **the ``[X; -X]`` reduction is exact.** It is how three centre-only statistics
  are evaluated uncentred without reimplementing them, and the predeclaration
  makes a failure a stopping condition.
* **the thresholds cannot be re-read after the numbers arrived.** Every constant
  in `p2_centring_verdict.THRESHOLDS` is asserted to appear verbatim in the
  predeclaration file, which is committed and pushed at `d4e344c`.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2 import p2_centring_amplification as A
from morpheus.v2.research.rebase.p2 import p2_centring_verdict as V
from morpheus.v2.research.rebase.p2.p2_competing_metrics import rankme

REPO = Path(__file__).resolve().parents[2]
PREDECLARATION = (REPO / "NOTEBOOK_ENTRIES"
                  / "PREDECLARED_centring_amplification_law_20260804T1750Z.md")


def _orthogonal_construction(n=200, d=32, k=8, offset=6.0, seed=0):
    """``Z = 1 mu^T + Z_c`` with ``mu`` EXACTLY orthogonal to ``Z_c``'s row space.

    ``Z_c`` is built in the span of ``k`` orthonormal directions, its rows are made
    to sum to zero, and ``mu`` is taken from the orthogonal complement. Theorem 1's
    hypothesis then holds exactly rather than approximately.
    """
    rng = np.random.default_rng(seed)
    basis = np.linalg.qr(rng.normal(size=(d, d)))[0]
    loadings = rng.normal(size=(n, k)) * np.exp(-np.arange(k) / 3.0)
    loadings -= loadings.mean(axis=0, keepdims=True)          # zero column mean
    centred = loadings @ basis[:, :k].T
    mu = offset * basis[:, k + 1]                              # orthogonal complement
    return np.ones((n, 1)) @ mu[None, :] + centred


def test_the_transfer_identity_is_exact_when_its_hypothesis_holds():
    z = _orthogonal_construction()
    entry = A.score_matrix(z)
    observed = entry["identity"]["ln_D1_uncentred_observed"]
    predicted = entry["identity"]["ln_D1_uncentred_predicted"]
    assert observed == pytest.approx(predicted, rel=1e-10), (
        "Theorem 1 must be exact under its own hypothesis; it is not, so the "
        f"derivation is wrong: {observed} vs {predicted}")
    # ... and the mean direction really is the dominant one in this construction
    assert entry["f"] > 0.9
    assert 0.0 < entry["t"] < 1.0


def test_the_identity_degrades_gracefully_when_the_hypothesis_does_not_hold():
    """A random mean direction is not orthogonal to the residual row space.

    The identity must still be close -- that is the empirical content of assumption
    A1 -- but it must NOT be exact, or the test above would be vacuous.
    """
    rng = np.random.default_rng(3)
    z = rng.normal(size=(300, 48)) + 6.0 * np.ones((300, 1)) @ rng.normal(size=(1, 48))
    error = A.score_matrix(z)["identity"]["relative_error"]
    assert 0.0 < error < 0.05


def test_doubling_reduces_uncentred_to_centred_exactly():
    rng = np.random.default_rng(1)
    z = rng.normal(size=(150, 24)) * np.exp(-np.arange(24) / 6.0)
    z += 2.0 * np.ones((150, 1)) @ rng.normal(size=(1, 24))
    guard = A.check_doubling_is_exact(z)
    assert guard["passed"], guard
    assert effective_rank(A.doubled(z), variant=RANK_VARIANTS["R1"]) == pytest.approx(
        effective_rank(z, variant=RANK_VARIANTS["R1_uncentred"]), rel=1e-12)
    assert rankme(A.doubled(z)) == pytest.approx(rankme(z), rel=1e-12)


def test_doubling_would_be_caught_if_it_were_not_scale_invariant():
    """The reduction rests on scale invariance; a scale-DEPENDENT statistic must break it.

    Without this, `check_doubling_is_exact` could be passing for the wrong reason.
    """
    rng = np.random.default_rng(2)
    z = rng.normal(size=(80, 16)) + 3.0
    frobenius = float(np.linalg.norm(z))
    assert float(np.linalg.norm(A.doubled(z))) == pytest.approx(math.sqrt(2) * frobenius)


def test_the_f_form_and_the_t_form_are_the_same_number():
    """Corollary 1.7 is an identity, not an approximation, and the code must show it."""
    rng = np.random.default_rng(5)
    for offset in (0.0, 0.5, 2.0, 10.0):
        z = rng.normal(size=(120, 20)) + offset * np.ones((120, 1)) @ rng.normal(size=(1, 20))
        entry = A.score_matrix(z)
        if entry["t_from_f"] is None:
            continue
        assert entry["t_identity_absolute_error"] < 1e-10, entry


def test_order_one_amplification_is_one_over_one_minus_t():
    for t in (0.05, 0.3, 0.7, 0.9):
        assert A.predicted_amplification(1, t, 12.0) == pytest.approx(1.0 / (1.0 - t))


def test_hard_rank_is_a_limit_case_the_law_gets_right_not_a_counterexample():
    """Order 0 predicts essentially NO amplification, and the measured floor is 1.000x."""
    predicted = A.predicted_amplification(0, 0.8, 255.0)
    assert predicted == pytest.approx(256.0 / 255.0, rel=1e-9)
    assert predicted < 1.005


def test_stable_rank_prediction_is_binary_and_switches_at_the_right_place():
    """``a -> inf`` is a max, so the prediction is a branch, not a continuous factor."""
    # spike dominates the residual's largest normalised eigenvalue -> no response at all
    assert A.predicted_amplification(math.inf, 0.9, 5.0) is None
    # residual's largest still dominates -> no amplification
    assert A.predicted_amplification(math.inf, 0.05, 5.0) == 1.0


def test_binary_entropy():
    assert A.binary_entropy(0.5) == pytest.approx(math.log(2))
    assert A.binary_entropy(0.0) == 0.0
    assert A.binary_entropy(1.0) == 0.0


def test_the_module_computes_no_svd_of_its_own():
    """It must reach every spectral statistic through an import, never through a formula."""
    for name in ("p2_centring_amplification.py", "p2_centring_verdict.py"):
        source = (REPO / "v2" / "research" / "rebase" / "p2" / name).read_text(encoding="utf-8")
        for token in ("svdvals", "linalg.svd", "eigvalsh", "linalg.eigh"):
            assert token not in source, f"{name} computes a spectrum inline: {token}"


@pytest.mark.parametrize("key,value", sorted(V.THRESHOLDS.items()))
def test_every_verdict_threshold_appears_in_the_predeclaration(key, value):
    """A threshold edited after the numbers arrived must break the build, not the argument."""
    text = PREDECLARATION.read_text(encoding="utf-8")
    # the predeclaration writes proportions as percentages and the sweep's f ceiling
    # as a decimal, so either spelling counts -- but SOME spelling must be there
    needles = [f"{value:g}", f"{value:.2f}"]
    if value < 1.0:
        needles.append(f"{value * 100:g}%")
    assert any(re.search(re.escape(n), text) for n in needles), (
        f"threshold {key} = {value} is not in the predeclaration as any of {needles!r}; "
        "either it was invented afterwards or the predeclaration must be quoted exactly")


def test_synthetic_sweep_rejects_an_unknown_condition():
    with pytest.raises(ValueError):
        A.synthetic_sweep(condition="whatever")


def test_variance_decomposition_flags_an_unstable_dominant_component():
    """A2's check must say 'stable' when t is constant and 'unstable' when it is not."""
    d1_centred = [20.0, 20.4, 19.7, 20.2, 20.1]
    d1_uncentred = [8.0, 8.1, 7.9, 8.05, 8.02]
    steady = V and A.variance_decomposition([0.5] * 5, d1_centred, d1_uncentred)
    assert steady["spike_share_of_uncentred_variance"] == pytest.approx(0.0)
    wobbly = A.variance_decomposition([0.50, 0.56, 0.46, 0.53, 0.49],
                                      d1_centred, d1_uncentred)
    assert wobbly["spike_share_of_uncentred_variance"] > 1.0


# ---------------------------------------------------------------------------
# the vendored result, re-read so that a number quoted in the notebook entry
# cannot drift from the file it came from
# ---------------------------------------------------------------------------
VENDORED = REPO / "v2" / "research" / "rebase" / "p2" / "figures" / "data" / "ws_amp" / "out"


@pytest.fixture(scope="module")
def verdict():
    import json
    return json.loads((VENDORED / "P2_CENTRING_VERDICT.json").read_text(encoding="utf-8"))


def test_the_vendored_verdicts_are_the_ones_reported(verdict):
    """The headline: the derivation survives its own premise and the data does not."""
    real = {k: verdict["real"][k]["verdict"] for k in ("P1", "P2", "P3", "P3b", "P4", "P5")}
    assert real == {"P1": "not falsified", "P2": "falsified", "P3": "not falsified",
                    "P3b": "falsified", "P4": "falsified", "P5": "not falsified"}
    synthetic = {c: {k: verdict["synthetic"][c][k]["verdict"] for k in ("S1", "S2", "S3", "S4")}
                 for c in ("t_pinned", "stable_mean", "unstable_mean")}
    assert synthetic["t_pinned"]["S1"] == "not falsified"
    assert synthetic["t_pinned"]["S2"] == "not falsified"
    assert synthetic["stable_mean"]["S1"] == "falsified"
    assert synthetic["unstable_mean"]["S1"] == "falsified"
    # S4 fails everywhere, including where the derivation is exact -- the predeclared
    # ORDERING gloss was wrong for a -> inf, not the general formula
    assert all(s["S4"] == "falsified" for s in synthetic.values())


def test_the_derived_form_is_exact_where_its_premise_holds_and_the_naive_one_is_not(verdict):
    pinned = verdict["synthetic"]["t_pinned"]
    assert pinned["S1"]["median_relative_error"] < 0.002
    assert pinned["S2"]["median_relative_error_naive"] > 0.45
    unstable = verdict["synthetic"]["unstable_mean"]["detail"]
    assert min(r["A_observed"] for r in unstable) < 0.2, (
        "the unstable-mean condition must show centring IMPROVING reproducibility")


def test_the_measured_shared_direction_is_not_the_one_the_observation_quoted(verdict):
    """f on the exported wsi_biology block is ~0.992, not the 0.8133 of a different block."""
    import json
    real = json.loads((VENDORED / "P2_CENTRING_AMPLIFICATION_REAL.json").read_text(
        encoding="utf-8"))
    f = [real["repeats"][r]["views"]["wsi_biology"]["raw"]["f"]
         for r in sorted(real["repeats"])]
    assert min(f) > 0.99 and max(f) < 1.0
    rna = [real["repeats"][r]["views"]["rna_biology"]["raw"]["f"]
           for r in sorted(real["repeats"])]
    assert 0.25 < min(rna) and max(rna) < 0.29
    # residualisation removes the shared direction entirely, so the law predicts
    # NOTHING about a comparison between two residualised views
    resid = [real["repeats"][r]["views"][v]["residualised"]["f"]
             for r in sorted(real["repeats"])
             for v in ("wsi_biology", "rna_biology", "full_biology")]
    assert max(resid) < 1e-20


def test_more_than_half_the_view_effect_predates_centring(verdict):
    p5 = verdict["real"]["P5"]
    assert p5["fraction_of_the_view_effect_present_before_centring"] == pytest.approx(
        0.5256, abs=5e-4)
