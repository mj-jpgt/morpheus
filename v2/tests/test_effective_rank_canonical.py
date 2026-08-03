"""One definition of effective rank, everywhere, checked three ways.

This repository shipped THREE mutually incompatible statistics under the name
`effective_rank` across nine call sites (`paper/P2_RANK_DRAFT.md` §3.1 counted
three of them; a full scan while writing these tests found six more). The
historical values quoted in that paper were therefore not all measured with the
same statistic. These tests exist so that cannot recur:

1. `test_hand_computed_*` pin the canonical function to a value derivable with
   pencil and paper from a matrix whose spectrum is known exactly.
2. `test_every_call_site_resolves_to_the_one_function` asserts object identity
   across the importable call sites.
3. `test_no_second_definition_exists_in_the_tree` walks the source tree and
   fails on any SVD-based rank computed outside `calibra/spectral.py`.

The canonical definition and the reasons for each of its three choices are in
`v2/calibra/spectral.py`; do not restate them here, restate them in the paper.
"""
from __future__ import annotations

import ast
import math
from pathlib import Path

import numpy as np
import pytest
import torch

from morpheus.v2.calibra import spectral
from morpheus.v2.calibra.spectral import CANONICAL, RANK_VARIANTS, effective_rank

REPO = Path(spectral.__file__).resolve().parents[2]  # .../v2/calibra/spectral.py -> repo root

# A 6x3 matrix with exactly zero column means whose columns are mutually
# orthogonal with norms 2*sqrt(2), sqrt(2), sqrt(2) -- so its singular values are
# exactly proportional to (2, 1, 1).
CENTRED_KNOWN_SPECTRUM = np.array(
    [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0],
     [0.0, 1.0, 0.0], [0.0, -1.0, 0.0],
     [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])

# For sigma proportional to (2, 1, 1): p = (1/2, 1/4, 1/4),
#   H = -(1/2)ln(1/2) - 2*(1/4)ln(1/4) = (1/2)ln2 + (1/2)ln4 = (1/2)ln8
#   erank = exp(H) = sqrt(8) = 2*sqrt(2)
# and the order-2 Hill number is 1/sum(p^2) = 1/(1/4 + 1/16 + 1/16) = 8/3.
ERANK_211 = 2.0 * math.sqrt(2.0)          # 2.8284271247461903
ORDER2_211 = 8.0 / 3.0                    # 2.6666666666666665


def test_hand_computed_value_on_a_known_spectrum_uncentred():
    """diag(2, 1, 1) has singular values exactly (2, 1, 1)."""
    matrix = np.diag([2.0, 1.0, 1.0])
    assert effective_rank(matrix, centre=False) == pytest.approx(ERANK_211, rel=1e-12)
    assert effective_rank(matrix, centre=False, order=2) == pytest.approx(ORDER2_211, rel=1e-12)


def test_hand_computed_value_on_a_known_spectrum_centred():
    """The canonical path recovers the same value through a non-zero column mean.

    Adding a constant row offset changes the raw spectrum but not the centred one,
    so this simultaneously pins the value AND proves centring is actually applied.
    """
    offset = np.array([5.0, -3.0, 7.0])
    shifted = CENTRED_KNOWN_SPECTRUM + offset
    assert effective_rank(shifted) == pytest.approx(ERANK_211, rel=1e-12)
    assert effective_rank(shifted, variant=RANK_VARIANTS["R2"]) == pytest.approx(ORDER2_211, rel=1e-12)
    # ... and the uncentred value is a long way off (1.68 against 2.83 on this
    # matrix), which is the whole point of having to state which one was used.
    assert effective_rank(shifted, centre=False) < 0.7 * ERANK_211


def test_roy_vetterli_property_1_bounds():
    """1 <= erank(A) <= rank(A) <= min(M, N), with both equality cases."""
    rng = np.random.default_rng(0)
    for shape in [(40, 12), (12, 40), (200, 200)]:
        x = rng.normal(size=shape)
        value = effective_rank(x, centre=False)
        assert 1.0 <= value <= np.linalg.matrix_rank(x) <= min(shape)
    # equality below: all mass on one singular value
    rank_one = np.outer(np.arange(1.0, 9.0), np.arange(1.0, 5.0))
    assert effective_rank(rank_one, centre=False) == pytest.approx(1.0, abs=1e-9)
    # equality above: a flat spectrum reaches the algebraic rank exactly
    assert effective_rank(np.eye(7), centre=False) == pytest.approx(7.0, rel=1e-12)


def test_order_two_is_never_above_order_one():
    """R2/R3 are order-2 Hill numbers and are <= R1 for EVERY matrix.

    This is why the historical R2 and R3 values cannot be compared with the R1
    values in the same table: the difference has a fixed sign.
    """
    rng = np.random.default_rng(7)
    for _ in range(25):
        x = rng.normal(size=(60, 24)) @ np.diag(rng.random(24) ** 3)
        assert effective_rank(x, order=2) <= effective_rank(x, order=1) + 1e-9


def test_scale_invariance():
    """erank(cA) == erank(A) for c != 0 -- Roy & Vetterli Property 2.

    The previous implementation used an ABSOLUTE 1e-12 singular-value cut and
    failed this: scaling a matrix down by 1e-9 silently emptied the spectrum.
    """
    rng = np.random.default_rng(11)
    x = rng.normal(size=(80, 30))
    base = effective_rank(x)
    for scale in (1e-9, 1e-3, 1e3, 1e9):
        assert effective_rank(x * scale) == pytest.approx(base, rel=1e-9)


def test_degenerate_inputs_score_zero_at_any_scale():
    assert effective_rank(np.ones((50, 8))) == 0.0
    assert effective_rank(np.ones((50, 8)) * 1e12) == 0.0
    assert effective_rank(torch.ones(16, 256)) == 0.0
    assert effective_rank(np.zeros((5, 5))) == 0.0


def test_torch_and_numpy_backends_agree():
    """The training tripwire uses the on-device path; it must be the same number."""
    rng = np.random.default_rng(3)
    x = rng.normal(size=(64, 96)) @ np.diag(rng.random(96) ** 2)
    for name, variant in RANK_VARIANTS.items():
        numpy_value = effective_rank(x, variant=variant)
        torch_value = effective_rank(torch.from_numpy(x), variant=variant)
        assert torch_value == pytest.approx(numpy_value, rel=1e-9), name


def test_row_normalisation_is_a_real_difference_not_a_rescaling():
    """R3's row L2-normalisation changes the answer; it is not a no-op."""
    rng = np.random.default_rng(5)
    x = rng.normal(size=(50, 20)) * rng.lognormal(sigma=2.0, size=(50, 1))
    assert effective_rank(x) != pytest.approx(effective_rank(x, normalise_rows=True), rel=1e-3)


def test_uncentred_rank_is_a_function_of_the_mean_and_centred_rank_is_not():
    """The reason the canonical variant deviates from Roy & Vetterli by centring.

    Uncentred effective rank moves with the magnitude of the column mean relative
    to the spread, in BOTH directions, so it is not a property of the spread at
    all. Centred effective rank is exactly invariant to a shift.
    """
    rng = np.random.default_rng(13)
    spread = rng.normal(size=(128, 64))
    mean = rng.normal(size=64)
    centred_value = effective_rank(spread)

    # (a) a large shared offset makes an isotropic, near-full-rank representation
    #     read as collapsed when it is not: the mean direction takes almost all
    #     the spectral mass.
    assert effective_rank(spread + 2000.0 * mean, centre=False) < 1.05
    assert centred_value > 40.0

    # (b) on the collapse family z_i = m + a_i*u with m NOT parallel to u and of
    #     comparable magnitude, uncentred rank reads ~2 where the representation
    #     has exactly one direction of variation.
    direction = rng.normal(size=64)
    direction /= np.linalg.norm(direction)
    offset = rng.normal(size=64)
    offset -= offset @ direction * direction          # orthogonal to u
    offset /= np.linalg.norm(offset)
    collapsed = offset + np.outer(rng.normal(size=128), direction)
    assert effective_rank(collapsed) == pytest.approx(1.0, abs=1e-6)
    assert effective_rank(collapsed, centre=False) > 1.4

    # (c) the canonical value does not move under any shift; the uncentred one does.
    for shift in (0.0, 1.0, 50.0, 5000.0):
        assert effective_rank(spread + shift * mean) == pytest.approx(centred_value, rel=1e-9)
    assert effective_rank(spread + 5000.0 * mean, centre=False) != pytest.approx(
        effective_rank(spread, centre=False), rel=1e-2)


def test_variant_registry_matches_the_three_historical_statistics():
    assert RANK_VARIANTS["R1"] == CANONICAL == spectral.RankVariant(True, False, 1)
    assert RANK_VARIANTS["R2"] == spectral.RankVariant(True, False, 2)   # d1_audit.py
    assert RANK_VARIANTS["R3"] == spectral.RankVariant(True, True, 2)    # geometry probes
    assert CANONICAL.label == "centred|order1"


def test_every_call_site_resolves_to_the_one_function():
    """Object identity, not "looks the same"."""
    import morpheus.v2.calibra as calibra_package
    import morpheus.v2.calibra.e1_rank_information as e1
    import morpheus.v2.calibra.e2_expressible_intersection as e2
    import morpheus.v2.calibra.run_calibra as run_calibra
    import morpheus.v2.runner as runner
    import morpheus.v2.tests.test_stress_collapse as stress
    import morpheus.v2.training as training

    for module in (calibra_package, e1, e2, run_calibra, runner, stress, training):
        assert module.effective_rank is effective_rank, module.__name__

    e0 = pytest.importorskip("morpheus.v2.calibra.e0_basis_transfer")
    assert e0.effective_rank is effective_rank
    # the alias must delegate, not reimplement
    rng = np.random.default_rng(17)
    sample = torch.from_numpy(rng.normal(size=(40, 16)))
    assert e0._effective_rank_torch(sample) == pytest.approx(effective_rank(sample), rel=1e-12)


#: Files allowed to compute an SVD. Anything else that does is either a new
#: duplicate of effective rank or a statistic that needs its own name and test.
#: `spectral.py` is the one rank implementation; the rest use SVDs for CCA,
#: whitening, PCA bases, subspace overlap or a hard numerical rank.
SVD_ALLOWLIST = {
    "v2/calibra/spectral.py",                    # THE effective-rank implementation
    "v2/calibra/e0_basis_transfer.py",           # subspace overlap + numerical rank
    "v2/calibra/e1_rank_information.py",         # stable_rank (named differently), CCA
    "v2/calibra/e2_expressible_intersection.py", # _svd_rank (hard rank), CCA
    "v2/calibra/induced_correlation_sweep.py",   # design-matrix conditioning
    "v2/discovery_targets.py",                   # PCA basis
    "v2/pbs.py",                                 # perturbation-basis construction
    "v2/research/rebase/d1_geometry_probe.py",   # torch.linalg.matrix_rank, a HARD rank
    "v2/calibra/claim_guards.py",                # the token appears only in prose
    "tests/test_claim_guards.py",                # ditto
    "v2/tests/test_effective_rank_canonical.py", # this file, which names the tokens
}

RANK_SHAPED = ("svdvals", "linalg.svd")


def test_no_second_definition_exists_in_the_tree():
    """No SVD-based rank statistic outside `spectral.py`, and no `def effective_rank`.

    The failure this guards against is not hypothetical: it is exactly what
    `paper/P2_RANK_DRAFT.md` §3.1 reports as a finding about this repository.
    """
    offenders_svd: list[str] = []
    offenders_def: list[str] = []
    for path in sorted(REPO.rglob("*.py")):
        relative = path.relative_to(REPO).as_posix()
        if "__pycache__" in relative or relative.startswith(("runs/", "pytmp/")):
            continue
        source = path.read_text(encoding="utf-8", errors="ignore")
        if any(token in source for token in RANK_SHAPED) and relative not in SVD_ALLOWLIST:
            offenders_svd.append(relative)
        try:
            tree = ast.parse(source, filename=relative)
        except SyntaxError:  # vendored / generated sources are not ours to police
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef) or "effective_rank" not in node.name:
                continue
            if node.name.startswith("test_"):  # tests OF the function, not copies of it
                continue
            if relative == "v2/calibra/spectral.py":
                continue
            body = [s for s in node.body if not (isinstance(s, ast.Expr)
                                                 and isinstance(s.value, ast.Constant))]
            delegates = (len(body) == 1 and isinstance(body[0], ast.Return)
                         and isinstance(body[0].value, ast.Call)
                         and getattr(body[0].value.func, "id", None) == "effective_rank")
            if not delegates:
                offenders_def.append(f"{relative}:{node.lineno}:{node.name}")

    assert not offenders_svd, (
        "SVD-based rank computed outside calibra/spectral.py -- add it to the "
        f"allowlist only after confirming it is not effective rank: {offenders_svd}")
    assert not offenders_def, (
        "a second definition of effective_rank has reappeared; the whole point of "
        f"this module is that there is exactly one: {offenders_def}")


def test_the_deviating_call_sites_name_their_deviation():
    """The tripwire and the gate probe use R3. That must be visible in the source."""
    for relative in ("v2/training.py", "v2/runner.py"):
        source = (REPO / relative).read_text(encoding="utf-8")
        assert 'RANK_VARIANTS["R3"]' in source, (
            f"{relative} calls effective_rank for a threshold calibrated on R3; the "
            "variant must be named at the call site, never left to the default")
