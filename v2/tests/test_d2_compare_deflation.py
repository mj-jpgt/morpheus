"""Guards for the proliferation-deflation path in d2_compare.

The failure mode being defended against is the one this project keeps hitting: a
term that is silently off. A deflation that quietly does nothing would return
"the gap survives proliferation removal" for *any* data, which is a result that
cannot be wrong and is therefore worthless.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.residualise import confound_design, cross_fitted_residuals
from morpheus.v2.calibra.spectral import top_canonical_correlation
from morpheus.v2.research.rebase.d2_compare import load_deflation_block


def _axis_npz(tmp_path, n=120, k=6, seed=0):
    rng = np.random.default_rng(seed)
    path = tmp_path / "axes.npz"
    scores = rng.normal(size=(n, k))
    np.savez(path, patient_ids=np.asarray([f"TCGA-AA-{i:04d}" for i in range(n)]),
             target_names=np.asarray([f"PBS_{j:03d}" for j in range(k)]), scores=scores)
    return path, scores


def test_deflation_block_rejects_the_ways_a_cut_can_be_silently_wrong(tmp_path):
    path, _ = _axis_npz(tmp_path)
    good = tmp_path / "good.txt"
    good.write_text("PBS_000\nPBS_002\n")
    index, columns, names = load_deflation_block(str(path), str(good))
    assert names == ["PBS_000", "PBS_002"] and columns.shape == (120, 2)
    assert index["TCGA-AA-0000"] == 0

    empty = tmp_path / "empty.txt"; empty.write_text("\n  \n")
    with pytest.raises(ValueError, match="is empty"):
        load_deflation_block(str(path), str(empty))

    dupe = tmp_path / "dupe.txt"; dupe.write_text("PBS_000\nPBS_000\n")
    with pytest.raises(ValueError, match="lists an axis twice"):
        load_deflation_block(str(path), str(dupe))

    absent = tmp_path / "absent.txt"; absent.write_text("PBS_000\nPBS_999\n")
    with pytest.raises(ValueError, match="absent from"):
        load_deflation_block(str(path), str(absent))


def test_a_constant_axis_is_rejected_rather_than_residualised_as_a_no_op(tmp_path):
    """Regressing out a constant column removes nothing but reports a cut was made."""
    path = tmp_path / "axes.npz"
    scores = np.random.default_rng(1).normal(size=(80, 3))
    scores[:, 1] = 4.0
    np.savez(path, patient_ids=np.asarray([f"TCGA-AA-{i:04d}" for i in range(80)]),
             target_names=np.asarray(["PBS_000", "PBS_001", "PBS_002"]), scores=scores)
    axes = tmp_path / "a.txt"; axes.write_text("PBS_001\n")
    with pytest.raises(ValueError, match="non-finite or constant"):
        load_deflation_block(str(path), str(axes))


def test_deflation_is_live_and_symmetric_across_both_arms(tmp_path):
    """THE test. Build a signal that lives ENTIRELY in the deflated subspace and
    check the canonical correlation collapses for both arms -- if the deflation
    were off, or applied to one arm only, this would not happen."""
    rng = np.random.default_rng(7)
    n = 400
    path, axis_scores = _axis_npz(tmp_path, n=n, k=6, seed=7)
    axes = tmp_path / "cut.txt"; axes.write_text("PBS_000\nPBS_001\n")
    _, columns, names = load_deflation_block(str(path), str(axes))

    # Targets and both arms are driven by axes 0 and 1 and nothing else.
    driver = axis_scores[:, :2]
    y = driver @ rng.normal(size=(2, 5)) + 0.01 * rng.normal(size=(n, 5))
    arm_a = driver @ rng.normal(size=(2, 8)) + 0.01 * rng.normal(size=(n, 8))
    arm_b = driver @ rng.normal(size=(2, 8)) + 0.01 * rng.normal(size=(n, 8))

    import pandas as pd
    cancers = np.asarray(["A"] * (n // 2) + ["B"] * (n // 2))
    base = confound_design(pd.DataFrame({"cancer": cancers}), ["cancer"])
    frame = pd.DataFrame({"cancer": cancers})
    for offset, name in enumerate(names):
        frame[name] = columns[:, offset]
    deflated = confound_design(frame, ["cancer", *names])
    assert deflated.shape[1] == base.shape[1] + 2

    def channel(design, arm):
        return top_canonical_correlation(cross_fitted_residuals(arm, design, seed=0),
                                         cross_fitted_residuals(y, design, seed=0), n_components=3)

    for arm in (arm_a, arm_b):
        assert channel(base, arm) > 0.95, "signal must be present before deflation"
        assert channel(deflated, arm) < 0.5, "deflation is off: the planted subspace survived"


def test_deflation_leaves_an_orthogonal_signal_alone(tmp_path):
    """The mirror of the above: a cut must not be a blunt instrument that removes
    everything. Signal outside the deflated subspace has to survive, or a
    'surviving gap' could never be observed."""
    rng = np.random.default_rng(11)
    n = 400
    path, axis_scores = _axis_npz(tmp_path, n=n, k=6, seed=11)
    axes = tmp_path / "cut.txt"; axes.write_text("PBS_000\nPBS_001\n")
    _, columns, names = load_deflation_block(str(path), str(axes))

    other = axis_scores[:, 3:5]          # untouched by the cut
    y = other @ rng.normal(size=(2, 5)) + 0.01 * rng.normal(size=(n, 5))
    arm = other @ rng.normal(size=(2, 8)) + 0.01 * rng.normal(size=(n, 8))

    import pandas as pd
    cancers = np.asarray(["A"] * (n // 2) + ["B"] * (n // 2))
    frame = pd.DataFrame({"cancer": cancers})
    for offset, name in enumerate(names):
        frame[name] = columns[:, offset]
    deflated = confound_design(frame, ["cancer", *names])
    survived = top_canonical_correlation(cross_fitted_residuals(arm, deflated, seed=0),
                                         cross_fitted_residuals(y, deflated, seed=0), n_components=3)
    assert survived > 0.95, "the cut removed signal it does not span"
