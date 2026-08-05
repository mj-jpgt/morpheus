"""Ledger-accounting and stage-order tests for the P5 pilot funnel.

These are mechanics tests only -- they check the pipeline's bookkeeping (a cell can
never survive a later stage without having survived every earlier one; BH-FDR is
applied once across the whole stage-2 test set; stage 4 only runs on stage-3
survivors) on tiny hand-built arrays, not on `build_synthetic_data`'s full-size
ladder (which is intentionally too slow for a unit test). See
`NOTEBOOK_ENTRIES/p5_pilot_funnel_synthetic_dry_run_20260805T...` for the actual
(synthetic, non-biological) pilot run this module supports.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.p5_pilot_funnel import (
    enumerate_candidate_space, ledger, ledger_rows, run_funnel, stage1_coarse_filter,
    stage2_certify, stage3_multiplicity_correction, synthetic_confusion,
)


def _toy_data(seed: int = 0, n: int = 60, n_axes: int = 4, n_sites: int = 3):
    """Two strata, three targets each, small enough to run in well under a second."""
    rng = np.random.default_rng(seed)
    strata = ["A", "B"]
    target_ids = ["T0", "T1", "T2"]
    data = {}
    for stratum in strata:
        site = rng.integers(0, n_sites, size=n).astype(str)
        x = rng.normal(size=(n, n_axes))
        half = n // 2
        order = rng.permutation(n)
        disc, repl = order[:half], order[half:]
        y = {t: rng.normal(size=n) for t in target_ids}
        data[stratum] = {
            "x_discovery": x[disc], "site_discovery": site[disc],
            "y_discovery": {t: v[disc] for t, v in y.items()},
            "x_replication": x[repl], "site_replication": site[repl],
            "y_replication": {t: v[repl] for t, v in y.items()},
        }
    return strata, target_ids, data


def test_stage_order_is_a_strict_funnel() -> None:
    """No cell may reach a later stage without clearing every earlier one."""
    strata, target_ids, data = _toy_data()
    results = run_funnel(strata, target_ids, data, seed=0)
    for r in results:
        if r.stage2_attempted:
            assert r.stage1_pass
        if r.stage3_survivor:
            assert r.stage2_attempted
        if r.stage4_attempted:
            assert r.stage3_survivor
        if r.stage4_replicated:
            assert r.stage4_attempted


def test_ledger_counts_are_monotone_non_increasing() -> None:
    strata, target_ids, data = _toy_data()
    results = run_funnel(strata, target_ids, data, seed=0)
    summary = ledger(results)
    assert summary["n_candidates_entered_stage0"] == len(strata) * len(target_ids)
    assert summary["n_survived_stage1_coarse_filter"] <= summary["n_candidates_entered_stage0"]
    assert summary["n_reached_stage2_certify"] <= summary["n_survived_stage1_coarse_filter"]
    assert summary["n_cleared_stage3_bh_fdr"] <= summary["n_reached_stage2_certify"]
    assert summary["n_attempted_stage4_replication"] <= summary["n_cleared_stage3_bh_fdr"]
    assert summary["n_replicated_stage4"] <= summary["n_attempted_stage4_replication"]


def test_stage1_keep_fraction_is_honoured() -> None:
    strata, target_ids, data = _toy_data()
    cells = enumerate_candidate_space(strata, target_ids)
    results = stage1_coarse_filter(cells, data, keep_fraction=0.5)
    assert sum(r.stage1_pass for r in results) == round(len(cells) * 0.5)


def test_bh_fdr_never_more_liberal_than_the_raw_p_values() -> None:
    """A regression guard on the BH wiring: with all-null (shuffled) data at a
    demanding q, the survivor set must never exceed the raw-p<=q count (BH is a
    strictly more conservative filter than an uncorrected threshold at the same q)."""
    strata, target_ids, data = _toy_data(seed=3, n=80)
    cells = enumerate_candidate_space(strata, target_ids)
    results = stage1_coarse_filter(cells, data, keep_fraction=1.0)  # force every cell to stage 2
    for r in results:
        stage2_certify(r, data, seed=3)
    stage3_multiplicity_correction(results, q=0.10)
    raw_below_q = sum(1 for r in results if r.stage2_attempted and np.isfinite(r.permutation_p)
                      and r.permutation_p <= 0.10)
    bh_survivors = sum(r.stage3_survivor for r in results)
    assert bh_survivors <= raw_below_q


def test_ledger_rows_frame_has_one_row_per_candidate_cell() -> None:
    strata, target_ids, data = _toy_data()
    results = run_funnel(strata, target_ids, data, seed=0)
    frame = ledger_rows(results)
    assert len(frame) == len(strata) * len(target_ids)
    assert set(frame["stratum"]) <= set(strata)


def test_synthetic_confusion_reports_planted_and_null_counts_separately() -> None:
    strata, target_ids, data = _toy_data()
    planted = {("A", "T0")}
    cells = enumerate_candidate_space(strata, target_ids, planted=planted)
    results = stage1_coarse_filter(cells, data, keep_fraction=1.0)
    confusion = synthetic_confusion(results)
    assert confusion["n_planted"] == 1
    assert confusion["n_null"] == len(cells) - 1


def test_main_refuses_real_data_flags_with_not_implemented(monkeypatch) -> None:
    """The real-data path is specified but was never exercised (no reachable
    artifact) -- it must fail loudly through the actual CLI, not silently fall
    through to synthetic data or crash on a missing file first.

    Uses a plain string output path rather than the ``tmp_path`` fixture: this
    checkout's default pytest temp root is not writable in this environment
    (``PermissionError`` on ``AppData\\Local\\Temp\\pytest-of-*``), a known,
    environment-specific hazard unrelated to this test's own logic. ``main``
    raises before ever touching ``--output``, so no directory needs to exist.
    """
    from morpheus.v2 import p5_pilot_funnel

    monkeypatch.setattr(
        "sys.argv",
        ["p5_pilot_funnel", "--artifact", "does_not_exist.npz", "--targets", "also_missing.npz",
         "--output", "unused_output_dir"],
    )
    with pytest.raises(NotImplementedError):
        p5_pilot_funnel.main()


def test_main_requires_a_data_source(monkeypatch) -> None:
    """Neither --synthetic-dry-run nor a real-data pair supplied -> argparse error, not a crash."""
    from morpheus.v2 import p5_pilot_funnel

    monkeypatch.setattr("sys.argv", ["p5_pilot_funnel", "--output", "unused_output_dir"])
    with pytest.raises(SystemExit):
        p5_pilot_funnel.main()
