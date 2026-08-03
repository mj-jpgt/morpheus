"""D2.3 guards. The failure mode being defended against is a design that cannot
return the falsifier — i.e. one that reports "legibility is not proliferation"
whatever the data say."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.research.rebase.d2_axis_proliferation import (axis_gene_statistics,
                                                               bimodality_report,
                                                               per_axis_legibility,
                                                               spearman_with_bootstrap_ci,
                                                               verdict)


def _annotations(genes):
    return pd.DataFrame({"gene": genes,
                         "proliferation_loading": [1.0 if i < 10 else 0.0 for i in range(len(genes))],
                         "essentiality_loading": np.linspace(0.0, 1.0, len(genes))})


def test_top_k_statistic_separates_a_proliferation_axis_that_the_weighted_mean_dilutes():
    """The reason both statistics are carried: the weighted mean hides the axis."""
    n_genes = 400
    genes = [f"G{i}" for i in range(n_genes)]
    # A DENSE basis, as an SVD gene basis actually is: the top genes lead, but
    # only just, so |loading| weights are close to uniform. This is the regime
    # the real k=128 basis sits in (its top-100 of 7,072 genes carries ~7% of
    # the |loading| mass against ~1.4% for uniform).
    basis = np.full((n_genes, 2), 0.90)
    basis[:10, 0] = 1.0          # axis 0 leads on the 10 proliferation genes
    basis[200:210, 1] = 1.0      # axis 1 leads on 10 non-proliferation genes
    table = axis_gene_statistics(basis, np.asarray(genes), _annotations(genes), top_k=10)
    # The concentrated statistic separates the two axes completely...
    assert table.loc[0, "prol_top10"] == pytest.approx(1.0)
    assert table.loc[1, "prol_top10"] == pytest.approx(0.0)
    # ...while the |loading|-weighted mean that build_pbs_targets emits leaves
    # BOTH axes pinned to the 10/400 = 0.025 background, i.e. it cannot see the
    # difference at all.
    background = 10 / n_genes
    assert abs(table.loc[0, "prol_wmean"] - background) < 0.01
    assert abs(table.loc[1, "prol_wmean"] - background) < 0.01


def test_permutation_null_is_within_stratum_and_is_not_zero():
    """Capacity alone gives a positive per-axis null; the null must reflect that."""
    rng = np.random.default_rng(0)
    n, p, k = 160, 24, 4
    strata = np.asarray(["A"] * 80 + ["B"] * 80)
    x = rng.normal(size=(n, p))
    y = rng.normal(size=(n, k))          # independent of x by construction
    design = np.zeros((n, 0))
    result = per_axis_legibility(x, y, design, strata, n_permutations=30, seed=1)
    assert result["null_p95"].shape == (k,)
    # A fitted 24-dim direction on 160 rows sees something even against noise.
    assert np.median(result["null_p95"]) > 0.0
    # Unrelated data must not be called legible more often than chance.
    assert int(np.sum(result["observed"] > result["null_p95"])) <= 1


def test_verdict_can_actually_fire_the_falsifier():
    """A design that cannot return the negative result is worthless."""
    frame = pd.DataFrame({
        "legibility": np.concatenate([np.linspace(0.6, 0.9, 32), np.linspace(0.0, 0.1, 96)]),
        "prol_top100": np.concatenate([np.linspace(0.6, 0.9, 32), np.linspace(0.0, 0.1, 96)]),
        "is_legible": np.concatenate([np.ones(32, bool), np.zeros(96, bool)]),
    })
    fired = verdict(frame, "prol_top100")
    assert fired["verdict"] == "FALSIFIER_FIRES_legibility_is_proliferation"
    assert fired["share_of_legible_that_are_loaded"] == pytest.approx(1.0)

    # Same legibility, proliferation deliberately anti-aligned with it.
    clean = frame.copy()
    clean["prol_top100"] = np.concatenate([np.linspace(0.0, 0.1, 32), np.linspace(0.6, 0.9, 96)])
    assert verdict(clean, "prol_top100")["verdict"] == "DISCHARGED_legibility_is_not_proliferation"

    # Too few legible axes is its own outcome, never a discharge.
    thin = frame.copy()
    thin["is_legible"] = np.concatenate([np.ones(5, bool), np.zeros(123, bool)])
    assert verdict(thin, "prol_top100")["verdict"] == "no_verdict_too_few_legible_axes"


def test_bimodality_report_distinguishes_one_cluster_from_two():
    rng = np.random.default_rng(3)
    one = rng.normal(0.0, 1.0, 300)
    two = np.concatenate([rng.normal(-6.0, 0.4, 150), rng.normal(6.0, 0.4, 150)])
    assert bimodality_report(one)["best_n_components"] == 1
    assert bimodality_report(two)["bimodal_by_bic"] is True


def test_spearman_ci_is_resampled_over_axes_and_covers_the_null():
    rng = np.random.default_rng(5)
    a = rng.normal(size=128)
    result = spearman_with_bootstrap_ci(a, rng.normal(size=128), n_boot=400)
    assert result["n_axes"] == 128
    assert result["ci95_low"] < 0.0 < result["ci95_high"]
    assert result["ci_excludes_zero"] is False
    assert spearman_with_bootstrap_ci(a, a, n_boot=400)["ci_excludes_zero"] is True
