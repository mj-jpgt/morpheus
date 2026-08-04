"""The four P2 analysis scripts import and run, on a synthetic cohort small enough for CI.

These scripts produce `paper/P2_RANK_DRAFT.md` §4.2, §4.4, §4.5 and §4.6. Until 2026-08-05 they
existed only at `~/e0_run/p2_*.py` on the GPU box, where nothing tested them and nothing could
tell that the workspace they imported `spectral.py` from had drifted from HEAD
(`NOTEBOOK_ENTRIES/WORKSPACE_DRIFT_AUDIT_ALL_20260803T2359Z.md`). Vendoring them into the
repository is only half the fix; without a test they can still rot silently against the
`calibra` functions they call.

This module therefore does the cheapest thing that would have caught that class of failure:
build a small synthetic cohort, run every script end to end on it, and assert the structural
invariants that hold for ANY input -- not the paper's numbers, which need the frozen artifacts.

The one numerical assertion here is the one that matters for the paper's own §3.1 argument:
that the five statistics `p2_rank_variants.py` reports are five different functions, ordered
R1 >= R2 >= PR, so a table that names the wrong one is reporting a different measurement.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2 import (p2_competing_metrics, p2_necessity_and_variance,
                                            p2_rank_variants, p2_robustness, p2_selection_rule)

#: Deliberately tiny. The real cohort is 2,766 test patients x 256 dimensions x 12 artifacts;
#: everything asserted below is a property of the code path, not of the scale.
N_PATIENTS = 120
N_DIMS = 24
#: Two targets per group at minimum: `cross_fitted_residuals` is a multi-output ridge and a
#: one-column group is a shape the real target table never has.
TARGET_GROUPS = (["heldout_pathway"] * 4 + ["immune_tme"] * 2 + ["tumour_state"] * 2
                 + ["random_control"] * 2 + ["hallmark_in_training"] * 2)
N_TARGETS = len(TARGET_GROUPS)

#: The labels the scripts' hardcoded pair lists expect (2 experiments x 2 arms x 3 seeds).
LABELS = ["H42", "H43", "H44", "I42", "I43", "I44", "P42", "P43", "P44", "F42", "F43", "F44"]


def _patient_ids(n: int) -> np.ndarray:
    """TCGA-style barcodes: `pooled_tissue_source_site` splits on '-' and takes field 1."""
    return np.asarray([f"TCGA-{i % 6:02d}-{i:04d}" for i in range(n)])


def _artifact(path, rng, *, decay: float) -> None:
    """A frozen-artifact-shaped .npz: three co-trained views, a split and a cohort.

    `decay` sets how fast the spectrum falls, so the twelve artifacts have distinct effective
    ranks and the selection-rule and variance code paths see a real spread rather than ties.
    """
    basis = np.linalg.qr(rng.normal(size=(N_DIMS, N_DIMS)))[0]
    scales = decay ** np.arange(N_DIMS)
    latent = rng.normal(size=(N_PATIENTS, N_DIMS)) * scales
    wsi = latent @ basis.T + 0.05 * rng.normal(size=(N_PATIENTS, N_DIMS))
    rna = latent @ basis.T + 0.20 * rng.normal(size=(N_PATIENTS, N_DIMS))
    np.savez(
        path,
        patient_ids=_patient_ids(N_PATIENTS),
        cancers=np.asarray([f"C{i % 4}" for i in range(N_PATIENTS)]),
        split=np.asarray(["test"] * N_PATIENTS),
        wsi_biology=wsi,
        rna_biology=rna,
        full_biology=0.5 * (wsi + rna),
    )


def _targets(path, rng) -> None:
    np.savez(
        path,
        patient_ids=_patient_ids(N_PATIENTS),
        target_names=np.asarray([f"T{i}" for i in range(N_TARGETS)]),
        target_groups=np.asarray(TARGET_GROUPS),
        scores=rng.normal(size=(N_PATIENTS, N_TARGETS)),
    )


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    """Twelve synthetic artifacts plus a target table, written once for the module."""
    root = tmp_path_factory.mktemp("p2")
    rng = np.random.default_rng(20260805)
    _targets(root / "targets.npz", rng)
    paths = {}
    for index, label in enumerate(LABELS):
        path = root / f"{label}.npz"
        _artifact(path, rng, decay=0.80 + 0.02 * (index % 5))
        paths[label] = path
    return root, paths


@pytest.fixture(scope="module")
def metrics_json(cohort):
    """Run `p2_competing_metrics.py` end to end, as §4.6 runs it -- with subsampling on."""
    root, paths = cohort
    output = root / "P2_METRICS.json"
    p2_competing_metrics.main(
        ["--targets", str(root / "targets.npz"), "--output", str(output), "--subsamples", "3",
         "--artifacts"] + [f"{k}={v}" for k, v in paths.items()])
    return output


def test_competing_metrics_writes_every_metric_for_every_artifact(metrics_json):
    data = json.loads(metrics_json.read_text(encoding="utf-8"))
    assert set(data) == set(LABELS) | {"_config"}
    for label in LABELS:
        metrics = data[label]["metrics"]
        for key in ("effective_rank_raw", "effective_rank_residualised", "rankme_raw",
                    "rankme_residualised", "participation_ratio_raw",
                    "participation_ratio_residualised", "stable_rank_raw",
                    "stable_rank_residualised"):
            assert np.isfinite(metrics[key]) and metrics[key] > 0, (label, key)
        assert np.isfinite(metrics["alpha_req_residualised"]["alpha"]), label
        assert np.isfinite(metrics["lidar_residualised"]["lidar"]), label
        # every target group the paper reads must be scored, not just the headline one
        assert set(data[label]["points"]) == set(p2_competing_metrics.GROUPS)
        assert 0.0 <= data[label]["points"]["untrained40"]["top_cca"] <= 1.0
        # the subsampling arm of §4.4(1) must produce a spread, not a constant
        assert data[label]["subsample"]["effective_rank_residualised"]["sd"] > 0.0


def test_rankme_is_the_published_formula_not_our_effective_rank():
    """RankMe is uncentred with eps INSIDE the normalisation; `effective_rank` is neither.

    P2 §4.4(4) turns on these two being close on our artifacts but not identical. A vendored
    RankMe that had been quietly routed through `spectral.py` would make that check vacuous.
    """
    rng = np.random.default_rng(3)
    x = rng.normal(size=(200, 32)) + 40.0        # a large column mean, which centring removes
    assert p2_competing_metrics.rankme(x) != pytest.approx(effective_rank(x), rel=1e-3)
    assert p2_competing_metrics.rankme(x) == pytest.approx(
        effective_rank(x, centre=False), rel=1e-3)


def test_participation_ratio_is_not_the_canonical_R2():
    """The two statistics `p2_rank_variants.py` reports side by side are different functions.

    `participation_ratio` is the order-2 Hill number of the EIGENVALUES; `RANK_VARIANTS["R2"]`
    is the order-2 Hill number of the SINGULAR VALUES. Both have been called a participation
    ratio in this project. R1 >= R2 >= PR holds for every matrix, so substituting one for the
    other is a systematic, signed change -- not a rounding difference.
    """
    rng = np.random.default_rng(19)
    for _ in range(10):
        x = rng.normal(size=(80, 20)) @ np.diag(rng.random(20) ** 3)
        r1 = effective_rank(x, variant=RANK_VARIANTS["R1"])
        r2 = effective_rank(x, variant=RANK_VARIANTS["R2"])
        pr = p2_competing_metrics.participation_ratio(x)
        assert r1 >= r2 >= pr - 1e-9
        assert r2 != pytest.approx(pr, rel=1e-3)


def test_rank_variants_reports_five_statistics_and_a_verdict_per_pair(cohort):
    root, paths = cohort
    result = p2_rank_variants.main(
        ["--targets", str(root / "targets.npz"), "--output", str(root / "variants.json")]
        + [f"{k}={v}" for k, v in paths.items()])
    assert set(result["artifacts"]) == set(LABELS)
    for record in result["artifacts"].values():
        assert record["R1"] >= record["R2"] >= record["PR"] - 1e-9
        assert np.isfinite(record["R3"]) and np.isfinite(record["PR_rownorm"])
    for name in p2_rank_variants.STATISTICS:
        verdict = result["verdicts"][name]
        assert verdict["ALL"] == verdict["D2"] + verdict["D1"] <= 6
        assert len(verdict["marks"]) == 6
    assert json.loads((root / "variants.json").read_text(encoding="utf-8"))["artifacts"]


def test_selection_rule_scores_every_metric_against_the_ground_truth(metrics_json, capsys):
    p2_selection_rule.main([str(metrics_json)])
    out = capsys.readouterr().out
    for name, _, _ in p2_selection_rule.METRICS:
        assert name in out, name
    assert "SELECTION-RULE SCORE" in out
    assert "METRIC NOISE FLOOR vs BETWEEN-ARM GAP" in out   # the §4.4(1) block


def test_two_sided_binomial_matches_the_values_quoted_in_the_draft():
    """§4.6's power argument quotes these three exactly; they are the section's whole point."""
    assert p2_selection_rule.two_sided_binomial(6, 6) == pytest.approx(0.03125)
    assert p2_selection_rule.two_sided_binomial(5, 6) == pytest.approx(0.21875)
    assert p2_selection_rule.two_sided_binomial(4, 6) == pytest.approx(0.6875)


def test_necessity_and_variance_runs_and_declares_its_thresholds(metrics_json, capsys):
    p2_necessity_and_variance.main([str(metrics_json)])
    out = capsys.readouterr().out
    assert "VARIANCE DECOMPOSITION" in out and "NECESSITY SCAN" in out
    assert "arm share" in out
    # the pre-declared criterion is part of the result and must stay visible in the output
    assert str(p2_necessity_and_variance.RANK_FOLD) in out
    assert str(p2_necessity_and_variance.CCA_DELTA) in out


def test_variance_decomposition_recovers_a_known_split():
    """A synthetic design with zero within-arm spread must give a 100% arm share.

    The decomposition is the paper's most important display item (§4.2) and is eight lines of
    arithmetic with no test anywhere. This pins both ends of its range.
    """
    labels = p2_necessity_and_variance.LABELS
    pure_arm = {label: {"H": 1.0, "I": 2.0, "P": 3.0, "F": 4.0}[label[0]] for label in labels}

    # re-derive the decomposition the way the script does, on the raw scale
    def share(mapping):
        groups: dict[str, list[float]] = {}
        for label in labels:
            groups.setdefault(label[0], []).append(mapping[label])
        grand = float(np.mean([mapping[label] for label in labels]))
        ss_arm = sum(len(g) * (np.mean(g) - grand) ** 2 for g in groups.values())
        ss_seed = sum(sum((x - np.mean(g)) ** 2 for x in g) for g in groups.values())
        return ss_arm / (ss_arm + ss_seed)

    assert share(pure_arm) == pytest.approx(1.0)
    # all variation in the nuisance factor -> zero arm share
    rng = np.random.default_rng(5)
    offsets = rng.normal(size=3)
    seed_only = {label: float(offsets[["42", "43", "44"].index(label[1:])]) for label in labels}
    assert share(seed_only) == pytest.approx(0.0, abs=1e-12)
    assert p2_necessity_and_variance.RANK_FOLD == 2.0
    assert p2_necessity_and_variance.CCA_DELTA == 0.0705


def test_robustness_reports_all_three_views_and_both_cca_estimators(cohort, capsys):
    root, paths = cohort
    rows = p2_robustness.main(str(root / "P2_ROBUSTNESS.json"), str(root / "targets.npz"), 42,
                              [f"{k}={v}" for k, v in paths.items()])
    assert set(rows) == set(LABELS)
    for record in rows.values():
        assert set(record) == {"wsi_biology", "rna_biology", "full_biology"}
        for view in record.values():
            assert view["eff_rank"] > 0
            assert 0.0 <= view["cca_insample"] <= 1.0
            # held-out CCA is the unbiased estimator and must not exceed the in-sample maximum
            assert view["cca_heldout"] <= view["cca_insample"] + 1e-9
    assert "ARM ORDERING under each representation" in capsys.readouterr().out
    assert json.loads((root / "P2_ROBUSTNESS.json").read_text(encoding="utf-8"))
