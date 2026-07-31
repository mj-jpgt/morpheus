"""Tests for E0 basis transfer. These ARE the crux experiment's self-test.

E0 decides whether an entire research direction survives, so the thing that has
to be tested is not that E0 produces numbers -- it always did -- but that it can
produce the number *zero*. The audit that motivated this file showed it could
not: the only floor was a Haar-random subspace sitting at ``k / n_genes``, while
any two real gene-expression matrices clear that floor from generic
co-expression alone. Every test below that mentions an "arm" exists because of
that defect, and several of them fail against the pre-fix implementation.

All fixtures are small, synthetic and CPU-only; no GPU and no data files.
"""
from __future__ import annotations

import numpy as np
import torch

from morpheus.v2.calibra.e0_basis_transfer import (
    ROTATION_TOLERANCE,
    TransferConfig,
    _arm_result,
    _bootstrap_overlap,
    _decision,
    _energy_arms,
    _full_dictionary_rank,
    _matched_spectrum_rotated_null,
    _overlap,
    _rezscore,
    _right_svd,
    MatrixBundle,
)

DEVICE = torch.device("cpu")
# The production config is q=150/min_q=101 because k reaches 100. The estimator
# is identical here; only the rank budget is shrunk so it fits synthetic data.
CFG = TransferConfig(ks=(5, 10, 20), offsets=(0, 1, 2, 5), primary_offset=1, q=60, min_q=8,
                     draws=40, bootstrap_draws=40)


def _tensor(x: np.ndarray) -> torch.Tensor:
    return torch.as_tensor(np.ascontiguousarray(x, dtype=np.float32), device=DEVICE)


def _module_fixture(*, planted: bool, seed: int = 0, n_arm: int = 90, n_patients: int = 400,
                    n_genes: int = 240, n_modules: int = 60):
    """Two perturbation arms and a tumour matrix over one shared gene-module basis.

    ``planted=False`` is the audit's reproduction case: the arms and the tumour
    matrix share a gene-covariance basis (the synthetic stand-in for ribosomal /
    cell-cycle / MHC co-expression) but NO row-level biology -- every loading
    matrix is drawn independently. Subspace overlap is a function of the gene
    covariance alone, so this construction scores high while containing nothing
    a perturbation did.

    ``planted=True`` gives the responsive arm the tumour-aligned modules and the
    control arm a partly private set. Note that a *weaker* loading on the same
    modules would not be detectable and should not be: principal-angle overlap
    is scale-invariant, so only a difference in direction is a difference at all.
    """
    rng = np.random.default_rng(seed)
    basis = np.linalg.qr(rng.normal(size=(n_genes, n_modules)))[0]
    half = n_modules // 2
    tumour = np.concatenate([np.linspace(3.0, 1.5, half), np.full(n_modules - half, 0.08)])
    generic = tumour.copy()
    # Partly private: the control arm keeps ten of the tumour-aligned modules,
    # so it still clears the Haar floor the way real non-responsive
    # perturbations do, but its leading directions are its own.
    private = np.concatenate([np.linspace(2.9, 2.5, 10), np.zeros(half - 10), np.linspace(3.0, 1.8, n_modules - half)])

    def draw(rows: int, scale: np.ndarray) -> np.ndarray:
        return ((rng.normal(size=(rows, n_modules)) * scale) @ basis.T
                + .05 * rng.normal(size=(rows, n_genes))).astype(np.float32)

    return draw(n_arm, generic), draw(n_arm, private if planted else generic), draw(n_patients, tumour)


def _tumour_reference(t: np.ndarray, cfg: TransferConfig = CFG, seed: int = 7):
    """TCGA-side right singular vectors plus the split-half ceiling, as in `_context`."""
    tt = _rezscore(_tensor(t))
    vt, _ = _right_svd(tt, cfg.q, seed, min_q=cfg.min_q)
    cut = tt.shape[0] // 2
    va, _ = _right_svd(tt[:cut], cfg.q, seed + 1, min_q=cfg.min_q)
    vb, _ = _right_svd(tt[cut:], cfg.q, seed + 2, min_q=cfg.min_q)
    return vt, {k: _overlap(va, vb, k, cfg.primary_offset) for k in cfg.ks}


def _run_arms(responsive: np.ndarray, control: np.ndarray, t: np.ndarray, *, cfg: TransferConfig = CFG, seed: int = 11):
    """Both arms through the identical pipeline, exactly as `_context` runs them."""
    vt, ceiling = _tumour_reference(t, cfg)
    n_genes = t.shape[1]
    responsive_arm, _, _ = _arm_result(_rezscore(_tensor(responsive)), vt, ceiling=ceiling, n_genes=n_genes,
                                       cfg=cfg, draws=cfg.draws, seed=seed, arm="responsive_matched")
    control_arm, _, _ = _arm_result(_rezscore(_tensor(control)), vt, ceiling=ceiling, n_genes=n_genes,
                                    cfg=cfg, draws=cfg.draws, seed=seed + 1_000, arm="nonresponsive")
    return responsive_arm, control_arm


# --- THE AUDIT REPRODUCTION -----------------------------------------------

def test_zero_shared_biology_construction_does_not_pass():
    """THE regression test. Two matrices with no shared row-level biology but a
    shared generic gene-module covariance must NOT be certified as transfer.

    The pre-fix experiment compared the observed overlap only against a Haar
    floor at k/n_genes, so this construction returned "proceed". The fixture
    asserts that old criterion still fires (otherwise the test is vacuous) and
    that the decision rule nonetheless declines."""
    for seed in (3, 5, 13, 17):
        responsive, control, tumour = _module_fixture(planted=False, seed=seed)
        responsive_arm, control_arm = _run_arms(responsive, control, tumour)
        _assert_not_certified(responsive_arm, control_arm, seed)


def _assert_not_certified(responsive_arm, control_arm, seed):
    for k in CFG.ks:
        r = responsive_arm[f"k{k}"]
        # Non-vacuity: the construction is exactly the false positive. It clears
        # the Haar floor by a wide margin while containing no perturbation.
        assert r["beats_haar_floor"], f"seed={seed}: fixture does not reproduce the defect at k={k}"
        assert r["pc1_removed_overlap"] > max(.2, 3 * r["theoretical_random_overlap"]), r["pc1_removed_overlap"]
        # It even reaches the TCGA-self split-half ceiling. Normalising against
        # the ceiling is necessary for interpretation but is still not a control.
        assert r["normalised_alignment"] > .8, r["normalised_alignment"]
        decision = _decision(responsive_arm, control_arm, k)
        assert decision["decision_status"] == "available"
        assert decision["responsive_beats_haar_floor_only"] is True
        assert decision["responsive_exceeds_nonresponsive"] is False, (
            f"seed={seed} k={k}: certified transfer on a construction with zero shared biology "
            f"({decision['responsive_ci95']} vs {decision['nonresponsive_ci95']})")


def test_nonresponsive_control_arm_matches_responsive_when_signal_is_absent():
    """With no planted perturbation structure the arms must be indistinguishable."""
    responsive, control, tumour = _module_fixture(planted=False, seed=5)
    responsive_arm, control_arm = _run_arms(responsive, control, tumour, seed=23)
    for k in CFG.ks:
        r, c = responsive_arm[f"k{k}"], control_arm[f"k{k}"]
        assert max(r["bootstrap_ci95"][0], c["bootstrap_ci95"][0]) <= min(r["bootstrap_ci95"][1], c["bootstrap_ci95"][1]), (
            f"k={k}: bootstrap intervals separated with no signal present: {r['bootstrap_ci95']} vs {c['bootstrap_ci95']}")
        assert _decision(responsive_arm, control_arm, k)["responsive_exceeds_nonresponsive"] is False


def test_responsive_arm_exceeds_control_when_signal_is_planted():
    """The converse: a rule that always says "no" is not a rule. When the
    responsive arm really does carry the tumour-aligned directions and the
    control arm does not, the decision must come back positive -- even though
    the control arm itself still clears the Haar floor, as real data does."""
    for seed in (3, 5, 13, 17):
        responsive, control, tumour = _module_fixture(planted=True, seed=seed)
        responsive_arm, control_arm = _run_arms(responsive, control, tumour, seed=23)
        for k in CFG.ks:
            decision = _decision(responsive_arm, control_arm, k)
            assert decision["overlap_gap"] > 0, (seed, k)
            # k=5 is deliberately excluded from the CI assertion: a 5-dimensional
            # subspace estimated from 90 rows is unstable enough that the two
            # bootstrap intervals can touch even with a real effect. That is a
            # true property of the rule at small k, not a fixture weakness, and
            # it is why the run reports every k rather than one.
            if k >= 10:
                assert decision["responsive_exceeds_nonresponsive"] is True, (
                    f"seed={seed} k={k}: missed a planted effect "
                    f"{decision['responsive_ci95']} vs {decision['nonresponsive_ci95']}")
            # The control arm is a biological floor, not a Haar floor: it sits
            # well above k/n_genes and the rule still separates the arms.
            assert control_arm[f"k{k}"]["pc1_removed_overlap"] > 2 * control_arm[f"k{k}"]["theoretical_random_overlap"]


def test_arms_are_n_matched():
    """The responsive arm must be subsampled to the control arm's row count,
    deterministically. Subspace estimation quality depends on n, so an unmatched
    comparison measures sample size rather than biology."""
    energy = np.concatenate([np.full(400, 1e-6), np.full(37, .9), np.full(50, .2), [np.nan]])
    rows, meta = _energy_arms(energy, cfg=CFG, seed=17)
    assert meta["n_responsive_full"] == 400 and meta["n_nonresponsive"] == 37
    assert meta["n_responsive_matched"] == meta["n_nonresponsive"] == len(rows["responsive_matched"])
    assert meta["arms_are_n_matched"] is True
    assert set(rows["responsive_matched"].tolist()) <= set(rows["responsive_full"].tolist())
    assert not set(rows["responsive_matched"].tolist()) & set(rows["nonresponsive"].tolist())
    assert meta["n_energy_p_nonfinite"] == 1
    repeat, _ = _energy_arms(energy, cfg=CFG, seed=17)
    np.testing.assert_array_equal(rows["responsive_matched"], repeat["responsive_matched"])
    moved, _ = _energy_arms(energy, cfg=CFG, seed=18)
    assert not np.array_equal(rows["responsive_matched"], moved["responsive_matched"])


def test_unmatched_arms_cannot_certify_even_with_a_large_gap():
    """n-matching is part of the rule, not a footnote on it: a responsive arm
    estimated from more rows has a better-estimated subspace, so an unmatched
    win measures sample size."""
    responsive, control, tumour = _module_fixture(planted=True, seed=5, n_arm=180)
    responsive_arm, control_arm = _run_arms(responsive, control[:90], tumour, seed=23)
    for k in CFG.ks:
        decision = _decision(responsive_arm, control_arm, k)
        assert decision["n_rows"] == [180, 90] and decision["arms_are_n_matched"] is False
        assert decision["overlap_gap"] > 0, "fixture has no gap; the test would be vacuous"
        assert decision["responsive_exceeds_nonresponsive"] is False, (
            f"k={k}: certified transfer from arms that are not n-matched")


def test_insufficient_rows_emits_unavailable_not_crash():
    """RPE1's non-responsive arm has 80 rows and cannot support the rank budget.
    That must surface as the repo's unavailable status -- not an exception, and
    above all not a silent pass."""
    _, _, tumour = _module_fixture(planted=False, seed=9)
    vt, ceiling = _tumour_reference(tumour)
    tiny = _rezscore(_tensor(np.random.default_rng(0).normal(size=(20, tumour.shape[1]))))
    arm, vp, sp = _arm_result(tiny, vt, ceiling=ceiling, n_genes=tumour.shape[1], cfg=CFG,
                              draws=CFG.draws, seed=1, arm="nonresponsive")
    assert arm["available"] is False and vp is None and sp is None
    assert arm["metric"] == "status" and np.isnan(arm["value"])
    assert arm["note"].startswith("unavailable_") and "insufficient_rows" in arm["note"]
    assert f"k{CFG.ks[0]}" not in arm
    responsive, control, _ = _module_fixture(planted=True, seed=9)
    responsive_arm, _ = _run_arms(responsive, control, tumour)
    decision = _decision(responsive_arm, arm, CFG.ks[0])
    assert decision["responsive_exceeds_nonresponsive"] is False
    assert decision["decision_status"].startswith("unavailable_")
    assert decision["control_available"] is False


# --- PC1: the question G5/E0.b is named for --------------------------------

def test_pc1_share_detects_trivial_axis_alignment():
    """If all the alignment lives in PC1 you have discovered library size. The
    offset=0 statistic was never computed, so the run could not answer this."""
    rng = np.random.default_rng(4)
    n_genes, n_arm, n_patients = 240, 120, 400
    axis = rng.normal(size=n_genes); axis /= np.linalg.norm(axis)
    p_basis = np.linalg.qr(rng.normal(size=(n_genes, 40)))[0]
    t_basis = np.linalg.qr(rng.normal(size=(n_genes, 40)))[0]
    scale = np.linspace(1.0, .5, 40)
    perturbation = (30.0 * rng.normal(size=(n_arm, 1)) * axis
                    + (rng.normal(size=(n_arm, 40)) * scale) @ p_basis.T).astype(np.float32)
    tumour = (30.0 * rng.normal(size=(n_patients, 1)) * axis
              + (rng.normal(size=(n_patients, 40)) * scale) @ t_basis.T).astype(np.float32)
    cfg = TransferConfig(ks=(3,), offsets=(0, 1, 2), primary_offset=1, q=40, min_q=8, draws=30, bootstrap_draws=20)
    vt, ceiling = _tumour_reference(tumour, cfg)
    arm, _, _ = _arm_result(_rezscore(_tensor(perturbation)), vt, ceiling=ceiling, n_genes=n_genes,
                            cfg=cfg, draws=cfg.draws, seed=5, arm="all_perturbations")
    result = arm["k3"]
    assert set(result["overlap_by_offset"]) == {"0", "1", "2"}
    assert result["overlap_by_offset"]["0"] > 5 * result["overlap_by_offset"]["1"]
    assert result["pc1_share"] > .8, result["pc1_share"]
    # Contrast: structure spread across many shared directions is not a PC1 artefact.
    shared = np.linalg.qr(rng.normal(size=(n_genes, 40)))[0]
    spread_p = ((rng.normal(size=(n_arm, 40)) * scale) @ shared.T).astype(np.float32)
    spread_t = ((rng.normal(size=(n_patients, 40)) * scale) @ shared.T).astype(np.float32)
    vt2, ceiling2 = _tumour_reference(spread_t, cfg)
    spread = _arm_result(_rezscore(_tensor(spread_p)), vt2, ceiling=ceiling2, n_genes=n_genes,
                         cfg=cfg, draws=cfg.draws, seed=5, arm="all_perturbations")[0]["k3"]
    assert spread["pc1_share"] < .5, spread["pc1_share"]


# --- E0b algebraic rank ----------------------------------------------------

def test_algebraic_rank_recovers_known_rank():
    """A matrix of true rank 50 must report ~50, not min(n, p).

    The pre-fix cut was ``s > s[0] * 1e-10`` on a float32 spectrum whose eps is
    ~1.2e-7, so it never fired and every dictionary came back at full algebraic
    rank -- the exact overclaim E0b exists to correct."""
    rng = np.random.default_rng(2)
    low = (rng.normal(size=(200, 50)) @ rng.normal(size=(50, 300))).astype(np.float32)
    bundle = MatrixBundle(low, [f"G{i}" for i in range(300)], [f"r{i}" for i in range(200)], {})
    report = _full_dictionary_rank(bundle, DEVICE)
    assert report["max_possible_rank"] == 200
    assert 48 <= report["numerical_rank"] <= 51, report["numerical_rank"]
    assert report["numerical_rank"] < report["max_possible_rank"]
    assert report["effective_rank"] <= report["numerical_rank"] + 1e-6
    # A genuinely full-rank block must still report full rank (column centring
    # costs exactly one degree of freedom), so the cut is not simply too harsh.
    full = MatrixBundle(rng.normal(size=(120, 200)).astype(np.float32), [f"G{i}" for i in range(200)],
                        [f"r{i}" for i in range(120)], {})
    assert _full_dictionary_rank(full, DEVICE)["numerical_rank"] >= 119


# --- the statistic itself --------------------------------------------------

def test_overlap_statistic_is_bounded_and_self_overlap_is_one():
    generator = torch.Generator().manual_seed(19)
    basis, _ = torch.linalg.qr(torch.randn(200, 60, generator=generator), mode="reduced")
    other, _ = torch.linalg.qr(torch.randn(200, 60, generator=generator), mode="reduced")
    for k in (1, 5, 20, 50):
        assert abs(_overlap(basis, basis, k, 0) - 1.0) < 1e-5, k
        assert abs(_overlap(basis, basis, k, 1) - 1.0) < 1e-5, k
        for offset in (0, 1, 5):
            value = _overlap(basis, other, k, offset)
            assert 0.0 <= value <= 1.0 + 1e-9, (k, offset, value)
    # Two independent Haar subspaces sit at k/n_genes, which is why this floor
    # cannot separate real biology from generic co-expression.
    assert _overlap(basis, other, 20, 1) < 8 * 20 / 200


def test_matched_spectrum_invariant_is_measured_not_restated():
    """The reported spectrum error must be a real measurement. It used to be
    ``abs(x - x)``, so the E0.e gate could not fail for any rotation at all."""
    generator = torch.Generator().manual_seed(21)
    target, _ = torch.linalg.qr(torch.randn(128, 12, generator=generator), mode="reduced")
    singular = torch.arange(12, 0, -1, dtype=torch.float32)
    null, invariant = _matched_spectrum_rotated_null(singular, 128, target, k=5, draws=30, seed=3, offset=1)
    assert invariant["spectrum_checked_draws"] >= 1
    assert invariant["spectrum_max_relative_singular_error"] < ROTATION_TOLERANCE
    assert invariant["rotation_orthonormality_max_abs_error"] < ROTATION_TOLERANCE
    assert abs(invariant["null_retained_energy"] - invariant["source_retained_energy"]) < ROTATION_TOLERANCE * invariant["source_retained_energy"]
    assert null.std() > 0
    # A real measurement of a float32 QR is small but strictly positive.
    # `abs(x - x)` is exactly 0.0 for every input, which is what the old
    # `spectrum_energy_abs_error` reported and why its gate could not fail.
    assert invariant["rotation_orthonormality_max_abs_error"] > 0.0
    assert invariant["spectrum_max_relative_singular_error"] > 0.0
    assert "spectrum_energy_abs_error" not in invariant


# --- uncertainty on the observed statistic ---------------------------------

def test_bootstrap_ci_is_not_the_null_interval():
    """`effect_ci95` was the null's dispersion mirrored around the observation:
    algebraically the same test as E0.b, counted twice. The interval E0 now
    reports for the observed statistic must come from resampling P's rows."""
    responsive, control, tumour = _module_fixture(planted=True, seed=13)
    responsive_arm, _ = _run_arms(responsive, control, tumour, seed=31)
    for k in CFG.ks:
        r = responsive_arm[f"k{k}"]
        observed, null_median = r["pc1_removed_overlap"], r["null_median"]
        # `observed - null_effect_ci95` recovers the null's own quantiles, which
        # is all the old interval ever contained.
        null_quantiles = [observed - r["null_effect_ci95"][1], observed - r["null_effect_ci95"][0]]
        assert r["bootstrap_draws"] == CFG.bootstrap_draws
        assert r["bootstrap_std"] > 0 and r["bootstrap_ci95"][1] > r["bootstrap_ci95"][0]
        # Not the null interval, and not the null interval shifted onto the
        # observation either.
        assert abs(r["bootstrap_ci95"][0] - r["null_effect_ci95"][0]) > 1e-6
        assert abs(r["bootstrap_ci95"][1] - r["null_effect_ci95"][1]) > 1e-6
        assert abs(r["bootstrap_ci95"][0] - null_quantiles[0]) > 1e-6
        # The bootstrap is an interval for the OBSERVED statistic; the null is not.
        assert abs(r["bootstrap_median"] - observed) < abs(null_median - observed)
        assert r["bootstrap_ci95"][0] > null_median


def test_bootstrap_resamples_rows_and_is_seed_deterministic():
    responsive, _, tumour = _module_fixture(planted=True, seed=15)
    vt, _ = _tumour_reference(tumour)
    x = _rezscore(_tensor(responsive))
    first = _bootstrap_overlap(x, vt, cfg=CFG, seed=101)
    again = _bootstrap_overlap(x, vt, cfg=CFG, seed=101)
    moved = _bootstrap_overlap(x, vt, cfg=CFG, seed=202)
    for k in CFG.ks:
        assert first[k].size == CFG.bootstrap_draws
        np.testing.assert_allclose(first[k], again[k], rtol=0, atol=0)
        assert first[k].std() > 0
        assert not np.array_equal(first[k], moved[k])


# --- normalisation against the ceiling -------------------------------------

def test_normalised_alignment_is_reported_against_the_ceiling():
    """`_split_half` always computed a valid ceiling and nothing ever compared
    the observation to it, so "alignment = 0.4" stayed uninterpretable."""
    responsive, control, tumour = _module_fixture(planted=True, seed=17)
    responsive_arm, control_arm = _run_arms(responsive, control, tumour, seed=41)
    for k in CFG.ks:
        r, c = responsive_arm[f"k{k}"], control_arm[f"k{k}"]
        assert r["normalised_alignment_status"] == "available"
        expected = (r["pc1_removed_overlap"] - r["null_median"]) / (r["pc1_removed_ceiling"] - r["null_median"])
        assert abs(r["normalised_alignment"] - expected) < 1e-9
        assert r["normalised_alignment"] > c["normalised_alignment"]


def test_normalised_alignment_guards_a_degenerate_denominator():
    """A ceiling that collapses onto the null makes the ratio meaningless; it
    must be declared unavailable rather than divided into."""
    responsive, _, tumour = _module_fixture(planted=True, seed=19)
    vt, ceiling = _tumour_reference(tumour)
    collapsed = {k: 0.0 for k in CFG.ks}
    arm, _, _ = _arm_result(_rezscore(_tensor(responsive)), vt, ceiling=collapsed, n_genes=tumour.shape[1],
                            cfg=CFG, draws=CFG.draws, seed=3, arm="all_perturbations")
    for k in CFG.ks:
        assert not np.isfinite(arm[f"k{k}"]["normalised_alignment"])
        assert arm[f"k{k}"]["normalised_alignment_status"] == "unavailable_degenerate_ceiling_null_gap"
    assert ceiling[CFG.ks[0]] > 0  # the real ceiling is not degenerate


# --- machinery coverage -----------------------------------------------------
# Mutation testing showed the suite stopped at `_decision`: breaking observe(),
# _verdict's control-arm requirement, main's FAILED_ path, or z-scoring one arm
# differently from the other all passed 14/14. These three close that gap.

def test_observed_rows_never_fail_a_run():
    """A scientific outcome must not decide whether the run is valid. If observe()
    is reverted to a FAIL-able gate, a true negative gets quarantined into
    FAILED_*.json and reads as a crash."""
    import tempfile
    from morpheus.v2.calibra.gates import GateLedger
    with tempfile.TemporaryDirectory() as tmp:
        ledger = GateLedger(tmp, "T")
        ledger.add("G_valid", 1, "ok", True)
        ledger.observe("E0.g_responsive_exceeds_nonresponsive_k10", "false", "recorded finding")
        assert ledger.write() is True, "an OBSERVED row must not fail the run"
        assert [r["status"] for r in ledger.rows] == ["PASS", "OBSERVED"]
        ledger2 = GateLedger(tmp, "T")
        ledger2.add("G_valid", 1, "ok", False)
        assert ledger2.write() is False, "a FAIL row must still fail the run"


def test_verdict_never_supports_without_a_control_arm():
    """The false-positive guard. 'Supported' requires the biological control arm to
    have actually been compared -- that is the entire point of the redesign."""
    from morpheus.v2.calibra.e0_basis_transfer import _verdict
    ks = (10,)
    unavailable = {"k10": {"decision_n_matched": {"decision_status": "unavailable_insufficient_rows_50_lt_151",
                                                  "responsive_exceeds_nonresponsive": False}}}
    out = _verdict({"contexts": {"RPE1_signed_log1p": unavailable}}, ks)
    assert out["verdict"] == "indeterminate"
    assert out["transfer_supported"] is False

    good = {"k10": {"decision_n_matched": {"decision_status": "available",
                                           "responsive_exceeds_nonresponsive": True,
                                           "responsive_beats_haar_floor_only": True}}}
    mixed = _verdict({"contexts": {"K562_signed_log1p": good, "RPE1_signed_log1p": unavailable}}, ks)
    # B2: an undecidable context must not veto a decidable one, but must stay visible.
    assert mixed["verdict"] == "supported", mixed
    assert mixed["contexts_undecidable"] == ["RPE1_signed_log1p"]
    assert mixed["n_contexts_decided"] == 1


def test_both_arms_are_processed_identically():
    """Arm asymmetry silently invalidates the whole comparison: any difference in
    z-scoring, offset, k, q or null between the arms makes 'responsive > control'
    meaningless. Identical rows under different arm names must give identical numbers."""
    import torch
    from morpheus.v2.calibra.e0_basis_transfer import _arm_result, TransferConfig, _rezscore, _right_svd
    rng = np.random.default_rng(5)
    cfg = TransferConfig(q=60, min_q=8, ks=(10,), bootstrap_draws=20)
    rows = rng.normal(size=(120, 200)).astype(np.float32)
    target = rng.normal(size=(150, 200)).astype(np.float32)
    x = torch.as_tensor(rows); vt, _ = _right_svd(_rezscore(torch.as_tensor(target)), cfg.q, 1, min_q=cfg.min_q)
    ceiling = {10: 0.5}
    a, _, _ = _arm_result(x, vt, ceiling=ceiling, n_genes=200, cfg=cfg, draws=100, seed=7, arm="responsive_matched")
    b, _, _ = _arm_result(x, vt, ceiling=ceiling, n_genes=200, cfg=cfg, draws=100, seed=7, arm="nonresponsive")
    assert a["k10"]["pc1_removed_overlap"] == b["k10"]["pc1_removed_overlap"]
    assert a["k10"]["bootstrap_ci95"] == b["k10"]["bootstrap_ci95"]


def test_n_matching_is_symmetric():
    """Whichever arm is larger must be downsampled. Subsampling only the responsive
    arm silently sets n_matched=False whenever the control arm is bigger -- which
    happens as soon as nonresponsive_p is relaxed (K562 at p>0.05: 4251 control vs
    3132 responsive) -- auto-failing a context for bookkeeping, not science."""
    from morpheus.v2.calibra.e0_basis_transfer import _energy_arms, TransferConfig
    cfg = TransferConfig(responsive_p=0.01, nonresponsive_p=0.05)
    # control arm deliberately LARGER than the responsive arm
    energy = np.concatenate([np.full(30, 0.001), np.full(120, 0.5)])
    rows, meta = _energy_arms(energy, cfg=cfg, seed=0)
    assert len(rows["responsive_matched"]) == len(rows["nonresponsive"]) == 30, meta
    assert meta["arms_are_n_matched"] is True
    # and the usual direction still works
    energy2 = np.concatenate([np.full(200, 0.001), np.full(40, 0.9)])
    rows2, meta2 = _energy_arms(energy2, cfg=cfg, seed=0)
    assert len(rows2["responsive_matched"]) == len(rows2["nonresponsive"]) == 40, meta2
    # subsampled rows must be a genuine subset of the source arm
    assert set(rows2["responsive_matched"]).issubset(set(rows2["responsive_full"]))
