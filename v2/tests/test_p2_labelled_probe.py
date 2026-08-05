"""`p2_labelled_probe.py` imports and runs end to end on a synthetic cohort.

The script produces the labelled linear probe `paper/P2_RANK_DRAFT.md` §2.5 and §6.2 both mark
as never run. It is tested for the same reason the other five P2 scripts are
(`test_p2_analysis_scripts.py`): a vendored script with no test can rot silently against the
`calibra` functions it calls, and every statistic in this one is imported rather than written.

The assertions are structural invariants that hold for ANY input, plus three that are the point
of the module:

* the probe's estimator is a **linear classifier on frozen features**, and its null is produced
  by the **same estimator on permuted labels** -- so a probe at chance is legible as such;
* the five-repeat floor and the between-arm delta are computed for **every** statistic, so the
  probe is judged against its own floor exactly as §4.1 judges rank against rank's;
* wiring the probe through `_sklearn_oof` must give the identical number as calling the same
  estimator by hand -- the check that would catch a re-implementation drifting from the canonical
  fold loop.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from morpheus.v2.calibra.confound_certificate import _encode, _stratified_folds
from morpheus.v2.research.rebase.p2 import p2_labelled_probe

N_PATIENTS = 160
N_DIMS = 24
TARGET_GROUPS = ["heldout_pathway"] * 4 + ["immune_tme"] * 2 + ["tumour_state"] * 2
N_TARGETS = len(TARGET_GROUPS)
ARMS = ["H42", "I42", "H43", "I43", "H44", "I44", "P42", "F42", "P43", "F43", "P44", "F44"]
LABELS = ARMS + p2_labelled_probe.ENVELOPE


def _patient_ids(n: int) -> np.ndarray:
    return np.asarray([f"TCGA-{i % 6:02d}-{i:04d}" for i in range(n)])


def _cancers(n: int) -> np.ndarray:
    return np.asarray([f"C{i % 4}" for i in range(n)])


def _artifact(path, rng, *, decay: float, signal: float) -> None:
    """A frozen-artifact-shaped .npz: three co-trained views, a split, a cohort.

    `signal` controls how much of the cancer label is linearly recoverable, so the twelve arms
    give the probe a real spread rather than ties.
    """
    basis = np.linalg.qr(rng.normal(size=(N_DIMS, N_DIMS)))[0]
    latent = rng.normal(size=(N_PATIENTS, N_DIMS)) * decay ** np.arange(N_DIMS)
    class_index = _encode(_cancers(N_PATIENTS))[1]
    onehot = np.zeros((N_PATIENTS, N_DIMS))
    onehot[np.arange(N_PATIENTS), class_index] = 1.0
    wsi = (latent + signal * onehot) @ basis.T + 0.05 * rng.normal(size=(N_PATIENTS, N_DIMS))
    rna = (latent + signal * onehot) @ basis.T + 0.20 * rng.normal(size=(N_PATIENTS, N_DIMS))
    np.savez(path, patient_ids=_patient_ids(N_PATIENTS), cancers=_cancers(N_PATIENTS),
             split=np.asarray(["test"] * N_PATIENTS), wsi_biology=wsi, rna_biology=rna,
             full_biology=0.5 * (wsi + rna))


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    root = tmp_path_factory.mktemp("p2probe")
    rng = np.random.default_rng(20260805)
    np.savez(root / "targets.npz", patient_ids=_patient_ids(N_PATIENTS),
             target_names=np.asarray([f"T{i}" for i in range(N_TARGETS)]),
             target_groups=np.asarray(TARGET_GROUPS),
             scores=rng.normal(size=(N_PATIENTS, N_TARGETS)))
    paths = {}
    for index, label in enumerate(LABELS):
        path = root / f"{label}.npz"
        _artifact(path, rng, decay=0.80 + 0.02 * (index % 5), signal=0.6 + 0.1 * (index % 3))
        paths[label] = path
    # A label table shaped like E1's parquet, written as CSV so the test needs no pyarrow.
    # Balanced enough that `evaluate_known_covariate`'s >=30-per-class and
    # `within_cancer_auroc`'s >=5-minority-per-cancer guards are all satisfied.
    labels = pd.DataFrame({
        "patient_id": _patient_ids(N_PATIENTS), "cancer": _cancers(N_PATIENTS),
        "mut_TP53": ((np.arange(N_PATIENTS) % 3) == 0).astype(float),
    })
    labels.to_csv(root / "labels.csv", index=False)
    return root, paths


@pytest.fixture(scope="module")
def result(cohort):
    root, paths = cohort
    return p2_labelled_probe.main(
        ["--labels", str(root / "labels.csv"), "--targets", str(root / "targets.npz"),
         "--output", str(root / "P2_LABELLED_PROBE.json"), "--views", "wsi_biology",
         "--n-permutations-logistic", "3", "--n-permutations-lda", "5",
         "--n-boot", "40", "--n-permutations-auroc", "40",
         "--artifacts"] + [f"{k}={v}" for k, v in paths.items()])


def test_every_artifact_carries_both_probes_the_rank_and_the_channel(result):
    assert set(result) == set(LABELS) | {"_config", "_summary"}
    for label in LABELS:
        block = result[label]["views"]["wsi_biology"]
        assert np.isfinite(block["effective_rank_residualised"]) > 0
        assert 0.0 <= block["channel_untrained40"] <= 1.0
        raw = block["probe_A_cancer_type_raw"]
        for key in ("logistic_balanced_accuracy", "lda_balanced_accuracy"):
            assert 0.0 <= raw[key] <= 1.0, (label, key)
        assert raw["chance_rate"] == pytest.approx(1.0 / raw["n_classes"])
        # the null must be produced by the same estimator, and must sit near chance
        assert raw["lda_null"]["n_permutations"] == 5
        assert raw["logistic_null"]["n_permutations"] == 3
        tp53 = block["probe_B"]["mut_TP53"]
        assert tp53["status"] == "scored"
        assert tp53["adjustment"] == "cancer+pooled_tss_cross_fitted"
        assert np.isfinite(tp53["within_cancer_auroc"])
        assert np.isfinite(tp53["null_p95"])


def test_probe_beats_its_own_permutation_null_when_the_label_is_in_the_features(result):
    """The positive control. A probe that cannot clear its own null cannot order anything.

    The synthetic cohort puts the cancer label into the features on purpose, so if the observed
    value does not exceed the permutation null the wiring is wrong -- not the representation.
    """
    for label in LABELS:
        raw = result[label]["views"]["wsi_biology"]["probe_A_cancer_type_raw"]
        assert raw["lda_balanced_accuracy"] > raw["lda_null"]["null_max"], label
        assert raw["lda_null"]["permutation_p"] < 0.5, label


def test_adjusted_probe_is_the_must_fail_control_and_is_reported_not_asserted(result):
    """Cancer type after cancer+TSS residualisation must collapse towards chance.

    Asserted loosely here because the synthetic design is not the real confound design; what
    the test pins is that the control is COMPUTED for every artifact and every view, since a
    must-fail control that is silently skipped is worse than one that fails.
    """
    for label in LABELS:
        adjusted = result[label]["views"]["wsi_biology"]["probe_A_cancer_type_adjusted"]
        raw = result[label]["views"]["wsi_biology"]["probe_A_cancer_type_raw"]
        assert np.isfinite(adjusted["lda_balanced_accuracy"])
        assert adjusted["lda_balanced_accuracy"] < raw["lda_balanced_accuracy"], label


def test_summary_judges_every_statistic_against_its_own_five_repeat_floor(result):
    summary = result["_summary"]
    floors = summary["floors"]["wsi_biology"]
    assert set(floors) == set(p2_labelled_probe.STATISTICS)
    for name, floor in floors.items():
        assert floor["n"] == len(p2_labelled_probe.ENVELOPE), name
        assert floor["spread"] >= 0.0 and floor["ratio"] >= 1.0 - 1e-12, name
    rows = summary["pairs"]["wsi_biology"]
    assert len(rows) == len(p2_labelled_probe.PAIRS)
    for row in rows:
        for name in p2_labelled_probe.STATISTICS:
            cell = row["statistics"][name]
            assert cell["winner"] in (row["arm_a"], row["arm_b"], "tie")
            assert cell["delta"] == pytest.approx(cell["value_a"] - cell["value_b"])
            # the floor a difference is judged against must be the SAME statistic's own floor
            assert cell["floor_spread"] == pytest.approx(floors[name]["spread"])
            assert cell["resolvable_by_spread"] == (abs(cell["delta"]) > floors[name]["spread"])


def test_agreement_is_scored_against_the_channel_and_conflicts_name_rank(result):
    summary = result["_summary"]
    for name, record in summary["agreement"]["wsi_biology"].items():
        assert name != "channel_untrained40"
        assert record["agrees_with_channel"] + record["disagrees_with_channel"] == record["n_pairs"]
        assert record["n_pairs"] == len(p2_labelled_probe.PAIRS)
    for conflict in summary["channel_probe_conflicts"]:
        assert conflict["probe_winner"] != conflict["channel_winner"]
        assert conflict["rank_sides_with"] in ("probe", "channel")


def test_the_logistic_wrapper_is_the_canonical_fold_loop_not_a_reimplementation(cohort):
    """`logistic_balanced_accuracy_oof` must be `_sklearn_oof` plus an estimator, nothing else.

    If someone later inlines the fold loop "for clarity", this catches it: the wrapper has to
    reproduce, exactly, the number obtained by handing the same estimator factory to the
    canonical helper with the same folds.
    """
    from sklearn.linear_model import LogisticRegression

    from morpheus.v2.calibra.nonlinear_confound_probe import _sklearn_oof

    rng = np.random.default_rng(7)
    class_index = _encode(_cancers(N_PATIENTS))[1]
    onehot = np.zeros((N_PATIENTS, N_DIMS))
    onehot[np.arange(N_PATIENTS), class_index] = 1.0
    features = onehot + rng.normal(size=(N_PATIENTS, N_DIMS))
    folds = _stratified_folds(class_index, 5, 42)

    def _make():
        return LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced", solver="lbfgs",
                                  random_state=42, n_jobs=1)

    expected = _sklearn_oof(_make, features, class_index, 4, n_splits=5, seed=42, folds=folds,
                            standardise=True)
    assert p2_labelled_probe.logistic_balanced_accuracy_oof(
        features, class_index, 4, folds=folds, seed=42) == pytest.approx(expected)


def test_output_json_round_trips(cohort, result):
    root, _ = cohort
    on_disk = json.loads((root / "P2_LABELLED_PROBE.json").read_text(encoding="utf-8"))
    assert set(on_disk) == set(result)
    assert on_disk["_summary"]["floors"]["wsi_biology"]["probeA_lda"]["n"] == 5
