"""Focused E1 tests: provenance refusal and rank/information measurement wiring."""
from __future__ import annotations

import json

import numpy as np
import pytest

from morpheus.v2.calibra.e1_rank_information import (
    heldout_directional_cca,
    run_e1,
    stable_rank,
    validate_matched_artifacts,
)


def _artifact(path, ids, cancers, split, state, *, decorrelation, extra=None):
    manifest = {"seed": 42, "epochs": 8, "token_budget": 16384,
                "config": {"decorrelation_weight": decorrelation, "hidden_dim": 256}}
    if extra:
        manifest.update(extra)
    np.savez_compressed(path, patient_ids=np.asarray(ids), cancers=np.asarray(cancers), split=np.asarray(split),
                        trained_states=np.asarray(["wsi_biology"]), manifest_json=np.asarray(json.dumps(manifest)),
                        wsi_biology=np.asarray(state, dtype=np.float32))


def _fixture(tmp_path, n=80):
    rng = np.random.default_rng(8)
    ids = [f"TCGA-{i % 8:02d}-{i:04d}" for i in range(n)]
    cancers = np.asarray([f"C{i % 4}" for i in range(n)])
    split = np.asarray(["test"] * n)
    latent = rng.normal(size=(n, 5))
    y = latent @ rng.normal(size=(5, 12)) + 0.2 * rng.normal(size=(n, 12))
    low = latent[:, :2] @ rng.normal(size=(2, 16))
    high = latent @ rng.normal(size=(5, 16)) + 0.1 * rng.normal(size=(n, 16))
    before, after, targets = tmp_path / "before.npz", tmp_path / "after.npz", tmp_path / "targets.npz"
    _artifact(before, ids, cancers, split, low, decorrelation=0.0)
    _artifact(after, ids, cancers, split, high, decorrelation=0.04)
    np.savez_compressed(targets, patient_ids=np.asarray(ids), scores=y.astype(np.float32),
                        target_names=np.asarray([f"T{i}" for i in range(y.shape[1])]))
    return before, after, targets


def test_stable_rank_distinguishes_collapsed_and_broad_states():
    rng = np.random.default_rng(0)
    low = rng.normal(size=(100, 2)) @ rng.normal(size=(2, 20))
    high = rng.normal(size=(100, 20))
    assert stable_rank(low) < 3
    assert stable_rank(high) > 5


def test_heldout_directional_cca_is_directionwise_and_non_circular():
    rng = np.random.default_rng(1)
    z = rng.normal(size=(300, 4))
    x = z @ rng.normal(size=(4, 16)) + .1 * rng.normal(size=(300, 16))
    y = z @ rng.normal(size=(4, 11)) + .1 * rng.normal(size=(300, 11))
    values = heldout_directional_cca(x, y, n_components=8, seed=2)
    assert len(values) == 8
    assert np.isfinite(values).all()
    assert values[0] > 0.7


def test_matched_artifacts_refuse_nonintervention_manifest_difference(tmp_path):
    before, after, _ = _fixture(tmp_path)
    with np.load(after, allow_pickle=True) as raw:
        ids, cancers, split, state = raw["patient_ids"], raw["cancers"], raw["split"], raw["wsi_biology"]
    _artifact(after, ids, cancers, split, state, decorrelation=.04, extra={"epochs": 9})
    with pytest.raises(ValueError, match="non-intervention manifest"):
        validate_matched_artifacts(before, after)


def test_run_e1_emits_paired_rows_protocol_and_arrays(tmp_path):
    before, after, targets = _fixture(tmp_path)
    out = tmp_path / "out"
    rows = run_e1(before, after, targets, out, n_components=4, levels=(0.0, .1, .3),
                  n_draws=3, n_permutations=3, min_site_count=2, seed=4)
    assert {"before", "after", "paired"}.issubset(set(rows["arm"]))
    paired = rows[(rows.arm == "paired") & (rows.metric == "delta_effective_rank")]
    assert len(paired) == 1 and paired.value.iloc[0] > 0
    assert (out / "task_rows.csv").exists()
    assert (out / "e1_protocol.json").exists()
    assert (out / "e1_manifest.json").exists()
    assert (out / "e1_arrays.npz").exists()
    protocol = json.loads((out / "e1_protocol.json").read_text())
    assert protocol["floor_metric"] == "paired_spike_recovery_one_direction_cca"
