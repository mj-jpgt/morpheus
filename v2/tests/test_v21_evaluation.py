from __future__ import annotations

import numpy as np
import pandas as pd

from morpheus.v2.v21_evaluation import evaluate


def _artifact(path, ids, cancers, split, state) -> None:
    np.savez_compressed(
        path,
        patient_ids=np.asarray(ids),
        cancers=np.asarray(cancers),
        split=np.asarray(split),
        trained_states=np.asarray(["wsi_identity"]),
        artifact_version=np.asarray(4),
        manifest_json=np.asarray("{}"),
        wsi_identity=np.asarray(state, dtype=np.float32),
    )


def test_strict_evaluator_writes_paired_predictions(tmp_path) -> None:
    ids = [f"TCGA-AA-{index:04d}" for index in range(16)]
    cancers = ["A"] * 4 + ["B"] * 4 + ["C"] * 4 + ["D"] * 4
    split = ["train"] * 4 + ["val"] * 4 + ["test"] * 8
    rng = np.random.default_rng(7)
    feature = rng.normal(size=(16, 4))
    teacher = tmp_path / "teacher.npz"
    challenger = tmp_path / "challenger.npz"
    _artifact(teacher, ids, cancers, split, feature)
    _artifact(challenger, ids, cancers, split, feature + rng.normal(scale=.01, size=feature.shape))
    targets = pd.DataFrame({
        "patient_id": ids,
        "programme": feature[:, 0] + feature[:, 1],
        "cytolytic_activity": feature[:, 0],
        "hypoxia": feature[:, 1],
        "emt": feature[:, 2],
        "proliferation": feature[:, 3],
        "RANDOM_CONTROL__programme__00": feature[:, 2],
    })
    target_path = tmp_path / "targets.parquet"
    targets.to_parquet(target_path, index=False)
    teacher_map = tmp_path / "teachers.json"
    teacher_map.write_text('{"challenger": "teacher"}', encoding="utf-8")
    output = tmp_path / "out"
    evaluate([str(teacher), str(challenger)], str(target_path), str(output), teacher_map_path=str(teacher_map), repeats=100)
    rows = pd.read_csv(output / "task_rows.csv")
    prediction = pd.read_parquet(output / "molecular_predictions.parquet")
    assert set(prediction.method) == {"teacher", "challenger"}
    assert ((rows.task == "molecular_prompting") & (rows.metric == "pearson")).any()
    assert (rows.task == "paired_vs_teacher").any()
    assert (rows.task == "programme_neighbour_retrieval").any()
    assert (rows.task == "molecular_phenotype").any()
    assert (rows.task == "random_signature_control").any()
    validation_output = tmp_path / "validation"
    evaluate([str(teacher)], str(target_path), str(validation_output), evaluation_partition="val", repeats=100)
    assert pd.read_parquet(validation_output / "molecular_predictions.parquet").patient_id.isin(ids[4:8]).all()
