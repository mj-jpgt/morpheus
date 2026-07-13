import json
from pathlib import Path

import pandas as pd

from morpheus.src.data.create_splits import create_splits


def test_create_splits_has_no_patient_overlap(tmp_path: Path):
    processed = tmp_path / "processed"
    processed.mkdir()
    master = pd.DataFrame(
        {
            "patient_id": [f"TCGA-AA-{i:04d}" for i in range(30)],
            "cancer_type": ["A"] * 15 + ["B"] * 15,
            "has_clinical": True,
            "has_rna": True,
            "has_wsi": True,
        }
    )
    master.to_parquet(processed / "master_patient_table.parquet", index=False)
    cfg = {
        "project_root": str(tmp_path),
        "seed": 1,
        "paths": {"processed_dir": "processed", "outputs_dir": "outputs"},
        "splits": {
            "val_fraction": 0.2,
            "test_fraction": 0.2,
            "heldout_cancers": ["A"],
            "domain_adaptation_cancers": ["A"],
            "adaptation_fractions": [0.5],
        },
    }
    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
    outputs = create_splits(cfg_path)
    split = json.load(open(outputs["stratified"], encoding="utf-8"))["patient_ids"]
    train, val, test = map(set, (split["train"], split["val"], split["test"]))
    assert not (train & val)
    assert not (train & test)
    assert not (val & test)
