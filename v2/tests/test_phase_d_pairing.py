import argparse
import hashlib
import json

import pytest

from morpheus.v2.runner import _validate_d2_pair


def _args(manifest, arm="H", targets=""):
    return argparse.Namespace(
        d2_pair_manifest=str(manifest), d2_arm=arm, programme_targets=targets,
        seed=42, d2_analysis_role="primary", d2_pbs_components=128,
        fit_programme_legibility=True, objective_profile="programme_only",
        fixed_final_epoch=True, epochs=40, token_budget=32768, hidden_dim=512,
        layers=4, heads=8, learning_rate=2e-4, weight_decay=1e-2,
        decorrelation_weight=0.04, loss_warmup_epochs=4, mlp_clip_anchor="", mlp_clip_teacher="",
        data_config=str(manifest.parent / "config.json"), split_file=str(manifest.parent / "split.json"), device="cuda",
        teacher_warmup_epochs=0, gradient_diagnostics_every=25, programme_warmup_weight=.5, programme_weight=1.,
        programme_neighbourhood_weight=.2, programme_supcon_weight=.2, separation_weight=.01, variance_weight=.01,
        programme_head_dim=256,
        pretrain_epochs=0, pretrain_checkpoint="", pretrain_learning_rate=2e-4, pretrain_mask_fraction=.3,
        pretrain_view_keep_fraction=.7, pretrain_target_dim=128, snv_features="", cnv_features="", plip_teacher="",
        include_clinical=False, resume="", fit_development=True, expected_development_cancers=11, expected_heldout_cancers=22,
    )


def test_d2_pair_manifest_rejects_arm_asymmetry(tmp_path):
    target = tmp_path / "pbs.npz"; target.touch()
    config, split = tmp_path / "config.json", tmp_path / "split.json"
    config.write_text("{}"); split.write_text("{}")
    common = {"data_config": str(config.resolve()), "data_config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
              "split_file": str(split.resolve()), "split_file_sha256": hashlib.sha256(split.read_bytes()).hexdigest(),
              "objective_profile": "programme_only", "fixed_final_epoch": True, "epochs": 40,
              "token_budget": 32768, "hidden_dim": 512, "layers": 4, "heads": 8,
              "learning_rate": 2e-4, "weight_decay": 1e-2, "device": "cuda", "decorrelation_weight": .04,
              "loss_warmup_epochs": 4, "mlp_clip_anchor": "", "mlp_clip_teacher": "", "teacher_warmup_epochs": 0,
              "gradient_diagnostics_every": 25, "programme_warmup_weight": .5, "programme_weight": 1.,
              "programme_neighbourhood_weight": .2, "programme_supcon_weight": .2, "separation_weight": .01,
              "variance_weight": .01, "pretrain_epochs": 0, "pretrain_checkpoint": "", "pretrain_learning_rate": 2e-4,
              "programme_head_dim": 256,
              "pretrain_mask_fraction": .3, "pretrain_view_keep_fraction": .7, "pretrain_target_dim": 128,
              "snv_features": "", "cnv_features": "", "plip_teacher": "", "include_clinical": False,
              "resume": "", "fit_development": True, "expected_development_cancers": 11,
              "expected_heldout_cancers": 22, "fit_programme_legibility": True,
              "d2_analysis_role": "primary", "d2_pbs_components": 128}
    manifest = tmp_path / "pair.json"
    manifest.write_text(json.dumps({"schema_version": 1, "experiment": "D2_H_vs_I", "common_args": common,
                                     "fit_programme_legibility": True,
                                     "common_config_sha256": hashlib.sha256(json.dumps(common, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                                     "pbs_target_sha256": hashlib.sha256(target.read_bytes()).hexdigest(), "seeds": [42],
                                     "analysis_role": "primary", "pbs_components": 128,
                                     "targets": {"H": "", "I": str(target)}}))
    assert _validate_d2_pair(_args(manifest)) is not None
    bad = _args(manifest); bad.token_budget = 16384
    with pytest.raises(ValueError, match="asymmetry"):
        _validate_d2_pair(bad)
    bad = _args(manifest); bad.fit_programme_legibility = False
    with pytest.raises(ValueError, match="legibility"):
        _validate_d2_pair(bad)


def test_pbs_exclusion_must_be_requested_bounded_and_recorded():
    """The 14 paired patients with no PanCan RNA row are excluded, not silently dropped.
    Silently changing the cohort would invalidate comparisons against every previously
    published E0/CALIBRA number with no trace that the people differ."""
    import inspect
    from morpheus.v2 import build_pbs_targets as mod
    source = inspect.getsource(mod.build_pbs_targets)
    signature = inspect.signature(mod.build_pbs_targets)
    # opt-in, and off by default
    assert signature.parameters["allow_missing_rna"].default is False
    assert signature.parameters["max_missing_rna_fraction"].default == 0.01
    # refuses without the flag, and refuses beyond the ceiling even with it
    assert "allow_missing_rna=True to exclude them" in source
    assert "exceeds the" in source and "cohort change, not an exclusion" in source
    # every excluded patient is named in the manifest, plus count/fraction/ceiling
    for field in ("rna_missing_excluded_patients", "rna_missing_excluded_count",
                  "rna_missing_excluded_fraction", "rna_missing_exclusion_ceiling",
                  "paired_patient_count_before_exclusion", "cohort_deviation"):
        assert field in source, f"manifest must record {field}"
    # the retained view, not the full split, must reach the saved arrays and the fit rows
    assert "patient_ids=patient_ids, split=split_labels" in source
    assert "_digest_strings(patient_ids[fit_rows])" in source
