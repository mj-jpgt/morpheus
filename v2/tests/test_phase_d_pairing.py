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
