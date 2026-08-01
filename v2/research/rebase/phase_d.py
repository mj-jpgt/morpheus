"""Phase-D reproducible launch planning and paired-arm enforcement.

The runner validates its own invocation against the manifest created here, so
the H/PBS comparison cannot be made asymmetric by accidentally changing one
flag in a notebook or an overnight shell command.  This controller is designed
to be copied unchanged to Lambda and may be dry-run locally.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _require_axis_annotation_report(target: Path) -> None:
    """D2.3 is mandatory, so an unannotated PBS axis is not a completed D2."""
    with np.load(target, allow_pickle=False) as raw:  # imported lazily; controller otherwise stdlib-only
        manifest = json.loads(str(np.asarray(raw["manifest_json"]).item()))
    if not str(manifest.get("cohort_fit_population", "")).startswith("development"):
        raise RuntimeError("D2 blocked: final PBS target transform was not fit on development patients only")
    report = target.with_suffix(target.suffix + ".axis_annotations.csv")
    declared = manifest.get("axis_annotations", {})
    if declared.get("status") != "annotated" or not report.is_file() or _sha256(report) != declared.get("sha256"):
        raise RuntimeError("D2.3 blocked: PBS axis annotation report is absent, unannotated, or provenance-mismatched")
    with report.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"axis", "proliferation_loading", "essentiality_loading", "annotation_status"}
    if not rows or not required.issubset(rows[0]) or any(row.get("annotation_status") != "annotated" for row in rows):
        raise RuntimeError("D2.3 blocked: PBS axis report lacks complete proliferation/essentiality annotations")
    try:
        if any(not (math.isfinite(float(row["proliferation_loading"])) and
                    math.isfinite(float(row["essentiality_loading"]))) for row in rows):
            raise ValueError
    except (TypeError, ValueError):
        raise RuntimeError("D2.3 blocked: PBS axis report contains unavailable proliferation/essentiality values")


def _validate_pbs_role(target: Path, components: int, role: str) -> None:
    with np.load(target, allow_pickle=False) as raw:
        manifest = json.loads(str(np.asarray(raw["manifest_json"]).item()))
    if int(manifest.get("n_components", -1)) != components:
        raise RuntimeError("PBS artifact component count differs from requested D2 role")
    if role == "primary" and components != 128:
        raise RuntimeError("D2 primary is pre-registered at 128 PBS components; 64/256 are sensitivity only")
    if role == "sensitivity" and components not in {64, 256}:
        raise RuntimeError("D2 sensitivity role is reserved for the predeclared 64/256 component runs")


def d2_pair_manifest(args: argparse.Namespace) -> dict:
    target = Path(args.pbs_targets).resolve()
    if not target.is_file():
        raise FileNotFoundError(f"PBS target artifact does not exist: {target}")
    # These are precisely the knobs that can change the learned representation.
    # The runner checks them independently for both arms.  It is intentional
    # that output directory and target source are absent: they are the sole
    # permitted differences.
    digest = lambda path: _sha256(Path(path).resolve())
    common = {
        "data_config": str(Path(args.data_config).resolve()), "data_config_sha256": digest(args.data_config),
        "split_file": str(Path(args.split_file).resolve()), "split_file_sha256": digest(args.split_file),
        "epochs": args.epochs, "token_budget": args.token_budget, "hidden_dim": args.hidden_dim,
        "layers": args.layers, "heads": args.heads, "learning_rate": args.learning_rate,
        "weight_decay": args.weight_decay, "device": args.device, "mlp_clip_teacher": "",
        "mlp_clip_anchor": "", "teacher_warmup_epochs": 0,
        "gradient_diagnostics_every": args.gradient_diagnostics_every,
        "objective_profile": "programme_only", "decorrelation_weight": args.decorrelation_weight,
        "loss_warmup_epochs": args.loss_warmup_epochs, "programme_warmup_weight": args.programme_warmup_weight,
        "programme_weight": args.programme_weight, "programme_neighbourhood_weight": args.programme_neighbourhood_weight,
        "programme_supcon_weight": args.programme_supcon_weight, "separation_weight": args.separation_weight,
        "variance_weight": args.variance_weight, "programme_head_dim": args.programme_head_dim,
        "pretrain_epochs": 0, "pretrain_checkpoint": "",
        "pretrain_learning_rate": args.pretrain_learning_rate, "pretrain_mask_fraction": args.pretrain_mask_fraction,
        "pretrain_view_keep_fraction": args.pretrain_view_keep_fraction, "pretrain_target_dim": args.pretrain_target_dim,
        "snv_features": "", "cnv_features": "", "plip_teacher": "", "include_clinical": False,
        "resume": "", "fit_development": True, "fixed_final_epoch": True,
        "expected_development_cancers": args.expected_development_cancers,
        "expected_heldout_cancers": args.expected_heldout_cancers, "fit_programme_legibility": True,
        "d2_analysis_role": args.analysis_role, "d2_pbs_components": args.pbs_components,
    }
    return {
        "schema_version": 1, "experiment": "D2_H_vs_I",
        "target_only_difference": True, "fit_programme_legibility": True,
        "common_args": common, "common_config_sha256": hashlib.sha256(json.dumps(common, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest(),
        "targets": {"H": "", "I": str(target)},
        "pbs_target_sha256": _sha256(target), "seeds": args.seeds,
        "pbs_components": args.pbs_components, "analysis_role": args.analysis_role,
        "data_config": str(Path(args.data_config).resolve()), "split_file": str(Path(args.split_file).resolve()),
        "predeclared_primary": "legibility_weighted_grouped_cv",
        "notes": "H and I use the identical model, schedule, seed, split, and target-axis weighting procedure; only supervision coordinates differ.",
    }


def _runner_command(args: argparse.Namespace, manifest: Path, arm: str, seed: int) -> list[str]:
    output = Path(args.run_root).resolve() / f"d2_{arm.lower()}_seed{seed}"
    command = [sys.executable, "-m", "morpheus.v2.runner", "--data-config", str(Path(args.data_config).resolve()),
               "--split-file", str(Path(args.split_file).resolve()), "--output-dir", str(output),
               "--objective-profile", "programme_only", "--epochs", str(args.epochs),
               "--token-budget", str(args.token_budget), "--hidden-dim", str(args.hidden_dim),
               "--layers", str(args.layers), "--heads", str(args.heads), "--learning-rate", str(args.learning_rate),
               "--weight-decay", str(args.weight_decay), "--decorrelation-weight", str(args.decorrelation_weight),
               "--loss-warmup-epochs", str(args.loss_warmup_epochs), "--seed", str(seed), "--device", args.device,
               "--fit-development", "--fixed-final-epoch", "--fit-programme-legibility",
               "--programme-warmup-weight", str(args.programme_warmup_weight), "--programme-weight", str(args.programme_weight),
               "--programme-neighbourhood-weight", str(args.programme_neighbourhood_weight), "--programme-supcon-weight", str(args.programme_supcon_weight),
               "--separation-weight", str(args.separation_weight), "--variance-weight", str(args.variance_weight),
               "--programme-head-dim", str(args.programme_head_dim),
               "--gradient-diagnostics-every", str(args.gradient_diagnostics_every),
               "--pretrain-learning-rate", str(args.pretrain_learning_rate), "--pretrain-mask-fraction", str(args.pretrain_mask_fraction),
               "--pretrain-view-keep-fraction", str(args.pretrain_view_keep_fraction), "--pretrain-target-dim", str(args.pretrain_target_dim),
               "--expected-development-cancers", str(args.expected_development_cancers), "--expected-heldout-cancers", str(args.expected_heldout_cancers),
               "--d2-pair-manifest", str(manifest), "--d2-arm", arm, "--d2-analysis-role", args.analysis_role,
               "--d2-pbs-components", str(args.pbs_components)]
    if arm == "I":
        command += ["--programme-targets", str(Path(args.pbs_targets).resolve())]
    return command


def _append_gate_log(path: Path, experiment: str, gate: str, value: str, threshold: str, status: str) -> None:
    if not path.exists():
        path.write_text("| experiment | gate | value | threshold | status |\n|---|---|---|---|---|\n", encoding="utf-8")
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"| {experiment} | {gate} | {value} | {threshold} | {status} |\n")


def run_d2(args: argparse.Namespace) -> None:
    root = Path(args.run_root).resolve(); root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "D2_PAIR_MANIFEST.json"
    manifest = d2_pair_manifest(args); _write_json(manifest_path, manifest)
    commands = [{"arm": arm, "seed": seed, "argv": _runner_command(args, manifest_path, arm, seed)}
                for seed in args.seeds for arm in ("H", "I")]
    _write_json(root / "D2_LAUNCH_PLAN.json", {"pair_manifest": str(manifest_path), "commands": commands})
    if not args.execute:
        print(json.dumps({"status": "planned", "pair_manifest": str(manifest_path), "n_commands": len(commands)}, indent=2))
        return
    repo = Path(__file__).resolve().parents[3]
    if subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain"], text=True).strip():
        raise RuntimeError("G0.2 failed: Phase-D execution requires a clean committed worktree")
    if (root / "SUCCESS.json").exists():
        raise RuntimeError(f"refusing an already-complete D2 run root {root}")
    _require_axis_annotation_report(Path(args.pbs_targets).resolve())
    _validate_pbs_role(Path(args.pbs_targets).resolve(), args.pbs_components, args.analysis_role)
    if not args.calibra_targets:
        raise RuntimeError("D2 execution requires --calibra-targets to run the same-run G4 controls and paired comparison")
    if not Path(args.calibra_targets).is_file():
        raise FileNotFoundError(args.calibra_targets)
    gate_log = Path(__file__).with_name("nature") / "GATE_LOG.md"
    for item in commands:
        run_dir = root / f"d2_{item['arm'].lower()}_seed{item['seed']}"
        if run_dir.exists():
            raise RuntimeError(f"refusing stale D2 output directory {run_dir}; use a new run root")
        result = subprocess.run(item["argv"], cwd=Path(__file__).resolve().parents[3])
        status = "PASS" if result.returncode == 0 else "FAIL"
        _append_gate_log(gate_log, f"D2_{item['arm']}_seed{item['seed']}", "runner_exit",
                         str(result.returncode), "0", status)
        if result.returncode:
            raise RuntimeError(f"D2 {item['arm']} seed {item['seed']} failed; do not compare incomplete arms")
        liveness_path = root / f"d2_{item['arm'].lower()}_seed{item['seed']}" / "liveness.json"
        liveness = json.loads(liveness_path.read_text(encoding="utf-8"))
        # The runner has already fail-closed on its profile-specific G2.6
        # criteria before it writes this file.  We still independently assert
        # that the recorded real-model overfit and step-one gradients exist.
        overfit = liveness.get("overfit_one_batch", {})
        passed = bool(overfit) and all(
            float(value) > 0.0 for value in liveness.get("gradient_norms_first", {}).values())
        _append_gate_log(gate_log, f"D2_{item['arm']}_seed{item['seed']}", "G2_liveness",
                         json.dumps(liveness.get("parameter_relative_delta", {}), sort_keys=True), 
                         "overfit+nonzero_grads", "PASS" if passed else "FAIL")
        if not passed:
            raise RuntimeError("D2 liveness gate failed; no scientific result may be logged")
        _write_json(run_dir / "TRAIN_SUCCESS.json", {"pair_manifest_sha256": _sha256(manifest_path), "arm": item["arm"],
                                                         "seed": item["seed"], "liveness_passed": True})
    exports: dict[str, list[str]] = {"H": [], "I": []}
    for seed in args.seeds:
        for arm in ("H", "I"):
            run_dir = root / f"d2_{arm.lower()}_seed{seed}"
            output = root / "artifacts" / f"d2_{arm.lower()}_seed{seed}.npz"
            command = [sys.executable, "-m", "morpheus.v2.export", "--data-config", str(Path(args.data_config).resolve()),
                       "--split-file", str(Path(args.split_file).resolve()), "--checkpoint", str(run_dir / "last.pt"),
                       "--output", str(output), "--token-budget", str(args.token_budget), "--hidden-dim", str(args.hidden_dim),
                       "--layers", str(args.layers), "--heads", str(args.heads), "--device", args.device]
            if arm == "I":
                command += ["--programme-targets", str(Path(args.pbs_targets).resolve())]
            result = subprocess.run(command, cwd=repo)
            _append_gate_log(gate_log, f"D2_{arm}_seed{seed}", "artifact_export", str(result.returncode), "0",
                             "PASS" if result.returncode == 0 else "FAIL")
            if result.returncode or not output.is_file():
                raise RuntimeError("D2 artifact export failed; no comparison may use the checkpoint directly")
            exports[arm].append(str(output))
    calibra_output = root / "calibra"
    calibra = [sys.executable, "-m", "morpheus.v2.calibra.run_calibra", "--artifacts", *exports["H"], *exports["I"],
               "--targets", str(Path(args.calibra_targets).resolve()), "--output", str(calibra_output), "--n-draws", "40",
               "--n-components", "16", "--n-permutations", "2000", "--n-jobs", str(args.calibra_jobs),
               "--require-rna-positive-control", "--require-channel-gates"]
    result = subprocess.run(calibra, cwd=repo)
    _append_gate_log(gate_log, "D2", "G4_CALIBRA_controls", str(result.returncode), "0", "PASS" if result.returncode == 0 else "FAIL")
    if result.returncode:
        raise RuntimeError("D2 CALIBRA G4 controls failed; no D2 result is valid")
    comparison = [sys.executable, "-m", "morpheus.v2.research.rebase.d2_compare", "--hallmark-artifacts", *exports["H"],
                  "--pbs-artifacts", *exports["I"], "--targets", str(Path(args.calibra_targets).resolve()),
                  "--output", str(root / "D2_PAIRED_BOOTSTRAP.json"), "--repeats", str(args.bootstrap_repeats)]
    result = subprocess.run(comparison, cwd=repo)
    _append_gate_log(gate_log, "D2", "paired_patient_and_cancer_bootstrap", str(result.returncode), "0", "PASS" if result.returncode == 0 else "FAIL")
    if result.returncode:
        raise RuntimeError("D2 paired bootstrap failed")
    _write_json(root / "SUCCESS.json", {"pair_manifest_sha256": _sha256(manifest_path),
                                          "calibra_targets_sha256": _sha256(Path(args.calibra_targets).resolve()),
                                          "pbs_target_sha256": _sha256(Path(args.pbs_targets).resolve()),
                                          "analysis_role": args.analysis_role, "pbs_components": args.pbs_components,
                                          "artifacts": exports, "g4_calibra_passed": True,
                                          "paired_bootstrap": str(root / "D2_PAIRED_BOOTSTRAP.json")})


def run_d3(args: argparse.Namespace) -> None:
    """Run the purity sensitivity with mandatory same-run G4 controls."""
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, "-m", "morpheus.v2.calibra.run_calibra", "--artifacts", *args.artifacts,
               "--targets", str(Path(args.targets).resolve()), "--output", str(output),
               "--purity-table", str(Path(args.purity_table).resolve()), "--purity-source", args.purity_source,
               "--purity-reference", args.purity_reference, "--purity-units", args.purity_units,
               "--n-draws", str(args.n_draws), "--n-components", str(args.n_components),
               "--n-permutations", str(args.n_permutations), "--n-jobs", str(args.n_jobs),
               "--require-rna-positive-control", "--require-channel-gates"]
    if args.purity_patient_column:
        command += ["--purity-patient-column", args.purity_patient_column]
    if args.purity_column:
        command += ["--purity-column", args.purity_column]
    _write_json(output / "D3_LAUNCH_PLAN.json", {"command": command,
                "invariant": "complete-case before/after uses identical patients and the full CALIBRA channel procedure"})
    if not args.execute:
        print(json.dumps({"status": "planned", "output": str(output)}, indent=2)); return
    repo = Path(__file__).resolve().parents[3]
    if subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain"], text=True).strip():
        raise RuntimeError("G0.2 failed: Phase-D execution requires a clean committed worktree")
    if any((output / name).exists() for name in ("task_rows.csv", "SUCCESS.json", "FAILED.json")):
        raise RuntimeError(f"refusing stale D3 output directory {output}; use a new output path")
    _write_json(output / "D3_INPUT_MANIFEST.json", {"artifacts": [{"path": str(Path(path).resolve()), "sha256": _sha256(Path(path).resolve())} for path in args.artifacts],
                                                       "targets": {"path": str(Path(args.targets).resolve()), "sha256": _sha256(Path(args.targets).resolve())},
                                                       "purity": {"path": str(Path(args.purity_table).resolve()), "sha256": _sha256(Path(args.purity_table).resolve())}})
    result = subprocess.run(command, cwd=repo)
    status = "PASS" if result.returncode == 0 else "FAIL"
    gate_log = Path(__file__).with_name("nature") / "GATE_LOG.md"
    _append_gate_log(gate_log, "D3", "G4_same_run_rna_positive_and_channel_gates", str(result.returncode), "0", status)
    if result.returncode:
        raise RuntimeError("D3 failed a required CALIBRA gate; do not interpret purity attenuation")
    _write_json(output / "SUCCESS.json", {"g4_controls_required": True, "purity_complete_case_paired": True})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="task", required=True)
    d2 = sub.add_parser("d2", help="write or execute a target-only paired H-vs-PBS D2 run")
    d2.add_argument("--data-config", default="morpheus/configs/v1.json")
    d2.add_argument("--split-file", required=True)
    d2.add_argument("--pbs-targets", required=True)
    d2.add_argument("--pbs-components", type=int, default=128, choices=(64, 128, 256))
    d2.add_argument("--analysis-role", default="primary", choices=("primary", "sensitivity"))
    d2.add_argument("--run-root", required=True)
    d2.add_argument("--seeds", default="42,43,44")
    d2.add_argument("--epochs", type=int, default=40); d2.add_argument("--token-budget", type=int, default=32768)
    d2.add_argument("--hidden-dim", type=int, default=512); d2.add_argument("--layers", type=int, default=4); d2.add_argument("--heads", type=int, default=8)
    d2.add_argument("--learning-rate", type=float, default=2e-4); d2.add_argument("--weight-decay", type=float, default=1e-2)
    d2.add_argument("--decorrelation-weight", type=float, default=.04); d2.add_argument("--loss-warmup-epochs", type=int, default=4)
    d2.add_argument("--programme-warmup-weight", type=float, default=.50); d2.add_argument("--programme-weight", type=float, default=1.0)
    d2.add_argument("--programme-neighbourhood-weight", type=float, default=.20); d2.add_argument("--programme-supcon-weight", type=float, default=.20)
    d2.add_argument("--programme-head-dim", type=int, default=256, help="fixed masked programme-head width shared by H and PBS")
    d2.add_argument("--separation-weight", type=float, default=.01); d2.add_argument("--variance-weight", type=float, default=.01)
    d2.add_argument("--gradient-diagnostics-every", type=int, default=25)
    d2.add_argument("--pretrain-learning-rate", type=float, default=2e-4); d2.add_argument("--pretrain-mask-fraction", type=float, default=.30)
    d2.add_argument("--pretrain-view-keep-fraction", type=float, default=.70); d2.add_argument("--pretrain-target-dim", type=int, default=128)
    d2.add_argument("--expected-development-cancers", type=int, default=11); d2.add_argument("--expected-heldout-cancers", type=int, default=22)
    d2.add_argument("--calibra-targets", default="", help="frozen RNA target artifact for mandatory G4 CALIBRA and paired bootstrap")
    d2.add_argument("--calibra-jobs", type=int, default=1); d2.add_argument("--bootstrap-repeats", type=int, default=2000)
    d2.add_argument("--device", default="cuda"); d2.add_argument("--execute", action="store_true")
    d3 = sub.add_parser("d3", help="run the complete-case purity sensitivity with mandatory G4 gates")
    d3.add_argument("--artifacts", nargs="+", required=True); d3.add_argument("--targets", required=True); d3.add_argument("--output", required=True)
    d3.add_argument("--purity-table", required=True); d3.add_argument("--purity-source", required=True,
                    choices=("published_consensus", "absolute", "expression_derived"))
    d3.add_argument("--purity-reference", required=True); d3.add_argument("--purity-units", default="fraction", choices=("fraction", "percent"))
    d3.add_argument("--purity-patient-column", default=""); d3.add_argument("--purity-column", default="")
    d3.add_argument("--n-draws", type=int, default=40); d3.add_argument("--n-components", type=int, default=16)
    d3.add_argument("--n-permutations", type=int, default=2000); d3.add_argument("--n-jobs", type=int, default=1)
    d3.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.task == "d2":
        args.seeds = [int(seed) for seed in str(args.seeds).split(",") if seed]
        if not args.seeds or len(set(args.seeds)) != len(args.seeds):
            raise ValueError("--seeds must be a nonempty unique comma-separated list")
        run_d2(args)
    else:
        run_d3(args)


if __name__ == "__main__":
    main()
