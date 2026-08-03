"""D1 Stage-1 audit + the two readouts this run's pipeline will not produce.

Written 2026-08-03 while D1 was still training (epoch 17/40), so the checks are
fixed before any number exists. Mirrors the Post-D2 Stage-1 checklist in
NOTEBOOK.md, adapted to D1.

This run's `run_d1` predates the three-readout change (commit a4c74da), so it
emits only the UNRESTRICTED bootstrap -- the one in which 50 of 90 targets are
`programme_only`'s own Hallmark supervision. This script runs the missing
stratified and random_control readouts and writes D1_READOUT_INDEX.json.

usage:
  python -m morpheus.v2.research.rebase.d1_audit <run_root> --targets <frozen_rna_targets.npz>
"""
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np

from morpheus.v2.calibra.spectral import CANONICAL, RANK_VARIANTS, effective_rank

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run_root")
parser.add_argument("--targets", required=True, help="frozen_rna_targets.npz")
parser.add_argument("--repo", default=str(Path(__file__).resolve().parents[3].parent),
                    help="cwd for the d2_compare subprocess (the dir containing the morpheus package)")
parser.add_argument("--seeds", default="42,43,44")
parser.add_argument("--repeats", type=int, default=2000)
parser.add_argument("--skip-readouts", action="store_true")
args = parser.parse_args()

ROOT = Path(args.run_root).resolve()
SKIP = args.skip_readouts
TARGETS = Path(args.targets).resolve()
REPO = Path(args.repo).resolve()
SEEDS = [int(x) for x in str(args.seeds).replace(" ", ",").split(",") if x]
STRATIFIED = ["heldout_pathway", "immune_tme", "tumour_state"]
OUT: dict[str, object] = {"run_root": str(ROOT)}
FAIL: list[str] = []


def check(name: str, passed: bool, detail: object) -> None:
    OUT[name] = {"pass": bool(passed), "detail": detail}
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)
    if not passed:
        FAIL.append(name)


# ---------------------------------------------------------------- A1 completeness
runs = {f"d1_{a}_seed{s}": ROOT / f"d1_{a}_seed{s}" for s in SEEDS for a in ("p", "f")}
present, gates = {}, {}
for name, d in runs.items():
    ok = (d / "TRAIN_SUCCESS.json").is_file()
    present[name] = ok
    if ok and (d / "liveness.json").is_file():
        lv = json.loads((d / "liveness.json").read_text())
        grads = lv.get("gradient_norms_first", {})
        gates[name] = {"overfit_present": bool(lv.get("overfit_one_batch")),
                       "all_grads_positive": bool(grads) and all(float(v) > 0 for v in grads.values())}
check("A1_all_six_runs_complete", all(present.values()) and all(
    g["overfit_present"] and g["all_grads_positive"] for g in gates.values()),
      {"train_success": present, "liveness": gates})

if FAIL:
    (ROOT / "D1_AUDIT.json").write_text(json.dumps({**OUT, "verdict": f"FAIL: {FAIL}"}, indent=2),
                                        encoding="utf-8")
    raise SystemExit("A1 failed: the run is incomplete, so no readout may be computed. "
                     "Do not compare incomplete arms.")

# ---------------------------------------------------------------- A2/A3 readouts
arts = ROOT / "artifacts"
P = [str(arts / f"d1_p_seed{s}.npz") for s in SEEDS]
F = [str(arts / f"d1_f_seed{s}.npz") for s in SEEDS]
readouts = {
    "stratified": (["--target-groups", *STRATIFIED], ROOT / "D1_PAIRED_BOOTSTRAP_STRATIFIED.json"),
    "random_control": (["--target-groups", "random_control"], ROOT / "D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json"),
    "unrestricted": ([], ROOT / "D1_PAIRED_BOOTSTRAP.json"),
}
for name, (groups, path) in readouts.items():
    if SKIP or path.is_file():
        continue
    cmd = [sys.executable, "-m", "morpheus.v2.research.rebase.d2_compare",
           "--hallmark-artifacts", *P, "--pbs-artifacts", *F, "--targets", str(TARGETS),
           "--output", str(path), "--repeats", str(args.repeats),
           "--label-a", "programme_only", "--label-b", "programme_free", "--experiment", "D1", *groups]
    print(f"running {name} readout ...", flush=True)
    if subprocess.run(cmd, cwd=REPO).returncode:
        raise SystemExit(f"{name} readout failed")


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def deltas(res: dict) -> list[dict]:
    key = "programme_free_minus_programme_only"
    return [{"seed": SEEDS[i],
             "point_programme_only": p["point_programme_only"],
             "point_programme_free": p["point_programme_free"],
             "patient": p[key]["patient"], "cancer": p[key]["cancer"]}
            for i, p in enumerate(res["pairs"])]


strat = load(readouts["stratified"][1]); unres = load(readouts["unrestricted"][1])
ctrl = load(readouts["random_control"][1])
ds, du = deltas(strat), deltas(unres)
sign = lambda xs: {int(np.sign(x)) for x in xs}
strat_signs = sign([d["patient"]["point_delta"] for d in ds])
unres_signs = sign([d["patient"]["point_delta"] for d in du])
check("A2_stratified_direction_matches_unrestricted",
      strat_signs == unres_signs,
      {"n_targets_stratified": strat["n_targets"], "n_targets_unrestricted": unres["n_targets"],
       "stratified": [{"seed": d["seed"], "point_delta": d["patient"]["point_delta"],
                       "ci95": [d["patient"]["ci95_low"], d["patient"]["ci95_high"]]} for d in ds],
       "unrestricted": [{"seed": d["seed"], "point_delta": d["patient"]["point_delta"]} for d in du],
       "note": "a direction that flips here means the unrestricted number was measuring "
               "whose training targets were on the exam"})

dc = deltas(ctrl)
# The null for a 16-component canonical correlation is NOT zero -- CCA is biased
# upward, more so as components grow -- so "at chance" cannot be tested against 0.
# The testable statement is that RANDOM targets must score BELOW REAL ones, for
# both arms and every seed. If noise scores as well as biology, the instrument is
# manufacturing signal and every number on the project is void, not just D1's.
per_seed_ctrl = []
for i, d in enumerate(dc):
    for arm in ("programme_only", "programme_free"):
        per_seed_ctrl.append({"seed": d["seed"], "arm": arm,
                              "random_control": d[f"point_{arm}"],
                              "stratified": ds[i][f"point_{arm}"],
                              "below_real": bool(d[f"point_{arm}"] < ds[i][f"point_{arm}"])})
check("A3_random_control_scores_below_real_targets",
      all(r["below_real"] for r in per_seed_ctrl),
      {"n_targets": ctrl["n_targets"], "per_seed_arm": per_seed_ctrl,
       "delta_ci95": [[d["patient"]["ci95_low"], d["patient"]["ci95_high"]] for d in dc],
       "note": "REQUIRES HUMAN READ even when PASS: 'below real' is necessary, not sufficient. "
               "Judge the MARGIN. A control close to the real score is an alarm."})

# ---------------------------------------------------------------- A4 seed agreement
overlap = []
for i in range(len(ds)):
    for j in range(i + 1, len(ds)):
        a, b = ds[i]["patient"], ds[j]["patient"]
        overlap.append(bool(a["ci95_low"] <= b["ci95_high"] and b["ci95_low"] <= a["ci95_high"]))
check("A4_seed_agreement", len(strat_signs) == 1 and all(overlap),
      {"same_sign": len(strat_signs) == 1, "pairwise_ci_overlap": overlap,
       "per_seed": [{"seed": d["seed"], "point_delta": d["patient"]["point_delta"],
                     "ci95": [d["patient"]["ci95_low"], d["patient"]["ci95_high"]],
                     "p_improve": d["patient"]["p_improve"]} for d in ds]})

# ---------------------------------------------------------------- A5 effective rank
# This check previously carried its own inline definition -- the order-2 participation
# ratio (Sum sigma)^2 / Sum sigma^2, statistic R2 of paper/P2_RANK_DRAFT.md 3.1 -- which
# is NOT Roy & Vetterli's effective rank and is systematically LOWER than it. It now
# calls the one canonical implementation. R2 is reported alongside so that this table
# can be read against the historical D1-A numbers, but the headline column is R1.
def _rank_of_artifact(path: str) -> dict[str, float]:
    z = np.load(path, allow_pickle=True)
    x = np.asarray(z["wsi_biology"], float)[np.asarray(z["split"]).astype(str) == "test"]
    return {name: effective_rank(x, variant=RANK_VARIANTS[name]) for name in ("R1", "R2", "R3")}


ranks = {f"seed{s}": {"programme_only": _rank_of_artifact(P[i]),
                      "programme_free": _rank_of_artifact(F[i])}
         for i, s in enumerate(SEEDS)}
check("A5_effective_rank_reported", True,
      {"effective_rank": ranks,
       "canonical": "R1 = " + CANONICAL.label + " (Roy & Vetterli 2007 Definition 1, column-centred). "
                    "R2/R3 are order-2 Hill numbers of the same spectrum, retained only for "
                    "comparability with pre-2026-08-05 numbers; they are always <= R1.",
       "note": "REPORTED, NOT INTERPRETED (blocker 5). WSI states are ~0.80 collinear at init, so a "
               "narrower rank may reflect resistance to an already-collapsed view rather than "
               "dictionary content."})

# ---------------------------------------------------------------- A6 provenance
man = ROOT / "D1_PAIR_MANIFEST.json"
manifest = load(man) if man.is_file() else {}
check("A6_pair_manifest_records_one_difference",
      bool(manifest.get("objective_only_difference")) and manifest.get("experiment") == "D1_supervision_ablation",
      {"common_config_sha256": manifest.get("common_config_sha256"),
       "seeds": manifest.get("seeds"), "arms": manifest.get("arms"),
       "preregistered_prediction": manifest.get("preregistered_prediction"),
       "note": "D1 matching holds BY CONSTRUCTION (run_d1 asserts the argv differ only in profile); "
               "the manifest is a record, not a runner-side enforcement like D2's."})

# ---------------------------------------------------------------- index + verdict
(ROOT / "D1_READOUT_INDEX.json").write_text(json.dumps({
    "headline": "D1_PAIRED_BOOTSTRAP_STRATIFIED.json",
    "headline_rationale": f"the {strat['n_targets']} targets neither arm trained on: {', '.join(STRATIFIED)}",
    "negative_control": "D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json",
    "negative_control_requirement": "both arms at chance; a channel here voids every number on the project",
    "secondary_do_not_headline": "D1_PAIRED_BOOTSTRAP.json",
    "secondary_caveat": f"scores all {unres['n_targets']} non-control targets, of which 50 are "
                        "hallmark_in_training -- programme_only's own supervision.",
    "preregistered": "NOTEBOOK_ENTRIES/d1_readout_preregistration_20260803T1700Z.md",
    "audit": "D1_AUDIT.json",
}, indent=2), encoding="utf-8")
OUT["verdict"] = "PASS" if not FAIL else f"FAIL: {FAIL}"
(ROOT / "D1_AUDIT.json").write_text(json.dumps(OUT, indent=2), encoding="utf-8")
print(f"\nverdict: {OUT['verdict']}")
print("A3 and A5 require a human read even when marked PASS.")
