# P2 analysis scripts — vendored 2026-08-05

The five scripts in this directory produce `paper/P2_RANK_DRAFT.md` §4.2, §4.4, §4.5, §4.6 and
§4.7. Until this commit they existed **only** at `~/e0_run/p2_*.py` on `ubuntu@150.136.45.194`.
§6.2 of the draft lists vendoring them as a pre-submission requirement, because the paper's own
standard is that every quoted number traces to a file in the repository.

| script | draft section | reads | writes |
|---|---|---|---|
| `p2_competing_metrics.py` | §4.2, §4.4(1), §4.6 | frozen `.npz` artifacts + `frozen_rna_targets.npz` | `P2_METRICS_D2.json`, `P2_METRICS_D1.json` |
| `p2_selection_rule.py` | §4.6, §4.4(1) | the JSON above | stdout |
| `p2_necessity_and_variance.py` | **§4.2**, §4.7 | the JSON above | stdout |
| `p2_robustness.py` | §4.5(c) | artifacts + targets | `P2_ROBUSTNESS.json` |
| `p2_rank_variants.py` | §4.5(a) | artifacts + targets | optional JSON |

`v2/tests/test_p2_analysis_scripts.py` imports all five and runs them on a synthetic cohort, so
they cannot rot silently against the `calibra` functions they call.

## Running them

CPU only. **Thread caps are not optional** — multithreaded BLAS is ~23× slower on this stack, and
uncapped joblib workers each spawn a full BLAS pool
(`NOTEBOOK_ENTRIES/operational_shared_box_rules_20260804T0730Z.md`).

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
D2=~/e0_run/d2_v3; D1=~/e0_run/d1_v2/artifacts; T=~/e0_run/data/frozen_rna_targets.npz

python3 -m morpheus.v2.research.rebase.p2.p2_competing_metrics \
  --targets $T --output P2_METRICS_D2.json --subsamples 40 --artifacts \
  H42=$D2/d2_v3_s42/artifacts/d2_h_seed42.npz I42=$D2/d2_v3_s42/artifacts/d2_i_seed42.npz ...

python3 -m morpheus.v2.research.rebase.p2.p2_selection_rule P2_METRICS_D2.json P2_METRICS_D1.json
python3 -m morpheus.v2.research.rebase.p2.p2_necessity_and_variance P2_METRICS_D2.json P2_METRICS_D1.json
```

Run them from a workspace **verified equal to HEAD**, per-file, before quoting anything. The
reason is in `NOTEBOOK_ENTRIES/WORKSPACE_DRIFT_AUDIT_ALL_20260803T2359Z.md`: the workspace these
scripts originally ran in (`~/ws`) carried a `spectral.py` predating the rank canonicalisation
entirely, and the numbers survived only because the consolidation happened not to move the
default.

## Changes made at vendoring, and why

Recorded in full in `NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260805T0130Z.md`.

1. `p2_selection_rule.py`, `p2_necessity_and_variance.py` — **unchanged** apart from a header.
2. `p2_competing_metrics.py` — `main()` takes `argv` so it is testable; `participation_ratio_rownorm`
   added, used by no published number.
3. `p2_robustness.py` — its output path was hardcoded to `/home/ubuntu/e0_run/P2_ROBUSTNESS.json`
   and is now an argument.
4. `p2_rank_variants.py` — **rewritten**. The original began
   `sys.path.insert(0, "/home/ubuntu/ws")` and carried its own inline `R1`/`R2`/`R3`. It now calls
   `effective_rank(..., variant=RANK_VARIANTS[...])`. Its inline `R2` and `R3` were **not** R2 and
   R3: they were the order-2 Hill number of the *eigenvalue* distribution, `(Σσ²)²/Σσ⁴`, where
   `d1_audit.py`'s R2 — and therefore `RANK_VARIANTS["R2"]` — is `(Σσ)²/Σσ²`. Both statistics are
   now reported side by side (`R2`/`R3` against `PR`/`PR_rownorm`) so the difference is measured
   rather than assumed.
