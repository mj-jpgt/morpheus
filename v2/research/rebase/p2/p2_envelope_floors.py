"""The retraining floor, measured for EVERY statistic, view and block the paper uses.

Draft §4.1 measures one floor: canonical R1, on the exported `wsi_biology` block,
raw and residualised, across five identical `programme_only` retrains at seed 42
(`~/e0_run/d1_envelope/rep{1..5}.npz`). §4.1a then applies that one floor to
comparisons on eight other statistics and on four other blocks, which it flags
with a ‡ but cannot repair.

**This module repairs what is repairable.** The five repeats are already exported
and hold all three co-trained views, so every statistic in T1 and every view in
§4.5(c) can be given a floor of its own at no additional training cost. It is a
re-derivation from artifacts that exist, not a re-run.

WHAT IS AND IS NOT RECOVERABLE FROM THESE FIVE FILES.

  RECOVERABLE, and measured here: every statistic (R1, R2, R3, PR, PR_rownorm,
  RankMe, stable rank, α-ReQ, LiDAR, hard numerical rank) on every exported view
  (`wsi_biology`, `rna_biology`, `full_biology`), raw and residualised, plus the
  top-CCA channel on the same five runs so the asymmetry §4.1 rests on can be
  re-read under each statistic.

  NOT RECOVERABLE, and named as absent rather than left ambiguous: the fixed
  held-out probe, the in-run training batch and the 16-patient gate batch. The
  repeats were **exported**, not probed: `d1_envelope` holds `rep{n}.npz`
  artifacts and per-run `train_metrics.jsonl`, and neither carries a probe
  forward pass, a training batch's activations, or a gate batch. A floor on any
  of those three blocks requires re-running the five repeats with the probe
  attached, on a GPU. `ABSENT_BLOCKS` below states this per block, and the audit
  reads it: a comparison on one of those blocks is **unjudgeable**, which is
  neither a pass nor a failure.

WHAT A FLOOR HERE IS, AND IS NOT. It is `max / min` over five runs that differ
only in GPU non-determinism. It is a FLOOR TWICE OVER, exactly as §4.1 says of
the one it replaces: `programme_only` is this project's stable arm, and same-seed
repeats exclude seed variation entirely. n = 5, no interval, one arm, one seed,
one stack. A table of these numbers is not a table of estimated distributions and
must not be read as one.

THE SHAPE MATTERS AS MUCH AS THE NUMBER. §4.1's floor is bimodal — four runs
agree to 2% and rep2 lands at a third of them. Whether that shape survives a
change of statistic is a claim about the paper's thesis, not a detail: a
statistic in which rep2 is *not* the outlier would mean the divergence between
those five runs is visible to some measures and not others. `_shape` below
records it under a rule fixed here rather than chosen per statistic:

  * `outlier`  — the run whose removal minimises the fold of the remaining four;
  * `rest_fold`     — that remaining four's own fold;
  * `separation`    — full fold / rest fold;
  * `concordant`    — `rest_fold <= 1.05`;
  * `bimodal`       — `concordant` AND `separation >= 2.0`.

NOTHING IS COMPUTED INLINE. R1/R2/R3 and the channel come from
`v2.calibra.spectral`; RankMe, the participation ratios, stable rank, α-ReQ and
LiDAR from `p2_competing_metrics` (the four published alternatives, which must
NOT be routed through `spectral.py` — see §4.6). The hard numerical rank is
`numpy.linalg.matrix_rank`, a library call to the same routine
`d1_geometry_probe.py` uses, and is not an effective rank. Four inline-formula
substitutions have been caught in this paper and the tell was identical each
time: an arithmetic expression where an import belonged.

usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 python3 p2_envelope_floors.py \
        --reps '/home/ubuntu/e0_run/d1_envelope/rep*.npz' \
        --targets /home/ubuntu/e0_run/data/frozen_rna_targets.npz \
        --output ~/ws_floor/out/P2_ENVELOPE_FLOORS.json
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import (RANK_VARIANTS, effective_rank,
                                          top_canonical_correlation)
from morpheus.v2.research.rebase.p2.p2_competing_metrics import (alpha_req, lidar,
                                                                 participation_ratio,
                                                                 participation_ratio_rownorm,
                                                                 rankme, stable_rank)

#: The 40 targets neither D1 arm trained on — identical to `d1_envelope_readout.py`
#: and to `d2_compare`, so the channel measured here is the channel §4.1 quotes.
UNTRAINED_40 = ("heldout_pathway", "immune_tme", "tumour_state")

#: The three co-trained views the export carries. §4.5(c) reads rank on all three
#: and has a floor for none of them.
VIEWS = ("wsi_biology", "rna_biology", "full_biology")

#: α-ReQ's eigenvalue index range and LiDAR's ridge, at the values T1 was scored
#: with (`P2_METRICS_D2.json` `_config`). A floor measured at other settings would
#: not be a floor for the numbers in T1.
ALPHA_INDEX_RANGE = (11, 50)
LIDAR_DELTA = 1e-4

#: Statistic label -> callable on one (n, d) matrix. R1/R2/R3 are the canonical
#: variants of the one implementation; PR/PR_rownorm are the eigenvalue-distribution
#: participation ratios `p2_rank_variants.py` reports beside them; the rest are the
#: published alternatives of §4.6, at their published definitions.
STATISTICS = {
    "R1": lambda x: effective_rank(x, variant=RANK_VARIANTS["R1"]),
    "R2": lambda x: effective_rank(x, variant=RANK_VARIANTS["R2"]),
    "R3": lambda x: effective_rank(x, variant=RANK_VARIANTS["R3"]),
    "PR": participation_ratio,
    "PR_rownorm": participation_ratio_rownorm,
    "RankMe": rankme,
    "stable_rank": stable_rank,
    "alpha_req": lambda x: alpha_req(x, *ALPHA_INDEX_RANGE)["alpha"],
    # T1 scores α-ReQ as |α − 1|, which is this project's operationalisation and
    # not the authors' — §4.6 says so. Both are given a floor because the
    # selection rule is applied to the second.
    "alpha_req_abs_dev": lambda x: abs(alpha_req(x, *ALPHA_INDEX_RANGE)["alpha"] - 1.0),
    # NOT an effective rank. `numpy.linalg.matrix_rank`, the same routine
    # `d1_geometry_probe.py` calls through torch; §4.9's "16/16" instance is this
    # statistic, on a different block.
    "hard_rank": lambda x: float(np.linalg.matrix_rank(np.asarray(x, dtype=np.float64))),
}

#: LiDAR needs a positive PAIR, so it is a per-(rep, block) quantity over two views
#: rather than a per-view one, exactly as `p2_competing_metrics.score_artifact`
#: computes it. Its pseudo-view name is recorded so the floor cannot be matched to
#: a single-view comparison by accident.
LIDAR_VIEW = "wsi_rna_paired"

#: Blocks on which a floor CANNOT be measured from these exports, and what each
#: would cost. This is a result, not a caveat: it is exactly the list of blocks
#: draft §5's rank numbers live on.
ABSENT_BLOCKS = {
    "fixed held-out probe": (
        "The repeats were EXPORTED, not probed. `~/e0_run/d1_envelope/` holds "
        "`rep{n}.npz` and per-run `train_metrics.jsonl` and neither carries a probe "
        "forward pass. A floor here needs the five repeats re-run with "
        "`d1_momentum_probe.py` attached — GPU, five runs. Every rank number in "
        "draft §5 is on this block."),
    "training batch, in-run": (
        "The in-run tripwire reads R3 on the live training batch; those activations "
        "are never saved (`F3_TRIPWIRE_STEP200_R3_n5.json`: 'states never saved, so "
        "[NOT RECOMPUTABLE]'). A floor here needs the five repeats re-run with the "
        "tripwire logged at a fixed step — GPU."),
    "16-patient gate batch": (
        "The liveness gate's batch is constructed inside the gate run and is not an "
        "export. A floor here needs five identical gate runs — GPU. §4.9's '16/16' "
        "instance is on this block, and it is a hard numerical rank pinned at the "
        "batch size in every arm."),
    "282 held-out patients, live checkpoint": (
        "Read off a live D1-A checkpoint that was never exported; §4.9 records it as "
        "[NOT RECOMPUTED] under R1. A floor here needs the checkpoint re-created."),
}


# --------------------------------------------------------------------------
def _shape(values: dict[str, float]) -> dict:
    """The bimodality descriptor, under the rule fixed in this module's docstring.

    Reported per statistic because whether the shape survives a change of
    statistic is on the paper's thesis. Not a distributional claim: n = 5.
    """
    names = sorted(values)
    numbers = [values[n] for n in names]
    finite = [v for v in numbers if np.isfinite(v)]
    if len(finite) != len(numbers) or min(finite) <= 0:
        return {"defined": False,
                "why": "a non-finite or non-positive value — a fold is not defined"}
    full = max(numbers) / min(numbers)
    best = None
    for i in range(len(numbers)):
        rest = numbers[:i] + numbers[i + 1:]
        fold = max(rest) / min(rest)
        if best is None or fold < best[0]:
            best = (fold, names[i])
    rest_fold, outlier = best
    separation = full / rest_fold if rest_fold > 0 else float("inf")
    concordant = rest_fold <= 1.05
    others = [values[n] for n in names if n != outlier]
    return {
        "defined": True,
        "outlier": outlier,
        "outlier_is_low": values[outlier] < min(others),
        "rest_fold": rest_fold,
        "separation": separation,
        "concordant": bool(concordant),
        "bimodal": bool(concordant and separation >= 2.0),
    }


def _fold(values: dict[str, float]) -> dict:
    numbers = [values[n] for n in sorted(values)]
    finite = [v for v in numbers if np.isfinite(v)]
    entry = {"values": {k: float(v) for k, v in sorted(values.items())}, "n": len(numbers)}
    if len(finite) != len(numbers):
        entry.update({"min": None, "max": None, "fold": None,
                      "why": "a non-finite value — no fold"})
        entry["shape"] = _shape(values)
        return entry
    lo, hi = min(numbers), max(numbers)
    entry["min"], entry["max"] = float(lo), float(hi)
    entry["fold"] = float(hi / lo) if lo > 0 else None
    if entry["fold"] is None:
        entry["why"] = "min <= 0 — a fold is not defined; read min and max"
    entry["shape"] = _shape(values)
    return entry


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_targets(path: str):
    raw = np.load(path, allow_pickle=False)
    index = {p: i for i, p in enumerate(np.asarray(raw["patient_ids"]).astype(str))}
    keep = np.isin(np.asarray(raw["target_groups"]).astype(str), list(UNTRAINED_40))
    return index, np.asarray(raw["scores"], dtype=np.float64)[:, keep]


def score_repeat(path: str, index, target_values, seed: int = 42) -> dict:
    """Every statistic on every view, raw and residualised, for one exported repeat."""
    raw = np.load(path, allow_pickle=False)
    test = np.asarray(raw["split"]).astype(str) == "test"
    ids = np.asarray(raw["patient_ids"]).astype(str)[test]
    cancers = np.asarray(raw["cancers"]).astype(str)[test]
    tss, _ = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    ry = cross_fitted_residuals(target_values[[index[i] for i in ids]], design, seed=seed)

    blocks: dict[str, dict[str, np.ndarray]] = {"raw": {}, "residualised": {}}
    record: dict = {"path": str(Path(path).resolve()), "n_test": int(test.sum()),
                    "sha256": _sha256(path), "views": {}}
    for view in VIEWS:
        x = np.asarray(raw[view], dtype=np.float64)[test]
        blocks["raw"][view] = x
        blocks["residualised"][view] = cross_fitted_residuals(x, design, seed=seed)
        record["views"][view] = {
            block: {name: float(fn(blocks[block][view])) for name, fn in STATISTICS.items()}
            for block in ("raw", "residualised")}

    record["views"][LIDAR_VIEW] = {
        block: {"LiDAR": float(lidar([blocks[block]["wsi_biology"],
                                      blocks[block]["rna_biology"]],
                                     delta=LIDAR_DELTA)["lidar"])}
        for block in ("raw", "residualised")}

    # The co-measured comparator. §4.1's argument is the ASYMMETRY between these
    # two quantities on the same five runs, so the channel is measured here too
    # rather than quoted from the other table.
    record["channel_top_cca_16"] = float(
        top_canonical_correlation(blocks["residualised"]["wsi_biology"], ry, n_components=16))
    return record


def collect(reps: dict[str, dict]) -> dict:
    """{view: {block: {statistic: fold entry}}} across the repeats."""
    out: dict[str, dict] = {}
    for view in list(VIEWS) + [LIDAR_VIEW]:
        out[view] = {}
        for block in ("raw", "residualised"):
            names = sorted({s for r in reps.values() for s in r["views"][view][block]})
            out[view][block] = {
                stat: _fold({rep: r["views"][view][block][stat] for rep, r in reps.items()})
                for stat in names}
    out["channel"] = {"residualised": {"top_cca_16": _fold(
        {rep: r["channel_top_cca_16"] for rep, r in reps.items()})}}
    return out


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reps", default="/home/ubuntu/e0_run/d1_envelope/rep*.npz")
    ap.add_argument("--targets", default="/home/ubuntu/e0_run/data/frozen_rna_targets.npz")
    ap.add_argument("--output", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    paths = sorted(glob.glob(args.reps))
    if not paths:
        raise SystemExit(f"no repeats matched {args.reps!r}")
    index, target_values = load_targets(args.targets)

    reps = {}
    for path in paths:
        name = Path(path).stem
        reps[name] = score_repeat(path, index, target_values, seed=args.seed)
        m = reps[name]["views"]["wsi_biology"]
        print(f"{name:>6}  R1_raw={m['raw']['R1']:8.3f}  R1_res={m['residualised']['R1']:8.3f}  "
              f"R3_res={m['residualised']['R3']:8.3f}  "
              f"RankMe_raw={m['raw']['RankMe']:8.3f}  "
              f"channel={reps[name]['channel_top_cca_16']:.4f}", flush=True)

    out = {
        "_": ("Retraining floors for every statistic, view and block recoverable from the "
              "five exported same-seed repeats. A FLOOR TWICE OVER (stable arm; same-seed "
              "repeats exclude seed variation entirely), n = 5, no interval. Not a "
              "distribution."),
        "config": {"seed": args.seed, "reps_glob": args.reps, "targets": args.targets,
                   "alpha_index_range": list(ALPHA_INDEX_RANGE), "lidar_delta": LIDAR_DELTA,
                   "targets_groups": list(UNTRAINED_40),
                   "statistic_sites": {
                       "R1/R2/R3": "v2/calibra/spectral.py RANK_VARIANTS",
                       "channel": "v2/calibra/spectral.py top_canonical_correlation",
                       "PR/PR_rownorm/RankMe/stable_rank/alpha_req/LiDAR":
                           "v2/research/rebase/p2/p2_competing_metrics.py",
                       "hard_rank": "numpy.linalg.matrix_rank — NOT an effective rank"}},
        "shape_rule": ("outlier = the run whose removal minimises the remaining four's fold; "
                       "concordant = that fold <= 1.05; bimodal = concordant and "
                       "full fold / rest fold >= 2.0"),
        "absent_blocks": ABSENT_BLOCKS,
        "reps": reps,
        "floors": collect(reps),
    }
    path = Path(args.output).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")

    print("\nfloor (max/min over the five repeats)")
    print(f"{'view':>16} {'block':>13} {'statistic':>18} {'min':>10} {'max':>10} "
          f"{'fold':>8}  shape")
    for view, blocks in out["floors"].items():
        for block, stats in blocks.items():
            for stat, entry in stats.items():
                fold = "n/a" if entry["fold"] is None else f"{entry['fold']:.3f}x"
                shape = entry["shape"]
                tag = ("undefined" if not shape.get("defined") else
                       (f"BIMODAL outlier={shape['outlier']}" if shape["bimodal"]
                        else f"not bimodal (outlier={shape['outlier']}, "
                             f"rest={shape['rest_fold']:.3f}x, sep={shape['separation']:.2f})"))
                print(f"{view:>16} {block:>13} {stat:>18} "
                      f"{entry['min']:>10.4f} {entry['max']:>10.4f} {fold:>8}  {tag}")
    print("\nwrote", path)
    return out


if __name__ == "__main__":
    main()
