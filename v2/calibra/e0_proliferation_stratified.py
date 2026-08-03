"""E0 sensitivity analysis: stratify the responsive arm by proliferation loading.

This is `claim_guards.proliferation_deflation`'s remedy #2, run on the claim the
blocker is actually attached to. The blocker says the responsive arm is selected
on *having a detectable transcriptional effect*, which enriches for essential,
core-machinery, ribosome and cell-cycle genes -- so E0's responsive-minus-
nonresponsive gap could be proliferation-in-cell-lines matching proliferation-in-
tumours and would look identical to a real transfer result.

WHAT THIS IS NOT. It is not a re-certification of E0. It does not run E0's gate
ledger, and it must not be quoted as though it did. It imports and calls E0's own
``_arm_result`` and ``_decision`` verbatim so that the statistic, the PC1
stripping, the gene universe, the Haar null and the bootstrap are identical to
the run it is a sensitivity analysis of; the gated runner itself is untouched.

THE PLACEBO IS NOT OPTIONAL. Dropping proliferation-targeting perturbations
shrinks the responsive arm and re-draws the n-match against the control, so a
moving gap could be sample size rather than proliferation. ``responsive_placebo``
drops the same NUMBER of perturbations at random. Read it as:

  nonprolif gap ~= placebo gap  ->  the alignment is not proliferation-specific
  nonprolif gap <<  placebo gap  ->  the alignment WAS proliferation
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from .e0_basis_transfer import (TRANSFORM_SEED_OFFSET, TransferConfig, _align, _arm_result, _decision,
                                _energy_arms, _load_perturbation, _load_tcga, _overlap,
                                _restrict_tcga_to_registry, _rezscore, _right_svd, _split_half,
                                _unavailable_arm)

# Distinct offsets so the Haar draws and bootstrap resamples of the new arms are
# independent of E0's own arms and of each other. The two reference arms keep the
# offsets E0 itself uses, so `responsive_matched` and `nonresponsive` reproduce
# E0's own numbers rather than merely resembling them.
STRATIFIED_SEED_OFFSET = {"responsive_nonprolif": 52_000, "responsive_placebo": 65_000,
                          "responsive_prolif_only": 78_000}
# The control re-matched to the stratified arms' row count. Without it every
# stratified comparison is 869-vs-956 rows, and E0's `_decision` correctly
# refuses to certify an arm pair that is not n-matched -- so the decision would
# read False for a bookkeeping reason and look exactly like a scientific failure.
ARM_SEED_OFFSET = {"responsive_matched": 26_000, "nonresponsive": 39_000,
                   "nonresponsive_strata_matched": 91_000, **STRATIFIED_SEED_OFFSET}


def proliferation_gene_set(annotation_path: Path) -> set[str]:
    """Genes flagged proliferation by ``build_gene_annotations`` (MSigDB Hallmark union)."""
    table = (pd.read_parquet(annotation_path) if str(annotation_path).endswith(".parquet")
             else pd.read_csv(annotation_path))
    gene_column = next(c for c in ("gene", "gene_symbol", "symbol") if c in table.columns)
    flag = pd.to_numeric(table["proliferation_loading"], errors="coerce")
    return {str(g).strip().upper() for g, v in zip(table[gene_column], flag) if np.isfinite(v) and v > 0}


def stratify_responsive(responsive_rows: np.ndarray, targets: list[str | None] | None,
                        proliferation: set[str], *, seed: int) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    """Split the responsive arm by whether its TARGET gene is a proliferation gene.

    ``responsive_rows`` indexes the perturbation matrix; ``targets`` is row-aligned
    with that matrix. A perturbation whose target could not be parsed is treated as
    UNKNOWN and excluded from both strata rather than assumed non-proliferation --
    quietly banking unknowns as non-proliferation would inflate exactly the arm
    whose survival is the result.
    """
    if not targets:
        return {}, {"stratification_status": "unavailable_no_parsed_targets"}
    labels = np.asarray([targets[int(i)] if targets[int(i)] is not None else "" for i in responsive_rows])
    known = labels != ""
    is_prolif = np.asarray([label.upper() in proliferation for label in labels])
    nonprolif = responsive_rows[known & ~is_prolif]
    prolif_only = responsive_rows[known & is_prolif]
    n_dropped = int(len(responsive_rows) - len(nonprolif))
    # SIZE-matched placebo: drop the same NUMBER of perturbations, chosen at random
    # from the whole responsive arm. Same n, same n-match against the control, and
    # no relationship to proliferation.
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(len(responsive_rows), size=len(responsive_rows) - n_dropped, replace=False))
    placebo = responsive_rows[keep]
    if len(nonprolif) != len(placebo):                     # defensive; the draw above fixes the size
        raise ValueError(f"placebo {len(placebo)} is not size-matched to nonprolif {len(nonprolif)}")
    meta = {"stratification_status": "available", "n_responsive_in": int(len(responsive_rows)),
            "n_target_unknown": int((~known).sum()), "n_proliferation_targets": int(len(prolif_only)),
            "n_nonproliferation_targets": int(len(nonprolif)),
            "proliferation_target_fraction": float(len(prolif_only) / max(1, int(known.sum()))),
            "n_dropped_for_placebo": n_dropped, "placebo_seed": int(seed),
            "n_placebo": int(len(placebo)), "placebo_is_size_matched": bool(len(placebo) == len(nonprolif))}
    return {"responsive_nonprolif": nonprolif, "responsive_placebo": placebo,
            "responsive_prolif_only": prolif_only}, meta


def _gap(arm: dict[str, object], control: dict[str, object], k: int) -> float:
    key = f"k{k}"
    if not arm.get("available") or not control.get("available"):
        return float("nan")
    return float(arm[key]["pc1_removed_overlap"]) - float(control[key]["pc1_removed_overlap"])


def run_context(name: str, p, t, *, device: torch.device, cfg: TransferConfig, draws: int,
                seed: int, proliferation: set[str]) -> dict[str, object]:
    """Every arm through the identical E0 pipeline; only the rows of P differ."""
    px, tx, genes, join = _align(p, t)
    pt, tt = torch.as_tensor(px, device=device), torch.as_tensor(tx, device=device)
    vt, _ = _right_svd(tt, cfg.q, seed + 1, min_q=cfg.min_q)
    # Same split-half positive ceiling E0 uses, computed once with E0's seed
    # convention so `normalised_alignment` is on E0's scale.
    va, vb, split_meta = _split_half(tt, t.groups or [], cfg.q, seed + 2, min_q=cfg.min_q)
    ceiling = {k: _overlap(va, vb, k, cfg.primary_offset) for k in cfg.ks}

    arm_rows, arm_meta = _energy_arms(p.energy_p, cfg=cfg, seed=seed + 11)
    if not arm_rows:
        return {"context": name, "status": arm_meta.get("arm_status", "unavailable"), "arm_split": arm_meta}
    strata, strat_meta = stratify_responsive(arm_rows["responsive_matched"], p.targets, proliferation,
                                             seed=20260803)
    # Re-match the control to the stratified arms. E0's `_decision` requires
    # `arms_are_n_matched`; comparing an 869-row stratified arm against the
    # 956-row control returns False for bookkeeping, not biology.
    control_rows = arm_rows["nonresponsive"]
    strata_n = len(strata.get("responsive_nonprolif", []))
    if strata_n and len(control_rows) > strata_n:
        control_rows = np.sort(np.random.default_rng(20260803 + 1).choice(
            control_rows, size=strata_n, replace=False))
    strat_meta["n_nonresponsive_strata_matched"] = int(len(control_rows))
    arms: dict[str, object] = {}
    for arm_name, rows in [("responsive_matched", arm_rows["responsive_matched"]),
                           ("nonresponsive", arm_rows["nonresponsive"]),
                           ("nonresponsive_strata_matched", control_rows),
                           *[(n, strata.get(n)) for n in STRATIFIED_SEED_OFFSET]]:
        if rows is None or len(rows) == 0:
            arms[arm_name] = _unavailable_arm(arm_name, str(strat_meta.get("stratification_status",
                                                                           "unavailable_empty_arm")),
                                              0 if rows is None else len(rows))
            continue
        offset = ARM_SEED_OFFSET[arm_name]
        arm_x = _rezscore(pt[torch.as_tensor(np.asarray(rows), device=device)])
        arms[arm_name], _, _ = _arm_result(arm_x, vt, ceiling=ceiling, n_genes=len(genes), cfg=cfg,
                                           draws=draws, seed=seed + offset, arm=arm_name)
        print(f"  [{name}] arm {arm_name}: n={len(rows)} "
              f"k10_overlap={arms[arm_name].get('k10', {}).get('pc1_removed_overlap', float('nan')):.4f}",
              flush=True)

    # E0's own reference pair stays at its own n (reproducing E0_RESULT exactly);
    # every stratified arm is judged against the control re-matched to ITS n.
    control = arms["nonresponsive"]
    strata_control = arms["nonresponsive_strata_matched"]
    per_k: dict[str, object] = {}
    for k in cfg.ks:
        reference = _gap(arms["responsive_matched"], control, k)
        entry = {"gap_responsive_matched": reference,
                 "decision_responsive_matched": _decision(arms["responsive_matched"], control, k)}
        for arm_name in STRATIFIED_SEED_OFFSET:
            gap = _gap(arms[arm_name], strata_control, k)
            entry[f"gap_{arm_name}"] = gap
            entry[f"retention_{arm_name}"] = float(gap / reference) if np.isfinite(reference) and abs(reference) > 1e-12 else float("nan")
            entry[f"decision_{arm_name}"] = _decision(arms[arm_name], strata_control, k)
        per_k[f"k{k}"] = entry
    return {"context": name, "status": "available", "n_shared_genes": len(genes), "join": join,
            "arm_split": arm_meta, "stratification": strat_meta, "split": split_meta,
            "config": {"ks": list(cfg.ks), "q": cfg.q, "primary_offset": cfg.primary_offset,
                       "draws": draws, "bootstrap_draws": cfg.bootstrap_draws},
            "arms": arms, "by_k": per_k}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--perturbation", required=True, help="K562 gwps normalized bulk h5ad")
    parser.add_argument("--tcga", required=True)
    parser.add_argument("--tcga-registry", required=True)
    parser.add_argument("--annotations", required=True, help="per-gene annotation parquet")
    parser.add_argument("--output", required=True)
    parser.add_argument("--transforms", default="signed_log1p,clip_log1p")
    parser.add_argument("--draws", type=int, default=100)
    parser.add_argument("--bootstrap-draws", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    if args.draws < 100 or args.bootstrap_draws < 200:
        parser.error("E0's predeclared floors: --draws >=100 and --bootstrap-draws >=200")

    device = torch.device(args.device if torch.cuda.is_available() or args.device == "cpu" else "cpu")
    proliferation = proliferation_gene_set(Path(args.annotations))
    print(f"[e0-prolif] device={device} proliferation_genes={len(proliferation)}", flush=True)
    perturbation = _load_perturbation(Path(args.perturbation))
    print(f"[e0-prolif] perturbation rows={perturbation.x.shape[0]} "
          f"targets={perturbation.meta.get('guide_target_status')}", flush=True)
    cfg = TransferConfig(draws=args.draws, bootstrap_draws=args.bootstrap_draws)
    results: dict[str, object] = {"perturbation": perturbation.meta,
                                  "n_proliferation_genes": len(proliferation),
                                  "note": ("sensitivity analysis on E0's decision rule; imports E0's own "
                                           "_arm_result/_decision; does NOT re-run E0's gate ledger"),
                                  "contexts": {}}
    for transform in args.transforms.split(","):
        tcga = _restrict_tcga_to_registry(_load_tcga(Path(args.tcga), transform), Path(args.tcga_registry))
        name = f"K562_{transform}"
        print(f"[e0-prolif] context {name} patients={tcga.x.shape[0]}", flush=True)
        results["contexts"][name] = run_context(name, perturbation, tcga, device=device, cfg=cfg,
                                                draws=args.draws,
                                                seed=args.seed + TRANSFORM_SEED_OFFSET[transform],
                                                proliferation=proliferation)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
    print(f"[e0-prolif] done -> {out}", flush=True)


if __name__ == "__main__":
    main()
