"""P2: the LABELLED LINEAR PROBE that §2.5 and §6.2 both mark as never run.

Why this file exists
--------------------
`paper/P2_RANK_DRAFT.md` §2.5 ends with

    [STILL NOT MEASURED] -- a labelled linear probe on every artifact, which is the ground
    truth LiDAR and RankMe were validated against. Ours is a held-out canonical correlation
    against unsupervised molecular targets (§3.2), which is a different reference standard
    and is the one Zaiem et al. would attack.

§6.2 carries the same row and §6.3 concedes Zaiem et al. against us. A paper about whether a
label-free proxy tracks usable information has therefore never measured the labelled thing the
proxy is a proxy *for*. This module measures it, on the same twelve artifacts, the same three
co-trained views and the same cohort rule §4.1--§4.7 use, and scores it as a SELECTION RULE on
the same six matched pairs `p2_selection_rule.py` scores rank on.

Pre-registered in `NOTEBOOK_ENTRIES/PREDECLARED_p2_labelled_linear_probe_20260805T0040Z.md`,
committed before any statistic here was computed.

Two probes, and the reason there are two
----------------------------------------
* **Probe A -- cancer type, RAW representation.** The literal analogue of what RankMe and LiDAR
  validated against: a linear classifier fitted on frozen features, reported as out-of-fold
  accuracy. Labels come from the artifact's own `cancers` array; nothing is invented.
* **Probe B -- molecular / clinical endpoints, CONFOUND-ADJUSTED.** Matched to the construct
  this paper argues about. Labels come from E1's `e1_endpoint_labels.parquet`, built from the
  MC3 public MAF (`e1_clinical_endpoints_from_morphology_20260804T1930Z.md`).

They are different reference standards and they are allowed to disagree. When they do, both are
reported and neither is nominated the winner -- that is §4.6a's finding restated on the labelled
side.

Nothing statistical is implemented here
---------------------------------------
`_sklearn_oof`, `lda_oof_balanced_accuracy`, `balanced_accuracy`, `_stratified_folds`, `_encode`,
`global_permutations` and `null_summary` are imported from `calibra.confound_certificate` /
`calibra.nonlinear_confound_probe`; `evaluate_known_covariate` from
`calibra.known_covariate_control`; `confound_design` / `cross_fitted_residuals` /
`pooled_tissue_source_site` from `calibra.residualise`; `effective_rank` /
`top_canonical_correlation` from `calibra.spectral`. Re-deriving a statistic where an import
belonged is the most repeated defect on this project and it is what produced the withdrawn P1
cancer-type pair (`p1_cancer_type_pair_withdrawn_20260805T0230Z.md`).

Example
-------
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
    PYTHONPATH=$WS python -m morpheus.v2.research.rebase.p2.p2_labelled_probe \
        --artifacts H42=... I42=... P42=... F42=... REP1=... \
        --labels .../e1_endpoint_labels.parquet \
        --targets .../frozen_rna_targets.npz \
        --output P2_LABELLED_PROBE.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.confound_certificate import (_encode, _stratified_folds,
                                                      lda_oof_balanced_accuracy)
from morpheus.v2.calibra.known_covariate_control import evaluate_known_covariate
from morpheus.v2.calibra.nonlinear_confound_probe import (_sklearn_oof, global_permutations,
                                                          null_summary)
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import effective_rank, top_canonical_correlation

__all__ = ["logistic_balanced_accuracy_oof", "probe_cancer_type", "score_artifact", "summarise"]

VIEWS = ("wsi_biology", "rna_biology", "full_biology")

# The six matched pairs §4.6 scores rank on, in the same order and with the same labels.
PAIRS = [("D2", "s42", "H42", "I42"), ("D2", "s43", "H43", "I43"), ("D2", "s44", "H44", "I44"),
         ("D1", "s42", "P42", "F42"), ("D1", "s43", "P43", "F43"), ("D1", "s44", "P44", "F44")]

# The five same-seed retrains §4.1's 3.295x rank floor and 1.055x channel spread are measured on.
ENVELOPE = ["REP1", "REP2", "REP3", "REP4", "REP5"]

# Probe B endpoints. The band for TP53 is the one E1 pre-registered from the published
# literature (Kather 0.60-0.78; Noorbakhsh 0.65-0.80); the rest are EXPLORATORY -- they have no
# band, [0, 1] is a placeholder, and `passed` must not be read for them.
ENDPOINTS = {
    "mut_TP53": {"band": (0.60, 0.80), "exploratory": False},
    "grade_high": {"band": (0.0, 1.0), "exploratory": True},
    "stage_late": {"band": (0.0, 1.0), "exploratory": True},
    "mut_ATM": {"band": (0.0, 1.0), "exploratory": True},
    "mut_KMT2D": {"band": (0.0, 1.0), "exploratory": True},
    "mut_ARID1A": {"band": (0.0, 1.0), "exploratory": True},
}
PRIMARY_ENDPOINT = "mut_TP53"

# §3.2's molecular ground truth: the 40 targets neither arm was supervised on.
UNTRAINED40 = ("heldout_pathway", "immune_tme", "tumour_state")


# ---------------------------------------------------------------------------
# Probe A -- the RankMe/LiDAR analogue
# ---------------------------------------------------------------------------
def logistic_balanced_accuracy_oof(features, class_index, n_classes, *, n_splits: int = 5,
                                   seed: int = 42, folds=None, max_iter: int = 2000,
                                   C: float = 1.0) -> float:
    """Out-of-fold multinomial logistic-regression balanced accuracy.

    A thin `make_model` for `nonlinear_confound_probe._sklearn_oof`, which owns the fold loop,
    the train-fold-only standardisation and the call to `balanced_accuracy`. Nothing about the
    statistic is re-implemented here; the only new content is the estimator, which is the one
    RankMe and LiDAR's own validation used (a linear classifier on frozen features).

    `class_weight="balanced"` because the reported loss is balanced accuracy over 21 very
    unequal TCGA lineages; without it the classifier optimises plain accuracy and the probe
    reads low for a reason that has nothing to do with the representation. This mirrors the
    choice already made in `forest_balanced_accuracy_oof` and `rbf_svm_balanced_accuracy_oof`.
    """
    from sklearn.linear_model import LogisticRegression

    def _make():
        return LogisticRegression(C=float(C), max_iter=int(max_iter), class_weight="balanced",
                                  solver="lbfgs", random_state=int(seed), n_jobs=1)

    return _sklearn_oof(_make, features, class_index, n_classes, n_splits=n_splits, seed=seed,
                        folds=folds, standardise=True)


def probe_cancer_type(features, class_index, n_classes, *, folds, seed: int,
                      n_permutations_logistic: int, n_permutations_lda: int) -> dict:
    """Probe A on one feature block: two linear estimators, each against its own measured null.

    Two estimators on IDENTICAL folds, declared in advance, because
    `p1_cancer_type_pair_withdrawn_20260805T0230Z.md` measured this same statistic moving from
    0.16 to 0.73 across estimators on one cohort. A selection verdict carried by one estimator
    is not a verdict.
    """
    out = {"n_classes": int(n_classes), "chance_rate": 1.0 / float(n_classes)}

    t0 = time.time()
    out["logistic_balanced_accuracy"] = logistic_balanced_accuracy_oof(
        features, class_index, n_classes, folds=folds, seed=seed)
    out["logistic_seconds"] = round(time.time() - t0, 2)

    t0 = time.time()
    out["lda_balanced_accuracy"] = lda_oof_balanced_accuracy(
        features, class_index, n_classes, folds=folds, seed=seed)
    out["lda_seconds"] = round(time.time() - t0, 2)

    # Nulls: label permutation with NO stratification, which is what MEASURES chance. Same
    # estimator, same folds, same code path as the observed value.
    if n_permutations_lda:
        draws = global_permutations(class_index, n_permutations=n_permutations_lda, seed=seed)
        null = np.asarray([lda_oof_balanced_accuracy(features, d, n_classes, folds=folds, seed=seed)
                           for d in draws])
        out["lda_null"] = null_summary(null, out["lda_balanced_accuracy"])
    if n_permutations_logistic:
        draws = global_permutations(class_index, n_permutations=n_permutations_logistic, seed=seed)
        null = np.asarray([logistic_balanced_accuracy_oof(features, d, n_classes, folds=folds,
                                                          seed=seed) for d in draws])
        out["logistic_null"] = null_summary(null, out["logistic_balanced_accuracy"])
    return out


# ---------------------------------------------------------------------------
# plumbing
# ---------------------------------------------------------------------------
def load_targets(path, groups):
    """`p2_competing_metrics.load_targets`, verbatim in behaviour, so the two agree by rule."""
    raw = np.load(path, allow_pickle=False)
    ids = np.asarray(raw["patient_ids"]).astype(str)
    names = np.asarray(raw["target_names"]).astype(str)
    values = np.asarray(raw["scores"], dtype=np.float64)
    available = np.asarray(raw["target_groups"]).astype(str)
    keep = np.isin(available, list(groups)) if groups else ~np.char.startswith(names, "RANDOM_CONTROL__")
    return {p: i for i, p in enumerate(ids)}, values[:, keep]


def score_artifact(path, *, labels_table, targets_path, seed, views, endpoints,
                   n_permutations_logistic, n_permutations_lda, n_boot, n_permutations_auroc,
                   partition="test", null_views=None) -> dict:
    raw = np.load(path, allow_pickle=False)
    test = raw["split"].astype(str) == partition
    ids = raw["patient_ids"].astype(str)[test]
    cancers = raw["cancers"].astype(str)[test]

    site, _ = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": site}), ["cancer", "tss"])

    classes, class_index = _encode(cancers)
    n_classes = int(len(classes))
    folds = _stratified_folds(class_index, 5, seed)

    index, target_values = load_targets(targets_path, UNTRAINED40)
    rows = np.asarray([index[p] for p in ids])
    ry = cross_fitted_residuals(target_values[rows], design, seed=seed)

    label_rows = labels_table.set_index("patient_id").reindex(ids)

    entry = {"path": str(Path(path).resolve()), "n_test": int(test.sum()),
             "n_cancer_classes": n_classes, "n_confound_columns": int(design.shape[1]),
             "views": {}}
    for view in views:
        x = raw[view].astype(np.float64)[test]
        rx = cross_fitted_residuals(x, design, seed=seed)
        # The logistic permutation null costs ~5 s per draw against the LDA null's 0.04 s, so
        # it is computed only on the views named in `--null-views` (the readout view by
        # default) and its absence elsewhere is recorded rather than silently implied. The LDA
        # null is computed on every view and every block.
        logistic_draws = (n_permutations_logistic
                          if (null_views is None or view in null_views) else 0)
        block = {
            # re-quoted with the canonical functions so every column of the agreement table
            # comes off one cohort and one code path
            "effective_rank_raw": effective_rank(x),
            "effective_rank_residualised": effective_rank(rx),
            "channel_untrained40": top_canonical_correlation(rx, ry, n_components=16),
            # Probe A on the raw block -- the RankMe analogue
            "probe_A_cancer_type_raw": probe_cancer_type(
                x, class_index, n_classes, folds=folds, seed=seed,
                n_permutations_logistic=logistic_draws,
                n_permutations_lda=n_permutations_lda),
            # Probe A on the adjusted block -- the MUST-FAIL control. If cancer type survives
            # cancer+TSS residualisation well above chance, the adjustment is not doing what
            # §3.2 says and every downstream number in this run is void.
            "probe_A_cancer_type_adjusted": probe_cancer_type(
                rx, class_index, n_classes, folds=folds, seed=seed,
                n_permutations_logistic=0, n_permutations_lda=n_permutations_lda),
            "probe_B": {},
        }
        for endpoint in endpoints:
            spec = ENDPOINTS[endpoint]
            values = pd.to_numeric(label_rows[endpoint], errors="coerce").to_numpy()
            keep = np.isfinite(values)
            if keep.sum() < 60:
                block["probe_B"][endpoint] = {"status": "unavailable_insufficient_labelled_cases",
                                              "n_labelled": int(keep.sum())}
                continue
            result = evaluate_known_covariate(
                x[keep], values[keep] > 0.5, cancers[keep], ids[keep],
                expected_low=spec["band"][0], expected_high=spec["band"][1],
                seed=seed, n_boot=n_boot, n_permutations=n_permutations_auroc,
                residualise=True, covariate_name=endpoint)
            result["exploratory_no_preregistered_band"] = bool(spec["exploratory"])
            result.pop("per_cancer", None)
            block["probe_B"][endpoint] = result
        entry["views"][view] = block
        print(f"    {view:14s} rank_res={block['effective_rank_residualised']:8.3f} "
              f"channel={block['channel_untrained40']:.4f} "
              f"probeA_logit={block['probe_A_cancer_type_raw']['logistic_balanced_accuracy']:.4f} "
              f"probeA_lda={block['probe_A_cancer_type_raw']['lda_balanced_accuracy']:.4f} "
              f"probeA_adj_lda={block['probe_A_cancer_type_adjusted']['lda_balanced_accuracy']:.4f} "
              f"TP53={block['probe_B'].get(PRIMARY_ENDPOINT, {}).get('within_cancer_auroc', float('nan')):.4f}",
              flush=True)
    return entry


# ---------------------------------------------------------------------------
# the agreement arithmetic -- the question the paper needs answered
# ---------------------------------------------------------------------------
def _floor(values):
    """Spread over the five same-seed retrains, both ways §4.1 could be read.

    `ratio` is max/min, the form §4.1 quotes for rank (multiplicative, spanning 6.4-34.1).
    `spread` is max-min, the form appropriate to a bounded rate such as balanced accuracy or
    AUROC. Both are reported for every statistic and the verdict is stated under both, because
    which one applies to a bounded rate is a choice and the paper's own lesson is that a choice
    made silently is a defect.
    """
    values = np.asarray([v for v in values if np.isfinite(v)], dtype=np.float64)
    if len(values) < 2:
        return {"n": int(len(values)), "ratio": float("nan"), "spread": float("nan")}
    return {"n": int(len(values)), "min": float(values.min()), "max": float(values.max()),
            "ratio": float(values.max() / values.min()) if values.min() > 0 else float("nan"),
            "spread": float(values.max() - values.min())}


def _get(entry, view, path):
    node = entry["views"][view]
    for key in path:
        if node is None or key not in node:
            return float("nan")
        node = node[key]
    return float(node) if np.isscalar(node) or isinstance(node, (int, float)) else float("nan")


STATISTICS = {
    "effective_rank_residualised": ("effective_rank_residualised",),
    "channel_untrained40": ("channel_untrained40",),
    "probeA_logistic": ("probe_A_cancer_type_raw", "logistic_balanced_accuracy"),
    "probeA_lda": ("probe_A_cancer_type_raw", "lda_balanced_accuracy"),
    "probeB_TP53_within_cancer": ("probe_B", "mut_TP53", "within_cancer_auroc"),
    "probeB_TP53_pooled": ("probe_B", "mut_TP53", "pooled_auroc"),
}


def summarise(data, views=VIEWS) -> dict:
    """For every statistic, every view and every pair: which arm does it pick, and is it resolvable?

    "Resolvable" applies §4.1's own criterion to each statistic in turn -- the between-arm
    difference must exceed the SAME statistic's five-repeat same-seed floor on the SAME view.
    A difference inside the floor is UNRESOLVED and is scored as neither an agreement nor a
    disagreement, exactly as the predeclaration fixed.
    """
    floors = {}
    for view in views:
        floors[view] = {}
        for name, path in STATISTICS.items():
            floors[view][name] = _floor([_get(data[label], view, path)
                                         for label in ENVELOPE if label in data])

    table = {}
    for view in views:
        rows = []
        for experiment, seed, a, b in PAIRS:
            if a not in data or b not in data:
                continue
            row = {"experiment": experiment, "seed": seed, "arm_a": a, "arm_b": b, "statistics": {}}
            for name, path in STATISTICS.items():
                va, vb = _get(data[a], view, path), _get(data[b], view, path)
                delta = va - vb
                floor = floors[view][name]
                ratio = max(va, vb) / min(va, vb) if min(va, vb) > 0 else float("nan")
                resolvable_spread = bool(np.isfinite(floor["spread"]) and abs(delta) > floor["spread"])
                resolvable_ratio = bool(np.isfinite(floor["ratio"]) and np.isfinite(ratio)
                                        and ratio > floor["ratio"])
                row["statistics"][name] = {
                    "value_a": va, "value_b": vb, "delta": delta,
                    "winner": (a if delta > 0 else b) if np.isfinite(delta) and delta != 0 else "tie",
                    "between_arm_ratio": ratio,
                    "floor_spread": floor["spread"], "floor_ratio": floor["ratio"],
                    "resolvable_by_spread": resolvable_spread,
                    "resolvable_by_ratio": resolvable_ratio,
                }
            rows.append(row)
        table[view] = rows

    # Agreement counts. Reference = the channel (the paper's published ground truth); the
    # question §2.5 asks is whether the LABELLED probe says the same thing.
    agreement = {}
    for view in views:
        agreement[view] = {}
        for name in STATISTICS:
            if name == "channel_untrained40":
                continue
            agree = disagree = unresolved_self = 0
            detail = []
            for row in table[view]:
                ref = row["statistics"]["channel_untrained40"]["winner"]
                mine = row["statistics"][name]["winner"]
                res = row["statistics"][name]["resolvable_by_spread"]
                mark = "OK" if mine == ref else "MISS"
                if mine == ref:
                    agree += 1
                else:
                    disagree += 1
                if not res:
                    unresolved_self += 1
                detail.append({"pair": f"{row['experiment']} {row['seed']}", "winner": mine,
                               "channel_winner": ref, "mark": mark, "resolvable": res})
            agreement[view][name] = {"agrees_with_channel": agree, "disagrees_with_channel": disagree,
                                     "n_pairs": agree + disagree,
                                     "n_unresolved_against_own_floor": unresolved_self,
                                     "detail": detail}

    # The pairwise question the predeclaration's outcome D turns on: where the probe and the
    # channel disagree, which one does RANK side with?
    conflicts = []
    for view in views:
        for row in table[view]:
            ref = row["statistics"]["channel_untrained40"]["winner"]
            rank = row["statistics"]["effective_rank_residualised"]["winner"]
            for name in ("probeA_logistic", "probeA_lda", "probeB_TP53_within_cancer"):
                mine = row["statistics"][name]["winner"]
                if mine != ref:
                    conflicts.append({
                        "view": view, "pair": f"{row['experiment']} {row['seed']}",
                        "probe": name, "probe_winner": mine, "channel_winner": ref,
                        "rank_winner": rank, "rank_sides_with": "probe" if rank == mine else "channel",
                        "probe_resolvable_by_spread": row["statistics"][name]["resolvable_by_spread"],
                        "probe_delta": row["statistics"][name]["delta"],
                        "probe_floor_spread": row["statistics"][name]["floor_spread"],
                    })
    return {"floors": floors, "pairs": table, "agreement": agreement,
            "channel_probe_conflicts": conflicts}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifacts", nargs="+", required=True, help="LABEL=path")
    ap.add_argument("--labels", required=True, help="E1 endpoint label parquet")
    ap.add_argument("--targets", required=True, help="frozen_rna_targets.npz")
    ap.add_argument("--output", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--views", nargs="+", default=list(VIEWS))
    ap.add_argument("--endpoints", nargs="+", default=[PRIMARY_ENDPOINT])
    ap.add_argument("--wsi-only-endpoints", nargs="+", default=[],
                    help="endpoints scored on the FIRST view only, to bound cost")
    ap.add_argument("--n-permutations-logistic", type=int, default=25)
    ap.add_argument("--null-views", nargs="+", default=None,
                    help="views that get the (expensive) logistic permutation null; "
                         "default = every view named by --views")
    ap.add_argument("--n-permutations-lda", type=int, default=200)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--n-permutations-auroc", type=int, default=1000)
    a = ap.parse_args(argv)

    labels_table = (pd.read_parquet(a.labels) if str(a.labels).lower().endswith(".parquet")
                    else pd.read_csv(a.labels))
    out = {"_config": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(a).items()}}
    out["_config"]["targets"] = str(Path(a.targets).resolve())
    out["_config"]["labels"] = str(Path(a.labels).resolve())

    for item in a.artifacts:
        label, path = item.split("=", 1)
        print(f"[{label}] {path}", flush=True)
        entry = score_artifact(
            path, labels_table=labels_table, targets_path=a.targets, seed=a.seed,
            views=a.views, endpoints=list(a.endpoints),
            n_permutations_logistic=a.n_permutations_logistic,
            n_permutations_lda=a.n_permutations_lda, n_boot=a.n_boot,
            n_permutations_auroc=a.n_permutations_auroc, null_views=a.null_views)
        # extra endpoints, first view only, to bound cost -- declared, not silent
        if a.wsi_only_endpoints:
            extra = score_artifact(path, labels_table=labels_table, targets_path=a.targets,
                                   seed=a.seed, views=a.views[:1],
                                   endpoints=list(a.wsi_only_endpoints),
                                   n_permutations_logistic=0, n_permutations_lda=0,
                                   n_boot=a.n_boot, n_permutations_auroc=a.n_permutations_auroc)
            entry["views"][a.views[0]]["probe_B"].update(extra["views"][a.views[0]]["probe_B"])
        out[label] = entry

    out["_summary"] = summarise({k: v for k, v in out.items() if not k.startswith("_")},
                                views=a.views)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    Path(a.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("wrote", a.output)
    return out


if __name__ == "__main__":
    main()
