"""Leave-sites-out generalisation: does the channel survive on sites never seen?

Predeclared in ``NOTEBOOK_ENTRIES/PREDECLARED_leave_sites_out_20260804T1745Z.md``
before the first statistic was run.

Why this test exists
--------------------
``claim_guards.no_external_cohort`` blocks the two strongest claim types. Its stated
mechanism is that "a site-specific artifact would survive the current checks intact".
The *purpose* of an external cohort is to ask whether a result survives different
acquisition conditions, and that question can be asked today by holding out whole
tissue source sites. **This does not discharge the blocker** — same institution
network, same protocols, same era, same encoder — and nothing here edits
``claim_guards.py``. It is evidence about acquisition-condition generalisation.

Site is nested within cancer, and that dictates the folds
---------------------------------------------------------
TCGA assigns a TSS code per submitting site *per disease*. Measured on the frozen
artifact: 0 of 352 ``test``-partition sites contribute more than one cancer, and 0
of them appear in ``train``/``val``. Two consequences:

* Because the existing split holds out *cancers*, the encoder has already never seen
  any ``test`` site — so the existing channel number is already on encoder-unseen
  sites, but confounded with unseen cancers and unreadable as a site result.
* A naive whole-cohort site holdout would remove whole cancers and silently re-run
  the held-out-cancer test under a new name.

Folds are therefore built **within cancer**: every fold holds out whole sites, and
each cancer contributes sites to as many folds as it has sites. The two tests are
named ``held-out-cancer`` and ``held-out-site`` throughout and are never merged.

Why the adjusted arm residualises each block separately
-------------------------------------------------------
A held-out site has no column in a confound design fit on the training sites, so it
cannot be adjusted from there. Each block is adjusted by its own ``cancer + pooled
TSS`` structure, which is exactly what one does with a real external cohort. Because
site is nested in cancer the TSS dummies span the cancer dummies, so this arm removes
cancer too and measures within-site, within-cancer covariation only — the most
stringent reading available, and reported as such.

Why the matched random comparator is not optional
-------------------------------------------------
A held-out-site channel below the pooled in-sample number confounds two things: the
cost of *site shift*, and the cost of fitting the read-out out-of-sample at all. The
comparator re-runs the identical estimator on a random patient split with identical
per-(cancer, fold) counts. Without it a drop cannot be attributed to sites.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from .spectral import heldout_cca_projection, paired_absolute_correlation

__all__ = ["site_folds", "matched_random_folds", "fold_composition",
           "evaluate_fold", "leave_sites_out_channel"]


def site_folds(sites: np.ndarray, cancers: np.ndarray, *, n_folds: int = 5) -> np.ndarray:
    """Assign WHOLE sites to folds, balanced within cancer. Deterministic and seed-free.

    For each cancer independently, its sites are sorted by descending patient count
    (ties broken by site code, so the result does not depend on input order) and each
    site is handed to the fold currently holding fewest patients *of that cancer*. A
    site never straddles two folds, which is the whole point: holding out patients
    would leave the same scanner and the same protocol on both sides of the split.

    Balancing within cancer rather than globally is what keeps this a *site* test.
    A globally balanced assignment would put whole cancers in one fold, because TSS
    is nested in cancer, and would re-run the held-out-cancer test instead.

    Seed-free by design: there is no random choice to make, so there is no seed a
    later run could differ on and no split a caller could re-roll until it liked it.
    """
    sites = np.asarray(sites).astype(str)
    cancers = np.asarray(cancers).astype(str)
    if len(sites) != len(cancers):
        raise ValueError("sites and cancers must be aligned")
    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")
    fold = np.full(len(sites), -1, dtype=np.int64)
    for cancer in sorted(set(cancers.tolist())):
        in_cancer = np.flatnonzero(cancers == cancer)
        counts: dict[str, int] = {}
        for site in sites[in_cancer]:
            counts[site] = counts.get(site, 0) + 1
        ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        load = np.zeros(n_folds, dtype=np.int64)
        for site, size in ordered:
            target = int(np.argmin(load))
            load[target] += size
            fold[in_cancer[sites[in_cancer] == site]] = target
    if (fold < 0).any():
        raise AssertionError("every patient must receive a fold")
    return fold


def matched_random_folds(cancers: np.ndarray, fold: np.ndarray, *, seed: int) -> np.ndarray:
    """A random patient split with the SAME per-(cancer, fold) counts as ``fold``.

    Patients are permuted within cancer and re-dealt into the fold sizes the site
    assignment produced, so cancer composition and every fold size are matched
    exactly and the only thing that changes is whether a fold respects site
    boundaries. That is the contrast the comparator has to isolate.
    """
    cancers = np.asarray(cancers).astype(str)
    fold = np.asarray(fold, dtype=np.int64)
    rng = np.random.default_rng(seed)
    matched = np.full(len(fold), -1, dtype=np.int64)
    for cancer in sorted(set(cancers.tolist())):
        rows = np.flatnonzero(cancers == cancer)
        shuffled = rng.permutation(rows)
        matched[shuffled] = fold[rows]
    if (matched < 0).any():
        raise AssertionError("every patient must receive a matched fold")
    return matched


def fold_composition(fold: np.ndarray, sites: np.ndarray, cancers: np.ndarray) -> list[dict]:
    """Per-fold site, patient and cancer counts — TSS is heavily imbalanced, so report it."""
    fold = np.asarray(fold, dtype=np.int64)
    sites = np.asarray(sites).astype(str)
    cancers = np.asarray(cancers).astype(str)
    rows = []
    for f in sorted(set(fold.tolist())):
        held = fold == f
        held_sites = sorted(set(sites[held].tolist()))
        train_cancers = set(cancers[~held].tolist())
        shared = sorted(set(cancers[held].tolist()) & train_cancers)
        sizes = sorted((int((sites[held] == s).sum()) for s in held_sites), reverse=True)
        rows.append({
            "fold": int(f),
            "n_heldout_patients": int(held.sum()),
            "n_heldout_sites": len(held_sites),
            "n_training_patients": int((~held).sum()),
            "n_training_sites": len(set(sites[~held].tolist())),
            "n_heldout_cancers": len(set(cancers[held].tolist())),
            "n_cancers_shared_with_training": len(shared),
            "largest_heldout_site_n": sizes[0] if sizes else 0,
            "largest_heldout_site_share": float(sizes[0] / max(int(held.sum()), 1)) if sizes else 0.0,
            "heldout_site_size_median": float(np.median(sizes)) if sizes else 0.0,
        })
    return rows


def _adjusted(x: np.ndarray, y: np.ndarray, cancers: np.ndarray, sites: np.ndarray,
              rows: np.ndarray, *, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """CALIBRA's own ``cancer + pooled TSS`` adjustment, applied WITHIN one block."""
    import pandas as pd

    frame = pd.DataFrame({"cancer": cancers[rows], "tss": sites[rows]})
    design = confound_design(frame, ["cancer", "tss"])
    return (cross_fitted_residuals(x[rows], design, seed=seed),
            cross_fitted_residuals(y[rows], design, seed=seed))


def evaluate_fold(x: np.ndarray, y: np.ndarray, cancers: np.ndarray, sites: np.ndarray,
                  pooled_sites: np.ndarray, held: np.ndarray, *, n_components: int,
                  adjust: bool, n_permutations: int, n_boot: int, seed: int) -> dict:
    """Fit the channel on the training sites, score it on the held-out sites.

    The null permutes the cross-modal pairing WITHIN CANCER among the held-out rows
    and rescores through the *same* fitted directions — refitting per permutation
    would test a different hypothesis. The primary interval resamples held-out
    **sites**, not patients, because the question is whether another draw of sites
    would give this; a patient bootstrap treats two patients from one scanner as two
    independent facts and is reported only alongside.
    """
    train_rows = np.flatnonzero(~held)
    test_rows = np.flatnonzero(held)
    if len(train_rows) < 10 or len(test_rows) < 10:
        return {"status": "unavailable_fold_too_small",
                "n_heldout": int(len(test_rows)), "n_training": int(len(train_rows))}

    if adjust:
        x_train, y_train = _adjusted(x, y, cancers, pooled_sites, train_rows, seed=seed)
        x_test, y_test = _adjusted(x, y, cancers, pooled_sites, test_rows, seed=seed)
        # The two blocks are residualised separately but stay in the same coordinate
        # space, so the train-fitted projection applies to the test block unchanged.
        stacked_x = np.zeros_like(np.asarray(x, dtype=np.float64))
        stacked_y = np.zeros((len(x), np.asarray(y).shape[1]), dtype=np.float64)
        stacked_x[train_rows], stacked_x[test_rows] = x_train, x_test
        stacked_y[train_rows], stacked_y[test_rows] = y_train, y_test
        x_use, y_use = stacked_x, stacked_y
    else:
        x_use, y_use = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)

    px, py = heldout_cca_projection(x_use, y_use, train_rows, test_rows,
                                    n_components=n_components)
    if px.size == 0:
        return {"status": "unavailable_degenerate_whitening",
                "n_heldout": int(len(test_rows)), "n_training": int(len(train_rows))}
    observed = paired_absolute_correlation(px, py)

    held_cancers = cancers[test_rows]
    held_sites = sites[test_rows]
    rng = np.random.default_rng(seed)

    groups = [np.flatnonzero(held_cancers == level) for level in np.unique(held_cancers)]
    null = np.empty(n_permutations, dtype=np.float64)
    for i in range(n_permutations):
        order = np.arange(len(py))
        for rows in groups:
            order[rows] = rng.permutation(rows)
        null[i] = paired_absolute_correlation(px, py[order])
    null_p95 = float(np.percentile(null, 95))
    permutation_p = float((int(np.sum(null >= observed)) + 1.0) / (n_permutations + 1.0))

    unique_sites = np.unique(held_sites)
    site_rows = {s: np.flatnonzero(held_sites == s) for s in unique_sites}
    site_boot = np.full(n_boot, np.nan)
    for i in range(n_boot):
        drawn = rng.choice(len(unique_sites), size=len(unique_sites), replace=True)
        rows = np.concatenate([site_rows[unique_sites[j]] for j in drawn])
        site_boot[i] = paired_absolute_correlation(px[rows], py[rows])
    patient_boot = np.full(n_boot, np.nan)
    for i in range(n_boot):
        rows = rng.integers(0, len(py), size=len(py))
        patient_boot[i] = paired_absolute_correlation(px[rows], py[rows])

    def _ci(values):
        finite = values[np.isfinite(values)]
        if finite.size < 2:
            return [float("nan"), float("nan")]
        return [float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5))]

    site_ci = _ci(site_boot)
    return {
        "status": "scored",
        "adjusted": bool(adjust),
        "n_heldout": int(len(test_rows)), "n_training": int(len(train_rows)),
        "n_heldout_sites": int(len(unique_sites)),
        "observed": observed,
        "null_median": float(np.median(null)), "null_p95": null_p95,
        "permutation_p": permutation_p,
        "permutation_resolution": 1.0 / (n_permutations + 1.0),
        "site_cluster_ci95": site_ci,
        "patient_ci95": _ci(patient_boot),
        # The predeclared per-fold bar: clears the null AND the site-cluster interval
        # keeps it clear. The interval is the primary one because a site bootstrap is
        # the only resample that reflects "another draw of sites".
        "clears_null": bool(np.isfinite(observed) and observed > null_p95),
        "ci_clears_null": bool(np.isfinite(site_ci[0]) and site_ci[0] > null_p95),
    }


def leave_sites_out_channel(x: np.ndarray, y: np.ndarray, patient_ids: np.ndarray,
                            cancers: np.ndarray, *, n_folds: int = 5, n_components: int = 16,
                            adjust: bool = False, min_site_count: int = 10,
                            n_permutations: int = 1000, n_boot: int = 1000,
                            seed: int = 42) -> dict:
    """The whole predeclared test for one representation state and one arm."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    cancers = np.asarray(cancers).astype(str)
    # Raw TSS defines the folds: "OTHER" is not a site, it is a bag of rare sites,
    # and holding it out as a unit would hold out nothing in particular. Pooled TSS
    # is what the *adjustment* uses, exactly as every other CALIBRA analysis does.
    sites, _ = pooled_tissue_source_site(np.asarray(patient_ids).astype(str), min_site_count=1)
    pooled, frequent = pooled_tissue_source_site(np.asarray(patient_ids).astype(str),
                                                 min_site_count=min_site_count)
    fold = site_folds(sites, cancers, n_folds=n_folds)
    matched = matched_random_folds(cancers, fold, seed=seed)

    folds, comparator = [], []
    for f in range(n_folds):
        held = fold == f
        if not held.any():
            folds.append({"status": "unavailable_empty_fold", "fold": f})
            comparator.append({"status": "unavailable_empty_fold", "fold": f})
            continue
        result = evaluate_fold(x, y, cancers, sites, pooled, held, n_components=n_components,
                               adjust=adjust, n_permutations=n_permutations, n_boot=n_boot,
                               seed=seed + f)
        result["fold"] = f
        folds.append(result)
        reference = evaluate_fold(x, y, cancers, sites, pooled, matched == f,
                                  n_components=n_components, adjust=adjust,
                                  n_permutations=n_permutations, n_boot=n_boot, seed=seed + f)
        reference["fold"] = f
        comparator.append(reference)

    scored = [f for f in folds if f.get("status") == "scored"]
    matched_scored = [f for f in comparator if f.get("status") == "scored"]
    observed = np.asarray([f["observed"] for f in scored], dtype=np.float64)
    matched_observed = np.asarray([f["observed"] for f in matched_scored], dtype=np.float64)
    n_survive = sum(1 for f in scored if f["clears_null"] and f["ci_clears_null"])
    median_site = float(np.nanmedian(observed)) if observed.size else float("nan")
    median_matched = float(np.nanmedian(matched_observed)) if matched_observed.size else float("nan")
    ratio = float(median_site / median_matched) if np.isfinite(median_matched) and median_matched > 1e-12 else float("nan")

    # The verdict is the one fixed in the predeclaration; it is computed here rather
    # than read off by eye so that "it survived" cannot drift after the fact.
    majority = (len(scored) // 2) + 1
    if not scored:
        verdict = "unavailable_no_scored_folds"
    elif n_survive >= majority:
        verdict = "attenuated_but_present" if np.isfinite(ratio) and ratio < 0.5 else "survives"
    else:
        verdict = "collapses"

    return {
        "status": "scored" if scored else "unavailable_no_scored_folds",
        "test_name": "held-out-site",
        "not_this_test": "held-out-cancer",
        "discharges_no_external_cohort": False,
        "adjusted": bool(adjust),
        "adjustment": "cancer+pooled_tss_cross_fitted_within_block" if adjust else "none_raw_state",
        "n_folds": int(n_folds), "n_components": int(n_components),
        "n_patients": int(len(x)), "n_raw_sites": int(len(set(sites.tolist()))),
        "n_pooled_site_classes": int(len(set(pooled.tolist()))),
        "n_frequent_sites": int(len(frequent)), "min_site_count": int(min_site_count),
        "n_permutations": int(n_permutations), "n_boot": int(n_boot), "seed": int(seed),
        "composition": fold_composition(fold, sites, cancers),
        "folds": folds,
        "matched_random_folds": comparator,
        "n_folds_scored": len(scored),
        "n_folds_surviving": int(n_survive),
        "median_heldout_site_channel": median_site,
        "median_matched_random_channel": median_matched,
        "site_to_random_ratio": ratio,
        "verdict": verdict,
    }


def _load_targets(path: str, target_group: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The RNA target block, selected exactly as ``run_calibra.main`` selects it.

    ``--target-group`` empty is the default and is what reproduces the published
    channel: the block is every target whose name does NOT carry the
    ``RANDOM_CONTROL__`` prefix. Note that ``all_non_control`` is the *label*
    ``calibra_protocol.json`` records for that selection, not a value present in
    ``target_groups`` — passing it as a group would silently select nothing, so an
    empty selection raises rather than proceeding to a confusing failure downstream.
    """
    raw = np.load(path, allow_pickle=True)
    ids = np.asarray([str(p) for p in raw["patient_ids"]])
    names = np.asarray([str(t) for t in raw["target_names"]])
    groups = np.asarray([str(g) for g in raw["target_groups"]])
    scores = np.asarray(raw["scores"], dtype=np.float64)
    keep = np.ones(len(names), dtype=bool)
    if target_group:
        keep = groups == target_group
        if not keep.any():
            raise ValueError(f"--target-group {target_group!r} matches no target; "
                             f"available groups are {sorted(set(groups.tolist()))}")
    keep &= ~np.char.startswith(names, "RANDOM_CONTROL__")
    if not keep.any():
        raise ValueError("target selection is empty after excluding RANDOM_CONTROL__ columns")
    return ids, names[keep], scores[:, keep]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", nargs="+", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--partition", default="test", choices=("test", "val", "all"))
    parser.add_argument("--target-group", default="",
                        help="restrict to one target_group; empty (the default) is the "
                             "published channel's block: everything but RANDOM_CONTROL__")
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--min-site-count", type=int, default=10)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--states", nargs="*", default=[])
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    target_ids, target_names, scores = _load_targets(args.targets, args.target_group)
    index = {pid: i for i, pid in enumerate(target_ids)}

    results: dict[str, dict] = {}
    rows: list[dict] = []
    for artifact_path in args.artifacts:
        path = Path(artifact_path)
        method = path.stem
        raw = np.load(path, allow_pickle=True)
        declared = {str(s) for s in raw["trained_states"]} if "trained_states" in raw.files else set()
        patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
        cancers = np.asarray([str(c) for c in raw["cancers"]])
        split = np.asarray([str(s) for s in raw["split"]])
        mask = np.ones(len(patient_ids), dtype=bool) if args.partition == "all" else (split == args.partition)
        aligned = np.asarray([index.get(pid, -1) for pid in patient_ids])
        mask &= aligned >= 0
        y = scores[aligned[mask]]
        states = sorted(declared if not args.states else (declared & set(args.states)))
        for state in states:
            if state not in raw.files:
                results[f"{method}::{state}"] = {"status": "unavailable_state_declared_but_absent"}
                continue
            x = np.asarray(raw[state], dtype=np.float64)[mask]
            for adjust in (False, True):
                arm = "adjusted" if adjust else "unadjusted"
                key = f"{method}::{state}::{arm}"
                result = leave_sites_out_channel(
                    x, y, patient_ids[mask], cancers[mask], n_folds=args.n_folds,
                    n_components=args.n_components, adjust=adjust,
                    min_site_count=args.min_site_count, n_permutations=args.n_permutations,
                    n_boot=args.n_boot, seed=args.seed)
                results[key] = result
                for metric in ("median_heldout_site_channel", "median_matched_random_channel",
                               "site_to_random_ratio", "n_folds_surviving", "n_folds_scored"):
                    rows.append({"method": method, "representation_state": state,
                                 "task": "leave_sites_out", "target": arm, "metric": metric,
                                 "value": float(result.get(metric, float("nan"))),
                                 "note": result.get("verdict", result["status"])})
                print(f"[{key}] verdict={result.get('verdict')} "
                      f"site={result.get('median_heldout_site_channel', float('nan')):.4f} "
                      f"random={result.get('median_matched_random_channel', float('nan')):.4f} "
                      f"ratio={result.get('site_to_random_ratio', float('nan')):.3f} "
                      f"survive={result.get('n_folds_surviving')}/{result.get('n_folds_scored')}",
                      flush=True)
    (output / "leave_sites_out.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    import pandas as pd
    pd.DataFrame(rows).to_csv(output / "task_rows.csv", index=False)
    print(f"[done] {len(rows)} rows -> {output}", flush=True)


if __name__ == "__main__":
    main()
