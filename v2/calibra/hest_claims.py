"""Spot-level replication of five project claims on the HEST-1k spatial cohort.

Every reading this harness produces was written down first, in
``NOTEBOOK_ENTRIES/PREDECLARED_spatial_claim_replication_20260804T1800Z.md``.  The harness
exists so that the claims are re-measured by the *same* imported instrument that produced
the bulk-RNA numbers, with only the cohort and the confound changed:

  claim 1  confound adjustment works two-sided ....... ``confound_certificate.certify_axes``
                                                       + ``calibration.spike_recovery_curve``
  claim 2  chance is not zero ......................... ``calibration.permutation_null``
  claim 3  random gene sets reproduce the signal ...... ``run_calibra.score_target_block_per_column``
  claim 4  per-axis certification is insufficient ..... ``confound_certificate.certify_axes``
  claim 5  the zero-parameter baseline ................ ``hest.{pooled_r, within_slide_r}``

Three things in here are new, and each exists for a stated reason.

1. SUBSAMPLING IS DECLARED AND SLIDE-BALANCED.  ``permutation_null`` re-whitens the image
   block *inside* its permutation loop, so its cost is O(n * d^2) per permutation and a
   53,217-spot run is a multi-hour job (measured, 2026-08-03).  ``slide_stratified_subsample``
   takes a seeded ``min(m, n_slide)`` draw from every slide, so slide -- the confound under
   test -- stays balanced and the sample size is a declared parameter rather than a silent
   truncation.

2. A NON-MEAN-BASED CONFOUND PROBE.  ``certify_axes`` scores nearest-class-mean per axis and
   shrunk-covariance LDA jointly.  Both are mean-based, and cross-fitted residualisation on a
   one-hot slide design removes slide *means* by construction, so "the adjusted state cannot
   predict slide" is close to arithmetic.  ``knn_balanced_accuracy_oof`` asks the same
   question with a classifier that reads local neighbourhood structure instead of class
   means.  If it still recovers slide after adjustment, the certificate's pass is a statement
   about linear-mean classifiers, not about the representation.  This applies verbatim to the
   TCGA site certificate the claim comes from.

3. GENUINE RANDOM GENE SETS.  The 16 ``RANDOM_CONTROL__`` columns shipped in the target file
   are row-permutations of real columns: a *pairing* null, not a random gene set, and it
   should read ~0.  TCGA's T1.4 control was a covariate-matched random gene set, which keeps
   real cross-gene covariance.  ``matched_random_gene_sets`` reproduces that object here by
   drawing non-panel genes matched to each panel gene on train-spot (mean, log-variance).
   Both arms are reported; conflating them would answer a different question than the claim.

Normalisation note, load-bearing for claim 3.  ``hest_build.stage_assemble`` builds the real
targets as ``normalise_expression(panel_counts)`` -- the library size is the sum over the 50
PANEL genes, not over the transcriptome.  A random gene set is therefore normalised by its
own 50-gene total, so the two blocks differ in gene identity and in nothing else.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from .calibration import permutation_null, spike_recovery_curve
from .confound_certificate import _encode, _stratified_folds, balanced_accuracy, certify_axes
from .hest import (HestAdapterError, normalise_expression, per_slide_mean_baseline, pooled_r,
                   within_slide_r)
from .residualise import confound_design, cross_fitted_residuals, pooled_tissue_source_site
from .run_calibra import random_direction_column_correlation, score_target_block_per_column
from .spectral import cca_spectrum, effective_rank, heldout_top_cca

__all__ = ["slide_stratified_subsample", "knn_balanced_accuracy_oof", "matched_random_gene_sets",
           "spatial_confound_design", "capacity_floor_prediction", "slide_bootstrap_interval",
           "between_slide_variance_share"]

RANDOM_CONTROL_PREFIX = "RANDOM_CONTROL__"


# --- declared subsampling --------------------------------------------------------------

def slide_stratified_subsample(slide_ids, base_mask, per_slide: int, *, seed: int = 42) -> np.ndarray:
    """Row indices of a seeded, slide-balanced draw of ``min(per_slide, n_slide)`` spots.

    Slide-balanced rather than uniform because slide is the confound under test: a uniform
    draw would give the 9,079-spot slide 5.5x the weight of the 1,653-spot one and any
    slide-adjusted statistic would then be dominated by a single tissue block.  ``per_slide``
    <= 0 means "every row in ``base_mask``", so the full-cohort arm goes through the identical
    code path as the subsampled arms.
    """
    slides = np.asarray(slide_ids).astype(str)
    base = np.asarray(base_mask, dtype=bool)
    if len(slides) != len(base):
        raise HestAdapterError("slide ids and mask must be aligned")
    rows = np.flatnonzero(base)
    if per_slide is None or int(per_slide) <= 0:
        return rows
    rng = np.random.default_rng(seed)
    picked = []
    for slide in np.unique(slides[rows]):
        candidates = rows[slides[rows] == slide]
        take = int(min(int(per_slide), len(candidates)))
        picked.append(rng.choice(candidates, size=take, replace=False))
    return np.sort(np.concatenate(picked)) if picked else np.zeros(0, dtype=np.int64)


def spatial_confound_design(patient_ids, cancers, *, min_site_count: int = 10):
    """CALIBRA's own design, with the site term resolved on spot keys.

    Returned alongside the derived site labels so a caller can assert that the site field
    really is the slide.  ``pooled_tissue_source_site`` splits on ``-`` and takes field 1;
    ``hest.spot_key`` puts the slide there precisely so that this function needs no spatial
    special case.
    """
    import pandas as pd

    site, frequent = pooled_tissue_source_site(np.asarray(patient_ids).astype(str),
                                               min_site_count=min_site_count)
    frame = pd.DataFrame({"cancer": np.asarray(cancers).astype(str), "tss": site})
    return confound_design(frame, ["cancer", "tss"]), site, frequent


# --- claim 1c: a confound probe that is not mean-based ----------------------------------

def knn_balanced_accuracy_oof(features: np.ndarray, class_index: np.ndarray, n_classes: int, *,
                              k: int = 15, n_splits: int = 5, seed: int = 42,
                              folds: list[np.ndarray] | None = None, block: int = 512) -> float:
    """Out-of-fold k-nearest-neighbour balanced accuracy for a categorical confound.

    The companion to ``confound_certificate.certify_axes`` and the reason it is here: that
    certificate's two classifiers (nearest-class-mean per axis, shrunk LDA jointly) are both
    functions of class MEANS, and the adjustment being certified is a linear regression that
    removes class means.  Passing it is therefore nearly tautological.  A kNN vote depends on
    local neighbourhood structure, which mean-removal does not touch, so it can distinguish
    "the confound is gone" from "the confound is no longer in the first moment".

    Plain majority vote, no prior correction, because the intended caller passes a
    slide-balanced subsample in which the class priors are uniform by construction; the
    balance is asserted by the caller, not assumed here.
    """
    features = np.asarray(features, dtype=np.float64)
    class_index = np.asarray(class_index, dtype=np.int64)
    if features.ndim != 2 or len(features) != len(class_index):
        raise HestAdapterError("kNN probe needs an aligned [row, feature] matrix and label vector")
    folds = folds if folds is not None else _stratified_folds(class_index, n_splits, seed)
    predicted = np.empty(len(class_index), dtype=np.int64)
    square = (features ** 2).sum(axis=1)
    for test_rows in folds:
        train_rows = np.setdiff1d(np.arange(len(class_index)), test_rows, assume_unique=False)
        if len(train_rows) <= k:
            raise HestAdapterError("not enough training rows for the requested k")
        train_features = features[train_rows]
        train_square = square[train_rows]
        train_labels = class_index[train_rows]
        for start in range(0, len(test_rows), block):
            rows = test_rows[start:start + block]
            # squared euclidean distance without materialising differences
            distance = (square[rows][:, None] + train_square[None, :]
                        - 2.0 * (features[rows] @ train_features.T))
            nearest = np.argpartition(distance, kth=k - 1, axis=1)[:, :k]
            votes = train_labels[nearest]
            counts = np.zeros((len(rows), n_classes), dtype=np.int64)
            np.add.at(counts, (np.repeat(np.arange(len(rows)), k), votes.ravel()), 1)
            predicted[rows] = np.argmax(counts, axis=1)
    return balanced_accuracy(class_index, predicted, n_classes)


# --- claim 2: the predeclared capacity-floor law ----------------------------------------

def capacity_floor_prediction(n: int, n_components: int, design_rank: int = 0) -> float:
    """2*sqrt(k/n_eff): the Wachter/Marchenko-Pastur edge for two independent whitened blocks.

    This is the formula claim 2 was predeclared against, kept next to the measurement so the
    two can never drift apart.  ``design_rank`` is subtracted from n because both blocks are
    cross-fitted residuals of the same confound design and the residual space is that much
    smaller.  It is NOT a measured statistic and is never reported as one.
    """
    effective = max(int(n) - int(design_rank), 1)
    return float(2.0 * np.sqrt(float(n_components) / float(effective)))


# --- claim 5: intervals whose resampling unit is the slide ------------------------------

def slide_bootstrap_interval(per_slide_values: dict, *, n_boot: int = 2000, seed: int = 42) -> dict:
    """Percentile CI for a per-slide statistic, resampling SLIDES rather than spots.

    Spots inside a slide are neither independent nor even non-overlapping (a 128 um window on
    a 100 um pitch overlaps its neighbours), so a spot bootstrap would report an interval
    several times too narrow.  The bootstrap universe is 13 test slides, which is small, and
    the resulting interval is honestly wide.
    """
    slides = sorted(per_slide_values)
    values = np.asarray([np.nanmean(per_slide_values[s]) for s in slides], dtype=np.float64)
    if values.size < 2:
        return {"point": float(np.nanmean(values)) if values.size else float("nan"),
                "ci95": [float("nan"), float("nan")], "n_slides": int(values.size)}
    rng = np.random.default_rng(seed)
    draws = values[rng.integers(0, len(values), size=(int(n_boot), len(values)))].mean(axis=1)
    return {"point": float(np.nanmean(values)), "n_slides": int(values.size),
            "ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))]}


def between_slide_variance_share(targets: np.ndarray, slide_ids, mask=None) -> np.ndarray:
    """Per-gene share of total variance that lies BETWEEN slides rather than within them.

    The plain statement claim 5 asks for: if this is large, a pooled spot-level correlation
    is mostly a between-slide correlation and a per-slide mean can win it with no parameters.
    """
    targets = np.asarray(targets, dtype=np.float64)
    slides = np.asarray(slide_ids).astype(str)
    rows = np.ones(len(targets), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    values, groups = targets[rows], slides[rows]
    grand = values.mean(axis=0, keepdims=True)
    total = ((values - grand) ** 2).sum(axis=0)
    between = np.zeros(values.shape[1])
    for slide in np.unique(groups):
        sel = groups == slide
        between += sel.sum() * (values[sel].mean(axis=0) - grand[0]) ** 2
    return np.where(total > 1e-12, between / np.maximum(total, 1e-12), np.nan)


# --- claim 3: genuine random gene sets --------------------------------------------------

def matched_random_gene_sets(train_mean: np.ndarray, train_logvar: np.ndarray,
                             gene_names, panel_genes, *, n_sets: int = 5,
                             seed: int = 101) -> list[np.ndarray]:
    """Non-panel gene sets matched to the panel on (mean, log-variance), one draw per set.

    TCGA's T1.4 control matched training-only per-gene mean, variance and PC1 loading.  The
    first two are reproduced here; PC1 loading is not, and that is a declared weakening of
    the control in the direction that makes the control EASIER to beat, i.e. it biases
    against the claim rather than for it.

    Matching is a nearest-neighbour draw in standardised (mean, log-variance) space among a
    candidate pool of the 32 closest non-panel genes to each panel gene, sampled without
    replacement across the set so no gene appears twice.
    """
    names = np.asarray(gene_names).astype(str)
    panel = np.asarray(panel_genes).astype(str)
    mean = np.asarray(train_mean, dtype=np.float64)
    logvar = np.asarray(train_logvar, dtype=np.float64)
    if not (len(names) == len(mean) == len(logvar)):
        raise HestAdapterError("gene names, means and log-variances must be aligned")
    index = {g: i for i, g in enumerate(names.tolist())}
    missing = [g for g in panel.tolist() if g not in index]
    if missing:
        raise HestAdapterError(f"panel genes absent from the statistics table: {missing[:5]}")
    panel_rows = np.asarray([index[g] for g in panel.tolist()])
    is_panel = np.zeros(len(names), dtype=bool)
    is_panel[panel_rows] = True
    candidates = np.flatnonzero(~is_panel)
    if len(candidates) < 32 * len(panel_rows) // 8:
        raise HestAdapterError("candidate pool is too small to match the panel")
    scale_m = max(float(np.std(mean)), 1e-12)
    scale_v = max(float(np.std(logvar)), 1e-12)
    coords = np.column_stack([mean / scale_m, logvar / scale_v])
    sets = []
    for s in range(int(n_sets)):
        rng = np.random.default_rng(int(seed) + s)
        taken: set[int] = set()
        chosen = []
        for row in panel_rows:
            distance = ((coords[candidates] - coords[row]) ** 2).sum(axis=1)
            order = candidates[np.argsort(distance, kind="stable")]
            pool = [int(c) for c in order[:64] if int(c) not in taken][:32]
            if not pool:
                raise HestAdapterError("exhausted the matched candidate pool")
            pick = int(rng.choice(pool))
            taken.add(pick)
            chosen.append(pick)
        sets.append(np.asarray(sorted(chosen), dtype=np.int64))
    return sets


# --- stages -----------------------------------------------------------------------------

def _load(artifact: str, targets: str, partition: str):
    art = np.load(artifact, allow_pickle=True)
    tgt = np.load(targets, allow_pickle=True)
    ids = art["patient_ids"].astype(str)
    slides = art["slide_ids"].astype(str)
    cancers = art["cancers"].astype(str)
    split = art["split"].astype(str)
    x = art["wsi_identity"]
    names = tgt["target_names"].astype(str)
    order = {p: i for i, p in enumerate(tgt["patient_ids"].astype(str).tolist())}
    y_all = tgt["scores"][[order[p] for p in ids.tolist()]]
    is_control = np.char.startswith(names, RANDOM_CONTROL_PREFIX)
    mask = np.ones(len(ids), dtype=bool) if partition == "all" else (split == partition)
    # The fix this cohort depends on: the TSS field of a spot key must be the slide.
    site_field = np.asarray([p.split("-")[1] if len(p.split("-")) > 1 else "NA" for p in ids])
    if not (site_field == slides).all():
        raise HestAdapterError("spot keys do not put the slide in the TSS field; the slide "
                               "confound would silently drop out of every design")
    return {"ids": ids, "slides": slides, "cancers": cancers, "split": split, "x": x,
            "y": y_all[:, ~is_control], "y_names": names[~is_control],
            "controls": y_all[:, is_control], "control_names": names[is_control], "mask": mask}


def _rows_for(data, per_slide: int, seed: int) -> np.ndarray:
    return slide_stratified_subsample(data["slides"], data["mask"], per_slide, seed=seed)


def _blocks(data, rows):
    x = np.asarray(data["x"], dtype=np.float64)[rows]
    y = np.asarray(data["y"], dtype=np.float64)[rows]
    design, site, frequent = spatial_confound_design(data["ids"][rows], data["cancers"][rows])
    return x, y, design, site, frequent


def stage_null_sweep(args) -> None:
    """Claim 2: the spot-scale permutation null, against the predeclared 2*sqrt(k/n)."""
    data = _load(args.artifact, args.targets, args.partition)
    grid = [int(v) for v in args.per_slide_grid.split(",")]
    perms = [int(v) for v in args.perm_grid.split(",")]
    if len(perms) != len(grid):
        raise HestAdapterError("perm-grid must have one entry per per-slide-grid entry")
    results = []
    for per_slide, n_permutations in zip(grid, perms):
        rows = _rows_for(data, per_slide, args.seed)
        x, y, design, site, _ = _blocks(data, rows)
        strata_arms = {"cancer": data["cancers"][rows], "slide": site}
        record = {"per_slide": per_slide, "n": int(len(rows)),
                  "n_slides": int(np.unique(data["slides"][rows]).size),
                  "n_design_columns": int(design.shape[1]),
                  "design_rank": int(np.linalg.matrix_rank(design)),
                  "n_permutations": n_permutations, "n_components": args.n_components,
                  "predicted_2sqrt_k_over_n": capacity_floor_prediction(
                      len(rows), args.n_components,
                      design_rank=int(np.linalg.matrix_rank(design)))}
        for arm, strata in strata_arms.items():
            start = time.time()
            null = permutation_null(x, y, design, strata=strata, n_permutations=n_permutations,
                                    n_components=args.n_components, seed=args.seed,
                                    n_jobs=args.n_jobs)
            elapsed = time.time() - start
            record[f"strata_{arm}"] = {**{k: (float(v) if isinstance(v, float) else v)
                                          for k, v in null.items()},
                                       "wall_seconds": round(elapsed, 2),
                                       "seconds_per_permutation": round(elapsed / max(n_permutations, 1), 3)}
            print(f"[null] n={len(rows):6d} strata={arm:6s} observed={null['observed_top_cca']:.4f} "
                  f"null_median={null['null_median']:.4f} p95={null['null_p95']:.4f} "
                  f"predicted={record['predicted_2sqrt_k_over_n']:.4f} "
                  f"({elapsed/max(n_permutations,1):.2f} s/perm)", flush=True)
        results.append(record)
        _write(args.output, "claim2_null_sweep.json", {"grid": results, "seed": args.seed,
                                                       "partition": args.partition})


def stage_certificate(args) -> None:
    """Claims 1a, 1c and 4: raw vs adjusted slide predictability, three classifiers."""
    data = _load(args.artifact, args.targets, args.partition)
    rows = _rows_for(data, args.per_slide, args.seed)
    x = np.asarray(data["x"], dtype=np.float64)[rows]
    ids, cancers, slides = data["ids"][rows], data["cancers"][rows], data["slides"][rows]
    counts = np.unique(slides, return_counts=True)[1]
    balanced = bool(counts.max() == counts.min())
    out = {"n": int(len(rows)), "n_slides": int(len(counts)), "slide_balanced": balanced,
           "per_slide": args.per_slide, "n_permutations": args.n_permutations,
           "n_boot": args.n_boot, "seed": args.seed}
    classes, class_index = _encode(pooled_tissue_source_site(ids)[0])
    n_classes = int(len(classes))
    out["n_classes"] = n_classes
    out["chance_rate"] = 1.0 / n_classes
    design, _, _ = spatial_confound_design(ids, cancers)
    adjusted = cross_fitted_residuals(x, design, seed=args.seed)
    for arm, features in (("raw", x), ("adjusted", adjusted)):
        start = time.time()
        out[f"knn_{arm}"] = knn_balanced_accuracy_oof(features, class_index, n_classes,
                                                      k=args.knn_k, seed=args.seed)
        print(f"[knn ] {arm:8s} balanced_accuracy={out[f'knn_{arm}']:.4f} "
              f"chance={1.0/n_classes:.4f} ({time.time()-start:.0f} s)", flush=True)
    for residualise in (False, True):
        arm = "adjusted" if residualise else "raw"
        start = time.time()
        result = certify_axes(x, ids, cancers, n_permutations=args.n_permutations,
                              n_boot=args.n_boot, n_boot_axes=args.n_boot_axes,
                              seed=args.seed, residualise=residualise, n_jobs=args.n_jobs)
        result["wall_seconds"] = round(time.time() - start, 1)
        # the whole per-axis vector is large and is the evidence for claim 4; keep the tail
        for key in ("per_axis_balanced_accuracy", "per_axis_null_p95", "per_axis_permutation_p"):
            values = np.asarray(result.pop(key), dtype=np.float64)
            result[f"{key}_top20"] = np.sort(values)[::-1][:20].tolist()
            result[f"{key}_median"] = float(np.median(values))
        out[f"certificate_{arm}"] = result
        print(f"[cert] {arm:8s} per_axis_max={result['observed_max']:.4f} "
              f"joint={result['joint_lda_balanced_accuracy']:.4f} chance={result['chance_rate']:.4f} "
              f"breach={result['n_breaching_axes']}/{result['n_axes']} "
              f"({result['wall_seconds']:.0f} s)", flush=True)
    _write(args.output, "claim1_claim4_certificate.json", out)


def stage_channel(args) -> None:
    """Claim 1b: what the slide adjustment costs the channel, and the induced level-0 floor."""
    data = _load(args.artifact, args.targets, args.partition)
    rows = _rows_for(data, args.per_slide, args.seed)
    x, y, design, site, _ = _blocks(data, rows)
    result = spike_recovery_curve(x, y, design, levels=tuple(float(v) for v in args.levels.split(",")),
                                  n_draws=args.n_draws, n_components=args.n_components,
                                  seed=args.seed, n_jobs=args.n_jobs)
    summary = result.summary()
    unadjusted = cca_spectrum(x, y, n_components=args.n_components)
    summary.update({
        "n": int(len(rows)), "per_slide": args.per_slide,
        "n_design_columns": int(design.shape[1]),
        "design_rank": int(np.linalg.matrix_rank(design)),
        "effective_rank_x": effective_rank(x),
        "unadjusted_top_cca": float(unadjusted[0]) if unadjusted.size else float("nan"),
        "heldout_top_cca": heldout_top_cca(cross_fitted_residuals(x, design, seed=args.seed),
                                           cross_fitted_residuals(y, design, seed=args.seed),
                                           n_components=args.n_components, seed=args.seed),
    })
    print(f"[chan] n={summary['n']} unadjusted={summary['unadjusted_top_cca']:.4f} "
          f"adjusted={summary['observed']:.4f} heldout={summary['heldout_top_cca']:.4f} "
          f"attenuation={summary['attenuation_slope']:.4f} "
          f"induced_baseline={summary['baseline_recovered_median']:.4f} "
          f"transmission_floor={summary['transmission_floor']} "
          f"detection_floor={summary['detection_floor']}", flush=True)
    _write(args.output, "claim1b_channel.json", summary)


def stage_gene_sets(args) -> None:
    """Claim 3: permuted-column arm and genuine matched-random-gene-set arm, same rows."""
    data = _load(args.artifact, args.targets, args.partition)
    rows = _rows_for(data, args.per_slide, args.seed)
    x, y, design, _, _ = _blocks(data, rows)
    x_res = cross_fitted_residuals(x, design, seed=args.seed)
    y_res = cross_fitted_residuals(y, design, seed=args.seed)
    real_fitted = score_target_block_per_column(x_res, y_res, seed=args.seed)
    real_matched = random_direction_column_correlation(x_res, y_res, n_draws=args.n_draws,
                                                       seed=args.seed)
    out = {"n": int(len(rows)), "per_slide": args.per_slide, "seed": args.seed,
           "real_fitted_median": float(np.nanmedian(real_fitted)),
           "real_fitted_p95": float(np.nanpercentile(real_fitted, 95)),
           "real_matched_median": float(np.nanmedian(real_matched)),
           "n_real_targets": int(y.shape[1])}

    control = np.asarray(data["controls"], dtype=np.float64)[rows]
    control_res = cross_fitted_residuals(control, design, seed=args.seed)
    permuted_fitted = score_target_block_per_column(x_res, control_res, seed=args.seed)
    out["permuted_arm"] = {
        "n_columns": int(control.shape[1]),
        "fitted_median": float(np.nanmedian(permuted_fitted)),
        "fitted_p95": float(np.nanpercentile(permuted_fitted, 95)),
        "ratio_to_real": float(np.nanmedian(permuted_fitted) / np.nanmedian(real_fitted))
        if abs(float(np.nanmedian(real_fitted))) > 1e-9 else float("nan"),
        "what_this_is": "row-permutation of real columns: a pairing null, NOT a random gene set"}
    print(f"[gset] real_fitted_median={out['real_fitted_median']:.4f} "
          f"permuted={out['permuted_arm']['fitted_median']:.4f}", flush=True)

    if args.random_sets:
        blob = np.load(args.random_sets, allow_pickle=True)
        block_ids = blob["spot_ids"].astype(str)
        position = {p: i for i, p in enumerate(block_ids.tolist())}
        take = np.asarray([position[p] for p in data["ids"][rows].tolist()])
        arms = []
        for name in [str(k) for k in blob["set_names"]]:
            block = np.asarray(blob[name], dtype=np.float64)[take]
            block_res = cross_fitted_residuals(block, design, seed=args.seed)
            fitted = score_target_block_per_column(x_res, block_res, seed=args.seed)
            arms.append({"set": name, "n_columns": int(block.shape[1]),
                         "fitted_median": float(np.nanmedian(fitted)),
                         "fitted_p95": float(np.nanpercentile(fitted, 95)),
                         "ratio_to_real": float(np.nanmedian(fitted) / np.nanmedian(real_fitted))})
            print(f"[gset] {name}: fitted_median={arms[-1]['fitted_median']:.4f} "
                  f"ratio={arms[-1]['ratio_to_real']:.4f}", flush=True)
        ratios = np.asarray([a["ratio_to_real"] for a in arms])
        out["random_gene_set_arm"] = {
            "sets": arms, "n_sets": len(arms),
            "ratio_median": float(np.median(ratios)), "ratio_min": float(ratios.min()),
            "ratio_max": float(ratios.max()),
            "matched_on": str(blob["matched_on"]) if "matched_on" in blob.files else "",
            "what_this_is": "non-panel genes matched to the panel on train (mean, log-variance)"}
    _write(args.output, "claim3_gene_sets.json", out)


def stage_random_sets(args) -> None:
    """Build the matched random-gene-set target blocks from the raw HEST h5ad store."""
    plan = json.load(open(args.plan))
    panel = np.asarray(plan["panel_genes"], dtype=str)
    root = Path(args.root)
    slides = plan["slides"]
    train_ids = [s["id"] for s in slides if s["split"] == "train"]
    eval_ids = [s["id"] for s in slides if s["split"] == args.partition]

    common = None
    for meta in slides:
        genes = _h5ad_genes(root / "st" / f"{meta['id']}.h5ad")
        common = genes if common is None else np.intersect1d(common, genes)
    common = np.asarray(sorted(common.tolist()), dtype=str)
    print(f"[rsets] gene intersection = {len(common)}", flush=True)

    # Streaming train statistics on the SAME normalisation gene selection used, so the match
    # is made in the space the panel was chosen in.
    total = np.zeros(len(common))
    total_sq = np.zeros(len(common))
    n_rows = 0
    for sid in train_ids:
        expr = _h5ad_expression(root / "st" / f"{sid}.h5ad", common, plan)
        total += expr.sum(axis=0)
        total_sq += (expr ** 2).sum(axis=0)
        n_rows += len(expr)
        print(f"[rsets] train {sid}: {len(expr)} spots", flush=True)
    mean = total / max(n_rows, 1)
    var = np.maximum(total_sq / max(n_rows, 1) - mean ** 2, 1e-12)
    sets = matched_random_gene_sets(mean, np.log(var), common, panel,
                                    n_sets=args.n_sets, seed=args.seed)

    blocks: dict[str, list[np.ndarray]] = {f"random_set_{i}": [] for i in range(len(sets))}
    ids: list[np.ndarray] = []
    for sid in eval_ids:
        raw, barcodes = _h5ad_counts(root / "st" / f"{sid}.h5ad", common, plan)
        ids.append(np.asarray([f"HEST-{sid}-{b}" for b in barcodes], dtype=str))
        for i, gene_rows in enumerate(sets):
            # own-set library size, exactly as stage_assemble did for the 50-gene panel
            blocks[f"random_set_{i}"].append(normalise_expression(raw[:, gene_rows]))
        print(f"[rsets] {args.partition} {sid}: {len(barcodes)} spots", flush=True)
    arrays = {name: np.vstack(parts).astype(np.float32) for name, parts in blocks.items()}
    arrays["spot_ids"] = np.concatenate(ids)
    arrays["set_names"] = np.asarray(sorted(blocks), dtype=str)
    arrays["genes"] = np.asarray([",".join(common[s].tolist()) for s in sets], dtype=str)
    arrays["matched_on"] = np.asarray("train_spot_mean_and_log_variance_nearest_neighbour_pool32")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **arrays)
    print(f"[rsets] wrote {args.output}: {arrays['spot_ids'].shape[0]} spots, "
          f"{len(sets)} sets of {len(panel)} genes", flush=True)


def _h5ad_genes(path: Path) -> np.ndarray:
    import anndata

    return np.asarray(anndata.read_h5ad(path, backed="r").var_names).astype(str)


def _h5ad_counts(path: Path, common: np.ndarray, plan: dict):
    """Usable-spot counts restricted to ``common``, with the build's own filter applied.

    ``usable_spots`` is evaluated on the FULL transcriptome, exactly as ``hest_build``'s tile
    stage did, and only then are the genes restricted.  Filtering on the restricted matrix
    would apply a stricter library-size cut than the artifact used and would silently drop
    spots the artifact kept, so the two stores would no longer be row-comparable.
    """
    import anndata

    from .hest import usable_spots

    adata = anndata.read_h5ad(path)
    counts = adata.X
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    counts = np.asarray(counts, dtype=np.float32)
    keep = usable_spots(counts, plan["min_counts"], plan["min_genes"])
    order = {g: i for i, g in enumerate(np.asarray(adata.var_names).astype(str).tolist())}
    counts = counts[keep][:, [order[g] for g in common.tolist()]]
    barcodes = np.asarray(adata.obs_names).astype(str)[keep]
    return counts, barcodes


def _h5ad_expression(path: Path, common: np.ndarray, plan: dict) -> np.ndarray:
    counts, _ = _h5ad_counts(path, common, plan)
    return normalise_expression(counts).astype(np.float64)


def stage_baselines(args) -> None:
    """Claim 5: the zero-parameter baseline with slide-bootstrap intervals."""
    from sklearn.linear_model import Ridge

    data = _load(args.artifact, args.targets, "all")
    slides, split = data["slides"], data["split"]
    y = np.asarray(data["y"], dtype=np.float64)
    x = np.asarray(data["x"], dtype=np.float64)
    train, test = split == "train", split == args.partition
    predictions = {
        "per_slide_mean": per_slide_mean_baseline(y, slides, test),
        "global_mean": np.tile(y[train].mean(axis=0), (len(y), 1)),
        "ridge_hoptimus": Ridge(alpha=args.ridge_alpha).fit(x[train], y[train]).predict(x),
    }
    out = {"n_train_spots": int(train.sum()), "n_test_spots": int(test.sum()),
           "n_test_slides": int(np.unique(slides[test]).size), "n_genes": int(y.shape[1]),
           "ridge_alpha": args.ridge_alpha, "n_boot": args.n_boot,
           "bootstrap_unit": "slide", "partition": args.partition, "results": {}}
    share = between_slide_variance_share(y, slides, test)
    out["between_slide_variance_share"] = {
        "median": float(np.nanmedian(share)), "mean": float(np.nanmean(share)),
        "p10": float(np.nanpercentile(share, 10)), "p90": float(np.nanpercentile(share, 90))}
    for name, prediction in predictions.items():
        per_slide_within = {}
        for slide in np.unique(slides[test]):
            sel = test & (slides == slide)
            per_slide_within[slide] = within_slide_r(y, prediction, slides, sel)
        out["results"][name] = {
            "pooled_r_all_spots": float(np.nanmean(pooled_r(y, prediction, test))),
            "pooled_r_all_spots_median_gene": float(np.nanmedian(pooled_r(y, prediction, test))),
            "within_slide_r": slide_bootstrap_interval(per_slide_within, n_boot=args.n_boot,
                                                       seed=args.seed)}
        print(f"[base] {name:16s} pooled={out['results'][name]['pooled_r_all_spots']:.4f} "
              f"within={out['results'][name]['within_slide_r']['point']:.4f} "
              f"CI={out['results'][name]['within_slide_r']['ci95']}", flush=True)
    # The pooled statistic has no per-slide decomposition (it is defined across slides), so
    # its interval comes from a slide bootstrap of the pooled computation itself.
    rng = np.random.default_rng(args.seed)
    unique_slides = np.unique(slides[test])
    pooled_boot = {name: [] for name in predictions}
    for _ in range(args.n_boot // 4):
        draw = rng.choice(unique_slides, size=len(unique_slides), replace=True)
        sel = np.concatenate([np.flatnonzero(test & (slides == s)) for s in draw])
        for name, prediction in predictions.items():
            pooled_boot[name].append(float(np.nanmean(pooled_r(y[sel], prediction[sel]))))
    for name in predictions:
        values = np.asarray(pooled_boot[name])
        out["results"][name]["pooled_r_ci95"] = [float(np.percentile(values, 2.5)),
                                                 float(np.percentile(values, 97.5))]
        print(f"[base] {name:16s} pooled CI={out['results'][name]['pooled_r_ci95']}", flush=True)
    _write(args.output, "claim5_baselines.json", out)


def _write(output: str, name: str, payload: dict) -> None:
    directory = Path(output)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / name).write_text(json.dumps(payload, indent=1, default=float), encoding="utf-8")
    print(f"[done] wrote {directory / name}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)

    def common(p, *, partition_default="test"):
        p.add_argument("--artifact", required=True)
        p.add_argument("--targets", required=True)
        p.add_argument("--output", required=True)
        p.add_argument("--partition", default=partition_default, choices=("test", "val", "all"))
        p.add_argument("--seed", type=int, default=42)
        p.add_argument("--n-jobs", type=int, default=1)
        return p

    n = common(sub.add_parser("null-sweep"))
    n.add_argument("--per-slide-grid", default="40,100,213,400,800,1600,0")
    n.add_argument("--perm-grid", default="100,100,100,100,100,50,20")
    n.add_argument("--n-components", type=int, default=16)
    n.set_defaults(func=stage_null_sweep)

    c = common(sub.add_parser("certificate"))
    c.add_argument("--per-slide", type=int, default=213)
    c.add_argument("--n-permutations", type=int, default=200)
    c.add_argument("--n-boot", type=int, default=100)
    c.add_argument("--n-boot-axes", type=int, default=8)
    c.add_argument("--knn-k", type=int, default=15)
    c.set_defaults(func=stage_certificate)

    h = common(sub.add_parser("channel"))
    h.add_argument("--per-slide", type=int, default=213)
    h.add_argument("--levels", default="0.0,0.01,0.02,0.05,0.10,0.20,0.40")
    h.add_argument("--n-draws", type=int, default=25)
    h.add_argument("--n-components", type=int, default=16)
    h.set_defaults(func=stage_channel)

    g = common(sub.add_parser("gene-sets"))
    g.add_argument("--per-slide", type=int, default=213)
    g.add_argument("--n-draws", type=int, default=40)
    g.add_argument("--random-sets", default="", help="npz written by the random-sets stage")
    g.set_defaults(func=stage_gene_sets)

    r = sub.add_parser("random-sets")
    r.add_argument("--plan", required=True)
    r.add_argument("--root", required=True)
    r.add_argument("--output", required=True)
    r.add_argument("--partition", default="test")
    r.add_argument("--n-sets", type=int, default=5)
    r.add_argument("--seed", type=int, default=101)
    r.set_defaults(func=stage_random_sets)

    b = common(sub.add_parser("baselines"))
    b.add_argument("--ridge-alpha", type=float, default=1000.0)
    b.add_argument("--n-boot", type=int, default=2000)
    b.set_defaults(func=stage_baselines)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
