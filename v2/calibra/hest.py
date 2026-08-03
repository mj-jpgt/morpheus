"""HEST-1k spatial-transcriptomics adapter: a second target modality in the V2 artifact contract.

Every molecular target in this project so far is *bulk* RNA -- expression averaged over a
whole slide, while the image is sampled from tumour regions.  That mismatch is a standing
confound.  Spot-level spatial transcriptomics localises the target to the tissue it came
from, so the morphology->molecular channel can be measured without bulk averaging.

Three things in here are load-bearing and easy to get wrong:

1. PROTOCOL MATCH.  Our TCGA patches are 256x256 px resampled from a fixed 128 um x 128 um
   field of view (crop ``128/mpp`` native px, BICUBIC to 256, JPEG q75 4:2:0, no colour
   normalisation).  That geometry is magnification-invariant by construction and it is what
   H-Optimus-0 saw for all 271,710 existing TCGA patches.  Spatial tiles are cut to the SAME
   128 um field -- *not* to the spot diameter and *not* to a convenient pixel size.  The
   renderer is imported from the TCGA extractor rather than reimplemented, so the two stores
   cannot silently diverge.

2. THE SPOT -> 128 um WINDOW COSTS SOMETHING, AND THE COST IS RECORDED.  A Visium spot is
   55 um across on a 100 um pitch: 2376 um^2 of assayed tissue.  A 128 um window is
   16384 um^2 -- ``window_area_ratio()`` = 6.9x more tissue than the transcriptome came from,
   and the window's corners reach 90.5 um, so it also contains tissue assayed by the six
   neighbouring spots.  The image therefore sees strictly more than the target does.  This is
   an honest, quantified cost of protocol matching, not a bug; ``--fov-microns`` makes the
   spot-matched ablation a one-flag change.

3. THE UNIT OF ANALYSIS IS A SPOT, NOT A PATIENT.  The contract wants unique ``patient_ids``;
   here each row is one spot, keyed ``<sample_id>__<barcode>``.  ``slide_ids`` is carried as
   an extra array so splits stay slide-grouped (no slide straddles two partitions) and so the
   per-slide-mean baseline is reproducible from the artifact alone.

The baseline that matters is ``per_slide_mean_baseline``: predict every spot by the mean
expression of its own slide.  Zero parameters.  It is absent from the HEST leaderboard and
from the 46+ models the 2026 *Brief Bioinform* survey catalogues, and it scores well on the
pooled correlation those papers report -- because pooled correlation is dominated by
between-slide variation that slide identity alone already explains.  ``within_slide_r`` is
the metric that cannot be won this way.
"""
from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path
from typing import Sequence

import numpy as np

from ..contracts import ARTIFACT_VERSION, trainable_state_names

# Pinned empirically from the TCGA-UT store; duplicated here only as a guard so a drift in
# the extractor is caught by test_hest.py rather than by a silently incomparable embedding.
FOV_MICRONS = 128.0
OUTPUT_PX = 256
VISIUM_SPOT_DIAMETER_UM = 55.0
VISIUM_SPOT_PITCH_UM = 100.0
EMBED_DIM = 1536
RANDOM_CONTROL_PREFIX = "RANDOM_CONTROL__"

# H-Optimus-0's own timm ``pretrained_cfg``: 224 px input at crop_pct 0.875, crop_mode
# "center".  ``create_transform`` therefore resizes our 256 px patch (a no-op) and then
# CENTRE-CROPS it to 224.  The encoder consequently sees 128 * 224/256 = 112 um, not the
# 128 um we cut.  This applies identically to all 271,710 TCGA patches -- the same transform
# is resolved from the same cfg -- so comparability is unaffected, but the physical field
# under analysis is 112 um and saying "128 um" without this footnote overstates it by 14%.
HOPTIMUS_CROP_PCT = 0.875
HOPTIMUS_INPUT_PX = 224


class HestAdapterError(ValueError):
    pass


# --- geometry: spot -> protocol-matched window ----------------------------------------

def effective_fov_microns(fov_microns: float = FOV_MICRONS,
                          crop_pct: float = HOPTIMUS_CROP_PCT) -> float:
    """The field the encoder actually sees, after its own centre crop: 128 -> 112 um."""
    if not 0 < crop_pct <= 1:
        raise HestAdapterError("crop_pct must lie in (0, 1]")
    return float(fov_microns * crop_pct)


def window_area_ratio(fov_microns: float = FOV_MICRONS,
                      spot_diameter_um: float = VISIUM_SPOT_DIAMETER_UM) -> float:
    """How many times more tissue the image window sees than the assay did.

    The square 128 um window against the circular 55 um spot: 16384 / 2375.8 = 6.90.
    After H-Optimus-0's own centre crop the encoder sees 112 um, i.e. 5.28x -- use
    ``effective_fov_microns()`` for the number that describes what the model was shown.
    Both are reported in the artifact manifest so no downstream reader can forget them.
    """
    if fov_microns <= 0 or spot_diameter_um <= 0:
        raise HestAdapterError("field of view and spot diameter must be positive")
    return float(fov_microns ** 2 / (np.pi * (spot_diameter_um / 2.0) ** 2))


def crop_pixels(mpp: float, fov_microns: float = FOV_MICRONS) -> int:
    """Native level-0 pixels spanning the fixed physical field of view."""
    if not np.isfinite(mpp) or mpp <= 0:
        raise HestAdapterError(f"non-positive microns-per-pixel: {mpp}")
    return int(round(fov_microns / float(mpp)))


def spot_windows(coords: np.ndarray, mpp: float, slide_dims: tuple[int, int],
                 fov_microns: float = FOV_MICRONS) -> tuple[np.ndarray, np.ndarray]:
    """Level-0 top-left corners for windows centred on each spot, plus a keep mask.

    ``coords`` is [spot, 2] of level-0 (x, y) spot centroids.  A spot whose 128 um window
    would run off the slide is dropped rather than clipped or padded -- a clipped window is
    no longer a 128 um field, and silently keeping it would break the protocol match.
    """
    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise HestAdapterError("spot coordinates must be an [n_spots, 2] array")
    crop = crop_pixels(mpp, fov_microns)
    half = crop / 2.0
    corners = np.rint(coords - half).astype(np.int64)
    width, height = int(slide_dims[0]), int(slide_dims[1])
    keep = ((corners[:, 0] >= 0) & (corners[:, 1] >= 0)
            & (corners[:, 0] + crop <= width) & (corners[:, 1] + crop <= height))
    return corners, keep


# --- expression: counts -> targets ----------------------------------------------------

def normalise_expression(counts: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """CPM-style library-size normalisation then log1p, the field-standard spot transform.

    Spots with zero total counts survive as all-zero rows; they are removed upstream by
    ``usable_spots`` rather than silently producing NaNs here.
    """
    counts = np.asarray(counts, dtype=np.float64)
    if counts.ndim != 2:
        raise HestAdapterError("counts must be a [spot, gene] matrix")
    totals = counts.sum(axis=1, keepdims=True)
    scaled = counts * (target_sum / np.maximum(totals, 1.0))
    return np.log1p(scaled).astype(np.float32)


def usable_spots(counts: np.ndarray, min_counts: float = 100.0, min_genes: int = 50) -> np.ndarray:
    """Mask of spots with enough library to carry a target."""
    counts = np.asarray(counts, dtype=np.float64)
    return (counts.sum(axis=1) >= min_counts) & ((counts > 0).sum(axis=1) >= min_genes)


def select_target_genes(expr: np.ndarray, gene_names: Sequence[str], n_genes: int,
                        fit_mask: np.ndarray | None = None) -> np.ndarray:
    """Indices of the most variable genes, chosen on TRAIN spots only.

    Choosing on all spots would leak test-set variance structure into the target definition;
    ``fit_mask`` is therefore not optional in practice, only in signature.
    """
    expr = np.asarray(expr, dtype=np.float32)
    if expr.ndim != 2 or expr.shape[1] != len(gene_names):
        raise HestAdapterError("expression matrix and gene names disagree")
    rows = expr if fit_mask is None else expr[np.asarray(fit_mask, dtype=bool)]
    if rows.shape[0] < 2:
        raise HestAdapterError("need at least two spots to rank gene variance")
    variance = rows.var(axis=0)
    variance[~np.isfinite(variance)] = -np.inf
    order = np.argsort(-variance, kind="stable")
    return np.sort(order[:int(n_genes)])


# --- splits ---------------------------------------------------------------------------

def slide_grouped_split(slide_ids: Sequence[str], seed: int = 42,
                        val_fraction: float = 0.2, test_fraction: float = 0.3) -> np.ndarray:
    """Assign train/val/test so that no slide straddles two partitions.

    Spot-level random splitting would put neighbouring spots from the same slide -- often
    overlapping tissue -- on both sides of the split, and every reported number would be
    contaminated.  Slides are the smallest safe unit here because HEST's ``patient`` field is
    unpopulated for most public 10x samples (6 distinct labels across 44 slides).
    """
    slides = np.asarray(slide_ids).astype(str)
    unique = np.unique(slides)
    if unique.size < 3:
        raise HestAdapterError("need at least three slides to form train/val/test")
    order = np.random.default_rng(seed).permutation(unique.size)
    n_test = max(1, int(round(test_fraction * unique.size)))
    n_val = max(1, int(round(val_fraction * unique.size)))
    if n_test + n_val >= unique.size:
        raise HestAdapterError("val and test fractions leave no training slides")
    label = {}
    for rank, idx in enumerate(order):
        name = unique[idx]
        label[name] = "test" if rank < n_test else ("val" if rank < n_test + n_val else "train")
    return np.asarray([label[s] for s in slides], dtype=str)


# --- baselines ------------------------------------------------------------------------

def per_slide_mean_baseline(targets: np.ndarray, slide_ids: Sequence[str],
                            eval_mask: np.ndarray) -> np.ndarray:
    """The zero-parameter baseline: predict every spot by the mean of its OWN slide.

    Nothing is fitted and nothing is learned -- the prediction for a spot is a statistic of
    the very slide it belongs to, computed on the evaluation partition itself.  That is
    deliberately the most generous possible version of the baseline: if a trained model
    cannot beat a constant-per-slide predictor that was handed the answer key at slide
    granularity, the model has not demonstrated a morphology->molecular channel at all.

    Constant within a slide, so ``within_slide_r`` against it is 0 by construction while
    ``pooled_r`` can be high.  The gap between those two numbers is the whole point.
    """
    targets = np.asarray(targets, dtype=np.float64)
    slides = np.asarray(slide_ids).astype(str)
    eval_mask = np.asarray(eval_mask, dtype=bool)
    if not (len(targets) == len(slides) == len(eval_mask)):
        raise HestAdapterError("targets, slide ids and mask must be aligned")
    prediction = np.full(targets.shape, np.nan, dtype=np.float64)
    for slide in np.unique(slides[eval_mask]):
        rows = eval_mask & (slides == slide)
        prediction[rows] = targets[rows].mean(axis=0, keepdims=True)
    return prediction


def pooled_r(truth: np.ndarray, prediction: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """Per-gene Pearson r over all evaluation spots pooled together.

    This is the number the HEST leaderboard and most of the literature report.  It is
    inflated by between-slide variation, which is exactly what the per-slide-mean baseline
    exploits.
    """
    truth, prediction = np.asarray(truth, dtype=np.float64), np.asarray(prediction, dtype=np.float64)
    if mask is not None:
        rows = np.asarray(mask, dtype=bool)
        truth, prediction = truth[rows], prediction[rows]
    return _column_pearson(truth, prediction)


def within_slide_r(truth: np.ndarray, prediction: np.ndarray, slide_ids: Sequence[str],
                   mask: np.ndarray | None = None) -> np.ndarray:
    """Per-gene Pearson r computed inside each slide, then averaged over slides.

    Between-slide variation is removed, so slide identity buys nothing.  A model that only
    recognises which slide it is looking at scores 0 here.
    """
    truth, prediction = np.asarray(truth, dtype=np.float64), np.asarray(prediction, dtype=np.float64)
    slides = np.asarray(slide_ids).astype(str)
    rows = np.ones(len(truth), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    per_slide = []
    for slide in np.unique(slides[rows]):
        sel = rows & (slides == slide)
        if sel.sum() < 3:
            continue
        per_slide.append(_column_pearson(truth[sel], prediction[sel]))
    if not per_slide:
        raise HestAdapterError("no slide had enough evaluation spots for a within-slide correlation")
    return np.nanmean(np.vstack(per_slide), axis=0)


def _column_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pearson r per column; a constant column yields 0.0, never NaN-by-surprise."""
    a = a - a.mean(axis=0, keepdims=True)
    b = b - b.mean(axis=0, keepdims=True)
    denom = np.sqrt((a ** 2).sum(axis=0) * (b ** 2).sum(axis=0))
    out = np.zeros(a.shape[1], dtype=np.float64)
    good = denom > 1e-12
    out[good] = (a[:, good] * b[:, good]).sum(axis=0) / denom[good]
    return out


def cohort_classifier_auc(a: np.ndarray, b: np.ndarray, seed: int = 0, n_splits: int = 5) -> float:
    """Cross-validated AUC separating two embedding cohorts (e.g. TCGA vs HEST).

    If morphology embeddings from the two cohorts are linearly separable with high AUC, then
    residual batch signal is available to every downstream result that mixes them, and that
    must be declared rather than assumed away.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler

    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1]:
        raise HestAdapterError("cohort matrices must be 2-D with matching feature width")
    x = np.vstack([a, b])
    y = np.concatenate([np.zeros(len(a), dtype=int), np.ones(len(b), dtype=int)])
    scores = np.zeros(len(y), dtype=np.float64)
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for train, test in splitter.split(x, y):
        scaler = StandardScaler().fit(x[train])
        model = LogisticRegression(max_iter=2000, C=1.0)
        model.fit(scaler.transform(x[train]), y[train])
        scores[test] = model.predict_proba(scaler.transform(x[test]))[:, 1]
    return float(roc_auc_score(y, scores))


# --- artifact writers ------------------------------------------------------------------

def _normalise_rows(value: np.ndarray) -> np.ndarray:
    value = np.asarray(value, dtype=np.float32)
    return value / np.maximum(np.linalg.norm(value, axis=1, keepdims=True), 1e-8)


def _atomic_savez(target: Path, **arrays) -> Path:
    import os
    import tempfile

    target.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(dir=str(target.parent), suffix=".npz.tmp", delete=False)
    handle.close()
    np.savez_compressed(handle.name, **arrays)
    written = handle.name if handle.name.endswith(".npz") else handle.name + ".npz"
    os.replace(written, target)
    return target


def write_spatial_artifact(output: str | Path, *, spot_ids: Sequence[str], cancers: Sequence[str],
                           split: Sequence[str], slide_ids: Sequence[str], embeddings: np.ndarray,
                           state_name: str = "wsi_identity", config: dict | None = None,
                           source_provenance: dict | None = None) -> Path:
    """Write spot-level embeddings in the V2 artifact contract.

    ``patient_ids`` holds spot keys, because the pipeline's unit-of-analysis slot is what a
    row means, and here a row is a spot.  ``slide_ids`` rides along so the grouping structure
    survives into every downstream reader.
    """
    ids = np.asarray(spot_ids).astype(str)
    cancer = np.asarray(cancers).astype(str)
    parts = np.asarray(split).astype(str)
    slides = np.asarray(slide_ids).astype(str)
    x = np.asarray(embeddings, dtype=np.float32)
    if len(ids) == 0 or len(set(ids.tolist())) != len(ids):
        raise HestAdapterError("spot ids must be non-empty and unique")
    if not (len(ids) == len(cancer) == len(parts) == len(slides) == len(x)):
        raise HestAdapterError("spot ids, cancers, splits, slides and embeddings must be aligned")
    if x.ndim != 2 or not np.isfinite(x).all():
        raise HestAdapterError("embeddings must be a finite [spot, feature] matrix")
    if "train" not in set(parts.tolist()) or "test" not in set(parts.tolist()):
        raise HestAdapterError("artifact needs explicit train and test spots")
    for slide in np.unique(slides):
        if np.unique(parts[slides == slide]).size != 1:
            raise HestAdapterError(f"slide {slide} straddles more than one partition")

    trained = np.asarray(trainable_state_names(state_name), dtype=str)
    cfg = dict(config or {})
    cfg.setdefault("fov_microns", FOV_MICRONS)
    cfg.setdefault("output_px", OUTPUT_PX)
    cfg.setdefault("spot_diameter_um", VISIUM_SPOT_DIAMETER_UM)
    cfg.setdefault("window_area_ratio", round(window_area_ratio(), 4))
    cfg.setdefault("effective_fov_microns", round(effective_fov_microns(), 4))
    cfg.setdefault("effective_window_area_ratio",
                   round(window_area_ratio(effective_fov_microns()), 4))
    cfg.setdefault("unit_of_analysis", "visium_spot")
    manifest = {
        "artifact_version": ARTIFACT_VERSION, "method": "hest_spatial", "trained_states": trained.tolist(),
        "fit_population": "none_no_fitted_representation", "config": cfg,
        "source_provenance": source_provenance or {},
        "n_slides": int(np.unique(slides).size),
        "cohort_digest": sha256(json.dumps(
            {"patient_ids": ids.tolist(), "cancers": cancer.tolist(), "split": parts.tolist()},
            sort_keys=True).encode()).hexdigest(),
    }
    return _atomic_savez(
        Path(output), patient_ids=ids, cancers=cancer, split=parts, slide_ids=slides,
        trained_states=trained, artifact_version=np.asarray(ARTIFACT_VERSION),
        manifest_json=np.asarray(json.dumps(manifest, sort_keys=True)),
        **{state_name: _normalise_rows(x)})


def write_spatial_targets(output: str | Path, *, spot_ids: Sequence[str], scores: np.ndarray,
                          target_names: Sequence[str], target_groups: Sequence[str] | str,
                          slide_ids: Sequence[str] | None = None, n_random_controls: int = 16,
                          seed: int = 0, config: dict | None = None) -> Path:
    """Write spot-level expression targets in the CALIBRA target contract.

    Adds ``RANDOM_CONTROL__`` columns -- permutations of real target columns, which preserve
    the marginal distribution while destroying the spot correspondence -- because the readers
    split those off into a separate control block and a target file without them offers no
    null to compare against.
    """
    ids = np.asarray(spot_ids).astype(str)
    y = np.asarray(scores, dtype=np.float32)
    names = list(map(str, target_names))
    if len(set(ids.tolist())) != len(ids):
        raise HestAdapterError("target spot ids must be unique")
    if y.ndim != 2 or y.shape != (len(ids), len(names)):
        raise HestAdapterError("scores must be [spot, target] and aligned to ids and names")
    if not np.isfinite(y).all():
        raise HestAdapterError("target scores must be finite")
    stds = y.std(axis=0)
    if (stds < 1e-8).any():
        keep = stds >= 1e-8
        y, names = y[:, keep], [n for n, k in zip(names, keep) if k]
        if not names:
            raise HestAdapterError("every target column was constant; CALIBRA would refuse this file")
    groups = ([str(target_groups)] * len(names) if isinstance(target_groups, str)
              else list(map(str, target_groups)))
    if len(groups) != len(names):
        raise HestAdapterError("target groups must align to target names")

    rng = np.random.default_rng(seed)
    control_cols, control_names = [], []
    for i in range(int(n_random_controls)):
        source = i % y.shape[1]
        control_cols.append(y[rng.permutation(len(ids)), source])
        control_names.append(f"{RANDOM_CONTROL_PREFIX}{names[source]}__perm{i}")
    if control_cols:
        y = np.hstack([y, np.column_stack(control_cols).astype(np.float32)])
        names = names + control_names
        groups = groups + [f"{RANDOM_CONTROL_PREFIX}PERMUTED"] * len(control_names)

    manifest = {"target_kind": "hest_spatial_expression", "config": dict(config or {}),
                "n_random_controls": int(n_random_controls), "seed": int(seed),
                "scores_sha256": sha256(np.ascontiguousarray(y).tobytes()).hexdigest(),
                "patient_id_digest": sha256(json.dumps(ids.tolist(), sort_keys=True).encode()).hexdigest()}
    arrays = {"patient_ids": ids, "scores": y.astype(np.float32),
              "target_names": np.asarray(names, dtype=str), "target_groups": np.asarray(groups, dtype=str),
              "manifest_json": np.asarray(json.dumps(manifest, sort_keys=True))}
    if slide_ids is not None:
        arrays["slide_ids"] = np.asarray(slide_ids).astype(str)
    return _atomic_savez(Path(output), **arrays)
