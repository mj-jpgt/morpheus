"""Tests for the HEST spatial adapter.

The geometry tests ARE the protocol-match self-test: if the spot -> 128 um window mapping
drifts, the spatial embeddings stop being comparable to the 271,710 TCGA patches and every
cross-modality number becomes meaningless without anything failing loudly.

The per-slide-mean tests encode the finding the cohort exists to demonstrate: a
zero-parameter predictor scores well on the pooled correlation the literature reports, and
exactly zero on the within-slide correlation it does not.
"""
from __future__ import annotations

import numpy as np

from morpheus.v2.calibra.hest import (FOV_MICRONS, OUTPUT_PX, VISIUM_SPOT_DIAMETER_UM,
                                      HestAdapterError, cohort_classifier_auc, crop_pixels,
                                      effective_fov_microns, normalise_expression,
                                      per_slide_mean_baseline, pooled_r,
                                      select_target_genes, slide_grouped_split, spot_windows,
                                      usable_spots, within_slide_r, window_area_ratio,
                                      write_spatial_artifact, write_spatial_targets)
from morpheus.v2.contracts import declared_states, validate_artifact


def _spot_grid(n_side=6, pitch_um=100.0, mpp=0.25, origin=5000.0):
    """A square Visium-like lattice in level-0 pixels."""
    step = pitch_um / mpp
    xs = origin + step * np.arange(n_side)
    grid = np.stack(np.meshgrid(xs, xs, indexing="ij"), axis=-1).reshape(-1, 2)
    return grid


def _cohort(n_slides=6, per_slide=40, n_genes=12, seed=0, slide_shift=3.0, image_signal=1.0):
    """Spots carrying a between-slide offset plus a within-slide morphology->expression link."""
    rng = np.random.default_rng(seed)
    slides, x, y = [], [], []
    for s in range(n_slides):
        offset = rng.normal(scale=slide_shift, size=n_genes)
        latent = rng.normal(size=(per_slide, 3))
        emb = np.hstack([latent, rng.normal(size=(per_slide, 5))])
        weights = rng.normal(size=(3, n_genes))
        y.append(offset + image_signal * (latent @ weights) + rng.normal(scale=0.5, size=(per_slide, n_genes)))
        x.append(emb)
        slides.extend([f"S{s}"] * per_slide)
    return np.vstack(x), np.vstack(y), np.asarray(slides)


# --- geometry: the protocol match -----------------------------------------------------

def test_window_area_ratio_is_the_documented_cost():
    """A 128 um window sees 6.9x the tissue a 55 um Visium spot assayed."""
    ratio = window_area_ratio()
    assert abs(ratio - 6.90) < 0.01, f"documented cost drifted: {ratio}"
    # Matching the spot instead of the protocol would be 1.0 by definition.
    assert abs(window_area_ratio(VISIUM_SPOT_DIAMETER_UM, VISIUM_SPOT_DIAMETER_UM) - 4 / np.pi) < 1e-9


def test_effective_field_accounts_for_the_encoder_centre_crop():
    """H-Optimus-0's own pretrained_cfg is 224 px at crop_pct 0.875, so create_transform
    centre-crops our 256 px patch to 224 and the encoder sees 112 um, not 128.  Identical for
    the TCGA store, so comparability holds -- but the analysed field is 112 um."""
    assert abs(effective_fov_microns() - 112.0) < 1e-9
    assert abs(window_area_ratio(effective_fov_microns()) - 5.28) < 0.01
    assert effective_fov_microns(crop_pct=1.0) == FOV_MICRONS
    for bad in (0.0, -0.5, 1.5):
        try:
            effective_fov_microns(crop_pct=bad)
        except HestAdapterError:
            continue
        raise AssertionError(f"crop_pct={bad} should have raised")


def test_crop_pixels_is_magnification_invariant_in_microns():
    """The physical field is fixed; only the native pixel count moves with mpp."""
    assert crop_pixels(0.25) == 512      # 40x
    assert crop_pixels(0.50) == 256      # 20x, no resampling at all
    assert crop_pixels(0.171674) == 746  # the finest slide in the selected cohort
    for bad in (0.0, -0.25, np.nan, np.inf):
        try:
            crop_pixels(bad)
        except HestAdapterError:
            continue
        raise AssertionError(f"mpp={bad} should have raised")


def test_spot_windows_centre_on_the_spot_and_drop_off_slide():
    mpp = 0.25
    coords = _spot_grid(mpp=mpp)
    crop = crop_pixels(mpp)
    corners, keep = spot_windows(coords, mpp, (20000, 20000))
    centres = corners + crop / 2.0
    assert np.abs(centres - coords).max() <= 0.5, "window is not centred on the spot"
    assert keep.all(), "fixture should sit entirely on the slide"

    # A spot near the edge cannot yield a full 128 um field, so it must be dropped, not clipped.
    edge = np.array([[10.0, 10.0], [9000.0, 9000.0], [19995.0, 19995.0]])
    _, keep_edge = spot_windows(edge, mpp, (20000, 20000))
    assert keep_edge.tolist() == [False, True, False], keep_edge.tolist()


def test_spot_windows_rejects_malformed_coordinates():
    for bad in (np.zeros(5), np.zeros((5, 3))):
        try:
            spot_windows(bad, 0.25, (100, 100))
        except HestAdapterError:
            continue
        raise AssertionError(f"shape {bad.shape} should have raised")


def test_adapter_constants_match_the_tcga_extractor():
    """REGRESSION GUARD: the whole cohort is only comparable because the field of view is
    identical to the TCGA store.  If either side is edited alone, this fails."""
    from morpheus.v2.research.dilution import extract_normal_patches as tcga

    assert float(tcga.FOV_MICRONS) == FOV_MICRONS, (tcga.FOV_MICRONS, FOV_MICRONS)
    assert int(tcga.OUTPUT_PX) == OUTPUT_PX, (tcga.OUTPUT_PX, OUTPUT_PX)
    assert int(tcga.JPEG_QUALITY) == 75 and int(tcga.JPEG_SUBSAMPLING) == 2


# --- expression -----------------------------------------------------------------------

def test_normalise_expression_removes_library_size():
    rng = np.random.default_rng(0)
    profile = rng.integers(0, 50, size=20).astype(float)
    counts = np.vstack([profile, 10.0 * profile])  # same biology, 10x the depth
    assert counts[1].sum() > counts[0].sum() * 5, "fixture has no depth difference; test is vacuous"
    out = normalise_expression(counts)
    assert np.abs(out[0] - out[1]).max() < 1e-4, "depth survived normalisation"
    assert np.isfinite(out).all()


def test_usable_spots_drops_empty_spots():
    counts = np.zeros((4, 200))
    counts[0, :100] = 5.0    # deep and broad -> keep
    counts[1, :10] = 1.0     # too few counts and genes -> drop
    counts[2, :60] = 3.0     # keep
    keep = usable_spots(counts, min_counts=100.0, min_genes=50)
    assert keep.tolist() == [True, False, True, False], keep.tolist()


def test_select_target_genes_uses_only_the_fit_mask():
    """REGRESSION GUARD: selecting the panel on all spots leaks test-set variance structure
    into the definition of the target itself."""
    rng = np.random.default_rng(1)
    expr = rng.normal(scale=0.01, size=(100, 8))
    fit = np.zeros(100, dtype=bool)
    fit[:50] = True
    expr[:50, 2] += rng.normal(scale=5.0, size=50)    # variable in TRAIN only
    expr[50:, 6] += rng.normal(scale=5.0, size=50)    # variable in HELD-OUT only
    picked = select_target_genes(expr, [f"g{i}" for i in range(8)], 1, fit_mask=fit)
    assert picked.tolist() == [2], picked.tolist()
    leaked = select_target_genes(expr, [f"g{i}" for i in range(8)], 1)
    assert 6 in leaked.tolist() or 2 in leaked.tolist(), "fixture is not discriminating; test is vacuous"


# --- splits ---------------------------------------------------------------------------

def test_slide_grouped_split_never_straddles_a_slide():
    slides = np.repeat([f"S{i}" for i in range(10)], 7)
    split = slide_grouped_split(slides, seed=3)
    for slide in np.unique(slides):
        assert np.unique(split[slides == slide]).size == 1, f"{slide} straddles partitions"
    assert set(split.tolist()) == {"train", "val", "test"}, set(split.tolist())


def test_slide_grouped_split_refuses_impossible_geometry():
    try:
        slide_grouped_split(["a", "a", "b"], seed=0)
    except HestAdapterError:
        pass
    else:
        raise AssertionError("two slides cannot make three partitions")
    try:
        slide_grouped_split(np.repeat([f"S{i}" for i in range(4)], 3), val_fraction=0.5, test_fraction=0.5)
    except HestAdapterError:
        return
    raise AssertionError("val+test consuming every slide should have raised")


# --- the baseline that matters --------------------------------------------------------

def test_per_slide_mean_is_constant_within_a_slide():
    x, y, slides = _cohort()
    mask = np.ones(len(y), dtype=bool)
    prediction = per_slide_mean_baseline(y, slides, mask)
    for slide in np.unique(slides):
        rows = prediction[slides == slide]
        assert np.abs(rows - rows[0]).max() < 1e-9, "per-slide mean varies inside a slide"


def test_per_slide_mean_scores_high_pooled_and_zero_within_slide():
    """THE POINT OF THE COHORT: a zero-parameter predictor wins on the metric the HEST
    leaderboard reports, and scores exactly 0 on the one that removes slide identity."""
    x, y, slides = _cohort(slide_shift=6.0, image_signal=0.3)
    mask = np.ones(len(y), dtype=bool)
    prediction = per_slide_mean_baseline(y, slides, mask)
    pooled = np.nanmean(pooled_r(y, prediction, mask))
    within = np.nanmean(within_slide_r(y, prediction, slides, mask))
    assert pooled > 0.7, f"fixture lacks between-slide variation; test is vacuous ({pooled})"
    assert abs(within) < 1e-9, f"a constant-per-slide predictor cannot correlate within a slide ({within})"


def test_within_slide_r_rewards_a_predictor_the_pooled_metric_cannot_distinguish():
    """A model that tracks within-slide structure beats the slide mean only on within_slide_r."""
    x, y, slides = _cohort(slide_shift=6.0, image_signal=1.0, seed=5)
    mask = np.ones(len(y), dtype=bool)
    oracle = y + np.random.default_rng(0).normal(scale=0.3, size=y.shape)
    slide_mean = per_slide_mean_baseline(y, slides, mask)
    assert np.nanmean(within_slide_r(y, oracle, slides, mask)) > 0.8
    assert np.nanmean(within_slide_r(y, slide_mean, slides, mask)) < 1e-9
    # Both look respectable pooled, which is exactly why pooled alone is not evidence.
    assert np.nanmean(pooled_r(y, slide_mean, mask)) > 0.6


def test_pooled_r_matches_numpy_on_known_data():
    rng = np.random.default_rng(2)
    a = rng.normal(size=(200, 3))
    b = 0.7 * a + rng.normal(scale=0.5, size=(200, 3))
    expected = np.array([np.corrcoef(a[:, i], b[:, i])[0, 1] for i in range(3)])
    np.testing.assert_allclose(pooled_r(a, b), expected, rtol=1e-10, atol=1e-12)
    # A constant column yields 0.0 rather than a NaN surprise.
    assert pooled_r(a, np.ones_like(b))[0] == 0.0


def test_within_slide_r_ignores_between_slide_offsets():
    rng = np.random.default_rng(4)
    slides = np.repeat(["a", "b", "c"], 60)
    truth = rng.normal(size=(180, 2))
    shifted = truth + np.repeat(np.array([[10.0, -8.0], [0.0, 0.0], [-5.0, 7.0]]), 60, axis=0)
    within = within_slide_r(truth, shifted, slides)
    np.testing.assert_allclose(within, np.ones(2), rtol=1e-9, atol=1e-9)
    assert np.nanmean(pooled_r(truth, shifted)) < 0.9, "fixture has no offset; test is vacuous"


def test_per_slide_mean_rejects_misaligned_inputs():
    x, y, slides = _cohort(n_slides=3, per_slide=5)
    try:
        per_slide_mean_baseline(y, slides[:-1], np.ones(len(y), dtype=bool))
    except HestAdapterError:
        return
    raise AssertionError("misaligned slide ids should have raised")


# --- cohort classifier control --------------------------------------------------------

def test_cohort_classifier_auc_is_chance_on_one_distribution():
    rng = np.random.default_rng(7)
    a, b = rng.normal(size=(150, 10)), rng.normal(size=(150, 10))
    auc = cohort_classifier_auc(a, b, seed=0)
    assert 0.35 < auc < 0.65, f"same distribution should be near chance, got {auc}"


def test_cohort_classifier_auc_detects_a_batch_shift():
    rng = np.random.default_rng(8)
    a = rng.normal(size=(150, 10))
    b = rng.normal(size=(150, 10)) + 2.0
    auc = cohort_classifier_auc(a, b, seed=0)
    assert auc > 0.95, f"a two-sigma shift must be detectable, got {auc}"


# --- artifact contract ----------------------------------------------------------------

def _artifact_inputs(n_slides=4, per_slide=10, dim=6):
    rng = np.random.default_rng(0)
    slides = np.repeat([f"S{i}" for i in range(n_slides)], per_slide)
    ids = np.asarray([f"{s}__spot{i}" for i, s in enumerate(slides)])
    split = np.where(slides == "S0", "test", np.where(slides == "S1", "val", "train"))
    cancers = np.repeat(["COAD", "IDC", "COAD", "PRAD"], per_slide)
    return ids, cancers, split, slides, rng.normal(size=(len(ids), dim)).astype(np.float32)


def test_artifact_round_trips_through_the_real_validator(tmp_path):
    ids, cancers, split, slides, emb = _artifact_inputs()
    path = write_spatial_artifact(tmp_path / "a.npz", spot_ids=ids, cancers=cancers, split=split,
                                  slide_ids=slides, embeddings=emb)
    report = validate_artifact(path)
    assert report["n_patients"] == len(ids)
    assert report["trained_states"] == ["wsi_identity"]
    assert report["has_manifest"]
    with np.load(path, allow_pickle=False) as raw:
        assert declared_states(raw) == frozenset({"wsi_identity"})
        assert raw["split"].astype(str).tolist() == split.tolist()
        assert raw["slide_ids"].astype(str).tolist() == slides.tolist()
        rows = raw["wsi_identity"]
        np.testing.assert_allclose(np.linalg.norm(rows, axis=1), 1.0, rtol=1e-5, atol=1e-5)
        manifest = __import__("json").loads(str(raw["manifest_json"]))
        assert manifest["config"]["fov_microns"] == FOV_MICRONS
        assert abs(manifest["config"]["window_area_ratio"] - 6.90) < 0.01
        assert manifest["n_slides"] == 4


def test_atomic_write_leaves_no_temporary_files(tmp_path):
    """np.savez appends .npz when the name lacks it, so the NamedTemporaryFile placeholder is
    not the file that gets renamed and was being left behind as 0-byte litter."""
    ids, cancers, split, slides, emb = _artifact_inputs()
    write_spatial_artifact(tmp_path / "a.npz", spot_ids=ids, cancers=cancers, split=split,
                           slide_ids=slides, embeddings=emb)
    left = sorted(p.name for p in tmp_path.iterdir())
    assert left == ["a.npz"], left


def test_artifact_refuses_a_slide_that_straddles_partitions(tmp_path):
    """REGRESSION GUARD: spot-level splitting puts neighbouring -- often overlapping -- tissue
    on both sides of the split and silently contaminates every reported number."""
    ids, cancers, split, slides, emb = _artifact_inputs()
    split = split.copy()
    split[0] = "train"  # S0 is otherwise entirely test
    assert np.unique(split[slides == "S0"]).size == 2, "fixture does not straddle; test is vacuous"
    try:
        write_spatial_artifact(tmp_path / "b.npz", spot_ids=ids, cancers=cancers, split=split,
                               slide_ids=slides, embeddings=emb)
    except HestAdapterError:
        return
    raise AssertionError("a straddling slide should have raised")


def test_artifact_refuses_duplicate_or_misaligned_spots(tmp_path):
    ids, cancers, split, slides, emb = _artifact_inputs()
    dup = ids.copy()
    dup[1] = dup[0]
    for kwargs in ({"spot_ids": dup}, {"embeddings": emb[:-1]}, {"split": np.full(len(ids), "train")}):
        base = {"spot_ids": ids, "cancers": cancers, "split": split, "slide_ids": slides,
                "embeddings": emb}
        base.update(kwargs)
        try:
            write_spatial_artifact(tmp_path / "c.npz", **base)
        except HestAdapterError:
            continue
        raise AssertionError(f"{list(kwargs)} should have raised")


def test_targets_carry_a_random_control_block(tmp_path):
    ids, _, _, slides, _ = _artifact_inputs()
    rng = np.random.default_rng(0)
    scores = rng.normal(size=(len(ids), 5)).astype(np.float32)
    path = write_spatial_targets(tmp_path / "t.npz", spot_ids=ids, scores=scores,
                                 target_names=[f"GENE{i}" for i in range(5)],
                                 target_groups="HEST_SPOT_EXPRESSION", slide_ids=slides,
                                 n_random_controls=4, seed=0)
    with np.load(path, allow_pickle=False) as raw:
        names = raw["target_names"].astype(str)
        groups = raw["target_groups"].astype(str)
        control = np.char.startswith(names, "RANDOM_CONTROL__")
        assert control.sum() == 4, names.tolist()
        assert raw["scores"].shape == (len(ids), 9)
        assert set(groups[~control].tolist()) == {"HEST_SPOT_EXPRESSION"}
        # A permuted control keeps the marginal distribution and destroys the correspondence.
        real = raw["scores"][:, ~control]
        ctrl = raw["scores"][:, control]
        np.testing.assert_allclose(np.sort(ctrl[:, 0]), np.sort(real[:, 0]), rtol=1e-6, atol=1e-6)
        assert raw["scores"].std(axis=0).min() > 1e-8, "CALIBRA would refuse a constant column"


def test_cancer_label_falls_back_to_organ():
    """REGRESSION GUARD: 14 of 44 selected slides carry no oncotree code.  Letting those
    become the string "nan" invents a cancer type that pools unrelated tumours and then gets
    used as a confound covariate downstream."""
    from morpheus.v2.calibra.hest_build import _cancer_label

    assert _cancer_label({"oncotree_code": "COAD", "organ": "Bowel"}) == "COAD"
    for missing in ("nan", "", "None", "UNKNOWN", "TODO"):
        label = _cancer_label({"oncotree_code": missing, "organ": "Lymph node"})
        assert label == "ORGAN_LYMPH_NODE", (missing, label)
    assert _cancer_label({"oncotree_code": "nan", "organ": ""}) == "UNSPECIFIED"


def test_targets_drop_constant_columns(tmp_path):
    ids, _, _, _, _ = _artifact_inputs()
    scores = np.random.default_rng(0).normal(size=(len(ids), 3)).astype(np.float32)
    scores[:, 1] = 4.0
    path = write_spatial_targets(tmp_path / "t2.npz", spot_ids=ids, scores=scores,
                                 target_names=["A", "B", "C"], target_groups="G",
                                 n_random_controls=0, seed=0)
    with np.load(path, allow_pickle=False) as raw:
        assert raw["target_names"].astype(str).tolist() == ["A", "C"]
