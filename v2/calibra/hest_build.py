"""Build the HEST spatial cohort into V2 artifacts: scan -> tile+embed -> assemble -> baselines.

Stages are separate subcommands because the expensive one (``tile``) is GPU-bound and
resumable per slide, and because gene selection must happen AFTER the slide-grouped split is
fixed so the target definition never sees test spots.

    python -m morpheus.v2.calibra.hest_build scan     --root ... --metadata ... --output ...
    python -m morpheus.v2.calibra.hest_build tile     --root ... --plan ... --output ...
    python -m morpheus.v2.calibra.hest_build assemble --plan ... --staging ... --output ...
    python -m morpheus.v2.calibra.hest_build baselines --artifact ... --targets ... --output ...

The renderer is imported from the TCGA extractor on purpose: protocol identity with the
271,710 existing TCGA patches is the reason this cohort is comparable at all, and a
copy-pasted renderer would be free to drift.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .hest import (EMBED_DIM, FOV_MICRONS, OUTPUT_PX, HestAdapterError, crop_pixels,
                   normalise_expression, per_slide_mean_baseline, pooled_r, select_target_genes,
                   slide_grouped_split, spot_windows, usable_spots, within_slide_r,
                   window_area_ratio, write_spatial_artifact, write_spatial_targets)


def _render_patch():
    """The exact TCGA renderer -- 128 um field, BICUBIC to 256, JPEG q75 4:2:0 round trip."""
    from ..research.dilution.extract_normal_patches import (JPEG_QUALITY, JPEG_SUBSAMPLING,
                                                            render_patch)
    from ..research.dilution.extract_normal_patches import FOV_MICRONS as TCGA_FOV
    from ..research.dilution.extract_normal_patches import OUTPUT_PX as TCGA_PX
    if float(TCGA_FOV) != FOV_MICRONS or int(TCGA_PX) != OUTPUT_PX:
        raise HestAdapterError(
            f"TCGA patch spec drifted (fov={TCGA_FOV}, px={TCGA_PX}); spatial tiles would no "
            f"longer be protocol-matched to the existing H-Optimus store")
    return render_patch, int(JPEG_QUALITY), int(JPEG_SUBSAMPLING)


def _read_st(path: Path):
    import anndata

    adata = anndata.read_h5ad(path)
    counts = adata.X
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    genes = np.asarray(adata.var_names).astype(str)
    barcodes = np.asarray(adata.obs_names).astype(str)
    return np.asarray(counts, dtype=np.float32), coords, genes, barcodes


# --- stage: scan ----------------------------------------------------------------------

def stage_scan(args) -> None:
    """Fix the slide-grouped split and the gene panel, using train slides only."""
    import pandas as pd

    root = Path(args.root)
    meta = pd.read_csv(args.metadata)
    rows, gene_sets = [], {}
    for r in meta.itertuples():
        st_path = root / "st" / f"{r.id}.h5ad"
        wsi_path = root / "wsis" / f"{r.id}.tif"
        if not st_path.exists() or not wsi_path.exists():
            print(f"[scan] SKIP {r.id}: missing {'st' if not st_path.exists() else 'wsi'}", flush=True)
            continue
        counts, coords, genes, barcodes = _read_st(st_path)
        keep = usable_spots(counts, args.min_counts, args.min_genes)
        gene_sets[r.id] = set(genes.tolist())
        rows.append({"id": r.id, "organ": r.organ, "oncotree_code": r.oncotree_code,
                     "mpp": float(r.mpp), "n_spots": int(counts.shape[0]),
                     "n_usable": int(keep.sum()), "n_genes": int(len(genes)),
                     "crop_px": crop_pixels(float(r.mpp))})
        print(f"[scan] {r.id}: {keep.sum()}/{counts.shape[0]} usable spots, "
              f"{len(genes)} genes, crop {crop_pixels(float(r.mpp))} px", flush=True)
    if len(rows) < 3:
        raise HestAdapterError("need at least three usable slides")
    frame = pd.DataFrame(rows)
    frame["split"] = slide_grouped_split(frame.id.tolist(), seed=args.seed,
                                         val_fraction=args.val_fraction,
                                         test_fraction=args.test_fraction)

    common = set.intersection(*(gene_sets[i] for i in frame.id))
    common = np.asarray(sorted(common), dtype=str)
    print(f"[scan] {len(frame)} slides, gene intersection = {len(common)}", flush=True)
    if len(common) < args.n_target_genes:
        raise HestAdapterError(f"gene intersection {len(common)} < requested panel {args.n_target_genes}")

    train_ids = frame.id[frame.split == "train"].tolist()
    pooled = []
    for sid in train_ids:
        counts, _, genes, _ = _read_st(root / "st" / f"{sid}.h5ad")
        order = {g: i for i, g in enumerate(genes.tolist())}
        counts = counts[:, [order[g] for g in common.tolist()]]
        keep = usable_spots(counts, args.min_counts, args.min_genes)
        if keep.sum() == 0:
            continue
        pooled.append(normalise_expression(counts[keep]))
    if not pooled:
        raise HestAdapterError("no usable training spots for gene selection")
    pooled = np.vstack(pooled)
    panel_idx = select_target_genes(pooled, common.tolist(), args.n_target_genes)
    panel = common[panel_idx]
    print(f"[scan] panel = {len(panel)} genes selected on {len(pooled)} TRAIN spots", flush=True)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out.with_suffix(".slides.csv"), index=False)
    json.dump({"panel_genes": panel.tolist(), "common_genes": int(len(common)),
               "seed": args.seed, "n_target_genes": args.n_target_genes,
               "min_counts": args.min_counts, "min_genes": args.min_genes,
               "fov_microns": FOV_MICRONS, "output_px": OUTPUT_PX,
               "window_area_ratio": round(window_area_ratio(), 4),
               "slides": frame.to_dict(orient="records")}, open(out, "w"), indent=1)
    print(f"[scan] wrote {out}", flush=True)
    print(frame.groupby("split").agg(slides=("id", "size"), spots=("n_usable", "sum")).to_string(), flush=True)


# --- stage: tile ----------------------------------------------------------------------

def stage_tile(args) -> None:
    import openslide

    plan = json.load(open(args.plan))
    panel = np.asarray(plan["panel_genes"], dtype=str)
    root, staging = Path(args.root), Path(args.output)
    staging.mkdir(parents=True, exist_ok=True)
    render_patch, jpeg_q, jpeg_sub = _render_patch()

    pending = [s for s in plan["slides"] if not (staging / f"{s['id']}.npz").exists()]
    print(f"[tile] {len(plan['slides'])} slides, {len(pending)} pending", flush=True)
    if not pending:
        return
    model, transform, torch = _load_encoder(args.device)

    for slide_meta in pending:
        sid, mpp = slide_meta["id"], float(slide_meta["mpp"])
        try:
            counts, coords, genes, barcodes = _read_st(root / "st" / f"{sid}.h5ad")
            order = {g: i for i, g in enumerate(genes.tolist())}
            if not set(panel.tolist()).issubset(order):
                raise HestAdapterError(f"{sid} is missing panel genes")
            panel_counts = counts[:, [order[g] for g in panel.tolist()]]
            keep = usable_spots(counts, plan["min_counts"], plan["min_genes"])
            slide = openslide.OpenSlide(str(root / "wsis" / f"{sid}.tif"))
            corners, on_slide = spot_windows(coords, mpp, slide.dimensions, FOV_MICRONS)
            keep = keep & on_slide
            idx = np.flatnonzero(keep)
            if idx.size == 0:
                raise HestAdapterError("no spot survived the usable/on-slide filters")
            crop = crop_pixels(mpp, FOV_MICRONS)
            vectors = _embed_windows(slide, corners[idx], crop, render_patch, model, transform,
                                     torch, args.device, args.batch_size, args.readers)
            slide.close()
            if vectors.shape != (idx.size, EMBED_DIM) or not np.isfinite(vectors).all():
                raise HestAdapterError(f"bad embedding block {vectors.shape}")
            np.savez_compressed(
                staging / f"{sid}.npz", embeddings=vectors,
                spot_ids=np.asarray([f"{sid}__{b}" for b in barcodes[idx]], dtype=str),
                panel_counts=panel_counts[idx].astype(np.float32),
                total_counts=counts[idx].sum(axis=1).astype(np.float32),
                coords=coords[idx].astype(np.float64),
                meta=np.asarray(json.dumps({
                    "id": sid, "organ": slide_meta["organ"], "oncotree_code": slide_meta["oncotree_code"],
                    "split": slide_meta["split"], "mpp": mpp, "crop_px": crop,
                    "fov_microns": FOV_MICRONS, "output_px": OUTPUT_PX,
                    "jpeg_quality": jpeg_q, "jpeg_subsampling": "4:2:0",
                    "n_spots_total": int(counts.shape[0]), "n_spots_kept": int(idx.size),
                    "n_dropped_off_slide": int((~on_slide).sum())})))
            print(f"[tile] {sid}: {idx.size} spots embedded "
                  f"({(~on_slide).sum()} dropped off-slide), crop {crop} px", flush=True)
        except Exception as exc:  # one bad slide must not abort the sweep
            json.dump({"id": sid, "error": repr(exc)[:500]},
                      open(staging / f"{sid}.failed.json", "w"), indent=1)
            print(f"[tile] FAIL {sid}: {exc!r}", flush=True)


def _load_encoder(device: str):
    import timm
    import torch
    from timm.data import create_transform, resolve_data_config

    model = timm.create_model("hf-hub:bioptimus/H-optimus-0", pretrained=True).eval().to(device)
    transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model),
                                 is_training=False)
    return model, transform, torch


def _embed_windows(slide, corners, crop_px, render_patch, model, transform, torch,
                   device: str, batch_size: int, readers: int = 12) -> np.ndarray:
    """Embed one slide's windows, overlapping slide I/O with GPU work.

    Reading a 744x744 level-0 region off NFS costs ~155 ms cold, so a serial loop is
    I/O-bound at ~6 spots/s -- six hours for this cohort.  ``openslide_read_region`` is
    thread-safe and releases the GIL, so the reads are issued from a small thread pool while
    the GPU works on the previous batch.  Patch order is preserved, so the embedding rows
    still line up with ``corners``.
    """
    from concurrent.futures import ThreadPoolExecutor

    out = []
    with ThreadPoolExecutor(max_workers=readers) as pool:
        pending = None
        for start in range(0, len(corners), batch_size):
            chunk = corners[start:start + batch_size]
            future = pool.map(lambda c: render_patch(slide, int(c[0]), int(c[1]), crop_px), chunk)
            if pending is not None:
                out.append(_forward(model, transform, torch, list(pending), device))
            pending = future
        if pending is not None:
            out.append(_forward(model, transform, torch, list(pending), device))
    return np.concatenate(out, axis=0).astype(np.float32) if out else np.zeros((0, EMBED_DIM), np.float32)


def _forward(model, transform, torch, images, device: str) -> np.ndarray:
    """One H-Optimus-0 batch, in exactly the precision the TCGA store used."""
    batch = torch.stack([transform(im) for im in images]).to(device)
    with torch.inference_mode(), torch.autocast("cuda", dtype=torch.float16,
                                                enabled=device.startswith("cuda")):
        values = model(batch)
    if isinstance(values, (tuple, list)):
        values = values[0]
    return values.detach().float().cpu().numpy()


# --- stage: assemble ------------------------------------------------------------------

def stage_assemble(args) -> None:
    plan = json.load(open(args.plan))
    panel = np.asarray(plan["panel_genes"], dtype=str)
    staging, out = Path(args.staging), Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    ids, cancers, splits, slides, embeds, counts = [], [], [], [], [], []
    for slide_meta in plan["slides"]:
        path = staging / f"{slide_meta['id']}.npz"
        if not path.exists():
            print(f"[assemble] SKIP {slide_meta['id']}: not tiled", flush=True)
            continue
        blob = np.load(path, allow_pickle=False)
        meta = json.loads(str(blob["meta"]))
        n = len(blob["spot_ids"])
        ids.append(blob["spot_ids"].astype(str))
        embeds.append(blob["embeddings"])
        counts.append(blob["panel_counts"])
        cancers.append(np.full(n, str(meta["oncotree_code"]), dtype=object))
        splits.append(np.full(n, str(meta["split"]), dtype=object))
        slides.append(np.full(n, str(meta["id"]), dtype=object))
    if not ids:
        raise HestAdapterError("nothing tiled yet")
    ids = np.concatenate(ids)
    embeds = np.vstack(embeds)
    counts = np.vstack(counts)
    cancers = np.concatenate(cancers).astype(str)
    splits = np.concatenate(splits).astype(str)
    slides = np.concatenate(slides).astype(str)
    print(f"[assemble] {len(ids)} spots across {len(np.unique(slides))} slides, "
          f"embeddings {embeds.shape}", flush=True)

    expression = normalise_expression(counts)
    artifact = write_spatial_artifact(
        out / "hest_spatial_hoptimus.npz", spot_ids=ids, cancers=cancers, split=splits,
        slide_ids=slides, embeddings=embeds, state_name="wsi_identity",
        config={"cohort": "HEST-1k human cancer Visium", "encoder": "H-optimus-0",
                "n_slides": int(np.unique(slides).size), "n_panel_genes": int(len(panel)),
                "licence": "per-sample CC BY 4.0 subset of HEST-1k (CC BY-NC-SA 4.0 collection)"},
        source_provenance={"dataset": "MahmoodLab/hest", "metadata": "HEST_v1_3_0.csv"})
    targets = write_spatial_targets(
        out / "hest_spatial_targets.npz", spot_ids=ids, scores=expression,
        target_names=panel.tolist(), target_groups="HEST_SPOT_EXPRESSION", slide_ids=slides,
        n_random_controls=args.n_random_controls, seed=plan.get("seed", 42),
        config={"normalisation": "log1p CPM(1e4)", "n_panel_genes": int(len(panel))})
    print(f"[assemble] wrote {artifact}\n[assemble] wrote {targets}", flush=True)


# --- stage: baselines -----------------------------------------------------------------

def stage_baselines(args) -> None:
    from sklearn.linear_model import Ridge

    art = np.load(args.artifact, allow_pickle=False)
    tgt = np.load(args.targets, allow_pickle=False)
    ids, slides, splits = art["patient_ids"].astype(str), art["slide_ids"].astype(str), art["split"].astype(str)
    x = art["wsi_identity"].astype(np.float64)
    names = tgt["target_names"].astype(str)
    real = ~np.char.startswith(names, "RANDOM_CONTROL__")
    order = {p: i for i, p in enumerate(tgt["patient_ids"].astype(str).tolist())}
    y = tgt["scores"].astype(np.float64)[[order[p] for p in ids.tolist()]][:, real]
    names = names[real]
    train, test = splits == "train", splits == "test"
    print(f"[baselines] {train.sum()} train / {test.sum()} test spots, {y.shape[1]} genes", flush=True)

    results = {}

    slide_mean = per_slide_mean_baseline(y, slides, test)
    results["per_slide_mean"] = _score(y, slide_mean, slides, test)

    global_mean = np.tile(y[train].mean(axis=0), (len(y), 1))
    results["global_mean"] = _score(y, global_mean, slides, test)

    ridge = Ridge(alpha=args.ridge_alpha).fit(x[train], y[train])
    results["ridge_hoptimus"] = _score(y, ridge.predict(x), slides, test)

    for name, res in results.items():
        print(f"[baselines] {name:20s} pooled_r={res['pooled_r_mean']:.4f}  "
              f"within_slide_r={res['within_slide_r_mean']:.4f}", flush=True)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump({"results": results, "n_train_spots": int(train.sum()), "n_test_spots": int(test.sum()),
               "n_genes": int(y.shape[1]), "n_test_slides": int(np.unique(slides[test]).size),
               "ridge_alpha": args.ridge_alpha,
               "note": "pooled_r is the metric the HEST leaderboard reports; within_slide_r "
                       "removes between-slide variation, which the zero-parameter per-slide "
                       "mean baseline exploits."},
              open(out, "w"), indent=1)
    print(f"[baselines] wrote {out}", flush=True)


def _score(y, prediction, slides, mask) -> dict:
    pooled = pooled_r(y, prediction, mask)
    within = within_slide_r(y, prediction, slides, mask)
    return {"pooled_r_mean": float(np.nanmean(pooled)), "pooled_r_median": float(np.nanmedian(pooled)),
            "within_slide_r_mean": float(np.nanmean(within)),
            "within_slide_r_median": float(np.nanmedian(within))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)

    s = sub.add_parser("scan")
    s.add_argument("--root", required=True); s.add_argument("--metadata", required=True)
    s.add_argument("--output", required=True); s.add_argument("--seed", type=int, default=42)
    s.add_argument("--val-fraction", type=float, default=0.2)
    s.add_argument("--test-fraction", type=float, default=0.3)
    s.add_argument("--n-target-genes", type=int, default=50)
    s.add_argument("--min-counts", type=float, default=100.0)
    s.add_argument("--min-genes", type=int, default=50)
    s.set_defaults(func=stage_scan)

    t = sub.add_parser("tile")
    t.add_argument("--root", required=True); t.add_argument("--plan", required=True)
    t.add_argument("--output", required=True); t.add_argument("--device", default="cuda")
    t.add_argument("--batch-size", type=int, default=32)
    t.add_argument("--readers", type=int, default=12,
                   help="threads issuing slide reads; the stage is NFS-latency bound, not GPU bound")
    t.set_defaults(func=stage_tile)

    a = sub.add_parser("assemble")
    a.add_argument("--plan", required=True); a.add_argument("--staging", required=True)
    a.add_argument("--output", required=True)
    a.add_argument("--n-random-controls", type=int, default=16)
    a.set_defaults(func=stage_assemble)

    b = sub.add_parser("baselines")
    b.add_argument("--artifact", required=True); b.add_argument("--targets", required=True)
    b.add_argument("--output", required=True); b.add_argument("--ridge-alpha", type=float, default=1000.0)
    b.set_defaults(func=stage_baselines)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
