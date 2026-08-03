"""Patch + embed TCGA solid-tissue-NORMAL slides to the pinned TCGA-UT patch spec.

The tumour side of the study (271,710 H-Optimus-0 patch tokens) was sampled ONLY
inside pathologist-drawn tumour polygons.  This script produces the *other* half
of a dilution experiment: genuinely non-tumour tissue from the same cohort,
scanners and stain distribution, patched to the identical physical spec so the
two token populations are directly mixable.

Spec (pinned empirically from the TCGA-UT store, see DILUTION_CURVE.md):
  * 256x256 px output from a fixed 128 um x 128 um field of view at native
    resolution, i.e. crop ``round(128 / mpp)`` native pixels then resample to
    256x256.  Magnification-invariant by construction.
  * JPEG quality 75, 4:2:0 chroma subsampling, encoded and re-decoded so the
    encoder sees the same compression artefacts the tumour patches carry.
  * H-Optimus-0 via ``timm``, exactly the transform the tumour store used
    (``create_transform(**resolve_data_config(cfg), is_training=False)``).

Difference from the tumour side, and it is the point of the experiment: there is
no polygon here, so patches are drawn from a whole-slide tissue mask -- which is
precisely what an external cohort without annotations would give us.

Resumable: one ``.npz`` per slide under ``--output-dir/staging``; the ``.svs`` is
deleted after a slide is embedded so peak disk stays bounded.
"""
from __future__ import annotations

import argparse
import io
import json
import os
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

FOV_MICRONS = 128.0
OUTPUT_PX = 256
JPEG_QUALITY = 75
JPEG_SUBSAMPLING = 2  # 4:2:0
TISSUE_FRACTION = 0.50


def slide_mpp(slide) -> float:
    """Microns-per-pixel at level 0, refusing a slide that cannot state it."""
    for key in ("openslide.mpp-x", "aperio.MPP", "aperio.AppMag"):
        value = slide.properties.get(key)
        if key == "aperio.AppMag":
            break
        if value:
            return float(value)
    power = slide.properties.get("openslide.objective-power") or slide.properties.get("aperio.AppMag")
    if power:
        # TCGA Aperio: 40x ~ 0.25 um/px, 20x ~ 0.50 um/px.
        return 10.0 / float(power)
    raise ValueError("slide does not declare mpp or objective power")


def tissue_mask(slide, downsample_level: int) -> np.ndarray:
    """Saturation/brightness tissue mask on a low-resolution level.

    A whole-tissue pipeline has no polygon, so this stands in for one.  It is
    deliberately the plain thing an external-cohort pipeline would do.
    """
    thumb = slide.read_region((0, 0), downsample_level, slide.level_dimensions[downsample_level]).convert("RGB")
    hsv = np.asarray(thumb.convert("HSV"), dtype=np.uint8)
    saturation, value = hsv[..., 1], hsv[..., 2]
    return (saturation > 30) & (value > 40) & (value < 235)


def sample_patch_grid(slide, crop_px: int, mask: np.ndarray, mask_downsample: float,
                      n_patches: int, rng: np.random.Generator) -> list[tuple[int, int]]:
    """Non-overlapping grid positions whose tissue fraction clears the threshold."""
    width, height = slide.dimensions
    step = crop_px
    xs = np.arange(0, max(width - crop_px, 1), step)
    ys = np.arange(0, max(height - crop_px, 1), step)
    cell = max(int(round(crop_px / mask_downsample)), 1)
    candidates = []
    for y in ys:
        my = int(round(y / mask_downsample))
        row = mask[my:my + cell]
        if row.size == 0:
            continue
        for x in xs:
            mx = int(round(x / mask_downsample))
            block = row[:, mx:mx + cell]
            if block.size and block.mean() >= TISSUE_FRACTION:
                candidates.append((int(x), int(y)))
    if not candidates:
        return []
    if len(candidates) <= n_patches:
        return candidates
    index = rng.choice(len(candidates), size=n_patches, replace=False)
    return [candidates[i] for i in sorted(index)]


def render_patch(slide, x: int, y: int, crop_px: int) -> Image.Image:
    """Read the native field of view, resample to 256, round-trip through JPEG."""
    region = slide.read_region((x, y), 0, (crop_px, crop_px)).convert("RGB")
    if region.size != (OUTPUT_PX, OUTPUT_PX):
        region = region.resize((OUTPUT_PX, OUTPUT_PX), Image.BICUBIC)
    buffer = io.BytesIO()
    region.save(buffer, format="JPEG", quality=JPEG_QUALITY, subsampling=JPEG_SUBSAMPLING)
    buffer.seek(0)
    return Image.open(buffer).convert("RGB")


def load_encoder(device: str):
    import timm
    import torch
    from timm.data import create_transform, resolve_data_config

    model = timm.create_model("hf-hub:bioptimus/H-optimus-0", pretrained=True).eval().to(device)
    transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model), is_training=False)
    return model, transform, torch


def embed(model, transform, torch, images: list[Image.Image], device: str, batch_size: int) -> np.ndarray:
    out = []
    for start in range(0, len(images), batch_size):
        batch = torch.stack([transform(image) for image in images[start:start + batch_size]]).to(device)
        with torch.inference_mode(), torch.autocast("cuda", dtype=torch.float16, enabled=device.startswith("cuda")):
            values = model(batch)
        if isinstance(values, (tuple, list)):
            values = values[0]
        out.append(values.detach().float().cpu().numpy())
    return np.concatenate(out, axis=0).astype(np.float32) if out else np.zeros((0, 1536), np.float32)


def fetch(repo_file: str, slide_dir: Path) -> Path:
    from huggingface_hub import hf_hub_download

    return Path(hf_hub_download("mj0jpgg/tcga_cptac_cgga", repo_file, repo_type="dataset",
                                local_dir=str(slide_dir)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, help="CSV with file,size,patient,tss,cancer,slide")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--patches-per-slide", type=int, default=120)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--prefetch", type=int, default=3)
    args = parser.parse_args()

    output = Path(args.output_dir)
    staging = output / "staging"
    staging.mkdir(parents=True, exist_ok=True)
    slide_dir = output / "slides"
    slide_dir.mkdir(parents=True, exist_ok=True)

    catalog = pd.read_csv(args.catalog)
    pending = [row for row in catalog.itertuples()
               if not (staging / f"{row.slide}.npz").exists()]
    print(f"[plan] {len(catalog)} selected slides, {len(pending)} still to do", flush=True)
    if not pending:
        return

    import openslide

    model, transform, torch = load_encoder(args.device)
    pool = ThreadPoolExecutor(max_workers=args.prefetch)
    futures = {row.slide: pool.submit(fetch, row.file, slide_dir) for row in pending[:args.prefetch]}
    queue = list(pending)

    for position, row in enumerate(queue):
        try:
            path = futures.pop(row.slide).result()
        except Exception as exc:  # noqa: BLE001 - a broken slide must not kill the sweep
            print(f"[skip] {row.slide} download failed: {type(exc).__name__} {exc}", flush=True)
            (staging / f"{row.slide}.failed.json").write_text(json.dumps({"stage": "download", "error": str(exc)}))
            continue
        finally:
            nxt = position + args.prefetch
            if nxt < len(queue):
                futures[queue[nxt].slide] = pool.submit(fetch, queue[nxt].file, slide_dir)
        try:
            slide = openslide.OpenSlide(str(path))
            mpp = slide_mpp(slide)
            crop_px = int(round(FOV_MICRONS / mpp))
            level = min(len(slide.level_dimensions) - 1, 2)
            mask = tissue_mask(slide, level)
            mask_ds = float(slide.level_downsamples[level])
            rng = np.random.default_rng(args.seed + abs(hash(row.slide)) % 100000)
            positions = sample_patch_grid(slide, crop_px, mask, mask_ds, args.patches_per_slide, rng)
            if len(positions) < 10:
                raise ValueError(f"only {len(positions)} tissue patches found")
            images = [render_patch(slide, x, y, crop_px) for x, y in positions]
            vectors = embed(model, transform, torch, images, args.device, args.batch_size)
            if vectors.shape != (len(positions), 1536) or not np.isfinite(vectors).all():
                raise ValueError(f"bad embedding block {vectors.shape}")
            np.savez_compressed(
                staging / f"{row.slide}.npz", embeddings=vectors,
                positions=np.asarray(positions, dtype=np.int64),
                meta=np.asarray(json.dumps({
                    "slide": row.slide, "patient": row.patient, "cancer": str(row.cancer),
                    "tss": str(row.tss), "mpp": mpp, "crop_px": crop_px,
                    "fov_microns": FOV_MICRONS, "output_px": OUTPUT_PX,
                    "jpeg_quality": JPEG_QUALITY, "jpeg_subsampling": "4:2:0",
                    "tissue_fraction": TISSUE_FRACTION, "n_patches": len(positions),
                    "level0_dimensions": list(slide.dimensions),
                })))
            slide.close()
            print(f"[ok] {position + 1}/{len(queue)} {row.slide} {row.cancer} mpp={mpp:.4f} "
                  f"crop={crop_px} n={len(positions)}", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[skip] {row.slide} {type(exc).__name__} {exc}", flush=True)
            (staging / f"{row.slide}.failed.json").write_text(json.dumps({"stage": "patch", "error": str(exc)}))
        finally:
            try:
                os.remove(path)
            except OSError:
                pass
    pool.shutdown(wait=False)
    shutil.rmtree(slide_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
