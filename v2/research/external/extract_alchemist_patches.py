"""Patch + embed the ALCHEMIST external cohort to the *identical* TCGA patch spec.

Everything that touches a pixel is imported from
``v2/research/dilution/extract_normal_patches.py`` rather than reimplemented.
That module's constants are already asserted against the spatial adapter's by a
test, because the two drifted once; importing keeps ALCHEMIST inside that guard
instead of opening a third place for the spec to drift.

What is genuinely different here, and only this:

  * the slide arrives from the GDC data endpoint rather than the HF hub;
  * there is no polygon, so patches are drawn from the whole-slide tissue mask.
    This is the declared, irreducible deviation: TCGA patches are restricted to
    pathologist-drawn tumour polygons and no WSI tumour polygons exist for
    ALCHEMIST -- nor for any non-TCGA GDC project.
  * Aperio header fields are recorded per slide.  ALCHEMIST publishes no
    ``tissue_source_site``, so the scanner identity is the only batch-like
    variable that exists at all.  It is recorded, not silently used as a site.

Sampling rule: ``--patches-per-slide 30``, matched to the TCGA store, whose
per-slide count is 30 at the median and at both the 10th and 90th percentiles.

Resumable: one ``.npz`` per slide under ``--output-dir/staging``; the ``.svs`` is
deleted the moment its patches are embedded, so peak disk stays bounded at
roughly ``prefetch`` slides rather than the 1.7 TB the cohort would otherwise
occupy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3].parent))

from morpheus.v2.research.dilution.extract_normal_patches import (  # noqa: E402
    FOV_MICRONS,
    JPEG_QUALITY,
    JPEG_SUBSAMPLING,
    OUTPUT_PX,
    TISSUE_FRACTION,
    embed,
    load_encoder,
    render_patch,
    sample_patch_grid,
    slide_mpp,
    tissue_mask,
)

GDC_DATA = "https://api.gdc.cancer.gov/data"
# Header keys worth keeping: the scanner identity is the only batch-like variable
# ALCHEMIST exposes, and the magnification/mpp pair is what makes the fixed-micron
# crop magnification-invariant.
APERIO_KEYS = ("aperio.ScanScope ID", "aperio.AppMag", "aperio.MPP", "aperio.Date",
               "aperio.ScanScopeID", "openslide.vendor", "openslide.objective-power",
               "openslide.mpp-x", "openslide.mpp-y")


def download(file_id: str, destination: Path, expected_md5: str | None = None,
             retries: int = 4) -> Path:
    """Stream one open-access GDC file to disk, verifying its published md5."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        digest = hashlib.md5()
        try:
            request = urllib.request.Request(f"{GDC_DATA}/{file_id}")
            with urllib.request.urlopen(request, timeout=600) as response, \
                    open(destination, "wb") as handle:
                while True:
                    chunk = response.read(1 << 22)
                    if not chunk:
                        break
                    digest.update(chunk)
                    handle.write(chunk)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
            continue
        if expected_md5 and digest.hexdigest() != expected_md5:
            if attempt == retries - 1:
                raise ValueError(f"md5 mismatch for {file_id}: {digest.hexdigest()} != {expected_md5}")
            time.sleep(2 ** attempt)
            continue
        return destination
    raise RuntimeError("unreachable")


def slide_properties(slide) -> dict:
    return {key: str(slide.properties[key]) for key in APERIO_KEYS if key in slide.properties}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True,
                        help="CSV from build_alchemist_manifest.py: file,size,patient,tss,cancer,slide,case_id,md5sum")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--patches-per-slide", type=int, default=30,
                        help="30 matches the TCGA store's per-slide median and its 10th/90th percentiles")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--prefetch", type=int, default=6,
                        help="concurrent GDC downloads; 6 saturates the link at ~170 MB/s")
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--num-shards", type=int, default=1)
    args = parser.parse_args()

    output = Path(args.output_dir)
    staging = output / "staging"
    staging.mkdir(parents=True, exist_ok=True)
    slide_dir = output / "slides"
    slide_dir.mkdir(parents=True, exist_ok=True)

    catalog = pd.read_csv(args.catalog)
    if args.num_shards > 1:
        catalog = catalog.iloc[args.shard::args.num_shards]
    pending = [row for row in catalog.itertuples()
               if not (staging / f"{row.slide}.npz").exists()]
    print(f"[plan] {len(catalog)} slides in shard {args.shard}/{args.num_shards}, "
          f"{len(pending)} still to do, "
          f"{sum(int(r.size) for r in pending) / 1e12:.3f} TB to fetch", flush=True)
    if not pending:
        return

    import openslide

    model, transform, torch = load_encoder(args.device)
    pool = ThreadPoolExecutor(max_workers=args.prefetch)

    def submit(row):
        return pool.submit(download, row.file, slide_dir / f"{row.slide}.svs", str(row.md5sum))

    queue = list(pending)
    futures = {row.slide: submit(row) for row in queue[:args.prefetch]}
    started = time.time()

    for position, row in enumerate(queue):
        try:
            path = futures.pop(row.slide).result()
        except Exception as exc:  # noqa: BLE001 - one broken slide must not kill the cohort
            print(f"[skip] {row.slide} download failed: {type(exc).__name__} {exc}", flush=True)
            (staging / f"{row.slide}.failed.json").write_text(
                json.dumps({"stage": "download", "error": str(exc)}))
            continue
        finally:
            nxt = position + args.prefetch
            if nxt < len(queue):
                futures[queue[nxt].slide] = submit(queue[nxt])
        try:
            slide = openslide.OpenSlide(str(path))
            mpp = slide_mpp(slide)
            crop_px = int(round(FOV_MICRONS / mpp))
            level = min(len(slide.level_dimensions) - 1, 2)
            mask = tissue_mask(slide, level)
            mask_ds = float(slide.level_downsamples[level])
            rng = np.random.default_rng(args.seed + int(hashlib.sha256(
                str(row.slide).encode()).hexdigest()[:8], 16) % 100000)
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
                    "cohort": "ALCHEMIST-ALCH",
                    "slide": row.slide, "patient": row.patient, "case_id": row.case_id,
                    "cancer": str(row.cancer), "tss": str(row.tss),
                    "gdc_file_id": row.file, "md5sum": str(row.md5sum),
                    "mpp": mpp, "crop_px": crop_px,
                    "fov_microns": FOV_MICRONS, "output_px": OUTPUT_PX,
                    "jpeg_quality": JPEG_QUALITY, "jpeg_subsampling": "4:2:0",
                    "tissue_fraction": TISSUE_FRACTION, "n_patches": len(positions),
                    "sampling": "whole-slide tissue mask; NO tumour polygon exists for ALCHEMIST",
                    "level0_dimensions": list(slide.dimensions),
                    "slide_properties": slide_properties(slide),
                })))
            slide.close()
            rate = (position + 1) / max(time.time() - started, 1e-6) * 3600
            print(f"[ok] {position + 1}/{len(queue)} {row.slide} mpp={mpp:.4f} crop={crop_px} "
                  f"n={len(positions)} ({rate:.0f} slides/h)", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[skip] {row.slide} {type(exc).__name__} {exc}", flush=True)
            (staging / f"{row.slide}.failed.json").write_text(
                json.dumps({"stage": "patch", "error": str(exc)}))
        finally:
            try:
                os.remove(path)
            except OSError:
                pass
    pool.shutdown(wait=False)
    shutil.rmtree(slide_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
