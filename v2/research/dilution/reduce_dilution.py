"""Reduce every dilution level to 256 dimensions, identically, before CALIBRA.

WHY. The dilution artifact is `concat(mean, std)` over 1536-d H-Optimus tokens, so
each state is 3,072-wide. CALIBRA's permutation null re-whitens the representation
block once per permutation, and an SVD of a 2,766 x 3,072 matrix a thousand times
over, for each of seven levels, does not finish in the time available. Reducing to
256 makes the run affordable AND makes the dilution arm capacity-matched to the D2
states it is being read beside, which is the fairer comparison anyway.

WHAT IS AND IS NOT HELD FIXED. The reduction is refit **per level, on train rows
only**. Refitting is the honest choice: each dilution level is a different
representation, and freezing the level-0 basis would measure "how far the diluted
bags drift off the level-0 axes" rather than "how much of the channel survives" --
it would build the answer into the transform. Fitting on train rows only keeps the
held-out cancers out of the basis. The component count, the fit population, the
seed and the whitening convention are identical across levels, so nothing but the
dilution changes.

The number of components is capped at 256 and reported per level, along with the
retained variance, so a level whose spectrum collapsed cannot pass unnoticed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-components", type=int, default=256)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    from sklearn.decomposition import PCA

    raw = np.load(args.artifact, allow_pickle=True)
    manifest = json.loads(str(raw["manifest_json"]))
    split = np.asarray([str(s) for s in raw["split"]])
    train = split == "train"
    states = [str(s) for s in raw["trained_states"]]
    reduced: dict[str, np.ndarray] = {}
    detail = {}
    for state in states:
        values = np.asarray(raw[state], dtype=np.float64)
        width = int(min(args.n_components, train.sum() - 1, values.shape[1]))
        model = PCA(n_components=width, random_state=args.seed).fit(values[train])
        scores = model.transform(values).astype(np.float32)
        if not np.isfinite(scores).all() or float(scores.std(axis=0).min()) < 1e-8:
            raise RuntimeError(f"reduced state {state} is degenerate")
        reduced[state] = scores
        detail[state] = {"n_components": width,
                         "explained_variance_ratio_sum": float(model.explained_variance_ratio_.sum()),
                         "explained_variance_ratio_first": float(model.explained_variance_ratio_[0]),
                         "source_width": int(values.shape[1])}
        print(f"[reduce] {state}: {values.shape[1]} -> {width}, "
              f"variance retained {detail[state]['explained_variance_ratio_sum']:.4f}, "
              f"PC1 {detail[state]['explained_variance_ratio_first']:.4f}", flush=True)

    manifest = dict(manifest)
    manifest["method"] = str(manifest.get("method", "dilution")) + "_pca256"
    manifest["reduction"] = {"kind": "pca", "n_components": int(args.n_components),
                             "fit_population": "train_rows_only", "refit_per_level": True,
                             "seed": int(args.seed), "per_state": detail,
                             "rationale": "capacity-match to the 256-d D2 states and make the "
                                          "permutation null affordable; refit per level so the "
                                          "transform cannot encode the answer"}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, patient_ids=raw["patient_ids"], cancers=raw["cancers"], split=raw["split"],
                        trained_states=np.asarray(states, dtype=str),
                        artifact_version=raw["artifact_version"],
                        manifest_json=np.asarray(json.dumps(manifest, sort_keys=True)), **reduced)
    output.with_suffix(output.suffix + ".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[done] -> {output}", flush=True)


if __name__ == "__main__":
    main()
