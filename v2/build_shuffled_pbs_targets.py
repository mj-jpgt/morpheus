"""T1.5 -- build gene-label-shuffled PBS target blocks from the FROZEN dictionary.

WHY NOT ``build_pbs_targets --shuffle-gene-labels``. That flag exists and works,
but it refits the dictionary from the Perturb-seq reference, and doing so on this
machine now fails: the data config whose digest the frozen target file records
(``data_config_sha256 = 76927870...``) is no longer on disk, and every surviving
config declares a cohort that is missing 249 of the split's patients. Refitting
against a different config would silently produce a shuffled block bound to a
different cohort than the block it is being compared with, which is exactly the
comparison this control cannot afford to get wrong.

So the dictionary is not refit at all. It is REBOUND: the frozen
``gene_basis`` is read straight out of ``pbs_targets_k128_v2.npz`` and its ROWS
are permuted, which is precisely what the ``--shuffle-gene-labels`` path does
(``build_pbs_targets.py`` permutes ``gene_mean``/``gene_basis`` after the fit).
The directions, their spectrum and their mutual geometry are untouched; only the
gene each loading names has moved.

FIDELITY CHECK, RUN EVERY TIME. The frozen npz stores ``gene_basis`` but not
``gene_mean``, so the reconstruction here reproduces the stored scores only up to
a per-column additive constant. That constant is annihilated by every statistic
CALIBRA computes (all of them are centred), but "it should not matter" is not
evidence, so the builder recomputes the UNSHUFFLED scores from the same path and
refuses to write anything unless every column matches the frozen scores at
Pearson r >= 0.9999 after centring. A failure there means the expression
transform has drifted and the shuffled block would not be comparable.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .baseline_target_common import development_expression, load_reference_targets

__all__ = ["reconstruct_scores", "build_shuffled_pbs_targets"]


def reconstruct_scores(expression: np.ndarray, gene_basis: np.ndarray) -> np.ndarray:
    """``(E - gene_mean) @ B`` up to the per-column constant ``-gene_mean @ B``."""
    return np.asarray(expression, dtype=np.float64) @ np.asarray(gene_basis, dtype=np.float64)


def _column_agreement(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64) - np.mean(a, axis=0, keepdims=True)
    b = np.asarray(b, dtype=np.float64) - np.mean(b, axis=0, keepdims=True)
    denominator = np.linalg.norm(a, axis=0) * np.linalg.norm(b, axis=0)
    return np.where(denominator > 1e-30, np.sum(a * b, axis=0) / np.maximum(denominator, 1e-30), np.nan)


def build_shuffled_pbs_targets(*, pbs_targets: str, rna_table: str, output: str, seed: int,
                               fidelity_floor: float = 0.9999) -> dict:
    reference = load_reference_targets(pbs_targets)
    frozen = np.load(pbs_targets, allow_pickle=True)
    frozen_scores = np.asarray(frozen["scores"], dtype=np.float64)
    gene_basis = np.asarray(frozen["gene_basis"], dtype=np.float64)
    expression, transform = development_expression(rna_table, reference)

    reconstructed = reconstruct_scores(expression, gene_basis)
    agreement = _column_agreement(reconstructed, frozen_scores)
    worst = float(np.nanmin(agreement))
    if not np.isfinite(worst) or worst < fidelity_floor:
        raise RuntimeError(
            f"reconstruction of the UNSHUFFLED PBS scores agrees with the frozen block at only "
            f"r = {worst:.6f} (floor {fidelity_floor}); the expression transform has drifted and a "
            f"shuffled block built through it would not be comparable")

    order = np.random.default_rng(seed).permutation(len(gene_basis))
    fixed_points = int((order == np.arange(len(order))).sum())
    shuffled_scores = reconstruct_scores(expression, gene_basis[order]).astype(np.float32)
    if not np.isfinite(shuffled_scores).all() or float(shuffled_scores.std(axis=0).min()) < 1e-8:
        raise RuntimeError("shuffled PBS block is non-finite or has a constant column")

    names = np.asarray([f"PBSSHUF_{i:03d}" for i in range(shuffled_scores.shape[1])])
    manifest = {
        "schema_version": "1.0",
        "target_kind": "gene_label_shuffled_external_perturbation_dictionary_coordinates",
        "control": "T1.5_must_fail_shuffled_gene_labels",
        "construction": "frozen gene_basis rows permuted; directions, spectrum and geometry unchanged",
        "rebound_not_refit": True,
        "shuffle_seed": int(seed), "fixed_points": fixed_points,
        "permutation_digest": sha256(order.astype(np.int64).tobytes()).hexdigest(),
        "unshuffled_reconstruction_min_column_r": worst,
        "unshuffled_reconstruction_median_column_r": float(np.nanmedian(agreement)),
        "fidelity_floor": float(fidelity_floor),
        "gene_mean_offset_dropped": "scores differ from the frozen block by a per-column additive "
                                    "constant only; every CALIBRA statistic is centred",
        "pbs_targets": str(Path(pbs_targets).resolve()),
        "pbs_manifest_digest": reference["manifest_digest"],
        "gene_count": int(len(reference["genes"])), "gene_digest": reference["gene_digest"],
        "expression_transform": transform,
        "n_components": int(shuffled_scores.shape[1]),
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path, patient_ids=reference["patient_ids"], split=reference["split"], cancers=reference["cancers"],
        scores=shuffled_scores, target_names=names,
        target_groups=np.asarray(["PBS_GENE_LABEL_SHUFFLED"] * len(names)),
        genes=reference["genes"], gene_basis=gene_basis[order].astype(np.float32),
        singular_values=reference["singular_values"].astype(np.float32),
        manifest_json=np.asarray(json.dumps(manifest, sort_keys=True)))
    path.with_suffix(path.suffix + ".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pbs-targets", required=True)
    parser.add_argument("--rna-table", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--fidelity-floor", type=float, default=0.9999)
    args = parser.parse_args()
    print(json.dumps(build_shuffled_pbs_targets(
        pbs_targets=args.pbs_targets, rna_table=args.rna_table, output=args.output,
        seed=args.seed, fidelity_floor=args.fidelity_floor), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
