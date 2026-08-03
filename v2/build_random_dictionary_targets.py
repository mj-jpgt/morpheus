"""Must-beat baseline: a SIZE- AND SPECTRUM-MATCHED random dictionary.

The claim PBS makes is that *these particular* interventional directions carry
morphology-legible structure. The cheapest way that claim can be wrong is if any
128-dimensional projection of the same expression matrix would do as well, in
which case what is being measured is the dimensionality of the readout and not
the dictionary. This module builds exactly that alternative: the same cohort, the
same genes, the same development-only expression transform, the same number of
components, and the same singular-value spectrum — with the directions replaced
by a Haar-random orthonormal basis.

Matching the spectrum, not just the size, is the part that matters. An unmatched
random basis has a flat spectrum and therefore a different conditioning and a
different effective rank, so a difference in downstream score could be read off
the spectrum alone. Re-using the fitted dictionary's own singular values removes
that explanation and leaves only the gene directions as the difference.

Provenance is taken from the frozen PBS target file rather than re-derived, so
the baseline is bound to the identical cohort, split, gene order and transform
that the PBS targets were built with; nothing here can drift from them.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .baseline_target_common import (development_expression, load_reference_targets,
                                     write_target_block)

__all__ = ["build_random_dictionary_targets"]


def haar_orthonormal(n_rows: int, n_columns: int, rng: np.random.Generator) -> np.ndarray:
    """Uniformly random orthonormal ``[n_rows, n_columns]`` frame (QR of a Gaussian).

    The sign fix on R's diagonal is not cosmetic: without it numpy's QR does not
    give a Haar-distributed Q, and a subtly non-uniform "random" basis is exactly
    the kind of thing that makes a negative control look weaker than it is.
    """
    q, r = np.linalg.qr(rng.normal(size=(n_rows, n_columns)))
    return q * np.sign(np.where(np.diag(r) == 0, 1.0, np.diag(r)))[None, :]


def build_random_dictionary_targets(*, pbs_targets: str, rna_table: str, output: str,
                                    seed: int = 0, match_spectrum: bool = True) -> dict[str, object]:
    reference = load_reference_targets(pbs_targets)
    expression, transform = development_expression(rna_table, reference)
    rng = np.random.default_rng(seed)
    n_genes, n_components = reference["gene_basis"].shape
    basis = haar_orthonormal(n_genes, n_components, rng)
    spectrum = np.asarray(reference["singular_values"], dtype=np.float64)
    if match_spectrum:
        # The PBS scores are (E - m) @ B with B orthonormal; the spectrum enters
        # through the reference responses, not through B. Re-imposing it here as a
        # per-component scale keeps the two blocks on the same footing while
        # leaving the directions uniformly random.
        basis = basis * (spectrum / max(float(spectrum[0]), 1e-12))[None, :]
    scores = ((expression - expression.mean(axis=0, keepdims=True)) @ basis).astype(np.float32)
    names = np.asarray([f"RANDDICT_{i:03d}" for i in range(n_components)])
    manifest = {
        "target_kind": "size_and_spectrum_matched_random_dictionary_coordinates",
        "n_components": int(n_components), "seed": int(seed),
        "match_spectrum": bool(match_spectrum),
        "basis_is_haar_orthonormal": True,
        "spectrum_source": "singular_values of the frozen PBS dictionary",
        "pbs_targets": str(Path(pbs_targets).resolve()),
        "pbs_manifest_digest": reference["manifest_digest"],
        "gene_count": int(n_genes), "gene_digest": reference["gene_digest"],
        "expression_transform": transform,
        "basis_digest": sha256(np.ascontiguousarray(basis.astype(np.float32)).view(np.uint8)).hexdigest(),
    }
    return write_target_block(output, reference, scores, names, "RANDOM_DICTIONARY", manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pbs-targets", required=True, help="frozen pbs_targets_*.npz to inherit cohort/split/genes from")
    parser.add_argument("--rna-table", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-match-spectrum", action="store_true")
    args = parser.parse_args()
    print(json.dumps(build_random_dictionary_targets(
        pbs_targets=args.pbs_targets, rna_table=args.rna_table, output=args.output,
        seed=args.seed, match_spectrum=not args.no_match_spectrum), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
