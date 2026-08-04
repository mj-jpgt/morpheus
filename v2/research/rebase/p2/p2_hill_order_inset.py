"""F5 inset: how far apart `(sum s)^2 / sum s^2` and `(sum s^2)^2 / sum s^4` actually are.

`paper/P2_FIGURES.md` F5 asks for this and marks it `NEEDS EXTRACTION`:

    A small inset showing `(Ss)^2/Ss^2` against `(Ss^2)^2/Ss^4` on one synthetic spectrum family
    would make the size of the substitution concrete rather than asserted. That is NEEDS EXTRACTION
    -- cheap, CPU-only, deterministic new computation, which must be written into
    `v2/research/rebase/p2/` with a test rather than a scratchpad file.

Those are the two statistics whose confusion produced the SS4.5(a) correction: the draft's second and
third rows were labelled `R2`/`R3` but were computed as the order-2 Hill number of the EIGENVALUE
distribution, and are now relabelled `PR`/`PR_rownorm`. The correction moved a published count from
3 of 6 to 2 of 6.

**No rank is computed here.** `R1` and `R2` come from `spectral.effective_rank` under the named
`RANK_VARIANTS`; `PR` comes from `p2_competing_metrics.participation_ratio`. This module only builds
matrices with a prescribed spectrum and tabulates what those two imported functions return -- which is
the discipline whose absence produced the error the inset is about.

**How the matrices are built.** Each spectrum `s` is realised exactly, as `X = Q diag(s)` where `Q`
has orthonormal columns lying in the mean-zero subspace. Centring is then a no-op and the singular
values of the centred matrix are `s` to machine precision, so the closed forms are checkable against
the imported functions rather than approximated.

Run::

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
      python3 -m morpheus.v2.research.rebase.p2.p2_hill_order_inset --output P2_HILL_ORDER_INSET.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2.p2_competing_metrics import participation_ratio

#: Rows and spectrum length of the synthetic family. 64 components on 256 rows is the same
#: rows >> components regime the frozen artifacts sit in (2,766 x 256), so the inset is not
#: demonstrating a small-sample artefact.
N_ROWS = 256
N_COMPONENTS = 64

#: Power-law spectra `s_k = k^-a`. `a = 0` is flat (every statistic equals the component count and
#: the two agree exactly); increasing `a` concentrates the spectrum, which is the regime a real
#: representation lives in and the regime where the two statistics separate.
DECAYS = tuple(round(0.1 * i, 1) for i in range(21))

#: The hand-checkable anchor already used by `v2/tests/test_effective_rank_canonical.py`:
#: `s ~ (2, 1, 1)` gives R1 = 2*sqrt(2), R2 = 8/3 and PR = 2 exactly.
ANCHOR = (2.0, 1.0, 1.0)

SEED = 20260804


def centred_basis(n_rows: int, n_components: int, seed: int = SEED) -> np.ndarray:
    """Orthonormal columns spanning a subspace of the mean-zero hyperplane.

    Centring first and orthonormalising second is what makes the construction exact: a basis
    orthonormalised before centring has a non-zero column mean, and the centred spectrum is then
    not the one that was asked for.
    """
    raw = np.random.default_rng(seed).normal(size=(n_rows, n_components))
    basis = np.linalg.qr(raw - raw.mean(axis=0, keepdims=True))[0]
    return basis[:, :n_components]


def matrix_with_spectrum(singular_values, *, n_rows: int = N_ROWS, seed: int = SEED) -> np.ndarray:
    """A matrix whose CENTRED singular values are exactly `singular_values`."""
    values = np.asarray(singular_values, dtype=np.float64)
    return centred_basis(n_rows, len(values), seed) * values


def power_law_spectrum(decay: float, n_components: int = N_COMPONENTS) -> np.ndarray:
    return np.power(np.arange(1, n_components + 1, dtype=np.float64), -float(decay))


def score(singular_values, *, n_rows: int = N_ROWS, seed: int = SEED) -> dict[str, float]:
    """R1, R2 and PR on one spectrum, each from its own imported implementation."""
    x = matrix_with_spectrum(singular_values, n_rows=n_rows, seed=seed)
    r1 = effective_rank(x, variant=RANK_VARIANTS["R1"])
    r2 = effective_rank(x, variant=RANK_VARIANTS["R2"])
    pr = participation_ratio(x)
    return {"R1": r1, "R2": r2, "PR": pr,
            "R2_over_PR": r2 / pr if pr > 0 else float("nan"),
            "R1_over_PR": r1 / pr if pr > 0 else float("nan")}


def build(decays=DECAYS, *, n_rows: int = N_ROWS, n_components: int = N_COMPONENTS,
          seed: int = SEED) -> dict:
    rows = []
    for decay in decays:
        record = {"decay": float(decay)}
        record.update(score(power_law_spectrum(decay, n_components), n_rows=n_rows, seed=seed))
        rows.append(record)
    ratios = [row["R2_over_PR"] for row in rows]
    return {
        "figure": "P2 F5 inset",
        "question": "how large is the (sum s)^2/sum s^2 versus (sum s^2)^2/sum s^4 substitution",
        "statistics": {
            "R1": "spectral.CANONICAL -- Roy & Vetterli order 1, column-centred, rows at own norms",
            "R2": "RANK_VARIANTS['R2'] -- (sum s)^2 / sum s^2, order 2 on the SINGULAR VALUES",
            "PR": "p2_competing_metrics.participation_ratio -- (sum s^2)^2 / sum s^4, order 2 on the EIGENVALUES",
        },
        "block": "synthetic; the matrix is constructed to have exactly the stated centred spectrum",
        "construction": {"n_rows": n_rows, "n_components": n_components, "seed": seed,
                         "family": "power law s_k = k^-a"},
        "anchor": {"spectrum": list(ANCHOR), **score(ANCHOR, n_rows=n_rows, seed=seed),
                   "closed_form": {"R1": 2 * float(np.sqrt(2)), "R2": 8 / 3, "PR": 2.0}},
        "rows": rows,
        # A closed-form fingerprint that makes the table self-checking: PR is the order-2 Hill number
        # of the SQUARED spectrum, and squaring a power law doubles its exponent, so on this family
        # PR at decay `a` equals R2 at decay `2a` exactly. It is visible in the printed table
        # (PR(0.5) = 13.811 = R2(1.0)) and asserted in `v2/tests/test_p2_hill_order_inset.py`.
        "identity": "PR(s_k = k^-a) == R2(s_k = k^-2a), exactly",
        "summary": {
            "R2_over_PR_min": float(min(ratios)),
            "R2_over_PR_max": float(max(ratios)),
            "R2_over_PR_at_flat_spectrum": float(rows[0]["R2_over_PR"]),
            "R2_over_PR_at_steepest": float(rows[-1]["R2_over_PR"]),
        },
    }


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--n-rows", type=int, default=N_ROWS)
    parser.add_argument("--n-components", type=int, default=N_COMPONENTS)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)
    payload = build(n_rows=args.n_rows, n_components=args.n_components, seed=args.seed)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    anchor = payload["anchor"]
    print(f"anchor s ~ (2,1,1):  R1 {anchor['R1']:.10f}  R2 {anchor['R2']:.10f}  PR {anchor['PR']:.10f}")
    print(f"{'decay a':>8}{'R1':>10}{'R2':>10}{'PR':>10}{'R2/PR':>9}{'R1/PR':>9}")
    for row in payload["rows"]:
        print(f"{row['decay']:>8.1f}{row['R1']:>10.3f}{row['R2']:>10.3f}{row['PR']:>10.3f}"
              f"{row['R2_over_PR']:>9.3f}{row['R1_over_PR']:>9.3f}")
    summary = payload["summary"]
    print(f"R2/PR spans {summary['R2_over_PR_min']:.3f} to {summary['R2_over_PR_max']:.3f} "
          f"over a in [{DECAYS[0]}, {DECAYS[-1]}]")
    return payload


if __name__ == "__main__":
    main()
