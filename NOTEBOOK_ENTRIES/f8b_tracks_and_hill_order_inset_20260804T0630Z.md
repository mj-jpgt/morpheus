# 2026-08-04 06:30 UTC — Both `NEEDS EXTRACTION` items discharged, and the extraction found a third statistic substitution

**Logged:** 2026-08-04 06:30 UTC. **How obtained:** CPU only, thread-capped
(`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`); the A100 was at
99% on another chain and was read from but never computed on. Sources pulled from
`ubuntu@150.136.45.194` and vendored byte-identical (sha256 verified in both directions).

Two items were marked `NEEDS EXTRACTION` in `paper/P2_FIGURES.md`: per-step arrays for F8(b), and a
synthetic inset for F5. Both are now closed. The first one closed in two opposite directions, and on
the way it turned up something that belongs in §4.5(a)'s family.

---

## 1. Bad news first — `12.88 → 1.00` is **R3**, and three documents call it something else

`~/ws_d1/diag_d.py`, the A100 scratchpad script that produced F8(b)'s rank track, computes its
`eff-rank` column at lines 50–51 as:

```python
sv = torch.linalg.svdvals(wn.float() - wn.float().mean(0))
erank = float((sv.sum() ** 2) / (sv ** 2).sum())
```

with `wn = F.normalize(w, dim=-1)`. Rows L2-normalised, matrix column-centred, **order 2**. That is

```python
RANK_VARIANTS["R3"] == RankVariant(centre=True, normalise_rows=True, order=2)
```

It is not `spectral.CANONICAL` (order 1, rows at own norms). Yet:

| document | how it describes the number |
|---|---|
| `paper/P2_FIGURES.md` F8(a) | *"diagnostic-script centred effective rank"* — and labels the two rows beneath it **R3**, implying this one is something else |
| `paper/P2_FIGURES.md` F8(b) | *"the centred effective rank of the same objective"* |
| `paper/P1_CALIBRA_DRAFT.md` §4.11 | *"the centred effective rank of the same objective falls 12.88 → 1.00"* |
| `paper/P1_FIGURES.md` (F11 note) | same wording |
| `paper/P2_RANK_DRAFT.md` §4.9 table | *"centred eff-rank (diagnostic script)"* — under-specified, but does not assert order 1 |

**Same failure mode as §4.5(a), same tell: an inline formula where an import should have been.** It is
the third instance this project has caught, and like the first two it was invisible to review — the
label was plausible and the number was plausible. The same inline formula appears in
`diag_e/f/g/h.py` and `geom_probe.py` on the box, so **every** `eff-rank` column under
`~/e0_run/d1_diag/` is R3. `~/ws_d1/momentum_test.py:102-104`, written later, does import
`effective_rank` and reports R3 beside the canonical value; the `d1_diag` family predates that fix.

**Damage assessment — small, and it does not move a claim.** The instance is used only in §4.9/§4.10,
to argue that rank *does* catch total collapse. At the collapsed endpoint the label is nearly
irrelevant: a rank-1 matrix scores 1.00 under every variant in `RANK_VARIANTS`, so "R3 fell to 1.00"
and "R1 fell to ~1" say the same thing about the collapse. The starting value 12.88 is the part that
is unconstrained — row normalisation changes the spectrum, so R1 at that point is not recoverable by
inequality, and the diagnostic kept **no checkpoint and wrote no state**, so it cannot be recomputed
without re-running the training on a GPU. It belongs in the `[NOT RECOMPUTED — needs a GPU re-run]`
category the drafts already use.

**What I changed and what I did not.** `paper/P2_FIGURES.md` is corrected (F8(a) row label, F8(b)
prose, the F8 status block, and the P1-corrections list now names the two P1 files that need the same
fix). `paper/P2_RANK_DRAFT.md` and `NOTEBOOK.md` are **not** edited, per the standing instruction;
§4.9's row is under-specified rather than wrong, but should be relabelled by whoever owns the draft.
`P1_CALIBRA_DRAFT.md` and `P1_FIGURES.md` are also left alone and flagged instead.

The source script is vendored verbatim at `v2/research/rebase/p2/collapse_tracks/diag_d.py`
(sha256 `7448d7c5…`) and allowlisted in `v2/tests/test_effective_rank_canonical.py` with the reason:
it must **not** be rewritten to call `spectral.py`, because keeping the inline formula is what makes
the identification checkable instead of a claim in a notebook entry.

---

## 2. F8(b) — the extraction resolved in two opposite directions

`paper/P2_FIGURES.md` asked for per-step arrays, *"else before/after paired markers"*. The answer is
one of each.

**The rank track HAS per-step values.** `~/e0_run/d1_diag/diag_d.log` (sha256 `d2b49035…`) carries six
recorded steps of the clean in-batch InfoNCE arm — the arm with every other loss term zeroed:

| step | rank (R3) | acc@1 | pos cos | worst-neg cos | min margin | wsi-wsi cos |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | **12.88** | 0.188 | 0.1235 | 0.2186 | −0.2190 | 0.8008 |
| 25 | 1.21 | 0.125 | 0.9344 | 0.9872 | −0.1103 | 0.8702 |
| 50 | **1.00** | 0.188 | 0.9993 | 0.9993 | −0.0001 | 0.6119 |
| 100 | 1.00 | 0.125 | 0.9995 | 0.9995 | −0.0001 | 0.6103 |
| 200 | 1.00 | 0.062 | 0.9996 | 0.9996 | −0.0001 | 0.6165 |
| 400 | 1.00 | 0.125 | 0.9999 | 0.9999 | −0.0000 | 0.5904 |

That half of F8(b) can now be a curve rather than two markers.

**The collapse-evidence quantities have ONLY endpoints, and that is now confirmed rather than
assumed.** `~/e0_run/collapse_diag.log` (sha256 `bedadc23…`) is the sole surviving source of the
16/16 instance, and its script probes before and after training, not on a schedule. Arm A, every
value F8(b) quotes, reproduced exactly:

| quantity | before → after |
|---|---|
| in-batch InfoNCE (chance 2.7726) | 3.4681 → 2.7734 |
| retrieval acc@1 (chance 0.062) | 0.062 → **0.000** |
| cross-modal positive cosine | 0.0538 → 0.9959 |
| cross-modal negative cosine | 0.0816 → **0.9960** |
| WSI within-modality off-diagonal cosine | 0.7089 → 0.9999 |
| `z_biology` **hard** matrix rank | 16 → 16 (max 16) |

The negatives finishing marginally *above* the positives is preserved and asserted by the test — that
inversion is what makes the instance evidence rather than an anomaly.

**A discrepancy against `NOTEBOOK.md`, reported not worked around.** `NOTEBOOK.md:1327-1328`
tabulates arms B and C of this same experiment with three cells that the surviving log does not
support:

| cell | `NOTEBOOK.md` | `collapse_diag.log` |
|---|---:|---:|
| B, cross-modal negative cosine | 0.4922 | **0.4895** |
| B, WSI within-modality cosine | 0.4946 | **0.4918** |
| C, cross-modal positive cosine | 0.9998 | **0.9997** |

Arm A — the only arm F8(b) draws — matches to every digit. A search of `~/e0_run` finds 0.4922 and
0.4946 nowhere else (the 0.4922 in `g26_variance_floor_fix_20260803T0210Z.md` is arm D of the
*variance-floor-fix* run, a different experiment), so the most likely explanation is a second,
unretained run of the same nondeterministic diagnostic. Either way the notebook's B/C cells are not
reproducible from anything now on disk, and should not be quoted at four digits. `NOTEBOOK.md` is not
edited here.

**Vendored, with a parser and a test.** `v2/research/rebase/p2/collapse_tracks/` holds both logs and
the script; `v2/research/rebase/p2/p2_f8b_tracks.py` parses them and refuses to report if a log's
sha256 has drifted; `v2/tests/test_p2_f8b_tracks.py` pins the arrays, the endpoints, the 16/16
pinning, the negatives-above-positives inversion, and the R3 identification.

One parser detail worth recording because it nearly ate the panel's own subject: the log aligns its
columns with a variable number of spaces, and two of the six quantities — including the `16 -> 16`
rank pinning that F8(b) is *about* — use a single space. A fixed-width separator parsed four of six
quantities and dropped that one silently.

---

## 3. F5 inset — the size of the §4.5(a) substitution, on a synthetic family

`v2/research/rebase/p2/p2_hill_order_inset.py`, deterministic (seed 20260804), CPU-only, with
`v2/tests/test_p2_hill_order_inset.py`.

Power-law spectra `σ_k = k^−a`, 64 components on 256 rows, realised **exactly**: the matrix is an
orthonormal basis of the mean-zero subspace scaled by the requested spectrum, so centring is a no-op
and the centred singular values are the ones asked for (asserted in the test through the Gram
eigenvalues, an independent route).

| `a` | R1 | R2 `(Σσ)²/Σσ²` | PR `(Σσ²)²/Σσ⁴` | R2/PR | R1/PR |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 64.000 | 64.000 | 64.000 | 1.000 | 1.000 |
| 0.2 | 62.788 | 61.394 | 51.995 | 1.181 | 1.208 |
| 0.5 | 54.990 | 44.946 | 13.811 | 3.254 | 3.982 |
| 0.9 | 34.517 | 17.922 | 3.026 | **5.923** | 11.407 |
| 1.0 | 29.118 | 13.811 | 2.453 | 5.630 | **11.870** |
| 1.5 | 10.817 | 4.647 | 1.420 | 3.272 | 7.618 |
| 2.0 | 4.668 | 2.453 | 1.167 | 2.103 | 4.001 |

**R2/PR spans 1.000 to 5.923** over `a ∈ [0, 2]`; R1/PR reaches **11.87**. The ratio is non-monotone —
it peaks near `a ≈ 0.9` and falls again as the spectrum concentrates towards rank 1, where all three
statistics converge on 1. The two agree **only** on a flat spectrum.

**The hand anchor**, the same one `test_effective_rank_canonical.py` uses: `σ ∝ (2,1,1)` gives
**R1 = 2√2 = 2.8284271247**, **R2 = 8/3 = 2.6666666667**, **PR = 2.0000000000**, all three reproduced
to ten digits by the imported implementations. Three different numbers on one spectrum, which is the
whole argument of §3.1 in one line.

**A closed-form fingerprint fell out of the family and is now asserted.** On a power law, `PR` at
decay `a` equals `R2` at decay `2a` **exactly** — visible in the table (`PR(0.5) = 13.811 = R2(1.0)`,
`PR(0.2) = 51.995 = R2(0.4)`) — because PR is the order-2 Hill number of the *squared* spectrum and
squaring a power law doubles its exponent. It is a check a referee can run in one line, and it would
break immediately if either statistic were swapped for the other.

**No statistic is implemented in the script.** R1 and R2 come from `spectral.effective_rank` under the
named `RANK_VARIANTS`; PR comes from `p2_competing_metrics.participation_ratio`. The test asserts the
script contains no decomposition of its own — the discipline whose absence produced the error the
inset is about.

---

## Files / commits

New: `v2/research/rebase/p2/p2_f8b_tracks.py`, `v2/research/rebase/p2/p2_hill_order_inset.py`,
`v2/research/rebase/p2/collapse_tracks/{README.md,diag_d.log,collapse_diag.log,diag_d.py}`,
`v2/tests/test_p2_f8b_tracks.py`, `v2/tests/test_p2_hill_order_inset.py`.
Edited: `paper/P2_FIGURES.md` (F5 inset spec, F8 status and statistic labels, S6's two corrections,
the two pending-dependency rows, the would-be-figure row for E1),
`v2/tests/test_effective_rank_canonical.py` (one allowlist entry, with its reason).
Not edited: `NOTEBOOK.md`, `paper/P2_RANK_DRAFT.md`, `paper/P1_*`, anything under `v2/calibra/`.

**Suite: 277 passed in `v2/tests` (266 before; +11), thread-capped.** One pre-existing failure,
`test_paper_paths_resolve.py::test_no_box_output_basename_is_actually_in_the_repository`, is caused by
an untracked `v2/research/rebase/p2/figures/data/` tree belonging to concurrent work in this working
copy; it is not mine and is not touched.
