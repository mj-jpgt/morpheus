# F8(b) collapse tracks — the source logs, vendored

`paper/P2_FIGURES.md` F8 marks panel (b) **`NEEDS EXTRACTION`**: *"the collapse-evidence values are
recorded as endpoint pairs ... not per-step arrays. Either extract the per-step arrays from
`~/e0_run/d1_diag/`, or draw them as before/after paired markers labelled 'endpoint values as
recorded; per-step array not retained'. Do not interpolate."*

Both source logs and the script that wrote one of them are vendored here, byte-identical to the box,
so the answer is checkable rather than asserted. `p2_f8b_tracks.py` parses them into
`F8B_TRACKS.json`; `v2/tests/test_p2_f8b_tracks.py` pins the parse and the endpoint values.

| file | source on `ubuntu@150.136.45.194` | sha256 |
|---|---|---|
| `diag_d.log` | `~/e0_run/d1_diag/diag_d.log` | `d2b49035e174db9b62dca95d3bed163cd6fc5e0aa8da41ec6af12231babb40f7` |
| `collapse_diag.log` | `~/e0_run/collapse_diag.log` | `bedadc2300d68e798d93c575b05eeb7e1a575a2362327411b104cdc4069601dc` |
| `diag_d.py` | `~/ws_d1/diag_d.py` | `7448d7c5b587ea8aa70f4a247d292ed4c3e4fff37de26f463c11104753a6977c` |

## The answer to the extraction, in two halves

**Half one — the rank track: per-step arrays EXIST and are recovered.** `diag_d.log` carries six
recorded steps (0, 25, 50, 100, 200, 400) of the clean in-batch InfoNCE arm, each with loss,
retrieval acc@1, positive cosine, worst-negative cosine, minimum margin, WSI-WSI off-diagonal cosine,
rank and feature std. This is the `12.88 → 1.00 by step 50` track that F8(b) stacks under the 16/16
pinning, and it no longer has to be drawn as two markers.

**Half two — the collapse-evidence track: per-step arrays DO NOT EXIST.** `collapse_diag.log` is the
sole surviving source of `0.7089 → 0.9999`, `0.0538 → 0.9959`, `0.0816 → 0.9960` and
`acc@1 0.062 → 0.000`, and it records **only** the endpoint pair for each quantity — the script probes
before and after, not on a schedule. F8(b)'s fallback instruction therefore stands for that half of
the panel: draw before/after paired markers labelled *"endpoint values as recorded; per-step array not
retained"*, and do not interpolate.

## The statistic, which was not what the figure spec says it is

`diag_d.py:50-51` computes its `eff-rank` column as:

```python
sv = torch.linalg.svdvals(wn.float() - wn.float().mean(0))
erank = float((sv.sum() ** 2) / (sv ** 2).sum())
```

where `wn = F.normalize(w, dim=-1)`. Rows are L2-normalised, the matrix is column-centred, and the
statistic is the **order-2** Hill number of the singular values. In the vocabulary of
`v2/calibra/spectral.py` that is exactly

```python
RANK_VARIANTS["R3"] == RankVariant(centre=True, normalise_rows=True, order=2)
```

and it is **not** `spectral.CANONICAL` (order 1, rows at own norms). `12.88 → 1.00` is an **R3**
number. The same inline formula appears in `diag_e/f/g/h.py` and `geom_probe.py` on the box, so every
`eff-rank` column under `~/e0_run/d1_diag/` is R3.

`~/ws_d1/momentum_test.py:102-104` — written later — imports `effective_rank` and reports R3 and the
canonical value side by side. The `d1_diag` family predates that fix.

The canonical value of this instance **cannot be recovered**: the diagnostic held no checkpoint and
wrote no state, so recomputing R1 requires re-running the training. It belongs in the same
`[NOT RECOMPUTED — needs a GPU forward pass]` category the draft already uses, and should not be
quoted as "the centred effective rank" without the R3 label. At the collapsed endpoint the label
matters little (a rank-1 matrix scores 1.00 under every variant in `RANK_VARIANTS`); at the 12.88
starting value it is unconstrained, because row normalisation changes the spectrum and the order-1
Hill number is `>=` the order-2 one on any fixed spectrum.
