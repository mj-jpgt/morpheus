# P2 figures

The nine display items of `paper/P2_FIGURES.md`, drawn with matplotlib and nothing
else. Each is a standalone script; each writes PDF, SVG and a 300 dpi PNG into
`out/`.

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
python make_all.py            # everything, in the draft's order
python make_all.py --only F2  # one item
python make_all.py --strict   # a pending measurement becomes a hard failure
python fig_f2_variance.py     # or run any script on its own
```

| item | script | carries |
|---|---|---|
| F1 | `fig_f1_envelope.py` | the retraining floor **at n = 5** and the seven arm ratios inside it |
| F2 | `fig_f2_variance.py` | the variance decomposition — the paper's most important display item |
| F3 | `fig_f3_floor.py` | the per-arm reproducibility floor, **at n = 5** |
| F4 | `fig_f4_defeaters.py` | the four-way defeater check |
| F5 | `fig_f5_verdict.py` | verdict instability across statistic, block and view |
| F6 | `fig_f6_necessity.py` | the necessity test, which went against us |
| F7 | `fig_f7_dilution.py` | the dilution dose–response and its miscalibration factor |
| F8 | `fig_f8_collapse.py` | the collapse boundary and the 16/16 withdrawal |
| F9 | `fig_f9_decorrelation.py` | rank rises while a co-measured collapse measure rises with it |
| T1 | `tab_t1_selection_rule.py` | the selection rule — with the underpowering inside the table body |

## Where the numbers come from

**No script carries a figure value as a literal.** Every number is read either
from a vendored box artifact under `data/` or from a markdown source in this
repository, parsed at draw time.

`data/` is written only by `extract_from_box.py`, which pulls the vendored files
in a single `tar` stream (never per-file `scp`), distils the two summaries that
are too large to vendor whole, and records every resulting file in
`data/MANIFEST.json` with its box path, size and SHA-256. `test_p2_figures.py`
re-hashes the tree against that manifest. The figures never talk to the box.

Refresh it with:

```
python extract_from_box.py
```

Two files are **deliberately not vendored**. `D1_PAIRED_BOOTSTRAP.json`
(unsuffixed) scores all 90 non-control targets, 50 of which are `programme_only`'s
own supervision, and F6 must not be able to reach it. `P2_METRICS_ALL_SUBSAMPLED.json`
is superseded by the `subsample` block inside the two vendored metrics files.

## What the scripts refuse to do

Each script cross-checks its sources against one another and **stops rather than
plotting either** when they disagree. In practice:

* F2 recomputes the variance decomposition from the twelve per-artifact metrics
  and asserts it equals the printed decomposition in the verified workspace's run
  log.
* F3 cross-checks the five-seed tripwire extraction against the three seeds that
  the canonicalisation run wrote independently.
* F4 recovers `gap/sd = gap / hypot(sd_A, sd_B)` and asserts all six pairs
  against the printed table.
* F5 asserts every selection mark against the printed table, and asserts
  `R2 == R3` on the raw exported artifacts before saying so on the panel.
* F6 reads its pre-declared thresholds from the script that declares them and
  cross-checks them against the notebook entry that records them as fixed before
  the pair list was inspected.
* F7 asserts the recomputed raw-R1 dilution curve equals the rank column printed
  in `DILUTION_LOWER_BOUND.md`.
* T1 recomputes every mark the metrics files can support, and recomputes the
  exact two-sided binomial *p* from `math.comb` rather than copying it.

## Pending measurements

`p2fig.pending_panel` draws a hatched placeholder that names the blocking
measurement, the path its data will arrive at, and its predeclaration — and
records it so `make_all.py` reports it at the end. Under `--strict` it raises
`PendingMeasurement` instead. It never silently plots nothing.

**None is open.** F1(d)'s five same-seed retraining repeats under
`~/e0_run/d1_envelope/` reported on 2026-08-04, and F1 refuses to draw the
placeholder now that they have, so the figure cannot go stale unnoticed.
`make_all.py --strict` returns 0.

## Conventions

* **Axis labels name the statistic and the block.** `p2fig.axis_label` will not
  build a rank axis without both. Four statistics have been called
  `effective_rank` in this repository and the raw-versus-residualised
  distinction flips arm orderings, so an unqualified "effective rank" axis is not
  a well-defined object.
* **No panel mixes two rank statistics on one axis.**
* **Palette: Okabe & Ito**, validated with the dataviz palette checker (all
  checks pass). Every series also carries a distinct marker and linestyle or
  hatch, so identity never rests on hue and the figures survive greyscale
  printing.
* **Binding statements are rendered into the artwork**, not left to a manuscript
  caption — F8's withdrawal, T1's underpowering, F6's negative, F1's floor-twice-over
  caveat, F9's one-seed-per-level.
* **A floor drawn beside a ratio is the floor of THAT ratio's own block,
  statistic and reading step.** F9 draws the fixed held-out probe's step-400
  floors (1.4489x R3, 1.5702x R1) and must not draw F1's 3.295x, which is
  canonical R1 on the residualised exported block of a different arm at 40
  epochs. Drawing the wrong floor is the error draft 4.1a exists to catch, and
  F9 committed it until 2026-08-05.
