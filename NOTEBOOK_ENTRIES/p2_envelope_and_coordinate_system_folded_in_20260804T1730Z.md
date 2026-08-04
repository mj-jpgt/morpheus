## 2026-08-04 17:30 UTC — Both landed results folded into the papers; and the selection-rule verdict is NOT stable across ground-truth blocks — 12 of 12 metric rows move, and our own statistic can be handed a p = 0.031 by choosing the exam's coordinate system

**Logged:** 2026-08-04 17:30 UTC. **Writing and figure work only, CPU only; the GPU was not touched.**
Two results that had reported but had not reached the manuscripts were folded in, the figures they
change were regenerated, and one consequence that nobody had chased was computed. Thread caps
`OMP/OPENBLAS/MKL/NUMEXPR=1` throughout. **No statistic is computed inline anywhere in this work** —
every rank and channel value is read from the log written by
`v2/research/rebase/d1_envelope_readout.py`, which imports `calibra.spectral` and
`calibra.residualise`, or from `EXAM_PANEL.json`, written on the box by the same `d2_compare`
machinery.

### 0. Bad news first: the selection-rule verdict is not stable, and it is unstable in our favour as well as against us

`paper/P2_RANK_DRAFT.md` §4.6 scores twelve label-free metrics as selection rules against **one**
ground truth: the held-out channel onto 40 gene-set targets. The coordinate-system result established
that that arm contrast is a property of the target block. Nobody had chased what that does to §4.6.
Holding all twelve metrics fixed and swapping only the truth, once per block
(`v2/research/rebase/p2/p2_selection_rule_blocks.py`):

| exam block used as ground truth | Δ s42 | Δ s43 | Δ s44 | arm ordering |
|---|---:|---:|---:|:---:|
| gene sets, untrained 40 — **the published truth** | −0.1325 | −0.1089 | −0.1226 | H H H |
| PBS codes 128 (arm I's own supervision) | −0.0098 | **+0.0088** | −0.0026 | H **I** H |
| PCA basis 128 | −0.0201 | **+0.0049** | −0.0284 | H **I** H |
| gene-label-shuffled PBS 128 | −0.0359 | −0.0057 | −0.0175 | H H H |
| `random_control` gene sets 90 | −0.0099 | −0.0280 | −0.0268 | H H H |
| random dictionary 128 | −0.0597 | −0.0132 | −0.0454 | H H H |

**Two distinct winner patterns, and everything turns on seed 43** — the one seed whose arm ordering
the coordinate system flips.

| metric | gene sets 40 | PBS codes | PCA basis | shuffled PBS | rand ctrl | rand dict | stable? |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--|
| canonical effective rank (raw) | 2/3 | 3/3 | 3/3 | 2/3 | 2/3 | 2/3 | **no** |
| canonical effective rank (resid.) | 2/3 | 3/3 | 3/3 | 2/3 | 2/3 | 2/3 | **no** |
| RankMe (raw, as published) | 3/3 | 2/3 | 2/3 | 3/3 | 3/3 | 3/3 | **no** |
| RankMe (residualised) | 2/3 | 3/3 | 3/3 | 2/3 | 2/3 | 2/3 | **no** |
| participation ratio (raw / resid.) | 2/3 | 1/3 | 1/3 | 2/3 | 2/3 | 2/3 | **no** |
| stable rank (raw) | 2/3 | 1/3 | 1/3 | 2/3 | 2/3 | 2/3 | **no** |
| stable rank (resid.) | 1/3 | 0/3 | 0/3 | 1/3 | 1/3 | 1/3 | **no** |
| α-ReQ \|α−1\| (raw / resid.) | 2/3 | 3/3 | 3/3 | 2/3 | 2/3 | 2/3 | **no** |
| LiDAR (raw) | 0/3 | 1/3 | 1/3 | 0/3 | 0/3 | 0/3 | **no** |
| LiDAR (residualised) | 1/3 | 0/3 | 0/3 | 1/3 | 1/3 | 1/3 | **no** |

**Answer to the question as asked: the verdict is not stable. All twelve metric rows change their D2
count when the ground-truth block changes; not one is stable.** Three consequences, and the second is
the one that matters most:

1. **The ordering between the two rows §4.6 quotes against one another reverses on 2 of 6 blocks.** On
   the published gene sets RankMe-as-published beats canonical effective rank on D2, 3/3 against 2/3 —
   a fact §4.6 reports precisely because it cuts against us. On the dictionary's own codes and on a
   plain PCA basis the two **swap**.
2. **Our own statistic can be handed a nominally significant result by choosing the block.** Canonical
   effective rank goes 5/6 → **6/6 overall, exact two-sided p = 0.031** on two of the six blocks —
   "significant" by §4.6's own stated bar, produced by nothing but the coordinate system the exam is
   written in. **This is reported because it disqualifies the favourable reading as firmly as the
   unfavourable one**, and it is the strongest available argument that no count in §4.6 may be quoted
   without its target block.
3. **Only the D2 half can be re-scored at all.** The coordinate-system work re-scored the two **D2**
   arms; the D1 arms were never scored against any block but the gene sets. The D1 column is therefore
   held fixed in every row above and the ALL counts inherit that. **That is an absent measurement, not
   evidence that the D1 half is block-stable**, and it is now a named row in §6.2's missing-measurement
   table.

What *is* stable across all six blocks is the direction: arm H is never behind anywhere. That is the
same thing `D2_RESULT.md` §6 now says, and it is all that survives.

*Validation of the counterfactual itself: the `geneset_untrained40` column reproduces §4.6's published
table exactly, in all twelve rows, which is what licenses reading the other five. Three tests in
`v2/tests/test_p2_analysis_scripts.py` cover it, including one that feeds a panel agreeing with the
published truth everywhere and requires every row to come back `stable? yes`.*

---

### 1. Result 1 — the retraining envelope, folded into §4.1, the abstract and §4.7

Five identical `programme_only` retrains at seed 42 (`~/e0_run/d1_envelope/rep{1..5}`), readout
`~/e0_run/d1_envelope_readout.log`. **Verified before use**: the box copy of
`~/ws_d1/morpheus/v2/research/rebase/d1_envelope_readout.py` is byte-identical to repository HEAD by
git blob SHA-1 (`4ba868ece8a7a1d84c39f758503413182a6b5216`), so the log's provenance is not a matter
of trust.

| repeat | rank raw | rank resid. | channel |
|---|---:|---:|---:|
| 1 | 24.481 | 28.320 | 0.6182 |
| **2** | **8.033** | **8.834** | **0.5859** |
| 3 | 24.504 | 28.348 | 0.6123 |
| 4 | 24.990 | 29.106 | 0.6110 |
| 5 | 24.912 | 28.959 | 0.6098 |
| spread | **3.111×** | **3.295×** | **1.055×** |

**What moved in §4.1.** Rewritten around the measured floor. The n = 1 2.69× estimate is retained only
where it is labelled superseded. Added: the bimodality as the primary framing ("reproducible ~80% of
the time and catastrophically not ~20% of the time"), repeat 2 as the paper's cleanest single
observation, the floor-twice-over argument, and the count moving **six of seven → seven of seven**
reported with the predeclared scepticism, including that the seventh point clears by **1.5%**.

**A block subtlety that had to be handled and is now handled explicitly.** The floor is 3.295×
residualised and **3.111× raw**. D1-B seed 43's ratio is 3.246× residualised and 3.09× raw. Both are
inside their own block's floor — but a residualised ratio judged against the *raw* floor would fall
outside. §4.1's table, F1(b) and the figure plan now all judge each ratio against the floor measured
on its own block, and say so. The predeclaration and the earlier notebook entry quote D1's raw-block
ratios (2.02 / 3.09 / 1.68) while §4.1 and §4.7 quote the residualised ones (2.190 / 3.246 / 1.738);
both are now stated, with their blocks named, so the two cannot be confused for a disagreement.

**What moved in the abstract.** Repeat 2 is in it. The long abstract's item (i) is rebuilt around the
measured floor, the bimodality and the 3.295× / 1.055× asymmetry on identical inputs; the short
abstract likewise. The necessity paragraph now carries that its rank differences sit inside the floor,
so the test is **unresolvable rather than refuted** — with the predeclaration named in the same
sentence, so a reader can check that the reading was fixed before the number existed.

**What moved in §4.7.** The section keeps its title, its position ahead of everything favourable, its
table and its intervals. A blockquote at its head applies the predeclared `> 3.09×` band: D1 is
uninformative about rank **in either direction**; the necessity result is **not refuted**; the
asymmetry is the finding and holds only because both halves are reported. §4.7.4 gains an explicit
statement of the strongest objection — *"you made your own inconvenient result disappear by measuring
a wide enough noise floor"* — and three answers, none of which is that the result went away.

**§4.3 corrected while in the area.** Its step-200 tripwire table showed three seeds (`programme_only`
spread 1.003×) while F3 has been drawn at five since seeds 45 and 46 landed. The table is now the
five-seed one: `programme_only` **1.18×**, `programme_free` **6.05×**. Text and figure agreed on the
6.05× and disagreed on the sibling by a factor of 60; they now agree. Every downstream quotation of
"1.003×" in the draft was updated or scoped to the three seeds it describes.

### 2. Result 2 — D2's coordinate-system dependence, folded into `D2_RESULT.md` §6

§6's verdict sentence narrowed from *"on the held-out molecular channel"* to **"on gene-set–valued
targets, on the confound-residualised block"**, with both qualifications spelled out as load-bearing:
the six-block panel (−0.12 on exactly one block, every other at or inside the negative control), and
the raw block (gap 3–5× smaller, **seed 43 reverses sign**, point estimates only — the 3/3 sign
consistency and the ~0.12 magnitude are properties of residualisation). The surviving claim is stated
at its real strength — Hallmark is never worse anywhere tested and much better on gene sets, PBS buys
~0.01–0.04 in its own neighbourhood — and it still refutes P3's hypothesis, but **not as a 0.12
general effect**. The survival result and the empty T1.7(b) coverage travel with it, and §0's
consequence for P2 is recorded there too.

### 3. Figures regenerated

All nine display items re-rendered from the vendored data (`make_all.py`); four changed materially.

- **F1 — rebuilt.** Panel (d) was a hatched `[RETRAINING ENVELOPE PENDING]` placeholder and the script
  was written to **raise** rather than emit it once the repeat reported; it now does something else
  instead. New (a): the five repeats, rank and channel as stacked strips on a shared repeat axis,
  plotted **individually and never as a mean or band** because the distribution is bimodal. New (b):
  the seven ratios against the measured floor, each judged against **its own block's** floor, 7/7
  inside. (c) unchanged. New (d): the measured floor against the superseded n = 1 estimate, with D1's
  three ratios beneath and the "uninformative in either direction / not refuted / in our favour"
  statement in text. The figure recomputes the three spreads from the per-repeat values and
  **asserts them against the spreads the readout log itself printed** before drawing anything.
- **T1 — a second band added.** Below the table: the six blocks as columns, their arm orderings, and
  the D2 counts for the two rows §4.6 quotes, each with its ALL/6 and exact binomial. Closing line
  states the reversal and the manufactured 6/6. The band's counts are recomputed in the figure from
  the same per-artifact metrics the table's marks are recomputed from, and the exam panel's gene-set
  column is asserted against the metrics JSON's own untrained-40 contrast first.
- **F6(b) — the measured floor drawn on it.** §4.7 now reads this panel against the floor, so the
  panel carries it: a dashed line at 3.295×, and a footnote stating that 3/3 are below it, that D1
  therefore resolves rank in neither direction, and that this does not refute the necessity result
  opposite. It reads the floor from the same extracted file F1 does, so the two cannot diverge.
- **F2 — nothing numerical moved**, correctly: the variance decomposition is over the same twelve
  artifacts and neither result touches it. Its docstring was updated to say that its independence
  from F1's floor mattered when the floor was n = 1 and still matters now that it is measured.

F3, F4, F5, F7, F8 unchanged.

### 4. Anything that disagreed with its source

**Nothing that forced a stop.** Three things were reconciled rather than plotted:

1. **D1's rank ratios appear in two forms** — 2.02 / 3.09 / 1.68 (raw) in the predeclaration and the
   envelope entry, 2.190 / 3.246 / 1.738 (residualised) in §4.1 and §4.7. Not a disagreement: two
   blocks. Both are now quoted with their block named, and each is judged against the floor for its
   own block. Recorded because the numbers look like a contradiction and are not.
2. **§4.3's tripwire table was three seeds where F3 draws five** — `programme_only` 1.003× against
   1.18×. Corrected in the draft; the figure was right.
3. **"3.30×" against "3.295×".** The readout log prints `spread=3.295x`; the earlier notebook entry
   and its commit message round it to 3.30×. Same number. The papers and figures quote **3.295×** at
   three decimals throughout so the ambiguity cannot recur.

### 5. Suite

**377 passed, 0 failed** (`v2/tests` + `tests`, thread-capped, matplotlib present locally). Baseline
before this work was 373; the four new tests are three covering the block counterfactual and one
asserting that **no display item is pending** — which replaces the F1-specific `--strict` test that
had become a permanent skip once the last placeholder was retired. The `--strict` mechanism itself is
still covered, as a unit test on `pending_panel`. Two failures were caught and fixed en route, both
from `test_paper_paths_resolve.py`: five repository paths cited in a shorthand form
(`data/extracted/…`, `d2_coordinate_system/out/…`) that does not resolve from the repository root.

*`v2/tests/test_p2_figures.py` still cannot run in the box's `~/venv` — no matplotlib — and nothing
was installed there.*

### In plain terms

We had two finished experiments sitting outside the papers. The first says that training the same
model five times with everything identical gives the same rank four times and a number three times
smaller once, while the information content barely moves; that is now the headline, and the cost is
that our own inconvenient experiment falls inside that noise and can no longer tell us anything about
rank in either direction — which we say rather than treat as the result going away.

The second says the big D2 effect only exists when the exam is written in one particular set of units.
Following that through to where nobody had — the table in the rank paper that uses the D2 effect as
its answer key — the answer key turns out to move too. Every row of that table changes when you change
the units, the two methods it compares swap places, and we can give our own preferred metric a
"significant" result just by picking the units. We report that, because being able to manufacture a
result in our own favour is the clearest possible reason for nobody to believe any number in that
table without being told which units it was scored in.

### Files

- `paper/P2_RANK_DRAFT.md` — §4.1 rewritten, §4.6a new, §4.6 and §4.7 reframed, §4.3 corrected to
  five seeds, abstracts, status banner, §6.2, conclusion, Appendix A and C
- `paper/P2_FIGURES.md` — F1 and T1 rows rewritten, binding constraint 2, S4 status, pending table
- `v2/research/rebase/nature/D2_RESULT.md` — §6 narrowed
- `v2/research/rebase/p2/p2_selection_rule_blocks.py` — new
- `v2/research/rebase/p2/figures/{fig_f1_envelope,fig_f2_variance,fig_f6_necessity,tab_t1_selection_rule,extract_from_box,make_all}.py`
- `v2/research/rebase/p2/figures/data/` — `e0_run/d1_envelope_readout.log` vendored,
  `extracted/F1_RETRAINING_REPEAT.json` and `MANIFEST.json` refreshed
- `v2/tests/{test_p2_analysis_scripts,test_p2_figures}.py`
- Sources: `NOTEBOOK_ENTRIES/retraining_envelope_and_momentum_seeds_20260804T1000Z.md`,
  `NOTEBOOK_ENTRIES/d2_coordinate_system_result_20260804T0800Z.md`,
  `NOTEBOOK_ENTRIES/PREDECLARED_retraining_envelope_20260804T0330Z.md`

### Left undone, named rather than left to be discovered

- **§5.2 and §5.3 are not rewritten around the momentum seed replication.** It landed in the same
  entry as the envelope and is outside this task's scope, but the status banner, §5.3, Appendix C and
  the figure plan's S4 row now all record that it reported, that it resolves §5.3's disjunction in
  favour of separation, **and that its worst-case separation is 3.29× against the 3.295× floor** — so
  by §4.1's own criterion the momentum fix's rank difference is not resolvable either. Stated rather
  than left for a referee.
- **The D1 arms have no per-block ground truth.** Named in §6.2.
