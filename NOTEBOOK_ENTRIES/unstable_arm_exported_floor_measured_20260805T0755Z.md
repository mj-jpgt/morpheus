## 2026-08-05 07:55 UTC — The unstable arm's exported floor is **1.4254×**, not the ~5× I predicted, so **3.295× stands**. What does not stand is the paper's own badge of self-criticism: under the both-arm floor **RankMe as published is the LESS reproducible statistic (3.5484× against 3.1110×)**, the reversal of an abstract-level claim, and the paper's only selection clearing a floor on the exported block now clears nothing. The "12 of 12 resolvable" on the usable views falls to **7 of 12**

**Logged:** 2026-08-05 07:55 UTC. **Outcome: RESULT.**
**Predeclared in** `NOTEBOOK_ENTRIES/PREDECLARED_unstable_arm_exported_floor_20260805T0045Z.md`
(commit `a392c0a`), which fixed the prediction, the four readings and the flip thresholds **before the
runs were launched**. Interim record and the seed-axis supplementary:
`unstable_arm_exported_floor_seed_axis_and_the_queued_same_seed_run_20260805T0100Z.md` (`3eefc5f`).

---

## 0. The awkward findings, in order of how much they cost us

**(i) My predeclared prediction was wrong, and wrong in the paper's favour.** I predicted the
`programme_free` exported `wsi_biology` floor at **≈ 5×**, in an interval of 2.5×–12×, on the strength
of the probe block's "the collapsed arm carries the floor by ~2×". It measures **1.4254×** — less than
half my lower bound, and **2.3× smaller** than the stable arm's 3.2947×. That is predeclared outcome
**(D)**: *"the '~2× flattering' inference from the probe was wrong on this block — reported as a failed
prediction of mine, with the probe's own numbers left standing."* **§4.1's 3.295× is not an
underestimate. It is the conservative reading, and it stays.** The completeness audit's highest-priority
open item is answered, and it is answered in the paper's favour.

**(ii) But the same measurement reverses the paper's flagship self-criticism, and that is the finding
that matters.** Abstract item 8, §1.4 contribution 4 and §4.1b all carry: *"RankMe as published has a
retraining floor of **1.811×** on the raw exported block against our own centred statistic's
**3.111×** on the same five runs ... **by our own criterion the metric we criticise is more
reproducible than the one we chose to measure with**."* On `programme_free`'s five runs the same two
numbers are **RankMe 3.5484× against canonical R1 1.4395×** — RankMe is **2.46× worse**. Under the
both-arm floor the project already uses (`p2_probe_floors.combine()`: max of the arms), the comparison
is **3.5484× (RankMe) against 3.1110× (R1)** and **the direction is reversed.**

  **The paper's most-quoted piece of self-criticism was itself measured on one arm — the same defect,
  in the same place, that §4.1a exists to catch — and correcting it makes the paper look better.**
  Per this project's standing rule that a result in our favour gets the scepticism a result against us
  would get, it is reported here first and at full size.

**(iii) The consequence is that the exported block now has zero clearing selections.** §4.1a row 30
(`4.6-rankme-d2`) is the paper's one selection that clears a floor its own statistic and block license:
RankMe on D2 seed 43, **3.3817×** against 1.811×. Against the both-arm RankMe floor of 3.5484× it
**fails**. It is the **only** row in the whole 62-row audit whose verdict flips, and the selection
counts go **12 clear / 13 fail → 11 clear / 14 fail**. Every remaining clearing selection is §5's, on
the probe block.

**(iv) And the "usable views" half of the central claim loses five of its twelve comparisons.**
`rna_biology` **6/6 → 4/6**, `full_biology` **6/6 → 3/6**; **"12 of 12" becomes "7 of 12"**. §1.3's
*"the floor is fifty times smaller and every between-arm difference clears it"* is dead as written.
**This is predeclared outcome (B): the central claim survives in substance, its arithmetic does not.**

**(v) The one thing that softens (iv), and it is real rather than a rescue.** All five lost
comparisons are **D2 pairs (arm H against arm I)** — and *neither D2 arm is `programme_only` or
`programme_free`*, so this floor never licensed them in the first place. **All six D1 comparisons —
the only ones a `programme_only`/`programme_free` floor actually licenses — still clear, 3/3 on
`rna_biology` and 3/3 on `full_biology`.** The honest statement is that the count falls to 7 of 12 and
that the losses are concentrated exactly where the floor was already being stretched across
experiments (§6 below).

---

## 1. What was run

Five identical `programme_free` retrains, **seed 42 in all five**, 40 epochs, exported.
`~/chain_unstable_envelope.sh` is `~/chain_retrain_envelope.sh` — the script that produced 3.295× —
with **exactly one flag changed**, `--objective-profile programme_only` → `programme_free`. Same data
config, split file, architecture, optimiser, learning rate, weight decay, decorrelation, warmup,
`--programme-head-dim 256`, `--biology-key-momentum 0.999`, tripwire, `--gate-repeats 0
--rank-probe-repeats 0`, `--fit-development --fixed-final-epoch --restrict-to-split`, and the identical
`morpheus.v2.export` invocation. Five runs concurrent on one A100, as the stable arm's five were.

`~/e0_run/d1_envelope_pf/` — launched 05:37:22 UTC after the card cleared, **5/5 checkpoints** at
07:36:38, **5/5 exports** at 07:46:12. **Every repeat completed 40 epochs**, so the predeclared
"the unstable arm does not admit a floor because it does not reach 40 epochs" outcome did not fire and
no run was selected on.

**Queue discipline.** The A100 was saturated from 00:00 to 05:37 by another agent's twenty-run
600-step probe sweep (`~/ws_j2`, `launch_pf10.sh`). Nothing was launched into it; the chain polled
`nvidia-smi --query-compute-apps` — not a process-name match, so it could not match itself — and
required three consecutive clear checks 120 s apart.

**Workspace.** `~/ws_uf/morpheus`, `git -c core.autocrlf=false -c core.eol=lf archive` from `a392c0a`,
verified file by file by git blob SHA-1 against `git ls-tree -r`: **750/750, 0 missing, 0 extra, 0
differing**; tree object `b92991ad1e110c98b2dcee4e7195226317c07c96`, identical to `a392c0a^{tree}`.

**The two arms are code-matched in the training path, verified rather than assumed.** The stable arm's
five ran at commit `9cf6c84`. Against this workspace, after CR normalisation, `v2/runner.py`,
`v2/model.py`, `v2/data.py`, `v2/training.py`, `v2/losses.py`, `v2/contracts.py`, `v2/preflight.py`,
`v2/provenance.py`, `v2/plip.py`, `v2/pbs.py`, `v2/slide_pretraining.py`, `v2/export.py` and all of
`src/training/` are **identical**. The one module the runner imports that did change is
`v2/calibra/spectral.py`, and the two symbols it imports — `effective_rank` and `RANK_VARIANTS` — are
**identical by AST**. So the arm swap is the only difference in the training path.

**Readout.** `d1_envelope_readout.py` and `p2_envelope_floors.py`, both unchanged, both run from the
verified workspace. **No statistic was written for this entry**: R1/R2/R3 and the top-CCA channel are
imported from `v2/calibra/spectral.py`, RankMe / PR / PR_rownorm / stable rank / α-ReQ / LiDAR from
`p2_competing_metrics.py`, the hard rank is `numpy.linalg.matrix_rank`.

**The stable arm was re-scored through the same invocation in the same session** so the two arms are
comparable by construction rather than by trust. It reproduces the published values exactly —
**3.2947× / 3.1110× / 1.0193× / 1.0200× / channel 1.0553×**, and every per-repeat value in draft §4.1's
table (28.3202, 8.8340, 28.3482, 29.1057, 28.9588). **3.295× now has three independent sources.**

## 2. The floor, per repeat, never a mean

`programme_free`, exported artifact, held-out `test` partition, n = 2,766, canonical R1:

| repeat | R1 raw | **R1 residualised** | rna_biology R1 res | full_biology R1 res | **channel (top-CCA 16)** |
|---|---:|---:|---:|---:|---:|
| 1 | 9.0977 | **9.9004** | 14.3149 | 9.9526 | 0.5543 |
| 2 | 13.0964 | **14.0681** | 17.6142 | 14.0872 | 0.5476 |
| 3 | 11.1824 | **12.1054** | 14.6007 | 12.1962 | 0.5571 |
| 4 | 11.4008 | **12.2879** | 15.8621 | 12.3361 | 0.5611 |
| 5 | 12.9929 | **14.1120** | 16.8866 | 14.1530 | 0.5637 |
| **floor (max/min)** | **1.4395×** | **1.4254×** | **1.2305×** | **1.4220×** | **1.0293×** |

**Not bimodal on any view.** Under the rule fixed in `p2_envelope_floors._shape`, the WSI view gives
outlier rep1 (low), rest-fold 1.166×, separation 1.22 — nowhere near the 1.05× a four-run agreement
requires. **The stable arm's four-agree-plus-one-at-a-third signature does not reproduce on the
unstable arm**, so §4.1's *"reproducible about 80% of the time and catastrophically not about 20% of
the time"* is a property of `programme_only` on this block, not of the block and not of the method.
That is the third scope condition the sentence needs and does not have.

**The rank-versus-channel asymmetry survives on both arms**, which is the thing §4.1 actually rests on:
rank 1.4254× against channel 1.0293× here, rank 3.2947× against channel 1.0553× there. It is a factor
of 14.5 in excess-over-1 on the unstable arm and 41.5 on the stable one.

## 3. Both arms, every statistic — and the combined floor

Combined per the convention already fixed in `p2_probe_floors.combine()`: **`max` of the two arms'
folds, with the carrying arm named.** Never pooled and never averaged; a pooled fold across ten runs of
two arms would score the arms' genuine difference as noise.

**`wsi_biology`, residualised** — the block §4.1 measures

| statistic | `programme_free` | `programme_only` | **combined** | carried by |
|---|---:|---:|---:|---|
| **R1** | 1.4254 | **3.2947** | **3.2947×** | PO |
| RankMe | 1.4254 | **3.2946** | 3.2946× | PO |
| R3 | 1.4466 | **2.2901** | 2.2901× | PO |
| R2 | 1.4629 | **2.2244** | 2.2244× | PO |
| α-ReQ \|α−1\| | 1.7303 | **2.2985** | 2.2985× | PO |
| α-ReQ α | 1.5343 | **1.6923** | 1.6923× | PO |
| PR_rownorm | 1.2768 | **1.4663** | 1.4663× | PO |
| PR | 1.2464 | **1.4194** | 1.4194× | PO |
| stable rank | 1.1800 | **1.2239** | 1.2239× | PO |
| hard rank | 1.0000 | 1.0000 | 1.0000× | — |

**Every statistic on this view is carried by the stable arm.** The statistic axis §4.1b quotes —
1.000× to 3.295× on one block — is **unchanged** under the both-arm reading.

**`wsi_biology`, raw** — and this is where it reverses

| statistic | `programme_free` | `programme_only` | **combined** | carried by |
|---|---:|---:|---:|---|
| **RankMe (published)** | **3.5484** | 1.8111 | **3.5484×** | **PF** |
| **R1 (ours, centred)** | 1.4395 | **3.1110** | **3.1110×** | PO |
| R2 = R3 | 1.4306 | **2.1402** | 2.1402× | PO |
| α-ReQ \|α−1\| | 1.7479 | **2.1401** | 2.1401× | PO |
| PR = PR_rownorm | 1.2018 | **1.4464** | 1.4464× | PO |
| stable rank | 1.1674 | **1.2455** | 1.2455× | PO |
| hard rank | 1.0000 | 1.0000 | 1.0000× | — |

**RankMe as published is the only statistic on any view or block whose floor is carried by the unstable
arm**, and it is carried by a factor of two. Its per-repeat values there are **3.1164, 9.6465, 9.7732,
4.7458, 11.0583** — no degenerate near-zero, graded rather than bimodal (rest 2.330×). The statistic
axis on the **raw** block therefore no longer tops out at our R1: it runs 1.000× → **3.5484× (RankMe)**.

**The mechanism the paper gives is right and its sign is arm-conditional.** §4.6's account is that
RankMe's uncentred normalisation retains the mean-offset direction, which is large and stable, and that
every exported row has L2 norm exactly 1. On `programme_only` that direction *is* stable — RankMe raw
spans 1.99–3.60 while R1 raw spans 8.03–24.99. On `programme_free` the same uncentred direction spans
**3.12–11.06** while R1 raw spans 9.10–13.10. **The mean-offset direction is not reliably the stable
part; whether it is depends on the arm.** So the mechanism stands and the conclusion drawn from it does
not generalise beyond the arm it was measured on.

**`rna_biology` and `full_biology`, residualised** — here the unstable arm carries everything

| statistic | view | `programme_free` | `programme_only` | **combined** | ratio |
|---|---|---:|---:|---:|---:|
| **R1** | `rna_biology` | **1.2305** | 1.0193 | **1.2305×** | 1.21 |
| **R1** | `full_biology` | **1.4220** | 1.0200 | **1.4220×** | 1.39 |
| R2 | `rna_biology` | **1.3613** | 1.0231 | 1.3613× | 1.33 |
| R3 | `full_biology` | **1.4430** | 1.0186 | 1.4430× | 1.42 |
| stable rank | `rna_biology` | **1.3806** | 1.0269 | 1.3806× | 1.34 |
| α-ReQ \|α−1\| | `rna_biology` | **2.1592** | 1.0226 | 2.1592× | 2.11 |
| **LiDAR** (`wsi`/`rna` pair) | — | **2.5588** | 1.0596 | **2.5588×** | 2.42 |

**The arm effect on the floor inverts with the view, and the inversion is the result.** On
`wsi_biology` the stable arm is 2.31× noisier; on `rna_biology` and `full_biology` the unstable arm is
1.21× and 1.39× noisier, and on LiDAR 2.42×. `programme_only`'s RNA view is pinned at 27.22–27.75
across five retrains **because programme supervision holds it there**; `programme_free`'s three views
all sit at R1 ≈ 10–18 and move together. **1.019× is not a property of the RNA view. It is a property
of an arm whose objective supervises that view, and the paper has generalised it to "the view".**

## 4. What this does to the audit — one row flips, and it is the flagship

Recomputed read-only from `floor_audit.json`'s recorded ratios against the combined floors; nothing was
written to `floor_audit.json` (see §7).

| floor key | published | **combined** | × |
|---|---:|---:|---:|
| `RankMe_published_raw_export` | 1.811 | **3.5484** | 1.96 |
| `LiDAR_residualised_export` | 1.060 | **2.5588** | 2.41 |
| `LiDAR_raw_export` | 1.034 | **2.1708** | 2.10 |
| `R1_residualised_rna_view` | 1.019 | **1.2305** | 1.21 |
| `R1_raw_rna_view` | 1.023 | **1.2376** | 1.21 |
| `R1_residualised_full_view` | 1.020 | **1.4220** | 1.39 |
| `R1_raw_full_view` | 1.014 | **1.4355** | 1.42 |
| *all sixteen `wsi_biology` export keys* | *unchanged* | *unchanged* | 1.00 |

**Exactly one of the 62 rows changes verdict:**

| row | comparison | ratio | floor | was | **is** |
|---|---|---:|---|---|---|
| 30 / `4.6-rankme-d2` | RankMe as published, D2 seed 43 | **3.3817×** | 1.811× → **3.5484×** | clears | **fails** |

Selections: **12 clear / 13 fail → 11 clear / 14 fail**, of 25, 0 unjudgeable. **The exported artifact
block goes from one clearing selection to none.** §4.1b's *"RankMe is the only selection in this paper
that clears a floor on the exported artifact block"* becomes *no selection clears a floor on the
exported artifact block*, and every one of the eleven that do clear is §5's, on the probe.

## 5. The counts §4.5(c) and §4.1b rest on

Per pair, against each view's own combined floor:

| pair | `wsi_biology` (3.2947×) | `rna_biology` (1.2305×) | `full_biology` (1.4220×) |
|---|---|---|---|
| D2 s42 (H/I) | 1.5730 fail | 1.2054 **fail** *(was clear)* | 1.2345 **fail** *(was clear)* |
| D2 s43 (H/I) | 1.1858 fail | 1.1160 **fail** *(was clear)* | 1.0424 **fail** *(was clear)* |
| D2 s44 (H/I) | 1.0041 fail | 1.2380 clear *(by 0.61%)* | 1.1404 **fail** *(was clear)* |
| D1 s42 (P/F) | 2.1896 fail | 1.7659 clear | 2.2484 clear |
| D1 s43 (P/F) | 3.2463 fail | 2.8522 clear | 3.6057 clear |
| D1 s44 (P/F) | 1.7384 fail | 3.0137 clear | 5.2498 clear |
| **count** | **0/6 → 0/6** | **6/6 → 4/6** | **6/6 → 3/6** |

- **"12 of 12 resolvable on the other two views" → "7 of 12".**
- **Restricted to the D1 pairs — the only ones this floor licenses — 6 of 6 → 6 of 6, unchanged.**
- **All five losses are D2 pairs.** Restricted to D2: 6 of 6 → **1 of 6**, and the survivor clears by
  **0.61%**, which is not a margin anything should be quoted on.
- The view axis §1.4 calls **3.2×** becomes **2.68×** (3.2947 / 1.2305). §1.3's *"fifty times
  smaller"* — the excess-over-1 ratio, (3.295−1)/(1.019−1) = 121 — becomes **(3.2947−1)/(1.2305−1) =
  10.0**, i.e. **ten times smaller**, not fifty.
- The arm axis gains its first exported-block measurement: **2.31× on `wsi_biology` (stable noisier),
  1.21× on `rna_biology` and 1.39× on `full_biology` (unstable noisier)**. The paper currently quotes
  5.1× for this axis from the step-200 tripwire block only.

## 6. The scope point this measurement makes unavoidable

Three of the six pairs in every §4.5(c) / §4.1b count are **D2** (arm H, Hallmark, against arm I, PBS).
**Neither D2 arm is `programme_only` or `programme_free`, and no retraining floor of any kind has ever
been measured on a D2 arm.** Before today that stretch was invisible because the floor was 1.019× and
everything cleared it. It is now the entire difference between "12 of 12" and "7 of 12". Two readings
are available and the paper must pick one in writing:

1. **Conservative** — quote **7 of 12**, and say the five losses are on pairs whose arms have no
   measured floor, so the number is a floor-transfer and not a like-for-like verdict.
2. **Scope-matched** — quote **6 of 6 on the D1 pairs** and record the three D2 pairs as
   **unjudgeable**, which is what §4.1a's own block-matching rule does everywhere else.

Reading 2 is the one consistent with the rule the paper already enforces, and it costs the headline
count six comparisons rather than five. **Either way "12 of 12" cannot be printed.** Closing it
properly is five same-seed retrains per D2 arm — GPU, ten runs — and it is not currently on §6.2's list.

## 7. Prose locations the main session must change

`paper/P2_RANK_DRAFT.md` was **not edited by this session**, nor was `floor_audit.json`. `3.295` occurs
**91** times in the draft and `3.111` **27** times. They fall into five classes:

**Class A — UNCHANGED, no edit needed (the large majority).** Every use of 3.295× / 3.111× as the
`wsi_biology` exported floor, every "seven of seven fall inside a floor of 3.295×", every §4.1a row on
`wsi_biology`, the §4.1 repeat table, the 1.000×→3.295× statistic axis, and every historical sentence
("it once read as failing 3.295× by 0.3%"). The value is confirmed to four decimal places by a third
independent source. **What each of these should gain is the two words the measurement now licenses:
3.295× is the floor *of both arms*, carried by `programme_only`.** That is a strengthening, not a
correction.

*Line numbers are as of this commit against a 4,482-line draft and **will shift** — another session is
editing the same file concurrently. Each row therefore also gives the literal string to search for,
which is the durable locator.*

**Class B — the RankMe reversal. Must change; it is at abstract prominence.**
Search `1.811`: lines **42** (Status block item 8), **112**, **272**, **378**, **493**, **547**
(§1.4 contribution 4), **1434**, **1492** (§4.1a floors table, the RankMe raw row), **1558**, **1748**
(audit row 30), **1798**, **2002** (§4.1b box), **3989** (§6.2 closed-item row), **4187**, **4303**
(evidence-chain table), **4478** (scope caveat). In `paper/P2_FIGURES.md`: lines **22**, **175–176**,
**696**, **775**, **798** (S4 note (4)).
The correction is *not* a deletion: 1.811× against 3.111× remains **true on `programme_only`'s five
runs**. The new sentence is that **the direction is arm-conditional and reverses under the both-arm
floor — 3.5484× (RankMe) against 3.1110× (ours)** — which is this paper's own thesis one level up,
applied to its own instrument comparison.

**Class C — "the only selection that clears a floor on the exported artifact block".**
Search `only selection`: lines **47**, **276**, **494**, **550**, **4191**; and `exactly 1 clears` at
**112**. Becomes **no** selection clears on the exported block; all eleven that clear are §5's, on the
probe. `P2_FIGURES.md` line **696** carries the same sentence.

**Class D — the usable-view floors, counts and the view axis.**
Search `1.019`: lines **34**, **117**, **264**, **306**, **375**, **538**, **949**, **1434**, **1507**
(§4.1a floors table, `rna_biology` row), **1566**, **1604** (the view-axis table), **1936** (§4.1b box),
**1950** (the three-axis table), **2046**, **2162**, **2368**, **3989**, **4169**.
Search `12 of 12`: **37**, **267**, **542**, **4170**, **4474**; `all twelve` at **1936**;
`twelve of twelve` at **2370**.
Search `fifty times`: **262** and **480** (§1.3's claim box — the phrase is split across a line break
there, so a `fifty times smaller` search misses it).
Search `between views` / `between the co-trained views`: **306**, **538**, **1950**.
Audit rows **1745** (row 27, `rna_biology`) and **1746** (row 28, `full_biology`).
`P2_FIGURES.md`: **112**, **175**, **770** (the per-view table's `6 of 6`), **798**.
New values: floors **1.2305×** (`rna_biology`) and **1.4220×** (`full_biology`); counts **4 of 6** and
**3 of 6**; **7 of 12**; view axis **2.68×**; *"fifty times smaller"* → **"ten times smaller"**.

**Class E — `floor_audit.json`.** Seven floor values (§4) plus row `4.6-rankme-d2`'s `clears`
true → false and its `rests_on`, plus rows `4.5c-rna` / `4.5c-full` `rests_on` counts, plus the
`a`/`b` sources for those keys repointed at
`e0_run/d1_envelope_pf/out/P2_ENVELOPE_FLOORS_{PF,PO_RECHECK}.json`. **Deliberately not edited here**:
the §4.1a tables are generated from it and `v2/tests/test_p2_floor_audit.py` pins several of these
values as protective invariants, so changing the file without changing the draft in the same commit
would leave the two disagreeing. Four tests will need updating **as a decision, not as a test repair** —
in particular any that assert RankMe is the one clearing selection on the exported block.

## 8. Does the central claim survive? **Yes — with revised arithmetic, not a rewrite**

The claim is *"effective rank's usability as a selection signal is conditional on the co-trained view
it is read from and on the statistic it is read with."* Under the both-arm floor:

- **`wsi_biology` unusable: intact, and now arm-conservative.** 0 of 6 at 3.2947×, and the floor is
  confirmed as the *larger* of the two arms rather than assumed to be.
- **`rna_biology` / `full_biology` usable: intact in direction, halved in size.** 4 of 6 and 3 of 6
  against 0 of 6 is still a categorical difference between views on the same artifacts with the same
  statistic; on the D1 pairs the floor licenses it is 6 of 6 against 0 of 6, unchanged.
- **The view-conditionality itself: intact but re-attributed.** The gap between views narrows from
  3.2× to 2.68×, and the measurement says part of what looked like a view effect is an arm effect —
  `programme_only`'s RNA view is stable because its objective supervises it. **That is a genuine
  weakening of the mechanism, not of the phenomenon**, and §4.1b should say so.
- **The statistic-conditionality: intact and unchanged** on the residualised block; on the raw block
  the spread now runs to RankMe's 3.5484× rather than our R1's 3.1110×.
- **§4.2's variance decomposition, which the paper says it rests on: untouched.** It uses no floor.

**So: no rewrite of §1, no change of title, no withdrawal.** What must change is (a) the RankMe
sentence in the abstract, (b) "the only selection that clears", (c) every "12 of 12" / "fifty times
smaller" / "3.2× between views", and (d) the two view floors. The one place the paper gets *stronger*
is 3.295× itself, and the one place it gets weaker is the sentence it was proudest of.

## 9. What this still is not

n = 5 per arm, **one seed**, one stack, one architecture, one cohort, no interval. It is a floor once
over for the arm axis now and still a floor for the seed axis: §4.2's seed term is separate and larger.
It says nothing about the D2 arms, the probe block, the in-run training batch, the 16-patient gate batch
or the 282-patient live checkpoint. `no_external_cohort` is undischarged. **"The floor is 1.4254×" must
never be written as "rank varies 1.4254×".**

## 10. Suite

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python -m pytest v2/tests tests -q --basetemp=./pytmp
```

Baseline at `288a124`, taken before this work began:

```
672 passed, 2 skipped, 445 warnings in 154.67s (0:02:34)
```

Final, at `b33b2dc` plus this entry's one code change:

```
696 passed, 1 skipped, 447 warnings in 133.17s (0:02:13)
```

Both verbatim. **Zero failures in either.** The counts differ because three other agents committed to
this tree between the two runs (P5 pilot funnel, the labelled linear probe, the P2 limit-2 stress
sweep); **no test was added, changed, skipped or deleted by this entry.** The repository is reached as
`morpheus/` through a directory junction outside the tree, so nothing in the repo was modified to make
the suite import.

## 11. Files

- **Box:** `~/e0_run/d1_envelope_pf/rep{1..5}/`, `rep{1..5}.npz`,
  `out/{P2_ENVELOPE_FLOORS_PF.json, P2_ENVELOPE_FLOORS_PO_RECHECK.json, d1_envelope_pf_readout.log,
  floors_pf_run.log, floors_po_recheck.log}`; `~/e0_run/pf_seedaxis/out/SEEDAXIS_{f,p}.json`;
  `~/chain_unstable_envelope.sh`; `~/ws_uf/morpheus` (750/750, tree `b92991ad`)
- **Vendored** through the only sanctioned path, `figures/extract_from_box.py` (additive block; every
  previously vendored file re-hashed byte-identical) → `figures/data/e0_run/d1_envelope_pf/out/`,
  `figures/data/e0_run/pf_seedaxis/out/`, `data/MANIFEST.json`
- **Run unchanged:** `v2/research/rebase/p2/p2_envelope_floors.py`,
  `v2/research/rebase/d1_envelope_readout.py`
- **Not touched:** `paper/P2_RANK_DRAFT.md`, `paper/P2_FIGURES.md`,
  `v2/research/rebase/p2/floor_audit.json`, `v2/calibra/claim_guards.py`, `claim_evidence.json`, any
  other agent's `PREDECLARED_*`, and every file another agent had uncommitted in the shared tree.
