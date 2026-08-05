## 2026-08-05 01:50 UTC — The labelled linear probe, on all 12 artifacts and all 3 views. **It has rank's reproducibility failure, not the channel's** — and on §4.1's cleanest observation it moves with rank, 1.38× against the channel's 1.055×. The molecular half of the probe corroborates the channel, 24 of 27

**Predeclaration:** `NOTEBOOK_ENTRIES/PREDECLARED_p2_labelled_linear_probe_20260805T0040Z.md`, committed
at `6b3d8e7` **before any statistic in this entry existed**. Outcomes **A, B, C and D all occur**, on
different views — which the predeclaration did not anticipate and which is stated as such below.

**How obtained:** box `150.136.45.194`, workspace `~/ws_p2probe`, shipped as
`git -c core.autocrlf=false archive HEAD` and **verified byte-equal to HEAD, 753/753 tracked files by
per-file git blob SHA-1**, before anything ran. `~/venv`, threads capped to 1. CPU only, no GPU, no
retraining — every probe is fitted on frozen exported embeddings. Runner
`v2/research/rebase/p2/p2_labelled_probe.py`, vendored and tested (`v2/tests/test_p2_labelled_probe.py`)
**before** it was run. Every statistic imported from `v2/calibra/`; nothing computed inline. Output
`~/ws_p2probe/out/P2_LABELLED_PROBE.json` (merged from four shards + the four shard files).

---

## 0. The gap, and what it now says

§2.5 ended `[STILL NOT MEASURED] — a labelled linear probe on every artifact, which is the ground truth
LiDAR and RankMe were validated against`; §6.2 carried the same row; §6.3 concedes Zaiem et al. against
us on exactly that; `paper/P2_FIGURES.md` records it twice; and the completeness pass of
2026-08-04 called it *"the single most valuable missing measurement in the paper"*. It is now measured.

**Two probes, both predeclared.** **Probe A** — cancer type (21 held-out lineages), out-of-fold linear
classifier on the **raw** block, balanced accuracy: the literal analogue of RankMe's and LiDAR's own
validation. **Probe B** — six molecular / clinical endpoints from E1's MC3-derived label table,
confound-adjusted, size-weighted within-cancer AUROC: matched to the construct *this paper* argues
about. Neither is nominated the winner over the other.

---

## 1. The awkward finding, first: the labelled probe has rank's reproducibility failure, view for view — and the channel does not

Five identical `programme_only` retrains at seed 42 — the same runs §4.1's **3.295×** rank floor and
**1.055×** channel spread are measured on. Every statistic re-read on them:

| statistic, five same-seed retrains | `wsi_biology` | `rna_biology` | `full_biology` |
|---|---:|---:|---:|
| canonical R1, residualised (§4.1a) | **3.295×** | 1.019× | 1.020× |
| held-out channel, untrained-40 (§4.1) | 1.055× | 1.007× | 1.005× |
| **Probe A — cancer type, logistic** | **1.380×** (Δ 0.2076 bal. acc.) | 1.011× | 1.012× |
| **Probe A — cancer type, LDA** | **1.511×** (Δ 0.2270) | 1.010× | 1.012× |
| Probe B — TP53 within-cancer AUROC | 1.007× (Δ 0.0040) | 1.015× | 1.011× |

**§4.1b's central claim was that rank's usability is conditional on the co-trained view. That
conditionality is not a property of rank. The labelled reference standard has it too, in the same
direction, on the same runs.** The view where rank is irreproducible is the view where a labelled
linear probe is irreproducible; the two views where rank is reproducible to 1.02× are the views where
the probe is reproducible to 1.01×. What §4.1b reads as *"the difference between 'unusable' and
'usable' is a choice of co-trained view … that no paper we have read reports making"* is, on this
evidence, a real property of what the WSI encoder retrains to — visible to the metric under criticism
and to the reference standard alike.

### 1a. Repeat 2 — §4.1 calls it "the cleanest single observation this project has produced", and the probe changes what it observes

| repeat | canonical R1 | **Probe A, logistic** | Probe A, LDA | channel | Probe B, TP53 |
|---|---:|---:|---:|---:|---:|
| 1 | 28.320 | 0.7413 | 0.6571 | 0.6182 | 0.5811 |
| **2** | **8.834** | **0.5459** | **0.4440** | **0.5859** | **0.5839** |
| 3 | 28.348 | 0.7498 | 0.6709 | 0.6123 | 0.5851 |
| 4 | 29.106 | 0.7535 | 0.6675 | 0.6110 | 0.5821 |
| 5 | 28.959 | 0.7442 | 0.6620 | 0.6098 | 0.5824 |

§4.1 currently reads this row as: *"effective rank falls 3.3×; the channel falls 5%. One quantity moves
by a factor, the other barely moves at all, and there is no difference between the two runs to
attribute either move to."*

**The measurement that sentence was missing is now in it. Two quantities move by a factor and one does
not.** On repeat 2 the labelled cancer-type probe falls **27%** (0.7413–0.7535 → 0.5459), i.e. **1.38×**,
against the channel's 5%. The rank collapse on that run corresponds to a real, large, measurable loss
of **linearly decodable content**, seen by the exact reference standard RankMe and LiDAR were validated
against — and *not* seen by this paper's ground truth.

**The reconciliation, and it narrows the paper's claim rather than destroying it.** What the probe
loses on repeat 2 is *cancer type* — the lineage direction the channel residualises **out** by
construction. The two measurements are therefore not in contradiction: on `wsi_biology`, canonical R1
appears to track the lineage direction, which is precisely the direction §3.2's readout is built to
remove. **"Rank moved and nothing else did" is not a sentence this paper can keep. "Rank moved with the
confound the channel removes, and the channel did not move" is a sentence it can keep, and it is a
better one** — it says what rank is tracking instead of only what it is not.

---

## 2. Where the probe can resolve a between-arm difference, it sides with RANK against the channel more often than with the channel

Applying §4.1's own criterion to the probe — a between-arm difference counts only if it exceeds the
**same statistic's** five-repeat floor on the **same view**:

| view | Probe A resolvable | Probe A agrees with **channel** | Probe A agrees with **rank** |
|---|:--:|:--:|:--:|
| `wsi_biology`, logistic | **0 / 6** | 3/6 | 4/6 |
| `wsi_biology`, LDA | 1 / 6 | 5/6 | **6/6** |
| `rna_biology`, logistic | 6 / 6 | **0/6** | 3/6 |
| `rna_biology`, LDA | 6 / 6 | 3/6 | **6/6** |
| `full_biology`, logistic | 6 / 6 | 3/6 | **6/6** |
| `full_biology`, LDA | 6 / 6 | 3/6 | **6/6** |

**On `wsi_biology` the labelled probe is exactly as unusable as rank: 0 of 6 differences clear its own
retraining floor.** The reference standard the proxy is supposed to substitute for cannot adjudicate
the readout view either. §4.1's practitioner rule — *measure the floor on the view and with the
statistic you are actually going to compare with* — applies verbatim to a labelled linear probe, which
is not a caveat on the rule but its strongest support.

**On the two views where everything resolves, the probe picks rank's arm, not the channel's.** All
three D2 pairs go to arm **I** under Probe A on `rna_biology` and `full_biology` — the same three pairs
§4.5(c) and §4.1b score as rank's *wrong* answers. Aggregated over every (probe × view × pair) conflict
between a probe and the channel, **rank sides with the probe 14 times and with the channel 9**.

**This is the predeclaration's outcome D and it must be reported first, as it was fixed in advance.**
§4.1b's sharpest sentence — *"the view that makes rank usable is the view on which it is most often
wrong"* — is true **only relative to the channel**. Relative to a labelled linear probe on those same
views, rank is right on 6 of 6 pairs. The sentence needs the words "against our unsupervised readout"
in it, and §4.5(c)'s 3/6 and 2/9 counts need the same qualification.

---

## 3. The favourable half, reported with its own distrust checks applied

**Probe B — the molecular / clinical endpoints — corroborates the channel.** On `wsi_biology`, over the
six-endpoint panel × six pairs (36 comparisons): **29 agree with the channel; 27 clear the endpoint's
own five-repeat floor; of those 27, 24 agree and 3 disagree.**

| endpoint, `wsi_biology` | own floor (Δ) | resolvable | agrees with channel |
|---|---:|:--:|:--:|
| `mut_TP53` | 0.0040 | 6/6 | **3/6** |
| `grade_high` | 0.0138 | 5/6 | 5/6 |
| `stage_late` | 0.0524 | **0/6** | (5/6, none resolvable) |
| `mut_ATM` | 0.0248 | 6/6 | **6/6** |
| `mut_KMT2D` | 0.0271 | 6/6 | **6/6** |
| `mut_ARID1A` | 0.0104 | 4/6 | 4/6 |

**All three resolvable disagreements are the same endpoint on the same experiment**: TP53 picks
`programme_free` over `programme_only` in all three D1 seeds (−0.0107 / −0.0212 / −0.0262 against a
0.0040 floor), where the channel, rank and every other statistic pick `programme_only`. Taken alone
that would say the necessity result of §4.7 is scored against a ground truth the labelled standard
contradicts 3/3.

**It does not survive the panel, and distrust check 4 was predeclared for exactly this.** On the same
three D1 pairs, `mut_ATM` and `mut_KMT2D` pick `programme_only` **3/3 each**, `grade_high` 2/3 and
`mut_ARID1A` 2/3. **One endpoint out of six is a coordinate choice, which is §4.6a's own finding
restated on the labelled side.** The TP53 result is reported, and it is not evidence that D1's ordering
is wrong.

---

## 4. Controls, nulls, and every distrust check the predeclaration named

| check | result |
|---|---|
| **Must-fail control** — cancer type must collapse under cancer + pooled-TSS adjustment | **PASSES, 17/17.** Adjusted LDA balanced accuracy **0.0340–0.0598** against chance 1/21 = 0.0476 and a *measured* 200-draw null p95 of 0.0578–0.0620; permutation p 0.0746–0.9950. The closest to breaching is F44 at 0.0598, p = 0.0746 — inside, and named rather than rounded past. |
| 2. Probe above its own measured null | **17/17 for Probe A** (LDA raw 0.361–0.723 against null p95 ≈ 0.059, permutation p = 0.0050 = the 1/201 resolution floor; logistic p = 0.0476 = the 1/21 floor at 20 draws, **resolution-limited and disclosed**). **17/17 for TP53** (within-cancer 0.5640–0.6017 against a measured within-cancer null p95 of 0.5286–0.5324, p ≤ 0.002). |
| 1. Difference inside the probe's own floor | **Fires hard on `wsi_biology`:** 0/6 for Probe A. Those six pairs are `UNRESOLVED` and are scored as neither agreement nor disagreement, as predeclared. |
| 3. Agreement carried by pooled rather than within-cancer AUROC | Reported side by side. TP53 pooled agrees with the channel 3/6, the same as within-cancer, and the pooled column is a lineage proxy (E1's BRAF/KRAS/APC lesson) — no conclusion here rests on it. |
| 4. Agreement on one endpoint only | **Fires — see §3.** The only resolvable disagreements are TP53's, and the rest of the panel contradicts them. |
| 5. Ordering flips between estimators | **Fires.** Logistic and LDA disagree on the winner in **2/6** pairs on `wsi_biology`, **3/6** on `rna_biology`, 0/6 on `full_biology`. **We have reproduced Zaiem et al.'s critique against ourselves**, on our own probe, and §6.3 should say so. |
| 6. `rna_biology` / `full_biology` circularity | Declared before the run. Probe B on those views is reported and is **not** used to adjudicate any morphology claim. The panel endpoints were scored on `wsi_biology` only; their rna/full columns are absent, not null. |
| 7. Adjustment not working | Covered by the must-fail control, which passes. |

**Reproduction, as a validity check on this run's code path.** Every published number this run
recomputes reproduces **to every published digit**: all twelve artifacts' canonical R1 and untrained-40
channel (§4.1, §4.6), the five retrains' rank and channel (§4.1), the three per-view floors
3.295× / 1.019× / 1.020× and the channel spread 1.055× (§4.1a), and E1's six `d2_h_seed42` endpoint
AUROCs (0.5912 / 0.6777 / 0.5469 / 0.6135 / 0.5997 / 0.5981). A redundant second pass over the five
retrains, run by accident, reproduced the panel bit-for-bit.

---

## 5. What the main session should change, and where

**Do not treat this as closing §2.5 favourably. It closes it, and three of the four things it says cut
against the current text.**

| section | change |
|---|---|
| **§2.5** | Replace `[STILL NOT MEASURED] — a labelled linear probe on every artifact` with the measurement. It is no longer a gap; it is a **result**, and the headline of that result is §1 above, not §3. |
| **§6.2** | Strike the row *"A labelled linear probe on every artifact — Not run"*. Replace with a **CLOSED** row stating: probe floors 1.380×/1.511× on `wsi_biology` against 1.011×/1.012× on the other two views; 0/6 resolvable on the readout view; and that what is *still* absent is a paired bootstrap on the between-arm probe difference and any probe on a second cohort. |
| **§6.3** | Zaiem et al. is now conceded **twice over and with our own numbers**: our two linear estimators disagree on the winner in 2/6, 3/6 and 0/6 pairs across the three views. Add it. Also: the second-largest exposure named there — that the ground truth is one readout — is now **partly answered** (Probe B corroborates it 24/27 on the readout view) and **partly worsened** (Probe A contradicts it wherever it can resolve). Both belong. |
| **§4.1** | The repeat-2 paragraph must gain the probe column. *"One quantity moves by a factor, the other barely moves at all"* is **false as written** and must become "two quantities move by a factor — rank and the labelled linear probe — and the molecular channel does not." |
| **§4.1a** | The floor table gains three rows per view: Probe A logistic, Probe A LDA, Probe B TP53. |
| **§4.1b** | The claim *"rank's usability is conditional on the co-trained view"* survives but loses its implied attribution to rank: the labelled probe is conditional on the same view in the same direction. And *"the view that makes rank usable is the view on which it is most often wrong"* needs "against our unsupervised readout" written into it. |
| **§4.5(c)** | The 3/6 and 2/9 counts are relative to the channel. Against a labelled cancer-type probe on the same two views, rank is 6/6. Add the row; do not delete the existing one. |
| **§4.6 / §4.6a** | §4.6a showed every count moves when the target **block** changes. This adds that every count also moves when the **reference standard** changes — and the standard with the strongest external warrant is the one that disagrees. It is a seventh coordinate system, and it belongs in §4.6a's table logic rather than as a new section. |
| **§4.7** | TP53 picks `programme_free` 3/3 on D1 against a resolvable floor. Report it, and report immediately that ATM and KMT2D pick `programme_only` 3/3 on the same pairs. Do **not** let it stand as a contradiction of §4.7 and do not omit it. |
| `paper/P2_FIGURES.md` | Both rows recording the absence (the "figures the paper does NOT have" table and the status table) must be rewritten to **CLOSED**; the floor table is drawable as a fourth panel of F1. |

**Not touched by this work, deliberately:** `paper/P2_RANK_DRAFT.md` (other agents are reading and
editing it concurrently), `claim_guards.py`, `claim_evidence.json`, and every other agent's
`PREDECLARED_*` file.

---

## 6. What this cannot do

* **One stack, one cohort, one architecture family.** `claim_guards.no_external_cohort` is untouched.
* **The probe floor inherits every limitation §4.1 states of the rank floor**: five repeats, one arm
  (`programme_only`, the *stable* arm), one seed, one configuration, one stack. It is a **floor twice
  over** and is quoted as "floor", never "envelope".
* **No paired bootstrap on the between-arm probe difference.** The predeclared criterion was the
  five-repeat floor, and that is what is used. The bootstrap CIs reported are on the **level** (TP53:
  ~0.07 wide, far wider than any arm delta), which is the same level-versus-paired-difference asymmetry
  §3.5 already forces on the channel. Adding a paired probe bootstrap is the obvious next measurement
  and it is **not** run here.
* **Six pairs cannot support a rate** (§3.6 rule 3). No count above is quoted as one, in either
  direction, including the ones that favour us.
* **Probe A is not a molecular measurement** and is never described as one. The panel endpoints were
  scored on `wsi_biology` only.
* The logistic permutation null is **20 draws** (p floored at 0.0476) because each draw costs ~5 s
  against the LDA null's 0.04 s; the LDA null is 200 draws on every view and both blocks. Where a
  logistic p reads 0.0476 that is the resolution, not the evidence.

---

## Files / commits

- Predeclaration `NOTEBOOK_ENTRIES/PREDECLARED_p2_labelled_linear_probe_20260805T0040Z.md` — `6b3d8e7`
- Runner `v2/research/rebase/p2/p2_labelled_probe.py` — `288a124`, `2a3bd78`, `e772316`, `0aef953`
- Test `v2/tests/test_p2_labelled_probe.py` (8 tests)
- Box: `~/ws_p2probe/` (workspace verified 753/753 byte-equal to HEAD), `out/P2_LABELLED_PROBE.json`,
  `out/shard_s{1..4}.json`, `out/shard_reps_panel.json` (the redundant pass), logs `s{1..4}.log`,
  `panel_reps.log`
- Labels reused, not invented: each artifact's own `cancers` array;
  `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/e1_endpoints/inputs/e1_endpoint_labels.parquet`
  (E1, from the MC3 public MAF)
