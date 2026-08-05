## 2026-08-05 03:30 UTC — §5.4 limit 2 does NOT clear at ten repeats. The pass was a five-repeat draw, one run in ten at m = 0.99 is what breaks it, and the independent second five reaches the same verdict

**Logged:** 2026-08-05 03:30 UTC. **How obtained:** ten further GPU runs on the A100
(`150.136.45.194`) from `~/ws_j2` — the *same* workspace the original fifteen step-600 repeats were
run from, verified file-by-file against a manifest generated from the local checkout at `HEAD`, with
every file on the training path and the scoring path byte-identical after LF normalisation and the
whole `class RankVariant` → `def cca_spectrum` block of `spectral.py` (10,050 characters, both sides)
identical. Launched 00:00:02 UTC, 10 concurrent on an otherwise idle card, finished 02:46:23 UTC.
Scored by `v2/research/rebase/p2/p2_limit2_stress.py`, which imports its statistic table, its fold,
its sha and its shape rule from `p2_envelope_floors.py`, its rank variants from
`calibra.spectral.RANK_VARIANTS`, and its **direction convention** from
`p2_selection_rule.METRICS`. **Predeclared in full before anything ran** at
`NOTEBOOK_ENTRIES/PREDECLARED_probe_floor_n10_and_momentum_grid_20260805T0200Z.md`.

### 0. The awkward finding first

**The predeclared primary falsifier fired.** §5.4 limit 2 — *"the value this project actually
runs"*, `m = 0.999` over `m = 0.99` on the fixed held-out probe at step 600 — **no longer clears its
own floor once the floor is measured from ten same-seed repeats per arm instead of five.**

| | n = 5 (published 2026-08-05 00:00) | **n = 10** | verdict |
|---|---:|---:|---|
| row's ratio (fixed — two specific runs) | 1.262× | 1.262× | — |
| **R3 floor, its own two arms, its own step** | **1.195×** | **1.326×** | — |
| **Test A, `p2_floor_audit.check`'s rule** | **clears by 5.6%** | **DOES NOT CLEAR** | **broken** |
| which arm carries the floor | m = 0.999 | **m = 0.99** | changed hands |

This is not a new criterion and no threshold moved. It is `ratio > floor` — the paper's own rule, the
same checker, the same arms, the same step, the same statistic, the same block, the same seed. The
only thing that changed is that the floor was measured from twice as many runs.

**It is also not one unlucky draw.** Repeats 6–10 scored **alone**, as an independent n = 5 never
pooled with the first five, give an R3 floor of **1.279×** — also above 1.262×. Both halves of the
twenty runs put the floor above the row's ratio. **The original five were the favourable draw.**

And it is not confined to the block §5 quotes: on the `rna_biology` view the n = 10 R3 floor is
**1.276×**, also above 1.262×. Test A fails on both views.

### 1. What breaks it is one run, and that run is not excludable

The entire increase is on the **m = 0.99** arm, and it is a single repeat:

| arm | statistic | reps 1–5 | **reps 1–10** | fold 1–5 | **fold 1–10** |
|---|---|---|---|---:|---:|
| m = 0.999 | R3 | 6.741–8.054 | 6.741–8.054 | 1.195× | **1.195×** |
| m = 0.999 | R1 | 10.353–11.955 | 10.353–11.955 | 1.155× | **1.155×** |
| **m = 0.99** | **R3** | 5.142–5.922 | **4.465**–5.922 | 1.152× | **1.326×** |
| **m = 0.99** | **R1** | 6.462–7.124 | **5.520**–7.124 | 1.102× | **1.291×** |

**The m = 0.999 arm did not move at all** — all five new repeats landed inside the old range. Every
bit of the change is `m0.99 rep6`, which reads **R1 5.520 / R3 4.465** where the other nine span
6.462–7.124 / 5.142–5.922.

**Three reasons it stays in, each checked rather than asserted**, against the exclusion rule written
down before the runs:

1. **It trained.** `biology_contrastive` 7.575 at step 600 against chance `ln 80 = 4.382`, beside
   rep7's 7.853. Its RNA-view mutual cosine (0.9629) sits with the rest of the arm.
2. **It is not collapsed.** R1 5.520 against the m = 0 arm's 1.564–2.738. The predeclared exclusion
   bar was R1 < 4.
3. **The project's own shape rule does not flag it.** `outlier = rep6`, `rest_fold = 1.152`,
   `concordant = False` (1.152 > 1.05), so `bimodal = False`. Under
   `p2_envelope_floors._shape` — imported, not restated — this arm is **not bimodal**, exactly as all
   sixteen probe floors were not.
4. **It is not a batch effect.** m = 0.99 repeats 6–10 span 5.520–7.029 against repeats 1–5's
   6.462–7.124: the two ranges **overlap**, so they are two samples of one spread, not two levels.
   The predeclared distrust condition was non-overlap, and it did not occur. Checked at step 100 too,
   where both arms' batches overlap fully.

So the honest reading is: **same-seed GPU non-determinism at m = 0.99 produces a run about 20% lower
in rank roughly one time in ten, and five repeats did not see it.** The 00:00 entry's description of
this floor — *"the smallest of the sixteen probe floors at 1.195×, and it is carried by m = 0.999"*,
offered as evidence that the rule "measure both sides" returns the stable arm here — is overtaken in
both clauses. At ten repeats it is 1.326× and carried by **m = 0.99**.

### 2. Against the predeclaration, item by item

| predeclared | outcome |
|---|---|
| **B1** — n = 10 R3 floor ≥ 1.262×, the primary falsifier | **FIRED.** 1.326×. |
| **B2** — R1 fails Test B at n = 10 | did **not** fire. R1 still clears, 1.453× against 1.291×. |
| **B3** — ≥ 6 of 11 statistics fail Test A | **NOT EVALUABLE, and that is itself a result** — see §3. |
| **H1** — R1 survives Test B at n = 10 | **holds.** |
| **H2** — R3 still passes Test A | **fails.** |
| **H3** — a second, independent statistic passes Test B | **fails.** `R1_uncentred` and `RankMe` both passed at n = 5 and both **fail** at n = 10. R1 is now alone. |

**The verdict is BREAKS, not HOLDS and not AMBIGUOUS.** Two of the three "holds" conditions failed and
the primary falsifier fired. The one thing that could have broken and did not is R1's ratio test.

### 3. Test A exists only under R3 — the statistic the row fails the strong test on

Predeclaration item B3 asked whether a majority of statistics fail Test A. **It cannot be asked.**
The two runs the row is a claim *about* — `~/e0_run/d1_diag/long_m0.999.log` and `long_m0.99.log` —
were written by the version of `d1_momentum_probe.py` that printed four columns (`eff-rank`,
`feat-std`, `rna-rna`, `contrastive`). **There is no `CANONICAL` column in either file and their
states were never exported**, so no R1 ratio, and no PR / RankMe / stable-rank / α-ReQ ratio, can be
recovered for this row at all.

**So the only statistic under which §5.4 limit 2 can be given the audit's own test is R3 — the one
statistic under which the same ten runs already undercut it.** That is recorded in the module's
`ABSENT` rather than left as a silence. R1's support for this row was never Test A and cannot be
made into Test A without re-running the two runs the row quotes.

### 4. The disagreement is the Hill order, and it is narrower than "R1 says yes, R3 says no"

Every key of `RANK_VARIANTS` was scored beside the published alternatives, and **duplicates were
detected numerically rather than assumed**, because counting a statistic twice would inflate any
"N statistics agree":

| pair | worst relative difference over all 20 states | why |
|---|---:|---|
| R3 == R2 | 3.8e-08 | `z_biology` is L2-normalised at output, so `normalise_rows` is a no-op |
| R1_rownorm == R1 | 5.9e-08 | same reason |
| PR_rownorm == PR | 7.1e-09 | same reason |
| RankMe ≈ R1_uncentred | **3.0e-04** | both are exp(Shannon entropy) of the *uncentred* spectrum; they differ only by the `eps` RankMe adds **outside** the normalisation |

Thirteen labels are **seven** statistics. At n = 5 the ratio test split them exactly by **Hill order** —
every order-1 variant passed, every order-2 variant failed — across both centring settings and both
row-normalisation settings, so the split tracked the order and nothing else. At n = 10 only R1
survives.

**And α-ReQ needed its sign.** `|α − 1|` is smaller-is-better (`p2_selection_rule.METRICS` gives it
`direction = -1`), and grading it on raw values printed its **perfectly ordered** arms as INVERTED — a
sign convention that reads as a statistic contradicting the others. The direction table is imported
from the selection rule rather than restated, and a test asserts they agree.

### 5. What survives, and it survives more strongly: the arms are completely separated

Test B and Test A are both statements about a **ratio**. Test C asks about **order**: is every
`m = 0.999` repeat above every `m = 0.99` repeat? Unlike Test B — whose two sides both move against
the pass as repeats accumulate — **Test C gets more surprising the more repeats it survives.**

**At n = 10, all ten m = 0.999 repeats are above all ten m = 0.99 repeats**, under R1, R2, R3,
R1_rownorm, R1_uncentred, RankMe and (with its sign applied) α-ReQ `|α − 1|`. Exact one-sided
permutation probability of that arrangement under exchangeability of the twenty runs:
**1/184,756 = 5.4e-06**, up from 1/252 at five per arm.

**Scope, stated in the same words it will be reported in:** the repeats differ only in GPU
non-determinism — same seed, same workspace, same stack, one card — so exchangeability holds over
*that* noise source and nothing else. **This is a statement about run-to-run reproducibility and is
NOT a p-value for the momentum effect**, which would need the seed varied, and §4.2 measures the seed
as the dominant term.

Two statistics do **not** order the arms even at n = 10, and they are recorded rather than dropped:
**PR** (worst-case ratio 0.754×, ranges overlapping) and **stable rank** (0.724×, overlapping). The
uncentred order-2 variant `R2_uncentred` **loses** its n = 5 separation at n = 10.

### 6. The sharpest thing these runs say, and it cost no GPU: the verdict flips with the reading step

The row is quoted at step 600. `d1_momentum_probe.py` exports the probe states at **every** step it
reads, so the same two arms and the same repeats can be given the same three tests at 100–600 for
nothing. **At n = 5 the ratio verdict was already unstable; at n = 10 Test A fails almost everywhere.**

**n = 5** (B = worst-case ratio against that step's floor; A = the audit's rule):

| step | B: R1 | B: R3 | B: RankMe | **A (R3)** | C, all variants |
|---:|---|---|---|---|---|
| 100 | FAIL | FAIL | PASS | PASS | SEP |
| 200 | FAIL | FAIL | FAIL | **FAIL** | SEP |
| 300 | PASS | PASS | PASS | PASS | SEP |
| 400 | PASS | **FAIL** | FAIL | PASS | SEP |
| 500 | PASS | PASS | FAIL | PASS | SEP |
| 600 | PASS | **FAIL** | PASS | PASS | SEP |

R3's ratio verdict **passes at two of six readings and changes its answer four times in six
consecutive steps**. RankMe changes four times too. Canonical R1 changes once and settles from 300.

**n = 10:**

| step | B: R1 | B: R3 | B: RankMe | **A (R3)** |
|---:|---|---|---|---|
| 100 | FAIL | FAIL | PASS | PASS |
| 200 | FAIL | FAIL | FAIL | **FAIL** |
| 300 | FAIL | FAIL | FAIL | **FAIL** |
| 400 | PASS | FAIL | FAIL | **FAIL** |
| 500 | PASS | FAIL | FAIL | **FAIL** |
| 600 | PASS | FAIL | FAIL | **FAIL** |

**Test A now fails at five of the six steps the states were saved at**, passing only at step 100 where
both arms are still near their initialisation. Test C is SEP at every step at both repeat counts.

**So the statistics never disagreed about which arm is higher. They disagreed about whether the gap
clears the noise — and that answer depends on the reading step about as much as on the statistic.**

### 7. The arithmetic that decides what any of this is worth, and it cuts both ways

The floor is `max/min` over an arm's repeats, so it is **non-decreasing** in the repeat count; the
Test B separation is `min(high)/max(low)`, so it is **non-increasing**. **Adding repeats can only ever
make Test B harder.** Consequences, predeclared before the runs and honoured here:

1. **"R3's Test B got worse at n = 10" is not reported as a finding.** It could not have got better.
   A test asserts the monotonicity directly against the two files, so if the n = 5 and n = 10 outputs
   are ever not nested samples of one experiment the suite says so.
2. **"R1's Test B survived at n = 10" IS a finding** — it had ten repeats per arm to break it and did
   not.
3. **Test A is the one that can genuinely flip**, because the ratio is fixed while the floor grows.
   **It flipped.** That is why B1 was written as the primary falsifier and why the break is real
   rather than an artifact of asking a harder question.

### 8. What the paper should now say — flagged, NOT edited

`paper/P2_RANK_DRAFT.md` is being edited concurrently by a completeness-audit agent and has not been
touched by this work. The following are the exact locations that are now wrong.

1. **`floor_audit.json`, row `5.4-m0999-over-m099`** — `"clears": true` is **false at n = 10**. The
   floor entry `R3_probe_step600_m0999_vs_m099` reads `1.1947` and the measured value at ten repeats
   is **1.3263**. Its `rests_on` text (*"a pass by 5.6%, and the ten runs that measured the floor
   undercut it"*) should become "does not clear once its floor is measured from ten repeats per arm
   rather than five". **This file is owned by another agent and was not modified.**
2. **§4.1a's audit table and counting sentence** — this selection moves from *clear* to *do not
   clear*. The counts become **11 clear / 14 do not clear** of 25 if nothing else changes. The blocks
   are generated by `p2_floor_audit.py`, so regenerating is the fix, not hand-editing.
3. **§5.4 limit 2 prose** — *"The rule now licenses it — by 5.6%"* is withdrawn. So is the 00:00
   entry's §4 framing. The replacement is §0–§1 above.
4. **§"Status" item 12** — the fragility recorded there resolved, and it resolved **against** the
   pass. It should say the row was pushed to ten repeats and failed.
5. **The 2026-08-05 00:00 entry §2** — *"m = 0.999 against m = 0.99 has no collapsed arm — both train
   — so its floor is the smallest of the sixteen probe floors at 1.195×, and it is carried by
   m = 0.999. The rule is not 'the collapsed arm is noisier'. It is 'measure both sides', and here
   that returns the stable one."* **The last clause is now wrong**: at ten repeats the m = 0.99 arm is
   the noisier one and carries the floor. The rule "measure both sides" is *vindicated*, but the
   example chosen to illustrate it inverts.
6. **§3.1** should still record that R2 and R3 are the same statistic on this block; §4 above adds
   that R1 = R1_rownorm and PR = PR_rownorm for the same reason, and that RankMe is R1_uncentred plus
   its own eps — so §3.1's "three statistics travel under one name" is, on this block, **seven labels
   for four**.
7. **A new limitation worth stating explicitly**: §5.4 limit 2 can only be judged under **R3**,
   because the two runs it quotes predate the canonical column and the state export. The statistic
   the row fails the strong test on is the only one that can rule on it.

### 9. What is still open

The momentum grid (five same-seed repeats each at **m = 0.995** and **m = 0.98**, so the shipped
comparison sits inside {0, 0.98, 0.99, 0.995, 0.999} rather than being a two-point contrast) was
launched at 02:46 UTC as wave 2 and is running. Its predeclared reading rule — STEP / SMOOTH /
NEITHER — is §6 of the predeclaration. **It does not bear on §0: the row is already broken at its own
two arms.** Note the draft's §5.2 table is already non-monotone at this step (m = 0.9 reads 2.23
against m = 0's 2.81), which the draft's *"the effect is monotone in m"* does not accommodate.

### 10. Files

* `v2/research/rebase/p2/p2_limit2_stress.py` — Tests A, B, C at any repeat count, any step, any view,
  any pair or grid; duplicates detected; direction imported from `p2_selection_rule`.
* `v2/tests/test_p2_limit2_stress.py` — 19 tests, including the monotonicity assertion that makes the
  "what is this worth" labelling true, and the n = 10 break.
* `~/e0_run/d1_probefloor600/out/P2_LIMIT2_STRESS_{N5,N10,LATE5,N10_RNA}.json` and
  `limit2_stress_run.log`, vendored under `v2/research/rebase/p2/figures/data/`; each carries the
  sha256 of every probe state it read and its own `absent` list.
* The ten new run logs, `pf600_m{0.999,0.99}_rep{6..10}.log`.
