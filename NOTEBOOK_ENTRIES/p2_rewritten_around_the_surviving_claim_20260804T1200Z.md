## 2026-08-04 12:00 UTC — P2 rewritten around the claim that survived; the old claim is withdrawn, and P1 is de-duplicated

**Logged:** 2026-08-04 12:00 UTC. **How obtained:** writing and restructuring only. No experiment was
run, no GPU was touched, no artifact was read. Every number in the rewritten draft is quoted from
`NOTEBOOK_ENTRIES/p2_competing_metrics_and_necessity_test_20260803T2326Z.md`,
`p2_prior_art_citation_graph_sweep_20260803T2326Z.md`,
`effective_rank_canonicalised_and_every_instance_recomputed_20260804T0005Z.md`,
`rank_probe_repeat_variance_20260804T0900Z.md`, `g26_is_not_reproducible_20260804T0700Z.md`,
`PREDECLARED_D1_necessity_test_20260803T2300Z.md`, or the evidence files those entries cite. Suite
re-run to confirm green: **317 passed in 50.9 s**, thread-capped
(`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=1`), across
`v2/tests` + `tests`. Files touched are markdown only.

---

### 0. Bad news first

1. **The claim P2 was built on is dead and has been deleted from the draft.** "Effective rank does not
   track information content" is not supported by our own best-matched three-seed experiment. The
   necessity test confirmed RankMe 3/3 (`programme_only`, the higher-rank arm, carries the larger
   channel in every seed). Under `PREDECLARED_D1_necessity_test_20260803T2300Z.md` this is outcome
   **O2**, whose predeclared reading is "report as a limitation with the same prominence a confirmation
   would have had". It is now abstract paragraph 1 and §4.7, placed ahead of every result that favours
   us.
2. **§2.2's central sentence was false and is deleted.** Aldeneh et al. (ICASSP 2025) publish a
   within-method matched-arm selection-rule failure *and* a low-rank/high-information instance. §2.2
   now leads with it. The novelty claim is re-scoped to the within-method matched-arm evaluation in
   pathology/transcriptomics plus the reproducibility floor and the variance decomposition.
3. **§5 (the integrated liveness-gate and queue material) contains a real tension with §4.1 and the
   draft says so.** The momentum decision at step 600 is 2.81 against 7.42 — a **2.64× ratio, inside
   the 2.69× retraining envelope §4.1 uses to disqualify six of seven arm comparisons** — and the
   sweep is one seed per momentum value. §5.3 states this, gives the three reasons the fix is still
   reported as real (monotone across four m values, flat over 400 steps, readings at the collapse
   floor with cosine corroboration, tight fixed-seed repeats with an empty band), and states that
   **no seed replication of the momentum sweep exists**. That measurement is now named as the one that
   would close the gap between §4 and §5.
4. **The D1 channel numbers still have no intervals.** `D1_PAIRED_BOOTSTRAP_STRATIFIED.json` does not
   exist — a stale absolute path in the D1 audit chain. Those cells carry
   `[D1 PAIRED BOOTSTRAP PENDING]`. The one-line path bug is still unfixed.

---

### 1. The claim, before and after

| | before | after |
|---|---|---|
| headline | effective rank does not track information content | **effective rank is unusable as a selection signal, because its between-arm differences are smaller than its own within-arm reproducibility floor — inside the regime its proponents reserve for it** |
| carried by | six "dissociations", four of which do not contradict RankMe as written | §4.1 (6 of 7 ratios, 1.004×–2.19×, inside a 2.69× envelope); §4.2 (34.5% arm / 65.5% seed, F(3,8)=1.41 n.s., against the channel's 98.0% arm, F=128.2); §4.3 (6.05× vs 1.003× at matched step); §4.4 (the defeater check) |
| survives the necessity test? | no — the test was the claim's crux | **yes** — rank can be right on average and still be unusable for a single comparison, because the ratios it got right (1.74×–3.25×) are themselves inside the within-arm seed range (2.10×–3.75×) |
| depends on a sign count? | yes | **no** — §3.6 rule 3 states that at n=6 a flawless 6/6 gives p=0.031, and no count in the paper carries weight |

### 2. Section structure of the rewritten draft

Modelled on Leavitt & Morcos (ICLR 2021) — claim grammar as a property rather than a performance
result, and a defeater section in the same position as their §4.2.

1 Introduction (1.1 the appeal, taken seriously; 1.2 prior art and RankMe's hedges; **1.3 the claim**;
1.4 contributions; 1.5 relationship to P1) · 2 Related work (2.1 the proposal + canonical definition +
the ε; **2.2 prior negatives, led by Aldeneh**; **2.3 the strongest defences, new**; 2.4 collapse
literature; 2.5 alternatives computed on our artifacts; 2.6 verification table) · 3 Methods (3.1
canonical definition + ten call sites + the fourth degree of freedom; 3.2 measures and nulls; 3.3
cohorts and the 12 artifacts; 3.4 matching; 3.5 seed non-reproducibility; 3.6 reporting rules) ·
**4 Results (4.1 the envelope; 4.2 variance decomposition; 4.3 the floor belongs to the arm; 4.4
defeater check; 4.5 verdict instability; 4.6 selection rules, underpowered; 4.7 the necessity test,
which went against us; 4.8 dose–response; 4.9 historical instances; 4.10 the surviving use)** ·
**5 Worked example: gate + queue, graded in rank (5.1, 5.2, 5.3 where it strains)** · 6 Limitations as
named objections · 7 Conclusion · Appendices A (provenance), B (code), C (per-number caveats).

2,029 lines, up from 1,451. Every table carries a provenance line; every rank number names its
statistic **and** its block.

### 3. The seven corrections, and their status

| # | correction | status |
|---|---|---|
| 1 | §2.2's central sentence is false; lead with Aldeneh et al. (ICASSP 2025, DOI 10.1109/ICASSP49660.2025.10889651, arXiv:2409.10787); re-scope the contribution | **DONE.** Sentence deleted, Aldeneh leads §2.2 with the verbatim quote and an explicit "what is still ours after this" paragraph that calls the remaining distinction narrow |
| 2 | cite the defences (Deng arXiv:2510.10948; Awasthi *BMC Genomics* 26(1):710, DOI 10.1186/s12864-025-11913-2) and argue against the strongest form | **DONE.** New §2.3 with both in their own words plus Ruan and six further defences; the Awasthi rebuttal is confined to "it is cross-method, which RankMe disclaims", and §2.3 closes by conceding that the published cases where rank works are predominantly cross-method/cross-capacity |
| 3 | §4.2 dilution: −3.10%, 21.5×, range 1.95×–21.5×; note the correction moved it in our favour | **DONE** (now §4.8), including the sentence that the flattering number is not used anywhere |
| 4 | §4.3 seed 43 restated as implementation-dependent | **DONE** (now §4.9): "the rank ordering is wrong in 1 of 3 seeds under every combination, and never the same seed" |
| 5 | state the canonical definition explicitly, and RankMe's ε outside the division | **DONE** (§3.1, §2.1): order 1, column-centred, rows at own norms, LAPACK relative cut; RankMe's `p_k` sum to `1+min(N,K)ε` and its statistic is not the exponential of a Shannon entropy |
| 6 | selection-rule counts must not carry rhetorical weight | **DONE.** §3.6 rule 3 and §4.6 both state p=0.031 for a flawless 6/6, and §4.6 opens with "it is underpowered and must not carry rhetorical weight" |
| 7 | mark `[NOT RECOMPUTABLE]` instances rather than carrying numbers forward | **DONE.** Decorrelation and "16/16" are `[NOT RECOMPUTABLE — artifact never existed]`; D1-A's 9.81/1.71 is `[NOT RECOMPUTED]`; all three are excluded from every count |

**None of the seven could not be made.** Two are partial in a way worth recording: correction 2's
rebuttal of Deng et al. rests on an argument (they sweep capacity-like variables, we match on capacity)
and **not on a measurement** — we have not measured rank at capacity scale, and §6.2 says so. And
correction 4 leaves the D2 instance materially weaker than it was quoted as being, which is why §4.9
rather than §4.1 now carries it.

### 4. Structural requirements

- **Liveness gate + queue anchoring integrated as §5**, per the author's decision, framed as a worked
  example of rank used in the regime §4.10 defends. `paper/LIVENESS_GATE_DESIGN.md` and
  `paper/QUEUE_ANCHORING.md` are retained as working sources and both now carry a **SUPERSEDED AS A
  SUBMISSION UNIT** header pointing at P2 §5, so neither can be submitted separately and reuse P2's
  results. The three MoCo corrections §2.4 of the old draft demanded (identifier; advanced twice as a
  hypothesis; queue-specific scope) are made in the integrated text **and** in both source files.
- **P1 de-duplication — all four edits executed.** (i) `P1_FIGURES.md` F11 replaced by a `DELETED`
  stub stating both reasons; `P1_CALIBRA_DRAFT.md` §4.11 reduced to a two-paragraph pointer with no
  table and no rank numbers. (ii) The instance-3 description is gone with the table, and the pointer
  states explicitly that it was a hard `matrix_rank` at a ceiling of 16 and that the centred effective
  rank falls 12.88 → 1.00. (iii) P1 §2.6 rewritten: RankMe named as the canonical proposal with its
  full verified citation, Roy & Vetterli credited with the statistic and explicitly *not* with the
  quality claim, Jing et al. mis-grouping removed and the correction recorded; §2.7's
  `[CITATION NEEDED]` for §2.6 closed and the remaining two `[CITATION NEEDED]` items restated. (iv)
  `P1_FIGURES.md` F10(b) keeps the twin axis but drops the rank-versus-information claim and must now
  name the block each curve is computed on. P1 §4.12(iv) retained per the deconfliction note.
- `paper/P2_FIGURES.md` carries a **STALE** banner listing what the rewrite requires: a new headline
  figure built from §4.1's ratio table and §4.2's variance decomposition (the single most important
  display item, previously with no row at all), a required defeater-check figure, a required
  verdict-instability figure, F2 demoted, a necessity-test figure placed before the favourable ones,
  and F8 no longer pending on the rank side.

### 5. Everything in the evidence that still contradicts itself

Recorded because the report asked for it and because none of these is resolved by writing.

1. **§5.2's momentum ratio is inside §4.1's disqualifying envelope** (2.64× against 2.69×), single
   seed. Stated in §5.3 and in the `QUEUE_ANCHORING.md` header. Unresolved; needs a seed replication.
2. **D1 is 3/3 under R1 and R3 but 2/3 under R2**, the statistic an earlier draft nominated for that
   exact table. The draft quotes the honest form ("2–3 of 3 depending on which function you call
   effective rank") rather than a clean 3/3.
3. **The two P2 notebook entries' R3 rows disagree** because one is computed on the raw block and the
   other on the residualised block. Neither is wrong; the block was not stated. §4.5 names the block in
   every row and the draft cites both entries with the reconciliation.
4. **Two collinearity values for the same-sounding quantity** — 0.7362 (std 0.0314) on the full cohort
   versus 0.80 on the 16-patient gate batch — are not the same measurement. Appendix C says so; nothing
   depends on which is used.
5. **The 2.69× envelope rests on one retraining pair.** §4.2 exists because of that and reaches the
   same conclusion from 8 within-arm degrees of freedom, but §4.1's headline number is thin and §6.2
   now names the controlled repeat design as **the** most valuable missing measurement.
6. **The 2/9 view-restricted count in §4.5c inherits the `rna_biology` circularity caveat** and is not
   quoted as a rate.
7. **D1's preregistered escalation is flagged and unresolved.** `D1_PAIR_MANIFEST.json` says "if
   programme_only wins, the collapse story is wrong — escalate, do not proceed to D2". It won 3/3. This
   paper flags it and does not resolve it.

### 6. Files / commits

- `paper/P2_RANK_DRAFT.md` — rewritten (1,451 → 2,029 lines)
- `paper/P1_CALIBRA_DRAFT.md` — §2.6, §2.7, §4.11
- `paper/P1_FIGURES.md` — F10(b) caption, F11 deleted
- `paper/P2_FIGURES.md` — stale banner, deconfliction marked executed
- `paper/LIVENESS_GATE_DESIGN.md`, `paper/QUEUE_ANCHORING.md` — superseded headers + MoCo corrections
- Suite: **317 passed**, thread-capped. No source file changed.

### 7. Addendum — the workspace-drift audit, landed concurrently

`WORKSPACE_DRIFT_AUDIT_ALL_20260803T2359Z.md` (commit `0d6eed8`) arrived on the branch between this
work's two pushes and bears directly on §4.2 of the draft. Two things folded in:

1. **The §4.2 numbers were computed by a `spectral.py` predating the canonicalisation** — no
   `CANONICAL`, no `RANK_VARIANTS`. Recomputed against a verified-current workspace, all five reproduce
   exactly (34.5%, F 1.41, 29.1%, 98.0%, F 128.20). They reproduce because the consolidation added
   variants without moving the default and `residualise.py` is byte-identical everywhere. **That is
   luck, not design**, and §4.2 now says so, and names the object-identity test as the reason it cannot
   recur.
2. **`~/e0_run/p2_*.py` are not vendored into the repository**, unlike the rank-recomputation scripts.
   §6.2 now carries a row making vendoring them a pre-submission requirement, because until then the
   "every number traces to a file in the repository" rule is unsatisfied for §4.2, §4.4, §4.5 and §4.6.

Nothing in the audit changes a number in the draft.

---

### 8. Addendum 2 — the D1 bootstrap exists, A3's verdict is qualified, and our own margin is now in the argument

Four coordinator updates folded into `paper/P2_RANK_DRAFT.md` after the rewrite.

**1. `[D1 PAIRED BOOTSTRAP PENDING]` is closed — and the answer is less favourable than the point
estimates were.** The stratified bootstrap existed all along; the audit chain's stale absolute path hid
it. §4.7.2 now carries both estimators, in the predeclared direction
`Δ = channel(programme_free) − channel(programme_only)`:

| seed | Δ | patient CI₉₅ | cancer-cluster CI₉₅ |
|---|---:|---|---|
| 42 | −0.0705 | [−0.0938, −0.0444] | [−0.0957, −0.0180] |
| 43 | −0.0863 | [−0.1186, −0.0522] | **[−0.1386, +0.0006]** |
| 44 | −0.0961 | [−0.1314, −0.0618] | [−0.1535, −0.0016] |

**Decisive 3/3 on the patient bootstrap, 2/3 on the cancer-cluster bootstrap.** The cluster estimator
resamples whole cancer types and is the conservative one, so it is the one weighted — and the draft says
so in the abstract, the short abstract, the status block, §4.7.2, §4.7.3, §6.2, the conclusion and
Appendix C. The reason for that redundancy is stated in §4.7.2: **§4.6 refuses to let a 5/6 sign count
carry weight, so quoting the favourable one of two estimators on the paper's most load-bearing negative
would be incoherent.** Note the sign convention changed: the table previously used `P − F` (positive);
it now uses the predeclaration's `F − P` so the CIs and the Δ they bound agree.

**2. Audit check A3 is recorded with its qualification, not as a pass.** It passes **on arm difference**
(`random_control` gaps −0.022 / −0.007 / −0.032, all CIs spanning zero). But the **absolute** level is
0.44–0.48 against real targets' 0.51–0.62 and a within-cancer permutation null of **0.147** — random
gene sets carry roughly **three times the floor**, agreeing with T1.4's independent 76–82%. The verdict
now reads *"the instrument does not manufacture an arm difference; the absolute level is high and
separately explained"* in §4.7.3, §6.3 and Appendix C. It does not weaken D1, which is a paired arm
difference — the exact quantity the control clears — but it forbids reading any absolute channel level
against an assumed null of zero.

**3. Our own margin is now in the argument rather than waiting to be used against us.** The smallest arm
difference we quote (0.0705) is **about half** the real-versus-random-control margin (0.139). That is
§4.1's envelope argument arriving from the channel side and pointing at us. New paragraph in §4.7.4,
promoted into the short abstract and made §6.3's second exposure. It sets out the three respects in
which our paired, interval-backed differences sit differently from rank's uninterval-backed levels
(paired-vs-level, intervals-vs-none, and the variance decomposition depending on no margin at all) and
**states that none of the three makes our own margin large.** This pre-empts the sharpest available
objection — that we hold rank to a strict standard and ourselves to a loose one.

**4. §5.2's single seed was a defect, not a design choice.** The momentum harness had its **seed
hardcoded**; the sweep could not have varied seeds had it been asked to. §5.3 and the
`QUEUE_ANCHORING.md` header now say so, and both record that a seed-replicated sweep is **armed**
(`m ∈ {0, 0.999} × 3 seeds`, canonical statistic alongside the participation ratio). §5.3 states the
disjunction in advance: **if the distributions separate across seeds the tension with §4.1 disappears;
if they overlap, the momentum choice is a rank comparison this paper's own rule disqualifies and §5.2
must be rewritten to rest on downstream behaviour rather than on rank.** Committed before the result so
the reading cannot be chosen afterwards.

**Still open after this addendum:** the stale path in the D1 audit chain is unfixed (it hid a result for
a day); the seed-replicated momentum sweep has not reported; `~/e0_run/p2_*.py` are still unvendored.

Suite re-run after these edits: **317 passed**. Markdown only.

---

### 9. Addendum 3 — the vendoring caught a mislabelled statistic in our own analysis code

Commit `7b37dce` (concurrent, another agent) vendored the five P2 analysis scripts to
`v2/research/rebase/p2/` with an end-to-end test, closing the §6.2 row this rewrite had opened. Two
consequences for the draft, one of them a correction that moves a published number.

**What reproduced.** Re-run from a checkout verified byte-equal to HEAD (402/402 files by git blob
SHA-1): **§4.2, §4.4(1), §4.5(b), §4.5(c), §4.6 and §4.7 reproduce to every published digit.**

**What did not.** `p2_rank_variants.py` — the script behind **§4.5(a)** — began
`sys.path.insert(0, "/home/ubuntu/ws")`, the workspace the drift audit found most stale, and carried
its own inline `R1`/`R2`/`R3`. **Its "R2" and "R3" were not R2 and R3.** They were `(Σσ²)²/Σσ⁴`, the
order-2 Hill number of the **eigenvalue** distribution, where `d1_audit.py`'s R2 — and therefore
`RANK_VARIANTS["R2"]` — is `(Σσ)²/Σσ²`, the order-2 Hill number of the **singular-value**
distribution. Different statistics.

Folded into the draft as follows:

- §4.5(a)'s rows 2 and 3 are **relabelled `PR` and `PR_rownorm`**, which is what they are; the numbers
  are unchanged and reproduce cell for cell.
- **The headline count falls from 3 of 6 pairs to 2 of 6** (D2 s43, D2 s44). Propagated to the
  abstract, §1.4 contribution 5 and §4.4's forward reference.
- **§4.7.3's "D1 is only 2/3 under R2" qualification is WITHDRAWN.** Canonical R2 scores D1 **3/3**.
  This *removes* a qualification from the necessity result, i.e. **it costs us** — D1 now confirms
  RankMe under every canonical statistic we compute, and the only surviving qualification on D1 is the
  interval one (3/3 patient, 2/3 cluster). Propagated to the status block and Appendix C.
- §4.5(b)'s R2/R3 **levels** are unaffected; they were computed with the canonical function all along.

**The slot is closed.** `NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260804T0255Z.md` landed while
this was being written and carries the corrected per-cell rows, so §4.5(a) now tabulates all five
statistics — canonical R1/R2/R3 plus PR/PR_rownorm under their true names — rather than carrying a
pending marker. Canonical R2 reads `OK OK MISS OK OK OK` (5/6) and canonical R3 reads
`OK MISS OK OK OK OK` (5/6, identical to R1). Nothing was reconstructed by hand.

**A second, smaller error the same entry caught.** §4.5's provenance note had explained the two source
entries' differing R3 rows as raw-versus-residualised block. **That was wrong** — both are
residualised; the difference was the statistic. Corrected in §4.5(a). The R3 *levels* in §4.5(b) were
always canonical and reproduce exactly.

**And a duplication now flagged:** the `PR` row is identical cell for cell to §4.6's "participation
ratio" row. That identity is the arithmetic signature of the substitution, and §4.5(a) now says the two
must not be counted as independent evidence.

**The point worth keeping.** This is §3.1's finding — three functions under one name — recurring in
*our own analysis code*, and it survived until the "every number traces to a file in the repository"
rule was actually enforced rather than merely stated. §4.5(a) now says so at the same prominence as
the finding it qualifies.

**Working-tree note.** Commit `7b37dce`'s files carried further uncommitted modifications by a
concurrent agent when this addendum was written. Only `paper/P2_RANK_DRAFT.md`,
`paper/QUEUE_ANCHORING.md` and this entry were staged; nothing under `v2/` was touched or committed by
this work.

---

### 10. Addendum 4 — the null comparator was wrong, and the mislabelled statistic is promoted into the argument

**1. §4.7.3 compared D2's random controls against the wrong permutation null, and it is the exact error
this paper is about.** The comparator used was **0.147**, the dilution sweep's **within-cancer**
permutation at a different *n*. The correct comparator for D2/D1 is **0.140**, a **200-draw row-shuffle
of the residualised target matrix** (`NOTEBOOK_ENTRIES/D2_stratified_result_20260803T1210Z.md`).
Corrected in §4.7.3, §6.3 and Appendix C. **The conclusion is unchanged**: 0.4425–0.4810 against 0.140
is **≈3.2×** the floor (3.1–3.4× across the six control readings), where 0.147 gave 3.0×.

A **footnote added to §3.2** now records that this project quotes **at least three** permutation nulls —
**0.140** (D2 row-shuffle), **0.145–0.147** (dilution within-cancer), **0.151–0.158** elsewhere — that
they differ in *n*, component count and, for 0.140, **in the permutation procedure itself**, and that
P1's own audit says they **must not be carried across**
(`p1_submission_draft_20260803T1230Z.md` §5). The footnote records the error rather than only the fix,
because it is §3.1's failure mode in a second quantity.

**And it flags a label conflict inside our own notes**, unresolved: `d1_a3_verdict_and_effect_vs_floor`
calls 0.140 a *"within-cancer"* null where `p1_submission_draft` identifies it as a **row-shuffle**.
**The value is agreed; the procedure label is not.** The draft follows the latter and flags the former
rather than quietly picking one.

**2. The mislabelled statistic is promoted from erratum to evidence.** §4.5(a) now carries a named
subsection arguing that a *fourth* statistic under the name `effective_rank`, inside the analysis code
for the section that argues the name is unreliable, is the paper's thesis demonstrating itself. The
framing does the work the coordinator asked for: **the substitution was invisible to review, invisible
to the test suite, and invisible to the authors across two drafts and a full recomputation pass, and
was found only when the traceability rule was enforced by vendoring the producing code.** The
recommendation is therefore mechanical rather than moral: *against this class of error the defence is
mechanical provenance, not diligence.* Two supporting details are included so it reads as usable rather
than confessional — the **`PR` row is cell-for-cell identical to §4.6's participation-ratio row**, which
is the arithmetic fingerprint a reader can look for elsewhere (two differently-named rows agreeing
exactly are one statistic reported twice); and the correction moved one count **against** us and
withdrew one qualification **for** us, **the scatter signature of an error rather than a bias**.
Promoted to §1.4 as contribution 9 and to a new paragraph in §7.

**3. Sign convention verified end to end.** §4.6's table uses `Δ = A − B` (better arm first, all
positive) and §4.7.2 uses the predeclared `Δ = F − P` (negative). Both are correct for their table and
the two are now **explicitly reconciled in §4.6's caption**, which states that the same three D1 gaps
read −0.0705 / −0.0863 / −0.0961 in §4.7.2. Every other mention checked: §1 abstract, §4.7.1's
predeclared threshold, §4.7.2's narrative (−0.0863), §4.7.4, §6.3 and Appendix C use the magnitude
unsigned or the negative convention consistently. No figure caption in `P2_FIGURES.md` quotes a D1
delta.

Suite after these edits: **326 passed**, thread-capped. Markdown only; nothing under `v2/` touched.
