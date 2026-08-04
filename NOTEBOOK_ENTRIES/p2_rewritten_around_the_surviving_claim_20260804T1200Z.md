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
