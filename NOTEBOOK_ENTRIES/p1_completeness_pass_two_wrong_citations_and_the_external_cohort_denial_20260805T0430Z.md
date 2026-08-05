# P1 completeness / QA pass — two citations were wrong, and the paper was still denying an external cohort that exists

**2026-08-05 04:30 UTC.** Outcome: **FIX** (no measurement, no compute beyond the test suite and
bibliographic API calls).

**Predeclaration status: none required, and none was made.** Nothing here is a new number. Every
value written into `paper/P1_CALIBRA_DRAFT.md` or `paper/P1_FIGURES.md` in this pass was already on
disk in a notebook entry, a result file or a CSV ledger, and is re-read from its own recorded source.
The two counts I did compute — the `GATE_LOG.md` row census and the six raw ÷ adjusted ratios of
§4.2's site table — are arithmetic over data printed in the draft itself, and both are shown below so
they can be re-derived rather than trusted.

**Scope discipline.** `v2/calibra/claim_guards.py`, `v2/research/rebase/nature/claim_evidence.json`,
other agents' `PREDECLARED_*` files, `paper/P2_RANK_DRAFT.md`, `v2/attributable_basis.py` and
`v2/tests/test_attributable_basis.py` (another agent's in-flight work, uncommitted at the time) were
**not edited**. The GPU was never touched; this was CPU and HTTP only.

---

## 0. The awkward findings first

### 0.1 Two citations in P1 were wrong, and both were author-name errors on real papers

This is the failure mode a title-and-year spot-check cannot catch, and it is the closest thing yet to
the project's three historical fabrications.

| where | as it stood | what the record says |
|---|---|---|
| §2.5, the random-signature negative control | **"Venet, Dhanasekaran & Sotiriou** (*PLoS Comput. Biol.* 2011)" | Venet D, **Dumont JE, Detours V**, *PLoS Comput. Biol.* 7(10):e1002240, DOI 10.1371/journal.pcbi.1002240. Title, journal and year were correct. **Two of the three author names did not belong to this paper.** |
| §2.5, ComBat | "Johnson, Li & **Rabinowitz** 2007" | Johnson WE, Li C, **Rabinovic A**, *Biostatistics* 8(1):118–127, DOI 10.1093/biostatistics/kxj037 |

Both are now corrected in the draft with a visible `*(Correction, 2026-08-05: …)*` note beside each,
per the append-don't-rewrite convention. §2.7 records both in its status table and §5 limitation 15
now leads with them rather than with a generic "reference verification is incomplete".

**The transferable lesson, written into §2.7 so the next pass inherits it:** *check author lists, not
just identifiers.* Every fabricated-or-wrong citation this project has produced has had a plausible
title and a real venue. An arXiv-ID sweep of the kind P2's audit ran would have cleared both of these
without noticing anything, because neither is an arXiv paper and neither identifier is fake.

### 0.2 P1 asserted, in its abstract, that no external cohort had been through the instrument — seven hours after one had

`paper/P1_CALIBRA_DRAFT.md` said, in three places (abstract, §1.3, §5.1):

> No external cohort has been through the instrument; every number here is TCGA.

`NOTEBOOK_ENTRIES/alchemist_external_replication_RESULT_20260804T2115Z.md` (2026-08-04 21:15 UTC)
records 1,106 paired ALCHEMIST-ALCH NSCLC patients measured through the same instrument, channel
**R = 1.110** at matched *n* against the 841-patient TCGA-NSCLC comparator, *p* = 0.0033 at the 1/301
resolution floor, on the primary 59-target block. `PROJECT_GUIDE.md` §1 already describes P1 as
"replicated in an external cohort (ALCHEMIST)". The paper did not.

**This is a "closed state still written as open" of the same shape P2's completeness pass found, but
pointing the other way: the paper was under-claiming to the point of stating something false.**

The correction is deliberately narrow, because the *reason* the `no_external_cohort` blocker stays up
is real and is not the same as the sentence being true:

* What is still true, and is what the three passages now say: **no injection, transmission floor,
  detection floor or attenuation slope has ever been computed outside TCGA.** Every number in P1
  remains TCGA. The floors are what this paper is about.
* What was false: that no external cohort had been *through the instrument*. One has; it replicated
  the **channel**, which is a different quantity from the floors.
* Why the blocker nonetheless stays undischarged, now stated in the paper (§1.3, §5.1) rather than
  only in `NOTEBOOK_ENTRIES/decision_external_cohort_blocker_stays_20260804T2200Z.md`: the blocker
  gates only the two **per-axis** claim kinds, and an aggregate-channel replication is evidence at
  the wrong granularity to discharge a per-axis claim. `_is_discharged` counts external cohorts and
  cannot tell the difference — recorded there as a latent defect in the guard.
* The cohort classifier line that must travel with any external figure — **AUC 0.99906** TCGA vs
  ALCHEMIST against a within-TCGA control of 0.50016 — is carried into §5.1 with it.

**No ALCHEMIST number enters any results section, table or figure of P1.** It appears only in §1.3
and §5.1, as scope.

### 0.3 A drop range that no row of its own table produces

The abstract, contribution 1 and §4.2(3) all said adjusted joint site accuracy "falls **21–45×**".
Recomputed from the six rows of §4.2's own site table:

| state | raw | adjusted | ratio |
|---|---:|---:|---:|
| d2_h wsi_biology | 0.3633 | 0.0118 | 30.8× |
| d2_h full_biology | 0.2630 | 0.0101 | 26.0× |
| d2_h rna_biology | 0.2563 | 0.0074 | 34.6× |
| d2_i wsi_biology | 0.2348 | 0.0052 | 45.2× |
| d2_i full_biology | 0.2689 | 0.0085 | 31.6× |
| d2_i rna_biology | 0.2744 | 0.0079 | 34.7× |

The range is **26–45×**. **No row produces 21×.** Corrected in all three places, with the arithmetic
shown in §4.2 so a reader can re-derive it.

**It was not P1's error.** `v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md` §T1.3 says
"21–45×" directly above the same table, and P1 inherited it verbatim. **That file still says 21–45×
and should be corrected there too** — it is a result file owned by the Track 1 work, not touched
here, and is listed as a punch-list item below.

### 0.4 A ledger count that had been correct and quietly stopped being correct

§4.14 and `P1_FIGURES.md` S7 said `GATE_LOG.md` "current contents: **101 rows** — 62 gates, 39
observations, 7 failed gates". Census of the file as it stands (`csv.DictReader`, 2026-08-05):

```
total rows 126   PASS 74   OBSERVED 39   FAIL 13
rows tagged P1_track1_negative_control_battery: 101   PASS 55   OBSERVED 39   FAIL 7
```

**The 101/62/39/7 figures are exactly right for this paper's contribution and exactly wrong for the
file.** `GATE_LOG.md` is append-only and shared: 25 further rows have since been written by
`D3_purity_sensitivity`, `D2_proliferation_deflation`, `E0_proliferation_stratification` and
`D2_3_per_axis_proliferation`, six of them FAIL, none of them P1's. Both passages now name the
experiment tag and say explicitly that the file's line count is not the number to quote. This is the
same class of defect as P2's item 6 — a count that was true when written, in a file that keeps
growing underneath it.

### 0.5 The figure plan asserted a bound the draft had already qualified twice over

`P1_FIGURES.md` F2 panel (e) told the artist to annotate that "the labels alone account for
**6.0–11.2%** of the channel's excess over its own null", full stop. By 2026-08-05 the draft body had
qualified that number **twice** and the figure plan had neither qualification:

1. 6.0–11.2% holds **only** when the labels are encoded in the same design columns the adjustment
   residualises against. Re-encoded on the operator's own frozen design the same comparison reads
   −0.3% (*p* = 0.605), 6.4% (*p* = 0.086), −19.3% (*p* = 1.000) — none distinguishable from zero
   (`inductive_channel_and_ceiling_result_20260804T2345Z.md`).
2. It is a **single-partition** reading. At identical encoding and identical *n*, twelve partitions
   of the same cohort move the additive/transductive share from +0.0223 to +0.3777 — a **16.9×**
   spread (`inductive_channel_split_stability_20260805T0110Z.md`, "Bad news first" item 4).

A caption drawn from the un-qualified row would have printed a fixed bound for a quantity the project
has established moves by an order of magnitude across partitions. Fixed; both qualifications are now
required annotations, and the row says a bound may be drawn only as the twelve-partition range with
the number of partitions stated.

---

## 1. Citation verification — every marker, individually

Eleven markers were found by `grep -n "\[UNVERIFIED\]\|\[CITATION NEEDED\]"` (the brief predicted
"roughly 11"; the count was exact, but the marker at §1.1 covers *two* references and the §4.9 marker
covers *four*, so sixteen distinct works were checked). Every field below was read off a live
Crossref record, an arXiv abstract page or a publisher chapter-level record in this pass. **Nothing
below is from recollection.**

| # | as written in P1 | verdict | record consulted |
|---|---|---|---|
| 1 | Kather et al., *Nature Cancer* 2020 | **VERIFIED** | Crossref 10.1038/s43018-020-0087-6 — "Pan-cancer image-based detection of clinically actionable genetic alterations", *Nature Cancer* **1(8):789–799**, 2020-07-27, Kather JN, Heij LR, Grabsch HI, Loeffler C, Echle A, … An Author Correction also exists (10.1038/s43018-020-00149-6) and is noted in §2.7. |
| 2 | CHIEF, Wang et al., *Nature* 2024 | **VERIFIED** | Crossref 10.1038/s41586-024-07894-z — "A pathology foundation model for cancer diagnosis and prognosis prediction", *Nature* **634(8035):970–978**, Wang X, Zhao J, Marostica E, … |
| 3 | Prov-GigaPath, Xu et al., *Nature* 2024 | **VERIFIED** | Crossref 10.1038/s41586-024-07441-w — "A whole-slide foundation model for digital pathology from real-world data", *Nature* **630:181–188**, Xu H, Usuyama N, Bagga J, … |
| 4 | Kömen et al., arXiv:2411.05489, 2024 | **VERIFIED**, content checked | arXiv abstract page — "Do Histopathological Foundation Models Eliminate Batch Effects? A Comparative Study", Kömen J, Marienwald H, Dippel J, Hense J, submitted 2024-11-08. The abstract's "the impact of batch effects, e.g., systematic technical data differences across hospitals" supports exactly the use §1.2 makes of it. |
| 5 | Carloni et al., arXiv:2507.22092, 2025 | **VERIFIED**, and it gained a venue | arXiv abstract page — "Pathology Foundation Models are Scanner Sensitive: Benchmark and Mitigation with Contrastive ScanGen Loss", Carloni G, Brattoli B, Keum S, Park J, Lee T, Ahn CH, Pereira S, 2025-07-29. **Accepted as an oral at the MedAGI 2025 workshop, MICCAI** — a peer-reviewed venue the prose did not state, which strengthens the prior art against us. Content is scanner-induced variation, as §1.2 claims. |
| 6 | Schmitt et al., *JMIR* 2021 | **VERIFIED** | Crossref 10.2196/23436 — "Hidden Variables in Deep Learning Digital Pathology and Their Potential to Cause Batch Effects: Prediction Model Study", *J. Med. Internet Res.* **23(2):e23436**, 2021-02-02, Schmitt M, Maron RC, Hekler A, … |
| 7 | Venet, **Dhanasekaran & Sotiriou**, *PLoS Comput. Biol.* 2011 | **WRONG — CORRECTED** | Crossref 10.1371/journal.pcbi.1002240 — "Most Random Gene Expression Signatures Are Significantly Associated with Breast Cancer Outcome", *PLoS Comput. Biol.* **7(10):e1002240**, 2011-10-20, authors **Venet David; Dumont Jacques E.; Detours Vincent**. |
| 8 | Leek & Storey 2007 | **VERIFIED** | Crossref 10.1371/journal.pgen.0030161 — "Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis", *PLoS Genetics* **3(9):e161**, Leek JT, Storey JD. |
| 9 | Johnson, Li & **Rabinowitz** 2007 | **WRONG — CORRECTED** | Crossref 10.1093/biostatistics/kxj037 — "Adjusting batch effects in microarray expression data using empirical Bayes methods", *Biostatistics* **8(1):118–127**, authors Johnson W. Evan; Li Cheng; **Rabinovic Ariel**. (Crossref's `issued` is the 2006 advance-access date; the issue is January 2007 and the conventional citation year 2007 is kept.) |
| 10 | Muirhead (1982), *Aspects of Multivariate Statistical Theory* | **BOOK AND CHAPTER VERIFIED; PAGE NOT** | Crossref monograph record 10.1002/9780470316559 — Wiley 1982, ISBN 978-0-471-09442-5 / 978-0-470-31655-9, Muirhead Robb J. Chapter-level DOIs list **ch. 5, "Correlation Coefficients", pp. 144–195** as the relevant chapter. The section and page carrying the N − R statement have not been read. |
| 11 | Anderson (2003), *An Introduction to Multivariate Statistical Analysis* | **EDITION VERIFIED; CHAPTER AND PAGE NOT** | OpenLibrary — T. W. Anderson, **3rd edition, Wiley, 2003**, ISBN 0-471-36091-0 / 978-0-471-36091-9. Wiley's own page returned HTTP 402 and no chapter-level record was reachable, so the chapter is **not** identified for this book. This is the **one residual `[UNVERIFIED]` in P1**, and it is scoped in the draft to exactly that: the section and page, in books that need to be read rather than queried. |
| 12 | Naik et al. *Nat Commun* 2020;11:5727 (§3.9/§4.9, was `[CITATION NEEDED]`) | **VERIFIED exactly as written** | Crossref 10.1038/s41467-020-19334-3 — vol **11**, article **5727**, 2020-11-16, Naik N, Madani A, Esteva A, Keskar NS, Press MF, Ruderman D, Agus DB, Socher R. |
| 13 | Rawat et al. *Sci Rep* 2020;10:7275 | **VERIFIED exactly as written** | Crossref 10.1038/s41598-020-64156-4 — vol **10**, article **7275**, 2020-04-29. |
| 14 | Shamai et al. *JAMA Netw Open* 2019;2:e197700 | **VERIFIED exactly as written** | Crossref 10.1001/jamanetworkopen.2019.7700 — vol **2**, page **e197700**, 2019. |
| 15 | Couture et al. *npj Breast Cancer* 2018;4:30 | **VERIFIED exactly as written** | Crossref 10.1038/s41523-018-0079-1 — vol **4**, article **30**, 2018-09-03. |

**Eight further DOI-carrying references were re-verified opportunistically** (they were not marked,
but they are load-bearing and each cost one API call). All eight match what P1 says about them:
Howard 10.1038/s41467-021-24698-1 (*Nat Commun* **12**, art. 4423 — the draft's "12:4423" is right);
Murchan 10.1016/j.jpi.2024.100396 (*J. Pathol. Inform.* **15**:100396);
Winkler 10.1016/j.neuroimage.2020.117065 (*NeuroImage* **220**:117065, "Permutation inference for
canonical correlation analysis", Winkler AM, Renaud O, Smith SM, …);
Jiang 10.1101/gr.121095.111 (*Genome Res.* **21**:1543–1551);
Munro 10.1038/ncomms6125 (*Nat Commun* **5**, art. 5125);
Gerard 10.1186/s12859-020-3450-9 (*BMC Bioinformatics* **21**, art. 206 — note the paper's own title
is "Data-based RNA-seq simulations by binomial thinning"; `seqgendiff` is its **package** name, which
the draft uses correctly but which a copy-editor could mistake for a title);
Biwer 10.1103/PhysRevD.95.062002; Marek 10.1038/s41586-022-04492-9 (*Nature* **603**:654–660).

**What could not be done and why.** A batched arXiv Atom query over the twelve arXiv IDs §2.7 records
as "spot-check verified" (2501.16652, 2406.16192, 2501.18055, 2509.15482, 2510.23807, 2309.07778,
2409.09173, 2411.19666, 2307.00369, 1804.06788, 2105.04906, 2210.02885) returned **`Rate exceeded`**
on every attempt across this session — the arXiv API was saturated, plausibly by this project's own
concurrent agents. Those twelve were not marked `[UNVERIFIED]` and are out of this pass's declared
scope, but **they have not been re-checked here and I am not claiming they have**. The two arXiv
items that *were* marked (Kömen, Carloni) were verified through `arxiv.org/abs/` pages instead, which
were reachable. Re-running that batched query is a punch-list item, not a result.

**No fabricated identifier was found.** Every DOI and every arXiv ID in P1 resolves to a real record
whose title matches. The two defects were both in author lists.

---

## 2. `paper/P1_CALIBRA_DRAFT.md` — the full read-through, and what it found

Read end to end (1,711 lines before this pass). Beyond §0's findings:

| # | where | said | is |
|---|---|---|---|
| 1 | abstract + contribution 1 | the labels-only ceiling quoted as "6.0–11.2%" and "−0.3% to 6.4%" with no partition caveat, while §4.2 (edited earlier today) already says every ceiling figure there is "a single partition's reading, not a stable property of the cohort" | both passages now carry the 16.9× twelve-partition spread and say the ceiling is quotable as a direction, not a percentage |
| 2 | abstract + contribution 1 | **silent on the twelve-partition inductive validation entirely** — the strongest new evidence in §4.2 and the paper's only out-of-sample test appeared nowhere in the paper's own summary of itself | both now carry retention 0.987–1.052 (`d2_h`) / 0.922–1.079 (`d2_i`) over twelve partitions with every arm at the permutation floor |
| 3 | §4.2, inductive paragraph | reports the ranges correctly but **never says a predeclared bar fired**, in a paper whose §4.6.6, §4.10 and Appendix A all report failed predeclarations prominently | now states that `d2_i`'s spread is **0.1573** against the predeclared ≤ 0.10 bar (already 0.1573 on the predeclared eight, so not an artefact of the four-partition extension), that "retention 0.9710" may therefore not be written as a number, and why (a small-denominator effect: `d2_i`'s excess is 0.262–0.345 against `d2_h`'s 0.364–0.426) |
| 4 | §3.10 | "Test suite: **275 passed** in 44 s … verified at the time of writing" | 275 was the count when §3.10 was first written; the tree now collects 696. Replaced with the 2026-08-05 run and its one unrelated failure, named. |
| 5 | §2.7 status list | three bullets describing which references had and had not been checked, all now stale | replaced with a per-citation verification table (the one in §1 above), plus the "check author lists, not identifiers" lesson |
| 6 | §5 limitation 15 | "Reference verification is incomplete … every reference must be verified before submission" — generic, and it survived the discovery of two wrong citations without changing | now leads with the two specific errors and what they teach, and names the single scoped residual |

Plus §0.2 (external cohort, four passages), §0.3 (26–45×, three passages) and §0.4 (ledger census,
one passage).

**Not one measured value in this draft was wrong.** Every table, every floor, every attenuation
slope, every retention range and every permutation *p* checked out against its cited source. The
defects were: two author lists, one arithmetic range inherited from a result file, one count of a
file that grew, one stale test-suite line, and three places where a state had closed and the prose
had not been updated.

**Machine-generated tables: there are none in P1.** I checked for the analogue of P2's
`p2_floor_audit.py` generator/test pair. P1 has `v2/tests/test_paper_paths_resolve.py` (cited paths
must exist) and `v2/tests/test_paper_artifact_digests.py` (a quoted number must be hash-identified),
which are *checkers*, not generators — no table in P1 is emitted by a script, so nothing here needed
regenerating rather than editing. Both tests pass on the edited drafts (**19 passed, 1 skipped**),
including on the four new file paths this pass introduced into `P1_FIGURES.md`.

---

## 3. `paper/P1_FIGURES.md` — the same sweep

| # | row | fixed |
|---|---|---|
| 1 | **F2 (e)** | the un-qualified 6.0–11.2% bound — §0.5 above. Both encodings and the twelve-partition spread are now required annotations, and the status is qualified: "`PLOTTABLE` **provided both encodings and the twelve-partition spread are drawn** — the shared-column figure alone is not a plottable bound." |
| 2 | **F2 Data block** | still said "Panel (d) from `v2/research/rebase/nature/PHASE1_RESULT.md`" — **the source of the withdrawn 0.463 → 0.035 pair**, four rows below a header block that had already been corrected to say panel (d) comes from `P1_CANCER_TYPE_CERTIFICATE.json`. The plan contradicted itself about which file the corrected panel is drawn from, and the stale half named the withdrawn source. Rewritten, with an explicit "**Not from `PHASE1_RESULT.md`**". *(P2's audit found the identical shape of defect in F9 — a caption corrected in one place and not in the other. It recurs because a plan file has more than one place to say the same thing.)* |
| 3 | **F2 claim + panels** | no panel anywhere in the plan carried the twelve-partition inductive validation, added to §4.2 today. New **panel (f)** specified: retention strip plot, 24 points, seed-42 partition marked as rank 2/12 on both artifacts, `none`-arm overlay for scale (medians 0.2602 / 0.3521), and a caption required to state the breached ≤ 0.10 spread bar. Status `NEEDS EXTRACTION` — the values exist only in a notebook-entry table. F2's claim line updated to match. |
| 4 | **S7** | the 101-row ledger count, unfiltered — §0.4 above. Now says to filter on the `P1_track1_negative_control_battery` tag and that the file holds 126 rows. |
| 5 | **F9** | the four ER papers' `[CITATION NEEDED]` closed, with the verification date. |
| 6 | **"Figures the paper does NOT have"** | the external-cohort row said "No external cohort has been through the instrument", the same false statement as the draft. Rewritten to the true and narrower form, with an explicit instruction that the ALCHEMIST channel replication **must not be drawn as a floor figure**. |

---

## 4. The dangling claim — `morphology_to_pbs_axis_legibility`

**Finding: no paper draft depends on this claim record, and none references it by name.**
`grep -rn "morphology_to_pbs_axis_legibility"` over the whole tree returns exactly four hits — the
record itself in `v2/research/rebase/nature/claim_evidence.json`, and three notebook entries
(`claim_guards_read_evidence_20260803T1545Z.md` ×2,
`decision_external_cohort_blocker_stays_20260804T2200Z.md`). **Zero hits in `paper/`.**

What the drafts *do* reference is the claim **kind** `legible_axis`, which is a different thing and
stays: P1 §1.3 and §5.1 say `legible_axis`/`gene_attribution` claims are inadmissible and that none is
made; `P3_P4_PLAN.md` §7 and §246/§289/§328 discuss the kind's requirements. None of that depends on
this particular record existing, being open, or being pursued.

**The claim is moot, and P1 already prints the evidence that moots it.** §4.13's sixteen-row
must-beat table — the same table as `t11_t12_must_beat_baselines_20260803T0440Z.md` — reads, in the
draft's own words: *"against ordinary PCA of the same expression matrix, PBS loses 3/4 with a CI
excluding zero and ties the fourth — it never wins."* `PROJECT_GUIDE.md` §3 records P3 as having
pivoted away from PBS-as-competing-basis across five attempts.

**Recommendation (not executed here, per the brief): mark the record withdrawn/superseded rather than
delete it or leave it open.** The record is evidence — it documents a purity check that was really
run and really passed — and deleting it would erase that.

**A safe mechanism, checked against the guards rather than guessed:**
`claim_guards.load_claim_evidence` reads only `record["kind"]` and `record["evidence"]`, and
`evidence_digest` is computed over the `evidence` dict alone. A **claim-level** key — a sibling of
`kind`, `description` and `evidence`, e.g.

```json
"status": "withdrawn",
"superseded_by": "NOTEBOOK_ENTRIES/t11_t12_must_beat_baselines_20260803T0440Z.md"
```

— is therefore ignored by the loader and **not** covered by the digest, so it will not trip
`tests/test_claim_guards.py::test_the_shipped_evidence_file_is_internally_consistent`. Putting it
*inside* `evidence` would break that test, and would also make the guard treat it as an evidence
field. This is stated so whoever executes it does not have to rediscover it.

**Where a corresponding note is needed in the papers: nowhere in P1, and one place in P3.**

* **P1 — no change needed.** P1 never claims PBS-axis legibility, prints the losing table, and its
  `legible_axis` statements are about the claim kind. Adding a note about a record P1 does not cite
  would be noise.
* **`paper/P3_P4_PLAN.md` §1 — one line.** That section already declares the working title and
  hypothesis "(dead)". It does not say that the **registered claim record** carrying that hypothesis
  is still sitting open in `claim_evidence.json`. A reader of the plan would reasonably conclude the
  registry was cleaned up when the hypothesis died. It was not. One sentence there, naming the record
  and its withdrawn status, closes the loop. **Not written by me** — `P3_P4_PLAN.md` is another
  agent's live document today and the withdrawal has not actually been executed yet, so a note saying
  it has would be false.

---

## 5. Tests

Run per the project convention — repository reachable as `morpheus/`, thread-capped, `--basetemp`
because the Windows default temp root is not writable. **No environment was created or installed
into**; the default interpreter (`C:\Python313`, Python 3.13.5) already has what the suite needs, and
`~/venv` on the box was not touched.

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python -m pytest morpheus/v2/tests morpheus/tests -q --basetemp=<scratch>
```

**Final, verbatim: `1 failed, 694 passed, 1 skipped, 447 warnings in 117.63s (0:01:57)`.**

The one failure is
`morpheus/v2/tests/test_attributable_basis.py::test_cross_line_rotation_recovers_a_planted_agreeing_direction`.
**It is not mine and I did not touch it.** `v2/attributable_basis.py` and its test landed in commit
`a3937e7` during this pass and were both still modified in the working tree
(`git status -sb` shows ` M v2/attributable_basis.py`, ` M v2/tests/test_attributable_basis.py`) —
another agent's in-flight work. My changes in this pass are two markdown files under `paper/` plus a
notebook entry and a `PROJECT_GUIDE.md` §3 edit; none of them can reach that module. The two tests
that *do* police the drafts pass:
`test_paper_paths_resolve.py` + `test_paper_artifact_digests.py` → **19 passed, 1 skipped**.

**Two things about the test run that are worth recording because they cost time and will cost the
next agent the same.**

1. **`pytest v2/tests tests -q` from the repository root fails with 71 collection errors**
   (`ModuleNotFoundError: No module named 'morpheus'`). The suite requires the repo to be importable
   *as* `morpheus`, which the directory name `morpheus-rebase` is not. A junction plus `PYTHONPATH`
   is needed. The convention is recorded in P2's audit entry but not in `PROJECT_GUIDE.md` §4; it
   probably should be.
2. **An earlier run of mine reported `517 passed` and it was wrong.** A pre-existing
   `scratchpad/ws/morpheus` directory from an earlier session turned out to be a **stale copy**, not a
   junction (47 test files against the live tree's 58), and pytest happily ran it. It was caught only
   because 517 exactly matched a count recorded on 2026-08-04 and the tree has grown a lot since.
   **A workspace that is a copy rather than a link will silently report yesterday's suite as today's**
   — the same class of hazard as methodology rule 11 and rule 12, in a place neither of them covers.
   A fresh junction (`New-Item -ItemType Junction`) to the live repo was created and used for the
   reported run.

An intermediate run also caught `test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree`
failing on `v2/attributable_basis.py`. That was a race — the module landed before the allowlist entry
for it did, within the same minute — and it had cleared by the reported run. Recorded because it is
the third time this scan has been observed failing transiently on a legitimate file.

---

## 6. Punch list — what is still open on P1

**Nothing on this list is closeable by prose. Each needs a measurement, a retrieval, or an edit to a
file this pass deliberately did not own.**

*Needs retrieval, not compute:*

1. **Muirhead ch. 5 / Anderson ch. 4 at the page level.** The one residual `[UNVERIFIED]` in P1. The
   books are identified (Muirhead's relevant chapter is pinned to pp. 144–195); the section and page
   carrying the N − R statement need the books themselves. Neither is open-access and Wiley's
   chapter pages return HTTP 402.
2. **Re-run the batched arXiv Atom query over the twelve IDs §2.7 lists as "spot-check verified."**
   Rate-limited to failure throughout this session. They are not marked `[UNVERIFIED]` and were not in
   scope, but §2.7 asserts a verification state for them that this pass could not independently
   confirm. Cheap; one call when the API is quiet.
3. **Re-verify the venue/award labels flagged as lower-confidence by the 2026-07-29 quality audits.**
   Carloni turned out to carry an unstated peer-reviewed venue (MedAGI @ MICCAI 2025); the same is
   likely true of others, and it strengthens the prior art against us, so finding them is in the
   paper's interest.

*Needs an edit to a file this pass did not own:*

4. **`TRACK1_NEGATIVE_CONTROLS.md` §T1.3: "21–45×" → "26–45×."** §0.3 above. P1 is corrected; the
   result file it inherited the error from is not. Owned by the Track 1 work.
5. **`claim_evidence.json`: mark `morphology_to_pbs_axis_legibility` withdrawn/superseded.** §4 above,
   including the digest-safe mechanism.
6. **`paper/P3_P4_PLAN.md` §1: one sentence** recording that the registered claim record is (once it
   is) withdrawn. §4 above. Must follow 5, not precede it.

*Needs new measurement (named, not attempted — the GPU was saturated by another agent's run and this
was a QA pass regardless):*

7. **The floors on an external cohort.** ALCHEMIST replicated the **channel**. No injection,
   transmission floor, detection floor or attenuation slope has ever been computed outside TCGA — and
   the floors are what P1 is about. This is now stated as limitation 1 rather than implied by "no
   external cohort", which makes it a concrete, scopeable experiment rather than a blanket caveat.
8. **A per-axis external replication**, if `no_external_cohort` is ever to be discharged for a
   `legible_axis` claim. An aggregate-channel result cannot do it and `_is_discharged` cannot tell.
9. **A plot-ready export of the 24 inductive-retention values** (twelve partitions × two artifacts)
   so F2 (f) can move from `NEEDS EXTRACTION` to `PLOTTABLE`. Not a measurement — the runs exist —
   but it needs the per-run JSONs pulled off the box into the figure data tree.
10. **The `dx_normal` dilution arm.** Named as "the clearly indicated next experiment" by the
    ALCHEMIST entry, because ALCHEMIST's contaminant is adjacent normal tissue from the *same*
    patient and the dilution curve prices only `foreign_tumour`. Needs GPU re-embedding of TCGA
    normal slides — machinery that is now deployed.
11. **Reproduce `sim.py` / `sim3.py` into the repository** (limitation 14, unchanged by this pass).
    The 7.4 × 10⁻¹⁶ standalone verification is still scratchpad-only.

---

## 7. Files changed

* `paper/P1_CALIBRA_DRAFT.md` — citations (§1.1, §1.2, §2.2, §2.4, §2.5, §2.7, §4.9), the external
  cohort scope (abstract ×2, §1.3, §5.1), 26–45× (abstract, contribution 1, §4.2), the ceiling's
  partition caveat (abstract, contribution 1), the inductive validation and its breached spread bar
  (abstract, contribution 1, §4.2), the ledger census (§4.14), the test-suite line (§3.10), and
  limitation 15.
* `paper/P1_FIGURES.md` — F2 claim/artifacts/(e)/(f)/Data/Status, F9, S7, and the
  external-cohort row of "Figures the paper does NOT have".
* `PROJECT_GUIDE.md` §3 — P1's open list, per that file's §5. *(Untracked and authored by another
  agent as of this writing; edited in place, deliberately **not** staged, so its author commits it.)*
* This entry.

Related: [[p2_completeness_pass_stale_open_states_and_the_figure_that_drew_the_wrong_floor_20260804T2030Z]],
[[alchemist_external_replication_RESULT_20260804T2115Z]],
[[decision_external_cohort_blocker_stays_20260804T2200Z]],
[[inductive_channel_split_stability_20260805T0110Z]],
[[inductive_channel_and_ceiling_result_20260804T2345Z]],
[[t11_t12_must_beat_baselines_20260803T0440Z]]
