## 2026-08-05 03:30 UTC — The labelled linear probe folded into `paper/P2_RANK_DRAFT.md` and `paper/P2_FIGURES.md`. Nine locations edited, nine located; one machine-generated table flagged rather than hand-edited

**Prose only. No statistic in this entry was measured by this work.** Every number folded in comes from
`NOTEBOOK_ENTRIES/p2_labelled_linear_probe_result_20260805T0150Z.md`, predeclared in
`NOTEBOOK_ENTRIES/PREDECLARED_p2_labelled_linear_probe_20260805T0040Z.md`. Nothing was recomputed,
re-derived or rounded here.

---

## 1. What changed, and the two edits that are corrections rather than additions

**§4.1's repeat-2 paragraph was carrying a false sentence and it is withdrawn in the paragraph's own
words.** *"One quantity moves by a factor, the other barely moves at all"* is not true of repeat 2: the
labelled cancer-type probe falls **27%** on that run (0.7413–0.7535 → 0.5459, a fold of **1.380×**;
LDA 0.6571–0.6709 → 0.4440, **1.511×**) against the channel's 5%. The replacement says two quantities
move by a factor and one does not, and then gives the reconciliation — what the probe loses on repeat 2
is *cancer type*, which is the direction §3.2's channel residualises out by construction, so rank and
the unadjusted probe both partly track the confound on `wsi_biology` and the channel does not because
it is built not to. The paragraph is written as a self-correction ("we called this the cleanest single
observation and the description was wrong"), because the corrected reading locates a mechanism the old
one did not have.

**§4.1b's attribution was too narrow and that is the most consequential edit.** The view-conditional
reproducibility failure was written there as a property of the rank metric. It is a property of the
**view**: the labelled reference standard has floors of 1.380×/1.511× on `wsi_biology` against
1.011×/1.012× on `rna_biology`/`full_biology`, on the same five runs. The edit does **not** rescue rank
— its 3.295× floor stands, 0 of 6 differences clear it, and the subsection now says so explicitly in
the same list that states the correction.

**Both 14 and 9 are stated wherever the agreement split appears** (§4.1b, §4.5(c), §4.6a, §6.2, §6.3,
`P2_FIGURES.md`), never rounded into "the probe mostly agrees with X". Likewise 24-of-27 and its three
TP53 exceptions always travel together.

## 2. The nine target locations, one line each

| location | state | what it now says |
|---|---|---|
| **§2.5** | **edited cleanly** | `[STILL NOT MEASURED]` replaced with the result. Headline stated first and it is unfavourable: the probe has rank's reproducibility failure, not the channel's. Both probes described, must-fail control (17/17) recorded. |
| **§4.1** | **edited cleanly** | Repeat-2 paragraph rewritten as above, with the probe column and its provenance note. |
| **§4.1a** | **edited — but the generated table was NOT touched; see §3** | New subsection *"And the labelled reference standard has a floor on the same five runs, with the same view shape"* immediately after finding 5, carrying the three probe rows × three views as a **separate** table. |
| **§4.1b** | **edited cleanly — the largest edit** | Blockquote gains the attribution correction; new "the first of those three axes is not rank's" passage with a three-item list (the rule gets stronger / rank is still unusable / finding 4 loses its strongest reading); the tail-effect paragraph loses its implication that the redistribution was harmless; *"the view that makes rank usable is the view on which it is most often wrong"* gains **"against our own unsupervised readout"** plus the 6/6-against-probe and 14-vs-9 counts. |
| **§4.5(c)** | **edited cleanly** | Second-standard row added, not substituted, after the twelve-of-twelve paragraph; the "most often wrong" sentence there gains the same qualifier; 0/6-vs-6/6 on the six D2 (pair × view) comparisons, 14-vs-9, and the estimator disagreement (2/6, 3/6, 0/6). The existing circularity caveat is extended to the added row. |
| **§4.6** | **edited cleanly** | "Fourth, and largest" now also names the reference standard as a seventh coordinate system and points at §4.6a. |
| **§4.6a** | **edited cleanly** | "Three things follow" → "Four"; new fourth item, in §4.6a's table logic rather than as a new section, with both counts and an explicit statement that the probe rows are *not* folded into the six-block table because that table varies the exam's basis with the standard held fixed. |
| **§4.7** | **edited cleanly** | New deflation bullet in §4.7.3: TP53 picks `programme_free` 3/3 on D1 against a resolvable 0.0040 floor, reported first and in those words; then ATM 3/3, KMT2D 3/3, `grade_high` 2/3, `mut_ARID1A` 2/3, `stage_late` unresolvable, and the 24/27 panel total. Not allowed to stand as a contradiction of §4.7; not omitted. |
| **§6.2** | **edited cleanly** | *"A labelled linear probe on every artifact — Not run"* struck and replaced with a **CLOSED** row carrying all four findings and naming what is still absent: a **paired bootstrap on the between-arm probe difference**, and any probe on a second cohort. |
| **§6.3** | **edited cleanly** | Zaiem et al. conceded a second time with our own numbers (2/6, 3/6, 0/6 estimator disagreements). The "one readout" exposure is now recorded as **partly answered** (24/27 molecular panel) and **partly worsened** (cancer-type probe contradicts it wherever it resolves; 14-vs-9), with the note that it does not retire the objection because the probe resolves 0/6 on the readout view. |
| **`paper/P2_FIGURES.md`, "figures the paper does NOT have"** | **edited cleanly** | Row struck through and rewritten to **CLOSED**, drawable as a fourth panel of F1, with the instruction that F1's repeat-2 callout must gain the probe column and that both counts travel together. |
| **`paper/P2_FIGURES.md`, pending-dependencies table** | **edited cleanly** | *"Not run. No figure depends on it"* replaced with **RESOLVED**, F1 (new panel) + T1, and T1's seventh-coordinate-system requirement. |

**All nine draft sections and both figure rows contained the claim expected of them.** Nothing had
drifted out from under the result entry's edit list; no target had to be searched for elsewhere and
none was skipped.

## 3. The one thing deliberately not done, and why

**§4.1a's floor table is machine-generated and was not hand-edited.** It is rendered by
`v2/research/rebase/p2/p2_floor_audit.py` and asserted verbatim against the draft by
`v2/tests/test_p2_floor_audit.py::test_the_draft_prints_the_rendered_floor_table` (a literal
`render_floors_markdown(audit) in draft`). The result entry asks for three rows per view to be added to
it; adding them by hand would have broken that test and, worse, would have put values into a generated
table that the generator cannot reproduce. **The probe floors are therefore printed as a separate,
clearly-labelled table immediately below it, and the draft's own provenance note says why.**

**Outstanding, and it is a generator change rather than a prose one:** `p2_floor_audit.py` needs
`P2_LABELLED_PROBE.json` as an input before Probe A (logistic), Probe A (LDA) and Probe B (TP53) can
appear as rows of the rendered table. That is out of scope for a prose task and is recorded here rather
than worked around.

**Not touched:** `claim_guards.py`, `claim_evidence.json`, every `PREDECLARED_*` file, every other
agent's working file, and the abstract. **The abstract is flagged, not edited.** §4.1's repeat-2 row
appears there (draft lines ~253–258) as *"its rank down 3.3× and its channel down 5%"*. That is not
false as written — it does not assert that nothing else moved — but it is now incomplete in the same
way §4.1's paragraph was, and §4.1 now says the row "may not be quoted anywhere without its probe
column". **Whoever next owns the abstract should add the probe column to it.**

## 4. Suite state, verbatim

Run as `python -m pytest v2/tests -q -p no:randomly --basetemp=<scratch>` with the repository reachable
on `PYTHONPATH` under the package name `morpheus` (the tests import `morpheus.v2.*`; without both of
those the run reports 54 collection errors and 80 `tmp_path` `PermissionError`s from a concurrent
agent's `pytest-of-mobar` lock, neither of which is a test failure).

* **Before any edit: `605 passed, 2 skipped, 445 warnings in 136.89s`.**
* **After all edits: `605 passed, 2 skipped, 445 warnings in 138.26s`.**

No test changed state. `test_p2_floor_audit.py`'s two draft-versus-generator assertions
(`test_the_draft_prints_the_rendered_table`, `test_the_draft_prints_the_rendered_floor_table`) pass,
which is the check that the generated tables were left alone.

*(The 641/6 recorded at `da5c6a5` is a different run: that suite predates several test modules added
since, and this run's counts are quoted as measured rather than reconciled to it.)*

## 5. Files

- `paper/P2_RANK_DRAFT.md` — §2.5, §4.1, §4.1a, §4.1b, §4.5(c), §4.6, §4.6a, §4.7.3, §6.2, §6.3
- `paper/P2_FIGURES.md` — the "figures the paper does NOT have" row and the pending-dependencies row
- Source of every number: `NOTEBOOK_ENTRIES/p2_labelled_linear_probe_result_20260805T0150Z.md`
- Its predeclaration: `NOTEBOOK_ENTRIES/PREDECLARED_p2_labelled_linear_probe_20260805T0040Z.md`
