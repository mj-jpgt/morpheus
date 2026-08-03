## 2026-08-03 12:30 UTC — P1 written to a submission draft; six contradictions in the evidence base found and reconciled in text, three claims narrowed

**Logged:** 2026-08-03 12:30 UTC. **How obtained:** synthesis only, no new experiments. Read
`paper/P1_CALIBRA_DRAFT.md`, all of `v2/research/rebase/nature/`, `NOTEBOOK.md` (publication plan,
P1 completion plan, running log), all 33 files in `NOTEBOOK_ENTRIES/`, `HANDOFF_PHASE_D.md` §0,
`v2/calibra/claim_guards.py`, `v2/calibra/hest.py`, and `GATE_LOG.md` (101 rows, parsed for
`induced_baseline=`, `detection_floor=` and `heldout=` fields). Local machine; test suite re-run with
thread caps to confirm a green baseline.

### Technical

`paper/P1_CALIBRA_DRAFT.md` rewritten from front matter into a full submission draft: abstract (long
+ 200-word variant), introduction, related work, methods, results with real numbers and 22 tables,
limitations, conclusion, three appendices (predeclaration grades, provenance index, code index).
Every table carries a `provenance` line naming the evidence file, the ledger rows and the box path.
`paper/P1_FIGURES.md` created: 11 main figures, 4 main tables, 8 supplementary items, each with the
exact data file, the single claim it carries, and a `PLOTTABLE` / `NEEDS EXTRACTION` / `NOT MEASURED`
status; plus a closing table of nine figures the paper deliberately does not have, each pointing at
the draft section that says so in prose.

**Framing changed.** The old draft's headline was the induced correlation presented as a novel
methodological observation ("we have not seen reported"). That is withdrawn per `NOVELTY_SEARCH.md`.
The paper now leads with confound-adjustment verification and attenuation ≈ 1, and the induced
correlation appears as a *magnitude* claim with the Yule/FWL identity conceded in the introduction
and a seven-row prior-art table in Related Work.

**Six contradictions found in the evidence base.** All are reconciled in the draft text rather than
silently smoothed:

1. **Induced correlation, 0.067–0.140 vs 0.0748.** Phase 1b's range is the *spread across draws*
   within one run; Track 2's number is the *median over 40 draws*. Different statistics, not
   conflicting measurements. Draft quotes the median with its spread and says so (§4.6.7).
2. **Design column count: 99 vs 108 vs 109.** Three different designs at three different cohort
   sizes, all with the same pooling rule; the difference is which sites survive `min_site_count = 10`
   at each n. Draft carries a design table in Methods §3.2 and labels every result table with which
   design produced it.
3. **Detection floor: 0.20 / 0.25–0.30 / 0.30 / 0.40.** Not a single number. `GATE_LOG.md` shows it
   varies by *target block* on the same artifact (0.2 for random dictionary, 0.4 for PBS and PCA, 0.3
   for curated pathway on `d2_h::wsi_biology`). Draft states it is a property of
   (representation × target block × design × n) and reports the range (§4.5).
4. **Attenuation: 0.974–1.039 vs 1.07–1.12 vs 0.944–1.228 vs 0.855–1.130.** Four runs on different
   cohorts and grids. Draft tabulates all four separately, refuses to pool, and nominates the Track 1
   range as the one to quote (§4.3).
5. **Permutation null median: 0.147 vs 0.140 vs 0.151–0.158.** Different n, different component
   counts and — for D2's 0.140 — a *different permutation procedure* (row-shuffle of the residualised
   target matrix, not within-cancer permutation). Draft says the median is cohort- and
   capacity-specific and must not be carried across (§4.4).
6. **`DILUTION_LOWER_BOUND.md`'s filename contradicts its own §4**, which withdraws the phrase "lower
   bound". Draft uses "the cost of preparation-matched, information-free contamination" throughout and
   states the withdrawal explicitly (§4.10); `P1_FIGURES.md` binds the phrasing into the F10 caption.

**Two source defects found and corrected in the draft, not in the source files:**

* `TRACK1_NEGATIVE_CONTROLS.md`'s must-beat table (§T1.1/T1.2) has a mis-labelled header — it prints
  "baseline" twice and the true column order is (baseline name, PBS value, baseline value,
  difference). Verified against the `heldout=` field of the corresponding `T1.2_baseline_block::*`
  ledger rows (e.g. `pbs::d2_h::wsi_biology heldout=0.5032`, `pca::d2_h::wsi_biology heldout=0.5520`).
  The draft prints the corrected header and records the defect in a footnote.
* The `+107%` effective-rank instance (49.9 → 103.3) comes from an earlier codebase generation and a
  different benchmark statistic, and the `−17%` instance's two arms were never verified matched on
  epochs/LR/step budget (gate G0.4). Both caveats now travel with the four-instance table (§4.11) and
  with figure F11.

**Three claims narrowed.** (i) Induced correlation: from "novel observation" to "magnitude in
correlation units under correct cross-fitted residualisation of exactly orthogonal signals,
cross-modal, at TCGA scale" — nothing broader. (ii) Spike certification: from "first" to "has not
previously been applied to a confound-adjusted cross-modal biological analysis, to our knowledge",
with ERCC/Munro/Gerard/LIGO cited pre-emptively. (iii) The D2 supervision ablation: framed as an
ablation demonstrating the instrument's utility, explicitly *not* as a biological claim.

**One gap surfaced that no evidence file states plainly.** `observed_above_floor = 0` in every state
is correct, because the floor is in single-random-direction units and the headline channel is a
16-component multivariate maximum. The consequence is that **this repository contains no measurement
of the real channel, in the floor's own units, that exceeds the floor.** The instrument certifies
pipeline sensitivity; the significance of the channel rests entirely on the permutation null. The
draft states this in §4.5 and §5.7 and F4 has a dedicated scale panel for it.

**Baseline.** Full suite **275 passed** in 43.65 s with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`. Note for future agents: on this Windows checkout the tests
import `morpheus.*`, so they only collect from a directory containing a `morpheus` → repo symlink
(one exists at `<scratchpad>/ws/morpheus`), and `--basetemp` must be overridden because
`C:\Users\...\AppData\Local\Temp\pytest-of-mobar` is not writable. Running `pytest` from the repo root
gives 47 collection errors that are **not** test failures.

### In plain terms

The paper is written. It is now honestly a methods paper: how to check whether a
morphology-to-molecular analysis could have seen anything, and what you find when you check. The
strongest results are that our correction removes what it claims to remove, that it costs the signal
essentially nothing, and that there is a hard sensitivity limit that more patients cannot lower.

The thing the old draft treated as its discovery — that removing background variables from two
unrelated measurements makes the leftovers look related — turns out to be a formula from 1907 that
five different fields have already published warnings about. All that is ours is the size of it in a
setting nobody had measured. The draft says that in the introduction rather than waiting for a
referee to say it.

Reading everything side by side turned up six places where our own files disagree with each other —
mostly the same quantity measured on different cohorts and then quoted as if it were one number. None
of them is a wrong measurement; all of them would have been an easy referee hit. Each is now
reconciled in the text with both numbers shown.

### Meaning for the claim

* **P1 is now a complete draft, not a plan.** Every remaining gap is named in a Limitations section
  rather than an "Open gaps" list, because they are properties of the paper rather than tasks
  blocking it. The one that still blocks a *stronger* paper is the external cohort.
* **The phase gate in `NOTEBOOK.md` is only half met.** Item B (negative-control battery executed and
  written up, including the losses) is done. Item A (external / second-dataset demonstration of both
  floors) is not, and the draft submits without it as a declared scope decision rather than claiming
  it.
* **Three predeclared predictions failed and are reported as failures** (P0, P6, D2/D4 in the dilution
  arm). The draft's Appendix A tabulates all eleven with grades so the failure rate is visible rather
  than inferable.
* **Nothing in the draft requires a GPU run.** The remaining pre-submission work is bibliographic
  (resolve every `[UNVERIFIED]` and `[CITATION NEEDED]` against a live API), reproducing the two
  scratchpad verification scripts into the repo, and extracting four `NEEDS EXTRACTION` arrays from
  box run outputs for figures F1a, F2b, F4a and F8a.

### Files / commits

`paper/P1_CALIBRA_DRAFT.md` (rewritten), `paper/P1_FIGURES.md` (new), this entry.
