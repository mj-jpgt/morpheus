# P2 completeness / QA pass — seven closed items still written as open, and the figure that drew the wrong floor

**2026-08-04 20:30 UTC.** Outcome: **FIX** (no measurement).

**Predeclaration status: none required, and none was made.** Nothing here is a new number. Every value
written into the draft, the figure plan or a figure in this pass was already on disk, already
predeclared where a predeclaration applied, and is re-read from its own recorded source. The one code
change (`fig_f9_decorrelation.py`) switches which *already-measured* floor a band is drawn from; it
computes nothing. If a later reader finds a quantity in this entry that does not resolve to a file
that predates it, that is a defect in this entry.

---

## 0. The awkward finding first

**`v2/research/rebase/p2/figures/fig_f9_decorrelation.py` committed, in artwork, the exact error
§4.1a exists to catch.** It drew §4.1's **3.295×** retraining floor as a shaded band on both panels
and printed, in bold, inside the figure:

> `ONE SEED PER LEVEL. The rank change is ×1.854 under R3 and ×1.940 under the canonical statistic —
> BOTH INSIDE §4.1's ×3.295 same-seed retraining floor, the grey band in (a) and (b).`

3.295× is canonical R1 on the residualised **exported** `wsi_biology` block of a **different arm**
(`programme_only`) at **40 epochs**. F9's three runs are `programme_free` on the **fixed held-out
probe at step 400**. Different block, different arm, different duration — the mismatch the audit's
block-matching rule was written to refuse. The probe block's own step-400 floors have been measured
since 2026-08-04 16:20 (**1.4489× R3 / 1.5702× R1**,
`the_probe_block_has_a_floor_at_last_20260804T1620Z.md`), and against those **both changes clear**;
`floor_audit.json` rows 54–55 have said so for hours while the figure said the opposite.

The draft body had the same defect in prose (§4.9a, below). `paper/P2_FIGURES.md`'s F9 caption and
Appendix C had already been corrected; the figure and §4.9a had not. **A caption fixed in the plan
and not in the script is not a fixed caption** — the script is what ships.

F9 now reads its band **per statistic** from
`v2/research/rebase/p2/figures/data/e0_run/d1_probefloor/out/P2_PROBE_FLOORS.json` at step 400, and
prints that clearing carries nothing at one seed per level. The floor is read per statistic and not
once, because §4.1a's own measurement is that the floor is a property of the statistic (1.000× to
3.295× on one block); drawing R1's floor on the R3 panel would have been the second half of the same
error.

---

## 1. What this pass was, and what it was not

A completeness/QA sweep over `paper/P2_RANK_DRAFT.md` (4,100 lines, read end to end),
`paper/P2_FIGURES.md` (1,309 lines), §2's citations, and the figure pipeline. **No experiment was
run and no quantity was computed.** The floor audit's counts were regenerated from
`p2_floor_audit.summary()` rather than typed.

Everything found was one shape of error: **a state that had closed, still written as open.** Not one
rank value in this draft was wrong.

---

## 2. The seven stale open-states, and the evidence that closed each

| # | where | said | is |
|---|---|---|---|
| 1 | §4.1b | "**Eight** further selections clear, and every one of them is §5's" | **eleven** — `p2_floor_audit.summary()` gives 13 fail / 12 clear / 0 unjudgeable of 25 selections; §4.1a finding 1 and §"Status" item 12 already said eleven |
| 2 | §4.9a | the decorrelation rank change is inside §4.1's 3.295× floor, "and on a third block again — the fixed held-out probe, **which has no floor of its own**; §4.1a **rows 48–49**" | the probe has step-400 floors of 1.4489×/1.5702× and **both rows clear**; the rows are **54–55**, not 48–49 (48–49 are §5.2's staleness measurements) |
| 3 | Appendix C, floor-audit bullet | "**19** of the 62 rows still sit on a block for which no floor has been measured"; "for **four** blocks does not exist" | **16** (generated) and **three** (§4.1a) |
| 4 | Appendix C, §5's momentum sweep | the 2.64× step-600 reading "is **inside** §4.1's 3.295× disqualifying floor", full stop | 3.295× never licensed it; its own step-600 R3 floor is 1.749× and it **clears by 1.51×** (audit row 43, §5.4 row 3) |
| 5 | Appendix C, momentum value | "The choice of `m = 0.999` over `m = 0.99` (1.26×) **is not supported by anything**" | it clears an R3 floor of 1.195× built from its own two arms at its own step — **by 5.6%**, statistic-conditionally (under R3 the same ten runs separate the arms by only 1.138×, inside the floor; under R1 cleanly, 1.453× against 1.155×). What is unchanged is what the *value* rests on |
| 6 | `floor_audit.json`, row `5.4-row2-seedvaried` | `rests_on` ended "**It fails the floor by 0.3%**" while the row printed `clears: yes` against 1.3674× | rewritten to the three-state route (fails by 0.3% against the wrong floor → unjudgeable → clears by 2.4×), and the §4.1a table row **regenerated** from `p2_floor_audit.py --markdown` rather than hand-edited |
| 7 | §7 conclusion | a mid-sentence capital: "…and RankMe is / **The** only selection…" | typographic only |

**Item 6 is the one worth dwelling on.** The audit is machine-checked on every *value*, every *ratio*
and every *verdict*, and the test caught none of this, because the stale sentence lived in the row's
free-text `rests_on` field. A checker that validates numbers and not the prose beside them will print
"clears — yes" and "it fails the floor by 0.3%" in the same table row indefinitely. That is worth
recording as a limit of the mechanism §4.1a recommends, not just as a typo.

## 3. `paper/P2_FIGURES.md` — the same sweep, and it was worse there

The figure plan was rewritten 2026-08-04/05 against the narrowed claim, but eight rows still carried
the pre-2026-08-05 state. Corrected, each citing its closing entry:

- header box and **T8**: 13 / 9 / 3 → **13 / 12 / 0**; "one of the nine" → one of the twelve; "the
  other eight" → eleven; "19 of the 62 unjudgeable" → 16; counting history extended. A note now says
  every count in T8 is generated by `p2_floor_audit.summary()` and must be re-read, not edited.
- **T8**'s "the 13 / 11 / 1 split ... is unchanged" for the four new §5.2a rows → 13 / 12 / 0, and the
  reason those four are unjudgeable restated as *the floor exists but at `lr = 2e-4`*, not *no floor
  exists*.
- **T10**: gained the three later floors it does not contain (step 600 at m = 0.999/m = 0 → R3
  1.749×; step 600 at m = 0.999/m = 0.99 → R3 1.195× / R1 1.155×; capacity 64 at step 150 → R3
  1.705×), and its caption's "which is why one §5.4 row remains unjudgeable" is closed.
- **F9** required-annotation 2 told the artist to draw §4.1's 3.295× band — **directly contradicting
  F9's own caption three paragraphs below, which says it must not be drawn at all.** Rewritten to
  specify this block's own floors, and the data row now names `P2_PROBE_FLOORS.json`.
- **S4**: "the best-agreeing arm does not have the best rank" — the ordering read off a 1.036×
  difference that §5.2 and audit row 48 already withdrew — restated as the **equality** it is;
  the m = 0.999-over-m = 0.99 row and the "no like-for-like floor has been measured" clause updated
  (that clause sat two paragraphs after the same row cited T10, which measured it).
- **S9** caption: "on a block with no measured floor" → the floor exists but at the wrong learning
  rate; and the **rank/cosine dissociation is withdrawn** — across-arm cosine spread 0.223 against
  within-arm 0.250 at n = 3 per arm, so no arm difference may be read off the cosine there either,
  and the centring account offered as the alternative was measured and is degenerate.
- "**Figures the paper does NOT have**": the seed-replication row (now clears by 2.4×), the centred
  RNA-cosine row (**measured, and it settles nothing** — seven of nine readings within 0.04 of zero),
  and the like-for-like-floor row (**measured; it is T10**).
- **T6**: `INCOMPLETE` census row → `DOWNGRADED FROM INCOMPLETE`, plus the §2.6 state below.

## 4. §2 citations — three markers cleared, and no fabricated identifier

Verified by retrieval, not recollection. Every field below comes from a named record.

| item | was | is | source |
|---|---|---|---|
| Barlow Twins, Zbontar/Jing/Misra/LeCun/Deny | `[UNVERIFIED]`, "PMLR URL only" | **ICML 2021, PMLR 139:12310–12320**, arXiv:2103.03230; **no DOI exists** (PMLR mints none) | PMLR page + BibTeX `pmlr-v139-zbontar21a`, DBLP `conf/icml/ZbontarJMLD21`, arXiv Atom |
| LDReg, Huang et al. | `[UNVERIFIED]` | **ICLR 2024 (poster)**, arXiv:2401.10474v2 — a conference paper, not a preprint | DBLP `conf/iclr/HuangCEMH024`, OpenReview `oZyAqjAjJW`, arXiv Atom |
| Ruan/Zhang/Wang/Zhang, "Muon…" | `[NOT INDEPENDENTLY VERIFIED — S2 record only]` | **arXiv:2606.09658v1, 2026-06-08**; title and all four authors match exactly; **no venue, no DOI, v1 only** | arXiv Atom + `abs` page |

And the twelve `[S2 RECORD ONLY]` items plus 2404.10947 were resolved in **one batched arXiv Atom
query: 13 of 13 resolve**, every first-author surname and every spelled-out author list matching how
§2.2–§2.3 attribute them. **No identifier in §2 is fabricated.** The marker stays where the *content*
is abstract-level, which is the distinction the project's three historical fabrications actually
turned on. Five of them carry peer-reviewed venues the prose did not state (Cheng 2607.13432 = ICML
2026 Spotlight; Adilova 2606.21593 = ECML PKDD 2026; Dai 2510.17299 = NeurIPS 2025; Zhang/Deidda
2502.04591 = ICLR 2026; Zhang/Jiang 2404.10947 = CVPR 2026), which strengthens the prior art against
us. §2.6 now carries a state paragraph: **no bare `[UNVERIFIED]` remains**, and every residual flag
names an action.

## 5. The figures were rendered, and `~/venv` was not touched

`v2/tests/test_p2_figures.py` has been failing on `ModuleNotFoundError: matplotlib` in the shared
`~/venv` on the box. **That is a box-environment problem and not a figure problem.** On this Windows
workstation the default interpreter (`C:\Python313`, Python 3.13.5) already has matplotlib 3.10.8,
numpy 2.4.3 and scipy 1.16.0, so **no environment was created, installed into, or modified** — not
`~/venv`, not a new venv, not `pip install --user`. Nothing was installed at all.

```
cd v2/research/rebase/p2/figures
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python make_all.py --strict
→ 10 display item(s) written as PDF, SVG and 300 dpi PNG.   (exit 0)
```

**Nine of the ten re-rendered byte-identically in PNG.** Only F9 changed, which is the check that the
other nine were already current despite their file mtimes predating the last data refresh — the
PDF/SVG differences are embedded timestamps and clip-path ids only, and were reverted to keep the
diff honest. The three F9 files are committed.

The test's synthetic corpus gained a `P2_PROBE_FLOORS.json` stub, with **deliberately different folds
for R1 and R3**, so a script that read one floor and drew it on both panels fails the corpus.

## 6. Tests

Run as the project's convention has it — the repository reachable as `morpheus/`, thread-capped,
`--basetemp` because the Windows default temp root is not writable:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  python -m pytest v2/tests tests -q --basetemp=./pytmp
```

Baseline before this pass (HEAD `6f28814`): **638 passed, 1 skipped**, verbatim.
Final, after this pass and after a concurrent agent's commits landed in the same tree
(`80957e3` and predecessors): **664 passed, 2 skipped**, verbatim. Zero failures.

*One transient failure is recorded because it was observed and is not ours.* Between the two runs,
`test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree` failed on
`v2/tests/test_p2_limit2_stress.py`, added by commit `80957e3` — the AST/SVD scan matched the literal
string `np.linalg.svd` inside that file's own *assertion that no SVD is present*. It was not touched
here and had cleared by the final run. It is noted because a scan that cannot distinguish "computes
an SVD" from "asserts the absence of one" will keep doing this, and the allowlist will keep growing
for the wrong reason.

*One failure in this pass was ours and is fixed.* `test_paper_paths_resolve.py` rejected a bare
`extracted/F1_RETRAINING_REPEAT.json` introduced into `P2_FIGURES.md`; rewritten to the full
repository path. The test is doing exactly its job.

## 7. Punch list — what is still open on P2

**Nothing on this list is closeable by prose or by re-rendering. Each needs a measurement or a
retrieval.**

*Needs new measurement (out of scope for a QA pass; named, not attempted):*

1. **A floor on the three remaining blocks** — the in-run training batch, the 16-patient gate batch,
   the 282-patient live checkpoint. §4.3's headline (6.05× against 1.18×) and §4.9's two historical
   instances sit on these, and 16 of 62 audit rows are unjudgeable for want of them. Cost recorded
   per block in `P2_ENVELOPE_FLOORS.json`'s `absent_blocks`; needs a GPU.
2. **A probe floor at `lr = 1e-3` and `4e-5`** — §5.2a's four `direction` rows are unjudgeable
   because the probe floor is at `2e-4`, and §5.2a's own result is that the learning rate moves rank
   more than anything else. Five same-seed repeats per arm at each rate.
3. **A floor measured on the UNSTABLE arm of the exported block.** Every exported-block floor is
   `programme_only`. The probe measurement showed the collapsed arm carries the floor by ~2×, so
   §4.1's 3.295× is very likely an underestimate by about that factor — and the direction of that
   error is *against* the paper's own headline count.
4. **A labelled linear probe on every artifact** — the reference standard RankMe and LiDAR were
   validated against. This is the single most valuable missing measurement in the paper (§6.2).
5. **A per-block ground truth for the D1 arms.** §4.6a re-scores only D2; the D1 column is held fixed
   in every row and its block-stability is unmeasured, not established.
6. **Anything on a second architecture, cohort or modality pair.** `no_external_cohort` undischarged.
7. **D1-A's 9.81 / 1.71 under R1** — a GPU forward pass from surviving checkpoints.

*Needs retrieval, not compute:*

8. **Luisto et al., arXiv:2604.14815, at full text** — the one item in the 453-work census that could
   scoop the framing (Finnish histopathology + rank-geometry selection).
9. **Zhang/Jiang 2404.10947's sign, confirmed in the body**, before it may be cited as a necessity
   violation. Its abstract reports the *opposite* direction from the usual framing.
10. **Fang et al. (ICLR 2024, arXiv:2403.00642) at full text** — the paper a referee will say this one
    duplicates.
11. **Aldeneh et al. (ICASSP 2025) at full text** — it leads §2.2 and only its abstract has been read.

*Housekeeping:*

12. The stale absolute path in the D1 audit chain (`/lambda/nfs/.../d1_v2/artifacts/` against the real
    `/home/ubuntu/e0_run/...`) is **still unfixed**; §6.2 says so and it hid a bootstrap for a day.
13. `~/venv` on the box still has no matplotlib, so `test_p2_figures.py` still fails **there**. It is
    not fixed here, deliberately: installing into the shared environment is against convention. The
    fix is a box-side isolated venv for figure rendering, or accepting that figures render locally.
