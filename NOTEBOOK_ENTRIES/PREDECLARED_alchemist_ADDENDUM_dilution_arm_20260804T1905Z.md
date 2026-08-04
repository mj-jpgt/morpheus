# ADDENDUM to the ALCHEMIST predeclaration — the dilution curve does not price this contaminant

**UTC** 2026-08-04T19:05Z
**Status** written **before** any ALCHEMIST channel number exists. Extraction is at 349 of 1,106
slides; the channel has not been run on the full cohort. This weakens the predeclaration it
amends; it is filed separately rather than edited into it so the original bar stays legible.

**Amends** `NOTEBOOK_ENTRIES/PREDECLARED_alchemist_external_replication_20260804T1830Z.md` §3.

---

## The discrepancy

The task brief states, and my own §3 repeated:

> We have a calibration for this: the dilution curve measured what contamination costs
> (channel retained 1.000/0.999/0.968/0.905/0.804/0.607/0.333 at 0→80% foreign tissue,
> half-loss at ~68%). Use it to state the expected penalty in advance.

The seven numbers are correct and I verified them against
`v2/research/rebase/nature/DILUTION_LOWER_BOUND.md` §2 line by line. **What is not correct is
the implied applicability.** That source says, in its own §4 and §5:

- the measured arm is **`foreign_tumour`**: the bag is diluted with **same-cancer tumour
  patches from *other patients***;
- the three normal-tissue arms — `pooled`, `matched`, `dx_normal` — were **never run**.
  *"all three need GPU re-embedding of normal slides, which was out of scope here"*;
- and the document explicitly forbids the reading I was about to make:
  *"Until those arms run, the correct phrasing is 'the cost of preparation-matched,
  information-free contamination', not 'a lower bound'. The number should not be quoted with
  the word 'lower bound' attached unless this paragraph travels with it."*
- it further records a mechanism arguing the **opposite sign**: normal tissue is off-manifold
  similarly across all patients, so it adds a near-constant offset to a mean-pooled bag, which
  damages between-patient variance *less* than a patient-specific foreign tumour shift does.
  If that dominates, `foreign_tumour` is an **upper** bound on the cost, not a lower one.

**ALCHEMIST's contamination is the untested case.** Tissue-level sampling of an FFPE NSCLC
resection admits adjacent normal lung, stroma, vessels and immune aggregates **from the same
patient**. That is not other patients' tumour. It is much closer to `dx_normal`, which is
exactly the arm nobody has measured.

## What this changes

§3's numeric band stands as *an* anchor, but it is **not a calibration of this deviation**, and
the phrase "the dilution-predicted penalty" must not be used unqualified for ALCHEMIST. The bar
in §4 is **unchanged** — 0.60 and 0.30 were already chosen with slack for uncalibrated
cross-institution shift, and moving them now, with the extraction still running, would be
fitting the bar to an argument.

What changes is the **interpretation** attached to each verdict:

- **REPLICATES (R >= 0.60).** Reads the same. If anything it is now a slightly stronger result,
  because it clears a band derived from a contaminant that may well be harsher than the real one.
- **ATTENUATED BUT PRESENT (0.30 <= R < 0.60).** The shortfall may **not** be attributed to
  contamination by pointing at the curve. The curve prices a different contaminant, in an
  unknown direction. The shortfall is unattributed.
- **FAILS (R < 0.30).** Likewise: "it was just the missing polygons" is **not** an available
  explanation, because no measurement in this project supports it.

## The experiment that would fix this, named so it is not quietly skipped

Run the `dx_normal` arm. `v2/research/dilution/extract_normal_patches.py` and
`select_normal_slides.py` already exist and are the intended path; the blocker recorded in
`DILUTION_LOWER_BOUND.md` §5 is GPU re-embedding of TCGA solid-tissue-normal slides, and that
machinery is now deployed and working on this box for ALCHEMIST. Until it runs, this project has
**no** measured price for same-patient non-tumour contamination, and every statement about what
tissue-level sampling costs — including any made in the ALCHEMIST result — is an extrapolation
from a different experiment.
