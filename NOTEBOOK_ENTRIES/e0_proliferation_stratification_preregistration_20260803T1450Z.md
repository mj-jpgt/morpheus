# PRE-REGISTRATION — E0 responsive-arm stratification by proliferation loading (`claim_guards` remedy #2), and a real evidence file for the claim guards

**Logged:** 2026-08-03 14:50 UTC. **Written and committed BEFORE the E0 stratification runs.**
Follows `3e6333f` and `f35de66`.

**How obtained:** planning only. Inputs inspected: `v2/calibra/e0_basis_transfer.py`,
`v2/research/rebase/nature/E0_RESULT.md`, `v2/calibra/claim_guards.py`, `tests/test_claim_guards.py`,
`~/e0_run/data/{K562_gwps_normalized_bulk_01.h5ad,tcga_pancan_rna.parquet,tcga_cancer_registry.parquet,gene_annotations.parquet}`.

---

## 1. What is being run, and what it is not

`claim_guards.proliferation_deflation` names two remedies. Remedy #1 (regress the programme out) was
executed on D2 and returned a **placebo-controlled null for proliferation**, not a certificate
(`d2_proliferation_deflation_result_20260803T1440Z.md`). This runs **remedy #2 — stratify the
responsive arm by proliferation loading — on E0 itself**, which is the claim the blocker is attached to.

E0's statistic: principal-angle overlap between a Replogle perturbation-response subspace and TCGA
bulk RNA, PC1 stripped (`primary_offset=1`), at k ∈ {10, 25, 50, 100}, q=150. E0's decision rule is
**not** beating the Haar floor; it is the n-matched **responsive arm's bootstrap interval lying
entirely above the non-responsive arm's**. `E0_RESULT.md` §0: K562 signed, gap **+0.0727 / +0.0532 /
+0.0405 / +0.0325** at k = 10/25/50/100, ~10–11% of the split-half ceiling.

**Implementation is a separate module, `v2/calibra/e0_proliferation_stratified.py`, which imports and
calls E0's own `_arm_result` and `_decision` verbatim.** The gated E0 runner is not edited. This is a
sensitivity analysis on E0's decision rule, not a re-certification of E0, and it does not re-run E0's
gate ledger — stated plainly so nobody reads it as one.

## 2. The arms

All arms are **n-matched to a common row count** (E0 already does this; subspace-estimation quality is
a function of n, so an unmatched comparison measures sample size as much as biology).

| arm | definition |
|---|---|
| `nonresponsive` | `energy_test_p_value > 0.5` — E0's biological negative control, unchanged |
| `responsive_matched` | `energy_test_p_value < 0.01`, n-matched — E0's own responsive arm, the reference |
| **`responsive_nonprolif`** | responsive perturbations whose **target gene is NOT** in the MSigDB Hallmark proliferation union (E2F_TARGETS, G2M_CHECKPOINT, MYC_TARGETS_V1/V2, MITOTIC_SPINDLE) — the same set `build_gene_annotations.py` uses |
| **`responsive_placebo`** | responsive perturbations with an equal **number** of randomly chosen perturbations dropped, seed 20260803 — the size-matched control |

Target genes come from the validated `gene_transcript` index (`_parse_targets`). If that index is
unavailable the stratified arm is emitted as `unavailable_*` per repo convention and **no verdict is
issued** — it is not silently replaced by a proxy.

## 3. Why the placebo is mandatory here, and how to read it

Dropping proliferation-targeting perturbations makes the responsive arm smaller and re-draws the
n-match. A shrinking gap could be that, not proliferation. The size-matched placebo removes the same
number of perturbations at random. Fixed reading, identical in logic to the D2 run:

- **`responsive_nonprolif` gap ≈ `responsive_placebo` gap** ⇒ removing proliferation targets does
  nothing beyond removing that many perturbations ⇒ **the alignment is not proliferation-specific**.
- **`responsive_nonprolif` gap ≪ `responsive_placebo` gap** ⇒ the alignment *was* proliferation ⇒
  **the falsifier fires and the blocker stands.**

Note the asymmetry with D2 and take it seriously: in D2, placebo ≈ real meant "nothing special about
proliferation", which was a *null*. Here the same equality is a *positive* discharge, because here the
quantity being defended (E0's gap) is the claim itself rather than a difference between two arms.

## 4. What counts as pass / fail

Primary context: **K562, both transforms** (`signed_log1p`, `clip_log1p`). RPE1 is not run: it was
`UNDECIDABLE` in E0 for want of non-responsive rows, so it can neither confirm nor refute.

Reference quantity: gap = `responsive_arm.pc1_removed_overlap − nonresponsive.pc1_removed_overlap`,
per k, per transform. Retention = gap(`responsive_nonprolif`) ÷ gap(`responsive_matched`), **within
the same run**, never across runs.

| verdict | condition |
|---|---|
| **DISCHARGED** | `responsive_nonprolif` clears E0's own non-overlapping-bootstrap decision against `nonresponsive` at **4/4 k in both transforms**, retention **≥ 70%** at every k, **and** retention is within ±15 percentage points of the placebo's retention |
| **PARTIAL** | retention 40–70%, or the decision fails at 1–2 of the 8 (k × transform) cells |
| **FALSIFIER FIRES — the alignment was proliferation** | retention **< 40%**, or the decision fails at ≥ 3 of the 8 cells, **and** the placebo retains materially more |
| **NO VERDICT** | target index unavailable, or fewer rows than `cfg.min_rows` (=q+1=151) in the stratified arm |

If DISCHARGED is reached, `proliferation_deflation` is discharged for E0's `transfer` claim and the
guard state is updated deliberately (see §5). If anything else is reached, it is not, and the reason
is recorded.

**Cost/contention.** ~15 min GPU. `nvidia-smi` checked before launch; D1 training holds ~30 GB of
82 GB and E0 peaks at 4.5 GB, so memory is not contended, but SM time is — the run is kept to the two
K562 transforms and launched once, not swept. Everything around it stays thread-capped and CPU.

## 5. Second task — making the claim guards actually guard something

**The defect.** `validate_claim` reads evidence from **nowhere**. Nothing in production builds a claim
dict. The project's record of E0's admissibility exists only as a hardcoded fixture at
`tests/test_claim_guards.py:135`. So the guard system cannot catch anything; it can only restate a
decision someone already made, and "discharging a blocker" means editing a test to say so. That is the
same failure mode this project keeps finding elsewhere: **a check that looks like it ran and didn't.**

**Pre-declared requirements, written before the code:**

1. A checksummed JSON **evidence file** under `v2/research/rebase/nature/`, one record per claim,
   each carrying: claim `kind`, the evidence fields `_DISCHARGED_BY` names, and for every field a
   **provenance triple** (the run that produced it, the notebook entry, the commit).
2. `validate_claim` gains a path-taking loader that reads that file. **Unreadable, missing, or
   checksum-mismatched evidence must be inadmissible**, never permissive — the same default that
   already governs an unknown claim `kind`.
3. **Evidence with no provenance is not evidence.** A field set to `True` without a resolvable run,
   entry and commit must be treated as absent, so that "someone typed True" cannot discharge anything.
4. The fixture test is rewritten to assert the guard's **logic** against synthetic claims — each
   blocker fires when its field is absent/False, clears when properly evidenced, and every failure
   mode above is exercised — instead of pinning one hardcoded verdict.
5. The real E0 record moves out of the test and into the evidence file, with its **current** state
   (`proliferation_controlled` set by §4's outcome, `single_platform` still undischarged), so the
   project's claim state is data with provenance rather than a line of test code.
6. A test asserting that **the evidence file's E0 record is inadmissible while `single_platform`
   stands** — i.e. the guard keeps biting after this change, and this refactor is not a quiet
   loosening.

**Falsifier for the refactor itself:** if after the change a claim can be made admissible by editing
only the evidence file, with no analysis run and no provenance, the refactor has failed and must be
reverted. A test asserts exactly that.

### Files / commits

- To be written: `v2/calibra/e0_proliferation_stratified.py`, output `~/e0_run/d3/e0_prolif/`.
- This pre-registration, committed before the E0 run.
