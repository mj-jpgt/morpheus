# E0 responsive-arm stratification — the blocker's premise is TRUE (responsive perturbations are 3.0× enriched for proliferation targets, p=1.1e-13) and its conclusion still does not follow: removing every proliferation-targeting perturbation costs 2–9% of the gap, indistinguishable from removing the same number at random. `proliferation_deflation` is DISCHARGED for E0.

**Logged:** 2026-08-03 15:30 UTC. Pre-registered in
`NOTEBOOK_ENTRIES/e0_proliferation_stratification_preregistration_20260803T1450Z.md`, committed
`5f331f6`, **before** the run.

**How obtained:** Lambda A100 `150.136.45.194`, `~/ws_d3`. GPU slot taken deliberately with
coordinator authorisation; `nvidia-smi` checked first (D1 training holding ~31 GB of 82 GB, E0 peaks
at ~4.5 GB, so memory uncontended; SM time shared for ~18 min). Everything around it thread-capped.
`morpheus.v2.calibra.e0_proliferation_stratified --perturbation K562_gwps_normalized_bulk_01.h5ad
--tcga tcga_pancan_rna.parquet --tcga-registry tcga_cancer_registry.parquet --annotations
gene_annotations.parquet --transforms signed_log1p,clip_log1p --draws 200 --bootstrap-draws 1000
--seed 42`. Output `~/e0_run/d3/e0_prolif/e0_proliferation_stratified.json`.

---

### Technical

#### Provenance: the reference arms reproduce `E0_RESULT.md` exactly

The module imports E0's own `_arm_result` and `_decision` verbatim. Its `responsive_matched` and
`nonresponsive` arms return **0.0822 / 0.0095** (signed, k=10) and **0.0833 / 0.0095** (clip, k=10) —
the exact values in `E0_RESULT.md` §0. So the stratified arms are measured on E0's instrument, not on
a lookalike. **This is a sensitivity analysis on E0's decision rule; it does not re-run E0's gate
ledger and must not be quoted as a re-certification.**

#### The blocker's premise, tested for the first time — and it is correct

`claim_guards.proliferation_deflation` asserts the responsive arm "is selected on HAVING a detectable
transcriptional effect, which enriches for essential, core-machinery, ribosome and cell-cycle genes."
Nobody had checked. Over all 10,744 non-control K562 perturbations:

| arm | n | targets a proliferation gene | enrichment vs all |
|---|---:|---:|---:|
| all non-control | 10,744 | 6.09% | 1.00× |
| **responsive** (`energy p < 0.01`) | 3,973 | **9.59%** | **1.58×** |
| non-responsive (`energy p > 0.5`) | 1,206 | 3.40% | 0.56× |

Responsive vs non-responsive odds ratio **3.01**, Fisher exact **p = 1.1 × 10⁻¹³**. **The premise is
real and strongly supported.** The blocker was a good hypothesis, not a formality.

#### The result — the premise holds and the conclusion still fails

Every arm n-matched; the control is re-matched to each stratified arm's row count. Gap =
responsive − non-responsive overlap, PC1 stripped, E0's own statistic.

**K562 signed_log1p** — 87 of 956 responsive perturbations (9.1%) dropped:

| k | gap, E0's arm | gap, non-proliferation | gap, **placebo** | retention | **placebo retention** | E0 decision (matched / nonprolif / placebo) |
|---|---:|---:|---:|---:|---:|---|
| 10 | 0.0727 | 0.0704 | 0.0701 | 96.9% | 96.5% | True / **True** / True |
| 25 | 0.0532 | 0.0502 | 0.0498 | 94.3% | 93.5% | True / **True** / True |
| 50 | 0.0405 | 0.0394 | 0.0387 | 97.1% | 95.5% | True / **True** / True |
| 100 | 0.0325 | 0.0312 | 0.0314 | 96.2% | 96.8% | True / **True** / True |

**K562 clip_log1p** — 107 of 956 (11.2%) dropped:

| k | gap, E0's arm | gap, non-proliferation | gap, **placebo** | retention | **placebo retention** | E0 decision |
|---|---:|---:|---:|---:|---:|---|
| 10 | 0.0738 | 0.0726 | 0.0768 | 98.3% | 104.1% | True / **True** / True |
| 25 | 0.0558 | 0.0538 | 0.0548 | 96.4% | 98.3% | True / **True** / True |
| 50 | 0.0442 | 0.0404 | 0.0439 | 91.4% | 99.3% | True / **True** / True |
| 100 | 0.0331 | 0.0310 | 0.0329 | 93.6% | 99.5% | True / **True** / True |

Paired `bootstrap_diff_ci95` — the interval E0's decision actually uses — is entirely above zero for
the non-proliferation arm in **8/8** cells, e.g. signed k=10 **[+0.0464, +0.0664]** at n=[869, 869].

**Against the pre-declared bar:**

| condition | required | observed | |
|---|---|---|---|
| clears E0's decision | 4/4 k, both transforms | **8/8** | ✓ |
| retention | ≥ 70% at every k | **91.4–98.3%** | ✓ |
| within ±15 pp of placebo | yes | **0.4–7.9 pp** | ✓ |

**→ DISCHARGED.** The falsifier — "the alignment was proliferation", requiring retention < 40% or the
placebo retaining materially more — did not fire in any cell. In `clip` the placebo actually retains
*more* than the real cut at 3 of 4 k, i.e. dropping proliferation-targeting perturbations costs no
more than dropping arbitrary ones.

#### A dead end, logged: the first run was defective and its answer was wrong

The first execution reported `responsive_exceeds_nonresponsive: False` at **8/8** cells and looked
like a clean falsification. It was a bug of mine. I stratified the responsive arm *after*
`_energy_arms` had already n-matched it, so every stratified comparison was 869 rows against a
956-row control, and E0's `_decision` correctly refuses to certify an arm pair with
`arms_are_n_matched: False`. The decision was returning False for **bookkeeping, not biology** — the
paired `bootstrap_diff_ci95` was [+0.0457, +0.0661], entirely above zero, in the same output.

What caught it was reading the decision's internals rather than its boolean. Had I reported the
boolean I would have filed a confident false negative and left the blocker standing on an artefact.
The defective output is retained at `~/e0_run/d3/e0_prolif/UNMATCHED_control_DEFECTIVE.json`. The fix
adds a `nonresponsive_strata_matched` control and a regression test that every arm has a distinct
seed offset (a second bug, an eagerly-evaluated `dict.get` default, was caught the same way).

### In plain terms

The objection was that E0's result might be nothing but cell division: the perturbations that "did
something" are disproportionately ones that hit growth machinery, and tumours differ from each other
in growth rate, so the two could line up for the most boring reason in cancer biology.

The first half of that objection is simply true, and nobody had checked it before. Perturbations with
a detectable effect are three times more likely to target a proliferation gene than perturbations
without one.

The second half is not. Throw out every single proliferation-targeting perturbation and the result
keeps 91–98% of its strength and still passes E0's own test in all eight cells. Throw out the same
number of randomly chosen perturbations instead and you lose just as much — sometimes more. So what
was lost was the perturbations, not the proliferation.

### Meaning for the claim

**`claim_guards.proliferation_deflation` is DISCHARGED for E0's `transfer` claim**, on the remedy the
blocker itself names, with a size-matched placebo. This is recorded in the new evidence file rather
than by editing a test fixture (see the companion entry).

**E0 remains INADMISSIBLE.** `single_platform` stands and is untouched: K562 and RPE1 are two lineages
on one Perturb-seq protocol, and nothing here changes that. The guard still bites, and a test asserts
that it does.

Scope, stated precisely:
- Discharged **for K562, both transforms, k ∈ {10,25,50,100}**. RPE1 was `UNDECIDABLE` in E0 for want
  of non-responsive rows and was not run; it can neither confirm nor refute.
- The stratification is by **target gene identity**. It does not remove proliferation *signal* from
  the response profiles of the surviving perturbations — a knockdown of a non-proliferation gene can
  still produce a cell-cycle response. A stronger test would residualise the proliferation programme
  out of the gene space on both sides; that is remedy #1 applied to E0, and it was not run.
- 667 proliferation genes (MSigDB Hallmark union, as `build_gene_annotations.py` defines it). A
  different definition would move the cut; `prol_top100`-style sensitivity was not repeated here.

Together with the two earlier placebo-controlled results — D2.3 (85–95 of 95 non-proliferation PBS
axes legible at ~90% of median) and the deflated D2 gap (proliferation cut indistinguishable from a
rank-matched arbitrary cut) — the proliferation confound has now been tested three independent ways
and survived all three.

### Limitations

1. `responsive_prolif_only` (n = 87 / 107) is below `cfg.min_rows` = 151 and is correctly reported
   `unavailable`. So the complementary statement — "the proliferation perturbations *alone* also
   align" — is **not** measured, and must not be asserted.
2. `--draws 200` for the Haar null rather than E0's 1,000. The Haar floor is explicitly not used by
   E0's decision rule, and bootstrap draws were kept at 1,000 where the decision does depend on them.
3. One seed (42). The stratification and placebo draws are deterministic at that seed.

### Files / commits

- `~/e0_run/d3/e0_prolif/e0_proliferation_stratified.json` (result),
  `~/e0_run/d3/e0_prolif/UNMATCHED_control_DEFECTIVE.json` (the logged dead end),
  `~/e0_run/d3/run_e0_prolif.sh`, `~/e0_run/d3/logs/e0_prolif{,2}.log`
- Code: `v2/calibra/e0_proliferation_stratified.py`, `v2/tests/test_e0_proliferation_stratified.py`
