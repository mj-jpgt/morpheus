## 2026-08-04 08:00 UTC — D2's headline is coordinate-system dependent: changing the exam's basis moves the arm contrast by +0.12, which is the whole effect. The PCA comparison is NOT circular.

**Logged:** 2026-08-04 08:00 UTC. **Predeclared in**
`NOTEBOOK_ENTRIES/PREDECLARED_D2_coordinate_system_confound_20260804T0640Z.md`, committed
(`f51f36e`) before any number below was computed.

**How obtained.** Workspace `~/ws_d2sym/morpheus` on the A100 (`150.136.45.194`), built from
`git -c core.autocrlf=false archive HEAD` at commit `80a14d6` and verified **452/452 files by git
blob SHA-1** (0 mismatched, 0 missing, 0 extra). Every module used is byte-identical at `064baea`
(the HEAD a concurrent agent moved to mid-build):
`git diff --name-only 80a14d6 064baea -- v2/research/rebase/d2_compare.py v2/calibra/
v2/paired_bootstrap.py v2/survival_evaluation.py v2/tests/` is empty. CPU only,
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`; the GPU chain
was not contended for. **No statistic is computed inline anywhere in this work** — every channel
number comes from `calibra.spectral.top_canonical_correlation` / `heldout_top_cca`, every
residualisation from `calibra.residualise`, every interval from `paired_bootstrap` or
`survival_evaluation`.

### 0. The anchor: the workspace reproduces published D2 exactly

Before anything new, `d2_compare` was re-run on this workspace against `frozen_rna_targets.npz`
with the published invocation. Point estimates, **residualised** block (cancer + pooled TSS, 108
design columns, 84 sites, `n_test = 2,766`), top-CCA at 16 components:

| exam | seed 42 | seed 43 | seed 44 | published |
|---|---:|---:|---:|---|
| untrained 40, Δ (PBS−H) | **−0.1325** | **−0.1089** | **−0.1226** | −0.1325 / −0.1089 / −0.1226 ✓ |
| all 90 non-control, Δ | −0.1359 | −0.1077 | −0.1192 | −0.1359 / −0.1077 / −0.1192 ✓ |
| `random_control`, Δ | −0.0099 | −0.0280 | −0.0268 | −0.0099 / −0.0280 / −0.0268 ✓ |

Every point estimate reproduces to four decimals. The full 2,000-repeat bootstrap on the untrained
40 (`ANCHOR_geneset40.json`):

| seed | Δ | patient CI₉₅ | p | cancer CI₉₅ | p | published patient / cancer CI |
|---|---:|:---:|---:|:---:|---:|---|
| 42 | −0.1325 | [−0.1605, −0.0993] | 0.0000 | [−0.1792, −0.0632] | 0.0010 | identical to 4 dp |
| 43 | −0.1089 | [−0.1459, −0.0733] | 0.0000 | [−0.1604, −0.0108] | 0.0125 | [−0.1460,−0.0749] / [−0.1623,−0.0118] |
| 44 | −0.1226 | [−0.1483, −0.0867] | 0.0000 | [−0.1643, −0.0427] | 0.0010 | [−0.1502,−0.0866] / [−0.1653,−0.0411] |

Seed 42 reproduces exactly; seeds 43 and 44 agree to ~0.002 on the interval bounds and are identical
on the point estimate. That residue is expected and benign: `d2_compare` seeds each pair's bootstrap
at `args.seed + pair_index`, and the published readout was run as separate per-pair jobs
("13 jobs, ~13 minutes", `D2_RESULT.md` §5), so pairs 2 and 3 drew a different resample sequence
there than in this single three-pair invocation. Both CIs exclude zero in 3/3 either way. The
instrument and the workspace are sound, so what follows is about the claim, not the plumbing.

---

### 1. THE SYMMETRIC EXAM — both arms scored on the 128 PBS dictionary codes

Arm I's supervision file is `pbs_targets_k128_v2.npz`, sha256
`4f7d6f409988a8191bd41a84cea0a2e12096ac50f41bfe287e72ef6f7e40fd40`, recorded as `pbs_target_sha256`
in `D2_PAIR_MANIFEST.json` for all three seeds and byte-identical to the file used as the target
block here. **This exam is arm I's own supervision**, exactly as the gene-set exam is 50/90 arm H's
own supervision. Same `d2_compare`, unmodified; only `--targets` changes.

**Residualised block, top-CCA at 16 components, n_test = 2,766:**

| seed | exam | Hallmark (arm H) | PBS (arm I) | Δ (PBS−H) |
|---|---|---:|---:|---:|
| 42 | gene sets, untrained 40 | 0.6126 | 0.4800 | **−0.1325** |
| 42 | **PBS codes (128)** | 0.5560 | 0.5462 | **−0.0098** |
| 43 | gene sets, untrained 40 | 0.5970 | 0.4882 | **−0.1089** |
| 43 | **PBS codes (128)** | 0.5202 | 0.5291 | **+0.0088** |
| 44 | gene sets, untrained 40 | 0.5983 | 0.4757 | **−0.1226** |
| 44 | **PBS codes (128)** | 0.5335 | 0.5309 | **−0.0026** |

**Full 2,000-repeat paired bootstrap on the PBS-code exam** (`T1_pbs_codes.json`), same design,
same partition, same seed:

| seed | Δ (PBS−H) | patient CI₉₅ | p_improve | cancer CI₉₅ | p_improve |
|---|---:|:---:|---:|:---:|---:|
| 42 | −0.0098 | [−0.0422, **+0.0232**] | 0.3080 | [−0.0498, **+0.0317**] | 0.2990 |
| 43 | +0.0088 | [**−0.0177**, +0.0426] | 0.7380 | [**−0.0237**, +0.0503] | 0.7080 |
| 44 | −0.0026 | [−0.0342, **+0.0230**] | 0.3810 | [−0.0454, **+0.0373**] | 0.4115 |

**Reading, against the predeclared rule.** The predeclared third branch is taken, and it is met in
its strict form: **all six intervals cover zero** — patient and cancer, in 3/3 seeds — with
`p_improve` between 0.30 and 0.74. Δ is −0.0098 / +0.0088 / −0.0026, the sign is not consistent
across seeds, and every value sits inside the −0.0099…−0.0280 band the published `random_control`
negative control already occupies. Per the predeclaration this is **uninformative about which
supervision is better in that space**, and it is reported as such and not read either way. Arm I
does **not** win on its own codes, so the clean "each arm is better in its own space" refutation is
*not* what happened.

For contrast, on the gene-set exam the same instrument put every one of the six intervals *away*
from zero at p ≤ 0.0125. Same arms, same patients, same design, same statistic, same number of
repeats — only the coordinate system of the target block differs.

**But the contrast between the two exams is the finding, and it is large and stable:**

| seed | Δ on gene sets | Δ on PBS codes | shift from changing the exam's basis |
|---|---:|---:|---:|
| 42 | −0.1325 | −0.0098 | **+0.1227** |
| 43 | −0.1089 | +0.0088 | **+0.1177** |
| 44 | −0.1226 | −0.0026 | **+0.1200** |

**Rotating the exam from Hallmark coordinates into the dictionary's own coordinates moves the arm
contrast by +0.118 to +0.120 — i.e. by the entire size of the published effect**, in all three
seeds, with a spread of 0.005. Decomposed per arm, each arm is relatively better in its own space
and both effects are real: arm I gains ~+0.05 moving from gene sets to PBS codes (0.480→0.546,
0.488→0.529, 0.476→0.531) while arm H loses ~−0.06 (0.613→0.556, 0.597→0.520, 0.598→0.534).

**What this does to the published claim.** D2_RESULT.md §6 currently reads "PBS underperforms
Hallmark supervision on the held-out molecular channel, by ~0.11–0.13 top-CCA". *The held-out
molecular channel* is not what was measured. What was measured is the channel **onto gene-set–valued
targets**. On the only other molecular target space on disk, the two arms are indistinguishable.
The sentence must be narrowed to its actual scope; as written it generalises over a coordinate
choice that carries the entire effect.

**What survives, stated at its real strength.** Arm I never beats arm H anywhere tested — not on
gene sets, not on random-control gene sets, and not on its own 128 supervision codes. On an exam
maximally generous to arm I, arm I can only draw. That is a genuine and non-vacuous residual
finding, but it is a claim of **no disadvantage for Hallmark supervision**, not a claim of a 0.12
advantage for it. P3's headline hypothesis — that interventional coordinates are a *better*
supervision target — remains unsupported; the *magnitude* of its refutation does not.

#### 1a. The −0.12 is unique to gene-set targets — it is absent on every other block on disk

Two readings of the tie were still open after Test 1: either arm I's training bought it exactly
enough to draw on its own coordinates (specific), or the arms tie on *any* expression-derived
128-column block and the gene-set exam is the outlier (generic). Four such blocks already exist on
disk and neither arm trained on three of them. Same `d2_compare` machinery, point estimates only
(`EXAM_PANEL.json`), residualised block, top-CCA at 16, Δ = arm I − arm H:

| exam block | seed 42 | seed 43 | seed 44 | mean |
|---|---:|---:|---:|---:|
| gene sets, untrained 40 | −0.1325 | −0.1089 | −0.1226 | **−0.1213** |
| **PBS codes 128 (arm I's own supervision)** | −0.0098 | +0.0088 | −0.0026 | **−0.0012** |
| PCA basis 128 | −0.0201 | +0.0049 | −0.0284 | −0.0145 |
| gene-label-shuffled PBS 128 (s1) | −0.0359 | −0.0057 | −0.0175 | −0.0197 |
| `random_control` gene sets (90) | −0.0099 | −0.0280 | −0.0268 | −0.0216 |
| size/spectrum-matched random dictionary 128 | −0.0597 | −0.0132 | −0.0454 | −0.0394 |

**The −0.12 appears on exactly one block: the gene sets.** On all four 128-column
expression-derived blocks the arm gap is −0.001 to −0.04 — at or inside the −0.0216 mean of the
published `random_control` negative control. The published effect is 3–100× larger than the gap on
any other molecular target space that exists here.

This also settles the specificity question, and modestly in arm I's favour: arm I does best,
relative to arm H, on precisely the block it was trained on (−0.0012), and progressively worse on
blocks further from its supervision (PCA −0.0145, shuffled-PBS −0.0197, random dictionary −0.0394).
Its training did buy it something, on the order of 0.01–0.04, in its own neighbourhood. Arm H is
still never behind by more than noise anywhere.

#### 1b. The two exams really are different spaces, so the shift is not a relabelling

A +0.12 shift from swapping target blocks would be inexplicable if the blocks were near-identical.
They are not. Canonical correlations between the 40 gene-set targets and the 128 PBS codes on the
same 2,766 test patients (`cca_spectrum` at 16 components), and Roy–Vetterli effective rank of each
block (`spectral.effective_rank`, CANONICAL = centred, order 1):

| block | erank(gene sets 40) | erank(PBS 128) | top CCA | mean of 16 | >0.9 | >0.5 |
|---|---:|---:|---:|---:|---:|---:|
| raw | 26.45 | 67.22 | 0.9584 | 0.4792 | 1/16 | 8/16 |
| residualised | 28.37 | 74.39 | 0.8846 | 0.3507 | 0/16 | 4/16 |

On the residualised block — the one every D2 number is computed on — the two target spaces share
three strong directions (0.885 / 0.850 / 0.817) and then diverge fast: the 4th canonical correlation
is already 0.521 and the mean over 16 is 0.35. The exams overlap in part and differ in most, which
is exactly the configuration in which "which basis is the exam written in" is a live confound rather
than a pedantic one.

#### 1c. The raw block, which has never been reported, is weaker still

The task brief flagged that raw-vs-residualised flips arm orderings. It does. **Unresidualised**
(no cancer/TSS removal), same statistic, same partition:

| exam | seed 42 | seed 43 | seed 44 |
|---|---:|---:|---:|
| gene sets, untrained 40 | −0.0453 | **+0.0043** | −0.0224 |
| gene sets, all 90 | −0.0441 | **+0.0092** | −0.0135 |
| `random_control` | −0.0252 | **+0.0237** | **+0.0255** |
| PBS codes (128) | −0.0306 | **+0.0100** | −0.0098 |

These are **point estimates only** — `d2_compare` residualises unconditionally, so no bootstrap
interval was produced for the raw block and none is claimed. On that understanding: the raw gap is
3–5× smaller than the residualised one, and its **sign reverses in seed 43** on both the untrained-40
and the all-90 exam. The 3/3 sign consistency and the ~0.12 magnitude are properties of the
residualised block specifically. Removing cancer and tissue-source-site is defensible — it is the whole point of
the design — but it must be said that it does not merely clean the estimate, it produces most of it.
Any statement of the D2 effect must carry its block.

---

### 2. A space neither arm trained on

Reported as found rather than forced, per the predeclaration.

**The pre-registered clinical endpoint cannot be run on this contrast at all.**
`p1_evidence/inputs/tcga_clinical_covariates.parquet` labels 690 ER-status and 688 PR-status
patients. **Every one is BRCA, and BRCA is a development cancer in the maximal split: 585 train /
105 val / 0 test.** The D2 comparison lives entirely on the test partition. P1's T1.7(b) covariate
already declared itself in-distribution and therefore weak; on *this* comparison its coverage is not
weak but **empty**. It cannot arbitrate D2 and no amount of care makes it able to.

**No drug-response data** on either machine ties to these patients (the GDSC/Tahoe material on the
box is cell-line, not TCGA).

**Overall survival is available and is the one genuinely neutral endpoint.**
`runs/v22_a10_11v21_20260725/discovery_inputs/tcga_cdr_outcomes.parquet` carries `os_time_days` /
`os_event` for **2,765 of the 2,766 test patients** (PFI for 2,159). It is a clinical endpoint, in
neither training coordinate system, on exactly the D2 partition. Instruments:
`survival_evaluation.evaluate_coxnet_endpoint` (Coxnet fit on development rows only, cancer-grouped
inner CV for alpha, risks exported on held-out test rows) and
`survival_evaluation.paired_cindex_bootstrap` (paired Harrell C-index difference, patient and
cancer-cluster). Top-CCA is deliberately not used on a one-column target — it is a maximum over 256
in-sample directions and `spectral.heldout_single_direction_correlation`'s docstring already records
that this manufactures a channel out of noise.

**Declared deviation from the predeclaration:** 400 bootstrap repeats, not 2,000.
`tasks.harrell_cindex` is an O(n²) pure-Python double loop, measured at **2.52 s per call** at
n = 2,765; the paired bootstrap needs `repeats × 2 modes × 2` calls, so 2,000 repeats is ~5.6 h per
cell against ~67 min at 400. 400 is also the repeat count the published T1.1 baseline bootstrap
used. Recorded here rather than quietly substituted.

**The probe is live, not degenerate.** Coxnet fit on the 8 development cancers and scored on the 21
held-out ones reaches Harrell C-index **0.551–0.571**, comfortably above 0.5. A null below would
have been uninterpretable; it is not one.

**Result — paired C-index, arm I minus arm H, 2,765 test patients, 400 repeats:**

| block | seed | C-index H | C-index I | Δ (I−H) | patient CI₉₅ | p_improve | cancer CI₉₅ | p_improve |
|---|---|---:|---:|---:|:---:|---:|:---:|---:|
| **residualised** | 42 | 0.5587 | 0.5452 | −0.0136 | [−0.0345, **+0.0099**] | 0.1425 | [−0.0464, **+0.0211**] | 0.2450 |
| **residualised** | 43 | 0.5663 | 0.5470 | −0.0193 | [−0.0386, **+0.0010**] | 0.0350 | [−0.0426, **+0.0070**] | 0.0700 |
| **residualised** | 44 | 0.5532 | 0.5361 | −0.0171 | [−0.0362, **+0.0046**] | 0.0725 | [−0.0387, **+0.0074**] | 0.1050 |
| raw | 42 | 0.5592 | 0.5593 | +0.0001 | [−0.0176, +0.0200] | 0.5575 | [−0.0399, +0.0484] | 0.5150 |
| raw | 43 | 0.5512 | 0.5712 | **+0.0200** | **[+0.0014, +0.0377]** | 0.9875 | [−0.0372, +0.0934] | 0.7375 |
| raw | 44 | 0.5705 | 0.5588 | −0.0117 | [−0.0272, +0.0087] | 0.1100 | [−0.0541, +0.0287] | 0.2900 |

**Reading, against the predeclared rule. The predeclared "neither separates" branch is taken, and it
is reported as uninformative rather than read either way.** On the residualised block — the D2
design — **no interval excludes zero in any seed**, patient or cancer. Seed 43's patient interval
comes closest and still covers zero at +0.0010.

What can honestly be said beside that: the *direction* on the residualised block is consistent with
D2 in **3/3 seeds** (arm H ahead by 0.0136 / 0.0193 / 0.0171 C-index, `p_improve` 0.14 / 0.035 /
0.073). That is a real hint that the D2 direction is not purely a gene-set artifact. It is a hint
only. The predeclared bar was a CI excluding zero and it was not cleared, and **no cross-seed
combination was predeclared, so none is performed** — combining three `p_improve` values after
seeing them is exactly the manoeuvre this project's predeclaration discipline exists to prevent.

On the raw block the picture is mixed and one cell points the other way: seed 43's patient interval
**excludes zero in arm I's favour** (+0.0200, [+0.0014, +0.0377]). The raw/residualised split flips
the survival ordering on seed 43 just as it flips the top-CCA ordering there.

**What Test 2 establishes, and what it does not.** It does not rescue the general claim and it does
not refute it. Its firm contributions are two: the *direction* survives on a neutral endpoint at 3/3
but below the predeclared bar, and — separately and more solidly — **the one pre-registered clinical
control this project has cannot be run on the partition its headline results live on.** That second
point is a finding about the evidence base, not about either arm, and it holds regardless of what
the survival numbers had said.

---

### 3. Is the PCA comparison circular? No — and the brief's premise needed correcting first

**Correction to the premise, recorded before the test was run.** The brief supposed PCA is a
*predictor* handed the gene-set targets' basis. In T1.1 it is not: the PCA block and the PBS block
are both **targets**, scored against a fixed frozen representation. So "PCA is handed the targets'
basis" is not an available failure mode. The available one is different, and sharper:
`spectral.cca_spectrum` PCA-whitens **both** sides to `n_components` before taking singular values,
so top-CCA at 16 components depends on a 128-column block **only through that block's top-16
principal subspace**. If that drove the win, "capacity-matched at 128" would be false — 112 columns
would never enter — and the result would report where expression variance sits rather than anything
about interventional content.

**That structural claim is verified, not argued.** `top_canonical_correlation(x, y, k=16)` on the
residualised 128-column PBS block, arm H seed 42 (`SUBSPACE_CHECK.json`):

| transform of the 128-column block | top-CCA (k=16) | Δ vs original |
|---|---:|---:|
| original, all 128 columns | 0.555987077425 | — |
| its own rank-16 PCA reconstruction | 0.555987077425 | −4.4e−16 |
| that rank-16 block under a random invertible 16×16 map | 0.555987077425 | +1.1e−16 |
| top-16 PCs only, each rescaled by a random factor | 0.555987077425 | −5.6e−16 |
| **control:** PCs 17–32 substituted for PCs 1–16 | 0.331745972308 | **−0.224** |
| **control:** spectrum reversed | 0.160398071602 | **−0.396** |

Invariant to 1e−16 under all three transforms that preserve the top-16 principal subspace, and moved
by −0.22 / −0.40 by the two that change it. (A rank-32 reconstruction was tried first as the control
and is *also* invariant — it preserves the top-16 subspace, so it is not a control at all. Recorded
because an inert control reads exactly like a passing one, which is this project's own standing
lesson.) **So at k = 16, 112 of the 128 columns never enter the statistic**, and "capacity-matched at
128" describes how the blocks were built, not what was compared. That is a real defect in how the
T1.1 result is worded — but, as the sweep below shows, it is not what decided it.

**Reproduction first.** The T1.1 main table's statistic is `heldout_top_cca` (its header says so);
the block-level table in the same entry quotes `adjusted_top_cca`. Both reproduce exactly on
`runs/d2_final/artifacts/d2_{h,i}_seed42.npz`, 108 design columns, 84 sites, n_test = 2,766:

| cell (k=16) | PBS | PCA | RANDDICT | published |
|---|---:|---:|---:|---|
| d2_h `wsi_biology`, heldout | 0.5032 | 0.5520 | 0.4551 | 0.5032 / 0.5520 / 0.4551 ✓ |
| d2_h `full_biology`, heldout | 0.8417 | 0.8776 | 0.8102 | 0.8417 / 0.8776 / 0.8102 ✓ |
| d2_i `wsi_biology`, heldout | 0.4605 | 0.4905 | 0.4108 | 0.4605 / 0.4905 / 0.4108 ✓ |
| d2_i `full_biology`, heldout | 0.8634 | 0.8714 | 0.8487 | 0.8634 / 0.8714 / 0.8487 ✓ |
| d2_h `wsi_biology`, adjusted | 0.5504 | 0.5901 | 0.5226 | 0.5504 / — / 0.5226 ✓ |

**The budget sweep — PCA minus PBS at k ∈ {8,16,32,64,128}** (positive = PCA beats the dictionary):

| cell | k=8 | k=16 | k=32 | k=64 | k=128 |
|---|---:|---:|---:|---:|---:|
| d2_h `wsi_biology`, heldout | −0.0106 | +0.0488 | +0.0534 | +0.0487 | +0.0151 |
| d2_h `full_biology`, heldout | +0.0186 | +0.0359 | +0.0189 | +0.0162 | +0.0228 |
| d2_i `wsi_biology`, heldout | +0.0280 | +0.0300 | +0.0243 | +0.0271 | +0.0461 |
| d2_i `full_biology`, heldout | −0.0027 | +0.0080 | +0.0188 | +0.0107 | +0.0224 |
| d2_h `wsi_biology`, adjusted | +0.0302 | +0.0396 | +0.0445 | +0.0364 | +0.0285 |
| d2_h `full_biology`, adjusted | +0.0298 | +0.0324 | +0.0182 | +0.0130 | +0.0156 |
| d2_i `wsi_biology`, adjusted | +0.0162 | +0.0184 | +0.0179 | +0.0283 | +0.0306 |
| d2_i `full_biology`, adjusted | +0.0085 | +0.0128 | +0.0126 | +0.0137 | +0.0108 |

PCA is ahead in **18 of 20** held-out cells and **20 of 20** adjusted cells, and it is ahead at
**k = 128**, where the whole 128-column capacity of both blocks enters the statistic. The predeclared
"advantage is present at 16 and gone by 64/128" branch is **not** taken.

**Paired bootstrap on the block difference** (`wsi_biology`, 400 repeats to match the published
`--n-boot 400`; representation held fixed, the two target blocks resampled together — legitimate
because `top_canonical_correlation(a,b) == top_canonical_correlation(b,a)`, asserted on the real
matrices at 1e-10 inside the script):

| cell | Δ (PCA−PBS) | patient CI₉₅ | cancer CI₉₅ | p_improve (patient) |
|---|---:|:---:|:---:|---:|
| d2_h, k=16 | +0.0396 | [+0.0264, +0.0536] | [+0.0126, +0.0609] | 1.0000 |
| d2_h, k=128 | +0.0285 | [+0.0122, +0.0381] | [−0.0119, +0.0384] | 1.0000 |
| d2_i, k=16 | +0.0184 | [+0.0048, +0.0308] | [−0.0058, +0.0410] | 0.9925 |
| d2_i, k=128 | +0.0306 | [+0.0116, +0.0379] | [−0.0129, +0.0395] | 1.0000 |

The patient CI excludes zero in **4/4**, at both budgets. The cancer CI excludes zero in 1/4 and has
p_improve ≥ 0.91 in the rest.

**The variance-concentration explanation fails, and fails backwards.** Roy–Vetterli effective rank
(`spectral.effective_rank`, CANONICAL) of each **residualised** target block — low means the block's
variance sits in few directions, which is what would make a truncated readout easy:

| block | effective rank (residualised) |
|---|---:|
| PBS 128 | **74.39** |
| random dictionary 128 | 96.75 |
| PCA basis 128 | **97.92** |

PCA has the **highest** effective rank of the three — the least concentrated spectrum — and still
wins at every budget. It is not being handed an easy low-rank subspace; if anything the dictionary
is, and the dictionary loses anyway.

*(Method note, recorded rather than silently corrected: the first version of this script computed
"share of block variance in the top 16 PCs" from an inline `np.linalg.svd` — 0.859 / 0.726 / 0.651
for PBS / RANDDICT / PCA, the same ordering. That is a spectral statistic written by hand, and
`v2/tests/test_effective_rank_canonical.py::test_no_second_definition_exists_in_the_tree` failed on
it the moment the script was vendored. It has been replaced by the canonical imported statistic. The
guard did exactly what it exists to do, on this work, and the number it objected to was mine.)*

**Verdict on Task 3: the comparison is not circular and the T1.1 result stands.** Ordinary PCA of
the expression matrix beats the interventional dictionary at every readout budget from 8 to 128,
with the patient CI excluding zero at both ends of the sweep, and the one mechanism that could have
made it an artifact (spectral concentration handed to the 16-component readout) points the other
way. The deflationary conclusion of `t11_t12_must_beat_baselines_20260803T0440Z.md` should be read
as **stronger**, not weaker, than it was written.

**One incidental correction to that entry.** It reported PBS vs the random dictionary as "wins 2/4,
ties 2/4, never loses". Across the full sweep PBS − RANDDICT is positive in **20/20** cells on both
statistics, and the paired bootstrap gives patient CIs excluding zero in 3/4 of the cells tested
(d2_h k=16 +0.0278 [+0.0042,+0.0515]; d2_i k=16 +0.0802 [+0.0374,+0.1013]; d2_i k=128 +0.0234
[+0.0011,+0.0306]; d2_h k=128 +0.0090 [−0.0092,+0.0236]). The narrow claim P3 retains —
interventional coordinates beat random projections of the same matrix — is better supported than
that entry allowed. It remains a claim PCA also satisfies, and by more.

---

### In plain terms

We asked whether the perturbation dictionary lost its headline contest because it is worse, or
because the exam was written in its opponent's units. Mostly the units. Re-marking the same two
models on the dictionary's own 128 questions moves the score gap by about 0.12 — the whole gap — and
on those questions the two models are statistically tied. The gap of 0.12 turns up on exactly one
kind of exam out of the six we could set, and that is the kind the published result used. The
dictionary still never *wins* anywhere, so this is not a reversal; but the published sentence claims
a general advantage across "the molecular channel", and that generality is not there.

We then looked for an exam in neither method's units. Only one exists: how long the patients lived.
The pre-registered alternative — breast-cancer hormone-receptor status — turns out to have no
coverage at all on the patients these results are computed on, which is worth knowing on its own.
On survival, the Hallmark method is ahead in all three repeats, but by an amount too small for the
error bars to separate, so we report it as undecided rather than as a win.

Separately, we checked whether the other bad result for the dictionary — that ordinary PCA of the
same data beats it — was rigged by how the scoring works. It was not. PCA wins at every setting we
tried, including the one that uses all 128 dimensions, and the obvious way the scoring could have
favoured PCA turns out to favour the dictionary instead.

### Meaning for the claim

* **D2_RESULT.md §6 must be narrowed.** "PBS underperforms Hallmark supervision on the held-out
  molecular channel by ~0.11–0.13" is supported only for **gene-set–valued targets on the
  residualised block**. On the dictionary's own codes all six 2,000-repeat intervals cover zero, and
  the gap on every other 128-column block on disk sits at or inside the negative control. Both the
  coordinate system and the block must travel with the number.
* **The surviving form** is: *Hallmark supervision is never worse than PBS supervision on any target
  space tested, and is much better on gene-set targets; PBS supervision buys ~0.01–0.04 even on its
  own coordinates, which is not enough to draw level plus a margin.* That still refutes P3's
  hypothesis — interventional coordinates are not a better supervision target — but it does not
  support a 0.12 general effect, and the published verdict should not be quoted as one.
* **Survival is directional support that does not reach the bar.** Arm H ahead 3/3 by 0.014–0.019
  C-index on the residualised block, every interval covering zero. Quotable only as "the direction
  reproduces on a neutral endpoint; the separation does not."
* **T1.1 is unaffected and strengthened.** The PCA result is not a readout artifact, at any budget.
  Its wording needs one repair: "capacity-matched at 128" is not what the k=16 statistic compares.
* **The evidence base has a real hole**, independent of any arm: the only pre-registered clinical
  control on this project has **zero** coverage on the partition every headline result is computed
  on. Any future claim that leans on T1.7(b) as validation for a test-partition result is leaning on
  nothing.

### Files / provenance

Scripts: `~/ws_d2sym/{preview.py,exam_geometry.py,task3_sweep.py,task3_heldout.py,task3_boot.py,task2_cell.py}`.
Outputs: `~/ws_d2sym/out/{PREVIEW_points.json, EXAM_GEOMETRY.json, ANCHOR_geneset40.json,
T1_pbs_codes.json, T3_budget_sweep.json, T3_heldout_sweep.json, T3_boot_d2*_wsi_biology_k*.json,
T2_os_s*_*.json}`, copied into `v2/research/rebase/nature/d2_coordinate_system/`.

`scikit-survival` 0.25.0 is not in `~/venv` and was installed into a **separate** `~/venv_surv`
rather than into the shared environment, because other agents have jobs running against `~/venv` and
a dependency upgrade underneath them would be exactly the kind of silent breakage this project keeps
finding. `torch` and the rest resolve from `~/venv` via `PYTHONPATH`; numpy 2.2.6 / sklearn 1.7.2 /
pandas 2.3.3 are identical in both, checked before use.

**Suite.** `v2/tests` on the verified workspace: **278 passed, 1 failed in 46.37 s.** The one failure
is `test_paper_paths_resolve.py::test_no_box_output_basename_is_actually_in_the_repository` and is
**pre-existing and not from this work**: it fires because commit `84b1bf9` (another agent's P2 figure
work) vendored `P2_METRICS_D*.json`, `RANK_RECOMPUTE*.json`, `D1_PAIR*.json` and
`D2_PER_ARTIFACT_READOUT.json` into `v2/research/rebase/p2/figures/data/`, and the test requires that
newly-vendored files be removed from `BOX_OUTPUT_BASENAMES` so their citations start being checked.
This workspace is a pristine 452/452-verified archive with zero edits, so the failure cannot originate
here. Left for the P2 agent whose files trigger it rather than edited underneath them.
