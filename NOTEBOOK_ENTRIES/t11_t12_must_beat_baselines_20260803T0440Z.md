## 2026-08-03 04:40 UTC — Must-beat table with paired bootstrap CIs: ordinary PCA of the expression matrix BEATS the interventional dictionary in 3 of 4 cells and never loses to it

**Logged:** 2026-08-03 04:40 UTC. **How obtained:** `run_calibra` on four target blocks through the identical instrument (`--partition test --levels 0.0,0.05,0.10,0.20,0.40 --n-draws 16 --n-components 16 --n-permutations 500 --seed 42`), plus `p1_evidence/baseline_paired_bootstrap.py --n-boot 400` for the CI on every difference. Lambda box `~/ws_p1`, 2,766 held-out patients, 108-column cancer + pooled-TSS design.

### Technical

A table of point estimates is not a comparison, so every difference below is a **paired** bootstrap:
PBS and the baseline are scored on the *same* resample of patients, 400 draws, and the CI is on the
difference. All blocks share the identical residualisation, confound design, partition and seed; only
the target block changes.

**Held-out top-CCA, PBS minus baseline (n_components = 16):**

| artifact | state | baseline | PBS | baseline | difference | CI95 of difference | verdict |
|---|---|---|---:|---:|---:|---|---|
| d2_h | **wsi_biology** | random dictionary | 0.5032 | 0.4551 | +0.0481 | [−0.0177, +0.0693] | TIE |
| d2_h | **wsi_biology** | **PCA basis** | 0.5032 | **0.5520** | **−0.0488** | **[−0.0734, −0.0183]** | **BASELINE WINS** |
| d2_h | **wsi_biology** | gene-label-shuffled (s1) | 0.5032 | 0.5146 | −0.0114 | [−0.0601, +0.0445] | TIE |
| d2_h | **wsi_biology** | gene-label-shuffled (s2) | 0.5032 | 0.5187 | −0.0155 | [−0.0561, +0.0318] | TIE |
| d2_h | full_biology | random dictionary | 0.8417 | 0.8102 | +0.0315 | [+0.0241, +0.0653] | PBS wins |
| d2_h | full_biology | **PCA basis** | 0.8417 | **0.8776** | **−0.0359** | **[−0.0483, −0.0236]** | **BASELINE WINS** |
| d2_h | full_biology | gene-label-shuffled (s1) | 0.8417 | 0.8140 | +0.0277 | [+0.0197, +0.0660] | PBS wins |
| d2_h | full_biology | gene-label-shuffled (s2) | 0.8417 | 0.8085 | +0.0332 | [+0.0116, +0.0632] | PBS wins |
| d2_i | **wsi_biology** | random dictionary | 0.4605 | 0.4108 | +0.0497 | [+0.0251, +0.1372] | PBS wins |
| d2_i | **wsi_biology** | **PCA basis** | 0.4605 | 0.4905 | −0.0300 | [−0.0429, +0.0053] | TIE |
| d2_i | **wsi_biology** | gene-label-shuffled (s1) | 0.4605 | 0.4245 | +0.0360 | [+0.0111, +0.0926] | PBS wins |
| d2_i | **wsi_biology** | gene-label-shuffled (s2) | 0.4605 | 0.4317 | +0.0288 | [−0.0030, +0.0887] | TIE |
| d2_i | full_biology | random dictionary | 0.8634 | 0.8487 | +0.0147 | [−0.0004, +0.0283] | TIE |
| d2_i | full_biology | **PCA basis** | 0.8634 | **0.8714** | **−0.0080** | **[−0.0233, −0.0001]** | **BASELINE WINS (marginal)** |
| d2_i | full_biology | gene-label-shuffled (s1) | 0.8634 | 0.8408 | +0.0227 | [−0.0029, +0.0288] | TIE |
| d2_i | full_biology | gene-label-shuffled (s2) | 0.8634 | 0.8374 | +0.0260 | [+0.0057, +0.0383] | PBS wins |

Read down the three opponents:

* **vs. size- and spectrum-matched random dictionary:** PBS wins 2/4 cells, ties 2/4, **never loses.**
* **vs. gene-label-shuffled PBS:** PBS wins 3/8 draws, ties 5/8, **never loses.**
* **vs. ordinary PCA of the same expression matrix, fit on development rows only:**
  **PBS loses 3/4 cells with a CI excluding zero and ties the fourth. It never wins.**

**That is the deflationary result the T1.1 spec was designed to force**, verbatim from the module
docstring: "if the top 128 principal components of the same expression matrix are as morphology-legible
as the interventional dictionary, then the dictionary has contributed nothing that ordinary variance
decomposition does not already supply, and the interventional framing is decoration." They are not as
legible — **they are more legible.**

The leak discipline is not the escape hatch: the PCA basis is fit on development rows only
(`split != "test"`), through the same `fit_development_expression_transform` PBS uses, so no held-out
cancer influenced it. Both blocks are 128 columns, so this is capacity-matched.

**A readout-dependence that corrects my earlier T1.5 entry.** At `n_components = 32` with a 500-draw
patient bootstrap, PBS and the gene-label-shuffled block were indistinguishable on all three shuffle
draws (01:50 UTC entry). At `n_components = 16` here, PBS beats the shuffled block on 3 of 8 draws with
a CI excluding zero. **So "the shuffled dictionary is indistinguishable from the real one" is
readout-dependent and was overstated.** The defensible version is: *the gene-label shuffle costs the
dictionary little and sometimes nothing, and never costs it more than a few hundredths of a CCA* — which
is still enough to void gene-level attribution, but is not the flat tie the first entry implied. The
earlier entry should be read with this correction attached.

Block-level supporting numbers from the four `run_calibra` runs (adjusted top-CCA, and the same value
minus the permutation null median of ~0.147, which is the honest above-chance quantity):

| block | d2_h wsi_biology | above null | d2_i wsi_biology | above null |
|---|---:|---:|---:|---:|
| frozen curated pathway (82 cols) | 0.6052 | 0.4569 | 0.4703 | 0.3231 |
| PBS (128) | 0.5504 | 0.4050 | 0.5217 | 0.3730 |
| random dictionary (128) | 0.5226 | 0.3761 | 0.4415 | 0.2919 |

The induced (level-0) baseline differs by block — 0.106 for the curated-pathway block, 0.071 for PBS,
0.081 for the random dictionary on d2_h `wsi_biology` — which is itself a Track 2 confirmation: blocks
whose coordinates the confound design explains less produce less induced correlation.

**Not built:** the text-prior (GenePT-style) and capacity-matched cell-composition blocks. Both need an
external resource (a gene text-embedding table; a deconvolution signature matrix) that is on neither
machine. The plan named these two as the most likely to slip and they slipped;
`claim_guards.composition_attribution` therefore stays undischarged. The zero-parameter cancer-mean
baseline is degenerate on this split by construction — the maximal split holds out whole cancers, so
every test patient gets the global training mean and the state is constant on the evaluated partition.

### In plain terms

We tested our interventional gene dictionary against the cheap alternatives a reviewer would reach for,
scoring all of them with exactly the same instrument on exactly the same patients, and putting an error
bar on every comparison rather than just reporting two numbers.

Against a randomly generated dictionary of the same size and shape, ours wins or draws — good. Against
scrambling which gene each coordinate refers to, ours wins or draws — good.

Against plain old principal-component analysis of the very same expression data, ours **loses** in three
of four comparisons and draws in the fourth. It never wins. PCA is the most obvious thing anyone would
try, it takes one line of code, and on this evidence it reads the images better than the dictionary we
built.

### Meaning for the claim

* **P3's central claim is not supported and must be withdrawn in its current form.** "Interventional
  coordinates beat curated pathway scores / are legible in a way generic decompositions are not" is
  contradicted by the strongest available comparison, with a CI. What P3 can still say is the much
  narrower "interventional coordinates beat *random* projections of the same expression matrix" — true
  here, but a far weaker claim, and one PCA also satisfies.
* **P4 inherits the block:** an axis whose legibility PCA reproduces or exceeds cannot be certified as
  carrying interventional meaning.
* **P1 is strengthened by all of this.** The instrument produced a clean, CI-backed negative against the
  project's own preferred method, on a protocol fixed in advance, and the add/observe separation
  ensured that negative was recorded as a finding rather than quarantining the run
  (`test_a_losing_baseline_is_an_observation_and_cannot_move_the_verdict`). That is precisely the
  behaviour P1 claims to enforce, demonstrated on a result we did not want.
* **A correction to my own 01:50 entry** is recorded above: the T1.5 "indistinguishable" statement is
  readout-dependent and was overstated.

### Files / commits

`p1_evidence/baseline_paired_bootstrap.py`, `v2/build_random_dictionary_targets.py`,
`v2/build_pca_basis_targets.py`, `v2/build_shuffled_pbs_targets.py`.
Results: `p1_evidence/track1/baseline_paired_bootstrap.{json,csv}` and
`p1_evidence/track1/calibra_{pbs,randdict,pca,pbsshuf1}/` under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/`.
