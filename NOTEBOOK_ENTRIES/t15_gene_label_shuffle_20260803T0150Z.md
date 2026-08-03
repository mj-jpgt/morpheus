## 2026-08-03 01:50 UTC — T1.5 must-FAIL control 3: the subspace survives a gene-label shuffle so completely that the fitted dictionary is indistinguishable from a random projection

**Logged:** 2026-08-03 01:50 UTC. **How obtained:** `python -m morpheus.v2.build_shuffled_pbs_targets` (seeds 1,2,3) then `python -m morpheus.v2.calibra.gene_label_shuffle_control --artifact d2_h_seed42.npz --state wsi_biology --partition test --n-boot 500`, Lambda box `~/ws_p1`, 2,766 held-out patients, 99-column cancer+pooled-TSS design, seed 42.

### Technical

**Build note first, because it changes provenance.** `build_pbs_targets --shuffle-gene-labels` exists
and is correct, but it **cannot be run on this machine any more**: it refits the dictionary from the
Perturb-seq reference through a data config, and the config whose digest the frozen target file
records (`data_config_sha256 = 76927870…`) is no longer on disk. Every surviving config declares a
cohort missing 249 of the split's patients, so `restrict_cohort_to_split` refuses. Refitting through a
different config would have produced a shuffled block bound to a *different cohort* than the block it
is compared against — silently. New module `v2/build_shuffled_pbs_targets.py` therefore **rebinds
rather than refits**: it reads `gene_basis` straight out of `pbs_targets_k128_v2.npz` and permutes its
rows, which is exactly what the `--shuffle-gene-labels` path does after the fit. It refuses to write
unless the *unshuffled* reconstruction reproduces the frozen scores at Pearson r ≥ 0.9999 per column.
Measured: min column r = **0.99999999999999**, median **0.9999999999999998**. The frozen dictionary is
reproduced exactly.

**Part (ii), attribution must collapse — PASSES, but it is a build check, not a finding.**

| shuffle seed | median \|Spearman\| by axis index | max by index | median best-match \|Spearman\| | max best-match |
|---|---:|---:|---:|---:|
| 1 | 0.0069 | 0.0307 | 0.0330 | 0.0466 |
| 2 | 0.0073 | — | 0.0327 | — |
| 3 | 0.0077 | — | 0.0334 | — |

All ≤ the 0.05 bar, CIs tight and far from it (e.g. seed 1 by-index CI95 [0.0049, 0.0088]). The
best-match statistic — for each true axis, the largest \|Spearman\| against *any* of the 128 shuffled
axes, which is the harder test a reviewer would ask for — also stays at 0.033. **But this is true by
construction**: the shuffle permutes `gene_basis` rows after the fit, so a permuted loading vector is
uncorrelated with itself. A non-null here would have meant the shuffle failed to take effect. It is
reported as a build-integrity check and may not be quoted as evidence that our gene attribution is
shuffle-sensitive in any deeper sense.

**Part (i), the subspace must persist — PASSES, and that is the problem.**

| shuffle seed | held-out top-CCA, true dictionary | shuffled | paired difference (true − shuffled) | 95% CI of difference | inside CI95 of true? |
|---|---:|---:|---:|---|---|
| 1 | 0.5411 | **0.5600** | −0.0189 | [−0.0489, +0.0384] | yes |
| 2 | 0.5411 | **0.5360** | +0.0051 | [−0.0564, +0.0418] | yes |
| 3 | 0.5411 | **0.4771** | +0.0640 | [−0.0260, +0.0988] | no (0.4771 vs CI95_true [0.4874, 0.5962]) |

CI95 on the true value: [0.4874, 0.5962], 500 patient bootstrap draws.

By the letter of the T1.5 criterion the control passes on 2 of 3 shuffle draws and fails containment
on the third. By the statistic that actually decides whether the two differ — the paired bootstrap
difference — **all three CIs cover zero**, and on one draw the shuffled dictionary scores *higher*
than the real one. Median across draws: shuffled 0.5360 vs true 0.5411.

**What that means mechanically.** After the row permutation the target block is a spectrum-matched
random projection of the same expression matrix. So "the subspace persists" is not the reassuring
result it reads as. The honest statement is: **any spectrum-matched 128-dimensional projection of this
expression matrix is as legible to `wsi_biology` as the fitted interventional dictionary is.** The
readout is reading expression variance, not the dictionary.

This is the T1.1 random-dictionary must-beat baseline arriving by a second, independent road, and it
arrives with the answer "we do not beat it". The direct T1.1 comparison (`randdict_targets.npz`
through the same instrument) is still running and will either corroborate or contradict this; it is
the decisive check and it is not being waited for before this is written down.

### In plain terms

We scrambled which gene each coordinate of our dictionary refers to and then asked whether the images
could still predict the scrambled coordinates. They could — just as well as before. Two things follow.
The good one is the one we asked for: the *overall pattern* survives label scrambling, so a claim of
the form "there is molecular structure here that morphology can see" is safe. The bad one is the one we
did not ask for: if scrambling the gene labels does not hurt, then the specific biological directions
we spent the dictionary on are not what the images are reading. A randomly chosen set of 128
directions through the same expression data would have done the same job.

### Meaning for the claim

* **P1 (this paper).** The control behaves as a control should, and both halves are reportable. The
  attribution half must carry its by-construction caveat, exactly as the E0-side `_gene_label_shuffle_null`
  already does.
* **P3 — this is the damaging one.** P3's claim is "interventional coordinates beat curated pathway
  scores / are legible in a way generic decompositions are not". On this evidence, at this readout, on
  this representation, **that claim is not supported**: the dictionary is statistically
  indistinguishable from a spectrum-matched random projection. P3 may not assert the dictionary is
  doing work until it can show a comparison in which the fitted basis wins with a CI excluding zero.
* **P4.** `inspect_gene` cannot ship on this evidence. A fluent gene-level answer built over an
  attribution that is invariant to the thing it names is precisely the "launders uncertainty into
  prose" failure P4 forbids.
* **Method note for P1's own credibility:** the containment test in the T1.5 spec ("shuffled value
  inside the bootstrap CI of the unshuffled value") is a weak test — a wide CI passes it by accident.
  The paired difference CI is what should be quoted. Both are reported here so the difference is
  visible rather than chosen.

### Files / commits

`v2/build_shuffled_pbs_targets.py`, `v2/calibra/gene_label_shuffle_control.py` (commit 1c4b4b5).
Results: `p1_evidence/track1/gene_label_shuffle/gene_label_shuffle_control.json`,
inputs `p1_evidence/inputs/pbs_shuffled_seed{1,2,3}.npz` — all under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/` (persistent).
