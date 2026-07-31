# HANDOFF — HEALTH GATES (run these or your results mean nothing)

**Purpose: stop us concluding "the science says no" when the truth is "the implementation was broken."**

Every failure mode below has *already happened once on this project*. A null result produced by a
scrambled join, a constant column, a silently no-op loss term, or a stale file is indistinguishable
from a real null **unless you check** — and we have twice nearly written up a technical fault as a
scientific finding.

## THE GOVERNING RULE

> **A negative result is reportable only if the positive control passed in the same run, on the same
> data, through the same code path.**
>
> If the positive control did not pass: **you have not measured anything.** Do not log a finding,
> do not proceed, fix the pipeline. If a gate fails and you cannot explain why within one hour of
> work, **escalate — do not work around it.**

Record every gate outcome in `v2/research/rebase/nature/GATE_LOG.md`, one row per
`(experiment, gate, value, threshold, PASS/FAIL)`. **A missing gate row is treated as a FAIL.**

---

## G0 — Provenance (before *any* experiment; ~1 min)

| Gate | Check | Fail action |
|---|---|---|
| G0.1 Workspace is real | working dir is a **junction**, not a copy (`Get-Item ... .LinkType -eq 'Junction'`) | STOP. This silently validated 62-line stale code against a 107-line source once already (mistake #2). |
| G0.2 Code identity | record `git rev-parse HEAD` + dirty-file list in the gate log | STOP if dirty and unexplained |
| G0.3 Artifact identity | SHA-256 (first 16 hex) + mtime + size of every input `.npz`/`.h5ad`/`.tsv` | STOP on mismatch with a prior run |
| G0.4 Manifest read | for each artifact, print `manifest_json` epochs / LR / seed / token budget | if runs being compared differ → the comparison is **suggestive, not causal**. Say so in the log. |

## G1 — Data integrity (per data source; catches the silent-nan class)

| Gate | Check | Threshold |
|---|---|---|
| **G1.1 No constant columns** | per-column std on **the split actually being evaluated** | **zero columns with std < 1e-8.** *This exact bug produced `nan` specificity: `data.hallmark` is train-fold-only and is constant on the held-out split.* |
| G1.2 Join is a real join | match on **IDs**, never positional index; log `n_left, n_right, n_matched` | `n_matched ≥ 0.8·min(n_left,n_right)`, else investigate |
| G1.3 Join is not scrambled | shuffle-detection: re-run G4.1 positive control | see G4 |
| G1.4 Scale sanity | max value, fraction of exact zeros, skew | flag if TCGA looks un-logged (max > 1e4) — PC1 will be library size |
| G1.5 No all-zero / all-NaN rows | count and drop explicitly, log the count | never drop silently |
| G1.6 Gene symbol mapping | overlap count after `SYMBOL\|ENTREZ` split; assert a housekeeping panel (`ACTB, GAPDH, TUBB, RPL13A, B2M`) is present in the intersection | if housekeepers are missing the mapping is broken, not the biology |

## G2 — Liveness (anything that trains: E2, Milestone D)

**These catch the dead-implementation class. `feature_decorrelation` once contributed ~1e-6 to the loss
with a parameter delta of 4e-5 — it was "working" in the sense that it ran, and doing nothing.**

| Gate | Check | Threshold |
|---|---|---|
| **G2.1 Params actually moved** | `‖θ_final − θ_init‖ / ‖θ_init‖` | **> 1e-2.** Below that the module is a no-op regardless of what the loss curve says. |
| **G2.2 Every loss term is live** | log each term's **magnitude separately**, per epoch | any term < 1e-4 × the largest term is **effectively off** — flag it loudly, do not average it into a total |
| G2.3 Gradients reach every group | per-parameter-group grad norm at step 1 and at the end | no group at exactly 0.0 (that means detached or unregistered) |
| G2.4 Loss actually decreased | train loss final vs init | ≥ 20% reduction, else undertrained — a "collapse" observed here is an artifact |
| G2.5 Not undertrained | loss still falling at the last epoch? | if yes, extend before drawing any conclusion about the converged state |
| **G2.6 Overfit-one-batch** | can the head drive loss → ~0 on a single fixed batch? | **must succeed.** If it cannot memorise 1 batch, it cannot learn anything, and every downstream number is meaningless. This is the cheapest and most decisive liveness test there is. |
| G2.7 Guard clauses fire on real shapes | assert any `min_batch`-style guard is exercised at the *actual* batch size | uncapped H-Optimus batches hold B≈1–3; an `min_batch=8` guard silently disabled the term on every real step |

## G3 — Representation health (every latent, every experiment)

| Gate | Check | Threshold |
|---|---|---|
| **G3.1 Effective rank** | `calibra.spectral.effective_rank` — **singular values, NOT covariance eigenvalues** | report always. σ² inflates collapse ~6× and that error reached a paper draft (mistake #1). |
| G3.2 Dead dimensions | fraction with variance < 1e-6 × mean variance | report; > 0.5 means half the head is off |
| G3.3 Duplicate dimensions | fraction of pairs with \|corr\| > 0.99 | high ⇒ apparent width is fake width |
| G3.4 Not sample-constant | does the latent vary *across samples* at all? | std over samples > 0 per dim, else the encoder is returning a bias vector |
| G3.5 Not batch/site-degenerate | R² of latent on site/scanner one-hot | report; if ≈1 the "representation" is a site code |
| G3.6 Norm sanity | mean/median L2 norm; NaN/Inf count | **zero** NaN/Inf. Any NaN = stop. |

## G4 — Statistical validity (the false-conclusion guards)

**This is the section that decides whether a null is real.**

| Gate | Check | Required outcome |
|---|---|---|
| **G4.1 POSITIVE CONTROL** | `rna_*` states → RNA-derived targets (circular by construction) | **MUST be strong.** If RNA→RNA is weak, the join/pipeline is broken, not the model. **This is the single most important gate in this file.** |
| **G4.2 NEGATIVE CONTROL** | shuffle the patient pairing between X and Y, re-measure | **MUST collapse to null.** If a shuffled pairing still scores, you are measuring capacity, not signal. |
| G4.3 Permutation null is sane | null distribution centred near the theoretical chance level, not degenerate; report its spread | a null with ~0 variance means the permutation isn't permuting |
| G4.4 Held-out ≈ in-sample | report both; large gap = capacity inflation | in-sample CCA on pure noise reads > 0.3 — never report it alone |
| G4.5 Null resolution | permutation p is floored at 1/(n+1) | raise n for any headline number; never write "p < 0.001" from 100 draws |
| G4.6 Effect size, not just significance | report the effect and its CI alongside every p | — |

## G5 — Experiment-specific gates

### E0 (basis transfer) — needs a false-**positive** guard as much as a false-negative one

| Gate | Check | Why |
|---|---|---|
| **E0.a Delta, not absolute** | confirm whether the stored matrix is already control-subtracted **before** subtracting | if it isn't a delta, "perturbation directions" are dominated by mean expression level and the alignment is trivially high |
| **E0.b Strip the trivial axis** | recompute principal angles after removing the mean-expression axis / PC1 from **both** matrices | **if all alignment lives in PC1, you have discovered library size.** This is the most likely false positive. |
| E0.c Orientation | assert `P` is (perturbations × genes) and TCGA is (samples × genes) after transpose; print both shapes | a transpose error compares the wrong axes and reads as null |
| **E0.d Ceiling control** | principal angles of TCGA vs **itself, split in half** | gives the achievable maximum. Without it, "alignment = 0.4" is uninterpretable. |
| E0.e Floor control | the matched-spectrum random-rotation null, ≥100 draws | already specified; verify the rotation preserves singular values |
| E0.f Cross-context replication | repeat with RPE1 | agreement across two unrelated lineages ⇒ biological, not line-specific |

**E0 verdict is only valid with (d) ceiling, (e) floor, and (b) trivial-axis stripping all present.**
A number between an unmeasured floor and an unmeasured ceiling is not a result.

### E1 (is the added rank empty?)
- Report effective rank **and** the above-floor-direction count **from the same run** — the whole claim
  is that these two move differently. Computing them from different runs invalidates it.
- Detection floor must be finite; if `nan`, the spike design is unpaired (fixed once already — the
  paired design measures per-draw increments over that draw's own level-0 baseline).

### E2 (expressible intersection) — needs a capacity control
- **E2.a:** at `k = full expressible dimension`, the head **must** reach high effective rank.
  If it cannot, the sweep measures optimisation failure, not the hypothesis. **Run this first.**
- **E2.b:** weight-decay control — repeat the sweep at wd=0. If collapse tracks wd rather than *k*,
  the mechanism is regularisation, and the H2 claim is dead. Report that.
- G2.1–G2.6 apply in full at every *k*.

### E3 / E4 (channel measurement on frozen artifacts)
- G1.1 (constant columns) is **mandatory** here — this is precisely where the `nan` came from.
- G0.4: if `programme_only` and `identity_only` differ in epochs/LR/budget, the ablation is confounded.
  **Report it; do not quietly compare them anyway.**
- E4 must report the **whole ladder** including `raw_hoptimus_meanstd`, not just the winner.

### E5 (multivariate vs univariate)
- The frozen null file has **1 draw**. Regenerate at 20 (`build_matched_random_controls`).
  Using the frozen file is a silent invalidation.
- Multivariate CCA is a **maximum**, per-target Pearson is a **mean**. Report both; never conflate.

---

## Escalate immediately (do not work around, do not proceed)

1. **G4.1 positive control fails** — pipeline broken; every other number is void.
2. **G4.2 negative control passes** — you are measuring capacity, not signal.
3. **G2.6 overfit-one-batch fails** — the model cannot learn; no result is interpretable.
4. Any NaN/Inf in a latent (G3.6).
5. **E0 at null** — the engine dies; we pivot to reporting the transfer failure.
6. **E3 shows `programme_only` has a strong biology channel** — F2 was an anchoring artifact and PBS
   loses its primary evidence.
7. Any gate you cannot explain within ~1 hour.

## Self-audit protocol (**maximum 3 agents total, including yourself**)

You are agent 1. After each experiment produces a result **and its gates pass**:

- **Agent 2 — adversarial auditor.** Spawn ONE. Its job is to **REFUTE the finding, not confirm it.**
  Give it: the result, the gate log, the code path, and this instruction —
  > *"Assume this result is an artifact. Find the mechanism. Check specifically: stale inputs,
  > constant/degenerate columns, scrambled joins, transposed matrices, a loss term that is
  > effectively off, in-sample inflation, a degenerate null, and whether the positive and negative
  > controls actually ran on this exact data. Default to 'artifact' if you are uncertain."*
  It must return a verdict of **ARTIFACT / INCONCLUSIVE / SOUND** with the evidence for it.
- **Agent 3 — tiebreak only.** Spawn only if agent 2 returns ARTIFACT/INCONCLUSIVE and you disagree.
  Give it both positions and no hint of which is yours.

**Never exceed 3 concurrent agents.** If agent 2 says ARTIFACT and agent 3 agrees, the result does
**not** go in the experiment log as a finding — it goes in as a defect, and you fix it.

## GPU

- **E0 / E0b: use it.** SVD of 11,258×8,248 plus ≥100 matched-spectrum null draws is hours on CPU and
  minutes on GPU (`torch.linalg.svd` on device; randomized SVD is fine for top-*k* subspaces —
  state which you used).
- **E1 / E3 / E4 / E5: CPU is genuinely fine** (2,530 × 256 frozen embeddings — seconds).
- **E2: CPU fine**, but GPU makes the *k*-sweep × seeds grid cheap enough to run properly.
- **Milestone D (the PBS retrain): GPU required.** That is where the paper's win lives.

Record device, wall-clock, and library versions per experiment. A result that cannot be reproduced on
a stated device with stated versions is not yet a result.
