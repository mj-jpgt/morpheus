> # ⚠ ACTIVE ENTRY POINT IS NOW `PROJECT_GUIDE.md`
> This file's build-phase content is historical — folded into the paper drafts long since. Read
> `PROJECT_GUIDE.md` for current project state, the per-paper build plan, and the methodology rules.
> Unlike this chain, `PROJECT_GUIDE.md` is not superseded by a newer handoff file — it is maintained
> in place.

# HANDOFF — PHASE D: build. **This supersedes the E-series as the active task list.**

**Read in this order:** this file → `HANDOFF_GATES.md` (mandatory; G2 *liveness* matters most here because
Phase D trains) → `v2/research/rebase/ENGINE_CLD.md` (why). Branch `research/rebase-vision`, currently
`e33a853`. `HANDOFF_EXPERIMENTS_NOW.md` (E0–E5) is **closed**.

> **Before any GPU run you MUST complete the recursive audit loop in §2.** Two prior audits each caught a
> defect that would have produced a confident, wrong answer, and both defects passed their own test suites.

---

## 0. CONTEXT UPDATE — what changed since you last had context

### The instrument was broken and is now fixed
CALIBRA's spike readout scored recovery with `top_canonical_correlation`, a **maximum over 16
components**, while the spike lives on **one** known direction. Ambient sits at ~0.97, so every
detection floor came back `NaN` on real data while all 11 synthetic self-tests passed. Three nested
defects (max-readout; partial replacement because `a` was standardised and `y·v` was not; absolute value
taken *before* pairing, which destroys the paired comparison since induced correlation has random sign).

**Consequences you must know:**
- **The confound adjustment does NOT destroy signal** — attenuation 0.94–1.23, i.e. ≈1. This retires the
  objection that killed three earlier theses.
- **Residualising two orthogonal signals through a shared confound design INDUCES correlation between
  them** (0.067–0.140 for the 99-column cancer+TSS design at n=2,530). Nobody reports this. It is why the
  floor is read with a paired test.
- Two floors now exist and are **not interchangeable**: `transmission_floor` (paired, near-noiseless,
  **never quote as a detection limit**) and `detection_floor` (unpaired, conservative, quotable ≈0.2 WSI).

### F2 IS WITHDRAWN
E3 ran. `wsi_identity` changes by **2.6e-04** between `full` and `identity_only`, against **1.4e-01** for
the biology head. The identity head is the **frozen MLP-CLIP teacher passed through**. So "the head
trained for biology is worse at biology" restated *"the frozen teacher beats our biology head"* — a
distillation observation, not a claim about objectives. **Do not cite F2.** Any draft text asserting
"molecular supervision degrades the molecular channel" must be removed until Task D1 provides the arm
that was never trained.

### What replaced it is stronger
Effective rank **−17%** (38.48 → 32.06) with the molecular channel **unchanged** (0.4768 → 0.4748).
Combined with the earlier +107% rank at flat specificity, **rank ≠ information is now demonstrated in
both directions** across independent experiments.

### E0 / E0b are done; the E-series is closed
- E0 as originally written **could not return a negative**: its "matched-spectrum null" ignored the
  spectrum and was a Haar-uniform random subspace (~k/n_genes) while any two real expression matrices
  score 0.4–0.8 from generic co-expression. A zero-shared-biology construction passed every gate. Fixed
  with a **biological control arm** (`energy_test_p_value`: responsive p<0.01 vs non-responsive p>0.5).
- **E0 result:** `supported`, K562, **~10% of the achievable ceiling**, `pc1_share` 0.09–0.21 (not library
  size). The control arm absorbs **55%** of the raw overlap at k=100 — half of what a random-floor rule
  would have credited to perturbation biology is generic expression structure.
- **Replication:** cross-lineage `supported`, 4/4 contexts. At k=25, K562 +0.0387 vs RPE1 +0.0394 — **within
  2%** when n-matched. The n-matching was load-bearing: uncapped, K562 reads +0.0671 vs +0.0496, a **35%
  inflation from sample size alone** that would have been reported as lineage specificity.
- **E0b:** 8,403 perturbations have **effective rank 132.1** (RPE1 113.9), stable rank 17.4, coherence
  0.85. *"11,000 independent causal directions" is dead.* **`n_equivalence_classes` returned n and is
  mis-specified — do not use it.**

### Decision: stop validating, start building
Further E0-phase work has diminishing returns because **E0 does not measure the operation we are
building.** It is a feasibility screen ("can this dictionary express tumour biology at all?" — yes,
partially, replicated). The causal claim is carried by **per-axis certification downstream**, not by E0.
GTEx is **dropped** — it interrogates E0's statistic, not the system.
**The proliferation confound is answered better by Task D2 than by any further E0 control** (see §3).

### Claims are now enforced in code
`v2/calibra/claim_guards.py` + `tests/test_claim_guards.py` (15 tests). Run `validate_claim()` before
any write-up. **E0 is currently an INADMISSIBLE transfer claim**, blocked on `proliferation_deflation`
and `single_platform`, pinned by `test_current_e0_result_is_not_yet_an_admissible_transfer_claim`. If
you discharge one, that test fails and you must update it **deliberately**.

---

## 1. THE TASKS — follow these literally

> **Nothing below runs on GPU until §2 (recursive self-audit) returns clean.** Every step here is
> CPU-verifiable first. A GPU run started on unaudited code is the exact failure this handoff exists to
> prevent.

---

### D1 — The missing arm. *Does molecular supervision help or hurt?*

#### D1.0 — Why you must retrain BOTH arms, not reuse the existing one

The three diagnostic artifacts record only this in `source_manifest`:
`configuration_sha256`, `git_commit`, `git_dirty`, `package`, `package_root`, `source_tree_sha256`.

**No epochs. No learning rate. No token budget. No seed.** And measured:
- `configuration_sha256` **differs** across all three arms (expected — different profile — but it means
  the hyperparameters are unverifiable from a hash);
- `git_dirty` is **True** for all three.

So **G0.4 cannot be discharged from disk.** Do not reuse `diagnostic_programme_only_seed42.npz` as D1's
baseline. **Train both arms yourself from one command that differs only in `--objective-profile`**, so
matching holds by construction rather than by hope.

#### D1.1 — Add the `programme_free` objective profile

Six edit sites. All six are required; five will pass tests while the sixth silently breaks the run.

| # | file:line | change |
|---|---|---|
| 1 | `v2/training.py:47-48` | add `"programme_free"` to the allowed set in `__post_init__` |
| 2 | `v2/training.py:69-72` | new branch in `weights()` (contract below) |
| 3 | `v2/training.py:158-161` | `liveness_parameter_groups` returns **the same groups as `programme_only`**: `("wsi", "rna", "shared", "biology_programme")` |
| 4 | `v2/runner.py:344-354` | `_trained_states_for_profile` returns **exactly** `["wsi_biology", "rna_biology", "full_biology"]` |
| 5 | `v2/runner.py:390` | extend `use_anchor` so `programme_free` **also** skips the MLP-CLIP anchor |
| 6 | `v2/runner.py:563` | add `"programme_free"` to the `--objective-profile` choices tuple |

**Weight contract for `programme_free`** — exact, and asserted in a test before any GPU time:

```
programme            == 0.0     # no Hallmark Gaussian-NLL regression
neighbourhood        == 0.0     # the diagnosed collapse mechanism
supcon               == 0.0     # the diagnosed collapse mechanism
identity             == 0.0
fusion_identity      == 0.0
patient_consistency  == 0.0
semantic             == 0.0
decorrelation        == the SAME value programme_only uses
                        (training.py:69-72 keeps it there; if the arms differ on decorrelation
                         you are measuring decorrelation, not supervision)
<new biology contrastive term>  > 0.0
```

**The replacement objective.** `v2/losses.py:13` already provides
`symmetric_infonce(left, right, temperature=0.07)`. Wire it as **RNA-paired InfoNCE on the biology
view**: positives are `(wsi_biology, rna_biology)` of the same patient, negatives are other patients in
the batch. That is the point — a **rank-preserving contrastive signal** instead of regression onto a
~50-D target, which is what flattened the head.

**Batch-size trap (G2.7).** Uncapped H-Optimus batches hold **B ≈ 1–3 patients**. InfoNCE at B=1 has no
negatives and is a silent no-op — exactly how `feature_decorrelation`'s `min_batch=8` guard disabled
itself on every real step. **Reuse the detached feature-bank ring buffer already built for
`feature_decorrelation`**, or assert the effective negative count per step and fail loudly below 8.

#### D1.2 — CPU verification, before any GPU

```python
from morpheus.v2.training import <ScheduleClass>
free = <ScheduleClass>(objective_profile="programme_free").weights(epoch=99)   # post-warmup
prog = <ScheduleClass>(objective_profile="programme_only").weights(epoch=99)
assert free["programme"] == free["neighbourhood"] == free["supcon"] == 0.0
assert free["decorrelation"] == prog["decorrelation"]      # the arms must differ in ONE thing only
assert free["identity"] == free["fusion_identity"] == 0.0
assert free[<new_key>] > 0.0
from morpheus.v2.runner import _trained_states_for_profile
assert _trained_states_for_profile("programme_free") == ["wsi_biology", "rna_biology", "full_biology"]
```

Then, still cheap:
- **G2.6 overfit-one-batch** on `programme_free` — loss must reach ~0 on a single fixed batch.
  **If this fails, stop.** Nothing downstream is interpretable.
- **G2.3** gradient norm of the `biology_programme` group **non-zero at step 1**. Zero means the new loss
  is detached or its parameters are unregistered.
- **G2.2** log every loss term separately for 3 epochs; the new contrastive term must be
  **> 1e-4 x the largest term**. Below that it is off, whatever the total says.
- **G2.1** `||dtheta||/||theta|| > 1e-2` after a short run.

#### D1.3 — Run

Identical command for both arms, differing **only** in `--objective-profile`. Seeds `42, 43, 44`.
Record the full argv in the run log.

```
python -m morpheus.v2.runner --split-file <SPLIT> --output-dir runs/d1_<arm>_seed<S> \
  --objective-profile {programme_only|programme_free} \
  --epochs 40 --token-budget 32768 --hidden-dim 512 --layers 4 --heads 8 \
  --learning-rate 2e-4 --weight-decay 1e-2 --decorrelation-weight 0.04 \
  --loss-warmup-epochs 4 --seed {42|43|44} --device cuda
```

**Pass no `--mlp-clip-anchor` and no `--mlp-clip-teacher` to either arm.** With an anchor, D1 inherits
the exact artifact that killed F2 and answers nothing.

#### D1.4 — Measure

```
python -m morpheus.v2.calibra.run_calibra --artifacts <the 6 npz> \
  --targets <frozen_rna_targets.npz> --output runs/d1_calibra \
  --n-draws 40 --n-components 16 --n-permutations 2000 --n-jobs 30 --seed 42
```

`--targets` **must** be `frozen_rna_targets.npz`. `data.hallmark` is train-fold-only and **constant on
the test split** — that is where the original `nan` came from (G1.1).

Report per arm: held-out CCA, `effective_rank`, within-cancer specificity
(`honest_metrics.macro_group_pearson`, `control_adjusted_specificity`), and a **paired bootstrap on the
between-arm difference** (`v2/paired_bootstrap.py`). F2's headline gap had no CI, and that is partly how
it survived into a paper draft.

**Pre-register before running:** *`programme_free` >= `programme_only` on the held-out molecular
channel.* If so, molecular supervision is not helping -> F2 is restored as an objective claim and PBS is
motivated. If `programme_only` wins, the collapse story is wrong — **escalate, do not proceed to D2.**

---

### D2 — PBS head-to-head. *Interventional coordinates vs curated pathways.* **Highest value; run this first if choosing one.**

**Why this subsumes the proliferation confound:** Hallmark already contains proliferation programmes. If
our dictionary's only real content is proliferation, it **cannot** beat Hallmark — there is nothing left
to win on. A win means it carries something curated pathways lack.

#### D2.1 — Build the dictionary. Reuse `v2/pbs.py`; do not reimplement.

```python
from morpheus.v2.pbs import ReferenceDictionary, LegibilityOperator, weighted_code_loss
D = ReferenceDictionary.fit(responses=P, genes=genes, atom_ids=atom_ids, n_components=128)
codes = D.encode_expression(expression=tcga_rna, genes=tcga_genes)   # exact gene identity is mandatory
```

- `P` is the K562 perturbation matrix built **exactly as E0 builds it** — reuse
  `v2.calibra.e0_basis_transfer._load_perturbation`. Do **not** write a second loader.
- **`n_components = 128`**, with a sensitivity run at 64 and 256. Rationale: E0b measured **effective
  rank 132.1**. **Do NOT use 8,403 atoms** and **do NOT use `n_equivalence_classes`** — it returned *n*
  and is mis-specified.
- Fit `D` on the **dev cancers only**, then encode all patients. Fitting on everything leaks the test
  split into the supervision target.

#### D2.2 — Two arms, identical but for the supervision target

| arm | target |
|---|---|
| **H** | ~50 Hallmark scores (the current baseline; `programme_only` reproduces it) |
| **I** | the 128-D interventional codes from D2.1 |

Same architecture, epochs, LR, token budget, seeds `42, 43, 44`. **Only the regression target differs.**

#### D2.3 — Report, including the field that answers proliferation for free

Held-out CALIBRA channel, `effective_rank`, within-cancer specificity, standard benchmark — **plus, for
every axis, its proliferation / essentiality loading.** You will have per-axis gene loadings anyway.
*If every legible axis comes back proliferation-loaded, that is the deflation, visible without a
separate experiment.*

**Escalate immediately if I ~= H** (overlapping paired-bootstrap CIs). That means the interventional
dictionary's content already sits inside curated pathways, and the rebase premise is in trouble.

---

### D3 — Purity into the adjustment set. CPU. Run while D1/D2 train.

Open since Phase 1, and a **complete alternative explanation** for the morphology<->molecular channel.

**No TCGA purity table is on disk.** CPTAC has `tumor_purity__washu.parquet` per cohort under
`data/raw/hf_tcga_cptac_cgga/cptac/tables/<cohort>/` — that is for the external cohort later, **not** for
TCGA.

1. **Preferred:** ask for the published TCGA consensus purity (ABSOLUTE / PanCanAtlas). Open access.
2. **Fallback:** compute an ESTIMATE-style stromal+immune score from the expression we hold. **Emit
   `purity_source="expression_derived"`** and state the partial circularity in the write-up: you are
   residualising expression on a quantity computed from that same expression.

Add it to `confound_design` in `v2/calibra/run_calibra.py`, re-run the Phase-1 channel measurement, and
report the channel **before and after**. **If the channel dies when purity enters, that is a finding —
report it, do not bury it.**

---

## 2. RECURSIVE SELF-AUDIT — mandatory, before any GPU run

**No GPU job starts until this loop returns clean.** Two prior audits on this project each found a defect
that would have produced a confident, wrong answer — one where the experiment could **never return a
negative**, one where a **true negative would have been filed as a crash**. Both passed their own test
suites at the time.

**The loop.** Max **2** auditor agents at a time; you are the third.

1. Write the code and get your own tests green.
2. Spawn **one adversarial auditor**. Give it the diff, the tests, this handoff and `HANDOFF_GATES.md`.
   Instruct it verbatim:
   > *"Assume this is broken and find the mechanism. Check specifically: a loss term that is silently
   > off; a guard that no-ops at the real batch size; an arm asymmetry that invalidates the comparison;
   > leakage of the test split into a fit; a scientific outcome wired into a pass/fail gate; and whether
   > this design can return a NEGATIVE result at all. Default to BROKEN if uncertain. Return GO or NO-GO
   > with file:line and a minimal fix."*
3. **If NO-GO: fix, then return to step 2 with a FRESH auditor.** Repeat until an auditor returns GO.
   **Do not argue an auditor into a GO — fix the code.**
4. Optionally spawn a **second** auditor for an independent GO on the fixed code. If the two disagree,
   treat it as NO-GO and iterate.
5. **Only then** run on GPU.

**Termination:** stop when a fresh auditor, seeing the code for the first time, returns GO with no
blockers. If you reach **4 rounds** without a clean GO, **escalate to the mastermind** — that means the
design is wrong, not the implementation.

**"Sizeable and conceptual" means:** anything that changes a number, a verdict, or what a claim licenses.
Ignore style. A dead loss term, a leaked split, an unmatched arm, an ungated NaN, or a design that cannot
produce a negative result are all blockers.

---

## 3. Gates that bite hardest here

Phase D **trains**, so the G2 liveness family applies in full, and every one of these has bitten us:
- **G2.6 overfit-one-batch** — cheapest decisive test. Run first on every new profile.
- **G2.1** `||dtheta||/||theta|| > 1e-2` — `feature_decorrelation` once contributed ~1e-6 with a 4e-5
  parameter delta.
- **G2.2** every loss term logged **separately**; below 1e-4 x the largest it is effectively off.
- **G2.7** guards must fire at the **real** batch size (B ~ 1-3, not 8).
- **G4.1 positive control** (`rna_*` -> RNA targets) must pass **in the same run**, or you measured
  nothing.
- **G1.1** no constant columns **on the evaluated split**.
- **G3.1** `calibra.spectral.effective_rank` — **singular values**, never covariance eigenvalues (a 6x
  error that reached a paper draft).
- **G0.2** clean worktree before launching, or provenance fails and the output is quarantined.

## 4. Do not
- Do not cite F2, or quote the marginal `bootstrap_ci95` as a 95% CI (biased low by 0.04-0.09, flagged
  `bootstrap_ci95_is_biased_low`; decisions use the **paired** difference).
- Do not reuse `diagnostic_programme_only_seed42.npz` as D1's baseline (see D1.0).
- Do not pass an MLP-CLIP anchor to either D1 arm.
- Do not use `n_equivalence_classes`, 8,403 atoms, or the phrase "11,000 independent causal directions".
- Do not fit the dictionary on all patients — that leaks the test split.
- Do not compare gap magnitudes across E0 runs with different n or q.
- Do not claim E0 shows the alignment is *biological* — both lineages share one platform and one
  effect-size-monotone statistic.
- Do not run further E0 controls (GTEx, second platform) — superseded by D2.

## 5. Escalate immediately
1. G4.1 positive control fails, or G2.6 overfit-one-batch fails.
2. **D2 returns I ~= H** — the premise is in trouble.
3. **D1 returns `programme_only` > `programme_free`** — the collapse story is wrong; do not proceed to D2.
4. The new contrastive term is below 1e-4 x the largest loss term and you cannot fix it.
5. The channel vanishes when purity enters (D3) — a finding; stop and report before building on it.
6. **4 audit rounds without a clean GO.**

## 6. Logging and provenance
Gates -> `v2/research/rebase/nature/GATE_LOG.md`. Results **including negatives** ->
`v2/research/rebase/nature/EXPERIMENT_LOG.md`. Record device, wall-clock, library versions and the full
argv per run. Commit and push per task once gates **and** audit both pass. A missing gate row is a FAIL.

**GPU:** required for D1/D2 — the first thing on this project that genuinely needs it. D3 is CPU. Ask for
the SSH login when you need it. Launch only from a **clean worktree**.
