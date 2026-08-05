## 2026-08-05 00:40 UTC — PREDECLARATION: the labelled linear probe P2 §2.5 marks `[STILL NOT MEASURED]`, and how its agreement or disagreement with rank will be read

**Logged:** 2026-08-05 00:40 UTC (clock checked against the box, `150.136.45.194`, which agrees).
**Committed before any probe statistic is computed.** No probe number for any D1 or D2 artifact
exists at the time of writing except E1's four TP53 readings
(`e1_clinical_endpoints_from_morphology_20260804T1930Z.md` §3b), which are quoted below **as the
prior that fixes my expectation**, not as a result of this run.

---

### 1. The gap this closes, in the draft's own words

`paper/P2_RANK_DRAFT.md` §2.5 ends:

> `[STILL NOT MEASURED]` — a labelled linear probe on every artifact, which is the ground truth LiDAR
> and RankMe were validated against. Ours is a held-out canonical correlation against unsupervised
> molecular targets (§3.2), which is a different reference standard and is the one Zaiem et al. would
> attack.

§6.2 carries the same row (*"A labelled linear probe on every artifact | Not run"*), and §6.3 concedes
Zaiem et al. against us. The referee-facing form of the objection is short: **a paper about whether a
label-free proxy tracks usable information has never measured the labelled thing the proxy is a proxy
for.** This run measures it.

### 2. What will be run — fixed here so it cannot be chosen afterwards

**Compute.** CPU only, on `150.136.45.194` because the artifacts and the label table live there and
nowhere else (no `.npz` and no `.parquet` exists in the local checkout). **No retraining, no GPU** —
every probe is fitted on frozen exported embeddings. Threads capped to 1 per
`operational_shared_box_rules_20260804T0730Z.md`.

**Artifacts.** The same twelve §4.1–§4.7 rest on — D2 `{H,I}` × seeds 42/43/44, D1-B `{P,F}` ×
seeds 42/43/44 — **plus** the five `~/e0_run/d1_envelope/rep{1..5}.npz` same-seed retrains, which are
what §4.1's 3.295× rank floor and 1.055× channel spread are measured on.

**Views.** All three co-trained views §4.5(c) exploits: `wsi_biology`, `rna_biology`, `full_biology`.

**Cohort.** `split == "test"` in the artifact (2,766 patients), intersected with label availability —
`run_calibra.main`'s own cohort rule, the same one `p2_competing_metrics.py` and
`p1_cancer_type_certificate.select_cohort` use.

**Labels — reused, not invented.** Two tables that already exist on this project:

| probe | label | source | n held out |
|---|---|---|---|
| **A (primary)** | cancer type, 21 test classes | `cancers` inside each artifact | 2,766 |
| **B (secondary)** | `mut_TP53` | `…/e1_endpoints/inputs/e1_endpoint_labels.parquet`, built by E1 from the MC3 public MAF | 2,686 |
| B panel | `grade_high`, `stage_late`, `mut_ATM`, `mut_KMT2D`, `mut_ARID1A` | same parquet | 1,209 / 2,061 / 2,686 |

**Probe A is nominated as the primary because it is the literal analogue of the reference standard
RankMe and LiDAR validated against**: a linear classifier fitted on frozen features, top-1 /
balanced accuracy, on the **raw** representation with no confound adjustment — RankMe's ImageNet
linear probe, transposed. Probe B is the secondary because it is matched to the construct *this
paper* argues about (the molecular channel) and is confound-adjusted the way §3.2's channel is.
**When A and B disagree, both are reported and neither is nominated the winner** — that disagreement
is itself the §4.6a lesson repeated on the labelled side, and it will be written that way.

**Estimators — canonical functions, imported, nothing reimplemented.**

* Probe A: multinomial logistic regression, out-of-fold, folds from
  `calibra.confound_certificate._stratified_folds`, classes from `_encode`, scored with
  `calibra.confound_certificate.balanced_accuracy`. **Secondary estimator on the identical folds:
  `calibra.confound_certificate.lda_oof_balanced_accuracy`**, because
  `p1_cancer_type_pair_withdrawn_20260805T0230Z.md` measured this statistic moving from 0.16 to 0.73
  across estimators on one cohort, and a single estimator cannot be trusted to carry an ordering.
* Probe B: `calibra.known_covariate_control.evaluate_known_covariate`, **unmodified**, with
  `residualise=True` — out-of-fold ridge scores, size-weighted within-cancer AUROC as the primary,
  pooled AUROC reported beside it, 1,000-draw bootstrap CI, 1,000 within-cancer label permutations.
* Residualisation, where used: `calibra.residualise.confound_design` /
  `cross_fitted_residuals` / `pooled_tissue_source_site`, `min_site_count=10` — the identical
  99-column cancer + pooled-TSS design.
* Rank and channel, where re-quoted: `calibra.spectral.effective_rank` /
  `top_canonical_correlation`. **No statistic is computed inline.** Re-deriving a statistic where an
  import belonged is the repeated defect on this project and it produced the withdrawn P1 pair.

**Nothing is compared against zero.** Every probe number is read against three things it must clear:
(i) its **measured** permutation null (`null_p95`; label permutation for A, within-cancer label
permutation for B), (ii) chance = 1/21 for A's balanced accuracy, and (iii) **its own same-seed
retraining floor**, computed on the five `d1_envelope` repeats with the identical estimator, view and
statistic — the exact quantity §4.1 demands of rank, now demanded of the probe.

**A between-arm probe difference counts as `resolvable` only if it exceeds that five-repeat floor.**
A difference inside the floor is `UNRESOLVED` and may not be scored as an agreement *or* a
disagreement, in either direction.

**Must-fail control, run in the same command.** Probe A on the **adjusted** block must collapse to
near chance. If cancer-type balanced accuracy survives cancer + TSS residualisation at well above
1/21, the adjustment is not doing what §3.2 says it does and **every downstream number in this run is
void** — reported as such rather than worked around.

### 3. What I expect, written before the numbers exist

**Probe B (molecular) will most likely resolve nothing.** E1 already measured within-cancer TP53
AUROC on four of the twelve artifacts: `d2_h_seed42` 0.5912 [0.5618, 0.6274], `d2_i_seed42` 0.5827
[0.5566, 0.6307], `d2_h_seed43` 0.6017, `d2_h_seed44` 0.5847. The one between-arm pair already
visible there is **H42 − I42 = +0.0085 against CIs about 0.065 wide**. Extrapolating, I expect the
six between-arm TP53 differences to be of order 0.01 and the probe's own retraining floor to be of
the same order or larger. **My modal prediction is that Probe B agrees with the channel in
*direction* on most pairs and resolves none of them.**

**Probe A (cancer type) will resolve differences and I expect it to agree with the channel less
often than rank does.** Balanced accuracy over 21 classes at n = 2,766 has small sampling noise, so
arm differences should clear their own floor. But cancer type is precisely the direction the channel
residualises **out**, so there is no reason for the two orderings to coincide, and a disagreement
here is the outcome I think most likely.

**On rank, the comparison is fixed by the published tables.** §4.5(c): the channel orders every pair
H/H/H and P/P/P on all three views; canonical R1 orders them H,I,H / P,P,P on `wsi_biology` (5/6
agreement) and I,I,I / P,P,P on `rna_biology` and `full_biology` (3/6). Those are the numbers the
probe's own record will be set against, unchanged.

### 4. The four outcomes and how each is read

| outcome | reading — fixed now |
|---|---|
| **A. Probe agrees with the channel, and resolves the pairs.** | The paper's ground truth is corroborated by the standard RankMe and LiDAR were validated against; §2.5's `[STILL NOT MEASURED]` closes **favourably**; rank's record against the probe is essentially its record against the channel. **This is the paper-supporting outcome and §5's distrust checks apply to it in full.** |
| **B. Probe resolves nothing.** | The reference standard is **no more decisive than the proxy** on this stack. §2.5 closes as *"measured, and it does not adjudicate"*. The paper may **not** then claim that "the channel is resolvable where rank is not" is corroborated by a labelled probe, and §4.1's asymmetry sentence must say the probe was measured and could not confirm it. Modal expectation for Probe B. |
| **C. Probe disagrees with the channel on ≥ 1 pair.** | **HEADLINE, reported first.** §4.6a showed every OK/MISS count moves when the target *block* changes; this would show it also moves when the *reference standard* changes, and moves against the standard with the strongest external warrant. §4.6, §4.6a and §6.3 all need it. |
| **D. Probe agrees with RANK where it disagrees with the channel.** | **Worst case for the paper, reported first and most prominently.** On any such pair the paper has scored rank a MISS against a ground truth the labelled standard contradicts, i.e. rank was right and the exam was wrong. This must be stated in §4.6 and §6.3 in those words, not as a footnote. |

### 5. What would make me distrust a favourable result

A "favourable" result here means one that supports P2's negative conclusion about rank — outcome A
with rank scoring poorly, or any outcome in which the probe vindicates the channel. **Seven checks,
all of which will be reported whether or not they bite:**

1. **The probe's between-arm difference is inside its own five-repeat retraining floor.** Then the
   agreement is *unmeasured*, not confirmed, and must be written `UNRESOLVED`. This is the same trap
   §4.1 caught rank in and there is no reason the probe is exempt.
2. **The probe's observed value is inside its measured permutation null.** A probe at chance orders
   arms by noise; its agreement with anything is worth nothing.
3. **The agreement is carried by the pooled AUROC rather than the within-cancer AUROC.** E1 measured
   pooled 0.73–0.75 against within-cancer 0.47–0.56 for BRAF/KRAS/APC — lineage identification
   wearing a molecular label. Any probe ordering that survives only pooled is reported as a lineage
   artifact.
4. **The agreement holds on one endpoint of the panel and not the others.** One endpoint out of six
   is a coordinate choice, which is exactly §4.6a's finding restated.
5. **The ordering flips between the logistic and the LDA estimator on the same folds.** Zaiem et al.
   is cited in §6.3 for precisely this; if our own two estimators disagree we have reproduced the
   critique against ourselves and must say so.
6. **`rna_biology` / `full_biology` probe values on molecular endpoints are partly circular** — an
   RNA-derived view predicting a molecular endpoint is not an image→molecular measurement. Declared
   now, exactly as §4.5(c) declares it for the channel; those two views' Probe B numbers are
   reported and are **not** used to adjudicate any morphology claim.
7. **The must-fail control does not fail** (cancer-type accuracy surviving adjustment). Then nothing
   in this run is interpretable.

**And one asymmetry I am fixing now so it cannot be chosen later:** if the probe *disagrees* with the
channel, I do not get to retreat to "the probe is a worse standard". Probe A's estimator, folds, null
and floor are all declared above, before any number; if it disagrees, the disagreement is the result.

### 6. What this run cannot do

* It is one stack, one cohort (TCGA), one architecture family — `claim_guards.no_external_cohort`
  stays undischarged and nothing here discharges it.
* The five-repeat probe floor inherits every limitation §4.1 states of the rank floor: one arm
  (`programme_only`, the stable arm), one seed, one configuration, one stack. **It is a floor twice
  over, and will be quoted as "floor", never "envelope".**
* Probe A is not a molecular measurement and will never be described as one.
* Six pairs cannot support a rate (§3.6 rule 3). No count produced here is quoted as one, in either
  direction.

### Files

* Runner (to be written): `v2/research/rebase/p2/p2_labelled_probe.py`, vendored in the repository
  before it is run, per §6.2's traceability requirement.
* Output: `P2_LABELLED_PROBE.json`, brought back and summarised in a `NOTEBOOK_ENTRIES/` result entry.
* Nothing in `paper/P2_RANK_DRAFT.md`, `claim_guards.py`, `claim_evidence.json` or any other agent's
  `PREDECLARED_*` file is touched by this work.
