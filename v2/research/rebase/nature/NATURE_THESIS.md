# NATURE_THESIS — the decisive proposal

**Chief scientist verdict · 2026-07-29**
**Inputs read:** `thesis_T1_spatial_target.md`, `thesis_T2_virtual_perturbation.md`, `thesis_T3_denovo_discovery.md`, `thesis_T5_clinical_actionability.md`, `kill_priorart_{T1,T2,T3}.md`, `kill_feas_{T1,T2,T3}.md`, `data_spatial.md`, `data_external_wsi_outcome.md`, `data_perturb_drug.md`, `nature_bar.md`, `discovery_standard.md`.

---

## 0. THE HONEST HEADLINE, STATED FIRST

**No thesis in the set clears *Nature* main, and none of the four survives both assassins as written.** T1 and T2 and T3 were each killed on feasibility; T5 was killed on data availability and scooped three times over. Anyone who tells you otherwise is reading the scout reports and skipping the kill reports.

But the kill reports are not uniform. Read them as a group and a pattern falls out: **every kill landed on the *deliverable*, and none landed on the *instrument*.** Three separate assassins independently converged on the same missing object:

- `kill_feas_T3.md` §4 — *"there is by construction no signal of known ground-truth strength left with which to demonstrate that the pipeline still works… the authors have not demonstrated that their adjustment preserves signal known to be present; their null is therefore uninterpretable."*
- `kill_feas_T2.md` §1 — *"a rigorous null requires a demonstrated-sensitive instrument; there is no instrument."*
- `kill_feas_T1.md` §1 — *"the experiment cannot produce evidence either way."*

Three theses died of the same disease: **an uncalibrated measuring instrument.** Nobody in this field can tell you what effect size their confound adjustment would have missed. That is the hole. It is a real hole, it is cheap to fill, filling it is a methods contribution in its own right, and — decisively — **it converts the exact experiments that were just killed into experiments that work.**

That is the thesis. One primary direction, below. Ranked fallbacks in §6.

**Reachable venue, stated plainly and up front: *Nature Methods* or *Nature Biomedical Engineering* for the guaranteed core; *Nature Cancer* / *Nature Medicine* if and only if the biology leg (§3, Component 3) lands. *Nature* main is not reachable from here and should not be the target.** Justification against the calibrated bar in §2. The user's own criteria list "new evaluation paradigm" as Nature-tier-acceptable; this satisfies that criterion, and Nature Methods / Nature BME is where that currency is actually spent (`nature_bar.md` D1, C2, C3).

---

## 1. THE CHOSEN THESIS

### 1.1 The name

**CALIBRA** — *Calibrated Cross-Modal Biological Recovery Audit.*
Released artefact: an open harness, `morpheus-cal`.

### 1.2 The claim, as a single falsifiable proposition about biology (G0)

> **In human tumours, the fraction of measured H&E→molecular performance attributable to a non-biological identity term — cohort and submitting-site identity for bulk targets, per-slide identity for spot-level targets — exceeds the fraction attributable to patient- or spot-specific molecular state; and the residual biological capacity of the channel is small, encoder-invariant, and target-modality-invariant, at a magnitude we report against a calibrated detection floor.**

Note what is in the sentence and what is not. There is no model noun. It survives the G0 test in `discovery_standard.md`: it would be equally interesting produced by a competitor's model. It is a statement about *the data relationship*, which is what G3 demands once method-invariance is observed.

**The single falsifier, pre-registered:** any encoder that, under the fixed adjustment protocol, exceeds the calibrated tie band by a margin whose 95% CI excludes zero, on ≥2 independent cohorts. If that happens the paper inverts to a positive result — *"here is the encoder property that carries biological capacity"* — which is a better paper. **Both branches publish. That is the point of choosing this one.**

### 1.3 Why it survived the prior-art assassin

Each assassin returned `killed = false` on novelty for the object CALIBRA measures, and the negatives are hard, not fuzzy:

- **No recovery-calibrated cross-modal claim exists.** `kill_priorart_T3.md` §3: arXiv `histology AND transcriptome AND confounder` = **0 results**; PubMed `histology + deep learning + gene expression prediction + tumor purity + confounding` = **0 results**; PubMed title-restricted (pitfall|caveat|confound|critical assessment|limitations|benchmark) × (gene expression|transcriptome) × (histology|pathology|WSI|image) = **exactly 1 hit** (Wang et al., Nat Commun 2025, PMID 39934114 — verified real, already cited, clean). No paper reports a confound-adjusted scalar effect size for the tumour morphology↔transcriptome channel; no paper reports a dimension *k* against a permutation null.
- **No mean-subtracted, composition-adjusted, or random-panel-nulled spot-level leaderboard exists.** `kill_priorart_T1.md` §1: all five kill conditions K1–K5 unfound; K5 (bulk-vs-spatial ceiling comparison) returned **zero hits on every API**. The clincher is §2.7 — the *Briefings in Bioinformatics* 2026 survey (DOI 10.1093/bib/bbag255) catalogues **46+ H&E→ST models 2020–2025 and discusses none of the four controls**. A whole-subfield survey missing all four is close to proof.
- **The strongest collision, DECAT (arXiv:2605.31504, TCGA n=8,979), is null-referenced but not recovery-calibrated.** It classifies whether a *representation* for a *given task* is confounded. It does not, and cannot, tell you what effect size it would have missed. CALIBRA's recovery curve is precisely the thing DECAT lacks, and DECAT becomes a cheap post-hoc baseline we run and pass rather than a competitor.
- **The field is circling the confound without converting it.** HESCAPE (arXiv:2508.01490) named batch effects; CHRep (arXiv:2604.21573) named "slide-level appearance shifts and regression-driven over-smoothing"; COAST, HEXST and HiST all ship *training-time fixes* for mean-dominance. Five independent groups building workarounds for a quantity nobody has measured is the strongest possible evidence that measuring it is both novel and wanted. **We must cite all five up front as corroboration, not pretend the field is blind.** That is the difference between a landed audit and a desk reject.

### 1.4 Why it survived the feasibility assassin

Take each fatal obstacle in turn and show it is *designed out*, not argued away.

| Kill | Where it landed | How CALIBRA removes it |
|---|---|---|
| **"Full confound control and power are mutually exclusive"** (`kill_feas_T3` §1). Grade+stage exist for only **2,077 cases / 8 cancer types**; RIN missing for **9,800/11,428 (85.8%)**; held-out slice r_min = 0.081 > target 0.07. | T3's certificate demanded grade+stage+RIN. | **We do not require them.** The kill report itself concedes (§8) the full-coverage stratum — cancer type + TSS site + purity + deconvolution on n=6,192 — is **powered at r_min = 0.036** and "costs hours." CALIBRA lives in that stratum by design. Grade/stage/RIN become a *reported coverage limitation with a sensitivity analysis on the 2,077-case subset*, not a load-bearing certificate. |
| **"Uncalibratable instrument"** (`kill_feas_T3` §4) — the deepest objection, explicitly n-independent and "unfixable." Every known-strong axis is inside the residualisation set, so a flat spectrum is indistinguishable from over-residualisation. | T3 had no positive control that survived its own adjustment. | **This is what CALIBRA is.** Two mechanisms: (a) synthetic **spike-recovery** at known r_true through the *identical* pipeline including residualisation, yielding an empirical recovery curve and detection floor; (b) **held-out-confound positive controls** — deliberately exclude one morphologically-legible covariate (MSI, TP53, consensus subtype) from the adjustment set and require the pipeline to recover it at its known strength, rotating which is held out. The reviewer's one-sentence rejection is answered with a figure. |
| **"Residual scored against DepMap is identically zero by construction"** (`kill_feas_T2` §1) — mean-zero-in-patient quantity vs constant-in-patient target. | T2 tried to score a *per-patient* residual against a *lineage-constant* target. | CALIBRA never does this. Perturb-seq/DepMap are used **only** for G8 pre-specified *directional* concordance on *programs* — "perturbing X moves module M in direction d" — with a random-gene-set null, exactly as `discovery_standard.md` G8 specifies. That estimand is not mean-zero-in-patient and is scoreable. |
| **"Data-processing inequality makes it information-free"** (`kill_feas_T2` §2) — RNA→dependency is a fixed patient-independent map, so I(WSI;dep\|type) ≤ I(WSI;RNA\|type) = the already-known +0.07. | T2's headline number was a coordinate change on a number already in hand. | CALIBRA's headline number is **not** the +0.07. It is the **decomposition ratio** N:I:C:B and the **calibrated floor** — quantities the project does not have and nobody has published. The DPI argument does not apply because we are not passing the signal through a lossy fixed map; we are partitioning the signal we already measure. |
| **"Composition is unmeasurable; CellViT vs 1536-d FM is capacity-confounded; deconvolution is circular"** (`kill_feas_T1` §1). | T1 compared a ~5–15-d nuclei vector against a 1536-d embedding and called the gap biology. | **Capacity-matched baselines.** Report B as a *curve in effective dimension d*: PCA-truncate the FM embedding down to d, and random-feature-lift the composition vector up to d, and compare at matched d. Add the **per-slide-mean predictor (zero free parameters)**, which is immune to capacity confounding and is the most damning baseline available. Deconvolution-based composition is reported **separately and labelled as an upper bound**, never as the primary adjustment, precisely because it is a linear functional of the target. |
| **"Target-invariance is a curve-fit with four analyst knobs"** (`kill_feas_T1` §2). | T1 asserted numerical equality between two statistics with different sampling units and no common denominator. | **We assert no numerical equality.** We assert three *structural* properties, each measured in its own native metric space: (i) after the appropriate zero-parameter naive baseline is subtracted, the encoder ranking collapses inside the calibrated tie band; (ii) the identity term exceeds the residual biological term; (iii) both hold in both target modalities. One fixed, pre-registered adjustment protocol, declared before any HEST number is computed. The garden of forking paths is closed by pre-registration, not by argument. |
| **"Unmeasured attenuation — TCGA slide provenance is uncontrolled"** (`kill_feas_T3` §5). 30,326 slides = **18,425 Tissue (frozen, adjacent to the RNA aliquot) + 11,901 Diagnostic (FFPE, different block)**; repo grep finds no `slide_type` field anywhere. | T3 sold a calibrated number without the attenuation constant. | **This kill converts into a measurement.** TS slides bracket the RNA aliquot; DX slides come from a different block. **The TS-vs-DX delta, within TCGA, on the same patients, *is* the tissue-mismatch attenuation constant.** It is a GDC metadata join onto an existing store — hours, no re-extraction. GTEx (same-sample image+RNA, 25,306 samples) provides the independent upper anchor. Nobody has published this constant. It is fallback paper F3. |
| **"Disk: 48.2 GB free of 927 GB, 95% full, OneDrive-synced"** (`kill_feas_T1` §4). | Operational. | **Verified again this session: `C: 927G, 879G used, 48G avail, 95%`.** Hard blocker. Buy ≥8 TB external before anything else (§7, item 0). Not a scientific issue; do not let it become one. |

**What no assassin disputed, and what CALIBRA is built on:** the 6,192-patient paired cohort is real and on disk; H-Optimus patch features are extracted and uncapped; compute is "a rounding error on one A100" (`kill_feas_T3` §6, stated verbatim); purity and deconvolution are locally computable; HEST-1k's 1,276 profiles are aligned to WSIs at <1.15 µm/px and openly downloadable; and — from the *T2* assassin, unprompted — *"the pairing T2 needs does not exist; the pairing T1 needs already exists, aligned"* with an instrument that is *"demonstrably sensitive (spot Pearson 0.230→0.431)."*

---

## 2. IS IT NATURE-TIER? — scored against `nature_bar.md`, honestly

The bar: **a capability that did not previously exist or a finding that overturns a prevailing assumption, proven on independently-collected data from institutions the model never saw. Pay in ≥2 of four currencies; one must be external generalization.**

| Currency | CALIBRA pays | Evidence it counts |
|---|---|---|
| **1. Novelty** (≥1) | **Two items.** (a) *A capability that did not exist*: recovery-calibrated certification of a cross-modal biological claim — the first cross-modal audit that reports what effect size it would have missed. (b) *A finding that changes what the field believes*: encoder rank on H&E→molecular benchmarks is uninformative about biological capacity, because the metric is dominated by an identity term. | This is precisely the **D1 slot** — Nat Biomed Eng 2026, PMID 41034516, which trained **no new model** and cleared the bar on three transferable findings ("diversity > volume"; "FMs are complementary, so ranking them is the wrong frame"). Also the **C2/C3 slot** — Ahlmann-Eltze *Nat Methods* 2025 and Kernfeld *Genome Biol* 2025, both pure "you did not run the baseline" corrections. |
| **2. External generalization** (**MANDATORY**, observed floor ≥5 independent cohorts) | **Comfortably cleared, and this is CALIBRA's strongest suit.** TCGA (32 cancers, 241 TSS sites) · CPTAC (14 collections, ~2,400 WSI subjects, CC BY 4.0) · HANCOCK (763 patients, Erlangen, 14-yr OS, CC BY 4.0) · SurGen (843 CRC, NHS Scotland) · GTEx (970 donors, 29 tissues) · CAMELYON17 (5 Dutch centres, CC0) · HEST-1k (**180 source cohorts, 26 organs**) · STimage-1K4M (1,149 slides). | Floor is 5 (ENLIGHT-DeepPT). CHIEF used 32 slide sets / 24 hospitals; Hetairos 11 centres. HEST alone brings 180 cohorts and 26 organs. **Multi-continent, multi-scanner, multi-platform, and — uniquely — multi-*target-modality*.** |
| **3. Breadth-or-depth** (≥1) | **Breadth.** ≥19 frozen encoders × 3 target modalities (bulk RNA, spot ST, CPTAC proteomics) × ≥6 cohorts × the N/I/C/B decomposition, plus survival. | Reference frame is CONCH 14 / MUSK 23 / Prov-GigaPath 26 / UNI 34 / GPFM 72 tasks. D1's own frame is 19 models — `nature_bar.md` states explicitly: *"19 models is now the reference frame."* We meet it. |
| **4. Credibility** (≥1, increasingly ≥2) | **Three.** (a) benchmarked against the actual field standard (the live HEST leaderboard and the 19-FM Nat BME frame), not legacy CNNs; (b) **biological mechanism recovered and interpreted** — what the identity term *is* (site, scanner, block provenance, composition) and what survives it; (c) an open released harness others can run on their own claims. Plus DECAT run as an external post-hoc diagnostic. | HEX and ENLIGHT-DeepPT both close on mechanism. `nature_bar.md`: *"The interpretation step costs no GPU and is what converts a Methods-tier result into a Nature Medicine / Nature Cancer-tier one."* |

**Currencies paid: 3–4 of 4, with the mandatory one paid heavily.** That is above the floor of every accepted paper in the corpus.

### Scored against `discovery_standard.md` — the no-wet-lab checklist

Mandatory **G0–G4**: G0 typed (§1.2, no model noun) · G1 the identity term **is** the confound-only null, printed beside ours · G2 zero-parameter naive predictor + ridge + ResNet50 + EfficientNet-B0 + composition, all numbers printed with CIs · G3 ≥19 structurally different encoders, and **the invariance is the finding, correctly re-typed as a data-level claim** as G3 requires · G4 ≥3 seeds, ≥2 preprocessing pipelines, hyperparameter sweep, effect-size distribution reported not best run. **All five cleared, and G1/G2/G3 are already substantially in hand from the existing +0.07 work — which `discovery_standard.md` §6 calls out as "rare."**

Independent-evidence gates, need ≥2, **we carry five**: **G5** CPTAC/HANCOCK/SurGen/GTEx/CAMELYON17/HEST · **G6** CPTAC proteomics+phospho (the cheapest available orthogonal modality, inventoried and unwired — `discovery_standard.md` names wiring it as the #1 credibility purchase per unit effort) · **G7** held-out-confound recovery *plus* registered above-floor programs — both halves, as required · **G8** Replogle GWPS directional concordance with a random-gene-set null, stating the Squires-yes / Varici-no identifiability boundary explicitly · **G9** HANCOCK 14-yr OS, SurGen n=426, CPTAC-3 640 deaths, covariate-only C-index printed alongside.

Tier-raisers, **all three**: **G10** the named-artifact section with a stated artifact fraction *is the paper* · **G11** analysis frozen and registered before HANCOCK/SurGen are downloaded · **G12** the harness is a public standing bet.

`discovery_standard.md`'s own rubric: *"G0–G4 + G8 + two others + G10 → Mechanistic discovery without wet lab. Nature / Science tier."* **CALIBRA lands on that top row.** I am reporting that faithfully, and I am also telling you it is the rubric's ceiling, not a prediction.

### The honest ceiling, stated without hedging

**Guaranteed core (Phases 1–2, cannot fail to produce a result): *Nature Methods* or *Nature Biomedical Engineering*.** The precedents are exact: D1 (Nat BME 2026, a third-party benchmark of 19 FMs that trained nothing), C2 (Nat Methods 2025, "does not yet outperform simple linear baselines"), D2 (Nat Med 2025 editorial, "A benchmarking crisis in biomedical machine learning" — the venue announcing appetite for this genre).

**Upside branch (Phase 3 lands): *Nature Cancer* or *Nature Medicine*.** Requires above-floor programs with a CPTAC protein shadow, Perturb-seq directional concordance, and a ΔC-index increment in HANCOCK/SurGen. That is the HEX shape (Nat Med 2026, trained on **382 samples**, validated on **2,298 patients across 6 cohorts**, closed on mechanism).

**Not reachable: *Nature* main.** Nature-main in this space in 2024–26 went to Prov-GigaPath (1.3B tiles), CHIEF (60,530 WSIs / 24 hospitals), MUSK (50M images + 1B tokens), PathChat. Those are scale-currency or new-interaction-modality papers. We have neither and cannot buy either on one A100 with open data — and `nature_bar.md` §Part 3 is explicit that the scale slot is both unreachable *and* now devalued. **Do not aim at it. Aiming at it is how you get a Nature Communications paper instead of a Nature Methods paper.**

---

## 3. THE METHOD

Frozen encoders throughout. **Train nothing large.** The A100 does inference-once-and-cache plus small heads — the exact profile `nature_bar.md` §Part 3 identifies as the reachable pattern.

### Component 1 — The instrument (NEW; this is the methods contribution)

For any tuple (encoder *E*, target *Y*, confound set *S*, estimator *f*, cohort):

1. **Cross-fitted residualisation.** Regress out *S* from both modalities, fold-safe — nuisance models never see the evaluation fold. *S* is restricted to **full-coverage covariates only**: cancer type, TSS site (or site-cluster where df is prohibitive), tumour purity/ESTIMATE, cell-type deconvolution proportions. Grade/stage/RIN enter only as a sensitivity analysis on the subsets where they exist (2,077 cases / 8 cancers; 1,628 cases respectively), reported as coverage-limited.
2. **Spectrum + permutation null.** Cross-fitted canonical correlation spectrum {ρ₁…ρ_k}; null by permuting patients within cancer-type × site strata. Report the **whole spectrum**, never the top component.
3. **Spike-recovery calibration — the new object.** Draw a random unit direction *u* in image-embedding space and *v* in molecular space; inject `Y' = Y + α·(Xu)vᵀ` with α set so the true induced canonical correlation is `r_true`. Sweep `r_true ∈ {0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40}`, ≥100 random (u,v) draws per level. Run the **entire** pipeline — residualisation included — on spiked data. Output:
   - the **recovery curve** r̂ vs r_true,
   - the **empirical detection floor** (smallest r_true recovered in ≥80% of draws),
   - the **attenuation slope**,
   - and the observed real-data r̂ plotted **on that curve**.
   Additionally spike along *structured* directions (a Hallmark program's loading vector) as well as random ones, since a random direction is a favourable case.
4. **Held-out-confound positive control.** Rotate one morphologically-legible covariate out of *S* (MSI status, TP53 mutation, PAM50/CMS consensus subtype) and require recovery at its independently-known strength. This is G7a executed as an instrument check rather than as a finding.
5. **Attenuation anchor.** Estimate the image↔RNA tissue-mismatch constant two ways: the **TS-vs-DX slide-type contrast within TCGA** (18,425 vs 11,901 slides; a GDC metadata join) and the **same-sample GTEx** cohort. Report the observed channel both raw and attenuation-corrected, with the correction's own uncertainty.

**Why this is the contribution.** Spike-recovery calibration is routine in genomics (ERCC spike-ins) and in Bayesian methodology (simulation-based calibration). It has **never** been applied to certify a cross-modal morphology–molecular claim. It is the difference between "we adjusted and found nothing" and "we adjusted, our floor is 0.031, and we found 0.068." One of those is a null. The other is a measurement.

### Component 2 — The identity decomposition (the finding)

For every (encoder × target × cohort), partition the reported metric into four nested terms, each with a bootstrap CI:

- **N — naive term.** Zero free parameters. Cancer-type mean for bulk; **per-slide mean for spots**. This single baseline is missing from the entire HEST leaderboard and from the 46+ models the 2026 *Brief Bioinform* survey catalogues.
- **I — identity term.** What an oracle on {TSS site, scanner/vendor, slide ID, source cohort} alone achieves beyond N. Includes the CAMELYON17 measurement of *I* in isolation (5 centres, no molecular target at all — the cleanest possible read on the pure identity channel).
- **C — composition term.** Capacity-matched: CellViT nuclei counts/types (ships with HEST-1k) lifted by random features to dimension d, versus the FM embedding PCA-truncated to the same d. Report **B(d) as a curve**, not a point. Deconvolution-derived composition reported separately and explicitly labelled an upper bound (it is a linear functional of the target).
- **B — residual biological capacity.** Everything else, read against the Component-1 floor.

**The claim is N + I + C ≫ B, and B is flat across encoders.** Encoders: ≥19, spanning H-Optimus-0/1, UNI/UNI2-h, CONCH, Virchow/Virchow2, Prov-GigaPath, TITAN, GPFM, CTransPath, plus ResNet50 (2015) and EfficientNet-B0 (5.3M params) as the deliberately unfair floor.

### Component 3 — The capacity spectrum (the upside; this is what buys Nature Cancer)

With the instrument calibrated, ask *what is inside B*:

- Score ~10,000 MSigDB / Reactome / Hallmark programs, consensus subtypes, and de-novo cross-fitted components. **Plot every program's calibrated effect size against the detection floor on the same axis.** Programs above floor are the channel's contents; programs below floor are honestly reported as indeterminate rather than absent.
- **Irreducibility screen** at |r| < **0.5**, not 0.8 — `kill_feas_T3` §7.2 is right that 0.8 is an accept-the-null threshold with no power analysis, and a component at r=0.79 with a proliferation signature would be falsely called novel. Fix this now regardless of what else happens.
- **Protein shadow** (G6): do above-floor programs have CPTAC proteomic/phospho correlates? Morphology→RNA→protein agreement is far harder to explain as artefact than morphology→RNA alone.
- **Spatial localisation** (HEST): where in the tissue does an above-floor program live? Converts a number into a picture a pathologist can adjudicate. This is *localisation*, not predictability — the only defensible use of ST here, exactly as `thesis_T3` §2.1 concluded.
- **Perturbational concordance** (G8): for any program claimed mechanistic, a pre-specified directional prediction against Replogle GWPS with a random-gene-set null. State the Squires ≈1-intervention/node identifiability regime as satisfied and the Varici interaction regime as **not** satisfied (no doubles) — reviewers reward that and punish its absence.
- **Prognostic increment** (G9): ΔC-index over stage+grade+subtype+known signatures in HANCOCK / SurGen / CPTAC-3, bootstrap CI, REMARK-compliant (JNCI 2018, doi:10.1093/jnci/djy088).

### What is borrowed vs new

**Borrowed:** frozen public FM encoders; cross-fitted CCA; permutation nulls; ESTIMATE / xCell / quanTIseq / CIBERSORTx; CellViT nuclei (ships with HEST-1k); Fisher-z power arithmetic; Celligner-style alignment logic; DECAT as a post-hoc baseline.

**New:** (1) the spike-recovery calibration curve and empirical detection floor for a cross-modal channel; (2) the held-out-confound positive control that survives its own adjustment; (3) the N/I/C/B decomposition with capacity-matched composition baselines; (4) the TS-vs-DX attenuation constant; (5) replication of the *structure* (not a number) across two target modalities and ≥6 cohorts; (6) the released harness.

### How it uses MORPHEUS's existing assets

- **6,192 paired TCGA WSI+RNA, H-Optimus features extracted and uncapped** → the Phase-1 discovery cohort. Zero acquisition cost. (Resolve the 6,192 vs 6,443 discrepancy first — `kill_feas_T2` §0 recomputed 6,443 from disk; both are internally consistent, likely a QC filter. A pre-registered paper cannot ship with two cohort sizes.)
- **The leakage-controlled held-out-cancer split** → the transfer arm. **Fix the count before anything is registered:** the thesis text says 14 dev / 21 test (=35, against 32 cancers), `build_paired_split.py` and `preflight.py` default to **11/22 (=33)**, and `refit_mlp_clip.py` defaults to **11/21 (=32)**. Three different numbers, one of them internal to the repo. This is a one-hour fix and a fatal embarrassment if it reaches a reviewer.
- **Established finding (i)** — +0.07 method-invariant, 46–49% cohort structure → the anchor observation the whole paper explains. G1/G2/G3 already substantially satisfied.
- **Established finding (ii)** — effective-rank collapse ~40 vs ~180 sibling, covariance term recovers +53 (2.1×, 3 seeds) with **no change in benchmark score** → **this is where finding (ii) finally gets typed correctly.** Currently it fails G0 ("a representation has higher rank" is not a proposition about biology). Inside CALIBRA it becomes evidence: *if the recovered +53 dimensions are the ones carrying above-floor programs in the capacity spectrum, then rank collapse with no benchmark movement is a direct demonstration that the benchmark is blind to the dimensions that carry biology.* Two independent results converging on "the benchmark is blind" is a stronger thesis than either alone — `discovery_standard.md` §6.3 says exactly this. **This is the internal coherence that makes CALIBRA a MORPHEUS paper rather than a generic audit.**
- **Replogle GWPS (11,258 × 8,248), CCLE, GDSC, DepMap scripts** → G8, used correctly (directional program-level concordance), not the way T2 used them.
- **CPTAC proteomics inventory** → G6, wired at last.

---

## 4. VALIDATION PLAN

Meeting the no-wet-lab discovery checklist. Every cohort below is open-access; every one has a verified route in the scout files.

| # | Cohort | URL / accession | Licence & size | What it validates | Falsification criterion |
|---|---|---|---|---|---|
| 1 | **TCGA** (on disk) | — | — | Discovery: recovery curve, N/I/C/B, capacity spectrum. n=6,192; r_min=0.036 in the full-coverage stratum. | Detection floor > 0.05 after reducing *S* to full-coverage covariates → the instrument cannot resolve the channel; report the floor-vs-covariate tradeoff curve as the result and stop. |
| 2 | **CPTAC** — 14 TCIA collections | `cancerimagingarchive.net/browse-collections` ; pathology browser `pathdb.cancerimagingarchive.net/eaglescope/dist/` ; Aspera `faspex.cancerimagingarchive.net` ; RNA via GDC `api.gdc.cancer.gov` ; proteomics via PDC `pdc.cancer.gov` | **CC BY 4.0**; ~2,400 WSI subjects; LUAD 431.5 GB, CCRCC 190 GB, UCEC 154 GB, GBM 151.5 GB → budget ~1.5–2 TB | **G5** external replication of the *decomposition ratio* + **G6** orthogonal proteomics/phospho shadow. | The identity fraction differs from TCGA by more than its bootstrap CI → the phenomenon is TCGA-specific, and the paper narrows to "TCGA benchmarks are broken." Still publishable, lower tier. |
| 3 | **HANCOCK** | `doi.org/10.7937/rcty-5h16` → `cancerimagingarchive.net/collection/HANCOCK` (Aspera) | **CC BY 4.0**; 763 HNSC patients; 4.5 TB; H&E + LN slides + IHC TMA + text; **OS to 14 years** | **G5 + G9 + G11.** Registered analysis frozen *before* download. Best long-follow-up outcome cohort in open data. | ΔC-index CI over stage+subtype includes 0 **and** the identity-only model matches the full model → the clinical arm reports null, which is the audit's point. |
| 4 | **SurGen** | BioImage Archive **S-BIAD1285**, `doi.org/10.6019/S-BIAD1285` (Globus) | 1,020 WSIs / 843 CRC; survival on 426; KRAS/NRAS/BRAF/MMR. **Licence field COULD-NOT-VERIFY — check the landing page first.** | **G5 + G9** second outcome cohort, second nation, molecular labels for the held-out-confound control. | Same as #3. |
| 5 | **GTEx v8/v10** | `gtexportal.org/home/histologyPage` + `gtexportal.org/home/downloads/adult-gtex` | Open; 25,306 samples / 970 donors / 29 tissues; **image and RNA from the same sample** | **The attenuation anchor** (the only same-sample cohort that exists) **and a specificity negative control**: tumour-intrinsic components must *not* replicate in post-mortem normal tissue; composition components should. | Loadings replicate *fully* in GTEx → the component is a tissue-composition/pre-analytic axis, not tumour biology. **This is failure mode #1 and it is the modal outcome — design for it.** |
| 6 | **HEST-1k + HEST-Bench** | `huggingface.co/datasets/MahmoodLab/hest` (2.01 TB, subset by oncotree) and `.../hest-bench` (42.2 GB); code `github.com/mahmoodlab/HEST` | **CC BY-NC-SA 4.0**, gated (HF account + accept terms; `raw/main/README.md` returns 401). Need `cellvit_seg` + full `st` matrices from the **parent**, not just bench. | **Second target modality.** 1,276 profiles, 180 cohorts, 26 organs. Per-slide-mean baseline, composition-only delta, random-gene-panel null. | Per-task bootstrap CIs exceed the ~0.008 inter-encoder spread → "the encoders tie" is unprovable from HEST alone; pool with #7. **Check this first (§5, Risk 4).** |
| 7 | **STimage-1K4M** | `huggingface.co/datasets/jiawennnn/STimage-1K4M` ; `github.com/JiawenChenn/STimage-1K4M` | **MIT** (audit per-source GEO terms), ungated; 1,149 slides / **4,293,195 pairs**; full transcriptome per spot | **Power backup for #6** and the escape from the top-50-HVG artefact. HF viewer is broken (17,138 inconsistent gene columns) — `snapshot_download` and write your own loader. | If HEST and STimage disagree on the decomposition → platform-specific artefact; report both, claim neither. |
| 8 | **CAMELYON17** | `s3://camelyon-dataset` (`aws s3 ls --no-sign-request`) | **CC0**; 1,000 WSIs / 200 patients / **5 Dutch centres** | **The pure identity term in isolation** — 5 centres, no molecular target. Measures *I* with no biology to confuse it. | Centre is *not* decodable above chance → contradicts Howard 2021 and undermines the identity thesis. Very unlikely. |
| 9 | **DepMap 26Q1 + Sanger Project Score** | figshare `ndownloader.figshare.com/files/62677015` (431 MB, **CC BY 4.0**); Zenodo `10.5281/zenodo.20355477`; `cog.sanger.ac.uk/cmp/download/essentiality_matrices.zip` (241 MB) | Open, scriptable (**the DepMap portal itself is behind Cloudflare Turnstile — no API**) | **G8 pipeline-independence.** Any concordance that does not survive both Broad *and* Sanger pipelines is a Broad artefact. | Concordance present in DepMap, absent in Project Score → drop the G8 leg, keep G5/G6/G7/G9. |
| 10 | **Replogle GWPS** (on disk) | — | 11,258 × 8,248 pseudobulk | **G8** pre-specified directional program concordance + random-gene-set null. | Above-floor programs show no directional concordance beyond the random-gene null → report as such; the audit does not depend on it. |
| 11 | **CMap LINCS 2020** | `s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/` | Open. **clue.io retired 2026-01-31; the LINCS-2020 rebuild has no GEO accession. MIRROR IT NOW — this is a genuine data-loss risk.** | Secondary G8 (genetic↔chemical concordance in the 27 `trt_xpr` lines). | — |
| 12 | **MSigDB** (Hallmark, Reactome, C2, C5) | `gsea-msigdb.org/gsea/msigdb` | Open | The >10,000-set irreducibility screen at |r| < 0.5. | Leading component correlates ≥0.5 with any catalogued set → restatement, not discovery. Report it and move on. |
| 13 | *(optional, clinical arm)* HER2-TUMOR-ROIS `10.7937/E65C-AM96` (n=85, CC BY 4.0); Ovarian Bevacizumab Response `10.7937/TCIA.985G-EY35` (n=78); Post-NAT-BRCA (n=64); IMPRESS (n≈126, **email the authors — no Data Availability statement exists**) | TCIA | CC BY 4.0 | The complete open universe of WSI + treatment + response: **~350 patients across four unrelated contexts.** Enumerating it is itself unpublished (fallback F4). | Any ΔAUC claim below 0.10 here is undetectable. **Do not make one.** |

**Not to be used, and why (state this in the paper so reviewers do not ask):** AACR GENIE has **no images at all** — PRISSMM "Pathology/Imaging" are abstractions of *reports*. No open immunotherapy cohort with slides exists at scale: TCGA has ipilimumab 20 / pembrolizumab 4 / nivolumab 1 / durvalumab 1. CPTAC is dead for response (PR=9, PD=6, SD=33). MOSAIC and POSEIDON are not open. Hartwig requires an application.

---

## 5. REMAINING KILL RISKS AND THE CHEAPEST GO/NO-GO FOR EACH

Ordered by *(probability × damage) / cost to resolve*. **Everything in the top four is resolved in week 1 with zero downloads.**

| # | Risk | Earliest, cheapest resolver | Cost | Decision rule |
|---|---|---|---|---|
| **1** | **The residualisation destroys the signal** — the recovery curve shows the detection floor exceeds any plausible effect, and the instrument is as uninterpretable as T3's was. *This is the risk that decides the programme.* | **Run the spike-recovery calibration on the on-disk 6,192 with *S* = {cancer type, TSS site, purity, deconvolution}.** No downloads. | **2–3 days.** | Floor ≤ 0.03 → proceed. 0.03–0.05 → shrink *S* (site→site-cluster; drop deconvolution) and report the floor-vs-covariate tradeoff curve, which is itself a publishable object. > 0.05 at every reduction → **stop and pivot to fallback F3 (the attenuation constant)**, which does not depend on residualisation at all. |
| **2** | **DECAT (arXiv:2605.31504) already did it.** Same cohort family, same modality pair, null-referenced, 8,979 TCGA patients. | **Read DECAT in full.** Then run it as a post-hoc baseline on MORPHEUS embeddings — it is model-agnostic, so it is hours. | **1 day + hours.** | If DECAT reports a *recovery* curve (not just a null reference), the methods novelty is gone → pivot the contribution to multi-cohort × multi-target-modality replication + the capacity spectrum, and cite DECAT as the instrument. If not — the likely case — cite it, differentiate in the abstract, pass it, and use passing it as evidence. |
| **3** | **Scooped on the control-adjusted ST leaderboard.** Live threats: HistoPrism (ICLR 2026, arXiv:2601.21560, explicitly claims to fix "prior variance-based assessment limitations"), CPNN (arXiv:2603.18461, learns cell-type compositional weights from images), sCellST (Nat Commun 2026, doi:10.1038/s41467-025-67965-1, nuclei model *matching* patch-FM at spot level), HESCAPE, CHRep. | **Read all five in full.** Also recheck bioRxiv `10.1101/2023.09.20.558624` (403 this session — the one unresolved residual risk from the prior-art sweep). | **2 days.** | If HistoPrism reports the mean/composition decomposition → demote HEST Tier 1 to a confirmation and lead with bulk + CPTAC. If sCellST's result already is the composition delta → the composition baseline must be counts/types **only**, no learned morphology embedding, and the deliverable is the FM-minus-composition delta on HEST-Bench specifically. |
| **4** | **HEST per-task n = 2–24 makes "the encoders tie" statistically unprovable.** Seven of nine tasks have 2–4 patients; PRAD is 23 samples from 2 patients. | **Bootstrap the per-task CI width on HEST-Bench.** (`thesis_T1` F4, still unresolved.) Note the leaderboard publishes **point estimates with no standard deviations**, so the "tie within 0.008" premise cannot be assessed from published numbers at all. | **1 afternoon** (after H-Optimus feature extraction — `hest-bench/fm_v1` ships **only ctranspath**, so budget single-digit A100-hours first). | CI width > 0.008 → pool tasks and add STimage-1K4M + SpaRED for power, and never claim a per-task tie. |
| **5** | **Slide provenance is uncontrolled** (18,425 Tissue vs 11,901 Diagnostic; no `slide_type` field anywhere in the repo). | **GDC metadata join onto the existing store.** No re-extraction. | **1 day.** | If the store is DX-dominated, the TS-vs-DX contrast is thin → fall back to GTEx same-sample for the attenuation anchor. Either way this becomes fallback paper F3. |
| **6** | **Disk.** Verified again this session: `C: 927G total, 879G used, 48G free, 95%`, on a OneDrive-synced path. HEST-Bench alone is 88% of remaining free space. | **Buy ≥8 TB external, non-OneDrive.** | **1 day, ~$150.** | Non-negotiable, blocks everything in §7. Do it in week 1. |
| **7** | **B > 0 and encoder-dependent** — the thesis is wrong. | Falls out of Phase 1 automatically. | Free. | **Not a kill.** The paper inverts to a positive result: "here is the encoder property that carries biological capacity." Better paper, harder to write. This is the high floor working as designed. |
| **8** | **The leading component is tumour purity / composition in new coordinates.** The modal outcome, per Fu 2020 ("reflect tumor composition") and Jones/Engelhardt (13,360 GTEx samples → cell-type heterogeneity, ischemic time, demographics, mechanical ventilation). | Regress the component on ESTIMATE/xCell/CIBERSORTx; and the GTEx specificity check (#5 in §4). | Days, Phase 1–2. | R² > 0.6 → the component is composition. **Report it, with the artifact fraction, as G10.** Barrio-Hernandez conceded 4% in their abstract and that is *why* the rest was believed. |
| **9** | **Venue lands below Nature Methods.** | Nothing resolves this early; it is decided by whether Component 3 produces above-floor programs with a protein shadow. | — | Accept it. The floor is a strong specialist paper; there is no version of this programme where four months of work produces nothing publishable. |

---

## 6. THREE-PHASE EXECUTION PLAN

### Phase 0 — Week 1. Reads, fixes, and the disk. **No downloads, no GPU.**

- Buy storage (Risk 6). Move the working tree off the OneDrive-synced path.
- Read DECAT, HistoPrism, CPNN, sCellST, HESCAPE in full (Risks 2, 3).
- **Fix the split count.** Thesis says 14/21 (=35); `build_paired_split.py`/`preflight.py` default to 11/22 (=33); `refit_mlp_clip.py` defaults to 11/21 (=32). Pick one, make the repo internally consistent, register it.
- **Resolve 6,192 vs 6,443.** One cohort size, documented QC filter.
- GDC slide-type join → TS/DX label on every slide in the store (Risk 5).
- Run the TCGA RECIST×WSI×RNA join — the 20-minute query `thesis_T5` flags as "the single most decision-relevant number in this document." It decides whether the clinical arm exists at all.
- Register the analysis plan and the fixed adjustment protocol. **Before any HEST or HANCOCK byte is fetched.** This is G11 and it is free.

### Phase 1 — Weeks 2–6. **FIRST RUNNABLE MILESTONE.**

> **Deliverable: the recovery curve and the N/I/C/B decomposition on the 6,192-patient on-disk TCGA cohort, for ≥19 frozen encoders, with a permutation null and an empirical detection floor.**

Everything is on disk. Compute is hours-to-days on one A100 (`kill_feas_T3` §6: *"compute is a rounding error"*). This produces the paper's Figure 1 and Figure 2 and **cannot fail to produce a number** — the only question is which number.

Milestone acceptance: a plot of r̂ vs r_true with an 80%-recovery floor marked, the observed channel placed on it, and N/I/C/B with bootstrap CIs for every encoder.

### Phase 2 — Months 2–4. Second target modality, first external cohorts.

- HEST-Bench + HEST-1k parent subset (`cellvit_seg`, full `st` matrices) → per-slide-mean baseline, capacity-matched composition delta, random-gene-panel null. Pool STimage-1K4M if Risk 4 bites.
- CPTAC WSI + GDC RNA + PDC proteomics → G5 replication of the decomposition ratio, G6 protein shadow.
- CAMELYON17 → *I* measured in isolation across 5 centres.
- GTEx → attenuation anchor + specificity negative control.

Deliverable: the structure replicates across target modality and cohort — or it does not, which is equally a finding and equally reported.

### Phase 3 — Months 4–7. The biology leg and the registered outcome cohorts.

- Capacity spectrum with the floor drawn on the same axis; irreducibility screen at |r| < 0.5 against >10,000 sets.
- Replogle directional concordance with a random-gene-set null; Sanger Project Score as the pipeline-independence check.
- **Download HANCOCK and SurGen only now**, with the analysis already frozen. ΔC-index over stage+subtype+known signatures, REMARK-compliant, covariate-only C-index printed alongside.
- Write the named-artifact section with a **stated artifact fraction** (G10).

### "Meaningful outcomes along the way" — the fallback papers, ranked

Every one of these ships from infrastructure the primary thesis needs anyway. **There is no branch of this programme that produces nothing.**

| Rank | Fallback | Ships after | Honest venue |
|---|---|---|---|
| **F1** | **The calibrated cross-modal audit harness + the identity decomposition on TCGA and CPTAC.** The instrument alone, with two cohorts. Guaranteed. | Phase 1 + early Phase 2 | *Nature Methods* Brief Communication / *Nature Communications* |
| **F2** | **The control-adjusted HEST leaderboard + the composition-only delta.** T1 Tier 1+2, correctly scoped as a component rather than a thesis. Both assassins agreed this is worth doing and not worth being a thesis. | Phase 2 | NeurIPS D&B / *Nature Methods* Matters Arising / *Brief Bioinform* |
| **F3** | **The image↔RNA tissue-mismatch attenuation constant.** TS-vs-DX within TCGA; same-sample GTEx as the anchor. Unpublished, needed by everyone in this field, and derivable from a metadata join. **Independent of Risk 1 — this is the escape hatch if the residualisation kill lands.** | Phase 0–1 | *Bioinformatics* / *Genome Biology* / *GigaScience* |
| **F4** | **The open-data-universe audit for clinical actionability.** `thesis_T5` §3 enumerated it with live GDC/TCIA counts and nobody has published it: ~1,400–1,500 TCGA chemo-era patients + ~350 external across four unrelated contexts, plus the power floor the field sits ~10× below. | Phase 0 (the join) | Correspondence / data resource / *Nature Medicine* comment |
| **F5** | **Finding (ii) re-typed.** Rank collapse ~40 vs ~180, +53 recovery (2.1×, 3 seeds), **with the recovered dimensions shown to carry above-floor programs.** This is the only route by which finding (ii) becomes biology rather than trivia. | Phase 3 | *Nature Methods* if the biology link lands; ICML/NeurIPS otherwise |

---

## 7. WHAT THE USER MUST GO DOWNLOAD

**0. FIRST, BEFORE ANYTHING ELSE — buy ≥8 TB external storage, not on a OneDrive-synced path.** Verified this session: `C: 927 GB total, 879 GB used, 48 GB free, 95% full`. HEST-Bench alone (42.2 GB) eats 88% of what is left. Nothing below can start until this is done. 4 TB is the aggressive-subsetting minimum; 8 TB is comfortable.

**Tier 1 — needed for Phase 1–2 (get these first):**

1. **Pathology FM weights** (HuggingFace; most are "gated: auto" = instant on form submission) — `bioptimus/H-optimus-0` and `H-optimus-1` (Apache-2.0), `MahmoodLab/UNI` + `UNI2-h`, `MahmoodLab/CONCH`, `paige-ai/Virchow` + `Virchow2`, `prov-gigapath/prov-gigapath`, `MahmoodLab/TITAN`, GPFM, CTransPath, plus torchvision ResNet50 and EfficientNet-B0. **Target ≥19 to meet the Nat BME reference frame.**
2. **HEST-Bench** — `huggingface.co/datasets/MahmoodLab/hest-bench` — 42.2 GB, gated (HF account + accept terms), CC BY-NC-SA 4.0. *Note: `fm_v1` ships **only ctranspath** — you must extract H-Optimus features yourself.*
3. **HEST-1k parent, oncology subset** — `huggingface.co/datasets/MahmoodLab/hest` — full repo is 2.01 TB; use `snapshot_download(allow_patterns=...)` filtered by oncotree code. **You need `cellvit_seg` and the full `st` matrices, which are in the parent, not in bench.** Budget 500 GB–1 TB. Library: `github.com/mahmoodlab/HEST`.
4. **STimage-1K4M** — `huggingface.co/datasets/jiawennnn/STimage-1K4M` — MIT, ungated, 1,149 slides / 4,293,195 pairs, full transcriptome per spot. **The HF viewer is broken; use `snapshot_download` on the `ST/`, `Visium/`, `meta/` folders and write your own loader.**
5. **CAMELYON17** — `aws s3 sync --no-sign-request s3://camelyon-dataset .` (us-west-2). CC0, 1,000 WSIs, 5 Dutch centres. Instant, free, no registration.
6. **DepMap 26Q1 Chronos gene effect** — `https://ndownloader.figshare.com/files/62677015` (431 MB, CC BY 4.0, DOI 10.6084/m9.figshare.31660582). **Do not scrape the DepMap portal — it is behind Cloudflare Turnstile.** Plus DepMap 25Q3 omics/`Model.csv` from Zenodo `10.5281/zenodo.20355477`.
7. **Sanger Project Score** — `https://cog.sanger.ac.uk/cmp/download/essentiality_matrices.zip` (241 MB, direct, no gate). The independent-pipeline check.
8. **MSigDB** — `gsea-msigdb.org/gsea/msigdb` — Hallmark + Reactome + C2 + C5.
9. **CMap LINCS 2020 — MIRROR THIS WEEK.** `https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/`. clue.io retired 2026-01-31 and **the LINCS-2020 rebuild has no GEO accession**. At minimum: `siginfo_beta.txt` (465 MB), `cellinfo_beta.txt`, `compoundinfo_beta.txt`. Ideally the level-5 GCTX too. This is the one item with a real deadline.

**Tier 2 — needed for Phase 2–3:**

10. **CPTAC pathology, all 14 TCIA collections** — browse `cancerimagingarchive.net/browse-collections` and `pathdb.cancerimagingarchive.net/eaglescope/dist/`; download via **IBM Aspera Faspex** (`faspex.cancerimagingarchive.net/aspera/faspex/public/package`) — use the `ascp` CLI, not the browser plugin. CC BY 4.0, ~2,400 subjects, budget 1.5–2 TB. **The "Limited access" badge on CPTAC-GBM/HNSCC applies to the radiology only; the slides remain CC BY 4.0.**
11. **CPTAC RNA-seq** from GDC (`api.gdc.cancer.gov`, CPTAC-3 = 1,683 cases incl. 640 deaths; CPTAC-2 = 342) and **proteomics/phospho** from PDC (`pdc.cancer.gov`). **Note: GDC returns ZERO slide images for CPTAC-3 — imaging is TCIA-only, so this is a three-portal join across three ID conventions. Compute the true triple-intersection n; do not quote 2,400 as the paired n.**
12. **GTEx v8/v10** — expression matrices from `gtexportal.org/home/downloads/adult-gtex`; histology images from `gtexportal.org/home/histologyPage`. Fully open.

**Tier 3 — download only after the analysis is registered (G11):**

13. **HANCOCK** — `doi.org/10.7937/rcty-5h16` → `cancerimagingarchive.net/collection/HANCOCK`, Aspera, CC BY 4.0, **4.5 TB**. (The project site `hancock.research.fau.edu` did not resolve; use the TCIA mirror.)
14. **SurGen** — BioImage Archive **S-BIAD1285**, `doi.org/10.6019/S-BIAD1285`, over **Globus**. **Check the licence field on the landing page before use — it was COULD-NOT-VERIFY.**

**Tier 4 — optional clinical arm:**

15. TCIA **HER2-TUMOR-ROIS** (`10.7937/E65C-AM96`, n=85, 40 GB, CC BY 4.0) · **Ovarian Bevacizumab Response** (`10.7937/TCIA.985G-EY35`, n=78, 253.8 GB) · **Post-NAT-BRCA** (n=64) · **IMPRESS** (n≈126 — **no Data Availability statement exists; email Kun Huang / Zaibo Li**).

**Already on disk, no action:** TCGA paired WSI+RNA (6,192/6,443), Replogle GWPS + K562-essential + RPE1, CCLE (436 MB), GDSC (131 MB).

**Do not bother with:** AACR GENIE (no images at all), MOSAIC, POSEIDON, ORIEN/AVATAR (all closed), any open ICI cohort with slides (does not exist — TCGA has pembrolizumab n=4). PAIP is unreachable (`wisepaip.org` refused connection). HTAN needs a browser session or the Synapse API, not WebFetch — worth one re-scout, not worth planning around.

---

## 8. THE ONE-PARAGRAPH VERSION, FOR THE ABSTRACT

*Deep-learning models predict molecular state from routine H&E at correlations that appear substantial and that rank hundreds of pathology foundation models against one another. We show that this ranking is uninformative. Using a spike-recovery calibration procedure that reports, for the first time, the effect size a cross-modal confound adjustment would have missed, we decompose measured H&E→molecular performance into a zero-parameter naive term, a non-biological identity term, a cell-composition term, and a residual biological capacity. Across ≥19 frozen encoders, three target modalities and ≥6 independently collected cohorts spanning four continents and two target technologies, the identity term exceeds the residual biological term, and the residual is invariant to encoder choice. We report the channel's calibrated capacity, the programs that lie above the detection floor, their proteomic shadow, and their prognostic increment — and we release the harness, so that any cross-modal biological claim can state what it would have missed.*

---

## 9. CITATION HYGIENE INHERITED INTO THIS DOCUMENT

Every number in this proposal is carried from a file in this directory whose author verified it against a live API this session, or was re-verified by me. **Nothing here was written from memory.** Items still flagged **COULD-NOT-VERIFY** upstream and therefore not to be asserted in a manuscript: the CPTAC WSI∩RNA∩proteomics triple-intersection n; the SurGen licence; the IMPRESS access route and cohort composition; the HEST per-technology split; PANDA's exact counts and licence; HTAN's entire contents; TANGLE (Jaume et al., CVPR 2024); Coladan-human3K's release channel; any citation count (Semantic Scholar returned HTTP 429 across every session).

**Three fabricated citations have already contaminated this project, and the T2 sweep caught a fourth in formation** — "Coladan, Genome Medicine 2026" is a **method name, not an author**; the real citation is Wang Z, Yang C, Tang X, Yin E, Yao Y, Luo Y, He J, Sun N, *Genome Medicine* 2026, PMID 42449400, DOI 10.1186/s13073-026-01713-y. Fix it wherever it appears. Verify every reference against a live API before it enters a manuscript.
