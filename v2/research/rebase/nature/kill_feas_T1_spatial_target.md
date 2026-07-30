# KILL REPORT — T1_spatial_target (inverted: "target-invariance of the +0.07–0.10 ceiling")

**Role: feasibility assassin. Constraints assumed: no wet lab, open-access data only, one A100-40GB.**
**Verdict: KILLED.**
Date: 2026-07-29. All external facts below were re-verified this session via WebFetch (arXiv API,
PubMed E-utilities, HuggingFace API, raw GitHub). WebSearch was unavailable. Anything unverified is
marked COULD-NOT-VERIFY and is not used as evidence.

---

## 0. What I am killing

Not the strawman. I am killing the **inverted** thesis the scout recommends
(`thesis_T1_spatial_target.md` §4.2), i.e.:

> The ~+0.07–0.10 method-invariant ceiling is TARGET-INVARIANT — it holds identically for bulk RNA
> over 6,192 TCGA patients and for spot-level ST over 1,276 HEST samples; the apparent r≈0.42 is
> per-slide mean + cell-type composition and collapses to the same +0.07–0.10 after adjustment.

The scout rates this **Feasibility 5/5** and calls it "very likely the most executable thesis in the
set." **That rating is wrong, and it is wrong for a specific reason: it is computed against the wrong
bottleneck.** Compute is not the bottleneck. *Measurability of the load-bearing quantity* is, and it
is not fixable with money, GPU-hours, or an external drive.

I concede up front, because an honest kill has to:

- **Compute is fine.** H-Optimus-0 is a ViT-G/14 (~1.1B params). HEST-Bench is 259 samples; at Visium
  densities that is order 10⁶ spot patches. Forward-only inference at 224² on an A100-40GB is
  single-digit GPU-hours. Nothing in Tiers 1–2 is compute-bound. **Do not kill this on compute.**
- **Access is fine.** HEST is gated-but-free (HF account + accept terms), CC BY-NC-SA 4.0 — acceptable
  for academic publication. `bioptimus/H-optimus-0` is `"gated":"auto"` (instant approval on form
  submission) under **Apache-2.0**. Verified via `huggingface.co/api/models/bioptimus/H-optimus-0`.
  **Do not kill this on access.**
- **The cell-line→tumor gap does not apply.** T1 uses no CCLE/GDSC/Perturb-seq. That standard attack
  is inapplicable here and I am not going to manufacture it.

The kill is on three other grounds, in descending order of finality.

---

## 1. FATAL — "cell-type composition" is not measurable under the constraints, so Tier 2 is uninterpretable in *both* directions

Tier 2 is the load-bearing deliverable. The scout's own words: *"the damning number… if that delta is
≈0, the field's premier molecular benchmark is a nuclei-counting benchmark."* The entire mechanism of
the inverted thesis — "the r≈0.42 is composition" — rests on being able to **hold composition fixed
and look at the residual**. Under no-wet-lab + open-access, there are exactly three ways to
operationalise "composition." **All three fail, and they fail in different directions, which is why
the experiment cannot produce evidence either way.**

### (a) CellViT nuclei counts/types — what the thesis actually proposes → *capacity-confounded, non-identifiable*

`cellvit_seg` does exist in the HEST parent repo (verified: top-level dirs are `cellvit_seg`,
`metadata`, `patches`, `patches_vis`, `pixel_size_vis`, `spatial_plots`, `st`, `thumbnails`,
`tissue_seg`, `transcripts`, `wsis`, `xenium_seg`). So the data is there. The problem is what it *is*:
a PanNuke-class nucleus detector's **prediction of composition from the very same H&E patch**, yielding
a ~5–15-dimensional count vector.

Now run the comparison the thesis specifies: 1536-d H-Optimus embedding vs ~5-d nuclei-count vector,
both predicting the same 50 genes. **The FM wins essentially by construction** — on representational
capacity, before any biology enters. And then the result is non-identifiable between:

- H1 (thesis is wrong): the FM encodes molecular state beyond composition; or
- H2 (thesis is right but untestable this way): CellViT is a lossy composition estimator and the FM's
  advantage is a *better estimate of the same composition*.

There is no way to separate H1 from H2 with image-derived composition, because both hypotheses predict
the identical observable. Worse: this outcome — FM ≫ nuclei baseline — is the scout's own
**pre-registered falsifier F2**. So the modal outcome of the load-bearing experiment is *falsification
of the thesis for a reason that has nothing to do with biology*. An experiment whose most likely result
is a methodological artifact that happens to trip your own falsifier is not a 5/5-feasible experiment.

### (b) Deconvolution of spot expression (RCTD / cell2location / CARD) → *circular*

The obvious upgrade — get composition from the transcriptome instead of the image — is worse.
Deconvolved composition is a **linear functional of the target**. Regressing the target on a linear
functional of itself removes target variance by construction, so the post-adjustment residual is small
for *every* predictor and the FM-vs-baseline delta compresses toward zero **whether or not the
mechanism is real**. That is not a test; it is an identity. Any referee sees this in one reading, and
it is precisely the direction that would manufacture the thesis's desired "+0.07–0.10 collapse."

### (c) Independent ground-truth composition → *does not exist where the number lives*

The only non-circular, non-capacity-confounded route is a per-cell ground truth on the same tissue:
mIF/CODEX/IMC (wet lab — excluded) or single-cell-resolution ST. HEST does contain Xenium (verified:
"27 new Xenium" Sep-2024 plus "18 new Xenium including Xenium 5k" Feb-2026; `xenium_seg` and
`transcripts` dirs exist). But this route closes on itself:

1. **Xenium panels are 300–500 genes selected to discriminate cell types.** On the one platform where
   composition is measurable, the target is composition *by panel design*. Finding "expression is
   explained by composition" there is a property of Xenium's panel-selection procedure, not a finding
   about H&E morphology. A referee will say exactly this.
2. **The headline number does not live there.** The 0.4229 leaderboard is computed on the Visium/legacy
   per-organ tasks (IDC, PRAD, PAAD, SKCM, COAD, READ, CCRCC, LUNG, LYMPH_IDC, HCC — verified: 10 task
   directories in `MahmoodLab/hest-bench`, not 9). **Where the number is, composition is unmeasurable;
   where composition is measurable, the number does not exist.** There is no overlap to run the
   experiment on.

### The mandate question, answered directly

> *"Does the required data actually exist at adequate scale with the needed pairing?"*

**No.** The pairing this thesis requires is **H&E + spot expression + composition ground truth
independent of both, on the platform carrying the headline result**. That triple does not exist in
open data and cannot be created without wet lab. Every substitute is either circular (b) or
capacity-confounded (a), and the one honest measurement (c) sits on a platform whose target is
composition by construction.

### Novelty-adjacent corroboration (flagging; not my lane)

**CPNN — arXiv:2603.18461v1**, "Cell-Type Prototype-Informed Neural Network for Gene Expression
Estimation from Pathology Images," Nishimura, Bise, Matsuo, Hirose, Kojima, 19 Mar 2026. VERIFIED
(arXiv API, full abstract retrieved). It "learns cell-type compositional weights directly from images
and models the relationship between prototypes and observed **bulk or spatial** expression," evaluated
on **three slide-level and three patch-level ST datasets**, reporting the **highest Spearman across all
settings** plus interpretable per-cell-type weights. Two consequences for T1: the composition mechanism
is already built and published, and a composition-structured model is already known to be *competitive
or better* — which drains Tier 2 of its news value even in the unlikely branch where it works.

---

## 2. FATAL — "target-invariance" is not a measurable quantity; it is a curve-fit with the answer pre-specified

The headline claim asserts that two numbers are *the same*. They are not the same kind of number.

| | Bulk (MORPHEUS) | Spot (HEST-Bench) |
|---|---|---|
| Sampling unit of the correlation | **patients**, within cancer | **spots**, pooled within held-out slides |
| Gene universe | genome-scale (BulkFormer) | **top-50 HVGs, re-selected per organ task** |
| Gene selection | none | data-driven HVG on the target |
| Independent replicates for inference | ~190/cancer (6,192 / 32) | **2–24 slides per task** |
| Cohort / platform | TCGA FFPE, one pipeline | 180 cohorts, legacy-ST / Visium / VisiumHD / Xenium |

Nothing constrains these two statistics to coincide. Their variance denominators are different objects.
So "**it holds identically**" is not a prediction that can be confirmed — it is a coincidence that would
have to be *produced*, and the thesis hands the analyst **four independent knobs** to produce it with
the target value (+0.07) known in advance:

1. which mean to subtract (per-slide / per-cohort / per-gene / per-organ);
2. what counts as a "spatially-patterned" gene (Moran's I? SPARK? threshold?);
3. how the random-gene-panel null is constructed (matched on expression level? on variance? on dropout?);
4. how composition is adjusted out (§1, itself undefined).

Each of these moves the resulting delta by **more than the ±0.015 band the claim asserts**. This is a
textbook garden of forking paths, and it is the first thing a statistically literate referee at any
serious venue will say. The thesis is not falsifiable as written — which is the mirror image of a
feasibility failure: you cannot build an experiment whose outcome is informative.

### The one principled fix is underpowered *and* measures a different thing

The right way to make the comparison commensurable is to **pseudo-bulk each HEST slide** and predict
that from the WSI — a genuine matched bulk-vs-spot contrast on identical tissue. Two problems:

- **Power.** 1,276 samples / 26 organs ≈ 49 per organ, heavily skewed, and 2–24 per HEST-Bench task,
  against ~190/cancer in TCGA. At n = 20, the Fisher-z standard error of a Pearson r is
  1/√(n−3) ≈ **0.24** — roughly an order of magnitude wider than the 0.03 band ("+0.07 vs +0.10") the
  claim is trying to resolve. Pooling the 10 tasks to buy power destroys the within-organ conditioning
  that the whole bulk comparison depends on. This is the scout's own F4 gate, unresolved, and it does
  not pass on the pseudo-bulk route.
- **Construct validity.** Visium pseudo-bulk is not bulk RNA-seq: partial tissue capture (one section,
  spot grid, non-exhaustive), 3′ bias, capture-efficiency and per-slide batch. Even a powered version
  would compare bulk-RNA prediction to a *different measurement* and re-import the exact per-slide
  confound the thesis is trying to subtract out.

**No matched large bulk+ST cohort exists to escape this.** Verified: PubMed query for
("spatial transcriptomics" AND "bulk RNA" AND matched/same-patients AND cohort) returned **5 records
total** (PMIDs 42426811, 42126225, 41876493, 41738622, 39058036) — single-disease atlases and biomarker
studies (cervical SCC, MIBC, stage-I LUAD, HCC, Merkel cell), none a hundreds-of-patients matched
bulk+ST resource. Nothing at the scale required. This confirms the scout's limitation #2 and elevates
it from a caveat to a structural blocker for the headline claim.

### And the published tie cannot be assessed at all from the published numbers

Verified on `raw.githubusercontent.com/mahmoodlab/HEST/main/README.md`: the 25-model leaderboard reports
**point estimates only — no standard deviations**. So "25 FMs tie within 0.008," the fact that makes the
method-invariance parallel look compelling, **cannot be established from the leaderboard**. It must be
re-derived from scratch, and the re-derivation is exactly the underpowered inference above.

---

## 3. FATAL (against the stated goal) — success does not clear the bar, and the tier that would is contradicted by evidence already in hand

User constraint: Nature-tier (new architecture / SOTA / new evaluation paradigm / real biological
discovery). Score Tiers 1–3 against it:

- **Tier 1 + Tier 2 fully successful** = a control-adjusted leaderboard audit. Scout's own ceiling: 2/5,
  NeurIPS D&B / Nature Methods correspondence / Nat Comms at the absolute top. **Does not meet the
  stated goal.**
- **Tier 3** — isolate the composition-independent gene/pathway subset and show it transfers to held-out
  cancers and to bulk TCGA — is the only Nature-tier route, and it fails three ways at once:
  1. It **inherits §1**: "composition-independent" is undefinable without a composition measurement.
     The dependent variable of Tier 3 cannot be constructed.
  2. It is **directly contradicted by evidence the project already holds**: *"Benchmarking the
     translational potential of spatial gene expression prediction from histology,"* Wang C., Chan A.S.,
     Fu X. et al., **Nat Commun 2025, PMID 39934114, DOI 10.1038/s41467-025-56618-y** — 11 image→ST
     methods, 5 ST datasets, **external TCGA validation**, finding within-image spot accuracy is a poor
     predictor of cross-study transfer. That is Tier 3's exact transfer experiment, already run,
     already negative.
  3. The scout itself prices it "unlikely."

**A ladder whose achievable rungs miss the goal and whose goal-meeting rung is undefinable and
pre-refuted is not a research programme.**

### What a referee will demand that this project cannot produce

1. An **independent composition ground truth** on the Visium tasks carrying the headline number (§1 —
   impossible without wet lab).
2. Proof the conclusion is not an artifact of the **top-50-HVG-per-task target**. HEST-Bench ships
   top-50 HVG `.h5ad` only; establishing that the effect holds genome-wide requires the full `st`
   matrices from the 2 TB parent (§4).
3. A head-to-head against **HistoPrism (arXiv:2601.21560, ICLR 2026)** and any HEST v2 — both of which
   the scout flags as live F3 falsifiers and neither of which has been read.

---

## 4. Supporting resource findings — the 5/5 feasibility rating rests on two statements that are false

These are schedule risks, not kills, but they matter because they are the literal justification for the
5/5 rating and both are wrong.

**(i) "using the H-Optimus features already on disk" — FALSE, twice over.**
- MORPHEUS's store is **hard-scoped to TCGA**. `src/data/hoptimus_patch_store.py` line 1:
  `"""Strict, versioned H-Optimus-0 patch-store primitives for TCGA-UT."""`, and its ID normaliser
  rejects anything else:
  ```python
  def normalize_tcga_patient_id(value: str) -> str:
      parts = str(value).strip().upper().split("-")
      if len(parts) < 3 or parts[0] != "TCGA" or len(parts[1]) != 2 or len(parts[2]) != 4:
          raise ValueError(f"Not a TCGA patient/slide identifier: {value!r}")
  ```
  Required metadata columns are `patient_id, slide_id, cancer_type, internal_split, external_split` —
  a TCGA schema. HEST samples cannot enter this store without a rewrite.
- **HEST-Bench does not ship H-Optimus features either.** Verified via
  `huggingface.co/api/datasets/MahmoodLab/hest-bench/tree/main/fm_v1`: `fm_v1` contains
  **`ctranspath` and nothing else**. The leaderboard's H-Optimus numbers are not reproducible from
  distributed features; you must extract them yourself for every bench patch.

  Consequence: not fatal (single-digit A100-hours), but the sentence carrying the 5/5 rating is wrong,
  and "weeks" should be read as "weeks *after* a feature-extraction and schema-porting task that was
  scoped at zero."

**(ii) "the entry point is 42 GB, not 2 TB" — FALSE for this specific protocol.**
- HEST-Bench (42.2 GB) ships **top-50 HVG per task**. But adjustment knob #2 —
  *"restrict to spatially-patterned genes"* — requires computing spatial autocorrelation over the
  **full** gene set to select on. Selecting spatially-patterned genes from within 50 already-HVG-selected
  genes is a near-null operation and will not survive review. That control needs the full `st` matrices.
- `cellvit_seg`, the entire Tier 2 dependency, is in the **2 TB `MahmoodLab/hest` parent**, not in
  `hest-bench`.
- **Hardware, measured this session: `C:` is 927 GB, 877.9 GB used, 48.2 GB free — 95% full**, and the
  working tree is a **OneDrive-synced path**. HEST-Bench alone (42.2 GB) consumes 88% of remaining free
  space, before extracted features (~10⁶ patches × 1536 × fp32 ≈ 6 GB per encoder). Tier 2 cannot begin.
  Solvable with an external drive — flagged as an operational blocker, not the kill.

**(iii) License.** All HEST assets are **CC BY-NC-SA 4.0**. Share-alike propagates to any derived
control-adjusted leaderboard or dataset you publish. Fine academically; note it.

---

## 5. Falsifier for *this kill* (what would revive T1)

I would withdraw the kill if **all** of the following became true:

1. A public cohort appears with **H&E + spot-level ST + independent per-cell composition ground truth**
   (mIF/CODEX/IMC or untargeted single-cell-resolution ST) on the **Visium-class tasks** where the
   r≈0.42 number lives, at n sufficient for slide-clustered inference — dissolving §1; **and**
2. A matched **bulk-RNA + ST same-patient cohort** at n ≳ 200/indication appears, making
   target-invariance a directly measurable contrast rather than a cross-cohort coincidence —
   dissolving §2; **and**
3. The claim is restated as a **pre-registered directional prediction with a single fixed adjustment
   protocol** ("the adjusted spot delta is < X"), not as a numerical coincidence with four free knobs.

(1) and (2) require wet lab or a resource that does not exist. (3) is free and should be adopted
regardless — but on its own it converts the thesis into a modest, correctly-scoped benchmark note,
which is §3's verdict anyway.

---

## 6. Bottom line

**KILLED — but salvage the afternoon, not the thesis.**

The scout's operational recommendation is right for the wrong reason. It says "run the inverted version
as the cheap control experiment for whichever thesis you pick." That is correct, and this report does
not contradict it: **computing a per-slide-mean baseline and a random-gene-panel null on HEST-Bench is
a legitimate, cheap, one-to-two-week robustness check** that hardens the existing bulk +0.07 finding by
showing the confound structure recurs on a second modality. Do that. It costs an afternoon of thinking
and a few A100-hours (after a feature-extraction task that was scoped at zero, and after freeing disk).

What must **not** happen is promoting it to a thesis. As a thesis it requires measuring a quantity
("composition, independent of the image and of the target") that no open-access asset provides and no
amount of compute manufactures; it asserts an equality between two statistics with no common
denominator, no matched cohort to compare them on, and four analyst knobs pointed at a pre-known
answer; and its only Nature-tier rung is both undefinable and already reported negative in
Nat Commun 2025. Feasibility is not 5/5. On the deliverable that carries the claim, it is 1/5.

---

### Verification ledger (this session)

| Fact | Source | Status |
|---|---|---|
| HEST leaderboard = point estimates, **no std devs**; 2 TB; CC BY-NC-SA; ≥45 Xenium samples added 2024+2026 | `raw.githubusercontent.com/mahmoodlab/HEST/main/README.md` | VERIFIED |
| HEST parent repo dirs incl. `cellvit_seg`, `st`, `xenium_seg`, `transcripts` | `huggingface.co/api/datasets/MahmoodLab/hest/tree/main` | VERIFIED |
| `hest-bench` = **10** task dirs + `fm_v1`; `fm_v1` contains **only `ctranspath`** | `huggingface.co/api/datasets/MahmoodLab/hest-bench/tree/main[/fm_v1]` | VERIFIED |
| H-optimus-0: `"gated":"auto"`, Apache-2.0, affiliation form | `huggingface.co/api/models/bioptimus/H-optimus-0` | VERIFIED |
| HEST-1k: 1,229 profiles, 153 cohorts, 26 organs, 25 cancer types, 2.1M pairs, 76M nuclei; NeurIPS'24 | arXiv API; S2 API (164 citations) | VERIFIED |
| CPNN arXiv:2603.18461v1, 19 Mar 2026, Nishimura et al., full abstract | arXiv API | VERIFIED |
| STFlow ICML 2025, ">18% relative improvement"; SPADE (Med Image Anal 2025); HistoGPA 2026 | arXiv API; S2 API | VERIFIED |
| No large matched bulk+ST same-patient cohort: 5 PubMed hits, all single-disease atlases | E-utilities esearch + esummary | VERIFIED |
| MORPHEUS H-Optimus store is TCGA-only | `src/data/hoptimus_patch_store.py` (read) | VERIFIED |
| `C:` 927 GB, 48.2 GB free, 95% used, OneDrive-synced | `df -h` + `Get-PSDrive` | VERIFIED |
| Nat Commun 2025 translational benchmark PMID 39934114 | carried from scout report; **not independently re-verified this session** | INHERITED |
| HistoPrism arXiv:2601.21560 ICLR 2026 | carried from scout report; **not independently re-verified this session** | INHERITED |
| Semantic Scholar composition/matched-cohort queries | HTTP 429 on both | COULD-NOT-VERIFY |
| TANGLE (Jaume et al., CVPR 2024) | still unverified | COULD-NOT-VERIFY |
