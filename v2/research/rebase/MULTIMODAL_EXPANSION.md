# Multiscale Causal Observability System — the multimodal expansion

*Roadmap document. Nothing here is claimable until the Tier-1 experiments in
`HANDOFF_EXPERIMENTS_NOW.md` return. Build order is deliberate and staged; do not jump ahead.*

---

## 0. The organising question

Not *"can we fuse more modalities into a better embedding?"* — that space is crowded (OmiCLIP, STORM,
HEX, ROSIE, TITAN, PathChat, HONeYBEE). Instead:

> **Which intervention-derived biological effects are observable at each biological scale, through
> which clinical modality, and with what level of evidence?**

The output is a **multimodal causal observability atlas**, not a multimodal embedding.

## 1. THE PREREQUISITE RULE (answers "do we certify before we build promptable?")

> **You cannot prompt what you cannot certify. Certification is not a delay before the interesting
> work — it is the thing that makes the interesting work trustworthy.**

A natural-language interface that reports "this tumour shows immune exclusion" is only as sound as the
certificate behind that axis. Without one you have a confident interface over an uncertified
representation, which is **worse than no interface** — it launders uncertainty into fluent prose.

**This rule applies at every stage, including the H&E + RNA/ST system being built now.** For each
modality `m`, the gate before exposing anything through the copilot is:

1. the observability operator `Π_m` is estimated **on a discovery fold only**;
2. its axes clear the **CALIBRA detection floor**;
3. they pass the **confound certificate** (must fail to predict site/scanner/batch);
4. they **replicate** in untouched patients and ≥1 external cohort;
5. failures are recorded and exposed alongside successes.

Only certified axes become promptable. Uncertified ones are visible in the atlas **marked as
uncertified** — never returned as an answer.

*This is where the original "promptable multimodal representation" vision lands: it was right, but it
needed a foundation. The atlas is that foundation; the copilot is a query layer over it.*

## 2. Build order (do not reorder)

| Stage | Modalities | Gate to pass before proceeding |
|---|---|---|
| **S1 (now)** | H&E + bulk RNA | Tier-1 experiments E0–E3; certified legible subspace exists |
| **S2** | H&E + spatial transcriptomics | Spatial replication of S1 axes; cell-of-origin resolved |
| **S3** | + spatial proteomics (mIF/CODEX/IMC) | RNA→protein propagation certified |
| **S4** | + DNA/CNV/methylation **as context**, not targets | Context-transport model validated |
| **S5** | + radiology | Macro-scale observability certified |
| **S6** | copilot / promptable layer | **All of §1 satisfied for every exposed axis** |

## 3. Modality-specific observability operators

Each modality observes a *different consequence* of perturbation; they are not interchangeable channels.

From the intervention dictionary `D = [δ₁ … δ_G]` and patient code `y_i ≈ D a_i`, learn a separately
certified operator per modality:

```
Π_m : C → L_m          (C = intervention-response space, L_m = recoverable subspace of modality m)
```

Giving every gene a **modality observability fingerprint**:

```
M_g = [ ‖Π_H&E δ_g‖ , ‖Π_sRNA δ_g‖ , ‖Π_prot δ_g‖ , ‖Π_radio δ_g‖ ]
```

The question upgrades from *"is this gene H&E-visible?"* to *"at what scale does this perturbation
become visible, where does it manifest, and which modalities independently confirm it?"*

| Modality | Scientific role |
|---|---|
| H&E / WSI | tissue architecture, cytology, necrosis, differentiation, immune morphology |
| Spatial transcriptomics | **where** the perturbational programme is active |
| Spatial proteomics | whether RNA-level effects **propagate to protein**, cell states, neighbourhoods |
| Bulk RNA / proteomics | higher-powered patient-level estimates |
| DNA / SNV / CNV / methylation | **context** — upstream constraints and natural perturbations |
| Radiology | whole-lesion scale, invasion, multifocality, longitudinal change |
| Clinical / path text | human interpretation, specimen context |
| Perturbation assays | the causal dictionary itself |

## 4. Shared, private and synergistic effects

- **Shared:** `L_shared = L_H&E ∩ L_prot ∩ L_sRNA`. An axis independently recoverable from three
  modalities has a far stronger mechanistic interpretation than one inferred from a slide alone.
- **Modality-private:** transcriptionally strong but morphologically invisible; protein-visible but
  RNA-weak; radiologically visible but not resolvable on one biopsy. **These are not failures — they
  define the observational limits of each assay**, and they are publishable results.
- **Synergy:** `Δ_syn(k; m,n) = R_{m,n}(k) − max(R_m(k), R_n(k))`.
  **Certify synergy only when the joint model beats the strongest unimodal model by a pre-registered
  margin, in untouched patients and an external cohort.** This is the guard against the standard
  multimodal-paper failure: a 1% gain reported as "discovering cross-modal biology."

## 5. Spatial transcriptomics — four roles, not just validation

1. **Localise causal coordinates.** Per spot, solve `a_s = argmin ‖y_s − D a‖² + λΩ(a)` → a spatial map
   for every certified axis. (Not per-gene spot prediction — that is now commodity.)
2. **Test morphological colocalisation.** H&E-predicted axis map vs directly measured axis map: spatial
   correlation, region overlap, boundary agreement, neighbourhood agreement, patient-clustered
   uncertainty, cross-platform replication.
3. **Resolve cell-of-origin — critical.** A bulk axis may arise from tumour-cell state change, immune
   infiltration, stromal composition, necrosis/purity, or multi-population interaction. Label each axis
   as *tumour-autonomous / immune-mediated / stromal / composition-driven / neighbourhood-dependent*.
   **Without this, "legible axis" may just mean "cell composition" — see the baseline in §9.**
4. **Detect spatially conditional effects.** Does the same perturbation have different visible
   consequences at invasive margins vs hypoxic cores vs immune-rich vs necrotic regions, or
   primary vs metastatic?

## 6. Spatial proteomics — prioritise over another RNA modality

RNA measures *intended* programmes; protein and neighbourhood show whether they became **functional
tissue phenotypes**. Per axis: does the transcriptional effect propagate to protein? which cell types?
does it alter immune–tumour contact, exclusion, compartmentalisation? is it colocalised with the
H&E-visible regions? *Do not compete on virtual staining (HEX, ROSIE, VirTues exist) — use measured or
predicted protein maps to **certify downstream consequences** of causal axes.*

## 7. DNA / CNV as context, not targets

A mutation is **not** a controlled perturbation: it may be subclonal, co-occurring with other drivers,
lineage-dependent, LoF vs GoF, and CNV affects neighbouring genes. Use DNA/CNV/methylation to define the
context-dependent operator `δ_g(c_i)` where `c_i` = lineage, co-mutations, copy state, epigenetic state,
composition, treatment. Then ask: **in which genomic contexts is the reference effect transportable to
this patient?** Catalogue fields: reference effect, context-adjusted effect, transportability
confidence, natural-genetic support, **contradictory contexts**.

## 8. The promptable layer — a scientific query engine, not a chatbot

PathChat and TITAN already do image+language. The differentiator:

> **Every natural-language answer compiles into deterministic atlas operations and returns with
> certificates, provenance and visible evidence. The LLM never decides whether a gene is causal — it
> translates a question into validated tool calls, and never emits a number it did not retrieve.**

Tools: `search_axes` · `score_case` · `show_axis_regions` · `compare_modalities` · `inspect_gene` ·
`inspect_certificate` · `test_confound` · `compare_cohorts` · `find_counterexamples` ·
`simulate_perturbation` · `export_report`.

Every response carries: slide overlays, spatial maps, gene loadings, effect sizes + CIs, cohort
replication, confound results, **negative evidence**, provenance, model/certificate versions.

## 9. BASELINES — how we know we are not fooling ourselves

**This section is mandatory. Each baseline targets a specific way the result could be fake.**

### Baselines the method must BEAT
| Baseline | What it kills if it wins |
|---|---|
| **Curated pathway supervision** (Hallmark/Reactome; i.e. MoPE's answer) | the entire premise — if curated does as well, interventional coordinates buy nothing |
| **Random dictionary** (size/spectrum-matched random directions replacing `δ_g`) | the "causal" content — if random works equally, the dictionary is just a rotation |
| **PCA / NMF basis of tumour expression** | the need for perturbation data at all |
| **Text-prior (GenePT-style LLM gene embeddings)** | devastating if it wins — means we recovered literature, not measurement |
| **Cell-composition** (deconvolution, CellViT nuclei counts) — **capacity-matched** | *the most likely true explanation*: "legible" may just mean "cell composition is visible" |
| **Zero-parameter naive** (cancer-type mean; **per-slide mean** for spots) | everything — this baseline is missing from the HEST leaderboard and must be reported |
| **Best unimodal model** (for any multimodal claim) | all synergy claims |

### Controls that must FAIL (if they pass, the pipeline is broken)
- **Site/scanner/batch prediction** from certified axes → must fail (the confound certificate).
- **Random gene sets** → must not clear the detection floor.
- **Shuffled gene labels** on the dictionary (same directions, permuted labels) → gene-level
  attribution must collapse while the subspace persists. *This separates "the subspace is real" from
  "we can name which gene it belongs to" — they are different claims.*
- **Modality-shuffled pairing** (modality *m* of patient *i* with patient *j*) → all cross-modal
  agreement must vanish.

### Positive controls that must PASS
- RNA-input states predicting RNA-derived targets (circular by construction) → strong channel.
- A **held-out known-legible covariate** (e.g. MSI, TP53, consensus subtype excluded from the
  adjustment set) → must be recovered at its independently known strength.
- Synthetic spike at a level above the detection floor → recovered.

**Interpretation rule:** a result is reportable only when it beats *every* "must beat" baseline **and**
every "must fail" control actually failed. Report all of them, including the ones we lose.

## 10. Architecture — modular, open, separately certifiable

`causal-dictionary` (build/harmonise `D`, geometry, equivalence groups, effective rank) →
`tumor-coder` (project bulk/sc/spatial into `D`; structured sparse coding, uncertainty) →
`observability/{histology,spatial_rna,spatial_protein,radiology,multimodal}` (learn `Π_m`) →
`certify` (spike recovery, confounds, held-out, external, spatial, stability, falsification) →
`atlas` (Gene · Perturbation · Axis · EquivalenceGroup · Modality · CancerContext · Cohort ·
Certificate · Falsification · Case · Region) → `copilot` (NL → tool calls only) →
`viewer` (browser dashboard + **QuPath extension**, following the WSInfer deployment pattern).

**Standards — adopt, don't invent:** SpatialData + AnnData (spatial), OME-pyramids (WSI/multiplex),
DICOM (radiology), Zarr (arrays), Parquet (tables), versioned JSON (certificates), immutable run
manifests (provenance). Public data via CRDC / HTAN.

## 11. Flagship deliverable — the Gene × Modality Causal Observability Matrix

Rows: genes or **perturbation-equivalence groups** (not genes alone — see the identifiability limit).
Columns: H&E · spatial RNA · spatial protein · radiology · bulk RNA · proteomics · cross-modal
combinations. Each cell: certified effect magnitude, confidence, #cohorts, cancer contexts, evidence
tier, failure status. Click-through to source perturbation evidence, patient-level distributions,
spatial maps, representative slide regions, protein/RNA consequences, confound results, and
**counterexamples**.

## 12. Operating modes
**Atlas** (explore certified biology across cancers/cohorts/modalities) · **Case** (one patient: active
axes, localisation, cross-modal agreement, similar cases, uncertainty, *failed* certificates) ·
**Discovery** (external researchers fit candidate axes, get detection limits + confound checks +
replication, and submit successful axes back).

## 13. Honest risks
- **Cell composition explains everything** (§9). The single most likely deflation. Test it early.
- Each added modality multiplies confounds; the certificate must be re-earned per modality, never inherited.
- Radiology↔pathology pairing is scarce; expect low power, report `n` and abstain rather than certify.
- Context-transport (`δ_g(c)`) is unsolved; without it the atlas describes *reference* effects, not
  patient-specific ones. Say so explicitly in every gene record.
- Scope discipline: **S1 must clear its gates before S2 begins.** The failure mode of this project has
  been breadth substituting for depth.
