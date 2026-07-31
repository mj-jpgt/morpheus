# Causal Legibility Decomposition (CLD) — engine, hypothesis, architecture, discovery plan

*Status: proposed thesis. The load-bearing hypotheses are NOT yet tested; the experiments that can
kill them are cheap and specified in `HANDOFF_EXPERIMENTS_NOW.md`. Nothing here should be cited as a
result.*

---

## 1. The observation that motivates everything

Two measurements from our own runs, which no existing story explains together:

**(O1) We manufactured capacity with no information in it.** Adding a covariance-decorrelation term
raised the biology head's effective rank **49.9 → 103.3** (2.1×, 3 seeds) and moved the task metric
**0.1366 → 0.1367**. The ~53 added dimensions were *empty*. This cuts against the implicit premise of
the anti-collapse literature (VICReg/Barlow-style decorrelation), which treats rank as a proxy for
representational health. **Rank is manipulable independently of information; it is not the quantity.**

**(O2) The channel is strong in aggregate and weak in every named direction.** Held-out canonical
correlation **0.477** (permutation chance 0.151) against per-target within-cancer, control-adjusted
specificity of **+0.07 — identical across every method including baselines.** If morphology cleanly
encoded a few pathways we would see high per-target scores for those pathways. We do not. The
information is **distributed and misaligned with the curated basis.**

Together: the morphology→molecular channel is not rank-limited, it is *content*-limited, and its
content does not align with the vocabulary the field tests it against (Hallmark/Reactome — built from
bulk transcriptomics and literature, not from tissue architecture).

## 2. The hypothesis

> **H1 (engine hypothesis).** The morphology-visible subspace of tumour expression is spanned by the
> transcriptional consequences of perturbing a small set of **tissue-architecture-controlling** genes
> (adhesion, ECM, cytoskeleton, proliferation, immune recruitment) — **not** by metabolic or signalling
> programmes.

- **If true:** we have stated what histology *is*, as a measurement device, in causal gene terms.
- **If false** (the visible subspace is spanned by metabolic or apparently arbitrary perturbations):
  more surprising, equally reportable — it would mean tissue morphology reflects something other than
  the genes that control tissue architecture.

Two supporting hypotheses, each independently testable:

> **H2 (mechanism of collapse).** A representation supervised onto targets collapses to the rank of the
> **expressible intersection** (target span ∩ what the input modality can carry), *not* to the target's
> own rank. NRC1 (Andriopoulos, NeurIPS 2024) predicts collapse to the **target** rank; our biology head
> sits at ~38 against a Hallmark target of effective rank ~92. **That gap is the anomaly.** If H2 holds
> it is a refinement of published theory; if the head tracks the nominal target rank instead, it is
> plain NRC1 and we drop the claim.

> **H3 (basis transfer).** Perturbation-response structure learned in cell lines (K562, RPE1) aligns
> with expression structure in solid tumours well enough to serve as a basis. **This is the crux risk**
> — it is the same cell-line→tumour gap that killed our earlier virtual-perturbation direction. Here it
> is tested *first* rather than assumed, and a negative result is publishable in its own right: a great
> deal of current work implicitly assumes cell-line perturbation atlases transfer to tumours, and
> nobody has measured it this directly.

## 3. The engine

**Existing framing (what everyone does):** treat perturbation data as a *prediction target* — "predict
the response to perturbing gene G." This is crowded (GEARS, CPA, scGPT-perturb, image-based virtual
perturbation) and, for the image→dependency direction, information-theoretically capped.

**CLD (ours):** treat the 11,258 Perturb-seq perturbations as **11,258 directions in expression space
with known causal meaning — a causal dictionary** — and use it to *decompose an observationally
measured cross-modal channel into causal coordinates*.

```
CALIBRA  ──►  V = the certified morphology-visible subspace of tumour expression
                        │
Perturb-seq ──►  P = causal dictionary (11,258 labelled directions)
                        │
                        ▼
        decompose V in the basis P  ──►  causal coordinates of what morphology sees
```

The output is a statement of the form *"the morphology channel is spanned by the transcriptional
consequences of perturbing genes {G₁…G_k}"* — an interpretable, causally-labelled description of an
otherwise anonymous statistical subspace.

**Why this is not the scooped idea.** We are not predicting perturbation response from an image. We
are choosing a *coordinate system* for a channel measured entirely within TCGA (paired WSI+RNA). The
perturbation data supplies only the basis vectors, so the data-processing argument that capped the
prediction framing does not apply — no signal is passed through a lossy fixed map.

## 4. The architecture: Perturbation-Basis Supervision (PBS)

**What it replaces.** The current biology head is supervised by regression onto ~50 curated Hallmark
scores plus a neighbour-KL and supervised-contrastive term. That is the diagnosed collapse mechanism,
and F2 is the evidence that it is actively harmful (the head trained *for* biology carries less
molecular signal than the head trained for retrieval).

**PBS.** Supervise the biology head on **coordinates in the causal dictionary**, restricted to the
*legible* subspace. Two properties hold by construction: the coordinates are **interventional**, not
curated; and **we never supervise on directions the input cannot express**.

> ### CORRECTION — it is a DICTIONARY, not a basis, and "no collapse" was wrong
> `D = [δ₁ … δ_G] ∈ ℝ^{p×G}` with G ≈ 11,000 in p ≈ 8,000 gene dimensions is **overcomplete and
> highly correlated**. `rank(D) ≤ min(p, G)`, and the *effective* rank is far lower because many
> perturbations converge on shared stress / cell-cycle / apoptosis / interferon / ribosomal responses.
> **11,258 perturbations is not 11,258 recoverable directions**, and the earlier claim that PBS is
> "high-rank so it cannot collapse" was simply false. Claim instead: *a genome-scale, non-curated
> dictionary of ~11,000 measured perturbations rather than dozens of predefined pathways* — and
> **report** algebraic rank, effective rank, dictionary coherence, perturbation-retrieval accuracy,
> and the number of distinguishable perturbation equivalence classes.
>
> **Biological convergence is not model collapse.** If disrupting several different mitotic genes all
> produce the same necrosis/pleomorphism phenotype, that is a fact about the measurement, not a
> failure. The correct object is the **quotient** `C / ker(T)`, where `C` is intervention-response
> space and `T` maps molecular perturbation effects to observable morphology; perturbations `g,h` are
> equivalent when `Tδ_g ≈ Tδ_h` even though `δ_g ≠ δ_h`. We want **no artificial collapse in the source
> dictionary, but permissible biological convergence in the morphologically distinguishable quotient.**

> ### CORRECTION — the dictionary is causal only in its experimental context
> A CRISPRi signature in K562 is an interventional effect **in K562**. Write it honestly as
> `δ_g(c) = E[Y | do(g), c] − E[Y | do(control), c]`, where `c` is cell type, basal state, genotype and
> environment. Without a context-transport model the catalogue measures *the H&E-legible component of a
> reference perturbation signature* — **not** how perturbing that gene would move this patient's tumour.
> Do not write the stronger sentence. The eventual patient-specific form is `S_ig = ‖Π_L(c_i) δ_g(c_i)‖`.

**The loss — legibility as an OPERATOR, not per-gene weights.** A diagonal vector of per-gene weights
wrongly assumes each atom is individually visible or invisible. Morphology may express *combinations*:
neither A nor B individually predictable while a shared A+B programme is highly visible; or two
individually-predictable atoms mutually indistinguishable. So define a PSD **legibility operator**
`W ⪰ 0` on intervention-coordinate space, and use the noise-whitened dictionary Gram
`G_D = Dᵀ Σ⁻¹ D`:

```
L_leg = Σ_i (â_i − a_i)ᵀ W^{1/2} G_D W^{1/2} (â_i − a_i)  +  λ Ω(â_i)
```
with `a_i` the tumour's sparse intervention code, `â_i = f_θ(X_i)` the H&E-predicted code, `W`
retaining only certified recoverable combinations, and `G_D` preventing over-penalisation of
coordinate errors between near-equivalent atoms. **This operator form is the novel component**; the
sparse-coding backbone is standard (cite Webster, CellCap) and must not be claimed.

**Ω must not be plain LASSO.** In a highly correlated dictionary LASSO arbitrarily selects one of
several near-identical atoms, so a gene could be labelled "morphologically consequential" while an
almost-equivalent gene is dropped on sampling noise. Use elastic-net / group sparsity / graph-fused
penalties over similar perturbations, plus **stability selection**, and **report perturbation
equivalence groups whenever individual attribution is not identifiable.**

> ### CRITICAL — avoid circular certification
> If legibility enters the sparse prior, the pipeline becomes self-fulfilling: estimate g is legible →
> bias the decomposition toward g → train the image model to predict g → "discover" that g is legible.
> **This would invalidate the entire catalogue.** Mandatory nested protocol:
> **A** estimate and *freeze* `D` on perturbation experiments containing no pathology cohort →
> **B** encode tumour molecular profiles into `a_i` **without using their images** →
> **C** estimate `W` / the legible subspace on a **discovery fold only**, after confound control →
> **D** train the final H&E model with `D` and `W` **frozen** →
> **E** certify on **untouched** patients + CPTAC + spatial.
> The certificate must rest on data that determined *neither* the dictionary, *nor* the legibility
> operator, *nor* the prior.

**Prediction (falsifiable):** PBS should (a) not collapse, (b) beat curated-target supervision on the
calibrated channel, and (c) leave the confounded benchmark roughly unchanged — the last being the point,
since it is the metric that could not see the difference.

## 5. Discovery experiments (what the engine produces)

**D1 — genome-scale legibility catalogue.** For each gene × cancer, report **two** quantities against
the certified legible projector `Π_L`, in the noise-whitened metric:

```
M_g = ‖Π_L δ_g‖_{Σ⁻¹}                    (absolute legible effect magnitude)
F_g = ‖Π_L δ_g‖²_{Σ⁻¹} / ‖δ_g‖²_{Σ⁻¹}     (fraction of the total effect that is legible)
```
The four quadrants are separately meaningful (large M/large F = strong and mostly visible; large
M/small F = big molecular effect, small visible component; etc.). **Note:** `M_g` needs only `δ_g` and
`Π_L` — **no sparse coding**. Sparse coding is required for *representing tumours and training the H&E
model in those coordinates*, not for ranking genes. Keep those two claims separate.

> **Mandatory nuisance control — otherwise the catalogue is a list of essential genes.** Knockouts
> causing death, cycle arrest, DNA damage, mitochondrial failure or broad stress produce real but
> uninteresting morphology. Report both `M_g^total` and
> `M_g^specific = ‖Π_{L ∩ N^⊥} δ_g‖`, where `N` spans generic nuisance responses (viability,
> proliferation, global stress — definable from DepMap essentiality, which is on disk).
> A catalogue dominated by essential genes is a negative result about the method, not a discovery.

**Negative entries are first-class** — "gene G is not legible above floor 0.03" is a claim nobody can
currently make, and it says which biomarkers will never be readable from a slide.

**Per-gene record must include:** source assay + biological context `c`; `M_g^total`, `M_g^specific`,
`F_g`; the visible modes it loads on and its **equivalence group**; CI and certified lower bound;
guide-efficiency / off-target sensitivity; sensitivity to viability & proliferation; internal, CPTAC
and spatial replication status; and **explicitly listed failed tests.**

**D2 — the composition of the visible subspace (tests H1).** Are the legible directions enriched for
architecture-controlling gene classes (adhesion/ECM/cytoskeleton/proliferation/immune) relative to a
size- and expression-matched null? This is the direct test of H1 and it is a one-figure result.

**D3 — novel morphology–biology links.** Genes that are legible but have **no known morphological
phenotype**. Each is a falsifiable hypothesis — *"perturbing G should change tissue architecture in way
W"* — testable by anyone with a lab, and by us in spatial data.

**D4 — replication.** Spot-level (HEST-1k) and independent cohort (CPTAC) replication of D1/D2. This is
what separates "a property of the channel" from "an artifact of bulk-RNA averaging in TCGA."

## 6. Honest accounting of novelty

*Two claims below were externally challenged and then **verified real by direct API lookup** — treat
them as hard constraints, not opinions.*

| Component | Status |
|---|---|
| Learn coordinates from perturbations, not curated pathways | **Established.** sVAE+ / sparse-mechanism-shift. Cite. |
| Sparse dictionary learning over correlated gene effects | **Established.** **Webster** (graph-regularised, explicitly handles pleiotropy/correlated effects) and **CellCap** (sparse perturbation→programme with context attention). Cite, do not claim. |
| "Don't supervise H&E on molecular variation it cannot observe" | **PRIOR ART — VERIFIED.** **MoPE**, arXiv **2606.02877v2**, Jun 2026, *Pathway-Structured Privileged Distillation for Deployable Computational Pathology*. Explicitly states the partial-observability problem. **Its solution is to retreat to ~50 curated Hallmark pathways.** Our problem statement is therefore *not* new — but MoPE **validates** it, and the solution space is open. |
| Genome-wide gene→morphology catalogue | **PRIOR ART — VERIFIED.** **PERISCOPE**, *A genome-wide atlas of human cell morphology*, **Nature Methods 2025**, >20,000 genes / >30M cells. **But confirmed: cultured cells only — no patient tumour H&E.** So "quantify how much perturbing a gene changes morphology" is *not* novel on its own. |
| Perturbation→tissue morphology in vivo | **Exists** (Perturb-map, CRISPRmap, Perturb-DBiT). **Validation assets, not competitors** — none learns an H&E supervision objective or a certified pan-cancer legibility catalogue. |
| Perturbation atlases as prediction targets | **Crowded.** GEARS/CPA/scGPT-perturb; not our framing. |
| "Virtual perturbation from histology" | **Taken.** Retired in `NEAR_COLLISIONS.md`. |
| Collapse onto low-rank targets | **Published** (NRC1). Only the expressible-intersection refinement (H2) could be new. |
| Rank ↑ with information flat (O1) | **Ours**, unpublished; cuts against anti-collapse practice. |
| **Learning the H&E-legible subspace of an interventional dictionary** | **No precedent found** (ours + external search agree). |
| **Legibility-*operator* interventional sparse supervision** | **Appears genuinely new.** |
| **Certified gene-level projection into patient-tumour-visible causal modes** | **Appears highly novel** — but requires the careful causal wording in §4. |
| "All ~11,000 directions are independent and causal in tumours" | **NOT DEFENSIBLE. Never claim this.** |

**The delta, stated precisely.** MoPE solves partial observability by *retreating to a small curated
vocabulary*. PERISCOPE builds a genome-wide morphology atlas *in cultured cells with fluorescent
phenotypes*. Neither derives coordinates from measured interventions, estimates which combinations are
legible from **routine patient H&E**, restricts supervision to that certified span, and certifies each
mode and gene across patients, an external cohort and spatial tissue. **That composition is the claim.**

**One-sentence contribution:** *We replace curated molecular supervision with a genome-scale
intervention dictionary and learn the certified quotient of that dictionary that routine tumour
morphology can distinguish.*

## 7. Kill conditions (pre-registered)

- **H3 fails** (cell-line perturbation structure does not align with tumour expression structure above
  null) → the causal-basis engine dies. Report the transfer failure; fall back to a data-driven legible
  basis without causal labels (the discovery catalogue survives, the causal interpretation does not).
- **H2 fails** (head tracks nominal target rank) → drop the mechanism claim, keep the engine.
- **H1 fails** (legible directions not enriched for architecture genes) → report it; it is the more
  surprising result and reframes rather than kills the paper.
- **F2 deflates** (the identity>biology gap is explained by MLP-CLIP anchoring) → PBS loses its
  motivating evidence and must be justified on the collapse mechanism alone.

## 8. Venue framing (honest)
The discovery catalogue + certification is a *Nature Methods / Nature BME* shape. H1 landing — a causal
account of what histology measures — is what would lift it to *Nature Cancer / Nature Medicine*. The
general ML contribution (legibility-matched supervision; rank-is-not-information) is the ICLR/NeurIPS
main-track component. None of this is claimable until §7 has been run.
