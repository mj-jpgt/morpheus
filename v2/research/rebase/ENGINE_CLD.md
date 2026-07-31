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
*legible* subset (those directions that exceed the CALIBRA detection floor and pass the confound
certificate). Three properties hold **by construction**:

1. the basis is **causal**, not curated;
2. it is **high-rank** (~11k directions, so no collapse to a ~50-D manifold);
3. **we never supervise on directions the input cannot express** — the legibility filter enforces it.

**The loss.** Supervision over a large, highly correlated dictionary requires sparse coding, not dense
regression. Proposed objective, per patient:

```
L_PBS = || Πᴠ(y) − D·a ||²  +  λ₁·‖a‖₁  +  λ₂·Σⱼ wⱼ·|aⱼ|
```
where `D` is the (legible) perturbation dictionary, `a` the sparse causal coordinates, `Πᴠ` projection
onto the certified-visible subspace, and `wⱼ` a **legibility weight** — the inverse of direction *j*'s
calibrated detection floor, so the model is penalised for spending capacity on directions the channel
demonstrably cannot carry. The legibility-weighted prior is the novel component; the sparse-coding
backbone is standard and should be cited as such.

**Prediction (falsifiable):** PBS should (a) not collapse, (b) beat curated-target supervision on the
calibrated channel, and (c) leave the confounded benchmark roughly unchanged — the last being the point,
since it is the metric that could not see the difference.

## 5. Discovery experiments (what the engine produces)

**D1 — genome-scale legibility catalogue.** For each of ~20,000 genes × 21 held-out cancers: is the
transcriptional consequence of perturbing it morphologically legible, at what calibrated effect size,
and does it survive the confound certificate. ≈420,000 certified measurements. **Negative entries are
first-class results** — "gene G is *not* legible above floor 0.03" is a claim nobody in this field can
currently make, and it tells you which biomarkers will never be readable from a slide.

**D2 — the composition of the visible subspace (tests H1).** Are the legible directions enriched for
architecture-controlling gene classes (adhesion/ECM/cytoskeleton/proliferation/immune) relative to a
size- and expression-matched null? This is the direct test of H1 and it is a one-figure result.

**D3 — novel morphology–biology links.** Genes that are legible but have **no known morphological
phenotype**. Each is a falsifiable hypothesis — *"perturbing G should change tissue architecture in way
W"* — testable by anyone with a lab, and by us in spatial data.

**D4 — replication.** Spot-level (HEST-1k) and independent cohort (CPTAC) replication of D1/D2. This is
what separates "a property of the channel" from "an artifact of bulk-RNA averaging in TCGA."

## 6. Honest accounting of novelty

| Component | Status |
|---|---|
| Sparse coding / dictionary learning | **Borrowed.** Cite, do not claim. |
| Perturbation atlases as prediction targets | **Crowded.** GEARS/CPA/scGPT-perturb; explicitly not our framing. |
| "Virtual perturbation from histology" | **Taken.** Retired in `NEAR_COLLISIONS.md`; do not reuse the phrase. |
| Collapse onto low-rank targets | **Published** (NRC1). Only the *expressible-intersection* refinement (H2) could be new. |
| Rank ↑ with information flat (O1) | **Ours**, unpublished, and cuts against anti-collapse practice. |
| Perturbation atlas as a **basis** for decomposing an observational cross-modal channel | **No instance found** in our 15-lane corpus. The candidate novel move. |
| Legibility-weighted sparse supervision (PBS) | **Novel as far as we know**; depends on the CALIBRA floor existing. |
| Certified genome-scale legibility catalogue | **Novel**; enabled by CALIBRA, and the reason `kill_feas_T3.md`'s "uncalibratable instrument" objection no longer applies. |

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
