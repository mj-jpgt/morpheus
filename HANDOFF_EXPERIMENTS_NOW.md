# HANDOFF — the four experiments to run NOW

**Read `HANDOFF_BUILD_AGENT.md` first** (orientation, the 14-item mistake list, run commands), then
`v2/research/rebase/ENGINE_CLD.md` (what these experiments are for). Branch: `research/rebase-vision`.

**All four run on data already on disk. No GPU. No downloads.** Each one is decisive and each outcome
is reportable — including the failures. Run them in the order given; E0 can kill the whole engine, so
it goes first.

**Report every result, including negatives, in `v2/research/rebase/nature/EXPERIMENT_LOG.md`.** A
negative here is a finding, not a setback — that is the entire point of running them first.

---

## E0 — Basis transfer: do cell-line perturbation directions align with tumour expression? **[CRUX]**

**Tests H3.** If this fails, the causal-basis engine dies today and we save months.

**Data**
- `PRISM/data/perturbseq/K562_gwps_normalized_bulk_01.h5ad` — 11,258 perturbations × 8,248 genes
  (also `K562_essential_...` 2,285 and `rpe1_...` 2,679 as replication contexts)
- TCGA bulk expression: `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp (1).tsv` (1.88 GB,
  gene-rows × sample-cols, `SYMBOL|ENTREZ` row labels, quoted TSV). `v2/prepare_pancan_rna.py` was
  written for exactly this layout — reuse it, do not write a new parser.

**Method**
1. Build the perturbation matrix `P` (perturbations × genes) as **change vs control** (the `obs` block
   carries `control_expr`; confirm whether the stored matrix is already a delta before subtracting —
   check, do not assume).
2. Intersect gene symbols with TCGA; report the overlap count (expect a few thousand).
3. **Do not test "does P span the space"** — with 11,258 rows in ~8k dimensions it trivially will.
   Test **structural alignment**: compute principal angles between the top-*k* principal subspaces of
   `P` and of centred TCGA expression, for k ∈ {10, 25, 50, 100}.
4. **Null:** random rotations of `P` with matched spectrum (preserve singular values, randomise
   directions), ≥100 draws. Report the observed alignment against that null distribution.
5. Repeat with RPE1 as a second cell context — agreement between two unrelated lineages is evidence
   the alignment is biological rather than line-specific.

**Outcomes**
- **Aligned above null** → H3 holds, the causal dictionary is usable, proceed to E1–E3.
- **At null** → **the engine dies.** Report it: *cell-line perturbation atlases do not transfer to
  tumour expression structure* — a directly useful negative that much current work assumes away.
  Fall back to a data-driven legible basis without causal labels (the discovery catalogue survives;
  the causal interpretation does not).

---

## E1 — Is the added rank empty? (formalise O1)

**Our strongest owned observation, currently a single anecdote.** Make it a measurement.

**Data:** the three diagnostic artifacts in
`discovery_evidence_v2/runs/v21_release_20260720_retry3_resume_safe/artifacts/` plus the CALIBRA
results in `runs/calibra_v2_local/`.

**Method.** For a representation before and after the decorrelation term (rank 49.9 → 103.3), measure
**information**, not rank: (a) held-out canonical correlation with the molecular targets
(`calibra.spectral.heldout_top_cca`); (b) the *number of canonical directions that individually exceed
the CALIBRA detection floor* — the honest count of informative dimensions; (c) per-dimension
information density = (b) / effective rank.

**Prediction:** effective rank rises ~2×, the count of above-floor directions does **not**. If so,
state it plainly: **decorrelation-based anti-collapse regularisation adds capacity without adding
information, and effective rank is not a valid proxy for representational health.** That is a direct,
evidenced challenge to standard practice.

*If the informative count also rises*, O1 was a metric artifact — say so, and the anti-collapse story
needs rewriting.

---

## E2 — Expressible-intersection vs plain NRC1 (tests H2)

**Method (synthetic first, CPU, hours).** Construct targets by projecting real molecular targets onto
subspaces of the **image-expressible** space of known dimension *k* ∈ {5, 10, 20, 50, 100}, holding the
nominal target dimensionality fixed. Train the same small head identically at each *k*, then measure
the trained head's effective rank (`calibra.spectral.effective_rank` — **singular values, see mistake
#1**).

**Three-way falsification**
- final rank tracks **k** → H2 holds; a genuine refinement of NRC1 (collapse to the *expressible*
  intersection, not the target rank). This is the ICLR/NeurIPS-track claim.
- final rank tracks the **nominal target rank** → plain NRC1. **Drop the mechanism claim entirely.**
- tracks neither → collapse is driven by the KL/supcon geometry or weight decay; report that instead.

---

## E3 — The objective ablation, with the anchoring control (tests F2)

As specified in `HANDOFF_BUILD_AGENT.md` §5-A1. Restated here because it gates the PBS architecture.

**Method.** Run the CALIBRA channel measurement on `diagnostic_identity_only_seed42.npz`,
`diagnostic_programme_only_seed42.npz`, `diagnostic_full_seed42.npz`.
**First** read each artifact's `manifest_json` and confirm epochs / LR / seed / token budget match; if
they do not, the ablation is **suggestive, not causal** — report that rather than hiding it.

**The anchoring control is the key question.** `z_identity` is anchored on the frozen MLP-CLIP teacher
with a residual of ≈0, so F2 may partly restate "MLP-CLIP beats our biology head."
`programme_only` has **no anchor** — if its biology channel is still weak, the effect is about the
**objective** (F2 holds, PBS is motivated); if it is strong, F2 was an **anchoring artifact** and PBS
loses its primary evidence. **Escalate immediately in that case.**

---

## Decision table — what to build next

| E0 | E2 | E3 | Next |
|---|---|---|---|
| pass | pass | pass | Build the full CLD engine + PBS. Strongest case. |
| pass | fail | pass | Build CLD + PBS; drop the H2 mechanism claim. |
| pass | any | **fail** | Build CLD/discovery catalogue only; PBS unmotivated — escalate. |
| **fail** | any | any | **Stop. Report the transfer failure.** Fall back to a data-driven legible basis; the discovery catalogue survives without causal labels. |

## Rules
- Emit the repo's unavailability convention everywhere: `metric="status"`, `value=NaN`,
  `note="unavailable_<reason>"`. Never drop a row silently.
- Every headline number needs its permutation null and its held-out counterpart. No absolute CCA alone.
- `rna_*` states are the positive control (RNA→RNA is circular). If they don't show a strong channel,
  the pipeline is broken, not the model.
- Verify the workspace is a real junction, not a stale copy, before trusting any run (mistake #2).
- Do not retrain anything. E0–E3 are all frozen-artifact / on-disk analyses.
