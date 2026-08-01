# Lane 2 — Programme-target / neighbour-KL redesign to avoid the low-rank manifold

Remit: keep programme/biology supervision WITHOUT pinning the 256-D biology head onto a
~50-D (effectively rank-5-6) target manifold. Ground every claim in the actual code.

## Queries/searches run

- WebSearch: "BLEEP spot expression contrastive histology gene target dimensionality highly variable genes"
- WebSearch: "mclSTExp spatial transcriptomics histology contrastive learning gene expression prediction"
- WebSearch: "HEST-bench gene expression prediction benchmark 50 highly variable genes target dimensionality Pearson"
- WebSearch: "contrastive learning target similarity matrix distillation rank collapse effective rank low-dimensional target manifold"
- WebFetch arxiv 2510.03455 (PEaRL, pathway-enhanced representation learning) — target dimensionality + collapse discussion
- Read: morpheus/v2/losses.py, model.py, training.py, export.py, runner.py (target construction), tests/test_v21_model.py

## Sources

Web:
- BLEEP — Spatially Resolved Gene Expression Prediction via Bi-modal Contrastive Learning, NeurIPS 2023. https://arxiv.org/abs/2306.01859 , https://github.com/bowang-lab/BLEEP . Targets: full expression profile aligned by contrastive learning; evaluated on marker genes, top-50 HEG, top-50 HVG; prediction via kNN imputation in the JOINT EMBEDDING space (not per-gene regression, not a similarity-Gram KL).
- mclSTExp — Multimodal contrastive learning for spatial gene expression, Brief. Bioinform. 2024. https://academic.oup.com/bib/article/25/6/bbae551/7848899 , https://arxiv.org/abs/2407.08216 . Targets: highly variable genes; contrastive maximizes cosine similarity of truly-paired image and gene-expression embeddings (embedding-to-embedding, learns its own rank).
- HEST-Benchmark / HEST-1k, Mahmood lab, NeurIPS 2024. https://github.com/mahmoodlab/HEST . Target: 50 highly variable genes as PER-GENE raw expression; ridge regression; metric = mean per-gene Pearson over top-50 HVG. Direct regression, NOT a similarity-matrix KL.
- PEaRL — Pathway-Enhanced Representation Learning. https://arxiv.org/html/2510.03455 . When it DOES use pathway scores it uses 609-1,100 pathways (Breast 775 / Skin 609 / Lymph 1,100) and 1,000 HVGs for the gene head. Aligns modalities in a 256-D CONTRASTIVE space; treats pathway reduction as feature engineering. Does NOT KL-match a target Gram matrix and does not report/analyze collapse.
- "Breaking the Geometric Bottleneck: Contrastive Expansion in Asymmetric Cross-Modal Distillation." https://arxiv.org/pdf/2603.06698 . Cosine distillation to a target collapses the student to low effective rank (student ~16 vs teacher ~88); adding an InfoNCE contrastive expands the intrinsic manifold ~2.4x.
- WERank (rank-degradation prevention). https://arxiv.org/pdf/2402.09586 . Matrix Information Theory for SSL. https://arxiv.org/pdf/2305.17326 . Corroborate: per-dimension variance floors do not prevent rank collapse; only cross-dimension (covariance/uniformity) terms do.

Code (file:line):
- losses.py:46-53 — `programme_neighbourhood_loss`: KL(row-softmax(z_biology sim / 0.20) || row-softmax(target sim)). This is soft-target distillation of the 256-D biology Gram matrix onto the target Gram matrix.
- losses.py:29-34 — `variance_floor`: per-dimension `relu(1.0 - std)`, does NOT constrain cross-dimension covariance, so does not prevent rank collapse.
- losses.py:37-43 — `whitened_cross_covariance` (separation, weight 0.01): only decorrelates identity vs biology; does not raise biology's own rank.
- losses.py:1, 56-66 — `supervised_programme_contrastive`: positives = programme-similar patients; still uses z_biology-vs-z_biology similarity (same low-rank target graph), weight 0.20.
- model.py:283-290 — `biology_state = self.biology(biology)`; z_biology = normalize(biology_state); programme_mean/log_variance are Linear(256 -> programme_dim=50) downstream of biology_state.
- training.py:194-196 — identity gets symmetric InfoNCE (z_identity WSI vs RNA); biology gets NO paired contrastive.
- training.py:157-160 — neighbourhood term added at weight `weights["neighbourhood"]` (0.20 after warmup, schedule at training.py:29).
- runner.py:33-40 — `residualise_programmes`: cancer-residualized + per-dim standardized 50-D Hallmark, train-fold fit.
- runner.py:172-188 — 50-D programme target, top-8 neighbour graph (`_v2_neighbour_indices`), positive mask.
- runner.py:138-141 — batch exposes `programme_target` (50-D), `programme_positive_mask`, `programme_neighbor_indices`, `programme_target_mask`.
- export.py:54-55 — widths dict: only *_identity/*_biology/full_patient exported; z_context, z_uncertainty, 5x residual NOT exported (dead wiring, confirmed).
- tests/test_v21_model.py:19-34 — reusable `_batch()` fixture (programme_target width 3 in tests).

## Findings

1. Mechanism is a distillation collapse, not a capacity limit. `programme_neighbourhood_loss`
   (losses.py:46-53) forces the 256-D z_biology pairwise-similarity distribution to EQUAL the
   pairwise-similarity distribution of the 50-D cancer-residualized Hallmark target. A Gram
   matrix of 50-D vectors has rank <= 50; cancer-residualized Hallmark scores are strongly
   cross-correlated, so the target's EFFECTIVE rank is ~5-6 — precisely the observed z_biology
   effective rank. KL-matching a normalized-similarity (cosine) structure is textbook
   "cosine distillation onto a low-rank teacher," which the geometric-bottleneck paper shows
   drives student effective rank down to the teacher's. Identity avoids this because it has a
   paired InfoNCE (training.py:194) whose uniformity term expands rank (observed ~84).

2. No WSI->expression method in the literature uses a target-Gram-KL. They use either
   (a) paired image<->expression InfoNCE that learns its OWN high-rank embedding structure
   (BLEEP, mclSTExp, PEaRL), or (b) direct per-target regression against 50-1000 raw values
   (HEST-bench, PEaRL gene head). The MORPHEUS neighbour-KL is the outlier and is the specific
   loss that transfers the target's low rank into the head.

3. Target dimensionality in the field is not the problem; the LOSS FORM is. HEST-bench uses only
   50 targets and does not collapse, because it regresses each gene independently (rank of the
   50-D OUTPUT does not constrain the rank of the input EMBEDDING). MORPHEUS's Gaussian-NLL head
   (losses.py:69-85, model.py:290) is likewise a benign per-target regression — it is NOT the
   collapse driver. The neighbour-KL is.

4. `variance_floor` (losses.py:29-34) and per-dim standardization (runner.py:35) are per-dimension
   and provably cannot prevent rank collapse (WERank / matrix-information-theory sources): a rank-6
   representation can have unit per-dimension std. Only a cross-dimension covariance/uniformity
   term raises rank.

5. supcon (losses.py:56-66) does NOT fix this even though it is "contrastive": its positives are
   the low-rank programme-neighbour graph and its logits are z_biology-vs-z_biology, so it pushes
   toward the SAME 50-D geometry. It is not an independent high-rank signal.

## Recommended change (file:line, exact) — smallest change that fixes the mechanism

Give biology the SAME manifold-expanding objective identity already has (a paired InfoNCE),
against a molecular positive, and STOP forcing 256-D similarity to equal 50-D similarity. This is
exactly the BLEEP/mclSTExp/PEaRL pattern and the geometric-bottleneck fix (InfoNCE expands rank
~2.4x). No new modality, no target-dimensionality change to exports, no architecture change.

Concretely (two edits, one weight retarget):

A. Replace the neighbour-KL with a symmetric InfoNCE between z_biology (WSI) and a projection of
   the raw RNA embedding, using the EXISTING top-8 programme neighbours as extra positives.
   - training.py:157-160: replace the `programme_neighbourhood_loss(...)` call with a call to a
     new `biology_paired_contrastive(z_biology_wsi, z_biology_rna, positive_mask)` that reuses
     `symmetric_infonce` (losses.py:13-18) semantics (paired diagonal + programme-similar
     off-diagonal positives from `programme_positive_mask`, runner.py:141). Keep weight
     `neighbourhood_after_warmup` (training.py:29) but rename its role; 0.20 is a sound start.
   - The RNA view is already computed every step (out_rna, training.py:183) and z_biology exists
     for RNA (model.py:286), so the positive pair is z_biology[WSI] <-> z_biology[RNA] with NO new
     forward pass. This is a true paired contrastive on the biology head — currently absent by
     design (losses.py:1) — and is the single missing high-rank signal.

B. Add ONE cross-dimension anti-collapse term on z_biology so it cannot re-collapse if positives
   are sparse in a batch. Reuse the machinery already in the file: add a whitened-covariance
   OFF-DIAGONAL penalty on z_biology alone (a VICReg-style covariance term). Minimal form:
   in losses.py add `def covariance_offdiag(state)` returning the mean squared off-diagonal of
   the centered covariance of `state`, and add it in training.py next to the existing
   `variance` term (training.py:232-237) at the existing `variance` weight 0.01. This directly
   targets rank (per the WERank / matrix-info sources) whereas the current `variance_floor` cannot.

C. Keep the Gaussian-NLL programme regression head unchanged (model.py:290, losses.py:69-85,
   weight 1.0). It is not the collapse driver and it preserves interpretable programme readout for
   molecular prompting / the pearson selection metric (runner.py:224).

Optional (only if A+B under-fix): raise the regression target rank from 50 Hallmark to ~500-1000
HVG per-gene (BLEEP/PEaRL scale) by swapping the parquet at runner.py:172 / v2_tcga.yaml:7 and
widening `programme_dim` (model.py:189, 211-212; export.py:52). This is a larger change (retrains,
re-fits residualization, changes export/eval schema) and is NOT needed to fix the MECHANISM — do it
only if higher-rank REGRESSION is independently wanted. The mechanism fix is A+B.

Do NOT: keep `programme_neighbourhood_loss` and merely reweight it — any nonzero weight re-imposes
the low-rank Gram. Removing it and relying on supcon alone (losses.py:56-66) also fails, because
supcon uses z_biology-vs-z_biology on the same neighbour graph and does not add an independent
high-rank axis; the RNA-paired positive in (A) is what breaks the symmetry.

## Risks & scaling

- Batch-size sensitivity: InfoNCE rank on biology depends on negatives per batch. With the uncapped
  H-Optimus token budget the effective patient count per batch is small; the existing programme
  memory bank (training.py:70-116) stores DETACHED embeddings so it gives NO gradient — it can
  supply extra *negatives* for the new contrastive but positives must stay in-batch/attached.
  Mitigation: keep the top-8 in-batch positive mask (runner.py:141) as the primary signal; treat
  the bank as negatives only. Watch batch composition so a batch is not all-same-cancer (would
  make paired biology InfoNCE trivial and re-collapse).
- Leakage: unchanged. The RNA embedding is a real paired modality already used for identity
  InfoNCE (training.py:194); pairing biology to it introduces no new label and the residualization
  fit stays train-fold-only (runner.py:169-172).
- Separation term interaction: `whitened_cross_covariance` identity<->biology (losses.py:37-43)
  still applies; raising biology rank can only help it satisfy separation honestly.
- Test coverage gap: NO existing test catches rank collapse (confirmed). Add a CPU test using the
  `_batch()` fixture (tests/test_v21_model.py:19) asserting `torch.linalg.matrix_rank` (or an
  effective-rank/entropy-of-singular-values proxy) of z_biology over a synthetic multi-cancer batch
  is > a floor after a few steps, and a NaN guard on the new contrastive when a batch has zero
  positives (mirror the `< 3` / `not positive_mask.any()` guards at losses.py:48,57).
- Scaling to Lambda A10+: no added forward passes (A reuses out_rna); (B) is O(256^2) per batch,
  negligible. Optional HVG swap would increase target memory (1000-D vs 50-D) and re-fit cost but
  is out of scope for the mechanism fix.
- Metric continuity: `programme_mean_pearson` selection (runner.py:227) and export widths
  (export.py:54-55) are UNCHANGED by A+B, so downstream artifacts/contracts stay valid
  (contracts.validate_artifact isfinite still holds).
