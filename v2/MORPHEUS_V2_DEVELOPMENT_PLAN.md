# MORPHEUS V2 Local-First Implementation Specification

## Decision

V2 is a task-specialized tumour-state model, not a single latent vector forced to solve retrieval, biology and clinical prediction simultaneously. The all-patch MLP-CLIP result is the retrieval anchor. V2 succeeds only if it retains held-out-cancer retrieval within 5% relative of that anchor and reproducibly improves at least two non-retrieval tasks, calibration, or missing-modality robustness.

All code is written and tested locally in `morpheus/v2/` before transfer to Lambda. Remote scripts are execution wrappers, never the source of truth.

## Data And Protocol

The primary claim is Protocol C: 11 development cancers, three cancer-grouped inner folds, then one 22-cancer outer test. Random patient splits are diagnostic only; within-cancer retrieval is a separate shortcut-control protocol. All transforms, thresholds, programme graphs, prototype mixing, calibrators and reference libraries fit only on the relevant inner-training population.

The registry is patient-level and versioned. It records patient, cancer, split, source/sample/slide/image provenance, patch count, modality availability, survival labels and source digests. All slides and patches from one patient remain in one split. Every V2, V1 and baseline export must include patient ID, cancer, split, availability mask, view-specific states, residuals, uncertainty, programme outputs, and code/config/split/fit-population digests.

Regenerate Hallmark scores from raw RNA for every eligible WSI-RNA patient before programme evaluation. The prior table lacks outer-test coverage, so it cannot support an outer-test programme claim.

## Architecture

`morpheus/v2/model.py` is the authoritative model.

1. H-Optimus 1536-D patch vectors plus coordinates and slide IDs map to 512-D tokens.
2. Orthogonally initialized, normalized soft-slot pooling performs patch-to-local (8 slots), local-to-slide (8 slots), and slide-to-patient (16 slots) aggregation. Softmax is over valid token scores. Dynamic token batching controls step size; patches are never capped and every valid patch is visited once per epoch.
3. BulkFormer RNA and optional clinical/SNV/CNV inputs use typed adapters with explicit availability masks. The WSI-only path physically excludes non-WSI tokens.
4. A 512-D, four-layer, eight-head Query Former produces `z_identity` (256-D), `z_biology` (256-D), `z_context` (128-D), five modality residual states, uncertainty, and `z_patient`.
5. The identity projection may warm-start from an inner-fold-only all-patch MLP-CLIP anchor. Biology, context and residual heads initialize independently.

Clinical, SNV and CNV adapters activate only after coverage, source-alignment and TESSERA contamination gates pass. Protein/phosphoprotein, DepMap, LINCS, spatial and organoid tasks are separate external-confirmation stages.

## Training

`morpheus/v2/losses.py` and `morpheus/v2/training.py` define the sole supported objective and checkpoint contract.

| Stage | Epochs | Weights |
|---|---:|---|
| Warm-up | 1-4 | identity InfoNCE 1.0; programme Gaussian NLL 0.5 |
| Specialization | 5+ | identity 0.25; programme 1.0; programme-neighbour KL 0.20; programme supervised contrastive 0.20; identity-biology cross-covariance 0.01; RNA reconstruction at most 0.03 |
| Robustness | after stable validation retrieval | natural-mask subset consistency, uncertainty calibration and selective risk |

`z_biology` never receives paired-patient CLIP/InfoNCE. Programme positives are derived from train-fold programme targets only. Use AdamW, BF16, TF32, cosine decay, gradient clipping at 1.0, and checkpoints containing model, optimizer, scheduler, RNG, sampler state, manifest and selected metric. Save last, best-retrieval, best-programme and inner-CV Pareto checkpoints; never select on the outer test.

Before final training, calibrate total-token budgets 16,384, 32,768 and 65,536. Select the largest measured tokens/sec below 90% GPU memory; never use artificial memory occupancy as a performance target.

## Eight Mandatory Tasks

Every task runs on identical patient rows, masks, candidate sets, seeds and patient-bootstrap intervals for V2 and all baselines.

1. Bidirectional WSI-RNA retrieval: R@1/5/10, MRR, median rank, matched/unmatched cosine and within-cancer retrieval.
2. WSI-only molecular prompting: train-library Soft-kNN and direct ridge Hallmark prediction; Pearson, Spearman, R2, calibration, programme panels and matched random gene-set controls.
3. Global and train-only residual programme-state prediction: macro/per-cancer correlation and programme-neighbour retrieval.
4. Zero-shot cancer-label classification using public prompted descriptions and frozen patient states.
5. Few-shot unseen-cancer classification using frozen encoders and disjoint support prototypes for k=1,2,5,10,20.
6. Coverage-gated molecular phenotype prediction: high TMB, high aneuploidy, stage group and recurrent alterations; AUROC, AUPRC, balanced accuracy, macro-F1, Brier and ECE.
7. Survival/prognosis: clinical-only, WSI-only, RNA-only, available-genomic, early/late fusion and V2 heads; Harrell/IPCW C-index, time AUC, integrated Brier and calibration.
8. Missing-modality/next-best-test plus immune-discordance: selective risk, uncertainty reduction, recommendation utility, purity/cancer-controlled enrichment and hypothesis cards.

Sparse programme stability is a required control for task 8: compare V2 programmes to PCA, NMF, ICA, clustering and random sparse features using seed/fold stability, cancer leakage and pathway grounding.

## Zero- And Few-Shot Cancer Protocol

`morpheus/v2/text_prototypes.py` implements these tasks. Build one reviewed description per seen and unseen TCGA cancer using public knowledge only: tissue, broad histology, morphology, canonical molecular features and broad immune/stromal context. Prohibit TCGA counts, patient examples, test-derived statistics or data-derived wording. Record source URLs and review status.

Embed descriptions once with pinned `microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext`, attention-mask mean pooling and L2 normalization. Train a small mapper from seen-cancer patient states to this frozen space only on development cancers. Evaluate WSI-only, RNA-only, full available-modality, `z_identity`, `z_biology`, `z_context` and `z_patient` views in both candidate sets: 22 unseen cancers and all 33 cancers. Report top-1/3/5, macro-F1, balanced accuracy, calibration, per-cancer accuracy and confusion matrices.

Zero-shot controls are random descriptions, cancer-name-only prompts, full prompts, raw H-Optimus, raw BulkFormer, raw concatenation, CCA/MLP-CLIP and V1. Few-shot evaluation freezes all encoders, samples identical disjoint support/query episodes for every method, and optionally blends text/support prototypes with alpha tuned only in seen-cancer validation episodes.

## Baselines, Gates And Local Verification

Required retrieval baselines: no alignment, Ridge, CCA, Procrustes, all-patch MIL regression, MLP-CLIP, SigLIP, hard-negative CLIP and V1. Required fusion baselines: WSI-only, RNA-only, clinical-only, genomic-only when eligible, early concatenation, parameter-matched late fusion and a masked transformer without typed slots. V2 ablations: retrieval anchor only, identity-only, no biology separation, no residual/uncertainty, and full V2.

The word `disentangled` is allowed only if all final seeds pass CKA <= 0.80, mean absolute identity-biology cosine <= 0.70, 10-NN overlap <= 0.65, a ten-point biology cancer-probe reduction where appropriate, retained programme value and cross-head ablations. Otherwise call the states specialized.

Local tests must cover split leakage, fit populations, duplicate patient/slide membership, complete patch coverage, deterministic resume, WSI-only isolation, slot entropy/diversity, shared candidate universes, support/query disjointness, description provenance and end-to-end artifact schema production. Add `transformers` and `safetensors` before generating text embeddings. Transfer only a passing local source manifest, test report and smoke DAG to Lambda.
