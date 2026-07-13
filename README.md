# MORPHEUS

MORPHEUS is a leakage-controlled multimodal cancer representation framework for integrating fine-grained H&E pathology, RNA, genomics and clinical context to retrieve related patients, infer tumour states and support falsifiable biological hypotheses across held-out cancer types.

## Repository layout

- `src/data/` — patient registry construction, source auditing, cancer-held-out splits and H-Optimus patch-store utilities.
- `src/encoders/` — WSI, RNA, clinical and genomic adapters.
- `src/models/` — original MORPHEUS/BioQueryFormer model components.
- `src/training/` — training and baseline runners.
- `src/eval/` — shared retrieval, molecular prompting and diagnostic evaluation code.
- `tests/` — unit and protocol tests.
- `v2/` — the separate V2 implementation, design specification, uncapped patch runner and discovery-task evaluation layer. This directory is present on the `v2` branch.

## Scientific protocol

The primary protocol is cancer-held-out evaluation: development occurs only in the designated development cancer types, while all preprocessing, residualisation, model selection and reference libraries are fit within the relevant training population. All slides and all patches belonging to a patient remain in one split. Generated data, feature stores, model weights and run outputs are intentionally excluded from Git.

V2 uses hierarchical, uncapped H-Optimus patch tokens, isolated WSI/RNA identity alignment, programme-specialized biological states, and shared artifact-based evaluation against V1 and baseline representations. See [`v2/MORPHEUS_V2_DEVELOPMENT_PLAN.md`](v2/MORPHEUS_V2_DEVELOPMENT_PLAN.md) on the `v2` branch for the detailed protocol.

## Scope

MORPHEUS is a research codebase. Results should be interpreted under the documented split, coverage, contamination and task-specific claim gates; it is not a clinical decision system.

## Data and model access

The repository contains no patient-level data or model weights. Reproduction requires separately authorized access to the relevant TCGA-UT assets, H-Optimus embeddings/model access, and RNA/genomic source data. Run configuration paths are defined in `configs/` and are expected to point to local or persistent compute storage.
