# MORPHEUS-Rebase Literature Fleet — Lane Scoping

## Goal
Assemble **>=400 real, verifiable papers**, each **mapped to one or more of the five rebase axes**, converging on **>=3 novelty-vetted research paths** for MORPHEUS.

### Five rebase axes
- **A1** — Promptable unified representation + NL task auto-detection (infer/route the task vs hard-coded probes).
- **A2** — Identified, pathway-addressable slots enabling reliable prompting (identifiability, per-programme addressability).
- **A3** — NL<->biology grounding + emergent-knowledge elicitation **and its evaluation**.
- **A4** — Multimodal prompting: when to ENCODE a modality (proteomics/phospho/CNV/SNV) vs treat it as RAG context; frozen-trunk plug-in.
- **A5** — Interventional/causal-geometry queries (counterfactual perturbation/drug as a query, not a retrained classifier).

## Lanes (15, disjoint)

**l01_multimodal_repr — Multimodal representation learning & fusion.** General-domain joint embedding, contrastive/fusion architectures. NOT: domain-specific pathology or omics encoders (l02/l03), nor alignment *theory* (l09).

**l02_pathology_fm — Pathology / WSI foundation models.** Histology/WSI pretraining, tissue encoders, slide-level trunks. NOT: molecular/omics FMs (l03) or generic fusion methods (l01).

**l03_omics_fm — Omics & multi-omics foundation models.** Genomic/transcriptomic/proteomic pretrained models and multi-omics integration. NOT: pathology imaging (l02); NOT perturbation/causal omics (l12).

**l04_molecular_nl_prompting — Molecular & NL-prompting of scientific models.** Text-conditioned molecular/biological models, prompt formats for science. NOT: unified promptable *interfaces* with task routing (l05); NOT agentic tool-use (l13).

**l05_promptable_unified — Promptable / task-general unified interfaces + task auto-detection.** One-model-many-tasks routing, implicit task inference. NOT: the prompt-encoding of molecules (l04); NOT compositional instruction-following (l14).

**l06_emergence_eval — Emergent-capability & elicitation evaluation.** Measuring emergent knowledge, elicitation/probing methodology. NOT: general benchmark construction (l11); NOT discovery pipelines (l07).

**l07_ai_discovery — AI for biological discovery / hypothesis generation.** Model-driven hypothesis proposal and scientific finding. NOT: agentic multi-tool orchestration (l13); NOT emergence *measurement* (l06).

**l08_multimodal_rag — Multimodal RAG vs. encoded-modality integration.** Retrieval-augmented multimodal context, encode-vs-retrieve trade-offs. NOT: fusion-encoder architectures (l01); NOT agent tool-calling loops (l13).

**l09_alignment_identifiability — Cross-modal alignment, identifiability & disentanglement theory.** Formal identifiability, disentanglement, alignment guarantees. NOT: empirical fusion models (l01); NOT causal-representation interventions (l12).

**l10_decision_support — Decision-support methods (retrospective only).** Retrospective predictive/decision models and their methodology. NOT: prospective/clinical-trial claims; NOT benchmark/confound design (l11).

**l11_benchmarks_confound — Benchmarks & confound-aware evaluation.** Benchmark design, leakage/confound control, evaluation rigor. NOT: emergence-specific elicitation eval (l06); NOT decision-support modeling (l10).

**l12_interventional_causal — Interventional / causal representation & perturbation modeling.** Perturbation prediction, causal representation, counterfactual geometry. NOT: identifiability *theory* alone (l09); NOT descriptive omics FMs (l03).

**l13_agentic_science — Agentic LLM scientific workflows & tool-use.** Autonomous agents, planning, tool orchestration for science. NOT: single-model discovery (l07); NOT retrieval-only augmentation (l08).

**l14_compositional_instruction — Compositional & instruction-following multimodal.** Compositional generalization, multimodal instruction-following. NOT: task auto-detection/routing (l05); NOT molecular prompt formats (l04).

**l15_steelman_prior_art — Adversarial prior-art / steelman-the-null.** Strongest existing systems that already do what MORPHEUS claims; null hypotheses. NOT: constructive method surveys (all other lanes); it argues *against* novelty.
