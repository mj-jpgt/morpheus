## Omics & multi-omics foundation models

Lane `l03_omics_fm`. Remit: RNA/cell/gene foundation models, proteomics/phospho representation learning, multi-omics integration; representation of molecular state. Every entry maps to MORPHEUS rebase axes A1-A5:
- **A1** promptable unified representation + NL task auto-detection
- **A2** identified, pathway-addressable slots (identifiability / per-programme addressability)
- **A3** NL<->biology grounding + emergent-knowledge elicitation AND its evaluation
- **A4** multimodal prompting: when to ENCODE a modality vs treat as RAG context; frozen-trunk plug-in
- **A5** interventional/causal-geometry queries (perturbation/drug as a query, not a retrained classifier)

---

### 1. scGPT: toward building a foundation model for single-cell multi-omics using generative AI
Cui, Wang, Maan, Pang, Luo, Wang. **Nature Methods**, 2024. https://www.nature.com/articles/s41592-024-02201-0 (bioRxiv 2023.04.30.538439)

**Takeaway:** A single generative transformer pretrained on 33M cells that is fine-tuned into cell annotation, batch/multi-omic integration, perturbation prediction and GRN inference.

**Technical summary:** scGPT treats genes as tokens and cells as "sentences," pretrained with a masked/generative expression-prediction objective over 33M cells, using specialized attention for the non-sequential bag-of-genes structure. A single backbone transfers via fine-tuning to cell-type annotation, multi-batch and multi-omic (RNA+ATAC+protein) integration, genetic perturbation-response prediction, and gene-network inference. Attention maps are mined post hoc as candidate gene programmes.

**Plain-English:** Like GPT but for cells: it learns the "grammar" of which genes fire together, then can be nudged toward many different single-cell tasks with a little extra training.

**Applicability:** A1 (one backbone, many tasks — but tasks are hard-coded fine-tuning heads, not NL-routed, exactly the gap MORPHEUS targets), A4 (demonstrates encoding ATAC/protein into the same token space), A5 (perturbation head is a *fine-tuned classifier*, not a query — the thing MORPHEUS reframes). Design implication: scGPT is the canonical baseline; MORPHEUS must show that task routing + counterfactual queries beat per-task fine-tuning heads.

**Novelty implication:** Pre-empts any bare "unified single-cell representation" claim. MORPHEUS novelty must live in *promptable auto-routing (A1)* and *interventional queries (A5)*, not in multi-task capability per se.

---

### 2. Transfer learning enables predictions in network biology (Geneformer)
Theodoris, Xiao, Chopra, Chaffin, Al Sayed, Hill, Mannion, Ellinor. **Nature**, 2023. https://www.nature.com/articles/s41586-023-06139-9

**Takeaway:** Context-aware attention model pretrained on ~30M transcriptomes that encodes network hierarchy in its attention weights and transfers to low-data disease settings.

**Technical summary:** Cells are represented as rank-value-encoded gene tokens (genes ordered by expression rank), pretrained self-supervised on ~30M single-cell transcriptomes (later scaled to ~104M). Fine-tuning with limited labels enables chromatin/network-dynamic predictions and in-silico deletion analyses; applied to cardiomyopathy it nominated therapeutic targets validated in iPSC models.

**Plain-English:** By reading tens of millions of cells, the model builds an internal map of how genes regulate each other, which it can reuse to make predictions in diseases where data is scarce.

**Applicability:** A2 (rank encoding gives a defined, addressable gene-ordering — a partial identifiability story), A3 (attention weights claimed to encode network hierarchy = emergent-knowledge elicitation, but evaluation is anecdotal), A5 (in-silico gene deletion is a genuine interventional query in embedding space). Design implication: MORPHEUS's A5 counterfactual machinery should be benchmarked directly against Geneformer's in-silico deletion.

**Novelty implication:** Geneformer already frames *in-silico perturbation as elicitation from a frozen model* — strengthens the premise but pre-empts naive "first to do in-silico deletion" claims. MORPHEUS must differentiate on *rigorous evaluation of emergent knowledge (A3)* and *NL-addressable programmes (A2)*.

---

### 3. Scaling and quantization of a large-scale foundation model enables resource-efficient predictions in network biology
Geneformer team (Theodoris lab et al.). **Nature Computational Science**, 2026. https://www.nature.com/articles/s43588-026-00972-4

**Takeaway:** Scales Geneformer to ~104M-cell pretraining with quantization, showing compute-efficient inference without losing network-biology predictive power.

**Technical summary:** Expands the pretraining corpus to ~104M human single-cell transcriptomes across broader tissue/disease contexts and applies quantization to reduce memory/compute for downstream fine-tuning and in-silico perturbation. Reports that resource-efficient (quantized) variants retain accuracy on network-biology tasks, addressing deployability of large single-cell FMs.

**Plain-English:** A leaner, bigger Geneformer — trained on far more cells but compressed so it can run cheaply while still making good predictions.

**Applicability:** A4 (frozen-trunk plug-in economics: quantized frozen backbones make MORPHEUS's "frozen trunk + light adapters" pattern practical), A2 (larger corpus tests whether programme addressability survives scaling). Design implication: supports a frozen-quantized-trunk deployment story for MORPHEUS adapters.

**Novelty implication:** Neutral/strengthening — validates the frozen-trunk economics MORPHEUS assumes; not a novelty threat.

---

### 4. Large-scale foundation model on single-cell transcriptomics (scFoundation / xTrimoscFoundationα)
Hao, Gong, Cheng, Zhang, Zhou, Xie, Song. **Nature Methods**, 2024. https://www.nature.com/articles/s41592-024-02305-7 (bioRxiv 2023.05.29.542705)

**Takeaway:** 100M-parameter asymmetric transformer over ~50M cells and 20k genes with a read-depth-aware objective, SOTA on drug-response and perturbation tasks.

**Technical summary:** Uses an asymmetric encoder-decoder ("xTrimoGene") to model all ~20,000 genes per cell efficiently, with a read-depth-aware pretraining task that couples low- and high-depth cells to learn co-expression. Achieves SOTA on expression enhancement/imputation, bulk and single-cell drug-response prediction, perturbation prediction, cell annotation and gene-module inference; embeddings feed downstream models like GEARS and DeepCDR.

**Plain-English:** A very large cell model that reads every gene at once and corrects for how deeply each cell was sequenced, making it strong at predicting drug responses.

**Applicability:** A4 (embeddings are designed as a plug-in feature bank for external models — direct evidence for MORPHEUS's frozen-trunk-plus-adapter route), A5 (drug/perturbation response prediction). Design implication: scFoundation's "embedding-as-plug-in" is the RAG-context-vs-encode boundary MORPHEUS must formalize under A4.

**Novelty implication:** Pre-empts "large FM embeddings improve perturbation/drug prediction." MORPHEUS must move from *feeding embeddings to a fixed head* to *drug-as-a-query in causal geometry (A5)*.

---

### 5. Universal cell embedding provides a foundation model for cell biology (UCE)
Rosen, Roohani, Agrawal, Samotorčan, TSP Consortium, Quake, Leskovec. **Nature**, 2026 (bioRxiv 2023.11.28.568918). https://www.biorxiv.org/content/10.1101/2023.11.28.568918v3.full

**Takeaway:** Zero-shot, cross-species cell embeddings via ESM2 protein-based gene tokenization — no fine-tuning, robust to batch, generalizes to unseen species.

**Technical summary:** UCE tokenizes genes through their ESM2 protein-language-model embeddings, so any gene from any species maps into one shared latent space; trained self-supervised on 36M cells across 8 species. It embeds entirely new datasets/species zero-shot (no retraining), outperforms Geneformer/scGPT on integration/batch metrics, and shows emergent cell-lineage organization.

**Plain-English:** It represents every cell — from any animal — in one universal space by describing genes through their proteins, so it can place a brand-new species' cells on the map without extra training.

**Applicability:** A1 (true zero-shot embedding without task-specific fine-tuning — closest existing thing to a promptable universal representation), A3 (emergent lineage structure = emergent-knowledge claim), A4 (protein-sequence grounding of gene tokens is a multimodal-encode design pattern). Design implication: UCE proves zero-shot universality is possible; MORPHEUS's advance is making that space *NL-promptable and interventional*, not merely universal.

**Novelty implication:** Strong pre-empt on "universal zero-shot cell representation." MORPHEUS must not claim universality as novelty — the wedge is A1 NL-routing + A2 addressability on top of a universal space.

---

### 6. scBERT as a large-scale pretrained deep language model for cell-type annotation of single-cell RNA-seq data
Yang, Wang, Chen, Zhang, Cui, Yu, Liu, Wang, Yao. **Nature Machine Intelligence**, 2022. https://www.nature.com/articles/s42256-022-00534-z (bioRxiv 2021.12.05.471261)

**Takeaway:** The first BERT-style single-cell model — pretrain on unlabeled scRNA-seq, fine-tune for annotation with interpretable gene-attention.

**Technical summary:** Adapts BERT's masked-pretrain/fine-tune recipe with gene2vec-style embeddings and a Performer attention backbone to handle >16k genes. Pretrained on unlabeled scRNA-seq, then fine-tuned for cell-type annotation, showing robustness to batch effects, novel-cell-type discovery, and gene-attention interpretability.

**Plain-English:** The idea that started the field: treat a cell like a sentence of genes, pretrain like a language model, then specialize it to label cell types.

**Applicability:** A1 (early pretrain-then-finetune template MORPHEUS supersedes), A3 (gene-attention interpretability as a proto-elicitation method). Design implication: historical baseline; establishes that interpretability of attention is expected but under-evaluated.

**Novelty implication:** Foundational prior art; no direct threat, but anchors the "pretrain+finetune head" paradigm MORPHEUS argues against.

---

### 7. CellPLM: Pre-training of Cell Language Model Beyond Single Cells
Wen, Tang, Dai, Wang, Ding, Jin, Tang. **ICLR**, 2024 (bioRxiv 2023.10.03.560734). https://openreview.net/forum?id=BKXvPDekud

**Takeaway:** Models cells as tokens and tissues as sentences, adding cell-cell relations and a spatial-informed Gaussian-mixture latent prior.

**Technical summary:** CellPLM aggregates gene embeddings into cell tokens, then models inter-cell relations within tissues (using spatial transcriptomics in pretraining) via a transformer with a Gaussian-mixture prior as inductive bias against data scarcity/noise. It surpasses prior FMs on annotation, denoising/imputation, zero-shot embedding, and gene-perturbation prediction.

**Plain-English:** Instead of treating each cell in isolation, it also learns how neighboring cells relate within a tissue, using spatial data to teach it that context.

**Applicability:** A2 (structured latent prior is a step toward identifiable/addressable latent structure), A4 (spatial modality folded into pretraining), A5 (perturbation prediction). Design implication: the Gaussian-mixture prior is a concrete mechanism for MORPHEUS's A2 "identified slots" — a latent-structure design worth comparing against.

**Novelty implication:** Partially pre-empts A2 (structured latent) — MORPHEUS must show its slots are *pathway-addressable and promptable*, not just mixture components.

---

### 8. Simple and effective embedding model for single-cell biology built from ChatGPT (GenePT)
Chen, Zou. **Nature Biomedical Engineering**, 2024 (bioRxiv 2023.10.16.562533). https://www.nature.com/articles/s41551-024-01284-6

**Takeaway:** Gene/cell embeddings built purely from GPT text descriptions of genes rival Geneformer/scGPT without any expression pretraining.

**Technical summary:** GenePT embeds each gene from its NCBI text description via GPT-3.5, then forms cell embeddings by expression-weighted gene-embedding averages or by GPT sentence-embeddings of expression-ordered gene names. With no expression-data pretraining or curation, it matches or beats Geneformer/scGPT on gene-property and cell-type classification.

**Plain-English:** Surprisingly, describing genes in English (via ChatGPT) and combining those descriptions gives cell representations as good as models trained on millions of cells.

**Applicability:** A3 (direct evidence that NL/literature knowledge is a *substitute* for learned biology — core to NL<->biology grounding), A1 (language-native representation), A4 (text-as-modality vs expression-as-modality tradeoff). Design implication: GenePT is the strongest argument that MORPHEUS should treat literature/NL as first-class grounding (and sometimes RAG context) rather than always encoding expression.

**Novelty implication:** Reframes MORPHEUS's A3 — if plain text embeddings are "hard to beat," MORPHEUS must show that *grounded, interventional, expression-conditioned* NL beats text-only priors, or risks the "PCA/GPT-text is enough" critique.

---

### 9. Nicheformer: a foundation model for single-cell and spatial omics
Schaar, Tejada-Lapuerta, Palla, Gutgesell, Halle, Minaeva, Vornholz, Dony, Drummer, Bahrami, Theis. **Nature Methods**, 2025 (bioRxiv 2024.04.15.589472). https://www.nature.com/articles/s41592-025-02814-z

**Takeaway:** Trained on 110M dissociated + spatial cells (SpatialCorpus-110M), it transfers spatial context onto dissociated scRNA-seq.

**Technical summary:** A transformer pretrained jointly on ~57M dissociated and ~53M spatially resolved cells across 73 tissues, learning representations that encode niche/spatial context. Excels at spatial-composition and spatial-label prediction under linear probing and fine-tuning, and can predict the likely spatial context of dissociated cells.

**Plain-English:** By learning from both isolated cells and cells with known tissue locations, it can guess where a cell "belongs" in a tissue even when that spatial info was never measured.

**Applicability:** A4 (spatial modality as an encode target; linear-probing = frozen-trunk evaluation), A5 (predicting spatial context is a form of counterfactual "where would this cell sit?"). Design implication: Nicheformer's linear-probe protocol is a clean template for MORPHEUS's frozen-trunk adapter evaluation.

**Novelty implication:** Pre-empts "spatial context transfer." MORPHEUS's A4 must be about *deciding* when to encode spatial vs use it as retrieval context, not just encoding it.

---

### 10. Multi-omics single-cell data integration and regulatory inference with graph-linked embedding (GLUE / scGLUE)
Cao, Gao. **Nature Biotechnology**, 2022 (bioRxiv 2021.08.22.457275). https://www.nature.com/articles/s41587-022-01284-4

**Takeaway:** Aligns unpaired scRNA/scATAC/scmethylation into one embedding using a prior knowledge graph of feature relations, at million-cell scale.

**Technical summary:** GLUE uses modality-specific variational autoencoders tied by a guidance graph encoding prior regulatory relations (e.g., peak-to-gene), enabling unpaired multi-omics integration and simultaneous regulatory inference. It outperforms contemporaneous integration tools on accuracy/robustness/scalability and supports triple-omics integration for atlas building.

**Plain-English:** It stitches together different measurement types from different cells by using known biology (which regions regulate which genes) as the glue.

**Applicability:** A4 (the canonical "encode multiple omics into a shared space via prior graph" — a structured alternative to RAG context), A2 (guidance graph injects addressable prior structure), A5 (regulatory inference). Design implication: GLUE's guidance graph is a template for grounding MORPHEUS's A2 slots in known pathway/regulatory priors rather than unsupervised factors.

**Novelty implication:** Pre-empts unsupervised multi-omics integration novelty. MORPHEUS differentiates by making the integrated space *promptable and interventional*, not just aligned.

---

### 11. CellFM: a large-scale foundation model pre-trained on transcriptomics of 100 million human cells
Zeng et al. **Nature Communications**, 2025. https://www.nature.com/articles/s41467-025-59926-5

**Takeaway:** An 800M-parameter single-species (human) FM pretrained on ~100M cells, arguing depth-on-human beats broad cross-species dilution.

**Technical summary:** CellFM scales to ~100M human cells with a large (hundreds-of-millions-parameter) transformer using an efficient attention/low-rank design to cover the full gene panel. It reports gains on annotation, imputation, perturbation and gene-function tasks, positioning human-only scale as competitive with or better than multi-species FMs.

**Plain-English:** A very large model trained only on human cells, betting that going deep on one species beats spreading across many.

**Applicability:** A4/A2 (scale-vs-breadth tradeoff informs whether MORPHEUS should specialize or universalize its trunk), A5 (perturbation tasks). Design implication: evidence for a human-specialized trunk if MORPHEUS prioritizes clinical/pathway addressability over cross-species universality.

**Novelty implication:** Neutral — a scale point in the design space; not a direct claim threat.

---

### 12. A large-scale foundation model for bulk transcriptomes (BulkFormer)
Kang, Bo, et al. **Cell Systems**, 2026 (bioRxiv 2025.06.11.659222). https://www.cell.com/cell-systems/abstract/S2405-4712(26)00139-0 · https://github.com/KangBoming/BulkFormer

**Takeaway:** A ~150M-parameter bulk-RNA-seq FM (GNN + Performer hybrid) on ~500k bulk profiles that beats scRNA-trained FMs on bulk tasks at far lower cost.

**Technical summary:** BulkFormer combines a graph neural network capturing explicit gene-gene interactions with a Performer module for global expression dependencies, pretrained on >500k human bulk transcriptomes covering ~20k protein-coding genes. It outperforms scRNA-seq FMs across imputation, disease annotation, prognosis, drug-response, compound-perturbation simulation, and gene-essentiality scoring — highlighting that bulk needs its own FM.

**Plain-English:** Most cell models are trained on single-cell data; this one is built for bulk RNA (whole-tissue) data and does clinical/pharma tasks better and cheaper.

**Applicability:** A4 (bulk vs single-cell as distinct modalities — when to encode which; hybrid GNN grounds gene-gene structure), A5 (compound-perturbation simulation, gene essentiality as interventional queries). Design implication: MORPHEUS's A4 modality-routing should explicitly cover bulk RNA (clinical reality) not just scRNA; BulkFormer is a strong bulk baseline.

**Novelty implication:** Pre-empts "one FM for all transcriptomics." Supports MORPHEUS's argument that modality/assay identity matters and must be routed (A4).

---

### 13. Benchmarking Transcriptomics Foundation Models for Perturbation Analysis: one PCA still rules them all
Bendidi, Whitfield, Kenyon-Dean, Ben Yedder, El Mesbahi, Noutahi, Denton. **arXiv:2410.13956**, 2024 (NeurIPS workshop). https://arxiv.org/abs/2410.13956

**Takeaway:** Under fair evaluation, scVI and plain PCA beat single-cell foundation-model embeddings for understanding biological perturbations.

**Technical summary:** The authors argue current perturbation-analysis benchmarks are inconsistent/leaky, propose a more rigorous protocol, and find that classical baselines (PCA) and scVI outperform popular transcriptomics FMs at capturing perturbation effects. The result questions whether large-scale pretraining currently buys anything for the flagship perturbation use case.

**Plain-English:** When you test them fairly, the fancy foundation models don't beat a decades-old dimensionality-reduction trick at predicting how cells respond to being perturbed.

**Applicability:** A5 (direct evidence that FM embeddings are *not yet* good at perturbation geometry — the exact gap MORPHEUS's causal-query axis must close), A3 (benchmarking rigor). Design implication: MORPHEUS MUST include PCA/scVI baselines and a leakage-controlled protocol, and demonstrate its A5 machinery beats them — otherwise the whole premise is vulnerable.

**Novelty implication:** Reframes the field — the biggest risk AND opportunity for MORPHEUS. Strengthens the *motivation* (existing FMs fail at causal/perturbation queries) but sets a hard bar: beat PCA/scVI or the A5 claim collapses.

---

### 14. BMFM-RNA: An Open Framework for Building and Evaluating Transcriptomic Foundation Models
IBM Research (Dandala et al.). **arXiv:2506.14861**, 2025. https://arxiv.org/abs/2506.14861

**Takeaway:** An open, modular framework unifying training objectives/architectures so transcriptomic FMs can be compared apples-to-apples.

**Technical summary:** BMFM-RNA provides a reproducible pipeline covering data processing, multiple pretraining objectives (masked, contrastive, generative), and standardized downstream evaluation, letting researchers ablate design choices under one roof. It reproduces and stress-tests existing FMs and exposes how much reported gains depend on evaluation setup.

**Plain-English:** A common testbench so people stop comparing cell models unfairly and can see which design choices actually matter.

**Applicability:** A3 (evaluation infrastructure for emergent-knowledge and task claims), A1 (multi-objective backbone comparison). Design implication: adopt BMFM-RNA-style standardized eval so MORPHEUS's A1/A3 claims are credible and reproducible.

**Novelty implication:** Neutral tooling; raises the evaluation bar MORPHEUS must clear.

---

### 15. Predicting transcriptional outcomes of novel multigene perturbations with GEARS
Roohani, Huang, Leskovec. **Nature Biotechnology**, 2024 (bioRxiv 2022). https://www.nature.com/articles/s41587-023-01905-6

**Takeaway:** A GNN over a gene co-expression + GO knowledge graph predicts expression outcomes of *unseen* single and combinatorial genetic perturbations.

**Technical summary:** GEARS represents each gene with a knowledge-graph-informed embedding and predicts post-perturbation expression by composing perturbation embeddings, enabling zero-shot prediction for genes/combinations never perturbed in training, including non-additive (epistatic) effects. It outperforms prior perturbation-response models and flags genetic-interaction types.

**Plain-English:** Given a gene knockout it has never seen — even two at once — it predicts how the whole cell's expression will shift, using a graph of known gene relationships.

**Applicability:** A5 (the reference model for perturbation-as-prediction, including unseen combinations — MORPHEUS's causal-query axis must engage it), A2 (knowledge-graph gene embeddings = addressable prior structure). Design implication: GEARS is the A5 target to beat/subsume; MORPHEUS should express GEARS-style queries as NL-promptable counterfactuals over a frozen trunk.

**Novelty implication:** Strong pre-empt on "predict novel perturbation outcomes." MORPHEUS must reframe from a *dedicated perturbation model* to *perturbation as a query on a general representation* to be novel.

---

### 16. scGen predicts single-cell perturbation responses
Lotfollahi, Wolf, Theis. **Nature Methods**, 2019. https://www.nature.com/articles/s41592-019-0494-8

**Takeaway:** A VAE + latent-space vector arithmetic generalizes perturbation/stimulation responses to unseen cell types.

**Technical summary:** scGen learns a VAE latent space where a perturbation is a difference vector; adding that vector transports control cells of an unseen cell type to their predicted perturbed state (out-of-distribution generalization). It predicted infection/stimulation responses across species and cell types.

**Plain-English:** It learns "the direction of a perturbation" in a compressed space, then applies that same shift to new cell types to guess their response.

**Applicability:** A5 (latent vector arithmetic = the original "counterfactual as geometry" idea MORPHEUS generalizes), A2 (perturbation as an interpretable latent direction). Design implication: scGen's latent-direction abstraction is the conceptual seed of MORPHEUS's causal-geometry axis; MORPHEUS scales it to a promptable, multi-programme setting.

**Novelty implication:** Prior art for "perturbation = latent vector." MORPHEUS novelty is making these directions *identified, addressable (A2) and NL-queried (A1)* rather than one hand-defined vector.

---

### 17. Predicting cellular responses to complex perturbations at scale with the Compositional Perturbation Autoencoder (CPA)
Lotfollahi, Klimovskaia Susmelj, De Donno, ... Theis. **Molecular Systems Biology**, 2023 (bioRxiv 2021). https://www.embopress.org/doi/full/10.15252/msb.202211517

**Takeaway:** Disentangles cell state from perturbation and covariate effects into composable latent factors, enabling dose/combination interpolation.

**Technical summary:** CPA is an autoencoder with adversarially disentangled latent spaces for basal state, drug/genetic perturbation, dose, and covariates, so responses to unseen drug combinations and doses can be composed and interpolated. It generalizes across perturbations, doses, and cell lines and yields interpretable perturbation embeddings.

**Plain-English:** It separates "what kind of cell this is" from "what was done to it," so you can mix-and-match perturbations and doses to predict new combinations.

**Applicability:** A2 (explicit disentanglement into addressable factors — a direct identifiability/addressability precedent), A5 (compositional counterfactual queries over dose/combination). Design implication: CPA's disentanglement is the strongest existing model of MORPHEUS's A2 "pathway-addressable slots"; MORPHEUS must show its slots are more general and NL-addressable.

**Novelty implication:** Significant pre-empt on A2+A5 combined. MORPHEUS must clearly exceed CPA on *number/granularity of identified programmes* and *NL prompting*, or A2 reads as incremental.

---

### 18. Joint probabilistic modeling of single-cell multi-omic data with totalVI
Gayoso, Steier, Lopez, Regier, Nazor, Streets, Yosef. **Nature Methods**, 2021. https://www.nature.com/articles/s41592-020-01050-x

**Takeaway:** A deep generative model jointly embedding CITE-seq RNA + surface-protein counts with calibrated uncertainty and denoising.

**Technical summary:** totalVI is a VAE that jointly models paired transcript and antibody-derived-tag (protein) counts, handling protein background/ambient noise and enabling integrated embedding, imputation of missing proteins, and differential expression with posterior uncertainty. It is a scvi-tools reference for RNA+protein integration.

**Plain-English:** For assays that measure both RNA and surface proteins in the same cell, it learns one combined representation and cleans up the noisy protein signal.

**Applicability:** A4 (canonical "encode protein modality jointly with RNA" — the encode side of MORPHEUS's encode-vs-RAG decision), A2 (probabilistic latent with uncertainty). Design implication: totalVI is the baseline for *encoding* proteomics; MORPHEUS's A4 must justify when a frozen-trunk RAG treatment beats joint encoding like this.

**Novelty implication:** Pre-empts "joint RNA+protein encoding." Supports MORPHEUS's A4 framing that encoding is one option among several, with a decision rule.

---

### 19. MultiVI: deep generative model for the integration of multimodal data
Ashuach, Gabitto, Koodli, Saldi, Jordan, Yosef. **Nature Methods**, 2023 (bioRxiv 2021). https://www.nature.com/articles/s41592-023-01909-9

**Takeaway:** Probabilistically integrates paired and unpaired RNA + ATAC (+ protein) into a shared latent, imputing missing modalities.

**Technical summary:** MultiVI builds modality-specific encoders sharing a joint latent, aligning paired multiome cells while also placing unpaired scRNA-only or scATAC-only cells in the same space, and imputes the missing modality with uncertainty. It extends the scvi-tools family to mosaic multi-omic integration.

**Plain-English:** It merges cells measured in different ways — some with RNA, some with chromatin, some with both — into one space and fills in what each cell is missing.

**Applicability:** A4 (mosaic/partial-modality integration — realistic clinical setting where MORPHEUS must handle missing modalities), A2 (shared latent). Design implication: MORPHEUS's A4 must handle mosaic availability (not every sample has every omic); MultiVI is the encode-based baseline.

**Novelty implication:** Pre-empts mosaic multi-omics integration. MORPHEUS differentiates on promptability/queries, not on integration itself.

---

### 20. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data
Argelaguet, Arnol, Bredikhin, Deloro, Velten, Marioni, Stegle. **Genome Biology**, 2020. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1

**Takeaway:** Unsupervised factor analysis that decomposes multi-omic/multi-group data into interpretable, sparse latent factors.

**Technical summary:** MOFA+ extends group factor analysis with stochastic variational inference to scale to single-cell multi-omics, learning sparse factors that are annotatable to specific molecular programmes and can capture group/condition structure. Factors are interpretable via their feature loadings across omics layers.

**Plain-English:** A principled way to find a handful of hidden "axes of variation" that span RNA, chromatin, protein, etc., each interpretable as a biological programme.

**Applicability:** A2 (the classical identifiable/interpretable-factor model — the statistical benchmark for MORPHEUS's "identified pathway-addressable slots"), A4 (multi-omic decomposition). Design implication: MORPHEUS's A2 identifiability claims should be benchmarked against MOFA+ factors and their interpretability, and argue what a learned FM adds over linear factor analysis.

**Novelty implication:** Reframes A2 — if linear MOFA+ factors are already interpretable and addressable, MORPHEUS must show its slots capture *nonlinear, prompt-conditioned, causal* structure MOFA+ cannot.

---

### 21. Prosit: proteome-wide prediction of peptide tandem mass spectra by deep learning
Gessulat, Schmidt, Zolg, ... Wilhelm, Kuster. **Nature Methods**, 2019. https://www.nature.com/articles/s41592-019-0426-7

**Takeaway:** Deep learning predicts peptide fragmentation spectra and retention time, transforming proteomics identification/quantification.

**Technical summary:** Prosit is an encoder-decoder (bidirectional recurrent) network trained on >500k tryptic peptides to predict fragment-ion intensities and iRT retention times as a function of sequence, charge, and collision energy. Predicted spectral libraries boost identification rates and enable data-independent-acquisition workflows.

**Plain-English:** It learns the "chemical fingerprint" a peptide will produce in a mass spectrometer, so proteomics software can identify proteins far more reliably.

**Applicability:** A4 (representation learning for the proteomics modality — grounding how peptide/protein signals should be encoded before integration), A3 (sequence->measurement grounding). Design implication: informs MORPHEUS's proteomics-encoding path (A4): peptide-level DL representations are mature and could feed a frozen adapter rather than being re-learned.

**Novelty implication:** Neutral; establishes proteomics DL representations exist and are strong, supporting an "encode via existing proteomics models" option in A4.

---

### 22. DeepPhospho accelerates DIA phosphoproteome profiling through in-silico library generation
Lou, Shui, Li, ... Yang. **Nature Communications**, 2021. https://www.nature.com/articles/s41467-021-26979-1

**Takeaway:** A deep model predicts phosphopeptide spectra and retention time, enabling large-scale DIA phosphoproteomics.

**Technical summary:** DeepPhospho uses transformer/LSTM sequence models trained on phosphopeptide libraries to predict fragment intensities and iRT specifically for phosphorylated peptides, generating in-silico spectral libraries that substantially increase phosphosite identifications in DIA data. It addresses the modality-specific difficulty that phospho-modification shifts fragmentation.

**Plain-English:** Phosphorylation changes how peptides shatter in the mass spec; this model learns those altered patterns so far more phosphorylation sites can be detected.

**Applicability:** A4 (phospho as a distinct modality needing dedicated representation — directly informs MORPHEUS's phospho encode-vs-RAG decision), A2 (site-level addressability). Design implication: phospho signals are modality-specific enough that MORPHEUS should either encode via a phospho-specialized adapter or treat phospho-site facts as RAG context — DeepPhospho evidences the "specialized encoder" pole.

**Novelty implication:** Neutral; supports A4's claim that phospho cannot be naively folded into an RNA trunk.

---

### 23. Phosformer: an explainable transformer model for protein-kinase-specific phosphorylation predictions
Zhou, Yeung, Kannan, et al. **Bioinformatics**, 2023. https://academic.oup.com/bioinformatics/article/39/2/btad046/7000004

**Takeaway:** A transformer predicts kinase-specific phosphorylation sites with attention-based explainability of substrate motifs.

**Technical summary:** Phosformer jointly encodes kinase domain and substrate peptide sequences in a transformer to predict whether a given kinase phosphorylates a given site, generalizing across the kinome, and exposes attention over sequence positions for motif interpretability. It provides kinase-substrate relationship representations learnable from sequence alone.

**Plain-English:** Given a kinase and a candidate site, it predicts whether that kinase would tag it — and shows which residues drove the call.

**Applicability:** A2 (kinase-substrate = addressable signaling-programme structure — a phospho-pathway addressability precedent), A3 (attention-as-explanation, an elicitation-evaluation angle), A4 (sequence-based phospho representation). Design implication: kinase-substrate graphs are candidate priors for MORPHEUS's A2 signaling-pathway slots; Phosformer shows they are learnable and addressable.

**Novelty implication:** Neutral/strengthening for A2 in the signaling domain; not a whole-representation threat.

---

### 24. Evolutionary-scale prediction of atomic-level protein structure with a language model (ESM-2 / ESMFold)
Lin, Akin, Rao, ... Rives. **Science**, 2023. https://www.science.org/doi/10.1126/science.ade2574

**Takeaway:** A 15B-parameter protein language model learns structure-predictive representations from sequence alone, enabling fast structure prediction.

**Technical summary:** ESM-2 is a masked protein language model scaled to 15B parameters whose internal representations encode tertiary structure, powering ESMFold for single-sequence structure prediction without MSAs. The embeddings are widely reused as frozen protein features (as in UCE's gene tokenization).

**Plain-English:** By reading hundreds of millions of protein sequences like text, it learns enough to predict 3D structure and produces protein "fingerprints" others reuse.

**Applicability:** A4 (the reference frozen protein-representation used to ground gene/protein tokens — the plug-in modality pattern), A3 (emergent structural knowledge from sequence). Design implication: MORPHEUS's protein/proteomics grounding (A4) can lean on frozen ESM-2 embeddings rather than learning protein representations from scratch, as UCE already does.

**Novelty implication:** Neutral enabling tech; reinforces A4's frozen-plug-in strategy.

---

### 25. Simulating 500 million years of evolution with a language model (ESM3)
Hayes, Rao, ... Rives (EvolutionaryScale). **Science**, 2025 (bioRxiv 2024). https://www.science.org/doi/10.1126/science.ads0018

**Takeaway:** A multimodal generative protein model over sequence, structure, and function that can *design* novel functional proteins as a generative query.

**Technical summary:** ESM3 is a masked generative transformer trained jointly over protein sequence, structure, and function tokens, enabling prompted generation/infilling across modalities; it generated a novel fluorescent protein far from known sequences. It demonstrates promptable multimodal generation and function-conditioned design.

**Plain-English:** A protein model you can prompt across sequence, shape, and function to invent new proteins — it created a glowing protein unlike any in nature.

**Applicability:** A1 (promptable multimodal generation — the clearest existing "prompt the model to do the task" precedent, in proteins), A4 (joint sequence/structure/function tokens), A5 (generative design = a constructive counterfactual query). Design implication: ESM3 is the proof-of-concept that MORPHEUS's A1 promptability and A5 generative-query goals are achievable — but in proteins, not cell state; MORPHEUS ports the idea to molecular-state/pathway space.

**Novelty implication:** Strong conceptual pre-empt on A1/A5 *mechanism* (promptable multimodal generation). MORPHEUS's novelty must be the *domain transfer* to cell-state/pathway prompting + NL task auto-detection, not the promptable-generation idea itself.

---

### 26. Cell2Sentence: Teaching Large Language Models the Language of Biology (C2S / C2S-Scale)
Levine, Lévy, Rizvi, ... van Dijk. **ICML**, 2024 (bioRxiv 2023.09.11.557287); C2S-Scale (Yale/Google) 2025 scale-up. https://openreview.net/forum?id=EWt5wsEdvc

**Takeaway:** Represents each cell as a "sentence" of rank-ordered gene names so a standard LLM can be fine-tuned to do single-cell tasks in natural language.

**Technical summary:** C2S converts expression profiles into text sequences (genes ordered by expression) and fine-tunes pretrained LLMs to generate cells, annotate types, and answer questions, unifying generation and NL interaction; C2S-Scale extends this to multi-billion-parameter models and reports emergent single-cell reasoning. It bridges cell state and natural language directly.

**Plain-English:** It turns a cell into a sentence of gene names and teaches ChatGPT-style models to read and write cells, so you can talk to your single-cell data in English.

**Applicability:** A1 (NL-native single-cell interaction — the closest existing approach to MORPHEUS's promptable/task-auto-detect goal), A3 (NL<->biology grounding and emergent reasoning claims), A5 (generative cell prediction). Design implication: C2S is the most direct competitor to MORPHEUS's A1/A3 thesis; MORPHEUS must show its identified-slot/interventional approach beats "just tokenize cells as text for an LLM."

**Novelty implication:** Biggest A1/A3 pre-empt in the lane. MORPHEUS must differentiate on A2 (identifiable addressable slots) and A5 (causal queries) — pure NL-tokenization of cells is already claimed.

---

### 27. scPRINT: pre-training on 50 million cells allows robust gene-network predictions
Kalfon, Samaran, Peyré, Cantini. **Nature Communications**, 2025 (bioRxiv 2024.07.29.605556). https://www.nature.com/articles/s41467-025-58699-1

**Takeaway:** A denoising/generative FM whose attention yields cell-type-specific gene regulatory networks in zero shot.

**Technical summary:** scPRINT is pretrained with denoising and bottleneck objectives on ~50M cells and extracts gene-gene regulatory networks directly from attention without task-specific supervision, benchmarked against GRN-inference methods. It targets network inference as a first-class output rather than a downstream fine-tune.

**Plain-English:** Trained to clean up noisy cells, it also learns which genes regulate which — producing regulatory network maps for a given cell type without extra training.

**Applicability:** A3 (emergent GRN knowledge + its explicit evaluation against GRN benchmarks — a model for MORPHEUS's A3 elicitation-evaluation), A2 (gene-network structure as addressable). Design implication: scPRINT's GRN-extraction-with-benchmark is a template for how MORPHEUS should *evaluate* emergent biological knowledge (A3) rather than assert it.

**Novelty implication:** Pre-empts "emergent GRN from attention." Strengthens A3's *evaluation* demand; MORPHEUS must go beyond GRN to broader prompt-elicited knowledge with quantitative eval.

---

### 28. Toward universal cell embeddings: integrating single-cell RNA-seq datasets across species with SATURN
Rosen, Brbić, Roohani, ... Leskovec. **Nature Methods**, 2024 (bioRxiv 2023). https://www.nature.com/articles/s41592-024-02191-z

**Takeaway:** Cross-species integration via protein-embedding "macrogenes," coupling homology to a shared cell space.

**Technical summary:** SATURN maps genes from different species into shared "macrogenes" defined by protein-language-model (ESM) embedding similarity, then learns a common cell embedding across species. It integrates datasets lacking one-to-one orthologs and enables cross-species cell-type comparison.

**Plain-English:** It groups genes across species by protein similarity so cells from, say, frog and human can live in the same map even when genes don't line up one-to-one.

**Applicability:** A4 (protein-grounded gene bridging — a multimodal grounding pattern), A2 (macrogenes as addressable meta-features). Design implication: macrogenes are a concrete construction for MORPHEUS's A2 addressable units that are biologically grounded (via protein), not arbitrary latent dims.

**Novelty implication:** Pre-empts protein-grounded cross-species universality (with UCE). MORPHEUS should treat protein grounding as accepted machinery, not novelty.

---

### 29. Contextual AI models for single-cell protein biology (PINNACLE)
Li, Zitnik, et al. **Nature Methods**, 2024 (bioRxiv 2023). https://www.nature.com/articles/s41592-024-02341-3

**Takeaway:** Learns context-specific protein representations across 156 cell types via protein-interaction networks + single-cell data, improving therapeutic-target prediction.

**Technical summary:** PINNACLE is a graph neural network producing cell-type-contextualized protein embeddings from multi-organ single-cell data and protein-protein interaction networks, yielding hundreds of thousands of context-aware protein representations. It improves 3D-structure-based binding-site prediction and nominates cell-type-specific therapeutic targets.

**Plain-English:** The same protein behaves differently in different cell types; this model gives each protein a representation tailored to its cellular context, sharpening drug-target predictions.

**Applicability:** A4 (protein biology encoded with cellular context — the proteomics side of MORPHEUS's multimodal representation), A2 (per-cell-type protein addressability), A5 (target nomination as an actionable query). Design implication: PINNACLE shows context-specific protein slots are learnable — a template for MORPHEUS's A2 proteomics-programme addressability and A4 protein integration.

**Novelty implication:** Pre-empts "context-specific protein representation for targets." MORPHEUS differentiates by unifying protein context with promptable NL queries rather than a fixed GNN.

---

### 30. scMulan: a multitask generative pre-trained language model for single-cell analysis
Bian, Xu, Chen, ... Zhang, Jiang. **RECOMB / bioRxiv 2024.01.25.577152**, 2024. https://www.biorxiv.org/content/10.1101/2024.01.25.577152

**Takeaway:** Encodes cells as multi-attribute "cell-language" sentences (expression + metadata) so one generative model does annotation, batch integration, and conditional cell generation.

**Technical summary:** scMulan tokenizes not just genes but cell metadata (tissue, cell type, technology) into a unified "cell sentence," pretrained generatively on ~10M cells, enabling multitask inference (annotation, integration) and conditional generation prompted by attributes. Tasks are specified by prompt-like attribute conditioning rather than separate fine-tuned heads.

**Plain-English:** It writes cells as sentences that include both their genes and their labels, so you can prompt it (e.g., "generate a liver T cell") and get many tasks from one model.

**Applicability:** A1 (attribute-conditioned prompting = a partial version of MORPHEUS's task auto-detection/routing), A5 (conditional cell generation as a query). Design implication: scMulan's attribute-conditioning is a stepping stone toward A1; MORPHEUS must add free-form NL task *inference* (not fixed attribute slots) and interventional queries.

**Novelty implication:** Partial pre-empt on A1 (prompt-conditioned multitask). MORPHEUS's edge is NL task auto-detection and identified pathway slots beyond metadata attributes.

---

### 31. LangCell: Language-Cell Pre-training for Cell Identity Understanding
Zhao, Yang, Sun, ... Yao, Wang. **ICML**, 2024 (arXiv 2405.06708). https://arxiv.org/abs/2405.06708

**Takeaway:** Jointly pretrains cell-expression and natural-language descriptions so cells can be understood/annotated in zero shot from text.

**Technical summary:** LangCell contrastively/generatively aligns single-cell expression encoders with textual descriptions of cell identity (type, pathway, disease), enabling zero-shot cell-type annotation and text-based retrieval by grounding cell embeddings in language. It targets the cell-text grounding gap directly.

**Plain-English:** It teaches a cell model and a language model to agree, so you can describe a cell type in words and have it recognized without labeled examples.

**Applicability:** A3 (explicit NL<->cell grounding with zero-shot evaluation — core to MORPHEUS's grounding axis), A1 (text-queried annotation). Design implication: LangCell is a direct A3 baseline; MORPHEUS must show richer grounding (pathway/programme-level, interventional) than cell-identity text alignment.

**Novelty implication:** Pre-empts "align cells with natural-language identity for zero-shot." MORPHEUS differentiates on eliciting/evaluating *emergent* knowledge and on interventional grounding, not static identity matching.

---

### 32. AIDO.Cell: scaling dense representations for single cells with transcriptome-scale context
GenBio AI (Ho, Ellington, et al.). **bioRxiv / ICML 2024 workshop**, 2024. https://www.biorxiv.org/content/10.1101/2024.11.28.625303

**Takeaway:** A modular "AI-driven digital organism" cell FM scaled to hundreds of millions of parameters with full-transcriptome context, positioned as a composable module.

**Technical summary:** AIDO.Cell scales bidirectional transformers to model all ~20k genes per cell with dense (non-sparsified) attention over transcriptome-scale context, part of GenBio's modular multi-scale foundation-model stack meant to interoperate with DNA/protein modules. It reports competitive annotation, perturbation, and integration performance and emphasizes composability across biological scales.

**Plain-English:** A big cell model built as one Lego brick in a larger system meant to connect DNA, RNA, protein, and cell models together.

**Applicability:** A4 (explicit multi-scale/multimodal composability — the architectural vision closest to MORPHEUS's frozen-trunk + modality plug-ins), A2 (dense full-transcriptome context). Design implication: AIDO's modular stack is a competing blueprint for MORPHEUS's A4 plug-in architecture; MORPHEUS must justify NL-promptable routing over static module composition.

**Novelty implication:** Pre-empts "modular multimodal biological FM stack." MORPHEUS's wedge is A1 NL task-routing + A5 causal queries layered on such a stack.

---

### 33. Generative pretraining from pan-cancer/pan-tissue transcriptomes for cell representation (tGPT)
Shen, Liu, Yang, ... Chen. **iScience**, 2023. https://www.cell.com/iscience/fulltext/S2589-0042(23)00854-7

**Takeaway:** An autoregressive transformer over gene-expression rankings pretrained on ~22M samples yields transferable cell/bulk representations.

**Technical summary:** tGPT applies GPT-style autoregressive next-token prediction over rank-ordered gene tokens across ~22M single-cell and bulk samples, producing embeddings transferable to cell clustering, bulk tissue analysis, and survival/tumor characterization. It is an early demonstration that pure generative pretraining transfers across single-cell and bulk regimes.

**Plain-English:** A GPT for expression rankings that works on both single cells and bulk tissue, giving reusable representations for clustering and cancer analysis.

**Applicability:** A4 (bridges single-cell and bulk in one generative model — relevant to MORPHEUS's cross-assay routing), A1 (generative pretraining backbone). Design implication: supports a shared trunk spanning bulk+single-cell (relevant to clinical bulk data) under A4.

**Novelty implication:** Neutral early prior art; reinforces feasibility of a unified single-cell/bulk trunk.

---

## Cross-cutting synthesis for MORPHEUS

- **A1 (promptable + task auto-detection):** Closest prior art is C2S/C2S-Scale (#26), scMulan (#30), ESM3 (#25) — all show *prompt-conditioned* generation, but none do free-form NL *task auto-detection/routing*. This is the clearest open lane.
- **A2 (identified, pathway-addressable slots):** CPA (#17), MOFA+ (#20), PINNACLE (#29), Phosformer (#23) collectively occupy the identifiability/addressability space with disentanglement, linear factors, and context-specific protein/kinase slots. MORPHEUS must show slots that are simultaneously *identified, nonlinear, pathway-grounded, and NL-addressable* — no single prior does all four.
- **A3 (grounding + emergent-knowledge eval):** GenePT (#8) and the perturbation benchmark (#13) are warnings — text priors and PCA are hard to beat, and emergent-knowledge claims are often under-evaluated. scPRINT (#27) and LangCell (#31) model good evaluation practice. MORPHEUS must ship *quantitative emergent-knowledge evaluation*, not attention anecdotes.
- **A4 (encode vs RAG; frozen plug-in):** totalVI/MultiVI/GLUE (#18/#19/#10) = encode side; ESM-2/SATURN/PINNACLE (#24/#28/#29) = protein grounding; BulkFormer (#12) and DeepPhospho/Phosformer (#22/#23) show modality-specific encoders. No prior formalizes a *decision rule* for encode-vs-RAG — MORPHEUS's opportunity.
- **A5 (interventional/causal queries):** Geneformer in-silico deletion (#2), scGen (#16), CPA (#17), GEARS (#15), ESM3 design (#25). These treat perturbation as prediction/design, but #13 shows FMs currently *lose to PCA/scVI* on perturbation. MORPHEUS's A5 is the highest-risk, highest-reward axis: the motivation is validated by the failure results, but the bar (beat PCA/scVI under leakage control) is concrete and hard.
