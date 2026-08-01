## Pathology / WSI foundation models

Lane remit: WSI/H&E foundation models and vision-language pathology (UNI, Virchow(2), CONCH, TITAN, PLIP, MUSK, Prov-GigaPath, H-Optimus, CHIEF, Phikon and relatives) — their pretraining, prompting/zero-shot, and molecular-prediction capabilities. Every entry is mapped to MORPHEUS's five rebase axes: A1 (promptable unified representation + NL task auto-detection), A2 (identified, pathway-addressable slots / identifiability), A3 (NL<->biology grounding + emergent-knowledge elicitation and its evaluation), A4 (multimodal prompting: encode-vs-RAG a modality; frozen-trunk plug-in), A5 (interventional/causal-geometry / counterfactual-as-query).

Note on grounding: axis mappings below argue by analogy to MORPHEUS (a transcriptomics/cell-state programme model) from the pathology-imaging modality. The recurring cross-modal lesson is that pathology FMs have largely NOT built the promptable/identifiable/interventional capabilities MORPHEUS claims, which is where the novelty risk (or opportunity) lies.

---

### Tile / patch encoders (self-supervised)

1. **Towards a general-purpose foundation model for computational pathology (UNI)** — Chen, Ding, Lu, Williamson, ... Mahmood. *Nature Medicine*, 2024. https://www.nature.com/articles/s41591-024-02857-3
- Takeaway: A ViT-L/16 DINOv2 tile encoder pretrained on ~100M patches from >100k WSIs (Mass-100K) that set the reference bar for label-efficient, general-purpose pathology features.
- Technical summary: UNI is self-supervised (DINOv2) on ~100M tissue patches from >100,000 diagnostic H&E slides spanning 20 tissue types, with no labels. Evaluated on 34 clinical tasks, it excels at few-shot linear probing and enables resolution-agnostic ROI classification and up-to-108-way OrganNet subtyping, outperforming prior encoders (CTransPath, REMEDIS) especially in low-label regimes.
- Plain-English: A vision model that learned what tissue "looks like" from a hundred million unlabeled image tiles, so it can be adapted to new diagnostic tasks with very few examples.
- Applicability: A4 (canonical frozen-trunk tile encoder — the plug-in substrate a MORPHEUS-style system would consume, not retrain); A1 (its "general-purpose" claim is embedding-generality, NOT task auto-detection — probes are still hand-built per task). Design implication: adopt UNI-class embeddings as a frozen input modality; the promptable-routing layer is exactly what UNI lacks and MORPHEUS must add.
- Novelty implication: Pre-empts any MORPHEUS claim to "general-purpose pathology representation" at the embedding level, but does NOT pre-empt promptable NL task auto-detection (A1) — UNI needs a bespoke linear head per task. Reframes MORPHEUS's contribution as the routing/prompting layer above a UNI-like trunk.

2. **Virchow: A Million-Slide Digital Pathology Foundation Model** — Vorontsov, Bozkurt, Casson, Shaikovski, ... (Paige). arXiv:2309.07778, 2023. https://arxiv.org/abs/2309.07778
- Takeaway: A 632M-param ViT-H DINOv2 encoder trained on 1.5M MSK H&E slides; scaling data+model lifts pan-cancer detection (0.949 specimen AUC across 17 cancers; 0.937 on 7 rare cancers).
- Technical summary: Self-supervised DINOv2 at scale on 1.5M WSIs. Tile features feed a downstream pan-cancer aggregator; the paper is an explicit data/model scaling study showing gains on tile benchmarks and slide-level biomarker prediction, particularly under limited labels.
- Plain-English: Training a bigger image model on a million+ slides makes cancer detection markedly better, including for rare cancers with little data.
- Applicability: A4 (frozen tile encoder feeding slide aggregators — the encode-once/reuse pattern); weak A2 (features are entangled, not per-programme addressable). Design implication: Virchow-scale trunks are commodity inputs; MORPHEUS should treat them as frozen and invest elsewhere.
- Novelty implication: Strengthens the "scale of unlabeled pretraining wins" narrative but says nothing about identifiability or intervention — leaves A2/A5 open for MORPHEUS.

3. **Virchow2: Scaling Self-Supervised Mixed Magnification Models in Pathology** — Zimmermann, Vorontsov, Viret, ... (Paige). arXiv:2408.00738, 2024. https://arxiv.org/abs/2408.00738
- Takeaway: 632M (and 1.85B "Virchow2G") ViT trained on 3.1M multi-magnification WSIs (5x-40x, H&E+IHC) with pathology-specific DINOv2 modifications.
- Technical summary: Extends Virchow to mixed magnifications (2.0/1.0/0.5/0.25 microns per pixel) and multiple stains, introducing domain-tuned augmentation and regularization to the DINOv2 recipe. Achieves SOTA tile features across a broad benchmark suite; demonstrates that magnification/stain diversity in pretraining improves transfer.
- Plain-English: A follow-up that trains on slides at many zoom levels and stain types at once, producing more robust tissue features.
- Applicability: A4 (multi-magnification/multi-stain encoder — informs which raw modalities are worth encoding into the trunk vs supplying as context). Design implication: MORPHEUS's "encode-vs-RAG" decision should note that magnification and stain are cheap to fold INTO the encoder; molecular modalities are the hard case.
- Novelty implication: Pre-empts naive claims that mixing modalities always needs a new architecture — for imaging axes, folding into pretraining suffices. Sharpens MORPHEUS's A4 to be specifically about molecular (proteomics/CNV/SNV) modalities where encoding vs RAG is genuinely contested.

4. **Phikon-v2, A large and public feature extractor for biomarker prediction** — Filiot, Jacob, Mac Kain, Saillard (Owkin). arXiv:2409.09173, 2024. https://arxiv.org/abs/2409.09173
- Takeaway: Public ViT-L DINOv2 on 460M tiles from >100 public cohorts; matches proprietary FMs on 8 slide-level biomarker tasks — but shows top-model margins are statistically non-significant.
- Technical summary: Trained on 460M tiles / >30 cancer sites from public data. On 8 external biomarker-prediction cohorts it ties leading proprietary encoders; a simple ensembling gives +1.75 AUC, and a smaller internal model beats Phikon-v2 on MSI — evidencing diminishing returns from scale alone.
- Plain-English: A freely available model that does about as well as secret commercial ones for predicting molecular biomarkers from H&E, and the differences between all the top models are within noise.
- Applicability: A4 (biomarker/molecular prediction from H&E — the core evidence base for what H&E can and cannot reveal about molecular state); A2 (the "margins are non-significant" finding suggests current trunks encode similar, entangled features — not identifiable programmes). Design implication: MORPHEUS should not expect trunk choice to be decisive; the differentiator is the promptable/identifiable layer, not the encoder.
- Novelty implication: REFRAMES the field's scaling race — if encoders are interchangeable, MORPHEUS's value proposition must be structural (identifiability, prompting, intervention), not "a better trunk." Strong support for MORPHEUS positioning above the encoder.

5. **Scaling Self-Supervised Learning for Histopathology with Masked Image Modeling (Phikon)** — Filiot, Ghermi, Olivier, ... Saillard (Owkin). medRxiv 2023.07.21.23292757, 2023. https://www.medrxiv.org/content/10.1101/2023.07.21.23292757
- Takeaway: The original Phikon — a ViT-B iBOT (masked image modeling) encoder trained on ~40M TCGA tiles; one of the first open pathology FMs and a widely used baseline.
- Technical summary: Uses iBOT/masked-image-modeling self-supervision on TCGA pan-cancer tiles, showing MIM beats prior SSL for histology and enabling weakly-supervised slide classification via attention-based MIL on frozen features.
- Plain-English: An early open model that learned tissue features by masking and reconstructing image patches, becoming a standard comparison point.
- Applicability: A4 (frozen-feature MIL is the dominant consumption pattern MORPHEUS inherits). Design implication: confirms attention-MIL over frozen tiles as the default slide-level readout; MORPHEUS's routing must interoperate with this.
- Novelty implication: Establishes the frozen-trunk + task-head paradigm MORPHEUS explicitly aims to transcend with promptability — supports the "everyone hard-codes probes" gap MORPHEUS targets (A1).

6. **Transformer-based unsupervised contrastive learning for histopathological image classification (CTransPath)** — Wang, Yang, Zhang, ... Han. *Medical Image Analysis*, 2022. https://doi.org/10.1016/j.media.2022.102559
- Takeaway: A CNN-transformer hybrid with semantically-relevant contrastive learning (SRCL) on ~15M pan-cancer tiles — the pre-DINOv2 baseline that most FM papers benchmark against.
- Technical summary: SRCL augments MoCo v3 by pulling in pseudo-positives that are semantically similar across the dataset, mitigating the false-negative problem of instance discrimination on repetitive tissue. Yields transferable features for classification, retrieval, and MIL.
- Plain-English: An earlier self-supervised model that learned to treat different-but-similar-looking tissue patches as related, improving on standard contrastive learning.
- Applicability: A4 (legacy frozen encoder; historical baseline). Design implication: mainly context — establishes the trajectory from contrastive to DINOv2/iBOT trunks.
- Novelty implication: Neutral; useful as the "before" state showing the field never built promptable or identifiable representations.

7. **Hibou: A Family of Foundational Vision Transformers for Pathology** — Nechaev, Pchelnikov, Ivanova (HistAI). arXiv:2406.05074, 2024. https://arxiv.org/abs/2406.05074
- Takeaway: Open Hibou-B/L DINOv2 encoders on >1M proprietary WSIs across diverse stains; Hibou-L leads several patch/slide benchmarks.
- Technical summary: DINOv2 pretraining on a >1M-WSI multi-stain proprietary corpus, open-sourced. Competitive tile/slide performance, reinforcing that the DINOv2-on-large-corpus recipe is now commoditized.
- Plain-English: Another strong, freely released tissue-feature model trained on over a million slides.
- Applicability: A4 (interchangeable frozen trunk option). Design implication: expands the menu of frozen encoders; nothing changes MORPHEUS's above-trunk strategy.
- Novelty implication: Further evidence that trunks are a crowded commodity — strengthens the case for MORPHEUS to differentiate on A1/A2/A5.

8. **RudolfV: A Foundation Model by Pathologists for Pathologists** — Dippel, Feulner, Winterhoff, ... Alber (Aignostics/Charite). arXiv:2401.04079, 2024. https://arxiv.org/abs/2401.04079
- Takeaway: Encoder trained with pathologist-guided data curation across 58 tissue types and 129 stains from 15+ labs, emphasizing domain-knowledge-shaped pretraining.
- Technical summary: Injects expert structure into data selection/curation (stain and tissue diversity, semantic clustering) rather than pure scale; reports SOTA on tumor-microenvironment profiling, biomarker eval, and reference-case search with strong robustness.
- Plain-English: A model built with pathologists deciding what data to feed it, showing curation quality can rival raw scale.
- Applicability: A2 (expert-guided structure is a proxy for imposing addressable semantics — but stops short of identifiable slots); A3 (encodes domain knowledge implicitly). Design implication: MORPHEUS's identifiability goal is a principled version of RudolfV's ad hoc expert curation.
- Novelty implication: Partially pre-empts "domain-knowledge-shaped representation" but only via data curation, not architecturally identifiable programmes — leaves A2's per-programme addressability novel.

9. **H-Optimus-0 (and H-Optimus-1)** — Bioptimus. Open model release, 2024-2025. https://huggingface.co/bioptimus/H-optimus-0
- Takeaway: A ~1.1B-param ViT-g DINOv2 open-weights tile encoder trained on hundreds of millions of tiles from >500k slides; one of the strongest open trunks for biomarker/mutation tasks.
- Technical summary: Released as open weights (model card / GitHub, no full paper for H-Optimus-0) with a giant ViT trained on a large proprietary corpus; positioned for feature extraction in mutation and biomarker prediction pipelines.
- Plain-English: A very large freely available tissue-feature model tuned to help predict genetic and molecular readouts from H&E.
- Applicability: A4 (heavyweight frozen trunk explicitly aimed at molecular prediction). Design implication: strong candidate encoder if MORPHEUS ingests H&E; validates that molecular signal is partially recoverable from morphology.
- Novelty implication: Strengthens the premise that morphology carries molecular information (A4 rationale) but, like all trunks, offers no promptability or intervention — A1/A5 stay open.

---

### Slide-level / whole-slide encoders

10. **A whole-slide foundation model for digital pathology from real-world data (Prov-GigaPath)** — Xu, Usuyama, Bagga, ... Poon, Wang (Microsoft/Providence/UW). *Nature*, 2024. https://www.nature.com/articles/s41586-024-07441-w
- Takeaway: Tile encoder + LongNet-based slide aggregator pretrained on 1.3B tiles / 171k real-world slides (Providence), enabling gigapixel slide-level representations.
- Technical summary: Two-stage GigaPath: DINOv2 tile encoder, then a LongNet dilated-attention slide encoder over tens of thousands of tiles to produce one slide embedding. SOTA on 26 tasks including several mutation/biomarker predictions; real-world (non-TCGA) data is a key ingredient.
- Plain-English: A model that summarizes an entire gigapixel slide into one representation by using an attention method that scales to enormous numbers of tiles.
- Applicability: A4 (slide-level encoder — the object MORPHEUS-style prompting would operate over); A1 (still task-specific heads, no NL routing). Design implication: MORPHEUS can consume slide embeddings as one modality token; the aggregator itself is not promptable.
- Novelty implication: Pre-empts "whole-slide unified representation" at the embedding level; does NOT pre-empt NL-promptable task auto-detection. Reinforces MORPHEUS's A1 gap.

11. **PRISM: A Multi-Modal Generative Foundation Model for Slide-Level Histopathology** — Shaikovski, Casson, Severson, ... Fuchs (Paige). arXiv:2405.10254, 2024. https://arxiv.org/abs/2405.10254
- Takeaway: Slide encoder over Virchow tile embeddings + clinical-report text, giving zero-shot cancer detection/subtyping and report generation from one slide representation.
- Technical summary: Aggregates Virchow tiles into a slide embedding aligned to clinical reports (multi-modal generative pretraining). Zero-shot detection approaches supervised performance; slide-embedding + linear classifier beats supervised MIL aggregators; label-efficient biomarker prediction (10% data beats fully-trained baselines).
- Plain-English: A model that turns a whole slide into a text-aligned summary, letting it detect and describe cancers without task-specific training.
- Applicability: A1 (zero-shot via text alignment is a partial promptability); A3 (report generation is NL<->biology grounding); A4 (text-as-context vs encoded). Design implication: PRISM shows slide<->report alignment yields zero-shot behavior — MORPHEUS's NL prompting should build on this but add explicit task auto-detection and pathway addressability.
- Novelty implication: Partially pre-empts A1/A3 (zero-shot text-conditioned inference exists), but PRISM's "prompts" are free-text retrieval-style, not identified task routing over addressable programmes. Reframes MORPHEUS A1 as principled task-inference beyond CLIP-style zero-shot.

12. **A multimodal whole-slide foundation model for pathology (TITAN)** — Ding, Wagner, Song, ... Mahmood. *Nature Medicine*, 2025 (arXiv:2411.19666, 2024). https://arxiv.org/abs/2411.19666
- Takeaway: Slide encoder pretrained by (1) visual SSL over 335k WSIs then (2) vision-language alignment to reports + 423k synthetic captions — strong zero/few-shot, rare-cancer retrieval, and report generation without clinical labels.
- Technical summary: TITAN (Transformer image+text alignment) uses iBOT-style slide SSL followed by CoCa-style vision-language alignment. Outperforms ROI and slide FMs across linear probe, few-/zero-shot classification, cross-modal retrieval, and report generation; synthetic captions from a pathology copilot expand supervision.
- Plain-English: A whole-slide model that first learns tissue structure, then aligns it to pathology reports, so it can classify, retrieve, and write reports about slides it was never explicitly trained to label.
- Applicability: A1 (zero-shot task transfer via language); A3 (report generation + caption grounding and its evaluation); A4 (reports as encoded-vs-context). Design implication: TITAN is the closest existing "promptable slide model" — MORPHEUS should benchmark against it and target the gaps: no explicit task auto-detection, no pathway-addressable slots, no intervention.
- Novelty implication: The strongest pre-emption risk for MORPHEUS A1/A3 — a slide FM already does language-conditioned zero-shot and report generation. MORPHEUS must clearly differentiate: identifiable/addressable programmes (A2) and interventional queries (A5), which TITAN lacks.

---

### Vision-language pathology models

13. **A visual-language foundation model for computational pathology (CONCH)** — Lu, Chen, Williamson, ... Mahmood. *Nature Medicine*, 2024. https://www.nature.com/articles/s41591-024-02856-4
- Takeaway: CoCa-style image-text model on 1.17M histopathology image-caption pairs; SOTA zero-shot classification, segmentation, captioning, and cross-modal retrieval across 14 benchmarks.
- Technical summary: Contrastive-captioning pretraining on curated educational/PubMed image-caption pairs. Enables text-prompted zero-shot classification and image<->text retrieval; the vision encoder also serves as a strong tile feature extractor.
- Plain-English: A CLIP-like model for pathology that connects tissue images to text, so you can classify or search slides by describing them in words.
- Applicability: A1 (text-prompted zero-shot classification is early promptability); A3 (NL<->morphology grounding); A4 (captions as encoded modality). Design implication: CONCH's text encoder is a natural NL front-end; MORPHEUS's task auto-detection can sit atop CONCH-style embeddings rather than reinvent NL grounding.
- Novelty implication: Pre-empts "connect natural language to pathology" broadly, but zero-shot label matching is not task auto-detection nor pathway addressability. Sharpens MORPHEUS's A1/A2 claims to be about routing and identifiable slots, not CLIP retrieval.

14. **A visual-language foundation model for pathology image analysis using medical Twitter (PLIP)** — Huang, Bianchi, Yuksekgonul, Montine, Zou. *Nature Medicine*, 2023. https://www.nature.com/articles/s41591-023-02504-3
- Takeaway: First pathology vision-language FM — CLIP fine-tuned on 208k image-text pairs (OpenPath, from medical Twitter); strong zero-shot tissue classification.
- Technical summary: Fine-tunes CLIP on OpenPath (Twitter-sourced pathology image-caption pairs), achieving zero-shot F1 0.565-0.832 vs 0.030-0.481 for base CLIP across 4 external datasets, plus image retrieval.
- Plain-English: The first model to teach an image-text AI about pathology using captioned images from doctors' tweets, enabling classify-by-description without labels.
- Applicability: A1 (text-prompted zero-shot); A3 (NL grounding). Design implication: establishes feasibility of NL-driven pathology inference; MORPHEUS extends from "match a label string" to "infer and route the task."
- Novelty implication: Historically pre-empts the vision-language framing; but its promptability is single-shot label matching, leaving A1 task auto-detection, A2 addressability, A5 intervention fully open.

15. **A vision-language foundation model for precision oncology (MUSK)** — Xiang, Zhang, Chen, ... Zou (Stanford). *Nature*, 2025. https://www.nature.com/articles/s41586-024-08378-w
- Takeaway: Unified masked modeling on 50M unpaired pathology images + 1B pathology text tokens, then aligned on 1M image-text pairs; strong across 23 patch/slide tasks incl. molecular biomarkers and outcome prediction.
- Technical summary: MUSK (multimodal transformer, unified masked modeling) leverages large-scale UNPAIRED image and text pretraining before contrastive alignment. Predicts melanoma relapse, pan-cancer prognosis, and immunotherapy response (lung, gastroesophageal) from >8,000 patients, plus VQA, retrieval, and molecular biomarker prediction.
- Plain-English: A model that learns from huge amounts of pathology pictures and text separately, then links them, and does especially well at predicting patient outcomes and treatment response.
- Applicability: A3 (NL<->biology grounding, VQA-style knowledge elicitation); A4 (molecular biomarker + outcome prediction from morphology); A5-adjacent (immunotherapy-response prediction is a treatment-conditioned readout, though correlational not interventional). Design implication: MUSK is the closest to MORPHEUS's "predict response" ambition — but its response prediction is a trained classifier, not a counterfactual drug query. MORPHEUS's A5 is precisely this gap.
- Novelty implication: Strong pre-emption risk for "predict immunotherapy response from a foundation model," but MUSK does it as supervised outcome regression, NOT as an interventional/counterfactual query. MORPHEUS A5 (drug-as-query on a causal geometry) remains distinct and novel.

16. **Pathology Language and Image Pre-training / OpenPath ecosystem — Quilt-1M: One Million Image-Text Pairs for Histopathology** — Ikezogwo, Seyfioglu, Ghezloo, ... Shapiro, Krishna. NeurIPS 2023, arXiv:2306.11207. https://arxiv.org/abs/2306.11207
- Takeaway: Largest open pathology image-text dataset (1M pairs, 802k from educational YouTube via ASR+LLM curation); QuiltNet (CLIP fine-tune) beats SOTA zero-shot/linear across 13 datasets.
- Technical summary: Assembles pairs from histopathology lecture videos (ASR-aligned frames+narration) plus Twitter/papers; releases QuiltNet, which improves zero-shot and linear-probe over PLIP across 8 sub-pathologies and cross-modal retrieval.
- Plain-English: A giant dataset of pathology images paired with expert narration mined from teaching videos, used to train a better image-text model.
- Applicability: A3 (the data substrate for NL<->pathology grounding and its evaluation). Design implication: candidate corpus if MORPHEUS needs pathology-language grounding or an evaluation set for emergent-knowledge probes.
- Novelty implication: Dataset-level enabler; neutral on MORPHEUS claims but underpins the A3 evaluation infrastructure MORPHEUS needs.

17. **Knowledge-enhanced Pretraining for Vision-language Pathology Foundation Model on Cancer Diagnosis (KEEP)** — Zhou, Sun, He, ... Zhang, Sun, Wang, Xie. arXiv:2412.13126, 2024. https://arxiv.org/abs/2412.13126
- Takeaway: Structures millions of image-text pairs with a disease knowledge graph (11,454 diseases, 139,143 attributes) into 143k ontology-aligned groups; improves zero-shot, especially rare subtypes, across 18 benchmarks.
- Technical summary: KEEP (Knowledge-Enhanced Pathology) injects a disease ontology hierarchy into vision-language pretraining so alignment respects semantic structure, yielding gains on rare disease subtyping and 14k+ WSI evaluation.
- Plain-English: A pathology image-text model that is taught the medical hierarchy of diseases, so its language understanding is organized like a real diagnostic ontology.
- Applicability: A2 (ontology hierarchy is a concrete route toward addressable, semantically-organized slots); A3 (structured NL<->biology grounding). Design implication: KEEP's knowledge-graph structuring is a template for MORPHEUS's pathway-addressable slots — impose external biological structure on the representation, don't hope it emerges.
- Novelty implication: Partially pre-empts A2's "impose knowledge structure" idea via disease ontology, but operates on the text/label side, not on identifiable per-programme latent slots. MORPHEUS's identifiability-in-the-representation claim remains distinct.

18. **BiomedCLIP: a multimodal biomedical foundation model pretrained from fifteen million scientific image-text pairs** — Zhang, Xu, Usuyama, ... Poon, Wang, Gao (Microsoft). arXiv:2303.00915, 2023. https://arxiv.org/abs/2303.00915
- Takeaway: General biomedical CLIP on PMC-15M (15M figure-caption pairs, incl. pathology); a broad cross-domain baseline for medical image-text tasks.
- Technical summary: Domain-adapted CLIP (PubMedBERT text encoder, tuned ViT) on 15M PMC pairs across radiology, pathology, microscopy; SOTA on retrieval, classification, VQA vs prior medical CLIPs.
- Plain-English: A general medical image-text model trained on 15 million figures from scientific papers, usable as a baseline for pathology too.
- Applicability: A3 (broad biomedical NL grounding); A4 (cross-modal encoding). Design implication: general-domain baseline; MORPHEUS should show a pathology-specialized promptable model beats a general biomedical CLIP.
- Novelty implication: Neutral baseline; underscores that generic image-text pretraining does not deliver identifiability or intervention.

---

### Molecular-driven / multimodal (H&E + omics)

19. **A pathology foundation model for cancer diagnosis and prognosis prediction (CHIEF)** — Wang, Zhao, Marostica, ... Yu (Harvard). *Nature*, 2024. https://www.nature.com/articles/s41586-024-07894-z
- Takeaway: Weakly-supervised slide framework combining unsupervised tile pretraining + weakly-supervised whole-slide pretraining; predicts diagnosis, molecular profiles, and validated survival across international cohorts.
- Technical summary: CHIEF pretrains tile features then a weakly-supervised slide model with text-anatomy conditioning, generalizing across digitization protocols and populations; predicts genomic/molecular profiles and patient outcomes with cross-cohort validation.
- Plain-English: A model that reads tumor slides to predict both the tumor's molecular makeup and how the patient will fare, tested across countries.
- Applicability: A4 (molecular-profile prediction from H&E — encode-vs-RAG evidence); A5-adjacent (prognosis, but correlational). Design implication: CHIEF shows morphology predicts molecular state cross-population; MORPHEUS can treat H&E as a molecular-informative modality but must add intervention beyond CHIEF's static prognosis.
- Novelty implication: Pre-empts "predict molecular profile + outcome from a pathology FM," but again as supervised prediction. MORPHEUS's counterfactual/interventional framing (A5) is the untouched axis.

20. **Molecular-driven Foundation Model for Oncologic Pathology (THREADS)** — Vaidya, Zhang, Jaume, Song, ... Mahmood. arXiv:2501.16652, 2025. https://arxiv.org/abs/2501.16652
- Takeaway: Slide FM pretrained on 47,171 H&E sections paired with genomic AND transcriptomic profiles — the largest paired morphology-molecular corpus — leading 54 oncology tasks incl. mutation/IHC/treatment-response prediction.
- Technical summary: Multimodal alignment of slide embeddings to paired genomics+transcriptomics during pretraining, so the slide encoder internalizes molecular composition. Outperforms baselines on subtyping, grading, mutation prediction, IHC status, treatment-response, and survival with strong label efficiency and rare-event strength.
- Plain-English: A slide model taught alongside each slide's DNA and RNA data, so it learns the molecular meaning of tissue appearance and predicts genetics and treatment response well.
- Applicability: A4 (the definitive "encode molecular modalities into the trunk during pretraining" datapoint — molecular as encoded, not RAG); A5-adjacent (treatment-response prediction). Design implication: THREADS answers part of MORPHEUS A4 — encoding genomics/transcriptomics INTO pretraining works. MORPHEUS's contribution shifts to (a) proteomics/phospho/CNV/SNV modalities THREADS omits and (b) making the molecular alignment interventional, not just predictive.
- Novelty implication: Strong pre-emption of "align H&E to omics in a foundation model" for DNA/RNA. Reframes MORPHEUS A4 toward the modalities THREADS did NOT encode (proteomics/phospho) and toward the encode-vs-RAG DECISION criterion; A5 remains fully open.

21. **A Multimodal Knowledge-enhanced Whole-slide Pathology Foundation Model (mSTAR)** — Xu, Wang, Zhou, ... Chan, Liang, Wang, Han, Chen. arXiv:2407.15362, 2024. https://arxiv.org/abs/2407.15362
- Takeaway: Injects whole-slide multimodal context (reports + gene expression) into patch representations via "Multimodal Self-TAught pRetraining"; trained on 26k slide-level pairs, evaluated on 97 oncology tasks / 15 cancers.
- Technical summary: mSTAR uses report + gene-expression signals to teach a slide-aware context back into patch features, so downstream MIL inherits multimodal context without needing molecular data at inference.
- Plain-English: A model that uses each slide's report and gene expression during training to make its image features "molecular-aware," even though only the image is needed later.
- Applicability: A4 (molecular+text as training-time encoded context, absent at inference — a concrete "encode not RAG" pattern where the modality is distilled into the trunk). Design implication: supports a MORPHEUS design where scarce molecular modalities are encoded at train time and dropped at inference; contrast with RAG when the modality is available per-query.
- Novelty implication: Pre-empts "distill molecular context into image features." MORPHEUS must articulate WHEN to distill (encode) vs retrieve (RAG) per modality availability/reliability — the decision rule, not the mechanism, is the novel A4 contribution.

22. **CARE: A Molecular-Guided Foundation Model with Adaptive Region Modeling for Whole Slide Image Analysis** — Zhang, Gong, Pang, ... Crispin-Ortuzar, Yu, Li, Gao. CVPR 2026, arXiv:2602.21637. https://arxiv.org/abs/2602.21637
- Takeaway: Two-stage FM — SSL on 34k WSIs, then cross-modal alignment to RNA AND protein profiles to guide morphology-based region construction; leads 33 tasks incl. molecular prediction and survival.
- Technical summary: CARE learns to partition slides into morphologically/molecularly coherent regions guided by RNA+protein alignment, selecting representative ROIs without segmentation labels; molecular guidance improves both region quality and downstream molecular/survival prediction with less pretraining data.
- Plain-English: A model that lets RNA and protein data teach it which tissue regions matter, then focuses on those, improving molecular and survival predictions.
- Applicability: A4 (RNA + PROTEIN alignment — extends THREADS toward proteomics, directly relevant to MORPHEUS's proteomics/phospho question); A2 (molecular-guided regions are a step toward biologically-addressable structure). Design implication: CARE shows protein-level alignment is feasible and improves region identifiability — MORPHEUS's proteomics-encoding case has precedent; the addressable-region idea prefigures pathway-addressable slots.
- Novelty implication: Partially pre-empts "align morphology to proteomics" AND "molecular-guided addressable regions." MORPHEUS's per-PATHWAY (not per-region) identifiability and interventional querying remain the differentiators; monitor CARE as the closest A2/A4 competitor.

23. **Pan-cancer image-based detection of clinically actionable genetic alterations** — Kather, Heij, Grabsch, ... (deep learning from H&E). *Nature Cancer*, 2020. https://www.nature.com/articles/s43018-020-0087-6
- Takeaway: Seminal demonstration that mutations, MSI, and molecular subtypes are predictable directly from H&E across many cancers with deep learning — the empirical basis for "morphology encodes molecular state."
- Technical summary: Trains CNNs on TCGA H&E to predict point mutations, MSI, and molecular subtypes pan-cancer, establishing which alterations are morphology-visible and setting the biomarker-from-H&E research program.
- Plain-English: Early proof that a computer can guess a tumor's key genetic mutations just from a stained tissue image.
- Applicability: A4 (defines the ceiling/floor of what H&E alone reveals molecularly — informs when a molecular modality MUST be supplied vs can be inferred). Design implication: MORPHEUS's encode-vs-RAG rule should be calibrated by which alterations are morphology-recoverable (encode/infer) vs morphology-invisible (must RAG/supply).
- Novelty implication: Grounds the A4 premise; neutral on MORPHEUS's structural claims but sets the predictability prior MORPHEUS should cite.

---

### Slide aggregation architecture (context)

24. **Scaling Vision Transformers to Gigapixel Images via Hierarchical Self-Supervised Learning (HIPT)** — Chen, Chen, Li, Chen, Trister, Krishnan, Mahmood. CVPR 2022, arXiv:2206.02647. https://arxiv.org/abs/2206.02647
- Takeaway: Two-level hierarchical ViT (16x16 -> 256 -> 4096 px) with self-supervised pretraining on 10.6k WSIs; early principled whole-slide representation learning.
- Technical summary: Exploits WSI hierarchy with nested DINO pretraining at cell and tissue-microenvironment scales, aggregating to slide-level for subtyping and survival; a precursor to modern slide encoders.
- Plain-English: A model that builds up understanding of a slide from cells to tissue regions to the whole slide, in stages.
- Applicability: A4 (hierarchical aggregation architecture MORPHEUS may reuse). Design implication: hierarchy is the standard way to reach slide level; MORPHEUS's promptable layer sits above whichever aggregator is used.
- Novelty implication: Architectural context; neutral on promptability/identifiability.

---

### Conversational / generative pathology (NL grounding)

25. **A Foundational Multimodal Vision Language AI Assistant for Human Pathology (PathChat)** — Lu, Chen, Williamson, Chen, ... Parwani, Mahmood. *Nature*, 2024 (arXiv:2312.07814, 2023). https://arxiv.org/abs/2312.07814
- Takeaway: Vision encoder (CONCH-lineage, 100M images) + LLM, instruction-tuned on 250k+ pathology visual-language instructions; 87% MCQ accuracy with clinical context, beats GPT-4V.
- Technical summary: Couples a pathology vision encoder to an LLM and fine-tunes on diverse disease-agnostic instructions, yielding an interactive assistant that answers open-ended and multiple-choice pathology questions and is pathologist-preferred over GPT-4V/LLaVA-Med.
- Plain-English: A ChatGPT-like assistant for pathology that looks at a tissue image and answers questions about it in natural language.
- Applicability: A1 (interactive NL interface — closest existing thing to "describe the task in language"); A3 (NL<->biology grounding + emergent-knowledge elicitation via QA, and its MCQ evaluation). Design implication: PathChat is the NL front-end paradigm MORPHEUS's task auto-detection resembles; but it answers about an image, it does not ROUTE to identified quantitative programme readouts or run interventions.
- Novelty implication: Pre-emption risk for A1/A3 "talk to your pathology model." MORPHEUS must differentiate: PathChat generates language ABOUT a slide; MORPHEUS infers a task and routes to addressable, quantitative, interventionable programme slots. The distinction (conversation vs identified-slot routing + intervention) is the defensible novelty.

---

### Benchmarks, evaluation, and critical analysis (informs A2/A3 evaluation)

26. **Towards A Generalizable Pathology Foundation Model via Unified Knowledge Distillation (GPFM + benchmark)** — Ma, Guo, Zhou, Wang, ... Chen. arXiv:2407.18449, 2024/2025. https://arxiv.org/abs/2407.18449
- Takeaway: Distills multiple expert FMs into one GPFM (190M images / ~72k slides / 34 tissues) AND contributes a 72-task benchmark across 6 clinical categories incl. survival and VQA.
- Technical summary: Unified expert+self knowledge distillation produces a single strong encoder; the paper's lasting value is a broad, standardized 72-task evaluation exposing where FMs still fail.
- Plain-English: Combines several existing pathology models into one and tests all of them on 72 tasks to see who really generalizes.
- Applicability: A3 (evaluation infrastructure — how to MEASURE generality/emergent knowledge); A1 (still per-task heads). Design implication: adopt GPFM's task suite as a MORPHEUS evaluation harness; extend it with promptability and intervention benchmarks it lacks.
- Novelty implication: Highlights that "general" is measured only by task coverage, never by task auto-detection or intervention — supports MORPHEUS's claim that A1/A5 are unmeasured gaps.

27. **Towards Large-Scale Training of Pathology Foundation Models (kaiko + eva)** — kaiko.ai (Aben, de Jong, Gatopoulos, ... Tang). arXiv:2404.15217, 2024. https://arxiv.org/abs/2404.15217
- Takeaway: Open TCGA-trained FMs plus the open-source "eva" evaluation framework standardizing how pathology FMs are compared.
- Technical summary: Scalable SSL pipeline on public WSIs with strong results (breast subtyping, colorectal nuclei segmentation), and — most usefully — a released, reproducible evaluation harness (eva) for frozen-feature benchmarking.
- Plain-English: Open models trained on public slides, plus a standard scoreboard tool so everyone tests pathology models the same way.
- Applicability: A3/A2 (standardized frozen-feature evaluation — the measurement substrate). Design implication: use eva for reproducible frozen-trunk comparisons when selecting MORPHEUS's encoder inputs.
- Novelty implication: Neutral tooling; reinforces that current evaluation measures embedding quality, not promptability/identifiability/intervention — MORPHEUS must define new metrics.

28. **Comparing Computational Pathology Foundation Models using Representational Similarity Analysis** — Mishra, Lotter. arXiv:2509.15482, 2025. https://arxiv.org/abs/2509.15482
- Takeaway: RSA across 6 FMs finds strong sensitivity to SLIDE-SPECIFIC (batch) features but WEAK disease-related structure; stain normalization cuts slide-dependence 5.5-20.5%.
- Technical summary: Applies representational similarity analysis; UNI2/Virchow2 are most distinct, Prov-GigaPath most average; vision-language models encode more compact representations than vision-only; shared training paradigms do not yield shared representations. Crucially, disease signal is under-represented relative to site/batch signal.
- Plain-English: A neuroscience-style analysis showing these pathology models mostly encode which slide/scanner an image came from, and only weakly encode the actual disease.
- Applicability: A2 (direct evidence AGAINST current identifiability — features are entangled with batch, not organized by biological programme); A4 (batch confound informs encode-vs-RAG). Design implication: this is the empirical motivation for MORPHEUS A2 — if trunks encode site more than disease, identifiable pathway-addressable slots are needed, not assumed.
- Novelty implication: STRENGTHENS MORPHEUS's core A2 premise: existing FMs are NOT identifiable/disease-organized. Use as primary citation that per-programme addressability is an open, needed problem.

29. **HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis** — Jaume, Doucet, Song, Lu, ... Kim, Mahmood. NeurIPS 2024 Spotlight, arXiv:2406.16192. https://arxiv.org/abs/2406.16192
- Takeaway: 1,229 spatial-transcriptomics profiles paired with H&E (153 cohorts, 26 organs; 2.1M expression-morphology pairs), with a benchmark for predicting gene expression from morphology.
- Technical summary: Releases paired ST+H&E data and HEST-Library; benchmarks FMs on gene-expression prediction from tiles, biomarker discovery, and multimodal learning — quantifying how much transcriptomic signal is morphology-recoverable per gene/organ.
- Plain-English: A large dataset linking tissue images to spatially-resolved gene activity, used to test how well models can read gene expression off the image.
- Applicability: A4 (the empirical map of which genes/pathways are morphology-encodable vs not — the exact evidence needed for MORPHEUS's encode-vs-RAG decision at the pathway/gene level); A3 (grounding morphology to molecular readouts). Design implication: use HEST-1k to calibrate MORPHEUS's per-pathway encode-vs-RAG rule and to evaluate emergent molecular knowledge.
- Novelty implication: Provides the measurement substrate for A4's decision criterion; strengthens the feasibility of pathway-level molecular grounding while showing many genes are weakly predictable (justifying RAG for those).

30. **A Survey on Computational Pathology Foundation Models: Datasets, Adaptation Strategies, and Evaluation Tasks** — Li, Wan, Wu, ... Sorger, Semenov, Zhao. arXiv:2501.15724, 2025. https://arxiv.org/abs/2501.15724
- Takeaway: Systematic review of CPathFMs (uni- vs multi-modal), adaptation strategies, and evaluation, explicitly flagging the absence of standardized benchmarks and limited data access as field-wide gaps.
- Technical summary: Taxonomizes pretraining data, SSL/contrastive/multimodal strategies, and downstream evaluation (segmentation, classification, biomarker discovery); identifies open problems in reproducibility and generalization.
- Plain-English: A map of the whole pathology-foundation-model field: what data and methods exist, how they're tested, and what's still missing.
- Applicability: A1-A5 (landscape orientation). Design implication: positions MORPHEUS relative to the surveyed axis of "adaptation strategies" — none of which are NL task auto-detection, identifiable slots, or intervention. Design context for the whole rebase.
- Novelty implication: Confirms via a neutral third-party survey that the field frames advances as data/adaptation/evaluation — NOT as promptability (A1), identifiability (A2), emergent-knowledge elicitation (A3), encode-vs-RAG decisions (A4), or intervention (A5). Strongest external evidence that MORPHEUS's five axes are genuinely under-addressed.
