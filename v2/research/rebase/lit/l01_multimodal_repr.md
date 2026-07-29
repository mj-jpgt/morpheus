# Multimodal representation learning & fusion

Lane id: `l01_multimodal_repr`. Remit: contrastive/aligned multimodal representation learning; fusion architectures (early/late/attention/PoE-MoE); modality dropout; cross-modal alignment (CLIP/ImageBind/Perceiver-IO lineage). Foundations a promptable multimodal biomedical representation builds on.

MORPHEUS rebase axes referenced below: **A1** promptable unified rep + NL task auto-detection; **A2** identified, pathway-addressable slots; **A3** NL<->biology grounding + emergent-knowledge elicitation & its evaluation; **A4** multimodal prompting — when to ENCODE a modality vs treat as RAG context, frozen-trunk plug-in; **A5** interventional/causal-geometry queries.

---

### 1. Learning Transferable Visual Models From Natural Language Supervision (CLIP)
Radford, Kim, Hallacy, Ramesh, Goh, Agarwal, Sastry, Askell, Mishkin, Clark, Krueger, Sutskever. ICML 2021. arXiv:2103.00020.
- **Takeaway:** A dual-encoder trained with a symmetric InfoNCE contrastive loss on 400M web image-text pairs yields a joint embedding space that supports zero-shot transfer by phrasing tasks as natural-language prompts.
- **Technical summary:** Separate image and text encoders are projected to a shared space; a temperature-scaled contrastive objective pulls matched pairs together and pushes mismatched pairs apart across the batch. Zero-shot classification is done by embedding class-name prompts ("a photo of a {label}") and taking nearest text — a ResNet-50x64/ViT variant matches supervised ImageNet ResNet-50 with no fine-tuning. Prompt engineering and ensembling of templates materially move accuracy.
- **Plain-English:** Instead of training on fixed labels, it learns which pictures and which sentences go together, so you can later "ask" it about a new category just by describing it in words.
- **Applicability:** A1 (the founding demonstration that a task can be *specified in natural language at inference* against a frozen representation rather than hard-coded as a probe head — the direct ancestor of MORPHEUS's promptable rep), A3 (text is the grounding interface to a learned space). Design implication: MORPHEUS's NL task auto-detection is a generalization of CLIP prompt-based zero-shot to biological programmes.
- **Novelty implication:** Pre-empts any bare claim of "prompt a representation in natural language." MORPHEUS must differentiate on *task auto-detection/routing* and *pathway-addressable structure*, not on prompting per se, which CLIP owns.

### 2. Scaling Up Visual and Vision-Language Representation Learning With Noisy Text Supervision (ALIGN)
Jia, Yang, Xia, Chen, Parekh, Pham, Le, Sung, Li, Duerig. ICML 2021. arXiv:2102.05918.
- **Takeaway:** A simple dual-encoder contrastive model trained on 1.8B *unfiltered* noisy alt-text pairs matches or beats curated-data models, showing scale can substitute for curation.
- **Technical summary:** EfficientNet image encoder + BERT text encoder aligned by normalized softmax contrastive loss; no dataset cleaning beyond frequency-based filtering. Sets SOTA on Flickr30K/MSCOCO retrieval and strong zero-shot ImageNet, outperforming cross-attention models at retrieval.
- **Plain-English:** If you have enough messy image-caption pairs from the web, you don't need to clean them — sheer volume teaches a good shared space.
- **Applicability:** A4 (encode-vs-context: evidence that noisy paired supervision at scale is enough to *encode* a modality into a shared trunk), A3. Implication: for MORPHEUS, modalities with abundant but noisy NL pairings (e.g., gene descriptions) may be worth encoding rather than RAG-ing.
- **Novelty implication:** Reframes — supports the position that encoding a modality only pays off at data scale; MORPHEUS's "when to encode vs RAG" decision rule is the novel contribution ALIGN does not address.

### 3. Perceiver: General Perception with Iterative Attention
Jaegle, Gimenez, Brock, Zisserman, Vinyals, Carreira. ICML 2021. arXiv:2103.03206.
- **Takeaway:** A fixed-size latent bottleneck cross-attends into arbitrarily large/heterogeneous input arrays, decoupling compute from input size and modality.
- **Technical summary:** Inputs (pixels, audio, point clouds) are projected via cross-attention onto a small set of latent vectors, then processed by a latent self-attention tower; the same architecture handles images, audio, video, and point clouds with only positional/Fourier feature changes. Scales linearly, not quadratically, with input size.
- **Plain-English:** Rather than let every input token talk to every other (too expensive), it funnels any kind of raw data through a small shared "workspace" of latent slots.
- **Applicability:** A2 (a small set of latent slots is the architectural seed of "identified slots"), A4 (a single modality-agnostic trunk that new modalities plug into). Implication: MORPHEUS's addressable per-programme slots can be realized as identified Perceiver-style latents.
- **Novelty implication:** Pre-empts "latent slots for multimodal fusion" as architecture; MORPHEUS must claim *identifiability/addressability* of those slots (which Perceiver latents lack — they are anonymous and permutation-free), not their existence.

### 4. Perceiver IO: A General Architecture for Structured Inputs & Outputs
Jaegle, Borgeaud, Alayrac, Doersch, Ionescu, Ding, et al. ICLR 2022. arXiv:2107.14795.
- **Takeaway:** Adds a query-based decoder to Perceiver so arbitrary structured outputs are produced by learned output queries, giving one architecture for language, vision, multimodal, and RL.
- **Technical summary:** Read-process-write: encode inputs to latents via cross-attention, refine with latent self-attention, then decode by cross-attending output queries (whose semantics define each output) against the latents. Matches BERT on GLUE with no tokenizer, does optical flow, StarCraft II, and multimodal autoencoding under one design.
- **Plain-English:** A general model that reads any input and, by asking the right "output questions," produces any shaped answer — words, pixels, or actions.
- **Applicability:** A1 (output queries are a mechanism to *route* a specified task through a shared trunk — a template for MORPHEUS task routing), A2, A4. Implication: MORPHEUS task auto-detection could compile an NL request into a Perceiver-IO-style output query set.
- **Novelty implication:** Strengthens A1 feasibility but also pre-empts "query the representation for a task"; MORPHEUS's differentiator is *inferring* the query from NL, not being handed a query.

### 5. Flamingo: a Visual Language Model for Few-Shot Learning
Alayrac, Donahue, Luc, Miech, Barr, Hasson, et al. NeurIPS 2022. arXiv:2204.14198.
- **Takeaway:** Bridges a frozen vision encoder and a frozen LLM with trainable gated cross-attention + a Perceiver Resampler, enabling in-context few-shot multimodal tasks from interleaved image-text prompts.
- **Technical summary:** A Perceiver Resampler maps variable visual features to fixed tokens; gated cross-attention-dense layers interleaved into a frozen Chinchilla LM inject vision while preserving language priors. Trained on interleaved web corpora; outperforms task-specific fine-tuned models on many benchmarks with only a handful of examples.
- **Plain-English:** Keep a strong language model and a strong vision model frozen, glue them with a small trainable adapter, and the system learns new visual tasks from a few examples shown in the prompt.
- **Applicability:** A4 (canonical *frozen-trunk plug-in*: freeze powerful unimodal models, train only the bridge — the exact pattern MORPHEUS proposes for plugging in proteomics/CNV), A1 (in-context task specification). Implication: adopt gated cross-attention adapters so encoded biomedical modalities attach to a frozen NL trunk without catastrophic forgetting.
- **Novelty implication:** Pre-empts "frozen-trunk multimodal plug-in via adapters." MORPHEUS must claim the biological/identifiability and encode-vs-RAG decision layers, not the frozen-adapter mechanism.

### 6. CoCa: Contrastive Captioners are Image-Text Foundation Models
Yu, Wang, Vasudevan, Yeung, Seyedhosseini, Wu. TMLR 2022. arXiv:2205.01917.
- **Takeaway:** Jointly trains a contrastive (alignment) and a captioning (generative) objective in one encoder-decoder, getting both a retrieval-quality embedding and a generative head.
- **Technical summary:** A split decoder — unimodal layers produce a contrastive text embedding, cross-attention layers produce captions — lets one forward pass serve both objectives; yields SOTA on ImageNet, retrieval, VQA, and captioning with a single pretrained model.
- **Plain-English:** One model learns both to *match* images and text and to *describe* images, so it is good at retrieval and generation at once.
- **Applicability:** A1/A3 (unifying an alignment objective with a generative NL head is how MORPHEUS could both align modalities and *emit* NL answers), A4. Implication: a MORPHEUS trunk should carry both a contrastive alignment head and a generative NL head to support prompting and explanation.
- **Novelty implication:** Reframes MORPHEUS's "grounding + generation" as an established recipe; the novelty must be in *emergent-knowledge elicitation and its evaluation* (A3), not in combining contrastive+generative losses.

### 7. FLAVA: A Foundational Language And Vision Alignment Model
Singh, Hu, Goswami, Couairon, Galuba, Rohrbach, Kiela. CVPR 2022. arXiv:2112.04482.
- **Takeaway:** A single model with unimodal image, unimodal text, and a multimodal fusion encoder trained jointly on contrastive, masked, and multimodal objectives to be good at vision, language, AND cross-modal tasks.
- **Technical summary:** ViT image encoder + text encoder feed a transformer multimodal encoder; losses combine global contrastive, masked-image/masked-language modeling, and image-text matching. Evaluated across 35 vision, language, and multimodal tasks from public data only.
- **Plain-English:** One foundation model built to handle pictures alone, text alone, and the two together, rather than specializing in just one.
- **Applicability:** A4 (explicit unimodal + fusion branches — informs when MORPHEUS keeps a modality in its own encoder vs a shared fusion encoder), A1. Implication: a hybrid design with per-modality encoders plus a shared fusion trunk supports selective encoding.
- **Novelty implication:** Pre-empts "one model for unimodal and multimodal tasks." MORPHEUS differentiates on biological addressability and interventional queries, which FLAVA does not attempt.

### 8. BLIP: Bootstrapping Language-Image Pre-training for Unified Vision-Language Understanding and Generation
Li, Li, Xiong, Hoi. ICML 2022. arXiv:2201.12086.
- **Takeaway:** Unifies understanding and generation with a multimodal mixture-of-encoder-decoder plus a caption-and-filter bootstrapping step that cleans noisy web captions.
- **Technical summary:** A shared transformer operates in three modes (unimodal, image-grounded encoder, image-grounded decoder) trained with contrastive, matching, and language-modeling losses; a captioner generates synthetic captions and a filter removes noisy ones to self-improve the training set.
- **Plain-English:** The model writes better captions for its own training images and throws out bad ones, then retrains on the improved data.
- **Applicability:** A3 (bootstrapping/self-filtering of NL annotations is a template for cleaning noisy biological text used for grounding), A4. Implication: MORPHEUS could bootstrap NL descriptions of programmes and filter them by consistency.
- **Novelty implication:** Neutral; mainly informs data pipeline. No direct MORPHEUS claim collision.

### 9. BLIP-2: Bootstrapping Language-Image Pre-training with Frozen Image Encoders and Large Language Models
Li, Li, Savarese, Hoi. ICML 2023. arXiv:2301.12597.
- **Takeaway:** A lightweight Querying Transformer (Q-Former) with learned query tokens bridges a *frozen* image encoder and a *frozen* LLM in two stages, beating Flamingo-80B on zero-shot VQA with ~54x fewer trainable params.
- **Technical summary:** Stage 1 trains Q-Former against the frozen vision encoder with contrastive/matching/generation losses to extract the most text-relevant visual tokens; stage 2 feeds those query outputs as soft prompts into a frozen LLM. Only the Q-Former is trained.
- **Plain-English:** A small translator module learns to hand a frozen language model exactly the visual details it needs, so almost nothing big has to be retrained.
- **Applicability:** A4 (the strongest evidence for MORPHEUS's frozen-trunk-plus-tiny-bridge thesis; Q-Former's learned queries = a candidate mechanism for encoding a biomedical modality into a few soft-prompt tokens), A1. Implication: encode proteomics/phospho as a handful of Q-Former-style query tokens injected into a frozen NL trunk.
- **Novelty implication:** Strongly pre-empts "encode a modality as soft-prompt tokens for a frozen LLM." MORPHEUS's decision rule for *which* modalities merit a Q-Former vs RAG, and pathway-addressability of the query tokens, is the defensible novelty.

### 10. ImageBind: One Embedding Space To Bind Them All
Girdhar, El-Nouby, Liu, Singh, Alwala, Joulin, Misra. CVPR 2023. arXiv:2305.05665.
- **Takeaway:** Binds six modalities (image, text, audio, depth, thermal, IMU) into one space using *only* image-paired data per modality, yielding emergent cross-modal alignment between never-co-observed pairs.
- **Technical summary:** Each modality encoder is contrastively aligned to a frozen image/text (CLIP-style) space using naturally occurring image-modality pairs; alignment is transitive, so audio-text retrieval emerges without audio-text training pairs. Enables cross-modal retrieval, arithmetic, and generation zero-shot.
- **Plain-English:** By anchoring every sensor type to images, all the sensors end up sharing one space — and can be compared even for combinations never seen together in training.
- **Applicability:** A4 (blueprint for adding N biomedical modalities by aligning each to a shared anchor without all-pairs data — critical since proteomics-phospho-CNV paired data is scarce), A2. Implication: MORPHEUS can bind sparse modalities to a common NL/expression anchor and get emergent cross-modal transfer.
- **Novelty implication:** Pre-empts "bind many modalities via a common anchor with emergent cross-modal alignment." MORPHEUS must claim *biological* emergent knowledge and its evaluation (A3), and interventional geometry (A5), not the binding recipe.

### 11. LanguageBind: Extending Video-Language Pretraining to N-modality by Language-based Semantic Alignment
Zhu, Lin, Ning, Yan, Cui, Wang, et al. ICLR 2024. arXiv:2310.01852.
- **Takeaway:** Uses *language* (not image) as the binding anchor across video, infrared, depth, audio, so the shared space is directly NL-addressable.
- **Technical summary:** Freezes a CLIP-style language encoder and contrastively aligns each modality encoder to it via modality-language pairs, with a large multimodal dataset (VIDAL); language-centered alignment improves zero-shot performance over image-centered binding on several modalities.
- **Plain-English:** Instead of anchoring everything to images, anchor everything to text — making the joint space something you can query in words directly.
- **Applicability:** A1/A3 (language-as-anchor is exactly MORPHEUS's premise that the shared space should be promptable in NL), A4. Implication: prefer an NL-anchored binding over an expression-anchored one so biological queries land natively in the space.
- **Novelty implication:** Pre-empts "language-anchored N-modality binding." Reinforces that MORPHEUS's edge is task auto-detection + pathway slots + causal queries, not the NL-anchor idea.

### 12. Meta-Transformer: A Unified Framework for Multimodal Learning
Zhang, Gong, Zhang, Li, Qiao, Ouyang, Yue. arXiv:2307.10802 (2023).
- **Takeaway:** A single *frozen* transformer encoder (pretrained on images) processes 12 modalities via per-modality tokenizers, showing a shared backbone can generalize across radically different data with no modality-specific backbone.
- **Technical summary:** Data-to-sequence tokenizers map each modality (text, image, point cloud, audio, graph, hyperspectral, IMU, etc.) into a common token space; a frozen shared encoder plus lightweight task heads handle all 12. Competitive across many benchmarks without training the backbone.
- **Plain-English:** One frozen brain, twelve different "input adapters" — the same core network handles everything from graphs to audio.
- **Applicability:** A4 (extreme frozen-trunk plug-in: only tokenizer + head are trained per modality — directly informs MORPHEUS's frozen-trunk plug-in cost model), A2. Implication: encoding a new biomedical modality may cost only a tokenizer + head against a frozen MORPHEUS trunk.
- **Novelty implication:** Pre-empts "one frozen backbone for many modalities via tokenizers." MORPHEUS's contribution is the biology-specific decision of what to encode and identifiable slots, not the frozen-shared-encoder claim.

### 13. ONE-PEACE: Exploring One General Representation Model Toward Unlimited Modalities
Wang, Wang, Lin, Bai, Zhou, Zhou, Wang, Zhou. arXiv:2305.11172 (2023).
- **Takeaway:** A scalable architecture with shared self-attention and modality-specific adapters/FFNs, extensible to new modalities without retraining existing ones, aligned by cross-modal contrastive + intra-modal denoising.
- **Technical summary:** Modality adapters + a shared transformer with modality-specific feed-forward "experts"; pretraining combines cross-modal contrastive and masked/denoising intra-modal objectives across vision, audio, language. Designed so a new modality adds an adapter+FFN branch, preserving prior alignment.
- **Plain-English:** A model built so you can keep bolting on new senses without breaking the ones it already learned.
- **Applicability:** A4 (modality-specific FFN experts + shared attention is a concrete design for MORPHEUS to add proteomics/phospho branches incrementally — MoE-style per-modality routing), A2. Implication: use per-modality expert FFNs so added modalities do not disturb existing programme slots.
- **Novelty implication:** Pre-empts "extensible per-modality experts on a shared trunk." MORPHEUS must claim identifiable/pathway-addressable experts and causal queries.

### 14. Attention Bottlenecks for Multimodal Fusion (MBT)
Nagrani, Yang, Arnab, Jansen, Schmid, Sun. NeurIPS 2021. arXiv:2107.00135.
- **Takeaway:** Forcing cross-modal information through a few "fusion bottleneck" tokens (and only in later layers) improves accuracy over free all-to-all cross-attention while cutting compute ~50%.
- **Technical summary:** Modalities attend freely within themselves but exchange information only through a small shared set of bottleneck latents; restricting fusion to later layers lets early layers specialize. SOTA on audio-visual video classification at lower FLOPs.
- **Plain-English:** Don't let two data streams gossip freely — make them pass only the essentials through a narrow shared channel, which is both cheaper and better.
- **Applicability:** A2 (bottleneck tokens are candidate *identified fusion slots*; late-fusion specialization argues for keeping modality-specific programme structure before mixing), A4. Implication: MORPHEUS should fuse modalities through a small, addressable bottleneck rather than dense cross-attention, aiding both identifiability and cost.
- **Novelty implication:** Strengthens A2 feasibility (bottleneck = slot substrate) but pre-empts "fusion through shared latent bottleneck"; MORPHEUS's claim is that the bottleneck slots are *pathway-identified and promptable*.

### 15. LiT: Zero-Shot Transfer with Locked-image Tuning
Zhai, Wang, Mustafa, Steiner, Keysers, Kolesnikov, Beyer. CVPR 2022. arXiv:2111.07991.
- **Takeaway:** Contrastively tuning only the *text* tower while keeping a strong pretrained image tower *locked* dramatically improves zero-shot transfer over training both from scratch.
- **Technical summary:** A pre-trained (even supervised) image encoder is frozen; only the text encoder is contrastively trained to align to it. This "locked-image, tuned-text" setup beats fully-trained and fully-fine-tuned CLIP-style alternatives on zero-shot ImageNet and retrieval.
- **Plain-English:** Freeze the vision half, only teach the text half to speak its language — you get better zero-shot recognition than training both together.
- **Applicability:** A4 (direct evidence that MORPHEUS can *freeze a strong modality encoder* — e.g., an existing single-cell foundation trunk — and only train the NL/alignment side), A1. Implication: lock the biological trunk, tune the NL interface, to preserve learned biology while gaining promptability.
- **Novelty implication:** Pre-empts "freeze the representation, tune only the language interface." MORPHEUS's defensible novelty is the biological addressability and causal query layer atop the locked trunk.

### 16. Multimodal Few-Shot Learning with Frozen Language Models (Frozen)
Tsimpoukelli, Menick, Cabi, Eslami, Vinyals, Hill. NeurIPS 2021. arXiv:2106.13884.
- **Takeaway:** Train only a vision encoder to emit soft prompts for a *frozen* LLM; the system inherits the LLM's few-shot, in-context abilities for multimodal tasks.
- **Technical summary:** A vision encoder maps images to a sequence of embeddings prepended to a frozen autoregressive LM; trained only on image captioning, it generalizes to VQA and few-shot visual tasks by in-context prompting, transferring the LM's rapid-learning behavior to vision.
- **Plain-English:** Teach a camera to whisper to a language model in the model's own vocabulary, and the model can suddenly answer questions about pictures using its existing few-shot talent.
- **Applicability:** A4 (the seminal "encode a modality as a soft prompt for a frozen LM" result — the mechanistic core of MORPHEUS's encode-as-prompt option), A1. Implication: encoded biomedical modalities can be injected as prefix soft prompts so the frozen NL trunk's reasoning transfers.
- **Novelty implication:** Pre-empts the *mechanism* of encode-as-soft-prompt. MORPHEUS must own the *policy* (encode vs RAG) and the biological grounding/eval.

### 17. Are Multimodal Transformers Robust to Missing Modality?
Ma, Ren, Zhao, Testuggine, Peng. CVPR 2022. arXiv:2204.05454.
- **Takeaway:** Multimodal transformers degrade sharply when a modality is missing at test time, and the optimal fusion strategy is *dataset-dependent*; a searched/robust fusion and multi-task modality-dropout training improve robustness.
- **Technical summary:** Systematically shows accuracy collapse under missing modalities, finds no single fusion (early/late/mid) is universally best, and proposes automatic fusion-architecture search plus training that optimizes for missing-modality cases.
- **Plain-English:** These models often lean too hard on one input; take it away and they fail — so you must train and design deliberately for the case where an input is absent.
- **Applicability:** A4 (biomedical modalities are chronically missing — phospho/CNV present for some patients only; this directly motivates modality-dropout training and a fusion that tolerates absence — a core MORPHEUS robustness requirement), A2. Implication: MORPHEUS must train with modality dropout and validate graceful degradation, and the encode-vs-RAG choice should factor missingness.
- **Novelty implication:** Reframes — supplies the problem statement MORPHEUS's "when to encode vs RAG" answers; strengthens the case that a plug-in/optional-modality design is necessary, not just convenient.

### 18. Multimodal Generative Models for Scalable Weakly-Supervised Learning (MVAE, Product-of-Experts)
Wu, Goodman. NeurIPS 2018. arXiv:1802.05335.
- **Takeaway:** Models the joint multimodal posterior as a *product of experts* over per-modality posteriors, so any subset of modalities can be encoded and missing modalities are handled naturally.
- **Technical summary:** Each modality has its own inference network; the joint posterior is their product (with a prior expert), giving a principled way to marginalize absent modalities and train with sub-sampled modality subsets. Enables cross-modal generation and weak supervision.
- **Plain-English:** Combine each sensor's "opinion" by multiplying them, which lets the model work even when some sensors are switched off.
- **Applicability:** A4 (the PoE formulation is the classical answer to *variable/missing modality fusion* — foundational for MORPHEUS's optional-modality encoding), A2 (per-modality experts = addressable contributions). Implication: a PoE combination rule lets MORPHEUS accept any subset of {RNA, protein, phospho, CNV} and infer the rest.
- **Novelty implication:** Pre-empts PoE fusion as a mechanism. MORPHEUS should cite as the fusion primitive and claim novelty in identifiability + prompting + causal use, not in PoE itself.

### 19. Variational Mixture-of-Experts Autoencoders for Multi-Modal Deep Generative Models (MMVAE)
Shi, Siddharth, Paige, Torr. NeurIPS 2019. arXiv:1911.03393.
- **Takeaway:** Replaces PoE with a *mixture of experts* joint posterior to better satisfy four desiderata (coherent joint/cross generation, latent factorization, synergy), fixing PoE's tendency to be dominated by overconfident modalities.
- **Technical summary:** Joint posterior = average of per-modality posteriors; a stratified/IWAE objective yields both coherent cross-modal generation and informative latents. Analyzes when MoE beats PoE on generative coherence.
- **Plain-English:** Instead of multiplying sensor opinions (which lets a loud one dominate), average them, giving fairer, more coherent cross-modal generation.
- **Applicability:** A4 (PoE-vs-MoE tradeoff is exactly the fusion-rule choice MORPHEUS must make per modality set), A2. Implication: MORPHEUS's fusion rule should be chosen (or learned) rather than defaulted; MoE guards against a dominant modality drowning weaker programmes.
- **Novelty implication:** Neutral/foundational — frames the PoE/MoE design axis MORPHEUS operates within.

### 20. Generalized Multimodal ELBO (MoPoE-VAE)
Sutter, Daunhawer, Vogt. ICLR 2021. arXiv:2105.02470.
- **Takeaway:** A mixture-of-products-of-experts objective evaluates all modality *subsets* in the powerset, unifying PoE and MoE and improving both joint and conditional generation without extra training objectives.
- **Technical summary:** Defines subset posteriors as products of the experts in each subset, and the joint posterior as their mixture over the powerset; a single ELBO covers all self- and cross-modal terms, giving better scalability across modality subsets than PoE or MoE alone.
- **Plain-English:** Consider every possible combination of available sensors and blend them, getting the best of both the "multiply" and "average" strategies.
- **Applicability:** A4 (principled handling of *any subset* of biomedical modalities — the general case MORPHEUS faces with irregular multi-omic coverage), A2 (subset-addressable structure). Implication: a MoPoE-style objective lets MORPHEUS train once and serve any observed modality subset.
- **Novelty implication:** Pre-empts subset-generalizing fusion. MORPHEUS differentiates via NL-promptability and interventional queries layered on top.

### 21. Contrastive Multiview Coding (CMC)
Tian, Krishnan, Isola. ECCV 2020. arXiv:1906.05849.
- **Takeaway:** Treats multiple views/modalities of a scene as different "views" to be made mutually predictive by contrastive learning; more views yield better representations.
- **Technical summary:** Maximizes mutual information between views via NCE over positive (same-scene) vs negative pairs, scaling to many views; shows representation quality improves as views are added and that view choice matters.
- **Plain-English:** Show a model several angles on the same thing and train it so the angles agree — the more angles, the richer what it learns.
- **Applicability:** A4 (theoretical grounding for treating each omic layer as a "view" and aligning them contrastively — informs multi-omic alignment), A3. Implication: MORPHEUS can add omic layers as additional contrastive views, expecting monotone gains with informative views.
- **Novelty implication:** Foundational; establishes the multi-view contrastive premise MORPHEUS builds on. No direct claim collision.

### 22. VATT: Transformers for Multimodal Self-Supervised Learning from Raw Video, Audio and Text
Akbari, Yuan, Qian, Chuang, Chang, Cui, Gong. NeurIPS 2021. arXiv:2104.11178.
- **Takeaway:** A single modality-agnostic transformer trained by hierarchical contrastive loss on raw video/audio/text, with weights shared across modalities ("one transformer for all").
- **Technical summary:** Linear tokenization of each raw modality feeds a shared transformer; a hierarchical contrastive objective aligns video-audio at a fine level and video-text at a coarser level (different granularities). Modality-agnostic weight sharing works competitively with modality-specific models.
- **Plain-English:** The same transformer, with shared weights, learns from video, sound, and text at once by matching them at appropriate levels of detail.
- **Applicability:** A4 (hierarchical/granularity-aware alignment — biological modalities align at different granularities, e.g., pathway vs gene), A2. Implication: MORPHEUS alignment losses should be hierarchical (programme-level vs feature-level) rather than flat.
- **Novelty implication:** Reframes alignment as multi-granular — useful design input; MORPHEUS's novelty is biological granularity + addressability, not the shared-transformer idea.

### 23. MultiMAE: Multi-modal Multi-task Masked Autoencoders
Bachmann, Mizrahi, Atanov, Zamir. ECCV 2022. arXiv:2204.01678.
- **Takeaway:** A masked-autoencoding pretext across multiple modalities (RGB, depth, semantics) with heavy cross-modal masking learns representations that transfer even when only one modality is available at downstream time.
- **Technical summary:** Randomly masks tokens across all modalities and reconstructs each from the shared visible set, forcing cross-modal predictive coding; a single pretrained encoder then fine-tunes with any subset of the modalities present.
- **Plain-English:** Hide most of the data across all sensors and make the model fill in the gaps, teaching it to infer one modality from another.
- **Applicability:** A4 (masked cross-modal reconstruction is a self-supervised route to *impute missing omic layers* and to learn encode-worthy shared structure), A5 (reconstruction-from-partial is a mild counterfactual: "what would phospho look like given RNA?"). Implication: pretrain MORPHEUS with cross-modal masking so it can impute absent modalities and support conditional queries.
- **Novelty implication:** Strengthens A4/A5 feasibility; pre-empts "cross-modal masked imputation." MORPHEUS's interventional-geometry claim (A5) must go beyond passive imputation to *causal/counterfactual* perturbation.

### 24. Image as a Foreign Language: BEiT Pretraining for Vision and Vision-Language Tasks (BEiT-3)
Wang, Bao, Dong, Bjorck, Peng, Liu, et al. CVPR 2023. arXiv:2208.10442.
- **Takeaway:** Treats images as "Imglish" tokens and uses a single Multiway Transformer with shared self-attention and modality-expert FFNs, unifying masked "language" modeling across text, image, and image-text.
- **Technical summary:** One masked-modeling objective over image tokens, text tokens, and pairs; Multiway Transformer routes each token type to a modality-expert FFN while sharing attention. SOTA across vision and vision-language benchmarks with a single pretraining task.
- **Plain-English:** Turn images into a kind of language and train one "fill in the blank" model over words and image-words together.
- **Applicability:** A2/A4 (Multiway modality-expert FFNs with shared attention = a clean template for MORPHEUS's per-programme/per-modality experts), A1. Implication: adopt Multiway routing so each biomedical modality has an expert while the trunk stays shared and promptable.
- **Novelty implication:** Pre-empts "unified masked modeling with modality-expert routing." MORPHEUS claims biological identifiability of experts + causal queries.

### 25. Uni-Perceiver: Pre-training Unified Architecture for Generic Perception for Zero-shot and Few-shot Tasks
Zhu, Zhu, Chen, Shen, Wang, Dai, et al. CVPR 2022. arXiv:2112.01522.
- **Takeaway:** Casts *every* task (across modalities) as maximizing similarity between an input and a candidate-target representation in one shared space, enabling zero/few-shot generalization to *unseen tasks* without task-specific heads.
- **Technical summary:** A single encoder maps inputs and candidate targets (labels, captions, next tokens) into a joint space; the training objective is a unified maximum-likelihood over input-target similarity, so novel tasks are handled by supplying new candidate targets — no new head.
- **Plain-English:** Frame all tasks as "which candidate answer best matches this input?" so the model can tackle brand-new tasks just by offering it new candidate answers.
- **Applicability:** A1 (the closest prior to MORPHEUS's *task auto-detection/routing without hard-coded probes* — tasks are specified by their candidate-target set, not a fixed head), A2. Implication: MORPHEUS can implement NL task auto-detection as choosing/generating the candidate-target space that an NL request implies.
- **Novelty implication:** **Pre-empts the core of A1** ("no task-specific heads; specify the task and route it"). MORPHEUS must differentiate on *inferring the task from free-form NL* (Uni-Perceiver still needs the task's candidate set supplied) and on biological/pathway addressability + causal queries.

### 26. data2vec: A General Framework for Self-supervised Learning in Speech, Vision and Language
Baevski, Hsu, Xu, Babu, Gu, Auli. ICML 2022. arXiv:2202.03555.
- **Takeaway:** One self-supervised recipe — predict latent target representations of the full input from a masked view via a teacher-EMA — works across speech, vision, and text with the same objective.
- **Technical summary:** A student predicts contextualized latent targets produced by an EMA teacher on the unmasked input; because targets are *learned latents* (not modality-specific tokens), the identical objective spans modalities and captures higher-order features.
- **Plain-English:** Rather than predicting raw pixels or words, predict the model's own internal summary of the hidden parts — a trick that works the same way for sound, images, and text.
- **Applicability:** A4 (a single modality-agnostic SSL objective that could pretrain heterogeneous biomedical modalities uniformly), A2. Implication: latent-target prediction gives MORPHEUS a common pretraining objective across omics without hand-designed per-modality targets.
- **Novelty implication:** Foundational; supports the feasibility of a unified biomedical SSL objective. No direct claim collision.

### 27. Joint probabilistic modeling of single-cell multi-omic data with totalVI
Gayoso, Steier, Lopez, Regier, Nazor, Streets, Yosef. Nature Methods 18, 2021.
- **Takeaway:** A conditional VAE jointly models paired RNA + surface protein (CITE-seq), decomposing biological vs technical (protein background, batch) factors into one latent cell-state representation.
- **Technical summary:** Per-cell latent variable generates both mRNA (negative binomial) and protein (mixture accounting for background) via neural decoders; scalable stochastic variational inference handles millions of cells, supports imputation, denoising, differential expression, and batch integration.
- **Plain-English:** For cells measured by both their genes and their surface proteins, it learns one unified "cell state" that explains both while stripping out technical noise.
- **Applicability:** A4 (canonical *encoded* biomedical multimodal fusion — RNA+protein into a shared latent; a reference for what "encode this modality" concretely buys), A2 (latent cell state as addressable representation). Implication: MORPHEUS's protein/phospho encoding path can follow totalVI-style modality-specific likelihoods into a shared latent.
- **Novelty implication:** Strengthens A4 as prior art for encoding proteomics; MORPHEUS must add promptability, pathway addressability, and causal queries that totalVI lacks.

### 28. MultiVI: Probabilistic joint analysis of single-cell multimodal data (Nature Methods)
Ashuach, Gabitto, Koodli, Saldi, Jordan, Yosef. Nature Methods 20, 2023.
- **Takeaway:** Extends the VAE framework to jointly model RNA + ATAC (+ optionally protein), and crucially *integrates paired and unpaired* cells by aligning modality-specific latents into a shared space.
- **Technical summary:** Modality-specific encoders/decoders (NB for RNA, Bernoulli for ATAC) map to a shared latent with a penalty aligning per-modality latents; missing modalities are imputed, and unpaired single-modality cells are embedded jointly with paired cells.
- **Plain-English:** It puts cells measured by genes, by chromatin, or by both into a single shared map — even when most cells only have one measurement.
- **Applicability:** A4 (directly addresses MORPHEUS's core problem: fuse modalities when most samples are *partially observed*; alignment-of-latents is a concrete encode-vs-impute mechanism), A2. Implication: MORPHEUS should align modality-specific latents so unpaired/partial multi-omic samples still land in one space.
- **Novelty implication:** Pre-empts "probabilistic fusion of partially paired multi-omics." MORPHEUS differentiates via NL promptability, task auto-detection, and interventional queries.

### 29. Multi-omics single-cell data integration and regulatory inference with graph-linked embedding (scGLUE)
Cao, Gao. Nature Biotechnology 40, 2022.
- **Takeaway:** Integrates *unpaired* multi-omics (RNA, ATAC, methylation) by embedding cells into a shared space linked through a prior *knowledge graph* of regulatory relationships across feature spaces.
- **Technical summary:** Per-omics VAEs share a latent space; a guidance graph encoding cross-omics feature relations (e.g., gene-peak links) plus adversarial alignment ties the modalities, enabling integration without any paired cells and yielding regulatory inferences.
- **Plain-English:** It aligns different molecular measurements from *different* cells by using known biological links between their features as a bridge.
- **Applicability:** A4 (fusion via a *knowledge graph prior* — a bridge between encoding and RAG: structured biological knowledge guides alignment), A2/A3 (graph-linked features are inherently addressable and grounded). Implication: MORPHEUS can use a regulatory-graph prior to align modalities and to make slots pathway-addressable.
- **Novelty implication:** Reframes the encode-vs-RAG boundary — knowledge-graph-guided alignment is a hybrid MORPHEUS should position against; supports A2/A3 while pre-empting "graph-guided multi-omic alignment."

### 30. Contrastive Learning of Medical Visual Representations from Paired Images and Text (ConVIRT)
Zhang, Jiang, Miura, Manning, Langlotz. MLHC 2022 (arXiv:2010.00747, 2020).
- **Takeaway:** The pre-CLIP demonstration that bidirectional image-text contrastive learning on paired radiology reports yields label-efficient medical image representations.
- **Technical summary:** Paired chest X-ray + report; a bidirectional contrastive objective aligns image and text encoders; downstream classification/retrieval needs far fewer labels than ImageNet-transfer baselines. Direct ancestor of biomedical CLIP variants.
- **Plain-English:** Pairing medical images with their written reports and training them to match produces strong image features without lots of manual labels.
- **Applicability:** A3/A4 (foundational proof that *biomedical* NL (reports) can supervise a representation via contrastive alignment — the biomedical grounding precedent), A1. Implication: MORPHEUS's NL grounding for biology has direct lineage; report-style text can supervise representation alignment.
- **Novelty implication:** Foundational for biomedical NL grounding; MORPHEUS must claim emergent-knowledge elicitation + evaluation beyond contrastive report alignment.

### 31. Making the Most of Text Semantics to Improve Biomedical Vision-Language Processing (BioViL)
Boecking, Usuyama, Bannur, Castro, Schwaighofer, Hyland, Wetscherek, Naumann, Nori, Alvarez-Valle, Poon, Oktay. ECCV 2022. arXiv:2204.09817.
- **Takeaway:** Shows that a domain-specific text model and text-aware losses (exploiting radiology report semantics, incl. a new CXR-BERT and local alignment) substantially improve biomedical vision-language grounding and phrase grounding.
- **Technical summary:** Pretrains a specialized clinical text encoder (CXR-BERT) and adds a local/global image-text contrastive plus masked-language objective; releases a phrase-grounding benchmark (MS-CXR). Text-quality gains transfer to zero-shot and grounding tasks.
- **Plain-English:** Investing in a language model that truly understands clinical text — not generic text — makes the whole image-text system markedly better and more precise about *where* findings are.
- **Applicability:** A3 (domain-specialized NL encoder is essential for biological grounding; generic NL under-grounds — a direct design constraint for MORPHEUS's NL interface), A2 (local alignment ~ addressable regions/slots). Implication: MORPHEUS needs a biology-specialized text encoder, not off-the-shelf NL, and local (slot-level) alignment objectives.
- **Novelty implication:** Reframes — establishes that biomedical NL grounding quality gates everything; strengthens MORPHEUS's A3 emphasis while pre-empting "domain text encoder improves grounding."

### 32. A visual-language foundation model for pathology image analysis using medical Twitter (PLIP)
Huang, Bianchi, Yuksekgonul, Montine, Zou. Nature Medicine 29, 2023.
- **Takeaway:** Fine-tunes CLIP on ~200K pathology image-text pairs harvested from medical Twitter (OpenPath), yielding a pathology foundation model with strong zero-shot classification and image-text retrieval.
- **Technical summary:** Curates OpenPath from public pathology social-media posts with hashtags; contrastive CLIP fine-tuning produces PLIP, which does zero-shot tissue classification and text-to-image / image-to-image retrieval, outperforming generic CLIP on pathology.
- **Plain-English:** By collecting pathologists' publicly shared annotated images, they taught a CLIP-style model the language of pathology, so it can recognize tissue types just from text descriptions.
- **Applicability:** A3/A4 (proof that *biological* zero-shot prompting works when the contrastive space is trained on domain image-text; informs whether to encode a modality by pairing it with NL), A1. Implication: for modalities with harvestable NL pairings, contrastive encoding gives promptable zero-shot biology.
- **Novelty implication:** Pre-empts "domain CLIP gives biological zero-shot prompting." MORPHEUS must claim multimodal (non-image) omics + task auto-detection + causal queries beyond image-text zero-shot.

### 33. A visual-language foundation model for computational pathology (CONCH)
Lu, Chen, Williamson, Chen, Liang, Mahmood. Nature Medicine 30, 2024.
- **Takeaway:** A pathology vision-language foundation model pretrained on ~1.17M image-caption pairs with an iterative CoCa-style contrastive+captioning objective, SOTA across 14 diverse pathology benchmarks including zero-shot subtyping and text-to-image retrieval.
- **Technical summary:** Curates a large histopathology image-caption corpus; trains a CoCa-like model (contrastive alignment + generative captioning) with task-agnostic pretraining, enabling zero-shot classification, segmentation-via-retrieval, and captioning without task-specific heads.
- **Plain-English:** A general pathology model that both matches images to text and describes them, and can be pointed at many new pathology tasks just by prompting.
- **Applicability:** A1/A3 (a task-agnostic biomedical foundation model promptable in NL across many tasks — the closest biomedical analogue to MORPHEUS's promptable-rep ambition, in the image-text sub-domain), A4. Implication: benchmark MORPHEUS's promptability against CONCH-style zero-shot task breadth, extended to molecular modalities.
- **Novelty implication:** **Pre-empts "task-agnostic promptable biomedical foundation model"** in pathology. MORPHEUS's defensible novelty is molecular/omic modalities + task auto-detection from free NL + identifiable pathway slots + interventional queries — none of which CONCH addresses.

### 34. PaLI: A Jointly-Scaled Multilingual Language-Image Model
Chen, Wang, Changpinyo, Piergiovanni, Padlewski, Salz, et al. ICLR 2023. arXiv:2209.06794.
- **Takeaway:** Scales a unified image+text-to-text model where *all* vision-language tasks are cast as text generation, showing joint scaling of vision and language components drives broad multitask/multilingual gains.
- **Technical summary:** A large ViT feeds an encoder-decoder text transformer; every task (captioning, VQA, detection-as-text, OCR) is a text-generation target, trained on a huge multilingual image-text mixture. Establishes benefits of scaling the vision component alongside language.
- **Plain-English:** Turn every image task into "generate the right text," then scale both halves — the model gets broadly capable across languages and tasks.
- **Applicability:** A1 (unifying tasks as NL generation is a route to task auto-detection: the model reads the request and generates the answer), A3. Implication: MORPHEUS can express heterogeneous biological queries as a single text-generation interface over the shared rep.
- **Novelty implication:** Pre-empts "all tasks as text-in/text-out over a multimodal trunk." MORPHEUS differentiates on biological grounding, addressable slots, and causal-geometry queries rather than generation-unification.
