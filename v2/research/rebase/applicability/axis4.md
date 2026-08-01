# Axis A4 — Multimodal prompting: ENCODE vs RETRIEVE, and the frozen-trunk plug-in

**Axis remit.** For each auxiliary molecular modality (proteomics, phospho, CNV, SNV, and by extension bulk RNA, ATAC, spatial), decide *when to ENCODE it* into the shared trunk (a learned adapter/latent) vs *when to treat it as RAG context* (retrieved evidence around a frozen trunk). Includes the frozen-trunk plug-in cost model and the mosaic/missing-modality reality of clinical multi-omics.

This synthesis reads across all 15 harvested lanes and extracts what actually bears on A4. The dominant sources are `l08_multimodal_rag` (the encode-vs-retrieve backbone), `l01_multimodal_repr` (fusion/frozen-plug-in mechanisms), `l03_omics_fm` (biological encode/retrieve baselines), and `l15_steelman_prior_art` (what is already claimed). Supporting hits come from `l04_molecular_nl_prompting` and `l05_promptable_unified`.

---

## Key findings (with citations)

### The tradeoff itself is settled prior art — the mechanism is not the contribution
- **kNN-LM** (Khandelwal et al., ICLR 2020) is the founding evidence that a *frozen trunk + external datastore* beats pushing everything into parameters, and gains concentrate on the **rare / long-tail / factual** slice. A datastore swap = per-cohort adaptation with no retraining. (`l08` #1)
- **RETRO** (Borgeaud et al., ICML 2022): a 7.5B model with chunked cross-attention over 2T retrieved tokens matches GPT-3 (175B) — "encode the reasoning core small, retrieve the knowledge." (`l08` #2)
- **RAG** (Lewis et al., NeurIPS 2020) and **In-Context RALM** (Ram et al., TACL 2023) define "modality as retrieved context"; In-Context RALM is the *zero-surgery baseline* — prepend retrieved evidence to a frozen model, no training — that any "encode this modality" proposal must beat. (`l08` #4, #7)
- **Fine-Tuning or Retrieval?** (Ovadia et al., EMNLP 2024): for injecting *new* knowledge, RAG consistently beats unsupervised fine-tuning and avoids catastrophic forgetting. Direct decision rule for a late-arriving cohort (e.g., a new CPTAC batch) added *after* the trunk is frozen. (`l08` #10)
- **Implication:** MORPHEUS cannot claim the encode-vs-retrieve mechanism, the frozen-trunk-plus-datastore pattern, or "retrieval beats fine-tuning for new knowledge." These are 4–6 years old in NLP. A4's defensible novelty must be a *molecular-modality-selection criterion*.

### Long-tail / data-scarce modalities are exactly where retrieval wins — the sharpest quantitative case
- **NPM** (Min et al., ACL 2023 Findings) *replaces* the parametric softmax with a nonparametric distribution over corpus phrases; largest wins on rare words, rare senses, rare facts. Reframes A4 from "add RAG context" to "retrieve *instead of* encode a modality." (`l08` #8)
- **Atlas** (Izacard et al., JMLR 2023): retrieval buys *few-shot* competence (42% NQ with 64 examples, beating a 540B model) where an encoder would overfit. (`l08` #5)
- **Implication:** sparsely-observed, heavy-tailed molecular modalities — CPTAC proteomics (inventory-only), rare phospho-sites, rare mutation combinations — map onto the exact regime where retrieval out-generalizes an encoder that must spend parameters on the tail.

### But the retrieval advantage may be an artifact of under-trained encoders — the evidentiary bar
- **Why do Nearest Neighbor LMs Work?** (Xu, Alon, Neubig 2023): kNN-LM gains decompose into three mechanistic factors that can be **largely internalized into the parametric model** via better training. (`l08` #9)
- **Implication (critical devil's advocate):** before MORPHEUS declares a modality "retrieve-only," it must rule out that a better-trained encoder captures the same signal. A genuine "must-retrieve" claim needs a modality where the gap *survives distillation*. This is the honest bar A4's study must clear.

### Frozen-trunk plug-in mechanisms are mature and diverse — pick, don't invent
- **BLIP-2 Q-Former** (Li et al., ICML 2023): a lightweight querying transformer with learned query tokens bridges a frozen encoder and a frozen LLM, beating Flamingo-80B with ~54× fewer trainable params. Strongest evidence for "encode a modality as a handful of soft-prompt tokens." (`l01` #9; echoed in `l04` — the Q-Former is a concrete template for the `(batch, n_pathways, D)` slots the TQI needs)
- **Frozen** (Tsimpoukelli et al., NeurIPS 2021): the seminal "encode a modality as a soft prompt for a frozen LM." **LiT** (Zhai et al., CVPR 2022): lock the strong encoder, tune only the text/alignment side. **Flamingo** (Alayrac et al., 2022): gated cross-attention adapters into a frozen LM. **Meta-Transformer** (2023): only tokenizer + head trained per modality against a frozen backbone. **BEiT-3 / ONE-PEACE**: shared attention + per-modality expert FFNs, extensible without disturbing prior alignment. (`l01` #16, #15, #5, #12, #24, #13)
- **Implication:** the plug-in mechanism itself is prior art. Encoding a new molecular modality can cost only a tokenizer/Q-Former + adapter against a frozen (even quantized — Geneformer 2026) oncology trunk. Novelty is *which* modality merits this vs RAG, plus pathway-addressability of the query tokens.

### Missing / mosaic modalities are the chronic clinical reality — the design must factor absence
- **Are Multimodal Transformers Robust to Missing Modality?** (Ma et al., CVPR 2022): sharp degradation under missing modalities; no single fusion is universally best; modality-dropout training is required. (`l01` #17)
- **PoE / MMVAE / MoPoE** (Wu & Goodman 2018; Shi et al. 2019; Sutter et al., ICLR 2021): principled subset-marginalizing fusion — train once, serve *any observed subset* of {RNA, protein, phospho, CNV, SNV}. MoE guards against a dominant modality drowning weaker programmes. (`l01` #18–20)
- **MultiVI / totalVI / scGLUE** (Nature Methods 2023 / 2021; Nature Biotech 2022): the biology-side *encode* baselines for partially-paired/mosaic multi-omics; scGLUE uses a regulatory *knowledge-graph prior* as a bridge between pure encoding and RAG. (`l03` #19, #18, #10)
- **Implication:** clinical samples rarely have every omic. A4's choice must be posed as "encode / retrieve / marginalize" per available subset, trained with modality dropout, with validated graceful degradation.

### The biomedical-multimodal-RAG surface is crowded — but imaging/text-centric
- **RA-CM3** (ICML 2023), **MedRAG/MIRAGE** (2024), **MMed-RAG** (ICLR 2025), **RULE** (EMNLP 2024, +47.4% factual accuracy via *calibrated retrieval count*), **HeteroRAG** (heterogeneous per-source retrievers), **FactMM-RAG** (retrieval keys encode *biological fact structure* via RadGraph), **RA-RRG** (retrieve fine-grained *phrases*, not documents). (`l08` #14–20)
- The biomedical-RAG survey (He et al. 2025) explicitly names **"multimodal RAG misalignment"** as unsolved and frames a reasoning/latency/privacy trilemma. (`l08` #22)
- **Implication:** all of this is imaging+text. **No cited work retrieves over *structured molecular* modalities that have no natural text/image form.** That is the open lane, named as an open problem by the survey. RULE/MMed-RAG's calibrated "how much to trust retrieved vs internal" is prior art to *build on* for a safety gate, not to claim.

### Biology-side "retrieve instead of encode" already has a devastating null
- **GenePT** (Chen & Zou, Nat. Biomed. Eng. 2024): representing genes/cells by *retrieved* LLM text-embeddings of gene descriptions — with **no expression pretraining** — rivals or beats Geneformer/scGPT. A live null hypothesis: if a text/RAG baseline matches an encoded latent, the "encode" side of A4 is unjustified for that modality. (`l03` #8; `l15` #19)
- **scRAG** (ACL 2025 Findings): retrieve reference cells + KG triples rather than encode, beating trained classifiers on *unseen tissues*. **BioBridge** (ICLR 2024): parameter-efficient KG bridges connect *frozen* unimodal biomedical FMs without fine-tuning them. **RetMol** (ICLR 2023): steer a frozen generator with retrieved exemplars, no fine-tuning. (`l08` #25, #23, #24)
- **Implication:** biology-side retrieval exists but keys on **raw similarity or generic external KGs**, never on an *identified causal tumor-state*. Retrieval keyed by A2 pathway slots is the unclaimed bridge between A2 and A4.

### It is a spectrum, not a binary
- **Memorizing Transformers** (Wu et al., ICLR 2022): a modality can be *memory read via attention* — the middle ground between fused-encode and retrieved-text. Reframes A4 as **fuse ↔ attend-to-memory ↔ retrieve-text**, and MORPHEUS can claim novelty in *placing biological modalities on this spectrum*. (`l08` #6)

### Modality-specific molecular encoders are mature (the "encode via existing model" pole)
- **Prosit** (Nat. Methods 2019), **DeepPhospho** (Nat. Commun. 2021), **Phosformer** (Bioinformatics 2023), **ESM-2** (Science 2023), **PINNACLE** (Nat. Methods 2024): proteomics/phospho/protein DL representations are strong and reusable. Phospho fragmentation is modality-specific enough that it **cannot be naively folded into an RNA trunk**. (`l03` #21–24, #29)
- **Implication:** the "encode" path can lean on *frozen existing* proteomics/phospho encoders feeding a light adapter, rather than relearning them — cheapening the encode option and sharpening the comparison against retrieval.

### In-domain oncology pre-emptions to differentiate against
- **SurvPath** (CVPR 2024): already tokenizes transcriptomics into **named biological pathway tokens** fused with histology — but tokens are *defined a priori by gene sets*, not *identified*, and are neither NL-promptable nor counterfactual. **PORPOISE** (Cancer Cell 2022): the "encode both modalities and fuse" default. **Med-Gemini** (2024): a production system that already *blends encoded modalities with retrieved context* — but never asks *which molecular modality* to encode vs retrieve. The **Boehm et al. (Nat. Rev. Cancer 2022)** multimodal-oncology fusion taxonomy **omits RAG-as-a-modality-choice entirely**. (`l15` #31, #32, #10, #33)

---

## What this implies for MORPHEUS design

**Build ON (established, low-risk, cite-don't-claim):**
1. **Frozen-trunk plug-in via Q-Former / gated cross-attention / locked-encoder tuning** (BLIP-2, Flamingo, LiT, Frozen). Encode a candidate modality as a few pathway-addressable soft-prompt tokens against a frozen (quantized) trunk.
2. **Subset-marginalizing fusion (PoE / MoPoE) + modality-dropout training** (Wu&Goodman, Sutter, Ma et al.) so any observed multi-omic subset is served by one model with graceful degradation.
3. **Modality-specific frozen molecular encoders** (ESM-2, DeepPhospho, Phosformer, Prosit, PINNACLE) feeding light adapters — don't relearn proteomics/phospho.
4. **Calibrated retrieval + safety gate** (RULE, MMed-RAG, Self-RAG) for the closed-RAG card renderer.

**NEW-DESIGN (the genuine, defensible A4 gap):**
1. **A principled, pre-registered molecular modality-selection rule** — no cited work formalizes *when* proteomics/phospho/CNV/SNV should be encoded vs retrieved. The theory (long-tail favors retrieval: kNN-LM, NPM, Atlas) and the null (GenePT) exist in fragments; the *biological instantiation with a decision criterion* does not.
2. **Retrieval over structured molecular memory with no text/image form** — the survey-named "multimodal RAG misalignment" gap. All medical RAG is imaging/text.
3. **Retrieval keyed by identified A2 pathway/causal slots** rather than raw similarity (differentiates from scRAG, BioBridge, GenePT).

**Hard evidentiary bars A4 must clear (or the axis is hype):**
- Beat the **In-Context RALM zero-surgery baseline** — if prepending retrieved proteomic evidence matches an encoded adapter, the adapter hasn't earned its place.
- Beat the **GenePT text-embedding null** for any modality it recommends *encoding*.
- Show any "must-retrieve" gap **survives distillation** into the trunk (Xu/Alon/Neubig control arm) — otherwise the gap is just an under-trained encoder.

---

## Candidate research directions

### D1. The Modality Encodability Score (MES): a pre-registered encode-vs-retrieve decision rule for molecular modalities
**Claim.** There exists a computable, per-modality score — a function of (a) cohort coverage/sparsity, (b) tail-heaviness of the signal distribution, (c) the *survives-distillation* gap (does an encoder capture signal a retrieval head cannot after the retrieval trick is distilled into the trunk, à la Xu/Alon/Neubig), and (d) margin over the In-Context RALM zero-surgery baseline — that *predicts* whether {proteomics, phospho, CNV, SNV, bulk-RNA} should be encoded into the frozen oncology trunk or served as retrieved context, and this rule generalizes to held-out modalities/cohorts.
**Why novel.** The encode-vs-retrieve theory (kNN-LM, NPM, Atlas), the distillation caveat (Xu/Alon/Neubig), the fine-tune-vs-RAG result (Ovadia), and the GenePT null all exist — but *no* cited work turns them into a quantitative, testable *modality-selection criterion for molecular biology*. The steelman lane and both RAG/omics syntheses independently flag this exact rule as the single unoccupied A4 contribution. It converts A4 from an engineering choice into a falsifiable, pre-registerable claim with concrete baselines (In-Context RALM, GenePT) built in.

### D2. Pathway-slot-keyed molecular memory: placing modalities on the fuse ↔ attend ↔ retrieve spectrum
**Claim.** Structured molecular modalities with no natural text/image form (phospho, CNV, SNV) are best integrated as an **attention-readable external memory keyed by identified A2 pathway/programme slots** — a Memorizing-Transformers-style middle ground — rather than either dense fusion (SurvPath's a-priori pathway tokens) or free-text RAG. Retrieval keys are *identified causal tumor-state slots*, not raw expression similarity.
**Why novel.** Directly targets the biomedical-RAG survey's named-open "multimodal RAG misalignment": every existing multimodal/medical RAG system (RA-CM3, MMed-RAG, HeteroRAG, scRAG, BioBridge) retrieves images/text or keys on raw similarity / generic KGs. None retrieve over structured molecular memory addressed by an *identified* representation. This is the concrete bridge between A2 (addressable slots) and A4 that `l08`'s synthesis says "none pre-empt, but depends on A2 delivering identifiability" — making it the highest-leverage A2×A4 claim, and it reframes A4 as a spectrum rather than a binary.

### D3. Mosaic frozen-trunk plug-in with an adaptive per-query encode/retrieve/marginalize gate
**Claim.** A single frozen-trunk plug-in accepts *any observed subset* of molecular modalities via MoPoE subset-marginalization, is trained with modality dropout, and carries a **calibrated per-query gate** (Self-RAG / RULE-style) that decides, for each *absent* modality and each NL task, whether to retrieve it as context, marginalize it, or abstain — with validated graceful degradation and no catastrophic reliance on any single modality.
**Why novel.** The mosaic-fusion VAEs (MultiVI, totalVI, MoPoE) handle *encoded* missing modalities but have no retrieval path and no NL-conditioned gate; the medical-RAG calibration work (RULE, MMed-RAG, Self-RAG) handles adaptive retrieval but only for imaging/text and never for molecular mosaics. Combining subset-marginalizing fusion + an adaptive retrieve-or-marginalize gate + calibrated missing-modality behavior on a frozen oncology trunk is unclaimed, and it operationalizes the clinical reality (Ma et al.'s missing-modality collapse) that the "encode everything" oncology default (PORPOISE) and the RAG-free fusion taxonomy (Boehm et al.) both ignore.
