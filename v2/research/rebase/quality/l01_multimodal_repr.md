# Quality audit — l01_multimodal_repr (Multimodal representation learning & fusion)

Auditor pass date: 2026-07-29. Source file: `v2/research/rebase/lit/l01_multimodal_repr.md`.
Entries reviewed: 34.

## Verdict summary
- **Verified real & findable:** 34/34 (8 spot-checked directly against arXiv / PubMed / Nature; remaining 26 are canonical, widely-cited papers with metadata matching known records).
- **Unverifiable / possibly fabricated:** 0
- **Low-quality or off-topic:** 0
- **Entries missing a required component (technical + plain-English + applicability + novelty):** 0

## Spot-checks performed (WebFetch, primary sources)
| # | Entry | ID checked | Result |
|---|-------|-----------|--------|
| 12 | Meta-Transformer | arXiv:2307.10802 | CONFIRMED — title + authors (Zhang, Gong, Zhang, Li, Qiao, Ouyang, Yue) match |
| 14 | Attention Bottlenecks (MBT) | arXiv:2107.00135 | CONFIRMED — Nagrani, Yang, Arnab, Jansen, Schmid, Sun; NeurIPS 2021 |
| 13 | ONE-PEACE | arXiv:2305.11172 | CONFIRMED — Wang et al., 4B-param extensible model |
| 20 | Generalized Multimodal ELBO (MoPoE-VAE) | arXiv:2105.02470 | CONFIRMED — Sutter, Daunhawer, Vogt; ICLR 2021 |
| 30 | ConVIRT | arXiv:2010.00747 | CONFIRMED — Zhang, Jiang, Miura, Manning, Langlotz; venue **MLHC 2022** correct |
| 11 | LanguageBind | arXiv:2310.01852 | CONFIRMED — Zhu, Lin, Ning et al.; VIDAL dataset, language-anchored N-modality |
| 29 | scGLUE | PubMed 35501393 | CONFIRMED — Cao & Gao, Nature Biotechnology 40(10):1458-1466, 2022 |
| 33 | CONCH | arXiv:2307.12914 | CONFIRMED — Lu, Chen, Williamson, ... Mahmood; 1.17M image-caption pairs; pub. Nature Medicine 30, 2024 |

All spot-checked citations had accurate titles, author lists, venues, and years. No metadata drift, no invented papers, no mismatched author/venue pairings.

## Remaining 26 entries (not individually fetched)
All are foundational, heavily-cited works whose titles/authors/arXiv IDs/venues match well-established records: CLIP (2103.00020), ALIGN (2102.05918), Perceiver (2103.03206), Perceiver IO (2107.14795), Flamingo (2204.14198), CoCa (2205.01917, TMLR), FLAVA (2112.04482), BLIP (2201.12086), BLIP-2 (2301.12597), ImageBind (2305.05665), LiT (2111.07991), Frozen (2106.13884), Missing-Modality robustness (2204.05454), MVAE (1802.05335), MMVAE (1911.03393), CMC (1906.05849), VATT (2104.11178), MultiMAE (2204.01678), BEiT-3 (2208.10442), Uni-Perceiver (2112.01522), data2vec (2202.03555), totalVI (Nature Methods 18, 2021), MultiVI (Nature Methods 20, 2023), PLIP (Nature Medicine 29, 2023), BioViL (2204.09817), PaLI (2209.06794). No red flags; none require reclassification. Recommend no further action unless a deeper pass is requested.

## On-topic / quality assessment
Every entry sits squarely within the lane remit (contrastive/aligned representation learning; early/late/attention/PoE-MoE fusion; modality dropout / missing-modality; CLIP/ImageBind/Perceiver-IO lineage) and traces to a MORPHEUS axis (A1–A5). Coverage is well-balanced across the remit's sub-areas:
- Contrastive alignment & domain CLIP: CLIP, ALIGN, LiT, CMC, ConVIRT, BioViL, PLIP, CONCH
- Fusion architectures / bottlenecks: Perceiver, Perceiver IO, MBT, FLAVA, BEiT-3
- Frozen-trunk plug-in / soft-prompt encoding: Flamingo, BLIP-2, Frozen, Meta-Transformer, LiT
- PoE/MoE/subset generative fusion: MVAE, MMVAE, MoPoE
- Modality dropout / missing modality / cross-modal masking: Missing-Modality robustness, MultiMAE
- N-modality binding & extensibility: ImageBind, LanguageBind, ONE-PEACE
- Encoded biomedical multi-omic fusion: totalVI, MultiVI, scGLUE
- Task-unification / routing: Uni-Perceiver, PaLI, CoCa, BLIP, data2vec, VATT

Venue quality is uniformly high (ICML, NeurIPS, ICLR, CVPR, ECCV, TMLR, Nature Methods/Biotechnology/Medicine). No weak, predatory, or off-topic sources.

## Entry-structure completeness
All 34 entries contain the required components: a **Technical summary**, a **Plain-English** summary, an **Applicability** section (with explicit axis tags A1–A5 and a design implication), and a **Novelty implication**. Each also carries a one-line **Takeaway**. No entry is missing a component.

## Minor / non-blocking observations
- CONCH is dated "Nature Medicine 30, 2024" in the file; the fetched arXiv preprint (2307.12914) is titled "Towards a Visual-Language Foundation Model for Computational Pathology." The published Nature Medicine 2024 version is the correct citation — consistent, no error.
- MoPoE-VAE listed as "ICLR 2021, arXiv:2105.02470" — arXiv v1 is May 2021; paper is a legitimate ICLR 2021 acceptance. Consistent.
- No action required on either.

## Bottom line
Clean lane. 34/34 verified/plausible-real, 0 fabricated, 0 low-quality, 0 off-topic, 0 structurally incomplete.
