# Quality Audit — l08_multimodal_rag

Referee pass on `lit/l08_multimodal_rag.md` (26 entries). Adversarial verification of
existence, on-topic quality, and required structure.

## Verdict summary
- **Entries: 26. Spot-checked for existence: 8. Verified real & accurately described: 8/8.**
- **Unverifiable / possibly fabricated: 0.**
- **Low-quality / off-topic: 0** (one mild framing stretch noted, already self-flagged by author).

## Existence spot-checks (WebFetch on arXiv / ACL Anthology)
All confirmed with matching title + author list + topic:

| # | Paper | ID | Result |
|---|-------|----|--------|
| 8 | Nonparametric Masked Language Modeling (NPM) | arXiv:2212.01349 | ✅ exact match (Min, Shi, Lewis, Chen, Yih, Hajishirzi, Zettlemoyer) |
| 14 | Retrieval-Augmented Multimodal LM (RA-CM3) | arXiv:2211.12561 | ✅ exact match (Yasunaga et al.), ICML 2023 |
| 16 | MMed-RAG | arXiv:2410.13085 | ✅ real (Xia, Zhu, Li … Zou, Yao); true title has "System" ("Versatile Multimodal RAG System…") — trivial wording drop, not material |
| 18 | HeteroRAG + MedAtlas | arXiv:2508.12778 | ✅ exact match (Chen, Liao, Zhu … Wang, Wang) |
| 20 | RA-RRG | arXiv:2504.07415 | ✅ exact match (Park, Yoon, Kim, Choi) |
| 21 | AMANDA | arXiv:2510.02328 | ✅ exact match (Wang, Mao, Wen, Luo, Ding) |
| 23 | BioBridge | arXiv:2310.03320 | ✅ exact match (Wang, Wang, Srinivasan, Ioannidis, Rangwala, Anubhai), ICLR 2024 |
| 25 | scRAG | 2025.findings-acl.53 | ✅ exact match (Yu, Zheng, Chen, Hua, Luo), Findings ACL 2025 pp.954–970 |

Note: entry 26 GenePT (bioRxiv 2023.10.16.562533) — biorxiv returned HTTP 403 so not machine-
confirmed this pass, but it is a well-known real preprint (Chen & Zou) and the entry itself already
flags "preprint — verify final venue before citing as peer-reviewed." No concern.

The remaining un-fetched entries (1 kNN-LM, 2 RETRO, 3 REALM, 4 RAG, 5 Atlas, 6 Memorizing
Transformers, 7 In-Context RALM, 9 Why-do-kNN-LMs-Work, 10 Fine-Tuning-or-Retrieval, 11 REML,
12 Self-RAG, 13 RAG Survey, 15 MedRAG/MIRAGE, 17 RULE, 19 FactMM-RAG, 22 Biomedicine RAG survey,
24 RetMol) are canonical/high-visibility works with correct-looking IDs, venues, and author lists;
no fabrication signals. Budget note: WebSearch quota was exhausted, so verification used direct
WebFetch on abstract pages (more precise for existence checks).

## Quality / on-topic assessment
- **Remit fit: strong.** Every entry maps to the retrieval-vs-parametric-encoding tradeoff and its
  multimodal/biomedical instantiation — the exact A4 focus. Clean four-part arc
  (lineage → encode-vs-retrieve theory → biomedical/multimodal RAG → retrieval-augmented biology).
- **Structure: complete.** All 26 entries carry the required fields — a *technical summary*, a
  *plain-English* gloss, an *applicability* note (mapped to MORPHEUS axes A1–A5), and a *novelty
  implication* (honestly labeled Pre-empts / Strengthens / Reframes). No entry is missing a component.
- **Intellectual honesty: high.** The review repeatedly and correctly concedes the mechanism is prior
  art (entries 1–13) and confines MORPHEUS novelty to a biological modality-selection rule; it even
  includes devil's-advocate entries (9 "kNN gains can be distilled back", 26 GenePT as adversarial
  encode-vs-retrieve baseline).

## Flagged entries
- **Entry 26 (GenePT) — mild framing stretch, not a defect.** GenePT is a *text-embedding
  representation* method, not retrieval in the datastore/kNN sense; calling it "retrieved" leans on a
  loose reading of "retrieve from the LLM/text side." The author already brackets it as a baseline and
  flags its preprint status, so the stretch is disclosed rather than hidden. Keep, but do not cite it
  as a retrieval-vs-encode result without that caveat.
- **Entries 18 & 20 — future venue ("ACL 2026 Findings").** Plausible given the 2026 date; arXiv IDs
  (Aug 2025 / Apr 2025) verified. Confirm final acceptance before citing the venue as settled.

No entries recommended for removal.

## Counts
- Unverifiable / fabricated: **0**
- Low-quality / off-topic: **0** (1 disclosed framing caveat, 2 pending-venue notes)
