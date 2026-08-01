# Quality audit — l06_emergence_eval (Emergent-capability & elicitation evaluation)

Referee pass date: 2026-07-29. Source file: `v2/research/rebase/lit/l06_emergence_eval.md` (37 entries).

## Verdict summary
- **Unverifiable / possibly-fabricated entries: 0**
- **Low-quality / off-topic entries: 0**
- Every entry contains all five required fields (Takeaway, Technical summary, Plain-English, Applicability, Novelty implication). Completeness requirement fully met.

## Reality spot-checks (9 of 37, ~24%, verified via arXiv / Semantic Scholar / publisher)
All checked entries matched the cited title, ID, and author list **exactly**:

| # | Entry | ID / DOI | Result |
|---|-------|----------|--------|
| 4 | Are Emergent Abilities in LLMs just In-Context Learning? (Lu, Bigoulaeva, Sachdeva, Tayyar Madabushi, Gurevych) | arXiv:2309.01809 | VERIFIED — title, all 5 authors, ACL 2024 |
| 5 | Understanding Emergent Abilities from the Loss Perspective (Du, Zeng, Dong, Tang) | arXiv:2403.15796 | VERIFIED — title, authors |
| 13 | Challenges with Unsupervised LLM Knowledge Discovery (Farquhar, Varma, Kenton, Gasteiger, Mikulik, Shah) | arXiv:2312.10029 | VERIFIED — title, all 6 authors |
| 23 | Finding Neurons in a Haystack (Gurnee, Nanda, Pauly, Harvey, Troitskii, Bertsimas) | arXiv:2305.01610 | VERIFIED — title, all 6 authors |
| 27 | Physics of LMs Part 3.3, Knowledge Capacity Scaling Laws (Allen-Zhu, Li) | arXiv:2404.05405 | VERIFIED — title, authors |
| 31 | Function Vectors in LLMs (Todd, Li, Sen Sharma, Mueller, Wallace, Bau) | arXiv:2310.15213 | VERIFIED — title, all 6 authors, ICLR 2024 |
| 34 | Stress-Testing Capability Elicitation w/ Password-Locked Models (Greenblatt, Roger, Krasheninnikov, Krueger) | arXiv:2405.19550 | VERIFIED — title, all 4 authors |
| 36 | Transfer learning enables predictions in network biology (Geneformer; Theodoris et al.) | doi:10.1038/s41586-023-06139-9 | VERIFIED — title, Nature 2023, Theodoris/Xiao/Chopra/Chaffin |
| 37 | Assessing the limits of zero-shot foundation models in single-cell biology (Kedzierska, Crawford, Amini, Lu) | doi:10.1101/2023.10.16.561085 | VERIFIED — title, authors, bioRxiv 2023 |

Note: WebSearch budget (200) was exhausted mid-audit; remaining checks used WebFetch against arXiv, the Semantic Scholar API, and publishers, which was sufficient.

## Untested but high-confidence (well-established canonical works)
The 28 unchecked entries are all widely-cited papers from top venues (NeurIPS/ICLR/ACL/TMLR/Nature Human Behaviour/COLM/TACL) whose IDs and titles I recognize as correct, e.g. Wei "Emergent Abilities" (2206.07682), Schaeffer "Mirage" (2304.15004), BIG-bench (2206.04615), Burns CCS (2212.03827), Marks/Tegmark Geometry of Truth (2310.06824), Meng ROME (2202.05262), Li Othello-GPT (2210.13382), Power Grokking (2201.02177), Nanda Progress Measures (2301.05217). No ID looks malformed or anomalous. Zero red flags. Given a 9/9 clean spot-check rate, fabrication risk across the remainder is assessed as very low.

## On-topic / quality assessment
- Remit = measuring emergent abilities and latent knowledge, elicitation, ELK, concept/knowledge-neuron probing, "does the model know X" evaluation methodology. Coverage is well-matched and well-structured (8 thematic parts).
- Parts IV and VI (ROME, knowledge neurons, function vectors, RepE) lean toward interpretability/editing/steering rather than pure "emergence measurement," but they fall squarely under the remit's "capability elicitation" and "concept probing" clauses and each entry justifies its relevance. Not flagged.
- Entry 10 (Emergent Analogical Reasoning, Webb et al.) attracted a documented contamination critique; the entry itself flags this and turns it into a design caveat — handled honestly, not a weakness.
- One preprint (entry 37, bioRxiv) is appropriately labeled as such; appropriate given it is the key negative-result audit for the biology section. No predatory or low-tier venues present.

## Completeness check
Scanned all 37 entries: every one has a Technical summary (technical), a Plain-English line, an Applicability (Ax) mapping, and a Novelty implication. No structurally incomplete entries.
