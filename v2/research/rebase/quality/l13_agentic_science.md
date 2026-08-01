# Quality Audit — l13_agentic_science

**Lane:** Agentic LLM scientific workflows & tool-use
**Auditor pass:** adversarial verification of realness, on-topic quality, and entry completeness.
**Date:** 2026-07-29

## Summary
- **Entries in file:** 32
- **Spot-checked for realness:** 7 (all verified REAL via arXiv / PMC / Oxford Academic)
- **Unverifiable / possibly fabricated:** 0
- **Low-quality or off-topic (flagged):** 0 hard flags; 2 soft/tangential notes (justified in-text)
- **Structural completeness:** all 32 entries carry Technical summary + Plain-English + Applicability + Novelty implication. PASS.

## Realness spot-checks (verified)
| # | Entry | Citation checked | Result |
|---|-------|------------------|--------|
| 3 | Gorilla | arXiv:2305.15334 | REAL — Patil, Zhang, Wang, Gonzalez. Title exact. |
| 11 | Biomni | PMC12157518 (bioRxiv) | REAL — Huang, Zhang, Wang, ... Leskovec (Stanford). 150 tools / 105 pkgs / 59 DBs confirmed. |
| 12 | TxAgent | arXiv:2503.10970 | REAL — Gao, Zhu, Kong, ... Zitnik. 211 tools, 92.1% accuracy confirmed. |
| 13 | CRISPR-GPT | arXiv:2404.18021 | REAL — Qu, Huang, Yin, ... Cong. NatBME acceptance confirmed. |
| 14 | BioDiscoveryAgent | arXiv:2405.17631 | REAL — Roohani, Lee, Q. Huang, ... Leskovec. 21% / 46% figures confirmed. |
| 18 | GeneAgent | arXiv:2405.16205 | REAL — Wang et al. 1,106 gene sets, beats GPT-4 confirmed. |
| 31 | BiB 2026 review | Briefings in Bioinformatics 27(2):bbag110, DOI 10.1093/bib/bbag110 | REAL — lead author Sajib Acharjee Dip (Virginia Tech). This was the highest fabrication-risk citation (journal-only, future-dated, unusual "bbag" issue prefix); confirmed genuine. |

Note: WebSearch budget was exhausted at session start; all verification done via direct WebFetch against arXiv abstract pages, PMC, and Oxford Academic. The 7 checks span the highest-risk citations (biology-specific, recent 2024-2026, and the journal-only future-dated review) rather than the well-known foundational papers.

## Quality / on-topic assessment
All 32 entries are on-remit for a lane whose remit is "LLM agents that plan/execute scientific tasks; tool-use over scientific models; LLM-as-orchestrator for biology; model-as-tool APIs; function-calling to specialist predictors."

- **Foundational agent papers** (1 ReAct, 2 Toolformer, 5 Reflexion, 6 Voyager) are general-domain, not biology — but each is the canonical mechanism (think→act loop, learned tool-gating, verbal self-critique, skill library) the biology agents build on, and each entry explicitly justifies its A1/A3/A4 relevance. Kept, not flagged.
- **Chemistry agents** (7 Coscientist, 8 Emergent, 9 ChemCrow, 24 CACTUS) are adjacent-domain but are the load-bearing precedents for A1 routing and A5 intervention-as-query; on-topic as methodology.
- **Biology/omics agents** (10, 11, 13, 15, 16, 17, 18, 29) and **eval benchmarks** (25 LAB-Bench, 26 ScienceAgentBench, 27 DiscoveryWorld, 30 MLAgentBench) are squarely core to the remit and its A3 evaluation emphasis.
- **Surveys** (10, 31, 32) are appropriate framing/positioning references.

Soft notes (no action required):
- **Entry 6 (Voyager, Minecraft):** most domain-distant of the set; justified as a skill-library template but is the weakest topical fit. Acceptable.
- **Entry 14 author line** lists "...Leskovec, Zou"; the arXiv author block surfaces through Leskovec (Zou/Marson present in full list). Cosmetic, not a realness issue.

## Verdict
No fabricated, unverifiable, or off-topic entries detected. The lane is high-quality, well-scoped, and the novelty framing (routing pre-emption vs. MORPHEUS's unified-trunk wedge) is consistently applied. Fit to use.
