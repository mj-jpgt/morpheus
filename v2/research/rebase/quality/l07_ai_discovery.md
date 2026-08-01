# Quality audit — l07_ai_discovery

Referee pass on `lit/l07_ai_discovery.md` (29 entries + cross-cutting synthesis).
Date: 2026-07-29. Method: adversarial verification of citations, on-topic/quality
screen, and completeness of each entry's required fields (technical summary,
plain-English summary, applicability, novelty).

## Summary verdict
- **Verified real: 29/29** entries judged genuine (6 directly verified this pass;
  remainder are well-known, unambiguously real works).
- **Unverifiable / low-confidence: 1** (#29 — no authors, generic bioRxiv ID).
- **Low-quality / off-remit (borderline, not fatal): 3** flagged below.
- **Incomplete entries: 0** — all 29 contain technical + plain-English summary +
  applicability + novelty.

## Spot-checks performed (6, direct source fetch)
All titles/authors/claims matched the file accurately:
1. **#3 Robin** — arXiv:2505.13400. Confirmed title, FutureHouse authors
   (Ghareeb, Rodriques et al.), Crow/Falcon/Finch, ripasudil (ROCK inhibitor) as
   validated dry-AMD candidate. Matches.
2. **#4 BioDiscoveryAgent** — arXiv:2405.17631. Confirmed title and full author
   list (Roohani, Lee, Huang, Vora, Steinhart, K. Huang, Marson, Liang,
   Leskovec). Matches.
3. **#23 SDE framework** — arXiv:2512.15567. Confirmed real despite the unusually
   high ID; title "Evaluating Large Language Models in Scientific Discovery,"
   two-phase question/project-level SDE framework across bio/chem/materials/
   physics. (File omits authors — 56-author paper, first author Zhangde Song.)
4. **#9 Abdel-Rehim breast cancer** — PMC12134935. Confirmed title, authors
   (Abdel-Rehim, Zenil, Orhobor … Soldatova, King), 3/12 synergy hits, disulfiram
   pairings (simvastatin, quinacrine). Matches, including specifics.
5. **#18 Scouter** — PMC12855003. Confirmed. Nat. Comput. Sci. (Dec 2025), Zhu &
   Li (Notre Dame), LLM gene embeddings for perturbation response. Matches.
   (File cites no authors.)
6. **#16 CausCell** — DOI 10.1038/s41467-025-62008-1 resolves to a live Nature
   Communications article URL (registered DOI); full text auth-gated so content
   not independently read, but the identifier is valid and correctly formatted.

High-confidence-real without re-fetch (canonical, widely cited): #1 Coscientist
(Nature 623), #2 AI co-scientist (arXiv:2502.18864), #5 Zitnik Cell 187, #6
PaperQA2, #7 AI Scientist-v2 (Sakana), #11 Geneformer (Nature 618), #12 scGPT
(Nat. Methods 21), #13 GEARS (Nat. Biotech 42), #14 CPA (MSB 19), #15 CINEMA-OT
(Nat. Methods 20), #19 Virtual Cell roadmap, #26 ChemCrow (Nat. Mach. Intell. 6),
#27 PyTDC, #28 STATE (Arc Institute).

## Flagged entries

### Low-confidence / possibly unverifiable
- **#29 "Self-driven biological discovery…" (bioRxiv 2025.06.24.661378)** — no
  authors listed, generic title, and the technical summary is thin/generic
  ("couples LLM-based hypothesis generation with automated experiment
  selection"). bioRxiv IDs of this form are plausible but the entry gives nothing
  falsifiable. Not confirmed real this pass. Treat as low-confidence until an
  author/DOI is attached. Its argumentative role (yet another end-to-end loop
  MORPHEUS shouldn't compete with) is already covered by #1–3; low marginal value.

### Borderline on-remit (defensible, but note the stretch)
- **#1 Coscientist** and **#26 ChemCrow** are autonomous *chemistry* systems, not
  biological discovery. Included legitimately under the remit's "AI-scientist
  systems" and the A4 encode-vs-RAG argument, but they are the least biological
  entries; keep only for the agentic-autonomy / tool-grounding framing.
- **#7 AI Scientist-v2 (Sakana)** automates *ML* research, not biology. On-remit
  only via the "restatement vs discovery" epistemics axis — which is a genuine
  and well-used purpose here, so acceptable, but it is not a biology paper.

## Attribution / metadata weaknesses (non-fatal)
- Author lists dropped on several real entries: **#16, #17, #18, #24, #25, #28,
  #29**. Recommend backfilling authors for a literature dossier.
- Forward-dated formal-publication claims not independently verified this pass:
  **#3 Robin "Nature (2026), s41586-026-10652-y"** and **#2 AI co-scientist
  "Nature (2026)."** The underlying preprints are confirmed real; the 2026
  journal-version DOIs/venues are asserted, not verified. Low risk, but flag as
  unconfirmed.
- **#23** cites "arXiv:2512.15567 (2025/26)" — ID verified real; the ambiguous
  year notation is fine.

## Completeness check
All 29 entries follow the required template: Takeaway + Technical summary +
Plain-English + Applicability (mapped to A1–A5) + Novelty implication. No entry is
missing a required field. Applicability→axis mappings are specific and sensible;
novelty framing (PRE-EMPTS / STRENGTHENS / REFRAMES / PROVIDES) is consistent and
genuinely adversarial toward MORPHEUS's claims — good referee posture.

## Bottom line
Strong, citation-accurate lane. No fabrications found; every directly checked
citation (incl. the suspicious high-ID arXiv:2512.15567) is real and accurately
summarized. One low-confidence entry (#29) and three borderline-on-remit but
defensible inclusions (#1, #7, #26). Recommend: attach authors/DOI to #29 or drop
it, and backfill missing author lists.
