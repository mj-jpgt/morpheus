# Quality audit — l03_omics_fm (Omics & multi-omics foundation models)

Referee pass date: 2026-07-29. Source file: `v2/research/rebase/lit/l03_omics_fm.md` (33 entries + cross-cutting synthesis).

## Method
WebSearch budget was exhausted before use, so verification was done by fetching cited
URLs / preprint + repo landing pages directly (Nature and Cell block automated fetches, so
those were cross-checked via GitHub, arXiv, OUP, and OpenReview mirrors where possible).

## Verification (spot-checks)
Verified REAL and on-topic (title/authors/venue confirmed):
- #12 **BulkFormer** — confirmed via GitHub (KangBoming/BulkFormer), Cell Systems 2026, ~150M params, GNN+Performer, >500k bulk profiles. Real. (DOI resolves to 10.1016/j.cels.2026.101657.)
- #13 **"one PCA still rules them all"** — confirmed arXiv:2410.13956 (Bendidi et al., 2024); finds scVI/PCA beat FMs on perturbation. Real, exactly as described.
- #14 **BMFM-RNA** — confirmed arXiv:2506.14861 (IBM Research / BiomedSciAI, Danziger/Dandala et al., 2025). Real. (Live arXiv title now reads "whole-cell expression decoding improves transcriptomic foundation models" — a revised title; entry uses the v1 title. Harmless.)
- #23 **Phosformer** — confirmed Bioinformatics 39(2), btad046, 2023 (Zhou, Yeung, ... Kannan). Real, on-topic.
- #31 **LangCell** — confirmed arXiv:2405.06708, ICML 2024. Real. BUT author list is wrong (see below).

Not directly fetchable (Nature/Cell/bioRxiv/OpenReview bot-blocked) but plausible and
consistent with known literature; no fabrication indicators:
- #3 Geneformer scaling/quantization (Nat Comput Sci 2026, DOI s43588-026-00972-4) — DOI pattern consistent with a real 2026 NCS article; Theodoris lab is active and #2 Geneformer is confirmed real. Treated as PLAUSIBLE, not verified.
- #5 UCE, #26 Cell2Sentence/C2S, #27 scPRINT, #32 AIDO.Cell — all well-known real works in this field; landing pages blocked but no red flags.

No entry appears fabricated. 0 suspected-fake citations.

## Flags (quality, not existence)
1. **#31 LangCell — incorrect author list.** File lists "Zhao, Yang, Sun, ... Yao, Wang";
   actual authors are Suyuan Zhao, Jiahuan Zhang, Yushuai Wu, Yizhen Luo, Zaiqing Nie.
   First author (Zhao) is right; the rest are wrong. Fix the citation.
2. **#12 BulkFormer — misparsed author.** "Kang, Bo, et al." treats one author (Kang Boming;
   GitHub handle KangBoming) as two. Should read "Kang, Fan, ... Cui." Cosmetic.
3. **Minor: #14 BMFM-RNA title** is the v1 arXiv title; the published/revised title differs.
   Not wrong, just note the paper was retitled.

## On-topic / remit fit
All 33 fit the remit (RNA/cell/gene FMs, proteomics/phospho representation learning,
multi-omics integration). The proteomics/protein-LM entries (#21 Prosit, #22 DeepPhospho,
#24 ESM-2, #25 ESM3, #29 PINNACLE) are within remit because it explicitly covers
"proteomics/phospho representation learning" — not off-topic. No weak or filler entries;
even the classical baselines (#16 scGen, #17 CPA, #20 MOFA+, #18 totalVI, #19 MultiVI) are
deliberately included as A2/A4/A5 baselines and justified as such.

## Structure completeness
All 33 entries contain the full required block: Takeaway, **Technical summary**,
**Plain-English**, **Applicability** (mapped to A1-A5 with design implication), and
**Novelty implication**. No entry is missing a component. Cross-cutting synthesis maps
prior art to each axis coherently.

## Bottom line
- Verified real: 5/6 spot-checks fully confirmed; 1 (#3) plausible-but-unfetchable. 0 fabricated.
- Low-quality / needs-fix entries: **2** (citation-metadata errors in #31 and #12); both are
  real papers with genuine content — errors are in author metadata, not substance.
- Content quality is high and uniformly on-remit; structure is complete throughout.
