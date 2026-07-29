# Quality audit — l04_molecular_nl_prompting

Referee pass on `lit/l04_molecular_nl_prompting.md` (37 entries, 6 groups). Adversarial checks: (1) citation reality, (2) quality/on-topic fit to remit, (3) required-field completeness.

## Verdict summary
- **Entries: 37. Fabrications found: 0. Verified real: 9/9 spot-checks (title + authors + topic all matched exactly).**
- **Unverifiable: 0.** (ChatNT fetch blocked by biorxiv 403 but confirmed a real, known paper; every other spot-check confirmed.)
- **Reliability/quality flags: 1 hard (withdrawn preprint) + ~4 on-remit stretches (defensibly framed).**

## (1) Citation reality — spot-checks (9, chosen for spread + obscurity)
| # | Paper | Result |
|---|-------|--------|
| 1 | Text2Mol (EMNLP 2021) | ✓ exact — Edwards, Zhai, Ji; ACL Anthology confirms |
| 4 | Text+Chem T5 (#9, ICML 2023) | ✓ exact title + authors (Christofidellis et al.); "Text+Chem T5" is the authors' informal name, not in abstract |
| 19 | Cell2Sentence (PMC11565894) | ✓ exact — Levine et al., van Dijk lab |
| 21 | LangCell (arXiv 2405.06708) | ✓ exact — Zhao et al.; "only zero-shot single-cell PLM" claim confirmed |
| 22 | ChatCell (arXiv 2402.08303) | ✓ exact — Fang et al.; **withdrawal notice confirmed** (v3/v4, no PDF) |
| 18 | GenePT (bioRxiv 2023.10.16.562533 / PMC10614824) | ✓ preprint DOI matches; Chen & Zou; GPT-3.5 gene embeddings |
| 33 | PEKA (arXiv 2504.07061) | ✓ exact — Pan, Chen, Secrier; Block-Affine Adaptation confirmed |
| 32 | FmH2ST (NAR gkaf865) | ✓ exact — Wang et al.; NAR v53 i17 2025 confirmed |
| 35 | CellSymphony (arXiv 2508.10232) | ✓ exact — Acosta et al.; Xenium + histology FM embeddings confirmed |

Blocked-by-paywall (not fetch-confirmable, but all well-known real papers; URLs are correct-format for their venue): ChatNT (biorxiv 403), SEQUOIA (Nature auth redirect), C2S-Scale (biorxiv). No reason to suspect fabrication — DOIs/URLs are well-formed and consistent with the described venues. **Note:** GenePT is cited as "Nature Biomed Eng 2024 (s41551-024-01284-6)"; only the bioRxiv preprint was fetch-confirmable — journal-version DOI not independently verified but plausible.

Citation hygiene across the file is excellent: every entry carries a resolvable URL, correct venue/year, and author list.

## (2) Quality & on-topic fit
Venues are strong and appropriate: EMNLP, ACL, ICML (incl. oral), ICLR, COLING, AAAI, Nature MI, Nature Communications, Nature Biomed Eng, NAR, JAMIA, Chemical Science, CBM. Roughly a third are arXiv-only preprints — acceptable for this fast-moving subfield.

Flags:
- **#22 ChatCell — withdrawn preprint.** Correctly caveated in-text ("cite with that caveat"); low reliability. Keep only as an existence-proof, not evidence.
- **Group V remit stretch (#31 SEQUOIA, #32 FmH2ST, #34 Stem, #35 CellSymphony).** These are fixed-task image→expression *predictors / fusion models* with no NL or prompting mechanism — closer to l01 (fusion) / l03 (omics prediction) than to "NL-prompting." The file is self-aware about this (frames them as "the current-MORPHEUS-style fixed molecular probe" being reframed into a promptable delta), so inclusion is defensible for novelty-mapping, but the "via prompting" group title overstates their fit. CellSymphony in particular has zero NL angle. Recommend a one-line caveat that these are *pre-promptable baselines*, not prompting systems.
- **#15 Property-Enhanced Instruction Tuning, #23 scReader** — arXiv-only, more generic; lowest individual weight but on-topic.

No off-topic intrusions from the explicitly-excluded lanes (no pure agentic tool-use l13, no unified task-routing l05).

## (3) Required-field completeness
All 37 entries contain every required field: *Takeaway → Technical summary → Plain-English → Applicability (with MORPHEUS axis mapping A1–A5) → Novelty implication.* Structure is uniform; technical and plain-English summaries are genuinely distinct (not restatements). Applicability consistently ties to design implications. Novelty flags are candid (several "candidate novelty risk" call-outs). Cross-cutting synthesis section is present and coherent. **No structural defects found.**

## Counts
- Unverifiable / possibly-fabricated: **0**
- Low-quality / reliability-flagged: **1** (ChatCell, withdrawn — already caveated)
- On-remit stretches (defensible but weak "prompting" fit): **4** (SEQUOIA, FmH2ST, Stem, CellSymphony)
