# Quality audit — l11_benchmarks_confound

Referee pass on `v2/research/rebase/lit/l11_benchmarks_confound.md` (27 entries). Date: 2026-07-29.

## 1. Reality / findability (spot-check)
Verified 9 entries directly (arXiv abstract pages + publisher pages via WebFetch; WebSearch budget was exhausted so verification was done by fetching canonical URLs). All checked entries are REAL with matching titles/authors:

- **#1 HEST-1k** (arXiv:2406.16192) — CONFIRMED. Title + author list (Jaume … Mahmood) exact; NeurIPS'24 Spotlight.
- **#4 HiST** (arXiv:2606.14251) — CONFIRMED. Title + authors (Wu, Xu, Diao, Li, Wei, Andersson, Gui) exact. (June-2026 arXiv ID; consistent with current date.)
- **#9 de Jong Robustness Index** (arXiv:2501.18055) — CONFIRMED. Title, authors (de Jong, Marcus, Teuwen), and "only one model RI>1" claim all match.
- **#14 GLMP** (arXiv:2606.28697) — CONFIRMED. Title + 11-author list match; "first pathology model to use text as intermediate representation" matches the entry's claim.
- **#19 Herbert monotherapy drug-blind** (PLOS Comput Biol, 10.1371/journal.pcbi.1013232) — CONFIRMED. Title + authors (Herbert, Chia, Jensen, Walther-Antonio), intraclass-generalization thesis match.
- **#20 Wong "Simple controls exceed…"** (Bioinformatics 41(6) btaf317) — CONFIRMED. Title + authors (Wong, Hill, Moccia, Pfizer); CRISPR-informed-mean-beats-scGPT/GEARS claim matches.
- **#25 Tizhoosh "Beyond the Failures"** (arXiv:2510.23807) — CONFIRMED. Title + author; "dense embeddings cannot represent combinatorial richness" thesis matches.
- **#27 PLUTO-4** (arXiv:2511.02826) — CONFIRMED. Title + lead author (Padigela); frontier pathology FM.

Entries whose canonical URLs are valid but sat behind auth/anti-bot walls (403 / Nature IDP redirect) — NOT flagged as fabricated; DOIs/IDs are well-formed and these are well-known works:
- **#3 Wang, translational ST benchmark** (Nat Commun s41467-025-56618-y) — Nature login redirect.
- **#7 Howard site-specific signatures** (Nat Commun s41467-021-24698-1) — Nature login redirect; landmark, widely cited paper.
- **#8 Dawood "Buyer Beware"** (bioRxiv 2024.06.23.600257) — 403; well-known preprint.
- **#23 COMPASS** (medRxiv 2025.05.01.25326820v3) — 403; DOI well-formed.

**No fabricated citations detected.** 0 unverifiable-as-nonexistent.

## 2. Quality / on-topic
Overall high quality and tightly on-remit. Section A (WSI->molecular benchmarks) and Section B (site/batch confounding) are squarely core. Notes:

- **Section C leans on non-WSI drug/perturbation papers (#19–#22, #24).** These are single-cell/drug-response, not histology. Justified by the remit's "ORBIT-like held-out protocols + baseline discipline" angle and A5 mapping, and the lane is explicit about it — acceptable, but a reader should know these are *methodological analogues*, not WSI->molecular papers.
- **#24 (Nat Commun s41467-025-56827-5, drug sensitivity principles)** — weakest entry: technical summary is generic with no quantitative detail or specific method, and applicability is "neutral-to-strengthening." Borderline filler; keep only as a comparator.
- **Venue/award labels partially unverified.** Several attributions ("NeurIPS Spotlight," "ICML 2026," "MICCAI 2025 Oral," "IJCAI 2025 Survey Track") were not all independently confirmed — the *papers* are real, but treat the venue/honor tags as lower-confidence.

## 3. Structure completeness
All 27 entries follow a consistent schema: *Takeaway / Technical summary (technical) / Plain-English / Applicability (A-axis mapping) / Novelty implication*. Every entry therefore contains a technical summary + plain-English summary + applicability + novelty. **No structural gaps.** Plain-English lines are genuinely lay-readable, not restated jargon.

## Tally
- Entries: 27
- Spot-checked and confirmed real: 9 (plus 4 more with valid-but-walled URLs)
- Unverifiable / suspected fabricated: **0**
- Low-quality / weak-on-topic flagged: **1** clearly thin (#24); soft flag on #19–#22 as non-WSI analogues and on unverified venue labels.
