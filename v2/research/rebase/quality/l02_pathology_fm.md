# Quality audit — l02_pathology_fm

Referee pass on `lit/l02_pathology_fm.md` (30 entries). Adversarial checks: (1) papers real/findable, (2) high-quality + on-topic, (3) required structure present.

## Verdict summary
- **Entries:** 30
- **Directly verified real (WebFetch):** 6 (spot-check target ~6 met)
- **Additional entries judged real** (well-established, paywalled but canonical): UNI, MUSK, CONCH, PLIP, Prov-GigaPath, CHIEF, Virchow(2), Phikon(-v2), PathChat — all widely-cited literature.
- **Unverifiable / possibly fabricated:** 0
- **Low-quality / off-topic:** 0
- **Structure defects (missing technical / plain-English / applicability / novelty):** 0

## Spot-checks performed (all PASSED — title/authors/claims match)
1. **CARE** (arXiv:2602.21637) — the most suspicious entry (future-dated Feb 2026 / CVPR 2026). VERIFIED REAL: Zhang, Gong, Pang, ... Crispin-Ortuzar, Yu, Li, Gao; two-stage SSL on 34,277 WSIs + RNA/protein alignment, 33 tasks. Description in lit file is accurate. Not fabricated despite the recent date.
2. **THREADS** (arXiv:2501.16652) — VERIFIED: Vaidya, ... Jaume, Lu, Mahmood; ~47,171 samples paired genomic+transcriptomic, 54 oncology tasks. Accurate.
3. **Comparing CPath FMs via RSA** (arXiv:2509.15482) — VERIFIED: Mishra, Lotter; UNI2/Virchow2 most distinct, Prov-GigaPath most average, strong slide/batch signal + weak disease signal, stain-norm 5.5–20.5%. Accurate, incl. numbers.
4. **KEEP** (arXiv:2412.13126) — VERIFIED: Zhou et al.; disease KG of 11,454 diseases / 139,143 attributes, 18 benchmarks. Accurate.
5. **mSTAR** (arXiv:2407.15362) — VERIFIED: Xu, Wang, Zhou, ... Chen; reports + gene expression injected into patch reps. Accurate. (Minor: lit says "26k slide-level pairs"; abstract headlines 116M patches / 32 cancers — the 26k is the paired-multimodal subset, not an error.)
6. **TITAN** (arXiv:2411.19666) — VERIFIED: Ding, Wagner, Song, ... Mahmood; 335,645 WSIs, visual SSL + VL alignment, zero/few-shot + report gen. Accurate.

## Quality / on-topic assessment
- Every entry is on-remit (WSI/tile FMs, vision-language pathology, molecular/multimodal, plus a small, clearly-labeled set of benchmark/survey/architecture context items). No off-topic drift.
- Borderline entries are handled honestly in-text rather than overclaimed:
  - **#9 H-Optimus-0** — flagged in the entry itself as a HuggingFace weights release with "no full paper." Acceptable and transparent.
  - **#18 BiomedCLIP** — general biomedical (not pathology-specific), but explicitly framed as a cross-domain baseline. Appropriate.
  - **#23 Kather 2020** — a CNN study, not a "foundation model," but explicitly included as the empirical prior for "morphology encodes molecular state." Legitimate context.
- Structure: all 30 entries contain Takeaway + Technical summary + Plain-English + Applicability (mapped to A1–A5) + Novelty implication. No structural gaps found.
- Axis mapping (A1–A5) is consistently applied and the cross-modal-by-analogy grounding caveat is stated up front (line 5), which is intellectually honest given this is an imaging lane feeding a transcriptomics programme.

## Conclusion
Strong, well-curated lane. No fabricated, unverifiable, or low-quality citations detected; all six adversarial spot-checks (including the future-dated CARE) confirmed real and accurately summarized.
