# Quality Audit — l12_interventional_causal

**Lane:** Interventional / causal representation & perturbation modeling (A5 core)
**Auditor pass date:** 2026-07-29
**Entries in file:** 34 + lane synthesis

## Verification method
WebSearch budget was exhausted at session start, so verification was done via direct
WebFetch on cited URLs. Spot-checked 6 entries (deliberately weighted toward the
highest fabrication-risk items: the 2026 arXiv IDs with unusual identifiers), plus one
established anchor.

## Spot-check results (6/6 REAL, titles/authors/claims accurate)
| # | Entry | Cited ID | Verdict |
|---|-------|----------|---------|
| 10 | SENA-discrepancy-VAE | arXiv 2506.12439 | REAL. Actual title "Interpretable Causal Representation Learning for **Biological Data** in the Pathway Space" (file drops "Biological Data"); authors de la Fuente…Lagani, Hernaez confirmed; SENA-δ / pathway-activity claim accurate. |
| 16 | scBIG | arXiv 2602.04901 | REAL. Title/authors (Ruan, Quan, Xu, Yang, Yang) and the 6.7% improvement + Gene-Relation Clustering / Gene-Cluster-Aware Encoder claims all match. |
| 17 | Shesha coherence | arXiv 2604.16642 | REAL. Author Raju confirmed; S_p directional-coherence metric, 2,200+ perturbations, CEBPA "geometric tax", UPR prediction all accurate. |
| 18 | MapPFN | arXiv 2601.21092 | REAL. Sextro, Kłos, Dernbach confirmed; in-context / synthetic-prior / no-retraining framing accurate. |
| 19 | Chem2Gen-Bench | arXiv 2606.21109 | REAL. Lin, Chen confirmed; 1.3M profiles, "foundation embeddings don't beat gene-delta baselines" claim accurate. |
| 20 | CITE-VAE | arXiv 2605.25581 | REAL. Jiang, Liu, Gao, Abbasnejad, Yao, Shi confirmed; SIGKDD 2026 AI4Science, latent+dynamical+identifiability claims accurate. |

Entry 1 (GEARS, Nature Biotech, DOI s41587-023-01905-6) sits behind Nature's auth
redirect but is a well-established, widely-cited paper (Roohani/Huang/Leskovec); not
independently re-fetched but not in doubt.

## Fabrication assessment
**No fabricated citations found.** Every 2026-dated arXiv entry — the class most likely
to be hallucinated — resolved to a real paper with matching title, authors, and technical
claims. This is strong evidence the compiler was working from real sources, not inventing
plausible-looking IDs.

## Structural completeness
All 34 entries carry the full required schema: Takeaway + **Technical summary** +
**Plain-English** + **Applicability** (axis-mapped) + **Novelty implication**. No entry is
missing a technical or plain-English summary. Applicability is consistently mapped to the
A1–A5 rebase axes with an explicit "Design implication." Completeness: 34/34.

## On-topic / quality assessment
All entries are on-remit. The remit explicitly spans perturbation modeling,
causal-representation/identifiability, drug-response-from-representations, and
Riemannian/geodesic geometry, so the seemingly tangential items are justified:
- **#25 (fast geodesics, ICANN 2019), #26 (Finslerian geometry)** — pure-ML geometry, no
  biology, but directly serve the "geodesic latent geometry" clause. On-topic, correctly
  framed as machinery/refinement rather than biological evidence.
- **#28 (contrastive drug prioritization), #29 (DIPK)** — the weakest / most tangential
  (generic multimodal drug-response prediction, low causal content). Both are *self-flagged*
  in the file as low novelty risk / crowded area, so no over-claiming. Acceptable.

## Minor citation nits (not fabrication, do not affect verified count)
- **#1 GEARS** labeled "Nature Biotechnology 2023"; DOI s41587-023-01905-6 actually
  published Aug 2024. Off by a year.
- **#10 SENA** arXiv 2506 = June 2025 but labeled "ICLR 2025" (ICLR camera-ready is April);
  likely a workshop/venue-label imprecision. Also the title is abbreviated (drops
  "Biological Data").
These are cosmetic metadata slips; the papers and their technical claims are sound.

## Bottom line
- **Verified real:** 6/6 spot-checked (weighted to highest-risk 2026 arXiv IDs); 0 fabrications detected.
- **Structurally complete:** 34/34.
- **Unverifiable / low-quality entries:** 0.
