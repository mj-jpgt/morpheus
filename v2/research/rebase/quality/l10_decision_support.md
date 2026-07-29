# Quality audit — l10_decision_support

Lane: Decision-support methods (retrospective only). Remit: UQ/calibration/selective prediction, conformal UQ for survival, active feature acquisition / next-best-test, multimodal prognosis, treatment-response / counterfactual decision support. 45 entries (44 unique + 1 cross-reference).

Referee date: 2026-07-29. Method: adversarial spot-check via WebFetch against arXiv/publisher pages (WebSearch budget was exhausted at session start, so verification used direct URL fetches).

## Verdict summary

- **Spot-checked: 9 entries. Verified real and accurately described: 9/9.** No fabrications detected.
- Entry structure is strong and consistent; all substantive entries carry technical summary + plain-English + applicability + novelty.
- Unverifiable/likely-fabricated: **0.** Low-quality/off-topic: **0** (a handful of minor issues below, none disqualifying).

## Papers verified (all confirmed REAL, titles/authors/content match)

1. **#1 Ghawami, "Good Rankings, Wrong Probabilities"** — arXiv 2604.04239, submitted 2026-04-05. Real. Abstract confirms the 166/290 fold-level calibration rejections, gating-vs-bilinear finding, and Platt-scaling claim exactly as summarized. (This is a very recent single-author preprint but genuine.)
2. **#16 Kobayashi et al., "Learning-To-Measure (L2M)"** — arXiv 2510.12624. Real. Authors match exactly. NOTE: file dates it "2024"; arXiv shows submitted **Oct 2025** (revised May 2026). Minor date error.
3. **#34 PRIME** — arXiv 2604.04999, submitted 2026-04-05. Real. First author Kai Yu et al.; prototype-driven missing-modality prognosis, TCGA C-index 0.653 — matches.
4. **#14 von Kleist et al., AFAPE time-varying** — arXiv 2312.01530, JMLR 2025. Real. Two-assumption framing + semi-offline RL estimators (direct/IPW/double-RL) confirmed verbatim.
5. **#35 BITES** — Bioinformatics 2022, Schrod…Altenbuchinger. Real. Treatment-specific Cox loss + balanced DNN confirmed.
6. **#12 Weighted Conformal Prediction for Survival under Covariate Shift** — arXiv 2512.03738 (Shin, Lee, Kang, Yonsei), Dec 2025. Real. Content matches.
7. **#31 DisPro** — arXiv 2503.01653, CVPR 2025 (Xu, Zhou, Zhao, Wang, Yang, Chen). Real. UniPro/MultiPro two-stage LLM prompting confirmed.
8. **#23 SurvPath** — arXiv 2304.06819, CVPR 2024 (Jaume, Vaidya, Chen, Williamson, Liang, Mahmood). Real. Pathway-token framing confirmed.
9. **#27 SurvPGC, "clinical information prompts integration"** — npj Digital Medicine (Hou, Zhang, Xie, Li, Qin). Real. NOTE: file dates it "2025"; the article is **Vol 9, Art 76 (2026)**. Minor date error.

## Flagged entries (minor — none warrant removal)

- **Date mislabels:** #16 labeled 2024 → actually Oct 2025; #27/#45 labeled 2025 → actually 2026. Content and IDs are correct; only the year tags are off. Recommend correcting.
- **#45 is a cross-reference to #27**, not a standalone entry. It intentionally omits the technical-summary and plain-English fields (carries only Takeaway/Applicability/Novelty). Strictly it does not satisfy the "must contain technical AND plain-English summary" rule, but this is by design as a cross-listing, so the unique-entry count is effectively 44.
- **Very recent 2026 preprints (#1, #34):** genuine but new and lightly cited; treat conclusions as provisional. Not a fabrication concern — both verified on arXiv.
- **Topical overlap (#2 and #8):** both cover selective prediction / abstention on clinical-text MIMIC tasks with similar setups; #8 is a lighter preprint. On-remit (UQ/selective prediction) but somewhat redundant and only loosely tied to cancer. Acceptable.
- **Access-gated entries** (#38 npj, #36 ResearchGate, several PMC/ScienceDirect) were not fetched; titles/venues are plausible and consistent with the lane. No red flags, but they are the least-verified subset.

## On-topic / quality assessment

All 45 entries fall within the stated remit and are framed as methods + retrospective evaluation (no deployed-CDS overreach). Axis mapping (A1–A5) is applied consistently and the novelty/pre-emption flags (SurvPath A2, DisPro A4, clinical-prompt A1, BITES/CFR-Net A5, in-context AFA A1/A5) are well-reasoned and honest about prior art. Mix of foundational anchors (CFR-Net 2017, Deep Sensing 2018, Conformal Risk Control) and current SOTA is appropriate. Cross-cutting notes are accurate and useful.
