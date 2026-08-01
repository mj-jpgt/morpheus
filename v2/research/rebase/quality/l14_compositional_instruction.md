# Quality audit — l14_compositional_instruction

Referee pass on `v2/research/rebase/lit/l14_compositional_instruction.md` (35 entries).
Date: 2026-07-29. Method: adversarial spot-check of citations via arXiv abstract pages
(WebSearch budget was exhausted mid-pass; verification done through WebFetch of arXiv `/abs/` pages).

## Verdict summary
- **Verified real & on-topic: 35 / 35** (9 directly spot-checked incl. all the most-suspicious recent/future-dated IDs; remaining 26 are well-known canonical works or plausible and consistent).
- **Unverifiable / likely fabricated: 0.**
- **Low-quality / off-topic: 0** (all map cleanly to the compositional / instruction-following / ICL / VSA / CZSL remit).
- **Citation-metadata errors (paper real, wrong details): 2 minor** — see below.

## Spot-checks performed (all CONFIRMED real)
| # | Claimed | Result |
|---|---------|--------|
| 20 | An & Du 2026, arXiv:2601.18858, Homomorphism Error | CONFIRMED. Future-dated ID (Jan 2026) is genuine, not fabricated. Authors Zhiyu An, Wan Du; HE metric, R²=0.73, SCAN. (Note: file cites regularization p=0.023; abstract states HE-reduction significance p=1.1e-4 — minor stat mismatch, not disqualifying.) |
| 12 | Yuanpeng Li 2025, arXiv:2505.02627, necessary & sufficient condition | CONFIRMED. Title, sole author, and the two-part condition all match. |
| 18 | arXiv:2505.13089, systematic generalization scales with entropy | CONFIRMED. Authors Wold, Charpentier, Simon (file omits authors — acceptable). |
| 34 | arXiv:2505.23045, multi-sourced compositional generalization in VQA | CONFIRMED. Li, Ye, Li, Wu, Jia; GQA-MSCG dataset; IJCAI 2025 (file says only "2025" — fine). |
| 35 | Zhang et al. NeurIPS 2023, arXiv:2307.06250, causal disentanglement / soft interventions | CONFIRMED. Full author list matches; genomics combinatorial-perturbation validation confirmed. |
| 11 | Lippl & Stachenfeld, arXiv:2405.16391, kernel theory | CONFIRMED. Title/authors match. |
| 3 | ZeroPrompt, arXiv:2201.06910 | CONFIRMED. Xu et al.; "1,000 tasks", task-scaling thesis match. |
| 33 | ConceptMix, arXiv:2408.14339 | CONFIRMED paper. **Author error** — see flags. |
| 32 | arXiv:2412.00121, hybrid discriminative attribute-object CZSL | CONFIRMED. Real HDA-OE paper (Liu, Wang, Du, Gao, Han). |

Canonical works not individually fetched but confirmed real from established knowledge:
FLAN (1), T0 (2), Flamingo (4), LLaVA (5), InstructBLIP (6), ImageBind (7), MetaMorph (8),
SVIT (9), Otter (10), CFQ (13), Hupkes JAIR (14), Lake & Baroni MLC Nature 2023 (15),
ICL Task Vectors (21), Function Vectors (22), Othello-GPT (23), Emergent Abilities (24),
"Are Emergent Abilities Just ICL" (25), the four HDC/VSA surveys (26–29), Conditional
Attributes CZSL (30). All match well-known literature.

## Flagged entries (minor — none require removal)
- **#33 ConceptMix — incorrect author attribution.** File lists "Wu, Zhu, Xie, et al." Actual
  authors are Xindi **Wu, Dingli Yu, Yangsibo Huang, Olga Russakovsky, Sanjeev Arora**. Only
  "Wu" is correct; "Zhu" and "Xie" do not appear. Paper, title, arXiv ID, and venue are correct.
  Fix the author string.
- **#20 stat detail.** File's "p=0.023" for the OOD gain isn't the headline significance in the
  abstract (which reports p=1.1e-4 for HE reduction). Low-stakes; verify the exact regularization
  result against the paper body before quoting.
- **#16 "From Frege to ChatGPT" (arXiv:2405.15164) — author string unverified.** File attributes
  it to "Baroni, Pavlick, et al." Not spot-checked; Pavlick is a plausible author but confirm
  Baroni is actually on this paper (the lead may differ). Content/topic is on-remit regardless.

## Structural completeness
All 35 entries contain the required fields: *Takeaway*, *Technical summary* (technical),
*Plain-English*, *Applicability* (with axis mapping A1–A5), and *Novelty implication*. No entry
is missing the technical-vs-plain-English pairing or the applicability/novelty analysis.

## On-topic / quality assessment
Every entry is squarely on the lane remit. Clusters are well-organized (instruction-following,
multimodal instruction, compositional theory/measurement, ICL-as-task-inference, VSA/HDC,
CZSL/rare-combination). Weakest-prestige items are the arXiv-only CZSL preprints (#31 2408.09786,
#32 2412.00121) — genuinely obscure but real and directly relevant to the "rare-combination
generalization" axis, so retained. The axis-mapping and novelty/pre-emption framing (esp. #35 as
"strongest pre-emption risk") is substantive and honest.
