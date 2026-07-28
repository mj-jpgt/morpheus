# MORPHEUS Rank-Collapse Paper — Build & Assembly

This directory holds the source for *The Effective-Rank Fingerprint: Diagnosing
Cohort-Confounded WSI→Molecular Alignment from a Single Trained Model*.

## Branch

Work is on **`paper/rank-collapse-diagnostic`** (branched off `main`).

## Layout

The paper is authored as one Markdown file per section, then stitched into a
single `main.md`:

| File | Section in `main.md` |
|---|---|
| `01_abstract_intro.md` | Abstract + §1 Introduction |
| `02_related_work.md` | §2 Related Work |
| `03_method.md` | §3 Method (architecture + the fix) |
| `04_metrics.md` | §4 Metrics (confound-aware evaluation) |
| `05_experiments_results.md` | §5 Experiments & Results |
| `06_discussion_limitations.md` | §6 Discussion, Risks, Limitations |
| `references.bib` | Bibliography (BibTeX) |
| **`main.md`** | **Assembled single-file draft (generated)** |
| `README.md` | This file |

`main.md` is the stitched deliverable: title, author placeholder, abstract, then
the six sections in order (intro → related work → method → metrics → experiments
→ discussion), followed by a References note pointing at `references.bib`.

## Tables status

Empirical tables carry **verified numbers** and must not be edited when
re-stitching. Completeness:

| Table | Content | Status |
|---|---|---|
| T1 | Dual-head effective-rank spectrum (identity 84.3 vs biology 5.3–6.0) | **complete** |
| T2 | Confound decomposition ladder (global → within-cancer, −46 to −49%) | **complete** |
| T3 | Method-invariance of control-adjusted specificity (~+0.07) + head-to-head | **complete** |
| Retrieval | Recall@k, baseline leads | **complete** |
| **T4** | VICReg decorrelation ablation (rank recovery vs. specificity) | **[queued: Lambda]** — pre-registered hypothesis only, no numbers reported |

T4 is a pre-registered ablation whose full multi-seed λ run is queued on the
Lambda cluster; the draft states the hypothesis and the exact metric it will
report, and deliberately reports **no** T4 numbers.

## Verified numbers (do not alter without re-deriving)

- Effective rank: `wsi_identity` 84.3, `rna_identity` 37.5, `wsi_biology` 6.0,
  `rna_biology` 4.4, `full_biology` 5.3, `full_patient` 8.0.
- Modality gaps: identity 0.296, biology 0.475.
- Anchor: residual_scale −0.0011, correction_norm 0.00073, gate_mean 0.646.
- Held-out split: n = 2530, d = 256. Seeds {42, 43, 44}; 180 Hallmark targets.
- Global → within-cancer drop: ~46–49%. Random-gene null ~0.30–0.32 pooled /
  ~0.15 within. Method-invariant control-adjusted specificity ~+0.07.

## Regenerate `main.md`

`main.md` is produced by concatenating the section files in order with the
title/author/abstract front matter and the References note. To rebuild it after
editing a section file, re-stitch in this order:

```
01_abstract_intro.md → 02_related_work.md → 03_method.md →
04_metrics.md → 05_experiments_results.md → 06_discussion_limitations.md
```

then append the References note pointing to `references.bib`. Keep all verified
numbers and the **[queued: Lambda]** marker on T4 intact. When T4 completes,
drop the Lambda run's numbers into §5.5 and flip its status here to *complete*.

## Render to PDF / HTML

`main.md` uses `\citep{...}` keys and LaTeX math. Render with a citation-aware
toolchain, e.g.:

```
pandoc --citeproc main.md --bibliography=references.bib -o main.pdf
```

Note: `references.bib` currently omits two sources cited only in prose — Roy &
Vetterli (2007, effective rank) and Seurat/Tirosh (2016, `AddModuleScore`);
add them before a citation-complete render.
