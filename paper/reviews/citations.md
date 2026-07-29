# Citation Verification — MORPHEUS rank-collapse paper

**Scope:** every entry in `paper/references.bib` (17 BibTeX entries + 2 prose-only citations), checked against `main.md` for (a) existence with the stated title/venue/year and (b) support for the claim it is attached to. Verified via WebSearch/WebFetch, July 2026.

**Headline:** No hallucinated references — every cited work exists. But **two entries are mis-attributed**: the claim they support belongs to a *different, real* paper. Both must be fixed before submission.

---

## MAJOR — mis-attributions (2)

### M1. `wang2020understanding` — wrong paper for the claim (and garbled in-text label)

- **Bib entry holds:** Wang & Isola, *"Understanding Contrastive Representation Learning through Alignment and Uniformity on the Hypersphere,"* ICML 2020. **This paper is real** (ICML 2020, PMLR v119, pp. 9929–9939) — but it is about the alignment and uniformity properties of the contrastive loss.
- **Claim it is attached to** (§2.4, also §5.4, §6.4, abstract/intro concession list): *"Wang et al. showed the minimal-sufficient representation from contrastive learning discards task-relevant, non-shared information and is therefore insufficient downstream."* Labeled in prose as **"Wang & Isola, CVPR 2022" / "Wang, CVPR 2022" / "Wang CVPR22."**
- **Problem:** The alignment/uniformity paper makes no such minimal-sufficient-representation claim. That claim is the thesis of a **different** paper: **Haoqing Wang, Xun Guo, Zhi-Hong Deng, Yan Lu, "Rethinking Minimal Sufficient Representation in Contrastive Learning," CVPR 2022 (Oral), arXiv:2203.07004.** The prose label "Wang & Isola, CVPR 2022" is a garbled hybrid — it pairs the *authors* of the ICML 2020 paper (Tongzhou Wang & Phillip Isola) with the *venue/year* of the CVPR 2022 paper (whose authors are Haoqing Wang et al., not Isola).
- **Correction:** Replace the bib entry with the CVPR 2022 paper:
  ```bibtex
  @inproceedings{wang2022rethinking,
    title     = {Rethinking Minimal Sufficient Representation in Contrastive Learning},
    author    = {Wang, Haoqing and Guo, Xun and Deng, Zhi-Hong and Lu, Yan},
    booktitle = {Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
    year      = {2022},
    eprint    = {2203.07004},
    archivePrefix = {arXiv},
    primaryClass  = {cs.CV}
  }
  ```
  Update the in-text label from "Wang & Isola, CVPR 2022" to "Wang et al., CVPR 2022" everywhere (§2.4, §5.4, §6.4, abstract/intro, T-tables). If the alignment/uniformity result is *also* wanted (it is not currently used for any claim), keep the ICML 2020 entry under a distinct key — but as written, no claim in the paper needs it.

### M2. `tizhoosh2026rethinking` — "Buyer Beware" claim belongs to a different paper

- **Bib entry holds:** Tizhoosh, *"Rethinking foundation models in pathology,"* Nature Biomedical Engineering 2026, DOI 10.1038/s41551-026-01696-6. **This paper is real** — but it is a comment about *foundation-model architecture* (dense embeddings can't represent tissue's combinatorial richness; flaws in self-supervision, patch design, noise-fragile pretraining). It does **not** discuss confounders inflating omics-from-histology biomarkers.
- **Claim it is attached to** (§2.2, also abstract/intro, §5 concession, R2 in the risks table): *"'Buyer Beware' catalogued how co-dependencies and confounders (e.g., tumour mutational burden) inflate omics-from-histology biomarker estimates."* Cited throughout as **"Buyer Beware, NBME 2026."**
- **Problem:** That claim precisely describes a **different** paper actually titled *"Buyer Beware"*: **Muhammad Dawood, Kim Branson, Sabine Tejpar, Nasir Rajpoot, Fayyaz ul Amir Afsar Minhas, "Buyer Beware: confounding factors and biases abound when predicting omics-based biomarkers from histological images," bioRxiv 2024, doi:10.1101/2024.06.23.600257.** That paper documents that mutation predictions are confounded by overall tumour mutational burden and that biomarker co-dependencies drive apparent performance — exactly the cited claim. The bib note (`Comment / "buyer beware" on pathology foundation models`) shows the two papers were conflated: the authors labeled the Tizhoosh comment "Buyer Beware," but the real "Buyer Beware" paper is Dawood et al.
- **Correction:** Add the Dawood et al. entry and attach the confounder/TMB claim to it:
  ```bibtex
  @article{dawood2024buyer,
    title   = {Buyer Beware: confounding factors and biases abound when predicting omics-based biomarkers from histological images},
    author  = {Dawood, Muhammad and Branson, Kim and Tejpar, Sabine and Rajpoot, Nasir and Minhas, Fayyaz ul Amir Afsar},
    journal = {bioRxiv},
    year    = {2024},
    doi     = {10.1101/2024.06.23.600257}
  }
  ```
  Retarget the "Buyer Beware" in-text citations (abstract/intro, §2.2, §5, R2) to `dawood2024buyer`. The Tizhoosh NBME 2026 entry may be kept **only** if used for a foundation-model-architecture point — but as currently written every use of it is the confounder claim, which it does not support. Also drop or fix the misleading "NBME 2026" venue label wherever "Buyer Beware" appears.

---

## VERIFIED CORRECT (15)

All of the following exist with the stated title, venue, and year, and support the claim they are attached to in `main.md`:

| key | verdict |
|---|---|
| `jaume2024hest` | OK. HEST-1k, NeurIPS 2024 Datasets & Benchmarks Track, arXiv:2406.16192. Supports "1,229 ST/WSI pairs + HEST-Benchmark" (§2.1). |
| `schmauch2020he2rna` | OK. Nat. Commun. 11:3877 (2020), DOI 10.1038/s41467-020-17678-4. Supports HE2RNA / random-gene-null lineage (§2.1, §4.4, §4.7, §5). |
| `fu2020pancancer` | OK. Nat. Cancer 1(8):800–810 (2020), DOI 10.1038/s43018-020-0085-8. Supports "signal reflects tumour type/composition; ~50% halving" (§2.1, §4.3, §5.2). |
| `howard2021impact` | OK. Nat. Commun. 12:4423 (2021), DOI 10.1038/s41467-021-24698-1. Supports "site signatures survive normalization, bias predictions" (§2.2). |
| `dehkharghanian2023biased` | OK. Diagnostic Pathology 18:67 (2023), DOI 10.1186/s13000-023-01355-3. Supports "networks recover TCGA acquisition site" (§2.2). |
| `steiner2026decat` | OK. DECAT, arXiv:2605.31504 (2026). Supports "null-referenced, confounder-label-free decision procedure" (§2.2, R1). |
| `gindra2025hescape` | OK. HESCAPE, arXiv:2508.01490 (2025). Supports "contrastive pretraining degrades expression prediction; batch effects" (§2.2). |
| `jing2022understanding` | OK. ICLR 2022 (arXiv:2110.09348). Supports "dimensional collapse" (§2.3, §3.2). Minor: verify the `url` OpenReview forum id and consider adding `eprint=2110.09348`. |
| `andriopoulos2024prevalence` | OK. NeurIPS 2024 (arXiv:2409.04180). NRC1 = last-layer features collapse to target-subspace span; matches claim (§2.3, §3.2). |
| `bardes2022vicreg` | OK. VICReg, ICLR 2022 (arXiv:2105.04906). Variance hinge + off-diagonal covariance term; matches §2.3/§3.3 fix framing. |
| `liang2023factorized` | OK. FactorCL, NeurIPS 2023 (arXiv:2306.05268). Shared/unique factorization; matches §2.4. |
| `venet2011most` | OK. PLoS Comput. Biol. 7(10):e1002240 (2011), DOI 10.1371/journal.pcbi.1002240. Random signatures predict outcome; matches §4.4/§4.7. |
| `zhai2023sigmoid` | OK. SigLIP, ICCV 2023, pp. 11975–11986 (arXiv:2303.15343). Used as the SigLIP baseline (§5). |

---

## MINOR / prose-only citations to add

- **Roy & Vetterli (2007)** — cited in §4.6 for the effective-rank definition; flagged in `main.md` as pending. Real: O. Roy & M. Vetterli, *"The effective rank: A measure of effective dimensionality,"* EUSIPCO 2007, pp. 606–610. Add to bib.
- **Tirosh 2016 / Seurat `AddModuleScore`** — cited in §4.4/§4.7 for the background-subtracted module-score control; flagged as pending. Real: Tirosh et al., *"Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq,"* Science 352(6282):189–196 (2016), DOI 10.1126/science.aad0501 (method operationalized in Seurat's `AddModuleScore`). Add to bib.
- **`tizhoosh2026rethinking` note field** — remove the `"buyer beware"` label from the note regardless of M2 resolution; it drove the conflation.

---

## Bottom line

No fabricated references — all 17 bib entries and both prose citations point to real works. But two are mis-attributed to claims made by *other* real papers (M1: Haoqing Wang CVPR 2022, not Wang & Isola ICML 2020; M2: Dawood et al. "Buyer Beware" bioRxiv 2024, not Tizhoosh NBME 2026). Both errors propagate across the abstract, related work, results-concession list, and risks table, so the fix touches multiple sections. These are correctness-of-attribution defects, not citation-formatting nits.
