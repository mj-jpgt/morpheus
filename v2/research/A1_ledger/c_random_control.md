# Lane c_random_control — Random / null-model controls in gene-signature prediction

Remit: methodology on random-gene-set / size-matched null controls for scoring gene
signatures / pathway activity; is "specificity over random gene sets" standard; is a
model predicting random signatures ~as well as real ones a recognized failure mode /
needed control.

## Queries run
- random gene set null control specificity scoring gene signature pathway activity
- size-matched random gene sets control AUCell VISION signature scoring null model
- deep learning predicts random gene signatures as well as real signatures failure mode benchmark expression
- histopathology image gene expression prediction confounded cross-cancer cohort structure benchmark spurious
- GSEA permutation null model random gene set size bias enrichment score standard control
- most random gene expression signatures associated breast cancer survival prognosis Venet Delorenzi
- pan-cancer batch effect tissue of origin confounds molecular prediction TCGA cohort structure inflates correlation
- benchmark single cell signature scoring specificity random control gene set 2023 2024 false positive rate
- foundation model pathology gene signature prediction random baseline control specificity effective rank collapse 2024 2025
- dimensional collapse contrastive learning VICReg variance regularization embedding rank collapse per-dimension insufficient
- gene signature specificity control reporting standard practice recommendation permutation null bioinformatics best practice

## Sources
- Venet, Dhanasekaran & Sotiriou 2011, "Most Random Gene Expression Signatures Are
  Significantly Associated with Breast Cancer Outcome," PLoS Comput Biol
  (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002240) —
  foundational: random gene sets predict outcome ~as well as real signatures; established
  random-signature null as a required control.
- Domany group 2018, "Association between expression of random gene sets and survival is
  evident in multiple cancer types and may be explained by sub-classification," PMC5839591
  (https://pmc.ncbi.nlm.nih.gov/articles/PMC5839591/) — extends random-bias to 24/34 TCGA
  cancers; attributes it to cohort sub-structure (exactly our cross-cancer confound).
- Zhang et al. 2022, "Signature-scoring methods developed for bulk samples are not
  adequate for cancer scRNA-seq," eLife 71994 (https://elifesciences.org/articles/71994) —
  quantifies false-positive rates on null (downsampled-identity) signatures per method.
- STAR Protocols 2022, "Protocol to benchmark gene expression signature scoring techniques
  for scRNA-seq in cancer" (https://www.cell.com/star-protocols/fulltext/S2666-1667(22)00757-2)
  — formal protocol using random gene sets to compute specificity/false-positive rate;
  shows random-control specificity IS a codified standard practice.
- NARGAB 2024, "Comparative analysis of single-cell pathway scoring methods and a novel
  approach" (https://academic.oup.com/nargab/article/6/3/lqae124/7770961) — random gene
  sets of sizes 50-300 at 5 noise levels as null; AUCell/JASMINE still fire on random sets.
- AUCell / AddModuleScore docs (https://www.bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html)
  — AddModuleScore subtracts a randomly-selected size/expression-matched control feature
  set; size-matched random background is built into standard tooling.
- Ahlmann-Eltze, Huber & Anders 2025, "Simple controls exceed best deep learning
  algorithms ... foundation model effectiveness," PMC12202205
  (https://pmc.ncbi.nlm.nih.gov/articles/PMC12202205/) — a CRISPR-informed mean baseline
  beats GEARS/scGPT; random-init ~= pretrained; field "lacks simple but appropriate
  benchmarks." Direct analogue to our "every method ~+0.07 incl. baseline."
- Domain-adversarial / tissue-of-origin confound work
  (https://arxiv.org/pdf/2504.10343 ; https://www.nature.com/articles/s41698-022-00302-7) —
  tissue-of-origin dominates pan-cancer expression embeddings; standard fix is orthogonal
  projection out of the tissue subspace (= our within-cancer stratification).

## Findings
- Random/size-matched gene-set nulls are an ESTABLISHED, codified control for signature
  scoring. AddModuleScore (Seurat) subtracts a randomly selected, size- and
  expression-bin-matched control set by construction; AUCell/UCell/VISION benchmarks and a
  dedicated STAR Protocol (2022) generate random gene sets (n=50-300, multiple noise levels,
  ~200-5000 sets) purely to measure specificity / false-positive rate. [STAR Protocols 2022;
  NARGAB 2024; AUCell docs]
- "A model scores random signatures about as well as real ones" is a RECOGNIZED failure
  mode, not a novel observation. Venet 2011 showed any random ≥100-gene set has ~90% chance
  of significant breast-cancer survival association; signatures for unrelated phenomena
  (postprandial laughter, mouse social defeat, skin fibroblast localization) were all
  significantly prognostic. [Venet 2011]
- The failure mode is driven by COHORT STRUCTURE, matching our cross-cancer confound
  directly. The 2018 multi-cancer paper found random-gene-set survival bias in 24/34 TCGA
  datasets (up to ~99% of random sets significant in GBMLGG), and showed it is explained by
  latent sub-classification: clustering removed the bias in 65/106 clusters; proliferation
  (PCNA) adjustment alone was insufficient (only 2/17 positively-biased sets lost
  significance). This is the transcriptomic analogue of our finding that ~46-49% of Pearson
  is cross-cancer cohort structure and that within-cancer, random-adjusted specificity is
  tiny. [Domany 2018]
- Nominal significance overstates real signal; empirical/permutation nulls are the corrective
  and are standard. GSEA's own model normalizes ES against gene-permutation nulls to remove
  gene-set-SIZE bias; empirical-p reanalysis of the 48k random signatures left only ~2%
  significant vs real signatures at p<1e-15. So requiring a random-control-ADJUSTED metric is
  established methodology. [GSEA docs; BMC Med Genomics PMC6842262]
- Signature-scoring methods themselves emit substantial false positives on null signatures,
  and specificity over random/null sets is explicitly reported: eLife 2022 gives AUCell 7.5%,
  JASMINE 4%, ssGSEA up to 46% false-positive rates on identity/null comparisons; NARGAB 2024
  reports AUCell/JASMINE still fire on random gene sets while AddModuleScore/scPS are
  cleaner. [eLife 2022; NARGAB 2024]
- The "simple baseline matches every fancy model" pattern is independently recognized in
  genomics DL: Ahlmann-Eltze et al. 2025 show a mean-expression baseline beats GEARS/scGPT
  (Pearson-delta +0.08/+0.11), pretrained ~= random init (0.004, p=0.89), and explicitly
  criticize the field for lacking appropriate baseline controls. This mirrors our result that
  the random-adjusted within-cancer gain is ~+0.07 for EVERY method including the MLP-CLIP
  baseline. [Ahlmann-Eltze 2025]
- Tissue-of-origin / cohort confounding of molecular prediction is well documented; the
  recommended remedy (orthogonal projection out of the tissue subspace, or within-stratum
  evaluation) is exactly our within-cancer decomposition, though it is applied to
  expression-embedding and survival tasks, not WSI->Hallmark prompting. [arXiv 2504.10343;
  npj Precision Oncology 2022]

## Novelty verdict
PARTIALLY NOVEL — the CONTROL and FAILURE MODE are prior art; the SPECIFIC APPLICATION to
WSI->molecular-programme prompting with a quantified decomposition is not something we found
reported.

Already established (our finding is NOT novel on these points):
- Random / size-matched gene-set nulls as a specificity control — standard, codified
  (STAR Protocols 2022, AddModuleScore, GSEA gene-permutation null).
- "Predicts random signatures ~as well as real ones" as a diagnosed failure mode — Venet
  2011 is the canonical citation; reviewers in this area will expect it.
- The mechanism (cohort/tissue sub-structure inflating apparent signal, proliferation not
  fully explaining it, within-stratum evaluation as the fix) — Domany 2018 and the
  tissue-of-origin confounding literature.
- Simple/random baselines matching sophisticated models, and the field under-using baseline
  controls — Ahlmann-Eltze 2025 (perturbation prediction).

Appears novel (nothing found reporting this exact result):
- Bringing the random-gene-set null control specifically into the WSI (histology image) ->
  Hallmark molecular-programme PROMPTING benchmark. Prior random-null work is on
  transcriptome-based survival/scoring, not image->pathway prompting.
- The quantified decomposition on that task: ~46-49% of global Pearson = cross-cancer cohort
  structure (0.348 -> 0.188 within-cancer), random-gene-set control ~0.30-0.32 global, and a
  genuine within-cancer, random-control-adjusted specificity of only ~+0.07 for every method
  including the baseline. The magnitudes and the "random control scores nearly as high as the
  real target" quantification in this modality are our own contribution.
- SigLIP marginally beating MLP-CLIP (+0.005 within-cancer, 62% of targets) under this
  corrected metric — an application-specific empirical result with no external analogue found.

Bottom line: our finding is a correct and well-grounded APPLICATION of an established control
philosophy to a benchmark where it had not been applied. The reusable contribution is the
quantified "random-control-adjusted within-cancer specificity" recipe for WSI->pathway
prompting, not the general idea that random controls are needed (which is standard and should
be cited to Venet 2011 / Domany 2018 / STAR Protocols 2022 to preempt reviewer pushback).
