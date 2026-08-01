# A1: Research Lane Scoping

## Master Question
Is the MORPHEUS multi-objective head-collapse + confounded WSI->molecular benchmark finding a real, novel, usable contribution, and what is the strongest paper?

## Grounding Results (seed-42 MORPHEUS V2)
- Biology head collapsed to effective rank ~5-6 of 256; identity head healthy ~84; anchor residual ~= 0 (identity == frozen MLP-CLIP teacher).
- WSI->Hallmark prompting benchmark confounded: ~46-49% of Pearson is cross-cancer cohort structure (global 0.348 -> within-cancer 0.188 for MLP-CLIP); random-gene-set control ~0.30-0.32 global; genuine within-cancer, random-adjusted specificity ~+0.07 Pearson for EVERY method incl. baseline.
- SigLIP is the only variant beating MLP-CLIP on molecular prompting (+0.005 within-cancer, wins 62% of targets).
- Mechanism: 256-D biology head aligned to a ~50-D Hallmark target via neighbour-KL + supcon; per-dimension VICReg variance floor (0.01) did NOT prevent rank collapse.

## Lanes

**a_prior_benchmarks — Prior WSI->molecular/pathway prediction benchmarks.** Find existing WSI->expression/pathway/Hallmark benchmarks (HEST-bench, ST-Net, Hist2ST, BLEEP, mclSTExp, TITAN/CONCH); assess whether they report the cross-cancer confound, within-cohort controls, or random-gene-set nulls. MUST NOT: cover generic batch-effect critiques or random-null methodology (b, c).

**b_batch_confound — Batch/cohort-confound critiques in comp-path & multimodal genomics.** Find work showing performance driven by site/batch/cancer-type shortcuts (Howard, Dehkharghanian, site-preserved splits). MUST NOT: enumerate specific WSI->molecular benchmarks or their metrics (a).

**c_random_control — Random/null-model controls in gene-signature prediction.** Find size-matched random-gene-set / permutation null methodology (GSEA, VISION, AUCell); is "specificity over random signatures" standard? MUST NOT: address cohort/batch confounds or head-collapse geometry (b, d).

**d_head_collapse — Multi-objective/dimensional head collapse as a named phenomenon.** Find literature on rank/dimensional/neural collapse from low-rank-target alignment or multi-task interference (DirectCLR, low-rank simplicity bias). MUST NOT: cover the retrieval-vs-regression trade-off theory (f).

**e_benchmark_packaging — How ML benchmarks get adopted.** Research task cards, standardized splits, leaderboards, eval harness, Croissant/HF hosting, held-out servers; define minimum viable packaging. MUST NOT: assess scientific novelty or closest papers (g).

**f_theory — Why one shared latent can't be optimal for retrieval + regression.** Find identifiability/geometry/rate-distortion results on alignment-uniformity vs regression capacity (FactorCL, modality gap, sufficient-vs-minimal reps). MUST NOT: cover empirical collapse naming (d).

**g_venue_novelty — Venue/novelty scan + 3 closest papers.** Scan MICCAI, CVPR/ICCV-medical, Nature Methods/Comms, NeurIPS D&B, ICLR; name 3 closest papers and articulate our delta. MUST NOT: argue the null case against contribution (h).

**h_steelman_null — Steelman that this is NOT a contribution.** Build the strongest adversarial case (confound already known, collapse a trivial bug, ~0.07 signal uninteresting, site-preserved splits solve it); cite counter-evidence. MUST NOT: advocate for novelty or select the strongest positive paper (g).
