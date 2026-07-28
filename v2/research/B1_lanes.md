# MORPHEUS Fix-Research: Lane Scoping (B1)

## Master Question
How do we implement the biology-collapse fixes in THIS codebase, what else is broken, and how do we catch breakage cheaply?

Context: the L2-normalized 256-D biology head has effective rank ~5-6 of 256 (identity ~84). Programme NLL, neighbour-KL, supcon, separation, and a per-dim variance floor all coexist but none prevents rank collapse. Substantial wiring (z_context, z_uncertainty, 5 modality residuals, coordinate path, detached programme memory bank) is computed but never trained or exported. Tests catch none of the observed failure modes.

## Lanes

**1_collapse_remedy — Anti-collapse remedy for the L2-normalized 256-D biology head**
Remit: concrete fixes to raise biology effective rank (whitening/decorrelation, rank-aware loss, normalization change) implemented against model.py/losses.py.
NOT: do not touch the programme-target construction or neighbour-KL formulation (Lane 2).

**2_target_redesign — Programme-target / neighbour-KL redesign to avoid low-rank manifold**
Remit: whether the 50-D Hallmark cancer-residualized target and top-8 neighbour-KL structurally force a low-rank manifold, and how to redesign the target/graph.
NOT: do not propose head-side normalization or decorrelation remedies (Lane 1).

**3_siglip_swap — SigLIP sigmoid loss swap for the identity contrastive**
Remit: replace the symmetric InfoNCE (temp 0.07) on z_identity with a SigLIP sigmoid loss; integration and hyperparameters.
NOT: do not add any contrastive to the biology head or alter biology losses.

**4_grad_conflict — Gradient-conflict instrumentation + mitigation**
Remit: measure and mitigate inter-loss gradient conflict across the active loss set (NLL, KL, supcon, separation, identity).
NOT: do not redefine any individual loss's mathematical form (Lanes 1-3); only weighting/projection/logging.

**5_honest_metric — Confound-aware / honest evaluation metric harness**
Remit: evaluation metrics that expose the cross-cancer molecular-prompting confound (~50% cross-cancer, +0.07 specific over random) under the held-out protocol.
NOT: do not build failing pre-run guards or CI stress tests (Lane 6).

**6_stress_tests — Fast stress tests that fail before an expensive full run**
Remit: cheap CPU tests catching rank collapse, KL/NLL NaN, log_variance mode collapse, no-op loss, anchor-gate saturation; reuse existing fixtures.
NOT: do not define the scientific evaluation metrics themselves (Lane 5).

**7_dead_wiring — Dead-wiring audit and remediation priority**
Remit: audit z_context, z_uncertainty, 20 residual slots, coordinate path, semantic head, detached memory bank; prioritize wire-up vs delete.
NOT: do not implement biology-head or loss redesigns (Lanes 1-4); only classify and sequence.
