# Decision: `no_external_cohort` stays undischarged, and the reason is granularity

2026-08-04 · main session · decides the item `alchemist_external_replication_RESULT_20260804T2115Z.md`
deliberately left open.

## The decision

ALCHEMIST is acquired, the channel replicates (R = 1.110 at n matched, p = 0.0033), and
`no_external_cohort` **stays up**. This is not caution for its own sake.

## Why not the AUC

The obvious reason to hesitate is the cohort classifier: **AUC 0.99906** separating TCGA from
ALCHEMIST, within-TCGA control 0.50016. That number is real and must be printed beside every
external figure. But it does not by itself block anything, because **each cohort is adjusted and
measured entirely inside itself and the two are never pooled** — a between-cohort direction is not
available to either measurement. What the AUC forbids is a *cohort-invariance* claim, and no such
claim exists on this project. Replication of a finding is not transfer of a model.

## The reason that does bind

`_REQUIREMENTS` gates exactly two claim kinds on this blocker:

    legible_axis     -> composition_attribution, purity_confound, no_external_cohort
    gene_attribution -> composition_attribution, purity_confound, sign_blind, no_external_cohort

Both are **per-axis** claims. `claim_evidence.json` registers one of them,
`morphology_to_pbs_axis_legibility` (kind `legible_axis`).

ALCHEMIST replicated the **aggregate channel across 59 targets**. It did not replicate any
individual axis, and `_is_discharged` cannot tell the difference — it tests
`len(external_cohorts) >= 1`. Writing `["ALCHEMIST"]` into that field would therefore let a
per-axis claim pass on evidence measured at a different granularity. That is precisely the failure
mode the guards exist to prevent, and it would be a self-inflicted one.

**A per-axis external result would discharge it. A pooled-channel external result does not.**

## What this costs, and what it does not

Nothing that was admissible becomes inadmissible, and P1's headline is unaffected: the channel
claim is not of kind `legible_axis` or `gene_attribution`, so it was **never gated by this blocker
at all**. The external cohort strengthens it voluntarily rather than unblocking it. The claim that
stays blocked — morphology reading named PBS axes — was already blocked by
`composition_attribution` and `purity_confound` independently, so this decision changes no verdict
today. It only prevents a future discharge from being taken on the wrong evidence.

## The latent defect this exposed

`_is_discharged` treats external evidence as a count, with no record of *what* was replicated. Any
future external result, at any granularity, discharges every claim kind equally. Logged as a real
weakness in the guard rather than repaired in the same breath as the result that revealed it —
repairing it now would mean the fix and its motivating case land in one undiffable commit.

Related: [[alchemist_external_replication_RESULT_20260804T2115Z]],
[[PREDECLARED_alchemist_ADDENDUM_cancer_labels_20260804T2015Z]]
