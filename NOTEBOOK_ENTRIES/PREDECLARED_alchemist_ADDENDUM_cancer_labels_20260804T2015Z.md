# ADDENDUM 2 to the ALCHEMIST predeclaration — the cancer-label decision, declared before measuring

**UTC** 2026-08-04T20:15Z
**Status** written **before** any ALCHEMIST channel number exists. Extraction stood at 1,029 of
1,106 slides; `RESULT*` directories were empty; the chained analysis job was stopped by PID
before it fired precisely so this could be declared first.

**Amends** `PREDECLARED_alchemist_external_replication_20260804T1830Z.md` §1 (Adjustment).
**Prompted by** a coordinator report that ALCHEMIST's cancer labels are disjoint from TCGA's, so
every row refuses with `UnseenLevelError` under the persisted TCGA operator.

---

## 1. The reported blocker does not bind on this measurement, and the reason matters

The coordinator's diagnosis of `on_unseen_level="zero"` is correct and I am not adopting it: a
zero one-hot row is the design's *reference* level, so zeroing would silently file every
ALCHEMIST patient under whichever TCGA cancer happens to be the dropped category.

But this comparison never applies a TCGA-fitted operator to ALCHEMIST rows. **Each cohort's
confound design is fitted within the cohort being measured**, and each cohort is residualised
against its own design and correlated with its own targets. The cohorts are never pooled, so no
cross-cohort operator application occurs and `UnseenLevelError` cannot arise. This is not an
oversight being defended after the fact; it is what "measure the channel on ALCHEMIST with the
same instrument" has to mean for a *replication* question rather than a *transfer* question.

Applying the persisted TCGA operator to ALCHEMIST would actively introduce the error it is meant
to prevent. A confound design exists to remove the cancer-type variance *present in the rows
being adjusted*. A TCGA-fitted operator encodes TCGA's cancer means; subtracting those from
ALCHEMIST patients does not remove ALCHEMIST's histology structure, it adds a constant offset
and leaves the histology confound sitting **inside** the ALCHEMIST channel, inflating it.

## 2. The decisive fact, verified from the operator's own provenance

`runs_misc/tcga_operators/tcga_operator__wsi_identity__cancer_only__test_seed42.provenance.json`
records `n_reference_rows: 2530` and 22 `design_columns`:

```
cancer_ACC BLCA CESC CHOL COAD DLBC ESCA KICH KIRP LIHC MESO PAAD PCPG READ
       SARC SKCM STAD TGCT THYM UCS UVM  +  cancer_nan
```

**There is no `cancer_LUAD` and no `cancer_LUSC`.** Cross-checked against the frozen cohort:
the TCGA test split is 2,766 patients over exactly those 21 cancers, holding **0 LUAD and
0 LUSC** — all 846 TCGA NSCLC patients sit in train (748) and val (98), by the 11v21 holdout
design.

So the persisted `cancer_only` operator **cannot adjust either arm** of a lung-versus-lung
comparison. TCGA's own LUAD and LUSC patients raise `UnseenLevelError` against it exactly as
ALCHEMIST does. The blocker is therefore not an ALCHEMIST label-mapping problem; it is that the
operator was fitted on a patient set from which lung cancer is entirely absent.

This also makes the coordinator's option 1 as literally worded — "map to the TCGA codes the
operator was fitted on" — impossible: the operator was not fitted on LUAD or LUSC. And option 2,
"map the whole cohort to a single TCGA lung code", would map onto a column that does not exist.

## 3. What the ALCHEMIST arm therefore uses, and the answer to the concrete question

**Does each ALCHEMIST case carry a histology distinguishing adenocarcinoma from squamous?
Yes — verified, not assumed.** `api.gdc.cancer.gov/cases` faceted on
`diagnoses.primary_diagnosis` is populated for **all 1,176** cases (unlike
`tissue_source_site`, which is `_missing` for all 1,176, and `ajcc_pathologic_stage`, likewise).
In the 1,106-patient paired cohort: Adenocarcinoma NOS 808, Squamous cell carcinoma NOS 209,
adenosquamous 21, adenocarcinoma in situ 14, large cell 13, carcinoma NOS 11, and 8 rarer types.

**Primary arm:** ALCHEMIST's own `primary_diagnosis`, rare levels pooled to `OTHER` at the same
`min_site_count = 10` rule the project uses elsewhere — 7 levels. This is the design that
actually removes ALCHEMIST's histology variance, which is the job of the adjustment.

**Declared sensitivity arm, and this is the coordinator's option 1 run as a check rather than as
the primary:** the same histology collapsed onto the TCGA lung code space via an explicit,
reviewable mapping (`ALCHEMIST_TO_TCGA_LUNG` in `v2/research/external/alchemist_channel.py`):

| ALCHEMIST histology | mapped |
|---|---|
| adenocarcinoma NOS / in situ / acinar / mucinous / invasive mucinous / mixed mucinous / lepidic / mixed subtypes | **LUAD** |
| squamous cell carcinoma NOS | **LUSC** |
| adenosquamous, large cell, carcinoma NOS, non-small cell, pleomorphic, spindle cell, large cell neuroendocrine | **OTHER** |

Adenosquamous carcinoma is genuinely both lineages, so it goes to `OTHER` rather than being
forced to whichever of LUAD/LUSC would flatter the result.

**What it assumes:** that ICD-O morphology families map onto TCGA's LUAD/LUSC project codes.
That is a histology-to-project-code identification, not a claim that an ALCHEMIST adenocarcinoma
and a TCGA LUAD case are exchangeable. It is used only to test label granularity.

**The predeclared reading, fixed now:** if the ALCHEMIST channel is materially different under
the 7-level and 3-level designs, the cancer-label choice is setting the answer and I will report
it as a finding rather than pick the better number. If the two agree, the choice is not
load-bearing and both are reported anyway.

## 4. Two points from the coordinator adopted verbatim

- **The comparison is cancer-adjusted only, on both arms.** Already the case since the original
  predeclaration §1, for the reason given there: ALCHEMIST has no site field, and an adjustment
  applied to one arm only is not the same instrument. TCGA is additionally reported with
  `cancer + pooled TSS` so the cost of dropping the site term is visible.
- **The dilution curve cannot price this contaminant**, per
  `PREDECLARED_alchemist_ADDENDUM_dilution_arm_20260804T1905Z.md`. The external channel arrives
  without a predicted penalty. The expected direction is stated qualitatively only, and a weak
  result will not be converted into "tissue-level sampling explains it" — there is no number
  behind that sentence and this project has not measured one.

## 5. Unchanged

The bar in §4 of the original predeclaration — REPLICATES at p < 0.01 and R >= 0.60, ATTENUATED
at 0.30 <= R < 0.60, FAILS otherwise — is **not** modified. R is computed from the primary arm.
The cohort-classifier AUC is printed before any channel number, by construction.
