# PREDECLARED — P5 §4.3 novelty-filter feasibility scoping, and a stages-0-4 pilot funnel

**UTC** 2026-08-05T07:50Z
**Task:** `paper/P5_DISCOVERY_PLAN.md` §4.3 (novelty-filter feasibility scoping) and §5 (a small
pilot of the discovery funnel's stages 0-4, novelty tier-1 only). Written by the agent sent per
the plan's closing line ("An agent is being sent on the §4.3 novelty-filter feasibility scoping
and a first draft of the stage 0-4 pilot pipeline now").

---

## 0. A process deviation, disclosed before anything else

**Part A's core connectivity checks were already run before this file was committed.** Before
realising the full predeclare-then-measure discipline applied to environment-scoping as well as
to scientific measurement, four requests were made to live external APIs as basic environment
reconnaissance: (1) the exact PubMed `esearch` query given in the task brief
(`TP53[Title/Abstract] AND lung adenocarcinoma[Title/Abstract] AND (histology OR morphology)`),
(2) one `esummary` call, (3) a 10-request rapid-fire burst against `esearch` to see whether a 429
appears, (4) one Open Targets GraphQL query for TP53's associated diseases. All four are
reachability/behaviour facts (HTTP status code, JSON shape, timing), not results with a
favourable/unfavourable direction a predeclaration is designed to keep honest — there is no way
to unconsciously shade "did the socket connect" — but the rule is "predeclare before you measure"
without a carve-out for that judgment, so this is recorded as a deviation rather than quietly
absorbed. Everything from here on, including the remaining Part A characterisation (sustained
rate-limit behaviour, Open Targets burst behaviour, a deliberate zero-hit query) and all of Part
B, is predeclared before it is run.

---

## 1. Part A — novelty-filter feasibility scoping: what would count as what

**Already known before running anything further** (from the four calls above, reported in full
regardless of which way they cut): both `eutils.ncbi.nlm.nih.gov` and
`api.platform.opentargets.org` are reachable over HTTPS from this sandbox and returned real,
parseable data for the exact query given in the brief. That is a favourable-looking start and is
exactly the kind of result rule 3 (push until it breaks) says not to bank — the rest of this
section is what would break it.

### What is being measured next
1. Sustained-rate behaviour of PubMed E-utilities without an API key: repeated `esearch` calls at
   decreasing intervals to find where `429` starts, and whether a ~3 req/s pace avoids it.
2. Whether Open Targets shows any rate limiting under a rapid burst (it has no documented
   anonymous-key requirement, unlike PubMed).
3. A deliberately obscure/unlikely gene-cancer-morphology combination, to see the *shape* of a
   genuine "no hit" response (empty `idlist`, not an error) — needed because a "no hit" that is
   actually a malformed query would be indistinguishable from a real negative without checking.

### What would make the "usable for a novelty check" conclusion favourable but wrong
- If the exact query from the brief returned 0 hits (it does not — 413, confirmed above) that
  would have suggested the query syntax itself was broken, not that PubMed is reachable.
- If rate-limiting is *silent* (e.g. empty-but-200 responses under load) rather than an explicit
  `429`, an automated tier-3 checker would misread "was throttled" as "found nothing", which is
  the single most dangerous failure mode for a "no obvious prior hit" claim. This is checked for
  explicitly, not assumed absent because a 429 was seen elsewhere.
- Any result here says only that the **query mechanics** work. It says nothing about whether the
  free-text/MeSH mapping actually surfaces the kind of paper this pipeline would need to find
  (a gene-morphology correlational finding, not a gene-disease association) — see §1.5 below,
  predeclared as the honest limitation regardless of how the connectivity checks land.

### False-negative judgment call, stated as a judgment, not computed
A "no hit" from either API will be reported as differently trustworthy depending on tier:
- **Open Targets tier-2** (curated gene-disease association scores aggregating genetic,
  somatic-mutation, drug and literature evidence): a "no hit" is reasonably strong evidence of
  absence *for the gene-disease relationship itself*, because the resource is curated and
  score-based rather than a raw text index. It says nothing about a gene-**morphology** claim.
- **PubMed tier-3** (title/abstract text search): a "no hit" is weak evidence. It misses gene
  alias mismatches (a search on the HGNC symbol will not find a paper using only "p53"-style
  historical naming unless MeSH translation catches it), non-English literature, preprints not
  yet indexed, and — most importantly for this specific pipeline — the tier-3 query as specified
  in the plan (`gene AND cancer-type AND (histology OR morphology OR imaging)`) is a **generic**
  co-occurrence query, not a search for "this gene's expression correlates with this specific
  morphological/imaging axis", which is a narrow publication niche that plain title/abstract
  MeSH search is not well shaped to find even when such a paper exists. Expect this tier's
  false-negative rate for the pipeline's *actual* claim type to be meaningfully higher than for a
  generic "is this gene studied in this cancer at all" query — a double-digit-percent miss rate
  is a defensible expectation, not a number this pass will try to compute precisely.

---

## 2. Part B — pilot funnel, stages 0-4, novelty tier-1 only

### 2.1 What is expected before running: the data-access question first

**Prediction, stated before checking:** given this is a Windows checkout mirrored from a Lambda
box per the task brief, and the last three notebook entries that reference a live box give three
different IPs across four days (`132.145.196.200`, `150.136.45.194`, and the two addresses in
`~/.ssh/config`, `150.136.95.220` / `150.136.215.164`), the actual current box is expected to be
unreachable or to require credentials this checkout does not hold — Lambda instances are commonly
re-provisioned with a new address between sessions and nothing in the repo pins "the current one".
**If that prediction is wrong** (a box is reachable and holds a current, hash-verifiable
`frozen_rna_targets.npz` and at least one representation-state artifact) the real pilot in §2.2
runs on real data. **If it is right**, §2.3's synthetic-mechanics fallback runs instead, and is
reported as exactly that — a code-path validation, not a biological result — regardless of how
favourable its numbers look. No number produced by §2.3 will be described anywhere as a TCGA or
biological finding.

### 2.2 If real data is reachable: the candidate space, fixed now

- **Candidate targets:** every column of `frozen_rna_targets.npz` whose `target_groups` value is
  **not** `hallmark_in_training` (that group is the D2 arm-H supervision target itself — scoring
  it is measuring training-signal recovery, not discovery, per the plan's tier-1 novelty filter).
  From the evidence-ledger counts already on record (`24 heldout_pathway`, `8 immune_tme`,
  `8 tumour_state`), that is **40 candidate target columns**, fixed by the artifact's own schema,
  not chosen after looking at results.
- **Cancer strata:** the 5 largest cancers by patient count in the loaded artifact's paired split
  (ties broken alphabetically by TCGA code), each cancer's cohort restricted to that cancer's own
  patients only (a within-stratum analysis, not the pooled cancer-type confound `run_calibra.py`
  normally adjusts for — restricting to one cancer removes the need to adjust for it and lets a
  finding be about *that* cancer specifically, matching the plan's "cross the target space with
  cancer strata").
- **Candidate space size: 40 x 5 = 200 cells**, fixed before stage 1 runs.
- **Stage 1 (coarse pre-filter):** `spectral.top_canonical_correlation(x, y_column, n_components=1)`,
  in-sample, unadjusted — cheap and canonical. **Predeclared keep rule:** the top 40% of the 200
  cells by `|r|` advance to stage 2 (≈80 cells), fixed as a fraction rather than a magnitude
  threshold so the stage-2 test count is knowable before any real number is seen.
- **Stage 2 (certify):** `residualise.confound_design` (pooled tissue-source-site only, since
  cancer type is constant within a stratum), `residualise.cross_fitted_residuals`,
  `calibration.permutation_null`, `calibration.spike_recovery_curve`,
  `spectral.heldout_single_direction_correlation` as the channel statistic (the single-target
  analogue of `run_calibra._channel_measurement`'s `heldout_top_cca`), graded with
  `calibration.channel_clears_floor` — unchanged imports, no reimplementation.
- **Stage 3 (FDR):** `scipy.stats.false_discovery_control(method="bh")` — a canonical library
  implementation, not a hand-rolled Benjamini-Hochberg formula, applied **once**, across all ≈80
  stage-2 permutation p-values (the whole predeclared stage-2 test set, not just the ones that
  look good), at **q = 0.10**, fixed now.
- **Stage 4 (replicate):** a further split of each stratum's patients into two halves at stage 0
  (discovery / replication), with the certified statistic re-measured on the replication half for
  every stage-3 survivor. ALCHEMIST is not attempted in this pilot even if TCGA data is reachable
  — it needs its own patch extraction and cohort assembly (see
  `NOTEBOOK_ENTRIES/PREDECLARED_alchemist_external_replication_20260804T1830Z.md`), which is out
  of scope for a same-day pilot.

### 2.3 If real data is not reachable: the synthetic mechanics fallback, fixed now

A controlled synthetic ladder, explicitly **not** a biological measurement, whose only purpose is
to check the funnel's stages behave the way the plan says they should — an engineering
counterfactual (rule 6: every claim ships a must-fail and a must-pass control), not a discovery.

- **200 synthetic cells**, matching the real design's 40 x 5 shape: 5 synthetic "strata"
  (labelled with real TCGA cancer codes for readability — BRCA, LUAD, KIRC, HNSC, THCA — with **no
  claim that the underlying values are those cancers' real biology**) x 40 synthetic target
  columns, ~250 synthetic patients per stratum (chosen to be of the right order for a large TCGA
  cancer type), 32 synthetic representation axes (reduced from the plan's 256 for CPU wall-clock
  in a same-day pilot — stated as a deviation, not hidden).
- **20 of the 200 cells (10%) carry a planted signal** at `r_true = 0.18`, built with
  `calibration.spike_targets` against pure-noise axes and a pure-noise target column, so the axes
  and target for those 20 cells are otherwise indistinguishable from the 180 null cells except for
  the planted direction. The remaining 180 are pure Gaussian noise plus the same confound
  structure (a synthetic site covariate) as the planted cells, so the null cells are a real null,
  not an easy one.
- **Same stage 1/2/3/4 procedure as §2.2**, called on the synthetic arrays through the identical
  code path (one pilot module, one `--synthetic-dry-run` flag switching only the data source).
- **Predeclared pass/fail for the funnel itself:**
  - *Must-pass:* a majority of the 20 planted cells should survive stage 1 (their true |r| is
    constructed to sit above the median of the null distribution) and a nontrivial fraction should
    still be BH-FDR significant at q=0.10 after stage 2/3 — a funnel that loses every planted
    signal by stage 3 has a shape problem worth reporting, not hiding.
  - *Must-fail:* the **BH-FDR survivor count among the 180 pure-null cells should be small** —
    under the null, BH at q=0.10 should keep the expected false-discovery proportion at or below
    10% of whatever it calls significant; a null-cell survival rate wildly above that is a defect
    in this pilot's stage-2/3 wiring, not a permissible finding, and will be reported as a
    methodology-repair item, not smoothed over.
  - **What would make this pilot's "the funnel works" conclusion distrusted even if the counts
    look right:** if the planted-cell recovery rate is high only because the coarse stage-1 filter
    or the permutation null is somehow leaking the planted direction (e.g. a shared random seed
    reused where it should not be) rather than genuinely detecting it — checked by confirming the
    180 null cells use independently drawn noise from the 20 planted cells, not the same draws
    with a spike added on top of a shared base.
- **Explicitly out of scope for this fallback:** no claim about real TCGA attrition rates, no
  claim about how many *real* candidates would survive stage 3 — that number can only come from
  §2.2, and if §2.2 does not run, it is not estimated, guessed, or interpolated from the synthetic
  ladder's counts.

### 2.2/2.3 shared reporting requirement
Whichever path runs, the deliverable is the ledger the plan asks for: candidates entered, stage-1
survivors, stage-2-certified, stage-3 BH-FDR survivors, stage-4 replicated — reported in full even
if the terminal count is zero, per rule 2 (bad news first) and the plan's own §6 ("report the
attrition at every stage rather than only the final shortlist").

---

## 3. What would make either part's headline claim distrusted on its own terms

- Part A: a "reachable" verdict that turns out to depend on a cached DNS/proxy answer rather than
  a live round trip — checked by requiring a real HTTP status code and a real, query-specific JSON
  body (not a cached 304 or an error page with a 200 status) on every call quoted.
- Part B (§2.2 path): any per-cell statistic computed by an inline formula instead of the imported
  `calibra` function — checked by import identity, the same discipline
  `test_effective_rank_canonical.py` already enforces for rank.
- Part B (§2.3 path): reported as if it were evidence about real tumours. It is not, and every
  place it is written up will say so in the same sentence as the numbers, not once at the top and
  then dropped.
