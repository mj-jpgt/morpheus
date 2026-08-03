## 2026-08-03 01:55 UTC — T1.7(b) positive control PASSES: image-only ER status recovered at 0.878 against a pre-registered literature band of [0.78, 0.92]

**Logged:** 2026-08-03 01:55 UTC. **How obtained:** `python -m morpheus.v2.calibra.known_covariate_control --covariate-column er_positive --expected-low 0.78 --expected-high 0.92 --preregistration PREREG_known_covariate.json --partition all --n-boot 1000 --n-permutations 1000`, Lambda box `~/ws_p1`, artifacts `d2_h_seed42.npz` and `d2_i_seed42.npz`.

### Technical

**Substitution, recorded not silent.** The plan named MSI, TP53 or consensus subtype. None is usable
from the data on either machine. The TCGA PanCan clinical mirror carries `microsatellite_instability`
only as an assay-performed flag — 7 `YES`, 74 `NO`, 6,171 `NONE` — with **no MSI-H/MSI-L/MSS calls
whatsoever**, and there is no mutation table and no consensus-subtype table on disk. Substituted:
**breast carcinoma ER status by IHC** (`breast_carcinoma_estrogen_receptor_status`, 690 labelled
patients, 528 positive / 162 negative), with progesterone receptor as a second, weaker anchor.

**Pre-registered before the run** (`p1_evidence/inputs/PREREG_known_covariate.json`, written 01:45 UTC,
run started 01:47 UTC):

* ER expected band **[0.78, 0.92]**, point estimate 0.86 — Naik et al. *Nat Commun* 2020;11:5727 (ER
  AUC 0.92 internal ABCTB, **0.86 on TCGA external validation**); Rawat et al. *Sci Rep*
  2020;10:7275 (~0.86–0.89); Shamai et al. *JAMA Netw Open* 2019;2:e197700 (0.80); Couture et al.
  *npj Breast Cancer* 2018;4:30 (84% accuracy).
* PR expected band **[0.70, 0.85]**, point estimate 0.78 — Naik 2020 (0.81 internal, 0.75 TCGA).
* Grading rule fixed in advance: PASS only if the bootstrap CI of the **within-cancer** AUROC overlaps
  the band. Entirely below = under-recovery, FAIL. Entirely above = leak, FAIL.

**Result — ER, within-cancer AUROC, 1,000-draw bootstrap CI, 1,000 within-cancer label permutations:**

| artifact | state | adjustment | within-cancer AUROC | CI95 | null p95 | verdict |
|---|---|---|---:|---|---:|---|
| d2_h | **wsi_biology** | raw | **0.8781** | [0.8457, 0.9115] | 0.546 | **PASS** |
| d2_h | wsi_biology | cancer+TSS | 0.8714 | [0.8379, 0.9055] | 0.546 | PASS |
| d2_h | full_biology | raw | 0.9195 | [0.8901, 0.9455] | 0.544 | PASS |
| d2_h | rna_biology | raw | 0.9198 | [0.8911, 0.9459] | 0.545 | PASS |
| d2_i | **wsi_biology** | raw | **0.8667** | [0.8360, 0.8971] | 0.542 | **PASS** |
| d2_i | wsi_biology | cancer+TSS | 0.8644 | [0.8340, 0.8946] | 0.544 | PASS |
| d2_i | full_biology | raw | 0.9401 | [0.9127, 0.9640] | 0.545 | PASS (marginal) |
| d2_i | rna_biology | raw | 0.9395 | [0.9123, 0.9633] | 0.544 | PASS (marginal) |

PR, same protocol: `full_biology` 0.8341 [0.8009, 0.8684], `rna_biology` 0.8349 [0.8016, 0.8689] —
both inside [0.70, 0.85], PASS.

The two image-only numbers, **0.878 and 0.867**, sit essentially on the literature point estimate of
0.86 for TCGA external validation. The measured within-cancer chance level is 0.542–0.546, not 0.5, so
grading against an assumed 0.5 would have been wrong by four points.

**Two honest qualifications, both declared in the pre-registration before the run:**

1. **BRCA is a development cancer in the maximal split.** Its 690 labelled patients are train/val, so
   the control had to be run on `--partition all` and is *in-distribution* for the representation.
   This is a weaker control than the plan intended.
2. **Only one cancer carries the label**, so the within-cancer AUROC collapses to the BRCA AUROC and
   the lineage-guessing protection that statistic normally provides is not exercised (pooled and
   within-cancer are identical here for exactly that reason).

One thing to watch rather than celebrate: `full_biology` and `rna_biology` on d2_i reach 0.940, with
CI lower bounds of 0.913 — they clear the band only because the CI just overlaps 0.92. Those states
take RNA as input and ER status is close to a monotone function of *ESR1* expression, so scoring above
an H&E-derived literature band is expected by construction there and is not evidence about morphology.
The image-only `wsi_biology` row is the one that carries the control.

### In plain terms

Before running anything we wrote down, from four published papers, how well someone should be able to
guess a breast tumour's oestrogen-receptor status from the slide alone: somewhere between 0.78 and
0.92 on the usual scale, most likely about 0.86. Our image representation got 0.878. It landed almost
exactly where the literature said it should. That matters because a positive control that merely
"finds something" proves nothing — the pipeline could be finding an artefact. Landing on a
pre-committed external number is much harder to fake.

The caveats are that breast cancer was in our training set (so this is not a held-out test), and it is
the only cancer with the label, so this is one control rather than a battery.

### Meaning for the claim

* `HANDOFF_GATES.md`'s governing rule — a negative result is reportable only if the positive control
  passed in the same run, on the same data, through the same code path — is **discharged for the
  image-only channel** by this row. Every negative result in this session (T1.5's random-projection
  finding, the Track 2 nulls, the dilution decline) is therefore admissible rather than
  indistinguishable from a broken pipeline.
* The confound adjustment costs almost nothing here: 0.8781 → 0.8714 raw → adjusted. Consistent with
  the attenuation ≈ 1 finding, on an entirely different statistic.
* **P2 depends on this specifically**, being a negative-result paper.

### Files / commits

`v2/calibra/known_covariate_control.py` (inherited, audited, unchanged).
Pre-registration `p1_evidence/inputs/PREREG_known_covariate.json`; results
`p1_evidence/track1/covariate_er/`, `p1_evidence/track1/covariate_pr/`; covariate table
`p1_evidence/inputs/tcga_clinical_covariates.parquet` — all under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/`.
