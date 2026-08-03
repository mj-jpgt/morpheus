# D3 — the morphology→molecular channel SURVIVES real ABSOLUTE purity entering the adjustment set, losing 2–6% of its excess over null; a rank-matched placebo moves it by +0.2%

**Logged:** 2026-08-03 13:30 UTC. Pre-registered in
`NOTEBOOK_ENTRIES/d3_d2p3_preregistration_20260803T1300Z.md`, committed `cd9b056`, **before** the run.

**How obtained:** Lambda A100 box `150.136.45.194`, workspace `~/ws_d3` (CPU only; the GPU training
queue was at 98% utilisation throughout and was not touched). All commands under
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`.

```
python -m morpheus.v2.calibra.run_calibra \
  --artifacts d2_h_seed42.npz d2_i_seed42.npz --targets frozen_rna_targets.npz \
  --partition test --levels 0.0,0.01,0.02,0.05,0.10,0.20,0.30,0.40,0.50 \
  --n-draws 40 --n-components 16 --n-permutations 2000 --seed 42 --n-jobs 10 \
  --purity-table tcga_absolute_purity_by_patient.csv --purity-source absolute \
  --require-rna-positive-control
```
Outputs on persistent NFS: `~/e0_run/d3/{main_absolute,placebo_rank_matched,purity,d2_3}/`.
Readout settings are deliberately identical to the 2026-08-03 02:30 UTC T1.2 ledger.

---

### Technical

#### The blocker is resolved with a real table, not a proxy

The ledger's "no TCGA consensus purity table on disk" is **correct** — an exhaustive
`find /lambda/nfs/geeg -iname "*purity*"` returns only the ten CPTAC
`cptac/tables/<cohort>/tumor_purity__washu.parquet` files, which are the external cohort, exactly as
`HANDOFF_PHASE_D.md` §D3 says. The clinical table `tcga_pancan_clinical.parquet` (6,429 × 746) has
no purity column either; all 746 were enumerated.

**But the box has outbound network, and the PanCanAtlas ABSOLUTE calls are openly downloadable.**

| field | value |
|---|---|
| file | `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` |
| url | `https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5` |
| sha256 | `f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd` |
| rows | 10,786 samples; 10,642 finite purity; range 0.08–1.00, median 0.65 |
| `purity_source` | **`absolute`** |

**No expression-derived fallback was used, and the circularity caveat does not apply.** ABSOLUTE
infers purity from SNP-array copy number and allelic fractions — from **DNA**, not from the RNA the
targets are built from. The previous agent's dx_year/age + TSS substitution is superseded and was
not used.

Purity is per-*sample*; the cohort is keyed by 12-character patient. Each patient was matched to
**the exact RNA source sample it was built from** via `tcga_pancan_rna.sample_map.parquet`, joined on
the 15-character `TCGA-XX-XXXX-NN` key, rather than guessing a sample type. 59 patients (0.6%) still
carried more than one purity-bearing sample and were resolved primary-first (01 before 03/06/02/05),
then lexicographically; **their within-patient purity spread has median 0.16, which is a real
limitation and is recorded, not smoothed away.** Coverage: **test 95.81%**, train 96.18%, val 97.05%.
The analysis is complete-case, `n = 2,650` of 2,766, and before/after run on the *identical* patient
set, so cohort composition cannot masquerade as a purity effect.

#### The result — quoted as paired within-run deltas, never absolute levels

`n = 2,650`; design 105 → **106** columns; `gates_pass: true`, `channel_gate_failures: []`,
`rna_positive_control_passed: true` (the RNA→RNA control clears its own pairing null in **both** the
before and after arm). Primary statistic is `excess_over_null_median`, because the 16-component
permutation null is **not zero** — it measures 0.1487–0.1504 here, consistent with the 0.1463–0.1483
on the T1.2 ledger.

| artifact | state | excess over null, before → after | retained | `adjusted_top_cca` Δ | `heldout_top_cca` Δ |
|---|---|---|---:|---:|---:|
| d2_h_seed42 | **wsi_biology** | 0.4729 → 0.4456 | **94.2%** | −0.0273 | −0.0343 |
| d2_i_seed42 | **wsi_biology** | 0.3374 → 0.3304 | **97.9%** | −0.0069 | −0.0195 |
| d2_h_seed42 | full_biology | 0.7419 → 0.6911 | 93.2% | −0.0505 | −0.0568 |
| d2_i_seed42 | full_biology | 0.6942 → 0.6739 | 97.1% | −0.0202 | −0.0220 |
| d2_h_seed42 | rna_biology * | 0.7388 → 0.6887 | 93.2% | −0.0497 | −0.0544 |
| d2_i_seed42 | rna_biology * | 0.7039 → 0.6765 | 96.1% | −0.0271 | −0.0292 |

\* RNA→RNA, circular by construction: the same-run positive control.

`permutation_p` stays at the 1/2001 resolution floor in **every** state, before and after. The
observed value stays 3.6× above `null_p95` in the worst case (0.4791 vs 0.1674).

**Verdict: SURVIVES.** The pre-declared bar was ≥80% of excess over null retained, observed still
above `null_p95`, and `permutation_p` still at floor. All three hold, with margin, in both artifacts.
**The falsifier did not fire.**

#### The rank-matched placebo — this is what makes the 94% readable

Adding purity adds a design column, so a channel drop could be the design gaining rank rather than
tumour content being removed. An otherwise identical run used a **placebo table: same patients, same
106-column design, ABSOLUTE purity values permuted across patients** (seed 20260803), declared
`--purity-source placebo_rank_matched`.

| artifact | state | real Δ | **placebo Δ** |
|---|---|---:|---:|
| d2_h_seed42 | wsi_biology | −0.0273 | **+0.0018** |
| d2_i_seed42 | wsi_biology | −0.0070 | **+0.0012** |
| d2_h_seed42 | full_biology | −0.0507 | **+0.0000** |
| d2_i_seed42 | full_biology | −0.0204 | **+0.0001** |
| d2_h_seed42 | rna_biology | −0.0501 | **+0.0002** |
| d2_i_seed42 | rna_biology | −0.0274 | **+0.0003** |

The placebo moves the channel by at most +0.002 — i.e. the one extra design column costs nothing.
**So the small loss is genuinely purity, and the surviving 94–98% is genuinely surviving.**

All **72** before-arm rows are identical to 0.0 between the two runs, which incidentally proves the
`run_calibra` edit made between them (the placebo source label) did not touch the numeric path.

#### The detail that matters most, and it is not the headline

**Purity does not preferentially explain the morphology channel.** On d2_h the morphology→molecular
channel loses 5.8% of its excess over null; the RNA→RNA circular control loses 6.8% and
`full_biology` 6.8%, on the same patients in the same run. On d2_i, morphology loses 2.1% against
3.9% for RNA→RNA. In both artifacts **the image channel is the *least* affected of the three.**
Purity adjustment shaves a common few percent off everything — which is what one expects if purity
is a mild global covariate of expression — rather than selectively deflating the image result.
That is the opposite of the "it was reading tumour content" hypothesis.

### In plain terms

The worry was that when the model looks at a slide and predicts molecular state, it might just be
counting how much tumour is in the picture — more tumour cells, more tumour-like expression, and
nothing biological in between. We found the real published measurement of how much tumour is in each
sample (from DNA copy number, so it is not derived from the expression we are trying to predict),
put it into the list of things we correct for, and measured the channel again on exactly the same
patients.

The channel lost about 2–6% of its strength and kept the rest. To check that even that small loss
was real, we repeated the whole thing with the purity numbers shuffled between patients — same
arithmetic, same number of correction terms, meaningless values. That version lost nothing. So the
2–6% is really purity, and the remaining 94–98% is really not.

The most useful detail: the image channel lost *less* than the RNA-to-RNA control did. If the image
result had been tumour-content in disguise, it should have lost the most.

### Meaning for the claim

- The D3 ledger row moves from NOT STARTED to **done, with a real `purity_source="absolute"`**, and
  its stated falsifier — "the channel dying when purity enters the adjustment set" — **did not fire**.
- `claim_guards.purity_confound` is discharged **for the states measured here** by
  `purity_in_adjustment_set`. It gates `legible_axis` and `gene_attribution` claims. It is **not**
  discharged for anything measured on a different cohort, partition or target block; this run is
  TCGA test split, 2,650 patients, 90 non-control RNA targets, 16 components.
- Escalation item 5 in `HANDOFF_PHASE_D.md` §5 ("the channel vanishes when purity enters") is **not**
  triggered.
- What this does **not** license: `no_external_cohort` and `composition_attribution` are untouched.
  Purity is not composition — knowing the tumour fraction says nothing about *which* non-tumour cells
  make up the rest, and T1.4 already showed 76–82% of the per-target channel is reproduced by
  covariate-matched random gene sets. D3 removes one alternative explanation; it does not promote the
  channel to a biological claim.

### Limitations, stated plainly

1. **116 of 2,766 test patients (4.2%) have no ABSOLUTE call** and are dropped. The complete-case
   cohort is therefore mildly non-random — ABSOLUTE fails more often on low-purity samples, which is
   exactly the tail that matters most for this question. This biases towards *under*-detecting a
   purity effect and should be read that way.
2. The 59 multi-sample patients were resolved by rule, and their within-patient purity spread has
   median 0.16.
3. `detection_floor` moved 0.30 → 0.40 in four of six state/artifact cells. It is measured on a
   coarse level grid, so a one-step move is within its resolution and is not interpreted.
4. Seed 42 only. Seeds 43/44 artifacts exist and this run would take ~7 CPU-minutes each; not done.

### Files / commits

- `~/e0_run/d3/purity/{abs.txt,tcga_absolute_purity_by_patient.csv,placebo_shuffled_purity_by_patient.csv}`
- `~/e0_run/d3/main_absolute/{task_rows.csv,calibra_summary.json,calibra_gates.json,calibra_protocol.json}`
- `~/e0_run/d3/placebo_rank_matched/` (same four files)
- `~/e0_run/d3/{run_main.sh,run_placebo.sh,logs/}`
- Code: `9bc7085` — `run_calibra --purity-source placebo_rank_matched` + test.
