## 2026-08-04 00:05 UTC — One effective rank, ten call sites: the canonical definition, and every historical instance recomputed under it

**Logged:** 2026-08-04 00:05 UTC. **How obtained:** source audit of the repository at
`research/rebase-vision`; full-text re-verification of Roy & Vetterli (EUSIPCO 2007) and RankMe
(arXiv:2210.02885v3) from primary PDFs; recomputation on frozen `.npz` artifacts on the A100
(`150.136.45.194`), **CPU only, thread-capped, no GPU** — a D1-m099 training pair was occupying the
GPU throughout (PIDs 217116/217117, 41 min elapsed at the time of measurement) and was not touched.
Scripts `~/ws_rank/recompute_rank.py` and `~/ws_rank/recompute_p1b.py`; outputs
`~/ws_rank/RANK_RECOMPUTE.json` and `~/ws_rank/RANK_RECOMPUTE_P1B.json`. The canonical implementation
used on the box is byte-identical to the committed blob
`9a3368a0d9d0c3c14a203353e806ff02227710af` = `85c0fa8:v2/calibra/spectral.py`
(sha256 `c403e3f3f7604af7ad5fe6f330eaf837fa0f969a1327cdf18bacd491d0b6e944`).

**Relationship to the concurrent entry.** `p2_competing_metrics_and_necessity_test_20260803T2326Z.md`
scores R1/R2/R3 as *selection rules* using the taxonomy this work established. §4 below adds a
**fourth** unstated choice that entry does not vary and which flips two of its cells; §7 states the
reconciliation required. Nothing here contradicts it.

---

### 1. The blocker, restated with the true count

`paper/P2_RANK_DRAFT.md` §3.1 records three mutually incompatible statistics named `effective_rank`.
A full AST + SVD scan of the tree finds **ten call sites carrying three statistics**, not three:

| statistic | definition | sites |
|---|---|---|
| **R1** | Roy & Vetterli: `exp(-Σ pᵢ ln pᵢ)`, `p = σ/Σσ`, column-centred | `v2/calibra/spectral.py:14`; torch duplicates at `v2/run_rank_ablation.py:35`, `v2/tests/test_stress_collapse.py:23`, `v2/calibra/e0_basis_transfer.py:432`; a fifth inline at `e0_basis_transfer.py:480` with a *different* (float32-eps) tolerance |
| **R2** | order-2 participation ratio `(Σσ)²/Σσ²`, centred | `v2/research/rebase/d1_audit.py:149` |
| **R3** | R2 on **L2-normalised rows** | `v2/research/rebase/d1_geometry_probe.py:50`; `v2/training.py:569`; `v2/runner.py:942`; `v2/research/rebase/d1_collapse_causal_test.py:75` |

Two of the R3 sites were not in the draft's census and both are **live abort thresholds**: the in-run
rank tripwire (`--rank-tripwire-minimum 4.0`) and the gate probe (`_require_rank_probe`, same bar).

**R1 and R2 are not variants of one statistic — they are Hill numbers of orders 1 and 2 of the same
singular-value distribution.** `(Σσ)²/Σσ² = 1/Σpᵢ²` exactly. Hill numbers are non-increasing in the
order, so **R2 ≤ R1 and R3 ≤ R1 for every matrix**, with equality only on a flat spectrum. The
disagreement therefore has a fixed sign; it is not noise. This is now asserted in
`v2/tests/test_effective_rank_canonical.py::test_order_two_is_never_above_order_one`.

---

### 2. The canonical definition, and why each choice was made

**CANONICAL = Roy & Vetterli Definition 1, order 1, column-centred, rows at their own norms.**
`v2/calibra/spectral.py`, `CANONICAL.label == "centred|order1"`.

Both source papers were re-read at full text before choosing, and the reading **changes what can be
claimed**:

| | Roy & Vetterli 2007 | RankMe (ICML 2023) |
|---|---|---|
| Definition | `erank(A) = exp{H(p)}`, `H = −Σ p_k log p_k`, `p_k = σ_k/‖σ‖₁`, natural log | adopts it, with `p_k = σ_k(Z)/‖σ(Z)‖₁ + ε` |
| Matrix centred before the SVD | **No** — "a complex-valued non-all-zero matrix A of size M × N whose singular value decomposition is given by A = UDV*"; no preprocessing appears anywhere in the paper | **No** — the paper is silent on centring and applies none |
| Rows L2-normalised | **No** | **No** |
| Near-zero σ | convention `0 log 0 = 0`; **no ε, no tolerance** | additive ε **outside** the division |
| ε value | n/a | **never stated for their Eq. (2)**; `10⁻⁷ for float32` is given only for the *contrasting* threshold-rank definition |
| Bounds | Property 1: `1 ≤ erank(A) ≤ rank(A) ≤ Q` | — |
| Quality/downstream claim | **None.** Closest is "assess the loss incurred by dimensionality reduction methods, such as PCA" | yes, that is the paper's thesis |

The draft's §2.1 statement that `spectral.py` "implements exactly this" is **not accurate and must be
amended**: `spectral.py` centres, and neither source paper does. Three consequences.

**(i) Order 1, not order 2.** R1 is the published statistic and the only one comparable to any number
outside this repository. R2 and R3 are a different quantity that happens to share the name.

**(ii) Centred — a deliberate deviation, to be stated wherever a value is quoted.** The reason is not
convention but that *uncentred effective rank is not a property of the representation's spread at
all*: it is a function of the column mean's magnitude relative to that spread, and it moves in **both**
directions.

- A large shared offset on an isotropic, near-full-rank representation drives uncentred erank to ~1 —
  it reads as total collapse when nothing has collapsed.
- On the collapse family `zᵢ = m + aᵢ·u` documented at
  `g26_centring_fix_20260803T0730Z.md`, with `m` comparable to the spread and not parallel to `u`,
  uncentred erank reads ~2 where there is exactly one direction of variation.
- Centred erank is *exactly* invariant to a shift; uncentred is not.

All four statements are pinned in
`test_effective_rank_canonical.py::test_uncentred_rank_is_a_function_of_the_mean_and_centred_rank_is_not`.
Centring is also what the majority of this project's quoted values used, so it is the choice that
keeps the historical record recomputable — see §5.

**(iii) No row normalisation.** It is in neither paper, it discards norm variation that is part of the
representation, and — see §4 — it is one of the choices that flips a published verdict.

**A defect fixed in passing.** The singular-value cut was **absolute** at `1e-12`, which breaks the
scale invariance Roy & Vetterli prove (Property 2: `erank(cA) = erank(A)`): scaling a matrix by 1e-9
silently emptied the spectrum and returned 0. It is now the standard LAPACK relative cut
`σ > σ_max·max(n,p)·eps`, with a separate *scale-relative* degeneracy floor at 1e-12 of the input's own
norm so that a representation collapsed to within float noise still scores 0 rather than the
dimensionality of the noise. **This change moved no historical number**: the maximum relative
difference between the old absolute cut and the new relative one, over all 68 recomputed
artifact × block combinations, is **0.000e+00**.

**Not reconciled, and stated as such.** RankMe's ε sits *outside* the division (verified at glyph
coordinates in the v3 PDF), so its `p_k` sum to `1 + min(N,K)·ε` and its statistic is not the
exponential of a Shannon entropy. No number in this repository is comparable to a published RankMe
value. (The concurrent entry measures the practical size of that gap on D2: `rankme_residualised`
23.391 against `effective_rank_residualised` 23.387.)

---

### 3. One implementation, everywhere

`v2/calibra/spectral.py` now holds the only definition. `effective_rank(x, *, centre, normalise_rows,
order, tolerance, variant)` defaults to `CANONICAL`; torch inputs are computed on-device in float64
without a host transfer, because the tripwire evaluates it every optimisation step, and the two
backends are asserted equal. `RANK_VARIANTS` names R1/R2/R3 and three further combinations so the
comparison in §4 runs through the same code as the production path.

All ten sites now import it. The two **live abort thresholds keep R3**, named explicitly in the call
(`RANK_VARIANTS["R3"]`), because the 4.0 bar was calibrated against R3 readings (6.92–7.25 healthy,
1.46–1.98 collapsed) and silently changing the statistic under an abort threshold would change which
runs are killed. Both now log the canonical value beside it (`biology_effective_rank_canonical`) so
the bar can be recalibrated on evidence.

`v2/tests/test_effective_rank_canonical.py` (13 tests) pins:

- hand-computed values on a matrix with known spectrum — `σ ∝ (2,1,1)` gives `erank = 2√2 =
  2.8284271247461903` and order-2 `= 8/3`, checked both uncentred on `diag(2,1,1)` and through a
  non-zero column mean, which simultaneously proves centring is applied;
- Roy & Vetterli Property 1 including both equality cases;
- `R2 ≤ R1` on 25 random spectra; scale invariance; torch ≡ numpy on all six variants;
- **object identity** across every importable call site;
- an **AST + SVD scan of the whole tree** that fails if a second definition, or any unallowlisted
  SVD-based rank, reappears.

Suite: **317 passed in 48.0 s** with thread caps (304 before; +13).

---

### 4. How much the choice matters — and a fourth choice nobody has been stating

The three statistics do not differ by a few percent. Over all 68 recomputed artifact × block
combinations:

| ratio | min | median | max |
|---|---:|---:|---:|
| R2 / R1 | 0.338 | 0.629 | 0.813 |
| R3 / R1 | 0.351 | 0.655 | 0.826 |
| R1 uncentred / R1 | 0.116 | 0.995 | 1.000 |

**R2 and R3 read between one third and four fifths of R1.** A rank quoted without its statistic is
uncertain by a factor of up to three.

**The fourth choice: which matrix.** Every CALIBRA readout computes `effective_rank(x)` on the **raw**
representation block (`run_calibra.py:140`) while the channel it is compared against is a
`heldout_top_cca` on the **confound-residualised** block. `d2_readout.py` reports both and the draft
quotes the **residualised** one for D2 and the **raw** one for the dilution and Phase-1b tables.
*Nothing in the draft says so, and the two are different representations.* On D2 seed 43 and seed 44
this choice flips the R3 verdict:

| seed | statistic | H (raw) | I (raw) | higher | H (resid.) | I (resid.) | higher |
|---|---|---:|---:|:--:|---:|---:|:--:|
| 43 | R1 | 24.674 | 28.800 | I | 28.772 | 34.117 | I |
| 43 | R2 | 11.720 | 11.111 | **H** | 13.227 | 12.972 | **H** |
| 43 | **R3** | 11.720 | 11.111 | **H** | 14.746 | 15.915 | **I** |
| 44 | R1 | 8.447 | 8.313 | H | 9.143 | 9.105 | H |
| 44 | R2 | 5.385 | 5.449 | I | 5.733 | 5.815 | I |
| 44 | **R3** | 5.385 | 5.449 | **I** | 6.564 | 6.302 | **H** |

(Arm H wins the channel in all three seeds, so "higher = H" is the ordering a rank rule would need.)

Note also that on the **raw** exported artifacts `R2 ≡ R3` exactly, because the exported `z_biology`
states are already L2-normalised by the model — the R2/R3 distinction only exists after
residualisation. That is why the distinction was invisible for so long.

**The headline.** The draft's single strongest sentence — §4.3, seed 43, *"the arm with the higher
effective rank loses (34.12 vs 28.77)"* — is **true under the canonical R1 and false under R2**, the
statistic §4.8.3 nominates for the pending D1 table (R2 gives H 13.23 > I 12.97, ordering correct).
The claim survives under the canonical definition. It does not survive a referee who runs the other
function in the same repository.

---

### 5. Every historical instance, recomputed

Canonical R1. "Original implementation" is what actually produced the published number, established by
reproducing it. **Every surviving instance reproduces to the digits published** — the canonicalisation
changed no value, only the label and the guarantee.

| # | instance | original value | original implementation | artifact | recomputed, canonical R1 | conclusion survives? |
|---|---|---|---|---|---|---|
| **6** | D2 seed 42, arm H / arm I | 23.39 / 14.87 | R1, **residualised** block | `~/e0_run/d2_v3/d2_v3_s42/artifacts/` | **23.3868 / 14.8675** | **Yes**, exactly |
| **6** | D2 seed 43, arm H / arm I | 28.77 / 34.12 | R1, residualised | `d2_v3_s43/artifacts/` | **28.7715 / 34.1168** | **Yes** under R1 — **but the inversion is statistic-specific** (§4): R2 orders H above I |
| **6** | D2 seed 44, arm H / arm I | 9.14 / 9.11 | R1, residualised | `d2_v3_s44/artifacts/` | **9.1426 / 9.1052** | **Yes** under R1; R2 reverses the sign of the 0.4% gap |
| **6** | re-export of arm H seed 42 | 8.68 | R1, residualised | `d2_v3/recovered_artifacts/` | **8.6809** | **Yes** |
| **4** | dilution, level 0.00 → 0.80 | 196.2 → 161.2, **−18%** | R1, **raw** block | `~/p1_out/dilution/dilution_foreign_tumour_pca256.npz` | **196.187 → 161.226**, −17.82% | **Value yes; the −18% and the "3.7×" do NOT survive as numbers.** See below |
| **2** | Phase 1b, `full` → `programme_only` | 38.48 → 32.06, **−17%** | R1, raw block | `/lambda/nfs/.../runs/v21_release_20260720_retry3_resume_safe/artifacts/` | **38.4834 → 32.0594**, −16.7% | **Yes.** Direction robust across every variant (−16.7% R1 raw, −18.9% R2/R3 raw, −18.4% R1 residualised, −22.6% R1 uncentred) |
| **3** | "rank pinned at 16/16" | 16/16 | `torch.linalg.matrix_rank`, **not an effective rank** | **`[NOT RECOMPUTABLE — artifact lost]`** — a train-time batch of 16, never exported | withdrawal in §4.6 **stands**, and is now stronger: the column was never any of R1/R2/R3 |
| **1** | decorrelation term, 49.9 → 103.3 | +107% | **unknown** | **`[NOT RECOMPUTABLE — artifact lost]`**; the cited `paper/.../RESULTS.md` does not exist | remains history only; must stay excluded from every count |
| **5** | D1-A, `programme_only` 9.81 / 10.47, `programme_free` 1.71 | R3 on live checkpoints | R3, raw, GPU forward | checkpoints survive at `~/e0_run/d1_v1/d1_{p42,p43,p44,f42}/last.pt` (`f43`, `f44` were never written — the gate refused them) | **`[NOT RECOMPUTED — requires a GPU forward pass; the GPU was running D1-m099 and was not contended for]`** | pending; the source's own prohibition (*"the contrastive arm never trained"*) is unaffected |

**The dilution result is the one that does not survive as quoted.** The published −18% is R1 on the
**raw** block. On the **residualised** block — the one the channel is actually read from — the same
canonical statistic falls only **210.179 → 203.667, −3.10%**. The abstract's *"miscalibrated by a
factor of 3.7"* therefore ranges over:

| block | statistic | level 0.00 → 0.80 | rank lost | miscalibration vs channel −66.7% |
|---|---|---|---:|---:|
| raw | **R1 (canonical, as published)** | 196.187 → 161.226 | 17.82% | **3.74×** |
| raw | R2 | 147.039 → 96.854 | 34.13% | 1.95× |
| raw | R3 | 150.493 → 105.215 | 30.09% | 2.22× |
| raw | R1 uncentred | 193.752 → 155.116 | 19.94% | 3.34× |
| **residualised** | **R1 (canonical, matched to the channel)** | **210.179 → 203.667** | **3.10%** | **21.53×** |
| residualised | R2 | 170.674 → 160.946 | 5.70% | 11.70× |
| residualised | R3 | 173.348 → 165.545 | 4.50% | 14.82× |

**The direction survives everywhere — rank always under-reports the information loss, in every
variant. The magnitude does not: it spans 1.95× to 21.5×.** The honest statement is that on matched
preprocessing the under-reporting is *larger* than published, which strengthens the paper's point and
destroys its quoted constant. `§4.2` and the abstract must be rewritten to quote a range with the
block named, or to quote the matched-preprocessing figure with the raw one as a sensitivity.

**Instance 5 is no longer pending on the rank side.** D1-B finished; all six artifacts exist at
`~/e0_run/d1_v2/artifacts/`. Canonical R1 on the residualised block:

| seed | `programme_only` | `programme_free` | ratio |
|---|---:|---:|---:|
| 42 | 29.3813 | 13.4184 | 2.19× |
| 43 | 24.6730 | 7.6003 | 3.25× |
| 44 | 11.1148 | 6.3937 | 1.74× |

The draft's §4.8.3 states the D1 rank column "will be **statistic R2**". It is now R1, and R2 would
have given roughly 60% of these values. The channel side is being computed concurrently
(`d2_compare` on `D1_PAIRED_BOOTSTRAP_STRATIFIED.json`) and belongs to that work, not this one.

---

### 6. Rank's own instability, assembled under one definition

Canonical R1 on the residualised held-out `wsi_biology` block unless marked. **Paired within-run
differences only for the channel; rank is quoted as a level because every practice this paper is about
quotes it as a level.**

| what varies | values | spread | statistic | provenance |
|---|---|---:|---|---|
| **retraining, same seed, same configuration** (D2 arm H, seed 42) | re-export of the surviving checkpoint **8.6809** vs retrained **23.3868** | **2.69×** | R1 | `~/e0_run/d2_v3/recovered_artifacts/` vs `d2_v3_s42/artifacts/`; `D2_RESULT.md` §4 |
| seeds 42/43/44 of one configuration (D2 arm H) | 23.3868 / 28.7715 / 9.1426 | 3.15× | R1 | `d2_v3_s{42,43,44}/artifacts/` |
| seeds 42/43/44 of one configuration (D2 arm I) | 14.8675 / 34.1168 / 9.1052 | **3.75×** | R1 | as above |
| seeds 42/43/44 (D1-B `programme_only`) | 29.3813 / 24.6730 / 11.1148 | 2.64× | R1 | `~/e0_run/d1_v2/artifacts/` |
| seeds 42/43/44 (D1-B `programme_free`) | 13.4184 / 7.6003 / 6.3937 | 2.10× | R1 | as above |
| **seeds 42/43/44 in flight, global step 200** (D1-B `programme_free`) | **7.545 / 45.646 / 12.194** | **6.05×** | **R3** (the tripwire statistic; not recomputable under R1 — the in-training states were never saved) | `~/e0_run/d1_v2/d1_f_seed{42,43,44}/train_metrics.jsonl`, key `train_rank_tripwire_observed`, epoch 11 |
| same, `programme_only` | 110.765 / 110.879 / 111.078 | **1.003×** | R3 | as above |

**The last two rows are the sharpest thing in this table.** At the *same* global step, on the *same*
configuration, with only the seed differing, the supervised arm reproduces to **0.3%** and the
contrastive arm spans **6.05×**. Rank's reproducibility is not a property of the statistic — it is a
property of the arm being measured, and it is worst exactly where the arm is interesting. A
reproducibility envelope measured on a well-behaved arm does not transfer.

**Against the differences a rank-based selection rule would act on:**

| comparison | rank ratio between arms | inside the 2.69× retraining envelope? |
|---|---:|---|
| D2 seed 42 | 1.573× | **yes** |
| D2 seed 43 | 1.186× | **yes** |
| D2 seed 44 | 1.004× | **yes** |
| Phase 1b (instance 2, single seed) | 1.200× | **yes** |
| D1-B seed 42 | 2.190× | **yes** |
| D1-B seed 44 | 1.738× | **yes** |
| D1-B seed 43 | 3.246× | no (only just) |

**Six of the seven between-arm rank differences this project has ever measured are smaller than the
spread of the same statistic when one configuration is retrained with the same seed.** Across those
same D2 runs the paired channel difference is −0.1325 / −0.1089 / −0.1226 — same sign 3/3, spread
0.024 on a mean of −0.121.

---

### 7. Required edits, and one reconciliation

To `paper/P2_RANK_DRAFT.md`:

1. §2.1 — *"Our `v2/calibra/spectral.py:14-29` implements exactly this"* is **wrong**: `spectral.py`
   centres and Roy & Vetterli do not. Replace with the deviation as stated in §2 above.
2. §3.1 — the census is ten sites, not three; add R2/R3 = order-2 Hill numbers and the fixed-sign
   inequality; add the **raw vs residualised** choice as a fourth unstated degree of freedom.
3. §4.2 and the abstract — the "−18%" and "3.7×" must carry the block and the range (§5).
4. §4.3 — add that seed 43's inversion holds under R1 and R3-on-residualised and **not** under R2.
5. §4.8.3 — the D1 rank column is R1, not R2; the D1-B rank values are in §5 above.
6. §5.2 — the row *"Rank and channel measured with one statistic across all instances | … the
   historical numbers cannot be recomputed"* is now **false** for instances 2, 4 and 6 and must be
   replaced by this entry's table. It remains true for instances 1 and 3.
7. Appendix B — one `effective_rank`, at `v2/calibra/spectral.py`.

**Reconciliation with `p2_competing_metrics_and_necessity_test_20260803T2326Z.md` §4.3.** That table's
R3 row (D2: OK / OK / MISS) is R3 on the **raw** block. On the **residualised** block R3 reads
OK / MISS / OK — both D2 cells that differ between R1 and R3 flip. Neither is wrong; the block was not
stated. That entry's §4.3 should name the block, and its conclusion — that the three statistics return
different verdicts on 3 of 6 pairs — is if anything understated.

### In plain terms

The project had been measuring "effective rank" three different ways, in ten places, without saying
which. Two of those places are switches that kill training runs. All three are now the same function,
which defaults to the definition in the published paper, and a test fails if anybody adds a fourth.

Every number the paper already published recomputes exactly, so nothing has to be retracted for being
mis-measured. But two things came out of the recomputation that do change the paper. First, the choice
of statistic is worth up to a factor of three, and on the paper's single strongest example — one seed
where the higher-rank arm loses — the other function in the same repository gives the opposite
ordering. Second, the dilution result's headline number, "rank falls 18% while the information falls
67%", was measured on a different version of the representation from the one the information was
measured on. Measured on the same version, rank falls 3%. The finding gets stronger and the number
quoted for it is wrong.

Two historical instances cannot be recomputed at all because their artifacts never existed, and one
more needs a GPU that is busy training.

### Meaning for the claim

The paper's central negative is intact and is now defensible against the obvious referee attack.
What must change is precision, not direction: every rank in the paper needs a statistic and a block
attached, the dilution magnitude must be quoted as a range, and the seed-43 inversion must be labelled
as holding under the canonical definition specifically. The instability contribution is strengthened —
the in-flight step-200 numbers show 6.05× against 1.003× at the same step, which is a cleaner
demonstration of the reproducibility floor than anything previously in §4.7.

### Files / commits

- `85c0fa8` — `v2/calibra/spectral.py` (canonical), `v2/tests/test_effective_rank_canonical.py`, and
  the nine call sites; suite 317 passed in 48.0 s
- `~/ws_rank/recompute_rank.py`, `~/ws_rank/recompute_p1b.py`,
  `~/ws_rank/RANK_RECOMPUTE.json`, `~/ws_rank/RANK_RECOMPUTE_P1B.json` on the A100
- Artifacts read: `~/e0_run/d2_v3/{d2_v3_s42,s43,s44}/artifacts/`, `~/e0_run/d2_v3/recovered_artifacts/`,
  `~/e0_run/d1_v2/artifacts/`, `~/e0_run/d1_v2/*/train_metrics.jsonl`,
  `~/p1_out/dilution/dilution_foreign_tumour_pca256.npz`,
  `/lambda/nfs/geeg/biorag3_persistent_20260711/runs/v21_release_20260720_retry3_resume_safe/artifacts/`
- Prior: `p2_competing_metrics_and_necessity_test_20260803T2326Z.md`,
  `g26_centring_fix_20260803T0730Z.md`, `rank_probe_repeat_variance_20260804T0900Z.md`

**One thing deliberately not done.** The `v22_a10_11v21_20260725` copies of the Phase-1b artifacts give
different values for the same artifact names (`full/wsi_biology` 42.0065 against 38.4834). They are
**not** a retraining-variance data point: `configuration_sha256` and `source_tree_sha256` both differ
(`ca8b57bd…`/`4630f4d9…` against `9f191bae…`/`1e2c4b20…`), so they are different configurations from
different code snapshots. Recording this so nobody quotes it as a repeat.
