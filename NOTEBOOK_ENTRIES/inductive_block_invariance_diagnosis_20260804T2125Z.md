## The inductive operator is correct and the assertion was too strong: one-row-vs-block agrees bit-for-bit in everything the operator owns, and to 4 eps (8.882e-16) in what BLAS owns — and to 0.000e+00 on all 14 deployed operators

**Logged:** 2026-08-04 21:25 UTC.
**Question:** `v2/tests/test_inductive_adjustment.py::test_one_row_at_a_time_equals_the_whole_block`
asserts `np.allclose(single, block[i:i+1], atol=0, rtol=0)` — bit-exactness — and fails.
Introduced at `e071d6c`; reported failing by an agent that did not touch it. Is the
implementation wrong, or is the assertion unattainable for a legitimate implementation?
**How obtained:** box `150.136.45.194`, `~/ws_ind` from
`git -c core.autocrlf=false archive HEAD` at `f8277e2`, per-file blob SHA-1 verified against
the local repo for `v2/calibra/inductive_adjustment.py` (`73c94aa1…`),
`v2/tests/test_inductive_adjustment.py` (`79dab7d3…`), `v2/calibra/residualise.py`
(`b2204ab5…`). `~/venv` (numpy 2.2.6, sklearn 1.7.2, scipy-openblas), threads capped to 1.
Cross-checked on the Windows workstation (numpy 2.4.3, scipy-openblas).

---

### 1. The verdict: the implementation is not wrong, and the assertion could not have held

The failure is real and reproducible, and it is **not** a leak of a cohort statistic. It is
`numpy` dispatching `(1, d) @ (d, k)` to a BLAS **GEMV** and `(m, d) @ (d, k)` to a blocked
**GEMM**, which sum the `d` products in different orders. Isolated away from this module
entirely, on dense random operands:

| operands | rows where `A[i:i+1] @ B != (A @ B)[i:i+1]` | max abs |
|---|---:|---:|
| `(12,19) @ (19,12)`, numpy 2.2.6 | 12 / 12 | 3.553e-15 |
| `(12,19) @ (19,12)`, numpy 2.4.3 | 12 / 12 | 2.665e-15 |

So the property `atol=0, rtol=0` was unattainable on the day it was written, on both numpy
versions available to this project. The test did not regress; it was introduced red.

### 2. What *is* bit-exact — and it is the property the operator is actually for

Every part of the map that this module owns is exactly row-wise, and this is asserted with
`array_equal`, not with a tolerance:

* **the design encoding.** `design(block)[i]` equals `design(row i)` bit-for-bit. `DesignSpec`
  freezes the one-hot levels and the numeric mean/s.d. against the reference, so nothing is
  recomputed from the incoming rows.
* **invariance to the rest of the block.** Hold row 0 fixed, replace all 11 companions with
  different cancers, different sites and values scaled by 100: row 0's adjusted coordinates are
  bit-identical (`max|diff| = 0.0`). A leaked mean, s.d. or centre cannot survive this; a BLAS
  blocking difference cannot be detected by it, because GEMM blocking depends on the shape and
  not on the values.
* **invariance to row order.** A permutation of the block returns every row unchanged, bit-for-bit.

Only the block's *size* moves the answer, and only in the last bit.

### 3. The measured deviation

Sweep: 3 synthetic cohorts × 3 residualiser settings `((5,1.0,42), (3,1.0,7), (10,2.0,1))` ×
block sizes `{2, 12, 64, 240}`, on the test's 19-column design, which contains a **continuous**
covariate (`purity`) and its missingness indicator.

| quantity | worst over the sweep |
|---|---:|
| `max abs` block vs single | **8.882e-16** (= 4 × eps) |
| `max rel` | 1.690e-13 (attained on a near-zero entry) |
| `max ULP` distance | 1024 (same reason) |
| max abs adjusted value | 7.237 |
| entries differing | ≈ 1.5–3% |

On the **14 persisted operators** at
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/runs_misc/tcga_operators/`
(design widths 22 and 99, 256 output columns, block sizes 8 and 64):

| quantity | value |
|---|---:|
| `max abs` block vs single | **0.000e+00** |
| `max ULP` | 0 |
| entries differing | 0 |
| design + composition bit-exact | True (14/14) |

The reason is structural, not luck: those designs are `cancer` + pooled `tss` and nothing else,
so every design row is one-hot. The products are exactly representable and every other term is
a zero, and re-associating a sum of exact values and zeros is exact. **Every deployed operator
— including the ones ALCHEMIST's external comparison runs through — satisfies the property
bit-for-bit, and no published number rests on the tolerance below.**

### 4. The tolerance, and why it cannot absorb a leak

A tolerance is defensible only against the thing it must still catch. Two *real* cohort-statistic
leaks, constructed on the same 12-row block the deviation was measured on:

| counterfactual leak | size |
|---|---:|
| re-centre on the new block (what transductive `cross_fitted_residuals` does) | 2.018e+00 |
| standardise `purity` by the new cohort's own mean/s.d. instead of the reference's | **1.285e-01** |
| BLAS block-size deviation (measured, above) | 8.882e-16 |

Set `NEW_ROW_BLOCK_INVARIANCE_ATOL = NEW_ROW_BLOCK_INVARIANCE_RTOL = 1e-12`. That is
**1.1e3 above** the measured deviation — room for a different BLAS — and **1.3e11 below** the
smaller of the two leaks. There is no cohort statistic that fits in the gap.

### 5. Considered and rejected: making it exact by construction

Adjusting rows one at a time inside `transform` would make the equality hold by construction and
would change *no* published number (the deployed designs already agree exactly). It was rejected
because it buys nothing that is not already guaranteed — §2 is the property, and it is already
bit-exact — at the cost of `n × K` Python-level `Ridge.predict` calls, and because touching the
numeric path of an operator whose whole purpose is that its numbers do not move is the wrong
trade for a last-bit cosmetic.

### 6. Do the 14 persisted operators need refitting? No.

The implementation is unchanged in its numerics — this commit adds constants, a measurement
function and docstrings, and rewrites tests. Independently, `fit.log` in the operator directory
records `identity=True maxabs=0.000e+00 round_trip=True` for all 14, i.e. the load-bearing
property (`adjust_reference` ≡ `cross_fitted_residuals`, bit-for-bit) holds as fitted, and §3
shows the new-row path is bit-exact for all 14 as well. Nothing to refit.

### 7. What changed

`v2/calibra/inductive_adjustment.py`
* module docstring: a new section stating the two-tier guarantee and its mechanism.
* `NEW_ROW_BLOCK_INVARIANCE_MEASURED_ATOL = 8.882e-16`,
  `NEW_ROW_BLOCK_INVARIANCE_MEASURED_ATOL_DEPLOYED = 0.0`,
  `COHORT_STATISTIC_LEAK_SMALLEST_MEASURED = 1.285e-01`,
  `NEW_ROW_BLOCK_INVARIANCE_ATOL = NEW_ROW_BLOCK_INVARIANCE_RTOL = 1e-12` — so no test and no
  document re-derives them inline.
* `block_invariance_deviation(operator, matrix, frame, ...)` — returns the deviation and the two
  bit-exactness flags, never a verdict.

`v2/tests/test_inductive_adjustment.py`
* `test_one_row_at_a_time_equals_the_whole_block` — keeps its name and its subject; now asserts
  the bit-exact half with `array_equal` and the BLAS half against both a scale-free bound
  (`64 · eps · max|value|`) and the module tolerance.
* `test_a_row_is_unmoved_by_whatever_else_is_in_its_block` — new, bit-exact.
* `test_a_categorical_only_design_is_bit_exact_block_versus_single` — new; reproduces the
  deployed operators' exactness without the box.
* `test_the_block_size_tolerance_cannot_absorb_a_cohort_statistic` — new; asserts the separation
  in §4 rather than asserting it in prose.

`v2/tests/test_inductive_adjustment.py`: 20 passed.

### Honest constraints

The 8.882e-16 is a worst case over the sweep in §3, not a proof of a bound; a different BLAS,
a wider design or a representation on a much larger scale could exceed it, which is why the
scale-free `64 · eps · max|value|` assertion sits alongside the flat tolerance. The exactness in
§3 for the deployed operators is a consequence of those designs being purely categorical and
would not survive adding a continuous covariate (e.g. purity or mRNAsi) to the adjustment set —
if that happens, the deviation moves from 0.000e+00 to the §3 regime, still ~1e11 below any leak.
